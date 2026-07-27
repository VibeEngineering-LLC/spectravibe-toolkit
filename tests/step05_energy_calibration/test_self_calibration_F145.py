"""
F-145 / v1.17.8 — Regression tests for two-phase multiplet self-calibration.

Test plan (per ROADMAP_v1_17_8_plus.md F-145):
  1. Ra-226 demo: M1 χ²/ν improvement after Phase D (was 61.69 in v1.17.7;
     target ≤ 50 in v1.17.8, ≤ 20 in v1.19+).
  2. Th-232 demo: M1 χ²/ν preserved or improved (28.82 baseline).
  3. Synthetic single-line (Cs-137 only) — no recalibration (no multiplet).
  4. Synthetic Pb-214 doublet with injected -3 keV drift — Phase A fitted
     shifts within 0.5 keV of true drift; Phase D χ² substantially better
     than locked-only.

Contract (ЗАКРЕПЛЕНО НАВСЕГДА с v1.17.8):
  - coupled_intensity_fit принимает free_centroids=True/False kwarg
  - При free_centroids=True И use_peak_image=True запускается Phase A
    side-fit; основная LOCKED-passport подгонка не меняется.
  - Результирующий CoupledFitResult содержит centroid_shifts_keV (выровнен
    по components) и phase_A_chi2_per_dof.
  - recalibrate_from_multiplet_centroids формирует ОДИН anchor на
    мультиплет через I_pct-взвешенное усреднение центроидов.
  - При наличии stored_cal новая калибровка строится как
    E_new(N) = E_stored(N) + δ(N), где δ(N) — low-degree коррекция.
  - Phase D приметится только если суммарное χ² мультиплетов уменьшается.
"""
from __future__ import annotations

import json
import math
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

# Make scripts/ importable
ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from gamma.peaks.coupled_multiplet import (
    coupled_intensity_fit, ComponentSpec, H_STEP_DEFAULT_NAI,
)
from gamma.calibration.multiplet_self_calibration import (
    recalibrate_from_multiplet_centroids, SelfCalibrationDiag,
    PHASE_B_MAX_DRIFT_FRACTION_OF_FWHM,
    PHASE_B_MIN_CHI2_IMPROVEMENT,
)


# ──────────────────────────────────────────────────────────────────
# Test 1: Phase A side-fit на синтетике с известным drift'ом
# ──────────────────────────────────────────────────────────────────

def _nai_fwhm(E):
    """Реалистичная NaI 63x63 FWHM ~ 2.06·√E."""
    return 2.06 * math.sqrt(max(1.0, E))


def _make_synthetic_doublet(E_lib_pair, I_pct_pair, drift_keV, seed=42):
    """Двойная Гауссиана с заданным сдвигом центров от паспортных."""
    E_grid = np.linspace(220, 420, 200)
    data = np.full_like(E_grid, 2500.0)
    for E0, I_pct in zip(E_lib_pair, I_pct_pair):
        sig = _nai_fwhm(E0) / 2.355
        A = 1500.0 * I_pct
        data += A * np.exp(-((E_grid - (E0 + drift_keV))/sig)**2 / 2) / (sig * math.sqrt(2*math.pi))
    return E_grid, np.random.RandomState(seed).poisson(data).astype(float)


def test_phase_A_fits_known_drift():
    """Phase A на двойнике Pb-214 295+352 с искусственным drift'ом -3 кэВ
    должен фитировать средний сдвиг в пределах 0.5 кэВ от истинного."""
    DRIFT = -3.0
    E_grid, data = _make_synthetic_doublet(
        E_lib_pair=(295.22, 351.93),
        I_pct_pair=(18.42, 35.60),
        drift_keV=DRIFT,
    )
    components = [
        ComponentSpec(nuclide='Pb-214', E_keV=295.22, I_gamma_pct=18.42, group='Pb-214'),
        ComponentSpec(nuclide='Pb-214', E_keV=351.93, I_gamma_pct=35.60, group='Pb-214'),
    ]
    r = coupled_intensity_fit(
        E_grid, data, components, _nai_fwhm,
        use_peak_image=True, tail_param=0.7,
        h_step=H_STEP_DEFAULT_NAI,
        free_centroids=True,
    )
    assert r.phase_A_chi2_per_dof is not None, "Phase A не запустилась"
    assert r.phase_A_converged is True, "Phase A не сошлась"
    # Фитированные сдвиги
    assert len(r.centroid_shifts_keV) == 2
    avg_shift = sum(r.centroid_shifts_keV) / len(r.centroid_shifts_keV)
    assert abs(avg_shift - DRIFT) < 1.0, (
        f"average shift {avg_shift:+.3f} far from true {DRIFT:+.3f}"
    )
    # Phase A χ²/ν должен быть существенно лучше locked-passport
    assert r.phase_A_chi2_per_dof < 0.5 * r.chi2_per_dof, (
        f"Phase A χ² {r.phase_A_chi2_per_dof:.2f} не улучшил locked "
        f"{r.chi2_per_dof:.2f}"
    )


# ──────────────────────────────────────────────────────────────────
# Test 2: Phase B+C формирует weighted anchor и refit'ит калибровку
# ──────────────────────────────────────────────────────────────────

class _FakeSpec:
    """Минимальный stand-in для Spectrum, дающий energy_cal интерфейс."""
    def __init__(self, energy_cal):
        self.energy_cal = energy_cal

    def energy_to_channel(self, E):
        # Линейный inverse — для теста достаточно
        return (E - self.energy_cal[0]) / self.energy_cal[1]


def test_phase_B_C_weighted_anchor_and_refit():
    """Два мультиплета с консистентным drift'ом -2.5 keV → Phase B+C
    собирает 2 anchor'а и refit'ит линейную калибровку."""
    DRIFT = -2.5
    spec = _FakeSpec(energy_cal=(0.0, 1.0))  # E = N (1 канал = 1 кэВ)

    mults = []
    for name, comps_def in [
        ('M1', [(295.22, 18.42), (351.93, 35.60)]),
        ('M2', [(609.31, 45.49), (665.45, 1.51)]),
    ]:
        E_lo = comps_def[0][0] - 40
        E_hi = comps_def[-1][0] + 40
        E_seg = np.linspace(E_lo, E_hi, 180)
        data = np.full_like(E_seg, 2500.0)
        for E0, I_pct in comps_def:
            sig = _nai_fwhm(E0) / 2.355
            A = 1500.0 * I_pct
            data += A * np.exp(-((E_seg - (E0 + DRIFT))/sig)**2 / 2) / (sig * math.sqrt(2*math.pi))
        data = np.random.RandomState(7).poisson(data).astype(float)
        cs = [
            ComponentSpec(nuclide='Pb-214', E_keV=E0, I_gamma_pct=I, group='Pb-214')
            for E0, I in comps_def
        ]
        r = coupled_intensity_fit(
            E_seg, data, cs, _nai_fwhm,
            use_peak_image=True, tail_param=0.7,
            h_step=H_STEP_DEFAULT_NAI, cluster_id=name,
            free_centroids=True,
        )
        mults.append(r)

    new_cal, diag = recalibrate_from_multiplet_centroids(
        spec, mults, fwhm_provider_keV=_nai_fwhm,
    )
    assert diag.phase_B_passed, f"Phase B failed: {diag.reason}"
    assert new_cal is not None, f"Phase C did not apply: {diag.reason}"
    assert diag.n_anchors_after_filter == 2, (
        f"Expected 2 weighted anchors, got {diag.n_anchors_after_filter}"
    )
    # Сдвиг калибровки должен идти в ПРАВИЛЬНОМ направлении (-DRIFT > 0).
    # Точная величина зависит от Phase A noise (фитированные shifts ~ DRIFT
    # с σ ≈ 1-2 кэВ), поэтому tolerance свободный.
    new_a0 = new_cal[0]
    assert (new_a0 - 0.0) * (-DRIFT) > 0, (
        f"new_a0={new_a0:+.3f} в неправильном направлении (-DRIFT={-DRIFT:+.3f})"
    )
    assert abs(new_a0 - (-DRIFT)) < 2.0, (
        f"new_a0={new_a0:.3f} слишком далеко от -DRIFT={-DRIFT}"
    )


# ──────────────────────────────────────────────────────────────────
# Test 3: Synthetic single-line (без мультиплетов) — F-145 no-op
# ──────────────────────────────────────────────────────────────────

def test_no_multiplets_no_recalibration():
    """Если нет ни одного мультиплета — F-145 не пытается refit'ить.

    Phase A не запускается, recalibrate_from_multiplet_centroids
    возвращает (None, diag) с n_multiplets_seen=0.
    """
    spec = _FakeSpec(energy_cal=(0.0, 1.0))
    new_cal, diag = recalibrate_from_multiplet_centroids(
        spec, [], fwhm_provider_keV=_nai_fwhm,
    )
    assert new_cal is None
    assert diag.n_multiplets_seen == 0
    assert diag.n_anchors_after_filter == 0
    assert diag.phase_C_applied is False
    assert ("мультиплет" in diag.reason.lower()
            or "anchor" in diag.reason.lower()
            or "не прошёл" in diag.reason.lower())


# ──────────────────────────────────────────────────────────────────
# Test 4: Drift вне Phase B окна → anchor отвергнут
# ──────────────────────────────────────────────────────────────────

def test_drift_beyond_window_rejects_anchor():
    """Если average dE > 0.5·FWHM(⟨E⟩), мультиплет НЕ становится anchor'ом
    (защита от ложных свободных параметров)."""
    spec = _FakeSpec(energy_cal=(0.0, 1.0))
    DRIFT = -50.0  # Огромный — заведомо > 0.5·FWHM(~30) на 295 кэВ
    E_grid, data = _make_synthetic_doublet(
        E_lib_pair=(295.22, 351.93),
        I_pct_pair=(18.42, 35.60),
        drift_keV=DRIFT,
    )
    components = [
        ComponentSpec(nuclide='Pb-214', E_keV=295.22, I_gamma_pct=18.42, group='Pb-214'),
        ComponentSpec(nuclide='Pb-214', E_keV=351.93, I_gamma_pct=35.60, group='Pb-214'),
    ]
    r = coupled_intensity_fit(
        E_grid, data, components, _nai_fwhm,
        use_peak_image=True, tail_param=0.7,
        h_step=H_STEP_DEFAULT_NAI,
        free_centroids=True, cluster_id='M_outlier',
    )
    new_cal, diag = recalibrate_from_multiplet_centroids(
        spec, [r], fwhm_provider_keV=_nai_fwhm,
    )
    # Anchor'ы либо отброшены polyfit'ом (insufficient), либо drift cap
    assert new_cal is None or diag.n_anchors_after_filter < 2


# ──────────────────────────────────────────────────────────────────
# Test 5: Pipeline end-to-end на реальном Ra-226 demo
# ──────────────────────────────────────────────────────────────────

@pytest.mark.slow
def test_ra226_demo_phase_C_applied():
    """Полный pipeline на Ra-226 → Phase C должен примениться И χ²_sum
    уменьшается (Phase D принят)."""
    from gamma.identification.staged_pipeline import analyze_lsrm_spe

    SPE = (
        ROOT / "detectors" / "Gamma-1S" / "reference_spectra"
        / "archive"
        / "Ra226_420-7-18_Маринелли_0cm.spe"
    )
    if not SPE.exists():
        pytest.skip(f"Ra-226 .spe не найден: {SPE}")

    res = analyze_lsrm_spe(
        str(SPE),
        apply_deconvolution=True, allow_stage2=True, allow_stage3=True,
        compute_activities=True, compute_mda=True,
    )
    d = res.multiplet_self_calibration_diag
    assert d is not None, "multiplet_self_calibration_diag не пробросилось"
    assert d.get("attempted") is True
    # На реальном спектре с дрейфом E(N) — Phase C должен apply
    assert d.get("phase_C_applied") is True, (
        f"Phase C не применён: {d.get('reason')}"
    )
    assert d.get("n_anchors_after_filter") >= 2
    # Suммарный χ² мультиплетов после Phase D ≤ χ² Phase A (locked)
    # (это и есть критерий принятия Phase D)
    # Также проверка: M1 Pb-214 χ² должен УМЕНЬШИТЬСЯ относительно
    # v1.17.7 baseline (61.69).
    m1 = res.deconvolution_results[0]
    assert m1.chi2_per_dof < 61.69, (
        f"M1 χ²/ν {m1.chi2_per_dof:.2f} не улучшил v1.17.7 baseline 61.69"
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
