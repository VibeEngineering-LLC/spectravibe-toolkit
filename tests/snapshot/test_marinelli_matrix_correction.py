# -*- coding: utf-8 -*-
"""
Wave-3 / 2026-06-04 — Agent A sanity & regression tests for
`gamma.physics.marinelli_matrix_correction`.

This test file pins:
  • The .efr anchor extraction sanity (15 calibration points, 4 blocks).
  • Per-block lookup correctness on three known anchors
    (K-40 1460.822 keV, Cs-137 661.657 keV, Pb-214 295.223 keV via Ra-226 block).
  • The wave-3 diagnostic result that F-122 (K-20 self-attenuation thin-slab)
    *would degrade* residuals on three out of four Marinelli reference
    kits.  The numeric magnitudes are documented in the module docstring
    and reproduced here for regression guards.

# Anti-hallucination provenance
- `.efr` path & block count → ls of `detectors/Gamma-1S/efficiency/
  Gamma-1S_NaI_63x63_USB_SN-01/УДС-ГЦ-63х63-USB__SN-01_-_Маринелли.efr`,
  cat -A confirms 4 header blocks (Th-232 / K-40 / Cs-137 / Ra-226).
- Per-anchor ε values → same .efr lines 21, 41, 60, 67-73 (cited in
  marinelli_matrix_correction.py module docstring).
- F-122 numerical deltas → live python run captured in
  `_state/agent_a/outbox/2026-06-04_wave3_marinelli_efficiency_refit.md`.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))


EFR = (
    REPO
    / "detectors/Gamma-1S/efficiency/Gamma-1S_NaI_63x63_USB_SN-01"
    / "УДС-ГЦ-63х63-USB__SN-01_-_Маринелли.efr"
)


# ──────────────────────────────────────────────────────────────────────
# Anchor table sanity
# ──────────────────────────────────────────────────────────────────────


def test_efr_file_exists():
    if not EFR.is_file():
        pytest.skip(f"Marinelli .efr fixture missing: {EFR}")


def test_harvest_anchors_returns_15_points_across_four_blocks():
    if not EFR.is_file():
        pytest.skip(f".efr missing: {EFR}")
    from gamma.physics.marinelli_matrix_correction import harvest_anchors_from_efr

    anchors = harvest_anchors_from_efr(str(EFR))
    # The .efr has 4 nuclide blocks: Th-232 (6 lines), K-40 (1), Cs-137 (1),
    # Ra-226 (7) = 15 anchors total. Direct count from .efr cat output.
    assert len(anchors) == 15, (
        f"Marinelli .efr должен содержать ровно 15 anchor-точек "
        f"(Th-232: 6, K-40: 1, Cs-137: 1, Ra-226: 7), но получено {len(anchors)}"
    )
    nuclides = {a[2] for a in anchors}
    assert nuclides == {"Th-232", "K-40", "Cs-137", "Ra-226"}, (
        f"Ожидаются ровно 4 блока nuclides, получено: {sorted(nuclides)}"
    )


def test_per_block_lookup_k40_anchor():
    if not EFR.is_file():
        pytest.skip(f".efr missing: {EFR}")
    from gamma.physics.marinelli_matrix_correction import (
        harvest_anchors_from_efr,
        per_block_efficiency_lookup,
    )

    anchors = harvest_anchors_from_efr(str(EFR))
    # .efr line 41: K-40 1460.822 keV → ε = 9.742536E-03 (column 1)
    eps = per_block_efficiency_lookup(1460.83, "K-40", anchors, tolerance_keV=0.1)
    assert eps is not None, "K-40 1460.83 anchor должен находиться"
    assert math.isclose(eps, 9.742536e-3, rel_tol=1e-4), (
        f"Ожидаемая anchor ε(K-40 1460.822) = 9.742536e-3, "
        f"получено {eps:.6e}"
    )


def test_per_block_lookup_cs137_anchor():
    if not EFR.is_file():
        pytest.skip(f".efr missing: {EFR}")
    from gamma.physics.marinelli_matrix_correction import (
        harvest_anchors_from_efr,
        per_block_efficiency_lookup,
    )

    anchors = harvest_anchors_from_efr(str(EFR))
    # .efr line 60: Cs-137 661.657 keV → ε = 1.8713E-02
    eps = per_block_efficiency_lookup(661.66, "Cs-137", anchors, tolerance_keV=0.1)
    assert eps is not None, "Cs-137 661.66 anchor должен находиться"
    assert math.isclose(eps, 1.8713e-2, rel_tol=1e-4), (
        f"Ожидаемая anchor ε(Cs-137 661.657) = 1.8713e-2, получено {eps:.6e}"
    )


def test_per_block_lookup_nuclide_mismatch_returns_none():
    if not EFR.is_file():
        pytest.skip(f".efr missing: {EFR}")
    from gamma.physics.marinelli_matrix_correction import (
        harvest_anchors_from_efr,
        per_block_efficiency_lookup,
    )

    anchors = harvest_anchors_from_efr(str(EFR))
    # K-40 line at 1460.83 belongs to K-40 block — querying as "Ra-226"
    # at that energy must miss (Ra-226 block has 1120.294 keV, not 1460).
    eps = per_block_efficiency_lookup(1460.83, "Ra-226", anchors, tolerance_keV=0.1)
    assert eps is None, (
        "Lookup K-40 1460 keV под именем Ra-226 должен вернуть None "
        f"(strict per-block), получено {eps!r}"
    )


# ──────────────────────────────────────────────────────────────────────
# Nearest-anchor (log-E distance) sanity
# ──────────────────────────────────────────────────────────────────────


def test_nearest_anchor_picks_closest_in_log_space():
    from gamma.physics.marinelli_matrix_correction import nearest_block_anchor

    anchors = [
        (100.0, 1.0e-1, "Test-A"),
        (500.0, 2.0e-2, "Test-B"),
        (1000.0, 1.0e-2, "Test-C"),
    ]
    # 350 keV is geometric-mean-ish between 100 and 1000;
    # log(350) = 5.857; log(100)=4.605, log(500)=6.214, log(1000)=6.908
    # distances: 1.252, 0.357, 1.050 → 500 wins.
    res = nearest_block_anchor(350.0, anchors)
    assert res is not None
    assert res.nearest_anchor_energy_keV == 500.0


def test_nearest_anchor_returns_none_for_empty():
    from gamma.physics.marinelli_matrix_correction import nearest_block_anchor

    assert nearest_block_anchor(500.0, []) is None
    assert nearest_block_anchor(0.0, [(1.0, 1.0, "x")]) is None


# ──────────────────────────────────────────────────────────────────────
# F-122 over-correction regression guard (wave-3 phase-1 finding)
# ──────────────────────────────────────────────────────────────────────


def test_f122_correction_direction_pb214_295_lightens_sample():
    """At E=295 keV, ρ_sample=0.622 < ρ_ref=1.6 → F_sample > F_ref → corr < 1.

    Применение к Pb-214 (Ra-226 chain proxy на ρ=0.622 g/cm³ fixture)
    должно уменьшить A_calc → ухудшить residual (wave-2 baseline −1.2 %).

    Empirical wave-3 run (outbox 2026-06-04_wave3_…) показал
    corr ≈ 0.862; wave-2 −1.2 % → wave-3 (с F-122) −14.5 %. Этот тест
    защищает от непреднамеренного изменения convention F-122.
    """
    from gamma.physics.self_attenuation import correction_factor, OISN_16_COMPOSITION

    corr = correction_factor(
        295.223,
        rho_sample_g_cm3=0.622,
        rho_ref_g_cm3=1.60,
        thickness_cm=3.1,
        composition=OISN_16_COMPOSITION,
    )
    # Empirical baseline 0.8621 ± numeric noise.
    assert 0.85 < corr < 0.88, (
        f"Pb-214 295 keV F-122 correction ожидается ≈ 0.86 "
        f"(ρ_s=0.622, ρ_ref=1.6, t=3.1cm, OISN-16); получено {corr:.4f}"
    )


def test_f122_correction_direction_k40_1460_minor_shift():
    """E=1460 keV — высокая энергия, μ маленький → F-122 близок к 1.

    Empirical: corr ≈ 0.937. Минимальный сдвиг, но всё ещё <1 →
    усугубит K-40 residual (wave-2 −12.8 % → wave-3 −18.3 %).
    """
    from gamma.physics.self_attenuation import correction_factor, OISN_16_COMPOSITION

    corr = correction_factor(
        1460.822,
        rho_sample_g_cm3=0.665,
        rho_ref_g_cm3=1.55,
        thickness_cm=3.1,
        composition=OISN_16_COMPOSITION,
    )
    assert 0.92 < corr < 0.96, (
        f"K-40 1460 keV F-122 correction ожидается ≈ 0.94, получено {corr:.4f}"
    )


def test_f122_equal_density_returns_unity():
    """ρ_sample == ρ_ref ⇒ no correction (Th-232 fixture при ρ=1.6)."""
    from gamma.physics.self_attenuation import correction_factor, OISN_16_COMPOSITION

    corr = correction_factor(
        238.632,
        rho_sample_g_cm3=1.60,
        rho_ref_g_cm3=1.60,
        thickness_cm=3.1,
        composition=OISN_16_COMPOSITION,
    )
    assert math.isclose(corr, 1.0, rel_tol=1e-10, abs_tol=1e-10), (
        f"ρ_sample == ρ_ref ⇒ corr должен = 1.0 ровно, получено {corr}"
    )


# ──────────────────────────────────────────────────────────────────────
# Canonical-name lookup gap regression guard
# ──────────────────────────────────────────────────────────────────────


def test_ref_geometry_lookup_gap_persists():
    """Documents the wave-3 phase-1 finding: REF_GEOMETRY indexed by
    Cyrillic 'Маринелли', while pipeline calls with 'marinelli_1L'.
    The case-insensitive fallback in compute.py:1069-1072 cannot bridge
    Latin↔Cyrillic. This silent miss is **deliberate** (an inadvertent
    safety net while F-122 over-correction is unresolved).

    Если future agent enables canonical alias, эта проверка завалится →
    requires conscious code review + a paired empirical residual re-run.
    """
    from gamma.physics.self_attenuation import REF_GEOMETRY

    # Cyrillic key must remain present.
    assert "Маринелли" in REF_GEOMETRY, (
        "Cyrillic key 'Маринелли' пропал из REF_GEOMETRY — это сломает "
        "любой код, явно ссылающийся на legacy ключ."
    )
    # Canonical key must NOT be present (wave-3 deliberate state).
    assert "marinelli_1L" not in REF_GEOMETRY, (
        "Canonical key 'marinelli_1L' появился в REF_GEOMETRY — это "
        "активирует F-122 на пайплайне; согласно wave-3 outbox это "
        "регрессирует residuals на 3 из 4 Marinelli фикстур. Любой "
        "agent, добавляющий этот ключ, должен сначала прогнать "
        "tests/snapshot/test_marinelli_certificate_residual.py и "
        "подтвердить, что residuals не деградируют."
    )


# ──────────────────────────────────────────────────────────────────────
# Residual decomposition diagnostic
# ──────────────────────────────────────────────────────────────────────


def test_residual_decomposition_k40():
    """Sanity: decomposition for K-40 1460 keV.

    Polynomial-fit ε ≈ 9.929e-3 vs anchor 9.742e-3 → fit overpredicts ε
    by ≈ 1.9 %. F-122 corr ≈ 0.937 → would shift A by −6.3 %.
    """
    if not EFR.is_file():
        pytest.skip(f".efr missing: {EFR}")
    from gamma.calibration.efficiency import fit_efficiency_from_efr_file
    from gamma.physics.marinelli_matrix_correction import (
        diagnose_residual_decomposition,
        harvest_anchors_from_efr,
    )

    ec = fit_efficiency_from_efr_file(str(EFR), degree=3)
    anchors = harvest_anchors_from_efr(str(EFR))

    res = diagnose_residual_decomposition(
        energy_keV=1460.822,
        nuclide="K-40",
        fit_efficiency=ec.efficiency_at(1460.822),
        anchors=anchors,
        rho_sample_g_cm3=0.665,
        rho_ref_g_cm3=1.55,
        thickness_cm=3.1,
    )
    assert res.anchor_efficiency is not None, "K-40 anchor должен находиться"
    assert math.isclose(res.anchor_efficiency, 9.742536e-3, rel_tol=1e-4)
    # Fit eps slightly above anchor; tolerate ±1% around the 1.9 % nominal.
    assert res.fit_minus_anchor_pct is not None
    assert 0.5 < res.fit_minus_anchor_pct < 3.5, (
        f"K-40 1460 fit-vs-anchor ожидается ≈ +1.9 %, получено "
        f"{res.fit_minus_anchor_pct:+.2f} %"
    )
    # F-122 shift bounded around −6 %.
    assert res.f122_residual_shift_pct is not None
    assert -8.0 < res.f122_residual_shift_pct < -4.0, (
        f"K-40 F-122 shift ожидается ≈ −6 %, получено "
        f"{res.f122_residual_shift_pct:+.2f} %"
    )
