"""Regression — F-117 связанная подгонка мультиплета M1 (Th-232 fixture).

Контракт зафиксирован в references/demo_contract_v1_17_2/multiplet_M1_coupled.json:
  • ROI 754-1114 keV (каналы 260-380)
  • χ²/ν ≈ 17.02, closure ≈ -0.92%
  • Площади: Ac-228 911 ≈ 116996, 964.77 ≈ 22628, 969 ≈ 71649
  • Связь по интенсивностям: A(Ac-228) общий

Допуск на χ²/ν ослаблен до 50 (из-за вариаций FWHM-модели между
v1.17.2 hand-crafted и v1.17.5 default_NaI). Допуск на площади
ослаблен до 15% по той же причине.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

from gamma.io.readers import read_spectrum
from gamma.identification.staged_pipeline import (
    build_fwhm_model, fwhm_keV_at_energy,
)
from gamma.peaks.coupled_multiplet import (
    coupled_intensity_fit, ComponentSpec,
)


_FIXTURE = (
    "detectors/Gamma-1S/reference_spectra/"
    "archive/"
    "Th232_420-7-17_Маринелли_0cm.spe"
)


def _maybe_skip():
    if not Path(_FIXTURE).is_file():
        print(f"  ⚠ skipping (fixture missing): {_FIXTURE}")
        return True
    return False


def test_coupled_fit_M1_chi2_and_closure():
    if _maybe_skip():
        return
    spec = read_spectrum(_FIXTURE)
    fwhm_model, _ = build_fwhm_model(spec)
    fwhm_at = lambda E: fwhm_keV_at_energy(fwhm_model, E)
    lo, hi = 260, 380
    energies = np.array([spec.channel_to_energy(c) for c in range(lo, hi)])
    counts = spec.counts[lo:hi].astype(np.float64)
    comps = [
        ComponentSpec("Ac-228", 911.204, 25.8, group="Ac-228"),
        ComponentSpec("Ac-228", 964.77,  4.99, group="Ac-228"),
        ComponentSpec("Ac-228", 968.971, 15.8, group="Ac-228"),
        ComponentSpec("Tl-208", 860.6,   4.5,  group=""),
    ]
    res = coupled_intensity_fit(
        energies, counts, comps, fwhm_at,
        continuum="step_linear", roi_low_ch=lo, cluster_id="M1c",
    )
    # F-117 contract: χ²/ν conservatively bounded
    assert res.chi2_per_dof < 50.0, (
        f"M1 χ²/ν={res.chi2_per_dof:.2f} > 50 (gold=17.02)"
    )
    # |closure| < 5 %
    assert abs(res.closure_pct) < 5.0, (
        f"M1 closure {res.closure_pct:.2f}% не в допуске ±5%"
    )
    print(f"  ✓ test_coupled_fit_M1_chi2_and_closure "
          f"(χ²/ν={res.chi2_per_dof:.2f}, closure={res.closure_pct:.2f}%)")


def test_coupled_fit_M1_areas_within_15pct():
    if _maybe_skip():
        return
    spec = read_spectrum(_FIXTURE)
    fwhm_model, _ = build_fwhm_model(spec)
    fwhm_at = lambda E: fwhm_keV_at_energy(fwhm_model, E)
    lo, hi = 260, 380
    energies = np.array([spec.channel_to_energy(c) for c in range(lo, hi)])
    counts = spec.counts[lo:hi].astype(np.float64)
    comps = [
        ComponentSpec("Ac-228", 911.204, 25.8, group="Ac-228"),
        ComponentSpec("Ac-228", 964.77,  4.99, group="Ac-228"),
        ComponentSpec("Ac-228", 968.971, 15.8, group="Ac-228"),
        ComponentSpec("Tl-208", 860.6,   4.5,  group=""),
    ]
    res = coupled_intensity_fit(
        energies, counts, comps, fwhm_at,
        continuum="step_linear", roi_low_ch=lo, cluster_id="M1c",
    )
    # Gold targets from references/demo_contract_v1_17_2/multiplet_M1_coupled.json
    GOLD = {
        911.204: 116995.79,
        964.77:   22628.25,
        968.971:  71648.58,
    }
    for cf in res.components:
        E_key = round(cf.E_keV, 3)
        if E_key not in GOLD:
            continue  # Tl-208 → 0 (skip)
        gold = GOLD[E_key]
        rel = abs(cf.area - gold) / gold
        assert rel < 0.15, (
            f"area({cf.E_keV:.2f})={cf.area:.0f} vs gold {gold:.0f}, "
            f"Δ={rel:.1%} > 15 %"
        )
    # Связь по интенсивностям: соотношения площадей должны примерно
    # равняться соотношениям I_pct (доказательство coupling)
    a_911 = next(c.area for c in res.components if abs(c.E_keV - 911.204) < 0.5)
    a_964 = next(c.area for c in res.components if abs(c.E_keV - 964.77) < 0.5)
    a_969 = next(c.area for c in res.components if abs(c.E_keV - 968.971) < 0.5)
    r_964 = a_964 / a_911
    r_969 = a_969 / a_911
    r_964_lib = 4.99 / 25.8
    r_969_lib = 15.8 / 25.8
    assert abs(r_964 - r_964_lib) / r_964_lib < 0.01, (
        f"coupling violated: a(964)/a(911)={r_964:.4f} vs lib {r_964_lib:.4f}"
    )
    assert abs(r_969 - r_969_lib) / r_969_lib < 0.01, (
        f"coupling violated: a(969)/a(911)={r_969:.4f} vs lib {r_969_lib:.4f}"
    )
    print(f"  ✓ test_coupled_fit_M1_areas_within_15pct "
          f"(a(911)={a_911:.0f}, ratios match libs)")


def test_phantom_inclusive_recovers_absorbed_flux():
    """BUG-32ζ / task #82 — phantom-inclusion-in-fit.

    Synthetic spectrum: 3 неблендяющихся Gauss на 800/810/820 keV с
    «инжектированными» площадями 10000 / 3000 / 5000. Сценарий моделирует
    F-387.1 top-K cap'инг: «kept» компоненты на 800 + 820, «phantom» на
    810.

    Контракт:
      (a) phantom_inclusive=False (default): phantom 810 НЕ в fit'е;
          его flux ~3000 counts абсорбируется в kept-компоненты →
          area(800) или area(820) значительно (>10%) выше реальной.
      (b) phantom_inclusive=True: phantom получает свободную площадь
          с Tikhonov zero-prior penalty; восстанавливает свою площадь
          в пределах ~3000 ± 2σ_Poisson, а kept-площади возвращаются
          к инжектированным (≤10% bias).
    """
    from gamma.peaks.coupled_multiplet import (
        coupled_intensity_fit, ComponentSpec,
    )

    # ─── synthetic spectrum: 3 unit-area Gaussians × площади ─────────
    rng = np.random.default_rng(seed=20260603)
    SQRT_2PI = np.sqrt(2.0 * np.pi)

    # FWHM(E) NaI-подобный: R(662) ≈ 7%, FWHM(E) ∝ sqrt(E)
    def fwhm_at(E):
        return 0.07 * np.sqrt(E * 662.0)

    # Energy grid: 700-950 keV, 1 keV/channel.
    # FWHM(800) ≈ 51 keV; линии на 800/820/870 → min separation 20 кэВ
    # (sep/FWHM ≈ 0.39) — типичный unresolved multiplet (Rayleigh-near,
    # F-387.1 unresolved criterion как раз срабатывает при sep < ~1·FWHM).
    # При sep/FWHM < 0.3 deconvolution становится near-singular и
    # phantom-flux фундаментально не разрешим даже без penalty; при
    # sep/FWHM > 0.4 — phantom recovers ≥95%.
    E = np.linspace(700.0, 950.0, 251)

    E_lines = [800.0, 820.0, 870.0]
    A_inject = [10000.0, 3000.0, 5000.0]   # injected areas (counts)
    baseline = 200.0                        # flat continuum

    counts_clean = np.full_like(E, baseline)
    for E0, A in zip(E_lines, A_inject):
        sigma = fwhm_at(E0) / 2.355
        bin_w = E[1] - E[0]
        gauss = np.exp(-0.5 * ((E - E0) / sigma) ** 2) / (sigma * SQRT_2PI)
        counts_clean += A * gauss * bin_w

    # Poisson noise (realistic)
    counts = rng.poisson(np.maximum(counts_clean, 0.0)).astype(np.float64)

    # ─── kept components (top-K=2 by intensity: 800 (I=100) и 870 (I=50)) ─
    # 820 (I=30) — phantom, демотирован top-K cap'ом.
    # Group="" → каждый получает собственную свободную площадь
    # (моделируем «kept» как независимые fit-target'ы, не coupled через
    # один нуклид — чтобы изолировать BUG-32ζ механизм от F-117 coupling).
    kept = [
        ComponentSpec("Nuc-A", 800.0, 100.0, group=""),
        ComponentSpec("Nuc-A", 870.0,  50.0, group=""),
    ]
    phantom = [
        ComponentSpec("Nuc-A", 820.0,  30.0, group=""),
    ]

    # ─── (a) phantom_inclusive=False (current behaviour) ──────────────
    res_off = coupled_intensity_fit(
        E, counts, kept, fwhm_at,
        continuum="linear",
        cluster_id="phantom_off",
        # phantom_components=() (default), lambda_phantom_rel=0 (default)
    )
    a800_off = next(c.area for c in res_off.components if abs(c.E_keV - 800.0) < 0.1)
    a870_off = next(c.area for c in res_off.components if abs(c.E_keV - 870.0) < 0.1)
    bias_off = max(
        abs(a800_off - A_inject[0]) / A_inject[0],
        abs(a870_off - A_inject[2]) / A_inject[2],
    )

    # ─── (b) phantom_inclusive=True ───────────────────────────────────
    res_on = coupled_intensity_fit(
        E, counts, kept, fwhm_at,
        continuum="linear",
        cluster_id="phantom_on",
        phantom_components=phantom,
        lambda_phantom_rel=1e-3,
    )
    a800_on = next(c.area for c in res_on.components if abs(c.E_keV - 800.0) < 0.1)
    a870_on = next(c.area for c in res_on.components if abs(c.E_keV - 870.0) < 0.1)
    a820_on = next(c.area for c in res_on.components if abs(c.E_keV - 820.0) < 0.1)
    bias_on = max(
        abs(a800_on - A_inject[0]) / A_inject[0],
        abs(a870_on - A_inject[2]) / A_inject[2],
    )

    # ─── контрактные проверки ─────────────────────────────────────────
    # 1) Phantom-on восстанавливает phantom-площадь в разумных пределах.
    #    2σ Poisson ≈ 2·√3000 ≈ 110 → 3.7% от инжектированных. Даём
    #    более широкий 25% допуск из-за остаточных корреляций с соседями
    #    при sep/FWHM ≈ 0.4 (унresolved multiplet с deconvolution-noise).
    rel_820 = abs(a820_on - A_inject[1]) / A_inject[1]
    assert rel_820 < 0.25, (
        f"phantom area 820 recovery failed: got {a820_on:.0f}, "
        f"injected {A_inject[1]:.0f}, Δ={rel_820:.1%}"
    )
    # 2) Phantom-on уменьшает bias на kept-компонентах vs phantom-off
    #    (phantom flux больше не абсорбируется в соседей)
    assert bias_on < bias_off, (
        f"phantom_inclusive did NOT reduce kept-component bias: "
        f"bias_off={bias_off:.1%}, bias_on={bias_on:.1%}"
    )
    # 3) Bias_on ≤ 5% — phantom-on восстанавливает kept к инжектированным
    #    в пределах ~Poisson-σ (для 5000-counts линии σ ≈ 1.4%).
    assert bias_on < 0.05, (
        f"phantom-on bias on kept still > 5%: {bias_on:.1%}"
    )
    # 4) Default (phantom_inclusive=False) — bias_off демонстрирует
    #    значительное phantom flux absorption (>10%); если меньше — наш
    #    syntheticgeometry не triggerит BUG-32 симптом и test не
    #    валидирует contract.
    assert bias_off > 0.10, (
        f"phantom_off bias too small ({bias_off:.1%}) — synthetic "
        f"geometry не triggerит BUG-32 absorption"
    )
    print(
        f"  ✓ test_phantom_inclusive_recovers_absorbed_flux "
        f"(bias_off={bias_off:.1%}, bias_on={bias_on:.1%}, "
        f"a(820|on)={a820_on:.0f} vs injected {A_inject[1]:.0f})"
    )


def test_phantom_inclusive_default_off_is_backcompat():
    """BUG-32ζ regression-guard: default kwargs (phantom_components=(),
    lambda_phantom_rel=0.0) дают результат БИТНО-идентичный pre-BUG-32ζ
    path. Если user не передал phantom-параметры — control flow тот же.

    Проверяется через unmodified M1 fixture path.
    """
    if _maybe_skip():
        return
    spec = read_spectrum(_FIXTURE)
    fwhm_model, _ = build_fwhm_model(spec)
    fwhm_at = lambda E: fwhm_keV_at_energy(fwhm_model, E)
    lo, hi = 260, 380
    energies = np.array([spec.channel_to_energy(c) for c in range(lo, hi)])
    counts = spec.counts[lo:hi].astype(np.float64)
    comps = [
        ComponentSpec("Ac-228", 911.204, 25.8, group="Ac-228"),
        ComponentSpec("Ac-228", 964.77,  4.99, group="Ac-228"),
        ComponentSpec("Ac-228", 968.971, 15.8, group="Ac-228"),
        ComponentSpec("Tl-208", 860.6,   4.5,  group=""),
    ]
    # Pre-BUG-32ζ call (без новых kwargs)
    res_orig = coupled_intensity_fit(
        energies, counts, comps, fwhm_at,
        continuum="step_linear", roi_low_ch=lo, cluster_id="M1c",
    )
    # New call с явными defaults (phantom_components=(), lambda=0.0)
    res_default = coupled_intensity_fit(
        energies, counts, comps, fwhm_at,
        continuum="step_linear", roi_low_ch=lo, cluster_id="M1c",
        phantom_components=(),
        lambda_phantom_rel=0.0,
    )
    # χ²/ν и площади должны совпадать ДО последнего бита
    assert abs(res_orig.chi2_per_dof - res_default.chi2_per_dof) < 1e-9, (
        f"back-compat violated: chi2 differs "
        f"{res_orig.chi2_per_dof} vs {res_default.chi2_per_dof}"
    )
    for c1, c2 in zip(res_orig.components, res_default.components):
        assert abs(c1.area - c2.area) < 1e-9, (
            f"back-compat violated: area({c1.E_keV})="
            f"{c1.area} vs {c2.area}"
        )
    print(f"  ✓ test_phantom_inclusive_default_off_is_backcompat")


if __name__ == "__main__":
    test_coupled_fit_M1_chi2_and_closure()
    test_coupled_fit_M1_areas_within_15pct()
    test_phantom_inclusive_recovers_absorbed_flux()
    test_phantom_inclusive_default_off_is_backcompat()
    print("All M1 coupled-fit tests passed.")
