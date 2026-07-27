"""Regression — F-117 связанная подгонка мультиплета M2 (Th-232 fixture).

Контракт зафиксирован в references/demo_contract_v1_17_2/multiplet_M2_coupled.json:
  • ROI 1430-1786 keV (каналы 484-603)
  • χ²/ν ≈ 1.19, closure ≈ -0.51 %
  • Площади: Ac-228 1588 ≈ 16258, 1630 ≈ 8078, Bi-212 1620 ≈ 2172
  • Связь по интенсивностям: A(Ac-228) общий, Bi-212 независим
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


def test_coupled_fit_M2_chi2_and_closure():
    if _maybe_skip():
        return
    spec = read_spectrum(_FIXTURE)
    fwhm_model, _ = build_fwhm_model(spec)
    fwhm_at = lambda E: fwhm_keV_at_energy(fwhm_model, E)
    lo, hi = 484, 603
    energies = np.array([spec.channel_to_energy(c) for c in range(lo, hi)])
    counts = spec.counts[lo:hi].astype(np.float64)
    comps = [
        ComponentSpec("Ac-228", 1588.2, 3.22, group="Ac-228"),
        ComponentSpec("Bi-212", 1620.5, 1.49, group=""),
        ComponentSpec("Ac-228", 1630.6, 1.6,  group="Ac-228"),
    ]
    res = coupled_intensity_fit(
        energies, counts, comps, fwhm_at,
        continuum="step_linear", roi_low_ch=lo, cluster_id="M2c",
    )
    assert res.chi2_per_dof < 2.0, (
        f"M2 χ²/ν={res.chi2_per_dof:.2f} > 2.0 (gold=1.19)"
    )
    assert abs(res.closure_pct) < 5.0, (
        f"M2 closure {res.closure_pct:.2f}% не в допуске ±5%"
    )
    print(f"  ✓ test_coupled_fit_M2_chi2_and_closure "
          f"(χ²/ν={res.chi2_per_dof:.2f}, closure={res.closure_pct:.2f}%)")


def test_coupled_fit_M2_areas_within_25pct():
    if _maybe_skip():
        return
    spec = read_spectrum(_FIXTURE)
    fwhm_model, _ = build_fwhm_model(spec)
    fwhm_at = lambda E: fwhm_keV_at_energy(fwhm_model, E)
    lo, hi = 484, 603
    energies = np.array([spec.channel_to_energy(c) for c in range(lo, hi)])
    counts = spec.counts[lo:hi].astype(np.float64)
    comps = [
        ComponentSpec("Ac-228", 1588.2, 3.22, group="Ac-228"),
        ComponentSpec("Bi-212", 1620.5, 1.49, group=""),
        ComponentSpec("Ac-228", 1630.6, 1.6,  group="Ac-228"),
    ]
    res = coupled_intensity_fit(
        energies, counts, comps, fwhm_at,
        continuum="step_linear", roi_low_ch=lo, cluster_id="M2c",
    )
    GOLD = {
        1588.2: 16257.55,
        1620.5:  2171.53,
        1630.6:  8078.28,
    }
    for cf in res.components:
        E_key = round(cf.E_keV, 1)
        if E_key not in GOLD:
            continue
        gold = GOLD[E_key]
        # 25 % допуск для слабого Bi-212 1620 (низкое S/N), 10 % для Ac-228
        tol = 0.25 if abs(cf.E_keV - 1620.5) < 0.5 else 0.10
        rel = abs(cf.area - gold) / gold
        assert rel < tol, (
            f"area({cf.E_keV:.2f})={cf.area:.0f} vs gold {gold:.0f}, "
            f"Δ={rel:.1%} > {tol*100:.0f} %"
        )
    print(f"  ✓ test_coupled_fit_M2_areas_within_25pct")


if __name__ == "__main__":
    test_coupled_fit_M2_chi2_and_closure()
    test_coupled_fit_M2_areas_within_25pct()
    print("All M2 coupled-fit tests passed.")
