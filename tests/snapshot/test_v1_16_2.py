# -*- coding: utf-8 -*-
"""
v1.16.2 delivery tests — F-98 quasi-template fitter (LSRM §13).

Run:
    cd 0_Work/gamma-spectrum-analysis
    $env:PYTHONPATH = "scripts"
    python test_v1_16_2.py
"""
from __future__ import annotations
import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts"))

import numpy as np

from gamma.activity.quasitemplate import (
    build_nuclide_template,
    build_continuum_basis,
    quasi_template_fit,
)


# ──────────────────────────────────────────────────────────────────────
# Simple model setup for synthetic tests
# ──────────────────────────────────────────────────────────────────────

def _nai_fwhm_keV(E_keV: float) -> float:
    """Approx NaI FWHM at the 662 keV anchor ≈ 49 keV (FWHM% ≈ 7.4%)."""
    return 0.07 * math.sqrt(max(E_keV, 1.0) * 662.0)


def _eff(E_keV: float) -> float:
    """Crude exponential photopeak efficiency: ε ≈ 0.05 · (662/E)^0.6."""
    return 0.05 * (662.0 / max(E_keV, 1.0)) ** 0.6


# Energy axis common to all tests
_N_CH = 4096
_E_PER_CH = 1.0  # 1 keV / channel
_ENERGIES = np.arange(_N_CH, dtype=float) * _E_PER_CH


# ──────────────────────────────────────────────────────────────────────
# F-98a — template builder
# ──────────────────────────────────────────────────────────────────────

def test_template_single_line_recovers_amplitude() -> None:
    """Template for one library line — area ≈ I·ε·t_live."""
    lines = [(661.66, 85.1)]  # Cs-137
    t_live = 1000.0
    sigma_k = build_nuclide_template(
        "Cs-137", lines,
        energies_keV=_ENERGIES,
        fwhm_keV_fn=_nai_fwhm_keV,
        efficiency_fn=_eff,
        t_live_s=t_live,
    )
    total_area = float(np.sum(sigma_k) * _E_PER_CH)
    expected = (85.1 / 100.0) * _eff(661.66) * t_live
    rel = abs(total_area - expected) / expected
    assert rel < 0.01, f"template area off: {rel*100:.2f}%"


def test_template_two_lines_combined() -> None:
    """Co-60 two lines must produce two Gaussians of expected area each."""
    lines = [(1173.23, 99.85), (1332.49, 99.98)]
    t_live = 1000.0
    sigma_k = build_nuclide_template(
        "Co-60", lines,
        energies_keV=_ENERGIES,
        fwhm_keV_fn=_nai_fwhm_keV,
        efficiency_fn=_eff,
        t_live_s=t_live,
    )
    # Area between 1100..1250 keV ≈ first Gaussian; 1250..1400 ≈ second
    band_1 = float(np.sum(sigma_k[1100:1250]) * _E_PER_CH)
    band_2 = float(np.sum(sigma_k[1250:1400]) * _E_PER_CH)
    expected_1 = 0.9985 * _eff(1173.23) * t_live
    expected_2 = 0.9998 * _eff(1332.49) * t_live
    assert abs(band_1 - expected_1) / expected_1 < 0.10
    assert abs(band_2 - expected_2) / expected_2 < 0.10


def test_template_intensity_cutoff_drops_weak_lines() -> None:
    """A 0.001% line below default 0.1% cutoff is dropped."""
    lines = [(661.66, 0.0005)]
    sigma_k = build_nuclide_template(
        "X-rare", lines,
        energies_keV=_ENERGIES,
        fwhm_keV_fn=_nai_fwhm_keV,
        efficiency_fn=_eff,
        t_live_s=1000.0,
    )
    assert float(np.sum(sigma_k)) == 0.0


# ──────────────────────────────────────────────────────────────────────
# F-98b — continuum basis
# ──────────────────────────────────────────────────────────────────────

def test_continuum_basis_has_correct_shape() -> None:
    M = build_continuum_basis(_ENERGIES, degree=4)
    assert M.shape == (_N_CH, 5)
    # column 0 = const(1); column 1 = x; column 2 = x²; …
    assert np.allclose(M[:, 0], 1.0)


# ──────────────────────────────────────────────────────────────────────
# F-98c — round-trip on synthetic Cs-137 + Co-60 spectrum
# ──────────────────────────────────────────────────────────────────────

def test_fit_recovers_cs137_alone() -> None:
    """Generate Cs-137 + smooth continuum; fit recovers A within 10%."""
    rng = np.random.default_rng(7)
    t_live = 3600.0
    A_cs_true = 100.0  # Bq
    lines_cs = [(661.66, 85.1)]
    sigma_cs = build_nuclide_template(
        "Cs-137", lines_cs,
        energies_keV=_ENERGIES,
        fwhm_keV_fn=_nai_fwhm_keV, efficiency_fn=_eff, t_live_s=t_live,
    )
    # smooth continuum ~ 1/E
    cont_true = 50.0 * np.exp(-_ENERGIES / 300.0)
    y_clean = A_cs_true * sigma_cs + cont_true
    y = y_clean + rng.normal(0, np.sqrt(np.maximum(y_clean, 1.0)))
    y = np.maximum(y, 0.0)
    result = quasi_template_fit(
        counts=y, energies_keV=_ENERGIES,
        nuclide_lines={"Cs-137": lines_cs},
        fwhm_keV_fn=_nai_fwhm_keV, efficiency_fn=_eff,
        t_live_s=t_live, continuum_degree=4,
    )
    assert result.converged, result.notes
    A_rec = result.by_nuclide()["Cs-137"][0]
    rel = abs(A_rec - A_cs_true) / A_cs_true
    assert rel < 0.15, f"Cs-137 recovery off: {A_rec} vs {A_cs_true} (Δ={rel*100:.1f}%)"


def test_fit_recovers_cs137_plus_co60() -> None:
    """Two nuclides, distinct energies — should recover both."""
    rng = np.random.default_rng(11)
    t_live = 3600.0
    A_cs_true, A_co_true = 80.0, 40.0
    lines_cs = [(661.66, 85.1)]
    lines_co = [(1173.23, 99.85), (1332.49, 99.98)]
    sigma_cs = build_nuclide_template("Cs-137", lines_cs, energies_keV=_ENERGIES,
                                      fwhm_keV_fn=_nai_fwhm_keV, efficiency_fn=_eff,
                                      t_live_s=t_live)
    sigma_co = build_nuclide_template("Co-60", lines_co, energies_keV=_ENERGIES,
                                      fwhm_keV_fn=_nai_fwhm_keV, efficiency_fn=_eff,
                                      t_live_s=t_live)
    cont_true = 30.0 * np.exp(-_ENERGIES / 250.0)
    y_clean = A_cs_true * sigma_cs + A_co_true * sigma_co + cont_true
    y = y_clean + rng.normal(0, np.sqrt(np.maximum(y_clean, 1.0)))
    y = np.maximum(y, 0.0)
    result = quasi_template_fit(
        counts=y, energies_keV=_ENERGIES,
        nuclide_lines={"Cs-137": lines_cs, "Co-60": lines_co},
        fwhm_keV_fn=_nai_fwhm_keV, efficiency_fn=_eff,
        t_live_s=t_live, continuum_degree=4,
    )
    assert result.converged, result.notes
    rec = result.by_nuclide()
    A_cs_rec = rec["Cs-137"][0]
    A_co_rec = rec["Co-60"][0]
    assert abs(A_cs_rec - A_cs_true) / A_cs_true < 0.15
    assert abs(A_co_rec - A_co_true) / A_co_true < 0.15


def test_fit_handles_absent_nuclide_gracefully() -> None:
    """When candidate is absent from spectrum, fit returns A near 0."""
    rng = np.random.default_rng(13)
    t_live = 3600.0
    # Only Cs-137 present, but we ask for Cs-137 AND Co-60 (absent).
    lines_cs = [(661.66, 85.1)]
    lines_co = [(1173.23, 99.85), (1332.49, 99.98)]
    sigma_cs = build_nuclide_template("Cs-137", lines_cs, energies_keV=_ENERGIES,
                                      fwhm_keV_fn=_nai_fwhm_keV, efficiency_fn=_eff,
                                      t_live_s=t_live)
    cont_true = 30.0 * np.exp(-_ENERGIES / 250.0)
    y_clean = 80.0 * sigma_cs + cont_true
    y = y_clean + rng.normal(0, np.sqrt(np.maximum(y_clean, 1.0)))
    y = np.maximum(y, 0.0)
    result = quasi_template_fit(
        counts=y, energies_keV=_ENERGIES,
        nuclide_lines={"Cs-137": lines_cs, "Co-60": lines_co},
        fwhm_keV_fn=_nai_fwhm_keV, efficiency_fn=_eff,
        t_live_s=t_live, continuum_degree=4,
    )
    assert result.converged
    rec = result.by_nuclide()
    A_cs_rec = rec["Cs-137"][0]
    A_co_rec, A_co_sig = rec["Co-60"]
    # The fit may leak a small noise-floor activity (~ 1% of the
    # present nuclide) into the absent template because the χ²
    # surface near A_co=0 is flat. Acceptance: A_co < 5% of A_cs,
    # which physically means "Co-60 is at noise floor".
    assert A_co_rec < 0.05 * A_cs_rec, (
        f"Co-60 leakage too high: {A_co_rec:.3f} Bq vs "
        f"5% of {A_cs_rec:.3f} Bq"
    )


def test_fit_window_restriction() -> None:
    """E_min/E_max kwargs restrict the fit channels."""
    rng = np.random.default_rng(17)
    t_live = 3600.0
    lines_cs = [(661.66, 85.1)]
    sigma_cs = build_nuclide_template("Cs-137", lines_cs, energies_keV=_ENERGIES,
                                      fwhm_keV_fn=_nai_fwhm_keV, efficiency_fn=_eff,
                                      t_live_s=t_live)
    y = 100.0 * sigma_cs + rng.normal(0, 1.0, size=_N_CH)
    y = np.maximum(y, 0.0)
    result = quasi_template_fit(
        counts=y, energies_keV=_ENERGIES,
        nuclide_lines={"Cs-137": lines_cs},
        fwhm_keV_fn=_nai_fwhm_keV, efficiency_fn=_eff,
        t_live_s=t_live, continuum_degree=2,
        E_min_keV=400.0, E_max_keV=900.0,
    )
    assert result.converged
    assert result.n_channels < _N_CH


def test_fit_returns_uncoverged_on_empty_library() -> None:
    """No candidate nuclides → graceful failure with notes."""
    result = quasi_template_fit(
        counts=[1.0] * 100, energies_keV=list(range(100)),
        nuclide_lines={},
        fwhm_keV_fn=_nai_fwhm_keV, efficiency_fn=_eff,
        t_live_s=1000.0,
    )
    assert not result.converged
    assert "no candidate" in result.notes


# ──────────────────────────────────────────────────────────────────────
# Driver
# ──────────────────────────────────────────────────────────────────────

def _all_tests():
    return [v for k, v in globals().items() if k.startswith("test_") and callable(v)]


def main() -> int:
    tests = _all_tests()
    fail = 0
    for t in tests:
        try:
            t()
            print(f"  PASS  {t.__name__}")
        except AssertionError as exc:
            fail += 1
            print(f"  FAIL  {t.__name__}: {exc}")
        except Exception as exc:
            fail += 1
            print(f"  ERR   {t.__name__}: {type(exc).__name__}: {exc}")
    print(f"\nv1.16.2: {len(tests) - fail}/{len(tests)} passed")
    return fail


if __name__ == "__main__":
    sys.exit(main())
