"""
Tests for Phase 2.1c efficiency calibration:
  - Efficiency curve fitting (log-log polynomial)
  - Activity-ratio (efficiency-corrected) proportionality check
  - Validates K-12 closure for close-by line nuclides
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))


@pytest.fixture(scope="module", autouse=True)
def _restore_nuclide_library_after_module():
    """BUG-21 / v1.18.32 — TD-2 test isolation fix (module-scoped).

    ``test_check_intensity_proportionality_co60`` (and friends) load the
    Gamma-1S library with ``merge_mode="override"`` into the
    module-global ``_CACHE`` in ``gamma.data.nuclide_library``. Without
    teardown, downstream tests (notably
    ``tests/snapshot/test_f389_v2_activity_parity.py``) inherit the
    override and drift — see BUG-21 / TD-2 commentary in
    ``tests/step08_multiplets/test_deconvolve.py`` for full mechanism.
    """
    yield
    from gamma.data.nuclide_library import reset_cache
    reset_cache()

from gamma.calibration.efficiency import (
    EfficiencyCurve, fit_efficiency_curve, fit_efficiency_from_efr_file,
    REFERENCE_ENERGY_KEV,
)
from gamma.identification import check_intensity_proportionality


EFR_PATH = ("detectors/Gamma-1S/efficiency/Gamma-1S_NaI_63x63_USB_SN-01/"
            "УДС-ГЦ-63х63-USB__SN-01_-_Точечная-5см.efr")


# --- Efficiency curve fitting ---

def test_efficiency_fit_synthetic():
    """Fit a known ε(E) curve to synthetic data — should recover the
    underlying coefficients."""
    # Generate data with known coefficients a0+a1*ln(E)+a2*ln(E)²
    # Pick something physical: rising at low E, falling at high E
    a0_true = 1.0
    a1_true = -0.5
    a2_true = -0.1
    E_vals = [60, 100, 200, 400, 700, 1000, 1500, 2000, 2600]
    eps_vals = [math.exp(a0_true + a1_true*math.log(E) + a2_true*math.log(E)**2)
                for E in E_vals]
    dpct = [2.0] * len(E_vals)
    curve = fit_efficiency_curve(E_vals, eps_vals, dpct, degree=2)
    # Coefficients should match within tight tolerance
    assert abs(curve.coefficients[0] - a0_true) < 0.05, f"a0={curve.coefficients[0]} vs {a0_true}"
    assert abs(curve.coefficients[1] - a1_true) < 0.05
    assert abs(curve.coefficients[2] - a2_true) < 0.05
    # χ²/dof should be very small (noiseless data)
    assert curve.chi2_per_dof < 0.001
    print(f"  ✓ test_efficiency_fit_synthetic (coefs recovered to <5%)")


def test_efficiency_fit_real_5cm():
    """Fit the Gamma-1S 5cm point source efficiency from .efr file."""
    curve = fit_efficiency_from_efr_file(EFR_PATH, degree=3)
    # Range check
    assert 50 < curve.E_min_keV < 80, f"E_min={curve.E_min_keV}"
    assert curve.E_max_keV > 2000, f"E_max={curve.E_max_keV}"
    # Should have 24 data points
    assert curve.n_points_used == 24
    # Cs-137 661.66 keV: library says ε=1.3742e-02 — fit should predict this
    # very accurately (it's one of the data points)
    eps_661 = curve(661.66)
    assert abs(eps_661 - 1.3742e-02) / 1.3742e-02 < 0.10, \
        f"ε(661) predicted {eps_661}, library 1.3742e-02"
    # Curve should be monotonically decreasing above ~100 keV (NaI behaviour)
    eps_200 = curve(200)
    eps_500 = curve(500)
    eps_1000 = curve(1000)
    eps_2000 = curve(2000)
    assert eps_200 > eps_500 > eps_1000 > eps_2000, \
        f"Monotonic decrease violated: {eps_200}, {eps_500}, {eps_1000}, {eps_2000}"
    print(f"  ✓ test_efficiency_fit_real_5cm "
          f"(ε(661)={eps_661:.4e}, χ²/dof={curve.chi2_per_dof:.2f})")


def test_efficiency_extrapolation_flag():
    """Curve should report when called outside its calibrated range."""
    curve = fit_efficiency_from_efr_file(EFR_PATH, degree=3)
    assert not curve.is_extrapolating(661.66)  # Cs-137 — well inside
    assert curve.is_extrapolating(10.0)        # below E_min
    assert curve.is_extrapolating(5000.0)      # above E_max
    print(f"  ✓ test_efficiency_extrapolation_flag")


def test_efficiency_fit_input_validation():
    """fit_efficiency_curve should validate inputs."""
    # Mismatched array lengths
    try:
        fit_efficiency_curve([100, 200], [0.01], [2.0])
        assert False, "Should have raised ValueError"
    except ValueError:
        pass
    # Too few points for degree
    try:
        fit_efficiency_curve([100, 200], [0.01, 0.005], [2.0, 2.0], degree=3)
        assert False, "Should have raised ValueError"
    except ValueError:
        pass
    print(f"  ✓ test_efficiency_fit_input_validation")


def test_efficiency_curve_repr():
    """EfficiencyCurve __repr__ should be informative."""
    curve = fit_efficiency_from_efr_file(EFR_PATH, degree=3)
    s = repr(curve)
    assert "EfficiencyCurve" in s
    assert "degree" in s
    assert "χ²" in s or "chi2" in s.lower()
    print(f"  ✓ test_efficiency_curve_repr: {s}")


# --- Proportionality with efficiency correction ---

def make_match(E, I_pct, area, sigma=10.0):
    """Build a synthetic LineMatch."""
    from gamma.identification.identify import LineMatch
    return LineMatch(
        nuclide="X",
        library_E_keV=E,
        library_I_pct=I_pct,
        peak_channel=int(E * 5),
        peak_E_keV=E,
        peak_sigma=sigma,
        residual_keV=0.5,
        is_characteristic=False,
        peak_area=area,
    )


def test_proportionality_efficiency_corrected_passes():
    """When efficiency curve is provided, activity ratios should equal
    1.0 for a single nuclide."""
    # Build a fake curve: ε(E) = 0.1 * exp(-E/500)
    def fake_eff(E):
        return 0.1 * math.exp(-E / 500.0)

    # Construct matches where (area / ε / I) is the same constant for all
    # Activity = 1000 Bq (notional), so area = activity * ε * I
    activity = 1000.0
    matches = [
        make_match(200, 50.0, activity * fake_eff(200) * 0.5),
        make_match(500, 30.0, activity * fake_eff(500) * 0.3),
        make_match(1500, 40.0, activity * fake_eff(1500) * 0.4),
    ]
    result = check_intensity_proportionality(
        "X", matches, efficiency_curve=fake_eff,
    )
    assert result.passed, f"Should pass with consistent activity: {result.reason}"
    assert "activity ratio" in result.reason
    print(f"  ✓ test_proportionality_efficiency_corrected_passes")


def test_proportionality_efficiency_corrected_catches_wide_separation():
    """Without efficiency correction, widely-separated lines fail.
    WITH correction, they pass."""
    # ε(E) varies by factor 10 between 200 and 2000 keV
    def fake_eff(E):
        return 0.1 * math.exp(-E / 500.0)

    activity = 1000.0
    # Two lines: 200 keV and 2000 keV — factor 10 in energy.
    # Without ε correction: area ratio ≠ I ratio (off by ε ratio)
    matches = [
        make_match(200, 50.0, activity * fake_eff(200) * 0.5),
        make_match(2000, 50.0, activity * fake_eff(2000) * 0.5),
    ]
    # With efficiency: activity ratio = 1.0 → PASS
    result_with_eff = check_intensity_proportionality(
        "X", matches, efficiency_curve=fake_eff,
    )
    assert result_with_eff.passed
    # Without efficiency: area ratio is ε(200)/ε(2000) ≈ exp(1800/500) ≈ 36×
    # The wide tolerance is 5×, so this FAILS
    result_no_eff = check_intensity_proportionality(
        "X", matches, efficiency_curve=None,
    )
    assert not result_no_eff.passed, \
        f"Without ε should fail for wide-separation: {result_no_eff.reason}"
    print(f"  ✓ test_proportionality_efficiency_corrected_catches_wide_separation "
          f"(K-12 closure verified)")


def test_proportionality_co60_real_data():
    """On real Co-60 spectrum (5cm geometry, same detector as efficiency),
    activity ratio of 1173 vs 1332 should be ~1.0 within counting stats."""
    from gamma.io.readers import read_spectrum
    from gamma.peaks.search import mariscotti_search
    from gamma.calibration.fwhm_provider import make_fwhm_at_channel_provider
    from gamma.data.nuclide_library import (
        load_external_library, reset_cache,
    )
    from gamma.identification import (
        build_identification_window, identify_nuclides,
        disambiguate_identifications,
    )
    reset_cache()
    load_external_library(
        "references/nuclide_libraries/Gamma-1S_NaI_63x63_USB_SN-01_lsrm_v2.lib",
        merge_mode="override", split_chains=True,
    )
    eff = fit_efficiency_from_efr_file(EFR_PATH, degree=3)
    spec = read_spectrum(
        "detectors/Gamma-1S/reference_spectra/archive/"
        "Co-60__043_02_2019_Точечная-5см_5cm.spe"
    )
    fwhm_at = make_fwhm_at_channel_provider(spec, fallback_channels=15.0)
    peaks = mariscotti_search(spec.counts, fwhm_channels=fwhm_at, sigma_threshold=5.0)
    window = build_identification_window("NaI", delta_E0_keV=15.0)
    raw = identify_nuclides(found_peaks=peaks, spec=spec, window=window,
                             fwhm_at_channel=fwhm_at)
    refined = disambiguate_identifications(raw)
    co60 = next((ni for ni in refined.detected_nuclides if ni.nuclide == "Co-60"), None)
    assert co60 is not None, "Co-60 should be detected"
    matches = list(co60.matched_lines)
    prop = check_intensity_proportionality(
        "Co-60", matches, efficiency_curve=eff,
    )
    assert prop.passed, f"Co-60 activity ratio should pass: {prop.reason}"
    print(f"  ✓ test_proportionality_co60_real_data: "
          f"{prop.n_ratios_passed}/{prop.n_ratios_passed+prop.n_ratios_failed} pairs PASS")


if __name__ == "__main__":
    print("Running Phase 2.1c efficiency calibration tests...\n")
    test_efficiency_fit_synthetic()
    test_efficiency_fit_real_5cm()
    test_efficiency_extrapolation_flag()
    test_efficiency_fit_input_validation()
    test_efficiency_curve_repr()
    test_proportionality_efficiency_corrected_passes()
    test_proportionality_efficiency_corrected_catches_wide_separation()
    test_proportionality_co60_real_data()
    print("\n✓ All Phase 2.1c efficiency calibration tests passed.")
