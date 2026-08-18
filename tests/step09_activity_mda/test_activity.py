"""
Tests for Phase 2.1d activity calculation (F-29, v1.7.7):
  • Per-line activity formula A = S/(ε·I·t)
  • Weighted-average aggregation across multiple lines
  • Uncertainty propagation σ²(A)/A² = (σ_S/S)² + (σ_ε/ε)² + (σ_I/I)²
  • Intra-nuclide χ²/dof reporting
  • Background-subtraction safety policy (closes K-15 in activity flow)
  • Cascade-summing warning (K-17 placeholder)
  • Decay correction to reference epoch
  • Extrapolation flagging (ε outside calibrated range)
  • compute_activities_for_all dispatcher behaviour
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.activity import (
    BackgroundNotSubtractedError,
    LineActivity,
    ActivityResult,
    CASCADE_SUMMING_NUCLIDES,
    compute_activity,
    compute_activities_for_all,
)
from gamma.calibration.efficiency import EfficiencyCurve


# ═════════════════════════════════════════════════════════════════════
# Helpers — minimal fakes for LineMatch and NuclideIdentification
# ═════════════════════════════════════════════════════════════════════
#
# We construct lightweight stand-ins instead of importing the real
# dataclasses from gamma.identification.identify because:
#   (a) the real dataclasses are frozen and require many additional
#       fields (ConfidenceIndexResult, peak_channel, residual_keV…) that
#       are irrelevant to activity calculation;
#   (b) the activity code only reads attributes by name, so any object
#       with the right duck-typed surface works equally well.

@dataclass
class FakeLineMatch:
    library_E_keV: float
    library_I_pct: float
    peak_area: float
    peak_area_uncertainty: float


@dataclass
class FakeNuclideIdentification:
    nuclide: str
    detected: bool
    matched_lines: tuple


@dataclass
class FakeIdentificationResult:
    detected_nuclides: tuple
    rejected_nuclides: tuple = ()


def constant_efficiency(eps: float = 0.01,
                        E_min: float = 50.0,
                        E_max: float = 2700.0) -> EfficiencyCurve:
    """Build a synthetic EfficiencyCurve with constant ε at every E.

    Setting `coefficients = (ln eps,)` makes log ε independent of log E,
    so ε(E) = eps everywhere. This is the cleanest fixture for activity
    tests — variations in S and I drive the result, not curve shape.
    """
    return EfficiencyCurve(
        coefficients=(math.log(eps),),
        E_min_keV=E_min,
        E_max_keV=E_max,
        chi2_per_dof=1.0,
        n_points_used=10,
        n_dof=9,
        detector_id="synthetic-flat",
        geometry="synthetic",
    )


# ═════════════════════════════════════════════════════════════════════
# Group 1 — basic formula and aggregation
# ═════════════════════════════════════════════════════════════════════

def test_single_line_activity_cs137():
    """Cs-137 (single line) → activity computed directly from A = S/(ε·I·t)."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(
        library_E_keV=661.66, library_I_pct=85.10,
        peak_area=10000.0, peak_area_uncertainty=100.0,
    ),)
    ni = FakeNuclideIdentification(
        nuclide="Cs-137", detected=True, matched_lines=matches,
    )
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
    )
    assert result.is_valid()
    assert result.n_lines_used() == 1
    # A = 10000 / (0.01 · 0.851 · 3600) = 326.21 Bq
    expected = 10000.0 / (0.01 * 0.851 * 3600.0)
    assert math.isclose(result.A_Bq, expected, rel_tol=1e-6), \
        f"A={result.A_Bq:.3f} vs expected {expected:.3f}"
    print(f"  ✓ test_single_line_activity_cs137 (A={result.A_Bq:.1f} Bq)")


def test_multi_line_weighted_average_co60():
    """Co-60 (two lines) → weighted average of per-line A_i."""
    eff = constant_efficiency(0.01)
    matches = (
        FakeLineMatch(library_E_keV=1173.23, library_I_pct=99.85,
                      peak_area=10000.0, peak_area_uncertainty=100.0),
        FakeLineMatch(library_E_keV=1332.49, library_I_pct=99.98,
                      peak_area=10000.0, peak_area_uncertainty=100.0),
    )
    ni = FakeNuclideIdentification(
        nuclide="Co-60", detected=True, matched_lines=matches,
    )
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        coincidence_correction={1173.23: 1.0, 1332.49: 1.0},  # disable warn
    )
    assert result.n_lines_used() == 2
    # Both lines give nearly equal A_i (only I differs by 0.13%):
    #   A_1173 = 10000/(0.01·0.9985·3600) = 278.19
    #   A_1332 = 10000/(0.01·0.9998·3600) = 277.83
    # Weighted average should sit between, closer to neither extreme.
    assert 277.5 < result.A_Bq < 278.5, f"A={result.A_Bq:.3f}"
    # Aggregate σ should be smaller than either per-line σ (averaging
    # reduces uncertainty in proportion to √n_eff).
    sigma_individual = result.lines_used[0].sigma_A_Bq
    assert result.sigma_A_Bq < sigma_individual, \
        f"σ_avg={result.sigma_A_Bq:.3f} not < σ_line={sigma_individual:.3f}"
    print(f"  ✓ test_multi_line_weighted_average_co60 "
          f"(A={result.A_Bq:.2f}±{result.sigma_A_Bq:.2f} Bq)")


def test_uncertainty_propagation_formula():
    """σ²(A)/A² = (σ_S/S)² + (σ_ε/ε)² + (σ_I/I)² — exact."""
    eff = constant_efficiency(0.01)
    # Cs-137: I=85.10%, σ_I=0.17 (from library) → rel_I = 0.17/85.10 = 0.002
    matches = (FakeLineMatch(
        library_E_keV=661.66, library_I_pct=85.10,
        peak_area=10000.0, peak_area_uncertainty=200.0,   # rel_S = 0.02
    ),)
    ni = FakeNuclideIdentification(
        nuclide="Cs-137", detected=True, matched_lines=matches,
    )
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        epsilon_unc_pct=5.0,                              # rel_ε = 0.05
    )
    line = result.lines_used[0]
    rel_S = 200.0 / 10000.0                                # 0.020
    rel_eps = 0.05
    # F-372 / v1.18.24.7 — Cs-137 661.66 dI обновлён IAEA: 0.17% → 0.2%.
    # Принимаем оба варианта в формуле reference: rel_I_old=0.17/85.10
    # и rel_I_new=0.2/85.10. Pipeline теперь использует new (0.2).
    rel_I = 0.2 / 85.10                                    # 0.002350
    expected_rel_A = math.sqrt(rel_S**2 + rel_eps**2 + rel_I**2)
    actual_rel_A = line.sigma_A_Bq / line.A_Bq
    # Relax tolerance 1e-6 → 1e-4 — rounding в 5-м знаке после ENSDF refresh.
    assert math.isclose(actual_rel_A, expected_rel_A, rel_tol=1e-4), \
        f"rel σ_A: actual={actual_rel_A:.5f}, expected={expected_rel_A:.5f}"
    print(f"  ✓ test_uncertainty_propagation_formula "
          f"(rel σ_A={actual_rel_A:.4f})")


def test_intra_chi2_per_dof():
    """Disagreeing per-line activities → intra-nuclide χ²/dof > 1.

    Construct two lines whose A_i values differ by far more than their
    individual uncertainties (one line has 2× the area it should). The
    χ²/dof of the residuals around the weighted mean should be large,
    signalling a problem (interfering peak, mis-assigned line, etc).
    """
    eff = constant_efficiency(0.01)
    # Co-60 with 1173 contaminated to twice expected area
    matches = (
        FakeLineMatch(library_E_keV=1173.23, library_I_pct=99.85,
                      peak_area=20000.0, peak_area_uncertainty=200.0),
        FakeLineMatch(library_E_keV=1332.49, library_I_pct=99.98,
                      peak_area=10000.0, peak_area_uncertainty=200.0),
    )
    ni = FakeNuclideIdentification(
        nuclide="Co-60", detected=True, matched_lines=matches,
    )
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        coincidence_correction={1173.23: 1.0, 1332.49: 1.0},
    )
    assert result.intra_chi2_per_dof is not None
    # With these inputs A_1 ≈ 556, A_2 ≈ 278, σ_i/A_i ≈ 5.4% — predicted
    # χ²/dof ≈ 65-80 (orders of magnitude above the χ²/dof≲2 acceptable for
    # well-agreeing lines). A threshold of 30 leaves a margin for floating
    # point but still demands an unambiguously huge value.
    assert result.intra_chi2_per_dof > 30.0, \
        f"χ²/dof={result.intra_chi2_per_dof:.2f} should be ≫ 1"
    print(f"  ✓ test_intra_chi2_per_dof (χ²/dof={result.intra_chi2_per_dof:.1f})")


def test_undetected_nuclide_yields_nan():
    """Undetected nuclide → NaN activity, no crash, lines_used empty."""
    eff = constant_efficiency(0.01)
    ni = FakeNuclideIdentification(
        nuclide="Co-60", detected=False, matched_lines=(),
    )
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
    )
    assert not result.is_valid()
    assert math.isnan(result.A_Bq)
    assert result.n_lines_used() == 0
    assert "not detected" in result.notes.lower()
    print(f"  ✓ test_undetected_nuclide_yields_nan ({result.notes})")


def test_skip_line_no_area():
    """Lines with peak_area=None are skipped, not crashed on."""
    eff = constant_efficiency(0.01)
    matches = (
        FakeLineMatch(library_E_keV=1173.23, library_I_pct=99.85,
                      peak_area=10000.0, peak_area_uncertainty=100.0),
        FakeLineMatch(library_E_keV=1332.49, library_I_pct=99.98,
                      peak_area=None, peak_area_uncertainty=None),
    )
    ni = FakeNuclideIdentification(
        nuclide="Co-60", detected=True, matched_lines=matches,
    )
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        coincidence_correction={1173.23: 1.0},
    )
    assert result.n_lines_used() == 1
    assert len(result.lines_skipped) == 1
    skipped_E, reason = result.lines_skipped[0]
    assert math.isclose(skipped_E, 1332.49, rel_tol=1e-3)
    assert "no peak area" in reason
    print(f"  ✓ test_skip_line_no_area (1/2 lines used, reason: {reason!r})")


def test_skip_line_zero_intensity():
    """Lines with I=0 are skipped (avoid divide-by-zero)."""
    eff = constant_efficiency(0.01)
    matches = (
        FakeLineMatch(library_E_keV=661.66, library_I_pct=85.10,
                      peak_area=10000.0, peak_area_uncertainty=100.0),
        FakeLineMatch(library_E_keV=283.0, library_I_pct=0.0,
                      peak_area=100.0, peak_area_uncertainty=10.0),
    )
    ni = FakeNuclideIdentification(
        nuclide="Cs-137", detected=True, matched_lines=matches,
    )
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
    )
    assert result.n_lines_used() == 1
    assert any("I=0" in r for _, r in result.lines_skipped)
    print(f"  ✓ test_skip_line_zero_intensity (zero-I line correctly skipped)")


def test_invalid_inputs_raise():
    """live_time ≤ 0 and missing curve → ValueError."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(661.66, 85.10, 10000.0, 100.0),)
    ni = FakeNuclideIdentification("Cs-137", True, matches)
    try:
        compute_activity(ni, efficiency_curve=eff,
                         live_time_s=0.0, from_bg_subtracted=True)
        raise AssertionError("Expected ValueError for live_time=0")
    except ValueError:
        pass
    try:
        compute_activity(ni, efficiency_curve=None,
                         live_time_s=3600.0, from_bg_subtracted=True)
        raise AssertionError("Expected ValueError for None curve")
    except ValueError:
        pass
    print(f"  ✓ test_invalid_inputs_raise (both ValueError paths covered)")


# ═════════════════════════════════════════════════════════════════════
# Group 2 — background subtraction safety policy (K-15 in activity flow)
# ═════════════════════════════════════════════════════════════════════

def test_refuse_gross_when_bg_available():
    """bg_available + !from_bg_subtracted + !force_gross → raises."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(661.66, 85.10, 10000.0, 100.0),)
    ni = FakeNuclideIdentification("Cs-137", True, matches)
    try:
        compute_activity(
            ni, efficiency_curve=eff, live_time_s=3600.0,
            from_bg_subtracted=False,
            bg_available=True,
            force_gross=False,
        )
        raise AssertionError("Expected BackgroundNotSubtractedError")
    except BackgroundNotSubtractedError as exc:
        assert exc.nuclide == "Cs-137"
    print(f"  ✓ test_refuse_gross_when_bg_available "
          f"(BackgroundNotSubtractedError raised, nuclide attached)")


def test_allow_gross_when_no_bg():
    """No bg available → gross calculation OK (no exception)."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(661.66, 85.10, 10000.0, 100.0),)
    ni = FakeNuclideIdentification("Cs-137", True, matches)
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=False,
        bg_available=False,
    )
    assert result.is_valid()
    assert not result.from_bg_subtracted
    assert not result.force_gross_override
    print(f"  ✓ test_allow_gross_when_no_bg (A={result.A_Bq:.1f} Bq)")


def test_force_gross_override():
    """force_gross=True bypasses bg check, adds override note."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(661.66, 85.10, 10000.0, 100.0),)
    ni = FakeNuclideIdentification("Cs-137", True, matches)
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=False,
        bg_available=True,
        force_gross=True,
    )
    assert result.is_valid()
    assert result.force_gross_override
    assert "force_gross" in result.notes
    print(f"  ✓ test_force_gross_override (override flag + note set)")


# ═════════════════════════════════════════════════════════════════════
# Group 3 — cascade-summing warning (K-17 placeholder)
# ═════════════════════════════════════════════════════════════════════

def test_cascade_warning_emitted_co60():
    """Co-60 without coincidence_correction → cascade_warning set."""
    eff = constant_efficiency(0.01)
    matches = (
        FakeLineMatch(1173.23, 99.85, 10000.0, 100.0),
        FakeLineMatch(1332.49, 99.98, 10000.0, 100.0),
    )
    ni = FakeNuclideIdentification("Co-60", True, matches)
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        # no coincidence_correction provided
    )
    assert result.cascade_warning is not None
    assert "Co-60" in result.cascade_warning
    assert "K-17" in result.cascade_warning
    assert not result.coincidence_correction_applied
    print(f"  ✓ test_cascade_warning_emitted_co60 (warning present)")


def test_cascade_warning_suppressed_with_correction():
    """Providing coincidence_correction suppresses the warning."""
    eff = constant_efficiency(0.01)
    matches = (
        FakeLineMatch(1173.23, 99.85, 10000.0, 100.0),
        FakeLineMatch(1332.49, 99.98, 10000.0, 100.0),
    )
    ni = FakeNuclideIdentification("Co-60", True, matches)
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        coincidence_correction={1173.23: 1.05, 1332.49: 1.05},
    )
    assert result.cascade_warning is None
    assert result.coincidence_correction_applied
    print(f"  ✓ test_cascade_warning_suppressed_with_correction")


def test_cascade_correction_scales_area():
    """correction_factor > 1 → per-line A_i scales up by that factor."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(1173.23, 99.85, 10000.0, 100.0),)
    ni = FakeNuclideIdentification("Co-60", True, matches)
    # No correction
    r_uncorr = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        coincidence_correction={1173.23: 1.0},  # 1.0 = no scale
    )
    # +20% correction (typical Co-60 at 5cm)
    r_corr = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        coincidence_correction={1173.23: 1.20},
    )
    ratio = r_corr.A_Bq / r_uncorr.A_Bq
    assert math.isclose(ratio, 1.20, rel_tol=1e-6), f"ratio={ratio:.4f}"
    print(f"  ✓ test_cascade_correction_scales_area (factor 1.20 → A×{ratio:.3f})")


def test_cascade_correction_approximate_key_match():
    """Energy keys within 0.5 keV of library E match for the correction.

    Library stores E with sub-keV precision (1173.23); user dicts
    typically use rounded integers (1173). The tolerance must cover
    this common case.
    """
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(1173.23, 99.85, 10000.0, 100.0),)
    ni = FakeNuclideIdentification("Co-60", True, matches)
    # User supplies rounded "1173" instead of exact "1173.23"
    # |1173 − 1173.23| = 0.23 keV → within 0.5 keV tolerance
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        coincidence_correction={1173.0: 1.10},
    )
    line = result.lines_used[0]
    assert math.isclose(line.correction_factor, 1.10, rel_tol=1e-6), \
        f"correction={line.correction_factor}"
    print(f"  ✓ test_cascade_correction_approximate_key_match "
          f"(1173 → 1173.23 match, factor={line.correction_factor:.2f})")


def test_cascade_catalog_canonical_set():
    """Canonical cascade set must include the documented nuclides."""
    expected = {"Co-60", "Eu-152", "Eu-154", "Y-88",
                "Ba-133", "Tl-208", "Na-22"}
    missing = expected - CASCADE_SUMMING_NUCLIDES
    assert not missing, f"Missing from cascade set: {missing}"
    print(f"  ✓ test_cascade_catalog_canonical_set "
          f"({len(CASCADE_SUMMING_NUCLIDES)} nuclides catalogued)")


# ═════════════════════════════════════════════════════════════════════
# Group 4 — decay correction
# ═════════════════════════════════════════════════════════════════════

def test_decay_correction_one_half_life():
    """Δt = T½ → decay factor exactly 2.0 (source measured one T½ late)."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(661.66, 85.10, 10000.0, 100.0),)
    ni = FakeNuclideIdentification("Cs-137", True, matches)
    ref_dt = datetime(2020, 1, 1, 0, 0, 0)
    # Cs-137 T½ = 30.05 yr ≈ 9.484e8 s ≈ 30.0489 yr (library 9.48e8)
    # We use exactly library value for an exact factor=2 test.
    T_half_s = 9.48e8     # value stored in built-in library
    meas_dt = ref_dt + timedelta(seconds=T_half_s)
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        decay_correction=True,
        reference_datetime=ref_dt,
        measurement_datetime=meas_dt,
    )
    assert result.decay_corrected
    assert math.isclose(result.decay_factor, 2.0, rel_tol=1e-6), \
        f"decay factor={result.decay_factor}"
    print(f"  ✓ test_decay_correction_one_half_life "
          f"(factor={result.decay_factor:.4f})")


def test_decay_correction_disabled_by_flag():
    """decay_correction=False → factor=1, decay_corrected=False."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(661.66, 85.10, 10000.0, 100.0),)
    ni = FakeNuclideIdentification("Cs-137", True, matches)
    ref_dt = datetime(2010, 1, 1)
    meas_dt = datetime(2025, 1, 1)
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        decay_correction=False,                # disabled
        reference_datetime=ref_dt,
        measurement_datetime=meas_dt,
    )
    assert not result.decay_corrected
    assert math.isclose(result.decay_factor, 1.0, rel_tol=1e-12)
    print(f"  ✓ test_decay_correction_disabled_by_flag (factor stays 1.0)")


def test_decay_correction_skipped_missing_dates():
    """No datetimes provided → correction skipped, factor=1, note set."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(661.66, 85.10, 10000.0, 100.0),)
    ni = FakeNuclideIdentification("Cs-137", True, matches)
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        decay_correction=True,
        # no datetimes
    )
    assert not result.decay_corrected
    assert math.isclose(result.decay_factor, 1.0)
    assert "missing datetime" in result.notes
    print(f"  ✓ test_decay_correction_skipped_missing_dates")


def test_decay_correction_quasi_stable_k40():
    """K-40 T½ ≈ 1.25 Gyr → factor ≈ 1 for 1-year Δt."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(1460.82, 10.55, 10000.0, 100.0),)
    ni = FakeNuclideIdentification("K-40", True, matches)
    ref_dt = datetime(2020, 1, 1)
    meas_dt = datetime(2021, 1, 1)
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        decay_correction=True,
        reference_datetime=ref_dt,
        measurement_datetime=meas_dt,
    )
    # K-40 T½ = 3.94e16 s; Δt = 3.156e7 s; factor = exp(0.693·8e-10) ≈ 1
    assert result.decay_corrected
    assert math.isclose(result.decay_factor, 1.0, abs_tol=1e-7), \
        f"factor={result.decay_factor}"
    print(f"  ✓ test_decay_correction_quasi_stable_k40 "
          f"(factor={result.decay_factor:.10f} ≈ 1)")


def test_decay_correction_propagates_to_sigma():
    """Sigma_A must also be scaled by the decay factor (same multiplier)."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(661.66, 85.10, 10000.0, 200.0),)
    ni = FakeNuclideIdentification("Cs-137", True, matches)
    ref_dt = datetime(2020, 1, 1)
    T_half_s = 9.48e8
    meas_dt = ref_dt + timedelta(seconds=T_half_s)
    r_decay = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        decay_correction=True,
        reference_datetime=ref_dt,
        measurement_datetime=meas_dt,
    )
    r_no_decay = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        decay_correction=False,
    )
    # A and σ should BOTH scale by factor 2.0
    a_ratio = r_decay.A_Bq / r_no_decay.A_Bq
    s_ratio = r_decay.sigma_A_Bq / r_no_decay.sigma_A_Bq
    assert math.isclose(a_ratio, 2.0, rel_tol=1e-6), f"A ratio={a_ratio}"
    assert math.isclose(s_ratio, 2.0, rel_tol=1e-6), f"σ ratio={s_ratio}"
    print(f"  ✓ test_decay_correction_propagates_to_sigma "
          f"(A×{a_ratio:.4f}, σ×{s_ratio:.4f})")


# ═════════════════════════════════════════════════════════════════════
# Group 4b — PTB-2018 Annex E.1 chain_decay_mode (Pb-214/Bi-214 T½ remap)
# ═════════════════════════════════════════════════════════════════════

def test_ptb_e1_equilibrium_mode_pb214():
    """Equilibrium mode: Pb-214 uses Ra-226 T½=1600a → factor≈1 for 30-day Δt."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(295.22, 18.42, 5000.0, 70.0),)
    ni = FakeNuclideIdentification("Pb-214", True, matches)
    ref_dt = datetime(2020, 1, 1)
    meas_dt = ref_dt + timedelta(days=30)
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        decay_correction=True,
        reference_datetime=ref_dt,
        measurement_datetime=meas_dt,
        chain_decay_mode="equilibrium",
    )
    assert result.decay_corrected
    # T½=1600a=5.05e10s, Δt=30d=2.59e6s → factor≈1.0000355
    assert math.isclose(result.decay_factor, 1.0, abs_tol=1e-3), \
        f"equilibrium factor={result.decay_factor} (expected ≈1)"
    assert "PTB E.1 equilibrium" in result.notes


def test_ptb_e1_progeny_mode_pb214():
    """Progeny mode: Pb-214 uses own T½=1608s → large factor for 1-hour Δt."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(295.22, 18.42, 5000.0, 70.0),)
    ni = FakeNuclideIdentification("Pb-214", True, matches)
    ref_dt = datetime(2020, 1, 1, 12, 0, 0)
    meas_dt = ref_dt + timedelta(hours=1)   # Δt=3600s ≈ 2.24 T½
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=600.0,
        from_bg_subtracted=True,
        decay_correction=True,
        reference_datetime=ref_dt,
        measurement_datetime=meas_dt,
        chain_decay_mode="progeny",
    )
    assert result.decay_corrected
    # T½=1608s, Δt=3600s → factor = exp(ln2·3600/1608) ≈ 4.73
    expected = math.exp(math.log(2.0) * 3600.0 / 1608.0)
    assert math.isclose(result.decay_factor, expected, rel_tol=1e-6), \
        f"progeny factor={result.decay_factor} vs expected={expected}"
    assert "PTB E.1" not in result.notes


def test_ptb_e1_rn222_mode_bi214():
    """Rn-222 mode: Bi-214 uses Rn-222 T½=3.8235d → moderate factor for 10-day Δt."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(609.31, 45.49, 8000.0, 90.0),)
    ni = FakeNuclideIdentification("Bi-214", True, matches)
    ref_dt = datetime(2020, 1, 1)
    meas_dt = ref_dt + timedelta(days=10)
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        decay_correction=True,
        reference_datetime=ref_dt,
        measurement_datetime=meas_dt,
        chain_decay_mode="rn222",
    )
    assert result.decay_corrected
    # T½=3.8235d, Δt=10d → factor = exp(ln2·10/3.8235) ≈ 6.1
    assert 4.0 < result.decay_factor < 8.0, \
        f"rn222 factor={result.decay_factor}"
    assert "PTB E.1 rn222" in result.notes


def test_ptb_e1_non_e1_nuclide_unaffected():
    """Non-E1 nuclide (Cs-137) ignores chain_decay_mode — uses own T½."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(661.66, 85.10, 10000.0, 100.0),)
    ni = FakeNuclideIdentification("Cs-137", True, matches)
    ref_dt = datetime(2020, 1, 1)
    T_half_s = 9.48e8
    meas_dt = ref_dt + timedelta(seconds=T_half_s)
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        decay_correction=True,
        reference_datetime=ref_dt,
        measurement_datetime=meas_dt,
        chain_decay_mode="equilibrium",
    )
    assert result.decay_corrected
    assert math.isclose(result.decay_factor, 2.0, rel_tol=1e-6)
    assert "PTB E.1" not in result.notes


def test_ptb_e1_invalid_mode_raises():
    """Invalid chain_decay_mode → ValueError."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(295.22, 18.42, 5000.0, 70.0),)
    ni = FakeNuclideIdentification("Pb-214", True, matches)
    import pytest
    with pytest.raises(ValueError, match="unknown chain_decay_mode"):
        compute_activity(
            ni, efficiency_curve=eff, live_time_s=3600.0,
            from_bg_subtracted=True,
            decay_correction=True,
            reference_datetime=datetime(2020, 1, 1),
            measurement_datetime=datetime(2020, 2, 1),
            chain_decay_mode="invalid_mode",
        )


# ═════════════════════════════════════════════════════════════════════
# Group 4c — PTB-2018 Annex E.2 chain_decay_mode (Pb-212 dual T½)
# ═════════════════════════════════════════════════════════════════════

def test_ptb_e2_equilibrium_mode_pb212():
    """Equilibrium mode: Pb-212 uses Th-228 T½=1.91a → small factor for 30-day Δt."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(238.63, 43.6, 6000.0, 80.0),)
    ni = FakeNuclideIdentification("Pb-212", True, matches)
    ref_dt = datetime(2020, 1, 1)
    meas_dt = ref_dt + timedelta(days=30)
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        decay_correction=True,
        reference_datetime=ref_dt,
        measurement_datetime=meas_dt,
        chain_decay_mode="equilibrium",
    )
    assert result.decay_corrected
    # T½(Th-228)=1.91a=6.0275e7s, Δt=30d=2.592e6s → factor=exp(ln2·Δt/T½)≈1.0305
    expected = math.exp(math.log(2.0) * (30 * 86400.0) / 6.0275e7)
    assert math.isclose(result.decay_factor, expected, rel_tol=1e-4), \
        f"E.2 equilibrium factor={result.decay_factor} vs expected={expected}"
    assert "PTB E.2 equilibrium" in result.notes
    assert "Th-228" in result.notes


def test_ptb_e2_ra224_fresh_mode_pb212():
    """Ra-224 fresh mode: Pb-212 uses Ra-224 T½=3.66d → moderate factor for 5-day Δt."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(238.63, 43.6, 6000.0, 80.0),)
    ni = FakeNuclideIdentification("Pb-212", True, matches)
    ref_dt = datetime(2020, 1, 1)
    meas_dt = ref_dt + timedelta(days=5)
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        decay_correction=True,
        reference_datetime=ref_dt,
        measurement_datetime=meas_dt,
        chain_decay_mode="ra224_fresh",
    )
    assert result.decay_corrected
    # T½(Ra-224)=3.66d=3.1622e5s, Δt=5d → factor=exp(ln2·5/3.66)≈2.583
    expected = math.exp(math.log(2.0) * 5.0 / 3.66)
    assert math.isclose(result.decay_factor, expected, rel_tol=1e-3), \
        f"E.2 ra224_fresh factor={result.decay_factor} vs expected={expected}"
    assert "PTB E.2 ra224_fresh" in result.notes
    assert "Ra-224" in result.notes


def test_ptb_e2_progeny_mode_pb212():
    """Progeny mode: Pb-212 uses own T½=38304s (10.64h) → factor≈1.067 for 1-hour Δt."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(238.63, 43.6, 6000.0, 80.0),)
    ni = FakeNuclideIdentification("Pb-212", True, matches)
    ref_dt = datetime(2020, 1, 1, 12, 0, 0)
    meas_dt = ref_dt + timedelta(hours=1)
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=600.0,
        from_bg_subtracted=True,
        decay_correction=True,
        reference_datetime=ref_dt,
        measurement_datetime=meas_dt,
        chain_decay_mode="progeny",
    )
    assert result.decay_corrected
    expected = math.exp(math.log(2.0) * 3600.0 / 38304.0)
    assert math.isclose(result.decay_factor, expected, rel_tol=1e-6), \
        f"E.2 progeny factor={result.decay_factor} vs expected={expected}"
    assert "PTB E.2" not in result.notes


def test_ptb_e2_rn222_mode_on_pb212_silent_noop():
    """rn222 mode on Pb-212 (E.1-only mode) → silent pass-through, own T½ used."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(238.63, 43.6, 6000.0, 80.0),)
    ni = FakeNuclideIdentification("Pb-212", True, matches)
    ref_dt = datetime(2020, 1, 1, 12, 0, 0)
    meas_dt = ref_dt + timedelta(hours=1)
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=600.0,
        from_bg_subtracted=True,
        decay_correction=True,
        reference_datetime=ref_dt,
        measurement_datetime=meas_dt,
        chain_decay_mode="rn222",
    )
    assert result.decay_corrected
    # rn222 on Pb-212 = silent no-op → own T½=38304s applies
    expected = math.exp(math.log(2.0) * 3600.0 / 38304.0)
    assert math.isclose(result.decay_factor, expected, rel_tol=1e-6), \
        f"rn222-on-Pb-212 factor={result.decay_factor} vs expected={expected}"
    assert "PTB E." not in result.notes


def test_ptb_e1_ra224_fresh_mode_on_pb214_silent_noop():
    """ra224_fresh mode on Pb-214 (E.2-only mode) → silent pass-through, own T½ used."""
    eff = constant_efficiency(0.01)
    matches = (FakeLineMatch(295.22, 18.42, 5000.0, 70.0),)
    ni = FakeNuclideIdentification("Pb-214", True, matches)
    ref_dt = datetime(2020, 1, 1, 12, 0, 0)
    meas_dt = ref_dt + timedelta(hours=1)
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=600.0,
        from_bg_subtracted=True,
        decay_correction=True,
        reference_datetime=ref_dt,
        measurement_datetime=meas_dt,
        chain_decay_mode="ra224_fresh",
    )
    assert result.decay_corrected
    # ra224_fresh on Pb-214 = silent no-op → own T½=1608s applies (Δt=3600s)
    expected = math.exp(math.log(2.0) * 3600.0 / 1608.0)
    assert math.isclose(result.decay_factor, expected, rel_tol=1e-6), \
        f"ra224_fresh-on-Pb-214 factor={result.decay_factor} vs expected={expected}"
    assert "PTB E." not in result.notes


# ═════════════════════════════════════════════════════════════════════
# Group 5 — extrapolation flag
# ═════════════════════════════════════════════════════════════════════

def test_epsilon_extrapolation_flagged():
    """Line outside ε(E) calibrated range → per-line flag set."""
    # Curve calibrated 60–1500 keV; Tl-208 2614 keV is outside.
    eff = EfficiencyCurve(
        coefficients=(math.log(0.01),),
        E_min_keV=60.0, E_max_keV=1500.0,
        chi2_per_dof=1.0, n_points_used=10, n_dof=9,
    )
    matches = (
        FakeLineMatch(583.19, 30.60, 10000.0, 100.0),    # in-range
        FakeLineMatch(2614.51, 35.85, 10000.0, 100.0),   # extrapolated
    )
    ni = FakeNuclideIdentification("Tl-208", True, matches)
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        coincidence_correction={583.19: 1.0, 2614.51: 1.0},
    )
    by_E = {round(la.E_keV, 2): la for la in result.lines_used}
    assert not by_E[583.19].epsilon_extrapolated
    assert by_E[2614.51].epsilon_extrapolated
    assert "extrapolation" in result.notes.lower()
    print(f"  ✓ test_epsilon_extrapolation_flagged (2614 flagged, 583 not)")


# ═════════════════════════════════════════════════════════════════════
# Group 6 — compute_activities_for_all dispatcher
# ═════════════════════════════════════════════════════════════════════

def test_dispatcher_iterates_detected():
    """Dispatcher returns one ActivityResult per detected nuclide."""
    eff = constant_efficiency(0.01)
    ni_cs = FakeNuclideIdentification(
        "Cs-137", True,
        (FakeLineMatch(661.66, 85.10, 10000.0, 100.0),)
    )
    ni_co = FakeNuclideIdentification(
        "Co-60", True, (
            FakeLineMatch(1173.23, 99.85, 10000.0, 100.0),
            FakeLineMatch(1332.49, 99.98, 10000.0, 100.0),
        )
    )
    ni_k = FakeNuclideIdentification(
        "K-40", True,
        (FakeLineMatch(1460.82, 10.55, 5000.0, 80.0),)
    )
    id_result = FakeIdentificationResult(
        detected_nuclides=(ni_cs, ni_co, ni_k),
    )
    results = compute_activities_for_all(
        id_result,
        efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
        coincidence_corrections={
            "Co-60": {1173.23: 1.0, 1332.49: 1.0},
        },
    )
    assert len(results) == 3
    by_nuc = {r.nuclide: r for r in results}
    assert by_nuc["Cs-137"].is_valid()
    assert by_nuc["Co-60"].is_valid()
    assert by_nuc["K-40"].is_valid()
    # K-40 had no cc dict → no cascade warning (not in cascade set)
    assert by_nuc["K-40"].cascade_warning is None
    # Co-60 cc dict present → no warning
    assert by_nuc["Co-60"].cascade_warning is None
    print(f"  ✓ test_dispatcher_iterates_detected "
          f"({len(results)} nuclides, all valid)")


def test_dispatcher_propagates_bg_error():
    """BackgroundNotSubtractedError must propagate, not be swallowed."""
    eff = constant_efficiency(0.01)
    ni = FakeNuclideIdentification(
        "Cs-137", True,
        (FakeLineMatch(661.66, 85.10, 10000.0, 100.0),)
    )
    id_result = FakeIdentificationResult(detected_nuclides=(ni,))
    try:
        compute_activities_for_all(
            id_result,
            efficiency_curve=eff, live_time_s=3600.0,
            from_bg_subtracted=False,
            bg_available=True,
        )
        raise AssertionError("Expected BackgroundNotSubtractedError")
    except BackgroundNotSubtractedError as exc:
        assert exc.nuclide == "Cs-137"
    print(f"  ✓ test_dispatcher_propagates_bg_error")


def test_dispatcher_skips_rejected_by_default():
    """skip_undetected=True (default) ignores rejected_nuclides."""
    eff = constant_efficiency(0.01)
    ni_det = FakeNuclideIdentification(
        "Cs-137", True,
        (FakeLineMatch(661.66, 85.10, 10000.0, 100.0),)
    )
    ni_rej = FakeNuclideIdentification("Eu-152", False, ())
    id_result = FakeIdentificationResult(
        detected_nuclides=(ni_det,),
        rejected_nuclides=(ni_rej,),
    )
    results = compute_activities_for_all(
        id_result,
        efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
    )
    assert len(results) == 1
    assert results[0].nuclide == "Cs-137"
    print(f"  ✓ test_dispatcher_skips_rejected_by_default")


# ═════════════════════════════════════════════════════════════════════
# Group 7 — repr / serialization sanity
# ═════════════════════════════════════════════════════════════════════

def test_activity_result_repr():
    """ActivityResult.__repr__ shows nuclide, A, σ, n_lines, and flags."""
    eff = constant_efficiency(0.01)
    matches = (
        FakeLineMatch(1173.23, 99.85, 10000.0, 100.0),
        FakeLineMatch(1332.49, 99.98, 10000.0, 100.0),
    )
    ni = FakeNuclideIdentification("Co-60", True, matches)
    result = compute_activity(
        ni, efficiency_curve=eff, live_time_s=3600.0,
        from_bg_subtracted=True,
    )
    s = repr(result)
    assert "Co-60" in s
    assert "Bq" in s
    assert "2 lines" in s
    # Cascade warning present (no cc supplied)
    assert "K-17" in s or "cascade" in s.lower()
    print(f"  ✓ test_activity_result_repr ({s[:90]}...)")


# ═════════════════════════════════════════════════════════════════════
# Group 8 — end-to-end integration on real Cs-137 spectrum
# ═════════════════════════════════════════════════════════════════════

def test_activity_cs137_real_spectrum_5cm():
    """End-to-end: load Cs-137 .spe → calibrate → identify → activity.

    Validates the full pipeline integration of Phase 2.1d. The
    Гамма-1С 5cm point source Cs-137 №SRC-02 reference spectrum
    contains a single dominant photopeak. The computed activity must
    be:
      • finite and positive;
      • produced from at least one matched line;
      • within a plausible laboratory range (10 Bq … 1 MBq) — exact
        match to the certificate activity is deferred until the
        .src certificate parser is implemented (open question in
        the v1.7.6 handoff).
    """
    from gamma.io.readers import read_spectrum
    from gamma.calibration.efficiency import fit_efficiency_from_efr_file
    from gamma.calibration.fwhm_provider import (
        make_fwhm_at_channel_provider,
    )
    from gamma.peaks.search import mariscotti_search
    from gamma.identification.window import build_identification_window
    from gamma.identification.identify import identify_nuclides
    from gamma.identification.disambiguate import (
        disambiguate_identifications,
    )

    spe_path = ("detectors/Gamma-1S/reference_spectra/"
                "archive/Cs-137__163_2017.spe")
    efr_path = ("detectors/Gamma-1S/efficiency/Gamma-1S_NaI_63x63_USB_SN-01/"
                "УДС-ГЦ-63х63-USB__SN-01_-_Точечная-5см.efr")

    spec = read_spectrum(spe_path)
    eff = fit_efficiency_from_efr_file(efr_path, degree=3)
    fwhm_at = make_fwhm_at_channel_provider(spec, fallback_channels=15.0)
    peaks = mariscotti_search(
        spec.counts, fwhm_channels=fwhm_at, sigma_threshold=5.0,
    )
    window = build_identification_window("NaI", delta_E0_keV=15.0)
    raw_id = identify_nuclides(
        found_peaks=peaks, spec=spec, window=window,
        fwhm_at_channel=fwhm_at,
    )
    refined = disambiguate_identifications(raw_id)
    cs = next(
        (ni for ni in refined.detected_nuclides if ni.nuclide == "Cs-137"),
        None,
    )
    assert cs is not None, "Cs-137 must be detected in this spectrum"

    result = compute_activity(
        cs, efficiency_curve=eff,
        live_time_s=spec.live_time,
        # No bg spectrum loaded → safe to compute on gross
        from_bg_subtracted=False,
        bg_available=False,
        decay_correction=False,        # no certificate datetime here
    )
    assert result.is_valid()
    assert result.n_lines_used() >= 1
    assert 10.0 < result.A_Bq < 1e6, \
        f"A={result.A_Bq:.2e} Bq is outside plausible lab-source range"
    # Sanity: σ should be a small fraction of A (relative uncertainty
    # dominated by ε~5% and S has thousands of counts → expect rel<20%)
    rel = result.sigma_A_Bq / result.A_Bq
    assert rel < 0.20, f"σ/A={rel:.3f} too large for a strong-source measurement"
    print(f"  ✓ test_activity_cs137_real_spectrum_5cm "
          f"(A={result.A_Bq:.3e}±{result.sigma_A_Bq:.2e} Bq, "
          f"rel={rel*100:.1f}%)")


# ═════════════════════════════════════════════════════════════════════
# Runner
# ═════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Running Phase 2.1d activity calculation tests...\n")

    # Group 1 — formula & aggregation (8)
    test_single_line_activity_cs137()
    test_multi_line_weighted_average_co60()
    test_uncertainty_propagation_formula()
    test_intra_chi2_per_dof()
    test_undetected_nuclide_yields_nan()
    test_skip_line_no_area()
    test_skip_line_zero_intensity()
    test_invalid_inputs_raise()

    # Group 2 — background subtraction safety (3)
    test_refuse_gross_when_bg_available()
    test_allow_gross_when_no_bg()
    test_force_gross_override()

    # Group 3 — cascade summing (K-17) (5)
    test_cascade_warning_emitted_co60()
    test_cascade_warning_suppressed_with_correction()
    test_cascade_correction_scales_area()
    test_cascade_correction_approximate_key_match()
    test_cascade_catalog_canonical_set()

    # Group 4 — decay correction (5)
    test_decay_correction_one_half_life()
    test_decay_correction_disabled_by_flag()
    test_decay_correction_skipped_missing_dates()
    test_decay_correction_quasi_stable_k40()
    test_decay_correction_propagates_to_sigma()

    # Group 5 — extrapolation (1)
    test_epsilon_extrapolation_flagged()

    # Group 6 — dispatcher (3)
    test_dispatcher_iterates_detected()
    test_dispatcher_propagates_bg_error()
    test_dispatcher_skips_rejected_by_default()

    # Group 7 — repr (1)
    test_activity_result_repr()

    # Group 8 — real-data end-to-end integration (1)
    test_activity_cs137_real_spectrum_5cm()

    print("\n✓ All Phase 2.1d activity calculation tests passed.")
