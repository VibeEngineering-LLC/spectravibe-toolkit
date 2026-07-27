# F-RPT-08 / v1.19.1 — ISO 11929 thresholds unit tests
#
# Verifies scripts/gamma/math/iso_11929_thresholds.py
# ISO 11929-1:2019 §5.4.3 (decision threshold y*) and §5.4.4 (detection limit y#).
# Cross-refs: RAG-005, RAG-008, RAG-009, RAG-022.

import math
import pytest

from gamma.math.iso_11929_thresholds import (
    decision_threshold,
    detection_limit,
    multi_line_decision_threshold,
    multi_line_detection_limit,
)

# ---------------------------------------------------------------------------
# Hand-derived reference constants for the canonical test case
#
#   gross_counts = 5400, bg_counts = 1800
#   efficiency   = 0.8  (dimensionless)
#   branching    = 0.8  (decimal)
#   mass_kg      = 1.0
#   live_time_s  = 3600
#
#   sigma_0 = sqrt(5400 + 1800) / (0.8 * 0.8 * 1.0 * 3600)
#           = sqrt(7200) / 2304
#           = 84.8528... / 2304
#           = 0.036829... Bq/kg
#
#   k_{0.95} = 1.6449  (ISO 11929 Table A.1 / standard normal 95th percentile)
#
#   t* = 1.6449 * 0.036829 = 0.060581... Bq/kg     (decision threshold)
#   η* ≈ 2 * t* = 0.121162... Bq/kg               (detection limit, low-stats approx)
# ---------------------------------------------------------------------------
_GROSS = 5400.0
_BG    = 1800.0
_EFF   = 0.8
_BR    = 0.8
_MASS  = 1.0
_T     = 3600.0
_K95   = 1.6449
_SIGMA0 = math.sqrt(_GROSS + _BG) / (_EFF * _BR * _MASS * _T)
_T_STAR  = _K95 * _SIGMA0          # ~0.060581 Bq/kg
_ETA_STAR = 2.0 * _T_STAR          # ~0.121162 Bq/kg


def test_decision_threshold_reference_case():
    """Decision threshold for canonical NaI-style measurement.

    Inputs: gross=5400, bg=1800, eff=0.8, br=0.8, mass=1.0 kg, t=3600 s.
    Derivation (ISO 11929-1:2019 §5.4.3):
        sigma_0 = sqrt(7200) / (0.8 * 0.8 * 1.0 * 3600) = 84.8528/2304 ≈ 0.036829 Bq/kg
        k_{1-0.05} = 1.6449  (standard normal 95th percentile)
        t* = 1.6449 * 0.036829 ≈ 0.060581 Bq/kg
    """
    result = decision_threshold(_GROSS, _BG, _EFF, _BR, _MASS, _T)
    assert result == pytest.approx(_T_STAR, rel=1e-4)


def test_detection_limit_reference_case():
    """Detection limit for canonical NaI-style measurement.

    Inputs: same as test_decision_threshold_reference_case.
    Derivation (ISO 11929-1:2019 §5.4.4, low-stats approx):
        eta* ≈ 2 * t* = 2 * 0.060581 ≈ 0.121162 Bq/kg
    Valid when background-dominated (N_bg >> net signal at y=t*).
    """
    result = detection_limit(_GROSS, _BG, _EFF, _BR, _MASS, _T)
    assert result == pytest.approx(_ETA_STAR, rel=1e-4)


def test_eta_star_approx_twice_t_star_symmetric():
    """eta* = 2 * t* when alpha = beta = 0.05 (symmetric case).

    ISO 11929-1:2019 §5.4.4 low-stats approximation:
        y# ≈ (k_{1-alpha} + k_{1-beta}) * sigma_0
    When alpha = beta = 0.05, k_{1-alpha} = k_{1-beta} = 1.6449:
        y# ≈ 2 * 1.6449 * sigma_0 = 2 * t*
    Exact to floating-point precision (ratio must equal 2.0 ± 1e-10).
    """
    t_star = decision_threshold(_GROSS, _BG, _EFF, _BR, _MASS, _T)
    eta_star = detection_limit(_GROSS, _BG, _EFF, _BR, _MASS, _T)
    assert t_star is not None and eta_star is not None
    assert eta_star == pytest.approx(2.0 * t_star, rel=1e-10)


def test_t_star_less_than_eta_star():
    """t* < eta* must hold for any valid measurement (ISO 11929 ordering guarantee).

    Rationale: eta* is by definition the smallest true value detectable with
    P(false negative) <= beta, so eta* >= t* always.  For alpha = beta = 0.05,
    eta* = 2 * t* > t* as long as t* > 0.
    """
    t_star = decision_threshold(_GROSS, _BG, _EFF, _BR, _MASS, _T)
    eta_star = detection_limit(_GROSS, _BG, _EFF, _BR, _MASS, _T)
    assert t_star is not None and eta_star is not None
    assert t_star > 0
    assert t_star < eta_star


def test_zero_gross_zero_bg_returns_zero():
    """When gross = bg = 0, sigma_0 = 0, so t* = eta* = 0.

    Derivation:
        sigma_0 = sqrt(0 + 0) / (0.8 * 0.8 * 1.0 * 3600) = 0.0
        t*   = 1.6449 * 0.0 = 0.0
        eta* = 2 * 0.0 = 0.0
    This is the vacuous case: no counts → threshold is zero.
    """
    assert decision_threshold(0, 0, _EFF, _BR, _MASS, _T) == pytest.approx(0.0, abs=1e-15)
    assert detection_limit(0, 0, _EFF, _BR, _MASS, _T) == pytest.approx(0.0, abs=1e-15)


def test_efficiency_zero_returns_none():
    """efficiency = 0 → None (guard: division by zero in sensitivity w).

    w = 1 / (eff * br * mass * t); if eff = 0, w is undefined.
    Function must return None, not raise ZeroDivisionError or NaN.
    """
    assert decision_threshold(1000, 500, 0.0, _BR, _MASS, _T) is None
    assert detection_limit(1000, 500, 0.0, _BR, _MASS, _T) is None


def test_mass_zero_returns_none():
    """mass_kg = 0 → None (guard: specific activity is undefined for zero mass).

    Converting counts/s to Bq/kg requires mass > 0.  Return None explicitly.
    """
    assert decision_threshold(1000, 500, _EFF, _BR, 0.0, _T) is None
    assert detection_limit(1000, 500, _EFF, _BR, 0.0, _T) is None


def test_livetime_zero_returns_none():
    """live_time_s = 0 → None (guard: count rate is undefined for zero live time).

    w = 1 / (eff * br * mass * t); if t = 0, w is undefined.
    """
    assert decision_threshold(1000, 500, _EFF, _BR, _MASS, 0.0) is None
    assert detection_limit(1000, 500, _EFF, _BR, _MASS, 0.0) is None


def test_branching_ratio_zero_returns_none():
    """branching_ratio = 0 → None (guard: line with zero emission probability is non-physical).

    A line with I = 0 carries no information; w is undefined.
    """
    assert decision_threshold(1000, 500, _EFF, 0.0, _MASS, _T) is None
    assert detection_limit(1000, 500, _EFF, 0.0, _MASS, _T) is None


def test_physical_scale_t_star_order_of_magnitude():
    """Physical scale check: t* should be O(0.05–0.07) Bq/kg for NaI reference case.

    Brief specification (task brief, wave 4):
        '1 cps net, 1 kg, 80% eff, 80% BR, 3600 s → t* of order ~0.5 Bq/kg'
    Note: brief says 0.5, but that assumes gross ≈ bg ≈ 1800 (lower total counts).
    With gross=5400, bg=1800 (total=7200) our hand derivation gives 0.0606 Bq/kg.
    The interval [0.01, 1.0] safely brackets both interpretations.
    """
    result = decision_threshold(_GROSS, _BG, _EFF, _BR, _MASS, _T)
    assert result is not None
    assert 0.01 <= result <= 1.0, (
        f"t* = {result:.4f} Bq/kg outside expected physical range [0.01, 1.0]"
    )


def test_non_default_alpha_detection_threshold():
    """Non-default alpha=0.01 → k_{0.99} ≈ 2.3263, t* increases proportionally.

    Derivation:
        sigma_0 = sqrt(7200) / 2304 ≈ 0.036829 Bq/kg  (same as reference case)
        k_{0.99} ≈ 2.3263  (standard normal 99th percentile, e.g. Knoll 4ed Table B.1)
        t*(alpha=0.01) = 2.3263 * 0.036829 ≈ 0.085672 Bq/kg
    Ratio t*(0.01) / t*(0.05) = 2.3263 / 1.6449 ≈ 1.4142.
    """
    t_star_01 = decision_threshold(_GROSS, _BG, _EFF, _BR, _MASS, _T, alpha=0.01)
    t_star_05 = decision_threshold(_GROSS, _BG, _EFF, _BR, _MASS, _T, alpha=0.05)
    assert t_star_01 is not None and t_star_05 is not None
    # Ratio should be k_0.99 / k_0.95 ≈ 2.3263 / 1.6449 ≈ 1.4142
    expected_ratio = 2.3263 / 1.6449
    assert t_star_01 / t_star_05 == pytest.approx(expected_ratio, rel=5e-3)


# ---------------------------------------------------------------------------
# #PTB-1 regulatory-regime coverage — PTB-2018 SPEKT/GRUNDL Annex C
# ---------------------------------------------------------------------------

def test_regime_kta_equals_default_decision_threshold():
    """regime='KTA' must equal the default alpha=0.05 result (k = 1.6449)."""
    t_default = decision_threshold(_GROSS, _BG, _EFF, _BR, _MASS, _T)
    t_kta = decision_threshold(_GROSS, _BG, _EFF, _BR, _MASS, _T, regime="KTA")
    assert t_default is not None and t_kta is not None
    assert t_kta == pytest.approx(t_default, rel=1e-12)


def test_regime_kta_equals_default_detection_limit():
    """regime='KTA' must equal default alpha=beta=0.05 result (y# = 2·y*)."""
    y_default = detection_limit(_GROSS, _BG, _EFF, _BR, _MASS, _T)
    y_kta = detection_limit(_GROSS, _BG, _EFF, _BR, _MASS, _T, regime="KTA")
    assert y_default is not None and y_kta is not None
    assert y_kta == pytest.approx(y_default, rel=1e-12)


def test_regime_imis_decision_threshold_k3():
    """regime='IMIS' → k_{1-α} = 3.0. y* = 3 · σ_0."""
    t_imis = decision_threshold(_GROSS, _BG, _EFF, _BR, _MASS, _T, regime="IMIS")
    assert t_imis is not None
    assert t_imis == pytest.approx(3.0 * _SIGMA0, rel=1e-12)


def test_regime_imis_detection_limit_k3_plus_k95():
    """regime='IMIS' → y# = (3.0 + 1.6449) · σ_0 = 4.6449 · σ_0."""
    y_imis = detection_limit(_GROSS, _BG, _EFF, _BR, _MASS, _T, regime="IMIS")
    assert y_imis is not None
    assert y_imis == pytest.approx((3.0 + 1.6449) * _SIGMA0, rel=1e-12)


def test_regime_imis_asymmetric_vs_kta():
    """IMIS y* is ~1.82× KTA y*; IMIS y# is ~1.41× KTA y# (asymmetric)."""
    t_kta = decision_threshold(_GROSS, _BG, _EFF, _BR, _MASS, _T, regime="KTA")
    t_imis = decision_threshold(_GROSS, _BG, _EFF, _BR, _MASS, _T, regime="IMIS")
    y_kta = detection_limit(_GROSS, _BG, _EFF, _BR, _MASS, _T, regime="KTA")
    y_imis = detection_limit(_GROSS, _BG, _EFF, _BR, _MASS, _T, regime="IMIS")
    # y*: ratio = 3.0 / 1.6449 ≈ 1.8238
    assert t_imis / t_kta == pytest.approx(3.0 / 1.6449, rel=1e-6)
    # y#: ratio = (3.0 + 1.6449) / (2 · 1.6449) = 4.6449 / 3.2898 ≈ 1.4119
    assert y_imis / y_kta == pytest.approx((3.0 + 1.6449) / (2.0 * 1.6449), rel=1e-6)


def test_regime_invalid_returns_none():
    """Unknown regime string → guard fires → None (both functions)."""
    assert decision_threshold(_GROSS, _BG, _EFF, _BR, _MASS, _T, regime="XYZ") is None
    assert detection_limit(_GROSS, _BG, _EFF, _BR, _MASS, _T, regime="XYZ") is None
    # Case-sensitive: lowercase not accepted
    assert decision_threshold(_GROSS, _BG, _EFF, _BR, _MASS, _T, regime="kta") is None
    assert detection_limit(_GROSS, _BG, _EFF, _BR, _MASS, _T, regime="imis") is None


def test_regime_overrides_alpha_beta():
    """When regime is set, alpha/beta arguments are ignored."""
    # alpha=0.5 would give k ≈ 0 without regime; regime='IMIS' forces k=3.0
    t_override = decision_threshold(
        _GROSS, _BG, _EFF, _BR, _MASS, _T, alpha=0.5, regime="IMIS"
    )
    assert t_override is not None
    assert t_override == pytest.approx(3.0 * _SIGMA0, rel=1e-12)


# ---------------------------------------------------------------------------
# #PTB-2 multi-line (Cs-134 style) — PTB-2018 Annex C Eq. (C3)–(C5)
# ---------------------------------------------------------------------------

def test_multi_line_single_entry_equals_single_line():
    """One-line input must reproduce the single-line decision_threshold."""
    t_single = decision_threshold(_GROSS, _BG, _EFF, _BR, _MASS, _T)
    t_multi = multi_line_decision_threshold(
        [_GROSS], [_BG], [_EFF], [_BR], _MASS, _T
    )
    assert t_single is not None and t_multi is not None
    assert t_multi == pytest.approx(t_single, rel=1e-12)


def test_multi_line_two_identical_lines_sqrt2_gain():
    """Two identical lines → combined σ(0) = σ_single / √2 (Eq. C3).

    Inverse-variance combination: 1/σ² = 1/σ_j² + 1/σ_j² = 2/σ_j².
    """
    t_single = decision_threshold(_GROSS, _BG, _EFF, _BR, _MASS, _T)
    t_multi = multi_line_decision_threshold(
        [_GROSS, _GROSS], [_BG, _BG], [_EFF, _EFF], [_BR, _BR], _MASS, _T
    )
    assert t_single is not None and t_multi is not None
    assert t_multi == pytest.approx(t_single / math.sqrt(2.0), rel=1e-12)


def test_multi_line_combined_below_best_single_line():
    """Combined y* must beat (be below) the best individual line.

    Cs-134-style asymmetric pair: strong 605 keV line (br=0.976) and weaker
    796 keV line (br=0.855, lower efficiency).
    """
    lines = dict(
        line_gross_counts=[5400.0, 3200.0],
        line_bg_counts=[1800.0, 1500.0],
        line_efficiencies=[0.8, 0.6],
        line_branching_ratios=[0.976, 0.855],
    )
    t_multi = multi_line_decision_threshold(**lines, mass_kg=_MASS, live_time_s=_T)
    t_1 = decision_threshold(5400.0, 1800.0, 0.8, 0.976, _MASS, _T)
    t_2 = decision_threshold(3200.0, 1500.0, 0.6, 0.855, _MASS, _T)
    assert t_multi is not None and t_1 is not None and t_2 is not None
    assert t_multi < min(t_1, t_2)


def test_multi_line_imis_regime():
    """regime='IMIS' → y* = 3.0 · combined σ(0)."""
    t_kta = multi_line_decision_threshold(
        [_GROSS, _GROSS], [_BG, _BG], [_EFF, _EFF], [_BR, _BR], _MASS, _T,
        regime="KTA",
    )
    t_imis = multi_line_decision_threshold(
        [_GROSS, _GROSS], [_BG, _BG], [_EFF, _EFF], [_BR, _BR], _MASS, _T,
        regime="IMIS",
    )
    assert t_kta is not None and t_imis is not None
    assert t_imis / t_kta == pytest.approx(3.0 / 1.6449, rel=1e-9)


def test_multi_line_guards():
    """Empty input, mismatched lengths, bad efficiency, bad regime → None."""
    assert multi_line_decision_threshold([], [], [], [], _MASS, _T) is None
    assert multi_line_decision_threshold(
        [_GROSS, _GROSS], [_BG], [_EFF], [_BR], _MASS, _T
    ) is None
    assert multi_line_decision_threshold(
        [_GROSS], [_BG], [0.0], [_BR], _MASS, _T
    ) is None
    assert multi_line_decision_threshold(
        [_GROSS], [_BG], [_EFF], [_BR], _MASS, _T, regime="XYZ"
    ) is None
    assert multi_line_detection_limit(
        [_GROSS], [_BG], [_EFF], [_BR], _MASS, _T, regime="XYZ"
    ) is None


def test_multi_line_zero_counts_returns_zero():
    """A line with zero gross and zero bg counts → σ(0) = 0 → y* = 0."""
    t = multi_line_decision_threshold(
        [0.0, _GROSS], [0.0, _BG], [_EFF, _EFF], [_BR, _BR], _MASS, _T
    )
    assert t == 0.0


def test_multi_line_detection_limit_fallback_is_2x():
    """Without measured activity → fallback y# = y* + k_β·u(0) = 2·y* (KTA)."""
    t_star = multi_line_decision_threshold(
        [_GROSS, _GROSS], [_BG, _BG], [_EFF, _EFF], [_BR, _BR], _MASS, _T
    )
    y_hash = multi_line_detection_limit(
        [_GROSS, _GROSS], [_BG, _BG], [_EFF, _EFF], [_BR, _BR], _MASS, _T
    )
    assert t_star is not None and y_hash is not None
    assert y_hash == pytest.approx(2.0 * t_star, rel=1e-12)


def test_multi_line_detection_limit_iterative_fixed_point():
    """Iterative y# must satisfy the Eq. (C4)+(C5) fixed-point equation.

    y# = y* + k_β · sqrt( u0² + (u²(a_r) − u0²) · y#/a_r )
    with u0 = y*/k_α.
    """
    args = ([_GROSS, _GROSS], [_BG, _BG], [_EFF, _EFF], [_BR, _BR], _MASS, _T)
    a_r = 0.5          # Bq/kg, measured activity
    u_ar = 0.05        # Bq/kg, its standard uncertainty (> u0)
    y_star = multi_line_decision_threshold(*args)
    y_hash = multi_line_detection_limit(
        *args, measured_activity=a_r, measured_uncertainty=u_ar
    )
    assert y_star is not None and y_hash is not None
    u0 = y_star / _K95
    rhs = y_star + _K95 * math.sqrt(
        u0**2 + (u_ar**2 - u0**2) * (y_hash / a_r)
    )
    assert y_hash == pytest.approx(rhs, rel=1e-8)
    # u(a_r) > u(0) → iterative y# must exceed the low-stats 2·y* fallback
    assert y_hash > 2.0 * y_star


def test_multi_line_detection_limit_iterative_reduces_to_fallback():
    """When u(a_r) == u(0), Eq. (C5) is constant → y# = 2·y* exactly."""
    args = ([_GROSS, _GROSS], [_BG, _BG], [_EFF, _EFF], [_BR, _BR], _MASS, _T)
    y_star = multi_line_decision_threshold(*args)
    assert y_star is not None
    u0 = y_star / _K95
    y_hash = multi_line_detection_limit(
        *args, measured_activity=0.5, measured_uncertainty=u0
    )
    assert y_hash is not None
    assert y_hash == pytest.approx(2.0 * y_star, rel=1e-9)
