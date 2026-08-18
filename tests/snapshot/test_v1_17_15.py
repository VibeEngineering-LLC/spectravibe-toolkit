# -*- coding: utf-8 -*-
"""
v1.17.15 delivery tests — Calibration robustness.

Covers F-284..F-287 (T-029, T-059, T-060, T-031).
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest

SCRIPTS = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts")
)
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)


# ─────────────────────────────────────────────────────────────────────
# F-284 / T-029 — 2-zone efficiency
# ─────────────────────────────────────────────────────────────────────

def test_F284_two_zone_basic_fit():
    from gamma.calibration.two_zone_efficiency import (
        fit_two_zone_efficiency, evaluate_two_zone_efficiency,
    )
    # Smoothly varying ε(E) — synthetic, no scatter
    E = np.array([60.0, 122.0, 186.0, 250.0, 661.0, 1173.0, 1332.0, 1461.0, 2614.0])
    eps = np.array([0.025, 0.022, 0.018, 0.015, 0.008, 0.0055, 0.0050, 0.0045, 0.0025])
    fit = fit_two_zone_efficiency(E, eps)
    assert fit.n_low_anchors >= 3
    assert fit.n_high_anchors >= 4
    # Eps at 661 should be ~0.008
    eps_661 = evaluate_two_zone_efficiency(fit, 661.0)
    assert 0.005 < eps_661 < 0.012


def test_F284_two_zone_c1_continuity_at_split():
    """Значение ε на split должно совпадать в обеих зонах."""
    from gamma.calibration.two_zone_efficiency import (
        fit_two_zone_efficiency, evaluate_two_zone_efficiency,
    )
    E = np.array([60.0, 122.0, 186.0, 250.0, 661.0, 1173.0, 1332.0, 1461.0, 2614.0])
    eps = np.array([0.025, 0.022, 0.018, 0.015, 0.008, 0.0055, 0.0050, 0.0045, 0.0025])
    fit = fit_two_zone_efficiency(E, eps, split_E_keV=250.0)
    # Approach from below and above
    eps_below = evaluate_two_zone_efficiency(fit, 249.99)
    eps_above = evaluate_two_zone_efficiency(fit, 250.01)
    rel_diff = abs(eps_above - eps_below) / max(eps_below, 1e-9)
    assert rel_diff < 0.05, f"C¹ discontinuity {rel_diff:.2%}"


def test_F284_too_few_points_raises():
    from gamma.calibration.two_zone_efficiency import fit_two_zone_efficiency
    with pytest.raises(ValueError):
        fit_two_zone_efficiency([100.0, 661.0], [0.02, 0.008])


def test_F284_negative_eps_raises():
    from gamma.calibration.two_zone_efficiency import fit_two_zone_efficiency
    with pytest.raises(ValueError):
        fit_two_zone_efficiency([100.0, 300.0, 661.0], [0.02, -0.01, 0.008])


# ─────────────────────────────────────────────────────────────────────
# F-285 / T-059 — GOST INL / NL metrics
# ─────────────────────────────────────────────────────────────────────

def test_F285_perfect_calibration_zero_inl():
    from gamma.calibration.gost_metrics import compute_gost_linearity_metrics
    E = [122.0, 661.0, 1173.0, 1332.0, 1461.0, 2614.0]
    metrics = compute_gost_linearity_metrics(E, E)
    assert metrics.INL_pct_of_full_scale == 0.0
    assert metrics.accepts_inl


def test_F285_inl_calculation():
    """Residual 2 кэВ на full_scale 2614 = ~0.077 %."""
    from gamma.calibration.gost_metrics import compute_gost_linearity_metrics
    anchors = [661.0, 1461.0, 2614.0]
    fitted = [661.0, 1461.0, 2616.0]   # +2 кэВ residual at 2614
    metrics = compute_gost_linearity_metrics(anchors, fitted)
    expected_inl = 100.0 * 2.0 / 2614.0
    assert metrics.INL_pct_of_full_scale == pytest.approx(expected_inl, abs=0.001)


def test_F285_inl_fails_at_2_pct():
    from gamma.calibration.gost_metrics import compute_gost_linearity_metrics
    anchors = [661.0, 2614.0]
    fitted = [661.0, 2700.0]    # ≈ 3.3 % INL
    metrics = compute_gost_linearity_metrics(anchors, fitted)
    assert metrics.INL_pct_of_full_scale > 2.0
    assert not metrics.accepts_inl


# ─────────────────────────────────────────────────────────────────────
# F-286 / T-060 — Calibration verification loop
# ─────────────────────────────────────────────────────────────────────

def test_F286_perfect_match_accepts():
    from gamma.calibration.verification_loop import (
        make_standard_comparison, verify_calibration_against_standards,
    )
    stds = [
        make_standard_comparison(
            nuclide="Cs-137", line_keV=661.66,
            A_passport_Bq=1000.0, A_passport_unc_Bq=10.0,
            A_computed_Bq=1005.0, A_computed_unc_Bq=15.0,
        ),
        make_standard_comparison(
            nuclide="Co-60", line_keV=1173.23,
            A_passport_Bq=500.0, A_passport_unc_Bq=5.0,
            A_computed_Bq=498.0, A_computed_unc_Bq=8.0,
        ),
    ]
    out = verify_calibration_against_standards(stds)
    assert out.acceptable
    assert out.n_within_2_sigma == 2


def test_F286_failing_standard_marked():
    from gamma.calibration.verification_loop import (
        make_standard_comparison, verify_calibration_against_standards,
    )
    # 5-sigma deviation
    s = make_standard_comparison(
        nuclide="Cs-137", line_keV=661.66,
        A_passport_Bq=1000.0, A_passport_unc_Bq=5.0,
        A_computed_Bq=1100.0, A_computed_unc_Bq=10.0,
    )
    assert not s.within_2_sigma
    out = verify_calibration_against_standards([s])
    assert not out.acceptable


def test_F286_threshold_20_percent():
    """Если 1 из 5 эталонов выпал → fraction=0.2 → ровно на пороге → accept."""
    from gamma.calibration.verification_loop import (
        make_standard_comparison, verify_calibration_against_standards,
    )
    good = make_standard_comparison(
        nuclide="X", line_keV=100.0,
        A_passport_Bq=100.0, A_passport_unc_Bq=1.0,
        A_computed_Bq=100.0, A_computed_unc_Bq=1.0,
    )
    bad = make_standard_comparison(
        nuclide="Y", line_keV=200.0,
        A_passport_Bq=100.0, A_passport_unc_Bq=1.0,
        A_computed_Bq=200.0, A_computed_unc_Bq=1.0,
    )
    stds = [good, good, good, good, bad]
    out = verify_calibration_against_standards(stds)
    assert out.fraction_failing == 0.2
    assert out.acceptable   # ровно на пороге → accept


# ─────────────────────────────────────────────────────────────────────
# F-287 / T-031 — bg subtraction strategy
# ─────────────────────────────────────────────────────────────────────

def test_F287_nai_strategy_is_channel_by_channel():
    from gamma.calibration.bg_strategy import recommend_bg_strategy
    assert recommend_bg_strategy("NaI").strategy == "channel_by_channel"
    assert recommend_bg_strategy("CsI").strategy == "channel_by_channel"


def test_F287_hpge_strategy_is_per_peak():
    from gamma.calibration.bg_strategy import recommend_bg_strategy
    assert recommend_bg_strategy("HPGe").strategy == "per_peak"
    assert recommend_bg_strategy("LaBr3").strategy == "per_peak"


def test_F287_subtract_bg_channel_by_channel():
    from gamma.calibration.bg_strategy import subtract_bg
    counts = np.array([100.0, 200.0, 300.0])
    bg = np.array([20.0, 30.0, 40.0])
    net = subtract_bg(counts, bg, strategy="channel_by_channel")
    assert np.allclose(net, [80.0, 170.0, 260.0])


def test_F287_subtract_bg_with_time_ratio():
    from gamma.calibration.bg_strategy import subtract_bg
    counts = np.array([100.0, 200.0])
    bg = np.array([50.0, 100.0])
    # sample 2× longer than bg → bg должен умножиться на 2
    net = subtract_bg(counts, bg, bg_live_time_ratio=2.0)
    assert np.allclose(net, [0.0, 0.0])
