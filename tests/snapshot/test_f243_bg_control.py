"""F-243 / v1.18.29 — BG control 10% gate (sum_Y, t_live, rate-ratio).

Validates the three independent gates of
``gamma.io.bg_control.validate_background``:

  1. ``|F - F_ref| / F_ref < rate_tolerance``  (default 10 %)
  2. ``sum_Y >= min_sum_counts``                (default 1000)
  3. ``t_live >= min_live_time_s``              (default 600 s)

When ``F_ref <= 0`` the rate-ratio check is skipped; the other two gates
must still operate.
"""
from __future__ import annotations

import math
import sys

import pytest

sys.path.insert(0, "scripts")

from gamma.io.bg_control import (  # noqa: E402
    validate_background,
    validate_background_roi,
    BgControlResult,
    check_bg_quality,
    BgQualityReport,
    bg_z_test,
    BgZTestResult,
    Z_TIER_STABLE_MAX,
    Z_TIER_BORDERLINE_MAX,
)
from types import SimpleNamespace  # noqa: E402


def _mk_spec(counts, live_time):
    """Lightweight Spectrum stand-in (duck-typed: counts + live_time)."""
    return SimpleNamespace(counts=list(counts), live_time=float(live_time))


def test_all_gates_pass_happy_path():
    """All three gates pass cleanly (F~F_ref, sum>=1000, t_live>=600)."""
    r = validate_background(F=10.0, F_ref=9.5, sum_Y=5000.0, t_live=1200.0)
    assert isinstance(r, BgControlResult)
    assert r.ok is True
    assert r.failures == ()
    # rate_ratio == 10/9.5
    assert math.isclose(r.rate_ratio, 10.0 / 9.5, rel_tol=1e-9)
    assert r.sum_counts == 5000.0
    assert r.live_time_s == 1200.0


def test_rate_mismatch_fails():
    """|F-F0|/F0 = 0.15 > 0.1 → rate-mismatch failure."""
    r = validate_background(F=11.5, F_ref=10.0, sum_Y=5000.0, t_live=1200.0)
    assert r.ok is False
    assert len(r.failures) == 1
    assert "rate" in r.failures[0].lower()
    # Other gates passed
    assert r.sum_counts == 5000.0
    assert r.live_time_s == 1200.0


def test_low_counts_fails():
    """sum_Y = 500 < 1000 → counts failure."""
    r = validate_background(F=10.0, F_ref=10.0, sum_Y=500.0, t_live=1200.0)
    assert r.ok is False
    assert len(r.failures) == 1
    assert "sum_counts" in r.failures[0]
    assert "500" in r.failures[0]


def test_short_live_time_fails():
    """t_live = 300 < 600 → live-time failure."""
    r = validate_background(F=10.0, F_ref=10.0, sum_Y=5000.0, t_live=300.0)
    assert r.ok is False
    assert len(r.failures) == 1
    assert "t_live" in r.failures[0]
    assert "300" in r.failures[0]


def test_multiple_failures():
    """All three gates fail → three failure strings."""
    r = validate_background(F=20.0, F_ref=10.0, sum_Y=500.0, t_live=100.0)
    assert r.ok is False
    assert len(r.failures) == 3
    joined = " | ".join(r.failures).lower()
    assert "rate" in joined
    assert "sum_counts" in joined
    assert "t_live" in joined


def test_zero_ref_skips_ratio():
    """F_ref = 0 → ratio check is skipped, other gates still apply."""
    # Sum and live_time both OK → overall ok=True, ratio = NaN
    r = validate_background(F=10.0, F_ref=0.0, sum_Y=5000.0, t_live=1200.0)
    assert r.ok is True
    assert math.isnan(r.rate_ratio)
    assert r.failures == ()

    # Low counts still trips the counts gate even with F_ref=0
    r2 = validate_background(F=10.0, F_ref=0.0, sum_Y=400.0, t_live=1200.0)
    assert r2.ok is False
    assert len(r2.failures) == 1
    assert "sum_counts" in r2.failures[0]
    assert math.isnan(r2.rate_ratio)


def test_custom_thresholds_loosened():
    """Loosening thresholds turns previously-failing inputs into pass."""
    # Default would fail (rate diff 0.15, counts 800, t_live 500),
    # but with relaxed thresholds everything passes.
    r = validate_background(
        F=11.5, F_ref=10.0,
        sum_Y=800.0, t_live=500.0,
        rate_tolerance=0.2,
        min_sum_counts=500.0,
        min_live_time_s=300.0,
    )
    assert r.ok is True
    assert r.failures == ()


def test_frozen_dataclass_immutability():
    """BgControlResult is frozen — mutation must raise."""
    r = validate_background(F=10.0, F_ref=9.5, sum_Y=5000.0, t_live=1200.0)
    with pytest.raises((AttributeError, Exception)):
        r.ok = False  # type: ignore[misc]


# ════════════════════════════════════════════════════════════════════
# High-level `check_bg_quality(sample_spec, bg_spec)` — Hybrid layer
# Per AGENT_A_BRIEF_F-243.md (Agent D decision 2026-06-02)
# ════════════════════════════════════════════════════════════════════


def test_three_gates_pass_on_clean_pair():
    """Clean Th-232-like sample+bg pair: all three gates pass.

    Setup:
      sample.live_time = bg.live_time = 1200 s ⇒ ratio 1.0  (Gate 3 ✔)
      sum_bg = 5000 counts ⇒ σ = 1/sqrt(5000) ≈ 1.4 %      (Gate 2 ✔)
      rate_sample ≈ rate_bg (5% drift) ⇒ flux_drift 5%      (Gate 1 ✔)
    """
    sample = _mk_spec([10] * 600, 1200.0)   # sum=6000, rate=5.0
    bg     = _mk_spec([8.75] * 600, 1200.0) # sum=5250, rate≈4.375
    # flux_drift = |5.0 - 4.375| / 4.375 ≈ 0.143 — too high, redo:
    # actually we want clean ⇒ make rates equal-ish.
    sample = _mk_spec([5] * 1000, 1200.0)   # sum=5000, rate≈4.166
    bg     = _mk_spec([5] * 1000, 1200.0)   # sum=5000, rate≈4.166 ⇒ drift 0
    r = check_bg_quality(sample, bg)
    assert isinstance(r, BgQualityReport)
    assert r.passed is True
    assert all(g["passed"] for g in r.gates.values()), r.gates
    assert r.notes == []
    assert set(r.gates.keys()) == {"flux_drift", "sum_y_stat", "live_time_min"}


def test_flux_drift_triggers_fail():
    """Synthetic bg with +20% rate ⇒ Gate 1 (flux_drift) trips, others OK."""
    sample = _mk_spec([10] * 1000, 1000.0)  # sum=10000, rate=10.0
    # rate_bg = 8.0 ⇒ flux_drift = |10-8|/8 = 0.25 > 0.10
    bg = _mk_spec([8] * 1000, 1000.0)       # sum=8000,  rate=8.0
    r = check_bg_quality(sample, bg)
    assert r.passed is False
    assert r.gates["flux_drift"]["passed"] is False
    assert r.gates["sum_y_stat"]["passed"] is True   # N=8000, σ=1.1%
    assert r.gates["live_time_min"]["passed"] is True  # ratio = 1.0
    assert r.gates["flux_drift"]["value"] > 0.10
    assert any("drift" in n.lower() for n in r.notes)


def test_short_bg_triggers_gate3():
    """bg.live_time = 10% of sample ⇒ Gate 3 trips, others OK."""
    sample = _mk_spec([10] * 1000, 1000.0)  # sum=10000, rate=10.0
    # Same rate (10.0) so flux_drift = 0; large enough N for σ ≤ 10%.
    bg = _mk_spec([10] * 100, 100.0)        # sum=1000, rate=10.0
    r = check_bg_quality(sample, bg)
    assert r.passed is False
    assert r.gates["live_time_min"]["passed"] is False
    assert r.gates["flux_drift"]["passed"] is True
    assert r.gates["sum_y_stat"]["passed"] is True
    assert r.gates["live_time_min"]["value"] == pytest.approx(0.10, abs=1e-9)
    assert any("live_time" in n.lower() for n in r.notes)


def test_low_statistics_triggers_gate2():
    """bg with N=50 ⇒ σ = 1/√50 ≈ 14.1% > 10% ⇒ Gate 2 trips."""
    # Match rates to keep Gate 1 happy; same live_time to keep Gate 3 happy.
    sample = _mk_spec([1] * 50, 1000.0)     # sum=50, rate=0.05
    bg     = _mk_spec([1] * 50, 1000.0)     # sum=50, rate=0.05 ⇒ drift 0
    r = check_bg_quality(sample, bg)
    assert r.passed is False
    assert r.gates["sum_y_stat"]["passed"] is False
    assert r.gates["flux_drift"]["passed"] is True
    assert r.gates["live_time_min"]["passed"] is True
    expected_sigma = 1.0 / math.sqrt(50)
    assert r.gates["sum_y_stat"]["value"] == pytest.approx(expected_sigma)
    assert any(("σ" in n) or ("poisson" in n.lower()) for n in r.notes)


# ════════════════════════════════════════════════════════════════════
# F-243.1 — ROI-windowed flux drift (v1.18.30)
# ════════════════════════════════════════════════════════════════════


def _mk_spec_with_cal(counts, live_time, energy_cal=(0.0, 1.0)):
    """Spectrum stand-in with energy calibration: E(ch) = a0 + a1*ch + ..."""
    return SimpleNamespace(
        counts=list(counts),
        live_time=float(live_time),
        energy_cal=energy_cal,
    )


def test_f243_1_roi_pass():
    """ROI [600,700 keV], ROI drift ≈5.3% < 10% → ROI gate and overall PASS.

    Linear calibration E(ch)=ch keV; 1024 channels.
    sample ROI: 101 ch × 100 cts = 10100; bg ROI: 101 ch × 95 cts = 9595.
    drift = |10100−9595| / 9595 ≈ 5.3 % < 10 % → PASS.
    Gross gates also pass (drift ≈4.8%, σ≈1.0%, ratio=1.0).
    """
    sample_counts = [100 if 600 <= i <= 700 else 1 for i in range(1024)]
    bg_counts     = [95  if 600 <= i <= 700 else 1 for i in range(1024)]
    sample = _mk_spec_with_cal(sample_counts, 1200.0)
    bg     = _mk_spec_with_cal(bg_counts,     1200.0)
    r = check_bg_quality(sample, bg, roi_keV_list=[(600.0, 700.0)])
    assert isinstance(r, BgQualityReport)
    assert r.passed is True
    assert "roi_600_700" in r.gates
    assert r.gates["roi_600_700"]["passed"] is True
    assert r.notes == []


def test_f243_1_roi_fail():
    """ROI [600,700 keV], ROI drift ≈25% > 10% → ROI gate FAIL, report FAIL.

    bg ROI: 101 ch × 80 cts = 8080 vs sample 10100.
    drift = |10100−8080| / 8080 ≈ 25 % > 10 % → FAIL.
    """
    sample_counts = [100 if 600 <= i <= 700 else 1 for i in range(1024)]
    bg_counts     = [80  if 600 <= i <= 700 else 1 for i in range(1024)]
    sample = _mk_spec_with_cal(sample_counts, 1200.0)
    bg     = _mk_spec_with_cal(bg_counts,     1200.0)
    r = check_bg_quality(sample, bg, roi_keV_list=[(600.0, 700.0)])
    assert r.passed is False
    assert "roi_600_700" in r.gates
    assert r.gates["roi_600_700"]["passed"] is False
    assert any("ROI" in n for n in r.notes)


def test_f243_1_empty_roi_fallback():
    """roi_keV_list=[] and roi_keV_list=None → bit-identical to v1.18.29 path.

    Uses same clean pair as test_three_gates_pass_on_clean_pair to verify
    the result is unchanged when no ROI windows are requested.
    """
    sample = _mk_spec([5] * 1000, 1200.0)
    bg     = _mk_spec([5] * 1000, 1200.0)
    r_base  = check_bg_quality(sample, bg)
    r_empty = check_bg_quality(sample, bg, roi_keV_list=[])
    r_none  = check_bg_quality(sample, bg, roi_keV_list=None)
    assert r_empty.passed == r_base.passed
    assert r_empty.gates  == r_base.gates
    assert r_empty.notes  == r_base.notes
    assert r_none.passed  == r_base.passed
    assert r_none.gates   == r_base.gates
    assert r_none.notes   == r_base.notes


def test_f243_1_f_ref_zero_guard():
    """F_ref_roi=0 → rate check skipped (NaN ratio), no DivisionByZero.

    sum_Y_roi=2000 ≥ 1000 and t_live=1200 ≥ 600 → ok=True despite NaN ratio.
    """
    r = validate_background_roi(
        F_roi=5.0,
        F_ref_roi=0.0,
        sum_Y_roi=2000.0,
        t_live=1200.0,
    )
    assert isinstance(r, BgControlResult)
    assert r.ok is True
    assert math.isnan(r.rate_ratio)
    assert r.failures == ()


def test_f243_1_boundary_10pct():
    """Drift exactly 0.10 → PASS (convention: diff ≤ rate_tolerance passes).

    |F_roi − F_ref_roi| / F_ref_roi = |11 − 10| / 10 = 0.10 = rate_tolerance.
    Gate condition: diff > tolerance → fail; 0.10 > 0.10 is False → PASS.
    """
    r = validate_background_roi(
        F_roi=11.0,
        F_ref_roi=10.0,
        sum_Y_roi=2000.0,
        t_live=1200.0,
    )
    assert r.ok is True
    assert r.failures == ()
    assert math.isclose(
        abs(11.0 - 10.0) / 10.0, 0.10, rel_tol=1e-12
    ), "pre-condition: drift must be exactly 0.10"


# ════════════════════════════════════════════════════════════════════
# BUG-35 — Poisson |z|-test for bg stability (RAG-022 placeholder)
# Per ISO 11929-2:2019 §6: z = (B1 − B2) / √(B1 + B2)
# Methodology v2 §criterion 5 tiers: <2 stable, [2,3) borderline, ≥3 reject
# ════════════════════════════════════════════════════════════════════


def test_bug35_z_test_stable_equal_counts():
    """B1 == B2 → z=0 → tier='stable', PASS.

    Trivial case: two identical Poisson realisations differ by 0.
    """
    r = bg_z_test(1000.0, 1000.0)
    assert isinstance(r, BgZTestResult)
    assert r.z == 0.0
    assert r.abs_z == 0.0
    assert r.tier == "stable"
    assert r.passed is True
    assert r.B1 == 1000
    assert r.B2 == 1000
    assert r.note == ""


def test_bug35_z_test_stable_small_drift():
    """B1=1050, B2=1000 → z = 50/√2050 ≈ 1.10 < 2 → stable.

    Per ISO 11929-2:2019 §6 formula:
        z = (B1 - B2) / sqrt(B1 + B2)
        = 50 / sqrt(2050)
        ≈ 1.1043
    Within |z| < 2.0 stable tier per methodology v2 criterion 5.
    """
    r = bg_z_test(1050.0, 1000.0)
    expected_z = 50.0 / math.sqrt(2050.0)
    assert math.isclose(r.z, expected_z, rel_tol=1e-12)
    assert r.abs_z < Z_TIER_STABLE_MAX
    assert r.tier == "stable"
    assert r.passed is True


def test_bug35_z_test_borderline_tier():
    """|z| in [2, 3) → tier='borderline', soft PASS.

    Construct B1=1100, B2=1000 → z=100/√2100 ≈ 2.18 ∈ [2, 3).
    Per methodology v2 criterion 5: borderline = recount recommended,
    but still PASS (not a hard reject).
    """
    r = bg_z_test(1100.0, 1000.0)
    expected_z = 100.0 / math.sqrt(2100.0)
    assert math.isclose(r.z, expected_z, rel_tol=1e-12)
    assert Z_TIER_STABLE_MAX <= r.abs_z < Z_TIER_BORDERLINE_MAX
    assert r.tier == "borderline"
    assert r.passed is True  # soft PASS per methodology v2
    assert "borderline" in r.note.lower()


def test_bug35_z_test_reject_tier():
    """|z| ≥ 3 → tier='reject', FAIL.

    Construct B1=1200, B2=1000 → z=200/√2200 ≈ 4.26 > 3.
    Per ISO 11929-2:2019 §6: 'exceedance >3·σ indicates non-stationarity'.
    """
    r = bg_z_test(1200.0, 1000.0)
    expected_z = 200.0 / math.sqrt(2200.0)
    assert math.isclose(r.z, expected_z, rel_tol=1e-12)
    assert r.abs_z >= Z_TIER_BORDERLINE_MAX
    assert r.tier == "reject"
    assert r.passed is False
    assert "non-stationary" in r.note.lower() or "reject" in r.note.lower()


def test_bug35_z_test_signed_z_negative_when_b1_lt_b2():
    """B1 < B2 → signed z is negative; abs_z drives tier decision."""
    r = bg_z_test(900.0, 1000.0)
    expected_z = -100.0 / math.sqrt(1900.0)
    assert math.isclose(r.z, expected_z, rel_tol=1e-12)
    assert r.z < 0.0
    assert r.abs_z == abs(r.z)
    # |z| ≈ 2.29 → borderline
    assert r.tier == "borderline"


def test_bug35_z_test_zero_counts_undefined():
    """B1=B2=0 → undefined, FAIL (do not silently treat as stable).

    Empty spectrum or zeroed bg is a programming error upstream — z-test
    is mathematically undefined when B1+B2 = 0.
    """
    r = bg_z_test(0.0, 0.0)
    assert math.isnan(r.z)
    assert math.isnan(r.abs_z)
    assert r.tier == "undefined"
    assert r.passed is False
    assert "undefined" in r.note.lower()


def test_bug35_z_test_self_correcting_strong_bg():
    """At strong bg, large absolute drift can still be tier='stable'.

    B1+B2 = 2_000_000; absolute diff = 1000 counts.
    z = 1000 / √2_000_000 ≈ 0.707 < 2 → stable.

    The same absolute diff (1000) at weak bg (B1+B2=2000) gives
    z ≈ 22.4 → reject. This is the self-correction property that
    motivates the supersession of the 10% engineering rule.
    """
    # Strong bg path
    r_strong = bg_z_test(1_000_500.0, 999_500.0)
    expected_z_strong = 1000.0 / math.sqrt(2_000_000.0)
    assert math.isclose(r_strong.z, expected_z_strong, rel_tol=1e-12)
    assert r_strong.tier == "stable"
    assert r_strong.passed is True

    # Weak bg path — same absolute diff, very different z
    r_weak = bg_z_test(1500.0, 500.0)
    # z = 1000 / sqrt(2000) ≈ 22.36 → reject
    assert r_weak.abs_z > Z_TIER_BORDERLINE_MAX
    assert r_weak.tier == "reject"
    assert r_weak.passed is False


def test_bug35_z_test_boundary_z_equals_2_is_borderline():
    """|z| == 2.0 exactly → borderline tier (boundary inclusive on lower).

    Convention: tier = 'stable' if abs_z < 2.0 (strict), borderline if
    2.0 ≤ abs_z < 3.0. Mirrors the conservative 'recount recommended'
    reading of methodology v2 §criterion 5.

    Construct: B1 + B2 = 100; |z| = 2 → |B1 − B2| = 20.
    Use B1 = 60, B2 = 40 → z = 20/√100 = 2.0 exactly.
    """
    r = bg_z_test(60.0, 40.0)
    assert math.isclose(r.z, 2.0, rel_tol=1e-12)
    assert r.tier == "borderline"
    assert r.passed is True


def test_bug35_check_bg_quality_no_z_arg_backwards_compat():
    """peak_z_roi_keV_list=None / [] → no z-gates added, identical result.

    Critical backwards-compat guarantee: existing call sites without the
    new parameter must produce bit-identical reports to v1.18.30.
    """
    sample = _mk_spec([5] * 1000, 1200.0)
    bg     = _mk_spec([5] * 1000, 1200.0)
    r_base  = check_bg_quality(sample, bg)
    r_none  = check_bg_quality(sample, bg, peak_z_roi_keV_list=None)
    r_empty = check_bg_quality(sample, bg, peak_z_roi_keV_list=[])
    assert r_none.passed == r_base.passed
    assert r_none.gates  == r_base.gates
    assert r_none.notes  == r_base.notes
    assert r_empty.passed == r_base.passed
    assert r_empty.gates  == r_base.gates
    assert r_empty.notes  == r_base.notes
    # No z_roi_* keys should appear
    assert not any(k.startswith("z_roi_") for k in r_base.gates)


def test_bug35_check_bg_quality_z_roi_pass():
    """ROI [600,700 keV] with B1≈B2 → z≈0 → z_roi gate PASS, overall PASS.

    Linear calibration E(ch)=ch keV; 1024 channels.
    Sample ROI counts = bg ROI counts = 101 × 100 = 10100 each.
    z = 0/√20200 = 0 → tier='stable' → PASS.
    """
    sample_counts = [100 if 600 <= i <= 700 else 1 for i in range(1024)]
    bg_counts     = [100 if 600 <= i <= 700 else 1 for i in range(1024)]
    sample = _mk_spec_with_cal(sample_counts, 1200.0)
    bg     = _mk_spec_with_cal(bg_counts,     1200.0)
    r = check_bg_quality(sample, bg, peak_z_roi_keV_list=[(600.0, 700.0)])
    assert r.passed is True
    assert "z_roi_600_700" in r.gates
    g = r.gates["z_roi_600_700"]
    assert g["passed"] is True
    assert g["tier"] == "stable"
    assert math.isclose(g["value"], 0.0, abs_tol=1e-12)
    assert g["B1"] == g["B2"] == 10100


def test_bug35_check_bg_quality_z_roi_reject():
    """ROI [600,700 keV] with strong drift → z ≥ 3 → reject FAIL.

    B1 = 101 × 100 = 10100 (sample), B2 = 101 × 50 = 5050 (bg).
    z = 5050 / √15150 ≈ 41.04 → tier='reject' → overall FAIL.
    Note: this case the 10% engineering rule also fails (drift 100%),
    but the z-test gives a self-correcting σ-quantified diagnostic.
    """
    sample_counts = [100 if 600 <= i <= 700 else 1 for i in range(1024)]
    bg_counts     = [50  if 600 <= i <= 700 else 1 for i in range(1024)]
    sample = _mk_spec_with_cal(sample_counts, 1200.0)
    bg     = _mk_spec_with_cal(bg_counts,     1200.0)
    r = check_bg_quality(sample, bg, peak_z_roi_keV_list=[(600.0, 700.0)])
    assert r.passed is False
    g = r.gates["z_roi_600_700"]
    assert g["passed"] is False
    assert g["tier"] == "reject"
    assert g["abs_z"] >= Z_TIER_BORDERLINE_MAX
    # Note explicitly references ISO §6
    assert any(
        "z-test" in n.lower() and "ISO 11929-2".lower() in n.lower()
        for n in r.notes
    )


def test_bug35_check_bg_quality_z_requires_energy_cal():
    """peak_z_roi_keV_list without energy_cal on specs → ValueError.

    Mirrors the existing roi_keV_list contract (RAG-008): per-window
    operations require calibrated energy axis to map keV → channels.
    """
    sample = _mk_spec([5] * 1000, 1200.0)   # no energy_cal
    bg     = _mk_spec([5] * 1000, 1200.0)
    with pytest.raises(ValueError, match="peak_z_roi_keV_list requires energy_cal"):
        check_bg_quality(sample, bg, peak_z_roi_keV_list=[(600.0, 700.0)])
