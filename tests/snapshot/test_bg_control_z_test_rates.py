"""Wave 4 (2026-06-04) — A territory.

Regression guard for ``bg_z_test_rates(c1, t1, c2, t2)`` — rate-normalised
two-sample z-test for unequal Poisson live-times.

Source: BUG-35 z-test followup #3 (deferred from wave 1) —
``_state/agent_a/outbox/2026-06-04_backlog_top1_bug35_z_test.md``.

Formula (Gilmore & Joss §5.5):
    R_i = c_i / t_i,   Var(R_i) = c_i / t_i²
    z   = (R1 − R2) / √( Var(R1) + Var(R2) )
        = (c1/t1 − c2/t2) / √( c1/t1² + c2/t2² )

When t1 == t2 == t the rate-form reduces algebraically to the
integer-count form already in use:
    z = (c1 − c2) / √(c1 + c2)

Cite-list:
    * F-157 (LSRM > Будыка > Gilmore — Gilmore canonical for §5.5)
    * Gilmore & Joss §5.5 (eq. 5.21 — comparison of two count rates)
    * ISO 11929-2:2019 §6 (Poisson propagation foundation, same as
      ``bg_z_test``)
    * Methodology v2 §criterion 5 (3-sigma reject tier shared with
      integer-count form)
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from gamma.io.bg_control import (  # noqa: E402
    Z_TIER_BORDERLINE_MAX,
    bg_z_test,
    bg_z_test_rates,
)


# ─────────────────────────────────────────────────────────────────────
# 1) Equal-time invariant: matches the integer-count form exactly
# ─────────────────────────────────────────────────────────────────────

def test_equal_time_matches_integer_form():
    """When t1 == t2 == t, the rate-form reduces algebraically to
    ``(c1 − c2) / √(c1 + c2)``. Pin numerical equivalence with the
    existing ``bg_z_test`` kernel — guards against silent divergence
    of the two APIs."""
    cases = [(1000, 1100), (500, 600), (10_000, 9_500), (100, 150)]
    t = 3600.0
    for c1, c2 in cases:
        z_rates, _ = bg_z_test_rates(c1, t, c2, t)
        z_int = bg_z_test(c1, c2).z
        assert math.isfinite(z_rates) and math.isfinite(z_int)
        assert abs(z_rates - z_int) < 1e-12, (
            f"equal-time divergence at (c1={c1}, c2={c2}): "
            f"rates={z_rates}, integer={z_int}"
        )


# ─────────────────────────────────────────────────────────────────────
# 2) Unequal time: zero z when rates are identical
# ─────────────────────────────────────────────────────────────────────

def test_unequal_time_identical_rates_yields_zero_z():
    """1000 counts / 100 s = 10 cps; 10 000 counts / 1000 s = 10 cps.
    Same rate, vastly different live-times → z = 0. This is the
    pathology the legacy integer-form misdiagnoses: it would yield
    ``(1000 − 10 000) / √11 000`` ≈ −85.8 (massive false reject)
    because it ignores live-time."""
    z, sig = bg_z_test_rates(1000.0, 100.0, 10_000.0, 1000.0)
    assert math.isfinite(z)
    assert abs(z) < 1e-9, f"identical rates should give z≈0, got {z}"
    assert sig is False


# ─────────────────────────────────────────────────────────────────────
# 3) Unequal time: known-significant case
# ─────────────────────────────────────────────────────────────────────

def test_unequal_time_clearly_significant_case():
    """R1=10 cps from 1000 c / 100 s; R2=20 cps from 20 000 c / 1000 s.
    Var(R1) = 1000 / 10 000 = 0.1; Var(R2) = 20 000 / 1e6 = 0.02.
    z = (10 − 20) / √0.12 ≈ −28.87 → ``|z|`` ≫ 3 → significant."""
    z, sig = bg_z_test_rates(1000.0, 100.0, 20_000.0, 1000.0)
    assert math.isfinite(z)
    expected = (10.0 - 20.0) / math.sqrt(1000.0 / 100.0**2 + 20_000.0 / 1000.0**2)
    assert abs(z - expected) < 1e-9, f"formula drift: got {z}, expected {expected}"
    assert sig is True
    assert abs(z) > Z_TIER_BORDERLINE_MAX


# ─────────────────────────────────────────────────────────────────────
# 4) Sign direction is preserved
# ─────────────────────────────────────────────────────────────────────

def test_sign_direction_preserved():
    """Swapping (c1,t1) ↔ (c2,t2) must flip the sign of z exactly,
    not change its magnitude. (Pin to catch accidental abs() inside
    the kernel — the signed form is needed for one-sided BG-rising
    vs BG-falling diagnostics.)"""
    z_pos, _ = bg_z_test_rates(2000.0, 100.0, 1000.0, 100.0)
    z_neg, _ = bg_z_test_rates(1000.0, 100.0, 2000.0, 100.0)
    assert z_pos > 0 and z_neg < 0
    assert abs(z_pos + z_neg) < 1e-12, f"sign asymmetry: {z_pos} + {z_neg} != 0"


# ─────────────────────────────────────────────────────────────────────
# 5) Three-sigma boundary uses ``> 3.0`` (strict)
# ─────────────────────────────────────────────────────────────────────

def test_three_sigma_boundary_strict():
    """``is_significant`` is True iff ``|z| > 3.0`` (strict). Pin the
    boundary direction so it cannot silently flip to ``≥``."""
    # Construct a case where z is well below 3 → not significant.
    z_low, sig_low = bg_z_test_rates(1000.0, 100.0, 1100.0, 100.0)
    assert math.isfinite(z_low)
    # |z| = 100 / sqrt(2100) ≈ 2.18 → below 3.0
    assert abs(z_low) < Z_TIER_BORDERLINE_MAX
    assert sig_low is False

    # And a case where z is well above 3 → significant.
    z_hi, sig_hi = bg_z_test_rates(1000.0, 100.0, 1500.0, 100.0)
    # |z| = 500 / sqrt(2500) = 10 → far above 3
    assert abs(z_hi) > Z_TIER_BORDERLINE_MAX
    assert sig_hi is True


# ─────────────────────────────────────────────────────────────────────
# 6) Degenerate live-time → NaN, not significant
# ─────────────────────────────────────────────────────────────────────

def test_nonpositive_live_time_returns_nan():
    """t ≤ 0 is degenerate. Return (NaN, False) — never raise; the
    caller may legitimately receive zero-live-time spectra from a
    badly truncated file and we want graceful downstream."""
    for t1, t2 in [(0.0, 100.0), (100.0, 0.0), (-1.0, 100.0), (100.0, -1.0)]:
        z, sig = bg_z_test_rates(1000.0, t1, 1000.0, t2)
        assert math.isnan(z), f"t1={t1}, t2={t2}: expected NaN, got {z}"
        assert sig is False


# ─────────────────────────────────────────────────────────────────────
# 7) Zero counts on both sides → NaN
# ─────────────────────────────────────────────────────────────────────

def test_zero_counts_returns_nan():
    """c1 + c2 ≤ 0 → variance undefined → NaN, not significant.
    Matches ``bg_z_test`` behaviour for the ``total <= 0`` branch."""
    z, sig = bg_z_test_rates(0.0, 100.0, 0.0, 100.0)
    assert math.isnan(z)
    assert sig is False


# ─────────────────────────────────────────────────────────────────────
# 8) Negative counts → fail-closed NaN
# ─────────────────────────────────────────────────────────────────────

def test_negative_counts_fail_closed():
    """Negative counts are nonsensical for Poisson and would make
    the variance term go negative under the sqrt. Fail closed: NaN,
    not significant, no exception."""
    z, sig = bg_z_test_rates(-10.0, 100.0, 100.0, 100.0)
    assert math.isnan(z)
    assert sig is False


# ─────────────────────────────────────────────────────────────────────
# 9) Float c_i accepted (ROI sums may be fractional after slicing)
# ─────────────────────────────────────────────────────────────────────

def test_float_counts_accepted_no_rounding():
    """Unlike ``bg_z_test`` (rounds to int for Poisson semantics),
    the rate form is variance-correct for any non-negative real
    ``c_i``. Pin: do NOT round."""
    # 1000.5 c / 100 s vs 1000.0 c / 100 s → small but nonzero z
    z_frac, _ = bg_z_test_rates(1000.5, 100.0, 1000.0, 100.0)
    z_int, _ = bg_z_test_rates(1000.0, 100.0, 1000.0, 100.0)
    assert math.isfinite(z_frac) and math.isfinite(z_int)
    assert z_int == 0.0
    assert z_frac != 0.0, "fractional input must produce a nonzero z"
    assert z_frac > 0


# ─────────────────────────────────────────────────────────────────────
# 10) Return shape is (float, bool)
# ─────────────────────────────────────────────────────────────────────

def test_return_shape_is_float_bool_tuple():
    """Signature contract: ``-> tuple[float, bool]``. Caller code
    unpacks as ``z, sig = bg_z_test_rates(...)``; a regression to
    ``(np.float64, np.bool_)`` would break operator-friendly JSON
    serialisation downstream."""
    z, sig = bg_z_test_rates(100.0, 10.0, 150.0, 10.0)
    assert isinstance(z, float)
    assert isinstance(sig, bool)
