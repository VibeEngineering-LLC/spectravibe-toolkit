"""Wave 4 (2026-06-04) — A territory.

Regression guard for the zero-point (a0) safety check added in
``scripts/gamma/calibration/bg_subtract_dual_mode.py``.

Background:
    Agent B wave 3 documented a silent BG-subtraction bug
    (``_state/agent_b/outbox/2026-06-04_wave3_k40_regen_and_bg_subtraction_investigation.md``,
    lines 112–211): source spectrum with energy_cal a0=+47.669 keV and
    background spectrum with a0=−25.814 keV → Δa0 = 73.5 keV. The two
    spectra had nearly identical gain (Δa1/a1 = 0.28% < 0.5%) so the
    pre-existing logic selected ``rate_normalized_channel`` mode, which
    subtracted bg-channel-i from src-channel-i even though those
    channels mapped to different energies. At the K-40 line this
    suppressed BG subtraction by ~50 %.

These tests pin the new ``ZERO_POINT_MATCH_THRESHOLD_KEV = 30 keV`` gate
that forces ``energy_aligned`` mode when ``|Δa0| > 30 keV``, regardless
of gain match.

Cite-list:
    * F-243 (BG subtraction safety family)
    * F-157 (LSRM > Будыка > Gilmore precedence)
    * Gilmore & Joss §6.4 (NaI(Tl) FWHM model — defends 30 keV threshold
      choice as ≈ 0.3 × NaI FWHM at 1460 keV ≈ 10 × HPGe FWHM)
    * ISO 11929-2:2019 §6 (counting-statistics propagation)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from gamma.calibration.bg_subtract_dual_mode import (  # noqa: E402
    GAIN_MATCH_THRESHOLD,
    ZERO_POINT_MATCH_THRESHOLD_KEV,
    subtract_background,
)
from gamma.spectrum import Spectrum  # noqa: E402


# ─────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────

def _make_spec(counts, live_time, *, a0=0.0, a1=1.0, source_path="test"):
    """Build a minimal Spectrum suitable for subtract_background()."""
    spec = Spectrum(
        counts=np.array(counts, dtype=np.int64),
        live_time=float(live_time),
        real_time=float(live_time),
        source_path=source_path,
        source_format="test",
    )
    spec.energy_cal = (float(a0), float(a1))
    spec.n_channels = len(counts)
    spec.n_channels_raw = len(counts)
    return spec


# 1024-channel synthetic spectrum, flat counts so the choice of mode is
# visible only in the result.mode label (not in the numeric output).
_FLAT_COUNTS = [10] * 1024


# ─────────────────────────────────────────────────────────────────────
# 1) Sentinel: zero a0 delta + gain match → rate_normalized_channel
# ─────────────────────────────────────────────────────────────────────

def test_zero_a0_delta_gain_match_uses_rate_normalized():
    """No a0 delta, no gain delta → preserve legacy fast path."""
    src = _make_spec(_FLAT_COUNTS, 100.0, a0=0.0, a1=3.0)
    bg = _make_spec(_FLAT_COUNTS, 100.0, a0=0.0, a1=3.0)
    res = subtract_background(src, bg, user_confirmed_applicable=True)
    assert res.mode == "rate_normalized_channel"


# ─────────────────────────────────────────────────────────────────────
# 2) THE B-wave-3 scenario: tiny gain delta, huge a0 delta → forced align
# ─────────────────────────────────────────────────────────────────────

def test_b_wave3_scenario_a0_73keV_forces_energy_aligned():
    """B wave 3 numbers verbatim.

    a0_src=+47.669, a1_src=2.99926; a0_bg=-25.814, a1_bg=2.990863.
    Δa1/a1 = |2.99926-2.990863|/2.99926 ≈ 0.28 % < 0.5 % → would have
    selected ``rate_normalized_channel`` before the fix. Δa0 = 73.5 keV
    ≫ 30 keV → safety gate must force ``energy_aligned``.
    """
    src = _make_spec(_FLAT_COUNTS, 100.0, a0=47.669, a1=2.99926)
    bg = _make_spec(_FLAT_COUNTS, 100.0, a0=-25.814, a1=2.990863)
    # Sanity: confirm gain match alone would have classified as
    # rate-normalised (regression guard against future GAIN_MATCH_THRESHOLD
    # tightening that would mask this test).
    gain_delta_rel = abs(src.energy_cal[1] - bg.energy_cal[1]) / src.energy_cal[1]
    assert gain_delta_rel < GAIN_MATCH_THRESHOLD, (
        f"Test premise broken: gain delta {gain_delta_rel:.4f} no longer "
        f"< GAIN_MATCH_THRESHOLD={GAIN_MATCH_THRESHOLD}"
    )
    # Sanity: confirm Δa0 is in the regime the gate is meant to catch.
    a0_delta = abs(src.energy_cal[0] - bg.energy_cal[0])
    assert a0_delta > ZERO_POINT_MATCH_THRESHOLD_KEV, (
        f"Test premise broken: Δa0={a0_delta:.2f} keV no longer above "
        f"ZERO_POINT_MATCH_THRESHOLD_KEV={ZERO_POINT_MATCH_THRESHOLD_KEV}"
    )
    res = subtract_background(src, bg, user_confirmed_applicable=True)
    assert res.mode == "energy_aligned", (
        f"B wave 3 regression: Δa0={a0_delta:.1f} keV must force "
        f"energy_aligned mode, got {res.mode!r}"
    )


# ─────────────────────────────────────────────────────────────────────
# 3) Below-threshold a0 drift (normal recal noise) → keep fast path
# ─────────────────────────────────────────────────────────────────────

def test_small_a0_drift_below_threshold_keeps_rate_normalized():
    """Δa0 = 5 keV (typical day-to-day NaI recal drift) ≪ 30 keV
    threshold — should NOT trigger the safety gate."""
    src = _make_spec(_FLAT_COUNTS, 100.0, a0=0.0, a1=3.0)
    bg = _make_spec(_FLAT_COUNTS, 100.0, a0=5.0, a1=3.0)
    res = subtract_background(src, bg, user_confirmed_applicable=True)
    assert res.mode == "rate_normalized_channel"


# ─────────────────────────────────────────────────────────────────────
# 4) Boundary: exactly at the threshold → still rate_normalized
# ─────────────────────────────────────────────────────────────────────

def test_a0_delta_exactly_at_threshold_keeps_rate_normalized():
    """Boundary semantics: the gate triggers on ``> 30 keV`` (strict),
    so Δa0 = 30.0 keV exactly stays in the fast path. This pins the
    chosen boundary direction so future refactors do not flip it."""
    src = _make_spec(_FLAT_COUNTS, 100.0, a0=0.0, a1=3.0)
    bg = _make_spec(_FLAT_COUNTS, 100.0,
                    a0=ZERO_POINT_MATCH_THRESHOLD_KEV, a1=3.0)
    res = subtract_background(src, bg, user_confirmed_applicable=True)
    assert res.mode == "rate_normalized_channel"


# ─────────────────────────────────────────────────────────────────────
# 5) Just above threshold → triggers safety gate
# ─────────────────────────────────────────────────────────────────────

def test_a0_delta_just_above_threshold_forces_energy_aligned():
    """Δa0 = 30.001 keV → safety gate fires."""
    src = _make_spec(_FLAT_COUNTS, 100.0, a0=0.0, a1=3.0)
    bg = _make_spec(_FLAT_COUNTS, 100.0,
                    a0=ZERO_POINT_MATCH_THRESHOLD_KEV + 1e-3, a1=3.0)
    res = subtract_background(src, bg, user_confirmed_applicable=True)
    assert res.mode == "energy_aligned"


# ─────────────────────────────────────────────────────────────────────
# 6) Sign-independence: negative a0 delta of the same magnitude triggers
# ─────────────────────────────────────────────────────────────────────

def test_a0_delta_sign_independent():
    """The gate uses ``abs(a0_src - a0_bg)``; sign of the difference
    must not matter. Swap which spectrum has the high a0."""
    src = _make_spec(_FLAT_COUNTS, 100.0, a0=-30.0, a1=3.0)
    bg = _make_spec(_FLAT_COUNTS, 100.0, a0=+30.0, a1=3.0)  # Δ=60 keV
    res = subtract_background(src, bg, user_confirmed_applicable=True)
    assert res.mode == "energy_aligned"


# ─────────────────────────────────────────────────────────────────────
# 7) force_mode overrides the safety gate (operator escape hatch)
# ─────────────────────────────────────────────────────────────────────

def test_force_mode_overrides_safety_gate():
    """``force_mode`` is documented as bypassing the auto-select.
    The safety gate must also be bypassable when the operator
    explicitly requests rate_normalized_channel — needed for
    differential debugging / regression bisection."""
    src = _make_spec(_FLAT_COUNTS, 100.0, a0=+47.669, a1=2.99926)
    bg = _make_spec(_FLAT_COUNTS, 100.0, a0=-25.814, a1=2.990863)
    res = subtract_background(
        src, bg,
        user_confirmed_applicable=True,
        force_mode="rate_normalized_channel",
    )
    assert res.mode == "rate_normalized_channel"


# ─────────────────────────────────────────────────────────────────────
# 8) Notes string surfaces Δa0 for diagnostics
# ─────────────────────────────────────────────────────────────────────

def test_notes_string_contains_a0_delta():
    """Notes string carries diagnostic Δa0 value (operator visibility).

    Pinned because B wave 3 had to compute this delta manually from
    raw energy_cal tuples — having it in the result.notes string makes
    the next anomaly easy to spot in JSON reports.
    """
    src = _make_spec(_FLAT_COUNTS, 100.0, a0=+47.669, a1=2.99926)
    bg = _make_spec(_FLAT_COUNTS, 100.0, a0=-25.814, a1=2.990863)
    res = subtract_background(src, bg, user_confirmed_applicable=True)
    assert "Δa₀" in res.notes
    assert "73.4" in res.notes or "73.5" in res.notes, (
        f"Notes string should include Δa0=73.48 keV: {res.notes!r}"
    )


# ─────────────────────────────────────────────────────────────────────
# 9) Threshold constant is exported and stable
# ─────────────────────────────────────────────────────────────────────

def test_threshold_constant_is_30_keV():
    """The threshold literal is API-visible (operators may want to
    inspect it). Pin to 30.0 — changes require updating the F-243
    documentation and the B wave 3 outbox cross-references."""
    assert ZERO_POINT_MATCH_THRESHOLD_KEV == 30.0
