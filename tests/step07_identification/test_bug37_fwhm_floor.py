"""BUG-37 / RAG-043 — FWHM² floor reduction from 1.0 → 0.01 keV² in
``fwhm_keV_at_energy`` (effective FWHM floor 1.0 → 0.1 keV).

Rationale (physics):
  The old `max(val, 1.0)` floor on FWHM² silently clamped FWHM to 1.0 keV
  whenever a poorly-constrained `lsrm_peaks_table` quadratic fit
  extrapolated to FWHM² < 1 keV² at low E. This masked physically
  unreasonable model output (any real NaI(Tl) scintillator gives
  FWHM ≥ 6 keV at E ≥ 50 keV; HPGe gives FWHM ≈ 0.5-1.5 keV) and could
  cause peak-search windows to under-cover Am-241 @ 59.5 keV when the
  identification window scales with FWHM.

Lowered floor 0.01 keV² (= 0.1 keV FWHM) is a pure numerical-safety net
for degenerate model fits (val ≤ 0). No physical detector hits this floor
under any realistic configuration.

Tests:
  1. Am-241 @ 59.5 keV: post-floor FWHM matches model output (not clamped
     to 1 keV) for a synthetic NaI-like model and for the default model.
  2. Numerical safety: with degenerate model (val=0), floor returns
     finite positive FWHM (no NaN/zero, no division-by-zero downstream).
  3. Cs-137 @ 661.7 keV: default model FWHM unchanged by floor change.
  4. Floor activation: with an extreme degenerate model that produces
     val < 0.01 keV², output is exactly sqrt(0.01) = 0.1 keV.
"""
from __future__ import annotations

import math
import os
import sys

import pytest

SCRIPTS = os.path.join(
    os.path.dirname(__file__), os.pardir, os.pardir, "scripts",
)
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)

from gamma.identification.staged_pipeline import (  # noqa: E402
    _DEFAULT_NAI_FWHM_MODEL,
    fwhm_keV_at_energy,
)


# ─────────────────────────────────────────────────────────────────────
# Test 1 — Am-241 @ 59.5 keV: model output not clamped to 1.0 keV
# ─────────────────────────────────────────────────────────────────────

def test_fwhm_floor_does_not_kill_am241_59keV_default_NaI():
    """Default NaI 63×63 model at 59.5 keV: FWHM ≈ 13 keV, far above
    floor=0.1 keV. The 1.0 floor never triggered here, but we confirm
    the fix preserves correct output.
    """
    f = fwhm_keV_at_energy(_DEFAULT_NAI_FWHM_MODEL, 59.5)
    assert 10.0 < f < 18.0, (
        f"Default NaI FWHM(59.5)={f:.3f} keV outside expected 10-18 keV. "
        f"Either model regressed or floor masked output."
    )


def test_fwhm_floor_does_not_kill_am241_59keV_lsrm_alpha_sqrt_E():
    """α·√E model anchored at single low-E point. Replicates the
    `lsrm_peaks_table` 1-anchor fallback path. At 59.5 keV the floor
    must NOT clamp output to 1.0 keV.
    """
    # α=1.5 → FWHM(E) = 1.5·√E; at E=59.5 → FWHM ≈ 11.57 keV.
    alpha = 1.5
    model = (0.0, alpha ** 2, 0.0)
    f = fwhm_keV_at_energy(model, 59.5)
    expected = alpha * math.sqrt(59.5)
    assert f == pytest.approx(expected, rel=1e-6), (
        f"α·√E model with α={alpha}: FWHM(59.5)={f:.4f}, "
        f"expected {expected:.4f}. Floor leaking into normal output."
    )
    assert f > 5.0, (
        f"α·√E FWHM(59.5)={f:.3f} ≤ 5 keV — floor masking output."
    )


# ─────────────────────────────────────────────────────────────────────
# Test 2 — numerical safety: degenerate model → finite positive
# ─────────────────────────────────────────────────────────────────────

def test_fwhm_floor_numerical_safety_zero_model():
    """Degenerate (0,0,0) model → val=0 → without floor sqrt(0)=0
    → division-by-zero downstream. Floor must return finite >0.
    """
    f = fwhm_keV_at_energy((0.0, 0.0, 0.0), 100.0)
    assert math.isfinite(f), f"FWHM not finite for degenerate model: {f}"
    assert f > 0, f"FWHM={f} ≤ 0 for degenerate model"
    # Floor active → expect exactly 0.1 keV (sqrt(0.01))
    assert f == pytest.approx(0.1, abs=1e-9), (
        f"Degenerate model should hit floor 0.1 keV, got {f:.6f}"
    )


def test_fwhm_floor_numerical_safety_negative_model():
    """Negative-output model (unphysical extrapolation) → val<0
    must still return floor, not NaN."""
    # Model that yields negative at low E
    f = fwhm_keV_at_energy((-100.0, 0.0, 0.0), 10.0)
    assert math.isfinite(f), f"FWHM not finite for negative model: {f}"
    assert f == pytest.approx(0.1, abs=1e-9), (
        f"Negative-val model should hit floor 0.1 keV, got {f:.6f}"
    )


# ─────────────────────────────────────────────────────────────────────
# Test 3 — Cs-137 @ 661.7 keV calibration unchanged
# ─────────────────────────────────────────────────────────────────────

def test_fwhm_floor_preserves_661_calibration():
    """Default NaI model at Cs-137 line: FWHM ≈ 47 keV, no floor effect.
    Floor change from 1.0 → 0.01 must not perturb anchor-line output.
    """
    f = fwhm_keV_at_energy(_DEFAULT_NAI_FWHM_MODEL, 661.66)
    assert 45.0 <= f <= 50.0, (
        f"FWHM(661.66)={f:.2f} keV outside 45-50 keV band — "
        f"floor change broke calibration."
    )


# ─────────────────────────────────────────────────────────────────────
# Test 4 — floor activation: val < 0.01 keV² returns exactly 0.1 keV
# ─────────────────────────────────────────────────────────────────────

def test_fwhm_floor_low_E_floor_active_at_0p01_keV2():
    """Construct model that explicitly gives val < 0.01 keV² to confirm
    new floor value. Old floor 1.0 keV² → FWHM=1.0; new 0.01 → FWHM=0.1.
    """
    # Constant val=0.005 keV² (below new floor 0.01, well below old 1.0).
    # Use model (a=0.005, b=0, c=0). At any E ≥ 5 keV, val = 0.005.
    model = (0.005, 0.0, 0.0)
    f = fwhm_keV_at_energy(model, 100.0)
    # New floor active → returns sqrt(0.01) = 0.1 (NOT sqrt(0.005) ≈ 0.0707)
    assert f == pytest.approx(0.1, abs=1e-9), (
        f"Expected floor-clamped 0.1 keV, got {f:.6f}. "
        f"Either floor not active or floor value wrong."
    )
    # And confirm it's also NOT the old 1.0 floor (regression guard)
    assert f < 0.5, (
        f"FWHM={f} ≥ 0.5 keV — old 1.0 floor may have leaked back in."
    )


def test_fwhm_floor_above_floor_passes_through():
    """val just above floor (0.05 keV²) → output = sqrt(0.05) ≈ 0.2236,
    NOT clamped. Confirms floor is bottom-only, not always-applied.
    """
    model = (0.05, 0.0, 0.0)
    f = fwhm_keV_at_energy(model, 100.0)
    expected = math.sqrt(0.05)
    assert f == pytest.approx(expected, rel=1e-9), (
        f"val=0.05 (above floor 0.01) should pass through unchanged. "
        f"Got {f:.6f}, expected {expected:.6f}."
    )
