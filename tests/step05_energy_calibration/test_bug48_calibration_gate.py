# -*- coding: utf-8 -*-
"""BUG-48 — Hard calibration self-consistency gate.

Unit-tests for ``gamma.calibration.calibration_gate.evaluate_calibration_gate``.

Each test constructs a minimal Spectrum-like object with only the fields the
gate inspects (``energy_cal``, ``n_channels``, ``stored_fwhm_calibration``,
``extras``). No real fixtures are read — this keeps the gate test deterministic
and decoupled from the actual reader / detector library.

Coverage matrix:
    G1 — energy-cal monotonicity:
        * test_G1_pass_linear_increasing
        * test_G1_hard_fail_negative_slope
        * test_G1_hard_fail_no_calibration
    G2 — channel-range sanity:
        * test_G2_hard_fail_truncated_range
    G3 — FWHM monotonicity (soft warning):
        * test_G3_soft_warn_decreasing_simplesqrt
    G4 — FWHM plausibility band (soft warning):
        * test_G4_soft_warn_nai_out_of_band
        * test_G4_pass_hpge_in_band
    Pipeline integration:
        * test_extras_populated_when_gate_runs_via_pipeline_hook

Citations: gate criteria sourced in module docstring; see RAG-043
(Gilmore §6.4 + LSRM-9.4 §3.2) for FWHM-band evidence.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from gamma.calibration.calibration_gate import (
    evaluate_calibration_gate,
    CalibrationGateResult,
    HPGE_BAND_R,
    SCINT_BAND_R,
    E_LOW_KEV_MIN,
    E_HIGH_KEV_MIN,
    STANDARD_FWHM_ANCHORS_KEV,
    MEASURED_PEAK_THRESHOLD,
    ANCHOR_OVERSHOOT_FACTOR,
)
from gamma.spectrum import StoredFwhmCalibration, FwhmCalPeak


# ---------------------------------------------------------------------------
# Test scaffolding — minimal Spectrum stand-in
# ---------------------------------------------------------------------------

class _SpecStub:
    """Just the fields the gate reads. Avoids importing the heavy reader."""

    def __init__(
        self,
        *,
        energy_cal=None,
        n_channels=8192,
        stored_fwhm_calibration=None,
    ):
        self.energy_cal = energy_cal
        self.n_channels = n_channels
        self.stored_fwhm_calibration = stored_fwhm_calibration
        self.extras = {}


# ---------------------------------------------------------------------------
# G1 — Energy-cal monotonicity
# ---------------------------------------------------------------------------

def test_G1_pass_linear_increasing():
    """Healthy NaI-like cal (a0=0, a1=0.4 keV/ch over 8192 ch → 0–3276 keV)."""
    spec = _SpecStub(energy_cal=(0.0, 0.4), n_channels=8192)
    r = evaluate_calibration_gate(spec)
    assert isinstance(r, CalibrationGateResult)
    assert r.passed, f"Expected pass, got hard_failures={r.hard_failures}"
    assert not r.hard_failures
    assert "G1" in r.criteria_evaluated
    assert "G2" in r.criteria_evaluated


def test_G1_hard_fail_negative_slope():
    """Inverted axis (a1 < 0) → must hard-fail G1."""
    spec = _SpecStub(energy_cal=(3000.0, -0.4), n_channels=8192)
    r = evaluate_calibration_gate(spec)
    assert not r.passed
    codes = [h["code"] for h in r.hard_failures]
    assert "G1" in codes
    # G2 will also flag because E(n_ch-1) = 3000 - 0.4*8191 ≈ -276
    # That's expected; the test only requires G1 to be raised.


def test_G1_hard_fail_no_calibration():
    """Missing energy_cal → hard fail with explicit G1 entry."""
    spec = _SpecStub(energy_cal=None, n_channels=4096)
    r = evaluate_calibration_gate(spec)
    assert not r.passed
    assert r.hard_failures[0]["code"] == "G1"
    assert "missing energy calibration" in r.reason.lower()


# ---------------------------------------------------------------------------
# G2 — Channel-range sanity
# ---------------------------------------------------------------------------

def test_G2_hard_fail_truncated_range():
    """Spectrum whose max E = 50 keV (below E_HIGH_KEV_MIN=100) → G2 fail."""
    # Linear 0–50 keV across 1024 channels: a1 = 50/1023
    spec = _SpecStub(energy_cal=(0.0, 50.0 / 1023.0), n_channels=1024)
    r = evaluate_calibration_gate(spec)
    assert not r.passed
    codes = [h["code"] for h in r.hard_failures]
    assert "G2" in codes
    # G1 should still pass (slope > 0)
    g1_failures = [h for h in r.hard_failures if h["code"] == "G1"]
    assert not g1_failures


# ---------------------------------------------------------------------------
# G3 — FWHM monotonicity (soft warning)
# ---------------------------------------------------------------------------

def test_G3_soft_warn_decreasing_simplesqrt():
    """SimpleSqrtFwhm with c1 < 0 (FWHM decreases with channel) → G3 soft."""
    sf = StoredFwhmCalibration(
        coefficients=(100.0, -0.05),     # c0=100, c1=-0.05 → non-physical
        model="SimpleSqrtFwhm",
        calibration_peaks=[
            FwhmCalPeak(channel=500, energy_keV=200.0, fwhm_channels=20.0),
        ],
    )
    spec = _SpecStub(
        energy_cal=(0.0, 0.4),
        n_channels=8192,
        stored_fwhm_calibration=sf,
    )
    r = evaluate_calibration_gate(spec)
    assert r.passed, "G3 is soft — must not block pass/fail verdict"
    soft_codes = [w["code"] for w in r.soft_warnings]
    assert "G3" in soft_codes
    # Ensure message references the non-physical coefficient
    g3 = next(w for w in r.soft_warnings if w["code"] == "G3")
    assert "c1=" in g3["message"]


# ---------------------------------------------------------------------------
# G4 — FWHM physical plausibility band (soft warning)
# ---------------------------------------------------------------------------

def test_G4_soft_warn_nai_out_of_band():
    """NaI-class cal peak with FWHM/E = 0.5% (below 2% scintillator floor)."""
    # E(ch=1500) = 1500 * 0.4 = 600 keV
    # FWHM 7.5 channels * 0.4 keV/ch = 3 keV → R = 0.5% (HPGe-band value)
    # But the second peak at 200 keV with FWHM 25 ch (10 keV → R=5%) makes
    # the gate classify as scintillator → first peak then flags.
    sf = StoredFwhmCalibration(
        coefficients=(),
        model="",
        calibration_peaks=[
            FwhmCalPeak(channel=500, energy_keV=200.0, fwhm_channels=25.0),
            FwhmCalPeak(channel=1500, energy_keV=600.0, fwhm_channels=7.5),
        ],
    )
    spec = _SpecStub(
        energy_cal=(0.0, 0.4),
        n_channels=8192,
        stored_fwhm_calibration=sf,
    )
    r = evaluate_calibration_gate(spec)
    assert r.passed
    # Class hint chosen from peak nearest 662 keV → 600 keV peak (R=0.5%)
    # which falls inside HPGe band → classified HPGe. Then the 200 keV
    # peak (R=5%, scintillator-range) falls OUTSIDE the HPGe band [0.05, 2]%.
    # Either direction, at least one G4 warning must fire.
    g4_warnings = [w for w in r.soft_warnings if w["code"] == "G4"]
    assert g4_warnings, (
        f"Expected G4 band warning; soft_warnings={r.soft_warnings}, "
        f"class hint={r.detector_class_hint}"
    )


def test_G4_pass_hpge_in_band():
    """HPGe-class cal peaks all inside [0.05%, 2%] band → no G4 warnings."""
    # HPGe: 0.3 keV/ch, peak at 1332.5 keV with FWHM 2 keV → R=0.15% (in band)
    # peak at 661.7 keV with FWHM 1.3 keV → R=0.20% (in band)
    sf = StoredFwhmCalibration(
        coefficients=(),
        model="",
        calibration_peaks=[
            FwhmCalPeak(
                channel=int(661.7 / 0.3),
                energy_keV=661.7,
                fwhm_channels=1.3 / 0.3,
            ),
            FwhmCalPeak(
                channel=int(1332.5 / 0.3),
                energy_keV=1332.5,
                fwhm_channels=2.0 / 0.3,
            ),
        ],
    )
    spec = _SpecStub(
        energy_cal=(0.0, 0.3),
        n_channels=8192,
        stored_fwhm_calibration=sf,
    )
    r = evaluate_calibration_gate(spec)
    assert r.passed
    assert r.detector_class_hint == "HPGe"
    g4_warnings = [w for w in r.soft_warnings if w["code"] == "G4"]
    assert not g4_warnings, (
        f"Unexpected G4 warnings for healthy HPGe cal: {g4_warnings}"
    )


# ---------------------------------------------------------------------------
# Integration — extras populated on pipeline read path
# ---------------------------------------------------------------------------

def test_extras_populated_when_gate_runs_via_pipeline_hook():
    """Direct invocation mirrors the hook in staged_pipeline.py.

    Verifies the contract that downstream report.json consumers rely on:
    ``spec.extras['calibration_gate']`` is a JSON-friendly dict with the
    expected top-level keys.
    """
    spec = _SpecStub(energy_cal=(0.0, 0.4), n_channels=8192)
    r = evaluate_calibration_gate(spec)
    spec.extras["calibration_gate"] = r.as_dict()
    payload = spec.extras["calibration_gate"]
    expected_keys = {
        "passed", "hard_failures", "soft_warnings",
        "criteria_evaluated", "detector_class_hint",
        "E_low_keV", "E_high_keV", "reason",
    }
    assert expected_keys.issubset(payload.keys())
    assert isinstance(payload["passed"], bool)
    assert isinstance(payload["hard_failures"], list)
    assert isinstance(payload["soft_warnings"], list)


# ---------------------------------------------------------------------------
# Sanity — constants
# ---------------------------------------------------------------------------

def test_band_constants_sanity():
    """Band thresholds must be ordered and physically reasonable."""
    assert HPGE_BAND_R[0] < HPGE_BAND_R[1]
    assert SCINT_BAND_R[0] < SCINT_BAND_R[1]
    # HPGe band must be entirely below scintillator band (no overlap)
    assert HPGE_BAND_R[1] <= SCINT_BAND_R[0]
    # Channel-range floors physically reasonable
    assert E_LOW_KEV_MIN < 0       # small negative intercept tolerated
    assert E_HIGH_KEV_MIN > 0      # spectrum must cover positive range


# ---------------------------------------------------------------------------
# G3_anchor / G4_anchor — anchor-mode evaluators (BUG-50)
# ---------------------------------------------------------------------------
#
# Anchor-mode activates on the LSRM corpus where the stored FWHM model is
# ``lsrm_fwhm_polynomial_in_E`` (coefficients-only, no measured cal peaks).
# The gate falls back to evaluating the polynomial at every standard anchor
# energy that lies inside the spectrum's calibrated range. Tests below
# build a synthetic LSRM-like ``StoredFwhmCalibration`` so that the
# polynomial argument is z = √E_keV per LSRM Algorithmic Foundations §8.3.
#
# FWHM_keV(E) = Σ_k c_k · z^k,  z = √E_keV
# For a healthy NaI-like fit on the Co-60 / Cs-137 / K-40 region, the
# polynomial passes through ~46 keV @ 662 keV and ~72 keV @ 1332 keV
# (from validate_certs.py:465-467 documentation).

import math as _math


def _lsrm_fwhm_at(coefs, E_keV):
    """Mirror the production model: FWHM_keV(E) = Σ c_k · (√E)^k."""
    z = _math.sqrt(E_keV)
    return sum(c * (z ** k) for k, c in enumerate(coefs))


# Healthy NaI-like FWHM coefficients that reproduce the validate_certs.py
# §465-467 reference points: ~46 keV @ 662 (≈7% R) and ~72 keV @ 1332
# (≈5.4% R). Constructed analytically: solve a·z = FWHM at z=√E for two
# anchors, then add a zero linear term so the curve is monotone in √E.
# z1 = √661.657 ≈ 25.72; z2 = √1332.492 ≈ 36.50.
# Use simple a·z form: FWHM(E) = a·√E with a = 46/25.72 ≈ 1.789
# → @ 1332 keV: 1.789 · 36.50 ≈ 65.3 keV (5.0% R). Both NaI-plausible.
_HEALTHY_NAI_COEFS = (0.0, 1.789)
# Sanity at construction time — anchors land where the docstring claims
assert 40.0 < _lsrm_fwhm_at(_HEALTHY_NAI_COEFS, 661.657) < 60.0
assert 60.0 < _lsrm_fwhm_at(_HEALTHY_NAI_COEFS, 1332.492) < 80.0


def test_g3_anchor_mode_triggers_on_sparse_peaks():
    """LSRM file with empty calibration_peaks → anchor-mode runs.

    Healthy NaI polynomial → monotone increasing FWHM across the anchor
    set → no G3_anchor warning emitted.
    """
    sf = StoredFwhmCalibration(
        coefficients=_HEALTHY_NAI_COEFS,
        model="lsrm_fwhm_polynomial_in_E",
        calibration_peaks=[],   # the canonical LSRM-corpus shape
    )
    # NaI-class energy axis: 0.4 keV/ch over 8192 ch → 0–3276 keV
    spec = _SpecStub(
        energy_cal=(0.0, 0.4),
        n_channels=8192,
        stored_fwhm_calibration=sf,
    )
    r = evaluate_calibration_gate(spec)
    assert r.passed, f"Expected pass, soft={r.soft_warnings}"
    # criteria must show the _anchor suffix because cal_peaks=[] (<4)
    assert "G3_anchor" in r.criteria_evaluated
    assert "G4_anchor" in r.criteria_evaluated
    # No regular G3/G4 entries should appear when anchor-mode ran
    assert "G3" not in r.criteria_evaluated
    assert "G4" not in r.criteria_evaluated
    # No G3_anchor soft warning on a monotone polynomial
    g3a = [w for w in r.soft_warnings if w["code"] == "G3_anchor"]
    assert not g3a, f"Unexpected G3_anchor warning: {g3a}"


def test_g3_anchor_mode_warns_on_non_monotone_fwhm():
    """Polynomial that *decreases* with energy → G3_anchor soft warning.

    Construct coefficients so FWHM(E) = 100 − 0.05·z; over the in-range
    anchor sequence the FWHM strictly decreases as z grows.
    """
    coefs = (100.0, -2.5)
    # Sanity at construction time — the model decreases over the anchor
    # range used by the test (Co-57 122 keV → Tl-208 2614 keV).
    f_lo = _lsrm_fwhm_at(coefs, 122.0607)
    f_hi = _lsrm_fwhm_at(coefs, 2614.511)
    assert f_hi < f_lo, (f_lo, f_hi)

    sf = StoredFwhmCalibration(
        coefficients=coefs,
        model="lsrm_fwhm_polynomial_in_E",
        calibration_peaks=[],
    )
    spec = _SpecStub(
        energy_cal=(0.0, 0.4),
        n_channels=8192,
        stored_fwhm_calibration=sf,
    )
    r = evaluate_calibration_gate(spec)
    # Soft warning (G3 is non-blocking) — verdict is PASS_WITH_SOFT
    assert r.passed
    g3a = [w for w in r.soft_warnings if w["code"] == "G3_anchor"]
    assert g3a, f"Expected G3_anchor warning; soft={r.soft_warnings}"
    # Message must reference monotone violation
    assert "non-monotone" in g3a[0]["message"].lower()


def test_g4_anchor_mode_warns_on_excessive_fwhm():
    """FWHM at a standard anchor far above NaI band (20% × 2 overshoot)
    → G4_anchor soft warning.

    Construct coefficients so FWHM(1460.82 keV) ≈ 800 keV → R ≈ 55%,
    well above the scintillator band even after the 2× overshoot
    margin (scintillator hi=20%, hi_eff=40%).
    """
    # Pick coefs that yield ~800 keV @ K-40 anchor: FWHM = a·√E with
    # a = 800/√1460.82 ≈ 20.93 → R(1460.82) ≈ 55%. Far above hi_eff=40%.
    coefs = (0.0, 20.93)
    f_k40 = _lsrm_fwhm_at(coefs, 1460.822)
    R_k40 = f_k40 / 1460.822
    assert R_k40 > 0.40, (f_k40, R_k40)   # exceeds hi_eff for scintillator

    sf = StoredFwhmCalibration(
        coefficients=coefs,
        model="lsrm_fwhm_polynomial_in_E",
        calibration_peaks=[],
    )
    spec = _SpecStub(
        energy_cal=(0.0, 0.4),
        n_channels=8192,
        stored_fwhm_calibration=sf,
    )
    r = evaluate_calibration_gate(spec)
    assert r.passed   # G4 soft, not blocking
    g4a = [w for w in r.soft_warnings if w["code"] == "G4_anchor"]
    assert g4a, f"Expected G4_anchor warning; soft={r.soft_warnings}"
    # The K-40 anchor (1460.82) must appear in at least one warning
    k40_msg = [w for w in g4a if "K-40" in w["message"]]
    assert k40_msg, f"K-40 anchor missing from G4_anchor list: {g4a}"
    # Class hint must be scintillator (R @ 662 ≈ 21%)
    assert r.detector_class_hint == "scintillator"


def test_anchor_mode_skipped_when_enough_measured_peaks():
    """≥4 usable measured cal peaks → measured-mode runs (no _anchor tags).

    Reproduces a vendor file that carries an explicit calibration_peaks
    list (AtomSpectra-style behaviour, where measured cal peaks are
    always present). The gate stays in the original measured-mode
    G3/G4 branch and ``criteria_evaluated`` lists plain "G3" / "G4".
    """
    # 4 cal peaks across 200–1500 keV with monotone FWHM @ HPGe band
    a1 = 0.3
    peaks = [
        FwhmCalPeak(channel=int(E / a1), energy_keV=E, fwhm_channels=(E * 0.0015) / a1)
        for E in (200.0, 500.0, 1000.0, 1500.0)
    ]
    sf = StoredFwhmCalibration(
        coefficients=(),
        model="",   # AtomSpectra-style — measured peaks dominate
        calibration_peaks=peaks,
    )
    spec = _SpecStub(
        energy_cal=(0.0, a1),
        n_channels=8192,
        stored_fwhm_calibration=sf,
    )
    r = evaluate_calibration_gate(spec)
    assert r.passed
    # No _anchor suffix in this regime
    assert "G3_anchor" not in r.criteria_evaluated
    assert "G4_anchor" not in r.criteria_evaluated
    assert "G3" in r.criteria_evaluated
    assert "G4" in r.criteria_evaluated


def test_anchor_mode_skips_anchors_outside_energy_range():
    """Truncated spectrum (50–1500 keV) → 2614.51 anchor silently skipped.

    The gate should not crash on out-of-range anchors and should not
    emit a warning when the polynomial cannot be evaluated at them.
    """
    sf = StoredFwhmCalibration(
        coefficients=_HEALTHY_NAI_COEFS,
        model="lsrm_fwhm_polynomial_in_E",
        calibration_peaks=[],
    )
    # 1024 ch × ~1.466 keV/ch ≈ 1500 keV at the top, well below
    # Tl-208 anchor 2614.51. Cs-137 (661.7) and Co-60 (1173, 1332)
    # are still in range, plus K-40 (1460.8) and Co-57 (122.06).
    a1 = 1500.0 / 1023.0
    spec = _SpecStub(
        energy_cal=(50.0, a1),
        n_channels=1024,
        stored_fwhm_calibration=sf,
    )
    r = evaluate_calibration_gate(spec)
    assert r.passed
    # Polynomial is monotone → no G3_anchor warning
    g3a = [w for w in r.soft_warnings if w["code"] == "G3_anchor"]
    assert not g3a, f"Unexpected G3_anchor on healthy poly: {g3a}"
    # No G4_anchor warning should mention Tl-208 (out of range)
    g4a = [w for w in r.soft_warnings if w["code"] == "G4_anchor"]
    for w in g4a:
        assert "Tl-208" not in w["message"], (
            f"Out-of-range Tl-208 anchor leaked into G4_anchor: {w}"
        )


def test_standard_fwhm_anchors_sanity():
    """STANDARD_FWHM_ANCHORS_KEV must be sorted, non-empty, well-sourced."""
    assert STANDARD_FWHM_ANCHORS_KEV
    energies = [E for (E, _label, _src) in STANDARD_FWHM_ANCHORS_KEV]
    assert energies == sorted(energies), "Anchor energies must be ascending"
    # Each anchor carries a nuclide label and a citation tag
    for (E, label, src) in STANDARD_FWHM_ANCHORS_KEV:
        assert isinstance(E, float) and E > 0
        assert isinstance(label, str) and label
        assert isinstance(src, str) and src
    # Constants used by the orchestrator must have sane shapes
    assert MEASURED_PEAK_THRESHOLD >= 1
    assert ANCHOR_OVERSHOOT_FACTOR >= 1.0
