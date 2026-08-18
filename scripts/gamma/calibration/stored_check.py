"""
Stored-calibration verification.

Implements the skip rule from SKILL.md §5: if the file's stored
calibration produces residuals < 0.3 · FWHM at every confirmed library
anchor line, reuse it. Otherwise rebuild from scratch.

The verification needs:
  - A spectrum with at least a stored energy calibration
  - A list of found peaks (from `gamma.peaks.search`)
  - A library of expected anchor lines (from `gamma.data.anchors`)
  - An estimate of the local FWHM (from stored FWHM cal or fallback)

This module returns a structured result, not a boolean — the AI needs
to see *why* a check passed or failed (which lines matched, what the
residuals are) to make the autonomous-staging decision.

Token economy: result dict has fixed size O(N_anchors). Counts array
is not in the result.

────────────────────────────────────────────────────────────────────
v1.6 — Adaptive FWHM matching window
────────────────────────────────────────────────────────────────────
Match window is now a true 1·FWHM(E_anchor) at every anchor, sourced
from the vendor's stored FWHM model via
`make_fwhm_at_channel_provider`. The old `fwhm_fallback_channels=10`
became too tight for high-energy NaI peaks (Tl-208 2614 keV has
FWHM ≈ 90 keV on a 50×50 crystal, three orders of magnitude wider
than the previous fixed 10-channel ≈ 6.5 keV window). With the fix,
Tl-208 now matches consistently on natural-background spectra.

A bug was corrected at the same time: the previous code interpreted
`SimpleSqrtFwhm` as FWHM = c₀ + c₁·√N. Direct check against the file's
embedded `<CalibrationPeaks>` shows the actual convention is FWHM² =
c₀ + c₁·N (the “Sqrt” names the outer square root, not √N). The fix
lives in `gamma.calibration.fwhm_provider`; the previous code in this
module is removed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from gamma.spectrum import Spectrum
from gamma.peaks.search import FoundPeak
from gamma.calibration.fwhm_provider import make_fwhm_at_channel_provider


@dataclass
class StoredCheckResult:
    """
    Outcome of stored-calibration verification.

    Attributes:
        passed: True if stored calibration meets the skip criterion.
        n_anchors_tested: How many library energies we attempted to match.
        n_anchors_matched: How many we successfully matched to a peak.
        max_residual_keV: Largest |E_observed − E_expected| over matched
                          anchors.
        max_residual_over_fwhm: Same residual divided by local FWHM (the
                                actual skip-criterion value; must be
                                < 0.3 for `passed=True`).
        matches: List of dicts (one per matched anchor) with details.
        unmatched: List of energies we expected but couldn't find.
        threshold: The fwhm fraction threshold used (default 0.3).
        reason: Human-readable summary.
        fwhm_source: "stored_model" | "calibration_peaks" | "fallback"
                     — for diagnostics, which source produced the
                     matching window.
    """
    passed: bool
    n_anchors_tested: int
    n_anchors_matched: int
    max_residual_keV: float
    max_residual_over_fwhm: float
    matches: list = field(default_factory=list)
    unmatched: list = field(default_factory=list)
    threshold: float = 0.3
    reason: str = ""
    fwhm_source: str = ""


def check_stored_calibration(
    spec: Spectrum,
    found_peaks: list,
    *,
    threshold_fwhm: float = 0.3,
    fwhm_fallback_channels: float = 10.0,
    anchor_priorities: Optional[list] = None,
    match_window_fwhm: float = 1.0,
) -> StoredCheckResult:
    """
    Verify the stored calibration against found peaks using anchor lines.

    Args:
        spec: parsed Spectrum (must have spec.energy_cal set)
        found_peaks: output from mariscotti_search()
        threshold_fwhm: skip-criterion threshold (default 0.3·FWHM)
        fwhm_fallback_channels: used only when neither stored FWHM
            model nor calibration peaks are available (default 10).
        anchor_priorities: list of priorities to try (default [1, 2]).
                          Lower priority = stronger anchor.
        match_window_fwhm: matching window size in multiples of local
                          FWHM. Default 1.0 — a peak must be within 1
                          FWHM of an anchor energy to be considered a
                          match. The skip rule then evaluates whether
                          the residual is within `threshold_fwhm`.

    Returns:
        StoredCheckResult with passed=True if all matched anchors fall
        within threshold_fwhm of their library energy.
    """
    if spec.energy_cal is None:
        return StoredCheckResult(
            passed=False,
            n_anchors_tested=0,
            n_anchors_matched=0,
            max_residual_keV=float("inf"),
            max_residual_over_fwhm=float("inf"),
            reason="No stored calibration present",
        )

    if not found_peaks:
        return StoredCheckResult(
            passed=False,
            n_anchors_tested=0,
            n_anchors_matched=0,
            max_residual_keV=float("inf"),
            max_residual_over_fwhm=float("inf"),
            reason="No peaks found — cannot verify stored calibration",
        )

    # Late import to avoid circular load issues
    from gamma.data.anchors import patterns_by_priority

    # Collect candidate anchor energies, deduplicated
    if anchor_priorities is None:
        anchor_priorities = [1, 2]
    expected_energies = set()
    for prio in anchor_priorities:
        for pat in patterns_by_priority(max_priority=prio):
            for E in pat.get("lines", []):
                expected_energies.add(float(E))
    expected_energies = sorted(expected_energies)

    # Build the FWHM-at-channel provider (single source of truth) and
    # remember which source it ended up using for diagnostics
    fwhm_at_ch = make_fwhm_at_channel_provider(
        spec, fallback_channels=fwhm_fallback_channels,
    )
    fwhm_source = _diagnose_fwhm_source(spec)

    # Local FWHM estimator at energy E:
    #   1. translate E → channel via stored energy cal
    #   2. call fwhm_at_ch(channel) → fwhm_channels
    #   3. convert fwhm_channels → keV via local dE/dN at that channel
    def fwhm_at_E_keV(E: float) -> float:
        ch = spec.energy_to_channel(E)
        if ch is None:
            return _channels_to_keV(spec, fwhm_fallback_channels, E)
        ch_int = int(round(float(ch)))
        fwhm_ch = fwhm_at_ch(ch_int)
        return _channels_to_keV(spec, fwhm_ch, E)

    # Peak energies in keV under stored cal
    peak_energies = np.array(
        [spec.channel_to_energy(p.channel) for p in found_peaks],
        dtype=np.float64,
    )

    matches = []
    unmatched = []
    max_residual_keV = 0.0
    max_residual_over_fwhm = 0.0

    for E_expected in expected_energies:
        fwhm = fwhm_at_E_keV(E_expected)
        if fwhm <= 0:
            continue
        deltas = np.abs(peak_energies - E_expected)
        idx = int(np.argmin(deltas))
        delta = float(peak_energies[idx] - E_expected)

        if abs(delta) <= match_window_fwhm * fwhm:
            residual_keV = abs(delta)
            residual_over_fwhm = residual_keV / fwhm
            matches.append({
                "anchor_keV": E_expected,
                "matched_peak_channel": found_peaks[idx].channel,
                "matched_peak_energy_keV": float(peak_energies[idx]),
                "residual_keV": residual_keV,
                "fwhm_keV": fwhm,
                "residual_over_fwhm": residual_over_fwhm,
            })
            max_residual_keV = max(max_residual_keV, residual_keV)
            max_residual_over_fwhm = max(max_residual_over_fwhm, residual_over_fwhm)
        else:
            unmatched.append(E_expected)

    if not matches:
        return StoredCheckResult(
            passed=False,
            n_anchors_tested=len(expected_energies),
            n_anchors_matched=0,
            max_residual_keV=float("inf"),
            max_residual_over_fwhm=float("inf"),
            unmatched=unmatched,
            threshold=threshold_fwhm,
            reason=(f"None of {len(expected_energies)} anchor energies "
                    f"matched any found peak"),
            fwhm_source=fwhm_source,
        )

    passed = max_residual_over_fwhm < threshold_fwhm
    if passed:
        reason = (f"Stored calibration OK: max residual "
                  f"{max_residual_keV:.2f} keV = "
                  f"{max_residual_over_fwhm:.2f}·FWHM "
                  f"(< {threshold_fwhm}·FWHM) over "
                  f"{len(matches)} anchor lines")
    else:
        reason = (f"Stored calibration fails skip rule: max residual "
                  f"{max_residual_keV:.2f} keV = "
                  f"{max_residual_over_fwhm:.2f}·FWHM "
                  f"(needs < {threshold_fwhm}·FWHM)")

    return StoredCheckResult(
        passed=passed,
        n_anchors_tested=len(expected_energies),
        n_anchors_matched=len(matches),
        max_residual_keV=max_residual_keV,
        max_residual_over_fwhm=max_residual_over_fwhm,
        matches=matches,
        unmatched=unmatched,
        threshold=threshold_fwhm,
        reason=reason,
        fwhm_source=fwhm_source,
    )


# ============================================================================
# Diagnostics
# ============================================================================

def _diagnose_fwhm_source(spec: Spectrum) -> str:
    """Return a short tag identifying which FWHM source will be used."""
    sf = getattr(spec, "stored_fwhm_calibration", None)
    if sf is None:
        return "fallback"
    # Mirror the resolution order of make_fwhm_at_channel_provider.
    # If SimpleSqrtFwhm model gives sensible values at every cal peak
    # (within 30%), the model wins. Otherwise we drop to cal peaks.
    if sf.model == "SimpleSqrtFwhm" and len(sf.coefficients) >= 2:
        import math as _m
        c0, c1 = float(sf.coefficients[0]), float(sf.coefficients[1])
        cal = list(sf.calibration_peaks or [])
        if cal:
            ok = True
            for cp in cal:
                if cp.fwhm_channels and cp.fwhm_channels > 0:
                    arg = c0 + c1 * cp.channel
                    pred = _m.sqrt(arg) if arg > 0 else 0.0
                    if pred <= 0 or abs(pred - cp.fwhm_channels) / cp.fwhm_channels > 0.30:
                        ok = False
                        break
            if ok:
                return "stored_model"
            else:
                return "calibration_peaks"
        return "stored_model"
    if sf.calibration_peaks:
        return "calibration_peaks"
    return "fallback"


# ============================================================================
# Helpers — channel ↔ keV conversion at a given energy
# ============================================================================

def _channels_to_keV(spec: Spectrum, fwhm_channels: float, E_keV: float) -> float:
    """
    Convert an FWHM in channels to keV using the local slope of the
    stored energy calibration evaluated at E_keV.
    """
    if spec.energy_cal is None or len(spec.energy_cal) < 2:
        return float(fwhm_channels)
    ch = spec.energy_to_channel(E_keV)
    if ch is None:
        return float(fwhm_channels)
    dE_dN = sum(
        i * a * (ch ** (i - 1))
        for i, a in enumerate(spec.energy_cal) if i > 0
    )
    if dE_dN <= 0:
        return float(fwhm_channels)
    return float(fwhm_channels * dE_dN)
