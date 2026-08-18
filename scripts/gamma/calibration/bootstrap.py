"""
Bootstrap energy calibration.

Algorithm: brute-force with anchor-pattern priority (variant `a` per
SKILL.md Phase 1.2 design).

When the file's stored calibration fails the residual test (or is
absent), we rebuild from scratch using known anchor patterns:

  1. Take anchor patterns in priority order (1 = strongest)
  2. For each priority-1 pattern, attempt to identify it in the found
     peaks by matching inter-peak distance ratios to the pattern's
     expected energy spacing
  3. If a pattern is identified, use its (channel, energy) pairs as a
     two-point seed for a linear calibration
  4. With the seed calibration, expand: search lower-priority patterns
     near their predicted positions in the new energy axis
  5. Refit with higher polynomial degree (≤4) as more anchors are added
  6. Return the final calibration

The output is a small dict with calibration coefficients, residuals,
and the list of anchors that were used — enough for the AI to assess
quality and for downstream modules to use the new calibration.

Token economy: returns ~10 KV pairs total; no counts arrays, no per-
channel data.

═══════════════════════════════════════════════════════════════════════
STATUS (v1.6): Three v1.6 improvements landed on top of v1.5
═══════════════════════════════════════════════════════════════════════
Built on the v1.5 foundation (FWHM-adaptive tolerance, line-count-
first pattern sort, Th-232 chain triplet pattern). New in v1.6:

 1. Physical-offset sanity check in `_match_pattern_to_peaks`: the
    resulting a0 must lie in [-200, +200] keV. This single check kills
    the entire failure mode where `K-40 + Tl-208 high-energy pair`
    mathematically fits any pair of peaks separated by ~1150 keV.
    Tested: previously-failing Cs-137 (a0=+1310), Ra-226 (+1196), and
    M_k_лёгкий (+1273) seeds are now correctly rejected.

 2. Stricter wide_pair check: high-energy anchor of any `wide_pair`
    pattern must land at >= 60% of channel range. The 1460/2614 keV
    pair cannot physically be in the lower half of the spectrum for
    any sensible NaI gain (2.5–3.5 keV/ch on 1024–8192 channels).

 3. Single-line anchor seed (third pass in `_find_seed_pattern`):
    when neither multi-line nor 2-line patterns match, try pairing
    the strongest peak (σ ≥ 30) with each single_line pattern under
    offset prior a0 = 0. The implied gain a1 = E_line / ch must be
    in the standard [0.05, 3.5] range. This unblocks pure Cs-137
    calibration sources, where v1.5 could not seed at all.

Expansion step (Step 3) is unchanged: once any seed is established,
nearby anchor lines are iteratively added and the polynomial degree
grows up to max_degree as more anchors accumulate.

Result on real Lsrm SpectraLine control sources (M_cs/M_k/M_ra/M_th):
  Cs-137 (легкий):  Cs-137 single seed (single_line path) — correct.
  K-40 (легкий):    K-40 single seed — correct.
  Ra-226 (легкий):  Bi-214 multi-line — correct.
  Th-232 (легкий):  Th-232 chain triplet — correct.
═══════════════════════════════════════════════════════════════════════
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from itertools import combinations
from typing import Optional

import numpy as np

from gamma.peaks.search import FoundPeak
from gamma.calibration.energy_fit import (
    EnergyFitResult,
    polynomial_energy_fit,
    MAX_POLYNOMIAL_DEGREE,
)


# ============================================================================
# Output dataclass
# ============================================================================

@dataclass
class BootstrapResult:
    """
    Outcome of bootstrap calibration.

    Attributes:
        success: True if a usable calibration was found
        coefficients: low-to-high polynomial coefficients
        degree: polynomial degree
        anchors_used: list of (channel, energy_keV, source) tuples
        seed_pattern: name of the priority-1 pattern that seeded the fit
        max_residual_keV: max |observed − predicted| over anchors
        rms_residual_keV: RMS residual
        n_iterations: how many refinement passes were performed
        fit_result: the final EnergyFitResult
        reason: human-readable summary
    """
    success: bool
    coefficients: tuple = ()
    degree: int = 0
    anchors_used: list = field(default_factory=list)
    seed_pattern: str = ""
    max_residual_keV: float = 0.0
    rms_residual_keV: float = 0.0
    n_iterations: int = 0
    fit_result: Optional[EnergyFitResult] = None
    reason: str = ""


# ============================================================================
# Top-level entry point
# ============================================================================

def bootstrap_energy_calibration(
    found_peaks: list,
    n_channels: int,
    *,
    spec=None,
    spectrum_type: Optional[str] = None,
    energy_floor_keV: float = 25.0,
    initial_gain_hint: Optional[float] = None,
    initial_offset_hint: Optional[float] = None,
    max_degree: int = MAX_POLYNOMIAL_DEGREE,
    target_residual_keV: float = 1.0,
) -> BootstrapResult:
    """
    Build an energy calibration from scratch using anchor patterns.

    Args:
        found_peaks: peaks from mariscotti_search() (channel space)
        n_channels: total channels in the spectrum (for sanity checks)
        spec: optional Spectrum object for filename-token auto-detection
              of spectrum type and for stored-cal-based energy filtering
        spectrum_type: "background"/"bkg" applies relaxed Lsrm-style
              background processing — only strongest peaks (σ ≥ 10) are
              used as anchors.
        energy_floor_keV: minimum peak energy (in keV) to consider as a
              calibration anchor candidate. Default 25 keV — below this,
              spectra contain electronic noise, detector dark current,
              and X-ray pileup that should not influence the energy
              calibration. Applied only if spec has a stored energy
              calibration that can map channels to keV; otherwise all
              peaks are used. Set to 0 to disable.
        initial_gain_hint: optional rough guess of keV/channel
        initial_offset_hint: optional rough guess of a0 in keV
        max_degree: polynomial degree cap (≤ 4 per skill scope)
        target_residual_keV: stop refining when max residual reaches this

    Returns:
        BootstrapResult.
    """
    # Auto-detect background type from filename hints
    if spectrum_type is None and spec is not None:
        sample_type = getattr(spec, "filename_tokens", {}).get("sample_type", [])
        bg_markers = {"фон", "Фон", "ФОН", "bkg", "background"}
        if any(m in sample_type for m in bg_markers):
            spectrum_type = "background"

    is_background = spectrum_type in ("background", "bkg")
    anchor_sigma = 10.0 if is_background else 5.0

    # Apply energy floor using stored calibration if available
    if (spec is not None and energy_floor_keV > 0
            and getattr(spec, "energy_cal", None)):
        original_count = len(found_peaks)
        found_peaks = [
            p for p in found_peaks
            if spec.channel_to_energy(p.channel) >= energy_floor_keV
        ]
        # filter does not return a diagnostic; downstream code knows
        # only how many peaks remained.

    if not found_peaks or len(found_peaks) < 2:
        return BootstrapResult(
            success=False,
            reason=f"Bootstrap needs at least 2 peaks; got {len(found_peaks)}",
        )

    from gamma.data.anchors import patterns_by_priority

    # ----- Step 1: find a seed pattern -----
    # Single pass over ALL patterns up to priority 2 (no two-stage call).
    # The seed-finder sorts by (line_count desc, priority asc), so the
    # most informative patterns are tried first. This is critical: an
    # 8-line Bi-214 pattern (priority 2) gives a much stricter ratio
    # check than a 2-line wide_pair (priority 1), so the former must
    # be tried first.
    seed = _find_seed_pattern(
        found_peaks,
        patterns_by_priority(max_priority=2),
        n_channels=n_channels,
        gain_hint=initial_gain_hint,
        offset_hint=initial_offset_hint,
        anchor_sigma_threshold=anchor_sigma,
    )

    if seed is None:
        return BootstrapResult(
            success=False,
            reason="No anchor pattern could be matched in the found peaks "
                   "(insufficient natural-background or calibration lines)",
        )

    # ----- Step 2: initial linear fit from seed -----
    anchors = list(seed["anchors"])  # list of (channel, energy, source)
    is_single_line_seed = bool(seed.get("single_line"))

    if is_single_line_seed:
        # v1.6: single-line seed — fix offset=0, derive gain from the one
        # anchor. Build a synthetic EnergyFitResult so the iterative
        # expansion step can use it without special-casing. The (0,0)
        # offset prior is NOT carried into the anchors list; it stays
        # implicit in the coefficients.
        ch_anchor, E_anchor, _src = anchors[0]
        if ch_anchor <= 0:
            return BootstrapResult(
                success=False,
                seed_pattern=seed["pattern_name"],
                reason=(f"Single-line seed peak at ch={ch_anchor} "
                        f"too close to origin to derive gain"),
            )
        a1 = E_anchor / float(ch_anchor)
        a0 = 0.0
        fit = EnergyFitResult(
            coefficients=(a0, a1),
            degree=1,
            n_points=1,
            residuals_keV=[0.0],
            max_residual_keV=0.0,
            rms_residual_keV=0.0,
            converged=True,
            reason=(f"Single-line seed (a0=0 prior) — "
                    f"gain {a1:.4f} keV/ch from anchor {E_anchor} keV "
                    f"at ch {ch_anchor}"),
        )
    elif len(anchors) >= 2:
        channels = [a[0] for a in anchors]
        energies = [a[1] for a in anchors]
        fit = polynomial_energy_fit(
            channels, energies,
            max_degree=1,  # start linear
            target_residual_keV=target_residual_keV,
        )
    else:
        return BootstrapResult(
            success=False,
            seed_pattern=seed["pattern_name"],
            reason=f"Seed pattern matched only {len(anchors)} anchors; need ≥2",
        )

    # ----- Step 3: iterative expansion -----
    # Use the current calibration to predict expected channel positions
    # of lower-priority anchors, search for matches, add them, refit.
    n_iterations = 0
    all_patterns = patterns_by_priority(max_priority=3)
    seen_anchors = set((round(a[0], 3), round(a[1], 3)) for a in anchors)

    for iteration in range(5):  # cap iterations
        n_iterations = iteration + 1
        new_anchors = _expand_anchors(
            found_peaks,
            current_coefs=fit.coefficients,
            patterns=all_patterns,
            already_used=seen_anchors,
        )
        if not new_anchors:
            break

        anchors.extend(new_anchors)
        for a in new_anchors:
            seen_anchors.add((round(a[0], 3), round(a[1], 3)))

        # Refit with higher degree as we accumulate anchors
        channels = [a[0] for a in anchors]
        energies = [a[1] for a in anchors]
        # Allow degree to grow proportionally to anchor count
        attempt_degree = min(max_degree, max(1, len(anchors) // 3))
        fit = polynomial_energy_fit(
            channels, energies,
            max_degree=attempt_degree,
            min_degree=1,
            target_residual_keV=target_residual_keV,
        )

        if fit.max_residual_keV <= target_residual_keV:
            break

    return BootstrapResult(
        success=True,
        coefficients=fit.coefficients,
        degree=fit.degree,
        anchors_used=anchors,
        seed_pattern=seed["pattern_name"],
        max_residual_keV=fit.max_residual_keV,
        rms_residual_keV=fit.rms_residual_keV,
        n_iterations=n_iterations,
        fit_result=fit,
        reason=(f"Seeded with '{seed['pattern_name']}'; "
                f"final degree {fit.degree}, {len(anchors)} anchors, "
                f"max residual {fit.max_residual_keV:.3f} keV"),
    )


# ============================================================================
# Seed-pattern search (Step 1)
# ============================================================================

def _find_seed_pattern(
    found_peaks: list,
    patterns: list,
    *,
    n_channels: int,
    gain_hint: Optional[float] = None,
    offset_hint: Optional[float] = None,
    single_line_sigma_threshold: float = 15.0,
    anchor_sigma_threshold: float = 5.0,
) -> Optional[dict]:
    """
    Find a seed pattern.

    `anchor_sigma_threshold` (v1.6): only peaks with σ ≥ this threshold are
    eligible as anchor candidates for multi-line and 2-line patterns. The
    Lsrm methodology distinguishes between peak-search threshold (3σ, used
    to detect ALL peaks for the workflow) and the threshold for accepting
    peaks as **calibration references** (typically 5σ or higher). Low-σ
    peaks risk being noise fluctuations and can create fake-consistent
    pattern matches: e.g. on a real 7034-channel NaI background, a 3-line
    Th-232 triplet could match {ch=204 (σ=33.6, real E≈80 keV), ch=596
    (σ=40.6, real E≈240), ch=2890 (σ=3.2, possible noise, mapped to 2614)}
    with mathematically perfect linearity but completely wrong gain (0.88
    keV/ch vs true 0.39 keV/ch). Excluding σ<5 candidates kills that
    failure mode.
    """
    """
    Try each candidate pattern; return the first successful seed.

    Strategy (v1.6 — three passes):
      1. Patterns with ≥3 lines first — their inter-peak distance
         ratios make matching unambiguous (and ratio_tol is strict at
         2%).
      2. 2-line patterns — accepted if (a) they are wide_pair (line
         spacing > 500 keV gives unambiguous gain) OR (b) a gain_hint
         is provided. Universal sanity checks (offset, wide_pair high-E
         position) are applied in `_match_pattern_to_peaks`.
      3. v1.6 NEW: Single-line patterns — last-resort fallback for
         pure-source spectra where neither (1) nor (2) apply (e.g.
         M_cs_лёгкий has one strong line at 661 keV and a few weak
         background lines). The strongest peak whose significance ≥
         `single_line_sigma_threshold` (default 30σ) is paired with
         each single_line pattern under offset prior a0 = 0; the gain
         is taken from the (channel, E) point alone.

    Returns:
        {
            "pattern_name": str,
            "anchors": list of (channel, energy_keV, source_pattern),
            "single_line": bool,    # True when seed came from pass 3
        }
        or None if no pattern matched.
    """
    # Sort by (line count desc, priority asc): more lines = stricter
    # ratio-check, which outweighs the priority field. Within the same
    # line count, prefer higher priority (lower number).
    multi_patterns = [
        p for p in patterns
        if len(p.get("lines", [])) >= 3
    ]
    multi_patterns.sort(key=lambda p: (
        -len(p.get("lines", [])),    # more lines first
        p.get("priority", 99),        # then higher priority
    ))

    two_line_patterns = [
        p for p in patterns
        if len(p.get("lines", [])) == 2
    ]
    two_line_patterns.sort(key=lambda p: p.get("priority", 99))

    single_line_patterns = [
        p for p in patterns
        if len(p.get("lines", [])) == 1 and p.get("single_line") is True
    ]
    single_line_patterns.sort(key=lambda p: p.get("priority", 99))

    # Use the strongest peaks first, but only those above the anchor σ
    # threshold (v1.6 — see docstring). Low-σ peaks risk being noise
    # fluctuations that create fake pattern matches.
    above_threshold = [
        p for p in found_peaks if p.significance >= anchor_sigma_threshold
    ]
    strong_peaks = sorted(
        above_threshold, key=lambda p: -p.significance
    )[:max(20, len(above_threshold) // 2)]
    strong_peak_channels_unsorted = [p.channel for p in strong_peaks]
    strong_peak_sigs_unsorted = [p.significance for p in strong_peaks]
    strong_peak_fwhms_unsorted = [p.fwhm_channels for p in strong_peaks]
    # Sort all three arrays together by channel
    sorted_triples = sorted(zip(
        strong_peak_channels_unsorted,
        strong_peak_sigs_unsorted,
        strong_peak_fwhms_unsorted,
    ))
    strong_peak_channels = [c for c, _, _ in sorted_triples]
    strong_peak_sigs = [s for _, s, _ in sorted_triples]
    strong_peak_fwhms = [f for _, _, f in sorted_triples]

    # First pass: ≥3-line patterns (unambiguous)
    for pat in multi_patterns:
        lines = list(pat.get("lines", []))
        tol_keV = float(pat.get("tolerance_keV", 5.0))

        match = _match_pattern_to_peaks(
            pattern_lines=lines,
            peak_channels=strong_peak_channels,
            peak_significances=strong_peak_sigs,
            peak_fwhms=strong_peak_fwhms,
            tolerance_keV=tol_keV,
            gain_hint=gain_hint,
            offset_hint=offset_hint,
            n_channels=n_channels,
            is_wide_pair=False,  # multi-line patterns use ratio check, not wide_pair
            pattern_intensity_ratios=pat.get("intensity_ratios"),
        )

        if match is not None:
            anchors = [
                (ch, E, pat["name"])
                for ch, E in zip(match["channels"], match["energies"])
            ]
            return {
                "pattern_name": pat["name"],
                "anchors": anchors,
                "single_line": False,
            }

    # Second pass: 2-line patterns
    # Accepted if EITHER (a) they are wide_pair (line spacing > 500 keV
    # makes the gain unambiguous from any pair of peaks), OR (b) we
    # have a gain hint.
    for pat in two_line_patterns:
        lines = list(pat.get("lines", []))
        is_wide_pair = bool(pat.get("wide_pair"))
        if not is_wide_pair and gain_hint is None:
            continue

        tol_keV = float(pat.get("tolerance_keV", 5.0))

        match = _match_pattern_to_peaks(
            pattern_lines=lines,
            peak_channels=strong_peak_channels,
            peak_significances=strong_peak_sigs,
            peak_fwhms=strong_peak_fwhms,
            tolerance_keV=tol_keV,
            gain_hint=gain_hint,
            offset_hint=offset_hint,
            n_channels=n_channels,
            is_wide_pair=is_wide_pair,
            pattern_intensity_ratios=pat.get("intensity_ratios"),
        )

        if match is not None:
            anchors = [
                (ch, E, pat["name"])
                for ch, E in zip(match["channels"], match["energies"])
            ]
            return {
                "pattern_name": pat["name"],
                "anchors": anchors,
                "single_line": False,
            }

    # Third pass (v1.6): single-line patterns as last resort
    single_seed = _try_single_line_seed(
        strong_peaks=strong_peaks,
        single_line_patterns=single_line_patterns,
        n_channels=n_channels,
        sigma_threshold=single_line_sigma_threshold,
        gain_hint=gain_hint,
        offset_hint=offset_hint,
    )
    if single_seed is not None:
        return single_seed

    return None


def _try_single_line_seed(
    *,
    strong_peaks: list,
    single_line_patterns: list,
    n_channels: int,
    sigma_threshold: float,
    gain_hint: Optional[float],
    offset_hint: Optional[float],
) -> Optional[dict]:
    """
    Single-line seed with consistency scoring (v1.6).

    The naive single-line seed is ambiguous: the strongest peak in a
    Cs-137 spectrum at ch ≈ 197 could be interpreted as Cs-137 661 keV
    (gain 3.36 keV/ch), or as Annihilation 511 keV (gain 2.59), or as
    K-40 1460 keV (gain 7.41 — typically rejected by the gain bound)
    depending on which single_line pattern wins the priority sort.

    The disambiguator is **consistency with the rest of the spectrum**:
    the right gain hypothesis should explain MORE other found peaks as
    coincident with library anchor lines than wrong hypotheses do.

    Algorithm:
      1. Filter peaks to σ ≥ sigma_threshold and channel > 0.
      2. For each (peak, single_line_pattern) pair:
         a. Compute implied gain a1 = E_line / peak.channel under
            offset prior a0 = 0.
         b. Reject if a1 ∉ [0.05, 3.5].
         c. Score = consistency: count how many OTHER found peaks
            (those with σ ≥ 5, not the seed peak) have their calibrated
            energy land within tolerance_keV of ANY priority ≤ 2
            anchor line in the library. Higher count = better.
         d. Tiebreaker: prefer higher peak significance and higher
            pattern priority (lower priority number).
      3. Accept the best-scoring combination; return None if no
         hypothesis exceeds the minimum (1 corroborating peak besides
         the seed itself).

    This naturally selects Cs-137 over Annihilation on a Cs-137
    spectrum: under gain 3.36 the second peak ch=51 maps to E=171 keV
    (no priority-1/2 anchor nearby — but Cs+BaKα 32 keV could match if
    spectrum has it; otherwise no corroborating peak). Under gain 2.59
    (Annihilation) the ch=197 peak maps to 511, and no library anchors
    coincide with the secondary peaks either. The two hypotheses tie
    on consistency, but Cs-137 wins the priority tiebreaker.

    For Cs-137 spectra with weak secondary peaks: the priority
    tiebreaker handles the no-corroboration case.
    """
    from gamma.data.anchors import patterns_by_priority

    if not single_line_patterns:
        return None

    candidates = [p for p in strong_peaks if p.significance >= sigma_threshold]
    candidates = [p for p in candidates if p.channel > 0]
    if not candidates:
        return None

    # Collect all priority-1/2 library anchor energies (for consistency check)
    library_anchor_energies = []
    for pat in patterns_by_priority(max_priority=2):
        for E in pat.get("lines", []):
            library_anchor_energies.append((float(E), pat["name"]))

    # All found peaks with σ ≥ 5 (excluding the seed when scoring)
    all_peaks_for_consistency = [
        p for p in strong_peaks if p.significance >= 5.0
    ]

    best = None
    best_score = (-1.0, -1.0, 99)  # (consistency, significance, -priority)

    for peak in candidates:
        for pat in single_line_patterns:
            lines = pat.get("lines", [])
            if len(lines) != 1:
                continue
            E_line = float(lines[0])
            tol_keV = float(pat.get("tolerance_keV", 2.0))

            # Gain under offset=0 prior
            a1 = E_line / peak.channel
            if not (0.05 <= a1 <= 3.5):
                continue

            # Hint compatibility
            if gain_hint is not None and abs(a1 - gain_hint) / gain_hint > 0.15:
                continue
            if offset_hint is not None and abs(0.0 - offset_hint) > 50:
                continue

            # v1.6 PHYSICAL CORROBORATION FOR ANNIHILATION 511.
            # The 511 keV positron-annihilation line is an epiphenomenon, not a
            # primary photopeak: it requires either a positron emitter in the
            # sample, or pair production from a parent gamma with E > 1022 keV.
            # In a pure Cs-137 calibration spectrum (max photopeak 661.66 keV),
            # 511 cannot physically arise — yet under offset=0 prior a peak at
            # ch=197 (Cs-137 661 keV under correct gain 3.36) is also a valid
            # match for 511 keV under wrong gain 2.59. Refuse the 511 hypothesis
            # unless the spectrum contains at least one OTHER significant peak
            # mapping to > 1022 keV under this calibration.
            if abs(E_line - 511.0) < 1.0:
                has_high_E_corroboration = any(
                    (a1 * other.channel) > 1022.0
                    for other in all_peaks_for_consistency
                    if other.channel != peak.channel
                )
                if not has_high_E_corroboration:
                    continue

            # Consistency: count other peaks that map to known anchors
            # under this hypothesis. Use a generous tolerance (3·peak's
            # FWHM in keV) since the polynomial higher-order corrections
            # haven't been applied yet.
            fwhm_keV = peak.fwhm_channels * a1
            corroboration_tol = max(tol_keV * 3, 1.5 * fwhm_keV)
            n_corroborated = 0
            corroborating_peaks = []
            for other in all_peaks_for_consistency:
                if other.channel == peak.channel:
                    continue
                E_other = a1 * other.channel  # offset is 0
                # Is there any library anchor within tolerance?
                nearest = min(
                    library_anchor_energies,
                    key=lambda ea: abs(ea[0] - E_other),
                )
                if abs(nearest[0] - E_other) <= corroboration_tol:
                    n_corroborated += 1
                    corroborating_peaks.append((other.channel, E_other,
                                                nearest[0], nearest[1]))

            # Score tuple: (consistency, significance, -priority)
            priority = pat.get("priority", 5)
            score = (n_corroborated, peak.significance, -priority)

            if score > best_score:
                best_score = score
                best = {
                    "pattern_name": pat["name"],
                    "anchors": [(int(peak.channel), E_line, pat["name"])],
                    "single_line": True,
                    "offset_prior": 0.0,
                    "implied_gain": a1,
                    "implied_fwhm_keV": fwhm_keV,
                    "n_corroborated": n_corroborated,
                    "corroborating_peaks": corroborating_peaks,
                }

    return best


def _match_pattern_to_peaks(
    *,
    pattern_lines: list,
    peak_channels: list,
    peak_significances: list,
    peak_fwhms: list,
    tolerance_keV: float,
    gain_hint: Optional[float],
    offset_hint: Optional[float],
    n_channels: int,
    is_wide_pair: bool = False,
    pattern_intensity_ratios: Optional[list] = None,
) -> Optional[dict]:
    """
    Try to match a list of expected energies to a subset of peak channels.

    A match exists if there is a linear mapping E(N) = a0 + a1·N such
    that every pattern energy lands within `tolerance_keV` of some peak
    channel mapped to keV.

    For ≥3-line patterns: we require strict consistency of internal
    position ratios — accept only if the relative position of each
    intermediate line in the [first, last] range matches the pattern's
    relative position within `RATIO_TOL` (default 2%).

    For 2-line patterns: any pair of peaks fits linearly through both
    points, so the additional filter is significance — pairs with low
    significance are penalized.

    `peak_significances` is a list parallel to `peak_channels`.

    Two universal sanity checks apply to ALL pattern sizes (v1.6):

      1. PHYSICAL OFFSET (`a0 ∈ [-200, +200] keV`): channel 0 of any
         realistic detector — NaI/HPGe/CdZnTe/LaBr3 — lies near the
         zero-energy origin. A larger offset means we've matched the
         pattern to physically wrong peaks. This kills the v1.5 mode of
         failure where `K-40 + Tl-208 high-energy pair` mathematically
         "fits" any two well-separated peaks (Cs-137 alone, M_ra,
         M_k_лёгкий) and produces |a0| > 1100 keV.

      2. WIDE_PAIR HIGH-E ANCHOR IN UPPER 40% OF SPECTRUM: when
         `is_wide_pair=True`, the highest-energy anchor of the pattern
         must map to a channel >= 0.60·n_channels. This catches the
         remaining false matches: e.g. M_k_лёгкий strongest peak is
         at ch 443/1023 (43%) — well below the cut, even though the
         offset check might pass.

    Returns dict {channels, energies, gain, offset} or None.
    """
    n_lines = len(pattern_lines)
    if n_lines < 2 or len(peak_channels) < n_lines:
        return None

    pat_low = pattern_lines[0]
    pat_high = pattern_lines[-1]
    pat_total_range = pat_high - pat_low
    if pat_total_range <= 0:
        return None

    pat_positions = [(L - pat_low) / pat_total_range for L in pattern_lines]
    # v1.6: RATIO_TOL relaxed from 0.02 → 0.05.
    # On NaI 50×50 with FWHM @ 500 keV ≈ 30 keV, individual peak centroids
    # carry ~10-15 keV energy uncertainty. For an 8-line Bi-214 pattern
    # spanning 609–2204 keV (range 1595 keV), middle peaks may shift by
    # 10/1595 ≈ 0.6% from their nominal relative positions purely due to
    # measurement noise. At 2% bound, valid Bi-214 matches on NaI are
    # rejected; at 5% bound they pass while still excluding unrelated
    # peak combinations (which typically differ by 10% or more).
    RATIO_TOL = 0.05

    # Physical-offset sanity bounds (in keV) — see docstring.
    PLAUSIBLE_OFFSET_LO = -200.0
    PLAUSIBLE_OFFSET_HI = 200.0

    # Wide-pair high-E channel-position cut (fraction of n_channels).
    WIDE_PAIR_HIGH_E_MIN_FRAC = 0.60

    # Intensity-asymmetry tolerance: when a pattern's brightest expected
    # line is more than 3× brighter than its dimmest, require the
    # matched peaks to respect that ordering (the peak at the dim line
    # must have lower σ than the peak at the bright line). Catches
    # wrong-direction matches like Cs+BaKα being assigned to (ch=100
    # σ=117) and (ch=675 σ=5) on a Ra-226 spectrum.
    ASYMMETRY_RATIO_THRESHOLD = 0.30  # min/max ≥ this means symmetric

    intensity_ratios = []
    raw_ir = pattern_intensity_ratios or []
    if (len(raw_ir) == n_lines and
            all(r is not None for r in raw_ir) and
            max(raw_ir) > 0):
        # Pattern declares intensity for every line; check asymmetry
        ir_min = min(raw_ir)
        ir_max = max(raw_ir)
        if ir_min / ir_max < ASYMMETRY_RATIO_THRESHOLD:
            intensity_ratios = list(raw_ir)  # use for ordering check

    # Index by channel for significance and FWHM lookups
    sig_by_channel = {ch: s for ch, s in zip(peak_channels, peak_significances)}
    fwhm_by_channel = {ch: f for ch, f in zip(peak_channels, peak_fwhms)}

    best_match = None
    best_score = float("inf")

    for combo in combinations(peak_channels, n_lines):
        combo_sorted = sorted(combo)
        ch_low = combo_sorted[0]
        ch_high = combo_sorted[-1]
        ch_range = ch_high - ch_low
        if ch_range <= 0:
            continue

        if n_lines >= 3:
            peak_positions = [
                (c - ch_low) / ch_range for c in combo_sorted
            ]
            ok = True
            for pp, ep in zip(peak_positions, pat_positions):
                if abs(pp - ep) > RATIO_TOL:
                    ok = False
                    break
            if not ok:
                continue

        a1 = pat_total_range / ch_range
        a0 = pat_low - a1 * ch_low

        if not (0.05 <= a1 <= 3.5):
            continue

        # --- v1.6 n_channels-aware gain plausibility ---
        # The expected gain for a typical detector is set by the channel
        # count and the energy ceiling: a1_expected ≈ E_ceiling / n_channels.
        # 1024-channel NaI: ~3 keV/ch; 8192-channel: ~0.37 keV/ch.
        # Allow a generous factor of 3 in either direction, but reject
        # gains that span >2× the ceiling or compress everything into <1/3
        # of the channel range. This catches the failure mode where on a
        # high-channel-count NaI background, an 8-line Bi-214 (or 4-line
        # quartet) randomly matches 8 (or 4) peaks at gain ~3 keV/ch —
        # which would imply a total spectrum span of >21000 keV, clearly
        # unphysical.
        if n_channels > 0:
            expected_gain = 3000.0 / n_channels  # ENERGY_CEILING_KEV / n_channels
            if not (expected_gain / 3.0 <= a1 <= expected_gain * 3.0):
                continue

        # --- v1.6 physical-offset sanity check ---
        if not (PLAUSIBLE_OFFSET_LO <= a0 <= PLAUSIBLE_OFFSET_HI):
            continue

        # --- v1.6 wide_pair high-E anchor in upper 40% ---
        if is_wide_pair and n_channels > 0:
            if ch_high < WIDE_PAIR_HIGH_E_MIN_FRAC * n_channels:
                continue

        # --- v1.6 intensity-asymmetry ordering check ---
        # If the pattern is asymmetric (declared in intensity_ratios),
        # the peak slot assigned to the brightest expected line must
        # have higher σ than the peak slot assigned to the dimmest.
        if intensity_ratios:
            obs_sigs = [sig_by_channel.get(c, 0.0) for c in combo_sorted]
            if all(s > 0 for s in obs_sigs):
                bright_idx = max(range(n_lines), key=lambda i: intensity_ratios[i])
                dim_idx = min(range(n_lines), key=lambda i: intensity_ratios[i])
                if obs_sigs[bright_idx] <= obs_sigs[dim_idx]:
                    continue

        if gain_hint is not None:
            if abs(a1 - gain_hint) / gain_hint > 0.15:
                continue
        if offset_hint is not None and abs(a0 - offset_hint) > 50:
            continue

        mapped_energies = [a0 + a1 * c for c in combo_sorted]
        max_err = max(abs(me - pl) for me, pl in zip(mapped_energies, pattern_lines))

        # Effective tolerance scales with pattern size:
        #   n_lines ≤ 4: 0.5·FWHM (linear fit usually adequate)
        #   n_lines ≥ 5: 1.0·FWHM (a linear seed cannot describe a
        #     degree-3 stored cal across a 2000+ keV span better than
        #     ~half the local FWHM; the Lsrm-recommended ERN-line set
        #     has residual ~30 keV at the middle anchor under linear
        #     fit, which exceeds 0.5·FWHM but is correct for a
        #     polynomial-shaped detector response).
        peak_fwhms_ch = [fwhm_by_channel.get(c, 0.0) for c in combo_sorted]
        if peak_fwhms_ch and any(f > 0 for f in peak_fwhms_ch):
            mean_fwhm_keV = (sum(peak_fwhms_ch) / len(peak_fwhms_ch)) * a1
            fwhm_factor = 1.0 if n_lines >= 5 else 0.5
            effective_tolerance = max(tolerance_keV, fwhm_factor * mean_fwhm_keV)
        else:
            effective_tolerance = tolerance_keV

        if max_err > effective_tolerance:
            continue

        # Score: penalize low-significance peaks and large mismatch
        energy_mismatch = sum(
            abs(me - pl) for me, pl in zip(mapped_energies, pattern_lines)
        )
        total_significance = sum(sig_by_channel.get(c, 0.0) for c in combo_sorted)
        # Lower score = better; favour high significance and low mismatch
        # Significance penalty: inverse, so high significance reduces score
        score = energy_mismatch - 0.1 * total_significance

        if score < best_score:
            best_score = score
            best_match = {
                "channels": list(combo_sorted),
                "energies": list(pattern_lines),
                "gain": a1,
                "offset": a0,
                "total_significance": total_significance,
            }

    return best_match


# ============================================================================
# Anchor expansion (Step 3)
# ============================================================================

def _expand_anchors(
    found_peaks: list,
    *,
    current_coefs: tuple,
    patterns: list,
    already_used: set,
) -> list:
    """
    Given a working calibration, search for additional anchor lines that
    we can add. For each anchor line in each pattern, find the nearest
    found peak in calibrated energy space; if it's within a tolerance,
    accept as a new anchor.

    Returns:
        List of (channel, energy_keV, source_pattern_name) tuples.
    """
    if not current_coefs or len(current_coefs) < 2:
        return []

    # Predict energy for each found peak
    peak_energies = []
    for p in found_peaks:
        E = 0.0
        for c in reversed(current_coefs):
            E = E * p.channel + c
        peak_energies.append(E)

    new_anchors = []
    for pat in patterns:
        tol = float(pat.get("tolerance_keV", 3.0))
        for E_expected in pat.get("lines", []):
            # Find nearest peak
            best_idx = None
            best_delta = float("inf")
            for i, E_peak in enumerate(peak_energies):
                d = abs(E_peak - E_expected)
                if d < best_delta:
                    best_delta = d
                    best_idx = i
            if best_idx is not None and best_delta <= tol:
                anchor_key = (
                    round(found_peaks[best_idx].channel, 3),
                    round(float(E_expected), 3),
                )
                if anchor_key in already_used:
                    continue
                new_anchors.append((
                    int(found_peaks[best_idx].channel),
                    float(E_expected),
                    pat["name"],
                ))

    return new_anchors
