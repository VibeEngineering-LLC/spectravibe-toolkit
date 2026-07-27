"""
Multi-line intensity-ratio proportionality check and rare-isotope
priors for identification disambiguation.

For any nuclide with multiple γ-lines, the relative intensities of
its lines are fixed by nuclear physics (branching ratios from
ENSDF). When a candidate nuclide is matched to spectrum peaks, the
ratios of those peak areas should match the library intensity ratios
within a tolerance set by detector efficiency variation across the
energy range.

If the proportionality check fails for a multi-line nuclide, the
identification is suspect: either the peaks don't really belong to
this nuclide, or there is an interfering nuclide contributing to one
or more of the peaks.

Expert methodology (Lsrm + practitioner consensus):

  • For close-by lines (within ~factor 2 in energy), ε(E) is nearly
    constant, so peak-σ ratios should match library I ratios within
    factor ~2-3.

  • For widely separated lines (>factor 5 in energy), ε(E) variation
    on NaI is significant; tolerance increases.

  • Rare/anthropogenic isotopes (Zn-65, Co-60 in environmental
    samples, Cs-134) have very low prior probability — they need
    MULTIPLE corroborating lines with proportional intensities
    before being accepted.

  • Natural U is always a mixture of U-238 (99.27%) and U-235
    (0.72%). When U is genuinely present, characteristic lines from
    BOTH isotopes appear with proportionality determined by isotopic
    abundance (enrichment level varies — depleted/natural/enriched
    U have different U-235/U-238 ratios).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple, Callable


# Rare-isotope prior: multiplicative penalty applied to the
# Confidence Index when this isotope is identified, reflecting its
# low natural occurrence. A penalty of 0.1 effectively requires the
# identification to be 10x more confident than for common isotopes
# before being accepted.
#
# Calibration:
#   1.0  = common, no penalty (e.g. K-40, Cs-137, Ra-chain, Th-chain)
#   0.5  = uncommon (e.g. Co-60 may appear from industrial/medical)
#   0.2  = rare (Ba-133 contamination, Eu-152 calibration sources)
#   0.05 = very rare (Zn-65 fallout product — almost never seen)
RARE_ISOTOPE_PRIOR = {
    "Zn-65": 0.05,       # very rare — fallout-era only
    "Cs-134": 0.3,       # reactor fission product, rare in env. samples
    "Co-58": 0.1,        # activation product, very rare environmentally
    "Mn-54": 0.1,        # activation product, very rare
    "Co-60": 0.5,        # industrial/medical, occasional
    "Eu-152": 0.5,       # calibration source; wide E-range (121-1408 keV) means
                         # raw-sigma proportionality test unreliable without ε-correction
    "Eu-154": 0.2,       # calibration source isotope
    "Ba-133": 0.2,       # calibration source isotope
    "Na-22": 0.1,        # cosmogenic, very rare in normal environment
    "Am-241": 0.3,       # smoke-detector / industrial
    "Be-7":   0.5,       # cosmogenic, occasionally on filters
    "Sb-125": 0.2,       # rare fission product
    "Ru-106": 0.2,       # rare fission product
    "Ag-110m": 0.2,      # rare activation/fission product
    "I-131":  0.2,       # short-lived medical
    "I-133":  0.1,       # very short-lived medical
    "Co-57":  0.2,       # calibration source; essentially absent in env/LSRM samples
    "In-111": 0.1,       # short-lived medical diagnostic; never in env samples
    # natural & common: default 1.0
}

# Default prior for nuclides not in the table
DEFAULT_PRIOR = 1.0

# Multi-line intensity-ratio tolerance multiplicative factor.
# Observed ratio must satisfy:
#   library_ratio / ratio_tol  ≤  observed_ratio  ≤  library_ratio · ratio_tol
RATIO_TOLERANCE_FACTOR = 3.0

# For lines separated by more than this energy factor, ε(E) variation
# matters; widen the tolerance further.
WIDE_ENERGY_FACTOR = 5.0
WIDE_RATIO_TOLERANCE_FACTOR = 5.0


@dataclass(frozen=True)
class ProportionalityCheckResult:
    """Result of an intensity-ratio proportionality check."""
    nuclide: str
    n_lines_checked: int
    n_ratios_passed: int
    n_ratios_failed: int
    passed: bool   # True if ≥ minimum fraction of ratios are consistent
    reason: str
    failed_pairs: tuple = ()   # tuple of (E1, E2, observed, expected) for failures


def get_prior(nuclide: str) -> float:
    """Return the rare-isotope prior weight (multiplier on CI)."""
    return RARE_ISOTOPE_PRIOR.get(nuclide, DEFAULT_PRIOR)


def check_intensity_proportionality(
    nuclide: str,
    matched_lines: list,
    *,
    min_lines_required: int = 2,
    ratio_tolerance: float = RATIO_TOLERANCE_FACTOR,
    min_intensity_threshold_pct: float = 1.0,
    use_peak_area_when_available: bool = True,
    efficiency_curve: Optional[Callable] = None,
) -> ProportionalityCheckResult:
    """
    Check whether matched peaks of a candidate nuclide have intensities
    proportional to library values.

    Three modes of operation (in priority order):

      1. EFFICIENCY-CORRECTED (most accurate): if `efficiency_curve` is
         provided AND all matched lines have valid peak_area, compare
         (area / (ε(E) · I)) ratios instead of bare amplitude ratios.
         This is the activity ratio — should equal 1.0 ± counting
         statistics for a single nuclide. Tolerance: factor 2 (tight).
         **Closes K-12** (widely-separated lines).

      2. AREA-BASED (good for close-by lines): if all matched lines
         have peak_area but no efficiency_curve, use (area / I) ratios.
         Tolerance: factor 3 default. May fail for lines >5× apart in
         energy due to ε(E) variation.

      3. SIGMA-PROXY (fallback): use significance_currie as amplitude proxy.
         Same tolerance scaling as area mode but biased by baseline
         continuum shape.

    Args:
        nuclide: nuclide name (for diagnostics)
        matched_lines: list of LineMatch with library_E_keV,
            library_I_pct, significance_currie, optionally peak_area
        min_lines_required: minimum matched lines for the check
        ratio_tolerance: multiplicative deviation tolerance (default 3)
        min_intensity_threshold_pct: skip lines below this I%
        use_peak_area_when_available: if True (default), prefer area
            over σ when available
        efficiency_curve: callable(E_keV) → ε(E) for efficiency
            correction. If None, no correction applied.

    Returns:
        ProportionalityCheckResult.
    """
    # Determine which amplitude metric to use
    use_area = False
    if use_peak_area_when_available:
        all_areas_sane = all(
            getattr(m, "peak_area", None) is not None
            and (getattr(m, "peak_area", 0) or 0) > 0
            for m in matched_lines
        )
        use_area = all_areas_sane and len(matched_lines) > 0

    # Efficiency correction is only meaningful when we have areas
    use_eff = use_area and efficiency_curve is not None

    def _normalised_amplitude(m):
        """
        Return amplitude normalised by intensity (and efficiency if
        available). For a single nuclide, this should be the same value
        across all lines (= nuclide activity × live_time, up to constants).
        """
        if use_area:
            amp = float(m.peak_area)
        else:
            amp = float(m.significance_currie or 0.0)
        if m.library_I_pct <= 0:
            return None
        # Normalise by intensity
        amp_per_I = amp / m.library_I_pct
        if use_eff:
            eps = efficiency_curve(m.library_E_keV)
            if eps is None or eps <= 0:
                return None
            return amp_per_I / eps
        return amp_per_I

    # Filter to lines with library intensity above threshold and positive amplitude
    relevant = [
        m for m in matched_lines
        if m.library_I_pct >= min_intensity_threshold_pct
        and (float(m.peak_area) if use_area else float(m.significance_currie or 0.0)) > 0
    ]

    if len(relevant) < min_lines_required:
        return ProportionalityCheckResult(
            nuclide=nuclide,
            n_lines_checked=len(relevant),
            n_ratios_passed=0,
            n_ratios_failed=0,
            passed=True,
            reason=f"Only {len(relevant)} matched lines with I≥{min_intensity_threshold_pct}%; "
                   f"need ≥{min_lines_required} for proportionality check. Defer judgement.",
        )

    # When efficiency-corrected, the metric is "activity-equivalent" and
    # should be constant across lines. Use tighter tolerance.
    eff_corrected_tolerance = 2.0  # factor 2 for activity-ratio mode

    n_passed = 0
    n_failed = 0
    failed_pairs = []

    for i in range(len(relevant)):
        for j in range(i + 1, len(relevant)):
            m_i = relevant[i]
            m_j = relevant[j]
            amp_i = _normalised_amplitude(m_i)
            amp_j = _normalised_amplitude(m_j)
            if amp_i is None or amp_j is None or amp_j <= 0:
                continue

            if use_eff:
                # Activity-ratio mode: amp_i and amp_j are both
                # activity-equivalent; their ratio should be 1.0
                observed_ratio = amp_i / amp_j
                expected_ratio = 1.0
                tol = eff_corrected_tolerance
            else:
                # Legacy mode: amp_i = (raw_amp_i / I_i), amp_j similar.
                # Their ratio equals raw amplitude ratio (since I drops
                # out — but we've still normalised explicitly for clarity)
                observed_ratio = amp_i / amp_j
                expected_ratio = 1.0  # since we normalised by I
                # Use wide tolerance for distant energies
                E_ratio = max(m_i.library_E_keV, m_j.library_E_keV) / \
                          min(m_i.library_E_keV, m_j.library_E_keV)
                tol = (WIDE_RATIO_TOLERANCE_FACTOR
                       if E_ratio > WIDE_ENERGY_FACTOR
                       else ratio_tolerance)

            ratio_min = expected_ratio / tol
            ratio_max = expected_ratio * tol

            if ratio_min <= observed_ratio <= ratio_max:
                n_passed += 1
            else:
                n_failed += 1
                failed_pairs.append((
                    m_i.library_E_keV, m_j.library_E_keV,
                    observed_ratio, expected_ratio,
                ))

    total = n_passed + n_failed
    if total == 0:
        return ProportionalityCheckResult(
            nuclide=nuclide,
            n_lines_checked=len(relevant),
            n_ratios_passed=0,
            n_ratios_failed=0,
            passed=True,
            reason="No usable ratio pairs.",
        )
    pass_fraction = n_passed / total
    passed_overall = pass_fraction >= 0.60

    proxy_str = "area/ε/I (activity ratio)" if use_eff else \
                ("area/I" if use_area else "σ/I (σ proxy)")
    return ProportionalityCheckResult(
        nuclide=nuclide,
        n_lines_checked=len(relevant),
        n_ratios_passed=n_passed,
        n_ratios_failed=n_failed,
        passed=passed_overall,
        reason=(f"{n_passed}/{total} ratio pairs proportional via "
                f"{proxy_str} (tol {eff_corrected_tolerance if use_eff else ratio_tolerance:.0f}×); "
                f"{'PASS' if passed_overall else 'FAIL'} (need ≥60%)"),
        failed_pairs=tuple(failed_pairs),
    )


__all__ = [
    "RARE_ISOTOPE_PRIOR", "DEFAULT_PRIOR",
    "RATIO_TOLERANCE_FACTOR",
    "ProportionalityCheckResult",
    "get_prior",
    "check_intensity_proportionality",
]
