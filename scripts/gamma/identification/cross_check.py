"""
Secondary-feature cross-check for nuclide identification.

For each candidate photopeak, the spectrum should contain corroborating
secondary features at kinematically determined positions:

  - Compton edge at E_γ · (2·E_γ) / (m_e·c² + 2·E_γ)
  - Backscatter peak at E_γ − E_Compton_edge
  - For E_γ > 1022 keV: SEP at E_γ − 511, DEP at E_γ − 1022

If these features ARE present at the right energies, confidence in the
photopeak identification is boosted. If they are ABSENT (with sufficient
spectrum statistics to have seen them), the identification is suspect.

This complements the Confidence Index (Lsrm §14.3) which counts library
lines but does not look at physics-based secondaries.

The cross-check works on already-identified nuclides — it does NOT
discover new ones. Its role is to PROMOTE or DEMOTE confidence in the
identifications made by `identify_nuclides`.

Implementation considerations:
  - Compton edges and backscatter peaks are broad (typically 2-3× FWHM
    at the edge energy) because the kinematic edge is smeared by
    multiple-scattering and detector response.
  - The Compton continuum is high under these features → significance
    threshold for "feature present" is lower than for photopeaks.
  - Only strong photopeaks (σ > 30) have detectable secondaries on
    typical detectors.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from gamma.physics.secondary import (
    compton_edge_keV, backscatter_keV,
    single_escape_keV, double_escape_keV,
    PAIR_PRODUCTION_THRESHOLD_KEV,
)


@dataclass(frozen=True)
class SecondaryFeatureCheck:
    """Result of checking expected secondary features for a photopeak."""

    photopeak_E_keV: float
    photopeak_sigma: float

    # Expected positions
    expected_compton_edge_keV: float
    expected_backscatter_keV: float
    expected_SEP_keV: Optional[float]
    expected_DEP_keV: Optional[float]

    # Detection: did we find a peak/excess at each expected position?
    found_compton_edge: bool
    found_backscatter: bool
    found_SEP: bool
    found_DEP: bool

    # Boost: how much the secondary features add to confidence (0..1)
    confidence_boost: float

    notes: str = ""


def check_secondary_features(
    photopeak_E_keV: float,
    photopeak_sigma: float,
    found_peaks: list,
    spec,
    *,
    feature_window_keV: float = 30.0,
    min_sigma_for_check: float = 30.0,
) -> SecondaryFeatureCheck:
    """
    For one candidate photopeak, verify presence of expected secondaries.

    Args:
        photopeak_E_keV: library energy of the candidate photopeak
        photopeak_sigma: significance of the matched peak in the spectrum
        found_peaks: list of all FoundPeaks in the spectrum
        spec: Spectrum object (for channel-to-energy mapping)
        feature_window_keV: half-width of search window for each
            secondary feature, in keV. Larger than typical
            identification windows because Compton features are broad
            and SEP/DEP are spread by multiple scattering.
        min_sigma_for_check: only run this check if photopeak σ exceeds
            this threshold (weaker photopeaks don't produce detectable
            secondaries).

    Returns:
        SecondaryFeatureCheck.
    """
    # Compute expected positions
    expected_compton = compton_edge_keV(photopeak_E_keV)
    expected_backscatter = backscatter_keV(photopeak_E_keV)
    expected_SEP = single_escape_keV(photopeak_E_keV)
    expected_DEP = double_escape_keV(photopeak_E_keV)

    # If photopeak is too weak, secondaries are too dim to see
    if photopeak_sigma < min_sigma_for_check:
        return SecondaryFeatureCheck(
            photopeak_E_keV=photopeak_E_keV,
            photopeak_sigma=photopeak_sigma,
            expected_compton_edge_keV=expected_compton,
            expected_backscatter_keV=expected_backscatter,
            expected_SEP_keV=expected_SEP,
            expected_DEP_keV=expected_DEP,
            found_compton_edge=False,
            found_backscatter=False,
            found_SEP=False,
            found_DEP=False,
            confidence_boost=0.0,
            notes=f"Photopeak σ={photopeak_sigma:.1f} < {min_sigma_for_check} threshold; "
                  f"secondaries not searched.",
        )

    # Map found peaks to energies (use spec's stored calibration)
    peak_energies = [
        (p, spec.channel_to_energy(p.channel))
        for p in found_peaks
    ]

    def find_peak_near(E_target_keV: float) -> bool:
        """True if a peak exists within feature_window_keV of E_target."""
        for p, E in peak_energies:
            if abs(E - E_target_keV) <= feature_window_keV:
                # Significance threshold for "feature present" is lower
                # — Compton/backscatter peaks usually have σ ≈ 5–20
                if p.significance >= 3.0:
                    return True
        return False

    found_compton = find_peak_near(expected_compton)
    found_backscatter = find_peak_near(expected_backscatter)
    found_SEP = (expected_SEP is not None and find_peak_near(expected_SEP))
    found_DEP = (expected_DEP is not None and find_peak_near(expected_DEP))

    # Confidence boost: each detected feature adds +0.25 (max 1.0).
    # Compton edge and backscatter are the strongest indicators.
    boost = 0.0
    if found_compton:
        boost += 0.30
    if found_backscatter:
        boost += 0.30
    if found_SEP:
        boost += 0.20
    if found_DEP:
        boost += 0.20
    boost = min(1.0, boost)

    return SecondaryFeatureCheck(
        photopeak_E_keV=photopeak_E_keV,
        photopeak_sigma=photopeak_sigma,
        expected_compton_edge_keV=expected_compton,
        expected_backscatter_keV=expected_backscatter,
        expected_SEP_keV=expected_SEP,
        expected_DEP_keV=expected_DEP,
        found_compton_edge=found_compton,
        found_backscatter=found_backscatter,
        found_SEP=found_SEP,
        found_DEP=found_DEP,
        confidence_boost=boost,
    )


def cross_check_identification(
    identification_result,
    found_peaks: list,
    spec,
    *,
    feature_window_keV: float = 30.0,
) -> dict:
    """
    Run secondary-feature cross-check on all detected nuclides.

    For each detected nuclide, check its characteristic line for
    secondary features. Returns a dict mapping nuclide name to a
    SecondaryFeatureCheck.

    This adds physics-based corroboration to the library-line-only
    Confidence Index from `identify_nuclides`.

    Args:
        identification_result: IdentificationResult from identify_nuclides
        found_peaks: all peaks from mariscotti_search
        spec: Spectrum object
        feature_window_keV: search window for secondaries

    Returns:
        dict {nuclide_name: SecondaryFeatureCheck}
    """
    results = {}
    for ni in identification_result.detected_nuclides:
        if not ni.matched_lines:
            continue
        # Check the characteristic line (the strongest matched)
        char_match = None
        for m in ni.matched_lines:
            if m.is_characteristic:
                char_match = m
                break
        if char_match is None:
            continue
        check = check_secondary_features(
            photopeak_E_keV=char_match.library_E_keV,
            photopeak_sigma=char_match.significance_currie or 0.0,
            found_peaks=found_peaks,
            spec=spec,
            feature_window_keV=feature_window_keV,
        )
        results[ni.nuclide] = check
    return results


# G2 / v1.31.2 -- intensity-ratio chi^2 gate (annotation, no decision impact).
#
# Spec: AUDIT_F-419_skill_vs_canon.md / G2 -- per-nuclide consistency check
# across matched lines. For each detected nuclide with >=2 matched lines the
# gate computes
#
#     Q_i        = peak_area_i / I_gamma_i   (counts per unit %, proxy of A*t*eps)
#     sigma_Q_i  = peak_area_unc_i / I_gamma_i
#     Q_mean     = sum(w_i * Q_i) / sum(w_i),   w_i = 1 / sigma_Q_i^2
#     chi^2     = sum(w_i * (Q_i - Q_mean)^2)
#     ndof       = n_lines - 1
#     ratio      = chi^2 / ndof
#
# Verdict (thresholds locked by audit):
#     ratio <= 1.5   -> "pass_strict"
#     ratio <= 3.0   -> "pass_lenient"
#     ratio  > 3.0   -> "fail"
#
# This is an annotation only: callers do not down-weight identification on
# "fail" -- they surface the value in the JSON / HTML report for operator
# triage. Without efficiency-curve normalisation the gate is biased toward
# "fail" when the matched lines span a wide energy range on a detector with
# steep eps(E) (NaI below 200 keV). Documented limitation -- operator chooses
# the mode.
INTENSITY_RATIO_CHI2_STRICT_THRESHOLD = 1.5
INTENSITY_RATIO_CHI2_LENIENT_THRESHOLD = 3.0


def _line_area_and_sigma(line_match, fallback_sigma_floor: float = 1.0):
    """Return (area, sigma_area) for one LineMatch, or (None, None).

    Prefer explicit peak_area + peak_area_uncertainty. When uncertainty is
    missing fall back to area / significance_currie (Currie-style proxy).
    """
    area = getattr(line_match, "peak_area", None)
    try:
        area_f = float(area) if area is not None else None
    except (TypeError, ValueError):
        area_f = None
    if area_f is None or area_f <= 0:
        return None, None
    unc = getattr(line_match, "peak_area_uncertainty", None)
    try:
        unc_f = float(unc) if unc is not None else None
    except (TypeError, ValueError):
        unc_f = None
    if unc_f is not None and unc_f > 0:
        return area_f, unc_f
    sig = getattr(line_match, "significance_currie", None)
    try:
        sig_f = float(sig) if sig is not None else None
    except (TypeError, ValueError):
        sig_f = None
    if sig_f is not None and sig_f > 0:
        # area / sigma ~ Currie significance -> sigma = area / significance
        return area_f, area_f / sig_f
    return area_f, max(fallback_sigma_floor, area_f ** 0.5)


def intensity_ratio_chi2_gate(
    identification_result,
    *,
    min_lines: int = 2,
    min_intensity_pct: float = 1.0,
) -> dict:
    """Compute per-nuclide intensity-ratio chi^2 gate.

    Args:
        identification_result: IdentificationResult from identify_nuclides.
        min_lines: minimum number of usable lines required to compute the
            gate (default 2 -- below that ndof <= 0).
        min_intensity_pct: drop lines with library_I_pct below this floor;
            tiny intensities blow up Q_i = area/I and dominate chi^2.

    Returns:
        dict mapping nuclide name -> dict with keys
            chi2          : float
            ndof          : int
            ratio         : float (chi2 / ndof)
            verdict       : "pass_strict" / "pass_lenient" / "fail"
            n_lines_used  : int
            n_lines_skipped : int
            q_mean        : float
            lines         : list of {E_keV, I_pct, area, sigma_area, q, q_residual_sigma}
        Nuclides without enough lines are omitted (caller can detect absence).
    """
    out: dict = {}
    if identification_result is None:
        return out
    detected = getattr(identification_result, "detected_nuclides", ()) or ()
    for ni in detected:
        nuc = getattr(ni, "nuclide", None)
        if not nuc:
            continue
        matched = list(getattr(ni, "matched_lines", ()) or ())
        usable = []
        n_skipped = 0
        for m in matched:
            I_pct = float(getattr(m, "library_I_pct", 0.0) or 0.0)
            if I_pct < min_intensity_pct:
                n_skipped += 1
                continue
            area, sigma_area = _line_area_and_sigma(m)
            if area is None or sigma_area is None or sigma_area <= 0:
                n_skipped += 1
                continue
            E_lib = float(getattr(m, "library_E_keV", 0.0) or 0.0)
            q_i = area / I_pct
            sigma_q = sigma_area / I_pct
            usable.append({
                "E_keV": E_lib,
                "I_pct": I_pct,
                "area": area,
                "sigma_area": sigma_area,
                "q": q_i,
                "sigma_q": sigma_q,
            })
        n_used = len(usable)
        if n_used < min_lines:
            continue
        weights = [1.0 / (u["sigma_q"] ** 2) for u in usable]
        w_sum = sum(weights)
        if w_sum <= 0:
            continue
        q_mean = sum(w * u["q"] for w, u in zip(weights, usable)) / w_sum
        chi2 = sum(w * (u["q"] - q_mean) ** 2 for w, u in zip(weights, usable))
        ndof = n_used - 1
        ratio = chi2 / ndof if ndof > 0 else float("inf")
        if ratio <= INTENSITY_RATIO_CHI2_STRICT_THRESHOLD:
            verdict = "pass_strict"
        elif ratio <= INTENSITY_RATIO_CHI2_LENIENT_THRESHOLD:
            verdict = "pass_lenient"
        else:
            verdict = "fail"
        for u in usable:
            u["q_residual_sigma"] = (u["q"] - q_mean) / u["sigma_q"] if u["sigma_q"] > 0 else None
        out[nuc] = {
            "chi2": chi2,
            "ndof": ndof,
            "ratio": ratio,
            "verdict": verdict,
            "n_lines_used": n_used,
            "n_lines_skipped": n_skipped,
            "q_mean": q_mean,
            "lines": usable,
        }
    return out


# Extend __all__ -- preserved original names + new gate.
__all__ = ["SecondaryFeatureCheck", "check_secondary_features",
           "cross_check_identification",
           "intensity_ratio_chi2_gate",
           "INTENSITY_RATIO_CHI2_STRICT_THRESHOLD",
           "INTENSITY_RATIO_CHI2_LENIENT_THRESHOLD"]