"""
Detector type classification.

Per SKILL.md §4, classify the detector by relative resolution R = FWHM/E,
combined with stored peak-shape priors:

  HPGe:        R(662) ~ 0.3-0.5 % (≪ 1%); symmetric Gaussian
  CdZnTe:      R(662) ~ 1.5-2.5 %; left-tail asymmetric (Hypermet)
  LaBr3(Ce):   R(662) ~ 2.5-3 %; symmetric, intrinsic peaks present
  CeBr3:       R(662) ~ 3.5-4 %; symmetric
  NaI(Tl):     R(662) ~ 6-8 %; mildly asymmetric (right-tail)

The classification uses TWO inputs (when available):
  - Measured R at any reference energy (typically ~600-800 keV).
  - Stored peak_type from SimpleSqrtFwhmCalibration:
      0 = pure Gaussian (HPGe, LaBr3, CeBr3)
      1 = Hypermet (NaI with right-tail asymmetry; CdZnTe)

If the spectrum carries the AtomSpectra stored FWHM model with peak_type
and tails, that's a strong prior — but we still verify against the
measured R.

Token economy: returns a compact result dict; no spectrum data.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class DetectorTypeResult:
    """Outcome of detector-type classification."""
    detector_type: str                            # "HPGe", "NaI", "LaBr3", "CeBr3", "CdZnTe", "unknown"
    confidence: float                             # 0.0–1.0
    R_at_662: Optional[float] = None              # R = FWHM/E at 662 keV (or nearest)
    R_reference_energy_keV: Optional[float] = None
    candidates_considered: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    reason: str = ""

    # Suggested FWHM model for this detector (input to fit_fwhm_*)
    suggested_fwhm_model: str = ""                # "hpge" | "scintillator"


# Resolution-range table (R = FWHM/E in fraction)
# Order matters: matched from finest to coarsest resolution
_R_RANGES = [
    ("HPGe",       0.0,   0.010,  "hpge"),
    ("CdZnTe",     0.010, 0.030,  "hpge"),       # use HPGe model (semiconductor)
    ("LaBr3",      0.020, 0.035,  "scintillator"),
    ("CeBr3",      0.030, 0.045,  "scintillator"),
    ("NaI",        0.040, 0.100,  "scintillator"),
]


def classify_detector(
    *,
    fwhm_keV: Optional[float] = None,
    at_energy_keV: Optional[float] = None,
    peak_type_hint: Optional[int] = None,
    stored_fwhm_chi2_per_dof: Optional[float] = None,
) -> DetectorTypeResult:
    """
    Classify the detector type from one (energy, FWHM) measurement plus
    optional priors from the file's stored calibration.

    Args:
        fwhm_keV: measured FWHM in keV at some reference energy
                  (preferably near 662 keV)
        at_energy_keV: the energy where FWHM was measured
        peak_type_hint: AtomSpectra SimpleSqrtFwhmCalibration peak_type
                        (0 = Gaussian, 1 = Hypermet)
        stored_fwhm_chi2_per_dof: quality of the file's FWHM fit

    Returns:
        DetectorTypeResult with the most likely detector type and
        confidence.
    """
    if fwhm_keV is None or at_energy_keV is None or at_energy_keV <= 0:
        return DetectorTypeResult(
            detector_type="unknown",
            confidence=0.0,
            reason="No FWHM measurement provided",
        )

    R = fwhm_keV / at_energy_keV

    # Scale R to its equivalent at 662 keV using approximate √E law
    # (this is rough — full FWHM(E) fit would be better, but here we
    # need a single number for classification)
    R_at_662_approx = R * (at_energy_keV / 662.0) ** 0.5

    # Find candidate detector types whose R range contains R_at_662
    candidates = []
    for name, R_lo, R_hi, fwhm_model in _R_RANGES:
        # Slight expansion to make the brackets soft
        margin = 0.0025
        if R_lo - margin <= R_at_662_approx <= R_hi + margin:
            # Distance from the midpoint
            mid = 0.5 * (R_lo + R_hi)
            distance = abs(R_at_662_approx - mid) / mid if mid > 0 else 1.0
            candidates.append({
                "type": name,
                "R_range": (R_lo, R_hi),
                "distance_to_midpoint": distance,
                "fwhm_model": fwhm_model,
            })

    if not candidates:
        return DetectorTypeResult(
            detector_type="unknown",
            confidence=0.0,
            R_at_662=R_at_662_approx,
            R_reference_energy_keV=at_energy_keV,
            reason=(f"R(662)≈{R_at_662_approx*100:.2f}% doesn't fit any "
                    f"known detector range"),
        )

    # Pick the best candidate by distance to range midpoint
    best = min(candidates, key=lambda c: c["distance_to_midpoint"])

    # Adjust confidence based on:
    #   - distance to midpoint (closer = higher confidence)
    #   - agreement with peak_type_hint:
    #       Hypermet (1) supports NaI / CdZnTe
    #       Gaussian (0) supports HPGe / LaBr3 / CeBr3
    base_confidence = 1.0 - min(1.0, best["distance_to_midpoint"])

    notes = []
    if peak_type_hint is not None:
        if peak_type_hint == 1:
            asymmetric_types = ("NaI", "CdZnTe")
            if best["type"] in asymmetric_types:
                base_confidence = min(1.0, base_confidence + 0.15)
                notes.append("stored peak_type=1 (Hypermet) consistent with "
                             f"{best['type']} asymmetric peak shape")
            else:
                base_confidence *= 0.7
                notes.append(f"stored peak_type=1 (Hypermet) unexpected for "
                             f"{best['type']} — verify")
        elif peak_type_hint == 0:
            symmetric_types = ("HPGe", "LaBr3", "CeBr3")
            if best["type"] in symmetric_types:
                base_confidence = min(1.0, base_confidence + 0.1)
                notes.append("stored peak_type=0 (Gaussian) consistent with "
                             f"{best['type']}")

    if stored_fwhm_chi2_per_dof is not None and stored_fwhm_chi2_per_dof > 5.0:
        notes.append(f"stored FWHM fit χ²/ν = {stored_fwhm_chi2_per_dof:.1f} "
                     f"is high — measured FWHM may be unreliable")

    return DetectorTypeResult(
        detector_type=best["type"],
        confidence=float(base_confidence),
        R_at_662=R_at_662_approx,
        R_reference_energy_keV=at_energy_keV,
        candidates_considered=[c["type"] for c in candidates],
        notes=notes,
        suggested_fwhm_model=best["fwhm_model"],
        reason=(f"R(662)≈{R_at_662_approx*100:.2f}% → {best['type']} "
                f"(range {best['R_range'][0]*100:.1f}–"
                f"{best['R_range'][1]*100:.1f}%)"),
    )
