"""
Secondary peak energies from Compton kinematics and pair production.

For any photopeak at energy E_γ, the detector spectrum also contains
secondary features whose positions are determined by physics — not by
the nuclide library. These features are **diagnostic**: their presence
together with a candidate photopeak corroborates that the photopeak is
real and correctly identified. Their absence (when expected) is a
warning sign.

This module provides closed-form formulas for:

  • Compton edge: maximum energy transferred to an electron by a γ-ray
    scattering at 180° from the detector active volume. Appears as a
    sharp drop on the high-E side of the Compton continuum, ~200 keV
    below the photopeak for typical γ energies.

  • Backscatter peak: energy of the γ-ray that scattered at 180° in
    the source surroundings (collimator, shielding, sample matrix) and
    returned to the detector. E_back + E_compton_edge = E_γ exactly,
    so backscatter sits ~200 keV above zero for typical γ.

  • Single-escape peak (SEP): E_γ − 511 keV. Appears when a γ produces
    an e⁺e⁻ pair inside the detector and one of the two 511-keV
    annihilation γ's escapes. Only present for E_γ > 1022 keV (pair-
    production threshold) and is most prominent for E_γ > 1500 keV.

  • Double-escape peak (DEP): E_γ − 1022 keV. Both annihilation γ's
    escape. Only for E_γ > 1022 keV, prominent above 2 MeV.

  • Annihilation peak: always at 511.00 keV when positron-emitter or
    pair production occurs (epiphenomenon).

  • Sum peaks: E_γ1 + E_γ2 for coincidence-summed γ pairs. Listed only
    for the most-likely pairs (cascades from the same decay).

On NaI 50×50 (FWHM @ 500 keV ≈ 30 keV), Compton-edge and backscatter
features are broad and shifted from the kinematic limit — see
`compton_edge_keV` docstring for the appropriate width. On HPGe (FWHM
@ 500 keV ≈ 1.5 keV), the kinematic limits are sharp and these
features can be located precisely.

Methodology references:
  - Compton scattering kinematics: A.H. Compton, Phys. Rev. 21 (1923) 483
  - Gilmore & Joss "Practical Gamma-ray Spectrometry" 3rd Ed., §2.3.2
  - Knoll "Radiation Detection and Measurement" 4th Ed., §10.II.D

The 511 keV used internally is the electron rest mass energy
m_e·c² = 510.998950 keV (CODATA 2018).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

# Electron rest mass energy — CODATA 2018 value
M_E_C2_KEV = 510.998950

# Pair production threshold = 2·m_e·c²
PAIR_PRODUCTION_THRESHOLD_KEV = 2.0 * M_E_C2_KEV  # 1021.998 keV


@dataclass(frozen=True)
class SecondaryFeatures:
    """Expected positions of all secondary spectrum features for one γ-line."""
    photopeak_keV: float
    compton_edge_keV: float
    backscatter_keV: float
    single_escape_keV: Optional[float]  # None if E_γ < 1022 keV
    double_escape_keV: Optional[float]  # None if E_γ < 1022 keV

    def summary(self) -> str:
        """Human-readable one-line summary."""
        parts = [
            f"photopeak={self.photopeak_keV:.1f}",
            f"Compton_edge={self.compton_edge_keV:.1f}",
            f"backscatter={self.backscatter_keV:.1f}",
        ]
        if self.single_escape_keV is not None:
            parts.append(f"SEP={self.single_escape_keV:.1f}")
        if self.double_escape_keV is not None:
            parts.append(f"DEP={self.double_escape_keV:.1f}")
        return "  ".join(parts)


def compton_edge_keV(photopeak_keV: float) -> float:
    """
    Compton edge energy: maximum energy a γ can transfer to an electron
    by a 180° backscatter inside the detector.

    Formula (from Compton scattering kinematics at θ=180°):
        E_edge = E_γ · (2·E_γ) / (m_e·c² + 2·E_γ)

    Equivalently: E_edge = E_γ − E_backscatter.

    On HPGe this edge is sharp (drop in counts within ~1 FWHM). On NaI
    the broad FWHM smears it into a "Compton shoulder" rather than an
    edge, but the kinematic position is still useful as the upper bound
    of the Compton continuum.

    Args:
        photopeak_keV: E_γ in keV (must be > 0)

    Returns:
        Compton edge energy in keV.

    Examples:
        >>> round(compton_edge_keV(661.66), 1)  # Cs-137
        477.3
        >>> round(compton_edge_keV(1460.82), 1)  # K-40
        1243.4
        >>> round(compton_edge_keV(2614.51), 1)  # Tl-208
        2381.8
    """
    if photopeak_keV <= 0:
        raise ValueError(f"photopeak_keV must be > 0, got {photopeak_keV}")
    E = float(photopeak_keV)
    return E * (2.0 * E) / (M_E_C2_KEV + 2.0 * E)


def backscatter_keV(photopeak_keV: float) -> float:
    """
    Backscatter peak energy: energy of a γ that scattered at 180° in
    the source surroundings and returned to the detector.

    Formula:
        E_back = E_γ − E_Compton_edge
               = E_γ · m_e·c² / (m_e·c² + 2·E_γ)

    For high-E γ (E_γ ≫ m_e·c²), E_back → m_e·c²/2 = 256 keV — that is,
    the backscatter peak from any high-energy line converges around
    250 keV, which is why a broad "backscatter bump" near 200–250 keV
    is a generic feature of γ-ray spectra in real environments.

    Args:
        photopeak_keV: E_γ in keV

    Returns:
        Backscatter peak energy in keV.

    Examples:
        >>> round(backscatter_keV(661.66), 1)  # Cs-137
        184.3
        >>> round(backscatter_keV(1460.82), 1)  # K-40
        217.5
        >>> round(backscatter_keV(2614.51), 1)  # Tl-208
        232.8
    """
    if photopeak_keV <= 0:
        raise ValueError(f"photopeak_keV must be > 0, got {photopeak_keV}")
    E = float(photopeak_keV)
    return E * M_E_C2_KEV / (M_E_C2_KEV + 2.0 * E)


def single_escape_keV(photopeak_keV: float) -> Optional[float]:
    """
    Single-escape peak (SEP): E_γ − 511 keV.

    Appears when E_γ > 1022 keV (pair production threshold) and one of
    the two 511-keV annihilation γ's escapes from the detector before
    being absorbed. SEP becomes prominent for E_γ above ~1500 keV.

    Args:
        photopeak_keV: E_γ in keV

    Returns:
        SEP energy in keV, or None if E_γ ≤ pair production threshold.

    Examples:
        >>> single_escape_keV(661.66) is None  # below threshold
        True
        >>> round(single_escape_keV(2614.51), 1)  # Tl-208 SEP
        2103.5
    """
    if photopeak_keV <= PAIR_PRODUCTION_THRESHOLD_KEV:
        return None
    return float(photopeak_keV) - M_E_C2_KEV


def double_escape_keV(photopeak_keV: float) -> Optional[float]:
    """
    Double-escape peak (DEP): E_γ − 1022 keV.

    Appears when E_γ > 1022 keV and BOTH annihilation γ's escape from
    the detector. DEP becomes prominent for E_γ above ~2 MeV.

    Args:
        photopeak_keV: E_γ in keV

    Returns:
        DEP energy in keV, or None if E_γ ≤ pair production threshold.

    Examples:
        >>> double_escape_keV(661.66) is None  # below threshold
        True
        >>> round(double_escape_keV(2614.51), 1)  # Tl-208 DEP
        1592.5
    """
    if photopeak_keV <= PAIR_PRODUCTION_THRESHOLD_KEV:
        return None
    return float(photopeak_keV) - PAIR_PRODUCTION_THRESHOLD_KEV


def secondary_features(photopeak_keV: float) -> SecondaryFeatures:
    """
    Compute all secondary spectrum features for one photopeak.

    Useful for cross-checking a candidate identification: if the
    photopeak is real, the spectrum should contain at least the
    Compton edge and backscatter at their kinematic positions (the
    SEP/DEP appear only above 1022 keV).

    Args:
        photopeak_keV: E_γ in keV

    Returns:
        SecondaryFeatures with kinematic positions of all features.

    Example:
        >>> sf = secondary_features(1460.82)  # K-40
        >>> round(sf.compton_edge_keV, 1)
        1243.4
        >>> round(sf.backscatter_keV, 1)
        217.5
        >>> round(sf.single_escape_keV, 1)
        949.8
        >>> round(sf.double_escape_keV, 1)
        438.8
    """
    return SecondaryFeatures(
        photopeak_keV=float(photopeak_keV),
        compton_edge_keV=compton_edge_keV(photopeak_keV),
        backscatter_keV=backscatter_keV(photopeak_keV),
        single_escape_keV=single_escape_keV(photopeak_keV),
        double_escape_keV=double_escape_keV(photopeak_keV),
    )


def sum_peak_keV(E1_keV: float, E2_keV: float) -> float:
    """
    Sum (coincidence-summing) peak energy: E1 + E2.

    Appears when two cascade γ's from the same decay are detected
    simultaneously within the resolving time of the spectrometer. Sum
    peaks are most prominent for short-lived intermediate nuclear
    levels and at close source-to-detector geometries.

    Args:
        E1_keV, E2_keV: energies of the two cascading γ's

    Returns:
        Sum peak energy in keV.

    Examples:
        >>> round(sum_peak_keV(1173.23, 1332.49), 2)  # Co-60 sum
        2505.72
    """
    return float(E1_keV) + float(E2_keV)


__all__ = [
    "M_E_C2_KEV", "PAIR_PRODUCTION_THRESHOLD_KEV",
    "SecondaryFeatures",
    "compton_edge_keV", "backscatter_keV",
    "single_escape_keV", "double_escape_keV",
    "secondary_features", "sum_peak_keV",
]
