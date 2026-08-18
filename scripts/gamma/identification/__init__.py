"""
Nuclide identification per Lsrm Algorithmic Foundations §6, §14.

Modules:
  - window: identification window δE(E) for matching peaks to library
    lines, with detector-appropriate scaling (sqrt_E for scintillators,
    linear for HPGe).
  - mda: Minimum Detectable Activity, decision threshold, and detection
    limit per ISO 11929:2019.
  - confidence: Confidence Index per Lsrm §14.3 — discriminates strong
    multi-line identifications from weak single-line guesses.
  - identify: main Lsrm-style algorithm "find peaks for nuclide lines",
    starting from the characteristic line of each candidate.
"""

from gamma.identification.window import (
    REFERENCE_ENERGY_KEV, DEFAULT_DELTA_E0_KEV,
    IdentificationWindow,
    build_identification_window,
    identification_window_from_fwhm,
    build_id_window_k_fwhm,
)
# F-167 — canonical k·FWHM(E) ID window constants & helpers.
from gamma.identification.id_window import (
    DetectorClass,
    ID_WINDOW_K_FWHM,
    id_window_keV,
    normalize_detector_class,
)
from gamma.identification.mda import (
    K_ALPHA_95, K_ALPHA_99, K_ALPHA_999,
    MdaResult, mda_for_peak, characteristic_line_of_nuclide,
)
from gamma.identification.confidence import (
    ConfidenceIndexResult, confidence_index,
)
from gamma.identification.identify import (
    LineMatch, NuclideIdentification, IdentificationResult,
    identify_nuclides,
)
from gamma.identification.cross_check import (
    SecondaryFeatureCheck, check_secondary_features,
    cross_check_identification,
)
from gamma.identification.disambiguate import (
    NATURAL_CHAIN_NUCLIDES, POSITRON_EMITTERS_NEAR_511, NAI_CONFUSION_MAP,
    disambiguate_identifications,
)
from gamma.identification.proportionality import (
    RARE_ISOTOPE_PRIOR, DEFAULT_PRIOR, RATIO_TOLERANCE_FACTOR,
    ProportionalityCheckResult,
    get_prior, check_intensity_proportionality,
)
from gamma.identification.chain_equilibrium import (
    RA_226_CHAIN_GROUPS, TH_232_CHAIN_GROUPS,
    ChainEquilibriumResult, check_ra226_chain_equilibrium,
)

__all__ = [
    "REFERENCE_ENERGY_KEV", "DEFAULT_DELTA_E0_KEV",
    "IdentificationWindow",
    "build_identification_window",
    "identification_window_from_fwhm",
    "build_id_window_k_fwhm",
    "DetectorClass",
    "ID_WINDOW_K_FWHM",
    "id_window_keV",
    "normalize_detector_class",
    "K_ALPHA_95", "K_ALPHA_99", "K_ALPHA_999",
    "MdaResult", "mda_for_peak", "characteristic_line_of_nuclide",
    "ConfidenceIndexResult", "confidence_index",
    "LineMatch", "NuclideIdentification", "IdentificationResult",
    "identify_nuclides",
    "SecondaryFeatureCheck", "check_secondary_features",
    "cross_check_identification",
    "NATURAL_CHAIN_NUCLIDES", "POSITRON_EMITTERS_NEAR_511", "NAI_CONFUSION_MAP",
    "disambiguate_identifications",
    "RARE_ISOTOPE_PRIOR", "DEFAULT_PRIOR", "RATIO_TOLERANCE_FACTOR",
    "ProportionalityCheckResult",
    "get_prior", "check_intensity_proportionality",
    "RA_226_CHAIN_GROUPS", "TH_232_CHAIN_GROUPS",
    "ChainEquilibriumResult", "check_ra226_chain_equilibrium",
]
