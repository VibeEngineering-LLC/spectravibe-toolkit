"""
Physics models for γ-ray spectrometry:

- secondary: Compton edge, backscatter, escape peaks, sum peaks
  (kinematic positions of secondary spectrum features used for
  identification cross-checks)
- pileup: random pile-up (rate ∝ cps²) and cascade-coincidence
  summing (rate ∝ activity, geometry-dependent) detection
"""
from gamma.physics.secondary import (
    M_E_C2_KEV, PAIR_PRODUCTION_THRESHOLD_KEV,
    SecondaryFeatures,
    compton_edge_keV, backscatter_keV,
    single_escape_keV, double_escape_keV,
    secondary_features, sum_peak_keV,
)
from gamma.physics.pileup import (
    KNOWN_CASCADE_PAIRS, PileupCandidate, detect_pileup_peaks,
)

__all__ = [
    "M_E_C2_KEV", "PAIR_PRODUCTION_THRESHOLD_KEV",
    "SecondaryFeatures",
    "compton_edge_keV", "backscatter_keV",
    "single_escape_keV", "double_escape_keV",
    "secondary_features", "sum_peak_keV",
    "KNOWN_CASCADE_PAIRS", "PileupCandidate", "detect_pileup_peaks",
]
