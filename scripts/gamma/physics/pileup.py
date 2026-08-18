"""
Detection of pile-up (random summing) and cascade (true coincidence)
summing peaks in γ-ray spectra.

Two physically distinct phenomena produce "extra" peaks at energies
that are sums of other photopeak energies. Misidentifying them as
new nuclides is a common pitfall.

═══════════════════════════════════════════════════════════════════════
1. Random pile-up (pulse pile-up)
═══════════════════════════════════════════════════════════════════════

Two INDEPENDENT γ-rays arrive at the detector within the resolving
time τ of the electronics; their pulses overlap and the system
records them as a single event with summed amplitude.

  • Rate proportional to (total_cps)² · τ
  • Affects ALL energy pairs in the spectrum (regardless of origin)
  • Most prominent at high count rates (>several kcps total)
  • Independent of source nuclide structure
  • Position: E_pileup ≈ E1 + E2 for any two coincident γ's

Diagnostic: pile-up peaks appear at sums of the strongest photopeaks
when the total count rate is high. A peak at exactly 2·E (twice the
photopeak energy) — "self-pile-up" — is the textbook signature.

Typical resolving time for NaI + commodity electronics: τ ≈ 1-3 μs.
For τ=2 μs and total_cps=1000, random pile-up rate per pair is
about (1000)² · 2e-6 = 2 cps — non-trivial.

═══════════════════════════════════════════════════════════════════════
2. Cascade (true coincidence) summing
═══════════════════════════════════════════════════════════════════════

Two γ-rays from the SAME decay event (cascade transitions of one
nucleus) are detected in coincidence because the intermediate
nuclear level has T½ ≪ detector resolving time.

  • Rate proportional to activity × (Ω/4π)² for two-detector cases,
    or just (Ω/4π) for self-coincidence in one detector
  • Affects only γ pairs from cascade transitions
  • Most prominent at large solid angles (Marinelli, close geom.)
  • Independent of total spectrum count rate
  • Position: E_cascade = E1 + E2 for specific known cascades

Examples on common reference sources:
  • Co-60: 1173.23 + 1332.49 → 2505.72 keV (always present in close
    geometry; absent at far geometry)
  • Y-88: 898.04 + 1836.06 → 2734.10 keV
  • Th-232 chain Tl-208: 583.19 + 2614.51 → 3197.70 keV
  • Cs-134: 569 + 605 → 1174; 605 + 796 → 1401; 605 + 802 → 1407
  • Eu-152: many cascades; 122 + 344 → 466; 244 + 1408 → 1652 etc.

Cascade-sum effect REDUCES the apparent area of the contributing
photopeaks (their counts redistribute into the sum peak). Above ~1%
solid angle this becomes significant for activity calculations.

═══════════════════════════════════════════════════════════════════════
Disambiguation
═══════════════════════════════════════════════════════════════════════

If a peak at E_sum = E1+E2 is observed and both E1 and E2 photopeaks
are present:

  • At LOW count rate (total cps < 200) and CLOSE geometry → most
    likely cascade summing
  • At HIGH count rate (total cps > 1000) and ANY geometry → pile-up
    dominates; cascade contribution also possible if cascade pair
  • At HIGH count rate AND far geometry → primarily pile-up (cascade
    contribution suppressed by small solid angle)

The two effects can coexist. They are physically separable by:
  • Reducing source-detector distance (boosts cascade, weak effect on
    pile-up)
  • Reducing activity (proportional drop for cascade, quadratic drop
    for pile-up)

Reference: Knoll "Radiation Detection" 4th Ed., Ch. 17 (NaI summing);
Lsrm Algorithmic Foundations §10 (coincidence corrections).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


# Known cascade-summing pairs in commonly-used nuclides.
# (nuclide, E1_keV, E2_keV, E_sum_keV)
KNOWN_CASCADE_PAIRS = [
    # Co-60: classic example
    ("Co-60", 1173.23, 1332.49, 2505.72),
    # Y-88
    ("Y-88", 898.04, 1836.06, 2734.10),
    # Th-232 chain (Tl-208 has β-cascade through 2614 + 583)
    ("Tl-208", 583.19, 2614.51, 3197.70),
    ("Tl-208", 510.77, 2614.51, 3125.28),
    ("Tl-208", 510.77, 583.19, 1093.96),    # often seen in Th-232 Marinelli
    ("Tl-208", 277.36, 583.19, 860.55),
    ("Tl-208", 583.19, 763.13, 1346.32),
    # Th-228 chain Bi-212 has cascade pairs
    ("Bi-212", 727.33, 1620.50, 2347.83),
    # Ac-228 has many cascade pairs in Th chain
    ("Ac-228", 338.32, 911.20, 1249.52),
    ("Ac-228", 338.32, 968.97, 1307.29),
    ("Ac-228", 911.20, 968.97, 1880.17),
    # Cs-134 (β-decay populated cascade)
    ("Cs-134", 569.32, 604.72, 1174.04),
    ("Cs-134", 604.72, 795.86, 1400.58),
    ("Cs-134", 604.72, 802.04, 1406.76),
    ("Cs-134", 569.32, 795.86, 1365.18),
    # Eu-152 (many cascades; only the dominant ones listed)
    ("Eu-152", 121.78, 344.28, 466.06),
    ("Eu-152", 121.78, 1408.01, 1529.79),
    ("Eu-152", 244.70, 1408.01, 1652.71),
    ("Eu-152", 344.28, 1408.01, 1752.29),
    # Na-22 (positron + 1274 cascade — actual coincidence of 2 annihilation γ
    # and the 1274 prompt)
    ("Na-22", 511.0, 1274.54, 1785.54),
    # Ra-226 chain Bi-214 has multiple cascade pairs
    ("Bi-214", 609.31, 1120.29, 1729.60),
    ("Bi-214", 609.31, 1764.49, 2373.80),
    # Th-228 chain Pb-212 → Bi-212 cascade
    ("Pb-212", 238.63, 300.09, 538.72),
]


# Map chain-member nuclide name (in Lsrm library) to the actual daughter
# nuclide that emits the cascade γ. E.g. Th-232 library entry includes
# Tl-208 lines because Tl-208 is a daughter in the chain; the cascade
# pair "Tl-208 583+2614" should be detected even when our library lists
# only "Th-232" with the combined chain.
CHAIN_DAUGHTER_NAMES = {
    "Th-232": ["Tl-208", "Pb-212", "Bi-212", "Ac-228"],
    "Th-228": ["Tl-208", "Pb-212", "Bi-212"],
    "Ra-226": ["Pb-214", "Bi-214", "Pb-210"],
    "U-238": ["Th-234", "Pa-234m", "Ra-226", "Pb-214", "Bi-214"],
    "U-235": ["Th-231"],
}


@dataclass(frozen=True)
class PileupCandidate:
    """Detected potential pile-up or sum peak."""
    type: str                       # "random_pileup" or "cascade_sum" or "ambiguous"
    candidate_E_keV: float          # observed peak energy
    parent_E1_keV: float            # first contributing γ
    parent_E2_keV: float            # second contributing γ
    parent_nuclide: Optional[str]   # for cascade only; None for random
    expected_E_sum_keV: float       # E1 + E2 (or 2·E1 for self-pile)
    energy_residual_keV: float      # |observed - expected|
    cps_total: Optional[float]      # total cps if computed
    reason: str

    @property
    def is_self_pileup(self) -> bool:
        return abs(self.parent_E1_keV - self.parent_E2_keV) < 1.0


def detect_pileup_peaks(
    *,
    found_peaks: list,
    spec,
    detected_nuclides: Optional[list] = None,
    pileup_window_keV: float = 30.0,
    high_cps_threshold: float = 500.0,
    min_parent_sigma: float = 10.0,
) -> list:
    """
    Examine detected peaks for potential pile-up / cascade-sum signatures.

    Args:
        found_peaks: peaks from mariscotti_search (with channel + significance)
        spec: Spectrum (for channel→keV mapping and live_time)
        detected_nuclides: list of NuclideIdentification (from
            identify_nuclides). If provided, cascade summing is checked
            against known cascades of these nuclides.
        pileup_window_keV: tolerance for matching observed peak to a
            predicted sum. On NaI 50×50 FWHM @ 2500 keV ≈ 80 keV, so
            30 keV is appropriate; tighter for HPGe.
        high_cps_threshold: total cps above which random pile-up
            becomes a serious candidate explanation. Below this,
            random pile-up is unlikely.
        min_parent_sigma: only consider photopeaks with σ above this
            as parents (weak peaks don't generate detectable pile-up).

    Returns:
        list of PileupCandidate.

    Algorithm:
      1. Compute total cps from spectrum
      2. Build set of "parent" photopeaks (σ ≥ min_parent_sigma)
      3. For each pair (E1, E2) of parent peaks, predict E_sum = E1+E2
      4. Scan found_peaks for any near E_sum within pileup_window_keV
      5. Classify each match:
         - cascade_sum: if (parent_nuclide, E1, E2) matches a known
           cascade pair from KNOWN_CASCADE_PAIRS
         - random_pileup: if total cps ≥ high_cps_threshold
         - ambiguous: cannot distinguish
    """
    total_cps = spec.counts.sum() / spec.live_time if spec.live_time > 0 else 0.0

    # Build parent peaks with energy
    parents = []
    for p in found_peaks:
        if p.significance < min_parent_sigma:
            continue
        E = spec.channel_to_energy(p.channel)
        if E > 0:
            parents.append((E, p))

    if len(parents) < 1:
        return []

    # Build set of found peak energies for sum matching
    found_E = [(spec.channel_to_energy(p.channel), p) for p in found_peaks]

    # Known cascade pairs keyed by detected nuclide (case-insensitive set)
    detected_names = set()
    if detected_nuclides is not None:
        detected_names = {ni.nuclide for ni in detected_nuclides}
        # Also include chain daughters: if "Th-232" detected, treat
        # Tl-208/Pb-212/Bi-212 as if they were detected too (they're
        # in the chain and contribute lines to the same spectrum)
        for parent, daughters in CHAIN_DAUGHTER_NAMES.items():
            if parent in detected_names:
                detected_names.update(daughters)

    candidates = []

    # --- Cascade-sum candidates ---
    for nuc, E1, E2, E_sum in KNOWN_CASCADE_PAIRS:
        if detected_names and nuc not in detected_names:
            continue
        # Both parents must be present in found_peaks (within FWHM)
        parent1_present = any(
            abs(E - E1) < pileup_window_keV / 2 for E, _ in found_E
        )
        parent2_present = any(
            abs(E - E2) < pileup_window_keV / 2 for E, _ in found_E
        )
        if not (parent1_present and parent2_present):
            continue
        # Look for a peak near E_sum
        for E_obs, p in found_E:
            if abs(E_obs - E_sum) <= pileup_window_keV:
                candidates.append(PileupCandidate(
                    type="cascade_sum",
                    candidate_E_keV=E_obs,
                    parent_E1_keV=E1,
                    parent_E2_keV=E2,
                    parent_nuclide=nuc,
                    expected_E_sum_keV=E_sum,
                    energy_residual_keV=abs(E_obs - E_sum),
                    cps_total=total_cps,
                    reason=(f"Cascade pair {nuc} {E1:.1f}+{E2:.1f} keV → "
                            f"sum at {E_sum:.1f} keV; observed peak at "
                            f"{E_obs:.1f} keV (Δ={abs(E_obs - E_sum):.1f})"),
                ))
                break

    # --- Random pile-up candidates ---
    # Only relevant when total cps is significant
    if total_cps >= high_cps_threshold:
        # Sort parents by significance descending — pile-up is dominated
        # by the strongest peaks (rate ∝ σ_i · σ_j)
        parents_sorted = sorted(parents, key=lambda x: -x[1].significance)
        # Consider top-5 strongest as candidates for pile-up parents
        for i, (E1, p1) in enumerate(parents_sorted[:5]):
            for j, (E2, p2) in enumerate(parents_sorted[:5]):
                if j < i:
                    continue
                E_sum_expected = E1 + E2
                # Look for a found peak near this sum
                for E_obs, p_obs in found_E:
                    if abs(E_obs - E_sum_expected) <= pileup_window_keV:
                        # Skip if this was already classified as cascade
                        if any(c.candidate_E_keV == E_obs and c.type == "cascade_sum"
                               for c in candidates):
                            continue
                        candidates.append(PileupCandidate(
                            type="random_pileup",
                            candidate_E_keV=E_obs,
                            parent_E1_keV=E1,
                            parent_E2_keV=E2,
                            parent_nuclide=None,
                            expected_E_sum_keV=E_sum_expected,
                            energy_residual_keV=abs(E_obs - E_sum_expected),
                            cps_total=total_cps,
                            reason=(f"High cps={total_cps:.0f} (>{high_cps_threshold}); "
                                    f"sum of strong photopeaks {E1:.1f}+{E2:.1f} → "
                                    f"{E_sum_expected:.1f} keV matches observed "
                                    f"peak at {E_obs:.1f} keV"),
                        ))
                        break

    return candidates


__all__ = [
    "KNOWN_CASCADE_PAIRS",
    "PileupCandidate",
    "detect_pileup_peaks",
]
