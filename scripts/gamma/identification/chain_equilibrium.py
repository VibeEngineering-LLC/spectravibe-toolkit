"""
Decay-chain equilibrium model for natural-decay chains.

For natural-decay chains (Ra-226, Th-232, U-238, U-235), parent and
daughter nuclides reach secular equilibrium only if the chain is
sealed. In practice, equilibrium is often broken:

  • Ra-226 → Rn-222 (gaseous!) → Po-218 → Pb-214 → Bi-214 → Pb-210
    Rn-222 (T½ = 3.8 d) can escape from porous samples, fluids, soils.
    When Rn escapes, the Pb-214 / Bi-214 / Po-214 daughters that
    follow it in the chain are depleted; Ra-226 (186 keV) itself is
    NOT affected by Rn loss.
    Pb-210 (T½ = 22.3 y) is a long-lived terminal-region daughter
    that integrates over decades of past Rn emanation. It does NOT
    track current Ra-226 activity in environmental samples.

  • Th-232 chain: no gaseous intermediate of consequence on usual
    time scales; secular equilibrium normally holds except in fresh
    chemical separations.

  • U-238 chain: Ra-226 itself is a chain member. The Pb-210 problem
    above applies.

  • Pb-210 has another, completely independent source: **the lead
    metal of detector shielding**. Standard Pb has 0.1-100 Bq/kg of
    Pb-210 contamination depending on its age and ore source. The
    46.5 keV line seen in low-background spectra is usually from
    the shielding, NOT from the sample.

Therefore, for identification of Ra-226 in a sample, the relevant
proportionality groups are:

  GROUP A — directly tied to Ra-226 activity:
    Ra-226 186.21 keV (intrinsic)

  GROUP B — short-lived Rn-222 daughters, affected by Rn loss:
    Pb-214 295.22, 351.93 keV
    Bi-214 609.31, 768.36, 934.06, 1120.29, 1238.11, 1377.67,
           1764.49, 2204.21 keV

  GROUP C — long-lived Pb-210 (NOT usable for sample Ra-226):
    Pb-210 46.5 keV  ← typically dominated by shielding/contamination

The proportionality check must be done WITHIN each group, not across
groups. Ratios within Group B remain fixed (those nuclides are in
short-term secular equilibrium with each other even if Rn escapes —
once Rn decays in place its daughters track each other). But the
ratio of Group A to Group B reflects Rn retention, which varies.

This module:
  • Defines the chain groups
  • Provides a per-group proportionality check
  • Reports the inferred Rn-retention fraction (if both Ra-226 186
    and Bi-214 lines are detected with proportionality holding
    within each group)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from gamma.identification.proportionality import (
    check_intensity_proportionality, ProportionalityCheckResult,
)


# Ra-226 chain groups for proportionality. Each group is a list of
# (nuclide_name, line_energy_keV).
RA_226_CHAIN_GROUPS = {
    # Group A: Ra-226 intrinsic
    "Ra-226_intrinsic": [
        ("Ra-226", 186.21),
    ],
    # Group B: short-lived Rn-222 daughters
    # Within this group, lines are always in mutual equilibrium because
    # half-lives are short (min–hours) vs the integration time of a
    # gamma measurement (hours–days). Rn loss reduces all of Group B
    # uniformly relative to Group A but preserves intra-group ratios.
    "Rn222_daughters": [
        ("Pb-214", 295.22),
        ("Pb-214", 351.93),
        ("Bi-214", 609.31),
        ("Bi-214", 1120.29),
        ("Bi-214", 1764.49),
    ],
    # Group C: Pb-210 — has independent source (shielding) and is NOT
    # used for Ra-226 identification in samples. Listed for awareness.
    "Pb210_long_lived": [
        ("Pb-210", 46.54),
    ],
}

TH_232_CHAIN_GROUPS = {
    # Th-232 chain typically reaches secular equilibrium; treat as one
    # group, but Tl-208 is a branched daughter (Bi-212 → Tl-208 36%) so
    # its lines have a known branching factor applied.
    "Th232_daughters": [
        ("Ac-228", 911.16),
        ("Ac-228", 968.97),
        ("Pb-212", 238.63),
        ("Bi-212", 727.33),
        ("Tl-208", 583.19),
        ("Tl-208", 2614.51),
    ],
}


@dataclass(frozen=True)
class ChainEquilibriumResult:
    """Result of a decay-chain proportionality analysis."""
    chain_name: str
    group_results: dict   # group_name → ProportionalityCheckResult
    # Inferred Rn-retention fraction (Group A / expected from Group B).
    # 1.0 = full equilibrium; < 1.0 = Rn escape; > 1.0 = Ra excess
    # (recent introduction).
    rn_retention_ratio: Optional[float]
    # Diagnostic: is the chain identification CONSISTENT (each group
    # proportional internally, even if cross-group equilibrium broken)
    chain_consistent: bool
    notes: str = ""


def check_ra226_chain_equilibrium(
    matched_lines_by_nuclide: dict,
) -> ChainEquilibriumResult:
    """
    Check Ra-226 chain identification with disequilibrium awareness.

    Args:
        matched_lines_by_nuclide: dict {nuclide_name: list_of_LineMatch}

    Returns:
        ChainEquilibriumResult with per-group proportionality and the
        inferred Rn retention.

    Algorithm:
      1. For each group (intrinsic / Rn222_daughters / Pb210),
         collect the matched lines that belong to it.
      2. Check intensity-ratio proportionality WITHIN each group.
      3. Group A is just 1 line (Ra-226 186) — always trivially passes.
         Group B is the substantive proportionality check.
      4. If both Group A and Group B are populated, compute the
         Rn-retention ratio = (Group A peak σ) / (Group B mean σ
         normalised by library intensity).
      5. Overall chain is CONSISTENT iff each populated group's
         internal proportionality passes.

    Note: Pb-210 is NOT used to validate or invalidate Ra-226
    identification because Pb-210 in real spectra is dominated by
    shielding contamination, not by sample Ra-226 daughters.
    """
    group_results = {}

    # --- Group A: Ra-226 intrinsic ---
    ra226_matches = matched_lines_by_nuclide.get("Ra-226", [])
    group_results["Ra-226_intrinsic"] = check_intensity_proportionality(
        "Ra-226_intrinsic", ra226_matches,
        min_lines_required=2,  # only 1 line in group — defer
    )

    # --- Group B: Rn-222 daughters (Pb-214, Bi-214) ---
    # Aggregate matched lines from Pb-214 and Bi-214 into a single list
    # for the proportionality check — they're in mutual equilibrium.
    rn_daughter_matches = []
    for nuc in ("Pb-214", "Bi-214"):
        rn_daughter_matches.extend(matched_lines_by_nuclide.get(nuc, []))
    group_results["Rn222_daughters"] = check_intensity_proportionality(
        "Rn222_daughters", rn_daughter_matches,
        min_lines_required=2,
    )

    # --- Group C: Pb-210 — informational only ---
    # We don't validate the chain by Pb-210 because its source is
    # typically the lead shielding around the detector.
    pb210_matches = matched_lines_by_nuclide.get("Pb-210", [])
    if pb210_matches:
        group_results["Pb210_long_lived"] = ProportionalityCheckResult(
            nuclide="Pb-210_long_lived",
            n_lines_checked=len(pb210_matches),
            n_ratios_passed=0,
            n_ratios_failed=0,
            passed=True,  # not used for chain validation
            reason=("Pb-210 detected but NOT used for chain validation: "
                    "its source in low-background spectra is typically "
                    "shielding contamination, not sample Ra-226 daughters."),
        )

    # --- Cross-group: Rn-retention diagnostic ---
    # Use peak_area when populated (Phase 2.1a), else significance_currie (BUG-34 Phase 3a').
    def _amplitude(m):
        a = getattr(m, "peak_area", None)
        if a is not None and a > 0:
            return float(a)
        return float(m.significance_currie or 0.0)

    rn_retention = None
    if ra226_matches and rn_daughter_matches:
        # Take Ra-226 186 (Group A)
        ra_amp_normalised = None
        for m in ra226_matches:
            if abs(m.library_E_keV - 186.21) < 2.0:
                if m.library_I_pct > 0:
                    ra_amp_normalised = _amplitude(m) / m.library_I_pct
                break

        # Take Bi-214 609 (typical Group B reference)
        bi214_amp_normalised = None
        for m in rn_daughter_matches:
            if abs(m.library_E_keV - 609.31) < 4.0:
                if m.library_I_pct > 0:
                    bi214_amp_normalised = _amplitude(m) / m.library_I_pct
                break

        if ra_amp_normalised is not None and bi214_amp_normalised is not None \
                and bi214_amp_normalised > 0:
            rn_retention = ra_amp_normalised / bi214_amp_normalised

    # Chain is consistent if all populated groups pass internally
    chain_consistent = all(
        r.passed for r in group_results.values()
        if r.n_lines_checked >= 2  # only check groups with enough data
    )

    return ChainEquilibriumResult(
        chain_name="Ra-226",
        group_results=group_results,
        rn_retention_ratio=rn_retention,
        chain_consistent=chain_consistent,
        notes=(
            "Ra-226 chain validated by within-group proportionality. "
            "Cross-group ratios may be broken due to Rn-222 escape "
            "(reducing Group B relative to Group A). Pb-210 is NOT used "
            "for validation because of independent shielding-contamination "
            "source."
        ),
    )


__all__ = [
    "RA_226_CHAIN_GROUPS", "TH_232_CHAIN_GROUPS",
    "ChainEquilibriumResult",
    "check_ra226_chain_equilibrium",
]
