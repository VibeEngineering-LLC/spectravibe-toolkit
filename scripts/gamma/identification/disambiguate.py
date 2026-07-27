"""
Identification disambiguation per Lsrm Algorithmic Foundations §14.4
(Принцип наименьшего действия).

After `identify_nuclides` runs, the result often contains spurious
detections where library lines of different nuclides fall within the
same identification window on a low-resolution detector. Classic
examples on NaI 50×50:

  • 511 keV peak: matches both Na-22 (real annihilation) and Tl-208
    510.77 (Th chain). Tl-208 is much more common in natural
    background; Na-22 would require positron-emitting contamination.

  • 604.7 keV peak: matches both Cs-134 (604.72 keV main line) and
    Bi-214 (609.31 keV, Ra chain). On NaI these are 4.5 keV apart vs
    ~30 keV FWHM — unresolved.

  • 1115.5 keV peak: matches both Zn-65 (1115.55 keV) and Bi-214
    (1120.29 keV, Ra chain). Same problem — 4.7 keV apart on NaI.

  • 605 keV in Cs-134 also has secondary lines at 569, 796, 802 keV;
    if only the 605 line is matched, the identification is suspect.

Lsrm §14.4 "Principle of least action": prefer the identification
that requires the FEWEST assumptions and the SIMPLEST nuclide
inventory. If natural background nuclides (Ra, Th chains) already
explain the spectrum, do NOT add a positron emitter / artificial
nuclide to explain the same peaks.

Disambiguation rules implemented here:

  Rule 1 — "Strong chain trumps single line"
  If a multi-line nuclide of a natural-decay chain is detected with
  ≥2 matched lines, and a single-line nuclide claims the SAME peak,
  the single-line claim is removed.

  Rule 2 — "Natural background trumps positron emitter"
  If both Tl-208 (Th chain) and Na-22 claim the 511 keV region, and
  Tl-208 has additional matched lines (583 or 2614), the Na-22 claim
  is removed (the 511 ROI counts come from Tl-208 510.77 plus
  annihilation, indistinguishable on NaI).

  Rule 3 — "Higher-CI nuclide wins shared peaks"
  When two nuclides match the same peak, the one with higher CI keeps
  it; the other has that match removed. If it then has zero matched
  lines remaining, the nuclide is rejected entirely.

  Rule 4 — "Cross-check boost as tiebreaker"
  If two nuclides both have moderate-high CI and share a peak, the one
  with higher secondary-feature confirmation boost keeps the peak.

  Rule 5 — "Secondary-feature anti-misidentification" (F-40 / v1.7.18)
  If every matched line of a candidate nuclide falls inside the empirical
  position range (p10..p90) of a Compton edge, backscatter, escape peak,
  or other secondary feature of an ALREADY-DETECTED parent — the candidate
  is most likely a misattribution of the parent's continuum. Demote it.
  Classic case: Ac-228 911.20 keV ↔ Co-60 1173 keV Compton edge
  (observed range [906.9..912.5] on NaI 63×63). Catalog comes from
  `detectors/Gamma-1S/data/secondary_peaks_v2.json` (9 problem isotopes,
  Gamma-1S-specific) and is loaded via
  `gamma.physics.secondary_peaks.matches_secondary`.

This module operates on `IdentificationResult` to produce a refined
result.

Reference: Lsrm Algorithmic Foundations 2022 §14.3, §14.4;
Knoll §10 (Compton scattering), §11.A.5 (backscatter peaks).
"""

from __future__ import annotations

from dataclasses import replace
from typing import Optional

from gamma.identification.identify import (
    IdentificationResult, NuclideIdentification, LineMatch,
)
from gamma.identification.proportionality import (
    check_intensity_proportionality, get_prior,
    RARE_ISOTOPE_PRIOR,
)
from gamma.identification.chain_equilibrium import (
    check_ra226_chain_equilibrium,
)
from gamma.physics.secondary_peaks import (
    matches_secondary, load_catalog_v2,
)


# Nuclides that are part of natural decay chains and should be
# proportionality-checked AS A CHAIN rather than individually. The
# Ra-226 chain in particular has the radon-loss disequilibrium
# problem and needs group-wise checking via check_ra226_chain_equilibrium.
CHAIN_NUCLIDES_EXEMPT_FROM_INDIVIDUAL_PROP = {
    "Ra-226", "Pb-214", "Bi-214", "Pb-210",
}


# Sets of nuclides that are members of natural decay chains. When a
# natural-chain nuclide is detected with ≥2 lines, single-line claims
# on the same peaks from non-chain nuclides are suspect.
NATURAL_CHAIN_NUCLIDES = {
    "Ra-226", "Pb-214", "Bi-214", "Pb-210",
    "Th-232", "Ac-228", "Pb-212", "Bi-212", "Tl-208",
    "U-235", "U-238",
    "K-40",  # primordial
}

# Anthropogenic / positron-emitter nuclides that are commonly
# mis-identified as Tl-208 510.77 (annihilation peak).
POSITRON_EMITTERS_NEAR_511 = {
    "Na-22",  # also has 1274.5 keV secondary line — check that too
    "F-18",
    "Co-58",
    "Ga-68",
}

# Nuclides commonly confused with Ra-chain on NaI 50×50.
# Key: spurious nuclide, Value: list of natural-chain nuclides whose
# detection (≥2 lines) suggests the spurious one is the misattribution.
NAI_CONFUSION_MAP = {
    "Cs-134": ["Bi-214"],  # 604.72 vs 609.31
    "Zn-65":  ["Bi-214"],  # 1115.55 vs 1120.29
    "Na-22":  ["Tl-208"],  # 511.00 vs 510.77
    "Ba-133": ["Eu-152"],  # 356.01 vs Eu-152 344 keV region (overlap with cal drift)
}

# Pairs where the spurious nuclide's CHARACTERISTIC LINE shares the SAME
# physical peak channel as the chain nuclide's matched peak. Classic case:
# Co-57 (122.06 keV) and Eu-152 (121.78 keV) — 0.28 keV apart, same NaI
# channel. The chain nuclide's direct detection explains the spurious
# char peak even when the spurious has >1 matched lines (extra matches
# are secondary-feature coincidences or noise peaks with no real source).
#
# Conditions to reject spurious:
#   • chain has ≥ MIN_CHAIN_LINES_CHAR_OVERLAP matched lines
#   • chain CI > spurious CI
#   • spurious char-line peak channel ∈ chain's matched peak channels
NAI_CHAR_OVERLAP_PAIRS: dict = {
    "Co-57":  ["Eu-152"],  # 122.06 keV ≈ Eu-152 121.78 keV (same NaI channel)
    "In-111": ["Eu-152"],  # 245.35 keV ≈ Eu-152 244.70 keV; 171 keV = secondary/noise
}
MIN_CHAIN_LINES_CHAR_OVERLAP = 4  # chain must be strongly multi-line to override

# Indicators of the Ra-226 chain (any 2 of these detected → Ra chain
# is present, so 186 keV ROI is Ra-226 not U-235).
RA_CHAIN_INDICATORS = {"Pb-214", "Bi-214", "Pb-210"}


# Library intensities for the Ra-226/U-235 disambiguation.
# Per Lsrm methodology, U-235 must show 143.76 keV (10.96%) together
# with 185.72 keV (57.0%) in proportional intensity. Ra-226 186.21 keV
# (3.59%) must be proportional to its chain partners 351.93 keV
# (Pb-214, 35.60%) and 609.31 keV (Bi-214, 45.49%).
#
# Expected intensity ratios (in counts/s, eff factors approximately
# cancel for close energies):
#   U-235:        I(143.76) / I(185.72)  = 10.96 / 57.0   ≈ 0.192
#   Ra-226 chain: I(186.21) / I(351.93)  =  3.59 / 35.6  ≈ 0.101
#   Ra-226 chain: I(186.21) / I(609.31)  =  3.59 / 45.49 ≈ 0.079
RATIO_U235_143_OVER_186 = 10.96 / 57.0   # ≈ 0.192
RATIO_RA226_186_OVER_351 = 3.59 / 35.60  # ≈ 0.101
RATIO_RA226_186_OVER_609 = 3.59 / 45.49  # ≈ 0.079

# Acceptable ratio deviation (multiplicative): observed ratio must be
# within [expected/RATIO_TOL, expected·RATIO_TOL]. Factor 3 is generous
# but reflects ε(E) variation on NaI between e.g. 144 and 186 keV.
RATIO_TOLERANCE_FACTOR = 3.0


def disambiguate_identifications(
    result: IdentificationResult,
    *,
    min_chain_lines_to_override: int = 2,
    min_ci_difference_to_resolve: float = 1.0,
    apply_secondary_feature_rule: bool = True,
    secondary_max_lines: int = 2,
) -> IdentificationResult:
    """
    Apply Lsrm-style disambiguation rules to refine identification.

    Args:
        result: IdentificationResult from identify_nuclides
        min_chain_lines_to_override: a natural-chain nuclide needs at
            least this many matched lines (default 2) to override a
            single-line claim on shared peaks.
        min_ci_difference_to_resolve: when two non-chain nuclides
            share a peak, the higher-CI one wins only if its CI
            exceeds the other by at least this amount. Below this
            threshold both keep the match (ambiguous).
        apply_secondary_feature_rule: enable Rule 5 / F-40 — demote a
            candidate when ALL of its matched lines fall inside the
            empirical secondary-feature ranges (Compton edge,
            backscatter, escape) of an already-detected parent in the
            `secondary_peaks_v2` catalog. Default True.
        secondary_max_lines: the secondary-feature rule only acts on
            candidates with ≤ this many matched lines (default 2) —
            stronger multi-line evidence overrides the rule.

    Returns:
        Refined IdentificationResult with spurious nuclides moved to
        the rejected list.
    """
    detected = list(result.detected_nuclides)
    rejected = list(result.rejected_nuclides)
    notes_lines = []

    # Build a map: peak_channel → list of (nuclide_name, LineMatch) using it
    peak_users: dict = {}
    for ni in detected:
        for m in ni.matched_lines:
            peak_users.setdefault(m.peak_channel, []).append((ni.nuclide, m))

    # Build name → NuclideIdentification map for quick lookup
    ni_by_name = {ni.nuclide: ni for ni in detected}

    # --- Rule 2: Na-22 / positron-emitter ↔ Tl-208 ---
    # If Tl-208 detected with multiple lines and Na-22 claims 511,
    # remove Na-22.
    if "Tl-208" in ni_by_name:
        tl208 = ni_by_name["Tl-208"]
        if len(tl208.matched_lines) >= min_chain_lines_to_override:
            for emitter in POSITRON_EMITTERS_NEAR_511:
                if emitter in ni_by_name:
                    em = ni_by_name[emitter]
                    em_ci = em.confidence.CI if em.confidence else 0
                    # Only override if emitter has ≤1 matched line and low CI
                    if len(em.matched_lines) <= 1 and em_ci < 3.0:
                        notes_lines.append(
                            f"  Removed {emitter}: 511 ROI explained by "
                            f"Tl-208 510.77 ({len(tl208.matched_lines)} chain lines)"
                        )
                        # Move to rejected
                        new_rejected = NuclideIdentification(
                            nuclide=emitter, detected=False,
                            reason=(f"Superseded by Tl-208 chain identification "
                                    f"({len(tl208.matched_lines)} matched lines): "
                                    f"Lsrm §14.4 principle of least action."),
                            characteristic_line_keV=em.characteristic_line_keV,
                            matched_lines=(),
                        )
                        rejected.append(new_rejected)
                        detected = [d for d in detected if d.nuclide != emitter]
                        del ni_by_name[emitter]

    # --- Rule 1+2 generalised: NAI_CONFUSION_MAP ---
    for spurious, chain_nuclides in NAI_CONFUSION_MAP.items():
        if spurious not in ni_by_name:
            continue
        # Check if any chain nuclide is detected with ≥ threshold lines
        for chain in chain_nuclides:
            if chain in ni_by_name:
                chain_ni = ni_by_name[chain]
                if len(chain_ni.matched_lines) >= min_chain_lines_to_override:
                    spur = ni_by_name[spurious]
                    spur_ci = spur.confidence.CI if spur.confidence else 0
                    chain_ci = chain_ni.confidence.CI if chain_ni.confidence else 0
                    if (len(spur.matched_lines) <= 1
                            and spur_ci < chain_ci):
                        notes_lines.append(
                            f"  Removed {spurious} ({len(spur.matched_lines)} line, CI={spur_ci:.2f}): "
                            f"superseded by {chain} ({len(chain_ni.matched_lines)} lines, CI={chain_ci:.2f})"
                        )
                        new_rejected = NuclideIdentification(
                            nuclide=spurious, detected=False,
                            reason=(f"Superseded by {chain} natural-chain identification "
                                    f"({len(chain_ni.matched_lines)} matched lines, "
                                    f"CI={chain_ci:.2f} > {spur_ci:.2f}): "
                                    f"Lsrm §14.4 principle of least action."),
                            characteristic_line_keV=spur.characteristic_line_keV,
                            matched_lines=(),
                        )
                        rejected.append(new_rejected)
                        detected = [d for d in detected if d.nuclide != spurious]
                        if spurious in ni_by_name:
                            del ni_by_name[spurious]
                        break  # done with this spurious

    # --- Rule 1c: NAI_CHAR_OVERLAP_PAIRS ---
    # Spurious nuclide's characteristic line physically shares the SAME peak
    # channel as the chain nuclide's matched peak.  Because NaI cannot resolve
    # lines < ~2 keV apart, the chain nuclide directly explains the spurious
    # char peak even if the spurious has accumulated >1 matched lines (those
    # extra matches are assumed to be secondary-feature coincidences or noise
    # peaks with no real source in the sample).  We require the chain nuclide
    # to be strongly multi-line (≥ MIN_CHAIN_LINES_CHAR_OVERLAP) and to have
    # higher CI than the spurious nuclide.
    for _spur_name, _chain_names in NAI_CHAR_OVERLAP_PAIRS.items():
        if _spur_name not in ni_by_name:
            continue
        for _chain_name in _chain_names:
            if _chain_name not in ni_by_name:
                continue
            _chain_ni = ni_by_name[_chain_name]
            if len(_chain_ni.matched_lines) < MIN_CHAIN_LINES_CHAR_OVERLAP:
                continue
            _spur_ni = ni_by_name[_spur_name]
            _spur_ci = _spur_ni.confidence.CI if _spur_ni.confidence else 0
            _chain_ci = _chain_ni.confidence.CI if _chain_ni.confidence else 0
            if _spur_ci >= _chain_ci:
                continue
            _spur_char_keV = _spur_ni.characteristic_line_keV or 0.0
            _chain_ch_set = {m.peak_channel for m in _chain_ni.matched_lines}
            _char_shared = any(
                m.peak_channel in _chain_ch_set
                for m in _spur_ni.matched_lines
                if abs(m.library_E_keV - _spur_char_keV) < 5.0
            )
            if _char_shared:
                notes_lines.append(
                    f"  Removed {_spur_name} (char {_spur_char_keV:.2f} keV, "
                    f"{len(_spur_ni.matched_lines)} lines, CI={_spur_ci:.2f}): "
                    f"char peak explained by {_chain_name} "
                    f"({len(_chain_ni.matched_lines)} lines, CI={_chain_ci:.2f})"
                )
                _rej = NuclideIdentification(
                    nuclide=_spur_name, detected=False,
                    reason=(
                        f"Characteristic peak {_spur_char_keV:.2f} keV directly "
                        f"explained by {_chain_name} identification "
                        f"({len(_chain_ni.matched_lines)} matched lines, "
                        f"CI={_chain_ci:.2f} > {_spur_ci:.2f}): "
                        f"Lsrm §14.4 — principle of least action."
                    ),
                    characteristic_line_keV=_spur_ni.characteristic_line_keV,
                    matched_lines=(),
                )
                rejected.append(_rej)
                detected = [d for d in detected if d.nuclide != _spur_name]
                if _spur_name in ni_by_name:
                    del ni_by_name[_spur_name]
                break  # done with this spurious nuclide

    # --- Rule 4: Universal intensity-ratio proportionality check ---
    # For every multi-line nuclide that is NOT part of a natural decay
    # chain (those are checked by chain-specific equilibrium logic),
    # verify that observed peak σ ratios match library intensity ratios.
    nuclides_failing_proportionality = []
    for ni in list(detected):
        if len(ni.matched_lines) < 2:
            continue
        if ni.nuclide in CHAIN_NUCLIDES_EXEMPT_FROM_INDIVIDUAL_PROP:
            continue  # checked by chain equilibrium below
        prop = check_intensity_proportionality(
            ni.nuclide, list(ni.matched_lines),
            min_lines_required=2,
            min_intensity_threshold_pct=1.0,
        )
        if not prop.passed:
            nuclides_failing_proportionality.append((ni, prop))

    for ni, prop in nuclides_failing_proportionality:
        prior = get_prior(ni.nuclide)
        if prior <= 0.2:
            notes_lines.append(
                f"  Removed {ni.nuclide}: rare-isotope (prior={prior:.2f}) "
                f"failed proportionality ({prop.reason})"
            )
            new_rejected = NuclideIdentification(
                nuclide=ni.nuclide, detected=False,
                reason=(f"Rare-isotope identification failed intensity-ratio "
                        f"proportionality: {prop.reason}. "
                        f"Failed pairs: {prop.failed_pairs}. "
                        f"Lsrm §14.4 — without proportional multi-line "
                        f"evidence, the matches are likely chance overlaps "
                        f"with interfering nuclides."),
                characteristic_line_keV=ni.characteristic_line_keV,
                matched_lines=(),
            )
            rejected.append(new_rejected)
            detected = [d for d in detected if d.nuclide != ni.nuclide]
            if ni.nuclide in ni_by_name:
                del ni_by_name[ni.nuclide]
        else:
            notes_lines.append(
                f"  Warning {ni.nuclide}: failed proportionality "
                f"({prop.reason}); possible interference. Kept (common prior).")

    # --- Rule 5 / F-40: Secondary-feature anti-misidentification ---
    # When EVERY matched line of a candidate falls inside the empirical
    # secondary-feature position range (p10..p90) of an already-detected
    # parent in the v2 catalog, the candidate is most likely a
    # misattribution of the parent's Compton edge / backscatter / escape.
    # We exclude the parent's `photopeak` feature: true photopeak overlap
    # is handled by NAI_CONFUSION_MAP / Rule 3 (CI tiebreaker).
    if apply_secondary_feature_rule:
        cat_v2 = load_catalog_v2()
        catalog_nuclides = set(cat_v2.get("nuclides", {}).keys())
        detected_parents_in_cat = [
            ni.nuclide for ni in detected
            if ni.nuclide in catalog_nuclides
        ]

        # F-73 (v1.11.1): chain-membership exemption.
        # Nuclides of the SAME natural-decay chain as the parent must not
        # be demoted by the parent's secondaries. Pb-212 photopeak (238.6
        # keV, Th-232 chain) accidentally falls inside Tl-208 backscatter
        # window [234.7..235.9] on NaI 63×63 — but Pb-212 is a legitimate
        # Th-232 chain daughter that ALWAYS accompanies Tl-208 in natural
        # background, not a misattribution. Same logic for U-238 (Ra-226)
        # chain: Pb-214 / Bi-214 / Ra-226 must not exclude each other.
        from gamma.data.nuclide_library import get_nuclide as _gn

        def _chain_of(name: str) -> str:
            nuc = _gn(name)
            if nuc is None:
                return ""
            return str(nuc.get("chain", "") or "")

        if detected_parents_in_cat:
            for ni in list(detected):
                if not ni.matched_lines:
                    continue
                # Multi-line strong evidence overrides the rule.
                if len(ni.matched_lines) > secondary_max_lines:
                    continue
                # Parents available to explain this candidate's lines.
                # Excludes (a) the candidate itself and (b) parents that
                # share a natural-decay chain with the candidate.
                cand_chain = _chain_of(ni.nuclide)
                explaining_parents = []
                for p in detected_parents_in_cat:
                    if p == ni.nuclide:
                        continue
                    if cand_chain and cand_chain != "-" \
                            and _chain_of(p) == cand_chain:
                        # Same chain — both are legitimate ЕРН daughters,
                        # skip Rule 5 between them. (F-73)
                        continue
                    explaining_parents.append(p)
                if not explaining_parents:
                    continue

                explanations = []  # list of (parent, LineMatch, hit_dict)
                all_explained = True
                for m in ni.matched_lines:
                    line_hit = None
                    for parent in explaining_parents:
                        hits = matches_secondary(
                            parent, m.peak_E_keV, span="p10p90",
                        )
                        # Exclude the parent's own photopeak: true
                        # photopeak overlaps belong to other rules.
                        non_pp = [h for h in hits
                                  if h.get("feature") != "photopeak"]
                        if non_pp:
                            line_hit = (parent, m, non_pp[0])
                            break
                    if line_hit is None:
                        all_explained = False
                        break
                    explanations.append(line_hit)

                if all_explained and explanations:
                    descrs = []
                    for parent_nuc, mtch, hit in explanations:
                        lo, hi = hit["range"]
                        descrs.append(
                            f"{mtch.peak_E_keV:.1f} keV ↔ {parent_nuc} "
                            f"{hit['feature']} "
                            f"[{lo:.1f}..{hi:.1f}]"
                        )
                    notes_lines.append(
                        f"  Removed {ni.nuclide}: all "
                        f"{len(explanations)} matched line(s) explained as "
                        f"secondary features — {'; '.join(descrs)}"
                    )
                    new_rejected = NuclideIdentification(
                        nuclide=ni.nuclide, detected=False,
                        reason=(f"All matched lines fall inside observed "
                                f"secondary-feature ranges (p10..p90) of "
                                f"detected parent(s): {'; '.join(descrs)}. "
                                f"Positions consistent with Compton edge, "
                                f"backscatter, or escape continuum of the "
                                f"parent γ-emission, not a new nuclide line. "
                                f"Lsrm §14.4 + Knoll §10 — Gamma-1S NaI 63×63 "
                                f"secondary_peaks_v2 catalog (F-40)."),
                        characteristic_line_keV=ni.characteristic_line_keV,
                        matched_lines=(),
                    )
                    rejected.append(new_rejected)
                    detected = [
                        d for d in detected if d.nuclide != ni.nuclide
                    ]
                    if ni.nuclide in ni_by_name:
                        del ni_by_name[ni.nuclide]

    # --- Rule 4b: Ra-226 chain equilibrium analysis ---
    # Ra-chain nuclides have the radon-loss disequilibrium problem:
    # ratios WITHIN Group A (Ra-226 186) and WITHIN Group B (Rn-222
    # daughters Pb-214/Bi-214) are preserved, but the A/B ratio can be
    # arbitrary. So we check each group separately. Pb-210 has its own
    # independent source (shielding lead contamination) and is NOT used
    # to validate the chain.
    ra_chain_matches_by_nuc = {}
    for nuc in ("Ra-226", "Pb-214", "Bi-214", "Pb-210"):
        if nuc in ni_by_name:
            ra_chain_matches_by_nuc[nuc] = list(ni_by_name[nuc].matched_lines)

    if ra_chain_matches_by_nuc:
        ra_chain = check_ra226_chain_equilibrium(ra_chain_matches_by_nuc)
        if not ra_chain.chain_consistent:
            # Find which group failed and demote the relevant nuclides
            for group_name, group_result in ra_chain.group_results.items():
                if group_name == "Pb210_long_lived":
                    continue  # informational only
                if (group_result.n_lines_checked >= 2
                        and not group_result.passed):
                    notes_lines.append(
                        f"  Ra-chain group {group_name} failed: "
                        f"{group_result.reason}"
                    )
                    # Don't reject — chain disequilibrium is a real
                    # physical phenomenon (Rn loss). Just flag.
        else:
            if ra_chain.rn_retention_ratio is not None:
                notes_lines.append(
                    f"  Ra-chain consistent. Rn retention indicator "
                    f"(σ_186/I_186)/(σ_609/I_609) = "
                    f"{ra_chain.rn_retention_ratio:.3f} "
                    f"(equilibrium ≈ ε(186)/ε(609) on NaI ≈ 2)."
                )

    # --- Rule 4c: U-235 vs Ra-226 at 186 keV ROI ---
    # The 185.72 (U-235) and 186.21 (Ra-226) lines are unresolvable on
    # any scintillator (Δ = 0.5 keV << FWHM). Decision rule:
    #   • If Ra-chain (Pb-214 and/or Bi-214) is detected with ≥2 lines
    #     and consistent within Group B, the 186 ROI → Ra-226.
    #   • If U-235 shows BOTH 143.76 and 185.72 in proportion
    #     (universal Rule 4 above verified this), AND Ra-chain is NOT
    #     present, the 186 ROI → U-235.
    #   • If both proportionalities hold, it's a mixture — keep both.
    # For now we apply the simpler rule: Ra-chain presence wins.
    if "U-235" in ni_by_name:
        u235 = ni_by_name["U-235"]
        # Check Ra-chain presence (excluding Pb-210 which is unreliable)
        ra_chain_indicators_present = [
            n for n in ("Pb-214", "Bi-214")
            if n in ni_by_name and len(ni_by_name[n].matched_lines) >= 2
        ]
        if ra_chain_indicators_present:
            # Ra-chain wins the 186 ROI.
            notes_lines.append(
                f"  Removed U-235: Ra-chain present with multi-line evidence "
                f"({', '.join(ra_chain_indicators_present)} ≥2 lines each); "
                f"186 keV ROI attributed to Ra-226 186.21 keV."
            )
            new_rejected = NuclideIdentification(
                nuclide="U-235", detected=False,
                reason=(f"Ra-chain disambiguation: Pb-214/Bi-214 detected "
                        f"with proportional Group B → 186 keV ROI = "
                        f"Ra-226 186.21 keV, not U-235 185.72 keV. "
                        f"Lsrm §14.4 principle of least action."),
                characteristic_line_keV=u235.characteristic_line_keV,
                matched_lines=(),
            )
            rejected.append(new_rejected)
            detected = [d for d in detected if d.nuclide != "U-235"]
            del ni_by_name["U-235"]

    # --- Rule 6 / TD-3: Single-peak-collapse demotion ---
    # When EVERY matched line of a (rare-prior, multi-window-library)
    # nuclide hits the same single peak_channel, the "multi-line"
    # confirmation is illusory: it's just N library lines that all fall
    # inside one ID window on a low-resolution detector. The peak is
    # actually a single feature and should count as 1-line evidence,
    # not as multi-line confirmation. This guards against the
    # Cs-137 → U-235 false positive introduced when the IAEA refresh
    # expanded U-235 from 4 to 13 lines, so the 661.66 keV continuum +
    # 185.71-keV window now traps up to 5 nearby U-235 library lines
    # on a single peak. Same logic applies to any rare-prior nuclide
    # with dense library lines.
    REQUIRE_DISTINCT_PEAKS_NUCLIDES = {
        "U-235", "U-238", "Th-234", "Pa-234m", "Pb-210",
        "Cs-134",  # 4 close lines around 605 keV → same problem
        "Eu-152", "Eu-154", "Eu-155",  # dense L-α x-ray multiplets
        "Am-241",
    }
    for ni in list(detected):
        if ni.nuclide not in REQUIRE_DISTINCT_PEAKS_NUCLIDES:
            continue
        if len(ni.matched_lines) < 2:
            continue
        distinct_peaks = {m.peak_channel for m in ni.matched_lines}
        if len(distinct_peaks) >= 2:
            continue  # multi-peak — OK
        # All lines collapsed onto one peak — single-line evidence.
        # For the "dense library / window-trap" nuclide set, this is
        # treated as insufficient evidence per methodology regardless
        # of prior: post-IAEA-refresh these nuclides have ≥4 library
        # lines within a single NaI window, so multi-line confirmation
        # MUST come from distinct peaks, not from N library lines all
        # hitting the same feature.
        ch = next(iter(distinct_peaks))
        n_lines = len(ni.matched_lines)
        notes_lines.append(
            f"  Removed {ni.nuclide}: {n_lines} matched library lines "
            f"all collapsed onto a single peak (ch={ch}); single-peak "
            f"evidence is insufficient for a dense-library nuclide. "
            f"TD-3 multi-peak confirmation rule."
        )
        new_rejected = NuclideIdentification(
            nuclide=ni.nuclide, detected=False,
            reason=(f"Multi-peak confirmation required for dense-library "
                    f"nuclide: all {n_lines} matched library lines "
                    f"collapsed onto single peak ch={ch}. "
                    f"Post-IAEA-refresh library has dense neighbouring "
                    f"lines that all fall in one ID window on NaI; "
                    f"this is not multi-line evidence — it's a single "
                    f"feature claimed by N nearby library lines. "
                    f"Lsrm §14.4 principle of least action (TD-3)."),
            characteristic_line_keV=ni.characteristic_line_keV,
            matched_lines=(),
        )
        rejected.append(new_rejected)
        detected = [d for d in detected if d.nuclide != ni.nuclide]
        if ni.nuclide in ni_by_name:
            del ni_by_name[ni.nuclide]

    # --- Rule 7 / BUG-51: Nuisance-line suppression ---
    # When a peak is claimed by multiple nuclides, and ONE claimant uses a
    # very weak nuisance line (library_I_pct < 5% of that nuclide's max
    # library line intensity) while ANOTHER claimant uses a strong line
    # (library_I_pct > 50% of ITS nuclide's max), the weak claim is most
    # likely a chance overlap and should be suppressed.
    #
    # Motivating case (A4 §4 diagnostic 2026-06-04):
    #   • 508.38 keV peak (511 annihilation, σ=302) was being attributed to
    #     Eu-152 503.467 (I=0.1524% → 0.53% of Eu-152 max=28.53% at 121.78).
    #   • 656.89 keV peak (Cs-137 661, σ=170) was being attributed to
    #     Eu-152 656.489 (I=0.1441% → 0.50% of Eu-152 max).
    # These nuisance Eu-152 lines diluted the Eu-152 multi-line aggregation
    # (BUG-39 weighted-mean residual −59.3% on AmTiCsEu fixture).
    #
    # Guard 1: skip the cut if the candidate nuclide has ALL library lines
    # < NUISANCE_THRESHOLD_PCT — for some rare isotopes ALL lines are weak.
    # We must not strip every claim from such a nuclide.
    # Guard 2: only apply when the OTHER claimant on the same peak has a
    # genuinely strong library line (I_other ≥ NUISANCE_STRONG_OTHER_PCT
    # of its own nuclide's max).
    # Guard 3: never strip a claim on the candidate's CHARACTERISTIC line —
    # that line was selected as the minimum-MDA detection criterion.
    NUISANCE_THRESHOLD_FRAC = 0.05   # weak if < 5% of nuclide's max line I
    NUISANCE_STRONG_OTHER_FRAC = 0.50  # rival is strong if ≥ 50% of its max

    # Build a cache of (nuclide → max library line intensity) for every
    # detected nuclide. Pull from data/nuclides.json via get_nuclide.
    from gamma.data.nuclide_library import get_nuclide as _gn_b51

    def _max_library_I(name: str) -> float:
        rec = _gn_b51(name)
        if rec is None:
            return 0.0
        lines = rec.get("lines", []) or []
        if not lines:
            return 0.0
        # line format: [E, I, sigma_I]
        return max((float(l[1]) for l in lines if len(l) >= 2), default=0.0)

    # Chain-membership cache for Guard 4 (TD-3 / F-73 precedent — same as
    # Rules 3 and 5). When the candidate and the strong-rival both belong
    # to the SAME natural decay chain (Th-232 / Ra-226 / U-238 / U-235),
    # both are legitimate daughters that always accompany each other in
    # secular equilibrium. Stripping a weak Tl-208 claim because Pb-212
    # has a strong rival line collapses Tl-208 activity to zero on Th-232
    # fixtures (regression seen 2026-06-05 on test_bi212_tcs,
    # test_multiplet_nnls_ac228, test_td2_cache_isolation — Tl-208 Bq
    # collapses to 1e-24 because all non-characteristic lines hit chain
    # rivals). Skip the cut for same-chain co-claimants.
    def _chain_of_b51(name: str) -> str:
        rec = _gn_b51(name)
        if rec is None:
            return ""
        return str(rec.get("chain", "") or "")

    max_I_by_nuclide = {n: _max_library_I(n) for n in ni_by_name}
    chain_by_nuclide = {n: _chain_of_b51(n) for n in ni_by_name}

    # Rebuild peak_users (Rules 1/2/4/5/6 may have removed nuclides).
    peak_users = {}
    for ni in detected:
        for m in ni.matched_lines:
            peak_users.setdefault(m.peak_channel, []).append((ni.nuclide, m))

    b51_to_clear: dict = {}  # nuclide → set of peak_channels to remove
    b51_notes = []
    for ch, users in peak_users.items():
        if len(users) <= 1:
            continue
        # For each candidate user on this peak, check if it is a nuisance.
        # (a) library_I_pct < NUISANCE_THRESHOLD_FRAC × nuclide_max_I.
        # (b) some OTHER user on this peak has library_I_pct ≥
        #     NUISANCE_STRONG_OTHER_FRAC × its nuclide_max_I.
        for name, match in users:
            # Guard 3 — never strip the characteristic line claim
            if getattr(match, "is_characteristic", False):
                continue
            nuc_max = max_I_by_nuclide.get(name, 0.0)
            if nuc_max <= 0.0:
                continue
            cand_I = float(getattr(match, "library_I_pct", 0.0) or 0.0)
            cand_frac = cand_I / nuc_max
            if cand_frac >= NUISANCE_THRESHOLD_FRAC:
                continue
            # Guard 1 — nuclide must have at least one strong (≥5%) line
            # in its library overall; if not, do not strip (rare-emitter
            # case where every line is weak).
            #
            # Equivalent here: the nuclide max IS its strongest line; if
            # nuc_max itself < 5 percent, the rule does not fire — every
            # line of this nuclide is weak by nature (rare-emitter case,
            # e.g. Pu/Am ic_xray composites). We treat the absolute
            # threshold as 5% — i.e. require the nuclide to have at
            # least one library line ≥ 5% absolute intensity.
            if nuc_max < 5.0:
                continue
            # Check if some other claimant on this peak is genuinely strong.
            # Guard 4 — same-chain exemption (TD-3 / F-73 precedent):
            # if the strong rival and the candidate share the same natural
            # decay chain, do NOT strip the candidate. Both are legitimate
            # secular-equilibrium daughters.
            cand_chain = chain_by_nuclide.get(name, "")
            other_strong = False
            other_name = ""
            other_I = 0.0
            for other, other_match in users:
                if other == name:
                    continue
                other_max = max_I_by_nuclide.get(other, 0.0)
                if other_max <= 0.0:
                    continue
                I_o = float(getattr(other_match, "library_I_pct", 0.0) or 0.0)
                if I_o / other_max >= NUISANCE_STRONG_OTHER_FRAC:
                    other_chain = chain_by_nuclide.get(other, "")
                    # Skip same-chain rivals — both are legitimate daughters.
                    if (cand_chain and other_chain
                            and cand_chain != "-" and other_chain != "-"
                            and cand_chain == other_chain):
                        continue
                    other_strong = True
                    other_name = other
                    other_I = I_o
                    break
            if not other_strong:
                continue
            b51_to_clear.setdefault(name, set()).add(ch)
            b51_notes.append(
                f"  Removed {name} nuisance-line claim on ch={ch} "
                f"(library_I={cand_I:.3f}% = {cand_frac*100:.1f}% of "
                f"nuclide max {nuc_max:.2f}%): peak better explained by "
                f"{other_name} strong line (I={other_I:.2f}%). BUG-51."
            )

    if b51_to_clear:
        new_detected = []
        for ni in detected:
            clears = b51_to_clear.get(ni.nuclide, set())
            if not clears:
                new_detected.append(ni)
                continue
            remaining = tuple(
                m for m in ni.matched_lines if m.peak_channel not in clears
            )
            if not remaining:
                # Should not happen given Guard 3 (characteristic preserved),
                # but defensively reject if it does.
                new_rejected = NuclideIdentification(
                    nuclide=ni.nuclide, detected=False,
                    reason=(f"All matched lines were nuisance-line claims on "
                            f"peaks better explained by other nuclides "
                            f"(BUG-51 cut)."),
                    characteristic_line_keV=ni.characteristic_line_keV,
                    matched_lines=(),
                )
                rejected.append(new_rejected)
                continue
            new_detected.append(replace(ni, matched_lines=remaining))
        detected = new_detected
        ni_by_name = {ni.nuclide: ni for ni in detected}
        notes_lines.extend(b51_notes)

    # --- Rule 3: shared peaks with CI tiebreaker ---
    # Rebuild peak map after removals
    peak_users = {}
    for ni in detected:
        for m in ni.matched_lines:
            peak_users.setdefault(m.peak_channel, []).append((ni.nuclide, m))

    # For peaks claimed by multiple nuclides, the lower-CI nuclide
    # loses its claim if the higher-CI nuclide exceeds it by threshold.
    # This applies regardless of how many lines the loser has — if all
    # of a nuclide's claimed peaks are explained by other higher-CI
    # nuclides, it should not survive (Lsrm §14.4 principle of least
    # action).
    #
    # TD-3 (v1.18.27.1): two exemptions added to prevent false negatives
    # on Pb-212 (238.6 keV) and K-40 (1460.8 keV) post-IAEA-refresh:
    #   (a) Natural-chain co-membership: if both claimants share the
    #       same natural decay chain (e.g. Tl-208 + Pb-212 both Th-232),
    #       neither should cancel the other on a shared peak — both are
    #       legitimate chain daughters that always accompany each other
    #       (analogous to F-73 exemption in Rule 5).
    #   (b) Dominant-intensity reversal: if the loser's claimed line on
    #       this peak has library intensity ≥ INTENSITY_DOMINANCE_FACTOR
    #       times higher than the winner's claimed line, AND the loser's
    #       line is the loser's characteristic line, the peak is in fact
    #       better explained by the loser's line (high CI of winner came
    #       from its OTHER lines, not from this peak). Keep both claims.
    from gamma.data.nuclide_library import get_nuclide as _gn_r3

    def _chain_of_r3(name: str) -> str:
        nuc = _gn_r3(name)
        if nuc is None:
            return ""
        return str(nuc.get("chain", "") or "")

    INTENSITY_DOMINANCE_FACTOR = 10.0  # loser's I ≥ 10× winner's I → keep

    nuclides_to_clear: dict = {}  # nuclide → set of peak_channels to remove
    for ch, users in peak_users.items():
        if len(users) <= 1:
            continue
        # Multiple claims on same channel
        users_with_ci = []
        for name, match in users:
            ni = ni_by_name.get(name)
            ci = ni.confidence.CI if (ni and ni.confidence) else 0
            n_lines = len(ni.matched_lines) if ni else 0
            users_with_ci.append((ci, n_lines, name, match))
        # Sort: higher CI first; on tie, more lines first
        users_with_ci.sort(key=lambda x: (-x[0], -x[1]))
        winner = users_with_ci[0]
        winner_ci = winner[0]
        winner_name = winner[2]
        winner_match = winner[3]
        winner_chain = _chain_of_r3(winner_name)
        for loser_ci, loser_n_lines, loser_name, loser_match in users_with_ci[1:]:
            if winner_ci - loser_ci < min_ci_difference_to_resolve:
                continue
            # (a) TD-3 chain-membership exemption — skip if same natural chain
            loser_chain = _chain_of_r3(loser_name)
            if (winner_chain and loser_chain
                    and winner_chain != "-" and loser_chain != "-"
                    and winner_chain == loser_chain):
                notes_lines.append(
                    f"  Kept {loser_name} claim on ch={ch}: same natural "
                    f"chain as {winner_name} ({winner_chain}); both are "
                    f"legitimate chain daughters (TD-3)."
                )
                continue
            # (b) TD-3 intensity-dominance reversal — the peak is more
            # naturally attributed to the high-intensity loser line.
            w_I = float(getattr(winner_match, "library_I_pct", 0.0) or 0.0)
            l_I = float(getattr(loser_match, "library_I_pct", 0.0) or 0.0)
            loser_char = bool(getattr(loser_match, "is_characteristic", False))
            if (loser_char and w_I > 0.0
                    and l_I >= INTENSITY_DOMINANCE_FACTOR * w_I):
                notes_lines.append(
                    f"  Kept {loser_name} claim on ch={ch}: characteristic "
                    f"line (I={l_I:.2f}%) dominates {winner_name}'s match "
                    f"on this peak (I={w_I:.2f}%); peak better explained "
                    f"by {loser_name} (TD-3 intensity dominance)."
                )
                continue
            nuclides_to_clear.setdefault(loser_name, set()).add(ch)
            notes_lines.append(
                f"  Removed {loser_name} claim on ch={ch}: "
                f"superseded by {winner_name} (CI {winner_ci:.2f} > {loser_ci:.2f})"
            )

    # Apply the clear: remove matches from nuclides; if a nuclide is
    # left with zero matched lines, move it to rejected.
    new_detected = []
    for ni in detected:
        clears = nuclides_to_clear.get(ni.nuclide, set())
        if not clears:
            new_detected.append(ni)
            continue
        remaining = tuple(m for m in ni.matched_lines if m.peak_channel not in clears)
        if not remaining:
            new_rejected = NuclideIdentification(
                nuclide=ni.nuclide, detected=False,
                reason=(f"All matched lines removed by disambiguation: "
                        f"every claimed peak is better explained by another nuclide. "
                        f"Lsrm §14.4 principle of least action."),
                characteristic_line_keV=ni.characteristic_line_keV,
                matched_lines=(),
            )
            rejected.append(new_rejected)
            continue
        # Re-flag characteristic if needed
        has_char = any(m.is_characteristic for m in remaining)
        if not has_char:
            # Promote highest-σ remaining line as new "characteristic"
            best = max(remaining, key=lambda m: (m.significance_currie or 0.0))
            promoted = LineMatch(
                nuclide=best.nuclide,
                library_E_keV=best.library_E_keV,
                library_I_pct=best.library_I_pct,
                peak_channel=best.peak_channel,
                peak_E_keV=best.peak_E_keV,
                peak_sigma=best.peak_sigma,
                residual_keV=best.residual_keV,
                is_characteristic=True,
                peak_area=best.peak_area,
                peak_area_uncertainty=best.peak_area_uncertainty,
                peak_area_source=best.peak_area_source,
                # BUG-34 Phase 1+2: propagate explicit successors
                significance_currie=best.significance_currie,
                gauss_sigma_keV=best.gauss_sigma_keV,
            )
            remaining = tuple(
                promoted if m.peak_channel == best.peak_channel else m
                for m in remaining
            )
        new_detected.append(replace(ni, matched_lines=remaining))

    detected = new_detected

    # Refresh notes
    new_notes = result.notes
    if notes_lines:
        new_notes = (new_notes + "\nDisambiguation applied:\n" +
                     "\n".join(notes_lines)).strip()

    return IdentificationResult(
        detector_type=result.detector_type,
        window=result.window,
        candidates_considered=result.candidates_considered,
        detected_nuclides=tuple(detected),
        rejected_nuclides=tuple(rejected),
        unmatched_peaks=result.unmatched_peaks,
        notes=new_notes,
    )


__all__ = [
    "NATURAL_CHAIN_NUCLIDES",
    "POSITRON_EMITTERS_NEAR_511",
    "NAI_CONFUSION_MAP",
    "NAI_CHAR_OVERLAP_PAIRS",
    "MIN_CHAIN_LINES_CHAR_OVERLAP",
    "disambiguate_identifications",
]
