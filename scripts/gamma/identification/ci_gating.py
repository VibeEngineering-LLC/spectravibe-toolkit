"""
Confidence-index gating of identified nuclides (F-60 / v1.11.1).

After identify_nuclides + disambiguate, every NuclideIdentification has
a `confidence.CI` value. Per Lsrm Algorithmic Foundations §14.3, CI is

  CI = log10[ 1 / (δE₁ · δE₂ · ... · δI₂ · δI₃ · ...) ]

where δE_i = ΔE_i/E_i are relative energy uncertainties of confirmed
lines and δI_i are relative intensity-ratio uncertainties. Higher CI
means more reliable identification.

Reference CI values per Lsrm methodology:
  Cs-137 (single line)       : ~2
  Co-60 (well-resolved pair) : ~5–7
  Eu-152 (multi-line)        : ~20–67

This module classifies each identification into three tiers based on
CI thresholds:

  confirmed  (CI ≥ 2.0)   — reliable; safe to report as primary result
  tentative  (1.0 ≤ CI<2) — present-but-weak; flag for operator review
  noise      (CI < 1.0)   — likely spurious; demote unless externally
                            corroborated (e.g. matches an anchor-rank
                            anchor or an express pattern)

Default thresholds can be overridden via `gate_identifications` kwargs.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence


# ──────────────────────────────────────────────────────────────────
# Default thresholds — per Lsrm methodology reference values
# ──────────────────────────────────────────────────────────────────

CI_CONFIRMED_THRESHOLD = 2.0   # >= → confirmed
CI_TENTATIVE_THRESHOLD = 1.0   # >= → tentative; below → noise


# ──────────────────────────────────────────────────────────────────
# Result
# ──────────────────────────────────────────────────────────────────

@dataclass
class CIGating:
    confirmed: List               # NuclideIdentification with CI ≥ confirmed
    tentative: List               # NuclideIdentification with tentative ≤ CI < confirmed
    noise: List                   # CI < tentative
    promoted_by_anchor: List = field(default_factory=list)
    promoted_by_pattern: List = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    @property
    def all_classified(self) -> List:
        return self.confirmed + self.tentative + self.noise


# ──────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────

def gate_identifications(
    detected_nuclides: List,
    *,
    confirmed_threshold: float = CI_CONFIRMED_THRESHOLD,
    tentative_threshold: float = CI_TENTATIVE_THRESHOLD,
    anchor_confirmed_nuclides: Optional[Sequence[str]] = None,
    pattern_confirmed_nuclides: Optional[Sequence[str]] = None,
) -> CIGating:
    """
    Partition the detected_nuclides list into confidence tiers.

    A nuclide is **promoted from `tentative` to `confirmed`** if it is
    independently corroborated by either:
      • an anchor-rank match (Pass A) with rank ≤ 7, or
      • an express pattern confirmation (Pass B)

    This handles the K-40 single-line case (intrinsic CI ≈ 1.5–1.6) —
    K-40 1461 keV is one of the most reliable identifications on NaI
    despite its low single-line CI, so anchor confirmation promotes it.

    A nuclide is **demoted from `noise` to `tentative`** if any anchor
    or pattern corroborates it — single-line low-CI detections that
    nonetheless match high-rank anchors deserve a second look.

    Args:
        detected_nuclides: list of NuclideIdentification (after disambiguate)
        confirmed_threshold: CI value ≥ this → confirmed (default 2.0)
        tentative_threshold: CI value ≥ this → tentative (default 1.0)
        anchor_confirmed_nuclides: nuclide names that have an anchor-rank
            match (Pass A) — used for cross-promotion. Optional.
        pattern_confirmed_nuclides: nuclide names confirmed by an express
            pattern (Pass B). Optional.

    Returns:
        CIGating with three lists (confirmed/tentative/noise) plus
        provenance tracking of promoted entries.
    """
    anchors = set(anchor_confirmed_nuclides or ())
    patterns = set(pattern_confirmed_nuclides or ())

    confirmed: List = []
    tentative: List = []
    noise: List = []
    promoted_by_anchor: List = []
    promoted_by_pattern: List = []

    for nid in detected_nuclides:
        ci = getattr(getattr(nid, "confidence", None), "CI", None)
        if ci is None:
            # No CI computed — treat as tentative
            tier = "tentative"
        elif ci >= confirmed_threshold:
            tier = "confirmed"
        elif ci >= tentative_threshold:
            tier = "tentative"
        else:
            tier = "noise"

        # Cross-promotion via independent corroboration
        if tier == "tentative" and nid.nuclide in anchors:
            tier = "confirmed"
            promoted_by_anchor.append(nid)
        elif tier == "tentative" and nid.nuclide in patterns:
            tier = "confirmed"
            promoted_by_pattern.append(nid)
        elif tier == "noise" and (nid.nuclide in anchors or nid.nuclide in patterns):
            tier = "tentative"
            # Track promotion regardless of source
            if nid.nuclide in anchors:
                promoted_by_anchor.append(nid)
            else:
                promoted_by_pattern.append(nid)

        if tier == "confirmed":
            confirmed.append(nid)
        elif tier == "tentative":
            tentative.append(nid)
        else:
            noise.append(nid)

    notes: List[str] = []
    if promoted_by_anchor:
        names = ", ".join(n.nuclide for n in promoted_by_anchor)
        notes.append(f"Повышены до confirmed по anchor-rank: {names}")
    if promoted_by_pattern:
        names = ", ".join(n.nuclide for n in promoted_by_pattern)
        notes.append(f"Повышены до confirmed по express pattern: {names}")
    if noise:
        names = ", ".join(n.nuclide for n in noise)
        notes.append(f"Низкий CI (< {tentative_threshold}) → noise: {names}")

    return CIGating(
        confirmed=confirmed,
        tentative=tentative,
        noise=noise,
        promoted_by_anchor=promoted_by_anchor,
        promoted_by_pattern=promoted_by_pattern,
        notes=notes,
    )


__all__ = [
    "CIGating",
    "gate_identifications",
    "CI_CONFIRMED_THRESHOLD",
    "CI_TENTATIVE_THRESHOLD",
]
