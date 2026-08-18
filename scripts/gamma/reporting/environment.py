"""
Measurement-environment classifier (Step 2 of SKILL.md, surfaced by the
report generators).

Two canonical environments per SKILL.md:

* **low_background** (shielded castle): Pb K-XRF at 73–87 keV present,
  ЕРН lines (⁴⁰K 1460, ²¹⁴Bi series, ²⁰⁸Tl 2614) strongly suppressed.
* **natural** (open lab): no Pb K-XRF, strong ЕРН lines from concrete
  walls.

Diagnostic signals (per SKILL.md Step 2):

1. Presence/absence of Pb K-XRF at 73–87 keV in residual_classifications
   (label = "xrf" with element = "Pb").
2. Filename hints (`bkg`, `lab`, `shield`, `castle`).
3. Strength of ЕРН lines (K-40 1460, Tl-208 2614, Bi-214 1764) — these
   surface either as detected nuclides (post Stage 1) or as anchor
   matches.

The classifier returns one of three labels:

* `"low_background"` — Pb K-XRF is present (shield fluorescence) AND
  ЕРН strength is consistent with shielding suppression.
* `"natural"` — no Pb K-XRF AND at least one ЕРН anchor at significant
  rate.
* `"unknown"` — neither signature is decisive (e.g. a pure-source
  measurement on a clean isolator, or sparse data).
"""
from __future__ import annotations

from typing import Optional


ENV_NATURAL = "natural"
ENV_LOW_BG = "low_background"
ENV_UNKNOWN = "unknown"
ENV_BACKGROUND_ONLY = "background_only"


# Filename tokens hinting at the environment. These are consulted as a
# tie-breaker when the spectral evidence is ambiguous.
_LOW_BG_HINT_TOKENS = ("shield", "castle", "защит", "защ_", "закр_кр", "пасс")
_NATURAL_HINT_TOKENS = ("открыт", "open_lid", "openlid", "лаба", "lab_open")


def _is_background_only(staged_result) -> bool:
    """True if the analysed file is a pure background spectrum (D-01).

    Detected by either:
      * ``result.is_background`` flag = True;
      * ``sample_type_hint`` / ``sample_type_canonical`` in
        {"background", "bg"};
      * filename basename starts with ``bg_`` / ``bkg_``.
    """
    if getattr(staged_result, "is_background", False):
        return True
    hint = (getattr(staged_result, "sample_type_hint", "") or "").lower()
    canon = (getattr(staged_result, "sample_type_canonical", "") or "").lower()
    if hint in ("background", "bg") or canon in ("background", "bg"):
        return True
    spec = getattr(staged_result, "spec", None)
    if spec is not None:
        sp = getattr(spec, "source_path", "") or ""
        if "\\" in sp:
            leaf = sp.rsplit("\\", 1)[-1]
        elif "/" in sp:
            leaf = sp.rsplit("/", 1)[-1]
        else:
            leaf = sp
        low = leaf.lower()
        if low.startswith("bg_") or low.startswith("bkg_"):
            return True
    return False


def _filename_hint(filename: str) -> Optional[str]:
    """Inspect filename for environment tokens; return ENV_* or None."""
    if not filename:
        return None
    low = filename.lower()
    for tok in _LOW_BG_HINT_TOKENS:
        if tok in low:
            return ENV_LOW_BG
    for tok in _NATURAL_HINT_TOKENS:
        if tok in low:
            return ENV_NATURAL
    return None


def _has_pb_kxrf(residual_classifications) -> bool:
    """True if at least one residual was classified as Pb K-XRF.

    The Pb fluorescence triplet is the canonical low-background castle
    signature — concrete walls give no shield-XRF fluorescence, only
    direct ЕРН emission.
    """
    if not residual_classifications:
        return False
    for rc in residual_classifications:
        if getattr(rc, "label", "") == "xrf" and getattr(rc, "element", "") == "Pb":
            return True
    return False


def _ern_strength_score(staged_result) -> float:
    """Return a 0..1 score of how strongly ЕРН lines are present.

    1.0 → several ЕРН anchors confirmed at significant σ (typical of
    natural-background concrete-wall labs).
    0.0 → no ЕРН anchors at all (typical of pure-source measurements
    or strongly shielded backgrounds).

    Uses `anchor_matches` and `pattern_confirmations` because they are
    populated for both background spectra (mode = background_7line) and
    sample spectra (mode = sample_anchor_rank).
    """
    score = 0.0
    ern_nuclides = {"K-40", "Tl-208", "Bi-214", "Ac-228", "Pb-212", "Pb-214"}

    am = list(getattr(staged_result, "anchor_matches", []) or [])
    n_confirmed_anchors = 0
    for m in am:
        anchor = getattr(m, "anchor", None)
        if anchor is None:
            continue
        nuc = getattr(anchor, "nuclide", "")
        if nuc in ern_nuclides and not getattr(m, "partner_required_but_missing", True):
            n_confirmed_anchors += 1
    # Each confirmed ЕРН anchor contributes 0.25, capped at 1.0
    score += min(1.0, 0.25 * n_confirmed_anchors)

    # Detected nuclides also count
    detected = list(getattr(staged_result, "final_detected", []) or [])
    n_ern_detected = sum(1 for n in detected
                         if getattr(n, "nuclide", "") in ern_nuclides)
    score = min(1.0, score + 0.20 * n_ern_detected)
    return score


def classify_environment(staged_result) -> str:
    """Return one of ENV_NATURAL / ENV_LOW_BG / ENV_UNKNOWN.

    Rule table (first match wins):

    | Pb K-XRF | ЕРН score | Filename hint | → Result          |
    |----------|-----------|---------------|--------------------|
    | yes      | any       | any           | low_background    |
    | no       | ≥ 0.5     | any           | natural           |
    | no       | < 0.5     | low_bg token  | low_background    |
    | no       | < 0.5     | natural token | natural           |
    | no       | < 0.5     | none          | unknown           |

    The first rule (Pb K-XRF wins outright) reflects that shield
    fluorescence is the strongest single-signal indicator — there is
    no natural source of 74–87 keV Pb XRF in a non-shielded environment.
    """
    # D-01 — when the input is a pure background spectrum, environment
    # is "background_only" and downstream narrative must NOT speak
    # about a sample.
    if _is_background_only(staged_result):
        return ENV_BACKGROUND_ONLY

    residuals = getattr(staged_result, "residual_classifications", []) or []
    if _has_pb_kxrf(residuals):
        return ENV_LOW_BG

    ern = _ern_strength_score(staged_result)
    spec = getattr(staged_result, "spec", None)
    filename = ""
    if spec is not None:
        sp = getattr(spec, "source_path", "") or ""
        # Take the leaf only (cross-platform)
        if "\\" in sp:
            filename = sp.rsplit("\\", 1)[-1]
        elif "/" in sp:
            filename = sp.rsplit("/", 1)[-1]
        else:
            filename = sp

    if ern >= 0.5:
        return ENV_NATURAL

    hint = _filename_hint(filename)
    if hint:
        return hint

    return ENV_UNKNOWN


# G3 / v1.31.2 -- continuum-level diagnostic (annotation, no decision impact).
#
# Spec: AUDIT_F-419_skill_vs_canon.md / G3 -- expose three integral metrics
# next to the existing categorical environment label so an operator can
# sanity-check classify_environment() at a glance:
#
#   * total_cps                 = total counts across spectrum / live_time
#   * bg_line_dominance_pct     = 100 * (sum of identified-ERN line areas
#                                  / total_counts), proxy of how much of
#                                  the total flux is naturally-occurring
#                                  background activity
#   * environment_hint_by_cps   = "low_bg" (<1 cps), "intermediate" (1-10),
#                                  "natural" (10-100), "high" (>100)
#
# Thresholds locked by audit. Hint is independent of classify_environment
# (which uses Pb K-XRF + ERN-anchor logic) -- divergence between them is the
# diagnostic signal the operator wants to surface.
ERN_NUCLIDES_FOR_DOMINANCE = frozenset({
    "K-40", "Tl-208", "Bi-214", "Ac-228", "Pb-212", "Pb-214",
    "Bi-212", "Ra-226", "Th-232", "U-238",
})

CPS_HINT_LOW_BG_MAX = 1.0
CPS_HINT_INTERMEDIATE_MAX = 10.0
CPS_HINT_NATURAL_MAX = 100.0


def _spec_total_counts_and_live_time(spec):
    """Return (total_counts, live_time) from a Spectrum; (None, None) on failure."""
    if spec is None:
        return None, None
    try:
        counts = getattr(spec, "counts", None)
        if counts is None:
            return None, None
        total = float(sum(counts))
    except (TypeError, ValueError):
        return None, None
    try:
        t_live = float(getattr(spec, "live_time", 0.0) or 0.0)
    except (TypeError, ValueError):
        return total, None
    if t_live <= 0:
        return total, None
    return total, t_live


def _ern_line_area_sum(staged_result) -> float:
    """Sum of peak_area across matched lines of ERN nuclides."""
    total = 0.0
    detected = list(getattr(staged_result, "final_detected", []) or [])
    for ni in detected:
        nuc = str(getattr(ni, "nuclide", "") or "")
        if nuc not in ERN_NUCLIDES_FOR_DOMINANCE:
            continue
        for m in getattr(ni, "matched_lines", ()) or ():
            area = getattr(m, "peak_area", None)
            try:
                a = float(area) if area is not None else 0.0
            except (TypeError, ValueError):
                a = 0.0
            if a > 0:
                total += a
    return total


def _cps_hint(total_cps: float) -> str:
    if total_cps < CPS_HINT_LOW_BG_MAX:
        return "low_bg"
    if total_cps < CPS_HINT_INTERMEDIATE_MAX:
        return "intermediate"
    if total_cps < CPS_HINT_NATURAL_MAX:
        return "natural"
    return "high"


def continuum_diagnostic(staged_result) -> dict:
    """Compute G3 continuum diagnostic block.

    Returns a dict with the three metrics. Missing inputs -> None values;
    callers serialize as-is. Never raises -- always returns a dict.
    """
    spec = getattr(staged_result, "spec", None)
    total_counts, t_live = _spec_total_counts_and_live_time(spec)
    if total_counts is None or t_live is None or t_live <= 0:
        return {
            "total_cps": None,
            "bg_line_dominance_pct": None,
            "environment_hint_by_cps": None,
            "ern_line_area_sum": None,
            "total_counts": total_counts,
            "live_time_s": t_live,
        }
    total_cps = total_counts / t_live
    ern_area = _ern_line_area_sum(staged_result)
    dominance_pct = (100.0 * ern_area / total_counts) if total_counts > 0 else 0.0
    return {
        "total_cps": total_cps,
        "bg_line_dominance_pct": dominance_pct,
        "environment_hint_by_cps": _cps_hint(total_cps),
        "ern_line_area_sum": ern_area,
        "total_counts": total_counts,
        "live_time_s": t_live,
    }

__all__ = [
    "classify_environment",
    "ENV_NATURAL", "ENV_LOW_BG", "ENV_UNKNOWN", "ENV_BACKGROUND_ONLY",
    "continuum_diagnostic", "ERN_NUCLIDES_FOR_DOMINANCE",
    "CPS_HINT_LOW_BG_MAX", "CPS_HINT_INTERMEDIATE_MAX", "CPS_HINT_NATURAL_MAX",
]