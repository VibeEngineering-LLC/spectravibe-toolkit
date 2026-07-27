"""
In-chat summary — 3–8 lines per references/06_report_format.md.

Template:

  Spectrum: <filename> (<live_time_s> s live, dead-time <dead_time_pct>%)
  Detector: <detector_type>, environment: <environment>
  Identified nuclides: <list with CI in parens, sorted by CI desc>
  Completeness DC: <dc_pct>%, <flag>
  Shielding: <inferred composition>
  MDA highlights: <Cs-137 / Co-60 / Ra-226 / Th-232 in Bq or cps>
  Warnings: <if any>
  Report: <path/to/report.json>

The summary is generated from the JSON report dict, so the chat layer
does not depend on the StagedAnalysisResult internals.
"""
from __future__ import annotations

from typing import List, Optional


_MDA_HIGHLIGHT_NUCLIDES = ["Cs-137", "Co-60", "K-40", "Bi-214", "Tl-208", "Ac-228"]


def _fmt_dead_time(d) -> str:
    if d is None:
        return "?"
    return f"{d:.1f}%"


def _fmt_live(t) -> str:
    if t is None or t <= 0:
        return "?"
    if t >= 3600:
        return f"{t / 3600.0:.2f} h"
    if t >= 60:
        return f"{t / 60.0:.1f} min"
    return f"{t:.0f} s"


def _identified_line(json_dict) -> str:
    """One line listing identified nuclides with CI in parens, CI-desc sorted."""
    nucs = list(json_dict.get("identified_nuclides", []) or [])
    if not nucs:
        return "Identified: none in default Stage-1 candidate list."
    # Sort by CI desc; missing CI goes last
    nucs.sort(
        key=lambda n: (-(n.get("confidence_index") or 0.0), n.get("nuclide", "")),
    )
    parts = []
    for n in nucs[:8]:
        nuc = n.get("nuclide", "?")
        ci = n.get("confidence_index")
        if ci is not None:
            parts.append(f"{nuc} (CI={ci:.1f})")
        else:
            parts.append(nuc)
    extra = ""
    if len(nucs) > 8:
        extra = f" (+{len(nucs) - 8} more)"
    return "Identified: " + ", ".join(parts) + extra


def _completeness_line(json_dict) -> Optional[str]:
    cmp = json_dict.get("completeness", {}) or {}
    dc = cmp.get("dc_pct")
    flag = cmp.get("flag") or "unknown"
    if dc is None:
        return None
    return f"Completeness: DC = {dc:.1f}% [{flag}]."


def _shielding_line(json_dict) -> Optional[str]:
    """Infer shielding composition from the elemental XRF table."""
    xrf = json_dict.get("elemental_xrf", []) or []
    if not xrf:
        return None
    elements = [e.get("element", "?") for e in xrf]
    return "Shielding signature: " + ", ".join(elements) + " K-XRF observed."


def _mda_highlights_line(json_dict) -> Optional[str]:
    mda_rows = json_dict.get("mda", []) or []
    if not mda_rows:
        return None
    # Group by nuclide; pick the characteristic-line MDA (lowest E for that nuc)
    by_nuc = {}
    for row in mda_rows:
        nuc = row.get("nuclide", "")
        mda = row.get("MDA_Bq")
        if mda is None:
            continue
        prev = by_nuc.get(nuc)
        if prev is None or mda < prev[0]:
            by_nuc[nuc] = (mda, row.get("line_E_keV"))
    parts = []
    for nuc in _MDA_HIGHLIGHT_NUCLIDES:
        if nuc in by_nuc:
            mda_val, E = by_nuc[nuc]
            parts.append(f"{nuc}@{E:.0f}={mda_val:.2g} Bq")
    if not parts:
        return None
    return "MDA: " + "; ".join(parts) + "."


def _warnings_line(json_dict) -> Optional[str]:
    w = json_dict.get("warnings", []) or []
    if not w:
        return None
    # P0-8: warnings may contain dicts (structured alarm entries) or plain strings.
    def _w_str(entry) -> str:
        if isinstance(entry, dict):
            return entry.get("message") or entry.get("code") or str(entry)
        return str(entry)
    return "Warnings: " + " | ".join(_w_str(x) for x in w[:3])


def _chain_dominance_line(json_dict) -> Optional[str]:
    """F-88/F-89d — surface chain dominance + K-40/Ac-228 overlap warning +
    chain suppression."""
    diag = json_dict.get("diagnostics", {}) or {}
    cd = diag.get("chain_dominance") or {}
    k40_warn = bool(diag.get("k40_ac228_overlap_warning"))
    bits = []
    if cd.get("th232_dominant"):
        bits.append("Th-232 chain DOMINANT")
    if cd.get("u238_dominant"):
        bits.append("U-238 chain DOMINANT")
    if k40_warn:
        bits.append("⚠ K-40/Ac-228 overlap")
    # F-89d
    sup = cd.get("suppressed_chains") or []
    if sup:
        bits.append(f"suppressed: {'+'.join(sup)}")
    if not bits:
        return None
    return "Express: " + ", ".join(bits) + "."


def _bg_filename_line(json_dict) -> Optional[str]:
    """F-89a/b — one line covering filename hint + background status."""
    h = json_dict.get("header", {}) or {}
    bits = []
    hints = h.get("filename_isotope_hints") or []
    if hints:
        bits.append(f"filename → {', '.join(hints)}")
    bg = h.get("background_status", "")
    if bg == "subtracted_from_external_file":
        bits.append("bg SUBTRACTED")
    elif bg == "auto_resolved_from_directory":
        # F-135 / v1.17.7 — фон автоматически найден и вычтен
        bits.append("bg AUTO-SUBTRACTED (F-131)")
    elif bg == "embedded_present_not_subtracted":
        bits.append("bg present, NOT subtracted")
    elif bg == "absent_no_subtraction":
        bits.append("bg NOT subtracted (none avail.)")
    if not bits:
        return None
    return "Inputs: " + " · ".join(bits) + "."


def build_chat_summary(
    result_or_json,
    *,
    json_dict=None,
    report_path: Optional[str] = None,
) -> str:
    """Generate a 3–8 line in-chat summary string.

    Accepts either:
    * a `StagedAnalysisResult` — internally calls `build_json_report`;
    * or a precomputed JSON dict via the `json_dict` keyword.

    Hard upper bound of 8 lines per references/06_report_format.md.
    """
    if json_dict is None:
        # Lazy import to avoid circular dep
        from gamma.reporting.json_report import build_json_report
        json_dict = build_json_report(result_or_json)
    elif json_dict is not None and result_or_json is not None:
        # both provided; trust the dict
        pass
    if json_dict is None:
        return "Spectrum: <empty result>"

    header = json_dict.get("header", {}) or {}
    fn = header.get("filename") or "<unknown>"
    live = header.get("live_time_s")
    dead = header.get("dead_time_pct")
    detector = header.get("detector_canonical") or header.get("detector_type") or "?"
    geometry = header.get("geometry_canonical") or header.get("geometry") or "?"
    env = header.get("environment") or "unknown"

    lines: List[str] = []
    lines.append(
        f"Spectrum: {fn}  ({_fmt_live(live)} live, dead-time {_fmt_dead_time(dead)})"
    )
    lines.append(f"Detector: {detector} / {geometry}  •  environment: {env}")
    # F-89a/b — inputs row: filename binding hypothesis + bg status
    inp_ln = _bg_filename_line(json_dict)
    if inp_ln:
        lines.append(inp_ln)
    # F-88 — chain dominance up front, before nuclide list, because
    # it sets the interpretive context.
    cd_ln = _chain_dominance_line(json_dict)
    if cd_ln:
        lines.append(cd_ln)
    lines.append(_identified_line(json_dict))

    cmp_ln = _completeness_line(json_dict)
    if cmp_ln:
        lines.append(cmp_ln)
    sh_ln = _shielding_line(json_dict)
    if sh_ln:
        lines.append(sh_ln)
    mda_ln = _mda_highlights_line(json_dict)
    if mda_ln:
        lines.append(mda_ln)
    wn_ln = _warnings_line(json_dict)
    if wn_ln:
        lines.append(wn_ln)
    if report_path:
        lines.append(f"Report: {report_path}")

    # Hard upper bound — 8 lines per reference 06.
    return "\n".join(lines[:8])


__all__ = ["build_chat_summary"]
