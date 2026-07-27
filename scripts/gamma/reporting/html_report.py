"""
HTML report renderer (F-86c / v1.15.0).

Renders the same 13 sections as the Markdown report, but as a
self-contained HTML page that can be emailed or opened in any
browser. PNG plots are embedded as base64 ``data:`` URIs so the file
travels as a single artefact (no broken image links).

CSS is inlined (no external dependencies). The styling is
deliberately minimal — readable monospace tables, conservative
colour coding for confirmed / tentative / noise tiers.

API
---

    build_html_report(result_or_json, *, json_dict=None,
                      plots=None) -> str

`plots` accepts the output dict of
:func:`gamma.reporting.plots.build_all_plots`. PNGs are read from
disk and base64-encoded into the HTML. Pass ``plots=None`` for a
text-only report.
"""
from __future__ import annotations

import base64
import html
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

# F-114 / v1.17.3 — canonical interactive report form. The default
# path now produces HTML in the demo skeleton; pass legacy=True to
# build_html_report() to fall back to the static document below.
from gamma.reporting.interactive_html import render_interactive_html


# ──────────────────────────────────────────────────────────────────
# CSS (inlined)
# ──────────────────────────────────────────────────────────────────

_CSS = """
:root {
  --fg: #1a1a1a;
  --muted: #666;
  --accent: #1f3a5f;
  --primary: #cc3322;
  --secondary: #118866;
  --warning: #cc8a00;
  --confirmed: #2e7d32;
  --tentative: #ed6c02;
  --noise: #b71c1c;
  --bg: #fafafa;
  --table-stripe: #f0f0f0;
}
body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
  color: var(--fg);
  background: white;
  max-width: 1100px;
  margin: 1.5em auto;
  padding: 0 1.5em;
  line-height: 1.45;
}
h1 { font-size: 1.6em; border-bottom: 2px solid var(--accent); padding-bottom: 0.3em; }
h2 { font-size: 1.2em; color: var(--accent); margin-top: 1.5em; border-bottom: 1px solid #ddd; padding-bottom: 0.2em; }
h3 { font-size: 1.05em; color: #333; margin-top: 1em; }
table {
  border-collapse: collapse;
  width: 100%;
  font-size: 0.88em;
  margin: 0.5em 0;
}
th, td {
  border: 1px solid #ddd;
  padding: 4px 8px;
  text-align: left;
  vertical-align: top;
}
th { background: var(--bg); font-weight: 600; }
tbody tr:nth-child(even) { background: var(--table-stripe); }
code, pre {
  font-family: "SF Mono", Monaco, Consolas, "Courier New", monospace;
  font-size: 0.92em;
  background: var(--bg);
  padding: 1px 4px;
  border-radius: 3px;
}
.tier-confirmed { color: var(--confirmed); font-weight: 600; }
.tier-tentative { color: var(--tentative); font-weight: 600; }
.tier-noise { color: var(--noise); }
.muted { color: var(--muted); font-style: italic; }
.warning { color: var(--warning); }
.label-row td:first-child { font-weight: 600; width: 30%; }
img.plot { max-width: 100%; height: auto; margin: 0.6em 0; border: 1px solid #ddd; padding: 4px; background: white; }
.footer { margin-top: 2em; padding-top: 0.8em; border-top: 1px solid #ddd; color: var(--muted); font-size: 0.85em; }
.placeholder { padding: 12px; background: var(--bg); border: 1px dashed #ccc; color: var(--muted); }
"""


# ──────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────

def _e(x) -> str:
    """HTML-escape (with None → em-dash)."""
    if x is None:
        return "—"
    return html.escape(str(x))


def _fmt(x, fmt: str = "{:.2f}") -> str:
    if x is None:
        return "—"
    try:
        return html.escape(fmt.format(x))
    except (TypeError, ValueError):
        return html.escape(str(x))


def _yn(b) -> str:
    return "yes" if b else "no"


def _img_data_uri(path: Optional[str]) -> Optional[str]:
    """Base64-encode a PNG (or other image) on disk into a data URI."""
    if not path:
        return None
    p = Path(path)
    if not p.exists():
        return None
    try:
        data = p.read_bytes()
    except OSError:
        return None
    b64 = base64.b64encode(data).decode("ascii")
    suffix = p.suffix.lower().lstrip(".") or "png"
    mime = "image/png" if suffix == "png" else f"image/{suffix}"
    return f"data:{mime};base64,{b64}"


def _kv_table(rows) -> str:
    """Render a list of (label, value) tuples as a 2-column table."""
    parts = ["<table class='label-row'>"]
    for k, v in rows:
        parts.append(f"<tr><td>{_e(k)}</td><td>{_e(v)}</td></tr>")
    parts.append("</table>")
    return "".join(parts)


def _data_table(rows: List[dict], columns: List[tuple]) -> str:
    """Generic table renderer; `columns` is [(header, key, fmt_or_None), ...]."""
    if not rows:
        return "<p class='muted'>(no rows)</p>"
    parts = ["<table><thead><tr>"]
    for col in columns:
        parts.append(f"<th>{_e(col[0])}</th>")
    parts.append("</tr></thead><tbody>")
    for r in rows:
        parts.append("<tr>")
        for header, key, fmt in columns:
            v = r.get(key)
            if v is None:
                cell = "—"
            elif fmt is None:
                cell = _e(v)
            else:
                try:
                    cell = html.escape(fmt.format(v))
                except (TypeError, ValueError):
                    cell = _e(v)
            parts.append(f"<td>{cell}</td>")
        parts.append("</tr>")
    parts.append("</tbody></table>")
    return "".join(parts)


# ──────────────────────────────────────────────────────────────────
# Section renderers (mirror markdown_report.py)
# ──────────────────────────────────────────────────────────────────

def _section_header(h: dict) -> str:
    # F-89a — explicit background status label
    bg_status = h.get("background_status", "")
    bg_label = _BACKGROUND_STATUS_LABEL.get(
        bg_status, f"({bg_status})" if bg_status else "—"
    )
    # F-89b — filename binding hypothesis
    isotope_hints = h.get("filename_isotope_hints") or []
    isotope_hint_str = ", ".join(isotope_hints) if isotope_hints else "—"
    # BUG-39 / BUG-40 — visible detector-fallback row when applicable.
    fb = h.get("detector_fallback") or {}
    fb_reason = fb.get("reason") if isinstance(fb, dict) else ""
    if fb_reason and fb_reason != "profile_loaded_no_fallback":
        fb_label = (
            f"⚠ {fb.get('requested', '?')} → {fb.get('actual', '?')} "
            f"({fb_reason})"
        )
    else:
        fb_label = "—"
    rows = [
        ("Filename",         h.get("filename")),
        ("Sample ID",        h.get("sample_id") or "—"),
        ("Operator",         h.get("operator") or "—"),
        ("Start datetime",   h.get("start_datetime") or "—"),
        ("Live time, s",     _fmt(h.get("live_time_s"), "{:.1f}")),
        ("Real time, s",     _fmt(h.get("real_time_s"), "{:.1f}")),
        ("Dead time, %",     _fmt(h.get("dead_time_pct"), "{:.2f}")),
        ("Detector",         h.get("detector_canonical") or h.get("detector_type") or "—"),
        ("Detector profile fallback", fb_label),
        ("Detector ID",      h.get("detector_id") or "—"),
        ("Geometry",         h.get("geometry") or "—"),
        ("Geometry canon.",  h.get("geometry_canonical") or "—"),
        ("Environment",      h.get("environment")),
        ("Background",       bg_label),
        ("Filename isotope hint", isotope_hint_str),
        ("Analysis mode",    h.get("analysis_mode") or "—"),
        ("Energy ceiling, keV", _fmt(h.get("energy_ceiling_keV"), "{:.0f}")),
        ("Energy max kept, keV", _fmt(h.get("energy_max_keV_kept"), "{:.1f}")),
        ("Channels (kept)",  h.get("n_channels")),
        ("Channels (raw)",   h.get("n_channels_raw")),
        ("Dropped (overflow)", h.get("dropped_high_energy_count")),
    ]
    return "<h2>1. Header</h2>" + _kv_table(rows)


def _section_detector(h: dict, calib: dict) -> str:
    fwhm = calib.get("fwhm_cal", {}) or {}
    fwhm661 = fwhm.get("fwhm_at_661_keV")
    coefs = fwhm.get("coefficients", [])
    coef_str = ", ".join(_fmt(c, "{:.4g}") for c in coefs[:3])
    body = (
        f"<p>Detector canonical: <strong>"
        f"{_e(h.get('detector_canonical') or h.get('detector_type') or '—')}</strong>.</p>"
        f"<p>FWHM(E) model: <code>{_e(fwhm.get('model','?'))}</code> with coefficients "
        f"<code>({coef_str})</code>. FWHM at 661.66 keV ≈ "
        f"<strong>{_fmt(fwhm661, '{:.1f}')} keV</strong> "
        f"(source: {_e(fwhm.get('source','?'))}).</p>"
    )
    return "<h2>2. Detector type</h2>" + body


def _section_calibration(calib: dict) -> str:
    e = calib.get("energy_cal", {}) or {}
    e_coefs = e.get("coefficients", [])
    e_coef_str = ", ".join(_fmt(c, "{:.6g}") for c in e_coefs)
    body = [
        f"<p><strong>Energy calibration</strong>: degree {_e(e.get('degree','?'))} polynomial, "
        f"coefficients <code>({e_coef_str})</code>, source: <code>{_e(e.get('source','?'))}</code>.</p>"
    ]
    slc = calib.get("seven_line_check")
    if slc is not None:
        body.append("<p><strong>7-line ЕРН calibration check (Lsrm §9):</strong></p>")
        body.append(_kv_table([
            ("Lines checked", slc.get("lines_total", 7)),
            ("Lines present", slc.get("lines_present", 0)),
            ("Max residual, keV", _fmt(slc.get("max_residual_keV"), "{:.2f}")),
            ("Max |Δ|/FWHM", _fmt(slc.get("max_residual_fwhm_fraction"), "{:.0%}")),
            ("Calibration quality", slc.get("quality")),
            ("Note", slc.get("quality_note")),
        ]))
    return "<h2>3. Calibration</h2>" + "".join(body)


def _section_identified(rows: List[dict]) -> str:
    if not rows:
        return ("<h2>7. Identified nuclides</h2>"
                "<p class='muted'>(no nuclides confirmed)</p>")
    parts = ["<h2>7. Identified nuclides</h2>"]
    for r in rows:
        tier = r.get("tier", "?")
        tier_class = {
            "confirmed": "tier-confirmed",
            "tentative": "tier-tentative",
            "noise": "tier-noise",
        }.get(tier, "")
        parts.append(
            f"<h3>{_e(r.get('nuclide','?'))}  "
            f"<span class='{tier_class}'>(tier: {_e(tier)})</span></h3>"
        )
        parts.append(_kv_table([
            ("Lines matched", r.get("n_matched_lines")),
            ("Characteristic line, keV", _fmt(r.get("characteristic_line_keV"), "{:.2f}")),
            ("CI", f"{_fmt(r.get('confidence_index'), '{:.2f}')}  ({_e(r.get('confidence_level','?'))})"),
            ("Peak rate, cps", _fmt(r.get("peak_rate_cps"), "{:.3g}")),
            ("Activity, Bq", f"{_fmt(r.get('activity_Bq'), '{:.3g}')} ± {_fmt(r.get('activity_sigma_Bq'), '{:.2g}')}"),
        ]))
        if r.get("specific_activity_Bq_per_kg") is not None:
            parts.append(
                f"<p>Specific activity: "
                f"{_fmt(r.get('specific_activity_Bq_per_kg'), '{:.3g}')} ± "
                f"{_fmt(r.get('specific_activity_sigma_Bq_per_kg'), '{:.2g}')} Bq/kg</p>"
            )
        if r.get("cascade_warning"):
            parts.append(f"<p class='warning'>⚠ Cascade warning: {_e(r['cascade_warning'])}</p>")
        parts.append(f"<p class='muted'>Reason: {_e(r.get('reason','?'))}</p>")
    return "".join(parts)


def _section_diagnostics(diag: dict) -> str:
    rows = [
        ("Measurement environment",      diag.get("measurement_environment", "?")),
        ("Dead time, %",                 _fmt(diag.get("dead_time_pct"), "{:.2f}")),
        ("Dead-time correction applied", _yn(diag.get("dead_time_correction_applied"))),
        ("TCS correction applied",       _yn(diag.get("tcs_correction_applied"))),
        ("Pile-up indicator",            _yn(diag.get("pile_up_indicator"))),
        ("Annihilation 511 observed",    _yn(diag.get("annihilation_511_observed"))),
        ("Escape peaks observed",        diag.get("n_escape_peaks", 0)),
        ("Sum peaks observed",           diag.get("n_sum_peaks", 0)),
        ("XRF residuals",                diag.get("n_xrf_residuals", 0)),
        ("Background subtracted",        _yn(diag.get("background_subtracted"))),
        ("Calibration quality",          diag.get("calibration_quality") or "—"),
        ("Completeness DC, %",           _fmt(diag.get("completeness_dc_pct"), "{:.2f}")),
        ("Completeness flag",            diag.get("completeness_flag") or "—"),
        ("Cascade nuclides flagged",     ", ".join(diag.get("cascade_warning_nuclides") or []) or "—"),
        ("FWHM model source",            diag.get("fwhm_model_source") or "—"),
        ("Efficiency source",            diag.get("efficiency_source") or "—"),
        ("Efficiency loaded",            _yn(diag.get("efficiency_loaded"))),
    ]
    body = ["<h2>12. Diagnostics</h2>", _kv_table(rows)]

    intr = diag.get("intrinsic_activity_signature") or {}
    if intr:
        body.append("<h3>Intrinsic detector activity</h3>")
        body.append(
            f"<p>Detector: <strong>{_e(intr.get('detector_canonical','?'))}</strong> — "
            f"intrinsic Bq/cm³: {_fmt(intr.get('Bq_per_cm3'), '{:.4g}')}.</p>"
        )
        if intr.get("expected_artefacts"):
            body.append("<p>Expected artefacts:</p><ul>")
            for a in intr["expected_artefacts"]:
                body.append(
                    f"<li><strong>{_e(a.get('kind','?'))}</strong> — "
                    f"rule: <code>{_e(a.get('rule','?'))}</code> "
                    f"({_e(a.get('note',''))})</li>"
                )
            body.append("</ul>")
        if intr.get("absent_signatures"):
            body.append(
                f"<p>Absent signatures (sanity check): "
                f"{_e(', '.join(intr['absent_signatures']))}</p>"
            )
    return "".join(body)


_VERSION_HISTORY = [
    ("2026-05", "v1.15.2", "F-89 filename binding hypothesis + chain suppression + bg status"),
    ("2026-05", "v1.15.1", "F-88 user-priority express order + chain-dominance hard-prior"),
    ("2026-05", "v1.15.0", "F-86 plots + HTML + CLI; F-87 anchor-seeding at Step 5α"),
    ("2026-05", "v1.14.0", "Step 11 reporting module (F-85): JSON + chat + Markdown"),
    ("2026-05", "v1.13.0", "Round 5 (F-84): activities, MDA, multiplet deconvolution"),
    ("2026-05", "v1.12.0", "Gamma-1S detector isolation (F-83)"),
]


_BACKGROUND_STATUS_LABEL = {
    "subtracted_from_external_file": "Background SUBTRACTED (external file)",
    # F-135 / v1.17.7 — auto-resolved bg
    "auto_resolved_from_directory":
        "Background SUBTRACTED (auto-resolved F-131: same detector, "
        "compatible geometry, |Δt| ≤ 90 d)",
    "embedded_present_not_subtracted":
        "Background NOT subtracted (embedded bg present but unused)",
    "absent_no_subtraction":
        "Background NOT subtracted (no background available — cps include "
        "natural contribution)",
}


def _section_priority_express(json_dict: dict) -> str:
    """F-88 v1.15.1 — priority express anchors + chain dominance."""
    findings = json_dict.get("priority_express_findings", []) or []
    diag = json_dict.get("diagnostics", {}) or {}
    cd = diag.get("chain_dominance") or {}
    k40_warn = bool(diag.get("k40_ac228_overlap_warning"))

    parts = ["<h2>3α. Priority express anchors (user methodology)</h2>"]

    if not findings:
        parts.append("<p class='muted'>(priority findings not populated — pre-v1.15.1 result)</p>")
        return "".join(parts)

    parts.append(
        "<table><thead><tr>"
        "<th>#</th><th>Signal</th><th>Status</th><th>σ</th><th>Note</th>"
        "</tr></thead><tbody>"
    )
    for pf in findings:
        if pf.get("matched"):
            cls = "tier-confirmed"
            mark = "✔ MATCHED"
        else:
            cls = "muted"
            mark = "✘ missing"
        sig = pf.get("max_significance_sigma")
        sig_str = _fmt(sig, "{:.1f}") if sig else "—"
        parts.append(
            f"<tr>"
            f"<td>{_e(pf.get('order','?'))}</td>"
            f"<td>{_e(pf.get('label','?'))}</td>"
            f"<td class='{cls}'>{mark}</td>"
            f"<td>{sig_str}</td>"
            f"<td>{_e(pf.get('note',''))}</td>"
            f"</tr>"
        )
    parts.append("</tbody></table>")

    # Chain dominance verdict
    parts.append("<h3>Chain dominance verdict</h3>")
    if cd:
        th_cls = "tier-confirmed" if cd.get("th232_dominant") else "muted"
        u_cls = "tier-confirmed" if cd.get("u238_dominant") else "muted"
        parts.append(
            "<ul>"
            f"<li><span class='{th_cls}'>Th-232 chain dominant: "
            f"{'YES' if cd.get('th232_dominant') else 'no'}</span> "
            f"(strength σ ≤ {_fmt(cd.get('th232_strength_sigma'), '{:.1f}')})</li>"
        )
        if cd.get("th232_evidence"):
            parts.append("<li>Th-232 evidence:<ul>")
            for e in cd["th232_evidence"]:
                parts.append(f"<li>{_e(e)}</li>")
            parts.append("</ul></li>")
        parts.append(
            f"<li><span class='{u_cls}'>U-238 chain dominant: "
            f"{'YES' if cd.get('u238_dominant') else 'no'}</span> "
            f"(strength σ ≤ {_fmt(cd.get('u238_strength_sigma'), '{:.1f}')})</li>"
        )
        if cd.get("u238_evidence"):
            parts.append("<li>U-238 evidence:<ul>")
            for e in cd["u238_evidence"]:
                parts.append(f"<li>{_e(e)}</li>")
            parts.append("</ul></li>")
        if cd.get("reason"):
            parts.append(f"<li class='muted'>Reason: {_e(cd['reason'])}</li>")
        parts.append("</ul>")
    else:
        parts.append("<p class='muted'>(chain dominance not computed)</p>")

    if k40_warn:
        parts.append(
            "<div class='placeholder warning' style='border-color:#cc8a00;color:#cc8a00;'>"
            "<strong>⚠ K-40 / Ac-228 overlap warning (F-88)</strong><br/>"
            "Th-232 chain is dominant AND the K-40 1460.82 keV priority "
            "signal matched. On NaI 63×63 the Ac-228 1459.20 keV line "
            "(I=0.85%) cannot be resolved from K-40 — the K-40 peak "
            "area is contaminated by Ac-228 contribution. K-40 activity "
            "should not be reported without deconvolution against the "
            "confirmed Tl-208 2614.51 anchor or a separate Ac-228 "
            "reference."
            "</div>"
        )

    # F-89d — chain suppression notice
    suppressed = cd.get("suppressed_chains") or []
    sup_reason = cd.get("suppression_reason") or ""
    dropped = cd.get("chain_filtered_out_nuclides") or []
    if suppressed:
        parts.append(
            "<div class='placeholder' style='border-color:#1f3a5f;color:#1f3a5f;'>"
            "<strong>ⓘ Chain suppression by filename binding (F-89d)</strong><br/>"
            f"Suppressed chain(s): <strong>{_e(', '.join(suppressed))}</strong><br/>"
        )
        if dropped:
            parts.append(
                f"Nuclides dropped from identifications: "
                f"<strong>{_e(', '.join(dropped))}</strong><br/>"
            )
        if sup_reason:
            parts.append(_e(sup_reason))
        parts.append("</div>")

    return "".join(parts)


def _section_version_history() -> str:
    rows = [
        f"<tr><td>{_e(d)}</td><td>{_e(v)}</td><td>{_e(n)}</td></tr>"
        for d, v, n in _VERSION_HISTORY
    ]
    return (
        "<h2>13. Version history</h2>"
        "<table><thead><tr><th>Date</th><th>Skill version</th><th>Changes</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table>"
    )


# ──────────────────────────────────────────────────────────────────
# Public
# ──────────────────────────────────────────────────────────────────

def build_html_report(
    result_or_json,
    *,
    json_dict: Optional[Dict[str, Any]] = None,
    plots: Optional[Dict[str, Any]] = None,
    legacy: bool = False,
    cost_estimate: Optional[Dict[str, Any]] = None,
    bundle_index: bool = False,
) -> str:
    """Render the JSON report as a self-contained HTML document.

    F-114 / v1.17.3: by default this returns the canonical interactive
    HTML — same skeleton as ``references/demo_contract_v1_17_2/report.html``. Pass
    ``legacy=True`` to fall back to the older static document below
    (kept for tests / debugging only).
    """
    if json_dict is None:
        from gamma.reporting.json_report import build_json_report
        json_dict = build_json_report(result_or_json)

    if not legacy:
        # Canonical form (F-114). Requires the StagedAnalysisResult —
        # if a plain JSON dict was passed instead, fall back to legacy.
        if hasattr(result_or_json, "spec") and hasattr(result_or_json, "fwhm_model"):
            return render_interactive_html(
                json_dict, result_or_json,
                cost_estimate=cost_estimate,
                bundle_index=bundle_index,
            )
        # else: legacy path is the only option

    h = json_dict.get("header", {}) or {}
    calib = json_dict.get("calibration", {}) or {}
    diag = json_dict.get("diagnostics", {}) or {}
    cmp = json_dict.get("completeness", {}) or {}

    title = f"Gamma report — {_e(h.get('filename','spectrum'))}"

    body: List[str] = []
    # F-UX-03 / 2026-06-04 — back-link «← На главную» убрана (supersedes
    # F-RPT-02). См. _state/agent_a/inbox/2026-06-04_correction_9_*.
    body.append(f"<h1>Gamma-spectrum analysis report</h1>")
    body.append(
        f"<p class='muted'>skill {_e(json_dict.get('skill_version','?'))}, "
        f"schema {_e(json_dict.get('schema_version','?'))}</p>"
    )

    body.append(_section_header(h))
    body.append(_section_detector(h, calib))
    body.append(_section_calibration(calib))
    body.append(_section_priority_express(json_dict))

    # Section 4 — Primary FEPs (Основные пики полного поглощения)
    body.append("<h2>4. Основные пики полного поглощения</h2>")
    pf_rows_html = []
    for row in (json_dict.get("primary_feps", []) or []):
        row = dict(row)
        if row.get("is_upper_limit_artifact"):
            sigma_val = row.get("peak_area_sigma")
            if isinstance(sigma_val, (int, float)) and sigma_val > 0:
                row["peak_area_counts"] = f"< {sigma_val:.3g} (<MDA)"
            else:
                row["peak_area_counts"] = "<MDA"
            row["rate_cps"] = "—"
        pf_rows_html.append(row)
    body.append(_data_table(
        pf_rows_html,
        [
            ("Nuclide",      "nuclide",           None),
            ("E_lib, keV",   "library_E_keV",     "{:.2f}"),
            ("I_lib, %",     "library_I_pct",     "{:.3g}"),
            ("E_peak, keV",  "peak_E_keV",        "{:.2f}"),
            ("FWHM, keV",    "fwhm_keV",          "{:.2f}"),
            ("Area",         "peak_area_counts",  "{:.3g}"),
            ("σ_area",       "peak_area_sigma",   "{:.3g}"),
            ("Rate, cps",    "rate_cps",          "{:.3g}"),
            ("Source",       "peak_area_source",  None),
            ("Char.",        "is_characteristic", None),
        ],
    ))

    # Section 5 — Secondary peaks (Вторичные пики)
    body.append("<h2>5. Вторичные пики</h2>")
    body.append(_data_table(
        json_dict.get("secondary_peaks", []) or [],
        [
            ("E, keV",        "energy_keV",      "{:.2f}"),
            ("σ",             "significance",    "{:.1f}"),
            ("Type",          "type",            None),
            ("Feature",       "feature_kind",    None),
            ("Parent",        "parent_nuclide",  None),
            ("Parent E",      "parent_line_keV", "{:.2f}"),
            ("Note",          "note",            None),
        ],
    ))

    # Section 6 — XRF
    body.append("<h2>6. Elemental XRF</h2>")
    xrf = json_dict.get("elemental_xrf", []) or []
    if not xrf:
        body.append("<p class='muted'>(no XRF residuals classified)</p>")
    else:
        for entry in xrf:
            body.append(
                f"<h3>{_e(entry.get('element','?'))}  "
                f"<span class='muted'>({_e(entry.get('mechanism','?'))})</span></h3>"
            )
            body.append(f"<p>Observed lines: {_e(entry.get('n_observed', 0))}</p>")
            if entry.get("observed_lines"):
                body.append("<ul>")
                for ln in entry["observed_lines"]:
                    body.append(
                        f"<li>{_fmt(ln.get('energy_keV'), '{:.2f}')} keV "
                        f"(σ={_fmt(ln.get('significance'), '{:.1f}')}, "
                        f"lib {_fmt(ln.get('library_E_keV'), '{:.2f}')}, "
                        f"Δ={_fmt(ln.get('delta_keV'), '{:.2f}')})</li>"
                    )
                body.append("</ul>")

    # Section 7 — Identified nuclides
    body.append(_section_identified(json_dict.get("identified_nuclides", []) or []))

    # Section 8 — Unidentified + DC%
    body.append("<h2>8. Unidentified significant peaks</h2>")
    body.append(_data_table(
        json_dict.get("unidentified_peaks", []) or [],
        [
            ("E, keV",       "energy_keV",   "{:.2f}"),
            ("σ",            "significance", "{:.1f}"),
            ("Label",        "label",        None),
            ("Note",         "note",         None),
        ],
    ))
    body.append(
        f"<p><strong>Dose Contribution (DC)</strong> = "
        f"{_fmt(cmp.get('dc_pct'), '{:.2f}')} %  "
        f"[{_e(cmp.get('flag') or '—')}]</p>"
    )

    # Section 9 — Spectrum plot
    body.append("<h2>9. Spectrum plot</h2>")
    spec_uri = _img_data_uri((plots or {}).get("spectrum") if plots else None)
    if spec_uri:
        body.append(f"<img class='plot' src='{spec_uri}' alt='Spectrum overlay'/>")
    else:
        body.append(
            "<p class='placeholder'>Plot generation deferred — "
            "call <code>build_report(..., write_plots=True)</code> to render PNGs.</p>"
        )

    # Section 10 — Multiplet deconvolution (Разложение мультиплетов)
    # BUG-5 / v1.18.30+ (Agent B): consistency с interactive_html и markdown.
    # Legacy путь рендерит только sample-multiplets.
    # BUG-14 / v1.18.30+ (Agent B): для чисто-фоновых спектров
    # (diagnostics.measurement_environment == "background_only") секция
    # пропускается — пики не значимы над шумовым континуумом. Симметрично
    # с interactive_html и markdown_report (BUG-14).
    _is_bg_only_legacy = (diag.get("measurement_environment") == "background_only")
    if not _is_bg_only_legacy:
        body.append("<h2>10. Мультиплеты — разложение в спектре образца</h2>")
        body.append(
            "<p class='muted'>Первичная подгонка по библиотечным интенсивностям. "
            "Площади взяты из спектра образца — основной источник данных для "
            "расчёта активности.</p>"
        )
        decons = json_dict.get("multiplet_deconvolutions", []) or []
        mult_paths = (plots or {}).get("multiplets") if plots else None
        if not decons:
            body.append("<p class='muted'>(no multiplet clusters processed in this run)</p>")
        else:
            for i, d in enumerate(decons, start=1):
                body.append(
                    f"<h3>Cluster {i}  "
                    f"(χ²/dof = {_fmt(d.get('chi2_per_dof'), '{:.3f}')}, "
                    f"converged: {_yn(d.get('converged'))})</h3>"
                )
                if mult_paths and i - 1 < len(mult_paths):
                    uri = _img_data_uri(mult_paths[i - 1])
                    if uri:
                        body.append(f"<img class='plot' src='{uri}' alt='Cluster {i}'/>")
                body.append(_data_table(
                    d.get("components", []) or [],
                    [
                        ("Nuclide",  "nuclide",                None),
                        ("E_lib",    "line_E_keV",             "{:.2f}"),
                        ("Area",     "deconvolved_area",       "{:.3g}"),
                        ("σ_area",   "deconvolved_area_sigma", "{:.3g}"),
                    ],
                ))

    # Section 11 — MDA table
    body.append("<h2>11. MDA table (ISO 11929 / Lsrm §6.3)</h2>")
    body.append(_data_table(
        json_dict.get("mda", []) or [],
        [
            ("Nuclide",       "nuclide",                   None),
            ("E, keV",        "line_E_keV",                "{:.2f}"),
            ("I, %",          "intensity_pct",             "{:.3g}"),
            ("ε",             "efficiency",                "{:.3g}"),
            ("L_C, counts",   "decision_threshold_counts", "{:.3g}"),
            ("L_D, counts",   "detection_limit_counts",    "{:.3g}"),
            ("MDA, Bq",       "MDA_Bq",                    "{:.3g}"),
        ],
    ))

    # Section 12 — Diagnostics
    body.append(_section_diagnostics(diag))

    # Warnings
    warns = json_dict.get("warnings", []) or []
    if warns:
        body.append("<h2>Warnings</h2><ul>")
        for w in warns:
            # P0-8: structured alarm entries (dict) rendered via message field.
            if isinstance(w, dict):
                w_text = w.get("message") or w.get("code") or str(w)
            else:
                w_text = w
            body.append(f"<li class='warning'>{_e(w_text)}</li>")
        body.append("</ul>")

    # Section 13 — Version history
    body.append(_section_version_history())
    body.append("<p class='footer'>Generated by gamma.reporting (v1.15.0).</p>")

    html_doc = (
        "<!DOCTYPE html>\n"
        "<html lang='en'>\n<head>\n"
        "<meta charset='utf-8'/>\n"
        f"<title>{title}</title>\n"
        f"<style>{_CSS}</style>\n"
        "</head>\n<body>\n"
        + "\n".join(body)
        + "\n</body>\n</html>\n"
    )
    return html_doc


__all__ = ["build_html_report"]
