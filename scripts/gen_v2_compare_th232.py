# -*- coding: utf-8 -*-
from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""v1.18.25.0 — Th-232 v2 vs production compare generator (post-hoc).

F-390 (v1.18.25.3+): 2-column layout.
  * Шапка (meta + KPI) — full width.
  * Основная часть — `display:grid; grid-template-columns:1fr 1fr`:
      колонка V2 (слева) | колонка Production (справа).
    Внутри каждой колонки: identified nuclides chips, primary FEP
    peaks table, multiplet clusters list, secondary peaks list.
  * Финальная diff-секция (matched, only_v2, only_prod) — full width.

Читает два готовых JSON-отчёта (production и V2) из demo_reports/<run>/
и собирает:
  * compare_data.json — KPI + peaks + clusters + nuclides diff
  * v2_compare_report.html — 2-column report (см. эталон
    `demo_reports/v1_18_24_th232_full/v2_compare/v2_compare_report.html`).

Использование:
    python scripts/gen_v2_compare_th232.py <demo_run_dir>

где <demo_run_dir> содержит sample/ + sample_v2/ с *_report.json.

F-380 contract (НЕ менять): primary_feps использует ключи
peak_channel, peak_E_keV, peak_area_counts, peak_area_sigma.
"""
import html as _html
import json
import sys
from pathlib import Path


def _load_report(report_dir: Path) -> dict:
    candidates = list(report_dir.glob("*_report.json"))
    if not candidates:
        raise FileNotFoundError(f"No *_report.json in {report_dir}")
    with candidates[0].open("r", encoding="utf-8") as f:
        return json.load(f)


def _peak_rows(peaks: list[dict]) -> list[dict]:
    """F-380 / v1.18.25.2 — primary_feps реальные ключи:
    peak_channel, peak_E_keV, peak_area_counts, peak_area_sigma.
    Старая версия читала channel/E_keV/significance → все 0/None.

    F-391 / v1.18.27 — filter phantom anchors (peak_area_source
    library_anchor / library_anchor_phantom) ДО построения rows. Они
    представляют library evidence без реального signal — место в
    multiplet evidence, не в primary peaks list.

    BUG-11 / v1.18.30+ (Agent B) — dedupe по (nuclide, peak_channel).
    Multiplet fitter может назначать несколько library-линий на один
    канал (e.g. Ac-228 ch=203 → lib_E 562.5 / 570.91 / 572.14 / 583.41,
    все с одинаковым peak_E_keV=583.01). Старый ключ
    (nuclide, round(E_keV, 0)) корректно сворачивал такие случаи только
    благодаря тому, что peak_E_keV одинаковый, но это совпадение, не
    инвариант. peak_channel — единственная стабильная identity-метка
    того, что речь идёт о ОДНОМ физическом канале. Tiebreak: max(S/σ),
    затем max(library_I_pct) — берём строку с реальным сигналом, а при
    равных σ — наиболее вероятную линию по библиотеке.
    """
    PHANTOM = {"library_anchor", "library_anchor_phantom"}
    # F-391 — filter phantoms
    visible = [
        p for p in peaks
        if str(p.get("peak_area_source") or "") not in PHANTOM
    ]
    rows = []
    for p in visible:
        ch = p.get("peak_channel", p.get("channel"))
        E = p.get("peak_E_keV", p.get("E_keV", 0.0))
        # σ = peak_area / peak_area_sigma как S/σ ratio
        area = p.get("peak_area_counts") or 0.0
        sigma_area = p.get("peak_area_sigma") or 1.0
        sig = (float(area) / float(sigma_area)) if sigma_area else 0.0
        rows.append({
            "channel": int(ch) if ch is not None else None,
            "E_keV": round(float(E or 0.0), 2),
            "sigma": round(float(sig), 1),
            "nuclide": str(p.get("nuclide") or ""),
            "library_E_keV": round(float(p.get("library_E_keV") or 0.0), 2),
            "library_I_pct": round(float(p.get("library_I_pct") or 0.0), 3),
        })
    # BUG-11 — dedupe по (nuclide, peak_channel). Tiebreak: max(S/σ),
    # затем max(library_I_pct). Ключ-fallback на round(E_keV) если channel
    # отсутствует (legacy primary_feps без peak_channel).
    dedup: dict = {}
    for r in rows:
        ch_key = r["channel"] if r["channel"] is not None else int(round(r["E_keV"]))
        key = (r["nuclide"], ch_key)
        prev = dedup.get(key)
        if prev is None:
            dedup[key] = r
            continue
        # Tiebreak: max(sigma), затем max(library_I_pct).
        if (r["sigma"], r["library_I_pct"]) > (prev["sigma"], prev["library_I_pct"]):
            dedup[key] = r
    return list(dedup.values())


def _match_peaks(prod_peaks: list[dict], v2_peaks: list[dict], tol_keV: float = 4.0):
    matched, only_prod, only_v2 = [], [], []
    used_v2 = set()
    for p in prod_peaks:
        Ep = float(p.get("E_keV", 0.0))
        best_j, best_dE = -1, tol_keV * 1.5
        for j, q in enumerate(v2_peaks):
            if j in used_v2:
                continue
            dE = abs(Ep - float(q.get("E_keV", 0.0)))
            if dE < best_dE:
                best_dE, best_j = dE, j
        if best_j >= 0 and best_dE <= tol_keV:
            matched.append({"prod": p, "v2": v2_peaks[best_j], "dE": round(best_dE, 2)})
            used_v2.add(best_j)
        else:
            only_prod.append(p)
    for j, q in enumerate(v2_peaks):
        if j not in used_v2:
            only_v2.append(q)
    return matched, only_prod, only_v2


def _nuclides_set(report: dict) -> list[str]:
    return sorted({(n.get("nuclide") or n.get("name") or "").strip()
                   for n in (report.get("identified_nuclides") or [])
                   if n.get("nuclide") or n.get("name")})


# ──────────────────────────────────────────────────────────────────
# F-390 HTML render helpers — single-column-of-column blocks
# ──────────────────────────────────────────────────────────────────

def _esc(x) -> str:
    return _html.escape(str(x), quote=True)


def _chips_html(items: list[str], badge_bg: str = "#37a") -> str:
    if not items:
        return "<span class='dim'>—</span>"
    return " ".join(
        f"<span class='chip' style='background:{badge_bg};'>{_esc(n)}</span>"
        for n in items
    )


def _primary_peaks_table_html(rows: list[dict]) -> str:
    """Compact primary FEP table for column-internal use.

    F-393 / v1.18.27 — pre-sorted by E asc + sortable headers (num/str hints).
    Default sort marker (.sorted.asc) на колонке E применяется JS-кодом
    в footer compare-HTML.
    """
    if not rows:
        return "<div class='dim'>нет пиков полного поглощения</div>"
    rows_sorted = sorted(rows, key=lambda r: float(r.get("E_keV") or 0.0))
    body = "\n".join(
        f"<tr><td>{r['channel']}</td>"
        f"<td>{r['E_keV']:.2f}</td>"
        f"<td>{_esc(r.get('nuclide') or '—')}</td>"
        f"<td>{r['library_E_keV']:.2f}</td>"
        f"<td>{r['sigma']:.1f}</td></tr>"
        for r in rows_sorted
    )
    return (
        "<table class='peaks' data-sortable='true'>"
        "<thead><tr>"
        "<th data-sort='num' data-sortable='true' onclick=\"sortTable(this)\">канал</th>"
        "<th data-sort='num' data-sortable='true' onclick=\"sortTable(this)\""
        " data-default-sort='asc'>E, кэВ</th>"
        "<th data-sort='str' data-sortable='true' onclick=\"sortTable(this)\">нуклид</th>"
        "<th data-sort='num' data-sortable='true' onclick=\"sortTable(this)\">E_lib</th>"
        "<th data-sort='num' data-sortable='true' onclick=\"sortTable(this)\">S/σ</th>"
        "</tr></thead>"
        f"<tbody>{body}</tbody></table>"
    )


def _multiplet_clusters_html(clusters: list[dict], badge_label: str,
                             badge_bg: str) -> str:
    """Render multiplet clusters block (один <div class='cluster'> на каждый)."""
    if not clusters:
        return "<div class='dim'>нет кластеров мультиплетов</div>"
    out = []
    for c in clusters:
        cid = _esc(c.get("cluster_id") or "")
        converged = bool(c.get("converged"))
        chi2 = c.get("chi2_per_dof")
        closure = c.get("closure_pct") or c.get("closure")
        E_lo = c.get("E_lo_keV") or c.get("e_lo_keV")
        E_hi = c.get("E_hi_keV") or c.get("e_hi_keV")
        range_str = (f"{E_lo:.1f} – {E_hi:.1f} кэВ"
                     if (E_lo is not None and E_hi is not None) else "")
        chi2_str = f"χ²/ν = <b>{chi2:.2f}</b>" if isinstance(chi2, (int, float)) else ""
        closure_str = (f"closure = <b>{closure}%</b>"
                       if closure not in (None, "") else "")
        conv_badge = ("<span class='ok'>✓ сошёлся</span>"
                      if converged else "<span class='bad'>✗ не сошёлся</span>")
        comp_rows = []
        # F-391 / v1.18.27 — phantom anchors помечаем visual class 'phantom'
        # (зачёркивается + dim) и оставляем для evidence: пользователь видит,
        # что library line была учтена как anchor, но не fit'ena как
        # отдельный пик. Для multi-active кластеров это полезно показать.
        PHANTOM_SRC = {"library_anchor", "library_anchor_phantom"}
        for comp in (c.get("components") or []):
            E_lib = comp.get("line_E_keV") or comp.get("E_keV") or 0.0
            area = comp.get("deconvolved_area") or comp.get("area") or 0.0
            sigma = (comp.get("deconvolved_area_sigma")
                     or comp.get("area_sigma") or comp.get("sigma_area") or 0.0)
            I_pct = comp.get("library_I_pct") or comp.get("I_pct")
            I_str = f"{I_pct:.2f}" if isinstance(I_pct, (int, float)) else "—"
            is_phantom = str(
                comp.get("peak_area_source") or ""
            ) in PHANTOM_SRC
            tr_cls = " class='phantom'" if is_phantom else ""
            nuc_label = _esc(comp.get('nuclide') or '?')
            if is_phantom:
                # F-397.3 / v1.18.28 (Agent B) — фикс title attribute:
                # старый вариант ("не fit'ena, evidence-only") содержал
                # апостроф внутри одинарных кавычек title='...', что
                # обрывало tooltip на "не fit" + ломало parsing атрибутов.
                # Использую двойные кавычки для title + чистый RU текст.
                nuc_label = (
                    f"{nuc_label} "
                    '<span class="dim" title="library anchor — '
                    'не подгоняется, evidence-only">(якорь)</span>'
                )
            comp_rows.append(
                f"<tr{tr_cls}><td>{nuc_label}</td>"
                f"<td>{float(E_lib):.2f}</td>"
                f"<td>{I_str}</td>"
                f"<td>{int(float(area))}</td>"
                f"<td>{int(float(sigma))}</td></tr>"
            )
        # F-393 / v1.18.27 — sortable component table
        # v1.27.9 lane A: extract empty-row fallback out of f-string expression
        # to avoid PEP-701 backslash-in-fstring (fails on Python 3.11 CI).
        empty_row = '<tr><td colspan="5" class="dim">—</td></tr>'
        comp_body = "".join(comp_rows) or empty_row
        comp_table = (
            "<table class='comp' data-sortable='true'>"
            "<thead><tr>"
            "<th data-sort='str' data-sortable='true' onclick=\"sortTable(this)\">Нуклид</th>"
            "<th data-sort='num' data-sortable='true' onclick=\"sortTable(this)\""
            " data-default-sort='asc'>E, кэВ</th>"
            "<th data-sort='num' data-sortable='true' onclick=\"sortTable(this)\">I, %</th>"
            "<th data-sort='num' data-sortable='true' onclick=\"sortTable(this)\">area</th>"
            "<th data-sort='num' data-sortable='true' onclick=\"sortTable(this)\">σ_area</th>"
            "</tr></thead>"
            f"<tbody>{comp_body}</tbody></table>"
        )
        # F-397.3 / v1.18.28 (Agent B) — suppress empty <b></b> tag when
        # cluster_id missing (V2 phantom-/wide-CC кластеры не получают
        # стабильного M-имени и отдавали `<b></b>` в выводе).
        cid_html = f"<b>{cid}</b> " if cid else ""
        out.append(
            f"<div class='cluster'>"
            f"<div class='cluster-head'>"
            f"<span class='badge' style='background:{badge_bg}'>{_esc(badge_label)}</span> "
            f"{cid_html}"
            f"<span class='dim'>{range_str}</span> "
            f"<span class='dim'>{chi2_str}</span> "
            f"<span class='dim'>{closure_str}</span> "
            f"{conv_badge}</div>"
            f"{comp_table}</div>"
        )
    return "\n".join(out)


def _secondary_peaks_html(items: list[dict]) -> str:
    """Render compact secondary peaks list (Compton / SE / DE / annihilation).

    F-393 / v1.18.27 — pre-sorted by E asc + sortable headers (num/str hints).
    """
    if not items:
        return "<div class='dim'>нет вторичных пиков</div>"
    items_sorted = sorted(
        items,
        key=lambda s: float(s.get("energy_keV") or s.get("E_keV") or 0.0),
    )
    rows = []
    for s in items_sorted:
        ch = s.get("channel")
        E = s.get("energy_keV") or s.get("E_keV") or 0.0
        sig = s.get("significance") or 0.0
        kind = s.get("feature_kind") or s.get("type") or "—"
        parent = s.get("parent_nuclide") or "—"
        parent_E = s.get("parent_line_keV")
        parent_str = (f"{_esc(parent)} ({parent_E:.1f})"
                      if isinstance(parent_E, (int, float)) else _esc(parent))
        rows.append(
            f"<tr><td>{ch if ch is not None else '—'}</td>"
            f"<td>{float(E):.2f}</td>"
            f"<td>{_esc(kind)}</td>"
            f"<td>{parent_str}</td>"
            f"<td>{float(sig):.1f}</td></tr>"
        )
    return (
        "<table class='peaks' data-sortable='true'>"
        "<thead><tr>"
        "<th data-sort='num' data-sortable='true' onclick=\"sortTable(this)\">канал</th>"
        "<th data-sort='num' data-sortable='true' onclick=\"sortTable(this)\""
        " data-default-sort='asc'>E, кэВ</th>"
        "<th data-sort='str' data-sortable='true' onclick=\"sortTable(this)\">тип</th>"
        "<th data-sort='str' data-sortable='true' onclick=\"sortTable(this)\">родитель</th>"
        "<th data-sort='num' data-sortable='true' onclick=\"sortTable(this)\">S/σ</th>"
        "</tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table>"
    )


# ──────────────────────────────────────────────────────────────────
# diff-section helpers (full-width footer block)
# ──────────────────────────────────────────────────────────────────

def _diff_peak_row_html(p: dict) -> str:
    nuc = _esc(p.get('nuclide') or '?')
    E_lib = (f"{p.get('library_E_keV', 0):.2f}"
             if p.get('library_E_keV') else "—")
    I_str = (f"{p.get('library_I_pct', 0):.2f}"
             if p.get('library_I_pct') else "—")
    return (
        f"<tr><td>{p['channel']}</td><td>{p['E_keV']:.2f}</td>"
        f"<td><span class='chip' style='background:#37a;font-size:10px;'>{nuc}</span></td>"
        f"<td>{E_lib}</td><td>{I_str}</td>"
        f"<td>{p['sigma']:.1f}</td></tr>"
    )


def main():
    if len(sys.argv) < 2:
        print("usage: gen_v2_compare_th232.py <demo_run_dir>")
        sys.exit(1)
    run = Path(sys.argv[1]).resolve()
    prod_dir = run / "sample"
    v2_dir = run / "sample_v2"
    out_dir = run / "v2_compare"
    out_dir.mkdir(parents=True, exist_ok=True)

    prod = _load_report(prod_dir)
    v2 = _load_report(v2_dir)

    prod_peaks_raw = prod.get("primary_feps") or []
    v2_peaks_raw = v2.get("primary_feps") or []
    prod_peaks = _peak_rows(prod_peaks_raw)
    v2_peaks = _peak_rows(v2_peaks_raw)

    matched, only_prod, only_v2 = _match_peaks(prod_peaks, v2_peaks)

    prod_nuc = _nuclides_set(prod)
    v2_nuc = _nuclides_set(v2)
    only_v2_nuc = sorted(set(v2_nuc) - set(prod_nuc))
    only_prod_nuc = sorted(set(prod_nuc) - set(v2_nuc))

    prod_mult_list = prod.get("multiplet_deconvolutions") or []
    v2_mult_list = v2.get("multiplet_deconvolutions") or []
    prod_sec_list = prod.get("secondary_peaks") or []
    v2_sec_list = v2.get("secondary_peaks") or []
    prod_mult = len(prod_mult_list)
    v2_mult = len(v2_mult_list)
    prod_sec = len(prod_sec_list)
    v2_sec = len(v2_sec_list)
    prod_un = len(prod.get("unidentified_peaks") or [])
    v2_un = len(v2.get("unidentified_peaks") or [])

    header = prod.get("header") or {}
    sample_name = header.get("sample_filename") or header.get("filename") or "—"
    bg_name = header.get("background_filename") or "—"
    mass_kg = (prod.get("calibration") or {}).get("sample_mass_kg") \
        or header.get("sample_mass_kg") \
        or v2.get("header", {}).get("sample_mass_kg") \
        or 1.6
    skill_v = prod.get("skill_version") or "v1.18.25.0"

    compare = {
        "meta": {
            "sample": sample_name,
            "background": bg_name,
            "mass_kg": float(mass_kg) if mass_kg else None,
            "skill_version": skill_v,
            "production_n_peaks": len(prod_peaks),
            "production_n_multiplets": prod_mult,
            "production_n_secondary": prod_sec,
            "production_n_unidentified": prod_un,
            "v2_n_peaks": len(v2_peaks),
            "v2_n_multiplets": v2_mult,
            "v2_n_secondary": v2_sec,
            "v2_n_unidentified": v2_un,
            "nuclides_production": prod_nuc,
            "nuclides_v2": v2_nuc,
            "nuclides_only_in_v2": only_v2_nuc,
            "nuclides_only_in_production": only_prod_nuc,
            "peak_match_tol_keV": 4.0,
            "n_matched": len(matched),
            "n_only_prod": len(only_prod),
            "n_only_v2": len(only_v2),
        },
        "production_peaks": prod_peaks,
        "v2_peaks": v2_peaks,
        "matched_peaks": matched,
        "only_in_production": only_prod,
        "only_in_v2": only_v2,
    }

    with (out_dir / "compare_data.json").open("w", encoding="utf-8") as f:
        json.dump(compare, f, ensure_ascii=False, indent=2)

    # ─────────────── HTML render ───────────────
    # F-393 / v1.18.27 — diff tables по умолчанию sorted by E asc, чтобы
    # совпадали с visual marker (.sorted.asc) на колонке E.
    matched_sorted = sorted(
        matched, key=lambda m: float(m.get("prod", {}).get("E_keV") or 0.0)
    )
    only_prod_sorted = sorted(only_prod, key=lambda p: float(p.get("E_keV") or 0.0))
    only_v2_sorted = sorted(only_v2, key=lambda p: float(p.get("E_keV") or 0.0))

    matched_rows = "\n".join(
        f"<tr><td>{m['prod']['E_keV']:.2f}</td>"
        f"<td>{m['v2']['E_keV']:.2f}</td>"
        f"<td>{m['dE']:.2f}</td>"
        f"<td>{m['prod']['sigma']:.1f}</td>"
        f"<td>{m['v2']['sigma']:.1f}</td></tr>"
        for m in matched_sorted
    )

    only_prod_html = "\n".join(_diff_peak_row_html(p) for p in only_prod_sorted) or \
        "<tr><td colspan='6' class='dim'>нет</td></tr>"
    only_v2_html = "\n".join(_diff_peak_row_html(p) for p in only_v2_sorted) or \
        "<tr><td colspan='6' class='dim'>нет</td></tr>"

    # BUG-12 / v1.18.30 — выравнивание секций «Сводка по конвейерам» по
    # горизонтали. Раньше V2 и Production рендерились как два независимых
    # <div class='col'>, и при разной длине списков (V2: 21 пик vs Prod: 12)
    # заголовок «Кластеры мультиплетов» в V2 уезжал значительно ниже, чем
    # в Production — карточки M1↔M1, M2↔M2 не сопоставлялись визуально.
    #
    # Новая структура: единый CSS Grid `pipeline-summary-grid` с 5 рядами,
    # где каждая из 4 секций (nuclides / primary peaks / multiplets /
    # secondary) занимает 1 ряд из 2 ячеек (V2 cell слева, Prod cell справа),
    # и над ними — ряд с заголовками колонок (V2 / Production). Так как
    # высота ряда = max(height(v2 cell), height(prod cell)), заголовок
    # каждой секции автоматически выровнен между колонками. См.
    # `.pipeline-summary-grid` в CSS ниже.
    #
    # F-395: все user-facing заголовки секций — RU.
    def _ps_cell(col_kind: str, section_id: str, header_html: str,
                 content_html: str) -> str:
        """Render one section-cell: header + content in a single grid cell.

        col_kind: 'v2' | 'prod'. Используется для accent-border и для
        data-col атрибута (тесты опираются на него).
        section_id: 'nuclides' | 'primary' | 'multiplets' | 'secondary'.
        """
        return (
            f"<div class='ps-cell ps-cell-{col_kind}' "
            f"data-section='{section_id}' data-col='{col_kind}'>"
            f"<h3 class='ps-section-header' "
            f"data-section='{section_id}' data-col='{col_kind}'>"
            f"{header_html}</h3>"
            f"<div class='ps-section-content'>{content_html}</div>"
            f"</div>"
        )

    # Column title cells (row 1)
    ps_title_v2 = (
        "<div class='ps-col-title ps-cell-v2' data-col='v2'>"
        "V2 (experimental)</div>"
    )
    ps_title_prod = (
        "<div class='ps-col-title ps-cell-prod' data-col='prod'>"
        "Production</div>"
    )

    # Rows 2-5: per-section cells
    ps_v2_nuclides = _ps_cell(
        "v2", "nuclides",
        f"Идентифицированные нуклиды ({len(v2_nuc)})",
        f"<div class='chips-row'>{_chips_html(v2_nuc, '#3a7')}</div>",
    )
    ps_prod_nuclides = _ps_cell(
        "prod", "nuclides",
        f"Идентифицированные нуклиды ({len(prod_nuc)})",
        f"<div class='chips-row'>{_chips_html(prod_nuc, '#37a')}</div>",
    )
    ps_v2_peaks = _ps_cell(
        "v2", "primary",
        f"Основные пики полного поглощения ({len(v2_peaks)})",
        _primary_peaks_table_html(v2_peaks),
    )
    ps_prod_peaks = _ps_cell(
        "prod", "primary",
        f"Основные пики полного поглощения ({len(prod_peaks)})",
        _primary_peaks_table_html(prod_peaks),
    )
    ps_v2_mult = _ps_cell(
        "v2", "multiplets",
        f"Кластеры мультиплетов ({v2_mult})",
        _multiplet_clusters_html(v2_mult_list, 'V2', '#3a7'),
    )
    ps_prod_mult = _ps_cell(
        "prod", "multiplets",
        f"Кластеры мультиплетов ({prod_mult})",
        _multiplet_clusters_html(prod_mult_list, 'PROD', '#37a'),
    )
    ps_v2_sec = _ps_cell(
        "v2", "secondary",
        f"Вторичные пики ({v2_sec})",
        _secondary_peaks_html(v2_sec_list),
    )
    ps_prod_sec = _ps_cell(
        "prod", "secondary",
        f"Вторичные пики ({prod_sec})",
        _secondary_peaks_html(prod_sec_list),
    )

    html = f"""<!doctype html>
<html lang="ru"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>v2 vs Production — Th-232 (F-390) — {skill_v}</title>
<style>
:root {{
  --bg:#fbfaf3; --bg-secondary:#f5f4ee;
  --text:#1a1a1a; --text-secondary:#5f5e5a;
  --border:rgba(0,0,0,.18); --border-secondary:rgba(0,0,0,.20);
  --radius-md:6px; --accent-v2:#3a7; --accent-prod:#37a;
}}
@media (prefers-color-scheme:dark) {{
  :root {{ --bg:#1c1c1b; --bg-secondary:#262625;
           --text:#ece9d8; --text-secondary:#b4b2a9;
           --border:rgba(255,255,255,.22); --border-secondary:rgba(255,255,255,.22); }}
}}
body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
     background:var(--bg);color:var(--text);margin:0;padding:0;}}
.page{{max-width:1280px;margin:0 auto;padding:24px;}}
h1{{font-size:22px;margin:0 0 6px;}}
.sub{{font-size:13px;color:var(--text-secondary);margin:0 0 24px;}}
.meta{{border:0.5px solid var(--border-secondary);border-radius:var(--radius-md);
       padding:14px 18px;margin-bottom:18px;background:var(--bg-secondary);}}
.meta table{{width:100%;border-collapse:collapse;font-size:13px;}}
.meta th{{text-align:left;font-weight:500;color:var(--text-secondary);
         padding:3px 12px 3px 0;width:30%;}}
.meta td{{padding:3px 0;}}
h2{{font-size:17px;margin:28px 0 12px;border-bottom:0.5px solid var(--border-secondary);padding-bottom:6px;text-align:left;padding-left:0;vertical-align:middle;}}
h3{{font-size:13px;margin:14px 0 6px;color:var(--text-secondary);font-weight:500;
    text-transform:uppercase;letter-spacing:.5px;text-align:left;padding-left:0;vertical-align:middle;}}
.grid{{display:grid;grid-template-columns:1fr 1fr;gap:18px;}}
.col{{border:0.5px solid var(--border-secondary);border-radius:var(--radius-md);
      padding:14px 16px;background:var(--bg-secondary);}}
.col.v2{{border-left:3px solid var(--accent-v2);}}
.col.prod{{border-left:3px solid var(--accent-prod);}}
.col-title{{font-size:16px;font-weight:600;margin:0 0 12px;padding-bottom:8px;
            border-bottom:0.5px solid var(--border-secondary);}}
.col.v2 .col-title{{color:var(--accent-v2);}}
.col.prod .col-title{{color:var(--accent-prod);}}
/* BUG-12 / v1.18.30 — pipeline-summary-grid: aligns V2/Prod section
   headers по горизонтали. Каждая секция (nuclides / primary peaks /
   multiplets / secondary peaks) — одна строка из двух ячеек. Высота
   строки = max(V2 cell, Prod cell), поэтому следующий заголовок секции
   стартует в обеих колонках синхронно даже при разной длине списков. */
.pipeline-summary-grid{{display:grid;grid-template-columns:1fr 1fr;
  gap:0 18px;margin-bottom:18px;}}
.pipeline-summary-grid > .ps-col-title{{font-size:16px;font-weight:600;
  padding:14px 16px 8px;background:var(--bg-secondary);
  border:0.5px solid var(--border-secondary);
  border-bottom:0.5px solid var(--border-secondary);
  border-radius:var(--radius-md) var(--radius-md) 0 0;
  margin:0;}}
.pipeline-summary-grid > .ps-cell-v2.ps-col-title{{
  color:var(--accent-v2);border-left:3px solid var(--accent-v2);}}
.pipeline-summary-grid > .ps-cell-prod.ps-col-title{{
  color:var(--accent-prod);border-left:3px solid var(--accent-prod);}}
.pipeline-summary-grid > .ps-cell{{padding:0 16px 14px;
  background:var(--bg-secondary);
  border-left:0.5px solid var(--border-secondary);
  border-right:0.5px solid var(--border-secondary);
  align-self:stretch;}}
.pipeline-summary-grid > .ps-cell-v2.ps-cell{{
  border-left:3px solid var(--accent-v2);}}
.pipeline-summary-grid > .ps-cell-prod.ps-cell{{
  border-left:3px solid var(--accent-prod);}}
/* last-row cells get bottom border + radius. Грид заполняется построчно;
   мы знаем, что 'secondary' — финальная секция, поэтому таргетим её. */
.pipeline-summary-grid > .ps-cell[data-section="secondary"]{{
  border-bottom:0.5px solid var(--border-secondary);
  border-radius:0 0 var(--radius-md) var(--radius-md);
  padding-bottom:14px;}}
.pipeline-summary-grid .ps-section-header{{font-size:13px;
  margin:14px 0 6px;color:var(--text-secondary);font-weight:500;
  text-transform:uppercase;letter-spacing:.5px;text-align:left;}}
.pipeline-summary-grid .ps-section-content{{}}
.box{{border:0.5px solid var(--border-secondary);border-radius:var(--radius-md);
      padding:14px;background:var(--bg-secondary);}}
.peaks{{border-collapse:collapse;width:100%;font-size:12px;}}
.peaks th,.peaks td{{border-bottom:0.5px solid var(--border-secondary);padding:4px 8px;text-align:left;}}
.peaks thead{{background:var(--bg);}}
.comp{{border-collapse:collapse;width:100%;font-size:11.5px;margin-top:6px;}}
.comp th,.comp td{{border-bottom:0.5px solid var(--border-secondary);padding:3px 6px;text-align:left;}}
.comp thead{{font-size:11px;color:var(--text-secondary);}}
.cluster{{border:0.5px solid var(--border-secondary);border-radius:var(--radius-md);
          padding:10px;margin-bottom:10px;background:var(--bg);}}
.cluster-head{{font-size:12.5px;margin-bottom:6px;}}
.badge{{display:inline-block;padding:2px 8px;border-radius:10px;color:#fff;
        font-size:10px;font-weight:600;letter-spacing:.5px;margin-right:6px;}}
.chip{{display:inline-block;padding:2px 8px;border-radius:10px;color:#fff;
       font-size:11px;font-weight:500;margin:0 4px 4px 0;}}
.chips-row{{margin:0 0 4px;}}
.dim{{color:var(--text-secondary);}}
/* F-391 / v1.18.27 — phantom anchor styling: library evidence без fit */
.comp tr.phantom td{{opacity:0.55;font-style:italic;}}
.ok{{color:#3a7;font-size:11px;margin-left:8px;}}
.bad{{color:#a23;font-size:11px;margin-left:8px;}}
.kpi{{display:flex;gap:18px;flex-wrap:wrap;margin:18px 0;}}
.kpi-card{{flex:1;min-width:180px;border:0.5px solid var(--border-secondary);
           border-radius:var(--radius-md);padding:12px 14px;background:var(--bg-secondary);}}
.kpi-label{{font-size:11px;color:var(--text-secondary);text-transform:uppercase;letter-spacing:.5px;}}
.kpi-value{{font-size:22px;font-weight:600;margin-top:4px;}}
.kpi-diff{{font-size:11px;color:var(--text-secondary);margin-top:2px;}}
.footer{{margin-top:36px;padding-top:14px;border-top:0.5px solid var(--border-secondary);
        font-size:11px;color:var(--text-secondary);}}
code{{background:rgba(0,0,0,.06);padding:1px 6px;border-radius:3px;font-size:11px;}}
@media (prefers-color-scheme:dark) {{
  code{{background:rgba(255,255,255,.1);}}
}}
@media (max-width:900px) {{
  .grid{{grid-template-columns:1fr;}}
  .pipeline-summary-grid{{grid-template-columns:1fr;}}
}}
/* F-393 / v1.18.27 — sortable table headers */
table[data-sortable] th[data-sortable]{{cursor:pointer;user-select:none;white-space:nowrap;}}
table[data-sortable] th[data-sortable]:hover{{color:var(--text);}}
table[data-sortable] th[data-sortable].sorted::after{{content:" ▾";color:var(--text-secondary);}}
table[data-sortable] th[data-sortable].sorted.asc::after{{content:" ▴";}}
/* F-RPT-02 / v1.18.29 — back-to-bundle nav link */
.back-nav{{display:block;margin:0 0 14px;font-size:12.5px;}}
.back-nav a{{display:inline-flex;align-items:center;gap:6px;
  padding:5px 10px;border:0.5px solid var(--border-secondary);
  border-radius:var(--radius-md);text-decoration:none;color:var(--text-secondary);
  background:transparent;transition:background .12s ease,color .12s ease;}}
.back-nav a:hover{{background:var(--bg-secondary);color:var(--text);}}
</style></head>
<body><div class="page">

<!-- F-UX-03 / 2026-06-04 — back-link убрана (supersedes F-RPT-02) -->

<h1>Сравнение конвейеров Production и Experimental V2 — Th-232 (F-390)</h1>
<p class="sub">
SpectraVibe {skill_v} · Маринелли 1 л · Gamma-1S NaI 63×63 ·
параллельный запуск <code>analyze_lsrm_spe</code> и V2 monkey-patched конвейера на одном спектре
</p>

<div class="meta">
<table>
<tr><th>Спектр пробы</th><td>{_esc(sample_name)}</td></tr>
<tr><th>Спектр фона</th><td>{_esc(bg_name)}</td></tr>
<tr><th>Масса пробы, кг</th><td>{mass_kg:.3f}</td></tr>
<tr><th>Версия скилла</th><td>{_esc(skill_v)}</td></tr>
<tr><th>Допуск сопоставления пиков</th><td>±4.0 кэВ</td></tr>
</table>
</div>

<div class="kpi">
  <div class="kpi-card">
    <div class="kpi-label">Пики ППП — V2</div>
    <div class="kpi-value">{len(v2_peaks)}</div>
    <div class="kpi-diff">двойной поиск (Мариcкотти ∪ согласованный фильтр)</div>
  </div>
  <div class="kpi-card">
    <div class="kpi-label">Пики ППП — production</div>
    <div class="kpi-value">{len(prod_peaks)}</div>
    <div class="kpi-diff">только Мариcкотти</div>
  </div>
  <div class="kpi-card">
    <div class="kpi-label">Совпадающие пики</div>
    <div class="kpi-value">{len(matched)}</div>
    <div class="kpi-diff">|ΔE| ≤ 4.0 кэВ</div>
  </div>
  <div class="kpi-card">
    <div class="kpi-label">Только в V2</div>
    <div class="kpi-value">+{len(only_v2)}</div>
    <div class="kpi-diff">не нашлось в production</div>
  </div>
  <div class="kpi-card">
    <div class="kpi-label">Только в production</div>
    <div class="kpi-value">+{len(only_prod)}</div>
    <div class="kpi-diff">не нашлось в V2</div>
  </div>
  <div class="kpi-card">
    <div class="kpi-label">Мультиплеты V2/P</div>
    <div class="kpi-value">{v2_mult}/{prod_mult}</div>
    <div class="kpi-diff">разложенные кластеры</div>
  </div>
</div>

<h2>Сводка по конвейерам</h2>
<!-- BUG-12 / v1.18.30 — pipeline-summary-grid: each section is one grid
     row (V2 cell + Prod cell), so headers align horizontally across columns
     regardless of column content length. -->
<div class="pipeline-summary-grid">
  {ps_title_v2}
  {ps_title_prod}
  {ps_v2_nuclides}
  {ps_prod_nuclides}
  {ps_v2_peaks}
  {ps_prod_peaks}
  {ps_v2_mult}
  {ps_prod_mult}
  {ps_v2_sec}
  {ps_prod_sec}
</div>

<h2>Различия пиков</h2>
<h3>Совпадающие пики (|ΔE| ≤ 4 кэВ) — {len(matched)} шт</h3>
<table class="peaks" data-sortable="true">
<thead><tr>
<th data-sort="num" data-sortable="true" onclick="sortTable(this)" data-default-sort="asc">E_prod, кэВ</th>
<th data-sort="num" data-sortable="true" onclick="sortTable(this)">E_V2, кэВ</th>
<th data-sort="num" data-sortable="true" onclick="sortTable(this)">ΔE, кэВ</th>
<th data-sort="num" data-sortable="true" onclick="sortTable(this)">σ_prod</th>
<th data-sort="num" data-sortable="true" onclick="sortTable(this)">σ_V2</th>
</tr></thead>
<tbody>
{matched_rows or '<tr><td colspan="5" class="dim">нет совпадений</td></tr>'}
</tbody>
</table>

<div class="grid" style="margin-top:18px;">
  <div class="box">
    <h3>Только в V2 ({len(only_v2)})</h3>
    <table class="peaks" data-sortable="true">
    <thead><tr>
    <th data-sort="num" data-sortable="true" onclick="sortTable(this)">канал</th>
    <th data-sort="num" data-sortable="true" onclick="sortTable(this)" data-default-sort="asc">E, кэВ</th>
    <th data-sort="str" data-sortable="true" onclick="sortTable(this)">нуклид</th>
    <th data-sort="num" data-sortable="true" onclick="sortTable(this)">E_lib</th>
    <th data-sort="num" data-sortable="true" onclick="sortTable(this)">I, %</th>
    <th data-sort="num" data-sortable="true" onclick="sortTable(this)">S/σ</th>
    </tr></thead>
    <tbody>{only_v2_html}</tbody>
    </table>
  </div>
  <div class="box">
    <h3>Только в production ({len(only_prod)})</h3>
    <table class="peaks" data-sortable="true">
    <thead><tr>
    <th data-sort="num" data-sortable="true" onclick="sortTable(this)">канал</th>
    <th data-sort="num" data-sortable="true" onclick="sortTable(this)" data-default-sort="asc">E, кэВ</th>
    <th data-sort="str" data-sortable="true" onclick="sortTable(this)">нуклид</th>
    <th data-sort="num" data-sortable="true" onclick="sortTable(this)">E_lib</th>
    <th data-sort="num" data-sortable="true" onclick="sortTable(this)">I, %</th>
    <th data-sort="num" data-sortable="true" onclick="sortTable(this)">S/σ</th>
    </tr></thead>
    <tbody>{only_prod_html}</tbody>
    </table>
  </div>
</div>

<h3>Нуклиды-различия</h3>
<div class="grid">
  <div class="box">
    <div style="font-size:12px;color:var(--text-secondary);margin-bottom:6px;">
      Только в V2 ({len(only_v2_nuc)})
    </div>
    {_chips_html(only_v2_nuc, '#3a7')}
  </div>
  <div class="box">
    <div style="font-size:12px;color:var(--text-secondary);margin-bottom:6px;">
      Только в production ({len(only_prod_nuc)})
    </div>
    {_chips_html(only_prod_nuc, '#37a')}
  </div>
</div>

<h2>Метрики сравнения</h2>
<table class="peaks">
<thead><tr><th>Метрика</th><th>V2</th><th>Production</th><th>Δ (V2 − P)</th></tr></thead>
<tbody>
<tr><td>Пиков полного поглощения</td><td>{len(v2_peaks)}</td><td>{len(prod_peaks)}</td><td>{len(v2_peaks)-len(prod_peaks):+d}</td></tr>
<tr><td>Идентифицированных нуклидов</td><td>{len(v2_nuc)}</td><td>{len(prod_nuc)}</td><td>{len(v2_nuc)-len(prod_nuc):+d}</td></tr>
<tr><td>Кластеров мультиплетов</td><td>{v2_mult}</td><td>{prod_mult}</td><td>{v2_mult-prod_mult:+d}</td></tr>
<tr><td>Вторичных пиков (комптон / пик вылета SE/DE)</td><td>{v2_sec}</td><td>{prod_sec}</td><td>{v2_sec-prod_sec:+d}</td></tr>
<tr><td>Неидентифицированных пиков</td><td>{v2_un}</td><td>{prod_un}</td><td>{v2_un-prod_un:+d}</td></tr>
</tbody>
</table>

<div class="footer">
SpectraVibe {skill_v} · F-354 peak_pipeline_v2 · F-367 V2 production-pipeline integration ·
F-390 двухколоночный лейаут · F-393 sortable peak tables · F-395 RU локализация ·
масса пробы {mass_kg:.3f} кг ·
сгенерировано <code>scripts/gen_v2_compare_th232.py</code>
</div>

</div>
<script>
// F-393 / v1.18.27 — vanilla sort util для peak tables.
// Каждый <th data-sort="num|str"> переключает порядок строк родительского
// <table>. Числовое значение парсится из cell.textContent (regex), либо
// из data-sort атрибута cell (если задан). Stable for strings (locale ru).
function sortTable(th) {{
  var table = th.closest('table');
  if (!table) return;
  var tbody = table.tBodies[0];
  if (!tbody) return;
  var ths = table.tHead ? table.tHead.querySelectorAll('th') : [];
  var colIdx = -1;
  for (var i = 0; i < ths.length; i++) {{ if (ths[i] === th) {{ colIdx = i; break; }} }}
  if (colIdx < 0) return;
  var type = th.getAttribute('data-sort') || 'str';
  var curDir = th.getAttribute('data-sort-dir');
  var asc = (curDir === 'asc') ? false : true;  // toggle
  var rows = Array.prototype.slice.call(tbody.rows);
  rows.sort(function (a, b) {{
    var ca = a.cells[colIdx], cb = b.cells[colIdx];
    if (!ca || !cb) return 0;
    var va = (ca.getAttribute('data-sort') !== null) ? ca.getAttribute('data-sort') : ca.textContent;
    var vb = (cb.getAttribute('data-sort') !== null) ? cb.getAttribute('data-sort') : cb.textContent;
    va = (va || '').trim(); vb = (vb || '').trim();
    var cmp;
    if (type === 'num') {{
      var na = parseFloat(va.replace(/[^0-9.\\-]/g, ''));
      var nb = parseFloat(vb.replace(/[^0-9.\\-]/g, ''));
      var aok = !isNaN(na), bok = !isNaN(nb);
      if (aok && bok) cmp = na - nb;
      else if (aok) cmp = -1;
      else if (bok) cmp = 1;
      else cmp = va.localeCompare(vb, 'ru');
    }} else {{
      cmp = va.localeCompare(vb, 'ru');
    }}
    return asc ? cmp : -cmp;
  }});
  // re-attach in order
  rows.forEach(function (r) {{ tbody.appendChild(r); }});
  // update visual indicators on all th of this header row
  for (var j = 0; j < ths.length; j++) {{
    ths[j].classList.remove('sorted', 'asc');
    ths[j].removeAttribute('data-sort-dir');
  }}
  th.classList.add('sorted');
  if (asc) th.classList.add('asc');
  th.setAttribute('data-sort-dir', asc ? 'asc' : 'desc');
}}
// Initial visual markers — set .sorted.asc on every th[data-default-sort="asc"].
(function () {{
  var ths = document.querySelectorAll('th[data-default-sort="asc"]');
  for (var i = 0; i < ths.length; i++) {{
    ths[i].classList.add('sorted', 'asc');
    ths[i].setAttribute('data-sort-dir', 'asc');
  }}
}})();
</script>
</body></html>
"""
    with (out_dir / "v2_compare_report.html").open("w", encoding="utf-8") as f:
        f.write(html)

    print(f"OK compare_data.json + v2_compare_report.html → {out_dir}")
    print(f"   production: peaks={len(prod_peaks)} nuc={len(prod_nuc)} mult={prod_mult} sec={prod_sec}")
    print(f"   V2:         peaks={len(v2_peaks)} nuc={len(v2_nuc)} mult={v2_mult} sec={v2_sec}")
    print(f"   matched={len(matched)}  only_prod={len(only_prod)}  only_v2={len(only_v2)}")
    print(f"   nuc only_v2={only_v2_nuc}  nuc only_prod={only_prod_nuc}")


if __name__ == "__main__":
    main()
