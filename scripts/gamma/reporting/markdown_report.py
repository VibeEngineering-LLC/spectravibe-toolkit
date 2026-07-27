"""
Markdown-отчёт (v1.17.4) — полностью на русском (D-04, F-108 RU glossary).

Генерируется по запросу пользователя; основной артефакт — JSON.
13 обязательных разделов спецификации сохранены, но названия и
табличные заголовки переведены на русский.
"""
from __future__ import annotations

import os
from typing import List, Optional

from gamma.reporting._kind_ru import FEATURE_KIND_RU as _FEATURE_KIND_RU


_VERSION_HISTORY = [
    ("2026-05", "v1.17.4", "F-115 анонимизация + русификация + F-110 композит + F-111 цепочка"),
    ("2026-05", "v1.17.3", "F-114 канонический интерактивный отчёт"),
    ("2026-05", "v1.17.0", "F-100 метод шаблонов"),
    ("2026-05", "v1.16.2", "F-98 квази-шаблонный подгонщик"),
    ("2026-05", "v1.16.1", "F-95 формула мёртвого времени, F-96 фон как априорная информация"),
    ("2026-05", "v1.16.0", "F-90..F-93 пик-образ, σ, CI, верхний предел"),
    ("2026-05", "v1.15.2", "F-89 привязка по имени файла, статус фона"),
    ("2026-05", "v1.15.1", "F-88 экспресс-порядок и цепочное доминирование"),
    ("2026-05", "v1.15.0", "F-86 графики, HTML, CLI"),
    ("2026-05", "v1.14.0", "Шаг 11 модуль отчётности"),
    ("2026-05", "v1.13.0", "Раунд 5: активности и MDA"),
]


_BACKGROUND_STATUS_LABEL = {
    "subtracted_from_external_file": "Фон ВЫЧТЕН (внешний файл)",
    # F-135 / v1.17.7 — авто-найденный + вычтенный фон
    "auto_resolved_from_directory":
        "Фон ВЫЧТЕН (авто-подбор F-131: тот же детектор, совместимая "
        "геометрия, |Δt| ≤ 90 дн.)",
    "embedded_present_not_subtracted":
        "Фон НЕ вычтен (встроенный фон присутствует, но не использован)",
    "absent_no_subtraction":
        "Фон НЕ вычтен (фоновый спектр отсутствует — cps включают "
        "природный вклад)",
}


def _rel_image_path(image_path: str, md_dir: Optional[str]) -> str:
    if md_dir is None:
        return str(image_path).replace("\\", "/")
    try:
        rel = os.path.relpath(image_path, md_dir)
    except ValueError:
        return str(image_path).replace("\\", "/")
    return rel.replace("\\", "/")


def _fmt_or_dash(x, fmt: str = "{:.2f}") -> str:
    if x is None:
        return "—"
    try:
        return fmt.format(x)
    except (TypeError, ValueError):
        return str(x)


def _h(level: int, text: str) -> str:
    return f"\n{'#' * level} {text}\n"


def _render_header(h: dict) -> str:
    out = [_h(2, "1. Заголовок")]
    bg_status = h.get("background_status", "")
    bg_label = _BACKGROUND_STATUS_LABEL.get(
        bg_status, f"({bg_status})" if bg_status else "—"
    )
    isotope_hints = h.get("filename_isotope_hints") or []
    isotope_hint_str = ", ".join(isotope_hints) if isotope_hints else "—"
    bg_fn = h.get("background_filename") or ""
    # BUG-39 / BUG-40 — visible detector-fallback row when applicable.
    fb = h.get("detector_fallback") or {}
    fb_reason = fb.get("reason") if isinstance(fb, dict) else ""
    if fb_reason and fb_reason != "profile_loaded_no_fallback":
        # Map machine-readable reason → Russian label (avoid EN-leak per F-386).
        # F2-A (2026-06-21): only "profile_not_on_disk" remains as a real
        # fallback reason; "efficiency_tbd_using_fallback_profile" was
        # retired together with the bogus Gamma-1S stub profile.
        _REASON_RU = {
            "profile_not_on_disk": "профиль не найден",
        }
        reason_ru = _REASON_RU.get(fb_reason, "подмена")
        fb_label = (
            f"⚠ {fb.get('requested', '?')} → {fb.get('actual', '?')} "
            f"({reason_ru})"
        )
    else:
        fb_label = "—"
    rows = [
        ("Имя файла образца",               h.get("sample_filename") or h.get("filename", "—")),
        # F-144 / v1.17.7 — имя фон-файла обязательно в отчёте
        ("Имя файла фона",                  bg_fn if bg_fn else "—"),
        ("Идентификатор образца",           h.get("sample_id", "—") or "—"),
        ("Оператор",                        h.get("operator") or "—"),
        ("Дата и время начала",             h.get("start_datetime") or "—"),
        ("t_живое, с",                      _fmt_or_dash(h.get("live_time_s"), "{:.1f}")),
        ("t_реальное, с",                   _fmt_or_dash(h.get("real_time_s"), "{:.1f}")),
        ("Мёртвое время, %",                _fmt_or_dash(h.get("dead_time_pct"), "{:.2f}")),
        ("Детектор",                        h.get("detector_canonical") or h.get("detector_type") or "—"),
        ("Подмена профиля детектора",       fb_label),
        ("Идентификатор детектора",         h.get("detector_id") or "—"),
        ("Геометрия",                       h.get("geometry") or "—"),
        ("Геометрия (канон.)",              _ru_cell(h.get("geometry_canonical") or "—")),
        ("Окружение",                       _ru_cell(h.get("environment", "—"))),
        ("**Фон**",                         f"**{bg_label}**"),
        ("Подсказка изотопа из имени",      isotope_hint_str),
        ("Режим анализа",                   _ru_cell(h.get("analysis_mode") or "—")),
        ("Верхняя граница, кэВ",            _fmt_or_dash(h.get("energy_ceiling_keV"), "{:.0f}")),
        ("Максимум использован, кэВ",       _fmt_or_dash(h.get("energy_max_keV_kept"), "{:.1f}")),
        ("Каналы (использовано)",           h.get("n_channels", "—")),
        ("Каналы (всего)",                  h.get("n_channels_raw", "—")),
        ("Отброшено (переполнение)",        h.get("dropped_high_energy_count", "—")),
    ]
    out.append("| Поле | Значение |")
    out.append("|---|---|")
    for k, v in rows:
        out.append(f"| {k} | {v} |")
    return "\n".join(out)


def _render_detector(h: dict, calib: dict) -> str:
    out = [_h(2, "2. Тип детектора")]
    fwhm = calib.get("fwhm_cal", {}) or {}
    fwhm661 = fwhm.get("fwhm_at_661_keV")
    coefs = fwhm.get("coefficients", [])
    coef_str = ", ".join(_fmt_or_dash(c, "{:.4g}") for c in coefs[:3])
    out.append(
        f"Канонический детектор: **{h.get('detector_canonical') or h.get('detector_type') or '—'}**.\n"
    )
    out.append(
        f"Модель ПШПВ(E): `{fwhm.get('model', '?')}` с коэффициентами `({coef_str})`. "
        f"ПШПВ при 661.66 кэВ ≈ **{_fmt_or_dash(fwhm661, '{:.1f}')} кэВ** "
        f"(источник: {_ru_cell(fwhm.get('source', '?'))})."
    )
    return "\n".join(out)


def _render_calibration(calib: dict) -> str:
    out = [_h(2, "3. Калибровка")]
    e = calib.get("energy_cal", {}) or {}
    e_coefs = e.get("coefficients", [])
    e_coef_str = ", ".join(_fmt_or_dash(c, "{:.6g}") for c in e_coefs)
    out.append(
        f"**Энергетическая калибровка**: полином степени {e.get('degree','?')}, "
        f"коэффициенты `({e_coef_str})`, источник: `{e.get('source','?')}`."
    )
    slc = calib.get("seven_line_check") or None
    if slc is not None:
        out.append("")
        out.append("**Проверка калибровки по 7 линиям ЕРН (Lsrm §9)**:")
        total = slc.get("lines_total", 7)
        present = slc.get("lines_present", 0)
        out.append(f"- проверено опорных линий: {total}")
        out.append(f"- найдено совпадений: {present}/{total}")
        out.append(f"- максимальная невязка: {_fmt_or_dash(slc.get('max_residual_keV'), '{:.2f}')} кэВ")
        max_frac = slc.get("max_residual_fwhm_fraction")
        if max_frac is not None:
            out.append(f"- max|Δ|/ПШПВ: {_fmt_or_dash(max_frac, '{:.0%}')}")
        cal_q = slc.get("quality") or "—"
        out.append(f"- качество калибровки: **{_ru_cell(cal_q)}**")
        note = slc.get("quality_note")
        if note:
            out.append(f"- пояснение: {note}")
    return "\n".join(out)


_CELL_TRANSLATIONS = {
    "True": "да",
    "False": "нет",
    "cowell": "Кауэлл",
    "stored": "сохранён",
    "failed": "не удалось",
    "default_NaI_63x63": "модель NaI 63×63",
    "natural": "природное",
    "low_background": "низкий фон",
    "background_only": "только фон",
    "unknown": "не определено",
    "sample_anchor_rank": "ранжирование по якорям",
    "background_7line": "фон по 7 линиям ЕРН",
    "marinelli_1L": "Маринелли 1 л",
    "subtracted_from_external_file": "вычтен из внешнего файла",
    "auto_resolved_from_directory": "вычтен (авто-подбор F-131)",
    "absent_no_subtraction": "отсутствует — без вычитания",
    "embedded_present_not_subtracted": "встроен, но не вычтен",
    "weighted_mean": "взвешенное среднее",
    "spread": "разброс",
    # TD-4: feature_kind translations imported from shared _kind_ru module.
    # Do NOT add new feature_kind entries here — extend _kind_ru.FEATURE_KIND_RU.
    **_FEATURE_KIND_RU,
    "IC_Cs-137": "ВК Cs-137",
    "cluster": "композит",
    "true_unmatched": "истинно неотождест.",
    "xrf": "рентген. флуор.",
    "edge_of_range": "край диапазона",
    "anchor_seeded": "посеяно по якорю",
    "auto_resolve_overlap": "авто-разрешение наложения",
    # F-117 / F-118 (v1.17.5)
    "scatter": "разброс",
    "deconvolved": "по деконволюции",
    "deconvolved_coupled": "по связанной деконволюции",
    "cowell": "Cowell",
    "lsrm_peaks_table": "табл. пиков LSRM",
    "failed": "не подгоняется",
    "complete": "полно",
    "moderate": "умеренный",
    "high": "высокий",
    "low": "низкий",
    "very_low": "очень низкий",
    "very_high": "очень высокий",
    "good": "хорошее",
    "poor": "слабое",
    "ok": "ok",
    "drift": "дрейф",
    "broken": "не валидна",
    "n/a": "—",
    "marginal": "пограничное",
    "noise": "шум",
    "confirmed": "подтверждён",
    "tentative": "предположительный",
    "detected": "обнаружен",
}


def _ru_text(s) -> str:
    """Translate free-form text with the cell map + ru_filter.

    Used for individual cells that may be embedded English phrases.
    """
    if s is None:
        return "—"
    text = str(s)
    for en, ru in _CELL_TRANSLATIONS.items():
        text = text.replace(en, ru)
    return text


def _ru_cell(v) -> str:
    if isinstance(v, bool):
        return "да" if v else "нет"
    s = str(v)
    return _CELL_TRANSLATIONS.get(s, s)


def _render_table(title: str, rows: List[dict], columns: List[tuple]) -> str:
    out = [_h(2, title)]
    if not rows:
        out.append("_(нет данных)_")
        return "\n".join(out)
    out.append("| " + " | ".join(c[0] for c in columns) + " |")
    out.append("|" + "|".join("---" for _ in columns) + "|")
    for r in rows:
        cells = []
        for _, key, fmt in columns:
            v = r.get(key)
            if v is None:
                cells.append("—")
            elif fmt is None:
                cells.append(_ru_cell(v))
            else:
                try:
                    cells.append(fmt.format(v))
                except (TypeError, ValueError):
                    cells.append(_ru_cell(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def _render_identified_nuclides(rows: List[dict]) -> str:
    out = [_h(2, "7. Идентифицированные нуклиды")]
    if not rows:
        out.append("_(нет нуклидов, подтверждённых Stage-1 + disambiguate)_")
        return "\n".join(out)
    for r in rows:
        out.append(f"\n### {r.get('nuclide','?')}  "
                   f"(уровень: {_ru_cell(r.get('tier','?'))})")
        out.append(f"- Совпавших линий: {r.get('n_matched_lines','?')}")
        out.append(f"- Характеристическая линия: {_fmt_or_dash(r.get('characteristic_line_keV'), '{:.2f}')} кэВ")
        out.append(f"- Индекс достоверности (CI): {_fmt_or_dash(r.get('confidence_index'), '{:.2f}')}  "
                   f"({_ru_cell(r.get('confidence_level','?'))})")
        out.append(f"- Скорость счёта: {_fmt_or_dash(r.get('peak_rate_cps'), '{:.3g}')} cps")
        if r.get("is_upper_limit"):
            ul = r.get("upper_limit_Bq")
            if ul is not None:
                out.append(f"- Активность: **≤ {_fmt_or_dash(ul, '{:.3g}')} Бк (верхний предел, L_D)** — "
                           f"σ/A = {_fmt_or_dash(r.get('activity_relative_sigma'), '{:.0%}')} > 50% "
                           f"по LSRM §11")
            else:
                out.append(f"- Активность: **верхний предел** — "
                           f"σ/A = {_fmt_or_dash(r.get('activity_relative_sigma'), '{:.0%}')} > 50%, "
                           f"L_D недоступен")
        else:
            out.append(f"- Активность: {_fmt_or_dash(r.get('activity_Bq'), '{:.3g}')} ± "
                       f"{_fmt_or_dash(r.get('activity_sigma_Bq'), '{:.2g}')} Бк "
                       f"(σ по методу {_ru_cell(r.get('activity_sigma_method', 'weighted_mean'))})")
        if r.get("specific_activity_Bq_per_kg") is not None:
            out.append(f"- Удельная активность: {_fmt_or_dash(r.get('specific_activity_Bq_per_kg'), '{:.3g}')} ± "
                       f"{_fmt_or_dash(r.get('specific_activity_sigma_Bq_per_kg'), '{:.2g}')} Бк/кг")
        out.append(f"- Поправка на распад: {'да' if r.get('decay_corrected') else 'нет'} "
                   f"(коэффициент {_fmt_or_dash(r.get('decay_factor'), '{:.4f}')})")
        if r.get("cascade_warning"):
            out.append(f"- Предупреждение о каскадном суммировании: "
                       f"{_ru_filter(r['cascade_warning']) or '—'}")
        reason_raw = r.get('reason', '?') or "?"
        # Try aggressive translation; if any 4+ letter ASCII word
        # remains, fall back to a safe stub.
        reason_ru = _ru_filter(reason_raw)
        if not reason_ru:
            reason_ru = "См. диагностику и таблицу характеристических линий."
        out.append(f"- Обоснование: {reason_ru}")
    return "\n".join(out)


def _render_priority_express(json_dict: dict) -> str:
    out = [_h(2, "3α. Приоритетные экспресс-опорные линии (методика пользователя)")]
    findings = json_dict.get("priority_express_findings", []) or []
    diag = json_dict.get("diagnostics", {}) or {}
    cd = diag.get("chain_dominance") or {}
    k40_warn = bool(diag.get("k40_ac228_overlap_warning"))

    if not findings:
        out.append("_(приоритетные результаты не заполнены — отчёт до v1.15.1)_")
        return "\n".join(out)

    out.append("| # | Сигнал | Статус | σ | Примечание |")
    out.append("|---|---|---|---|---|")
    for pf in findings:
        marker = "✔ ПОДТВЕРЖДЕНО" if pf.get("matched") else "✘ отсутствует"
        sig = pf.get("max_significance_sigma")
        sig_str = _fmt_or_dash(sig, "{:.1f}") if sig else "—"
        label_raw = pf.get('label', '?') or "?"
        label_ru = (label_raw
                    .replace("keV", "кэВ")
                    .replace("chain", "цепочка")
                    .replace("(pair)", "(пара)")
                    .replace(" pair", " пара")
                    .replace("series", "ряд")
                    .replace("Ra-pair", "Ra-пара"))
        note_raw = pf.get('note', '') or ""
        note_ru = (note_raw
                   .replace("MATCHED", "ПОДТВЕРЖДЕНО")
                   .replace("matched", "подтверждено")
                   .replace("missing", "отсутствует")
                   .replace("lines", "линий"))
        out.append(
            f"| {pf.get('order','?')} "
            f"| {label_ru} "
            f"| {marker} "
            f"| {sig_str} "
            f"| {note_ru} |"
        )

    out.append("")
    out.append("### Доминирование цепочки")
    if cd:
        th_str = "**ДА**" if cd.get("th232_dominant") else "нет"
        u_str = "**ДА**" if cd.get("u238_dominant") else "нет"
        out.append(
            f"- **Цепочка Th-232 доминирует**: {th_str}  "
            f"(σ ≤ {_fmt_or_dash(cd.get('th232_strength_sigma'), '{:.1f}')})"
        )
        if cd.get("th232_evidence"):
            out.append("  - Подтверждение:")
            for e in cd["th232_evidence"]:
                out.append(f"    - {e}")
        out.append(
            f"- **Цепочка U-238 доминирует**: {u_str}  "
            f"(σ ≤ {_fmt_or_dash(cd.get('u238_strength_sigma'), '{:.1f}')})"
        )
        if cd.get("u238_evidence"):
            out.append("  - Подтверждение:")
            for e in cd["u238_evidence"]:
                out.append(f"    - {e}")
        if cd.get("reason"):
            r = _ru_filter(cd['reason']) or cd['reason']
            out.append(f"- Обоснование: _{r}_")
    else:
        out.append("_(доминирование цепочки не вычислено)_")

    if k40_warn:
        out.append("")
        out.append("> ПРЕДУПРЕЖДЕНИЕ: наложение K-40 / Ac-228 (F-88)")
        out.append("> ")
        out.append("> Цепочка Th-232 доминирует И приоритетный сигнал K-40 1460.82 кэВ "
                   "подтверждён. На NaI 63×63 линия Ac-228 1459.20 кэВ (I=0.85%) "
                   "не разрешается от K-40 — площадь пика K-40 содержит вклад "
                   "Ac-228. Активность K-40 не следует приводить без разложения "
                   "относительно подтверждённой опорной линии Tl-208 2614.51 "
                   "или отдельной опоры Ac-228.")

    suppressed = cd.get("suppressed_chains") or []
    sup_reason = cd.get("suppression_reason") or ""
    dropped = cd.get("chain_filtered_out_nuclides") or []

    # F-388 / v1.18.26 — conditional render: блок подавления отображается
    # только если хотя бы один primary_fep принадлежит одной из подавляемых
    # цепочек. Иначе сообщение пользователю — шум: нет ни «pickup'ов»,
    # ни matched_lines для упоминаемой цепочки. Симметрично interactive_html
    # реализации.
    _U_CHAIN_MEMBERS = {
        "Pb-214", "Bi-214", "Pb-210", "Bi-210", "Po-214", "Po-218",
        "Ra-226", "Rn-222", "U-238", "U-234", "Th-234", "Pa-234",
    }
    _TH_CHAIN_MEMBERS = {
        "Tl-208", "Pb-212", "Ac-228", "Bi-212", "Th-228", "Th-232",
        "Ra-228", "Po-216", "Po-212",
    }

    if suppressed:
        feps = json_dict.get("primary_feps") or []
        present_nuclides = {fp.get("nuclide") or "" for fp in feps}

        def _chain_has_evidence(chain_name: str) -> bool:
            name = chain_name.strip()
            if "U-238" in name or "Ra-226" in name or "U_" in name:
                members = _U_CHAIN_MEMBERS
            elif "Th-232" in name or "Th_" in name:
                members = _TH_CHAIN_MEMBERS
            else:
                return True
            return bool(present_nuclides & members)

        relevant = [c for c in suppressed if _chain_has_evidence(c)]
        if relevant:
            out.append("")
            out.append(f"> Подавление цепочки по привязке имени файла (F-89d)")
            out.append("> ")
            out.append(f"> Подавлённые цепочки: **{', '.join(relevant)}**")
            if dropped:
                out.append(f"> ")
                out.append(f"> Нуклиды, исключённые из идентификации: "
                           f"**{', '.join(dropped)}**")
            if sup_reason:
                # BUG-6 / v1.18.30+ (Agent B): fallback переписан без bare
                # F-id (F-317 strip срезает F-89 и оставляет «правилу d.»).
                ru = _ru_filter(sup_reason) or (
                    "Цепочка подавлена согласно привязке имени файла "
                    "(filename binding): нет независимых линий, кроме хвостов "
                    "других цепочек."
                )
                out.append("> ")
                out.append("> " + ru)

    return "\n".join(out)


def _render_diagnostics(diag: dict, cmp: dict) -> str:
    out = [_h(2, "12. Диагностика")]
    rows = [
        ("Окружение измерения",            _ru_cell(diag.get("measurement_environment", "?"))),
        ("Мёртвое время, %",               _fmt_or_dash(diag.get("dead_time_pct"), "{:.2f}")),
        ("Поправка на мёртвое время",      "да" if diag.get("dead_time_correction_applied") else "нет"),
        ("Поправка TCS применена",         "да" if diag.get("tcs_correction_applied") else "нет"),
        ("Индикатор наложений",            "да" if diag.get("pile_up_indicator") else "нет"),
        ("Аннигиляция 511 обнаружена",     "да" if diag.get("annihilation_511_observed") else "нет"),
        # F-386.1 / v1.18.28 (Agent B) — терминология «пик вылета»,
        # не «ускользание» (закрепляет F-386 hardlock).
        ("Пиков вылета",                   diag.get("n_escape_peaks", 0)),
        ("Сумматорных пиков",              diag.get("n_sum_peaks", 0)),
        ("XRF-остатков",                   diag.get("n_xrf_residuals", 0)),
        ("Фон вычтен",                     "да" if diag.get("background_subtracted") else "нет"),
        ("Качество калибровки",            _ru_cell(diag.get("calibration_quality") or "—")),
        ("Полнота DC, %",                  _fmt_or_dash(diag.get("completeness_dc_pct"), "{:.2f}")),
        ("Флаг полноты",                   _ru_cell(diag.get("completeness_flag") or "—")),
        ("Каскадные нуклиды",              ", ".join(diag.get("cascade_warning_nuclides") or []) or "—"),
        ("Источник модели ПШПВ",           _ru_cell(diag.get("fwhm_model_source") or "—")),
        ("Источник кривой эффективности",  diag.get("efficiency_source") or "—"),
        ("Кривая эффективности загружена", "да" if diag.get("efficiency_loaded") else "нет"),
    ]
    out.append("| Поле | Значение |")
    out.append("|---|---|")
    for k, v in rows:
        out.append(f"| {k} | {v} |")

    intr = diag.get("intrinsic_activity_signature") or {}
    if intr:
        out.append("")
        out.append("### Собственная активность детектора")
        out.append(f"Детектор: **{intr.get('detector_canonical','?')}** — "
                   f"собственная активность Бк/см³: {_fmt_or_dash(intr.get('Bq_per_cm3'), '{:.4g}')}.")
        # Artefact / signature listings carry English text from the
        # JSON dataset; we render them as plain RU stubs to stay
        # F-108-compliant.  The detailed catalogue is available in
        # the JSON file for downstream consumers.
        if intr.get("expected_artefacts"):
            out.append("")
            out.append("Ожидаемые артефакты: см. JSON-отчёт (поле "
                       "`intrinsic_activity_signature.expected_artefacts`).")
        if intr.get("absent_signatures"):
            out.append("")
            out.append("Отсутствующие сигнатуры (проверка): см. JSON-отчёт.")
    return "\n".join(out)


def _render_version_history() -> str:
    out = [_h(2, "13. История версий")]
    out.append("| Дата | Версия | Изменения |")
    out.append("|---|---|---|")
    for date, ver, note in _VERSION_HISTORY:
        out.append(f"| {date} | {ver} | {note} |")
    return "\n".join(out)


def build_markdown_report(
    result_or_json,
    *,
    json_dict=None,
    plots: Optional[dict] = None,
    md_dir: Optional[str] = None,
) -> str:
    """Сборка Markdown-отчёта (полностью RU, v1.17.4).

    Принимает либо ``StagedAnalysisResult`` (внутри соберётся
    JSON), либо уже готовый ``json_dict``.
    """
    if json_dict is None:
        from gamma.reporting.json_report import build_json_report
        json_dict = build_json_report(result_or_json)

    h = json_dict.get("header", {}) or {}
    calib = json_dict.get("calibration", {}) or {}

    parts: List[str] = []
    parts.append(f"# Отчёт о гамма-спектрометрическом анализе\n")
    parts.append(f"_(сборка {json_dict.get('skill_version','?')}, "
                 f"схема {json_dict.get('schema_version','?')})_\n")

    parts.append(_render_header(h))
    parts.append(_render_detector(h, calib))
    parts.append(_render_calibration(calib))

    parts.append(_render_priority_express(json_dict))

    pf_rows = []
    for row in (json_dict.get("primary_feps", []) or []):
        row = dict(row)
        if row.get("is_upper_limit_artifact"):
            sigma_val = row.get("peak_area_sigma")
            if isinstance(sigma_val, (int, float)) and sigma_val > 0:
                row["peak_area_counts"] = f"< {sigma_val:.3g} (<MDA)"
            else:
                row["peak_area_counts"] = "<MDA"
            row["rate_cps"] = "—"
        pf_rows.append(row)
    parts.append(_render_table(
        "4. Основные пики полного поглощения",
        pf_rows,
        [
            ("Нуклид",            "nuclide",           None),
            ("E_библ., кэВ",      "library_E_keV",     "{:.2f}"),
            ("I_библ., %",        "library_I_pct",     "{:.3g}"),
            ("Канал",             "peak_channel",      None),
            ("E_пика, кэВ",       "peak_E_keV",        "{:.2f}"),
            ("ПШПВ, кэВ",         "fwhm_keV",          "{:.2f}"),
            ("Площадь, отсч.",    "peak_area_counts",  "{:.3g}"),
            ("σ_площ.",           "peak_area_sigma",   "{:.3g}"),
            ("Скор., cps",        "rate_cps",          "{:.3g}"),
            ("Источник",          "peak_area_source",  None),
            ("Характ.",           "is_characteristic", None),
        ],
    ))
    sp_rows = []
    for sp in (json_dict.get("secondary_peaks", []) or []):
        sp = dict(sp)
        # Translate the free-text `note` field; if translation leaves
        # any English ≥4 letters, replace with a generic stub.
        if sp.get("note"):
            note_ru = _ru_filter(sp["note"])
            if not note_ru:
                note_ru = _ru_text(sp["note"])
                # Final scrub: drop residual English ASCII words ≥4 letters
                import re as _re
                if _re.search(r"[A-Za-z]{4,}", note_ru):
                    # Strip residual English to keep field RU-only
                    note_ru = _re.sub(r"\b[A-Za-z]{4,}\b", "", note_ru).strip()
                    note_ru = _re.sub(r"\s+", " ", note_ru) or "—"
            sp["note"] = note_ru
        # Translate type / feature_kind via cell map
        for key in ("type", "feature_kind"):
            if sp.get(key):
                sp[key] = _ru_cell(sp[key])
        sp_rows.append(sp)
    parts.append(_render_table(
        "5. Вторичные пики",
        sp_rows,
        [
            ("E, кэВ",            "energy_keV",        "{:.2f}"),
            ("σ",                 "significance",      "{:.1f}"),
            ("Тип",               "type",              None),
            ("Признак",           "feature_kind",      None),
            ("Родитель",          "parent_nuclide",    None),
            ("E_родителя, кэВ",   "parent_line_keV",   "{:.2f}"),
            ("Примечание",        "note",              None),
        ],
    ))

    parts.append(_h(2, "6. Элементная XRF"))
    xrf = json_dict.get("elemental_xrf", []) or []
    if not xrf:
        parts.append("_(нет XRF-остатков классифицировано)_")
    else:
        for entry in xrf:
            parts.append(f"\n### {entry.get('element','?')}  "
                         f"({_ru_cell(entry.get('mechanism','?'))})")
            parts.append(f"- Линий: {entry.get('n_observed', 0)}")
            for ln in entry.get("observed_lines", []) or []:
                parts.append(
                    f"  - {_fmt_or_dash(ln.get('energy_keV'), '{:.2f}')} кэВ  "
                    f"(σ={_fmt_or_dash(ln.get('significance'), '{:.1f}')}, "
                    f"библ. {_fmt_or_dash(ln.get('library_E_keV'), '{:.2f}')}, "
                    f"Δ={_fmt_or_dash(ln.get('delta_keV'), '{:.2f}')})"
                )

    parts.append(_render_identified_nuclides(
        json_dict.get("identified_nuclides", []) or []
    ))

    parts.append(_render_table(
        "8. Неидентифицированные значимые пики",
        json_dict.get("unidentified_peaks", []) or [],
        [
            ("E, кэВ",       "energy_keV",   "{:.2f}"),
            ("σ",            "significance", "{:.1f}"),
            ("Метка",        "label",        None),
            ("Примечание",   "note",         None),
        ],
    ))
    cmp = json_dict.get("completeness", {}) or {}
    parts.append(
        f"\n**Полнота описания (DC)** = "
        f"{_fmt_or_dash(cmp.get('dc_pct'), '{:.2f}')} %  "
        f"[{_ru_cell(cmp.get('flag') or '—')}]\n"
    )

    parts.append(_h(2, "9. График спектра"))
    spectrum_path = (plots or {}).get("spectrum") if plots else None
    if spectrum_path:
        rel = _rel_image_path(spectrum_path, md_dir)
        parts.append(f"![График спектра]({rel})")
    else:
        parts.append(
            "_(Графики не сгенерированы — "
            "вызовите `build_report(..., write_plots=True)` для PNG.)_"
        )

    # BUG-5 / v1.18.30+ (Agent B): заголовок явно обозначает источник (образец).
    # Markdown-отчёт рендерит только sample-мультиплеты — bg-мультиплеты живут
    # только в интерактивном HTML (см. interactive_html.py::_build_multiplet_blocks).
    # BUG-14 / v1.18.30+ (Agent B): для чисто-фоновых спектров
    # (diagnostics.measurement_environment == "background_only") секция 10
    # «Мультиплеты — разложение в спектре образца» (и подсекция F-145
    # самокалибровки) пропускается полностью — пики не значимы над шумовым
    # континуумом, подгонка χ²/ν к нему вводит оператора в заблуждение.
    # Симметрично с interactive_html (BUG-14). multiplet_deconvolutions
    # остаётся доступным в JSON для аудита; просто не рендерится в MD.
    _diag_md = json_dict.get("diagnostics", {}) or {}
    _is_bg_only_md = (_diag_md.get("measurement_environment") == "background_only")
    if not _is_bg_only_md:
        parts.append(_h(2, "10. Мультиплеты — разложение в спектре образца"))
        parts.append(
            "_Первичная подгонка по библиотечным интенсивностям. Площади "
            "**S** взяты из спектра образца — основной источник данных для "
            "расчёта активности._"
        )
        decons = json_dict.get("multiplet_deconvolutions", []) or []
        multiplet_paths = (plots or {}).get("multiplets") if plots else None
        if not decons:
            parts.append("_(в этом прогоне мультиплеты не обрабатывались)_")
        else:
            for i, d in enumerate(decons, start=1):
                cluster_label = d.get("cluster_id") or f"M{i}"
                phase_A_chi2 = d.get("F145_phase_A_chi2_per_dof")
                f145_header = ""
                if phase_A_chi2 is not None:
                    f145_header = (
                        f"; F-145 фаза А χ²/ν = "
                        f"{_fmt_or_dash(phase_A_chi2, '{:.2f}')}"
                    )
                parts.append(f"\n### Кластер {cluster_label}  (χ²/ν = "
                             f"{_fmt_or_dash(d.get('chi2_per_dof'), '{:.3f}')}, "
                             f"сошлось: {'да' if d.get('converged') else 'нет'}"
                             f"{f145_header})")
                if multiplet_paths and i - 1 < len(multiplet_paths):
                    rel = _rel_image_path(multiplet_paths[i - 1], md_dir)
                    parts.append(f"![Кластер {cluster_label}]({rel})\n")
                # F-145: показать колонку «сдвиг dE» если есть данные Phase A
                has_shift = any(
                    c.get("F145_centroid_shift_keV") is not None
                    for c in d.get("components", []) or []
                )
                if has_shift:
                    parts.append("| Нуклид | E_библ., кэВ | dE_фаза_А, кэВ | Площадь | σ_площ. |")
                    parts.append("|---|---|---|---|---|")
                    for c in d.get("components", []) or []:
                        parts.append(
                            f"| {c.get('nuclide','?')} "
                            f"| {_fmt_or_dash(c.get('line_E_keV'), '{:.2f}')} "
                            f"| {_fmt_or_dash(c.get('F145_centroid_shift_keV'), '{:+.2f}')} "
                            f"| {_fmt_or_dash(c.get('deconvolved_area'), '{:.3g}')} "
                            f"| {_fmt_or_dash(c.get('deconvolved_area_sigma'), '{:.3g}')} |"
                        )
                else:
                    parts.append("| Нуклид | E_библ., кэВ | Площадь | σ_площ. |")
                    parts.append("|---|---|---|---|")
                    for c in d.get("components", []) or []:
                        parts.append(
                            f"| {c.get('nuclide','?')} "
                            f"| {_fmt_or_dash(c.get('line_E_keV'), '{:.2f}')} "
                            f"| {_fmt_or_dash(c.get('deconvolved_area'), '{:.3g}')} "
                            f"| {_fmt_or_dash(c.get('deconvolved_area_sigma'), '{:.3g}')} |"
                        )

            # F-145 / v1.17.8 — Раздел самокалибровки по итогам фаз А→D
            f145 = json_dict.get("multiplet_self_calibration") or {}
            if f145.get("attempted"):
                parts.append("\n### F-145: двухфазная самокалибровка по мультиплетам")
                if f145.get("phase_C_applied"):
                    parts.append(
                        f"- Шкала E(N) пересчитана по {f145.get('n_anchors_after_filter', 0)} "
                        f"опорным точкам (центроиды мультиплетов, степень полинома "
                        f"коррекции = {f145.get('degree_used')})."
                    )
                    old_r = f145.get("old_residual_max_keV")
                    new_r = f145.get("new_residual_max_keV")
                    if old_r is not None and new_r is not None:
                        parts.append(
                            f"- Максимальный остаточный сдвиг: "
                            f"{_fmt_or_dash(old_r, '{:.3f}')} → "
                            f"{_fmt_or_dash(new_r, '{:.3f}')} кэВ."
                        )
                    if f145.get("anchors_used"):
                        parts.append("- Опорные точки (взвешенно по I_γ внутри мультиплета):")
                        for a in f145["anchors_used"]:
                            # source: multiplet_<id>_I_pct_weighted → перевод
                            src_raw = str(a.get("source",""))
                            src_ru = src_raw.replace(
                                "multiplet_", "мультиплет "
                            ).replace(
                                "_I_pct_weighted", " (взвеш. по I_γ)"
                            )
                            parts.append(
                                f"  - {a.get('nuclide','?')} "
                                f"{_fmt_or_dash(a.get('E_passport_keV'), '{:.2f}')} кэВ → "
                                f"сдвиг {_fmt_or_dash(a.get('drift_keV'), '{:+.3f}')} кэВ "
                                f"(канал {_fmt_or_dash(a.get('channel_fitted'), '{:.2f}')}, "
                                f"источник: {src_ru})"
                            )
                else:
                    parts.append(
                        f"- Калибровка не пересчитана: {f145.get('reason','')}"
                    )

    parts.append(_render_table(
        "11. Таблица MDA (ISO 11929 / Lsrm §6.3)",
        json_dict.get("mda", []) or [],
        [
            ("Нуклид",       "nuclide",                None),
            ("E, кэВ",       "line_E_keV",             "{:.2f}"),
            ("I, %",         "intensity_pct",          "{:.3g}"),
            ("ε",            "efficiency",             "{:.3g}"),
            ("L_C, отсч.",   "decision_threshold_counts", "{:.3g}"),
            ("L_D, отсч.",   "detection_limit_counts", "{:.3g}"),
            ("MDA, Бк",      "MDA_Bq",                 "{:.3g}"),
        ],
    ))

    parts.append(_render_diagnostics(
        json_dict.get("diagnostics", {}) or {},
        cmp,
    ))

    warnings = json_dict.get("warnings", []) or []
    # BUG-39 / BUG-40 — when a detector_fallback warning is present in the
    # warnings list, swap the bilingual EN string for the RU-only variant
    # stored in `header.detector_fallback.human_ru`. This keeps the MD
    # report compliant with F-386 (no English ASCII words leak).
    fb_header = (json_dict.get("header", {}) or {}).get("detector_fallback") or {}
    fb_en = (fb_header.get("human_en") or "").strip() if isinstance(fb_header, dict) else ""
    fb_ru = (fb_header.get("human_ru") or "").strip() if isinstance(fb_header, dict) else ""
    if warnings:
        parts.append(_h(2, "Предупреждения"))
        for w in warnings:
            # BUG-40 — warnings may now be heterogeneous (str | dict). Dict
            # entries carry a machine-readable ``code`` field. Render a
            # hand-written RU equivalent for known codes; never emit the
            # English ``message`` field directly (F-386 EN-leak gate).
            if isinstance(w, dict):
                parts.append(f"- {_render_warning_dict_ru(w)}")
                continue
            if fb_en and fb_ru and fb_en in w:
                # Drop the EN portion; emit the RU message directly.
                parts.append(f"- {fb_ru}")
            else:
                parts.append(f"- {_translate_warning(w)}")

    # F-132 / v1.17.7 — обязательный раздел «Оценка стоимости анализа»
    # с поэтапной разбивкой + итогом «~N токенов или M% от бесплатной
    # 5-часовой сессии».
    parts.append(_render_cost_estimate(json_dict.get("cost_estimate")))

    parts.append(_render_version_history())

    return "\n".join(parts) + "\n"


def _render_cost_estimate(cost_block: Optional[dict]) -> str:
    """F-132 / v1.17.7 — раздел «Оценка стоимости анализа» (всегда).

    Структура:
      * Заголовок «Оценка стоимости анализа»
      * Таблица по 10 этапам Step 1..11 (base + complexity + итог + почему)
      * Строка-итог + % от 5-часовой сессии
    """
    out: List[str] = []
    out.append(_h(2, "Оценка стоимости анализа"))
    if not cost_block:
        out.append("_(оценка стоимости недоступна для этого прогона)_")
        return "\n".join(out)
    by_stage = cost_block.get("by_stage") or []
    total = int(cost_block.get("tokens_total") or 0)
    budget = int(cost_block.get("session_token_budget") or 0)
    pct = float(cost_block.get("session_pct") or 0.0)
    override = bool(cost_block.get("override_used"))
    detail = cost_block.get("detail") or ""

    out.append(
        "Поэтапная разбивка стоимости анализа в токенах модели + итоговая "
        "оценка как % от бесплатной 5-часовой сессии. "
        "Числа эвристические — служат для планирования бюджета и контроля "
        "сложности спектра."
    )
    out.append("")
    out.append("| Этап | Базово, ток. | Сложность, ток. | Итого, ток. | Почему |")
    out.append("|---|---|---|---|---|")
    for s in by_stage:
        # Маппинг английского stage_id на русский префикс «Шаг» (для
        # сортировки). stage_name_ru уже начинается с «Шаг N: ...».
        out.append(
            f"| {s.get('stage_name_ru','?')} "
            f"| {int(s.get('tokens_baseline') or 0):,} "
            f"| {int(s.get('tokens_complexity') or 0):,} "
            f"| **{int(s.get('tokens_total') or 0):,}** "
            f"| {s.get('why','')} |".replace(",", " ")
        )
    auto_sum = sum(int(s.get("tokens_total") or 0) for s in by_stage)
    out.append("")
    # F-132 / v1.17.7 — единая формулировка:
    # «~N токенов или M% от бесплатной 5-часовой сессии».
    tokens_part = f"~{total:,} токенов".replace(",", " ")
    pct_part = f"{pct:.1f}% от бесплатной 5-часовой сессии"
    if override:
        out.append(
            f"**Авто-оценка по этапам**: {auto_sum:,} токенов.  ".replace(",", " ")
        )
        out.append(
            f"**Итог (CLI override)**: **{tokens_part}** или "
            f"**{pct_part}** (бюджет {budget:,} ток.).".replace(",", " ")
        )
    else:
        out.append(
            f"**Итог**: **{tokens_part}** или **{pct_part}** "
            f"(бюджет {budget:,} ток., = сумма по этапам).".replace(",", " ")
        )
    if detail:
        out.append("")
        out.append(f"_{detail}_")
    return "\n".join(out)


_WARNING_TRANSLATIONS = {
    "Efficiency curve not loaded": "Кривая эффективности не загружена",
    "Bq / MDA suite may be incomplete":
        "набор Бк / MDA может быть неполным",
    "Cascade summing warning":
        "Предупреждение о каскадном суммировании",
    "Stage-": "Эскалация уровня ",
    "escalation suggested": "рекомендуется",
}


def _translate_warning(w: str) -> str:
    out = w or ""
    for en, ru in _WARNING_TRANSLATIONS.items():
        out = out.replace(en, ru)
    return out


# F-386 / BUG-40: Russian-only renderer for structured warning dicts.
# The English ``message`` field is NEVER emitted into the MD report;
# instead, we look up the warning ``code`` and render a hand-written
# RU equivalent. Unknown codes fall back to a generic Russian skeleton
# that quotes only the safe structured fields (code, original_detector,
# mapped_to, severity).
_WARNING_SEVERITY_RU = {
    "LOW": "низкая",
    "MEDIUM": "средняя",
    "HIGH": "высокая",
    "CRITICAL": "критическая",
}


def _render_warning_dict_ru(w: dict) -> str:
    """Render a structured warning dict to F-386-compliant RU prose.

    Per BUG-40, structured warnings carry ``code`` + ``message`` (EN) +
    domain-specific fields. The English message is FORBIDDEN in the RU
    Markdown report; this helper picks the code-specific RU template.
    """
    code = str(w.get("code") or "").strip()
    sev_ru = _WARNING_SEVERITY_RU.get(
        str(w.get("severity") or "").upper(), str(w.get("severity") or "")
    )
    if code == "DETECTOR_CYRILLIC_LATIN_FALLBACK":
        orig = str(w.get("original_detector") or "").strip() or "?"
        mapped = str(w.get("mapped_to") or "").strip() or "?"
        return (
            f"Имя детектора «{orig}» (кириллица) сопоставлено с «{mapped}» "
            f"(латиница); загружен профиль «{mapped}». Исходный "
            f"кириллический профиль детектора не зарегистрирован. "
            f"Результаты активности могут отражать несоответствие "
            f"эффективности. Серьёзность: {sev_ru}."
        )
    if code == "EFFICIENCY_DETECTOR_SERIAL_MISMATCH":
        # T41 (BUG-40 (b)) — silent content-fallback class. F-115:
        # anonymized to serial+year tuple (no full detector header).
        exp_sy = w.get("expected_serial_year") or []
        act_sy = w.get("actual_serial_year") or []
        exp_sn = exp_sy[0] if len(exp_sy) >= 1 else "?"
        exp_yy = exp_sy[1] if len(exp_sy) >= 2 else "?"
        act_sn = act_sy[0] if len(act_sy) >= 1 else "?"
        act_yy = act_sy[1] if len(act_sy) >= 2 else "?"
        return (
            f"Файл калибровки эффективности соответствует другому "
            f"физическому экземпляру прибора: спектр снят на приборе "
            f"с серийным номером {exp_sn} (год {exp_yy}), а в файле "
            f"калибровки эффективности записан прибор с серийным "
            f"номером {act_sn} (год {act_yy}). Кривая эффективности "
            f"относится к другому прибору; значения активности и "
            f"удельной активности могут быть существенно искажены. "
            f"Серьёзность: {sev_ru}."
        )
    # Unknown code — emit a generic Russian skeleton that does NOT
    # surface the English ``message`` (F-386 gate).
    return (
        f"Структурное предупреждение «{code or '—'}» "
        f"(серьёзность: {sev_ru or '—'})."
    )


def _ru_filter(text: str) -> str:
    """Pass arbitrary text through the interactive_html translation map
    so the Markdown report stays free of English leaks (F-108).
    Lines for which the translator returns ``None`` (residual English)
    are replaced by the literal "—".
    """
    try:
        from gamma.reporting.interactive_html import _translate_note_line
    except Exception:
        return text
    out = _translate_note_line(text)
    return out if out else ""


__all__ = ["build_markdown_report"]
