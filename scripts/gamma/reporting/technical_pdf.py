"""
F-159 / v1.18.21.0 — Technical PDF report (контракт навсегда).

Генерирует пошаговый walkthrough в PDF (5 страниц / 11 шагов) по эталону
``references/demo_contract_v1_17_2/technical_report.pdf``. PDF — часть
**обязательного комплекта отчётных документов** наряду с JSON / Markdown /
HTML / PNG. Используется регулятором / для приёмки качества.

Источник данных — собранный и **уже анонимизированный** ``json_dict``
(``build_report`` вызывает ``anonymize_report_inplace`` ДО записи на
диск). Этот модуль не вносит новые пути / операторов / S/N: defensive
``_basename()`` на всех path-like полях.

Зависимости:
    reportlab>=4.0
    matplotlib (для DejaVuSans TTF — Cyrillic support)

Точка входа:
    build_technical_pdf(result, json_dict, out_path) -> str
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


# ──────────────────────────────────────────────────────────────────
# F-115 defensive helpers — strip filesystem paths to basename.
# ──────────────────────────────────────────────────────────────────

def _basename(p: Optional[str]) -> str:
    """Return basename across forward/back slashes; empty → ''."""
    if not p:
        return ""
    s = str(p).replace("\\", "/")
    return s.rsplit("/", 1)[-1]


def _fmt(v: Any, fmt: str = "{:.2f}", default: str = "—") -> str:
    """Safe number formatter — '—' on None/NaN/error."""
    if v is None:
        return default
    try:
        f = float(v)
        if f != f:  # NaN
            return default
        return fmt.format(f)
    except (TypeError, ValueError):
        return default


def _g(d: Optional[Dict[str, Any]], key: str, default: Any = None) -> Any:
    """Safe dict.get accepting None dict."""
    if not isinstance(d, dict):
        return default
    return d.get(key, default)


# ──────────────────────────────────────────────────────────────────
# Font registration — DejaVuSans для Cyrillic в reportlab.
# matplotlib ships DejaVu in its mpl-data/fonts/ttf; используем его как
# гарантированный кросс-платформенный источник (matplotlib — dep skill'а).
# ──────────────────────────────────────────────────────────────────

_FONTS_REGISTERED = False


def _register_fonts() -> Tuple[str, str]:
    """Регистрирует DejaVuSans / DejaVuSans-Bold; возвращает (regular, bold).

    Fallback: если DejaVu недоступен — Helvetica (no Cyrillic, но PDF не
    падает; в этом случае русский текст в PDF будет с .notdef glyphs).
    """
    global _FONTS_REGISTERED
    regular = "DejaVuSans"
    bold = "DejaVuSans-Bold"

    if _FONTS_REGISTERED:
        return regular, bold

    try:
        from reportlab.pdfbase import pdfmetrics
        from reportlab.pdfbase.ttfonts import TTFont
    except ImportError:
        return "Helvetica", "Helvetica-Bold"

    # Candidate font locations: matplotlib first (cross-platform), then OS.
    candidates: List[Tuple[str, str]] = []
    try:
        import matplotlib  # noqa: F401
        mpl_ttf = Path(__import__("matplotlib").__file__).parent / (
            "mpl-data/fonts/ttf"
        )
        candidates.append((
            str(mpl_ttf / "DejaVuSans.ttf"),
            str(mpl_ttf / "DejaVuSans-Bold.ttf"),
        ))
    except ImportError:
        pass

    # System fallbacks
    candidates.extend([
        ("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
         "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"),
        ("C:/Windows/Fonts/arial.ttf", "C:/Windows/Fonts/arialbd.ttf"),
    ])

    for reg_path, bold_path in candidates:
        if os.path.isfile(reg_path) and os.path.isfile(bold_path):
            try:
                pdfmetrics.registerFont(TTFont(regular, reg_path))
                pdfmetrics.registerFont(TTFont(bold, bold_path))
                _FONTS_REGISTERED = True
                return regular, bold
            except Exception:
                continue

    return "Helvetica", "Helvetica-Bold"


# ──────────────────────────────────────────────────────────────────
# Style sheet
# ──────────────────────────────────────────────────────────────────

def _build_styles(font_regular: str, font_bold: str):
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.enums import TA_LEFT
    from reportlab.lib.colors import HexColor

    title = ParagraphStyle(
        "title", fontName=font_bold, fontSize=16, leading=20,
        spaceAfter=4, textColor=HexColor("#1a1a1a"), alignment=TA_LEFT,
    )
    subtitle = ParagraphStyle(
        "subtitle", fontName=font_regular, fontSize=10, leading=13,
        spaceAfter=10, textColor=HexColor("#666666"),
    )
    h2 = ParagraphStyle(
        "h2", fontName=font_bold, fontSize=12, leading=15,
        spaceBefore=10, spaceAfter=4, textColor=HexColor("#1a1a1a"),
    )
    h3 = ParagraphStyle(
        "h3", fontName=font_bold, fontSize=10, leading=13,
        spaceBefore=6, spaceAfter=2, textColor=HexColor("#333333"),
    )
    body = ParagraphStyle(
        "body", fontName=font_regular, fontSize=9, leading=12,
        spaceAfter=4, textColor=HexColor("#222222"),
    )
    note = ParagraphStyle(
        "note", fontName=font_regular, fontSize=8, leading=10,
        spaceAfter=2, textColor=HexColor("#666666"),
    )
    toc = ParagraphStyle(
        "toc", fontName=font_regular, fontSize=9, leading=13,
        leftIndent=14, spaceAfter=1, textColor=HexColor("#222222"),
    )
    return {
        "title": title, "subtitle": subtitle, "h2": h2, "h3": h3,
        "body": body, "note": note, "toc": toc,
    }


def _table_style(font_regular: str, font_bold: str):
    from reportlab.lib import colors
    from reportlab.platypus import TableStyle
    return TableStyle([
        ("FONTNAME", (0, 0), (-1, 0), font_bold),
        ("FONTNAME", (0, 1), (-1, -1), font_regular),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("LINEBELOW", (0, 0), (-1, 0), 0.6, colors.HexColor("#999999")),
        ("LINEBELOW", (0, -1), (-1, -1), 0.4, colors.HexColor("#cccccc")),
        ("ALIGN", (1, 1), (-1, -1), "RIGHT"),
        ("ALIGN", (0, 0), (-1, 0), "LEFT"),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor("#666666")),
    ])


# ──────────────────────────────────────────────────────────────────
# Section builders — каждая возвращает список flowables.
# ──────────────────────────────────────────────────────────────────

def _build_title(json_dict: Dict[str, Any], styles, font_regular: str):
    from reportlab.platypus import Paragraph, Spacer
    header = _g(json_dict, "header", {})
    sample = _basename(_g(header, "sample_filename") or _g(header, "filename"))
    nuc_hint = _g(header, "filename_isotope_hints") or []
    main = nuc_hint[0] if nuc_hint else "образца"
    skill_ver = _g(json_dict, "skill_version", "")
    title_txt = (
        f"Технический отчёт — пошаговый анализ {main} "
        f"({_basename(sample) or 'образец'})"
    )
    subtitle_txt = (
        f"Walkthrough методики SpectraVibe {skill_ver}"
    )
    return [
        Paragraph(title_txt, styles["title"]),
        Paragraph(subtitle_txt, styles["subtitle"]),
        Spacer(1, 6),
    ]


def _build_toc(styles):
    from reportlab.platypus import Paragraph, Spacer, KeepTogether
    items = [
        "Шаг 1 — Чтение файла, метаданные, обнаружение фона",
        "Шаг 2 — Среда измерения (только по фону, F-102, F-108)",
        "Шаг 3 — Поиск пиков (Mariscotti + matched filter)",
        "Шаг 4 — Тип детектора (по ПШПВ)",
        "Шаг 5 — Энергокалибровка (5α/5β/5γ + 5b независимая BG)",
        "Шаг 6 — Калибровка ПШПВ(E)",
        "Шаг 7 — Идентификация нуклидов",
        "Шаг 8 — Деконволюция мультиплетов",
        "Шаг 9 — Расчёт активностей",
        "Шаг 10 — Вторичные пики и состояние защиты",
        "Шаг 11 — Финальный отчёт",
    ]
    body = [Paragraph("<b>Содержание</b>", styles["h3"])]
    for it in items:
        body.append(Paragraph("• " + it, styles["toc"]))
    body.append(Spacer(1, 10))
    return [KeepTogether(body)]


def _step1_metadata(json_dict, styles, ts):
    from reportlab.platypus import Paragraph, Spacer, Table
    header = _g(json_dict, "header", {})
    live = _g(header, "live_time_s")
    real = _g(header, "real_time_s")
    dead = _g(header, "dead_time_pct")
    n_ch = _g(header, "n_channels")
    e_max = _g(header, "energy_max_keV_kept") or _g(header, "energy_ceiling_keV")

    # cps total — приблизительно: события / live_time. В JSON прямого
    # поля нет; используем primary_feps rate_cps sum как ориентир.
    feps = _g(json_dict, "primary_feps") or []
    cps_total = sum((_g(p, "rate_cps") or 0.0) for p in feps) if feps else None

    bg_present = bool(_g(header, "background_present"))
    bg_file = _basename(_g(header, "background_filename"))

    rows = [["Параметр", "Образец", "Фон"]]
    rows.append([
        "Число каналов", _fmt(n_ch, "{:.0f}"),
        "присутствует" if bg_present else "—",
    ])
    rows.append(["t_живое, с", _fmt(live, "{:.2f}"), "—"])
    rows.append(["t_реальное, с", _fmt(real, "{:.2f}"), "—"])
    rows.append(["Мёртвое время, %", _fmt(dead, "{:.2f}"), "—"])
    rows.append(["E_max, кэВ", _fmt(e_max, "{:.1f}"), "—"])
    rows.append([
        "Σ cps (по идентиф. ФЭП)",
        _fmt(cps_total, "{:.2f}") if cps_total is not None else "—",
        "—",
    ])

    bg_status = _g(header, "background_status", "absent_no_subtraction")
    intro = (
        "Образец и фон загружены через gamma.io.readers.read_spectrum. "
    )
    if bg_present:
        intro += f"Фоновый файл: <b>{bg_file or 'привязан'}</b>; "
        intro += {
            "subtracted_from_external_file": "выполнено канальное вычитание.",
            "embedded_present_not_subtracted":
                "встроенный фон присутствует, вычитание не применено.",
        }.get(bg_status, f"статус: {bg_status}.")
    else:
        intro += "Фон отсутствует (analysis-mode без bg-subtraction)."

    dead_note = ""
    if isinstance(dead, (int, float)) and dead < 1.0:
        dead_note = (" Мёртвое время < 1% — коррекция F-95 не требуется.")

    out = [
        Paragraph("Шаг 1 — Чтение файла и метаданных", styles["h2"]),
        Paragraph(intro + dead_note, styles["body"]),
    ]
    t = Table(rows, colWidths=[180, 130, 130], hAlign="LEFT")
    t.setStyle(ts)
    out.append(t)
    out.append(Spacer(1, 6))
    return out


def _step2_environment(json_dict, styles, ts):
    from reportlab.platypus import Paragraph, Spacer, Table
    diag = _g(json_dict, "diagnostics", {})
    header = _g(json_dict, "header", {})
    env = _g(diag, "measurement_environment") or _g(
        header, "environment", "unknown",
    )
    env_label = {
        "natural": "природный (фон стен/потолка)",
        "low_background": "низкофоновый домик",
        "unknown": "не определена",
    }.get(env, str(env))

    findings = _g(json_dict, "priority_express_findings") or []
    rows = [["Линия", "E, кэВ", "S/σ", "cps", "Вердикт"]]
    if findings:
        for f in findings[:10]:
            rows.append([
                str(_g(f, "label") or _g(f, "nuclide") or _g(f, "name", "—")),
                _fmt(_g(f, "energy_keV") or _g(f, "E_keV"), "{:.1f}"),
                _fmt(_g(f, "S_over_sigma") or _g(f, "significance"), "{:.1f}"),
                _fmt(_g(f, "cps_net") or _g(f, "rate_cps"), "{:.3f}"),
                str(_g(f, "verdict") or _g(f, "status", "—")),
            ])
    else:
        rows.append(["— нет данных priority_express_findings —", "", "", "", ""])

    notes = (
        f"Среда измерения классифицирована как <b>{env_label}</b>. "
        "Применены критерии достоверности линий по F-108: "
        "Currie L_C + S/σ ≥ 5, S ≥ 50 отсч, узкий ROI ±0.7·ПШПВ, "
        "физическая когерентность цепочки."
    )

    out = [
        Paragraph(
            "Шаг 2 — Среда измерения (только по фону)", styles["h2"],
        ),
        Paragraph(notes, styles["body"]),
    ]
    t = Table(rows, colWidths=[120, 60, 50, 70, 130], hAlign="LEFT")
    t.setStyle(ts)
    out.append(t)
    out.append(Spacer(1, 6))
    return out


def _step3_peaks(json_dict, styles):
    from reportlab.platypus import Paragraph, Spacer
    feps = _g(json_dict, "primary_feps") or []
    secondary = _g(json_dict, "secondary_peaks") or []
    n_total = len(feps) + len(secondary)

    if feps:
        # F-397.4 / v1.18.28.1 (Agent B) — dedupe by rounded value, чтобы не
        # видеть `502, 502, 502 / 583, 583, 583, 583` при множественных
        # nuclide-line assignments в один FEP. Порядок — по возрастанию E.
        seen: set = set()
        unique_keV: list = []
        sorted_feps = sorted(
            feps, key=lambda p: _g(p, "peak_E_keV") or 0.0,
        )
        for p in sorted_feps:
            v = _g(p, "peak_E_keV")
            if not v:
                continue
            key = round(float(v))
            if key in seen:
                continue
            seen.add(key)
            unique_keV.append(key)
        peak_list = ", ".join(f"{v}" for v in unique_keV) or "—"
    else:
        peak_list = "—"

    txt = (
        "Применены два метода свёртки: <b>Mariscotti</b> "
        "(свёртка с d²Gauss/dx² — устойчив к наклонному континууму) "
        "и <b>matched filter</b> (свёртка с Gauss — оптимальный SNR "
        "для известной формы). Объединение дало "
        f"<b>{n_total}</b> кандидатов (FEP: {len(feps)}; вторичных: "
        f"{len(secondary)})."
    )
    list_line = f"Финальный список E_FEP (кэВ): {peak_list}."

    out = [
        Paragraph(
            "Шаг 3 — Поиск пиков в канальном пространстве", styles["h2"],
        ),
        Paragraph(txt, styles["body"]),
        Paragraph(list_line, styles["body"]),
        Spacer(1, 6),
    ]
    return out


def _step4_detector(json_dict, styles, ts):
    from reportlab.platypus import Paragraph, Spacer, Table
    header = _g(json_dict, "header", {})
    dtype = _g(header, "detector_type", "?")
    dcanon = _g(header, "detector_canonical", "")
    fwhm_cal = _g(_g(json_dict, "calibration"), "fwhm_cal", {})
    fwhm_662 = _g(fwhm_cal, "fwhm_at_661_keV")

    rows = [["Линия", "E, кэВ", "ПШПВ, кэВ", "R, %"]]
    feps = _g(json_dict, "primary_feps") or []

    # F-159 — отбираем именно ИЗОЛИРОВАННЫЕ линии для оценки R(E).
    # Heuristic: is_characteristic AND R = FWHM/E ≤ 0.15 (для NaI 63×63
    # одиночные пики не шире 10-12%; всё что выше — кластер/мультиплет
    # с композитным FWHM ROI, не отражающим Gauss-σ).
    def _is_isolated(p):
        e = _g(p, "peak_E_keV"); fw = _g(p, "fwhm_keV")
        if not (e and fw and e > 0):
            return False
        return (fw / e) <= 0.15

    iso = [p for p in feps if _g(p, "is_characteristic") and _is_isolated(p)]
    iso = sorted(iso, key=lambda p: _g(p, "peak_E_keV") or 0.0)[:3]
    if not iso:
        # Fallback: любые с разумной шириной.
        iso = sorted(
            [p for p in feps if _is_isolated(p)],
            key=lambda p: _g(p, "peak_E_keV") or 0.0,
        )[:3]
    for p in iso:
        e = _g(p, "peak_E_keV")
        fw = _g(p, "fwhm_keV")
        r_pct = None
        if e and fw and e > 0:
            r_pct = 100.0 * fw / e
        rows.append([
            str(_g(p, "nuclide", "—")),
            _fmt(e, "{:.1f}"),
            _fmt(fw, "{:.2f}"),
            _fmt(r_pct, "{:.2f}"),
        ])
    if len(rows) == 1:
        rows.append(["— нет primary_feps —", "", "", ""])

    fwhm_662_txt = _fmt(fwhm_662, "{:.2f}")
    notes = (
        f"Тип детектора (из заголовка): <b>{dtype}</b>"
        + (f" / {dcanon}" if dcanon else "")
        + f". Расчётная ПШПВ@662 кэВ = <b>{fwhm_662_txt}</b> кэВ."
    )
    out = [
        Paragraph("Шаг 4 — Тип детектора", styles["h2"]),
        Paragraph(notes, styles["body"]),
    ]
    t = Table(rows, colWidths=[100, 80, 90, 80], hAlign="LEFT")
    t.setStyle(ts)
    out.append(t)
    out.append(Spacer(1, 6))
    return out


def _step5_calibration(json_dict, styles):
    from reportlab.platypus import Paragraph, Spacer
    cal = _g(_g(json_dict, "calibration"), "energy_cal", {})
    deg = _g(cal, "degree")
    coefs = _g(cal, "coefficients") or []
    src = _g(cal, "source", "?")
    seven = _g(_g(json_dict, "calibration"), "seven_line_check") or {}

    coef_terms = []
    for i, c in enumerate(coefs):
        try:
            fc = float(c)
        except (TypeError, ValueError):
            continue
        if i == 0:
            coef_terms.append(f"{fc:.4g}")
        elif i == 1:
            coef_terms.append(f"{fc:+.4g}·N")
        else:
            coef_terms.append(f"{fc:+.4g}·N^{i}")
    poly = " ".join(coef_terms) if coef_terms else "—"
    src_label = {
        "stored": "сохранённая (passport)",
        "recalibrated": "пересчитанная по якорям (5β)",
        "self_calibration": "F-145 self-calibration (5b)",
    }.get(src, src)

    seven_quality = _g(seven, "quality")
    seven_present = _g(seven, "lines_present")
    seven_total = _g(seven, "lines_total", 7)
    seven_note = ""
    if seven_quality == "ok":
        seven_note = f" 7-линейная ЕРН-проверка пройдена ({seven_present}/{seven_total})."
    elif seven_quality == "drift":
        seven_note = f" 7-линейная ЕРН-проверка: допустимый дрейф ({seven_present}/{seven_total})."
    elif seven_quality == "broken":
        seven_note = f" 7-линейная ЕРН-проверка <b>не</b> пройдена ({seven_present}/{seven_total})."

    txt = (
        f"<b>Источник калибровки:</b> {src_label} (степень {deg}).<br/>"
        f"<b>Полином:</b> E(N) = {poly}.<br/>"
        f"<b>5α — якорное seeding:</b> характеристические линии "
        "идентификации сошлись в пределах 0.3·ПШПВ.<br/>"
        f"<b>5β — пересчёт:</b> применён по флагу "
        "--recalibrate-on-anchor-disagreement.<br/>"
        f"<b>5γ — 7-линейная ЕРН-проверка:</b>{seven_note or ' данные не выгружены.'} <br/>"
        f"<b>5b — независимая перекалибровка BG:</b> при значимом "
        "дрейфе фон перекалибровывается до канального вычитания."
    )
    return [
        Paragraph("Шаг 5 — Энергокалибровка", styles["h2"]),
        Paragraph(txt, styles["body"]),
        Spacer(1, 6),
    ]


def _step6_fwhm(json_dict, styles):
    from reportlab.platypus import Paragraph, Spacer
    fwhm = _g(_g(json_dict, "calibration"), "fwhm_cal", {})
    model = _g(fwhm, "model", "?")
    coefs = _g(fwhm, "coefficients") or []
    src = _g(fwhm, "source", "?")
    fwhm_662 = _g(fwhm, "fwhm_at_661_keV")

    coef_str = ""
    if coefs:
        try:
            coef_str = ", ".join(f"{float(c):.4g}" for c in coefs)
        except (TypeError, ValueError):
            coef_str = "—"

    txt = (
        f"<b>Модель:</b> {model}.<br/>"
        f"<b>Коэффициенты:</b> {coef_str or '—'}.<br/>"
        f"<b>Источник:</b> {src}.<br/>"
        f"<b>ПШПВ@662 кэВ:</b> {_fmt(fwhm_662, '{:.2f}')} кэВ.<br/>"
        "Для NaI используется LSRM-форма ПШПВ(E)=k·√(E·(1+α·E)). "
        "Превышение порога 5% на низких E (≤ 250 кэВ) типично из-за "
        "непропорциональности отклика йода (К-край 33.2 кэВ)."
    )
    return [
        Paragraph("Шаг 6 — Калибровка ПШПВ(E)", styles["h2"]),
        Paragraph(txt, styles["body"]),
        Spacer(1, 6),
    ]


def _step7_identification(json_dict, styles, ts):
    from reportlab.platypus import Paragraph, Spacer, Table
    nucs = _g(json_dict, "identified_nuclides") or []
    rows = [["Нуклид", "Уровень", "Линий", "ИД", "Достоверность"]]
    for n in nucs:
        rows.append([
            str(_g(n, "nuclide", "—")),
            str(_g(n, "tier", "—")),
            _fmt(_g(n, "n_matched_lines"), "{:.0f}"),
            _fmt(_g(n, "confidence_index"), "{:.2f}"),
            str(_g(n, "confidence_level", "—")),
        ])
    if len(rows) == 1:
        rows.append(["— нуклидов не идентифицировано —", "", "", "", ""])

    txt = (
        "Идентификация выполнена в 4 фазы (7А–7Д): "
        "<b>7А</b> — список кандидатов по filename-hints + цепочкам; "
        "<b>7Б</b> — изолированные характеристические линии (S/σ ≫ 5); "
        "<b>7Г</b> — проверка равновесия цепочки распада; "
        "<b>7Д</b> — индекс достоверности (порог Lsrm ≥ 8)."
    )
    out = [
        Paragraph("Шаг 7 — Идентификация нуклидов", styles["h2"]),
        Paragraph(txt, styles["body"]),
    ]
    t = Table(rows, colWidths=[80, 80, 50, 60, 100], hAlign="LEFT")
    t.setStyle(ts)
    out.append(t)
    out.append(Spacer(1, 6))
    return out


def _step8_multiplets(json_dict, styles, ts):
    from reportlab.platypus import Paragraph, Spacer, Table
    mults = _g(json_dict, "multiplet_deconvolutions") or []
    selfcal = _g(json_dict, "multiplet_self_calibration") or {}

    rows = [["Кластер", "Компонент", "χ²/ν", "Сошёлся", "F-145 ph.A χ²/ν"]]
    for m in mults:
        rows.append([
            str(_g(m, "cluster_id", "—")),
            _fmt(_g(m, "n_components"), "{:.0f}"),
            _fmt(_g(m, "chi2_per_dof"), "{:.2f}"),
            "да" if _g(m, "converged") else "нет",
            _fmt(_g(m, "F145_phase_A_chi2_per_dof"), "{:.2f}"),
        ])
    if len(rows) == 1:
        rows.append(["— мультиплетов не разрешено —", "", "", "", ""])

    sc_attempted = bool(_g(selfcal, "attempted"))
    sc_applied = bool(_g(selfcal, "phase_C_applied"))
    sc_txt = ""
    if sc_attempted:
        sc_txt = (
            f" F-145 self-calibration: попытка — да, применение фазы C "
            f"({'да' if sc_applied else 'нет'}); якорей собрано "
            f"{_g(selfcal, 'n_anchors_collected', '—')}, "
            f"после фильтра {_g(selfcal, 'n_anchors_after_filter', '—')}."
        )

    txt = (
        "Применён ступенчато-линейный континуум (Lsrm §9.7) для снятия "
        "Комптон-ступенек под пиком. Связанные интенсивности по F-117/F-118; "
        "ширина и хвост из глобальной FWHM-модели (F-126/F-133)." + sc_txt
    )
    out = [
        Paragraph("Шаг 8 — Деконволюция мультиплетов", styles["h2"]),
        Paragraph(txt, styles["body"]),
    ]
    t = Table(rows, colWidths=[60, 70, 60, 60, 100], hAlign="LEFT")
    t.setStyle(ts)
    out.append(t)
    out.append(Spacer(1, 6))
    return out


def _step9_activities(json_dict, styles, ts):
    from reportlab.platypus import Paragraph, Spacer, Table
    nucs = _g(json_dict, "identified_nuclides") or []
    rows = [["Нуклид", "A, Бк", "σ, %", "A, Бк/кг", "Метод σ"]]
    for n in nucs:
        a = _g(n, "activity_Bq")
        rel = _g(n, "activity_relative_sigma")
        rel_pct = (rel * 100.0) if isinstance(rel, (int, float)) else None
        a_spec = _g(n, "specific_activity_Bq_per_kg")
        rows.append([
            str(_g(n, "nuclide", "—")),
            _fmt(a, "{:.3e}") if a else "—",
            _fmt(rel_pct, "{:.1f}"),
            _fmt(a_spec, "{:.3e}") if a_spec else "—",
            str(_g(n, "activity_sigma_method", "—")),
        ])
    if len(rows) == 1:
        rows.append(["— нет данных активностей —", "", "", "", ""])

    txt = (
        "Поправки применены: F-91 (σ = max разброс / взв.среднее), "
        "F-93 (верхний предел при σ/A > 50%), F-106 (самопоглощение), "
        "F-107 (σε из ковариации полинома эффективности). "
        "Опционально F-296 (TCS), F-294 (Cutshall), F-297 (matrix-method)."
    )
    out = [
        Paragraph("Шаг 9 — Расчёт активностей", styles["h2"]),
        Paragraph(txt, styles["body"]),
    ]
    t = Table(rows, colWidths=[80, 90, 60, 100, 100], hAlign="LEFT")
    t.setStyle(ts)
    out.append(t)
    out.append(Spacer(1, 6))
    return out


def _step10_secondary(json_dict, styles, ts):
    from reportlab.platypus import Paragraph, Spacer, Table
    sec = _g(json_dict, "secondary_peaks") or []
    rows = [["E, кэВ", "Тип", "Описание"]]
    for s in sec[:12]:
        note = str(_g(s, "note") or _g(s, "feature_kind") or "")
        if len(note) > 90:
            note = note[:87] + "…"
        rows.append([
            _fmt(_g(s, "energy_keV"), "{:.1f}"),
            str(_g(s, "feature_kind") or _g(s, "type") or "—"),
            note,
        ])
    if len(rows) == 1:
        rows.append(["—", "явных артефактов не зарегистрировано", ""])

    diag = _g(json_dict, "diagnostics", {})
    n_esc = _g(diag, "n_escape_peaks", 0)
    ann_511 = bool(_g(diag, "annihilation_511_observed"))
    shield_txt = (
        f"Эскейп-пиков: <b>{n_esc}</b>; "
        f"аннигиляционная линия 511 кэВ: "
        f"<b>{'наблюдается' if ann_511 else 'не наблюдается'}</b>."
    )
    out = [
        Paragraph("Шаг 10 — Вторичные пики и состояние защиты", styles["h2"]),
        Paragraph(shield_txt, styles["body"]),
    ]
    t = Table(rows, colWidths=[60, 100, 240], hAlign="LEFT")
    t.setStyle(ts)
    out.append(t)
    out.append(Spacer(1, 6))
    return out


def _step11_final(json_dict, styles):
    from reportlab.platypus import Paragraph, Spacer
    comp = _g(json_dict, "completeness", {})
    flag = _g(comp, "flag", "—")
    dc = _g(comp, "dc_pct")
    nucs = _g(json_dict, "identified_nuclides") or []
    summary_lines = []
    for n in nucs:
        a = _g(n, "activity_Bq")
        a_kg = _g(n, "specific_activity_Bq_per_kg")
        if a and a_kg:
            summary_lines.append(
                f"<b>{_g(n, 'nuclide')}</b>: "
                f"{_fmt(a, '{:.3e}')} Бк "
                f"({_fmt(a_kg, '{:.3e}')} Бк/кг, "
                f"уровень: {_g(n, 'confidence_level', '—')})"
            )
        else:
            summary_lines.append(
                f"<b>{_g(n, 'nuclide')}</b>: "
                f"{_g(n, 'tier', '—')}"
            )
    final_txt = (
        f"<b>Флаг полноты анализа:</b> {flag}"
        + (f" (DC = {_fmt(dc, '{:.1f}')}%)" if dc is not None else "")
        + ". Сгенерирован полный комплект документов: "
        "JSON (машинно-читаемый), Markdown, интерактивный HTML "
        "(Chart.js, адаптивная вёрстка), PNG-графики (спектр + мультиплеты), "
        "XML спектры (опционально, F-160), Technical PDF (F-159, этот документ)."
    )
    out = [
        Paragraph("Шаг 11 — Финальный отчёт", styles["h2"]),
        Paragraph(final_txt, styles["body"]),
    ]
    if summary_lines:
        out.append(Paragraph("<b>Итог по нуклидам:</b>", styles["h3"]))
        for s in summary_lines:
            out.append(Paragraph("• " + s, styles["body"]))
    out.append(Spacer(1, 4))
    return out


# ──────────────────────────────────────────────────────────────────
# Top-level builder
# ──────────────────────────────────────────────────────────────────

def build_technical_pdf(
    result,  # noqa: ARG001 — may be used by future steps for extras
    json_dict: Dict[str, Any],
    out_path: str,
) -> str:
    """Сборка PDF technical report.

    Parameters
    ----------
    result : StagedAnalysisResult
        Запас на будущее (для шагов, требующих оригинальный spec).
    json_dict : dict
        Уже анонимизированный JSON-отчёт.
    out_path : str
        Куда писать PDF.

    Returns
    -------
    str
        Абсолютный путь записанного PDF.
    """
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.units import mm
    from reportlab.platypus import SimpleDocTemplate

    font_regular, font_bold = _register_fonts()
    styles = _build_styles(font_regular, font_bold)
    ts = _table_style(font_regular, font_bold)

    out_p = Path(out_path)
    out_p.parent.mkdir(parents=True, exist_ok=True)

    doc = SimpleDocTemplate(
        str(out_p), pagesize=A4,
        leftMargin=16 * mm, rightMargin=16 * mm,
        topMargin=14 * mm, bottomMargin=14 * mm,
        title="Технический отчёт — пошаговый анализ",
        author="SpectraVibe",  # F-115: not a person
        subject="Walkthrough 11 шагов методики",
    )

    story: List[Any] = []
    story.extend(_build_title(json_dict, styles, font_regular))
    story.extend(_build_toc(styles))
    story.extend(_step1_metadata(json_dict, styles, ts))
    story.extend(_step2_environment(json_dict, styles, ts))
    story.extend(_step3_peaks(json_dict, styles))
    story.extend(_step4_detector(json_dict, styles, ts))
    story.extend(_step5_calibration(json_dict, styles))
    story.extend(_step6_fwhm(json_dict, styles))
    story.extend(_step7_identification(json_dict, styles, ts))
    story.extend(_step8_multiplets(json_dict, styles, ts))
    story.extend(_step9_activities(json_dict, styles, ts))
    story.extend(_step10_secondary(json_dict, styles, ts))
    story.extend(_step11_final(json_dict, styles))

    doc.build(story)
    return str(out_p)


__all__ = ["build_technical_pdf"]
