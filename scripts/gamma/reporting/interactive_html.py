"""
Interactive HTML report renderer — canonical form (F-114 / v1.17.3).

The hand-crafted demo at ``references/demo_contract_v1_17_2/report.html`` is the
locked canonical shape for every full-report run. This renderer fills
the template skeleton at
``scripts/gamma/reporting/templates/interactive_v1_17_2.html`` with
per-sample data taken from the JSON report and the
:class:`StagedAnalysisResult`.

Public API
----------

    render_interactive_html(report, analysis_result, *, cost_estimate=None) -> str

`report` is the JSON dict produced by
:func:`gamma.reporting.build_json_report`. `analysis_result` is the
:class:`StagedAnalysisResult`. `cost_estimate` is an optional dict
``{'tokens': int, 'session_pct': str, 'detail': str}``.

Substitution uses ``str.replace`` (not ``str.format``) because the
template CSS contains literal ``{...}`` braces.
"""
from __future__ import annotations

import html as _html
import json
import logging
import math
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# BUG-25 / v1.18.31+ (Agent B) — debug-level logger для прозрачности
# фильтрации вырожденных мультиплетов. Используется в
# `_filter_meaningful_multiplets` (см. ниже).
log = logging.getLogger("gamma.reporting.interactive_html")

# F-452 / v1.33.0 — FwhmModel polymorphic API.
from gamma.identification.staged_pipeline import (  # noqa: E402
    fwhm_keV_at_energy as _fwhm_keV_at_energy,
)


# ──────────────────────────────────────────────────────────────────
# Colour palette (matches the demo)
# ──────────────────────────────────────────────────────────────────

_COL_TH    = "#7F77DD"   # Th-232 chain
_COL_NAT   = "#1D9E75"   # K-40, Cs-137, Co-60, Ra-226 family, Am-241
_COL_PHYS  = "#888780"   # Compton edges, escape peaks, sum peaks, XRF
# BUG-20 / v1.18.30+ (Agent B): distinct color for SECONDARY artefact peaks.
# Symptom: secondary processes (compton_edge, SE, DE, backscatter, 511, sum)
# раньше получали parent's chain color (Th-232 SE/DE были синие, как primary
# FEPs Th-цепочки). На общем спектре оператор не отличал артефакт детектора
# от реальной γ-линии. Orange #E8884F даёт контраст к _COL_TH (#7F77DD,
# blue) и _COL_NAT (#1D9E75, green) — WCAG AA contrast > 4.5:1 against
# white background; на dark тоже видно (relative luminance 0.34 vs 0.22 bg).
_COL_SECONDARY = "#E8884F"  # SE/DE/sum/511/backscatter/compton — все detector artefacts


_TH_CHAIN_NUCLIDES = {
    "Tl-208", "Pb-212", "Ac-228", "Bi-212", "Th-228", "Th-232",
    "Ra-228", "Po-216", "Po-212",
}
_U_CHAIN_NUCLIDES = {
    "Pb-214", "Bi-214", "Pb-210", "Bi-210", "Po-214", "Po-218",
    "Ra-226", "Rn-222", "U-238", "U-234", "Th-234", "Pa-234",
}
_NATURAL_NUCLIDES = {
    "K-40", "Cs-137", "Cs-134", "Co-60", "Am-241", "Na-22", "Be-7",
    "I-131", "I-125", "Eu-152", "Eu-154", "Cs-137m", "Ba-137m",
}

_T_HALF_RU = {
    "Tl-208": "3.05 мин",
    "Pb-212": "10.6 ч",
    "Ac-228": "6.15 ч",
    "Bi-212": "60.6 мин",
    "Th-228": "1.91 года",
    "Th-232": "1.40×10¹⁰ лет",
    "Ra-228": "5.75 года",
    "Pb-214": "26.8 мин",
    "Bi-214": "19.9 мин",
    "Ra-226": "1600 лет",
    "K-40":   "1.25×10⁹ лет",
    "Cs-137": "30.08 года",
    "Cs-134": "2.06 года",
    "Co-60":  "5.27 года",
    "Am-241": "432.6 года",
    "Na-22":  "2.60 года",
    "I-131":  "8.02 сут",
    "Eu-152": "13.5 года",
}


# F-111 / F-111b — Th-232 chain library lines (I_γ ≥ 0.5%) within the
# typical Gamma-1S analyzed range [50..3000 keV]. Used to inject
# placeholder entries for lines confirmed by chain equilibrium but
# not directly fit in this run, so peaks/rows/detail satisfy the
# single-source-of-truth invariant.
_TH232_CHAIN_LIBRARY = [
    # (nuclide, E_keV, I_gamma_pct, t_half_ru)
    # Ac-228 129.065: low-E Th chain anchor (I_apparent=2.42% in equilibrium).
    # CLAUDE.md gotcha: real line, not Compton/NE-artefact. Key diagnostic anchor.
    ("Ac-228", 129.065, 2.42, "6.15 ч"),
    ("Pb-212", 238.63, 43.6, "10.6 ч"),
    ("Ac-228", 209.30, 3.89, "6.15 ч"),
    ("Ac-228", 270.20, 3.46, "6.15 ч"),
    ("Ac-228", 277.40, 2.33, "6.15 ч"),
    ("Ac-228", 338.32, 11.27, "6.15 ч"),
    ("Ac-228", 463.00, 4.40, "6.15 ч"),
    ("Ac-228", 911.20, 25.8, "6.15 ч"),
    ("Ac-228", 964.77, 4.99, "6.15 ч"),
    ("Ac-228", 968.97, 15.8, "6.15 ч"),
    ("Ac-228", 1588.2, 3.22, "6.15 ч"),
    ("Ac-228", 1630.6, 1.51, "6.15 ч"),
    ("Bi-212", 727.33, 6.67, "60.6 мин"),
    ("Tl-208", 510.77, 22.6, "3.05 мин"),
    ("Tl-208", 583.19, 85.0, "3.05 мин"),
    ("Tl-208", 763.13, 1.79, "3.05 мин"),
    ("Tl-208", 860.56, 12.5, "3.05 мин"),
    ("Tl-208", 2614.5, 99.75, "3.05 мин"),
]


# RU translations for English tokens commonly seen in pipeline_notes
# and warnings.  Lines that still carry untranslated English tokens
# (4+ letters) after substitution are filtered out (D-04, F-108).
# NB: order matters — longer phrases come first so they win.
_NOTE_RU_MAP = [
    # Longer phrases first
    ("filename isotope hints driving candidate list",
     "подсказки изотопа из имени файла формируют список кандидатов"),
    ("Step 7A.1: filename isotope hints driving candidate list",
     "Шаг 7А.1: подсказки изотопа из имени файла формируют список кандидатов"),
    ("Filename binds source to", "Имя файла привязывает источник к"),
    ("chain suppression", "подавление цепочки"),
    ("isotope hints", "подсказки по изотопам"),
    ("Multiplet deconvolution", "Разложение мультиплета"),
    ("overlap threshold", "порог наложения"),
    ("cluster(s)", "кластер(ов)"),
    ("clusters", "кластеров"),
    ("cluster", "кластер"),
    ("peak area(s) replaced", "площадей пиков обновлено"),
    ("peak area replaced", "площадь пика обновлена"),
    ("driving candidate list", "формирует список кандидатов"),
    ("is often Tl-208", "часто соответствует Tl-208"),
    ("Step 7A.1", "Шаг 7А.1"),
    ("Step 7A", "Шаг 7А"),
    ("Step 5α", "Шаг 5α"),
    ("Step 5", "Шаг 5"),
    ("Step 1", "Шаг 1"),
    ("Stage 1", "Шаг 1"),
    ("Stage 2", "Шаг 2"),
    ("Stage 3", "Шаг 3"),
    # F-108 / v1.17.6 — covers pipeline RECOMMENDATION notes like
    # "Stage-2 рекомендуется" (hyphenated form).
    ("Stage-1", "Шаг 1"),
    ("Stage-2", "Шаг 2"),
    ("Stage-3", "Шаг 3"),
    ("FWHM model source", "Источник модели ПШПВ"),
    ("default_NaI_63x63", "модель по умолчанию NaI 63×63"),
    ("WARNING", "ПРЕДУПРЕЖДЕНИЕ"),
    ("Warning", "Предупреждение"),
    ("warning", "предупреждение"),
    ("Cs-134 is a known cascade emitter",
     "Cs-134 — известный каскадный излучатель"),
    ("Cascade summing warning", "Предупреждение о каскадном суммировании"),
    ("cascade summing", "каскадное суммирование"),
    ("Background SUBTRACTED", "Фон ВЫЧТЕН"),
    ("Background NOT subtracted", "Фон НЕ вычтен"),
    ("external file", "внешний файл"),
    ("embedded bg", "встроенный фон"),
    ("from filename", "по имени файла"),
    ("dominant", "доминирует"),
    ("trump-card rule", "правило определяющей линии"),
    ("trump card", "определяющая линия"),
    ("gain drift", "дрейф усиления"),
    ("upper limit", "верхний предел"),
    ("Efficiency curve not loaded",
     "Кривая эффективности не загружена"),
    ("Bq / MDA suite may be incomplete",
     "набор Бк / MDA может быть неполным"),
    ("escalation suggested", "рекомендуется"),
    # F-89d chain suppression preamble
    ("chain DOMINANT", "цепочка ДОМИНИРУЕТ"),
    ("DOMINANT", "ДОМИНИРУЕТ"),
    ("dropped", "исключены"),
    ("from identifications", "из идентификации"),
    ("suppressed", "подавлены"),
    ("matched at", "совпала при"),
    ("evidence", "признаки"),
    ("would require", "потребует"),
    ("Bi-214 quartet", "квартет Bi-214"),
    ("not satisfied here", "здесь не выполняется"),
    ("on NaI 63×63 the 609 keV peak is often Tl-208 583 keV shifted by Compton overlap",
     "на NaI 63×63 пик 609 кэВ часто соответствует Tl-208 583, смещённому из-за наложения Комптона"),
    ("U-238 evidence (Bi-214 Ra-pair 609+1764) is suppressed because",
     "Признаки U-238 (Ra-пара Bi-214 609+1764) подавлены, поскольку"),
    ("U-238 dominance would require the Bi-214 quartet (≥3 of 609, 1120, 1764, 2204) — not satisfied here",
     "Доминирование U-238 потребовало бы квартет Bi-214 (≥3 из 609, 1120, 1764, 2204) — здесь не выполняется"),
    ("Tl-208 / Pb-212 / Ac-228 / Bi-212 confirmed as strong-prior candidates",
     "Tl-208 / Pb-212 / Ac-228 / Bi-212 подтверждены как кандидаты с сильным приоритетом"),
    ("confirmed as strong-prior candidates",
     "подтверждены как кандидаты с сильным приоритетом"),
    ("strong-prior candidates", "кандидаты с сильным приоритетом"),
    ("empirical range", "эмпирический диапазон"),
    ("single escape", "пик одиночного вылета"),
    ("double escape", "пик двойного вылета"),
    ("sum-peak", "сумматорный пик"),
    ("sum peak", "сумматорный пик"),
    ("Compton edge", "Комптон-край"),
    ("compton edge", "комптон-край"),
    ("compton", "комптон"),
    ("Compton", "Комптон"),
    ("of Ac-228", "от Ac-228"),
    ("of Tl-208", "от Tl-208"),
    ("of Pb-212", "от Pb-212"),
    ("of Bi-212", "от Bi-212"),
    ("from Tl-208", "от Tl-208"),
    # F-386.1 / v1.18.28 (Agent B) — «вылет», не «ускользание» (F-386 hardlock).
    ("escape", "вылет"),
    ("Ra-pair", "Ra-пара"),
    ("pair", "пара"),
    ("Th K-series", "Th K-серия"),
    ("K-series", "K-серия"),
    ("series", "серия"),
    ("Filename binds source to Th-232 chain only.",
     "Имя файла привязывает источник только к цепочке Th-232."),
    ("Filename binds source to U-238 chain only.",
     "Имя файла привязывает источник только к цепочке U-238."),
]


# Whitelist of ASCII tokens (≥4 letters) allowed in HTML/MD body text.
# Anything else triggers note-line dropping.  Case-insensitive match.
_RU_BODY_ALLOWED_TOKENS = {
    # units / abbreviations
    "kev", "keV", "cps", "bq", "mda", "fwhm", "roi", "tcs", "xrf",
    "nai", "hpge", "labr", "cebr", "cdznte", "labr3", "iso", "iaea",
    "lsrm", "stage", "ern", "html", "json", "pdf", "udp",
    # nuclide stems
    "cs", "co", "am", "eu", "pb", "bi", "tl", "ac", "ra", "rn", "th",
    "po", "ba", "na", "ce", "la", "br",
    # math/physics symbols (most are non-ASCII already)
    "alpha", "beta", "gamma", "sigma", "chi",
    # F-IDs
    # filenames that may leak (allow .spe/.efr/.efa)
    "spe", "efr", "efa", "src", "txt", "csv", "html", "xml",
}


_F_RULE_RE = re.compile(r"\bF-\d+[a-z]?\b", re.IGNORECASE)
_NUCLIDE_RE = re.compile(r"\b(?:Cs|Co|Am|Eu|Pb|Bi|Tl|Ac|Ra|Rn|Th|"
                         r"Po|Ba|Na|K|U|I|Be|Mn)-\d+m?\b")
_NUMBER_RE = re.compile(r"\b\d+(?:[\.,]\d+)?\b")
_ASCII_WORD_RE = re.compile(r"[A-Za-z]{4,}")


def _translate_note_line(line: str) -> Optional[str]:
    """Translate an English pipeline_notes line to Russian.

    Returns ``None`` if any English token of 4+ letters remains after
    mapping — the line is dropped silently (D-04).
    """
    if not line:
        return None
    s = str(line)
    for en, ru in _NOTE_RU_MAP:
        s = s.replace(en, ru)
    # Now strip protected tokens (F-IDs, nuclides, numbers)
    probe = _F_RULE_RE.sub(" ", s)
    probe = _NUCLIDE_RE.sub(" ", probe)
    probe = _NUMBER_RE.sub(" ", probe)
    # Find lingering English words ≥4 letters
    for m in _ASCII_WORD_RE.finditer(probe):
        tok = m.group(0).lower()
        if tok in _RU_BODY_ALLOWED_TOKENS:
            continue
        # Unmapped English token — drop the whole line (better silent
        # than English-leaked).
        return None
    return s

_GEOMETRY_RU = {
    "marinelli": "Маринелли",
    "denta": "Дента",
    "petri": "Петри",
    "petri60": "Петри 60 мл",
    "denta120": "Дента 120 мл",
    "point": "точечная",
    "point5cm": "Точ-5 см",
    "point25cm": "Точ-25 см",
}


# ──────────────────────────────────────────────────────────────────
# F-334.1 / v1.18.18.7 — Background library lines for bg-view markers
# ──────────────────────────────────────────────────────────────────
#
# Когда в HTML-отчёте toggle переключён в режим «Фон» (F-332), маркеры
# должны соответствовать характерным линиям ФОНА, а не sample-пикам.
# Берём найденные peak-search на bg энергии и сопоставляем со
# справочными линиями (NORM + аннигиляция + Pb-XRF от защиты).
# Tolerance 8 keV — соответствует FWHM NaI 63×63 на средних энергиях.

_BG_LINES_DICT: List[Tuple[float, str]] = [
    # (E_keV, isotope)  — отсортировано по возрастанию энергии
    (40.0,   "Pb K-XRF"),       # K-флуоресценция свинца защиты
    (74.8,   "Pb Kα1"),
    (84.9,   "Pb Kα2"),
    (87.3,   "Pb Kβ"),
    (185.7,  "U-235"),          # часто в Pb-защите как примесь
    (238.6,  "Pb-212"),         # Th-цепочка
    (295.2,  "Pb-214"),         # U-цепочка
    (338.3,  "Ac-228"),         # Th-цепочка
    (351.9,  "Pb-214"),         # U-цепочка
    (511.0,  "Аннигиляция"),    # e+/e-
    (583.2,  "Tl-208"),         # Th-цепочка
    (609.3,  "Bi-214"),         # U-цепочка
    (661.7,  "Cs-137"),         # антропогенное загрязнение
    (727.3,  "Bi-212"),
    (768.4,  "Bi-214"),
    (911.2,  "Ac-228"),         # Th-цепочка, доминантный
    (969.0,  "Ac-228"),
    (1120.3, "Bi-214"),
    (1238.1, "Bi-214"),
    (1377.7, "Bi-214"),
    (1460.8, "K-40"),           # K-40 — повсюду, доминирует в фоне
    (1764.5, "Bi-214"),         # U-цепочка
    (2204.1, "Bi-214"),
    (2614.5, "Tl-208"),         # Th-цепочка, верхняя линия
]


def _match_bg_isotope(
    E_keV: float,
    tolerance_keV: Optional[float] = None,
) -> Optional[Tuple[str, float]]:
    """F-334.1 — return (isotope_label, E_lib_keV) for closest bg-library line.

    Tolerance логика:
      • Если ``tolerance_keV`` задан явно — используется как absolute window.
      • Иначе FWHM-зависимая (NaI 63×63 model): δ_E = 0.5·FWHM где
        FWHM² ≈ 3.3·E (kэВ²) → FWHM(E) ≈ √(3.3·E). Для E=200: FWHM≈26→δ=13;
        для E=1461: FWHM≈69→δ=35. Это покрывает gain-drift между sample и
        bg измерениями (особенно когда bg-файл из другого года).

    Returns None если нет совпадений; иначе (isotope, E_lib_keV).
    """
    if E_keV is None or math.isnan(E_keV) or E_keV <= 0:
        return None
    if tolerance_keV is None:
        # Default NaI 63×63 FWHM model: FWHM² = 3.3·E (kэВ²)
        # BUG-37 / RAG-043 — FWHM² floor lowered from 1.0 → 0.01 (0.1 keV
        # FWHM floor) for consistency with `fwhm_keV_at_energy` (see
        # `staged_pipeline.py:530`). In practice this floor is unreachable
        # here because `3.3·E` at E ≥ 50 keV gives val ≥ 165 → FWHM ≥ 12.8
        # keV; the 8 keV outer `max(8.0, ...)` already bounds tolerance.
        # Keeping the safety net only for numerical edge cases.
        fwhm = math.sqrt(max(3.3 * float(E_keV), 0.01))
        tolerance_keV = max(8.0, 0.5 * fwhm)
    best: Optional[Tuple[str, float, float]] = None
    for E_lib, isotope in _BG_LINES_DICT:
        delta = abs(E_keV - E_lib)
        if delta <= tolerance_keV and (best is None or delta < best[2]):
            best = (isotope, E_lib, delta)
    if best is None:
        return None
    return (best[0], best[1])


def _detect_bg_peaks(
    analysis_result,
    top_n: int = 20,
    sigma_threshold: float = 3.0,
) -> List[Dict[str, Any]]:
    """F-334.1 — peak-search на bg спектре + matching против static dict.

    Возвращает список ``[{e, label, isotope, E_lib, intensity, is_top}]``
    отсортированный по убыванию intensity. ``top_n`` ограничивает кол-во
    точек на графике (читаемость). ``is_top`` помечает первые 5 для
    full-label (F-334.4-параллельно).

    Возвращает [] если:
      • фон не вычитался (`bg_counts_on_sample_grid` is None);
      • peak-search вернул 0 пиков (плохой/пустой фон);
      • импорт convolution_peak_search упал (defensive).

    Sigma threshold 3.0 (vs 2.5 для sample) — bg имеет худшую статистику,
    нужен более жёсткий cutoff чтобы не показывать шум.
    """
    bg_grid = getattr(analysis_result, "bg_counts_on_sample_grid", None)
    spec = getattr(analysis_result, "spec", None)
    fwhm_model = getattr(analysis_result, "fwhm_model", None)
    if bg_grid is None or spec is None or fwhm_model is None:
        return []
    try:
        import numpy as np
        bg_scaled = np.asarray(bg_grid, dtype=np.float64)
        if bg_scaled.size < 50:
            return []
        # F-58 contract: bg_counts_on_sample_grid = np.interp(...) * k где
        # k = t_sample/t_bg < 1. Для peak-search нужны native bg statistics
        # (иначе пики «утоплены» в шум после scaling вниз). Делим на k чтобы
        # восстановить эффективные исходные counts на sample-сетке.
        k = float(getattr(analysis_result, "bg_scale_factor", 0.0) or 0.0)
        if 0 < k < 1.0:
            bg_native = bg_scaled / k
        else:
            bg_native = bg_scaled
        if float(bg_native.sum()) < 100.0:
            return []
        from gamma.peaks.convolution_search import convolution_peak_search
        from gamma.identification.staged_pipeline import _make_fwhm_at_channel
        fwhm_at_ch = _make_fwhm_at_channel(spec, fwhm_model)
        peaks_raw = convolution_peak_search(
            bg_native,
            fwhm_channels=fwhm_at_ch,
            sigma_threshold=sigma_threshold,
            min_separation_factor=0.6,
            edge_margin=10,
        )
    except Exception:
        return []
    if not peaks_raw:
        return []

    enriched: List[Dict[str, Any]] = []
    for p in peaks_raw:
        try:
            e_kev = spec.channel_to_energy(int(p.channel))
            if e_kev is None:
                continue
            e_kev = float(e_kev)
        except Exception:
            continue
        intensity = float(getattr(p, "height", 0.0) or 0.0)
        if intensity <= 0:
            continue
        match = _match_bg_isotope(e_kev)  # FWHM-derived tolerance (default)
        if match is not None:
            isotope, E_lib = match
            label = f"{isotope} {E_lib:.1f}"
        else:
            isotope = ""
            E_lib = e_kev
            label = f"{e_kev:.1f} кэВ"
        enriched.append({
            "id": f"bg{int(round(e_kev))}",
            "e": round(e_kev, 1),
            "label": label,
            "isotope": isotope,
            "E_lib": round(E_lib, 1),
            "intensity": round(intensity, 1),
            "color": _COL_PHYS if not isotope else _chain_color(isotope, isotope),
        })

    # Sort by intensity desc, take top_n
    enriched.sort(key=lambda d: -d["intensity"])
    enriched = enriched[: top_n]
    # Mark top-5 for full labels (F-334.4 parity)
    for i, d in enumerate(enriched):
        d["is_top"] = i < 5
    # F-370 / v1.18.24.5 — FORCED top labels для диагностически важных
    # линий цепочек: даже слабые пики Tl-208 2614 / K-40 1461 / Bi-214 609
    # должны быть подписаны на view «Фон», иначе теряется источник
    # сигнала. Решение по аналогии с F-335.8 (boost SE/DE/sum для sample).
    # 2614 кэВ в фоне типично слабый (intensity ~5-10) и не попадает в
    # top-5 по чистому intensity-ranking, но его НАЛИЧИЕ — единственный
    # надёжный маркер Th-232 цепочки в pollution-free фоне.
    for d in enriched:
        iso = d.get("isotope") or ""
        E_lib = d.get("E_lib") or 0.0
        # Tl-208 2614 — top Th-232 line (99.75% I_γ, чистая ROI)
        if iso == "Tl-208" and 2600 <= E_lib <= 2625:
            d["is_top"] = True
        # K-40 1461 — одиночник, всегда диагностически важен
        elif iso == "K-40" and 1450 <= E_lib <= 1475:
            d["is_top"] = True
        # Bi-214 609 — top Ra-226 line
        elif iso == "Bi-214" and 605 <= E_lib <= 615:
            d["is_top"] = True
        # Ac-228 911 / 969 — Th-232 anchors
        elif iso == "Ac-228" and (905 <= E_lib <= 920 or 960 <= E_lib <= 975):
            d["is_top"] = True
        # Pb-212 238 — top Th-232 low-E line
        elif iso == "Pb-212" and 235 <= E_lib <= 245:
            d["is_top"] = True
        # Cs-137 661.6 — техногенный, всегда важен если он в фоне
        elif iso == "Cs-137" and 655 <= E_lib <= 668:
            d["is_top"] = True
    # BUG-43 / 2026-06-04 (Agent B) — 511 keV ID rule: если Tl-208 цепочка
    # подтверждена (есть пик у 583.2 или 2614.5 keV), атрибутировать 511 keV →
    # Tl-208 510.77 keV (Th-232 chain daughter), а не «Аннигиляция».
    # Физика: Tl-208 имеет γ-линию 510.77 keV (BR ~8.1% по nuclides.json).
    # На NaI 7% FWHM линии 510.77 и 511.00 неразрешимы (ΔE=0.23 keV << FWHM~36).
    # Когда Th-232 chain присутствует в фоне (маркер: Tl-208 583.2 / 2614.5),
    # атрибуция к Tl-208 методологически достовернее (ЛСРМ §14.4 наименьшего
    # действия). Источник: пользователь 2026-06-04.
    tl208_char_E_ranges = [(575, 592), (2600, 2625), (850, 870)]  # 583.2, 2614.5, 860.6
    tl208_confirmed = any(
        d.get("isotope") == "Tl-208"
        and any(lo <= (d.get("E_lib") or 0.0) <= hi for lo, hi in tl208_char_E_ranges)
        for d in enriched
    )
    if tl208_confirmed:
        for d in enriched:
            if d.get("isotope") == "Аннигиляция" and 500 <= d.get("E_lib", 0) <= 520:
                d["isotope"] = "Tl-208"
                d["E_lib"] = 510.77
                d["label"] = "Tl-208 510.8"
                d["color"] = _chain_color("Tl-208", "Tl-208")
                d["_511_tl208_override"] = True  # sentinel for augmenter
    # Re-sort by energy ascending for final list (matches sample peaks ordering)
    enriched.sort(key=lambda d: d["e"])
    return enriched


# ──────────────────────────────────────────────────────────────────
# F-334.4 / v1.18.18.7 — Mark top-N sample peaks for full labelling
# ──────────────────────────────────────────────────────────────────

def _mark_top_peaks(
    peaks: List[Dict[str, Any]],
    report: Dict[str, Any],
    top_n: int = 5,
) -> None:
    """F-334.4 / F-335.8 — set `is_top=True` для top-N пиков + special features.

    Правила (после F-335.8):
      • Primary FEPs (real fits) — score = peak_area_counts.
      • Secondary peaks — score = significance × 10.
      • F-335.8 BOOST: special methodological features (SE, DE, sum_peak,
        composite_cluster для Th/U цепочек) получают score = 1e8 →
        ВСЕГДА попадают в top независимо от area. Это критично для
        Th-232 спектров где SE 2103 / DE 1593 / 87-кластер — обязательные
        образовательные labels.
      • F-111 chain-completeness placeholders → score=0 (non-top).

    Top-N: сначала отбираются BOOST-features (forced); затем добавляются
    top-N остальных по area. Итоговое кол-во labelled может превысить
    top_n при наличии boost.
    """
    # F-335.8 — feature kinds которые получают unconditional boost.
    BOOST_KINDS = frozenset({
        "single_escape", "double_escape", "sum_peak",
        "composite_cluster", "cluster", "annihilation_511",
    })

    area_by_id: Dict[str, float] = {}
    for fp in (report.get("primary_feps") or []):
        e = _safe_float(fp.get("peak_E_keV"))
        if e is None:
            continue
        area = _safe_float(fp.get("peak_area_counts")) or 0.0
        area_by_id[_peak_id(e)] = float(area)
    for sp in (report.get("secondary_peaks") or []):
        e = _safe_float(sp.get("energy_keV"))
        if e is None:
            continue
        sig = _safe_float(sp.get("significance")) or 0.0
        existing = area_by_id.get(_peak_id(e), 0.0)
        area_by_id[_peak_id(e)] = max(existing, float(sig) * 10.0)

    # F-335.8 — forced top from BOOST_KINDS (always labelled)
    forced_top: set = set()
    for p in peaks:
        kind = p.get("feature_kind", "")
        if kind in BOOST_KINDS:
            forced_top.add(p.get("id", ""))

    # Top-N from area-scored peaks (excluding forced ones to avoid double-count)
    scored = [
        (p, area_by_id.get(p.get("id", ""), 0.0))
        for p in peaks
        if p.get("id") not in forced_top
    ]
    scored.sort(key=lambda x: -x[1])
    top_ids = set(forced_top)
    for i in range(min(top_n, len(scored))):
        if scored[i][1] > 0:
            top_ids.add(scored[i][0]["id"])

    # Always force-label max-E and min-E primary FEPs (key physics anchors:
    # Tl-208 2614 = highest-E calibration anchor; Ac-228 129 = low-E Th chain anchor).
    _pfe_all = [(float(p.get("e") or 0.0), p.get("id", ""))
                for p in peaks
                if p.get("feature_kind") in ("primary_fep", "chain_completeness")]
    _pfe_area = [(e, pid) for (e, pid) in _pfe_all
                 if area_by_id.get(pid, 0.0) > 0]
    if _pfe_area:
        top_ids.add(max(_pfe_area)[1])  # max-E with real area (e.g. Tl-208 2614)
    if _pfe_all:
        top_ids.add(min(_pfe_all)[1])   # min-E unconditionally (e.g. Ac-228 129)

    # F-385 / v1.18.26 — annotation tiebreaker: nuclide-FEP > secondary.
    #
    # При коллизии в пределах ±1·FWHM по энергии (NaI rough FWHM ≈ 6% E + 5 кэВ)
    # secondary peak (SE/DE/sum/annihilation/backscatter, is_secondary=True)
    # уступает primary FEP'у, даже если secondary был forced через BOOST_KINDS.
    # Это устраняет случаи когда label DE 1593 (Tl-208) перекрывает соседний
    # Ac-228 1588 nuclide-FEP на общем спектре Th-232.
    #
    # Edge case: если в окне нет primary_fep — secondary остаётся top (нет
    # конфликта, маркер нужен для образовательной ценности).
    primary_by_e: List[Tuple[float, str]] = [
        (float(p["e"]), p["id"]) for p in peaks
        if p.get("feature_kind") == "primary_fep"
    ]
    if primary_by_e:
        for p in peaks:
            if not p.get("is_secondary"):
                continue
            if p.get("id") not in top_ids:
                continue
            e_sec = float(p["e"])
            # Rough NaI FWHM in keV; tight enough to catch SE/DE vs FEP labels
            # на NaI 63×63 (FWHM ~50 кэВ @ 1500 кэВ) и не слишком тесный.
            fwhm_kev = 0.06 * e_sec + 5.0
            collide = [
                pid for (e_p, pid) in primary_by_e
                if abs(e_p - e_sec) <= fwhm_kev
            ]
            if collide:
                # secondary уступает primary slot; primary форсится в top_ids
                top_ids.discard(p.get("id", ""))
                for pid in collide:
                    top_ids.add(pid)

    # BUG-17 / v1.18.30+ (Agent B) — collision detection on top-labels.
    # Symptom: Ac-228 имеет 3 близкие линии 911.2 / 964.8 / 968.97 кэВ. При
    # `top_n=5` все три попадают в top_ids → их full-labels на Chart.js
    # рендерятся overlapping и нечитаемо (964 vs 969 разнесены на 4 кэВ).
    #
    # Fix: после выбора top_ids делаем pass с порогом MIN_DX_KEV: если
    # энергия пика < MIN_DX_KEV от уже принятой top-метки — отбрасываем
    # этот pid из top_ids (он остаётся в peaks для tooltip-only маркера,
    # т.е. line+dot всё ещё рисуется, но без full label). Порядок:
    # forced_top (BOOST_KINDS) идут ПЕРВЫМИ — образовательные SE/DE/sum
    # никогда не уступают; затем primary FEP по убыванию area_score.
    #
    # MIN_DX_KEV=15 выбрано как компромисс: Ac-228 911↔964 (Δ=53) keep both;
    # 964↔969 (Δ=4) → отбрасываем второй (слабее по I_γ). На zoom ~ 800-1000
    # кэВ метки занимают ~30 px = ~10-15 кэВ, поэтому 15 — минимум разделения.
    MIN_DX_KEV = 15.0
    if len(top_ids) > 1:
        peaks_by_id = {p.get("id", ""): p for p in peaks}

        def _label_score(pid: str) -> Tuple[float, float]:
            """Sort key for ordering top labels by priority.

            forced_top (BOOST features) — приоритет absolute (1e18);
            затем по area_by_id (выше area → выше приоритет — приходит
            раньше и резервирует slot). E_kev — secondary tiebreaker.
            """
            base = 1e18 if pid in forced_top else area_by_id.get(pid, 0.0)
            p = peaks_by_id.get(pid, {})
            return (-base, float(p.get("e") or 0.0))

        ordered = sorted(top_ids, key=_label_score)
        placed: List[float] = []
        survivors: set = set()
        for pid in ordered:
            p = peaks_by_id.get(pid)
            if p is None:
                continue
            e_kev = float(p.get("e") or 0.0)
            if any(abs(e_kev - e_p) < MIN_DX_KEV for e_p in placed):
                # Collision — этот label не получит full text (line+dot
                # остаются для tooltip), уступает уже принятому соседу.
                continue
            placed.append(e_kev)
            survivors.add(pid)
        top_ids = survivors

    for p in peaks:
        p["is_top"] = bool(p.get("id") in top_ids)
        # Cleanup intermediate field
        p.pop("_area_score", None)


# ──────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────

def _esc(x: Any) -> str:
    if x is None:
        return ""
    return _html.escape(str(x), quote=True)


def _round_e(e: Optional[float], digits: int = 1) -> float:
    if e is None:
        return 0.0
    try:
        return round(float(e), digits)
    except (TypeError, ValueError):
        return 0.0


def _safe_float(x) -> Optional[float]:
    if x is None:
        return None
    try:
        v = float(x)
        if math.isnan(v) or math.isinf(v):
            return None
        return v
    except (TypeError, ValueError):
        return None


def _chain_color(nuclide: str, label: str = "") -> str:
    if not nuclide and not label:
        return _COL_PHYS
    nl = (nuclide or "").strip()
    if nl in _TH_CHAIN_NUCLIDES or nl in _U_CHAIN_NUCLIDES:
        return _COL_TH
    if nl in _NATURAL_NUCLIDES:
        return _COL_NAT
    lab = (label or "").lower()
    if any(tok in lab for tok in (
        "комптон", "compton", "se ", "de ", "single esc", "double esc",
        "сумм", "sum peak", "pile", "k-ри", "xrf", "флуор", "anneal",
        "annihil", "511", "tcs", "intrinsic"
    )):
        return _COL_PHYS
    # Unknown — treat as physical artefact
    return _COL_PHYS


def _peak_id(E_keV: float) -> str:
    return f"p{int(round(E_keV))}"


def _load_template() -> str:
    """Load the interactive HTML template — prefer importlib.resources."""
    try:
        from importlib.resources import files  # py3.9+
        ref = files("gamma.reporting.templates").joinpath(
            "interactive_v1_17_2.html"
        )
        return ref.read_text(encoding="utf-8")
    except Exception:
        # Fallback: relative path next to this module
        here = Path(__file__).resolve().parent
        path = here / "templates" / "interactive_v1_17_2.html"
        return path.read_text(encoding="utf-8")


# ──────────────────────────────────────────────────────────────────
# Header block helpers
# ──────────────────────────────────────────────────────────────────

def _geometry_ru(canonical: str, fallback: str = "") -> str:
    key = (canonical or "").strip().lower()
    if key in _GEOMETRY_RU:
        return _GEOMETRY_RU[key]
    # Heuristic: substring match
    for k, v in _GEOMETRY_RU.items():
        if k in key:
            return v
    return canonical or fallback or "—"


def _build_title(header: Dict[str, Any]) -> str:
    fname = header.get("filename") or "образец"
    # Try to extract a nuclide hint from filename
    hints = header.get("filename_isotope_hints") or []
    nuclide = hints[0] if hints else None
    if not nuclide:
        # Try to derive from filename token like "Th232"
        m = re.search(r"\b([A-Za-z]{1,2})[\-_]?(\d{2,3})\b", fname)
        if m:
            nuclide = f"{m.group(1).capitalize()}-{m.group(2)}"
        else:
            nuclide = fname.rsplit(".", 1)[0][:32]
    geom = _geometry_ru(header.get("geometry_canonical") or "",
                       header.get("geometry") or "")
    return f"Гамма-спектр образца {nuclide}, геометрия {geom}"


def _build_subtitle(
    header: Dict[str, Any],
    diag: Dict[str, Any],
    f145: Optional[Dict[str, Any]] = None,
    energy_cal: Optional[Dict[str, Any]] = None,
) -> str:
    detector = header.get("detector_canonical") or header.get("detector_type") or "детектор"
    n_ch = header.get("n_channels") or 0
    live = float(header.get("live_time_s") or 0.0)
    live_h = live / 3600.0
    bg_status = header.get("background_status") or ""
    # F-325 / v1.18.18.1 — bg-state на графике должен быть однозначен.
    # На chart показывается то, что в spec.counts ПОСЛЕ всех преобразований
    # (т.е. если bg_subtracted=True → net spectrum, иначе gross).
    if bg_status == "auto_resolved_from_directory":
        # F-135 / v1.17.7 — фон автоматически найден и вычтен
        bg_txt = "график — net спектр (фон вычтен, авто-подбор F-131)"
    elif header.get("background_subtracted"):
        bg_txt = "график — net спектр (фон вычтен)"
    elif bg_status.startswith("embedded"):
        bg_txt = "график — gross спектр (фон встроен, не вычтен)"
    else:
        bg_txt = "график — gross спектр (без вычитания фона)"
    # F-144 / v1.17.7 — имена sample-файла и фон-файла в subtitle.
    sample_fn = header.get("sample_filename") or header.get("filename") or ""
    bg_fn = header.get("background_filename") or ""
    files_part = ""
    if sample_fn:
        files_part = f" · файл образца: {sample_fn}"
    if bg_fn:
        files_part += f" · файл фона: {bg_fn}"
    # G5 / v1.31.1 — провенанс энергетической калибровки (stored vs rebuilt).
    cal_part = ""
    if energy_cal:
        label = str(energy_cal.get("source_label") or "").strip()
        if label and label != "источник не указан":
            cal_part = f" · калибровка: {label}"
    # F-145 / v1.17.8 — пометка о пересчёте калибровки
    f145_part = ""
    if f145 and f145.get("phase_C_applied"):
        n_anch = f145.get("n_anchors_after_filter", 0)
        new_r = f145.get("new_residual_max_keV")
        old_r = f145.get("old_residual_max_keV")
        if new_r is not None and old_r is not None:
            f145_part = (
                f" · F-145: E(N) пересчитана по {n_anch} anchor'ам "
                f"(residual {old_r:.2f} → {new_r:.2f} кэВ)"
            )
    return (
        f"Детектор {detector} · {n_ch} каналов · "
        f"t_живое {live:.0f} с ({live_h:.2f} ч) · {bg_txt}"
        f"{files_part}"
        f"{cal_part}"
        f"{f145_part}"
    )


def _build_grid_cards(header: Dict[str, Any],
                      sample_mass_kg: Optional[float]) -> str:
    geom_ru = _geometry_ru(header.get("geometry_canonical") or "",
                          header.get("geometry") or "")
    mass_str = "—"
    if sample_mass_kg is not None:
        try:
            mass_str = f"{float(sample_mass_kg):.3f} кг"
        except Exception:
            pass
    shielding = "Pb 50 мм (свинцовая защита)"
    cards = [
        f'    <div><div class="k">Геометрия</div><div class="v">{_esc(geom_ru)}</div></div>',
        f'    <div><div class="k">Масса образца</div><div class="v">{_esc(mass_str)}</div></div>',
        f'    <div><div class="k">Защита</div><div class="v">{_esc(shielding)}</div></div>',
    ]
    return "\n".join(cards)


def _build_legend_items() -> str:
    items = [
        f'    <span><span class="fp-sw" style="background:{_COL_TH}"></span>Цепочка Th/U</span>',
        f'    <span><span class="fp-sw" style="background:{_COL_NAT}"></span>Натуральные / искусств.</span>',
        f'    <span><span class="fp-sw" style="background:{_COL_PHYS}"></span>Физ. артефакты</span>',
        # BUG-20 / v1.18.30+ (Agent B) — distinct orange swatch для secondary
        # artefacts (SE/DE/sum/511/backscatter/compton). Раньше эти пики
        # отображались parent's chain color и сливались с primary FEPs.
        f'    <span><span class="fp-sw" style="background:{_COL_SECONDARY}"></span>Вторичные процессы</span>',
    ]
    return "\n".join(items)


def _build_controls_info(header: Dict[str, Any]) -> str:
    detector = header.get("detector_canonical") or header.get("detector_type") or ""
    live = float(header.get("live_time_s") or 0.0)
    live_h = live / 3600.0
    return f"{detector} · {live_h:.2f} ч"


def _build_table_headers() -> str:
    # F-393 / v1.18.27 — sortable columns: data-sort hint (num|str) + title
    # tooltip; default sort = `line` asc (visual marker set client-side).
    # F-TBL-01 / v1.19.1 (Agent A, correction #5 2026-06-03):
    # column «T½» removed (it's a reference constant, not a per-peak diagnostic);
    # replaced with two per-peak diagnostic columns FWHM (kEV) and Z (σ).
    # FWHM source: primary_feps[].fwhm_keV (Mariscotti peak-fit). Z source:
    # primary_feps[].peak_significance_z (S/√B per Currie/ISO 11929) — distinct
    # from rate_sigma_cps (uncertainty of cps), которое идёт в Комментарий.
    return (
        '<tr>'
        '<th data-col="iso" data-sort="str" data-sortable="true"'
        ' title="Кликните для сортировки по нуклиду (А→Я)">Изотоп</th>'
        '<th data-col="line" data-sort="num" data-sortable="true"'
        ' title="Кликните для сортировки по энергии (по умолчанию)">Линия (изм/спр), кэВ</th>'
        '<th data-col="fwhm" data-sort="num" data-sortable="true"'
        ' title="Кликните для сортировки по ПШПВ (FWHM_measured); badge ⚠ если |ΔFWHM|/FWHM_cal > 10%">FWHM, кэВ</th>'
        '<th data-col="z" data-sort="num" data-sortable="true"'
        ' title="Кликните для сортировки по значимости пика Z = S/√B (Currie/ISO 11929)">Z, σ</th>'
        '<th data-col="a" data-sort="num" data-sortable="true"'
        ' title="Кликните для сортировки по активности">A, Бк/кг</th>'
        '<th data-col="cmt" data-sort="str" data-sortable="true"'
        ' title="Кликните для сортировки по комментарию">Комментарий</th>'
        '</tr>'
    )


# ──────────────────────────────────────────────────────────────────
# Data array builders
# ──────────────────────────────────────────────────────────────────

def _build_E_C(analysis_result) -> Tuple[List[float], List[int]]:
    """Sample spectrum into an E (keV) / C (counts) array (every 4 channels)."""
    spec = getattr(analysis_result, "spec", None)
    if spec is None or spec.counts is None:
        return [], []
    try:
        import numpy as np
        counts = np.asarray(spec.counts, dtype=np.int64)
    except Exception:
        counts = list(spec.counts)
    n = len(counts)
    # Energy axis via channel_to_energy
    E: List[float] = []
    C: List[int] = []
    # Stride to keep total points around ~1000 max
    stride = max(1, n // 1000)
    for ch in range(0, n, stride):
        try:
            e = spec.channel_to_energy(int(ch))
            if e is None:
                continue
            E.append(round(float(e), 1))
            C.append(int(counts[ch]))
        except Exception:
            continue
    return E, C


def _build_chart_payload(analysis_result) -> Dict[str, Any]:
    """F-332 / v1.18.18.5 — build payload для 4-way chart toggle.

    Возвращает dict с:
      * `E`         — энергетическая ось (кэВ), общая для всех 4 видов.
      * `C_net`     — net counts (отображается по умолчанию когда фон
                      вычтен; равно spec.counts ПОСЛЕ subtract_background).
      * `C_gross`   — gross sample counts (sample.counts до вычитания
                      фона); None если фон не применялся.
      * `C_bg`      — фон, ре-биннирован на сетку образца и масштабирован
                      по live-time (т.е. напрямую сопоставим с C_gross
                      для overlay/subtract). None если фон не применялся.
      * `has_background`  — bool, есть ли валидные C_gross + C_bg для
                           toggle (когда False — UI показывает только
                           «Образец» без переключателя).
      * `t_sample`, `t_bg`, `bg_scale` — метаданные для tooltip.

    Калибровка: оба массива выводятся на ОДНОЙ энергетической оси
    (sample's channel→keV mapping). Фон уже интерполирован на эту
    сетку в `subtract_background()` (F-58 / bg_subtract_energy.py:89).
    Live-time scaling также применён там (k = t_sample / t_bg, см.
    BackgroundSubtractionResult.scale_factor).
    """
    spec = getattr(analysis_result, "spec", None)
    if spec is None or spec.counts is None:
        return {"E": [], "C_net": [], "C_gross": None, "C_bg": None,
                "has_background": False}
    try:
        import numpy as np
        counts_net = np.asarray(spec.counts, dtype=np.float64)
    except Exception:
        counts_net = list(spec.counts)
        np = None

    n = len(counts_net)
    stride = max(1, n // 1000)

    gross = getattr(analysis_result, "gross_counts", None)
    bg_grid = getattr(analysis_result, "bg_counts_on_sample_grid", None)
    has_bg = (
        gross is not None
        and bg_grid is not None
        and len(gross) == n
        and len(bg_grid) == n
    )

    E: List[float] = []
    C_net: List[int] = []
    C_gross: List[float] = []
    C_bg: List[float] = []

    for ch in range(0, n, stride):
        try:
            e = spec.channel_to_energy(int(ch))
            if e is None:
                continue
            E.append(round(float(e), 1))
            C_net.append(int(round(float(counts_net[ch]))))
            if has_bg:
                C_gross.append(round(float(gross[ch]), 2))
                C_bg.append(round(float(bg_grid[ch]), 2))
        except Exception:
            continue

    # F-335.3 / v1.18.18.8 — pre-compute Savitzky-Golay smoothing для каждого
    # массива. Window=5 polyorder=2 (стандарт в гамма-спектрометрии, PyMca/
    # GADRAS). Default toggle OFF в JS — оператор включает по кнопке.
    # mode='nearest' для краёв (не дополнять нулями — иначе перекос на краях).
    def _savgol(arr: List[float]) -> List[float]:
        if not arr or len(arr) < 5:
            return list(arr) if arr else []
        try:
            from scipy.signal import savgol_filter
            import numpy as _np
            sm = savgol_filter(_np.asarray(arr, dtype=_np.float64),
                               window_length=5, polyorder=2, mode='nearest')
            return [round(float(v), 4) for v in sm]
        except Exception:
            return list(arr)

    payload: Dict[str, Any] = {
        "E": E,
        "C_net": C_net,
        "C_gross": C_gross if has_bg else None,
        "C_bg": C_bg if has_bg else None,
        "has_background": bool(has_bg),
        # F-334.2 / v1.18.18.7 — default Y units: cps (counts-per-second).
        # JS toggle setYUnits('cps'|'counts') реализуется в template.
        # t_sample используется для sample/clean/overlay; t_bg для bg view.
        "default_y_units": "cps",
        "t_sample": float(getattr(spec, "live_time", 0.0) or 0.0),
        # F-335.3 / v1.18.18.8 — smoothed arrays (Savitzky-Golay w=5 p=2)
        "C_net_smooth": _savgol(C_net),
        "smooth_method": "savitzky_golay",
        "smooth_window": 5,
        "smooth_polyorder": 2,
    }
    if has_bg:
        payload["t_bg"] = float(getattr(analysis_result, "bg_live_time", 0.0) or 0.0)
        payload["bg_scale"] = float(getattr(analysis_result, "bg_scale_factor", 0.0) or 0.0)
        # F-334.1 / v1.18.18.7 — peaks найденные на bg-спектре + matching
        # против static-dict для isotope label. Используется в setView('bg')
        # для подмены аннотаций (вместо sample-пиков, которые в bg-view
        # бессмысленны).
        payload["bg_peaks"] = _detect_bg_peaks(analysis_result, top_n=20)
        # F-335.3 / v1.18.18.8 — smoothed companion arrays
        payload["C_gross_smooth"] = _savgol(C_gross)
        payload["C_bg_smooth"] = _savgol(C_bg)
    else:
        payload["bg_peaks"] = []
        payload["C_gross_smooth"] = None
        payload["C_bg_smooth"] = None
    return payload


def _build_peaks(report: Dict[str, Any]) -> List[Dict[str, Any]]:
    """List of {id,e,label,color}, sorted by energy ascending.

    F-111 / F-111b (D-17, D-18): when Th-232 chain is dominant, library
    lines with I_γ ≥ 0.5% are injected as placeholders so every legend
    entry appears on the chart.
    """
    seen: Dict[str, Dict[str, Any]] = {}
    # Primary FEPs
    for fp in (report.get("primary_feps") or []):
        e = _safe_float(fp.get("peak_E_keV"))
        if e is None:
            continue
        # F-335.8 / v1.18.18.14 — drop zero-area FEPs (artefact attempts
        # типа Tl-208 510 которые collide с аннигиляцией 511 → area=0
        # после fit). Не показываем такие на чарте.
        area = _safe_float(fp.get("peak_area_counts")) or 0.0
        if area <= 0:
            continue
        nuc = fp.get("nuclide") or ""
        lib_e = _safe_float(fp.get("library_E_keV"))
        # F-335.8 — label transparency: показать «изм/спр» когда drift > 1.5 кэВ
        # (типичная NaI калибровка-нестабильность). Иначе оставляем компактную
        # форму «Pb-212 239».
        if lib_e is None:
            label = f"{nuc} {e:.0f}"
        elif abs(e - lib_e) > 1.5:
            label = f"{nuc} {e:.0f}/{lib_e:.0f}"
        else:
            label = f"{nuc} {lib_e:.0f}"
        pid = _peak_id(e)
        if pid in seen:
            continue
        seen[pid] = {
            "id": pid, "e": round(e, 1), "label": label,
            "color": _chain_color(nuc, label),
            # F-335.8 — feature_kind для _mark_top_peaks boost: primary FEPs всегда top по area
            "feature_kind": "primary_fep",
            "is_secondary": False,  # F-147 — primary FEPs всегда видимы
            "_area_score": float(area),
        }
    # Secondary peaks
    # F-147 / v1.18.22.0 — feature_kinds, которые на чарте получают
    # ОТДЕЛЬНЫЙ визуальный класс (точечный dash + collapsible toggle).
    # Эти процессы — артефакты физики детектора, не FEP'ы реальных γ-линий
    # нуклидов в образце; их полезно различать визуально и иметь
    # возможность скрыть.
    F147_SECONDARY_PHYS_KINDS = frozenset({
        "compton_edge", "backscatter", "annihilation_511",
        "single_escape", "double_escape", "sum_peak",
    })
    for sp in (report.get("secondary_peaks") or []):
        e = _safe_float(sp.get("energy_keV"))
        if e is None:
            continue
        pid = _peak_id(e)
        if pid in seen:
            continue
        # F-141 / v1.17.7 — feature_kind важнее type для labelling, потому
        # что computed features имеют type="расчётная_сигнатура" но точный
        # feature_kind (compton_edge/backscatter/single_escape/double_escape).
        kind = sp.get("feature_kind") or sp.get("type") or ""
        parent = sp.get("parent_nuclide") or ""
        parent_E = _safe_float(sp.get("parent_line_keV"))
        # F-335.8 — informative labels: include parent isotope + theoretical E
        # для SE/DE; явная привязка к цепочке для composite_cluster.
        parent_tag = f" {parent}" if parent else ""
        if kind == "annihilation_511":
            nice = "Аннигиляция 511"
        elif kind == "single_escape" and parent_E:
            nice = f"SE{parent_tag} {parent_E - 511:.0f}"
        elif kind == "single_escape":
            nice = f"SE{parent_tag} {int(round(e))}"
        elif kind == "double_escape" and parent_E:
            nice = f"DE{parent_tag} {parent_E - 1022:.0f}"
        elif kind == "double_escape":
            nice = f"DE{parent_tag} {int(round(e))}"
        elif kind == "sum_peak":
            nice = f"Сумм{parent_tag} {int(round(e))}"
        elif kind == "compton_edge":
            nice = f"K_C{parent_tag} {int(round(e))}"
        elif kind == "backscatter":
            nice = f"Обр.расс.{parent_tag} {int(round(e))}"
        elif kind in ("cluster", "composite_cluster"):
            # F-335.8 — для Th-232 кластера явно называем «Th-кластер»
            if parent == "Th-232" or "Th-232" in parent:
                nice = f"Th-кластер {int(round(e))} кэВ"
            else:
                nice = f"Кластер {int(round(e))}"
        elif kind == "расчётная_сигнатура":
            nice = f"Расчёт{parent_tag} {int(round(e))}"
        else:
            nice = kind or f"физ. {int(round(e))}"
        # BUG-20 / v1.18.30+ (Agent B) fix #1 — distinct orange for
        # secondary artefacts (SE/DE/sum/511/backscatter/compton_edge),
        # regardless of parent nuclide. Раньше parent's chain color делал
        # их визуально неотличимыми от primary FEPs Th/U цепочек.
        # cluster/composite_cluster — это группа РЕАЛЬНЫХ γ-линий, не
        # артефакт, оставляем им chain color.
        if kind in F147_SECONDARY_PHYS_KINDS:
            color = _COL_SECONDARY
        else:
            color = _chain_color(parent, nice) if parent else _COL_PHYS
        seen[pid] = {
            "id": pid, "e": round(e, 1), "label": nice,
            "color": color,
            # F-335.8 — feature_kind boost flag для top selection
            "feature_kind": kind,
            # F-147 / v1.18.22.0 — флаг для UI toggle «Вторичные процессы»:
            # позволяет JS отдельно рендерить annotation (dotted dash) и
            # массово прятать/показывать эти артефакты детектора.
            "is_secondary": kind in F147_SECONDARY_PHYS_KINDS,
            "_area_score": 0.0,
        }

    # F-111 / D-17: inject chain-completeness placeholders into peaks
    # F-335.8 / v1.18.18.14 — дедуп placeholder если та же (nuc, lib_E) уже
    # представлена primary FEP (в пределах ±10 кэВ для NaI calibration drift).
    # Иначе на чарте появляются двойные маркеры (Pb-212 235/239 + Pb-212 239).
    cd = (report.get("diagnostics") or {}).get("chain_dominance") or {}
    if cd.get("th232_dominant"):
        # Index existing primary FEPs by (nuclide, library_E_keV)
        fep_lib_set = set()
        for fp in (report.get("primary_feps") or []):
            nuc_fp = fp.get("nuclide") or ""
            lib_fp = _safe_float(fp.get("library_E_keV"))
            area_fp = _safe_float(fp.get("peak_area_counts")) or 0.0
            if nuc_fp and lib_fp is not None and area_fp > 0:
                fep_lib_set.add((nuc_fp, round(lib_fp, 1)))
        for nuc, e_lib, _I, _th in _TH232_CHAIN_LIBRARY:
            pid = _peak_id(e_lib)
            if pid in seen:
                continue
            # Skip if a primary FEP already covers this (nuc, lib_E)
            if (nuc, round(e_lib, 1)) in fep_lib_set:
                continue
            seen[pid] = {
                "id": pid, "e": round(e_lib, 1),
                "label": f"{nuc} {e_lib:.0f}",
                "color": _COL_TH,
                "feature_kind": "chain_completeness",
                "is_secondary": False,  # F-147 — chain placeholders — это γ-линии нуклидов, не артефакты
                "_area_score": 0.0,
            }

    peaks = sorted(seen.values(), key=lambda p: p["e"])
    return peaks


def _build_detail(report: Dict[str, Any],
                  peaks: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """For each peak id: {title, color, desc}."""
    detail: Dict[str, Dict[str, Any]] = {}
    # Index by id from primary feps
    fep_by_id: Dict[str, Dict[str, Any]] = {}
    for fp in (report.get("primary_feps") or []):
        e = _safe_float(fp.get("peak_E_keV"))
        if e is None:
            continue
        fep_by_id[_peak_id(e)] = fp
    sp_by_id: Dict[str, Dict[str, Any]] = {}
    for sp in (report.get("secondary_peaks") or []):
        e = _safe_float(sp.get("energy_keV"))
        if e is None:
            continue
        sp_by_id[_peak_id(e)] = sp

    # F-111 chain placeholder set (only used when Th dominant)
    cd = (report.get("diagnostics") or {}).get("chain_dominance") or {}
    chain_placeholder: Dict[str, str] = {}
    if cd.get("th232_dominant"):
        for nuc, e_lib, _I, _th in _TH232_CHAIN_LIBRARY:
            chain_placeholder[_peak_id(e_lib)] = nuc

    for p in peaks:
        pid = p["id"]
        title = p["label"]
        color = p["color"]
        if pid in fep_by_id:
            fp = fep_by_id[pid]
            nuc = fp.get("nuclide") or ""
            I_g = _safe_float(fp.get("library_I_pct"))
            area = _safe_float(fp.get("peak_area_counts"))
            sigma = _safe_float(fp.get("rate_sigma_cps")) or 0.0
            parts = []
            if nuc:
                parts.append(f"Линия {nuc}.")
            if I_g is not None:
                parts.append(f"Библиотечная интенсивность I_γ = {I_g:.2f}%.")
            if area is not None:
                parts.append(f"Площадь пика {area:.0f} отсч.")
            if not parts:
                parts.append("Идентифицированная характеристическая линия.")
            desc = " ".join(parts)
        elif pid in sp_by_id:
            sp = sp_by_id[pid]
            note = sp.get("note") or ""
            kind = sp.get("type") or ""
            parent = sp.get("parent_nuclide") or ""
            desc = note or f"Вторичный пик типа «{kind}»" + (
                f" от {parent}." if parent else "."
            )
        elif pid in chain_placeholder:
            desc = "Подтверждена по равновесию цепочки Th-232 (F-111)"
        else:
            desc = title
        detail[pid] = {"title": title, "color": color, "desc": desc}
    return detail


def _build_rows(report: Dict[str, Any],
                analysis_result,
                is_bg: bool = False) -> List[Dict[str, Any]]:
    """Table rows.

    BUG-19 / v1.18.30+ (Agent B): каждая строка получает `section`-метку для
    разделения таблицы на 3 подсекции в JS:
      • "primary_detected" — найденные FEP с матчем в библиотеке (есть peak_E)
      • "weak_candidate"   — слабые библиотечные кандидаты (если бы они были,
                              сейчас сюда же попадают совсем низко-I_γ FEP с
                              σ < 3 — но в текущей реализации ID-pipeline
                              такие отфильтровываются раньше, секция остаётся
                              в фреймворке для будущих расширений)
      • "secondary"        — backscatter / SE / DE / 511 / sum / compton_edge

    BUG-23 / v1.18.31+ (Agent B): убрана подсекция "chain_expected" — F-111
    placeholder-строки из цепочки (которых физически нет в спектре) больше
    НЕ добавляются в таблицу. Пользовательская обратная связь: «перечислять
    не найденные линии не нужно». Линии цепочки, реально подтверждённые в
    primary_feps, остаются в "primary_detected" — секции по равновесию более
    нет ни в Sample, ни в Bg представлении.

    Также: для cepочечных нуклидов (Ac-228, Tl-208, Pb-212, Bi-212, Bi-214,
    Pb-214) per-line активность `A_i ± σ_A_i` вытаскиваем из
    `analysis_result.activities[*].lines_used` (matching по library_E_keV),
    вместо общего weighted-mean. Weighted-mean остаётся в SUMMARY_CARD.
    """
    rows: List[Dict[str, Any]] = []
    # Identified nuclides → activity per nuclide (specific Bq/kg → fallback)
    act_by_nuc: Dict[str, Dict[str, Any]] = {}
    for n in (report.get("identified_nuclides") or []):
        act_by_nuc[n.get("nuclide") or ""] = n

    # BUG-19: per-line A_i ± σ_A_i lookup from analysis_result.activities.
    # Key: (nuclide, round(library_E_keV, 1)) → (A_Bq, sigma_A_Bq)
    line_act_by_key: Dict[Tuple[str, float], Tuple[float, float]] = {}
    # BUG-26: nuclides that have at least one LineActivity in analysis_result.
    # If a row's nuclide IS in this set but (nuc, lib_E) is NOT in
    # line_act_by_key, that row represents a line that compute.py
    # dedup-skipped (BUG-15) — drop it from the table. If the nuclide is
    # NOT in this set at all (no activity computed), keep the row and let
    # the weighted-mean fallback apply (BUG-19 backward-compat).
    nuclides_with_line_acts: set = set()
    sample_mass_kg: Optional[float] = None
    try:
        ars = getattr(analysis_result, "activities", None) or ()
        for ar in ars:
            nuc = getattr(ar, "nuclide", None)
            if not nuc:
                continue
            any_la_for_nuc = False
            for la in (getattr(ar, "lines_used", ()) or ()):
                ek = _safe_float(getattr(la, "E_keV", None))
                if ek is None:
                    continue
                a = _safe_float(getattr(la, "A_Bq", None))
                s = _safe_float(getattr(la, "sigma_A_Bq", None))
                if a is None:
                    continue
                line_act_by_key[(nuc, round(ek, 1))] = (a, s if s is not None else 0.0)
                any_la_for_nuc = True
            if any_la_for_nuc:
                nuclides_with_line_acts.add(nuc)
        sample_mass_kg = _safe_float(
            getattr(analysis_result, "sample_mass_kg", None)
        )
    except Exception:
        # Defensive: rows must still render even if analysis_result shape changes.
        line_act_by_key = {}
        nuclides_with_line_acts = set()
        sample_mass_kg = None

    def _fmt_line_activity(nuc: str, lib_e: float,
                           fallback_specific: str) -> str:
        """Per-line A_i if available, else fall back to nuclide weighted mean.

        Returns Bq/kg formatted string if mass is known; otherwise Bq.
        """
        la = line_act_by_key.get((nuc, round(lib_e, 1)))
        if la is None:
            return fallback_specific
        A, sA = la
        if sample_mass_kg and sample_mass_kg > 0:
            A_per_kg = A / sample_mass_kg
            sA_per_kg = sA / sample_mass_kg
            return f"{A_per_kg:.0f} ± {sA_per_kg:.0f}"
        return f"{A:.0f} ± {sA:.0f} Бк"

    # ─── BUG-26 / v1.18.31+ (Agent B) — row-layer dedup + phantom-zero filter ──
    # Mirrors BUG-15 logic from scripts/gamma/activity/compute.py at the
    # reporting layer. compute.py drops duplicate library lines (same
    # peak_channel, same nuclide → keep highest-I; same peak_channel,
    # cross-nuclide → keep characteristic owner) when computing per-line
    # activity. The reporting layer ignored this and rendered ALL matched
    # lines, producing phantom Ac-228/Tl-208 rows with A=0±σ or with the
    # nuclide weighted-mean A spread across every duplicate.
    #
    # Pass (a) — within-nuclide: for each (nuclide, peak_channel) group,
    # keep the entry with highest library_I_pct. Tie-break: highest
    # peak_area_counts, then smallest |library_E_keV − peak_E_keV|.
    # Pass (b) — cross-nuclide: build channel→characteristic-owners map.
    # For channels claimed by ≥2 nuclides where ≥1 is is_characteristic,
    # drop the non-owner entries.
    # Pass (c) — phantom-zero: rows with peak_area_counts == 0 (or missing)
    # provide no information and are dropped from the primary table.
    primary_feps_raw = list(report.get("primary_feps") or [])

    # Pass (c) — phantom-zero filter (S=0 or missing).
    def _S_of(fp: Dict[str, Any]) -> float:
        s = _safe_float(fp.get("peak_area_counts"))
        return float(s) if s is not None else 0.0

    primary_feps_filtered = [fp for fp in primary_feps_raw if _S_of(fp) > 0]

    # Pass (a) — within-nuclide dedup by (nuclide, peak_channel).
    by_nuc_ch: Dict[Tuple[str, int], List[Dict[str, Any]]] = {}
    for fp in primary_feps_filtered:
        nuc_key = fp.get("nuclide") or ""
        ch_raw = fp.get("peak_channel")
        ch_key = int(ch_raw) if ch_raw is not None else -1
        if ch_key < 0:
            # Without a channel we cannot dedup; keep as-is via unique sentinel.
            ch_key = -id(fp)
        by_nuc_ch.setdefault((nuc_key, ch_key), []).append(fp)
    keep_ids: set = set()
    for (nuc_key, ch_key), group in by_nuc_ch.items():
        if len(group) == 1:
            keep_ids.add(id(group[0]))
            continue
        def _within_key(fp: Dict[str, Any]) -> Tuple[float, float, float]:
            I = _safe_float(fp.get("library_I_pct")) or 0.0
            S = _safe_float(fp.get("peak_area_counts")) or 0.0
            lib_e = _safe_float(fp.get("library_E_keV")) or 0.0
            obs_e = _safe_float(fp.get("peak_E_keV")) or lib_e
            # Smallest |Δ| → highest negative for max() — invert sign.
            return (float(I), float(S), -abs(float(lib_e) - float(obs_e)))
        winner = max(group, key=_within_key)
        keep_ids.add(id(winner))
    primary_feps_a = [fp for fp in primary_feps_filtered if id(fp) in keep_ids]

    # Pass (b) — cross-nuclide dedup by characteristic owner.
    channel_owners: Dict[int, set] = {}
    channel_members: Dict[int, List[Dict[str, Any]]] = {}
    for fp in primary_feps_a:
        ch_raw = fp.get("peak_channel")
        if ch_raw is None:
            continue
        ch = int(ch_raw)
        channel_members.setdefault(ch, []).append(fp)
        if bool(fp.get("is_characteristic", False)):
            channel_owners.setdefault(ch, set()).add(fp.get("nuclide") or "")
    drop_ids: set = set()
    for ch, members in channel_members.items():
        nuclides_here = {m.get("nuclide") or "" for m in members}
        if len(nuclides_here) < 2:
            continue
        owners = channel_owners.get(ch, set())
        if not owners:
            continue
        for m in members:
            if (m.get("nuclide") or "") not in owners:
                drop_ids.add(id(m))
    primary_feps_dedup = [fp for fp in primary_feps_a if id(fp) not in drop_ids]

    seen_lines: set = set()
    for fp in primary_feps_dedup:
        e = _safe_float(fp.get("peak_E_keV"))
        if e is None:
            continue
        pid = _peak_id(e)
        nuc = fp.get("nuclide") or ""
        lib_e = _safe_float(fp.get("library_E_keV")) or e
        # BUG-26: if compute.py computed activity for this nuclide but did
        # NOT include this specific line in lines_used, it means the line
        # was dedup-skipped in compute_activity. Drop the row.
        if nuc in nuclides_with_line_acts:
            if (nuc, round(lib_e, 1)) not in line_act_by_key:
                continue
        thalf = _T_HALF_RU.get(nuc, "—")
        # Pill class
        if nuc in _TH_CHAIN_NUCLIDES or nuc in _U_CHAIN_NUCLIDES:
            kls = "fp-long"
        elif nuc in _NATURAL_NUCLIDES:
            kls = "fp-nat"
        else:
            kls = "fp-phys"
        line_str = f"{e:.1f} / {lib_e:.1f}"
        # F-TBL-01 / v1.19.1 (Agent A, correction #5 2026-06-03):
        # FWHM measured (kEV) — Mariscotti peak-fit output. Primary source is
        # `fwhm_keV`; fallback derive from `gauss_sigma_keV * 2.355`. If neither
        # available — show "—".
        fwhm_val = _safe_float(fp.get("fwhm_keV"))
        if fwhm_val is None:
            gs = _safe_float(fp.get("gauss_sigma_keV"))
            if gs is not None:
                fwhm_val = gs * 2.354820045
        fwhm_str = f"{fwhm_val:.1f}" if fwhm_val is not None else "—"
        # Z, σ — peak significance (S/√B per Currie / ISO 11929). Distinct from
        # rate_sigma_cps (count rate uncertainty, goes into Комментарий).
        # Pipeline emits this as `peak_significance_z`; fall back to alternate
        # naming if upstream changes.
        z_val = _safe_float(
            fp.get("peak_significance_z")
            or fp.get("peak_significance")
            or fp.get("significance_sigma")
            or fp.get("z_score")
        )
        z_str = f"{z_val:.1f}" if z_val is not None else "—"
        # Default: weighted-mean specific activity (fallback path).
        a_str_specific = "—"
        act = act_by_nuc.get(nuc)
        if act is not None:
            sa = _safe_float(act.get("specific_activity_Bq_per_kg"))
            sa_s = _safe_float(act.get("specific_activity_sigma_Bq_per_kg"))
            if act.get("is_upper_limit"):
                ul = _safe_float(act.get("upper_limit_Bq"))
                a_str_specific = f"< {ul:.1f}" if ul else "верхний предел"
            elif sa is not None:
                if sa_s is not None:
                    a_str_specific = f"{sa:.0f} ± {sa_s:.0f}"
                else:
                    a_str_specific = f"{sa:.0f}"
        # BUG-19: per-line activity for chain nuclides (Ac-228, Tl-208, etc.)
        a_str = _fmt_line_activity(nuc, lib_e, a_str_specific)
        if is_bg:
            a_str = "—"
        elif a_str.startswith("0 ±") or a_str == "0":
            continue  # skip phantom-zero activity rows (unconfirmed area)
        cmt_parts = []
        I_g = _safe_float(fp.get("library_I_pct"))
        if I_g is not None:
            cmt_parts.append(f"I_γ={I_g:.2f}%")
        sig = _safe_float(fp.get("rate_sigma_cps"))
        if sig is not None and sig > 0:
            cmt_parts.append(f"σ={sig:.3f} cps")
        cmt = ", ".join(cmt_parts) or "характеристическая линия"
        # G4 / v1.31.0 — «фон» pill for sample-side peaks that mirror a bg
        # identification. Annotation comes from json_report bg_carryover pass.
        # Rendered as <span class="fp-pill fp-bg"> via tmpl innerHTML.
        bg_co = fp.get("bg_carryover") or {}
        if bg_co.get("matched"):
            d_keV = bg_co.get("delta_E_keV")
            d_txt = f" (Δ={d_keV:.1f})" if isinstance(d_keV, (int, float)) else ""
            cmt = (
                f'<span class="fp-pill fp-bg" title="линия также найдена в '
                f'фоновом спектре{d_txt} — возможный остаток вычитания">'
                f'фон{d_txt}</span> {cmt}'
            )
        rows.append({
            "peak": pid, "iso": nuc, "kls": kls,
            "line": line_str,
            # F-TBL-01 / v1.19.1: `t` kept for backward-compat (sortable JS
            # may still address it); FWHM and Z added as new columns. The
            # template emits `fwhm` and `z` instead of `t` for visual cells.
            "t": thalf, "fwhm": fwhm_str, "z": z_str,
            "a": a_str, "cmt": cmt,
            "section": "primary_detected",
        })
        seen_lines.add((nuc, round(lib_e, 1)))

    # Secondary physical peaks
    for sp in (report.get("secondary_peaks") or []):
        e = _safe_float(sp.get("energy_keV"))
        if e is None:
            continue
        pid = _peak_id(e)
        kind = sp.get("type") or ""
        parent = sp.get("parent_nuclide") or "—"
        sig = _safe_float(sp.get("significance"))
        note = sp.get("note") or ""
        cmt = note or kind
        # F-TBL-01: secondary peaks rarely have measured FWHM (they are
        # placeholder XRF clusters / sum / SE / DE); Z is `significance`
        # field when present.
        sec_z_str = f"{sig:.1f}" if sig else "—"
        rows.append({
            "peak": pid,
            "iso": parent or "артефакт",
            "kls": "fp-phys",
            "line": f"{e:.1f} / —",
            "t": "—",
            "fwhm": "—",
            "z": sec_z_str,
            "a": "—",
            "cmt": cmt,
            "section": "secondary",
        })

    # BUG-23 / v1.18.31+ (Agent B): F-111 chain-completeness placeholder rows
    # удалены из таблицы. Раньше при th232_dominant / u238_dominant сюда
    # добавлялись строки с library_E без peak_E ("присутствует по цепочке"),
    # что пользователь расценил как шум — «перечислять не найденные линии не
    # нужно». Placeholder-маркеры на чарте (если есть) остаются под управлением
    # `_build_peaks`/`_build_detail` — bug касается только таблицы пиков.
    # `seen_lines` сохранён выше для возможных будущих weak_candidate logic.

    # D-03 — sort rows by ascending energy. We parse the energy out of
    # the `line` field; missing energies sort last.
    def _row_energy(r: Dict[str, Any]) -> float:
        s = r.get("line", "") or ""
        # `line` is like "238.6 / 238.6" or "— / 463.0" — pick the
        # first numeric value we see.
        m = re.search(r"\d+(?:\.\d+)?", s)
        if m:
            try:
                return float(m.group(0))
            except ValueError:
                return float("inf")
        return float("inf")

    rows.sort(key=_row_energy)
    return rows


# BUG-13 / v1.18.30+ (Agent B) — sync background-peak table with plot annotations.
# Symptom: на view «Фон» plot аннотирует K-40 1461 keV (label "K-40 1460.8"
# через _detect_bg_peaks → _match_bg_isotope против _BG_LINES_DICT), но в
# таблице «Пики фона» этой строки нет. _build_rows() для bg_report_view
# берёт `primary_feps` из `background_primary_feps`, которые проходят через
# pipeline chain-filtering — K-40 / Cs-137 / Pb K-XRF / U-235 как «не в
# доминирующей цепочке» туда не попадают, хотя их пики физически видны.
# Plot и таблица должны быть согласованы (одно — overlay другого).
#
# Fix: после `rows_bg = _build_rows(bg_report_view, ...)` вызываем эту
# функцию, которая зовёт ТУ ЖЕ `_detect_bg_peaks(analysis_result)` что и
# plot-аннотация (это единственный «оракул» истины для bg-вид) и добавляет
# в таблицу строки для каждого matched isotope, ещё не покрытого
# bg_report_view (по сопоставлению (isotope, round(E_lib, 1))).
_NATURAL_BG_T_HALF_RU: Dict[str, str] = {
    # K-40 + антропогенные + Pb K-XRF + U-235 → лебедь typical bg-only lines.
    "K-40":      "1.25×10⁹ лет",
    "Cs-137":    "30.08 года",
    "Pb K-XRF":  "—",                     # X-ray fluorescence — не radioactive decay
    "Pb Kα1":    "—",
    "Pb Kα2":    "—",
    "Pb Kβ":     "—",
    "U-235":     "7.04×10⁸ лет",
    "Аннигиляция": "—",                   # 511 keV e+e-
}


def _augment_bg_rows_with_natural_lines(
    rows_bg: List[Dict[str, Any]],
    analysis_result,
) -> List[Dict[str, Any]]:
    """BUG-13: добавить строки для K-40 / Cs-137 / Pb K-XRF / U-235 / 511,
    которые plot уже подписал, но pipeline в `background_primary_feps` не
    вернул (они не в Th/U цепочке, отсекаются на стадии identification).

    Возвращает НОВЫЙ список (не мутирует входной). Если detect_bg_peaks
    вернул пусто — возвращает rows_bg as-is (fail-safe для bg_grid=None
    случая).
    """
    try:
        bg_detected = _detect_bg_peaks(analysis_result, top_n=20, sigma_threshold=3.0)
    except Exception:
        # _detect_bg_peaks уже defensive, но дублируем чтобы row-builder не падал.
        return rows_bg
    if not bg_detected:
        return rows_bg

    # Index already-present rows by approximate (iso, lib_E).
    existing_keys: set = set()
    for r in rows_bg:
        iso = (r.get("iso") or "").strip()
        # rows используют формат "238.6 / 238.6" или "— / 463.0".
        line = r.get("line") or ""
        m = re.search(r"(\d+(?:\.\d+)?)\s*$", line)  # last number = lib_E
        if iso and m:
            try:
                lib_e = float(m.group(1))
                existing_keys.add((iso, round(lib_e, 0)))   # ±0.5 keV tolerance
            except ValueError:
                pass

    out = list(rows_bg)
    for d in bg_detected:
        iso = (d.get("isotope") or "").strip()
        if not iso:
            continue
        # BUG-43 / 2026-06-04 — Tl-208 override: если _detect_bg_peaks
        # переопределил 511 keV «Аннигиляция» → «Tl-208» (sentinel _511_tl208_override),
        # это не обычный Tl-208 char-line (583/2614) из primary_feps — это
        # особый случай 510.77 keV Tl-208 line. Пропускаем chain-filter
        # и добавляем строку явно с правильным pill-классом.
        is_511_tl208_override = bool(d.get("_511_tl208_override"))
        # Только natural-BG / pollution / X-ray isotopes — Th/U цепочка
        # уже идёт через primary_feps path (если в фоне есть Pb-212 — он там
        # уже есть). Это страховка от двойных строк для Pb-212 238.6 / Tl-208
        # 583.2 которые matched оба раза.
        # Исключение: Tl-208 511 override — это отдельная линия 510.77,
        # не дублирует primary_feps 583/2614.
        if (iso in _TH_CHAIN_NUCLIDES or iso in _U_CHAIN_NUCLIDES) and not is_511_tl208_override:
            continue
        E_lib = float(d.get("E_lib") or d.get("e") or 0.0)
        E_meas = float(d.get("e") or 0.0)
        if (iso, round(E_lib, 0)) in existing_keys:
            continue
        # Build natural row matching _build_rows output shape exactly.
        if is_511_tl208_override:
            # Tl-208 510.77 — chain nuclide pill (fp-long), not phys artefact
            kls = "fp-long"
        elif iso in _NATURAL_NUCLIDES:
            kls = "fp-nat"
        else:
            # Pb K-XRF, U-235, Аннигиляция → phys (artefact/non-chain)
            kls = "fp-phys"
        # BUG-43 — для 511→Tl-208 override используем _T_HALF_RU (цепочка),
        # для остальных — _NATURAL_BG_T_HALF_RU (K-40, Cs-137, XRF, U-235).
        if is_511_tl208_override:
            thalf = _T_HALF_RU.get(iso, "—")
        else:
            thalf = _NATURAL_BG_T_HALF_RU.get(iso, "—")
        intensity = float(d.get("intensity") or 0.0)
        cmt_parts = [f"intensity={intensity:.0f}"]
        if iso == "Аннигиляция":
            cmt_parts.append("e⁺e⁻ → 511 кэВ (космика + пара)")
        elif is_511_tl208_override:
            cmt_parts.append("Tl-208 510.77 кэВ (Th-232 chain); + возможен вклад от аннигиляции (e⁺e⁻)")
        elif iso.startswith("Pb K"):
            cmt_parts.append("K-флуоресценция Pb-защиты")
        elif iso == "K-40":
            cmt_parts.append("естественный калий (доминирует в bg)")
        elif iso == "U-235":
            cmt_parts.append("обычно примесь в Pb-защите")
        elif iso == "Cs-137":
            cmt_parts.append("антропогенное загрязнение")
        cmt = ", ".join(cmt_parts)
        out.append({
            "peak": d.get("id") or f"bg{int(round(E_meas))}",
            "iso": iso,
            "kls": kls,
            "line": f"{E_meas:.1f} / {E_lib:.1f}",
            "t": thalf,
            "a": "—",
            "cmt": cmt,
            # BUG-19 / BUG-13 — natural bg lines = primary detected (они РЕАЛЬНО
            # видны в фоне, plot их подписал) — не chain_expected.
            "section": "primary_detected",
        })
        existing_keys.add((iso, round(E_lib, 0)))

    # Re-sort by ascending energy (same key as _build_rows).
    def _row_energy(r: Dict[str, Any]) -> float:
        s = r.get("line", "") or ""
        m = re.search(r"\d+(?:\.\d+)?", s)
        if m:
            try:
                return float(m.group(0))
            except ValueError:
                return float("inf")
        return float("inf")

    out.sort(key=_row_energy)
    return out


def _sync_peaks_rows_detail(peaks: List[Dict[str, Any]],
                             rows: List[Dict[str, Any]],
                             detail: Dict[str, Any]) -> None:
    """F-111b: enforce set(peaks ids) == set(rows[].peak) == set(detail keys).

    Mutates the three structures in place so every id appears in all three.
    Orphan peaks gain a chain-completeness row + detail entry; orphan rows
    gain a peak entry; orphan detail entries gain peak + row.
    """
    pk_ids = {p.get("id") for p in peaks if p.get("id")}
    rw_ids = {r.get("peak") for r in rows if r.get("peak")}
    dt_ids = set(detail.keys())
    all_ids = pk_ids | rw_ids | dt_ids

    pk_by_id = {p["id"]: p for p in peaks if p.get("id")}
    rw_by_id = {r["peak"]: r for r in rows if r.get("peak")}

    for pid in all_ids:
        # Ensure peak exists
        if pid not in pk_ids:
            r = rw_by_id.get(pid)
            d = detail.get(pid)
            e_kev = 0.0
            label = pid
            color = _COL_PHYS
            if r:
                m = re.search(r"\d+(?:\.\d+)?", r.get("line", ""))
                if m:
                    try:
                        e_kev = float(m.group(0))
                    except ValueError:
                        pass
                label = f'{r.get("iso", "—")} {e_kev:g}'
                color = _chain_color(r.get("iso", ""), r.get("iso", ""))
            elif d:
                label = d.get("title", pid)
                color = d.get("color", _COL_PHYS)
                m = re.search(r"\d+(?:\.\d+)?", label)
                if m:
                    try:
                        e_kev = float(m.group(0))
                    except ValueError:
                        pass
            peaks.append({"id": pid, "e": e_kev, "label": label, "color": color})

        # Ensure row exists
        if pid not in rw_ids:
            p = pk_by_id.get(pid) or next((q for q in peaks if q.get("id") == pid), None)
            d = detail.get(pid)
            # BUG-23 fix: skip chain_completeness placeholders only (no real FEP).
            # For all other orphan peaks (double_escape, annihilation, secondary etc.)
            # create a row to preserve the peaks<->rows<->detail invariant (F-111b).
            if (p or {}).get("feature_kind") == "chain_completeness":
                pass  # chain_completeness: chart marker only, no table row
            else:
                e_kev = float((p or {}).get("e", 0.0) or 0.0)
                label = (p or {}).get("label") or (d or {}).get("title") or pid
                rows.append({
                    "peak": pid,
                    "iso": label,
                    "kls": "fp-phys",
                    "line": f"{e_kev:.1f} / вЂ”",
                    "t": "вЂ”",
                    "fwhm": "вЂ”",
                    "z": "вЂ”",
                    "a": "вЂ”",
                    "cmt": (d or {}).get("desc", ""),
                    "section": "secondary",
                })


        # Ensure detail exists
        if pid not in dt_ids:
            p = pk_by_id.get(pid) or next((q for q in peaks if q.get("id") == pid), None)
            r = rw_by_id.get(pid) or next((q for q in rows if q.get("peak") == pid), None)
            title = (p or {}).get("label") or (r or {}).get("iso") or pid
            color = (p or {}).get("color", _COL_PHYS)
            desc = ((r or {}).get("cmt")
                    or "Линия указана для согласованности peaks↔rows↔detail (F-111b).")
            detail[pid] = {"title": title, "color": color, "desc": desc}

    # Resort peaks by energy
    peaks.sort(key=lambda p: float(p.get("e", 0.0) or 0.0))
    # Resort rows by energy
    def _row_e(r):
        m = re.search(r"\d+(?:\.\d+)?", r.get("line", "") or "")
        try:
            return float(m.group(0)) if m else float("inf")
        except ValueError:
            return float("inf")
    rows.sort(key=_row_e)


def _gauss(x: float, E0: float, area: float, sigma: float) -> float:
    if sigma <= 0:
        return 0.0
    return (area / (sigma * math.sqrt(2.0 * math.pi))
            * math.exp(-((x - E0) ** 2) / (2.0 * sigma * sigma)))


def _build_multiplets_data(report: Dict[str, Any],
                           analysis_result) -> List[Dict[str, Any]]:
    """Build the multiplets[] JS array.

    Required keys per multiplet:
      id, title, chi2_per_dof, closure_pct, roi_low_ch, roi_high_ch,
      n_channels, E_keV[], data[], continuum[], total[], components[].
    Each component: {nuclide, E_keV, I_pct, area, g_plus_cont[], g_base[]}.
    """
    decons = analysis_result.deconvolution_results or []
    if not decons:
        return []

    spec = analysis_result.spec
    try:
        import numpy as np
        counts = np.asarray(spec.counts, dtype=np.float64) if spec.counts is not None else None
    except Exception:
        return []
    if counts is None or counts.size == 0:
        return []

    # FWHM model — F-452: polymorphic FwhmModel callable | legacy 3-tuple.
    _fwhm_model_local = analysis_result.fwhm_model
    def fwhm_kev(E: float) -> float:
        return _fwhm_keV_at_energy(_fwhm_model_local, float(E))

    out = []
    for i, d in enumerate(decons, start=1):
        roi_lo = int(getattr(d, "roi_low_ch", 0))
        roi_hi = int(getattr(d, "roi_high_ch", 0))
        if not (0 <= roi_lo < roi_hi <= len(counts)):
            continue

        # F-134 / v1.17.7 — overlay-данные из CoupledFitResult.
        # Если F-117 связанный fit заполнил overlay_* поля, рендер берёт
        # их точно. Это устраняет рассинхрон с настоящим fit'ом, где
        # форма Гаусс+tail+step не реконструируется канальными формулами.
        # F-373 — теперь legacy deconvolve_multiplet тоже заполняет
        # overlay_data/continuum/total/components (в channel space). Когда
        # overlay_E_keV отсутствует, мы строим ось energy из roi-каналов
        # через spec.channel_to_energy. Это исправляет рендер step_linear
        # континуума для free-NNLS подгонок (раньше шёл legacy fallback
        # с неверной формулой cp[0]+cp[1]·ch без step-члена).
        has_overlay = bool(
            getattr(d, "overlay_data", None)
            and getattr(d, "overlay_continuum", None)
            and getattr(d, "overlay_total", None)
        )

        if has_overlay:
            data_arr = [float(v) for v in d.overlay_data]
            cont_arr = [float(v) for v in d.overlay_continuum]
            total_arr = [float(v) for v in d.overlay_total]
            n_pts = len(data_arr)
            if getattr(d, "overlay_E_keV", None):
                E_arr = [round(float(v), 2) for v in d.overlay_E_keV]
            else:
                # F-373 — derive E_keV from spec at ROI channels
                E_arr = [
                    round(float(spec.channel_to_energy(int(ch))), 2)
                    for ch in range(roi_lo, roi_lo + n_pts)
                ]
            # Per-component g_plus_cont — прямо из overlay.
            # g_base = continuum (база для рендера компоненты).
            overlays = d.overlay_components or ()
            comp_payload = []
            for k, (comp, area) in enumerate(
                zip(getattr(d, "components", []) or [],
                    getattr(d, "areas", []) or [])
            ):
                nuc = getattr(comp, "nuclide", "") or ""
                E_line = float(getattr(comp, "line_E_keV", 0.0) or 0.0)
                I_pct = float(getattr(comp, "library_I_pct", 0.0) or 0.0)
                area_f = float(area or 0.0)
                if k < len(overlays):
                    g_plus_cont = [float(v) for v in overlays[k]]
                else:
                    g_plus_cont = list(cont_arr)
                g_base = list(cont_arr)
                # BUG-24 / v1.18.31+ (Agent B) — если overlay-компонент
                # отсутствует или вырожден (g_plus_cont ≈ cont_arr) при
                # значимой площади S > 0, синтезируем Гауссиану из
                # FWHM-модели чтобы заливка была видна. Иначе компонент
                # числится в легенде (`comp_payload` всегда полный), но
                # на холсте — пусто (peakY - continuum = 0). Это и есть
                # симптом «вертикальная линия + label есть, заливки нет».
                if area_f > 0:
                    delta_max = max(
                        (g - b for g, b in zip(g_plus_cont, g_base)),
                        default=0.0,
                    )
                    if delta_max <= 1e-9:
                        sigma_e = fwhm_kev(E_line) / 2.355
                        g_plus_cont = [
                            g_base[j] + _gauss(E_arr[j], E_line, area_f, sigma_e)
                            for j in range(n_pts)
                        ]
                comp_payload.append({
                    "nuclide": nuc,
                    "E_keV": round(E_line, 2),
                    "I_pct": round(I_pct, 3),
                    "area": round(area_f, 1),
                    "g_plus_cont": [round(v, 3) for v in g_plus_cont],
                    "g_base": [round(v, 3) for v in g_base],
                })
            sum_data = sum(data_arr) or 1.0
            sum_total = sum(total_arr)
            closure_pct = (sum_total - sum_data) / sum_data * 100.0
        else:
            # Legacy fallback (свободный NNLS-fit без overlay): канальный
            # рендер. Используется только для свободно подобранных
            # одиночных пиков из apply_multiplet_deconvolution.
            n_pts = roi_hi - roi_lo
            E_arr = []
            data_arr = []
            for ch in range(roi_lo, roi_hi):
                e = spec.channel_to_energy(int(ch))
                E_arr.append(round(float(e), 2))
                data_arr.append(float(counts[ch]))
            cp = list(getattr(d, "continuum_params", ()) or ())
            # F-373 — match deconvolve_multiplet's model: B(ch) = β₀ + β₁·(ch - x_mid)
            # [+ β_step·0.5·erfc((ch - x_step)/(σ_step·√2))]. Old formula
            # B(ch) = cp[0] + cp[1]·ch was wrong (raw ch instead of (ch -
            # x_mid)) and dropped the step term entirely on step_linear,
            # producing artefacts like a continuum that climbs through the
            # data or sinks below zero. Clamp ≥0 (counts can't be negative).
            cmodel = str(getattr(d, "continuum_model", "linear") or "linear")
            x_mid = 0.5 * (roi_lo + roi_hi - 1)
            comps = getattr(d, "components", ()) or ()
            if comps:
                fwhms = [float(getattr(c, "fwhm_channels", 1.0) or 1.0) for c in comps]
                weights = [max(w, 1e-6) for w in fwhms]
                centers = [float(getattr(c, "center_channel", 0.0) or 0.0) for c in comps]
                if sum(weights) > 0:
                    x_step = sum(cc * w for cc, w in zip(centers, weights)) / sum(weights)
                else:
                    x_step = x_mid
                sigma_step_ch = max(fwhms) / 2.355
            else:
                x_step = x_mid
                sigma_step_ch = 1.0
            import math as _math
            cont_arr = []
            for j, ch in enumerate(range(roi_lo, roi_hi)):
                c0 = cp[0] if len(cp) > 0 else 0.0
                c1 = cp[1] if len(cp) > 1 else 0.0
                val = float(c0) + float(c1) * (ch - x_mid)
                if cmodel == "step_linear" and len(cp) >= 3:
                    # 0.5·erfc((ch - x_step)/(σ_step·√2))
                    arg = (ch - x_step) / (sigma_step_ch * _math.sqrt(2.0))
                    val += float(cp[2]) * 0.5 * _math.erfc(arg)
                # F-373 clamp: continuum must be ≥ 0
                cont_arr.append(max(val, 0.0))
            comp_gauss: List[List[float]] = []
            for comp, area in zip(getattr(d, "components", []) or [],
                                  getattr(d, "areas", []) or []):
                E_line = float(getattr(comp, "line_E_keV", 0.0) or 0.0)
                fwhm = fwhm_kev(E_line)
                sigma = fwhm / 2.355
                curve = [_gauss(e, E_line, float(area or 0.0), sigma)
                         for e in E_arr]
                comp_gauss.append(curve)
            total_arr = [cont_arr[j] + sum(cg[j] for cg in comp_gauss)
                         for j in range(n_pts)]
            sum_data = sum(data_arr) or 1.0
            sum_total = sum(total_arr)
            closure_pct = (sum_total - sum_data) / sum_data * 100.0
            comp_payload = []
            for k, (comp, area) in enumerate(
                zip(getattr(d, "components", []) or [],
                    getattr(d, "areas", []) or [])
            ):
                nuc = getattr(comp, "nuclide", "") or ""
                E_line = float(getattr(comp, "line_E_keV", 0.0) or 0.0)
                I_pct = float(getattr(comp, "library_I_pct", 0.0) or 0.0)
                this_curve = comp_gauss[k]
                others = [sum(comp_gauss[kk][j]
                              for kk in range(len(comp_gauss))
                              if kk != k) for j in range(n_pts)]
                g_base = [cont_arr[j] + others[j] for j in range(n_pts)]
                g_plus_cont = [g_base[j] + this_curve[j] for j in range(n_pts)]
                comp_payload.append({
                    "nuclide": nuc,
                    "E_keV": round(E_line, 2),
                    "I_pct": round(I_pct, 3),
                    "area": round(float(area or 0.0), 1),
                    "g_plus_cont": [round(v, 3) for v in g_plus_cont],
                    "g_base": [round(v, 3) for v in g_base],
                })

        # Title — list nuclides + energies
        title_parts = []
        for comp in getattr(d, "components", []) or []:
            nuc = getattr(comp, "nuclide", "") or ""
            E_line = float(getattr(comp, "line_E_keV", 0.0) or 0.0)
            title_parts.append(f"{nuc} {E_line:.1f}")
        title = (f"Мультиплет M{i} — связанная подгонка "
                 + " + ".join(title_parts) + " кэВ")

        # F-145 / v1.17.8 — Phase A χ² для отображения рядом с locked χ²
        pA_chi2_val = getattr(d, "phase_A_chi2_per_dof", None)

        # F-378 / v1.18.25.1 — mark top-K компонент по площади для
        # selective annotation. Маркёры на холсте рисуются только для
        # is_top компонент (≤TOP_K=3 максимальных по deconvolved_area);
        # остальные компоненты остаются в datasets (Chart.js покажет их
        # в hover-tooltip как обычно). Без этого M3 V2 с 7 компонентами
        # превращался в визуальную кашу — 5 «degenerate» Ac-228 маркёров
        # перекрывались в одном окне.
        TOP_K = 3
        if comp_payload:
            ranked = sorted(
                enumerate(comp_payload),
                key=lambda kv: abs(float(kv[1].get("area") or 0.0)),
                reverse=True,
            )
            top_idx = {k for k, _ in ranked[:TOP_K]}
            for k, cp in enumerate(comp_payload):
                cp["is_top"] = bool(k in top_idx)

        out.append({
            "id": f"M{i}c",
            "cluster_id": str(getattr(d, "cluster_id", "") or ""),
            "title": title,
            "chi2_per_dof": round(float(getattr(d, "chi2_per_dof", 0.0) or 0.0), 3),
            "closure_pct": round(closure_pct, 2),
            "roi_low_ch": roi_lo,
            "roi_high_ch": roi_hi,
            "n_channels": n_pts,
            "E_keV": E_arr,
            "data": data_arr,
            "continuum": [round(v, 3) for v in cont_arr],
            "total": [round(v, 3) for v in total_arr],
            "components": comp_payload,
            # F-145 опциональные поля для caption блока
            "phase_A_chi2_per_dof": (
                round(float(pA_chi2_val), 2) if pA_chi2_val is not None else None
            ),
            # F-392 / v1.18.27 + F-392.1 / v1.18.28 — выбранная continuum
            # модель: "linear" | "step_linear" | "step_linear_multi" |
            # "quadratic". Используется в meta-бейдже для отличия multi-
            # anchor step от обычного step (см. _build_multiplet_blocks).
            "continuum_model": str(getattr(d, "continuum_model", "") or ""),
        })
    return out


# BUG-25 / v1.18.31+ (Agent B) — пороги фильтрации вырожденных
# мультиплетов в HTML/JS render слое. Применяются ПОСЛЕ деконволюции
# (`coupled_multiplet.py` владеет fit-результатами и НЕ трогается):
#
#   • MIN_COMPONENT_SNR — минимальный S/σ_S для «значимого» компонента;
#     при S=0 (явный flat NNLS solution) σ_S обычно тоже 0 → автоматом
#     не значимо. При σ_S=0 и S>0 — считаем значимым (legacy fit без
#     uncertainty).
#   • MAX_ACCEPTABLE_CHI2_DOF — катастрофические fit'ы (~ много σ выше
#     ожидаемого ~1-5) указывают на отсутствие реального сигнала или
#     неправильный континуум.
#
# Эмпирически подобрано: Th-232 demo M1 (χ²/ν≈29), M3 (χ²/ν≈931 с
# Pb-212 238 S=570k), M4 (χ²/ν≈180 с Tl-208 583 S=205k) ВСЕ имеют
# хотя бы один компонент с S>0 — проходят. Чисто-noise мультиплет с
# всеми S=0 и χ²/ν~500-1000+ отсекается.
MIN_COMPONENT_SNR = 3.0
MAX_ACCEPTABLE_CHI2_DOF = 1000.0


def _is_meaningful_multiplet(mp: Dict[str, Any]) -> bool:
    """BUG-25 / v1.18.31+ (Agent B) — гейтер вырожденных мультиплетов.

    Returns False если ни один компонент не несёт значимого сигнала
    ИЛИ если fit катастрофически плох (χ²/ν выше потолка). Применяется
    ТОЛЬКО в reporting layer — в JSON multiplet_deconvolutions всегда
    остаётся полный список для аудита.

    Контракт «meaningful»:
      • Хотя бы один компонент с S > 0 (NNLS-solved, не fixed-zero)
        ИЛИ значимый по S/σ_S ≥ MIN_COMPONENT_SNR;
      • chi2_per_dof < MAX_ACCEPTABLE_CHI2_DOF.

    Если ВСЕ компоненты дают S=0 — fit сошёлся к нулевому вкладу пиков
    (только континуум) — для оператора это шумовое окно, не реальный
    мультиплет. Отображать его в HTML отчёте только запутывает.
    """
    chi2 = float(mp.get("chi2_per_dof") or 0.0)
    if chi2 >= MAX_ACCEPTABLE_CHI2_DOF:
        return False
    has_signal = False
    for c in mp.get("components", []) or []:
        S = float(c.get("area") or 0.0)
        if S <= 0:
            continue
        sigma_S = c.get("area_sigma")  # optional; not always present in payload
        try:
            sigma_f = float(sigma_S) if sigma_S is not None else 0.0
        except (TypeError, ValueError):
            sigma_f = 0.0
        # Если σ_S не задан (или 0) — S>0 считаем значимым.
        # Если σ_S>0 — требуем S/σ ≥ MIN_COMPONENT_SNR.
        if sigma_f <= 0 or (S / sigma_f) >= MIN_COMPONENT_SNR:
            has_signal = True
            break
    return has_signal


def _filter_meaningful_multiplets(
    multiplets: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """BUG-25 helper: применить `_is_meaningful_multiplet` и логировать
    отсев. Возвращает новый список (исходный не мутируется)."""
    kept: List[Dict[str, Any]] = []
    for mp in multiplets:
        if _is_meaningful_multiplet(mp):
            kept.append(mp)
        else:
            log.debug(
                "BUG-25: multiplet %s skipped "
                "(all components S=0 or chi2/dof=%.2f above %.1f cap)",
                mp.get("id", "?"),
                float(mp.get("chi2_per_dof") or 0.0),
                MAX_ACCEPTABLE_CHI2_DOF,
            )
    return kept


def _build_multiplet_blocks(
    multiplets: List[Dict[str, Any]],
    *,
    bg: bool = False,
) -> str:
    """Render HTML blocks for each multiplet (before the summary).

    F-397.2 / v1.18.28 (Agent B): when ``bg=True``, suffix every canvas id
    with ``-bg`` AND suffix the section heading с пометкой "(фон)", чтобы
    при печати/без JS-toggle оба блока были визуально различимы. Sample-
    блок (default bg=False) сохраняет старый id-namespace для обратной
    совместимости с тестами и canvas-таргетингом в JS.

    BUG-5 / v1.18.30+ (Agent B): два набора мультиплетов (sample vs bg)
    подгоняются на разных спектрах — это НЕ дублирование, а разные
    исходные данные:
      • sample → анализ образца (S = площади в спектре образца)
      • bg     → анализ референсного фонового спектра (S = площади в фоне)
    Заголовки явно разделяют роли, чтобы оператор при печати/PDF/копировании
    не считал два блока одним и тем же расчётом.
    """
    if bg:
        title_label = (
            "Мультиплеты — разложение в фоновом спектре "
            "(референс для вычитания, связанная подгонка по библиотечным "
            "интенсивностям)"
        )
        subhead = (
            'Подгонка применена к спектру фонового измерения отдельно от '
            'образца — приведена для контроля чистоты фона. Площади '
            '<b>S</b> ниже относятся к фоновому спектру, не к образцу.'
        )
    else:
        title_label = (
            "Мультиплеты — разложение в спектре образца "
            "(первичная подгонка, связанная по библиотечным интенсивностям)"
        )
        subhead = (
            'Площади <b>S</b> ниже взяты из спектра образца — это '
            'основной источник характеристик пиков для расчёта активности.'
        )

    if not multiplets:
        empty_title = (
            "Мультиплеты — разложение в фоновом спектре"
            if bg else
            "Мультиплеты — разложение в спектре образца"
        )
        return (
            f'  <h2 style="font-size:17px;font-weight:500;margin:28px 0 10px;">'
            f'{empty_title}</h2>\n'
            '  <div class="fp-mp-block">\n'
            '    <p class="fp-mp-note" style="color:var(--text-tertiary);font-style:italic;">'
            'Мультиплеты не разрешались в этом прогоне '
            '(нет перекрывающихся характеристических линий выше порога).'
            '</p>\n'
            '  </div>'
        )
    parts = [
        f'  <h2 style="font-size:17px;font-weight:500;margin:28px 0 10px;">'
        f'{title_label}</h2>',
        f'  <p class="fp-mp-section-note" '
        f'style="margin:0 0 12px;font-size:12.5px;color:var(--text-secondary);">'
        f'{subhead}</p>'
    ]
    id_suffix = "-bg" if bg else ""
    for m in multiplets:
        title = _esc(m.get("title") or "Мультиплет")
        chi = m.get("chi2_per_dof", 0.0)
        n_ch = m.get("n_channels", 0)
        closure = m.get("closure_pct", 0.0)
        # F-145 / v1.17.8: Phase A χ² для caption — если был запущен
        pA_chi2 = m.get("phase_A_chi2_per_dof")
        pA_part = ""
        if pA_chi2 is not None:
            pA_part = f" · фаза А (F-145) χ²/ν = {pA_chi2:.2f}"
        # F-392 / v1.18.27 + F-392.1 / v1.18.28 — badge зависит от выбранной
        # continuum модели. Для step_linear_multi показываем расширенный
        # бейдж (multi-anchor step). Для legacy / unknown — старый текст
        # "ступенчато-линейный континуум" (back-compat).
        cmodel = str(m.get("continuum_model") or "")
        if cmodel == "step_linear_multi":
            cont_badge = (
                '<span class="fp-mp-cmodel fp-mp-cmodel-multi" '
                'style="color:var(--accent,#1565c0);font-weight:500;">'
                'ступ.-лин. континуум (multi-anchor, F-392)'
                '</span>'
            )
        elif cmodel == "step_linear":
            cont_badge = "ступенчато-линейный континуум"
        elif cmodel == "quadratic":
            cont_badge = "квадратичный континуум"
        elif cmodel == "linear":
            cont_badge = "линейный континуум"
        else:
            cont_badge = "ступенчато-линейный континуум"
        # Compose component summary.
        # BUG-24 / v1.18.31+ (Agent B) — оставляем в текстовой легенде
        # только компоненты с area > 0. Это симметрично JS-рендеру
        # (interactive_v1_17_2.html:1366 пропускает area === 0 при
        # построении datasets и lineAnnots). Раньше HTML-нота показывала
        # ВСЕ компоненты (включая S=0 phantom-fit), а на холсте их не было
        # → оператор не находил заявленный «компонент» на графике.
        comp_strs = []
        for c in m.get("components", []):
            area = c.get("area") or 0.0
            if area <= 0:
                continue
            comp_strs.append(
                f"<b>{_esc(c.get('nuclide',''))} {c.get('E_keV',0):.1f}</b> — "
                f"I={c.get('I_pct',0):.2f}%, S={area:.0f}"
            )
        note = "; ".join(comp_strs) or "—"
        canvas_id = f'mp-{_esc(m.get("id","M"))}{id_suffix}'
        parts.append(
            f'  <div class="fp-mp-block">\n'
            f'    <h3>{title}</h3>\n'
            f'    <div class="fp-mp-meta">'
            f'χ²/ν = {chi:.2f} · ROI {n_ch} кан · '
            f'{cont_badge} · '
            f'закрытие Δ = {closure:.1f}%'
            f'{pA_part}'
            f'</div>\n'
            f'    <div class="fp-chart-mp">'
            f'<canvas id="{canvas_id}" width="860" height="340"></canvas>'
            f'</div>\n'
            f'    <p class="fp-mp-note">{note}.</p>\n'
            f'  </div>'
        )
    return "\n".join(parts)


# ──────────────────────────────────────────────────────────────────
# Summary / notes
# ──────────────────────────────────────────────────────────────────

# F-368 / v1.18.24.3 — natural-chain parent mapping для summary card.
# Контракт: в заголовке «Итоговая активность» показываем РОДИТЕЛЬСКИЙ
# нуклид цепочки (Th-232 / Ra-226), а лучшую line-carrier нуклид-дочку
# (Ac-228 / Bi-214 / ...) уносим в note. В природных цепочках при
# секулярном равновесии активность дочки = активности родителя, и
# источник всегда — родитель.
#
# Mapping:
#   Th-232 chain — daughters Ac-228, Tl-208, Pb-212, Bi-212 → parent Th-232
#   U-238/Ra-226 chain — Rn-222 короткоживущие daughters
#     Pb-214, Bi-214 → parent Ra-226
#     (Ra-226 сам — parent, если найден по 186 кэВ)
#   U-235 chain — обычно через 185.7 / 143.8 кэВ → parent U-235
#   Одиночные/техногенные (Cs-137, Co-60, K-40, I-131, Eu-152, …) —
#     карта возвращает None, заголовок остаётся как есть.
_CHAIN_DAUGHTER_TO_PARENT: Dict[str, str] = {
    # Th-232 chain
    "Ac-228": "Th-232",
    "Tl-208": "Th-232",
    "Pb-212": "Th-232",
    "Bi-212": "Th-232",
    # Ra-226 chain (U-238 chain, post-Rn-222 segment).
    # Ra-226 сам как parent для Bi-214/Pb-214 — это distinction между
    # «измеряем Ra-226 по 186 кэВ» (parent сам) и «измеряем Ra-226 по
    # daughters» (через Bi-214). В обоих случаях reported parent = Ra-226.
    "Bi-214": "Ra-226",
    "Pb-214": "Ra-226",
}


def _chain_parent_for(nuclide: str) -> Optional[str]:
    """F-368: вернуть родительский нуклид цепочки, если ``nuclide`` —
    короткоживущая daughter природной цепочки в секулярном равновесии.
    Иначе None (одиночные / техногенные / сам parent)."""
    if not nuclide:
        return None
    return _CHAIN_DAUGHTER_TO_PARENT.get(nuclide)


def _build_summary_card(report: Dict[str, Any]) -> str:
    ids = report.get("identified_nuclides") or []
    # Pick the nuclide with highest CI / lowest sigma_rel
    if not ids:
        return (
            '    <div class="label">Удельная активность</div>\n'
            '    <div class="val">Удельная активность не определена</div>\n'
            '    <div class="note">См. диагностику.</div>'
        )
    # Sort: prefer not upper-limit, then by lowest activity_relative_sigma
    def score(n):
        rel = n.get("activity_relative_sigma")
        return (1 if n.get("is_upper_limit") else 0,
                rel if rel is not None else 99.0)
    ids_sorted = sorted(ids, key=score)
    best = ids_sorted[0]
    nuc = best.get("nuclide") or "?"
    sa = _safe_float(best.get("specific_activity_Bq_per_kg"))
    sa_s = _safe_float(best.get("specific_activity_sigma_Bq_per_kg"))
    if sa is not None and sa_s is not None:
        val = f"{sa:.0f} ± {sa_s:.0f} Бк/кг"
        is_specific = True
    elif sa is not None:
        val = f"{sa:.0f} Бк/кг"
        is_specific = True
    else:
        a = _safe_float(best.get("activity_Bq"))
        val = (f"{a:.0f} Бк" if a is not None else "верхний предел")
        # «верхний предел» без массы — не знаем kind, default «активность».
        is_specific = False
    n_lines = best.get("n_matched_lines") or 0

    # F-368 / v1.18.24.3 — chain-parent в заголовке для природных
    # цепочек (Th-232, Ra-226). Daughter best-line carrier уходит в note.
    # Лейбл активности зависит от единиц измерения: если значение в
    # Бк/кг (есть масса пробы) — это удельная активность; если в Бк
    # (точечный источник, нет массы) — общая активность.
    activity_kind = "удельная активность" if is_specific else "активность"
    chain_parent = _chain_parent_for(nuclide=nuc)
    # BUG-6 / v1.18.30+ (Agent B): расписана формула метода σ вместо
    # ссылки «по правилу F-NN.» — bare F-id срезается F-317 pipeline в
    # `build.py::_f317_apply_user_facing_compliance`, оставляя обрубленную
    # фразу «Метод σ по правилу .». Дескриптор пережёвывается strip-ом
    # без потери смысла.
    sigma_method_desc = (
        "σ итоговой активности — max(σ_взвешенное среднее, σ_разброс) "
        "по совпадающим линиям (консервативная агрегация ЛСРМ §7)."
    )
    if chain_parent:
        title = f"Итоговая {activity_kind} {chain_parent}"
        note = (
            f"Цепочка {chain_parent} в секулярном равновесии. "
            f"Измерено по дочке {nuc} ({n_lines} совпадающих линий). "
            f"{sigma_method_desc}"
        )
    else:
        title = f"Итоговая {activity_kind} {nuc}"
        note = (
            f"Лучшая нуклид-линия: {nuc} ({n_lines} совпадающих линий). "
            f"{sigma_method_desc}"
        )

    return (
        f'    <div class="label">{_esc(title)}</div>\n'
        f'    <div class="val">{_esc(val)}</div>\n'
        f'    <div class="note">{_esc(note)}</div>'
    )


def _build_notes_blocks(report: Dict[str, Any], analysis_result) -> str:
    """Notes — assembled from pipeline_notes / warnings / diagnostics."""
    diag = report.get("diagnostics") or {}
    cd = diag.get("chain_dominance") or {}
    multiplets_present = bool(report.get("multiplet_deconvolutions"))
    has_tl208_2614 = False
    for fp in (report.get("primary_feps") or []):
        if (fp.get("nuclide") == "Tl-208"
                and abs((_safe_float(fp.get("library_E_keV")) or 0) - 2614.5) < 5):
            has_tl208_2614 = True
            break
    th_dom = bool(cd.get("th232_dominant"))
    u_dom = bool(cd.get("u238_dominant"))
    suppressed = cd.get("suppressed_chains") or []

    blocks = []

    if th_dom or u_dom:
        chain_name = "Th-232" if th_dom else "U-238 / Ra-226"
        evidence = cd.get("th232_evidence") if th_dom else cd.get("u238_evidence")
        ev_text = ", ".join(evidence or []) if evidence else "ключевые линии цепочки"
        blocks.append(
            f'    <h3>Что определяет цепочку {chain_name}</h3>\n'
            f'    <p>Подтверждение цепочки {chain_name} строится на присутствии '
            f'{_esc(ev_text)}. В условиях секулярного равновесия активности '
            f'дочерних совпадают с материнским нуклидом.</p>'
        )

    if multiplets_present:
        blocks.append(
            '    <h3>Связанная подгонка мультиплетов</h3>\n'
            '    <p>Линии одного нуклида в мультиплете связываются через <b>одну '
            'свободную активность</b>; их площади распределяются по библиотечным '
            'интенсивностям. Решается как линейная задача наименьших квадратов с '
            'неотрицательными переменными (NNLS / lsq_linear).</p>'
        )

    # F-388 / v1.18.26 — conditional render: блок «Почему … не определяется»
    # имеет смысл только когда у нуклидов соответствующей цепочки есть
    # фактические присутствия в спектре (primary FEPs). Если ни одного
    # peak от U_chain / Th_chain не обнаружено, блок просто шум для
    # пользователя — render skip.
    if suppressed:
        feps = report.get("primary_feps") or []
        present_nuclides = {fp.get("nuclide") or "" for fp in feps}

        def _chain_has_evidence(chain_name: str) -> bool:
            name = chain_name.strip()
            if "U-238" in name or "Ra-226" in name or "U_" in name:
                members = _U_CHAIN_NUCLIDES
            elif "Th-232" in name or "Th_" in name:
                members = _TH_CHAIN_NUCLIDES
            else:
                # Неизвестная цепочка — оставляем прежнее поведение (render).
                return True
            return bool(present_nuclides & members)

        relevant = [c for c in suppressed if _chain_has_evidence(c)]
        if relevant:
            other = ", ".join(relevant)
            # BUG-6 / v1.18.30+ (Agent B): fallback-фраза переписана без
            # bare F-id — F-317 strip пропускает «F-89» и оставляет
            # «по правилу d.». Описательный fallback переживает strip.
            _suppress_fallback = (
                "Цепочка подавлена согласно привязке имени файла "
                "(filename binding): нет независимых линий, кроме хвостов "
                "других цепочек."
            )
            reason_raw = cd.get("suppression_reason") or _suppress_fallback
            reason_ru = _translate_note_line(reason_raw) or _suppress_fallback
            blocks.append(
                f'    <h3>Почему {_esc(other)} цепочка не определяется</h3>\n'
                f'    <p>{_esc(reason_ru)}</p>'
            )

    if has_tl208_2614:
        # BUG-6 / v1.18.30+ (Agent B): расписан смысл σε вместо ссылки
        # «по правилу F-NN» — bare F-id срезается F-317 pipeline в
        # `build.py::_f317_apply_user_facing_compliance`. Inline-определение
        # σε нужно потому что термин используется только здесь и не
        # раскрывается в других местах user-facing вывода.
        blocks.append(
            '    <h3>Линия 2614 кэВ — на границе кривой эффективности</h3>\n'
            '    <p>Tl-208 2614 кэВ сидит на верхнем краю Зоны 2 полинома EFA. '
            'На границе зоны σε (стандартное отклонение энергии пика, '
            'вычисленное из ковариационной матрицы полинома эффективности) '
            'может быть значительно больше табличного значения.</p>'
        )

    if th_dom:
        blocks.append(
            '    <h3>Зона 70-90 кэВ — без покомпонентного разложения</h3>\n'
            '    <p>Кластер 73-90 кэВ — композит Pb K-РИ от внутренней '
            'конверсии Pb-212 239 + Tl-208/Bi-212 + Th-228 84.4 + флуоресценция Pb. '
            'На NaI c ПШПВ ≈ 12 кэВ при шаге 2-3 кэВ покомпонентное разделение '
            'физически невозможно (F-110). Кластер не входит в расчёт активности.</p>'
        )

    # Warnings — pass through RU-translation filter (D-04, F-108).
    # BUG-40: warnings may be heterogeneous (str | dict). For dict
    # entries we render an F-386-compliant RU string up front, then
    # send the string through the existing _translate_note_line filter
    # so the rest of the pipeline is unchanged.
    warns = report.get("warnings") or []
    warns_strs = []
    for w in warns:
        if isinstance(w, dict):
            from gamma.reporting.markdown_report import _render_warning_dict_ru
            warns_strs.append(_render_warning_dict_ru(w))
        else:
            warns_strs.append(w)
    warns_ru = [_translate_note_line(w) for w in warns_strs]
    warns_ru = [w for w in warns_ru if w]
    if warns_ru:
        blocks.append(
            '    <h3>Предупреждения</h3>\n'
            '    <ul>\n' +
            "\n".join(f"      <li>{_esc(w)}</li>" for w in warns_ru) +
            '\n    </ul>'
        )

    # Pipeline notes — translate and drop English-leaked lines.
    pn = report.get("pipeline_notes") or []
    pn_ru = [_translate_note_line(n) for n in pn[:30]]
    pn_ru = [n for n in pn_ru if n]
    if pn_ru:
        blocks.append(
            '    <h3>Заметки пайплайна</h3>\n'
            '    <ul>\n' +
            "\n".join(f"      <li>{_esc(n)}</li>" for n in pn_ru[:20]) +
            '\n    </ul>'
        )

    # Conclusion — adapt narrative for pure background spectra (D-01).
    diag = report.get("diagnostics") or {}
    env = diag.get("measurement_environment", "")
    ids = report.get("identified_nuclides") or []
    detected_names = ", ".join(n.get("nuclide", "") for n in ids if n.get("nuclide"))
    if env == "background_only":
        if detected_names:
            conclusion = (
                f"<b>Фоновый спектр.</b> В фоне зарегистрированы линии: "
                f"{detected_names}. Никаких выводов про образец не делается — "
                f"это измерение фона."
            )
        else:
            conclusion = (
                "<b>Фоновый спектр.</b> Над порогом значимости линий не "
                "обнаружено. Никаких выводов про образец не делается — "
                "это измерение фона."
            )
    elif detected_names:
        conclusion = (
            f"<b>Идентифицированы нуклиды:</b> {detected_names}. "
            f"Подробности — в таблице характеристических линий выше."
        )
    else:
        conclusion = (
            "Характеристические линии выше порога значимости не найдены."
        )
    blocks.append(
        f'    <h3>Заключение</h3>\n'
        f'    <p>{conclusion}</p>'
    )
    return "\n\n".join(blocks)


def _build_cost_footer(cost_estimate: Optional[Dict[str, Any]]) -> str:
    """F-132 / v1.17.7 — ВСЕГДА выводит footer стоимости в HTML.

    Раньше (до F-132) при `cost_estimate=None` возвращалась пустая
    строка и footer вообще не рисовался. Теперь footer обязателен;
    при отсутствии данных от вызывающего показываем явный placeholder
    «оценка недоступна» — но НИКОГДА не молчим.
    """
    if not cost_estimate:
        cost_estimate = {
            "tokens": 0,
            "session_pct": "оценка стоимости недоступна",
            "detail": "Авто-оценка (F-132) не передана сборщику HTML. "
                      "Поэтапная разбивка должна приводиться в текстовом "
                      "отчёте (раздел «Оценка стоимости анализа»).",
        }
    tokens = cost_estimate.get("tokens") or 0
    session = cost_estimate.get("session_pct") or ""
    detail = cost_estimate.get("detail") or ""
    # F-132 формулировка (закреплена): «~N токенов или M% от бесплатной
    # 5-часовой сессии». Если CLI/build передал готовую session-строку
    # вида "8.5% от бесплатной 5-часовой сессии" — берём её целиком;
    # иначе показываем только токены без процентной части.
    tokens_str = f"~{int(tokens):,} токенов".replace(",", " ")
    if session:
        cost_line = f"{tokens_str} или {_esc(session)}"
    else:
        cost_line = tokens_str
    return (
        '  <div style="margin:36px 0 12px;padding:14px 16px;'
        'background:var(--bg-secondary);border-radius:var(--radius-md);'
        'font-size:12px;color:var(--text-secondary);line-height:1.55;">\n'
        '    <div style="font-size:11.5px;text-transform:uppercase;'
        'letter-spacing:0.04em;color:var(--text-tertiary);margin-bottom:4px;">'
        'Стоимость анализа (оценка) — F-132</div>\n'
        f'    <div style="color:var(--text-primary);font-weight:500;">'
        f'{cost_line}</div>\n'
        f'    <div style="margin-top:8px;color:var(--text-tertiary);'
        f'font-size:11px;font-style:italic;">{_esc(detail)}</div>\n'
        '  </div>'
    )


# ──────────────────────────────────────────────────────────────────
# F-FIT-VIEW / v1.22.1 — fit overlay payload for interactive HTML toggle
# ──────────────────────────────────────────────────────────────────

def _build_fit_overlay_payload(report: Dict[str, Any]) -> Dict[str, Any]:
    """F-FIT-VIEW / v1.22.1 — extract fit_overlay section from JSON report.

    The fit_overlay section is produced by json_report._build_fit_overlay().
    This function prepares a JS-friendly payload dict. Returns a safe
    empty structure when fit_overlay is absent (backward compat).

    Returns dict with:
      ``peaks``              — list of {peak_id, nuclide, energy_keV, amp_counts,
                               sigma_keV, source, label}
      ``multiplet_continua`` — list of {cluster_id, E_keV[], continuum[], total[],
                               components[{nuclide, energy_keV, sigma_keV,
                               amp_counts, label, g_curve[]}]}
    """
    fo = report.get("fit_overlay") or {}
    peaks = fo.get("peaks") or []
    multiplet_continua = fo.get("multiplet_continua") or []
    return {
        "peaks": peaks,
        "multiplet_continua": multiplet_continua,
    }


def _build_visual_similarity_card(report: Dict[str, Any], is_background_only: bool) -> str:
    """F-070 W3 / v1.24.0 — render F-VISUAL-SIMILARITY HTML card.

    Card placement: after F-FIT-VIEW, before Decision summary.

    Returns empty string for background-only spectra (consistent with bg-aware UI).
    Returns a collapsed skeleton card with reason message when enabled=False.
    Returns full top-3 match table otherwise.

    Cite: _state/agent_b/inbox/2026-06-04_F-070-W3_html_card_json_wiring.md §HTML card design
    """
    # Background-only: omit card entirely.
    if is_background_only:
        return ""

    vs = report.get("visual_similarity") or {}
    enabled = vs.get("enabled", False)

    # Verdict label localisation (RU) — F-108 compliance: no EN leak in body text.
    _VERDICT_RU = {
        "match":     "совпадение",
        "ambiguous": "неоднозначно",
        "mismatch":  "несовпадение",
    }

    if not enabled:
        reason = vs.get("reason") or "unknown"
        reason_esc = _esc(reason)
        return (
            '<div class="fp-vs-card" style="display:none;">\n'
            '  <div class="fp-vs-header">'
            '<span class="fp-vs-title">Визуальное сопоставление спектров</span>'
            '</div>\n'
            f'  <div class="fp-vs-disabled">Шаблоны недоступны (причина: {reason_esc})</div>\n'
            '</div>'
        )

    # Extract fields from the visual_similarity block.
    query_geometry = vs.get("query_geometry")
    matches = vs.get("matches") or []
    n_templates_total = 24  # W2 canonical count (F-070 brief)
    n_matches = len(matches)
    geometry_display = _esc(query_geometry or "все геометрии")

    subtitle = (
        f"Топ-{n_matches} совпадений из {n_templates_total} канонических шаблонов "
        f"(геометрия: {geometry_display})"
    )
    # Link text stays short (RU) so the file path in href doesn't leak as body text.
    policy_link = "audit/_rag/visual_templates/SIMILARITY_POLICY.md"
    policy_link_text = "RAG-047 → Политика сходства"

    # Table rows.
    rows_html = ""
    for rank, m in enumerate(matches, start=1):
        nuclide = _esc(m.get("nuclide") or "?")
        template_id = _esc(m.get("template_id") or "")
        cosine_raw = m.get("cosine_raw")
        cosine_adj = m.get("cosine_adjusted")
        verdict = (m.get("verdict") or "mismatch").lower()
        tier = (m.get("tier") or "").upper()
        stale = m.get("stale_reference", False)
        cert_dates = m.get("cert_reference_dates") or []

        # Cosine cell: show adjusted; tooltip if downweighted (F-070 W3 brief §Body).
        if (cosine_raw is not None and cosine_adj is not None
                and abs(cosine_raw - cosine_adj) > 0.001):
            # Tooltip text uses RU-friendly abbreviations, no forbidden EN words.
            tooltip_text = f"исх. {cosine_raw:.4f}, ×0.7 → {cosine_adj:.4f}"
            cosine_cell = (
                f'<span title="{_esc(tooltip_text)}" '
                f'style="cursor:help;text-decoration:underline dotted;">'
                f'{cosine_adj:.4f}</span>'
            )
        else:
            val = cosine_adj if cosine_adj is not None else cosine_raw
            cosine_cell = f"{val:.4f}" if val is not None else "—"

        # Verdict badge — RU label text, CSS class for colour coding.
        verdict_class = {
            "match":     "verdict-match",
            "ambiguous": "verdict-ambiguous",
            "mismatch":  "verdict-mismatch",
        }.get(verdict, "verdict-mismatch")
        verdict_label_ru = _VERDICT_RU.get(verdict, verdict)
        verdict_badge = f'<span class="fp-vs-badge {verdict_class}">{_esc(verdict_label_ru)}</span>'

        # Tier badge — just the letter (A/B/C); stale in Russian.
        tier_class = {"A": "tier-a-badge", "B": "tier-b-badge", "C": "tier-c-badge"}.get(tier, "tier-c-badge")
        tier_badge = f'<span class="fp-vs-badge {tier_class}">{_esc(tier)}</span>' if tier else ""
        if stale:
            tier_badge += ' <span class="fp-vs-badge stale-badge">устаревший</span>'

        # Certificate dates cell.
        cert_cell = ""
        if cert_dates:
            if len(cert_dates) == 1:
                cert_cell = f"Серт.: {_esc(cert_dates[0])}"
            else:
                years = sorted(set(d[:4] for d in cert_dates if len(d) >= 4))
                if len(years) >= 2:
                    cert_cell = f"Серт.: {_esc(years[0])}–{_esc(years[-1])}"
                else:
                    cert_cell = f"Серт.: {_esc(cert_dates[0])}"

        rows_html += (
            f"  <tr>\n"
            f"    <td>{rank}</td>\n"
            f"    <td>{nuclide}</td>\n"
            f"    <td><code style='font-family:var(--font-mono);font-size:10.5px;'>{template_id}</code></td>\n"
            f"    <td class='fp-num'>{cosine_cell}</td>\n"
            f"    <td>{verdict_badge}</td>\n"
            f"    <td>{tier_badge}</td>\n"
            f"    <td style='font-size:11px;color:var(--text-tertiary);'>{cert_cell}</td>\n"
            f"  </tr>\n"
        )

    # Verdict summary banner — fully in Russian.
    verdict_summary = (vs.get("verdict_summary") or "mismatch").lower()
    verdict_nuclide = vs.get("verdict_summary_nuclide")
    banner_class = {
        "match":     "match-banner",
        "ambiguous": "ambiguous-banner",
        "mismatch":  "mismatch-banner",
    }.get(verdict_summary, "mismatch-banner")
    verdict_summary_ru = _VERDICT_RU.get(verdict_summary, verdict_summary)
    if verdict_nuclide:
        banner_text = f"Лучшее: {_esc(verdict_nuclide)} ({_esc(verdict_summary_ru)})"
    elif verdict_summary == "ambiguous":
        banner_text = "Лучшее: неоднозначно → см. итоговое решение (F-157 ratio check)"
    else:
        banner_text = f"Лучшее: {_esc(verdict_summary_ru)}"

    # Disclaimer — no EN leak: decay_age описывается по-русски, TIER → «уровень C».
    disclaimer = (
        "Шаблоны NaI 63×63 7%-FWHM, возраст ≥ 15 лет помечены уровнем C ×0.7. "
        "Не заменяет F-157/F-330 — комбинируется в финальном итоговом решении."
    )

    return (
        '<div class="fp-vs-card">\n'
        '  <div class="fp-vs-header">\n'
        '    <span class="fp-vs-title">Визуальное сопоставление спектров</span>\n'
        f'    <a href="{_esc(policy_link)}" style="font-size:11px;color:var(--text-tertiary);" target="_blank">{policy_link_text}</a>\n'
        '  </div>\n'
        f'  <div class="fp-vs-subtitle">{_esc(subtitle)}</div>\n'
        '  <div style="overflow-x:auto;">\n'
        '  <table class="fp-vs-tbl">\n'
        '    <thead><tr>\n'
        '      <th>#</th><th>Нуклид</th><th>Шаблон</th>'
        '<th>Сходство</th><th>Вердикт</th><th>Ур.</th><th>Серт. дата</th>\n'
        '    </tr></thead>\n'
        f'    <tbody>\n{rows_html}    </tbody>\n'
        '  </table>\n'
        '  </div>\n'
        f'  <div class="fp-vs-verdict-banner {banner_class}">{banner_text}</div>\n'
        f'  <div class="fp-vs-disclaimer">{_esc(disclaimer)}</div>\n'
        '</div>'
    )


# ──────────────────────────────────────────────────────────────────
# Public entry point
# ──────────────────────────────────────────────────────────────────

def render_interactive_html(
    report: Dict[str, Any],
    analysis_result,
    *,
    cost_estimate: Optional[Dict[str, Any]] = None,
    bundle_index: bool = False,
) -> str:
    """Render the canonical v1.17.2 interactive HTML report.

    Parameters
    ----------
    bundle_index : bool, default False
        Mirrors ``cfg.reports.bundle_index`` (gate из #118). When True —
        emits ``<nav class="fp-back-nav">`` link to ``../index.html``;
        when False (default post-#118) — placeholder substitutes to empty
        string so the sub-report stands alone without a broken back-link.
    """
    template = _load_template()

    header = report.get("header") or {}
    diag = report.get("diagnostics") or {}

    # Sample mass
    sample_mass_kg = None
    for n in (report.get("identified_nuclides") or []):
        sa = n.get("specific_activity_Bq_per_kg")
        ab = n.get("activity_Bq")
        if sa is not None and ab and ab > 0:
            sample_mass_kg = ab / (sa if sa > 0 else 1.0)
            break
    # If wrapper passed it via spec.extras (best effort)
    if sample_mass_kg is None:
        try:
            sm = (analysis_result.spec.extras or {}).get("sample_mass_kg")
            if sm:
                sample_mass_kg = float(sm)
        except Exception:
            pass

    # Build text substitutions
    subs: Dict[str, str] = {}
    subs["TITLE"] = _build_title(header)
    subs["SUBTITLE"] = _build_subtitle(
        header, diag,
        f145=report.get("multiplet_self_calibration"),
        energy_cal=(report.get("calibration") or {}).get("energy_cal"),
    )
    subs["GRID_CARDS"] = _build_grid_cards(header, sample_mass_kg)
    subs["LEGEND_ITEMS"] = _build_legend_items()
    subs["CONTROLS_INFO"] = _build_controls_info(header)
    subs["TABLE_HEADERS"] = _build_table_headers()
    # F-UX-03 / 2026-06-04 — back-nav «← На главную» убрана по UX-correction
    # пользователя (см. _state/agent_a/inbox/2026-06-04_correction_9_*).
    # Supersedes F-RPT-02 (back-link к bundle index). Placeholder сохранён
    # пустым для совместимости с template ключом BACK_NAV.
    subs["BACK_NAV"] = ""

    # Data arrays
    E, C = _build_E_C(analysis_result)
    peaks = _build_peaks(report)
    detail = _build_detail(report, peaks)
    rows = _build_rows(report, analysis_result)
    multiplets = _build_multiplets_data(report, analysis_result)

    # F-111b: enforce set(peaks ids) == set(rows[].peak) == set(detail keys)
    _sync_peaks_rows_detail(peaks, rows, detail)

    # F-334.4 / v1.18.18.7 — пометить top-5 пиков по net_area для full
    # label; остальные получат tooltip-only маркеры (Chart.js annotation
    # с `display:false` на label, но видимой line+dot). Снижает clutter
    # на NORM-спектрах с 15+ линиями.
    _mark_top_peaks(peaks, report, top_n=5)

    # BUG-14 / v1.18.30+ (Agent B) — для чисто-фоновых спектров
    # (env == "background_only") разложение мультиплетов не имеет смысла:
    # пики не значимы над континуумом, χ²/ν подгонка к шумовому фону
    # вводит оператора в заблуждение. Полностью скрываем секцию (sample
    # и bg-suffix варианты) — оставляем data_* массивы пустыми, чтобы
    # canvas-таргетинг в JS не пытался рисовать пустые мультиплеты.
    # Источник данных (`multiplet_deconvolutions` /
    # `background_multiplet_deconvolutions`) не трогаем — деконволюция
    # остаётся в JSON для аудита, просто не рендерится в HTML.
    is_background_only = (diag.get("measurement_environment") == "background_only")

    if is_background_only:
        subs["MULTIPLET_BLOCKS"] = ""
        subs["MULTIPLET_BLOCKS_BG"] = ""
        rows_bg: List[Dict[str, Any]] = []
        multiplets_bg: List[Dict[str, Any]] = []
        # Sample-серия мультиплетов тоже не нужна — обнуляем DATA_MULTIPLETS
        # чтобы JS не пытался рисовать на отсутствующих canvas (блоки выше
        # пустые, canvas-элементы не создаются).
        multiplets = []
    else:
        # BUG-25 / v1.18.31+ (Agent B) — отсечь вырожденные мультиплеты
        # (все компоненты S=0 или катастрофический χ²/ν) ДО рендеринга.
        # Reporting layer; coupled_multiplet.py / deconvolve.py не
        # трогаются — фильтрация только на представлении. JSON
        # multiplet_deconvolutions сохраняет полный список для аудита.
        multiplets = _filter_meaningful_multiplets(multiplets)
        subs["MULTIPLET_BLOCKS"] = _build_multiplet_blocks(multiplets)
        # F-397 / v1.18.27 — bg-only peak block (rows + multiplet blocks). Когда
        # фон не анализировался (`background_*` keys пустые) — rows_bg=[] и
        # multiplet_blocks_bg = placeholder.
        bg_report_view = {
            "primary_feps": report.get("background_primary_feps") or [],
            "secondary_peaks": report.get("background_secondary_peaks") or [],
            "multiplet_deconvolutions":
                report.get("background_multiplet_deconvolutions") or [],
            "identified_nuclides": [],   # фон без активностей (compute_activities=False)
            "diagnostics": {},           # без chain_dominance — bg-таблица не получает chain placeholders
        }
        rows_bg = _build_rows(bg_report_view, analysis_result, is_bg=True)
        # BUG-13 / v1.18.30+ (Agent B) — sync bg-table с plot annotations.
        # Plot подписывает K-40 1461 / Cs-137 661 / Pb K-XRF из
        # _detect_bg_peaks (matching against _BG_LINES_DICT), но в
        # bg_report_view.primary_feps этих линий НЕТ (они отсекаются
        # chain-filtering в identification). Augmenter добавляет недостающие
        # строки используя тот же detect-источник что и plot.
        rows_bg = _augment_bg_rows_with_natural_lines(rows_bg, analysis_result)
        # BUG-43 / 2026-06-04 (Agent B) — НЕ выводить мультиплет-блок в BG-view.
        # Мультиплеты всегда относятся к образцу: разложение Eu-152/Pb-214
        # и аналогичные являются результатом анализа sample-спектра.
        # Показывать их под секцией «Найденные пики (фон)» вводит в заблуждение.
        # _build_multiplets_data читает analysis_result.deconvolution_results
        # (sample) вне зависимости от переданного report — данные bg-секции
        # совпадают с sample-данными, порождая ложный «мультиплет фона».
        # JSON background_multiplet_deconvolutions сохраняется для аудита.
        # Данные DATA_MULTIPLETS_BG обнуляются: canvas'ов нет → JS не рисует.
        multiplets_bg: List[Dict[str, Any]] = []
        subs["MULTIPLET_BLOCKS_BG"] = ""
    # F-UX-04 / 2026-06-04 — для чисто-фоновых спектров (env == "background_only")
    # блок «Итоговая удельная активность» не выводим: фоновый спектр — control,
    # активности нуклидов не измеряются и с паспортом не сравниваются. Источник:
    # _state/agent_a/inbox/2026-06-04_correction_10_no_activity_for_background.md
    # (создано пользователем на основе screenshot отчёта branch_b_background).
    if is_background_only:
        subs["SUMMARY_CARD"] = (
            '    <div class="label">Фоновый спектр</div>\n'
            '    <div class="val">Активность не измеряется</div>\n'
            '    <div class="note">Спектр используется как контроль фона. '
            'Удельная активность нуклидов и сравнение с паспортом '
            'неприменимы.</div>'
        )
    else:
        subs["SUMMARY_CARD"] = _build_summary_card(report)
    # F-070 W3 / v1.24.0 — visual similarity card (after F-FIT-VIEW, before Decision summary).
    subs["VISUAL_SIMILARITY_CARD"] = _build_visual_similarity_card(report, is_background_only)
    subs["NOTES_BLOCKS"] = _build_notes_blocks(report, analysis_result)
    subs["COST_FOOTER"] = _build_cost_footer(cost_estimate)

    # Serialize JS data — compact JSON
    subs["DATA_E"] = json.dumps(E, ensure_ascii=False)
    subs["DATA_C"] = json.dumps(C, ensure_ascii=False)
    subs["DATA_PEAKS"] = json.dumps(peaks, ensure_ascii=False)
    subs["DATA_DETAIL"] = json.dumps(detail, ensure_ascii=False)
    subs["DATA_ROWS"] = json.dumps(rows, ensure_ascii=False)
    subs["DATA_MULTIPLETS"] = json.dumps(multiplets, ensure_ascii=False)
    # F-397 / v1.18.27 — bg rows array (sample's view-toggle handler swaps in)
    subs["DATA_ROWS_BG"] = json.dumps(rows_bg, ensure_ascii=False)
    # F-397.2 / v1.18.28 (Agent B) — bg multiplets array. Параллельная
    # серия canvases (id="mp-<X>-bg") должна получать данные отдельно от
    # sample-серии (id="mp-<X>"), иначе bg multiplet charts остаются
    # пустыми при mode-bg. См. interactive_v1_17_2.html — multiplets_bg
    # forEach в конце инициализации Chart.js блока.
    subs["DATA_MULTIPLETS_BG"] = json.dumps(multiplets_bg, ensure_ascii=False)
    # F-332 / v1.18.18.5 — chart toggle payload (gross + bg + net на
    # общей energy axis с применённой live-time-scaling). Если фон не
    # вычитался — `has_background=False`, JS прячет toggle UI.
    chart_payload = _build_chart_payload(analysis_result)
    subs["DATA_CHART"] = json.dumps(chart_payload, ensure_ascii=False)
    # F-FIT-VIEW / v1.22.1 — fit overlay payload (per-peak Gaussians + multiplet
    # continua). Taken from report.fit_overlay (built by json_report._build_fit_overlay).
    fit_overlay_payload = _build_fit_overlay_payload(report)
    subs["DATA_FIT_OVERLAY"] = json.dumps(fit_overlay_payload, ensure_ascii=False)

    out = template
    for key, val in subs.items():
        out = out.replace("{{" + key + "}}", val)
    return out


__all__ = ["render_interactive_html"]
