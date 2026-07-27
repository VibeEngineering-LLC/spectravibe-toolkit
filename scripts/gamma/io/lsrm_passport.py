# -*- coding: utf-8 -*-
"""F-330 / v1.18.18.4 — LSRM .spe COMMENT passport parser.

LSRM SpectraLine calibration-source файлы хранят паспортные данные
эталона прямо в поле COMMENT, например:

    COMMENT=(Sum/T =    70.67 Sum*I/T =  8219.98)
            Cs-137 - 1890 Бк/кг (5%) от 30.05.97г.

или (К-40 пишется кириллической «К», уровень = только проценты без
референсной даты):

    COMMENT=(Sum/T =   21.56 Sum*I/T = 3770.29) К-40 - 2540 Бк/кг (10%)

Парсер извлекает структурированную информацию (LsrmPassportEntry) для
последующей передачи в `passport_activity_Bq` через wrapper
(F-326/F-330 auto-routing). Это закрывает разрыв «оператор должен
вручную ввести данные паспорта», когда сами данные уже физически
присутствуют в файле.

Поддерживаемые форматы:
- Nuclide: латиница (Cs, Co, Ra, Th, Na, K, ...) ИЛИ кириллица одно-двух
  буквенный символ + дефис + массовое число (с опциональным «m» для
  метастабильных, e.g. Tc-99m, Ba-137m)
- Активность: десятичная (с «.» или «,»), опционально с e-нотацией
- Единицы: «Бк/кг», «Бк·кг⁻¹», «Bq/kg», «Бк» (всеобщая активность)
- Погрешность: «(X%)» или «X%» — расширённая неопределённость U(P=0.95)
  если контекст не указан, считаем именно её (LSRM-конвенция)
- Дата: «от DD.MM.YY[YY]г?.?» или «from DD.MM.YY[YY]»
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from datetime import date, datetime
from typing import List, Optional


# ─── Cyrillic → Latin element-symbol map ────────────────────────────
# LSRM operators routinely type «К-40» (cyrillic К) instead of «K-40».
# Также встречается «С-137» (cyrillic С), «В-12» (cyrillic В = vit B12
# placeholder, не нуклид), «Н-3» (cyrillic Н = Latin H = tritium), и т.д.
# Маппинг по visually-identical Cyrillic→Latin glyphs.
_CYR_TO_LAT = {
    "А": "A", "В": "B", "Е": "E", "К": "K", "М": "M",
    "Н": "H", "О": "O", "Р": "P", "С": "C", "Т": "T",
    "Х": "X", "У": "Y",
    # lowercase fallbacks (rarely used in nuclide names but defensive)
    "а": "a", "в": "b", "е": "e", "к": "k", "м": "m",
    "н": "h", "о": "o", "р": "p", "с": "c", "т": "t",
    "х": "x", "у": "y",
}


def _normalize_nuclide_symbol(raw: str) -> str:
    """Convert Cyrillic look-alike chars to Latin in a nuclide symbol.

    «К-40» → «K-40», «С-137» → «C-137», «Н-3» → «H-3».
    """
    return "".join(_CYR_TO_LAT.get(ch, ch) for ch in raw)


# ─── Passport regex ─────────────────────────────────────────────────
# Components (named groups for clarity):
#   nuclide  — element-symbol + mass-number, e.g. «Cs-137», «К-40»,
#              «Tc-99m», «Ba-137m». Element symbol — 1 or 2 letters
#              (Cyrillic or Latin), Mass — 1-3 digits with optional «m».
#   value    — decimal activity, optional e-notation. Accepts both
#              «.» and «,» as decimal separator.
#   unit     — «Бк/кг» / «Бк·кг⁻¹» / «Bq/kg» / «Бк» / «Bq»
#   uncert   — uncertainty percent (with or without «%», parens
#              optional)
#   ref_date — optional reference date «от DD.MM.YY[YY]г?.?»
_PASSPORT_RE = re.compile(
    r"""
    (?P<nuclide>[A-ZА-ЯЁ][a-zа-яё]?-?\d{1,3}m?)   # nuclide symbol
    \s*[-–—]\s*                                    # separator dash
    (?P<value>\d+(?:[.,]\d+)?(?:[eE][+-]?\d+)?)    # numeric value
    \s*
    (?P<unit>Бк(?:/кг|·кг⁻¹)?|Bq(?:/kg)?|кБк(?:/кг)?|kBq(?:/kg)?)
    \s*
    \(?\s*
    (?P<uncert>\d+(?:[.,]\d+)?)\s*%?
    \s*\)?
    # Date prefix: «от/from» желателен, но допустимы малформированные
    # legacy LSRM cases (cp1251 mojibake — soft-hyphen + NBSP вместо
    # «от»). Принимаем любой короткий non-digit gap (0-8 chars) перед
    # date, ИЛИ канонический «от»/«from» с whitespace.
    \s*
    (?:
        (?:(?:от|from)\s+|[^\d\n]{0,8})
        (?P<day>\d{1,2})[.\-/]
        (?P<month>\d{1,2})[.\-/]
        (?P<year>\d{2,4})
        г?\.?
    )?
    """,
    re.VERBOSE | re.IGNORECASE,
)


# ─── BUG-49: Поверка-2016 lab COMMENT format ────────────────────────
# Operator lab «Поверка 2016» (Gamma-1S detector, 22 Tier-1 fixtures
# verified 2026-06-04) writes calibration-source passport data в COMMENT
# в собственном format'е, отличающемся от canonical F-330:
#
#   «<Nuclide> A=<value> <Бк|Бк/кг> dA=<unc>% DD-MM-YYYY»
#
# Differences vs canonical F-330:
#   - no dash between nuclide и value
#   - explicit «A=» prefix перед value
#   - «dA=» prefix перед uncertainty, with «%» suffix instead of «(%)»
#   - no «от» before date
#   - date format DD-MM-YYYY (dashes) vs canonical DD.MM.YY (dots)
#
# Examples (verbatim from probe outbox 2026-06-04_lsrm_spe_descriptions.json):
#   «Am-241 A=118000 Бк dA=5% 03-12-2013»       (source #42.13, Бк)
#   «Cd-109 A=1.033E6 Бк dA=2% 01-10-2008»      (source #SRC-05, e-notation)
#   «Cs-137 A=1760 Бк/кг dA=5% 24-05-2002»      (Marinelli 420-7-14, Бк/кг)
#
# Probe ran F-330 regex against all 22 files → 0/22 matches
# (silent corruption). New regex below matches 22/22.
_POVERKA_2016_RE = re.compile(
    r"""
    (?P<nuclide>[A-ZА-ЯЁ][a-zа-яё]?-?\d{1,3}m?)        # nuclide symbol
    \s+A\s*=\s*
    (?P<value>\d+(?:[.,]\d+)?(?:[eE][+-]?\d+)?)         # numeric value
    \s*
    (?P<unit>Бк/кг|Бк·кг⁻¹|Бк|Bq/kg|Bq|кБк/кг|кБк|kBq/kg|kBq)
    \s+
    dA\s*=\s*(?P<uncert>\d+(?:[.,]\d+)?)\s*%
    \s+
    (?P<day>\d{1,2})[-./](?P<month>\d{1,2})[-./](?P<year>\d{2,4})
    """,
    re.VERBOSE | re.IGNORECASE,
)


@dataclass
class LsrmPassportEntry:
    """One nuclide entry extracted from LSRM .spe COMMENT field.

    Attributes
    ----------
    nuclide : str
        Normalised nuclide symbol (Latin element, e.g. «Cs-137»,
        «K-40»). Cyrillic look-alikes converted to Latin.
    value : float
        Numeric activity as written в COMMENT (no unit conversion).
    unit : str
        Unit string as parsed: «Бк/кг», «Бк», «kBq», etc.
    uncertainty_pct : float
        Reported uncertainty (typically U(P=0.95) per LSRM convention).
    reference_date : datetime.date | None
        Reference date for decay correction; None if not given в COMMENT.
    raw_match : str
        Substring matched (для отладки и аудита).
    """

    nuclide: str
    value: float
    unit: str
    uncertainty_pct: float
    reference_date: Optional[date]
    raw_match: str

    @property
    def is_specific_activity(self) -> bool:
        """True если единица — Бк/кг (массовая активность)."""
        return "/" in self.unit or "·" in self.unit

    @property
    def is_kilo(self) -> bool:
        """True если префикс kBq / кБк (множитель ×1000)."""
        return self.unit.lower().startswith(("kbq", "кбк"))

    def value_Bq(self) -> Optional[float]:
        """Return activity in Бк (kBq → Bq), если is_specific_activity=False."""
        if self.is_specific_activity:
            return None
        return float(self.value) * (1000.0 if self.is_kilo else 1.0)

    def value_Bq_per_kg(self) -> Optional[float]:
        """Return specific activity in Бк/кг (kBq/kg → Bq/kg),
        если is_specific_activity=True."""
        if not self.is_specific_activity:
            return None
        return float(self.value) * (1000.0 if self.is_kilo else 1.0)


def _parse_value(raw: str) -> float:
    """Convert «1890» or «1,89e3» or «1.05E3» → float."""
    return float(raw.replace(",", "."))


def _parse_ref_date(
    day: Optional[str], month: Optional[str], year: Optional[str]
) -> Optional[date]:
    """Convert (D, M, Y) strings → date. Year may be 2 or 4 digits.

    For 2-digit years: <30 → 20XX (e.g. 25 → 2025); ≥30 → 19XX
    (e.g. 97 → 1997). LSRM legacy fixtures cover 1990s–2010s, so the
    cutoff at 30 catches both eras correctly through ~2029.
    """
    if not (day and month and year):
        return None
    try:
        d = int(day)
        m = int(month)
        y = int(year)
        if y < 100:
            y = 2000 + y if y < 30 else 1900 + y
        return date(y, m, d)
    except (ValueError, TypeError):
        return None


def parse_lsrm_passport_comment(
    comment: str,
) -> List[LsrmPassportEntry]:
    """Extract passport entries from an LSRM .spe COMMENT field.

    Returns list of LsrmPassportEntry — может быть пустой (e.g. COMMENT
    содержит только Sum/T бухгалтерию без passport-данных) или
    multi-element (некоторые лабораторные эталоны — multinuclide).

    Examples
    --------
    >>> entries = parse_lsrm_passport_comment(
    ...     "(Sum/T = 70.67) Cs-137 - 1890 Бк/кг (5%) от 30.05.97г."
    ... )
    >>> entries[0].nuclide, entries[0].value_Bq_per_kg(), entries[0].reference_date
    ('Cs-137', 1890.0, datetime.date(1997, 5, 30))

    >>> entries = parse_lsrm_passport_comment(
    ...     "(Sum/T = 21.56) К-40 - 2540 Бк/кг (10%)"
    ... )
    >>> entries[0].nuclide, entries[0].uncertainty_pct
    ('K-40', 10.0)
    """
    if not comment:
        return []

    out: List[LsrmPassportEntry] = []
    matched_spans: List[tuple] = []
    for m in _PASSPORT_RE.finditer(comment):
        try:
            raw_nuc = m.group("nuclide")
            nuc = _normalize_nuclide_symbol(raw_nuc).upper()
            # Promote first-letter-only-upper format (К-40 → K-40): the
            # normaliser keeps case, so re-title for canonical form.
            # «Cs-137» stays «CS-137» after .upper(), restore «Cs-137».
            nuc = _canonicalize_nuclide(nuc)
            value = _parse_value(m.group("value"))
            unit = m.group("unit")
            uncert = _parse_value(m.group("uncert"))
            ref = _parse_ref_date(
                m.group("day"), m.group("month"), m.group("year")
            )
            out.append(
                LsrmPassportEntry(
                    nuclide=nuc,
                    value=value,
                    unit=unit,
                    uncertainty_pct=uncert,
                    reference_date=ref,
                    raw_match=m.group(0).strip(),
                )
            )
            matched_spans.append(m.span())
        except (ValueError, AttributeError):
            # Malformed entry — skip silently, keep parsing the rest.
            continue

    # BUG-49: Fallback to Поверка-2016 lab format
    # «<Nuc> A=<v> <Бк|Бк/кг> dA=<u>% DD-MM-YYYY». Only matches that do
    # NOT overlap with canonical F-330 matches are added — this keeps
    # mixed-format COMMENTs (canonical + Поверка-2016) working without
    # double-counting.
    for m in _POVERKA_2016_RE.finditer(comment):
        ms, me = m.span()
        # Skip if this match overlaps with an F-330 canonical match.
        if any(not (me <= s or ms >= e) for (s, e) in matched_spans):
            continue
        try:
            raw_nuc = m.group("nuclide")
            nuc = _canonicalize_nuclide(
                _normalize_nuclide_symbol(raw_nuc).upper()
            )
            value = _parse_value(m.group("value"))
            unit = m.group("unit")
            uncert = _parse_value(m.group("uncert"))
            ref = _parse_ref_date(
                m.group("day"), m.group("month"), m.group("year")
            )
            out.append(
                LsrmPassportEntry(
                    nuclide=nuc,
                    value=value,
                    unit=unit,
                    uncertainty_pct=uncert,
                    reference_date=ref,
                    raw_match=m.group(0).strip(),
                )
            )
        except (ValueError, AttributeError):
            continue

    return out


def _canonicalize_nuclide(symbol: str) -> str:
    """Convert «CS-137» → «Cs-137», «CO-60» → «Co-60», «K-40» → «K-40».

    Element symbol: first letter upper, rest lower; mass number
    unchanged; metastable «m» suffix lower.
    """
    # Match: optional letters + dash + digits + optional 'm'
    m = re.match(r"^([A-Z]+)-?(\d+)([Mm]?)$", symbol)
    if not m:
        return symbol
    el = m.group(1)
    mass = m.group(2)
    meta = m.group(3).lower()
    # Title-case the element symbol
    el = el[:1].upper() + el[1:].lower()
    return f"{el}-{mass}{meta}"


# ─── Decay correction helpers ───────────────────────────────────────

# Half-lives in seconds for nuclides commonly seen в калибровочных
# источниках LSRM. Источник: NNDC / ENSDF + DDEP. Только те, у которых
# t½ соизмерим со срокам хранения эталона (~years).
# Для долгоживущих (K-40, Th-232, U-238) decay correction негligible.
_HALF_LIFE_S = {
    "H-3":    3.8878e8,    # 12.32 y
    "Na-22":  8.2090e7,    # 2.6018 y
    "Mn-54":  2.6960e7,    # 312.05 d
    "Co-57":  2.3489e7,    # 271.74 d
    "Co-60":  1.6633e8,    # 5.2711 y
    "Sr-90":  9.0876e8,    # 28.79 y (parent of Y-90)
    "Y-90":   2.3045e5,    # 64.0 h (transient eq with Sr-90)
    "Cd-109": 4.0061e7,    # 463.5 d
    "Sb-125": 8.7066e7,    # 2.7586 y
    "Cs-134": 6.5114e7,    # 2.0648 y
    "Cs-137": 9.4894e8,    # 30.07 y
    "Ba-133": 3.3137e8,    # 10.51 y
    "Eu-152": 4.2718e8,    # 13.537 y
    "Eu-154": 2.7137e8,    # 8.601 y
    "Eu-155": 1.5023e8,    # 4.753 y
    "Pb-210": 7.0224e8,    # 22.26 y (Ra-226 daughter)
    "Po-210": 1.1957e7,    # 138.376 d
    "Am-241": 1.3633e10,   # 432.2 y
    # Long-lived NORM / cosmogenic — decay correction <0.01% per century:
    "K-40":   4.0337e16,   # 1.2778e9 y
    "Ra-226": 5.0492e10,   # 1600 y
    "Th-232": 4.4334e17,   # 1.405e10 y
    "U-235":  2.2210e16,   # 7.04e8 y
    "U-238":  1.4099e17,   # 4.468e9 y
}


def half_life_seconds(nuclide: str) -> Optional[float]:
    """Return half-life в секундах для известных нуклидов; иначе None."""
    return _HALF_LIFE_S.get(nuclide)


def decay_correct(
    A0: float,
    nuclide: str,
    ref_date: date,
    meas_date: date,
) -> Optional[float]:
    """Apply exponential decay A(t) = A0 · exp(-λ·Δt).

    Args:
        A0: активность на reference date (Bq или Bq/kg — единицы не
            изменяются).
        nuclide: канонический symbol (e.g. «Cs-137»).
        ref_date: дата паспорта.
        meas_date: дата измерения.

    Returns:
        Корректированная активность на meas_date, или None если
        t½ не известна.
    """
    T = half_life_seconds(nuclide)
    if T is None or T <= 0:
        return None
    dt_s = (meas_date - ref_date).total_seconds() if isinstance(meas_date, datetime) else (
        (meas_date - ref_date).days * 86400.0
    )
    if dt_s == 0:
        return float(A0)
    import math
    lam = math.log(2.0) / T
    return float(A0) * math.exp(-lam * dt_s)


__all__ = [
    "LsrmPassportEntry",
    "parse_lsrm_passport_comment",
    "half_life_seconds",
    "decay_correct",
]
