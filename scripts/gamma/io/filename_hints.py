"""
Filename token parsing.

The filename is metadata. It often contains the sample identifier, nuclide
hints (Cs137, Co60), geometry (Marinelli, 1L), date, and the operator's
shorthand. Per SKILL.md: these tokens are PRIORS for candidate-list
building in step 7, not substitutes for spectral evidence.

Russian terms are recognised alongside English (Фон = background, Дски,
лаба = lab, etc.) because real-world filenames are mixed-language.
"""

from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path

# Nuclide tokens — matched case-insensitively. Order doesn't matter.
NUCLIDE_TOKENS = [
    "Cs137", "Cs134", "Co60", "Co57", "Eu152", "Am241", "Ra226", "Ba133",
    "Na22", "K40", "Th232", "U235", "U238", "Mn54", "Zn65", "Be7", "Mo99",
    "Tc99m", "I131", "F18", "Na24",
    # F-89b / v1.15.2 — extended for Th-chain / U-chain members so that
    # filenames mentioning chain proxies directly are recognised
    "Tl208", "Pb212", "Ac228", "Bi212", "Bi214", "Pb214", "Pb210",
    # BUG-36 (Wave 3 C1, 2026-06-05) — Ti-44 / Sc-44 calibration set.
    # AmTiCsEu Marinelli fixture filenames carry "Ti44" or "Sc44" tokens;
    # without these entries `extract_isotope_hints` silently dropped them,
    # so even when the filename signalled Ti-44 the candidate list never
    # bound the chain. Tokens are case-insensitive via `_match_token`.
    "Ti44", "Sc44",
]


# F-89b / v1.15.2 — canonical isotope label mapping.
# Filenames use a compact form like "Cs137"; the rest of the codebase
# (anchor_ranks, identify, nuclide library) uses canonical "Cs-137"
# with a hyphen between element and mass number. The mapping below is
# the single source of truth for that transform.
_NUCLIDE_CANONICAL: dict = {
    "Cs137": "Cs-137", "Cs134": "Cs-134",
    "Co60":  "Co-60",  "Co57":  "Co-57",
    "Eu152": "Eu-152",
    "Am241": "Am-241",
    "Ra226": "Ra-226",
    "Ba133": "Ba-133",
    "Na22":  "Na-22",  "Na24": "Na-24",
    "K40":   "K-40",
    "Th232": "Th-232",
    "U235":  "U-235",  "U238": "U-238",
    "Mn54":  "Mn-54",
    "Zn65":  "Zn-65",
    "Be7":   "Be-7",
    "Mo99":  "Mo-99",
    "Tc99m": "Tc-99m",
    "I131":  "I-131",
    "F18":   "F-18",
    # Th-chain / U-chain proxies in the extended NUCLIDE_TOKENS set
    "Tl208": "Tl-208", "Pb212": "Pb-212", "Ac228": "Ac-228",
    "Bi212": "Bi-212",
    "Bi214": "Bi-214", "Pb214": "Pb-214", "Pb210": "Pb-210",
    # BUG-36 (Wave 3 C1, 2026-06-05) — Ti-44 / Sc-44 calibration set.
    "Ti44": "Ti-44", "Sc44": "Sc-44",
}


# F-89b / v1.15.2 — chain membership table.
# Used by F-89d to decide whether a Th-only filename should suppress
# the U-238 chain (and vice versa). When a filename hint canonicalises
# to a nuclide listed here, the corresponding chain is "claimed" by
# the filename and competing chains are suppressed (unless they have
# strong independent evidence).
_NUCLIDE_TO_CHAIN: dict = {
    # Th-232 chain (Ra-228 → Ac-228 → Th-228 → Ra-224 → Rn-220 → Po-216 →
    # Pb-212 → Bi-212 → Po-212 / Tl-208 → Pb-208)
    "Th-232": "Th-232",
    "Tl-208": "Th-232", "Pb-212": "Th-232", "Ac-228": "Th-232",
    "Bi-212": "Th-232",
    # U-238 chain via Ra-226 → Rn-222 → Po-218 → Pb-214 → Bi-214 → Po-214
    "U-238": "U-238", "Ra-226": "U-238",
    "Bi-214": "U-238", "Pb-214": "U-238", "Pb-210": "U-238",
    "U-235": "U-235",
    # BUG-36 (Wave 3 C1, 2026-06-05) — Ti-44 / Sc-44 belong to the
    # Ti-44 → Sc-44 → Ca-44 (stable) decay chain. Both nuclides bind
    # the same chain so a "Ti44" filename also claims Sc-44 and vice
    # versa (secular equilibrium in any aged source).
    "Ti-44": "Ti-44", "Sc-44": "Ti-44",
}

# Geometry tokens (EN + RU — F-52a v1.11.1, added Russian geometry vocabulary
# used by LSRM SpectraLine on Gamma-1S, СЕГ-1КП, etc.).
# Ordered so that container-type tokens (Маринелли/Дента/Петри/Точ) match
# BEFORE distance-only tokens (0cm/5cm/25cm), because the container is the
# primary geometry identifier — distance is secondary.
GEOMETRY_TOKENS = [
    # Container types — Russian (LSRM-native), case-insensitive
    "Маринелли 1л", "Маринелли 0.5л",
    "Дента-120мл", "Петри-60мл",
    "Точ-25см", "Точ-10см", "Точ-5см",
    "Маринелли", "маринелли",
    "Дента", "дента",
    "Петри", "петри",
    "Точечный", "точечный", "Точ",
    "Колодец", "колодец",
    # Container types — English
    "Marinelli", "Petri", "Denta", "point",
    "1L", "500mL", "100mL", "GS5050", "GS3030", "GS7070",
    # Distance markers (fallback when no container token)
    "25cm", "20cm", "10cm", "5cm", "0cm",
]

# Detector tokens
DETECTOR_TOKENS = [
    "HPGe", "NaI", "LaBr", "LaBr3", "CeBr", "CeBr3", "CZT", "CdZnTe",
    "AtomSpectra",
]

# Sample-type and background tokens (incl. Russian)
SAMPLE_TYPE_TOKENS_EN = [
    "soil", "water", "air", "filter", "food", "concrete",
    "metal", "calib", "check_source", "bkg", "background",
]
SAMPLE_TYPE_TOKENS_RU = [
    # Background / measurement-type markers
    "фон", "Фон", "ФОН",
    "проба", "Проба", "ПРОБА",
    "лаба", "Лаба", "ЛАБА",
    "образец", "Образец",
    "Дски", "дски",
    "калиб", "Калиб",
    # F-52a v1.11.1: matrix-type markers commonly used by LSRM laboratories
    "Грунт", "грунт", "ГРУНТ",
    "Почва", "почва",
    "Вода", "вода", "ВОДА",
    "Воздух", "воздух",
    "Молоко", "молоко",
    "Зола", "зола",
    "Бетон", "бетон",
    "Пища", "пища",
    "Рыба", "рыба",
    "Хлеб", "хлеб",
    "Лес", "лес", "Лесная подстилка",
    "Песок", "песок",
    "Глина", "глина",
    "Торф", "торф",
    "Фильтр", "фильтр",
]


# A separator-class boundary that, unlike \b, treats underscores and
# CJK/Cyrillic non-word transitions as separators.
_SEP = r"(?:(?<=^)|(?<=[\s\-_\.\(\)\[\]])|(?<=[\u0400-\u04FF]))"
_END = r"(?:(?=$)|(?=[\s\-_\.\(\)\[\]])|(?=[\u0400-\u04FF]))"


def _match_token(token: str, text: str) -> bool:
    """True if token appears in text bounded by recognized separators."""
    pat = _SEP + re.escape(token) + _END
    return re.search(pat, text, re.IGNORECASE) is not None


def extract_isotope_hints(filename: str) -> list:
    """F-89b / v1.15.2 — extract canonical isotope labels from a filename.

    Returns a list of canonical labels (e.g. ``["Cs-137", "Co-60"]``).
    Empty list when no nuclide token is recognised. Order matches the
    raw NUCLIDE_TOKENS scan; duplicates collapsed.

    Use case: filenames like ``Cs137_420-7-14_Маринелли_0cm.spe`` carry
    a binding hypothesis — the user labelled the spectrum as Cs-137,
    so Cs-137 MUST enter the Step-7 candidate list per SKILL.md §7A.1
    ("Nuclides suggested by filename/metadata — highest priority").
    """
    name = Path(filename).stem
    out: list = []
    for tok in NUCLIDE_TOKENS:
        if _match_token(tok, name):
            canon = _NUCLIDE_CANONICAL.get(tok)
            if canon and canon not in out:
                out.append(canon)
    return out


def chains_claimed_by_isotope_hints(isotope_hints: list) -> set:
    """F-89b / v1.15.2 — which chains the filename binds.

    Returns a set like ``{"Th-232"}`` for ``["Th-232"]`` or
    ``{"Th-232", "Cs-137"}`` for a mixed Th+Cs source. Single-line
    nuclides without a chain (Cs-137, K-40, Co-60, Am-241) are
    returned as their own label — callers use the set as a chain
    restriction filter for F-89d (U-238 suppression rule).
    """
    out: set = set()
    for nuc in isotope_hints or ():
        chain = _NUCLIDE_TO_CHAIN.get(nuc, nuc)
        out.add(chain)
    return out


def parse_filename(filename: str) -> dict:
    """
    Extract recognized tokens from a filename.

    Returns a dict with these keys (all optional, only set if found):
      - nuclides: list[str]            — token-matched candidates (raw tokens)
      - isotope_hints: list[str]       — F-89b canonical labels (Cs-137, Th-232…)
      - chains_claimed: list[str]      — F-89b chains bound by the filename
      - geometry: str                  — first matched geometry token (raw)
      - geometry_canonical: str        — canonical name via gamma.data.aliases (F-78)
      - detector: str                  — first matched detector token (raw)
      - detector_canonical: str        — canonical (F-78)
      - sample_type: str               — first matched sample-type token (raw)
      - sample_type_canonical: str     — canonical (F-78)
      - is_background_hint: bool       — true if a background-marker token matched
      - date: datetime | None
      - raw: str                       — original stem
    """
    name = Path(filename).stem
    out = {
        "nuclides": [],
        "isotope_hints": [],
        "chains_claimed": [],
        "geometry": "",
        "geometry_canonical": "",
        "detector": "",
        "detector_canonical": "",
        "sample_type": "",
        "sample_type_canonical": "",
        "is_background_hint": False,
        "date": None,
        "raw": name,
    }

    # Nuclides — raw tokens + F-89b canonical labels
    for tok in NUCLIDE_TOKENS:
        if _match_token(tok, name):
            out["nuclides"].append(tok)
            canon = _NUCLIDE_CANONICAL.get(tok)
            if canon and canon not in out["isotope_hints"]:
                out["isotope_hints"].append(canon)
    out["chains_claimed"] = sorted(chains_claimed_by_isotope_hints(
        out["isotope_hints"]
    ))

    # Geometry
    for tok in GEOMETRY_TOKENS:
        if _match_token(tok, name):
            out["geometry"] = tok
            break

    # Detector
    for tok in DETECTOR_TOKENS:
        if _match_token(tok, name):
            out["detector"] = tok
            break

    # Sample type / background hints (EN + RU)
    for tok in SAMPLE_TYPE_TOKENS_EN:
        if _match_token(tok, name):
            out["sample_type"] = tok
            if tok.lower() in ("bkg", "background"):
                out["is_background_hint"] = True
            break
    if not out["sample_type"]:
        # Russian tokens — case-sensitive because Russian doesn't use \b reliably
        for tok in SAMPLE_TYPE_TOKENS_RU:
            if tok in name:
                out["sample_type"] = tok
                if tok.lower().startswith("фон"):
                    out["is_background_hint"] = True
                break

    # F-78: canonicalise via central registry. The per-token search above
    # has already picked the most-specific match for each category from
    # its ordered token list — we prefer that result over a whole-filename
    # containment scan (which can pick the first-by-length synonym from
    # any category, e.g. "Вода" → water in a "Фон Вода ..." filename
    # where "Фон" → background is the right answer for sample_type).
    try:
        from gamma.data.aliases import canonicalize
        for cat, raw_key, canon_key in (
            ("geometry", "geometry", "geometry_canonical"),
            ("detector", "detector", "detector_canonical"),
            ("sample_type", "sample_type", "sample_type_canonical"),
        ):
            c = ""
            # First — canonicalise the per-token raw match (authoritative)
            if out[raw_key]:
                c = canonicalize(cat, out[raw_key]) or ""
            # Fallback — try the whole filename (catches multi-word forms
            # like "Маринелли 1л" if per-token picked just "Маринелли")
            if not c:
                c = canonicalize(cat, name) or ""
            if c:
                out[canon_key] = c
    except Exception:
        # Alias registry is optional infrastructure; never break parse.
        pass

    # Date — try YYYYMMDD, YYYY-MM-DD, DD.MM.YYYY
    for pattern, fmt in (
        (r"(\d{4})(\d{2})(\d{2})",
         lambda m: datetime(int(m[1]), int(m[2]), int(m[3]))),
        (r"(\d{4})-(\d{2})-(\d{2})",
         lambda m: datetime(int(m[1]), int(m[2]), int(m[3]))),
        (r"(\d{2})\.(\d{2})\.(\d{4})",
         lambda m: datetime(int(m[3]), int(m[2]), int(m[1]))),
    ):
        m = re.search(pattern, name)
        if m:
            try:
                out["date"] = fmt(m)
                break
            except (ValueError, KeyError):
                pass

    return out
