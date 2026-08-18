"""
Canonical-name + synonym registry (F-78 / v1.12.0).

v1.12.0 (F2-A, 2026-06-21): canonical 'Gamma-1' + 'C' (Latin)  →  'Gamma-1' + 'S' (ASCII transliteration of cyrillic 'С')
— cyrillic «Гамма-1С» transliterates to ASCII 'Gamma-1S' (S = correct
transliteration of cyrillic С, U+0421 → ASCII S). The legacy ASCII-C
form (historical homoglyph typo) is kept as a synonym for backward
compatibility — all three forms (`Гамма-1С`, ASCII-C legacy, `Gamma-1S`)
now resolve to canonical `Gamma-1S`. The cyrillic_to_latin_collision()
predicate still returns True for cyrillic-raw → ASCII-canonical pairs;
post-merge warnings only fire when the profile-load path itself
fell back (gate in json_report._build_warnings), which no longer
happens for the canonical detector.

Translates raw text tokens — geometry, detector, software, sample matrix
— from any of their Russian/English/informal/with-or-without-diacritics
forms to a single stable canonical name (snake_case ASCII).

Design principles
-----------------
* The set of canonicals is small, stable, and ASCII; downstream code
  (algorithms, file paths, JSON keys, .efr lookup) MUST use canonicals.
* Synonyms are matched case-INsensitively with whitespace + punctuation
  normalisation, so `"Маринелли 1 л"` ≡ `"маринелли_1l"` ≡ `"Marinelli 1L"`.
* When multiple canonicals match a raw token, the FIRST in `aliases.json`
  wins — order in the file is the authoritative tiebreaker.
* The registry lives in `data/aliases.json` for ease of editing and
  user-approval workflow; this module just loads + indexes + queries.

Public API
----------
* :func:`canonicalize(category, raw)`  →  canonical name | None
* :func:`is_known(category, raw)`      →  bool
* :func:`synonyms_of(category, canonical)`  →  list[str]
* :func:`list_canonicals(category)`    →  list[str]
* :data:`CATEGORIES`                   →  ("geometry", "detector", "software", "sample_type")
"""

from __future__ import annotations

import json
import re
import unicodedata
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Tuple


CATEGORIES: Tuple[str, ...] = ("geometry", "detector", "software", "sample_type")


# Path to the registry JSON. Resolved once at import time; tests can monkey-patch.
_REGISTRY_PATH = (
    Path(__file__).resolve().parents[3] / "data" / "aliases.json"
)


# ──────────────────────────────────────────────────────────────────
# Normalisation
# ──────────────────────────────────────────────────────────────────

# Strip combining marks (NFD-decomposed Unicode). Also fold Cyrillic
# letters with stress accents to base letters. Leave the underlying
# Cyrillic codepoints unchanged.
def _strip_combining(s: str) -> str:
    return "".join(
        c for c in unicodedata.normalize("NFD", s)
        if unicodedata.category(c) != "Mn"
    )


# Map all whitespace and punctuation runs to a single ASCII space.
_PUNCT_RE = re.compile(r"[\s\-_./×x×]+")


def _normalise(s: str) -> str:
    """Canonical normalisation for synonym lookup.

    Strip diacritics → casefold → unify separators → collapse spaces.
    Note: Cyrillic letters are kept (Cyrillic is a primary script for
    LSRM filenames).
    """
    if s is None:
        return ""
    s = str(s)
    s = _strip_combining(s)
    s = s.casefold()
    s = _PUNCT_RE.sub(" ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


# ──────────────────────────────────────────────────────────────────
# Registry loading
# ──────────────────────────────────────────────────────────────────

@lru_cache(maxsize=1)
def _load_registry() -> Dict[str, Dict[str, List[str]]]:
    """Parse aliases.json and strip _meta. Cached for the process lifetime."""
    with _REGISTRY_PATH.open("r", encoding="utf-8") as f:
        raw = json.load(f)
    return {k: v for k, v in raw.items() if not k.startswith("_")}


@lru_cache(maxsize=1)
def _build_lookup_index() -> Dict[str, List[Tuple[str, str]]]:
    """
    Build: category → list of (normalised_synonym, canonical_name).

    The list preserves the source order of canonicals in aliases.json so
    that the first canonical wins on tie.
    """
    out: Dict[str, List[Tuple[str, str]]] = {}
    reg = _load_registry()
    for cat in CATEGORIES:
        cat_map = reg.get(cat, {})
        flat: List[Tuple[str, str]] = []
        for canonical, syns in cat_map.items():
            seen_local = set()
            # Always include the canonical itself as a synonym.
            for s in [canonical] + list(syns):
                n = _normalise(s)
                if n and n not in seen_local:
                    flat.append((n, canonical))
                    seen_local.add(n)
        out[cat] = flat
    return out


# ──────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────

def canonicalize(category: str, raw: str) -> Optional[str]:
    """
    Resolve any synonym (or the canonical itself) to its canonical name.

    Matches longest-first within the category. Returns None when no
    synonym matches. Case-insensitive, punctuation-insensitive.

    Examples:
        canonicalize("geometry", "Маринелли 1л") → "marinelli_1L"
        canonicalize("detector", "БДЭГ-63×63")   → "NaI_63x63"
        canonicalize("software", "ЛСРМ SpectraLine") → "lsrm_spectraline"
        canonicalize("sample_type", "Грунт")     → "soil"
        canonicalize("geometry", "Xyzzy")        → None
    """
    if category not in CATEGORIES:
        raise ValueError(f"unknown alias category: {category!r}")
    if not raw:
        return None
    norm = _normalise(raw)
    if not norm:
        return None
    idx = _build_lookup_index()[category]
    # Exact match first (most common).
    for (syn, canonical) in idx:
        if syn == norm:
            return canonical
    # Containment match: if the raw text *contains* a known synonym as a
    # whole token (so that "soil_sample_2024" → "soil"). Try longest
    # synonyms first to avoid matching "soil" inside "soilage".
    for (syn, canonical) in sorted(idx, key=lambda x: -len(x[0])):
        if not syn:
            continue
        if re.search(r"(?<![\wЀ-ӿ])" + re.escape(syn) +
                     r"(?![\wЀ-ӿ])", norm):
            return canonical
    return None


def is_known(category: str, raw: str) -> bool:
    """True if `raw` resolves to any canonical in `category`."""
    return canonicalize(category, raw) is not None


def synonyms_of(category: str, canonical: str) -> List[str]:
    """Return the human-readable synonym list (incl. the canonical itself)."""
    reg = _load_registry()
    cat = reg.get(category, {})
    if canonical not in cat:
        return []
    return [canonical] + list(cat[canonical])


def list_canonicals(category: str) -> List[str]:
    """List of canonical names in registry order."""
    return list(_load_registry().get(category, {}).keys())


# ──────────────────────────────────────────────────────────────────
# BUG-40 — Cyrillic / Latin homoglyph detection
# ──────────────────────────────────────────────────────────────────
#
# The LSRM SpectraLine ``CONFIGNAME`` header frequently encodes the
# detector complex with a mix of Cyrillic + Latin letters that look
# identical on screen (А↔A, В↔B, С↔C, Е↔E, Н↔H, К↔K, М↔M, О↔O, Р↔P,
# Т↔T, Х↔X, У↔Y, …). When the canonicalizer maps such a header to a
# pure-ASCII canonical (e.g. ``"Гамма-1С"`` → ``"Gamma-1S"``), the
# substitution is invisible to the operator unless the report flags it.
#
# BUG-40 (KFI:1401-1422) requires an explicit warning whenever:
#
#   1. the raw header contains at least one Cyrillic letter, AND
#   2. the resolved canonical is pure ASCII, AND
#   3. the canonical profile triggers a fallback (since F2-A 2026-06-21
#      only ``profile_not_on_disk`` remains — the legacy
#      ``efficiency_tbd_using_fallback_profile`` reason was retired
#      together with the bogus Gamma-1S stub profile) — i.e. the
#      Cyrillic complex is not first-class in the detector registry.
#
# This module exposes the two cheap pure-string predicates; the
# pipeline composes them with the existing :mod:`gamma.detectors.profile`
# fallback record to decide whether a structured warning fires.

# Cyrillic alphabet blocks:
#   * U+0400–U+04FF: Cyrillic (basic + supplement)
#   * U+0500–U+052F: Cyrillic Supplement
# Punctuation, digits, ASCII letters and whitespace are NOT covered
# here — only "is this codepoint a Cyrillic letter" matters.

def contains_cyrillic_letters(raw: Optional[str]) -> bool:
    """True iff ``raw`` contains at least one Cyrillic letter.

    Punctuation, digits, ASCII letters and whitespace do NOT count.
    Returns False for ``None`` / empty string.

    Examples:
        contains_cyrillic_letters("Гамма-1С")     → True
        contains_cyrillic_letters("Gamma-1S")     → False
        contains_cyrillic_letters("Gamma-1С")     → True  (mixed)
        contains_cyrillic_letters("№SN-02")     → False (№ is U+2116, not Cyrillic)
        contains_cyrillic_letters("")             → False
        contains_cyrillic_letters(None)           → False
    """
    if not raw:
        return False
    for ch in raw:
        cp = ord(ch)
        if 0x0400 <= cp <= 0x04FF or 0x0500 <= cp <= 0x052F:
            return True
    return False


def cyrillic_to_latin_collision(
    raw: Optional[str], canonical: Optional[str]
) -> bool:
    """True iff ``raw`` is Cyrillic-bearing AND ``canonical`` is pure ASCII.

    This is the BUG-40 trigger predicate. Use it AFTER a successful
    canonicalize() round-trip to flag the homoglyph substitution.

    Examples:
        cyrillic_to_latin_collision("Гамма-1С", "Gamma-1S")   → True
        cyrillic_to_latin_collision("Gamma-1S", "Gamma-1S")   → False
        cyrillic_to_latin_collision("Гамма-1С", "Гамма-1С")   → False (canon also CYR)
        cyrillic_to_latin_collision("", "Gamma-1S")           → False (empty raw)
        cyrillic_to_latin_collision("Гамма-1С", "")           → False (no canonical)
    """
    if not raw or not canonical:
        return False
    if not contains_cyrillic_letters(raw):
        return False
    # Canonical is "pure ASCII" iff every codepoint < 0x80.
    return all(ord(c) < 0x80 for c in canonical)


__all__ = [
    "CATEGORIES",
    "canonicalize",
    "is_known",
    "synonyms_of",
    "list_canonicals",
    "contains_cyrillic_letters",
    "cyrillic_to_latin_collision",
]
