"""BUG-40 (P0-NEW Wave 3 C4, v1.25.0) — Cyrillic→Latin homoglyph
detector-name substitution emits a structured ``warnings`` entry with
machine-readable ``code = "DETECTOR_CYRILLIC_LATIN_FALLBACK"``.

KFI: ``KNOWN_AND_FIXED_ISSUES.md:1401-1422``.

Symptom (verbatim from KFI):
  Spectrum ``CONFIGNAME = "Гамма-1С №SN-02"`` (Cyrillic «С») →
  pipeline reports ``detector_canonical = "Gamma-1S"`` (Latin «C»)
  without any visible warning in ``report.json``. Efficiency file
  loaded: ``УДС-ГЦ-63х63-USB_-_Маринелли.efr`` from Поверка-2024 on
  Gamma-1S №0086, NOT the Gamma-1S Поверка-2016 cert. This drives
  Cs-137 +9.4% residual.

Fix verified here:
  * ``scripts/gamma/data/aliases.py`` exposes
    ``contains_cyrillic_letters()`` + ``cyrillic_to_latin_collision()``.
  * ``scripts/gamma/reporting/json_report.py::_build_warnings`` emits a
    dict with ``code = "DETECTOR_CYRILLIC_LATIN_FALLBACK"`` when:
      (a) the detector_fallback record carries a non-clean reason; AND
      (b) the winning canonicalization source raw string contained
          Cyrillic; AND
      (c) the resolved canonical is pure ASCII.

Tests cover:
  A. Predicate ``contains_cyrillic_letters`` over canonical examples.
  B. Predicate ``cyrillic_to_latin_collision`` over canonical examples.
  C. ``_build_warnings`` emits the structured DETECTOR_CYRILLIC_LATIN_FALLBACK
     dict for the brief's Cyrillic-«С» fixture.
  D. ``_build_warnings`` does NOT emit the structured dict when the
     detector header is pure ASCII (backward compat).
  E. ``_build_warnings`` does NOT emit the structured dict when the
     profile loaded cleanly (no fallback).
"""
from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.data.aliases import (                                       # noqa: E402
    contains_cyrillic_letters,
    cyrillic_to_latin_collision,
)
from gamma.reporting.json_report import _build_warnings                # noqa: E402


# ---------------------------------------------------------------------------
# Test A — contains_cyrillic_letters predicate
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "raw, expected",
    [
        # Pure Cyrillic word.
        ("Гамма", True),
        # Mixed Cyrillic + Latin + punctuation + digits.
        ("Гамма-1С №SN-02", True),
        # Single Cyrillic letter inside otherwise ASCII text.
        ("Gamma-1С", True),
        # Pure Latin — no Cyrillic.
        ("Gamma-1S", False),
        # Punctuation + digits only — no letters at all.
        ("№SN-02", False),
        # Empty string.
        ("", False),
        # None defensive.
        (None, False),
        # Mathematical Greek letter — NOT Cyrillic.
        ("αβγ", False),
    ],
    ids=[
        "pure-cyrillic",
        "configname-mixed",
        "single-cyrillic-in-ascii",
        "pure-latin",
        "no-letters",
        "empty",
        "none",
        "greek-not-cyrillic",
    ],
)
def test_contains_cyrillic_letters(raw, expected):
    assert contains_cyrillic_letters(raw) is expected


# ---------------------------------------------------------------------------
# Test B — cyrillic_to_latin_collision predicate
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "raw, canonical, expected",
    [
        # KFI fixture: Cyrillic header → ASCII canonical = collision.
        ("Гамма-1С №SN-02", "Gamma-1S", True),
        # Bare Cyrillic header → ASCII canonical Gamma-1S.
        ("Гамма-1С", "Gamma-1S", True),
        # Legacy ASCII-C target — predicate stays True for any pure-ASCII
        # canonical (preserves backward-compat after F2-A canonical rename).
        ("Гамма-1С", "Gamma-1" + chr(0x43), True),
        # Latin header → ASCII canonical = NO collision.
        ("Gamma-1S", "Gamma-1S", False),
        # Cyrillic header with Cyrillic canonical = NO collision
        # (theoretical — current registry does not produce this).
        ("Гамма-1С", "Гамма-1С", False),
        # Empty raw.
        ("", "Gamma-1S", False),
        # None raw.
        (None, "Gamma-1S", False),
        # Empty canonical.
        ("Гамма-1С", "", False),
        # None canonical.
        ("Гамма-1С", None, False),
    ],
    ids=[
        "configname-S-with-serial-to-Gamma-1S",
        "bare-cyrillic-to-Gamma-1S",
        "bare-cyrillic-to-legacy-ASCII-C",
        "latin-S-to-Gamma-1S-no-collision",
        "cyrillic-to-cyrillic-no-collision",
        "empty-raw",
        "none-raw",
        "empty-canonical",
        "none-canonical",
    ],
)
def test_cyrillic_to_latin_collision(raw, canonical, expected):
    assert cyrillic_to_latin_collision(raw, canonical) is expected


# ---------------------------------------------------------------------------
# Test C — _build_warnings emits DETECTOR_CYRILLIC_LATIN_FALLBACK
# ---------------------------------------------------------------------------
def _make_minimal_result(detector_fallback):
    """Construct a minimal SimpleNamespace stand-in for StagedAnalysisResult."""
    return SimpleNamespace(
        efficiency_curve=None,
        activities=None,
        mda_per_line=None,
        next_stage_recommended=None,
        next_stage_reason="",
        detector_fallback=detector_fallback,
    )


def test_build_warnings_emits_cyrillic_latin_fallback_dict():
    """Cyrillic detector header + missing profile → structured dict (BUG-40).

    F2-A note (2026-06-21): the real «Гамма-1С» CONFIGNAME now resolves
    to canonical Gamma-1S whose profile loads cleanly, so this test
    models a *hypothetical* future cyrillic complex («Гамма-XYZ») that
    has no profile JSON yet — the cyrillic→ASCII collision predicate
    still fires AND the fallback is real (profile_not_on_disk), so the
    structured warning must be emitted.
    """
    fallback = {
        "requested": "Gamma-XYZ",
        "actual": "Gamma-1S",
        "reason": "profile_not_on_disk",
        "human": (
            "Detector profile references/detectors/Gamma-XYZ.json not "
            "found on disk — pipeline fell back to Gamma-1S."
        ),
        "human_en": (
            "Detector profile references/detectors/Gamma-XYZ.json not "
            "found on disk — pipeline fell back to Gamma-1S."
        ),
        "human_ru": (
            "Профиль детектора Gamma-XYZ отсутствует на диске — "
            "применены параметры Gamma-1S."
        ),
        # BUG-40 fields appended by staged_pipeline.py:
        "original_raw": "Гамма-XYZ №0001-26",
        "cyrillic_to_latin_collision": True,
    }
    result = _make_minimal_result(fallback)
    warnings = _build_warnings(result)

    # At least one entry must be a dict with code = DETECTOR_CYRILLIC_LATIN_FALLBACK.
    code_entries = [
        w for w in warnings
        if isinstance(w, dict) and w.get("code") == "DETECTOR_CYRILLIC_LATIN_FALLBACK"
    ]
    assert len(code_entries) == 1, (
        f"Expected exactly one DETECTOR_CYRILLIC_LATIN_FALLBACK dict, "
        f"got warnings = {warnings!r}"
    )
    entry = code_entries[0]
    # Verify all required fields from KFI:1401-1422 brief structure.
    assert entry["code"] == "DETECTOR_CYRILLIC_LATIN_FALLBACK"
    assert "message" in entry and entry["message"]
    assert "original_detector" in entry
    assert "mapped_to" in entry
    assert entry["severity"] == "MEDIUM"
    # original_detector should reflect the Cyrillic original_raw.
    assert "Гамма-XYZ" in entry["original_detector"]
    # mapped_to should be the ASCII canonical.
    assert all(ord(c) < 0x80 for c in entry["mapped_to"])


def test_build_warnings_no_cyrillic_warning_for_ascii_detector():
    """ASCII detector header → no DETECTOR_CYRILLIC_LATIN_FALLBACK dict (backward compat)."""
    fallback = {
        "requested": "Gamma-XYZ",
        "actual": "Gamma-1S",
        "reason": "profile_not_on_disk",
        "human": (
            "Detector profile references/detectors/Gamma-XYZ.json not "
            "found on disk — pipeline fell back to Gamma-1S."
        ),
        "human_en": (
            "Detector profile references/detectors/Gamma-XYZ.json not "
            "found on disk — pipeline fell back to Gamma-1S."
        ),
        "human_ru": (
            "Профиль детектора Gamma-XYZ отсутствует на диске — "
            "применены параметры Gamma-1S."
        ),
        "original_raw": "Gamma-XYZ",
        "cyrillic_to_latin_collision": False,
    }
    result = _make_minimal_result(fallback)
    warnings = _build_warnings(result)

    code_entries = [
        w for w in warnings
        if isinstance(w, dict) and w.get("code") == "DETECTOR_CYRILLIC_LATIN_FALLBACK"
    ]
    assert code_entries == [], (
        f"Did NOT expect a DETECTOR_CYRILLIC_LATIN_FALLBACK dict for an "
        f"ASCII detector header, got: {code_entries!r}"
    )
    # The text warning for the fallback itself should still be present.
    assert any(
        isinstance(w, str) and "Gamma-1S" in w
        for w in warnings
    )


def test_build_warnings_no_cyrillic_warning_when_profile_loaded_clean():
    """Profile loaded cleanly → no warning at all, even with Cyrillic raw."""
    fallback = {
        "requested": "Gamma-1S",
        "actual": "Gamma-1S",
        "reason": "profile_loaded_no_fallback",
        "human": "",
        "human_en": "",
        "human_ru": "",
        "original_raw": "Гамма-1С",  # Cyrillic raw but no fallback path → no warning.
        "cyrillic_to_latin_collision": True,
    }
    result = _make_minimal_result(fallback)
    warnings = _build_warnings(result)

    code_entries = [
        w for w in warnings
        if isinstance(w, dict) and w.get("code") == "DETECTOR_CYRILLIC_LATIN_FALLBACK"
    ]
    assert code_entries == [], (
        f"Expected NO DETECTOR_CYRILLIC_LATIN_FALLBACK dict when profile "
        f"loaded cleanly (no fallback path), got: {code_entries!r}"
    )


def test_build_warnings_no_detector_fallback_at_all():
    """No detector_fallback attribute → no detector warnings of any kind."""
    result = _make_minimal_result(None)
    warnings = _build_warnings(result)
    assert all(
        not (isinstance(w, dict) and w.get("code") == "DETECTOR_CYRILLIC_LATIN_FALLBACK")
        for w in warnings
    )
