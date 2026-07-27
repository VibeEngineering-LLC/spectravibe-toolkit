"""F2-A (2026-06-21) — homoglyph-resolve invariants for the canonical
detector rename ``Gamma-1C`` (ASCII «C», legacy typo) →
``Gamma-1S`` (ASCII «S», correct transliteration of cyrillic «Гамма-1С»).

Three operator-facing invariants the rename MUST preserve:

  1. All three forms — cyrillic «Гамма-1С», legacy ASCII-C «Gamma-1C»,
     canonical ASCII-S «Gamma-1S» — resolve to a single canonical
     string equal to ``"Gamma-1S"`` via ``canonicalize("detector", ...)``.

  2. The cyrillic header with the LSRM serial suffix
     («Гамма-1С №SN-01») also resolves to ``"Gamma-1S"`` (containment
     match on the canonicalizer, leaving the serial untouched).

  3. The resolved canonical ``"Gamma-1S"`` has a first-class profile on
     disk: ``detect_silent_fallback("Gamma-1S")`` returns
     ``reason == "profile_loaded_no_fallback"`` (no stub, no fallback).
     This kills the historical ``efficiency_tbd_using_fallback_profile``
     branch retired in F2-A Phase 5.

Why this test is needed (Censor brief):
  Phase 6 of the F2-A decree mandates a *dedicated* homoglyph-resolve
  test next to the broader BUG-40 warning-emission test. BUG-40
  (``tests/step04_detector_type/test_bug40_cyrillic_latin_warning.py``)
  proves the WARNING fires when the homoglyph DOES cause a fallback;
  this file proves the canonical rename ELIMINATES that fallback for
  the real «Гамма-1С» complex.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.data.aliases import canonicalize                              # noqa: E402
from gamma.detectors.profile import detect_silent_fallback               # noqa: E402


# ---------------------------------------------------------------------------
# Invariant 1+2 — three forms + serial-suffix containment match
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "raw",
    [
        "Гамма-1С",                  # cyrillic canonical input
        "Gamma-1C",                  # legacy ASCII-C (historical typo)
        "Gamma-1S",                  # ASCII-S canonical (idempotent)
        "Гамма-1С №SN-01",         # cyrillic + LSRM serial suffix
        "гамма-1с",                  # cyrillic lowercase
        "GAMMA-1S",                  # ASCII uppercase
        "Gamma_1S",                  # underscore separator
    ],
    ids=[
        "cyrillic-bare",
        "legacy-ascii-c",
        "ascii-s-idempotent",
        "cyrillic-with-serial",
        "cyrillic-lowercase",
        "ascii-uppercase",
        "underscore-sep",
    ],
)
def test_f2a_all_forms_resolve_to_canonical_gamma_1s(raw):
    """Every legitimate homoglyph / synonym → canonical 'Gamma-1S'."""
    canonical = canonicalize("detector", raw)
    assert canonical == "Gamma-1S", (
        f"F2-A invariant violated: canonicalize('detector', {raw!r}) = "
        f"{canonical!r}, expected 'Gamma-1S'"
    )


# ---------------------------------------------------------------------------
# Invariant 3 — Gamma-1S profile loads clean (no fallback)
# ---------------------------------------------------------------------------
def test_f2a_gamma_1s_profile_loads_without_fallback():
    """The canonical Gamma-1S detector must have a first-class profile.

    F2-A Phase 5 retired the bogus Gamma-1S stub-profile branch — if this
    test starts seeing a non-clean reason, the disk profile has been
    deleted or moved and the regression must be addressed before the
    fallback path silently re-emerges.
    """
    record = detect_silent_fallback("Gamma-1S")
    fallback = record.as_dict()
    assert fallback.get("reason") == "profile_loaded_no_fallback", (
        f"F2-A invariant violated: Gamma-1S profile load reported "
        f"reason={fallback.get('reason')!r} (expected "
        f"'profile_loaded_no_fallback'); full record = {fallback!r}"
    )
    # Both 'requested' and 'actual' must equal the canonical name.
    assert fallback.get("requested") == "Gamma-1S"
    assert fallback.get("actual") == "Gamma-1S"


# ---------------------------------------------------------------------------
# Negative — unrelated text should NOT accidentally resolve to Gamma-1S
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "raw",
    [
        "HPGe",                      # different detector type
        "NaI_63x63",                 # bare crystal, no Gamma-1 prefix
        "Xyzzy",                     # nonsense
    ],
)
def test_f2a_unrelated_does_not_resolve_to_gamma_1s(raw):
    canonical = canonicalize("detector", raw)
    assert canonical != "Gamma-1S", (
        f"F2-A invariant violated: canonicalize('detector', {raw!r}) "
        f"unexpectedly resolved to 'Gamma-1S'"
    )