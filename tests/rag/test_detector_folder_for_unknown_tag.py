# -*- coding: utf-8 -*-
"""V126-04 — Unit tests for ``_detector_folder_for`` unknown-tag handling.

Pre-fix gap (brief 2026-06-05_P2_V126-04-05):
    ``scripts/rag/build_spectra_index.py:_detector_folder_for`` silently
    returned an empty string for any detector tag not present in
    ``_DETECTOR_TAG_TO_FOLDER``. ``test_detector_folder_field_in_all_records``
    in ``tests/detectors/test_detector_folder_structure.py`` asserted only
    ``"detector_folder" not in rec`` (key presence) — empty strings passed.

Post-fix contract:
    * Unmapped non-empty tags raise a ``UserWarning`` whose message starts
      with ``"_detector_folder_for: unmapped detector tag"`` (Option A in
      the brief).
    * Known mapped tags still return their canonical folder.
    * Empty / whitespace-only inputs still return ``""`` silently (no warn).
    * A new SPECTRA_INDEX-level test asserts that every record has a
      *non-empty* ``detector_folder`` (RED if any future record ships an
      empty string for the field).
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SPECTRA_INDEX_PATH = REPO_ROOT / "audit" / "_rag" / "SPECTRA_INDEX.json"


# ─── Fixture ─────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def spectra_index() -> dict:
    if not SPECTRA_INDEX_PATH.exists():
        pytest.skip(f"{SPECTRA_INDEX_PATH.name} not present — F-300/W3 not run yet.")
    return json.loads(SPECTRA_INDEX_PATH.read_text(encoding="utf-8"))


# ─── 1. Unmapped tag emits UserWarning (RED before Option A fix) ─────────────

def test_detector_folder_for_unknown_tag_warns() -> None:
    """An unmapped non-empty tag must surface as a ``UserWarning``.

    This is the V126-04 acceptance criterion #1/#3 from the brief: pre-fix
    behaviour was a silent empty-string return; post-fix, callers (and
    downstream tests) get a visible signal so that
    ``_DETECTOR_TAG_TO_FOLDER`` can be kept in sync with the data feed.
    """
    from scripts.rag.build_spectra_index import _detector_folder_for

    with pytest.warns(UserWarning, match=r"unmapped detector tag"):
        result = _detector_folder_for("NonExistentTag_XYZ")
    assert result == "", (
        "Caller-side contract: unknown tags still return '' so existing "
        "downstream code can choose to fall back gracefully; the warning "
        "is the visible signal."
    )


# ─── 2. Known mapped tag returns canonical folder (no warn) ──────────────────

def test_detector_folder_for_known_tag_no_warn() -> None:
    """Known tags must keep returning their canonical folder, no warning."""
    from scripts.rag.build_spectra_index import _detector_folder_for

    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)  # any UserWarning -> raises
        assert _detector_folder_for("Gamma-1S") == "detectors/Gamma-1S/"
        assert _detector_folder_for("Handy_LaBr") == "detectors/Handy_LaBr/"
        # Legacy Gamma-1S synonym still maps to Gamma-1S (taxonomy lock 2026-06-05).
        assert _detector_folder_for("Gamma-1S") == "detectors/Gamma-1S/"


# ─── 3. Empty / whitespace inputs return "" silently (no warn) ───────────────

def test_detector_folder_for_empty_input_silent() -> None:
    """Empty string / whitespace inputs are a legitimate "no detector tag"
    case (e.g. a record whose ingestion did not populate the field); they
    must NOT trigger the unmapped-tag warning."""
    from scripts.rag.build_spectra_index import _detector_folder_for

    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        assert _detector_folder_for("") == ""
        assert _detector_folder_for("   ") == ""


# ─── 4. detector_folder non-empty for every SPECTRA_INDEX record ─────────────

def test_detector_folder_field_nonempty_in_all_records(
    spectra_index: dict,
) -> None:
    """SPECTRA_INDEX schema 0.2: ``detector_folder`` must be non-empty for
    every record. Key-presence is asserted elsewhere; this test catches the
    "empty string slipped through" case that V126-04 makes RED.
    """
    empty: list[str] = []
    for i, rec in enumerate(spectra_index["spectra"]):
        if not rec.get("detector_folder"):
            sid = rec.get("spectrum_id", f"[index {i}]")
            empty.append(sid)
    assert not empty, (
        f"{len(empty)} records have empty/missing detector_folder "
        f"(first 5: {empty[:5]})"
    )
