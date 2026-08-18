# -*- coding: utf-8 -*-
"""Structural integrity tests for the F-300/W2 detector-class folder skeletons.

Each entry in ``audit/_rag/SPECTRA_INDEX.json`` ``indexes.by_detector`` must have
a matching ``detectors/<class>/`` subtree.  Non-Gamma-1S classes also require a
``data/SPECTRA_MANIFEST.json``.  AtomSpectra is explicitly isolated per F-307 and
must NOT appear in SPECTRA_INDEX.  Raw ``.spe`` copies are forbidden in these
skeleton folders per F-115 (operator-side raw stays in LSRM tree).

Coverage (≥ 6 tests required by the wave brief):

1. ``test_all_detector_classes_have_readme`` — every by_detector key has a README.md.
2. ``test_non_gamma1s_classes_have_manifest`` — 9 non-Gamma-1S classes have SPECTRA_MANIFEST.json.
3. ``test_non_gamma1s_no_raw_lsrm`` — none of the 9 skeleton folders contain a raw_lsrm/ subtree.
4. ``test_atomspectra_not_in_spectra_index`` — F-307 isolation: AtomSpectra absent from SPECTRA_INDEX.
5. ``test_detector_folder_field_in_all_records`` — schema 0.2: every spectra[] record has detector_folder.
6. ``test_spectra_index_schema_version_0_2`` — _meta.schema_version == "0.2".
7. ``test_gamma1s_record_count_preserved`` — Gamma-1S == 394 records invariant.
8. ``test_total_spectra_count`` — total SPECTRA_INDEX records == 556 (394 + 162 new classes).

F-300 wave: W2 (Agent C) + W3+W4 (Agent A).
Written by Agent B (F-300/W7).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SPECTRA_INDEX_PATH = REPO_ROOT / "audit" / "_rag" / "SPECTRA_INDEX.json"
DETECTORS_ROOT = REPO_ROOT / "detectors"

# The 9 new skeleton classes added in F-300/W2.  Gamma-1S pre-existed.
NON_GAMMA1C_CLASSES = [
    "Handy_LaBr",
    "Handy_HPGe",
    "Handy_NaI",
    "GP_HPGe20",
    "Simple_HPGe",
    "Simple_NaI",
    "Simple_TeCd",
    "Simple_SiLi",
    "Simple_Alpha",
]


# ─── Fixtures ────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def spectra_index() -> dict:
    if not SPECTRA_INDEX_PATH.exists():
        pytest.skip(f"{SPECTRA_INDEX_PATH} not present — F-300/W3+W4 not run yet.")
    return json.loads(SPECTRA_INDEX_PATH.read_text(encoding="utf-8"))


# ─── 1. README presence for all by_detector keys ────────────────────────────

def test_all_detector_classes_have_readme(spectra_index: dict) -> None:
    """Every detector class listed in SPECTRA_INDEX.indexes.by_detector must have
    a detectors/<class>/README.md on disk."""
    missing: list[str] = []
    for cls in spectra_index["indexes"]["by_detector"]:
        readme = DETECTORS_ROOT / cls / "README.md"
        if not readme.exists():
            missing.append(f"{cls} → {readme}")
    assert not missing, (
        f"Detector classes missing README.md in detectors/ tree: {missing}"
    )


# ─── 2. SPECTRA_MANIFEST.json for non-Gamma-1S classes ─────────────────────

def test_non_gamma1s_classes_have_manifest() -> None:
    """9 new skeleton classes from F-300/W2 must each have
    detectors/<class>/data/SPECTRA_MANIFEST.json."""
    missing: list[str] = []
    for cls in NON_GAMMA1C_CLASSES:
        manifest = DETECTORS_ROOT / cls / "data" / "SPECTRA_MANIFEST.json"
        if not manifest.exists():
            missing.append(f"{cls} → {manifest}")
    assert not missing, (
        f"Non-Gamma-1S classes missing SPECTRA_MANIFEST.json: {missing}"
    )


# ─── 3. raw_lsrm/ never committed for non-Gamma-1S (F-115 + F-070 W4) ──────

def test_non_gamma1s_no_raw_lsrm() -> None:
    """F-115 — raw .spe stays operator-side. F-070 W4 (cc126ff 2026-06-06,
    operator-authorized) broadened raw_lsrm/ ingest to all 11 detectors at
    filesystem level (.gitignore pattern `detectors/*/raw_lsrm/`). Contract
    enforced here: NO raw_lsrm/ content for non-Gamma-1S classes may be tracked
    by git (committed). On-disk presence is allowed and gitignored."""
    import subprocess

    violations: list[str] = []
    for cls in NON_GAMMA1C_CLASSES:
        path_prefix = f"detectors/{cls}/raw_lsrm/"
        result = subprocess.run(
            ["git", "ls-files", "--", path_prefix],
            capture_output=True,
            text=True,
            cwd=str(REPO_ROOT),
        )
        tracked = [line for line in result.stdout.splitlines() if line.strip()]
        if tracked:
            violations.append(
                f"{cls}/raw_lsrm/ has {len(tracked)} tracked file(s) — F-115 violation"
            )
    assert not violations, (
        f"raw_lsrm/ content committed for non-Gamma-1S: {violations}"
    )


# ─── 4. AtomSpectra F-307 isolation ─────────────────────────────────────────

def test_atomspectra_not_in_spectra_index(spectra_index: dict) -> None:
    """F-307 isolation: AtomSpectra is pre-Phase-1 and must NOT appear in
    SPECTRA_INDEX.indexes.by_detector (deferred to Gamma-1S Phase 3 GA)."""
    assert "AtomSpectra" not in spectra_index["indexes"]["by_detector"], (
        "AtomSpectra found in SPECTRA_INDEX.indexes.by_detector — "
        "F-307 isolation violated (should be deferred to Phase 3 GA)"
    )


# ─── 5. detector_folder field present on all records (schema 0.2) ───────────

def test_detector_folder_field_in_all_records(spectra_index: dict) -> None:
    """SPECTRA_INDEX schema 0.2 adds a 'detector_folder' field to every record.
    Absence means W3+W4 rebuild did not complete correctly."""
    missing: list[str] = []
    for i, rec in enumerate(spectra_index["spectra"]):
        if "detector_folder" not in rec:
            sid = rec.get("spectrum_id", f"[index {i}]")
            missing.append(sid)
    assert not missing, (
        f"{len(missing)} spectra[] records missing 'detector_folder' field "
        f"(first 5: {missing[:5]})"
    )


# ─── 6. Schema version == 0.2 ───────────────────────────────────────────────

def test_spectra_index_schema_version_0_2(spectra_index: dict) -> None:
    """F-300/W3+W4 bumped _meta.schema_version from '0.1' to '0.2'."""
    got = spectra_index["_meta"]["schema_version"]
    assert got == "0.2", (
        f"SPECTRA_INDEX _meta.schema_version expected '0.2', got {got!r}"
    )


# ─── 7. Gamma-1S record count invariant ─────────────────────────────────────

def test_gamma1s_record_count_preserved(spectra_index: dict) -> None:
    """The Gamma-1S physical station must retain its 394 spectra records
    after the W3+W4 rebuild that added the detector_folder field."""
    g1c = [
        r for r in spectra_index["spectra"]
        if r.get("detector_folder") == "detectors/Gamma-1S/"
    ]
    assert len(g1c) == 394, (
        f"Gamma-1S record count expected 394, got {len(g1c)} — "
        "W3+W4 rebuild may have dropped or duplicated Gamma-1S entries"
    )


# ─── 8. Total spectra count ──────────────────────────────────────────────────

def test_total_spectra_count(spectra_index: dict) -> None:
    """After F-300/W2+W3+W4, SPECTRA_INDEX must have exactly 556 records:
    394 (Gamma-1S) + 162 (9 new skeleton classes)."""
    total = len(spectra_index["spectra"])
    assert total == 556, (
        f"SPECTRA_INDEX total records expected 556, got {total}"
    )
