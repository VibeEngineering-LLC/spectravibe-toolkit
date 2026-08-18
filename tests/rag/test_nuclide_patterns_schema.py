# -*- coding: utf-8 -*-
"""Schema, integrity and cross-link tests for the F-300/W5 NPR layer.

NPR = Nuclide Pattern Record, written under
``audit/_rag/nuclide_patterns/<detector_class>/NPR-<isotope>-<class>.json``.
The registry lives in ``audit/_rag/NUCLIDE_PATTERNS_INDEX.json``; each NPR
is mirrored as a methodology-tier entry (RAG-048..RAG-071) in
``audit/_rag/RAG_INDEX.json`` and cross-links into
``audit/_rag/SPECTRA_INDEX.json`` via ``representative_spectrum_ids``.

Coverage (≥ 6 tests required by the wave brief):

1. ``test_index_top_level_keys`` — `_meta`, `patterns` present; schema_version "1.0".
2. ``test_patterns_count_at_least_20`` — index has ≥ 20 NPR entries.
3. ``test_per_pattern_file_resolves`` — every `patterns[i].file` exists on disk.
4. ``test_per_pattern_required_fields`` — each NPR file has the 12 schema fields.
5. ``test_rag_id_resolves_in_RAG_INDEX`` — each `rag_id` is a key in RAG_INDEX.
6. ``test_representative_spectrum_ids_resolve_in_SPECTRA_INDEX``.
7. ``test_f115_no_lsrm_leak_in_provenance``.
8. ``test_detector_class_matches_folder_layout``.

Generator: ``audit/_drafts/_ollama_helpers/_session_2026-06-05/gen_w5_npr_files.py``.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
NPR_INDEX_PATH = REPO_ROOT / "audit" / "_rag" / "NUCLIDE_PATTERNS_INDEX.json"
RAG_INDEX_PATH = REPO_ROOT / "audit" / "_rag" / "RAG_INDEX.json"
SPECTRA_INDEX_PATH = REPO_ROOT / "audit" / "_rag" / "SPECTRA_INDEX.json"

REQUIRED_NPR_FIELDS = {
    "pattern_id",
    "rag_id",
    "detector_class",
    "crystal_class",
    "isotope",
    "isotope_group",
    "expected_peaks_keV",
    "characteristic_features",
    "representative_spectrum_ids",
    "provenance",
    "schema_version",
    "created",
}

VALID_DETECTOR_CLASSES = {
    "Gamma-1S",
    "Handy_LaBr",
    "Handy_HPGe",
    "GP_HPGe20",
}


# ─── Fixtures ────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def npr_index() -> dict:
    if not NPR_INDEX_PATH.exists():
        pytest.skip(f"{NPR_INDEX_PATH} not present — F-300/W5 not run yet.")
    return json.loads(NPR_INDEX_PATH.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def rag_index() -> dict:
    return json.loads(RAG_INDEX_PATH.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def spectra_index() -> dict:
    if not SPECTRA_INDEX_PATH.exists():
        pytest.skip(f"{SPECTRA_INDEX_PATH} not present — F-300/W3 not run yet.")
    return json.loads(SPECTRA_INDEX_PATH.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def npr_records(npr_index: dict) -> list[tuple[Path, dict]]:
    """Load every NPR file referenced from the index. Returns (path, data)."""
    out: list[tuple[Path, dict]] = []
    for p in npr_index["patterns"]:
        fp = REPO_ROOT / p["file"]
        data = json.loads(fp.read_text(encoding="utf-8"))
        out.append((fp, data))
    return out


# ─── 1. Index top-level schema ───────────────────────────────────────────────

def test_index_top_level_keys(npr_index: dict) -> None:
    assert "_meta" in npr_index, "missing top-level '_meta'"
    assert "patterns" in npr_index, "missing top-level 'patterns'"
    meta = npr_index["_meta"]
    assert meta.get("schema_version") == "1.1", (
        f"_meta.schema_version expected '1.1', got {meta.get('schema_version')!r}"
    )
    assert meta.get("owner"), "_meta.owner must be set"
    assert isinstance(npr_index["patterns"], list)


# ─── 2. Population gate ──────────────────────────────────────────────────────

def test_patterns_count_at_least_20(npr_index: dict) -> None:
    n = len(npr_index["patterns"])
    assert n >= 20, f"expected ≥20 NPR entries, got {n}"


# ─── 3. Per-pattern file existence ───────────────────────────────────────────

def test_per_pattern_file_resolves(npr_index: dict) -> None:
    missing: list[str] = []
    for p in npr_index["patterns"]:
        fp = REPO_ROOT / p["file"]
        if not fp.exists():
            missing.append(f"{p['pattern_id']} → {p['file']}")
    assert not missing, f"NPR files missing: {missing}"


# ─── 4. Required fields per NPR ──────────────────────────────────────────────

def test_per_pattern_required_fields(npr_records: list[tuple[Path, dict]]) -> None:
    failures: dict[str, list[str]] = {}
    for fp, data in npr_records:
        missing = REQUIRED_NPR_FIELDS - set(data)
        if missing:
            failures[str(fp.name)] = sorted(missing)
        # Type sanity on the required fields
        assert isinstance(data.get("expected_peaks_keV"), list), (
            f"{fp.name}: expected_peaks_keV must be list")
        assert isinstance(data.get("characteristic_features"), list), (
            f"{fp.name}: characteristic_features must be list")
        assert isinstance(data.get("representative_spectrum_ids"), list), (
            f"{fp.name}: representative_spectrum_ids must be list")
        assert isinstance(data.get("provenance"), dict), (
            f"{fp.name}: provenance must be dict")
        assert data.get("schema_version") == "1.1", (
            f"{fp.name}: schema_version != '1.1'")
    assert not failures, f"NPRs with missing fields: {failures}"


# ─── 5. RAG_INDEX cross-link ─────────────────────────────────────────────────

def test_rag_id_resolves_in_RAG_INDEX(
    npr_index: dict, npr_records: list[tuple[Path, dict]], rag_index: dict
) -> None:
    entries = rag_index["entries"]
    unresolved: list[str] = []
    for fp, data in npr_records:
        rag_id = data["rag_id"]
        if rag_id not in entries:
            unresolved.append(f"{fp.name} → {rag_id}")
            continue
        # Cross-check: RAG_INDEX entry must mention the NPR file path
        ent = entries[rag_id]
        xrefs = ent.get("cross_references", [])
        expected_file = data.get("rag_id")  # noqa: F841 — semantics check below
        # The NPR file path under cross_references must match
        npr_path_posix = next(
            (p["file"] for p in npr_index["patterns"]
             if p["pattern_id"] == data["pattern_id"]), None,
        )
        assert npr_path_posix, f"{rag_id}: not in NUCLIDE_PATTERNS_INDEX"
        assert npr_path_posix in xrefs, (
            f"{rag_id}: RAG_INDEX cross_references missing NPR path "
            f"{npr_path_posix!r}"
        )
    assert not unresolved, f"NPR rag_ids unresolved in RAG_INDEX: {unresolved}"


# ─── 6. SPECTRA_INDEX cross-link ─────────────────────────────────────────────

def test_representative_spectrum_ids_resolve_in_SPECTRA_INDEX(
    npr_records: list[tuple[Path, dict]], spectra_index: dict
) -> None:
    valid_ids = {s["spectrum_id"] for s in spectra_index["spectra"]}
    unresolved: list[str] = []
    empty: list[str] = []
    for fp, data in npr_records:
        rep_ids = data["representative_spectrum_ids"]
        if not rep_ids:
            empty.append(fp.name)
            continue
        for sid in rep_ids:
            if sid not in valid_ids:
                unresolved.append(f"{fp.name} → {sid}")
    assert not empty, (
        f"NPRs with empty representative_spectrum_ids (must have ≥1): {empty}"
    )
    assert not unresolved, (
        f"NPR representative_spectrum_ids unresolved in SPECTRA_INDEX: "
        f"{unresolved}"
    )


# ─── 7. F-115 anonymization ──────────────────────────────────────────────────

def test_f115_no_lsrm_leak_in_provenance(
    npr_records: list[tuple[Path, dict]]
) -> None:
    leaks: list[str] = []
    for fp, data in npr_records:
        prov = data["provenance"]
        for key in ("spectra_source", "notes"):
            val = prov.get(key, "")
            if "<LSRM>" in str(val):
                leaks.append(f"{fp.name}::provenance.{key}")
    assert not leaks, f"F-115 violation — <LSRM> leaks in NPR provenance: {leaks}"


# ─── 8. Detector-class ↔ folder layout ───────────────────────────────────────

def test_detector_class_matches_folder_layout(
    npr_records: list[tuple[Path, dict]]
) -> None:
    """NPR file lives at
       audit/_rag/nuclide_patterns/<X>/...json  →  data.detector_class == X."""
    mismatches: list[str] = []
    for fp, data in npr_records:
        # parent folder name = detector_class
        folder_class = fp.parent.name
        if folder_class != data["detector_class"]:
            mismatches.append(
                f"{fp.name}: folder={folder_class!r} but "
                f"detector_class={data['detector_class']!r}"
            )
        if data["detector_class"] not in VALID_DETECTOR_CLASSES:
            mismatches.append(
                f"{fp.name}: detector_class {data['detector_class']!r} not in "
                f"{sorted(VALID_DETECTOR_CLASSES)}"
            )
    assert not mismatches, f"detector_class ↔ folder mismatches: {mismatches}"


# ─── 9. Disk ↔ index reconciliation: no orphan NPR files ────────────────────
#
# V126-05 (brief 2026-06-05_P2_V126-04-05): the existing 8 tests above are
# index-driven — every test iterates ``npr_index["patterns"]``. If a developer
# drops a new NPR JSON under ``audit/_rag/nuclide_patterns/<class>/`` and
# forgets to register it in NUCLIDE_PATTERNS_INDEX.json, the orphan file is
# silently ignored by all schema / cross-link tests. This reconciliation
# closes the disk-side blind spot.
#
# Direction caught here: file on disk but NOT in index (orphan).
# Reverse direction (index entry but file missing) is already caught by
# ``test_per_pattern_file_resolves``.

def test_no_orphan_npr_files_on_disk(npr_index: dict) -> None:
    """Every NPR JSON under ``audit/_rag/nuclide_patterns/**`` must be
    registered in ``NUCLIDE_PATTERNS_INDEX.json``.

    GREEN today (24 indexed == 24 on disk). RED the moment an unregistered
    NPR file appears under any detector-class subfolder. The index file
    itself is excluded from the disk-side set (it is not a pattern record).
    """
    indexed_files = {p["file"] for p in npr_index["patterns"]}
    npr_root = REPO_ROOT / "audit" / "_rag" / "nuclide_patterns"
    # Defensive: if the folder does not exist (skeleton-only checkout),
    # there can be no orphans by definition.
    if not npr_root.exists():
        pytest.skip(f"{npr_root.name} not present — F-300/W5 not run yet.")
    disk_files: set[str] = set()
    for f in npr_root.rglob("*.json"):
        # POSIX-style relative path matches how the index stores entries
        # (forward slashes regardless of OS).
        rel = f.relative_to(REPO_ROOT).as_posix()
        disk_files.add(rel)
    orphans = sorted(disk_files - indexed_files)
    assert not orphans, (
        f"{len(orphans)} NPR files on disk not registered in "
        f"NUCLIDE_PATTERNS_INDEX.json: {orphans}"
    )
