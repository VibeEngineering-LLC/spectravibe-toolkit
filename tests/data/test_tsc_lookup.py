"""Schema/sanity tests for data/tsc_lookup.json (v1.19.0 Phase 4).

These tests verify the TSCF lookup table is well-formed but DO NOT exercise
any pipeline integration -- per PLAN_v1_18_32_to_v1_19_0_TCS_INTEGRATION.md
§3 Phase 4, the data layer is decoupled from activity-computation until
Phase 5 (deferred to v1.19.1+).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

# Repository root: <root>/tests/data/test_tsc_lookup.py -> parents[2]
REPO_ROOT = Path(__file__).resolve().parents[2]
LOOKUP_PATH = REPO_ROOT / "data" / "tsc_lookup.json"
NUCLIDES_PATH = REPO_ROOT / "data" / "nuclides.json"
DOC_CORPUS_PATH = REPO_ROOT / "audit" / "_rag" / "DOC_CORPUS_INDEX.json"

# Required fields per entry (PLAN §3 Phase 4 schema)
REQUIRED_FIELDS = {
    "nuclide",
    "energy_keV",
    "geometry",
    "matrix_density_g_cm3",
    "detector_class",
    "detector_model_ref",
    "tscf",
    "tscf_uncertainty",
    "source_doc_id",
}


@pytest.fixture(scope="module")
def lookup() -> dict:
    return json.loads(LOOKUP_PATH.read_text(encoding="utf-8"))


def test_tsc_lookup_json_valid(lookup: dict) -> None:
    """File loads, schema version field present and matches v1 marker."""
    assert lookup.get("$schema") == "tsc_lookup_v1", (
        "expected $schema == 'tsc_lookup_v1' marker for v1.19.0 Phase 4 data layer"
    )
    assert isinstance(lookup.get("entries"), list), "entries must be a list"
    assert len(lookup["entries"]) > 0, "tsc_lookup.json must not be empty"
    # entries_count metadata consistency
    if "entries_count" in lookup:
        assert lookup["entries_count"] == len(lookup["entries"]), (
            "entries_count metadata must match len(entries)"
        )


def test_tsc_lookup_all_entries_have_required_fields(lookup: dict) -> None:
    """Every entry carries the 9 required fields per PLAN §3 Phase 4."""
    missing_report: list[tuple[int, set[str]]] = []
    for idx, entry in enumerate(lookup["entries"]):
        missing = REQUIRED_FIELDS - set(entry.keys())
        if missing:
            missing_report.append((idx, missing))
    assert not missing_report, (
        f"{len(missing_report)} entries missing required fields: "
        f"{missing_report[:3]}{'...' if len(missing_report) > 3 else ''}"
    )


def test_tsc_lookup_all_entries_hpge_coaxial(lookup: dict) -> None:
    """All v1.19.0 entries are detector_class == 'HPGe-coaxial' (PLAN §3 Phase 4
    coverage matrix). Later versions may add other classes; this guard catches
    accidental cross-class entries.
    """
    wrong = [
        (i, e["detector_class"])
        for i, e in enumerate(lookup["entries"])
        if e.get("detector_class") != "HPGe-coaxial"
    ]
    assert not wrong, f"non-HPGe-coaxial entries found in v1.19.0 lookup: {wrong[:5]}"
    # also assert the supported-list metadata mentions it
    assert "HPGe-coaxial" in lookup.get("detector_class_supported", []), (
        "detector_class_supported metadata must list 'HPGe-coaxial' for v1.19.0"
    )


def test_tsc_lookup_geometries_in_enum(lookup: dict) -> None:
    """Every entry.geometry is a member of GeometryClass."""
    from scripts.gamma.data.geometry_classes import GeometryClass

    allowed = {g.value for g in GeometryClass}
    bad = [
        (i, e["geometry"])
        for i, e in enumerate(lookup["entries"])
        if e.get("geometry") not in allowed
    ]
    assert not bad, (
        f"entries with geometry NOT in GeometryClass enum: {bad[:5]}; "
        f"allowed={sorted(allowed)}"
    )


def test_tsc_lookup_tscf_in_plausible_range(lookup: dict) -> None:
    """TSCF must be in [0.5, 2.0] -- sanity gate.

    All-cascade emitters in close geometry can reach TSCF ~1.27 (Cs-134 801.9
    keV PPAQ per Giubrone 2016 Table 7). 0.5 lower bound guards against
    accidental insertion of inverse correction; 2.0 upper bound guards against
    typos and unphysical entries.
    """
    bad = []
    for i, e in enumerate(lookup["entries"]):
        tscf = e.get("tscf")
        if tscf is None or not (0.5 <= float(tscf) <= 2.0):
            bad.append((i, e.get("nuclide"), e.get("energy_keV"), tscf))
    assert not bad, f"entries with tscf outside [0.5, 2.0]: {bad[:5]}"


def test_tsc_lookup_source_docs_in_corpus(lookup: dict) -> None:
    """Every source_doc_id resolves to an entry in DOC_CORPUS_INDEX.json."""
    corpus = json.loads(DOC_CORPUS_PATH.read_text(encoding="utf-8"))
    corpus_entries = corpus.get("entries", {})
    # entries may be dict (keyed by doc_id) or list of {doc_id: ...} dicts
    if isinstance(corpus_entries, dict):
        valid_ids = set(corpus_entries.keys())
    else:
        valid_ids = {
            e.get("doc_id") or e.get("id")
            for e in corpus_entries
            if isinstance(e, dict)
        }

    referenced = {e.get("source_doc_id") for e in lookup["entries"]}
    missing = referenced - valid_ids
    assert not missing, (
        f"source_doc_id values not registered in DOC_CORPUS_INDEX.json: {sorted(missing)}; "
        f"corpus has {len(valid_ids)} entries"
    )


def test_tsc_lookup_tscf_uncertainty_positive(lookup: dict) -> None:
    """Uncertainty must be a positive number for every entry (sanity gate)."""
    bad = [
        (i, e.get("nuclide"), e.get("tscf_uncertainty"))
        for i, e in enumerate(lookup["entries"])
        if e.get("tscf_uncertainty") is None or float(e["tscf_uncertainty"]) <= 0
    ]
    assert not bad, f"entries with non-positive tscf_uncertainty: {bad[:5]}"


def test_tsc_lookup_entry_count_at_least_target(lookup: dict) -> None:
    """v1.19.0 Phase 4 target: 47 entries (27 Giubrone + 20 Ordonez)."""
    count = len(lookup["entries"])
    assert count >= 47, (
        f"expected >= 47 entries in v1.19.0 Phase 4 baseline, got {count}"
    )
    # Source breakdown sanity
    from collections import Counter
    by_src = Counter(e.get("source_doc_id") for e in lookup["entries"])
    assert by_src.get("giubrone-2016", 0) >= 27, (
        f"expected >= 27 Giubrone entries, got {by_src.get('giubrone-2016')}"
    )
    assert by_src.get("ordonez-2019-rpc", 0) >= 20, (
        f"expected >= 20 Ordonez entries, got {by_src.get('ordonez-2019-rpc')}"
    )
