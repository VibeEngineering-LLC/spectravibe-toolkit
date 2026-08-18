# -*- coding: utf-8 -*-
"""Schema, integrity and stability tests for ``audit/_rag/SPECTRA_INDEX.json``.

Coverage:

1. Schema validation — every spectrum record has the required keys; tier
   is in {A, B, C, D}.
2. F-115 compliance — every ``rel_path`` starts with ``<LSRM>\\``.
3. ID stability — ``spectrum_id``s are unique, follow the
   ``RAG-SPEC-NNNN`` pattern and the assignment order matches sha256-ASC.
4. Dedup correctness — every id in ``duplicates[].spectrum_ids`` is
   present in ``spectra[]``; duplicates only appear for ≥ 2 distinct
   paths sharing the same sha256.
5. Index integrity — every id under ``indexes.*`` exists in
   ``spectra[]``.
6. Methodology cross-link sanity — every id in
   ``linked_rag_methodology`` (per spectrum) refers to an entry in
   ``audit/_rag/RAG_INDEX.json`` if non-empty.
7. Drift cohort sanity — cohorts have ≥ 3 distinct passport years AND
   all members are A-tier with single-nuclide passport.

Generator: ``scripts/rag/build_spectra_index.py``.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SPECTRA_INDEX_PATH = REPO_ROOT / "audit" / "_rag" / "SPECTRA_INDEX.json"
RAG_INDEX_PATH = REPO_ROOT / "audit" / "_rag" / "RAG_INDEX.json"

# Required keys on each ``spectra[*]`` record. Matches the design §2
# schema (with engineering additions: ``parse_error``, ``quality_flags``,
# ``priority`` for diagnostic continuity).
REQUIRED_RECORD_KEYS = {
    "spectrum_id", "sha256", "rel_path", "quality_tier",
    "detector_tag", "detector_id_header", "geometry",
    "live_time_s", "real_time_s", "dead_time_frac", "channels",
    "total_counts", "acq_started_at",
    "energy_cal", "energy_cal_source", "fwhm_cal", "fwhm_model",
    "energy_range_keV", "passport", "tags", "use_cases",
    "linked_rag_methodology", "quality_flags",
}

VALID_TIERS = {"A", "B", "C", "D"}
SPEC_ID_RE = re.compile(r"^RAG-SPEC-(\d{4,})$")


# ─── Module-level fixtures ───────────────────────────────────────────────────

@pytest.fixture(scope="module")
def spectra_index() -> dict:
    if not SPECTRA_INDEX_PATH.exists():
        pytest.skip(
            f"{SPECTRA_INDEX_PATH} not present — generator has not been run."
        )
    return json.loads(SPECTRA_INDEX_PATH.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def rag_index() -> dict:
    return json.loads(RAG_INDEX_PATH.read_text(encoding="utf-8"))


# ─── 1. Schema validation ────────────────────────────────────────────────────

def test_top_level_keys(spectra_index: dict) -> None:
    for k in (
        "schema_version", "generator", "lsrm_root_placeholder",
        "tier_definitions", "spectra", "indexes", "duplicates",
    ):
        assert k in spectra_index, f"top-level key missing: {k}"
    assert spectra_index["schema_version"] == "0.2"  # F-300 W3+W4: bumped from 0.1 (detector_folder field added)
    assert spectra_index["lsrm_root_placeholder"] == "<LSRM>"
    assert set(spectra_index["tier_definitions"]) == VALID_TIERS


def test_record_schema(spectra_index: dict) -> None:
    spectra = spectra_index["spectra"]
    assert isinstance(spectra, list)
    assert len(spectra) > 0, "SPECTRA_INDEX is empty"
    missing_summary: dict[str, int] = {}
    for r in spectra:
        missing = REQUIRED_RECORD_KEYS - set(r)
        if missing:
            for k in missing:
                missing_summary[k] = missing_summary.get(k, 0) + 1
        assert r["quality_tier"] in VALID_TIERS, (
            f"{r['spectrum_id']}: invalid tier {r['quality_tier']!r}"
        )
    assert not missing_summary, (
        f"records missing required keys (key→count): {missing_summary}"
    )


# ─── 2. F-115 compliance ─────────────────────────────────────────────────────

def test_f115_all_paths_use_placeholder(spectra_index: dict) -> None:
    leaks: list[str] = []
    for r in spectra_index["spectra"]:
        if not r["rel_path"].startswith("<LSRM>\\"):
            leaks.append(f"{r['spectrum_id']}: {r['rel_path']!r}")
    # Also check the dedup table.
    for d in spectra_index["duplicates"]:
        for p in d.get("paths", []):
            if not p.startswith("<LSRM>\\"):
                leaks.append(f"dup-{d['sha256'][:12]}: {p!r}")
    assert not leaks, (
        "F-115 violation — rel_path entries leak operator-absolute paths "
        f"(showing first 5): {leaks[:5]}"
    )


# ─── 3. ID stability & ordering ──────────────────────────────────────────────

def test_spectrum_ids_unique_and_well_formed(spectra_index: dict) -> None:
    ids = [r["spectrum_id"] for r in spectra_index["spectra"]]
    # Unique.
    assert len(ids) == len(set(ids)), "duplicate spectrum_id values present"
    # Well-formed.
    bad = [sid for sid in ids if not SPEC_ID_RE.match(sid)]
    assert not bad, f"malformed spectrum_ids: {bad[:5]}"


def test_id_assignment_matches_sha256_asc(spectra_index: dict) -> None:
    """`RAG-SPEC-NNNN` numeric order MUST match sha256-ASC sort of records.

    This is the spec contract from the design proposal §2 — re-runs
    preserve IDs by sha256 lookup. The smallest sha256 gets the smallest
    ID on the genesis run, and stays put thereafter.
    """
    spectra = spectra_index["spectra"]
    sorted_by_sha = sorted(spectra, key=lambda r: r["sha256"])
    sorted_by_id_num = sorted(
        spectra,
        key=lambda r: int(SPEC_ID_RE.match(r["spectrum_id"]).group(1)),
    )
    sha_order_ids = [r["spectrum_id"] for r in sorted_by_sha]
    id_order_ids = [r["spectrum_id"] for r in sorted_by_id_num]
    assert sha_order_ids == id_order_ids, (
        "spectrum_id numeric order does not follow sha256-ASC — ID "
        "stability contract broken."
    )


# ─── 4. Dedup correctness ────────────────────────────────────────────────────

def test_duplicates_consistent(spectra_index: dict) -> None:
    all_ids = {r["spectrum_id"] for r in spectra_index["spectra"]}
    all_shas = {r["sha256"]: r["spectrum_id"] for r in spectra_index["spectra"]}
    for d in spectra_index["duplicates"]:
        sha = d["sha256"]
        assert sha in all_shas, f"duplicate group references unknown sha {sha}"
        for sid in d["spectrum_ids"]:
            assert sid in all_ids, f"duplicate group references unknown id {sid}"
        # By construction the canonical winner is the only id in this list,
        # but the paths list MUST have ≥ 2 entries.
        assert len(d["paths"]) >= 2, (
            f"dup group {sha[:12]} has < 2 paths — should not be a dup"
        )
        assert d["spectrum_ids"] == [all_shas[sha]], (
            f"dup group {sha[:12]} canonical id mismatch"
        )


# ─── 5. Index integrity ──────────────────────────────────────────────────────

def test_indexes_all_ids_resolve(spectra_index: dict) -> None:
    all_ids = {r["spectrum_id"] for r in spectra_index["spectra"]}
    idx = spectra_index["indexes"]
    for axis in ("by_nuclide", "by_detector", "by_geometry",
                 "by_quality_tier", "by_passport_year"):
        for key, ids in idx[axis].items():
            unknown = [sid for sid in ids if sid not in all_ids]
            assert not unknown, (
                f"indexes.{axis}[{key!r}] references unknown ids: "
                f"{unknown[:5]}"
            )
            # And lists are sorted+unique (deterministic re-run contract).
            assert ids == sorted(set(ids)), (
                f"indexes.{axis}[{key!r}] is not sorted-unique"
            )
    for cohort in idx["drift_cohorts"]:
        unknown = [sid for sid in cohort["members"] if sid not in all_ids]
        assert not unknown, (
            f"drift cohort {cohort['cohort_id']} references unknown ids: "
            f"{unknown[:5]}"
        )


# ─── 6. Methodology cross-link sanity ────────────────────────────────────────

def test_methodology_cross_links_resolve(
    spectra_index: dict, rag_index: dict,
) -> None:
    """Every non-empty ``linked_rag_methodology`` entry must point at a
    real ``RAG-NNN`` (methodology side).

    Most records have empty lists in the genesis index — the cross-link
    is added manually as use cases activate. We only validate that
    *non-empty* lists are well-formed.
    """
    known_methodology = set(rag_index.get("entries", {}).keys())
    leaks: list[tuple[str, str]] = []
    for r in spectra_index["spectra"]:
        for mid in r.get("linked_rag_methodology") or []:
            if mid not in known_methodology:
                leaks.append((r["spectrum_id"], mid))
    assert not leaks, (
        f"linked_rag_methodology entries do not resolve in RAG_INDEX: "
        f"{leaks[:5]}"
    )


# ─── 7. Drift cohort sanity ──────────────────────────────────────────────────

def test_drift_cohort_invariants(spectra_index: dict) -> None:
    spec_by_id = {r["spectrum_id"]: r for r in spectra_index["spectra"]}
    for cohort in spectra_index["indexes"]["drift_cohorts"]:
        assert len(cohort["year_span"]) >= 3, (
            f"{cohort['cohort_id']}: year_span < 3 should not be a cohort"
        )
        for sid in cohort["members"]:
            rec = spec_by_id.get(sid)
            assert rec is not None, f"unknown member {sid}"
            assert rec["quality_tier"] == "A", (
                f"{cohort['cohort_id']} member {sid} is not A-tier"
            )
            nucs = {p.get("nuclide") for p in rec["passport"] if p.get("nuclide")}
            assert nucs == {cohort["nuclide"]}, (
                f"{cohort['cohort_id']} member {sid} passport nuclides "
                f"{nucs} != cohort nuclide {cohort['nuclide']}"
            )
            assert rec["detector_tag"] == cohort["detector_tag"], (
                f"{cohort['cohort_id']} member {sid} detector mismatch"
            )


# ─── 8. Tier definitions + tier counts non-zero invariant ────────────────────

def test_a_tier_population_nonzero(spectra_index: dict) -> None:
    """The Spectrum-RAG exists primarily to enable F-QC-01 sweep on
    real fixtures. The A-tier population must therefore be non-empty
    (otherwise the index has no use).

    The current LSRM tree yields 166 A-tier records per the v1.22.x
    eval; we use ≥ 50 as a loose, future-proof floor to detect
    catastrophic regression of the passport-extraction pipeline without
    being brittle to small library changes.
    """
    n_a = len(spectra_index["indexes"]["by_quality_tier"]["A"])
    assert n_a >= 50, (
        f"A-tier population dropped to {n_a} (< 50) — likely a regression "
        f"in the passport parser or scope-glob enumeration."
    )
