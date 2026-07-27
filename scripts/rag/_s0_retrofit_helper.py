# -*- coding: utf-8 -*-
"""F-070 W4 S0 — In-place schema 0.1→0.2 retrofit helper.

One-shot script. Performs:
  1. Snapshot 24 canonical VT-*.json + INDEX (sha256 + feature_vector.values).
  2. Append schema 0.2 fields per lookup table (no mutation of existing fields
     except `schema_version` `0.1`→`0.2`).
  3. Update VISUAL_TEMPLATES_INDEX.json entries with schema 0.2 fields +
     bump _meta.schema_version + add pending_review:true for denta_100ml.
  4. Verify feature_vector.values bit-for-bit identical, schema_version=0.2.
  5. Save post-edit checksums.

Idempotent: if a file already has `schema_version: "0.2"` it skips re-adding
the new fields (uses keys-present probe).

Anti-hallucination: every (template_id → vessel_class, etc.) mapping is in
the LOOKUP dict below; no inference. Sources cited inline.

Run:
    python scripts/rag/_s0_retrofit_helper.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
VT_ROOT = REPO_ROOT / "audit" / "_rag" / "visual_templates"
INDEX_PATH = REPO_ROOT / "audit" / "_rag" / "VISUAL_TEMPLATES_INDEX.json"
PRE_CHECKSUMS_PATH = VT_ROOT / "_s0_pre_edit_checksums.txt"
POST_CHECKSUMS_PATH = VT_ROOT / "_s0_post_edit_checksums.txt"
FEATURE_VECTOR_SNAPSHOT_PATH = VT_ROOT / "_s0_feature_vector_snapshot.json"

# ---------------------------------------------------------------------------
# Lookup table — per replan §3 / inbox brief table.
# Keys are basename minus .json (i.e. template_id).
# Sources: detectors/Gamma-1S/README.md §3 (vessel canonicalization) and
# audit/_plans/F-070_W4_S0_migration_replan.md §3 field-mapping table.
# ---------------------------------------------------------------------------

LOOKUP: dict[str, dict] = {
    # marinelli_0cm/ (4) — vessel marinelli_1L per README.md §3:97 + lab convention §3:108
    "VT-CS137-MARI0CM-2024": {
        "folder": "marinelli_0cm",
        "vessel_class": "marinelli_1L",
        "useful_sample_volume_ml": 1000,
        "placement_distance_cm": 0,
        "effective_thickness_mm": [26, 2],
        "__needs_operator_review": False,
    },
    "VT-K40-MARI0CM-2024": {
        "folder": "marinelli_0cm",
        "vessel_class": "marinelli_1L",
        "useful_sample_volume_ml": 1000,
        "placement_distance_cm": 0,
        "effective_thickness_mm": [26, 2],
        "__needs_operator_review": False,
    },
    "VT-RA226-MARI0CM-2024": {
        "folder": "marinelli_0cm",
        "vessel_class": "marinelli_1L",
        "useful_sample_volume_ml": 1000,
        "placement_distance_cm": 0,
        "effective_thickness_mm": [26, 2],
        "__needs_operator_review": False,
    },
    "VT-TH232-MARI0CM-2024": {
        "folder": "marinelli_0cm",
        "vessel_class": "marinelli_1L",
        "useful_sample_volume_ml": 1000,
        "placement_distance_cm": 0,
        "effective_thickness_mm": [26, 2],
        "__needs_operator_review": False,
    },
    # denta_120ml/ (4) — vessel denta_120ml per README.md §3:99
    "VT-CS137-DENTA120ML-2024": {
        "folder": "denta_120ml",
        "vessel_class": "denta_120ml",
        "useful_sample_volume_ml": 120,
        "placement_distance_cm": 0,
        "effective_thickness_mm": [36, 2],
        "__needs_operator_review": False,
    },
    "VT-K40-DENTA120ML-2024": {
        "folder": "denta_120ml",
        "vessel_class": "denta_120ml",
        "useful_sample_volume_ml": 120,
        "placement_distance_cm": 0,
        "effective_thickness_mm": [36, 2],
        "__needs_operator_review": False,
    },
    "VT-RA226-DENTA120ML-2024": {
        "folder": "denta_120ml",
        "vessel_class": "denta_120ml",
        "useful_sample_volume_ml": 120,
        "placement_distance_cm": 0,
        "effective_thickness_mm": [36, 2],
        "__needs_operator_review": False,
    },
    "VT-TH232-DENTA120ML-2024": {
        "folder": "denta_120ml",
        "vessel_class": "denta_120ml",
        "useful_sample_volume_ml": 120,
        "placement_distance_cm": 0,
        "effective_thickness_mm": [36, 2],
        "__needs_operator_review": False,
    },
    # petri_60ml/ (5) — vessel petri_75ml + useful 60ml user lock 2026-06-05 (README.md §3:115,126)
    "VT-CS137-PETRI60ML-2024": {
        "folder": "petri_60ml",
        "vessel_class": "petri_75ml",
        "useful_sample_volume_ml": 60,
        "placement_distance_cm": 0,
        "effective_thickness_mm": [15, 2],
        "__needs_operator_review": False,
        "__detector_scope": "Gamma-1S",
    },
    "VT-K40-PETRI60ML-2024": {
        "folder": "petri_60ml",
        "vessel_class": "petri_75ml",
        "useful_sample_volume_ml": 60,
        "placement_distance_cm": 0,
        "effective_thickness_mm": [15, 2],
        "__needs_operator_review": False,
        "__detector_scope": "Gamma-1S",
    },
    "VT-RA226-PETRI60ML-2024": {
        "folder": "petri_60ml",
        "vessel_class": "petri_75ml",
        "useful_sample_volume_ml": 60,
        "placement_distance_cm": 0,
        "effective_thickness_mm": [15, 2],
        "__needs_operator_review": False,
        "__detector_scope": "Gamma-1S",
    },
    "VT-TH232-PETRI60ML-MERGED": {
        "folder": "petri_60ml",
        "vessel_class": "petri_75ml",
        "useful_sample_volume_ml": 60,
        "placement_distance_cm": 0,
        "effective_thickness_mm": [15, 2],
        "__needs_operator_review": False,
        "__detector_scope": "Gamma-1S",
    },
    "VT-MIXAMCSEUTI-PETRI60ML-2016": {
        "folder": "petri_60ml",
        "vessel_class": "petri_75ml",
        "useful_sample_volume_ml": 60,
        "placement_distance_cm": 0,
        "effective_thickness_mm": [15, 2],
        "__needs_operator_review": False,
        "__detector_scope": "Gamma-1S",
    },
    # denta_100ml/ (2) — non-ЛСРМ fallback (README.md §3:117-124)
    "VT-MIXAMCSEUTI-DENTA100ML-2016": {
        "folder": "denta_100ml",
        "vessel_class": "other_denta_100ml",
        "useful_sample_volume_ml": None,
        "placement_distance_cm": 0,
        "effective_thickness_mm": None,
        "__needs_operator_review": True,
        "lsrm_standard": False,
        "__pending_review_reason": (
            "Дента-100 reconciliation pending: user lock "
            "2026-06-05 (detectors/Gamma-1S/README.md §3 vessel canonicalization). "
            "Three candidates: (a) non-ЛСРМ vessel; (b) partial-fill "
            "120ml; (c) typo. Awaiting explicit user instruction."
        ),
    },
    "VT-TH232-DENTA100ML-2016": {
        "folder": "denta_100ml",
        "vessel_class": "other_denta_100ml",
        "useful_sample_volume_ml": None,
        "placement_distance_cm": 0,
        "effective_thickness_mm": None,
        "__needs_operator_review": True,
        "lsrm_standard": False,
        "__pending_review_reason": (
            "Дента-100 reconciliation pending: user lock "
            "2026-06-05 (detectors/Gamma-1S/README.md §3 vessel canonicalization). "
            "Three candidates: (a) non-ЛСРМ vessel; (b) partial-fill "
            "120ml; (c) typo. Awaiting explicit user instruction."
        ),
    },
    # pointlike_5cm/ (9) — point_source + 5cm distance per GEOM_TAG POINT5CM
    "VT-AM241-POINT5CM-2024": {
        "folder": "pointlike_5cm",
        "vessel_class": "point_source",
        "useful_sample_volume_ml": None,
        "placement_distance_cm": 5,
        "effective_thickness_mm": None,
        "__needs_operator_review": False,
    },
    "VT-BA133-POINT5CM-2024": {
        "folder": "pointlike_5cm",
        "vessel_class": "point_source",
        "useful_sample_volume_ml": None,
        "placement_distance_cm": 5,
        "effective_thickness_mm": None,
        "__needs_operator_review": False,
    },
    "VT-CD109-POINT5CM-2024": {
        "folder": "pointlike_5cm",
        "vessel_class": "point_source",
        "useful_sample_volume_ml": None,
        "placement_distance_cm": 5,
        "effective_thickness_mm": None,
        "__needs_operator_review": False,
    },
    "VT-CO60-POINT5CM-2024": {
        "folder": "pointlike_5cm",
        "vessel_class": "point_source",
        "useful_sample_volume_ml": None,
        "placement_distance_cm": 5,
        "effective_thickness_mm": None,
        "__needs_operator_review": False,
    },
    "VT-CS137-POINT5CM-2024": {
        "folder": "pointlike_5cm",
        "vessel_class": "point_source",
        "useful_sample_volume_ml": None,
        "placement_distance_cm": 5,
        "effective_thickness_mm": None,
        "__needs_operator_review": False,
    },
    "VT-EU152-POINT5CM-2024": {
        "folder": "pointlike_5cm",
        "vessel_class": "point_source",
        "useful_sample_volume_ml": None,
        "placement_distance_cm": 5,
        "effective_thickness_mm": None,
        "__needs_operator_review": False,
    },
    "VT-NA22-POINT5CM-2024": {
        "folder": "pointlike_5cm",
        "vessel_class": "point_source",
        "useful_sample_volume_ml": None,
        "placement_distance_cm": 5,
        "effective_thickness_mm": None,
        "__needs_operator_review": False,
    },
    "VT-TH228-POINT5CM-2024": {
        "folder": "pointlike_5cm",
        "vessel_class": "point_source",
        "useful_sample_volume_ml": None,
        "placement_distance_cm": 5,
        "effective_thickness_mm": None,
        "__needs_operator_review": False,
    },
    "VT-Y88-POINT5CM-2024": {
        "folder": "pointlike_5cm",
        "vessel_class": "point_source",
        "useful_sample_volume_ml": None,
        "placement_distance_cm": 5,
        "effective_thickness_mm": None,
        "__needs_operator_review": False,
    },
}

# Constant fields applied to every retrofitted template.
CRYSTAL_CLASS = "NaI-63x63"
STATION_OBSERVED_ON = "Gamma-1S"
NEW_SCHEMA_VERSION = "0.2"


def _sha256_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _list_canonical_files() -> list[Path]:
    files: list[Path] = []
    for tid, meta in LOOKUP.items():
        files.append(VT_ROOT / meta["folder"] / f"{tid}.json")
    return sorted(files)


def snapshot_pre(files: list[Path]) -> None:
    fv_snap: dict[str, list] = {}
    lines = []
    for f in files:
        h = _sha256_file(f)
        rel = f.relative_to(REPO_ROOT).as_posix()
        lines.append(f"{h}  {rel}")
        with f.open(encoding="utf-8") as fh:
            t = json.load(fh)
        fv_snap[t["template_id"]] = list(t["feature_vector"]["values"])
    PRE_CHECKSUMS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    FEATURE_VECTOR_SNAPSHOT_PATH.write_text(
        json.dumps(fv_snap, indent=2), encoding="utf-8"
    )
    print(f"[snapshot] pre-edit checksums + feature vectors stored ({len(files)} files).")


def retrofit_template(path: Path) -> bool:
    """Apply schema 0.2 retrofit to one VT-*.json. Returns True if mutated."""
    with path.open(encoding="utf-8") as fh:
        t = json.load(fh)
    tid = t["template_id"]
    meta = LOOKUP.get(tid)
    if meta is None:
        print(f"[skip] no LOOKUP entry for {tid}", file=sys.stderr)
        return False

    # Build new field set (only ones missing get added — idempotent).
    additions: dict = {}
    if t.get("schema_version") != NEW_SCHEMA_VERSION:
        additions["schema_version"] = NEW_SCHEMA_VERSION
    if "crystal_class" not in t:
        additions["crystal_class"] = CRYSTAL_CLASS
    if "station_observed_on" not in t:
        additions["station_observed_on"] = STATION_OBSERVED_ON
    if "vessel_class" not in t:
        additions["vessel_class"] = meta["vessel_class"]
    if "useful_sample_volume_ml" not in t:
        additions["useful_sample_volume_ml"] = meta["useful_sample_volume_ml"]
    if "placement_distance_cm" not in t:
        additions["placement_distance_cm"] = meta["placement_distance_cm"]
    if "effective_thickness_mm" not in t:
        additions["effective_thickness_mm"] = meta["effective_thickness_mm"]
    if "__needs_operator_review" not in t:
        additions["__needs_operator_review"] = meta["__needs_operator_review"]

    # __detector_scope only for petri_60ml rows (per replan §3).
    if "__detector_scope" in meta and "__detector_scope" not in t:
        additions["__detector_scope"] = meta["__detector_scope"]

    # denta_100ml fallback flags.
    if "lsrm_standard" in meta and "lsrm_standard" not in t:
        additions["lsrm_standard"] = meta["lsrm_standard"]
    if "__pending_review_reason" in meta and "__pending_review_reason" not in t:
        additions["__pending_review_reason"] = meta["__pending_review_reason"]

    # __geometry_provenance (single-layer lsrm_source_default).
    if "__geometry_provenance" not in t:
        additions["__geometry_provenance"] = [
            {
                "layer": "lsrm_source_default",
                "field": "vessel_class",
                "value": meta["vessel_class"],
                "raw": (
                    "ЛСРМ source стр. 11, "
                    "detectors/Gamma-1S/README.md §3 Canonical vessel-classes"
                ),
            }
        ]
    if "__geometry_conflicts" not in t:
        additions["__geometry_conflicts"] = []

    if not additions:
        return False

    # In-place mutation preserving existing key order (additions appended at end).
    # Strategy: rebuild dict in order existing keys + new additions, except
    # schema_version (replace in place if present, NOT append duplicate).
    new_obj: dict = {}
    for k, v in t.items():
        if k == "schema_version" and "schema_version" in additions:
            new_obj[k] = additions.pop("schema_version")
        else:
            new_obj[k] = v
    for k, v in additions.items():
        new_obj[k] = v

    # Write with sort_keys=False, indent=2, ensure_ascii=False.
    path.write_text(
        json.dumps(new_obj, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return True


def retrofit_index() -> None:
    with INDEX_PATH.open(encoding="utf-8") as fh:
        idx = json.load(fh)
    idx["_meta"]["schema_version"] = NEW_SCHEMA_VERSION
    idx["_meta"]["wave"] = "F-070-W4-S0"
    idx["_meta"]["generated_at_utc"] = datetime.now(tz=timezone.utc).isoformat(
        timespec="seconds"
    )
    for entry in idx["entries"]:
        tid = entry.get("template_id")
        meta = LOOKUP.get(tid)
        if meta is None:
            continue
        entry["schema_version"] = NEW_SCHEMA_VERSION
        entry["crystal_class"] = CRYSTAL_CLASS
        entry["station_observed_on"] = STATION_OBSERVED_ON
        entry["vessel_class"] = meta["vessel_class"]
        if meta["folder"] == "denta_100ml":
            entry["pending_review"] = True
    INDEX_PATH.write_text(
        json.dumps(idx, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"[index] updated {len(idx['entries'])} entries + _meta")


def snapshot_post(files: list[Path]) -> None:
    lines = []
    for f in files:
        h = _sha256_file(f)
        rel = f.relative_to(REPO_ROOT).as_posix()
        lines.append(f"{h}  {rel}")
    POST_CHECKSUMS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[snapshot] post-edit checksums stored ({len(files)} files).")


def verify_feature_vectors_unchanged(files: list[Path]) -> tuple[bool, list[str]]:
    fv_pre = json.loads(FEATURE_VECTOR_SNAPSHOT_PATH.read_text(encoding="utf-8"))
    failures: list[str] = []
    for f in files:
        with f.open(encoding="utf-8") as fh:
            t = json.load(fh)
        tid = t["template_id"]
        post = list(t["feature_vector"]["values"])
        pre = fv_pre.get(tid)
        if pre is None:
            failures.append(f"{tid}: no pre-snapshot")
            continue
        if pre != post:
            failures.append(f"{tid}: feature_vector mutated (len pre={len(pre)} post={len(post)})")
    return (len(failures) == 0, failures)


def verify_schema_fields(files: list[Path]) -> tuple[bool, list[str]]:
    failures: list[str] = []
    required = {
        "crystal_class",
        "station_observed_on",
        "vessel_class",
        "useful_sample_volume_ml",
        "placement_distance_cm",
        "effective_thickness_mm",
        "__geometry_provenance",
        "__geometry_conflicts",
        "__needs_operator_review",
    }
    for f in files:
        with f.open(encoding="utf-8") as fh:
            t = json.load(fh)
        tid = t["template_id"]
        if t.get("schema_version") != NEW_SCHEMA_VERSION:
            failures.append(f"{tid}: schema_version={t.get('schema_version')!r}, expected {NEW_SCHEMA_VERSION}")
        miss = required - set(t)
        if miss:
            failures.append(f"{tid}: missing {sorted(miss)}")
    return (len(failures) == 0, failures)


def main() -> int:
    files = _list_canonical_files()
    print(f"[start] {len(files)} canonical templates targeted")
    assert len(files) == 24, f"Expected 24 canonical templates, got {len(files)}"

    if not PRE_CHECKSUMS_PATH.exists():
        snapshot_pre(files)
    else:
        print(f"[snapshot] pre-edit checksums already exist at {PRE_CHECKSUMS_PATH.name}; preserving.")

    n_mutated = 0
    for f in files:
        if retrofit_template(f):
            n_mutated += 1
    print(f"[retrofit] mutated {n_mutated}/{len(files)} template files")

    retrofit_index()

    snapshot_post(files)

    fv_ok, fv_fails = verify_feature_vectors_unchanged(files)
    schema_ok, schema_fails = verify_schema_fields(files)
    if not fv_ok:
        print("[FAIL] feature vector regression:")
        for f in fv_fails:
            print(" ", f)
    if not schema_ok:
        print("[FAIL] schema field gap:")
        for f in schema_fails:
            print(" ", f)
    if fv_ok and schema_ok:
        print("[OK] feature vectors intact, schema 0.2 fields present on all 24 templates")
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
