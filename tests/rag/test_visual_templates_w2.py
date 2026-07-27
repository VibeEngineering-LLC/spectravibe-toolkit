"""Tests for F-070 Wave 2 — canonical template extraction + VISUAL_TEMPLATES_INDEX.

Scope (per brief §"Exit criteria"):
- Schema validity of every canonical template (required fields, correct types).
- Cert validation gate: ≥3 nuclides incl. Cs-137 + K-40 must hit ≤30% residual.
- Feature vector dim invariant (declared dim = 128).
- Drift-isolation: NO Поверка-2016 Маринелли / Точка-25см records in canonical dirs.
- Tier-C presence: Ba-133 record must exist, tier="C", decay_age_years > 15.
- Geometry split correctness: denta_100ml vs denta_120ml separated per D2.
- Index skeleton structure: ≥10 entries with cross-link to constituent raw-ingest paths.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
VT_ROOT = REPO_ROOT / "audit" / "_rag" / "visual_templates"
INDEX_PATH = REPO_ROOT / "audit" / "_rag" / "VISUAL_TEMPLATES_INDEX.json"
LOG_PATH = VT_ROOT / "_extraction_log.json"

CANONICAL_GEOMETRIES = ("pointlike_5cm", "marinelli_0cm", "petri_60ml", "denta_100ml", "denta_120ml")
REQUIRED_TOP_FIELDS = {
    "template_id", "nuclide", "geometry_class", "detector_class", "source_epoch",
    "tier", "characteristic_lines", "compton_continuum", "background_flags",
    "feature_vector", "qc", "provenance", "schema_version",
}
REQUIRED_PROVENANCE_FIELDS = {
    "source_files", "source_dirs", "measurement_dates", "certificate_ref",
    "certificate_activity_Bq", "certificate_reference_date", "decay_age_years",
    "geometry_field_value",
}


def _load_all_canonical() -> list[dict]:
    out = []
    for geom in CANONICAL_GEOMETRIES:
        d = VT_ROOT / geom
        if not d.exists():
            continue
        for jp in sorted(d.glob("VT-*.json")):
            out.append(json.loads(jp.read_text(encoding="utf-8")))
    return out


def test_canonical_templates_min_count_and_schema():
    """≥10 canonical templates exist; all have required top-level + provenance fields."""
    templates = _load_all_canonical()
    assert len(templates) >= 10, f"expected ≥10 canonical templates, got {len(templates)}"
    for t in templates:
        missing = REQUIRED_TOP_FIELDS - set(t)
        assert not missing, f"{t.get('template_id')} missing top-level fields: {missing}"
        prov_missing = REQUIRED_PROVENANCE_FIELDS - set(t["provenance"])
        assert not prov_missing, f"{t['template_id']} missing provenance fields: {prov_missing}"
        assert t["schema_version"] == "0.2"  # F-070 W4 S0 bump 0.1→0.2 (2026-06-05)
        assert t["geometry_class"] in CANONICAL_GEOMETRIES
        assert t["detector_class"] == "Gamma-1S"
        assert t["tier"] in ("A", "B", "C")
        assert isinstance(t["characteristic_lines"], list) and len(t["characteristic_lines"]) >= 1
        for ln in t["characteristic_lines"]:
            for k in ("energy_keV", "intensity", "expected_fwhm_keV", "area_relative"):
                assert k in ln, f"{t['template_id']} line missing {k}"


def test_feature_vector_dim_128_and_l2_normalized():
    """Declared dim = 128; values L2-normalized (or zero); channel grid present."""
    templates = _load_all_canonical()
    for t in templates:
        fv = t["feature_vector"]
        assert fv["dim"] == 128, f"{t['template_id']} feature vector dim != 128 (got {fv['dim']})"
        assert len(fv["values"]) == 128
        assert len(fv["channel_grid_keV"]) == 128
        assert fv["normalization"] == "l2"
        assert fv["encoding"] == "peak_emphasis_log_continuum_suppressed"
        norm = math.sqrt(sum(v * v for v in fv["values"]))
        # Either exactly zero (degenerate) or unit norm within rounding (6-digit precision)
        assert norm == 0 or abs(norm - 1.0) < 5e-3, f"{t['template_id']} feature vector not L2 (norm={norm})"


def test_cert_validation_gate_cs137_k40_plus_one():
    """At least 3 nuclides incl. Cs-137 and K-40 pass the 30% residual gate."""
    log = json.loads(LOG_PATH.read_text(encoding="utf-8"))
    cert_rows = log["cert_validation_summary"]
    passed = [r for r in cert_rows if r.get("gate_pass_30pct") is True]
    passed_nuclides = {r["nuclide"] for r in passed}
    assert "Cs-137" in passed_nuclides, f"Cs-137 missing from cert-gate PASS set ({passed_nuclides})"
    assert "K-40" in passed_nuclides, f"K-40 missing from cert-gate PASS set ({passed_nuclides})"
    assert len(passed_nuclides) >= 3, f"need ≥3 distinct nuclides passing, got {passed_nuclides}"


def test_drift_isolation_zero_poverka_marinelli_or_t25cm_in_canonical():
    """ZERO canonical templates may trace to Поверка-2016 Маринелли / Точка-25см."""
    templates = _load_all_canonical()
    forbidden_markers = ("Поверка-2016/Маринелли", "Поверка-2016/Точка")
    for t in templates:
        for sf in t["provenance"]["source_files"]:
            for m in forbidden_markers:
                assert m not in sf, (
                    f"{t['template_id']} traces to drift-only source '{sf}' "
                    f"(marker '{m}') — must live under _drift_study/ only"
                )


def test_ba133_tier_c_with_decay_age():
    """Ba-133 canonical template must exist, tier='C', decay_age_years > 15."""
    templates = _load_all_canonical()
    ba = [t for t in templates if t["nuclide"] == "Ba-133"]
    assert ba, "Ba-133 canonical template missing"
    t = ba[0]
    assert t["tier"] == "C", f"Ba-133 tier expected 'C', got {t['tier']}"
    age = t["provenance"]["decay_age_years"]
    assert age is not None and age > 15.0, f"Ba-133 decay_age_years expected > 15, got {age}"
    rationale = t["provenance"].get("tier_rationale", "")
    assert "Ba-133" in rationale and "D3" in rationale


def test_denta_split_100ml_vs_120ml():
    """denta_100ml vs denta_120ml split — historical evidence vs canonical.

    2026-06-06: operator reclassified 2 Поверка-2016 «Дента-100» records as
    Дента-120 typo (task #12 «Treat as Дента-120 typo»). After reclassification,
    denta_100ml/ is empty by design — folder retained for future legitimate
    denta_100ml records if they ever arrive (non-ЛСРМ vessel scenario).

    For reclassified denta_120ml entries the provenance.geometry_field_value
    still legitimately reads "Дента-100" (historical evidence trail — the
    original cert label preserved as part of the provenance audit chain).
    Only the 2024 canonical denta_120ml records carry a "Дента-120" marker.
    """
    d100 = VT_ROOT / "denta_100ml"
    d120 = VT_ROOT / "denta_120ml"
    assert d100.exists() and d120.exists(), "both denta_100ml and denta_120ml dirs must exist"
    n100 = list(d100.glob("VT-*.json"))
    n120 = list(d120.glob("VT-*.json"))
    # denta_100ml: zero allowed after 2026-06-06 reclassification
    assert len(n100) >= 0, f"denta_100ml unexpected: {len(n100)}"
    # denta_120ml: must have ≥1 (original 4 from W2 + 2 reclassified from W4 task #12 = 6)
    assert len(n120) >= 1, f"denta_120ml expected ≥1 template, got {len(n120)}"
    # If anything is still in denta_100ml, it must legitimately mention Дента-100
    for jp in n100:
        t = json.loads(jp.read_text(encoding="utf-8"))
        assert t["geometry_class"] == "denta_100ml"
        geom_fields = t["provenance"]["geometry_field_value"]
        assert any("Дента-100" in g for g in geom_fields), (
            f"{t['template_id']} in denta_100ml dir but geometry_field={geom_fields}"
        )
    for jp in n120:
        t = json.loads(jp.read_text(encoding="utf-8"))
        assert t["geometry_class"] == "denta_120ml"
        geom_fields = t["provenance"]["geometry_field_value"]
        reclass = t.get("__operator_reclassification_resolved")
        if reclass:
            # 2026-06-06 task #12 reclassified records: original Дента-100
            # marker preserved as historical evidence in provenance.
            assert any("Дента-100" in g for g in geom_fields), (
                f"{t['template_id']} reclassified but geometry_field "
                f"missing original Дента-100 marker: {geom_fields}"
            )
        else:
            # Pure-canonical denta_120ml records (e.g. 2024 epoch) carry the
            # Дента-120 marker.
            assert any("Дента-120" in g for g in geom_fields), (
                f"{t['template_id']} in denta_120ml dir but geometry_field={geom_fields}"
            )


def test_visual_templates_index_skeleton():
    """Index has ≥10 entries cross-linking to constituent raw-ingest paths."""
    assert INDEX_PATH.exists(), f"VISUAL_TEMPLATES_INDEX.json missing at {INDEX_PATH}"
    idx = json.loads(INDEX_PATH.read_text(encoding="utf-8"))
    assert idx["_meta"]["schema_version"] == "0.2"  # F-070 W4 S0 bump 0.1→0.2 (2026-06-05)
    assert idx["_meta"]["rag_id"] == "RAG-047"
    entries = idx["entries"]
    assert len(entries) >= 10, f"expected ≥10 index entries, got {len(entries)}"
    for e in entries:
        for k in ("template_id", "path", "nuclide", "geometry_class", "detector_class",
                  "tier", "n_constituents", "constituent_raw_ingest_paths", "schema_version"):
            assert k in e, f"index entry {e.get('template_id')} missing {k}"
        assert isinstance(e["constituent_raw_ingest_paths"], list) and len(e["constituent_raw_ingest_paths"]) >= 1
        for p in e["constituent_raw_ingest_paths"]:
            assert p.startswith("audit/_rag/visual_templates/_raw_ingest/"), \
                f"constituent path '{p}' must be under _raw_ingest/"


def test_geometry_inconsistency_log_present():
    """`_extraction_log.json` records D5 geometry-field cross-check inconsistencies."""
    log = json.loads(LOG_PATH.read_text(encoding="utf-8"))
    # geometry_remap_log present in _meta
    assert "geometry_remap_log" in log["_meta"], "geometry_remap_log missing from _extraction_log _meta"
    # Each template has cert_validation block
    for t in log["templates"]:
        assert "cert_validation" in t
        assert "tier" in t and "tier_rationale" in t
        for c in t["constituents"]:
            for k in ("basename", "source_file", "source_file_sha256", "geometry_field"):
                assert k in c, f"constituent of {t['template_id']} missing {k}"


def test_no_canonical_template_under_drift_dir():
    """The _drift_study/ subtree must NEVER contain canonical VT-*.json
    (only its `_raw_ingest_poverka2016/` raw records)."""
    drift_dir = VT_ROOT / "_drift_study"
    if not drift_dir.exists():
        pytest.skip("_drift_study dir not present")
    for jp in drift_dir.rglob("VT-*.json"):
        # All matches here are mistakes
        pytest.fail(f"Canonical template found under _drift_study/: {jp}")
