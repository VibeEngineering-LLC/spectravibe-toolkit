# -*- coding: utf-8 -*-
"""F-070 W3 — tests for `scripts/rag/visual_similarity.py`.

Coverage matrix (see brief §"Tests" 1..12):

 1.  test_query_vector_dim_128_l2_normalized
 2.  test_query_vector_matches_canonical_encoding_cs137  (≥0.999 vs stored)
 3.  test_load_templates_all_24_records_present
 4.  test_load_templates_geometry_filter_pointlike (exactly 9)
 5.  test_score_top_k_descending_cosine
 6.  test_score_cs137_self_match_match_verdict (or 'ambiguous' for tier C)
 7.  test_score_cs137_query_vs_k40_template_not_match
 8.  test_tier_c_downweight_applied (Ba-133 adjusted = 0.70 × raw)
 9.  test_stale_reference_badge_set_for_age_15y_plus (Cs-137 MARI)
10.  test_k40_marinelli_dual_cert_dates_surface_both
11.  test_verdict_thresholds_boundaries
12.  test_empty_geometry_filter_returns_all (None filter = 24)

Plus auxiliary sanity tests that pin the public API and the stale-badge
boundary that's load-bearing for the HTML card.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from rag.visual_similarity import (  # noqa: E402
    FEATURE_VECTOR_DIM,
    STALE_REFERENCE_AGE_YEARS,
    THRESHOLD_AMBIGUOUS_LOWER,
    THRESHOLD_MATCH,
    TIER_C_DOWNWEIGHT,
    compute_query_vector,
    load_templates,
    score_against_templates,
    score_query,
)

VT_ROOT = REPO_ROOT / "audit" / "_rag" / "visual_templates"
INDEX_PATH = REPO_ROOT / "audit" / "_rag" / "VISUAL_TEMPLATES_INDEX.json"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def all_templates() -> list[dict]:
    return load_templates(index_path=INDEX_PATH)


@pytest.fixture(scope="module")
def cs137_point_raw() -> dict:
    return json.loads(
        (VT_ROOT / "_raw_ingest" / "pointlike_5cm" / "Cs-137__163_2017.spe.json").read_text(
            encoding="utf-8"
        )
    )


@pytest.fixture(scope="module")
def cs137_point_template() -> dict:
    return json.loads(
        (VT_ROOT / "pointlike_5cm" / "VT-CS137-POINT5CM-2024.json").read_text(encoding="utf-8")
    )


def _load_counts(raw_record: dict) -> tuple[np.ndarray, dict]:
    """Read the .spe file backing a raw-ingest record and return
    (counts_array, energy_calib_dict)."""
    sys.path.insert(0, str(REPO_ROOT / "scripts"))
    from gamma.io.lsrm_spe import read_lsrm_spe  # type: ignore

    spe_path = REPO_ROOT / raw_record["source_file"]
    spec = read_lsrm_spe(str(spe_path))
    counts = np.asarray(spec.counts, dtype=np.int64)
    e_cal = {"coefficients": list(spec.energy_cal)}
    return counts, e_cal


# ---------------------------------------------------------------------------
# Test 1 — encoder shape/normalisation
# ---------------------------------------------------------------------------

def test_query_vector_dim_128_l2_normalized(cs137_point_raw):
    counts, e_cal = _load_counts(cs137_point_raw)
    q = compute_query_vector(counts, e_cal)
    assert q.shape == (FEATURE_VECTOR_DIM,)
    assert q.dtype == np.float64
    assert abs(np.linalg.norm(q) - 1.0) < 1e-9


# ---------------------------------------------------------------------------
# Test 2 — bit-identical (≥0.999) regression against stored canonical
# ---------------------------------------------------------------------------

def test_query_vector_matches_canonical_encoding_cs137(cs137_point_raw, cs137_point_template):
    counts, e_cal = _load_counts(cs137_point_raw)
    q = compute_query_vector(counts, e_cal)
    stored = np.asarray(cs137_point_template["feature_vector"]["values"], dtype=np.float64)
    assert stored.shape == q.shape
    stored_unit = stored / max(float(np.linalg.norm(stored)), 1e-12)
    q_unit = q / max(float(np.linalg.norm(q)), 1e-12)
    cos = float(np.dot(stored_unit, q_unit))
    # ≥0.999 per brief; in practice should be ≥0.9999 modulo JSON 6-digit rounding.
    assert cos >= 0.999, f"encoder drift: cos={cos:.6f} vs W2 stored vector"


# ---------------------------------------------------------------------------
# Test 3 — index loads all 24
# ---------------------------------------------------------------------------

def test_load_templates_all_24_records_present(all_templates):
    # F-070 W4 task #12 (2026-06-06): operator chose «Treat as Дента-120 typo».
    # The 2 previously-pending denta_100ml templates were reclassified to
    # canonical denta_120ml/ and `pending_review` was set to False on both
    # index entries — so `load_templates()` now exposes all 24 canonical
    # templates to the similarity API. The `pending_review` filter in
    # load_templates() remains in place as a defensive guard for any future
    # USER-BLOCKED records that may be added.
    assert len(all_templates) == 24
    ids = {t["template_id"] for t in all_templates}
    # cross-check with index, minus any pending_review entries (currently 0)
    idx = json.loads(INDEX_PATH.read_text(encoding="utf-8"))
    expected = {e["template_id"] for e in idx["entries"] if not e.get("pending_review")}
    assert ids == expected
    # And confirm INDEX itself has 24 entries (no entries dropped on rename)
    assert len(idx["entries"]) == 24


# ---------------------------------------------------------------------------
# Test 4 — geometry filter (pointlike_5cm has 9 entries per W2)
# ---------------------------------------------------------------------------

def test_load_templates_geometry_filter_pointlike():
    point_only = load_templates(geometry_class="pointlike_5cm", index_path=INDEX_PATH)
    assert len(point_only) == 9
    assert all(t["geometry_class"] == "pointlike_5cm" for t in point_only)


# ---------------------------------------------------------------------------
# Test 5 — top-K sort order
# ---------------------------------------------------------------------------

def test_score_top_k_descending_cosine(all_templates, cs137_point_raw):
    counts, e_cal = _load_counts(cs137_point_raw)
    q = compute_query_vector(counts, e_cal)
    top5 = score_against_templates(q, all_templates, top_k=5)
    assert len(top5) == 5
    raws = [r["cosine_raw"] for r in top5]
    assert raws == sorted(raws, reverse=True)


# ---------------------------------------------------------------------------
# Test 6 — Cs-137 self-match → verdict is "match" (or "ambiguous" if tier C)
# ---------------------------------------------------------------------------

def test_score_cs137_self_match_top_is_cs137(all_templates, cs137_point_raw):
    counts, e_cal = _load_counts(cs137_point_raw)
    q = compute_query_vector(counts, e_cal)
    top = score_against_templates(q, all_templates, top_k=3)
    assert top[0]["template_id"] == "VT-CS137-POINT5CM-2024"
    assert top[0]["nuclide"] == "Cs-137"
    # The point Cs-137 template is tier B (decay_age 7.4y), so raw==adjusted.
    assert top[0]["cosine_raw"] >= THRESHOLD_MATCH, (
        f"self-match raw cosine {top[0]['cosine_raw']:.4f} below match threshold"
    )
    assert top[0]["verdict"] == "match"


# ---------------------------------------------------------------------------
# Test 7 — Cs-137 query vs K-40 templates: K-40 never wins, never "match"
# ---------------------------------------------------------------------------

def test_score_cs137_query_vs_k40_template_not_match(all_templates, cs137_point_raw):
    counts, e_cal = _load_counts(cs137_point_raw)
    q = compute_query_vector(counts, e_cal)
    results = score_against_templates(q, all_templates, top_k=24)
    # K-40 results must NOT carry "match" verdict for a pure Cs-137 query
    k40_results = [r for r in results if r["nuclide"] == "K-40"]
    assert k40_results, "expected at least one K-40 template in the library"
    for r in k40_results:
        assert r["verdict"] in {"ambiguous", "mismatch"}, (
            f"K-40 template {r['template_id']} scored 'match' against Cs-137 query "
            f"(adjusted={r['cosine_adjusted']:.4f})"
        )


# ---------------------------------------------------------------------------
# Test 8 — Tier-C downweight
# ---------------------------------------------------------------------------

def test_tier_c_downweight_applied_ba133(all_templates, cs137_point_raw):
    counts, e_cal = _load_counts(cs137_point_raw)
    q = compute_query_vector(counts, e_cal)
    results = score_against_templates(q, all_templates, top_k=24)
    ba133 = next((r for r in results if r["template_id"] == "VT-BA133-POINT5CM-2024"), None)
    assert ba133 is not None, "Ba-133 template missing from library"
    assert ba133["tier"] == "C"
    assert ba133["cosine_adjusted"] == pytest.approx(
        ba133["cosine_raw"] * TIER_C_DOWNWEIGHT, abs=1e-12
    )


# ---------------------------------------------------------------------------
# Test 9 — stale-reference badge for ≥15y decay age (Cs-137 Marinelli)
# ---------------------------------------------------------------------------

def test_stale_reference_badge_set_for_age_15y_plus(all_templates):
    cs_mari = next(
        (t for t in all_templates if t["template_id"] == "VT-CS137-MARI0CM-2024"), None
    )
    assert cs_mari is not None
    decay = cs_mari["provenance"]["decay_age_years"]
    assert decay >= STALE_REFERENCE_AGE_YEARS

    # Re-encode against itself: query vector = its own stored vector → cos=1
    q = np.asarray(cs_mari["feature_vector"]["values"], dtype=np.float64)
    results = score_against_templates(q, [cs_mari], top_k=1)
    assert results[0]["stale_reference"] is True
    assert results[0]["decay_age_years"] >= STALE_REFERENCE_AGE_YEARS


# ---------------------------------------------------------------------------
# Test 10 — K-40 Marinelli dual cert reference dates surface both
# ---------------------------------------------------------------------------

def test_k40_marinelli_dual_cert_dates_surface_both(all_templates):
    k40 = next(
        (t for t in all_templates if t["template_id"] == "VT-K40-MARI0CM-2024"), None
    )
    assert k40 is not None
    dates = k40["cert_reference_dates"]
    # W2 raw-ingest carries 2002-05-24 + 2007-09-17 across the two constituents.
    assert "2002-05-24" in dates, f"missing 2002-05-24 in {dates}"
    assert "2007-09-17" in dates, f"missing 2007-09-17 in {dates}"
    assert len(dates) >= 2

    # And the scored result echoes them
    q = np.asarray(k40["feature_vector"]["values"], dtype=np.float64)
    scored = score_against_templates(q, [k40], top_k=1)
    surfaced = scored[0]["cert_reference_dates"]
    assert "2002-05-24" in surfaced and "2007-09-17" in surfaced


# ---------------------------------------------------------------------------
# Test 11 — verdict threshold boundaries
# ---------------------------------------------------------------------------

def _synth_template(template_id: str, tier: str, dim: int = FEATURE_VECTOR_DIM) -> dict:
    """Build a single-spike unit vector template for boundary tests.

    Using two synthetic templates with controlled cosine produces deterministic
    pre-downweight scores irrespective of the W2 library.
    """
    values = [0.0] * dim
    values[0] = 1.0
    return {
        "template_id": template_id,
        "nuclide": "SYN",
        "geometry_class": "pointlike_5cm",
        "tier": tier,
        "feature_vector": {
            "values": values, "dim": dim,
            "normalization": "l2", "encoding": "synthetic",
        },
        "provenance": {"decay_age_years": 0.0},
        "cert_reference_dates": [],
    }


def test_verdict_thresholds_boundaries():
    """Synthesise query vectors yielding exact adjusted cosines at threshold edges.

    Using tier=A so adjusted == raw, then build a 2-D-effective vector
    (q[0]=cos_target, q[1]=sin_target, rest=0) against a [1,0,...,0] template.
    """
    tmpl = _synth_template("SYN-A", tier="A")

    def _query_with_cos(c: float) -> np.ndarray:
        v = np.zeros(FEATURE_VECTOR_DIM, dtype=np.float64)
        v[0] = c
        # build orthogonal residual so ||v||=1 → cosine == c
        if 1.0 - c * c > 0:
            v[1] = math_sqrt(1.0 - c * c)
        return v

    # boundary at 0.93 → match
    r = score_against_templates(_query_with_cos(0.93), [tmpl], top_k=1)[0]
    assert r["cosine_adjusted"] == pytest.approx(0.93, abs=1e-9)
    assert r["verdict"] == "match"

    # 0.9299 → ambiguous (just below match)
    r = score_against_templates(_query_with_cos(0.9299), [tmpl], top_k=1)[0]
    assert r["verdict"] == "ambiguous"

    # 0.85 → ambiguous (boundary inclusive on lower)
    r = score_against_templates(_query_with_cos(0.85), [tmpl], top_k=1)[0]
    assert r["cosine_adjusted"] == pytest.approx(0.85, abs=1e-9)
    assert r["verdict"] == "ambiguous"

    # 0.8499 → mismatch
    r = score_against_templates(_query_with_cos(0.8499), [tmpl], top_k=1)[0]
    assert r["verdict"] == "mismatch"


# (helper to avoid module-level math import in test fn)
def math_sqrt(x: float) -> float:
    import math
    return math.sqrt(x)


# ---------------------------------------------------------------------------
# Test 12 — None filter returns all 24
# ---------------------------------------------------------------------------

def test_empty_geometry_filter_returns_all():
    # F-070 W4 task #12 (2026-06-06): the 2 previously-pending denta_100ml
    # templates were reclassified to denta_120ml and `pending_review` cleared.
    # load_templates(geometry_class=None) now returns all 24 canonical entries.
    # See test_load_templates_all_24_records_present for full rationale.
    all_t = load_templates(geometry_class=None, index_path=INDEX_PATH)
    assert len(all_t) == 24


# ---------------------------------------------------------------------------
# Bonus 13 — end-to-end score_query convenience wrapper
# ---------------------------------------------------------------------------

def test_score_query_end_to_end_cs137(cs137_point_raw):
    counts, e_cal = _load_counts(cs137_point_raw)
    top = score_query(counts, e_cal, geometry_class="pointlike_5cm", top_k=3,
                      index_path=INDEX_PATH)
    assert len(top) == 3
    assert top[0]["template_id"] == "VT-CS137-POINT5CM-2024"
    assert top[0]["verdict"] == "match"


# ---------------------------------------------------------------------------
# Bonus 14 — stale boundary at exactly 15.0y is inclusive
# ---------------------------------------------------------------------------

def test_stale_reference_boundary_inclusive():
    """`decay_age_years >= 15.0` is the rule — exactly 15.0 must fire stale=True."""
    values = [0.0] * FEATURE_VECTOR_DIM
    values[0] = 1.0
    tmpl_15 = {
        "template_id": "SYN-15Y", "nuclide": "SYN", "geometry_class": "pointlike_5cm",
        "tier": "B",
        "feature_vector": {"values": values, "dim": FEATURE_VECTOR_DIM,
                            "normalization": "l2", "encoding": "synthetic"},
        "provenance": {"decay_age_years": 15.0},
        "cert_reference_dates": [],
    }
    tmpl_149 = dict(tmpl_15)
    tmpl_149["template_id"] = "SYN-14.9Y"
    tmpl_149["provenance"] = {"decay_age_years": 14.9}

    q = np.zeros(FEATURE_VECTOR_DIM, dtype=np.float64); q[0] = 1.0
    r15 = score_against_templates(q, [tmpl_15], top_k=1)[0]
    r149 = score_against_templates(q, [tmpl_149], top_k=1)[0]
    assert r15["stale_reference"] is True
    assert r149["stale_reference"] is False
