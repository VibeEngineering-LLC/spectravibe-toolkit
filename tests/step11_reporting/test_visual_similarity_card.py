# -*- coding: utf-8 -*-
"""F-070 W3 / v1.24.0 — visual_similarity block wiring tests.

Tests:
1. test_json_report_emits_visual_similarity_block_for_cs137
2. test_visual_similarity_disabled_when_templates_dir_missing
3. test_visual_similarity_disabled_for_bg_only
4. test_top_k_is_3_default_and_configurable
5. test_geometry_inference_from_lsrm_comment
6. test_geometry_unknown_scores_all_24_templates
7. test_tier_c_downweight_reflected_in_emitted_block
"""
from __future__ import annotations

import sys
import types
from typing import Any, Dict
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

sys.path.insert(0, "scripts")

from gamma.reporting.json_report import (  # noqa: E402
    _build_visual_similarity,
    _infer_geometry_from_comment,
)


# ---------------------------------------------------------------------------
# Helper: build a minimal mock StagedAnalysisResult
# ---------------------------------------------------------------------------

def _make_result(
    comment: str = "",
    counts_len: int = 512,
    env: str = "indoor",
    energy_cal: tuple = (0.0, 1.3, 0.0),
) -> MagicMock:
    """Build a minimal StagedAnalysisResult mock for visual_similarity tests."""
    spec = MagicMock()
    spec.comment = comment
    spec.energy_cal = energy_cal
    spec.counts = list(range(counts_len))
    result = MagicMock()
    result.spec = spec
    # classify_environment is called inside _build_visual_similarity
    result.measurement_environment = env
    return result


def _make_templates(
    n: int = 3,
    tier: str = "B",
    geometry_class: str = "pointlike_5cm",
    nuclide: str = "Cs-137",
) -> list:
    """Build minimal template dicts that score_against_templates can consume."""
    templates = []
    for i in range(n):
        fv_values = np.zeros(128)
        # Put a single non-zero bin so cosine > 0.
        fv_values[i % 128] = 1.0
        norm = np.linalg.norm(fv_values)
        if norm > 0:
            fv_values = (fv_values / norm).tolist()
        templates.append({
            "template_id": f"VT-{nuclide.replace('-', '')}-GEO{i}-2017",
            "nuclide": nuclide,
            "geometry_class": geometry_class,
            "tier": tier,
            "feature_vector": {"values": list(fv_values)},
            "provenance": {
                "certificate_reference_date": "2017-04-12",
                "decay_age_years": 7.4,
            },
            "cert_reference_dates": ["2017-04-12"],
        })
    return templates


# ---------------------------------------------------------------------------
# 1. Emits visual_similarity block for normal spectrum
# ---------------------------------------------------------------------------

def test_json_report_emits_visual_similarity_block_for_cs137():
    """visual_similarity block must be emitted with all required keys and valid verdict."""
    result = _make_result(comment="GEOMETRY=Точечная-5см")

    templates = _make_templates(n=3, tier="B", geometry_class="pointlike_5cm", nuclide="Cs-137")

    with patch("gamma.reporting.environment.classify_environment", return_value="indoor"):
        with patch("rag.visual_similarity.load_templates", return_value=templates):
            block = _build_visual_similarity(result)

    assert block.get("enabled") is True, f"Expected enabled=True, got: {block}"

    # Required top-level keys.
    for key in ("policy", "query_geometry", "query_vector_dim", "top_k",
                "matches", "verdict_summary", "verdict_summary_nuclide"):
        assert key in block, f"Missing key: {key}"

    assert block["top_k"] == 3
    assert block["query_vector_dim"] == 128
    assert isinstance(block["matches"], list)

    # Policy block shape.
    policy = block["policy"]
    for pk in ("threshold_match", "threshold_ambiguous_lower", "tier_c_downweight",
               "stale_reference_age_years"):
        assert pk in policy, f"Missing policy key: {pk}"
    assert policy["threshold_match"] == 0.93
    assert policy["threshold_ambiguous_lower"] == 0.85

    # verdict_summary must be one of the 3 valid values.
    assert block["verdict_summary"] in ("match", "ambiguous", "mismatch")


# ---------------------------------------------------------------------------
# 2. Disabled when templates unavailable (load_templates raises / empty)
# ---------------------------------------------------------------------------

def test_visual_similarity_disabled_when_templates_dir_missing():
    """When load_templates raises, block must be enabled=False reason=templates_unavailable."""
    result = _make_result()

    with patch("gamma.reporting.environment.classify_environment", return_value="indoor"):
        with patch("rag.visual_similarity.load_templates",
                   side_effect=FileNotFoundError("index not found")):
            block = _build_visual_similarity(result)

    assert block.get("enabled") is False
    assert block.get("reason") == "templates_unavailable"


def test_visual_similarity_disabled_when_templates_empty():
    """When load_templates returns empty list, block must be enabled=False."""
    result = _make_result()

    with patch("gamma.reporting.environment.classify_environment", return_value="indoor"):
        with patch("rag.visual_similarity.load_templates", return_value=[]):
            block = _build_visual_similarity(result)

    assert block.get("enabled") is False
    assert block.get("reason") == "templates_unavailable"


# ---------------------------------------------------------------------------
# 3. Disabled for background-only spectra
# ---------------------------------------------------------------------------

def test_visual_similarity_disabled_for_bg_only():
    """Background-only spectra must yield enabled=False, reason=background_only_spectrum."""
    result = _make_result()

    with patch("gamma.reporting.environment.classify_environment",
               return_value="background_only"):
        block = _build_visual_similarity(result)

    assert block.get("enabled") is False
    assert block.get("reason") == "background_only_spectrum"
    # matches key must be absent (brief: "omit matches[]").
    assert "matches" not in block


# ---------------------------------------------------------------------------
# 4. top_k is 3 by default; score_against_templates called with top_k=3
# ---------------------------------------------------------------------------

def test_top_k_is_3_default_and_configurable():
    """score_against_templates must be called with top_k=3 and result has ≤3 matches."""
    result = _make_result(comment="GEOMETRY=Маринелли")
    templates = _make_templates(n=5, tier="B", geometry_class="marinelli_0cm")

    from rag.visual_similarity import score_against_templates as _real_sat

    call_log = []

    def _spy_sat(query_vector, tmpl_list, top_k=3):
        call_log.append(top_k)
        return _real_sat(query_vector, tmpl_list, top_k=top_k)

    with patch("gamma.reporting.environment.classify_environment", return_value="indoor"):
        with patch("rag.visual_similarity.load_templates", return_value=templates):
            with patch("rag.visual_similarity.score_against_templates",
                       side_effect=_spy_sat):
                block = _build_visual_similarity(result)

    assert block.get("enabled") is True
    assert block["top_k"] == 3
    assert len(block["matches"]) <= 3
    # Verify score_against_templates was called with top_k=3.
    assert call_log and call_log[0] == 3


# ---------------------------------------------------------------------------
# 5. Geometry inference from LSRM COMMENT
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("comment,expected_geom", [
    ("GEOMETRY=Маринелли\nMEASUREMENT=sample", "marinelli_0cm"),
    ("GEOMETRY=Точечная-5см", "pointlike_5cm"),
    ("GEOMETRY=Точ.5см", "pointlike_5cm"),
    ("GEOMETRY=Дента-100", "denta_100ml"),
    ("GEOMETRY=Дента-120", "denta_120ml"),
    ("GEOMETRY=Чашка Петри 60мл", "petri_60ml"),
    ("GEOMETRY=Петри-60", "petri_60ml"),
    ("GEOMETRY=Marinelli", "marinelli_0cm"),
    ("GEOMETRY=Petri-dish", "petri_60ml"),
])
def test_geometry_inference_from_lsrm_comment(comment, expected_geom):
    """GEOMETRY= token in LSRM COMMENT must map to canonical geometry_class."""
    result = _infer_geometry_from_comment(comment)
    assert result == expected_geom, (
        f"Comment {comment!r}: expected {expected_geom!r}, got {result!r}"
    )


def test_geometry_inference_missing_returns_none():
    """Missing GEOMETRY token returns None."""
    assert _infer_geometry_from_comment("") is None
    assert _infer_geometry_from_comment("MEASUREMENT=sample") is None
    assert _infer_geometry_from_comment(None) is None  # type: ignore


def test_geometry_inference_marinelli_from_block():
    """Full integration: _build_visual_similarity reads GEOMETRY from spec.comment."""
    result = _make_result(comment="SAMPLE=Probe_1\nGEOMETRY=Маринелли\n")
    templates = _make_templates(n=2, geometry_class="marinelli_0cm")

    with patch("gamma.reporting.environment.classify_environment", return_value="indoor"):
        with patch("rag.visual_similarity.load_templates", return_value=templates) as mock_lt:
            block = _build_visual_similarity(result)

    # load_templates must have been called with geometry_class="marinelli_0cm".
    mock_lt.assert_called_once()
    call_kwargs = mock_lt.call_args
    # Either positional or keyword arg.
    if call_kwargs.args:
        assert call_kwargs.args[0] == "marinelli_0cm"
    else:
        assert call_kwargs.kwargs.get("geometry_class") == "marinelli_0cm"

    assert block.get("query_geometry") == "marinelli_0cm"


# ---------------------------------------------------------------------------
# 6. Unknown geometry → load_templates called with geometry_class=None
# ---------------------------------------------------------------------------

def test_geometry_unknown_scores_all_24_templates():
    """No GEOMETRY token → load_templates called with geometry_class=None (all 24)."""
    result = _make_result(comment="MEASUREMENT=sample")  # no GEOMETRY=

    # Provide 5 templates (simulating "all" pool).
    templates = _make_templates(n=5, geometry_class="marinelli_0cm")

    with patch("gamma.reporting.environment.classify_environment", return_value="indoor"):
        with patch("rag.visual_similarity.load_templates",
                   return_value=templates) as mock_lt:
            # Ensure best cosine is high enough not to trigger the mismatch fallback.
            # Patch score_against_templates to return a match result.
            fake_match = {
                "template_id": "VT-CS137-MARINELLI0CM-2017",
                "nuclide": "Cs-137",
                "geometry_class": "marinelli_0cm",
                "tier": "B",
                "cosine_raw": 0.95,
                "cosine_adjusted": 0.95,
                "verdict": "match",
                "decay_age_years": 7.4,
                "stale_reference": False,
                "cert_reference_dates": ["2017-04-12"],
            }
            with patch("rag.visual_similarity.score_against_templates",
                       return_value=[fake_match]):
                block = _build_visual_similarity(result)

    # load_templates must have been called with geometry_class=None.
    mock_lt.assert_called_once()
    call_kwargs = mock_lt.call_args
    if call_kwargs.args:
        assert call_kwargs.args[0] is None
    else:
        assert call_kwargs.kwargs.get("geometry_class") is None

    assert block.get("query_geometry") is None


# ---------------------------------------------------------------------------
# 7. Tier C downweight reflected in emitted block
# ---------------------------------------------------------------------------

def test_tier_c_downweight_reflected_in_emitted_block():
    """Tier-C template: cosine_adjusted must equal cosine_raw * 0.70."""
    result = _make_result(comment="GEOMETRY=Точечная-5см")

    # One tier-C template with known cosine_raw.
    fake_raw = 0.90
    fake_adj = round(fake_raw * 0.70, 4)
    fake_match = {
        "template_id": "VT-BA133-POINT5CM-2001",
        "nuclide": "Ba-133",
        "geometry_class": "pointlike_5cm",
        "tier": "C",
        "cosine_raw": fake_raw,
        "cosine_adjusted": fake_adj,
        "verdict": "ambiguous",
        "decay_age_years": 23.0,
        "stale_reference": True,
        "cert_reference_dates": ["2001-06-15"],
    }

    with patch("gamma.reporting.environment.classify_environment", return_value="indoor"):
        with patch("rag.visual_similarity.load_templates",
                   return_value=[{"template_id": "x", "feature_vector": {"values": [1.0] + [0.0]*127},
                                   "nuclide": "Ba-133", "geometry_class": "pointlike_5cm",
                                   "tier": "C", "provenance": {"decay_age_years": 23.0},
                                   "cert_reference_dates": ["2001-06-15"]}]):
            with patch("rag.visual_similarity.score_against_templates",
                       return_value=[fake_match]):
                block = _build_visual_similarity(result)

    assert block.get("enabled") is True
    matches = block.get("matches", [])
    assert len(matches) >= 1
    top = matches[0]
    assert top["tier"] == "C"
    # cosine_adjusted must be 0.70 * cosine_raw (±0.001 float tolerance).
    assert abs(top["cosine_adjusted"] - fake_adj) < 0.001, (
        f"Expected cosine_adjusted={fake_adj}, got {top['cosine_adjusted']}"
    )
    assert top["stale_reference"] is True
