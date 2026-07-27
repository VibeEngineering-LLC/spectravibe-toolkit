# -*- coding: utf-8 -*-
"""F-070 W3 / v1.24.0 — visual similarity HTML card snapshot tests.

Tests:
1. test_html_card_present_after_fit_view_before_decision
2. test_html_card_renders_verdict_badge_color_classes
"""
from __future__ import annotations

import sys
from unittest.mock import patch

import pytest

sys.path.insert(0, "scripts")

from gamma.reporting.interactive_html import _build_visual_similarity_card  # noqa: E402


# ---------------------------------------------------------------------------
# Helper: minimal visual_similarity block (enabled=True with one match)
# ---------------------------------------------------------------------------

def _make_vs_block(
    verdict: str = "match",
    nuclide: str = "Cs-137",
    tier: str = "B",
    stale: bool = False,
    cosine_raw: float = 0.946,
    cosine_adj: float = 0.946,
) -> dict:
    return {
        "enabled": True,
        "policy": {
            "threshold_match": 0.93,
            "threshold_ambiguous_lower": 0.85,
            "tier_c_downweight": 0.70,
            "stale_reference_age_years": 15.0,
        },
        "query_geometry": "pointlike_5cm",
        "query_vector_dim": 128,
        "top_k": 3,
        "matches": [
            {
                "template_id": f"VT-{nuclide.replace('-', '')}-POINT5CM-2017",
                "nuclide": nuclide,
                "geometry_class": "pointlike_5cm",
                "tier": tier,
                "cosine_raw": cosine_raw,
                "cosine_adjusted": cosine_adj,
                "verdict": verdict,
                "decay_age_years": 7.4,
                "stale_reference": stale,
                "cert_reference_dates": ["2017-04-12"],
            }
        ],
        "verdict_summary": verdict,
        "verdict_summary_nuclide": nuclide if verdict != "mismatch" else None,
    }


def _make_full_report(vs_block: dict) -> dict:
    """Minimal JSON report dict with only fields needed for visual similarity card."""
    return {
        "visual_similarity": vs_block,
        "identified_nuclides": [],
        "diagnostics": {"measurement_environment": "indoor"},
    }


# ---------------------------------------------------------------------------
# 1. Card structure: present after fit-view area, before decision summary
# ---------------------------------------------------------------------------

def test_html_card_present_after_fit_view_before_decision():
    """HTML card must render with fp-vs-card class and contain the card title."""
    report = _make_full_report(_make_vs_block())
    html = _build_visual_similarity_card(report, is_background_only=False)

    assert "fp-vs-card" in html, "Missing fp-vs-card class in rendered HTML"
    assert "Визуальное сопоставление спектров" in html, "Card title missing"
    assert "SIMILARITY_POLICY.md" in html, "Policy link missing"
    # Table must be present.
    assert "fp-vs-tbl" in html, "Missing fp-vs-tbl class — table not rendered"
    # At least one match row.
    assert "Cs-137" in html, "Nuclide name missing in table"
    # Verdict banner.
    assert "fp-vs-verdict-banner" in html, "Verdict summary banner missing"


def test_html_card_absent_for_background_only():
    """Card must be empty string for background-only spectra."""
    report = _make_full_report(_make_vs_block())
    html = _build_visual_similarity_card(report, is_background_only=True)
    assert html == "", f"Expected empty string for bg-only, got: {html[:100]!r}"


def test_html_card_disabled_skeleton_when_enabled_false():
    """Disabled block must render a skeleton with reason message."""
    report = _make_full_report({"enabled": False, "reason": "templates_unavailable"})
    html = _build_visual_similarity_card(report, is_background_only=False)

    assert "fp-vs-card" in html, "Skeleton card missing fp-vs-card class"
    assert "templates_unavailable" in html, "Reason text missing in skeleton card"
    # Skeleton should be hidden by default (display:none).
    assert "display:none" in html, "Disabled skeleton must have display:none"


# ---------------------------------------------------------------------------
# 2. Verdict badge CSS classes
# ---------------------------------------------------------------------------

def test_html_card_renders_verdict_badge_color_classes():
    """Each verdict value must produce the corresponding CSS badge class."""
    for verdict, expected_class in [
        ("match",     "verdict-match"),
        ("ambiguous", "verdict-ambiguous"),
        ("mismatch",  "verdict-mismatch"),
    ]:
        report = _make_full_report(_make_vs_block(verdict=verdict))
        html = _build_visual_similarity_card(report, is_background_only=False)
        assert expected_class in html, (
            f"Expected CSS class {expected_class!r} for verdict={verdict!r}"
        )


def test_html_card_tier_badge_classes_present():
    """Tier A/B/C must each produce the correct tier badge CSS class."""
    for tier, expected_class in [
        ("A", "tier-a-badge"),
        ("B", "tier-b-badge"),
        ("C", "tier-c-badge"),
    ]:
        report = _make_full_report(_make_vs_block(tier=tier))
        html = _build_visual_similarity_card(report, is_background_only=False)
        assert expected_class in html, (
            f"Expected CSS class {expected_class!r} for tier={tier!r}"
        )


def test_html_card_stale_badge_shown_for_tier_c_stale():
    """Tier-C stale template must show stale-badge in rendered HTML."""
    report = _make_full_report(_make_vs_block(tier="C", stale=True,
                                               cosine_raw=0.90, cosine_adj=0.63))
    html = _build_visual_similarity_card(report, is_background_only=False)
    assert "stale-badge" in html, "stale-badge class missing for stale=True Tier-C"
    # Label is in Russian per F-108 compliance: "устаревший" instead of "stale ref".
    assert "устаревший" in html, "RU stale label missing in badge"


def test_html_card_downweight_tooltip_shown_for_tier_c():
    """Tier-C downweighted cosine must show a tooltip with raw value."""
    # cosine_raw=0.90, cosine_adj=0.63 (= 0.90 * 0.70)
    report = _make_full_report(
        _make_vs_block(tier="C", cosine_raw=0.9000, cosine_adj=0.6300)
    )
    html = _build_visual_similarity_card(report, is_background_only=False)
    # The tooltip text should reference the raw value.
    assert "0.9000" in html or "raw" in html.lower(), (
        "Raw cosine value or 'raw' label missing in Tier-C tooltip"
    )


def test_html_card_k40_dual_cert_epoch():
    """K-40 with dual cert dates must render the range in the Cert. epoch column."""
    vs_block = {
        "enabled": True,
        "policy": {
            "threshold_match": 0.93,
            "threshold_ambiguous_lower": 0.85,
            "tier_c_downweight": 0.70,
            "stale_reference_age_years": 15.0,
        },
        "query_geometry": "marinelli_0cm",
        "query_vector_dim": 128,
        "top_k": 3,
        "matches": [
            {
                "template_id": "VT-K40-MARINELLI0CM-2004",
                "nuclide": "K-40",
                "geometry_class": "marinelli_0cm",
                "tier": "B",
                "cosine_raw": 0.94,
                "cosine_adjusted": 0.94,
                "verdict": "match",
                "decay_age_years": 22.0,
                "stale_reference": False,
                "cert_reference_dates": ["2002-03-15", "2007-08-20"],
            }
        ],
        "verdict_summary": "match",
        "verdict_summary_nuclide": "K-40",
    }
    report = _make_full_report(vs_block)
    html = _build_visual_similarity_card(report, is_background_only=False)

    # Should show "2002–2007" range (year range).
    assert "2002" in html and "2007" in html, (
        "Dual cert date year range (2002-2007) not present in card HTML"
    )
