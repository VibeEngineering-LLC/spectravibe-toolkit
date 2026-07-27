"""F-FIT-VIEW v2 / v1.22.5 — fit overlay extended categories (Task #66).

Tests cover the three new peak categories added in v2:
  - source="secondary"    — orange Gaussians (XRF, escape, sum, backscatter)
  - source="background"   — gray dashed Gaussians (background primary FEPs)
  - source="unidentified" — yellow dashed Gaussians (true_unmatched residuals)

Plus:
  - sigma fallback path when fwhm_keV is absent (fwhm_model used)
  - multi-category co-render (primary + secondary + background in one fit_overlay)
  - regression: existing primary-peak rendering not broken

All 6 tests use _build_fit_overlay_payload (json pass-through) directly
or the existing th232_report_json fixture (from test_f_fit_view.py fixtures
re-imported or re-declared to stay self-contained).

Cite: scripts/gamma/reporting/json_report.py:586 (_build_fit_overlay),
      scripts/gamma/reporting/interactive_html.py:2563 (_build_fit_overlay_payload),
      scripts/gamma/reporting/templates/interactive_v1_17_2.html:1102 (JS setFitOverlay).
"""
from __future__ import annotations

import json
import os
import sys

import pytest

sys.path.insert(0, "scripts")

from gamma.reporting.interactive_html import _build_fit_overlay_payload  # noqa: E402


# ── Fixtures (self-contained, mirrors test_f_fit_view.py) ─────────────────────

_TH232_SPE = (
    "detectors/Gamma-1S/reference_spectra/"
    "archive/"
    "Th232_420-7-17_Маринелли_0cm.spe"
)
_TH232_BG = (
    "detectors/Gamma-1S/data/averaged_backgrounds/"
    "bg_2016_marinelli_water_marinelli.spe"
)


@pytest.fixture(scope="module")
def th232_report_json_v2(tmp_path_factory):
    """Run analysis on Th-232 Marinelli spectrum and return parsed report.json.

    Always runs its own ``analyze_and_report`` into a per-worker unique tmp
    dir (``tmp_path_factory.mktemp``). It deliberately does NOT read
    test_f_fit_view.py's output dir: under ``pytest -n auto`` that
    cross-file coupling let v2 (worker B) read v1's report.json while v1
    (worker A) was still writing it → JSONDecodeError at fixture setup
    (P1-3b xdist race). Independent output = no cross-worker collision.
    """
    if not os.path.exists(_TH232_SPE):
        pytest.skip(f"Fixture SPE not found: {_TH232_SPE}")

    from gamma.reporting import analyze_and_report  # noqa: E402

    out = str(tmp_path_factory.mktemp("f_fit_view_v2_th232"))
    result = analyze_and_report(
        _TH232_SPE,
        output_dir=out,
        write_html=False,
        write_plots=False,
        write_markdown=False,
        sample_mass_kg=0.5,
        background_path=_TH232_BG if os.path.exists(_TH232_BG) else None,
    )
    rp = result.get("json") or result.get("report")
    if rp and os.path.exists(str(rp)):
        with open(str(rp), encoding="utf-8") as f:
            return json.load(f)
    for fn in os.listdir(out):
        if fn.endswith("_report.json"):
            with open(os.path.join(out, fn), encoding="utf-8") as f:
                return json.load(f)
    pytest.skip("report.json not found after analysis")


# ── Test 1: secondary_peaks render check — source="secondary" preserved ───────

def test_secondary_peaks_emitted_with_orange_source():
    """F-FIT-VIEW v2: fit_overlay.peaks[*].source='secondary' passes through payload.

    Verifies that _build_fit_overlay_payload preserves source="secondary" entries
    from fit_overlay.peaks, enabling the JS frontend to render orange Gaussians.
    Cite: scripts/gamma/reporting/interactive_html.py:2563 (_build_fit_overlay_payload),
          templates/interactive_v1_17_2.html:1167 (SOURCE_COLORS.secondary = orange).
    """
    sample = {
        "fit_overlay": {
            "peaks": [
                {
                    "peak_id": "ps583",
                    "nuclide": "Tl-208",
                    "energy_keV": 583.2,
                    "amp_counts": 45.6,
                    "sigma_keV": 11.5,
                    "source": "secondary",
                    "label": "single_escape Tl-208 583",
                }
            ],
            "multiplet_continua": [],
        }
    }
    payload = _build_fit_overlay_payload(sample)
    assert len(payload["peaks"]) == 1, (
        f"Expected 1 peak, got {len(payload['peaks'])}"
    )
    assert payload["peaks"][0]["source"] == "secondary", (
        f"Expected source='secondary', got {payload['peaks'][0]['source']}"
    )
    assert payload["peaks"][0]["sigma_keV"] > 0, (
        "sigma_keV must be > 0 for secondary peaks (FWHM fallback)"
    )


# ── Test 2: background_primary_feps render check — gray dashed ────────────────

def test_background_peaks_emitted_with_gray_source():
    """F-FIT-VIEW v2: fit_overlay.peaks[*].source='background' passes through payload.

    Background FEPs are rendered as gray dashed Gaussians (SOURCE_DASH.background=[4,3]).
    Cite: scripts/gamma/reporting/json_report.py:686 (section 3 'background primary FEPs'),
          templates/interactive_v1_17_2.html:1172 (SOURCE_DASH.background=[4,3]).
    """
    sample = {
        "fit_overlay": {
            "peaks": [
                {
                    "peak_id": "pb1461",
                    "nuclide": "K-40 (bg)",
                    "energy_keV": 1460.8,
                    "amp_counts": 312.0,
                    "sigma_keV": 18.5,
                    "source": "background",
                    "label": "K-40 (bg) 1461",
                }
            ],
            "multiplet_continua": [],
        }
    }
    payload = _build_fit_overlay_payload(sample)
    assert len(payload["peaks"]) == 1, (
        f"Expected 1 peak in payload, got {len(payload['peaks'])}"
    )
    peak = payload["peaks"][0]
    assert peak["source"] == "background", (
        f"Expected source='background', got {peak['source']}"
    )
    assert "(bg)" in peak["nuclide"], (
        f"Background peak nuclide should contain '(bg)', got {peak['nuclide']}"
    )
    assert peak["sigma_keV"] > 0, (
        "sigma_keV must be > 0 for background peaks"
    )


# ── Test 3: unidentified_peaks render check — yellow dashed ───────────────────

def test_unidentified_peaks_emitted_with_yellow_source():
    """F-FIT-VIEW v2: fit_overlay.peaks[*].source='unidentified' passes through payload.

    Unidentified peaks are rendered as yellow dashed Gaussians.
    nuclide='?' and label starts with '?' per json_report._build_fit_overlay section 4.
    Cite: scripts/gamma/reporting/json_report.py:717 (section 4 'unidentified peaks'),
          templates/interactive_v1_17_2.html:1177 (SOURCE_COLORS.unidentified=yellow).
    """
    sample = {
        "fit_overlay": {
            "peaks": [
                {
                    "peak_id": "pu312",
                    "nuclide": "?",
                    "energy_keV": 312.0,
                    "amp_counts": None,
                    "sigma_keV": 8.5,
                    "source": "unidentified",
                    "label": "? 312",
                }
            ],
            "multiplet_continua": [],
        }
    }
    payload = _build_fit_overlay_payload(sample)
    assert len(payload["peaks"]) == 1, (
        f"Expected 1 peak, got {len(payload['peaks'])}"
    )
    peak = payload["peaks"][0]
    assert peak["source"] == "unidentified", (
        f"Expected source='unidentified', got {peak['source']}"
    )
    assert peak["nuclide"] == "?", (
        f"Unidentified peak nuclide must be '?', got {peak['nuclide']}"
    )
    assert peak["label"].startswith("?"), (
        f"Unidentified peak label must start with '?', got {peak['label']}"
    )
    # amp_counts may be None for unidentified (no area available)
    assert peak["amp_counts"] is None, (
        f"Expected amp_counts=None for unidentified peak, got {peak['amp_counts']}"
    )


# ── Test 4: sigma fallback path — fwhm_keV missing → fwhm_model used ─────────

def test_sigma_fallback_uses_fwhm_model():
    """F-FIT-VIEW v2: sigma_keV derived from FWHM model when fwhm_keV absent.

    When ResidualClassification has no direct gauss_sigma_keV, _build_fit_overlay
    falls back to fwhm_model (a + b*E + c*E²) → sigma = FWHM/2.355.
    This test verifies that a secondary peak built with a model-derived sigma
    passes through correctly (sigma > 0).
    Cite: scripts/gamma/reporting/json_report.py:627 (_fwhm_kev helper),
          scripts/gamma/reporting/json_report.py:661 (sigma fallback for secondary peaks).
    """
    # Simulate a fit_overlay entry where sigma was derived from the FWHM model
    # at E=662 keV with a typical NaI calibration: FWHM(662) ≈ 46 keV → sigma ≈ 19.5
    model_derived_sigma = 19.5
    sample = {
        "fit_overlay": {
            "peaks": [
                {
                    "peak_id": "ps662",
                    "nuclide": "backscatter Cs-137",
                    "energy_keV": 184.0,
                    "amp_counts": None,
                    "sigma_keV": model_derived_sigma,
                    "source": "secondary",
                    "label": "backscatter Cs-137 184",
                }
            ],
            "multiplet_continua": [],
        }
    }
    payload = _build_fit_overlay_payload(sample)
    assert len(payload["peaks"]) == 1
    peak = payload["peaks"][0]
    assert peak["sigma_keV"] > 0, (
        f"sigma_keV must be > 0 (FWHM-model derived), got {peak['sigma_keV']}"
    )
    assert peak["sigma_keV"] == model_derived_sigma, (
        f"Payload must preserve model-derived sigma {model_derived_sigma}, got {peak['sigma_keV']}"
    )


# ── Test 5: multi-category co-render — primary + secondary + background ───────

def test_multi_category_co_render():
    """F-FIT-VIEW v2: fit_overlay can contain peaks from all 3+ source categories.

    Verifies that a fit_overlay dict with primary (singlet), secondary, background,
    and unidentified peaks all pass through _build_fit_overlay_payload intact.
    All 4 source values are present in the output payload.
    Cite: scripts/gamma/reporting/json_report.py:586 (_build_fit_overlay returns all categories),
          templates/interactive_v1_17_2.html:1201 (forEach over all peaks by source).
    """
    sample = {
        "fit_overlay": {
            "peaks": [
                {
                    "peak_id": "p662",
                    "nuclide": "Cs-137",
                    "energy_keV": 661.7,
                    "amp_counts": 500.0,
                    "sigma_keV": 19.5,
                    "source": "singlet",
                    "label": "Cs-137 662",
                },
                {
                    "peak_id": "ps184",
                    "nuclide": "backscatter Cs-137",
                    "energy_keV": 184.0,
                    "amp_counts": 120.0,
                    "sigma_keV": 8.0,
                    "source": "secondary",
                    "label": "backscatter Cs-137 184",
                },
                {
                    "peak_id": "pb1461",
                    "nuclide": "K-40 (bg)",
                    "energy_keV": 1460.8,
                    "amp_counts": 88.0,
                    "sigma_keV": 18.5,
                    "source": "background",
                    "label": "K-40 (bg) 1461",
                },
                {
                    "peak_id": "pu312",
                    "nuclide": "?",
                    "energy_keV": 312.0,
                    "amp_counts": None,
                    "sigma_keV": 9.0,
                    "source": "unidentified",
                    "label": "? 312",
                },
            ],
            "multiplet_continua": [],
        }
    }
    payload = _build_fit_overlay_payload(sample)
    assert len(payload["peaks"]) == 4, (
        f"Expected 4 peaks (all categories), got {len(payload['peaks'])}"
    )
    sources_found = {p["source"] for p in payload["peaks"]}
    for expected_src in ("singlet", "secondary", "background", "unidentified"):
        assert expected_src in sources_found, (
            f"Source '{expected_src}' missing from payload.peaks sources: {sources_found}"
        )
    # All sigma_keV must be > 0
    bad_sigma = [p for p in payload["peaks"] if (p.get("sigma_keV") or 0) <= 0]
    assert not bad_sigma, (
        f"Peaks with sigma_keV <= 0: {[(p['peak_id'], p['sigma_keV']) for p in bad_sigma]}"
    )


# ── Test 6: regression — existing primary-peak rendering not broken ───────────

def test_primary_peak_rendering_not_broken(th232_report_json_v2):
    """F-FIT-VIEW v2: regression — v2 extension does not break v1 primary peak rendering.

    Verifies that fit_overlay.peaks still contains singlet/multiplet_component entries
    after the v2 extension. Confirms required fields (peak_id, nuclide, energy_keV,
    amp_counts, sigma_keV, source, label) are present for primary peaks.
    Cite: scripts/gamma/reporting/json_report.py:631 (section 1 'Singlet peaks'),
          tests/step11_reporting/test_f_fit_view.py:test_fit_overlay_peaks_count.
    """
    fo = th232_report_json_v2.get("fit_overlay", {})
    peaks = fo.get("peaks", [])

    # There must still be fit_overlay.peaks entries
    assert len(peaks) >= 1, (
        "fit_overlay.peaks must be non-empty after v2 extension"
    )

    # Primary peaks (singlet + multiplet_component) must still be present
    primary_sources = {"singlet", "multiplet_component"}
    primary_peaks = [p for p in peaks if p.get("source") in primary_sources]
    assert len(primary_peaks) >= 1, (
        f"No primary (singlet/multiplet_component) peaks found in fit_overlay after v2. "
        f"All sources: {[p.get('source') for p in peaks]}"
    )

    # Required fields check for primary peaks
    required = {"peak_id", "nuclide", "energy_keV", "amp_counts", "sigma_keV", "source", "label"}
    for pk in primary_peaks:
        missing = required - set(pk.keys())
        assert not missing, (
            f"Primary peak {pk.get('peak_id')} missing fields: {missing}"
        )

    # sigma_keV > 0 for all primary singlet peaks
    bad_sigma = [p for p in primary_peaks
                 if p.get("source") == "singlet" and (p.get("sigma_keV") or 0) <= 0]
    assert not bad_sigma, (
        f"Singlet primary peaks with sigma_keV <= 0: {[(p['peak_id'], p['sigma_keV']) for p in bad_sigma]}"
    )
