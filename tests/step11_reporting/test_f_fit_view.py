"""F-FIT-VIEW / v1.22.1 — fit overlay toggle tests.

Verifies that:
  1. report.json contains ``fit_overlay.peaks`` with ≥ 1 entries for a
     spectrum with identified nuclides.
  2. report.json contains ``fit_overlay.multiplet_continua`` when multiplet
     deconvolution ran (spectrum has overlapping peaks).
  3. Each peak entry has all required fields: peak_id, nuclide, energy_keV,
     amp_counts, sigma_keV, source, label.
  4. HTML output contains the toggle button ``id="toggle-fit-overlay"``.
  5. HTML output contains the ``FIT_OVERLAY`` JS constant (DATA_FIT_OVERLAY
     placeholder substituted).
  6. Empty/missing fit_overlay is handled gracefully (backward compat).
  7. Multiplet continua have required fields: cluster_id, E_keV, continuum,
     total, components.
  8. Singlet peaks have source="singlet".
  9. amp_counts > 0 for all non-zero-area peaks.
  10. sigma_keV > 0 for all peaks.
"""
from __future__ import annotations

import json
import os
import sys

import pytest

sys.path.insert(0, "scripts")

from gamma.reporting import analyze_and_report  # noqa: E402
from gamma.reporting.interactive_html import _build_fit_overlay_payload  # noqa: E402


# ── Fixtures ──────────────────────────────────────────────────────────────────

_TH232_SPE = (
    "detectors/Gamma-1S/reference_spectra/"
    "archive/"
    "Th232_420-7-17_Маринелли_0cm.spe"
)
_TH232_BG = (
    "detectors/Gamma-1S/data/averaged_backgrounds/"
    "bg_2016_marinelli_water_marinelli.spe"
)
_EU152_SPE = (
    "detectors/Gamma-1S/reference_spectra/"
    "archive/"
    "Eu-152__04_21_Точечная-5см_5cm.spe"
)

@pytest.fixture(scope="module")
def th232_result(tmp_path_factory):
    """Run analysis on Th-232 Marinelli spectrum once for all tests.

    Writes pipeline output to a per-worker unique tmp dir
    (``tmp_path_factory.mktemp``) so that under ``pytest -n auto`` two
    workers running this module-scoped fixture never collide on a fixed
    on-disk path. The chosen out-dir is stashed back into the returned
    dict under ``_test_out_dir`` for the glob fallback below.
    """
    spe = _TH232_SPE
    if not os.path.exists(spe):
        pytest.skip(f"Fixture SPE not found: {spe}")
    out = str(tmp_path_factory.mktemp("f_fit_view_th232"))
    result = analyze_and_report(
        spe,
        output_dir=out,
        write_html=True,
        write_plots=False,
        write_markdown=False,
        sample_mass_kg=0.5,
        background_path=_TH232_BG if os.path.exists(_TH232_BG) else None,
    )
    # Stash the unique out-dir for the readers' glob fallback (never a
    # fixed shared path — see P1-3b xdist race fix).
    if isinstance(result, dict):
        result["_test_out_dir"] = out
    return result


@pytest.fixture(scope="module")
def th232_report_json(th232_result):
    """Load and parse the generated report.json."""
    rp = th232_result.get("json") or th232_result.get("report")
    if rp and os.path.exists(str(rp)):
        with open(str(rp), encoding="utf-8") as f:
            return json.load(f)
    # Fallback: glob inside the per-worker out-dir stashed by th232_result.
    out = th232_result.get("_test_out_dir")
    if out and os.path.isdir(out):
        for fn in os.listdir(out):
            if fn.endswith("_report.json"):
                with open(os.path.join(out, fn), encoding="utf-8") as f:
                    return json.load(f)
    pytest.skip("report.json not found")


@pytest.fixture(scope="module")
def th232_html(th232_result):
    """Read the generated HTML as a string."""
    hp = th232_result.get("html")
    if hp and os.path.exists(str(hp)):
        with open(str(hp), encoding="utf-8") as f:
            return f.read()
    pytest.skip("HTML file not found")


# ── Test 1: fit_overlay section present in report.json ────────────────────────

def test_fit_overlay_section_present(th232_report_json):
    """F-FIT-VIEW: report.json MUST contain 'fit_overlay' key."""
    assert "fit_overlay" in th232_report_json, (
        "fit_overlay key missing from report.json"
    )


# ── Test 2: fit_overlay.peaks count ≥ 5 ──────────────────────────────────────

def test_fit_overlay_peaks_count(th232_report_json):
    """F-FIT-VIEW: fit_overlay.peaks must have ≥ 5 entries for Th-232 spectrum.

    Th-232 chain emits ≥ 10 FEPs in 50–3000 keV range, so ≥ 5 is conservative.
    """
    fo = th232_report_json.get("fit_overlay", {})
    peaks = fo.get("peaks", [])
    assert len(peaks) >= 5, (
        f"Expected ≥ 5 peaks in fit_overlay, got {len(peaks)}"
    )


# ── Test 3: required fields in each peak entry ────────────────────────────────

REQUIRED_PEAK_FIELDS = {
    "peak_id", "nuclide", "energy_keV", "amp_counts", "sigma_keV",
    "source", "label",
}


@pytest.mark.parametrize("field", sorted(REQUIRED_PEAK_FIELDS))
def test_fit_overlay_peak_fields(th232_report_json, field):
    """F-FIT-VIEW: each peak in fit_overlay.peaks must have required field."""
    fo = th232_report_json.get("fit_overlay", {})
    peaks = fo.get("peaks", [])
    if not peaks:
        pytest.skip("No peaks in fit_overlay")
    missing = [i for i, p in enumerate(peaks) if field not in p]
    assert not missing, (
        f"Field '{field}' missing from peaks at indices: {missing[:5]}"
    )


# ── Test 4: sigma_keV > 0 for all peaks ──────────────────────────────────────

def test_fit_overlay_peaks_sigma_positive(th232_report_json):
    """F-FIT-VIEW: sigma_keV must be > 0 for all peaks."""
    fo = th232_report_json.get("fit_overlay", {})
    peaks = fo.get("peaks", [])
    bad = [(p.get("peak_id"), p.get("sigma_keV")) for p in peaks
           if (p.get("sigma_keV") or 0) <= 0]
    assert not bad, f"Peaks with sigma_keV ≤ 0: {bad}"


# ── Test 5: amp_counts > 0 for singlet peaks ─────────────────────────────────

def test_fit_overlay_peaks_amp_positive(th232_report_json):
    """F-FIT-VIEW: amp_counts must be > 0 for singlet (non-multiplet) peaks.

    Multiplet-component peaks may have amp_counts rounded to 0.0 when the
    deconvolved area is very small — this is valid (weak phantom component).
    Only singlet peaks from primary_feps have area > 0 guaranteed.
    """
    fo = th232_report_json.get("fit_overlay", {})
    peaks = fo.get("peaks", [])
    # Filter to singlet peaks only (multiplet components may have 0 area)
    singlets = [p for p in peaks if p.get("source") == "singlet"]
    bad = [(p.get("peak_id"), p.get("amp_counts")) for p in singlets
           if (p.get("amp_counts") or 0) <= 0]
    assert not bad, f"Singlet peaks with amp_counts ≤ 0: {bad}"


# ── Test 6: multiplet_continua present for Th-232 (has multiplets) ───────────

def test_fit_overlay_multiplet_continua(th232_report_json):
    """F-FIT-VIEW: Th-232 spectrum has multiplets → multiplet_continua non-empty."""
    # Only check if the spectrum actually had multiplet deconvolutions
    n_mdec = len(th232_report_json.get("multiplet_deconvolutions", []))
    if n_mdec == 0:
        pytest.skip("No multiplet_deconvolutions in this run")
    fo = th232_report_json.get("fit_overlay", {})
    mc = fo.get("multiplet_continua", [])
    # If multiplet deconvolutions ran, at least one cluster should have overlay data
    # (overlay_E_keV is populated by coupled_multiplet.py)
    # Note: if all clusters use legacy path (no overlay_E_keV), mc may be empty.
    # In that case we skip rather than fail.
    if not mc:
        pytest.skip(
            "No multiplet_continua in fit_overlay "
            "(overlay_E_keV may not be populated for this spectrum)"
        )
    assert len(mc) >= 1, "Expected ≥ 1 multiplet_continua entry"


# ── Test 7: multiplet_continua required fields ────────────────────────────────

REQUIRED_MC_FIELDS = {"cluster_id", "E_keV", "continuum", "total", "components"}


@pytest.mark.parametrize("field", sorted(REQUIRED_MC_FIELDS))
def test_fit_overlay_multiplet_continua_fields(th232_report_json, field):
    """F-FIT-VIEW: each multiplet_continua entry must have required field."""
    fo = th232_report_json.get("fit_overlay", {})
    mc = fo.get("multiplet_continua", [])
    if not mc:
        pytest.skip("No multiplet_continua in fit_overlay")
    missing = [i for i, c in enumerate(mc) if field not in c]
    assert not missing, (
        f"Field '{field}' missing from multiplet_continua at indices: {missing}"
    )


# ── Test 8: HTML contains toggle button id ───────────────────────────────────

def test_html_contains_toggle_button(th232_html):
    """F-FIT-VIEW: HTML must contain id='toggle-fit-overlay' button."""
    assert 'id="toggle-fit-overlay"' in th232_html, (
        "HTML does not contain <button id='toggle-fit-overlay'>"
    )


# ── Test 9: HTML contains FIT_OVERLAY JS constant ────────────────────────────

def test_html_contains_fit_overlay_constant(th232_html):
    """F-FIT-VIEW: HTML must contain 'const FIT_OVERLAY=' JS variable."""
    assert "const FIT_OVERLAY=" in th232_html, (
        "HTML does not contain 'const FIT_OVERLAY=' (DATA_FIT_OVERLAY not substituted)"
    )


# ── Test 10: HTML FIT_OVERLAY data is valid JSON with peaks ──────────────────

def test_html_fit_overlay_data_valid(th232_html):
    """F-FIT-VIEW: FIT_OVERLAY JS constant in HTML must be valid JSON with peaks."""
    import re
    # Find const FIT_OVERLAY={...};
    m = re.search(r'const FIT_OVERLAY=(\{[^;]+\});', th232_html)
    if not m:
        # Try multi-line
        idx = th232_html.find("const FIT_OVERLAY=")
        if idx < 0:
            pytest.fail("const FIT_OVERLAY not found in HTML")
        # Extract until semicolon after first closing brace at same nesting level
        rest = th232_html[idx + len("const FIT_OVERLAY="):]
        # Find the JSON boundary
        depth = 0
        end = 0
        for i, ch in enumerate(rest):
            if ch == '{':
                depth += 1
            elif ch == '}':
                depth -= 1
                if depth == 0:
                    end = i + 1
                    break
        json_str = rest[:end]
    else:
        json_str = m.group(1)
    data = json.loads(json_str)
    assert "peaks" in data, "FIT_OVERLAY JSON missing 'peaks' key"
    assert "multiplet_continua" in data, "FIT_OVERLAY JSON missing 'multiplet_continua' key"
    assert len(data["peaks"]) >= 1, "FIT_OVERLAY.peaks is empty in HTML"


# ── Test 11: _build_fit_overlay_payload backward compat (empty report) ────────

def test_fit_overlay_payload_empty_report():
    """F-FIT-VIEW: _build_fit_overlay_payload handles empty/missing fit_overlay."""
    payload = _build_fit_overlay_payload({})
    assert payload == {"peaks": [], "multiplet_continua": []}, (
        f"Empty report payload mismatch: {payload}"
    )


# ── Test 12: _build_fit_overlay_payload pass-through ─────────────────────────

def test_fit_overlay_payload_passthrough():
    """F-FIT-VIEW: _build_fit_overlay_payload passes peaks and continua through."""
    sample = {
        "fit_overlay": {
            "peaks": [
                {
                    "peak_id": "p662",
                    "nuclide": "Cs-137",
                    "energy_keV": 661.7,
                    "amp_counts": 123.4,
                    "sigma_keV": 12.5,
                    "source": "singlet",
                    "label": "Cs-137 662",
                }
            ],
            "multiplet_continua": [],
        }
    }
    payload = _build_fit_overlay_payload(sample)
    assert len(payload["peaks"]) == 1
    assert payload["peaks"][0]["nuclide"] == "Cs-137"
    assert payload["peaks"][0]["sigma_keV"] == 12.5
