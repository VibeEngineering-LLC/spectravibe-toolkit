"""F-129 / v1.17.7 — CLI флаг --peak-search-method для convolution / compare."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))


def test_staged_pipeline_accepts_peak_search_method():
    """analyze_lsrm_spe должен принимать peak_search_method kwarg."""
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    fixture = (Path(__file__).parent.parent.parent
               / "detectors" / "Gamma-1S" / "reference_spectra"
               / "archive"
               / "Cs137_420-7-14_Маринелли_0cm.spe")
    if not fixture.exists():
        pytest.skip(f"fixture missing: {fixture}")
    r = analyze_lsrm_spe(str(fixture), complete_workflow=True,
                         peak_search_method="convolution")
    assert r.peak_search_method == "convolution"


def test_compare_mode_emits_comparison_block():
    """Режим 'compare' должен заполнить peak_search_method_comparison."""
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    fixture = (Path(__file__).parent.parent.parent
               / "detectors" / "Gamma-1S" / "reference_spectra"
               / "archive"
               / "Cs137_420-7-14_Маринелли_0cm.spe")
    if not fixture.exists():
        pytest.skip(f"fixture missing: {fixture}")
    r = analyze_lsrm_spe(str(fixture), complete_workflow=True,
                         peak_search_method="compare")
    assert r.peak_search_method_comparison is not None
    cmp = r.peak_search_method_comparison
    assert "agreement_fraction" in cmp or "n_a" in cmp


def test_mariscotti_remains_default():
    """Без kwarg метод = 'mariscotti'."""
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    fixture = (Path(__file__).parent.parent.parent
               / "detectors" / "Gamma-1S" / "reference_spectra"
               / "archive"
               / "Cs137_420-7-14_Маринелли_0cm.spe")
    if not fixture.exists():
        pytest.skip(f"fixture missing: {fixture}")
    r = analyze_lsrm_spe(str(fixture), complete_workflow=True)
    assert r.peak_search_method == "mariscotti"
    assert r.peak_search_method_comparison is None


def test_unknown_method_falls_back():
    """Неизвестный метод → fallback к mariscotti, без падения."""
    from gamma.identification.staged_pipeline import _run_peak_search
    from gamma.io.readers import read_spectrum
    fixture = (Path(__file__).parent.parent.parent
               / "detectors" / "Gamma-1S" / "reference_spectra"
               / "archive"
               / "Cs137_420-7-14_Маринелли_0cm.spe")
    if not fixture.exists():
        pytest.skip(f"fixture missing: {fixture}")
    spec = read_spectrum(str(fixture))
    peaks, cmp = _run_peak_search(
        spec, lambda ch: 30.0, 3.0, method="bogus_method",
    )
    assert peaks is not None
    assert cmp is None


def test_json_report_includes_peak_search_method():
    """json_report должен сериализовать peak_search_method в diagnostics."""
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    from gamma.reporting.json_report import build_json_report
    fixture = (Path(__file__).parent.parent.parent
               / "detectors" / "Gamma-1S" / "reference_spectra"
               / "archive"
               / "Cs137_420-7-14_Маринелли_0cm.spe")
    if not fixture.exists():
        pytest.skip(f"fixture missing: {fixture}")
    r = analyze_lsrm_spe(str(fixture), complete_workflow=True,
                         peak_search_method="compare")
    rep = build_json_report(r)
    diag = rep.get("diagnostics", {})
    assert diag.get("peak_search_method") == "compare"
    assert diag.get("peak_search_method_comparison") is not None
