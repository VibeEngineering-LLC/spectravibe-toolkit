"""F-130 / v1.17.7 — auto-detect sample density из .spe метаданных."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))


def test_auto_extract_density_from_material_ro():
    """MATERIAL.Ro = 1.6 → ρ_sample=1.6, source='material_ro'."""
    from gamma.io.lsrm_spe import _auto_extract_density
    fields = {"MATERIAL": '{"Name":"OISN-16","Ro":1.6,"Compound":[]}'}
    result = _auto_extract_density(fields)
    assert result == (1.6, "material_ro")


def test_auto_extract_density_from_mass_volume():
    """SAMPLEMASS=1600;16, SAMPLEVOLUME=1000;10 → ρ=1.6."""
    from gamma.io.lsrm_spe import _auto_extract_density
    fields = {
        "SAMPLEMASS": "1600.0;16.0",
        "SAMPLEVOLUME": "1000.0;10.0",
    }
    result = _auto_extract_density(fields)
    assert result is not None
    assert abs(result[0] - 1.6) < 1e-6
    assert result[1] == "sample_mass_over_volume"


def test_auto_extract_density_probe_fallback():
    """Если SAMPLE* нет, используется PROBE*."""
    from gamma.io.lsrm_spe import _auto_extract_density
    fields = {
        "PROBEMASS": "500.0;5.0",
        "PROBEVOLUME": "250.0;2.5",
    }
    result = _auto_extract_density(fields)
    assert result is not None
    assert abs(result[0] - 2.0) < 1e-6
    assert result[1] == "probe_mass_over_volume"


def test_auto_extract_density_returns_none_when_missing():
    """Без полей → None."""
    from gamma.io.lsrm_spe import _auto_extract_density
    assert _auto_extract_density({}) is None


def test_auto_extract_density_rejects_pathological_values():
    """ρ < 0.1 или > 10 → отвергается (санитарный диапазон)."""
    from gamma.io.lsrm_spe import _auto_extract_density
    fields = {"MATERIAL": '{"Ro":100.0}'}
    assert _auto_extract_density(fields) is None
    fields = {"MATERIAL": '{"Ro":0.01}'}
    assert _auto_extract_density(fields) is None


def test_th232_fixture_auto_density():
    """E2E: Th-232 fixture сам подбирает ρ=1.6 в spec.extras."""
    from gamma.io.readers import read_spectrum
    fixture = (Path(__file__).parent.parent.parent
               / "detectors" / "Gamma-1S" / "reference_spectra"
               / "archive"
               / "Th232_420-7-17_Маринелли_0cm.spe")
    if not fixture.exists():
        pytest.skip(f"fixture missing: {fixture}")
    spec = read_spectrum(str(fixture))
    assert spec.extras.get("lsrm_sample_density_g_cm3") == 1.6
    assert spec.extras.get("lsrm_density_source") == "material_ro"


def test_staged_pipeline_uses_auto_density():
    """E2E: analyze_lsrm_spe без --sample-density-g-cm3 должен авто-найти ρ."""
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    fixture = (Path(__file__).parent.parent.parent
               / "detectors" / "Gamma-1S" / "reference_spectra"
               / "archive"
               / "Th232_420-7-17_Маринелли_0cm.spe")
    if not fixture.exists():
        pytest.skip(f"fixture missing: {fixture}")
    r = analyze_lsrm_spe(str(fixture), complete_workflow=True)
    # F-130 нарративная заметка должна упоминать material_ro
    notes_str = " ".join(r.notes or [])
    assert "F-122" in notes_str, (
        "F-122 self-attenuation note должна быть в pipeline notes "
        "(triggered by auto-detected density)"
    )


def test_explicit_density_overrides_auto():
    """Если CLI задал ρ, авто-detection пропускается."""
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    fixture = (Path(__file__).parent.parent.parent
               / "detectors" / "Gamma-1S" / "reference_spectra"
               / "archive"
               / "Th232_420-7-17_Маринелли_0cm.spe")
    if not fixture.exists():
        pytest.skip(f"fixture missing: {fixture}")
    r = analyze_lsrm_spe(
        str(fixture), complete_workflow=True,
        sample_density_g_cm3=1.95,
    )
    notes_str = " ".join(r.notes or [])
    # CLI override должен показываться как "CLI флаг" в нарративе
    assert "1.95" in notes_str or "CLI" in notes_str
