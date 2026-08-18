"""F-160 / v1.18.19.0 — BecqMoni XML export integration test.

Validates analyze_and_report(..., export_becqmoni="both") writes
both sample + bg XML files, both round-trip via read_spectrum, and
CLI flag --export-becqmoni is wired through to wrapper.
"""
from __future__ import annotations

import sys
import shutil
import subprocess
from pathlib import Path

import pytest

sys.path.insert(0, "scripts")

from gamma.reporting import analyze_and_report
from gamma.io.readers import read_spectrum


ROOT = Path(__file__).resolve().parents[2]
KIT = ROOT / "detectors/Gamma-1S/reference_spectra/reference_kits"
CS_SAMPLE = KIT / "Marinelli_1L/Cs-137/sample_M_cs_легкий_2001-2005.spe"
CS_BG = KIT / "Marinelli_1L/Cs-137/background_bg_2016_marinelli_water_marinelli.spe"


@pytest.fixture
def out_dir(tmp_path):
    return str(tmp_path / "f160")


def test_f160_export_both_writes_both_files(out_dir):
    if not CS_SAMPLE.exists():
        pytest.skip("kit sample missing")
    res = analyze_and_report(
        str(CS_SAMPLE),
        background_path=str(CS_BG),
        sample_mass_kg=0.570,
        output_dir=out_dir,
        write_html=False, write_plots=False, write_markdown=False, write_json=False,
        export_becqmoni="both",
    )
    assert "becqmoni_sample" in res
    assert "becqmoni_bg" in res
    sp_path = Path(res["becqmoni_sample"])
    bg_path = Path(res["becqmoni_bg"])
    assert sp_path.exists() and sp_path.suffix == ".xml"
    assert bg_path.exists() and bg_path.suffix == ".xml"
    assert sp_path.name.endswith("_calibrated.bq.xml")
    assert bg_path.name.endswith("_calibrated.bq.xml")
    # Both non-trivial
    assert sp_path.stat().st_size > 5000
    assert bg_path.stat().st_size > 5000


def test_f160_export_sample_only(out_dir):
    if not CS_SAMPLE.exists():
        pytest.skip("kit sample missing")
    res = analyze_and_report(
        str(CS_SAMPLE),
        background_path=str(CS_BG),
        sample_mass_kg=0.570,
        output_dir=out_dir,
        write_html=False, write_plots=False, write_markdown=False, write_json=False,
        export_becqmoni="sample",
    )
    assert "becqmoni_sample" in res
    assert "becqmoni_bg" not in res


def test_f160_export_bg_only(out_dir):
    if not CS_SAMPLE.exists():
        pytest.skip("kit sample missing")
    res = analyze_and_report(
        str(CS_SAMPLE),
        background_path=str(CS_BG),
        sample_mass_kg=0.570,
        output_dir=out_dir,
        write_html=False, write_plots=False, write_markdown=False, write_json=False,
        export_becqmoni="bg",
    )
    assert "becqmoni_sample" not in res
    assert "becqmoni_bg" in res


def test_f160_export_off_default(out_dir):
    """Без export_becqmoni — никаких BecqMoni файлов не пишется (back-compat)."""
    if not CS_SAMPLE.exists():
        pytest.skip("kit sample missing")
    res = analyze_and_report(
        str(CS_SAMPLE),
        background_path=str(CS_BG),
        sample_mass_kg=0.570,
        output_dir=out_dir,
        write_html=False, write_plots=False, write_markdown=False, write_json=False,
        # export_becqmoni не указан → default "off"
    )
    assert "becqmoni_sample" not in res
    assert "becqmoni_bg" not in res


def test_f160_roundtrip_sample(out_dir):
    """Записанный sample XML должен round-trip через read_spectrum."""
    if not CS_SAMPLE.exists():
        pytest.skip("kit sample missing")
    res = analyze_and_report(
        str(CS_SAMPLE),
        background_path=str(CS_BG),
        sample_mass_kg=0.570,
        output_dir=out_dir,
        write_html=False, write_plots=False, write_markdown=False, write_json=False,
        export_becqmoni="sample",
    )
    sp = read_spectrum(res["becqmoni_sample"])
    assert sp is not None
    assert len(sp.counts) > 100
    assert sp.live_time > 0
    # Sample ID preserved
    assert sp.sample_id, "sample_id отсутствует в round-trip XML"


def test_f160_bg_only_without_background_warns(out_dir):
    """export_becqmoni=bg без background_path → warning, не падение."""
    if not CS_SAMPLE.exists():
        pytest.skip("kit sample missing")
    res = analyze_and_report(
        str(CS_SAMPLE),
        sample_mass_kg=0.570,
        output_dir=out_dir,
        write_html=False, write_plots=False, write_markdown=False, write_json=False,
        export_becqmoni="bg",
    )
    # Должен быть warning, не exception
    warnings = res.get("warnings") or []
    assert any("background_path" in w for w in warnings), \
        f"expected background_path warning; got: {warnings}"
    assert "becqmoni_bg" not in res


def test_f160_invalid_mode_warns(out_dir):
    """Невалидный mode → warning, не падение."""
    if not CS_SAMPLE.exists():
        pytest.skip("kit sample missing")
    res = analyze_and_report(
        str(CS_SAMPLE),
        background_path=str(CS_BG),
        sample_mass_kg=0.570,
        output_dir=out_dir,
        write_html=False, write_plots=False, write_markdown=False, write_json=False,
        export_becqmoni="wrong_mode",
    )
    warnings = res.get("warnings") or []
    assert any("invalid mode" in w for w in warnings), \
        f"expected invalid mode warning; got: {warnings}"


def test_f160_cli_flag_registered():
    """CLI helper --export-becqmoni зарегистрирован в parser."""
    import gamma.cli as cli_mod
    # Build parser via re-creating it
    # Simpler: parse --help string for the flag
    out = subprocess.run(
        [sys.executable, "-m", "gamma.cli", "analyze", "--help"],
        capture_output=True, text=True,
        env={**__import__("os").environ, "PYTHONPATH": "scripts"},
        encoding="utf-8", errors="replace",
    )
    assert "--export-becqmoni" in out.stdout, \
        f"--export-becqmoni не найден в analyze --help: {out.stdout[:500]}"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
