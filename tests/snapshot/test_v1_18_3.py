# -*- coding: utf-8 -*-
"""v1.18.3 — Peak-image production: .cpt CLI tool + legacy-bridge helper."""
from __future__ import annotations
import json, math, os, subprocess, sys
from pathlib import Path
import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


# ──────────────────────────────────────────────────────────────────
# F-299 bridge: legacy peak_image params → anchor
# ──────────────────────────────────────────────────────────────────

def test_F299_bridge_pure_gaussian_no_tail():
    """T=0, h_step=0 → tail_fraction=0, slope=0."""
    from gamma.peaks.peak_image_tabulated import (
        anchor_from_legacy_peak_image_params,
    )
    a = anchor_from_legacy_peak_image_params(
        E_keV=662.0, fwhm_keV=46.5,
        tail_param_T=0.0, h_step_frac=0.0,
    )
    assert a.E_keV == 662.0
    assert a.fwhm_keV == 46.5
    assert a.tail_fraction == 0.0
    assert a.tail_slope_inv_keV == 0.0
    assert a.step_height_frac == 0.0


def test_F299_bridge_typical_nai_params():
    """NaI typical: T=0.7, h_step=0.03."""
    from gamma.peaks.peak_image_tabulated import (
        anchor_from_legacy_peak_image_params,
    )
    a = anchor_from_legacy_peak_image_params(
        E_keV=662.0, fwhm_keV=46.5,
        tail_param_T=0.7, h_step_frac=0.03,
    )
    # tail_fraction = exp(-T²/2) = exp(-0.245) ≈ 0.783
    assert a.tail_fraction == pytest.approx(math.exp(-0.245), abs=0.01)
    # sigma = 46.5 / 2.3548 ≈ 19.75; slope = T/sigma = 0.7/19.75 ≈ 0.0354
    sigma = 46.5 / 2.354820045
    assert a.tail_slope_inv_keV == pytest.approx(0.7 / sigma, rel=1e-3)
    assert a.step_height_frac == 0.03


def test_F299_bridge_validates_inputs():
    from gamma.peaks.peak_image_tabulated import (
        anchor_from_legacy_peak_image_params,
    )
    with pytest.raises(ValueError):
        anchor_from_legacy_peak_image_params(
            E_keV=0, fwhm_keV=46.5,
        )
    with pytest.raises(ValueError):
        anchor_from_legacy_peak_image_params(
            E_keV=662.0, fwhm_keV=-1,
        )


def test_F299_bridge_negative_T_treated_as_pure_gauss():
    from gamma.peaks.peak_image_tabulated import (
        anchor_from_legacy_peak_image_params,
    )
    a = anchor_from_legacy_peak_image_params(
        E_keV=662.0, fwhm_keV=46.5,
        tail_param_T=-0.5,    # negative → ignored (pure gaussian)
    )
    assert a.tail_fraction == 0.0


# ──────────────────────────────────────────────────────────────────
# F-301 CLI tool: build/read/inspect roundtrip
# ──────────────────────────────────────────────────────────────────

CPT_TOOL = SCRIPTS / "cpt_tool.py"


def _run_cli(*args, timeout=20):
    """Run scripts/cpt_tool.py with given args."""
    cmd = [sys.executable, str(CPT_TOOL), *map(str, args)]
    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"
    return subprocess.run(
        cmd, capture_output=True, text=True, timeout=timeout,
        env=env, encoding="utf-8", errors="replace",
    )


def test_cli_build_creates_valid_cpt(tmp_path):
    out = tmp_path / "test.cpt"
    r = _run_cli(
        "build",
        "--detector-id", "Gamma-1S",
        "--detector-class", "NaI",
        "--diameter-mm", "63",
        "--anchor", "122.0:12.5",
        "--anchor", "662.0:46.5",
        "--anchor", "1332.0:78.0",
        "--out", str(out),
    )
    assert r.returncode == 0, r.stderr
    assert out.exists()
    assert out.stat().st_size > 100
    xml = out.read_text(encoding="utf-8")
    assert "Gamma-1S" in xml
    assert "<anchor" in xml
    assert "662" in xml


def test_cli_read_roundtrip_via_build(tmp_path):
    out = tmp_path / "test.cpt"
    _run_cli(
        "build",
        "--detector-id", "Gamma-1S",
        "--detector-class", "NaI",
        "--diameter-mm", "63",
        "--anchor", "662.0:46.5",
        "--out", str(out),
    )
    r = _run_cli("read", str(out))
    assert r.returncode == 0, r.stderr
    payload = json.loads(r.stdout)
    assert payload["detector_id"] == "Gamma-1S"
    assert payload["detector_class"] == "NaI"
    assert payload["crystal_diameter_mm"] == pytest.approx(63.0)
    assert len(payload["anchors"]) == 1
    assert payload["anchors"][0]["E_keV"] == pytest.approx(662.0)


def test_cli_read_with_output_file(tmp_path):
    cpt = tmp_path / "test.cpt"
    js = tmp_path / "out.json"
    _run_cli(
        "build", "--detector-id", "X", "--detector-class", "NaI",
        "--diameter-mm", "63", "--anchor", "662.0:46.5",
        "--out", str(cpt),
    )
    r = _run_cli("read", str(cpt), "-o", str(js))
    assert r.returncode == 0
    assert js.exists()
    payload = json.loads(js.read_text(encoding="utf-8"))
    assert payload["detector_id"] == "X"


def test_cli_inspect_prints_summary(tmp_path):
    cpt = tmp_path / "test.cpt"
    _run_cli(
        "build", "--detector-id", "Gamma-1S", "--detector-class", "NaI",
        "--diameter-mm", "63",
        "--anchor", "122.0:12.5",
        "--anchor", "662.0:46.5",
        "--out", str(cpt),
    )
    r = _run_cli("inspect", str(cpt))
    assert r.returncode == 0, r.stderr
    assert "Gamma-1S" in r.stdout
    assert "NaI" in r.stdout
    assert "anchors" in r.stdout.lower()
    assert "FWHM%@662keV" in r.stdout


def test_cli_read_nonexistent_returns_error(tmp_path):
    r = _run_cli("read", str(tmp_path / "no-such.cpt"))
    assert r.returncode != 0
    assert "ERROR" in r.stderr


def test_cli_build_invalid_anchor_format_rejected(tmp_path):
    out = tmp_path / "test.cpt"
    r = _run_cli(
        "build", "--detector-id", "X", "--detector-class", "NaI",
        "--diameter-mm", "63",
        "--anchor", "not_a_pair_at_all",
        "--out", str(out),
    )
    assert r.returncode != 0


def test_cli_no_command_errors():
    r = _run_cli()
    assert r.returncode != 0
