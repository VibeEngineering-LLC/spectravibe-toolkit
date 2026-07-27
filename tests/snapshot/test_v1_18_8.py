# -*- coding: utf-8 -*-
"""v1.18.8 — CLI flags для opt-in activity integrations v1.18.1..v1.18.4.

Verifies:
- argparse принимает все 7 новых флагов без ошибок
- флаги корректно пробрасываются через CLI → wrapper → orchestrator
- --help содержит описание каждого флага
- back-compat: без флагов pipeline работает как раньше (default OFF)
"""
from __future__ import annotations
import os, subprocess, sys
from pathlib import Path
import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


def _run_cli(*args, timeout=30):
    """Run python -m gamma.cli with given args."""
    cmd = [sys.executable, "-m", "gamma.cli", *map(str, args)]
    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONPATH"] = str(SCRIPTS) + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        cmd, capture_output=True, text=True, timeout=timeout,
        env=env, encoding="utf-8", errors="replace",
        cwd=str(REPO),
    )


# ──────────────────────────────────────────────────────────────────
# CLI argparse acceptance
# ──────────────────────────────────────────────────────────────────

def test_cli_help_lists_v1_18_8_flags():
    r = _run_cli("analyze", "--help")
    assert r.returncode == 0, r.stderr
    for flag in [
        "--enable-tcs-correction",
        "--tcs-detector-id",
        "--enable-cutshall-self-abs",
        "--cutshall-path-cm",
        "--cutshall-calib-density-g-cm3",
        "--enable-matrix-method",
        "--matrix-method-energy-tolerance-keV",
    ]:
        assert flag in r.stdout, f"Missing flag in --help: {flag}"


def test_cli_help_mentions_f_ids_for_each_flag():
    r = _run_cli("analyze", "--help")
    assert r.returncode == 0
    for f_id in ["F-294", "F-295", "F-296", "F-297"]:
        assert f_id in r.stdout, f"Missing F-id reference: {f_id}"


def test_cli_help_mentions_default_off():
    """User должен знать что флаги default OFF (back-compat)."""
    r = _run_cli("analyze", "--help")
    assert "Default OFF" in r.stdout


def test_cli_rejects_invalid_tcs_detector_id():
    """tcs-detector-id ограничен choices."""
    r = _run_cli(
        "analyze", "fake.spe",
        "--tcs-detector-id", "NOT_A_VALID_DETECTOR",
    )
    assert r.returncode != 0
    assert "invalid choice" in r.stderr.lower() or "error" in r.stderr.lower()


def test_cli_accepts_valid_tcs_detector_id():
    """Все 4 валидных choices принимаются."""
    for det in ["Gamma-1S", "3in3", "4in4", "NaI_63x63"]:
        r = _run_cli(
            "analyze", "/non/existent/file.spe",
            "--tcs-detector-id", det,
        )
        # File doesn't exist → ERROR на read stage, но argparse ОК
        # (если argparse падает, мы получили бы "invalid choice")
        assert "invalid choice" not in r.stderr.lower(), (
            f"Detector {det} unexpectedly rejected"
        )


# ──────────────────────────────────────────────────────────────────
# Wrapper passthrough
# ──────────────────────────────────────────────────────────────────

def test_wrapper_orchestrator_keys_include_v1_18_8():
    """_ORCHESTRATOR_KEYS должен содержать все 7 новых ключей."""
    from gamma.reporting.wrapper import _ORCHESTRATOR_KEYS
    for key in [
        "enable_tcs_correction",
        "tcs_detector_id",
        "enable_cutshall_self_abs",
        "cutshall_path_cm",
        "cutshall_calib_density_g_cm3",
        "enable_matrix_method",
        "matrix_method_energy_tolerance_keV",
    ]:
        assert key in _ORCHESTRATOR_KEYS, (
            f"Missing key in _ORCHESTRATOR_KEYS: {key}"
        )


# ──────────────────────────────────────────────────────────────────
# Orchestrator signature acceptance
# ──────────────────────────────────────────────────────────────────

def test_orchestrator_accepts_v1_18_8_kwargs():
    """analyze_lsrm_spe сигнатура должна принимать все 7 kwargs."""
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    import inspect
    sig = inspect.signature(analyze_lsrm_spe)
    for key in [
        "enable_tcs_correction",
        "tcs_detector_id",
        "enable_cutshall_self_abs",
        "cutshall_path_cm",
        "cutshall_calib_density_g_cm3",
        "enable_matrix_method",
        "matrix_method_energy_tolerance_keV",
    ]:
        assert key in sig.parameters, (
            f"analyze_lsrm_spe missing parameter: {key}"
        )


def test_orchestrator_defaults_all_off():
    """All enable-flags default to False, identifiers/numbers — sane defaults."""
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    import inspect
    params = inspect.signature(analyze_lsrm_spe).parameters
    assert params["enable_tcs_correction"].default is False
    assert params["enable_cutshall_self_abs"].default is False
    assert params["enable_matrix_method"].default is False
    assert params["tcs_detector_id"].default == "Gamma-1S"
    assert params["cutshall_calib_density_g_cm3"].default == 1.0
    assert params["matrix_method_energy_tolerance_keV"].default == 1.0


# ──────────────────────────────────────────────────────────────────
# Compute layer signature (already tested in test_v1_18_1.py + test_v1_18_2.py,
# но повторяем здесь как safety-net)
# ──────────────────────────────────────────────────────────────────

def test_compute_activities_for_all_accepts_v1_18_8_kwargs():
    from gamma.activity.compute import compute_activities_for_all
    import inspect
    sig = inspect.signature(compute_activities_for_all)
    for key in [
        "enable_tcs_correction",
        "tcs_detector_id",
        "enable_cutshall_self_abs",
        "cutshall_path_cm",
        "cutshall_calib_density_g_cm3",
        "enable_matrix_method",
        "matrix_method_energy_tolerance_keV",
    ]:
        assert key in sig.parameters
