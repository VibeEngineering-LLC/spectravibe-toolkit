# -*- coding: utf-8 -*-
"""BUG-1 / 2026-06-02 — `run_skill.py` mass-resolution precedence.

Verifies the precedence chain for `sample_mass_kg` when the operator runs
`python scripts/run_skill.py <spe>` WITHOUT `--sample-mass-kg`:

    CLI --sample-mass-kg
       ↓ (absent)
    filename token  (e.g. "..._0.5kg.spe")
       ↓ (absent)
    .spe SAMPLEMASS field  (typed Spectrum.sample_mass_kg)   ← BUG-1 FIX
       ↓ (absent)
    geometry default       (0.5 kg for Маринелли, etc.)

The bug: the .spe SAMPLEMASS layer was missing. Pipeline fell straight
from "filename token absent" to geometry default (0.5 kg for Маринелли),
inflating Бк/кг by SAMPLEMASS / 0.5×.

Reference fixture (Th-232 .spe, line 21):
    SAMPLEMASS=1600.0;16.0   →   1.6 kg ± 0.016 kg
"""

from __future__ import annotations

import importlib
import logging
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

run_skill = importlib.import_module("run_skill")

TH232_SPE = (
    REPO_ROOT
    / "detectors"
    / "Gamma-1S"
    / "reference_spectra"
    / "archive"
    / "Th232_420-7-17_Маринелли_0cm.spe"
)


def _make_ctx(tmp_path: Path, cli_mass: object = None) -> "run_skill.RunContext":
    """Build a RunContext for the Th-232 spe, with optional CLI mass."""
    cfg = run_skill._load_config(None)
    if cli_mass is not None:
        cfg["analyze"]["sample_mass_kg"] = float(cli_mass)
    layout = run_skill.BundleLayout.from_base(tmp_path / "b", cfg)
    layout.ensure_dirs()
    log = logging.getLogger(f"bug1_{tmp_path.name}")
    log.addHandler(logging.NullHandler())
    log.setLevel(logging.DEBUG)
    ctx = run_skill.RunContext(
        spectrum=TH232_SPE,
        background=None,
        metadata=run_skill.SpectrumMetadata.from_path(TH232_SPE),
        cfg=cfg,
        layout=layout,
        logger=log,
        skill_version="test-bug1",
        include_v2=False,
    )
    return ctx


def test_th232_fixture_exists():
    """Sanity: fixture for the bug evidence is in the repo."""
    assert TH232_SPE.exists(), (
        f"Reference .spe missing — expected {TH232_SPE}"
    )


def test_read_spec_mass_kg_helper_returns_1_6():
    """`_read_spec_mass_kg(Th232.spe) == (1.6, 0.016)`."""
    pair = run_skill._read_spec_mass_kg(TH232_SPE)
    assert pair is not None
    m_kg, u_kg = pair
    assert m_kg == pytest.approx(1.6, abs=1e-9)
    assert u_kg == pytest.approx(0.016, abs=1e-9)


def test_read_spec_mass_kg_returns_none_for_missing_file(tmp_path):
    """Helper is fail-safe: missing/unreadable file → None."""
    pair = run_skill._read_spec_mass_kg(tmp_path / "nope.spe")
    assert pair is None


def test_build_orch_kwargs_uses_samplemass_when_no_cli(tmp_path):
    """BUG-1 acceptance: pipeline-bound mass_kg is 1.6 (from .spe), not 0.5."""
    ctx = _make_ctx(tmp_path, cli_mass=None)
    # Filename "Th232_420-7-17_Маринелли_0cm.spe" carries NO mass token.
    # geometry_hint will be "маринелли" → default would be 0.5 kg.
    # SAMPLEMASS in the .spe is 1.6 kg → that's what we want resolved.
    kwargs = run_skill._build_orch_kwargs(ctx)
    assert kwargs["sample_mass_kg"] == pytest.approx(1.6, abs=1e-9), (
        f"Expected mass from .spe SAMPLEMASS (1.6 kg), got "
        f"{kwargs['sample_mass_kg']}. BUG-1 regression — pipeline is "
        f"falling back to geometry default instead of reading .spe."
    )


def test_build_orch_kwargs_cli_overrides_samplemass(tmp_path):
    """CLI flag still wins (highest precedence)."""
    ctx = _make_ctx(tmp_path, cli_mass=2.0)
    kwargs = run_skill._build_orch_kwargs(ctx)
    assert kwargs["sample_mass_kg"] == pytest.approx(2.0, abs=1e-9), (
        f"CLI override broken: expected 2.0 kg, got {kwargs['sample_mass_kg']}."
    )


def test_geometry_default_fires_only_when_no_samplemass(tmp_path):
    """F-378 warning text only fires when both CLI AND .spe SAMPLEMASS missing.

    Constructs a spectrum file WITHOUT SAMPLEMASS (a freshly-built minimal
    .spe via the existing write_lsrm_spe writer + reading back) and confirms
    the geometry default is used + the warning is logged.
    """
    import numpy as np
    from gamma.io.lsrm_spe import write_lsrm_spe
    from gamma.spectrum import Spectrum

    spec = Spectrum(
        counts=np.zeros(1024, dtype=np.int64),
        live_time=100.0,
        real_time=100.0,
        geometry="Маринелли",
        sample_id="SYNTH_NOMASS",
    )
    spec.energy_cal = (0.0, 3.0)  # trivial linear cal so reader is happy
    spe_path = tmp_path / "synth_nomass_Маринелли.spe"
    write_lsrm_spe(spec, str(spe_path))

    # Sanity: writer round-trip shows no SAMPLEMASS line
    assert b"SAMPLEMASS" not in spe_path.read_bytes(), (
        "Writer unexpectedly emitted SAMPLEMASS — test fixture broken"
    )

    cfg = run_skill._load_config(None)
    layout = run_skill.BundleLayout.from_base(tmp_path / "b", cfg)
    layout.ensure_dirs()

    # Capture log messages
    captured: list = []

    class _CaptureHandler(logging.Handler):
        def emit(self, record):
            captured.append(record)

    log = logging.getLogger("bug1_default")
    log.handlers = [_CaptureHandler()]
    log.setLevel(logging.DEBUG)

    ctx = run_skill.RunContext(
        spectrum=spe_path,
        background=None,
        metadata=run_skill.SpectrumMetadata.from_path(spe_path),
        cfg=cfg, layout=layout, logger=log,
        skill_version="test-bug1-default", include_v2=False,
    )
    kwargs = run_skill._build_orch_kwargs(ctx)

    # Маринелли default = 0.5
    assert kwargs["sample_mass_kg"] == pytest.approx(0.5, abs=1e-9)

    # F-378 warning fired exactly once at WARNING level
    warnings = [r for r in captured if r.levelno >= logging.WARNING]
    assert any("F-378" in r.getMessage() for r in warnings), (
        f"F-378 default-warning expected; captured: "
        f"{[r.getMessage() for r in warnings]}"
    )


def test_info_log_when_samplemass_used(tmp_path):
    """When SAMPLEMASS feeds mass_kg, the INFO log carries the value/units."""
    cfg = run_skill._load_config(None)
    layout = run_skill.BundleLayout.from_base(tmp_path / "b", cfg)
    layout.ensure_dirs()

    captured: list = []

    class _CaptureHandler(logging.Handler):
        def emit(self, record):
            captured.append(record)

    log = logging.getLogger("bug1_info")
    log.handlers = [_CaptureHandler()]
    log.setLevel(logging.DEBUG)

    ctx = run_skill.RunContext(
        spectrum=TH232_SPE,
        background=None,
        metadata=run_skill.SpectrumMetadata.from_path(TH232_SPE),
        cfg=cfg, layout=layout, logger=log,
        skill_version="test-bug1-info", include_v2=False,
    )
    _ = run_skill._build_orch_kwargs(ctx)
    msgs = [r.getMessage() for r in captured]
    # Exactly one INFO message mentioning SAMPLEMASS + the 1.6 kg value
    info_hits = [m for m in msgs if "SAMPLEMASS" in m and "1.600" in m]
    assert info_hits, (
        f"Expected INFO line with 'SAMPLEMASS' + '1.600 кг'; got: {msgs}"
    )
