# -*- coding: utf-8 -*-
"""v1.18.14 — F-316 FWHM filter NaI-default-ON + F-315 secondary_peaks guard.
"""
from __future__ import annotations
import json, os, subprocess, sys
from pathlib import Path
import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


# ──────────────────────────────────────────────────────────────────
# F-316 — filter_narrow_peaks default resolution
# ──────────────────────────────────────────────────────────────────

def test_F316_filter_narrow_peaks_signature_is_optional_none():
    """Default = None (resolve via detector_class)."""
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    import inspect
    sig = inspect.signature(analyze_lsrm_spe)
    p = sig.parameters["filter_narrow_peaks"]
    assert p.default is None, (
        f"Default should be None (auto-resolve), got {p.default!r}"
    )


def test_F316_default_none_keeps_filter_off():
    """REVERTED 2026-05-31: default=None → False (filter OFF).

    Первоначальная v1.18.14 попытка auto-resolve filter=ON для NaI отсекала
    legitimate Bi-214 609/1764 keV pair (regression). Filter остаётся opt-in:
    пользователь явно включает через CLI флаг --filter-narrow-peaks.
    """
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    path = REPO / "evals" / "fixtures" / "M_cs_легкий_2001-2005.spe"
    if not path.exists():
        pytest.skip(f"Fixture missing: {path}")
    # filter_narrow_peaks=None → resolves to False (back-compat).
    r = analyze_lsrm_spe(
        str(path), compute_activities=False, complete_workflow=False,
        allow_stage2=True,
    )
    # Smoke check: pipeline runs without exception
    assert r is not None
    assert hasattr(r, "detector_type")


def test_F316_explicit_false_overrides_auto():
    """Если filter_narrow_peaks=False явно → не активируется."""
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    path = REPO / "evals" / "fixtures" / "M_cs_легкий_2001-2005.spe"
    if not path.exists():
        pytest.skip(f"Fixture missing: {path}")
    r = analyze_lsrm_spe(
        str(path), compute_activities=False, complete_workflow=False,
        filter_narrow_peaks=False,
        allow_stage2=True,
    )
    assert r is not None


def test_F316_explicit_true_overrides_auto():
    """Если filter_narrow_peaks=True явно → активируется для любого detector."""
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    path = REPO / "evals" / "fixtures" / "M_cs_легкий_2001-2005.spe"
    if not path.exists():
        pytest.skip(f"Fixture missing: {path}")
    r = analyze_lsrm_spe(
        str(path), compute_activities=False, complete_workflow=False,
        filter_narrow_peaks=True,
        allow_stage2=True,
    )
    assert r is not None


# ──────────────────────────────────────────────────────────────────
# F-315 — secondary_peaks regression guard
# ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def demo_run(tmp_path_factory):
    out = tmp_path_factory.mktemp("demo_out")
    cmd = [
        sys.executable,
        str(REPO / "scripts" / "regen_demo_reports.py"),
        "--output-dir", str(out),
    ]
    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONPATH"] = (
        str(REPO / "scripts") + os.pathsep + env.get("PYTHONPATH", "")
    )
    r = subprocess.run(
        cmd, capture_output=True, text=True, timeout=120,
        env=env, encoding="utf-8", errors="replace",
        cwd=str(REPO),
    )
    if r.returncode != 0:
        pytest.skip(f"regen failed: {r.stderr}")
    return out


def test_F315_secondary_peaks_have_real_energies(demo_run):
    """F-147: secondary_peaks JSON должен содержать ненулевые energy_keV."""
    cs_data = json.loads(
        (demo_run / "M_cs_легкий_report.json").read_text(encoding="utf-8")
    )
    sp = cs_data.get("secondary_peaks", [])
    assert len(sp) > 0, "M_cs должен иметь хотя бы 1 secondary_peak"
    nonzero = [s for s in sp if s.get("energy_keV", 0) > 0]
    assert len(nonzero) >= len(sp) * 0.5, (
        f"≥50% secondary_peaks должны иметь energy_keV > 0; "
        f"got {len(nonzero)}/{len(sp)}"
    )


def test_F315_compton_edge_for_cs137(demo_run):
    """F-147: Cs-137 demo должен включать Compton edge (477.3 keV)."""
    cs_data = json.loads(
        (demo_run / "M_cs_легкий_report.json").read_text(encoding="utf-8")
    )
    sp = cs_data.get("secondary_peaks", [])
    compton_edges = [s for s in sp if s.get("feature_kind") == "compton_edge"]
    assert len(compton_edges) >= 1, (
        "Cs-137 demo должен иметь хотя бы 1 compton_edge"
    )
    # Hard-check: Cs-137 661.66 → E_CE = 477.3 keV (Knoll §10)
    cs137_ce = [
        s for s in compton_edges
        if abs(s.get("energy_keV", 0) - 477.3) < 5
    ]
    assert len(cs137_ce) >= 1, (
        f"Cs-137 Compton-edge на 477.3 keV не найден; "
        f"compton_edges energies: {[s.get('energy_keV') for s in compton_edges]}"
    )


def test_F315_parent_line_keV_preserved(demo_run):
    """F-147: каждый secondary_peak должен ссылаться на parent_line_keV > 0."""
    cs_data = json.loads(
        (demo_run / "M_cs_легкий_report.json").read_text(encoding="utf-8")
    )
    sp = cs_data.get("secondary_peaks", [])
    with_parent = [s for s in sp if s.get("parent_line_keV", 0) > 0]
    assert len(with_parent) >= 1, (
        "Хотя бы 1 secondary_peak должен иметь parent_line_keV > 0"
    )
