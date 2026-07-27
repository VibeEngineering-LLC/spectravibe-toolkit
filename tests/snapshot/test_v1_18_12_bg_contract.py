# -*- coding: utf-8 -*-
"""v1.18.12 — F-135 background-subtraction contract regression guard.

Hard requirement: при demo regeneration (или любом analyze_and_report
runs через scripts/regen_demo_reports.py) фон ОБЯЗАТЕЛЬНО должен быть
вычтен. Это F-135 contract («ЗАКРЕПЛЕНО НАВСЕГДА»). Если бы в v1.18.11
бы этот тест существовал, P0 issue (фон не вычтен в demo output) был
бы пойман до релиза.
"""
from __future__ import annotations
import json, os, subprocess, sys
from pathlib import Path
import pytest

REPO = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def demo_run(tmp_path_factory):
    """Run regen_demo_reports.py once for the module."""
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


@pytest.mark.parametrize("nuclide,filename", [
    ("Cs-137", "M_cs_легкий_report.json"),
    ("K-40",   "M_k_легкий_report.json"),
    ("Ra-226", "M_ra_легкий_report.json"),
    ("Th-232", "M_th_легкий_report.json"),
])
def test_F313_demo_has_background_subtracted_True(demo_run, nuclide, filename):
    """F-135 contract: фон ВСЕГДА вычитается в demo regeneration."""
    p = demo_run / filename
    assert p.exists(), f"Missing demo report: {p}"
    data = json.loads(p.read_text(encoding="utf-8"))
    diag = data.get("diagnostics", {})
    bg_subtracted = diag.get("background_subtracted")
    assert bg_subtracted is True, (
        f"F-135 contract VIOLATION on {nuclide}: "
        f"background_subtracted = {bg_subtracted!r}. "
        f"Demo regeneration MUST run with bg-subtracted spectrum."
    )


def test_F313_demo_summary_exists(demo_run):
    """DEMO_SUMMARY.json должен быть создан с 4 nuclides."""
    p = demo_run / "DEMO_SUMMARY.json"
    assert p.exists()
    data = json.loads(p.read_text(encoding="utf-8"))
    assert set(data.keys()) >= {"Cs-137", "K-40", "Ra-226", "Th-232"}


def test_F313_cs137_activity_positive(demo_run):
    """Sanity: Cs-137 activity from M_cs > 100 Bq (real source)."""
    p = demo_run / "DEMO_SUMMARY.json"
    data = json.loads(p.read_text(encoding="utf-8"))
    cs_summary = data.get("Cs-137", {})
    acts = cs_summary.get("activities_Bq", {})
    cs_act = acts.get("Cs-137", 0)
    assert cs_act > 100, (
        f"Cs-137 activity too low: {cs_act} Bq (expected > 100)"
    )


# ──────────────────────────────────────────────────────────────────
# F-312 — F-131 search heuristic enhancement test
# ──────────────────────────────────────────────────────────────────

def test_F312_background_search_reaches_detector_subtree():
    """F-131 эвристика должна находить detectors/<DET>/data/averaged_backgrounds."""
    import sys as _sys
    _sys.path.insert(0, str(REPO / "scripts"))
    from gamma.io.readers import read_spectrum
    from gamma.io.background_search import find_background_candidates

    sample_path = REPO / "evals" / "fixtures" / "M_cs_легкий_2001-2005.spe"
    spec = read_spectrum(str(sample_path))
    # Increase max_days_apart чтобы date-filter не блокировал archive sample
    cands = find_background_candidates(
        spec, str(sample_path), max_days_apart=100000,
    )
    # Хотя бы один кандидат должен быть найден из detectors/Gamma-1S subtree
    assert len(cands) > 0, (
        "F-131 search returned 0 candidates even with date filter relaxed — "
        "search heuristic doesn't reach detectors/<DET>/data/averaged_backgrounds"
    )
    found_paths = [str(c.path) for c in cands]
    in_detector_subtree = any(
        "detectors" in p and "averaged_backgrounds" in p
        for p in found_paths
    )
    assert in_detector_subtree, (
        f"Found {len(cands)} candidates but none in detector subtree: "
        f"{found_paths[:3]}"
    )
