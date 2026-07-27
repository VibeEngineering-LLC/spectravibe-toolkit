# F-RPT-08 / v1.19.1 — Smoke test: ISO 11929 fields emitted in JSON report
#
# Verifies that build_json_report() emits
#   decision_threshold_Bq_per_kg  (float or null)
#   detection_limit_Bq_per_kg     (float or null)
# in every element of the identified_nuclides list.
#
# Test strategy: build a minimal StagedAnalysisResult using the
# staged pipeline on a real fixture (Cs-137), then check the JSON output.
# Skipped when evals/fixtures are absent (CI without large test data).

from __future__ import annotations

import importlib
import json
import pathlib
import sys
from typing import Optional

import pytest

# ---------------------------------------------------------------------------
# Path plumbing — same pattern as test_run_skill_orchestration.py
# ---------------------------------------------------------------------------
REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

FIXTURES = REPO_ROOT / "evals" / "fixtures"
CS137_FIXTURE = FIXTURES / "M_cs_легкий_2001-2005.spe"
DEFAULT_BG = (
    REPO_ROOT
    / "detectors"
    / "Gamma-1S"
    / "data"
    / "averaged_backgrounds"
    / "bg_2016_marinelli_water_marinelli.spe"
)


def _have_fixtures() -> bool:
    return CS137_FIXTURE.exists() and DEFAULT_BG.exists()


# ---------------------------------------------------------------------------
# Direct import of build_json_report + staged pipeline
# ---------------------------------------------------------------------------
from gamma.reporting.json_report import build_json_report  # noqa: E402
from gamma.identification.staged_pipeline import analyze_lsrm_spe  # noqa: E402


@pytest.fixture(scope="module")
def cs137_json_report():
    """Run the Cs-137 fixture through analyze_lsrm_spe + build_json_report.

    Returns the full report dict.  Fixture is module-scoped to avoid
    running the pipeline more than once per session.
    """
    result = analyze_lsrm_spe(
        str(CS137_FIXTURE),
        background_path=str(DEFAULT_BG),
        sample_mass_kg=0.570,
        compute_activities=True,
        compute_mda=True,
        allow_stage2=True,
        allow_stage3=True,
    )
    return build_json_report(result)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _have_fixtures(), reason="evals/fixtures missing")
def test_identified_nuclides_have_decision_threshold_key(cs137_json_report):
    """Every element in identified_nuclides must have decision_threshold_Bq_per_kg key.

    F-RPT-08 / v1.19.1 — the key must exist (even if the value is None when
    mass or MDA data are unavailable).  Key absence would indicate a wiring
    regression in json_report.py _build_identified_nuclides().
    """
    ids = cs137_json_report.get("identified_nuclides", [])
    assert len(ids) > 0, "Expected at least one identified nuclide in Cs-137 fixture"
    for entry in ids:
        assert "decision_threshold_Bq_per_kg" in entry, (
            f"Missing 'decision_threshold_Bq_per_kg' key in nuclide '{entry.get('nuclide')}'"
        )


@pytest.mark.skipif(not _have_fixtures(), reason="evals/fixtures missing")
def test_identified_nuclides_have_detection_limit_key(cs137_json_report):
    """Every element in identified_nuclides must have detection_limit_Bq_per_kg key.

    F-RPT-08 / v1.19.1 — analogous to decision_threshold check above.
    """
    ids = cs137_json_report.get("identified_nuclides", [])
    assert len(ids) > 0
    for entry in ids:
        assert "detection_limit_Bq_per_kg" in entry, (
            f"Missing 'detection_limit_Bq_per_kg' key in nuclide '{entry.get('nuclide')}'"
        )


@pytest.mark.skipif(not _have_fixtures(), reason="evals/fixtures missing")
def test_iso_threshold_values_are_float_or_none(cs137_json_report):
    """decision_threshold_Bq_per_kg and detection_limit_Bq_per_kg must be float or None.

    F-RPT-08 / v1.19.1 — JSON schema requirement: no NaN, no Inf, no string.
    Values may be None when mass or MDA data are absent.
    """
    ids = cs137_json_report.get("identified_nuclides", [])
    for entry in ids:
        nuc = entry.get("nuclide", "?")
        dt = entry.get("decision_threshold_Bq_per_kg")
        dl = entry.get("detection_limit_Bq_per_kg")
        assert dt is None or isinstance(dt, (int, float)), (
            f"[{nuc}] decision_threshold_Bq_per_kg={dt!r} is not float or None"
        )
        assert dl is None or isinstance(dl, (int, float)), (
            f"[{nuc}] detection_limit_Bq_per_kg={dl!r} is not float or None"
        )
        # NaN / Inf must never appear (json_report._safe_float guarantee)
        if isinstance(dt, float):
            import math
            assert not math.isnan(dt) and not math.isinf(dt), (
                f"[{nuc}] decision_threshold_Bq_per_kg is NaN/Inf"
            )
        if isinstance(dl, float):
            import math
            assert not math.isnan(dl) and not math.isinf(dl), (
                f"[{nuc}] detection_limit_Bq_per_kg is NaN/Inf"
            )


@pytest.mark.skipif(not _have_fixtures(), reason="evals/fixtures missing")
def test_iso_detection_limit_ge_decision_threshold(cs137_json_report):
    """When both fields are non-None floats, detection_limit >= decision_threshold.

    ISO 11929-1:2019 §5.4.4 guarantee: y# >= y* always.
    The approximation y# ≈ 2·y* ensures y# = 2·y* > y* when y* > 0.
    """
    ids = cs137_json_report.get("identified_nuclides", [])
    for entry in ids:
        nuc = entry.get("nuclide", "?")
        dt = entry.get("decision_threshold_Bq_per_kg")
        dl = entry.get("detection_limit_Bq_per_kg")
        if dt is not None and dl is not None:
            assert dl >= dt, (
                f"[{nuc}] detection_limit={dl:.4f} < decision_threshold={dt:.4f} "
                f"— violates ISO 11929 ordering"
            )


@pytest.mark.skipif(not _have_fixtures(), reason="evals/fixtures missing")
def test_json_serializable_after_iso_fields_added(cs137_json_report):
    """build_json_report() output must be fully JSON-serializable including new fields.

    Regression guard: NaN / Inf / non-serializable types in new fields would cause
    json.dumps() to raise TypeError or ValueError.
    """
    # Should not raise
    serialized = json.dumps(cs137_json_report, ensure_ascii=False)
    assert '"decision_threshold_Bq_per_kg"' in serialized
    assert '"detection_limit_Bq_per_kg"' in serialized


@pytest.mark.skipif(not _have_fixtures(), reason="evals/fixtures missing")
def test_cs137_iso_thresholds_physical_order_of_magnitude(cs137_json_report):
    """Cs-137 decision threshold should be non-zero when mass and MDA data are present.

    For a 570 g Cs-137 sample with 3600 s live time and NaI 50×50 mm detector
    (efficiency ~2% at 662 keV), the decision threshold should be in the range
    [0.001, 100] Bq/kg — a deliberately wide range that catches gross formula errors
    (e.g., units confusion kg vs g, efficiency in pct instead of decimal).
    """
    ids = cs137_json_report.get("identified_nuclides", [])
    cs137_entry = next(
        (e for e in ids if e.get("nuclide") == "Cs-137"), None
    )
    if cs137_entry is None:
        pytest.skip("Cs-137 not identified in this run — pipeline-level issue, not RPT-08")
    dt = cs137_entry.get("decision_threshold_Bq_per_kg")
    if dt is None:
        pytest.skip("decision_threshold_Bq_per_kg is None — no MDA data for Cs-137")
    assert 0.001 <= dt <= 100.0, (
        f"Cs-137 decision_threshold_Bq_per_kg={dt:.4f} outside physical range [0.001, 100]"
    )
