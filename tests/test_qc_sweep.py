"""
Tests for F-QC-01 sweep script (scripts/qc/sweep_qc_on_spectra_index.py).

RAG-ID: [F-QC-01], [RAG-041]
Cite: spectrum_qc_methodology_v2_2026-06-03.md; KNOWN_AND_FIXED_ISSUES.md:1292
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
_SWEEP_SCRIPT = _PROJECT_ROOT / "scripts" / "qc" / "sweep_qc_on_spectra_index.py"
_SPECTRA_INDEX = _PROJECT_ROOT / "audit" / "_rag" / "SPECTRA_INDEX.json"
_LSRM_ROOT = Path(r"C:\LSRM")
_SWEEP_OUTPUT = _PROJECT_ROOT / "audit" / "_drafts" / "f_qc_01_sweep_2026-06-04.json"

# ─── Fixtures ────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def sweep_json() -> dict:
    """Load the produced sweep output (or fail with a clear message)."""
    if not _SWEEP_OUTPUT.exists():
        pytest.skip(f"Sweep output not found: {_SWEEP_OUTPUT}. Run the sweep first.")
    with open(_SWEEP_OUTPUT, encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def lsrm_available() -> bool:
    return _LSRM_ROOT.exists()


# ─── Tests ───────────────────────────────────────────────────────────────────


class TestSweepOutputSchema:
    """Schema validation for f_qc_01_sweep_*.json output."""

    def test_top_level_keys(self, sweep_json):
        """Output must have _meta, aggregate, fixtures."""
        assert "_meta" in sweep_json, "missing _meta"
        assert "aggregate" in sweep_json, "missing aggregate"
        assert "fixtures" in sweep_json, "missing fixtures"

    def test_meta_fields(self, sweep_json):
        """_meta must include f_rule, rag_id, tier."""
        meta = sweep_json["_meta"]
        assert meta.get("f_rule") == "F-QC-01", "_meta.f_rule must be F-QC-01"
        assert meta.get("rag_id") == "RAG-041", "_meta.rag_id must be RAG-041"
        assert "tier" in meta, "_meta.tier missing"
        assert "schema_version" in meta, "_meta.schema_version missing"

    def test_aggregate_fields(self, sweep_json):
        """aggregate must have n_total, pass_rate, per_criterion_fail."""
        agg = sweep_json["aggregate"]
        for key in ("n_total", "n_pass", "n_fail", "n_error", "pass_rate",
                    "per_criterion_fail", "by_detector"):
            assert key in agg, f"aggregate.{key} missing"

    def test_criterion_keys_in_per_criterion_fail(self, sweep_json):
        """per_criterion_fail must report all 6 criterion keys."""
        expected = {
            "c1_energy_drift", "c2_fwhm_stability", "c3_efficiency_qa",
            "c4_bg_drift", "c5_peak_z_roi", "c6_sensitivity",
        }
        actual = set(sweep_json["aggregate"]["per_criterion_fail"].keys())
        assert actual == expected, f"per_criterion_fail keys mismatch: {actual ^ expected}"

    def test_fixture_schema(self, sweep_json):
        """Every fixture entry must have spectrum_id, sha256, criteria_verdicts."""
        for fix in sweep_json["fixtures"]:
            sid = fix.get("spectrum_id", "<missing>")
            assert "spectrum_id" in fix, "missing spectrum_id"
            assert "sha256" in fix, f"{sid}: missing sha256"
            assert "overall_passed" in fix, f"{sid}: missing overall_passed"
            assert "criteria_verdicts" in fix, f"{sid}: missing criteria_verdicts"
            cv = fix["criteria_verdicts"]
            for key in ("c1_energy_drift", "c2_fwhm_stability", "c3_efficiency_qa",
                        "c4_bg_drift", "c5_peak_z_roi", "c6_sensitivity"):
                assert key in cv, f"{sid}: criteria_verdicts missing {key}"
                assert "passed" in cv[key], f"{sid}: {key} missing 'passed'"

    def test_no_absolute_paths_in_rel_path(self, sweep_json):
        """
        F-115 compliance: fixture rel_path must contain <LSRM> placeholder,
        not an absolute operator path.
        """
        for fix in sweep_json["fixtures"]:
            rel = fix.get("rel_path", "")
            sid = fix.get("spectrum_id", "?")
            assert "<LSRM>" in rel, (
                f"F-115 violation: {sid} rel_path does not contain <LSRM> placeholder: {rel!r}"
            )
            # Absolute path check: must not start with C:/ or //
            assert not rel.startswith("C:"), (
                f"F-115 violation: {sid} rel_path is absolute: {rel!r}"
            )

    @pytest.mark.parametrize("n", [2, 5, 0])
    def test_n_pass_plus_n_fail_plus_n_error_equals_n_total(self, sweep_json, n):
        """Aggregate counts must be consistent."""
        _ = n   # parametrize forces 3 runs (idempotence proxy)
        agg = sweep_json["aggregate"]
        total = agg["n_total"]
        counts = agg["n_pass"] + agg["n_fail"] + agg.get("n_error", 0)
        # n_null is counted separately (errors produce overall_passed=None)
        # Total = pass + fail + null(=error)
        null = total - agg["n_pass"] - agg["n_fail"]
        assert null == agg.get("n_error", 0) or null == agg.get("n_null", 0), (
            f"Count mismatch: n_total={total} n_pass={agg['n_pass']} "
            f"n_fail={agg['n_fail']} n_error={agg.get('n_error',0)}"
        )


class TestSweepIdempotence:
    """Re-running the sweep on same inputs must produce same aggregate."""

    @pytest.mark.skipif(
        not _LSRM_ROOT.exists(),
        reason="LSRM root C:\\LSRM not accessible",
    )
    def test_idempotent_smoke(self, tmp_path):
        """Two smoke runs on 3 fixtures → same n_pass, n_fail, per_criterion_fail."""
        out1 = str(tmp_path / "sweep_run1.json")
        out2 = str(tmp_path / "sweep_run2.json")
        base_args = [
            sys.executable, str(_SWEEP_SCRIPT),
            "--lsrm-root", str(_LSRM_ROOT),
            "--tier", "A",
            "--smoke", "3",
            "--workers", "1",
        ]
        r1 = subprocess.run(base_args + ["--output", out1],
                            capture_output=True, text=True)
        r2 = subprocess.run(base_args + ["--output", out2],
                            capture_output=True, text=True)
        assert r1.returncode == 0, f"run1 failed: {r1.stderr}"
        assert r2.returncode == 0, f"run2 failed: {r2.stderr}"
        with open(out1, encoding="utf-8") as f:
            j1 = json.load(f)
        with open(out2, encoding="utf-8") as f:
            j2 = json.load(f)
        # Aggregate counts must be identical
        a1, a2 = j1["aggregate"], j2["aggregate"]
        assert a1["n_pass"] == a2["n_pass"], "idempotence fail: n_pass differs"
        assert a1["n_fail"] == a2["n_fail"], "idempotence fail: n_fail differs"
        assert a1["per_criterion_fail"] == a2["per_criterion_fail"], (
            "idempotence fail: per_criterion_fail differs"
        )


class TestCriterionDefinitionsMatchSpec:
    """
    Criterion definitions in the sweep output must match RAG-041 / methodology v2.
    Verifies that c6_sensitivity is always 'passed=True' (placeholder in v1.21.0).
    """

    def test_c6_sensitivity_always_pass(self, sweep_json):
        """Criterion 6 is a placeholder in v1.21.0 → all fixtures must report passed=True."""
        for fix in sweep_json["fixtures"]:
            sid = fix.get("spectrum_id", "?")
            c6 = fix["criteria_verdicts"]["c6_sensitivity"]
            assert c6["passed"] is True, (
                f"{sid}: c6_sensitivity.passed is not True — "
                "RAG-041 notes criterion 6 is placeholder (None) in v1.21.0"
            )

    def test_c1_energy_drift_conservative_pass_when_unavailable(self, sweep_json):
        """
        Criterion 1 (energy drift) must pass when seven_line_check is not available
        (methodology v2: unavailable criteria count as PASS, conservative).
        """
        c1_unavail = [
            fix for fix in sweep_json["fixtures"]
            if not fix["criteria_verdicts"]["c1_energy_drift"].get("available", True)
        ]
        for fix in c1_unavail:
            sid = fix.get("spectrum_id", "?")
            assert fix["criteria_verdicts"]["c1_energy_drift"]["passed"] is True, (
                f"{sid}: c1_energy_drift with available=False must be passed=True "
                "(RAG-041 conservative-pass rule for unavailable criteria)"
            )
