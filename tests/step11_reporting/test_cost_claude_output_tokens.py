"""P0-8 / v1.25.0 — CostEstimate.claude_output_tokens field + 20k alarm.

Contracts:
  * CostEstimate dataclass has field claude_output_tokens: int = 0.
  * report.json["cost_estimate"]["claude_output_tokens"] is present,
    is an int, and is >= 0.
  * When cost_estimate["output_tokens"] >= 20_000 is passed to build_report,
    report.json["warnings"] contains an entry with code COST_HIGH_OUTPUT_TOKENS.
  * When output_tokens < 20_000 (or 0), the alarm entry is absent.
"""
from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.reporting.cost_estimator import (
    CostEstimate,
    StageCostEstimate,
    estimate_total_cost,
    DEFAULT_SESSION_TOKEN_BUDGET,
)


FIXTURE_DIR = (
    Path(__file__).parent.parent.parent
    / "detectors" / "Gamma-1S" / "reference_spectra" / "archive"
)
CS137_FIXTURE = FIXTURE_DIR / "Cs137_420-7-14_Маринелли_0cm.spe"


def _need(p: Path):
    if not p.exists():
        pytest.skip(f"fixture missing: {p}")


# ─── Unit-level: CostEstimate dataclass ─────────────────────────────────────

def test_cost_estimate_has_claude_output_tokens_field():
    """CostEstimate dataclass must have claude_output_tokens field defaulting to 0."""
    ce = CostEstimate(
        tokens_total=1000,
        session_token_budget=DEFAULT_SESSION_TOKEN_BUDGET,
        session_pct=0.5,
    )
    assert hasattr(ce, "claude_output_tokens")
    assert ce.claude_output_tokens == 0


def test_cost_estimate_claude_output_tokens_type():
    """claude_output_tokens must be an int."""
    ce = CostEstimate(
        tokens_total=1000,
        session_token_budget=DEFAULT_SESSION_TOKEN_BUDGET,
        session_pct=0.5,
        claude_output_tokens=5000,
    )
    assert isinstance(ce.claude_output_tokens, int)
    assert ce.claude_output_tokens == 5000


def test_cost_estimate_claude_output_tokens_non_negative():
    """claude_output_tokens must be >= 0."""
    ce = CostEstimate(
        tokens_total=500,
        session_token_budget=DEFAULT_SESSION_TOKEN_BUDGET,
        session_pct=0.25,
        claude_output_tokens=0,
    )
    assert ce.claude_output_tokens >= 0


def test_to_dict_includes_claude_output_tokens():
    """CostEstimate.to_dict() must include claude_output_tokens as int."""
    ce = CostEstimate(
        tokens_total=8000,
        session_token_budget=DEFAULT_SESSION_TOKEN_BUDGET,
        session_pct=4.0,
        claude_output_tokens=7500,
    )
    d = ce.to_dict()
    assert "claude_output_tokens" in d
    assert isinstance(d["claude_output_tokens"], int)
    assert d["claude_output_tokens"] == 7500


def test_to_dict_claude_output_tokens_default_zero():
    """to_dict() claude_output_tokens defaults to 0 when not set."""
    ce = CostEstimate(
        tokens_total=3000,
        session_token_budget=DEFAULT_SESSION_TOKEN_BUDGET,
        session_pct=1.5,
    )
    d = ce.to_dict()
    assert d["claude_output_tokens"] == 0


def test_dataclasses_replace_on_frozen():
    """CostEstimate is frozen=True; dataclasses.replace must work."""
    ce = CostEstimate(
        tokens_total=2000,
        session_token_budget=DEFAULT_SESSION_TOKEN_BUDGET,
        session_pct=1.0,
        claude_output_tokens=0,
    )
    ce2 = dataclasses.replace(ce, claude_output_tokens=25000)
    assert ce2.claude_output_tokens == 25000
    # Original unchanged (frozen)
    assert ce.claude_output_tokens == 0


# ─── Integration: build_report wiring ────────────────────────────────────────

def _make_result():
    """Return a StagedAnalysisResult from CS-137 fixture for integration tests."""
    _need(CS137_FIXTURE)
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    return analyze_lsrm_spe(str(CS137_FIXTURE), complete_workflow=True)


def test_build_report_cost_estimate_has_claude_output_tokens():
    """build_report json_dict must include cost_estimate.claude_output_tokens (int, >= 0)."""
    result = _make_result()
    from gamma.reporting.build import build_report
    out = build_report(result)
    ce = out["json_dict"].get("cost_estimate")
    assert ce is not None, "cost_estimate missing from json_dict"
    assert "claude_output_tokens" in ce, "claude_output_tokens missing from cost_estimate"
    assert isinstance(ce["claude_output_tokens"], int)
    assert ce["claude_output_tokens"] >= 0


def test_build_report_default_claude_output_tokens_is_zero():
    """When no output_tokens is passed, claude_output_tokens in report defaults to 0."""
    result = _make_result()
    from gamma.reporting.build import build_report
    out = build_report(result)
    ce = out["json_dict"]["cost_estimate"]
    assert ce["claude_output_tokens"] == 0


def test_build_report_wires_output_tokens_from_cost_estimate_param():
    """When cost_estimate dict passes output_tokens=5000, report shows 5000."""
    result = _make_result()
    from gamma.reporting.build import build_report
    out = build_report(result, cost_estimate={"output_tokens": 5000})
    ce = out["json_dict"]["cost_estimate"]
    assert ce["claude_output_tokens"] == 5000


# ─── 20k alarm tests ─────────────────────────────────────────────────────────

def test_alarm_triggered_when_output_tokens_at_20k_threshold():
    """Alarm COST_HIGH_OUTPUT_TOKENS fires when output_tokens == 20000."""
    result = _make_result()
    from gamma.reporting.build import build_report
    out = build_report(result, cost_estimate={"output_tokens": 20000})
    warnings = out["json_dict"].get("warnings", [])
    alarm_entries = [
        w for w in warnings
        if isinstance(w, dict) and w.get("code") == "COST_HIGH_OUTPUT_TOKENS"
    ]
    assert alarm_entries, (
        "Expected COST_HIGH_OUTPUT_TOKENS warning at exactly 20000 tokens"
    )
    alarm = alarm_entries[0]
    assert alarm["claude_output_tokens"] == 20000
    assert alarm["threshold"] == 20_000
    assert alarm["severity"] == "INFO"


def test_alarm_triggered_at_25000_tokens():
    """Alarm COST_HIGH_OUTPUT_TOKENS fires when output_tokens == 25000."""
    result = _make_result()
    from gamma.reporting.build import build_report
    out = build_report(result, cost_estimate={"output_tokens": 25000})
    warnings = out["json_dict"].get("warnings", [])
    alarm_entries = [
        w for w in warnings
        if isinstance(w, dict) and w.get("code") == "COST_HIGH_OUTPUT_TOKENS"
    ]
    assert alarm_entries, (
        "Expected COST_HIGH_OUTPUT_TOKENS warning for 25000 output tokens"
    )
    alarm = alarm_entries[0]
    assert alarm["claude_output_tokens"] == 25000
    assert "20000" in alarm["message"] or "20_000" in alarm["message"] or "20000" in str(alarm["threshold"])


def test_alarm_absent_below_threshold():
    """No alarm when output_tokens < 20000."""
    result = _make_result()
    from gamma.reporting.build import build_report
    out = build_report(result, cost_estimate={"output_tokens": 19999})
    warnings = out["json_dict"].get("warnings", [])
    alarm_entries = [
        w for w in warnings
        if isinstance(w, dict) and w.get("code") == "COST_HIGH_OUTPUT_TOKENS"
    ]
    assert not alarm_entries, (
        "COST_HIGH_OUTPUT_TOKENS should NOT fire for 19999 < 20000"
    )


def test_alarm_absent_when_output_tokens_zero():
    """No alarm when output_tokens == 0 (default / unknown)."""
    result = _make_result()
    from gamma.reporting.build import build_report
    out = build_report(result)
    warnings = out["json_dict"].get("warnings", [])
    alarm_entries = [
        w for w in warnings
        if isinstance(w, dict) and w.get("code") == "COST_HIGH_OUTPUT_TOKENS"
    ]
    assert not alarm_entries, (
        "COST_HIGH_OUTPUT_TOKENS should NOT fire when claude_output_tokens == 0"
    )


def test_alarm_message_contains_token_count():
    """Alarm message includes the actual token count."""
    result = _make_result()
    from gamma.reporting.build import build_report
    out = build_report(result, cost_estimate={"output_tokens": 30000})
    warnings = out["json_dict"].get("warnings", [])
    alarm = next(
        (w for w in warnings
         if isinstance(w, dict) and w.get("code") == "COST_HIGH_OUTPUT_TOKENS"),
        None,
    )
    assert alarm is not None
    assert "30000" in alarm["message"]
