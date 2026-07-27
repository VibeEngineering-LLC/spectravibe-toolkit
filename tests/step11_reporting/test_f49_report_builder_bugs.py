"""F49 unit tests — report-builder logic bugs (F4 + F5 + F6).

F4: cert_zcheck verdict при A=None и sA=None ДОЛЖЕН быть "n/a", не
    "upper-limit" (последний зарезервирован за реальным is_upper_limit =
    fitted с S<L_D).
F5: confidence_level() для CI∈[3, 5) ДОЛЖЕН возвращать "low" (порог 5.0;
    раньше 3.0 → нуклиды Cs-134/Na-22 выпадали в "moderate" при единичной
    линии — жёсткий verdict без low-confidence пометки).
F6: completeness_dc_pct в JSON-отчёте ДОЛЖЕН быть числом, когда completeness
    посчитан (был баг: getattr(cmp, "DC_pct") вместо "dc_percent" →
    атрибут не существует → JSON всегда None → markdown рендерил DC="—"
    но flag="полно", что вводило оператора в заблуждение).
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

from gamma.identification.completeness import (
    CompletenessResult,
    compute_completeness,
)
from gamma.identification.confidence import ConfidenceIndexResult


# ──────────────────────────────────────────────────────────────────
# F5 — confidence_level() new low/moderate threshold = 5.0
# ──────────────────────────────────────────────────────────────────

def _make_ci(value: float, nuclide: str = "test") -> ConfidenceIndexResult:
    return ConfidenceIndexResult(
        nuclide=nuclide,
        n_lines_used=2,
        energy_uncertainty_factor=1e-4,
        intensity_uncertainty_factor=1e-1,
        CI=value,
    )


def test_f5_ci_below_3_is_low():
    assert _make_ci(1.8, "Cs-137").confidence_level() == "low"


def test_f5_ci_in_3_to_5_is_low_after_fix():
    """CI∈[3,5) теперь low (раньше moderate под порогом 3.0)."""
    assert _make_ci(3.5, "Na-22").confidence_level() == "low"
    assert _make_ci(4.4, "Cs-134").confidence_level() == "low"
    assert _make_ci(4.99).confidence_level() == "low"


def test_f5_ci_5_is_moderate():
    assert _make_ci(5.0, "boundary").confidence_level() == "moderate"
    assert _make_ci(5.9, "Co-60").confidence_level() == "moderate"
    assert _make_ci(8.5, "Ba-133").confidence_level() == "moderate"


def test_f5_ci_at_or_above_10_is_high():
    assert _make_ci(10.0, "boundary").confidence_level() == "high"
    assert _make_ci(18.3, "Eu-152").confidence_level() == "high"


# ──────────────────────────────────────────────────────────────────
# F6 — completeness.dc_percent attribute / JSON wiring
# ──────────────────────────────────────────────────────────────────

def test_f6_completeness_result_attribute_naming():
    """Канонический атрибут — dc_percent (lowercase). DC_pct НЕ существует."""
    r = compute_completeness(detected_nuclides=[], unmatched_peaks_E_area=[])
    assert hasattr(r, "dc_percent"), "dc_percent missing — атрибут переименован?"
    assert not hasattr(r, "DC_pct"), \
        "DC_pct появилось — старый bug-ловушка возвращается?"
    assert r.dc_percent == 0.0
    assert r.flag == "n/a"


def test_f6_completeness_dc_percent_non_zero():
    class FakeMatch:
        def __init__(self, area, E):
            self.peak_area = area
            self.library_E_keV = E

    class FakeNid:
        def __init__(self, matched):
            self.matched_lines = matched

    nid = FakeNid([FakeMatch(1000.0, 500.0)])
    unmatched = [(800.0, 200.0)]
    r = compute_completeness([nid], unmatched)
    assert r.dc_percent > 0.0
    assert getattr(r, "dc_percent", None) is not None
    assert getattr(r, "DC_pct", None) is None


# ──────────────────────────────────────────────────────────────────
# F4 — cert_zcheck verdict 'n/a' when A=None and sA=None
# ──────────────────────────────────────────────────────────────────

def _load_cert_zcheck():
    here = Path(__file__).resolve()
    cert_path = here.parents[2] / "scripts" / "cert_zcheck.py"
    spec = importlib.util.spec_from_file_location("cert_zcheck_mod", cert_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_f4_verdict_na_when_a_and_sigma_none(tmp_path: Path):
    mod = _load_cert_zcheck()
    report = tmp_path / "fake_report.json"
    report.write_text(json.dumps({
        "identified_nuclides": [
            {
                "nuclide": "Bi-212",
                "tier": "tentative",
                "characteristic_line_keV": 727.3,
                "specific_activity_Bq_per_kg": None,
                "specific_activity_sigma_Bq_per_kg": None,
            },
            {
                "nuclide": "Tl-208",
                "tier": "confirmed",
                "characteristic_line_keV": 2614.5,
                "specific_activity_Bq_per_kg": 1900.0,
                "specific_activity_sigma_Bq_per_kg": 60.0,
            },
        ]
    }), encoding="utf-8")
    out = mod.compute(report, A_cert=1940.0, rel_cert=0.06, cert_name="test")
    by_nuclide = {r["nuclide"]: r for r in out["nuclides"]}
    assert by_nuclide["Bi-212"]["verdict"] == "n/a"
    assert by_nuclide["Bi-212"]["verdict"] != "upper-limit"
    assert by_nuclide["Tl-208"]["verdict"] == "PASS"


def test_f4_summary_excludes_na_from_confirmed(tmp_path: Path):
    mod = _load_cert_zcheck()
    report = tmp_path / "fake_report.json"
    report.write_text(json.dumps({
        "identified_nuclides": [
            {
                "nuclide": "Bi-212",
                "tier": "tentative",
                "characteristic_line_keV": 727.3,
                "specific_activity_Bq_per_kg": None,
                "specific_activity_sigma_Bq_per_kg": None,
            },
        ]
    }), encoding="utf-8")
    out = mod.compute(report, A_cert=1940.0, rel_cert=0.06, cert_name="test")
    assert out["summary"]["n_confirmed"] == 0
    assert out["summary"]["all_pass"] is False