"""F-440 / v1.30.0 — unit-tests для two-phase weak-line completion.

Проверяет:
  * complete_weak_lines: разделение на phase1_fitted vs phase2_completed
  * Формула S_expected = A * I/100 * eps * t_live * f_self_abs * f_TCS
  * completeness_pct = 100 * S_fitted / (S_fitted + S_completed)
  * Phantom lines (peak_area_source содержит "phantom") → всегда weak
  * Intensity gate: I_pct < min_grouping_intensity_pct → weak
  * Contamination tracking через multiplet ROI
  * to_json_block: структура JSON-блока
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.activity.weak_line_completion import (  # type: ignore
    DEFAULT_MIN_GROUPING_INTENSITY_PCT,
    DEFAULT_MIN_GROUPING_SNR,
    CompletedLine,
    FittedLineSummary,
    NuclideCompletion,
    WeakContamination,
    complete_weak_lines,
    to_json_block,
)


@dataclass
class FakeActivityResult:
    nuclide: str
    A_Bq: float
    sigma_A_Bq: float


@dataclass
class FakeLineMatch:
    library_E_keV: float
    library_I_pct: float
    peak_area: Optional[float] = None
    peak_area_uncertainty: Optional[float] = None
    peak_area_source: str = ""
    significance_currie: Optional[float] = None


class FakeEfficiencyCurve:
    """Stub-efficiency: eps(E) = 0.05 на всех E, никогда не extrapolated."""

    def efficiency_at(self, E_keV: float) -> float:
        return 0.05

    def is_extrapolating(self, E_keV: float) -> bool:
        return False


def _tl208_library() -> dict:
    return {
        "Tl-208": {
            "lines": [
                # E_keV, I_pct, sigma_I_pct_abs
                (583.19, 30.55, 0.30),   # strong
                (2614.51, 35.85, 0.36),  # strong
                (277.37, 2.37, 0.05),    # weak (I < 3%)
                (763.13, 0.64, 0.02),    # very weak
                (860.56, 4.49, 0.10),    # mid (above 3% but no measured area in matches)
            ]
        }
    }


def _make_matches_tl208():
    """3 matches: 583 strong (cowell), 2614 strong, 277 phantom (Phase 1 gate)."""
    return [
        FakeLineMatch(
            library_E_keV=583.19, library_I_pct=30.55,
            peak_area=2.0e5, peak_area_uncertainty=2.0e3,
            peak_area_source="cowell",
        ),
        FakeLineMatch(
            library_E_keV=2614.51, library_I_pct=35.85,
            peak_area=5.0e4, peak_area_uncertainty=5.0e2,
            peak_area_source="cowell",
        ),
        FakeLineMatch(
            library_E_keV=277.37, library_I_pct=2.37,
            peak_area=None, peak_area_uncertainty=None,
            peak_area_source="library_anchor_phantom",
        ),
    ]


# ─────────────────────────────────────────────────────────────────
# Tests
# ─────────────────────────────────────────────────────────────────


def test_phase1_phase2_split_basic():
    """Strong lines с реальной площадью → phase1; phantom → phase2 (вместе с
    нефитированными library-линиями)."""
    matches = _make_matches_tl208()
    activities = [FakeActivityResult("Tl-208", A_Bq=3000.0, sigma_A_Bq=100.0)]
    out = complete_weak_lines(
        activities,
        matched_lines_by_nuclide={"Tl-208": matches},
        nuclide_library=_tl208_library(),
        efficiency_curve=FakeEfficiencyCurve(),
        t_live=1000.0,
        sample_mass_kg=1.6,
        min_grouping_snr=5.0,
        min_grouping_intensity_pct=3.0,
    )
    assert "Tl-208" in out
    c = out["Tl-208"]
    # Phase 1: только 583, 2614 (cowell c real area, I>=3)
    e_strong = sorted(round(f.E_keV, 1) for f in c.phase1_fitted_lines)
    assert e_strong == [583.2, 2614.5], f"phase1 mismatch: {e_strong}"
    # Phase 2: 277 (phantom; в библиотеке + phantom-match), 763, 860 — не в strong-set
    e_compl = sorted(round(cl.E_keV, 1) for cl in c.phase2_completed_lines)
    assert e_compl == [277.4, 763.1, 860.6], f"phase2 mismatch: {e_compl}"


def test_s_expected_formula():
    """S_expected должен соответствовать A * I/100 * eps * t_live."""
    matches = _make_matches_tl208()
    A = 2000.0
    activities = [FakeActivityResult("Tl-208", A_Bq=A, sigma_A_Bq=80.0)]
    out = complete_weak_lines(
        activities,
        matched_lines_by_nuclide={"Tl-208": matches},
        nuclide_library=_tl208_library(),
        efficiency_curve=FakeEfficiencyCurve(),
        t_live=500.0,
        sample_mass_kg=1.6,
    )
    c = out["Tl-208"]
    # Find the 763.13 line in phase2 (I_pct=0.64)
    line_763 = next((cl for cl in c.phase2_completed_lines if round(cl.E_keV, 2) == 763.13), None)
    assert line_763 is not None
    expected = A * (0.64 / 100.0) * 0.05 * 500.0 * 1.0 * 1.0  # f_self=f_tcs=1
    assert math.isclose(line_763.S_expected, expected, rel_tol=1e-6), \
        f"S_expected mismatch: got {line_763.S_expected}, expected {expected}"


def test_completeness_pct_formula():
    """completeness_pct = 100 * fitted / (fitted + completed)."""
    matches = _make_matches_tl208()
    activities = [FakeActivityResult("Tl-208", A_Bq=3000.0, sigma_A_Bq=100.0)]
    out = complete_weak_lines(
        activities,
        matched_lines_by_nuclide={"Tl-208": matches},
        nuclide_library=_tl208_library(),
        efficiency_curve=FakeEfficiencyCurve(),
        t_live=1000.0,
        sample_mass_kg=1.6,
    )
    c = out["Tl-208"]
    fitted_total = sum(f.S_measured for f in c.phase1_fitted_lines if f.S_measured is not None)
    completed_total = sum(cl.S_expected for cl in c.phase2_completed_lines)
    expected_pct = 100.0 * fitted_total / (fitted_total + completed_total)
    assert math.isclose(c.completeness_pct, expected_pct, rel_tol=1e-6)
    assert math.isclose(c.fitted_area_total, fitted_total, rel_tol=1e-6)
    assert math.isclose(c.completed_area_total, completed_total, rel_tol=1e-6)


def test_phantom_line_is_weak():
    """Линия с peak_area_source='library_anchor_phantom' всегда → weak,
    даже если I_pct высокий."""
    fake = FakeLineMatch(
        library_E_keV=583.19, library_I_pct=30.55,
        peak_area=1.0e5, peak_area_uncertainty=1.0e3,
        peak_area_source="library_anchor_phantom",
    )
    activities = [FakeActivityResult("Tl-208", A_Bq=3000.0, sigma_A_Bq=100.0)]
    out = complete_weak_lines(
        activities,
        matched_lines_by_nuclide={"Tl-208": [fake]},
        nuclide_library=_tl208_library(),
        efficiency_curve=FakeEfficiencyCurve(),
        t_live=1000.0,
        sample_mass_kg=1.6,
        min_grouping_intensity_pct=3.0,
    )
    c = out["Tl-208"]
    assert len(c.phase1_fitted_lines) == 0, "phantom-line должен быть weak, не fitted"


def test_intensity_gate_skips_low_I():
    """Линия с I_pct < min_grouping_intensity_pct и реальной площадью → weak."""
    fake = FakeLineMatch(
        library_E_keV=277.37, library_I_pct=2.37,
        peak_area=1.0e4, peak_area_uncertainty=1.0e3,  # snr=10 (above SNR gate)
        peak_area_source="cowell",
    )
    activities = [FakeActivityResult("Tl-208", A_Bq=3000.0, sigma_A_Bq=100.0)]
    out = complete_weak_lines(
        activities,
        matched_lines_by_nuclide={"Tl-208": [fake]},
        nuclide_library=_tl208_library(),
        efficiency_curve=FakeEfficiencyCurve(),
        t_live=1000.0,
        sample_mass_kg=1.6,
        min_grouping_snr=5.0,
        min_grouping_intensity_pct=3.0,
    )
    c = out["Tl-208"]
    assert len(c.phase1_fitted_lines) == 0, \
        "Line с I_pct=2.37<3.0 должна быть weak независимо от SNR"


def test_no_activity_returns_no_completed_lines():
    """A_Bq=None → phase2 list пуст (нечем считать S_expected)."""
    matches = _make_matches_tl208()
    activities = []  # no activity for Tl-208
    out = complete_weak_lines(
        activities,
        matched_lines_by_nuclide={"Tl-208": matches},
        nuclide_library=_tl208_library(),
        efficiency_curve=FakeEfficiencyCurve(),
        t_live=1000.0,
        sample_mass_kg=1.6,
    )
    c = out["Tl-208"]
    assert len(c.phase2_completed_lines) == 0
    # completeness=100% когда нет completed (бесполезно, но не crash)
    assert c.completeness_pct == 100.0 or len(c.phase1_fitted_lines) == 0


def test_contamination_into_strong_peaks():
    """Tl-208 contaminates Ac-228 в M1 multiplet (763 → ROI 700-800)."""
    tl208_matches = _make_matches_tl208()
    ac228_matches = [
        FakeLineMatch(
            library_E_keV=911.20, library_I_pct=25.8,
            peak_area=1.0e5, peak_area_uncertainty=1.0e3,
            peak_area_source="cowell",
        ),
    ]
    ac228_lib = {"Ac-228": {"lines": [(911.20, 25.8, 0.3)]}}
    library = {**_tl208_library(), **ac228_lib}
    activities = [
        FakeActivityResult("Tl-208", A_Bq=3000.0, sigma_A_Bq=100.0),
        FakeActivityResult("Ac-228", A_Bq=2000.0, sigma_A_Bq=80.0),
    ]
    multiplet_rois = [
        {"label": "M1", "E_lo_keV": 700.0, "E_hi_keV": 920.0},
    ]
    out = complete_weak_lines(
        activities,
        matched_lines_by_nuclide={
            "Tl-208": tl208_matches,
            "Ac-228": ac228_matches,
        },
        nuclide_library=library,
        efficiency_curve=FakeEfficiencyCurve(),
        t_live=1000.0,
        sample_mass_kg=1.6,
        multiplet_rois=multiplet_rois,
        contamination_threshold_pct=0.5,
    )
    c_tl = out["Tl-208"]
    # Tl-208's weak lines 763 + 860 попадают в M1 (700-920 keV)
    contam = [w for w in c_tl.weak_contamination_into_strong_peaks
              if w.strong_nuclide == "Ac-228" and w.multiplet_label == "M1"]
    assert len(contam) == 1, f"Expected contamination Tl→Ac M1, got {contam}"
    assert contam[0].S_contamination > 0
    assert contam[0].fraction_of_strong_pct > 0


def test_to_json_block_structure():
    """JSON-блок содержит обязательные ключи на per-nuclide уровне."""
    matches = _make_matches_tl208()
    activities = [FakeActivityResult("Tl-208", A_Bq=3000.0, sigma_A_Bq=100.0)]
    completions = complete_weak_lines(
        activities,
        matched_lines_by_nuclide={"Tl-208": matches},
        nuclide_library=_tl208_library(),
        efficiency_curve=FakeEfficiencyCurve(),
        t_live=1000.0,
        sample_mass_kg=1.6,
    )
    block = to_json_block(completions)
    assert "Tl-208" in block
    nuc_block = block["Tl-208"]
    for key in (
        "phase1_fitted_lines", "phase1_activity_Bq",
        "phase1_activity_Bq_kg", "phase2_completed_lines",
        "fitted_area_total", "completed_area_total",
        "completeness_pct", "weak_contamination_into_strong_peaks",
    ):
        assert key in nuc_block, f"missing key {key} in JSON block"
    # All numeric fields can be safely serialised
    import json
    s = json.dumps(block, ensure_ascii=False)
    assert "Tl-208" in s


def test_defaults_match_brief():
    """F-440 brief defaults: SNR=5.0, I_gamma=3.0%."""
    assert DEFAULT_MIN_GROUPING_SNR == 5.0
    assert DEFAULT_MIN_GROUPING_INTENSITY_PCT == 3.0