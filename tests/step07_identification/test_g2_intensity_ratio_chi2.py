"""G2 / v1.31.2 -- intensity-ratio chi^2 gate unit tests.

Exercises the annotation-only gate added to cross_check.py:
* coherent lines -> small chi^2, "pass_strict"
* divergent lines -> large chi^2, "fail"
* mid-range ratio -> "pass_lenient"
* per-nuclide skip when too few usable lines
* fallback area/significance estimator when peak_area_uncertainty missing
"""
from __future__ import annotations

import sys
import pathlib

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from gamma.identification.cross_check import (
    intensity_ratio_chi2_gate,
    INTENSITY_RATIO_CHI2_STRICT_THRESHOLD,
    INTENSITY_RATIO_CHI2_LENIENT_THRESHOLD,
)


class _LM:
    """Mimic LineMatch with the fields the gate reads."""

    def __init__(
        self,
        library_E_keV,
        library_I_pct,
        peak_area,
        peak_area_uncertainty=None,
        significance_currie=None,
    ):
        self.library_E_keV = library_E_keV
        self.library_I_pct = library_I_pct
        self.peak_area = peak_area
        self.peak_area_uncertainty = peak_area_uncertainty
        self.significance_currie = significance_currie


class _NI:
    def __init__(self, nuclide, matched_lines):
        self.nuclide = nuclide
        self.matched_lines = tuple(matched_lines)


class _IR:
    def __init__(self, detected):
        self.detected_nuclides = tuple(detected)


def test_coherent_lines_pass_strict():
    # Q = area / I should be ~constant across lines
    lines = [
        _LM(338.3, 11.27, 1127.0, peak_area_uncertainty=33.0),
        _LM(911.2, 25.8, 2580.0, peak_area_uncertainty=70.0),
        _LM(968.97, 15.8, 1580.0, peak_area_uncertainty=45.0),
    ]
    ident = _IR([_NI("Ac-228", lines)])
    gate = intensity_ratio_chi2_gate(ident)
    assert "Ac-228" in gate
    res = gate["Ac-228"]
    assert res["n_lines_used"] == 3
    assert res["ndof"] == 2
    # Q should all be ~100
    assert abs(res["q_mean"] - 100.0) < 1.0
    assert res["ratio"] <= INTENSITY_RATIO_CHI2_STRICT_THRESHOLD
    assert res["verdict"] == "pass_strict"


def test_divergent_lines_fail():
    # Line areas drift hard from the I-implied prediction
    lines = [
        _LM(338.3, 11.27, 1127.0, peak_area_uncertainty=10.0),  # Q=100
        _LM(911.2, 25.8, 5160.0, peak_area_uncertainty=10.0),   # Q=200
        _LM(968.97, 15.8, 4740.0, peak_area_uncertainty=10.0),  # Q=300
    ]
    ident = _IR([_NI("Ac-228", lines)])
    gate = intensity_ratio_chi2_gate(ident)
    res = gate["Ac-228"]
    assert res["ratio"] > INTENSITY_RATIO_CHI2_LENIENT_THRESHOLD
    assert res["verdict"] == "fail"


def test_skip_when_only_one_usable_line():
    # Only one line with intensity above the floor
    lines = [_LM(338.3, 11.27, 1127.0, peak_area_uncertainty=33.0)]
    ident = _IR([_NI("Ac-228", lines)])
    gate = intensity_ratio_chi2_gate(ident)
    assert "Ac-228" not in gate


def test_fallback_uncertainty_from_significance():
    # Significance present, uncertainty missing -> derive sigma = area / sig
    lines = [
        _LM(338.3, 11.27, 1127.0, peak_area_uncertainty=None, significance_currie=34.0),
        _LM(911.2, 25.8, 2580.0, peak_area_uncertainty=None, significance_currie=37.0),
    ]
    ident = _IR([_NI("Ac-228", lines)])
    gate = intensity_ratio_chi2_gate(ident)
    assert "Ac-228" in gate
    res = gate["Ac-228"]
    assert res["n_lines_used"] == 2
    # Coherent areas -> low ratio
    assert res["ratio"] < INTENSITY_RATIO_CHI2_LENIENT_THRESHOLD


def test_low_intensity_lines_skipped():
    # 0.5% intensity is below the 1% default floor
    lines = [
        _LM(338.3, 11.27, 1127.0, peak_area_uncertainty=33.0),
        _LM(911.2, 25.8, 2580.0, peak_area_uncertainty=70.0),
        _LM(99.0, 0.5, 50.0, peak_area_uncertainty=10.0),
    ]
    ident = _IR([_NI("Ac-228", lines)])
    gate = intensity_ratio_chi2_gate(ident)
    res = gate["Ac-228"]
    assert res["n_lines_used"] == 2
    assert res["n_lines_skipped"] == 1


def test_thresholds_locked():
    # Audit-locked numeric values
    assert INTENSITY_RATIO_CHI2_STRICT_THRESHOLD == 1.5
    assert INTENSITY_RATIO_CHI2_LENIENT_THRESHOLD == 3.0


def test_empty_identification_returns_empty_dict():
    assert intensity_ratio_chi2_gate(None) == {}
    assert intensity_ratio_chi2_gate(_IR([])) == {}


def test_zero_area_lines_skipped():
    lines = [
        _LM(338.3, 11.27, 0.0, peak_area_uncertainty=10.0),
        _LM(911.2, 25.8, 2580.0, peak_area_uncertainty=70.0),
    ]
    ident = _IR([_NI("Ac-228", lines)])
    gate = intensity_ratio_chi2_gate(ident)
    # Only one usable line -> nuclide omitted
    assert "Ac-228" not in gate