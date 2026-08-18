"""G3 / v1.31.2 -- continuum-level diagnostic unit tests.

Exercises continuum_diagnostic() added to environment.py:
* low-cps spectrum -> hint "low_bg"
* mid-cps spectrum -> hint "intermediate"/"natural"
* high-cps spectrum -> hint "high"
* ERN-line dominance computation
* missing spec / live_time -> graceful None
"""
from __future__ import annotations

import sys
import pathlib

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from gamma.reporting.environment import (
    continuum_diagnostic,
    ERN_NUCLIDES_FOR_DOMINANCE,
    CPS_HINT_LOW_BG_MAX,
    CPS_HINT_INTERMEDIATE_MAX,
    CPS_HINT_NATURAL_MAX,
)


class _Spec:
    def __init__(self, counts, live_time):
        self.counts = counts
        self.live_time = live_time


class _LM:
    def __init__(self, peak_area):
        self.peak_area = peak_area


class _NI:
    def __init__(self, nuclide, lines):
        self.nuclide = nuclide
        self.matched_lines = tuple(lines)


class _Result:
    def __init__(self, spec=None, detected=()):
        self.spec = spec
        self.final_detected = tuple(detected)


def test_low_cps_hint():
    # 100 counts over 200 s -> 0.5 cps -> low_bg
    spec = _Spec(counts=[100], live_time=200.0)
    res = _Result(spec=spec)
    d = continuum_diagnostic(res)
    assert d["total_cps"] == 0.5
    assert d["environment_hint_by_cps"] == "low_bg"
    assert d["bg_line_dominance_pct"] == 0.0


def test_intermediate_hint():
    # 1000 counts / 200 s -> 5 cps -> intermediate
    spec = _Spec(counts=[1000], live_time=200.0)
    res = _Result(spec=spec)
    d = continuum_diagnostic(res)
    assert d["total_cps"] == 5.0
    assert d["environment_hint_by_cps"] == "intermediate"


def test_natural_hint():
    # 10000 counts / 200 s -> 50 cps -> natural
    spec = _Spec(counts=[10000], live_time=200.0)
    res = _Result(spec=spec)
    d = continuum_diagnostic(res)
    assert d["total_cps"] == 50.0
    assert d["environment_hint_by_cps"] == "natural"


def test_high_hint():
    # 100000 counts / 200 s -> 500 cps -> high
    spec = _Spec(counts=[100000], live_time=200.0)
    res = _Result(spec=spec)
    d = continuum_diagnostic(res)
    assert d["total_cps"] == 500.0
    assert d["environment_hint_by_cps"] == "high"


def test_ern_line_dominance_computed():
    # 1000 total counts, 200 cps total from K-40+Tl-208
    spec = _Spec(counts=[1000], live_time=10.0)
    k40 = _NI("K-40", [_LM(150.0), _LM(50.0)])
    tl208 = _NI("Tl-208", [_LM(100.0)])
    # Non-ERN nuclide should be skipped
    cs137 = _NI("Cs-137", [_LM(500.0)])
    res = _Result(spec=spec, detected=[k40, tl208, cs137])
    d = continuum_diagnostic(res)
    assert d["ern_line_area_sum"] == 300.0
    # 300 / 1000 = 30%
    assert abs(d["bg_line_dominance_pct"] - 30.0) < 1e-9


def test_missing_spec_returns_none_values():
    res = _Result(spec=None)
    d = continuum_diagnostic(res)
    assert d["total_cps"] is None
    assert d["bg_line_dominance_pct"] is None
    assert d["environment_hint_by_cps"] is None


def test_zero_live_time_returns_none():
    spec = _Spec(counts=[100], live_time=0.0)
    res = _Result(spec=spec)
    d = continuum_diagnostic(res)
    assert d["total_cps"] is None
    assert d["environment_hint_by_cps"] is None


def test_thresholds_locked():
    assert CPS_HINT_LOW_BG_MAX == 1.0
    assert CPS_HINT_INTERMEDIATE_MAX == 10.0
    assert CPS_HINT_NATURAL_MAX == 100.0


def test_ern_set_contains_expected_nuclides():
    for nuc in ("K-40", "Tl-208", "Bi-214", "Ac-228", "Pb-212", "Pb-214"):
        assert nuc in ERN_NUCLIDES_FOR_DOMINANCE


def test_diagnostic_never_raises_on_garbage():
    # Verify graceful path with object that lacks attributes
    class Bare: pass
    res = Bare()
    d = continuum_diagnostic(res)
    assert d["total_cps"] is None
    assert isinstance(d, dict)