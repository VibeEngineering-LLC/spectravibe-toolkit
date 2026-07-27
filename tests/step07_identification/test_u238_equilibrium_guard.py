"""F-119 sub (v1.17.6) — chain-equilibrium guard для цепочки U-238.

В v1.17.5 ``chain_equilibrium_guard`` уже поддерживает U-238 в реестре
``CHAIN_MEMBERS`` (Pb-214 / Bi-214 / Po-214 / Pb-210 / Ra-226 / Th-234).
Однако регрессионного покрытия не было — этот тест закрепляет контракт
для v1.17.6.
"""
from __future__ import annotations

import math
import os
import sys

import pytest

SCRIPTS = os.path.join(os.path.dirname(__file__), "scripts")
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)

from gamma.activity.compute import (
    chain_equilibrium_guard, ActivityResult, CHAIN_MEMBERS,
)


def _ar(nuc: str, A: float) -> ActivityResult:
    """Minimal valid ActivityResult helper."""
    return ActivityResult(
        nuclide=nuc, A_Bq=float(A), sigma_A_Bq=float(A) * 0.05,
        lines_used=(),
    )


def test_chain_members_includes_u238_chain():
    """CHAIN_MEMBERS содержит U-238 цепочку с правильными нуклидами."""
    assert "U-238" in CHAIN_MEMBERS
    u_chain = CHAIN_MEMBERS["U-238"]
    # Ключевые нуклиды U-238 цепочки на NaI / HPGe
    for nuc in ("Bi-214", "Pb-214", "Pb-210", "Ra-226", "Po-214"):
        assert nuc in u_chain, f"{nuc} отсутствует в U-238 chain"


def test_u238_equilibrium_no_outliers_when_balanced():
    """При сбалансированной цепочке U-238 (ratio ≤ 5) — нет outlier-ов."""
    activities = [
        _ar("Bi-214", 100.0),
        _ar("Pb-214", 120.0),
        _ar("Pb-210", 80.0),
        _ar("Ra-226", 110.0),
    ]
    report = chain_equilibrium_guard(
        activities, chains=("U-238",), ratio_threshold=5.0,
    )
    assert "U-238" in report
    block = report["U-238"]
    assert block["in_equilibrium"]
    assert block["outliers"] == []
    assert block["ratio"] < 5.0


def test_u238_equilibrium_flags_outliers_when_ratio_exceeds_threshold():
    """При расхождении > 5× в цепочке U-238 — outlier помечается."""
    activities = [
        _ar("Bi-214", 100.0),
        _ar("Pb-214", 110.0),
        _ar("Pb-210", 5000.0),    # outlier: ×50 выше медианы
        _ar("Ra-226", 105.0),
    ]
    report = chain_equilibrium_guard(
        activities, chains=("U-238",), ratio_threshold=5.0,
    )
    block = report["U-238"]
    assert not block["in_equilibrium"]
    assert "Pb-210" in block["outliers"]
    # ratio = 5000 / 100 = 50
    assert block["ratio"] > 5.0


def test_chain_equilibrium_guard_runs_both_chains_in_one_call():
    """Один вызов покрывает обе цепочки — нет regression на Th-232 случай."""
    activities = [
        _ar("Tl-208", 50.0), _ar("Ac-228", 60.0), _ar("Pb-212", 55.0),
        _ar("Bi-214", 100.0), _ar("Pb-214", 110.0),
    ]
    report = chain_equilibrium_guard(activities)
    assert "Th-232" in report
    assert "U-238" in report
    assert report["Th-232"]["in_equilibrium"]
    assert report["U-238"]["in_equilibrium"]
