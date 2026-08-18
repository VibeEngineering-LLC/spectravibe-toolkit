# -*- coding: utf-8 -*-
"""F-329 / v1.18.18.3 (ROADMAP TD-2) — Compton-residual classification
для Co-60 и Na-22 (расширение F-143 whitelist).

Контракт:
  - Когда Co-60 в `detected_nuclides`, любой peak в 200-1200 кэВ →
    LBL_CHAIN_SECONDARY (compton_residual, parent=Co-60).
  - Когда Na-22 в `detected_nuclides`, peak в 150-400 кэВ →
    LBL_CHAIN_SECONDARY (compton_residual, parent=Na-22).
  - Эти rules дополняют существующие F-143 для Cs-137 (100-400) и
    K-40 (200-1300).

Цель: устранить FALSE-positive matches с U-235 / Ra-226 / Bi-214 в
spectrum'ах со strong calibration-source signal.
"""
from __future__ import annotations
import sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from gamma.identification.residual_classifier import (
    classify_residual,
    LBL_CHAIN_SECONDARY,
)


def _stub_classify(E_keV: float, detected: list[str]):
    """Wrap classify_residual с минимальными required-параметрами."""
    return classify_residual(
        peak_E_keV=E_keV,
        sigma=8.0,                       # arbitrary above threshold
        detected_nuclide_names=detected,
        fwhm_at_keV=15.0,                # NaI default ~15 keV @ 600 keV
        energy_max_keV=3000.0,
    )


# ─── Co-60 Compton residual (200-1200 keV) ──────────────────────────

@pytest.mark.parametrize("E_keV", [250, 500, 750, 1000, 1150])
def test_F329_co60_compton_residual_classified(E_keV):
    """Peak в 200-1200 кэВ при наличии Co-60 → compton_residual."""
    r = _stub_classify(E_keV, ["Co-60"])
    assert r.label == LBL_CHAIN_SECONDARY, (
        f"E={E_keV}: label={r.label} (ожидался chain_secondary)"
    )
    assert r.parent_nuclide == "Co-60"
    assert r.feature_kind == "compton_residual"


@pytest.mark.parametrize("E_keV", [50, 150, 199.9, 1200.1, 1500])
def test_F329_co60_outside_window_not_classified(E_keV):
    """Peak ВНЕ 200-1200 кэВ → НЕ classified как Co-60 residual."""
    r = _stub_classify(E_keV, ["Co-60"])
    # Может быть LBL_TRUE_UNMATCHED или LBL_EDGE_OF_RANGE,
    # но НЕ compton_residual от Co-60
    assert not (
        r.label == LBL_CHAIN_SECONDARY
        and r.parent_nuclide == "Co-60"
    )


def test_F329_co60_not_in_detected_no_classification():
    """Без Co-60 в detected_nuclides — peak не attributed Co-60."""
    r = _stub_classify(500.0, ["Cs-137"])
    assert r.parent_nuclide != "Co-60"


# ─── Na-22 Compton residual (150-400 keV) ───────────────────────────

# Note: 252.5 keV исключён — это double_escape от Na-22 1274.5 keV
# (1274.5-1022=252.5), которая срабатывает раньше F-329 в order of checks
# (correct precedence: более specific feature > generic compton).
@pytest.mark.parametrize("E_keV", [160, 200, 350, 399])
def test_F329_na22_compton_residual_classified(E_keV):
    """Peak в 150-400 кэВ при наличии Na-22 → compton_residual."""
    r = _stub_classify(E_keV, ["Na-22"])
    assert r.label == LBL_CHAIN_SECONDARY
    assert r.parent_nuclide == "Na-22"
    assert r.feature_kind == "compton_residual"


@pytest.mark.parametrize("E_keV", [100, 149.9, 400.1, 600])
def test_F329_na22_outside_window_not_classified(E_keV):
    """Peak ВНЕ 150-400 кэВ → НЕ classified как Na-22 residual."""
    r = _stub_classify(E_keV, ["Na-22"])
    assert not (
        r.label == LBL_CHAIN_SECONDARY
        and r.parent_nuclide == "Na-22"
    )


# ─── Backward compat: Cs-137 / K-40 rules должны продолжать работать ──

def test_F329_cs137_still_classified_in_100_400():
    """F-143 baseline: Cs-137 100-400 кэВ residual."""
    r = _stub_classify(250.0, ["Cs-137"])
    assert r.label == LBL_CHAIN_SECONDARY
    assert r.parent_nuclide == "Cs-137"


def test_F329_k40_still_classified_in_200_1300():
    """F-143 baseline: K-40 200-1300 кэВ residual."""
    r = _stub_classify(800.0, ["K-40"])
    assert r.label == LBL_CHAIN_SECONDARY
    assert r.parent_nuclide == "K-40"
