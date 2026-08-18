"""Regression — F-119 chain-equilibrium guard.

Контракт: для Th-232 фикстуры все члены цепи (Tl-208, Pb-212,
Ac-228, Bi-212) должны быть в пределах фактора 5× друг от друга
(секулярное равновесие); ни один не помечается как outlier.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

from gamma.activity.compute import chain_equilibrium_guard


class _FakeActivity:
    def __init__(self, nuclide, A_Bq, is_valid_=True):
        self.nuclide = nuclide
        self.A_Bq = A_Bq
        self.sigma_A_Bq = A_Bq * 0.1
        self._valid = is_valid_

    def is_valid(self):
        return self._valid and self.A_Bq > 0


def test_chain_equilibrium_no_outliers_for_th232():
    # Модель: Bi-212 ≈ Ac-228 ≈ Tl-208 ≈ Pb-212 ≈ 1500 Бк/кг — все в
    # пределах фактора 2.
    acts = [
        _FakeActivity("Tl-208", 1500.0),
        _FakeActivity("Pb-212", 1700.0),
        _FakeActivity("Ac-228", 1900.0),
        _FakeActivity("Bi-212", 1300.0),
    ]
    result = chain_equilibrium_guard(acts)
    assert "Th-232" in result, f"Th-232 chain absent in result: {result}"
    th = result["Th-232"]
    assert th["in_equilibrium"] is True, f"unexpected outliers: {th}"
    assert th["outliers"] == [], f"unexpected outliers: {th['outliers']}"
    assert th["ratio"] < 5.0, f"ratio {th['ratio']:.2f} should be < 5"
    print(f"  ✓ test_chain_equilibrium_no_outliers_for_th232 "
          f"(ratio={th['ratio']:.2f}, median={th['median_Bq']:.0f} Bq)")


def test_chain_equilibrium_flags_outliers_when_ratio_exceeds_5x():
    # Модель: Ac-228 = 50 Бк/кг (выпад), остальные ≈ 1500 Бк/кг → ratio ≈ 30
    acts = [
        _FakeActivity("Tl-208", 1500.0),
        _FakeActivity("Pb-212", 1700.0),
        _FakeActivity("Ac-228", 50.0),
        _FakeActivity("Bi-212", 1300.0),
    ]
    result = chain_equilibrium_guard(acts)
    th = result["Th-232"]
    assert th["in_equilibrium"] is False
    assert "Ac-228" in th["outliers"], f"expected Ac-228 outlier: {th}"
    print(f"  ✓ test_chain_equilibrium_flags_outliers_when_ratio_exceeds_5x "
          f"(ratio={th['ratio']:.1f}, outliers={th['outliers']})")


if __name__ == "__main__":
    test_chain_equilibrium_no_outliers_for_th232()
    test_chain_equilibrium_flags_outliers_when_ratio_exceeds_5x()
    print("Chain-equilibrium guard regression PASS.")
