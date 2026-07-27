"""F-125 (v1.17.6) — Refit NaI 63×63 FWHM model regression.

Закрепляет новые коэффициенты ``_DEFAULT_NAI_FWHM_MODEL``, полученные
по 26 анкер-точкам из калибровочных спектров Gamma-1S. Проверяет:
  1. Модель монотонно возрастает на [50, 3000] кэВ.
  2. Не уходит в NaN / отрицательные значения нигде в диапазоне.
  3. Совпадает с измеренными FWHM на ключевых линиях с rms ≤ 5 кэВ.
  4. На E=661.66 даёт FWHM = 47 ± 2 кэВ (R = 7.1 ± 0.3 %).
"""
from __future__ import annotations

import math
import os
import sys

import pytest

SCRIPTS = os.path.join(os.path.dirname(__file__), "scripts")
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)

from gamma.identification.staged_pipeline import (
    _DEFAULT_NAI_FWHM_MODEL, fwhm_keV_at_energy,
)


# Анкер-точки из реальных калибровочных спектров (E_keV, FWHM_keV)
ANCHOR_POINTS = [
    (84.373, 11.865), (121.782, 15.438), (238.632, 23.507),
    (240.986, 23.654), (244.697, 23.8), (277.371, 25.903),
    (300.087, 27.259), (344.279, 29.835), (443.961, 35.152),
    (510.77, 38.666), (583.187, 42.157), (661.657, 45.739),
    (727.33, 48.699), (763.13, 50.241), (778.905, 50.748),
    (785.37, 51.183), (860.557, 54.3), (964.057, 58.561),
    (1078.62, 63.504), (1085.837, 63.163), (1112.076, 64.123),
    (1173.228, 66.018), (1332.492, 71.53), (1408.013, 74.315),
    (1620.5, 80.391), (2614.511, 106.985),
]


def test_default_nai_fwhm_model_is_v1_17_6():
    """Коэффициенты обновлены в v1.17.6 (F-125)."""
    a, b, c = _DEFAULT_NAI_FWHM_MODEL
    # Refit anchor: a = 0 (constrained), b ≈ 2.95, c ≈ 5.76e-4
    assert a == pytest.approx(0.0, abs=1e-3)
    assert b == pytest.approx(2.95, abs=0.5)
    assert c == pytest.approx(5.76e-4, abs=2e-4)


def test_fwhm_monotonic_on_full_range():
    """FWHM(E) монотонно неубывает на [50, 3000] кэВ."""
    energies = list(range(50, 3001, 10))
    fwhms = [fwhm_keV_at_energy(_DEFAULT_NAI_FWHM_MODEL, e)
             for e in energies]
    for i in range(1, len(fwhms)):
        assert fwhms[i] >= fwhms[i - 1] - 1e-3, (
            f"FWHM not monotonic at E={energies[i]}: "
            f"{fwhms[i-1]:.3f} → {fwhms[i]:.3f}"
        )


def test_fwhm_no_nan_anywhere():
    """FWHM(E) определён и > 0 на [10, 5000] кэВ (граничные значения)."""
    for e in [10, 50, 100, 661, 1460, 2614, 3000, 5000]:
        f = fwhm_keV_at_energy(_DEFAULT_NAI_FWHM_MODEL, e)
        assert math.isfinite(f), f"FWHM NaN at E={e}"
        assert f > 0, f"FWHM={f} ≤ 0 at E={e}"


def test_fwhm_matches_anchor_points_within_rms_5kev():
    """rms между моделью и анкер-точками ≤ 5 кэВ."""
    residuals = []
    for e, f_measured in ANCHOR_POINTS:
        f_model = fwhm_keV_at_energy(_DEFAULT_NAI_FWHM_MODEL, e)
        residuals.append(f_model - f_measured)
    n = len(residuals)
    rms = math.sqrt(sum(r * r for r in residuals) / n)
    assert rms < 5.0, f"rms={rms:.3f} ≥ 5 keV"


def test_fwhm_at_cs137_within_tight_band():
    """E=661.66 кэВ → FWHM в полосе 45-50 кэВ (R=6.8-7.6%)."""
    f = fwhm_keV_at_energy(_DEFAULT_NAI_FWHM_MODEL, 661.66)
    assert 45.0 <= f <= 50.0, f"FWHM(661.66)={f:.2f} вне полосы 45-50 кэВ"
    R_pct = 100.0 * f / 661.66
    assert 6.8 <= R_pct <= 7.6, f"R={R_pct:.2f}% вне полосы 6.8-7.6%"


def test_fwhm_at_low_energy_nonzero():
    """E=80 кэВ → FWHM > 10 кэВ (не уходит в 1.0 floor как в старой модели)."""
    f = fwhm_keV_at_energy(_DEFAULT_NAI_FWHM_MODEL, 80.0)
    assert f > 10.0, (
        f"FWHM(80)={f:.2f} ≤ 10 кэВ — модель проваливается на низких "
        f"энергиях (бывшая проблема при v1.17.5)"
    )


def test_fwhm_at_tl208_2614_within_band():
    """E=2614.51 кэВ → FWHM в полосе 100-115 кэВ."""
    f = fwhm_keV_at_energy(_DEFAULT_NAI_FWHM_MODEL, 2614.51)
    assert 100.0 <= f <= 115.0, f"FWHM(2614.51)={f:.2f} вне 100-115 кэВ"
