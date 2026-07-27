"""
F-294 (v1.17.20, T-002 + T-024) — Marinelli self-absorption correction.

Когда **плотность образца** ρ_sample отличается от плотности эталона
ρ_calib, использованного в efficiency calibration, FEP-эффективность
ε(E) меняется из-за разного **самопоглощения** γ-квантов в матрице.

Это особенно критично на низких энергиях (E < 200 keV), где
сечение фотопоглощения резко растёт. Для Marinelli 0.5 L на NaI
63×63 эффект может достигать ±15 % при ρ 0.5 → 2.5 g/cm³.

Базовая формула (приближение «эквивалентный путь»)
--------------------------------------------------
Для одинаковой геометрии (Marinelli один объём) и матрицы,
близкой к water-equivalent, корректирующий фактор записывается как:

    f_abs(E) = ε_sample(E) / ε_calib(E)
             ≈ (1 - exp(-μ_sample(E)·x̄)) / (1 - exp(-μ_calib(E)·x̄))
                × exp(-(μ_sample - μ_calib)·d̄_outer)

где:

  • μ(E)   = ρ · (μ/ρ)(E)            — линейный коэф. поглощения, 1/cm
  • x̄     — средняя длина пути в образце (Marinelli ≈ 1.5–2.0 cm для 0.5 L)
  • d̄_outer — средняя толщина «крышки» между матрицей и детектором
              (для Marinelli стенки cup ≈ 0.1 cm полиэтилена → ничтожно)

Первый множитель — поправка на **внутреннее самопоглощение**
(больше образца → больше самопоглощение → меньше escape γ).
Второй — поправка на **прохождение от центра масс образца до детектора**
(почти 1 для Marinelli, т.к. cup тонкий).

Что НЕ покрыто этим модулем
---------------------------
- Для произвольной геометрии (vial, plate, soil deep-bed) нужна
  numerical integration по объёму с MC — выходит за scope T-002.
- Resonance K-edge effects (Pb, U) — требует full XCOM.
- Coherent scattering (Rayleigh) — отдельная корректировка ниже 100 keV.

Этот модуль даёт **first-order analytical** correction; для
production-grade рекомендуется EffCalcMC / GeantMC.

References
----------
- ЛСРМ Algorithmic Foundations 2022 § 9.2 «Самопоглощение в образце»
- Gilmore & Joss 3rd Ed. § 7.3 «Sample-related effects»
- NIST XCOM (mass attenuation coefficients) — μ/ρ tables
- Cutshall NH et al., NIM 206 (1983) 309 — direct-transmission method
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Sequence


# Marinelli 0.5 L geometry default (LSRM ОЭБ).
MARINELLI_05L_MEAN_PATH_CM = 1.75       # x̄ — average γ-path in active matrix
MARINELLI_05L_CUP_THICKNESS_CM = 0.10   # стенка cup (PE / PP)

# Reference μ/ρ table for water-equivalent matrix (low to high E).
# Источник: NIST XCOM, water Z_eff ≈ 7.42. Used as bootstrap fallback
# when user не передал явный mu_over_rho_func.
# Format: (E_keV, μ/ρ_cm2_per_g)
_WATER_MU_OVER_RHO_TABLE = [
    (30.0,  0.376),
    (40.0,  0.268),
    (50.0,  0.227),
    (60.0,  0.206),
    (80.0,  0.184),
    (100.0, 0.171),
    (150.0, 0.151),
    (200.0, 0.137),
    (300.0, 0.119),
    (400.0, 0.106),
    (500.0, 0.097),
    (600.0, 0.089),
    (800.0, 0.079),
    (1000.0, 0.071),
    (1250.0, 0.063),
    (1500.0, 0.057),
    (2000.0, 0.049),
    (3000.0, 0.040),
]


def mu_over_rho_water(E_keV: float) -> float:
    """Линейная интерполяция μ/ρ для water-equivalent matrix.

    Источник: NIST XCOM. Используется как fallback, когда user не
    передал tabulated coefficient для своей матрицы.

    Parameters
    ----------
    E_keV : float
        Энергия γ-кванта.

    Returns
    -------
    μ/ρ, cm² / g.
    """
    if E_keV <= _WATER_MU_OVER_RHO_TABLE[0][0]:
        return _WATER_MU_OVER_RHO_TABLE[0][1]
    if E_keV >= _WATER_MU_OVER_RHO_TABLE[-1][0]:
        return _WATER_MU_OVER_RHO_TABLE[-1][1]
    # Log-log interp (more accurate for power-law-like cross-section).
    for i in range(len(_WATER_MU_OVER_RHO_TABLE) - 1):
        E0, mu0 = _WATER_MU_OVER_RHO_TABLE[i]
        E1, mu1 = _WATER_MU_OVER_RHO_TABLE[i + 1]
        if E0 <= E_keV <= E1:
            ln_mu = (
                math.log(mu0)
                + (math.log(mu1) - math.log(mu0))
                * (math.log(E_keV) - math.log(E0))
                / (math.log(E1) - math.log(E0))
            )
            return math.exp(ln_mu)
    return _WATER_MU_OVER_RHO_TABLE[-1][1]


@dataclass(frozen=True)
class SelfAbsorptionInputs:
    """Параметры расчёта f_abs(E)."""

    E_keV: float
    rho_sample_g_cm3: float          # плотность измеряемого образца
    rho_calib_g_cm3: float           # плотность эталона калибровки
    mean_path_cm: float = MARINELLI_05L_MEAN_PATH_CM
    cup_thickness_cm: float = MARINELLI_05L_CUP_THICKNESS_CM
    mu_over_rho_sample: Optional[float] = None   # cm²/g; default = water
    mu_over_rho_calib: Optional[float] = None    # cm²/g; default = water


def self_absorption_factor(inputs: SelfAbsorptionInputs) -> float:
    r"""Compute self-absorption correction factor f_abs(E).

    Formula (Cutshall-type, 1st order analytic for Marinelli):

        f_abs = [(1 - exp(-μ_s·x̄)) / (μ_s·x̄)]
              / [(1 - exp(-μ_c·x̄)) / (μ_c·x̄)]
              × exp(-(μ_s - μ_c) · d̄_cup)

    where μ = ρ · (μ/ρ). Если плотности равны и μ/ρ одинаковы → f_abs = 1.

    Returns
    -------
    f_abs (decimal). Use as:
        A_corrected = A_apparent / f_abs
    or equivalently:
        ε_effective = ε_calib · f_abs

    Limitations
    -----------
    1st-order analytic for cylindrically-symmetric Marinelli. Не учитывает
    multiple scattering, build-up factor, edge effects. Для ρ ≥ 2.0 g/cm³
    погрешность может достигать 10 %. Для precision используйте MC.
    """
    if inputs.E_keV <= 0:
        raise ValueError("E_keV must be > 0")
    if inputs.rho_sample_g_cm3 <= 0 or inputs.rho_calib_g_cm3 <= 0:
        raise ValueError("densities must be > 0")
    if inputs.mean_path_cm <= 0:
        raise ValueError("mean_path_cm must be > 0")

    mu_rho_s = (
        inputs.mu_over_rho_sample
        if inputs.mu_over_rho_sample is not None
        else mu_over_rho_water(inputs.E_keV)
    )
    mu_rho_c = (
        inputs.mu_over_rho_calib
        if inputs.mu_over_rho_calib is not None
        else mu_over_rho_water(inputs.E_keV)
    )

    mu_s = inputs.rho_sample_g_cm3 * mu_rho_s
    mu_c = inputs.rho_calib_g_cm3 * mu_rho_c

    def _escape_fraction(mu: float, x: float) -> float:
        """(1 - e^{-μx}) / (μx) — escape probability normalized."""
        mux = mu * x
        if mux < 1e-9:
            return 1.0 - 0.5 * mux   # Taylor: 1 - μx/2 + (μx)²/6 - ...
        return (1.0 - math.exp(-mux)) / mux

    inner_ratio = (
        _escape_fraction(mu_s, inputs.mean_path_cm)
        / _escape_fraction(mu_c, inputs.mean_path_cm)
    )
    outer_correction = math.exp(
        -(mu_s - mu_c) * inputs.cup_thickness_cm
    )

    return inner_ratio * outer_correction


def correct_activity_for_self_absorption(
    A_apparent_Bq: float, f_abs: float,
) -> float:
    """Применить self-absorption correction к apparent activity.

        A_corrected = A_apparent / f_abs

    Если f_abs < 1 (больше самопоглощение чем эталон) →
    A_corrected > A_apparent (true activity ВЫШЕ).
    """
    if f_abs <= 0:
        raise ValueError(f"f_abs must be > 0, got {f_abs}")
    return A_apparent_Bq / f_abs


def batch_self_absorption_factors(
    energies_keV: Sequence[float],
    rho_sample: float,
    rho_calib: float,
    mean_path_cm: float = MARINELLI_05L_MEAN_PATH_CM,
    cup_thickness_cm: float = MARINELLI_05L_CUP_THICKNESS_CM,
) -> list[float]:
    """Бакетный расчёт f_abs для серии энергий (water matrix default)."""
    return [
        self_absorption_factor(
            SelfAbsorptionInputs(
                E_keV=E,
                rho_sample_g_cm3=rho_sample,
                rho_calib_g_cm3=rho_calib,
                mean_path_cm=mean_path_cm,
                cup_thickness_cm=cup_thickness_cm,
            )
        )
        for E in energies_keV
    ]


__all__ = [
    "MARINELLI_05L_MEAN_PATH_CM",
    "MARINELLI_05L_CUP_THICKNESS_CM",
    "mu_over_rho_water",
    "SelfAbsorptionInputs",
    "self_absorption_factor",
    "correct_activity_for_self_absorption",
    "batch_self_absorption_factors",
]
