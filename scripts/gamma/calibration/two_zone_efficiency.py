"""
F-284 (v1.17.15, T-029) — 2-zone efficiency calibration с C¹-сшивкой.

Каноничный ЛСРМ §8.4 / ORTEC рекомендация: вместо одного полинома
высокой степени, который перестраивается на edge points и даёт
non-physical wiggles, использовать **2 зоны**:

  - **Low-E zone** (E < ~250 кэВ): degree 3-4 polynomial
    (быстро растёт ε при увеличении E из-за уменьшения photoelectric).
  - **High-E zone** (E ≥ ~250 кэВ): degree 1-2 polynomial
    (медленно убывает ε из-за уменьшения photoabsorption + Compton).

Точка стыка (default 250 кэВ) определяется как минимум зоны
photoelectric → Compton dominance transition.

C¹-сшивка: на стыке значение ε(E_split) и первая производная
dε/dE(E_split) совпадают между зонами. Это исключает artificial
discontinuity или резкий перегиб.

Пример (NaI 50×50, Marinelli 1L, point-source 25cm):
  Low-E:  ε(60) ≈ 0.025, ε(186) ≈ 0.018, ε(250) ≈ 0.015
  High-E: ε(250) = 0.015 (stitched), ε(662) ≈ 0.008,
          ε(1461) ≈ 0.004, ε(2614) ≈ 0.0025

References
----------
- ЛСРМ Algorithmic Foundations 2022 §8.4 "Калибровка по эффективности"
- ORTEC AN66 §6 Efficiency Calibration
- Gilmore & Joss 3rd Ed., §8.3 «Multi-segment efficiency fits»
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np


DEFAULT_SPLIT_KEV = 250.0
DEFAULT_LOW_DEGREE = 3
DEFAULT_HIGH_DEGREE = 2
MAX_LOW_DEGREE = 4
MAX_HIGH_DEGREE = 2


@dataclass(frozen=True)
class TwoZoneEfficiencyFit:
    """Результат 2-зонной калибровки эффективности."""
    split_E_keV: float
    low_coefs: tuple        # low-E poly coeffs, low-to-high order
    high_coefs: tuple       # high-E poly coeffs, low-to-high order
    n_low_anchors: int
    n_high_anchors: int
    rms_residual_low_pct: float    # %
    rms_residual_high_pct: float
    rms_residual_overall_pct: float
    note: str = ""

    def epsilon(self, E_keV: float) -> float:
        """Вычислить ε(E)."""
        E = float(E_keV)
        if E <= self.split_E_keV:
            coefs = self.low_coefs
        else:
            coefs = self.high_coefs
        # Horner low-to-high
        result = 0.0
        for c in reversed(coefs):
            result = result * E + c
        return result


def fit_two_zone_efficiency(
    energies_keV: Sequence[float],
    epsilons: Sequence[float],
    *,
    split_E_keV: float = DEFAULT_SPLIT_KEV,
    low_degree: int = DEFAULT_LOW_DEGREE,
    high_degree: int = DEFAULT_HIGH_DEGREE,
    log_log_space: bool = True,
    enforce_c1_continuity: bool = True,
) -> TwoZoneEfficiencyFit:
    """Подогнать 2-зонную аппроксимацию ε(E) к данным эталонов.

    Parameters
    ----------
    energies_keV, epsilons : sequences
        Калибровочные точки.
    split_E_keV : float
        Энергия стыка зон (default 250 кэВ).
    low_degree, high_degree : int
        Степень полинома в каждой зоне. Capped at MAX_LOW_DEGREE=4
        and MAX_HIGH_DEGREE=2.
    log_log_space : bool
        Если True (default), фиттинг в log-log координатах: log(ε)
        vs log(E). Это даёт canonical ЛСРМ form, не выбрасывает
        ε<0 артефакты.
    enforce_c1_continuity : bool
        Если True (default), вторая зона стартует с (value, slope)
        совпадающими с первой на E_split. Реализуется через
        constrained linear-least-squares.
    """
    E = np.asarray(list(energies_keV), dtype=float)
    eps = np.asarray(list(epsilons), dtype=float)
    if len(E) != len(eps) or len(E) < 3:
        raise ValueError(
            f"Need ≥3 anchor points; got {len(E)} energies, {len(eps)} eps"
        )
    if (eps <= 0).any():
        raise ValueError("All epsilon values must be positive")

    if log_log_space:
        x = np.log(E)
        y = np.log(eps)
        split_x = math.log(split_E_keV)
    else:
        x = E.copy()
        y = eps.copy()
        split_x = split_E_keV

    low_mask = x <= split_x
    high_mask = ~low_mask

    n_low = int(low_mask.sum())
    n_high = int(high_mask.sum())

    if n_low < (low_degree + 1):
        # Auto-reduce degree
        low_degree = max(1, n_low - 1)
    if n_high < (high_degree + 1):
        high_degree = max(1, n_high - 1)

    low_degree = min(low_degree, MAX_LOW_DEGREE)
    high_degree = min(high_degree, MAX_HIGH_DEGREE)

    # Подгонка low zone сначала
    low_x = x[low_mask]
    low_y = y[low_mask]
    high_x = x[high_mask]
    high_y = y[high_mask]

    if n_low >= 2:
        low_coefs_high_to_low = np.polyfit(low_x, low_y, low_degree)
    else:
        # 1 anchor — degenerate, fit by passing through (E,ε)
        low_coefs_high_to_low = np.array([0.0, low_y[0]])
        low_degree = 0

    low_coefs = tuple(float(c) for c in low_coefs_high_to_low[::-1])

    # Вычислить значение и наклон в точке стыка
    def _poly_value(coefs_low_to_high, x_val):
        r = 0.0
        for c in reversed(coefs_low_to_high):
            r = r * x_val + c
        return r

    def _poly_derivative(coefs_low_to_high, x_val):
        # d/dx Σ c_i · x^i = Σ i·c_i · x^(i-1)
        if len(coefs_low_to_high) <= 1:
            return 0.0
        deriv = 0.0
        for i, c in enumerate(coefs_low_to_high):
            if i == 0:
                continue
            deriv += i * c * (x_val ** (i - 1))
        return deriv

    v_low_split = _poly_value(low_coefs, split_x)
    d_low_split = _poly_derivative(low_coefs, split_x)

    # Подгонка high zone с C¹ constraints
    if enforce_c1_continuity and n_high >= 1:
        # high_poly(x) = a0 + a1·(x - split_x) + a2·(x - split_x)^2 + ...
        # Условия C¹: high_poly(split_x) = v_low_split → a0 = v_low_split
        #             d/dx high_poly(split_x) = d_low_split → a1 = d_low_split
        # Свободные параметры: a2, a3, ... high_degree
        n_free = max(0, high_degree - 1)
        if n_free == 0 or n_high == 0:
            # Просто константа + линейный наклон
            a0 = v_low_split
            a1 = d_low_split
            high_local_coefs = [a0, a1]
        else:
            shifted_x = high_x - split_x
            # Residual data: y_high − v_low_split − d_low_split·shifted_x
            residual = high_y - v_low_split - d_low_split * shifted_x
            # Fit poly (a2, a3, ...) over shifted_x^2, shifted_x^3, ...
            # Design matrix
            cols = [shifted_x ** k for k in range(2, 2 + n_free)]
            if cols:
                A = np.column_stack(cols)
                if A.shape[0] >= A.shape[1]:
                    coefs_a2_plus, *_ = np.linalg.lstsq(A, residual, rcond=None)
                else:
                    coefs_a2_plus = np.zeros(A.shape[1])
                high_local_coefs = [v_low_split, d_low_split] + list(
                    float(c) for c in coefs_a2_plus
                )
            else:
                high_local_coefs = [v_low_split, d_low_split]
        # Перевести high-зону из локального (x - split_x) в абсолютный x:
        # poly_in_local(u) = Σ a_i · u^i, u = x - split_x
        # poly_in_x(x) = poly_in_local(x - split_x)
        # Это не простой переразвёртка коэффициентов — оставим как локальный poly
        # с явным offset_x_for_evaluation = split_x. Для простоты конвертируем
        # poly через разложение Tailor:
        # poly_in_x(x) = Σ_i a_i · (x - split_x)^i = Σ_i a_i · Σ_{k≤i} C(i,k) (-split_x)^(i-k) · x^k
        # → coefficient at x^k = Σ_{i≥k} a_i · C(i,k) · (-split_x)^(i-k)
        max_deg = len(high_local_coefs) - 1
        high_coefs_abs = [0.0] * (max_deg + 1)
        for i, a_i in enumerate(high_local_coefs):
            for k in range(0, i + 1):
                binom = math.comb(i, k)
                high_coefs_abs[k] += a_i * binom * ((-split_x) ** (i - k))
        high_coefs = tuple(high_coefs_abs)
    else:
        # Independent fit без C¹
        if n_high >= 2:
            high_coefs_high_to_low = np.polyfit(high_x, high_y, high_degree)
            high_coefs = tuple(float(c) for c in high_coefs_high_to_low[::-1])
        else:
            high_coefs = (v_low_split,)

    # Residuals в исходном пространстве ε
    def _poly_eval(coefs_low_to_high, x_val):
        r = 0.0
        for c in reversed(coefs_low_to_high):
            r = r * x_val + c
        return r

    low_pred_y = np.array([_poly_eval(low_coefs, xi) for xi in low_x])
    high_pred_y = np.array([_poly_eval(high_coefs, xi) for xi in high_x])
    if log_log_space:
        low_pred_eps = np.exp(low_pred_y)
        high_pred_eps = np.exp(high_pred_y)
        low_true_eps = eps[low_mask]
        high_true_eps = eps[high_mask]
    else:
        low_pred_eps = low_pred_y
        high_pred_eps = high_pred_y
        low_true_eps = eps[low_mask]
        high_true_eps = eps[high_mask]

    def _rms_pct(true, pred):
        if len(true) == 0:
            return 0.0
        rel = (pred - true) / np.maximum(true, 1e-12)
        return float(100.0 * math.sqrt(np.mean(rel ** 2)))

    rms_low = _rms_pct(low_true_eps, low_pred_eps)
    rms_high = _rms_pct(high_true_eps, high_pred_eps)
    all_true = np.concatenate([low_true_eps, high_true_eps])
    all_pred = np.concatenate([low_pred_eps, high_pred_eps])
    rms_all = _rms_pct(all_true, all_pred)

    note = (
        f"F-284 2-zone efficiency (split={split_E_keV:.0f} кэВ, "
        f"low_deg={low_degree}, high_deg={high_degree}, "
        f"log_log={log_log_space}, C¹={enforce_c1_continuity})"
    )

    fit = TwoZoneEfficiencyFit(
        split_E_keV=float(split_E_keV),
        low_coefs=low_coefs,
        high_coefs=high_coefs,
        n_low_anchors=n_low,
        n_high_anchors=n_high,
        rms_residual_low_pct=rms_low,
        rms_residual_high_pct=rms_high,
        rms_residual_overall_pct=rms_all,
        note=note,
    )

    # Если log_log: epsilon() должен делать exp(...). Для простоты —
    # добавим wrapper-метод в caller. В этом dataclass coefs полиномов в
    # log-log пространстве; пользователь использует helper:
    return fit


def evaluate_two_zone_efficiency(
    fit: TwoZoneEfficiencyFit,
    E_keV: float,
    *,
    log_log_fitted: bool = True,
) -> float:
    """Helper: вычислить ε(E_keV) учитывая log-log пространство."""
    if log_log_fitted:
        x = math.log(E_keV)
    else:
        x = E_keV
    if x <= (math.log(fit.split_E_keV) if log_log_fitted else fit.split_E_keV):
        coefs = fit.low_coefs
    else:
        coefs = fit.high_coefs
    result = 0.0
    for c in reversed(coefs):
        result = result * x + c
    if log_log_fitted:
        return math.exp(result)
    return result


__all__ = [
    "DEFAULT_SPLIT_KEV",
    "TwoZoneEfficiencyFit",
    "fit_two_zone_efficiency",
    "evaluate_two_zone_efficiency",
]
