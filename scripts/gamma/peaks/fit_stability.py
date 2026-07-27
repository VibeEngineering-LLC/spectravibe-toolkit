"""
F-281 (v1.17.13, T-052) — Sequential nonlinear-parameter exclusion.

ЛСРМ §15.2 / Будыка §7.7 — при подгонке мультиплета на NaI слабые
компоненты могут давать нестабильный нелинейный fit, когда число
свободных нелинейных параметров (FWHM, position, step) ≫ информация
в данных. Каноничный рецепт: **последовательно** исключать наименее
определённые параметры:

  1. dS/S > 0.05 для step → исключить step (зафиксировать h_step = 0)
  2. dS/S > 0.10 для FWHM → исключить FWHM (зафиксировать sigma на
     калиброванном значении)
  3. dS/S > 1.00 для всех нелинейных → fallback на линейный fit с
     зафиксированными centroids / FWHM / step
  4. Если даже линейный fit не сходится — вернуть NaN с notes.

Эта стратегия даёт robust convergence на слабых пиках вблизи MDA
(где статистика < 3-5σ над фоном) ценой потери одной-двух степеней
свободы. Для сильных пиков (S/dS < 0.05) изменений нет.

References
----------
- ЛСРМ Алгоритмические основы 2022 §15.2
- Будыка А.К. "Применение полупроводниковых детекторов..." 2021 §7.7
- F-281 contract: применяется как post-process после первого scipy.curve_fit
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


# Каноничные пороги ЛСРМ §15.2.
EXCLUDE_STEP_THRESHOLD = 0.05      # dS/S > 5 % → drop step
EXCLUDE_FWHM_THRESHOLD = 0.10      # dS/S > 10 % → drop FWHM
EXCLUDE_ALL_NL_THRESHOLD = 1.00    # dS/S > 100 % → fallback to linear


@dataclass(frozen=True)
class NonlinearExclusionDecision:
    """Решение пайплайна о том, какие нелинейные параметры исключить."""
    drop_step: bool
    drop_fwhm: bool
    drop_all_nl: bool
    reason: str           # человекочитаемое объяснение
    triggering_dS_over_S: float


def decide_nonlinear_exclusions(
    relative_S_uncertainty: float,
) -> NonlinearExclusionDecision:
    """Применить ЛСРМ §15.2 правила exclusion к фиту мультиплета.

    Parameters
    ----------
    relative_S_uncertainty : float
        dS/S, где S — суммарная площадь мультиплета, dS — её 1σ
        uncertainty из ковариационной матрицы. Безразмерная величина.

    Returns
    -------
    NonlinearExclusionDecision
        Какие группы нелинейных параметров исключить из следующей
        итерации. Caller сам реализует exclusion (фиксирование
        h_step=0, sigma=cal_value, и т.д.).
    """
    dS_S = float(relative_S_uncertainty)
    if dS_S > EXCLUDE_ALL_NL_THRESHOLD:
        return NonlinearExclusionDecision(
            drop_step=True, drop_fwhm=True, drop_all_nl=True,
            reason=(
                f"dS/S = {dS_S:.2f} > {EXCLUDE_ALL_NL_THRESHOLD:.2f} → "
                f"исключить ВСЕ нелинейные параметры (linear fit на "
                f"калиброванных centroid/FWHM/step)"
            ),
            triggering_dS_over_S=dS_S,
        )
    if dS_S > EXCLUDE_FWHM_THRESHOLD:
        return NonlinearExclusionDecision(
            drop_step=True, drop_fwhm=True, drop_all_nl=False,
            reason=(
                f"dS/S = {dS_S:.2%} > {EXCLUDE_FWHM_THRESHOLD:.0%} → "
                f"исключить FWHM (на калиброванном значении) и step"
            ),
            triggering_dS_over_S=dS_S,
        )
    if dS_S > EXCLUDE_STEP_THRESHOLD:
        return NonlinearExclusionDecision(
            drop_step=True, drop_fwhm=False, drop_all_nl=False,
            reason=(
                f"dS/S = {dS_S:.2%} > {EXCLUDE_STEP_THRESHOLD:.0%} → "
                f"исключить только step (h_step → 0)"
            ),
            triggering_dS_over_S=dS_S,
        )
    return NonlinearExclusionDecision(
        drop_step=False, drop_fwhm=False, drop_all_nl=False,
        reason=f"dS/S = {dS_S:.2%} ≤ {EXCLUDE_STEP_THRESHOLD:.0%} — fit стабилен",
        triggering_dS_over_S=dS_S,
    )


__all__ = [
    "EXCLUDE_STEP_THRESHOLD",
    "EXCLUDE_FWHM_THRESHOLD",
    "EXCLUDE_ALL_NL_THRESHOLD",
    "NonlinearExclusionDecision",
    "decide_nonlinear_exclusions",
]
