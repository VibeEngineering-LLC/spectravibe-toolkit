"""
F-300 (v1.17.21, T-021b) — Log-spline interpolation для tabulated peak shape
                            (ЛСРМ §8.4.3).

Между anchor-точками (см. F-299 `peak_image_tabulated.py`) параметры
формы пика интерполируются в **log-log space**, что соответствует
физическому поведению:

  • FWHM(E) ~ √E (statistics) + const (electronics): power-law хорошо
    приближается прямой в log-log.
  • tail_fraction(E) обычно убывает с E как E^(-0.5)..E^(-1).
  • step_height_frac(E) монотонно растёт с E (Compton scattering).

Алгоритм:
  1. Если anchor нашёлся точно (в пределах tolerance) — вернуть anchor.
  2. Иначе найти 2 ближайших anchor по обе стороны.
  3. Log-log линейная интерполяция (для FWHM / tail / step параметров).
  4. Linear interp для phase-параметров (asymmetry — не log).

Поведение на границах:
  • E < E_min(anchors) → extrapolate log-log (с warning, если задан).
  • E > E_max(anchors) → extrapolate log-log (с warning).

Альтернативы (deferred):
  • Cubic spline через scipy (доступно, но добавляет hard-dep).
  • Catmull-Rom smoothing для visually-smooth shapes.
  • Bezier для precision peak templates.

References
----------
- ЛСРМ Algorithmic Foundations 2022 § 8.4.3 «Интерполяция между
  anchor-точками»
- Numerical Recipes 3rd Ed. § 3.4 (linear interpolation)
- Press WH «Cubic splines» NR § 3.3 (если перейдём на cubic)
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

# Импорт TabulatedPeakImage делается lazy (избегаем circular deps).


@dataclass(frozen=True)
class InterpolatedPeakShape:
    """Результат интерполяции на конкретной энергии E."""

    E_keV: float
    fwhm_keV: float
    tail_fraction: float
    tail_slope_inv_keV: float
    step_height_frac: float
    asymmetry: float
    was_extrapolated: bool      # True если E вне диапазона anchors


def _log_log_interp(x: float, x0: float, y0: float, x1: float, y1: float,
                    ) -> float:
    """Линейная log-log интерполяция между (x0,y0) и (x1,y1)."""
    if y0 <= 0 or y1 <= 0:
        # Fallback: линейная если хоть один y неположителен
        t = (x - x0) / (x1 - x0) if x1 != x0 else 0.0
        return y0 + t * (y1 - y0)
    if x <= 0 or x0 <= 0 or x1 <= 0:
        t = (x - x0) / (x1 - x0) if x1 != x0 else 0.0
        return y0 + t * (y1 - y0)
    ln_x, ln_x0, ln_x1 = math.log(x), math.log(x0), math.log(x1)
    ln_y0, ln_y1 = math.log(y0), math.log(y1)
    if ln_x1 == ln_x0:
        return y0
    frac = (ln_x - ln_x0) / (ln_x1 - ln_x0)
    return math.exp(ln_y0 + frac * (ln_y1 - ln_y0))


def _linear_interp(x: float, x0: float, y0: float, x1: float, y1: float,
                   ) -> float:
    if x1 == x0:
        return y0
    t = (x - x0) / (x1 - x0)
    return y0 + t * (y1 - y0)


def interpolate_peak_shape(
    tabulated_peak_image,        # TabulatedPeakImage (lazy-typed)
    E_keV: float,
    exact_match_tolerance_keV: float = 0.5,
) -> InterpolatedPeakShape:
    """Интерполировать peak shape параметры на произвольную энергию.

    Parameters
    ----------
    tabulated_peak_image : TabulatedPeakImage
        Из `peak_image_tabulated`.
    E_keV : float
        Энергия, для которой нужна форма.
    exact_match_tolerance_keV : float
        Если anchor найден в пределах этого допуска — вернуть его без interp.

    Returns
    -------
    InterpolatedPeakShape.
    """
    if not tabulated_peak_image.anchors:
        raise ValueError("Empty anchors list — невозможно интерполировать")
    if E_keV <= 0:
        raise ValueError(f"E_keV={E_keV} must be > 0")

    sorted_a = sorted(
        tabulated_peak_image.anchors, key=lambda a: a.E_keV,
    )

    # Exact match
    exact = tabulated_peak_image.anchor_at_E(E_keV, exact_match_tolerance_keV)
    if exact is not None:
        return InterpolatedPeakShape(
            E_keV=E_keV,
            fwhm_keV=exact.fwhm_keV,
            tail_fraction=exact.tail_fraction,
            tail_slope_inv_keV=exact.tail_slope_inv_keV,
            step_height_frac=exact.step_height_frac,
            asymmetry=exact.asymmetry,
            was_extrapolated=False,
        )

    # Extrapolation cases
    if E_keV < sorted_a[0].E_keV:
        a0, a1 = sorted_a[0], sorted_a[1] if len(sorted_a) > 1 else sorted_a[0]
        extrap = True
    elif E_keV > sorted_a[-1].E_keV:
        a0 = sorted_a[-2] if len(sorted_a) > 1 else sorted_a[-1]
        a1 = sorted_a[-1]
        extrap = True
    else:
        # Найти bracket [a0, a1] вокруг E
        a0 = a1 = sorted_a[0]
        for i in range(len(sorted_a) - 1):
            if sorted_a[i].E_keV <= E_keV <= sorted_a[i + 1].E_keV:
                a0, a1 = sorted_a[i], sorted_a[i + 1]
                break
        extrap = False

    return InterpolatedPeakShape(
        E_keV=E_keV,
        fwhm_keV=_log_log_interp(
            E_keV, a0.E_keV, a0.fwhm_keV, a1.E_keV, a1.fwhm_keV,
        ),
        tail_fraction=_log_log_interp(
            E_keV, a0.E_keV, max(a0.tail_fraction, 1e-9),
            a1.E_keV, max(a1.tail_fraction, 1e-9),
        ) if a0.tail_fraction > 0 and a1.tail_fraction > 0
            else _linear_interp(
                E_keV, a0.E_keV, a0.tail_fraction,
                a1.E_keV, a1.tail_fraction,
            ),
        tail_slope_inv_keV=_linear_interp(
            E_keV, a0.E_keV, a0.tail_slope_inv_keV,
            a1.E_keV, a1.tail_slope_inv_keV,
        ),
        step_height_frac=_log_log_interp(
            E_keV, a0.E_keV, max(a0.step_height_frac, 1e-9),
            a1.E_keV, max(a1.step_height_frac, 1e-9),
        ) if a0.step_height_frac > 0 and a1.step_height_frac > 0
            else _linear_interp(
                E_keV, a0.E_keV, a0.step_height_frac,
                a1.E_keV, a1.step_height_frac,
            ),
        asymmetry=_linear_interp(
            E_keV, a0.E_keV, a0.asymmetry, a1.E_keV, a1.asymmetry,
        ),
        was_extrapolated=extrap,
    )


def fwhm_at_E(tabulated_peak_image, E_keV: float) -> float:
    """Convenience: только FWHM на E."""
    return interpolate_peak_shape(tabulated_peak_image, E_keV).fwhm_keV


def batch_interpolate(
    tabulated_peak_image,
    energies_keV: Sequence[float],
) -> List[InterpolatedPeakShape]:
    """Бакетная интерполяция для серии E."""
    return [
        interpolate_peak_shape(tabulated_peak_image, E)
        for E in energies_keV
    ]


__all__ = [
    "InterpolatedPeakShape",
    "interpolate_peak_shape",
    "fwhm_at_E",
    "batch_interpolate",
]
