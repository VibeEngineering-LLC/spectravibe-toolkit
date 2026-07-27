"""
F-280 (v1.17.13, T-040 + T-079) — Background-fit options under the FEP.

Two helpers for advanced background modelling in single-peak ROIs:

  1. **Smoothed-step background under FEP** (T-040 / Будыка eq.7.29)
     Pure-linear baseline under the FEP underestimates net area by ~1-3 %
     on broad NaI peaks sitting on Compton-step continuum (e.g. Cs-137
     661 on K-40 1461 Compton edge). Adding a sigmoid-shape "step"
     anchored at peak position recovers those counts.

  2. **Asymmetric background regions** (T-079)
     When an interfering peak sits just ABOVE the target peak (within
     ~3-4 FWHM), the standard symmetric bg-region [μ-Δ, μ-2Δ] ∪
     [μ+Δ, μ+2Δ] biases the upper estimate. Replace by asymmetric
     regions: e.g. 10 channels below, only 3 channels above.

Both helpers are **additive options** — callers explicitly opt in via
flags. Default coupled_multiplet behaviour is unchanged.

References
----------
- Будыка А.К. "Применение полупроводниковых детекторов..." 2021, eq.7.29
- Knoll "Radiation Detection and Measurement" 4th Ed., §10.5 Compton step
- Gilmore & Joss 3rd Ed., §6.4 Background subtraction strategies
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


# ──────────────────────────────────────────────────────────────────
# T-040 — Smoothed-step background under FEP
# ──────────────────────────────────────────────────────────────────

def smoothed_step_bg(
    E: np.ndarray,
    *,
    E_peak: float,
    sigma_peak: float,
    bg_left: float,
    bg_right: float,
    slope: float = 0.0,
) -> np.ndarray:
    """Будыка eq.7.29 — гладкая ступенька + линейный наклон.

    Модель фона под пиком:

        B(E) = bg_right + (bg_left − bg_right)·S(E)  +  slope·(E − E_peak)

    где S(E) = 0.5·erfc((E − E_peak)/(σ·√2)) — нормированная ступенька,
    переходящая от 1 (слева) к 0 (справа от пика).

    Parameters
    ----------
    E : np.ndarray
        Сетка энергий ROI (keV).
    E_peak : float
        Положение центра пика (keV).
    sigma_peak : float
        Сигма пика (keV) — определяет «ширину» перехода ступеньки.
    bg_left : float
        Уровень фона ЛЕВО от пика (counts).
    bg_right : float
        Уровень фона СПРАВА от пика (counts).
    slope : float
        Линейный наклон (counts / keV); default 0.

    Returns
    -------
    B(E) : np.ndarray
        Фон под пиком.
    """
    if sigma_peak <= 0:
        # Pathological → constant fallback
        return np.full_like(E, 0.5 * (bg_left + bg_right))
    sqrt_2 = math.sqrt(2.0)
    try:
        from scipy.special import erfc
        S = 0.5 * erfc((E - E_peak) / (sigma_peak * sqrt_2))
    except ImportError:
        # numpy-vectorize fallback
        S = np.array([0.5 * math.erfc((e - E_peak) / (sigma_peak * sqrt_2))
                      for e in E])
    return bg_right + (bg_left - bg_right) * S + slope * (E - E_peak)


# ──────────────────────────────────────────────────────────────────
# T-079 — Asymmetric background regions
# ──────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class AsymmetricBgRegions:
    """Параметры для asymmetric bg-окна вокруг пика."""
    below_lo_ch: int     # начало нижнего окна (channel index)
    below_hi_ch: int     # конец нижнего окна (exclusive)
    above_lo_ch: int     # начало верхнего окна
    above_hi_ch: int     # конец верхнего окна
    interferent_above: bool   # True если выше есть интерферирующий пик
    interferent_below: bool   # True если ниже есть интерферирующий пик


def asymmetric_bg_regions(
    peak_channel: int,
    fwhm_channels: float,
    *,
    interferent_above_channel: Optional[int] = None,
    interferent_below_channel: Optional[int] = None,
    default_width_fwhm: float = 1.0,
    narrow_width_fwhm: float = 0.3,
    gap_fwhm: float = 0.5,
) -> AsymmetricBgRegions:
    """Вычислить asymmetric bg-регионы вокруг пика.

    Если интерферирующий пик находится близко (в пределах 4·FWHM от
    target peak), bg-окно на ТУ сторону сужается до `narrow_width_fwhm`.
    Стандартная сторона использует `default_width_fwhm`.

    Geometry::

      peak  ←gap→  bg-window  ←... interferent | clean side: wider window

    Parameters
    ----------
    peak_channel : int
        Центральный канал target peak.
    fwhm_channels : float
        Полная ширина на полувысоте target peak (каналы).
    interferent_above_channel, interferent_below_channel : Optional[int]
        Каналы интерферирующих пиков (если есть).
    default_width_fwhm, narrow_width_fwhm : float
        Width of bg-windows expressed as fractions of FWHM.
    gap_fwhm : float
        Зазор от центра peak до ближайшего края bg-окна (FWHM).
    """
    fwhm = max(1.0, float(fwhm_channels))
    gap = max(1, int(round(gap_fwhm * fwhm)))
    default_w = max(2, int(round(default_width_fwhm * fwhm)))
    narrow_w = max(1, int(round(narrow_width_fwhm * fwhm)))

    # Below: interferent_below within 4·FWHM → narrow window
    below_window_w = default_w
    int_below = False
    if (interferent_below_channel is not None
            and (peak_channel - interferent_below_channel) <= 4.0 * fwhm
            and (peak_channel - interferent_below_channel) > 0):
        below_window_w = narrow_w
        int_below = True

    above_window_w = default_w
    int_above = False
    if (interferent_above_channel is not None
            and (interferent_above_channel - peak_channel) <= 4.0 * fwhm
            and (interferent_above_channel - peak_channel) > 0):
        above_window_w = narrow_w
        int_above = True

    return AsymmetricBgRegions(
        below_lo_ch=peak_channel - gap - below_window_w,
        below_hi_ch=peak_channel - gap,
        above_lo_ch=peak_channel + gap,
        above_hi_ch=peak_channel + gap + above_window_w,
        interferent_above=int_above,
        interferent_below=int_below,
    )


def asymmetric_bg_estimate(
    counts: np.ndarray,
    regions: AsymmetricBgRegions,
) -> Tuple[float, float]:
    """Оценить уровень фона на основе asymmetric регионов.

    Returns
    -------
    (bg_left, bg_right) : tuple of float
        Средний фон counts/channel в нижнем и верхнем окнах.
        Если регион выходит за границы массива — возвращает counts с
        clip.
    """
    n = len(counts)
    lo = max(0, regions.below_lo_ch)
    hi = min(n, regions.below_hi_ch)
    if hi > lo:
        bg_left = float(np.mean(counts[lo:hi]))
    else:
        bg_left = 0.0
    lo2 = max(0, regions.above_lo_ch)
    hi2 = min(n, regions.above_hi_ch)
    if hi2 > lo2:
        bg_right = float(np.mean(counts[lo2:hi2]))
    else:
        bg_right = 0.0
    return bg_left, bg_right


__all__ = [
    "smoothed_step_bg",
    "AsymmetricBgRegions",
    "asymmetric_bg_regions",
    "asymmetric_bg_estimate",
]
