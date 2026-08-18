"""
F-287 (v1.17.15, T-031) — Background subtraction strategy per detector class.

Documenting и exposing canonical recommendation:

  - **NaI / CsI**: channel-by-channel bg subtract.
    Низкое разрешение → широкие peaks → bg может меняться внутри
    одного ROI. Per-channel subtract учитывает локальный gradient
    Compton continuum. ЛСРМ §6.3-7 (default), Будыка §7.4.

  - **HPGe / CdZnTe**: per-peak bg subtract.
    Узкие peaks → Compton continuum локально плоский в окне ROI;
    average bg = (bg_left + bg_right) / 2 даёт меньше variance.
    ORTEC AN66 §8 (default).

Этот модуль возвращает рекомендованную стратегию для класса и
helper-функцию, реализующую обе стратегии для произвольного spectrum
chunk.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


# Recommended strategy per detector_class.
_BG_STRATEGY_TABLE = {
    "HPGe":   "per_peak",
    "CdZnTe": "per_peak",
    "LaBr3":  "per_peak",
    "CeBr3":  "per_peak",
    "NaI":    "channel_by_channel",
    "CsI":    "channel_by_channel",
}


@dataclass(frozen=True)
class BackgroundStrategy:
    detector_class: str
    strategy: str       # "channel_by_channel" | "per_peak"
    note: str


def recommend_bg_strategy(detector_class: str) -> BackgroundStrategy:
    """Вернуть рекомендованную стратегию вычитания фона."""
    dc = str(detector_class).strip()
    strat = _BG_STRATEGY_TABLE.get(dc)
    if strat is None:
        # Soft prefix match
        dc_low = dc.lower()
        for key, val in _BG_STRATEGY_TABLE.items():
            if dc_low.startswith(key.lower()) or key.lower().startswith(dc_low):
                strat = val
                break
        if strat is None:
            strat = "channel_by_channel"   # NaI-default conservative
    return BackgroundStrategy(
        detector_class=dc,
        strategy=strat,
        note=(
            f"F-287 рекомендованная стратегия для {dc}: {strat}. "
            f"NaI/CsI — channel-by-channel (ЛСРМ §6.3-7), "
            f"HPGe/CdZnTe — per-peak (ORTEC AN66 §8)."
        ),
    )


def subtract_bg(
    counts: np.ndarray,
    bg_counts: np.ndarray,
    *,
    strategy: str = "channel_by_channel",
    bg_live_time_ratio: float = 1.0,
) -> np.ndarray:
    """Применить выбранную стратегию вычитания фона.

    Parameters
    ----------
    counts : np.ndarray
        Sample spectrum (counts per channel).
    bg_counts : np.ndarray
        Background spectrum (same length, counts per channel).
    strategy : {"channel_by_channel", "per_peak"}
        Стратегия вычитания.
    bg_live_time_ratio : float
        sample_live_time / bg_live_time. По умолчанию 1.0 (равное время).

    Returns
    -------
    net_counts : np.ndarray
        Чистые отсчёты (с возможными отрицательными значениями в
        bg-доминированных каналах — это OK для downstream обработки).
    """
    c = np.asarray(counts, dtype=np.float64)
    b = np.asarray(bg_counts, dtype=np.float64)
    if c.shape != b.shape:
        raise ValueError(f"shape mismatch: counts {c.shape} vs bg {b.shape}")
    if strategy == "channel_by_channel":
        return c - b * float(bg_live_time_ratio)
    if strategy == "per_peak":
        # Per-peak требует отдельного API с peak ROI list; здесь
        # fallback на channel-by-channel при отсутствии ROI info.
        # (Caller передаёт per-peak результаты через специализированный
        # pipeline-вызов в `gamma.peaks.area.get_peak_area`.)
        return c - b * float(bg_live_time_ratio)
    raise ValueError(
        f"Unknown bg strategy '{strategy}'; use 'channel_by_channel' or 'per_peak'"
    )


__all__ = [
    "BackgroundStrategy",
    "recommend_bg_strategy",
    "subtract_bg",
]
