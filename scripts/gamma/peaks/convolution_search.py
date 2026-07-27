"""F-124 (v1.17.6) — Convolution-based peak search (D-20).

Альтернативный метод поиска пиков на основе свёртки спектра с
матчинг-фильтром «forward-Gaussian»:

    R(ch) = Σ_k counts[ch+k] · g(k)

где g(k) — нормированный единично-площадный гауссиан с σ = FWHM/2.355
в окне ±N·σ. Локальные максимумы R(ch) дают позиции пиков, локальный
уровень континуума — оценку фона.

Это близко к классическому matched-filter detection (Gilmore & Joss
§9.3) и работает несколько иначе, чем Mariscotti second-derivative:

  • Mariscotti чувствителен к «горбу» — реагирует на ВТОРУЮ
    производную, обнаруживая локальный максимум кривизны.
  • Convolution-search чувствителен к ФОРМЕ — реагирует на
    локальное «накопление» сигнала вокруг гауссиана.

Эти два подхода дают слегка разные ответы на близких дублетах
и широких пиках с наклоном континуума. Для контроля методологии
(D-20: «не будет ли лучше искать пики методом свёртки? Или прогнать
оба метода и сравнить результаты?») мы прогоняем оба и сравниваем.

API
---
``convolution_peak_search(counts, fwhm_channels, ...)`` → list[FoundPeak]
    Тот же ``FoundPeak`` dataclass, что у Mariscotti — для совместимости
    с downstream-кодом.

``compare_peak_methods(peaks_a, peaks_b, tolerance_channels=1.5)`` → dict
    Сравнение списков по позициям. Возвращает agreed / a_only / b_only.

Reference: Gilmore & Joss "Practical Gamma-ray Spectrometry" 3rd Ed.,
§9.3 "Peak Search Methods" — matched-filter detection.
"""
from __future__ import annotations

import math
from typing import Callable, List, Optional, Union

import numpy as np

from gamma.peaks.search import FoundPeak, FwhmSpec, SigmaSpec


SQRT_2PI = math.sqrt(2.0 * math.pi)


def _fwhm_at(fwhm_channels: FwhmSpec, ch: int) -> float:
    """Resolve scalar-or-callable FWHM spec to a float at given channel."""
    if callable(fwhm_channels):
        return float(fwhm_channels(ch))
    return float(fwhm_channels)


def _sigma_at(sigma_threshold: SigmaSpec, ch: int) -> float:
    if callable(sigma_threshold):
        return float(sigma_threshold(ch))
    return float(sigma_threshold)


def _gaussian_kernel(half_width_channels: int, sigma_ch: float) -> np.ndarray:
    """Unit-area discrete Gaussian kernel в диапазоне [-N, +N] каналов."""
    k = np.arange(-half_width_channels, half_width_channels + 1)
    g = np.exp(-(k ** 2) * 0.5 / (sigma_ch ** 2))
    s = float(g.sum())
    if s <= 0:
        return g
    return g / s


def _band_segments(
    n_channels: int, fwhm_channels: FwhmSpec, max_ratio: float = 1.3,
) -> List[tuple]:
    """Разбиение спектра на полосы, в пределах которых FWHM меняется
    не более чем в ``max_ratio`` раз. Возвращает [(lo, hi), ...].
    """
    if not callable(fwhm_channels):
        return [(0, n_channels)]
    segments: List[tuple] = []
    lo = 0
    fwhm_lo = _fwhm_at(fwhm_channels, lo)
    for ch in range(1, n_channels):
        fwhm_here = _fwhm_at(fwhm_channels, ch)
        if fwhm_here <= 0 or fwhm_lo <= 0:
            continue
        ratio = fwhm_here / fwhm_lo
        if ratio > max_ratio or ratio < (1.0 / max_ratio):
            segments.append((lo, ch))
            lo = ch
            fwhm_lo = fwhm_here
    segments.append((lo, n_channels))
    return segments


def convolution_peak_search(
    counts: np.ndarray,
    *,
    fwhm_channels: FwhmSpec,
    sigma_threshold: SigmaSpec = 2.5,
    min_separation_factor: float = 0.6,
    edge_margin: int = 10,
    kernel_half_widths_sigma: float = 3.0,
) -> List[FoundPeak]:
    """Найти пики методом свёртки с матчинг-фильтром-гауссианом.

    Алгоритм:
      1. Для каждой банды (где FWHM ≈ const) собрать ядро-гауссиан
         с σ = FWHM/2.355.
      2. Свернуть с counts через ``np.convolve(mode='same')`` —
         результат R(ch) — «отклик» фильтра.
      3. Оценить локальный континуум B(ch) скользящим минимумом
         в окне ±3·σ → ROI net = R(ch) − B(ch)·kernel_size.
      4. Найти локальные максимумы R(ch).
      5. Currie significance: significance = (R − B) / √(B·Σg²),
         где Σg² — норма ядра (F-269: убран лишний K-фактор).
      6. Отсеять пики с significance < sigma_threshold и
         расстоянием < min_separation_factor·FWHM от уже принятого.

    Параметры идентичны ``mariscotti_search`` для drop-in
    совместимости.
    """
    counts_arr = np.asarray(counts, dtype=np.float64)
    n = len(counts_arr)
    if n < 10:
        return []

    segments = _band_segments(n, fwhm_channels)
    significance = np.zeros(n, dtype=np.float64)
    background_per_ch = np.zeros(n, dtype=np.float64)
    height_per_ch = np.zeros(n, dtype=np.float64)
    fwhm_per_ch = np.zeros(n, dtype=np.float64)

    for lo, hi in segments:
        mid = (lo + hi) // 2
        fwhm_loc = max(1.0, _fwhm_at(fwhm_channels, mid))
        sigma_ch = fwhm_loc / 2.355
        half_w = max(2, int(round(kernel_half_widths_sigma * sigma_ch)))
        kernel = _gaussian_kernel(half_w, sigma_ch)
        K = len(kernel)
        # Свёртка R(ch) = Σ counts[ch+k]·g(k)
        # np.convolve(mode='same') вписывает результат в исходный
        # длину; на краях возможны искажения, поэтому учитываем
        # edge_margin при поиске максимумов.
        # Берём срез counts с буфером для свёртки.
        buf_lo = max(0, lo - half_w)
        buf_hi = min(n, hi + half_w)
        local_counts = counts_arr[buf_lo:buf_hi]
        if len(local_counts) < K:
            continue
        R = np.convolve(local_counts, kernel, mode="same")
        # R[i] соответствует local_counts[i] = counts[buf_lo + i]
        # Континуум: скользящее окно равной длины K — берём среднее
        # «крыльев» (без центральных ±sigma_ch каналов).
        window = K
        wing = max(1, int(round(0.5 * sigma_ch)))
        # Скользящее среднее
        if len(local_counts) >= window:
            kernel_uniform = np.ones(window) / window
            B = np.convolve(local_counts, kernel_uniform, mode="same")
        else:
            B = np.full_like(local_counts, float(np.mean(local_counts)))
        # Оценка сигнала: net = R − B·Σg
        # (Σg = 1 по нормировке, но для R = свёртка с unit-area
        # gauss-kernel net ≈ R − B.)
        net = R - B
        # F-269 (v1.17.11, T-018) — matched-filter variance.
        # Для R = Σ_k g(k)·X_k с независимыми Poisson X_k (Var=B per
        # channel) дисперсия отклика:
        #     Var(R) = Σ_k g(k)² · Var(X_k) = B · Σg²
        # Лишний множитель K (длина ядра), стоявший здесь до v1.17.10,
        # завышал σ_R в √K раз и недооценивал significance в √K раз —
        # из-за чего слабые пики не проходили sigma_threshold.
        # См. Gilmore & Joss §9.3 «Matched-filter detection».
        sum_g2 = float(np.sum(kernel * kernel))
        # Защита от деления на ноль
        sigma_R = np.sqrt(np.maximum(B * sum_g2, 1.0))
        sig = net / np.maximum(sigma_R, 1e-9)
        # Пишем в глобальный массив
        for i_local in range(len(local_counts)):
            i_global = buf_lo + i_local
            if i_global < lo or i_global >= hi:
                continue
            if sig[i_local] > significance[i_global]:
                significance[i_global] = float(sig[i_local])
                background_per_ch[i_global] = float(B[i_local])
                height_per_ch[i_global] = float(max(0.0, R[i_local] - B[i_local]))
                fwhm_per_ch[i_global] = float(fwhm_loc)

    # Локальные максимумы significance, отсев по edge / min_separation
    peaks: List[FoundPeak] = []
    candidates: List[tuple] = []  # (ch, significance)
    for ch in range(edge_margin, n - edge_margin):
        thr = _sigma_at(sigma_threshold, ch)
        if significance[ch] < thr:
            continue
        # Локальный максимум: > соседей в окне ±1
        if (significance[ch] >= significance[ch - 1]
                and significance[ch] >= significance[ch + 1]):
            candidates.append((ch, significance[ch]))
    # Отсев по min_separation_factor·FWHM
    candidates.sort(key=lambda x: -x[1])  # сильнейшие первые
    accepted_channels: List[int] = []
    for ch, sig in candidates:
        fwhm_here = fwhm_per_ch[ch] if fwhm_per_ch[ch] > 0 \
            else _fwhm_at(fwhm_channels, ch)
        min_sep = max(1.0, min_separation_factor * fwhm_here)
        if any(abs(ch - acc) < min_sep for acc in accepted_channels):
            continue
        accepted_channels.append(ch)
        # area_estimate: net·K ≈ Σ_k g(k) ≈ 1 → R − B
        height = height_per_ch[ch]
        # Площадь = height · σ · √(2π) для гауссиана
        sigma_ch = max(1.0, fwhm_here / 2.355)
        area = float(height * sigma_ch * SQRT_2PI)
        peaks.append(FoundPeak(
            channel=int(ch),
            height=float(height),
            fwhm_channels=float(fwhm_here),
            significance=float(sig),
            area_estimate=float(area),
            sigma_area_estimate=float(math.sqrt(max(area, 1.0))),
        ))
    # Сортировка по каналу для удобства downstream
    peaks.sort(key=lambda p: p.channel)
    return peaks


def compare_peak_methods(
    peaks_a: List[FoundPeak],
    peaks_b: List[FoundPeak],
    *,
    tolerance_channels: float = 1.5,
) -> dict:
    """Сравнить два списка пиков по позициям (channel).

    Возвращает словарь:
      • ``agreed``    — список пар (peak_a, peak_b), совпавших в пределах
                        ``tolerance_channels``.
      • ``a_only``    — пики только из метода A.
      • ``b_only``    — пики только из метода B.
      • ``agreement_fraction`` — |agreed| / max(|A|, |B|).
      • ``mean_residual_channels`` — средняя |ch_a - ch_b| для agreed.

    Используется как методологический контроль (D-20): если методы
    расходятся > 5% по числу пиков или > 1·FWHM по позициям, нужен
    разбор.
    """
    used_b: set = set()
    agreed: List[tuple] = []
    for pa in peaks_a:
        # Найти ближайший pb в пределах tolerance, ещё не привязанный
        best_i = -1
        best_d = float("inf")
        for i, pb in enumerate(peaks_b):
            if i in used_b:
                continue
            d = abs(pa.channel - pb.channel)
            if d < tolerance_channels and d < best_d:
                best_d = d
                best_i = i
        if best_i >= 0:
            used_b.add(best_i)
            agreed.append((pa, peaks_b[best_i]))
    a_only = [pa for pa in peaks_a
              if not any(pa is x for x, _ in agreed)]
    b_only = [pb for i, pb in enumerate(peaks_b) if i not in used_b]
    if peaks_a or peaks_b:
        agreement_fraction = (
            len(agreed) / max(len(peaks_a), len(peaks_b), 1)
        )
    else:
        agreement_fraction = 1.0
    if agreed:
        mean_res = sum(abs(pa.channel - pb.channel) for pa, pb in agreed)
        mean_res /= len(agreed)
    else:
        mean_res = 0.0
    return {
        "agreed": agreed,
        "a_only": a_only,
        "b_only": b_only,
        "agreement_fraction": float(agreement_fraction),
        "mean_residual_channels": float(mean_res),
        "n_a": len(peaks_a),
        "n_b": len(peaks_b),
        "tolerance_channels": float(tolerance_channels),
    }


__all__ = [
    "convolution_peak_search",
    "compare_peak_methods",
]
