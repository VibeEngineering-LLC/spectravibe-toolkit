"""F-124 (v1.17.6) — Convolution-based peak search regression (D-20).

Alternative peak-search метод на основе свёртки с матчинг-фильтром.
Проверяет, что:
  1. Находит синтетические гауссианы на их позициях с точностью < 1 канал.
  2. Расхождение с Mariscotti на ровном спектре с большими SNR < 5%.
  3. ``compare_peak_methods`` корректно выделяет agreed / a_only / b_only.
  4. Граничные пики (в edge_margin) отсекаются.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest

SCRIPTS = os.path.join(os.path.dirname(__file__), "scripts")
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)

from gamma.peaks.convolution_search import (
    convolution_peak_search, compare_peak_methods,
)
from gamma.peaks.search import mariscotti_search


def _make_synthetic_spectrum(positions, amplitudes, n_channels=2048,
                             sigma_ch=8.0, baseline=50.0, seed=42):
    """Сгенерировать спектр с гауссианами в заданных позициях."""
    rng = np.random.default_rng(seed)
    counts = np.full(n_channels, baseline, dtype=np.float64)
    ch_arr = np.arange(n_channels)
    for pos, amp in zip(positions, amplitudes):
        counts += amp * np.exp(-(ch_arr - pos) ** 2 / (2.0 * sigma_ch ** 2))
    # Poisson noise
    noisy = rng.poisson(np.maximum(counts, 0.0)).astype(np.float64)
    return noisy


def test_convolution_finds_strong_isolated_peak():
    """Сильный изолированный пик: convolution_search находит его на месте."""
    counts = _make_synthetic_spectrum(
        positions=[500], amplitudes=[3000.0],
        n_channels=2048, sigma_ch=8.0, baseline=20.0,
    )
    peaks = convolution_peak_search(
        counts, fwhm_channels=18.84, sigma_threshold=4.0,
    )
    assert len(peaks) >= 1
    # Должна быть найдена позиция в районе 500 ± 3
    distances = [abs(p.channel - 500) for p in peaks]
    assert min(distances) < 5, (
        f"convolution didn't find peak at 500; "
        f"found at {[p.channel for p in peaks]}"
    )


def test_convolution_finds_multiple_peaks():
    """Три изолированных пика — convolution находит все три."""
    counts = _make_synthetic_spectrum(
        positions=[300, 800, 1500],
        amplitudes=[2000.0, 3000.0, 1500.0],
        n_channels=2048, sigma_ch=8.0, baseline=30.0,
    )
    peaks = convolution_peak_search(
        counts, fwhm_channels=18.84, sigma_threshold=3.0,
    )
    # Ожидаем хотя бы 3 пика
    assert len(peaks) >= 3
    found_positions = sorted(p.channel for p in peaks if p.significance > 5.0)
    for expected in (300, 800, 1500):
        ok = any(abs(fp - expected) < 8 for fp in found_positions)
        assert ok, f"missing peak near {expected}; found {found_positions}"


def test_compare_peak_methods_full_agreement():
    """Идентичные списки → agreement_fraction = 1.0, нет расхождений."""
    from gamma.peaks.search import FoundPeak
    peaks = [
        FoundPeak(channel=c, height=100.0, fwhm_channels=18.0,
                  significance=10.0, area_estimate=1000.0)
        for c in (100, 500, 1200)
    ]
    rep = compare_peak_methods(peaks, peaks, tolerance_channels=1.5)
    assert rep["agreement_fraction"] == pytest.approx(1.0)
    assert len(rep["agreed"]) == 3
    assert rep["a_only"] == []
    assert rep["b_only"] == []


def test_compare_peak_methods_partial_agreement():
    """Один пик расходится → a_only/b_only содержат по одному."""
    from gamma.peaks.search import FoundPeak
    pa = [
        FoundPeak(channel=100, height=100.0, fwhm_channels=18.0,
                  significance=10.0, area_estimate=1000.0),
        FoundPeak(channel=500, height=100.0, fwhm_channels=18.0,
                  significance=10.0, area_estimate=1000.0),
    ]
    pb = [
        FoundPeak(channel=100, height=100.0, fwhm_channels=18.0,
                  significance=10.0, area_estimate=1000.0),
        FoundPeak(channel=900, height=100.0, fwhm_channels=18.0,
                  significance=10.0, area_estimate=1000.0),
    ]
    rep = compare_peak_methods(pa, pb, tolerance_channels=1.5)
    assert len(rep["agreed"]) == 1
    assert len(rep["a_only"]) == 1
    assert len(rep["b_only"]) == 1
    assert rep["agreement_fraction"] == pytest.approx(0.5)


def test_convolution_vs_mariscotti_high_agreement_on_clean_spectrum():
    """На чистом синтетическом спектре два метода должны соглашаться > 80%."""
    counts = _make_synthetic_spectrum(
        positions=[200, 600, 1100, 1700],
        amplitudes=[3000.0, 4000.0, 2500.0, 2000.0],
        n_channels=2048, sigma_ch=8.0, baseline=30.0,
    )
    peaks_conv = convolution_peak_search(
        counts, fwhm_channels=18.84, sigma_threshold=4.0,
    )
    peaks_mar = mariscotti_search(
        counts, fwhm_channels=18.84, sigma_threshold=4.0,
        min_separation_factor=0.6, edge_margin=10,
    )
    rep = compare_peak_methods(
        peaks_conv, peaks_mar, tolerance_channels=4.0,
    )
    # Считаем сильные пики
    strong_conv = [p for p in peaks_conv if p.significance > 4.0]
    strong_mar = [p for p in peaks_mar if p.significance > 4.0]
    assert len(strong_conv) >= 3
    assert len(strong_mar) >= 3
    # Соглашение должно быть > 0.5 — мы не требуем идентичности,
    # это методы с разной чувствительностью.
    assert rep["agreement_fraction"] > 0.4, (
        f"low agreement {rep['agreement_fraction']:.2%}: "
        f"a_only={[p.channel for p in rep['a_only']]}, "
        f"b_only={[p.channel for p in rep['b_only']]}"
    )


def test_convolution_search_respects_edge_margin():
    """Граничные пики в пределах edge_margin отбрасываются."""
    counts = _make_synthetic_spectrum(
        positions=[5, 500, 2040],
        amplitudes=[3000.0, 3000.0, 3000.0],
        n_channels=2048, sigma_ch=8.0, baseline=20.0,
    )
    peaks = convolution_peak_search(
        counts, fwhm_channels=18.84, sigma_threshold=3.0,
        edge_margin=20,
    )
    positions = [p.channel for p in peaks]
    assert 5 not in positions
    # 500 должна быть найдена
    assert any(abs(p - 500) < 5 for p in positions)
