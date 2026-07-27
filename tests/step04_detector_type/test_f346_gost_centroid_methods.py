"""F-346 / v1.18.23.0 — ГОСТ 26874-86 § 3.3.2-3.3.4 методы (unit-тесты).

Покрытие:
  * § 3.3.2 — pedestal (symmetric / asymmetric / threshold 2%)
  * § 3.3.3.1 — graphical centroid
  * § 3.3.3.2 — weighted_mean (формула 5)
  * § 3.3.3.3 — graphoanalytic (формулы 6-10)
  * § 3.3.4.2 — linear_interp FWHM (формула 12)
  * § 3.3.4.3 — graphoanalytic FWHM (формула 13 + corrected)
  * ROI gate (peak ± 0.5·FWHM)
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.peaks.centroid_gost import (  # noqa: E402
    PEDESTAL_THRESHOLD_FRAC,
    gost_centroid_graphical,
    gost_centroid_graphoanalytic,
    gost_centroid_weighted_mean,
    gost_fwhm_graphoanalytic,
    gost_fwhm_linear_interp,
    gost_pedestal_asymmetric,
    gost_pedestal_symmetric,
    gost_roi_from_fwhm,
    gost_select_pedestal_method,
    should_subtract_pedestal,
)


# ──────────────────────────────────────────────────────────────────
# Synthetic Gaussian helper
# ──────────────────────────────────────────────────────────────────

def _gauss_peak(
    n_ch: int = 1024,
    n_0: float = 500.0,
    sigma: float = 5.0,
    N_max: float = 10000.0,
    background: float = 0.0,
    slope: float = 0.0,
    poisson_seed: int | None = None,
) -> np.ndarray:
    """Синтетический Gauss-пик + опциональный наклонный фон + Пуассон-шум."""
    ch = np.arange(n_ch)
    peak = N_max * np.exp(-((ch - n_0) ** 2) / (2 * sigma ** 2))
    bg = background + slope * (ch - n_0)
    counts = peak + bg
    if poisson_seed is not None:
        rng = np.random.default_rng(poisson_seed)
        counts = rng.poisson(np.clip(counts, 0, None)).astype(float)
    return counts


# ──────────────────────────────────────────────────────────────────
# § 3.3.2 — Pedestal
# ──────────────────────────────────────────────────────────────────

class TestPedestalGate:
    def test_no_subtract_when_bg_below_2pct(self):
        # N_max = peak + bg = 10000 + 50 = 10050; фон = 50 → 50/10050 ≈ 0.5% < 2%
        counts = _gauss_peak(N_max=10000.0, background=50.0, sigma=5.0)
        do_sub, n_max, n_bg = should_subtract_pedestal(
            counts, peak_channel=500, fwhm_channels=11.77,
        )
        assert do_sub is False
        # n_max включает фон: 10000 + 50 = 10050
        assert n_max == pytest.approx(10050.0, rel=0.001)
        assert n_bg == pytest.approx(50.0, abs=2.0)

    def test_subtract_when_bg_above_2pct(self):
        # N_max = 10000, фон = 500 → 5% > 2%
        counts = _gauss_peak(N_max=10000.0, background=500.0, sigma=5.0)
        do_sub, _, _ = should_subtract_pedestal(
            counts, peak_channel=500, fwhm_channels=11.77,
        )
        assert do_sub is True

    def test_threshold_constant(self):
        assert PEDESTAL_THRESHOLD_FRAC == 0.02


class TestPedestalSymmetric:
    def test_horizontal_line_through_average(self):
        counts = _gauss_peak(N_max=10000.0, background=200.0, slope=0.0, sigma=5.0)
        res = gost_pedestal_symmetric(counts, peak_channel=500, fwhm_channels=11.77)
        assert res.method == "symmetric"
        # Горизонтальная линия — все значения близки к 200.
        assert res.pedestal.mean() == pytest.approx(200.0, abs=2.0)
        assert res.pedestal.std() < 1.0
        # counts_net ≈ pure peak в центре ROI.
        peak_idx_in_roi = 500 - res.roi_lo
        assert res.counts_net[peak_idx_in_roi] == pytest.approx(10000.0, rel=0.01)


class TestPedestalAsymmetric:
    def test_linear_slope_recovered(self):
        # Линейный наклон: background = 200 + 0.5·(ch − 500)
        counts = _gauss_peak(
            N_max=10000.0, background=200.0, slope=0.5, sigma=5.0,
        )
        res = gost_pedestal_asymmetric(
            counts, peak_channel=500, fwhm_channels=11.77,
            pedestal_gap_fwhm=4.0, pedestal_window=9,
        )
        assert res.method == "asymmetric"
        # gap ≥ 4·FWHM
        assert res.pedestal_gap_channels >= int(round(4 * 11.77))
        # N_l_avg в районе background(450) = 200 + 0.5·(450 − 500) = 175
        # N_h_avg в районе background(550) = 200 + 0.5·(550 − 500) = 225
        # (точные значения зависят от gap; проверяем порядок)
        assert res.N_l_avg < res.N_h_avg
        assert res.N_l_avg == pytest.approx(175, abs=5)
        assert res.N_h_avg == pytest.approx(225, abs=5)

    def test_pedestal_subtraction_recovers_peak(self):
        # После вычитания фона counts_net на пике ≈ N_max.
        counts = _gauss_peak(
            N_max=10000.0, background=200.0, slope=0.5, sigma=5.0,
        )
        res = gost_pedestal_asymmetric(
            counts, peak_channel=500, fwhm_channels=11.77,
        )
        peak_in_roi = 500 - res.roi_lo
        # Допускаем ≤ 1% relative error
        assert res.counts_net[peak_in_roi] == pytest.approx(10000.0, rel=0.01)


class TestPedestalAutoSelect:
    def test_auto_picks_symmetric_for_flat_bg(self):
        # bg=500 → 500/10500 ≈ 4.8% > 2% threshold → не отсеется
        counts = _gauss_peak(N_max=10000.0, background=500.0, slope=0.0)
        res = gost_select_pedestal_method(
            counts, peak_channel=500, fwhm_channels=11.77,
        )
        assert res.method == "symmetric"

    def test_auto_picks_asymmetric_for_steep_slope(self):
        # Сильный наклон: 50 counts/canal → N_l ≈ 200 − 50·47 = −2150, N_h ≈ 200 + 50·47 = 2550
        # diff = 4700, σ_l = √(N̅/9) ≈ 5 → diff/σ_l >> 1.0
        counts = _gauss_peak(N_max=50000.0, background=5000.0, slope=50.0)
        res = gost_select_pedestal_method(
            counts, peak_channel=500, fwhm_channels=11.77,
        )
        assert res.method == "asymmetric"

    def test_auto_skips_when_bg_below_threshold(self):
        counts = _gauss_peak(N_max=10000.0, background=50.0, slope=0.0)
        res = gost_select_pedestal_method(
            counts, peak_channel=500, fwhm_channels=11.77,
        )
        assert res.method == "none"
        assert "< 0.02" in res.skipped_reason


# ──────────────────────────────────────────────────────────────────
# § 3.3.3 — Centroid
# ──────────────────────────────────────────────────────────────────

class TestCentroidGraphical:
    def test_centroid_on_pure_gauss(self):
        counts = _gauss_peak(n_0=500.0, sigma=5.0, N_max=10000.0)
        ped = gost_pedestal_symmetric(counts, 500, 11.77)
        res = gost_centroid_graphical(ped.counts_net, channel_offset=ped.roi_lo)
        assert res.method == "graphical"
        assert res.n_c == pytest.approx(500.0, abs=0.1)

    def test_centroid_on_non_integer_position(self):
        counts = _gauss_peak(n_0=500.5, sigma=5.0, N_max=10000.0)
        ped = gost_pedestal_symmetric(counts, 500, 11.77)
        res = gost_centroid_graphical(ped.counts_net, channel_offset=ped.roi_lo)
        # Графический метод должен дать ~500.5 (с дискретизационной погрешностью)
        assert res.n_c == pytest.approx(500.5, abs=0.15)


class TestCentroidWeightedMean:
    def test_weighted_mean_on_pure_gauss(self):
        counts = _gauss_peak(n_0=500.0, sigma=5.0, N_max=10000.0)
        ped = gost_pedestal_symmetric(counts, 500, 11.77)
        res = gost_centroid_weighted_mean(ped.counts_net, channel_offset=ped.roi_lo)
        assert res.method == "weighted_mean"
        assert res.n_c == pytest.approx(500.0, abs=0.05)
        assert res.sigma_n_c > 0
        # n_points_used — только выше полувысоты, для σ=5 это ~2·1.177·5 ≈ 12 каналов
        assert 8 <= res.n_points_used <= 16

    def test_weighted_mean_on_shifted_gauss(self):
        # Дискретизационная погрешность weighted_mean на нецелочисленной
        # позиции n_0=512.3 (σ=5 ch) — ≤ 0.2 канала, фактически ~0.12.
        counts = _gauss_peak(n_0=512.3, sigma=5.0, N_max=10000.0)
        ped = gost_pedestal_symmetric(counts, 512, 11.77)
        res = gost_centroid_weighted_mean(ped.counts_net, channel_offset=ped.roi_lo)
        assert res.n_c == pytest.approx(512.3, abs=0.2)


class TestCentroidGraphoanalytic:
    def test_graphoanalytic_on_pure_gauss(self):
        counts = _gauss_peak(n_0=500.0, sigma=5.0, N_max=10000.0)
        ped = gost_pedestal_symmetric(counts, 500, 11.77)
        res = gost_centroid_graphoanalytic(ped.counts_net, channel_offset=ped.roi_lo)
        assert res.method == "graphoanalytic"
        assert res.n_c == pytest.approx(500.0, abs=0.05)
        # A = 1/σ² = 1/25 = 0.04
        assert res.A == pytest.approx(0.04, rel=0.01)
        # σ восстановленная: 1/√A = 5.0
        recovered_sigma = 1.0 / math.sqrt(res.A)
        assert recovered_sigma == pytest.approx(5.0, rel=0.01)

    def test_graphoanalytic_raises_for_too_few_pairs(self):
        # Очень узкий пик — мало точек выше половины
        counts = np.zeros(100)
        counts[50] = 1000.0
        counts[51] = 100.0
        with pytest.raises(ValueError, match="only|≥"):
            gost_centroid_graphoanalytic(counts)


# ──────────────────────────────────────────────────────────────────
# § 3.3.4 — FWHM
# ──────────────────────────────────────────────────────────────────

class TestFwhmLinearInterp:
    @pytest.mark.parametrize("sigma", [2.0, 5.0, 10.0, 20.0])
    def test_fwhm_matches_theoretical(self, sigma):
        # Theoretical FWHM = σ · 2·√(2·ln2) ≈ σ · 2.3548
        theoretical = sigma * 2.354820045
        counts = _gauss_peak(
            n_0=500.0, sigma=sigma, N_max=10000.0, n_ch=1024,
        )
        ped = gost_pedestal_symmetric(counts, 500, theoretical * 1.1)
        res = gost_fwhm_linear_interp(ped.counts_net)
        assert res.fwhm_channels == pytest.approx(theoretical, rel=0.01)


class TestFwhmGraphoanalytic:
    def test_corrected_formula_matches_theoretical(self):
        sigma = 5.0
        theoretical = sigma * 2.354820045
        # A = 1/σ² = 0.04
        A = 1.0 / sigma ** 2
        res = gost_fwhm_graphoanalytic(A, use_corrected_formula=True)
        assert res.fwhm_channels == pytest.approx(theoretical, rel=0.001)
        assert res.corrected_gauss_fwhm == pytest.approx(theoretical, rel=0.001)

    def test_raw_gost_formula_is_sqrt2_smaller(self):
        sigma = 5.0
        A = 1.0 / sigma ** 2
        res = gost_fwhm_graphoanalytic(A, use_corrected_formula=False)
        theoretical = sigma * 2.354820045
        # ГОСТ-формула (13) даёт theoretical / √2
        assert res.fwhm_channels == pytest.approx(theoretical / math.sqrt(2), rel=0.001)

    def test_ratio_is_sqrt2(self):
        A = 0.04
        res_c = gost_fwhm_graphoanalytic(A, use_corrected_formula=True)
        ratio = res_c.corrected_gauss_fwhm / res_c.raw_gost_fwhm
        assert ratio == pytest.approx(math.sqrt(2), rel=0.0001)


# ──────────────────────────────────────────────────────────────────
# ROI helper
# ──────────────────────────────────────────────────────────────────

class TestRoiSelection:
    def test_roi_around_peak_center(self):
        lo, hi = gost_roi_from_fwhm(
            peak_channel=500, fwhm_channels=10.0, half_fwhm=2.5, n_total=1024,
        )
        assert lo == 475
        assert hi == 526  # +25 +1 exclusive

    def test_roi_clamps_to_spectrum_edges(self):
        lo, hi = gost_roi_from_fwhm(
            peak_channel=5, fwhm_channels=10.0, half_fwhm=2.5, n_total=1024,
        )
        assert lo == 0
        assert hi == 31


# ──────────────────────────────────────────────────────────────────
# Cross-method consistency (все 3 метода центроид сходятся)
# ──────────────────────────────────────────────────────────────────

class TestCrossMethodConsistency:
    @pytest.mark.parametrize("n_0_exact", [500.0, 500.3, 500.7, 499.5])
    def test_three_methods_agree_within_quarter_channel(self, n_0_exact):
        counts = _gauss_peak(n_0=n_0_exact, sigma=5.0, N_max=10000.0)
        ped = gost_pedestal_symmetric(counts, int(round(n_0_exact)), 11.77)
        c_gr = gost_centroid_graphical(ped.counts_net, channel_offset=ped.roi_lo)
        c_wm = gost_centroid_weighted_mean(ped.counts_net, channel_offset=ped.roi_lo)
        c_ga = gost_centroid_graphoanalytic(ped.counts_net, channel_offset=ped.roi_lo)
        # Все три должны сойтись на чистом Гауссе в пределах 0.25 канала
        results = [c_gr.n_c, c_wm.n_c, c_ga.n_c]
        spread = max(results) - min(results)
        assert spread < 0.25, f"methods spread {spread:.4f} ≥ 0.25 ch"
        # И все близки к истинной позиции
        for r in results:
            assert abs(r - n_0_exact) < 0.25
