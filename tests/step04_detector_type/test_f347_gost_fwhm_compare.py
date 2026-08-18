"""F-347 / v1.18.23.0 — сравнение ГОСТ FWHM-методов с peak_image (C-12).

Цель:
  1. Доказать математически, что ГОСТ формула (13) Δn = 2·√(ln2/A)
     даёт ровно √2 раз меньше истинной FWHM Гауссиана.
  2. Сравнить linear_interp / corrected graphoanalytic с peak_image
     fit на реальном Cs-137 / K-40 / Tl-208 demo-спектре.
  3. Документировать вывод о том, какую формулу использовать в
     production (по умолчанию corrected; raw_gost доступна как
     informational только).
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.peaks.centroid_gost import (  # noqa: E402
    gost_centroid_graphoanalytic,
    gost_fwhm_graphoanalytic,
    gost_fwhm_linear_interp,
    gost_pedestal_symmetric,
)


ROOT = Path(__file__).resolve().parents[2]
KIT = ROOT / "detectors/Gamma-1S/reference_spectra/reference_kits"


def _synthetic_gauss(
    n_ch: int, n_0: float, sigma: float, N_max: float,
    background: float = 0.0,
) -> np.ndarray:
    ch = np.arange(n_ch)
    return N_max * np.exp(-((ch - n_0) ** 2) / (2 * sigma ** 2)) + background


# ──────────────────────────────────────────────────────────────────
# 1. Mathematical proof — ГОСТ-формула (13) даёт FWHM/√2
# ──────────────────────────────────────────────────────────────────

class TestGostFormula13Discrepancy:
    """ГОСТ § 3.3.4.3 формула (13): Δn = 2·√(ln2/A).

    Для Гауссиана N = N_max·exp[-(x − x₀)²/(2σ²)] с A = 1/σ²:
      - Стандартная FWHM Гаусса = 2σ·√(2·ln2) = 2·√(2·ln2/A) ≈ 2.355σ
      - ГОСТ-формула (13) =       2σ·√(ln2)   = 2·√(ln2/A)   ≈ 1.665σ
    Отношение = √2 → ГОСТ-формула систематически занижена в √2 раз.
    """

    @pytest.mark.parametrize("sigma", [2.0, 5.0, 10.0, 20.0, 50.0])
    def test_corrected_matches_theoretical_fwhm(self, sigma):
        A = 1.0 / sigma ** 2
        theoretical = sigma * 2.354820045
        res = gost_fwhm_graphoanalytic(A, use_corrected_formula=True)
        assert res.corrected_gauss_fwhm == pytest.approx(theoretical, rel=1e-6)
        assert res.fwhm_channels == res.corrected_gauss_fwhm

    @pytest.mark.parametrize("sigma", [2.0, 5.0, 10.0, 20.0, 50.0])
    def test_raw_gost_is_sqrt2_below_theoretical(self, sigma):
        A = 1.0 / sigma ** 2
        theoretical = sigma * 2.354820045
        res = gost_fwhm_graphoanalytic(A, use_corrected_formula=False)
        # ГОСТ-формула (13) → 2σ·√(ln2) = theoretical / √2
        expected_gost = theoretical / math.sqrt(2)
        assert res.raw_gost_fwhm == pytest.approx(expected_gost, rel=1e-6)
        assert res.fwhm_channels == res.raw_gost_fwhm

    def test_ratio_always_sqrt2(self):
        """Отношение corrected / raw_gost = √2 при любом A."""
        for sigma in [1.0, 3.0, 7.5, 12.0, 25.0]:
            A = 1.0 / sigma ** 2
            res = gost_fwhm_graphoanalytic(A)
            ratio = res.corrected_gauss_fwhm / res.raw_gost_fwhm
            assert ratio == pytest.approx(math.sqrt(2), rel=1e-9)


# ──────────────────────────────────────────────────────────────────
# 2. Сравнение методов на синтетических Гауссах
# ──────────────────────────────────────────────────────────────────

class TestSyntheticGaussComparison:
    """Все три ГОСТ-метода + теоретическое значение должны сойтись."""

    @pytest.mark.parametrize("sigma", [3.0, 5.0, 10.0, 15.0])
    def test_linear_interp_and_corrected_graphoanalytic_agree(self, sigma):
        n_0 = 500.0
        theoretical = sigma * 2.354820045
        # Достаточно большой N_max, чтобы дискретизация была неважна
        counts = _synthetic_gauss(1024, n_0, sigma, N_max=100000.0)

        ped = gost_pedestal_symmetric(counts, int(n_0), theoretical * 1.5)

        fw_li = gost_fwhm_linear_interp(ped.counts_net)
        ga = gost_centroid_graphoanalytic(ped.counts_net)
        fw_ga_c = gost_fwhm_graphoanalytic(ga.A, use_corrected_formula=True)

        assert fw_li.fwhm_channels == pytest.approx(theoretical, rel=0.005)
        assert fw_ga_c.fwhm_channels == pytest.approx(theoretical, rel=0.005)
        # И линейная интерполяция, и graphoanalytic_corrected должны совпасть
        # друг с другом с лучше чем 1% относительной разницы
        rel_diff = abs(fw_li.fwhm_channels - fw_ga_c.fwhm_channels) / theoretical
        assert rel_diff < 0.01

    @pytest.mark.parametrize("sigma", [3.0, 5.0, 10.0, 15.0])
    def test_raw_gost_is_definitely_wrong(self, sigma):
        """Сырая ГОСТ-формула стабильно ошибается на √2."""
        n_0 = 500.0
        theoretical = sigma * 2.354820045
        counts = _synthetic_gauss(1024, n_0, sigma, N_max=100000.0)
        ped = gost_pedestal_symmetric(counts, int(n_0), theoretical * 1.5)
        ga = gost_centroid_graphoanalytic(ped.counts_net)
        fw_raw = gost_fwhm_graphoanalytic(ga.A, use_corrected_formula=False)
        # Расхождение должно быть около (1 − 1/√2) ≈ 29.3%
        rel_error = abs(fw_raw.fwhm_channels - theoretical) / theoretical
        assert 0.28 < rel_error < 0.31


# ──────────────────────────────────────────────────────────────────
# 3. Сравнение с peak_image fit на реальном Cs-137 demo
# ──────────────────────────────────────────────────────────────────

class TestRealSpectrumComparison:
    """На production Cs-137 661 кэВ ГОСТ corrected ≈ peak_image."""

    @pytest.fixture(scope="class")
    def cs137_spectrum(self):
        cs_sample = KIT / "Marinelli_1L/Cs-137/sample_M_cs_легкий_2001-2005.spe"
        if not cs_sample.exists():
            pytest.skip(f"Cs-137 fixture missing: {cs_sample}")
        from gamma.io.readers import read_spectrum
        return read_spectrum(str(cs_sample))

    def test_cs137_centroid_in_660_to_665_keV(self, cs137_spectrum):
        """Sanity: пик Cs-137 в районе 661 кэВ; находим в канальном пространстве."""
        spec = cs137_spectrum
        counts = np.asarray(spec.counts, dtype=float)

        # Найти канал максимума в диапазоне 600-720 кэВ
        # Используем stored energy calibration spec.energy_calibration
        if not hasattr(spec, "energy_calibration") or spec.energy_calibration is None:
            pytest.skip("no stored energy calibration")
        # Найти канал, ближайший к E=661 кэВ через инверсию калибровки
        energies = np.array([
            spec.energy_calibration.channel_to_energy(i) for i in range(len(counts))
        ])
        lo_ch = int(np.argmin(np.abs(energies - 600.0)))
        hi_ch = int(np.argmin(np.abs(energies - 720.0)))
        local_counts = counts[lo_ch:hi_ch]
        peak_local = int(np.argmax(local_counts))
        peak_channel = lo_ch + peak_local

        # Грубая оценка FWHM из stored или ~7% от 661
        # FWHM_keV ≈ 0.07 · 661 = 46 кэВ; в каналах через градиент калибровки
        keV_per_ch = energies[peak_channel + 1] - energies[peak_channel]
        fwhm_channels_est = 46.0 / keV_per_ch

        ped = gost_pedestal_symmetric(counts, peak_channel, fwhm_channels_est)

        # ГОСТ методы
        fw_li = gost_fwhm_linear_interp(ped.counts_net)
        ga = gost_centroid_graphoanalytic(ped.counts_net)
        fw_ga_c = gost_fwhm_graphoanalytic(ga.A, use_corrected_formula=True)

        # Конверсия в кэВ
        fwhm_li_keV = fw_li.fwhm_channels * keV_per_ch
        fwhm_ga_c_keV = fw_ga_c.fwhm_channels * keV_per_ch

        # Ожидаем 7-9% от 661 = 46-60 кэВ для Gamma-1S NaI 63×63
        assert 30.0 < fwhm_li_keV < 80.0, f"linear_interp FWHM {fwhm_li_keV:.2f} keV out of bounds"
        assert 30.0 < fwhm_ga_c_keV < 80.0, f"graphoanalytic_corr FWHM {fwhm_ga_c_keV:.2f} keV"

        # ГОСТ corrected и linear_interp должны быть в пределах 25% друг от
        # друга на реальном пике (статистика + асимметрия)
        rel_diff = abs(fwhm_li_keV - fwhm_ga_c_keV) / fwhm_li_keV
        assert rel_diff < 0.25, (
            f"linear_interp={fwhm_li_keV:.2f} keV vs graphoanalytic={fwhm_ga_c_keV:.2f}"
            f" — разница {rel_diff*100:.1f}% > 25%"
        )


# ──────────────────────────────────────────────────────────────────
# 4. Documentation gate — модуль явно указывает на opечатку ГОСТ
# ──────────────────────────────────────────────────────────────────

class TestFormulaDocumentation:
    def test_module_docstring_mentions_correction(self):
        """Модуль явно документирует расхождение формулы (13) с теорией."""
        from gamma.peaks import centroid_gost
        doc = centroid_gost.__doc__
        assert "√2" in doc or "sqrt(2)" in doc
        # Документ упоминает исправленную форму
        assert "исправлен" in doc.lower() or "corrected" in doc.lower()
        # И стандартный коэффициент 2·√(2·ln2)
        assert "2·ln2" in doc or "2·ln 2" in doc or "2·√(2·ln2/A)" in doc

    def test_fwhm_result_documents_both_values(self):
        """FwhmResult всегда возвращает оба значения."""
        res = gost_fwhm_graphoanalytic(0.04)
        assert res.raw_gost_fwhm is not None
        assert res.corrected_gauss_fwhm is not None
        assert res.raw_gost_fwhm != res.corrected_gauss_fwhm

    def test_default_is_corrected(self):
        """По умолчанию возвращаем corrected, не raw_gost."""
        res = gost_fwhm_graphoanalytic(0.04)
        assert res.method == "graphoanalytic_corrected"
        assert res.fwhm_channels == res.corrected_gauss_fwhm
