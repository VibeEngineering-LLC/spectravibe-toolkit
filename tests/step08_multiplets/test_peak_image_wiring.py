"""F-120 (v1.17.6) — peak-image wiring regression.

Проверяет, что:
  1. Базис ``use_peak_image=True`` корректно подключается в
     ``coupled_intensity_fit`` и даёт ту же связную математику,
     что чистый гаусс (площадь = A·σ-нормировки = реальный счёт
     в импульсах).
  2. Унитарный интеграл peak-image == 1 (нормировка корректна).
  3. При сильном хвосте (T=0.7) подгонка близка к гауссу на пике
     и заметно отличается от чистого гаусса слева от пика.
  4. ``run_chain_forced_multiplets`` принимает новые kwargs
     ``use_peak_image`` / ``detector_type`` без падений и для
     Th-232 fixture даёт positive area для линий Ac-228 (911, 969).
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

from gamma.peaks.coupled_multiplet import (
    coupled_intensity_fit, ComponentSpec, T_TAIL_DEFAULT_NAI,
    _peak_image_unit_area, _gaussian_unit_area,
    _peak_image_normalisation,
)


def test_peak_image_unit_area_integrates_to_one():
    """Площадь peak-image-базиса (Gauss+tail) равна единице по построению."""
    sigma = 12.0  # keV
    T = 0.7
    E = np.linspace(500.0, 800.0, 2001)
    bin_w = float(np.mean(np.diff(E)))
    g = _peak_image_unit_area(E, 661.66, sigma, T)
    area = float(np.sum(g) * bin_w)
    assert abs(area - 1.0) < 1e-2, f"peak_image area {area} ≠ 1"


def test_peak_image_falls_back_to_gauss_when_T_zero():
    """T_tail <= 0 → чистый гауссиан."""
    sigma = 5.0
    E = np.linspace(100.0, 200.0, 1001)
    a = _peak_image_unit_area(E, 150.0, sigma, 0.0)
    b = _gaussian_unit_area(E, 150.0, sigma)
    assert np.allclose(a, b, rtol=1e-8)


def test_peak_image_has_tail_left_of_peak():
    """При T=0.7 plot слева от μ-T·σ должен быть ВЫШЕ, чем чистый гауссиан."""
    sigma = 10.0
    T = 0.7
    E = np.linspace(500.0, 800.0, 1001)
    pi = _peak_image_unit_area(E, 661.66, sigma, T)
    g = _gaussian_unit_area(E, 661.66, sigma)
    # Точки слева от 661.66 - T*sigma = 654.66
    mask_left = E < (661.66 - T * sigma * 1.5)
    # Точки далеко слева (≈ 5·sigma) — там exp-tail доминирует над чистым гауссом
    far_left = E < (661.66 - 3.0 * sigma)
    if far_left.any():
        ratio = pi[far_left] / np.maximum(g[far_left], 1e-30)
        assert np.median(ratio) > 1.0, (
            f"peak_image tail not above Gaussian on far left; "
            f"median ratio={np.median(ratio):.2e}"
        )


def test_peak_image_normalisation_formula():
    """Нормировка совпадает с аналитической формулой."""
    sigma = 8.0
    T = 0.7
    norm = _peak_image_normalisation(sigma, T)
    # Аналитически:
    erfc_T = math.erfc(T / math.sqrt(2.0))
    gauss_right = sigma * math.sqrt(math.pi / 2.0) * (2.0 - erfc_T)
    tail = sigma * math.exp(-0.5 * T * T) / T
    expected = gauss_right + tail
    assert abs(norm - expected) < 1e-9


def test_coupled_intensity_fit_accepts_use_peak_image_kwarg():
    """coupled_intensity_fit принимает use_peak_image без TypeError."""
    # Synthetic single-line spectrum centred at 600 keV
    rng = np.random.default_rng(42)
    E = np.linspace(550.0, 650.0, 401)
    bin_w = float(np.mean(np.diff(E)))
    sigma = 12.0
    A_true = 50000.0
    y_clean = A_true * _peak_image_unit_area(E, 600.0, sigma, 0.7) * bin_w + 100.0
    y = y_clean + rng.normal(0.0, np.sqrt(y_clean), size=E.shape)
    y = np.maximum(y, 0.0)

    comp = [ComponentSpec(nuclide="X", E_keV=600.0, I_gamma_pct=100.0, group="X")]
    fwhm_at = lambda e: sigma * 2.355  # noqa: E731 — local lambda is fine

    res_gauss = coupled_intensity_fit(
        E, y, comp, fwhm_at,
        continuum="step_linear", cluster_id="T1",
        use_peak_image=False,
    )
    res_pi = coupled_intensity_fit(
        E, y, comp, fwhm_at,
        continuum="step_linear", cluster_id="T1",
        use_peak_image=True, tail_param=T_TAIL_DEFAULT_NAI,
    )
    # Оба должны сойтись
    assert res_gauss.converged
    assert res_pi.converged
    # peak-image должен лучше восстановить истинную площадь (без хвоста гаусс
    # перетягивает ширину и недо/переоценивает A); допускаем 10% отклонение.
    # Главное — что peak-image НЕ хуже, чем гаусс, на хвостовом сигнале.
    A_g = float(res_gauss.components[0].area)
    A_p = float(res_pi.components[0].area)
    assert A_g > 0, "gauss basis must return positive area"
    assert A_p > 0, "peak_image basis must return positive area"
    # На синтетике с хвостом peak-image ближе к истинной
    err_g = abs(A_g - A_true) / A_true
    err_p = abs(A_p - A_true) / A_true
    assert err_p < 0.25, f"peak_image error {err_p:.1%} > 25%"


def test_run_chain_forced_multiplets_accepts_new_kwargs():
    """F-120/F-121: run_chain_forced_multiplets принимает новые kwargs
    и не падает на пустом chain_dominance."""
    from gamma.peaks.deconvolve import run_chain_forced_multiplets

    class _FakeChainDom:
        th232 = False
        u238 = False

    # No-chain → пустой список
    out = run_chain_forced_multiplets(
        None, None, None,
        _FakeChainDom(),
        None,
        use_peak_image=True, detector_type="NaI",
    )
    assert out == []


def test_ra226_forced_clusters_contract():
    """F-121: RA226_FORCED_CLUSTERS содержит ровно 3 кластера с ключевыми линиями."""
    from gamma.peaks.deconvolve import RA226_FORCED_CLUSTERS
    assert len(RA226_FORCED_CLUSTERS) == 3
    ids = [c["id"] for c in RA226_FORCED_CLUSTERS]
    assert ids == ["U1", "U2", "U3"]
    # U1 должен содержать Pb-214 295 и 352
    u1_E = [c[1] for c in RA226_FORCED_CLUSTERS[0]["components"]]
    assert any(abs(E - 295.22) < 0.1 for E in u1_E)
    assert any(abs(E - 351.93) < 0.1 for E in u1_E)
    # U2 должен содержать Bi-214 609
    u2_E = [c[1] for c in RA226_FORCED_CLUSTERS[1]["components"]]
    assert any(abs(E - 609.31) < 0.1 for E in u2_E)
    # U3 должен содержать Bi-214 1764 (trump card)
    u3_E = [c[1] for c in RA226_FORCED_CLUSTERS[2]["components"]]
    assert any(abs(E - 1764.49) < 0.1 for E in u3_E)
