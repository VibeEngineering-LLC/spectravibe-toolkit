"""F-133 / v1.17.7 — Per-line «ступенька под пиком» (ГОСТ форма).

Контракты, зафиксированные навсегда:
  • ВСЯКАЯ форма пика γ-спектрометрии на NaI = Гаусс + tail + per-line step
  • h_step (доля step от peak height) по умолчанию 0.03 для NaI (LSRM §8.4.4)
  • Центроиды пиков СТРОГО по паспортным энергиям (нет глобального dE сдвига)
  • В мультиплетах интенсивности связаны через библиотечные I_γ (F-117)
  • Глобальная β_step колонка подавляется, когда per-line step активен —
    иначе ступенька удваивается
  • Площадь Гаусс+tail = A_k (unit-area), step НЕ входит в A_k
  • На Th-232 demo M1 closure после F-133 близко к 0 (vs −78% до фикса)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.peaks.coupled_multiplet import (
    ComponentSpec, H_STEP_DEFAULT_NAI,
    _peak_image_with_step, _peak_image_unit_area,
    _gaussian_unit_area, _peak_image_normalisation,
    coupled_intensity_fit,
)


# ─── unit tests ──────────────────────────────────────────────────────

def test_h_step_default_is_0_03():
    """ГОСТ / LSRM §8.4.4 → h_step ≈ 0.03 для NaI 63×63."""
    assert H_STEP_DEFAULT_NAI == 0.03


def test_step_is_zero_above_peak():
    """Справа от пика (E >> E0) step → 0 (erfc → 0). Абсолютное
    превышение vs no_step должно быть пренебрежимо мало по сравнению
    с peak_height в центре."""
    sigma = 30.0
    E0 = 911.2
    T = 0.7
    h_step = 0.03
    E = np.array([E0 + 5 * sigma])    # 5σ справа
    with_step = _peak_image_with_step(E, E0, sigma, T, h_step)
    no_step = _peak_image_unit_area(E, E0, sigma, T)
    # peak_height в центре пика
    norm = _peak_image_normalisation(sigma, T)
    peak_height = 1.0 / norm
    # step справа от пика должен быть < 1e-5 от peak_height
    excess = abs(with_step[0] - no_step[0])
    assert excess < 1e-5 * peak_height


def test_step_is_finite_below_peak():
    """Слева от пика (E << E0) step → h_step·peak_height (постоянный фон)."""
    sigma = 30.0
    E0 = 911.2
    T = 0.7
    h_step = 0.03
    # 5σ слева — но это в tail region, не в чистом Гауссе.
    # Для чистоты теста возьмём 10σ слева (далеко за tail).
    E = np.array([E0 - 10 * sigma])
    with_step = _peak_image_with_step(E, E0, sigma, T, h_step)
    norm = _peak_image_normalisation(sigma, T)
    expected_step = h_step * (1.0 / norm)
    # Tail в этой точке практически 0
    assert abs(with_step[0] - expected_step) < 0.1 * expected_step


def test_step_is_half_at_centroid():
    """В центре пика (E = E0) step = h_step·peak_height/2 (erfc(0)=1)."""
    sigma = 30.0
    E0 = 911.2
    T = 0.7
    h_step = 0.03
    E = np.array([E0])
    with_step = _peak_image_with_step(E, E0, sigma, T, h_step)
    no_step = _peak_image_unit_area(E, E0, sigma, T)
    # Превышение — это step-вклад при E=E0 = h_step·peak_height·0.5
    norm = _peak_image_normalisation(sigma, T)
    expected_excess = h_step * (1.0 / norm) * 0.5
    actual_excess = with_step[0] - no_step[0]
    assert abs(actual_excess - expected_excess) < 1e-6 * expected_excess


def test_h_step_zero_equivalent_to_no_step():
    """h_step=0 → возвращает чистый Гаусс+tail (back-compat с v1.17.6)."""
    sigma = 30.0
    E0 = 911.2
    T = 0.7
    E = np.linspace(E0 - 100, E0 + 100, 50)
    with_step = _peak_image_with_step(E, E0, sigma, T, 0.0)
    no_step = _peak_image_unit_area(E, E0, sigma, T)
    assert np.allclose(with_step, no_step)


def test_coupled_intensity_fit_accepts_h_step():
    """coupled_intensity_fit принимает h_step kwarg без падения."""
    E = np.linspace(900.0, 1000.0, 200)
    counts = np.full_like(E, 50.0)
    res = coupled_intensity_fit(
        E, counts, [ComponentSpec("Test", 950.0, 100.0)],
        lambda x: 10.0,
        use_peak_image=True, h_step=0.03,
    )
    # basis_label должен содержать "+step=0.030"
    assert "step=0.030" in res.notes


def test_global_beta_step_suppressed_when_per_line_step_active():
    """При h_step > 0 глобальная β_step колонка не добавляется в матрицу."""
    E = np.linspace(900.0, 1000.0, 200)
    counts = np.full_like(E, 50.0)
    res_no_step = coupled_intensity_fit(
        E, counts, [ComponentSpec("Test", 950.0, 100.0)],
        lambda x: 10.0, use_peak_image=True, h_step=0.0,
        continuum="step_linear",
    )
    res_with_step = coupled_intensity_fit(
        E, counts, [ComponentSpec("Test", 950.0, 100.0)],
        lambda x: 10.0, use_peak_image=True, h_step=0.03,
        continuum="step_linear",
    )
    # continuum_params: [β₀, β₁, β_step]   vs   [β₀, β₁]
    assert len(res_no_step.continuum_params) == 3
    assert len(res_with_step.continuum_params) == 2


def test_centroids_strictly_from_library_E():
    """F-133 контракт: центроиды строго по паспортным E_keV.
    После nonlinear refine компоненты ДОЛЖНЫ остаться в своих E_keV.
    """
    E = np.linspace(900.0, 1000.0, 200)
    rng = np.random.default_rng(0)
    sigma = 12.0
    gauss = (np.exp(-((E - 952.0) / sigma) ** 2 * 0.5)
             / (sigma * np.sqrt(2 * np.pi))) * 50000.0
    counts = rng.poisson(np.maximum(100.0 + gauss, 1.0)).astype(float)
    res = coupled_intensity_fit(
        E, counts, [ComponentSpec("Test", 950.0, 100.0)],
        lambda x: 10.0,
        use_peak_image=True, h_step=0.03,
        nonlinear_refine=True,
    )
    # E_keV компоненты не сдвигается ни на что
    assert res.components[0].E_keV == 950.0


# ─── E2E на Th-232 fixture ────────────────────────────────────────────

TH232_FIXTURE = (Path(__file__).parent.parent.parent / "detectors" / "Gamma-1S"
                 / "reference_spectra"
                 / "archive"
                 / "Th232_420-7-17_Маринелли_0cm.spe")


def test_th232_m1_closure_near_zero():
    """E2E: M1 closure после F-133 должна быть |Δ| < 5 % (vs −78.8% до фикса).
    Эталон v1.17.2: −0.9 %."""
    if not TH232_FIXTURE.exists():
        pytest.skip("fixture missing")
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    r = analyze_lsrm_spe(str(TH232_FIXTURE), complete_workflow=True)
    m1 = None
    for d in r.deconvolution_results or []:
        if d.notes and "кластер M1" in d.notes:
            m1 = d
            break
    assert m1 is not None
    # closure_pct запихнут в notes как "closure Δ=-0.61%". Распарсим.
    import re
    m = re.search(r"closure Δ=([\-+]?\d+\.\d+)%", m1.notes)
    assert m, f"closure не найдено в notes: {m1.notes}"
    closure = float(m.group(1))
    assert abs(closure) < 5.0, (
        f"M1 closure |Δ|={abs(closure):.2f}% должно быть < 5 % после F-133 "
        f"(до фикса было ≈ −78.8%; эталон v1.17.2: −0.9%)"
    )


def test_th232_m1_areas_match_v1_17_2_reference():
    """E2E: M1 площади должны попасть в полосу ±5 % от эталона v1.17.2.

    Эталон: Ac-228 911=116996, 964.77=22628, 968.97=71649.
    """
    if not TH232_FIXTURE.exists():
        pytest.skip("fixture missing")
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    r = analyze_lsrm_spe(str(TH232_FIXTURE), complete_workflow=True)
    m1 = None
    for d in r.deconvolution_results or []:
        if d.notes and "кластер M1" in d.notes:
            m1 = d
            break
    assert m1 is not None
    expected = {
        911.2: 116996,
        964.77: 22628,
        968.97: 71649,
    }
    for comp, area in zip(m1.components, m1.areas):
        for E_ref, A_ref in expected.items():
            if abs(comp.line_E_keV - E_ref) < 0.1:
                rel_err = abs(area - A_ref) / A_ref
                assert rel_err < 0.05, (
                    f"Ac-228 @ {E_ref}: area={area:.0f}, "
                    f"эталон={A_ref}, отклонение {rel_err*100:.1f}% > 5 %"
                )
