# -*- coding: utf-8 -*-
"""v1.17.20 — Activity accuracy slice 1 (F-294..F-298)."""
from __future__ import annotations
import math, os, sys
import pytest

SCRIPTS = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts")
)
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)


# ──────────────────────────────────────────────────────────────────
# F-294 — Marinelli self-absorption (T-002 + T-024)
# ──────────────────────────────────────────────────────────────────

def test_F294_self_absorption_unity_when_equal_density():
    from gamma.activity.self_absorption import (
        SelfAbsorptionInputs, self_absorption_factor,
    )
    f = self_absorption_factor(SelfAbsorptionInputs(
        E_keV=662.0, rho_sample_g_cm3=1.0, rho_calib_g_cm3=1.0,
    ))
    assert f == pytest.approx(1.0, abs=1e-6)


def test_F294_self_absorption_less_than_unity_when_sample_denser():
    """Образец плотнее эталона → больше внутреннего самопоглощения
    → f_abs < 1 → A_corrected > A_apparent."""
    from gamma.activity.self_absorption import (
        SelfAbsorptionInputs, self_absorption_factor,
    )
    f = self_absorption_factor(SelfAbsorptionInputs(
        E_keV=100.0, rho_sample_g_cm3=2.0, rho_calib_g_cm3=1.0,
    ))
    assert 0.0 < f < 1.0


def test_F294_self_absorption_greater_than_unity_when_sample_lighter():
    """Образец легче эталона → меньше самопоглощения → f_abs > 1."""
    from gamma.activity.self_absorption import (
        SelfAbsorptionInputs, self_absorption_factor,
    )
    f = self_absorption_factor(SelfAbsorptionInputs(
        E_keV=100.0, rho_sample_g_cm3=0.5, rho_calib_g_cm3=1.0,
    ))
    assert f > 1.0


def test_F294_high_E_less_sensitive_than_low_E():
    """На 1500 keV эффект слабее чем на 100 keV (mu/rho меньше)."""
    from gamma.activity.self_absorption import (
        SelfAbsorptionInputs, self_absorption_factor,
    )
    f_low = self_absorption_factor(SelfAbsorptionInputs(
        E_keV=100.0, rho_sample_g_cm3=2.5, rho_calib_g_cm3=1.0,
    ))
    f_high = self_absorption_factor(SelfAbsorptionInputs(
        E_keV=1500.0, rho_sample_g_cm3=2.5, rho_calib_g_cm3=1.0,
    ))
    # |1 - f| на низкой E больше чем на высокой
    assert abs(1.0 - f_low) > abs(1.0 - f_high)


def test_F294_correct_activity_inverts_factor():
    from gamma.activity.self_absorption import (
        correct_activity_for_self_absorption,
    )
    A_corr = correct_activity_for_self_absorption(100.0, 0.8)
    assert A_corr == pytest.approx(125.0, rel=1e-6)


def test_F294_mu_over_rho_water_decreases_with_E():
    from gamma.activity.self_absorption import mu_over_rho_water
    assert mu_over_rho_water(50.0) > mu_over_rho_water(500.0)
    assert mu_over_rho_water(500.0) > mu_over_rho_water(2000.0)


# ──────────────────────────────────────────────────────────────────
# F-295 — P/T ratio NaI (T-011)
# ──────────────────────────────────────────────────────────────────

def test_F295_pt_ratio_decreases_with_E():
    """P/T уменьшается с E (Compton continuum растёт быстрее FEP)."""
    from gamma.activity.pt_ratio_nai import pt_ratio_nai
    assert pt_ratio_nai(100.0) > pt_ratio_nai(500.0)
    assert pt_ratio_nai(500.0) > pt_ratio_nai(1500.0)


def test_F295_pt_ratio_larger_for_bigger_crystal():
    """4"×4" даёт большее P/T чем Gamma-1S 63 mm."""
    from gamma.activity.pt_ratio_nai import pt_ratio_nai
    pt_small = pt_ratio_nai(1000.0, crystal_diameter_mm=63.0)
    pt_big = pt_ratio_nai(1000.0, crystal_diameter_mm=102.0)
    assert pt_big > pt_small


def test_F295_pt_ratio_gilmore_3in3_anchor():
    """Прямая проверка Gilmore Table 8.4 на 1000 keV для 3"×3"."""
    from gamma.activity.pt_ratio_nai import pt_ratio_nai
    pt = pt_ratio_nai(1000.0, crystal_diameter_mm=76.0)
    assert pt == pytest.approx(0.36, abs=0.02)


def test_F295_pt_ratio_for_detector_id_gamma_1c():
    from gamma.activity.pt_ratio_nai import pt_ratio_for_detector
    pt = pt_ratio_for_detector(662.0, "Gamma-1S")
    assert 0.1 < pt < 0.7


def test_F295_total_efficiency_inverse():
    from gamma.activity.pt_ratio_nai import (
        pt_ratio_nai, total_efficiency_from_fep,
    )
    pt = pt_ratio_nai(662.0, crystal_diameter_mm=63.0)
    eps_T = total_efficiency_from_fep(eps_fep=0.01, E_keV=662.0)
    assert eps_T == pytest.approx(0.01 / pt, rel=1e-6)


# ──────────────────────────────────────────────────────────────────
# F-296 — TCS correction close-geometry (T-007 + T-082)
# ──────────────────────────────────────────────────────────────────

def test_F296_tcs_co60_significant_close_geometry():
    """Co-60 каскад 1173+1332 в close geometry даёт ≥ 5 % correction."""
    from gamma.activity.tcs_close_geometry import (
        compute_tcs_correction_for_nuclide,
    )
    # Простая модель ε_T (close geometry — большая ε_T): возвращаем 0.10
    def eps_T(E): return 0.10
    res = compute_tcs_correction_for_nuclide(
        E_i_keV=1173.2,
        nuclide_id="Co-60",
        total_efficiency_func=eps_T,
    )
    assert res.is_significant
    assert res.correction_factor > 1.05


def test_F296_tcs_cs137_no_correction():
    """Cs-137 — одиночная линия, TCS correction = 1.0."""
    from gamma.activity.tcs_close_geometry import (
        compute_tcs_correction_for_nuclide,
    )
    def eps_T(E): return 0.10
    res = compute_tcs_correction_for_nuclide(
        E_i_keV=661.66,
        nuclide_id="Cs-137",
        total_efficiency_func=eps_T,
    )
    assert res.correction_factor == pytest.approx(1.0, abs=1e-6)
    assert not res.is_significant


def test_F296_tcs_geometry_heuristic():
    from gamma.activity.tcs_close_geometry import (
        is_tcs_significant_for_geometry,
    )
    assert is_tcs_significant_for_geometry(0.5)   # Marinelli surface
    assert is_tcs_significant_for_geometry(5.0)
    assert not is_tcs_significant_for_geometry(25.0)   # point source far


def test_F296_unknown_nuclide_raises():
    from gamma.activity.tcs_close_geometry import (
        compute_tcs_correction_for_nuclide,
    )
    with pytest.raises(KeyError):
        compute_tcs_correction_for_nuclide(
            E_i_keV=100.0, nuclide_id="Fictional-999",
            total_efficiency_func=lambda E: 0.05,
        )


# ──────────────────────────────────────────────────────────────────
# F-297 — Matrix method χ² (T-027)
# ──────────────────────────────────────────────────────────────────

def test_F297_matrix_method_single_nuclide_recovers_activity():
    """Один нуклид, один peak — A должно восстановиться точно."""
    from gamma.activity.matrix_method_chi2 import (
        PeakObservation, NuclideContribution, solve_matrix_method,
    )
    # Cs-137 1000 Bq → expected counts = ε·I·t·A = 0.015·0.851·1000·1000 = 12765
    A_true = 1000.0
    eps, I, t = 0.015, 0.851, 1000.0
    counts = eps * I * t * A_true
    peaks = [PeakObservation(E_keV=661.66, counts=counts)]
    contribs = {
        "Cs-137": [NuclideContribution(
            nuclide="Cs-137", E_keV=661.66, intensity_decimal=I,
            efficiency=eps, live_time_seconds=t,
        )]
    }
    res = solve_matrix_method(peaks, contribs)
    assert res.activities_Bq["Cs-137"] == pytest.approx(A_true, rel=0.02)


def test_F297_matrix_method_two_nuclides_two_peaks():
    """2 нуклида, 2 peak — система определима, решение единственно."""
    from gamma.activity.matrix_method_chi2 import (
        PeakObservation, NuclideContribution, solve_matrix_method,
    )
    A_cs, A_co = 500.0, 800.0
    t = 1000.0
    # Cs-137 в 662 → counts1
    counts1 = 0.015 * 0.851 * t * A_cs
    # Co-60 в 1173 → counts2
    counts2 = 0.010 * 0.9985 * t * A_co
    peaks = [
        PeakObservation(E_keV=661.66, counts=counts1),
        PeakObservation(E_keV=1173.2, counts=counts2),
    ]
    contribs = {
        "Cs-137": [NuclideContribution("Cs-137", 661.66, 0.851, 0.015, t)],
        "Co-60":  [NuclideContribution("Co-60",  1173.2, 0.9985, 0.010, t)],
    }
    res = solve_matrix_method(peaks, contribs)
    assert res.activities_Bq["Cs-137"] == pytest.approx(A_cs, rel=0.05)
    assert res.activities_Bq["Co-60"]  == pytest.approx(A_co, rel=0.05)
    assert res.is_acceptable


def test_F297_matrix_method_under_determined_raises():
    from gamma.activity.matrix_method_chi2 import (
        PeakObservation, NuclideContribution, solve_matrix_method,
    )
    peaks = [PeakObservation(E_keV=100.0, counts=500.0)]
    contribs = {
        "A": [NuclideContribution("A", 100.0, 1.0, 0.01, 1000.0)],
        "B": [NuclideContribution("B", 100.0, 1.0, 0.01, 1000.0)],
    }
    # 1 peak, 2 nuclides → under-determined
    with pytest.raises(ValueError):
        solve_matrix_method(peaks, contribs)


# ──────────────────────────────────────────────────────────────────
# F-298 — bg lines builder (T-013)
# ──────────────────────────────────────────────────────────────────

def test_F298_filter_k40_in_high_window():
    from gamma.activity.bg_lines_builder import filter_bg_lines_in_window
    lines = filter_bg_lines_in_window(1450.0, 1480.0)
    assert any(L.nuclide == "K-40" for L in lines)


def test_F298_filter_th232_chain_around_2614():
    from gamma.activity.bg_lines_builder import filter_bg_lines_in_window
    lines = filter_bg_lines_in_window(2600.0, 2640.0)
    assert any(L.nuclide == "Tl-208" and abs(L.E_keV - 2614.51) < 1.0
               for L in lines)


def test_F298_filter_by_parent_chain():
    from gamma.activity.bg_lines_builder import filter_bg_lines_in_window
    lines = filter_bg_lines_in_window(
        100.0, 3000.0, parent_chains=["Th-232"],
    )
    assert lines
    assert all(L.parent_chain == "Th-232" for L in lines)


def test_F298_build_f131_input_sorted_by_E():
    from gamma.activity.bg_lines_builder import build_f131_input
    inp = build_f131_input(roi_E_min_keV=100.0, roi_E_max_keV=2700.0)
    assert inp.expected_peaks_keV == sorted(inp.expected_peaks_keV)
    assert len(inp.expected_peak_labels) == len(inp.expected_peaks_keV)


def test_F298_anchor_candidates_include_K40_and_Tl208():
    from gamma.activity.bg_lines_builder import get_anchor_candidates
    anchors = get_anchor_candidates()
    nuclides = {a.nuclide for a in anchors}
    assert "K-40" in nuclides
    assert "Tl-208" in nuclides


def test_F298_chain_dominance_th_chain():
    """Если detected lines содержат 583+2614 → доминирует Th-232 chain."""
    from gamma.activity.bg_lines_builder import classify_chain_dominance
    counts = classify_chain_dominance([583.0, 2614.0, 911.0])
    assert counts.get("Th-232", 0) >= 3
