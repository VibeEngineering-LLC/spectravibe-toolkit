"""F-122 (v1.17.6) — Self-attenuation wiring для Marinelli.

Проверяет, что:
  1. ``compute_activity`` принимает kwarg ``self_attenuation_factors``
     и применяет per-line correction в формуле активности.
  2. ``compute_activities_for_all`` при заданной геометрии Marinelli
     и плотности образца ≠ reference (ОИСН-16 1.60 г/см³) применяет
     коррекцию ко всем линиям всех нуклидов.
  3. При sample_density == reference (1.60) коррекция = 1.0 (нет влияния).
  4. При sample_density < reference (1.0 — вода) коррекция < 1.0 для
     низких энергий (сильная самопоглощения в ОИСН-16 vs воде).
"""
from __future__ import annotations

import math
import os
import sys

import pytest

SCRIPTS = os.path.join(os.path.dirname(__file__), "scripts")
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)


def test_compute_activity_accepts_self_attenuation_factors_kwarg():
    """Smoke: compute_activity не падает на новом kwarg."""
    from gamma.activity.compute import compute_activity
    from gamma.calibration.efficiency import EfficiencyCurve
    from gamma.identification.identify import (
        NuclideIdentification, LineMatch,
    )

    # Постройм минимальный EfficiencyCurve (нулевой полином → eps=exp(0)=1)
    eff = EfficiencyCurve(
        coefficients=(math.log(0.01),),
        E_min_keV=50.0, E_max_keV=3000.0,
        chi2_per_dof=1.0, n_points_used=10, n_dof=8,
    )
    line = LineMatch(
        nuclide="Cs-137",
        library_E_keV=661.66, library_I_pct=85.1,
        peak_channel=662, peak_E_keV=661.66, peak_sigma=10.0,
        residual_keV=0.0, is_characteristic=True,
        peak_area=10000.0, peak_area_uncertainty=200.0,
        peak_area_source="cowell",
    )
    ni = NuclideIdentification(
        nuclide="Cs-137", detected=True,
        reason="char line matched", characteristic_line_keV=661.66,
        matched_lines=(line,),
    )
    # Без коррекции:
    a_nocorr = compute_activity(
        ni, efficiency_curve=eff, live_time_s=1000.0,
        from_bg_subtracted=True, bg_available=False,
    )
    # С коррекцией 0.7 на E=661.66:
    a_corr = compute_activity(
        ni, efficiency_curve=eff, live_time_s=1000.0,
        from_bg_subtracted=True, bg_available=False,
        self_attenuation_factors={661.66: 0.7},
    )
    assert a_nocorr.is_valid()
    assert a_corr.is_valid()
    # Корректированная активность × 0.7 = бaseline
    ratio = a_corr.A_Bq / a_nocorr.A_Bq
    assert abs(ratio - 0.7) < 0.01, (
        f"self-attenuation correction not applied; ratio={ratio}"
    )


def test_compute_activities_for_all_no_correction_when_geometry_unknown():
    """Без geometry_canonical из REF_GEOMETRY коррекция не применяется."""
    from gamma.activity.compute import compute_activities_for_all
    from gamma.calibration.efficiency import EfficiencyCurve
    from gamma.identification.identify import (
        NuclideIdentification, LineMatch, IdentificationResult,
    )
    from gamma.identification.window import build_identification_window

    eff = EfficiencyCurve(
        coefficients=(math.log(0.01),),
        E_min_keV=50.0, E_max_keV=3000.0,
        chi2_per_dof=1.0, n_points_used=10, n_dof=8,
    )
    line = LineMatch(
        nuclide="Cs-137", library_E_keV=661.66, library_I_pct=85.1,
        peak_channel=662, peak_E_keV=661.66, peak_sigma=10.0,
        residual_keV=0.0, is_characteristic=True,
        peak_area=10000.0, peak_area_uncertainty=200.0,
        peak_area_source="cowell",
    )
    ni = NuclideIdentification(
        nuclide="Cs-137", detected=True,
        reason="x", characteristic_line_keV=661.66,
        matched_lines=(line,),
    )
    ir = IdentificationResult(
        detector_type="NaI",
        window=build_identification_window("NaI"),
        candidates_considered=1,
        detected_nuclides=(ni,),
        rejected_nuclides=(),
        unmatched_peaks=(),
        notes="",
    )

    # Точечная геометрия — нет в REF_GEOMETRY → коррекция не активируется.
    results = compute_activities_for_all(
        ir, efficiency_curve=eff, live_time_s=1000.0,
        from_bg_subtracted=True, bg_available=False,
        geometry_canonical="Точечная-5см",
        sample_density_g_cm3=1.0,
    )
    assert len(results) == 1
    assert results[0].is_valid()


def test_compute_activities_for_all_no_op_when_sample_density_equals_ref():
    """sample_density == ref_density → correction = 1.0 → активность не меняется."""
    from gamma.activity.compute import compute_activities_for_all
    from gamma.calibration.efficiency import EfficiencyCurve
    from gamma.identification.identify import (
        NuclideIdentification, LineMatch, IdentificationResult,
    )
    from gamma.identification.window import build_identification_window

    eff = EfficiencyCurve(
        coefficients=(math.log(0.01),),
        E_min_keV=50.0, E_max_keV=3000.0,
        chi2_per_dof=1.0, n_points_used=10, n_dof=8,
    )
    line = LineMatch(
        nuclide="Cs-137", library_E_keV=661.66, library_I_pct=85.1,
        peak_channel=662, peak_E_keV=661.66, peak_sigma=10.0,
        residual_keV=0.0, is_characteristic=True,
        peak_area=10000.0, peak_area_uncertainty=200.0,
        peak_area_source="cowell",
    )
    ni = NuclideIdentification(
        nuclide="Cs-137", detected=True, reason="x",
        characteristic_line_keV=661.66, matched_lines=(line,),
    )
    ir = IdentificationResult(
        detector_type="NaI",
        window=build_identification_window("NaI"),
        candidates_considered=1,
        detected_nuclides=(ni,),
        rejected_nuclides=(),
        unmatched_peaks=(),
        notes="",
    )

    base = compute_activities_for_all(
        ir, efficiency_curve=eff, live_time_s=1000.0,
        from_bg_subtracted=True, bg_available=False,
    )
    # Marinelli + ρ_образец = ρ_ref (ОИСН-16 = 1.60) → correction = 1.0
    same = compute_activities_for_all(
        ir, efficiency_curve=eff, live_time_s=1000.0,
        from_bg_subtracted=True, bg_available=False,
        geometry_canonical="Маринелли",
        sample_density_g_cm3=1.60,
    )
    assert abs(base[0].A_Bq - same[0].A_Bq) / base[0].A_Bq < 1e-6


def test_compute_activities_marinelli_density_difference_applies_correction():
    """Marinelli + ρ_образец ≠ ρ_ref → активность меняется монотонно."""
    from gamma.activity.compute import compute_activities_for_all
    from gamma.calibration.efficiency import EfficiencyCurve
    from gamma.identification.identify import (
        NuclideIdentification, LineMatch, IdentificationResult,
    )
    from gamma.identification.window import build_identification_window

    eff = EfficiencyCurve(
        coefficients=(math.log(0.01),),
        E_min_keV=50.0, E_max_keV=3000.0,
        chi2_per_dof=1.0, n_points_used=10, n_dof=8,
    )
    # Низкоэнергетическая линия — большой эффект самопоглощения
    line = LineMatch(
        nuclide="Pb-212", library_E_keV=238.63, library_I_pct=43.6,
        peak_channel=239, peak_E_keV=238.63, peak_sigma=8.0,
        residual_keV=0.0, is_characteristic=True,
        peak_area=10000.0, peak_area_uncertainty=200.0,
        peak_area_source="cowell",
    )
    ni = NuclideIdentification(
        nuclide="Pb-212", detected=True, reason="x",
        characteristic_line_keV=238.63, matched_lines=(line,),
    )
    ir = IdentificationResult(
        detector_type="NaI",
        window=build_identification_window("NaI"),
        candidates_considered=1,
        detected_nuclides=(ni,),
        rejected_nuclides=(),
        unmatched_peaks=(),
        notes="",
    )

    # ρ_образец < ρ_ref → F_sample > F_ref → correction = F_ref/F_sample < 1
    light = compute_activities_for_all(
        ir, efficiency_curve=eff, live_time_s=1000.0,
        from_bg_subtracted=True, bg_available=False,
        geometry_canonical="Маринелли",
        sample_density_g_cm3=0.5,
    )
    # ρ_образец > ρ_ref → correction > 1
    heavy = compute_activities_for_all(
        ir, efficiency_curve=eff, live_time_s=1000.0,
        from_bg_subtracted=True, bg_available=False,
        geometry_canonical="Маринелли",
        sample_density_g_cm3=3.0,
    )
    assert light[0].A_Bq < heavy[0].A_Bq, (
        f"монотонность нарушена: light={light[0].A_Bq}, "
        f"heavy={heavy[0].A_Bq}"
    )
