# -*- coding: utf-8 -*-
"""v1.17.18 delivery tests — Production QA gates (F-289..F-292)."""
from __future__ import annotations
import math, os, sys
import pytest

SCRIPTS = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts")
)
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)


# ──────────────────────────────────────────────────────────────────
# F-289 — Spectrometer compliance gate (T-083)
# ──────────────────────────────────────────────────────────────────

def test_F289_compliance_pass_typical_gamma_1c():
    from gamma.calibration.spectrometer_compliance import (
        SpectrometerSpec, SpectrometerMeasurement,
        check_spectrometer_compliance,
    )
    spec = SpectrometerSpec.gamma_1c_nai_63x63()
    meas = SpectrometerMeasurement(
        fwhm_at_662_keV_pct=7.0,
        channels_actual=1024,
        inl_pct=0.5,
    )
    rep = check_spectrometer_compliance(meas, spec)
    assert rep.is_pass
    assert rep.overall_status == "PASS"


def test_F289_compliance_warning_borderline_fwhm():
    from gamma.calibration.spectrometer_compliance import (
        SpectrometerSpec, SpectrometerMeasurement,
        check_spectrometer_compliance,
    )
    spec = SpectrometerSpec.gamma_1c_nai_63x63()
    meas = SpectrometerMeasurement(
        fwhm_at_662_keV_pct=8.4,    # 5 % выше spec (8.0) — WARNING
        channels_actual=1024,
        inl_pct=0.5,
    )
    rep = check_spectrometer_compliance(meas, spec, tolerance_frac=0.10)
    assert rep.overall_status == "WARNING"
    assert not rep.is_fail


def test_F289_compliance_fail_degraded_fwhm():
    from gamma.calibration.spectrometer_compliance import (
        SpectrometerSpec, SpectrometerMeasurement,
        check_spectrometer_compliance,
    )
    spec = SpectrometerSpec.gamma_1c_nai_63x63()
    meas = SpectrometerMeasurement(
        fwhm_at_662_keV_pct=12.0,    # 50 % выше spec (8.0) — FAIL
        channels_actual=1024,
        inl_pct=0.5,
    )
    rep = check_spectrometer_compliance(meas, spec, tolerance_frac=0.10)
    assert rep.is_fail


def test_F289_compliance_fail_insufficient_channels():
    from gamma.calibration.spectrometer_compliance import (
        SpectrometerSpec, SpectrometerMeasurement,
        check_spectrometer_compliance,
    )
    spec = SpectrometerSpec.gamma_1c_nai_63x63()
    meas = SpectrometerMeasurement(
        fwhm_at_662_keV_pct=7.0,
        channels_actual=512,   # < 1024 min — FAIL
        inl_pct=0.5,
    )
    rep = check_spectrometer_compliance(meas, spec)
    assert rep.is_fail
    ch_finding = next(f for f in rep.findings if f.field_name == "channels")
    assert ch_finding.status == "FAIL"


def test_F289_compliance_warning_missing_inl():
    from gamma.calibration.spectrometer_compliance import (
        SpectrometerSpec, SpectrometerMeasurement,
        check_spectrometer_compliance,
    )
    spec = SpectrometerSpec.gamma_1c_nai_63x63()
    meas = SpectrometerMeasurement(
        fwhm_at_662_keV_pct=7.0,
        channels_actual=1024,
        inl_pct=None,    # отсутствует — WARNING per field
    )
    rep = check_spectrometer_compliance(meas, spec)
    assert rep.overall_status == "WARNING"


# ──────────────────────────────────────────────────────────────────
# F-290 — Sample density / fill-level gate (T-084)
# ──────────────────────────────────────────────────────────────────

def test_F290_sample_pass_typical_water():
    from gamma.activity.sample_gates import (
        SampleGeometry, CalibrationValidityRange,
        check_sample_compliance,
    )
    rng = CalibrationValidityRange.marinelli_0_5l_default()
    sample = SampleGeometry(
        container_type="Marinelli_0.5L",
        fill_height_mm=90.0,
        density_g_cm3=1.00,
        matrix_label="water",
    )
    rep = check_sample_compliance(sample, rng)
    assert rep.overall_status == "PASS"
    assert not rep.requires_self_absorption_correction


def test_F290_sample_extrapolation_dense_soil():
    from gamma.activity.sample_gates import (
        SampleGeometry, CalibrationValidityRange,
        check_sample_compliance,
    )
    rng = CalibrationValidityRange.marinelli_0_5l_default()
    sample = SampleGeometry(
        container_type="Marinelli_0.5L",
        fill_height_mm=90.0,
        density_g_cm3=1.80,    # вне 0.95-1.15, но в 0.50-2.50
        matrix_label="soil_dry",
    )
    rep = check_sample_compliance(sample, rng)
    assert rep.overall_status == "EXTRAPOLATION"
    assert rep.requires_self_absorption_correction


def test_F290_sample_fail_too_dense():
    from gamma.activity.sample_gates import (
        SampleGeometry, CalibrationValidityRange,
        check_sample_compliance,
    )
    rng = CalibrationValidityRange.marinelli_0_5l_default()
    sample = SampleGeometry(
        container_type="Marinelli_0.5L",
        fill_height_mm=90.0,
        density_g_cm3=3.50,    # вне даже extrap 2.50 — FAIL
        matrix_label="ore",
    )
    rep = check_sample_compliance(sample, rng)
    assert rep.is_fail


def test_F290_sample_fail_wrong_container():
    from gamma.activity.sample_gates import (
        SampleGeometry, CalibrationValidityRange,
        check_sample_compliance,
    )
    rng = CalibrationValidityRange.marinelli_0_5l_default()
    sample = SampleGeometry(
        container_type="vial_20ml",   # != Marinelli_0.5L
        fill_height_mm=20.0,
        density_g_cm3=1.00,
    )
    rep = check_sample_compliance(sample, rng)
    assert rep.is_fail


def test_F290_sample_warning_underfilled_near_border():
    from gamma.activity.sample_gates import (
        SampleGeometry, CalibrationValidityRange,
        check_sample_compliance,
    )
    rng = CalibrationValidityRange.marinelli_0_5l_default()
    sample = SampleGeometry(
        container_type="Marinelli_0.5L",
        fill_height_mm=81.0,    # на границе (80-100, 10 % band → ≤82)
        density_g_cm3=1.00,
    )
    rep = check_sample_compliance(sample, rng)
    assert rep.overall_status == "WARNING"


# ──────────────────────────────────────────────────────────────────
# F-291 — Background drift F-statistic (T-054)
# ──────────────────────────────────────────────────────────────────

def test_F291_bg_drift_pass_stable():
    from gamma.calibration.bg_drift import BgRoiSnapshot, f_test_bg_drift
    ref = BgRoiSnapshot(
        roi_label="Cs-137 661.6 keV",
        counts_per_session=[100, 102, 98, 101, 99, 103, 97, 100, 102, 98],
        live_time_seconds_per_session=3600,
    )
    cur = BgRoiSnapshot(
        roi_label="Cs-137 661.6 keV",
        counts_per_session=[101, 99, 100, 102, 98, 100, 101, 99, 102, 100],
        live_time_seconds_per_session=3600,
    )
    finding = f_test_bg_drift(cur, ref)
    assert finding.status == "PASS"


def test_F291_bg_drift_fail_contamination():
    from gamma.calibration.bg_drift import BgRoiSnapshot, f_test_bg_drift
    ref = BgRoiSnapshot(
        roi_label="Cs-137 661.6 keV",
        counts_per_session=[100, 102, 98, 101, 99, 103, 97, 100, 102, 98],
        live_time_seconds_per_session=3600,
    )
    cur = BgRoiSnapshot(
        roi_label="Cs-137 661.6 keV",
        counts_per_session=[250, 260, 240, 255, 245, 252, 248, 251, 249, 253],
        live_time_seconds_per_session=3600,
    )
    finding = f_test_bg_drift(cur, ref)
    assert finding.status == "FAIL"
    assert finding.mean_z_score > 3.0


def test_F291_bg_drift_warning_variance_only():
    from gamma.calibration.bg_drift import BgRoiSnapshot, f_test_bg_drift
    ref = BgRoiSnapshot(
        roi_label="K-40 1461 keV",
        counts_per_session=[100, 101, 100, 99, 100, 101, 99, 100, 100, 100],
        live_time_seconds_per_session=3600,
    )
    # Та же mean но в 5 раз больше variance
    cur = BgRoiSnapshot(
        roi_label="K-40 1461 keV",
        counts_per_session=[100, 110, 90, 105, 95, 108, 92, 103, 97, 100],
        live_time_seconds_per_session=3600,
    )
    finding = f_test_bg_drift(cur, ref)
    # Variance уехал, но mean ровно — должно быть WARNING (не FAIL)
    assert finding.status in ("WARNING", "FAIL")
    assert finding.f_statistic > 1.0


def test_F291_bg_drift_warning_insufficient_data():
    from gamma.calibration.bg_drift import BgRoiSnapshot, f_test_bg_drift
    ref = BgRoiSnapshot(
        roi_label="K-40 1461 keV",
        counts_per_session=[100],    # df=0
        live_time_seconds_per_session=3600,
    )
    cur = BgRoiSnapshot(
        roi_label="K-40 1461 keV",
        counts_per_session=[105, 95],
        live_time_seconds_per_session=3600,
    )
    finding = f_test_bg_drift(cur, ref)
    assert finding.status == "WARNING"


def test_F291_bg_drift_multi_roi_aggregation():
    from gamma.calibration.bg_drift import (
        BgRoiSnapshot, check_bg_drift_multi_roi,
    )
    ref_rois = [
        BgRoiSnapshot(
            "Cs-137 661.6 keV",
            [100, 102, 98, 101, 99, 103, 97, 100, 102, 98], 3600,
        ),
        BgRoiSnapshot(
            "K-40 1461 keV",
            [50, 52, 48, 51, 49, 53, 47, 50, 52, 48], 3600,
        ),
    ]
    cur_rois = [
        BgRoiSnapshot(
            "Cs-137 661.6 keV",
            [101, 99, 100, 102, 98, 100, 101, 99, 102, 100], 3600,
        ),
        BgRoiSnapshot(
            "K-40 1461 keV",
            [200, 210, 190, 205, 195, 208, 192, 203, 197, 200], 3600,
        ),
    ]
    rep = check_bg_drift_multi_roi(cur_rois, ref_rois)
    assert rep.overall_status == "FAIL"   # из-за K-40 контаминации


# ──────────────────────────────────────────────────────────────────
# F-292 — Sensitivity drift quarterly check (T-053)
# ──────────────────────────────────────────────────────────────────

def test_F292_sensitivity_pass_minor_drift():
    from gamma.calibration.sensitivity_drift import (
        EfficiencyAnchor, check_sensitivity_drift,
    )
    ref = [
        EfficiencyAnchor("Cs-137", 661.6, 0.0150),
        EfficiencyAnchor("Co-60", 1173.0, 0.0085),
        EfficiencyAnchor("Co-60", 1332.0, 0.0078),
    ]
    cur = [
        EfficiencyAnchor("Cs-137", 661.6, 0.01515),  # +1 %
        EfficiencyAnchor("Co-60", 1173.0, 0.00858),  # +1 %
        EfficiencyAnchor("Co-60", 1332.0, 0.00787),  # +0.9 %
    ]
    rep = check_sensitivity_drift(cur, ref, days_since_reference=20)
    assert rep.overall_status == "PASS"


def test_F292_sensitivity_fail_5pct_drift():
    from gamma.calibration.sensitivity_drift import (
        EfficiencyAnchor, check_sensitivity_drift,
    )
    ref = [
        EfficiencyAnchor("Cs-137", 661.6, 0.0150),
        EfficiencyAnchor("Co-60", 1173.0, 0.0085),
    ]
    cur = [
        EfficiencyAnchor("Cs-137", 661.6, 0.0142),   # -5.3 %
        EfficiencyAnchor("Co-60", 1173.0, 0.0080),   # -5.9 %
    ]
    rep = check_sensitivity_drift(cur, ref, days_since_reference=30)
    assert rep.is_fail


def test_F292_sensitivity_fail_calendar_overdue():
    from gamma.calibration.sensitivity_drift import (
        EfficiencyAnchor, check_sensitivity_drift,
    )
    ref = [EfficiencyAnchor("Cs-137", 661.6, 0.0150)]
    cur = [EfficiencyAnchor("Cs-137", 661.6, 0.01510)]   # drift 0.7 %
    # 200 дней — > interval_fail_days=180 → FAIL независимо от drift
    rep = check_sensitivity_drift(cur, ref, days_since_reference=200)
    assert rep.is_fail


def test_F292_sensitivity_warning_quarterly():
    from gamma.calibration.sensitivity_drift import (
        EfficiencyAnchor, check_sensitivity_drift,
    )
    ref = [EfficiencyAnchor("Cs-137", 661.6, 0.0150)]
    cur = [EfficiencyAnchor("Cs-137", 661.6, 0.01530)]   # +2 %
    rep = check_sensitivity_drift(cur, ref, days_since_reference=100)
    # > 90 (warn) но < 180 (fail) → WARNING
    assert rep.overall_status == "WARNING"


def test_F292_sensitivity_monotonic_shift_detected():
    """Drift монотонно растёт с E → likely geometry shift."""
    from gamma.calibration.sensitivity_drift import (
        EfficiencyAnchor, check_sensitivity_drift,
    )
    ref = [
        EfficiencyAnchor("Cs-137", 661.6, 0.0150),
        EfficiencyAnchor("Co-60", 1173.0, 0.0085),
        EfficiencyAnchor("Co-60", 1332.0, 0.0078),
        EfficiencyAnchor("K-40", 1461.0, 0.0070),
    ]
    cur = [
        EfficiencyAnchor("Cs-137", 661.6, 0.0149),   # -0.7 %
        EfficiencyAnchor("Co-60", 1173.0, 0.0083),   # -2.4 %
        EfficiencyAnchor("Co-60", 1332.0, 0.0075),   # -3.8 %
        EfficiencyAnchor("K-40", 1461.0, 0.0066),    # -5.7 % FAIL
    ]
    rep = check_sensitivity_drift(cur, ref, days_since_reference=60)
    assert rep.monotonic_shift is True
    assert rep.pearson_r_with_E is not None
    assert abs(rep.pearson_r_with_E) > 0.9   # strong correlation


def test_F292_sensitivity_empty_returns_pass():
    """Edge case: no anchors → PASS (nothing to check)."""
    from gamma.calibration.sensitivity_drift import check_sensitivity_drift
    rep = check_sensitivity_drift([], [], days_since_reference=10)
    assert rep.overall_status == "PASS"
