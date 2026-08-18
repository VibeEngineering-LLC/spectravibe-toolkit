"""F-446: Phase C adaptive polynomial degree tests.

n<3   -> deg=0 (constant shift, preserves spacing)
n==3  -> deg=1 (linear, 1 redundant point)
n>=4  -> deg=min(max_degree, n-1) (parabola for NaI cap)
"""
from __future__ import annotations
from types import SimpleNamespace
from typing import List

import numpy as np
import pytest

from gamma.calibration.multiplet_self_calibration import (
    CentroidAnchor,
    SelfCalibrationDiag,
    PHASE_C_ADAPTIVE_DEGREE,
    PHASE_C_MAX_DEGREE_NAI,
    PHASE_C_MIN_ANCHORS_FOR_LINEAR,
    PHASE_C_MIN_ANCHORS_FOR_PARABOLA,
    recalibrate_from_multiplet_centroids,
)
from gamma.calibration import multiplet_self_calibration as mscal


def _make_spec(stored_cal):
    a0, a1 = stored_cal[0], stored_cal[1]
    a2 = stored_cal[2] if len(stored_cal) > 2 else 0.0

    def e2c(E):
        return (E - a0) / a1

    spec = SimpleNamespace(
        energy_cal=tuple(stored_cal),
        energy_to_channel=e2c,
        channel_to_energy=lambda N: a0 + a1 * N + a2 * N * N,
    )
    return spec


def _fwhm_provider(E):
    return 0.05 * (E ** 0.5) * 10.0


def _build_anchors(E_passports, delta_targets, stored_cal, fwhm_provider):
    a0, a1 = stored_cal[0], stored_cal[1]
    anchors = []
    for idx in range(len(E_passports)):
        E = E_passports[idx]
        d = delta_targets[idx]
        ch = (E - d - a0) / a1
        fwhm = fwhm_provider(E)
        anchors.append(CentroidAnchor(
            nuclide="Nuc_" + str(idx),
            E_passport_keV=float(E),
            E_fitted_keV=float(E),
            channel_fitted=float(ch),
            fwhm_keV=float(fwhm),
            drift_fraction_of_fwhm=0.0,
            source="cluster_delta_" + str(idx),
        ))
    return anchors


def _eval_cal(N, cal):
    E = 0.0
    for c in reversed(cal):
        E = E * N + c
    return E


def test_phase_c_degree_n_lt_3_is_constant():
    deg1, reason1 = mscal._f446_choose_phase_c_degree(1, PHASE_C_MAX_DEGREE_NAI)
    deg2, reason2 = mscal._f446_choose_phase_c_degree(2, PHASE_C_MAX_DEGREE_NAI)
    assert deg1 == 0
    assert deg2 == 0
    assert "constant" in reason1
    assert "constant" in reason2


def test_phase_c_degree_n_eq_3_is_linear():
    deg, reason = mscal._f446_choose_phase_c_degree(3, PHASE_C_MAX_DEGREE_NAI)
    assert deg == 1
    assert "linear" in reason


def test_phase_c_degree_n_ge_4_is_capped_parabola():
    for n in (4, 5, 10):
        deg, _ = mscal._f446_choose_phase_c_degree(n, PHASE_C_MAX_DEGREE_NAI)
        assert deg == 2


def test_phase_c_degree_respects_lower_max_degree_cap():
    deg, _ = mscal._f446_choose_phase_c_degree(5, max_degree=1)
    assert deg == 1


def test_phase_c_adaptive_marker_constant_set():
    assert PHASE_C_ADAPTIVE_DEGREE is True


def test_phase_c_compute_constant_delta_uniform_mean():
    class _A:
        pass

    a1, a2, a3 = _A(), _A(), _A()
    out = mscal._f446_compute_constant_delta([a1, a2, a3], [1.0, 2.0, 3.0])
    assert abs(out - 2.0) < 1e-9


def test_phase_c_n_eq_2_uses_constant_shift_in_integration():
    stored_cal = (0.5, 3.0)
    spec = _make_spec(stored_cal)
    anchors = _build_anchors([941.0, 1591.0], [3.0, 3.5], stored_cal, _fwhm_provider)
    extras = []
    for a in anchors:
        extras.append((a.E_passport_keV, a.channel_fitted, a.source))
    new_cal, diag = recalibrate_from_multiplet_centroids(
        spec,
        fitted_multiplets=[],
        fwhm_provider_keV=_fwhm_provider,
        extra_anchors=extras,
        max_degree=PHASE_C_MAX_DEGREE_NAI,
        min_anchors=2,
        use_cluster_global=False,
    )
    assert new_cal is not None, "Phase C rejected: " + diag.reason
    assert diag.delta_degree_used == 0
    assert diag.delta_const_keV is not None
    assert abs(diag.delta_const_keV - 3.25) < 0.01


def test_phase_c_constant_preserves_spacing():
    stored_cal = (0.5, 3.0)
    spec = _make_spec(stored_cal)
    anchors = _build_anchors([600.0, 1500.0], [2.5, 2.5], stored_cal, _fwhm_provider)
    extras = []
    for a in anchors:
        extras.append((a.E_passport_keV, a.channel_fitted, a.source))
    new_cal, diag = recalibrate_from_multiplet_centroids(
        spec,
        fitted_multiplets=[],
        fwhm_provider_keV=_fwhm_provider,
        extra_anchors=extras,
        max_degree=PHASE_C_MAX_DEGREE_NAI,
        min_anchors=2,
        use_cluster_global=False,
    )
    assert new_cal is not None
    ch1 = anchors[0].channel_fitted
    ch2 = anchors[1].channel_fitted
    old_spacing = _eval_cal(ch2, stored_cal) - _eval_cal(ch1, stored_cal)
    new_spacing = _eval_cal(ch2, new_cal) - _eval_cal(ch1, new_cal)
    assert abs(new_spacing - old_spacing) < 1e-9


def test_phase_c_constant_no_extrapolation_blowup():
    stored_cal = (0.5, 3.0)
    spec = _make_spec(stored_cal)
    anchors = _build_anchors([900.5, 1800.5], [3.0, 3.6], stored_cal, _fwhm_provider)
    extras = []
    for a in anchors:
        extras.append((a.E_passport_keV, a.channel_fitted, a.source))
    new_cal, diag = recalibrate_from_multiplet_centroids(
        spec,
        fitted_multiplets=[],
        fwhm_provider_keV=_fwhm_provider,
        extra_anchors=extras,
        max_degree=PHASE_C_MAX_DEGREE_NAI,
        min_anchors=2,
        use_cluster_global=False,
    )
    assert new_cal is not None
    assert diag.delta_degree_used == 0
    max_anchor_abs = max(abs(3.0), abs(3.6))
    for ch_far in (80.0, 900.0):
        delta_applied = _eval_cal(ch_far, new_cal) - _eval_cal(ch_far, stored_cal)
        assert abs(delta_applied) <= max_anchor_abs * 1.1


def test_phase_c_n_eq_3_uses_linear_in_integration():
    stored_cal = (0.0, 3.0)
    spec = _make_spec(stored_cal)
    anchors = _build_anchors([300.0, 600.0, 900.0], [1.0, 2.0, 3.0], stored_cal, _fwhm_provider)
    extras = []
    for a in anchors:
        extras.append((a.E_passport_keV, a.channel_fitted, a.source))
    new_cal, diag = recalibrate_from_multiplet_centroids(
        spec,
        fitted_multiplets=[],
        fwhm_provider_keV=_fwhm_provider,
        extra_anchors=extras,
        max_degree=PHASE_C_MAX_DEGREE_NAI,
        min_anchors=3,
        use_cluster_global=False,
    )
    assert new_cal is not None
    assert diag.delta_degree_used == 1
    assert diag.delta_const_keV is None


def test_phase_c_n_eq_4_uses_parabola_capped_at_2_in_integration():
    stored_cal = (0.0, 3.0)
    spec = _make_spec(stored_cal)
    anchors = _build_anchors(
        [300.0, 600.0, 900.0, 1500.0],
        [1.0, 2.0, 3.0, 5.0],
        stored_cal,
        _fwhm_provider,
    )
    extras = []
    for a in anchors:
        extras.append((a.E_passport_keV, a.channel_fitted, a.source))
    new_cal, diag = recalibrate_from_multiplet_centroids(
        spec,
        fitted_multiplets=[],
        fwhm_provider_keV=_fwhm_provider,
        extra_anchors=extras,
        max_degree=PHASE_C_MAX_DEGREE_NAI,
        min_anchors=4,
        use_cluster_global=False,
    )
    assert new_cal is not None
    assert diag.delta_degree_used == 2
    assert diag.delta_const_keV is None


def test_phase_c_diag_choice_reason_recorded():
    stored_cal = (0.0, 3.0)
    spec = _make_spec(stored_cal)
    anchors = _build_anchors([300.0, 1500.0], [2.5, 2.5], stored_cal, _fwhm_provider)
    extras = []
    for a in anchors:
        extras.append((a.E_passport_keV, a.channel_fitted, a.source))
    _, diag = recalibrate_from_multiplet_centroids(
        spec,
        fitted_multiplets=[],
        fwhm_provider_keV=_fwhm_provider,
        extra_anchors=extras,
        max_degree=PHASE_C_MAX_DEGREE_NAI,
        min_anchors=2,
        use_cluster_global=False,
    )
    assert diag.degree_choice_reason
    assert "n_anchors=2" in diag.degree_choice_reason
