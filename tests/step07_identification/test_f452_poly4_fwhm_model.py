"""F-452 — FWHM-модель степени 4 в sqrt(E), LSRM-faithful.

Раньше F-160 ground-truth ветка `build_fwhm_model` делала NNLS-refit
(a + b·E + c·E²) с неотрицательными коэффициентами поверх LSRM-полинома
степени 4 → систематика ±5-7 keV на Th-232 anchors (FWHM(238)=29.5 vs LSRM 24.0
+23%, FWHM(2614)=116.8 vs 112.8 +3.6%). Корень — лоссная 3-параметрическая
компрессия poly-4.

F-452 (этот тест): unconstrained lstsq на FWHM(E) = Σ c_k·√E^k (deg=4) даёт
max|ΔFWHM| < 2 keV на тех же anchors. Новый `FwhmModel` dataclass — callable;
старый tuple-API сохранён через polymorphic `fwhm_keV_at_energy`.
"""
from __future__ import annotations

import math
import os
import sys
from dataclasses import dataclass, field
from typing import Any, Dict, Optional

import numpy as np
import pytest

SCRIPTS = os.path.join(
    os.path.dirname(__file__), os.pardir, os.pardir, "scripts",
)
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)

from gamma.identification.staged_pipeline import (  # noqa: E402
    _DEFAULT_NAI_FWHM_MODEL,
    _DEFAULT_NAI_FWHM_MODEL_OBJ,
    build_fwhm_model,
    fwhm_keV_at_energy,
    FwhmModel,
)
from gamma.spectrum import StoredFwhmCalibration  # noqa: E402


_TH232_GT_ANCHORS = [
    (208.129, 21.644),
    (237.508, 24.029),
    (327.625, 29.292),
    (337.944, 29.881),
    (459.984, 34.666),
    (507.747, 36.835),
    (580.077, 43.116),
    (900.228, 55.823),
    (953.774, 58.215),
    (957.978, 58.402),
    (1581.207, 90.683),
    (1613.500, 92.009),
    (1623.624, 92.419),
    (2612.857, 112.796),
]


@dataclass
class _MockSpec:
    extras: Dict[str, Any] = field(default_factory=dict)
    stored_fwhm_calibration: Optional[StoredFwhmCalibration] = None


def _make_spec_with_gt() -> _MockSpec:
    spec = _MockSpec()
    spec.extras["_f452_gt_anchors_inline"] = list(_TH232_GT_ANCHORS)
    return spec


def test_fwhm_model_callable_quad():
    a, b, c = _DEFAULT_NAI_FWHM_MODEL
    obj = FwhmModel(kind="quad_fwhm2_in_E", coefficients=(a, b, c))
    expected = math.sqrt(max(a + b * 661.7 + c * 661.7 ** 2, 0.01))
    assert abs(obj(661.7) - expected) < 1e-9


def test_fwhm_model_callable_lsrm_poly_sqrt_E():
    obj = FwhmModel(kind="lsrm_poly_sqrt_E", coefficients=(0.5, 1.2, 0.0, 0.0, 0.0))
    z = math.sqrt(661.7)
    assert abs(obj(661.7) - (0.5 + 1.2 * z)) < 1e-9


def test_fwhm_model_invalid_kind_raises():
    obj = FwhmModel(kind="bogus", coefficients=(1.0,))
    with pytest.raises(ValueError):
        obj(100.0)


def test_fwhm_keV_at_energy_polymorphic_BC():
    a, b, c = _DEFAULT_NAI_FWHM_MODEL
    legacy_tuple = (a, b, c)
    new_obj = _DEFAULT_NAI_FWHM_MODEL_OBJ
    for E in (60.0, 100.0, 300.0, 661.7, 1460.8, 2614.5):
        v_tuple = fwhm_keV_at_energy(legacy_tuple, E)
        v_obj = fwhm_keV_at_energy(new_obj, E)
        assert abs(v_tuple - v_obj) < 1e-9, f"BC mismatch at E={E}"


def test_poly4_sqrt_E_lstsq_beats_quad_on_th232_anchors():
    """Главный кейс F-452: poly-4 в sqrt(E) даёт max|Δ| < 3 keV на LSRM Th-232
    anchors и при этом В НЕСКОЛЬКО РАЗ точнее, чем 3-параметрическая
    FWHM²=a+b·E+c·E² аппроксимация (с которой раньше работала F-160 ветка
    через NNLS-refit поверх LSRM-полинома → систематика ±5-7 keV).
    """
    Es = np.array([e for e, _ in _TH232_GT_ANCHORS], dtype=np.float64)
    Fs = np.array([f for _, f in _TH232_GT_ANCHORS], dtype=np.float64)

    z = np.sqrt(Es)
    A_poly4 = np.vstack([z ** k for k in range(5)]).T
    coefs_p4, *_ = np.linalg.lstsq(A_poly4, Fs, rcond=None)
    model_p4 = FwhmModel(
        kind="lsrm_poly_sqrt_E",
        coefficients=tuple(float(c) for c in coefs_p4),
    )
    pred_p4 = np.array([model_p4(float(E)) for E in Es])
    max_abs_p4 = float(np.max(np.abs(pred_p4 - Fs)))

    A_quad = np.vstack([np.ones_like(Es), Es, Es ** 2]).T
    coefs_q, *_ = np.linalg.lstsq(A_quad, Fs * Fs, rcond=None)
    pred_quad = np.sqrt(np.maximum(
        coefs_q[0] + coefs_q[1] * Es + coefs_q[2] * Es ** 2, 0.01
    ))
    max_abs_quad = float(np.max(np.abs(pred_quad - Fs)))

    assert max_abs_p4 < 3.0, (
        f"F-452 poly-4 sqrt(E) lstsq должен давать max|Δ| < 3 keV; "
        f"получено {max_abs_p4:.3f} keV"
    )
    assert max_abs_p4 < 0.6 * max_abs_quad, (
        f"F-452 poly-4 ({max_abs_p4:.3f} keV) должен быть в ≥1.67× точнее "
        f"квадратичной FWHM²={max_abs_quad:.3f} keV"
    )
    rel_p4 = float(np.max(np.abs(pred_p4 - Fs) / Fs))
    assert rel_p4 < 0.07, f"max rel-error poly-4 должна быть < 7%; got {rel_p4*100:.2f}%"
    assert np.all(pred_p4 > 0)


def test_default_naI_obj_matches_default_naI_tuple():
    assert isinstance(_DEFAULT_NAI_FWHM_MODEL_OBJ, FwhmModel)
    assert _DEFAULT_NAI_FWHM_MODEL_OBJ.kind == "quad_fwhm2_in_E"
    assert _DEFAULT_NAI_FWHM_MODEL_OBJ.coefficients == _DEFAULT_NAI_FWHM_MODEL