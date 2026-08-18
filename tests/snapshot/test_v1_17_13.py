# -*- coding: utf-8 -*-
"""
v1.17.13 delivery tests — Deconvolution + peak shape.

Covers fixes F-280..F-283 (T-040 + T-079 + T-052 + T-003 + T-004).
T-073 (cross-correlation parallel search) is already present via
gamma.peaks.convolution_search; only its existence is asserted.

Run:
    cd 0_Work/gamma-spectrum-analysis
    $env:PYTHONPATH = "scripts"
    python -m pytest tests/snapshot/test_v1_17_13.py -v
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest

SCRIPTS = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts")
)
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)


# ─────────────────────────────────────────────────────────────────────
# F-280 / T-040 — smoothed-step background
# ─────────────────────────────────────────────────────────────────────

def test_F280_smoothed_step_bg_shape():
    """Step transitions from bg_left (слева) to bg_right (справа от пика)."""
    from gamma.peaks.background_options import smoothed_step_bg
    E = np.linspace(640.0, 690.0, 501)
    bg = smoothed_step_bg(E, E_peak=661.66, sigma_peak=10.0,
                          bg_left=200.0, bg_right=150.0)
    # Слева от пика → bg ≈ 200
    assert bg[0] == pytest.approx(200.0, abs=2.0)
    # Справа → bg ≈ 150
    assert bg[-1] == pytest.approx(150.0, abs=2.0)
    # На пике — среднее
    mid_idx = np.argmin(np.abs(E - 661.66))
    assert 170.0 < bg[mid_idx] < 180.0


def test_F280_smoothed_step_bg_with_slope():
    """Линейный наклон добавляется поверх ступеньки."""
    from gamma.peaks.background_options import smoothed_step_bg
    E = np.array([640.0, 661.66, 680.0])
    bg = smoothed_step_bg(E, E_peak=661.66, sigma_peak=5.0,
                          bg_left=100.0, bg_right=100.0,
                          slope=1.0)   # +1 counts/keV
    # На пике slope·0 = 0; слева slope·(-21.66) ≈ -22; справа +18.3
    assert bg[1] == pytest.approx(100.0, abs=1.0)
    assert bg[0] == pytest.approx(100.0 - 21.66, abs=1.0)
    assert bg[2] == pytest.approx(100.0 + 18.34, abs=1.0)


# ─────────────────────────────────────────────────────────────────────
# F-280 / T-079 — asymmetric bg regions
# ─────────────────────────────────────────────────────────────────────

def test_F280_symmetric_regions_when_no_interferent():
    """Без interferent — оба окна default ширины (1.0·FWHM)."""
    from gamma.peaks.background_options import asymmetric_bg_regions
    r = asymmetric_bg_regions(peak_channel=500, fwhm_channels=10.0)
    width_below = r.below_hi_ch - r.below_lo_ch
    width_above = r.above_hi_ch - r.above_lo_ch
    assert abs(width_below - width_above) <= 1  # симметрично


def test_F280_asymmetric_when_interferent_above():
    """Interferent above → above window narrow, below wide."""
    from gamma.peaks.background_options import asymmetric_bg_regions
    r = asymmetric_bg_regions(
        peak_channel=500,
        fwhm_channels=10.0,
        interferent_above_channel=525,    # 25 channels above = 2.5·FWHM
    )
    width_above = r.above_hi_ch - r.above_lo_ch
    width_below = r.below_hi_ch - r.below_lo_ch
    assert width_above < width_below
    assert r.interferent_above
    assert not r.interferent_below


def test_F280_asymmetric_bg_estimate():
    from gamma.peaks.background_options import (
        asymmetric_bg_regions, asymmetric_bg_estimate,
    )
    # Peak в районе 102 (каналы 100-104), вокруг чистый фон = 50.
    # FWHM=10 → gap=5, default_w=10. above_lo=107, above_hi=117 — все 50.
    counts = np.array([50.0] * 100 + [200.0] * 5 + [50.0] * 100)
    r = asymmetric_bg_regions(peak_channel=102, fwhm_channels=10.0)
    bg_left, bg_right = asymmetric_bg_estimate(counts, r)
    assert bg_left == pytest.approx(50.0, abs=1.0)
    assert bg_right == pytest.approx(50.0, abs=1.0)


# ─────────────────────────────────────────────────────────────────────
# F-281 / T-052 — sequential nonlinear-parameter exclusion
# ─────────────────────────────────────────────────────────────────────

def test_F281_low_dS_no_exclusion():
    from gamma.peaks.fit_stability import decide_nonlinear_exclusions
    d = decide_nonlinear_exclusions(0.02)
    assert not d.drop_step
    assert not d.drop_fwhm
    assert not d.drop_all_nl


def test_F281_medium_dS_drops_step():
    from gamma.peaks.fit_stability import decide_nonlinear_exclusions
    d = decide_nonlinear_exclusions(0.07)
    assert d.drop_step
    assert not d.drop_fwhm
    assert not d.drop_all_nl


def test_F281_high_dS_drops_step_and_fwhm():
    from gamma.peaks.fit_stability import decide_nonlinear_exclusions
    d = decide_nonlinear_exclusions(0.15)
    assert d.drop_step
    assert d.drop_fwhm
    assert not d.drop_all_nl


def test_F281_extreme_dS_drops_all():
    from gamma.peaks.fit_stability import decide_nonlinear_exclusions
    d = decide_nonlinear_exclusions(1.2)
    assert d.drop_step
    assert d.drop_fwhm
    assert d.drop_all_nl


# ─────────────────────────────────────────────────────────────────────
# F-282 / T-004 — NaI identification defaults
# ─────────────────────────────────────────────────────────────────────

def test_F282_nai_defaults_correct():
    from gamma.identification.nai_defaults import get_identification_defaults
    d = get_identification_defaults("NaI")
    assert d.lib_reduction is False
    assert d.peak_overlap_fwhm == 2.0
    assert d.use_tcc is False


def test_F282_hpge_defaults_differ():
    from gamma.identification.nai_defaults import get_identification_defaults
    nai = get_identification_defaults("NaI")
    hpge = get_identification_defaults("HPGe")
    assert hpge.lib_reduction is True
    assert hpge.peak_overlap_fwhm > nai.peak_overlap_fwhm
    assert hpge.use_tcc is True


def test_F282_unknown_falls_back_to_nai():
    from gamma.identification.nai_defaults import get_identification_defaults
    d = get_identification_defaults("MysteryDetector")
    # Falls back to NaI defaults (conservative)
    assert d.peak_overlap_fwhm == 2.0


# ─────────────────────────────────────────────────────────────────────
# F-283 / T-003 — pure-Gaussian policy guard
# ─────────────────────────────────────────────────────────────────────

def test_F283_nai_strict_returns_pure_gaussian():
    from gamma.peaks.pure_gaussian_policy import force_pure_gaussian_for
    p = force_pure_gaussian_for("NaI")
    assert p is not None
    assert p.tail_param == 0.0
    assert p.h_step == 0.0
    assert not p.use_T_E_model


def test_F283_hpge_no_strict_policy():
    from gamma.peaks.pure_gaussian_policy import force_pure_gaussian_for
    assert force_pure_gaussian_for("HPGe") is None


def test_F283_recommend_strict_lsrm_nai_pure():
    from gamma.peaks.pure_gaussian_policy import recommend_peak_shape
    tail, step, te = recommend_peak_shape("NaI", strict_lsrm=True)
    assert tail == 0.0
    assert step == 0.0
    assert not te


def test_F283_recommend_default_nai_has_tail():
    """Default (non-strict) — оставляет F-127 calibrated tail."""
    from gamma.peaks.pure_gaussian_policy import recommend_peak_shape
    tail, step, te = recommend_peak_shape("NaI", strict_lsrm=False)
    assert tail == 0.7
    assert step == 0.03
    assert te


# ─────────────────────────────────────────────────────────────────────
# T-073 — cross-correlation parallel peak search already present
# ─────────────────────────────────────────────────────────────────────

def test_T073_convolution_search_importable():
    """convolution_peak_search и compare_peak_methods доступны (F-124)."""
    from gamma.peaks.convolution_search import (
        convolution_peak_search, compare_peak_methods,
    )
    assert callable(convolution_peak_search)
    assert callable(compare_peak_methods)


def test_T073_staged_pipeline_compare_mode_available():
    """analyze_lsrm_spe принимает peak_search_method='compare'."""
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    import inspect
    sig = inspect.signature(analyze_lsrm_spe)
    assert "peak_search_method" in sig.parameters
