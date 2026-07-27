"""
Gauss + linear + erfc-step continuum fit (Gilmore section 9.7 / Lsrm section 9).

G1 fix (F-419, 2026-06-10): the skill `gamma-spectrum-analysis` Step 8/9
mandates "step-and-linear continuum (Gilmore section 9.7)" for isolated peak
integration. The previous Cowell-method baseline (linear wing-fit only)
systematically misses the Compton step under photopeaks sitting on a
strong continuum. Measured errors on Th-232 Marinelli poverka run
2026-06-03: Tl-208 2614 -38%, Tl-208 583 +45%, Ac-228 338 +39%,
Bi-212 727 missed entirely.

This module brings the same erfc-step continuum convention already used
in `deconvolve.py:_smooth_step` to the isolated single-peak path.

Returns `PeakAreaResult` from `area.py` for compatibility with
`integrate_peaks()` / `get_peak_area()` dispatch.
"""
from __future__ import annotations

import math
import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np

from .area import PeakAreaResult

SQRT_2PI = math.sqrt(2.0 * math.pi)
SQRT_2 = math.sqrt(2.0)
FWHM_TO_SIGMA = 2.355  # legacy project constant (≈ 2√(2·ln2))


def _free_sigma_enabled() -> bool:
    """F-449: free-σ opt-in toggle (analogue of LSRM «ПШПВ» flag ON)."""
    val = os.environ.get("GAMMA_FREE_SIGMA", "").strip().lower()
    return val in ("1", "true", "yes", "on")


def _erfc_vec(x):
    """Vectorised erfc — scipy if present, math.erfc fallback."""
    try:
        from scipy.special import erfc as _scipy_erfc

        return _scipy_erfc(x)
    except ImportError:
        from math import erfc as _math_erfc

        arr = np.asarray(x, dtype=np.float64)
        out = np.empty_like(arr)
        flat_in = arr.ravel()
        flat_out = out.ravel()
        for i in range(flat_in.size):
            flat_out[i] = _math_erfc(float(flat_in[i]))
        return out


def _model(x, H, c, sigma, a0, a1, h_step, x_mid):
    """Gauss + linear + erfc-step. Step centroid locked to peak centroid c."""
    gauss = H * np.exp(-((x - c) ** 2) / (2.0 * sigma * sigma))
    linear = a0 + a1 * (x - x_mid)
    step = h_step * 0.5 * _erfc_vec((x - c) / (sigma * SQRT_2))
    return gauss + linear + step


def gauss_erfc_step_fit(
    counts,
    *,
    peak_channel: int,
    fwhm_channels: float,
    window_factor: float = 2.5,
    sigma_bound_low: float = 0.6,
    sigma_bound_high: float = 1.6,
    fix_centroid: bool = False,
):
    """
    Compute peak net area with Gauss + linear + erfc-step continuum.
    G1 fix for skill Step 8/9 contract (Gilmore 9.7 / Lsrm 9).

    F-449 (operator-locked 2026-06-16): by default σ is HARD-LOCKED to
    sigma_cal = fwhm_channels/2.355 (LSRM «ПШПВ» flag OFF) and removed
    from the free parameters; reported FWHM == calibration FWHM. Set env
    GAMMA_FREE_SIGMA=1 to restore the legacy free-σ mode (σ bounded by
    [sigma_bound_low, sigma_bound_high]·sigma_cal); the sigma_bound_*
    args take effect ONLY in that mode.

    Free params (locked, default): [H, c, a0, a1, h_step]
    Free params (GAMMA_FREE_SIGMA=1): [H, c, sigma, a0, a1, h_step]
      H >= 0,  c in peak_channel +/- half_window/2 (or fixed),
      a0/a1 unconstrained,  h_step >= 0 (Compton step on low-E side).

    Area = H * sigma * sqrt(2*pi) — analytic Gaussian integral.
    Step is part of the continuum and NOT counted into net area.

    Returns PeakAreaResult with method="gauss_erfc_step".
    """
    n_ch = len(counts)
    if peak_channel < 0 or peak_channel >= n_ch:
        raise ValueError(f"peak_channel {peak_channel} outside [0, {n_ch})")
    if fwhm_channels <= 0:
        raise ValueError(f"fwhm_channels must be > 0, got {fwhm_channels}")

    sigma_cal = fwhm_channels / 2.355
    half_window = max(3, int(round(window_factor * fwhm_channels)))
    roi_low = max(0, peak_channel - half_window)
    roi_high = min(n_ch, peak_channel + half_window + 1)

    if roi_high - roi_low < 7:
        return PeakAreaResult(
            method="gauss_erfc_step",
            peak_channel=float(peak_channel),
            net_area_counts=0.0,
            net_area_uncertainty=0.0,
            gross_counts=0.0,
            baseline_counts=0.0,
            roi_low_ch=roi_low,
            roi_high_ch=roi_high,
            converged=False,
            notes=f"ROI too small for 6-param fit ({roi_high - roi_low} ch)",
        )

    roi_x = np.arange(roi_low, roi_high, dtype=np.float64)
    roi_y = np.asarray(counts[roi_low:roi_high], dtype=np.float64)
    x_mid = float(peak_channel)

    wing_size = max(2, (roi_high - roi_low) // 5)
    wing_x = np.concatenate([roi_x[:wing_size], roi_x[-wing_size:]])
    wing_y = np.concatenate([roi_y[:wing_size], roi_y[-wing_size:]])
    try:
        a1_init, a0_const_init = np.polyfit(wing_x, wing_y, 1)
        a0_init = float(a0_const_init + a1_init * x_mid)
        a1_init = float(a1_init)
    except np.linalg.LinAlgError:
        a0_init = float(wing_y.mean())
        a1_init = 0.0

    left_wing_mean = float(roi_y[:wing_size].mean())
    right_wing_mean = float(roi_y[-wing_size:].mean())
    h_step_init = max(0.0, left_wing_mean - right_wing_mean)

    H_init = max(1.0, float(roi_y.max()) - a0_init - 0.5 * h_step_init)

    try:
        from scipy.optimize import curve_fit
    except ImportError:
        return PeakAreaResult(
            method="gauss_erfc_step",
            peak_channel=float(peak_channel),
            net_area_counts=0.0,
            net_area_uncertainty=0.0,
            gross_counts=float(roi_y.sum()),
            baseline_counts=0.0,
            roi_low_ch=roi_low,
            roi_high_ch=roi_high,
            converged=False,
            notes="scipy.optimize unavailable; gauss_erfc_step needs scipy",
        )

    if fix_centroid:
        c_lo, c_hi = float(peak_channel) - 1e-6, float(peak_channel) + 1e-6
    else:
        c_lo = float(peak_channel) - half_window / 2.0
        c_hi = float(peak_channel) + half_window / 2.0

    # F-449: default σ-LOCK to calibration; free-σ only via GAMMA_FREE_SIGMA=1.
    free_sigma = _free_sigma_enabled()

    if free_sigma:
        sigma_lo = sigma_bound_low * sigma_cal
        sigma_hi = sigma_bound_high * sigma_cal
        if sigma_hi <= sigma_lo:
            sigma_hi = sigma_lo * 1.001
        p0 = [H_init, float(peak_channel), sigma_cal,
              a0_init, a1_init, h_step_init]
        bounds_lo = [0.0, c_lo, sigma_lo, -np.inf, -np.inf, 0.0]
        bounds_hi = [np.inf, c_hi, sigma_hi, np.inf, np.inf, np.inf]
    else:
        # σ hard-locked: free params [H, c, a0, a1, h_step].
        p0 = [H_init, float(peak_channel), a0_init, a1_init, h_step_init]
        bounds_lo = [0.0, c_lo, -np.inf, -np.inf, 0.0]
        bounds_hi = [np.inf, c_hi, np.inf, np.inf, np.inf]

    for i in range(len(p0)):
        if p0[i] < bounds_lo[i]:
            p0[i] = bounds_lo[i] + 1e-9
        elif p0[i] > bounds_hi[i]:
            p0[i] = bounds_hi[i] - 1e-9

    sigma_y = np.sqrt(np.maximum(roi_y, 1.0))

    def _model_partial(x, H, c, sig, a0, a1, h_s):
        return _model(x, H, c, sig, a0, a1, h_s, x_mid)

    def _fit_model_locked(x, H, c, a0, a1, h_s):
        return _model(x, H, c, sigma_cal, a0, a1, h_s, x_mid)

    _fit_model = _model_partial if free_sigma else _fit_model_locked

    try:
        popt, pcov = curve_fit(
            _fit_model,
            roi_x,
            roi_y,
            p0=p0,
            sigma=sigma_y,
            absolute_sigma=True,
            bounds=(bounds_lo, bounds_hi),
            maxfev=5000,
        )
    except Exception as e:
        return PeakAreaResult(
            method="gauss_erfc_step",
            peak_channel=float(peak_channel),
            net_area_counts=0.0,
            net_area_uncertainty=0.0,
            gross_counts=float(roi_y.sum()),
            baseline_counts=0.0,
            roi_low_ch=roi_low,
            roi_high_ch=roi_high,
            converged=False,
            notes=f"gauss_erfc_step fit failed: {type(e).__name__}: {e}",
        )

    if free_sigma:
        H_fit, c_fit, sigma_fit, a0_fit, a1_fit, h_step_fit = popt
        H_var = max(0.0, float(pcov[0, 0]))
        sig_var = max(0.0, float(pcov[2, 2]))
        H_sig_cov = float(pcov[0, 2])
        dA_dH = sigma_fit * SQRT_2PI
        dA_dsig = H_fit * SQRT_2PI
        area_var = (
            (dA_dH ** 2) * H_var
            + (dA_dsig ** 2) * sig_var
            + 2.0 * dA_dH * dA_dsig * H_sig_cov
        )
    else:
        # σ-locked: free params [H, c, a0, a1, h_step]; σ const = sigma_cal.
        H_fit, c_fit, a0_fit, a1_fit, h_step_fit = popt
        sigma_fit = sigma_cal
        H_var = max(0.0, float(pcov[0, 0]))
        # σ constant → area_var = (σ·√2π)²·H_var (F-449 contract).
        area_var = (sigma_fit * SQRT_2PI) ** 2 * H_var
    area_unc = math.sqrt(max(0.0, area_var))

    area = float(H_fit * sigma_fit * SQRT_2PI)

    continuum_per_ch = (
        a0_fit
        + a1_fit * (roi_x - x_mid)
        + h_step_fit * 0.5 * _erfc_vec((roi_x - c_fit) / (sigma_fit * SQRT_2))
    )
    baseline_counts = max(0.0, float(np.sum(continuum_per_ch)))
    gross_counts = float(roi_y.sum())

    y_model = _fit_model(roi_x, *popt)
    resid = (roi_y - y_model) / sigma_y
    n_params = 6 if free_sigma else 5
    dof = max(1, len(roi_y) - n_params)
    chi2_per_dof = float(np.sum(resid ** 2) / dof)

    # converged := fit produced positive height/area. We do NOT gate on
    # chi2/nu here: NaI photopeaks on heavy continua routinely produce
    # chi2/nu ~ 5-30 because the symmetric-Gauss model omits low-E
    # tailing. v2 hybrid step9 accepts such fits (chi2_red up to 30+).
    # Caller can inspect fit_residual_norm for a domain-specific gate.
    converged = bool(H_fit > 0 and area > 0)

    return PeakAreaResult(
        method="gauss_erfc_step",
        peak_channel=float(c_fit),
        net_area_counts=area,
        net_area_uncertainty=area_unc,
        gross_counts=gross_counts,
        baseline_counts=baseline_counts,
        roi_low_ch=roi_low,
        roi_high_ch=roi_high,
        fwhm_channels=float(sigma_fit * FWHM_TO_SIGMA),
        fit_height=float(H_fit),
        fit_residual_norm=chi2_per_dof,
        converged=converged,
        notes=(
            f"continuum=linear+erfc_step; h_step={h_step_fit:.2f}; "
            f"a0={a0_fit:.2f}, a1={a1_fit:.3f}; "
            f"sigma_cal={sigma_cal:.3f} fit_sigma={sigma_fit:.3f}; "
            f"sigma_mode={'free' if free_sigma else 'locked(F-449)'}"
        ),
    )


__all__ = ["gauss_erfc_step_fit"]