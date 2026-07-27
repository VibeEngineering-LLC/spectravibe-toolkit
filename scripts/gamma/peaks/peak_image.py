"""
LSRM-style peak-image model — Gaussian + low-energy exponential tail
+ Compton step (LSRM Algorithmic Foundations §8.4.2.1 and §8.4.4).

Background
----------
A pure Gaussian is a good first-order model for the FEP shape in NaI
but not for the low-energy side of the peak nor for the continuum
immediately under the peak. LSRM exposes a three-component description
of the peak-image:

    F(x) = peak(x) + step(x) + bg(x)

where:

    peak(x) = Gaussian + low-energy exponential tail
    step(x) = Compton-step (erfc-shaped) coupling left-side continuum
              to the peak height (charge loss, small-angle Compton
              scattering escape — LSRM §8.4.4)
    bg(x)   = polynomial background (handled outside this module)

The peak-image:

           ┌  A·exp(−(x−μ)²/(2σ²))                        if x ≥ μ − T·σ
    f(x) = ┤
           └  A·exp(T·(x−μ)/σ + T²/2)                     if x < μ − T·σ

The tail-coupling point is μ − T·σ; for x to the LEFT of this point
the curve switches from Gaussian to a left-going exponential whose
slope is set by σ and amplitude continuity at the join. T is the
"tail parameter" in fractions of FWHM (LSRM, "хвостовой параметр").
Smaller T → stronger / longer tail. For NaI a typical T is 0.5–0.9.

The Compton step:

    step(x) = (A_step / 2) · erfc((x − μ) / (σ·√2))

This is the LSRM-Gilmore convention — half the peak height multiplied
by the complementary error function centered on μ. It produces a smooth
"step down" from the high left-continuum to the low right-continuum
across the peak.

A_step is parameterised as a fraction of the peak height:

    A_step = h_step_frac · A      with typical h_step_frac ∈ [0, 0.05]

Why this matters for Gamma-1S (NaI 63×63)
-----------------------------------------
At E < 200 keV the Gaussian-only fit biases peak-area by ~5–8 %
because:

  * The low-energy tail is unmodelled → fit pulls the FWHM wider
    to compensate, mis-attributing tail counts to the Gaussian core.
  * The Compton step is unmodelled → polynomial background takes its
    place and bends the background fit, shifting the apparent area.

This module returns model values directly and exposes a curve_fit
wrapper `fit_peak_image()` that takes a ROI and returns area,
position, sigma, T, h_step plus their covariance.

Calibration (T(E), h_step(E))
-----------------------------
T is mildly energy-dependent — LSRM recommends calibrating it from
a reference source covering ≥3 energies. We expose
`calibrate_T_of_E()` which fits a linear model T(E) = T0 + T1·E.
For h_step a constant is usually sufficient unless K-edge effects
intrude (not relevant on NaI 63×63 above ~50 keV).

Reference: LSRM Algorithmic Foundations 2022, §8.4 «Калибровка по
форме линии. Пик-образ» (стр. 8-3 — 8-7), Figure 8.4 «Параметры
пика и хвостовой параметр», Figure 8.5 «Комптоновская ступенька».
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import Optional, Tuple, Sequence, Callable

import numpy as np

# -- 1/√(2π) and √2 reused constants
_SQRT_TWO = math.sqrt(2.0)
_ERFC_AVAILABLE = hasattr(math, "erfc")
_FWHM_TO_SIGMA = 2.3548  # 2√(2·ln2)


def _free_sigma_enabled() -> bool:
    """F-449: free-σ opt-in toggle (analogue of LSRM «ПШПВ» flag ON)."""
    val = os.environ.get("GAMMA_FREE_SIGMA", "").strip().lower()
    return val in ("1", "true", "yes", "on")


def _erfc(x: float) -> float:
    """erfc with fallback for environments without math.erfc."""
    if _ERFC_AVAILABLE:
        return math.erfc(x)
    # Abramowitz & Stegun 7.1.26 — error < 1.5e-7 for x ≥ 0
    sign = 1.0 if x >= 0 else -1.0
    x = abs(x)
    a1, a2, a3, a4, a5 = 0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429
    p = 0.3275911
    t = 1.0 / (1.0 + p * x)
    y = (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t
    return 1.0 - sign * (1.0 - y * math.exp(-x * x))


# ---------------------------------------------------------------------------
# Functional form
# ---------------------------------------------------------------------------

def gaussian(x: np.ndarray, A: float, mu: float, sigma: float) -> np.ndarray:
    """Pure Gaussian, area = A · σ · √(2π) when integrated over R."""
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def gaussian_with_tail(
    x: np.ndarray, A: float, mu: float, sigma: float, T: float
) -> np.ndarray:
    """
    Gaussian + exponential low-energy tail (LSRM Figure 8.4).

    Parameters
    ----------
    x : np.ndarray         channels or keV
    A : peak height at x = μ
    mu : peak position
    sigma : Gaussian σ
    T : tail join distance in σ-units (smaller = stronger tail)

    The Gaussian branch holds for x ≥ μ − T·σ; for x < μ − T·σ the
    function switches to a left-going exponential whose slope is
    chosen so the function and its first derivative are continuous
    at the join point.
    """
    if T <= 0:
        # Pathological — fall back to pure Gaussian (avoid div by zero)
        return gaussian(x, A, mu, sigma)
    z = (x - mu) / sigma
    # Clamp the tail-branch exponent to avoid overflow at extreme z<<−T
    tail_arg = np.clip(T * z + 0.5 * T * T, -700.0, 700.0)
    out = np.where(
        z >= -T,
        A * np.exp(-0.5 * z * z),
        A * np.exp(tail_arg),
    )
    return out


def compton_step(
    x: np.ndarray, A_step: float, mu: float, sigma: float
) -> np.ndarray:
    """
    Compton step (LSRM Figure 8.5) — erfc descent centred on μ.

    A_step is the FULL step height (left-continuum − right-continuum
    asymptotic difference). The function rises smoothly from 0 (right
    of the peak) to A_step (left of the peak):

        step(x) = (A_step / 2) · erfc((x − μ) / (σ · √2))
    """
    if sigma <= 0 or A_step == 0:
        return np.zeros_like(x)
    # vectorise erfc
    arg = (x - mu) / (sigma * _SQRT_TWO)
    return 0.5 * A_step * _vectorised_erfc(arg)


def _vectorised_erfc(arr: np.ndarray) -> np.ndarray:
    """numpy-friendly erfc via scipy if available, else math.erfc."""
    try:
        from scipy.special import erfc  # type: ignore
        return erfc(arr)
    except Exception:
        return np.vectorize(_erfc)(arr)


def peak_image(
    x: np.ndarray,
    A: float, mu: float, sigma: float,
    T: float = 0.7, h_step_frac: float = 0.0,
) -> np.ndarray:
    """
    Full LSRM peak-image = Gaussian + tail + Compton step.

    The step amplitude is parameterised as a fraction of the peak
    height: A_step = h_step_frac · A.
    """
    peak = gaussian_with_tail(x, A, mu, sigma, T)
    step = compton_step(x, h_step_frac * A, mu, sigma)
    return peak + step


# ---------------------------------------------------------------------------
# Area integration (closed-form for Gaussian; numeric for tail)
# ---------------------------------------------------------------------------

def integrated_area(
    A: float, sigma: float, T: float = 0.7, bin_w: float = 1.0,
) -> float:
    """
    Closed-form integral of the Gaussian-with-tail (no step).

    For the Gaussian branch:
        ∫ A·exp(−z²/2) dx = A·σ·√(2π)   (full Gaussian)
    For the exponential tail branch (x ∈ (−∞, μ − T·σ)):
        ∫ A·exp(T·z + T²/2)·σ dz = A·σ·exp(−T²/2) / T

    The combined area is obtained by replacing the truncated-Gaussian
    contribution on the left of −T with the tail integral.

    Units / ``bin_w`` (F-271 — v1.17.11, T-017)
    -------------------------------------------
    Curve-fit returns ``A`` in **the same units as ``y``** (counts per
    channel for SPE spectra) and ``sigma`` in **the same units as
    ``x``**. When the fit is performed on (channels, counts/channel),
    σ is in channels and the analytic integral ``A·σ·√(2π)`` is already
    in counts (channels cancel against per-channel A).

    BUT when the fit is performed on (keV, counts/channel) — which is
    the case for energy-axis multiplet wrappers — σ is in keV and the
    raw integral has units of [counts · keV / channel]. Dividing by
    ``bin_w`` (keV per channel) recovers true counts.

    Pass ``bin_w = 1.0`` (default) when σ is in channels — preserves
    legacy behaviour for callers already in channel space.
    Pass ``bin_w = mean(diff(E_keV))`` when σ is in keV.
    """
    if sigma <= 0 or A <= 0:
        return 0.0
    # Gaussian contribution on z ≥ −T:
    #     A · σ · ∫_{−T}^{+∞} exp(−z²/2) dz  =  A·σ·√(π/2)·(1 + erf(T/√2))
    # Note: 1 + erf(T/√2) = 2 − erfc(T/√2)
    gauss_right = A * sigma * math.sqrt(math.pi / 2.0) * (2.0 - _erfc(T / _SQRT_TWO))
    # Exponential tail on z < −T:  A·exp(T·z + T²/2)
    #     ∫_{−∞}^{−T} exp(T·z + T²/2) dz = exp(T·(−T) + T²/2) / T = exp(−T²/2) / T
    if T > 0:
        tail = A * sigma * math.exp(-0.5 * T * T) / T
    else:
        tail = 0.0
    raw = gauss_right + tail
    if bin_w and bin_w > 0 and bin_w != 1.0:
        return raw / float(bin_w)
    return raw


# ---------------------------------------------------------------------------
# Fit wrapper (scipy.optimize.curve_fit, soft-fail on absence)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PeakImageFitResult:
    """Outcome of fit_peak_image()."""
    A: float
    mu: float
    sigma: float
    T: float
    h_step_frac: float
    area: float                          # integrated Gaussian+tail
    area_sigma: float                    # σ propagated from cov
    cov: tuple = ()                      # 5×5 cov as tuple of tuples
    chi2_per_dof: Optional[float] = None
    converged: bool = True
    notes: str = ""

    def __repr__(self) -> str:
        return (
            f"PeakImageFit(μ={self.mu:.2f}, σ={self.sigma:.3f}, "
            f"T={self.T:.2f}, step={self.h_step_frac:.3f}, "
            f"area={self.area:.0f}±{self.area_sigma:.0f}, "
            f"χ²/dof={self.chi2_per_dof})"
        )


def fit_peak_image(
    x: Sequence[float],
    y: Sequence[float],
    *,
    A0: Optional[float] = None,
    mu0: Optional[float] = None,
    sigma0: Optional[float] = None,
    sigma_fixed: Optional[float] = None,
    fwhm_channels: Optional[float] = None,
    fit_sigma: Optional[bool] = None,
    T0: float = 0.7,
    h_step0: float = 0.0,
    fit_T: bool = True,
    fit_step: bool = True,
    sigma_y: Optional[Sequence[float]] = None,
) -> PeakImageFitResult:
    """
    Fit the LSRM peak-image to a single peak ROI.

    Inputs are channel-or-keV abscissa and counts. Background
    (polynomial baseline) is expected to be subtracted upstream — we
    only fit peak + Compton step on top of zero baseline. If the
    Compton step is requested without a separate baseline subtraction,
    the step amplitude will absorb a constant offset; that is the
    intended behaviour for the immediate vicinity of the peak.

    Initial guesses (A0, mu0, sigma0) come from the simple-moment
    estimator if not supplied. T0/h_step0 default to NaI-typical
    values.

    F-449 (operator-locked 2026-06-16): σ is DETERMINED by the FWHM(E)
    calibration, not fitted per-peak (LSRM «ПШПВ» flag OFF). Pass the
    calibration width via ``sigma_fixed`` (x-units) or ``fwhm_channels``
    (→ σ = fwhm/2.3548); σ is then HARD-LOCKED (free: A, μ, T, h_step).
    ``fit_sigma``: None → env GAMMA_FREE_SIGMA decides (default locked
    when a width is given); True → force free σ; False → force locked.
    No width + no override → σ seeded from half-max & fitted (pre-F-449).

    Returns PeakImageFitResult with `converged=False` and `notes`
    set on failure; never raises.
    """
    try:
        from scipy.optimize import curve_fit  # type: ignore
    except Exception:
        return PeakImageFitResult(
            A=float("nan"), mu=float("nan"), sigma=float("nan"),
            T=T0, h_step_frac=h_step0,
            area=float("nan"), area_sigma=float("nan"),
            converged=False,
            notes="scipy.optimize.curve_fit unavailable",
        )

    x_arr = np.asarray(list(x), dtype=float)
    y_arr = np.asarray(list(y), dtype=float)
    if len(x_arr) < 5 or len(x_arr) != len(y_arr):
        return PeakImageFitResult(
            A=float("nan"), mu=float("nan"), sigma=float("nan"),
            T=T0, h_step_frac=h_step0,
            area=float("nan"), area_sigma=float("nan"),
            converged=False, notes="ROI too small or x/y mismatch",
        )

    # ── Initial guesses ─────────────────────────────────────────────
    if A0 is None:
        A0 = float(np.max(y_arr))
    if mu0 is None:
        mu0 = float(x_arr[int(np.argmax(y_arr))])
    if sigma0 is None:
        # Estimate from FWHM by finding y > A0/2 channels
        above = y_arr >= 0.5 * A0
        if above.sum() >= 2:
            x_above = x_arr[above]
            fwhm = float(x_above[-1] - x_above[0])
            sigma0 = max(fwhm / 2.355, (x_arr[1] - x_arr[0]) if len(x_arr) > 1 else 1.0)
        else:
            sigma0 = max(float(np.std(x_arr)) / 3.0, 1.0)

    # ── F-449: resolve calibration σ and lock/free decision ─────────
    sigma_cal = None
    if sigma_fixed is not None and sigma_fixed > 0:
        sigma_cal = float(sigma_fixed)
    elif fwhm_channels is not None and fwhm_channels > 0:
        sigma_cal = float(fwhm_channels) / _FWHM_TO_SIGMA
    have_cal_width = sigma_cal is not None and sigma_cal > 0
    if fit_sigma is True:
        lock_sigma = False
    elif fit_sigma is False:
        lock_sigma = have_cal_width  # cannot lock without a width
    else:  # None → env-controlled
        lock_sigma = have_cal_width and not _free_sigma_enabled()
    if lock_sigma:
        sigma0 = float(sigma_cal)  # seed locked value

    # ── Build model with possibly-fixed T and/or step ──────────────
    def make_model(fit_T: bool, fit_step: bool):
        if fit_T and fit_step:
            def model(x, A, mu, sigma, T, h):
                return peak_image(x, A, mu, sigma, T, h)
            p0 = [A0, mu0, sigma0, T0, h_step0]
            return model, p0
        if fit_T and not fit_step:
            def model(x, A, mu, sigma, T):
                return peak_image(x, A, mu, sigma, T, h_step0)
            p0 = [A0, mu0, sigma0, T0]
            return model, p0
        if not fit_T and fit_step:
            def model(x, A, mu, sigma, h):
                return peak_image(x, A, mu, sigma, T0, h)
            p0 = [A0, mu0, sigma0, h_step0]
            return model, p0
        # neither fit
        def model(x, A, mu, sigma):
            return peak_image(x, A, mu, sigma, T0, h_step0)
        p0 = [A0, mu0, sigma0]
        return model, p0

    model, p0 = make_model(fit_T, fit_step)

    try:
        popt, pcov = curve_fit(
            model, x_arr, y_arr, p0=p0,
            sigma=(np.asarray(list(sigma_y), dtype=float)
                   if sigma_y is not None else None),
            absolute_sigma=(sigma_y is not None),
            maxfev=5000,
        )
    except Exception as exc:
        return PeakImageFitResult(
            A=float("nan"), mu=float("nan"), sigma=float("nan"),
            T=T0, h_step_frac=h_step0,
            area=float("nan"), area_sigma=float("nan"),
            converged=False, notes=f"curve_fit failed: {exc}",
        )

    # ── Unpack ──────────────────────────────────────────────────────
    if fit_T and fit_step:
        A, mu, sigma, T, h = popt
    elif fit_T and not fit_step:
        A, mu, sigma, T = popt
        h = h_step0
    elif not fit_T and fit_step:
        A, mu, sigma, h = popt
        T = T0
    else:
        A, mu, sigma = popt
        T, h = T0, h_step0

    # F-271 (v1.17.11, T-017) — авто-bin_w из x-сетки.
    # Если x — каналы (диффы ≈ 1), bin_w=1 → результат без изменений
    # (backward-compat). Если x — keV, bin_w = средний шаг → результат
    # в counts (а не counts·keV/channel).
    if len(x_arr) >= 2:
        dx = float(np.mean(np.diff(x_arr)))
        bin_w = dx if dx > 0 else 1.0
    else:
        bin_w = 1.0
    area = integrated_area(A, sigma, T, bin_w=bin_w)
    # σ(area) ≈ |∂area/∂A|·σ_A + |∂area/∂σ|·σ_σ ; we use diag(pcov) only
    try:
        sigma_A_par = math.sqrt(max(0.0, float(pcov[0, 0])))
        sigma_sigma = math.sqrt(max(0.0, float(pcov[2, 2])))
    except Exception:
        sigma_A_par = float("nan")
        sigma_sigma = float("nan")
    # Approximate area uncertainty by relative quadrature on A and σ
    if A > 0 and sigma > 0:
        rel_A = sigma_A_par / A if math.isfinite(sigma_A_par) else 0.0
        rel_sig = sigma_sigma / sigma if math.isfinite(sigma_sigma) else 0.0
        area_sigma = abs(area) * math.sqrt(rel_A ** 2 + rel_sig ** 2)
    else:
        area_sigma = float("nan")

    # ── χ²/dof ──────────────────────────────────────────────────────
    try:
        resid = y_arr - model(x_arr, *popt)
        if sigma_y is not None:
            sy = np.asarray(list(sigma_y), dtype=float)
            chi2 = float(np.sum((resid / np.where(sy > 0, sy, 1.0)) ** 2))
        else:
            # Poisson approximation: σ = √y on bins with y > 0
            sy = np.where(y_arr > 0, np.sqrt(y_arr), 1.0)
            chi2 = float(np.sum((resid / sy) ** 2))
        dof = max(1, len(x_arr) - len(popt))
        chi2_dof = chi2 / dof
    except Exception:
        chi2_dof = None

    return PeakImageFitResult(
        A=float(A), mu=float(mu), sigma=float(sigma),
        T=float(T), h_step_frac=float(h),
        area=float(area), area_sigma=float(area_sigma),
        cov=tuple(tuple(row) for row in pcov.tolist()) if pcov is not None else (),
        chi2_per_dof=chi2_dof,
        converged=True,
        notes="",
    )


# ---------------------------------------------------------------------------
# T(E) calibration
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TailCalibration:
    """T(E) = T0 + T1·E linear calibration of the tail parameter."""
    T0: float
    T1: float
    energies_keV: tuple = ()
    T_values: tuple = ()
    notes: str = ""

    def __call__(self, E_keV: float) -> float:
        return max(0.05, self.T0 + self.T1 * E_keV)


def calibrate_T_of_E(
    energies_keV: Sequence[float],
    T_values: Sequence[float],
) -> TailCalibration:
    """
    Fit T(E) = T0 + T1·E from ≥2 reference points.

    For NaI 63×63 a typical calibration uses isolated lines at 661.66
    keV (Cs-137) and 1460.82 keV (K-40). Add a low-energy anchor
    (e.g. Am-241 59.54 keV) for better leverage at the K-edge region.

    Returns a TailCalibration callable T(E).
    """
    if len(energies_keV) != len(T_values) or len(energies_keV) < 2:
        # Fall back to a NaI-typical constant T = 0.7
        return TailCalibration(
            T0=0.7, T1=0.0,
            energies_keV=tuple(energies_keV),
            T_values=tuple(T_values),
            notes="insufficient points, falling back to T=0.7",
        )
    E = np.asarray(list(energies_keV), dtype=float)
    T = np.asarray(list(T_values), dtype=float)
    # Simple linear regression
    A = np.column_stack([np.ones_like(E), E])
    coef, *_ = np.linalg.lstsq(A, T, rcond=None)
    return TailCalibration(
        T0=float(coef[0]), T1=float(coef[1]),
        energies_keV=tuple(E.tolist()),
        T_values=tuple(T.tolist()),
        notes="linear T(E) fit",
    )


__all__ = [
    "gaussian",
    "gaussian_with_tail",
    "compton_step",
    "peak_image",
    "integrated_area",
    "fit_peak_image",
    "PeakImageFitResult",
    "TailCalibration",
    "calibrate_T_of_E",
]
