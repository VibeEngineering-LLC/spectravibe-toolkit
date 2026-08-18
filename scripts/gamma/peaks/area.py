"""
Peak area integration — converts raw spectrum counts under a peak into
a net peak area S (counts above background), the fundamental quantity
for activity calculations.

Two methods are implemented:

  1. **Cowell method** (Lsrm §5.2.5) — simplest, robust for low-stat
     spectra. Integrate raw counts in a window around the peak, subtract
     polynomial baseline fit to the wings.

         S = N_total - B_polynomial
         dS = √N_total  (Poisson)

     Pros: works on any peak shape, no fit failures.
     Cons: includes systematic bias when peak shape differs from
     symmetric (low-E tails on HPGe, Compton step on scintillators).

  2. **Gaussian fit** (Lsrm §5.2, simplified) — fit a Gaussian peak
     plus a linear (or quadratic) baseline by least-squares.

         y_i = H · exp(-(x_i - c)² / (2σ²)) + (a₀ + a₁·x_i)

     Area = H · σ · √(2π).

     Pros: more accurate for resolved peaks with known shape.
     Cons: fails on multiplets (use multi_gaussian_fit instead),
     fails on highly asymmetric peaks (HPGe needs tail).

For multiplets (overlapping peaks), use `multi_gaussian_fit` from
the deconvolution module (Phase 2.1b).

References:
  - Lsrm Algorithmic Foundations 2022, §5.1 (simple methods) and §5.2
    (model-based fitting), §5.2.5 (Cowell)
  - Gilmore & Joss, Practical Gamma-ray Spectrometry, Ch. 5
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


@dataclass(frozen=True)
class PeakAreaResult:
    """Result of peak area integration for one peak."""

    method: str                     # "cowell" or "gaussian_fit"
    peak_channel: float             # final center
    net_area_counts: float          # S — net counts in the peak
    net_area_uncertainty: float     # dS — propagated uncertainty
    gross_counts: float             # N_total — sum of raw counts in ROI
    baseline_counts: float          # B — baseline contribution under peak

    # Window used (for diagnostics)
    roi_low_ch: int
    roi_high_ch: int

    # Method-specific
    fwhm_channels: Optional[float] = None  # fitted FWHM (gaussian only)
    fit_height: Optional[float] = None     # peak height (gaussian only)
    fit_residual_norm: Optional[float] = None  # χ²/dof
    converged: bool = True

    notes: str = ""

    def __repr__(self) -> str:
        return (f"PeakAreaResult({self.method}: ch={self.peak_channel:.1f}, "
                f"S={self.net_area_counts:.0f}±{self.net_area_uncertainty:.0f}, "
                f"gross={self.gross_counts:.0f})")


def cowell_area(
    counts: np.ndarray,
    *,
    peak_channel: int,
    fwhm_channels: float,
    window_factor: float = 2.5,
    baseline_polynomial_order: int = 1,
    wing_fraction: float = 0.3,
) -> PeakAreaResult:
    """
    Compute peak net area by Cowell's method.

    Algorithm (Lsrm §5.2.5):
      1. Define ROI: peak_channel ± window_factor·FWHM
      2. Define baseline wings: first `wing_fraction` and last
         `wing_fraction` of the ROI
      3. Fit polynomial of given order to baseline wings
      4. Integrate polynomial across full ROI → B (baseline counts)
      5. Sum raw counts in ROI → N_total (gross counts)
      6. Net area S = N_total - B
      7. Uncertainty dS = √N_total (Poisson, dominant for sample-bg)

    The Cowell method is robust because it makes no assumption about
    peak shape — it works for Gaussian, asymmetric, or even
    fractionally-resolved peaks. Its main limitation is that it
    cannot separate overlapping peaks (for those use multi_gaussian_fit).

    Args:
        counts: spectrum array indexed by channel
        peak_channel: integer channel of the peak centroid
        fwhm_channels: full-width at half-max in channels at this energy
        window_factor: ROI half-width = window_factor · FWHM. Default
            2.5 captures ≥99% of a Gaussian peak.
        baseline_polynomial_order: 0=constant, 1=linear (default),
            2=quadratic. Linear is appropriate for most Compton
            backgrounds in NaI spectra.
        wing_fraction: fraction of ROI used on each side for baseline
            fit. Default 0.3 (each wing) leaves 0.4 of ROI as the
            "peak zone" where baseline is extrapolated.

    Returns:
        PeakAreaResult.

    Example (synthetic Gaussian on flat background):
        >>> import numpy as np
        >>> N = 1000
        >>> ch = np.arange(N)
        >>> gauss = 5000 * np.exp(-((ch - 500) / 10)**2 / 2)
        >>> counts = (100 + gauss).astype(int)
        >>> r = cowell_area(counts, peak_channel=500, fwhm_channels=23.55)
        >>> # FWHM = 2.355 · σ, so for σ=10 FWHM ≈ 23.55
        >>> # Expected area = 5000 · 10 · √(2π) ≈ 125331
        >>> bool(abs(r.net_area_counts - 125331) / 125331 < 0.05)
        True
    """
    n_ch = len(counts)
    if peak_channel < 0 or peak_channel >= n_ch:
        raise ValueError(f"peak_channel {peak_channel} outside spectrum [0, {n_ch})")
    if fwhm_channels <= 0:
        raise ValueError(f"fwhm_channels must be > 0, got {fwhm_channels}")

    half_window = max(1, int(round(window_factor * fwhm_channels)))
    roi_low = max(0, peak_channel - half_window)
    roi_high = min(n_ch, peak_channel + half_window + 1)

    if roi_high - roi_low < 3:
        # ROI too small for baseline + integration
        return PeakAreaResult(
            method="cowell",
            peak_channel=float(peak_channel),
            net_area_counts=0.0,
            net_area_uncertainty=0.0,
            gross_counts=0.0,
            baseline_counts=0.0,
            roi_low_ch=roi_low,
            roi_high_ch=roi_high,
            converged=False,
            notes=f"ROI too small ({roi_high - roi_low} channels < 3)",
        )

    roi = counts[roi_low:roi_high].astype(np.float64)
    roi_x = np.arange(roi_low, roi_high, dtype=np.float64)
    roi_width = roi_high - roi_low

    # Identify baseline wings: first `wing_fraction` and last
    # `wing_fraction` of the ROI.
    wing_size = max(2, int(round(wing_fraction * roi_width)))
    left_wing_idx = np.arange(0, wing_size)
    right_wing_idx = np.arange(roi_width - wing_size, roi_width)
    wing_idx = np.concatenate([left_wing_idx, right_wing_idx])
    wing_x = roi_x[wing_idx]
    wing_y = roi[wing_idx]

    # Fit polynomial to wings
    try:
        coefs = np.polyfit(wing_x, wing_y, baseline_polynomial_order)
    except (np.linalg.LinAlgError, ValueError) as e:
        return PeakAreaResult(
            method="cowell",
            peak_channel=float(peak_channel),
            net_area_counts=0.0,
            net_area_uncertainty=0.0,
            gross_counts=float(roi.sum()),
            baseline_counts=0.0,
            roi_low_ch=roi_low,
            roi_high_ch=roi_high,
            converged=False,
            notes=f"Baseline polyfit failed: {e}",
        )

    # Evaluate baseline across full ROI and integrate (sum)
    baseline_at_roi = np.polyval(coefs, roi_x)
    baseline_total = float(baseline_at_roi.sum())
    baseline_total = max(0.0, baseline_total)  # don't allow negative

    gross_total = float(roi.sum())
    net_area = gross_total - baseline_total
    # Poisson uncertainty: dominant term is √(gross + B·(t/t_bg)²);
    # for in-spectrum baseline subtraction, dS ≈ √(N_total + B) ≈ √N_total
    # because B ~ N_total dominates when there's a continuum.
    # Conservative: include both:
    net_uncertainty = math.sqrt(gross_total + baseline_total)

    # Sanity check: detect when the baseline extrapolation has gone
    # wrong because the wings contain peak features (Compton step,
    # neighboring peak, escape peak). Symptoms:
    #   • net_area is negative (baseline higher than actual counts)
    #   • OR wings have higher mean than peak center (baseline slopes
    #     wrong direction relative to peak)
    converged = True
    notes_extra = ""
    if net_area < 0:
        converged = False
        notes_extra = (f"negative net area ({net_area:.0f}) — likely "
                       f"baseline extrapolation error from peak features "
                       f"in wings (neighboring peak, Compton step, escape)")
    else:
        # Check: mean of wings vs counts at peak center.
        # If peak center is BELOW the wings, this is not a real peak
        # (likely just a dip in continuum) and area is meaningless.
        peak_center_count = roi[len(roi) // 2]
        wing_mean = wing_y.mean()
        if peak_center_count < wing_mean * 0.7:
            converged = False
            notes_extra = (f"peak center counts ({peak_center_count:.0f}) much "
                           f"lower than wing mean ({wing_mean:.0f}); ROI may "
                           f"contain a continuum dip rather than peak")

    return PeakAreaResult(
        method="cowell",
        peak_channel=float(peak_channel),
        net_area_counts=net_area,
        net_area_uncertainty=net_uncertainty,
        gross_counts=gross_total,
        baseline_counts=baseline_total,
        roi_low_ch=roi_low,
        roi_high_ch=roi_high,
        converged=converged,
        notes=(f"baseline poly order={baseline_polynomial_order}, "
               f"wings={wing_size} each"
               + (f"; {notes_extra}" if notes_extra else "")),
    )


def _gaussian_plus_baseline(x, H, c, sigma, a0, a1):
    """Gaussian peak + linear baseline."""
    return H * np.exp(-(x - c)**2 / (2.0 * sigma**2)) + a0 + a1 * x


def gaussian_fit_area(
    counts: np.ndarray,
    *,
    peak_channel: int,
    fwhm_channels: float,
    window_factor: float = 3.0,
    fix_centroid: bool = False,
    fix_fwhm: bool = False,
) -> PeakAreaResult:
    """
    Compute peak area by Gaussian + linear-baseline fit.

    Algorithm:
      1. ROI = peak_channel ± window_factor·FWHM
      2. Initial guesses:
         - H = max(counts in ROI) - mean(wing counts)
         - c = peak_channel
         - σ = fwhm_channels / 2.355
         - a₀, a₁ from linear fit to wings
      3. Non-linear least-squares fit of Gaussian + linear baseline
         using scipy.optimize (or numpy fallback if unavailable)
      4. Area = H · σ · √(2π)
      5. Uncertainty from covariance matrix or Poisson approximation

    Args:
        counts: spectrum array indexed by channel
        peak_channel: integer channel of the peak centroid (initial)
        fwhm_channels: initial FWHM in channels
        window_factor: ROI half-width = factor · FWHM. Default 3.
        fix_centroid: if True, holds centroid at peak_channel
            (use when calibration is known to be reliable).
        fix_fwhm: if True, holds FWHM at fwhm_channels (use when
            FWHM calibration is reliable).

    Returns:
        PeakAreaResult with fitted parameters.

    Notes:
        If scipy is not available, falls back to a simple bracketed
        search for H and uses analytical formulas for σ, c via moments.
        This is less accurate but does not require scipy.

    Example:
        >>> import numpy as np
        >>> N = 1000
        >>> ch = np.arange(N)
        >>> gauss = 5000 * np.exp(-((ch - 500) / 10)**2 / 2)
        >>> counts = (100 + gauss).astype(int)
        >>> r = gaussian_fit_area(counts, peak_channel=500, fwhm_channels=23.55)
        >>> bool(abs(r.net_area_counts - 125331) / 125331 < 0.02)
        True
    """
    n_ch = len(counts)
    if peak_channel < 0 or peak_channel >= n_ch:
        raise ValueError(f"peak_channel {peak_channel} outside [0, {n_ch})")
    if fwhm_channels <= 0:
        raise ValueError(f"fwhm_channels must be > 0, got {fwhm_channels}")

    sigma_init = fwhm_channels / 2.355
    half_window = max(3, int(round(window_factor * fwhm_channels)))
    roi_low = max(0, peak_channel - half_window)
    roi_high = min(n_ch, peak_channel + half_window + 1)

    if roi_high - roi_low < 6:
        return PeakAreaResult(
            method="gaussian_fit",
            peak_channel=float(peak_channel),
            net_area_counts=0.0,
            net_area_uncertainty=0.0,
            gross_counts=0.0,
            baseline_counts=0.0,
            roi_low_ch=roi_low,
            roi_high_ch=roi_high,
            converged=False,
            notes=f"ROI too small for Gaussian fit ({roi_high - roi_low} ch)",
        )

    roi_x = np.arange(roi_low, roi_high, dtype=np.float64)
    roi_y = counts[roi_low:roi_high].astype(np.float64)

    # Initial baseline from wings
    wing_size = max(2, (roi_high - roi_low) // 5)
    wing_x = np.concatenate([roi_x[:wing_size], roi_x[-wing_size:]])
    wing_y = np.concatenate([roi_y[:wing_size], roi_y[-wing_size:]])
    try:
        a1_init, a0_init = np.polyfit(wing_x, wing_y, 1)
    except np.linalg.LinAlgError:
        a0_init = float(wing_y.mean())
        a1_init = 0.0

    # Initial height: max counts - baseline at center
    baseline_at_peak = a0_init + a1_init * peak_channel
    H_init = float(roi_y.max()) - baseline_at_peak
    H_init = max(1.0, H_init)

    # Try scipy.optimize.curve_fit first
    try:
        from scipy.optimize import curve_fit

        if fix_centroid and fix_fwhm:
            # Linear least-squares for H, a0, a1 with c, sigma fixed
            c_fixed = float(peak_channel)
            s_fixed = sigma_init
            G = np.exp(-(roi_x - c_fixed)**2 / (2 * s_fixed**2))
            A = np.column_stack([G, np.ones_like(roi_x), roi_x])
            # Weights: 1/sqrt(max(y, 1)) — Poisson statistics
            w = 1.0 / np.sqrt(np.maximum(roi_y, 1.0))
            Aw = A * w[:, None]
            yw = roi_y * w
            sol, residuals, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
            H_fit, a0_fit, a1_fit = sol
            c_fit = c_fixed
            sigma_fit = s_fixed
            converged = True
            # Approximate parameter uncertainties from residuals
            n_params = 3
            dof = max(1, len(roi_y) - n_params)
            chi2 = float(np.sum(((np.dot(A, sol) - roi_y) * w)**2))
            chi2_per_dof = chi2 / dof
            # H uncertainty from inverse of normal matrix
            try:
                cov = np.linalg.inv(Aw.T @ Aw) * chi2_per_dof
                H_unc = math.sqrt(max(0, cov[0, 0]))
            except np.linalg.LinAlgError:
                H_unc = H_fit * 0.05  # fallback 5%
        else:
            # Full nonlinear fit
            p0 = [H_init, float(peak_channel), sigma_init, a0_init, a1_init]
            bounds_lo = [0, peak_channel - half_window/2, sigma_init * 0.3,
                         -np.inf, -np.inf]
            bounds_hi = [np.inf, peak_channel + half_window/2,
                         sigma_init * 3.0, np.inf, np.inf]
            sigma_y = np.sqrt(np.maximum(roi_y, 1.0))  # Poisson errors
            popt, pcov = curve_fit(
                _gaussian_plus_baseline, roi_x, roi_y,
                p0=p0, sigma=sigma_y, absolute_sigma=True,
                bounds=(bounds_lo, bounds_hi), maxfev=2000,
            )
            H_fit, c_fit, sigma_fit, a0_fit, a1_fit = popt
            converged = True
            # Uncertainty in area: dA/dH · dH + dA/dσ · dσ (covariance terms ignored)
            H_unc = math.sqrt(max(0, pcov[0, 0]))
            chi2_per_dof = 1.0  # not directly available without manual recompute
            sigma_unc = math.sqrt(max(0, pcov[2, 2]))
            # Area = H·σ·√(2π); dA² = (σ·√(2π))²·dH² + (H·√(2π))²·dσ²
            # but ignore covariance for simplicity here
    except ImportError:
        # scipy fallback: use moments-based estimate
        # H = max - baseline at peak
        H_fit = H_init
        c_fit = float(peak_channel)
        sigma_fit = sigma_init
        a0_fit = a0_init
        a1_fit = a1_init
        H_unc = math.sqrt(max(1.0, H_fit))
        chi2_per_dof = float("nan")
        converged = False
    except Exception as e:
        return PeakAreaResult(
            method="gaussian_fit",
            peak_channel=float(peak_channel),
            net_area_counts=0.0,
            net_area_uncertainty=0.0,
            gross_counts=float(roi_y.sum()),
            baseline_counts=0.0,
            roi_low_ch=roi_low,
            roi_high_ch=roi_high,
            converged=False,
            notes=f"Gaussian fit failed: {e}",
        )

    # Area = H · σ · √(2π)
    sqrt_2pi = math.sqrt(2.0 * math.pi)
    area = H_fit * sigma_fit * sqrt_2pi
    # Uncertainty (Poisson-dominated for typical spectra):
    # dArea² ≈ (σ·√(2π))² · dH²  (assuming dσ small)
    area_unc = sigma_fit * sqrt_2pi * H_unc

    # Compute baseline contribution under peak
    baseline_counts = float(np.sum(a0_fit + a1_fit * roi_x))
    baseline_counts = max(0.0, baseline_counts)
    gross_counts = float(roi_y.sum())

    return PeakAreaResult(
        method="gaussian_fit",
        peak_channel=c_fit,
        net_area_counts=area,
        net_area_uncertainty=area_unc,
        gross_counts=gross_counts,
        baseline_counts=baseline_counts,
        roi_low_ch=roi_low,
        roi_high_ch=roi_high,
        fwhm_channels=sigma_fit * 2.355,
        fit_height=H_fit,
        fit_residual_norm=chi2_per_dof,
        converged=converged,
    )


def integrate_peaks(
    counts: np.ndarray,
    peaks: list,
    *,
    method: str = "cowell",
    fwhm_at_channel=None,
    **kwargs,
) -> list:
    """
    Convenience: compute peak areas for a list of FoundPeak objects.

    Args:
        counts: spectrum array
        peaks: list of FoundPeak objects (must have .channel, .fwhm_channels)
        method: "cowell", "gaussian_fit", or "gauss_erfc_step" (Gilmore §9.7)
        fwhm_at_channel: optional callable to override FWHM per peak
            (e.g. when peak.fwhm_channels is unreliable for low-σ peaks)
        **kwargs: passed to the integration function

    Returns:
        list of PeakAreaResult, one per input peak (in same order).
    """
    results = []
    for p in peaks:
        if fwhm_at_channel is not None:
            fwhm_ch = float(fwhm_at_channel(p.channel))
        else:
            fwhm_ch = float(getattr(p, "fwhm_channels", 0))
        if fwhm_ch <= 0:
            fwhm_ch = 5.0  # last-resort default
        try:
            if method == "cowell":
                r = cowell_area(
                    counts, peak_channel=int(p.channel),
                    fwhm_channels=fwhm_ch, **kwargs,
                )
            elif method == "gaussian_fit":
                r = gaussian_fit_area(
                    counts, peak_channel=int(p.channel),
                    fwhm_channels=fwhm_ch, **kwargs,
                )
            elif method == "gauss_erfc_step":
                from .area_step_continuum import gauss_erfc_step_fit
                r = gauss_erfc_step_fit(
                    counts, peak_channel=int(p.channel),
                    fwhm_channels=fwhm_ch, **kwargs,
                )
            else:
                raise ValueError(f"Unknown method: {method!r}")
        except Exception as e:
            r = PeakAreaResult(
                method=method,
                peak_channel=float(p.channel),
                net_area_counts=0.0,
                net_area_uncertainty=0.0,
                gross_counts=0.0,
                baseline_counts=0.0,
                roi_low_ch=p.channel,
                roi_high_ch=p.channel,
                converged=False,
                notes=f"Integration failed: {e}",
            )
        results.append(r)
    return results


# ============================================================================
# Lsrm built-in PEAKS table preferred — Cowell fallback
# ============================================================================

def get_peak_area(
    spec,
    peak_channel: int,
    fwhm_channels: float,
    *,
    prefer_lsrm_table: bool = True,
    match_tolerance_fwhm: float = 0.8,
    window_factor: float = 2.5,
):
    """
    Return `(area, uncertainty, source)` for one peak.

    Strategy:
      1. If `prefer_lsrm_table` and `spec.extras['lsrm_peaks_table']` is
         populated AND a row matches `peak_channel` within
         `match_tolerance_fwhm × FWHM`, return that Lsrm-fitted area.
         This is the **Lsrm software's own Gaussian-fit area**, written
         by Lsrm into the SPE file's `<START PEAKS>` table. It is
         typically much more accurate than Cowell on closely-spaced
         peaks (e.g. Co-60 1173/1332, Y-88 898/1836) because it uses
         a full Gaussian-on-step-baseline fit while Cowell uses a
         linear baseline that loses the wing counts to the
         neighbouring peak.
      2. Otherwise, fall back to Cowell integration.
      3. If Cowell fails (`converged=False`), return `(None, None, "failed")`.

    Args:
        spec: Spectrum object (must have `.counts` and optionally
            `.extras['lsrm_peaks_table']`).
        peak_channel: integer channel of the peak centroid.
        fwhm_channels: FWHM (in channels) at this peak's energy.
        prefer_lsrm_table: if False, skip Lsrm lookup and use Cowell.
        match_tolerance_fwhm: how close (in units of FWHM) an Lsrm-table
            entry must be to match this peak. Default 0.8 FWHM means
            "the Lsrm peak and our peak are essentially the same".
        window_factor: passed through to Cowell when falling back.

    Returns:
        Tuple of `(area, uncertainty, source)` where `source` ∈
        `{"lsrm_peaks_table", "cowell", "failed"}`.

    Background (F-31, v1.7.9):
        Discovered during F-31 cascade-summing validation work. On
        Co-60 spectra, Cowell with default window underestimates the
        1173/1332 photopeak areas by ~30% because the linear baseline
        runs through the neighbouring peak's wing rather than under
        the actual Compton continuum. Using the Lsrm-fitted areas
        instead restores cert-validated activity for Co-60 to within
        combined 1σ.
    """
    extras = getattr(spec, "extras", {}) or {}
    if prefer_lsrm_table:
        lsrm_tbl = extras.get("lsrm_peaks_table")
        if lsrm_tbl:
            tol = float(match_tolerance_fwhm) * float(fwhm_channels)
            best = None
            best_dist = float("inf")
            for entry in lsrm_tbl:
                ch_e = entry.get("position_ch")
                if ch_e is None:
                    continue
                d = abs(float(ch_e) - float(peak_channel))
                if d < tol and d < best_dist:
                    best_dist = d
                    best = entry
            if best is not None:
                area = best.get("area")
                d_area = best.get("d_area")
                if area is not None and area > 0:
                    return (float(area), float(d_area or 0.0),
                            "lsrm_peaks_table")

    # Cowell fallback
    try:
        r = cowell_area(
            spec.counts,
            peak_channel=int(peak_channel),
            fwhm_channels=float(fwhm_channels),
            window_factor=window_factor,
        )
    except Exception:
        return (None, None, "failed")
    if r.converged and r.net_area_counts > 0:
        return (r.net_area_counts, r.net_area_uncertainty, "cowell")
    return (None, None, "failed")


__all__ = [
    "PeakAreaResult",
    "cowell_area",
    "gaussian_fit_area",
    "integrate_peaks",
    "get_peak_area",
    "gauss_erfc_step_fit",
]


# Re-export so callers can import gauss_erfc_step_fit from scripts.gamma.peaks.area
from .area_step_continuum import gauss_erfc_step_fit  # noqa: E402
