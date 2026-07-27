"""
Subcalibration (Подкалибровка) — linear-only refit of an existing
energy calibration.

Per Lsrm Algorithmic Foundations §8.2.1, the nonlinear part of an
energy calibration is dominated by intrinsic detector properties (PMT
non-linearity, scintillator non-proportionality, ADC integral
non-linearity) and is stable over the operating life of the detector.
The linear part (a₀ = zero-point offset, a₁ = gain) is what drifts
under operating conditions — temperature, supply voltage variations,
PMT gain aging, mechanical shock during transport.

Subcalibration adjusts ONLY a₀ and a₁ while keeping a₂, a₃, ...
fixed from the stored polynomial. This is the appropriate response
when stored_check finds the matches systematically off-shifted from
their library positions but the residuals are *small relative to
FWHM* (i.e. the shape of the energy-channel curve is right, only the
overall placement is wrong).

When subcalibration applies:
  - stored_check returns matches with mean residual ≤ 1·FWHM (the
    calibration is close but not exact)
  - stored_check returns ≥ 2 matched anchors (we need two linear
    parameters to fit)
  - the residuals show a *linear* trend in channel (offset+slope),
    not a polynomial trend (which would indicate the nonlinear part
    has also changed)

When subcalibration does NOT apply (fall back to full bootstrap):
  - residuals exceed 1·FWHM at multiple anchors → nonlinear part has
    changed
  - too few matched anchors → cannot determine a₀, a₁ reliably
  - residuals show clear curvature → polynomial part is wrong, need
    full refit

Per Lsrm methodology, subcalibration is much faster and more
reliable than bootstrap when applicable, because:
  - Fewer free parameters → smaller statistical uncertainty
  - Preserves the carefully-calibrated nonlinear shape of the
    detector response
  - Cannot diverge from physical reality (linear-only)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


@dataclass(frozen=True)
class SubcalibrationResult:
    """Result of a subcalibration (linear-only) refit."""

    success: bool
    reason: str

    # Updated linear coefficients (a₀ and a₁ refit; a₂...aₙ from stored)
    coefficients: Tuple[float, ...] = ()
    degree: int = 0

    # Residuals at each used anchor, in keV
    residuals_keV: Tuple[float, ...] = ()
    max_residual_keV: float = 0.0
    mean_residual_keV: float = 0.0

    # Which anchors were used
    anchors_used: Tuple[Tuple[int, float, str], ...] = ()

    # Diagnostic: the shifts found
    a0_shift: float = 0.0  # Δa₀ from stored
    a1_relative_shift: float = 0.0  # (a₁_new - a₁_stored) / a₁_stored


def subcalibration_refit(
    spec,
    matched_anchors: list,
    *,
    max_acceptable_residual_keV: Optional[float] = None,
    fwhm_at_channel=None,
) -> SubcalibrationResult:
    """
    Linear-only refit of stored energy calibration.

    Algorithm:
      1. Take stored polynomial p(N) = a₀ + a₁·N + a₂·N² + ... + aₖ·Nᵏ
      2. Split as p(N) = (a₀ + a₁·N) + g(N), where g(N) = a₂·N² + ...
         is the "nonlinear residue" — kept fixed.
      3. Fit anchor data: E_i ≈ a₀' + a₁'·N_i + g(N_i)
         → solve linear system for (a₀', a₁') from observations
         (N_i, E_i - g(N_i))
      4. New calibration: (a₀', a₁', a₂, a₃, ...)
      5. Sanity-check residuals against max_acceptable threshold

    Args:
        spec: Spectrum with `energy_cal` (tuple of coefficients,
              lowest order first)
        matched_anchors: list of (channel, library_energy_keV,
                         source_label) tuples — typically from
                         stored_check
        max_acceptable_residual_keV: reject the refit if max residual
              exceeds this. If None, use 1·FWHM(highest_anchor) from
              `fwhm_at_channel` when supplied, else 5.0 keV.
        fwhm_at_channel: optional callable to compute FWHM(channel) for
              adaptive acceptance threshold. If None, falls back to
              5 keV.

    Returns:
        SubcalibrationResult.
    """
    if not spec.energy_cal or len(spec.energy_cal) < 2:
        return SubcalibrationResult(
            success=False,
            reason="Stored calibration must have ≥2 coefficients to "
                   "subcalibrate (need a₀ and a₁ to refit, plus "
                   "nonlinear part to preserve).",
        )

    if len(matched_anchors) < 2:
        return SubcalibrationResult(
            success=False,
            reason=f"Need ≥2 matched anchors to refit linear part; "
                   f"got {len(matched_anchors)}",
        )

    stored = list(spec.energy_cal)
    a0_stored = float(stored[0])
    a1_stored = float(stored[1])
    nonlin = stored[2:]  # a₂, a₃, ... — these stay fixed

    # Compute the nonlinear residue g(N) at each anchor channel
    channels = np.array([float(a[0]) for a in matched_anchors])
    library_energies = np.array([float(a[1]) for a in matched_anchors])

    def g_of_N(N: float) -> float:
        """Nonlinear part of the stored polynomial: a₂·N² + a₃·N³ + ..."""
        if not nonlin:
            return 0.0
        # Horner-like evaluation starting from highest degree
        result = 0.0
        for i, c in enumerate(nonlin, start=2):
            result += float(c) * (N ** i)
        return result

    g_values = np.array([g_of_N(N) for N in channels])

    # Target for linear fit: E_i - g(N_i) = a₀' + a₁'·N_i
    targets = library_energies - g_values

    # Solve linear least-squares [1, N] · [a₀', a₁']ᵀ = targets
    A = np.vstack([np.ones_like(channels), channels]).T
    sol, *_ = np.linalg.lstsq(A, targets, rcond=None)
    a0_new = float(sol[0])
    a1_new = float(sol[1])

    # Build new full polynomial
    new_coefs = tuple([a0_new, a1_new] + [float(c) for c in nonlin])

    # Compute residuals at each anchor with the new calibration
    def eval_poly(N: float) -> float:
        result = 0.0
        for i, c in enumerate(new_coefs):
            result += c * (N ** i)
        return result

    residuals = []
    for ch, E_lib, _src in matched_anchors:
        E_predicted = eval_poly(float(ch))
        residuals.append(abs(E_predicted - float(E_lib)))
    residuals_arr = np.array(residuals)
    max_resid = float(residuals_arr.max())
    mean_resid = float(residuals_arr.mean())

    # Acceptance check
    if max_acceptable_residual_keV is None:
        if fwhm_at_channel is not None:
            highest_ch = int(channels.max())
            try:
                fwhm_ch_at_high = float(fwhm_at_channel(highest_ch))
                max_acceptable_residual_keV = fwhm_ch_at_high * a1_new
            except Exception:
                max_acceptable_residual_keV = 5.0
        else:
            max_acceptable_residual_keV = 5.0

    if max_resid > max_acceptable_residual_keV:
        return SubcalibrationResult(
            success=False,
            reason=f"Max residual {max_resid:.2f} keV exceeds "
                   f"acceptance threshold {max_acceptable_residual_keV:.2f} keV "
                   f"— nonlinear part of stored calibration may have "
                   f"changed; recommend full bootstrap instead.",
            coefficients=new_coefs,
            degree=len(new_coefs) - 1,
            residuals_keV=tuple(float(r) for r in residuals_arr),
            max_residual_keV=max_resid,
            mean_residual_keV=mean_resid,
            anchors_used=tuple((int(a[0]), float(a[1]), str(a[2]))
                              for a in matched_anchors),
            a0_shift=a0_new - a0_stored,
            a1_relative_shift=(a1_new - a1_stored) / a1_stored if a1_stored else 0.0,
        )

    return SubcalibrationResult(
        success=True,
        reason=f"Subcalibration converged: {len(matched_anchors)} anchors, "
               f"max residual {max_resid:.2f} keV, mean {mean_resid:.2f} keV. "
               f"Δa₀={a0_new - a0_stored:+.3f} keV, "
               f"Δa₁/a₁={100*(a1_new - a1_stored)/a1_stored:+.2f}% "
               f"relative to stored.",
        coefficients=new_coefs,
        degree=len(new_coefs) - 1,
        residuals_keV=tuple(float(r) for r in residuals_arr),
        max_residual_keV=max_resid,
        mean_residual_keV=mean_resid,
        anchors_used=tuple((int(a[0]), float(a[1]), str(a[2]))
                          for a in matched_anchors),
        a0_shift=a0_new - a0_stored,
        a1_relative_shift=(a1_new - a1_stored) / a1_stored if a1_stored else 0.0,
    )


__all__ = ["SubcalibrationResult", "subcalibration_refit"]
