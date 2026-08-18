"""
Polynomial energy-calibration fit.

Implements the SKILL.md §5.3 rule:
  - Polynomial degree ≤ 4
  - Start at degree 1, raise only if residuals exceed 0.3·FWHM at any
    anchor point
  - When degree-4 still fails, the caller (bootstrap) is expected to
    segment the energy range and refit piecewise — that segmentation
    logic is in `bootstrap.py`

The fit itself is np.polyfit, returning low-to-high coefficients to
match the Spectrum.energy_cal convention.

Token economy: returns a small dict of coefficients + residuals.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np


# Hard cap from SKILL.md scope
MAX_POLYNOMIAL_DEGREE = 4


@dataclass
class EnergyFitResult:
    """Outcome of a single polynomial energy fit."""
    coefficients: tuple              # low-to-high: a0, a1, ...
    degree: int
    n_points: int
    residuals_keV: list = field(default_factory=list)  # per-anchor: observed - predicted
    max_residual_keV: float = 0.0
    rms_residual_keV: float = 0.0
    converged: bool = True
    reason: str = ""

    def predict(self, channels) -> np.ndarray:
        """Evaluate the polynomial at given channels."""
        ch = np.asarray(channels, dtype=np.float64)
        out = np.zeros_like(ch)
        for c in reversed(self.coefficients):
            out = out * ch + c
        return out


def polynomial_energy_fit(
    channels,
    energies,
    *,
    max_degree: int = MAX_POLYNOMIAL_DEGREE,
    target_residual_keV: Optional[float] = None,
    min_degree: int = 1,
) -> EnergyFitResult:
    """
    Fit E(N) = a0 + a1·N + a2·N² + ... with adaptive degree.

    Starts at min_degree, raises by 1 until either max_residual_keV
    falls below target_residual_keV OR degree reaches max_degree.

    Args:
        channels: list of peak channel positions (integers or floats)
        energies: list of expected energies in keV, same length as channels
        max_degree: ceiling on polynomial degree (default 4)
        target_residual_keV: stop raising degree when max_residual <= this.
                             If None, always fits at max_degree without
                             early stopping (useful for tests).
        min_degree: starting degree (default 1)

    Returns:
        EnergyFitResult with the lowest-degree fit that satisfies the
        target (or, if none does, the max_degree fit).
    """
    channels = np.asarray(channels, dtype=np.float64)
    energies = np.asarray(energies, dtype=np.float64)
    n = channels.size
    if n != energies.size:
        raise ValueError("channels and energies must have the same length")
    if n < 2:
        return EnergyFitResult(
            coefficients=(),
            degree=0,
            n_points=n,
            converged=False,
            reason=f"Need at least 2 points; got {n}",
        )

    # Cap degree at n-1 (no over-determined fit) and at max_degree
    effective_cap = min(max_degree, n - 1)
    if effective_cap < min_degree:
        # Force min_degree even if it means underdetermined; fall back
        # to whatever is achievable
        min_degree = max(0, effective_cap)

    best_result = None
    for d in range(min_degree, effective_cap + 1):
        # np.polyfit returns coefficients high-to-low; reverse
        coefs_high_to_low = np.polyfit(channels, energies, d)
        coefs = tuple(float(c) for c in coefs_high_to_low[::-1])

        # Evaluate at the fit points
        predicted = np.zeros_like(channels)
        for c in reversed(coefs):
            predicted = predicted * channels + c
        residuals = energies - predicted
        max_res = float(np.max(np.abs(residuals)))
        rms_res = float(np.sqrt(np.mean(residuals ** 2)))

        result = EnergyFitResult(
            coefficients=coefs,
            degree=d,
            n_points=n,
            residuals_keV=residuals.tolist(),
            max_residual_keV=max_res,
            rms_residual_keV=rms_res,
            converged=True,
            reason=f"Degree {d} fit, max residual {max_res:.3f} keV",
        )
        best_result = result

        if target_residual_keV is not None and max_res <= target_residual_keV:
            best_result.reason += f" — meets target {target_residual_keV:.3f} keV"
            return best_result

    if target_residual_keV is not None:
        best_result.reason += (
            f" — did NOT meet target {target_residual_keV:.3f} keV "
            f"at max degree {best_result.degree}; piecewise fit may be needed"
        )

    return best_result
