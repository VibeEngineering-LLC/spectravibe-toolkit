"""
Photopeak efficiency calibration ε(E) per Lsrm §8.5.

Efficiency is the ratio of counts detected in the full-energy peak to
the number of γ-rays emitted by the source:

    ε(E) = S(E) / (A · I(E) · t)

where:
    S(E)  — net photopeak area at energy E
    A     — source activity (Bq)
    I(E)  — gamma-line emission probability (decimal)
    t     — live time (s)

The fitted curve ε(E) is the central object of Phase 2.1c:
  • Identification: proportionality check uses (area/ε/I) ratios
    instead of bare area ratios — closes K-12.
  • MDA: ε(E) is the denominator in detection-limit formula
    (closes K-06).
  • Activity calculation (Phase 2.1d): A = S / (ε · I · t).

Model selection per Lsrm §8.5.1 and §8.5.2:

  • NaI / scintillators (50-3000 keV):
        log ε(E) = a₀ + a₁·log E + a₂·(log E)² + ... + aₙ·(log E)ⁿ
    Polynomial in log-log space — captures both the low-E rise
    (geometric ε increasing through threshold) and the high-E fall
    (decreasing photoelectric absorption). Typical n=3 or 4.

  • HPGe in similar range: usually a piecewise model with break at
    ~120 keV where K-edge transitions stop dominating.

  • Volumetric samples (Marinelli, Petri): same log-log polynomial
    but coefficients differ from point-source by ~factor 5-10 at
    low E (geometric area effect) and by smaller factor at high E.

Source data: ε(E) calibration points come from .efa (aggregated) or
.efr (per-source) files parsed by `gamma.io.lsrm_efficiency`. Each
calibration point has known E, ε(E), and uncertainty dε%.

The fit uses **weighted log-log polynomial regression**:
  • Convert all data points to (log E, log ε) space
  • Convert uncertainty: d(log ε) = dε/ε ≈ dε_pct/100 (for small %)
  • Weighted least squares with weights w_i = 1/d(log ε_i)²
  • Default polynomial degree: 3 (sufficient for NaI 50-3000 keV)
  • Report chi²/dof to flag underfit/overfit

Reference: Lsrm Algorithmic Foundations 2022, §8.5;
Knoll, Radiation Detection and Measurement, 4th Ed., §10.III.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Callable, List

import numpy as np


# Reference energy at which efficiency uncertainty is reported in the
# fit output (typically a well-measured line, e.g. Cs-137 661 keV)
REFERENCE_ENERGY_KEV = 661.66


@dataclass(frozen=True)
class EfficiencyCurve:
    """
    Fitted photopeak efficiency curve ε(E).

    The curve is a polynomial in log-log space:
        log ε(E) = Σ_k coefs[k] · (log E)^k

    Args:
        coefficients: polynomial coefficients in log-log space, low to
            high order
        E_min_keV, E_max_keV: validity range (outside which extrapolation
            should be flagged)
        chi2_per_dof: goodness-of-fit metric (≈1 ideal)
        n_points_used: number of data points used in fit
        n_dof: degrees of freedom (n_points - n_coefs)
        residuals_pct: array of (E, residual/ε in %) for diagnostic
        detector_id: detector identifier (for traceability)
        geometry: source geometry (for traceability)
        notes: free-text diagnostic notes
    """
    coefficients: tuple
    E_min_keV: float
    E_max_keV: float
    chi2_per_dof: float
    n_points_used: int
    n_dof: int
    residuals_pct: tuple = ()
    detector_id: str = ""
    geometry: str = ""
    notes: str = ""

    def __call__(self, E_keV: float) -> float:
        """Evaluate ε(E) at a given energy."""
        return self.efficiency_at(E_keV)

    def efficiency_at(self, E_keV: float) -> float:
        """
        Compute ε(E) by evaluating the log-log polynomial.

        Args:
            E_keV: energy (keV) at which to evaluate efficiency

        Returns:
            ε(E) — photopeak efficiency (dimensionless, typical 1e-4 .. 0.1)

        Notes:
            • Returns 0 for E ≤ 0 (unphysical)
            • Extrapolation outside [E_min, E_max] is permitted but
              should be used with caution (the log-log polynomial may
              behave wildly far from fit support)
        """
        if E_keV <= 0:
            return 0.0
        log_E = math.log(E_keV)
        log_eps = 0.0
        for k, a_k in enumerate(self.coefficients):
            log_eps += a_k * (log_E ** k)
        return math.exp(log_eps)

    def is_extrapolating(self, E_keV: float, margin_factor: float = 1.1) -> bool:
        """True if E lies outside the calibrated range."""
        return E_keV < self.E_min_keV / margin_factor or \
               E_keV > self.E_max_keV * margin_factor

    def __repr__(self) -> str:
        return (f"EfficiencyCurve(degree={len(self.coefficients)-1}, "
                f"E=[{self.E_min_keV:.0f}-{self.E_max_keV:.0f}] keV, "
                f"χ²/dof={self.chi2_per_dof:.2f}, "
                f"n_points={self.n_points_used})")


def fit_efficiency_curve(
    energies_keV: list,
    efficiencies: list,
    uncertainties_pct: list,
    *,
    degree: int = 3,
    detector_id: str = "",
    geometry: str = "",
) -> EfficiencyCurve:
    """
    Fit a photopeak efficiency curve to calibration data points.

    Uses weighted least-squares regression of log ε on log E with a
    polynomial of specified degree.

    Args:
        energies_keV: list of γ-line energies (keV) where ε was measured
        efficiencies: list of ε values (decimal) — same length as energies
        uncertainties_pct: list of relative uncertainties (%) — same length
        degree: polynomial order in log-log space (default 3 for NaI)
        detector_id: detector identifier (passed through to result)
        geometry: source geometry name (passed through to result)

    Returns:
        EfficiencyCurve.

    Algorithm:
        1. Convert (E_i, ε_i) → (x_i, y_i) = (log E_i, log ε_i)
        2. Convert relative uncertainty d(log ε_i) = dε_i / ε_i ≈
           uncertainties_pct[i] / 100 (Taylor approximation for small
           relative errors)
        3. Solve weighted normal equations:
               (A^T · W · A) · coefs = A^T · W · y
           where A_ij = x_i^j, W_ii = 1/d(log ε_i)²
        4. Compute χ² = Σ w_i · (y_i - Σ_j coefs_j · x_i^j)²
        5. Return EfficiencyCurve with chi²/dof and residuals

    Notes:
        • degree=3 (4 coefficients) is appropriate for NaI 50-3000 keV
        • Higher degree (5+) overfits typical 15-30 calibration points
        • Lower degree (1-2) underfits the K-edge curvature on HPGe
        • If chi²/dof > 5, the model is suspect — check data outliers
        • If chi²/dof < 0.5, uncertainties may be overestimated
    """
    if len(energies_keV) != len(efficiencies) or \
       len(energies_keV) != len(uncertainties_pct):
        raise ValueError("energies, efficiencies, uncertainties must have same length")
    if len(energies_keV) < degree + 1:
        raise ValueError(f"Need at least {degree+1} points for degree-{degree} fit, "
                         f"got {len(energies_keV)}")
    # Filter out invalid points
    valid = [
        (E, eps, dpct)
        for E, eps, dpct in zip(energies_keV, efficiencies, uncertainties_pct)
        if E > 0 and eps > 0 and dpct > 0
    ]
    if len(valid) < degree + 1:
        raise ValueError(f"After filtering, only {len(valid)} valid points "
                         f"(need ≥{degree+1})")

    E_arr = np.array([v[0] for v in valid], dtype=np.float64)
    eps_arr = np.array([v[1] for v in valid], dtype=np.float64)
    dpct_arr = np.array([v[2] for v in valid], dtype=np.float64)

    # Transform to log-log space
    x = np.log(E_arr)
    y = np.log(eps_arr)
    # d(log ε) for small relative error ≈ dε/ε = dpct/100
    d_log_eps = dpct_arr / 100.0
    w = 1.0 / (d_log_eps ** 2)  # weights

    # Build Vandermonde-like matrix: A[i,j] = x[i]^j
    n_pts = len(x)
    n_coefs = degree + 1
    A = np.zeros((n_pts, n_coefs))
    for j in range(n_coefs):
        A[:, j] = x ** j

    # Weighted normal equations: (A^T W A) coefs = A^T W y
    W = np.diag(w)
    ATW = A.T @ W
    AtWA = ATW @ A
    AtWy = ATW @ y
    try:
        coefs = np.linalg.solve(AtWA, AtWy)
    except np.linalg.LinAlgError:
        # Fall back to unweighted polyfit
        coefs = np.polyfit(x, y, degree)[::-1]  # polyfit returns high-to-low

    # Compute residuals + chi²
    y_fit = A @ coefs
    residuals_log = y - y_fit  # in log space
    # Convert to % deviation in ε space: dε/ε ≈ d(log ε)
    residuals_pct = tuple(
        (float(E_arr[i]), float(residuals_log[i] * 100.0))
        for i in range(n_pts)
    )
    chi2 = float(np.sum(w * residuals_log ** 2))
    dof = max(1, n_pts - n_coefs)
    chi2_per_dof = chi2 / dof

    return EfficiencyCurve(
        coefficients=tuple(float(c) for c in coefs),
        E_min_keV=float(E_arr.min()),
        E_max_keV=float(E_arr.max()),
        chi2_per_dof=chi2_per_dof,
        n_points_used=n_pts,
        n_dof=dof,
        residuals_pct=residuals_pct,
        detector_id=detector_id,
        geometry=geometry,
        notes=(f"Log-log polynomial fit, degree {degree}, "
               f"χ²={chi2:.2f}, χ²/dof={chi2_per_dof:.3f}"),
    )


def fit_efficiency_from_efr_file(
    efr_path: str,
    *,
    degree: int = 3,
) -> EfficiencyCurve:
    """
    Convenience: load an .efr or .efa file and fit an efficiency curve.

    Args:
        efr_path: path to Lsrm .efa or .efr file
        degree: polynomial degree for log-log fit (default 3)

    Returns:
        EfficiencyCurve fitted to all points across all blocks in the
        file. Detector and geometry attributed from the first block.
    """
    from gamma.io.lsrm_efficiency import read_efficiency_file
    eff_file = read_efficiency_file(efr_path)
    if not eff_file.blocks:
        raise ValueError(f"No blocks in efficiency file: {efr_path}")
    all_pts = []
    for b in eff_file.blocks:
        all_pts.extend(b.points)
    if not all_pts:
        raise ValueError(f"No efficiency points in: {efr_path}")
    # Sort by energy
    all_pts.sort(key=lambda p: p.energy_keV)
    detector = eff_file.blocks[0].detector
    geometry = eff_file.blocks[0].geometry
    return fit_efficiency_curve(
        energies_keV=[p.energy_keV for p in all_pts],
        efficiencies=[p.efficiency for p in all_pts],
        uncertainties_pct=[p.efficiency_uncertainty_pct for p in all_pts],
        degree=degree,
        detector_id=detector,
        geometry=geometry,
    )


__all__ = [
    "REFERENCE_ENERGY_KEV",
    "EfficiencyCurve",
    "fit_efficiency_curve",
    "fit_efficiency_from_efr_file",
]
