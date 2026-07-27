"""
FWHM(E) calibration.

Two models are supported:

  HPGe (and other semiconductors):
    FWHM²(E) = a + b·E + c·E²
    (Gilmore §6.2, three-component: electronic noise, charge collection,
    Fano factor)

  Scintillator (NaI, LaBr₃, CeBr₃, CdZnTe):
    FWHM(E) = k · √(E + α·E²)
    where α captures non-proportionality of light yield

The choice between models is informed by the detector_type step (which
in turn uses the relative resolution at e.g. 662 keV). Here we provide
both, and the caller picks the appropriate one.

Token economy: returns a small dataclass with coefficients and residuals.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np


# BUG-37 / Wave 7 / 2026-06-05 — numerical floor for FWHM^2(E) inside
# ``FwhmFitResult.fwhm_at``. The legacy clamp ``np.maximum(val, 0.0)``
# would return ``sqrt(0) = 0`` for any energy where the fitted quadratic
# went negative, which downstream produced a 0-width match window
# (anything-matches / nothing-matches depending on the consumer).
# Mirrors the 0.01 keV^2 floor used in
# ``scripts/gamma/identification/staged_pipeline.py:fwhm_keV_at_energy``
# (FWHM floor ~0.1 keV). Pure last-resort safety net for
# numerically-degenerate fits; the **real** BUG-41 fix is in
# ``staged_pipeline.build_fwhm_model`` (pathology detection + LSRM
# stored sqrt(E) polynomial fallback + default-NaI-63x63 backstop).
_FWHM2_FLOOR_keV2 = 0.01


@dataclass
class FwhmFitResult:
    """Outcome of a FWHM(E) fit."""
    model: str                                    # "hpge" | "scintillator"
    coefficients: tuple                           # model-dependent
    n_points: int
    residuals_keV: list = field(default_factory=list)
    max_residual_keV: float = 0.0
    rms_residual_keV: float = 0.0
    converged: bool = True
    reason: str = ""

    def fwhm_at(self, E_keV) -> float:
        """Evaluate FWHM (keV) at energy E.

        BUG-37 fix: clamp ``val`` (i.e. ``FWHM^2``) at
        ``_FWHM2_FLOOR_keV2`` (= 0.01 keV^2 -> FWHM = 0.1 keV) instead of
        ``0.0``. This avoids ``sqrt(0)`` -> 0 keV -> 0-width match
        windows when an upstream HPGe quadratic fit went pathological at
        low E (val < 0 region). The real fix for pathological quadratics
        lives in ``staged_pipeline.build_fwhm_model``; this floor is a
        last-resort numerical safety net only.
        """
        E = np.asarray(E_keV, dtype=np.float64)
        if self.model == "hpge":
            a, b, c = self.coefficients
            val = a + b * E + c * E * E
            return np.sqrt(np.maximum(val, _FWHM2_FLOOR_keV2))
        elif self.model == "scintillator":
            k, alpha = self.coefficients
            val = E + alpha * E * E
            return k * np.sqrt(np.maximum(val, _FWHM2_FLOOR_keV2))
        raise ValueError(f"Unknown FWHM model: {self.model}")


def fit_fwhm_hpge(energies, fwhms_keV) -> FwhmFitResult:
    """
    Fit FWHM²(E) = a + b·E + c·E² (HPGe model).

    Args:
        energies: list of calibration peak energies in keV
        fwhms_keV: list of measured FWHM at those energies, in keV

    Returns:
        FwhmFitResult with model="hpge".
    """
    energies = np.asarray(energies, dtype=np.float64)
    fwhms = np.asarray(fwhms_keV, dtype=np.float64)
    n = energies.size

    if n < 3:
        # With < 3 points we drop the quadratic term
        if n < 2:
            return FwhmFitResult(
                model="hpge",
                coefficients=(0.0, 0.0, 0.0),
                n_points=n,
                converged=False,
                reason=f"Need ≥2 points; got {n}",
            )
        # Linear in E: FWHM² = a + b·E
        A = np.vstack([np.ones_like(energies), energies]).T
        sol, *_ = np.linalg.lstsq(A, fwhms ** 2, rcond=None)
        coefs = (float(sol[0]), float(sol[1]), 0.0)
    else:
        # Quadratic in E for FWHM²
        A = np.vstack([np.ones_like(energies), energies, energies ** 2]).T
        sol, *_ = np.linalg.lstsq(A, fwhms ** 2, rcond=None)
        coefs = (float(sol[0]), float(sol[1]), float(sol[2]))

    a, b, c = coefs
    predicted_sq = a + b * energies + c * energies ** 2
    predicted = np.sqrt(np.maximum(predicted_sq, 0.0))
    residuals = fwhms - predicted
    max_res = float(np.max(np.abs(residuals)))
    rms_res = float(np.sqrt(np.mean(residuals ** 2)))

    return FwhmFitResult(
        model="hpge",
        coefficients=coefs,
        n_points=n,
        residuals_keV=residuals.tolist(),
        max_residual_keV=max_res,
        rms_residual_keV=rms_res,
        converged=True,
        reason=f"HPGe FWHM² model fit, max residual {max_res:.3f} keV",
    )


def fit_fwhm_scintillator(energies, fwhms_keV) -> FwhmFitResult:
    """
    Fit FWHM(E) = k · √(E + α·E²) (scintillator model).

    Args:
        energies: list of calibration peak energies in keV
        fwhms_keV: measured FWHM at those energies, in keV

    Returns:
        FwhmFitResult with model="scintillator".
    """
    energies = np.asarray(energies, dtype=np.float64)
    fwhms = np.asarray(fwhms_keV, dtype=np.float64)
    n = energies.size

    if n < 2:
        return FwhmFitResult(
            model="scintillator",
            coefficients=(0.0, 0.0),
            n_points=n,
            converged=False,
            reason=f"Need ≥2 points; got {n}",
        )

    # Linearise: FWHM² / E = k² · (1 + α·E)
    # So if we let Y = FWHM²/E and X = E, then Y = k² + k²·α·X
    # Linear in (k², k²·α). Then α = (k²·α) / k².
    mask = energies > 0
    Y = (fwhms[mask] ** 2) / energies[mask]
    X = energies[mask]
    if X.size < 2:
        # All energies zero or invalid; fall back to constant k
        k = float(np.mean(fwhms / np.sqrt(np.maximum(energies, 1e-6))))
        coefs = (k, 0.0)
    else:
        A = np.vstack([np.ones_like(X), X]).T
        sol, *_ = np.linalg.lstsq(A, Y, rcond=None)
        k_sq = float(sol[0])
        k_sq_alpha = float(sol[1])
        if k_sq <= 0:
            # Degenerate; use direct estimate
            k = float(np.mean(fwhms / np.sqrt(np.maximum(X, 1e-6))))
            alpha = 0.0
        else:
            k = math.sqrt(k_sq)
            alpha = k_sq_alpha / k_sq if k_sq > 0 else 0.0
        coefs = (k, alpha)

    k, alpha = coefs
    predicted = k * np.sqrt(np.maximum(energies + alpha * energies ** 2, 0.0))
    residuals = fwhms - predicted
    max_res = float(np.max(np.abs(residuals)))
    rms_res = float(np.sqrt(np.mean(residuals ** 2)))

    return FwhmFitResult(
        model="scintillator",
        coefficients=coefs,
        n_points=n,
        residuals_keV=residuals.tolist(),
        max_residual_keV=max_res,
        rms_residual_keV=rms_res,
        converged=True,
        reason=f"Scintillator FWHM model fit, max residual {max_res:.3f} keV",
    )
