"""
Effective dead-time correction (LSRM Algorithmic Foundations §15).

Why
---
The MCA-reported live time `t_live` already corrects for the
analyzer's basic dead time, but at high count rates (≥ 5·10⁴ cps for
typical Aspect electronics, or any rate where the σ-curve of an
isolated calibration peak starts to broaden against expected
calibration) it leaves residual losses from random coincidence and
pulse pile-up. LSRM §15 prescribes an additional **effective dead time**

    t_d  =  A · Σᵢ yᵢ  +  B · Σᵢ (yᵢ · i)     (LSRM Formula 15-1)

where yᵢ is the count in channel i and (A, B) are detector-specific
empirical coefficients. The corrected live time is

    t_live^corr  =  t_live − t_d.

For Gamma-1S (NaI 63×63 + Aspect spectrometer) the typical operating
range is below 5·10⁴ cps; for routine Marinelli / Дента / Петри
geometries dead-time is below 1 % and the correction is small.
Becomes critical for **5 cm point geometry** with strong sources
(Co-60 calibration, Cs-137 cert sources).

The (A, B) coefficients depend on the analyzer's input-stage
response (shaping time, pile-up rejection mode). LSRM §15 specifies a
three-source calibration protocol:

  1. Low-load source (Co-60 alone) at a distance giving ≤ 500 cps.
     Record area S₁ of the 1173 keV line; this is the reference
     count rate y_ref = S₁ / t_live.
  2. Add a second source (¹³³Ba or ²⁴¹Am low-energy emitter) without
     moving the Co-60; raise total load to ~5·10⁴ cps. Measure the
     same 1173 keV peak again → y_a.
  3. Replace the second source with ¹⁵²Eu or ¹³⁷Cs at similar high
     load → y_b.

The dead-time relative losses are

    δ_a = (y_ref − y_a) / y_ref
    δ_b = (y_ref − y_b) / y_ref

and (A, B) satisfy the linear system

    A · S_a + B · M_a  =  δ_a · t_live_a
    A · S_b + B · M_b  =  δ_b · t_live_b

where  Sᵢ = Σ y    (total counts)  and  Mᵢ = Σ (y · channel).

Solve the 2×2 system; store (A, B) in the detector profile.

Per-spectrum correction is then a one-line application of Formula
15-1 to get a corrected live time, which downstream activity formulas
use in place of the raw MCA live time.

Reference
---------
LSRM Algorithmic Foundations 2022, §15.1, стр. 15-1.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np


@dataclass(frozen=True)
class DeadTimeCoefficients:
    """A and B coefficients for the LSRM dead-time formula."""
    A_s_per_count: float          # seconds per count
    B_s_per_count_channel: float  # seconds per (count × channel)
    notes: str = ""
    calibration_sources: tuple = ()  # diagnostic — calibrating runs

    def __repr__(self) -> str:
        return (f"DeadTimeCoeff(A={self.A_s_per_count:.3e} s/count, "
                f"B={self.B_s_per_count_channel:.3e} s/(count·ch))")


@dataclass(frozen=True)
class DeadTimeResult:
    """Outcome of effective_dead_time()."""
    t_dead_s: float                    # the effective dead time
    t_live_corr_s: float               # t_live − t_dead
    fraction_lost: float               # t_dead / t_live (≥ 0)
    Sigma_y: float                     # Σᵢ yᵢ
    Sigma_yi: float                    # Σᵢ (yᵢ · i)
    applied: bool = True               # False if A, B not available
    notes: str = ""

    def __repr__(self) -> str:
        return (f"DeadTimeResult(t_d={self.t_dead_s:.4f} s, "
                f"loss={self.fraction_lost*100:.2f}%, applied={self.applied})")


def effective_dead_time(
    counts: Sequence[float],
    t_live_s: float,
    coeffs: Optional[DeadTimeCoefficients],
) -> DeadTimeResult:
    """
    Compute t_d = A·Σy + B·Σ(y·i) and the corrected live time.

    If `coeffs` is None (detector not calibrated for A, B), returns a
    DeadTimeResult with `applied=False` and `t_dead_s=0`, leaving the
    live time unchanged but surfacing the limitation in `notes`.
    """
    y = np.asarray(list(counts), dtype=float)
    Sigma_y = float(np.sum(y))
    # channel indices are 0-based — channel 0 contributes 0 to Σ(y·i)
    i = np.arange(len(y), dtype=float)
    Sigma_yi = float(np.sum(y * i))

    if coeffs is None:
        return DeadTimeResult(
            t_dead_s=0.0,
            t_live_corr_s=float(t_live_s),
            fraction_lost=0.0,
            Sigma_y=Sigma_y,
            Sigma_yi=Sigma_yi,
            applied=False,
            notes="A, B coefficients not calibrated for this detector — "
                  "LSRM §15 correction skipped; activities are "
                  "UNCORRECTED for residual pile-up / random coincidence",
        )

    t_d = coeffs.A_s_per_count * Sigma_y + coeffs.B_s_per_count_channel * Sigma_yi
    if t_d < 0:
        # Numerical artefact: clamp to zero
        t_d = 0.0
    t_live_corr = max(0.0, float(t_live_s) - t_d)
    frac = (t_d / float(t_live_s)) if t_live_s > 0 else 0.0
    note = ""
    if frac > 0.05:
        note = (f"dead-time loss = {frac*100:.1f}% (>5%) — correction "
                f"applied per LSRM §15")
    elif frac > 0.0:
        note = f"dead-time loss = {frac*100:.2f}% (negligible)"
    return DeadTimeResult(
        t_dead_s=t_d,
        t_live_corr_s=t_live_corr,
        fraction_lost=frac,
        Sigma_y=Sigma_y,
        Sigma_yi=Sigma_yi,
        applied=True,
        notes=note,
    )


# ---------------------------------------------------------------------------
# Calibration: solve 2×2 for (A, B) from three measurements
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class _CalibPoint:
    """One calibration spectrum's summary."""
    label: str
    Sigma_y: float
    Sigma_yi: float
    rate_at_ref_keV: float   # observed count rate at the chosen reference line
    t_live_s: float


def calibrate_AB(
    low_load: _CalibPoint,
    high_load_a: _CalibPoint,
    high_load_b: _CalibPoint,
) -> DeadTimeCoefficients:
    """
    Solve LSRM §15 calibration for (A, B).

    Per LSRM §15:
        y_ref = rate_at_ref_keV(low_load)       (negligible dead time)
        δ_a   = (y_ref − rate_at_ref_keV(high_a)) / y_ref
        δ_b   = (y_ref − rate_at_ref_keV(high_b)) / y_ref

    The two unknowns (A, B) satisfy:
        A · S_a + B · M_a = δ_a · t_live_a
        A · S_b + B · M_b = δ_b · t_live_b

    where S = Σy, M = Σ(y·i). Solve via 2×2 Cramer's rule.

    Parameters
    ----------
    low_load    : ≤ 500 cps reference spectrum (e.g. Co-60 alone)
    high_load_a : high-load combination with low-energy emitter
    high_load_b : high-load combination with high-energy emitter
    """
    y_ref = low_load.rate_at_ref_keV
    if y_ref <= 0:
        raise ValueError("low-load reference rate must be > 0")

    delta_a = (y_ref - high_load_a.rate_at_ref_keV) / y_ref
    delta_b = (y_ref - high_load_b.rate_at_ref_keV) / y_ref

    Sa, Ma = high_load_a.Sigma_y, high_load_a.Sigma_yi
    Sb, Mb = high_load_b.Sigma_y, high_load_b.Sigma_yi

    rhs_a = delta_a * high_load_a.t_live_s
    rhs_b = delta_b * high_load_b.t_live_s

    det = Sa * Mb - Sb * Ma
    if abs(det) < 1e-30:
        raise ValueError(
            "calibration system singular — high-load points are not "
            "independent enough (Σy and Σ(y·i) collinear)"
        )

    A = (rhs_a * Mb - rhs_b * Ma) / det
    B = (Sa * rhs_b - Sb * rhs_a) / det
    return DeadTimeCoefficients(
        A_s_per_count=float(A),
        B_s_per_count_channel=float(B),
        notes=(f"calibrated from 3 sources; δ_a={delta_a:.3f}, "
               f"δ_b={delta_b:.3f}"),
        calibration_sources=(low_load.label, high_load_a.label, high_load_b.label),
    )


def make_calib_point(
    label: str, counts: Sequence[float], t_live_s: float,
    rate_at_ref_keV: float,
) -> _CalibPoint:
    """Helper to build a _CalibPoint from raw counts + live time + ref rate."""
    y = np.asarray(list(counts), dtype=float)
    return _CalibPoint(
        label=label,
        Sigma_y=float(np.sum(y)),
        Sigma_yi=float(np.sum(y * np.arange(len(y), dtype=float))),
        rate_at_ref_keV=float(rate_at_ref_keV),
        t_live_s=float(t_live_s),
    )


# ---------------------------------------------------------------------------
# Gamma-1S default profile (uncalibrated stub)
# ---------------------------------------------------------------------------

GAMMA_1C_DEFAULT_COEFFS: Optional[DeadTimeCoefficients] = None
"""Set to None until a calibration session yields (A, B). When None,
`effective_dead_time()` returns applied=False and emits a notes warning."""


def get_detector_coeffs(detector_name: str) -> Optional[DeadTimeCoefficients]:
    """Lookup table for per-detector (A, B) coefficients."""
    if detector_name == "Gamma-1S":
        return GAMMA_1C_DEFAULT_COEFFS
    return None


__all__ = [
    "DeadTimeCoefficients",
    "DeadTimeResult",
    "effective_dead_time",
    "calibrate_AB",
    "make_calib_point",
    "get_detector_coeffs",
    "GAMMA_1C_DEFAULT_COEFFS",
]
