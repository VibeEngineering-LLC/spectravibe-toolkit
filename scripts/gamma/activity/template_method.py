"""
Template method — LSRM Algorithmic Foundations §12.

Concept
-------
The «шаблонный метод» (template method) is LSRM's preferred routine-
assay scheme for scintillator detectors (NaI, CeBr₃) measuring samples
of known nuclide composition (e.g. routine ¹³⁷Cs + ⁴⁰K + ²³²Th content
in environmental Marinelli). The detector + geometry is first
calibrated with one reference source per candidate nuclide, producing
a *sensitivity matrix* R of shape (M_windows × K_nuclides). Subsequent
sample measurements are then reduced to solving a linear system

    n_i  =  Σ_k  R_ki · A_k  +  B_i

where n_i is the count rate in the i-th energy window, B_i is the
corresponding bg window rate, and A_k is the specific activity of
nuclide k in the sample. The solution is found by weighted LSQ:

    χ²  =  Σ_i  w_i · (n_i / t_p  −  Σ_k R_ki · A_k  −  B_i / t_b)²

with weights w_i = 1 / σ²_total. The system matrix `R^T · W · R` is
solved by Cholesky (here we use the more robust `numpy.linalg.lstsq`).

Calibration formulas (LSRM §12.2):

    R_ki  =  (S_kr_i / t_kr  −  B_i / t_b)  /  A_k

where S_kr_i is the gross counts in window i of the k-th calibration
spectrum, t_kr its live time, and A_k the certified specific activity
of the calibration source.

For routine assay the activity of each candidate nuclide is solved
simultaneously across all windows. The uncertainty σ²(A_k) is the
k-th diagonal of (Rᵀ·W·R)⁻¹.

Why this is the right tool for Gamma-1S routine assay
-----------------------------------------------------
For the project's bread-and-butter Marinelli measurements (¹³⁷Cs +
⁴⁰K + ²³²Th + ²³⁸U + ²²⁶Ra), the template method:

* Builds R once per geometry from cert sources we ALREADY have
  (`detectors/Gamma-1S/data/cert_fixtures/` — Cs-137, K-40, Th-232,
  Ra-226 Marinelli, Дента, Петри). No need to find peaks for every
  routine sample.
* Reduces analysis time per spectrum to a single matrix solve —
  ~seconds vs. ~minute for the full 11-step pipeline.
* Naturally handles the ~7% NaI FWHM by working on broad energy
  windows instead of resolving overlapping FEPs.
* Activities and σ(A) come directly from the cert-source
  uncertainties — no efficiency curve refit per sample.

Scope of this module (v1.17.0)
------------------------------
**Fitter + matrix builder only** — additive. The CLI hook
`--routine-assay` and the calibration-batch utility
`build_sensitivity_matrix.py` are deferred to follow-up tasks.

Inputs:
    calibration_specs : dict {nuclide: CalibrationSpec(...)}
                        Each CalibrationSpec carries the cert
                        spectrum counts, t_live, and certified A_k.
    bg_spec           : BackgroundSpec(counts, t_live)
    windows_keV       : list of (E_low, E_high) tuples
    energies_per_ch   : per-channel energy

Then `build_R_matrix(...)` returns the sensitivity matrix; and
`template_assay(sample_counts, t_sample, R, bg_window_rates, windows)`
solves for A.

Reference
---------
LSRM Algorithmic Foundations 2022, §12 «Шаблонный метод», стр. 12-1.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Window-sum helper
# ---------------------------------------------------------------------------

def sum_in_windows(
    counts: Sequence[float],
    energies_per_ch: Sequence[float],
    windows_keV: Sequence[Tuple[float, float]],
) -> np.ndarray:
    """
    Return total counts in each energy window.

    Channels whose centre energy falls within [E_low, E_high) are
    included. Out-of-range bins contribute 0.
    """
    counts = np.asarray(list(counts), dtype=float)
    E = np.asarray(list(energies_per_ch), dtype=float)
    out = np.zeros(len(windows_keV), dtype=float)
    for i, (lo, hi) in enumerate(windows_keV):
        mask = (E >= lo) & (E < hi)
        out[i] = float(np.sum(counts[mask]))
    return out


# ---------------------------------------------------------------------------
# Calibration / background spec containers
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CalibrationSpec:
    """One calibration-source spectrum for the matrix R."""
    nuclide: str
    counts: tuple                  # tuple[float, ...]
    t_live_s: float
    A_certified_Bq_per_kg: float  # or Bq if point source
    sigma_A_rel: float = 0.05      # relative cert uncertainty
    sample_mass_kg: float = 1.0    # for specific-activity calibrations


@dataclass(frozen=True)
class BackgroundSpec:
    """Background spectrum for window-rate subtraction."""
    counts: tuple
    t_live_s: float


# ---------------------------------------------------------------------------
# Build the sensitivity matrix R_ki
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SensitivityMatrix:
    """Sensitivity matrix R[k_window, k_nuclide]."""
    R: np.ndarray                # shape (M_windows, K_nuclides)
    sigma_R: np.ndarray          # 1σ per element
    nuclides: tuple              # ordered nuclide labels
    windows_keV: tuple           # ordered window tuples
    notes: str = ""

    def __repr__(self) -> str:
        return (f"SensitivityMatrix(R.shape={self.R.shape}, "
                f"nuclides={self.nuclides})")


def build_R_matrix(
    calibration_specs: Dict[str, CalibrationSpec],
    bg_spec: BackgroundSpec,
    *,
    windows_keV: Sequence[Tuple[float, float]],
    energies_per_ch: Sequence[float],
) -> SensitivityMatrix:
    """
    Compute R_ki = (S_kr_i / t_kr − B_i / t_b) / A_k  per LSRM §12.2.

    Parameters
    ----------
    calibration_specs : map nuclide name → CalibrationSpec
    bg_spec : the background spectrum
    windows_keV : list of (E_low, E_high) energy windows
    energies_per_ch : per-channel energy axis (shared by all specs in
                      this calibration set)

    Returns
    -------
    SensitivityMatrix R, σ(R) of shape (M, K).
    """
    nuclides = list(calibration_specs.keys())
    K = len(nuclides)
    M = len(windows_keV)
    R = np.zeros((M, K), dtype=float)
    sigma_R = np.zeros((M, K), dtype=float)

    # background window rates (cps)
    B_counts = sum_in_windows(bg_spec.counts, energies_per_ch, windows_keV)
    B_rate = B_counts / bg_spec.t_live_s
    # Poisson σ on B (counts ≥ 0)
    sigma_B_rate = np.sqrt(np.maximum(B_counts, 1.0)) / bg_spec.t_live_s

    for k, n in enumerate(nuclides):
        spec = calibration_specs[n]
        S_counts = sum_in_windows(spec.counts, energies_per_ch, windows_keV)
        S_rate = S_counts / spec.t_live_s
        # σ on S
        sigma_S_rate = np.sqrt(np.maximum(S_counts, 1.0)) / spec.t_live_s

        if spec.A_certified_Bq_per_kg <= 0:
            R[:, k] = np.nan
            sigma_R[:, k] = np.nan
            continue
        net_rate = S_rate - B_rate
        R[:, k] = net_rate / spec.A_certified_Bq_per_kg
        # Uncertainty: relative propagation
        rel_S = sigma_S_rate / np.where(S_rate > 0, S_rate, 1.0)
        rel_B = sigma_B_rate / np.where(B_rate > 0, B_rate, 1.0)
        sigma_net = np.sqrt(sigma_S_rate ** 2 + sigma_B_rate ** 2)
        # Cert uncertainty
        rel_cert = spec.sigma_A_rel
        sigma_R[:, k] = np.abs(R[:, k]) * np.sqrt(
            (sigma_net / np.where(net_rate > 0, net_rate, 1.0)) ** 2
            + rel_cert ** 2
        )
    return SensitivityMatrix(
        R=R, sigma_R=sigma_R,
        nuclides=tuple(nuclides),
        windows_keV=tuple(tuple(w) for w in windows_keV),
    )


# ---------------------------------------------------------------------------
# Routine sample assay via weighted LSQ
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TemplateAssayResult:
    """Outcome of template_assay()."""
    nuclides: tuple = ()
    activities_Bq_per_kg: tuple = ()
    sigma_activities_Bq_per_kg: tuple = ()
    chi2: float = float("nan")
    chi2_per_dof: float = float("nan")
    n_windows: int = 0
    converged: bool = False
    notes: str = ""

    def by_nuclide(self) -> Dict[str, Tuple[float, float]]:
        return {n: (a, s) for n, a, s in zip(
            self.nuclides, self.activities_Bq_per_kg, self.sigma_activities_Bq_per_kg,
        )}

    def __repr__(self) -> str:
        body = ", ".join(
            f"{n}={a:.3g}±{s:.2g} Bq/kg"
            for n, a, s in zip(self.nuclides, self.activities_Bq_per_kg,
                               self.sigma_activities_Bq_per_kg)
        )
        return f"TemplateAssay(χ²/dof={self.chi2_per_dof:.2f}, {body})"


def template_assay(
    sample_counts: Sequence[float],
    sample_t_live_s: float,
    sensitivity: SensitivityMatrix,
    bg_spec: BackgroundSpec,
    *,
    energies_per_ch: Sequence[float],
    poisson_floor: float = 1.0,
) -> TemplateAssayResult:
    """
    Solve  n_i / t_p  −  B_i / t_b  =  Σ_k R_ki · A_k   for A_k.

    Weighted LSQ with Poisson weights.
    """
    R = sensitivity.R
    M, K = R.shape
    nuclides = sensitivity.nuclides

    sample_window_counts = sum_in_windows(
        sample_counts, energies_per_ch, sensitivity.windows_keV,
    )
    bg_window_counts = sum_in_windows(
        bg_spec.counts, energies_per_ch, sensitivity.windows_keV,
    )
    n_rate = sample_window_counts / sample_t_live_s
    b_rate = bg_window_counts / bg_spec.t_live_s
    y = n_rate - b_rate

    # Variance: σ²(n_rate − b_rate) = n_counts/t_p² + b_counts/t_b²
    var = (
        np.maximum(sample_window_counts, poisson_floor) / sample_t_live_s ** 2
        + np.maximum(bg_window_counts, poisson_floor) / bg_spec.t_live_s ** 2
    )
    w = 1.0 / np.maximum(var, 1e-30)
    sqrt_w = np.sqrt(w)
    A_w = R * sqrt_w[:, None]
    b_w = y * sqrt_w

    try:
        beta, _, rank, _ = np.linalg.lstsq(A_w, b_w, rcond=None)
    except Exception as exc:
        return TemplateAssayResult(
            nuclides=nuclides, converged=False,
            notes=f"lstsq failed: {exc}",
        )

    # Covariance ≈ (Rᵀ·W·R)⁻¹
    try:
        RtWR = R.T @ (R * w[:, None])
        cov = np.linalg.pinv(RtWR)
        sigma_A = np.sqrt(np.maximum(0.0, np.diag(cov)))
    except Exception:
        sigma_A = np.full(K, float("nan"))

    # Negative-clip with notes (LSRM §11 50% gate will turn small
    # negatives into upper limits downstream)
    negative_flag = bool((beta < 0).any())
    activities = np.maximum(beta, 0.0)

    y_model = R @ beta
    chi2 = float(np.sum(w * (y - y_model) ** 2))
    dof = max(1, M - K)
    chi2_dof = chi2 / dof

    notes_parts: List[str] = []
    if rank < K:
        notes_parts.append(f"rank-deficient: rank={rank} < {K}")
    if negative_flag:
        notes_parts.append("some activities clipped to 0 (negative fit)")
    if chi2_dof > 2.0:
        notes_parts.append(f"poor fit (χ²/dof={chi2_dof:.2f}>2)")
    notes = "; ".join(notes_parts)

    return TemplateAssayResult(
        nuclides=nuclides,
        activities_Bq_per_kg=tuple(float(a) for a in activities),
        sigma_activities_Bq_per_kg=tuple(float(s) for s in sigma_A),
        chi2=chi2, chi2_per_dof=chi2_dof,
        n_windows=M,
        converged=True,
        notes=notes,
    )


# ---------------------------------------------------------------------------
# Default energy windows for natural-gamma routine assay (NaI)
# ---------------------------------------------------------------------------

DEFAULT_NAI_WINDOWS_KEV: Tuple[Tuple[float, float], ...] = (
    (50.0, 100.0),       # low-energy bremsstrahlung / Pb K-XRF
    (100.0, 200.0),      # Am-241 region
    (200.0, 350.0),      # Pb-212 238, low Th line cluster
    (350.0, 480.0),      # Pb-214 351
    (480.0, 580.0),      # Cs-134 region
    (580.0, 650.0),      # Tl-208 583, Bi-214 609 boundary
    (650.0, 760.0),      # Cs-137 661
    (760.0, 880.0),      # Tl-208 weak lines
    (880.0, 1000.0),     # Ac-228 911, 969
    (1000.0, 1200.0),    # Co-60 1173, Pb-214 1120
    (1200.0, 1400.0),    # Co-60 1332
    (1400.0, 1550.0),    # K-40 1461 + Ac-228 1459 overlap
    (1550.0, 1800.0),    # Bi-214 1764
    (1800.0, 2300.0),    # high natural background tail
    (2300.0, 2700.0),    # Tl-208 2614
)


__all__ = [
    "CalibrationSpec",
    "BackgroundSpec",
    "SensitivityMatrix",
    "build_R_matrix",
    "TemplateAssayResult",
    "template_assay",
    "sum_in_windows",
    "DEFAULT_NAI_WINDOWS_KEV",
]
