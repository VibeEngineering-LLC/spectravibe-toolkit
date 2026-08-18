"""
Quasi-template method — LSRM Algorithmic Foundations §13.

Concept
-------
Instead of finding peaks one-by-one and fitting each with its own
Gaussian, LSRM §13 (квазишаблонный метод) treats the WHOLE spectrum
as the sum of per-nuclide FEP templates plus one global Compton
continuum:

    S(E)  ≈  Σ_k  A_k · σ_k(E)  +  C(E) · ξ(E)

where:

* `A_k` — activity (Bq) of the k-th candidate nuclide (the free
  parameter we want to recover);
* `σ_k(E)` — per-nuclide *template* spectrum: a synthetic shape
  consisting of every library line of nuclide k convolved with the
  detector peak-image, weighted by intensity and efficiency. This is
  what would be observed if the source contained 1 Bq of nuclide k
  alone, in the same geometry, for the same live time.
* `C(E)` — common Compton/scatter continuum, fit as a smooth piecewise
  polynomial or a low-order spline.
* `ξ(E)` — optional fixed shape factor (e.g. an empirical pile-up
  envelope, deferred for v1.16.3+).

The solution is a linear least-squares problem in (A_1, …, A_K,
c_1, …, c_N) where c_n are the continuum parameters. Poisson weights
w_i = 1/max(y_i, 1) protect against zero-count bins. The matrix is
solved by `numpy.linalg.lstsq` (rank-revealing — handles collinear
templates by returning the minimum-norm solution).

Why this matters for Gamma-1S (NaI 63×63)
-----------------------------------------
The historical peak-by-peak workflow fails when:

* Two nuclides share peaks within FWHM (e.g. Pb-214 352 keV vs
  Bi-212 727 keV near the K-edge structure — not literally
  overlapping but close enough to bias each other on NaI 63×63 with
  ~50 keV FWHM at 700 keV);
* A weak nuclide hides under the Compton plateau of a strong one
  (Eu-152 121 keV under Cs-137 Compton plateau);
* The Cs-137 661 keV peak is the only feature on top of natural
  background — Mariscotti peak search may miss it if calibration
  drifts and the priority-express trump card doesn't fire.

The quasi-template fit uses the FULL spectrum and the FULL library
simultaneously — the activity of each nuclide is constrained by ALL
its lines, ratios included via the per-template structure.

Scope of this module (v1.16.2)
-----------------------------
This is the **fitter only** — additive, not yet wired into the
orchestrator. It is documented as an opt-in path; downstream
integration into the staged pipeline lands in v1.17/v1.18.

Inputs:
    counts        : per-channel counts of the sample spectrum
    energies      : energy at each channel (keV)
    nuclide_list  : ordered list of candidate nuclides
    library       : provides {nuclide: [(E_keV, I_pct, dI_pct), ...]}
    efficiency    : callable ε(E_keV) returning photopeak efficiency
    peak_image    : optional model factory for one peak σ_k(E_i)
                    around a library line — defaults to a unit-area
                    Gaussian with FWHM from `fwhm_keV(E)`
    fwhm_keV      : callable FWHM(E_keV)
    t_live_s      : live time (seconds)
    bg_continuum_degree : polynomial degree of the global continuum

Output:
    QuasiTemplateResult with per-nuclide A_k (Bq) and σ_k from the
    diagonal of the (X·W·X^T)^{−1} matrix.

Reference
---------
LSRM Algorithmic Foundations 2022, §13 «Квазишаблонный метод»,
стр. 13-1.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Template builder — one nuclide's synthetic shape
# ---------------------------------------------------------------------------

def build_nuclide_template(
    nuclide: str,
    library_lines: Sequence[Tuple[float, float]],   # (E_keV, I_pct) pairs
    *,
    energies_keV: np.ndarray,
    fwhm_keV_fn: Callable[[float], float],
    efficiency_fn: Callable[[float], float],
    t_live_s: float,
    intensity_cutoff_pct: float = 0.1,
) -> np.ndarray:
    """
    Build σ_k(E_i): the synthetic spectrum produced by 1 Bq of nuclide
    `nuclide` over `t_live_s` seconds, for ε(E) = efficiency_fn and
    FWHM(E) = fwhm_keV_fn.

    Each library line contributes a unit-area Gaussian centred at its
    energy, scaled by (I_pct/100 · ε(E) · t_live_s).

    The total area of σ_k over E ∈ [0, +∞] equals the expected total
    *photopeak* count rate per 1 Bq. The Compton continuum is NOT
    included in σ_k — it goes into the global C(E) term.

    Lines with I_pct < intensity_cutoff_pct are skipped (negligible
    contribution to the χ² fit and they only add collinearity).
    """
    out = np.zeros_like(energies_keV, dtype=float)
    for E_keV, I_pct in library_lines:
        if I_pct < intensity_cutoff_pct:
            continue
        eps = efficiency_fn(E_keV)
        if eps <= 0:
            continue
        fwhm = fwhm_keV_fn(E_keV)
        if fwhm <= 0:
            continue
        sigma = fwhm / 2.355
        amplitude_counts = (I_pct / 100.0) * eps * t_live_s
        # Unit-area Gaussian × amplitude_counts. The discretisation
        # error at FWHM ≥ 2 channels is negligible.
        z = (energies_keV - E_keV) / sigma
        out += amplitude_counts * np.exp(-0.5 * z * z) / (sigma * math.sqrt(2.0 * math.pi))
    return out


# ---------------------------------------------------------------------------
# Continuum basis — polynomial in E
# ---------------------------------------------------------------------------

def build_continuum_basis(
    energies_keV: np.ndarray, degree: int,
) -> np.ndarray:
    """
    Return (N_channels, degree+1) matrix of polynomial basis vectors
    in the normalised energy x = E / E_max (numerical conditioning).
    """
    E_max = float(np.max(energies_keV))
    if E_max <= 0:
        E_max = 1.0
    x = energies_keV / E_max
    cols = [x ** k for k in range(degree + 1)]
    return np.column_stack(cols)


# ---------------------------------------------------------------------------
# Quasi-template fit
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class QuasiTemplateResult:
    """Outcome of quasi_template_fit()."""
    nuclides: tuple = ()                  # tuple[str, ...]
    activities_Bq: tuple = ()             # tuple[float, ...] aligned to nuclides
    sigma_activities_Bq: tuple = ()       # 1σ on activities
    continuum_coeffs: tuple = ()          # polynomial continuum params
    chi2: float = float("nan")
    chi2_per_dof: float = float("nan")
    n_channels: int = 0
    n_nuclides: int = 0
    converged: bool = False
    notes: str = ""

    def by_nuclide(self) -> Dict[str, Tuple[float, float]]:
        """Map nuclide → (A_Bq, sigma_A_Bq)."""
        return {n: (a, s) for n, a, s in zip(
            self.nuclides, self.activities_Bq, self.sigma_activities_Bq,
        )}

    def __repr__(self) -> str:
        body = ", ".join(
            f"{n}={a:.2e}±{s:.1e}"
            for n, a, s in zip(self.nuclides, self.activities_Bq, self.sigma_activities_Bq)
        )
        return (f"QuasiTemplate(χ²/dof={self.chi2_per_dof:.2f}, "
                f"converged={self.converged}, {body})")


def quasi_template_fit(
    counts: Sequence[float],
    energies_keV: Sequence[float],
    nuclide_lines: Dict[str, Sequence[Tuple[float, float]]],
    *,
    fwhm_keV_fn: Callable[[float], float],
    efficiency_fn: Callable[[float], float],
    t_live_s: float,
    continuum_degree: int = 4,
    poisson_floor: float = 1.0,
    E_min_keV: Optional[float] = None,
    E_max_keV: Optional[float] = None,
) -> QuasiTemplateResult:
    """
    Run the LSRM §13 quasi-template fit.

    Parameters
    ----------
    counts        : per-channel counts (length N)
    energies_keV  : per-channel energy (length N)
    nuclide_lines : map nuclide → [(E_keV, I_pct), ...] from the
                    library, intensities in percent.
    fwhm_keV_fn   : detector FWHM(E) calibration
    efficiency_fn : detector efficiency ε(E) calibration
    t_live_s      : live time
    continuum_degree : polynomial degree of the C(E) continuum
                       (default 4 — LSRM-recommended for natural-bg
                       spectra; raise to 6 for very-low-stat cases).
    poisson_floor : minimum weight denominator (avoid /0 on empty bins)
    E_min_keV, E_max_keV : optional restriction of the fit window.

    Returns
    -------
    QuasiTemplateResult.
    """
    counts = np.asarray(list(counts), dtype=float)
    energies = np.asarray(list(energies_keV), dtype=float)
    if len(counts) != len(energies):
        return QuasiTemplateResult(
            converged=False,
            notes="counts/energies length mismatch",
        )
    # ── Restrict to fit window ─────────────────────────────────────
    mask = np.ones(len(energies), dtype=bool)
    if E_min_keV is not None:
        mask &= (energies >= E_min_keV)
    if E_max_keV is not None:
        mask &= (energies <= E_max_keV)
    if mask.sum() < 10:
        return QuasiTemplateResult(
            converged=False,
            notes=f"too few channels in fit window: {int(mask.sum())}",
        )
    E_fit = energies[mask]
    y_fit = counts[mask]

    nuclides = list(nuclide_lines.keys())
    if not nuclides:
        return QuasiTemplateResult(
            converged=False, notes="no candidate nuclides supplied",
        )

    # ── Build the model matrix ──────────────────────────────────────
    # Columns = K nuclide templates + (degree+1) continuum basis
    K = len(nuclides)
    sigma_columns: List[np.ndarray] = []
    for n in nuclides:
        col = build_nuclide_template(
            n, nuclide_lines[n],
            energies_keV=E_fit,
            fwhm_keV_fn=fwhm_keV_fn,
            efficiency_fn=efficiency_fn,
            t_live_s=t_live_s,
        )
        sigma_columns.append(col)
    cont = build_continuum_basis(E_fit, continuum_degree)

    X = np.column_stack(sigma_columns + [cont])
    # Poisson weights w_i = 1/max(y_i, poisson_floor)
    w = 1.0 / np.maximum(y_fit, poisson_floor)
    sqrt_w = np.sqrt(w)

    # Weighted LSQ: minimise || sqrt(w) · (y - X·β) ||²
    A_w = X * sqrt_w[:, None]
    b_w = y_fit * sqrt_w

    try:
        beta, residuals_lsq, rank, sv = np.linalg.lstsq(A_w, b_w, rcond=None)
    except Exception as exc:
        return QuasiTemplateResult(
            converged=False, notes=f"lstsq failed: {exc}",
        )

    # ── Activities and σ from covariance matrix ─────────────────────
    activities = beta[:K]
    continuum_coeffs = beta[K:]
    # Covariance ≈ (X^T W X)⁻¹ ; diag gives σ². Use Moore-Penrose on
    # rank-deficient cases.
    try:
        XtWX = X.T @ (X * w[:, None])
        cov = np.linalg.pinv(XtWX)
        sigma_A = np.sqrt(np.maximum(0.0, np.diag(cov)[:K]))
    except Exception:
        sigma_A = np.full(K, float("nan"))

    # Activities can come out negative if a template was forced into a
    # noisy bin — clip to 0 with a notes flag (LSRM §11 50% gate will
    # then turn them into upper limits downstream).
    negative_flag = bool((activities < 0).any())
    activities = np.maximum(activities, 0.0)

    # ── χ² statistics ───────────────────────────────────────────────
    y_model = X @ beta
    chi2 = float(np.sum(w * (y_fit - y_model) ** 2))
    dof = max(1, len(y_fit) - len(beta))
    chi2_dof = chi2 / dof

    notes_parts = []
    if rank < X.shape[1]:
        notes_parts.append(f"rank-deficient: rank={rank} < {X.shape[1]}")
    if negative_flag:
        notes_parts.append("some activities clipped to 0 (negative fit)")
    if chi2_dof > 2.0:
        notes_parts.append(f"poor fit (χ²/dof={chi2_dof:.2f}>2)")
    notes = "; ".join(notes_parts)

    return QuasiTemplateResult(
        nuclides=tuple(nuclides),
        activities_Bq=tuple(float(a) for a in activities),
        sigma_activities_Bq=tuple(float(s) for s in sigma_A),
        continuum_coeffs=tuple(float(c) for c in continuum_coeffs),
        chi2=chi2, chi2_per_dof=chi2_dof,
        n_channels=int(mask.sum()),
        n_nuclides=K,
        converged=True,
        notes=notes,
    )


__all__ = [
    "build_nuclide_template",
    "build_continuum_basis",
    "QuasiTemplateResult",
    "quasi_template_fit",
]
