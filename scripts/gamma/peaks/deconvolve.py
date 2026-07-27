"""
Multiplet deconvolution (Phase 2.1b, closes K-05).

Lsrm Algorithmic Foundations §9 / Gilmore §9.7 prescribe a strict
"identification-first" deconvolution:

  - Peak **positions** are FIXED from the library energies of confirmed
    nuclides (no centroid floating — they are known a priori).
  - Peak **widths** (σ) are FIXED from the calibrated FWHM(E) curve.
  - Only **areas** (one per component) and the **continuum** beneath
    the multiplet are free.

With positions and widths frozen, the fit model becomes **linear** in
its free parameters. This converts what would otherwise be a fragile
nonlinear Levenberg-Marquardt fit into a small linear non-negative
least-squares problem — well-posed whenever the components are
separated by more than about σ/2 in channel space.

Continuum model:

  *Linear*       : B(x) = β₀ + β₁·(x − x_mid)
  *step + linear*: B(x) = β₀ + β₁·(x − x_mid) + β_step·S(x)

  where S(x) = 0.5·erfc((x − x_step) / (σ_step·√2)) is a smooth step
  that goes from 1 well below `x_step` to 0 well above it. `x_step` is
  the area-weighted average of the component centroids (after a first
  pass) and `σ_step` is the largest component FWHM/2.355. β_step is
  constrained ≥ 0: the Compton continuum from a photopeak adds counts
  on the LOW-energy side, so the step is physically "down" as channel
  increases.

The free-parameter vector is

    p = [A₁, A₂, …, A_n, β₀, β₁]                      (linear continuum)
    p = [A₁, A₂, …, A_n, β₀, β₁, β_step]              (step + linear)

with A_k ≥ 0 and β_step ≥ 0 (when present). β₀ and β₁ are unconstrained
(continuum can rise OR fall with channel).

The Gaussian shape used is normalised so that the area parameter A_k
is literally the integral of the component over the entire real line:

    g_k(x) = (1 / (σ_k · √(2π))) · exp(-(x − c_k)² / (2 σ_k²))

so that A_k is what activity calculations downstream want.

Limitations and graceful fallbacks:

  - Components closer than 0.5·σ in channel space are physically
    unresolvable. We still attempt the fit but flag the result with
    `notes` and a `degenerate=True` marker on the offending pair so
    callers can decide to merge.
  - If `scipy.optimize.lsq_linear` is unavailable, falls back to
    `numpy.linalg.lstsq` followed by area-clipping to ≥ 0 and a single
    re-solve with the offending components removed.
  - Single-component "multiplet" is accepted and reduces to a regular
    single-Gaussian fit with a known shape. The result format stays
    the same so callers don't branch.

References:
  - Lsrm Algorithmic Foundations 2022, §9 (model-based multiplet
    deconvolution under fixed-position constraints)
  - Gilmore & Joss, Practical Gamma-ray Spectrometry, 3rd Ed., §9.7
    (continuum modelling — step + linear under the photopeak)
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence

import numpy as np


SQRT_2PI = math.sqrt(2.0 * math.pi)
SQRT_2 = math.sqrt(2.0)


# ============================================================================
# Data classes
# ============================================================================

@dataclass(frozen=True)
class MultipletComponent:
    """
    One component of a multiplet — a library line at a known energy
    mapped to its expected channel and FWHM.

    All four fields are fixed for the duration of the fit.
    """
    nuclide: str
    line_E_keV: float
    library_I_pct: float            # intensity from the library (informational)
    center_channel: float           # mapped from energy_cal
    fwhm_channels: float            # from fwhm_at_channel

    @property
    def sigma_channels(self) -> float:
        return self.fwhm_channels / 2.355


@dataclass(frozen=True)
class DeconvolutionResult:
    """
    Outcome of one multiplet deconvolution.

    `areas[k]` and `area_uncertainties[k]` correspond to
    `components[k]`. The 1-σ uncertainties come from the diagonal of
    the parameter covariance matrix; for correlated components callers
    should consult `covariance` (the full (n+m)×(n+m) matrix where m is
    the number of continuum parameters).
    """
    components: tuple                       # tuple[MultipletComponent, ...]
    areas: tuple                            # tuple[float, ...]
    area_uncertainties: tuple               # tuple[float, ...]
    continuum_params: tuple                 # tuple[float, ...] — β₀, β₁ [, β_step]
    continuum_model: str                    # "linear" or "step_linear"
    chi2_per_dof: float
    n_dof: int
    roi_low_ch: int
    roi_high_ch: int
    gross_counts: float                     # total counts in ROI
    converged: bool
    method: str                             # "lsq_linear" or "lstsq_fallback"
    degenerate_pairs: tuple = ()            # tuple[tuple[int, int], ...]
    covariance: Optional[np.ndarray] = None
    notes: str = ""
    # F-134 / v1.17.7 — pre-computed overlay arrays для PNG / HTML рендера.
    # Заполняются _coupled_to_deconv_result из CoupledFitResult; при их
    # наличии downstream-рендерер использует эти массивы напрямую вместо
    # перевычисления continuum/гауссов в канальных координатах (что для
    # формы Гаусс+tail+step делается неверно). None при заполнении через
    # старый канальный путь — тогда рендерер падает к legacy-логике.
    overlay_E_keV: Optional[tuple] = None             # tuple[float, ...]
    overlay_data: Optional[tuple] = None              # сырой y(E)
    overlay_continuum: Optional[tuple] = None         # B(E)
    overlay_total: Optional[tuple] = None             # B + Σ comp(E)
    overlay_components: Optional[tuple] = None        # tuple[tuple[float, ...]] — per comp B+only-this-comp
    # F-145 / v1.17.8 — Phase A free-centroid side-fit поля для self-calibration.
    # Заполняются _coupled_to_deconv_result из CoupledFitResult. Используются
    # outer pipeline'ом (staged_pipeline) для Phase B+C decision и JSON-report
    # блока `self_calibration`. None / пустой список — Phase A не запускалась.
    centroid_shifts_keV: tuple = ()                    # выровнено по components
    phase_A_chi2_per_dof: Optional[float] = None       # χ²/ν Phase A (free)
    phase_A_converged: bool = False
    cluster_id: str = ""                               # ⟵ из CoupledFitResult.id
    # F-387.1 / v1.18.26.1 — phantom anchors: компоненты, отрезанные
    # top-K cap (after Rayleigh-CC split). Они НЕ в `components` /
    # `areas` (исключены из fit'а), но сохраняются здесь как
    # `tuple[MultipletComponent, ...]` для downstream
    # reporting/evidence/identification. Reporting layer (json_report)
    # подмешивает их в JSON `components` массив с
    # `peak_area_source="library_anchor_phantom"`.
    phantom_components: tuple = ()
    # F-392.1 / v1.18.29 — propagate multi-step continuum diagnostics из
    # CoupledFitResult в downstream JSON reporter. multi_step_anchors —
    # tuple[(E_keV, σ_step), ...] для intense-anchor компонент при
    # continuum_model == "step_linear_multi"; пустой кортеж при других
    # моделях. multi_step_intensity_threshold_pct — порог library_I_pct
    # (typ. 4.0%); None для non-multi continuum.
    multi_step_anchors: tuple = ()
    multi_step_intensity_threshold_pct: Optional[float] = None

    def total_fit_area(self) -> float:
        return float(sum(self.areas))

    def area_by_nuclide(self) -> dict:
        """Sum areas of components belonging to the same nuclide."""
        out: dict = {}
        for comp, a in zip(self.components, self.areas):
            out[comp.nuclide] = out.get(comp.nuclide, 0.0) + a
        return out

    def __repr__(self) -> str:
        parts = ", ".join(
            f"{c.nuclide}@{c.line_E_keV:.1f}={a:.0f}±{u:.0f}"
            for c, a, u in zip(self.components, self.areas, self.area_uncertainties)
        )
        return (f"DeconvolutionResult(n={len(self.components)}, "
                f"χ²/ν={self.chi2_per_dof:.2f}, {parts})")


# ============================================================================
# Core fit
# ============================================================================

def _gaussian_normalised(x: np.ndarray, c: float, sigma: float) -> np.ndarray:
    """Unit-area Gaussian sampled at channel grid `x`."""
    return np.exp(-((x - c) / sigma) ** 2 * 0.5) / (sigma * SQRT_2PI)


def _smooth_step(x: np.ndarray, x_step: float, sigma_step: float) -> np.ndarray:
    """
    Smooth step: 1 well below `x_step`, 0 well above.

    Implemented as 0.5·erfc((x - x_step) / (σ·√2)).
    """
    from math import erfc as _erfc  # scalar fallback
    # numpy vectorised erfc via scipy if available, else loop
    try:
        from scipy.special import erfc as _vec_erfc  # type: ignore
        return 0.5 * _vec_erfc((x - x_step) / (sigma_step * SQRT_2))
    except ImportError:
        out = np.empty_like(x, dtype=np.float64)
        scale = 1.0 / (sigma_step * SQRT_2)
        for i, xi in enumerate(x):
            out[i] = 0.5 * _erfc((xi - x_step) * scale)
        return out


def _build_design_matrix(
    x: np.ndarray,
    components: Sequence[MultipletComponent],
    *,
    continuum: str,
    x_mid: float,
    x_step: float,
    sigma_step: float,
) -> tuple[np.ndarray, list]:
    """
    Build the (n_pixels × n_params) design matrix for the linear model.

    Columns: [g_1, g_2, …, g_n, 1, (x - x_mid)] for linear
             [..., S(x)] for step_linear

    Returns (A, param_names) where `param_names` is a list of strings
    naming each column for diagnostics.
    """
    cols = []
    names: list = []
    for c in components:
        cols.append(_gaussian_normalised(x, c.center_channel, c.sigma_channels))
        names.append(f"area:{c.nuclide}@{c.line_E_keV:.1f}")
    cols.append(np.ones_like(x))
    names.append("continuum:beta0")
    cols.append(x - x_mid)
    names.append("continuum:beta1_slope")
    if continuum == "step_linear":
        cols.append(_smooth_step(x, x_step, sigma_step))
        names.append("continuum:beta_step")
    A = np.column_stack(cols)
    return A, names


def _solve_constrained(
    A: np.ndarray,
    y: np.ndarray,
    sigma_y: np.ndarray,
    *,
    n_components: int,
    has_step: bool,
) -> tuple[np.ndarray, str, bool]:
    """
    Solve weighted LSQ with bounds:
      - area columns (first n_components) ≥ 0
      - step-continuum column (last, if has_step) ≥ 0
      - other columns unconstrained

    Returns (params, method_label, converged).
    """
    w = 1.0 / sigma_y
    Aw = A * w[:, None]
    yw = y * w

    n_params = A.shape[1]
    lb = np.full(n_params, -np.inf)
    ub = np.full(n_params, np.inf)
    lb[:n_components] = 0.0
    if has_step:
        lb[-1] = 0.0

    try:
        from scipy.optimize import lsq_linear  # type: ignore
        res = lsq_linear(Aw, yw, bounds=(lb, ub), method="trf",
                         tol=1e-10, max_iter=200)
        return res.x, "lsq_linear", bool(res.success)
    except ImportError:
        pass

    # Fallback: unconstrained lstsq; if any bound is violated, drop
    # that variable (set it to zero) and re-solve. One iteration is
    # usually enough for the small problems we have.
    sol, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
    violations = np.zeros(n_params, dtype=bool)
    violations[:n_components] = sol[:n_components] < 0
    if has_step:
        violations[-1] = sol[-1] < 0
    if not violations.any():
        return sol, "lstsq_fallback", True

    keep = ~violations
    sub = Aw[:, keep]
    if sub.shape[1] == 0:
        return np.zeros(n_params), "lstsq_fallback", False
    sub_sol, *_ = np.linalg.lstsq(sub, yw, rcond=None)
    full = np.zeros(n_params)
    full[keep] = sub_sol
    return full, "lstsq_fallback", True


def _parameter_covariance(
    A: np.ndarray,
    sigma_y: np.ndarray,
    residuals: np.ndarray,
    n_dof: int,
) -> Optional[np.ndarray]:
    """
    Estimate parameter covariance from the weighted normal matrix.
    Returns None if the matrix is singular.
    """
    w = 1.0 / sigma_y
    Aw = A * w[:, None]
    # AUDIT-F5 (2026-06-25): cov через SVD от Aw, не через inv(Aw.T @ Aw).
    # Эквивалентно матем., но cond(Aw) вместо cond²(Aw) → стабильные σ
    # для тесных мультиплетов. Центральные параметры здесь не считаются
    # этим путём (см. вызывающие места), поэтому правка влияет только
    # на репортируемую ковариацию.
    try:
        _U, _s, _Vt = np.linalg.svd(Aw, full_matrices=False)
    except np.linalg.LinAlgError:
        return None
    if _s.size == 0 or _s[-1] <= 0:
        return None
    normal_inv = (_Vt.T * (1.0 / (_s * _s))) @ _Vt
    chi2 = float(np.sum((residuals * w) ** 2))
    scale = chi2 / max(1, n_dof)
    return normal_inv * scale


def deconvolve_multiplet(
    counts,
    *,
    components: Sequence[MultipletComponent],
    continuum: str = "step_linear",
    roi_window_factor: float = 2.5,
    roi_low_ch: Optional[int] = None,
    roi_high_ch: Optional[int] = None,
    degenerate_separation_sigma: float = 0.5,
) -> DeconvolutionResult:
    """
    Fit a multiplet of fixed-position Gaussians plus a continuum.

    Args:
        counts: 1-D array of channel counts (full spectrum).
        components: ordered list of `MultipletComponent` — one per
            library line participating in the multiplet. Positions
            (`center_channel`) and widths (`fwhm_channels`) are
            treated as exact constraints; only the areas float.
        continuum: "linear" or "step_linear" (default).
        roi_window_factor: when `roi_low_ch`/`roi_high_ch` are not
            given, the ROI is set to
                [min(c_k) - F·FWHM_max, max(c_k) + F·FWHM_max]
            with F = `roi_window_factor`. Default 2.5 captures ≥99%
            of each Gaussian.
        roi_low_ch, roi_high_ch: explicit ROI overrides (half-open,
            like Python slicing).
        degenerate_separation_sigma: any pair of components closer than
            this fraction of the smaller σ is flagged as degenerate
            (still fit, but reported separately so callers can merge).

    Returns:
        `DeconvolutionResult` with one area per input component, the
        fitted continuum parameters, χ²/dof, and the full covariance.
    """
    if not components:
        raise ValueError("at least one MultipletComponent is required")
    if continuum not in ("linear", "step_linear"):
        raise ValueError(f"unknown continuum model: {continuum!r}")

    counts_arr = np.asarray(counts, dtype=np.float64)
    n_ch = len(counts_arr)

    # ----- ROI selection -----
    fwhm_max = max(c.fwhm_channels for c in components)
    if roi_low_ch is None:
        roi_low_ch = int(math.floor(
            min(c.center_channel for c in components) - roi_window_factor * fwhm_max
        ))
    if roi_high_ch is None:
        roi_high_ch = int(math.ceil(
            max(c.center_channel for c in components) + roi_window_factor * fwhm_max
        )) + 1  # half-open
    roi_low_ch = max(0, roi_low_ch)
    roi_high_ch = min(n_ch, roi_high_ch)
    if roi_high_ch - roi_low_ch < len(components) + 3:
        # not enough points to constrain n areas + continuum
        return DeconvolutionResult(
            components=tuple(components),
            areas=tuple([0.0] * len(components)),
            area_uncertainties=tuple([0.0] * len(components)),
            continuum_params=tuple(),
            continuum_model=continuum,
            chi2_per_dof=float("inf"),
            n_dof=0,
            roi_low_ch=roi_low_ch,
            roi_high_ch=roi_high_ch,
            gross_counts=float(counts_arr[roi_low_ch:roi_high_ch].sum())
                         if roi_high_ch > roi_low_ch else 0.0,
            converged=False,
            method="not_attempted",
            notes=f"ROI too narrow: {roi_high_ch - roi_low_ch} channels "
                  f"for {len(components)} components + continuum",
        )

    x = np.arange(roi_low_ch, roi_high_ch, dtype=np.float64)
    y = counts_arr[roi_low_ch:roi_high_ch]
    gross = float(y.sum())

    # Poisson error model with minimum floor (1 count) so empty bins
    # don't blow up the weighted residual.
    sigma_y = np.sqrt(np.maximum(y, 1.0))

    # ----- Step centre / width -----
    weights = np.maximum([c.fwhm_channels for c in components], 1e-6)
    x_step = float(sum(c.center_channel * w for c, w in zip(components, weights))
                   / sum(weights))
    sigma_step = fwhm_max / 2.355
    x_mid = 0.5 * (roi_low_ch + roi_high_ch - 1)

    # ----- Design matrix -----
    A, _names = _build_design_matrix(
        x, components,
        continuum=continuum, x_mid=x_mid,
        x_step=x_step, sigma_step=sigma_step,
    )

    has_step = (continuum == "step_linear")
    params, method_label, converged = _solve_constrained(
        A, y, sigma_y,
        n_components=len(components),
        has_step=has_step,
    )

    # ----- Residuals, χ², covariance -----
    model = A @ params
    residuals = y - model
    n_params = A.shape[1]
    n_dof = max(1, len(y) - n_params)
    chi2_per_dof = float(np.sum((residuals / sigma_y) ** 2)) / n_dof
    cov = _parameter_covariance(A, sigma_y, residuals, n_dof)

    # ----- Unpack -----
    n_comp = len(components)
    areas = tuple(float(v) for v in params[:n_comp])
    if cov is not None:
        area_uncs = tuple(
            float(math.sqrt(max(0.0, cov[i, i]))) for i in range(n_comp)
        )
    else:
        # Pessimistic Poisson scaling
        area_uncs = tuple(float(math.sqrt(max(a, 1.0))) for a in areas)
    continuum_params = tuple(float(v) for v in params[n_comp:])

    # ----- Degenerate-pair detection -----
    degenerate = []
    for i in range(n_comp):
        for j in range(i + 1, n_comp):
            c1, c2 = components[i], components[j]
            sep = abs(c1.center_channel - c2.center_channel)
            min_sigma = min(c1.sigma_channels, c2.sigma_channels)
            if sep < degenerate_separation_sigma * min_sigma:
                degenerate.append((i, j))

    notes_parts = []
    if degenerate:
        notes_parts.append(
            f"{len(degenerate)} degenerate pair(s) closer than "
            f"{degenerate_separation_sigma}·σ — areas of these are "
            f"strongly correlated, consult covariance"
        )
    if chi2_per_dof > 1.5:
        notes_parts.append(
            f"χ²/ν = {chi2_per_dof:.2f} > 1.5; possible missing "
            f"component or response-model inadequacy"
        )
    if chi2_per_dof < 0.5:
        notes_parts.append(
            f"χ²/ν = {chi2_per_dof:.2f} < 0.5; over-fitting or "
            f"σ_y over-estimated"
        )

    # F-373 — populate overlay_* arrays in channel space so the HTML/PNG
    # renderer can avoid the legacy fallback (which used cp[0] + cp[1]·ch
    # and ignored the step term, producing wrong continuum on step_linear).
    # We reconstruct continuum and per-component curves directly from the
    # solved design-matrix columns, then clamp continuum to ≥0 (physically
    # counts can't be negative; a negative-slope linear continuum could
    # otherwise dive below zero on long ROIs).
    try:
        n_cont = len(continuum_params)
        # continuum = β₀ * 1 + β₁ * (x - x_mid) [+ β_step * step(x)]
        cont_curve = (
            float(params[n_comp]) * np.ones_like(x, dtype=np.float64)
            + float(params[n_comp + 1]) * (x - x_mid)
        )
        if has_step and n_cont >= 3:
            cont_curve = cont_curve + float(params[n_comp + 2]) * _smooth_step(
                x, x_step, sigma_step,
            )
        # F-373 — clamp continuum to ≥0 (negative counts unphysical)
        cont_curve = np.maximum(cont_curve, 0.0)
        # Per-component overlays = continuum + only-this-Gaussian
        per_comp_overlays = []
        for k in range(n_comp):
            g_k = float(params[k]) * _gaussian_normalised(
                x, components[k].center_channel,
                components[k].sigma_channels,
            )
            per_comp_overlays.append(tuple(float(v) for v in cont_curve + g_k))
        # Total = continuum + Σ_k area_k · g_k. Use clamped continuum so
        # total is reconstructed self-consistently with what's shown.
        total_curve = cont_curve.copy()
        for k in range(n_comp):
            total_curve += float(params[k]) * _gaussian_normalised(
                x, components[k].center_channel,
                components[k].sigma_channels,
            )
        # NOTE: overlay_E_keV is **intentionally left unset**. deconvolve_multiplet
        # operates in channel space and has no Spectrum reference for the keV
        # mapping. The renderer (interactive_html._build_multiplets_data)
        # checks overlay_data/continuum/total/components — if those are
        # present, it builds the E_keV axis from spec.channel_to_energy on
        # [roi_low_ch, roi_high_ch), which is what we want.
        overlay_data_t = tuple(float(v) for v in y)
        overlay_cont_t = tuple(float(v) for v in cont_curve)
        overlay_tot_t = tuple(float(v) for v in total_curve)
        overlay_comps_t = tuple(per_comp_overlays)
    except Exception:
        overlay_data_t = None
        overlay_cont_t = None
        overlay_tot_t = None
        overlay_comps_t = None

    return DeconvolutionResult(
        components=tuple(components),
        areas=areas,
        area_uncertainties=area_uncs,
        continuum_params=continuum_params,
        continuum_model=continuum,
        chi2_per_dof=chi2_per_dof,
        n_dof=n_dof,
        roi_low_ch=roi_low_ch,
        roi_high_ch=roi_high_ch,
        gross_counts=gross,
        converged=converged,
        method=method_label,
        degenerate_pairs=tuple(degenerate),
        covariance=cov,
        notes="; ".join(notes_parts),
        overlay_data=overlay_data_t,
        overlay_continuum=overlay_cont_t,
        overlay_total=overlay_tot_t,
        overlay_components=overlay_comps_t,
    )


# ============================================================================
# F-392 / v1.18.27 — multi-step continuum auto-selection
# ============================================================================

def _f392_auto_select_continuum(
    components_with_I,
    base_continuum: str,
    *,
    roi_e_span_keV: float = 0.0,
    e_span_threshold_keV: float = 200.0,
    min_intense_anchors: int = 3,
    # F-392.1 / v1.18.27.1 — default понижен 5.0 → 4.0% (см. полный комментарий
    # в `coupled_intensity_fit` kwarg). Real-data: Ac-228 463 (4.4%) в Th-232
    # M3 PROD теперь quotient'ится как valid anchor вместе с Tl-208 510/583.
    intense_threshold_pct: float = 4.0,
    min_separation_keV: float = 40.0,
) -> str:
    """F-392 / v1.18.27 — выбрать continuum-модель для coupled-fit.

    Возвращает "step_linear_multi" когда выполнены ВСЕ:
      - base_continuum == "step_linear" (auto-promote ONLY from step_linear,
        чтобы не ломать back-compat для forced "linear" callers);
      - E_span ≥ e_span_threshold_keV. E_span = max(roi_e_span_keV,
        max(E_comp) - min(E_comp)) — ROI bounds доминируют, потому что
        широкий display-window expand'ит cluster в широкий ROI даже
        если внутренний spread компонент узкий (M4 PROD: components
        463-583 кэВ, но ROI 350-700 кэВ из-за display-window F-374);
      - ≥ min_intense_anchors компонент с library_I_pct ≥ intense_threshold_pct;
      - после слипания anchors-ближе-min_separation_keV остаётся
        ≥ min_intense_anchors anchors.

    Иначе возвращает base_continuum без изменения.

    `components_with_I` — iterable of (E_keV, I_pct).
    `roi_e_span_keV` — ширина ROI display-window в кэВ (если известна);
    при 0.0 проверяется только component span.

    Использование: в каждой call-точке coupled_intensity_fit (auto-cluster
    path в apply_multiplet_deconvolution и forced-cluster path в
    deconvolve_identified_multiplets).
    """
    if base_continuum != "step_linear":
        return base_continuum
    items = [(float(e), float(i)) for e, i in components_with_I]
    if not items:
        return base_continuum
    e_vals = [e for e, _ in items]
    comp_span = max(e_vals) - min(e_vals)
    e_span = max(float(roi_e_span_keV), comp_span)
    if e_span < e_span_threshold_keV:
        return base_continuum
    intense = sorted(
        ((e, i) for e, i in items if i >= intense_threshold_pct),
        key=lambda t: t[0],
    )
    if len(intense) < min_intense_anchors:
        return base_continuum
    # Схлопывание соседних anchors ближе min_separation_keV — оставить
    # сильнейшего из пары.
    merged = []
    for e_c, i_c in intense:
        if merged and (e_c - merged[-1][0]) < min_separation_keV:
            if i_c > merged[-1][1]:
                merged[-1] = (e_c, i_c)
        else:
            merged.append((e_c, i_c))
    if len(merged) < min_intense_anchors:
        return base_continuum
    return "step_linear_multi"


# ============================================================================
# F-391 / v1.18.27 — S/N significance gate for multiplet members
# ============================================================================

def _f391_peak_snr(line_match) -> float:
    """F-391 / v1.18.27 — compute peak Signal-to-Noise для одного LineMatch.

    Возвращает 0.0 для library-only anchors (F-381 ``peak_area_source ==
    "library_anchor"`` или F-387.1 ``"library_anchor_phantom"``) — у них
    нет измеренного peak_area, значит нет реального signal в этом
    energy bin. Phantom-status сохраняется через S/N gate.

    BUG-3 (2026-06-02): ``library_anchor_strong`` возвращает ``+inf`` —
    это библиотечные линии с intensity ≥ 5%·I_max_in_window, которые
    обязаны участвовать в fit как active components (не phantoms).
    Без этого strong lines (Ac-228 209/270 в Th-232 M3) поглощались
    в continuum, а fit натягивал амплитуды на weak соседей.

    Иерархия источников S/N (порядок предпочтения):
      1. ``peak_area / peak_area_uncertainty`` — канонический S/N
         «net area significance» (Gilmore & Joss 3rd Ed. §5.5, формула
         5.21 — Currie's L_C variant из net area + variance). Используется
         если оба значения присутствуют и uncertainty > 0.
      2. ``significance_currie`` — Currie L_C-style значимость S/sqrt(B) из peak
         finder (``gamma.peaks.search.FoundPeak.significance``). Fallback
         когда peak_area_uncertainty отсутствует / 0 (e.g. lsrm_peaks_table
         без uncertainty).
      3. Когда нет ни area+unc, ни significance_currie, ни явного «anchor» source-метки,
         возвращаем ``float('inf')`` — индикация «информация недоступна,
         НЕ применяй S/N gate» (back-compat для test fixtures и
         legacy code paths без populated значимости).

    References:
      - Gilmore & Joss, 3rd Ed., §5.5 «Peak detection limits», стр. 124:
        S/N ≡ net_area / σ_net_area для FEP-significance.
      - Currie 1968 — L_C критерий, эквивалентный 3σ при k=3.
    """
    src = str(getattr(line_match, "peak_area_source", "") or "")
    if src in ("library_anchor", "library_anchor_phantom"):
        return 0.0
    # BUG-3: strong library lines bypass S/N gate.
    if src == "library_anchor_strong":
        return float("inf")
    area = getattr(line_match, "peak_area", None)
    area_unc = getattr(line_match, "peak_area_uncertainty", None)
    if area is not None and area_unc is not None:
        try:
            a = float(area)
            u = float(area_unc)
            if u > 0:
                # area может быть отрицательной (continuum mismatch),
                # тогда S/N=0 — реального signal нет.
                if a <= 0:
                    return 0.0
                return a / u
        except (TypeError, ValueError):
            pass
    sig = getattr(line_match, "significance_currie", None)
    if sig is not None:
        try:
            sig_f = float(sig)
            if sig_f > 0:
                return sig_f
        except (TypeError, ValueError):
            pass
    # No measurement info available — return inf to disable gate
    # (back-compat для test fixtures и legacy paths).
    return float("inf")


def _f391_mark_phantom(line_match):
    """F-391 — пометить LineMatch как phantom (без peak_area), сохранив
    в cluster.components для evidence. Возвращает новый LineMatch с
    ``peak_area_source="library_anchor_phantom"``.

    Используется dataclass.replace для frozen LineMatch; для mutable
    тестовых fake'ов — прямая мутация атрибута.
    """
    from dataclasses import replace as _dc_replace_lm
    try:
        return _dc_replace_lm(
            line_match,
            peak_area=None,
            peak_area_uncertainty=None,
            peak_area_source="library_anchor_phantom",
        )
    except Exception:
        try:
            line_match.peak_area = None
            line_match.peak_area_uncertainty = None
            line_match.peak_area_source = "library_anchor_phantom"
        except Exception:
            pass
        return line_match


# ============================================================================
# F-447 helpers — guest-only Phase 1 weak-line filter (ROI-owner aware)
# ============================================================================

def _f447_identify_roi_owner(roi_matches):
    """F-447 — owner ROI = nuclide с max sum(I_gamma); tie → dom-line>30%.

    Возвращает идентификатор нуклида-владельца ROI либо None для
    multi-owner случая (равные суммы I_gamma без однозначного
    разрешения по dominant line I_gamma>30%). None caller трактует
    как «применить guest-фильтр ко всем линиям ROI».
    """
    if not roi_matches:
        return None
    nuc_I_sum: dict = {}
    nuc_has_dom: dict = {}
    for m in roi_matches:
        nuc = getattr(m, "nuclide", None)
        if nuc is None:
            continue
        I = float(getattr(m, "library_I_pct", 0.0) or 0.0)
        nuc_I_sum[nuc] = nuc_I_sum.get(nuc, 0.0) + I
        if I > 30.0:
            nuc_has_dom[nuc] = True
    if not nuc_I_sum:
        return None
    sorted_n = sorted(nuc_I_sum.items(), key=lambda kv: kv[1], reverse=True)
    top_nuc, top_sum = sorted_n[0]
    if len(sorted_n) >= 2 and abs(top_sum - sorted_n[1][1]) < 1e-6:
        tie = [n for n, s in sorted_n if abs(s - top_sum) < 1e-6]
        with_dom = [n for n in tie if nuc_has_dom.get(n, False)]
        return with_dom[0] if len(with_dom) == 1 else None
    return top_nuc


def _f447_proto_group_by_adjacency(matches, fwhm_at_channel,
                                   overlap_threshold_fwhm):
    """F-447 — sequential-adjacency прото-группы для ROI-owner анализа.

    matches сортируются по channel; соседи в пределах
    overlap_threshold_fwhm·FWHM_avg объединяются в одну прото-ROI.
    Грубее чем union-find, но достаточно для определения владельца ДО
    финального grouping в find_multiplet_regions.
    """
    if not matches:
        return []
    sm = sorted(matches, key=lambda m: m.peak_channel)
    groups = [[sm[0]]]
    for m in sm[1:]:
        prev = groups[-1][-1]
        try:
            f_avg = 0.5 * (
                float(fwhm_at_channel(prev.peak_channel))
                + float(fwhm_at_channel(m.peak_channel))
            )
        except Exception:
            f_avg = 0.0
        if (m.peak_channel - prev.peak_channel) <= (
                overlap_threshold_fwhm * f_avg):
            groups[-1].append(m)
        else:
            groups.append([m])
    return groups


# ============================================================================
# F-441 helpers - isolated-peak classifier (Rayleigh isolation)
# ============================================================================

def _f441_flatten_library_lines(nuclide_library, detected_nuclides):
    """F-441 - flatten library to [(nuclide, E_keV, I_pct), ...] for
    ALL confirmed detected nuclides.

    Used by ``_is_isolated_peak`` as the candidate-neighbour pool. Per
    the brief, the pool MUST include lines from all detected nuclides
    (not just the line's own nuclide) so that e.g. Ac-228 weak lines
    near a Tl-208 candidate are correctly considered as potential
    blockers of isolation.

    nuclide_library entries follow data/nuclides.json layout
    (``{"chain": ..., "lines": [(E, I), ...]}``). Both list/tuple and
    dict line records are accepted.
    """
    out = []
    if not nuclide_library or not detected_nuclides:
        return out
    for nuc in detected_nuclides:
        rec = nuclide_library.get(nuc, {}) or {}
        lines = rec.get("lines") or []
        for line in lines:
            if isinstance(line, (list, tuple)) and len(line) >= 2:
                E_lib = float(line[0])
                I_pct = float(line[1])
            elif isinstance(line, dict):
                E_lib = float(line.get("E_keV") or 0.0)
                I_pct = float(line.get("I_pct") or 0.0)
            else:
                continue
            if E_lib <= 0.0:
                continue
            out.append((nuc, E_lib, I_pct))
    return out


def _f441_fwhm_keV_at_channel(ch, fwhm_at_channel, spec):
    """F-441 - convert FWHM(channels) at channel ``ch`` to FWHM in keV.

    Uses the same pattern as F-381 phantom-sigma derivation
    (this file, ~lines 1227-1240, v1.21.0): evaluate
    ``spec.channel_to_energy`` at ch +- FWHM/2 channels and take the
    absolute difference. Returns 0.0 on any failure - caller treats
    0.0 as a degenerate case (no window -> classifier abstains).
    """
    if spec is None:
        return 0.0
    try:
        fwhm_ch = float(fwhm_at_channel(ch))
    except Exception:
        return 0.0
    if fwhm_ch <= 0.0:
        return 0.0
    try:
        half = fwhm_ch / 2.0
        E_lo = float(spec.channel_to_energy(float(ch) - half))
        E_hi = float(spec.channel_to_energy(float(ch) + half))
    except Exception:
        return 0.0
    return abs(E_hi - E_lo)


def _is_isolated_peak(
    m,
    fwhm_at_channel,
    spec,
    all_lib_lines,
    *,
    window_fwhm: float = 1.0,
    min_neighbor_I_pct: float = 3.0,
    self_match_keV: float = 0.1,
) -> bool:
    """F-441 - return True iff ``m`` has no strong library-line neighbour
    within ``+/- window_fwhm * FWHM_keV(E_m)``.

    A "strong neighbour" is any library line in the window (from the
    pre-flattened ``all_lib_lines`` pool covering all confirmed
    detected nuclides) with ``I_gamma_pct >= min_neighbor_I_pct``
    (default 3.0). The line itself is excluded by an
    ``|dE| < self_match_keV`` skip (default 0.1 keV) - this also
    drops trivial duplicates that may arise when the same nuclide
    appears twice in the detection list.

    Defaults (1.0 FWHM window, 3.0 I_pct threshold) calibrated on three
    reference NaI lines per F-441 brief table:

      * Tl-208 583.19 keV (FWHM ~47 keV on NaI 63x63) - strongest
        neighbour Ac-228 562.50 at 0.87% --> isolated.
      * K-40 1461 keV (FWHM ~64 keV) - no nuclide library neighbour
        in window --> isolated.
      * Tl-208 2614.51 keV (FWHM ~95 keV) - no nuclide library
        neighbour in window --> isolated (already routed to cowell
        in pre-F-441 baseline).

    When FWHM cannot be evaluated (degenerate spec / fwhm_at_channel),
    the classifier returns False (conservative - keep pre-F-441 routing
    so the change is no-op on degenerate calibration).

    Reference: Lord Rayleigh (1879), classical resolution criterion;
    LSRM Algorithmic Foundations 2025 Section 3.3 (peak grouping
    criteria); Gilmore & Joss 3rd Ed. Section 6.4 ("Peak overlap and
    the resolution criterion") - ``delta E >= FWHM`` as the boundary
    between isolated and multiplet behaviour on scintillators.
    """
    if not all_lib_lines:
        return True
    try:
        E_m = float(m.library_E_keV)
    except Exception:
        return True
    if E_m <= 0.0:
        return True
    fwhm_keV = _f441_fwhm_keV_at_channel(
        m.peak_channel, fwhm_at_channel, spec,
    )
    if fwhm_keV <= 0.0:
        # Conservative: pre-F-441 behaviour (let union-find decide).
        return False
    win = window_fwhm * fwhm_keV
    lo = E_m - win
    hi = E_m + win
    for _nuc, E_l, I_l in all_lib_lines:
        if not (lo <= E_l <= hi):
            continue
        if abs(E_l - E_m) < self_match_keV:
            continue
        if I_l >= min_neighbor_I_pct:
            return False
    return True


# ============================================================================
# Multiplet detection from identification results
# ============================================================================

def _split_zones_lzmax(
    clusters: list,
    spec,
    fwhm_at_channel: Callable[[float], float],
    *,
    max_zone_length_fwhm: float = 10.0,
    roi_window_factor: float = 2.5,
    _max_depth: int = 64,
) -> list:
    """Step 3 (LSRM Lzmax) — верхний предел длины зоны.

    Дробит чрезмерно длинную зону (участок под ОДНИМ полиномом фона) в точке
    минимума отсчётов между крайними пиками, чтобы снизить влияние не входящих
    в зону соседних пиков. Источник — Lsrm_algorithmic_foundations.pdf.md:
    482-485 (зона = совместно обрабатываемый участок, Lzmax — config-параметр)
    и 490-494 (при превышении длины — split «согласно критерию минимума по числу
    отсчётов»). Обоснование выбора определения длины зоны — CLAUDE.md, раздел
    «Ширина зоны / зонирование спектра».

    Длина зоны := размах пиков + 2·roi_window_factor·ПШПВ(центр) — полный участок
    под одним полиномом, т.е. интеграционная ROI ±roi_window_factor ПШПВ с каждого
    края (синхрон со Step 2). ПШПВ берётся в центральном канале зоны. Порог =
    max_zone_length_fwhm·ПШПВ(центр); дефолт 10 ПШПВ = значение Гамма-1С UI
    «макс. длина зоны 10 ПШПВ». Сравнение строгое («>»), разбиение рекурсивное.

    Graceful no-op (точное pre-Step-3 поведение): max_zone_length_fwhm <= 0.0,
    spec is None, пустой clusters или нечитаемые counts — вход возвращается без
    изменений.
    """
    if max_zone_length_fwhm <= 0.0 or spec is None or not clusters:
        return clusters
    try:
        counts = np.asarray(spec.counts, dtype=np.float64)
        n_ch = int(len(counts))
    except Exception:
        return clusters
    if n_ch <= 0:
        return clusters

    def _split_one(members: list, depth: int) -> list:
        if len(members) < 2:
            return [members]
        chans = sorted(float(m.peak_channel) for m in members)
        ch_lo, ch_hi = chans[0], chans[-1]
        center = 0.5 * (ch_lo + ch_hi)
        try:
            fwhm_c = float(fwhm_at_channel(center))
        except Exception:
            return [members]
        if fwhm_c <= 0.0:
            return [members]
        wing = roi_window_factor * fwhm_c
        zone_len = (ch_hi - ch_lo) + 2.0 * wing
        if zone_len <= max_zone_length_fwhm * fwhm_c:
            return [members]
        if depth <= 0:
            # Наблюдаемость: исчерпан лимит рекурсии при остающейся длинной
            # зоне. Без warn'а зона молча возвращалась бы как есть — это
            # ловится только косвенно (через размер кластера), что и стало
            # источником false-green для Step 3. Censor 2026-06-21.
            warnings.warn(
                f"_split_zones_lzmax: исчерпан _max_depth при "
                f"{len(members)} компонентах в зоне "
                f"[ch≈{ch_lo:.1f}..{ch_hi:.1f}, len/FWHM≈"
                f"{zone_len / fwhm_c:.2f}, порог={max_zone_length_fwhm}]; "
                f"зона возвращена без дальнейшего split.",
                RuntimeWarning,
                stacklevel=3,
            )
            return [members]
        # Зона длиннее Lzmax → режем в самой глубокой впадине СТРОГО между
        # крайними пиками (минимум отсчётов; pdf.md:492-493).
        lo_i = max(0, int(math.floor(ch_lo)) + 1)
        hi_i = min(n_ch - 1, int(math.ceil(ch_hi)) - 1)
        if hi_i <= lo_i:
            return [members]
        seg = counts[lo_i:hi_i + 1]
        if seg.size == 0:
            return [members]
        split_ch = lo_i + int(np.argmin(seg))
        left = [m for m in members if float(m.peak_channel) <= split_ch]
        right = [m for m in members if float(m.peak_channel) > split_ch]
        if not left or not right:
            return [members]
        return _split_one(left, depth - 1) + _split_one(right, depth - 1)

    out: list = []
    for cl in clusters:
        if not cl:
            out.append(cl)
            continue
        out.extend(_split_one(list(cl), _max_depth))
    out.sort(key=lambda g: g[0].peak_channel if g else 0.0)
    return out


def find_multiplet_regions(
    identification_result,
    fwhm_at_channel: Callable[[float], float],
    *,
    overlap_threshold_fwhm: float = 1.0,
    expand_to_display_window: bool = True,
    display_window_fwhm: float = 3.0,
    # F-381 / v1.18.25.2 — после expand_to_display_window, втягивать
    # ВСЕ library lines уже-detected нуклидов, попадающие в ROI кластера.
    # Без этого мультиплет не моделирует сильные линии нуклида, которые
    # не были matched как LineMatch (например, Tl-208 583 при
    # disambiguation gap'е). Требует spec для энергия→канал и
    # nuclide_library для library line lookup.
    spec=None,
    nuclide_library: Optional[dict] = None,
    min_library_intensity_pct: float = 0.5,
    # F-387 / v1.18.26 → F-387.1 / v1.18.26.1 — финальный фильтр
    # «неразрешённый мультиплет» по критерию Рэлея.
    #
    # Семантика **изменена** в v1.18.26.1:
    #   v1.18.26 (F-387):  pair (a,b) unresolved ⟺ |ΔE| < factor · FWHM_min
    #                       (factor=0.7, один глобальный min по cluster'у —
    #                        keep cluster as monolith if ANY pair unresolved)
    #   v1.18.26.1 (F-387.1): pair (a,b) unresolved ⟺
    #                       |ΔE| < factor · FWHM_avg(a, b)
    #                       где FWHM_avg(a, b) = (FWHM_a + FWHM_b) / 2.
    #                       factor=1.0 = классический критерий Рэлея:
    #                       пики резолвимы ⟺ ΔE ≥ FWHM_avg
    #                       (эквивалентно «10%-критерию полного разделения»
    #                        ΔE ≥ (FWHM_a + FWHM_b)/2 и впадине ~20%
    #                        от вершины пика).
    #
    # F-387.1 также заменяет «keep cluster as monolith» на:
    #   (a) построить граф вершин = components с edge (a, b) ⟺ unresolved;
    #   (b) разбить на connected components (BFS);
    #   (c) CC размера ≥ 2 = legitimate unresolved sub-cluster
    #       (multiplet), CC размера 1 = isolated singleton
    #       (cluster size=1 → single-component "multiplet" fit
    #        в deconvolve_multiplet, см. этого файла docstring);
    #   (d) для sub-cluster'а размера > max_components_per_cluster (default 3,
    #       по правилу пользователя «разрешать 2-3 выраженных в составе
    #       мультиплета»): top-K=3 по library_I_pct оставляем активными,
    #       остальные помечаются `peak_area_source="library_anchor_phantom"`
    #       (НЕ fit'ятся как отдельные Gaussians, но сохраняются в
    #        cluster.components для evidence/identification).
    #
    # Default factor 1.0 (Rayleigh). Старый F-387 default 0.7 — другая
    # семантика (FWHM_min, не FWHM_avg), поэтому численная разница
    # между 1.0·FWHM_avg и 0.7·FWHM_min не катастрофична на типичных
    # NaI clusters: для FWHM_a≈FWHM_b пары 0.7·FWHM ≈ 0.7·FWHM_avg vs
    # 1.0·FWHM_avg — Rayleigh пропускает чуть больше пар как unresolved
    # (что физически правильнее).
    #
    # References:
    #   - Lord Rayleigh (1879), "Investigations in optics" — критерий
    #     резолвимости двух точечных источников; для двух Гауссовых
    #     пиков эквивалентен ΔE ≥ FWHM (~20% впадина между пиками).
    #   - Lsrm Алгоритмические основы 2022, §9.4 — unresolved-multiplet
    #     fit при ΔE < FWHM.
    #   - Gilmore & Joss 3rd Ed., §6.4 «Peak overlap and the resolution
    #     critereon», стр. 158 — `δE = FWHM` как граница резолвимости
    #     для счётчика NaI.
    #
    # **NaI 63×63 calibration override (v1.18.26.1)**: default factor=1.1
    # (не 1.0). Эмпирически на Th-232 demo при factor=1.0 marginal pairs
    # (Ac-228 463 vs 504 factor=1.03; 911 vs 969 factor=1.01) попадают
    # в "resolved" зону по строгому Rayleigh, хотя физически на NaI с
    # валли-фил от интенсивных соседей они unresolved (Ac-228 911+969
    # не разделяются практически из-за weak peaks 860/965 между ними,
    # которые тонут на фоне доминантных). 10%-расширение порога к 1.1
    # покрывает borderline кейсы без излишнего слипания резолвимых.
    unresolved_separation_fwhm_factor: float = 1.1,
    # F-387.1 — top-K cap: после CC-split, sub-cluster крупнее max_K
    # урезается до top_K активных по library_I_pct; остальные
    # сохраняются как phantom anchors (см. выше).
    max_components_per_cluster: int = 3,
    # F-391 / v1.18.27 — S/N significance gate для multiplet members.
    # Per Gilmore & Joss 3rd Ed. §5.5, пик считается «реально измеренным»
    # при S/N ≥ 3 (Currie L_C, k=3). Компоненты с S/N ниже порога
    # помечаются как phantom anchors (peak_area_source=
    # "library_anchor_phantom") — НЕ участвуют в Rayleigh-CC build для
    # determining topology, НЕ становятся отдельными free-параметрами
    # в coupled fit, НО сохраняются в cluster для evidence/JSON.
    #
    # Default 3.0 = классический 3σ FEP-significance criterion. Гасит:
    #   - F-381 library-anchor lines, добавленные enrichment'ом, у
    #     которых нет measured peak (S/N=0);
    #   - matched lines с слабым signal (Ac-228 409 keV I=1.74% при
    #     отсутствии видимого пика — peak_area_uncertainty велика → S/N<3).
    #
    # 0.0 отключает gate (back-compat / diagnostic mode).
    min_significance_snr: float = 3.0,
    # BUG-3 / 2026-06-02 — strong-anchor library coverage (Fix #1).
    # OPT-IN: для каждого кластера запрашиваем у nuclide_library все линии
    # обнаруженных нуклидов И их chain-mates в окне [E_min, E_max] и
    # добавляем строкой `peak_area_source="library_anchor_strong"`, если
    # I_pct ≥ enable_strong_anchor_enrichment·I_max_in_window.
    #
    # 0.0 = ВЫКЛ. 0.05 = brief'овский порог 5% от I_max в окне (default).
    #
    # BUG-22 / 2026-06-02 — default flipped 0.0 → 0.05 to activate the
    # BUG-3 Fix #1 strong-anchor enrichment by default. The earlier
    # rationale for keeping it at 0.0 was a misdocumentation of the
    # LSRM FWHM polynomial domain (see lsrm_spe.py header): the
    # `lsrm_fwhm_polynomial_in_E` model is in fact a polynomial in
    # z = √E_keV (LSRM Algorithmic Foundations §8.3), and the
    # `validate_certs.py::make_lsrm_fwhm_provider` correctly evaluates
    # it on this domain. Once that domain is honored, FWHM at 238 keV
    # is ≈23 keV (not 58 keV from the broken polynomial-in-E reading),
    # and strong-anchor enrichment fits converge as designed on the
    # Th-232 demo fixture.
    enable_strong_anchor_enrichment: float = 0.05,
    # F-440 / v1.30.0 — Two-phase weak-line completion (Phase 1 gate).
    # Brief: _state/agent_a/inbox/2026-06-13_F-440_two_phase_weak_line_completion.md
    #
    # Линия включается в multiplet topology (grouping → Rayleigh-CC →
    # top-K cap → fit) только если ОБА:
    #   • S/N (peak_area / peak_area_uncertainty или significance_currie) ≥
    #     min_grouping_snr (default 5.0)
    #   • library_I_pct ≥ min_grouping_intensity_pct (default 3.0)
    #
    # Фильтр применяется в самом начале функции, ДО union-find. Слабые
    # matches перемечаются peak_area_source=library_anchor_phantom: они
    # автоматически вытесняются из active Rayleigh-CC graph и из top-K cap
    # (через _is_phantom-проверку), сохраняются в evidence для Phase 2
    # completion (scripts/gamma/activity/weak_line_completion.py), но НЕ
    # становятся free-параметрами fit'а.
    #
    # 0.0 для любого параметра отключает соответствующий gate
    # (back-compat / diagnostic mode). Для отключения F-440 целиком
    # передать min_grouping_snr=0.0 И min_grouping_intensity_pct=0.0.
    #
    # Reference: LSRM Algorithmic Foundations 2025 §6.2 «fit only what
    # you measure»; Gilmore & Joss 3rd Ed. §9.6.4 «library correction
    # afterwards».
    # F-440 / v1.30.0 — Phase 1 weak-line gate (S/N + I_pct). Default
    # OFF (0.0/0.0) для сохранения baseline behaviour. Включается
    # передачей 5.0/3.0 (бриф F-440 spec). Phase 2 weak-line completion
    # выполняется независимо (всегда), как информационный JSON-блок.
    min_grouping_snr: float = 0.0,
    min_grouping_intensity_pct: float = 0.0,
    # F-447 — guest-only ROI-owner aware фильтр. False (default) =
    # F-440 глобальный фильтр (как было). True = фильтруются ТОЛЬКО
    # guest-линии чужих нуклидов в ROI; library-линии нуклида-владельца
    # ROI сохраняются безусловно. Цель — не дропать phantom anchor'ы
    # соседей, которые удерживают continuum от reabsorption (см.
    # _f447_identify_roi_owner / _f447_proto_group_by_adjacency).
    f440_guest_only_filter: bool = False,
    # F-441 / v1.31.4 - isolated-peak classifier (Rayleigh isolation).
    # Brief: _state/agent_a/inbox/2026-06-13_F-441_isolated_peak_classifier.md
    #
    # Lines whose library energy has NO library-line neighbour (across
    # all confirmed detected nuclides) within +- isolated_window_fwhm *
    # FWHM_keV(E_m) at intensity >= isolated_min_neighbor_intensity_pct
    # are classified as "isolated" and removed from ``all_matches``
    # BEFORE union-find grouping. They then flow through the pipeline
    # unchanged (no cluster -> apply_multiplet_deconvolution leaves
    # peak_area_source untouched -> stays as the original cowell /
    # lsrm_peaks_table area from identify.py).
    #
    # Defaults 1.0 / 3.0 are calibrated on three reference NaI lines:
    #
    #   * Tl-208 583.19 keV (strongest neighbour Ac-228 562 @ 0.87%
    #     << 3%) -> isolated -> cowell;
    #   * K-40 1460.82 keV (no neighbour in window) -> isolated -> cowell;
    #   * Tl-208 2614.51 keV (no neighbour in window) -> isolated
    #     (already so in pre-F-441 baseline).
    #
    # Negative controls (must remain multiplets):
    #   * Ac-228 233/252/277 cluster M3 - three strong lines within
    #     +-1 FWHM, multiple I_pct >= 3% blockers -> multiplet.
    #   * Ac-228 911 / 969 region M4 - two strong lines within +-1 FWHM,
    #     each blocks the other -> multiplet.
    #
    # 0.0 for window OR 0.0 for intensity disables the classifier
    # (back-compat / diagnostic mode). spec is None or
    # nuclide_library is None also disables (same gating as F-381).
    #
    # Reference: Lord Rayleigh (1879), classical resolution criterion;
    # LSRM Algorithmic Foundations 2025 Section 3.3 (peak grouping
    # criteria); Gilmore & Joss 3rd Ed. Section 6.4 ("Peak overlap and
    # the resolution criterion") - delta E >= FWHM boundary on NaI.
    isolated_window_fwhm: float = 1.0,
    isolated_min_neighbor_intensity_pct: float = 3.0,
    # Step 3 (LSRM Lzmax) — zone-length cap (pdf.md:482-494). A
    # cluster whose one-polynomial informative segment (peak-span plus
    # ±roi_window_factor·FWHM baseline wings = integration ROI ±2.5 ПШПВ,
    # Step-2 synced) exceeds max_zone_length_fwhm·FWHM is split at the
    # interior minimum-counts channel into shorter zones. Гамма-1С UI
    # «макс. длина зоны 10 ПШПВ» → default 10.0. 0.0 disables (exact
    # pre-Step-3 behaviour); spec=None also disables. FWHM at zone center.
    max_zone_length_fwhm: float = 10.0,
    lzmax_roi_window_factor: float = 2.5,
) -> list:
    """
    Group `LineMatch` objects (across all detected nuclides) into
    multiplet candidate clusters.

    A cluster is a maximal set of `LineMatch` entries such that every
    entry is within `overlap_threshold_fwhm × FWHM` (in channel space)
    of at least one other entry in the same set. FWHM is evaluated at
    each line's expected channel via `fwhm_at_channel`. Clusters of
    size 1 are dropped — they're isolated peaks, not multiplets.

    The grouping uses transitive closure (single-linkage clustering),
    so three adjacent overlapping lines all land in the same cluster
    even if the outermost pair lies just outside the threshold.

    F-374 — when ``expand_to_display_window=True`` (default), each
    cluster's ROI is extended by ``display_window_fwhm × FWHM`` on
    each side (the typical chart display window) and ANY identified
    LineMatch falling inside this expanded ROI is added to the cluster,
    even if it didn't satisfy the initial overlap criterion. This
    ensures the deconvolution accounts for all peaks visible in the
    chart's energy range — preventing the renderer from showing a
    "fit" that misses the dominant lines and produces large closure
    residuals. Pre-F-374 behavior is preserved by passing
    ``expand_to_display_window=False``.

    Args:
        identification_result: `IdentificationResult` from
            `gamma.identification.identify.identify_nuclides`.
        fwhm_at_channel: callable(channel) → FWHM in channels.
        overlap_threshold_fwhm: pairs within this much × FWHM are
            considered overlapping. 1.0 = "the two FWHM circles
            touch" — a sensible default for "the peaks blur into
            each other on the detector".
        expand_to_display_window: when True (default), expand each
            cluster ROI to the chart display window and absorb any
            identified peaks inside that window.
        display_window_fwhm: half-width of the display-window
            expansion, in units of cluster-edge FWHM. Default 3.0
            (matches the typical ±3·FWHM chart axis range).
        unresolved_separation_fwhm_factor: F-387 → F-387.1. Per-pair
            Rayleigh criterion. Pair (a, b) is **unresolved** ⟺
            |Δch| < factor · FWHM_avg(a, b), где
            FWHM_avg = (FWHM_a + FWHM_b) / 2. Default 1.0 = Rayleigh:
            пики резолвимы ⟺ ΔE ≥ FWHM_avg. **Семантика изменена в
            v1.18.26.1**: было factor·FWHM_min (один глобальный
            minimum), стало factor·FWHM_avg (среднее пары) — это
            классический Rayleigh, а не приближение. factor=0.0
            отключает фильтр (back-compat / diagnostic mode).
        max_components_per_cluster: F-387.1 top-K cap. После CC-split,
            sub-cluster крупнее max_K → top-K по library_I_pct активны
            (fit'ятся), остальные → phantom anchors с
            `peak_area_source="library_anchor_phantom"`. Phantom-ы
            сохраняются в cluster.components для evidence, но НЕ
            создают отдельных свободных параметров в fit'е. Default 3
            (по правилу «разрешать 2-3 выраженных в составе мультиплета»).

    Returns:
        List of clusters. F-387.1: каждый исходный cluster может быть
        разбит на N sub-cluster'ов (CC-split). Cluster — это
        `list[LineMatch]`, отсортированный по `peak_channel`. Sub-cluster
        размера 1 (isolated singleton) разрешается на выходе — он
        приведёт к 1-component "multiplet" fit downstream (см.
        `deconvolve_multiplet` docstring: «Single-component multiplet
        is accepted and reduces to a regular single-Gaussian fit»).
        Список clusters отсортирован по первому channel внутри cluster.
    """
    all_matches = []
    for ni in identification_result.detected_nuclides:
        all_matches.extend(ni.matched_lines)
    if not all_matches:
        return []

    # F-440 / v1.30.0 — Phase 1 weak-line filter. Линии с I_pct <
    # min_grouping_intensity_pct ИЛИ S/N < min_grouping_snr исключаются
    # из grouping/Rayleigh-CC/topology, НО не модифицируются (LineMatch
    # immutable, оригинальные peak_area/source сохраняются для step 9 /
    # activity calc / Phase 2 weak-line completion).
    # Existing library_anchor / library_anchor_phantom — propagate as-is.
    if (min_grouping_snr > 0.0 or min_grouping_intensity_pct > 0.0):
        if f440_guest_only_filter:
            # F-447 — guest-only: owner лiний нуклида-владельца ROI
            # сохраняем безусловно; threshold-фильтр применяется ТОЛЬКО
            # к «гостевым» линиям чужих нуклидов в той же ROI. Мотивация:
            # глобальный фильтр Phase 1 (F-440 5.0/3.0) дропает phantom
            # anchor'ы соседей (Ac-228 weak lines в M5 Th-232 demo),
            # которые удерживают continuum от reabsorption Tl-208 583 →
            # area drops → F-91 §7 switches sigma method → activity
            # катастрофа. Guest-only сохраняет library_I_pct mass
            # внутри ROI владельца, дропая только реально лишние гостевые.
            proto_groups = _f447_proto_group_by_adjacency(
                all_matches, fwhm_at_channel, overlap_threshold_fwhm,
            )
            f440_kept = []
            for proto in proto_groups:
                owner = _f447_identify_roi_owner(proto)
                for m in proto:
                    src_m = str(getattr(m, "peak_area_source", "") or "")
                    if src_m in ("library_anchor", "library_anchor_phantom"):
                        f440_kept.append(m)
                        continue
                    if owner is not None and getattr(m, "nuclide", None) == owner:
                        f440_kept.append(m)
                        continue
                    I_pct = float(getattr(m, "library_I_pct", 0.0) or 0.0)
                    snr = _f391_peak_snr(m)
                    fails_I = (min_grouping_intensity_pct > 0.0
                               and I_pct < min_grouping_intensity_pct)
                    fails_snr = (min_grouping_snr > 0.0
                                 and snr < min_grouping_snr)
                    if not (fails_I or fails_snr):
                        f440_kept.append(m)
                    else:
                        # F-447 V2 — failed guest НЕ дропается, а помечается
                        # phantom anchor: сохраняет ROI topology / continuum
                        # convergence, но не получает свободный параметр в
                        # coupled fit. Drop сужает ROI → continuum reabsorbs
                        # owner area (Tl-208 catastrophe demo 2026-06-15).
                        f440_kept.append(_f391_mark_phantom(m))
            all_matches = f440_kept
        else:
            f440_kept = []
            for m in all_matches:
                src_m = str(getattr(m, "peak_area_source", "") or "")
                if src_m in ("library_anchor", "library_anchor_phantom"):
                    f440_kept.append(m)
                    continue
                I_pct = float(getattr(m, "library_I_pct", 0.0) or 0.0)
                snr = _f391_peak_snr(m)
                fails_I = (min_grouping_intensity_pct > 0.0
                           and I_pct < min_grouping_intensity_pct)
                fails_snr = (min_grouping_snr > 0.0
                             and snr < min_grouping_snr)
                if not (fails_I or fails_snr):
                    f440_kept.append(m)
                # else: weak — drop from grouping. LineMatch lives on in
                # identification.matched_lines (caller). Phase 2 picks it up.
            all_matches = f440_kept

    # F-441 / v1.31.4 - isolated-peak classifier. Lines with no strong
    # library-line neighbour within +/- isolated_window_fwhm * FWHM_keV
    # are removed from all_matches BEFORE union-find. They then flow
    # through identification_result.matched_lines unchanged - downstream
    # apply_multiplet_deconvolution leaves their peak_area_source intact
    # (typically "cowell" / "lsrm_peaks_table" from identify.py).
    #
    # Skip when classifier is mis-configured (window or intensity = 0)
    # OR when spec/library is unavailable (same gating as F-381).
    # library_anchor / library_anchor_phantom lines are NEVER tested -
    # they were not in all_matches at this point anyway (they are added
    # later via F-381 inside-cluster enrichment), but defensive guard.
    if (isolated_window_fwhm > 0.0
            and isolated_min_neighbor_intensity_pct > 0.0
            and spec is not None
            and nuclide_library is not None
            and all_matches):
        try:
            detected_nuc_names = [
                ni.nuclide for ni in identification_result.detected_nuclides
            ]
            f441_lib_pool = _f441_flatten_library_lines(
                nuclide_library, detected_nuc_names,
            )
            f441_kept = []
            for m in all_matches:
                src_m = str(getattr(m, "peak_area_source", "") or "")
                if src_m in ("library_anchor", "library_anchor_phantom"):
                    # Defensive: phantoms should not be in all_matches
                    # here, but keep for safety - they belong to F-381
                    # cluster enrichment, not isolation classification.
                    f441_kept.append(m)
                    continue
                if _is_isolated_peak(
                    m, fwhm_at_channel, spec, f441_lib_pool,
                    window_fwhm=isolated_window_fwhm,
                    min_neighbor_I_pct=isolated_min_neighbor_intensity_pct,
                ):
                    # Isolated - drop from all_matches so union-find never
                    # groups it. LineMatch survives untouched in
                    # identification_result.matched_lines for downstream.
                    continue
                f441_kept.append(m)
            all_matches = f441_kept
        except Exception:
            # Graceful: F-441 classifier failure must not break the
            # pipeline. Fall back to pre-F-441 grouping (all_matches
            # untouched).
            pass

    # Sort by channel for predictable output and faster grouping
    all_matches.sort(key=lambda m: m.peak_channel)

    # Union-find over indices
    parent = list(range(len(all_matches)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)

    # For every pair, check overlap. O(N²) is fine — N is small after
    # identification (rarely more than a few dozen matched lines).
    for i in range(len(all_matches)):
        ch_i = all_matches[i].peak_channel
        fwhm_i = float(fwhm_at_channel(ch_i))
        for j in range(i + 1, len(all_matches)):
            ch_j = all_matches[j].peak_channel
            fwhm_j = float(fwhm_at_channel(ch_j))
            threshold = overlap_threshold_fwhm * 0.5 * (fwhm_i + fwhm_j)
            if abs(ch_j - ch_i) <= threshold:
                union(i, j)
            else:
                # since matches are sorted by channel and FWHM is
                # bounded, anything further to the right won't reach
                # back to i either — early break
                if ch_j - ch_i > 2.0 * overlap_threshold_fwhm * max(fwhm_i, fwhm_j):
                    break

    groups: dict = {}
    for idx in range(len(all_matches)):
        root = find(idx)
        groups.setdefault(root, []).append(all_matches[idx])

    clusters = [g for g in groups.values() if len(g) >= 2]
    clusters.sort(key=lambda g: g[0].peak_channel)

    # F-374 — expand each cluster to its display window and absorb any
    # additional identified peaks that fall inside. This prevents the
    # deconvolution renderer from drawing a chart with strong peaks in
    # the spectrum that aren't accounted for by the fit (visible as a
    # "fit hugging the baseline while data has unmodelled peaks").
    if expand_to_display_window and clusters:
        # Sort each cluster internally first
        for cl in clusters:
            cl.sort(key=lambda m: m.peak_channel)
        # Build "in-some-cluster" set so we don't double-add
        clustered_keys = set()
        for cl in clusters:
            for m in cl:
                clustered_keys.add(id(m))
        # For every cluster, expand ROI and absorb any LineMatch inside
        for cl in clusters:
            ch_lo_cluster = min(m.peak_channel for m in cl)
            ch_hi_cluster = max(m.peak_channel for m in cl)
            fwhm_lo = float(fwhm_at_channel(ch_lo_cluster))
            fwhm_hi = float(fwhm_at_channel(ch_hi_cluster))
            ch_lo_disp = ch_lo_cluster - display_window_fwhm * fwhm_lo
            ch_hi_disp = ch_hi_cluster + display_window_fwhm * fwhm_hi
            for m in all_matches:
                if id(m) in clustered_keys:
                    continue
                if ch_lo_disp <= m.peak_channel <= ch_hi_disp:
                    cl.append(m)
                    clustered_keys.add(id(m))
            cl.sort(key=lambda m: m.peak_channel)

        # F-381 / v1.18.25.2 — pull ALL library lines of identified
        # nuclides (not only peak-matched LineMatch'es) that fall in
        # each cluster's expanded ROI. Without this:
        #   - Tl-208 583.187 кэВ (I=30.5%) если pipeline её не сматчил
        #     (disambiguation gap) → не моделируется → χ²↑;
        #   - Ac-228 463 кэВ (I=4.4%) если sub-threshold для S/σ → не
        #     попадает в кластер → fit видит unexplained пик.
        # Создаём "library-anchor" LineMatch с peak_area=None (downstream
        # узнает по peak_area_source="library_anchor"). nuclide_library
        # имеет формат data/nuclides.json (см. _DEFAULT_LIB ниже).
        if spec is not None and nuclide_library is not None:
            try:
                from gamma.identification.identify import LineMatch
                # Список detected нуклидов
                detected_nucs = [
                    ni.nuclide for ni in
                    identification_result.detected_nuclides
                ]
                # Сет уже-присутствующих (nuclide, round(E_lib,1)) ключей
                existing_keys = set()
                for cl in clusters:
                    for m in cl:
                        existing_keys.add(
                            (m.nuclide, round(float(m.library_E_keV), 1))
                        )
                for cl in clusters:
                    # F-381 — фантомы добавляются в окно display, НО их
                    # peak_channel clamp'ится в pre-F-381 [ch_lo_cluster,
                    # ch_hi_cluster] bounds, чтобы НЕ расширять cluster
                    # ROI bounds (это сорвало бы _overlaps_forced
                    # фильтр в staged_pipeline: phantoms доходили бы до
                    # forced ROI и весь auto-cluster отбрасывался).
                    ch_lo_cluster = min(m.peak_channel for m in cl)
                    ch_hi_cluster = max(m.peak_channel for m in cl)
                    fwhm_lo = float(fwhm_at_channel(ch_lo_cluster))
                    fwhm_hi = float(fwhm_at_channel(ch_hi_cluster))
                    ch_lo_disp = ch_lo_cluster - display_window_fwhm * fwhm_lo
                    ch_hi_disp = ch_hi_cluster + display_window_fwhm * fwhm_hi
                    try:
                        E_lo_disp = float(spec.channel_to_energy(int(ch_lo_disp)))
                        E_hi_disp = float(spec.channel_to_energy(int(ch_hi_disp)))
                    except Exception:
                        continue
                    if E_lo_disp > E_hi_disp:
                        E_lo_disp, E_hi_disp = E_hi_disp, E_lo_disp
                    added_here = []
                    for nuc in detected_nucs:
                        lines = nuclide_library.get(nuc, {}).get("lines") or []
                        for line in lines:
                            if isinstance(line, (list, tuple)) and len(line) >= 2:
                                E_lib = float(line[0])
                                I_pct = float(line[1])
                            elif isinstance(line, dict):
                                E_lib = float(line.get("E_keV") or 0.0)
                                I_pct = float(line.get("I_pct") or 0.0)
                            else:
                                continue
                            if I_pct < min_library_intensity_pct:
                                continue
                            if not (E_lo_disp <= E_lib <= E_hi_disp):
                                continue
                            key = (nuc, round(E_lib, 1))
                            if key in existing_keys:
                                continue
                            try:
                                ch_real = int(round(spec.energy_to_channel(E_lib)))
                            except Exception:
                                continue
                            # F-381: clamp phantom peak_channel в bounds
                            # cluster'а чтобы не расширять ROI. Реальная
                            # E_lib сохранена в library_E_keV; downstream
                            # coupled_intensity_fit использует именно E_lib
                            # (не peak_channel) для построения Гауссиана.
                            ch_clamped = max(ch_lo_cluster,
                                             min(ch_hi_cluster, ch_real))
                            # BUG-34 W3 phantom normalisation (v1.21.0):
                            # populate gauss_sigma_keV from FWHM at this
                            # channel so chain-sibling-blend gate in
                            # compute.py:798 sees a correct σ instead of
                            # 0.0 → it can correctly accept/reject blend
                            # at this phantom anchor.
                            try:
                                _fwhm_ch_p = float(
                                    fwhm_at_channel(ch_clamped))
                                if _fwhm_ch_p > 0:
                                    _half = _fwhm_ch_p / 2.0
                                    _E_lo = float(spec.channel_to_energy(
                                        float(ch_clamped) - _half))
                                    _E_hi = float(spec.channel_to_energy(
                                        float(ch_clamped) + _half))
                                    _fwhm_keV_p = abs(_E_hi - _E_lo)
                                    _sig_keV_p = (
                                        _fwhm_keV_p / 2.354820045
                                        if _fwhm_keV_p > 0 else None
                                    )
                                else:
                                    _sig_keV_p = None
                            except Exception:
                                _sig_keV_p = None
                            phantom = LineMatch(
                                nuclide=nuc,
                                library_E_keV=E_lib,
                                library_I_pct=I_pct,
                                peak_channel=ch_clamped,
                                peak_E_keV=E_lib,
                                peak_sigma=0.0,
                                residual_keV=0.0,
                                is_characteristic=False,
                                peak_area=None,
                                peak_area_uncertainty=None,
                                peak_area_source="library_anchor",
                                gauss_sigma_keV=_sig_keV_p,
                            )
                            added_here.append(phantom)
                            existing_keys.add(key)
                    if added_here:
                        cl.extend(added_here)
                        cl.sort(key=lambda m: m.peak_channel)
            except ImportError:
                pass
            except Exception:
                # graceful: F-381 enrichment is non-critical, fallback
                # к F-374-only behaviour при любой ошибке
                pass

    # ─────────────────────────────────────────────────────────────────
    # BUG-3 (2026-06-02) — library-driven multiplet coverage
    # ─────────────────────────────────────────────────────────────────
    # Корень проблемы: F-381 enrichment добавляет library lines с
    # absolute I_pct ≥ min_library_intensity_pct (default 0.5%). Но они
    # помечаются как ``library_anchor`` (peak_area=None) → F-391 S/N gate
    # сразу демотит их в ``library_anchor_phantom`` → они НЕ участвуют
    # в fit как active free parameters. Counts strong-but-unmeasured
    # линий поглощаются в continuum-baseline, искажая её, а fit натягивает
    # амплитуды weak соседей.
    #
    # Th-232 M3 (162-310 кэВ) catastrophic case: Ac-228 209 (I=3.89%) и
    # Ac-228 270 (I=3.46%) — physically intense, видны в спектре, но не
    # detect'ились peak finder'ом (sit на склоне continuum). F-381 их
    # добавляет как library_anchor → phantom → fit видит только Tl-208 233
    # (I=0.11%) и Pb-212 238 (I=43.6%). χ²/ν=155, физически перевёрнутые
    # ratios (Tl-208 233 area > Pb-212 238 area при I-ratio 1:400).
    #
    # Fix: relative intensity threshold per cluster window.
    #   1. Для каждого cluster вычислить I_max — максимум library I_pct
    #      по всем lines всех (detected + chain) нуклидов в окне.
    #   2. Любая library line с I_pct ≥ 0.05·I_max получает
    #      ``peak_area_source="library_anchor_strong"`` — strong anchor.
    #   3. _f391_peak_snr возвращает inf для strong → пропускает S/N gate.
    #   4. apply_multiplet_deconvolution._is_phantom_lm НЕ считает strong
    #      phantom'ом → strong становится active fit component.
    #
    # Chain awareness: nuclide library имеет атрибут "chain"
    # (e.g. Ac-228, Pb-212, Tl-208, Bi-212 все chain="Th-232"). Если
    # detected содержит хоть одного chain-member, expand'им enrichment
    # pool до всех members этой цепочки. Это покрывает случай когда
    # detected ID получил только Tl-208 (по 2614 anchor) но в кластере
    # M3 нужны строгие линии Ac-228/Pb-212 того же chain.
    if (clusters and spec is not None and nuclide_library is not None
            and enable_strong_anchor_enrichment > 0.0):
        try:
            from gamma.identification.identify import LineMatch
            # Chain expansion (per-cluster): для каждого nuclide В КЛАСТЕРЕ
            # добавить всех chain mates. nuclide_library record format:
            #   {nuc_name: {"chain": "Th-232" | None, "lines": [...]}}
            # CRITICAL (regression fix 2026-06-02): chain_pool вычисляется
            # ОТДЕЛЬНО ДЛЯ КАЖДОГО кластера, а не глобально по всем
            # detected_nuclides. Если глобальный pool, то в singleton-
            # cluster Co-60 1173 будет добавлена Eu-152 1112 (I=13.4%) как
            # strong-anchor (обе nuclides detected, но chain=None у обеих,
            # и pool их объединял через `if ch_name in detected_chains` —
            # пустой set, но cluster_nucs UNION забирал все). Per-cluster
            # pool гарантирует: enrichment добавляет только «свои» линии
            # того же nuclide или того же chain'а что member кластера.

            STRONG_REL_THRESHOLD = float(enable_strong_anchor_enrichment)

            for cl in clusters:
                if not cl:
                    continue
                # Per-cluster chain pool — собирается ТОЛЬКО из ACTIVE
                # nuclides (с измеренным peak_area_source ∈ {cowell,
                # lsrm_peaks_table, deconvolved, ...}). F-381-added
                # phantom-anchors (library_anchor / library_anchor_phantom)
                # НЕ включаются в pool: иначе Eu-152 1112 (F-381-добавленная
                # в Co-60 1173 кластер) превратилась бы в strong anchor,
                # сломав F-387.2 singleton-routing для Co-60 doublet'a.
                # Контракт: chain enrichment продлевает coverage только
                # ФИЗИЧЕСКИ детектированных линий, не speculative anchors.
                _PHANTOM_SRCS = (
                    "library_anchor", "library_anchor_phantom",
                )
                cluster_nucs = {
                    m.nuclide for m in cl
                    if m.nuclide
                    and str(getattr(m, "peak_area_source", "") or "")
                    not in _PHANTOM_SRCS
                }
                chain_pool = set(cluster_nucs)
                cluster_chains = set()
                for nuc in cluster_nucs:
                    rec = nuclide_library.get(nuc, {}) or {}
                    ch_name = rec.get("chain")
                    if ch_name:
                        cluster_chains.add(ch_name)
                if cluster_chains:
                    for nuc_name, rec in nuclide_library.items():
                        if (rec or {}).get("chain") in cluster_chains:
                            chain_pool.add(nuc_name)

                # Cluster window in energy (use existing peak_channel
                # bounds expanded by display_window_fwhm, same as F-381).
                ch_lo_cluster = min(m.peak_channel for m in cl)
                ch_hi_cluster = max(m.peak_channel for m in cl)
                fwhm_lo = float(fwhm_at_channel(ch_lo_cluster))
                fwhm_hi = float(fwhm_at_channel(ch_hi_cluster))
                ch_lo_disp = ch_lo_cluster - display_window_fwhm * fwhm_lo
                ch_hi_disp = ch_hi_cluster + display_window_fwhm * fwhm_hi
                try:
                    E_lo_disp = float(spec.channel_to_energy(int(ch_lo_disp)))
                    E_hi_disp = float(spec.channel_to_energy(int(ch_hi_disp)))
                except Exception:
                    continue
                if E_lo_disp > E_hi_disp:
                    E_lo_disp, E_hi_disp = E_hi_disp, E_lo_disp

                # Collect all in-window candidate (nuclide, E_lib, I_pct).
                in_window = []
                for nuc in chain_pool:
                    lines = (nuclide_library.get(nuc, {}) or {}).get(
                        "lines") or []
                    for line in lines:
                        if (isinstance(line, (list, tuple))
                                and len(line) >= 2):
                            E_lib = float(line[0])
                            I_pct = float(line[1])
                        elif isinstance(line, dict):
                            E_lib = float(line.get("E_keV") or 0.0)
                            I_pct = float(line.get("I_pct") or 0.0)
                        else:
                            continue
                        if not (E_lo_disp <= E_lib <= E_hi_disp):
                            continue
                        if I_pct <= 0:
                            continue
                        in_window.append((nuc, E_lib, I_pct))

                if not in_window:
                    continue
                I_max = max(I for _, _, I in in_window)
                if I_max <= 0:
                    continue
                rel_threshold_pct = STRONG_REL_THRESHOLD * I_max

                # Existing keys in cluster (avoid duplicates).
                existing_keys_local = set(
                    (m.nuclide, round(float(m.library_E_keV), 1))
                    for m in cl
                )
                # We also want to UPGRADE already-present library_anchor
                # (F-381 added but phantom) to strong if it qualifies.
                # Build map nuclide+rounded_E → cluster index.
                upgrade_map: dict = {}
                for idx_m, m in enumerate(cl):
                    key_m = (m.nuclide, round(float(m.library_E_keV), 1))
                    src_m = str(getattr(m, "peak_area_source", "") or "")
                    if src_m == "library_anchor":
                        upgrade_map[key_m] = idx_m

                added_strong = []
                upgraded_indices: list = []
                for nuc, E_lib, I_pct in in_window:
                    if I_pct < rel_threshold_pct:
                        continue
                    key = (nuc, round(E_lib, 1))
                    if key in existing_keys_local:
                        # Try to upgrade if it's a library_anchor.
                        if key in upgrade_map:
                            upgraded_indices.append((upgrade_map[key], I_pct))
                        continue
                    try:
                        ch_real = int(round(
                            spec.energy_to_channel(E_lib)))
                    except Exception:
                        continue
                    ch_clamped = max(
                        ch_lo_cluster, min(ch_hi_cluster, ch_real))
                    # BUG-34 W3 phantom (strong-anchor) normalisation
                    # (v1.21.0): populate gauss_sigma_keV — same
                    # rationale as the library_anchor case above.
                    try:
                        _fwhm_ch_s = float(fwhm_at_channel(ch_clamped))
                        if _fwhm_ch_s > 0:
                            _half_s = _fwhm_ch_s / 2.0
                            _E_lo_s = float(spec.channel_to_energy(
                                float(ch_clamped) - _half_s))
                            _E_hi_s = float(spec.channel_to_energy(
                                float(ch_clamped) + _half_s))
                            _fwhm_keV_s = abs(_E_hi_s - _E_lo_s)
                            _sig_keV_s = (
                                _fwhm_keV_s / 2.354820045
                                if _fwhm_keV_s > 0 else None
                            )
                        else:
                            _sig_keV_s = None
                    except Exception:
                        _sig_keV_s = None
                    strong = LineMatch(
                        nuclide=nuc,
                        library_E_keV=E_lib,
                        library_I_pct=I_pct,
                        peak_channel=ch_clamped,
                        peak_E_keV=E_lib,
                        peak_sigma=0.0,
                        residual_keV=0.0,
                        is_characteristic=False,
                        peak_area=None,
                        peak_area_uncertainty=None,
                        peak_area_source="library_anchor_strong",
                        gauss_sigma_keV=_sig_keV_s,
                    )
                    added_strong.append(strong)
                    existing_keys_local.add(key)
                if upgraded_indices:
                    from dataclasses import replace as _dc_replace_lm2
                    for idx_m, _Ip in upgraded_indices:
                        try:
                            cl[idx_m] = _dc_replace_lm2(
                                cl[idx_m],
                                peak_area_source="library_anchor_strong",
                            )
                        except Exception:
                            try:
                                cl[idx_m].peak_area_source = (
                                    "library_anchor_strong")
                            except Exception:
                                pass
                if added_strong:
                    cl.extend(added_strong)
                    cl.sort(key=lambda m: m.peak_channel)
        except ImportError:
            pass
        except Exception:
            # graceful: BUG-3 enrichment не critical, fallback к pre-fix
            pass

    # F-391 / v1.18.27 — S/N significance gate (LSRM-9.4 / Gilmore §5.5).
    # Помечаем low-S/N компоненты как phantom anchors ДО Rayleigh-CC build,
    # чтобы они не влияли на topology (CC graph) и не «достраивали»
    # multiplet'ы из чистого фона. F-381 library_anchor (peak_area=None)
    # сразу получают S/N=0 → phantom.
    if clusters and min_significance_snr > 0.0:
        new_clusters_snr: list = []
        for cl in clusters:
            rebuilt: list = []
            for m in cl:
                # Уже-phantom (F-387.1 top-K cap или F-381 anchor с None)
                # — оставляем как есть.
                src = str(getattr(m, "peak_area_source", "") or "")
                if src == "library_anchor_phantom":
                    rebuilt.append(m)
                    continue
                snr = _f391_peak_snr(m)
                if snr < min_significance_snr:
                    rebuilt.append(_f391_mark_phantom(m))
                else:
                    rebuilt.append(m)
            new_clusters_snr.append(rebuilt)
        clusters = new_clusters_snr

    # F-387 / v1.18.26 → F-387.1 / v1.18.26.1 — Rayleigh-CC split +
    # top-K cap (см. полное обоснование в docstring и в kwarg-комментарии
    # выше). Алгоритм:
    #   1. Для каждого cluster построить граф: vertices = LineMatch'и,
    #      edge(a, b) ⟺ |Δch| < factor · FWHM_avg(a, b).
    #   2. BFS → connected components.
    #   3. CC размера ≥ 2 → unresolved multiplet sub-cluster.
    #      CC размера 1 → isolated singleton (cluster-size=1 на выходе,
    #      downstream запустит single-Gaussian fit через
    #      `deconvolve_multiplet`).
    #   4. Sub-cluster крупнее `max_components_per_cluster` (default 3) →
    #      top-K по library_I_pct активны; остальные помечаются
    #      `peak_area_source="library_anchor_phantom"` и остаются
    #      в cluster для evidence, но НЕ fit'ятся (фильтруются
    #      в apply_multiplet_deconvolution при сборке
    #      MultipletComponent).
    #
    # Hard-locked retained на real demo:
    #   - M1 Ac-228 911/964.8/969: 964.8+969 ≈4 кэВ при FWHM 30 кэВ
    #     (FWHM_avg=30) → 0.14·FWHM_avg < 1.0 → unresolved.
    #   - M2 Ac-228 1588 + Bi-212 1620 + Ac-228 1630: 1620+1630 ≈10 кэВ
    #     при FWHM 37 (FWHM_avg=37) → 0.27·FWHM_avg < 1.0 → unresolved
    #     → CC {1620, 1630}. 1588 — изолирован
    #     (Δ=32 кэВ при FWHM_avg ≈37 → 0.86·FWHM_avg < 1.0 → тоже
    #      unresolved → CC {1588, 1620, 1630}, transitive через 1620).
    if clusters and unresolved_separation_fwhm_factor > 0.0:
        try:
            from gamma.identification.identify import LineMatch  # noqa
            _have_linematch = True
        except Exception:
            _have_linematch = False
        from dataclasses import replace as _dc_replace_lm
        from collections import deque as _deque

        new_clusters: list = []
        for cl in clusters:
            n = len(cl)
            if n < 2:
                # Singleton input — пропускаем (overlap-этап выше уже
                # отбрасывает size-1; добрался только при странных
                # пограничных случаях, e.g. enrichment добавил, но
                # overlap не сматчил). Сохраняем как singleton sub-cluster.
                new_clusters.append(list(cl))
                continue

            # F-391 / v1.18.27 — phantom anchors (низкий S/N или F-381
            # library_anchor) НЕ участвуют в Rayleigh-CC topology build:
            # они не имеют реального signal и не должны «прокладывать»
            # рёбра между real-component'ами или формировать собственные
            # CC. Они остаются в cluster как evidence и распределяются
            # в финальные sub-cluster'ы согласно ближайшей real-component
            # (или собираются в trailing phantom-only cluster и
            # отбрасываются downstream).
            _is_phantom = lambda m: str(  # noqa: E731
                getattr(m, "peak_area_source", "") or ""
            ) in ("library_anchor", "library_anchor_phantom")
            active_idxs = [i for i in range(n) if not _is_phantom(cl[i])]
            phantom_idxs = [i for i in range(n) if _is_phantom(cl[i])]

            # Шаг 1: построить adjacency list через Rayleigh-pair edges
            # ТОЛЬКО среди active. F-387.1 NB: используем `library_E_keV`
            # (через spec→channel если spec доступен), а НЕ `peak_channel`
            # — потому что F-381 phantom anchors имеют peak_channel
            # **clamp'ed в cluster bounds**, что искажает Δch при
            # graph-CC. Реальная энергия библиотечной линии — единственный
            # консистентный source.
            def _ch_for(m):
                if spec is not None:
                    try:
                        return float(spec.energy_to_channel(
                            float(m.library_E_keV)
                        ))
                    except Exception:
                        pass
                return float(m.peak_channel)
            channels = [_ch_for(m) for m in cl]
            fwhms = [float(fwhm_at_channel(c)) for c in channels]
            adj: list = [[] for _ in range(n)]
            # F-391 — edges только между active нодами
            for ai in range(len(active_idxs)):
                i = active_idxs[ai]
                for aj in range(ai + 1, len(active_idxs)):
                    j = active_idxs[aj]
                    fwhm_avg = 0.5 * (fwhms[i] + fwhms[j])
                    if fwhm_avg <= 0:
                        continue
                    if abs(channels[i] - channels[j]) < (
                        unresolved_separation_fwhm_factor * fwhm_avg
                    ):
                        adj[i].append(j)
                        adj[j].append(i)

            # Шаг 2: BFS → connected components ТОЛЬКО по active.
            # Phantoms не имеют edges → каждый бы стал собственным CC=1.
            # Вместо этого phantom'ы прикрепляются к ближайшему active CC
            # по library_E_keV (или, если active'ов нет, формируют один
            # phantom-only fallback CC, который downstream дропнет).
            visited = [False] * n
            components_idxs: list = []
            for start in active_idxs:
                if visited[start]:
                    continue
                queue = _deque([start])
                visited[start] = True
                cc: list = []
                while queue:
                    u = queue.popleft()
                    cc.append(u)
                    for v in adj[u]:
                        if not visited[v]:
                            visited[v] = True
                            queue.append(v)
                components_idxs.append(cc)
            # F-391 — присоединить каждый phantom к ближайшему CC по E_lib
            if components_idxs:
                # Pre-compute mean E per CC для быстрого nearest assignment
                cc_mean_E = [
                    sum(channels[i] for i in cc) / len(cc) for cc in components_idxs
                ]
                for p in phantom_idxs:
                    if visited[p]:
                        continue
                    # nearest CC by abs(channel difference)
                    best = min(
                        range(len(components_idxs)),
                        key=lambda k: abs(channels[p] - cc_mean_E[k]),
                    )
                    components_idxs[best].append(p)
                    visited[p] = True
            else:
                # All-phantom cluster — собираем в один dummy CC и
                # помечаем; downstream apply_multiplet_deconvolution
                # отбросит его (нет active components).
                cc_all_phantom: list = list(phantom_idxs)
                for p in phantom_idxs:
                    visited[p] = True
                if cc_all_phantom:
                    components_idxs.append(cc_all_phantom)

            # Шаг 3 + 4: для каждого CC собрать sub-cluster, применить
            # top-K cap для CC размера > max_components_per_cluster.
            # F-391 — top-K cap считаем ТОЛЬКО по active компонент:
            # phantoms (S/N<gate или F-381 anchors) уже отсеяны и
            # не должны учавствовать в ranking-е. Они присоединены к
            # ближайшему CC для evidence, но всегда остаются phantom.
            #
            # BUG-3 Fix #1 (2026-06-02) — library_anchor_strong BYPASS top-K
            # cap. Strong-anchors добавлены явно потому что физика требует
            # их присутствия для разрешения cluster'а: их demote обратно
            # в phantom (от cap) аннулирует весь смысл Fix #1. Top-K cap
            # применяется только к OBSERVED active'ам (cowell / lsrm /
            # deconvolved); strong-anchors учитываются ПОВЕРХ cap без
            # ranking'а. Если общее число AСTIVE > max после слияния —
            # это физически оправдано (n_active = n_observed + n_strong).
            for cc_idxs in components_idxs:
                sub_cluster = [cl[i] for i in cc_idxs]
                sub_cluster.sort(key=lambda m: m.peak_channel)
                actives = [m for m in sub_cluster if not _is_phantom(m)]
                phantoms_pre = [m for m in sub_cluster if _is_phantom(m)]
                # BUG-3 Fix #1: separate observed actives from strong anchors.
                strong_actives = [
                    m for m in actives
                    if str(getattr(m, "peak_area_source", "") or "")
                    == "library_anchor_strong"
                ]
                observed_actives = [
                    m for m in actives
                    if str(getattr(m, "peak_area_source", "") or "")
                    != "library_anchor_strong"
                ]
                if (len(observed_actives) > max_components_per_cluster
                        and max_components_per_cluster > 0):
                    # Sort observed actives by library_I_pct desc;
                    # top-K активны, остальные становятся phantom anchors.
                    # NB: strong_actives ALWAYS keep — bypass cap.
                    by_intensity = sorted(
                        observed_actives,
                        key=lambda m: float(
                            getattr(m, "library_I_pct", 0.0) or 0.0
                        ),
                        reverse=True,
                    )
                    top_K = by_intensity[:max_components_per_cluster]
                    top_K_ids = {id(m) for m in top_K}
                    strong_ids = {id(m) for m in strong_actives}
                    rebuilt: list = []
                    for m in sub_cluster:
                        if id(m) in top_K_ids or id(m) in strong_ids:
                            rebuilt.append(m)
                            continue
                        if _is_phantom(m):
                            # already phantom (F-391 or F-381)
                            rebuilt.append(m)
                            continue
                        # Demote excess observed active → phantom
                        try:
                            phantom = _dc_replace_lm(
                                m,
                                peak_area=None,
                                peak_area_uncertainty=None,
                                peak_area_source="library_anchor_phantom",
                            )
                        except Exception:
                            try:
                                phantom = m
                                m.peak_area = None
                                m.peak_area_uncertainty = None
                                m.peak_area_source = "library_anchor_phantom"
                            except Exception:
                                phantom = m
                        rebuilt.append(phantom)
                    rebuilt.sort(key=lambda m: m.peak_channel)
                    new_clusters.append(rebuilt)
                    _ = phantoms_pre  # suppress F-841
                else:
                    new_clusters.append(sub_cluster)

        new_clusters.sort(key=lambda g: g[0].peak_channel if g else 0.0)
        clusters = new_clusters

    # Step 3 (LSRM Lzmax): split any cluster whose one-polynomial zone
    # exceeds max_zone_length_fwhm·FWHM at the min-counts valley
    # (LSRM §4, pdf.md:482-494). No-op when disabled (0.0) or spec None.
    if clusters and max_zone_length_fwhm > 0.0 and spec is not None:
        clusters = _split_zones_lzmax(
            clusters,
            spec,
            fwhm_at_channel,
            max_zone_length_fwhm=max_zone_length_fwhm,
            roi_window_factor=lzmax_roi_window_factor,
        )

    return clusters


def _f298_inject_bg_anchors(
    components: list,
    cluster_E_min_keV: float,
    cluster_E_max_keV: float,
    spec,
    fwhm_at_channel: Callable[[float], float],
    *,
    min_intensity_pct: float = 1.0,
) -> list:
    """F-314 / v1.18.13 — inject F-96 канонические фоновые линии как
    additional MultipletComponent в окно энергии cluster.

    При активации (opt-in) ловит ситуации когда реальный fit-окно содержит
    bg-линии (K-40 1461, Tl-208 583/2614, Bi-214 609/1764, ...), которые
    не идентифицированы pipeline-ом (например, на чистом Cs-137 source
    с loose bg residual). Без bg-anchors deconvolve fit относит их counts
    к близким library-linям → artifact bias.

    Returns
    -------
    Расширенный список MultipletComponent (original + bg anchors).
    """
    try:
        from gamma.activity.bg_lines_builder import filter_bg_lines_in_window
    except ImportError:
        return list(components)
    # F-298 API: filter_bg_lines_in_window(E_min_keV, E_max_keV, min_intensity);
    # intensity_decimal (not pct). Convert local min_intensity_pct → decimal.
    bg_lines = filter_bg_lines_in_window(
        E_min_keV=cluster_E_min_keV,
        E_max_keV=cluster_E_max_keV,
        min_intensity=float(min_intensity_pct) / 100.0,
    )
    if not bg_lines:
        return list(components)
    # Dedupe tolerance: ±2 keV (типичное расстояние между library E_keV
    # и точной канонической энергией линии). Без этого зазора K-40 в
    # библиотеке (1461.0 округлённое) и в F-96 (1460.83 точное) считаются
    # разными → дубликат.
    existing_E = [c.line_E_keV for c in components]
    augmented = list(components)
    for bgl in bg_lines:
        E = float(bgl.E_keV)
        is_dup = any(abs(E - e_existing) < 2.0 for e_existing in existing_E)
        if is_dup:
            continue   # уже идентифицирован в кластере
        try:
            ch = float(spec.energy_to_channel(E))
            fwhm = float(fwhm_at_channel(ch))
        except (ValueError, AttributeError):
            continue
        augmented.append(MultipletComponent(
            nuclide=f"bg:{bgl.nuclide}",   # маркер "bg:" чтобы downstream
                                          # отличал anchor от real component
            line_E_keV=E,
            library_I_pct=float(bgl.intensity_decimal) * 100.0,
            center_channel=ch,
            fwhm_channels=fwhm,
        ))
    return augmented


def deconvolve_identified_multiplets(
    identification_result,
    spec,
    fwhm_at_channel: Callable[[float], float],
    *,
    overlap_threshold_fwhm: float = 1.0,
    continuum: str = "step_linear",
    enable_f96_bg_anchors: bool = False,
    bg_anchor_min_intensity_pct: float = 1.0,
    expand_to_display_window: bool = True,
    # F-378 / v1.18.25.1 — intensity-coupled fit вместо legacy
    # independent-Gaussians NNLS. По умолчанию ON: компоненты одного
    # нуклида связываются через ОДНУ свободную амплитуду A_nuc с
    # библиотечными интенсивностями I_k → устраняет singular matrix на
    # close-Е компонентах одного nuclide (M3 V2: χ²/ν 499 → ~5; M3
    # production: χ²/ν 24 → ~6). Legacy path сохранён под флагом
    # use_intensity_coupled=False (debug / regression baseline).
    use_intensity_coupled: bool = True,
    # F-387 / v1.18.26 → F-387.1 / v1.18.26.1 — Rayleigh-CC unresolved
    # criterion. **Semantic changed**: factor·FWHM_avg(per pair), не
    # factor·FWHM_min (global per cluster). Default 1.0 = Rayleigh.
    # См. `find_multiplet_regions` для полного обоснования.
    unresolved_separation_fwhm_factor: float = 1.1,
    # F-387.1 — top-K cap по library_I_pct: sub-cluster крупнее этого
    # → top-K активны, остальные → phantom anchors.
    max_components_per_cluster: int = 3,
    # F-391 / v1.18.27 — S/N significance gate. См. `find_multiplet_regions`.
    min_significance_snr: float = 3.0,
    min_significance_snr_singleton: float = 5.0,
    # F-440 / v1.30.0 — Phase 1 weak-line gate (passthrough to
    # find_multiplet_regions). См. там же для полного описания.
    # F-440 / v1.30.0 — Phase 1 weak-line gate (S/N + I_pct). Default
    # OFF (0.0/0.0) для сохранения baseline behaviour. Включается
    # передачей 5.0/3.0 (бриф F-440 spec). Phase 2 weak-line completion
    # выполняется независимо (всегда), как информационный JSON-блок.
    min_grouping_snr: float = 0.0,
    min_grouping_intensity_pct: float = 0.0,
    # F-447 passthrough — guest-only ROI-owner aware фильтр.
    f440_guest_only_filter: bool = False,
) -> list:
    """
    Convenience: detect every multiplet cluster in an identification
    result and run intensity-coupled deconvolution on each.

    Library line energies and intensities are taken from the
    `LineMatch` objects themselves; the component channel is the
    library energy mapped through `spec.channel_to_energy`-inverse via
    `spec.energy_to_channel`, and the component FWHM is read from
    `fwhm_at_channel(center_channel)`.

    F-378: with `use_intensity_coupled=True` (default), the fit routes
    through `coupled_intensity_fit` with `group=nuclide` — this is the
    same path that FORCED clusters (TH232/RA226) already use. Areas
    within a nuclide group are constrained as a_k = A_nuc · I_k / 100,
    which is the only physically-correct treatment when several lines of
    a single nuclide fall within a few FWHM (the typical case for
    Ac-228 around 503-583 keV or 321-341 keV).

    Returns:
        List of `DeconvolutionResult` in the same order as the
        multiplet clusters returned by `find_multiplet_regions`.
    """
    # F-381 / v1.18.25.2 — подтянуть library из nuclide_library для
    # enrichment кластера всеми library линиями identified нуклидов
    _f381_lib = None
    try:
        from gamma.data.nuclide_library import _load as _f381_load_lib
        _f381_lib = _f381_load_lib()
    except Exception:
        _f381_lib = None

    clusters = find_multiplet_regions(
        identification_result, fwhm_at_channel,
        overlap_threshold_fwhm=overlap_threshold_fwhm,
        expand_to_display_window=expand_to_display_window,
        spec=spec,
        nuclide_library=_f381_lib,
        unresolved_separation_fwhm_factor=unresolved_separation_fwhm_factor,
        max_components_per_cluster=max_components_per_cluster,
        min_significance_snr=min_significance_snr,
        # F-440 passthrough.
        min_grouping_snr=min_grouping_snr,
        min_grouping_intensity_pct=min_grouping_intensity_pct,
        # F-447 passthrough.
        f440_guest_only_filter=f440_guest_only_filter,
    )

    # F-387.2 / v1.18.27.1 — cluster acceptance: 0 actives → drop;
    # 1 active (singleton) → drop из multiplet array (routing к
    # primary_fep через LineMatch в identification_result, см.
    # apply_multiplet_deconvolution для полного обоснования);
    # ≥2 actives → keep.
    def _is_phantom_lm(m) -> bool:
        return str(getattr(m, "peak_area_source", "") or "") in (
            "library_anchor", "library_anchor_phantom",
        )
    filtered_clusters: list = []
    for cl in clusters:
        actives = [m for m in cl if not _is_phantom_lm(m)]
        if len(actives) < 2:
            # F-387.2: 0 active → phantom-only, нечего фитить.
            # 1 active → singleton, маршрут в primary_feps (не multiplet).
            _ = min_significance_snr_singleton  # explicit unused
            continue
        filtered_clusters.append(cl)
    clusters = filtered_clusters

    results = []

    # F-378 — local fwhm_at_energy helper для coupled path:
    # FWHM(E) [keV] = FWHM(ch) [channels] · bin_width_keV(ch).
    # При spec=None (unit-тесты) coupled path отключается → fallback legacy.
    if use_intensity_coupled and spec is not None \
            and getattr(spec, "counts", None) is not None:
        counts_arr = np.asarray(spec.counts, dtype=np.float64)
        n_ch_total = len(counts_arr)
    else:
        counts_arr = None
        n_ch_total = 0

    def _fwhm_at_energy_local(E_keV: float) -> float:
        try:
            ch = float(spec.energy_to_channel(float(E_keV)))
            fw_ch = float(fwhm_at_channel(ch))
            if fw_ch <= 0:
                return 1.0
            # bin width in keV at this channel
            try:
                e1 = float(spec.channel_to_energy(ch + 0.5))
                e0 = float(spec.channel_to_energy(ch - 0.5))
                bin_w = abs(e1 - e0)
            except Exception:
                bin_w = float(
                    getattr(spec, 'bin_width_keV', 1.0) or 1.0
                )
            if bin_w <= 0:
                bin_w = 1.0
            return max(0.1, fw_ch * bin_w)
        except Exception:
            return 1.0

    for cluster in clusters:
        components = []
        # F-387.1 — phantom anchors из top-K cap собираются в отдельный
        # список и прикрепляются к result для downstream
        # reporting/evidence (см. apply_multiplet_deconvolution).
        phantom_components_list: list = []
        for m in cluster:
            ch = float(spec.energy_to_channel(m.library_E_keV))
            fwhm = float(fwhm_at_channel(ch))
            mc = MultipletComponent(
                nuclide=m.nuclide,
                line_E_keV=m.library_E_keV,
                library_I_pct=m.library_I_pct,
                center_channel=ch,
                fwhm_channels=fwhm,
            )
            if (getattr(m, "peak_area_source", "")
                    == "library_anchor_phantom"):
                phantom_components_list.append(mc)
                continue
            components.append(mc)
        # F-314 / v1.18.13 — opt-in inject F-96 bg-anchors как добавочные
        # MultipletComponent чтобы fit не приписывал bg counts близким
        # real-linям.
        if enable_f96_bg_anchors and components:
            E_min = min(c.line_E_keV for c in components)
            E_max = max(c.line_E_keV for c in components)
            # Расширяем окно на ±3 FWHM (типичный fit window)
            try:
                avg_fwhm = sum(c.fwhm_channels for c in components) / len(components)
                fwhm_keV = avg_fwhm * float(getattr(spec, 'bin_width_keV', 1.0) or 1.0)
                E_min_window = E_min - 3.0 * fwhm_keV
                E_max_window = E_max + 3.0 * fwhm_keV
            except Exception:
                E_min_window = E_min - 50.0
                E_max_window = E_max + 50.0
            components = _f298_inject_bg_anchors(
                components, E_min_window, E_max_window,
                spec, fwhm_at_channel,
                min_intensity_pct=bg_anchor_min_intensity_pct,
            )

        # F-378 — intensity-coupled fit ТОЛЬКО при реальной degeneracy.
        # См. apply_multiplet_deconvolution для обоснования.
        # BUG-3 (2026-06-02): расширение trigger'а — coupled также при
        # ≥2 components одного nuclide независимо от separation
        # (strong-anchor enrichment). См. apply_multiplet_deconvolution
        # для полного обоснования.
        _is_degenerate = False
        if components and len(components) >= 2:
            from collections import Counter as _CounterBug3_LEG
            _nuc_counts_leg = _CounterBug3_LEG(
                c.nuclide for c in components if c.nuclide
            )
            if any(n >= 2 for n in _nuc_counts_leg.values()):
                _is_degenerate = True
            if not _is_degenerate:
                for _ic in range(len(components)):
                    for _jc in range(_ic + 1, len(components)):
                        _ci = components[_ic]
                        _cj = components[_jc]
                        if not (_ci.nuclide and _cj.nuclide
                                and _ci.nuclide == _cj.nuclide):
                            continue
                        _sep = abs(_ci.center_channel - _cj.center_channel)
                        _fw_min = min(_ci.fwhm_channels, _cj.fwhm_channels)
                        if _fw_min > 0 and _sep < 0.4 * _fw_min:
                            _is_degenerate = True
                            break
                    if _is_degenerate:
                        break

        # F-378 — preferred path: intensity-coupled fit (only when degenerate)
        if (use_intensity_coupled and _is_degenerate
                and components and counts_arr is not None):
            try:
                from gamma.peaks.coupled_multiplet import (
                    coupled_intensity_fit, ComponentSpec,
                )
                # ROI bounds в каналах
                ch_lo = int(math.floor(
                    min(c.center_channel for c in components)
                    - 2.5 * max(c.fwhm_channels for c in components)
                ))
                ch_hi = int(math.ceil(
                    max(c.center_channel for c in components)
                    + 2.5 * max(c.fwhm_channels for c in components)
                )) + 1
                ch_lo = max(0, ch_lo)
                ch_hi = min(n_ch_total, ch_hi)
                if ch_hi - ch_lo < len(components) + 3:
                    raise ValueError("ROI too narrow for coupled fit")
                roi_E = np.array(
                    [spec.channel_to_energy(c)
                     for c in range(ch_lo, ch_hi)],
                    dtype=np.float64,
                )
                roi_counts = counts_arr[ch_lo:ch_hi]
                # Группировка по nuclide: компоненты одного nuclide с
                # ненулевой I → group=nuclide; компоненты без nuclide
                # или с I≤0 → независимые (group="").
                comp_specs = []
                for c in components:
                    has_group = (
                        bool(c.nuclide)
                        and float(c.library_I_pct or 0.0) > 0.0
                    )
                    comp_specs.append(ComponentSpec(
                        nuclide=str(c.nuclide or ""),
                        E_keV=float(c.line_E_keV),
                        I_gamma_pct=float(c.library_I_pct or 0.0),
                        group=str(c.nuclide) if has_group else "",
                    ))
                coupled = coupled_intensity_fit(
                    energy_keV=roi_E,
                    counts=roi_counts,
                    components=comp_specs,
                    fwhm_at=_fwhm_at_energy_local,
                    continuum=continuum,
                    roi_low_ch=ch_lo,
                    cluster_id="",
                    title="",
                    use_peak_image=False,
                )
                _r = _coupled_to_deconv_result(
                    coupled, fwhm_at_channel,
                )
                if phantom_components_list:
                    from dataclasses import replace as _dc_replace_dr2
                    _r = _dc_replace_dr2(
                        _r,
                        phantom_components=tuple(phantom_components_list),
                    )
                results.append(_r)
                continue
            except Exception:
                # graceful fallback к legacy path при любой ошибке
                pass

        _r = deconvolve_multiplet(
            spec.counts,
            components=components,
            continuum=continuum,
        )
        if phantom_components_list:
            from dataclasses import replace as _dc_replace_dr3
            _r = _dc_replace_dr3(
                _r,
                phantom_components=tuple(phantom_components_list),
            )
        results.append(_r)
    return results


# F-272 (v1.17.11, T-058) — detector-class рекомендованные пороги χ²/ν
# для acceptance мультиплетных fit'ов. Различаются драматически между
# HPGe (~1) и сцинтилляторами (NaI 63x63 имеет естественную статистику
# Poisson с small-counts регионами, где χ²/ν ≤ 6 — типичное значение).
# Использовать жёсткий порог ~1 для NaI ведёт к false-rejection валидных
# fit'ов. См. ЛСРМ Алгоритмические основы 2022 §8 — рекомендуется
# χ²/ν ≤ 6 для NaI при INL <1%.
RECOMMENDED_CHI2_THRESHOLD = {
    "HPGe":   2.0,
    "CdZnTe": 3.0,
    "LaBr3":  4.0,
    "CeBr3":  4.0,
    "NaI":    6.0,
    "CsI":    6.0,
}


def recommended_chi2_threshold(detector_class: str) -> float:
    """F-272 — вернуть рекомендованный порог χ²/ν для acceptance fit'a
    мультиплета для данного класса детектора.

    Используется как default при вызове ``apply_multiplet_deconvolution``
    или ``staged_pipeline``. Возвращает 6.0 если класс неизвестен
    (NaI-conservative — не reject'ить лишнее).
    """
    return float(RECOMMENDED_CHI2_THRESHOLD.get(str(detector_class), 6.0))


def apply_multiplet_deconvolution(
    identification_result,
    spec,
    fwhm_at_channel: Callable[[float], float],
    *,
    overlap_threshold_fwhm: float = 1.0,
    continuum: str = "step_linear",
    max_chi2_per_dof: float = float("inf"),
    enable_f96_bg_anchors: bool = False,
    bg_anchor_min_intensity_pct: float = 1.0,
    expand_to_display_window: bool = True,
    # F-378 / v1.18.25.1 — intensity-coupled fit для auto-detected
    # кластеров. См. документацию в `deconvolve_identified_multiplets`.
    use_intensity_coupled: bool = True,
    # F-387 / v1.18.26 → F-387.1 / v1.18.26.1 — Rayleigh-CC criterion.
    # Default 1.0 = Rayleigh (factor·FWHM_avg per pair). Семантика
    # **изменена** в v1.18.26.1 — было factor·FWHM_min с keep-monolith,
    # стало factor·FWHM_avg с CC-split + top-K cap. См.
    # `find_multiplet_regions` для полного обоснования.
    unresolved_separation_fwhm_factor: float = 1.1,
    # F-387.1 — top-K cap по library_I_pct в sub-cluster крупнее этого
    # значения. Остальные → phantom anchors (не fit'ятся).
    max_components_per_cluster: int = 3,
    # F-391 / v1.18.27 — S/N significance gate. Multiplet member threshold
    # (Gilmore §5.5: пик «реально измерен» при S/N ≥ 3, Currie L_C k=3).
    # См. `find_multiplet_regions` для полного обоснования.
    min_significance_snr: float = 3.0,
    # F-391 — singleton acceptance threshold: cluster после CC split,
    # содержащий ТОЛЬКО 1 active компоненту (CC=1), принимается обратно
    # в primary_fep list если S/N этой компоненты ≥ singleton_threshold.
    # Иначе cluster дропается полностью. Default 5.0 — более строго чем
    # multiplet member (3.0), потому что singleton-фит против continuum'а
    # без anchor'ов от соседей легко даёт false positives на шуме.
    min_significance_snr_singleton: float = 5.0,
    # F-440 / v1.30.0 — Phase 1 weak-line gate passthrough.
    min_grouping_snr: float = 0.0,
    min_grouping_intensity_pct: float = 0.0,
    # F-447 passthrough — guest-only ROI-owner aware фильтр.
    f440_guest_only_filter: bool = False,
) -> tuple:
    """
    Post-pass: replace `LineMatch.peak_area` with deconvolved areas for
    every component that lies in a multiplet cluster, leaving isolated
    lines untouched.

    This is the pipeline-integration entry point for F-33: identification
    runs first and yields per-`LineMatch` peak areas integrated by Cowell
    or the Lsrm peak table; this function then revisits any `LineMatch`
    whose peak channel collides with another `LineMatch` (within
    `overlap_threshold_fwhm × FWHM`) and substitutes the deconvolved
    area.

    Behaviour:
      - Lines NOT in any cluster: returned `LineMatch` is unchanged.
      - Lines IN a cluster, fit converged and `chi²/ν ≤ max_chi2_per_dof`:
        `peak_area` and `peak_area_uncertainty` replaced with the
        deconvolved values, `peak_area_source = "deconvolved"`.
      - Lines IN a cluster, fit failed OR χ²/ν above the threshold: the
        `LineMatch` is left as-is and the cluster's
        `DeconvolutionResult` carries an informative `notes` string.

    Args:
        identification_result: result from `identify_nuclides` (already
            disambiguated if desired).
        spec: the `Spectrum` used to produce that identification — must
            have `.counts` and `.energy_to_channel`.
        fwhm_at_channel: callable mapping channel → FWHM (in channels).
        overlap_threshold_fwhm: same semantic as in
            `find_multiplet_regions`. Default 1.0 (FWHM-touching pairs
            count as a multiplet).
        continuum: "linear" or "step_linear" (default).
        max_chi2_per_dof: skip area replacement for clusters whose fit
            quality is worse than this. Default infinity (always
            replace when converged).

    Returns:
        (new_identification_result, deconvolutions) where
        `new_identification_result` is an `IdentificationResult` with
        updated `LineMatch` entries and `deconvolutions` is the list of
        `DeconvolutionResult` (one per cluster, in the same order as
        `find_multiplet_regions` would return them).
    """
    from dataclasses import replace as _dc_replace
    from gamma.identification.identify import (
        IdentificationResult, NuclideIdentification, LineMatch,
    )

    # F-381 / v1.18.25.2 — enrich clusters всеми library lines
    # identified нуклидов (см. find_multiplet_regions).
    _f381_lib = None
    try:
        from gamma.data.nuclide_library import _load as _f381_load_lib
        _f381_lib = _f381_load_lib()
    except Exception:
        _f381_lib = None

    clusters = find_multiplet_regions(
        identification_result, fwhm_at_channel,
        overlap_threshold_fwhm=overlap_threshold_fwhm,
        expand_to_display_window=expand_to_display_window,
        spec=spec,
        nuclide_library=_f381_lib,
        unresolved_separation_fwhm_factor=unresolved_separation_fwhm_factor,
        max_components_per_cluster=max_components_per_cluster,
        min_significance_snr=min_significance_snr,
        # F-440 passthrough.
        min_grouping_snr=min_grouping_snr,
        min_grouping_intensity_pct=min_grouping_intensity_pct,
        # F-447 passthrough.
        f440_guest_only_filter=f440_guest_only_filter,
    )

    # F-387.2 / v1.18.27.1 — singleton CC=1 routing к primary_fep.
    # После F-387.1 Rayleigh-CC split возможны sub-cluster'ы с
    # `n_active == 1` (один реальный component + произвольное число
    # phantom anchors из F-381 enrichment / F-387.1 top-K cap demote).
    # F-387.2 contract: «1-component multiplet» физически не имеет
    # смысла — это просто peak, который должен жить в primary_feps,
    # а не в multiplet_deconvolutions. Маршрутизация:
    #
    #   - 0 actives (phantom-only) → drop полностью (нечего фитить).
    #
    #   - 1 active с S/N < singleton_threshold → drop полностью.
    #     LineMatch уже помечен phantom через F-391 gate ИЛИ останется
    #     с peak_area_source="cowell" но с низким S/N в primary_feps.
    #     Multiplet entry не создаётся (нет multiplet evidence от
    #     соседей, single-component fit против continuum'а на шуме
    #     даёт false positives).
    #
    #   - 1 active с S/N ≥ singleton_threshold → **drop из multiplet array**
    #     (F-387.2 / v1.18.27.1, было «keep как 1-active multiplet»
    #     в v1.18.27 F-391 partial implementation). LineMatch для этого
    #     active component уже присутствует в
    #     identification_result.matched_lines с оригинальным
    #     peak_area_source (e.g. "cowell" / "lsrm_peaks_table") и
    #     measured peak_area. Через `_build_primary_feps` он попадает
    #     в JSON primary_feps без дополнительной обработки. Phantom
    #     anchors, attached к этому CC (F-381 library enrichment или
    #     F-387.1 demote), дропаются вместе с singleton — они имели
    #     смысл только как evidence в multiplet context (≥2 active
    #     compongnents), а singleton — это уже не multiplet.
    #
    #   - ≥2 actives → multiplet fit как раньше (F-381 evidence,
    #     F-387.1 split, F-391 phantom attach все сохранены).
    #
    # Обоснование «singleton → primary_feps»: исходно (v1.18.27 F-391)
    # high-S/N singleton сохранялся в multiplet_deconvolutions как
    # 1-component entry — это давало refined area через single-component
    # fit против continuum'а. На реальной Th-232 demo это приводило к
    # «singleton multiplet» entries (M3 V2: 6c но 1 active, M7 V2: 2c
    # но 1 active) — UX-нежелательные: пользователь видит «multiplet»
    # с одной реальной линией. F-387.2 убирает эту аномалию: либо
    # multiplet (≥2 active overlapping lines), либо primary FEP — нет
    # промежуточного «singleton multiplet» состояния.
    def _is_phantom_lm(m) -> bool:
        return str(getattr(m, "peak_area_source", "") or "") in (
            "library_anchor", "library_anchor_phantom",
        )

    filtered_clusters: list = []
    for cl in clusters:
        actives = [m for m in cl if not _is_phantom_lm(m)]
        n_active = len(actives)
        if n_active == 0:
            # Phantom-only cluster — drop полностью (нечего фитить).
            continue
        if n_active == 1:
            # F-387.2: singleton → drop из multiplet array независимо
            # от S/N. High-S/N singleton маршрутизируется к primary_feps
            # через LineMatch в identification_result. Low-S/N — также
            # дропается (нет multiplet evidence для refinement'а).
            # Параметр min_significance_snr_singleton сохранён в API
            # для back-compat / диагностических вызовов, но больше не
            # влияет на routing (всегда drop).
            _ = min_significance_snr_singleton  # F-387.2: explicit unused
            continue
        filtered_clusters.append(cl)
    clusters = filtered_clusters

    # (nuclide, library_E_keV_rounded) → (area, area_uncertainty)
    replacements: dict = {}
    deconvolutions = []

    # F-378 — fwhm_at_energy adapter для coupled path:
    # FWHM(E) [keV] = FWHM(ch) [channels] · bin_width_keV(ch)
    def _fwhm_at_energy_local(E_keV: float) -> float:
        try:
            ch = float(spec.energy_to_channel(float(E_keV)))
            fw_ch = float(fwhm_at_channel(ch))
            if fw_ch <= 0:
                return 1.0
            try:
                e1 = float(spec.channel_to_energy(ch + 0.5))
                e0 = float(spec.channel_to_energy(ch - 0.5))
                bin_w = abs(e1 - e0)
            except Exception:
                bin_w = float(
                    getattr(spec, 'bin_width_keV', 1.0) or 1.0
                )
            if bin_w <= 0:
                bin_w = 1.0
            return max(0.1, fw_ch * bin_w)
        except Exception:
            return 1.0

    # F-378 — counts_arr / n_ch_total нужны для coupled-path ROI вычислений.
    # При spec=None (некоторые unit-тесты) фоллбэк на legacy path.
    if spec is not None and getattr(spec, "counts", None) is not None:
        counts_arr = np.asarray(spec.counts, dtype=np.float64)
        n_ch_total = len(counts_arr)
    else:
        counts_arr = None
        n_ch_total = 0

    # ─────────────────────────────────────────────────────────────────
    # BUG-3 Fix #2 (2026-06-02) — pre-filter weak active candidates
    # ─────────────────────────────────────────────────────────────────
    # Корень: peak-detector иногда генерирует «шумовые» candidates на
    # слабые библиотечные линии (Tl-208 233 keV, I=0.11%, в Th-232 M3) —
    # они проходят S/N gate потому что peak-search нашёл noisy bump на
    # склоне сильной соседней линии. Fitter затем тратит DoF на эту
    # линию, амплитуду «занижает» (или искажает соседей).
    #
    # Fix: per-cluster relative intensity threshold. Среди active
    # candidates (НЕ phantoms) вычисляем I_max; всё что < 5%·I_max
    # демотим в library_anchor_phantom (так оно остаётся в evidence
    # но не участвует в fit). Strong-anchor (BUG-3 Fix #1) защищены —
    # они в active по построению (I_pct ≥ 5%·I_max_in_window).
    #
    # Контракт: filter применяется ПОСЛЕ BUG-3 Fix #1 strong-anchor
    # enrichment и ПОСЛЕ F-387.2 ≥2-actives validation, но ДО fit-loop.
    BUG3_WEAK_REL_THRESHOLD = 0.05  # 5% of I_max per cluster
    from dataclasses import replace as _dc_replace_bug3_f2
    for cl_idx, cl in enumerate(clusters):
        active_idxs_lm = [
            i for i, m in enumerate(cl)
            if str(getattr(m, "peak_area_source", "") or "") not in (
                "library_anchor", "library_anchor_phantom",
            )
        ]
        if len(active_idxs_lm) < 2:
            continue
        I_pcts = [
            float(getattr(cl[i], "library_I_pct", 0.0) or 0.0)
            for i in active_idxs_lm
        ]
        I_max_active = max(I_pcts) if I_pcts else 0.0
        if I_max_active <= 0:
            continue
        weak_threshold = BUG3_WEAK_REL_THRESHOLD * I_max_active
        for i in active_idxs_lm:
            m = cl[i]
            I_pct = float(getattr(m, "library_I_pct", 0.0) or 0.0)
            if I_pct >= weak_threshold:
                continue
            # Демотим weak в phantom — сохраняет evidence, убирает
            # из fit-list. NB: strong (library_anchor_strong) сюда
            # не попадает — он по построению ≥ 5%·I_max_in_window
            # и таким образом не «weak» в smaller-or-equal smysle.
            # Но на всякий случай проверим source — strong не демотим.
            src = str(getattr(m, "peak_area_source", "") or "")
            if src == "library_anchor_strong":
                continue
            try:
                cl[i] = _dc_replace_bug3_f2(
                    m, peak_area_source="library_anchor_phantom",
                )
            except Exception:
                try:
                    m.peak_area_source = "library_anchor_phantom"
                except Exception:
                    pass

    # F-387.2 re-validation: after Fix #2 demote, кластер может скатиться
    # ниже 2-actives → дропаем.
    filtered_after_bug3_f2: list = []
    for cl in clusters:
        actives_chk = [
            m for m in cl
            if str(getattr(m, "peak_area_source", "") or "") not in (
                "library_anchor", "library_anchor_phantom",
            )
        ]
        if len(actives_chk) >= 2:
            filtered_after_bug3_f2.append(cl)
    clusters = filtered_after_bug3_f2

    for cluster in clusters:
        components = []
        # F-387.1 — phantom anchors из top-K cap собираем в отдельный
        # список для прикрепления к DeconvolutionResult (downstream
        # reporting/evidence). Они НЕ становятся отдельными
        # MultipletComponent в fit-list.
        phantom_components_list: list = []
        for m in cluster:
            ch = float(spec.energy_to_channel(m.library_E_keV))
            fwhm = float(fwhm_at_channel(ch))
            mc = MultipletComponent(
                nuclide=m.nuclide,
                line_E_keV=m.library_E_keV,
                library_I_pct=m.library_I_pct,
                center_channel=ch,
                fwhm_channels=fwhm,
            )
            if (getattr(m, "peak_area_source", "")
                    == "library_anchor_phantom"):
                phantom_components_list.append(mc)
                continue
            components.append(mc)
        # F-322 / v1.18.16 — opt-in F-96 bg-anchors injection в cluster
        # для покрытия неучтённых фоновых линий (например, 511 keV
        # annihilation в M3 кластере Tl-208 510.77 + 583.19). Без bg-anchors
        # fit относит counts annihilation к ближайшей library line → χ²↑.
        if enable_f96_bg_anchors and components:
            E_min = min(c.line_E_keV for c in components)
            E_max = max(c.line_E_keV for c in components)
            try:
                avg_fwhm = sum(c.fwhm_channels for c in components) / len(components)
                fwhm_keV = avg_fwhm * float(getattr(spec, 'bin_width_keV', 1.0) or 1.0)
                E_min_w = E_min - 3.0 * fwhm_keV
                E_max_w = E_max + 3.0 * fwhm_keV
            except Exception:
                E_min_w = E_min - 50.0
                E_max_w = E_max + 50.0
            components = _f298_inject_bg_anchors(
                components, E_min_w, E_max_w,
                spec, fwhm_at_channel,
                min_intensity_pct=bg_anchor_min_intensity_pct,
            )

        # F-378 — intensity-coupled fit ТОЛЬКО когда кластер реально
        # degenerate: ≥1 пара компонент одного nuclide ближе чем
        # 0.4·FWHM_min. Иначе → legacy NNLS (back-compat для chain-ratio
        # тестов и production-pipeline activity вычислений).
        # Корень: на близко-E линиях одного nuclide design matrix
        # становится сингулярной → площади «слипаются» в одной
        # компоненте, остальные → ≈0. Intensity-coupling решает это,
        # связывая их через A_nuc·I_k/100.
        #
        # BUG-3 (2026-06-02) — расширение trigger'а: coupled-path также
        # запускается, если кластер содержит ≥2 компоненты одного
        # nuclide (любая separation, не только < 0.4·FWHM_min). Корень
        # на Th-232 M3: strong-anchor enrichment добавляет Ac-228 209/270
        # и Pb-212 238/300 в один cluster (Δ=61 keV для Ac-228; 62 keV
        # для Pb-212), при FWHM ≈ 25 keV в области 250 keV это
        # 2.5·FWHM — НЕ degenerate по старому критерию, но physically
        # это lines одного nuclide с фиксированным I-ratio (3.89:3.46
        # для Ac-228, 43.6:3.30 для Pb-212). Без coupling fit'er тратит
        # 4 free амплитуды и легко уходит в локальный минимум, где Ac-228
        # 270 «съедает» Pb-212 238 counts. Coupling reduces до 2 free
        # амплитуд (A_Ac228, A_Pb212) и forces физический I-ratio.
        is_degenerate = False
        if components and len(components) >= 2:
            # Per-nuclide multi-component check.
            from collections import Counter as _CounterBug3
            nuc_counts = _CounterBug3(
                c.nuclide for c in components if c.nuclide
            )
            if any(n >= 2 for n in nuc_counts.values()):
                is_degenerate = True
            if not is_degenerate:
                # Legacy proximity check (back-compat).
                for i_c in range(len(components)):
                    for j_c in range(i_c + 1, len(components)):
                        ci = components[i_c]
                        cj = components[j_c]
                        if not (ci.nuclide and cj.nuclide
                                and ci.nuclide == cj.nuclide):
                            continue
                        sep = abs(ci.center_channel - cj.center_channel)
                        fw_min = min(ci.fwhm_channels, cj.fwhm_channels)
                        if fw_min > 0 and sep < 0.4 * fw_min:
                            is_degenerate = True
                            break
                    if is_degenerate:
                        break

        res = None
        if (use_intensity_coupled and is_degenerate
                and components and counts_arr is not None):
            try:
                from gamma.peaks.coupled_multiplet import (
                    coupled_intensity_fit, ComponentSpec,
                )
                ch_lo = int(math.floor(
                    min(c.center_channel for c in components)
                    - 2.5 * max(c.fwhm_channels for c in components)
                ))
                ch_hi = int(math.ceil(
                    max(c.center_channel for c in components)
                    + 2.5 * max(c.fwhm_channels for c in components)
                )) + 1
                ch_lo = max(0, ch_lo)
                ch_hi = min(n_ch_total, ch_hi)
                if ch_hi - ch_lo < len(components) + 3:
                    raise ValueError("ROI too narrow for coupled fit")
                roi_E = np.array(
                    [spec.channel_to_energy(c)
                     for c in range(ch_lo, ch_hi)],
                    dtype=np.float64,
                )
                roi_counts = counts_arr[ch_lo:ch_hi]
                comp_specs = []
                for c in components:
                    has_group = (
                        bool(c.nuclide)
                        and float(c.library_I_pct or 0.0) > 0.0
                    )
                    comp_specs.append(ComponentSpec(
                        nuclide=str(c.nuclide or ""),
                        E_keV=float(c.line_E_keV),
                        I_gamma_pct=float(c.library_I_pct or 0.0),
                        group=str(c.nuclide) if has_group else "",
                    ))
                # F-392 / v1.18.27 — auto-select multi-step continuum
                # для широких multi-anchor coupled-fit ROI. ROI span
                # из roi_E (display-window expanded), а не span компонент.
                # F-392.1 / v1.18.27.1 — phantom anchors (F-387.1 top-K cap
                # demote и F-381 library_anchor) учитываются в anchor-list,
                # если их библиотечная энергия лежит ВНУТРИ fit-ROI. Это
                # физически корректно: continuum step jumps под Compton
                # edges of strong γ-lines существуют в спектре независимо
                # от того, фитится ли амплитуда этой линии. Пример Th-232
                # M3 PROD (ROI ≈ 480-590 кэВ active 503+509+562, phantom
                # Tl-208 510 I=8.1%, Tl-208 583 I=30.5%, Ac-228 463 I=4.4%):
                # без phantom anchors F-392 не активируется (max active
                # I=0.87%); с phantom anchors → 2-3 intense ≥5% in-ROI →
                # step_linear_multi даёт реалистичный continuum.
                roi_span_keV = (
                    float(roi_E[-1] - roi_E[0])
                    if len(roi_E) >= 2 else 0.0
                )
                # F-392.1: собрать anchor pool = active comp_specs +
                # in-ROI phantoms. ROI bounds — fit-ROI [roi_E[0], roi_E[-1]].
                _roi_lo_kev = float(roi_E[0]) if len(roi_E) else 0.0
                _roi_hi_kev = float(roi_E[-1]) if len(roi_E) else 0.0
                _f392_anchor_pool = [
                    (cs.E_keV, cs.I_gamma_pct) for cs in comp_specs
                ]
                for _ph in phantom_components_list:
                    _ph_E = float(_ph.line_E_keV)
                    if _roi_lo_kev <= _ph_E <= _roi_hi_kev:
                        _ph_I = float(_ph.library_I_pct or 0.0)
                        _f392_anchor_pool.append((_ph_E, _ph_I))
                continuum_coupled = _f392_auto_select_continuum(
                    _f392_anchor_pool,
                    continuum,
                    roi_e_span_keV=roi_span_keV,
                )
                coupled = coupled_intensity_fit(
                    energy_keV=roi_E,
                    counts=roi_counts,
                    components=comp_specs,
                    fwhm_at=_fwhm_at_energy_local,
                    continuum=continuum_coupled,
                    roi_low_ch=ch_lo,
                    cluster_id="",
                    title="",
                    use_peak_image=False,
                )
                res = _coupled_to_deconv_result(coupled, fwhm_at_channel)
            except Exception:
                res = None  # graceful fallback below

        if res is None:
            # F-392: legacy channel-space fit поддерживает только
            # "linear" / "step_linear" — multi не имеет здесь смысла
            # (отсутствует energy-aware anchor placement).
            legacy_continuum = (
                "step_linear" if continuum == "step_linear_multi"
                else continuum
            )
            res = deconvolve_multiplet(
                spec.counts,
                components=components,
                continuum=legacy_continuum,
            )
        # F-387.1 — прикрепить phantom_components к result. Frozen
        # dataclass → используем replace.
        if phantom_components_list:
            from dataclasses import replace as _dc_replace_dr
            res = _dc_replace_dr(
                res,
                phantom_components=tuple(phantom_components_list),
            )
        deconvolutions.append(res)
        if not res.converged or res.chi2_per_dof > max_chi2_per_dof:
            continue
        for comp, area, unc in zip(res.components, res.areas,
                                   res.area_uncertainties):
            key = (comp.nuclide, round(comp.line_E_keV, 3))
            replacements[key] = (area, unc)

    if not replacements:
        return identification_result, deconvolutions

    new_detected = []
    for ni in identification_result.detected_nuclides:
        new_matches = []
        any_change = False
        for m in ni.matched_lines:
            key = (m.nuclide, round(m.library_E_keV, 3))
            if key in replacements:
                area, unc = replacements[key]
                new_matches.append(_dc_replace(
                    m,
                    peak_area=float(area),
                    peak_area_uncertainty=float(unc),
                    peak_area_source="deconvolved",
                ))
                any_change = True
            else:
                new_matches.append(m)
        if any_change:
            new_detected.append(_dc_replace(
                ni, matched_lines=tuple(new_matches),
            ))
        else:
            new_detected.append(ni)

    notes = identification_result.notes
    n_replaced = sum(
        1 for ni in new_detected for m in ni.matched_lines
        if m.peak_area_source == "deconvolved"
    )
    note_line = (
        f"Multiplet deconvolution: {len(clusters)} cluster(s), "
        f"{n_replaced} peak area(s) replaced "
        f"(overlap threshold = {overlap_threshold_fwhm:.2f}·FWHM, "
        f"continuum = {continuum})"
    )
    new_notes = (notes + "\n" + note_line).strip() if notes else note_line

    new_result = IdentificationResult(
        detector_type=identification_result.detector_type,
        window=identification_result.window,
        candidates_considered=identification_result.candidates_considered,
        detected_nuclides=tuple(new_detected),
        rejected_nuclides=identification_result.rejected_nuclides,
        unmatched_peaks=identification_result.unmatched_peaks,
        notes=new_notes,
    )
    return new_result, deconvolutions


# ============================================================================
# F-118 / v1.17.5 — Chain-forced intensity-coupled multiplets
# ============================================================================

# Жёстко закреплённые ROI для Th-232 цепочки. Эмиттируются всегда,
# когда chain_dominance.dominant_chain == "Th-232" (методологический
# контракт из references/demo_contract_v1_17_2/multiplet_M{1,2}_coupled.json).
TH232_FORCED_CLUSTERS = (
    {
        "id": "M1",
        "E_lo_keV": 750.0,
        "E_hi_keV": 1115.0,
        "components": (
            ("Ac-228", 911.204, 25.8, "Ac-228"),
            ("Ac-228", 964.77,  4.99, "Ac-228"),
            ("Ac-228", 968.971, 15.8, "Ac-228"),
            ("Tl-208", 860.6,   4.5,  ""),
        ),
        "title": (
            "Мультиплет M1 (связанная подгонка) — Ac-228 911 + 964.8 + "
            "969 + Tl-208 860.6"
        ),
    },
    {
        "id": "M2",
        "E_lo_keV": 1430.0,
        "E_hi_keV": 1790.0,
        "components": (
            ("Ac-228", 1588.2,  3.22, "Ac-228"),
            ("Bi-212", 1620.5,  1.49, ""),
            ("Ac-228", 1630.6,  1.6,  "Ac-228"),
        ),
        "title": (
            "Мультиплет M2 (связанная подгонка) — Ac-228 1588 + 1630 + "
            "Bi-212 1620"
        ),
    },
)


# F-121 / v1.17.6 — зеркальный контракт для цепочки Ra-226 / U-238.
# При chain_dominance.u238 == True всегда эмиттируются три связанных
# подгонки кластеров:
#   U1 (Pb-214 295+352)       — самые яркие линии Pb-214 на NaI
#   U2 (Bi-214 609+665+806)   — кластер вокруг 609 кэВ (характеристический)
#   U3 (Bi-214 1120+1238+1378+1408+1764) — высокоэнергетический кластер;
#       1764 — изолированная линия, но включена в U3 как trump card,
#       поскольку при NaI 63×63 окружена двумя слабыми линиями Bi-214.
# I_γ из библиотечной базы (LNHB / IAEA).
RA226_FORCED_CLUSTERS = (
    {
        "id": "U1",
        "E_lo_keV": 270.0,
        "E_hi_keV": 380.0,
        "components": (
            ("Pb-214",  295.22, 18.42, "Pb-214"),
            ("Pb-214",  351.93, 35.6,  "Pb-214"),
        ),
        "title": (
            "Мультиплет U1 (связанная подгонка) — Pb-214 295 + 352"
        ),
    },
    {
        "id": "U2",
        "E_lo_keV": 580.0,
        "E_hi_keV": 830.0,
        "components": (
            ("Bi-214",  609.31, 45.49, "Bi-214"),
            ("Bi-214",  665.45,  1.531, "Bi-214"),
            ("Bi-214",  768.36,  4.892, "Bi-214"),
            ("Bi-214",  806.17,  1.262, "Bi-214"),
        ),
        "title": (
            "Мультиплет U2 (связанная подгонка) — Bi-214 609 + 665 + "
            "768 + 806"
        ),
    },
    {
        "id": "U3",
        "E_lo_keV": 1080.0,
        "E_hi_keV": 1830.0,
        "components": (
            ("Bi-214", 1120.29, 14.92, "Bi-214"),
            ("Bi-214", 1238.11,  5.834, "Bi-214"),
            ("Bi-214", 1377.67,  3.968, "Bi-214"),
            ("Bi-214", 1408.01,  2.389, "Bi-214"),
            ("Bi-214", 1729.60,  2.878, "Bi-214"),
            ("Bi-214", 1764.49, 15.31,  "Bi-214"),
        ),
        "title": (
            "Мультиплет U3 (связанная подгонка) — Bi-214 1120 + 1238 + "
            "1378 + 1408 + 1730 + 1764"
        ),
    },
)


def _coupled_to_deconv_result(coupled, fwhm_at_channel) -> DeconvolutionResult:
    """Преобразовать CoupledFitResult → DeconvolutionResult так, чтобы
    `apply_multiplet_deconvolution`-совместимая выдача и downstream
    JSON-репортер увидели M1/M2 как обычные deconvolution clusters.
    """
    components_tuple = []
    areas_tuple = []
    sigma_tuple = []
    for cf in coupled.components:
        ch = float(cf.E_keV)  # E in keV; centre_channel here is informative
        fwhm_ch = float(fwhm_at_channel(ch)) if fwhm_at_channel is not None else 1.0
        components_tuple.append(MultipletComponent(
            nuclide=cf.nuclide,
            line_E_keV=float(cf.E_keV),
            library_I_pct=float(cf.I_pct),
            center_channel=ch,
            fwhm_channels=fwhm_ch,
        ))
        areas_tuple.append(float(cf.area))
        sigma_tuple.append(float(cf.sigma_area))
    # F-134 / v1.17.7 — overlay arrays напрямую из CoupledFitResult.
    # Это даёт PNG/HTML рендеру точно те же кривые, что использовались
    # для расчёта χ²/ν и closure, без перевычисления в каналах.
    overlay_E = tuple(float(v) for v in (coupled.E_keV or ()))
    overlay_data = tuple(float(v) for v in (coupled.data or ()))
    overlay_cont = tuple(float(v) for v in (coupled.continuum or ()))
    overlay_total = tuple(float(v) for v in (coupled.total or ()))
    overlay_comps_list = []
    for row in (coupled.component_g_plus_cont or ()):
        overlay_comps_list.append(tuple(float(v) for v in row))
    overlay_comps = tuple(overlay_comps_list) if overlay_comps_list else None

    return DeconvolutionResult(
        components=tuple(components_tuple),
        areas=tuple(areas_tuple),
        area_uncertainties=tuple(sigma_tuple),
        continuum_params=tuple(coupled.continuum_params),
        continuum_model=coupled.continuum_model,
        chi2_per_dof=float(coupled.chi2_per_dof),
        n_dof=int(coupled.n_dof),
        roi_low_ch=int(coupled.roi_low_ch),
        roi_high_ch=int(coupled.roi_high_ch),
        gross_counts=float(sum(coupled.data)),
        converged=bool(coupled.converged),
        method=f"coupled_{coupled.method}",
        degenerate_pairs=(),
        covariance=None,
        notes=(
            f"F-117 связанная подгонка: closure Δ={coupled.closure_pct:.2f}%, "
            f"χ²/ν={coupled.chi2_per_dof:.2f}, метод={coupled.method}, "
            f"кластер {coupled.id or '?'}"
        ),
        # F-134
        overlay_E_keV=overlay_E or None,
        overlay_data=overlay_data or None,
        overlay_continuum=overlay_cont or None,
        overlay_total=overlay_total or None,
        overlay_components=overlay_comps,
        # F-145 / v1.17.8 — Phase A free-centroid side-fit
        centroid_shifts_keV=tuple(
            float(v) for v in (coupled.centroid_shifts_keV or ())
        ),
        phase_A_chi2_per_dof=(
            float(coupled.phase_A_chi2_per_dof)
            if coupled.phase_A_chi2_per_dof is not None else None
        ),
        phase_A_converged=bool(coupled.phase_A_converged),
        cluster_id=str(coupled.id or ""),
        # F-392.1 / v1.18.29 — propagate multi-step anchors / threshold
        # из CoupledFitResult в downstream JSON. getattr защищает legacy
        # CoupledFitResult без полей (старые pickled snapshots / mocks).
        multi_step_anchors=tuple(
            (float(e), float(s))
            for e, s in (getattr(coupled, "multi_step_anchors", ()) or ())
        ),
        multi_step_intensity_threshold_pct=(
            float(coupled.multi_step_intensity_threshold_pct)
            if getattr(coupled, "multi_step_intensity_threshold_pct", None) is not None
            else None
        ),
    )


def run_chain_forced_multiplets(
    spec,
    fwhm_at_channel: Callable[[float], float],
    fwhm_at_energy: Callable[[float], float],
    chain_dominance,
    filename_isotope_hints=None,
    *,
    use_peak_image: bool = False,
    detector_type: str = "NaI",
    use_T_E_model: bool = False,
    nonlinear_refine: bool = False,
    # F-133 / v1.17.7 — per-line ступенька под пиком (ГОСТ форма)
    h_step: Optional[float] = None,
    # F-145 / v1.17.8 — free-centroid Phase A side-fit (self-calibration)
    free_centroids: bool = False,
    # F-167 disambiguation: bounds для NLS-сдвига центроидов в Phase A
    # side-fit, **НЕ ID-окно** (F-167 ID-окно живёт в
    # `gamma.identification.id_window`). Передаётся в `coupled_intensity_fit`
    # верх по стеку.
    centroid_window_frac: float = 0.5,
) -> list:
    """F-118 / F-121: всегда эмиттировать связанные подгонки фиксированных
    мультиплетов для доминантной цепочки. Th-232 → M1+M2 (F-118),
    Ra-226/U-238 → U1+U2+U3 (F-121, v1.17.6). Возвращает список
    ``DeconvolutionResult`` отсортированный по нижней границе ROI.

    F-120 (v1.17.6): когда ``use_peak_image=True`` и детектор — NaI,
    подгонка использует базис peak-image (Gauss + low-E tail) вместо
    чистого гаусса. Это снижает смещение площадей на ~5-8% для линий
    с заметным низкоэнергетическим хвостом (NaI, E < 600 кэВ).

    F-127 (v1.17.7): когда ``use_T_E_model=True`` и use_peak_image=True
    и детектор — NaI, параметр хвоста T становится функцией энергии
    (``nai_tail_T_at(E_keV)``) вместо хардкода 0.7. Это улучшает
    fit-shape для широкого диапазона энергий M2 (1588–1630 кэВ).
    """
    from gamma.peaks.coupled_multiplet import (
        coupled_intensity_fit, ComponentSpec, T_TAIL_DEFAULT_NAI,
        nai_tail_T_at, H_STEP_DEFAULT_NAI,
    )
    if chain_dominance is None:
        return []
    is_th = bool(getattr(chain_dominance, "th232", False))
    is_u = bool(getattr(chain_dominance, "u238", False))
    if not (is_th or is_u):
        return []
    # F-121: выбор контракта кластеров по доминантной цепочке.
    if is_th:
        clusters_iter = TH232_FORCED_CLUSTERS
    else:
        clusters_iter = RA226_FORCED_CLUSTERS
    # F-120: базис только для NaI; для HPGe хвост пренебрежимо мал.
    is_nai = (str(detector_type or "").lower().startswith("nai"))
    enable_peak_image = bool(use_peak_image and is_nai)
    tail_T = T_TAIL_DEFAULT_NAI if enable_peak_image else 0.0
    # F-127 / v1.17.7: per-line T(E) модель vs хардкод.
    tail_T_at_cb = nai_tail_T_at if (enable_peak_image and use_T_E_model) else None
    # F-133 / v1.17.7: per-line ступенька под пиком (ГОСТ). По умолчанию
    # для NaI — H_STEP_DEFAULT_NAI=0.03; для HPGe (или если peak_image
    # отключён вообще) — 0.0.
    if h_step is None:
        h_step_eff = H_STEP_DEFAULT_NAI if enable_peak_image else 0.0
    else:
        h_step_eff = float(max(0.0, h_step))
    out = []
    counts = np.asarray(spec.counts, dtype=np.float64)
    n_ch = len(counts)
    for cluster in clusters_iter:
        try:
            ch_lo = int(round(spec.energy_to_channel(cluster["E_lo_keV"])))
            ch_hi = int(round(spec.energy_to_channel(cluster["E_hi_keV"])))
        except Exception:
            continue
        ch_lo = max(0, ch_lo)
        ch_hi = min(n_ch, ch_hi)
        if ch_hi - ch_lo < 4:
            continue
        roi_E = np.array(
            [spec.channel_to_energy(c) for c in range(ch_lo, ch_hi)],
            dtype=np.float64,
        )
        roi_counts = counts[ch_lo:ch_hi]
        comp_specs = [
            ComponentSpec(nuclide=n, E_keV=E, I_gamma_pct=I, group=g)
            for (n, E, I, g) in cluster["components"]
        ]
        try:
            # F-392 / v1.18.27 — auto-promote step_linear → step_linear_multi
            # для forced clusters с широким ROI и ≥3 intense anchors. Для
            # M1 (Th-232 750-1115) anchors {Ac-228 911, Ac-228 969} — только
            # 2, остаётся step_linear. Для U3 (1080-1830 Bi-214) anchors
            # {1120, 1238, 1378, 1764} — 4 anchor → multi.
            # ВАЖНО: forced multiplet path использует per-line step
            # (h_step_eff>0 для NaI), что подавляет глобальный β_step.
            # Auto-promote применяется только когда h_step_eff==0 (HPGe
            # или peak_image off) — иначе multi-step тоже подавлен и
            # промоут не имеет эффекта.
            roi_span_keV = (
                float(roi_E[-1] - roi_E[0])
                if len(roi_E) >= 2 else 0.0
            )
            continuum_coupled = (
                _f392_auto_select_continuum(
                    ((cs.E_keV, cs.I_gamma_pct) for cs in comp_specs),
                    "step_linear",
                    roi_e_span_keV=roi_span_keV,
                )
                if h_step_eff <= 0 else "step_linear"
            )
            coupled = coupled_intensity_fit(
                roi_E, roi_counts, comp_specs, fwhm_at_energy,
                continuum=continuum_coupled, roi_low_ch=ch_lo,
                cluster_id=cluster["id"], title=cluster["title"],
                use_peak_image=enable_peak_image,
                tail_param=tail_T,
                tail_T_at=tail_T_at_cb,
                nonlinear_refine=bool(nonlinear_refine and enable_peak_image),
                # F-133 / v1.17.7 — per-line step
                h_step=h_step_eff,
                # F-145 / v1.17.8 — Phase A side-fit (только когда NaI)
                free_centroids=bool(free_centroids and enable_peak_image),
                centroid_window_frac=float(centroid_window_frac),
            )
        except Exception:
            continue
        out.append(_coupled_to_deconv_result(coupled, fwhm_at_channel))
    return out


__all__ = [
    "MultipletComponent",
    "DeconvolutionResult",
    "deconvolve_multiplet",
    "find_multiplet_regions",
    "deconvolve_identified_multiplets",
    "apply_multiplet_deconvolution",
    # F-118
    "TH232_FORCED_CLUSTERS",
    "run_chain_forced_multiplets",
    # F-121 / v1.17.6
    "RA226_FORCED_CLUSTERS",
    # F-392 / v1.18.27 — multi-step continuum auto-selection
    "_f392_auto_select_continuum",
]
