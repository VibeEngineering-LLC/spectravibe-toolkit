"""
F-303 (v1.18.0, T-022b) — Quasi-template simultaneous WLS fit.

Solves the full-spectrum weighted least-squares problem:

    minimize  Σ_ch  ( S_obs[ch] − Σ_n A_n·T_n[ch] − B[ch] )² / σ²[ch]

over the unknown activities {A_n} of N library нуклидов.
T_n[ch] are per-nuclide PPP templates from F-302. B[ch] is the
externally-measured background spectrum (или нулевой если background
уже вычтен). Poisson weights w_ch = 1/max(S_obs[ch], 1).

This is a **stdlib-only** implementation (Gauss-Jordan для small N),
follows the same pattern as `matrix_method_chi2.py` (F-297).

Acceptance criterion: χ²_red ≤ 1.5 для valid fit.

References
----------
- ЛСРМ §13 (Algorithmic Foundations), quasi-template simultaneous
  fit
- ISO 11929:2019 §6.3 (Type-B uncertainty propagation)
- Press et al. "Numerical Recipes" §15.4 (general linear LSQ)
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Sequence


@dataclass
class QuasiTemplateFitResult:
    """Solution + diagnostics."""
    activities: dict[str, float]              # nuclide_id → A [Bq]
    sigma_activities: dict[str, float]        # 1-σ uncertainties [Bq]
    chi2: float
    dof: int
    chi2_red: float
    converged: bool
    nuclides_in_fit: list[str]
    notes: list[str] = field(default_factory=list)

    def is_accepted(self, chi2_red_max: float = 1.5) -> bool:
        return self.converged and self.chi2_red <= chi2_red_max


# ──────────────────────────────────────────────────────────────────
# Gauss-Jordan linear solver (skopied from matrix_method_chi2.py)
# ──────────────────────────────────────────────────────────────────

def _invert_matrix(M: list[list[float]]) -> Optional[list[list[float]]]:
    """In-place Gauss-Jordan with partial pivoting. None если singular."""
    N = len(M)
    if N == 0:
        return []
    # build augmented [M | I]
    A = [row[:] + [1.0 if j == i else 0.0 for j in range(N)]
         for i, row in enumerate(M)]
    for col in range(N):
        # partial pivot
        pivot_row = col
        max_abs = abs(A[col][col])
        for r in range(col + 1, N):
            v = abs(A[r][col])
            if v > max_abs:
                max_abs = v
                pivot_row = r
        if max_abs < 1e-15:
            return None
        if pivot_row != col:
            A[col], A[pivot_row] = A[pivot_row], A[col]
        # normalize pivot row
        piv = A[col][col]
        inv_piv = 1.0 / piv
        for j in range(2 * N):
            A[col][j] *= inv_piv
        # eliminate other rows
        for r in range(N):
            if r == col:
                continue
            factor = A[r][col]
            if factor == 0.0:
                continue
            for j in range(2 * N):
                A[r][j] -= factor * A[col][j]
    return [row[N:] for row in A]


def _mat_vec(M: list[list[float]], v: list[float]) -> list[float]:
    return [sum(M[i][j] * v[j] for j in range(len(v))) for i in range(len(M))]


# ──────────────────────────────────────────────────────────────────
# Public solver
# ──────────────────────────────────────────────────────────────────

def solve_quasi_template_fit(
    observed: Sequence[float],
    templates: Sequence,        # list of PPPTemplate
    live_time_s: float,
    background: Optional[Sequence[float]] = None,
    weights: Optional[Sequence[float]] = None,
    energy_window: Optional[tuple[int, int]] = None,
) -> QuasiTemplateFitResult:
    """Simultaneous WLS fit observed = Σ A_n · t·T_n + bg.

    Parameters
    ----------
    observed : Sequence[float]
        Observed counts per channel (gross spectrum).
    templates : Sequence[PPPTemplate]
        Per-nuclide PPP templates from F-302 (1 Bq · 1 sec per template).
    live_time_s : float
        Sample live time [s]. Used to scale templates: model[ch] = Σ A_n·t·T_n[ch].
    background : Optional[Sequence[float]]
        Background counts per channel (same n_channels as observed). Если None — нули.
    weights : Optional[Sequence[float]]
        Custom weights per channel. Если None — Poisson w = 1/max(observed, 1).
    energy_window : Optional[(ch_lo, ch_hi)]
        Restrict fit to this channel slice. Useful чтобы исключить
        мёртвые каналы LLD/ULD-области. Если None — весь спектр.

    Returns
    -------
    QuasiTemplateFitResult.
    """
    if live_time_s <= 0.0:
        raise ValueError(f"live_time_s must be >0, got {live_time_s}")
    if not templates:
        raise ValueError("Empty templates list")

    n_ch_full = len(observed)
    for i, t in enumerate(templates):
        if t.n_channels != n_ch_full:
            raise ValueError(
                f"Template[{i}] {t.nuclide_id}: n_channels={t.n_channels} "
                f"!= observed length {n_ch_full}"
            )

    if energy_window is not None:
        ch_lo, ch_hi = energy_window
        ch_lo = max(0, ch_lo)
        ch_hi = min(n_ch_full, ch_hi)
        if ch_hi <= ch_lo:
            raise ValueError(f"energy_window has zero width: {energy_window}")
    else:
        ch_lo, ch_hi = 0, n_ch_full

    N = len(templates)
    M_eff = ch_hi - ch_lo
    if M_eff < N:
        return QuasiTemplateFitResult(
            activities={t.nuclide_id: 0.0 for t in templates},
            sigma_activities={t.nuclide_id: float("inf") for t in templates},
            chi2=float("nan"), dof=0, chi2_red=float("nan"),
            converged=False,
            nuclides_in_fit=[t.nuclide_id for t in templates],
            notes=[f"Underdetermined: M={M_eff} channels < N={N} nuclides"],
        )

    if background is None:
        bg = [0.0] * n_ch_full
    else:
        if len(background) != n_ch_full:
            raise ValueError(
                f"background length {len(background)} != observed {n_ch_full}"
            )
        bg = list(background)

    if weights is None:
        # Poisson weights с гvardом против нулей
        w = [1.0 / max(float(observed[ch]), 1.0) for ch in range(n_ch_full)]
    else:
        if len(weights) != n_ch_full:
            raise ValueError(
                f"weights length {len(weights)} != observed {n_ch_full}"
            )
        w = list(weights)

    # Build normal equations: (Aᵀ·W·A)·x = Aᵀ·W·(y - bg)
    # где A[ch][n] = live_time · T_n[ch]
    # x_n = activities

    # Construct AᵀWA (N×N) и AᵀW·(y-bg) (N) только по active channels
    AtWA = [[0.0] * N for _ in range(N)]
    AtWy = [0.0] * N

    for ch in range(ch_lo, ch_hi):
        wi = w[ch]
        y_minus_bg = float(observed[ch]) - bg[ch]
        for n_i in range(N):
            ai = live_time_s * templates[n_i].counts[ch]
            AtWy[n_i] += wi * ai * y_minus_bg
            for n_j in range(n_i, N):
                aj = live_time_s * templates[n_j].counts[ch]
                AtWA[n_i][n_j] += wi * ai * aj
    # симметризация
    for i in range(N):
        for j in range(i + 1, N):
            AtWA[j][i] = AtWA[i][j]

    inv = _invert_matrix(AtWA)
    if inv is None:
        return QuasiTemplateFitResult(
            activities={t.nuclide_id: 0.0 for t in templates},
            sigma_activities={t.nuclide_id: float("inf") for t in templates},
            chi2=float("nan"), dof=M_eff - N, chi2_red=float("nan"),
            converged=False,
            nuclides_in_fit=[t.nuclide_id for t in templates],
            notes=["Singular AᵀWA — templates collinear or no signal"],
        )

    activities_vec = _mat_vec(inv, AtWy)
    sigmas_vec = [math.sqrt(max(inv[i][i], 0.0)) for i in range(N)]

    # χ²
    chi2 = 0.0
    for ch in range(ch_lo, ch_hi):
        model = bg[ch] + sum(
            activities_vec[n_i] * live_time_s * templates[n_i].counts[ch]
            for n_i in range(N)
        )
        chi2 += w[ch] * (float(observed[ch]) - model) ** 2
    dof = M_eff - N
    chi2_red = chi2 / dof if dof > 0 else float("nan")

    activities_dict: dict[str, float] = {}
    sigmas_dict: dict[str, float] = {}
    for n_i in range(N):
        nid = templates[n_i].nuclide_id
        activities_dict[nid] = activities_vec[n_i]
        sigmas_dict[nid] = sigmas_vec[n_i]

    notes: list[str] = []
    neg = [k for k, v in activities_dict.items() if v < 0.0]
    if neg:
        notes.append(
            f"Negative activities (consider non-negative LSQ): {neg}"
        )

    return QuasiTemplateFitResult(
        activities=activities_dict,
        sigma_activities=sigmas_dict,
        chi2=chi2,
        dof=dof,
        chi2_red=chi2_red,
        converged=True,
        nuclides_in_fit=[t.nuclide_id for t in templates],
        notes=notes,
    )


__all__ = [
    "QuasiTemplateFitResult",
    "solve_quasi_template_fit",
]
