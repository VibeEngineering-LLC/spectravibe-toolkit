"""
F-297 (v1.17.20, T-027) — Matrix method χ² minimization for multi-nuclide
                            activity deconvolution.

Когда наблюдаемые peak areas — линейная суперпозиция вкладов нескольких
нуклидов (например, Eu-152 line 121.78 keV + Cs-137 secondary + bg),
matrix method решает обратную задачу через **χ²-minimization**:

    y_i = Σ_n R_{i,n} · A_n + b_i + ε_i

где:
  • y_i — observed counts в peak i (после bg-subtraction)
  • R_{i,n} = ε(E_i) · I_n(E_i) · t — sensitivity matrix
              (peak i, nuclide n)
  • A_n — activity нуклида n (что мы хотим найти)
  • b_i — residual bg в i-м peak (если не вычтен)
  • ε_i — measurement noise, σ_i² = y_i + b_i (Poisson)

Решение:
    A = (R^T · W · R)^{-1} · R^T · W · y

где W = diag(1/σ_i²) — weight matrix.

χ² statistic:
    χ² = Σ_i (y_i - Σ_n R_{i,n} A_n)² / σ_i²
    χ²_reduced = χ² / (N_peaks - N_nuclides)

Acceptance criterion: χ²_red ≤ 1.5 — fit acceptable; > 3.0 — likely
missing nuclide или wrong library.

Зависимости
-----------
Использует только stdlib (нет numpy hard-dep) для добавочности;
работает на маленьких системах (≤ 30 peaks × 10 nuclides). Для
production пайплайна (большие системы) можно подключить numpy через
optional dependency.

References
----------
- ЛСРМ Algorithmic Foundations 2022 § 12 «Шаблонный метод»
- Gilmore & Joss «Practical Gamma-ray Spectrometry» 3rd Ed. § 9.6
- Trkov A. «Multi-channel analysis with χ² fitting» Nucl. Instrum.
  Methods A498 (2003) 425
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple


@dataclass(frozen=True)
class PeakObservation:
    """Один наблюдаемый peak: area, residual bg, energy."""
    E_keV: float
    counts: float
    counts_bg: float = 0.0     # residual bg counts in same ROI

    @property
    def variance(self) -> float:
        """Poisson variance: σ² ≈ counts + counts_bg."""
        return max(self.counts + self.counts_bg, 1.0)


@dataclass(frozen=True)
class NuclideContribution:
    """Линия n-го нуклида в i-м peak: вклад через R_{i,n}."""
    nuclide: str
    E_keV: float
    intensity_decimal: float    # I (per decay), decimal
    efficiency: float           # ε(E) at this energy
    live_time_seconds: float    # t (assumed common, но parametrize on per-line)

    @property
    def sensitivity(self) -> float:
        """R = ε · I · t — sensitivity (counts per Bq)."""
        return self.efficiency * self.intensity_decimal * self.live_time_seconds


@dataclass(frozen=True)
class MatrixMethodResult:
    """Результат χ²-минимизации."""
    activities_Bq: Dict[str, float]
    activity_uncertainties_Bq: Dict[str, float]
    chi2: float
    chi2_reduced: float
    n_peaks: int
    n_nuclides: int
    is_acceptable: bool                  # χ²_red ≤ 1.5
    needs_more_nuclides: bool            # χ²_red > 3.0
    residuals: List[float] = field(default_factory=list)


def _matrix_invert(M: List[List[float]]) -> List[List[float]]:
    """Gauss-Jordan inversion для small dense matrix (N ≤ 20).

    Не оптимизировано — для production используйте numpy. Здесь
    самодостаточная реализация чтобы не вводить hard-dep.
    """
    n = len(M)
    if any(len(row) != n for row in M):
        raise ValueError("Matrix must be square")
    # Augmented [M | I]
    A = [row[:] + [1.0 if i == j else 0.0 for j in range(n)]
         for i, row in enumerate(M)]
    for i in range(n):
        # Partial pivot
        pivot = abs(A[i][i])
        pivot_row = i
        for k in range(i + 1, n):
            if abs(A[k][i]) > pivot:
                pivot = abs(A[k][i])
                pivot_row = k
        if pivot < 1e-12:
            raise ValueError("Singular matrix (pivot ≈ 0)")
        if pivot_row != i:
            A[i], A[pivot_row] = A[pivot_row], A[i]
        # Normalize pivot row
        pv = A[i][i]
        A[i] = [x / pv for x in A[i]]
        # Eliminate other rows
        for k in range(n):
            if k == i:
                continue
            factor = A[k][i]
            A[k] = [A[k][j] - factor * A[i][j] for j in range(2 * n)]
    return [row[n:] for row in A]


def _mat_mul(A: List[List[float]], B: List[List[float]]) -> List[List[float]]:
    rows_A, cols_A = len(A), len(A[0])
    rows_B, cols_B = len(B), len(B[0])
    if cols_A != rows_B:
        raise ValueError(f"shape mismatch {rows_A}x{cols_A} × {rows_B}x{cols_B}")
    return [
        [sum(A[i][k] * B[k][j] for k in range(cols_A)) for j in range(cols_B)]
        for i in range(rows_A)
    ]


def _mat_vec(M: List[List[float]], v: List[float]) -> List[float]:
    return [sum(M[i][j] * v[j] for j in range(len(v))) for i in range(len(M))]


def solve_matrix_method(
    peaks: Sequence[PeakObservation],
    contributions: Dict[str, Sequence[NuclideContribution]],
    energy_tolerance_keV: float = 1.0,
) -> MatrixMethodResult:
    """Решить A = arg min χ²(R · A - y) с W = diag(1/σ²).

    Parameters
    ----------
    peaks : sequence of PeakObservation
        Наблюдаемые peak areas (counts) после fit.
    contributions : dict[nuclide → sequence[NuclideContribution]]
        Library: для каждого нуклида — list of line contributions.
    energy_tolerance_keV : float
        Допуск парирования peak.E_keV ↔ line.E_keV.

    Returns
    -------
    MatrixMethodResult.
    """
    nuclide_list = sorted(contributions.keys())
    N = len(peaks)
    M = len(nuclide_list)
    if N == 0 or M == 0:
        return MatrixMethodResult(
            activities_Bq={}, activity_uncertainties_Bq={},
            chi2=0.0, chi2_reduced=0.0,
            n_peaks=N, n_nuclides=M,
            is_acceptable=True, needs_more_nuclides=False,
        )
    if N < M:
        raise ValueError(
            f"Under-determined: peaks={N} < nuclides={M}"
        )

    # Build R matrix: N × M
    R = [[0.0] * M for _ in range(N)]
    for i, pk in enumerate(peaks):
        for j, nucl in enumerate(nuclide_list):
            for contrib in contributions[nucl]:
                if abs(contrib.E_keV - pk.E_keV) <= energy_tolerance_keV:
                    R[i][j] += contrib.sensitivity

    # Weight vector (1/σ²)
    w = [1.0 / pk.variance for pk in peaks]
    y = [pk.counts for pk in peaks]

    # Normal equations: (R^T W R) A = R^T W y
    RTW = [[R[i][j] * w[i] for i in range(N)] for j in range(M)]  # M×N
    RTWR = _mat_mul(RTW, R)            # M × M
    RTWy = _mat_vec(RTW, y)            # M

    try:
        inv = _matrix_invert(RTWR)
    except ValueError as e:
        raise ValueError(
            f"Sensitivity matrix non-invertible — "
            f"вероятно вырожденная library: {e}"
        )

    A = _mat_vec(inv, RTWy)            # M (activities)

    # Activity uncertainties: σ_A = sqrt(diag((R^T W R)^-1))
    sigma_A = [math.sqrt(max(inv[j][j], 0.0)) for j in range(M)]

    # χ² and residuals
    y_pred = _mat_vec(R, A)
    residuals = [y[i] - y_pred[i] for i in range(N)]
    chi2 = sum((residuals[i] ** 2) * w[i] for i in range(N))
    dof = max(N - M, 1)
    chi2_red = chi2 / dof

    return MatrixMethodResult(
        activities_Bq={nuclide_list[j]: A[j] for j in range(M)},
        activity_uncertainties_Bq={
            nuclide_list[j]: sigma_A[j] for j in range(M)
        },
        chi2=chi2,
        chi2_reduced=chi2_red,
        n_peaks=N,
        n_nuclides=M,
        is_acceptable=chi2_red <= 1.5,
        needs_more_nuclides=chi2_red > 3.0,
        residuals=residuals,
    )


__all__ = [
    "PeakObservation",
    "NuclideContribution",
    "MatrixMethodResult",
    "solve_matrix_method",
]
