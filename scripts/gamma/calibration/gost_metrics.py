"""
F-285 (v1.17.15, T-059) — GOST INL / NL metrics для energy calibration.

ГОСТ 26874-86 / МИ 1916-88 определяют две стандартные метрики
нелинейности энергетической калибровки:

  • **INL (ИНЛ)** — integral non-linearity. Maximum абсолютной
    разности между измеренными и предсказанными положениями anchor
    линий, **выраженная как доля полной шкалы**:

        INL = max|E_obs_i − E_pred_i| / E_full_scale

  • **NL (НЛ)** — differential non-linearity (point-to-point):

        NL = max|Δ(E_obs_i+1 − E_pred_i+1) − Δ(E_obs_i − E_pred_i)|
             / max|Δ_avg|

    Грубо — максимальное изменение residuals между соседними точками
    относительно среднего ширины канала.

Для NaI Gamma-1S 63×63 нормативные пороги:
  - INL ≤ 1.0 % (full scale) — годность к использованию
  - INL ≤ 0.5 %               — высшая точность
  - NL  ≤ 1.0 % per channel   — стабильность ADC

Если INL > 1.0 % — калибровка требует доработки (refit, перебор
anchor'ов, перерасчёт коэффициентов).

References
----------
- ГОСТ 26874-86 "Спектрометры энергий ионизирующих излучений"
- МИ 1916-88 "Методика поверки спектрометров"
- ЛСРМ Algorithmic Foundations §8.2
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence


@dataclass(frozen=True)
class GOSTLinearityMetrics:
    INL_pct_of_full_scale: float       # integral non-linearity, %
    NL_pct_of_full_scale: float        # point-to-point non-linearity, %
    max_residual_keV: float
    max_residual_E_keV: float          # энергия где наблюдался max residual
    n_anchors: int
    full_scale_keV: float
    accepts_inl: bool                  # INL ≤ INL_acceptable_threshold
    accepts_nl: bool                   # NL  ≤ NL_acceptable_threshold
    note: str = ""


# Каноничные пороги для Gamma-1S 63×63 (ГОСТ + ЛСРМ).
DEFAULT_INL_ACCEPTABLE_PCT = 1.0
DEFAULT_NL_ACCEPTABLE_PCT = 1.0


def compute_gost_linearity_metrics(
    anchor_energies_keV: Sequence[float],
    fitted_energies_keV: Sequence[float],
    *,
    full_scale_keV: Optional[float] = None,
    inl_acceptable_pct: float = DEFAULT_INL_ACCEPTABLE_PCT,
    nl_acceptable_pct: float = DEFAULT_NL_ACCEPTABLE_PCT,
) -> GOSTLinearityMetrics:
    """Compute ГОСТ INL / NL для калибровочного fit.

    Parameters
    ----------
    anchor_energies_keV : sequence
        Паспортные энергии калибровочных линий (keV), отсортировано
        по возрастанию.
    fitted_energies_keV : sequence
        Полученные из E(N) формулы значения (keV) на тех же каналах
        что anchor — т.е. E_fit(channel(anchor_i)).
    full_scale_keV : Optional[float]
        Полная шкала спектрометра (max paskali). Если None — берётся
        max(anchor_energies_keV).
    inl_acceptable_pct, nl_acceptable_pct : float
        Пороги для accept-флагов. Default 1.0 % (ГОСТ Gamma-1S).
    """
    if len(anchor_energies_keV) != len(fitted_energies_keV):
        raise ValueError("anchor and fitted arrays must be same length")
    n = len(anchor_energies_keV)
    if n < 2:
        return GOSTLinearityMetrics(
            INL_pct_of_full_scale=0.0,
            NL_pct_of_full_scale=0.0,
            max_residual_keV=0.0,
            max_residual_E_keV=0.0,
            n_anchors=n,
            full_scale_keV=float(full_scale_keV or 0.0),
            accepts_inl=False,
            accepts_nl=False,
            note=f"insufficient anchors ({n} < 2)",
        )

    residuals = [
        float(fitted_energies_keV[i]) - float(anchor_energies_keV[i])
        for i in range(n)
    ]
    abs_res = [abs(r) for r in residuals]
    max_res = max(abs_res)
    max_idx = abs_res.index(max_res)
    max_res_E = float(anchor_energies_keV[max_idx])

    fs = float(full_scale_keV) if full_scale_keV is not None else \
        max(float(e) for e in anchor_energies_keV)

    inl_pct = 100.0 * max_res / fs if fs > 0 else 0.0

    # NL — max(residual[i+1] - residual[i]) / fs
    if n >= 2:
        diffs = [abs(residuals[i + 1] - residuals[i]) for i in range(n - 1)]
        max_diff = max(diffs)
        nl_pct = 100.0 * max_diff / fs if fs > 0 else 0.0
    else:
        nl_pct = 0.0

    return GOSTLinearityMetrics(
        INL_pct_of_full_scale=inl_pct,
        NL_pct_of_full_scale=nl_pct,
        max_residual_keV=max_res,
        max_residual_E_keV=max_res_E,
        n_anchors=n,
        full_scale_keV=fs,
        accepts_inl=inl_pct <= inl_acceptable_pct,
        accepts_nl=nl_pct <= nl_acceptable_pct,
        note=(
            f"INL={inl_pct:.3f}% (≤{inl_acceptable_pct}%? "
            f"{'OK' if inl_pct <= inl_acceptable_pct else 'FAIL'}), "
            f"NL={nl_pct:.3f}% (≤{nl_acceptable_pct}%? "
            f"{'OK' if nl_pct <= nl_acceptable_pct else 'FAIL'})"
        ),
    )


__all__ = [
    "DEFAULT_INL_ACCEPTABLE_PCT",
    "DEFAULT_NL_ACCEPTABLE_PCT",
    "GOSTLinearityMetrics",
    "compute_gost_linearity_metrics",
]
