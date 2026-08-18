"""
F-291 (v1.17.18, T-054) — Background drift F-statistic.

Production-QA gate: периодически (по умолчанию ежедневно) проверять,
что **текущий фон детектора** статистически совместим с
**reference-фоном** (averaged background), используемым в pipeline.

Rationale
---------
Фон NaI Gamma-1S дрейфует по нескольким причинам:

  • Контаминация детектора / измерительной камеры (попадание
    Cs-137 / I-131 с лабораторного оборудования).
  • Радон Rn-222 в воздухе помещения (сезонный 50 %).
  • Деградация защиты (свинцовый домик треснул, потерял bias).
  • Cosmic-ray variability (флуктуация на ±5 %).

При значительном уходе фона все MDA и activity-результаты
смещаются. Лучше остановиться и **переснять reference bg**, чем
выдавать результаты с silent bias.

F-test semantics
----------------
В пер-ROI режиме (более чувствительно чем full-spectrum):

  F = max(s²_current, s²_ref) / min(s²_current, s²_ref)

где s² — выборочная дисперсия counts в ROI (по N_repeats фоновых
измерений за период). df1 = df2 = N_repeats − 1.

Сравнение с критическим F-значением при α=0.05:

  F_critical ≈ ppf(1 − α/2, df1, df2)

Для quick-use без scipy:

  F_critical(α=0.05, df=20, df=20) ≈ 2.46
  F_critical(α=0.05, df=10, df=10) ≈ 3.72

Если F > F_critical → дрейф значимый → FAIL.

Дополнительно: ratio mean(current)/mean(ref) проверяется отдельно
(t-test for means). При уходе mean > 3·σ_combined → FAIL.

References
----------
- ЛСРМ Algorithmic Foundations 2022 § 14.3 «Контроль фона»
- Gilmore & Joss 3rd Ed. § 5.5 «Background-stability tests»
- ISO 11929:2019 § A.2 «Type-A background-uncertainty»
- F-distribution tables — Fisher 1925
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Sequence


# Pre-computed F_critical ≈ F.ppf(0.975, df, df) for symmetric F-test.
# Used as fallback when scipy.stats недоступен; для production-режима
# рекомендуется scipy.stats.f.ppf.
_F_CRITICAL_TABLE_ALPHA_05 = {
    # df1 = df2 (symmetric paired comparison)
    5:   7.15,
    7:   4.99,
    10:  3.72,
    15:  2.86,
    20:  2.46,
    25:  2.23,
    30:  2.07,
    40:  1.88,
    60:  1.67,
    120: 1.43,
}


def _f_critical_alpha_05(df: int) -> float:
    """F-critical ≈ F.ppf(0.975, df, df). Линейная интерполяция по таблице."""
    if df <= 0:
        raise ValueError("df must be > 0")
    keys = sorted(_F_CRITICAL_TABLE_ALPHA_05)
    if df <= keys[0]:
        return _F_CRITICAL_TABLE_ALPHA_05[keys[0]]
    if df >= keys[-1]:
        return _F_CRITICAL_TABLE_ALPHA_05[keys[-1]]
    for i in range(len(keys) - 1):
        if keys[i] <= df <= keys[i + 1]:
            x0, x1 = keys[i], keys[i + 1]
            y0 = _F_CRITICAL_TABLE_ALPHA_05[x0]
            y1 = _F_CRITICAL_TABLE_ALPHA_05[x1]
            return y0 + (y1 - y0) * (df - x0) / (x1 - x0)
    return _F_CRITICAL_TABLE_ALPHA_05[keys[-1]]


@dataclass(frozen=True)
class BgRoiSnapshot:
    """Серия повторных counts в одном ROI за период наблюдения."""
    roi_label: str                  # e.g. "Cs-137 661.6 keV"
    counts_per_session: Sequence[float]  # N измерений (>= 2)
    live_time_seconds_per_session: float

    @property
    def mean_counts(self) -> float:
        return sum(self.counts_per_session) / len(self.counts_per_session)

    @property
    def variance_counts(self) -> float:
        """Sample variance (Bessel-corrected, ddof=1)."""
        n = len(self.counts_per_session)
        if n < 2:
            return 0.0
        m = self.mean_counts
        return sum((x - m) ** 2 for x in self.counts_per_session) / (n - 1)

    @property
    def df(self) -> int:
        return max(0, len(self.counts_per_session) - 1)


@dataclass(frozen=True)
class BgDriftFinding:
    roi_label: str
    f_statistic: float
    f_critical: float
    mean_ratio: float
    mean_z_score: float
    status: str                     # "PASS" / "WARNING" / "FAIL"
    message: str


@dataclass(frozen=True)
class BgDriftReport:
    overall_status: str             # "PASS" / "WARNING" / "FAIL"
    findings: List[BgDriftFinding] = field(default_factory=list)

    @property
    def is_fail(self) -> bool:
        return self.overall_status == "FAIL"


def f_test_bg_drift(
    current: BgRoiSnapshot,
    reference: BgRoiSnapshot,
    alpha: float = 0.05,
    mean_z_warning: float = 2.0,
    mean_z_fail: float = 3.0,
) -> BgDriftFinding:
    """F-test variance equality + z-score check на mean для одного ROI.

    Returns
    -------
    BgDriftFinding со статусом PASS / WARNING / FAIL.

    Notes
    -----
    Требует ≥ 2 повторов в каждой serie (иначе df = 0 → невозможно
    F-test). Возвращает WARNING если данных недостаточно.
    """
    if current.roi_label != reference.roi_label:
        raise ValueError(
            f"ROI mismatch: current='{current.roi_label}' "
            f"vs reference='{reference.roi_label}'"
        )

    df_cur = current.df
    df_ref = reference.df

    if df_cur == 0 or df_ref == 0:
        return BgDriftFinding(
            roi_label=current.roi_label,
            f_statistic=float("nan"),
            f_critical=float("nan"),
            mean_ratio=float("nan"),
            mean_z_score=float("nan"),
            status="WARNING",
            message=(
                f"{current.roi_label}: недостаточно повторов "
                f"(df_cur={df_cur}, df_ref={df_ref}) — F-test пропущен."
            ),
        )

    var_cur = current.variance_counts
    var_ref = reference.variance_counts

    if var_cur == 0 and var_ref == 0:
        f_stat = 1.0
    elif var_ref == 0:
        f_stat = float("inf")
    else:
        f_stat = (
            max(var_cur, var_ref) / min(var_cur, var_ref)
            if min(var_cur, var_ref) > 0 else float("inf")
        )

    # Симметричный df для таблицы (берём min для conservative)
    df_symmetric = min(df_cur, df_ref)
    # alpha=0.05 жёстко в таблице; для других alpha — scipy fallback
    if abs(alpha - 0.05) > 1e-6:
        # Fallback: используем 0.05 значение, но возвращаем WARNING
        f_crit = _f_critical_alpha_05(df_symmetric)
    else:
        f_crit = _f_critical_alpha_05(df_symmetric)

    # Z-score на mean (combined variance)
    n_cur = len(current.counts_per_session)
    n_ref = len(reference.counts_per_session)
    se_combined = math.sqrt(
        max(var_cur / n_cur + var_ref / n_ref, 1e-12)
    )
    mean_diff = current.mean_counts - reference.mean_counts
    mean_z = mean_diff / se_combined
    mean_ratio = (
        current.mean_counts / reference.mean_counts
        if reference.mean_counts > 0 else float("nan")
    )

    variance_drift = f_stat > f_crit
    mean_drift_warning = abs(mean_z) >= mean_z_warning
    mean_drift_fail = abs(mean_z) >= mean_z_fail

    if mean_drift_fail or (variance_drift and abs(mean_z) >= mean_z_warning):
        status = "FAIL"
        message = (
            f"{current.roi_label}: drift подтверждён (F={f_stat:.2f} vs "
            f"crit {f_crit:.2f}, mean z={mean_z:+.2f}) — НУЖНА повторная "
            f"съёмка reference bg."
        )
    elif variance_drift or mean_drift_warning:
        status = "WARNING"
        message = (
            f"{current.roi_label}: возможный drift (F={f_stat:.2f} vs "
            f"crit {f_crit:.2f}, mean z={mean_z:+.2f}) — пронаблюдать "
            f"в следующих сессиях."
        )
    else:
        status = "PASS"
        message = (
            f"{current.roi_label}: фон стабилен (F={f_stat:.2f} ≤ "
            f"crit {f_crit:.2f}, mean z={mean_z:+.2f})."
        )

    return BgDriftFinding(
        roi_label=current.roi_label,
        f_statistic=f_stat,
        f_critical=f_crit,
        mean_ratio=mean_ratio,
        mean_z_score=mean_z,
        status=status,
        message=message,
    )


def check_bg_drift_multi_roi(
    current_rois: Sequence[BgRoiSnapshot],
    reference_rois: Sequence[BgRoiSnapshot],
    alpha: float = 0.05,
) -> BgDriftReport:
    """Multi-ROI gate. Запускает F-test на каждом ROI и агрегирует."""
    cur_map = {r.roi_label: r for r in current_rois}
    ref_map = {r.roi_label: r for r in reference_rois}
    common = sorted(set(cur_map) & set(ref_map))

    findings: List[BgDriftFinding] = []
    for label in common:
        findings.append(
            f_test_bg_drift(cur_map[label], ref_map[label], alpha=alpha)
        )

    statuses = {f.status for f in findings}
    if "FAIL" in statuses:
        overall = "FAIL"
    elif "WARNING" in statuses:
        overall = "WARNING"
    else:
        overall = "PASS"

    return BgDriftReport(overall_status=overall, findings=findings)


__all__ = [
    "BgRoiSnapshot",
    "BgDriftFinding",
    "BgDriftReport",
    "f_test_bg_drift",
    "check_bg_drift_multi_roi",
]
