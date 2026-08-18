"""
F-292 (v1.17.18, T-053) — Sensitivity (efficiency) drift quarterly check.

Production-QA gate: ежеквартальная (default 90 дней) проверка
**отклонения текущей эффективности** ε(E) от reference-эфф-калибровки,
проведённой ранее. При уходе > 5 % per LSRM § 14.4 — переснять
calibration.

Rationale
---------
Эффективность ε(E) спектрометра дрейфует медленно, но систематически:

  • PMT gain — снижение усиления с возрастом фотокатода (NaI).
  • Optical coupling — отслоение силикона / усыхание optical grease.
  • NaI hygroscopy — поверхностная гигроскопичность ухудшает
    светосбор у низкоэнергетических лиц (< 100 keV).
  • Электронный тракт — деградация preamp / ADC.

Без периодической поверки calibration «уезжает» на 1-2 %/quarter,
что соответствует ~5 %/year — порог переснятия per ЛСРМ.

Gate semantics
--------------
Для каждого reference-energy E:

  drift_pct(E) = 100 · |ε_current(E) − ε_reference(E)| / ε_reference(E)

Status:

  • drift_pct ≤ 2.0 % AND ≤ 30 дней — PASS
  • drift_pct ≤ 5.0 % AND ≤ 90 дней — WARNING (наблюдать)
  • drift_pct > 5.0 % OR > 90 дней — FAIL (переснять calibration)

Дополнительно проверяется **shape consistency**: уход
монотонный по E указывает на geometry shift, а random — на counting
fluctuation. Используется Pearson r между ratio'ами.

References
----------
- ЛСРМ Algorithmic Foundations 2022 § 14.4 «Проверка
  эффективностной калибровки»
- ГОСТ Р 51086-97 § 6.4 «Периодическая поверка»
- IAEA TECDOC-1599 «QA in radioactivity testing»
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple


# Per-line drift thresholds, % (per ЛСРМ § 14.4)
DEFAULT_DRIFT_WARNING_PCT = 2.0
DEFAULT_DRIFT_FAIL_PCT = 5.0

# Calendar windows, days
DEFAULT_INTERVAL_OK_DAYS = 30
DEFAULT_INTERVAL_WARN_DAYS = 90    # quarterly check
DEFAULT_INTERVAL_FAIL_DAYS = 180


@dataclass(frozen=True)
class EfficiencyAnchor:
    """Точка ε(E) в reference или current calibration."""
    nuclide: str
    E_keV: float
    epsilon: float                 # decimal (0..1)
    sigma_epsilon: Optional[float] = None  # абсолютная (1-σ)


@dataclass(frozen=True)
class SensitivityDriftFinding:
    E_keV: float
    nuclide: str
    eps_current: float
    eps_reference: float
    drift_pct: float
    status: str                    # "PASS" / "WARNING" / "FAIL"
    message: str


@dataclass(frozen=True)
class SensitivityDriftReport:
    overall_status: str            # "PASS" / "WARNING" / "FAIL"
    days_since_reference: int
    monotonic_shift: Optional[bool]  # True если drift_pct монотонен по E
    pearson_r_with_E: Optional[float]
    findings: List[SensitivityDriftFinding] = field(default_factory=list)
    next_check_recommendation: str = ""

    @property
    def is_fail(self) -> bool:
        return self.overall_status == "FAIL"


def _pearson_r(xs: Sequence[float], ys: Sequence[float]) -> Optional[float]:
    """Plain Pearson correlation; None if degenerate."""
    n = len(xs)
    if n < 3 or n != len(ys):
        return None
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if dx == 0 or dy == 0:
        return None
    return num / (dx * dy)


def check_sensitivity_drift(
    current_anchors: Sequence[EfficiencyAnchor],
    reference_anchors: Sequence[EfficiencyAnchor],
    days_since_reference: int,
    drift_warning_pct: float = DEFAULT_DRIFT_WARNING_PCT,
    drift_fail_pct: float = DEFAULT_DRIFT_FAIL_PCT,
    interval_ok_days: int = DEFAULT_INTERVAL_OK_DAYS,
    interval_warn_days: int = DEFAULT_INTERVAL_WARN_DAYS,
    interval_fail_days: int = DEFAULT_INTERVAL_FAIL_DAYS,
    energy_tolerance_keV: float = 5.0,
) -> SensitivityDriftReport:
    """Главный gate ε-drift.

    Парует current и reference anchors по близости E_keV (± tolerance),
    считает per-line drift_pct, агрегирует overall_status.

    Parameters
    ----------
    days_since_reference : int
        Число дней с момента съёма reference calibration.
    drift_warning_pct, drift_fail_pct : float
        Порог в процентах. Default 2 / 5 (per ЛСРМ).
    interval_*_days : int
        Календарные окна. Default 30 / 90 / 180.
    energy_tolerance_keV : float
        Допуск парирования по E (default ± 5 keV).
    """
    if drift_warning_pct < 0 or drift_fail_pct < 0:
        raise ValueError("thresholds must be ≥ 0")
    if drift_fail_pct < drift_warning_pct:
        raise ValueError("drift_fail_pct must be ≥ drift_warning_pct")

    findings: List[SensitivityDriftFinding] = []
    drift_pcts_by_E: List[Tuple[float, float]] = []

    ref_sorted = sorted(reference_anchors, key=lambda a: a.E_keV)
    for cur in current_anchors:
        # Найти ближайший reference anchor по E
        best = None
        best_dE = float("inf")
        for ref in ref_sorted:
            dE = abs(ref.E_keV - cur.E_keV)
            if dE < best_dE:
                best_dE = dE
                best = ref
        if best is None or best_dE > energy_tolerance_keV:
            continue
        if best.epsilon <= 0:
            continue

        drift_pct = 100.0 * abs(cur.epsilon - best.epsilon) / best.epsilon

        if drift_pct <= drift_warning_pct:
            status = "PASS"
            msg = (
                f"E={cur.E_keV:.1f} keV ({cur.nuclide}): drift "
                f"{drift_pct:.2f}% ≤ warn {drift_warning_pct:.1f}% — OK."
            )
        elif drift_pct <= drift_fail_pct:
            status = "WARNING"
            msg = (
                f"E={cur.E_keV:.1f} keV ({cur.nuclide}): drift "
                f"{drift_pct:.2f}% в зоне {drift_warning_pct:.1f}–"
                f"{drift_fail_pct:.1f}% — наблюдать."
            )
        else:
            status = "FAIL"
            msg = (
                f"E={cur.E_keV:.1f} keV ({cur.nuclide}): drift "
                f"{drift_pct:.2f}% > {drift_fail_pct:.1f}% — переснять "
                f"efficiency calibration."
            )

        findings.append(
            SensitivityDriftFinding(
                E_keV=cur.E_keV,
                nuclide=cur.nuclide,
                eps_current=cur.epsilon,
                eps_reference=best.epsilon,
                drift_pct=drift_pct,
                status=status,
                message=msg,
            )
        )
        drift_pcts_by_E.append((cur.E_keV, drift_pct))

    statuses = {f.status for f in findings}

    # Calendar-based override
    if days_since_reference > interval_fail_days:
        cal_status = "FAIL"
    elif days_since_reference > interval_warn_days:
        cal_status = "WARNING"
    elif days_since_reference > interval_ok_days:
        cal_status = "WARNING"
    else:
        cal_status = "PASS"

    if "FAIL" in statuses or cal_status == "FAIL":
        overall = "FAIL"
    elif "WARNING" in statuses or cal_status == "WARNING":
        overall = "WARNING"
    else:
        overall = "PASS"

    # Pearson r: drift_pct vs E
    Es = [t[0] for t in drift_pcts_by_E]
    dps = [t[1] for t in drift_pcts_by_E]
    pr = _pearson_r(Es, dps)
    monotonic_shift = None
    if pr is not None:
        monotonic_shift = abs(pr) > 0.7

    # Next-check recommendation
    if overall == "FAIL":
        rec = "Переснять efficiency calibration НЕМЕДЛЕННО."
    elif overall == "WARNING":
        rec = (
            f"Запланировать переснятие в течение "
            f"{max(0, interval_warn_days - days_since_reference)} дней."
        )
    else:
        rec = (
            f"Следующая проверка через "
            f"{max(0, interval_warn_days - days_since_reference)} дней "
            f"(quarterly)."
        )

    return SensitivityDriftReport(
        overall_status=overall,
        days_since_reference=days_since_reference,
        monotonic_shift=monotonic_shift,
        pearson_r_with_E=pr,
        findings=findings,
        next_check_recommendation=rec,
    )


__all__ = [
    "DEFAULT_DRIFT_WARNING_PCT",
    "DEFAULT_DRIFT_FAIL_PCT",
    "DEFAULT_INTERVAL_OK_DAYS",
    "DEFAULT_INTERVAL_WARN_DAYS",
    "DEFAULT_INTERVAL_FAIL_DAYS",
    "EfficiencyAnchor",
    "SensitivityDriftFinding",
    "SensitivityDriftReport",
    "check_sensitivity_drift",
]
