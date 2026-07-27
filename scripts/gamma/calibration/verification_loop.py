"""
F-286 (v1.17.15, T-060) — Calibration verification loop.

После efficiency calibration на эталонных источниках необходимо
**пересчитать активность тех же эталонов** и убедиться, что
вычисленная активность совпадает с паспортной в пределах combined
uncertainty.

Это закрывает loop **calibration → verification → acceptance**:

  1. Пользователь измеряет эталоны Cs-137, Co-60, K-40 (ОСГИ),
     получает peak areas (counts).
  2. Подгоняется ε(E) поlinom (low+high zone, F-284).
  3. Verification: для каждого эталона A_computed = peak_area / (ε·I·t)
     должно лежать в [A_passport · (1 − k·σ), A_passport · (1 + k·σ)],
     где σ — combined uncertainty (efficiency + counting + I + …),
     k=2 (95 % CL).
  4. Если выпало больше R эталонов (R = 1 на 5 — порог ЛСРМ) →
     калибровка отклонена.

Этот модуль реализует Step 3-4 как pure-function helper.

References
----------
- ЛСРМ Algorithmic Foundations §8.4.5 "Проверка калибровки"
- ГОСТ Р 51086-97 §6.3 «Поверка спектрометра»
- ISO 17025:2017 §7.6 «Контроль качества результатов»
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Sequence


@dataclass(frozen=True)
class StandardComparison:
    """Сравнение паспортной и вычисленной активности для одного эталона."""
    nuclide: str
    line_keV: float
    A_passport_Bq: float
    A_passport_unc_Bq: float
    A_computed_Bq: float
    A_computed_unc_Bq: float
    combined_unc_Bq: float            # √(σ_passport² + σ_computed²)
    deviation_sigma: float            # |A_comp − A_pass| / combined_unc
    within_2_sigma: bool
    note: str = ""


@dataclass(frozen=True)
class CalibrationVerificationOutcome:
    """Результат проверки калибровки на наборе эталонов."""
    comparisons: List[StandardComparison]
    n_total: int
    n_within_2_sigma: int
    n_failing: int
    fraction_failing: float
    acceptable: bool                  # fraction_failing ≤ FAIL_TOLERANCE
    note: str = ""


# Если ≥ FAIL_TOLERANCE_FRACTION эталонов выпало за 2σ — калибровка отклонена.
FAIL_TOLERANCE_FRACTION = 0.20    # ≥20% выпавших → reject (ЛСРМ §8.4.5)


def verify_calibration_against_standards(
    standards: Sequence[StandardComparison],
    *,
    fail_tolerance_fraction: float = FAIL_TOLERANCE_FRACTION,
) -> CalibrationVerificationOutcome:
    """Проверить выборку эталонов против их паспортных активностей.

    Все StandardComparison уже должны быть посчитаны caller'ом
    (с правильным combined_unc); эта функция только агрегирует
    результат.
    """
    n = len(standards)
    if n == 0:
        return CalibrationVerificationOutcome(
            comparisons=[],
            n_total=0,
            n_within_2_sigma=0,
            n_failing=0,
            fraction_failing=0.0,
            acceptable=False,
            note="No standards to verify against",
        )
    within = sum(1 for s in standards if s.within_2_sigma)
    failing = n - within
    frac_fail = failing / n
    acceptable = frac_fail <= fail_tolerance_fraction
    return CalibrationVerificationOutcome(
        comparisons=list(standards),
        n_total=n,
        n_within_2_sigma=within,
        n_failing=failing,
        fraction_failing=frac_fail,
        acceptable=acceptable,
        note=(
            f"F-286: {within}/{n} эталонов внутри 2σ, "
            f"{failing}/{n} fail; fraction_failing={frac_fail:.0%} "
            f"(порог ≤{fail_tolerance_fraction:.0%}) — "
            f"{'ACCEPT' if acceptable else 'REJECT'} калибровку"
        ),
    )


def make_standard_comparison(
    *,
    nuclide: str,
    line_keV: float,
    A_passport_Bq: float,
    A_passport_unc_Bq: float,
    A_computed_Bq: float,
    A_computed_unc_Bq: float,
) -> StandardComparison:
    """Создать StandardComparison с правильным combined_unc и within_2σ."""
    combined = (float(A_passport_unc_Bq) ** 2
                + float(A_computed_unc_Bq) ** 2) ** 0.5
    if combined > 0:
        dev = abs(A_computed_Bq - A_passport_Bq) / combined
    else:
        dev = float("inf") if A_computed_Bq != A_passport_Bq else 0.0
    return StandardComparison(
        nuclide=str(nuclide),
        line_keV=float(line_keV),
        A_passport_Bq=float(A_passport_Bq),
        A_passport_unc_Bq=float(A_passport_unc_Bq),
        A_computed_Bq=float(A_computed_Bq),
        A_computed_unc_Bq=float(A_computed_unc_Bq),
        combined_unc_Bq=float(combined),
        deviation_sigma=float(dev),
        within_2_sigma=bool(dev <= 2.0),
    )


__all__ = [
    "FAIL_TOLERANCE_FRACTION",
    "StandardComparison",
    "CalibrationVerificationOutcome",
    "verify_calibration_against_standards",
    "make_standard_comparison",
]
