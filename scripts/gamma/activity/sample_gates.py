"""
F-290 (v1.17.18, T-084) — Sample density / fill-level gate.

Production-QA gate: BEFORE применять calibration к sample,
проверить, что **геометрия и плотность образца** находятся в пределах
валидности используемой efficiency calibration.

Rationale
---------
Efficiency calibration ε(E) на Gamma-1S Marinelli 0.5 L снимается
с эталонами **фиксированной геометрии и плотности**:

  • Marinelli 0.5 L: fill-height **80–100 mm** от дна стакана.
  • Плотность ρ: water-equivalent **0.95–1.15 g/cm³** (ОЭБ ЛСРМ).
  • При уходе из этого диапазона ε(E) меняется до ±15 % (low-E)
    из-за **самопоглощения** в образце (T-002, T-024).

Если оператор положил sample с ρ=2.5 г/см³ (грунт сухой плотный)
или fill-height=50 mm (недозагрузил) — calibration **неприменима**.
Любой результат A_Bq будет смещён в одну сторону, причём
систематически (не статистически).

Gate semantics
--------------
``check_sample_compliance(sample, calib_range)`` возвращает
``SampleComplianceReport`` со статусом и явным numeric flag
о необходимости либо:

  • применить **поправку самопоглощения** (Marinelli matrix method),
  • либо отклонить sample и потребовать перепаковку.

Этот модуль НЕ выполняет поправку — он только flag-гейт.
Поправка отдельным релизом (v1.17.16 deferred).

References
----------
- ЛСРМ Algorithmic Foundations 2022 § 9.2 «Самопоглощение
  в Marinelli-геометрии»
- ОЭБ Колибри/Gamma-1S tech-passport: Marinelli 0.5 L spec
- ISO 18589-3:2015 «Soil radioactivity measurement»
- Gilmore & Joss 3rd Ed. § 7.3 «Sample-source effects»
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional


# Default validity range для эталонной геометрии Marinelli 0.5 L на Gamma-1S.
# Источник: ОЭБ ЛСРМ для эталонов ОИСН-16 / IB-08 (water-equivalent matrix).
MARINELLI_0_5L_DEFAULT_RANGE = {
    "fill_height_mm_min": 80.0,
    "fill_height_mm_max": 100.0,
    "density_g_cm3_min": 0.95,
    "density_g_cm3_max": 1.15,
    # Расширенный диапазон (требует self-absorption correction):
    "fill_height_mm_extrap_min": 60.0,
    "fill_height_mm_extrap_max": 110.0,
    "density_g_cm3_extrap_min": 0.50,
    "density_g_cm3_extrap_max": 2.50,
}


@dataclass(frozen=True)
class SampleGeometry:
    """Фактическая геометрия и плотность измеряемого образца."""
    container_type: str               # "Marinelli_0.5L" / "Marinelli_1L" / "vial_20ml" / ...
    fill_height_mm: Optional[float] = None
    density_g_cm3: Optional[float] = None
    mass_kg: Optional[float] = None
    matrix_label: Optional[str] = None   # "water" / "soil_dry" / "ash" / ...


@dataclass(frozen=True)
class CalibrationValidityRange:
    """Диапазон, в котором efficiency-cal считается валидной без поправок."""
    container_type: str
    fill_height_mm_min: float
    fill_height_mm_max: float
    density_g_cm3_min: float
    density_g_cm3_max: float
    fill_height_mm_extrap_min: float
    fill_height_mm_extrap_max: float
    density_g_cm3_extrap_min: float
    density_g_cm3_extrap_max: float

    @classmethod
    def marinelli_0_5l_default(cls) -> "CalibrationValidityRange":
        """Default Marinelli 0.5 L range (ОЭБ ЛСРМ)."""
        return cls(
            container_type="Marinelli_0.5L",
            **MARINELLI_0_5L_DEFAULT_RANGE,
        )


@dataclass(frozen=True)
class SampleComplianceFinding:
    field_name: str
    status: str                    # "PASS" / "WARNING" / "EXTRAPOLATION" / "FAIL"
    measured_value: Optional[float]
    range_min: Optional[float]
    range_max: Optional[float]
    message: str


@dataclass(frozen=True)
class SampleComplianceReport:
    overall_status: str            # "PASS" / "WARNING" / "EXTRAPOLATION" / "FAIL"
    requires_self_absorption_correction: bool
    findings: List[SampleComplianceFinding] = field(default_factory=list)

    @property
    def is_fail(self) -> bool:
        return self.overall_status == "FAIL"


def _check_range(
    field_name: str,
    measured: Optional[float],
    nominal_min: float,
    nominal_max: float,
    extrap_min: float,
    extrap_max: float,
) -> SampleComplianceFinding:
    """4-уровневая проверка: PASS / WARNING (близко к границе)
    / EXTRAPOLATION (в расширенном диапазоне) / FAIL (вне всего)."""
    if measured is None:
        return SampleComplianceFinding(
            field_name=field_name,
            status="WARNING",
            measured_value=None,
            range_min=nominal_min,
            range_max=nominal_max,
            message=f"{field_name}: значение не указано — невозможно проверить.",
        )

    if nominal_min <= measured <= nominal_max:
        # Soft-warning при близости к границе (10 % от ширины диапазона)
        band = 0.10 * (nominal_max - nominal_min)
        if measured < nominal_min + band or measured > nominal_max - band:
            return SampleComplianceFinding(
                field_name=field_name,
                status="WARNING",
                measured_value=measured,
                range_min=nominal_min,
                range_max=nominal_max,
                message=(
                    f"{field_name}={measured:.3g} в нормальном диапазоне "
                    f"[{nominal_min:.3g}; {nominal_max:.3g}], но близко к границе."
                ),
            )
        return SampleComplianceFinding(
            field_name=field_name,
            status="PASS",
            measured_value=measured,
            range_min=nominal_min,
            range_max=nominal_max,
            message=(
                f"{field_name}={measured:.3g} в норме "
                f"[{nominal_min:.3g}; {nominal_max:.3g}]."
            ),
        )

    if extrap_min <= measured <= extrap_max:
        return SampleComplianceFinding(
            field_name=field_name,
            status="EXTRAPOLATION",
            measured_value=measured,
            range_min=nominal_min,
            range_max=nominal_max,
            message=(
                f"{field_name}={measured:.3g} вне нормального диапазона "
                f"[{nominal_min:.3g}; {nominal_max:.3g}], но в расширенном "
                f"[{extrap_min:.3g}; {extrap_max:.3g}] — требуется "
                f"self-absorption correction."
            ),
        )

    return SampleComplianceFinding(
        field_name=field_name,
        status="FAIL",
        measured_value=measured,
        range_min=nominal_min,
        range_max=nominal_max,
        message=(
            f"{field_name}={measured:.3g} вне допустимого диапазона "
            f"[{extrap_min:.3g}; {extrap_max:.3g}] — sample НЕ ПОДХОДИТ "
            f"для этой калибровки, перепакуйте или возьмите другой эталон."
        ),
    )


def check_sample_compliance(
    sample: SampleGeometry,
    calib_range: CalibrationValidityRange,
) -> SampleComplianceReport:
    """Главный gate — проверяет геометрию sample против validity range.

    Returns
    -------
    SampleComplianceReport. Используйте ``.is_fail`` чтобы остановить
    pipeline, ``.requires_self_absorption_correction`` чтобы поднять
    флаг для последующего matrix-method (T-027 deferred v1.17.16).
    """
    findings: List[SampleComplianceFinding] = []

    # Container-type guard
    if sample.container_type != calib_range.container_type:
        findings.append(
            SampleComplianceFinding(
                field_name="container_type",
                status="FAIL",
                measured_value=None,
                range_min=None,
                range_max=None,
                message=(
                    f"container_type: sample='{sample.container_type}' "
                    f"vs calibration='{calib_range.container_type}' — "
                    f"геометрии не совпадают."
                ),
            )
        )
    else:
        findings.append(
            SampleComplianceFinding(
                field_name="container_type",
                status="PASS",
                measured_value=None,
                range_min=None,
                range_max=None,
                message=f"container_type='{sample.container_type}' OK.",
            )
        )

    findings.append(
        _check_range(
            "fill_height_mm",
            sample.fill_height_mm,
            calib_range.fill_height_mm_min,
            calib_range.fill_height_mm_max,
            calib_range.fill_height_mm_extrap_min,
            calib_range.fill_height_mm_extrap_max,
        )
    )

    findings.append(
        _check_range(
            "density_g_cm3",
            sample.density_g_cm3,
            calib_range.density_g_cm3_min,
            calib_range.density_g_cm3_max,
            calib_range.density_g_cm3_extrap_min,
            calib_range.density_g_cm3_extrap_max,
        )
    )

    statuses = {f.status for f in findings}
    if "FAIL" in statuses:
        overall = "FAIL"
    elif "EXTRAPOLATION" in statuses:
        overall = "EXTRAPOLATION"
    elif "WARNING" in statuses:
        overall = "WARNING"
    else:
        overall = "PASS"

    requires_corr = "EXTRAPOLATION" in statuses

    return SampleComplianceReport(
        overall_status=overall,
        requires_self_absorption_correction=requires_corr,
        findings=findings,
    )


__all__ = [
    "MARINELLI_0_5L_DEFAULT_RANGE",
    "SampleGeometry",
    "CalibrationValidityRange",
    "SampleComplianceFinding",
    "SampleComplianceReport",
    "check_sample_compliance",
]
