"""
F-289 (v1.17.18, T-083) — Spectrometer compliance gate.

Production-QA gate: BEFORE accepting an analysis result, assert that
the **measured detector characteristics** (resolution, channel count,
efficiency at key energies) match the **certified specification** of
the spectrometer.

Rationale
---------
Gamma-1S NaI 63×63 is delivered with a factory certificate listing:

  • FWHM @ 662 keV (Cs-137) ≤ 8.0 % (typical 6.5–7.0 %)
  • Channels: 1024 (USB) или 2048
  • Integral non-linearity (ИНЛ) ≤ 1.0 %
  • Bias drift ≤ 0.5 % за 8 h

Если **измеренная** разрешающая способность ушла далеко от паспорта —
PMT деградировал, NaI кристалл потрескался, либо HV сбит. Любой
дальнейший анализ (peak search, ID, MDA) даст misleading результаты.
Лучше остановить pipeline с ясным сообщением, чем выдать invalid
отчёт.

Gate semantics
--------------
``check_spectrometer_compliance(measured, spec, tolerance_frac=0.10)``
возвращает ``ComplianceReport`` с status:

  • ``"PASS"`` — все измерения в пределах ``spec · (1 ± tolerance_frac)``
  • ``"WARNING"`` — одно отклонение в пределах 2·tolerance_frac
  • ``"FAIL"`` — любое отклонение > 2·tolerance_frac или критическое поле

FAIL означает: НЕ продолжать analysis, потребовать service-check.
WARNING означает: продолжить, но добавить flag в report.

References
----------
- ЛСРМ Algorithmic Foundations 2022 § 14.1 «Контроль работоспособности
  спектрометра»
- ГОСТ Р 51086-97 § 6 «Поверка»
- ISO 17025:2017 § 7.7 «Контроль качества результатов»
- Gamma-1S USB Certificate (factory spec sheet)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional


# Factory-spec defaults для Gamma-1S NaI 63×63 USB (LSRM/Aspect).
# Источник: certificate шаблон ЛСРМ + ОЭБ Колибри/Gamma-1S tech-passport.
GAMMA_1C_NAI_63X63_DEFAULT_SPEC = {
    "fwhm_at_662_keV_pct_max": 8.0,        # certificate ≤ 8 %
    "fwhm_at_662_keV_pct_typical": 7.0,    # typical Gamma-1S 6.5–7.0
    "channels_min": 1024,                  # USB вариант
    "inl_pct_max": 1.0,                    # integral non-linearity
    "bias_drift_8h_pct_max": 0.5,
}


@dataclass(frozen=True)
class SpectrometerSpec:
    """Сертифицированная спецификация спектрометра (фабричный паспорт)."""

    detector_class: str                 # "NaI" / "HPGe" / "LaBr" / "CeBr" / "CsI"
    model: str                          # e.g. "Gamma-1S NaI 63×63"
    fwhm_at_662_keV_pct_max: float      # ≤ это значение
    channels_min: int                   # минимально допустимое число каналов
    inl_pct_max: float = 1.0            # integral non-linearity
    efficiency_anchors: Dict[float, float] = field(default_factory=dict)
    # {E_keV: ε_passport} — опционально, для cross-check с efficiency-cal

    @classmethod
    def gamma_1c_nai_63x63(cls) -> "SpectrometerSpec":
        """Default factory spec for Gamma-1S NaI 63×63 USB."""
        return cls(
            detector_class="NaI",
            model="Gamma-1S NaI 63×63 USB",
            fwhm_at_662_keV_pct_max=GAMMA_1C_NAI_63X63_DEFAULT_SPEC["fwhm_at_662_keV_pct_max"],
            channels_min=GAMMA_1C_NAI_63X63_DEFAULT_SPEC["channels_min"],
            inl_pct_max=GAMMA_1C_NAI_63X63_DEFAULT_SPEC["inl_pct_max"],
        )


@dataclass(frozen=True)
class SpectrometerMeasurement:
    """Измеренные характеристики спектрометра (после fwhm-fit / energy-cal)."""

    fwhm_at_662_keV_pct: Optional[float] = None
    channels_actual: Optional[int] = None
    inl_pct: Optional[float] = None         # из gost_metrics.compute_gost_linearity_metrics
    bias_drift_8h_pct: Optional[float] = None


@dataclass(frozen=True)
class ComplianceFinding:
    field_name: str
    status: str                # "PASS" / "WARNING" / "FAIL"
    measured_value: Optional[float]
    spec_max: Optional[float]
    deviation_pct: Optional[float]
    message: str


@dataclass(frozen=True)
class ComplianceReport:
    overall_status: str              # "PASS" / "WARNING" / "FAIL"
    findings: List[ComplianceFinding] = field(default_factory=list)

    @property
    def is_fail(self) -> bool:
        return self.overall_status == "FAIL"

    @property
    def is_pass(self) -> bool:
        return self.overall_status == "PASS"


def _check_max(
    field_name: str,
    measured: Optional[float],
    spec_max: float,
    tolerance_frac: float,
) -> ComplianceFinding:
    """Универсальный «≤ spec» check с трёхуровневым статусом."""
    if measured is None:
        return ComplianceFinding(
            field_name=field_name,
            status="WARNING",
            measured_value=None,
            spec_max=spec_max,
            deviation_pct=None,
            message=f"{field_name}: измерение отсутствует — невозможно проверить compliance.",
        )

    if spec_max <= 0:
        return ComplianceFinding(
            field_name=field_name,
            status="WARNING",
            measured_value=measured,
            spec_max=spec_max,
            deviation_pct=None,
            message=f"{field_name}: spec_max ≤ 0, проверка пропущена.",
        )

    soft_limit = spec_max * (1.0 + tolerance_frac)
    hard_limit = spec_max * (1.0 + 2.0 * tolerance_frac)
    dev_pct = 100.0 * (measured - spec_max) / spec_max

    if measured <= spec_max:
        return ComplianceFinding(
            field_name=field_name,
            status="PASS",
            measured_value=measured,
            spec_max=spec_max,
            deviation_pct=dev_pct,
            message=f"{field_name}: {measured:.3g} ≤ spec {spec_max:.3g} → OK.",
        )
    if measured <= soft_limit:
        return ComplianceFinding(
            field_name=field_name,
            status="WARNING",
            measured_value=measured,
            spec_max=spec_max,
            deviation_pct=dev_pct,
            message=(
                f"{field_name}: {measured:.3g} > spec {spec_max:.3g} "
                f"на {dev_pct:+.1f}% (≤ tolerance {tolerance_frac*100:.0f}%)."
            ),
        )
    if measured <= hard_limit:
        return ComplianceFinding(
            field_name=field_name,
            status="WARNING",
            measured_value=measured,
            spec_max=spec_max,
            deviation_pct=dev_pct,
            message=(
                f"{field_name}: {measured:.3g} > spec {spec_max:.3g} "
                f"на {dev_pct:+.1f}% — рекомендуется service check."
            ),
        )
    return ComplianceFinding(
        field_name=field_name,
        status="FAIL",
        measured_value=measured,
        spec_max=spec_max,
        deviation_pct=dev_pct,
        message=(
            f"{field_name}: {measured:.3g} > spec {spec_max:.3g} "
            f"на {dev_pct:+.1f}% (> 2·tolerance) — анализ ОТКЛОНЁН, требуется service."
        ),
    )


def check_spectrometer_compliance(
    measurement: SpectrometerMeasurement,
    spec: SpectrometerSpec,
    tolerance_frac: float = 0.10,
) -> ComplianceReport:
    """Главная gate-функция: возвращает report со статусами всех проверок.

    Parameters
    ----------
    measurement : SpectrometerMeasurement
        Что фактически намерено в сегодняшней сессии / при последней калибровке.
    spec : SpectrometerSpec
        Сертифицированная факторская спецификация (создавайте через
        ``SpectrometerSpec.gamma_1c_nai_63x63()`` для default).
    tolerance_frac : float
        Допустимое относительное отклонение (default 10 % — практика ЛСРМ).
        WARNING при ≤ 2·tolerance_frac, FAIL при > 2·tolerance_frac.

    Returns
    -------
    ComplianceReport. Используйте ``.is_fail`` чтобы остановить pipeline.
    """
    if tolerance_frac < 0:
        raise ValueError("tolerance_frac must be ≥ 0")

    findings: List[ComplianceFinding] = []

    # FWHM @ 662 keV
    findings.append(
        _check_max(
            "fwhm_at_662_keV_pct",
            measurement.fwhm_at_662_keV_pct,
            spec.fwhm_at_662_keV_pct_max,
            tolerance_frac,
        )
    )

    # Channels — отдельная логика «≥ min»
    if measurement.channels_actual is None:
        findings.append(
            ComplianceFinding(
                field_name="channels",
                status="WARNING",
                measured_value=None,
                spec_max=float(spec.channels_min),
                deviation_pct=None,
                message="channels: число каналов не передано — пропущено.",
            )
        )
    elif measurement.channels_actual >= spec.channels_min:
        findings.append(
            ComplianceFinding(
                field_name="channels",
                status="PASS",
                measured_value=float(measurement.channels_actual),
                spec_max=float(spec.channels_min),
                deviation_pct=0.0,
                message=(
                    f"channels: {measurement.channels_actual} ≥ min "
                    f"{spec.channels_min} → OK."
                ),
            )
        )
    else:
        findings.append(
            ComplianceFinding(
                field_name="channels",
                status="FAIL",
                measured_value=float(measurement.channels_actual),
                spec_max=float(spec.channels_min),
                deviation_pct=100.0 * (
                    measurement.channels_actual - spec.channels_min
                ) / spec.channels_min,
                message=(
                    f"channels: {measurement.channels_actual} < min "
                    f"{spec.channels_min} — недостаточное разрешение спектра."
                ),
            )
        )

    # ИНЛ
    findings.append(
        _check_max(
            "inl_pct", measurement.inl_pct, spec.inl_pct_max, tolerance_frac,
        )
    )

    # Drift (если передан)
    if measurement.bias_drift_8h_pct is not None:
        findings.append(
            _check_max(
                "bias_drift_8h_pct",
                measurement.bias_drift_8h_pct,
                GAMMA_1C_NAI_63X63_DEFAULT_SPEC["bias_drift_8h_pct_max"],
                tolerance_frac,
            )
        )

    statuses = {f.status for f in findings}
    if "FAIL" in statuses:
        overall = "FAIL"
    elif "WARNING" in statuses:
        overall = "WARNING"
    else:
        overall = "PASS"

    return ComplianceReport(overall_status=overall, findings=findings)


__all__ = [
    "GAMMA_1C_NAI_63X63_DEFAULT_SPEC",
    "SpectrometerSpec",
    "SpectrometerMeasurement",
    "ComplianceFinding",
    "ComplianceReport",
    "check_spectrometer_compliance",
]
