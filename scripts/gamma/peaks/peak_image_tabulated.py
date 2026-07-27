"""
F-299 (v1.17.21, T-021a) — Tabulated peak image (ЛСРМ §8.4.1).

ЛСРМ алгоритм SpectraLine использует **табличное представление**
peak shape: вместо аналитической функции (Gauss/Voigt/peak_image
из F-90), peak хранится как набор anchor-точек по энергии, в каждой
из которых записаны параметры (FWHM, asymmetry, tail amplitude, step).

Преимущества:
  • Может описать **произвольную форму пика** (например, multi-Gaussian
    artifacts от electronic crosstalk, или asymmetric peaks на сцинтилляторах).
  • Естественно стыкуется с .cpt файлами (LSRM template format).
  • Smoothed-spline между anchor'ами даёт C² continuous shape.

Использование:
  • Calibrate-once на эталоне (Cs-137, Eu-152, Co-60) → store as
    TabulatedPeakImage → save as .cpt.
  • Apply-many на sample спектр: load .cpt → for each peak interpolate
    shape parameters → use as fit template.

Канонические anchor-параметры (per ЛСРМ §8.4.1)
-----------------------------------------------
В каждой anchor-точке E_i хранятся:
  • fwhm_keV(E_i)            — резолюция в точке
  • shape_amplitude_norm     — нормировка пик-нула (всегда 1.0)
  • tail_fraction(E_i)       — доля low-E tail компоненты
  • tail_slope(E_i)          — экспоненциальный спад tail (1/keV)
  • step_height(E_i)         — Compton step height (доля от FEP)
  • asymmetry(E_i)           — параметр Voigt-asymmetry

Для NaI обычно tail_fraction ≈ 0.02..0.05, step_height ≈ 0.03..0.10
(Будыка 7.29). Для HPGe: tail_fraction ≈ 0.001, step_height ≈ 0.005.

References
----------
- ЛСРМ Algorithmic Foundations 2022 § 8.4.1 «Табличное представление пика»
- Gilmore & Joss «Practical Gamma-ray Spectrometry» 3rd Ed. § 6.3
- Будыка А.К. «Спектрометрия ионизирующих излучений» 2021, § 7.4
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Sequence


# F-299 v1.18.3 — bridge from legacy peak_image (F-90) calibration
# to tabulated representation. Используется для конвертации
# существующих per-line T(E)/h_step(E) калибровок в .cpt.

def anchor_from_legacy_peak_image_params(
    E_keV: float,
    fwhm_keV: float,
    tail_param_T: float = 0.0,
    h_step_frac: float = 0.0,
    weight: float = 1.0,
) -> "PeakShapeAnchor":
    """Bridge: legacy peak_image (Gauss+tail+step) → tabulated anchor.

    Старая F-90 peak_image model хранит T-параметр в FWHM-единицах
    (точка перехода Gaussian → exponential tail в μ−T·σ). Для
    tabulated representation мы конвертируем:

      • tail_fraction ≈ exp(-T²/2)  — амплитуда tail на точке перехода
      • tail_slope_inv_keV ≈ T / σ  — экспоненциальный спад tail в 1/keV
      • step_height_frac = h_step_frac  — Compton step как доля от FEP

    Если T=0 — pure-Gaussian; tail_fraction=0.

    Parameters
    ----------
    E_keV : float
    fwhm_keV : float
    tail_param_T : float
        LSRM «хвостовой параметр» в долях FWHM (типично 0.5-0.9 для NaI).
    h_step_frac : float
        Compton step amplitude (доля FEP, типично 0.03 для NaI).
    weight : float

    Returns
    -------
    PeakShapeAnchor — готов к включению в TabulatedPeakImage.
    """
    if E_keV <= 0 or fwhm_keV <= 0:
        raise ValueError(
            f"E_keV and fwhm_keV must be > 0, got E={E_keV}, FWHM={fwhm_keV}"
        )
    sigma = fwhm_keV / 2.354820045
    if tail_param_T > 0.0 and sigma > 0.0:
        tail_fraction = math.exp(-0.5 * tail_param_T * tail_param_T)
        tail_slope_inv_keV = tail_param_T / sigma
    else:
        tail_fraction = 0.0
        tail_slope_inv_keV = 0.0
    return PeakShapeAnchor(
        E_keV=float(E_keV),
        fwhm_keV=float(fwhm_keV),
        tail_fraction=float(tail_fraction),
        tail_slope_inv_keV=float(tail_slope_inv_keV),
        step_height_frac=float(h_step_frac),
        asymmetry=0.0,
        weight=float(weight),
    )


@dataclass(frozen=True)
class PeakShapeAnchor:
    """Параметры формы пика в одной anchor-точке энергии."""

    E_keV: float
    fwhm_keV: float
    tail_fraction: float = 0.0      # 0 = pure Gaussian; для NaI ~0.03
    tail_slope_inv_keV: float = 0.0  # экспоненциальный спад tail (1/keV)
    step_height_frac: float = 0.0   # Compton step / FEP amplitude
    asymmetry: float = 0.0          # Voigt asymmetry (Lorentzian admixture)
    weight: float = 1.0             # вес anchor для spline-fit


@dataclass(frozen=True)
class TabulatedPeakImage:
    """Anchor-based peak shape definition (без spline-interpolation).

    Для interpolation используйте `peak_image_logspline.interpolate_peak_shape`
    с этим объектом как input.
    """

    detector_id: str                    # "Gamma-1S" / "HPGe-1" / etc.
    detector_class: str                 # "NaI" / "HPGe" / "CsI" / "LaBr" / "CeBr"
    crystal_diameter_mm: float
    anchors: List[PeakShapeAnchor] = field(default_factory=list)
    source_metadata: Optional[str] = None    # e.g. "Cs-137_2024-05-01"
    notes: Optional[str] = None

    def anchor_energies_sorted(self) -> List[float]:
        return sorted(a.E_keV for a in self.anchors)

    def anchor_at_E(self, E_keV: float, tolerance_keV: float = 1.0,
                    ) -> Optional[PeakShapeAnchor]:
        """Найти anchor с ближайшей E_keV в пределах tolerance."""
        if not self.anchors:
            return None
        best = min(self.anchors, key=lambda a: abs(a.E_keV - E_keV))
        if abs(best.E_keV - E_keV) <= tolerance_keV:
            return best
        return None

    def validate(self) -> List[str]:
        """Самопроверка: возвращает list of issue messages (empty = OK)."""
        issues: List[str] = []
        if not self.anchors:
            issues.append("no anchors defined")
            return issues
        # Все E > 0, FWHM > 0
        for a in self.anchors:
            if a.E_keV <= 0:
                issues.append(f"E_keV={a.E_keV} ≤ 0")
            if a.fwhm_keV <= 0:
                issues.append(f"FWHM={a.fwhm_keV} ≤ 0 at E={a.E_keV}")
            if not (0 <= a.tail_fraction <= 1):
                issues.append(
                    f"tail_fraction={a.tail_fraction} out of [0,1] "
                    f"at E={a.E_keV}"
                )
            if not (0 <= a.step_height_frac <= 1):
                issues.append(
                    f"step_height_frac={a.step_height_frac} out of [0,1] "
                    f"at E={a.E_keV}"
                )
        # FWHM монотонно растёт с E (физическая ожидание для всех детекторов)
        sorted_a = sorted(self.anchors, key=lambda a: a.E_keV)
        for i in range(len(sorted_a) - 1):
            if sorted_a[i + 1].fwhm_keV < sorted_a[i].fwhm_keV * 0.9:
                issues.append(
                    f"FWHM non-monotone: {sorted_a[i].fwhm_keV} @ "
                    f"{sorted_a[i].E_keV} → {sorted_a[i+1].fwhm_keV} @ "
                    f"{sorted_a[i+1].E_keV}"
                )
        # Detector class plausibility
        if self.detector_class == "NaI":
            for a in self.anchors:
                pct = 100.0 * a.fwhm_keV / a.E_keV
                if pct > 15 or pct < 3:
                    issues.append(
                        f"NaI FWHM% = {pct:.1f}% at E={a.E_keV} "
                        f"unusual (expected 3-15%)"
                    )
        return issues

    def estimate_fwhm_pct_at_662(self) -> Optional[float]:
        """FWHM @ 662 keV в процентах (canonical Cs-137 resolution metric)."""
        a = self.anchor_at_E(662.0, tolerance_keV=20.0)
        if a is None:
            return None
        return 100.0 * a.fwhm_keV / a.E_keV


def build_anchors_from_calibration(
    detector_id: str,
    detector_class: str,
    crystal_diameter_mm: float,
    calibration_pairs: Sequence[tuple],
) -> TabulatedPeakImage:
    """Convenience builder: из list of (E, FWHM) → TabulatedPeakImage.

    Для каждой пары создаётся anchor с default tail/step/asymmetry
    из preset detector_class (NaI: tail=0.03, step=0.05; HPGe: 0.001, 0.005).

    Parameters
    ----------
    calibration_pairs : sequence of (E_keV, fwhm_keV) | (E_keV, fwhm_keV, dict_overrides)
    """
    presets_by_class = {
        "NaI":  {"tail_fraction": 0.03, "tail_slope_inv_keV": 0.05,
                 "step_height_frac": 0.05, "asymmetry": 0.0},
        "CsI":  {"tail_fraction": 0.04, "tail_slope_inv_keV": 0.05,
                 "step_height_frac": 0.06, "asymmetry": 0.0},
        "HPGe": {"tail_fraction": 0.001, "tail_slope_inv_keV": 0.1,
                 "step_height_frac": 0.005, "asymmetry": 0.0},
        "LaBr": {"tail_fraction": 0.005, "tail_slope_inv_keV": 0.08,
                 "step_height_frac": 0.01, "asymmetry": 0.0},
        "CeBr": {"tail_fraction": 0.005, "tail_slope_inv_keV": 0.08,
                 "step_height_frac": 0.01, "asymmetry": 0.0},
    }
    defaults = presets_by_class.get(detector_class, presets_by_class["NaI"])

    anchors = []
    for pair in calibration_pairs:
        if len(pair) == 2:
            E, fwhm = pair
            overrides = {}
        elif len(pair) == 3:
            E, fwhm, overrides = pair
        else:
            raise ValueError(f"bad calibration pair: {pair}")
        merged = {**defaults, **(overrides or {})}
        anchors.append(PeakShapeAnchor(
            E_keV=float(E), fwhm_keV=float(fwhm), **merged,
        ))

    return TabulatedPeakImage(
        detector_id=detector_id,
        detector_class=detector_class,
        crystal_diameter_mm=crystal_diameter_mm,
        anchors=anchors,
    )


__all__ = [
    "PeakShapeAnchor",
    "TabulatedPeakImage",
    "build_anchors_from_calibration",
]
