"""
F-282 (v1.17.13, T-004) — NaI / CsI scintillator-specific identification
defaults.

ORTEC GammaVision V9 и LSRM SpectraLine оптимизированы под HPGe (высокое
разрешение, узкие peaks). На NaI 63×63 Gamma-1S те же defaults дают:

  - **LibReduction=ON** → отсекает множество нуклидов, чьи library lines
    «слились» в общий ROI из-за широких NaI peaks. Результат: false
    NEGATIVE для типичных смешанных проб.
  - **LibPeakCL=ON** → жёсткое CL=99% gating отбраковывает слабые но
    реальные tracers.
  - **PeakOverlap=3.5·FWHM** (HPGe-default) → объявляет multiplet
    избыточно часто на NaI, где FWHM ~30-50 кэВ. Каждый ID ROI
    тонет в overlap-кластеры.

Каноничные NaI defaults (ORTEC AN66 Appendix A для scintillators):

  | Parameter       | HPGe default | NaI/CsI default |
  |---|---|---|
  | LibReduction    | ON           | **OFF**         |
  | LibPeakCL       | 0.99         | **OFF / 0.80**  |
  | PeakOverlap_FWHM| 3.5          | **2.0**         |
  | UseTCC          | ON           | OFF (см. T-007) |

Этот модуль экспонирует **только** константы и dispatch-функцию. Wiring
в `staged_pipeline.py` / `identify.py` — отдельная задача (deferred).
По умолчанию ничего не меняется (back-compat).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class IdentificationDefaults:
    """Detector-class-specific identification defaults."""
    detector_class: str
    lib_reduction: bool
    lib_peak_cl: float           # 0.0 = off (use any line); 0.99 = strict HPGe
    peak_overlap_fwhm: float     # FWHM multiplier for multiplet detection
    use_tcc: bool                # Total Cascade Correction — see T-007/T-082
    note: str = ""


# Каноничные defaults per ORTEC AN66 + LSRM Algorithmic Foundations §6.
_DEFAULTS_TABLE = {
    "HPGe": IdentificationDefaults(
        detector_class="HPGe",
        lib_reduction=True,
        lib_peak_cl=0.99,
        peak_overlap_fwhm=3.5,
        use_tcc=True,
        note="ORTEC HPGe canonical (AN66 Appendix A)",
    ),
    "CdZnTe": IdentificationDefaults(
        detector_class="CdZnTe",
        lib_reduction=True,
        lib_peak_cl=0.95,
        peak_overlap_fwhm=3.0,
        use_tcc=False,   # plane geometry → cascade correction small
        note="Semi-conductor; mid-resolution",
    ),
    "LaBr3": IdentificationDefaults(
        detector_class="LaBr3",
        lib_reduction=False,
        lib_peak_cl=0.90,
        peak_overlap_fwhm=2.5,
        use_tcc=False,
        note="LaBr3 intrinsic activity 138La requires tail-aware fit",
    ),
    "CeBr3": IdentificationDefaults(
        detector_class="CeBr3",
        lib_reduction=False,
        lib_peak_cl=0.90,
        peak_overlap_fwhm=2.5,
        use_tcc=False,
        note="CeBr3 cleaner intrinsic than LaBr3",
    ),
    "NaI": IdentificationDefaults(
        detector_class="NaI",
        lib_reduction=False,
        lib_peak_cl=0.0,
        peak_overlap_fwhm=2.0,
        use_tcc=False,
        note=("LSRM-aligned NaI defaults (F-282/T-004): LibReduction OFF, "
              "PeakOverlap = 2.0·FWHM (vs HPGe 3.5), TCC OFF unless "
              "close geometry (см. F-007/F-082)"),
    ),
    "CsI": IdentificationDefaults(
        detector_class="CsI",
        lib_reduction=False,
        lib_peak_cl=0.0,
        peak_overlap_fwhm=2.0,
        use_tcc=False,
        note="As NaI",
    ),
}


def get_identification_defaults(detector_class: str) -> IdentificationDefaults:
    """Вернуть рекомендованные identification-defaults для класса.

    Если класс не в таблице — возвращает NaI defaults (наиболее
    consurvative — НЕ reject слабые ID).
    """
    dc = str(detector_class).strip()
    if dc in _DEFAULTS_TABLE:
        return _DEFAULTS_TABLE[dc]
    # Soft prefix match
    dc_low = dc.lower()
    for key, val in _DEFAULTS_TABLE.items():
        if dc_low.startswith(key.lower()) or key.lower().startswith(dc_low):
            return val
    return _DEFAULTS_TABLE["NaI"]


__all__ = [
    "IdentificationDefaults",
    "get_identification_defaults",
]
