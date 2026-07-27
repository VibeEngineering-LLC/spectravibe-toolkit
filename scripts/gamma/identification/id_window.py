# -*- coding: utf-8 -*-
"""F-167 — Канонические FWHM-множители для всех контекстов гамма-спектрометрии.

Этот модуль — единый справочник `k·FWHM(E)` множителей. Любое использование
"magic-number" множителя в коде (`0.5`, `1.5`, `3.08`, ...) должно идти ЧЕРЕЗ
эти именованные константы — иначе теряется привязка к первоисточнику и
возникает риск семантической путаницы между разными контекстами.

Контексты различаются принципиально:
- **ID matching** — окно поиска совпадения найденного пика с библиотечной
  линией нуклида (Phase D identification).
- **PBC match** — окно сопоставления пика с peaked-background-вычитанием.
- **Multiplet criterion** — порог "является ли следующий пик частью того же
  мультиплета".
- **Multiplet padding** — поля региона мультиплета для аппроксимации фона.
- **Lib↔unknown association** — окно ассоциации между unknown-пиком и
  библиотечной энергией в deconvolution.
- **Background search** — ширина окна поиска фона методами Auto-5/3/1.
- **ISO NORM peak window** — окно для интегрирования по ISO 11929 Method 18/19.
- **Directed Fit** — диапазон функции Directed Fit (ORTEC GV9 §6.3.2.2).

Все ORTEC-множители — для cross-validation, ЛСРМ-множители — primary
(per F-157: ЛСРМ > Будыка > Gilmore для Gamma-1S NaI 63×63).

Источники (Layer 1 RAG-ID):
    [LSRM-Algo-9]                    — primary для ID-окна k=1.5 NaI / 1.0 HPGe
    [ORTEC-GV9-Background-Auto]      — 6·MW·FWHM background search
    [ORTEC-GV9-Background-FWHM]      — X.X·FWHM background width
    [ORTEC-GV9-Deconvolution-Width]  — 3.3·FWHM lib↔unknown association
    [ORTEC-GV9-Multiplet-Region]     — 1.5·FWHM padding + 3.08·FWHM grouping
    [ORTEC-GV9-PBC]                  — 0.5·FWHM PBC match width (eq.86-87)
    [ORTEC-GV9-ISO-NORM]             — 1.2 / 2.5 ·FWHM ISO Singlet
    [ORTEC-GV9-Engines]              — 4.84·FWHM Directed Fit base

Layer 2 (для пользовательских отчётов):
    [7, §9]                          — ЛСРМ-Algorithmic-Foundations
    [22, §6.3.1.3 / §6.3.7 / §6.5.1 / §6.10.4 / §6.3.2.2]  — ORTEC GammaVision V9

F-167 plan: audit/_plans/F-167_id_window_k_fwhm.md
"""
from __future__ import annotations

from typing import Literal


DetectorClass = Literal["NaI", "CsI", "LaBr", "CeBr", "HPGe", "CdZnTe"]


# ===========================================================================
# § 1.  ID matching — окно матчинга найденного пика с библиотечной энергией
# ===========================================================================
# Source: [LSRM-Algo-9] (ЛСРМ-Algorithmic-Foundations, §9, PRIMARY per F-157).
# Контекст: Phase D — на основе уже-откалиброванной FWHM(E) проверяем
#           |E_found - E_library| <= k(detector) · FWHM(E_library).
# Для NaI/CsI k=1.5; для HPGe/LaBr/CeBr/CdZnTe k=1.0 (узкое разрешение).
ID_WINDOW_K_FWHM: dict[DetectorClass, float] = {
    "NaI":    1.5,   # [LSRM-Algo-9] — primary
    "CsI":    1.5,   # extrapolated from NaI per [SHENDRIK-2]
    "LaBr":   1.0,   # higher resolution → narrower
    "CeBr":   1.0,
    "HPGe":   1.0,   # [LSRM-Algo-9]
    "CdZnTe": 1.0,
}


# ===========================================================================
# § 2.  Background search & subtraction
# ===========================================================================
# Source: [ORTEC-GV9-Background-Auto] §6.3.1.3.
# Метод Automatic 5/3/1-point: ширина окна поиска фона = 6 × MatchWidth × FWHM.
# MatchWidth — конфигурируемый параметр (default ORTEC = 1.0).
BACKGROUND_SEARCH_K_MW_FWHM: float = 6.0

# Source: [ORTEC-GV9-Background-FWHM] §6.3.1.3.
# Альтернативный метод "X.X·FWHM" — пользователь задаёт коэффициент;
# default ORTEC = 1.0.
BACKGROUND_X_FWHM: float = 1.0


# ===========================================================================
# § 3.  PBC (Peaked Background Correction)
# ===========================================================================
# Source: [ORTEC-GV9-PBC] §6.10.4.1, eq.86-87.
# Match width при сопоставлении PBC-таблицы с найденным пиком.
# Важно: 0.5 — это PBC-специфичный множитель, НЕ ID-окно.
PBC_MATCH_K_FWHM: float = 0.5


# ===========================================================================
# § 4.  Multiplet grouping (criterion & padding)
# ===========================================================================
# Source: [ORTEC-GV9-Multiplet-Region] §6.5.1.
# "Является ли следующий пик частью того же мультиплета": расстояние между
# центроидами <= MULTIPLET_CRITERION_K_FWHM × FWHM (max of two).
MULTIPLET_CRITERION_K_FWHM: float = 3.08

# Source: [ORTEC-GV9-Multiplet-Region] §6.5.1.
# Padding каждой стороны мультиплет-региона для аппроксимации фона.
MULTIPLET_PADDING_K_FWHM: float = 1.5


# ===========================================================================
# § 5.  Library ↔ unknown peak association (deconvolution)
# ===========================================================================
# Source: [ORTEC-GV9-Deconvolution-Width] §6.3.7.
# Окно ассоциации unknown-пика с библиотечной линией ВНУТРИ deconvolution
# (это НЕ ID matching; ID matching работает с уже-decomposed пиками).
LIB_UNKNOWN_ASSOC_K_FWHM: float = 3.3


# ===========================================================================
# § 6.  ISO 11929 NORM — Singlet method
# ===========================================================================
# Source: [ORTEC-GV9-ISO-NORM] Ann.C / §6.10.4 (Method 18/19 ISO 11929:2010).
# Peak-dominant: ширина окна охватывает практически весь пик (~99.4%).
ISO_NORM_PEAK_DOMINANT_K_FWHM: float = 2.5
# Background-dominant: 1.2·FWHM ≈ 84% площади (требует f-correction).
ISO_NORM_BG_DOMINANT_K_FWHM: float = 1.2


# ===========================================================================
# § 7.  Directed Fit (ORTEC engine option)
# ===========================================================================
# Source: [ORTEC-GV9-Engines] §6.3.2.2.
# Базовая ширина диапазона функции Directed Fit; полная ширина =
# DIRECTED_FIT_K_FWHM_BASE × MatchWidth × DirectedFitFactor.
DIRECTED_FIT_K_FWHM_BASE: float = 4.84


# ===========================================================================
# § 8.  Phase D regularization (внутренний skill-параметр, НЕ ORTEC)
# ===========================================================================
# Source: F-145 / F-117 contracts.
# Мягкий зажим центроидов в Phase D self-calibration: после совмещения пика с
# библиотечной линией центроид может "плавать" в пределах
# ±PHASE_D_REGULARIZATION_K_FWHM × FWHM, чтобы дать пространство микро-дрейфу
# калибровки. Малый множитель (0.15) — НЕ матчинг-окно, а bounds регуляризации.
# Legacy-имя `PHASE_D_CENTROID_TOLERANCE_FRAC` сохраняется как алиас
# в `gamma.calibration.multiplet_self_calibration` (см. F-167 подзадачу 5).
PHASE_D_REGULARIZATION_K_FWHM: float = 0.15


# ===========================================================================
# Helpers
# ===========================================================================

def normalize_detector_class(s: str | None) -> DetectorClass:
    """Привести произвольную строку detector_type к каноничному DetectorClass.

    Поддерживает legacy-входы вида "NaI 63×63", "HPGe Canberra GR2018",
    "LaBr3", "CeBr3", "CdZnTe", регистро-независимо. По умолчанию (None
    или unknown) возвращает "NaI" (профиль Gamma-1S по F-157).

    Args:
        s: исходная строка с типом детектора.

    Returns:
        Один из: "NaI" | "CsI" | "LaBr" | "CeBr" | "HPGe" | "CdZnTe".
    """
    if not s:
        return "NaI"
    low = s.strip().lower()
    if "hpge" in low or low.startswith("ge"):
        return "HPGe"
    if "cdznte" in low or "czt" in low:
        return "CdZnTe"
    if "labr" in low:
        return "LaBr"
    if "cebr" in low:
        return "CeBr"
    if low.startswith("csi") or "csi" in low.replace("csi(", "csi"):
        return "CsI"
    # Default — NaI (Gamma-1S profile)
    return "NaI"


def id_window_keV(
    energy_keV: float,
    fwhm_keV: float,
    detector: DetectorClass = "NaI",
) -> float:
    """ID-окно ±k·FWHM(E) для матчинга найденного пика с библиотекой.

    Каноничный критерий по ЛСРМ-Algorithmic-Foundations §9:

        |E_found - E_library| < k(detector) · FWHM(E_library)

    где k = 1.5 для NaI/CsI и k = 1.0 для HPGe/LaBr/CeBr/CdZnTe.

    Args:
        energy_keV: энергия библиотечного пика (только для сигнатуры,
            формула не зависит от E — она через FWHM(E)).
        fwhm_keV: калиброванное FWHM на этой энергии (из F-168).
        detector: класс детектора, см. ``DetectorClass``. По умолчанию "NaI".

    Returns:
        Полная ширина окна (full width, не half-width) в кэВ.
        Для ±-сравнения используй ``window/2.0``.

    Источники:
        [LSRM-Algo-9] §9 (PRIMARY) — ЛСРМ-Algorithmic-Foundations.
        Cross-check: [ORTEC-GV9-Deconvolution-Width] — 3.3·FWHM ассоциация
                     lib↔unknown ВНУТРИ deconvolution (НЕ ID окно).

    Example:
        >>> round(id_window_keV(1461.0, 120.0, "NaI"), 1)
        360.0
        >>> round(id_window_keV(1461.0, 1.8, "HPGe"), 2)
        3.6
    """
    k = ID_WINDOW_K_FWHM[detector]
    return 2.0 * k * fwhm_keV


__all__ = [
    "DetectorClass",
    "ID_WINDOW_K_FWHM",
    "BACKGROUND_SEARCH_K_MW_FWHM",
    "BACKGROUND_X_FWHM",
    "PBC_MATCH_K_FWHM",
    "MULTIPLET_CRITERION_K_FWHM",
    "MULTIPLET_PADDING_K_FWHM",
    "LIB_UNKNOWN_ASSOC_K_FWHM",
    "ISO_NORM_PEAK_DOMINANT_K_FWHM",
    "ISO_NORM_BG_DOMINANT_K_FWHM",
    "DIRECTED_FIT_K_FWHM_BASE",
    "PHASE_D_REGULARIZATION_K_FWHM",
    "normalize_detector_class",
    "id_window_keV",
]
