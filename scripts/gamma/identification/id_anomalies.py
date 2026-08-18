"""
F-276 (v1.17.12, T-028 + T-030 + T-025) — Identification anomaly rules.

Detects three common ID pitfalls and emits structured annotations so
downstream code can flag/suppress without rewriting identification
logic itself:

  1. **511 keV ID veto** (T-028)
     The 511 keV line appears in **every** spectrum acquired in a
     shielded geometry due to cosmic-ray pair production interacting
     with the lead shield ("annihilation peak"). Promoting Na-22 /
     Cu-64 / Zn-65 SOLELY on the basis of a 511 keV match produces
     false-positives. Rule: if the only matched line is 511 keV ± 2σ,
     veto the ID.

  2. **Pb 72 keV fluorescence flag** (T-030)
     Graded shields (Pb-lined with Cd / Cu inner layers) emit
     characteristic Pb K-α / K-β X-rays around 72-75 keV when
     irradiated by γ-rays > 88 keV (Pb K-edge). On NaI 63×63 these
     fluorescence X-rays form a clear peak in every sample with even
     modest background γ activity. Misidentified as Eu-155 (86) or
     other low-energy γ-emitters.

  3. **Iodine K-escape annotation** (T-025)
     For NaI(Tl), incident γ-rays in 100-1000 keV range can ionize
     iodine K-shell → photoelectron escapes with K_α emission
     (~28.6 keV from I). This produces a small **escape peak** at
     E_γ − 28.6 keV. For Am-241 (59.5 keV photopeak → escape at
     ~30.9 keV), the escape peak can be larger than the photopeak.
     Annotation tells downstream code this is a satellite, not an
     independent nuclide.

All rules are **annotations only** — they don't auto-remove anything.
Caller decides whether to veto / warn / proceed based on the flags.

References
----------
- Knoll "Radiation Detection and Measurement" 4th Ed. §10 (NaI escape)
- Gilmore & Joss §11 (X-ray fluorescence in shielded chambers)
- ORTEC AN66 §8 (511 keV annihilation handling)
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional


# ──────────────────────────────────────────────────────────────────
# T-028 — 511 keV ID veto
# ──────────────────────────────────────────────────────────────────

# Nuclides that COULD legitimately contribute to 511 keV (β+ emitters)
# but should NEVER be identified solely on this line.
BETA_PLUS_NUCLIDES_REQUIRING_CORROBORATION = frozenset({
    "Na-22",   # 1274 keV pair line — REQUIRE this for ID
    "Cu-64",   # also 1346 keV
    "Zn-65",   # also 1115 keV
    "F-18",    # only 511 — never standalone IDable on γ spectrometry
    "Ge-68",   # daughter Ga-68 511 only
    "C-11",
    "N-13",
    "O-15",
})

# Tolerance (keV) around 511 to count as the annihilation line.
ANNIHILATION_LINE_TOLERANCE_KEV = 3.0


def is_only_511_line_match(
    matched_line_energies_keV: Iterable[float],
) -> bool:
    """True если единственная matched линия — это 511 кэВ ± tolerance.

    Использовать в Step-3 promotion: если nuclide ∈ BETA_PLUS_NUCLIDES
    и is_only_511_line_match(...) → veto (никакого confirmed).
    """
    matched = [float(e) for e in matched_line_energies_keV]
    if not matched:
        return False
    if len(matched) > 1:
        return False
    return abs(matched[0] - 511.0) <= ANNIHILATION_LINE_TOLERANCE_KEV


def annihilation_veto_reason(nuclide: str) -> Optional[str]:
    """Вернуть текст veto-причины если нуклид ∈ β+ списка и matched
    только 511, иначе None.

    Caller pattern:
        if (n := annihilation_veto_reason(nuclide)) and is_only_511_line_match(...):
            mark_as_rejected(nuc, reason=n)
    """
    if nuclide in BETA_PLUS_NUCLIDES_REQUIRING_CORROBORATION:
        return (
            f"F-276/T-028: {nuclide} требует corroborating γ-линии помимо "
            f"511 кэВ — annihilation peak присутствует во ВСЕХ спектрах из-за "
            f"cosmic pair-production в свинцовом домике. Без второй линии "
            f"идентификация considered false-positive."
        )
    return None


# ──────────────────────────────────────────────────────────────────
# T-030 — Pb fluorescence flag (72-75 keV)
# ──────────────────────────────────────────────────────────────────

# Pb K-α₁ = 74.97 keV, K-α₂ = 72.80 keV, K-β = 84.94/87.36 keV.
PB_FLUORESCENCE_LINES_KEV = (72.80, 74.97, 84.94, 87.36)
PB_FLUORESCENCE_WINDOW_KEV = 3.0   # NaI resolution at 75 keV ≈ ±8 keV;
                                   # narrow 3 keV для центра K-α доминанты


def is_pb_fluorescence_line(
    found_E_keV: float,
    window_keV: float = PB_FLUORESCENCE_WINDOW_KEV,
) -> bool:
    """True если найденная линия попадает в окно одной из Pb fluorescence
    K-лекций. Используется для предупреждения «возможно, X-ray из shield».
    """
    for pb in PB_FLUORESCENCE_LINES_KEV:
        if abs(found_E_keV - pb) <= window_keV:
            return True
    return False


def pb_fluorescence_note(found_E_keV: float) -> Optional[str]:
    """Вернуть текст-аннотацию если линия выглядит как Pb fluorescence."""
    if not is_pb_fluorescence_line(found_E_keV):
        return None
    return (
        f"F-276/T-030: пик ~{found_E_keV:.1f} кэВ совпадает с Pb K-X-ray "
        f"флуоресценции графитового/свинцового домика "
        f"(Pb Kα₁=74.97, Kα₂=72.80, Kβ=84.94/87.36 кэВ). Учитывать как "
        f"артефакт shield, а не нуклид."
    )


# ──────────────────────────────────────────────────────────────────
# T-025 — Iodine K-escape annotation (NaI only)
# ──────────────────────────────────────────────────────────────────

IODINE_K_ALPHA_ENERGY_KEV = 28.612    # I K-α₁; K-shell escape энергия
IODINE_K_ESCAPE_RELEVANT_MAX_KEV = 200.0  # выше escape pulse слаб
IODINE_K_ESCAPE_MIN_PARENT_KEV = 35.0     # ниже не имеет смысла escape


@dataclass(frozen=True)
class IodineEscapeAnnotation:
    parent_E_keV: float          # photopeak energy E_γ
    escape_E_keV: float          # E_γ − I_K_α (~28.6)
    note: str


def iodine_k_escape_for(
    parent_E_keV: float, detector_class: str = "NaI",
) -> Optional[IodineEscapeAnnotation]:
    """Вернуть annotation о ожидаемом K-escape peak для NaI/CsI ниже 200 кэВ.

    Для не-NaI/CsI или E_γ вне диапазона — возвращает None.

    Пример (Am-241 59.5 кэВ):
        >>> a = iodine_k_escape_for(59.54)
        >>> round(a.escape_E_keV, 2)
        30.93
    """
    if detector_class.strip() not in ("NaI", "CsI"):
        return None
    if not (IODINE_K_ESCAPE_MIN_PARENT_KEV <= parent_E_keV
            <= IODINE_K_ESCAPE_RELEVANT_MAX_KEV):
        return None
    esc_E = parent_E_keV - IODINE_K_ALPHA_ENERGY_KEV
    if esc_E <= 0:
        return None
    return IodineEscapeAnnotation(
        parent_E_keV=float(parent_E_keV),
        escape_E_keV=float(esc_E),
        note=(
            f"F-276/T-025: для photopeak {parent_E_keV:.1f} кэВ на "
            f"{detector_class} ожидается I K-escape satellite на "
            f"{esc_E:.1f} кэВ (= E_γ − {IODINE_K_ALPHA_ENERGY_KEV:.1f} кэВ). "
            f"Если в спектре найден пик ~{esc_E:.1f} кэВ — он, скорее всего, "
            f"K-escape, не отдельный нуклид."
        ),
    )


def find_iodine_escape_candidates(
    found_peak_energies_keV: Iterable[float],
    parent_photopeak_keV: float,
    detector_class: str = "NaI",
    tolerance_keV: float = 5.0,
) -> List[float]:
    """Найти среди found-peaks те, что попадают в K-escape окно
    относительно parent_photopeak."""
    ann = iodine_k_escape_for(parent_photopeak_keV, detector_class)
    if ann is None:
        return []
    matches = []
    for E in found_peak_energies_keV:
        if abs(float(E) - ann.escape_E_keV) <= tolerance_keV:
            matches.append(float(E))
    return matches


__all__ = [
    # T-028
    "BETA_PLUS_NUCLIDES_REQUIRING_CORROBORATION",
    "ANNIHILATION_LINE_TOLERANCE_KEV",
    "is_only_511_line_match",
    "annihilation_veto_reason",
    # T-030
    "PB_FLUORESCENCE_LINES_KEV",
    "PB_FLUORESCENCE_WINDOW_KEV",
    "is_pb_fluorescence_line",
    "pb_fluorescence_note",
    # T-025
    "IODINE_K_ALPHA_ENERGY_KEV",
    "IodineEscapeAnnotation",
    "iodine_k_escape_for",
    "find_iodine_escape_candidates",
]
