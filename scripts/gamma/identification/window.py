"""
Identification window δE(E) — the energy interval around a candidate
γ-line within which a found peak is considered to match.

Per Lsrm Algorithmic Foundations §6, identification searches for a peak
within an energy window centered on each library line. The window size
scales with energy following the detector's resolution curve:

  • Scintillators (NaI, CeBr, LaBr, CsI): δE(E) = δE₀ · √(E/661.66)
    — square-root scaling reflects the underlying FWHM(E) ∝ √E law
    of scintillator counting statistics. The reference δE₀ is the
    window at 661.66 keV (Cs-137 photopeak), a configuration parameter.

  • Semiconductor (HPGe): δE(E) ≈ const + small slope. Standard NaI
    formula over-shrinks at high E where HPGe FWHM grows slowly. For
    HPGe a constant ~1 keV window or linear-in-E model is used.

Typical δE₀ values (Lsrm-recommended):
  - HPGe         : 1.0 keV  (very narrow — characteristic of HPGe)
  - LaBr3/CeBr3  : 6.0 keV  (better than NaI)
  - NaI 50×50    : 15.0 keV (typical scintillator)
  - CsI/Tracor   : 20.0 keV (lower-quality scintillators)

The window must be **wide enough** to encompass calibration drift and
peak-fit uncertainty, yet **narrow enough** to avoid matching peaks
that belong to unrelated nearby lines. δE₀ in the 0.5–2× FWHM_at_661
range is typical.

Reference: Lsrm Algorithmic Foundations 2022, §6, formula
δE(E) = δE₀ · √(E/E_ref) with E_ref = 661.66 keV.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Callable, Union


# Reference energy at which the calibration parameter δE₀ is defined.
# This is the Cs-137 137mBa γ — the most common reference line.
REFERENCE_ENERGY_KEV = 661.66


# Default identification window at the reference energy, in keV.
# Sensible per-detector defaults from the Lsrm methodology.
DEFAULT_DELTA_E0_KEV = {
    "HPGe": 1.0,
    "LaBr3": 6.0,
    "CeBr3": 6.0,
    "NaI": 15.0,
    "CsI": 20.0,
    "unknown": 15.0,  # NaI default — most common
}


@dataclass(frozen=True)
class IdentificationWindow:
    """
    Bundle of an identification window function with diagnostic info.

    Supported scalings:
      - "sqrt_E"   : ЛСРМ §6 canonical δE(E) = δE₀ · √(E/E_ref)
      - "linear"   : HPGe-style δE(E) = δE₀ + slope · (E − E_ref)
      - "constant" : δE(E) = δE₀
      - "k_fwhm"   : F-167 canonical δE(E) = k · FWHM(E) per [LSRM-Algo-9]
                     — requires `k_fwhm` and `fwhm_provider_keV`.
      - "lsrm_sqrt_with_kfwhm_floor" : F-279 (v1.17.12, T-001) primary
                     ЛСРМ √E формула + floor k·FWHM(E). Берётся MAX:
                         max(δE₀·√(E/E_ref), k·FWHM(E))
                     Защищает high-E хвост от слишком узкого окна на
                     нелинейных FWHM-профилях. Требует и k_fwhm и
                     fwhm_provider_keV.
    """
    detector_type: str
    delta_E0_keV: float  # Window width at REFERENCE_ENERGY_KEV
    scaling: str  # "sqrt_E" | "linear" | "constant" | "k_fwhm" | "lsrm_sqrt_with_kfwhm_floor"
    # F-167 — populated only when scaling involves k_fwhm.
    k_fwhm: Optional[float] = None
    fwhm_provider_keV: Optional[Callable[[float], float]] = None

    def window_keV(self, E_keV: float) -> float:
        """Return half-window width in keV at the given energy."""
        if self.scaling == "k_fwhm":
            # F-167 canonical per [LSRM-Algo-9] §9: half-width = k · FWHM(E).
            # ID matching condition: |E_found - E_lib| <= k · FWHM(E_lib).
            if self.fwhm_provider_keV is None or self.k_fwhm is None:
                raise ValueError(
                    "F-167: scaling='k_fwhm' requires both `k_fwhm` and "
                    "`fwhm_provider_keV` to be set."
                )
            fwhm = float(self.fwhm_provider_keV(E_keV))
            if fwhm <= 0 or not math.isfinite(fwhm):
                # Защитный fallback: если FWHM-провайдер вернул бессмыслицу
                # (ноль/NaN на краях диапазона), используем δE₀ как страховку.
                return float(self.delta_E0_keV)
            return float(self.k_fwhm) * fwhm
        if self.scaling == "lsrm_sqrt_with_kfwhm_floor":
            # F-279 (v1.17.12, T-001) — primary ЛСРМ §6 √E формула с
            # k·FWHM(E) floor. Берётся MAX, чтобы:
            #   • На low-E (Am-241 60 кэВ) — sqrt_E даёт ~4.5 кэВ; для
            #     NaI FWHM ~7 кэВ × k=1.5 = 10.5 → floor спасает.
            #   • На high-E (Tl-208 2614) — sqrt_E даёт ~30 кэВ; FWHM
            #     ~120 × k=1.5 = 180 → kfwhm доминирует, как и должно.
            #   • На середине (Cs-137 661) — обе формулы совпадают.
            if self.fwhm_provider_keV is None or self.k_fwhm is None:
                raise ValueError(
                    "F-279: scaling='lsrm_sqrt_with_kfwhm_floor' requires "
                    "both `k_fwhm` and `fwhm_provider_keV`."
                )
            sqrt_e = self.delta_E0_keV * math.sqrt(E_keV / REFERENCE_ENERGY_KEV)
            try:
                fwhm = float(self.fwhm_provider_keV(E_keV))
                if fwhm <= 0 or not math.isfinite(fwhm):
                    return sqrt_e
                floor = float(self.k_fwhm) * fwhm
                return max(sqrt_e, floor)
            except Exception:
                return sqrt_e
        if self.scaling == "sqrt_E":
            return self.delta_E0_keV * math.sqrt(E_keV / REFERENCE_ENERGY_KEV)
        if self.scaling == "linear":
            # Linear approximation suitable for HPGe whose FWHM grows
            # slowly. δE(E) = δE₀ + (slope_per_keV)·(E − E_ref).
            # Use small slope so window doubles roughly every 5 MeV.
            slope_per_keV = self.delta_E0_keV / 5000.0
            return self.delta_E0_keV + slope_per_keV * (E_keV - REFERENCE_ENERGY_KEV)
        # constant
        return self.delta_E0_keV

    def in_window(self, E_query_keV: float, E_library_keV: float) -> bool:
        """True if a found peak at E_query is within the window around
        the library line E_library."""
        return abs(E_query_keV - E_library_keV) <= self.window_keV(E_library_keV)


def build_identification_window(
    detector_type: str,
    delta_E0_keV: Optional[float] = None,
) -> IdentificationWindow:
    """
    Build an IdentificationWindow appropriate for a detector type.

    Args:
        detector_type: one of "HPGe", "NaI", "LaBr3", "CeBr3", "CsI",
            "unknown". Case-insensitive prefix match.
        delta_E0_keV: explicit window at REFERENCE_ENERGY_KEV. If None,
            use detector-default from DEFAULT_DELTA_E0_KEV.

    Returns:
        IdentificationWindow with appropriate scaling.

    Examples:
        >>> w = build_identification_window("NaI")
        >>> round(w.window_keV(661.66), 2)
        15.0
        >>> round(w.window_keV(2614.51), 2)
        29.82
        >>> w.in_window(661.0, 661.66)
        True
        >>> w.in_window(680.0, 661.66)
        False
    """
    # Normalise detector type
    dt = (detector_type or "unknown").strip()
    # Case-insensitive prefix match
    dt_lower = dt.lower()
    canonical = "unknown"
    for key in DEFAULT_DELTA_E0_KEV:
        if dt_lower.startswith(key.lower()) or key.lower().startswith(dt_lower):
            canonical = key
            break
    # If not matched as prefix, also try substring
    if canonical == "unknown" and dt_lower != "unknown":
        for key in DEFAULT_DELTA_E0_KEV:
            if key.lower() in dt_lower or dt_lower in key.lower():
                canonical = key
                break

    if delta_E0_keV is None:
        delta_E0_keV = DEFAULT_DELTA_E0_KEV[canonical]

    # Scaling: HPGe = linear (resolution barely depends on E),
    # everything else (scintillators) = sqrt_E.
    if canonical == "HPGe":
        scaling = "linear"
    else:
        scaling = "sqrt_E"

    return IdentificationWindow(
        detector_type=canonical,
        delta_E0_keV=float(delta_E0_keV),
        scaling=scaling,
    )


def identification_window_from_fwhm(
    fwhm_at_661_keV: float,
    detector_type: str = "NaI",
    fwhm_multiple: float = 0.5,
) -> IdentificationWindow:
    """
    Derive an identification window from a measured FWHM at 661 keV.

    Lsrm-style: δE₀ ≈ 0.5·FWHM at 661.66 keV. This gives a window that
    matches the detector's actual resolution rather than a fixed
    nuclide-library tolerance.

    Args:
        fwhm_at_661_keV: measured FWHM at 661.66 keV in keV
        detector_type: detector type for scaling-law selection
        fwhm_multiple: factor multiplying FWHM to get δE₀ (default 0.5)

    Returns:
        IdentificationWindow.

    Example:
        >>> # Typical NaI 50×50: FWHM @ 661 = 50 keV, so δE₀ = 25 keV
        >>> w = identification_window_from_fwhm(50.0, "NaI")
        >>> round(w.delta_E0_keV, 1)
        25.0
        >>> round(w.window_keV(2614.51), 1)
        49.7
    """
    delta_E0 = fwhm_at_661_keV * fwhm_multiple
    return build_identification_window(detector_type, delta_E0_keV=delta_E0)


def build_id_window_k_fwhm(
    detector_class: str,
    fwhm_provider_keV: Callable[[float], float],
    *,
    k_override: Optional[float] = None,
) -> IdentificationWindow:
    """F-167 — Build canonical ID window per [LSRM-Algo-9] §9: ±k·FWHM(E).

    Returns a IdentificationWindow with scaling="k_fwhm" that, instead of
    the legacy `δE₀·√(E/E_ref)` heuristic, evaluates the half-window as
    `k(detector) · FWHM(E)` at each query energy. This is the canonical
    matching criterion from the LSRM Algorithmic Foundations methodology:

        |E_found − E_library| <= k(detector) · FWHM(E_library)

    with k = 1.5 for NaI/CsI and k = 1.0 for HPGe/LaBr/CeBr/CdZnTe.

    Args:
        detector_class: one of "NaI", "CsI", "LaBr", "CeBr", "HPGe", "CdZnTe".
            Free-form input is normalized via
            `gamma.identification.id_window.normalize_detector_class`.
        fwhm_provider_keV: callable(E_keV) → FWHM(E) in keV (from F-168
            scintillation FWHM model or similar).
        k_override: optional explicit k value overriding the per-detector
            default (use only for custom/research detectors).

    Returns:
        IdentificationWindow with scaling="k_fwhm".

    Examples:
        >>> w = build_id_window_k_fwhm("NaI", lambda E: 0.05 * math.sqrt(E))
        >>> round(w.window_keV(661.66), 2)
        1.93
        >>> w.scaling
        'k_fwhm'
        >>> w.k_fwhm
        1.5

    Source: [LSRM-Algo-9] §9; cross-check [ORTEC-GV9-Deconvolution-Width].
    """
    from gamma.identification.id_window import (
        ID_WINDOW_K_FWHM, normalize_detector_class,
    )
    canonical = normalize_detector_class(detector_class)
    k = float(k_override) if k_override is not None else float(ID_WINDOW_K_FWHM[canonical])
    # Pre-evaluate FWHM at reference energy for diagnostic δE₀ field;
    # never used in `window_keV` when scaling == "k_fwhm".
    try:
        fwhm_ref = float(fwhm_provider_keV(REFERENCE_ENERGY_KEV))
        if not math.isfinite(fwhm_ref) or fwhm_ref <= 0:
            fwhm_ref = float(DEFAULT_DELTA_E0_KEV.get(canonical, 15.0))
    except Exception:
        fwhm_ref = float(DEFAULT_DELTA_E0_KEV.get(canonical, 15.0))
    return IdentificationWindow(
        detector_type=canonical,
        delta_E0_keV=k * fwhm_ref,  # diagnostic-only — для __repr__/логов
        scaling="k_fwhm",
        k_fwhm=k,
        fwhm_provider_keV=fwhm_provider_keV,
    )


def build_id_window_lsrm_with_kfwhm_floor(
    detector_class: str,
    fwhm_provider_keV: Callable[[float], float],
    *,
    delta_E0_keV: Optional[float] = None,
    k_override: Optional[float] = None,
) -> IdentificationWindow:
    """F-279 (v1.17.12, T-001) — ЛСРМ √E primary с k·FWHM(E) floor.

    Returns IdentificationWindow с scaling="lsrm_sqrt_with_kfwhm_floor".
    Окно = max(δE₀·√(E/E_ref), k·FWHM(E)). Объединяет каноничную
    ЛСРМ §6 формулу и canonical F-167 k·FWHM(E) как floor.

    Default δE₀ — DEFAULT_DELTA_E0_KEV[detector_class] (NaI=15 кэВ).
    Default k — из ID_WINDOW_K_FWHM (NaI=1.5).
    """
    from gamma.identification.id_window import (
        ID_WINDOW_K_FWHM, normalize_detector_class,
    )
    canonical = normalize_detector_class(detector_class)
    de0 = float(delta_E0_keV) if delta_E0_keV is not None \
        else float(DEFAULT_DELTA_E0_KEV.get(canonical, 15.0))
    k = float(k_override) if k_override is not None \
        else float(ID_WINDOW_K_FWHM[canonical])
    return IdentificationWindow(
        detector_type=canonical,
        delta_E0_keV=de0,
        scaling="lsrm_sqrt_with_kfwhm_floor",
        k_fwhm=k,
        fwhm_provider_keV=fwhm_provider_keV,
    )


__all__ = [
    "REFERENCE_ENERGY_KEV",
    "DEFAULT_DELTA_E0_KEV",
    "IdentificationWindow",
    "build_identification_window",
    "identification_window_from_fwhm",
    "build_id_window_k_fwhm",
    "build_id_window_lsrm_with_kfwhm_floor",
]
