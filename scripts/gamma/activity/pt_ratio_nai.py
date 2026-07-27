"""
F-295 (v1.17.20, T-011) — Peak-to-total ratio P/T(E) for NaI detectors.

Peak-to-total ratio (P/T) — доля полного-энергетического пика (FEP)
в общем числе зарегистрированных событий для моноэнергетического
источника. Используется в:

  • **TCS correction** (T-007) — расчёт total efficiency ε_T = ε_FEP / P/T
  • **Matrix method** (T-027) — построение sensitivity matrix
  • **Compton-continuum bookkeeping** — нормировка фоновых событий

Источник данных
---------------
Gilmore & Joss «Practical Gamma-ray Spectrometry» 3rd Ed., **Table 8.4**
(p. 225) — приблизительные значения P/T для NaI(Tl) при
source-to-detector ≈ 25 cm. Различия с реальными значениями для
Marinelli-геометрии (близкая) могут достигать ±5 %, что приемлемо
для бутстрап-расчётов.

Поддерживаемые размеры детектора (с интерполяцией):

  • **3"×3"** (76×76 mm)
  • **4"×4"** (102×102 mm)
  • **Gamma-1S 63×63** — extrapolated (меньше 3"×3" на 17 %)

Для произвольного размера используется линейная интерполяция по
диаметру (small-effect для P/T в первом приближении).

References
----------
- Gilmore & Joss «Practical Gamma-ray Spectrometry» 3rd Ed., Table 8.4
- Heath R.L. ANCR-1000-2 (1974) — extended P/T tables for NaI
- ЛСРМ Algorithmic Foundations 2022 § 10 «Каскадное суммирование»
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Optional, Sequence, Tuple


# Gilmore Table 8.4 (Practical Gamma-ray Spectrometry 3rd ed., p. 225).
# (E_keV, P/T) для двух размеров детектора NaI(Tl). Стандартная geometry:
# source 25 cm от детектора (uncollimated point source).
GILMORE_TABLE_8_4_3IN3 = [
    (100.0, 0.92),
    (200.0, 0.84),
    (500.0, 0.55),
    (1000.0, 0.36),
    (1500.0, 0.27),
    (2000.0, 0.22),
    (2500.0, 0.18),
]

GILMORE_TABLE_8_4_4IN4 = [
    (100.0, 0.95),
    (200.0, 0.89),
    (500.0, 0.68),
    (1000.0, 0.50),
    (1500.0, 0.40),
    (2000.0, 0.34),
    (2500.0, 0.30),
]

# Размеры (диаметр кристалла в мм) для интерполяции.
_DETECTOR_DIAMETERS_MM = {
    "3in3": 76.0,    # 3"×3"
    "4in4": 102.0,   # 4"×4"
    "Gamma-1S": 63.0,
    "NaI_63x63": 63.0,
}


def _interp_loglog(
    table: Sequence[Tuple[float, float]], E_keV: float,
) -> float:
    """Log-log линейная интерполяция P/T(E)."""
    if E_keV <= table[0][0]:
        return table[0][1]
    if E_keV >= table[-1][0]:
        return table[-1][1]
    for i in range(len(table) - 1):
        E0, p0 = table[i]
        E1, p1 = table[i + 1]
        if E0 <= E_keV <= E1:
            ln_p = (
                math.log(p0)
                + (math.log(p1) - math.log(p0))
                * (math.log(E_keV) - math.log(E0))
                / (math.log(E1) - math.log(E0))
            )
            return math.exp(ln_p)
    return table[-1][1]


def pt_ratio_nai(
    E_keV: float,
    crystal_diameter_mm: float = 63.0,
) -> float:
    """Peak-to-total ratio P/T(E) для NaI(Tl) данного размера.

    Parameters
    ----------
    E_keV : float
        Энергия γ-кванта.
    crystal_diameter_mm : float
        Диаметр кристалла. Default 63 mm (Gamma-1S). Поддерживаются
        промежуточные значения через линейную интерполяцию между
        3"×3" (76 mm) и 4"×4" (102 mm).

    Returns
    -------
    P/T ratio (decimal, 0..1).

    Notes
    -----
    Для diameters < 76 mm (Gamma-1S 63 mm) используется linear extrapolation
    с поправкой ~ -2 % относительно 3"×3" — на основании эмпирического
    наблюдения снижения geometric efficiency для меньших кристаллов.
    """
    if E_keV <= 0:
        raise ValueError("E_keV must be > 0")

    pt_3in3 = _interp_loglog(GILMORE_TABLE_8_4_3IN3, E_keV)
    pt_4in4 = _interp_loglog(GILMORE_TABLE_8_4_4IN4, E_keV)

    d3, d4 = 76.0, 102.0

    if crystal_diameter_mm <= d3:
        # Extrapolation для меньших кристаллов: -2 % per cm меньше 76 mm
        delta_cm = (d3 - crystal_diameter_mm) / 10.0
        adjustment = -0.02 * delta_cm
        return max(0.05, min(1.0, pt_3in3 + adjustment))

    if crystal_diameter_mm >= d4:
        # Extrapolation для больших: cap at 4"×4" value (P/T saturates)
        return pt_4in4

    # Linear interpolation между 3"×3" и 4"×4"
    t = (crystal_diameter_mm - d3) / (d4 - d3)
    return pt_3in3 * (1.0 - t) + pt_4in4 * t


def pt_ratio_for_detector(E_keV: float, detector_id: str) -> float:
    """Convenience wrapper: P/T по preset-имени детектора."""
    if detector_id not in _DETECTOR_DIAMETERS_MM:
        raise KeyError(
            f"Unknown detector_id '{detector_id}'. "
            f"Known: {sorted(_DETECTOR_DIAMETERS_MM.keys())}"
        )
    return pt_ratio_nai(E_keV, _DETECTOR_DIAMETERS_MM[detector_id])


def total_efficiency_from_fep(
    eps_fep: float, E_keV: float, crystal_diameter_mm: float = 63.0,
) -> float:
    """Total efficiency ε_T = ε_FEP / P/T(E).

    Используется в TCS correction формуле:
        C_TCS = 1 / (1 - Σ_i p_i · ε_T(E_i))
    где p_i — coincidence probability.
    """
    pt = pt_ratio_nai(E_keV, crystal_diameter_mm)
    if pt <= 0:
        raise ValueError(f"P/T = {pt} (≤ 0) at E={E_keV}")
    return eps_fep / pt


__all__ = [
    "GILMORE_TABLE_8_4_3IN3",
    "GILMORE_TABLE_8_4_4IN4",
    "pt_ratio_nai",
    "pt_ratio_for_detector",
    "total_efficiency_from_fep",
]
