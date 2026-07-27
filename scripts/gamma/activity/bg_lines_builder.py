"""
F-298 (v1.17.20, T-013) — Background lines builder (F-96 ↔ F-131 bridge).

F-96 (ЛСРМ §9.4 «A-priori peak list of background lines») задаёт
канонический список ожидаемых background-линий: K-40, Th-232 chain
(Tl-208, Ac-228), U-238 chain (Bi-214, Pb-214), Bi-207, Co-60 etc.

F-131 (deconvolution coupled-multiplet structure) принимает список
ожидаемых peak.E_keV для constrained-fit.

Этот модуль строит **bridge**: преобразует F-96 a-priori library в
F-131-совместимый input для deconvolution фоновых ROI.

Use-cases
---------
1. **bg spectrum analysis** — выявить вклад K-40 1461 keV и Th-chain
   2614 keV без полной peak search.
2. **NORM-spectrum classification** — отделить natural radionuclide
   contribution (K + Th + U) от anthropogenic.
3. **Calibration anchor extraction** — bg линии K-40 / Tl-208 как
   anchors для bootstrap energy calibration в «slow-mode» (long bg
   acquisition).

References
----------
- ЛСРМ Algorithmic Foundations 2022 § 9.4 «Peak-list a-priori
  background line representation»
- Gilmore & Joss «Practical Gamma-ray Spectrometry» 3rd Ed. § 5
  «Background and shielding»
- ICRU Report 53 (1994) — environmental NORM library
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence


@dataclass(frozen=True)
class BgLine:
    """Одна background-линия из F-96 a-priori library."""
    nuclide: str
    E_keV: float
    intensity_decimal: float    # I per decay (decimal)
    parent_chain: str           # "Th-232" / "U-238" / "U-235" / "K-40" / "Bi-207" / ...
    is_anchor_candidate: bool = False     # подходит для E-calibration anchor


# F-96 canonical a-priori background library для NaI low-resolution.
# Источник: ЛСРМ §9.4, IAEA Reference materials, NORM dataset.
F96_BG_LIBRARY: List[BgLine] = [
    # K-40 — single line, always present
    BgLine("K-40", 1460.83, 0.1066, "K-40", is_anchor_candidate=True),

    # Th-232 chain (in radioactive equilibrium, common in soil/concrete)
    BgLine("Ac-228", 338.32, 0.1127, "Th-232"),
    BgLine("Ac-228", 911.20, 0.2585, "Th-232"),
    BgLine("Ac-228", 968.97, 0.1550, "Th-232"),
    BgLine("Pb-212", 238.63, 0.4330, "Th-232", is_anchor_candidate=True),
    BgLine("Tl-208", 583.19, 0.3060, "Th-232", is_anchor_candidate=True),
    BgLine("Tl-208", 860.56, 0.0450, "Th-232"),
    BgLine("Tl-208", 2614.51, 0.3585, "Th-232", is_anchor_candidate=True),
    BgLine("Bi-212", 727.33, 0.0667, "Th-232"),

    # U-238 chain (Ra-226 sub-chain)
    BgLine("Pb-214", 295.22, 0.1842, "U-238"),
    BgLine("Pb-214", 351.93, 0.3560, "U-238"),
    BgLine("Bi-214", 609.31, 0.4549, "U-238", is_anchor_candidate=True),
    BgLine("Bi-214", 768.36, 0.0489, "U-238"),
    BgLine("Bi-214", 1120.29, 0.1492, "U-238"),
    BgLine("Bi-214", 1764.49, 0.1531, "U-238", is_anchor_candidate=True),
    BgLine("Bi-214", 2204.21, 0.0489, "U-238"),
    BgLine("Ra-226", 186.21, 0.0359, "U-238"),

    # U-235 chain (presence ≈ 4.25 % natural U)
    BgLine("U-235", 143.76, 0.1096, "U-235"),
    BgLine("U-235", 163.33, 0.0508, "U-235"),
    BgLine("U-235", 185.72, 0.5740, "U-235"),    # F-005 NORM apportionment
    BgLine("U-235", 205.31, 0.0501, "U-235"),

    # Bi-207 (cosmic / Pb shield activation)
    BgLine("Bi-207", 569.69, 0.9774, "Bi-207"),
    BgLine("Bi-207", 1063.66, 0.7460, "Bi-207"),
    BgLine("Bi-207", 1770.23, 0.0686, "Bi-207"),

    # 511 keV annihilation (cosmic / pair production)
    BgLine("Annihilation", 511.00, 1.0, "annihilation"),

    # Co-60 contamination (lab background after spill)
    BgLine("Co-60", 1173.23, 0.9985, "Co-60_contamination"),
    BgLine("Co-60", 1332.49, 0.9998, "Co-60_contamination"),
]


@dataclass(frozen=True)
class F131DeconvolutionInput:
    """Input для F-131 coupled-multiplet fit на одной ROI."""
    roi_E_min_keV: float
    roi_E_max_keV: float
    expected_peaks_keV: List[float]
    expected_peak_labels: List[str] = field(default_factory=list)
    intensity_weights: List[float] = field(default_factory=list)


def filter_bg_lines_in_window(
    E_min_keV: float, E_max_keV: float,
    min_intensity: float = 0.005,
    parent_chains: Optional[Sequence[str]] = None,
    library: Sequence[BgLine] = F96_BG_LIBRARY,
) -> List[BgLine]:
    """Найти все F-96 bg линии в окне [E_min; E_max].

    Parameters
    ----------
    min_intensity : float
        Игнорировать линии слабее этого (decimal, default 0.005 = 0.5 %).
    parent_chains : sequence of str | None
        Фильтр по chain: только эти chains. None → все.
    """
    out = []
    for line in library:
        if not (E_min_keV <= line.E_keV <= E_max_keV):
            continue
        if line.intensity_decimal < min_intensity:
            continue
        if parent_chains is not None and line.parent_chain not in parent_chains:
            continue
        out.append(line)
    return out


def build_f131_input(
    roi_E_min_keV: float,
    roi_E_max_keV: float,
    min_intensity: float = 0.005,
    parent_chains: Optional[Sequence[str]] = None,
) -> F131DeconvolutionInput:
    """Построить F-131 deconvolution input из F-96 bg library.

    Returns
    -------
    F131DeconvolutionInput с expected_peaks_keV в порядке возрастания E.
    """
    lines = filter_bg_lines_in_window(
        roi_E_min_keV, roi_E_max_keV, min_intensity, parent_chains,
    )
    lines_sorted = sorted(lines, key=lambda L: L.E_keV)
    return F131DeconvolutionInput(
        roi_E_min_keV=roi_E_min_keV,
        roi_E_max_keV=roi_E_max_keV,
        expected_peaks_keV=[L.E_keV for L in lines_sorted],
        expected_peak_labels=[
            f"{L.nuclide} {L.E_keV:.2f}" for L in lines_sorted
        ],
        intensity_weights=[L.intensity_decimal for L in lines_sorted],
    )


def get_anchor_candidates(
    library: Sequence[BgLine] = F96_BG_LIBRARY,
) -> List[BgLine]:
    """Линии, помеченные как is_anchor_candidate (для bootstrap E-cal)."""
    return [L for L in library if L.is_anchor_candidate]


def classify_chain_dominance(
    detected_lines_keV: Sequence[float],
    tolerance_keV: float = 5.0,
) -> Dict[str, int]:
    """Подсчитать сколько detected lines принадлежит каждой parent_chain.

    Используется для quick «NORM vs anthropogenic» классификации.

    Returns
    -------
    dict[chain → count].
    """
    counts: Dict[str, int] = {}
    for E in detected_lines_keV:
        for line in F96_BG_LIBRARY:
            if abs(line.E_keV - E) <= tolerance_keV:
                counts[line.parent_chain] = counts.get(line.parent_chain, 0) + 1
                break
    return counts


__all__ = [
    "BgLine",
    "F96_BG_LIBRARY",
    "F131DeconvolutionInput",
    "filter_bg_lines_in_window",
    "build_f131_input",
    "get_anchor_candidates",
    "classify_chain_dominance",
]
