"""
Confidence Index (CI) for nuclide identification per Lsrm
Algorithmic Foundations §14.3.

CI quantifies **how unambiguous** a successful identification is. The
intuition: if a nuclide has a single line and the only thing checked
is its energy, identification is weak — many other nuclides could fit
in the same energy window. If a nuclide has many lines with known
relative intensities and all of them appear in the spectrum at the
right positions with the right ratios, identification is strong.

The formula (Lsrm §14.3):

   CI = log( 1 / (δE_1 · δE_2 · ... · δI_2 · δI_3 · ...) )

where:
   δE_i = ΔE_i / E_i — relative uncertainty in energy of line i
                       (essentially the fraction of the identification
                       window E occupies in the spectrum)
   δI_j = ΔI_j / I_j — relative uncertainty in intensity ratio
                       (for lines 2 and onwards; the first line is
                       used to fix the scale, so no δI_1 contribution)

Reading the formula intuitively:
   - Each δE_i is a small fraction (~10⁻²–10⁻³); product over many
     lines → very small.
   - Taking log reverses small → big, so MORE lines → BIGGER CI.
   - Adding intensity-ratio uncertainty further narrows the
     parameter space → bigger CI.

Calibration values for NaI (Lsrm Table 14-1):

   Cs-137: CI = 1.8   (single 661 keV line — low confidence)
   K-40:   CI = 2.2   (one line at 1460 keV)
   Na-22:  CI = 3.8   (511 + 1274.5 keV)
   Co-60:  CI = 5.9   (1173 + 1332 keV, known ratio)
   Cs-134: CI = 4.4
   Ba-133: CI = 8.5   (5+ lines)
   Eu-152: CI = 18.3  (very many lines — extremely confident)
   Th-232: CI = 16.6  (full chain with many lines)

For HPGe (much narrower windows) CI values are much higher. Lsrm
Table 14-1 also reports HPGe and LaBr3 values.

Use in identification:
   - CI < 5: low confidence — single line / few lines, easy to confuse
             (Cs-137 1.8, K-40 2.2, Na-22 3.8, Cs-134 4.4 on NaI — all low)
   - 5 ≤ CI < 10: moderate — multi-line, plausible (Co-60 5.9, Ba-133 8.5)
   - CI ≥ 10: high — many lines + ratios — robust ID (Eu-152 18.3, Th-232 16.6)

Reference: Lsrm Algorithmic Foundations 2022 §14.3.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class ConfidenceIndexResult:
    """Detailed CI computation for one candidate nuclide."""

    nuclide: str
    n_lines_used: int
    energy_uncertainty_factor: float    # Π δE_i (small)
    intensity_uncertainty_factor: float # Π δI_j for j≥2 (small or 1.0)
    CI: float                           # log(1/Π) — bigger = more confident

    # Per-line breakdown (energy, δE_i)
    per_line_dE: tuple = ()
    # Per-line breakdown (energy, δI_j) for j≥2
    per_line_dI: tuple = ()

    notes: str = ""

    def confidence_level(self) -> str:
        """Discretise CI into low/moderate/high per Lsrm intuition.

        Threshold 5.0 between low/moderate aligns with LSRM Table 14-1 NaI
        calibration: Cs-137 (1.8), K-40 (2.2), Na-22 (3.8), Cs-134 (4.4)
        are all single-line / few-line identifications that should read as
        low-confidence even when they pass the ci_gating CI≥2.0 confirmed
        threshold via anchor-rank promotion.
        """
        if self.CI < 5.0:
            return "low"
        if self.CI < 10.0:
            return "moderate"
        return "high"

    def __repr__(self) -> str:
        return (f"CI({self.nuclide}: {self.CI:.2f}, "
                f"{self.confidence_level()}, "
                f"{self.n_lines_used} lines)")


def confidence_index(
    nuclide: str,
    matched_lines: list,
    window_at: callable,
) -> ConfidenceIndexResult:
    """
    Compute Confidence Index for one candidate nuclide identification.

    Args:
        nuclide: name of the nuclide (for the result label)
        matched_lines: list of dicts with keys
                       {"E_keV": float, "I_pct": float,
                        "dI_pct": Optional[float]}
                       — the LIBRARY lines that were successfully
                       matched in the spectrum.
        window_at: callable(E_keV) → window half-width in keV (typically
                   from IdentificationWindow.window_keV)

    Returns:
        ConfidenceIndexResult.

    Examples (using a NaI-typical window of 15 keV at 661 + sqrt scaling):
        >>> from gamma.identification.window import build_identification_window
        >>> w = build_identification_window("NaI", delta_E0_keV=15.0)
        >>> # Cs-137: single line at 661.66 keV
        >>> result = confidence_index("Cs-137",
        ...                           [{"E_keV": 661.66, "I_pct": 85.1}],
        ...                           w.window_keV)
        >>> 1.5 < result.CI < 2.5
        True

        >>> # Co-60: two lines with known intensities
        >>> result = confidence_index("Co-60",
        ...     [{"E_keV": 1173.23, "I_pct": 99.85, "dI_pct": 0.03},
        ...      {"E_keV": 1332.49, "I_pct": 99.98, "dI_pct": 0.02}],
        ...     w.window_keV)
        >>> result.CI > 5
        True
    """
    if not matched_lines:
        return ConfidenceIndexResult(
            nuclide=nuclide, n_lines_used=0,
            energy_uncertainty_factor=1.0,
            intensity_uncertainty_factor=1.0,
            CI=0.0,
            notes="No matched lines",
        )

    # Energy uncertainty product: Π δE_i = Π (window(E_i) / E_i)
    energy_factors = []
    per_line_dE = []
    for line in matched_lines:
        E = float(line["E_keV"])
        if E <= 0:
            continue
        dE = float(window_at(E))
        delta_E_rel = dE / E
        if delta_E_rel <= 0 or delta_E_rel >= 1:
            continue
        energy_factors.append(delta_E_rel)
        per_line_dE.append((E, delta_E_rel))

    if not energy_factors:
        return ConfidenceIndexResult(
            nuclide=nuclide, n_lines_used=0,
            energy_uncertainty_factor=1.0,
            intensity_uncertainty_factor=1.0,
            CI=0.0,
            notes="No valid energies",
        )

    energy_prod = 1.0
    for f in energy_factors:
        energy_prod *= f

    # Intensity-ratio uncertainty product: Π δI_j for j≥2.
    # When relative intensities are known precisely, this factor
    # narrows the identification further. If dI_pct is not provided,
    # we use a default ~ 1% uncertainty (conservative).
    intensity_factors = []
    per_line_dI = []
    if len(matched_lines) > 1:
        # Reference line = first matched line (highest intensity by
        # convention is best, but order doesn't change the product)
        for line in matched_lines[1:]:
            I = float(line.get("I_pct", 0.0))
            if I <= 0:
                continue
            dI = float(line.get("dI_pct", 0.0))
            if dI <= 0:
                # Use a 1% default relative if unknown — Lsrm uses
                # this magnitude implicitly in their table
                delta_I_rel = 0.01
            else:
                delta_I_rel = dI / I
            if delta_I_rel <= 0 or delta_I_rel >= 1:
                continue
            intensity_factors.append(delta_I_rel)
            per_line_dI.append((float(line["E_keV"]), delta_I_rel))

    intensity_prod = 1.0
    for f in intensity_factors:
        intensity_prod *= f

    # CI = log10( 1 / (energy_prod · intensity_prod) )
    full_prod = energy_prod * intensity_prod
    CI = -math.log10(full_prod) if full_prod > 0 else 0.0

    return ConfidenceIndexResult(
        nuclide=nuclide,
        n_lines_used=len(energy_factors),
        energy_uncertainty_factor=energy_prod,
        intensity_uncertainty_factor=intensity_prod,
        CI=CI,
        per_line_dE=tuple(per_line_dE),
        per_line_dI=tuple(per_line_dI),
    )


__all__ = ["ConfidenceIndexResult", "confidence_index"]
