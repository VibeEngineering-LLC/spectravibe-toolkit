"""
Minimum Detectable Activity (MDA), decision threshold, and detection
limit per ISO 11929:2019 — formalised in Lsrm Algorithmic Foundations
§6.3.

For each γ-ray line of a candidate nuclide, MDA quantifies the
**smallest activity that could be detected** given the measured
background under the peak ROI. MDA is the natural threshold for
identification: a nuclide is "present" if its activity (computed from
peak areas) exceeds the MDA of its characteristic line.

Three quantities are computed:

  1. **Decision threshold L_C** — the count rate above which we
     conclude "a peak is present" (rejecting H_0: no signal). At the
     decision threshold, the probability of a false-positive is α
     (typically 5%).
        L_C = k_α · σ_0
     where σ_0 is the uncertainty of the net count rate under the
     assumption that no signal is present.

  2. **Detection limit L_D** — the count rate at which we would
     definitely detect the signal with probability (1-β). At the
     detection limit, β is the probability of a false-negative
     (typically 5%, so 1-β = 95%).
        L_D solves: σ(L_D) · k_β + L_C ≤ L_D

  3. **Minimum Detectable Activity A_MDA** — L_D converted to activity:
        A_MDA = L_D / (ε · I · t)
     where ε is the detection efficiency, I is the γ-line intensity,
     and t is the live time.

For γ-ray spectrometry with background subtraction, the standard
formula is:

  L_C = k_α · √(2·n_bg)      [counts above bg]
  L_D = k_α² + 2·k_α·√(n_bg + L_C²/4)   [approximation]

  Common: k_α = k_β = 1.645 for α = β = 5%.

Per Lsrm §6.3 (Formula 6.3-7), for a measurement of duration t with
background rate n_0 in the ROI:

  L_C = √( n_0/t  +  n_0_f/t_0_f ) · k_α

where n_0_f is the background rate observed in a separate background
measurement of duration t_0_f. When the background spectrum is the
same as the sample (no separate bg), the second term vanishes.

Methodology references:
  - ISO 11929:2019 — Determination of characteristic limits
  - Lsrm Algorithmic Foundations 2022, §6.3
  - Currie L.A., Anal. Chem. 40 (1968) 586
  - Gilmore & Joss, "Practical Gamma-ray Spectrometry" 3rd Ed., §5.6
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional


# k-quantiles of the standard normal distribution for the most common
# significance levels. k_α with α=0.05 gives one-sided 95% confidence.
# Lsrm default is k=1.645 for both decision and detection.
K_ALPHA_95 = 1.645   # P(Z > k) = 5%
K_ALPHA_99 = 2.326   # P(Z > k) = 1%
K_ALPHA_999 = 3.090  # P(Z > k) = 0.1%


@dataclass(frozen=True)
class MdaResult:
    """Result of an MDA / decision-threshold / detection-limit calculation."""

    # Input parameters (for diagnostics)
    line_energy_keV: float
    background_counts_in_ROI: float
    live_time_s: float
    efficiency: float           # detection efficiency at E_keV
    intensity: float            # gamma intensity (decimal: 0.851 for Cs-137)
    k_alpha: float              # quantile used

    # Computed quantities
    decision_threshold_counts: float    # L_C in counts above bg
    detection_limit_counts: float       # L_D in counts above bg
    decision_threshold_cps: float       # L_C as count rate
    detection_limit_cps: float          # L_D as count rate
    MDA_Bq: float                       # Minimum Detectable Activity in Bq

    # Quality / diagnostics
    notes: str = ""

    def __repr__(self) -> str:
        return (f"MdaResult(E={self.line_energy_keV:.1f} keV, "
                f"L_C={self.decision_threshold_counts:.1f} cnt, "
                f"L_D={self.detection_limit_counts:.1f} cnt, "
                f"MDA={self.MDA_Bq:.2e} Bq)")


def mda_for_peak(
    *,
    line_energy_keV: float,
    background_counts_in_ROI: float,
    live_time_s: float,
    efficiency: float,
    intensity_pct: float,
    background_time_s: Optional[float] = None,
    background_ROI_counts_separate: Optional[float] = None,
    k_alpha: float = K_ALPHA_95,
    k_beta: Optional[float] = None,
    wings_baseline_n_roi_channels: Optional[int] = None,
    wings_baseline_m_each_side_channels: Optional[int] = None,
) -> MdaResult:
    """
    Compute decision threshold, detection limit, and MDA for one γ-line.

    Per ISO 11929:2019 and Lsrm §6.3 (Formulas 6.3-7, 6.3-9).

    Args:
        line_energy_keV: library energy of the line
        background_counts_in_ROI: net background counts under the peak
            ROI (after subtracting the polynomial baseline at the peak
            position, integrated over the ROI window).
        live_time_s: sample-spectrum live time
        efficiency: absolute photopeak detection efficiency at this
            energy (counts in the photopeak per γ emitted by source).
            Typical NaI 50×50 at 10 cm: ε(661 keV) ≈ 0.02; HPGe BEGe:
            ε(661 keV) ≈ 0.001 depending on geometry.
        intensity_pct: γ-line emission probability per decay (in %).
            Examples: Cs-137 = 85.1, K-40 = 10.66, Co-60(1173) = 99.85.
        background_time_s: live time of a separate background
            measurement (None if not used; in that case the second
            variance term vanishes).
        background_ROI_counts_separate: background counts in the same
            ROI from a separate background measurement (None if
            unavailable).
        k_alpha: quantile for decision threshold (default 1.645 for 5%)
        k_beta: quantile for detection limit (default = k_alpha)

    Returns:
        MdaResult with decision threshold, detection limit, and MDA.

    Notes on the formula (Lsrm §6.3):
        σ_0 = √( n_0/t + n_0_f/t_0_f )   [counts/s]
        L_C = k_α · σ_0  [counts/s, gross above background]
        L_D ≈ k_α² + 2·k_α·σ_0  [for k_α = k_β; small-correction term]
        Convert to counts: L_C_counts = L_C · t

        MDA = L_D / (ε · I · t)

        where I is the gamma intensity in decimal form (e.g. 0.851
        for Cs-137).

    Example (typical Cs-137 measurement on NaI 50×50):
        >>> r = mda_for_peak(
        ...     line_energy_keV=661.66,
        ...     background_counts_in_ROI=2500,
        ...     live_time_s=3600,
        ...     efficiency=0.02,
        ...     intensity_pct=85.1,
        ... )
        >>> r.decision_threshold_counts > 0
        True
        >>> r.MDA_Bq > 0
        True
    """
    if k_beta is None:
        k_beta = k_alpha

    if live_time_s <= 0:
        return MdaResult(
            line_energy_keV=line_energy_keV,
            background_counts_in_ROI=background_counts_in_ROI,
            live_time_s=live_time_s,
            efficiency=efficiency,
            intensity=intensity_pct / 100.0,
            k_alpha=k_alpha,
            decision_threshold_counts=0.0,
            detection_limit_counts=0.0,
            decision_threshold_cps=0.0,
            detection_limit_cps=0.0,
            MDA_Bq=float("inf"),
            notes="Invalid live_time_s ≤ 0",
        )

    if efficiency <= 0 or intensity_pct <= 0:
        return MdaResult(
            line_energy_keV=line_energy_keV,
            background_counts_in_ROI=background_counts_in_ROI,
            live_time_s=live_time_s,
            efficiency=efficiency,
            intensity=intensity_pct / 100.0,
            k_alpha=k_alpha,
            decision_threshold_counts=0.0,
            detection_limit_counts=0.0,
            decision_threshold_cps=0.0,
            decision_limit_cps=0.0 if False else 0.0,  # type: ignore
            detection_limit_cps=0.0,
            MDA_Bq=float("inf"),
            notes=f"Invalid efficiency ({efficiency}) or intensity ({intensity_pct}%)",
        )

    n_0 = max(0.0, float(background_counts_in_ROI))
    t = float(live_time_s)
    intensity = intensity_pct / 100.0

    # Variance of net count rate at H_0 (no signal). The sample
    # measurement contributes n_0/t. A separate background contributes
    # n_0_f/t_0_f.
    var_rate = n_0 / (t * t)  # σ² of N/t where N has variance N
    if (background_time_s is not None and background_time_s > 0
            and background_ROI_counts_separate is not None
            and background_ROI_counts_separate > 0):
        var_rate += float(background_ROI_counts_separate) / (
            float(background_time_s) ** 2
        )

    # F-273 (v1.17.11, T-006) — wings-baseline вариант для случая,
    # когда фон под пиком оценивается из m каналов с каждой стороны
    # ROI шириной n каналов (Currie 1968, ЛСРМ § 6.3-7 случай 2).
    # Дополнительная дисперсия фона из плеч:
    #     σ²_bg_in_ROI = (n/2m)² · n · b_per_ch = (n²/(2m)·b)
    # Полная дисперсия net = b·n + (n²/(2m))·b = b·n·(1 + n/(2m))
    # → factor (1 + n/(2m)) к Var(n_0).
    # Используется ТОЛЬКО когда переданы оба параметра. Каноничная
    # форма Currie σ²_net = B·(1 + n_roi/(2m_each_side)).
    if (wings_baseline_n_roi_channels is not None
            and wings_baseline_m_each_side_channels is not None
            and wings_baseline_n_roi_channels > 0
            and wings_baseline_m_each_side_channels > 0):
        n_roi = float(wings_baseline_n_roi_channels)
        m_side = float(wings_baseline_m_each_side_channels)
        wings_factor = 1.0 + n_roi / (2.0 * m_side)
        var_rate *= wings_factor

    sigma_0_cps = math.sqrt(var_rate)

    # Decision threshold (in count rate)
    L_C_cps = k_alpha * sigma_0_cps

    # Detection limit (in count rate). For k_α = k_β = k, the standard
    # approximate solution to L_D = k·σ_0 + k·σ(L_D) is
    #   L_D ≈ 2·L_C + k²/t
    # (Lsrm Formula 6.3-6 in the simplified case where dε and dI are
    # negligible compared to counting statistics; this is the usual
    # case for relatively short measurements.)
    L_D_cps = 2.0 * L_C_cps + (k_alpha * k_alpha) / t

    # Convert to counts
    L_C_counts = L_C_cps * t
    L_D_counts = L_D_cps * t

    # MDA in Bq:  MDA = L_D_cps / (ε · I)
    MDA_Bq = L_D_cps / (efficiency * intensity)

    return MdaResult(
        line_energy_keV=line_energy_keV,
        background_counts_in_ROI=background_counts_in_ROI,
        live_time_s=live_time_s,
        efficiency=efficiency,
        intensity=intensity,
        k_alpha=k_alpha,
        decision_threshold_counts=L_C_counts,
        detection_limit_counts=L_D_counts,
        decision_threshold_cps=L_C_cps,
        detection_limit_cps=L_D_cps,
        MDA_Bq=MDA_Bq,
    )


def characteristic_line_of_nuclide(
    nuclide_lines: list,
    mda_per_line: dict,
) -> Optional[tuple]:
    """
    Identify the **characteristic line** of a nuclide per Lsrm §6:
    the line with the LOWEST MDA. This is the most informative line
    for confirming/refuting the nuclide's presence.

    Args:
        nuclide_lines: list of (E_keV, I_pct, ...) tuples — the
            library lines of the nuclide
        mda_per_line: dict {E_keV: MdaResult} with computed MDA per
            line. Lines not in the dict are skipped.

    Returns:
        (E_keV, MdaResult) of the characteristic line, or None if the
        nuclide has no usable lines.
    """
    candidates = []
    for line in nuclide_lines:
        E = line[0] if isinstance(line, (list, tuple)) else line
        if E in mda_per_line:
            candidates.append((E, mda_per_line[E]))
    if not candidates:
        return None
    candidates.sort(key=lambda x: x[1].MDA_Bq)
    return candidates[0]


__all__ = [
    "K_ALPHA_95", "K_ALPHA_99", "K_ALPHA_999",
    "MdaResult", "mda_for_peak", "characteristic_line_of_nuclide",
]
