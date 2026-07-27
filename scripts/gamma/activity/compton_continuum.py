"""
F-304 (v1.18.0, T-022c) — Analytic Compton continuum model для NaI.

Per-line Compton continuum представлен в виде:

  - Plateau на [E_Compton_min .. E_Compton_edge] с константной amplitude
  - Smoothed step-down за Compton-edge (gaussian smoothing с FWHM детектора)
  - Linear ramp down к E ≈ 0 (multiple-Compton + ackscatter tail)

Используется как B(ch) pedestal-составляющая в F-302 PPP-templates.
Total area под continuum = compton_area = peak_area·(1-P/T)/(P/T),
где P/T берётся из F-295 (Gilmore Table 8.4).

Klein-Nishina kinematics:
    E_Compton_edge(E) = E · 2αE / (1 + 2αE),  где α = E/m_e c² (m_e c² = 511 keV)

Это аналитическое приближение — для precision-mode желателен Monte-Carlo
response (deferred).

References
----------
- Gilmore §2.2.2 Compton scattering kinematics
- Knoll "Radiation Detection and Measurement" 4th ed §10.III
- ЛСРМ §13 quasi-template + LSRM §8 peak-image foundations
"""
from __future__ import annotations

import math
from typing import Callable, Optional

ELECTRON_REST_MASS_KEV = 510.9989461   # m_e c² [keV]


def compton_edge_keV(E_keV: float) -> float:
    """Compton-edge energy для incident photon E [keV].

    E_c = E · (2αE/(1+2αE)), α = 1/m_e c²
    Эквивалентная форма: E_c = 2E²/(m_e c² + 2E)
    """
    if E_keV <= 0.0:
        return 0.0
    return 2.0 * E_keV * E_keV / (ELECTRON_REST_MASS_KEV + 2.0 * E_keV)


def backscatter_peak_keV(E_keV: float) -> float:
    """Backscatter peak energy = E - E_compton_edge."""
    return max(0.0, E_keV - compton_edge_keV(E_keV))


# ──────────────────────────────────────────────────────────────────
# Continuum shape (per-line pedestal)
# ──────────────────────────────────────────────────────────────────

def compton_continuum_for_line(
    E_line_keV: float,
    compton_area: float,
    n_channels: int,
    channel_to_keV: Callable[[int], float],
    fwhm_keV_at_edge: Optional[float] = None,
    backscatter_fraction: float = 0.10,
) -> list[float]:
    """Build per-channel continuum counts для one gamma-line.

    Shape:
      - Plateau [E_backscatter .. E_edge] с roughly константной amplitude
        (Klein-Nishina dN/dE монотонно растёт к edge, но для NaI с σ ~ 50 keV
        smoothed → ≈ flat)
      - Smoothed step за edge (erf-smoothed by detector FWHM)
      - Backscatter peak ~10% area в районе E_backscatter
      - Линейное затухание от backscatter к E=0

    Parameters
    ----------
    E_line_keV : float
        Energy of the parent FEP line.
    compton_area : float
        Total counts to distribute в continuum (compton_area =
        peak_area·(1-P/T)/(P/T) from F-295).
    n_channels : int
    channel_to_keV : Callable[[int], float]
    fwhm_keV_at_edge : Optional[float]
        FWHM of detector at the Compton-edge energy. Если None — берём
        E_edge/14 (≈ 7% NaI default).
    backscatter_fraction : float
        Fraction of compton_area направляется в backscatter feature.
        Default 0.10 (Gilmore §2 typical NaI Marinelli).

    Returns
    -------
    list[float] длины n_channels.
    """
    if compton_area <= 0.0 or E_line_keV <= 0.0:
        return [0.0] * n_channels
    if n_channels <= 0:
        return []

    E_edge = compton_edge_keV(E_line_keV)
    E_back = backscatter_peak_keV(E_line_keV)
    if fwhm_keV_at_edge is None or fwhm_keV_at_edge <= 0.0:
        fwhm_edge = max(E_edge / 14.0, 1.0)
    else:
        fwhm_edge = fwhm_keV_at_edge
    sigma_edge = fwhm_edge / 2.354820045

    plateau_area = compton_area * (1.0 - backscatter_fraction)
    backscatter_area = compton_area * backscatter_fraction

    plateau_low = E_back
    plateau_high = E_edge
    plateau_width = max(plateau_high - plateau_low, 1.0)
    plateau_height_per_keV = plateau_area / plateau_width

    # Backscatter peak FWHM ≈ FWHM at E_back
    fwhm_back = max(fwhm_edge * 0.7, 0.5)
    sigma_back = fwhm_back / 2.354820045

    SQRT2 = math.sqrt(2.0)

    counts = [0.0] * n_channels
    for ch in range(n_channels):
        e_lo = channel_to_keV(ch - 0.5)
        e_hi = channel_to_keV(ch + 0.5)
        if e_hi < e_lo:
            e_lo, e_hi = e_hi, e_lo
        if e_hi <= 0.0:
            continue

        # 1. Plateau с erf-smoothed step за edge
        # smoothed plateau: integral{e_lo..e_hi} plateau_height·Φ((E_edge-E)/sigma_edge) dE
        # для constant plateau с smoothed step:
        if e_lo < plateau_high + 5.0 * sigma_edge:
            # Часть из plateau (sharp на low side, smoothed на high side)
            # Эффективная высота: plateau_height_per_keV · weight_in_plateau_region
            # Усредняем upper-half-CDF от Gaussian centered at E_edge
            mid = 0.5 * (e_lo + e_hi)
            if mid > plateau_low and mid < plateau_high + 5.0 * sigma_edge:
                bin_w = e_hi - e_lo
                # smoothing factor (1 - Φ((E - E_edge)/sigma)) ≈ 0.5·erfc((E-E_edge)/(sigma·sqrt2))
                z = (mid - E_edge) / (sigma_edge * SQRT2)
                smoothing = 0.5 * math.erfc(z)
                counts[ch] += plateau_height_per_keV * bin_w * smoothing

        # 2. Backscatter peak
        if E_back > 0.0 and backscatter_area > 0.0:
            z_hi = (e_hi - E_back) / (sigma_back * SQRT2)
            z_lo = (e_lo - E_back) / (sigma_back * SQRT2)
            cdf_diff = 0.5 * (math.erf(z_hi) - math.erf(z_lo))
            counts[ch] += backscatter_area * cdf_diff

    # Normalize: integral counts должен == compton_area ± rounding
    total_actual = sum(counts)
    if total_actual > 0.0:
        rescale = compton_area / total_actual
        if abs(rescale - 1.0) > 0.01:
            counts = [c * rescale for c in counts]

    return counts


def make_continuum_func(
    fwhm_keV_at: Optional[Callable[[float], float]] = None,
    backscatter_fraction: float = 0.10,
) -> Callable:
    """Factory для использования compton_continuum_for_line как `continuum_func`
    в F-302 `build_nuclide_template`.
    """
    def _func(E_line_keV, compton_area, n_channels, channel_to_keV):
        fwhm_edge = None
        if fwhm_keV_at is not None:
            try:
                fwhm_edge = fwhm_keV_at(compton_edge_keV(E_line_keV))
            except Exception:
                fwhm_edge = None
        return compton_continuum_for_line(
            E_line_keV=E_line_keV,
            compton_area=compton_area,
            n_channels=n_channels,
            channel_to_keV=channel_to_keV,
            fwhm_keV_at_edge=fwhm_edge,
            backscatter_fraction=backscatter_fraction,
        )
    return _func


__all__ = [
    "ELECTRON_REST_MASS_KEV",
    "compton_edge_keV",
    "backscatter_peak_keV",
    "compton_continuum_for_line",
    "make_continuum_func",
]
