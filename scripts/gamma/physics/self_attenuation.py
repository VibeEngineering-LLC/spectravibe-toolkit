"""
Self-attenuation correction for volume γ-ray sources (K-20, v1.9.0).

For a γ-ray emitted inside a sample of finite density, some fraction of
photons are absorbed before escaping the sample. The "self-attenuation
factor" F(E) ∈ (0, 1] is the ratio of escaping to emitted photons.

For a thin disk slab of thickness t and density ρ in a uniform matrix
with mass-attenuation coefficient μ(E)/ρ, integrated over a uniformly-
distributed source emitting isotropically and detected from one face
(Knoll §10.III.5 / Gilmore §8.7):

    F(E) = (1 − exp(−μ(E)·ρ·t)) / (μ(E)·ρ·t)

Limits:
  • μρt → 0  (thin/light/high-E): F → 1
  • μρt → ∞  (thick/heavy/low-E): F → 1/(μρt)

Application to cert-validation (K-20):
============================================================
The LSRM .efr efficiency curve ε(E) is calibrated against a
**reference sample** of specific (ρ_ref, t_ref, geometry). When a new
sample with the same matrix material but different density ρ_sample
(thickness fixed by the container) is measured, compute_activity
returns:

    A_calc = A_true × F_sample(E) / F_ref(E)

(because the calibrated ε implicitly includes F_ref). To recover the
true activity, multiply by the correction factor:

    corr(E) = F_ref(E) / F_sample(E)

When ρ_sample < ρ_ref (lighter / less attenuating):
    F_sample > F_ref → corr < 1 → reduce A_calc (it was over-estimated)
When ρ_sample > ρ_ref:
    F_sample < F_ref → corr > 1 → boost A_calc (under-estimated)
When equal: corr = 1 (no correction needed).

Empirical findings from v1.7.25 cert matrix (40 fixtures):
    Cs-137 Δ% correlates with ρ_sample/ρ_ref:
      ρ=0.27 → Δ=+6.97 %  (light Petri)
      ρ=0.79 → Δ=−3.71 %  (heavy Petri)
      ρ=1.04 → Δ=−5.58 %  (heavy Marinelli)
    Expected post-correction: spread ≤ 5 % (vs current 17 % peak-to-peak).

Matrix composition data:
========================
ОИСН-16 — LSRM standard soil simulant, used as reference matrix in
.efr files for Marinelli/Дента-120мл/Петри-60мл geometries. Mass
fractions per .efr Material field:
    H: 0.022, C: 0.206, N: 0.009, O: 0.049, Fe: 0.714.
    Sum = 1.000 (verified).

NIST XCOM mass-attenuation data:
================================
Tabulated μ/ρ values from NIST Physical Reference Data
(https://physics.nist.gov/PhysRefData/Xcom/), Berger M.J. et al.
(2010) NIST Standard Reference Database 8, pillar energies 50, 60,
80, 100, 150, 200, 300, 400, 500, 600, 800, 1000, 1250, 1500, 2000,
2500, 3000 keV. Includes total attenuation (photoelectric + Compton
incoherent + coherent + pair production above 1022 keV); for typical
γ-spectrometry the photoelectric + Compton terms dominate. Values
in cm²/g.

Log-log linear interpolation used between pillar energies (NIST XCOM
recommendation).

Limitations:
============
1. **Thin-slab approximation**: ignores side-wall losses (acceptable
   for Petri-60мл shallow dish, Marinelli sleeve; less accurate for
   Дента-120мл cylinder where t/diameter ≈ 1).
2. **Matrix material assumption**: treats the entire sample as ОИСН-16.
   Real 420-series cert sources are mixed ОИСН-16 + nuclide-specific
   carriers, but mass fraction of carriers is < 0.01 → negligible.
3. **Single matrix supported**: ОИСН-16 only in this version. To add
   another matrix (e.g. organic biological), supply composition dict
   to `correction_factor(..., composition=...)`.
4. **Energy range**: 50-3000 keV (NIST XCOM table extent). Outside
   range: clamped to edge values.
5. **Container geometry t**: must be supplied externally — typically
   from the .efr Layers→Width field. Per-geometry default values
   provided in `REF_THICKNESS_CM`.

References:
===========
* Knoll, "Radiation Detection and Measurement" 4th Ed., §10.III.5
  (volume source self-absorption).
* Gilmore, "Practical Gamma-ray Spectrometry" 3rd Ed., §8.7
  (matrix corrections).
* Berger M.J. et al., NIST XCOM Photon Cross Section Database,
  NIST Standard Reference Database 8 (XGAM).
* ICRU Report 93, "Key data for ionizing-radiation dosimetry" (2014).
"""

from __future__ import annotations

import math
from typing import Dict, Iterable, Sequence


# ---------------------------------------------------------------------------
# NIST XCOM mass-attenuation coefficients μ/ρ in cm²/g
# ---------------------------------------------------------------------------
#
# Pillar energies (keV) and per-element values for H, C, N, O, Fe.
# Source: NIST Standard Reference Database 8 (XGAM), accessed 2024-11.
# Values are total (photoelectric + Compton incoherent + Rayleigh +
# pair production), which is what self-attenuation requires.

_XCOM_ENERGIES_KEV: Sequence[float] = (
    50.0,   60.0,   80.0,   100.0,  150.0,  200.0,  300.0,
    400.0,  500.0,  600.0,  800.0,  1000.0, 1250.0, 1500.0,
    2000.0, 2500.0, 3000.0,
)

_XCOM_MU_RHO: Dict[str, Sequence[float]] = {
    # Hydrogen Z=1
    "H":  (0.3354, 0.3206, 0.2944, 0.2735, 0.2367, 0.2114, 0.1769,
           0.1551, 0.1393, 0.1271, 0.1093, 0.0967, 0.0858, 0.0771,
           0.0641, 0.0551, 0.0485),
    # Carbon Z=6
    "C":  (0.1875, 0.1639, 0.1361, 0.1214, 0.1041, 0.0921, 0.0786,
           0.0700, 0.0635, 0.0586, 0.0509, 0.0451, 0.0401, 0.0361,
           0.0301, 0.0260, 0.0231),
    # Nitrogen Z=7
    "N":  (0.1980, 0.1672, 0.1357, 0.1199, 0.1019, 0.0903, 0.0769,
           0.0683, 0.0623, 0.0573, 0.0498, 0.0442, 0.0394, 0.0357,
           0.0299, 0.0263, 0.0237),
    # Oxygen Z=8
    "O":  (0.2132, 0.1759, 0.1390, 0.1217, 0.1029, 0.0912, 0.0779,
           0.0687, 0.0626, 0.0577, 0.0501, 0.0445, 0.0395, 0.0358,
           0.0301, 0.0267, 0.0241),
    # Iron Z=26
    "Fe": (1.958,  1.205,  0.5952, 0.3717, 0.1964, 0.1460, 0.1099,
           0.0941, 0.0840, 0.0769, 0.0667, 0.0594, 0.0529, 0.0485,
           0.0427, 0.0394, 0.0379),
}


# ---------------------------------------------------------------------------
# Matrix compositions (mass fractions, sum = 1.0)
# ---------------------------------------------------------------------------

OISN_16_COMPOSITION: Dict[str, float] = {
    "H":  0.022,
    "C":  0.206,
    "N":  0.009,
    "O":  0.049,
    "Fe": 0.714,
}
"""LSRM ОИСН-16 standard soil simulant. Used as reference matrix
in .efr files for Marinelli/Дента-120мл/Петри-60мл geometries.
Mass fractions verified from .efr Material field (sum=1.000)."""


# Per-geometry reference (ρ, thickness) extracted from .efr Layers field
# of the LSRM Gamma-1С NaI 63×63 USB SN-01 detector calibration.
# Values verified 2024-11 by direct inspection of .efr metadata.
#
# Format: geometry_name → (volume_ml, ρ_ref_g_cm3, thickness_cm)
#
# thickness = effective sample column depth used in .efr calibration
# (Layers.Width field). For Marinelli the "thickness" is the wall
# thickness perpendicular to the detector axis; for Дента/Петри it's
# the dish/cylinder depth.
REF_GEOMETRY: Dict[str, tuple] = {
    # Only Маринелли is registered for K-20 external correction because
    # its .efr file has the outer Layers.Enable=false flag, meaning the
    # matrix self-attenuation is NOT already baked into the calibrated
    # ε(E). Дента-120мл and Петри-60мл .efr files have outer
    # Layers.Enable=true — matrix correction is already incorporated
    # into the calibration; applying K-20 externally would
    # double-correct and degrade accuracy.
    #
    # Format: geometry_name → (volume_ml, ρ_ref_g_cm3, thickness_cm).
    # Verified from .efr 2024-11.
    "Маринелли":   (1000.0, 1.60,  3.1),   # 1L water-equivalent, 31mm sleeve
    # "Дента-120мл": already corrected in .efr (Layers.Enable=true)
    # "Петри-60мл":  already corrected in .efr (Layers.Enable=true)
}


# ---------------------------------------------------------------------------
# Core physics
# ---------------------------------------------------------------------------

def _validate_energy(E_keV: float) -> None:
    if not math.isfinite(E_keV):
        raise ValueError(f"E_keV must be finite, got {E_keV!r}")
    if E_keV <= 0:
        raise ValueError(f"E_keV must be positive, got {E_keV}")


def element_mu_over_rho(E_keV: float, symbol: str) -> float:
    """Mass-attenuation coefficient μ/ρ (cm²/g) for a single element
    at energy E_keV, log-log interpolated from NIST XCOM tabulated
    values. Outside the table range [50, 3000] keV, clamps to nearest
    edge value.

    Args:
        E_keV: photon energy in keV
        symbol: element symbol (e.g. "H", "Fe"); must be in NIST XCOM table

    Returns:
        μ/ρ in cm²/g

    Raises:
        KeyError if `symbol` is not in the tabulated set
        ValueError if E_keV is non-positive or non-finite
    """
    _validate_energy(E_keV)
    if symbol not in _XCOM_MU_RHO:
        raise KeyError(
            f"element {symbol!r} not in XCOM table; "
            f"available: {sorted(_XCOM_MU_RHO)}"
        )
    energies = _XCOM_ENERGIES_KEV
    mu_rho_vals = _XCOM_MU_RHO[symbol]
    if E_keV <= energies[0]:
        return mu_rho_vals[0]
    if E_keV >= energies[-1]:
        return mu_rho_vals[-1]
    # Bracket: find i such that energies[i] <= E_keV < energies[i+1]
    for i in range(len(energies) - 1):
        if energies[i] <= E_keV <= energies[i + 1]:
            E1, E2 = energies[i], energies[i + 1]
            mu1, mu2 = mu_rho_vals[i], mu_rho_vals[i + 1]
            # log-log linear interpolation (NIST XCOM recommendation
            # for smooth photon cross sections)
            t = (math.log(E_keV) - math.log(E1)) / (math.log(E2) - math.log(E1))
            log_mu = math.log(mu1) + t * (math.log(mu2) - math.log(mu1))
            return math.exp(log_mu)
    raise RuntimeError(f"unreachable: E={E_keV} not bracketed by table")


def matrix_mu_over_rho(E_keV: float,
                       composition: Dict[str, float]) -> float:
    """Compute matrix-averaged μ/ρ at energy E_keV (cm²/g).

    For a multi-element matrix with mass fractions w_i,
        (μ/ρ)_matrix = Σ_i w_i × (μ/ρ)_i

    Args:
        E_keV: photon energy in keV
        composition: dict {element_symbol: mass_fraction_w_i}.
            Should sum to 1.0; caller responsible for normalization.

    Returns:
        Matrix-averaged μ/ρ in cm²/g
    """
    return sum(w * element_mu_over_rho(E_keV, sym)
               for sym, w in composition.items())


def slab_self_attenuation_factor(mu_over_rho_cm2_g: float,
                                 rho_g_cm3: float,
                                 thickness_cm: float) -> float:
    """Self-attenuation factor F(E) = (1 − exp(−μρt)) / (μρt) for a
    thin disk slab approximation.

    F ∈ (0, 1] for all physically meaningful inputs.
    F → 1 as μρt → 0 (thin/light/high-E; no attenuation).
    F → 1/(μρt) as μρt → ∞ (thick/heavy/low-E; small fraction escapes).

    Args:
        mu_over_rho_cm2_g: matrix μ/ρ in cm²/g (from matrix_mu_over_rho)
        rho_g_cm3: sample density in g/cm³
        thickness_cm: effective sample thickness in cm

    Returns:
        F ∈ (0, 1]
    """
    if rho_g_cm3 < 0 or thickness_cm < 0:
        raise ValueError(
            f"density and thickness must be non-negative; got "
            f"ρ={rho_g_cm3}, t={thickness_cm}"
        )
    x = mu_over_rho_cm2_g * rho_g_cm3 * thickness_cm
    if x < 1e-6:
        # Series expansion (1 - exp(-x))/x ≈ 1 - x/2 + x²/6 - ...
        # avoids loss-of-significance for very small x.
        return 1.0 - x / 2.0 + (x * x) / 6.0
    return (1.0 - math.exp(-x)) / x


def correction_factor(E_keV: float, *,
                      rho_sample_g_cm3: float,
                      rho_ref_g_cm3: float,
                      thickness_cm: float,
                      composition: Dict[str, float] = OISN_16_COMPOSITION,
                      ) -> float:
    """K-20 multiplicative correction factor F_ref/F_sample at energy
    E_keV. Apply to a `compute_activity` result to recover true
    activity:

        A_true = A_measured × correction_factor(E_keV, ...)

    Args:
        E_keV: line energy in keV
        rho_sample_g_cm3: actual sample density (g/cm³)
        rho_ref_g_cm3: .efr reference sample density (g/cm³)
        thickness_cm: effective sample thickness (cm), assumed same
            for both reference and sample (same container)
        composition: matrix mass fractions; defaults to ОИСН-16

    Returns:
        correction = F_ref(E) / F_sample(E)

        When ρ_sample < ρ_ref: correction < 1 (reduces over-estimate)
        When ρ_sample > ρ_ref: correction > 1 (boosts under-estimate)
        When ρ_sample = ρ_ref: correction = 1 (no correction)
    """
    if rho_ref_g_cm3 <= 0:
        # Point source / non-matrix calibration: no correction defined
        return 1.0
    mu_E_over_rho = matrix_mu_over_rho(E_keV, composition)
    F_sample = slab_self_attenuation_factor(
        mu_E_over_rho, rho_sample_g_cm3, thickness_cm)
    F_ref = slab_self_attenuation_factor(
        mu_E_over_rho, rho_ref_g_cm3, thickness_cm)
    if F_sample <= 0:
        # Pathological — sample so dense/thick that F → 0. Return 1
        # (no correction) and let caller flag the case.
        return 1.0
    return F_ref / F_sample


def weighted_mean_correction(E_keV_list: Iterable[float],
                             weights: Iterable[float], *,
                             rho_sample_g_cm3: float,
                             rho_ref_g_cm3: float,
                             thickness_cm: float,
                             composition: Dict[str, float] = OISN_16_COMPOSITION,
                             ) -> float:
    """Weighted mean of per-line correction factors. Useful when
    compute_activity produces an inverse-variance-weighted mean of
    line activities and the caller wants a single overall correction
    factor consistent with that weighting.

    Args:
        E_keV_list: line energies
        weights: weighting factors (typically 1/sigma² from
            compute_activity.lines_used); must be same length as
            E_keV_list

    Returns:
        Weighted mean correction. If sum(weights) == 0, returns 1.0
        (no correction).
    """
    E_keV_list = list(E_keV_list)
    weights = list(weights)
    if len(E_keV_list) != len(weights):
        raise ValueError(
            f"E_keV_list and weights length mismatch: "
            f"{len(E_keV_list)} vs {len(weights)}"
        )
    if not E_keV_list:
        return 1.0
    total_w = sum(weights)
    if total_w <= 0:
        return 1.0
    weighted_sum = sum(
        w * correction_factor(E, rho_sample_g_cm3=rho_sample_g_cm3,
                              rho_ref_g_cm3=rho_ref_g_cm3,
                              thickness_cm=thickness_cm,
                              composition=composition)
        for E, w in zip(E_keV_list, weights)
    )
    return weighted_sum / total_w
