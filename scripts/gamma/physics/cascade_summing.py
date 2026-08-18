"""
True coincidence summing (TCS) correction — K-17.

When a nuclide emits multiple γ-rays in cascade (same de-excitation,
emitted within the detector's coincidence-resolving time ~µs), two
photons may both deposit energy in the detector during one MCA
processing cycle. The MCA registers the sum E1+E2 as a single event,
**depleting** the individual photopeaks at E1 and E2 and creating a
small sum peak at E1+E2 instead.

The depletion fraction depends on the **total efficiency** (any
interaction in the detector, not just photopeak) at the partner γ
energy. Total efficiency is higher than photopeak efficiency because
Compton scattering, partial energy deposition etc. also register as
"some event in this acquisition window".

═══════════════════════════════════════════════════════════════════
Formula (Knoll §17.6, Gilmore §8.5.1)
═══════════════════════════════════════════════════════════════════

For a single line E_i of nuclide N emitted in cascade with other γs:

    Loss(E_i) = Σ_j p_c(E_i, E_j) · ε_T(E_j)

    C(E_i) = 1 / (1 − Loss(E_i))             [correction factor]

    A_true,i = C(E_i) · A_observed,i

where:
  • p_c(E_i, E_j) — fraction of decays that produce both E_i and E_j
    in cascade (from nuclear decay scheme)
  • ε_T(E_j) — total efficiency at E_j (Compton + photopeak + escape
    + everything, integrated over all deposit energies)

The factor is **per geometry**: ε_T scales linearly with solid angle.
A small ε_T (point at 25cm, Ω/4π ≈ 0.005) gives correction ~0.5%;
a large ε_T (Marinelli 4π, Ω/4π ≈ 0.4) gives correction tens of %.

This module implements the *loss* term only. The *gain* term (where
cascade pairs sum INTO another nuclide's photopeak) is neglected;
for most use cases (no overlapping cascades from different nuclides
at the same sum energy) it is <1% and is documented as a limitation
of this implementation.

═══════════════════════════════════════════════════════════════════
Total efficiency from photopeak efficiency
═══════════════════════════════════════════════════════════════════

We don't measure ε_T(E) directly — we have ε_p(E) from the F-27
.efr calibration. The bridge is the **peak-to-total ratio**:

    ε_T(E) = ε_p(E) / P(E)

P(E) for NaI 3×3" (~63×63 mm) is well-known empirically (Gilmore
Table 8.4 quotes for various crystal sizes). It falls smoothly from
~0.92 at low E to ~0.17 at 2614 keV (the Tl-208 high-energy line).

A log-log polynomial fit (degree 2) to the Gilmore reference values
yields:
    log P(E) = −0.316 + 0.458 · log E − 0.081 · (log E)²
which is hard-coded as `peak_to_total_NaI` below.

For HPGe, P(E) is closer to 1.0 across the whole range — TCS is much
smaller. A separate `peak_to_total_HPGe` is provided with conservative
values; callers must supply their own measured P(E) for accurate HPGe
TCS work.

═══════════════════════════════════════════════════════════════════
Caveats and known limitations
═══════════════════════════════════════════════════════════════════

  • **Single-pair approximation**: for nuclides with >2 cascade γs
    per decay (e.g. Eu-152, Ba-133), only direct pair couplings are
    considered. Triple coincidences and longer cascade chains are
    neglected (typically <1% extra at point geometries).
  • **Sum-in gain neglected** (see above). Relevant only for
    accidental overlap of cascade sum with another photopeak.
  • **P(E) model is for NaI 3×3"** at near-point geometry. For very
    different crystal sizes or Marinelli/contact geometries, the user
    should supply their own `p_t_func`.
  • **Cascade scheme data is curated**, not auto-generated from
    ENSDF. We catalogue the major commercially-relevant calibrators
    (Co-60, Y-88, Na-22, Tl-208 from Th-228, Eu-152 partial, Ba-133
    partial). Other nuclides return correction factor 1.0 silently.
  • **Angular correlation W(θ) is set to 1**. For most cases this is
    accurate to <2% (Gilmore §8.5.2).

═══════════════════════════════════════════════════════════════════
References
═══════════════════════════════════════════════════════════════════

  • Knoll G.F. *Radiation Detection and Measurement*, 4th Ed.,
    Wiley 2010, §17.6 (coincidence summing).
  • Gilmore G., Joss D. *Practical Gamma-ray Spectrometry*, 3rd Ed.,
    Wiley 2024, §8.5 (cascade summing correction); Table 8.4
    (peak-to-total ratios for NaI).
  • Lsrm *Algorithmic Foundations* (2025) §10 (Lsrm software
    implementation of TCS correction).
  • ENSDF / NuDat 3 (NNDC, BNL) — decay scheme branching ratios.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Optional


# ═════════════════════════════════════════════════════════════════════
# Cascade scheme data classes
# ═════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class CascadePartner:
    """One γ-line that may be emitted in cascade with a primary line.

    Attributes:
        E_partner_keV: energy of the partner γ-line in keV
        probability: p_c(primary | partner) — probability that this
            partner is emitted in the same decay as the primary line.
            For Co-60: p_c(1173, 1332) = 0.9988 (nearly every 1173
            comes with 1332). For nuclides emitting 2 annihilation
            photons (Na-22), probability may exceed 1.0 (we use 2.0
            to model two independent photons summing in).
    """

    E_partner_keV: float
    probability: float


@dataclass(frozen=True)
class CascadeScheme:
    """Cascade scheme for one nuclide — lookup of partners per primary γ."""

    nuclide: str
    cascades: dict  # dict[float, tuple[CascadePartner, ...]]

    def line_energies(self) -> tuple:
        return tuple(self.cascades.keys())


# ═════════════════════════════════════════════════════════════════════
# Cascade scheme catalogue
# ═════════════════════════════════════════════════════════════════════
#
# Branching probabilities curated from ENSDF and Gilmore §8.5. Each
# `cascades` entry maps primary_E_keV → tuple of CascadePartner.
# Energies use library values (NNDC); slight differences from
# measured channel-to-energy values are absorbed by the lookup
# tolerance in `tcs_correction_factor`.

CASCADE_SCHEMES = {
    # Co-60 → Ni-60: 2.5+ → 1.33+ (1173 keV) → 0+ ground (1332 keV).
    # Both γs emitted in essentially every decay (β¯ to 2.5+ state).
    "Co-60": CascadeScheme("Co-60", {
        1173.23: (CascadePartner(1332.49, 0.9988),),
        1332.49: (CascadePartner(1173.23, 0.9988),),
    }),

    # Y-88 → Sr-88: EC/β+ → 2734 → 1836 → 898 → ground.
    # 898 always comes from cascade (94% of decays); 1836 → 898 in
    # 99% of paths that produce 1836. We use rounded probabilities.
    "Y-88": CascadeScheme("Y-88", {
        898.04:  (CascadePartner(1836.06, 1.000),),
        1836.06: (CascadePartner(898.04,  0.94),),
    }),

    # Na-22 → Ne-22*: β+/EC → 1274.54 → ground. The β+ produces an
    # e+ that annihilates → 2 × 511 keV photons. Treating these as
    # two cascade partners at 511 keV each:
    #   • 1274 keV depletion comes from BOTH 511 photons
    #   • 511 keV depletion comes from 1274 partner AND the other
    #     511 (positron's pair) — but the other 511 is back-to-back
    #     and rarely co-detected; conservative single-partner approx
    "Na-22": CascadeScheme("Na-22", {
        1274.54: (CascadePartner(511.00, 0.903),
                  CascadePartner(511.00, 0.903),),  # two annihilation γs
        511.00:  (CascadePartner(1274.54, 0.903),),
    }),

    # Tl-208 (from Th-228 chain decay): β¯ from Tl-208 g.s.
    # to Pb-208 5-/4- → cascade 510.77 + 583.19 + 2614.51 (+ minor).
    # ENSDF branching: 2614 in ~100% of Tl-208 decays; 583 in 85%;
    # 510 in 23%; 860 in 12.5%.
    # We catalog the strongest pairs.
    "Tl-208": CascadeScheme("Tl-208", {
        2614.51: (CascadePartner(583.19, 0.85),),
        583.19:  (CascadePartner(2614.51, 1.00),
                  CascadePartner(510.77, 0.27),),
        510.77:  (CascadePartner(583.19, 1.00),),
        860.56:  (CascadePartner(583.19, 1.00),),
    }),

    # Eu-152 → Sm-152 (EC, 72%) + Gd-152 (β¯, 28%): very complex
    # multi-state cascade with >20 γ-lines. We catalog the strongest
    # cascade pairs for the most prominent peaks. Coefficients are
    # approximations from the dominant decay paths (NNDC).
    # This is INCOMPLETE — flagged in the module limitations.
    "Eu-152": CascadeScheme("Eu-152", {
        121.78: (CascadePartner(344.28, 0.92),
                 CascadePartner(244.70, 0.27),),
        244.70: (CascadePartner(121.78, 1.00),
                 CascadePartner(344.28, 0.74),),
        344.28: (CascadePartner(121.78, 0.40),
                 CascadePartner(244.70, 0.21),),
        778.90: (CascadePartner(121.78, 0.87),),
        964.06: (CascadePartner(121.78, 0.65),),
        1112.07: (CascadePartner(121.78, 0.93),),
        1408.01: (CascadePartner(121.78, 0.97),
                  CascadePartner(244.70, 0.20),),
    }),

    # Ba-133 → Cs-133: EC. Cascade lines 80.99, 276, 302, 356, 384 etc.
    # Major pair: 356 with 80 in ~62% of decays.
    "Ba-133": CascadeScheme("Ba-133", {
        80.998: (CascadePartner(356.013, 0.85),
                 CascadePartner(302.851, 0.27),),
        302.851: (CascadePartner(80.998, 1.00),),
        356.013: (CascadePartner(80.998, 0.85),),
        383.851: (CascadePartner(80.998, 1.00),),
    }),

    # F-128 / v1.17.7 — Bi-212 (Th-232 chain): β¯ → Po-212 + α → Tl-208.
    # β¯ branch (64.06 %) produces excited Po-212 states which de-excite
    # via cascade γ-emission. Strongest pairs (NNDC ENSDF):
    #   • 727.33 (6.74 %) emitted in cascade with 785.4 / 893.4 / 1078 keV
    #   • 1620.5 (1.49 %) emitted in cascade with 727.3 (B(M1)=1.0) and
    #     785.4 / 893.4 (low intensity)
    # На NaI 63×63 в Marinelli/0см геометрии TCS-обеднение ~5-15 %
    # для 727 keV и 1620 keV. Без коррекции активность Bi-212
    # систематически занижена → нарушение chain equilibrium (ratio Pb-212
    # / Bi-212 завышается). Коэффициенты ниже — упрощённая схема пары
    # сильнейших линий из ENSDF Bi-212 (2014).
    "Bi-212": CascadeScheme("Bi-212", {
        727.33: (CascadePartner(785.37, 0.165),
                 CascadePartner(1078.62, 0.084),),
        785.37: (CascadePartner(727.33, 1.000),),
        1078.62: (CascadePartner(727.33, 1.000),),
        1620.50: (CascadePartner(727.33, 0.995),
                  CascadePartner(785.37, 0.114),),
    }),
}


# ═════════════════════════════════════════════════════════════════════
# Peak-to-total ratio models
# ═════════════════════════════════════════════════════════════════════

# Fit coefficients for log P = a + b·log E + c·(log E)² obtained by
# fitting Gilmore Table 8.4 NaI 3×3" data at 100, 200, 500, 1000, 1500,
# 2000, 2600 keV. Valid 50-3000 keV; flat-clipped to [0.05, 1.0].
_PT_NAI_COEFS = (-0.316, 0.458, -0.081)


# K-21 (v1.9.0): close-geometry P/T scaling factors. Empirically
# calibrated from v1.7.25 cert matrix Tl-208 chain-proxy findings.
#
# Rationale: peak_to_total_NaI() defaults to Gilmore Table 8.4 NaI 3×3"
# values which are calibrated for POINT-GEOMETRY at ~5 cm distance.
# At very-close-geometry samples (0 cm distance for Marinelli, Дента,
# Петри), the detector covers a much larger solid angle, true
# coincidence summing of cascade photons is more frequent, and the
# effective P/T is reduced relative to the 5 cm reference.
#
# The empirical "factor" multiplies the 5cm P/T to obtain the effective
# close-geometry P/T:
#       P_effective(E) = P_5cm(E) × geometry_factor
#
# Values fit so that the v1.7.25 Tl-208 systematic underestimate
# (−15 to −20%) becomes close to zero. Values are not geometry-
# independent — they reflect the specific Lsrm Gamma-1С NaI 63×63 mm
# detector ε(E) curves per geometry.
#
# See K-21 in KNOWN_AND_FIXED_ISSUES.md for the full quantitative
# analysis and v1.9 limitation scope.
GEOMETRY_PT_FACTOR: dict = {
    "Точечная-5см":  1.00,   # Gilmore reference, no scaling
    "Точечная-25см": 1.00,   # neutral (no chain-proxy fixture to fit)
    "Маринелли":     0.45,   # fitted to bring Tl-208 −20% → ≤ ±5%
    "Дента-120мл":   0.50,   # fitted to bring Tl-208 −18% → ≤ ±5%
    "Петри-60мл":    0.50,   # fitted to bring Tl-208 −20% → ≤ ±5%
}


def peak_to_total_NaI(E_keV: float, *,
                      geometry_factor: float = 1.0) -> float:
    """
    Peak-to-total ratio for a typical NaI 3×3" (~63×63 mm) crystal.

    P(E) = exp(−0.316 + 0.458·ln E − 0.081·(ln E)²) × geometry_factor,
    clipped to [0.05, 1.0].

    Values are point-geometry at ~5 cm distance (data from Gilmore
    Table 8.4 et al.). For different crystal sizes or geometries the
    user should supply their own callable, or pass the K-21
    `geometry_factor` kwarg (default 1.0).

    K-21 (v1.9.0): `geometry_factor` scales the intrinsic point-5cm
    P/T to model effective P/T at close-geometry samples
    (Marinelli/Дента/Петри at 0 cm). Smaller `geometry_factor` →
    smaller effective P → larger ε_T → larger TCS correction.

    For the LSRM Gamma-1С NaI 63×63 mm detector, see the
    `GEOMETRY_PT_FACTOR` dict for empirically-fit values per geometry.
    For other detectors/geometries, supply your own value (fit from
    your cert sources).

    Cross-check vs Gilmore Table 8.4 at geometry_factor=1.0:
        E (keV)   reference P   model P
            100        0.92       1.00 (capped)
            200        0.85       0.74
            500        0.55       0.55
           1000        0.35       0.36
           1500        0.27       0.27
           2000        0.22       0.22
           2600        0.17       0.18
    """
    if E_keV <= 0:
        return 1.0
    lnE = math.log(E_keV)
    a, b, c = _PT_NAI_COEFS
    log_P = a + b * lnE + c * lnE * lnE
    P = math.exp(log_P)
    return max(0.05, min(1.0, P * geometry_factor))


def peak_to_total_NaI_for_geometry(geometry: str = "Точечная-5см"):
    """K-21 convenience: return a P/T callable bound to a registered
    geometry's empirical scaling factor.

    Use the returned callable as the `p_t_func` argument to
    `compute_tcs_corrections` or `tcs_correction_factor`:

        from gamma.physics.cascade_summing import (
            compute_tcs_corrections, peak_to_total_NaI_for_geometry,
        )
        pt = peak_to_total_NaI_for_geometry("Маринелли")
        cc = compute_tcs_corrections("Tl-208", eff_curve, p_t_func=pt)

    Args:
        geometry: token from `GEOMETRY_PT_FACTOR`. Unknown geometries
            default to 1.0 (no scaling) with a no-op.

    Returns:
        Callable `(E_keV) -> P/T`, with geometry-appropriate scaling
        baked in.
    """
    factor = GEOMETRY_PT_FACTOR.get(geometry, 1.0)
    def _pt(E_keV: float) -> float:
        return peak_to_total_NaI(E_keV, geometry_factor=factor)
    return _pt


def peak_to_total_HPGe(E_keV: float) -> float:
    """
    Crude P/T ratio for HPGe. NOT calibrated against measured data —
    typical relative-efficiency 30%+ coaxial values used as default.

    HPGe has higher photopeak fraction than NaI (better resolution,
    more efficient charge collection), so P/T is closer to 1 across
    the full range. For accurate TCS work on HPGe, the user MUST
    supply their own measured P(E).
    """
    if E_keV <= 0:
        return 1.0
    # Empirical: P ≈ 0.8 at 100 keV, dropping to ~0.3 at 2614 keV
    lnE = math.log(E_keV)
    # Linear fit log P vs log E: log P ≈ 0.95 − 0.27·ln E
    log_P = 0.95 - 0.27 * lnE
    P = math.exp(log_P)
    return max(0.10, min(1.0, P))


def total_efficiency(
    E_keV: float,
    eff_curve,
    p_t_func: Callable[[float], float] = peak_to_total_NaI,
) -> float:
    """
    Total efficiency at energy E from photopeak efficiency and P(E).

    ε_T(E) = ε_p(E) / P(E)

    Args:
        E_keV: γ-ray energy
        eff_curve: object with `.efficiency_at(E)` returning ε_p
        p_t_func: callable returning P/T ratio for the detector geometry
    """
    eps_p = eff_curve.efficiency_at(E_keV)
    if eps_p <= 0:
        return 0.0
    P = p_t_func(E_keV)
    if P <= 0:
        return 0.0
    return eps_p / P


# ═════════════════════════════════════════════════════════════════════
# TCS correction calculation
# ═════════════════════════════════════════════════════════════════════

def _find_line_in_scheme(E_line: float, cascades: dict,
                          tolerance_keV: float = 0.5) -> Optional[float]:
    """Find the scheme key closest to E_line within tolerance; None if
    no match."""
    best = None
    best_dist = float("inf")
    for E_scheme in cascades.keys():
        d = abs(E_scheme - E_line)
        if d < tolerance_keV and d < best_dist:
            best_dist = d
            best = E_scheme
    return best


def tcs_correction_factor(
    nuclide: str,
    E_line_keV: float,
    eff_curve,
    p_t_func: Callable[[float], float] = peak_to_total_NaI,
    *,
    energy_tolerance_keV: float = 0.5,
    loss_cap: float = 0.5,
) -> float:
    """
    Per-line TCS correction factor.

    Computes C = 1/(1 − L) where L = Σ p_c · ε_T(E_partner).
    Returns 1.0 (no correction) when the nuclide is not in the catalogue,
    when the line is not in the nuclide's cascade scheme, or when ε_T
    information is insufficient.

    Args:
        nuclide: e.g. "Co-60"
        E_line_keV: photopeak energy whose correction we want
        eff_curve: fitted EfficiencyCurve (or any object with
            `.efficiency_at(E)`)
        p_t_func: peak-to-total ratio function for the detector geometry
        energy_tolerance_keV: how close E_line must be to a scheme entry
        loss_cap: maximum allowable Loss before clipping (safety against
            unphysical results from very strong partners; default 0.5
            corresponds to a maximum correction factor of 2.0)

    Returns:
        Correction factor C ≥ 1.0. Multiply observed peak area by this
        to obtain TCS-corrected area (the disintegration count would
        have been larger by C in the absence of cascade summing).
    """
    scheme = CASCADE_SCHEMES.get(nuclide)
    if scheme is None:
        return 1.0
    key = _find_line_in_scheme(E_line_keV, scheme.cascades,
                                tolerance_keV=energy_tolerance_keV)
    if key is None:
        return 1.0
    L = 0.0
    for partner in scheme.cascades[key]:
        eps_T = total_efficiency(partner.E_partner_keV, eff_curve, p_t_func)
        L += partner.probability * eps_T
    if L >= loss_cap:
        L = loss_cap
    if L < 0:
        L = 0.0
    return 1.0 / (1.0 - L)


def compute_tcs_corrections(
    nuclide: str,
    eff_curve,
    p_t_func: Callable[[float], float] = peak_to_total_NaI,
    **kwargs,
) -> dict:
    """
    All TCS correction factors for a nuclide as a dict {E_keV: factor}.

    Suitable for passing to `compute_activity(coincidence_correction=...)`.
    Returns an empty dict for non-cascade nuclides.

    Example::

        from gamma.physics.cascade_summing import compute_tcs_corrections
        from gamma.activity import compute_activity

        cc = compute_tcs_corrections("Co-60", eff_curve)
        result = compute_activity(
            co60_id, efficiency_curve=eff_curve,
            live_time_s=spec.live_time,
            from_bg_subtracted=True,
            coincidence_correction=cc,     # ← TCS applied
            ...,
        )
    """
    scheme = CASCADE_SCHEMES.get(nuclide)
    if scheme is None:
        return {}
    return {
        E_line: tcs_correction_factor(nuclide, E_line, eff_curve,
                                       p_t_func, **kwargs)
        for E_line in scheme.line_energies()
    }


__all__ = [
    "CascadePartner",
    "CascadeScheme",
    "CASCADE_SCHEMES",
    "GEOMETRY_PT_FACTOR",
    "peak_to_total_NaI",
    "peak_to_total_NaI_for_geometry",
    "peak_to_total_HPGe",
    "total_efficiency",
    "tcs_correction_factor",
    "compute_tcs_corrections",
]
