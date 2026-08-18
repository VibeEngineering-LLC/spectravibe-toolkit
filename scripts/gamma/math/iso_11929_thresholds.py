"""
ISO 11929 decision threshold and detection limit — per-peak Bq/kg emit.

F-RPT-08 / v1.19.1 (wave 4).  #PTB-1 (2026-07-02): added IMIS/KTA regime.

Implements two functions for per-peak ISO 11929 characteristic limits
expressed in Bq/kg, suitable for direct inclusion in the JSON report.

Regulatory regimes (PTB-2018 SPEKT/GRUNDL Annex C)
--------------------------------------------------
* **KTA** (Kerntechnischer Ausschuss, German nuclear-safety authority) —
  symmetric confidence: k_{1-α} = k_{1-β} = 1.645 (α = β = 0.05).
  Default in this module (matches the wave-4 F-RPT-08 baseline).
* **IMIS** (Integriertes Mess- und Informationssystem, German environmental
  radioactivity monitoring) — asymmetric confidence:
  k_{1-α} = 3.0 (≈ 99.87 % → α ≈ 0.00135),
  k_{1-β} = 1.645 (β = 0.05).
  Selected via `regime="IMIS"`.

ISO 11929 references
--------------------
* ISO 11929-1:2019(E) §3.12 (Decision threshold, y*) — cited in RAG-005 source
  iso_11929_1_2019; §3.13 (Detection limit, y#).
* ISO 11929-1:2019(E) §5.4.3 (formula for decision threshold):
    y* = k_{1-α} · ũ(0)
  where ũ(0) is the standard uncertainty of the measurand estimator
  evaluated at the true value ỹ = 0.
* ISO 11929-1:2019(E) §5.4.4 (formula for detection limit):
    y# is the smallest y satisfying  y# = y* + k_{1-β} · ũ(y#)
  For k_{1-α} = k_{1-β} (symmetric, α = β = 0.05) and under the
  approximation ũ(y#) ≈ ũ(0) (low-counting-statistics regime where the
  Poisson variance at the detection limit is dominated by background,
  not signal), the iterative equation simplifies to:
    y# ≈ 2 · y*
  This approximation is explicitly noted as an approximation valid for
  low-count / background-dominated cases; the full iterative solution is
  out of scope for this wave (see tech_debt note below).

RAG cross-references
---------------------
* RAG-005 — Background quality control gates (Currie/ISO 11929 lineage)
* RAG-008 — ROI-windowed flux drift gate
* RAG-009 — Consolidated 13-source cite list for spectrum QC
* RAG-022 — Poisson |z|-test gate (BUG-35; ISO 11929-2:2019 §6)

Formula derivation (σ_0(0) in Bq/kg units)
---------------------------------------------
For a single-line gamma measurement, the measurand y is specific activity
(Bq/kg).  The sensitivity coefficient w that converts net counts → Bq/kg is:

    w = 1 / (ε · I · M · t)

where ε = efficiency (dimensionless), I = branching ratio (decimal),
M = sample mass (kg), t = live time (s).

Under the null hypothesis y = 0, the uncertainty of the net-counts
estimator is:

    u(n_net | y=0) = sqrt(N_gross + N_bg)          [counts, Poisson]

(ISO 11929-1:2019 §5.4.3: when subtracting a background measurement
of the same live time, the variance of the estimator at y=0 is
Var(N_gross) + Var(N_bg) = N_gross + N_bg by Poisson counting.)

Multiplying by w converts to Bq/kg:

    ũ(0) = sqrt(N_gross + N_bg) / (ε · I · M · t)   [Bq/kg]

Then:
    y* = k_{1-α} · ũ(0)
    y# ≈ 2 · y*    (low-stats approximation, α = β = 0.05)

Tech debt
----------
Full iterative y# solution per ISO 11929-1:2019 §5.4.4 is out of scope
for wave 4 (F-RPT-08).  The approximation y# ≈ 2·y* is exact only when
k_{1-α} = k_{1-β} AND ũ(y#) = ũ(0) (i.e., background-dominated regime).
For signal-dominated peaks (N_gross >> N_bg), the full iterative solution
can yield y# > 2·y* by up to ~30%.  Track in KNOWN_AND_FIXED_ISSUES.md
or roadmap as F-RPT-08.1.
"""

from __future__ import annotations

import math
from typing import Optional, Sequence


# ──────────────────────────────────────────────────────────────────
# k-quantile of the standard normal distribution
# (1-α quantile: P(Z > k) = α)
# ISO 11929 Table A.1 / ISO 11929-1:2019 §A.1; α = 0.05 → k = 1.6449
# Commonly rounded to 1.645 in applied spectrometry literature
# (Gilmore & Joss, Practical Gamma-ray Spectrometry, 3rd Ed §5.6;
# Knoll, Radiation Detection and Measurement, 4th Ed §3.II.B).
# ──────────────────────────────────────────────────────────────────
_K_95 = 1.6449  # one-sided 95% normal quantile, α = β = 0.05
_K_IMIS_ALPHA = 3.0  # PTB-2018 Annex C — IMIS k_{1-α} = 3 (≈ 99.87 %)
_K_IMIS_BETA = 1.6449  # PTB-2018 Annex C — IMIS k_{1-β} = 1.645
_VALID_REGIMES = ("KTA", "IMIS")


def _regime_k_pair(regime: str) -> Optional[tuple]:
    """Return (k_alpha, k_beta) for a regulatory regime, or None if invalid."""
    if regime == "KTA":
        return (_K_95, _K_95)
    if regime == "IMIS":
        return (_K_IMIS_ALPHA, _K_IMIS_BETA)
    return None


def decision_threshold(
    gross_counts: float,
    bg_counts: float,
    efficiency: float,
    branching_ratio: float,
    mass_kg: float,
    live_time_s: float,
    alpha: float = 0.05,
    regime: Optional[str] = None,
) -> Optional[float]:
    """Decision threshold y* in Bq/kg (ISO 11929-1:2019 §5.4.3).

    y* = k_{1-α} · ũ(0)

    where ũ(0) = sqrt(gross_counts + bg_counts) / (ε · I · M · t)

    All inputs are per-peak quantities integrated over the peak ROI
    (typically ±1 FWHM).

    Parameters
    ----------
    gross_counts : float
        Gross counts in the peak ROI from the sample spectrum (Poisson).
    bg_counts : float
        Background counts in the same ROI from the background spectrum,
        already scaled to the sample live time (Poisson).
    efficiency : float
        Absolute photopeak detection efficiency ε (dimensionless, 0–1).
        Typical NaI 50×50 mm at 10 cm: ε(661 keV) ≈ 0.02.
    branching_ratio : float
        Gamma-line emission probability I per nuclear decay (decimal, 0–1).
        Example: Cs-137 main line = 0.851, K-40 = 0.1066.
    mass_kg : float
        Sample mass in kg.
    live_time_s : float
        Sample live time in seconds.
    alpha : float, optional
        False-positive probability (default 0.05 → 5%).
        k_{1-α} is computed via scipy if available; for α = 0.05 the
        hardcoded constant _K_95 = 1.6449 is used for stability.
        Ignored when `regime` is set.
    regime : str, optional
        Regulatory regime override (PTB-2018 SPEKT/GRUNDL Annex C):
        - "KTA" → k_{1-α} = 1.6449 (equivalent to alpha=0.05, default).
        - "IMIS" → k_{1-α} = 3.0 (~99.87%, stricter environmental monitoring).
        When set, overrides `alpha`. Invalid regime → None.

    Returns
    -------
    float or None
        Decision threshold in Bq/kg, or None if any guard condition fires:
        - efficiency ≤ 0
        - branching_ratio ≤ 0
        - mass_kg ≤ 0
        - live_time_s ≤ 0
        - regime not in ("KTA", "IMIS")

    Notes
    -----
    Formula derivation: see module docstring.
    ISO 11929-1:2019 §5.4.3: y* = k_{1-α} · ũ(0).
    ũ(0) derived from Poisson counting variance at null hypothesis y = 0:
        u(n_net | y=0) = sqrt(N_gross + N_bg)   [ISO 11929 §A.3 example]
    Cross-references: RAG-005 (Currie/ISO lineage), RAG-009 (cite list),
    PTB-2018 SPEKT/GRUNDL Annex C (IMIS/KTA regimes).
    """
    # Guard: cannot compute sensible result if any denominator is zero.
    if efficiency <= 0 or branching_ratio <= 0 or mass_kg <= 0 or live_time_s <= 0:
        return None

    # Choose k quantile
    if regime is not None:
        k_pair = _regime_k_pair(regime)
        if k_pair is None:
            return None
        k = k_pair[0]
    elif alpha == 0.05:
        k = _K_95
    else:
        # For non-standard alpha: use approximation or caller-supplied value.
        # We compute an approximate k via the error function inverse.
        # Accurate to ±0.001 for alpha in [0.001, 0.5].
        k = _normal_quantile(1.0 - alpha)
        if k is None:
            return None

    # σ_0(0) in Bq/kg
    sigma_0 = math.sqrt(max(0.0, gross_counts) + max(0.0, bg_counts)) / (
        efficiency * branching_ratio * mass_kg * live_time_s
    )

    return k * sigma_0


def detection_limit(
    gross_counts: float,
    bg_counts: float,
    efficiency: float,
    branching_ratio: float,
    mass_kg: float,
    live_time_s: float,
    alpha: float = 0.05,
    beta: float = 0.05,
    regime: Optional[str] = None,
) -> Optional[float]:
    """Detection limit y# in Bq/kg (ISO 11929-1:2019 §5.4.4, low-stats approx).

    Under the low-counting-statistics approximation (background-dominated,
    ũ(y#) ≈ ũ(0)) and for symmetric decision probabilities (α = β = 0.05,
    k_{1-α} = k_{1-β} = 1.6449):

        y# ≈ 2 · y*

    where y* is the decision threshold from `decision_threshold()`.

    This approximation is exact when the Poisson variance at the detection
    limit is indistinguishable from the variance at zero signal — valid when
    N_bg >> expected net signal at y = y*. For signal-dominated peaks the
    full iterative solution per ISO 11929-1:2019 §5.4.4 is required (out of
    scope for F-RPT-08 wave 4; see tech_debt in module docstring).

    Parameters
    ----------
    gross_counts, bg_counts, efficiency, branching_ratio, mass_kg, live_time_s
        Same as `decision_threshold`.
    alpha : float, optional
        False-positive probability for the decision threshold (default 0.05).
        Ignored when `regime` is set.
    beta : float, optional
        False-negative probability for the detection limit (default 0.05).
        When alpha == beta, the symmetric approximation y# ≈ 2 · y* is used.
        When alpha != beta, k_{1-α} + k_{1-β} is used: y# ≈ (k_{1-α} + k_{1-β}) · ũ(0)
        (still under the ũ(y#) ≈ ũ(0) approximation).
        Ignored when `regime` is set.
    regime : str or None, optional
        Regulatory regime overriding (alpha, beta):
          * "KTA"  → k_{1-α} = k_{1-β} = 1.6449 (α = β = 0.05, symmetric).
          * "IMIS" → k_{1-α} = 3.0, k_{1-β} = 1.6449 (asymmetric per
            PTB-2018 SPEKT/GRUNDL Annex C).
        None (default) → use (alpha, beta) arguments.

    Returns
    -------
    float or None
        Detection limit in Bq/kg, or None if any guard condition fires:
        - efficiency, branching_ratio, mass_kg, live_time_s ≤ 0
        - regime not in ('KTA', 'IMIS')
        - alpha or beta not in (0, 1) when regime is None

    Notes
    -----
    ISO 11929-1:2019 §5.4.4: y# = smallest y satisfying
        y# = y* + k_{1-β} · ũ(y#)
    Approximation: ũ(y#) ≈ ũ(0) → y# ≈ y* + k_{1-β} · ũ(0)
    When k_{1-α} = k_{1-β}: y# ≈ 2 · k_{1-α} · ũ(0) = 2 · y*.
    Cross-references: RAG-005, RAG-008, RAG-009, RAG-022.
    """
    # Guard: same as decision_threshold
    if efficiency <= 0 or branching_ratio <= 0 or mass_kg <= 0 or live_time_s <= 0:
        return None

    # Compute k_{1-alpha} and k_{1-beta}
    if regime is not None:
        k_pair = _regime_k_pair(regime)
        if k_pair is None:
            return None
        k_alpha, k_beta = k_pair
    else:
        k_alpha = _K_95 if alpha == 0.05 else _normal_quantile(1.0 - alpha)
        k_beta = _K_95 if beta == 0.05 else _normal_quantile(1.0 - beta)
        if k_alpha is None or k_beta is None:
            return None

    # σ_0(0) in Bq/kg
    sigma_0 = math.sqrt(max(0.0, gross_counts) + max(0.0, bg_counts)) / (
        efficiency * branching_ratio * mass_kg * live_time_s
    )

    # y# ≈ (k_alpha + k_beta) · sigma_0   [ISO 11929 §5.4.4 low-stats approx]
    return (k_alpha + k_beta) * sigma_0


def _normal_quantile(p: float) -> Optional[float]:
    """Approximate upper quantile of N(0,1): P(Z ≤ z) = p.

    Uses rational approximation (Abramowitz & Stegun 26.2.17).
    Accurate to ±5×10⁻⁴ for 0.001 ≤ p ≤ 0.999.

    Returns None if p is out of (0, 1).
    """
    if not (0 < p < 1):
        return None
    # Use scipy.stats if available for better accuracy
    try:
        from scipy.stats import norm as _norm  # noqa: PLC0415
        return float(_norm.ppf(p))
    except ImportError:
        pass
    # Fallback: A&S 26.2.17 rational approximation for upper tail
    # t = sqrt(-2 * ln(1 - p)) for p > 0.5; mirror for p <= 0.5
    if p <= 0.5:
        q = p
        sign = -1.0
    else:
        q = 1.0 - p
        sign = 1.0
    t = math.sqrt(-2.0 * math.log(q))
    c = (2.515517, 0.802853, 0.010328)
    d = (1.432788, 0.189269, 0.001308)
    z = t - (c[0] + c[1]*t + c[2]*t*t) / (1 + d[0]*t + d[1]*t*t + d[2]*t*t*t)
    return sign * z


def multi_line_decision_threshold(
    line_gross_counts: Sequence[float],
    line_bg_counts: Sequence[float],
    line_efficiencies: Sequence[float],
    line_branching_ratios: Sequence[float],
    mass_kg: float,
    live_time_s: float,
    alpha: float = 0.05,
    regime: Optional[str] = None,
) -> Optional[float]:
    """Multi-line decision threshold y* in Bq/kg — PTB-2018 Annex C Eq. (C3).

    #PTB-2 (2026-07-02).  For emitters with several gamma peaks (e.g. Cs-134
    with 605 + 796 keV) the decision threshold combines all lines by
    inverse-variance weighting (PTB-2018 SPEKT/GRUNDL, p. γ-SPEKT/GRUNDL-54):

        y* = k_{1-α} · u(A = 0),   u(A = 0) = [ Σ_j 1/u_j²(0) ]^{-1/2}

    where the per-line zero-activity uncertainty uses the same counts-based
    variance model as `decision_threshold`:

        u_j(0) = sqrt(N_gross,j + N_bg,j) / (ε_j · p_j · m · t_live)

    (This is Eq. C3 with the PTB calibration factor φ_j = 1/(ε_j·p_j·m·t) and
    the background factor f_B absorbed into the explicit N_gross + N_bg
    variance — identical convention to the single-line functions above.)

    Parameters
    ----------
    line_gross_counts, line_bg_counts : sequence of float
        Per-line gross counts in the peak region and background (continuum)
        counts under the peak.  Equal length, one entry per gamma line.
    line_efficiencies, line_branching_ratios : sequence of float
        Per-line full-energy-peak efficiency ε_j and emission probability p_j.
    mass_kg, live_time_s : float
        Sample mass and live time (shared by all lines — same measurement).
    alpha : float, optional
        False-positive probability (default 0.05).  Ignored when `regime` set.
    regime : str or None, optional
        "KTA" | "IMIS" — see `decision_threshold`.

    Returns
    -------
    float or None
        Combined decision threshold in Bq/kg; 0.0 if any line has zero total
        counts (zero variance dominates); None on guard failure:
        - empty input or mismatched sequence lengths
        - mass_kg or live_time_s ≤ 0; any ε_j or p_j ≤ 0
        - invalid regime / alpha
    """
    n = len(line_gross_counts)
    if n == 0:
        return None
    if not (len(line_bg_counts) == len(line_efficiencies)
            == len(line_branching_ratios) == n):
        return None
    if mass_kg <= 0 or live_time_s <= 0:
        return None

    # Choose k quantile (same logic as decision_threshold)
    if regime is not None:
        k_pair = _regime_k_pair(regime)
        if k_pair is None:
            return None
        k = k_pair[0]
    elif alpha == 0.05:
        k = _K_95
    else:
        k = _normal_quantile(1.0 - alpha)
        if k is None:
            return None

    inv_var_sum = 0.0
    for gross, bg, eff, br in zip(
        line_gross_counts, line_bg_counts, line_efficiencies, line_branching_ratios
    ):
        if eff <= 0 or br <= 0:
            return None
        denom = eff * br * mass_kg * live_time_s
        var_j = (max(0.0, gross) + max(0.0, bg)) / (denom * denom)
        if var_j == 0.0:
            # Zero-variance line → combined σ(0) = 0 (consistent with the
            # single-line zero-counts case, which returns 0.0).
            return 0.0
        inv_var_sum += 1.0 / var_j

    sigma_0 = math.sqrt(1.0 / inv_var_sum)
    return k * sigma_0


def multi_line_detection_limit(
    line_gross_counts: Sequence[float],
    line_bg_counts: Sequence[float],
    line_efficiencies: Sequence[float],
    line_branching_ratios: Sequence[float],
    mass_kg: float,
    live_time_s: float,
    alpha: float = 0.05,
    beta: float = 0.05,
    regime: Optional[str] = None,
    measured_activity: Optional[float] = None,
    measured_uncertainty: Optional[float] = None,
    max_iterations: int = 100,
    rel_tolerance: float = 1e-9,
) -> Optional[float]:
    """Multi-line detection limit y# in Bq/kg — PTB-2018 Annex C Eq. (C4)+(C5).

    #PTB-2 (2026-07-02).  Iterative detection limit for multi-peak emitters
    (PTB-2018 SPEKT/GRUNDL, p. γ-SPEKT/GRUNDL-54):

        y# ≈ y* + k_{1-β} · u(y#')                                   (C4)
        u(y#') = sqrt( u²(0) + [u²(a_r) − u²(0)] · y#'/a_r )         (C5)

    with u(0) = y*/k_{1-α} the combined zero-activity uncertainty from
    Eq. (C3) and u(a_r) the standard uncertainty of the *measured* specific
    activity a_r.  Eq. (C5) linearly interpolates the variance between zero
    activity and the measured activity; the fixed point of (C4) is found by
    direct iteration starting from the ũ(y#) ≈ ũ(0) approximation.

    When `measured_activity` / `measured_uncertainty` are not available
    (activity not quantified), falls back to the non-iterative low-stats
    approximation  y# = y* + k_{1-β} · u(0)  — the same approximation used by
    the single-line `detection_limit`.

    Parameters
    ----------
    line_gross_counts, line_bg_counts, line_efficiencies,
    line_branching_ratios, mass_kg, live_time_s
        Same as `multi_line_decision_threshold`.
    alpha, beta : float, optional
        Error probabilities (defaults 0.05).  Ignored when `regime` set.
    regime : str or None, optional
        "KTA" | "IMIS" — see `decision_threshold`.
    measured_activity : float or None, optional
        Measured specific activity a_r in Bq/kg (must be > 0 to enable the
        iterative Eq. C5 variance interpolation).
    measured_uncertainty : float or None, optional
        Standard uncertainty u(a_r) in Bq/kg of the measured activity.
    max_iterations : int, optional
        Iteration cap for the (C4) fixed point (default 100).
    rel_tolerance : float, optional
        Relative convergence tolerance on y# (default 1e-9).

    Returns
    -------
    float or None
        Detection limit in Bq/kg, or None on guard failure (same guards as
        `multi_line_decision_threshold`) or if the iteration does not
        converge / the interpolated variance goes negative (pathological
        u(a_r) << u(0) at large extrapolation).
    """
    y_star = multi_line_decision_threshold(
        line_gross_counts, line_bg_counts,
        line_efficiencies, line_branching_ratios,
        mass_kg, live_time_s, alpha=alpha, regime=regime,
    )
    if y_star is None:
        return None

    # k_{1-alpha}, k_{1-beta}
    if regime is not None:
        k_pair = _regime_k_pair(regime)
        if k_pair is None:
            return None
        k_alpha, k_beta = k_pair
    else:
        k_alpha = _K_95 if alpha == 0.05 else _normal_quantile(1.0 - alpha)
        k_beta = _K_95 if beta == 0.05 else _normal_quantile(1.0 - beta)
        if k_alpha is None or k_beta is None:
            return None

    u_0 = y_star / k_alpha  # combined σ(0) from Eq. (C3)

    if (
        measured_activity is None
        or measured_uncertainty is None
        or measured_activity <= 0
        or measured_uncertainty < 0
    ):
        # Fallback: ũ(y#) ≈ ũ(0) low-stats approximation
        return y_star + k_beta * u_0

    # Iterative Eq. (C4) fixed point with Eq. (C5) variance interpolation
    u0_sq = u_0 * u_0
    um_sq = measured_uncertainty * measured_uncertainty
    y = y_star + k_beta * u_0  # start from the ũ(0) approximation
    for _ in range(max_iterations):
        u_sq = u0_sq + (um_sq - u0_sq) * (y / measured_activity)
        if u_sq < 0.0:
            return None
        y_next = y_star + k_beta * math.sqrt(u_sq)
        if abs(y_next - y) <= rel_tolerance * max(abs(y_next), 1e-30):
            return y_next
        y = y_next
    return None


__all__ = [
    "decision_threshold",
    "detection_limit",
    "multi_line_decision_threshold",
    "multi_line_detection_limit",
]
