"""
Per-nuclide activity calculation from identified γ-lines.

Implements the activity formula (Lsrm §8.4, Gilmore §5.7):

    A_i = (S_net,i · c_i) / (ε(E_i) · I_i · t_live)             [Bq]

per-line, where c_i is an optional correction factor (coincidence
summing, geometry, etc.). For multi-line nuclides the per-line
activities are combined by weighted averaging:

    A = Σ w_i · A_i / Σ w_i
    σ(A) = (Σ w_i)^(-1/2)
    w_i = 1 / σ²(A_i)

Uncertainty per line is propagated as (Gilmore §5.7.2):

    σ²(A_i) / A_i² = (σ_S/S)² + (σ_ε/ε)² + (σ_I/I)²

The (σ_t/t) term is negligible and dropped (live time is typically
known to better than 0.01% by the MCA clock).

═════════════════════════════════════════════════════════════════════
Background subtraction safety policy (closes K-15)
═════════════════════════════════════════════════════════════════════

Using gross peak areas for activity calculation systematically
overestimates the result by the background contribution under the
photopeak. Therefore:

  • If a background spectrum is available (caller signals via
    `bg_available=True`) AND the spectrum was NOT bg-subtracted,
    `compute_activity` raises `BackgroundNotSubtractedError`.
  • Caller can override with `force_gross=True` (e.g. for legacy
    debugging or when the bg contribution is known negligible).

═════════════════════════════════════════════════════════════════════
Cascade summing depletion (K-17 placeholder)
═════════════════════════════════════════════════════════════════════

For nuclides emitting γ-rays in coincidence (Tl-208 583+2614,
Co-60 1173+1332, Y-88 898+1836, Eu-152 multi-line, Ba-133, Eu-154),
two photons may sum into a single MCA event and be lost from the
individual photopeaks. This **depletes** the observed peak area,
biasing the activity result LOW by 5-30% depending on geometry
(worst at small source-detector distance). The correction is
γ-cascade and geometry specific.

This module accepts a `coincidence_correction` dict — `{E_keV: factor}`
where factor > 1 multiplies the observed area to recover the true
disintegration rate. When the dict is not supplied for a known-
cascade nuclide, a `cascade_warning` is emitted in the result so
the caller can flag the activity as biased-low. Full implementation
of the coincidence correction algorithm is deferred to a future
phase (see KNOWN_AND_FIXED_ISSUES.md K-17, Lsrm §10).

═════════════════════════════════════════════════════════════════════
Decay correction
═════════════════════════════════════════════════════════════════════

When a `reference_datetime` (e.g. source certificate date) is
provided alongside the spectrum's `measurement_datetime`, the
measured activity is corrected back to the reference epoch:

    A_ref = A_meas · exp(ln2 · Δt / T_½)
    Δt    = (t_meas − t_ref) seconds

For Δt > 0 (measurement after reference), the factor is > 1 — the
source has decayed since the certificate was issued, so the
measured activity is correspondingly lower. The formula reverses
sign symmetrically for Δt < 0.

Stable nuclides (T½ = ∞) yield factor = 1 with a note. Nuclides
without a known T½ in the library skip the correction and emit a
note (this is rare for the natural-gamma catalogue but can happen
for newly-loaded external nuclides).
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional, Dict, Iterable, List

from gamma.calibration.efficiency import EfficiencyCurve
from gamma.data.nuclide_library import get_nuclide

# BUG-27 / v1.18.32+: chain-sibling I_pct lookup table for in-blob
# multiplet correction. See `_chain_sibling_I_in_window` and the BUG-27
# fix in `compute_activity` (within-nuclide sum + chain-sibling sum).
try:
    from gamma.data.chain_decomposer import ALL_CHAIN_LINES as _ALL_CHAIN_LINES
except Exception:  # pragma: no cover — defensive: chain_decomposer is core
    _ALL_CHAIN_LINES = ()


# Nuclides that emit γ-rays in coincidence; their measured photopeak
# areas are biased LOW unless a cascade-summing correction is applied.
# This set drives the `cascade_warning` field in ActivityResult.
#
# This list is curated from Gilmore §8.5 and Lsrm §10. Maintenance:
# add a nuclide here when it has at least two γ-lines emitted within
# the typical detector coincidence-resolving time (~µs) from the same
# disintegration. The library's per-record `is_cascade` flag (where
# set) is consulted first; this set provides the canonical fallback.
CASCADE_SUMMING_NUCLIDES = frozenset({
    "Co-60", "Eu-152", "Eu-154", "Y-88", "Ba-133", "Tl-208",
    "Co-58", "Co-57", "Na-22", "Sc-46", "Bi-207",
    # F-128 / v1.17.7 — Bi-212 (Th-232 chain) has 727+1620 cascade pair.
    # Adding it triggers cascade_warning when no per-line coincidence_correction
    # is provided AND enables auto-TCS via gamma.physics.cascade_summing.
    "Bi-212",
})


# Area-method-aware TCS scaling (K-18, v1.7.13)
# ----------------------------------------------
# The analytic TCS correction `C(E) = 1/(1 − Σ p·ε_T)` derived in
# `gamma.physics.cascade_summing` assumes the photopeak area was
# integrated as a pure channel-sum / Cowell ROI — i.e., the only counts
# entering the area are those that landed under the gaussian and were
# NOT pulled away by coincidence summing.
#
# However, the Lsrm SpectraLine peak-table area is the result of a
# **wide-ROI Gaussian-on-step fit**. F-31a observed that this method
# already recovers most of the area deficit on close doublets, and
# v1.7.9 cert validation observed that Co-60 5 cm activity matched
# the certificate within +0.61% on Lsrm-table areas WITHOUT any TCS
# correction. Applying the full analytic TCS on top therefore
# overshoots by ~2.3% (it double-counts the recovered counts).
#
# This dict caps that double-counting. For each known area-source
# label (see `LineMatch.peak_area_source` from F-34), it gives the
# fraction of the analytic TCS effect to retain:
#
#   c_effective = 1 + (c_analytic - 1) * scale[area_source]
#
# Defaults:
#   - Cowell, deconvolved, failed, unknown → 1.0 (full TCS)
#   - Lsrm peak table → 0.0 (no TCS; Lsrm fit already recovers it)
#
# The defaults are conservative for everything except Lsrm-table.
# Callers can pass their own scale via `compute_activity(...,
# tcs_method_scale=...)` to override (e.g., 0.3 if you trust the
# Lsrm fit only partially in a particular geometry).
DEFAULT_TCS_METHOD_SCALE: Dict[str, float] = {
    "":                  1.0,   # unknown source → full TCS, safe default
    "cowell":            1.0,
    "deconvolved":       1.0,
    "failed":            1.0,
    "lsrm_peaks_table":  0.0,   # see docstring above
}


# ════════════════════════════════════════════════════════════════════════
# BUG-27 / v1.18.32+ — In-blob multiplet effective-I correction
# ════════════════════════════════════════════════════════════════════════
#
# Context. After BUG-15 (within-nuclide dedup at shared peak_channel), the
# survivor — the highest-I_pct line — receives the FULL peak area S_net
# of the unresolved blob and the dedup'd partners are dropped. The
# formula A_i = S/(ε·I·t) then assumes S contains ONLY photons of the
# survivor's library line. This is wrong when the blob is wide enough
# (NaI 63×63 at high energy: FWHM ≈ 100-120 keV) to capture multiple
# library transitions:
#
#   • Same-nuclide lines that were dedup'd by BUG-15 deposit their
#     photons in S but contribute zero to the I in the denominator
#     → A_i inflated by factor (I_dom + Σ I_dedup) / I_dom.
#   • Chain-equilibrium sibling lines (e.g. Bi-212 1620.50 keV in the
#     Ac-228 1588.20 keV blob) deposit photons in S whether or not
#     they were matched as separate library lines. Since chain
#     daughters in secular equilibrium have ≈ same activity, their
#     in-blob photons must be credited to the SAME I-sum.
#
# Closure case. Th-232 demo `peak_channel=539`, peak_E_obs=1597.35 keV,
# FWHM=115.5 keV, S=20660. Matched same-nuclide Ac-228 lines on this
# channel: 1580.53 (I=0.6%), 1588.20 (I=3.22%, winner), 1630.63
# (I=1.51%). Chain-sibling Bi-212 1620.50 (I=1.51%) inside ±FWHM/2
# but unmatched. True in-blob I_eff = 0.6+3.22+1.51+1.51 = 6.84%.
# Using I=3.22% inflates A_i by 6.84/3.22 = 2.12× → observed 3881 Bq
# vs nuclide-mean 1802 Bq (ratio 2.15). After BUG-27 fix:
# A_i(1588) ≈ 1827 Bq (within +1.4% of nuclide-mean).
#
# Cite: Gilmore §6.5 (unresolved multiplets), Lsrm §7 (multi-line
# activity combination), chain_decomposer.ALL_CHAIN_LINES for sibling
# library, KNOWN_AND_FIXED_ISSUES.md BUG-27 entry.

# Which chain a given nuclide belongs to (used by sibling lookup).
# Mirrors CHAIN_MEMBERS below but as a reverse-index. Hard-coded
# because CHAIN_MEMBERS is defined later in the same module — to
# avoid forward-reference complexity we duplicate the canonical
# membership here. Keep these in sync.
_NUCLIDE_TO_CHAIN: Dict[str, str] = {
    # Th-232 chain (Th-228 → Ra-224 → ... → Pb-208)
    "Th-232": "Th-232", "Ac-228": "Th-232", "Tl-208": "Th-232",
    "Pb-212": "Th-232", "Bi-212": "Th-232", "Th-228": "Th-232",
    "Ra-224": "Th-232",
    # U-238 chain (above and through Ra-226 sub-chain)
    "U-238": "U-238", "Bi-214": "U-238", "Pb-214": "U-238",
    "Pb-210": "U-238", "Ra-226": "U-238", "Po-214": "U-238",
    "Th-234": "U-238", "Pa-234m": "U-238",
    # U-235 chain
    "U-235": "U-235", "Th-231": "U-235",
}

# Minimum FWHM (keV) above which a peak is considered "wide enough"
# for the chain-sibling sum to apply. Narrow peaks (e.g. low-energy
# NaI, HPGe at any energy) already have FWHM small enough that
# library-line blending is rare; the same-nuclide BUG-15 dedup
# normally won't trigger there. Threshold chosen so that K-40 1460.8
# (Marinelli NaI FWHM ~ 60-70 keV) does NOT get the sibling boost
# (single isolated line, no siblings within ±30 keV anyway) but
# wide 1500-1700 keV multiplets DO get it.
_CHAIN_SIBLING_FWHM_MIN_KEV: float = 50.0

# Area-source prefixes for which the chain-sibling sum is enabled.
# Restriction: only deconvolved-coupled (multi-component fits) peaks
# physically merge unresolved lines into one S. Lsrm peak-table and
# Cowell ROI on isolated peaks should NOT get the sibling boost
# (they're applied to single-Gaussian peaks).
_CHAIN_SIBLING_SOURCES: tuple = ("deconvolved",)

# Minimum sibling library I_pct to include in the chain-sibling sum.
# A library line with I < 0.1% contributes <1% of the dominant in-blob
# photon population (Ac-228 1588 case: Bi-212 1620.50 at I=1.51%, Tl-208
# 1647.40 at I=0.077% — the latter is below detector resolution noise
# floor for any realistic S_blob with √S/S ≥ 0.5%). Excluding sub-0.1%
# lines also makes the effective-I deterministic against future library
# additions of trace emissions.
_CHAIN_SIBLING_I_MIN_PCT: float = 0.1


def _chain_sibling_I_in_window(
    *,
    nuclide: str,
    peak_E_keV: float,
    fwhm_keV: float,
    exclude_E_keV: Iterable[float] = (),
    tolerance_keV: float = 0.5,
    I_min_pct: float = _CHAIN_SIBLING_I_MIN_PCT,
) -> float:
    """BUG-27: sum of chain-sibling library I_pct within ±FWHM/2 of
    the observed peak centroid.

    Used to correct the effective emission probability of a BUG-15
    dedup survivor whose physical peak area S contains photons from
    chain-equilibrium siblings (e.g. Ac-228 1588.20 keV peak captures
    Bi-212 1620.50 keV photons at FWHM=115 keV on NaI 63×63).

    Args:
        nuclide: host nuclide name (e.g. "Ac-228"). Used to identify
            chain membership and to EXCLUDE same-nuclide lines (those
            are handled by the within-nuclide BUG-15 dedup sum, not
            here).
        peak_E_keV: observed peak centroid (keV). Window is centered
            on this.
        fwhm_keV: FWHM of the observed peak (keV). Window half-width.
        exclude_E_keV: library energies to skip (typically the
            survivor's own E and its dedup partners on the same
            channel — those are already in the same-nuclide sum).
        tolerance_keV: small ε for matching exclude_E_keV (avoid
            float equality issues).
        I_min_pct: minimum library I_pct to include in the sum.
            Defaults to _CHAIN_SIBLING_I_MIN_PCT (0.1%). Filters out
            trace siblings whose physical photon contribution is
            below detector counting-statistics noise.

    Returns:
        Sum of I_pct of chain-sibling library lines inside the window
        (0.0 if no siblings, nuclide not in known chain, FWHM < 0,
        or chain decomposer not loaded).
    """
    chain = _NUCLIDE_TO_CHAIN.get(nuclide)
    if chain is None:
        return 0.0
    if fwhm_keV <= 0 or peak_E_keV <= 0:
        return 0.0
    if not _ALL_CHAIN_LINES:
        return 0.0
    half = fwhm_keV / 2.0
    E_lo = peak_E_keV - half
    E_hi = peak_E_keV + half
    exclude_set = {round(float(e), 2) for e in exclude_E_keV}
    total = 0.0
    for E_ref, owner, I_pct, _tol in _ALL_CHAIN_LINES:
        if _NUCLIDE_TO_CHAIN.get(owner) != chain:
            continue
        if owner == nuclide:
            continue  # same-nuclide handled separately
        if not (E_lo <= E_ref <= E_hi):
            continue
        if round(float(E_ref), 2) in exclude_set:
            continue
        if float(I_pct) < I_min_pct:
            continue  # trace sibling — below counting-stat noise floor
        total += float(I_pct)
    return total


class BackgroundNotSubtractedError(RuntimeError):
    """Raised when activity is requested on a gross spectrum but a
    background is available for subtraction (caller must subtract
    first, or override with `force_gross=True`)."""

    def __init__(self, message: str, *, nuclide: str = ""):
        super().__init__(message)
        self.nuclide = nuclide


@dataclass(frozen=True)
class LineActivity:
    """Activity computed from one matched line of a nuclide."""

    E_keV: float
    I_pct: float                       # library intensity (%, not decimal)
    sigma_I_pct_relative: float        # relative uncertainty of I (% of I)
    S_net: float                       # net counts
    sigma_S: float                     # net counts uncertainty
    epsilon: float                     # photopeak efficiency at E
    epsilon_unc_pct: float             # relative ε uncertainty (% of ε)
    correction_factor: float           # multiplicative correction (≥ 1 usually)
    A_Bq: float                        # per-line activity
    sigma_A_Bq: float                  # per-line uncertainty
    weight: float                      # 1/σ² used in weighted average
    epsilon_extrapolated: bool         # True if outside ε(E) calibrated range

    def __repr__(self) -> str:
        return (f"LineActivity(E={self.E_keV:.2f} keV, "
                f"S={self.S_net:.0f}±{self.sigma_S:.0f}, "
                f"ε={self.epsilon:.3e}, "
                f"A={self.A_Bq:.3e}±{self.sigma_A_Bq:.2e} Bq)")


@dataclass(frozen=True)
class ActivityResult:
    """Aggregate activity result for one nuclide."""

    nuclide: str
    A_Bq: float                        # weighted average (NaN if no lines)
    sigma_A_Bq: float                  # uncertainty (NaN if no lines)
    lines_used: tuple                  # tuple[LineActivity, ...]
    lines_skipped: tuple = ()          # tuple[tuple[E_keV, reason], ...]

    # Quality of multi-line agreement: χ² / (n−1) of A_i around weighted mean
    intra_chi2_per_dof: Optional[float] = None

    # F-91: which σ estimate was chosen — "weighted_mean" or "scatter"
    # (LSRM §7: σ = max(scatter, weighted-mean) for multi-line nuclides)
    sigma_method: str = "weighted_mean"

    # Provenance and corrections applied
    from_bg_subtracted: bool = False
    force_gross_override: bool = False
    coincidence_correction_applied: bool = False
    cascade_warning: Optional[str] = None

    # Decay correction
    decay_corrected: bool = False
    decay_factor: float = 1.0          # multiplier A_meas → A_ref
    reference_datetime: Optional[datetime] = None
    measurement_datetime: Optional[datetime] = None
    T_half_s: Optional[float] = None

    notes: str = ""

    def n_lines_used(self) -> int:
        return len(self.lines_used)

    def is_valid(self) -> bool:
        """True if a finite activity was computed."""
        return (not math.isnan(self.A_Bq)) and self.A_Bq > 0

    def __repr__(self) -> str:
        if not self.is_valid():
            return f"ActivityResult({self.nuclide}: no activity — {self.notes})"
        warn = f" [{self.cascade_warning}]" if self.cascade_warning else ""
        decay_str = (f", decay×{self.decay_factor:.3f}"
                     if self.decay_corrected else "")
        chi2_str = (f", χ²/dof={self.intra_chi2_per_dof:.2f}"
                    if self.intra_chi2_per_dof is not None else "")
        return (f"ActivityResult({self.nuclide}: "
                f"A={self.A_Bq:.3e}±{self.sigma_A_Bq:.2e} Bq, "
                f"{self.n_lines_used()} lines"
                f"{chi2_str}{decay_str}{warn})")


def _is_cascade_nuclide(nuclide: str) -> bool:
    """Consult library record first, fall back to canonical set."""
    rec = get_nuclide(nuclide)
    if rec and rec.get("is_cascade"):
        return True
    return nuclide in CASCADE_SUMMING_NUCLIDES


def _decay_factor(
    T_half_s: Optional[float],
    measurement_datetime: Optional[datetime],
    reference_datetime: Optional[datetime],
) -> tuple:
    """
    Compute the multiplicative decay correction factor.

    Returns (factor, applied: bool, note: str).
    """
    if T_half_s is None:
        return (1.0, False, "T½ unknown — decay correction skipped")
    if measurement_datetime is None or reference_datetime is None:
        return (1.0, False, "missing datetime — decay correction skipped")
    if math.isinf(T_half_s):
        return (1.0, True, "stable nuclide (T½=∞)")
    if T_half_s <= 0:
        return (1.0, False, f"invalid T½={T_half_s}")
    dt_s = (measurement_datetime - reference_datetime).total_seconds()
    factor = math.exp(math.log(2.0) * dt_s / T_half_s)
    note = (f"Δt={dt_s:.0f}s (={dt_s/86400:.2f}d), "
            f"T½={T_half_s:.3g}s, factor={factor:.4f}")
    return (factor, True, note)


# ════════════════════════════════════════════════════════════════════════
# PTB-2018 Annex E — application-dependent T½ for Pb-214/Bi-214/Pb-212
# ════════════════════════════════════════════════════════════════════════
#
# γ-SPEKT/GRUNDL (ISSN 1865-8725, March 2018).
#
# Annex E.1 / Tab. E1 (page -67/-68) — Pb-214, Bi-214:
#   mode "equilibrium" (Tab. E1 row 1, DEFAULT): sealed environmental
#       sample (soil, sediment, construction material) measured >20 d
#       after container filling — Pb-214/Bi-214 supported by Ra-226,
#       T½(Ra-226) = 1600 a.
#   mode "rn222" (row 2): activated-carbon filter with sampled Rn-222 —
#       progenies supported by Rn-222, T½ = 3.8235 d.
#   mode "progeny" (row 3): air filter with unsupported radon progenies —
#       own half-life (Pb-214: 26.8 min, Bi-214: 19.8 min).
#
# Annex E.2 (page -68/-69) — Pb-212:
#   mode "equilibrium" (DEFAULT): solid samples (soil, sediments) OR
#       aged aqueous solution (>4 d after sampling) — Ra-224/Th-228
#       radioactive equilibrium re-established, activity determined via
#       Pb-212 with T½(Th-228) = 1.91 a.
#   mode "ra224_fresh": aqueous solution measured soon after sampling to
#       preserve short-lived Ra-224 (T½ = 3.66 d) that would otherwise
#       decay — Pb-212 tracks Ra-224, T½ := 3.66 d.
#   mode "progeny": unsupported Pb-212 (rare — separated source, decay
#       verification) — own T½ = 10.64 h.
#
# Cross-annex mode compatibility: "rn222" is E.1-only (applies only to
# Pb-214/Bi-214); "ra224_fresh" is E.2-only (applies only to Pb-212).
# For a non-matching (nuclide, mode) pair the library T½ passes through
# unchanged with note "" — silent no-op, backward-compatible.
#
# Without this remap, enabling decay correction on an equilibrium soil
# sample would apply the daughter's own short T½ (~20 min for Pb-214,
# ~10.6 h for Pb-212) over a days-scale Δt and blow the factor up by
# many orders of magnitude.
_PTB_E1_NUCLIDES = frozenset({"Pb-214", "Bi-214"})
_PTB_E2_NUCLIDES = frozenset({"Pb-212"})
_RA226_T_HALF_S = 5.0492e10   # 1600 a   (E.1 Tab. E1 row 1)
_RN222_T_HALF_S = 3.3035e5    # 3.8235 d (E.1 Tab. E1 row 2)
_TH228_T_HALF_S = 6.0275e7    # 1.91 a   (E.2 aged aqueous / soil)
_RA224_T_HALF_S = 3.1622e5    # 3.66 d   (E.2 fresh aqueous carve-out)

_CHAIN_DECAY_MODES = ("equilibrium", "rn222", "ra224_fresh", "progeny")


def _ptb_annex_e_half_life(
    nuclide: str,
    T_half_own: Optional[float],
    chain_decay_mode: str,
) -> tuple:
    """
    Return (T_half_s, note) for the decay correction of `nuclide`
    according to PTB-2018 Annex E.1 (Pb-214/Bi-214) / E.2 (Pb-212).

    Silent pass-through (unchanged T½, note "") for:
      - any nuclide outside {Pb-214, Bi-214, Pb-212};
      - mode "progeny" on any of the three;
      - E.1-only mode "rn222" on Pb-212;
      - E.2-only mode "ra224_fresh" on Pb-214/Bi-214.
    """
    if chain_decay_mode not in _CHAIN_DECAY_MODES:
        raise ValueError(
            f"unknown chain_decay_mode {chain_decay_mode!r}; "
            f"expected one of {_CHAIN_DECAY_MODES}"
        )
    if chain_decay_mode == "progeny":
        return (T_half_own, "")
    if nuclide in _PTB_E1_NUCLIDES:
        if chain_decay_mode == "equilibrium":
            return (
                _RA226_T_HALF_S,
                f"PTB E.1 equilibrium mode: {nuclide} supported by Ra-226 — "
                f"T½ := 1600 a (Tab. E1 row 1)",
            )
        if chain_decay_mode == "rn222":
            return (
                _RN222_T_HALF_S,
                f"PTB E.1 rn222 mode: {nuclide} supported by Rn-222 — "
                f"T½ := 3.8235 d (Tab. E1 row 2)",
            )
        # "ra224_fresh" on Pb-214/Bi-214 — E.2 mode not applicable, no-op.
        return (T_half_own, "")
    if nuclide in _PTB_E2_NUCLIDES:
        if chain_decay_mode == "equilibrium":
            return (
                _TH228_T_HALF_S,
                f"PTB E.2 equilibrium mode: {nuclide} supported by Th-228 — "
                f"T½ := 1.91 a (soil / aged aqueous >4 d, "
                f"Ra-224–Th-228 equilibrium)",
            )
        if chain_decay_mode == "ra224_fresh":
            return (
                _RA224_T_HALF_S,
                f"PTB E.2 ra224_fresh mode: {nuclide} tracks Ra-224 — "
                f"T½ := 3.66 d (aqueous solution measured soon after "
                f"sampling)",
            )
        # "rn222" on Pb-212 — E.1 mode not applicable, no-op.
        return (T_half_own, "")
    return (T_half_own, "")


# Back-compat alias for external callers (F-rule PTB-5 landed with the
# E.1-only name; keep pointing to the generalized function).
_ptb_e1_half_life = _ptb_annex_e_half_life


def compute_activity(
    nuclide_id,
    *,
    efficiency_curve: EfficiencyCurve,
    live_time_s: float,
    from_bg_subtracted: bool,
    bg_available: bool = False,
    force_gross: bool = False,
    coincidence_correction: Optional[Dict[float, float]] = None,
    tcs_method_scale: Optional[Dict[str, float]] = None,
    decay_correction: bool = True,
    reference_datetime: Optional[datetime] = None,
    measurement_datetime: Optional[datetime] = None,
    epsilon_unc_pct: float = 5.0,
    min_intensity_pct: float = 0.0,
    self_attenuation_factors: Optional[Dict[float, float]] = None,
    enable_tcs_correction: bool = False,
    tcs_detector_id: str = "Gamma-1S",
    cross_nuclide_skip_energies_keV: Optional[set] = None,
    chain_decay_mode: str = "equilibrium",
) -> ActivityResult:
    """
    Compute the activity of one identified nuclide from its matched
    γ-lines.

    Args:
        nuclide_id: a `NuclideIdentification` (from
            `gamma.identification.identify`). Must be `detected=True`
            and have at least one `LineMatch` with a positive
            `peak_area`.
        efficiency_curve: fitted `EfficiencyCurve` for the detector
            and geometry of the spectrum.
        live_time_s: live time of the source spectrum (s). Must be > 0.
        from_bg_subtracted: True if the spectrum on which the matched
            peak areas were integrated has been background-subtracted.
            False on a gross spectrum.
        bg_available: True if a background spectrum is available for
            this analysis (i.e. it should have been subtracted before
            calling). If True AND `from_bg_subtracted=False` AND
            `force_gross=False`, the call raises
            `BackgroundNotSubtractedError`. Default False (no bg
            available → safe to compute on gross).
        force_gross: bypass the bg-subtraction safety check. For
            legacy/debug paths only — adds a note to the result.
        coincidence_correction: optional dict mapping line energy
            (keV) → multiplicative correction factor (typically > 1
            for cascade-affected lines). If a known-cascade nuclide
            is identified WITHOUT this argument, a cascade_warning
            is emitted (the result is still computed but flagged as
            biased low).
        tcs_method_scale: optional dict mapping `LineMatch.peak_area_source`
            label → fraction of the analytic TCS effect to retain
            (K-18). When None, `DEFAULT_TCS_METHOD_SCALE` applies:
            the Lsrm peak-table source gets scale=0 (the wide-ROI
            Gaussian fit already recovers the summing-displaced
            counts that the analytic correction would otherwise add),
            every other source gets scale=1.0. The effective
            per-line correction becomes
            `c_eff = 1 + (c_analytic − 1) · scale[area_source]`.
            Unknown source labels fall back to scale=1.0 (full
            correction — safe).
        decay_correction: enable decay correction if reference and
            measurement datetimes both supplied.
        reference_datetime: certificate / reference date for the
            source. If supplied with measurement_datetime, the
            reported activity is corrected to this epoch.
        measurement_datetime: when the spectrum was acquired. If
            None, decay correction is skipped.
        epsilon_unc_pct: relative uncertainty of the efficiency
            curve (% of ε) used in error propagation. Default 5.0%
            (typical for log-log polynomial fit of NaI ε at the
            reference energy).
        min_intensity_pct: skip library lines below this absolute
            intensity (default 0 = use all matched lines).

    Returns:
        ActivityResult.

    Raises:
        ValueError: invalid live_time, no efficiency curve, etc.
        BackgroundNotSubtractedError: if bg_available and
            !from_bg_subtracted and !force_gross.
    """
    if live_time_s <= 0:
        raise ValueError(f"live_time_s must be > 0, got {live_time_s}")
    if efficiency_curve is None:
        raise ValueError("efficiency_curve is required")

    nuclide = getattr(nuclide_id, "nuclide", "?")

    # ── Safety: refuse gross calculation when bg is available ─────────
    if bg_available and not from_bg_subtracted and not force_gross:
        raise BackgroundNotSubtractedError(
            f"Cannot compute activity for {nuclide!r}: spectrum is "
            f"NOT background-subtracted but a background spectrum is "
            f"available. Subtract the background first, or pass "
            f"force_gross=True to override.",
            nuclide=nuclide,
        )

    if not getattr(nuclide_id, "detected", False):
        return ActivityResult(
            nuclide=nuclide,
            A_Bq=float("nan"),
            sigma_A_Bq=float("nan"),
            lines_used=(),
            lines_skipped=(),
            from_bg_subtracted=from_bg_subtracted,
            force_gross_override=force_gross,
            notes="Nuclide not detected — no activity computed",
        )

    matched = getattr(nuclide_id, "matched_lines", ())
    coincidence_correction = coincidence_correction or {}

    # F-296 / v1.18.1 — opt-in auto-TCS correction для cascade нуклидов
    # с preset cascades (Co-60, Eu-152, Ba-133). Активируется ТОЛЬКО когда:
    #   1) enable_tcs_correction=True
    #   2) coincidence_correction dict пустой (не дублируем user-supplied)
    #   3) nuclide есть в CASCADE_PRESETS
    # Использует F-295 total_efficiency_from_fep для ε_T(E_j).
    # Результат: per-line auto-correction добавляется в coincidence_correction
    # (тогда дальше работает существующая K-18 TCS-pipeline без изменений).
    tcs_auto_applied: Dict[float, float] = {}
    if enable_tcs_correction and not coincidence_correction:
        try:
            from gamma.activity.tcs_close_geometry import (
                CASCADE_PRESETS, compute_tcs_correction,
            )
            from gamma.activity.pt_ratio_nai import pt_ratio_for_detector
            pairs = CASCADE_PRESETS.get(nuclide, [])
            if pairs:
                def _eps_T(E_j: float) -> float:
                    eps_j = efficiency_curve.efficiency_at(E_j)
                    if eps_j is None or eps_j <= 0:
                        return 0.0
                    try:
                        return eps_j / pt_ratio_for_detector(
                            E_j, tcs_detector_id,
                        )
                    except Exception as exc:  # DEEP-06
                        warnings.warn(
                            f"[F296] TCS total-efficiency lookup failed at "
                            f"E={E_j:.2f} keV for detector_id={tcs_detector_id!r} "
                            f"({exc!r}). Unknown detector — cannot compute TCS; "
                            f"ε_T set to 0.0, C_TCS=1.0 (effective TCS-skip). "
                            f"To register this detector add an entry in "
                            f"'scripts/gamma/activity/pt_ratio_nai.py' "
                            f"(PT_RATIO_TABLE dict).",
                            UserWarning,
                            stacklevel=2,
                        )
                        return 0.0
                for m in matched:
                    E_i = float(getattr(m, "library_E_keV", 0.0))
                    if E_i <= 0:
                        continue
                    try:
                        res = compute_tcs_correction(
                            E_i, pairs, _eps_T,
                        )
                    except Exception as exc:  # DEEP-06
                        warnings.warn(
                            f"TCS per-line correction failed for "
                            f"nuclide={nuclide!r} E={E_i:.2f} keV "
                            f"({exc!r}); skipping this line — activity for "
                            f"this gamma may be uncorrected for true-coincidence.",
                            stacklevel=2,
                        )
                        continue
                    if res.is_significant:
                        tcs_auto_applied[E_i] = res.correction_factor
        except Exception as exc:  # DEEP-06
            warnings.warn(
                f"TCS auto-correction block failed entirely for "
                f"nuclide={nuclide!r} ({exc!r}); proceeding without "
                f"auto-TCS — delivered activities will not carry "
                f"true-coincidence-summing correction.",
                stacklevel=2,
            )
        if tcs_auto_applied:
            coincidence_correction = dict(tcs_auto_applied)

    # K-18: area-method-aware TCS scaling. None → defaults; dict → merge
    # with defaults so callers only need to override what differs.
    if tcs_method_scale is None:
        method_scale = dict(DEFAULT_TCS_METHOD_SCALE)
    else:
        method_scale = dict(DEFAULT_TCS_METHOD_SCALE)
        method_scale.update(tcs_method_scale)
    n_tcs_scaled = 0  # count lines where scale ≠ 1 was applied

    # σ_I per line is not stored on LineMatch (it's library data, not
    # match-result data). Look it up once from the library record;
    # the library stores `lines` as [E_keV, I_pct, sigma_I_pct].
    nuc_lib_record = get_nuclide(nuclide)
    nuc_lib_lines = nuc_lib_record.get("lines", []) if nuc_lib_record else []

    def _sigma_I_for(E: float) -> float:
        """Look up σ_I (in % units, same as I_pct) for a given energy.
        Returns 0.0 if the library record doesn't carry σ_I."""
        for ll in nuc_lib_lines:
            if abs(float(ll[0]) - E) < 0.05:
                return float(ll[2]) if len(ll) > 2 else 0.0
        return 0.0

    line_activities = []
    skipped = []

    # ─── BUG-15 / v1.18.31 — shared-peak deduplication ──────────────────
    # When multiple library lines of the SAME nuclide fall within the
    # same observed peak (typical for NaI's poor resolution: e.g. Ac-228
    # 911 + 904 keV unresolved, or 1580 + 1588 + 1630 keV unresolved),
    # the upstream coupled-multiplet fitter may assign the SAME peak
    # area to all components when the fit is degenerate (intensity ratios
    # not independently constrained). Crediting each library line with
    # the full peak area means the low-intensity lines yield A_i values
    # inflated by factors of 5-50× (because A_i ∝ 1/I_pct).
    #
    # Physical interpretation: when N library lines are unresolved into
    # one observed peak with one well-determined area S, the only
    # unbiased per-line estimator is to assign S to the line with the
    # highest I_pct (the dominant transition). The other components are
    # absorbed into its activity estimate. Reporting independent A_i
    # for each would be double-counting.
    #
    # This block builds a set of (peak_channel) keys that are safe to
    # use. For groups with >1 matched line sharing peak_channel and
    # S_net > 0, only the highest-I_pct line is kept; the rest are
    # added to `skipped` with reason "shared_peak_dedupe".
    #
    # See Lsrm §7 (combining multi-line activities) and Gilmore §6.5
    # (unresolved multiplet handling).
    _by_channel: Dict[int, list] = {}
    for _m in matched:
        _ch = int(getattr(_m, "peak_channel", -1) or -1)
        _S = getattr(_m, "peak_area", None)
        if _ch < 0 or _S is None or float(_S) <= 0:
            continue
        _by_channel.setdefault(_ch, []).append(_m)
    # For each shared-peak group, determine the survivor (highest I_pct).
    # Tie-break: highest peak_area, then smallest E.
    _shared_skip_keys: set = set()  # ids of matched_lines to skip
    # BUG-27 / v1.18.32+: in-blob I_eff boost map. For each survivor
    # of within-nuclide dedup, store the SUM of I_pct of the entire
    # dedup group (winner + all skipped same-nuclide siblings sharing
    # peak_channel). Used in A_i computation below to replace the
    # winner's bare library_I_pct with the physically correct
    # in-blob effective intensity. See module-level BUG-27 docstring
    # and chain_sibling_I_in_window for the cross-nuclide companion.
    _within_nuclide_I_eff: Dict[int, float] = {}
    # Also remember per-winner the set of dedup partner E_keVs so
    # the chain-sibling lookup can exclude them (same-nuclide handled
    # here, sibling sum handled at A_i time).
    _within_nuclide_partner_E: Dict[int, list] = {}
    for _ch, _group in _by_channel.items():
        if len(_group) <= 1:
            continue
        _winner = max(
            _group,
            key=lambda x: (
                float(getattr(x, "library_I_pct", 0.0) or 0.0),
                float(getattr(x, "peak_area", 0.0) or 0.0),
                -float(getattr(x, "library_E_keV", 0.0) or 0.0),
            ),
        )
        _I_sum = 0.0
        _partner_E = []
        for _m in _group:
            _I_sum += float(getattr(_m, "library_I_pct", 0.0) or 0.0)
            _partner_E.append(float(getattr(_m, "library_E_keV", 0.0) or 0.0))
            if _m is not _winner:
                _shared_skip_keys.add(id(_m))
        _within_nuclide_I_eff[id(_winner)] = _I_sum
        _within_nuclide_partner_E[id(_winner)] = _partner_E

    for m in matched:
        E = float(getattr(m, "library_E_keV", 0.0))
        I_pct = float(getattr(m, "library_I_pct", 0.0))
        S_net = getattr(m, "peak_area", None)
        sigma_S = getattr(m, "peak_area_uncertainty", None)
        # σ_I from library (in % units, e.g. 0.17 for 661.66 keV Cs-137)
        sigma_I_pct_rel = _sigma_I_for(E)

        # BUG-15: skip duplicate library lines sharing one observed peak
        # (only keep the dominant-I_pct line, see block above).
        if id(m) in _shared_skip_keys:
            ch = int(getattr(m, "peak_channel", -1) or -1)
            skipped.append(
                (E, f"shared_peak_dedupe ch={ch} (I={I_pct:.3f}% not dominant)")
            )
            continue

        # BUG-15: skip cross-nuclide false-match lines (peak owned by
        # the characteristic line of a different nuclide). See
        # compute_activities_for_all for the cross-nuclide ownership
        # rule that builds this skip set.
        if cross_nuclide_skip_energies_keV is not None:
            if round(E, 2) in cross_nuclide_skip_energies_keV:
                ch = int(getattr(m, "peak_channel", -1) or -1)
                skipped.append(
                    (E,
                     f"cross_nuclide_peak_owned ch={ch} (peak owned by "
                     f"another nuclide's characteristic line)")
                )
                continue

        if S_net is None or S_net <= 0:
            skipped.append((E, "no peak area"))
            continue
        if I_pct <= 0 or I_pct < min_intensity_pct:
            skipped.append((E, f"I={I_pct:.3f}% below threshold"))
            continue

        eps = efficiency_curve.efficiency_at(E)
        if eps is None or eps <= 0:
            skipped.append((E, f"ε({E:.1f})={eps}"))
            continue

        extrapolated = efficiency_curve.is_extrapolating(E)

        # Coincidence correction factor (default 1.0 = no correction)
        c = float(coincidence_correction.get(E, 1.0))
        # Allow approximate-key match within 0.5 keV. Library energies
        # are stored with sub-keV precision (e.g. 1173.23), but user
        # dicts often use rounded values (e.g. 1173). The 0.5 keV
        # tolerance is safe: cascade-relevant photopeaks from the same
        # nuclide are separated by ≥10 keV (typically 100+ keV), so
        # the wrong line cannot be picked up.
        if c == 1.0 and coincidence_correction:
            for key_E, key_c in coincidence_correction.items():
                if abs(float(key_E) - E) < 0.5:
                    c = float(key_c)
                    break

        # K-18: scale the TCS effect by area-method factor. Lsrm-table
        # areas already partially recover summing-displaced counts via
        # their wide-ROI Gaussian fit, so applying the full analytic
        # correction on top would double-count. See
        # `DEFAULT_TCS_METHOD_SCALE`.
        if c != 1.0:
            area_source = str(getattr(m, "peak_area_source", "") or "")
            alpha = method_scale.get(area_source, 1.0)
            if alpha != 1.0:
                c = 1.0 + (c - 1.0) * alpha
                n_tcs_scaled += 1

        # F-122 / v1.17.6 — self-attenuation per-line correction (для
        # Marinelli / Дента / Петри геометрий). self_attenuation_factors
        # — словарь {E_keV: F_ref/F_sample}. Применяется как
        # умножение в формуле активности (фактически — корректирует ε).
        # При factor=1.0 (по умолчанию) — не влияет на результат.
        if self_attenuation_factors:
            sa_factor = 1.0
            for key_E, key_f in self_attenuation_factors.items():
                if abs(float(key_E) - E) < 0.5:
                    sa_factor = float(key_f)
                    break
            if sa_factor > 0 and sa_factor != 1.0:
                c = c * sa_factor

        # ── BUG-27 / v1.18.32+ — effective-I in-blob correction ────────
        # When this line is the survivor of a BUG-15 within-nuclide
        # dedup group, the peak area S_net physically contains photons
        # from ALL dedup'd library lines (winner + skipped partners).
        # Using only the winner's library_I_pct in the denominator
        # inflates A_i by factor (sum I_dedup) / I_winner. Replace
        # I_pct with the in-blob sum stored in _within_nuclide_I_eff.
        # Additionally, for wide deconvolved-coupled multiplet peaks,
        # add chain-equilibrium-sibling library lines that fall in
        # ±FWHM/2 of the observed centroid (e.g. Ac-228 1588 blob
        # captures Bi-212 1620.50 photons in Th-232 secular eq.).
        # The peak's physical FWHM is approximated by 2.355·peak_sigma
        # (peak_sigma stored on LineMatch is the Gaussian σ in keV).
        I_eff_pct = I_pct  # default: no correction
        bug27_note_parts = []
        winner_I_sum = _within_nuclide_I_eff.get(id(m))
        is_dedup_winner = winner_I_sum is not None
        if is_dedup_winner and winner_I_sum > I_pct:
            I_eff_pct = winner_I_sum
            bug27_note_parts.append(
                f"within-nuclide ΣI={winner_I_sum:.3f}%"
            )
        # Chain-sibling sum (BUG-27 cross-nuclide blend contribution).
        # Gate: only apply to peaks where within-nuclide dedup actually
        # fired (multiple library lines collided on the same channel).
        # This is the only positive evidence the peak is an unresolved
        # multiplet blob — without dedup, a "single line" matched at
        # the centroid is presumed isolated even if its FWHM happens
        # to be large (low-resolution detector at low energy).
        peak_area_source = str(
            getattr(m, "peak_area_source", "") or ""
        )
        if is_dedup_winner and any(peak_area_source.startswith(p)
               for p in _CHAIN_SIBLING_SOURCES):
            peak_E_obs = float(getattr(m, "peak_E_keV", 0.0) or 0.0)
            # BUG-34 Phase 3c R1: prefer Semantic-B gauss_sigma_keV
            # (W3 writer), fall back to legacy peak_sigma for callers/
            # fixtures that have not migrated yet. Numerically identical
            # when only peak_sigma is set. See PLAN_v1_18_32_to_v1_19_0
            # §3 Phase 3.
            _sigma = getattr(m, "gauss_sigma_keV", None)
            if _sigma is None:
                _sigma = float(getattr(m, "peak_sigma", 0.0) or 0.0) or None
            # Convert σ → FWHM. gauss_sigma_keV / peak_sigma on
            # LineMatch is the Gaussian σ in keV (set during
            # identification phase). If both are missing/zero, derive
            # from a stored peak_fwhm if available; otherwise skip
            # the sibling boost (safe default).
            fwhm_obs = (2.354820045 * _sigma
                        if _sigma is not None and _sigma > 0
                        else float(getattr(m, "fwhm_keV", 0.0) or 0.0))
            if (fwhm_obs >= _CHAIN_SIBLING_FWHM_MIN_KEV
                    and peak_E_obs > 0):
                # Exclude same-nuclide partner E_keVs (already in
                # within-nuclide ΣI) and the winner's own E.
                exclude = list(
                    _within_nuclide_partner_E.get(id(m), [E])
                )
                if not exclude:
                    exclude = [E]
                sibling_I = _chain_sibling_I_in_window(
                    nuclide=nuclide,
                    peak_E_keV=peak_E_obs,
                    fwhm_keV=fwhm_obs,
                    exclude_E_keV=exclude,
                )
                if sibling_I > 0:
                    I_eff_pct = I_eff_pct + sibling_I
                    bug27_note_parts.append(
                        f"chain-siblings ΣI+={sibling_I:.3f}% "
                        f"(±FWHM/2={fwhm_obs/2:.1f} keV)"
                    )

        # Per-line activity:
        #    A = S · c / (ε · I_decimal · t)
        I_decimal = I_eff_pct / 100.0
        A_i = (S_net * c) / (eps * I_decimal * live_time_s)

        # Uncertainty propagation:
        #    (σ_A / A)² = (σ_S/S)² + (σ_ε/ε)² + (σ_I/I)²
        if sigma_S is None or sigma_S <= 0:
            # Use Poisson approximation as a floor
            sigma_S_eff = math.sqrt(max(1.0, S_net))
        else:
            sigma_S_eff = float(sigma_S)
        rel_S = sigma_S_eff / S_net
        rel_eps = epsilon_unc_pct / 100.0
        # sigma_I_pct_rel is the absolute σ on I in percent units
        # (e.g. for I=85.1%, sigma_I_pct_rel=0.17% means I = 85.1±0.17 %).
        # Relative uncertainty σ_I / I = sigma_I_pct_rel / I_pct_eff
        # (use the effective I in the denominator — the σ_I from the
        # library is on the dominant component; the in-blob ΣI carries
        # roughly the same relative library uncertainty).
        if I_eff_pct > 0:
            rel_I = sigma_I_pct_rel / I_eff_pct
        else:
            rel_I = 0.0
        rel_A_sq = rel_S * rel_S + rel_eps * rel_eps + rel_I * rel_I
        sigma_A_i = A_i * math.sqrt(rel_A_sq)
        if sigma_A_i <= 0:
            skipped.append((E, "zero uncertainty"))
            continue
        w_i = 1.0 / (sigma_A_i * sigma_A_i)

        # Record BUG-27 audit trail when correction fired (visible in
        # lines_skipped as "BUG27_audit" entries — not actually skipped,
        # but the only side-channel available without breaking the
        # frozen LineActivity contract).
        if bug27_note_parts:
            skipped.append(
                (E, "BUG27_audit: " + ", ".join(bug27_note_parts))
            )

        line_activities.append(LineActivity(
            E_keV=E, I_pct=I_eff_pct,
            sigma_I_pct_relative=sigma_I_pct_rel,
            S_net=float(S_net), sigma_S=sigma_S_eff,
            epsilon=eps, epsilon_unc_pct=epsilon_unc_pct,
            correction_factor=c,
            A_Bq=A_i, sigma_A_Bq=sigma_A_i, weight=w_i,
            epsilon_extrapolated=extrapolated,
        ))

    # ── Aggregate ────────────────────────────────────────────────────
    if not line_activities:
        return ActivityResult(
            nuclide=nuclide,
            A_Bq=float("nan"),
            sigma_A_Bq=float("nan"),
            lines_used=(),
            lines_skipped=tuple(skipped),
            from_bg_subtracted=from_bg_subtracted,
            force_gross_override=force_gross,
            notes=(f"No usable lines for activity "
                   f"({len(skipped)} skipped)"),
        )

    # F-452-FU2 (2026-06-22) — Currie L_C pre-filter ДО BUG-38/39 MAD-rejection.
    # Линии с A_i / σ_A_i < 1.0 — статистически НЕ задетектированы (Currie 1968,
    # Gilmore §5.5 detection limit). На Th-232 Marinelli fixture F-452 poly-4
    # FWHM приводит к более «честному» deconvolved_coupled fit, где
    # второстепенные линии Tl-208 (510, 860, 277 keV) получают peak_area ~ 1e-21
    # (numerical noise вокруг нуля), но БОЛЬШИНСТВО матченных линий нуклида
    # именно такие — тогда MAD-median падает в зону наножёлтых значений, и
    # РЕАЛЬНЫЕ линии (2614, 583 keV → A_i ~ 2700 Bq) выбрасываются как 19σ
    # outliers против консенсуса нанолиний → A_avg = 9.81e-24 Bq (наблюдалось
    # 2026-06-22, F-452 root cause investigation).
    #
    # Fix: убрать ДО aggregate любые линии с A_i ≤ 0 или ниже Currie L_C — они
    # не несут информации об активности, лишь искажают MAD-median. Сохраняем
    # provenance в lines_skipped с пометкой `below_Currie_LC`. Безопасно: если
    # все линии non-detection, line_activities остаётся пустым и обработается
    # стандартной веткой `if not line_activities` ниже.
    _PRE_MAD_CURRIE_LC = 1.0
    _significant: List[LineActivity] = []
    _below_LC: List[LineActivity] = []
    for la in line_activities:
        if la.A_Bq <= 0 or la.sigma_A_Bq <= 0:
            _below_LC.append(la)
            continue
        if la.A_Bq / la.sigma_A_Bq < _PRE_MAD_CURRIE_LC:
            _below_LC.append(la)
            continue
        _significant.append(la)
    if _below_LC and _significant:
        for la in _below_LC:
            ratio = (la.A_Bq / la.sigma_A_Bq) if la.sigma_A_Bq > 0 else float("nan")
            skipped.append(
                (la.E_keV,
                 f"below_Currie_LC A_i={la.A_Bq:.3g} Bq σ={la.sigma_A_Bq:.3g} "
                 f"A/σ={ratio:.3g} < {_PRE_MAD_CURRIE_LC:.1f} "
                 f"(pre-MAD non-detection filter, F-452-FU2)")
            )
        line_activities = _significant
    elif _below_LC and not _significant:
        # Все линии non-detection — оставляем как есть, ниже Currie L_C считается
        # «нет полезной информации»; обработается стандартной aggregate-веткой.
        pass

    if not line_activities:
        return ActivityResult(
            nuclide=nuclide,
            A_Bq=float("nan"),
            sigma_A_Bq=float("nan"),
            lines_used=(),
            lines_skipped=tuple(skipped),
            from_bg_subtracted=from_bg_subtracted,
            force_gross_override=force_gross,
            notes=(f"No usable lines for activity after Currie L_C filter "
                   f"({len(skipped)} skipped)"),
        )

    sum_w = sum(la.weight for la in line_activities)
    A_avg = sum(la.weight * la.A_Bq for la in line_activities) / sum_w
    sigma_A_weighted = 1.0 / math.sqrt(sum_w)

    # ── BUG-38/39 / v1.22.0 — Eu-152 multi-line outlier rejection ────
    # When a library line of a multi-line nuclide is matched to a peak
    # owned by an UN-modelled physical process (e.g. Eu-152 503.47 keV
    # picking up Ti-44 511 keV β+ annihilation, or 656.5 keV catching
    # Cs-137 661.66 keV in unresolved doublet), A_i can be wrong by
    # 3-4 orders of magnitude. Weights are inverse-variance, so a
    # huge A_i carries vanishing weight in the weighted mean — but it
    # bloats the χ²/dof term and forces sigma_method="scatter" with
    # σ_A inflated by 5-10×. The result is then often flagged as an
    # upper-limit even though the consistent lines (1408 keV alone
    # gives −11% from cert on Eu-152) yield a perfectly usable result.
    #
    # Fix: BEFORE finalising the aggregate, apply Chauvenet-style 3σ
    # outlier rejection on per-line A_i around the initial weighted
    # mean. Drop lines with |A_i − Ā_init| > 3·σ_A_i (i.e. lines that
    # disagree with the consensus by more than 3 sigma of their OWN
    # uncertainty). Recompute A_avg, σ_A_weighted, χ²/dof, σ_scatter
    # on the surviving lines.
    #
    # Safety rails:
    #   • Only apply when n_lines ≥ 3 (need at least 2 surviving for
    #     a meaningful intra-nuclide consistency assessment).
    #   • If rejection would leave < 1 line, REVERT (keep all).
    #   • Limit to ONE iteration (avoid runaway rejection on chains
    #     with truly broken secular equilibrium).
    #
    # Cite: Gilmore §5.7.3 (multi-line consistency check); Bevington &
    # Robinson "Data Reduction" §3.4 (Chauvenet's criterion); analogous
    # to ISO 11929 §6.3 consistency test before reporting combined A.
    # Robust-reference outlier rejection:
    #   Step 1: compute MEDIAN A_i as the initial consensus reference
    #           (robust to a single rogue line dragging the weighted
    #           mean far off — the very pathology we are fixing).
    #   Step 2: compute MAD = median(|A_i − median|) and convert to
    #           a Gaussian-equivalent σ via σ ≈ 1.4826·MAD (the standard
    #           consistency factor for the normal distribution).
    #   Step 3: reject lines whose deviation from the median exceeds
    #           MAX(N·σ_MAD, N·σ_A_i) where N = 3. We use the larger of
    #           the two scale estimates so a line is rejected ONLY if it
    #           disagrees BOTH with the consensus spread of the other
    #           lines AND with its own statistical uncertainty.
    #
    # Why MAD-based reference: the weighted mean is by construction
    # dominated by the line with the smallest σ_A_i. In the production
    # Eu-152 case, that is the 1408 keV line (A=1801, σ=94). The 503
    # keV outlier has A_i = 2.4e6 with huge σ (124 948), so its weight
    # is negligible — the weighted mean ends up at ≈ 825 Bq (closer to
    # 1408 line, pulled down a bit). Computing z = (A_i − Ā)/σ_A_i for
    # the 503 line gives z ≈ 19 (clearly an outlier), but ALSO for the
    # 121 keV line we get z ≈ 14 because the 121 line's σ_A_i is small
    # while its A_i (627 Bq) differs from Ā (825). Using the median A_i
    # as the reference fixes this: median([627, 1429, 5.6e6, 1856]) =
    # 1642 Bq → 121 keV deviation = 1015 Bq vs σ ≈ 400 → z ≈ 2.5 (kept);
    # 503 keV deviation = 5.6e6 vs σ ≈ 2.9e5 → z ≈ 19 (rejected).
    outliers_rejected: List[LineActivity] = []
    if len(line_activities) >= 3:
        OUTLIER_SIGMA_THRESHOLD = 3.0
        sorted_A = sorted(la.A_Bq for la in line_activities)
        n = len(sorted_A)
        if n % 2 == 1:
            A_median = sorted_A[n // 2]
        else:
            A_median = 0.5 * (sorted_A[n // 2 - 1] + sorted_A[n // 2])
        # MAD: median absolute deviation from the median
        mad = sorted(abs(la.A_Bq - A_median) for la in line_activities)
        if n % 2 == 1:
            MAD = mad[n // 2]
        else:
            MAD = 0.5 * (mad[n // 2 - 1] + mad[n // 2])
        # Gaussian-equivalent σ from MAD (Rousseeuw & Croux 1993).
        sigma_MAD = 1.4826 * MAD
        survivors: List[LineActivity] = []
        for la in line_activities:
            if la.sigma_A_Bq <= 0:
                survivors.append(la)
                continue
            dev = abs(la.A_Bq - A_median)
            # Reject only if deviation exceeds 3σ of BOTH per-line σ_A
            # AND the MAD-based consensus scale (or MAD-based when
            # σ_MAD == 0 — degenerate identical-A case).
            scale = max(la.sigma_A_Bq, sigma_MAD)
            z = dev / scale if scale > 0 else 0.0
            if z > OUTLIER_SIGMA_THRESHOLD:
                outliers_rejected.append(la)
            else:
                survivors.append(la)
        # Recompute on survivors only if rejection produced a non-empty
        # subset AND removed ≥1 line. Need ≥2 survivors so we still
        # have an aggregate worth combining; if only 1 survivor, we
        # could still report it, but the caller's confidence model
        # treats single-line as a special case — keep ≥1 here, the
        # weighted-mean math handles n=1 cleanly.
        if outliers_rejected and len(survivors) >= 1:
            line_activities = survivors
            sum_w = sum(la.weight for la in line_activities)
            A_avg = (
                sum(la.weight * la.A_Bq for la in line_activities) / sum_w
            )
            sigma_A_weighted = 1.0 / math.sqrt(sum_w)
            for la in outliers_rejected:
                # Quote the deviation against the post-rejection Ā
                # (more meaningful diagnostic for the operator).
                dev = abs(la.A_Bq - A_avg)
                z_report = (
                    dev / la.sigma_A_Bq if la.sigma_A_Bq > 0 else 0.0
                )
                skipped.append(
                    (la.E_keV,
                     f"outlier_rejected_3sigma A_i={la.A_Bq:.3g} Bq "
                     f"vs Ā={A_avg:.3g} Bq (z={z_report:.1f}σ; "
                     f"MAD-ref median={A_median:.3g})")
                )
        else:
            # Either no outliers, or rejection would empty the set.
            outliers_rejected = []

    # Intra-nuclide χ²/dof: how well do the per-line A_i agree?
    # Also compute scatter-based σ per LSRM §7 (стр. 7-6):
    #   σ_scatter = sqrt( Σ w_i (A_i − Ā)² / ((n−1) Σ w_i) )
    # The reported uncertainty is the MAX of the two — this protects
    # against overconfident σ when the per-line A_i disagree (e.g.
    # broken secular-equilibrium assumption in a DPR chain).
    sigma_A_scatter = None
    if len(line_activities) > 1:
        chi2 = sum(la.weight * (la.A_Bq - A_avg) ** 2 for la in line_activities)
        intra_chi2_per_dof = chi2 / (len(line_activities) - 1)
        # Scatter-based σ — analogue of Lsrm §7 spread term
        sigma_A_scatter = math.sqrt(chi2 / ((len(line_activities) - 1) * sum_w))
    else:
        intra_chi2_per_dof = None

    # F-91: LSRM §7 — σ = max(scatter, weighted-mean)
    if sigma_A_scatter is not None and sigma_A_scatter > sigma_A_weighted:
        sigma_A = sigma_A_scatter
        sigma_method = "scatter"
    else:
        sigma_A = sigma_A_weighted
        sigma_method = "weighted_mean"

    # ── Cascade-summing flag ─────────────────────────────────────────
    cascade_warning = None
    if _is_cascade_nuclide(nuclide) and not coincidence_correction:
        cascade_warning = (
            f"{nuclide} is a known cascade emitter (K-17); activity "
            f"is biased LOW without coincidence-summing correction. "
            f"Magnitude of bias is geometry-dependent (5-30%, worst "
            f"at small source-detector distance). Provide a "
            f"coincidence_correction dict to apply per-line factors."
        )

    # ── Decay correction ─────────────────────────────────────────────
    T_half = nuc_lib_record.get("T_half_s") if nuc_lib_record else None
    # PTB-2018 Annex E: remap T½ for Pb-214/Bi-214 (E.1) / Pb-212 (E.2)
    T_half_eff, ptb_e_note = _ptb_annex_e_half_life(nuclide, T_half, chain_decay_mode)
    decay_factor = 1.0
    decay_applied = False
    decay_note = ""
    if decay_correction:
        decay_factor, decay_applied, decay_note = _decay_factor(
            T_half_eff, measurement_datetime, reference_datetime,
        )
        if ptb_e_note and decay_applied:
            decay_note = f"{ptb_e_note}; {decay_note}"
        if decay_applied:
            A_avg = A_avg * decay_factor
            sigma_A = sigma_A * decay_factor

    # ── Assemble notes ───────────────────────────────────────────────
    note_parts = []
    if force_gross and not from_bg_subtracted:
        note_parts.append("force_gross=True (bg-subtraction bypassed)")
    if any(la.epsilon_extrapolated for la in line_activities):
        n_extrap = sum(1 for la in line_activities if la.epsilon_extrapolated)
        note_parts.append(
            f"{n_extrap} line(s) used ε(E) outside calibrated range — "
            f"extrapolation"
        )
    if coincidence_correction:
        if tcs_auto_applied:
            note_parts.append(
                f"F-296 / v1.18.1: auto-TCS correction applied to "
                f"{len(tcs_auto_applied)} line(s) "
                f"(detector={tcs_detector_id}, presets from "
                f"CASCADE_PRESETS)"
            )
        else:
            note_parts.append(
                f"coincidence correction applied to "
                f"{len(coincidence_correction)} line(s)"
            )
        if n_tcs_scaled > 0:
            note_parts.append(
                f"K-18: TCS scaled by area-method on "
                f"{n_tcs_scaled} line(s)"
            )
    if decay_note:
        note_parts.append(f"decay: {decay_note}")
    if outliers_rejected:
        rejected_E = ", ".join(f"{la.E_keV:.2f}" for la in outliers_rejected)
        note_parts.append(
            f"BUG-38/39: {len(outliers_rejected)} line(s) rejected as "
            f"3σ outliers [{rejected_E}] keV"
        )
    notes = "; ".join(note_parts)

    return ActivityResult(
        nuclide=nuclide,
        A_Bq=A_avg,
        sigma_A_Bq=sigma_A,
        lines_used=tuple(line_activities),
        lines_skipped=tuple(skipped),
        intra_chi2_per_dof=intra_chi2_per_dof,
        sigma_method=sigma_method,
        from_bg_subtracted=from_bg_subtracted,
        force_gross_override=(force_gross and not from_bg_subtracted),
        coincidence_correction_applied=bool(coincidence_correction),
        cascade_warning=cascade_warning,
        decay_corrected=decay_applied,
        decay_factor=decay_factor,
        reference_datetime=reference_datetime,
        measurement_datetime=measurement_datetime,
        T_half_s=T_half,
        notes=notes,
    )


def compute_activities_for_all(
    identification_result,
    *,
    efficiency_curve: EfficiencyCurve,
    live_time_s: float,
    from_bg_subtracted: bool,
    bg_available: bool = False,
    force_gross: bool = False,
    coincidence_corrections: Optional[Dict[str, Dict[float, float]]] = None,
    tcs_method_scale: Optional[Dict[str, float]] = None,
    decay_correction: bool = True,
    reference_datetime: Optional[datetime] = None,
    measurement_datetime: Optional[datetime] = None,
    epsilon_unc_pct: float = 5.0,
    skip_undetected: bool = True,
    # F-122 / v1.17.6 — self-attenuation wiring для volume-geometry образцов
    geometry_canonical: str = "",
    sample_density_g_cm3: Optional[float] = None,
    matrix_composition: Optional[Dict[str, float]] = None,
    # F-296 / v1.18.1 — opt-in auto-TCS correction (preset cascades)
    enable_tcs_correction: bool = False,
    tcs_detector_id: str = "Gamma-1S",
    # F-294 / v1.18.1 — opt-in Cutshall self-absorption fallback
    enable_cutshall_self_abs: bool = False,
    cutshall_path_cm: Optional[float] = None,
    cutshall_calib_density_g_cm3: float = 1.0,
    # F-297 / v1.18.2 — opt-in matrix-method simultaneous solver
    enable_matrix_method: bool = False,
    matrix_method_energy_tolerance_keV: float = 1.0,
    # PTB-2018 Annex E — Pb-214/Bi-214 (E.1) / Pb-212 (E.2) T½ selection mode
    chain_decay_mode: str = "equilibrium",
) -> list:
    """
    Compute activity for every detected nuclide in an IdentificationResult.

    Args:
        identification_result: result from
            `identify_nuclides(...)`.
        coincidence_corrections: dict mapping nuclide name →
            per-line correction dict (e.g.
            {"Co-60": {1173.23: 1.05, 1332.49: 1.05}}).
        skip_undetected: if True (default), only iterate over
            `detected_nuclides`; rejected nuclides yield no
            ActivityResult.
        Other args: passed through to `compute_activity`.

    Returns:
        list[ActivityResult], one per detected nuclide.
    """
    coincidence_corrections = coincidence_corrections or {}
    nuclides_to_process = list(
        getattr(identification_result, "detected_nuclides", ())
    )
    if not skip_undetected:
        nuclides_to_process += list(
            getattr(identification_result, "rejected_nuclides", ())
        )

    # F-122 / v1.17.6 — построить таблицу per-line self-attenuation
    # факторов один раз для всех нуклидов. Активируется только когда:
    #   1) geometry_canonical зарегистрирована в REF_GEOMETRY
    #   2) sample_density_g_cm3 передана
    #   3) sample_density != reference density (иначе корректно = 1.0)
    self_att_factors: Optional[Dict[float, float]] = None
    if geometry_canonical and sample_density_g_cm3 is not None:
        try:
            from gamma.physics.self_attenuation import (
                REF_GEOMETRY, correction_factor, OISN_16_COMPOSITION,
            )
            geom_info = REF_GEOMETRY.get(geometry_canonical)
            if geom_info is None:
                # Попробовать нормализованное имя (Маринелли с большой буквы)
                for key in REF_GEOMETRY:
                    if key.lower() == geometry_canonical.lower():
                        geom_info = REF_GEOMETRY[key]
                        break
            if geom_info is not None and sample_density_g_cm3 > 0:
                _, rho_ref, thickness = geom_info
                composition = matrix_composition or OISN_16_COMPOSITION
                # Собрать все уникальные E_keV из идентификации
                all_E: set = set()
                for ni in nuclides_to_process:
                    for m in getattr(ni, "matched_lines", ()):
                        E = float(getattr(m, "library_E_keV", 0.0))
                        if E > 0:
                            all_E.add(round(E, 2))
                self_att_factors = {}
                for E in all_E:
                    self_att_factors[E] = correction_factor(
                        E, rho_sample_g_cm3=sample_density_g_cm3,
                        rho_ref_g_cm3=rho_ref,
                        thickness_cm=thickness,
                        composition=composition,
                    )
        except Exception:
            self_att_factors = None

    # F-294 / v1.18.1 — Cutshall analytic self-absorption fallback.
    # Активируется ТОЛЬКО когда:
    #   1) enable_cutshall_self_abs=True
    #   2) self_att_factors не построены через F-122 (no REF_GEOMETRY match)
    #   3) sample_density_g_cm3 передана
    # Использует NIST XCOM water μ/ρ table — water-equivalent approximation.
    if (enable_cutshall_self_abs and not self_att_factors
            and sample_density_g_cm3 is not None):
        try:
            from gamma.activity.self_absorption import (
                batch_self_absorption_factors,
            )
            from gamma.activity.self_absorption import (
                MARINELLI_05L_MEAN_PATH_CM as _DEFAULT_PATH,
            )
            path_cm = (cutshall_path_cm
                       if cutshall_path_cm is not None else _DEFAULT_PATH)
            all_E: set = set()
            for ni in nuclides_to_process:
                for m in getattr(ni, "matched_lines", ()):
                    E = float(getattr(m, "library_E_keV", 0.0))
                    if E > 0:
                        all_E.add(round(E, 2))
            if all_E:
                E_list = sorted(all_E)
                f_abs_list = batch_self_absorption_factors(
                    energies_keV=E_list,
                    rho_sample=sample_density_g_cm3,
                    rho_calib=cutshall_calib_density_g_cm3,
                    mean_path_cm=path_cm,
                )
                # Cutshall: A_corrected = A_apparent / f_abs.
                # F-122 семантика: sa_factor = F_ref/F_sample → multiplicative.
                # Эквивалент: sa_factor = 1.0 / f_abs.
                self_att_factors = {
                    E: (1.0 / float(f)) if float(f) > 0 else 1.0
                    for E, f in zip(E_list, f_abs_list)
                }
        except Exception:
            self_att_factors = None

    # F-297 / v1.18.2 — opt-in matrix-method simultaneous solver.
    # Решает все нуклиды одновременно через WLS χ²-минимизацию
    # (alternative к per-nuclide weighted-mean). Требует ≥2 нуклидов
    # и ≥2 уникальных пиков (иначе singular matrix). На сбое — fallback
    # на стандартный путь.
    if enable_matrix_method and len(nuclides_to_process) >= 2:
        try:
            from gamma.activity.matrix_method_chi2 import (
                NuclideContribution, PeakObservation,
                solve_matrix_method,
            )
            peak_dict: Dict[float, PeakObservation] = {}
            contributions_dict: Dict[str, List[NuclideContribution]] = {}
            for ni in nuclides_to_process:
                nuc_name = getattr(ni, "nuclide", "")
                if not nuc_name:
                    continue
                lines_for_nuc: List[NuclideContribution] = []
                for m in getattr(ni, "matched_lines", ()):
                    E = float(getattr(m, "library_E_keV", 0.0))
                    I_pct = float(getattr(m, "library_I_pct", 0.0))
                    S = getattr(m, "peak_area", None)
                    if E <= 0 or I_pct <= 0 or S is None or S <= 0:
                        continue
                    eps = efficiency_curve.efficiency_at(E)
                    if eps is None or eps <= 0:
                        continue
                    key = round(E, 1)
                    if key not in peak_dict:
                        sigma_S_attr = getattr(
                            m, "peak_area_uncertainty", None,
                        )
                        if sigma_S_attr and sigma_S_attr > 0:
                            counts_bg_proxy = max(
                                (sigma_S_attr ** 2) - float(S), 0.0,
                            )
                        else:
                            counts_bg_proxy = 0.0
                        peak_dict[key] = PeakObservation(
                            E_keV=E, counts=float(S),
                            counts_bg=counts_bg_proxy,
                        )
                    lines_for_nuc.append(NuclideContribution(
                        nuclide=nuc_name, E_keV=E,
                        intensity_decimal=I_pct / 100.0,
                        efficiency=eps,
                        live_time_seconds=live_time_s,
                    ))
                if lines_for_nuc:
                    contributions_dict[nuc_name] = lines_for_nuc

            if (len(peak_dict) >= len(contributions_dict)
                    and len(contributions_dict) >= 2):
                mat = solve_matrix_method(
                    peaks=list(peak_dict.values()),
                    contributions=contributions_dict,
                    energy_tolerance_keV=matrix_method_energy_tolerance_keV,
                )
                matrix_results: list = []
                note_base = (
                    f"F-297 / v1.18.2 matrix_method "
                    f"χ²_red={mat.chi2_reduced:.2f}, "
                    f"is_acceptable={mat.is_acceptable}"
                )
                if mat.needs_more_nuclides:
                    note_base += " — needs_more_nuclides (χ²_red>3.0)"
                for nuc_name, A_Bq in mat.activities_Bq.items():
                    sigma = float(
                        mat.activity_uncertainties_Bq.get(nuc_name, 0.0)
                    )
                    matrix_results.append(ActivityResult(
                        nuclide=nuc_name,
                        A_Bq=float(A_Bq), sigma_A_Bq=sigma,
                        lines_used=(),
                        lines_skipped=(),
                        intra_chi2_per_dof=float(mat.chi2_reduced),
                        sigma_method="matrix_method",
                        from_bg_subtracted=from_bg_subtracted,
                        force_gross_override=False,
                        notes=note_base,
                    ))
                return matrix_results
        except Exception as exc:  # DEEP-06
            warnings.warn(
                f"matrix_method activity solve failed ({exc!r}); falling "
                f"back to per-nuclide weighted-mean — multi-nuclide "
                f"χ² coupling will not be applied to this result.",
                stacklevel=2,
            )

    # ─── BUG-15 / v1.18.31 — cross-nuclide peak-ownership map ──────────
    # When two different nuclides match the SAME observed peak (channel),
    # the photopeak counts physically belong to ONE transition. The
    # standard rule: the characteristic line of nuclide X at peak P
    # "owns" P; any non-characteristic line of nuclide Y at the same P
    # is a spurious energy-proximity match (typical when NaI FWHM is
    # wide enough that a low-I library line of Y falls within Δ_FWHM of
    # X's strong line — e.g. Ac-228 562 keV claiming the Tl-208 583 keV
    # peak, or Ac-228 145 keV claiming the Tl-208 233 keV neighbor).
    #
    # Rule: for each peak_channel, find the set of characteristic-line
    # owners (nuclides where the line on that channel has is_characteristic
    # = True). For every other matched line on that channel whose host
    # nuclide is NOT in that owner set, mark (nuclide, library_E_keV)
    # for skipping in compute_activity.
    #
    # Within-nuclide deduplication (multiple library lines of the same
    # nuclide sharing one channel) is still handled inside
    # compute_activity (see _shared_skip_keys block there).
    cross_nuclide_skip: Dict[str, set] = {}  # nuclide → set of rounded E_keV
    channel_characteristic_owners: Dict[int, set] = {}
    channel_matches: Dict[int, list] = {}
    for ni in nuclides_to_process:
        nuc_name = getattr(ni, "nuclide", "")
        for m in getattr(ni, "matched_lines", ()):
            ch = int(getattr(m, "peak_channel", -1) or -1)
            S = getattr(m, "peak_area", None)
            if ch < 0 or S is None or float(S) <= 0:
                continue
            channel_matches.setdefault(ch, []).append((nuc_name, m))
            if bool(getattr(m, "is_characteristic", False)):
                channel_characteristic_owners.setdefault(ch, set()).add(nuc_name)
    for ch, matches_here in channel_matches.items():
        nuclides_on_ch = {n for n, _ in matches_here}
        if len(nuclides_on_ch) < 2:
            continue
        owners = channel_characteristic_owners.get(ch, set())
        if not owners:
            # No characteristic owner — fall back to "highest I_pct nuclide
            # wins". Compute max I per nuclide on this channel.
            max_I_by_nuc: Dict[str, float] = {}
            for n, m in matches_here:
                cur = max_I_by_nuc.get(n, 0.0)
                I = float(getattr(m, "library_I_pct", 0.0) or 0.0)
                if I > cur:
                    max_I_by_nuc[n] = I
            if not max_I_by_nuc:
                continue
            top_nuc = max(max_I_by_nuc, key=lambda k: max_I_by_nuc[k])
            owners = {top_nuc}
        # Mark non-owner nuclides' lines on this channel for skip.
        for n, m in matches_here:
            if n in owners:
                continue
            E = round(float(getattr(m, "library_E_keV", 0.0) or 0.0), 2)
            cross_nuclide_skip.setdefault(n, set()).add(E)

    results = []
    for ni in nuclides_to_process:
        nuc_name = getattr(ni, "nuclide", "")
        cc = coincidence_corrections.get(nuc_name)
        try:
            res = compute_activity(
                ni,
                efficiency_curve=efficiency_curve,
                live_time_s=live_time_s,
                from_bg_subtracted=from_bg_subtracted,
                bg_available=bg_available,
                force_gross=force_gross,
                coincidence_correction=cc,
                tcs_method_scale=tcs_method_scale,
                decay_correction=decay_correction,
                reference_datetime=reference_datetime,
                measurement_datetime=measurement_datetime,
                epsilon_unc_pct=epsilon_unc_pct,
                self_attenuation_factors=self_att_factors,
                enable_tcs_correction=enable_tcs_correction,
                tcs_detector_id=tcs_detector_id,
                cross_nuclide_skip_energies_keV=(
                    cross_nuclide_skip.get(nuc_name)
                ),
                chain_decay_mode=chain_decay_mode,
            )
            results.append(res)
        except BackgroundNotSubtractedError:
            # Re-raise: this is a session-level error that must reach
            # the caller. Returning a partial list would silently
            # under-report.
            raise
    return results


# ============================================================================
# F-119 / v1.17.5 — Chain-equilibrium guard
# ============================================================================

# Канонические члены цепочек распада (для Th-232 и U-238). При
# проверке равновесия рассматриваются только нуклиды этих наборов;
# в равновесии все они дают примерно равную активность на одну атомную
# распад-цепь. Расхождение > 5× указывает на нарушение секулярного
# равновесия (молодой образец / выщелачивание / нарушенный матрикс
# и т.п.) и член-выпад исключается из взвешенного среднего по цепи.
CHAIN_MEMBERS: Dict[str, frozenset] = {
    "Th-232": frozenset({
        "Th-232", "Ac-228", "Tl-208", "Pb-212", "Bi-212", "Th-228", "Ra-224",
    }),
    "U-238": frozenset({
        "U-238", "Bi-214", "Pb-214", "Pb-210", "Ra-226", "Po-214", "Th-234",
    }),
}

# TD-2 / v1.18.29 — chain-head nuclides excluded from equilibrium ratio.
# Long-lived α-emitting chain parents (Th-232: T½=14 Gyr, U-238: T½=4.47 Gyr)
# emit NO direct γ-rays of practical intensity above the natural background.
# Any "activity" reported for them is an artifact of:
#   • Lsrm chain-bundled .lib entries (split_chains residue): the parent
#     keeps a long tail of low-I_pct lines that survive ENSDF reassignment
#     (X-rays, weak unassigned lines). When matched against real daughter
#     photopeaks (233, 583, 911 keV …) the weighted-mean yields a
#     spurious, very small A_Bq that becomes the chain min(A) and inflates
#     ratio = max/min.
#   • Library X-ray escape peaks routed under the parent.
# These activity numbers should never participate in the secular-
# equilibrium ratio test. They MAY remain in the `members` list as
# informational ('excluded_from_ratio': True) for diagnostic
# transparency, but neither min/max/median nor `in_equilibrium` are
# computed from them.
CHAIN_HEADS: frozenset = frozenset({"Th-232", "U-238"})


def chain_equilibrium_guard(
    activities: Iterable["ActivityResult"],
    *,
    chains: Iterable[str] = ("Th-232", "U-238"),
    ratio_threshold: float = 5.0,
) -> dict:
    """F-119: проверка секулярного равновесия по цепям распада.

    Для каждой указанной цепочки:
      • собирает активности её матеров (только валидных, A>0);
      • вычисляет ratio = max(A) / min(A) и median(A);
      • если ratio > ratio_threshold, помечает выбросы как
        equilibrium_flag = "outlier" — это нуклиды с A > 5·median(A)
        (либо A < median(A)/5, симметрично).

    Возвращает словарь, готовый к вставке в
    ``report["diagnostics"]["chain_equilibrium"]``:

      {
        "Th-232": {
          "ratio": ...,
          "median_Bq": ...,
          "members": [{"nuclide": ..., "A_Bq": ..., "is_outlier": bool}, ...],
          "outliers": [...],
        },
        ...
      }

    Также мутирует поле ``equilibrium_flag`` у переданных
    ActivityResult (используется downstream-кодом, если он его
    читает). Так как ActivityResult — frozen dataclass, мы делаем
    «soft tag» через возвращаемый отчёт (поле не существует на
    dataclass-е, чтобы не ломать contracts).

    По умолчанию проверяются Th-232 и U-238. ratio_threshold=5.0
    (типичный предел в γ-спектрометрии: больше пятикратного — это
    нарушение равновесия).
    """
    out: dict = {}
    by_name = {ar.nuclide: ar for ar in activities
               if ar is not None and ar.is_valid()}
    if not by_name:
        return out
    for chain in chains:
        members = CHAIN_MEMBERS.get(chain)
        if not members:
            continue
        # TD-2 / v1.18.29 — split members into ratio-contributing (daughters)
        # and informational (chain heads, see CHAIN_HEADS docstring). Heads
        # never produce direct γ-rays at practical intensity; any A_Bq
        # reported for them is a library artifact (Lsrm split_chains
        # residue / X-ray escape) and must NOT participate in min/max/
        # median.
        chain_acts: list = []           # (name, A_Bq) for daughters
        head_acts: list = []            # (name, A_Bq) for chain heads
        for name in members:
            ar = by_name.get(name)
            if ar is None or not ar.is_valid():
                continue
            if name in CHAIN_HEADS:
                head_acts.append((name, float(ar.A_Bq)))
            else:
                chain_acts.append((name, float(ar.A_Bq)))
        if len(chain_acts) < 2:
            continue
        values = sorted([v for _, v in chain_acts])
        median = values[len(values) // 2]
        if len(values) % 2 == 0:
            median = 0.5 * (values[len(values) // 2 - 1]
                            + values[len(values) // 2])
        if min(values) <= 0:
            continue
        ratio = max(values) / min(values)
        members_block = []
        outliers: List[str] = []
        for name, A in chain_acts:
            is_outlier = False
            if ratio > ratio_threshold and median > 0:
                if A > ratio_threshold * median or A < median / ratio_threshold:
                    is_outlier = True
                    outliers.append(name)
            members_block.append({
                "nuclide": name,
                "A_Bq": A,
                "is_outlier": bool(is_outlier),
                "excluded_from_ratio": False,
            })
        for name, A in head_acts:
            # Heads are reported for transparency but flagged as
            # excluded — they did NOT contribute to ratio/median.
            members_block.append({
                "nuclide": name,
                "A_Bq": A,
                "is_outlier": False,
                "excluded_from_ratio": True,
            })
        out[chain] = {
            "ratio": float(ratio),
            "median_Bq": float(median),
            "members": members_block,
            "outliers": outliers,
            "ratio_threshold": float(ratio_threshold),
            "in_equilibrium": bool(ratio <= ratio_threshold),
        }
    return out


__all__ = [
    "BackgroundNotSubtractedError",
    "LineActivity",
    "ActivityResult",
    "DEFAULT_TCS_METHOD_SCALE",
    "CASCADE_SUMMING_NUCLIDES",
    "CHAIN_MEMBERS",
    "CHAIN_HEADS",
    "compute_activity",
    "compute_activities_for_all",
    "chain_equilibrium_guard",
]
