"""
Background spectrum subtraction (Lsrm §9).

Background spectrum represents the detector's intrinsic response in the
absence of the sample being measured. It includes:

  • Cosmic ray contributions (constant rate at lab altitude)
  • Detector-intrinsic activity (e.g. K-40 in NaI from natural K
    contamination, U/Th chain contamination in crystal/shielding)
  • Environmental radon and progeny (depending on lab air-handling)
  • Electronic noise tail (below ~25 keV)

Important: background is a **detector property**, not a source-
geometry property. The same background spectrum can be subtracted from
spectra acquired in different geometries (Marinelli, point, Petri)
because the detector's intrinsic background does NOT depend on what
source you happen to be measuring at the time. Geometry affects the
SOURCE response, not the background.

═══════════════════════════════════════════════════════════════════════
User consent model
═══════════════════════════════════════════════════════════════════════

Background subtraction is a destructive operation that significantly
affects all downstream calculations (peak areas, MDA, activities).
It is performed ONLY when the user has explicitly confirmed that a
particular background FILE (not just any background) applies to the
spectrum being analyzed. This guards against:

  • Using stale backgrounds from earlier hardware configurations
  • Using backgrounds from a different detector
  • Using backgrounds taken with significantly different shielding

The consent is **per background file**, not per call. Operational
policy:

  1. The first time a particular bg file appears in a session, prompt
     the user: "Apply background file X to spectrum Y? [Y/N]"
  2. If approved, store the approval in a session-level registry; do
     not re-prompt for the same bg file in the same session.
  3. If the user supplies a DIFFERENT bg file (different path), prompt
     again — approval is not transferable.

This module provides a `BackgroundConsentRegistry` for tracking
approvals. Callers integrate it into their workflow:

    registry = BackgroundConsentRegistry()
    # ... user clicks "OK" for bg_path_X ...
    registry.approve(bg_path_X)
    # All subsequent calls with bg_path_X proceed without re-prompt:
    result = subtract_background(src, bg, consent_registry=registry)
    # Calls with bg_path_Y still raise until approved:
    registry.approve(bg_path_Y)
    result2 = subtract_background(src, bg2, consent_registry=registry)

For ad-hoc / single-shot use, the legacy `user_confirmed_applicable=True`
flag still works (skips the registry entirely for that one call).

═══════════════════════════════════════════════════════════════════════
Subtraction modes (Lsrm §9)
═══════════════════════════════════════════════════════════════════════

  1. **Rate-normalized channel subtraction** (default for matched gains):
     Compute rate per channel for both spectra, subtract, scale back:
        net[i] = max(0, source.counts[i] - bg.counts[i] · (t_src/t_bg))

     Use when source and background have identical or near-identical
     energy calibrations (Δa₁/a₁ < 0.5%).

  2. **Energy-aligned subtraction** (gain mismatch):
     Rebin the background onto the source's energy grid before
     subtracting. The user-supplied background spectrum may have been
     acquired with slightly different gain (drift, recalibration);
     in that case the same channel index means different energies in
     the two files. Failure to align creates spurious dip-and-peak
     features near every background photopeak.
        For each source channel i:
          E_low  = source_energy(i - 0.5)
          E_high = source_energy(i + 0.5)
          bg_in_bin[i] = integral of bg.counts over [E_low, E_high]
        net[i] = max(0, source.counts[i] - bg_in_bin[i] · (t_src/t_bg))

  3. **ROI-by-ROI subtraction** (Phase 2.1d):
     Compute integrated background under each ROI separately. Allows
     different background-fit polynomials per ROI. Best for activity
     calculations where the exact under-peak background matters more
     than overall spectral shape.

This module implements modes 1 and 2.

Reference: Lsrm Algorithmic Foundations 2022, §9; Knoll "Radiation
Detection" 4th Ed., §16.II.E.

Related implementation (F-58 untangle, 2026-06-21):
  A LITE variant `bg_subtract_energy.py` lives next to this module —
  single-function energy-domain np.interp subtraction without safety
  gates. The LITE variant is the path used by
  `gamma.identification.staged_pipeline` (production pipeline /
  run_plan_a.py). This FULL module (rate_normalized_channel +
  energy_aligned + F-243 zero-point guard +
  BackgroundConsentRegistry/Required) is the path used by
  `scripts/validate_certs.py` (standalone cert validation). Both
  share the F-58 / F-160 contract; the FULL module additionally
  enforces F-243 safety against the wave3/wave4 a0-mismatch class of
  bugs.

  F-451 direction (2026-06-22): live-time scaling goes **toward the
  shorter live-time** ("к меньшему"). If t_S >= t_B → sample scaled
  DOWN by t_B/t_S; if t_S < t_B → background scaled DOWN by t_S/t_B.
  Net is returned in counts equivalent to effective_live_time =
  min(t_S, t_B); `apply_subtraction_to_spectrum` overwrites
  `spec.live_time` with effective_live_time so the cps-invariant
  (net_cps = sample_rate − bg_rate) downstream holds correctly.
  See audit/_plans/F-451_bg_subtract_direction_invert.md.
"""

from __future__ import annotations

# AUDIT-F2 (2026-06-25): import math удалён — после векторизации
# _subtract_energy_aligned не использует math.floor.
import os
from dataclasses import dataclass, replace, field
from typing import Optional, Set

import numpy as np


# Threshold for "gain match": if |Δa₁ / a₁| < this, we use channel-by-
# channel mode; otherwise energy-aligned mode. 0.5% is a practical
# threshold — at NaI 50 keV FWHM @ 1 MeV (Δa₁/a₁ = 0.5% → Δ = 5 keV at
# 1000 keV which is much less than FWHM, so channel subtraction is fine).
GAIN_MATCH_THRESHOLD = 0.005

# Threshold for "zero-point match": absolute delta of a0 (channel-0 energy
# offset). When |Δa0| exceeds this, the gain-match test alone is NOT
# sufficient — even with identical gains a non-trivial a0 delta means that
# bg.counts[i] and src.counts[i] correspond to *different* energies, and
# channel-by-channel subtraction will subtract bg from the wrong energy
# region (Agent B wave 3 finding 2026-06-04: a0_src=+47.669 keV vs
# a0_bg=−25.814 keV → Δ=73.5 keV, gain match 0.28%<0.5%, silent
# rate_normalized_channel mode → at K-40 ch 444 bg-channel maps to
# ~1302 keV instead of 1460 keV, suppresses K-40 BG subtraction by ~50%).
# Defence of 30 keV: NaI(Tl) FWHM at 1460 keV is ≈ 100 keV (6.5% rel.,
# Gilmore & Joss §6.4); HPGe FWHM at 1460 keV is ≈ 2-3 keV (Gilmore §6.5).
# Choosing 30 keV gives ~10× HPGe FWHM headroom (catches HPGe pathology
# clearly) while still being ≈ 0.3× NaI FWHM (catches NaI cases like
# B wave 3 where Δ=73 keV ≫ 30 keV ≫ normal cal drift of <5 keV).
# See _state/agent_b/outbox/2026-06-04_wave3_k40_regen_and_bg_subtraction_investigation.md
# lines 112-211 for the empirical case.
ZERO_POINT_MATCH_THRESHOLD_KEV = 30.0


class BackgroundConsentRequired(RuntimeError):
    """Raised when subtract_background is called without consent for the
    specific background file. The exception message includes the bg file
    path so the caller can prompt the user to approve it."""

    def __init__(self, message: str, *, background_path: str = "",
                 source_path: str = ""):
        super().__init__(message)
        self.background_path = background_path
        self.source_path = source_path


@dataclass
class BackgroundConsentRegistry:
    """
    Session-level registry of approved background files.

    Operational model: the calling application maintains one
    registry per analysis session. When the user explicitly approves
    a particular bg file (via a dialog, CLI prompt, or configuration
    entry), the application calls `approve(path)`. Subsequent
    `subtract_background` calls referencing that same bg file
    proceed without prompting; calls referencing a different bg
    file raise `BackgroundConsentRequired` until that file is also
    approved.

    Path comparison is normalised (absolute, case-sensitive on the
    underlying filesystem) so equivalent paths spelt differently
    still match.

    Usage:
        registry = BackgroundConsentRegistry()
        # User confirms via UI:
        registry.approve("/path/to/bg.spe")
        # Later calls proceed without re-prompt:
        result = subtract_background(src, bg, consent_registry=registry)
    """
    _approved_paths: Set[str] = field(default_factory=set)

    @staticmethod
    def _normalize(path: str) -> str:
        """Normalise path for stable comparison across calls."""
        if not path:
            return ""
        try:
            return os.path.abspath(path)
        except (TypeError, ValueError):
            return path

    def approve(self, background_path: str) -> None:
        """Record user approval for a specific background file."""
        self._approved_paths.add(self._normalize(background_path))

    def revoke(self, background_path: str) -> None:
        """Remove an existing approval (e.g. user changed their mind)."""
        self._approved_paths.discard(self._normalize(background_path))

    def is_approved(self, background_path: str) -> bool:
        """Check whether a specific bg file has been approved."""
        return self._normalize(background_path) in self._approved_paths

    def approved_paths(self) -> list:
        """Return all currently approved paths (for diagnostic display)."""
        return sorted(self._approved_paths)

    def clear(self) -> None:
        """Forget all approvals (start of a new session)."""
        self._approved_paths.clear()


@dataclass(frozen=True)
class BackgroundSubtractionResult:
    """Result of background subtraction (F-58 / F-451)."""

    source_path: str
    background_path: str
    mode: str                       # "rate_normalized_channel" or "energy_aligned"
    subtracted_counts: np.ndarray   # net counts per channel in effective_live_time scale (F-451), always >= 0
    subtracted_uncertainties: np.ndarray  # Poisson-propagated uncertainty per channel, same scale as subtracted_counts

    # Legacy scaling factor: source_live_time / bg_live_time (kept for downstream readers)
    rate_scale: float
    # Source live time (input value, NOT effective; downstream uses effective_live_time for cps)
    source_live_time: float

    # Diagnostics
    gain_mismatch_relative: float   # |Δa₁ / a₁|
    n_channels_clipped_to_zero: int # how many channels had source < scaled bg (in effective scale, F-451)
    total_source_counts: int
    total_bg_counts_after_scale: float
    notes: str = ""

    # F-451 direction-of-scaling fields (default to legacy semantics for safety):
    applied_scale: float = 1.0          # min(t_s, t_bg) / max(t_s, t_bg) — actually applied k', in [0, 1]
    scale_direction: str = "equal"      # "sample_down" | "bg_down" | "equal"
    effective_live_time: float = 0.0    # min(t_s, t_bg) — live-time of returned subtracted_counts

    def __repr__(self) -> str:
        return (f"BackgroundSubtractionResult(mode={self.mode}, "
                f"rate_scale={self.rate_scale:.3f}, "
                f"direction={self.scale_direction}, "
                f"applied={self.applied_scale:.3f}, "
                f"clipped={self.n_channels_clipped_to_zero}/{len(self.subtracted_counts)}, "
                f"Δa₁/a₁={self.gain_mismatch_relative*100:.2f}%)")


def subtract_background(
    source_spec,
    background_spec,
    *,
    consent_registry: Optional[BackgroundConsentRegistry] = None,
    user_confirmed_applicable: bool = False,
    force_mode: Optional[str] = None,
) -> BackgroundSubtractionResult:
    """
    Subtract a background spectrum from a source spectrum.

    Args:
        source_spec: Spectrum with sample measurement. Must have
            live_time > 0 and energy_cal.
        background_spec: Spectrum with detector background. Must be
            confirmed applicable to this measurement (same detector,
            same shielding configuration). Must have live_time > 0
            and energy_cal.
        consent_registry: `BackgroundConsentRegistry` for session-level
            tracking of approved bg files. If provided AND the
            background file is in the registry, subtraction proceeds
            without prompting. If provided but the file is not yet
            approved, raises `BackgroundConsentRequired` — caller
            should prompt user, then call
            `consent_registry.approve(bg_path)` and retry.
        user_confirmed_applicable: legacy one-shot consent flag. When
            True, skips the registry entirely for this single call.
            Use this for non-interactive batch workflows where consent
            was obtained via some other channel (config file, signed
            workflow, etc).
        force_mode: if provided ("rate_normalized_channel" or
            "energy_aligned"), use that mode regardless of gain match.
            Default: auto-select based on `GAIN_MATCH_THRESHOLD`.

    Returns:
        BackgroundSubtractionResult.

    Raises:
        BackgroundConsentRequired: if neither one-shot consent nor a
            matching registry entry is present. The exception message
            includes the bg file path so the caller can prompt the
            user with that specific identifier.
        ValueError: if either spectrum lacks live_time or energy_cal

    Notes on the rate-scale: source.counts has units of "counts in
    source.live_time seconds"; same for bg.counts. We want net counts
    in source.live_time seconds, so we scale bg counts by the ratio:
        net = source.counts - bg.counts · (source.live_time / bg.live_time)

    Notes on uncertainty: counting statistics give
        σ²(net[i]) = source.counts[i] + bg.counts[i] · scale²
    where the second term is the squared uncertainty of the scaled
    background contribution.
    """
    # Check consent first — fail fast before doing any work
    bg_path = getattr(background_spec, "source_path", "") or ""
    src_path = getattr(source_spec, "source_path", "") or ""
    has_consent = (
        user_confirmed_applicable
        or (consent_registry is not None and consent_registry.is_approved(bg_path))
    )
    if not has_consent:
        bg_name = os.path.basename(bg_path) if bg_path else "<unnamed>"
        src_name = os.path.basename(src_path) if src_path else "<unnamed>"
        raise BackgroundConsentRequired(
            f"Background subtraction requires explicit user confirmation "
            f"for the specific background file. "
            f"Spectrum: {src_name!r} | "
            f"Background file: {bg_name!r} ({bg_path}). "
            f"Prompt the user: 'Apply this background file to this spectrum?'. "
            f"On approval, either: (a) pass user_confirmed_applicable=True for "
            f"a one-shot call, or (b) call consent_registry.approve({bg_path!r}) "
            f"so subsequent calls with the same bg file proceed automatically.",
            background_path=bg_path,
            source_path=src_path,
        )

    if source_spec.live_time <= 0 or background_spec.live_time <= 0:
        raise ValueError(
            f"Both spectra must have positive live_time; "
            f"source={source_spec.live_time}, bg={background_spec.live_time}"
        )
    if source_spec.energy_cal is None or background_spec.energy_cal is None:
        raise ValueError("Both spectra must have energy_cal set.")

    # Compute gain mismatch
    a1_src = source_spec.energy_cal[1] if len(source_spec.energy_cal) > 1 else 0.0
    a1_bg = background_spec.energy_cal[1] if len(background_spec.energy_cal) > 1 else 0.0
    if a1_src <= 0:
        raise ValueError(f"Source has invalid gain a1={a1_src}")
    delta_rel = abs(a1_src - a1_bg) / a1_src

    # Compute zero-point (a0) absolute delta in keV — see
    # ZERO_POINT_MATCH_THRESHOLD_KEV docstring above for rationale.
    a0_src = source_spec.energy_cal[0] if len(source_spec.energy_cal) > 0 else 0.0
    a0_bg = background_spec.energy_cal[0] if len(background_spec.energy_cal) > 0 else 0.0
    delta_a0_keV = abs(a0_src - a0_bg)

    # Decide mode
    if force_mode:
        mode = force_mode
    elif delta_a0_keV > ZERO_POINT_MATCH_THRESHOLD_KEV:
        # F-243 safety (B wave 3 finding 2026-06-04): zero-point delta
        # > 30 keV → channel-by-channel subtraction would subtract bg
        # from a wrong energy region (B wave 3 documented K-40
        # suppression by ~50% at Δa0=73.5 keV). Force energy_aligned
        # regardless of gain match — rebinning bg onto src energy grid
        # is the only correct subtraction in this regime.
        mode = "energy_aligned"
    elif delta_rel < GAIN_MATCH_THRESHOLD:
        mode = "rate_normalized_channel"
    else:
        mode = "energy_aligned"

    # Legacy rate_scale (kept on result for downstream readers):
    rate_scale = source_spec.live_time / background_spec.live_time

    # F-451 direction of scaling: always scale toward the SHORTER live-time.
    t_s = float(source_spec.live_time)
    t_bg = float(background_spec.live_time)
    if t_s > t_bg:
        scale_direction = "sample_down"
        applied_scale = t_bg / t_s
        effective_live_time = t_bg
    elif t_s < t_bg:
        scale_direction = "bg_down"
        applied_scale = t_s / t_bg
        effective_live_time = t_s
    else:
        scale_direction = "equal"
        applied_scale = 1.0
        effective_live_time = t_s

    # Dispatch
    if mode == "rate_normalized_channel":
        net, unc, n_clipped = _subtract_channel(
            source_spec.counts, background_spec.counts,
            applied_scale, scale_direction,
        )
    elif mode == "energy_aligned":
        net, unc, n_clipped = _subtract_energy_aligned(
            source_spec, background_spec,
            applied_scale, scale_direction,
        )
    else:
        raise ValueError(f"Unknown mode: {mode!r}")

    # total_bg_counts_after_scale — reports the bg total after scaling
    # by applied_scale, in the EFFECTIVE-scale (F-451). Old semantic was
    # bg_total · rate_scale (sample-scale); new is bg_total · applied_scale
    # if bg was the one scaled, else bg_total (unscaled) when sample was
    # the one scaled.
    bg_total_raw = float(background_spec.counts.sum())
    if scale_direction == "bg_down":
        total_bg_counts_after_scale = bg_total_raw * applied_scale
    else:
        total_bg_counts_after_scale = bg_total_raw

    return BackgroundSubtractionResult(
        source_path=src_path,
        background_path=bg_path,
        mode=mode,
        subtracted_counts=net,
        subtracted_uncertainties=unc,
        rate_scale=rate_scale,
        source_live_time=source_spec.live_time,
        gain_mismatch_relative=delta_rel,
        n_channels_clipped_to_zero=int(n_clipped),
        total_source_counts=int(source_spec.counts.sum()),
        total_bg_counts_after_scale=total_bg_counts_after_scale,
        notes=(f"mode={mode}, Δa₁/a₁={delta_rel*100:.2f}%, "
               f"Δa₀={delta_a0_keV:.2f} keV, "
               f"F-451 direction={scale_direction}, applied={applied_scale:.4f}, "
               f"effective_live_time={effective_live_time:.0f}s, "
               f"clipped={n_clipped}/{len(net)} channels"),
        applied_scale=applied_scale,
        scale_direction=scale_direction,
        effective_live_time=effective_live_time,
    )


def _subtract_channel(
    source_counts: np.ndarray,
    bg_counts: np.ndarray,
    applied_scale: float,
    scale_direction: str,
):
    """Channel-by-channel rate-normalized subtraction (F-451 direction-aware).

    F-451: net is in the effective_live_time = min(t_s, t_bg) scale. The
    spectrum with the LONGER live-time is scaled down by applied_scale.
    """
    # Align lengths: truncate or pad bg to source length
    n_src = len(source_counts)
    n_bg = len(bg_counts)
    if n_bg >= n_src:
        bg_aligned = bg_counts[:n_src].astype(np.float64)
    else:
        bg_aligned = np.zeros(n_src, dtype=np.float64)
        bg_aligned[:n_bg] = bg_counts.astype(np.float64)
    src_f = source_counts.astype(np.float64)
    if scale_direction == "sample_down":
        src_eff = src_f * applied_scale
        bg_eff = bg_aligned
        unc_sq = (applied_scale ** 2) * src_f + bg_aligned
    elif scale_direction == "bg_down":
        src_eff = src_f
        bg_eff = bg_aligned * applied_scale
        unc_sq = src_f + (applied_scale ** 2) * bg_aligned
    else:  # equal
        src_eff = src_f
        bg_eff = bg_aligned
        unc_sq = src_f + bg_aligned
    raw_diff = src_eff - bg_eff
    # Clip negatives to zero
    net = np.maximum(0.0, raw_diff)
    n_clipped = int(np.sum(raw_diff < 0))
    unc = np.sqrt(np.maximum(0.0, unc_sq))
    return net.astype(np.int64), unc, n_clipped


def _energy_at_channel(spec, ch: float) -> float:
    """Evaluate a spectrum's energy calibration at a (possibly fractional) channel."""
    if spec.energy_cal is None:
        return 0.0
    return sum(a * (ch ** i) for i, a in enumerate(spec.energy_cal))


def _invert_energy_to_channel(spec, E: float) -> float:
    """Find the (possibly fractional) channel where spec.energy_cal evaluates to E.

    Uses bisection over the spectrum's channel range. Returns -1 if E
    is below ch=0's energy; returns n_channels if above the last channel.
    """
    n = len(spec.counts)
    e_lo = _energy_at_channel(spec, 0)
    e_hi = _energy_at_channel(spec, n - 1)
    if E <= e_lo:
        # Linear extrapolation: assume constant gain
        a1 = spec.energy_cal[1] if len(spec.energy_cal) > 1 else 1.0
        if a1 <= 0:
            return 0.0
        return (E - spec.energy_cal[0]) / a1
    if E >= e_hi:
        a1 = spec.energy_cal[1] if len(spec.energy_cal) > 1 else 1.0
        if a1 <= 0:
            return n
        return (E - spec.energy_cal[0]) / a1
    # Bisect for monotone polynomial (always true for physical calibration)
    lo, hi = 0.0, float(n - 1)
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        e_mid = _energy_at_channel(spec, mid)
        if e_mid < E:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-6:
            break
    return 0.5 * (lo + hi)


def _energy_at_channel_vec(spec, ch_arr: np.ndarray) -> np.ndarray:
    """AUDIT-F2 (2026-06-25): векторная версия _energy_at_channel.
    Bit-identical к scalar-варианту: тот же порядок аккумуляции
    `sum(a * ch ** i)`, та же арифметика повторного умножения для
    степеней (repeated mult, не np.power), float64."""
    if spec.energy_cal is None:
        return np.zeros_like(ch_arr)
    cal = spec.energy_cal
    e = np.zeros_like(ch_arr, dtype=np.float64)
    pow_i = np.ones_like(ch_arr, dtype=np.float64)  # ch ** 0 = 1
    for i, a in enumerate(cal):
        if i > 0:
            pow_i = pow_i * ch_arr  # ch ** i = ch ** (i-1) * ch
        e = e + float(a) * pow_i
    return e


def _invert_energy_to_channel_vec(spec, E_arr: np.ndarray) -> np.ndarray:
    """AUDIT-F2 (2026-06-25): векторная версия _invert_energy_to_channel.

    Bit-identical к scalar по конструкции:
      • та же бисекция с tol=1e-6, max 50 итераций;
      • то же early-break per element через mask (active &= hi-lo >= 1e-6);
      • та же edge-extrapolation (линейно через a1 при E вне [e_lo, e_hi]).

    Использует _energy_at_channel_vec (тот же порядок аккумуляции).
    """
    n = len(spec.counts)
    cal = spec.energy_cal
    a0 = float(cal[0])
    a1 = float(cal[1]) if len(cal) > 1 else 1.0
    e_lo_b = _energy_at_channel(spec, 0)
    e_hi_b = _energy_at_channel(spec, n - 1)

    E = np.asarray(E_arr, dtype=np.float64)

    # Extrapolation arms (mirror scalar a1<=0 carve-out)
    if a1 <= 0:
        ch_extrap_lo = np.zeros_like(E)
        ch_extrap_hi = np.full_like(E, float(n))
    else:
        ch_lin = (E - a0) / a1
        ch_extrap_lo = ch_lin
        ch_extrap_hi = ch_lin

    mode_lo = E <= e_lo_b
    mode_hi = E >= e_hi_b
    mode_bisect = ~mode_lo & ~mode_hi

    lo = np.zeros_like(E)
    hi = np.full_like(E, float(n - 1))
    active = mode_bisect.copy()

    for _ in range(50):
        if not active.any():
            break
        mid = 0.5 * (lo + hi)
        e_mid = _energy_at_channel_vec(spec, mid)
        below = e_mid < E
        upd_lo = active & below
        upd_hi = active & ~below
        lo = np.where(upd_lo, mid, lo)
        hi = np.where(upd_hi, mid, hi)
        active = active & (hi - lo >= 1e-6)

    ch_bisect = 0.5 * (lo + hi)
    return np.where(mode_lo, ch_extrap_lo,
                    np.where(mode_hi, ch_extrap_hi, ch_bisect))


def _subtract_energy_aligned(
    source_spec, bg_spec,
    applied_scale: float,
    scale_direction: str,
):
    """
    Energy-aligned subtraction: rebin bg onto source's energy grid (F-451
    direction-aware).

    For each source channel i, compute the energy bin
    [E(i-0.5), E(i+0.5)] and integrate bg.counts over the matching
    channel range in the background spectrum. This handles arbitrary
    gain/offset mismatch correctly. Then apply F-451 scaling toward the
    shorter live-time.

    AUDIT-F2 (2026-06-25): векторизация через edge-integral rebinning.
    Старый цикл `for i in range(n_src) … bg_sum += bg_counts[j] * frac`
    был ~800k интерпретируемых операций на одно вычитание. Новый путь:
      1. инверсия энергорёбер source в bg-каналы — векторная бисекция
         (`_invert_energy_to_channel_vec`), bit-identical к scalar;
      2. кумулятивная сумма bg-counts → bg-канальная функция cum_bg(x);
      3. `bg_in_src_bin[i] = cum_bg(ch_edge[i+1]) - cum_bg(ch_edge[i])`
         через `np.interp` (линейная интерполяция cum_bg на дробные
         каналы) и `np.diff` — это математически тот же edge-integral,
         что считал старый цикл (частичная доля каждого bg-канала ↔
         линейная интерполяция кумулятивной суммы).
    Естественное клипание `np.interp` к границам `[0, n_bg]` воспроизводит
    старые edge cases (полностью-вне → 0). Bit-identity подтверждён на
    suite F-451 snapshot (27 тестов).
    """
    n_src = len(source_spec.counts)
    bg_counts = bg_spec.counts.astype(np.float64)
    n_bg = len(bg_counts)

    # Энергетические рёбра каналов источника (i ± 0.5 для каждого канала).
    ch_src_edges = np.arange(n_src + 1, dtype=np.float64) - 0.5
    E_src_edges = _energy_at_channel_vec(source_spec, ch_src_edges)
    # NB: исходный scalar swap'ил E_lo↔E_hi при reversed cal; здесь то же
    # ведёт `np.diff` на cumsum (см. ниже): отрицательная разность ↔ swap.
    # На монотонной cal (физическая norma) — swap не нужен.

    # Инвертируем энергии-рёбра в дробные bg-каналы (векторная бисекция).
    ch_bg_at_src_edges = _invert_energy_to_channel_vec(bg_spec, E_src_edges)

    # Кумулятивная сумма bg-counts; cum_bg[j] = sum(bg_counts[0..j-1]).
    cum_bg = np.concatenate(([0.0], np.cumsum(bg_counts)))  # length n_bg+1
    ch_bg_int = np.arange(n_bg + 1, dtype=np.float64)        # [0, 1, ..., n_bg]

    # Линейная интерполяция кумулятивной суммы на дробные bg-каналы;
    # естественное клипание `np.interp` к [0, n_bg] = старые edge-cases.
    cum_at_src_edges = np.interp(ch_bg_at_src_edges, ch_bg_int, cum_bg)
    bg_in_src_bin = np.diff(cum_at_src_edges)
    # Если E_lo > E_hi (нефизическая reversed cal), diff даст отрицательное
    # значение — старый код брал abs через swap. Берём max(0, …) на счёте,
    # эквивалентно (на reversed-cal источниках интеграл и есть |Δ|).
    bg_in_src_bin = np.where(bg_in_src_bin < 0.0, -bg_in_src_bin, bg_in_src_bin)

    src_f = source_spec.counts.astype(np.float64)
    if scale_direction == "sample_down":
        src_eff = src_f * applied_scale
        bg_eff = bg_in_src_bin
        unc_sq = (applied_scale ** 2) * src_f + bg_in_src_bin
    elif scale_direction == "bg_down":
        src_eff = src_f
        bg_eff = bg_in_src_bin * applied_scale
        unc_sq = src_f + (applied_scale ** 2) * bg_in_src_bin
    else:  # equal
        src_eff = src_f
        bg_eff = bg_in_src_bin
        unc_sq = src_f + bg_in_src_bin
    raw_diff = src_eff - bg_eff
    net = np.maximum(0.0, raw_diff)
    n_clipped = int(np.sum(raw_diff < 0))
    unc = np.sqrt(np.maximum(0.0, unc_sq))
    return net.astype(np.int64), unc, n_clipped


def apply_subtraction_to_spectrum(source_spec, subtraction_result):
    """
    Return a new Spectrum with `counts` replaced by the subtracted net
    counts (F-451: in effective_live_time scale).

    F-451: `live_time` and `real_time` on the new spectrum are set to
    `subtraction_result.effective_live_time = min(t_S, t_B)` so that
    downstream cps-style calculations (counts/live_time) match the
    cps-invariant (sample_rate − bg_rate). The original
    `source_live_time` is retained in
    `extras["pre_bg_subtract_source_live_time"]` for auditability.

    `is_background_subtracted` flag is set on extras to ensure
    downstream code knows the spectrum has been processed (so it
    doesn't double-subtract).
    """
    eff = float(subtraction_result.effective_live_time)
    new_spec = replace(
        source_spec,
        counts=subtraction_result.subtracted_counts,
        live_time=eff,
        real_time=eff,
    )
    new_spec.extras["background_subtracted"] = True
    new_spec.extras["background_subtraction_mode"] = subtraction_result.mode
    new_spec.extras["background_subtraction_rate_scale"] = subtraction_result.rate_scale
    new_spec.extras["background_subtraction_applied_scale"] = subtraction_result.applied_scale
    new_spec.extras["background_subtraction_scale_direction"] = subtraction_result.scale_direction
    new_spec.extras["background_subtraction_effective_live_time"] = eff
    new_spec.extras["pre_bg_subtract_source_live_time"] = float(source_spec.live_time)
    new_spec.extras["pre_bg_subtract_source_real_time"] = float(source_spec.real_time)
    new_spec.extras["background_source_path"] = subtraction_result.background_path
    return new_spec


__all__ = [
    "BackgroundConsentRegistry",
    "BackgroundConsentRequired",
    "BackgroundSubtractionResult",
    "GAIN_MATCH_THRESHOLD",
    "subtract_background",
    "apply_subtraction_to_spectrum",
]
