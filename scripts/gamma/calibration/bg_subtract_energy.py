"""
Sample-minus-background subtraction in energy space (F-58 / F-451).

Given a sample spectrum S and a paired background B taken with the
SAME geometry and detector (but possibly different energy calibration
and different live time), produce a net spectrum N on the sample's
channel grid.

F-451 direction (2026-06-22): the live-time scaling goes **toward the
shorter live-time** ("к меньшему"). Concretely:

  • If t_live_S >= t_live_B → SAMPLE is scaled DOWN by k' = t_B / t_S
    (sample_down direction). Net is in counts equivalent to t_B.
  • If t_live_S  < t_live_B → BACKGROUND is scaled DOWN by k' = t_S / t_B
    (bg_down direction). Net is in counts equivalent to t_S.

The resulting effective live-time of the net spectrum is min(t_S, t_B);
result.effective_live_time exposes this and downstream code (e.g.
`staged_pipeline`) overwrites `spec.live_time` with it so that
counts/live_time still yields the correct cps.

Rationale (operator, 2026-06-21): scaling UP to the longer measurement
artificially inflates Poisson σ from the noisier (shorter) spectrum;
scaling DOWN to the shorter measurement preserves "really measured"
counts and σ in their natural magnitude. The cps-invariant
(net_cps = sample_rate − bg_rate) holds in both directions; only the
underlying counts/σ magnitudes change.

Notes:
  • Energy-domain interpolation is the right choice when sample and
    background were taken with slightly different stored calibrations
    (the LSRM detector can drift between long measurements).
  • Channels in the sample whose energy falls outside the background's
    covered range are NOT subtracted (b_resampled = 0 there).
  • Negative net counts are clamped to 0 (statistically valid; they
    arise when noise pushes background above sample in a low-count
    bin).

Related implementation (F-58 untangle, 2026-06-21):
  A FULL variant `bg_subtract_dual_mode.py` lives next to this module
  with explicit two-mode subtraction (rate_normalized_channel +
  energy_aligned), the F-243 ZERO_POINT_MATCH_THRESHOLD_KEV=30 keV
  safety gate, and BackgroundConsentRegistry/Required. The FULL
  variant is the path used by `scripts/validate_certs.py` (standalone
  cert validation). This LITE variant — single energy-domain
  np.interp subtraction without safety gates — is the path used by
  `gamma.identification.staged_pipeline` (production pipeline /
  run_plan_a.py). LITE always energy-aligns, so FULL's gain-match
  safety gate is structurally not applicable here.

  F-451 plan: audit/_plans/F-451_bg_subtract_direction_invert.md.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from gamma.spectrum import Spectrum


@dataclass
class BackgroundSubtractionResult:
    """Output of `subtract_background`.

    F-451 B1 contract (2026-06-22): `background_counts_on_sample_grid`
    is ALWAYS the RAW bg counts re-binned on the sample energy grid (no
    live-time scaling applied). Legacy consumers that compute
    `net = sample − bg_on_sample_grid × scale_factor` keep working
    unchanged. The effective-scale variant (`bg_in_eff`, used for the
    F-451 σ-propagation and HTML 4-way charts) is exposed alongside in
    `background_counts_on_sample_grid_effective`.
    """
    net_counts: np.ndarray              # length = sample.n_channels, in effective_live_time scale (F-451)
    net_cps: np.ndarray                  # net_counts / effective_live_time (cps-invariant)
    background_counts_on_sample_grid: np.ndarray        # RAW bg on sample grid (no scaling). Legacy contract preserved.
    scale_factor: float                  # legacy: t_live_sample / t_live_background (kept for downstream readers)
    overlap_fraction: float              # fraction of sample channels covered
    notes: str = ""
    net_uncertainties: np.ndarray | None = None        # sqrt of Poisson-propagated variance in effective_live_time scale (F-451)
    gain_mismatch_relative: float = float("nan")       # |a1_sample - a1_bg| / max(|a1_sample|, |a1_bg|); nan if either energy_cal lacks index 1
    zero_point_mismatch_keV: float = float("nan")      # abs(a0_sample - a0_bg); nan if either energy_cal lacks index 0
    n_channels_clipped: int = 0                        # count of channels where raw_diff < 0 in effective_live_time scale (BEFORE clamp, F-451)
    # F-451 direction-of-scaling fields:
    applied_scale: float = 1.0                          # min(t_s, t_bg) / max(t_s, t_bg) — actually applied k', in [0, 1]
    scale_direction: str = "equal"                      # "sample_down" | "bg_down" | "equal"
    effective_live_time: float = 0.0                    # min(t_s, t_bg) — the live-time of the returned net counts
    # F-451 B1 (2026-06-22): effective-scale bg on sample grid (bg_in_eff).
    # For sample_down/equal directions: identical to RAW bg on sample grid.
    # For bg_down direction: RAW bg × applied_scale. Used by σ-propagation reconstruction
    # and the HTML 4-way chart toggle (sample / bg-effective / net / etc.).
    background_counts_on_sample_grid_effective: np.ndarray | None = None


def subtract_background(
    sample: Spectrum,
    background: Spectrum,
    *,
    clamp_negative_to_zero: bool = True,
) -> BackgroundSubtractionResult:
    """
    Subtract a background spectrum from a sample, returning net counts
    on the sample's channel grid.

    Args:
        sample: Spectrum to subtract from
        background: Spectrum to subtract (same detector & geometry, but
            independently calibrated)
        clamp_negative_to_zero: if True (default), clamp any negative
            net counts to 0. If False, returns the raw difference
            (useful for statistical residual analysis).

    Returns:
        BackgroundSubtractionResult with net_counts/net_cps arrays
        aligned to `sample`'s channel grid.
    """
    if sample.live_time <= 0 or background.live_time <= 0:
        raise ValueError("Both spectra must have positive live_time")

    # Build energy arrays for both spectra on their own grids
    n_s = sample.n_channels
    e_sample = np.array([sample.channel_to_energy(i) for i in range(n_s)])

    n_b = background.n_channels
    e_bg = np.array([background.channel_to_energy(i) for i in range(n_b)])
    bg_counts = background.counts.astype(float)

    # Live-time scaling (F-451: scale to the SHORTER live-time, "к меньшему").
    t_s = float(sample.live_time)
    t_bg = float(background.live_time)
    scale_factor_legacy = t_s / t_bg   # legacy semantic: kept on result.scale_factor for downstream readers
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

    # Energy-grid interpolation: where sample energy is outside the
    # background's covered range, set to 0 (no subtraction).
    # bg_on_sample is the RAW background counts (no live-time scaling yet),
    # re-binned onto the sample's energy grid.
    e_bg_lo = float(np.min(e_bg))
    e_bg_hi = float(np.max(e_bg))
    bg_on_sample = np.interp(
        e_sample, e_bg, bg_counts, left=0.0, right=0.0,
    )

    # Coverage diagnostic — what fraction of sample channels has
    # actual background data underneath (within the background's
    # energy range)?
    in_range = (e_sample >= e_bg_lo) & (e_sample <= e_bg_hi)
    overlap = float(in_range.sum()) / float(n_s) if n_s else 0.0

    sample_counts = sample.counts.astype(float)

    # F-451: scale the LONGER spectrum down by applied_scale.
    if scale_direction == "sample_down":
        sample_in_eff = sample_counts * applied_scale
        bg_in_eff = bg_on_sample
        # σ²(net) = applied_scale²·N_s_raw + N_bg_raw
        net_uncertainties = np.sqrt(
            (applied_scale ** 2) * sample_counts + bg_on_sample
        )
    elif scale_direction == "bg_down":
        sample_in_eff = sample_counts
        bg_in_eff = bg_on_sample * applied_scale
        # σ²(net) = N_s_raw + applied_scale²·N_bg_raw
        net_uncertainties = np.sqrt(
            sample_counts + (applied_scale ** 2) * bg_on_sample
        )
    else:  # equal
        sample_in_eff = sample_counts
        bg_in_eff = bg_on_sample
        net_uncertainties = np.sqrt(sample_counts + bg_on_sample)

    raw_diff = sample_in_eff - bg_in_eff
    n_channels_clipped = int(np.sum(raw_diff < 0))

    # Compute gain mismatch
    a1_s = sample.energy_cal[1] if len(sample.energy_cal) >= 2 else None
    a1_b = background.energy_cal[1] if len(background.energy_cal) >= 2 else None
    if a1_s is None or a1_b is None:
        gain_mismatch_relative = float("nan")
    else:
        max_a1 = max(abs(a1_s), abs(a1_b))
        if max_a1 == 0:
            gain_mismatch_relative = float("nan")
        else:
            gain_mismatch_relative = abs(a1_s - a1_b) / max_a1

    # Compute zero point mismatch
    a0_s = sample.energy_cal[0] if len(sample.energy_cal) >= 1 else None
    a0_b = background.energy_cal[0] if len(background.energy_cal) >= 1 else None
    if a0_s is None or a0_b is None:
        zero_point_mismatch_keV = float("nan")
    else:
        zero_point_mismatch_keV = abs(a0_s - a0_b)

    net = raw_diff.copy()
    if clamp_negative_to_zero:
        net = np.maximum(net, 0.0)

    # F-451 B1 (2026-06-22): contract preserves legacy semantics — `background_counts_on_sample_grid`
    # is RAW bg on sample grid (no scaling). The effective-scale variant goes in a separate field.
    bg_on_sample_raw = bg_on_sample
    bg_on_sample_effective = bg_in_eff

    note = (f"F-451 scale_direction={scale_direction}, applied={applied_scale:.4f}, "
            f"t_s={t_s:.0f}s, t_bg={t_bg:.0f}s, effective_live_time={effective_live_time:.0f}s. "
            f"scale_factor(legacy t_s/t_bg)={scale_factor_legacy:.3f}. "
            f"Покрытие фоновой сетки: {overlap*100:.1f}% каналов образца.")
    if overlap < 0.95:
        note += f" ⚠ Неполное покрытие — net spectrum в краях недостоверен."

    if not np.isnan(zero_point_mismatch_keV) and zero_point_mismatch_keV > 5.0:
        note += f" Δa₀={zero_point_mismatch_keV:.2f} keV (см. F-243 в bg_subtract_dual_mode.py для strict-mode subtraction)."

    return BackgroundSubtractionResult(
        net_counts=net,
        net_cps=net / effective_live_time,
        background_counts_on_sample_grid=bg_on_sample_raw,
        scale_factor=scale_factor_legacy,
        overlap_fraction=overlap,
        notes=note,
        net_uncertainties=net_uncertainties,
        gain_mismatch_relative=gain_mismatch_relative,
        zero_point_mismatch_keV=zero_point_mismatch_keV,
        n_channels_clipped=n_channels_clipped,
        applied_scale=applied_scale,
        scale_direction=scale_direction,
        effective_live_time=effective_live_time,
        background_counts_on_sample_grid_effective=bg_on_sample_effective,
    )


__all__ = [
    "BackgroundSubtractionResult",
    "subtract_background",
]