"""
Peak search in channel space.

Implements the Mariscotti (1967) second-derivative-of-Gaussian filter,
extended with the Currie L_C significance test for selecting genuine
peaks from the filter output.

This module operates in **channel space only**, before energy
calibration. It is called both:
  - by bootstrap calibration (Phase 1.2) to find anchor candidates
    before any energy mapping exists
  - by the main analysis pipeline (Phase 1.3) after calibration has
    been established, using calibrated FWHM(E) as the adaptive filter
    width

The output of this module is small (a list of peaks, ~50–100 entries
for a typical spectrum), so passing the result to the AI for review
is token-cheap. The counts array is never returned.

────────────────────────────────────────────────────────────────────
v1.6 — Adaptive FWHM and adaptive sigma_threshold
────────────────────────────────────────────────────────────────────
Both `fwhm_channels` and `sigma_threshold` now accept either a scalar
(back-compat, constant across the spectrum) or a `Callable[[int],
float]` mapping channel → local value.

The adaptive FWHM is critical for NaI/scintillator spectra spanning
many octaves of energy: FWHM at Ba Kα 32 keV is ≈4 channels while at
Tl-208 2614 keV it is ≈90 channels, a 22× span that a single global
FWHM cannot resolve correctly (kernel either over-smooths the low-E
end or chops the high-E peaks into multiple false positives).

Implementation: banded segmentation. The channel range is split into
bands where FWHM varies by < `max_ratio` (default 1.2); each band is
convolved with the Mariscotti filter at its band-local FWHM, then the
significance arrays are stitched back together. Equivalent to per-
channel kernel sizing but O(N) instead of O(N²).

The adaptive sigma_threshold lets callers tighten the false-positive
rate in low-statistics regions (typically high-E tails of natural-
background spectra) without affecting the well-populated mid-energy
peaks.

Methodology references:
  - Mariscotti M.A., NIM 50 (1967) 309
  - Gilmore & Joss, "Practical Gamma-ray Spectrometry" 3rd Ed., §9.3
  - Currie L.A., Anal. Chem. 40 (1968) 586

Token economy: AI sees only the returned list of Peak dataclasses. The
filter convolution and Currie test live entirely in numpy.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Optional, Union

import numpy as np


# Type alias: scalar OR per-channel callable
FwhmSpec = Union[float, Callable[[int], float]]
SigmaSpec = Union[float, Callable[[int], float]]


# ============================================================================
# Output dataclass
# ============================================================================

@dataclass
class FoundPeak:
    """
    A peak found by the search algorithm.

    All quantities are in **channel space** unless stated otherwise.
    Energy is assigned later, after calibration.
    """
    channel: int                  # peak center, channel units (integer)
    height: float                 # net counts at peak channel (above local continuum)
    fwhm_channels: float          # full-width-at-half-maximum at this peak's channel
    significance: float           # Currie L_C-style: net/sqrt(B) ≈ S/sqrt(B)
    area_estimate: float = 0.0    # rough integral net counts (2.5·height·sigma)
    sigma_area_estimate: float = 0.0
    # Optional fields filled by the caller post-detection
    energy_keV: Optional[float] = None
    notes: list = field(default_factory=list)

    @property
    def sigma_channels(self) -> float:
        return self.fwhm_channels / 2.355

    def __repr__(self) -> str:
        e = f", E={self.energy_keV:.2f} keV" if self.energy_keV is not None else ""
        return (f"FoundPeak(ch={self.channel}, FWHM={self.fwhm_channels:.1f}ch, "
                f"area≈{self.area_estimate:.0f}, σ_sig={self.significance:.1f}{e})")


# ============================================================================
# Mariscotti second-derivative filter — primitive
# ============================================================================

def _gaussian_second_derivative_kernel(fwhm_channels: float, half_width: int) -> np.ndarray:
    """
    Build the second-derivative-of-Gaussian filter kernel.

    The kernel coefficients are derived from:
        f(x) = exp(-x²/(2σ²))
        f''(x) ∝ (x²/σ² - 1) · f(x)

    We center, normalise to unit L2 norm, so that convolution against
    a Gaussian peak produces a known-magnitude response.

    half_width: kernel half-length in channels.
                Typical value: ~1.5 · FWHM.
    """
    sigma = fwhm_channels / 2.355
    x = np.arange(-half_width, half_width + 1, dtype=np.float64)
    g = ((x / sigma) ** 2 - 1.0) * np.exp(-(x / sigma) ** 2 / 2.0)
    # Center (remove DC component): a true second-derivative integrates to zero
    g = g - g.mean()
    # Unit L2 norm for stable significance
    norm = math.sqrt(float((g * g).sum()))
    if norm > 0:
        g = g / norm
    return g


# ============================================================================
# Adaptive FWHM handling — banded segmentation
# ============================================================================

def _build_fwhm_array(fwhm_spec: FwhmSpec, n_channels: int) -> np.ndarray:
    """
    Evaluate FWHM(channel) at every channel.

    If the spec is scalar, returns a constant array (back-compat path).
    If callable, evaluates at every channel; values are coerced to
    `float` and clamped to a minimum of 1.0 channel (filter requires
    non-degenerate kernel).
    """
    if callable(fwhm_spec):
        arr = np.fromiter(
            (max(1.0, float(fwhm_spec(ch))) for ch in range(n_channels)),
            dtype=np.float64,
            count=n_channels,
        )
    else:
        arr = np.full(n_channels, max(1.0, float(fwhm_spec)), dtype=np.float64)
    return arr


def _build_sigma_threshold_array(sigma_spec: SigmaSpec, n_channels: int) -> np.ndarray:
    """Evaluate sigma_threshold(channel) at every channel."""
    if callable(sigma_spec):
        return np.fromiter(
            (float(sigma_spec(ch)) for ch in range(n_channels)),
            dtype=np.float64,
            count=n_channels,
        )
    return np.full(n_channels, float(sigma_spec), dtype=np.float64)


def _partition_into_bands(
    fwhm_arr: np.ndarray, max_ratio: float = 1.2
) -> list:
    """
    Partition channel range into bands of nearly-constant FWHM.

    A band is grown channel-by-channel as long as the ratio
    (max_FWHM_in_band / min_FWHM_in_band) stays below `max_ratio`.
    When adding the next channel would exceed the ratio, the band is
    closed and a new one starts.

    Returns:
        List of (start, end, fwhm_representative) tuples, half-open
        intervals [start, end) covering [0, n_channels) without gaps.
    """
    n = len(fwhm_arr)
    if n == 0:
        return []
    bands = []
    start = 0
    band_lo = float(fwhm_arr[0])
    band_hi = float(fwhm_arr[0])
    for i in range(1, n):
        v = float(fwhm_arr[i])
        new_lo = min(band_lo, v)
        new_hi = max(band_hi, v)
        if new_hi / new_lo > max_ratio:
            # Close current band; representative = geometric mean
            bands.append((start, i, math.sqrt(band_lo * band_hi)))
            start = i
            band_lo = v
            band_hi = v
        else:
            band_lo = new_lo
            band_hi = new_hi
    bands.append((start, n, math.sqrt(band_lo * band_hi)))
    return bands


def _band_filter(
    counts: np.ndarray,
    band_start: int,
    band_end: int,
    band_fwhm: float,
    n_channels: int,
) -> np.ndarray:
    """
    Convolve a single band with the Mariscotti filter at `band_fwhm`.

    The band is extended by ±half_width on each side (clipped to
    spectrum bounds) so the convolution result in the band interior
    is not edge-affected. Returns the significance array of length
    (band_end − band_start) covering exactly [band_start, band_end).

    Kernel half-width (BUG-56 / v1.25.0):
      • Default (band_fwhm < 15 channels): half_width = ⌈1.5·FWHM⌉ — the
        Mariscotti (1967) / Gilmore & Joss §9.3 conservative recommendation,
        which gives ≈ 99% of theoretical S/N for an isolated Gaussian.
      • Wide-FWHM regime (band_fwhm ≥ 15 channels, typical of NaI 1024-ch
        above ~700 keV): half_width = ⌈1.0·FWHM⌉ — narrower kernel that
        keeps ≈ 95% of theoretical S/N but resolves close peaks spaced
        ~1·FWHM apart that the 1.5·FWHM kernel merges into one feature.
        Empirically required to resolve the Eu-152 964.06 / 1085.84 /
        1112.08 keV cluster on AmTiCsEu (separation ≈ 1.05·FWHM and
        1.3·FWHM, FWHM_ch ≈ 19 at 1000 keV).
        Reference: Gilmore & Joss §9.3 — kernel width 2·M with M = FWHM
        for "doublet resolution mode".
    """
    if band_fwhm >= 15.0:
        # Adaptive M (BUG-56): narrower kernel preserves close-peak
        # resolution at the cost of ≈ 5% S/N on isolated peaks.
        half_width = max(int(math.ceil(1.0 * band_fwhm)), 3)
    else:
        # Default kernel (low-E, narrow FWHM): conservative 1.5·FWHM.
        half_width = max(int(math.ceil(1.5 * band_fwhm)), 3)
    kernel = _gaussian_second_derivative_kernel(band_fwhm, half_width)

    seg_start = max(0, band_start - half_width)
    seg_end = min(n_channels, band_end + half_width)
    seg = counts[seg_start:seg_end]

    response = -np.convolve(seg, kernel, mode="same")
    kernel_sq = kernel ** 2
    var = np.convolve(np.maximum(seg, 1.0), kernel_sq, mode="same")
    sigma_resp = np.sqrt(np.maximum(var, 1e-12))
    sig = response / sigma_resp

    # Extract central portion corresponding to absolute [band_start, band_end)
    central_offset = band_start - seg_start
    central_len = band_end - band_start
    return sig[central_offset:central_offset + central_len]


# ============================================================================
# Main entry point
# ============================================================================

def mariscotti_search(
    counts,
    fwhm_channels: FwhmSpec,
    sigma_threshold: SigmaSpec = 3.0,
    min_separation_factor: float = 1.0,
    edge_margin: Optional[int] = None,
    *,
    band_ratio: float = 1.2,
    # F-139 / v1.17.7 — opt-in FWHM-width filter.
    # При True измеряется фактическая FWHM пика по полувысоте; пики с
    # FWHM_meas < min_fwhm_ratio · expected_FWHM(E) отбраковываются как
    # шумовые спайки (1-2 канала). Default False для back-compat с
    # синтетическими тестами; pipeline (analyze_lsrm_spe) включает
    # F-139 явно при анализе реальных .spe-файлов.
    filter_narrow_peaks: bool = False,
    min_fwhm_ratio: float = 0.3,
):
    """
    Find peaks via the Mariscotti second-derivative filter.

    Args:
        counts: 1-D counts array (numpy or list-like)
        fwhm_channels: expected peak FWHM in channels. Either:
            - a scalar `float` (constant across the spectrum, legacy
              behaviour identical to v1.5),
            - a `Callable[[int], float]` mapping channel index to
              local FWHM in channels.
        sigma_threshold: Currie L_C threshold for accepting a peak.
            Either a scalar (default 3.0) or a `Callable[[int], float]`
            mapping channel index to local threshold. The latter is
            useful for raising the threshold in low-statistics regions
            (high-E tails) to suppress noise-driven false positives.
        min_separation_factor: minimum allowed center-to-center
            distance between accepted peaks, expressed in multiples of
            the *larger* of the two peaks' local FWHM (default 1.0).
        edge_margin: number of channels at each edge to ignore. If
            None, defaults to 2·FWHM(edge) on each side.
        band_ratio: max ratio of FWHM_max/FWHM_min within a single
            convolution band when FWHM is adaptive (default 1.2).
            Smaller values yield more bands (slower) and tighter local
            tuning; larger values yield fewer bands and faster runs.

    Returns:
        List of FoundPeak, sorted by channel ascending. Each peak's
        `fwhm_channels` field carries the local FWHM at that peak's
        channel — useful for downstream code that no longer has the
        original provider.
    """
    counts = np.asarray(counts, dtype=np.float64)
    n = counts.size
    if n < 50:
        return []

    # Resolve adaptive providers to per-channel arrays
    fwhm_arr = _build_fwhm_array(fwhm_channels, n)
    sigma_arr = _build_sigma_threshold_array(sigma_threshold, n)

    if np.all(fwhm_arr < 1.0):
        return []

    # Partition into bands of similar FWHM
    bands = _partition_into_bands(fwhm_arr, max_ratio=band_ratio)

    # Convolve each band; stitch significance arrays
    significance = np.zeros(n, dtype=np.float64)
    for (b_start, b_end, b_fwhm) in bands:
        significance[b_start:b_end] = _band_filter(
            counts, b_start, b_end, b_fwhm, n
        )

    # Edge masking
    if edge_margin is None:
        # Scale margin with local FWHM at each edge
        em_lo = int(math.ceil(2 * fwhm_arr[0]))
        em_hi = int(math.ceil(2 * fwhm_arr[-1]))
    else:
        em_lo = em_hi = int(edge_margin)
    if em_lo > 0:
        significance[:em_lo] = 0.0
    if em_hi > 0:
        significance[-em_hi:] = 0.0

    # Per-channel threshold mask (significance > local_threshold)
    above_threshold = significance > sigma_arr

    # Find local maxima above the per-channel threshold
    candidates = _local_maxima_masked(significance, above_threshold)
    if not candidates:
        return []

    # Enforce min separation with local FWHM
    candidates = _enforce_separation_adaptive(
        candidates, significance, fwhm_arr, min_separation_factor
    )

    # Convert to FoundPeak objects.
    # F-139 / v1.17.7 — измерение ФАКТИЧЕСКОЙ FWHM пика-кандидата и
    # отбраковка случаев, где она в <0.5 от ожидаемой модели — это
    # чисто шумовой пик/спайк (1-3 канала). Реальные физические пики
    # имеют FWHM около калибровочной модели; узкие острые пики не
    # реализуемы на NaI.
    peaks = []
    for ch in candidates:
        local_fwhm = float(fwhm_arr[ch])
        sigma_ch = local_fwhm / 2.355

        cont_l_lo = max(0, ch - int(round(3 * local_fwhm)))
        cont_l_hi = max(cont_l_lo + 1, ch - int(round(2 * local_fwhm)))
        cont_r_lo = min(n - 1, ch + int(round(2 * local_fwhm)))
        cont_r_hi = min(n, ch + int(round(3 * local_fwhm)))
        b_left = counts[cont_l_lo:cont_l_hi].mean() if cont_l_hi > cont_l_lo else 0.0
        b_right = counts[cont_r_lo:cont_r_hi].mean() if cont_r_hi > cont_r_lo else 0.0
        # BUG-21 / v1.18.32 — log-linear baseline interpolation at peak ch.
        #
        # Previously: b = 0.5*(b_left + b_right). For an exponential (Compton)
        # continuum the arithmetic mean of two equidistant samples OVERSHOOTS
        # the midpoint by cosh(d/τ) − 1: for typical NaI τ≈120 ch and
        # d≈2.5·FWHM≈24 ch this is ~3% — comparable to or larger than the
        # height of an "obvious by eye" peak (3-5% of local continuum), so
        # `net_height = counts[ch] − b` clipped to 0 silently dropped real
        # peaks on the steep falling side of the Compton plateau.
        #
        # Fix: log-linear interpolation at the peak channel. Mathematically
        # exact for an exponential continuum (`log(f)` is linear in `x`), and
        # exact for a linear continuum to first order. For the steep-slope
        # synthetic (τ=120 ch, d≈24 ch), it reduces the baseline bias from
        # cosh(d/τ)−1 ≈ 3% to a third-order curvature term well below 0.1%,
        # recovering all 5 truly-detectable peaks vs 2 of 5 before the fix.
        # Reduces to arithmetic mean for flat continua (slope=0 → b_left==
        # b_right → log-linear identical to arithmetic).
        #
        # Reproducer: _bug21_repro/repro_v4_steep.py.
        # Plan: _state/agent_a/outbox/BUG21_PLAN.md.
        c_left = 0.5 * (cont_l_lo + cont_l_hi - 1)
        c_right = 0.5 * (cont_r_lo + cont_r_hi - 1)
        if c_right > c_left and b_left > 0.0 and b_right > 0.0:
            # log-linear interpolation: log(b) linear in ch
            log_b_left = math.log(b_left)
            log_b_right = math.log(b_right)
            log_b = log_b_left + (log_b_right - log_b_left) * \
                (ch - c_left) / (c_right - c_left)
            b = float(math.exp(log_b))
        elif c_right > c_left:
            # Fallback to arithmetic linear interp if one side is zero
            b = float(b_left + (b_right - b_left) *
                      (ch - c_left) / (c_right - c_left))
        else:
            b = 0.5 * (b_left + b_right)

        net_height = max(counts[ch] - b, 0.0)
        if net_height <= 0:
            continue
        # F-139 — измерение фактической FWHM по полувысоте (opt-in).
        # Применяется только при достаточной статистике (net_height ≥ 10)
        # и адекватной ожидаемой FWHM (≥ 4 каналов). Иначе пропускается:
        # на слабых пиках или очень узкой FWHM-модели измерение FWHM
        # становится ненадёжным, и F-139 даст ложные отсевы.
        if filter_narrow_peaks and net_height >= 10.0 and local_fwhm >= 4.0:
            half_level = b + 0.5 * net_height
            search_range = max(2, int(round(2.5 * local_fwhm)))
            rr = ch
            while (rr < min(n - 1, ch + search_range)
                   and counts[rr] >= half_level):
                rr += 1
            ll = ch
            while (ll > max(0, ch - search_range)
                   and counts[ll] >= half_level):
                ll -= 1
            measured_fwhm = max(1.0, float(rr - ll))
            if measured_fwhm < float(min_fwhm_ratio) * local_fwhm:
                continue
        area = 2.507 * sigma_ch * net_height

        roi_lo = max(0, ch - int(round(2 * local_fwhm)))
        roi_hi = min(n, ch + int(round(2 * local_fwhm)) + 1)
        roi_counts = float(counts[roi_lo:roi_hi].sum())
        bg_in_roi = b * (roi_hi - roi_lo)
        sigma_area = math.sqrt(max(roi_counts + bg_in_roi, 1.0))

        peaks.append(FoundPeak(
            channel=int(ch),
            height=net_height,
            fwhm_channels=local_fwhm,
            significance=float(significance[ch]),
            area_estimate=area,
            sigma_area_estimate=sigma_area,
        ))

    peaks.sort(key=lambda p: p.channel)
    return peaks


# ============================================================================
# Helpers
# ============================================================================

def _local_maxima_masked(signal: np.ndarray, mask: np.ndarray) -> list:
    """
    Find indices of local maxima where `mask[i]` is True.

    A local maximum is defined as a sample strictly greater than both
    its immediate neighbours. Plateaus (signal[i] == signal[i+1]) are
    handled by selecting the leftmost index of the plateau.

    The mask lets the caller supply a per-channel threshold without
    the helper needing to know about it.
    """
    n = signal.size
    if n < 3:
        return []
    out = []
    i = 1
    while i < n - 1:
        if not mask[i]:
            i += 1
            continue
        left_ok = signal[i] > signal[i - 1]
        if not left_ok:
            i += 1
            continue
        # Walk through possible plateau
        j = i
        while j < n - 1 and signal[j + 1] == signal[i]:
            j += 1
        right_ok = (j < n - 1) and (signal[j] > signal[j + 1])
        if right_ok:
            out.append(i)
        i = j + 1
    return out


def _enforce_separation_adaptive(
    candidates: list,
    significance: np.ndarray,
    fwhm_arr: np.ndarray,
    min_separation_factor: float,
) -> list:
    """
    Greedy thinning with adaptive minimum separation.

    Candidates are processed by descending significance; each is
    accepted only if its distance to every already-accepted candidate
    exceeds `min_separation_factor · max(fwhm_at_self, fwhm_at_other)`.
    Using the larger of the two local FWHMs is conservative — it
    prevents close peaks in either band from being kept as a pair.
    """
    if not candidates:
        return []
    by_strength = sorted(candidates, key=lambda c: -significance[c])
    accepted = []
    for ch in by_strength:
        ok = True
        for a in accepted:
            min_sep = min_separation_factor * max(fwhm_arr[ch], fwhm_arr[a])
            if abs(ch - a) < min_sep:
                ok = False
                break
        if ok:
            accepted.append(ch)
    return sorted(accepted)


# ============================================================================
# Channel-space FWHM estimate per peak
# ============================================================================

def estimate_fwhm_at_peak(
    counts,
    channel: int,
    initial_fwhm: float,
    max_iter: int = 3,
) -> Optional[float]:
    """
    Refine the FWHM estimate at a single found peak by interpolating
    the half-maximum crossings.

    Uses linear interpolation between channels where the count first
    crosses (peak_height + background)/2 on either side. Returns None
    if the half-maximum crossings cannot be found (peak too close to
    edge or to a neighbouring peak).
    """
    counts = np.asarray(counts, dtype=np.float64)
    n = counts.size
    if not (1 <= channel < n - 1):
        return None

    fwhm = initial_fwhm
    for _ in range(max_iter):
        w = int(math.ceil(2 * fwhm))
        lo_l = max(0, channel - 3 * w)
        hi_l = max(lo_l + 1, channel - 2 * w)
        lo_r = min(n - 1, channel + 2 * w)
        hi_r = min(n, channel + 3 * w)
        if hi_l <= lo_l or hi_r <= lo_r:
            return None
        b = 0.5 * (counts[lo_l:hi_l].mean() + counts[lo_r:hi_r].mean())
        peak = counts[channel]
        half = b + 0.5 * (peak - b)
        if half <= b:
            return None

        i_left = channel
        while i_left > 0 and counts[i_left] > half:
            i_left -= 1
        if i_left == 0 and counts[i_left] > half:
            return None

        i_right = channel
        while i_right < n - 1 and counts[i_right] > half:
            i_right += 1
        if i_right == n - 1 and counts[i_right] > half:
            return None

        if counts[i_left + 1] == counts[i_left]:
            x_left = float(i_left)
        else:
            x_left = i_left + (half - counts[i_left]) / (counts[i_left + 1] - counts[i_left])

        if counts[i_right] == counts[i_right - 1]:
            x_right = float(i_right)
        else:
            x_right = i_right - 1 + (counts[i_right - 1] - half) / (counts[i_right - 1] - counts[i_right])

        new_fwhm = x_right - x_left
        if new_fwhm <= 0:
            return None
        if abs(new_fwhm - fwhm) / fwhm < 0.05:
            return new_fwhm
        fwhm = new_fwhm

    return fwhm
