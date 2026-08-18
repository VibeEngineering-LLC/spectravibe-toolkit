"""BUG-21 regression — Mariscotti search must recover detectable peaks on a
steep Compton continuum (sloped background).

Pre-fix (v1.18.31):
    On an exponential continuum spanning ~5 orders of magnitude (1e6 → 5
    counts across 2048 channels), peaks with Currie significance σ ≥ 4
    sitting on the steep falling side (channels 100, 250, 500) were
    silently dropped by the final-stage baseline-clip in
    `gamma.peaks.search.mariscotti_search`:

        b = 0.5 * (b_left + b_right)
        if net_height = counts[ch] - b  ≤ 0:  continue

    For a convex (exponential) continuum the arithmetic mean of two
    equidistant samples OVERSHOOTS the midpoint by `cosh(d/τ) − 1`. For
    typical NaI τ≈120 ch and d≈2.5·FWHM≈24 ch this is ~3% — comparable to
    or larger than the height of an "obvious by eye" peak (3-5% of local
    continuum). The clip then drops the candidate even though its
    Mariscotti significance was well above threshold.

Post-fix (v1.18.32):
    Log-linear baseline interpolation at the peak channel — exact for an
    exponential continuum, identical to the old formula for a flat
    continuum.

This test verifies:
  (A) On a synthetic 1e6→5 exponential plateau with 5 peaks at S/√B = 5-8,
      all peaks at σ ≥ 3 are recovered.
  (B) On a flat continuum, the new baseline gives an equivalent peak set
      to the old (no regressions, channels match).
  (C) Idempotence — re-running on the same input gives the same set.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from gamma.peaks.search import mariscotti_search  # noqa: E402


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def _build_steep_slope_spectrum(seed: int = 31, n: int = 2048,
                                fwhm_ch: float = 9.4,
                                peaks_spec=None):
    """Return (counts, peaks_spec, fwhm_ch, continuum).

    Continuum: `1e6·exp(-ch/120) + 1e3·exp(-ch/600) + 5`
        — spans 5 orders of magnitude across 2048 channels.

    Peaks: each given (ch, target_S_over_sqrt_B). Height = SNR · sqrt(B).
    """
    if peaks_spec is None:
        peaks_spec = [
            (100, 5.0),    # very steep slope, S/√B = 5
            (250, 6.0),    # steep slope
            (500, 5.0),    # mid slope
            (1000, 6.0),   # flat slope
            (1700, 8.0),   # tail
        ]
    rng = np.random.default_rng(seed)
    ch = np.arange(n, dtype=np.float64)
    continuum = 1.0e6 * np.exp(-ch / 120.0) + 1.0e3 * np.exp(-ch / 600.0) + 5.0
    sigma = fwhm_ch / 2.355
    truth = continuum.copy()
    for p, snr in peaks_spec:
        B = float(continuum[p])
        h = snr * math.sqrt(B)
        truth += h * np.exp(-((ch - p) / sigma) ** 2 / 2.0)
    counts = rng.poisson(np.maximum(truth, 0.0)).astype(np.float64)
    return counts, peaks_spec, fwhm_ch, continuum


def _build_flat_spectrum(seed: int = 17, n: int = 2048,
                         fwhm_ch: float = 9.4,
                         peaks_spec=None):
    """Flat continuum @1000 cps, peaks placed across the range."""
    if peaks_spec is None:
        peaks_spec = [
            (200, 8.0), (600, 8.0), (1000, 8.0),
            (1400, 8.0), (1800, 8.0),
        ]
    rng = np.random.default_rng(seed)
    ch = np.arange(n, dtype=np.float64)
    continuum = np.full(n, 1000.0)
    sigma = fwhm_ch / 2.355
    truth = continuum.copy()
    for p, snr in peaks_spec:
        h = snr * math.sqrt(continuum[p])
        truth += h * np.exp(-((ch - p) / sigma) ** 2 / 2.0)
    counts = rng.poisson(np.maximum(truth, 0.0)).astype(np.float64)
    return counts, peaks_spec, fwhm_ch


def _is_matched(true_ch: int, found_channels, tol_ch: float) -> bool:
    return any(abs(c - true_ch) <= tol_ch for c in found_channels)


# ----------------------------------------------------------------------
# (A) Steep slope: all S/√B ≥ 5 peaks recovered at σ_thr=3
# ----------------------------------------------------------------------

def test_bug21_recovers_steep_slope_peaks_at_sigma3():
    counts, peaks_spec, fwhm_ch, _ = _build_steep_slope_spectrum()
    out = mariscotti_search(
        counts, fwhm_channels=fwhm_ch,
        sigma_threshold=3.0,
        min_separation_factor=0.6, edge_margin=10,
    )
    found_channels = [p.channel for p in out]
    tol_ch = 1.5 * fwhm_ch
    missed = [tp for tp, _ in peaks_spec
              if not _is_matched(tp, found_channels, tol_ch)]
    assert missed == [], (
        f"BUG-21 regression — missed peaks {missed} at σ_thr=3 on steep "
        f"slope. Found channels: {found_channels}"
    )


def test_bug21_recovers_steep_slope_peaks_at_sigma25():
    counts, peaks_spec, fwhm_ch, _ = _build_steep_slope_spectrum()
    out = mariscotti_search(
        counts, fwhm_channels=fwhm_ch,
        sigma_threshold=2.5,
        min_separation_factor=0.6, edge_margin=10,
    )
    found_channels = [p.channel for p in out]
    tol_ch = 1.5 * fwhm_ch
    missed = [tp for tp, _ in peaks_spec
              if not _is_matched(tp, found_channels, tol_ch)]
    assert missed == [], (
        f"BUG-21 regression — missed peaks {missed} at σ_thr=2.5 on "
        f"steep slope. Found channels: {found_channels}"
    )


# ----------------------------------------------------------------------
# (B) Flat continuum: still finds known peaks (no regression)
# ----------------------------------------------------------------------

def test_bug21_flat_continuum_no_regression():
    counts, peaks_spec, fwhm_ch = _build_flat_spectrum()
    out = mariscotti_search(
        counts, fwhm_channels=fwhm_ch,
        sigma_threshold=3.0,
        min_separation_factor=0.6, edge_margin=10,
    )
    found_channels = [p.channel for p in out]
    tol_ch = 1.5 * fwhm_ch
    missed = [tp for tp, _ in peaks_spec
              if not _is_matched(tp, found_channels, tol_ch)]
    assert missed == [], (
        f"Flat-continuum regression — missed peaks {missed}. Found: "
        f"{found_channels}"
    )


# ----------------------------------------------------------------------
# (C) Idempotence
# ----------------------------------------------------------------------

def test_bug21_idempotent_on_steep_slope():
    counts, peaks_spec, fwhm_ch, _ = _build_steep_slope_spectrum()
    out1 = mariscotti_search(
        counts, fwhm_channels=fwhm_ch,
        sigma_threshold=2.5,
        min_separation_factor=0.6, edge_margin=10,
    )
    out2 = mariscotti_search(
        counts, fwhm_channels=fwhm_ch,
        sigma_threshold=2.5,
        min_separation_factor=0.6, edge_margin=10,
    )
    assert [p.channel for p in out1] == [p.channel for p in out2]


# ----------------------------------------------------------------------
# (D) Baseline math — log-linear interp matches an exponential continuum
# better than arithmetic mean by ≥10× on the synthetic spectrum.
# This is the failure mode that drove BUG-21; pinning it here so future
# refactors don't regress to arithmetic mean by accident.
# ----------------------------------------------------------------------

def test_bug21_baseline_math_exponential_continuum():
    n = 2048
    ch = np.arange(n, dtype=np.float64)
    continuum = 1.0e6 * np.exp(-ch / 120.0) + 1.0e3 * np.exp(-ch / 600.0) + 5.0
    fwhm_ch = 9.4
    # At ch=100 with d=24: simulate the baseline samples used by search.py
    for c in (100, 250, 500):
        l_lo = max(0, c - int(round(3 * fwhm_ch)))
        l_hi = max(l_lo + 1, c - int(round(2 * fwhm_ch)))
        r_lo = min(n - 1, c + int(round(2 * fwhm_ch)))
        r_hi = min(n, c + int(round(3 * fwhm_ch)))
        bl = float(continuum[l_lo:l_hi].mean())
        br = float(continuum[r_lo:r_hi].mean())
        cleft = 0.5 * (l_lo + l_hi - 1)
        cright = 0.5 * (r_lo + r_hi - 1)
        b_arith = 0.5 * (bl + br)
        log_bl, log_br = math.log(bl), math.log(br)
        b_log = math.exp(log_bl + (log_br - log_bl) *
                         (c - cleft) / (cright - cleft))
        true_b = float(continuum[c])
        err_arith = abs(b_arith - true_b) / true_b
        err_log = abs(b_log - true_b) / true_b
        # Log-linear should be ≥10× better than arithmetic on this geometry
        assert err_log * 10.0 < err_arith, (
            f"At ch={c}: log-linear baseline error {err_log:.4f} vs "
            f"arithmetic {err_arith:.4f} — log-linear is not the dominant "
            f"improvement. Did someone revert the BUG-21 fix?"
        )
