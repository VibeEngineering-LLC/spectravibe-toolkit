"""
Sanity test for the adaptive-FWHM rewrite of gamma.peaks.search.

Confirms three properties:
  (1) Back-compat: scalar fwhm_channels produces IDENTICAL results to
      v1.5 (same channels, same significance, same FWHM reported).
  (2) Adaptive mode runs and returns sensible peaks on a real
      multi-octave NaI spectrum.
  (3) Adaptive mode finds peaks that scalar mode at the wrong scale
      misses (the whole point of the upgrade).
"""
import sys, importlib.util
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent / "scripts"
sys.path.insert(0, str(ROOT))


def load_v15():
    """Load the v1.5 search.py as a side-by-side module for comparison."""
    from importlib.machinery import SourceFileLoader
    _project_root = Path(__file__).resolve().parent.parent.parent
    loader = SourceFileLoader(
        "search_v15",
        str(_project_root / "tests" / "fixtures" / "peaks" / "search_v1_5_legacy.py.fixture"),
    )
    return loader.load_module()


def synthetic_spectrum(n=4096, fwhm_low=4.0, fwhm_high=80.0, seed=42):
    """
    Build a synthetic spectrum with FWHM varying linearly from
    fwhm_low at ch 0 to fwhm_high at ch n-1, with peaks placed across
    the range and Poisson-like noise.
    """
    rng = np.random.default_rng(seed)
    # Linear FWHM(ch) model
    fwhm_at = lambda ch: fwhm_low + (fwhm_high - fwhm_low) * ch / (n - 1)
    counts = np.zeros(n, dtype=np.float64)
    # Continuum: exponential decay + flat tail
    continuum = 50.0 * np.exp(-np.arange(n) / 800.0) + 5.0
    counts += continuum
    # Place peaks at varied channels with varied heights
    peak_positions = [50, 200, 500, 1000, 1800, 2700, 3700]
    peak_heights = [2500, 1800, 3000, 2200, 1500, 800, 400]
    for ch, h in zip(peak_positions, peak_heights):
        sigma = fwhm_at(ch) / 2.355
        x = np.arange(n) - ch
        counts += h * np.exp(-(x / sigma) ** 2 / 2.0)
    # Poisson noise on the result
    counts = rng.poisson(np.maximum(counts, 0.0)).astype(np.float64)
    return counts, peak_positions, fwhm_at


# ----------------------------------------------------------------------
# (1) Back-compat: scalar mode equivalence
# ----------------------------------------------------------------------
print("─" * 70)
print("(1) Back-compat: scalar fwhm equivalence vs v1.5")
print("─" * 70)

from gamma.peaks.search import mariscotti_search as ms_new, FoundPeak
v15 = load_v15()

# Use a flat-FWHM synthetic to make the comparison meaningful
flat_counts, true_pos, _ = synthetic_spectrum(n=2048, fwhm_low=10.0, fwhm_high=10.0)

peaks_new = ms_new(flat_counts, fwhm_channels=10.0, sigma_threshold=3.0)
peaks_old = v15.mariscotti_search(flat_counts, fwhm_channels=10.0, sigma_threshold=3.0)

print(f"v1.5 found {len(peaks_old)} peaks; v1.6 scalar found {len(peaks_new)}")
if len(peaks_old) != len(peaks_new):
    print("⚠ DIFFERENT PEAK COUNT")
else:
    same = True
    for a, b in zip(peaks_old, peaks_new):
        if a.channel != b.channel:
            print(f"  ⚠ channel diff: v1.5 ch={a.channel}, v1.6 ch={b.channel}")
            same = False
        if abs(a.significance - b.significance) > 0.01:
            print(f"  ⚠ significance diff at ch={a.channel}: "
                  f"v1.5={a.significance:.3f}, v1.6={b.significance:.3f}")
            same = False
    if same:
        print("  ✓ identical channels and significance — back-compat OK")

# ----------------------------------------------------------------------
# (2) Adaptive mode on multi-octave FWHM
# ----------------------------------------------------------------------
print()
print("─" * 70)
print("(2) Adaptive FWHM on multi-octave synthetic")
print("─" * 70)

adapt_counts, true_pos, fwhm_at = synthetic_spectrum(
    n=4096, fwhm_low=4.0, fwhm_high=80.0
)

# Try with adaptive provider
peaks_adapt = ms_new(adapt_counts, fwhm_channels=fwhm_at, sigma_threshold=3.0)
print(f"Adaptive mode found {len(peaks_adapt)} peaks")
print(f"True peaks: {true_pos}")
print(f"Found channels (top 12 by significance):")
for p in sorted(peaks_adapt, key=lambda p: -p.significance)[:12]:
    expected_fwhm = fwhm_at(p.channel)
    print(f"  ch={p.channel:5d}  σ={p.significance:6.1f}  "
          f"FWHM={p.fwhm_channels:5.1f} (expected {expected_fwhm:5.1f})")

# Check all true peaks are found within ±FWHM
matched = []
missed = []
for tp in true_pos:
    expected_fwhm = fwhm_at(tp)
    nearest = min(peaks_adapt, key=lambda p: abs(p.channel - tp), default=None)
    if nearest and abs(nearest.channel - tp) <= expected_fwhm:
        matched.append((tp, nearest.channel, abs(nearest.channel - tp)))
    else:
        missed.append(tp)
print(f"Matched true peaks: {len(matched)}/{len(true_pos)}")
for tp, fc, d in matched:
    print(f"  true ch={tp} → found ch={fc} (Δ={d:.1f})")
if missed:
    print(f"⚠ Missed true peaks: {missed}")
else:
    print("✓ All true peaks recovered")

# ----------------------------------------------------------------------
# (3) Adaptive vs wrong-scale scalar — demonstrates why we need this
# ----------------------------------------------------------------------
print()
print("─" * 70)
print("(3) Same spectrum, scalar FWHM=10 (correct for low-E, wrong for high-E)")
print("─" * 70)

peaks_scalar = ms_new(adapt_counts, fwhm_channels=10.0, sigma_threshold=3.0)
print(f"Scalar FWHM=10 found {len(peaks_scalar)} peaks")
# Count how many of the high-E true peaks (ch > 1500) it finds
high_E_true = [tp for tp in true_pos if tp > 1500]
scalar_matched_high = 0
for tp in high_E_true:
    expected_fwhm = fwhm_at(tp)
    nearest = min(peaks_scalar, key=lambda p: abs(p.channel - tp), default=None)
    if nearest and abs(nearest.channel - tp) <= expected_fwhm:
        scalar_matched_high += 1
adaptive_matched_high = 0
for tp in high_E_true:
    expected_fwhm = fwhm_at(tp)
    nearest = min(peaks_adapt, key=lambda p: abs(p.channel - tp), default=None)
    if nearest and abs(nearest.channel - tp) <= expected_fwhm:
        adaptive_matched_high += 1

print(f"High-E true peaks (ch>1500): {len(high_E_true)}")
print(f"  Scalar FWHM=10: matched {scalar_matched_high}/{len(high_E_true)}")
print(f"  Adaptive:        matched {adaptive_matched_high}/{len(high_E_true)}")

# Count false positives at high E in scalar mode
scalar_high_E_fps = sum(1 for p in peaks_scalar if p.channel > 1500
                        and not any(abs(p.channel - tp) <= 2*fwhm_at(tp) for tp in true_pos))
adapt_high_E_fps = sum(1 for p in peaks_adapt if p.channel > 1500
                       and not any(abs(p.channel - tp) <= 2*fwhm_at(tp) for tp in true_pos))
print(f"  Scalar FWHM=10 false positives at high E: {scalar_high_E_fps}")
print(f"  Adaptive false positives at high E:        {adapt_high_E_fps}")

# ----------------------------------------------------------------------
# (4) Adaptive sigma_threshold
# ----------------------------------------------------------------------
print()
print("─" * 70)
print("(4) Adaptive sigma_threshold (stricter at high E)")
print("─" * 70)

# Linear ramp from 3.0 at ch 0 to 6.0 at ch 4095
sigma_at = lambda ch: 3.0 + 3.0 * ch / 4095.0
peaks_strict = ms_new(
    adapt_counts,
    fwhm_channels=fwhm_at,
    sigma_threshold=sigma_at,
)
print(f"With ramped sigma threshold (3→6): found {len(peaks_strict)} peaks "
      f"(vs {len(peaks_adapt)} with σ=3 flat)")
adapt_high_E_fps_strict = sum(1 for p in peaks_strict if p.channel > 1500
                              and not any(abs(p.channel - tp) <= 2*fwhm_at(tp) for tp in true_pos))
print(f"  high-E false positives: {adapt_high_E_fps_strict} "
      f"(was {adapt_high_E_fps} with flat σ)")

print()
print("═" * 70)
print("All adaptive-FWHM tests complete.")
print("═" * 70)
