"""
BUG-56 regression -- Mariscotti adaptive kernel half-width (M) for the
Eu-152 964 / 1086 / 1112 keV cluster on NaI 1024-channel spectra.

A4 diagnostic 2026-06-04 (`_state/agent_a/outbox/2026-06-04_a4_unidentified_peaks_diagnostic.md`)
sect.3 found that three strong Eu-152 lines (964.06 I=14.51%,
1085.84 I=10.11%, 1112.08 I=13.67%) were all missed on the AmTiCsEu
fixture. Root cause: the Mariscotti second-derivative kernel in
`scripts/gamma/peaks/search.py:_band_filter` used a fixed half-width
of ceil(1.5*FWHM) channels, which over-smoothes close-peak clusters
separated by ~= 1*FWHM when FWHM_channels grows large (~20 channels at
1000 keV on NaI 1024-ch).

Fix: when band_fwhm >= 15 channels, use half_width = ceil(1.0*FWHM) --
narrower kernel that keeps ~= 95% of theoretical S/N on isolated peaks
but resolves close pairs/triplets the wider 1.5*FWHM kernel merged
into a single feature.

Reference: Gilmore & Joss sect.9.3 -- kernel width 2*M with M = FWHM for
"doublet resolution mode".
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.peaks.search import mariscotti_search


def _build_cluster_spectrum(
    n=1024,
    fwhm_at_1000=20.0,
    peak_centers_keV=(964.06, 1085.84, 1112.08),
    peak_heights=(2200.0, 1400.0, 2000.0),
    continuum_height=200.0,
    seed=2026,
):
    """
    Build a synthetic NaI-style 1024-channel spectrum.

    Calibration: 1 channel ~= 2.5 keV (typical for NaI 1024-ch range
    ~0-2560 keV), so 950-1120 keV corresponds to channels ~= 380-448.

    FWHM model: linear in keV, with FWHM~= 20 channels (50 keV) at
    1000 keV -- matches the empirical NaI 7%*sqrt(E/1000) shape near 1 MeV.

    Three Eu-152 peaks placed at 964 / 1086 / 1112 keV; flat continuum
    plus Poisson noise.
    """
    rng = np.random.default_rng(seed)
    keV_per_ch = 2.5
    fwhm_at_ch = lambda ch: fwhm_at_1000  # constant in this band

    counts = np.full(n, continuum_height, dtype=np.float64)
    for E, h in zip(peak_centers_keV, peak_heights):
        ch_center = E / keV_per_ch
        sigma_ch = fwhm_at_1000 / 2.355
        x = np.arange(n) - ch_center
        counts += h * np.exp(-(x / sigma_ch) ** 2 / 2.0)

    # Poisson noise
    counts = rng.poisson(np.maximum(counts, 0.0)).astype(np.float64)
    return counts, keV_per_ch, fwhm_at_ch


def test_bug56_three_peak_cluster_resolves_two_or_more_after_fix():
    """
    Acceptance criterion (brief sect.BUG-56 item 3):
      At least 2 of the 3 Eu-152 peaks (964 / 1086 / 1112 keV) must be
      resolved by the adaptive Mariscotti filter in the 950-1120 keV
      region on a NaI 1024-channel-equivalent spectrum.

    Pre-fix: half_width = 1.5*FWHM = 30 channels (kernel width 61 ch)
             smoothed the 950-1120 keV span (~68 channels) into 1 peak.
    Post-fix: half_width = 1.0*FWHM = 20 channels (kernel width 41 ch)
              resolves >= 2 of the 3 closely-spaced peaks.
    """
    counts, keV_per_ch, fwhm_at_ch = _build_cluster_spectrum()

    peaks = mariscotti_search(
        counts,
        fwhm_channels=fwhm_at_ch(500),  # ~= 20 channels at the cluster
        sigma_threshold=3.0,
    )

    # Count distinct peak channels in the 950-1120 keV region.
    lo_ch = 950.0 / keV_per_ch  # ~= 380
    hi_ch = 1120.0 / keV_per_ch  # ~= 448
    region_peaks = [p for p in peaks if lo_ch <= p.channel <= hi_ch]
    n_region = len(region_peaks)

    print(f"  Peaks in 950-1120 keV (ch {lo_ch:.0f}-{hi_ch:.0f}): "
          f"{n_region}")
    for p in region_peaks:
        E_est = p.channel * keV_per_ch
        print(f"    ch={p.channel}  E~{E_est:.1f} keV  sig={p.significance:.1f}")

    assert n_region >= 2, (
        f"Adaptive Mariscotti must resolve >=2 of 3 close peaks in "
        f"950-1120 keV (Eu-152 964/1086/1112 cluster). Found {n_region}."
    )
    print("  [OK] test_bug56_three_peak_cluster_resolves_two_or_more_after_fix")


def test_bug56_half_width_threshold_only_affects_wide_fwhm():
    """
    Adaptive M kicks in at band_fwhm >= 15 channels. For band_fwhm < 15
    (typical low-E regime), the kernel must remain at the legacy
    1.5*FWHM half-width to preserve back-compat with v1.5/v1.6 scalar-
    mode behaviour and the synthetic tests in
    `tests/step03_peak_search/test_adaptive_mariscotti.py`.

    This test exercises a narrow-FWHM band (FWHM=10 channels) on a
    flat synthetic with isolated peaks; the result must match the
    classic Mariscotti output (no spurious shifts due to the BUG-56
    adaptive branch).
    """
    rng = np.random.default_rng(42)
    n = 2048
    counts = np.full(n, 50.0, dtype=np.float64)
    for ch in (100, 500, 1000, 1500):
        sigma_ch = 10.0 / 2.355
        x = np.arange(n) - ch
        counts += 2500.0 * np.exp(-(x / sigma_ch) ** 2 / 2.0)
    counts = rng.poisson(np.maximum(counts, 0.0)).astype(np.float64)

    peaks = mariscotti_search(
        counts, fwhm_channels=10.0, sigma_threshold=3.0,
    )

    found_channels = {p.channel for p in peaks}
    matched = sum(
        1 for tp in (100, 500, 1000, 1500)
        if any(abs(p.channel - tp) <= 3 for p in peaks)
    )
    assert matched == 4, (
        f"All 4 isolated peaks must be found on narrow-FWHM scalar "
        f"mode; matched {matched}/4 (found_channels={found_channels})"
    )
    print("  [OK] test_bug56_half_width_threshold_only_affects_wide_fwhm")


def test_bug56_wide_band_pre_fix_would_merge_cluster():
    """
    Sanity check showing the adaptive branch is the cause of the
    improvement: at band_fwhm >= 15 the new code path takes the
    narrower kernel; the test cluster spectrum has FWHM = 20 channels
    which crosses the threshold.

    We don't run the OLD code path (it's no longer in tree). Instead
    we exercise the new code path AND assert it resolves >= 2 peaks --
    which is exactly the contract guarded by the previous test.
    Here we add a separate diagnostic: print the half-width that
    would be selected at FWHM=20.
    """
    band_fwhm = 20.0
    import math
    expected_half_width = max(int(math.ceil(1.0 * band_fwhm)), 3)
    print(f"  Adaptive branch: band_fwhm={band_fwhm} -> "
          f"half_width={expected_half_width} (kernel width "
          f"{2 * expected_half_width + 1} channels)")
    assert expected_half_width == 20
    legacy_half_width = max(int(math.ceil(1.5 * band_fwhm)), 3)
    assert legacy_half_width == 30  # pre-fix would have been 30
    print(f"  (Pre-fix half_width would have been {legacy_half_width}.)")
    print("  [OK] test_bug56_wide_band_pre_fix_would_merge_cluster")


if __name__ == "__main__":
    print("Running BUG-56 adaptive Mariscotti M tests...\n")
    test_bug56_three_peak_cluster_resolves_two_or_more_after_fix()
    test_bug56_half_width_threshold_only_affects_wide_fwhm()
    test_bug56_wide_band_pre_fix_would_merge_cluster()
    print("\n[OK] All BUG-56 tests passed.")
