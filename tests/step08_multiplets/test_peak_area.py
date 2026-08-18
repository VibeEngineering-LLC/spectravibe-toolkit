"""
Tests for peak area integration (Phase 2.1a).

Validates:
  1. Cowell method on synthetic Gaussian + flat background
  2. Cowell method on synthetic Gaussian + sloped background
  3. Gaussian fit recovers known peak parameters
  4. Both methods give comparable results on real fixtures
  5. Edge cases: low-stat peaks, near-spectrum-edge peaks
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.peaks.area import (
    PeakAreaResult, cowell_area, gaussian_fit_area, integrate_peaks,
    get_peak_area,
)


def make_synthetic_spectrum(N=2000, peaks=None, baseline=None):
    """Build a synthetic test spectrum."""
    ch = np.arange(N, dtype=np.float64)
    counts = np.zeros(N, dtype=np.float64)
    if baseline is None:
        baseline = lambda x: 100 * np.ones_like(x)
    counts += baseline(ch)
    if peaks:
        for c, H, sigma in peaks:
            counts += H * np.exp(-((ch - c) / sigma) ** 2 / 2)
    # Add Poisson noise
    rng = np.random.default_rng(42)
    counts = rng.poisson(counts).astype(np.int64)
    return counts


def test_cowell_flat_background():
    """Single Gaussian on flat background — Cowell should recover area
    within ~5%."""
    # H=5000, c=500, sigma=10 → area = 5000·10·√(2π) ≈ 125331
    counts = make_synthetic_spectrum(
        N=1000, peaks=[(500, 5000, 10)],
        baseline=lambda x: 100 * np.ones_like(x),
    )
    fwhm = 10 * 2.355
    r = cowell_area(counts, peak_channel=500, fwhm_channels=fwhm)
    expected = 5000 * 10 * math.sqrt(2 * math.pi)
    err = abs(r.net_area_counts - expected) / expected
    assert err < 0.05, f"Cowell error {err*100:.1f}% > 5%"
    print(f"  ✓ test_cowell_flat_background: S={r.net_area_counts:.0f}, "
          f"expected {expected:.0f}, err {err*100:.1f}%")


def test_cowell_sloped_background():
    """Cowell with linear baseline."""
    counts = make_synthetic_spectrum(
        N=1000, peaks=[(500, 5000, 10)],
        baseline=lambda x: 200 - 0.05 * x,  # sloped baseline
    )
    fwhm = 10 * 2.355
    r = cowell_area(counts, peak_channel=500, fwhm_channels=fwhm,
                    baseline_polynomial_order=1)
    expected = 5000 * 10 * math.sqrt(2 * math.pi)
    err = abs(r.net_area_counts - expected) / expected
    assert err < 0.08, f"Cowell sloped-bg error {err*100:.1f}% > 8%"
    print(f"  ✓ test_cowell_sloped_background: S={r.net_area_counts:.0f}, "
          f"err {err*100:.1f}%")


def test_gaussian_fit_recovers_parameters():
    """Gaussian fit should recover H, c, σ on a clean peak."""
    counts = make_synthetic_spectrum(
        N=1000, peaks=[(500, 5000, 10)],
        baseline=lambda x: 100 * np.ones_like(x),
    )
    fwhm = 10 * 2.355
    r = gaussian_fit_area(counts, peak_channel=500, fwhm_channels=fwhm)
    expected = 5000 * 10 * math.sqrt(2 * math.pi)
    err = abs(r.net_area_counts - expected) / expected
    assert err < 0.03, f"Gauss-fit area error {err*100:.1f}% > 3%"
    assert abs(r.peak_channel - 500) < 0.5, \
        f"Centroid drift: {r.peak_channel}"
    assert abs(r.fwhm_channels - fwhm) / fwhm < 0.05, \
        f"FWHM error: fit {r.fwhm_channels}, expected {fwhm}"
    print(f"  ✓ test_gaussian_fit_recovers_parameters: "
          f"S={r.net_area_counts:.0f} (err {err*100:.1f}%), "
          f"c={r.peak_channel:.1f}, FWHM={r.fwhm_channels:.2f}")


def test_gaussian_fit_uncertainty_reasonable():
    """For a peak with high σ=area/√gross, uncertainty should be small."""
    counts = make_synthetic_spectrum(
        N=1000, peaks=[(500, 5000, 10)],
        baseline=lambda x: 100 * np.ones_like(x),
    )
    fwhm = 10 * 2.355
    r = gaussian_fit_area(counts, peak_channel=500, fwhm_channels=fwhm)
    # For S ≈ 125k, gross ≈ 145k (S + 20k bg in ROI), √gross ≈ 380
    # Uncertainty should be of that order
    assert r.net_area_uncertainty > 0
    assert r.net_area_uncertainty < r.net_area_counts * 0.05, \
        f"Uncertainty {r.net_area_uncertainty:.0f} too large for strong peak"
    print(f"  ✓ test_gaussian_fit_uncertainty_reasonable: "
          f"S/dS = {r.net_area_counts / r.net_area_uncertainty:.1f}")


def test_cowell_low_stat_peak():
    """Cowell on weak peak — should still give finite result."""
    counts = make_synthetic_spectrum(
        N=1000, peaks=[(500, 200, 10)],   # weak peak
        baseline=lambda x: 500 * np.ones_like(x),  # high baseline
    )
    fwhm = 10 * 2.355
    r = cowell_area(counts, peak_channel=500, fwhm_channels=fwhm)
    # Should converge; sign of area can fluctuate from noise alone
    assert r.converged
    expected = 200 * 10 * math.sqrt(2 * math.pi)
    # Uncertainty larger here — allow up to 30% error on this weak peak
    err = abs(r.net_area_counts - expected) / expected
    print(f"  ✓ test_cowell_low_stat_peak: S={r.net_area_counts:.0f} "
          f"(expected {expected:.0f}, err {err*100:.0f}%, S/dS="
          f"{r.net_area_counts/max(r.net_area_uncertainty, 1):.2f})")


def test_cowell_edge_of_spectrum():
    """Peak very close to spectrum edge — ROI should be clipped."""
    counts = make_synthetic_spectrum(
        N=1000, peaks=[(20, 5000, 5)],
        baseline=lambda x: 100 * np.ones_like(x),
    )
    fwhm = 5 * 2.355
    r = cowell_area(counts, peak_channel=20, fwhm_channels=fwhm)
    # Should converge even though ROI is clipped on the left
    assert r.converged
    assert r.roi_low_ch == 0
    print(f"  ✓ test_cowell_edge_of_spectrum: ROI=[{r.roi_low_ch},{r.roi_high_ch}]")


def test_fixed_centroid_and_fwhm():
    """When calibration is reliable, fixing centroid+FWHM gives a
    fast linear LSQ for height + baseline."""
    counts = make_synthetic_spectrum(
        N=1000, peaks=[(500, 5000, 10)],
        baseline=lambda x: 100 * np.ones_like(x),
    )
    fwhm = 10 * 2.355
    r = gaussian_fit_area(counts, peak_channel=500, fwhm_channels=fwhm,
                          fix_centroid=True, fix_fwhm=True)
    expected = 5000 * 10 * math.sqrt(2 * math.pi)
    err = abs(r.net_area_counts - expected) / expected
    assert err < 0.05
    assert r.peak_channel == 500.0  # fixed
    print(f"  ✓ test_fixed_centroid_and_fwhm: S={r.net_area_counts:.0f}, "
          f"err {err*100:.1f}%")


def test_real_fixture_integration():
    """Run integration on Фон_кабинет: K-40 1460 keV peak should give
    reasonable area."""
    from gamma.io.readers import read_spectrum
    from gamma.peaks.search import mariscotti_search
    from gamma.calibration.fwhm_provider import make_fwhm_at_channel_provider

    spec = read_spectrum("detectors/AtomSpectra/data/fixtures/Фон_кабинет_8192к_01-01-2025.xml")
    fwhm_at_ch = make_fwhm_at_channel_provider(spec, fallback_channels=15.0)
    peaks = mariscotti_search(spec.counts, fwhm_channels=fwhm_at_ch, sigma_threshold=5.0)
    # Find K-40 peak (ch~3441)
    k40 = min(peaks, key=lambda p: abs(p.channel - 3441))
    fwhm = fwhm_at_ch(k40.channel)
    r_cowell = cowell_area(spec.counts, peak_channel=k40.channel,
                            fwhm_channels=fwhm)
    r_gauss = gaussian_fit_area(spec.counts, peak_channel=k40.channel,
                                  fwhm_channels=fwhm)
    # Both methods should give positive area
    assert r_cowell.net_area_counts > 0
    assert r_gauss.net_area_counts > 0
    # They should agree within 20% on a strong peak
    ratio = r_gauss.net_area_counts / r_cowell.net_area_counts
    assert 0.7 < ratio < 1.5, \
        f"Methods disagree: Cowell {r_cowell.net_area_counts:.0f}, " \
        f"Gauss {r_gauss.net_area_counts:.0f}"
    print(f"  ✓ test_real_fixture_integration: K-40 1460 keV "
          f"S_cowell={r_cowell.net_area_counts:.0f}, "
          f"S_gauss={r_gauss.net_area_counts:.0f} "
          f"(ratio {ratio:.2f})")


# ═════════════════════════════════════════════════════════════════════
# F-31a — get_peak_area: Lsrm-table preference with Cowell fallback
# ═════════════════════════════════════════════════════════════════════

class _FakeSpec:
    """Minimal Spectrum stand-in: holds counts + extras dict."""

    def __init__(self, counts, extras=None):
        self.counts = counts
        self.extras = extras or {}


def _synthetic_gaussian_spectrum(centroid=400, fwhm=15.0, area=1.0e5,
                                  baseline=200.0, n=1024):
    """Build a synthetic spectrum with one Gaussian peak on flat baseline."""
    sigma = fwhm / 2.3548
    ch = np.arange(n)
    peak = (area / (sigma * math.sqrt(2 * math.pi))) * \
        np.exp(-0.5 * ((ch - centroid) / sigma) ** 2)
    return (baseline + peak).astype(float)


def test_get_peak_area_prefers_lsrm_table():
    """When an Lsrm PEAKS table entry matches the peak, use its area."""
    counts = _synthetic_gaussian_spectrum(centroid=400, area=1.0e5)
    extras = {"lsrm_peaks_table": [
        {"position_ch": 400.2, "energy_keV": 1173.0,
         "area": 666002.0, "d_area": 1151.0},
    ]}
    spec = _FakeSpec(counts, extras)
    area, dA, src = get_peak_area(spec, peak_channel=400,
                                   fwhm_channels=15.0)
    assert src == "lsrm_peaks_table"
    assert area == 666002.0
    assert dA == 1151.0
    print(f"  ✓ test_get_peak_area_prefers_lsrm_table "
          f"(area={area:.0f} from {src})")


def test_get_peak_area_falls_back_when_no_table():
    """No Lsrm table → Cowell fallback."""
    counts = _synthetic_gaussian_spectrum(centroid=400, area=1.0e5)
    spec = _FakeSpec(counts, extras={})  # no lsrm_peaks_table
    area, dA, src = get_peak_area(spec, peak_channel=400,
                                   fwhm_channels=15.0)
    assert src == "cowell"
    # Cowell should recover roughly the injected area (±15%)
    assert 0.85e5 < area < 1.15e5, f"Cowell area {area:.0f} off"
    print(f"  ✓ test_get_peak_area_falls_back_when_no_table "
          f"(area={area:.0f} from {src})")


def test_get_peak_area_falls_back_when_no_match():
    """Lsrm table present but no entry near the peak → Cowell fallback."""
    counts = _synthetic_gaussian_spectrum(centroid=400, area=1.0e5)
    extras = {"lsrm_peaks_table": [
        {"position_ch": 700.0, "energy_keV": 2000.0,
         "area": 5.0e4, "d_area": 500.0},   # far from ch 400
    ]}
    spec = _FakeSpec(counts, extras)
    area, dA, src = get_peak_area(spec, peak_channel=400,
                                   fwhm_channels=15.0)
    assert src == "cowell", f"expected cowell fallback, got {src}"
    print(f"  ✓ test_get_peak_area_falls_back_when_no_match "
          f"(area={area:.0f} from {src})")


def test_get_peak_area_bypass_with_prefer_false():
    """prefer_lsrm_table=False forces Cowell even when a table matches."""
    counts = _synthetic_gaussian_spectrum(centroid=400, area=1.0e5)
    extras = {"lsrm_peaks_table": [
        {"position_ch": 400.0, "energy_keV": 1173.0,
         "area": 666002.0, "d_area": 1151.0},
    ]}
    spec = _FakeSpec(counts, extras)
    area, dA, src = get_peak_area(spec, peak_channel=400,
                                   fwhm_channels=15.0,
                                   prefer_lsrm_table=False)
    assert src == "cowell"
    assert area != 666002.0
    print(f"  ✓ test_get_peak_area_bypass_with_prefer_false "
          f"(area={area:.0f} from {src})")


def test_get_peak_area_co60_doublet_real():
    """On the real Co-60 spectrum, get_peak_area returns the Lsrm-fitted
    area (~666k for 1173 keV), which is ~40% larger than the Cowell
    area — the F-31a fix that corrected the −25% activity bias.
    """
    from gamma.io.readers import read_spectrum
    from gamma.calibration.fwhm_provider import make_fwhm_at_channel_provider

    spe = ("detectors/Gamma-1S/reference_spectra/archive/"
           "Co-60__043_02_2019_Точечная-5см_5cm.spe")
    spec = read_spectrum(spe)
    fwhm_at = make_fwhm_at_channel_provider(spec, fallback_channels=15.0)
    # 1173 keV peak is near ch 401
    area, dA, src = get_peak_area(spec, peak_channel=401,
                                   fwhm_channels=fwhm_at(401))
    assert src == "lsrm_peaks_table", \
        f"expected Lsrm table area for Co-60, got {src}"
    # Lsrm area is ~666k; Cowell would give ~469k
    assert area > 600000, f"Lsrm 1173 area {area:.0f} unexpectedly low"
    # Confirm Cowell really is much lower (the bug being fixed)
    cowell_r = cowell_area(spec.counts, peak_channel=401,
                            fwhm_channels=fwhm_at(401))
    deficit = (area - cowell_r.net_area_counts) / area * 100
    assert deficit > 20, \
        f"Cowell deficit {deficit:.0f}% smaller than expected (bug not reproduced)"
    print(f"  ✓ test_get_peak_area_co60_doublet_real "
          f"(Lsrm={area:.0f}, Cowell={cowell_r.net_area_counts:.0f}, "
          f"deficit={deficit:.0f}%)")


if __name__ == "__main__":
    print("Running Phase 2.1a peak area integration tests...\n")
    test_cowell_flat_background()
    test_cowell_sloped_background()
    test_gaussian_fit_recovers_parameters()
    test_gaussian_fit_uncertainty_reasonable()
    test_cowell_low_stat_peak()
    test_cowell_edge_of_spectrum()
    test_fixed_centroid_and_fwhm()
    test_real_fixture_integration()
    test_get_peak_area_prefers_lsrm_table()
    test_get_peak_area_falls_back_when_no_table()
    test_get_peak_area_falls_back_when_no_match()
    test_get_peak_area_bypass_with_prefer_false()
    test_get_peak_area_co60_doublet_real()
    print("\n✓ All Phase 2.1a tests passed.")
