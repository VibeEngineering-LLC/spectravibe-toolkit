"""
F-449 (operator-locked 2026-06-16) regression tests.

Two contracts under test:

A. ``gamma.peaks.area_step_continuum.gauss_erfc_step_fit`` sigma-lock:
   By default (LSRM "PSHPV" flag OFF) the peak width sigma is HARD-LOCKED
   to the calibration FWHM(E); the reported ``fwhm_channels`` equals the
   calibration FWHM, NOT the fitted width. Setting GAMMA_FREE_SIGMA=1
   restores the legacy free-sigma mode, in which the fitted width tracks
   the true peak width and therefore differs from the calibration FWHM.

B. ``gamma.calibration.fwhm_provider.make_fwhm_at_channel_provider``
   F-449 bootstrap: with ``bootstrap_from_peaks=True`` and no vendor
   FWHM model / calibration peaks, the provider fits FWHM(E) from the
   significant, isolated peaks of the spectrum (source-label
   "bootstrap_significant_peaks_*") instead of returning the generic
   constant. With ``bootstrap_from_peaks=False`` (default) the legacy
   constant fallback is preserved.

References: CLAUDE.md "FWHM lyubogo pika = kalibrovka FWHM(E)"
(operator-locked 2026-06-16); spectralinexx_2.0_basic_functions_rus.pdf
:1647-1650 (PSHPV flag binary semantics).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

# NB: import via gamma.peaks.area (which re-exports gauss_erfc_step_fit)
# to avoid the area <-> area_step_continuum circular-import trap that
# triggers when the submodule is imported first.
from gamma.peaks.area import gauss_erfc_step_fit
from gamma.calibration.fwhm_provider import make_fwhm_at_channel_provider


FWHM_TO_SIGMA = 2.355  # project constant used by the fit module


def _make_peak(center, true_sigma, height, n=400, baseline=80.0, slope=0.0):
    """Single Gaussian on a (optionally sloped) flat continuum, no noise.

    Noise-free so the locked vs free distinction is unambiguous: a free
    fit recovers ``true_sigma`` exactly, a locked fit reports the
    calibration sigma regardless.
    """
    x = np.arange(n, dtype=np.float64)
    counts = baseline + slope * (x - center)
    counts = counts + height * np.exp(-((x - center) ** 2) / (2.0 * true_sigma ** 2))
    return np.maximum(counts, 0.0)


# ----------------------------------------------------------------------
# A. sigma-lock contract on gauss_erfc_step_fit
# ----------------------------------------------------------------------

def test_locked_sigma_reports_calibration_fwhm(monkeypatch):
    """Default mode: reported fwhm_channels == calibration FWHM, even when
    the true peak is narrower/wider than the calibration says."""
    monkeypatch.delenv("GAMMA_FREE_SIGMA", raising=False)

    center = 200
    true_sigma = 6.0           # actual peak width
    cal_fwhm = 18.0            # calibration says wider (cal_sigma = 7.643)
    counts = _make_peak(center, true_sigma, height=5000.0)

    res = gauss_erfc_step_fit(
        counts, peak_channel=center, fwhm_channels=cal_fwhm,
    )
    assert res.converged, f"fit did not converge: {res.notes}"
    # Reported FWHM must equal the calibration FWHM (sigma locked).
    assert abs(res.fwhm_channels - cal_fwhm) < 1e-6, (
        f"locked mode must report calibration FWHM {cal_fwhm}, "
        f"got {res.fwhm_channels}"
    )
    assert "locked(F-449)" in res.notes


def test_free_sigma_recovers_true_width(monkeypatch):
    """GAMMA_FREE_SIGMA=1: fitted width tracks the true peak width and
    therefore differs from the calibration FWHM."""
    monkeypatch.setenv("GAMMA_FREE_SIGMA", "1")

    center = 200
    true_sigma = 6.0
    cal_fwhm = 18.0            # cal_sigma 7.643 != true_sigma 6.0
    true_fwhm = true_sigma * FWHM_TO_SIGMA
    counts = _make_peak(center, true_sigma, height=5000.0)

    res = gauss_erfc_step_fit(
        counts, peak_channel=center, fwhm_channels=cal_fwhm,
    )
    assert res.converged, f"fit did not converge: {res.notes}"
    # Free fit must recover the true width (bounds [0.6,1.6]*cal_sigma
    # admit true_sigma=6.0 since 0.6*7.643=4.59 <= 6.0 <= 12.23).
    assert abs(res.fwhm_channels - true_fwhm) / true_fwhm < 0.05, (
        f"free mode should recover true FWHM {true_fwhm:.2f}, "
        f"got {res.fwhm_channels:.2f}"
    )
    # And it must differ clearly from the calibration FWHM.
    assert abs(res.fwhm_channels - cal_fwhm) > 0.5, (
        "free fit width should differ from calibration FWHM"
    )
    assert "sigma_mode=free" in res.notes


def test_locked_and_free_differ(monkeypatch):
    """Same input, the two modes must produce different reported FWHM."""
    center, true_sigma, cal_fwhm = 200, 6.0, 18.0
    counts = _make_peak(center, true_sigma, height=5000.0)

    monkeypatch.delenv("GAMMA_FREE_SIGMA", raising=False)
    locked = gauss_erfc_step_fit(counts, peak_channel=center, fwhm_channels=cal_fwhm)

    monkeypatch.setenv("GAMMA_FREE_SIGMA", "1")
    free = gauss_erfc_step_fit(counts, peak_channel=center, fwhm_channels=cal_fwhm)

    assert abs(locked.fwhm_channels - free.fwhm_channels) > 0.5, (
        f"locked ({locked.fwhm_channels:.2f}) and free "
        f"({free.fwhm_channels:.2f}) FWHM should differ"
    )


def test_locked_area_consistent(monkeypatch):
    """Locked-mode net area = H*sigma_cal*sqrt(2pi); area_var uses the
    constant-sigma propagation (sigma*sqrt(2pi))^2 * H_var."""
    monkeypatch.delenv("GAMMA_FREE_SIGMA", raising=False)
    center, true_sigma, cal_fwhm = 200, 7.643, 18.0   # cal_sigma == true_sigma
    height = 5000.0
    counts = _make_peak(center, true_sigma, height=height)
    res = gauss_erfc_step_fit(counts, peak_channel=center, fwhm_channels=cal_fwhm)
    assert res.converged
    cal_sigma = cal_fwhm / FWHM_TO_SIGMA
    expected_area = height * cal_sigma * math.sqrt(2 * math.pi)
    err = abs(res.net_area_counts - expected_area) / expected_area
    assert err < 0.05, (
        f"locked area {res.net_area_counts:.0f} vs expected "
        f"{expected_area:.0f} ({err*100:.1f}%)"
    )
    assert res.net_area_uncertainty >= 0.0


# ----------------------------------------------------------------------
# B. fwhm_provider bootstrap-from-significant-peaks
# ----------------------------------------------------------------------

class _MiniSpec:
    """Minimal Spectrum-like object for the channel provider."""

    def __init__(self, counts, energy_cal):
        self.counts = counts
        self.energy_cal = energy_cal
        self.stored_fwhm_calibration = None

    def channel_to_energy(self, ch):
        return sum(a * ch ** i for i, a in enumerate(self.energy_cal))


def _make_nai_spectrum(a1=2.0, n=1024):
    """Synthetic NaI-like spectrum: linear E-cal, FWHM_keV = alpha*sqrt(E)."""
    counts = np.full(n, 20.0)

    def fwhm_keV_true(E):
        return 0.06 * math.sqrt(max(E, 1.0) * 100.0)

    for E, A in [(200.0, 4000.0), (600.0, 6000.0), (1000.0, 5000.0),
                 (1400.0, 3500.0), (1800.0, 3000.0)]:
        ch0 = E / a1
        sig = (fwhm_keV_true(E) / a1) / 2.3548
        x = np.arange(n)
        counts = counts + A * np.exp(-((x - ch0) ** 2) / (2 * sig * sig))
    return np.maximum(counts, 0.0), fwhm_keV_true


def test_bootstrap_off_is_constant_fallback():
    """Default (bootstrap_from_peaks=False) preserves constant fallback."""
    counts, _ = _make_nai_spectrum()
    spec = _MiniSpec(counts, (0.0, 2.0))
    prov = make_fwhm_at_channel_provider(spec, fallback_channels=15.0)
    assert getattr(prov, "fwhm_source", None) == "fallback"
    # Constant everywhere (above any physical floor).
    assert abs(prov(300) - 15.0) < 1e-6
    assert abs(prov(700) - 15.0) < 1e-6


def test_bootstrap_recovers_fwhm_curve():
    """bootstrap_from_peaks=True fits FWHM(E) from significant peaks."""
    counts, fwhm_keV_true = _make_nai_spectrum(a1=2.0)
    spec = _MiniSpec(counts, (0.0, 2.0))
    prov = make_fwhm_at_channel_provider(
        spec, fallback_channels=15.0, bootstrap_from_peaks=True,
    )
    src = getattr(prov, "fwhm_source", "")
    assert src.startswith("bootstrap_significant_peaks"), (
        f"expected bootstrap source, got {src!r}"
    )
    assert getattr(prov, "bootstrap_n_anchors", 0) >= 3
    # FWHM(E) recovered within 5% at several energies.
    a1 = 2.0
    for E in (200, 600, 1000, 1400, 1800):
        ch = int(round(E / a1))
        fwhm_keV = prov(ch) * a1
        tru = fwhm_keV_true(E)
        err = abs(fwhm_keV - tru) / tru
        assert err < 0.05, f"E={E}: boot {fwhm_keV:.2f} vs true {tru:.2f} ({err*100:.1f}%)"


def test_bootstrap_declines_without_energy_cal():
    """No energy_cal => bootstrap cannot run => constant fallback."""
    counts, _ = _make_nai_spectrum()
    spec = _MiniSpec(counts, None)
    prov = make_fwhm_at_channel_provider(
        spec, fallback_channels=15.0, bootstrap_from_peaks=True,
    )
    assert getattr(prov, "fwhm_source", None) == "fallback"


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-v"]))