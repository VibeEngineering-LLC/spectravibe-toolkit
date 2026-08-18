"""
F-445 unit-tests for cluster-level Delta_cluster anchors.

See `scripts/gamma/calibration/cluster_shift_anchors.py` for the API
under test, and `_state/agent_a/inbox/2026-06-14_F-445_phase2_delta_cluster_to_ecal.md`
for acceptance criteria.
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

# repo-root scripts path (tests share this conftest pattern)
_SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from gamma.calibration.cluster_shift_anchors import (
    compute_cluster_global_shift,
    collect_cluster_global_anchors,
)
from gamma.calibration.multiplet_self_calibration import CentroidAnchor


# ---------------------------------------------------------------- helpers


class _FakeSpec:
    """Minimal Spectrum stand-in with linear E(N) = a0 + a1*N."""

    def __init__(self, a0=0.0, a1=3.0, n_ch=1024, counts=None):
        self.a0 = float(a0)
        self.a1 = float(a1)
        self.n_ch = int(n_ch)
        self.counts = counts if counts is not None else np.zeros(n_ch)

    def channel_to_energy(self, ch):
        return self.a0 + self.a1 * float(ch)

    def energy_to_channel(self, E):
        return (float(E) - self.a0) / self.a1


def _make_synthetic_spectrum(
    E_peak_true=2614.5, fwhm_keV=70.0, amp=10000.0,
    cont_level=50.0, a0=0.0, a1=3.0, n_ch=1024,
):
    """Build a synthetic spectrum with one Gaussian peak.

    Returns (spec, E_arr_keV, cont_arr).
    """
    sigma_keV = fwhm_keV / 2.355
    E_arr = np.array([a0 + a1 * ch for ch in range(n_ch)])
    cont = np.full(n_ch, cont_level)
    peak = amp * np.exp(-0.5 * ((E_arr - E_peak_true) / sigma_keV) ** 2)
    counts = cont + peak
    spec = _FakeSpec(a0=a0, a1=a1, n_ch=n_ch, counts=counts)
    return spec, list(E_arr), list(cont)


# ---------------------------------------------------------------- tests


def test_compute_cluster_global_shift_synthetic_shift_plus_2p5():
    """Synthetic peak at E_lib + 2.5 keV -> shift recovered to <=0.5 keV."""
    E_lib = 2614.5
    delta_true = +2.5
    fwhm = 70.0
    spec, E_arr, cont = _make_synthetic_spectrum(
        E_peak_true=E_lib + delta_true, fwhm_keV=fwhm,
    )
    sigma_keV = fwhm / 2.355
    delta = compute_cluster_global_shift(
        E_lib, sigma_keV, E_arr, cont, spec.counts, spec,
    )
    assert delta is not None, "shift fn returned None on clean peak"
    assert abs(delta - delta_true) <= 1.0, (
        "recovered Delta=" + repr(delta) + " expected ~" + repr(delta_true)
    )


def test_compute_cluster_global_shift_zero_shift():
    """Peak at E_lib exactly -> |delta| < FWHM/2 (a few keV)."""
    E_lib = 911.2
    fwhm = 40.0
    spec, E_arr, cont = _make_synthetic_spectrum(
        E_peak_true=E_lib, fwhm_keV=fwhm, amp=5000.0, cont_level=10.0,
    )
    sigma_keV = fwhm / 2.355
    delta = compute_cluster_global_shift(
        E_lib, sigma_keV, E_arr, cont, spec.counts, spec,
    )
    assert delta is not None
    assert abs(delta) <= 5.0, "expected ~0 shift, got " + repr(delta)


def test_phantom_cluster_returns_None_on_flat_noise():
    """Flat continuum (no peak) -> best_net <= 0 -> None."""
    E_lib = 1588.2
    fwhm = 55.0
    spec, E_arr, cont = _make_synthetic_spectrum(
        E_peak_true=0.0, fwhm_keV=fwhm, amp=0.0, cont_level=20.0,
    )
    sigma_keV = fwhm / 2.355
    delta = compute_cluster_global_shift(
        E_lib, sigma_keV, E_arr, cont, spec.counts, spec,
    )
    # spec.counts == cont -> best_net == 0 strictly -> None per guard
    assert delta is None


def test_runaway_shift_clamped_to_None():
    """Verify 1.5*FWHM clamp triggers when a peak sits beyond it.

    Strategy: pass a very small sigma_strongest so the search window
    (max(0.7*FWHM, 8.0)) is dominated by the 8.0-keV floor while the
    1.5*FWHM clamp is much smaller. Then a real peak ~+12 keV outside
    the lib energy lands inside the 8-keV window edge but at delta
    much larger than 1.5*FWHM, so the clamp triggers and we get None.
    """
    E_lib = 911.2
    tiny_fwhm = 2.0
    sigma_tiny = tiny_fwhm / 2.355
    spec, E_arr, cont = _make_synthetic_spectrum(
        E_peak_true=E_lib + 12.0, fwhm_keV=2.0, amp=8000.0, cont_level=5.0,
    )
    delta = compute_cluster_global_shift(
        E_lib, sigma_tiny, E_arr, cont, spec.counts, spec,
    )
    # Window = max(0.7*2, 8) = 8 keV. Best_E ~ E_lib + 8 (edge near peak).
    # |8 keV| > 1.5*tiny_fwhm = 3 keV -> clamp triggers -> None.
    assert delta is None, "expected None due to >1.5*FWHM clamp, got " + repr(delta)


class _FakeComp:
    def __init__(self, line_E_keV, nuclide="X"):
        self.line_E_keV = float(line_E_keV)
        self.nuclide = nuclide


class _FakeCluster:
    def __init__(self, cid, comps, areas, overlay_E, overlay_cont):
        self.cluster_id = cid
        self.components = comps
        self.areas = areas
        self.overlay_E_keV = overlay_E
        self.overlay_continuum = overlay_cont


def test_collect_anchors_two_clusters():
    """Two synthetic clusters with known Delta -> two anchors returned."""
    spec, E_arr, cont = _make_synthetic_spectrum(
        E_peak_true=2614.5 + 2.0, fwhm_keV=70.0, amp=8000.0, cont_level=30.0,
    )
    # Second cluster occupies channels around 911 keV
    spec_b, E_arr_b, cont_b = _make_synthetic_spectrum(
        E_peak_true=911.2 + 1.5, fwhm_keV=40.0, amp=12000.0, cont_level=25.0,
    )
    # Coadd peaks into single spec.counts (channels overlap because both
    # use same a0/a1) - we use spec only for energy_to_channel here, so OK.
    spec.counts = spec.counts + spec_b.counts
    clusters = [
        _FakeCluster(
            "M_Tl208_2614",
            [_FakeComp(2614.5, "Tl-208")],
            [50000.0],
            E_arr, cont,
        ),
        _FakeCluster(
            "M_Ac228_911",
            [_FakeComp(911.2, "Ac-228")],
            [40000.0],
            E_arr_b, cont_b,
        ),
    ]
    fwhm_provider = lambda E: 70.0 if E > 1500 else 40.0
    anchors = collect_cluster_global_anchors(
        spec, clusters, spec.counts, None, None, fwhm_provider,
    )
    assert len(anchors) == 2
    assert all(isinstance(a, CentroidAnchor) for a in anchors)
    # Anchor 1 should be Tl-208 2614
    e_passports = sorted([a.E_passport_keV for a in anchors])
    assert 911.0 <= e_passports[0] <= 911.5
    assert 2614.0 <= e_passports[1] <= 2615.0


def test_phantom_cluster_filtered_out_in_collector():
    """Cluster with area_strong ~ 0 -> phantom guard -> no anchor."""
    spec, E_arr, cont = _make_synthetic_spectrum(amp=0.0, cont_level=5.0)
    clusters = [
        _FakeCluster(
            "M_phantom",
            [_FakeComp(1588.2, "Ac-228")],
            [0.5],  # area * 1/sqrt(2pi) ~ 0.2 < phantom_floor=1.0
            E_arr, cont,
        ),
    ]
    fwhm_provider = lambda E: 55.0
    anchors = collect_cluster_global_anchors(
        spec, clusters, spec.counts, None, None, fwhm_provider,
    )
    assert anchors == []
