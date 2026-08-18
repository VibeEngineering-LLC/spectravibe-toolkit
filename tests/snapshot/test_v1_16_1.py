# -*- coding: utf-8 -*-
"""
v1.16.1 delivery tests — F-95 dead-time + F-96 bg-lines a priori.

Run:
    cd 0_Work/gamma-spectrum-analysis
    $env:PYTHONPATH = "scripts"
    python test_v1_16_1.py
"""
from __future__ import annotations
import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts"))

import numpy as np

from gamma.physics.dead_time import (
    DeadTimeCoefficients, effective_dead_time,
    calibrate_AB, make_calib_point, get_detector_coeffs,
)
from gamma.physics.bg_lines_apriori import (
    BackgroundPeak, BackgroundLineLibrary,
    build_library_from_peak_list,
    inject_priors_into_roi,
)

# ──────────────────────────────────────────────────────────────────────
# F-95 — dead-time formula and calibration
# ──────────────────────────────────────────────────────────────────────

def test_dead_time_no_coeffs_returns_uncorrected() -> None:
    """When coeffs unknown, t_d=0 and applied=False with a note."""
    counts = [10.0] * 4096
    r = effective_dead_time(counts, t_live_s=1000.0, coeffs=None)
    assert not r.applied
    assert r.t_dead_s == 0.0
    assert r.t_live_corr_s == 1000.0
    assert "UNCORRECTED" in r.notes


def test_dead_time_applied_with_known_coeffs() -> None:
    """t_d = A·Σy + B·Σ(y·i); for uniform counts, both sums positive."""
    coeffs = DeadTimeCoefficients(
        A_s_per_count=1e-6, B_s_per_count_channel=1e-10,
    )
    n_ch = 4096
    counts = [10.0] * n_ch
    r = effective_dead_time(counts, t_live_s=1000.0, coeffs=coeffs)
    Sigma_y = 10.0 * n_ch
    # Σ(y·i) = 10 · Σi  for i in [0..n_ch-1] = 10 · n_ch·(n_ch-1)/2
    Sigma_yi = 10.0 * n_ch * (n_ch - 1) / 2
    expected_td = 1e-6 * Sigma_y + 1e-10 * Sigma_yi
    assert abs(r.t_dead_s - expected_td) < 1e-9
    assert abs(r.t_live_corr_s - (1000.0 - expected_td)) < 1e-9
    assert r.applied


def test_dead_time_calibrate_AB_round_trip() -> None:
    """
    Build three synthetic spectra with known A, B → calibrate_AB →
    recovered values agree with input.
    """
    A_true = 5e-6
    B_true = 1e-9
    n_ch = 1024
    rng = np.random.default_rng(0)

    def make_spec(scale: float):
        # rate around 'scale' per channel
        return scale + rng.normal(0, max(0.01, scale * 0.05), size=n_ch)

    # Low load: rate ≈ 500 cps · t_live = 500·1000 = 5e5 counts spread over channels
    low = make_spec(120.0)        # 120·1024 ≈ 1.2e5 total
    # High loads: ~50 × bigger counts
    high_a = make_spec(5000.0)
    high_b = make_spec(4500.0)

    # Reference "rate at line" — for synthetic test we just pick a
    # nominal y_ref and let dead-time apply self-consistently. To
    # match the formulation, we need actual rate-at-ref-line; here
    # we mock it: assume y_ref_low = 100 cps and observed reduces by
    # the dead-time multiplier.
    t_live = 1000.0
    y_ref_low = 100.0
    # apply dead time to compute the reduced rate at high loads
    S_a = float(np.sum(high_a))
    M_a = float(np.sum(high_a * np.arange(n_ch, dtype=float)))
    S_b = float(np.sum(high_b))
    M_b = float(np.sum(high_b * np.arange(n_ch, dtype=float)))
    t_d_a = A_true * S_a + B_true * M_a
    t_d_b = A_true * S_b + B_true * M_b
    # Reduced rates
    y_ref_a = y_ref_low * (1 - t_d_a / t_live)
    y_ref_b = y_ref_low * (1 - t_d_b / t_live)

    p_low = make_calib_point("low", low, t_live_s=t_live, rate_at_ref_keV=y_ref_low)
    p_a = make_calib_point("hiA", high_a, t_live_s=t_live, rate_at_ref_keV=y_ref_a)
    p_b = make_calib_point("hiB", high_b, t_live_s=t_live, rate_at_ref_keV=y_ref_b)

    rec = calibrate_AB(p_low, p_a, p_b)
    assert abs(rec.A_s_per_count - A_true) / A_true < 0.05, \
        f"A drift: {rec.A_s_per_count} vs {A_true}"
    assert abs(rec.B_s_per_count_channel - B_true) / B_true < 0.05, \
        f"B drift: {rec.B_s_per_count_channel} vs {B_true}"


def test_dead_time_singular_calibration_raises() -> None:
    """Identical high-load points → singular 2×2 → ValueError."""
    n_ch = 256
    counts = [1.0] * n_ch
    p_low = make_calib_point("low", counts, t_live_s=1000.0, rate_at_ref_keV=100.0)
    p_a = make_calib_point("a", counts, t_live_s=1000.0, rate_at_ref_keV=90.0)
    p_b = make_calib_point("b", counts, t_live_s=1000.0, rate_at_ref_keV=90.0)
    try:
        calibrate_AB(p_low, p_a, p_b)
    except ValueError as e:
        assert "singular" in str(e).lower()
    else:
        assert False, "expected ValueError for singular calibration"


def test_gamma_1c_default_uncalibrated() -> None:
    """Gamma-1S is uncalibrated by default; correction skipped."""
    coeffs = get_detector_coeffs("Gamma-1S")
    assert coeffs is None
    r = effective_dead_time([1.0] * 100, t_live_s=10.0, coeffs=coeffs)
    assert not r.applied


# ──────────────────────────────────────────────────────────────────────
# F-96 — background lines as a priori
# ──────────────────────────────────────────────────────────────────────

def _typical_lib() -> BackgroundLineLibrary:
    """Mock a typical natural-background line list for Gamma-1S."""
    return build_library_from_peak_list([
        {"E_keV": 1460.82, "rate_cps": 0.0080, "sigma_rate_cps": 0.0003, "nuclide_hint": "K-40"},
        {"E_keV": 609.32,  "rate_cps": 0.0025, "sigma_rate_cps": 0.0002, "nuclide_hint": "Bi-214"},
        {"E_keV": 2614.51, "rate_cps": 0.0011, "sigma_rate_cps": 0.0001, "nuclide_hint": "Tl-208"},
        {"E_keV": 911.20,  "rate_cps": 0.0009, "sigma_rate_cps": 0.0001, "nuclide_hint": "Ac-228"},
        {"E_keV": 511.0,   "rate_cps": 0.0030, "sigma_rate_cps": 0.0002, "nuclide_hint": "annihil"},
    ], source_label="natural-background-Marinelli", t_live_source_s=86400.0)


def test_bg_library_peaks_in_range() -> None:
    lib = _typical_lib()
    sub = lib.peaks_in_range(500.0, 700.0)
    energies = sorted(p.E_keV for p in sub)
    assert energies == [511.0, 609.32]


def test_bg_library_peak_near() -> None:
    lib = _typical_lib()
    p = lib.peak_near(1461.0, window_keV=5.0)
    assert p is not None and p.nuclide_hint == "K-40"
    none_p = lib.peak_near(700.0, window_keV=5.0)
    assert none_p is None


def test_bg_inject_priors_for_roi() -> None:
    lib = _typical_lib()
    inj = inject_priors_into_roi(
        lib,
        E_low_keV=400.0, E_high_keV=700.0,
        t_live_sample_s=1000.0,
    )
    # Two bg peaks (511, 609.32) in this ROI
    assert inj.n_extra() == 2
    # 511 keV at rate 0.0030 cps over 1000 s → 3 counts
    assert any(abs(amp - 3.0) < 0.01 for amp in inj.expected_amplitudes_counts)


def test_bg_inject_skips_overlapping_sample_peaks() -> None:
    """When the sample fit already has a peak at the same E, drop the prior."""
    lib = _typical_lib()
    inj = inject_priors_into_roi(
        lib,
        E_low_keV=400.0, E_high_keV=700.0,
        t_live_sample_s=1000.0,
        existing_peak_positions_keV=(511.5,),   # close to 511.0
        min_separation_keV=2.0,
    )
    # 511 dropped, 609.32 kept
    assert inj.n_extra() == 1
    assert abs(inj.extra_peak_positions_keV[0] - 609.32) < 0.01


def test_bg_amplitude_scales_with_sample_t_live() -> None:
    lib = _typical_lib()
    inj_short = inject_priors_into_roi(
        lib, E_low_keV=1450.0, E_high_keV=1470.0, t_live_sample_s=100.0,
    )
    inj_long = inject_priors_into_roi(
        lib, E_low_keV=1450.0, E_high_keV=1470.0, t_live_sample_s=10000.0,
    )
    # K-40 at 0.0080 cps: 0.80 vs 80 counts
    assert abs(inj_short.expected_amplitudes_counts[0] - 0.80) < 0.01
    assert abs(inj_long.expected_amplitudes_counts[0] - 80.0) < 0.01


# ──────────────────────────────────────────────────────────────────────
# Driver
# ──────────────────────────────────────────────────────────────────────

def _all_tests():
    return [v for k, v in globals().items() if k.startswith("test_") and callable(v)]


def main() -> int:
    tests = _all_tests()
    fail = 0
    for t in tests:
        try:
            t()
            print(f"  PASS  {t.__name__}")
        except AssertionError as exc:
            fail += 1
            print(f"  FAIL  {t.__name__}: {exc}")
        except Exception as exc:
            fail += 1
            print(f"  ERR   {t.__name__}: {type(exc).__name__}: {exc}")
    print(f"\nv1.16.1: {len(tests) - fail}/{len(tests)} passed")
    return fail


if __name__ == "__main__":
    sys.exit(main())
