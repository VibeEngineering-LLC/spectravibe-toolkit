"""
Tests for Phase 2.1e background subtraction:
  - Explicit consent requirement (safety policy)
  - Rate-normalized channel mode (matched gains)
  - Energy-aligned mode (gain mismatch)
  - Uncertainty propagation
  - Real Gamma-1S reference data validation
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

import numpy as np

from gamma.io.readers import read_spectrum
from gamma.calibration.bg_subtract_dual_mode import (
    subtract_background, apply_subtraction_to_spectrum,
    BackgroundConsentRequired, BackgroundConsentRegistry,
    GAIN_MATCH_THRESHOLD,
)
from gamma.spectrum import Spectrum


REF_DIR = "detectors/Gamma-1S/reference_spectra/archive"


# --- Synthetic spectrum helper ---

def make_synthetic(counts, live_time, a0=0, a1=1.0, source_path="test"):
    spec = Spectrum(
        counts=np.array(counts, dtype=np.int64),
        live_time=live_time,
        real_time=live_time,
        source_path=source_path,
        source_format="test",
    )
    spec.energy_cal = (a0, a1)
    spec.n_channels = len(counts)
    spec.n_channels_raw = len(counts)
    return spec


# --- Safety policy ---

def test_subtract_without_consent_raises():
    """Subtract should refuse without user_confirmed_applicable=True."""
    src = make_synthetic([100, 200, 300], 100, a1=1.0)
    bg = make_synthetic([10, 20, 30], 100, a1=1.0)
    try:
        subtract_background(src, bg)
        assert False, "Should have raised BackgroundConsentRequired"
    except BackgroundConsentRequired as e:
        assert "explicit user confirmation" in str(e)
    print(f"  ✓ test_subtract_without_consent_raises")


def test_subtract_with_consent_works():
    """With consent flag, subtraction proceeds."""
    src = make_synthetic([100, 200, 300], 100, a1=1.0)
    bg = make_synthetic([10, 20, 30], 100, a1=1.0)
    result = subtract_background(src, bg, user_confirmed_applicable=True)
    assert result.mode == "rate_normalized_channel"
    # net = source - bg * (100/100) = [90, 180, 270]
    np.testing.assert_array_equal(result.subtracted_counts, [90, 180, 270])
    print(f"  ✓ test_subtract_with_consent_works")


# --- Consent registry (session-level approval) ---

def test_consent_registry_empty_raises():
    """An empty registry behaves like no consent."""
    registry = BackgroundConsentRegistry()
    src = make_synthetic([100], 100, a1=1.0, source_path="/src.spe")
    bg = make_synthetic([10], 100, a1=1.0, source_path="/bg.spe")
    try:
        subtract_background(src, bg, consent_registry=registry)
        assert False, "Should have raised"
    except BackgroundConsentRequired as e:
        # Exception should carry the bg path so caller can prompt
        assert e.background_path.endswith("bg.spe"), f"got {e.background_path}"
    print(f"  ✓ test_consent_registry_empty_raises")


def test_consent_registry_after_approve_proceeds():
    """After approve(), subtraction proceeds without re-prompt."""
    registry = BackgroundConsentRegistry()
    src = make_synthetic([100], 100, a1=1.0, source_path="/src.spe")
    bg = make_synthetic([10], 100, a1=1.0, source_path="/bg.spe")
    registry.approve(bg.source_path)
    result = subtract_background(src, bg, consent_registry=registry)
    np.testing.assert_array_equal(result.subtracted_counts, [90])
    print(f"  ✓ test_consent_registry_after_approve_proceeds")


def test_consent_registry_per_file_isolation():
    """Approval for one bg file does not transfer to another bg file."""
    registry = BackgroundConsentRegistry()
    src = make_synthetic([100], 100, a1=1.0, source_path="/src.spe")
    bg1 = make_synthetic([10], 100, a1=1.0, source_path="/bg1.spe")
    bg2 = make_synthetic([5], 100, a1=1.0, source_path="/bg2.spe")
    registry.approve(bg1.source_path)
    # bg1 should proceed
    result1 = subtract_background(src, bg1, consent_registry=registry)
    assert result1.subtracted_counts[0] == 90
    # bg2 should raise (not approved)
    try:
        subtract_background(src, bg2, consent_registry=registry)
        assert False, "Should have raised for bg2"
    except BackgroundConsentRequired as e:
        assert e.background_path.endswith("bg2.spe")
    print(f"  ✓ test_consent_registry_per_file_isolation")


def test_consent_registry_multi_source_reuse():
    """One approval works for multiple source spectra (same bg)."""
    registry = BackgroundConsentRegistry()
    src1 = make_synthetic([100], 100, a1=1.0, source_path="/src1.spe")
    src2 = make_synthetic([200], 100, a1=1.0, source_path="/src2.spe")
    bg = make_synthetic([10], 100, a1=1.0, source_path="/bg.spe")
    registry.approve(bg.source_path)
    # Both source spectra should work without re-prompt
    r1 = subtract_background(src1, bg, consent_registry=registry)
    r2 = subtract_background(src2, bg, consent_registry=registry)
    assert r1.subtracted_counts[0] == 90
    assert r2.subtracted_counts[0] == 190
    print(f"  ✓ test_consent_registry_multi_source_reuse")


def test_consent_registry_revoke():
    """revoke() removes an existing approval."""
    registry = BackgroundConsentRegistry()
    src = make_synthetic([100], 100, a1=1.0, source_path="/src.spe")
    bg = make_synthetic([10], 100, a1=1.0, source_path="/bg.spe")
    registry.approve(bg.source_path)
    assert registry.is_approved(bg.source_path)
    # Subtraction works
    subtract_background(src, bg, consent_registry=registry)
    # Revoke
    registry.revoke(bg.source_path)
    assert not registry.is_approved(bg.source_path)
    # Now raises again
    try:
        subtract_background(src, bg, consent_registry=registry)
        assert False, "Should raise after revoke"
    except BackgroundConsentRequired:
        pass
    print(f"  ✓ test_consent_registry_revoke")


def test_consent_registry_clear():
    """clear() removes ALL approvals."""
    registry = BackgroundConsentRegistry()
    registry.approve("/bg1.spe")
    registry.approve("/bg2.spe")
    assert len(registry.approved_paths()) == 2
    registry.clear()
    assert len(registry.approved_paths()) == 0
    print(f"  ✓ test_consent_registry_clear")


# --- Mode selection ---

def test_mode_select_rate_normalized_when_matched():
    """When gains match, use channel mode."""
    src = make_synthetic([100], 100, a1=1.0)
    bg = make_synthetic([10], 100, a1=1.0001)   # 0.01% mismatch
    result = subtract_background(src, bg, user_confirmed_applicable=True)
    assert result.mode == "rate_normalized_channel"
    print(f"  ✓ test_mode_select_rate_normalized_when_matched")


def test_mode_select_energy_aligned_when_mismatched():
    """When gains differ by > threshold, use energy_aligned mode."""
    src = make_synthetic([100]*10, 100, a1=1.0)
    bg = make_synthetic([10]*10, 100, a1=1.05)   # 5% mismatch
    result = subtract_background(src, bg, user_confirmed_applicable=True)
    assert result.mode == "energy_aligned"
    assert result.gain_mismatch_relative > GAIN_MATCH_THRESHOLD
    print(f"  ✓ test_mode_select_energy_aligned_when_mismatched")


# --- Rate scaling ---

def test_rate_scaling_short_bg():
    """If bg is shorter than source, bg counts scaled UP."""
    # source: 1000s, bg: 100s — scale = 10x
    src = make_synthetic([100], 1000, a1=1.0)
    bg = make_synthetic([10], 100, a1=1.0)
    result = subtract_background(src, bg, user_confirmed_applicable=True)
    assert result.rate_scale == 10.0
    # net = 100 - 10*10 = 0
    np.testing.assert_array_equal(result.subtracted_counts, [0])
    print(f"  ✓ test_rate_scaling_short_bg (scale=10×)")


def test_rate_scaling_long_bg():
    """If bg is longer than source, bg counts scaled DOWN."""
    # source: 100s, bg: 1000s — scale = 0.1x
    src = make_synthetic([100], 100, a1=1.0)
    bg = make_synthetic([200], 1000, a1=1.0)
    result = subtract_background(src, bg, user_confirmed_applicable=True)
    assert result.rate_scale == 0.1
    # net = 100 - 200*0.1 = 80
    np.testing.assert_array_equal(result.subtracted_counts, [80])
    print(f"  ✓ test_rate_scaling_long_bg (scale=0.1×)")


# --- Clipping ---

def test_negative_clipped_to_zero():
    """Channels where bg > source should clip to zero, not go negative."""
    src = make_synthetic([5, 100, 50], 100, a1=1.0)
    bg = make_synthetic([10, 5, 100], 100, a1=1.0)   # bg[0]>src[0], bg[2]>src[2]
    result = subtract_background(src, bg, user_confirmed_applicable=True)
    np.testing.assert_array_equal(result.subtracted_counts, [0, 95, 0])
    assert result.n_channels_clipped_to_zero == 2
    print(f"  ✓ test_negative_clipped_to_zero (2 channels clipped)")


# --- Uncertainty propagation ---

def test_uncertainty_propagation():
    """σ²(net) = N_src + N_bg · scale² per Poisson propagation."""
    # src=100 (σ=10), bg=25 (σ=5), scale=1, expected σ²(net) = 100 + 25 = 125
    src = make_synthetic([100], 100, a1=1.0)
    bg = make_synthetic([25], 100, a1=1.0)
    result = subtract_background(src, bg, user_confirmed_applicable=True)
    expected_unc = np.sqrt(100 + 25)  # ≈ 11.18
    assert abs(result.subtracted_uncertainties[0] - expected_unc) < 0.01
    print(f"  ✓ test_uncertainty_propagation (σ={result.subtracted_uncertainties[0]:.2f})")


# --- Energy-aligned mode ---

def test_energy_aligned_preserves_total_after_scale():
    """When bg shape is uniform and gain matches, energy-aligned and
    channel modes give the same result for interior channels (modulo
    ±1 count from integer rounding in fractional-bin integration)."""
    src = make_synthetic([100]*10, 100, a1=1.0)
    bg = make_synthetic([20]*10, 100, a1=1.0)
    r_ch = subtract_background(src, bg, user_confirmed_applicable=True,
                                force_mode="rate_normalized_channel")
    r_ea = subtract_background(src, bg, user_confirmed_applicable=True,
                                force_mode="energy_aligned")
    # Compare interior channels — allow ±1 count rounding difference
    diff = np.abs(r_ch.subtracted_counts[1:-1].astype(int) -
                  r_ea.subtracted_counts[1:-1].astype(int))
    assert diff.max() <= 1, f"Max diff > 1: {diff}"
    print(f"  ✓ test_energy_aligned_preserves_total_after_scale "
          f"(interior matches within ±1)")


def test_energy_aligned_handles_offset():
    """Energy-aligned mode correctly handles offset shift (bg has
    different a0 than source)."""
    # bg shifted by 1 channel worth of energy in offset
    # src: ch 0 -> E=0, ch 1 -> E=1, ch 2 -> E=2 ...
    # bg:  ch 0 -> E=1, ch 1 -> E=2, ch 2 -> E=3 ... (offset by +1 keV)
    # When we subtract from source bin E=[0.5, 1.5], that maps to bg bin E=[0.5, 1.5]
    # which is bg fractional channels [-0.5, 0.5] — clipped to [0, 0.5] = half of bg[0]
    src = make_synthetic([100, 100, 100, 100, 100], 100, a0=0, a1=1.0)
    bg = make_synthetic([20, 20, 20, 20, 20], 100, a0=1.0, a1=1.0)
    result = subtract_background(src, bg, user_confirmed_applicable=True,
                                  force_mode="energy_aligned")
    # All bg counts uniformly 20; offset by 1 keV but gain matches.
    # For source channels well within bg's range, integrated bg should
    # still be ~20 per source channel.
    # Channels in the middle should give ~80 (100 - 20).
    middle = result.subtracted_counts[2]
    assert abs(middle - 80) < 5, f"Middle channel: {middle}"
    print(f"  ✓ test_energy_aligned_handles_offset (middle={middle})")


# --- Real reference data ---

def test_subtract_gamma1s_cd109():
    """Cd-109 real source vs Gamma-1S background.
    Cd-109 source has cps=41, bg has cps=7.6; Δa1=3% (mismatch)."""
    bg = read_spectrum(f"{REF_DIR}/Фон_закр_кр_вода_01.spe")
    src = read_spectrum(f"{REF_DIR}/Cd-109__175_04_2017_Точечная-5см_5cm.spe")
    result = subtract_background(src, bg, user_confirmed_applicable=True)
    # Should auto-select energy_aligned due to 3% gain mismatch
    assert result.mode == "energy_aligned"
    # rate_scale = src.live_time / bg.live_time = 1800/3600 = 0.5
    assert abs(result.rate_scale - 0.5) < 0.001
    # Some channels should be clipped (bg dominates low E)
    assert result.n_channels_clipped_to_zero > 0
    # Total net counts should be positive (source dominates)
    assert result.subtracted_counts.sum() > 0
    # And less than original
    assert result.subtracted_counts.sum() < src.counts.sum()
    print(f"  ✓ test_subtract_gamma1s_cd109 "
          f"(mode={result.mode}, scale={result.rate_scale:.2f}, "
          f"clipped={result.n_channels_clipped_to_zero})")


def test_apply_subtraction_to_spectrum():
    """apply_subtraction_to_spectrum returns a new Spectrum with bg-subtracted counts
    and the background_subtracted flag set."""
    src = make_synthetic([100, 200, 300], 100, a1=1.0)
    bg = make_synthetic([10, 20, 30], 100, a1=1.0)
    result = subtract_background(src, bg, user_confirmed_applicable=True)
    sub_spec = apply_subtraction_to_spectrum(src, result)
    np.testing.assert_array_equal(sub_spec.counts, [90, 180, 270])
    assert sub_spec.extras.get("background_subtracted") is True
    assert sub_spec.live_time == src.live_time   # preserved
    assert sub_spec.energy_cal == src.energy_cal  # preserved
    print(f"  ✓ test_apply_subtraction_to_spectrum")


if __name__ == "__main__":
    print("Running Phase 2.1e background subtraction tests...\n")
    test_subtract_without_consent_raises()
    test_subtract_with_consent_works()
    test_consent_registry_empty_raises()
    test_consent_registry_after_approve_proceeds()
    test_consent_registry_per_file_isolation()
    test_consent_registry_multi_source_reuse()
    test_consent_registry_revoke()
    test_consent_registry_clear()
    test_mode_select_rate_normalized_when_matched()
    test_mode_select_energy_aligned_when_mismatched()
    test_rate_scaling_short_bg()
    test_rate_scaling_long_bg()
    test_negative_clipped_to_zero()
    test_uncertainty_propagation()
    test_energy_aligned_preserves_total_after_scale()
    test_energy_aligned_handles_offset()
    test_subtract_gamma1s_cd109()
    test_apply_subtraction_to_spectrum()
    print("\n✓ All Phase 2.1e background subtraction tests passed.")
