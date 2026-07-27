"""
Snapshot tests for the 4 diagnostic fields added to BackgroundSubtractionResult in Phase 2 of the F-58 LITE/FULL untangle.
"""

from __future__ import annotations
import math
import numpy as np
import pytest
from gamma.spectrum import Spectrum
from gamma.calibration.bg_subtract_energy import subtract_background, BackgroundSubtractionResult


def _make_spec(counts, live_time, *, a0=0.0, a1=1.0, source_path="test"):
    spec = Spectrum(
        counts=np.array(counts, dtype=np.int64),
        live_time=float(live_time),
        real_time=float(live_time),
        source_path=source_path,
        source_format="test",
    )
    spec.energy_cal = (float(a0), float(a1))
    spec.n_channels = len(counts)
    spec.n_channels_raw = len(counts)
    return spec


def test_n_channels_clipped_counts_negative_diff():
    """Verify n_channels_clipped counts negative differences."""
    sample = _make_spec([10, 10, 10, 10, 10], 100, a0=0, a1=1)
    bg = _make_spec([5, 5, 20, 25, 30], 100, a0=0, a1=1)
    result = subtract_background(sample, bg)
    assert result.n_channels_clipped == 3


def test_n_channels_clipped_zero_when_sample_dominates():
    """Verify n_channels_clipped is zero when sample dominates."""
    sample = _make_spec([100] * 5, 100, a0=0, a1=1)
    bg = _make_spec([1] * 5, 100, a0=0, a1=1)
    result = subtract_background(sample, bg)
    assert result.n_channels_clipped == 0


def test_net_uncertainties_poisson_propagation_equal_live_times():
    """Verify net uncertainties with poisson propagation for equal live times."""
    sample = _make_spec([100] * 4, 100, a0=0, a1=1)
    bg = _make_spec([25] * 4, 100, a0=0, a1=1)
    result = subtract_background(sample, bg)
    assert result.net_uncertainties is not None
    assert result.net_uncertainties.shape == (4,)
    expected = np.sqrt(125.0)
    for val in result.net_uncertainties:
        assert val == pytest.approx(expected, rel=1e-9)


def test_net_uncertainties_with_scale_factor_k_neq_1():
    """Verify net uncertainties with scale factor k != 1 (F-451).

    sample (200 counts in 200 s) vs background (50 counts in 100 s).
    Legacy ratio t_s/t_bg = 2.0 is exposed as `scale_factor`.

    F-451 direction: t_s > t_bg → sample is scaled DOWN by
    applied_scale = t_bg/t_s = 0.5. In the effective (t_bg=100 s) scale:
        sample_eff = 200·0.5 = 100, bg_eff = 50, net = 50.
        σ²(net) = 0.5²·N_sample_raw + N_bg_raw = 0.25·200 + 50 = 100
        σ = 10.0   (legacy direction would have given 20.0; the
        cps-invariant 0.5 cps still holds: 50/100 s = 100/200 s).
    """
    sample = _make_spec([200] * 4, 200, a0=0, a1=1)
    bg = _make_spec([50] * 4, 100, a0=0, a1=1)
    result = subtract_background(sample, bg)
    assert result.scale_factor == pytest.approx(2.0, rel=1e-12)
    assert result.applied_scale == pytest.approx(0.5, rel=1e-12)
    assert result.scale_direction == "sample_down"
    assert result.effective_live_time == pytest.approx(100.0, rel=1e-12)
    expected = 10.0
    for val in result.net_uncertainties:
        assert val == pytest.approx(expected, rel=1e-9)


def test_gain_mismatch_relative_known_a1():
    """Verify gain mismatch relative when a1 is known."""
    sample = _make_spec([10] * 8, 100, a0=0.0, a1=3.0)
    bg = _make_spec([10] * 8, 100, a0=0.0, a1=2.97)
    result = subtract_background(sample, bg)
    expected = 0.01
    assert result.gain_mismatch_relative == pytest.approx(expected, rel=1e-9)


def test_gain_mismatch_relative_nan_when_a1_missing():
    """Verify gain mismatch relative is NaN when a1 is missing."""
    sample = _make_spec([10] * 8, 100, a0=0.0, a1=1.0)
    bg = _make_spec([10] * 8, 100, a0=0.0, a1=1.0)
    # Mutate bg to remove a1
    bg.energy_cal = (0.0,)
    result = subtract_background(sample, bg)
    assert math.isnan(result.gain_mismatch_relative)


def test_zero_point_mismatch_keV_known_a0():
    """Verify zero point mismatch keV when a0 is known and below threshold."""
    sample = _make_spec([10] * 8, 100, a0=0.0, a1=3.0)
    bg = _make_spec([10] * 8, 100, a0=4.5, a1=3.0)
    result = subtract_background(sample, bg)
    assert result.zero_point_mismatch_keV == pytest.approx(4.5, rel=1e-9)
    assert "F-243" not in result.notes


def test_notes_extension_above_5_keV_threshold():
    """Verify notes extension when zero point mismatch exceeds 5 keV."""
    sample = _make_spec([10] * 8, 100, a0=0.0, a1=3.0)
    bg = _make_spec([10] * 8, 100, a0=10.0, a1=3.0)
    result = subtract_background(sample, bg)
    assert result.zero_point_mismatch_keV == pytest.approx(10.0, rel=1e-9)
    assert "F-243" in result.notes
    assert "bg_subtract_dual_mode.py" in result.notes
    assert "keV" in result.notes


def test_clamp_negative_false_preserves_negative_net_but_count_unchanged():
    """Verify that clamp=False preserves negatives but counts are still tracked."""
    sample = _make_spec([10, 10, 10, 10], 100, a0=0, a1=1)
    bg = _make_spec([1, 50, 50, 1], 100, a0=0, a1=1)
    result = subtract_background(sample, bg, clamp_negative_to_zero=False)
    assert result.n_channels_clipped == 2
    assert any(result.net_counts < 0)
    # Now call with default (clamp=True)
    result2 = subtract_background(sample, bg, clamp_negative_to_zero=True)
    assert result2.n_channels_clipped == 2
    assert np.min(result2.net_counts) >= 0.0


# F-451: scale-direction inversion ("scale to smaller" / "к меньшему")

def test_f451_sample_down_when_sample_longer():
    """t_s > t_bg → sample is scaled DOWN to t_bg's scale."""
    sample = _make_spec([400] * 4, 400, a0=0, a1=1)
    bg = _make_spec([60] * 4, 100, a0=0, a1=1)
    result = subtract_background(sample, bg)
    assert result.scale_direction == "sample_down"
    assert result.applied_scale == pytest.approx(0.25, rel=1e-12)
    assert result.effective_live_time == pytest.approx(100.0, rel=1e-12)
    # net = 400·0.25 − 60 = 40 in effective (100 s) scale
    for v in result.net_counts:
        assert v == pytest.approx(40.0, rel=1e-9)
    # legacy ratio kept on scale_factor:
    assert result.scale_factor == pytest.approx(4.0, rel=1e-12)


def test_f451_bg_down_when_bg_longer():
    """t_s < t_bg → background is scaled DOWN to t_s's scale."""
    sample = _make_spec([100] * 4, 100, a0=0, a1=1)
    bg = _make_spec([200] * 4, 400, a0=0, a1=1)
    result = subtract_background(sample, bg)
    assert result.scale_direction == "bg_down"
    assert result.applied_scale == pytest.approx(0.25, rel=1e-12)
    assert result.effective_live_time == pytest.approx(100.0, rel=1e-12)
    # net = 100 − 200·0.25 = 50 in effective (100 s) scale
    for v in result.net_counts:
        assert v == pytest.approx(50.0, rel=1e-9)
    assert result.scale_factor == pytest.approx(0.25, rel=1e-12)


def test_f451_equal_live_time_no_scaling():
    """t_s == t_bg → equal direction, applied_scale=1, no scaling."""
    sample = _make_spec([100] * 4, 100, a0=0, a1=1)
    bg = _make_spec([25] * 4, 100, a0=0, a1=1)
    result = subtract_background(sample, bg)
    assert result.scale_direction == "equal"
    assert result.applied_scale == pytest.approx(1.0, rel=1e-12)
    assert result.effective_live_time == pytest.approx(100.0, rel=1e-12)
    for v in result.net_counts:
        assert v == pytest.approx(75.0, rel=1e-9)


def test_f451_cps_invariant_sample_down():
    """net_cps = sample_rate − bg_rate must hold under sample_down."""
    sample = _make_spec([400] * 4, 400, a0=0, a1=1)
    bg = _make_spec([60] * 4, 100, a0=0, a1=1)
    result = subtract_background(sample, bg)
    expected_cps = (400 / 400) - (60 / 100)  # 1.0 − 0.6 = 0.4
    for v in result.net_cps:
        assert v == pytest.approx(expected_cps, rel=1e-9)


def test_f451_cps_invariant_bg_down():
    """net_cps = sample_rate − bg_rate must hold under bg_down."""
    sample = _make_spec([100] * 4, 100, a0=0, a1=1)
    bg = _make_spec([200] * 4, 400, a0=0, a1=1)
    result = subtract_background(sample, bg)
    expected_cps = (100 / 100) - (200 / 400)  # 1.0 − 0.5 = 0.5
    for v in result.net_cps:
        assert v == pytest.approx(expected_cps, rel=1e-9)


def test_f451_notes_contain_direction_marker():
    """notes string carries F-451 + scale_direction + applied for downstream audit."""
    sample = _make_spec([400] * 4, 400, a0=0, a1=1)
    bg = _make_spec([60] * 4, 100, a0=0, a1=1)
    result = subtract_background(sample, bg)
    assert "F-451" in result.notes
    assert "sample_down" in result.notes
    assert "applied=0.2500" in result.notes