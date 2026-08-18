"""#PTB-3 (2026-07-02) — volume-source TCS averaging, PTB-2018 Annex D Eq. (D17).

Verifies scripts/gamma/activity/tcs_close_geometry.py::compute_tcs_correction_volume
against the point-source compute_tcs_correction and hand-derived two-voxel cases.

Discretization under test:
    C_TCS = 1 / (1 - sum_j p_ij * eps_T(E_j) * <s_p*s_t>/<s_p>)
where <.> is the volume-weighted mean over VolumeElement entries.
"""
import pytest

from gamma.activity.tcs_close_geometry import (
    VolumeElement,
    compute_tcs_correction,
    compute_tcs_correction_volume,
    CO60_PAIRS,
)


def _eps_t_010(E_keV):
    return 0.10


def test_volume_all_unit_scales_equals_point():
    """All s_p = s_t = 1 -> D17 average degenerates to the point formula."""
    point = compute_tcs_correction(
        E_i_keV=1173.2, nuclide_pairs=CO60_PAIRS,
        total_efficiency_func=_eps_t_010,
    )
    elements = [VolumeElement(weight=w, fep_efficiency_scale=1.0,
                              total_efficiency_scale=1.0)
                for w in (0.2, 0.5, 0.3)]
    vol = compute_tcs_correction_volume(
        E_i_keV=1173.2, nuclide_pairs=CO60_PAIRS,
        total_efficiency_func=_eps_t_010, volume_elements=elements,
    )
    assert vol.correction_factor == pytest.approx(point.correction_factor, rel=1e-12)
    assert vol.sum_L_ij == pytest.approx(point.sum_L_ij, rel=1e-12)
    assert vol.n_pairs_used == point.n_pairs_used


def test_volume_two_voxel_hand_derived():
    """Two voxels: near (s_p=1, s_t=1) and far (s_p=0.5, s_t=0.5), w=1 each.

    eff_factor = (1*1*1 + 1*0.5*0.5) / (1*1 + 1*0.5) = 1.25/1.5 = 5/6.
    Co-60 at 1173.2 keV, eps_T = 0.10: point sum_L = 0.998*0.10 = 0.0998.
    volume sum_L = 0.0998 * 5/6 = 0.08316666...
    C = 1/(1 - 0.0831667) = 1.0907...
    """
    elements = [
        VolumeElement(weight=1.0, fep_efficiency_scale=1.0,
                      total_efficiency_scale=1.0),
        VolumeElement(weight=1.0, fep_efficiency_scale=0.5,
                      total_efficiency_scale=0.5),
    ]
    vol = compute_tcs_correction_volume(
        E_i_keV=1173.2, nuclide_pairs=CO60_PAIRS,
        total_efficiency_func=_eps_t_010, volume_elements=elements,
    )
    expected_sum_l = 0.998 * 0.10 * (1.25 / 1.5)
    assert vol.sum_L_ij == pytest.approx(expected_sum_l, rel=1e-12)
    assert vol.correction_factor == pytest.approx(
        1.0 / (1.0 - expected_sum_l), rel=1e-12
    )


def test_volume_far_elements_reduce_correction():
    """Distant voxels (s_t < 1) -> less summing-out than the point estimate."""
    point = compute_tcs_correction(
        E_i_keV=1173.2, nuclide_pairs=CO60_PAIRS,
        total_efficiency_func=_eps_t_010,
    )
    elements = [
        VolumeElement(weight=1.0, fep_efficiency_scale=1.0,
                      total_efficiency_scale=1.0),
        VolumeElement(weight=2.0, fep_efficiency_scale=0.4,
                      total_efficiency_scale=0.3),
    ]
    vol = compute_tcs_correction_volume(
        E_i_keV=1173.2, nuclide_pairs=CO60_PAIRS,
        total_efficiency_func=_eps_t_010, volume_elements=elements,
    )
    assert vol.correction_factor < point.correction_factor
    assert vol.correction_factor > 1.0


def test_volume_weight_scale_invariance():
    """Multiplying all weights by a constant must not change the result."""
    base = [
        VolumeElement(1.0, 1.0, 1.0),
        VolumeElement(2.0, 0.6, 0.5),
    ]
    scaled = [
        VolumeElement(7.0, 1.0, 1.0),
        VolumeElement(14.0, 0.6, 0.5),
    ]
    a = compute_tcs_correction_volume(
        1173.2, CO60_PAIRS, _eps_t_010, volume_elements=base)
    b = compute_tcs_correction_volume(
        1173.2, CO60_PAIRS, _eps_t_010, volume_elements=scaled)
    assert a.correction_factor == pytest.approx(b.correction_factor, rel=1e-12)


def test_volume_guards_raise():
    """Empty elements, negative weight/scale, zero fep mass -> ValueError."""
    with pytest.raises(ValueError):
        compute_tcs_correction_volume(
            1173.2, CO60_PAIRS, _eps_t_010, volume_elements=[])
    with pytest.raises(ValueError):
        compute_tcs_correction_volume(
            1173.2, CO60_PAIRS, _eps_t_010,
            volume_elements=[VolumeElement(-1.0, 1.0, 1.0)])
    with pytest.raises(ValueError):
        compute_tcs_correction_volume(
            1173.2, CO60_PAIRS, _eps_t_010,
            volume_elements=[VolumeElement(1.0, -0.1, 1.0)])
    with pytest.raises(ValueError):
        compute_tcs_correction_volume(
            1173.2, CO60_PAIRS, _eps_t_010,
            volume_elements=[VolumeElement(1.0, 0.0, 1.0)])