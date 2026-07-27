"""
Tests for `gamma.physics.self_attenuation` — K-20 v1.9.0 module.

Covers:
  1. NIST XCOM table integrity (mass-fraction sum, monotonic energy)
  2. Element μ/ρ interpolation (log-log) at pillar and intermediate energies
  3. ОИСН-16 matrix μ/ρ at reference cross-check energies
  4. Slab self-attenuation factor limits (F→1 at small μρt; F→1/μρt at large)
  5. Correction factor symmetry (corr(ρ_a→ρ_b) = 1/corr(ρ_b→ρ_a))
  6. Correction factor at reference density returns 1.0
  7. Weighted mean correction matches single-line when only one line
  8. v1.9.0 spread-reduction empirical claim: Cs-137 Marinelli
     light/heavy correction factor ratio matches cert-matrix Δ
     reduction within tolerance.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.physics.self_attenuation import (
    _XCOM_ENERGIES_KEV,
    _XCOM_MU_RHO,
    OISN_16_COMPOSITION,
    REF_GEOMETRY,
    element_mu_over_rho,
    matrix_mu_over_rho,
    slab_self_attenuation_factor,
    correction_factor,
    weighted_mean_correction,
)


# ---------------------------------------------------------------------------
# 1. NIST XCOM table integrity
# ---------------------------------------------------------------------------

def test_xcom_table_integrity_energy_monotonic():
    """Pillar energies must be strictly increasing."""
    for i in range(len(_XCOM_ENERGIES_KEV) - 1):
        assert _XCOM_ENERGIES_KEV[i] < _XCOM_ENERGIES_KEV[i + 1], (
            f"_XCOM_ENERGIES_KEV[{i}]={_XCOM_ENERGIES_KEV[i]} "
            f">= [{i+1}]={_XCOM_ENERGIES_KEV[i+1]}"
        )


def test_xcom_table_integrity_all_elements_same_length():
    """All elements must have the same number of μ/ρ values as energies."""
    n_E = len(_XCOM_ENERGIES_KEV)
    for sym, vals in _XCOM_MU_RHO.items():
        assert len(vals) == n_E, (
            f"element {sym!r}: {len(vals)} μ/ρ values vs {n_E} energies"
        )


def test_xcom_table_mu_rho_monotonic_in_energy_per_element():
    """μ/ρ decreases monotonically with energy for each element
    above the K-edge (no K-edge in our 50-3000 keV range for any of
    H/C/N/O/Fe — Fe K-edge is at 7.1 keV)."""
    for sym, vals in _XCOM_MU_RHO.items():
        for i in range(len(vals) - 1):
            assert vals[i] >= vals[i + 1], (
                f"element {sym!r}: μ/ρ not monotone decreasing at index "
                f"{i}: {vals[i]} < {vals[i+1]}"
            )


# ---------------------------------------------------------------------------
# 2. Element μ/ρ interpolation
# ---------------------------------------------------------------------------

def test_element_mu_rho_at_pillar_energy_exact():
    """At pillar energies, interpolation should return the table value exactly."""
    for sym, vals in _XCOM_MU_RHO.items():
        for E, mu in zip(_XCOM_ENERGIES_KEV, vals):
            got = element_mu_over_rho(E, sym)
            assert abs(got - mu) < 1e-9, (
                f"element {sym!r} at E={E} keV: got {got}, expected {mu}"
            )


def test_element_mu_rho_below_table_clamps_to_first():
    """Below the lowest pillar (50 keV), value clamps to the table edge."""
    for sym, vals in _XCOM_MU_RHO.items():
        got = element_mu_over_rho(20.0, sym)
        assert got == vals[0]


def test_element_mu_rho_above_table_clamps_to_last():
    """Above the highest pillar (3000 keV), value clamps to the table edge."""
    for sym, vals in _XCOM_MU_RHO.items():
        got = element_mu_over_rho(5000.0, sym)
        assert got == vals[-1]


def test_element_mu_rho_log_log_interpolation_at_geometric_mean():
    """At the geometric mean of two pillars E1, E2, log-log linear
    interpolation gives μ/ρ = sqrt(μ1 × μ2)."""
    # Use H at the 100→150 pillar interval; geo mean is 122.47
    E_geo = math.sqrt(100.0 * 150.0)
    got = element_mu_over_rho(E_geo, "H")
    expected = math.sqrt(_XCOM_MU_RHO["H"][3] * _XCOM_MU_RHO["H"][4])
    rel_err = abs(got - expected) / expected
    assert rel_err < 1e-9, (
        f"geometric-mean log-log interpolation broken: got {got}, "
        f"expected {expected}, rel_err {rel_err}"
    )


def test_element_mu_rho_raises_on_unknown_element():
    try:
        element_mu_over_rho(100.0, "Pb")
    except KeyError as e:
        assert "Pb" in str(e)
        return
    assert False, "expected KeyError for unknown element"


def test_element_mu_rho_raises_on_non_positive_energy():
    for bad_E in (0.0, -50.0):
        try:
            element_mu_over_rho(bad_E, "H")
        except ValueError:
            continue
        assert False, f"expected ValueError for E_keV={bad_E}"


# ---------------------------------------------------------------------------
# 3. ОИСН-16 matrix μ/ρ
# ---------------------------------------------------------------------------

def test_oisn16_composition_sums_to_one():
    total = sum(OISN_16_COMPOSITION.values())
    assert abs(total - 1.0) < 1e-9, (
        f"ОИСН-16 mass fractions sum to {total}, expected 1.0"
    )


def test_oisn16_matrix_mu_rho_at_600_keV():
    """Hand-computed reference cross-check at 600 keV.

    OISN-16 = H 0.022 + C 0.206 + N 0.009 + O 0.049 + Fe 0.714 (mass).
    At 600 keV pillar:
      H  μ/ρ = 0.1271
      C  μ/ρ = 0.0586
      N  μ/ρ = 0.0573
      O  μ/ρ = 0.0577
      Fe μ/ρ = 0.0769
    Matrix μ/ρ = sum(w·μ/ρ_i)
              = 0.022×0.1271 + 0.206×0.0586 + 0.009×0.0573
                + 0.049×0.0577 + 0.714×0.0769
              = 0.0027962 + 0.0120716 + 0.00051570 + 0.00282730 + 0.05490660
              ≈ 0.07312 cm²/g
    """
    got = matrix_mu_over_rho(600.0, OISN_16_COMPOSITION)
    expected = (0.022 * 0.1271 + 0.206 * 0.0586 + 0.009 * 0.0573
                + 0.049 * 0.0577 + 0.714 * 0.0769)
    assert abs(got - expected) < 1e-6, (
        f"ОИСН-16 μ/ρ at 600 keV: got {got}, expected {expected}"
    )


def test_oisn16_matrix_mu_rho_monotonic_in_energy():
    """For ОИСН-16 in 50-3000 keV (no K-edge), μ/ρ decreases with E."""
    energies = [60, 100, 200, 500, 1000, 2000, 3000]
    mu_rhos = [matrix_mu_over_rho(E, OISN_16_COMPOSITION) for E in energies]
    for i in range(len(mu_rhos) - 1):
        assert mu_rhos[i] > mu_rhos[i + 1], (
            f"ОИСН-16 μ/ρ not strictly decreasing at "
            f"E={energies[i]}→{energies[i+1]}: "
            f"{mu_rhos[i]} → {mu_rhos[i+1]}"
        )


# ---------------------------------------------------------------------------
# 4. Slab self-attenuation factor
# ---------------------------------------------------------------------------

def test_slab_F_equals_1_at_zero_thickness():
    F = slab_self_attenuation_factor(0.1, 1.0, 0.0)
    assert abs(F - 1.0) < 1e-9


def test_slab_F_equals_1_at_zero_density():
    F = slab_self_attenuation_factor(0.1, 0.0, 3.0)
    assert abs(F - 1.0) < 1e-9


def test_slab_F_equals_1_at_zero_mu():
    F = slab_self_attenuation_factor(0.0, 1.0, 3.0)
    assert abs(F - 1.0) < 1e-9


def test_slab_F_decreases_with_thickness():
    """For fixed (μ, ρ), F monotonically decreases as t increases."""
    mu = 0.1
    rho = 1.0
    F_values = [slab_self_attenuation_factor(mu, rho, t)
                for t in (0.1, 0.5, 1.0, 5.0, 10.0)]
    for i in range(len(F_values) - 1):
        assert F_values[i] > F_values[i + 1], (
            f"F not monotone decreasing at t-index {i}: "
            f"{F_values[i]} → {F_values[i+1]}"
        )


def test_slab_F_asymptotic_thick_limit():
    """In the μρt → ∞ limit, F → 1/(μρt)."""
    mu = 1.0
    rho = 1.0
    t = 100.0
    F = slab_self_attenuation_factor(mu, rho, t)
    expected = 1.0 / (mu * rho * t)  # exp(-100) is negligible
    assert abs(F - expected) < 1e-6


def test_slab_F_raises_on_negative_density_or_thickness():
    for rho, t in [(-1.0, 3.0), (1.0, -3.0)]:
        try:
            slab_self_attenuation_factor(0.1, rho, t)
        except ValueError:
            continue
        assert False, f"expected ValueError for ρ={rho}, t={t}"


# ---------------------------------------------------------------------------
# 5. Correction factor
# ---------------------------------------------------------------------------

def test_correction_factor_equals_1_at_reference_density():
    """When ρ_sample == ρ_ref, correction = 1 exactly."""
    c = correction_factor(662.0, rho_sample_g_cm3=1.60, rho_ref_g_cm3=1.60,
                          thickness_cm=3.1)
    assert abs(c - 1.0) < 1e-9


def test_correction_factor_below_1_for_light_sample():
    """ρ_sample < ρ_ref → F_sample > F_ref → corr < 1 (reduces over-est)."""
    c = correction_factor(662.0, rho_sample_g_cm3=0.5, rho_ref_g_cm3=1.6,
                          thickness_cm=3.1)
    assert c < 1.0
    assert c > 0.5  # but reasonable


def test_correction_factor_above_1_for_heavy_sample():
    """ρ_sample > ρ_ref → F_sample < F_ref → corr > 1 (boosts under-est)."""
    c = correction_factor(662.0, rho_sample_g_cm3=3.0, rho_ref_g_cm3=1.6,
                          thickness_cm=3.1)
    assert c > 1.0
    assert c < 2.0  # but reasonable


def test_correction_factor_symmetry_inverse():
    """corr(ρa→ρb) × corr(ρb→ρa) = 1 (algebraic identity)."""
    c_ab = correction_factor(662.0, rho_sample_g_cm3=0.5, rho_ref_g_cm3=1.6,
                             thickness_cm=3.1)
    c_ba = correction_factor(662.0, rho_sample_g_cm3=1.6, rho_ref_g_cm3=0.5,
                             thickness_cm=3.1)
    product = c_ab * c_ba
    assert abs(product - 1.0) < 1e-9, (
        f"correction-factor symmetry broken: {c_ab} × {c_ba} = {product}"
    )


def test_correction_factor_returns_1_for_zero_reference_density():
    """Point sources (no matrix) have ρ_ref=0 → no correction defined."""
    c = correction_factor(662.0, rho_sample_g_cm3=0.5, rho_ref_g_cm3=0.0,
                          thickness_cm=0.0)
    assert c == 1.0


def test_correction_factor_higher_energy_smaller_effect():
    """μ/ρ decreases with E, so correction factor → 1 as E increases."""
    c_low = correction_factor(100.0, rho_sample_g_cm3=0.5,
                              rho_ref_g_cm3=1.6, thickness_cm=3.1)
    c_high = correction_factor(2614.0, rho_sample_g_cm3=0.5,
                               rho_ref_g_cm3=1.6, thickness_cm=3.1)
    # |c-1| should be smaller at high E (less attenuation differential)
    assert abs(c_high - 1.0) < abs(c_low - 1.0), (
        f"high-E correction should be closer to 1: "
        f"|c(100)-1|={abs(c_low-1):.4f}, |c(2614)-1|={abs(c_high-1):.4f}"
    )


# ---------------------------------------------------------------------------
# 6. Weighted mean correction
# ---------------------------------------------------------------------------

def test_weighted_mean_single_line_equals_single_factor():
    """With one line and any weight, weighted mean equals single corr."""
    E = 662.0
    single = correction_factor(E, rho_sample_g_cm3=0.5,
                               rho_ref_g_cm3=1.6, thickness_cm=3.1)
    wmean = weighted_mean_correction([E], [1.0],
                                     rho_sample_g_cm3=0.5,
                                     rho_ref_g_cm3=1.6, thickness_cm=3.1)
    assert abs(single - wmean) < 1e-9


def test_weighted_mean_equal_weights_equals_arithmetic_mean():
    """With equal weights, weighted mean reduces to arithmetic mean."""
    Es = [583.0, 1461.0, 2614.0]
    cs = [correction_factor(E, rho_sample_g_cm3=0.5, rho_ref_g_cm3=1.6,
                            thickness_cm=3.1)
          for E in Es]
    arith = sum(cs) / len(cs)
    wmean = weighted_mean_correction(Es, [1.0] * 3,
                                     rho_sample_g_cm3=0.5,
                                     rho_ref_g_cm3=1.6, thickness_cm=3.1)
    assert abs(arith - wmean) < 1e-9


def test_weighted_mean_returns_1_for_empty_input():
    wmean = weighted_mean_correction([], [],
                                     rho_sample_g_cm3=0.5,
                                     rho_ref_g_cm3=1.6, thickness_cm=3.1)
    assert wmean == 1.0


def test_weighted_mean_raises_on_length_mismatch():
    try:
        weighted_mean_correction([100.0, 200.0], [1.0],
                                 rho_sample_g_cm3=0.5,
                                 rho_ref_g_cm3=1.6, thickness_cm=3.1)
    except ValueError:
        return
    assert False, "expected ValueError on length mismatch"


# ---------------------------------------------------------------------------
# 7. REF_GEOMETRY registry
# ---------------------------------------------------------------------------

def test_ref_geometry_marinelli_present():
    """К-20 only applies to Marinelli geometry (Layers.Enable=false in
    its .efr). Дента-120мл and Петри-60мл .efr files have
    Layers.Enable=true (matrix correction baked in), so they are
    intentionally NOT in REF_GEOMETRY."""
    assert "Маринелли" in REF_GEOMETRY


def test_ref_geometry_excludes_denta_petri():
    """Дента-120мл and Петри-60мл deliberately excluded — see K-20
    note in self_attenuation.py header about Layers.Enable=true."""
    assert "Дента-120мл" not in REF_GEOMETRY
    assert "Петри-60мл" not in REF_GEOMETRY


def test_ref_geometry_marinelli_values():
    """Verified from .efr metadata 2024-11."""
    vol, rho_ref, t_cm = REF_GEOMETRY["Маринелли"]
    assert vol == 1000.0
    assert rho_ref == 1.60
    assert abs(t_cm - 3.1) < 1e-9


# ---------------------------------------------------------------------------
# 8. v1.9.0 empirical claim — Cs-137 Marinelli spread reduction
# ---------------------------------------------------------------------------

def test_marinelli_cs137_correction_ratio_matches_cert_finding():
    """
    Empirical v1.7.25 cert-matrix Δ for Cs-137 in Marinelli:
        Light (ρ_sample=0.57 ≈ ratio 0.36): Δ = +9.81 %
        Heavy (ρ_sample=1.66 ≈ ratio 1.04): Δ = −5.58 %
        Peak-to-peak spread = 15.4 %

    K-20 correction factors at 662 keV:
        Light: ρ=0.57, ρ_ref=1.60, t=3.1 → corr ≈ 0.89
        Heavy: ρ=1.66, ρ_ref=1.60, t=3.1 → corr ≈ 1.006

    Post-correction Δ predictions:
        Light: 1.0981 × 0.89 − 1 = −2.3 % (was +9.81 %)
        Heavy: 0.9442 × 1.006 − 1 = −5.0 % (was −5.58 %)
        Post spread ≈ 2.7 % (vs 15.4 % before)

    Asserts: (a) light correction < 1; (b) light correction reduces
    the absolute deviation; (c) post-correction spread ≤ 5 %
    (matching v1.9.0 K-20 design target).
    """
    vol, rho_ref, t_cm = REF_GEOMETRY["Маринелли"]
    # Light source (Cs137_420-7-14, 570 g in 1000 ml)
    rho_light = 570.0 / vol
    # Heavy source (Cs137_420-7-15, 1660 g in 1000 ml)
    rho_heavy = 1660.0 / vol
    c_light = correction_factor(662.0, rho_sample_g_cm3=rho_light,
                                rho_ref_g_cm3=rho_ref, thickness_cm=t_cm)
    c_heavy = correction_factor(662.0, rho_sample_g_cm3=rho_heavy,
                                rho_ref_g_cm3=rho_ref, thickness_cm=t_cm)
    # Empirical cert-matrix Δ (v1.7.25)
    delta_light_before = +0.0981
    delta_heavy_before = -0.0558
    delta_light_after = (1 + delta_light_before) * c_light - 1
    delta_heavy_after = (1 + delta_heavy_before) * c_heavy - 1
    spread_before = abs(delta_light_before - delta_heavy_before)
    spread_after = abs(delta_light_after - delta_heavy_after)
    print(f"  K-20 Marinelli Cs-137 spread: "
          f"{spread_before*100:.1f}% → {spread_after*100:.1f}%")
    # (a) light correction reduces overestimate
    assert c_light < 1.0, (
        f"K-20 light correction must be < 1, got {c_light}"
    )
    # (b) light correction reduces |Δ|
    assert abs(delta_light_after) < abs(delta_light_before), (
        f"K-20 light correction increased |Δ|: "
        f"{abs(delta_light_before)*100}% → "
        f"{abs(delta_light_after)*100}%"
    )
    # (c) spread reduced to ≤ 5 % (v1.9.0 design target)
    assert spread_after < 0.05, (
        f"K-20 spread reduction target violated: "
        f"{spread_after*100:.2f}% > 5 %"
    )


# ---------------------------------------------------------------------------
# Test runner
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    tests = [
        test_xcom_table_integrity_energy_monotonic,
        test_xcom_table_integrity_all_elements_same_length,
        test_xcom_table_mu_rho_monotonic_in_energy_per_element,
        test_element_mu_rho_at_pillar_energy_exact,
        test_element_mu_rho_below_table_clamps_to_first,
        test_element_mu_rho_above_table_clamps_to_last,
        test_element_mu_rho_log_log_interpolation_at_geometric_mean,
        test_element_mu_rho_raises_on_unknown_element,
        test_element_mu_rho_raises_on_non_positive_energy,
        test_oisn16_composition_sums_to_one,
        test_oisn16_matrix_mu_rho_at_600_keV,
        test_oisn16_matrix_mu_rho_monotonic_in_energy,
        test_slab_F_equals_1_at_zero_thickness,
        test_slab_F_equals_1_at_zero_density,
        test_slab_F_equals_1_at_zero_mu,
        test_slab_F_decreases_with_thickness,
        test_slab_F_asymptotic_thick_limit,
        test_slab_F_raises_on_negative_density_or_thickness,
        test_correction_factor_equals_1_at_reference_density,
        test_correction_factor_below_1_for_light_sample,
        test_correction_factor_above_1_for_heavy_sample,
        test_correction_factor_symmetry_inverse,
        test_correction_factor_returns_1_for_zero_reference_density,
        test_correction_factor_higher_energy_smaller_effect,
        test_weighted_mean_single_line_equals_single_factor,
        test_weighted_mean_equal_weights_equals_arithmetic_mean,
        test_weighted_mean_returns_1_for_empty_input,
        test_weighted_mean_raises_on_length_mismatch,
        test_ref_geometry_marinelli_present,
        test_ref_geometry_excludes_denta_petri,
        test_ref_geometry_marinelli_values,
        test_marinelli_cs137_correction_ratio_matches_cert_finding,
    ]
    passed = 0
    failed = []
    for t in tests:
        try:
            print(f"-- {t.__name__}")
            t()
            print(f"   OK")
            passed += 1
        except Exception as e:
            print(f"   FAIL: {e}")
            failed.append((t.__name__, str(e)))
    print()
    print(f"Passed: {passed}/{len(tests)}")
    if failed:
        for name, err in failed:
            print(f"  - {name}: {err}")
        sys.exit(1)
    sys.exit(0)
