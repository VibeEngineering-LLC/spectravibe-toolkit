"""
Tests for Phase 2.1g cascade-summing correction (F-31b, v1.7.9):
  • Cascade scheme catalogue data is internally consistent
  • Peak-to-total NaI model produces sensible values at anchor energies
  • Total efficiency formula ε_T = ε_p / P
  • TCS correction factor for Co-60 at point 5cm ≈ 2-3% (Knoll §17.6)
  • Non-cascade nuclides (Cs-137) return correction = 1.0 / empty dict
  • Energy tolerance for line matching
  • Loss cap behaviour (safety against unphysical extrapolation)
  • HPGe vs NaI P/T comparison
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.physics.cascade_summing import (
    CASCADE_SCHEMES,
    CascadeScheme,
    CascadePartner,
    peak_to_total_NaI,
    peak_to_total_HPGe,
    total_efficiency,
    tcs_correction_factor,
    compute_tcs_corrections,
)


# ═════════════════════════════════════════════════════════════════════
# Stub efficiency curve for unit tests independent of real data
# ═════════════════════════════════════════════════════════════════════

class StubEfficiency:
    """Stand-in for EfficiencyCurve — implements `.efficiency_at(E)`
    and __call__ via a simple power law ε(E) = a · E^b.

    Default coefficients chosen to roughly match a NaI 63×63 mm point
    5cm geometry (ε(661)≈1.4e-2, ε(1332)≈6e-3).
    """

    def __init__(self, a=12.6, b=-1.05):
        self.a = a
        self.b = b

    def efficiency_at(self, E_keV: float) -> float:
        if E_keV <= 0:
            return 0.0
        return self.a * (E_keV ** self.b)

    def is_extrapolating(self, E_keV: float, margin_factor: float = 1.1) -> bool:
        # Treat the stub as calibrated across the whole γ range used in
        # tests (50–3000 keV) so compute_activity never flags
        # extrapolation for the synthetic fixture.
        return not (50.0 <= E_keV <= 3000.0)

    def __call__(self, E_keV: float) -> float:
        return self.efficiency_at(E_keV)


class ZeroEfficiency:
    """ε returns 0 everywhere — used to verify edge-case handling."""

    def efficiency_at(self, E_keV: float) -> float:
        return 0.0


# ═════════════════════════════════════════════════════════════════════
# Group 1 — cascade scheme catalogue sanity
# ═════════════════════════════════════════════════════════════════════

def test_catalogue_contains_major_cascade_nuclides():
    """Co-60, Y-88, Na-22, Tl-208, Eu-152, Ba-133 must all be present."""
    required = {"Co-60", "Y-88", "Na-22", "Tl-208", "Eu-152", "Ba-133"}
    missing = required - set(CASCADE_SCHEMES.keys())
    assert not missing, f"Missing cascade schemes: {missing}"
    print(f"  ✓ test_catalogue_contains_major_cascade_nuclides "
          f"({len(CASCADE_SCHEMES)} nuclides)")


def test_cascade_scheme_data_well_formed():
    """Every CascadeScheme has cascades dict mapping E→tuple of partners
    with positive energies and probabilities."""
    for nuc, scheme in CASCADE_SCHEMES.items():
        assert isinstance(scheme, CascadeScheme)
        assert scheme.nuclide == nuc
        assert isinstance(scheme.cascades, dict)
        assert len(scheme.cascades) >= 1
        for E_primary, partners in scheme.cascades.items():
            assert E_primary > 0, f"{nuc}: non-positive primary E={E_primary}"
            assert isinstance(partners, tuple)
            assert len(partners) >= 1
            for p in partners:
                assert isinstance(p, CascadePartner)
                assert p.E_partner_keV > 0
                assert 0.0 < p.probability <= 2.01, \
                    f"{nuc} {E_primary}: bad probability {p.probability}"
    print(f"  ✓ test_cascade_scheme_data_well_formed")


def test_co60_scheme_specifics():
    """Co-60 has exactly two lines (1173, 1332), each cascading with
    the other at p≈1.0."""
    co60 = CASCADE_SCHEMES["Co-60"]
    assert set(co60.cascades.keys()) == {1173.23, 1332.49}
    p_1173 = co60.cascades[1173.23]
    assert len(p_1173) == 1
    assert abs(p_1173[0].E_partner_keV - 1332.49) < 0.1
    assert p_1173[0].probability > 0.99
    p_1332 = co60.cascades[1332.49]
    assert len(p_1332) == 1
    assert abs(p_1332[0].E_partner_keV - 1173.23) < 0.1
    assert p_1332[0].probability > 0.99
    print(f"  ✓ test_co60_scheme_specifics (mutual cascade, p>0.99)")


def test_na22_models_two_annihilation_photons():
    """Na-22 1274 cascade with 2× 511 keV (positron annihilation pair),
    modelled as two CascadePartner entries at 511 keV each.
    """
    na22 = CASCADE_SCHEMES["Na-22"]
    partners_1274 = na22.cascades[1274.54]
    # At least two entries pointing to ~511 keV
    e511_entries = [p for p in partners_1274
                    if abs(p.E_partner_keV - 511.0) < 1.0]
    assert len(e511_entries) >= 2, \
        f"Expected ≥2 511-keV partners for Na-22 1274, got {len(e511_entries)}"
    print(f"  ✓ test_na22_models_two_annihilation_photons "
          f"({len(e511_entries)} 511-keV partners)")


# ═════════════════════════════════════════════════════════════════════
# Group 2 — peak-to-total models
# ═════════════════════════════════════════════════════════════════════

def test_peak_to_total_NaI_anchor_values():
    """P(E) for NaI matches Gilmore Table 8.4 within ±20% at anchor energies.

    Gilmore reference (3×3" NaI, point geometry):
        100 keV: 0.92    1000 keV: 0.35
        200 keV: 0.85    1500 keV: 0.27
        500 keV: 0.55    2000 keV: 0.22
                          2600 keV: 0.17
    The fit is degree-2 in log-log space; ±20% tolerance allows for
    the smoothing across 7 anchor points.
    """
    expected = {
        200:  0.85,
        500:  0.55,
        1000: 0.35,
        1500: 0.27,
        2000: 0.22,
        2600: 0.17,
    }
    max_err_pct = 0.0
    for E, P_ref in expected.items():
        P_model = peak_to_total_NaI(E)
        err_pct = abs(P_model - P_ref) / P_ref * 100
        max_err_pct = max(max_err_pct, err_pct)
        assert err_pct < 20.0, \
            f"P/T NaI at {E} keV: model {P_model:.3f} vs ref {P_ref:.3f} "\
            f"deviation {err_pct:.1f}% > 20%"
    print(f"  ✓ test_peak_to_total_NaI_anchor_values "
          f"(max err {max_err_pct:.1f}%)")


def test_peak_to_total_monotonic_decrease():
    """P(E) for NaI must decrease monotonically (or stay flat) with E."""
    E_arr = [100, 200, 500, 1000, 1500, 2000, 2600]
    P_prev = peak_to_total_NaI(E_arr[0])
    for E in E_arr[1:]:
        P = peak_to_total_NaI(E)
        # Allow up to 1% increase due to round-off; reject larger
        assert P <= P_prev * 1.01, \
            f"P/T increased at E={E}: {P_prev:.4f} → {P:.4f}"
        P_prev = P
    print(f"  ✓ test_peak_to_total_monotonic_decrease")


def test_peak_to_total_clipped_to_unit_range():
    """P(E) always returns values in [0.05, 1.0] (low/high-E safety)."""
    for E in (10, 50, 100, 500, 1000, 5000, 10000):
        P = peak_to_total_NaI(E)
        assert 0.05 <= P <= 1.0, f"P({E}) = {P} outside [0.05, 1.0]"
    print(f"  ✓ test_peak_to_total_clipped_to_unit_range")


# ═════════════════════════════════════════════════════════════════════
# K-21 (v1.9.0) — close-geometry P/T scaling
# ═════════════════════════════════════════════════════════════════════

def test_k21_default_geometry_factor_is_one():
    """peak_to_total_NaI without `geometry_factor` returns the
    Gilmore-5cm reference value (geometry_factor=1.0)."""
    from gamma.physics.cascade_summing import peak_to_total_NaI
    for E in (200, 500, 1000, 1500, 2614):
        P_default = peak_to_total_NaI(E)
        P_explicit_one = peak_to_total_NaI(E, geometry_factor=1.0)
        assert P_default == P_explicit_one, (
            f"P_default({E}) = {P_default} != "
            f"P_explicit_one({E}) = {P_explicit_one}"
        )
    print(f"  ✓ test_k21_default_geometry_factor_is_one")


def test_k21_geometry_factor_scales_pt_linearly():
    """At a fixed E, P scales linearly with geometry_factor (subject
    to the [0.05, 1.0] clip)."""
    from gamma.physics.cascade_summing import peak_to_total_NaI
    # 1000 keV is well inside the unclipped range
    P_full = peak_to_total_NaI(1000.0, geometry_factor=1.0)
    P_half = peak_to_total_NaI(1000.0, geometry_factor=0.5)
    P_quarter = peak_to_total_NaI(1000.0, geometry_factor=0.25)
    # Allow tiny floating-point drift but must agree to 1e-9
    assert math.isclose(P_half, P_full * 0.5, rel_tol=1e-9)
    assert math.isclose(P_quarter, P_full * 0.25, rel_tol=1e-9)
    print(f"  ✓ test_k21_geometry_factor_scales_pt_linearly")


def test_k21_geometry_pt_factor_table_present():
    """The empirical GEOMETRY_PT_FACTOR registry must include all
    K-46 geometries and have non-negative factors ≤ 1.0."""
    from gamma.physics.cascade_summing import GEOMETRY_PT_FACTOR
    required = {"Точечная-5см", "Точечная-25см",
                "Маринелли", "Дента-120мл", "Петри-60мл"}
    missing = required - set(GEOMETRY_PT_FACTOR.keys())
    assert not missing, f"missing geometries: {missing}"
    for geom, factor in GEOMETRY_PT_FACTOR.items():
        assert 0.0 < factor <= 1.0, (
            f"GEOMETRY_PT_FACTOR[{geom!r}] = {factor} outside (0, 1]"
        )
    # Point-5cm must be 1.0 (Gilmore reference, no scaling)
    assert GEOMETRY_PT_FACTOR["Точечная-5см"] == 1.0
    print(f"  ✓ test_k21_geometry_pt_factor_table_present")


def test_k21_close_geometry_factors_smaller_than_one():
    """Close-geometry samples must have factor < 1 to model the
    larger cascade-coincidence loss at 0 cm distance."""
    from gamma.physics.cascade_summing import GEOMETRY_PT_FACTOR
    for geom in ("Маринелли", "Дента-120мл", "Петри-60мл"):
        assert GEOMETRY_PT_FACTOR[geom] < 1.0, (
            f"close-geometry {geom!r} must have factor < 1, "
            f"got {GEOMETRY_PT_FACTOR[geom]}"
        )
    print(f"  ✓ test_k21_close_geometry_factors_smaller_than_one")


def test_k21_peak_to_total_NaI_for_geometry_dispatcher():
    """peak_to_total_NaI_for_geometry returns a callable bound to the
    registered factor for the geometry token."""
    from gamma.physics.cascade_summing import (
        peak_to_total_NaI, peak_to_total_NaI_for_geometry,
        GEOMETRY_PT_FACTOR,
    )
    for geom, factor in GEOMETRY_PT_FACTOR.items():
        pt_func = peak_to_total_NaI_for_geometry(geom)
        for E in (500, 1000, 2614):
            P_bound = pt_func(E)
            P_manual = peak_to_total_NaI(E, geometry_factor=factor)
            assert math.isclose(P_bound, P_manual, rel_tol=1e-12), (
                f"[{geom}] at E={E}: bound={P_bound} vs manual={P_manual}"
            )
    # Unknown geometry → defaults to 1.0
    pt_unknown = peak_to_total_NaI_for_geometry("HPGe-far")
    assert math.isclose(pt_unknown(1000), peak_to_total_NaI(1000),
                        rel_tol=1e-12)
    print(f"  ✓ test_k21_peak_to_total_NaI_for_geometry_dispatcher")


def test_k21_close_geometry_increases_tcs_correction():
    """Lowering geometry_factor → smaller effective P → larger ε_T →
    larger TCS correction factor C. Verify for Tl-208 in a synthetic
    case (point-5cm vs Marinelli)."""
    from gamma.physics.cascade_summing import (
        peak_to_total_NaI_for_geometry, tcs_correction_factor,
    )
    eff = StubEfficiency()
    # Tl-208 583 keV is one of the cascade lines in CASCADE_SCHEMES.
    pt_5cm = peak_to_total_NaI_for_geometry("Точечная-5см")
    pt_marin = peak_to_total_NaI_for_geometry("Маринелли")
    c_5cm = tcs_correction_factor("Tl-208", 583.19, eff, p_t_func=pt_5cm)
    c_marin = tcs_correction_factor("Tl-208", 583.19, eff, p_t_func=pt_marin)
    assert c_marin > c_5cm, (
        f"K-21 close-geometry TCS must be larger: "
        f"5cm C={c_5cm:.4f}, Marinelli C={c_marin:.4f}"
    )
    print(f"  ✓ test_k21_close_geometry_increases_tcs_correction "
          f"(5cm C={c_5cm:.3f} vs Marinelli C={c_marin:.3f})")


def test_HPGe_PT_higher_than_NaI_at_high_E():
    """HPGe peak-to-total ratio exceeds NaI's at high E (≥ 1 MeV),
    because at high E NaI's smaller photopeak fraction dominates.
    At low E (≤500 keV) the larger NaI crystal volume can win.
    """
    for E in (1000, 1500, 2000):
        P_NaI = peak_to_total_NaI(E)
        P_HPGe = peak_to_total_HPGe(E)
        assert P_HPGe > P_NaI, \
            f"At {E} keV: HPGe P={P_HPGe:.3f} not > NaI P={P_NaI:.3f}"
    print(f"  ✓ test_HPGe_PT_higher_than_NaI_at_high_E")


# ═════════════════════════════════════════════════════════════════════
# Group 3 — total efficiency and TCS factor
# ═════════════════════════════════════════════════════════════════════

def test_total_efficiency_formula():
    """ε_T = ε_p / P (definition check)."""
    eff = StubEfficiency()
    for E in (200, 500, 1000, 1500):
        eps_p = eff.efficiency_at(E)
        P = peak_to_total_NaI(E)
        eps_T_expected = eps_p / P
        eps_T_actual = total_efficiency(E, eff)
        assert math.isclose(eps_T_actual, eps_T_expected, rel_tol=1e-9)
    print(f"  ✓ test_total_efficiency_formula")


def test_total_efficiency_zero_when_eps_p_zero():
    """If photopeak efficiency is 0, total efficiency is 0 (avoid div/0)."""
    assert total_efficiency(1000.0, ZeroEfficiency()) == 0.0
    print(f"  ✓ test_total_efficiency_zero_when_eps_p_zero")


def test_co60_tcs_correction_factor_range():
    """Co-60 TCS correction at point 5cm NaI 63×63 should be in
    [1.01, 1.10]. For our specific ε_p calibration the value is around
    1.022 (~2.2% correction).
    """
    eff = StubEfficiency()
    c_1173 = tcs_correction_factor("Co-60", 1173.23, eff)
    c_1332 = tcs_correction_factor("Co-60", 1332.49, eff)
    assert 1.01 < c_1173 < 1.10, \
        f"Co-60 1173 TCS factor {c_1173:.4f} outside reasonable range"
    assert 1.01 < c_1332 < 1.10, \
        f"Co-60 1332 TCS factor {c_1332:.4f} outside reasonable range"
    # Both lines should give similar corrections (symmetric scheme)
    assert abs(c_1173 - c_1332) < 0.005, \
        f"Asymmetric Co-60 corrections: {c_1173:.4f} vs {c_1332:.4f}"
    print(f"  ✓ test_co60_tcs_correction_factor_range "
          f"(C₁₁₇₃={c_1173:.4f}, C₁₃₃₂={c_1332:.4f})")


def test_cs137_tcs_correction_is_one():
    """Non-cascade nuclide (Cs-137) returns C=1.0 (no correction)."""
    eff = StubEfficiency()
    c = tcs_correction_factor("Cs-137", 661.66, eff)
    assert c == 1.0
    # And the dict version returns empty
    assert compute_tcs_corrections("Cs-137", eff) == {}
    print(f"  ✓ test_cs137_tcs_correction_is_one")


def test_unknown_nuclide_returns_one():
    """Unknown nuclide returns C=1.0 silently (no exception)."""
    eff = StubEfficiency()
    c = tcs_correction_factor("Xx-999", 500.0, eff)
    assert c == 1.0
    print(f"  ✓ test_unknown_nuclide_returns_one")


def test_unknown_line_returns_one():
    """Known nuclide, unknown line: C=1.0 (line not in scheme)."""
    eff = StubEfficiency()
    # 500 keV isn't a Co-60 line
    c = tcs_correction_factor("Co-60", 500.0, eff)
    assert c == 1.0
    print(f"  ✓ test_unknown_line_returns_one")


def test_energy_tolerance_matching():
    """tcs_correction_factor finds the scheme line within ±0.5 keV by default."""
    eff = StubEfficiency()
    # 1173.0 keV (close to scheme 1173.23)
    c = tcs_correction_factor("Co-60", 1173.0, eff)
    assert c > 1.01, "1173.0 should match 1173.23 within default tol"
    # 1180 keV (≈7 keV off) should NOT match
    c2 = tcs_correction_factor("Co-60", 1180.0, eff)
    assert c2 == 1.0
    print(f"  ✓ test_energy_tolerance_matching")


def test_loss_cap_safety():
    """Even with very high (artificial) ε_T, correction factor capped
    at 1/(1-loss_cap)."""
    # Use a fake ε_curve that returns absurdly high values
    class HugeEff:
        def efficiency_at(self, E):
            return 10.0  # nonphysical but tests the cap

    c = tcs_correction_factor("Co-60", 1173.23, HugeEff(),
                                loss_cap=0.5)
    # With loss_cap=0.5, max factor is 1/(1-0.5) = 2.0
    assert c <= 2.001, f"Cap violated: {c}"
    print(f"  ✓ test_loss_cap_safety (capped at {c:.3f})")


# ═════════════════════════════════════════════════════════════════════
# Group 4 — compute_tcs_corrections dispatcher
# ═════════════════════════════════════════════════════════════════════

def test_compute_tcs_corrections_returns_dict_with_all_lines():
    """compute_tcs_corrections returns one entry per primary line in
    the nuclide's scheme."""
    eff = StubEfficiency()
    for nuc, scheme in CASCADE_SCHEMES.items():
        cc = compute_tcs_corrections(nuc, eff)
        assert set(cc.keys()) == set(scheme.line_energies()), \
            f"{nuc}: dict keys mismatch scheme lines"
        for E, C in cc.items():
            assert C >= 1.0, f"{nuc} {E}: factor {C} < 1.0"
    print(f"  ✓ test_compute_tcs_corrections_returns_dict_with_all_lines")


def test_compute_tcs_corrections_compatible_with_compute_activity():
    """The returned dict can be passed as `coincidence_correction` arg
    to compute_activity. Run a synthetic Co-60 case and verify the
    activity scales up by the correction factor.
    """
    from gamma.activity import compute_activity
    from gamma.identification.identify import (
        NuclideIdentification, LineMatch,
    )

    eff = StubEfficiency()
    # Fabricate a NuclideIdentification with one Co-60 line and a
    # known peak area
    lm = LineMatch(
        nuclide="Co-60",
        library_E_keV=1173.23, library_I_pct=99.85,
        peak_channel=400, peak_E_keV=1173.0,
        peak_sigma=50.0, residual_keV=0.23,
        is_characteristic=True,
        peak_area=1.0e5, peak_area_uncertainty=1.0e3,
    )
    ni = NuclideIdentification(
        nuclide="Co-60", detected=True,
        reason="synthetic test fixture",
        characteristic_line_keV=1173.23,
        matched_lines=(lm,),
    )
    # Without TCS
    r_no = compute_activity(
        ni, efficiency_curve=eff, live_time_s=1800.0,
        from_bg_subtracted=True,
    )
    # With TCS
    tcs = compute_tcs_corrections("Co-60", eff)
    r_yes = compute_activity(
        ni, efficiency_curve=eff, live_time_s=1800.0,
        from_bg_subtracted=True,
        coincidence_correction=tcs,
    )
    # Activity must scale up exactly by C(1173.23)
    C = tcs[1173.23]
    ratio = r_yes.A_Bq / r_no.A_Bq
    assert math.isclose(ratio, C, rel_tol=1e-6), \
        f"A scaling {ratio:.6f} ≠ TCS factor {C:.6f}"
    print(f"  ✓ test_compute_tcs_corrections_compatible_with_compute_activity "
          f"(scaling matches C={C:.4f})")


# ═════════════════════════════════════════════════════════════════════
# Group 5 — geometry dependence and order-of-magnitude reality check
# ═════════════════════════════════════════════════════════════════════

def test_tcs_scales_with_efficiency():
    """Doubling ε_p doubles ε_T → larger TCS correction."""
    eff_low  = StubEfficiency(a=12.6, b=-1.05)   # baseline ~real
    eff_high = StubEfficiency(a=25.2, b=-1.05)   # 2× efficiency
    c_low  = tcs_correction_factor("Co-60", 1173.23, eff_low)
    c_high = tcs_correction_factor("Co-60", 1173.23, eff_high)
    # Loss doubles → factor changes more than 2× the original loss
    L_low  = 1.0 - 1.0 / c_low
    L_high = 1.0 - 1.0 / c_high
    assert L_high > L_low, "High-ε case must produce larger loss"
    assert math.isclose(L_high / L_low, 2.0, rel_tol=0.10), \
        f"Loss ratio {L_high/L_low:.3f} should be ≈ 2.0 (efficiency doubled)"
    print(f"  ✓ test_tcs_scales_with_efficiency "
          f"(loss_low={L_low*100:.2f}%, loss_high={L_high*100:.2f}%)")


def test_eu152_multiple_corrections_returned():
    """Eu-152 has many cascade lines → compute_tcs_corrections returns
    a dict with multiple entries."""
    eff = StubEfficiency()
    cc = compute_tcs_corrections("Eu-152", eff)
    assert len(cc) >= 5, f"Expected ≥5 Eu-152 corrections, got {len(cc)}"
    # All factors should be > 1
    assert all(c > 1.0 for c in cc.values())
    print(f"  ✓ test_eu152_multiple_corrections_returned "
          f"({len(cc)} lines, max factor={max(cc.values()):.3f})")


# ═════════════════════════════════════════════════════════════════════
# Runner
# ═════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Running Phase 2.1g cascade-summing tests (F-31b)...\n")

    # Group 1 — cascade scheme catalogue (4)
    test_catalogue_contains_major_cascade_nuclides()
    test_cascade_scheme_data_well_formed()
    test_co60_scheme_specifics()
    test_na22_models_two_annihilation_photons()

    # Group 2 — peak-to-total models (4)
    test_peak_to_total_NaI_anchor_values()
    test_peak_to_total_monotonic_decrease()
    test_peak_to_total_clipped_to_unit_range()
    test_HPGe_PT_higher_than_NaI_at_high_E()

    # Group 3 — TCS formula (7)
    test_total_efficiency_formula()
    test_total_efficiency_zero_when_eps_p_zero()
    test_co60_tcs_correction_factor_range()
    test_cs137_tcs_correction_is_one()
    test_unknown_nuclide_returns_one()
    test_unknown_line_returns_one()
    test_energy_tolerance_matching()
    test_loss_cap_safety()

    # Group 4 — dispatcher (2)
    test_compute_tcs_corrections_returns_dict_with_all_lines()
    test_compute_tcs_corrections_compatible_with_compute_activity()

    # Group 5 — geometry + reality (2)
    test_tcs_scales_with_efficiency()
    test_eu152_multiple_corrections_returned()

    # Group 6 — K-21 close-geometry P/T scaling (v1.9.0, 6)
    test_k21_default_geometry_factor_is_one()
    test_k21_geometry_factor_scales_pt_linearly()
    test_k21_geometry_pt_factor_table_present()
    test_k21_close_geometry_factors_smaller_than_one()
    test_k21_peak_to_total_NaI_for_geometry_dispatcher()
    test_k21_close_geometry_increases_tcs_correction()

    print("\n✓ All Phase 2.1g cascade-summing tests passed.")
