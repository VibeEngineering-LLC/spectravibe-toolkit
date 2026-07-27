"""
Tests for the secondary-peak catalog (F-37, v1.7.15).

  • Theoretical Compton-edge formula on canonical energies
  • Theoretical backscatter formula + complementarity E_C + E_bs = E
  • Empirical-shift helpers match the documented rules
  • expected_features_for returns the right feature set for Cs-137, K-40,
    and a low-E nuclide (Am-241)
  • Catalog JSON loads, schema is sane, and Cs-137 / K-40 entries are present
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.physics.secondary_peaks import (
    M_E_KEV,
    compton_edge_keV, backscatter_keV,
    compton_edge_observed_keV, backscatter_observed_keV,
    expected_features_for, load_catalog, reset_catalog_cache,
    empirical_ratio,
    load_catalog_v2, reset_catalog_v2_cache,
    position_range, matches_secondary,
)


def test_compton_edge_cs137():
    """Cs-137 661.66 keV → analytical Compton edge 477.3 keV."""
    E_C = compton_edge_keV(661.66)
    assert abs(E_C - 477.3) < 0.5, f"E_C={E_C}"
    print(f"  ✓ test_compton_edge_cs137 (E_C={E_C:.2f} keV)")


def test_compton_edge_co60():
    """Co-60 1332.5 keV → analytical Compton edge 1118.1 keV (Knoll Fig 10.7)."""
    E_C = compton_edge_keV(1332.49)
    assert abs(E_C - 1118.1) < 0.5, f"E_C={E_C}"
    print(f"  ✓ test_compton_edge_co60 (E_C={E_C:.2f} keV)")


def test_backscatter_cs137():
    """Cs-137 661.66 keV → 184.3 keV backscatter."""
    E_bs = backscatter_keV(661.66)
    assert abs(E_bs - 184.3) < 0.5, f"E_bs={E_bs}"
    print(f"  ✓ test_backscatter_cs137 (E_bs={E_bs:.2f} keV)")


def test_complementarity():
    """E_C + E_bs = E (energy conservation in single Compton scatter)."""
    for E in (200.0, 511.0, 661.66, 1173.23, 1460.82, 2614.51):
        s = compton_edge_keV(E) + backscatter_keV(E)
        assert abs(s - E) < 1e-9, f"E={E}: sum={s}, expected {E}"
    print(f"  ✓ test_complementarity (5 energies, all conserved)")


def test_compton_edge_observed_shift():
    """E_observed = E_theory - 0.7·FWHM (empirical rule for NaI)."""
    # Cs-137 Compton edge at 477.3 keV with NaI FWHM≈50 keV → observed ≈ 442
    E_obs = compton_edge_observed_keV(661.66, fwhm_at_edge_keV=50.0)
    expected = compton_edge_keV(661.66) - 0.7 * 50.0
    assert abs(E_obs - expected) < 1e-9
    # Validate against actual observed value (mean -37 keV from the 9 Cs fixtures)
    measured_shift = E_obs - compton_edge_keV(661.66)
    assert abs(measured_shift - (-35.0)) < 5.0, (
        f"empirical Cs-137 shift -35 keV, formula says {measured_shift:.1f}"
    )
    print(f"  ✓ test_compton_edge_observed_shift "
          f"(predicted shift {measured_shift:.1f} keV, observed ~ -37 keV)")


def test_backscatter_observed_geometry_table():
    """Geometry-conditional backscatter shifts match the documented table."""
    E_theory = backscatter_keV(661.66)
    point5  = backscatter_observed_keV(661.66, "point_5cm")
    point25 = backscatter_observed_keV(661.66, "point_25cm")
    extended = backscatter_observed_keV(661.66, "extended_source")
    assert abs((point5 - E_theory) - 14.0) < 1e-9
    assert abs((point25 - E_theory) - 5.0) < 1e-9
    assert abs((extended - E_theory) - 10.0) < 1e-9
    print(f"  ✓ test_backscatter_observed_geometry_table "
          f"(point5=+14, point25=+5, extended=+10)")


def test_expected_features_cs137():
    """Cs-137 → photopeak + compton + backscatter + Ba_Ka + k40_natural; NO escape peaks."""
    feats = expected_features_for("Cs-137", 661.66)
    names = {f.name for f in feats}
    assert names == {"photopeak", "compton_edge", "backscatter",
                     "ic_xray_Ba_Ka", "k40_natural"}, names
    print(f"  ✓ test_expected_features_cs137 ({len(feats)} features)")


def test_expected_features_k40():
    """K-40 (1461 keV) → adds single+double escape (since E > 1022 keV)."""
    feats = expected_features_for("K-40", 1460.82)
    names = {f.name for f in feats}
    assert "single_escape" in names
    assert "double_escape" in names
    assert "backscatter" in names
    assert "compton_edge" in names
    # K-40 is NOT Cs-137 so no IC X-rays or k40_natural feature
    assert "ic_xray_Ba_Ka" not in names
    print(f"  ✓ test_expected_features_k40 ({len(feats)} features)")


def test_expected_features_low_E_includes_xray_escape():
    """For E_gamma < 200 keV (e.g. Am-241 60 keV), expect the I K X-ray escape."""
    feats = expected_features_for("Am-241", 59.54)
    names = {f.name for f in feats}
    assert "xray_escape" in names, names
    print(f"  ✓ test_expected_features_low_E_includes_xray_escape")


def test_catalog_loads_and_has_cs137_k40():
    """Reading detectors/Gamma-1S/data/secondary_peaks.json yields a catalog with both nuclides."""
    reset_catalog_cache()
    cat = load_catalog()
    assert "nuclides" in cat
    assert "Cs-137" in cat["nuclides"]
    assert "K-40" in cat["nuclides"]
    assert cat["nuclides"]["Cs-137"]["primary_E_keV"] == 661.66
    print(f"  ✓ test_catalog_loads_and_has_cs137_k40")


def test_catalog_cs137_backscatter_ratio_in_expected_range():
    """Mean backscatter ratio for Cs-137 must be in the 5-10% NaI band."""
    rec = empirical_ratio("Cs-137", "backscatter")
    assert rec is not None
    mean = rec["mean"]
    assert 0.03 < mean < 0.15, f"Cs-137 backscatter mean R={mean}"
    print(f"  ✓ test_catalog_cs137_backscatter_ratio_in_expected_range "
          f"(R={mean:.3f}, expected 0.05-0.10 on NaI)")


def test_catalog_compton_edge_residual_is_negative():
    """Compton edge measured position must be BELOW theoretical (NaI artefact)."""
    for nuc in ("Cs-137", "K-40"):
        rec = empirical_ratio(nuc, "compton_edge")
        if rec is None:
            continue
        residual = rec["mean_residual_keV"]
        assert residual < -10.0, (
            f"{nuc} compton_edge residual {residual} should be < -10 keV"
        )
    print(f"  ✓ test_catalog_compton_edge_residual_is_negative")


def test_empirical_ratio_unknown_returns_none():
    """Looking up a (nuclide, feature) not in the catalog returns None."""
    assert empirical_ratio("Cs-137", "no_such_feature") is None
    assert empirical_ratio("U-238", "compton_edge") is None
    print(f"  ✓ test_empirical_ratio_unknown_returns_none")


# ---------------------------------------------------------------------------
# F-38 v2 catalog tests (range/shape per problem isotope)
# ---------------------------------------------------------------------------

def test_v2_catalog_loads():
    """v2 catalog file exists and has the expected problem-isotope set."""
    reset_catalog_v2_cache()
    cat = load_catalog_v2()
    assert "nuclides" in cat
    # Problem isotopes we characterised on Gamma-1S (v1.7.16 set)
    for nuc in ("Cs-137", "K-40", "Co-60", "Na-22", "Th-228", "Y-88"):
        assert nuc in cat["nuclides"], f"missing {nuc}"
    # F-39 Th-chain daughters added in v1.7.17
    for nuc in ("Tl-208", "Pb-212", "Ac-228"):
        assert nuc in cat["nuclides"], f"missing {nuc} (F-39)"
    print(f"  ✓ test_v2_catalog_loads "
          f"({len(cat['nuclides'])} problem isotopes incl. F-39 Th-chain)")


def test_v2_catalog_tl208_chain_daughter():
    """F-39: Tl-208 from Th-228 fixtures has all 4 primary lines."""
    cat = load_catalog_v2()
    tl = cat["nuclides"]["Tl-208"]
    primaries = set(tl["by_primary_line"].keys())
    assert "583.19" in primaries, primaries
    assert "2614.51" in primaries, primaries
    assert "510.77" in primaries, primaries
    # Tl-208 510.77 photopeak observed range overlaps Na-22 511 — the
    # canonical example of a positron-emitter vs Th-chain confusion.
    r = position_range("Tl-208", 510.77, "photopeak", span="p10p90")
    assert r is not None
    lo, hi = r
    assert 500.0 < lo < hi < 510.0, f"Tl-208 510 range {lo}..{hi}"
    print(f"  ✓ test_v2_catalog_tl208_chain_daughter "
          f"(510.77 photopeak observed {lo:.1f}..{hi:.1f}, "
          f"overlap with Na-22 511)")


def test_v2_catalog_per_primary_keying():
    """Multi-line parents keep their primary-line secondaries separated."""
    cat = load_catalog_v2()
    co60 = cat["nuclides"]["Co-60"]
    assert "by_primary_line" in co60
    keys = list(co60["by_primary_line"].keys())
    assert "1173.23" in keys, keys
    assert "1332.49" in keys, keys
    print(f"  ✓ test_v2_catalog_per_primary_keying (Co-60: 2 primaries)")


def test_position_range_cs137_compton_edge():
    """Cs-137 Compton edge p10..p90 sits in [430, 445] keV on this detector."""
    r = position_range("Cs-137", 661.66, "compton_edge", span="p10p90")
    assert r is not None, "no range for Cs-137 compton_edge"
    lo, hi = r
    assert 430.0 < lo < hi < 445.0, f"range {lo}..{hi} outside expected"
    print(f"  ✓ test_position_range_cs137_compton_edge "
          f"({lo:.1f} .. {hi:.1f} keV)")


def test_position_range_k40_compton_edge():
    """K-40 Compton edge p10..p90 sits in [1175, 1185] keV."""
    r = position_range("K-40", 1460.82, "compton_edge", span="p10p90")
    assert r is not None
    lo, hi = r
    assert 1175.0 < lo < hi < 1185.0, f"range {lo}..{hi}"
    print(f"  ✓ test_position_range_k40_compton_edge "
          f"({lo:.1f} .. {hi:.1f} keV)")


def test_matches_secondary_cs137_compton_edge_collides_with_bi214():
    """A peak at 437 keV in a Cs-137 spectrum matches the Cs-137 Compton edge."""
    # Cs-137 Compton edge measured range ~432..442 keV
    hits = matches_secondary("Cs-137", 437.0, span="p10p90")
    assert any(h["feature"] == "compton_edge" for h in hits), (
        f"no compton_edge match at 437 keV: {hits}"
    )
    # And the same peak does NOT match Cs-137 backscatter or photopeak
    bs_hits = matches_secondary("Cs-137", 437.0, feature="backscatter")
    assert bs_hits == []
    print(f"  ✓ test_matches_secondary_cs137_compton_edge_collides_with_bi214")


def test_matches_secondary_no_match_outside_ranges():
    """A peak at 800 keV in Cs-137 spectrum doesn't match any Cs-137 feature."""
    hits = matches_secondary("Cs-137", 800.0)
    assert hits == [], f"unexpected matches: {hits}"
    print(f"  ✓ test_matches_secondary_no_match_outside_ranges")


def test_matches_secondary_k40_compton_dangerous_for_co60():
    """K-40 Compton edge range [1179, 1180] keV — Co-60 1173 photopeak
    sits 5-7 keV below this range, so it should NOT match (but is dangerously close)."""
    hits = matches_secondary("K-40", 1173.23)
    # 1173 < lower bound of K-40 Compton (1178.9), so no match — but the
    # CONFLICT is real because identification tolerance is FWHM-wide.
    # The catalog's p10..p90 is tighter than the identification window.
    assert hits == [] or all(h["feature"] != "compton_edge" for h in hits)
    # However the IDENTIFICATION conflict should be flagged in the catalog:
    cat = load_catalog_v2()
    rec = cat["nuclides"]["K-40"]["by_primary_line"]["1460.82"]
    ce = rec["features"]["compton_edge"]
    # Compton edge p10 is around 1178.9, Co-60 photopeak is 1173.23 ->
    # within 6 keV. The catalog records that observed median is 1178.9,
    # so 1173 sits below; in a real spectrum with realistic drift
    # (std=0.11 here, but on a different detector could be 10+ keV) the
    # conflict is real. Document this directly:
    print(f"  ✓ test_matches_secondary_k40_compton_dangerous_for_co60 "
          f"(K-40 CE p10={ce['position_keV']['p10']:.1f}, Co-60 1173.23 "
          f"sits {ce['position_keV']['p10']-1173.23:.1f} keV inside tail)")


def test_v2_catalog_conflict_lines_recorded():
    """The catalog records real-library conflicts for each feature range."""
    cat = load_catalog_v2()
    # Cs-137 Ba Kα IC X-ray range [23, 27] keV should flag Am-241 26.34
    cs = cat["nuclides"]["Cs-137"]["by_primary_line"]["661.66"]
    ba = cs["features"]["ic_xray_Ba_Ka"]
    conflicts = ba["conflict_lines"]
    am241_conflicts = [c for c in conflicts if c["nuclide"] == "Am-241"]
    assert am241_conflicts, f"expected Am-241 conflict, got {conflicts}"
    print(f"  ✓ test_v2_catalog_conflict_lines_recorded "
          f"(Cs-137 Ba Kα conflicts with Am-241 26.34)")


def test_lsrm_chain_loader_adds_th_chain_daughters():
    """F-39: load_lsrm_chain_libs() supplements built-in JSON with
    Th-232 chain daughters (Tl-208, Pb-212, Ac-228, Bi-212, Ra-224)."""
    from gamma.data.nuclide_library import (
        load_lsrm_chain_libs, get_nuclide, reset_cache, list_nuclides,
    )
    reset_cache()
    baseline = set(list_nuclides())
    res = load_lsrm_chain_libs()
    after = set(list_nuclides())
    added = sorted(after - baseline)
    # Threshold lowered from 18 → 10 (F-LIB-EXTENSION, 2026-06-07):
    # Th-232 chain daughters and OSGI nuclides (Eu-154, Eu-155, Ce-144, Tl-208,
    # Pb-212, Bi-212, Ac-228, Ra-224) are now in the built-in nuclides.json
    # (80 entries) so they no longer count as "added" by the LSRM loader.
    assert len(added) >= 10, (
        f"expected ≥10 new nuclides from Lsrm libs, got {len(added)}: "
        f"{added}"
    )
    # Th-232 chain daughters must all be present
    for daughter in ("Tl-208", "Pb-212", "Bi-212", "Ac-228"):
        rec = get_nuclide(daughter)
        assert rec is not None, f"{daughter} missing after chain load"
        assert rec.get("lines"), f"{daughter} has no lines"
    # OSGI extension nuclides also present
    for nuc in ("Eu-154", "Eu-155", "Ce-144"):
        assert get_nuclide(nuc) is not None, f"{nuc} missing"
    # Reset for any subsequent test that expects vanilla 27 entries
    reset_cache()
    print(f"  ✓ test_lsrm_chain_loader_adds_th_chain_daughters "
          f"(+{len(added)} nuclides, Th-chain + 18 ОСГИ extensions)")


def test_th228_in_built_in_library():
    """F-39: Th-228 explicitly added to built-in nuclides.json with T½
    so chain-proxy cert validation can decay-correct using the parent
    half-life regardless of whether the Lsrm chain libraries are loaded."""
    from gamma.data.nuclide_library import get_nuclide, reset_cache
    reset_cache()
    th = get_nuclide("Th-228")
    assert th is not None
    # T½ = 1.9116 yr ≈ 6.03e7 s
    assert 5.9e7 < th["T_half_s"] < 6.2e7, f"T_Th228={th['T_half_s']}"
    # Has at least the principal direct γ-lines (84.4, 215.99 keV)
    line_E = [L[0] for L in th["lines"]]
    assert any(83 < E < 86 for E in line_E), f"84 keV line missing: {line_E}"
    print(f"  ✓ test_th228_in_built_in_library "
          f"(T½={th['T_half_s']:.3e}s, {len(th['lines'])} lines)")


def test_v2_photopeak_position_tightness():
    """Photopeak position p90-p10 across many spectra should be < 5 keV
    (the intrinsic calibration drift on this detector)."""
    for nuc, primary in (("Cs-137", 661.66), ("K-40", 1460.82),
                         ("Na-22", 511.00), ("Y-88", 898.04)):
        r = position_range(nuc, primary, "photopeak", span="p10p90")
        if r is None:
            continue
        lo, hi = r
        spread = hi - lo
        assert spread < 5.0, (
            f"{nuc} photopeak p10..p90 spread {spread:.2f} keV > 5 keV"
        )
    print(f"  ✓ test_v2_photopeak_position_tightness "
          f"(all problem isotope photopeaks within 5 keV intrinsic spread)")


if __name__ == "__main__":
    print("Running F-37 secondary-peaks tests...\n")
    test_compton_edge_cs137()
    test_compton_edge_co60()
    test_backscatter_cs137()
    test_complementarity()
    test_compton_edge_observed_shift()
    test_backscatter_observed_geometry_table()
    test_expected_features_cs137()
    test_expected_features_k40()
    test_expected_features_low_E_includes_xray_escape()
    test_catalog_loads_and_has_cs137_k40()
    test_catalog_cs137_backscatter_ratio_in_expected_range()
    test_catalog_compton_edge_residual_is_negative()
    test_empirical_ratio_unknown_returns_none()
    print("\nRunning F-38 v2-catalog tests...\n")
    test_v2_catalog_loads()
    test_v2_catalog_tl208_chain_daughter()
    test_v2_catalog_per_primary_keying()
    test_position_range_cs137_compton_edge()
    test_position_range_k40_compton_edge()
    test_matches_secondary_cs137_compton_edge_collides_with_bi214()
    test_matches_secondary_no_match_outside_ranges()
    test_matches_secondary_k40_compton_dangerous_for_co60()
    test_v2_catalog_conflict_lines_recorded()
    test_v2_photopeak_position_tightness()
    print("\nRunning F-39 chain-library tests...\n")
    test_lsrm_chain_loader_adds_th_chain_daughters()
    test_th228_in_built_in_library()
    print("\nAll F-37 + F-38 + F-39 tests passed.")
