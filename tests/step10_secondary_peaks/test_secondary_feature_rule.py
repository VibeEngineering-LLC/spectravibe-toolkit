"""
Tests for F-40 / Rule 5 — Secondary-feature anti-misidentification
in `disambiguate_identifications`.

Uses the v1.7.16 `secondary_peaks_v2.json` catalog. The catalog
characterises observed Compton-edge / backscatter / escape positions
for 9 problem isotopes (Cs-137, K-40, Co-60, Na-22, Y-88, Th-228,
Tl-208, Pb-212, Ac-228) on Gamma-1S NaI 63×63.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.identification.identify import (
    LineMatch, NuclideIdentification, IdentificationResult,
)
from gamma.identification.window import IdentificationWindow
from gamma.identification.disambiguate import disambiguate_identifications


def _window() -> IdentificationWindow:
    """Build a minimal NaI identification window for fixture results."""
    return IdentificationWindow(
        detector_type="NaI",
        delta_E0_keV=30.0,
        scaling="sqrt_E",
    )


def _make_match(nuc, E_lib, E_obs, I_pct=10.0, sigma=50.0,
                is_char=True, peak_ch=None):
    return LineMatch(
        nuclide=nuc,
        library_E_keV=E_lib,
        library_I_pct=I_pct,
        peak_channel=peak_ch if peak_ch is not None else int(E_obs * 5),
        peak_E_keV=E_obs,
        peak_sigma=sigma,
        residual_keV=abs(E_obs - E_lib),
        is_characteristic=is_char,
    )


def _make_ni(nuc, matches, char_E=None):
    if char_E is None and matches:
        char_E = matches[0].library_E_keV
    return NuclideIdentification(
        nuclide=nuc,
        detected=True,
        reason="fixture",
        characteristic_line_keV=char_E,
        matched_lines=tuple(matches),
    )


def _make_result(detected_list):
    return IdentificationResult(
        detector_type="NaI",
        window=_window(),
        candidates_considered=len(detected_list) + 5,
        detected_nuclides=tuple(detected_list),
        rejected_nuclides=(),
        unmatched_peaks=(),
        notes="",
    )


# --- Positive cases: candidate falls in parent's secondary range ---

def test_co60_compton_edge_demotes_single_line_candidate():
    """Co-60 1173 keV Compton edge [906.85..912.50] swallows a 910 keV
    single-line candidate (Ac-228 911.20 keV — the classic case)."""
    co60 = _make_ni("Co-60", [
        _make_match("Co-60", 1173.23, 1168.5),
        _make_match("Co-60", 1332.49, 1330.0, is_char=False),
    ])
    ac228 = _make_ni("Ac-228", [
        _make_match("Ac-228", 911.20, 910.0),  # in Co-60 1173 edge range
    ])
    result = _make_result([co60, ac228])
    refined = disambiguate_identifications(result)
    names_detected = [ni.nuclide for ni in refined.detected_nuclides]
    names_rejected = [ni.nuclide for ni in refined.rejected_nuclides]
    assert "Ac-228" not in names_detected, "Ac-228 must be demoted"
    assert "Ac-228" in names_rejected
    assert "Co-60" in names_detected, "Parent Co-60 must remain"
    # Reason should mention Compton edge and Co-60
    rej = [r for r in refined.rejected_nuclides if r.nuclide == "Ac-228"][0]
    assert "Co-60" in rej.reason and "compton_edge" in rej.reason
    print("  ✓ test_co60_compton_edge_demotes_single_line_candidate")


def test_cs137_backscatter_demotes_191_kev_candidate():
    """Cs-137 backscatter [186.76..196.68] catches a 191 keV candidate."""
    cs137 = _make_ni("Cs-137", [
        _make_match("Cs-137", 661.66, 657.5),
    ])
    # Imaginary candidate with single line at 191 keV (e.g. spurious Hg-203
    # 279 misidentification — fake fixture; the point is the rule fires).
    spur = _make_ni("FakeSpur", [
        _make_match("FakeSpur", 191.0, 191.0),
    ])
    result = _make_result([cs137, spur])
    refined = disambiguate_identifications(result)
    names = [ni.nuclide for ni in refined.detected_nuclides]
    assert "FakeSpur" not in names
    rej_names = [ni.nuclide for ni in refined.rejected_nuclides]
    assert "FakeSpur" in rej_names
    print("  ✓ test_cs137_backscatter_demotes_191_kev_candidate")


def test_k40_compton_edge_demotes_candidate_at_1179():
    """K-40 Compton edge [1178.9..1179.1] is a tight cluster — catches
    a candidate at 1179 keV (e.g. nominal 1173 Co-60 misplaced by drift)."""
    k40 = _make_ni("K-40", [
        _make_match("K-40", 1460.82, 1453.0),
    ])
    spur = _make_ni("FakeFoo", [
        _make_match("FakeFoo", 1179.0, 1179.0),
    ])
    result = _make_result([k40, spur])
    refined = disambiguate_identifications(result)
    names = [ni.nuclide for ni in refined.detected_nuclides]
    assert "FakeFoo" not in names
    print("  ✓ test_k40_compton_edge_demotes_candidate_at_1179")


# --- Negative cases: rule must NOT fire ---

def test_multi_line_candidate_not_demoted_even_if_one_is_secondary():
    """Candidate with line in secondary AND a line outside any secondary —
    secondary_max_lines=2 default tolerates up to 2 lines but requires ALL
    to be in secondaries. One outside → no demotion."""
    co60 = _make_ni("Co-60", [
        _make_match("Co-60", 1173.23, 1168.5),
        _make_match("Co-60", 1332.49, 1330.0, is_char=False),
    ])
    bi214 = _make_ni("Bi-214", [
        _make_match("Bi-214", 609.31, 609.0),  # NOT in any Co-60 secondary
        _make_match("Bi-214", 1120.29, 910.0, is_char=False),  # in Co-60 edge
    ])
    result = _make_result([co60, bi214])
    refined = disambiguate_identifications(result)
    names = [ni.nuclide for ni in refined.detected_nuclides]
    assert "Bi-214" in names, ("Bi-214 must survive — has a line outside "
                              "any parent's secondary range")
    print("  ✓ test_multi_line_candidate_not_demoted_even_if_one_is_secondary")


def test_strong_multi_line_evidence_not_demoted():
    """Candidate with >secondary_max_lines matched lines is never
    demoted by Rule 5 regardless of position."""
    co60 = _make_ni("Co-60", [
        _make_match("Co-60", 1173.23, 1168.5),
        _make_match("Co-60", 1332.49, 1330.0, is_char=False),
    ])
    # Suppose all three lines coincidentally fall in Co-60 secondaries
    # (artificial — point is the rule respects the line-count threshold).
    foo = _make_ni("Bar-XX", [
        _make_match("Bar-XX", 909.0, 909.0),  # in Co-60 1173 edge
        _make_match("Bar-XX", 910.0, 910.0, is_char=False),
        _make_match("Bar-XX", 911.0, 911.0, is_char=False),  # 3 lines > default 2
    ])
    result = _make_result([co60, foo])
    refined = disambiguate_identifications(result, secondary_max_lines=2)
    names = [ni.nuclide for ni in refined.detected_nuclides]
    assert "Bar-XX" in names, ("Multi-line candidate (3 > 2 threshold) "
                                "must not be demoted by Rule 5")
    print("  ✓ test_strong_multi_line_evidence_not_demoted")


def test_no_parent_in_catalog_rule_inert():
    """If no detected nuclide has a v2-catalog entry, rule is inert."""
    # Use a nuclide NOT in the catalog as the only "parent".
    # The v2 catalog covers Cs-137, K-40, Co-60, Na-22, Y-88,
    # Th-228, Tl-208, Pb-212, Ac-228. Use Eu-152 (not in catalog).
    eu152 = _make_ni("Eu-152", [
        _make_match("Eu-152", 121.78, 121.0),
        _make_match("Eu-152", 1408.0, 1408.0, is_char=False),
    ])
    candidate = _make_ni("Bi-214", [
        _make_match("Bi-214", 609.31, 909.0),  # would be Co-60 edge,
                                                # but no Co-60 detected
    ])
    result = _make_result([eu152, candidate])
    refined = disambiguate_identifications(result)
    names = [ni.nuclide for ni in refined.detected_nuclides]
    assert "Bi-214" in names, ("No catalog parent detected — Rule 5 must "
                                "not fire")
    print("  ✓ test_no_parent_in_catalog_rule_inert")


def test_photopeak_collision_not_handled_by_rule_5():
    """A candidate at the parent's photopeak energy is NOT demoted by
    Rule 5 (that's the domain of NAI_CONFUSION_MAP / Rule 3)."""
    co60 = _make_ni("Co-60", [
        _make_match("Co-60", 1173.23, 1168.5),
        _make_match("Co-60", 1332.49, 1330.0, is_char=False),
    ])
    # Candidate whose only line is at the Co-60 1332 photopeak — pure
    # photopeak collision. Rule 5 must NOT demote (no non-photopeak hit).
    foo = _make_ni("PhotoColl", [
        _make_match("PhotoColl", 1330.0, 1330.0),
    ])
    result = _make_result([co60, foo])
    refined = disambiguate_identifications(result)
    # PhotoColl may or may not be demoted by Rule 3 (CI tiebreaker), but
    # we assert Rule 5 alone does not handle it — examine notes for
    # the F-40 marker.
    notes = refined.notes or ""
    if "PhotoColl" in notes:
        # If Rule 3 picked it up, that's fine — but the F-40 marker
        # "secondary features" should NOT mention PhotoColl.
        f40_lines = [ln for ln in notes.split("\n")
                     if "PhotoColl" in ln and "secondary" in ln.lower()]
        assert not f40_lines, ("Rule 5 must not act on pure photopeak "
                                "collisions")
    print("  ✓ test_photopeak_collision_not_handled_by_rule_5")


def test_opt_out_disables_rule():
    """`apply_secondary_feature_rule=False` disables Rule 5."""
    co60 = _make_ni("Co-60", [
        _make_match("Co-60", 1173.23, 1168.5),
        _make_match("Co-60", 1332.49, 1330.0, is_char=False),
    ])
    ac228 = _make_ni("Ac-228", [
        _make_match("Ac-228", 911.20, 910.0),  # would be demoted by default
    ])
    result = _make_result([co60, ac228])
    refined = disambiguate_identifications(
        result, apply_secondary_feature_rule=False,
    )
    names = [ni.nuclide for ni in refined.detected_nuclides]
    assert "Ac-228" in names, ("Opt-out must keep Ac-228 detected "
                                "(rule disabled)")
    print("  ✓ test_opt_out_disables_rule")


def test_parent_itself_not_demoted_by_own_secondary():
    """Parent's own photopeak detection must NEVER be demoted by Rule 5.
    Even if its own matched line is in its own secondary range
    (pathological — should never happen — but be defensive)."""
    # Co-60 photopeak detection — should remain.
    co60 = _make_ni("Co-60", [
        _make_match("Co-60", 1173.23, 1168.5),
        _make_match("Co-60", 1332.49, 1330.0, is_char=False),
    ])
    result = _make_result([co60])
    refined = disambiguate_identifications(result)
    names = [ni.nuclide for ni in refined.detected_nuclides]
    assert "Co-60" in names
    print("  ✓ test_parent_itself_not_demoted_by_own_secondary")


if __name__ == "__main__":
    print("Running F-40 secondary-feature anti-misidentification tests...\n")
    test_co60_compton_edge_demotes_single_line_candidate()
    test_cs137_backscatter_demotes_191_kev_candidate()
    test_k40_compton_edge_demotes_candidate_at_1179()
    test_multi_line_candidate_not_demoted_even_if_one_is_secondary()
    test_strong_multi_line_evidence_not_demoted()
    test_no_parent_in_catalog_rule_inert()
    test_photopeak_collision_not_handled_by_rule_5()
    test_opt_out_disables_rule()
    test_parent_itself_not_demoted_by_own_secondary()
    print("\n✓ All 9 F-40 tests passed.")
