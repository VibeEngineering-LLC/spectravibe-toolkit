"""
Tests for proportionality check, rare-isotope prior, and decay-chain
equilibrium model — Phase 1.4.5 disambiguation enhancements.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.identification.proportionality import (
    check_intensity_proportionality, get_prior,
    RARE_ISOTOPE_PRIOR,
)
from gamma.identification.chain_equilibrium import (
    check_ra226_chain_equilibrium, RA_226_CHAIN_GROUPS,
)
from gamma.identification.identify import LineMatch


def make_match(nuclide, E, I_pct, sigma):
    """Build a LineMatch for testing."""
    return LineMatch(
        nuclide=nuclide,
        library_E_keV=E,
        library_I_pct=I_pct,
        peak_channel=int(E * 5),  # fake ch
        peak_E_keV=E,
        peak_sigma=sigma,
        residual_keV=0.5,
        is_characteristic=False,
        significance_currie=sigma,
    )


# --- Proportionality tests ---

def test_proportionality_passes_when_consistent():
    """Co-60: 1173 (99.85%) + 1332 (99.98%) with σ ≈ equal → pass."""
    matches = [
        make_match("Co-60", 1173.23, 99.85, 100.0),
        make_match("Co-60", 1332.49, 99.98, 100.0),
    ]
    result = check_intensity_proportionality("Co-60", matches)
    assert result.passed, f"Co-60 equal intensities should pass: {result.reason}"
    print(f"  ✓ test_proportionality_passes_when_consistent: {result.reason}")


def test_proportionality_fails_when_inconsistent():
    """If matches have wildly wrong ratio, proportionality fails."""
    # Library: 100 vs 50 (ratio 2). Observed σ: 10 vs 100 (ratio 0.1).
    # Off by factor 20 — should fail at factor-3 tolerance.
    matches = [
        make_match("X", 500.0, 100.0, 10.0),
        make_match("X", 1000.0, 50.0, 100.0),
    ]
    result = check_intensity_proportionality("X", matches)
    assert not result.passed, "Wildly inconsistent ratios should fail"
    print(f"  ✓ test_proportionality_fails_when_inconsistent: {result.reason}")


def test_proportionality_defers_single_line():
    """With only 1 matched line, can't reject — defer judgement."""
    matches = [make_match("X", 661.66, 85.1, 50.0)]
    result = check_intensity_proportionality("X", matches)
    assert result.passed, "Single-line should defer (pass)"
    assert "defer" in result.reason.lower() or "≥2" in result.reason
    print(f"  ✓ test_proportionality_defers_single_line")


def test_proportionality_low_intensity_filter():
    """Lines below threshold are filtered out — pretend all are weak."""
    matches = [
        make_match("X", 100.0, 0.5, 10.0),  # below 1% threshold
        make_match("X", 200.0, 0.3, 10.0),  # below 1% threshold
    ]
    result = check_intensity_proportionality(
        "X", matches, min_intensity_threshold_pct=1.0,
    )
    assert result.passed
    assert "Defer" in result.reason or "defer" in result.reason
    print(f"  ✓ test_proportionality_low_intensity_filter")


def test_rare_isotope_priors():
    """Verify prior values for known rare isotopes."""
    assert get_prior("Zn-65") == 0.05, "Zn-65 should be very rare"
    assert get_prior("Cs-134") == 0.3
    assert get_prior("Co-60") == 0.5
    assert get_prior("K-40") == 1.0  # default — common
    assert get_prior("Ra-226") == 1.0  # not in rare table — common
    assert get_prior("UnknownNuclide-99") == 1.0  # default
    print(f"  ✓ test_rare_isotope_priors")


# --- Chain equilibrium tests ---

def test_chain_groups_defined():
    """Chain groups should be sensibly defined."""
    assert "Ra-226_intrinsic" in RA_226_CHAIN_GROUPS
    assert "Rn222_daughters" in RA_226_CHAIN_GROUPS
    assert "Pb210_long_lived" in RA_226_CHAIN_GROUPS
    # Group A: just Ra-226 186
    assert RA_226_CHAIN_GROUPS["Ra-226_intrinsic"][0][1] == 186.21
    # Group B should have multiple lines
    assert len(RA_226_CHAIN_GROUPS["Rn222_daughters"]) >= 3
    print(f"  ✓ test_chain_groups_defined")


def test_chain_consistent_when_equilibrium_holds():
    """When all Group B lines have proportional intensities, chain is consistent."""
    # Bi-214 lines: 609 (45.49%) and 1120 (15.31%) — ratio 2.97
    matches_by_nuc = {
        "Ra-226": [make_match("Ra-226", 186.21, 3.59, 30.0)],
        "Bi-214": [
            make_match("Bi-214", 609.31, 45.49, 100.0),
            make_match("Bi-214", 1120.29, 15.31, 34.0),  # 34/100 ≈ 15.31/45.49
        ],
        "Pb-214": [
            make_match("Pb-214", 351.93, 35.60, 80.0),
        ],
    }
    result = check_ra226_chain_equilibrium(matches_by_nuc)
    assert result.chain_consistent, f"Chain should be consistent: {result.notes}"
    # Rn retention indicator should be non-None when both A and B groups
    # have data
    assert result.rn_retention_ratio is not None
    print(f"  ✓ test_chain_consistent_when_equilibrium_holds: Rn retention "
          f"{result.rn_retention_ratio:.2f}")


def test_pb210_excluded_from_chain_validation():
    """Pb-210 detection must not affect chain consistency check."""
    # Only Pb-210 + Pb-214 — no Ra-226, no Bi-214
    matches_by_nuc = {
        "Pb-210": [make_match("Pb-210", 46.54, 4.25, 50.0)],
        "Pb-214": [
            make_match("Pb-214", 295.22, 18.41, 40.0),
            make_match("Pb-214", 351.93, 35.60, 80.0),
        ],
    }
    result = check_ra226_chain_equilibrium(matches_by_nuc)
    # Pb-214 group internally proportional (40/80 ≈ 18.41/35.60)
    rn_group = result.group_results.get("Rn222_daughters")
    assert rn_group is not None
    assert rn_group.passed, f"Pb-214 group should pass: {rn_group.reason}"
    # Pb-210 group: marked as informational
    pb_group = result.group_results.get("Pb210_long_lived")
    assert pb_group is not None
    assert "shielding" in pb_group.reason.lower() or "not used" in pb_group.reason.lower()
    print(f"  ✓ test_pb210_excluded_from_chain_validation")


def test_chain_with_rn_loss_intra_group_still_passes():
    """When Rn-222 escapes, Group A (Ra-226) and Group B (Pb-214/Bi-214)
    have different ratios — but within Group B, intensities are still
    proportional. So chain validation should still PASS within groups."""
    # Group B internally consistent: Pb-214 352 σ=80, Bi-214 609 σ=100,
    # ratio 80/100 = 0.80, library 35.6/45.49 = 0.78 — close.
    # But Group A is depleted: Ra-226 186 σ=5 (very low) — would mean
    # most Rn has escaped.
    matches_by_nuc = {
        "Ra-226": [make_match("Ra-226", 186.21, 3.59, 5.0)],  # depleted
        "Bi-214": [
            make_match("Bi-214", 609.31, 45.49, 100.0),
            make_match("Bi-214", 1120.29, 15.31, 34.0),
        ],
        "Pb-214": [
            make_match("Pb-214", 351.93, 35.60, 80.0),
            make_match("Pb-214", 295.22, 18.41, 40.0),  # 40/80 ≈ 18.41/35.6
        ],
    }
    result = check_ra226_chain_equilibrium(matches_by_nuc)
    rn_group = result.group_results["Rn222_daughters"]
    assert rn_group.passed, f"Group B should pass even with Rn loss: {rn_group.reason}"
    assert result.chain_consistent, "Chain should be CONSISTENT (each group internally proportional)"
    # Rn retention should reflect the depletion
    print(f"  ✓ test_chain_with_rn_loss_intra_group_still_passes: "
          f"Rn retention ratio = {result.rn_retention_ratio:.3f}")


if __name__ == "__main__":
    print("Running proportionality & chain equilibrium tests...\n")
    test_proportionality_passes_when_consistent()
    test_proportionality_fails_when_inconsistent()
    test_proportionality_defers_single_line()
    test_proportionality_low_intensity_filter()
    test_rare_isotope_priors()
    test_chain_groups_defined()
    test_chain_consistent_when_equilibrium_holds()
    test_pb210_excluded_from_chain_validation()
    test_chain_with_rn_loss_intra_group_still_passes()
    print("\n✓ All disambiguation tests passed.")
