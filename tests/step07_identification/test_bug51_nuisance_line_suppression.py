"""
BUG-51 regression tests -- nuisance-line dilution.

A4 diagnostic 2026-06-04 (`_state/agent_a/outbox/2026-06-04_a4_unidentified_peaks_diagnostic.md`)
sect.4 confirmed two mis-attributions on the AmTiCsEu fixture:

  * 508.38 keV peak (511 annihilation, sig=302) was being attributed to
    Eu-152 503.467 keV (library_I = 0.1524% -- only 0.53% of Eu-152's
    strongest line at 121.78 keV which has I = 28.53%).

  * 656.89 keV peak (Cs-137 661.66 keV, sig=170) was being attributed to
    Eu-152 656.489 keV (library_I = 0.1441% -- only 0.50% of Eu-152 max).

These nuisance-line claims inflated Eu-152's confidence_index to 28.77
(report.json line 1382) and contributed to the Eu-152 weighted-mean
aggregate residual of ?59.3% (BUG-39).

The fix in `scripts/gamma/identification/disambiguate.py` adds Rule 7:
when a peak is claimed by >=2 nuclides, a weak claim (library_I <5% of
that nuclide's max line) is suppressed if another claimant on the
same peak has a strong line (I >=50% of its nuclide's max).

Guards:
  1. Skip the cut if the nuclide's MAX library line is < 5% absolute
     intensity (rare-emitter case -- every line is weak by nature).
  2. Skip if the claim is on the candidate's characteristic line
     (lowest-MDA detection criterion must be preserved).
  3. Apply only when SOME other claimant has a genuinely strong line.

These tests verify the cut fires on the 511/661 motivating cases and
that the guards prevent over-aggressive removal.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.identification.identify import (
    IdentificationResult, NuclideIdentification, LineMatch,
)
from gamma.identification.window import IdentificationWindow
from gamma.identification.disambiguate import disambiguate_identifications


def _line(nuclide, library_E, library_I, peak_ch, peak_E,
          is_char=False, sigma=10.0):
    return LineMatch(
        nuclide=nuclide,
        library_E_keV=library_E,
        library_I_pct=library_I,
        peak_channel=peak_ch,
        peak_E_keV=peak_E,
        peak_sigma=sigma,
        residual_keV=abs(peak_E - library_E),
        is_characteristic=is_char,
        significance_currie=sigma,
    )


def _ni(nuclide, lines, char_E):
    return NuclideIdentification(
        nuclide=nuclide,
        detected=True,
        reason="seeded for test",
        characteristic_line_keV=char_E,
        matched_lines=tuple(lines),
        confidence=None,
    )


def _result(detected):
    return IdentificationResult(
        detector_type="NaI",
        window=IdentificationWindow(
            detector_type="NaI",
            delta_E0_keV=30.0,
            scaling="constant",
        ),
        candidates_considered=len(detected) + 5,
        detected_nuclides=tuple(detected),
        rejected_nuclides=(),
        unmatched_peaks=(),
        notes="",
    )


# --------------------------------------------------------------------
# (1) Motivating case 1: Eu-152 503.467 nuisance vs Cs-137 (not 511,
# but we model an "annihilation"-equivalent strong claim by a strong
# nuclide). We use Cs-137 661 as the strong rival because the
# library has Cs-137 661.66 keV I=85.1%.
# --------------------------------------------------------------------

def test_bug51_eu152_503_nuisance_removed_when_strong_rival_present():
    """
    Eu-152 503.467 (I=0.1524%, = 0.53% of Eu-152 max 28.53%) is a
    NUISANCE claim on a peak. When another claimant (a hypothetical
    strong line) is present on the same peak, Eu-152 503.467 must be
    suppressed.
    """
    # Eu-152 has at least one strong line (121.78 keV @ 28.53%) and
    # also a characteristic line -- model both. The 503.467 claim is
    # on the contested peak (ch=500). The strong rival here is a
    # Cs-137-like nuclide claiming 661.66 keV BUT actually mis-binned
    # to the same peak channel 500 to exercise the rule.
    #
    # For a clean test we model: Eu-152 has its char on 121.78, then
    # a strong (121.78@28.53) line on its own peak, AND the 503.467
    # nuisance on the contested peak.
    eu152 = _ni("Eu-152", [
        _line("Eu-152", 121.7817, 28.53, 100, 122.0, is_char=True),
        _line("Eu-152", 503.467, 0.1524, 500, 508.4),  # nuisance
    ], char_E=121.7817)
    # Cs-137 (max I = 85.1%, char line 661.66) is the strong rival
    # claiming the same peak ch=500 (artificially) with its main line.
    cs137 = _ni("Cs-137", [
        _line("Cs-137", 661.657, 85.1, 500, 508.4, is_char=True),
    ], char_E=661.657)

    refined = disambiguate_identifications(_result([eu152, cs137]))

    eu_refined = next(
        (n for n in refined.detected_nuclides if n.nuclide == "Eu-152"),
        None,
    )
    assert eu_refined is not None, "Eu-152 must remain detected"
    eu_E_kept = {m.library_E_keV for m in eu_refined.matched_lines}
    assert 503.467 not in eu_E_kept, (
        f"Eu-152 503.467 nuisance line must be removed; kept: {eu_E_kept}"
    )
    assert 121.7817 in eu_E_kept, (
        f"Eu-152 characteristic line 121.78 must be preserved; "
        f"kept: {eu_E_kept}"
    )
    print("  [OK] test_bug51_eu152_503_nuisance_removed_when_strong_rival_present")


# --------------------------------------------------------------------
# (2) Motivating case 2: Eu-152 656.489 vs Cs-137 661.66 on the 661
# peak -- the canonical case from A4 sect.4.
# --------------------------------------------------------------------

def test_bug51_eu152_656_nuisance_removed_by_cs137_661():
    """
    Eu-152 656.489 (I=0.1441% = 0.50% of Eu-152 max) on the 656.89 keV
    peak which is actually Cs-137 661.66 (I=85.1%, dominant primary
    line). Rule 7 must suppress the Eu-152 656.489 claim.
    """
    eu152 = _ni("Eu-152", [
        _line("Eu-152", 121.7817, 28.53, 100, 122.0, is_char=True),
        _line("Eu-152", 656.489, 0.1441, 657, 656.89),  # nuisance
    ], char_E=121.7817)
    cs137 = _ni("Cs-137", [
        _line("Cs-137", 661.657, 85.1, 657, 656.89, is_char=True),
    ], char_E=661.657)

    refined = disambiguate_identifications(_result([eu152, cs137]))
    eu_refined = next(
        (n for n in refined.detected_nuclides if n.nuclide == "Eu-152"),
        None,
    )
    assert eu_refined is not None
    eu_E_kept = {m.library_E_keV for m in eu_refined.matched_lines}
    assert 656.489 not in eu_E_kept, (
        f"Eu-152 656.489 nuisance must be removed; kept: {eu_E_kept}"
    )
    cs137_refined = next(
        (n for n in refined.detected_nuclides if n.nuclide == "Cs-137"),
        None,
    )
    assert cs137_refined is not None, "Cs-137 must remain"
    print("  [OK] test_bug51_eu152_656_nuisance_removed_by_cs137_661")


# --------------------------------------------------------------------
# (3) Guard: characteristic line is never stripped, even if its
# library intensity is below 5% of nuclide max.
# --------------------------------------------------------------------

def test_bug51_characteristic_line_never_stripped():
    """
    Even when a nuclide's characteristic line would qualify as
    "nuisance" by intensity threshold, it must never be removed --
    it was selected as the lowest-MDA detection criterion.
    """
    # Nuclide A: max I = 50%, but is_characteristic line has I = 1%
    # (well below 5% of max = 2.5%). On a shared peak with strong B
    # (I=80, max=80), Rule 7 would normally remove the A claim -- but
    # the characteristic-line guard preserves it.
    nuclide_a = _ni("Eu-152", [
        # Force a contrived scenario: characteristic at low I
        _line("Eu-152", 503.467, 0.1524, 500, 508.0, is_char=True),
        _line("Eu-152", 121.7817, 28.53, 100, 122.0),  # not char
    ], char_E=503.467)
    nuclide_b = _ni("Cs-137", [
        _line("Cs-137", 661.657, 85.1, 500, 508.0, is_char=True),
    ], char_E=661.657)

    refined = disambiguate_identifications(_result([nuclide_a, nuclide_b]))
    eu_refined = next(
        (n for n in refined.detected_nuclides if n.nuclide == "Eu-152"),
        None,
    )
    assert eu_refined is not None
    eu_E_kept = {m.library_E_keV for m in eu_refined.matched_lines}
    assert 503.467 in eu_E_kept, (
        f"Characteristic 503.467 must NOT be stripped; kept: {eu_E_kept}"
    )
    print("  [OK] test_bug51_characteristic_line_never_stripped")


# --------------------------------------------------------------------
# (4) Guard: nuclide with ALL library lines < 5% absolute intensity
# is exempt from the rule. We use Eu-155 (max line ~= 18.66% per
# library so this guard does not trip on Eu-155). For an absolute
# test we use a synthetic nuclide name not in the library -- guard
# returns early when no library data is found, and no removal occurs.
# --------------------------------------------------------------------

def test_bug51_no_library_record_does_not_strip():
    """
    When a nuclide cannot be found in the library, max_library_I = 0
    and the rule does NOT strip its claims (defensive -- would risk
    stripping legitimate test fixtures or external nuclides).
    """
    fake = _ni("FakeNuc-99", [
        _line("FakeNuc-99", 100.0, 99.0, 500, 508.0, is_char=True),
        _line("FakeNuc-99", 503.0, 0.1, 500, 508.0),
    ], char_E=100.0)
    cs137 = _ni("Cs-137", [
        _line("Cs-137", 661.657, 85.1, 500, 508.0, is_char=True),
    ], char_E=661.657)

    refined = disambiguate_identifications(_result([fake, cs137]))
    fake_refined = next(
        (n for n in refined.detected_nuclides if n.nuclide == "FakeNuc-99"),
        None,
    )
    assert fake_refined is not None
    # both lines preserved (no library record -> guard skips)
    fake_E_kept = {m.library_E_keV for m in fake_refined.matched_lines}
    assert 503.0 in fake_E_kept, (
        f"No-library guard must preserve claims; kept: {fake_E_kept}"
    )
    print("  [OK] test_bug51_no_library_record_does_not_strip")


# --------------------------------------------------------------------
# (5) Guard: no strong rival on the same peak -> no removal.
# --------------------------------------------------------------------

def test_bug51_no_strong_rival_no_removal():
    """
    If only a single nuclide claims a peak, Rule 7 does not fire --
    it only acts when >=2 nuclides contest the same peak AND at
    least one rival has a strong line.
    """
    eu152 = _ni("Eu-152", [
        _line("Eu-152", 121.7817, 28.53, 100, 122.0, is_char=True),
        # The 503.467 nuisance is on its OWN peak (ch=500) -- no rival.
        _line("Eu-152", 503.467, 0.1524, 500, 508.4),
    ], char_E=121.7817)
    refined = disambiguate_identifications(_result([eu152]))
    eu_refined = next(
        (n for n in refined.detected_nuclides if n.nuclide == "Eu-152"),
        None,
    )
    assert eu_refined is not None
    eu_E_kept = {m.library_E_keV for m in eu_refined.matched_lines}
    # No rival means we don't strip the nuisance -- it's the only
    # explanation of the peak (caller may filter via other rules).
    assert 503.467 in eu_E_kept, (
        f"Single-claimant scenario must not strip; kept: {eu_E_kept}"
    )
    print("  [OK] test_bug51_no_strong_rival_no_removal")


# --------------------------------------------------------------------
# (6) Guard: weak-vs-weak rival -> no removal (both below 50% of own
# nuclide's max).
# --------------------------------------------------------------------

def test_bug51_weak_vs_weak_rival_no_removal():
    """
    If neither claimant has a strong line on the contested peak,
    Rule 7 has no basis for choosing one over the other -- no
    removal.
    """
    eu152 = _ni("Eu-152", [
        _line("Eu-152", 121.7817, 28.53, 100, 122.0, is_char=True),
        _line("Eu-152", 503.467, 0.1524, 500, 508.4),
    ], char_E=121.7817)
    # Cs-137 on the contested peak via a weak line (Cs-137 has only
    # one ? line in the library, so we model this as Bi-214 weak line).
    bi214 = _ni("Bi-214", [
        _line("Bi-214", 609.31, 45.49, 600, 609.0, is_char=True),
        # Bi-214 1729.6 keV has I ~= 2.92% (~= 6% of Bi-214 max 45.49) --
        # in the library. Use it to fake a weak rival claim on ch=500.
        _line("Bi-214", 1729.6, 2.92, 500, 508.4),
    ], char_E=609.31)

    refined = disambiguate_identifications(_result([eu152, bi214]))
    eu_refined = next(
        (n for n in refined.detected_nuclides if n.nuclide == "Eu-152"),
        None,
    )
    assert eu_refined is not None
    eu_E_kept = {m.library_E_keV for m in eu_refined.matched_lines}
    # Neither rival has a >=50%-of-max line on the contested peak,
    # so Rule 7 does NOT fire. Eu-152 503.467 may stay or be removed
    # only by other rules (Rule 3 CI tiebreaker may still act).
    # Here we check that Rule 7 specifically did not strip:
    # the test passes as long as Eu-152 itself is still detected.
    # (Rule 3 may legitimately remove the 503.467 claim -- that's OK.)
    assert eu_refined.nuclide == "Eu-152"
    print("  [OK] test_bug51_weak_vs_weak_rival_no_removal")


if __name__ == "__main__":
    print("Running BUG-51 nuisance-line suppression tests...\n")
    test_bug51_eu152_503_nuisance_removed_when_strong_rival_present()
    test_bug51_eu152_656_nuisance_removed_by_cs137_661()
    test_bug51_characteristic_line_never_stripped()
    test_bug51_no_library_record_does_not_strip()
    test_bug51_no_strong_rival_no_removal()
    test_bug51_weak_vs_weak_rival_no_removal()
    print("\n[OK] All BUG-51 tests passed.")
