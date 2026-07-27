"""BUG-26 / v1.18.31+ (Agent B) — peak-table dedup + phantom-zero filter.

Symptom (from screenshot of Th-232 sample peak table):
  Group A — phantom A=0 rows:
    Ac-228 209.2  → A = 0 ± 1207
    Ac-228 270.2  → A = 0 ± 1628
    Tl-208 277.4  → A = 0 ± 3194 (and 277 is actually Ac-228, cross-nuclide)
    Tl-208 860.6  → A = 0 ± 227
  Group B — multiple rows on same observed peak, all showing identical
  A = 1802 ± 217 (the Ac-228 weighted-mean):
    obs 502.5 → Tl-208 510.8 (real) + Ac-228 503.8 + 509.0 + 523.1 (phantoms)
    obs 914.0 → Ac-228 911.2 (real) + Ac-228 904.2 (phantom — same channel)

Root cause: BUG-15 (compute.py) introduced within-nuclide and cross-nuclide
dedup IN the activity-computation path. The reporting layer (_build_rows)
did NOT mirror this, so all matched_lines (including dedup-skipped ones)
rendered into the table — each with the nuclide weighted-mean A.

Fix: _build_rows now applies the same dedup before classification:
  (a) within-nuclide by (nuclide, peak_channel): keep max library_I_pct
  (b) cross-nuclide by characteristic owner
  (c) phantom-zero: drop rows with peak_area_counts == 0
  per-line: if compute.py has activities for the nuclide but the line is
  NOT in lines_used → drop (it was dedup-skipped upstream)
"""
from __future__ import annotations

import sys

sys.path.insert(0, "scripts")

from gamma.reporting.interactive_html import _build_rows  # noqa: E402


# ──── Minimal stubs (parallel to test_bug19_*) ────────────────────────


class _FakeLineActivity:
    def __init__(self, E_keV: float, A_Bq: float, sigma_A_Bq: float):
        self.E_keV = E_keV
        self.A_Bq = A_Bq
        self.sigma_A_Bq = sigma_A_Bq


class _FakeActivityResult:
    def __init__(self, nuclide: str, lines_used):
        self.nuclide = nuclide
        self.lines_used = tuple(lines_used)


class _FakeAnalysisResult:
    def __init__(self, activities, sample_mass_kg=1.6):
        self.activities = activities
        self.sample_mass_kg = sample_mass_kg


# ──── Tests ────────────────────────────────────────────────────────


def test_bug26_within_nuclide_dedup_same_channel():
    """(Ac-228, peak_channel=502) appears 3× in primary_feps (lib lines
    503.8 / 509.0 / 523.1), all mapped to the same observed peak. The
    table must show ONE row — the one with highest library_I_pct."""
    report = {
        "primary_feps": [
            # 3 Ac-228 lines collapsed into one observed peak at channel 502
            {"peak_E_keV": 502.5, "peak_area_counts": 1500,
             "peak_channel": 502, "nuclide": "Ac-228",
             "library_E_keV": 503.8, "library_I_pct": 0.7,
             "is_characteristic": False},
            {"peak_E_keV": 502.5, "peak_area_counts": 1500,
             "peak_channel": 502, "nuclide": "Ac-228",
             "library_E_keV": 509.0, "library_I_pct": 0.4,
             "is_characteristic": False},
            {"peak_E_keV": 502.5, "peak_area_counts": 1500,
             "peak_channel": 502, "nuclide": "Ac-228",
             "library_E_keV": 523.1, "library_I_pct": 1.2,  # winner
             "is_characteristic": False},
        ],
        "secondary_peaks": [],
        "diagnostics": {},
        "identified_nuclides": [{
            "nuclide": "Ac-228",
            "specific_activity_Bq_per_kg": 1802.0,
            "specific_activity_sigma_Bq_per_kg": 217.0,
        }],
    }
    rows = _build_rows(report, _FakeAnalysisResult(activities=[]))
    ac_rows = [r for r in rows if r.get("iso") == "Ac-228"]
    assert len(ac_rows) == 1, (
        "BUG-26 (a) within-nuclide dedup: (Ac-228, ch=502) must collapse to "
        "ONE row (highest library_I_pct); got {} rows: {!r}"
        .format(len(ac_rows), [r["line"] for r in ac_rows])
    )
    # Winner has library_E_keV = 523.1 (highest I_pct = 1.2)
    assert "523.1" in ac_rows[0]["line"], (
        "BUG-26 (a) tie-break: highest library_I_pct wins; got line={!r}"
        .format(ac_rows[0]["line"])
    )


def test_bug26_cross_nuclide_drops_non_characteristic_owner():
    """obs 277.4 channel: Tl-208 277.4 (non-characteristic) + Ac-228 277.4
    (is_characteristic=True). Tl-208 row must be dropped."""
    report = {
        "primary_feps": [
            {"peak_E_keV": 277.4, "peak_area_counts": 800,
             "peak_channel": 277, "nuclide": "Tl-208",
             "library_E_keV": 277.4, "library_I_pct": 6.6,
             "is_characteristic": False},
            {"peak_E_keV": 277.4, "peak_area_counts": 800,
             "peak_channel": 277, "nuclide": "Ac-228",
             "library_E_keV": 277.4, "library_I_pct": 2.2,
             "is_characteristic": True},  # characteristic owner
        ],
        "secondary_peaks": [],
        "diagnostics": {},
        "identified_nuclides": [],
    }
    rows = _build_rows(report, _FakeAnalysisResult(activities=[]))
    isos_at_277 = [r["iso"] for r in rows if "277" in r["line"]]
    assert "Tl-208" not in isos_at_277, (
        "BUG-26 (b) cross-nuclide: Tl-208 277.4 must be dropped — Ac-228 is "
        "the characteristic owner on this channel; got rows {!r}"
        .format(isos_at_277)
    )
    assert "Ac-228" in isos_at_277, (
        "BUG-26 (b): Ac-228 (owner) must survive on channel 277; got {!r}"
        .format(isos_at_277)
    )


def test_bug26_phantom_zero_filter_drops_S0_rows():
    """primary_feps with peak_area_counts == 0 carry no information and
    must NOT appear in the primary table (they cause A = 0 ± σ rows)."""
    report = {
        "primary_feps": [
            {"peak_E_keV": 209.2, "peak_area_counts": 0,
             "peak_channel": 209, "nuclide": "Ac-228",
             "library_E_keV": 209.2, "library_I_pct": 3.9},
            {"peak_E_keV": 270.2, "peak_area_counts": 0,
             "peak_channel": 270, "nuclide": "Ac-228",
             "library_E_keV": 270.2, "library_I_pct": 3.5},
            # one real line for the same nuclide
            {"peak_E_keV": 911.2, "peak_area_counts": 4500,
             "peak_channel": 911, "nuclide": "Ac-228",
             "library_E_keV": 911.2, "library_I_pct": 25.8},
        ],
        "secondary_peaks": [],
        "diagnostics": {},
        "identified_nuclides": [],
    }
    rows = _build_rows(report, _FakeAnalysisResult(activities=[]))
    energies = [r["line"] for r in rows if r.get("iso") == "Ac-228"]
    assert not any("209.2" in s for s in energies), (
        "BUG-26 (c) phantom-zero: S=0 row at 209.2 keV must be dropped; "
        "got {!r}".format(energies)
    )
    assert not any("270.2" in s for s in energies), (
        "BUG-26 (c) phantom-zero: S=0 row at 270.2 keV must be dropped; "
        "got {!r}".format(energies)
    )
    assert any("911.2" in s for s in energies), (
        "BUG-26: the real S>0 line at 911.2 keV must survive; got {!r}"
        .format(energies)
    )


def test_bug26_dedup_skipped_lines_dropped_when_activity_computed():
    """If compute.py produced an ActivityResult for Ac-228 with
    lines_used = [911.2 only] (because 904.2 was within-nuclide-dedup-
    skipped), the row for 904.2 must be dropped (not silently rendered
    with the weighted-mean activity)."""
    line_acts = [_FakeLineActivity(E_keV=911.2, A_Bq=2640.0, sigma_A_Bq=145.0)]
    ar_ac228 = _FakeActivityResult(nuclide="Ac-228", lines_used=line_acts)
    ar = _FakeAnalysisResult(activities=[ar_ac228], sample_mass_kg=1.6)

    report = {
        "primary_feps": [
            {"peak_E_keV": 914.0, "peak_area_counts": 4500,
             "peak_channel": 914, "nuclide": "Ac-228",
             "library_E_keV": 911.2, "library_I_pct": 25.8},
            # phantom: different channel but same nuclide, no LineActivity entry
            {"peak_E_keV": 904.0, "peak_area_counts": 3000,
             "peak_channel": 904, "nuclide": "Ac-228",
             "library_E_keV": 904.2, "library_I_pct": 0.85},
        ],
        "secondary_peaks": [],
        "diagnostics": {},
        "identified_nuclides": [{
            "nuclide": "Ac-228",
            "specific_activity_Bq_per_kg": 1650.0,
            "specific_activity_sigma_Bq_per_kg": 90.0,
        }],
    }
    rows = _build_rows(report, ar)
    ac_rows = [r for r in rows if r.get("iso") == "Ac-228"]
    libs = [r["line"] for r in ac_rows]
    assert any("911.2" in s for s in libs), (
        "BUG-26: real Ac-228 911.2 (with LineActivity entry) must survive; "
        "got {!r}".format(libs)
    )
    assert not any("904.2" in s for s in libs), (
        "BUG-26: Ac-228 904.2 was dedup-skipped in compute.py (no entry in "
        "lines_used) and must be dropped from the row table; got {!r}"
        .format(libs)
    )


if __name__ == "__main__":
    test_bug26_within_nuclide_dedup_same_channel()
    test_bug26_cross_nuclide_drops_non_characteristic_owner()
    test_bug26_phantom_zero_filter_drops_S0_rows()
    test_bug26_dedup_skipped_lines_dropped_when_activity_computed()
    print("OK")
