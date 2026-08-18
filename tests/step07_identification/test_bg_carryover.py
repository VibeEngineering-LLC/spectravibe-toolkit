"""bg_carryover — annotation pass tests (G4 / v1.31.0).

Pure-unit tests for `gamma.identification.bg_carryover`. No spectrum I/O,
no staged pipeline — we hand-craft inputs that mimic the shapes
`json_report._build_fit_overlay` and `_build_primary_feps` produce.

Covered:
- happy path: K-40 1461 in bg → sample K-40 1461.2 gets bg_carryover
- isolated sample line (no bg counterpart for that nuclide) → no FP
- energy outside window → no match
- nuclide mismatch (same energy, different nuclide) → no match
- empty bg catalog → returns 0
- primary_feps shape (peak_E_keV, sample_sources=None)
- S_bg_counts is emitted and matches the bg-side entry
- closest-by-energy wins when two bg candidates fall inside the window
- fwhm_e_keV=None → 3 keV flat fallback
"""
from __future__ import annotations

from gamma.identification.bg_carryover import (
    build_bg_peak_catalog,
    mark_bg_carryover,
)


# ── Fakes that mimic IdentifyResult.final_detected shape ────────────────


class _M:
    def __init__(self, peak_E_keV, peak_area=10.0):
        self.peak_E_keV = peak_E_keV
        self.peak_area = peak_area


class _NI:
    def __init__(self, nuclide, lines):
        self.nuclide = nuclide
        self.matched_lines = lines


class _Staged:
    def __init__(self, detected):
        self.final_detected = detected


def _bg(nuclides):
    """nuclides: list of (name, [(E, area), ...])"""
    detected = [_NI(n, [_M(e, a) for e, a in lines]) for n, lines in nuclides]
    return _Staged(detected)


def _fwhm_naI(E):
    """Toy NaI-class FWHM model: ~7% of E."""
    return max(0.07 * E, 1.0)


# ── build_bg_peak_catalog ───────────────────────────────────────────────


def test_build_catalog_empty_when_none():
    assert build_bg_peak_catalog(None) == []


def test_build_catalog_skips_zero_area():
    staged = _bg([("K-40", [(1461.0, 0.0), (1461.0, 12.5)])])
    cat = build_bg_peak_catalog(staged)
    assert len(cat) == 1
    assert cat[0]["S_counts"] == 12.5


def test_build_catalog_emits_expected_fields():
    staged = _bg([("Tl-208", [(2614.5, 7.0)])])
    cat = build_bg_peak_catalog(staged)
    assert cat == [{
        "nuclide": "Tl-208",
        "E_keV": 2614.5,
        "S_counts": 7.0,
        "peak_id": "pb2614",
    }]


# ── mark_bg_carryover: fit_overlay.peaks shape ──────────────────────────


def _sample_peak(E, nuc, source="singlet"):
    return {"energy_keV": E, "nuclide": nuc, "source": source}


def test_happy_path_k40_carryover():
    peaks = [_sample_peak(1461.2, "K-40")]
    cat = build_bg_peak_catalog(_bg([("K-40", [(1461.0, 12.5)])]))
    n = mark_bg_carryover(peaks, cat, _fwhm_naI)
    assert n == 1
    bg = peaks[0]["bg_carryover"]
    assert bg["matched"] is True
    assert bg["bg_peak_id"] == "pb1461"
    assert bg["E_bg_keV"] == 1461.0
    assert bg["delta_E_keV"] == 0.2
    assert bg["S_bg_counts"] == 12.5


def test_isolated_sample_line_no_annotation():
    # Cs-137 662 in sample but not in bg → no false positive.
    peaks = [_sample_peak(661.6, "Cs-137")]
    cat = build_bg_peak_catalog(_bg([("K-40", [(1461.0, 12.5)])]))
    n = mark_bg_carryover(peaks, cat, _fwhm_naI)
    assert n == 0
    assert "bg_carryover" not in peaks[0]


def test_energy_outside_window_no_match():
    # 1461 vs 1450 → Δ=11 keV, > 1.5·FWHM(~102) — wait, NaI FWHM at 1461 ≈ 102.
    # Use a tighter spread: 1500 vs 1461 → Δ=39, well outside window 1.5·102=153.
    # Actually NaI FWHM(1461) = 102 keV is huge; pick a NARROWER fwhm model.
    def tight(E):
        return 2.0  # 2 keV across the board → window 3 keV
    peaks = [_sample_peak(1465.0, "K-40")]
    cat = build_bg_peak_catalog(_bg([("K-40", [(1461.0, 12.5)])]))
    n = mark_bg_carryover(peaks, cat, tight)
    assert n == 0


def test_nuclide_mismatch_no_match():
    # Same energy, different nuclide → must not match.
    peaks = [_sample_peak(1461.0, "Eu-152")]
    cat = build_bg_peak_catalog(_bg([("K-40", [(1461.0, 12.5)])]))
    n = mark_bg_carryover(peaks, cat, _fwhm_naI)
    assert n == 0


def test_empty_catalog_returns_zero():
    peaks = [_sample_peak(1461.0, "K-40")]
    assert mark_bg_carryover(peaks, [], _fwhm_naI) == 0
    assert "bg_carryover" not in peaks[0]


def test_skips_non_sample_sources_by_default():
    # source=background/secondary/unidentified → not annotated.
    peaks = [
        {"energy_keV": 1461.0, "nuclide": "K-40", "source": "background"},
        {"energy_keV": 1461.0, "nuclide": "K-40", "source": "secondary"},
    ]
    cat = build_bg_peak_catalog(_bg([("K-40", [(1461.0, 12.5)])]))
    n = mark_bg_carryover(peaks, cat, _fwhm_naI)
    assert n == 0


def test_multiplet_component_is_annotated():
    peaks = [_sample_peak(2614.5, "Tl-208", source="multiplet_component")]
    cat = build_bg_peak_catalog(_bg([("Tl-208", [(2614.5, 7.0)])]))
    n = mark_bg_carryover(peaks, cat, _fwhm_naI)
    assert n == 1
    assert peaks[0]["bg_carryover"]["matched"] is True


def test_closest_candidate_wins_within_window():
    # Two bg-side K-40 candidates at 1460 and 1462; sample at 1461.3 → pick 1462 (Δ=0.7).
    def tight(E):
        return 2.0  # window 3 keV — both inside
    peaks = [_sample_peak(1461.3, "K-40")]
    staged = _bg([("K-40", [(1460.0, 5.0), (1462.0, 9.0)])])
    cat = build_bg_peak_catalog(staged)
    n = mark_bg_carryover(peaks, cat, tight)
    assert n == 1
    bg = peaks[0]["bg_carryover"]
    assert bg["E_bg_keV"] == 1462.0
    assert bg["S_bg_counts"] == 9.0


def test_fwhm_none_falls_back_to_flat_3_keV():
    # No fwhm model → 3 keV flat. 1463 vs 1461 → Δ=2 < window=4.5 (3·1.5) → match.
    peaks = [_sample_peak(1463.0, "K-40")]
    cat = build_bg_peak_catalog(_bg([("K-40", [(1461.0, 12.5)])]))
    n = mark_bg_carryover(peaks, cat, None)
    assert n == 1


# ── primary_feps shape (peak_E_keV, no source field) ────────────────────


def test_primary_feps_shape_annotated():
    feps = [{"peak_E_keV": 1461.2, "nuclide": "K-40"}]  # no 'source' key
    cat = build_bg_peak_catalog(_bg([("K-40", [(1461.0, 12.5)])]))
    n = mark_bg_carryover(
        feps, cat, _fwhm_naI,
        sample_sources=None, energy_field="peak_E_keV",
    )
    assert n == 1
    assert feps[0]["bg_carryover"]["matched"] is True


def test_primary_feps_no_match_no_annotation():
    feps = [{"peak_E_keV": 661.6, "nuclide": "Cs-137"}]
    cat = build_bg_peak_catalog(_bg([("K-40", [(1461.0, 12.5)])]))
    n = mark_bg_carryover(
        feps, cat, _fwhm_naI,
        sample_sources=None, energy_field="peak_E_keV",
    )
    assert n == 0
    assert "bg_carryover" not in feps[0]


def test_bg_nuclide_with_bg_suffix_still_matched():
    # bg catalog entry typed as "K-40 (bg)" → stem strips suffix → still match "K-40".
    cat = [{"nuclide": "K-40 (bg)", "E_keV": 1461.0,
            "S_counts": 12.5, "peak_id": "pb1461"}]
    peaks = [_sample_peak(1461.1, "K-40")]
    n = mark_bg_carryover(peaks, cat, _fwhm_naI)
    assert n == 1
    assert peaks[0]["bg_carryover"]["bg_peak_id"] == "pb1461"