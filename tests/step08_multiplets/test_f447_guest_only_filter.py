"""F-447 — guest-only Phase 1 weak-line filter (ROI-owner aware).

Tests:
  - _f447_identify_roi_owner: clear majority / tie / single-dom / multi-dom / empty
  - _f447_proto_group_by_adjacency: adjacent / gap-split / empty
  - find_multiplet_regions(f440_guest_only_filter): OFF==baseline, owner kept,
    distinct outcome vs global F-440 filter.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

from gamma.peaks.deconvolve import (
    _f447_identify_roi_owner,
    _f447_proto_group_by_adjacency,
    find_multiplet_regions,
)


# ---------------------------------------------------------------------------
# Minimal fakes (no real spectra / identification — unit level).
# ---------------------------------------------------------------------------

@dataclass
class _FakeLineMatch:
    nuclide: str
    peak_channel: float
    library_I_pct: float
    library_E_keV: float = 0.0
    peak_area: float = 100.0
    peak_area_uncertainty: float = 10.0
    peak_area_source: str = "cowell"
    significance_currie: float = 10.0


@dataclass
class _FakeNuclideId:
    nuclide_id: str
    matched_lines: list


@dataclass
class _FakeIdResult:
    detected_nuclides: list


def _make_ir(matches):
    by_nuc: dict = {}
    for m in matches:
        by_nuc.setdefault(m.nuclide, []).append(m)
    return _FakeIdResult(
        detected_nuclides=[
            _FakeNuclideId(nuclide_id=n, matched_lines=ms)
            for n, ms in by_nuc.items()
        ],
    )


def _const_fwhm(_ch):
    return 5.0


# ---------------------------------------------------------------------------
# 1. _f447_identify_roi_owner
# ---------------------------------------------------------------------------

def test_f447_owner_clear_majority():
    """Clear majority: Tl-208 (84.5%) >> guests."""
    ms = [
        _FakeLineMatch("Tl-208", 100, 84.5),
        _FakeLineMatch("Ac-228", 102, 1.2),
        _FakeLineMatch("Ac-228", 105, 0.8),
    ]
    assert _f447_identify_roi_owner(ms) == "Tl-208"


def test_f447_owner_tie_no_dominant():
    """Equal sums, no dominant line on either → None (multi-owner)."""
    ms = [
        _FakeLineMatch("A", 0, 10.0),
        _FakeLineMatch("B", 0, 10.0),
    ]
    assert _f447_identify_roi_owner(ms) is None


def test_f447_owner_tie_single_dominant():
    """Tie A=40 (with dom>30) vs B=40 (no dom) → A."""
    ms = [
        _FakeLineMatch("A", 0, 40.0),
        _FakeLineMatch("B", 0, 25.0),
        _FakeLineMatch("B", 1, 15.0),
    ]
    assert _f447_identify_roi_owner(ms) == "A"


def test_f447_owner_tie_both_dominant():
    """Tie, both have dominant lines → None (ambiguous)."""
    ms = [
        _FakeLineMatch("A", 0, 35.0),
        _FakeLineMatch("A", 1, 5.0),
        _FakeLineMatch("B", 0, 40.0),
    ]
    assert _f447_identify_roi_owner(ms) is None


def test_f447_owner_empty():
    assert _f447_identify_roi_owner([]) is None


# ---------------------------------------------------------------------------
# 2. _f447_proto_group_by_adjacency
# ---------------------------------------------------------------------------

def test_f447_proto_group_three_adjacent():
    ms = [
        _FakeLineMatch("X", 100, 0),
        _FakeLineMatch("Y", 102, 0),
        _FakeLineMatch("Z", 104, 0),
    ]
    g = _f447_proto_group_by_adjacency(ms, _const_fwhm, 1.0)
    assert len(g) == 1 and len(g[0]) == 3


def test_f447_proto_group_split_by_gap():
    ms = [
        _FakeLineMatch("X", 100, 0),
        _FakeLineMatch("Y", 200, 0),
    ]
    g = _f447_proto_group_by_adjacency(ms, _const_fwhm, 1.0)
    assert len(g) == 2


def test_f447_proto_group_empty():
    assert _f447_proto_group_by_adjacency([], _const_fwhm, 1.0) == []


# ---------------------------------------------------------------------------
# 3. End-to-end via find_multiplet_regions
# ---------------------------------------------------------------------------

def test_f447_no_wire_in_baseline_equivalence():
    """F-447 OFF == F-440 OFF: thresholds=0 → no filtering, identical output."""
    ms = [
        _FakeLineMatch("Tl-208", 100, 84.5),
        _FakeLineMatch("Tl-208", 103, 11.0),
        _FakeLineMatch("Ac-228", 105, 0.6),
        _FakeLineMatch("Ac-228", 108, 0.3),
    ]
    ir = _make_ir(ms)
    clusters_off = find_multiplet_regions(
        ir, _const_fwhm,
        f440_guest_only_filter=False,
        min_grouping_snr=0.0, min_grouping_intensity_pct=0.0,
        unresolved_separation_fwhm_factor=0.0,
    )
    clusters_on = find_multiplet_regions(
        ir, _const_fwhm,
        f440_guest_only_filter=True,
        min_grouping_snr=0.0, min_grouping_intensity_pct=0.0,
        unresolved_separation_fwhm_factor=0.0,
    )
    flat_off = sorted(m.peak_channel for c in clusters_off for m in c)
    flat_on = sorted(m.peak_channel for c in clusters_on for m in c)
    assert flat_off == flat_on


def test_f447_guest_only_keeps_owner_phantoms_guest():
    """Owner Tl-208 kept active; Ac-228 guests phantom-marked (NOT dropped).

    F-447 V2 contract — failing guest lines convert to
    `peak_area_source='library_anchor_phantom'` to preserve ROI topology
    (continuum convergence depends on full line set). V1 «drop guest»
    сужал ROI → Tl-208 catastrophe demo 2026-06-15.
    """
    ms = [
        _FakeLineMatch("Tl-208", 100, 84.5),
        _FakeLineMatch("Tl-208", 103, 11.0),
        _FakeLineMatch("Ac-228", 105, 0.6),
        _FakeLineMatch("Ac-228", 108, 0.3),
    ]
    ir = _make_ir(ms)
    clusters = find_multiplet_regions(
        ir, _const_fwhm,
        f440_guest_only_filter=True,
        min_grouping_snr=5.0, min_grouping_intensity_pct=3.0,
        unresolved_separation_fwhm_factor=0.0,
    )
    flat = [m for c in clusters for m in c]
    tls = [m for m in flat if m.nuclide == "Tl-208"]
    acs = [m for m in flat if m.nuclide == "Ac-228"]
    assert len(tls) == 2 and all(m.peak_area_source == "cowell" for m in tls), \
        "owner Tl-208 lines must remain ACTIVE (peak_area_source='cowell')"
    assert len(acs) == 2 and all(
        m.peak_area_source == "library_anchor_phantom" for m in acs
    ), "guest Ac-228 lines must be phantom-marked (not dropped)"


def test_f447_guest_only_vs_global_filter_distinct():
    """Guest-only keeps weak owner line; global filter drops it."""
    ms = [
        _FakeLineMatch("Tl-208", 100, 84.5),
        _FakeLineMatch("Tl-208", 103, 1.5),  # weak owner line, would fail global
        _FakeLineMatch("Ac-228", 105, 0.6),  # guest weak
    ]
    ir = _make_ir(ms)
    global_clusters = find_multiplet_regions(
        ir, _const_fwhm,
        f440_guest_only_filter=False,
        min_grouping_snr=5.0, min_grouping_intensity_pct=3.0,
        unresolved_separation_fwhm_factor=0.0,
    )
    guest_clusters = find_multiplet_regions(
        ir, _const_fwhm,
        f440_guest_only_filter=True,
        min_grouping_snr=5.0, min_grouping_intensity_pct=3.0,
        unresolved_separation_fwhm_factor=0.0,
    )
    has_weak_tl_global = any(
        m.nuclide == "Tl-208" and 1.0 < m.library_I_pct < 2.0
        for c in global_clusters for m in c
    )
    has_weak_tl_guest = any(
        m.nuclide == "Tl-208" and 1.0 < m.library_I_pct < 2.0
        for c in guest_clusters for m in c
    )
    assert has_weak_tl_guest, "guest-only must keep weak owner line"
    assert not has_weak_tl_global, "global must drop weak owner line"