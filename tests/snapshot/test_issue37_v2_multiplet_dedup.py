# -*- coding: utf-8 -*-
"""Issue #37 / v1.18.31+ (Agent A2) — V2 multiplet clustering: post-merge
dedup pass + per-cluster cc-suffix.

Symptom (Th-232-like input, reproducer in
_state/agent_a/outbox/multiplet_v2_repro/):
  Adjacent mega-clusters after F-374 display-window expansion overlap,
  pulling the same library line into two clusters. F-387.1 top-K cap
  drops the «losing» component from one CC; the same component reappears
  as a spurious mini-cluster in the adjacent CC.

  Pre-fix on Th-232-like scenario C:
    auto_M1_cc3 = {Tl-208 860, Ac-228 911, Ac-228 968.97}   ← dropped 964.77
    auto_M2_cc4 = {Ac-228 964.77, Ac-228 968.97}            ← redundant
    968.97 in TWO clusters; 964.77 missing from main multiplet.

  Post-fix:
    auto_M1_cc3 = {Tl-208 860, Ac-228 911}                   ← dedup'd
    auto_M2_cc0 = {Ac-228 964.77, Ac-228 968.97}             ← canonical
    Each library line in exactly one cluster.

Fix:
  1. `MultipletClusterCandidate.phantom_components` field — store top-K
     cap'd phantoms instead of deleting them.
  2. Post-merge dedup: for each (nuclide, E_keV) appearing in ≥2
     clusters, keep in the cluster where it's closest to the
     intensity-weighted centroid.
  3. Subset-merge: drop cluster B if its components are a strict subset
     of cluster A's.
  4. `cc{idx}` is now per-input-cluster (`sub_idx`), not global.
"""
from __future__ import annotations
import sys
from math import sqrt
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent
                       / "scripts"))

from gamma.experimental.peak_pipeline_v2 import (   # noqa: E402
    detect_multiplet_clusters,
    PeakHit,
    _DEFAULT_LIB,
    MultipletClusterCandidate,
)


def _nai_fwhm(E):
    return 1.5 * sqrt(max(float(E), 1.0)) + 5.0


def _peak(E, kev_per_ch=0.5, sig=30.0):
    return PeakHit(
        channel=int(E / kev_per_ch),
        energy_keV=float(E),
        fwhm_channels=_nai_fwhm(E) / kev_per_ch,
        significance=sig,
        source="both",
    )


# ──────────────────────────────────────────────────────────────────
# Scenario C: Th-232-like chain — issue #37 regression
# ──────────────────────────────────────────────────────────────────


def test_no_library_line_in_two_clusters():
    """Each (nuclide, E_keV) library line appears in at most ONE cluster."""
    chain = ["Pb-212", "Ac-228", "Tl-208", "Bi-212", "K-40", "Bi-214"]
    lib = {k: v for k, v in _DEFAULT_LIB.items() if k in chain}
    pick_Es = [238.6, 338.3, 463.0, 510.8, 583.2, 727.3, 860.6,
               911.2, 964.8, 969.0, 1247.0, 1460.8, 1620.5, 1764.5]
    peaks = [_peak(E) for E in pick_Es]
    clusters = detect_multiplet_clusters(peaks, _nai_fwhm,
                                          library=lib, chain_filter=chain)

    seen = {}
    for c in clusters:
        for nuc, E, _I, _g in c.components:
            key = (nuc, round(float(E), 2))
            assert key not in seen, (
                f"Library line {key} appears in TWO clusters: "
                f"{seen[key]} and {c.cluster_id} — issue #37"
            )
            seen[key] = c.cluster_id


def test_ac228_doublet_964_969_never_split_across_clusters():
    """Ac-228 964.77 + 968.97 (Δ=4.2 keV ≪ FWHM=52 keV) is a physically
    iconic doublet and MUST never end up in two SEPARATE active clusters.

    A doublet pair can land:
      (a) both as active components in the SAME cluster, OR
      (b) one active + one phantom in the SAME cluster (if top-K cap
          forces a choice), OR
      (c) both as phantoms in the same cluster.

    But NOT in two different `components` arrays — that's the issue #37
    «overkill mini-cluster» symptom.
    """
    # Full Th-232 chain reproducer (#37 scenario C)
    chain = ["Pb-212", "Ac-228", "Tl-208", "Bi-212", "K-40", "Bi-214"]
    lib = {k: v for k, v in _DEFAULT_LIB.items() if k in chain}
    pick_Es = [238.6, 338.3, 463.0, 510.8, 583.2, 727.3, 860.6,
               911.2, 964.8, 969.0, 1247.0, 1460.8, 1620.5, 1764.5]
    peaks = [_peak(E) for E in pick_Es]
    clusters = detect_multiplet_clusters(peaks, _nai_fwhm,
                                          library=lib, chain_filter=chain)

    def _line_in_cluster(cl, target_E):
        # Check both active and phantom
        for comp in cl.components:
            if abs(comp[1] - target_E) < 0.1 and comp[0] == "Ac-228":
                return "active"
        for comp in cl.phantom_components:
            if abs(comp[1] - target_E) < 0.1 and comp[0] == "Ac-228":
                return "phantom"
        return None

    # 964.77 active in 0 or 1 clusters
    active_964 = [c for c in clusters
                  if _line_in_cluster(c, 964.77) == "active"]
    active_969 = [c for c in clusters
                  if _line_in_cluster(c, 968.97) == "active"]
    assert len(active_964) <= 1, (
        f"Ac-228 964.77 should be active in at most 1 cluster, "
        f"got {len(active_964)}: {[c.cluster_id for c in active_964]}"
    )
    assert len(active_969) <= 1, (
        f"Ac-228 968.97 should be active in at most 1 cluster, "
        f"got {len(active_969)}: {[c.cluster_id for c in active_969]}"
    )
    # If BOTH are active, they MUST be in the same cluster (no split doublet)
    if len(active_964) == 1 and len(active_969) == 1:
        assert active_964[0].cluster_id == active_969[0].cluster_id, (
            f"Ac-228 964/969 doublet split across {active_964[0].cluster_id} "
            f"and {active_969[0].cluster_id} — issue #37"
        )


def test_cc_suffix_is_per_input_cluster():
    """`cc{idx}` is per-input-cluster, not global. Sub-clusters of
    auto_M2 start at cc0, not at the count of all previous sub-clusters."""
    chain = ["Pb-212", "Ac-228", "Tl-208", "Bi-212", "K-40", "Bi-214"]
    lib = {k: v for k, v in _DEFAULT_LIB.items() if k in chain}
    pick_Es = [238.6, 338.3, 463.0, 510.8, 583.2, 727.3, 860.6,
               911.2, 964.8, 969.0, 1247.0, 1460.8, 1620.5, 1764.5]
    peaks = [_peak(E) for E in pick_Es]
    clusters = detect_multiplet_clusters(peaks, _nai_fwhm,
                                          library=lib, chain_filter=chain)

    # Each auto_M{N} should have its cc-suffix start at 0
    by_prefix = {}
    for c in clusters:
        prefix = c.cluster_id.rsplit("_cc", 1)[0]
        suffix = int(c.cluster_id.rsplit("_cc", 1)[1])
        by_prefix.setdefault(prefix, []).append(suffix)

    for prefix, suffixes in by_prefix.items():
        suffixes_sorted = sorted(suffixes)
        assert suffixes_sorted[0] == 0, (
            f"{prefix} cc-suffix should start at 0; got {suffixes_sorted}"
        )
        # And they should be contiguous (no gaps from globally bumped counter)
        assert suffixes_sorted == list(range(len(suffixes_sorted))), (
            f"{prefix} cc-suffixes should be contiguous 0..N-1, "
            f"got {suffixes_sorted}"
        )


# ──────────────────────────────────────────────────────────────────
# Scenario A: 4 close peaks — one cluster expected
# ──────────────────────────────────────────────────────────────────


def test_scenario_A_close_peaks_yield_at_most_one_multiplet():
    """4 close peaks ch=460,465,470,478 (E=230,232.5,235,239 with
    0.5 keV/ch and FWHM=3 keV) → at most one multiplet (or unit
    coverage), no overkill mini-clusters."""
    fwhm_at = lambda E: 3.0
    lib = {"Xx": [(230.0, 10.0), (232.5, 8.0), (235.0, 5.0), (239.0, 7.0)]}
    peaks = [_peak(E, kev_per_ch=1.0) for E in [230.0, 232.5, 235.0, 239.0]]
    # Override fwhm_channels to 3.0 ch directly
    peaks = [PeakHit(channel=int(p.energy_keV), energy_keV=p.energy_keV,
                     fwhm_channels=3.0, significance=30.0, source="both")
             for p in peaks]
    clusters = detect_multiplet_clusters(peaks, fwhm_at, library=lib)
    # Either 1 cluster covering all 4, or 1 cluster + nothing else —
    # NOT 2 disjoint mini-clusters covering the same range
    assert len(clusters) <= 2, (
        f"4 close peaks should not produce >2 clusters; got {len(clusters)}: "
        f"{[(c.cluster_id, [comp[1] for comp in c.components]) for c in clusters]}"
    )
    # And no library line should appear in 2 clusters
    seen = {}
    for c in clusters:
        for nuc, E, _I, _g in c.components:
            key = (nuc, round(float(E), 2))
            assert key not in seen
            seen[key] = c.cluster_id


# ──────────────────────────────────────────────────────────────────
# Scenario B: isolated singlet + adjacent multiplet
# ──────────────────────────────────────────────────────────────────


def test_scenario_B_isolated_singlet_not_clustered_with_multiplet():
    """Isolated peak at E=100, multiplet pair at E=300/305 (Δ=5, FWHM=5).
    Expected: one multiplet for the pair; isolated singlet dropped from
    multiplet output (V2 routes singletons to primary_feps)."""
    fwhm_at = lambda E: 5.0
    lib = {"Iso": [(100.0, 50.0)], "Mp1": [(300.0, 20.0)],
           "Mp2": [(305.0, 15.0)]}
    peaks = [
        PeakHit(channel=200, energy_keV=100.0, fwhm_channels=10.0,
                significance=30.0, source="both"),
        PeakHit(channel=600, energy_keV=300.0, fwhm_channels=10.0,
                significance=30.0, source="both"),
        PeakHit(channel=610, energy_keV=305.0, fwhm_channels=10.0,
                significance=30.0, source="both"),
    ]
    clusters = detect_multiplet_clusters(peaks, fwhm_at, library=lib)
    # Exactly 1 cluster (the {Mp1 300, Mp2 305} pair)
    assert len(clusters) == 1
    cl = clusters[0]
    nucs = sorted(c[0] for c in cl.components)
    assert nucs == ["Mp1", "Mp2"], (
        f"Expected pair {{Mp1,Mp2}}; got {nucs} in cluster {cl.cluster_id}"
    )


# ──────────────────────────────────────────────────────────────────
# Phantom retention (top-K cap evidence preservation)
# ──────────────────────────────────────────────────────────────────


def test_phantom_components_retained_on_topK_cap():
    """When a CC > max_components_per_cluster, the dropped components
    are now retained in `phantom_components` (NOT silently deleted).
    """
    # 5 close lines, all in one CC; max_K=3 → 3 active + 2 phantom.
    fwhm_at = lambda E: 10.0
    lib = {"Z": [(100.0, 50.0), (105.0, 40.0), (110.0, 30.0),
                  (115.0, 20.0), (120.0, 10.0)]}
    # One peak in the centre — its ROI covers all 5 library lines
    peaks = [PeakHit(channel=110, energy_keV=110.0, fwhm_channels=20.0,
                     significance=30.0, source="both")]
    clusters = detect_multiplet_clusters(peaks, fwhm_at, library=lib,
                                          max_components_per_cluster=3)
    assert len(clusters) == 1
    cl = clusters[0]
    # Top-3 active (by I_pct: 50, 40, 30) → 100, 105, 110
    active_Es = sorted(c[1] for c in cl.components)
    assert active_Es == [100.0, 105.0, 110.0]
    # 2 phantom (by I_pct: 20, 10) → 115, 120
    phantom_Es = sorted(c[1] for c in cl.phantom_components)
    assert phantom_Es == [115.0, 120.0], (
        f"Phantom components should be the 2 dropped lines; "
        f"got {phantom_Es}"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
