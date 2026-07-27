"""F-387.1 / v1.18.26.1 — Rayleigh-CC split + top-K cap regression.

Эта таблица тестов покрывает F-387.1 семантические изменения:
  1. Rayleigh-CC граф (per-pair edges, factor·FWHM_avg).
  2. Connected-components BFS split: один input cluster → N sub-clusters.
  3. Top-K cap по library_I_pct (default K=3): остальные → phantom anchors
     с peak_area_source="library_anchor_phantom".

Empirical-data-driven кейсы построены из реальных Th-232 demo emergence
v1.18.26 (PROD M1/M2/M3 + V2 M3/M4/M5).

Synthetic кейсы покрывают boundary, transitive closure, top-K cap,
factor=0.0 (filter disabled), и Rayleigh-vs-F-387-old семантическое
изменение.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

from gamma.peaks.deconvolve import find_multiplet_regions


# ──────────────────────────────────────────────────────────────────
# Fixtures
# ──────────────────────────────────────────────────────────────────

class _FakeLineMatch:
    """Mutable stand-in for LineMatch."""

    def __init__(self, nuclide, library_E_keV, peak_channel,
                 library_I_pct=10.0, peak_area_source=""):
        self.nuclide = nuclide
        self.library_E_keV = library_E_keV
        self.library_I_pct = library_I_pct
        self.peak_channel = peak_channel
        self.peak_E_keV = library_E_keV
        self.is_characteristic = False
        self.peak_area = None
        self.peak_area_uncertainty = None
        self.peak_area_source = peak_area_source


class _FakeNuclideId:
    def __init__(self, nuclide, matched_lines):
        self.nuclide = nuclide
        self.matched_lines = matched_lines
        self.detected = True


class _FakeIdent:
    def __init__(self, detected):
        self.detected_nuclides = detected


def _all_components(clusters):
    """Flatten all clusters → list of LineMatch."""
    out = []
    for c in clusters:
        out.extend(c)
    return out


def _is_phantom(m):
    return getattr(m, "peak_area_source", "") == "library_anchor_phantom"


# ──────────────────────────────────────────────────────────────────
# Empirical: Th-232 demo M3 prod input (12 components)
# ──────────────────────────────────────────────────────────────────

class TestEmpiricalTh232M3Prod:
    """M3 prod (v1.18.26) ROI 400-700 кэВ, 12 components:
      Ac-228 409.5, 463.0
      Ac-228 503.8, 509.0, Tl-208 510.8, Ac-228 523.1
      Ac-228 562.5, 570.9, 572.1, Tl-208 583.2 (583.4 alias),
      Ac-228 583.4 — близкая пара
      Ac-228 674.8
    FWHM~24-27 кэВ.

    Expected F-387.1 split:
      CC1 {504, 509, 510, 523}  — все пары Δ < FWHM_avg=25
      CC2 {562, 571, 572, 583, 583} → top-3 active + 2 phantom
      Singletons: {409, 463, 674} (Δ > FWHM от ближайшего)
    """

    def _build(self):
        # Mock FWHM model: 25 кэВ flat
        fwhm_at = lambda ch: 25.0
        # peak_channel ≈ E_keV для удобства (bin_width=1)
        lines = [
            ("Ac-228", 409.5,  1.74),
            ("Ac-228", 463.0,  4.40),
            ("Ac-228", 503.7,  0.10),  # I_pct низкий
            ("Ac-228", 509.0,  0.13),
            ("Tl-208", 510.77, 22.6),
            ("Ac-228", 523.1,  0.10),
            ("Ac-228", 562.5,  0.85),
            ("Ac-228", 570.9,  0.18),
            ("Ac-228", 572.1,  0.13),
            ("Tl-208", 583.19, 30.6),  # strongest in 562-583
            ("Ac-228", 583.4,  0.15),
            ("Ac-228", 674.8,  0.20),
        ]
        # Группируем по nuclide для _FakeNuclideId structure
        from collections import defaultdict
        by_nuc = defaultdict(list)
        for nuc, E, I in lines:
            by_nuc[nuc].append(_FakeLineMatch(nuc, E, E, I))
        detected = tuple(
            _FakeNuclideId(nuc, ms) for nuc, ms in by_nuc.items()
        )
        return _FakeIdent(detected), fwhm_at

    def test_m3_prod_split_into_two_subclusters_plus_three_singletons(self):
        ident, fwhm_at = self._build()
        clusters = find_multiplet_regions(
            ident, fwhm_at,
            overlap_threshold_fwhm=4.0,  # широкий, чтобы overlap собрал всё
            expand_to_display_window=False,
        )
        # Singletons определим как cluster size=1.
        # CC pair/multiplet — size>=2.
        sizes = sorted(len(c) for c in clusters)
        singletons = [c for c in clusters if len(c) == 1]
        multiplets = [c for c in clusters if len(c) >= 2]

        # Ожидаемые singletons: 409, 463, 674
        singleton_Es = sorted(c[0].library_E_keV for c in singletons)
        assert any(abs(e - 409.5) < 0.1 for e in singleton_Es), (
            f"409.5 keV должна быть singleton (Δ до 463 = 53.5 > 25); "
            f"got singletons={singleton_Es}"
        )
        assert any(abs(e - 463.0) < 0.1 for e in singleton_Es), (
            f"463.0 keV должна быть singleton (Δ до 503.7 = 40.7 > 25); "
            f"got singletons={singleton_Es}"
        )
        assert any(abs(e - 674.8) < 0.1 for e in singleton_Es), (
            f"674.8 keV должна быть singleton; got singletons={singleton_Es}"
        )

        # Ожидаемые multiplet'ы:
        # CC1 вокруг 504-523 (4 компонента — все пары Δ < 25)
        # CC2 вокруг 562-583 (5 компонент → top-3 + 2 phantom)
        assert len(multiplets) == 2, (
            f"M3 prod должен дать 2 multiplet CCs; got {len(multiplets)}: "
            f"{[[round(m.library_E_keV, 1) for m in c] for c in multiplets]}"
        )

        # CC1 — components в [500, 525] кэВ
        cc1 = next(
            c for c in multiplets
            if all(500 <= m.library_E_keV <= 530 for m in c)
        )
        cc2 = next(
            c for c in multiplets
            if all(555 <= m.library_E_keV <= 590 for m in c)
        )
        # CC1: 4 components, max_K=3 → top-3 + 1 phantom
        assert len(cc1) == 4, f"CC1 has {len(cc1)} components"
        cc1_active = [m for m in cc1 if not _is_phantom(m)]
        cc1_phantom = [m for m in cc1 if _is_phantom(m)]
        assert len(cc1_active) == 3, (
            f"CC1 active should be top-3, got {len(cc1_active)}"
        )
        assert len(cc1_phantom) == 1
        # Top by intensity: 510.77 (22.6), 509 (0.13), 503.7 (0.1)
        # Wait — 503.7=0.10, 509=0.13, 523.1=0.10 — top-3 by intensity
        # are 510.77, 509.0, then tie 503.7/523.1; phantom must be either.
        # Most important: 510.77 (Tl-208 strongest) must be active.
        assert any(abs(m.library_E_keV - 510.77) < 0.1 for m in cc1_active), \
            "Tl-208 510.77 (top intensity) должен остаться active"

        # CC2: 5 components, max_K=3 → top-3 active + 2 phantom
        assert len(cc2) == 5, f"CC2 has {len(cc2)} components"
        cc2_active = [m for m in cc2 if not _is_phantom(m)]
        cc2_phantom = [m for m in cc2 if _is_phantom(m)]
        assert len(cc2_active) == 3, (
            f"CC2 active should be top-3; got {len(cc2_active)}"
        )
        assert len(cc2_phantom) == 2
        # Top intensity: 583.19 (30.6) обязан быть active
        assert any(abs(m.library_E_keV - 583.19) < 0.1 for m in cc2_active), \
            "Tl-208 583.19 (top intensity) должен остаться active"

        # Total cluster count: 3 singletons + 2 multiplets = 5
        assert len(clusters) == 5, (
            f"M3 prod → 3 singletons + 2 multiplets = 5; got {len(clusters)}"
        )


# ──────────────────────────────────────────────────────────────────
# Empirical: PROD M1 input (Tl-208 860 + Ac-228 911 + 965 + 969)
# ──────────────────────────────────────────────────────────────────

class TestEmpiricalProdM1:
    """PROD M1 input (v1.18.26): 4 components {Tl-208 860, Ac-228 911,
    965, 969}, FWHM=38 кэВ.

    Edges (Rayleigh, FWHM_avg=38):
      860-911: Δ=51 → 1.34·FWHM_avg → RESOLVED
      911-965: Δ=54 → 1.42·FWHM_avg → RESOLVED
      965-969: Δ=4 → 0.11·FWHM_avg → UNRESOLVED
      860-965: Δ=105 → resolved
      860-969: resolved
      911-969: Δ=58 → 1.53·FWHM_avg → RESOLVED

    Expected: 1 CC {965, 969} + 2 singletons {860, 911}.
    """

    def test_prod_m1_yields_cc_plus_two_singletons(self):
        fwhm_at = lambda ch: 38.0
        detected = (
            _FakeNuclideId("Ac-228", [
                _FakeLineMatch("Ac-228", 911.2, 911.2, 25.8),
                _FakeLineMatch("Ac-228", 964.77, 964.77, 4.99),
                _FakeLineMatch("Ac-228", 968.97, 968.97, 15.8),
            ]),
            _FakeNuclideId("Tl-208", [
                _FakeLineMatch("Tl-208", 860.6, 860.6, 12.5),
            ]),
        )
        clusters = find_multiplet_regions(
            _FakeIdent(detected), fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
        )
        sizes = sorted(len(c) for c in clusters)
        assert sizes == [1, 1, 2], (
            f"PROD M1: 1 CC + 2 singletons; got sizes {sizes}: "
            f"{[[m.library_E_keV for m in c] for c in clusters]}"
        )
        # Singletons: 860 и 911
        singleton_Es = sorted(c[0].library_E_keV for c in clusters if len(c) == 1)
        assert singleton_Es == [860.6, 911.2]
        # CC: {964.77, 968.97}
        cc = next(c for c in clusters if len(c) == 2)
        cc_Es = sorted(m.library_E_keV for m in cc)
        assert cc_Es == [964.77, 968.97]


# ──────────────────────────────────────────────────────────────────
# Empirical: V2 M5 input (Bi-214 1095 + 1110 + K-40 1247 + 1459)
# ──────────────────────────────────────────────────────────────────

class TestEmpiricalV2M5:
    """V2 M5 input (v1.18.26): {1095, 1110, 1247, 1459} keV, FWHM~45.

    Edges (Rayleigh, FWHM_avg=45):
      1095-1110: Δ=15 → 0.33·FWHM → UNRESOLVED → CC
      1110-1247: Δ=137 → 3.04 → RESOLVED
      1247-1459: Δ=212 → 4.71 → RESOLVED

    Architecture note: find_multiplet_regions overlap-этап отбрасывает
    line, не попадающие в overlap-кластер (size=1). 1247 и 1459
    изолированы от {1095, 1110} (Δ > overlap_threshold) → не входят
    в clusters output → их peak_area остаётся из upstream Cowell.

    F-387.1 split применяется ВНУТРИ collected clusters: для всех 4
    компонент cluster'а (overlap_threshold достаточно большой, чтобы
    собрать всех в один) разделение даст: CC1 {1095, 1110} + 2
    singletons {1247, 1459}.

    Этот тест проверяет именно тот сценарий — overlap_threshold large
    enough чтобы собрать все 4 в один cluster (моделирует пост-F-374
    enrichment expansion).
    """

    def test_v2_m5_yields_cc_plus_two_singletons(self):
        fwhm_at = lambda ch: 45.0
        detected = (
            _FakeNuclideId("Bi-214", [
                _FakeLineMatch("Bi-214", 1094.65, 1094.65, 0.34),
                _FakeLineMatch("Bi-214", 1109.95, 1109.95, 0.50),
            ]),
            _FakeNuclideId("K-40", [
                _FakeLineMatch("K-40", 1247.0, 1247.0, 0.10),
            ]),
            _FakeNuclideId("K-40b", [
                _FakeLineMatch("K-40b", 1459.0, 1459.0, 10.66),
            ]),
        )
        # overlap_threshold=10.0 → threshold=10·0.5·90=450 кэВ.
        # Δ_max между 1094 и 1459 = 365 < 450 → все 4 в один cluster.
        # F-387.1 split: CC {1094, 1110} (Δ=15<45) + singletons {1247, 1459}.
        clusters = find_multiplet_regions(
            _FakeIdent(detected), fwhm_at,
            overlap_threshold_fwhm=10.0,
            expand_to_display_window=False,
        )
        sizes = sorted(len(c) for c in clusters)
        assert sizes == [1, 1, 2], (
            f"V2 M5: 1 CC + 2 singletons; got sizes {sizes}: "
            f"{[[round(m.library_E_keV,1) for m in c] for c in clusters]}"
        )
        cc = next(c for c in clusters if len(c) == 2)
        cc_Es = sorted(m.library_E_keV for m in cc)
        assert cc_Es == [1094.65, 1109.95]


# ──────────────────────────────────────────────────────────────────
# Synthetic: algorithmic edge cases
# ──────────────────────────────────────────────────────────────────

class TestSyntheticAlgorithmic:
    """Synthetic кейсы — алгоритмическое поведение CC + top-K cap."""

    def test_all_resolved_three_singletons(self):
        """3 widely-separated peaks (ΔE >> FWHM_avg) → 3 singletons."""
        fwhm_at = lambda ch: 10.0
        detected = (
            _FakeNuclideId("A", [_FakeLineMatch("A", 100.0, 100.0)]),
            _FakeNuclideId("B", [_FakeLineMatch("B", 200.0, 200.0)]),
            _FakeNuclideId("C", [_FakeLineMatch("C", 300.0, 300.0)]),
        )
        clusters = find_multiplet_regions(
            _FakeIdent(detected), fwhm_at,
            overlap_threshold_fwhm=20.0,  # большой, чтобы overlap не отбросил
            expand_to_display_window=False,
        )
        assert len(clusters) == 3
        for c in clusters:
            assert len(c) == 1

    def test_all_unresolved_chain_5_components_top3_cap(self):
        """5 компонент с overlapping FWHM (Δ=0.5·FWHM каждый) →
        transitive CC из всех 5 → top-3 active + 2 phantom.
        """
        fwhm_at = lambda ch: 20.0
        # Δ = 10 = 0.5·FWHM_avg < 1.0 → каждый сосед unresolved
        # Positions: 100, 110, 120, 130, 140 → 5 vertices, цепочка edges
        # I_pct: 100→50, 200→40, 300→30, 400→20, 500→10
        detected = (
            _FakeNuclideId("N1", [_FakeLineMatch("N1", 100.0, 100.0, 50.0)]),
            _FakeNuclideId("N2", [_FakeLineMatch("N2", 110.0, 110.0, 40.0)]),
            _FakeNuclideId("N3", [_FakeLineMatch("N3", 120.0, 120.0, 30.0)]),
            _FakeNuclideId("N4", [_FakeLineMatch("N4", 130.0, 130.0, 20.0)]),
            _FakeNuclideId("N5", [_FakeLineMatch("N5", 140.0, 140.0, 10.0)]),
        )
        clusters = find_multiplet_regions(
            _FakeIdent(detected), fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
        )
        # 1 CC из всех 5 (transitive)
        assert len(clusters) == 1
        assert len(clusters[0]) == 5
        active = [m for m in clusters[0] if not _is_phantom(m)]
        phantoms = [m for m in clusters[0] if _is_phantom(m)]
        assert len(active) == 3, f"top-3 active; got {len(active)}"
        assert len(phantoms) == 2
        # Top-3 by I_pct: 100 (50), 110 (40), 120 (30)
        active_Es = sorted(m.library_E_keV for m in active)
        assert active_Es == [100.0, 110.0, 120.0]
        # Phantom: 130, 140
        phantom_Es = sorted(m.library_E_keV for m in phantoms)
        assert phantom_Es == [130.0, 140.0]

    def test_transitive_closure_a_c_resolved(self):
        """A-B unresolved, B-C unresolved, A-C resolved → CC {A,B,C}."""
        fwhm_at = lambda ch: 10.0
        # FWHM_avg=10. Set Δ_AB=5, Δ_BC=5, Δ_AC=10 (=FWHM_avg, boundary →
        # strict < → resolved direct, но transitive через B).
        detected = (
            _FakeNuclideId("A", [_FakeLineMatch("A", 100.0, 100.0, 30.0)]),
            _FakeNuclideId("B", [_FakeLineMatch("B", 105.0, 105.0, 20.0)]),
            _FakeNuclideId("C", [_FakeLineMatch("C", 110.0, 110.0, 10.0)]),
        )
        clusters = find_multiplet_regions(
            _FakeIdent(detected), fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
        )
        assert len(clusters) == 1
        assert len(clusters[0]) == 3

    def test_boundary_strict_less_than(self):
        """ΔE = FWHM_avg exactly → resolved (strict-<)."""
        fwhm_at = lambda ch: 20.0
        detected = (
            _FakeNuclideId("A", [_FakeLineMatch("A", 100.0, 100.0)]),
            _FakeNuclideId("B", [_FakeLineMatch("B", 120.0, 120.0)]),
        )
        clusters = find_multiplet_regions(
            _FakeIdent(detected), fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
            unresolved_separation_fwhm_factor=1.0,  # explicit Rayleigh
        )
        # Δ=20 = FWHM_avg=20 → НЕ <  → resolved → 2 singletons
        assert sorted(len(c) for c in clusters) == [1, 1]

    def test_rayleigh_vs_f387_old_semantic_change(self):
        """Rayleigh-vs-F-387-old proof:
        FWHM_min=20, FWHM_max=30, ΔE=22.
        Old F-387:    Δ/FWHM_min = 22/20 = 1.1 > 0.7 → resolved
        New Rayleigh: Δ/FWHM_avg = 22/25 = 0.88 < 1.0 → unresolved → CC
        """
        # Asymmetric FWHM: ch 100 → 20, ch 122 → 30
        def fwhm_at(ch):
            if ch < 111:
                return 20.0
            return 30.0
        detected = (
            _FakeNuclideId("A", [_FakeLineMatch("A", 100.0, 100.0, 30.0)]),
            _FakeNuclideId("B", [_FakeLineMatch("B", 122.0, 122.0, 20.0)]),
        )
        clusters = find_multiplet_regions(
            _FakeIdent(detected), fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
        )
        # New Rayleigh → CC {A, B}
        assert len(clusters) == 1 and len(clusters[0]) == 2, (
            f"Rayleigh Δ=22 < FWHM_avg=25 → CC; got {clusters!r}"
        )

    def test_top_k_cap_intensity_sort(self):
        """5-component CC c intensities [100, 80, 60, 40, 20] →
        top-3 active {100, 80, 60}, phantom {40, 20}."""
        fwhm_at = lambda ch: 20.0
        # Все компоненты в Δ=5 → CC из 5
        detected = (
            _FakeNuclideId("N1", [_FakeLineMatch("N1", 100.0, 100.0, 100.0)]),
            _FakeNuclideId("N2", [_FakeLineMatch("N2", 105.0, 105.0, 80.0)]),
            _FakeNuclideId("N3", [_FakeLineMatch("N3", 110.0, 110.0, 60.0)]),
            _FakeNuclideId("N4", [_FakeLineMatch("N4", 115.0, 115.0, 40.0)]),
            _FakeNuclideId("N5", [_FakeLineMatch("N5", 120.0, 120.0, 20.0)]),
        )
        clusters = find_multiplet_regions(
            _FakeIdent(detected), fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
        )
        assert len(clusters) == 1
        cl = clusters[0]
        assert len(cl) == 5
        active = sorted(
            (m.library_I_pct for m in cl if not _is_phantom(m)),
            reverse=True,
        )
        phantom = sorted(
            (m.library_I_pct for m in cl if _is_phantom(m)),
            reverse=True,
        )
        assert active == [100.0, 80.0, 60.0], (
            f"active top-3 by I: expected [100,80,60], got {active}"
        )
        assert phantom == [40.0, 20.0]

    def test_factor_zero_no_split(self):
        """factor=0.0 disables filter — все компоненты остаются в одном
        cluster без split."""
        fwhm_at = lambda ch: 10.0
        detected = (
            _FakeNuclideId("A", [_FakeLineMatch("A", 100.0, 100.0)]),
            _FakeNuclideId("B", [_FakeLineMatch("B", 200.0, 200.0)]),
            _FakeNuclideId("C", [_FakeLineMatch("C", 300.0, 300.0)]),
        )
        # overlap_threshold_fwhm=200 → все 3 → один cluster, factor=0.0 → не split
        clusters = find_multiplet_regions(
            _FakeIdent(detected), fwhm_at,
            overlap_threshold_fwhm=25.0,
            expand_to_display_window=False,
            unresolved_separation_fwhm_factor=0.0,
        )
        # overlap собрал {A,B} cluster, {C} > A на 200=20·FWHM сильно отделён
        # → 1 cluster {A,B,C} только если overlap-threshold ≥ 200/10=20.
        # Реально overlap_threshold_fwhm=25 → пары: A-B(100, ≤250) yes,
        # B-C(100, ≤250) yes → CC {A,B,C}.
        # factor=0.0 → CC-split не запускается → 1 cluster size=3.
        assert len(clusters) == 1
        assert len(clusters[0]) == 3
        # Все компоненты НЕ phantom (top-K cap НЕ применён при factor=0.0)
        assert all(not _is_phantom(m) for m in clusters[0])

    def test_max_K_zero_disables_topk_cap(self):
        """max_components_per_cluster=0 → top-K cap НЕ применяется
        (back-compat / diagnostic). Все компоненты остаются active."""
        fwhm_at = lambda ch: 20.0
        detected = (
            _FakeNuclideId("N1", [_FakeLineMatch("N1", 100.0, 100.0, 100.0)]),
            _FakeNuclideId("N2", [_FakeLineMatch("N2", 105.0, 105.0, 80.0)]),
            _FakeNuclideId("N3", [_FakeLineMatch("N3", 110.0, 110.0, 60.0)]),
            _FakeNuclideId("N4", [_FakeLineMatch("N4", 115.0, 115.0, 40.0)]),
            _FakeNuclideId("N5", [_FakeLineMatch("N5", 120.0, 120.0, 20.0)]),
        )
        clusters = find_multiplet_regions(
            _FakeIdent(detected), fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
            max_components_per_cluster=0,
        )
        assert len(clusters) == 1
        assert len(clusters[0]) == 5
        # Все active (top-K cap отключён)
        assert all(not _is_phantom(m) for m in clusters[0])


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
