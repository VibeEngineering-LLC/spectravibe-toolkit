"""F-387 / v1.18.26 → F-387.1 / v1.18.26.1 — Rayleigh-CC unresolved
criterion. Регрессионный snapshot-набор.

Контекст:
    После F-374 (expand_to_display_window) + F-381 (library-anchor
    enrichment) auto-clusters содержат 7-13 компонент с separations
    5-20 кэВ при FWHM 22-30 кэВ.

    **Семантика v1.18.26.1 (F-387.1):**
    Pair (a, b) unresolved ⟺ |ΔE| < factor · FWHM_avg(a, b),
    где FWHM_avg = (FWHM_a + FWHM_b)/2. Default factor=1.0 = Rayleigh.
    Старая v1.18.26 (F-387) использовала factor · FWHM_min — другая
    физика, более грубое приближение.

    F-387.1 также разбивает cluster на connected components через граф:
    unresolved pairs образуют edges, BFS-CC даёт sub-cluster'ы.
    Изолированные vertices становятся **singleton sub-cluster'ами**
    (cluster size=1 разрешён — downstream `deconvolve_multiplet`
    обрабатывает как trivial 1-component fit).

Покрытие:
    - PASS: M1 Ac-228 911 + 964.8 + 969 (964.8+969 пара 4 кэВ при
      FWHM 30 → 0.13·FWHM_avg < 1.0) → unresolved sub-cluster {964.8, 969}.
      911 изолирована: Δ=53 кэВ → 1.78·FWHM_avg > 1.0 → singleton.
    - SPLIT: synthetic «M3-prod» Ac-228 503.7 + Tl-208 583.19
      (separation 79.5 кэВ при FWHM 25 → 3.18·FWHM_avg) → 2 singletons.
    - SPLIT: synthetic 100+150 кэВ FWHM 15 → 2 singletons.
    - factor=0.0 отключает фильтр (back-compat diagnostic mode).
    - Mirror v2: detect_multiplet_clusters симметрично.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

from gamma.peaks.deconvolve import find_multiplet_regions


# ──────────────────────────────────────────────────────────────────
# Test fixtures — minimal LineMatch / IdentificationResult stand-ins
# ──────────────────────────────────────────────────────────────────

class _FakeLineMatch:
    """Stand-in for LineMatch — поля, читаемые find_multiplet_regions."""
    # NOTE: НЕ slots — тест-fake должен быть mutable (F-387.1 пишет
    # peak_area_source при top-K cap через setattr fallback).
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


# ──────────────────────────────────────────────────────────────────
# F-387.1 — find_multiplet_regions Rayleigh-CC + top-K
# ──────────────────────────────────────────────────────────────────

class TestF387FindMultipletRegionsCriterion:
    """Rayleigh: unresolved ⟺ ∃ пара компонент с |Δ| < factor·FWHM_avg.
    F-387.1 разбивает cluster на CCs; default factor=1.0."""

    def test_M1_ac228_911_965_969_yields_unresolved_pair_plus_singleton(self):
        """M1 hard-locked: Ac-228 911 + 964.8 + 969.
        964.8+969 пара = 4.2 кэВ при FWHM 30 → 0.14·FWHM_avg < 1.0
        → unresolved → CC {964.8, 969}. 911 отделена от 964.8 на
        53.6 кэВ (1.79·FWHM_avg > 1.0) → singleton CC {911}.

        Expected: 2 sub-cluster'а (1 unresolved-пара + 1 singleton).
        """
        fwhm_at = lambda ch: 30.0
        detected = (
            _FakeNuclideId("Ac-228", [
                _FakeLineMatch("Ac-228", 911.2, 911.2),
                _FakeLineMatch("Ac-228", 964.8, 964.8),
                _FakeLineMatch("Ac-228", 969.0, 969.0),
            ]),
        )
        clusters = find_multiplet_regions(
            _FakeIdent(detected), fwhm_at,
            overlap_threshold_fwhm=4.0,  # принудительно соберём в overlap-этапе
            expand_to_display_window=False,
        )
        # 1 cluster from overlap-stage → 2 CCs after F-387.1 split:
        # {911}, {964.8, 969}
        assert len(clusters) == 2, (
            f"M1 must split into 2 sub-clusters (CC singleton 911 + "
            f"CC pair 964.8+969); got {len(clusters)}: "
            f"{[[m.library_E_keV for m in c] for c in clusters]}"
        )
        # Найдём pair-CC
        pair_cc = next(
            (c for c in clusters if len(c) == 2), None,
        )
        singleton_cc = next(
            (c for c in clusters if len(c) == 1), None,
        )
        assert pair_cc is not None, "expected pair-CC"
        assert singleton_cc is not None, "expected singleton-CC for 911"
        pair_Es = sorted(m.library_E_keV for m in pair_cc)
        assert pair_Es == [964.8, 969.0], (
            f"pair-CC must be {964.8, 969.0}; got {pair_Es}"
        )
        assert singleton_cc[0].library_E_keV == 911.2

    def test_M3_production_504_583_splits_into_two_singletons(self):
        """Synthetic «M3-prod»: Ac-228 503.7 + Tl-208 583.19.
        Separation 79.5 кэВ при FWHM 25 → 3.18·FWHM_avg > 1.0.
        F-387.1 → 2 singleton sub-cluster'а (раньше F-387 reject'ил).
        """
        fwhm_at = lambda ch: 25.0
        detected = (
            _FakeNuclideId("Ac-228", [
                _FakeLineMatch("Ac-228", 503.7, 503.7),
            ]),
            _FakeNuclideId("Tl-208", [
                _FakeLineMatch("Tl-208", 583.19, 583.19),
            ]),
        )
        clusters = find_multiplet_regions(
            _FakeIdent(detected), fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
        )
        # F-387.1 SPLIT: 2 singletons (вместо F-387 reject)
        assert len(clusters) == 2, (
            f"M3-prod должен split на 2 singletons (sep 79.5 > "
            f"1.0·FWHM_avg=25); got {len(clusters)}: "
            f"{[[m.library_E_keV for m in c] for c in clusters]}"
        )
        sizes = sorted(len(c) for c in clusters)
        assert sizes == [1, 1], f"expected [1,1], got {sizes}"

    def test_synthetic_100_150_splits_into_two_singletons(self):
        """Synthetic пара: 100 + 150 кэВ при FWHM 15 → sep 50 кэВ
        = 3.33·FWHM_avg → 2 singletons.
        """
        fwhm_at = lambda ch: 15.0
        detected = (
            _FakeNuclideId("X", [_FakeLineMatch("X", 100.0, 100.0)]),
            _FakeNuclideId("Y", [_FakeLineMatch("Y", 150.0, 150.0)]),
        )
        clusters = find_multiplet_regions(
            _FakeIdent(detected), fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
        )
        assert len(clusters) == 2, (
            f"100+150 → 2 singletons; got {clusters!r}"
        )
        assert sorted(len(c) for c in clusters) == [1, 1]

    def test_factor_zero_disables_filter(self):
        """factor=0.0 → CC-split не запускается; cluster выходит
        cohesively как одна группа всех overlap-собранных компонент.
        """
        fwhm_at = lambda ch: 15.0
        detected = (
            _FakeNuclideId("X", [_FakeLineMatch("X", 100.0, 100.0)]),
            _FakeNuclideId("Y", [_FakeLineMatch("Y", 150.0, 150.0)]),
        )
        clusters = find_multiplet_regions(
            _FakeIdent(detected), fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
            unresolved_separation_fwhm_factor=0.0,
        )
        # factor=0.0 → F-387.1 фильтр выключен, возвращается original cluster
        # (overlap-этап собрал оба → 1 cluster size=2).
        assert len(clusters) == 1, (
            f"factor=0.0 → фильтр отключён, оригинальный cluster; "
            f"got {clusters!r}"
        )
        assert len(clusters[0]) == 2

    def test_M2_ac228_1588_bi212_1620_ac228_1630_survives_via_transitive(self):
        """M2: Ac-228 1588 + Bi-212 1620 + Ac-228 1630, FWHM=37.
        Edges (Rayleigh, FWHM_avg=37):
          1588-1620: Δ=32 → 0.86·FWHM < 1.0 → unresolved
          1620-1630: Δ=10 → 0.27·FWHM < 1.0 → unresolved
          1588-1630: Δ=42 → 1.14·FWHM > 1.0 → resolved
        Transitive через 1620 → CC {1588, 1620, 1630}.
        """
        fwhm_at = lambda ch: 37.0
        detected = (
            _FakeNuclideId("Ac-228", [
                _FakeLineMatch("Ac-228", 1588.2, 1588.2),
                _FakeLineMatch("Ac-228", 1630.6, 1630.6),
            ]),
            _FakeNuclideId("Bi-212", [
                _FakeLineMatch("Bi-212", 1620.5, 1620.5),
            ]),
        )
        clusters = find_multiplet_regions(
            _FakeIdent(detected), fwhm_at,
            overlap_threshold_fwhm=1.5,
            expand_to_display_window=False,
        )
        assert len(clusters) == 1, (
            f"M2 transitive CC {{1588, 1620, 1630}}; got {clusters!r}"
        )
        assert len(clusters[0]) == 3

    def test_default_factor_is_1_0_rayleigh_boundary(self):
        """Sanity: default factor=1.0. Boundary `<` strict.
        FWHM 20 для обоих компонент → FWHM_avg=20.
        Δ=20 кэВ = 1.0·FWHM_avg → строго resolved.
        Δ=18 кэВ = 0.9·FWHM_avg → unresolved.
        """
        fwhm_at = lambda ch: 20.0
        # boundary: Δ=20 = FWHM_avg → resolved → singletons
        detected_at_boundary = (
            _FakeNuclideId("A", [_FakeLineMatch("A", 100.0, 100.0)]),
            _FakeNuclideId("B", [_FakeLineMatch("B", 120.0, 120.0)]),
        )
        clusters = find_multiplet_regions(
            _FakeIdent(detected_at_boundary), fwhm_at,
            overlap_threshold_fwhm=2.0,
            expand_to_display_window=False,
            unresolved_separation_fwhm_factor=1.0,  # explicit Rayleigh boundary test
        )
        # 2 singleton CCs (strict-`<`)
        assert sorted(len(c) for c in clusters) == [1, 1], (
            f"boundary Δ=FWHM_avg → 2 singletons (strict-<); got {clusters!r}"
        )
        # below boundary
        detected_below = (
            _FakeNuclideId("A", [_FakeLineMatch("A", 100.0, 100.0)]),
            _FakeNuclideId("B", [_FakeLineMatch("B", 118.0, 118.0)]),
        )
        clusters = find_multiplet_regions(
            _FakeIdent(detected_below), fwhm_at,
            overlap_threshold_fwhm=2.0,
            expand_to_display_window=False,
            unresolved_separation_fwhm_factor=1.0,  # explicit Rayleigh boundary test
        )
        assert len(clusters) == 1 and len(clusters[0]) == 2, (
            f"Δ=18 < FWHM_avg=20 → unresolved CC; got {clusters!r}"
        )


# ──────────────────────────────────────────────────────────────────
# F-387.1 — peak_pipeline_v2.detect_multiplet_clusters mirror
# ──────────────────────────────────────────────────────────────────

class TestF387V2DetectMultipletClustersMirror:
    """Зеркальное поведение в experimental v2 pipeline."""

    def _make_hit(self, ch, E_keV, fwhm=20.0):
        from gamma.experimental.peak_pipeline_v2 import PeakHit
        return PeakHit(
            channel=ch, energy_keV=E_keV, fwhm_channels=fwhm,
            significance=20.0, source="both",
            mari_significance=20.0, conv_significance=20.0,
        )

    def test_v2_resolved_doublet_splits_into_singletons(self):
        """Synthetic 100+150 кэВ FWHM 15 → F-387.1 split на singletons.

        F-387.2 / v1.18.27.1 update: после Rayleigh-CC split резолвимые
        sub-cluster'ы size=1 (singletons) дропаются из multiplet output
        — они маршрутизируются в primary_feps через standard LineMatch
        path. До v1.18.27.1 возвращалось 2 singleton clusters,
        теперь — пустой список (singletons НЕ multiplets).

        NB: V2 upstream pick-grouping строит cluster только когда ≥2
        library lines попадают в ROI одного peak hit. Используем hit
        в центре (125) с roi_extend=2.0 → ROI [95..155] захватит оба
        library lines → cluster size=2 → F-387.1 split на singletons →
        F-387.2 drop.
        """
        from gamma.experimental.peak_pipeline_v2 import detect_multiplet_clusters
        lib = {
            "X": [(100.0, 50.0)],
            "Y": [(150.0, 50.0)],
        }
        fwhm_keV_at = lambda E: 15.0
        hits = [self._make_hit(125, 125.0)]
        clusters = detect_multiplet_clusters(
            hits, fwhm_keV_at,
            library=lib,
            roi_extend_fwhm=2.0,
        )
        # F-387.2: singletons → drop. Output should be empty
        # (1 cluster → CC split → 2 singletons → both dropped).
        assert len(clusters) == 0, (
            f"F-387.2: resolved doublet → 2 singletons → drop; "
            f"got {clusters!r}"
        )

    def test_v2_unresolved_pair_survives(self):
        """Synthetic 964.8+969 кэВ FWHM 30 → cluster выживает."""
        from gamma.experimental.peak_pipeline_v2 import detect_multiplet_clusters
        lib = {
            "Ac-228": [(964.77, 4.99), (968.97, 15.8)],
        }
        fwhm_keV_at = lambda E: 30.0
        hits = [self._make_hit(967, 967.0)]
        clusters = detect_multiplet_clusters(
            hits, fwhm_keV_at,
            library=lib,
            roi_extend_fwhm=2.5,
        )
        assert len(clusters) == 1 and clusters[0].n_components == 2, (
            f"v2 mirror 964.8+969 (sep 4.2 < FWHM_avg=30) → CC; "
            f"got {clusters!r}"
        )

    def test_v2_factor_zero_disables_filter(self):
        """factor=0.0 → back-compat: v2 не фильтрует."""
        from gamma.experimental.peak_pipeline_v2 import detect_multiplet_clusters
        lib = {
            "X": [(100.0, 50.0)],
            "Y": [(150.0, 50.0)],
        }
        fwhm_keV_at = lambda E: 15.0
        hits = [self._make_hit(125, 125.0)]
        clusters_off = detect_multiplet_clusters(
            hits, fwhm_keV_at,
            library=lib,
            roi_extend_fwhm=2.0,
            unresolved_separation_fwhm_factor=0.0,
        )
        clusters_on = detect_multiplet_clusters(
            hits, fwhm_keV_at,
            library=lib,
            roi_extend_fwhm=2.0,
        )
        assert len(clusters_off) == 1 and clusters_off[0].n_components == 2, (
            f"factor=0.0 должно отключать F-387.1 в v2; got {clusters_off!r}"
        )
        # F-387.1 split → 2 singletons → F-387.2 drop → empty.
        # (Раньше: 2 singleton clusters в output; v1.18.27.1: singletons
        # маршрутизируются в primary_feps, не в multiplet array).
        assert len(clusters_on) == 0, (
            f"F-387.2: default factor=1.1 split на singletons → drop; "
            f"got {clusters_on!r}"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
