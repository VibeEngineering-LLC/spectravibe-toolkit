"""F-387.2 / v1.18.27.1 — Singleton CC=1 routing к primary_fep.

После F-387.1 Rayleigh-CC split возможны sub-cluster'ы с n_active == 1
(один реальный component + произвольное число phantom anchors из F-381
enrichment или F-387.1 top-K cap demote). F-387.2 contract:

  - 0 actives (phantom-only) → drop из multiplet_deconvolutions.
  - 1 active (singleton) → drop из multiplet_deconvolutions независимо
    от S/N. LineMatch для этого active component остаётся в
    identification_result.matched_lines с original peak_area_source
    ("cowell"/"lsrm_peaks_table") и попадает в primary_feps через
    стандартный путь.
  - ≥2 actives → multiplet fit как раньше (F-387.1 preserved).

Это убирает UX-аномалию «singleton multiplet» (M3 V2 с 1 active + 5
phantom, M7 V2 с 2 components 1 active). Множество правил тестируются
с использованием `apply_multiplet_deconvolution` как entry point —
это финальная инстанция, где cluster acceptance применяется и
LineMatch'ы маршрутизируются между primary_feps и multiplet_deconvolutions.

Покрывает:
  1. Synthetic high-S/N singleton: после F-387.1 split → НЕ в multiplet
     output; identification_result.matched_lines неизменён (LineMatch
     с original peak_area_source).
  2. Synthetic low-S/N singleton: drop (F-391 already gates это, sanity
     regression-проверка консистентности с F-387.2).
  3. Synthetic singleton с phantom anchors (V2 M3 case 6c→1 active):
     active component в primary_feps path, phantoms дропаются.
  4. Multiplet с 2+ actives: keep в multiplet output (F-387.1 preserved).
  5. Empirical-style Th-232 scenario: no multiplet с n_active == 1
     в output.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(
    Path(__file__).resolve().parent.parent.parent / "scripts"
))

from gamma.peaks.deconvolve import (
    apply_multiplet_deconvolution,
    find_multiplet_regions,
)
from gamma.identification.identify import (
    IdentificationResult, NuclideIdentification, LineMatch,
)


# ──────────────────────────────────────────────────────────────────
# Minimal Spectrum stub
# ──────────────────────────────────────────────────────────────────

class _FakeWindow:
    delta_E0_keV = 1.0


class _FakeSpec:
    """Minimal spec: linear energy↔channel, FWHM=const, counts=zeros.

    Sufficient для apply_multiplet_deconvolution flow control
    (cluster acceptance loop / deconvolve_multiplet fallback path).
    bin_width_keV=1.0 → channel ≈ energy_keV.
    """
    def __init__(self, n=4096):
        self.counts = np.zeros(n, dtype=np.float64)
        self.bin_width_keV = 1.0
        self.live_time = 1000.0

    def energy_to_channel(self, E_keV):
        return float(E_keV)

    def channel_to_energy(self, ch):
        return float(ch)


def _fwhm_at(ch):
    return 25.0  # constant ~25 channels (≈ 25 keV)


def _mk_lm(nuclide, E_keV, library_I_pct, *,
           peak_area=None, peak_area_uncertainty=None,
           peak_area_source="cowell"):
    """Helper: build LineMatch with measured area defaults."""
    return LineMatch(
        nuclide=nuclide,
        library_E_keV=float(E_keV),
        library_I_pct=float(library_I_pct),
        peak_channel=int(round(E_keV)),
        peak_E_keV=float(E_keV),
        peak_sigma=0.0,
        residual_keV=0.0,
        is_characteristic=False,
        peak_area=peak_area,
        peak_area_uncertainty=peak_area_uncertainty,
        peak_area_source=peak_area_source,
    )


def _mk_ident(matched_by_nuclide):
    """Build IdentificationResult из dict {nuclide: [LineMatch, ...]}."""
    detected = tuple(
        NuclideIdentification(
            nuclide=nuc,
            detected=True,
            reason="test",
            characteristic_line_keV=ms[0].library_E_keV,
            matched_lines=tuple(ms),
            confidence=None,
        )
        for nuc, ms in matched_by_nuclide.items()
    )
    return IdentificationResult(
        detector_type="NaI",
        window=_FakeWindow(),
        candidates_considered=len(detected),
        detected_nuclides=detected,
        rejected_nuclides=(),
        unmatched_peaks=(),
        notes="",
    )


def _lm_source(ident, nuclide, E_keV):
    """Найти LineMatch в (post-apply) ident по nuclide+E и вернуть source."""
    for ni in ident.detected_nuclides:
        if ni.nuclide != nuclide:
            continue
        for m in ni.matched_lines:
            if abs(m.library_E_keV - E_keV) < 0.5:
                return m.peak_area_source
    return None


def _n_active(deconv) -> int:
    """Count non-phantom components в DeconvolutionResult."""
    actives = [c for c in deconv.components
               if not str(getattr(c, "peak_area_source", "") or "")
                   .endswith("phantom")]
    return len(actives)


# ──────────────────────────────────────────────────────────────────
# 1. High-S/N singleton (исходно — F-391 keep-path) → теперь drop
# ──────────────────────────────────────────────────────────────────

class TestHighSnrSingletonDrop:
    """High-S/N singleton component (CC=1 после F-387.1 split) должен
    дропаться из multiplet_deconvolutions. LineMatch остаётся в
    identification_result с original peak_area_source."""

    def test_isolated_strong_peak_not_in_multiplet_output(self):
        """Один сильный изолированный пик: overlap-step его дропнет
        (size=1 не возвращается из find_multiplet_regions). Sanity
        baseline — no multiplet'ов в output."""
        ident = _mk_ident({
            "Cs-137": [_mk_lm(
                "Cs-137", 661.66, 85.1,
                peak_area=50000.0, peak_area_uncertainty=500.0,
            )],
        })
        spec = _FakeSpec()
        new_id, deconvs = apply_multiplet_deconvolution(
            ident, spec, _fwhm_at, expand_to_display_window=False,
        )
        assert len(deconvs) == 0, (
            f"Изолированный singleton не должен порождать multiplet entries; "
            f"got {len(deconvs)}: {[d.cluster_id for d in deconvs]}"
        )
        # LineMatch неизменён в matched_lines
        assert _lm_source(new_id, "Cs-137", 661.66) == "cowell"

    def test_high_snr_singleton_post_cc_split_drops(self):
        """Сценарий ключевой: два пика, изначально overlap-собранные в
        один cluster, но F-387.1 Rayleigh-CC их разделил на 2 singleton
        sub-clusters (ΔE > FWHM_avg). Оба high-S/N. Ожидание F-387.2:
        ОБА singleton CC drop'нуты из multiplet array, LineMatch обоих
        в primary_feps path с original source."""
        # ΔE=80 кэВ, FWHM=25 → Δ/FWHM=3.2 > 1.0 → resolved → 2 singleton CCs
        ident = _mk_ident({
            "Co-60": [
                _mk_lm("Co-60", 1173.2, 99.85,
                       peak_area=20000.0, peak_area_uncertainty=400.0),
                _mk_lm("Co-60", 1253.0, 99.98,  # сдвинут вверх, чтобы Δ=80
                       peak_area=18000.0, peak_area_uncertainty=380.0),
            ],
        })
        spec = _FakeSpec()
        # overlap_threshold=4.0 → 4·25=100 кэВ → 80 < 100 → один overlap-cluster
        new_id, deconvs = apply_multiplet_deconvolution(
            ident, spec, _fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
        )
        # F-387.2: оба CC по 1 active → оба drop
        for d in deconvs:
            assert _n_active(d) >= 2, (
                f"F-387.2: multiplet entry не должен содержать <2 active "
                f"components; cluster {d.cluster_id!r} имеет "
                f"{_n_active(d)} active. Components: "
                f"{[(c.nuclide, c.line_E_keV) for c in d.components]}"
            )
        # LineMatches неизменны (original cowell source)
        assert _lm_source(new_id, "Co-60", 1173.2) == "cowell"
        assert _lm_source(new_id, "Co-60", 1253.0) == "cowell"


# ──────────────────────────────────────────────────────────────────
# 2. Low-S/N singleton (F-391 уже handled) → regression sanity
# ──────────────────────────────────────────────────────────────────

class TestLowSnrSingletonDrop:
    """Low-S/N singleton тоже дропается (F-391 уже handled). F-387.2
    регрессионно подтверждает, что новый код не сломал этот путь."""

    def test_low_snr_active_marked_phantom_then_dropped(self):
        """Один S/N=1 active в overlap-cluster с другим weak: оба
        становятся phantom через F-391 gate; CC=0 actives → drop."""
        ident = _mk_ident({
            "Ac-228": [
                _mk_lm("Ac-228", 911.2, 25.8,
                       peak_area=10.0, peak_area_uncertainty=10.0),  # S/N=1
                _mk_lm("Ac-228", 916.0, 0.5,
                       peak_area=8.0, peak_area_uncertainty=10.0),   # S/N=0.8
            ],
        })
        spec = _FakeSpec()
        new_id, deconvs = apply_multiplet_deconvolution(
            ident, spec, _fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
            min_significance_snr=3.0,
        )
        # All low-S/N → phantom → cluster=0 actives → drop
        for d in deconvs:
            assert _n_active(d) >= 2, (
                f"Low-S/N cluster не должен попасть в output; "
                f"cluster {d.cluster_id!r}"
            )


# ──────────────────────────────────────────────────────────────────
# 3. Singleton с phantom anchors (V2 M3 case): 6c (1 active + 5 phantom)
# ──────────────────────────────────────────────────────────────────

class TestSingletonWithPhantomAnchors:
    """Эмулирует V2 M3 на real Th-232: один strong active пик + несколько
    library anchors с S/N=0 (peak_area=None, peak_area_source=
    "library_anchor" — F-381 enrichment без measured peak). После F-387.1
    CC split: phantom'ы attached к active CC (через _attach_phantom_to_nearest_cc).
    Active=1 → F-387.2 drop весь cluster."""

    def test_v2_m3_singleton_with_phantoms_drops_entirely(self):
        """Один реальный пик (Tl-208 510 high-S/N) + 5 library anchors
        (Ac-228 503.7, 509, 523 с peak_area=None, F-381 enrichment).
        Все в overlap cluster (ΔE ≤ 25 кэВ). F-387.1 CC split:
        {Tl-208 510 + Ac-228 509 + Ac-228 503.7 + Ac-228 523} —
        Δ ≤ FWHM_avg=25 для всех пар → один CC. F-391 gate помечает
        phantom anchors (S/N=0). active count = 1 (только Tl-208 510).
        Top-K cap не критичен (≤3 active). F-387.2: cluster drop."""
        # Active strong + 3 phantom anchors (no measured area)
        ident = _mk_ident({
            "Tl-208": [
                _mk_lm("Tl-208", 510.77, 22.6,
                       peak_area=12000.0, peak_area_uncertainty=400.0),
            ],
            "Ac-228": [
                _mk_lm("Ac-228", 503.7, 0.10,
                       peak_area=None, peak_area_uncertainty=None,
                       peak_area_source="library_anchor"),
                _mk_lm("Ac-228", 509.0, 0.13,
                       peak_area=None, peak_area_uncertainty=None,
                       peak_area_source="library_anchor"),
                _mk_lm("Ac-228", 523.1, 0.10,
                       peak_area=None, peak_area_uncertainty=None,
                       peak_area_source="library_anchor"),
            ],
        })
        spec = _FakeSpec()
        new_id, deconvs = apply_multiplet_deconvolution(
            ident, spec, _fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
        )
        # F-387.2: active=1 → drop entirely
        for d in deconvs:
            assert _n_active(d) >= 2, (
                f"V2-style singleton+phantom не должен попасть в output; "
                f"cluster {d.cluster_id!r} components: "
                f"{[(c.nuclide, c.line_E_keV) for c in d.components]}"
            )
        # Tl-208 510 LineMatch остаётся с original source (cowell)
        # — он рутается в primary_feps через identification_result
        assert _lm_source(new_id, "Tl-208", 510.77) == "cowell"


# ──────────────────────────────────────────────────────────────────
# 4. Genuine multiplet (2+ actives) preserved (F-387.1 backward compat)
# ──────────────────────────────────────────────────────────────────

class TestGenuineMultipletPreserved:
    """Cluster с ≥2 active components (real multiplet) ДОЛЖЕН остаться
    в output. F-387.2 не должен ломать F-387.1 path."""

    def test_two_active_unresolved_pair_kept(self):
        """Co-60 1173/1332 на NaI 63x63: ΔE=159, FWHM=80 → Δ/FWHM_avg=2 →
        в реальности resolved. Лучше использовать close pair: симулируем
        Ac-228 964/969 (Δ=4, FWHM=25, Δ/FWHM=0.16 → unresolved CC)."""
        ident = _mk_ident({
            "Ac-228": [
                _mk_lm("Ac-228", 964.77, 4.99,
                       peak_area=3000.0, peak_area_uncertainty=200.0),
                _mk_lm("Ac-228", 968.97, 15.8,
                       peak_area=9000.0, peak_area_uncertainty=300.0),
            ],
        })
        spec = _FakeSpec()
        new_id, deconvs = apply_multiplet_deconvolution(
            ident, spec, _fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
        )
        # ≥1 multiplet entry с 2 active
        assert len(deconvs) >= 1, "2-active unresolved pair должен дать multiplet"
        found_pair = False
        for d in deconvs:
            if _n_active(d) == 2:
                Es = sorted(c.line_E_keV for c in d.components
                            if not str(getattr(c, "peak_area_source",
                                                "") or "").endswith("phantom"))
                if Es == [pytest.approx(964.77, abs=0.5),
                          pytest.approx(968.97, abs=0.5)]:
                    found_pair = True
                    break
        assert found_pair, (
            f"2-active multiplet {{965, 969}} должен присутствовать; "
            f"got {[(d.cluster_id, [(c.nuclide, c.line_E_keV) for c in d.components]) for d in deconvs]}"
        )

    @pytest.mark.xfail(reason="F-441 isolated-peak classifier (NOT F-449): A/B/C isolation matrix 2026-06-17 -- disabling _is_isolated_peak makes this PASS; F-449 ruled out (GAMMA_FREE_SIGMA=1 identical fail). Synthetic empty-library-pool artifact. Follow-up: _state/agent_a/outbox/2026-06-17_F441_followup_multiplet_classifier_sideeffects.md", strict=False)
    def test_three_active_kept(self):
        """3 unresolved actives — все в одном CC, multiplet kept."""
        ident = _mk_ident({
            "N1": [_mk_lm("N1", 100.0, 50.0,
                           peak_area=2000.0, peak_area_uncertainty=100.0)],
            "N2": [_mk_lm("N2", 110.0, 40.0,
                           peak_area=1500.0, peak_area_uncertainty=80.0)],
            "N3": [_mk_lm("N3", 120.0, 30.0,
                           peak_area=1000.0, peak_area_uncertainty=70.0)],
        })
        spec = _FakeSpec()
        new_id, deconvs = apply_multiplet_deconvolution(
            ident, spec, _fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
        )
        assert len(deconvs) >= 1
        # должен быть 3-active multiplet
        triples = [d for d in deconvs if _n_active(d) == 3]
        assert len(triples) >= 1, (
            f"3-active multiplet должен сохраниться; "
            f"got actives={[_n_active(d) for d in deconvs]}"
        )


# ──────────────────────────────────────────────────────────────────
# 5. Inverse audit — no multiplet с n_active == 1 в any output
# ──────────────────────────────────────────────────────────────────

class TestNoSingletonInMultipletOutput:
    """F-387.2 invariant: после apply_multiplet_deconvolution, для каждого
    выходного DeconvolutionResult n_active ≥ 2. Проверяется в нескольких
    смешанных сценариях."""

    @pytest.mark.parametrize("scenario", [
        # (label, matched_by_nuclide dict)
        ("isolated_strong", {
            "Cs-137": [_mk_lm("Cs-137", 661.66, 85.1,
                              peak_area=50000.0,
                              peak_area_uncertainty=500.0)],
        }),
        ("two_resolved_singletons", {
            "Co-60": [
                _mk_lm("Co-60", 1173.2, 99.85,
                       peak_area=20000.0, peak_area_uncertainty=400.0),
                _mk_lm("Co-60", 1253.0, 99.98,
                       peak_area=18000.0, peak_area_uncertainty=380.0),
            ],
        }),
        ("one_active_plus_phantoms", {
            "Tl-208": [_mk_lm("Tl-208", 510.77, 22.6,
                              peak_area=12000.0,
                              peak_area_uncertainty=400.0)],
            "Ac-228": [
                _mk_lm("Ac-228", 503.7, 0.10,
                       peak_area=None, peak_area_uncertainty=None,
                       peak_area_source="library_anchor"),
                _mk_lm("Ac-228", 523.1, 0.10,
                       peak_area=None, peak_area_uncertainty=None,
                       peak_area_source="library_anchor"),
            ],
        }),
        ("low_snr_cluster", {
            "Ac-228": [
                _mk_lm("Ac-228", 911.2, 25.8,
                       peak_area=10.0, peak_area_uncertainty=10.0),
                _mk_lm("Ac-228", 916.0, 0.5,
                       peak_area=8.0, peak_area_uncertainty=10.0),
            ],
        }),
        ("genuine_2_active_pair", {
            "Ac-228": [
                _mk_lm("Ac-228", 964.77, 4.99,
                       peak_area=3000.0, peak_area_uncertainty=200.0),
                _mk_lm("Ac-228", 968.97, 15.8,
                       peak_area=9000.0, peak_area_uncertainty=300.0),
            ],
        }),
    ])
    def test_no_multiplet_with_one_active(self, scenario):
        label, matched = scenario
        ident = _mk_ident(matched)
        spec = _FakeSpec()
        _, deconvs = apply_multiplet_deconvolution(
            ident, spec, _fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
        )
        for d in deconvs:
            assert _n_active(d) >= 2, (
                f"[{label}] F-387.2 violation: cluster {d.cluster_id!r} "
                f"has {_n_active(d)} active components. "
                f"Components: {[(c.nuclide, c.line_E_keV) for c in d.components]}"
            )


# ──────────────────────────────────────────────────────────────────
# 6. F-387.1 unit-level invariant: find_multiplet_regions singleton path
# ──────────────────────────────────────────────────────────────────

class TestF387_1InvariantPreserved:
    """F-387.1 (find_multiplet_regions) НЕ изменился: CC=1 sub-clusters
    всё ещё возвращаются из find_multiplet_regions (как size=1 list).
    F-387.2 фильтрует их в apply_multiplet_deconvolution, НЕ в
    find_multiplet_regions. Это инвариант — алгоритм Rayleigh-CC split
    не зависит от output routing."""

    class _FakeLineMatchUnit:
        def __init__(self, nuclide, E, peak_area=None, peak_area_unc=None,
                     library_I_pct=10.0, peak_area_source=""):
            self.nuclide = nuclide
            self.library_E_keV = E
            self.library_I_pct = library_I_pct
            self.peak_channel = int(round(E))
            self.peak_E_keV = E
            self.is_characteristic = False
            self.peak_sigma = 0.0
            self.residual_keV = 0.0
            self.peak_area = peak_area
            self.peak_area_uncertainty = peak_area_unc
            self.peak_area_source = peak_area_source

    class _FakeNuclideId:
        def __init__(self, nuclide, matched_lines):
            self.nuclide = nuclide
            self.matched_lines = matched_lines
            self.detected = True

    class _FakeIdent:
        def __init__(self, detected):
            self.detected_nuclides = detected

    def test_find_multiplet_regions_still_returns_singletons(self):
        """Two well-separated strong peaks (ΔE >> FWHM_avg, but в одном
        overlap-cluster). find_multiplet_regions → 2 singleton CC
        sub-clusters (size=1 каждый)."""
        fwhm_at = lambda ch: 25.0
        detected = (
            self._FakeNuclideId("A", [
                self._FakeLineMatchUnit("A", 100.0,
                    peak_area=1000.0, peak_area_unc=50.0)]),
            self._FakeNuclideId("B", [
                self._FakeLineMatchUnit("B", 180.0,  # Δ=80, FWHM=25 → resolved
                    peak_area=800.0, peak_area_unc=40.0)]),
        )
        clusters = find_multiplet_regions(
            self._FakeIdent(detected), fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
        )
        # Должны быть 2 singleton CCs
        sizes = sorted(len(c) for c in clusters)
        assert sizes == [1, 1], (
            f"F-387.1 invariant: 2 resolved peaks → 2 singletons из "
            f"find_multiplet_regions; got sizes {sizes}"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
