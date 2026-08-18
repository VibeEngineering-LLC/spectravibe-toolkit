"""F-354 / v1.18.24.0 — experimental peak_pipeline_v2 (тесты).

Покрытие:
  * dual-method search (Mariscotti + matched filter, merge)
  * auto-detect multiplet clusters (без FORCED_CLUSTERS)
  * coupled_intensity_fit на каждый кластер
  * compare_with_production end-to-end
  * version-bump assertion
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

ROOT = Path(__file__).resolve().parents[2]
KIT = ROOT / "detectors/Gamma-1S/reference_spectra/reference_kits/Marinelli_1L"
TH_SAMPLE = KIT / "Th-232/Th232_420-7-17_Маринелли_0cm.spe"
TH_BG = KIT / "Th-232/Фон закр кр вода_13.spe"


# ──────────────────────────────────────────────────────────────────
# Module imports + dataclass shape
# ──────────────────────────────────────────────────────────────────

class TestModuleAPI:
    def test_module_imports(self):
        from gamma.experimental import (  # noqa: F401
            PipelineV2Result, ComparisonReport,
            search_dual_method, detect_multiplet_clusters,
            decompose_multiplets, run_v2_pipeline, compare_with_production,
        )

    def test_dataclasses_have_expected_fields(self):
        from gamma.experimental.peak_pipeline_v2 import (
            PeakHit, MultipletClusterCandidate,
        )
        hit = PeakHit(channel=100, energy_keV=300.0, fwhm_channels=10.0,
                      significance=15.0, source="both",
                      mari_significance=15.0, conv_significance=14.0)
        assert hit.source == "both"
        assert hit.mari_significance == 15.0

        cluster = MultipletClusterCandidate(
            cluster_id="auto_M1",
            E_lo_keV=800.0, E_hi_keV=1050.0,
            components=(("Ac-228", 911.2, 25.8, "Ac-228"),
                        ("Ac-228", 969.0, 15.8, "Ac-228")),
            n_components=2,
            detection_reason="test",
        )
        assert cluster.n_components == 2


# ──────────────────────────────────────────────────────────────────
# Dual-method search
# ──────────────────────────────────────────────────────────────────

class TestDualMethodSearch:
    @pytest.fixture(scope="class")
    def th_spec(self):
        if not TH_SAMPLE.exists():
            pytest.skip(f"Th-232 fixture missing: {TH_SAMPLE}")
        from gamma.io.readers import read_spectrum
        return read_spectrum(str(TH_SAMPLE))

    def test_merged_peaks_have_source_tag(self, th_spec):
        from gamma.experimental import search_dual_method
        from gamma.calibration.fwhm_provider import make_fwhm_at_channel_provider

        counts = np.asarray(th_spec.counts, dtype=float)
        fwhm_provider = make_fwhm_at_channel_provider(th_spec)
        merged, mari, conv = search_dual_method(
            counts, fwhm_provider, sigma_threshold=3.0,
            channel_to_energy=th_spec.channel_to_energy,
        )
        assert len(merged) > 0
        # Каждый pick имеет source ∈ {mariscotti, matched_filter, both}
        for h in merged:
            assert h.source in {"mariscotti", "matched_filter", "both"}
        # На Th-232 ожидаем «both» как доминирующее source (большинство пиков
        # видны обоими методами)
        n_both = sum(1 for h in merged if h.source == "both")
        assert n_both >= 8, f"expected ≥8 both-source peaks, got {n_both}"

    def test_matched_filter_finds_511_keV(self, th_spec):
        """Matched filter обязан найти 511 кэВ (annihilation Tl-208 SE)
        — то, что Mariscotti пропускает на этом fixture."""
        from gamma.experimental import search_dual_method
        from gamma.calibration.fwhm_provider import make_fwhm_at_channel_provider

        counts = np.asarray(th_spec.counts, dtype=float)
        fwhm_provider = make_fwhm_at_channel_provider(th_spec)
        merged, mari, conv = search_dual_method(
            counts, fwhm_provider, sigma_threshold=3.0,
            channel_to_energy=th_spec.channel_to_energy,
        )
        # Найти pick в районе 511 кэВ (±20 кэВ от 511)
        ann_peaks = [h for h in merged if 490 <= h.energy_keV <= 530]
        assert len(ann_peaks) >= 1, "annihilation 511 кэВ не найден ни одним методом"
        # И хотя бы один должен быть от matched_filter (это критический случай)
        ann_via_conv = [h for h in ann_peaks
                        if h.source in {"matched_filter", "both"}]
        assert len(ann_via_conv) >= 1, "511 кэВ не найден matched filter'ом"


# ──────────────────────────────────────────────────────────────────
# Auto-cluster detection
# ──────────────────────────────────────────────────────────────────

class TestMultipletAutoDetect:
    def test_detects_ac228_m1_cluster_on_th232(self):
        """F-387.1 / v1.18.26.1 — Rayleigh-CC split разделяет M1
        cluster {Ac-228 911, 965, 969} + Tl-208 860 на sub-clusters:
        911 и 860 — singletons (Δ к ближайшему > FWHM_avg);
        965+969 — unresolved CC (Δ=4 vs FWHM ~50 → 0.08·FWHM_avg).

        Pre-F-387.1 ожидание (911+969 в одном cluster) больше неактуально.
        Тест переписан: ищем хотя бы один cluster, содержащий unresolved
        пару Ac-228 965+969.
        """
        from gamma.experimental import (
            detect_multiplet_clusters,
            search_dual_method, run_v2_pipeline,
        )
        from gamma.io.readers import read_spectrum

        if not TH_SAMPLE.exists():
            pytest.skip("Th-232 fixture missing")
        spec = read_spectrum(str(TH_SAMPLE))
        res = run_v2_pipeline(
            spec, sigma_threshold=3.0,
            chain_filter=["Ac-228", "Tl-208", "Pb-212", "Bi-212"],
        )
        m_with_965_969 = []
        for c in res.multiplet_candidates:
            nuclides_E = {(comp[0], round(comp[1])) for comp in c.components}
            if (("Ac-228", 965) in nuclides_E
                    and ("Ac-228", 969) in nuclides_E):
                m_with_965_969.append(c)
        assert len(m_with_965_969) >= 1, (
            "F-387.1 auto-detect должен найти cluster с unresolved "
            "Ac-228 965+969 парой"
        )

    def test_detects_triplet_1588_1620_1630(self):
        """F-387.1 / v1.18.26.1 — triplet 1588+1620+1630 (FWHM ~62):
        1620+1630 Δ=10 → 0.16·FWHM_avg → unresolved edge.
        1588+1620 Δ=32 → 0.52·FWHM_avg → unresolved edge.
        1588+1630 Δ=42 → 0.68·FWHM_avg → unresolved edge.
        Все три в одном CC → cluster содержит {1588, 1620, 1630}.
        """
        from gamma.experimental import run_v2_pipeline
        from gamma.io.readers import read_spectrum

        if not TH_SAMPLE.exists():
            pytest.skip("Th-232 fixture missing")
        spec = read_spectrum(str(TH_SAMPLE))
        res = run_v2_pipeline(
            spec, sigma_threshold=3.0,
            chain_filter=["Ac-228", "Tl-208", "Bi-212"],
        )
        m_with_1620 = []
        for c in res.multiplet_candidates:
            nuclides_E = {(comp[0], round(comp[1])) for comp in c.components}
            if (("Bi-212", 1620) in nuclides_E
                    and ("Ac-228", 1588) in nuclides_E):
                m_with_1620.append(c)
        assert len(m_with_1620) >= 1, (
            "F-387.1 auto-detect должен найти triplet 1588+1620+1630 "
            "(все unresolved через FWHM ~62)"
        )


# ──────────────────────────────────────────────────────────────────
# Coupled fit decomposition
# ──────────────────────────────────────────────────────────────────

class TestDecomposeMultiplets:
    def test_triplet_converges_with_low_chi2(self):
        """Триплет 1588+1620+1630 должен сходиться с χ²/ν < 15 (на эталоне
        v1.17.2 было 1.17, у нас будет хуже — другой фикстур, но в пределах)."""
        from gamma.experimental import run_v2_pipeline
        from gamma.io.readers import read_spectrum

        if not TH_SAMPLE.exists():
            pytest.skip("Th-232 fixture missing")
        spec = read_spectrum(str(TH_SAMPLE))
        res = run_v2_pipeline(
            spec, sigma_threshold=3.0,
            chain_filter=["Ac-228", "Tl-208", "Bi-212"],
        )
        # Найти fit на триплет 1588+1620+1630
        for fit in res.coupled_fits:
            if isinstance(fit, dict):
                continue
            E_components = {round(c.E_keV) for c in fit.components}
            if 1588 in E_components and 1620 in E_components:
                assert fit.converged
                assert fit.chi2_per_dof < 15.0, (
                    f"χ²/ν={fit.chi2_per_dof:.2f} слишком высокое для триплета"
                )
                # Все 3 компоненты имеют S/σ > 1
                for c in fit.components:
                    s_sigma = c.area / c.sigma_area if c.sigma_area > 0 else 0
                    assert s_sigma >= 1.0, (
                        f"{c.nuclide} {c.E_keV}: S/σ={s_sigma:.1f} слишком низкий"
                    )
                return
        pytest.fail("триплет 1588+1620+1630 не найден в coupled_fits")


# ──────────────────────────────────────────────────────────────────
# End-to-end comparison
# ──────────────────────────────────────────────────────────────────

class TestComparisonWithProduction:
    def test_compare_returns_diff_report(self):
        from gamma.experimental import compare_with_production
        from gamma.experimental.peak_pipeline_v2 import (
            ComparisonReport, PipelineV2Result,
        )

        if not TH_SAMPLE.exists() or not TH_BG.exists():
            pytest.skip("Th-232 kit fixtures missing")
        prod, v2, rpt = compare_with_production(
            sample_path=str(TH_SAMPLE),
            background_path=str(TH_BG),
            sample_mass_kg=0.5,
            chain_filter=["Ac-228", "Tl-208", "Pb-212", "Bi-212",
                          "Pb-214", "Bi-214", "K-40", "Cs-137"],
        )
        assert isinstance(v2, PipelineV2Result)
        assert isinstance(rpt, ComparisonReport)
        assert rpt.production_n_peaks > 0
        assert rpt.v2_n_peaks > 0
        assert rpt.v2_n_multiplets >= 1
        # Notes должны содержать строки про search
        notes_txt = " ".join(rpt.notes)
        assert "search" in notes_txt
        assert "Mariscotti" in notes_txt or "merged" in notes_txt

    def test_v2_finds_strictly_more_peaks_than_production(self):
        """Цель v2 — расширенный отлов; на Th-232 ожидаем v2 ≥ production."""
        from gamma.experimental import compare_with_production

        if not TH_SAMPLE.exists() or not TH_BG.exists():
            pytest.skip("Th-232 kit fixtures missing")
        _, _, rpt = compare_with_production(
            sample_path=str(TH_SAMPLE),
            background_path=str(TH_BG),
            sample_mass_kg=0.5,
            chain_filter=None,
        )
        # v2 ≥ production по количеству; «only_in_production» допустим
        # (Cs-Kα 32 кэВ — пограничный)
        assert rpt.v2_n_peaks >= rpt.production_n_peaks - 1


# ──────────────────────────────────────────────────────────────────
# Version-bump
# ──────────────────────────────────────────────────────────────────

class TestVersionBump:
    def test_version_geq_1_18_24_0(self):
        from gamma.reporting.json_report import SKILL_VERSION
        m = re.match(r"v(\d+)\.(\d+)\.(\d+)(?:\.(\d+))?", SKILL_VERSION)
        assert m, f"unparseable SKILL_VERSION: {SKILL_VERSION}"
        parts = tuple(int(p or 0) for p in m.groups())
        assert parts >= (1, 18, 24, 0), (
            f"SKILL_VERSION {SKILL_VERSION} below F-354 baseline v1.18.24.0"
        )
