# -*- coding: utf-8 -*-
"""F-398 / v1.18.28 — tests for scripts/run_skill.py orchestration.

Two tiers:
  • Fast unit tests (no pipeline) — config merge, metadata extraction,
    bundle layout, phase markers, index rendering, exit-code contract on
    bad input. ≈ runs in <1s.
  • Slow integration tests (`@pytest.mark.slow`) — full Th-232 pipeline,
    --resume, --batch, idempotency. Run with `pytest -m slow`.

Покрывает §5.6 AGENTS.md (smoke, resume, config override, batch, error
path, idempotency).
"""
from __future__ import annotations

import hashlib
import importlib
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict

import pytest

# ──────────────────────────────────────────────────────────────────
# Make scripts/ importable (conftest already does, but be defensive
# in case the test is invoked stand-alone).
# ──────────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

# Import the module under test. `scripts/run_skill.py` is at top of scripts/,
# so a plain import after sys.path tweak resolves it.
run_skill = importlib.import_module("run_skill")

FIXTURES = REPO_ROOT / "evals" / "fixtures"
CS137_FIXTURE = FIXTURES / "M_cs_легкий_2001-2005.spe"
K40_FIXTURE = FIXTURES / "M_k_легкий_2001-2005.spe"
TH232_FIXTURE = FIXTURES / "M_th_легкий_2001-2005.spe"
RA226_FIXTURE = FIXTURES / "M_ra_легкий_2001-2007.spe"

# F-397.1 — AtomSpectra fixture с embedded <BackgroundEnergySpectrum>.
ATOMSPECTRA_FIXTURES = (
    REPO_ROOT / "detectors" / "AtomSpectra" / "data" / "fixtures"
)
ATOM_CS137_WITH_EMBED_BG = ATOMSPECTRA_FIXTURES / "Cs137_0_см.xml"

# Default BG used by regen_demo_reports for archive fixtures (F-313).
DEFAULT_BG = (
    REPO_ROOT
    / "detectors"
    / "Gamma-1S"
    / "data"
    / "averaged_backgrounds"
    / "bg_2016_marinelli_water_marinelli.spe"
)


# ──────────────────────────────────────────────────────────────────
# Fast unit tests — no pipeline
# ──────────────────────────────────────────────────────────────────


class TestConfig:
    def test_default_config_has_all_sections(self):
        cfg = run_skill._load_config(None)
        assert "analyze" in cfg
        assert "multiplet" in cfg
        assert "v2" in cfg
        assert "artefacts" in cfg
        assert "output" in cfg

    def test_default_multiplet_factor_is_1_1(self):
        """F-387.1 NaI calibration default — owned by Agent A; the wrapper
        must mirror it, not override."""
        cfg = run_skill._load_config(None)
        assert cfg["multiplet"]["unresolved_separation_fwhm_factor"] == 1.1
        assert cfg["multiplet"]["max_components_per_cluster"] == 3

    def test_default_export_becqmoni_is_both(self):
        """F-RPT-04 / v1.18.29 — BecqMoni export default flipped 'both' → 'off'.

        Тест ловит регрессию default — если случайно вернётся 'both' (или
        что-то отличное от 'off'), это сломает оператора чьи fixtures
        генерируют невалидный BecqMoni XML. Имя теста сохранено для git
        blame trace; смысл инвертирован."""
        cfg = run_skill._load_config(None)
        assert cfg["analyze"]["export_becqmoni"] == "off"

    def test_deep_merge_preserves_unspecified_keys(self, tmp_path):
        user_cfg = {"analyze": {"sample_mass_kg": 0.5}}
        p = tmp_path / "cfg.json"
        p.write_text(json.dumps(user_cfg), encoding="utf-8")
        cfg = run_skill._load_config(p)
        assert cfg["analyze"]["sample_mass_kg"] == 0.5
        # Untouched keys preserved (F-RPT-04 default 'off' v1.18.29).
        assert cfg["analyze"]["export_becqmoni"] == "off"
        assert cfg["multiplet"]["unresolved_separation_fwhm_factor"] == 1.1

    def test_config_override_factor(self, tmp_path):
        user_cfg = {"multiplet": {"unresolved_separation_fwhm_factor": 1.0}}
        p = tmp_path / "cfg.json"
        p.write_text(json.dumps(user_cfg), encoding="utf-8")
        cfg = run_skill._load_config(p)
        assert cfg["multiplet"]["unresolved_separation_fwhm_factor"] == 1.0

    def test_missing_config_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            run_skill._load_config(tmp_path / "nope.json")


class TestMetadata:
    def test_extract_marinelli_geometry(self):
        meta = run_skill.SpectrumMetadata.from_path(
            Path("Th232_Маринелли_0cm.spe")
        )
        assert "маринелли" in (meta.geometry_hint or "").lower()
        assert meta.distance_cm == 0

    def test_extract_mass_kg_token(self):
        meta = run_skill.SpectrumMetadata.from_path(Path("sample_0.5kg.spe"))
        assert meta.sample_mass_kg == 0.5

    def test_extract_mass_g_token(self):
        meta = run_skill.SpectrumMetadata.from_path(Path("sample_500g.spe"))
        assert meta.sample_mass_kg == pytest.approx(0.5)

    def test_extract_background_token(self):
        meta = run_skill.SpectrumMetadata.from_path(Path("Фон_закр_кр_вода.spe"))
        assert meta.is_background_candidate is True

    def test_extract_detector_hint(self):
        meta = run_skill.SpectrumMetadata.from_path(Path("Gamma-1S_test.spe"))
        assert meta.detector_hint is not None
        assert "gamma" in meta.detector_hint

    def test_default_mass_marinelli(self):
        assert run_skill._default_mass_for_geometry("Маринелли") == 0.5

    def test_default_mass_unknown(self):
        assert run_skill._default_mass_for_geometry(None) == 1.0


class TestBundleLayout:
    def test_layout_paths_match_spec(self, tmp_path):
        cfg = run_skill._load_config(None)
        layout = run_skill.BundleLayout.from_base(tmp_path / "bundle", cfg)
        assert layout.sample.name == "sample"
        assert layout.sample_v2.name == "sample_v2"
        assert layout.compare.name == "v2_compare"
        assert layout.index_html.name == "index.html"
        assert layout.summary_json.name == "run_skill_summary.json"
        assert layout.log_file.name == "run_skill.log"

    def test_ensure_dirs_creates_phases(self, tmp_path):
        cfg = run_skill._load_config(None)
        layout = run_skill.BundleLayout.from_base(tmp_path / "b", cfg)
        layout.ensure_dirs()
        assert layout.base.exists()
        assert layout.phases_dir.exists()

    def test_marker_path(self, tmp_path):
        cfg = run_skill._load_config(None)
        layout = run_skill.BundleLayout.from_base(tmp_path / "b", cfg)
        assert layout.marker(3).name == "phase_3.done"


class TestAutoDetectBackground:
    def test_returns_none_when_no_siblings(self, tmp_path):
        sample = tmp_path / "Sample.spe"
        sample.write_text("", encoding="utf-8")
        import logging
        log = logging.getLogger("test")
        log.addHandler(logging.NullHandler())
        assert run_skill._auto_detect_background(sample, log) is None

    def test_finds_sibling_with_bg_token(self, tmp_path):
        sample = tmp_path / "Sample.spe"
        sample.write_text("", encoding="utf-8")
        bg = tmp_path / "Фон_закр.spe"
        bg.write_text("", encoding="utf-8")
        import logging
        log = logging.getLogger("test")
        log.addHandler(logging.NullHandler())
        result = run_skill._auto_detect_background(sample, log)
        assert result is not None
        assert result.name == "Фон_закр.spe"


class TestExitCodes:
    def test_missing_spectrum_returns_2(self, tmp_path):
        code = run_skill.main([
            str(tmp_path / "no_such.spe"),
            "--output-dir", str(tmp_path / "bundle"),
            "--quiet",
        ])
        assert code == 2

    def test_invalid_resume_path_returns_2(self, tmp_path):
        code = run_skill.main([
            "--resume", str(tmp_path / "no_bundle"),
            "--quiet",
        ])
        assert code == 2

    def test_resume_without_marker_returns_2(self, tmp_path):
        bundle = tmp_path / "bundle"
        bundle.mkdir()
        code = run_skill.main([
            "--resume", str(bundle),
            "--quiet",
        ])
        assert code == 2

    def test_no_args_returns_2(self):
        code = run_skill.main(["--quiet"])
        assert code == 2

    def test_version_flag(self, capsys):
        code = run_skill.main(["--version"])
        captured = capsys.readouterr()
        assert code == 0
        assert captured.out.strip()  # non-empty version string


class TestPhaseMarkers:
    def test_phase_not_done_initially(self, tmp_path):
        cfg = run_skill._load_config(None)
        layout = run_skill.BundleLayout.from_base(tmp_path / "b", cfg)
        layout.ensure_dirs()
        import logging
        log = logging.getLogger("test_markers")
        log.addHandler(logging.NullHandler())
        ctx = run_skill.RunContext(
            spectrum=tmp_path / "fake.spe",
            background=None,
            metadata=run_skill.SpectrumMetadata.from_path(Path("fake.spe")),
            cfg=cfg,
            layout=layout,
            logger=log,
            skill_version="test",
            include_v2=False,
        )
        assert not run_skill._phase_already_done(ctx, 2)

    def test_record_phase_creates_marker(self, tmp_path):
        cfg = run_skill._load_config(None)
        layout = run_skill.BundleLayout.from_base(tmp_path / "b", cfg)
        layout.ensure_dirs()
        import logging
        log = logging.getLogger("test_record")
        log.addHandler(logging.NullHandler())
        ctx = run_skill.RunContext(
            spectrum=tmp_path / "fake.spe",
            background=None,
            metadata=run_skill.SpectrumMetadata.from_path(Path("fake.spe")),
            cfg=cfg, layout=layout, logger=log,
            skill_version="test", include_v2=False,
        )
        res = run_skill.PhaseResult(
            phase=2, name="prod_analyze", status="ok", elapsed_s=1.0,
            detail={"x": 1},
        )
        run_skill._record_phase(ctx, 2, res)
        assert run_skill._phase_already_done(ctx, 2)
        loaded = run_skill._load_prior_phase(ctx, 2)
        assert loaded is not None
        assert loaded["detail"]["x"] == 1


class TestIndexRendering:
    def test_renders_html_with_stem(self, tmp_path):
        cfg = run_skill._load_config(None)
        layout = run_skill.BundleLayout.from_base(tmp_path / "b", cfg)
        layout.ensure_dirs()
        import logging
        log = logging.getLogger("test_index")
        log.addHandler(logging.NullHandler())
        ctx = run_skill.RunContext(
            spectrum=tmp_path / "Th232_Маринелли_0cm.spe",
            background=None,
            metadata=run_skill.SpectrumMetadata.from_path(
                Path("Th232_Маринелли_0cm.spe")
            ),
            cfg=cfg, layout=layout, logger=log,
            skill_version="test", include_v2=False,
        )
        html = run_skill._render_index_html(ctx)
        assert "<!doctype html" in html.lower()
        assert "Th232_Маринелли_0cm" in html
        assert "Полный отчёт образца (production)" in html  # F-RPT-01 cards label

    def test_v2_section_only_when_enabled(self, tmp_path):
        cfg = run_skill._load_config(None)
        layout = run_skill.BundleLayout.from_base(tmp_path / "b", cfg)
        layout.ensure_dirs()
        import logging
        log = logging.getLogger("test_index_v2")
        log.addHandler(logging.NullHandler())
        # include_v2 = False
        ctx_off = run_skill.RunContext(
            spectrum=tmp_path / "s.spe", background=None,
            metadata=run_skill.SpectrumMetadata.from_path(Path("s.spe")),
            cfg=cfg, layout=layout, logger=log,
            skill_version="test", include_v2=False,
        )
        html_off = run_skill._render_index_html(ctx_off)
        assert "V2-метод" not in html_off  # F-RPT-01: V2 card absent when include_v2=False

        # include_v2 = True
        ctx_on = run_skill.RunContext(
            spectrum=tmp_path / "s.spe", background=None,
            metadata=run_skill.SpectrumMetadata.from_path(Path("s.spe")),
            cfg=cfg, layout=layout, logger=log,
            skill_version="test", include_v2=True,
        )
        html_on = run_skill._render_index_html(ctx_on)
        assert "V2-метод" in html_on  # F-RPT-01: V2 card present when include_v2=True

    def test_no_pikvylet_terminology(self, tmp_path):
        """F-386: «пик вылета», not «ускользание»."""
        cfg = run_skill._load_config(None)
        layout = run_skill.BundleLayout.from_base(tmp_path / "b", cfg)
        layout.ensure_dirs()
        import logging
        log = logging.getLogger("test_terms")
        log.addHandler(logging.NullHandler())
        ctx = run_skill.RunContext(
            spectrum=tmp_path / "s.spe", background=None,
            metadata=run_skill.SpectrumMetadata.from_path(Path("s.spe")),
            cfg=cfg, layout=layout, logger=log,
            skill_version="test", include_v2=True,
        )
        html = run_skill._render_index_html(ctx)
        assert "ускользан" not in html.lower(), (
            "F-386: terminology lock — must not use «ускользание»"
        )


class TestF397_1_EmbeddedBgExtraction:
    """F-397.1: extract embedded bg from sample file → use as bg path."""

    def _make_fake_spec_with_embedded(self):
        """Минимальный Spectrum-like object для unit-теста extract'а."""
        from gamma.spectrum import Spectrum
        import numpy as np
        sample_counts = np.array([10, 20, 30, 25, 15], dtype=np.uint32)
        sample = Spectrum(
            counts=sample_counts,
            live_time=100.0,
            real_time=110.0,
            sample_id="TEST_SAMPLE",
            n_channels=len(sample_counts),
            n_channels_raw=len(sample_counts),
            energy_cal=(0.0, 1.0),
        )
        bg_counts = np.array([1, 2, 3, 2, 1], dtype=np.uint32)
        bg = Spectrum(
            counts=bg_counts,
            live_time=200.0,
            real_time=205.0,
            sample_id="TEST_BG_EMBEDDED",
            n_channels=len(bg_counts),
            n_channels_raw=len(bg_counts),
            energy_cal=(0.0, 1.0),
        )
        sample.background_embedded = bg
        return sample

    def test_extract_writes_file(self, tmp_path):
        """write_lsrm_spe + readback round-trips embedded bg."""
        import logging
        log = logging.getLogger("test_extract")
        log.addHandler(logging.NullHandler())
        cfg = run_skill._load_config(None)
        layout = run_skill.BundleLayout.from_base(tmp_path / "b", cfg)
        layout.ensure_dirs()
        spec = self._make_fake_spec_with_embedded()
        out = run_skill._extract_embedded_background(spec, layout, "Sample", log)
        assert out is not None
        assert out.exists()
        assert out.name == "Sample_embedded_bg.spe"
        # Round-trip — reader should accept it
        from gamma.io.readers import read_spectrum
        bg_read = read_spectrum(str(out))
        assert bg_read.live_time == 200.0

    def test_extract_returns_none_when_no_embedded(self, tmp_path):
        import logging
        log = logging.getLogger("test_extract_none")
        log.addHandler(logging.NullHandler())
        cfg = run_skill._load_config(None)
        layout = run_skill.BundleLayout.from_base(tmp_path / "b", cfg)
        layout.ensure_dirs()
        from gamma.spectrum import Spectrum
        import numpy as np
        c = np.array([1, 2, 3], dtype=np.uint32)
        spec = Spectrum(
            counts=c,
            live_time=10.0,
            real_time=11.0,
            sample_id="NO_EMBED",
            n_channels=len(c),
            n_channels_raw=len(c),
            energy_cal=(0.0, 1.0),
        )
        # No background_embedded set
        out = run_skill._extract_embedded_background(spec, layout, "Sample", log)
        assert out is None

    def test_extract_handles_None_spec(self, tmp_path):
        import logging
        log = logging.getLogger("test_extract_none_spec")
        log.addHandler(logging.NullHandler())
        cfg = run_skill._load_config(None)
        layout = run_skill.BundleLayout.from_base(tmp_path / "b", cfg)
        layout.ensure_dirs()
        out = run_skill._extract_embedded_background(None, layout, "X", log)
        assert out is None


class TestArtefactsCheck:
    def test_empty_dir_returns_missing(self, tmp_path):
        cfg = run_skill._load_config(None)
        d = tmp_path / "out"
        d.mkdir()
        missing = run_skill._check_required_artefacts(d, cfg)
        # When no artifacts present, all enabled categories report missing.
        assert "json" in missing
        assert "html" in missing
        assert "markdown" in missing

    def test_complete_dir_returns_empty(self, tmp_path):
        cfg = run_skill._load_config(None)
        # Disable XML/PDF/plots to keep the fixture minimal
        cfg["artefacts"]["xml_bq"] = False
        cfg["artefacts"]["technical_pdf"] = False
        cfg["artefacts"]["plots"] = False
        d = tmp_path / "out"
        d.mkdir()
        (d / "Sample_report.json").write_text("{}", encoding="utf-8")
        (d / "Sample_report.md").write_text("#", encoding="utf-8")
        (d / "Sample_report.html").write_text("<html>", encoding="utf-8")
        assert run_skill._check_required_artefacts(d, cfg) == []


# ──────────────────────────────────────────────────────────────────
# Slow integration tests (full pipeline)
# ──────────────────────────────────────────────────────────────────

# Pick the lightest reliable fixture: Cs-137 archive (single FEP, fast).
INTEGRATION_FIXTURE = CS137_FIXTURE


def _have_fixture() -> bool:
    return INTEGRATION_FIXTURE.exists() and DEFAULT_BG.exists()


@pytest.fixture
def smoke_bundle(tmp_path: Path) -> Path:
    """Per-test bundle directory under tmp."""
    return tmp_path / "bundle"


@pytest.mark.slow
@pytest.mark.skipif(not _have_fixture(), reason="archive fixture missing")
def test_smoke_full_pipeline(smoke_bundle):
    """End-to-end на Cs-137. Без V2 — production-only, экономия времени."""
    code = run_skill.main([
        str(INTEGRATION_FIXTURE),
        "--background", str(DEFAULT_BG),
        "--mass", "0.570",
        "--output-dir", str(smoke_bundle),
        "--quiet",
    ])
    assert code == 0, f"run_skill exit={code}"
    # Standard artefacts present?
    sample = smoke_bundle / "sample"
    assert sample.exists()
    assert any(sample.glob("*_report.json"))
    assert any(sample.glob("*_report.html"))
    assert any(sample.glob("*_report.md"))
    # Index + summary + log
    # F-243-wave3 fix: index.html is only generated when cfg.reports.bundle_index=true
    # (compare-mode). Default single-variant runs skip phase 7 (bundle_index=False gate
    # added in wave 2/3 per user request 2026-06-03). index.html absent is correct.
    assert (smoke_bundle / "run_skill_summary.json").exists()
    assert (smoke_bundle / "run_skill.log").exists()
    # Phase markers — at least 0, 1, 2, 3, 8 must exist.
    # F-243-wave3 fix: phase 7 (bundle_index) skipped by default → no phase_7.done marker.
    for phase in (0, 1, 2, 3, 8):
        assert (smoke_bundle / ".phases" / f"phase_{phase}.done").exists(), (
            f"marker phase_{phase}.done missing"
        )


@pytest.mark.slow
@pytest.mark.skipif(not _have_fixture(), reason="archive fixture missing")
def test_resume_skips_completed_phases(smoke_bundle):
    """После первого прогона повторный --resume должен пропустить phases."""
    # First run — full
    code1 = run_skill.main([
        str(INTEGRATION_FIXTURE),
        "--background", str(DEFAULT_BG),
        "--mass", "0.570",
        "--output-dir", str(smoke_bundle),
        "--quiet",
    ])
    assert code1 == 0

    # Capture summary.json before second run
    summary_path = smoke_bundle / "run_skill_summary.json"
    summary1 = json.loads(summary_path.read_text(encoding="utf-8"))

    # Second run with --resume — should be fast (resumed status everywhere)
    code2 = run_skill.main([
        "--resume", str(smoke_bundle),
        "--quiet",
    ])
    assert code2 == 0
    summary2 = json.loads(summary_path.read_text(encoding="utf-8"))

    # In resumed run, phases 0-3 must be marked "resumed".
    # F-243-wave3 fix: phase 7 (bundle_index) is gated by cfg.reports.bundle_index=False
    # by default (wave 2/3 change). Since it returns "skipped" on first run, no
    # phase_7.done marker is written, so on resume it returns "skipped" again — not
    # "resumed". Phase 7 resume is only exercised in compare-mode (bundle_index=true).
    for phase in (0, 1, 2, 3):
        status = summary2["phases"][str(phase)]["status"]
        assert status == "resumed", (
            f"phase {phase} expected 'resumed', got {status}"
        )
    # Phase 7 should be "skipped" (gated off) in default single-variant runs.
    status_7 = summary2["phases"]["7"]["status"]
    assert status_7 == "skipped", (
        f"phase 7 expected 'skipped' (bundle_index disabled), got {status_7}"
    )


@pytest.mark.slow
@pytest.mark.skipif(not _have_fixture(), reason="archive fixture missing")
def test_idempotency_json_stable(smoke_bundle, tmp_path):
    """Два независимых прогона на тот же спектр → одинаковый JSON SHA256.

    Note: technical_pdf может варьироваться (метаданные reportlab), потому
    проверяем только _report.json. В нём всё что важно для пайплайна.
    """
    bundle_a = tmp_path / "a"
    bundle_b = tmp_path / "b"
    args_common = [
        str(INTEGRATION_FIXTURE),
        "--background", str(DEFAULT_BG),
        "--mass", "0.570",
        "--no-pdf", "--no-plots", "--no-xml",  # speed up, focus on JSON
        "--quiet",
    ]
    assert run_skill.main(args_common + ["--output-dir", str(bundle_a)]) == 0
    assert run_skill.main(args_common + ["--output-dir", str(bundle_b)]) == 0

    json_a = next((bundle_a / "sample").glob("*_report.json"))
    json_b = next((bundle_b / "sample").glob("*_report.json"))

    # Identical bytes (or at minimum identical content modulo whitespace).
    h_a = _sha256_normalized_json(json_a)
    h_b = _sha256_normalized_json(json_b)
    assert h_a == h_b, (
        f"JSON differs between runs:\n  {json_a}\n  {json_b}"
    )


def _sha256_normalized_json(path: Path) -> str:
    """SHA256 после canonical serialization — стабилен к key-order и whitespace."""
    data = json.loads(path.read_text(encoding="utf-8"))
    # Drop timestamps / paths that legitimately vary between runs.
    data = _strip_volatile_keys(data)
    canon = json.dumps(data, sort_keys=True, ensure_ascii=False)
    return hashlib.sha256(canon.encode("utf-8")).hexdigest()


def _strip_volatile_keys(obj):
    """Удалить ключи, которые меняются между прогонами (даты, пути)."""
    VOLATILE = {
        "generated_at", "generated_at_unix", "report_generated_at",
        "input_path", "absolute_path", "output_dir", "cost_estimate",
    }
    if isinstance(obj, dict):
        return {k: _strip_volatile_keys(v) for k, v in obj.items()
                if k not in VOLATILE}
    if isinstance(obj, list):
        return [_strip_volatile_keys(x) for x in obj]
    return obj


@pytest.mark.slow
@pytest.mark.skipif(
    not ATOM_CS137_WITH_EMBED_BG.exists(),
    reason="AtomSpectra Cs137_0_см.xml fixture missing",
)
def test_f397_1_embedded_bg_auto_extracted(tmp_path: Path):
    """F-397.1: запуск на AtomSpectra-файле с embedded <BackgroundEnergySpectrum>
    БЕЗ --background → run_skill сам извлекает embedded bg и пайплайн
    обрабатывает его штатно (background_source=embedded_extracted в summary)."""
    bundle = tmp_path / "atom_bundle"
    code = run_skill.main([
        str(ATOM_CS137_WITH_EMBED_BG),
        "--output-dir", str(bundle),
        "--no-pdf", "--no-plots", "--no-xml",  # speed
        "--quiet",
    ])
    assert code == 0, f"run_skill exit={code}"
    # Embedded bg file должен быть создан в .embedded_bg/
    embed_dir = bundle / ".embedded_bg"
    assert embed_dir.exists(), "missing .embedded_bg/ marker dir"
    extracted = list(embed_dir.glob("*_embedded_bg.spe"))
    assert len(extracted) == 1, (
        f"expected exactly one extracted bg, got {[p.name for p in extracted]}"
    )
    # Summary должен зафиксировать source
    summary = json.loads(
        (bundle / "run_skill_summary.json").read_text(encoding="utf-8")
    )
    p1_detail = summary["phases"]["1"]["detail"]
    assert p1_detail["background_embedded_present"] is True
    assert p1_detail["background_source"] == "embedded_extracted"
    assert "background_extracted_path" in p1_detail


@pytest.mark.slow
@pytest.mark.skipif(not _have_fixture(), reason="archive fixture missing")
def test_batch_three_fixtures(tmp_path):
    """3 .spe → 3 bundles + manifest.csv."""
    fixtures_dir = tmp_path / "inputs"
    fixtures_dir.mkdir()
    # Copy 3 fixtures into a flat dir for the glob.
    import shutil
    candidates = [CS137_FIXTURE, K40_FIXTURE, TH232_FIXTURE]
    available = [f for f in candidates if f.exists()]
    if len(available) < 3:
        pytest.skip(f"only {len(available)} archive fixtures available")
    for f in available[:3]:
        shutil.copy(f, fixtures_dir / f.name)

    out_root = tmp_path / "batch_out"
    code = run_skill.main([
        "--batch", str(fixtures_dir),
        "--output-dir", str(out_root),
        "--no-pdf", "--no-plots", "--no-xml",  # speed
        "--quiet",
    ])
    # Exit code can be 0 or 1 (some bundles may have non-fatal partial fails);
    # we just need the manifest + 3 bundles.
    assert code in (0, 1)
    assert (out_root / "manifest.csv").exists()
    bundles = [p for p in out_root.iterdir() if p.is_dir()]
    assert len(bundles) == 3
    for b in bundles:
        assert (b / "run_skill_summary.json").exists()
