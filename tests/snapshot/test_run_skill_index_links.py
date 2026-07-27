# -*- coding: utf-8 -*-
"""BUG-10 / 2026-06-02 — `<a href="…">` integrity in run_skill index.html.

Regression: `_render_index_html` joined `ctx.metadata.stem + "_report.html"`
naïvely, but the writer (`gamma.reporting.build._safe_filename_stem`) scrubs
embedded S/N tokens (F-115). For files like `Th232_420-7-17_Маринелли_0cm
.spe` the bundle directory keeps the raw name, while the report file is
written as `Th232_Маринелли_0cm_report.html`. The cards-style `index.html`
therefore pointed at non-existent paths and every link 404'd.

These tests exercise `_render_index_html` end-to-end and parse the rendered
HTML to ensure every `<a class="card" href="…">` resolves to a file that
actually exists on disk inside the bundle.
"""
from __future__ import annotations

import importlib
import logging
import re
import sys
from pathlib import Path
from typing import Any, Dict, List

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

run_skill = importlib.import_module("run_skill")


_HREF_RE = re.compile(r'<a\s+class="card"\s+href="([^"]+)"', re.IGNORECASE)


def _hrefs(html: str) -> List[str]:
    """Extract every active card href (skipping `<div class="card dim">`)."""
    return _HREF_RE.findall(html)


def _make_ctx(
    *,
    bundle: Path,
    spectrum_name: str,
    include_v2: bool,
    phase2_html: Path | None = None,
    phase4_html: Path | None = None,
) -> Any:
    """Build a minimal `RunContext` wired to render an index page.

    `phase2_html` / `phase4_html`, if provided, are wired into
    `ctx.phase_results[N]["detail"]["html"]` so the resolver picks up
    artefact paths the way it would after a real or resumed pipeline run.
    """
    cfg = run_skill._load_config(None)
    layout = run_skill.BundleLayout.from_base(bundle, cfg)
    layout.ensure_dirs()
    layout.sample.mkdir(parents=True, exist_ok=True)
    layout.sample_v2.mkdir(parents=True, exist_ok=True)
    layout.compare.mkdir(parents=True, exist_ok=True)

    log = logging.getLogger(f"test_index_links_{bundle.name}")
    log.addHandler(logging.NullHandler())

    ctx = run_skill.RunContext(
        spectrum=bundle.parent / spectrum_name,
        background=None,
        metadata=run_skill.SpectrumMetadata.from_path(Path(spectrum_name)),
        cfg=cfg,
        layout=layout,
        logger=log,
        skill_version="test",
        include_v2=include_v2,
    )
    # Wire phase artefact details (mirrors `_record_phase` after
    # analyze_and_report returns).
    if phase2_html is not None:
        ctx.phase_results[2] = {"detail": {"html": str(phase2_html)}}
    if phase4_html is not None:
        ctx.phase_results[4] = {"detail": {"html": str(phase4_html)}}
    return ctx


# ──────────────────────────────────────────────────────────────────
# Core regression: BUG-10 — S/N-bearing filename → scrubbed report stem
# ──────────────────────────────────────────────────────────────────


class TestIndexLinkIntegrity:
    """Each `<a class="card" href="…">` in the rendered index must resolve
    to an existing file under `ctx.layout.base`."""

    def test_sn_stripped_stem_resolves(self, tmp_path: Path) -> None:
        """Reproducer for BUG-10.

        Bundle directory carries an S/N token (`420-7-17`) but the report
        writer scrubs it. Pre-fix, the index linked at
        `sample/Th232_420-7-17_Маринелли_0cm_report.html` (404).
        """
        bundle = tmp_path / "Th232_420-7-17_Маринелли_0cm"
        spectrum_name = "Th232_420-7-17_Маринелли_0cm.spe"

        # Lay down the actual files the writer would emit (S/N scrubbed).
        sample_report = bundle / "sample" / "Th232_Маринелли_0cm_report.html"
        sample_v2_report = (
            bundle / "sample_v2" / "Th232_Маринелли_0cm_report.html"
        )
        compare_report = bundle / "v2_compare" / "v2_compare_report.html"

        ctx = _make_ctx(
            bundle=bundle,
            spectrum_name=spectrum_name,
            include_v2=True,
            phase2_html=sample_report,
            phase4_html=sample_v2_report,
        )
        # Materialize the artefacts so `_exists` and link-target assertions
        # both pass.
        for f in (sample_report, sample_v2_report, compare_report):
            f.parent.mkdir(parents=True, exist_ok=True)
            f.write_text("<html><body>x</body></html>", encoding="utf-8")

        html = run_skill._render_index_html(ctx)
        hrefs = _hrefs(html)

        # Three active cards expected.
        assert len(hrefs) == 3, f"expected 3 card links, got {len(hrefs)}: {hrefs}"

        # None of them should still embed the raw S/N stem.
        for href in hrefs:
            assert "420-7-17" not in href, (
                f"BUG-10 regression: href still carries S/N token: {href!r}"
            )

        # Every href must resolve to an existing file.
        for href in hrefs:
            target = bundle / href
            assert target.exists(), (
                f"broken link: <a href={href!r}> → {target} does not exist"
            )

    def test_glob_fallback_when_phase_detail_missing(
        self, tmp_path: Path
    ) -> None:
        """`ctx.phase_results` empty (e.g. interrupted run, re-rendered
        from disk). Resolver must fall back to a `*_report.html` glob and
        still produce live links."""
        bundle = tmp_path / "Sample_Bundle_500g"
        spectrum_name = "Sample_500g.spe"

        sample_report = bundle / "sample" / "Sample_500g_report.html"
        sample_v2_report = bundle / "sample_v2" / "Sample_500g_report.html"
        compare_report = bundle / "v2_compare" / "v2_compare_report.html"
        for f in (sample_report, sample_v2_report, compare_report):
            f.parent.mkdir(parents=True, exist_ok=True)
            f.write_text("<html></html>", encoding="utf-8")

        # NB: no phase2_html / phase4_html — purely glob-driven.
        ctx = _make_ctx(
            bundle=bundle,
            spectrum_name=spectrum_name,
            include_v2=True,
        )

        html = run_skill._render_index_html(ctx)
        hrefs = _hrefs(html)
        assert len(hrefs) == 3
        for href in hrefs:
            assert (bundle / href).exists(), (
                f"glob-fallback returned non-existent target: {href}"
            )

    def test_v2_off_yields_single_card(self, tmp_path: Path) -> None:
        """With `include_v2=False` only the production card is rendered;
        that link must still resolve."""
        bundle = tmp_path / "Cs137_only"
        spectrum_name = "Cs137_2024.spe"
        sample_report = bundle / "sample" / "Cs137_2024_report.html"
        sample_report.parent.mkdir(parents=True, exist_ok=True)
        sample_report.write_text("<html></html>", encoding="utf-8")

        ctx = _make_ctx(
            bundle=bundle,
            spectrum_name=spectrum_name,
            include_v2=False,
            phase2_html=sample_report,
        )
        html = run_skill._render_index_html(ctx)
        hrefs = _hrefs(html)
        assert len(hrefs) == 1
        assert (bundle / hrefs[0]).exists()

    def test_demo_bundle_links_resolve(self) -> None:
        """Acceptance smoke against the committed Th-232 demo bundle.

        Skips when the bundle is absent (clean checkout / CI without
        demo regen). When present, every active `<a class="card" href>`
        must point at an existing file — protecting against any future
        regression in `_render_index_html` link generation.
        """
        demo_root = (
            REPO_ROOT
            / "demo_reports"
            / "Th232_420-7-17_Маринелли_0cm"
        )
        index_html = demo_root / "index.html"
        if not index_html.exists():
            pytest.skip(f"demo bundle absent: {index_html}")

        html = index_html.read_text(encoding="utf-8")
        hrefs = _hrefs(html)
        assert hrefs, "no active card links found in demo index.html"
        for href in hrefs:
            target = demo_root / href
            assert target.exists(), (
                f"demo index link broken: <a href={href!r}> → "
                f"{target} missing"
            )
