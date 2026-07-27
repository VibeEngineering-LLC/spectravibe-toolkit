# -*- coding: utf-8 -*-
"""F-397.2 / v1.18.28 (Agent B) — bg multiplet canvas rendering.

Контракт: при наличии bg-мультиплетов отрисованный HTML обязан содержать
УНИКАЛЬНЫЕ canvas-IDs для sample-блока (`mp-<X>`) и bg-блока (`mp-<X>-bg`).
Без уникальности два HTML-id collide → нарушение спецификации, плюс
`document.getElementById` возвращает только первый канвас, что оставляет
bg-multiplet charts всегда пустыми.

Также проверяем семантическое отличие H2 заголовков: bg-блок имеет
суффикс "(фон)" для безJS-доступности (print, reader-mode).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))


def _make_fake_multiplet(mid: str = "M1c", title: str = "Test"):
    return {
        "id": mid,
        "title": title,
        "chi2_per_dof": 1.0,
        "n_channels": 50,
        "closure_pct": 0.0,
        "phase_A_chi2_per_dof": None,
        "components": [{"nuclide": "Cs-137", "E_keV": 661.66,
                        "I_pct": 85.0, "area": 1000.0}],
    }


class TestF397_2_BgMultipletBlocks:
    """Unit tests for _build_multiplet_blocks(bg=True/False)."""

    def test_sample_block_uses_default_id_namespace(self):
        from gamma.reporting.interactive_html import _build_multiplet_blocks
        out = _build_multiplet_blocks([_make_fake_multiplet("M1c")], bg=False)
        assert 'id="mp-M1c"' in out
        assert 'id="mp-M1c-bg"' not in out

    def test_bg_block_uses_bg_suffix(self):
        from gamma.reporting.interactive_html import _build_multiplet_blocks
        out = _build_multiplet_blocks([_make_fake_multiplet("M1c")], bg=True)
        assert 'id="mp-M1c-bg"' in out
        # The "bg" suffix prevents id-collision with the sample block.

    def test_bg_block_header_has_phon_suffix(self):
        """BUG-5 / v1.18.30+ (Agent B): bg-блок имеет явный label «в фоновом
        спектре» вместо краткого «(фон)» — отличает от sample по роли."""
        from gamma.reporting.interactive_html import _build_multiplet_blocks
        out = _build_multiplet_blocks([_make_fake_multiplet()], bg=True)
        assert "в фоновом спектре" in out
        assert "референс для вычитания" in out

    def test_sample_block_header_has_no_phon_suffix(self):
        """BUG-5: sample-блок помечен «в спектре образца», без «фон»/«bg»."""
        from gamma.reporting.interactive_html import _build_multiplet_blocks
        out = _build_multiplet_blocks([_make_fake_multiplet()], bg=False)
        assert "в спектре образца" in out
        assert "первичная подгонка" in out
        assert "фоновом" not in out

    def test_empty_bg_block_has_phon_suffix_too(self):
        """BUG-5: empty-state bg тоже имеет label «в фоновом спектре»."""
        from gamma.reporting.interactive_html import _build_multiplet_blocks
        out = _build_multiplet_blocks([], bg=True)
        assert "в фоновом спектре" in out

    def test_empty_sample_block_no_phon_suffix(self):
        """BUG-5: empty-state sample → label «в спектре образца», без «фон»."""
        from gamma.reporting.interactive_html import _build_multiplet_blocks
        out = _build_multiplet_blocks([], bg=False)
        assert "в спектре образца" in out
        assert "фоновом" not in out

    def test_canvas_ids_unique_across_sample_and_bg(self):
        """The core anti-regression: combined sample+bg HTML must
        contain no duplicate `id="mp-..."` value."""
        from gamma.reporting.interactive_html import _build_multiplet_blocks
        ms = [_make_fake_multiplet("M1c"), _make_fake_multiplet("M2c")]
        sample = _build_multiplet_blocks(ms, bg=False)
        bg = _build_multiplet_blocks(ms, bg=True)
        combined = sample + bg
        ids = re.findall(r'id="(mp-[^"]+)"', combined)
        assert len(ids) == 4, f"expected 4 canvas ids, got {ids}"
        assert len(ids) == len(set(ids)), (
            f"duplicate canvas ids: {ids} (set={set(ids)})"
        )


# ──────────────────────────────────────────────────────────────────
# Integration test — slow, full pipeline + HTML render check
# ──────────────────────────────────────────────────────────────────

ARCHIVE_TH = REPO_ROOT / "evals" / "fixtures" / "M_th_легкий_2001-2005.spe"
DEFAULT_BG = (
    REPO_ROOT / "detectors" / "Gamma-1S" / "data"
    / "averaged_backgrounds" / "bg_2016_marinelli_water_marinelli.spe"
)


@pytest.mark.slow
@pytest.mark.skipif(
    not (ARCHIVE_TH.exists() and DEFAULT_BG.exists()),
    reason="archive Th-232 fixture or default bg missing",
)
def test_rendered_html_has_unique_canvas_ids(tmp_path):
    """End-to-end: реальный pipeline + reporter → HTML с уникальными canvas IDs."""
    import importlib
    run_skill = importlib.import_module("run_skill")
    bundle = tmp_path / "th_bundle"
    code = run_skill.main([
        str(ARCHIVE_TH),
        "--background", str(DEFAULT_BG),
        "--mass", "0.550",
        "--output-dir", str(bundle),
        "--no-pdf", "--no-plots", "--no-xml",  # speed
        "--quiet",
    ])
    assert code == 0
    html_files = list((bundle / "sample").glob("*_report.html"))
    assert len(html_files) == 1
    html = html_files[0].read_text(encoding="utf-8")
    # Извлекаем все mp-... id и проверяем их уникальность.
    ids = re.findall(r'id="(mp-[^"]+)"', html)
    # Archive Th-232 в эту fixture обычно не имеет bg-multiplets, но
    # отсутствие duplicates обязано выполняться всегда.
    assert len(ids) == len(set(ids)), (
        f"duplicate canvas ids in rendered HTML: "
        f"{[i for i in ids if ids.count(i) > 1]}"
    )
