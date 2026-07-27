# -*- coding: utf-8 -*-
"""F-388 / v1.18.26 — conditional U-238 chain-rejection block.

Контракт: блок «Почему … цепочка не определяется» (interactive HTML
notes block + markdown report express section) рендерится только
если у нуклидов соответствующей подавлённой цепочки есть хотя бы один
`primary_fep` в спектре. Иначе шум для пользователя.

Раньше: блок всегда рендерился из cd.suppressed_chains.
Сейчас: фильтрация по `_chain_has_evidence`.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))


# ──────────────────────────────────────────────────────────────────
# Helpers — minimal report fixture
# ──────────────────────────────────────────────────────────────────

def _make_report(primary_feps, suppressed_chains):
    """Минимальный JSON-отчёт shape для _build_notes_blocks / markdown."""
    return {
        "primary_feps": primary_feps,
        "secondary_peaks": [],
        "warnings": [],
        "pipeline_notes": [],
        "priority_express_findings": [],
        "diagnostics": {
            "chain_dominance": {
                "th232_dominant": False,
                "u238_dominant": False,
                "suppressed_chains": suppressed_chains,
                "suppression_reason": (
                    "Filename binds source to Th-232 chain only. "
                    "U-238 chain suppression."
                ),
                "chain_filtered_out_nuclides": [],
            },
        },
    }


# ──────────────────────────────────────────────────────────────────
# interactive_html — _build_notes_blocks
# ──────────────────────────────────────────────────────────────────

def test_f388_html_skips_block_when_no_u238_fep():
    """U-238 suppressed, но НЕТ ни одного nuclide из U-цепочки в FEP'ах →
    блок отсутствует."""
    from gamma.reporting.interactive_html import _build_notes_blocks

    report = _make_report(
        primary_feps=[
            {"nuclide": "Tl-208", "peak_E_keV": 2614.5,
             "library_E_keV": 2614.5, "peak_area_counts": 50000.0},
            {"nuclide": "Ac-228", "peak_E_keV": 911.2,
             "library_E_keV": 911.2, "peak_area_counts": 12000.0},
        ],
        suppressed_chains=["U-238"],
    )
    html = _build_notes_blocks(report, analysis_result=None)
    assert "не определяется" not in html, (
        "Блок не должен рендериться: U-238 не присутствует в primary_feps"
    )


def test_f388_html_renders_block_when_u238_fep_present():
    """Bi-214 примем как доказательство U-238 присутствия → блок рендерится."""
    from gamma.reporting.interactive_html import _build_notes_blocks

    report = _make_report(
        primary_feps=[
            {"nuclide": "Bi-214", "peak_E_keV": 609.3,
             "library_E_keV": 609.3, "peak_area_counts": 5000.0},
            {"nuclide": "Tl-208", "peak_E_keV": 2614.5,
             "library_E_keV": 2614.5, "peak_area_counts": 50000.0},
        ],
        suppressed_chains=["U-238"],
    )
    html = _build_notes_blocks(report, analysis_result=None)
    assert "не определяется" in html, (
        "Блок должен рендериться: Bi-214 — член U-238 цепочки"
    )
    assert "U-238" in html


def test_f388_html_keeps_th232_block_when_th_fep_present():
    """Симметрия: Th-232 suppressed + Ac-228 в FEP'ах → блок рендерится."""
    from gamma.reporting.interactive_html import _build_notes_blocks

    report = _make_report(
        primary_feps=[
            {"nuclide": "Ac-228", "peak_E_keV": 911.2,
             "library_E_keV": 911.2, "peak_area_counts": 12000.0},
        ],
        suppressed_chains=["Th-232"],
    )
    html = _build_notes_blocks(report, analysis_result=None)
    assert "не определяется" in html
    assert "Th-232" in html


def test_f388_html_partial_filter_when_some_chains_relevant():
    """Suppressed = [U-238, Th-232], но в FEP'ах только Bi-214 (U) →
    рендерим блок только для U-238."""
    from gamma.reporting.interactive_html import _build_notes_blocks

    report = _make_report(
        primary_feps=[
            {"nuclide": "Bi-214", "peak_E_keV": 609.3,
             "library_E_keV": 609.3, "peak_area_counts": 5000.0},
        ],
        suppressed_chains=["U-238", "Th-232"],
    )
    html = _build_notes_blocks(report, analysis_result=None)
    assert "не определяется" in html
    # Th-232 не должен попасть в финальный заголовок
    # (Bi-214 — только U-член)
    # Block heading формат: "Почему U-238 цепочка не определяется"
    # Поэтому проверим по тексту что Th-232 в заголовке отсутствует
    # (но может быть в других местах HTML).
    # Локализуем check внутрь самого блока через раздел <h3>...не определяется
    import re
    match = re.search(
        r"<h3>Почему (.+?) цепочка не определяется</h3>", html
    )
    assert match is not None, "Заголовок блока не найден"
    chains_in_title = match.group(1)
    assert "U-238" in chains_in_title
    assert "Th-232" not in chains_in_title


# ──────────────────────────────────────────────────────────────────
# markdown_report — _render_priority_express
# ──────────────────────────────────────────────────────────────────

def test_f388_markdown_skips_block_when_no_u238_fep():
    """Markdown симметрично: U-238 suppressed, но НЕТ U-членов в FEP'ах →
    блок «Подавление цепочки …» отсутствует."""
    from gamma.reporting.markdown_report import _render_priority_express

    report = _make_report(
        primary_feps=[
            {"nuclide": "Tl-208", "peak_E_keV": 2614.5,
             "library_E_keV": 2614.5, "peak_area_counts": 50000.0},
        ],
        suppressed_chains=["U-238"],
    )
    # markdown_report использует priority_express_findings → добавим один
    report["priority_express_findings"] = [
        {"order": 1, "label": "Tl-208 2614", "matched": True,
         "max_significance_sigma": 10.0, "note": "MATCHED"},
    ]
    md = _render_priority_express(report)
    assert "Подавление цепочки" not in md, (
        "Блок не должен рендериться: U-238 не присутствует в primary_feps"
    )


def test_f388_markdown_renders_block_when_u238_fep_present():
    """Bi-214 в FEP'ах → markdown блок рендерится."""
    from gamma.reporting.markdown_report import _render_priority_express

    report = _make_report(
        primary_feps=[
            {"nuclide": "Bi-214", "peak_E_keV": 609.3,
             "library_E_keV": 609.3, "peak_area_counts": 5000.0},
        ],
        suppressed_chains=["U-238"],
    )
    report["priority_express_findings"] = [
        {"order": 1, "label": "Bi-214 609", "matched": True,
         "max_significance_sigma": 4.0, "note": "MATCHED"},
    ]
    md = _render_priority_express(report)
    assert "Подавление цепочки" in md, (
        "Блок должен рендериться: Bi-214 — член U-238 цепочки"
    )


# ──────────────────────────────────────────────────────────────────
# Regression: пустой suppressed → блок отсутствует (как до F-388)
# ──────────────────────────────────────────────────────────────────

def test_f388_empty_suppressed_no_block_html():
    """suppressed_chains=[] → блок отсутствует (это было и до F-388)."""
    from gamma.reporting.interactive_html import _build_notes_blocks

    report = _make_report(
        primary_feps=[
            {"nuclide": "Bi-214", "peak_E_keV": 609.3,
             "library_E_keV": 609.3, "peak_area_counts": 5000.0},
        ],
        suppressed_chains=[],
    )
    html = _build_notes_blocks(report, analysis_result=None)
    assert "не определяется" not in html


def test_f388_empty_suppressed_no_block_markdown():
    """Markdown: suppressed_chains=[] → блок отсутствует."""
    from gamma.reporting.markdown_report import _render_priority_express

    report = _make_report(
        primary_feps=[
            {"nuclide": "Bi-214", "peak_E_keV": 609.3,
             "library_E_keV": 609.3, "peak_area_counts": 5000.0},
        ],
        suppressed_chains=[],
    )
    report["priority_express_findings"] = [
        {"order": 1, "label": "Bi-214 609", "matched": True,
         "max_significance_sigma": 4.0, "note": "MATCHED"},
    ]
    md = _render_priority_express(report)
    assert "Подавление цепочки" not in md
