# -*- coding: utf-8 -*-
"""BUG-16 — ГОСТ 26874-86 должен присутствовать в bibliography subset
демо-отчётов (Th-232 и др.).

User-reported symptom (screenshot): bibliography panel Th-232 demo report
НЕ содержит ГОСТ 26874-86 «Спектрометры энергий ионизирующих излучений.
Методы измерения основных параметров», хотя PDF присутствует в репо
(``books_library/01_methodology_pdf/GOST_26874-86_спектрометры_методы_измерения.pdf``)
и затекстовая запись есть в ``references/REFERENCES.md`` §1 запись № 24.

Корневая причина: ``scripts/gamma/reporting/build.py``
``_load_references_map()`` не содержал запись 24, а
``baseline_refs`` в ``_f318_append_gost_references()`` ограничивался
``{2, 7, 12, 19}`` — следовательно даже при автоматическом text-scan
ГОСТ 26874-86 не появлялся ни в каком отчёте.

Контракт после BUG-16:
  • ``_load_references_map()`` СОДЕРЖИТ запись 24 с точным названием
    титульного листа PDF.
  • ``baseline_refs`` СОДЕРЖИТ 24 (foundational normative document
    для γ-спектрометрии — должен быть в любом regulated отчёте).
  • PDF файл существует по канонической дорожке
    ``books_library/01_methodology_pdf/``.

Если будущая правка нарушит контракт — test fails.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from gamma.reporting.build import (  # noqa: E402
    _f318_append_gost_references,
    _load_references_map,
)


def test_bug16_entry_24_in_references_map():
    """``_load_references_map()`` должен содержать запись 24
    с упоминанием ГОСТ 26874 и точного титульного названия."""
    refs = _load_references_map()
    assert 24 in refs, (
        "BUG-16: запись 24 (ГОСТ 26874-86) отсутствует в "
        "_load_references_map(). Bibliography subset never renders her."
    )
    citation = refs[24]
    # Точное название титульного листа PDF (страница 1–2):
    # «Спектрометры энергий ионизирующих излучений.
    #  Методы измерения основных параметров»
    assert "26874" in citation, (
        f"BUG-16: citation для № 24 не содержит '26874'. Got: {citation!r}"
    )
    assert "Спектрометры энергий ионизирующих излучений" in citation, (
        "BUG-16: citation для № 24 не содержит точное название "
        "титульного листа ГОСТ 26874-86. Got: " + repr(citation)
    )
    assert "Методы измерения основных параметров" in citation, (
        "BUG-16: citation для № 24 не содержит подзаголовок "
        "«Методы измерения основных параметров». Got: " + repr(citation)
    )


def test_bug16_pdf_exists_in_books_library():
    """PDF файл должен существовать в books_library/01_methodology_pdf/.
    Antihallucination guard: ссылка в bibliography subset должна
    указывать на реально существующий артефакт в репозитории.
    """
    pdf_dir = REPO / "books_library" / "01_methodology_pdf"
    assert pdf_dir.is_dir(), (
        f"BUG-16: directory {pdf_dir} не существует — "
        "bibliography subset ссылается на отсутствующий артефакт."
    )
    pdf_path = pdf_dir / "GOST_26874-86_спектрометры_методы_измерения.pdf"
    assert pdf_path.exists(), (
        f"BUG-16: PDF файл {pdf_path.name} не найден в {pdf_dir}. "
        "Bibliography subset ссылается на отсутствующий артефакт."
    )
    # Sanity-size check: реальный ГОСТ — около 1.8 MB, не пустой stub.
    size = pdf_path.stat().st_size
    assert size > 100_000, (
        f"BUG-16: PDF файл {pdf_path.name} слишком мал "
        f"({size} bytes), вероятно corrupted / placeholder."
    )


def test_bug16_entry_24_present_in_md_baseline_subset():
    """``_f318_append_gost_references(text, format='md')`` без явных
    ссылок [N, локатор] в тексте должен всё равно включить запись 24
    в Список использованной литературы (через baseline_refs subset)."""
    # Минимальный input — markdown без явных Layer-2 цитат.
    text = "# Th-232 demo report\n\nResult body without citations.\n"
    out = _f318_append_gost_references(text, format="md")
    assert "Список использованной литературы" in out, (
        "BUG-16: _f318_append_gost_references не вставил appendix-секцию."
    )
    # Поиск строки "24. <citation>" в конце.
    import re
    m = re.search(r"^24\.\s+(.+)$", out, re.MULTILINE)
    assert m is not None, (
        "BUG-16: запись 24 (ГОСТ 26874-86) не появилась в MD bibliography "
        "subset, хотя должна быть включена через baseline_refs."
    )
    line = m.group(1)
    assert "26874" in line and "Спектрометры энергий" in line, (
        f"BUG-16: запись 24 в MD subset не соответствует ожиданию. "
        f"Got: {line!r}"
    )


def test_bug16_entry_24_present_in_html_baseline_subset():
    """Аналогичная проверка для HTML формата."""
    # Минимальный HTML wrapper c <div class="page"> (как interactive
    # template) — гарантия что appendix-injection пройдёт по основному
    # code-path, а не fallback'у.
    html = (
        "<html><body><div class=\"page\">"
        "<h1>Th-232 demo</h1>"
        "</div></body></html>"
    )
    out = _f318_append_gost_references(html, format="html")
    assert 'class="gost-references"' in out, (
        "BUG-16: HTML appendix-секция отсутствует."
    )
    # <li value="24">...</li>
    import re
    m = re.search(r'<li\s+value="24">([^<]+)</li>', out)
    assert m is not None, (
        "BUG-16: <li value=\"24\"> (ГОСТ 26874-86) не появился в HTML "
        "bibliography subset, хотя должен быть включён через baseline_refs."
    )
    line = m.group(1)
    assert "26874" in line and "Спектрометры энергий" in line, (
        f"BUG-16: <li value=\"24\"> в HTML subset не соответствует "
        f"ожиданию. Got: {line!r}"
    )
