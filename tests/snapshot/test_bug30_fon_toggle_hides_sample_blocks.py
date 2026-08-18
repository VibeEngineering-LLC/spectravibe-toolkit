# -*- coding: utf-8 -*-
"""BUG-30 / v1.22.5 — «Фон» toggle скрывает блоки, актуальные только для образца.

Контракт: когда оператор переключается в view=bg, CSS-правила под body.mode-bg
скрывают:
    .fp-summary          — activity card (итоговая удельная активность)
    .fp-notes            — notes/warnings/conclusion блок
    .passport-comparison — сравнение с паспортом (добавляется build.py)

Механизм: body.mode-bg добавляется JS-функцией setView() (уже была).
Новые CSS-правила (строки с BUG-30) дополняют существующие view-sample/view-bg
правила, добавленные в F-397.

Симптом (до фикса): оператор видел activity card + passport comparison в режиме
«Фон» — визуальный шум, источник недоразумений.
"""
from __future__ import annotations

from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]


def _read_template() -> str:
    tmpl = (
        REPO / "scripts" / "gamma" / "reporting" / "templates"
        / "interactive_v1_17_2.html"
    )
    if not tmpl.is_file():
        pytest.skip(f"template missing: {tmpl}")
    return tmpl.read_text(encoding="utf-8")


# ──────────────────────────────────────────────────────────────────
# 1) CSS rules присутствуют в template
# ──────────────────────────────────────────────────────────────────

def test_bug30_css_hides_fp_summary_in_mode_bg():
    """Template содержит CSS-правило, скрывающее .fp-summary при body.mode-bg."""
    tmpl = _read_template()
    # BUG-30 rule: body.mode-bg .fp-summary { display: none; }
    assert "body.mode-bg .fp-summary" in tmpl, (
        "BUG-30 CSS rule для .fp-summary не найден в template"
    )


def test_bug30_css_hides_fp_notes_in_mode_bg():
    """Template содержит CSS-правило, скрывающее .fp-notes при body.mode-bg."""
    tmpl = _read_template()
    assert "body.mode-bg .fp-notes" in tmpl, (
        "BUG-30 CSS rule для .fp-notes не найден в template"
    )


def test_bug30_css_hides_passport_comparison_in_mode_bg():
    """Template содержит CSS-правило, скрывающее .passport-comparison при body.mode-bg."""
    tmpl = _read_template()
    assert "body.mode-bg .passport-comparison" in tmpl, (
        "BUG-30 CSS rule для .passport-comparison не найден в template"
    )


def test_bug30_css_display_none():
    """CSS-правила BUG-30 используют display:none (не visibility:hidden или другое)."""
    tmpl = _read_template()
    # Найдём блок BUG-30
    assert "BUG-30" in tmpl, "BUG-30 marker не найден в template"
    idx = tmpl.find("BUG-30")
    # Ищем display: none в следующих 400 символах после маркера
    block = tmpl[idx: idx + 400]
    assert "display: none" in block or "display:none" in block, (
        "BUG-30 CSS block должен использовать display: none"
    )


# ──────────────────────────────────────────────────────────────────
# 2) setView() по-прежнему управляет mode-bg (регрессия на JS-механизм)
# ──────────────────────────────────────────────────────────────────

def test_bug30_setview_still_adds_mode_bg_class():
    """setView(view=bg) должен добавлять класс mode-bg на body.
    Проверяем что JS-механизм не сломан после добавления CSS-правил."""
    tmpl = _read_template()
    assert "classList.add('mode-bg')" in tmpl or 'classList.add("mode-bg")' in tmpl, (
        "setView JS handler не добавляет mode-bg — регрессия после BUG-30"
    )
    assert "classList.remove('mode-bg')" in tmpl or 'classList.remove("mode-bg")' in tmpl, (
        "setView JS handler не снимает mode-bg — регрессия после BUG-30"
    )


# ──────────────────────────────────────────────────────────────────
# 3) HTML содержит блоки .fp-summary и .fp-notes (не убраны из template)
# ──────────────────────────────────────────────────────────────────

def test_bug30_template_still_has_fp_summary_div():
    """div.fp-summary должен присутствовать в template (не удалён)."""
    tmpl = _read_template()
    assert 'class="fp-summary"' in tmpl, (
        "div.fp-summary отсутствует в template — неожиданное удаление"
    )


def test_bug30_template_still_has_fp_notes_div():
    """div.fp-notes должен присутствовать в template (не удалён)."""
    tmpl = _read_template()
    assert 'class="fp-notes"' in tmpl, (
        "div.fp-notes отсутствует в template — неожиданное удаление"
    )
