# -*- coding: utf-8 -*-
"""F-386.1 / v1.18.28 (Agent B) — anti-regression guard для F-386 hardlock.

F-386 (v1.18.24+) фиксирует терминологию gamma-escape peaks: «пик вылета»,
не «ускользание». При sweep по reporting/* обнаружились 3 рудимента до
v1.18.28, проскочившие mass-rename F-386:

  • markdown_report.py:436   — «Пиков ускользания»
  • cost_estimator.py:270    — «(рентген.флуор. / ускользание / сумма /
                                 край)»
  • interactive_html.py:201  — `("escape", "ускользание")` translation map

Все три исправлены в этом цикле. Guard-тест предотвращает повторное
введение «ускользан*» в строки reporting-pipeline (Agent B zone). XRF
K-escape («I K-ускользание (Kα/Kβ)») — semantically different (electron
escape from K-shell, not gamma-escape) — оставлено как есть и явно
вычитается из проверки.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
REPORTING = REPO_ROOT / "scripts" / "gamma" / "reporting"


# F-386 scope: gamma-escape peak terminology only.
# XRF K-escape is a distinct physics concept (K-shell electron escape) —
# its label "I K-ускользание" is intentionally outside F-386 scope.
XRF_K_ESCAPE_OK_FRAGMENTS = (
    "I K-ускользание (Kα)",
    "I K-ускользание (Kβ)",
    "K-ускользание (Kα)",
    "K-ускользание (Kβ)",
)


def _strip_xrf_k_escape(text: str) -> str:
    """Удалить XRF K-escape подстроки, чтобы они не давали false-positive
    в F-386 guard."""
    out = text
    for s in XRF_K_ESCAPE_OK_FRAGMENTS:
        out = out.replace(s, "")
    return out


REPORTING_PY_FILES = sorted(REPORTING.rglob("*.py"))


@pytest.mark.parametrize("path", REPORTING_PY_FILES, ids=lambda p: p.name)
def test_no_stale_uskolzania_in_reporting(path: Path):
    """Ни в одном .py файле reporting/ не должно встречаться «ускользан»
    (кроме XRF K-escape, которое explicitly выводится из scope)."""
    text = path.read_text(encoding="utf-8", errors="replace")
    cleaned = _strip_xrf_k_escape(text)
    # Также вычитаем комментарии с F-386* — там слово может упоминаться
    # как объяснение чего НЕ надо делать.
    cleaned = re.sub(r"#[^\n]*F-386[^\n]*\n", "\n", cleaned)
    # Docstring-mentions of F-386 (multiline тоже)
    # Простая эвристика: убираем строки, начинающиеся с "F-386"
    cleaned = re.sub(r"F-386[^\n]*", "", cleaned)
    assert "ускользан" not in cleaned, (
        f"F-386 violation: stale «ускользан*» in {path.name}:\n"
        + "\n".join(f"  {i+1}: {ln}" for i, ln in enumerate(text.splitlines())
                    if "ускользан" in ln
                    and not any(x in ln for x in XRF_K_ESCAPE_OK_FRAGMENTS)
                    and "F-386" not in ln)
    )


def test_no_stale_uskolzania_in_templates():
    """Templates HTML (interactive_v1_17_2.html) тоже под F-386 guard."""
    template = REPORTING / "templates" / "interactive_v1_17_2.html"
    if not template.exists():
        pytest.skip("template missing")
    text = template.read_text(encoding="utf-8", errors="replace")
    cleaned = _strip_xrf_k_escape(text)
    assert "ускользан" not in cleaned


def test_markdown_diagnostic_table_uses_vylet():
    """Sanity: в markdown_report.py диагностика строит «Пиков вылета»."""
    md = REPORTING / "markdown_report.py"
    text = md.read_text(encoding="utf-8")
    assert '"Пиков вылета"' in text or "'Пиков вылета'" in text
    assert "Пиков ускользан" not in text  # старая форма должна быть удалена


def test_interactive_html_escape_label_is_vylet():
    """Sanity: translation map (escape → вылет)."""
    ih = REPORTING / "interactive_html.py"
    text = ih.read_text(encoding="utf-8")
    assert '("escape", "вылет")' in text
    assert '("escape", "ускользание")' not in text


def test_cost_estimator_step10_uses_vylet():
    """Sanity: stage_name_ru использует «вылет»."""
    ce = REPORTING / "cost_estimator.py"
    text = ce.read_text(encoding="utf-8")
    # Стадия Step-10 классификация остаточных пиков
    assert "(рентген.флуор. / вылет / сумма / край)" in text
    assert "(рентген.флуор. / ускользание / сумма / край)" not in text
