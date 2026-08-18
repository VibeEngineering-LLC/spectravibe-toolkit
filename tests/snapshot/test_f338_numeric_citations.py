"""F-338 / v1.18.20.0 — numeric citation format (затекстовые ссылки).

Per-source locator: LSRM/standards → §, ISO → п., textbooks → с.
Auto-merge соседних bare-numeric [N][M] → [N, M] ascending unique.
Mixed runs (с локаторами) НЕ мержатся.
"""
from __future__ import annotations

import sys

import pytest

sys.path.insert(0, "scripts")

from gamma.reporting.citation_translator import (
    translate_text,
    _merge_adjacent_citations,
)


# ── Per-source locator format ──────────────────────────────────────

def test_f338_lsrm_algo_uses_section_marker():
    out, _ = translate_text("[LSRM-Algo-9.3]")
    assert out == "[7, § 9.3]"


def test_f338_gilmore_uses_page_marker():
    """Gilmore — textbook, локатор «с.» (страница)."""
    out, _ = translate_text("[GILMORE-225]")
    assert out == "[19, с. 225]"


def test_f338_budyka_uses_page_marker():
    """Будыка — учебник, локатор «с.»."""
    out, _ = translate_text("[BUDYKA-47]")
    assert out == "[12, с. 47]"


def test_f338_shendrik_uses_page_marker():
    """Шендрик — учебник, локатор «с.»."""
    out, _ = translate_text("[SHENDRIK-1-12]")
    assert "[14, с. 12]" == out or out.startswith("[14, с.")


def test_f338_iso_11929_uses_clause_marker():
    """ISO 11929 — clause notation «п.»."""
    out, _ = translate_text("[ISO-11929-5.3]")
    assert out == "[2, п. 5.3]"


def test_f338_gost_uses_section_marker():
    """ГОСТ — стандарт, локатор «§»."""
    out, _ = translate_text("[GOST-7.0.5-2008-7.3]")
    assert out == "[1, § 7.3]"


def test_f338_database_no_locator():
    """Базы данных (ENSDF, NIST-XCOM, IAEA-LC) — без локатора."""
    out, _ = translate_text("[ENSDF]")
    assert out == "[4]"


def test_f338_space_between_marker_and_number():
    """Per ГОСТ Р 7.0.5–2008 §7.3: между marker и числом — пробел.

    [7, §9.3] (без пробела) — WRONG.
    [7, § 9.3] (с пробелом) — CORRECT.
    """
    out, _ = translate_text("[LSRM-Algo-9.3]")
    assert "§ 9.3" in out
    assert "§9.3" not in out


# ── Auto-merge соседних ссылок ─────────────────────────────────────

def test_f338_merge_bare_adjacent():
    """Соседние bare-numeric ссылки → один блок с запятыми."""
    out, _ = translate_text("[ENSDF][NIST-XCOM][IAEA-LC]")
    assert out == "[4, 20, 21]"


def test_f338_merge_dedup_and_sort():
    """Дубликаты убираются; порядок — ascending."""
    assert _merge_adjacent_citations("[12][7][12]") == "[7, 12]"


def test_f338_merge_keep_separate_when_locators():
    """Mixed run (с локаторами) НЕ мержится."""
    out, _ = translate_text("[LSRM-Algo-9][BUDYKA-7.5]")
    # Оба с локаторами — не сольются
    assert out == "[7, § 9][12, с. 7.5]"


def test_f338_merge_skip_when_one_has_locator():
    """Первый с локатором, второй без — НЕ сольются (per Q2 answer)."""
    out, _ = translate_text("[LSRM-Algo-9][GILMORE]")
    # GILMORE без -N → no locator → bare [19]; LSRM-Algo-9 → [7, § 9].
    # Mixed → keep separate.
    assert out == "[7, § 9][19]"


def test_f338_merge_with_optional_whitespace():
    """Соседние с пробелом тоже мержатся (но не с запятой между)."""
    assert _merge_adjacent_citations("[5] [7] [12]") == "[5, 7, 12]"
    # С запятой между — это уже manual list, не мержим
    assert _merge_adjacent_citations("[5], [7]") == "[5], [7]"


def test_f338_no_op_single_citation():
    """Одиночная [N] остаётся как есть."""
    assert _merge_adjacent_citations("Один [7].") == "Один [7]."
    out, _ = translate_text("Один [GILMORE-225].")
    assert out == "Один [19, с. 225]."


# ── Cross-cutting ──────────────────────────────────────────────────

def test_f338_no_op_text_without_citations():
    """Текст без ссылок — без изменений."""
    src = "Просто текст без ссылок 123 и абв."
    out, _ = translate_text(src)
    assert out == src


def test_f338_unresolved_prefix_preserved():
    """Неизвестный RAG-prefix → fallback на исходный текст."""
    out, unr = translate_text("[UNKNOWN-PREFIX-1]")
    assert out == "[UNKNOWN-PREFIX-1]"
    assert "UNKNOWN-PREFIX-1" in unr


def test_f338_version_bump():
    """SKILL_VERSION ≥ v1.18.20.0 (F-338 introduced bump)."""
    from gamma.reporting.json_report import SKILL_VERSION
    import re as _re
    m = _re.match(r"^v(\d+)\.(\d+)\.(\d+)(?:\.(\d+))?$", SKILL_VERSION)
    assert m
    parts = tuple(int(g) if g else 0 for g in m.groups())
    assert parts >= (1, 18, 20, 0)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
