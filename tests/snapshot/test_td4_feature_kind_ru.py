# -*- coding: utf-8 -*-
"""TD-4 / v1.18.30 (Agent B) — shared feature_kind_ru() utility guard.

Проверяет что:
1. _kind_ru.FEATURE_KIND_RU и FEATURE_KIND_RU_SHORT содержат идентичные ключи.
2. feature_kind_ru() делегирует на правильную таблицу по флагу short=.
3. Неизвестный ключ возвращается as-is (fallback).
4. F-386 hard-lock: ни одна запись не содержит слово «ускользание».
5. markdown_report._CELL_TRANSLATIONS содержит feature_kind ключи через _kind_ru.
6. plots._PLOT_KIND_RU_SHARED указывает на FEATURE_KIND_RU_SHORT (тот же объект).
"""
from __future__ import annotations
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))


class TestKindRuModule:
    def test_keys_identical_in_both_tables(self):
        from gamma.reporting._kind_ru import FEATURE_KIND_RU, FEATURE_KIND_RU_SHORT
        assert set(FEATURE_KIND_RU.keys()) == set(FEATURE_KIND_RU_SHORT.keys()), (
            "TD-4: FEATURE_KIND_RU и FEATURE_KIND_RU_SHORT должны иметь одинаковые ключи"
        )

    def test_feature_kind_ru_long_default(self):
        from gamma.reporting._kind_ru import feature_kind_ru, FEATURE_KIND_RU
        for k, v in FEATURE_KIND_RU.items():
            assert feature_kind_ru(k) == v, f"feature_kind_ru({k!r}) != {v!r}"

    def test_feature_kind_ru_short_flag(self):
        from gamma.reporting._kind_ru import feature_kind_ru, FEATURE_KIND_RU_SHORT
        for k, v in FEATURE_KIND_RU_SHORT.items():
            assert feature_kind_ru(k, short=True) == v, (
                f"feature_kind_ru({k!r}, short=True) != {v!r}"
            )

    def test_unknown_key_passthrough(self):
        from gamma.reporting._kind_ru import feature_kind_ru
        assert feature_kind_ru("nonexistent_kind") == "nonexistent_kind"
        assert feature_kind_ru("nonexistent_kind", short=True) == "nonexistent_kind"

    def test_f386_no_uskol_in_translations(self):
        """F-386 hard-lock: «вылет», не «ускользание»."""
        from gamma.reporting._kind_ru import FEATURE_KIND_RU, FEATURE_KIND_RU_SHORT
        for table_name, table in [("LONG", FEATURE_KIND_RU), ("SHORT", FEATURE_KIND_RU_SHORT)]:
            for k, v in table.items():
                assert "ускол" not in v.lower(), (
                    f"F-386 violation in _kind_ru.{table_name}[{k!r}] = {v!r}"
                )

    def test_single_escape_present(self):
        from gamma.reporting._kind_ru import FEATURE_KIND_RU
        assert "single_escape" in FEATURE_KIND_RU
        assert "вылет" in FEATURE_KIND_RU["single_escape"].lower()

    def test_double_escape_present(self):
        from gamma.reporting._kind_ru import FEATURE_KIND_RU
        assert "double_escape" in FEATURE_KIND_RU
        assert "вылет" in FEATURE_KIND_RU["double_escape"].lower()


class TestMarkdownReportIntegration:
    def test_cell_translations_has_single_escape(self):
        """markdown_report._CELL_TRANSLATIONS должен содержать single_escape из _kind_ru."""
        from gamma.reporting.markdown_report import _CELL_TRANSLATIONS
        assert "single_escape" in _CELL_TRANSLATIONS
        assert "вылет" in _CELL_TRANSLATIONS["single_escape"].lower()

    def test_cell_translations_has_compton_edge(self):
        from gamma.reporting.markdown_report import _CELL_TRANSLATIONS
        assert "compton_edge" in _CELL_TRANSLATIONS

    def test_cell_translations_has_fluorescence_shield(self):
        from gamma.reporting.markdown_report import _CELL_TRANSLATIONS
        assert "fluorescence_shield" in _CELL_TRANSLATIONS

    def test_cell_translations_f386_no_uskol(self):
        """F-386: ни в одной ячейке нет «ускол»."""
        from gamma.reporting.markdown_report import _CELL_TRANSLATIONS
        for k, v in _CELL_TRANSLATIONS.items():
            assert "ускол" not in v.lower(), (
                f"F-386 violation in _CELL_TRANSLATIONS[{k!r}] = {v!r}"
            )


class TestPlotsIntegration:
    def test_plots_shared_map_is_feature_kind_ru_short(self):
        """plots._PLOT_KIND_RU_SHARED должен быть FEATURE_KIND_RU_SHORT из _kind_ru."""
        from gamma.reporting._kind_ru import FEATURE_KIND_RU_SHORT
        from gamma.reporting import plots
        assert plots._PLOT_KIND_RU_SHARED is FEATURE_KIND_RU_SHORT, (
            "TD-4: plots._PLOT_KIND_RU_SHARED должен ссылаться на _kind_ru.FEATURE_KIND_RU_SHORT"
        )

    def test_plots_py_no_inline_plot_kind_ru(self):
        """После TD-4 inline _PLOT_KIND_RU не должен определяться в plots.py."""
        src = (REPO_ROOT / "scripts" / "gamma" / "reporting" / "plots.py").read_text(
            encoding="utf-8"
        )
        assert "_PLOT_KIND_RU = {" not in src, (
            "TD-4: inline _PLOT_KIND_RU dict найден в plots.py — должен быть импортирован из _kind_ru"
        )
