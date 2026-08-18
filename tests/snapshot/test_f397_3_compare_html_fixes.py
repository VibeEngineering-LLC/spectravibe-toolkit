# -*- coding: utf-8 -*-
"""F-397.3 / v1.18.28 (Agent B) — v2_compare HTML defect fixes.

Закрывает 2 реальных дефекта рендера compare HTML, найденных при
визуальной диагностике v1_18_26_1_th232 demo bundle (run_skill bundle):

Defect 1 — broken tooltip apostrophe
    title='library anchor — не fit'ena, evidence-only'
    Апостроф в `fit'ena` СОВПАДАЕТ с single-quote-обёрткой `title='...'`,
    обрывая значение атрибута на "не fit". Остаток "ena, evidence-only'>"
    парсится как additional HTML attributes — недетерминированный
    рендеринг + XSS-shape risk. Исправлено: апостроф удалён, текст —
    чистый RU «не подгоняется», обёртка переключена на двойные кавычки.

Defect 2 — empty <b></b> for cluster без cluster_id
    V2 phantom-/wide-CC кластеры не имеют стабильного M-имени; раньше
    рендерилось `<b></b>` (пустой bold-тег). Исправлено: пропускаем
    <b>-обёртку, если cluster_id пустой.
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


def _import_compare():
    import importlib.util
    p = SCRIPTS_DIR / "gen_v2_compare_th232.py"
    spec = importlib.util.spec_from_file_location("gen_v2_compare_th232", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _make_phantom_cluster(cid: str = "M1"):
    return {
        "cluster_id": cid,
        "converged": True,
        "chi2_per_dof": 2.5,
        "closure_pct": 0.5,
        "E_lo_keV": 500.0,
        "E_hi_keV": 600.0,
        "components": [
            {
                "nuclide": "Ac-228",
                "line_E_keV": 503.82,
                "deconvolved_area": 1000.0,
                "deconvolved_area_sigma": 50.0,
                "library_I_pct": 0.18,
                "peak_area_source": "fitted",
            },
            {
                "nuclide": "Ac-228",
                "line_E_keV": 523.13,
                "deconvolved_area": 0.0,
                "deconvolved_area_sigma": 0.0,
                "library_I_pct": 0.1,
                "peak_area_source": "library_anchor_phantom",
            },
        ],
    }


class TestF397_3_TooltipFix:
    def test_phantom_tooltip_has_no_broken_apostrophe(self):
        mod = _import_compare()
        html = mod._multiplet_clusters_html(
            [_make_phantom_cluster("M1")],
            badge_label="V2",
            badge_bg="#3a7",
        )
        # Старый bug: `fit'ena` рвал атрибут — больше не должен встречаться
        assert "fit'ena" not in html
        assert "fitена" not in html
        assert "fitena" not in html.lower()

    def test_phantom_tooltip_uses_clean_ru(self):
        mod = _import_compare()
        html = mod._multiplet_clusters_html(
            [_make_phantom_cluster("M1")],
            badge_label="V2",
            badge_bg="#3a7",
        )
        # Новый текст — чистый RU
        assert "не подгоняется" in html
        assert "evidence-only" in html

    def test_phantom_tooltip_attribute_uses_double_quotes(self):
        """С двойными кавычками `title="...не подгоняется..."` апостроф
        в значении НЕ обрывает атрибут (см. F-397.3 docstring)."""
        mod = _import_compare()
        html = mod._multiplet_clusters_html(
            [_make_phantom_cluster("M1")],
            badge_label="V2",
            badge_bg="#3a7",
        )
        # Должен быть title="..." (double quotes), не title='...' с
        # апострофом внутри. Проверяем что title значение начинается
        # с "library anchor" и корректно закрывается двойной кавычкой.
        m = re.search(r'title="(library anchor[^"]*)"', html)
        assert m is not None, "title attribute не найден или сломан"
        assert "evidence-only" in m.group(1), (
            f"title значение оборвано: {m.group(1)!r}"
        )


class TestF397_3_EmptyClusterIdFix:
    def test_empty_cluster_id_suppresses_b_tag(self):
        mod = _import_compare()
        cluster = _make_phantom_cluster(cid="")  # пустой ID
        html = mod._multiplet_clusters_html(
            [cluster], badge_label="V2", badge_bg="#3a7"
        )
        # Без cid не должно быть <b></b>
        assert "<b></b>" not in html, "пустой <b></b> для безымянного кластера"

    def test_present_cluster_id_renders_b_tag(self):
        mod = _import_compare()
        cluster = _make_phantom_cluster(cid="M2")
        html = mod._multiplet_clusters_html(
            [cluster], badge_label="V2", badge_bg="#3a7"
        )
        assert "<b>M2</b>" in html

    def test_combined_no_empty_b_anywhere(self):
        """Сводный анти-регресс: список с named и unnamed кластерами →
        ни одного `<b></b>` в финальном HTML."""
        mod = _import_compare()
        named = _make_phantom_cluster("M1")
        unnamed = _make_phantom_cluster(cid="")
        html = mod._multiplet_clusters_html(
            [named, unnamed, named, unnamed],
            badge_label="V2", badge_bg="#3a7"
        )
        assert "<b></b>" not in html
        # Но <b>M1</b> должно быть
        assert html.count("<b>M1</b>") == 2
