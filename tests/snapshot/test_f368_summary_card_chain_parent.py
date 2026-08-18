# -*- coding: utf-8 -*-
"""F-368 / v1.18.24.3 — Summary card должна показывать ПАРЕНТ-нуклид
природной цепочки, а не daughter-line-carrier.

Регрессия: пользователь сообщил «опять потерян источник» (2026-06-01)
— заголовок секции «Итоговая удельная активность» в HTML отображал «Ac-228»
вместо «Th-232», теряя информацию о parent-нуклиде. В природной
Th-232 цепочке Ac-228 — короткоживущая дочка в секулярном равновесии;
парент Th-232 не имеет γ-линий в Gamma-1S диапазоне и измеряется
ВСЕГДА через дочки (Ac-228 / Tl-208 / Pb-212 / Bi-212). Аналогично
для Ra-226 цепочки (Pb-214 / Bi-214 → Ra-226).
"""
from __future__ import annotations
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))


def test_F368_chain_parent_mapping():
    """`_chain_parent_for(daughter)` returns the natural-chain parent."""
    from gamma.reporting.interactive_html import _chain_parent_for
    # Th-232 chain daughters
    assert _chain_parent_for("Ac-228") == "Th-232"
    assert _chain_parent_for("Tl-208") == "Th-232"
    assert _chain_parent_for("Pb-212") == "Th-232"
    assert _chain_parent_for("Bi-212") == "Th-232"
    # Ra-226 chain (post-Rn-222) daughters
    assert _chain_parent_for("Bi-214") == "Ra-226"
    assert _chain_parent_for("Pb-214") == "Ra-226"
    # Singletons / techno — no mapping (parent = self)
    assert _chain_parent_for("Cs-137") is None
    assert _chain_parent_for("Co-60") is None
    assert _chain_parent_for("K-40") is None
    assert _chain_parent_for("I-131") is None
    # Self-parent — Ra-226 measured via 186 кэВ → still Ra-226 in title,
    # but не через mapping (mapping для daughters, не self)
    assert _chain_parent_for("Ra-226") is None
    # Empty / None safety
    assert _chain_parent_for("") is None
    assert _chain_parent_for(None) is None


def test_F368_summary_card_th232_chain_shows_parent():
    """Когда best identification — Ac-228, summary card label =
    «Итоговая удельная активность Th-232» с note про Ac-228 как line carrier."""
    from gamma.reporting.interactive_html import _build_summary_card
    report = {
        "identified_nuclides": [{
            "nuclide": "Ac-228",
            "specific_activity_Bq_per_kg": 5400.0,
            "specific_activity_sigma_Bq_per_kg": 213.0,
            "activity_relative_sigma": 0.04,
            "n_matched_lines": 5,
            "is_upper_limit": False,
        }],
    }
    card = _build_summary_card(report)
    assert "Итоговая удельная активность Th-232" in card, (
        f"label должна быть Th-232 (parent цепочки), а не Ac-228: {card!r}"
    )
    assert "Итоговая удельная активность Ac-228" not in card, (
        "регрессия: показан daughter-нуклид вместо parent"
    )
    # Daughter-name carrier попадает в note
    assert "Ac-228" in card, "best line carrier (Ac-228) должна быть в note"
    assert "цепочк" in card.lower() or "цепочка" in card.lower(), (
        "должно быть упоминание цепочки/равновесия"
    )
    # Численные значения сохранены
    assert "5400" in card
    assert "213" in card


def test_F368_summary_card_ra226_chain_via_bi214():
    """Best ident по Bi-214 — должен показать Ra-226."""
    from gamma.reporting.interactive_html import _build_summary_card
    report = {
        "identified_nuclides": [{
            "nuclide": "Bi-214",
            "specific_activity_Bq_per_kg": 120.0,
            "specific_activity_sigma_Bq_per_kg": 11.0,
            "activity_relative_sigma": 0.09,
            "n_matched_lines": 4,
            "is_upper_limit": False,
        }],
    }
    card = _build_summary_card(report)
    assert "Итоговая удельная активность Ra-226" in card, (
        f"для Bi-214 должна показаться Ra-226: {card!r}"
    )


def test_F368_summary_card_singleton_keeps_nuclide_in_title():
    """Для одиночного нуклида (Cs-137) parent-mapping не применяется
    — title остаётся «Итоговая удельная активность Cs-137»."""
    from gamma.reporting.interactive_html import _build_summary_card
    report = {
        "identified_nuclides": [{
            "nuclide": "Cs-137",
            "specific_activity_Bq_per_kg": 1050.0,
            "specific_activity_sigma_Bq_per_kg": 42.0,
            "activity_relative_sigma": 0.04,
            "n_matched_lines": 1,
            "is_upper_limit": False,
        }],
    }
    card = _build_summary_card(report)
    assert "Итоговая удельная активность Cs-137" in card, (
        f"для одиночного Cs-137 title должен оставаться: {card!r}"
    )
    # «цепочка» НЕ должна упоминаться для одиночных нуклидов
    assert "секулярном равновесии" not in card, (
        "одиночный нуклид не имеет цепочки"
    )


def test_F368_th232_demo_html_shows_th232_in_summary():
    """End-to-end guard: Th-232 demo report содержит «Итоговая удельная активность
    Th-232» в HTML, НЕ «Ac-228» в этом заголовке."""
    demo = (
        REPO / "demo_reports" / "v1_18_24_th232_full" / "sample"
        / "Th232_Маринелли_0cm_report.html"
    )
    if not demo.is_file():
        pytest.skip(f"demo missing: {demo}")
    html = demo.read_text(encoding="utf-8")
    # В summary card должен быть Th-232 title
    assert "Итоговая удельная активность Th-232" in html, (
        "regression: production-demo не показывает Th-232 как parent"
    )
    # Ac-228 как daughter ОК в общем тексте (он действительно best line),
    # но НЕ должен фигурировать в качестве title итоговой активности
    assert "Итоговая удельная активность Ac-228" not in html, (
        "regression: production-demo показывает Ac-228 вместо Th-232 в title"
    )
