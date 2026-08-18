"""F-132 / v1.17.7 — Обязательная оценка стоимости анализа по этапам.

Контракты:
  • HTML отчёт ВСЕГДА содержит footer с итогом (даже если CLI не передал
    --cost-tokens).
  • Markdown отчёт содержит обязательный раздел «Оценка стоимости анализа»
    с поэтапной таблицей по 10 шагам Step 1..11.
  • JSON-отчёт содержит блок `cost_estimate` с полями:
      tokens_total, session_token_budget, session_pct,
      session_pct_formatted, by_stage (list), override_used, detail.
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.reporting.cost_estimator import (
    CostEstimate, StageCostEstimate,
    estimate_cost_per_stage, estimate_total_cost,
    DEFAULT_SESSION_TOKEN_BUDGET,
)


FIXTURE_DIR = (Path(__file__).parent.parent.parent / "detectors" / "Gamma-1S"
               / "reference_spectra" / "archive")
TH232_FIXTURE = FIXTURE_DIR / "Th232_420-7-17_Маринелли_0cm.spe"
CS137_FIXTURE = FIXTURE_DIR / "Cs137_420-7-14_Маринелли_0cm.spe"
K40_FIXTURE = FIXTURE_DIR / "K40_420-7-20_Маринелли_0cm.spe"


def _need(p: Path):
    if not p.exists():
        pytest.skip(f"fixture missing: {p}")


# ─── Module-level invariants ─────────────────────────────────────────

def test_default_session_token_budget_is_200000():
    """Default бюджет 5-часовой сессии = 200_000 токенов."""
    assert DEFAULT_SESSION_TOKEN_BUDGET == 200_000


def test_per_stage_returns_10_stages():
    """estimate_cost_per_stage возвращает ровно 10 этапов."""
    _need(TH232_FIXTURE)
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    r = analyze_lsrm_spe(str(TH232_FIXTURE), complete_workflow=True)
    stages = estimate_cost_per_stage(r)
    assert len(stages) == 10
    expected_ids = {
        "step_1", "step_2", "step_3_4_5", "step_5_express",
        "step_6_peak_search", "step_7_identification",
        "step_8_deconvolution", "step_9_activities",
        "step_10_residuals", "step_11_reporting",
    }
    actual = {s.stage_id for s in stages}
    assert actual == expected_ids


def test_per_stage_all_positive():
    """Каждый этап имеет tokens_total > 0 (есть базовая стоимость)."""
    _need(TH232_FIXTURE)
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    r = analyze_lsrm_spe(str(TH232_FIXTURE), complete_workflow=True)
    for s in estimate_cost_per_stage(r):
        assert s.tokens_total > 0
        assert s.tokens_baseline > 0
        assert s.tokens_complexity >= 0
        assert s.tokens_total == s.tokens_baseline + s.tokens_complexity


def test_per_stage_russian_names():
    """Все имена этапов на русском (содержат «Шаг N»)."""
    _need(TH232_FIXTURE)
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    r = analyze_lsrm_spe(str(TH232_FIXTURE), complete_workflow=True)
    for s in estimate_cost_per_stage(r):
        assert "Шаг" in s.stage_name_ru


def test_total_cost_equals_sum_of_stages():
    """Итог по умолчанию = сумма этапов."""
    _need(TH232_FIXTURE)
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    r = analyze_lsrm_spe(str(TH232_FIXTURE), complete_workflow=True)
    ce = estimate_total_cost(r)
    auto_sum = sum(s.tokens_total for s in ce.by_stage)
    assert ce.tokens_total == auto_sum
    assert not ce.override_used


def test_total_cost_override():
    """cost_tokens_override переопределяет итог, но per-stage сохраняется."""
    _need(TH232_FIXTURE)
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    r = analyze_lsrm_spe(str(TH232_FIXTURE), complete_workflow=True)
    ce = estimate_total_cost(r, cost_tokens_override=50000)
    assert ce.tokens_total == 50000
    assert ce.override_used
    auto_sum = sum(s.tokens_total for s in ce.by_stage)
    # per-stage таблица всё равно отражает авто-оценку
    assert auto_sum != ce.tokens_total


def test_session_pct_computed():
    """session_pct = 100 · tokens_total / budget."""
    _need(TH232_FIXTURE)
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    r = analyze_lsrm_spe(str(TH232_FIXTURE), complete_workflow=True)
    ce = estimate_total_cost(r, session_token_budget=100_000)
    expected = 100.0 * ce.tokens_total / 100_000
    assert abs(ce.session_pct - expected) < 1e-6


def test_complexity_grows_with_spectrum_complexity():
    """Th-232 (chain + multiplets + activities) дороже Cs-137 (моно-нуклид)."""
    _need(TH232_FIXTURE)
    _need(CS137_FIXTURE)
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    r_th = analyze_lsrm_spe(str(TH232_FIXTURE), complete_workflow=True)
    r_cs = analyze_lsrm_spe(str(CS137_FIXTURE), complete_workflow=True)
    ce_th = estimate_total_cost(r_th)
    ce_cs = estimate_total_cost(r_cs)
    assert ce_th.tokens_total > ce_cs.tokens_total, (
        f"Th-232 ({ce_th.tokens_total}) должен стоить дороже "
        f"Cs-137 ({ce_cs.tokens_total}) — больше нуклидов + мультиплеты"
    )


def test_to_dict_json_serializable():
    """CostEstimate.to_dict() даёт JSON-сериализуемый dict."""
    _need(TH232_FIXTURE)
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    r = analyze_lsrm_spe(str(TH232_FIXTURE), complete_workflow=True)
    ce = estimate_total_cost(r)
    d = ce.to_dict()
    s = json.dumps(d, ensure_ascii=False)
    assert "tokens_total" in s
    assert "by_stage" in s
    assert "session_pct_formatted" in s


# ─── Integration: JSON / Markdown / HTML отчёты ────────────────────

@pytest.fixture
def th232_outdir(tmp_path):
    """Сгенерировать report bundle на Th-232 fixture в tmp."""
    _need(TH232_FIXTURE)
    out = tmp_path / "_f132_out"
    out.mkdir()
    from gamma.reporting import analyze_and_report
    analyze_and_report(
        str(TH232_FIXTURE), output_dir=str(out),
        write_pdf=False, plot_dpi=110,
    )
    return out


def test_json_report_contains_cost_estimate(th232_outdir):
    """JSON-отчёт обязан содержать блок cost_estimate."""
    json_files = list(th232_outdir.glob("*.json"))
    assert json_files
    data = json.loads(json_files[0].read_text(encoding="utf-8"))
    ce = data.get("cost_estimate")
    assert ce is not None
    assert ce["tokens_total"] > 0
    assert ce["session_token_budget"] == DEFAULT_SESSION_TOKEN_BUDGET
    assert 0 < ce["session_pct"] < 100
    assert len(ce["by_stage"]) == 10


def test_markdown_report_has_cost_section(th232_outdir):
    """Markdown-отчёт обязан содержать раздел «Оценка стоимости анализа»."""
    md_files = list(th232_outdir.glob("*.md"))
    assert md_files
    md = md_files[0].read_text(encoding="utf-8")
    assert "Оценка стоимости анализа" in md
    # Таблица по этапам
    assert "Шаг 1" in md
    assert "Шаг 11" in md
    # Итог + единая формулировка «~N токенов или M% от бесплатной 5-часовой сессии»
    assert "Итог" in md
    assert "токенов" in md
    assert "или" in md
    assert "от бесплатной 5-часовой сессии" in md


def test_html_report_always_shows_cost_footer(th232_outdir):
    """HTML-отчёт обязан содержать footer стоимости.

    F-317 / v1.18.15: F-id strip из user-facing body убирает 'F-132' маркер.
    Test проверяет semantic content footer'a, не F-id метку.
    """
    html_files = list(th232_outdir.glob("*.html"))
    assert html_files
    html = html_files[0].read_text(encoding="utf-8")
    assert "Стоимость анализа" in html
    # Токены должны быть отрендерены (число + слово "токенов")
    import re
    m = re.search(r"~[\d ]+\s*токенов", html)
    assert m, "HTML footer должен показать число токенов"


def test_html_footer_present_even_without_cli_flags():
    """F-132 контракт: даже без --cost-tokens CLI флагов HTML footer присутствует.

    Здесь проверяется raw _build_cost_footer() output (до user-facing
    compliance pipeline) — F-id маркер ещё присутствует.
    """
    from gamma.reporting.interactive_html import _build_cost_footer
    footer = _build_cost_footer(None)
    # Раньше возвращалось "" — теперь должен быть не пустой
    assert footer != ""
    assert "F-132" in footer
    assert "Стоимость анализа" in footer


def test_html_footer_uses_provided_estimate_when_available():
    """Если cost_estimate передан, используется он."""
    from gamma.reporting.interactive_html import _build_cost_footer
    footer = _build_cost_footer({
        "tokens": 12345,
        "session_pct": "6.2% от бесплатной 5-часовой сессии",
        "detail": "manual override",
    })
    assert "12 345" in footer or "12345" in footer
    assert "6.2%" in footer
    assert "или" in footer  # F-132 формат «~N токенов или M% ...»
    assert "manual override" in footer


def test_per_stage_table_has_all_columns(th232_outdir):
    """Markdown таблица должна содержать все 6 колонок: Шаг, Этап, Базово, Сложность, Итого, Почему."""
    md_files = list(th232_outdir.glob("*.md"))
    md = md_files[0].read_text(encoding="utf-8")
    # Найти таблицу
    idx = md.find("Оценка стоимости анализа")
    section = md[idx:idx + 3000]
    for col in ["Шаг", "Этап", "Базово", "Сложность", "Итого", "Почему"]:
        assert col in section, f"missing column '{col}' in cost table"
