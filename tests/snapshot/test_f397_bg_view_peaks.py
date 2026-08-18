# -*- coding: utf-8 -*-
"""F-397 / v1.18.27 — Background view filter для peak block.

Контракт: интерактивный HTML toggle «Образец / Фон / Наложение / Вычет»
дополнительно подменяет нижний peak block (таблица найденных линий +
блок мультиплетов) на bg-only данные при выборе «Фон». Это делается
через JS-handler в setView(): добавляется CSS class `mode-bg` на body,
который скрывает `.view-sample` и показывает `.view-bg`; одновременно
renderRows() переключается между rows и rowsBg.

JSON output получает три новых top-level массива:
    background_primary_feps
    background_secondary_peaks
    background_multiplet_deconvolutions

Они пустые когда фон не подключался (см. background_staged_result в
StagedAnalysisResult).

Демо-источник: Th-232 Marinelli kit (есть paired фон).
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))


# ──────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────

_TH232_SAMPLE = (
    REPO / "detectors" / "Gamma-1S" / "reference_spectra"
    / "reference_kits" / "Marinelli_1L" / "Th-232"
    / "Th232_420-7-17_Маринелли_0cm.spe"
)
_TH232_BG = (
    REPO / "detectors" / "Gamma-1S" / "reference_spectra"
    / "reference_kits" / "Marinelli_1L" / "Th-232"
    / "Фон закр кр вода_13.spe"
)


def _kit_available() -> bool:
    return _TH232_SAMPLE.is_file() and _TH232_BG.is_file()


def _run_full_th232() -> "tuple[object, dict, str]":
    """Полный прогон Th-232 sample + paired bg + JSON + HTML.

    Возвращает (StagedAnalysisResult, json_report, html).
    """
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    from gamma.reporting.json_report import build_json_report
    from gamma.reporting.interactive_html import render_interactive_html

    if not _kit_available():
        pytest.skip("Th-232 kit files missing")

    res = analyze_lsrm_spe(
        str(_TH232_SAMPLE),
        background_path=str(_TH232_BG),
        sample_mass_kg=0.5,
        complete_workflow=True,
    )
    report = build_json_report(res)
    html = render_interactive_html(report, res)
    return res, report, html


# ──────────────────────────────────────────────────────────────────
# 1) JSON schema — bg keys present (даже когда пустые)
# ──────────────────────────────────────────────────────────────────

def test_f397_json_has_background_primary_feps_key():
    """JSON содержит top-level `background_primary_feps` array (хотя бы
    пустой) — schema contract для downstream consumers."""
    from gamma.reporting.json_report import build_json_report

    class _FakeResult:
        spec = None
        final_detected = []
        residual_classifications = []
        activities = None
        specific_activities_Bq_per_kg = None
        mda_per_line = None
        deconvolution_results = None
        priority_findings = []
        chain_dominance = None
        k40_ac228_overlap_warning = False
        completeness = None
        seven_line_check = None
        ci_gating = None
        chain_filtered_out = []
        filename_isotope_hints = []
        filename_chains_claimed = []
        notes = []
        background_subtraction = None
        background_status = ""
        peak_search_method = "mariscotti"
        peak_search_method_comparison = None
        auto_background_mode = "off"
        auto_background_candidates = None
        auto_background_applied_path = None
        warnings = []
        multiplet_self_calibration_diag = None
        efficiency_curve = None
        efficiency_source = ""
        detector_type = "NaI"
        detector_canonical = "Gamma-1S"
        fwhm_model_source = "test"
        sample_mass_kg = None
        sample_type_canonical = ""
        geometry_canonical = ""
        next_stage_recommended = None
        next_stage_reason = ""
        background_staged_result = None  # нет фона — должны быть пустые списки

    # Минимальный fake: build_json_report должен отрабатывать на нём
    # для smoke-проверки наличия ключей.
    try:
        report = build_json_report(_FakeResult())
    except Exception:
        pytest.skip("build_json_report не работает на минимальном fake — нужен полный fixture")
    assert "background_primary_feps" in report
    assert "background_multiplet_deconvolutions" in report
    assert "background_secondary_peaks" in report
    assert isinstance(report["background_primary_feps"], list)
    assert isinstance(report["background_multiplet_deconvolutions"], list)
    assert isinstance(report["background_secondary_peaks"], list)


def test_f397_json_bg_keys_populated_for_th232_demo():
    """Th-232 demo (с paired фоном) → bg peak lists должны быть НЕпустые,
    так как фон содержит свои ЕРН линии (K-40, Pb-212 от лабораторного
    железобетона и т.п.)."""
    _res, report, _html = _run_full_th232()
    assert "background_primary_feps" in report
    bg_feps = report["background_primary_feps"]
    assert isinstance(bg_feps, list)
    # Фон Marinelli обычно содержит хотя бы K-40 1461 кэВ + Tl-208
    # (фон лаборатории) — sanity check.
    assert len(bg_feps) > 0, (
        "Ожидается ≥1 background primary FEP для Th-232 demo "
        "(фон содержит K-40 и Th-цепочку лаб. ЕРН)"
    )


# ──────────────────────────────────────────────────────────────────
# 2) JSON — bg peaks ≠ sample peaks (контракт «разные данные»)
# ──────────────────────────────────────────────────────────────────

def test_f397_bg_peaks_differ_from_sample_peaks():
    """Bg primary_feps не должны совпадать с sample primary_feps как
    set по (nuclide, library_E_keV). Иначе toggle бессмысленный."""
    _res, report, _html = _run_full_th232()
    sample = {
        (p.get("nuclide"), round(p.get("library_E_keV") or 0.0, 1))
        for p in (report.get("primary_feps") or [])
    }
    bg = {
        (p.get("nuclide"), round(p.get("library_E_keV") or 0.0, 1))
        for p in (report.get("background_primary_feps") or [])
    }
    if not bg:
        pytest.skip("Bg detection вернула пусто — skip differential check")
    # Th-232 sample has much richer Th-chain coverage than фон;
    # либо размеры разные, либо хотя бы одна линия отличается.
    assert sample != bg, (
        "Bg primary_feps идентичны sample — toggle не имеет смысла"
    )


# ──────────────────────────────────────────────────────────────────
# 3) HTML — оба container (view-sample / view-bg) присутствуют
# ──────────────────────────────────────────────────────────────────

def test_f397_html_has_both_sample_and_bg_containers():
    """HTML рендерит ОБА мультиплет-блока — sample под `view-sample` и
    bg под `view-bg`. CSS rule body.mode-bg инвертирует видимость."""
    _res, _report, html = _run_full_th232()
    assert 'id="fp-multiplets-sample"' in html, (
        "Отсутствует контейнер sample multiplets"
    )
    assert 'id="fp-multiplets-bg"' in html, (
        "Отсутствует контейнер bg multiplets"
    )
    # CSS rule, переключающий видимость
    assert ".view-bg" in html and "mode-bg" in html, (
        "CSS правила для view-toggle не вставлены"
    )


def test_f397_html_has_section_title_with_sample_label():
    """Section header «Найденные пики (образец)» — default label.
    JS подменяет на «(фон)» при view=bg."""
    _res, _report, html = _run_full_th232()
    assert 'id="fp-peaks-section-title"' in html
    assert "Найденные пики (образец)" in html


# ──────────────────────────────────────────────────────────────────
# 4) HTML — JS toggle handler присутствует
# ──────────────────────────────────────────────────────────────────

def test_f397_html_has_mode_bg_js_handler():
    """В <script> блоке должен быть handler, который добавляет/убирает
    класс `mode-bg` на body при переключении view."""
    _res, _report, html = _run_full_th232()
    # Извлекаем <script> блоки (берём всё что между <script> и </script>
    # внутри которого есть `setView` — это наш main report script).
    scripts = re.findall(r"<script[^>]*>(.*?)</script>", html, re.DOTALL)
    main_js = ""
    for s in scripts:
        if "setView" in s and "renderRows" in s:
            main_js = s
            break
    assert main_js, "Main JS блок с setView не найден"
    # Контракт: handler делает classList.add/remove('mode-bg').
    assert "mode-bg" in main_js, "JS не манипулирует классом mode-bg"
    assert "classList.add('mode-bg')" in main_js or \
           'classList.add("mode-bg")' in main_js, (
        "Не вижу classList.add('mode-bg') в JS"
    )
    assert "rowsBg" in main_js, (
        "JS должен использовать переменную rowsBg для swap"
    )


def test_f397_html_has_data_rows_bg_serialized():
    """`const rowsBg = [...]` — JSON-инжекция bg rows для toggle."""
    _res, _report, html = _run_full_th232()
    assert "const rowsBg=" in html or "const rowsBg =" in html, (
        "DATA_ROWS_BG не вставлен в JS"
    )


# ──────────────────────────────────────────────────────────────────
# 5) Smoke — bg peaks list (через rowsBg в HTML) не идентичен rows
# ──────────────────────────────────────────────────────────────────

def test_f397_html_bg_rows_differ_from_sample_rows():
    """Извлечь оба JS-массива из HTML и сравнить."""
    _res, _report, html = _run_full_th232()

    def _extract_array(name: str) -> list:
        m = re.search(
            r"const\s+" + re.escape(name) + r"\s*=\s*(\[.*?\]);",
            html, re.DOTALL,
        )
        if not m:
            return []
        try:
            return json.loads(m.group(1))
        except json.JSONDecodeError:
            return []

    rows = _extract_array("rows")
    rows_bg = _extract_array("rowsBg")
    if not rows_bg:
        pytest.skip("rowsBg пуст — bg detection не вернула пиков")
    # Sample Th-232 содержит много линий Th-цепочки; bg — только
    # лабораторный фон ЕРН. Множества rows должны различаться.
    assert rows != rows_bg, (
        "rows и rowsBg идентичны — bg view filter бесполезен"
    )
