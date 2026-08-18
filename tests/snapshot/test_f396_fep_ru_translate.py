"""F-396 — RU translate of "Primary FEP" / "Secondary peaks" / "Multiplet
deconvolutions" display strings in production HTML / Markdown reports.

User feedback (п.7): "Primary FEP пиков не понятен человеку. Перевести все
надписи разделов и комментариев."

Semantics:
* Full form (section headers): «Основные пики полного поглощения»
* Compact (tables, chart legends): «Основные ФЭП»
* «Secondary peaks» (display) → «Вторичные пики»
* «Multiplet deconvolutions» (display) → «Разложение мультиплетов»

JSON keys (`primary_feps`, `secondary_peaks`, `multiplet_deconvolutions`)
и variable / function names — НЕ трогаем (API contract).
"""

from __future__ import annotations

import os
import re
import sys

sys.path.insert(0, "scripts")

from gamma.reporting import analyze_and_report  # noqa: E402


def _strip_inline_data_arrays(html: str) -> str:
    """Удаляем inline <script> блоки — там могут жить EN identifiers /
    JSON keys, которые НЕ user-facing display."""
    return re.sub(
        r"<script\b[^>]*>.*?</script>",
        "",
        html,
        flags=re.IGNORECASE | re.DOTALL,
    )


def _run_pipeline_and_get_outputs(out_dir):
    """Th-232 → analyze_and_report → (html_text, md_text).

    P1-3c: caller provides ``out_dir`` (per-test ``tmp_path``) so 9 concurrent
    xdist workers do not clobber each other's ``report.json`` mid-write
    (this was the original f396 flake — partial read missing ``primary_feps``).
    The historical fixed dir was ``demo_reports/_test_f396_fep_ru``.
    """
    out = str(out_dir)
    sp = (
        "detectors/Gamma-1S/reference_spectra/"
        "archive/"
        "Th232_420-7-17_Маринелли_0cm.spe"
    )
    bg = (
        "detectors/Gamma-1S/data/averaged_backgrounds/"
        "bg_2016_marinelli_water_marinelli.spe"
    )
    res = analyze_and_report(
        sp,
        output_dir=out,
        write_html=True,
        write_markdown=True,
        write_plots=False,
        sample_mass_kg=0.5,
        background_path=bg,
    )

    html_text = open(res["html"], encoding="utf-8").read()
    md_text = open(res["markdown"], encoding="utf-8").read()
    return html_text, md_text


# ──────────────────────────────────────────────────────────────────
# HTML production output
# ──────────────────────────────────────────────────────────────────

def test_f396_html_has_ru_peak_terminology(tmp_path):
    """В HTML body (вне <script>) должна присутствовать RU терминология
    для FEP. Интерактивный HTML использует «Изотоп / Линия / Комментарий»
    в table headers; legacy HTML — «Основные пики полного поглощения».

    Хотя бы один из этих RU-маркеров обязан быть в body.
    """
    html, _ = _run_pipeline_and_get_outputs(tmp_path)
    body = _strip_inline_data_arrays(html)
    has_ru = (
        "Основные пики полного поглощения" in body
        or "Основные ФЭП" in body
        or "Изотоп" in body                  # interactive table header
        or "характеристическ" in body.lower()  # «характеристические линии»
    )
    assert has_ru, (
        "F-396: ожидается RU терминология для FEP-пиков в HTML body "
        "(«Основные пики полного поглощения», «Основные ФЭП», "
        "«Изотоп» или «характеристические линии»). Найдено: ни одно."
    )


def test_f396_html_no_primary_fep_english_in_body(tmp_path):
    """«Primary FEP» / «Primary FEPs» как USER-FACING display НЕ должны
    появляться в HTML body (script/data blobs исключаем).

    Допускается «primary_fep» (snake_case JSON key) в data-* атрибутах,
    но НЕ «Primary FEP» с пробелом — это display фраза.
    """
    html, _ = _run_pipeline_and_get_outputs(tmp_path)
    body = _strip_inline_data_arrays(html)
    # Find display occurrences (capital P + space + capital FEP — это явно
    # user-visible header / chip text, не identifier).
    matches = re.findall(r"Primary\s+FEPs?", body)
    assert not matches, (
        f"F-396: устаревшая EN фраза «Primary FEP(s)» найдена в HTML body "
        f"({len(matches)} раз). Должна быть заменена на RU."
    )


def test_f396_html_multiplet_section_ru(tmp_path):
    """Раздел мультиплетов должен содержать RU stem «мультиплет» в
    заголовке секции. BUG-5 / v1.18.30+ — H2 переписан в форме
    «Мультиплеты — разложение в спектре образца / в фоновом спектре»;
    допускается как новая, так и устаревшая форма «Разложение мультиплетов»
    (legacy build-pipeline pass-through, html_report.py)."""
    html, _ = _run_pipeline_and_get_outputs(tmp_path)
    body = _strip_inline_data_arrays(html)
    has_ru = (
        "Мультиплеты — разложение" in body  # новая форма (BUG-5)
        or "Разложение мультиплет" in body  # legacy форма
    )
    assert has_ru, (
        "F-396: ожидается «Мультиплеты — разложение …» или "
        "«Разложение мультиплет(а|ов)» в HTML body."
    )


def test_f396_html_secondary_peaks_ru(tmp_path):
    """RU «Вторичн...» (пики/процессы) должны присутствовать как display
    в HTML. Интерактивный HTML использует «Вторичные процессы» (toggle
    UI), legacy — «Вторичные пики» (section header)."""
    html, _ = _run_pipeline_and_get_outputs(tmp_path)
    body = _strip_inline_data_arrays(html)
    has_ru = (
        "Вторичные пики" in body
        or "Вторичные процессы" in body
    )
    assert has_ru, (
        "F-396: ожидается «Вторичные пики» либо «Вторичные процессы» "
        "в HTML body."
    )


# ──────────────────────────────────────────────────────────────────
# Markdown production output
# ──────────────────────────────────────────────────────────────────

def test_f396_md_has_ru_primary_fep_header(tmp_path):
    """Markdown section 4 = «Основные пики полного поглощения»."""
    _, md = _run_pipeline_and_get_outputs(tmp_path)
    assert "Основные пики полного поглощения" in md, (
        "F-396: ожидается «4. Основные пики полного поглощения» в Markdown."
    )


def test_f396_md_no_primary_fep_english_in_body(tmp_path):
    """Markdown НЕ должен содержать EN display «Primary FEP(s)»."""
    _, md = _run_pipeline_and_get_outputs(tmp_path)
    matches = re.findall(r"Primary\s+FEPs?", md)
    assert not matches, (
        f"F-396: «Primary FEP(s)» найдено в Markdown ({len(matches)} раз)."
    )


def test_f396_md_multiplet_section_ru(tmp_path):
    """Markdown раздел 10. BUG-5 / v1.18.30+ — H2 переписан в форме
    «10. Мультиплеты — разложение в спектре образца»; допускается также
    устаревшая форма «Разложение мультиплетов» для legacy."""
    _, md = _run_pipeline_and_get_outputs(tmp_path)
    has_ru = (
        "Мультиплеты — разложение" in md  # новая форма (BUG-5)
        or "Разложение мультиплет" in md  # legacy
    )
    assert has_ru, (
        "F-396: ожидается «Мультиплеты — разложение …» или "
        "«Разложение мультиплет(а|ов)» в Markdown."
    )


def test_f396_md_secondary_peaks_ru(tmp_path):
    """Markdown раздел 5: «Вторичные пики»."""
    _, md = _run_pipeline_and_get_outputs(tmp_path)
    assert "Вторичные пики" in md, (
        "F-396: ожидается «5. Вторичные пики» в Markdown."
    )


# ──────────────────────────────────────────────────────────────────
# JSON keys preserved (API contract — must NOT be translated).
# Проверяем через JSON артефакт пайплайна, а не HTML inline data
# (production interactive HTML встраивает данные как JS-объект, ключи
# не появляются как литералы в HTML body — это deliberate optimization).
# ──────────────────────────────────────────────────────────────────

def test_f396_json_keys_still_english(tmp_path):
    """JSON keys остаются английскими (API contract).

    Проверяем через JSON-артефакт (build_json_report), а не HTML body.
    """
    import json
    from gamma.reporting import analyze_and_report

    # P1-3c: per-test tmp_path prevents xdist concurrent-write race on
    # fixed shared dir (was "demo_reports/_test_f396_fep_ru"). This is the
    # ORIGINAL flake site: 9 f396 tests writing the same dir → partial
    # report.json read missing "primary_feps".
    out = str(tmp_path)
    sp = (
        "detectors/Gamma-1S/reference_spectra/"
        "archive/"
        "Th232_420-7-17_Маринелли_0cm.spe"
    )
    bg = (
        "detectors/Gamma-1S/data/averaged_backgrounds/"
        "bg_2016_marinelli_water_marinelli.spe"
    )
    res = analyze_and_report(
        sp, output_dir=out, write_json=True,
        write_html=False, write_markdown=False, write_plots=False,
        sample_mass_kg=0.5, background_path=bg,
    )
    with open(res["json"], encoding="utf-8") as f:
        report = json.load(f)
    assert "primary_feps" in report, (
        "F-396: JSON key `primary_feps` должен сохраниться (API contract)."
    )
    assert "secondary_peaks" in report, (
        "F-396: JSON key `secondary_peaks` должен сохраниться."
    )
    assert "multiplet_deconvolutions" in report, (
        "F-396: JSON key `multiplet_deconvolutions` должен сохраниться."
    )


if __name__ == "__main__":
    # P1-3c: __main__ standalone path — pass real Path tmp dirs (use stdlib).
    import tempfile, pathlib
    _mk = lambda n: pathlib.Path(tempfile.mkdtemp(prefix=f"_test_f396_{n}_"))
    test_f396_html_has_ru_peak_terminology(_mk("a"))
    test_f396_html_no_primary_fep_english_in_body(_mk("b"))
    test_f396_html_multiplet_section_ru(_mk("c"))
    test_f396_html_secondary_peaks_ru(_mk("d"))
    test_f396_md_has_ru_primary_fep_header(_mk("e"))
    test_f396_md_no_primary_fep_english_in_body(_mk("f"))
    test_f396_md_multiplet_section_ru(_mk("g"))
    test_f396_md_secondary_peaks_ru(_mk("h"))
    test_f396_json_keys_still_english(_mk("i"))
    print("OK")
