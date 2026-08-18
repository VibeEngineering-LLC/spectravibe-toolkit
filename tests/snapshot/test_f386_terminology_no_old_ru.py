"""F-386 — terminology rename: «двойное/одиночное ускользание» → «пик двойного/одиночного вылета».

Guard-тест проверяет, что в production-выводе (interactive HTML + Markdown отчёт)
больше НЕ встречаются устаревшие RU-термины «двойное ускользание» / «одиночное
ускользание» (применимо к pair-production escape peaks).

Семантика теста:

* Запрет (always): «двойное ускользание», «двойного ускользания»,
  «одиночное ускользание», «одиночного ускызания» — устаревший термин,
  должен быть стёрт из всех reporting-путей (build.py, interactive_html.py,
  markdown_report.py, templates/*.html).
* Conditional present: «пик двойного вылета» (соответственно «пик одиночного
  вылета») — должно появляться в выводе, когда detected secondary peak
  типа double_escape / single_escape присутствует. Если detector такой
  пик не нашёл — отсутствие normal, не fail.

Generic «escape», «ускользание» (XRF K-escape, общий счётчик
«Пиков ускользания») — НЕ часть F-386 scope и НЕ проверяется.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, "scripts")

from gamma.reporting import analyze_and_report  # noqa: E402


# Старые термины (pair-production escape peaks) — должны быть полностью удалены
_FORBIDDEN_OLD = (
    "двойное ускользание",
    "двойного ускользания",
    "одиночное ускользание",
    "одиночного ускользания",
)

# Новые термины — конкретное наличие зависит от того, классифицировал ли
# pipeline хотя бы один SE/DE peak на тестовом синтетике (Th-232 chain
# содержит линию 2614.5 кэВ → DE 1592.5 + SE 2103.5)
_NEW_DE_FULL = "пик двойного вылета"
_NEW_SE_FULL = "пик одиночного вылета"


def _strip_inline_data_arrays(html: str) -> str:
    """Удаляем inline JSON-данные в <script>, чтобы не натыкаться на
    EN-only ключи / value strings (peak_type='double_escape' etc.).

    F-386 интересует только user-facing rendered text, не data blobs.
    """
    import re
    html = re.sub(
        r"<script\b[^>]*>.*?</script>",
        "",
        html,
        flags=re.IGNORECASE | re.DOTALL,
    )
    return html


def _run_pipeline_and_get_outputs(out_dir):
    """Th-232 → analyze_and_report → (html_text, md_text).

    P1-3c: caller provides ``out_dir`` (per-test ``tmp_path``) so concurrent
    xdist workers do not clobber each other's ``report.json`` mid-write.
    The historical fixed dir was ``demo_reports/_test_f386_terminology``.
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


def test_f386_no_old_terminology_in_html(tmp_path):
    """Старые «двойное/одиночное ускользание» НЕ должны быть в HTML body."""
    html, _ = _run_pipeline_and_get_outputs(tmp_path)
    html_body = _strip_inline_data_arrays(html)
    failing = [term for term in _FORBIDDEN_OLD if term in html_body]
    assert not failing, (
        f"F-386: устаревшие RU-термины присутствуют в HTML body: {failing}. "
        f"Должны быть заменены на «пик двойного вылета» / «пик одиночного вылета»."
    )


def test_f386_no_old_terminology_in_markdown(tmp_path):
    """Старые «двойное/одиночное ускользание» НЕ должны быть в Markdown отчёте."""
    _, md = _run_pipeline_and_get_outputs(tmp_path)
    failing = [term for term in _FORBIDDEN_OLD if term in md]
    assert not failing, (
        f"F-386: устаревшие RU-термины присутствуют в Markdown: {failing}. "
        f"Должны быть заменены на «пик двойного вылета» / «пик одиночного вылета»."
    )


def test_f386_new_terminology_present_when_de_detected(tmp_path):
    """Когда DE/SE peak detected, новый термин должен появиться хотя бы в одном
    из отчётов (HTML или MD).

    На Th-232 fixture есть линия 2614.5 кэВ → ожидаем DE/SE secondary peak
    в interactive_v1_17_2.html info-span (постоянная строка) ИЛИ в data
    rows если pipeline их классифицировал.
    """
    html, md = _run_pipeline_and_get_outputs(tmp_path)
    combined = html + "\n" + md
    # Info-span в template всегда отрисовывается (постоянная подпись
    # «комптоновский край · пик одиночного / двойного вылета · …»),
    # поэтому хотя бы один из новых терминов обязан присутствовать.
    has_new = (
        _NEW_DE_FULL in combined
        or _NEW_SE_FULL in combined
        or "пик одиночного / двойного вылета" in combined
        or "пик одиночного/двойного вылета" in combined
    )
    assert has_new, (
        "F-386: ни «пик двойного вылета», ни «пик одиночного вылета», "
        "ни комбинированная info-span строка не найдены в отчётах. "
        "Ожидается хотя бы статичная подпись info-span в interactive HTML."
    )


if __name__ == "__main__":
    # P1-3c: __main__ standalone path — pass real Path tmp dirs (use stdlib).
    import tempfile, pathlib
    _a = pathlib.Path(tempfile.mkdtemp(prefix="_test_f386_a_"))
    _b = pathlib.Path(tempfile.mkdtemp(prefix="_test_f386_b_"))
    _c = pathlib.Path(tempfile.mkdtemp(prefix="_test_f386_c_"))
    test_f386_no_old_terminology_in_html(_a)
    test_f386_no_old_terminology_in_markdown(_b)
    test_f386_new_terminology_present_when_de_detected(_c)
    print("OK")
