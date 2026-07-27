"""F-114 / v1.17.4 — canonical interactive report form lock.

Verifies that `python -m gamma.cli analyze ... --full-report` (i.e.
`gamma.reporting.analyze_and_report`) emits HTML in the exact
skeleton from `references/demo_contract_v1_17_2/report.html`:

* Chart.js + chartjs-plugin-annotation script tags
* `<canvas id="fp-sp">` + `.fp-tbl` table + at least one `.fp-mp-block`
* F-113 iOS Telegram WebView fixes: `__initReport`, the
  `setTimeout(__initReport, 200)` re-entry, the Chart-missing
  fallback message, and the 1500 ms slow-CDN re-entry.
* Russian column labels (F-108 glossary).
* JSON schema/skill version v1.17.4.
* D-11: ``fp-mp-block`` appears BEFORE the summary card in HTML.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, "scripts")

from gamma.reporting import analyze_and_report  # noqa: E402


def test_th232_interactive_report_anchors(tmp_path):
    # P1-3c: per-test tmp_path prevents xdist concurrent-write race on
    # fixed shared dir (was "demo_reports/_test_interactive").
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
        sp,
        output_dir=out,
        write_html=True,
        write_plots=False,
        write_markdown=False,
        sample_mass_kg=0.5,
        background_path=bg,
    )
    html_path = res["html"]
    assert html_path and os.path.exists(html_path), f"HTML missing: {html_path}"

    with open(html_path, encoding="utf-8") as f:
        html = f.read()

    # ── form skeleton anchors ────────────────────────────────────
    form_anchors = [
        '<canvas id="fp-sp"',
        '<table class="fp-tbl">',
        'fp-mp-block',
        '__initReport',
        'setTimeout(__initReport, 200)',
        'Chart.js не загрузился',
        'chartjs-plugin-annotation',
    ]
    for anchor in form_anchors:
        assert anchor in html, f"missing form anchor: {anchor!r}"

    # ── RU column / button labels (F-108 glossary) ──────────────
    ru_labels = ["Изотоп", "Линия", "Комментарий", "Лог Y", "Сбросить выбор"]
    for ru in ru_labels:
        assert ru in html, f"missing RU label: {ru!r}"

    # ── schema/skill version bump ───────────────────────────────
    json_path = res["json"]
    with open(json_path, encoding="utf-8") as f:
        json_text = f.read()
    # v1.26.0 bump: F-300 detector-class expansion + Nuclide Pattern RAG.
    # SKILL_VERSION constant — see `scripts/gamma/reporting/json_report.py`.
    # Bump on every release in lockstep with that constant.
    from gamma.reporting.json_report import SKILL_VERSION as _SKILL_VERSION
    assert _SKILL_VERSION in json_text, f"expected {_SKILL_VERSION!r} in JSON report"

    # ── D-11: fp-mp-block must appear BEFORE the rendered summary
    # card.  "fp-summary" first appears in the <style> block; locate
    # the rendered card body marker instead (the literal label text
    # only emitted by _build_summary_card).
    idx_mp = html.find('class="fp-mp-block"')
    if idx_mp == -1:
        idx_mp = html.find("fp-mp-block")
    # Anchor the summary body by the "label" wrapper rendered into
    # the <div class="fp-summary"> container.
    idx_sum_div = html.find('<div class="fp-summary">')
    assert idx_mp != -1, "fp-mp-block missing"
    assert idx_sum_div != -1, "summary card <div> missing"
    assert idx_mp < idx_sum_div, (
        f"D-11: fp-mp-block must precede summary div "
        f"(got mp={idx_mp}, sum_div={idx_sum_div})"
    )

    # ── BUG-6: truncated rule citations must not leak as «по правилу .»
    # F-317 strip срезает bare F-ids ("F-91", "F-107", "F-89d") и оставляет
    # обрубленные фразы типа «Метод σ по правилу .». Источники этих
    # обрубков переписаны в дескриптивные формулировки (без bare F-id).
    import re as _re
    bad_pattern = _re.compile(r"[Пп]о правилу\s*\.")
    matches = bad_pattern.findall(html)
    assert not matches, (
        f"BUG-6: найдена обрубленная фраза «по правилу .» "
        f"({len(matches)} раз). F-317 strip срезал bare F-id; "
        f"источник нужно переписать дескриптивно "
        f"(см. interactive_html.py BUG-6 comments)."
    )
    # Также проверяем что fallback «по правилу d.» из F-89d не утекает.
    assert "по правилу d." not in html, (
        "BUG-6: F-89d truncation regression — bare 'F-89' срезан, "
        "оставлен суффикс 'd.'. Источник fallback в interactive_html.py / "
        "markdown_report.py должен быть дескриптивным."
    )

    # ── BUG-5: multiplet sections must have distinct role labels
    # Sample-блок размечен как «в спектре образца», bg-блок — «в фоновом
    # спектре» (если он рендерится — присутствует, когда есть bg-данные).
    assert "в спектре образца" in html, (
        "BUG-5: sample multiplet section должен явно помечен "
        "«в спектре образца» — заголовок H2 переписан в interactive_html.py."
    )
    # bg-блок присутствует при наличии bg-spectrum analysis; для Th-232
    # demo с bg-файлом он должен быть. Если bg пуст — fallback empty-block
    # тоже содержит «в фоновом спектре».
    assert "в фоновом спектре" in html, (
        "BUG-5: bg multiplet section должен явно помечен "
        "«в фоновом спектре» — заголовок H2 переписан в interactive_html.py."
    )


def test_bug6_build_multiplet_blocks_no_rule_truncation_unit():
    """BUG-6 unit-level: _build_summary_card не должен содержать bare
    F-id, который срезает strip-pipeline."""
    sys.path.insert(0, "scripts")
    from gamma.reporting.interactive_html import _build_summary_card
    fake_report = {
        "identified_nuclides": [
            {
                "nuclide": "Ac-228",
                "specific_activity_Bq_per_kg": 6167.0,
                "specific_activity_sigma_Bq_per_kg": 1026.0,
                "n_matched_lines": 14,
            }
        ]
    }
    out = _build_summary_card(fake_report)
    import re as _re
    assert not _re.search(r"[Пп]о правилу\s*F-\d", out), (
        "BUG-6: _build_summary_card до сих пор содержит bare F-id "
        "в фразе «по правилу F-NN» — будет срезано strip-pipeline."
    )
    assert not _re.search(r"[Пп]о правилу\s*\.", out), (
        "BUG-6: _build_summary_card содержит обрубленную «по правилу .»."
    )
    # Дескриптор метода σ должен быть на месте
    assert "max(" in out and ("σ_взвешенное" in out or "взвешенное среднее" in out), (
        "BUG-6: ожидается описательный дескриптор метода σ "
        "(«max(σ_взвешенное среднее, σ_разброс) …»)."
    )


if __name__ == "__main__":
    import tempfile, pathlib
    test_th232_interactive_report_anchors(pathlib.Path(tempfile.mkdtemp(prefix="_test_interactive_")))
    test_bug6_build_multiplet_blocks_no_rule_truncation_unit()
    print("OK")
