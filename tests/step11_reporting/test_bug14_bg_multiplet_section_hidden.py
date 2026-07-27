"""BUG-14 / v1.18.30+ (Agent B) — interactive HTML / markdown / legacy HTML
должны СКРЫВАТЬ секцию «Разложение мультиплетов» для чисто-фоновых
спектров (env == "background_only").

Symptom: при анализе фонового спектра (`bg_2016_marinelli_water_marinelli.spe`)
в интерактивном HTML-отчёте рендерилась секция:
  «Мультиплеты — разложение в фоновом спектре
   (референс для вычитания, связанная подгонка по библиотечным
    интенсивностям)»
с χ²/ν ≈ 1.00 — это статистическая подгонка к шумовому континууму
(пики не значимы), вводит оператора в заблуждение.

Fix: в `interactive_html.render_interactive_html`,
`markdown_report.build_markdown_report`,
`html_report.build_html_report` (legacy путь) — при
`diagnostics.measurement_environment == "background_only"` секция
«Мультиплеты — разложение в спектре образца» И парная bg-suffix-секция
«Мультиплеты — разложение в фоновом спектре» полностью пропускаются.

Для sample-spectrum поведение не регрессирует — оба заголовка остаются
(см. BUG-5 test в `test_interactive_report.py`).
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, "scripts")

from gamma.reporting import analyze_and_report  # noqa: E402


# Подсекции/заголовки, которые должны исчезнуть из bg-only отчёта.
# Берём подстроки из канонических заголовков _build_multiplet_blocks
# (interactive_html.py:1600-1624) + h2 из markdown_report.py:614 +
# h2 из html_report.py:625.
FORBIDDEN_HEADERS_BG_ONLY = [
    "Мультиплеты — разложение в спектре образца",
    "Мультиплеты — разложение в фоновом спектре",
    # Заголовок старого «(фон)» формата (если где-то остался)
    "Разложение мультиплетов (фон)",
]


def test_bug14_bg_only_interactive_html_omits_multiplet_section(tmp_path):
    """Background-only spectrum → interactive HTML должен НЕ содержать
    заголовков секции разложения мультиплетов (ни sample, ни bg suffix)."""
    # P1-3c: per-test tmp_path prevents xdist concurrent-write race on
    # fixed shared dir (was "demo_reports/_test_bug14_bg_only").
    out = str(tmp_path)
    bg = (
        "detectors/Gamma-1S/data/averaged_backgrounds/"
        "bg_2016_marinelli_water_marinelli.spe"
    )
    assert os.path.exists(bg), f"fixture missing: {bg}"

    res = analyze_and_report(
        bg,
        output_dir=out,
        write_html=True,
        write_markdown=True,
        write_plots=False,
    )
    html_path = res["html"]
    md_path = res["markdown"]
    assert html_path and os.path.exists(html_path), f"HTML missing: {html_path}"
    assert md_path and os.path.exists(md_path), f"MD missing: {md_path}"

    html = open(html_path, encoding="utf-8").read()
    md = open(md_path, encoding="utf-8").read()

    # JSON должен подтверждать background_only env (anti-flake предохранитель).
    json_text = open(res["json"], encoding="utf-8").read()
    assert '"measurement_environment": "background_only"' in json_text, (
        "Test fixture invariant: bg spectrum должен классифицироваться как "
        "background_only — иначе тест не валиден. См. environment.py."
    )

    # ── Interactive HTML — ни одного заголовка секции 10 ────────────
    for needle in FORBIDDEN_HEADERS_BG_ONLY:
        assert needle not in html, (
            f"BUG-14: forbidden header {needle!r} present в interactive HTML "
            f"для background_only спектра. Источник — interactive_html.py "
            f"_build_multiplet_blocks (sample/bg). Должен быть скрыт целиком."
        )

    # ── Markdown — секция 10 также пропущена ───────────────────────
    for needle in FORBIDDEN_HEADERS_BG_ONLY:
        assert needle not in md, (
            f"BUG-14: forbidden header {needle!r} present в Markdown отчёте "
            f"для background_only спектра. Источник — markdown_report.py "
            f"section 10. Должен быть скрыт целиком."
        )


def test_bug14_sample_interactive_html_keeps_multiplet_section(tmp_path):
    """Sample-spectrum (Th-232 с background-вычитанием) → секции мультиплетов
    остаются на месте — fix НЕ должен регрессировать sample report
    (см. BUG-5 в test_interactive_report.py)."""
    # P1-3c: per-test tmp_path prevents xdist concurrent-write race on
    # fixed shared dir (was "demo_reports/_test_bug14_sample").
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
    assert os.path.exists(sp), f"fixture missing: {sp}"
    assert os.path.exists(bg), f"fixture missing: {bg}"

    res = analyze_and_report(
        sp,
        output_dir=out,
        write_html=True,
        write_markdown=True,
        write_plots=False,
        sample_mass_kg=1.6,
        background_path=bg,
    )
    html_path = res["html"]
    md_path = res["markdown"]
    html = open(html_path, encoding="utf-8").read()
    md = open(md_path, encoding="utf-8").read()

    # JSON должен подтверждать НЕ background_only env (anti-flake).
    json_text = open(res["json"], encoding="utf-8").read()
    assert '"measurement_environment": "background_only"' not in json_text, (
        "Test fixture invariant: Th-232 sample не должен классифицироваться как "
        "background_only — иначе тест не валиден."
    )

    # ── Sample HTML — sample-блок ДОЛЖЕН присутствовать ─────────────
    # BUG-43 / 2026-06-04 (Agent B): bg-suffix multiplet block теперь
    # СКРЫТ во всех sample-runs (MULTIPLET_BLOCKS_BG = ""). Мультиплеты
    # относятся к анализу образца и должны показываться только в
    # sample-view. BG-view не нуждается в блоке разложения мультиплетов.
    assert "в спектре образца" in html, (
        "BUG-14 regression: sample multiplet section должна оставаться "
        "«в спектре образца» для sample report — fix не должен трогать "
        "non-background-only путь. Источник — interactive_html.py "
        "_build_multiplet_blocks(bg=False)."
    )
    # BUG-43: bg-suffix section БОЛЬШЕ НЕ рендерится в sample run —
    # suppress intentional (user request 2026-06-04).
    # Не assert-им «в фоновом спектре» in html — это теперь верное поведение.

    # ── Sample MD — секция 10 ДОЛЖНА присутствовать ────────────────
    assert "10. Мультиплеты — разложение в спектре образца" in md, (
        "BUG-14 regression: sample MD должна содержать заголовок section 10. "
        "Источник — markdown_report.py."
    )


if __name__ == "__main__":
    import tempfile, pathlib
    test_bug14_bg_only_interactive_html_omits_multiplet_section(
        pathlib.Path(tempfile.mkdtemp(prefix="_test_bug14_bg_only_"))
    )
    test_bug14_sample_interactive_html_keeps_multiplet_section(
        pathlib.Path(tempfile.mkdtemp(prefix="_test_bug14_sample_"))
    )
    print("OK")
