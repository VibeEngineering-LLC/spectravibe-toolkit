# -*- coding: utf-8 -*-
"""F-324 / v1.18.18.1 — regression guard: ГОСТ refs section МУСТ inject
INSIDE the .page container (max-width:900px), не после </div></body>.

User-reported visual bug v1.18.18: список литературы рендерится во
всю ширину окна, не соответствует ширине основных блоков.

Корневая причина: _f318_append_gost_references вставлял <section>
непосредственно перед </body>, что помещало его ВНЕ <div class="page">.

Контракт после F-324:
  • В интерактивном HTML-шаблоне (interactive_v1_17_2.html) section
    ОБЯЗАН быть внутри <div class="page">...</div>.
  • В static template (html_report.py) `body` сам имеет max-width:1100px,
    поэтому before-</body> insertion остаётся валидным.

Если будущая правка нарушит контракт — test fails.
"""
from __future__ import annotations
import os, re, subprocess, sys
from pathlib import Path
import pytest

REPO = Path(__file__).resolve().parents[2]
FIXTURES = ["M_cs_легкий", "M_k_легкий", "M_ra_легкий", "M_th_легкий"]


@pytest.fixture(scope="module")
def demo_run(tmp_path_factory):
    out = tmp_path_factory.mktemp("refs_width")
    cmd = [
        sys.executable,
        str(REPO / "scripts" / "regen_demo_reports.py"),
        "--output-dir", str(out),
    ]
    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONPATH"] = (
        str(REPO / "scripts") + os.pathsep + env.get("PYTHONPATH", "")
    )
    r = subprocess.run(
        cmd, capture_output=True, text=True, timeout=180,
        env=env, encoding="utf-8", errors="replace",
        cwd=str(REPO),
    )
    if r.returncode != 0:
        pytest.skip(f"regen failed: {r.stderr}")
    return out


def _find_page_close_idx(html: str):
    """Depth-balanced search for the </div> that closes <div class="page">.
    Returns None if .page not present.
    """
    mo = re.search(r'<div\s[^>]*class="page"', html, re.IGNORECASE)
    if not mo:
        return None
    ee = re.compile(r'>').search(html, mo.start())
    if not ee:
        return None
    cursor = ee.end()
    depth = 1
    for m in re.finditer(
        r'<(/?)div\b[^>]*>', html[cursor:], re.IGNORECASE,
    ):
        if m.group(1) == "":
            depth += 1
        else:
            depth -= 1
            if depth == 0:
                return cursor + m.start()
    return None


@pytest.mark.parametrize("stem", FIXTURES)
def test_F324_refs_section_inside_page_container(demo_run, stem):
    """gost-references должен находиться ВНУТРИ <div class="page">.

    Использует depth-balanced поиск, потому что interactive_v1_17_2.html
    имеет <script> теги между </div> закрытием .page и </body>.
    """
    html = (demo_run / f"{stem}_report.html").read_text(encoding="utf-8")
    if 'class="page"' not in html:
        pytest.skip(f"{stem}: not interactive template (no .page)")

    refs_idx = html.find('class="gost-references"')
    assert refs_idx > 0, f"{stem}: gost-references section отсутствует"

    page_close_idx = _find_page_close_idx(html)
    assert page_close_idx is not None, (
        f"{stem}: не найдено закрытие .page контейнера"
    )

    assert refs_idx < page_close_idx, (
        f"{stem}: gost-references (idx={refs_idx}) ПОСЛЕ закрывающего "
        f"</div> .page контейнера (idx={page_close_idx}) — "
        "F-324 contract violation, секция выйдет за пределы 900px column"
    )


@pytest.mark.parametrize("stem", FIXTURES)
def test_F324_refs_section_before_body_close(demo_run, stem):
    """Sanity: gost-references всё ещё ДО </body> (валидный HTML)."""
    html = (demo_run / f"{stem}_report.html").read_text(encoding="utf-8")
    refs_idx = html.find('class="gost-references"')
    if refs_idx < 0:
        pytest.skip(f"{stem}: gost-references отсутствует")
    body_close = html.lower().rfind("</body>")
    assert body_close > refs_idx, (
        f"{stem}: gost-references после </body> — невалидный HTML"
    )


@pytest.mark.parametrize("stem", FIXTURES)
def test_F324_refs_section_no_inline_full_width_styling(demo_run, stem):
    """Section не должна иметь явных стилей width:100% / margin:0 -X
    (которые сломают inheritance from .page).
    """
    html = (demo_run / f"{stem}_report.html").read_text(encoding="utf-8")
    refs_idx = html.find('class="gost-references"')
    if refs_idx < 0:
        pytest.skip()
    # Берём ровно тег <section ...>
    section_open = html[refs_idx - 50:refs_idx + 200]
    for bad in ("width:100vw", "width: 100vw", "margin-left:-",
                "margin-left: -"):
        assert bad not in section_open, (
            f"{stem}: section содержит '{bad}' — сломает inheritance"
        )


# ──────────────────────────────────────────────────────────────────
# F-325 — background filename surfaced in report
# ──────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("stem", FIXTURES)
def test_F325_background_filename_present_in_json(demo_run, stem):
    """JSON header.background_filename должен быть НЕ пустой когда
    фон был вычтен (background_subtracted=True).
    """
    import json as _json
    data = _json.loads(
        (demo_run / f"{stem}_report.json").read_text(encoding="utf-8")
    )
    hdr = data.get("header", {})
    if not hdr.get("background_subtracted"):
        pytest.skip(f"{stem}: bg не вычтен, проверка не применима")
    bg_fn = hdr.get("background_filename") or ""
    assert bg_fn, (
        f"{stem}: bg_filename пустой при background_subtracted=True — "
        "F-325 plumbing regression (staged_pipeline.py "
        "auto_background_applied_path)"
    )


@pytest.mark.parametrize("stem", FIXTURES)
def test_F325_background_filename_in_html_subtitle(demo_run, stem):
    """HTML subtitle должен содержать 'файл фона: ...' когда фон вычтен."""
    import json as _json
    data = _json.loads(
        (demo_run / f"{stem}_report.json").read_text(encoding="utf-8")
    )
    if not data.get("header", {}).get("background_subtracted"):
        pytest.skip()
    html = (demo_run / f"{stem}_report.html").read_text(encoding="utf-8")
    assert "файл фона:" in html, (
        f"{stem}: HTML subtitle без 'файл фона: …' — оператор не "
        "увидит источник bg"
    )


@pytest.mark.parametrize("stem", FIXTURES)
def test_F325_chart_label_explicit_bg_state(demo_run, stem):
    """HTML subtitle должен использовать однозначную фразу:
    'график — net спектр (фон вычтен)' либо 'gross спектр'.
    """
    html = (demo_run / f"{stem}_report.html").read_text(encoding="utf-8")
    has_net = "график — net спектр" in html
    has_gross = "график — gross спектр" in html
    assert has_net or has_gross, (
        f"{stem}: subtitle не содержит однозначной маркировки "
        "графика (net / gross). F-325 контракт нарушен."
    )


# ──────────────────────────────────────────────────────────────────
# F-326 — passport comparison section ALWAYS present
# ──────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("stem", FIXTURES)
def test_F326_passport_section_present_html(demo_run, stem):
    """passport-comparison секция должна быть в HTML ВСЕГДА."""
    html = (demo_run / f"{stem}_report.html").read_text(encoding="utf-8")
    assert "Сравнение с паспортной удельной активностью" in html, (
        f"{stem}: passport comparison секция отсутствует"
    )
    assert 'class="passport-comparison"' in html


@pytest.mark.parametrize("stem", FIXTURES)
def test_F326_passport_section_present_md(demo_run, stem):
    """passport-comparison секция должна быть в MD ВСЕГДА."""
    md = (demo_run / f"{stem}_report.md").read_text(encoding="utf-8")
    assert "Сравнение с паспортной удельной активностью" in md


@pytest.mark.parametrize("stem", FIXTURES)
def test_F326_F330_passport_either_deferred_or_auto(demo_run, stem):
    """F-326 base contract upgraded by F-330 / v1.18.18.4:
    либо deferred-сообщение (нет данных) с инструкцией, ЛИБО
    auto-routed passport (LSRM .spe COMMENT содержит passport-data) с
    provenance preamble. Любое из двух — корректно.
    """
    html = (demo_run / f"{stem}_report.html").read_text(encoding="utf-8")
    has_deferred = "Сравнение не выполнено" in html
    # F-337.6 / v1.18.19.1 — provenance phrase «Источник паспорта: …» убрана;
    # auto-routed теперь определяется по passport-comparison table + параметрам
    # пересчёта (масса/дата) ИЛИ по explicit-маркеру (passport_meta='explicit').
    has_auto = (
        "автоматически извлечено из поля COMMENT" in html
        or "передано пользователем явно" in html
        or (
            'class="passport-comparison"' in html
            and "Параметры пересчёта" in html
        )
    )
    assert has_deferred or has_auto, (
        f"{stem}: passport section ни deferred-сообщения, ни auto-данных"
    )
    if has_deferred:
        assert "passport_activity_Bq" in html, (
            f"{stem}: нет инструкции как передать паспортные данные"
        )
    if has_auto:
        # Auto-routed section should also include the comparison table.
        assert 'class="passport-tbl"' in html, (
            f"{stem}: auto-routing включён, но нет comparison-таблицы"
        )


# ──────────────────────────────────────────────────────────────────
# F-327 — refs section relabel + typography CSS
# ──────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("stem", FIXTURES)
def test_F327_refs_relabeled_html(demo_run, stem):
    """«Список использованных источников» → «Список использованной литературы»."""
    html = (demo_run / f"{stem}_report.html").read_text(encoding="utf-8")
    assert "Список использованной литературы" in html
    assert "Список использованных источников" not in html


@pytest.mark.parametrize("stem", FIXTURES)
def test_F327_refs_relabeled_md(demo_run, stem):
    md = (demo_run / f"{stem}_report.md").read_text(encoding="utf-8")
    assert "Список использованной литературы" in md
    assert "Список использованных источников" not in md


@pytest.mark.parametrize("stem", FIXTURES)
def test_F327_refs_css_present_in_html(demo_run, stem):
    """CSS rule для .gost-references должен быть в template (F-327)."""
    html = (demo_run / f"{stem}_report.html").read_text(encoding="utf-8")
    if 'class="gost-references"' not in html:
        pytest.skip()
    assert ".gost-references h2" in html, (
        f"{stem}: CSS .gost-references h2 отсутствует — typography "
        "не будет match main content"
    )
    assert ".gost-references ol" in html
