"""F-395 — RU translate of compare report (v2_compare_report.html).

User feedback (п.6): «Блок «Diff пиков» не понятен человеку. Перевести все
надписи разделов и комментариев.»

Семантика:
* «Diff пиков» → «Различия пиков» (section header)
* «matched» (как display) → «Совпадающие»
* «only_v2» (как display) → «Только в V2»
* «only_prod» (как display) → «Только в production»
* «Identified nuclides» / «Primary FEP peaks» / «Multiplet clusters» /
  «Secondary peaks» — переведены на RU.
* «Прочие метрики» → «Метрики сравнения».

JSON keys в `compare_data.json` (`matched_peaks`, `only_in_v2`,
`only_in_production`, `peak_E_keV`, `peak_channel`) — НЕ трогаем
(data contract / F-380 hard-lock).

F-390 layout regression (`grid-template-columns:1fr 1fr`) — должен
сохраняться.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "gen_v2_compare_th232.py"


def _resolve_demo_root() -> Path | None:
    env = os.environ.get("GAMMA_DEMO_REPORTS_DIR")
    if env:
        p = Path(env)
        if p.is_dir():
            return p
    # F-384 fallback: соседняя папка рядом со скиллом
    sibling = REPO_ROOT.parent / "demo_reports"
    if sibling.is_dir():
        return sibling
    return None


def _find_demo_run(demo_root: Path) -> Path | None:
    """Берём первый run с sample/*_report.json и sample_v2/*_report.json."""
    for cand in sorted(demo_root.iterdir(), reverse=True):
        if not cand.is_dir():
            continue
        s = cand / "sample"
        v2 = cand / "sample_v2"
        if not (s.is_dir() and v2.is_dir()):
            continue
        if list(s.glob("*_report.json")) and list(v2.glob("*_report.json")):
            return cand
    return None


@pytest.fixture(scope="module")
def compare_run(tmp_path_factory) -> Path:
    """Копируем demo-run во временную папку и запускаем генератор там,
    чтобы не мутировать реальный demo_reports/."""
    demo_root = _resolve_demo_root()
    if demo_root is None:
        pytest.skip("demo_reports/ не найден (set GAMMA_DEMO_REPORTS_DIR)")
    src_run = _find_demo_run(demo_root)
    if src_run is None:
        pytest.skip(f"в {demo_root} нет run с sample/+sample_v2/*_report.json")

    tmp_root = tmp_path_factory.mktemp("f395_compare_ru")
    dst_run = tmp_root / src_run.name
    (dst_run).mkdir()
    for sub in ("sample", "sample_v2"):
        shutil.copytree(src_run / sub, dst_run / sub)

    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), str(dst_run)],
        capture_output=True, text=True, env=env, encoding="utf-8",
    )
    assert proc.returncode == 0, (
        f"gen_v2_compare_th232.py failed:\nstdout={proc.stdout}\nstderr={proc.stderr}"
    )
    return dst_run


def _read_html(run: Path) -> str:
    return (run / "v2_compare" / "v2_compare_report.html").read_text(encoding="utf-8")


# ──────────────────────────────────────────────────────────────────
# RU display strings present
# ──────────────────────────────────────────────────────────────────

def test_f395_diff_section_ru(compare_run: Path):
    """«Различия пиков» должно появляться вместо «Diff пиков»."""
    html = _read_html(compare_run)
    assert "Различия пиков" in html, (
        "F-395: ожидается «Различия пиков» в HTML body."
    )


def test_f395_matched_label_ru(compare_run: Path):
    """«Совпадающие» (как display label) должно появляться в HTML."""
    html = _read_html(compare_run)
    assert "Совпадающие" in html, (
        "F-395: ожидается «Совпадающие пики» (RU display) в HTML body."
    )


def test_f395_only_labels_ru(compare_run: Path):
    """«Только в V2» и «Только в production» (как display) в HTML."""
    html = _read_html(compare_run)
    assert "Только в V2" in html, "F-395: ожидается «Только в V2» в HTML body."
    assert "Только в production" in html, (
        "F-395: ожидается «Только в production» в HTML body."
    )


def test_f395_column_blocks_ru(compare_run: Path):
    """4 блока внутри каждой колонки — на RU."""
    html = _read_html(compare_run)
    for block in (
        "Идентифицированные нуклиды",
        "Основные пики полного поглощения",
        "Кластеры мультиплетов",
        "Вторичные пики",
    ):
        assert block in html, f"F-395: блок «{block}» не найден в HTML body."


def test_f395_metrics_section_ru(compare_run: Path):
    """«Метрики сравнения» должно появляться (вместо «Прочие метрики»)."""
    html = _read_html(compare_run)
    assert "Метрики сравнения" in html, (
        "F-395: ожидается заголовок «Метрики сравнения» в HTML body."
    )


# ──────────────────────────────────────────────────────────────────
# EN display strings absent (data keys в JSON — это OK,
# но в HTML body как display они появляться не должны)
# ──────────────────────────────────────────────────────────────────

def test_f395_no_diff_peaks_en_in_html(compare_run: Path):
    """«Diff пиков» (старая фраза) НЕ должно быть в HTML body."""
    html = _read_html(compare_run)
    assert "Diff пиков" not in html, (
        "F-395: устаревшая фраза «Diff пиков» найдена в HTML."
    )


def test_f395_no_en_block_headers_in_html(compare_run: Path):
    """EN блоки «Identified nuclides» / «Primary FEP peaks» /
    «Multiplet clusters» / «Secondary peaks» — НЕ в HTML body
    (это display, JSON keys не появляются как user-visible текст)."""
    html = _read_html(compare_run)
    for en in (
        "Identified nuclides",
        "Primary FEP peaks",
        "Multiplet clusters",
        "Secondary peaks",
    ):
        assert en not in html, (
            f"F-395: устаревший EN заголовок «{en}» найден в HTML body."
        )


def test_f395_no_only_v2_only_prod_keys_as_display(compare_run: Path):
    """`only_v2` / `only_prod` (snake_case JSON keys) НЕ должны
    появляться как user-visible display в HTML.

    Иначе говоря: «only_v2» как **literal** строка не должна
    встречаться в HTML body (это data contract key,
    а не display label)."""
    html = _read_html(compare_run)
    # Allow в комментариях/JSON-inline нет — current HTML их и не содержит.
    assert "only_v2" not in html, (
        "F-395: «only_v2» (JSON key) утёк в HTML display."
    )
    assert "only_prod" not in html, (
        "F-395: «only_prod» (JSON key) утёк в HTML display."
    )


def test_f395_no_matched_en_display_in_html(compare_run: Path):
    """«matched» (как literal EN display) НЕ должно встречаться.

    Note: позволяем «matched_peaks» (JSON key) — это data contract,
    но JSON в HTML не встраивается, так что просто проверяем,
    что слово «matched» в latin shape не появляется в body как
    display, например в «Совпавшие пики (matched, |ΔE|...)».
    """
    html = _read_html(compare_run)
    matches = re.findall(r"\bmatched\b", html)
    assert not matches, (
        f"F-395: «matched» (EN display) найдено в HTML body "
        f"{len(matches)} раз."
    )


# ──────────────────────────────────────────────────────────────────
# F-390 regression: 2-col grid сохраняется
# ──────────────────────────────────────────────────────────────────

def test_f395_f390_grid_layout_preserved(compare_run: Path):
    """F-390 hard-lock: `grid-template-columns:1fr 1fr` остаётся в CSS."""
    html = _read_html(compare_run)
    assert "grid-template-columns:1fr 1fr" in html, (
        "F-395/F-390: 2-col grid layout сломан."
    )


# ──────────────────────────────────────────────────────────────────
# F-380 keys сохранены в compare_data.json
# ──────────────────────────────────────────────────────────────────

def test_f395_f380_keys_preserved_in_compare_data(compare_run: Path):
    """F-380 hard-lock: ключи `peak_E_keV`, `peak_channel` остаются
    в peaks источника (primary_feps). Здесь проверяем структуру
    compare_data.json и наличие matched_peaks / only_in_v2 /
    only_in_production (data contract — НЕ переименованы)."""
    data = json.loads(
        (compare_run / "v2_compare" / "compare_data.json").read_text(
            encoding="utf-8"
        )
    )
    # Compare data structure preserved (data contract)
    assert "meta" in data
    assert "matched_peaks" in data, (
        "F-395: data key `matched_peaks` пропал — это data contract."
    )
    assert "only_in_v2" in data, (
        "F-395: data key `only_in_v2` пропал — это data contract."
    )
    assert "only_in_production" in data, (
        "F-395: data key `only_in_production` пропал — это data contract."
    )
    # meta cells preserved
    assert "n_matched" in data["meta"]
    assert "n_only_v2" in data["meta"]
    assert "n_only_prod" in data["meta"]


if __name__ == "__main__":
    print("OK (run via pytest)")
