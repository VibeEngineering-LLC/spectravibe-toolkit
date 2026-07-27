"""F-390 / v1.18.25.3+ — v2_compare report 2-column layout (smoke).

Контракт: `scripts/gen_v2_compare_th232.py` рендерит compare-отчёт в
виде шапка + 2-column grid (V2 слева, Production справа) + full-width
diff-секция внизу. Тест запускает скрипт на готовом demo-снимке и
проверяет, что HTML содержит ключевые layout-маркеры.

Зависит от наличия demo-снимка `v1_18_25_0_th232` (с sample/ +
sample_v2/ *_report.json). Если снимок отсутствует — skip.

F-384: demo_reports/ живёт вне скилла, путь резолвится через env-var
`GAMMA_DEMO_REPORTS_DIR` или fallback на сосeда `../demo_reports`.
"""
from __future__ import annotations

import json
import os
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

    tmp_root = tmp_path_factory.mktemp("f390_compare")
    dst_run = tmp_root / src_run.name
    # копируем только sample/ и sample_v2/, без v2_compare/ — он
    # должен быть сгенерирован «с нуля» нашим вызовом
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


def test_compare_html_exists_and_nonempty(compare_run: Path):
    out = compare_run / "v2_compare" / "v2_compare_report.html"
    assert out.is_file(), f"HTML не сгенерирован: {out}"
    assert out.stat().st_size > 1000, "HTML подозрительно мал"


def test_compare_html_2col_layout(compare_run: Path):
    """F-390 acceptance: 2-column grid + правильные заголовки колонок."""
    out = compare_run / "v2_compare" / "v2_compare_report.html"
    html = out.read_text(encoding="utf-8")

    # 1) CSS grid 2 колонки (V2 слева, prod справа)
    assert "grid-template-columns:1fr 1fr" in html, (
        "CSS .grid должен задавать grid-template-columns:1fr 1fr"
    )

    # 2) Заголовки колонок
    assert "V2 (experimental)" in html, "Нет заголовка «V2 (experimental)»"
    assert ">Production<" in html, "Нет заголовка «Production»"

    # 3) Каждая колонка содержит все 4 блока данных
    # F-395: блоки теперь имеют RU заголовки.
    for block in (
        "Идентифицированные нуклиды",
        "Основные пики полного поглощения",
        "Кластеры мультиплетов",
        "Вторичные пики",
    ):
        # должно встречаться минимум дважды (по разу в V2- и в Prod-колонке)
        assert html.count(block) >= 2, (
            f"Блок «{block}» должен появляться в обеих колонках, "
            f"найдено: {html.count(block)}"
        )


def test_compare_html_diff_section_present(compare_run: Path):
    """Финальная diff-секция (matched / only_v2 / only_prod) — full width."""
    html = (compare_run / "v2_compare" / "v2_compare_report.html").read_text(
        encoding="utf-8"
    )
    # F-395: «Diff пиков» переведено в «Различия пиков», «Совпавшие пики»
    # → «Совпадающие пики».
    assert "Различия пиков" in html
    assert "Совпадающие пики" in html
    # Both only-v2 и only-prod таблицы
    assert "Только в V2" in html
    assert "Только в production" in html


def test_compare_html_kpi_header_full_width(compare_run: Path):
    """Шапка с KPI-карточками должна остаться вверху как row, не в гриде."""
    html = (compare_run / "v2_compare" / "v2_compare_report.html").read_text(
        encoding="utf-8"
    )
    # kpi блок должен идти ДО первого <h2> «Сводка по конвейерам»
    kpi_idx = html.find('class="kpi"')
    grid_h2 = html.find("Сводка по конвейерам")
    assert kpi_idx > 0 and grid_h2 > kpi_idx, (
        "KPI cards должны быть в шапке, до 2-column grid"
    )


def test_compare_data_json_preserves_f380_keys(compare_run: Path):
    """F-380 hard-lock: peak rows используют ключи channel/E_keV/sigma,
    но source-данные из primary_feps остаются по контракту F-380.
    Здесь — sanity-check, что compare_data.json генерируется и не пуст."""
    data = json.loads(
        (compare_run / "v2_compare" / "compare_data.json").read_text(encoding="utf-8")
    )
    assert "meta" in data
    assert "v2_peaks" in data and "production_peaks" in data
    # хотя бы один из конвейеров что-то нашёл
    assert data["meta"]["v2_n_peaks"] + data["meta"]["production_n_peaks"] > 0
