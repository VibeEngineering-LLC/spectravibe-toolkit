"""BUG-12 / v1.18.30 — pipeline-summary-grid выравнивание секций.

Контракт: `scripts/gen_v2_compare_th232.py` рендерит блок «Сводка по
конвейерам» как CSS Grid с явными рядами, где каждая секция (nuclides /
primary peaks / multiplets / secondary peaks) занимает один ряд из двух
ячеек (V2 + Production). Это гарантирует, что заголовок одной секции
располагается в обеих колонках на одной горизонтальной линии независимо
от длины списков (V2 = 21 пик, Prod = 12 — характерный демо-случай).

Тест строит JSON-фикстуры с заведомо разной длиной столбцов
(5 пиков в V2, 2 пика в Prod), запускает генератор, и проверяет
структуру HTML:
  * существует контейнер `.pipeline-summary-grid`
  * для каждой из 4 секций в нём есть две ячейки `.ps-cell`
    с одинаковым `data-section` и противоположным `data-col`
  * заголовки секций (`<h3 class='ps-section-header'>`) парные
  * порядок ячеек в DOM построчный (V2-секция-K идёт прямо
    перед Prod-секцией-K, не блоком по колонкам) — что даёт CSS Grid
    auto-flow row-первое поведение и выравнивает заголовки.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "gen_v2_compare_th232.py"


SECTIONS = ("nuclides", "primary", "multiplets", "secondary")


def _make_report(*, n_peaks: int, n_nuc: int, n_mult: int, n_sec: int,
                 nuclide_prefix: str = "Tl") -> dict:
    """Minimal JSON-report payload с управляемыми длинами списков.

    F-380 contract: primary_feps использует ключи peak_channel /
    peak_E_keV / peak_area_counts / peak_area_sigma.
    """
    feps = []
    for k in range(n_peaks):
        feps.append({
            "peak_channel": 100 + k * 20,
            "peak_E_keV": 200.0 + k * 30.0,
            "peak_area_counts": 1000.0 * (k + 1),
            "peak_area_sigma": 50.0,
            "nuclide": f"{nuclide_prefix}-{208 + k}",
            "library_E_keV": 200.0 + k * 30.0,
            "library_I_pct": 5.0,
        })
    nuclides = [{"nuclide": f"{nuclide_prefix}-{208 + k}"} for k in range(n_nuc)]
    mults = []
    for k in range(n_mult):
        mults.append({
            "cluster_id": f"M{k + 1}",
            "converged": True,
            "chi2_per_dof": 1.2,
            "closure_pct": 95,
            "E_lo_keV": 200.0 + k * 50.0,
            "E_hi_keV": 240.0 + k * 50.0,
            "components": [{
                "nuclide": f"{nuclide_prefix}-{208 + k}",
                "line_E_keV": 215.0 + k * 50.0,
                "deconvolved_area": 500.0,
                "deconvolved_area_sigma": 25.0,
                "library_I_pct": 10.0,
            }],
        })
    secs = []
    for k in range(n_sec):
        secs.append({
            "channel": 800 + k * 30,
            "energy_keV": 1500.0 + k * 100.0,
            "significance": 4.5,
            "feature_kind": "compton_edge",
            "parent_nuclide": f"{nuclide_prefix}-{208 + k}",
            "parent_line_keV": 1800.0 + k * 100.0,
        })
    return {
        "skill_version": "v1.18.30-bug12-test",
        "header": {
            "sample_filename": "fixture.spe",
            "background_filename": "fixture_bg.spe",
            "sample_mass_kg": 1.6,
        },
        "calibration": {"sample_mass_kg": 1.6},
        "identified_nuclides": nuclides,
        "primary_feps": feps,
        "multiplet_deconvolutions": mults,
        "secondary_peaks": secs,
        "unidentified_peaks": [],
    }


@pytest.fixture(scope="module")
def mismatched_compare_html(tmp_path_factory) -> str:
    """Запускает gen_v2_compare_th232.py на фикстурах, где V2 заметно
    длиннее Production во всех секциях. Возвращает HTML."""
    root = tmp_path_factory.mktemp("bug12_align")
    run = root / "run"
    (run / "sample").mkdir(parents=True)
    (run / "sample_v2").mkdir(parents=True)

    # V2 — длинная колонка, Prod — короткая (BUG-12 reproduction)
    v2_report = _make_report(n_peaks=5, n_nuc=4, n_mult=3, n_sec=4)
    prod_report = _make_report(n_peaks=2, n_nuc=2, n_mult=1, n_sec=1)

    (run / "sample" / "fixture_report.json").write_text(
        json.dumps(prod_report, ensure_ascii=False), encoding="utf-8"
    )
    (run / "sample_v2" / "fixture_v2_report.json").write_text(
        json.dumps(v2_report, ensure_ascii=False), encoding="utf-8"
    )

    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), str(run)],
        capture_output=True, text=True, env=env, encoding="utf-8",
    )
    assert proc.returncode == 0, (
        f"gen_v2_compare_th232.py failed:\n"
        f"stdout={proc.stdout}\nstderr={proc.stderr}"
    )
    html = (run / "v2_compare" / "v2_compare_report.html").read_text(
        encoding="utf-8"
    )
    return html


def test_pipeline_summary_grid_container_present(mismatched_compare_html: str):
    """Контейнер `.pipeline-summary-grid` должен присутствовать в HTML."""
    html = mismatched_compare_html
    assert 'class="pipeline-summary-grid"' in html, (
        "Ожидается <div class='pipeline-summary-grid'> для выравнивания "
        "секций V2/Prod по горизонтали"
    )
    # CSS rule для grid-template-columns:1fr 1fr тоже должна быть
    assert ".pipeline-summary-grid{display:grid" in html.replace(" ", ""), (
        "CSS правило `.pipeline-summary-grid{display:grid; ...}` отсутствует"
    )


def test_each_section_has_paired_cells(mismatched_compare_html: str):
    """Для каждой из 4 секций должна быть ровно одна V2-ячейка и одна
    Prod-ячейка с совпадающим data-section."""
    html = mismatched_compare_html
    for section in SECTIONS:
        v2_pat = re.compile(
            r"<div\s+class='ps-cell ps-cell-v2'\s+"
            r"data-section='" + re.escape(section) + r"'\s+"
            r"data-col='v2'>"
        )
        prod_pat = re.compile(
            r"<div\s+class='ps-cell ps-cell-prod'\s+"
            r"data-section='" + re.escape(section) + r"'\s+"
            r"data-col='prod'>"
        )
        n_v2 = len(v2_pat.findall(html))
        n_prod = len(prod_pat.findall(html))
        assert n_v2 == 1, (
            f"Ожидается ровно 1 V2 ячейка для секции '{section}', найдено {n_v2}"
        )
        assert n_prod == 1, (
            f"Ожидается ровно 1 Prod ячейка для секции '{section}', "
            f"найдено {n_prod}"
        )


def test_section_cell_pairs_are_adjacent_in_dom(mismatched_compare_html: str):
    """V2-ячейка секции K должна идти в DOM непосредственно перед Prod-ячейкой
    той же секции — это row-первый порядок auto-flow grid, выравнивающий
    заголовки горизонтально.

    Проверяем по позиции вхождения data-section атрибутов в общем HTML."""
    html = mismatched_compare_html
    # Извлекаем все ячейки в порядке их появления
    cell_pat = re.compile(
        r"<div\s+class='ps-cell ps-cell-(v2|prod)'\s+"
        r"data-section='(\w+)'\s+data-col='(v2|prod)'>"
    )
    matches = cell_pat.findall(html)
    # Ожидаем последовательность: (v2, sec1), (prod, sec1), (v2, sec2), ...
    assert len(matches) == 2 * len(SECTIONS), (
        f"Ожидается {2*len(SECTIONS)} ячеек, найдено {len(matches)}"
    )
    for i, section in enumerate(SECTIONS):
        v2_cell = matches[2 * i]
        prod_cell = matches[2 * i + 1]
        assert v2_cell[0] == "v2" and v2_cell[1] == section, (
            f"Позиция {2*i}: ожидается (v2, {section}), получено {v2_cell}"
        )
        assert prod_cell[0] == "prod" and prod_cell[1] == section, (
            f"Позиция {2*i+1}: ожидается (prod, {section}), получено {prod_cell}"
        )


def test_section_headers_marked_with_class(mismatched_compare_html: str):
    """Заголовок каждой секции — <h3 class='ps-section-header'> с
    data-section атрибутом. Это позволяет JS/CSS целевой работы с
    заголовками без зависимости от текста."""
    html = mismatched_compare_html
    header_pat = re.compile(
        r"<h3\s+class='ps-section-header'\s+"
        r"data-section='(\w+)'\s+data-col='(v2|prod)'>"
    )
    headers = header_pat.findall(html)
    # 4 секции × 2 колонки = 8 заголовков
    assert len(headers) == 8, (
        f"Ожидается 8 заголовков секций (4 × V2/Prod), найдено {len(headers)}"
    )
    # Каждая секция представлена ровно в обеих колонках
    by_section: dict[str, set[str]] = {}
    for sec, col in headers:
        by_section.setdefault(sec, set()).add(col)
    for sec in SECTIONS:
        assert by_section.get(sec) == {"v2", "prod"}, (
            f"Секция '{sec}' должна иметь заголовки в обеих колонках, "
            f"найдено: {by_section.get(sec)}"
        )


def test_column_titles_v2_and_prod(mismatched_compare_html: str):
    """В первой строке grid — два column-title cell (V2 / Production)."""
    html = mismatched_compare_html
    assert re.search(
        r"<div\s+class='ps-col-title ps-cell-v2'\s+data-col='v2'>"
        r"V2 \(experimental\)</div>",
        html,
    ), "Ожидается V2 column title cell"
    assert re.search(
        r"<div\s+class='ps-col-title ps-cell-prod'\s+data-col='prod'>"
        r"Production</div>",
        html,
    ), "Ожидается Production column title cell"
