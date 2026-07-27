"""F-393 / v1.18.27 — sortable peak tables в interactive HTML + compare HTML.

Контракт:
  1) Заголовки колонок таблиц помечены `data-sortable="true"` (либо
     `onclick="sortTable"`) — кликабельны для сортировки.
  2) E (энергия) и Nuclide (изотоп) sortable.
  3) Default order — по возрастанию энергии (rows pre-sorted Python-кодом,
     визуальный маркер `.sorted.asc` ставится JS на колонке E).
  4) JS sort utility (vanilla, без библиотек) присутствует в HTML.

Проверяемые отчёты:
  • interactive HTML (scripts/gamma/reporting/interactive_html.py)
  • compare HTML  (scripts/gen_v2_compare_th232.py) — best-effort, skip
    если demo-snapshot отсутствует.
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

sys.path.insert(0, "scripts")

from gamma.reporting import analyze_and_report  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "gen_v2_compare_th232.py"


# ──────────────────────────────────────────────────────────────────
# (1) interactive HTML — основной peak-таблица
# ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def interactive_html(tmp_path_factory) -> str:
    # Per-worker unique out-dir: under ``pytest -n auto`` a fixed shared
    # path (formerly ``demo_reports/_test_f393_sortable``) was written
    # once per worker per module → concurrent writers / partial reads
    # (P1-3b xdist race). mktemp gives each worker a unique dir and also
    # stops polluting ``demo_reports/``.
    out = str(tmp_path_factory.mktemp("f393_sortable"))
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
        write_markdown=False,
        write_plots=False,
        sample_mass_kg=0.5,
        background_path=bg,
    )
    return Path(res["html"]).read_text(encoding="utf-8")


def test_interactive_th_data_sortable_attr(interactive_html: str):
    """F-393 — каждое th в peak-table имеет data-sortable="true"."""
    html = interactive_html
    # Грубая выборка: th с data-col="iso"/"line"/"a" должен быть sortable
    for col in ("iso", "line", "a"):
        # Pattern: <th data-col="<col>" ... data-sortable="true">
        pat = (r'<th[^>]*data-col="' + re.escape(col)
               + r'"[^>]*data-sortable="true"')
        assert re.search(pat, html), (
            f"<th data-col=\"{col}\"> must carry data-sortable=\"true\""
        )


def test_interactive_th_sort_type_hints(interactive_html: str):
    """F-393 — числовые колонки помечены data-sort="num",
    строковые — data-sort="str"."""
    html = interactive_html
    # E (line) — num
    assert re.search(
        r'<th[^>]*data-col="line"[^>]*data-sort="num"', html
    ), "E column (data-col=\"line\") must have data-sort=\"num\""
    # A — num
    assert re.search(
        r'<th[^>]*data-col="a"[^>]*data-sort="num"', html
    ), "A column (data-col=\"a\") must have data-sort=\"num\""
    # Iso — str
    assert re.search(
        r'<th[^>]*data-col="iso"[^>]*data-sort="str"', html
    ), "Iso column (data-col=\"iso\") must have data-sort=\"str\""


def test_interactive_initial_rows_sorted_by_energy_asc(interactive_html: str):
    """F-393 — initial rows ordered by E asc (regex extract row data)."""
    html = interactive_html
    m = re.search(r"const\s+rows\s*=\s*(\[)", html)
    assert m, "could not locate `const rows = [...]` in HTML"
    start = m.end() - 1
    depth, i = 0, start
    while i < len(html):
        c = html[i]
        if c == "[":
            depth += 1
        elif c == "]":
            depth -= 1
            if depth == 0:
                break
        i += 1
    rows = json.loads(html[start:i + 1])
    assert rows, "rows array is empty"

    energies = []
    for r in rows:
        line = r.get("line") or ""
        em = re.search(r"\d+(?:\.\d+)?", line)
        if em:
            energies.append(float(em.group(0)))
    assert len(energies) >= 5, f"too few energies parsed: {energies}"

    for j in range(1, len(energies)):
        assert energies[j] >= energies[j - 1] - 1e-6, (
            f"rows not sorted ascending at i={j}: "
            f"{energies[j - 1]} > {energies[j]}"
        )


def test_interactive_sort_js_handler_present(interactive_html: str):
    """F-393 — JS-функция сортировки (vanilla) внедрена в template."""
    html = interactive_html
    # Pattern: применяем sort к rows массиву + classList.add('sorted')
    assert "applySort" in html or "sortTable" in html, (
        "JS sort handler (applySort or sortTable) must be inlined"
    )
    # Visual indicator class
    assert ".sorted" in html, ".sorted CSS class missing for sort indicator"
    # Direction toggle markers
    assert "asc" in html and "sortState" in html, (
        "asc-toggle logic (sortState or data-sort-dir) must be present"
    )


def test_interactive_default_sort_marker(interactive_html: str):
    """F-393 — default visual marker (.sorted.asc) ставится на колонке line (E)."""
    html = interactive_html
    # JS-код должен явно помечать колонку "line" (E) как sorted asc на init.
    # Допускаем разные эквивалентные формулировки:
    #   • `x.dataset.col === 'line'` (или `==`) — runtime equality check
    #   • `data-default-sort="asc"` — declarative HTML hint
    has_marker = (
        re.search(r"\.dataset\.col\s*===?\s*['\"]line['\"]", html)
        or "data-default-sort" in html
        or re.search(r"sortState\s*=\s*\{\s*col\s*:\s*['\"]line['\"]", html)
    )
    assert has_marker, (
        "JS must apply default .sorted.asc to E column on init "
        "(via dataset.col==='line' check, sortState init, or data-default-sort attr)"
    )


# ──────────────────────────────────────────────────────────────────
# (2) compare HTML — best-effort smoke test (skip if demo missing)
# ──────────────────────────────────────────────────────────────────

def _resolve_demo_root() -> Path | None:
    env = os.environ.get("GAMMA_DEMO_REPORTS_DIR")
    if env:
        p = Path(env)
        if p.is_dir():
            return p
    sibling = REPO_ROOT.parent / "demo_reports"
    if sibling.is_dir():
        return sibling
    return None


def _find_demo_run(demo_root: Path) -> Path | None:
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
def compare_html(tmp_path_factory) -> str:
    demo_root = _resolve_demo_root()
    if demo_root is None:
        pytest.skip("demo_reports/ не найден (set GAMMA_DEMO_REPORTS_DIR)")
    src_run = _find_demo_run(demo_root)
    if src_run is None:
        pytest.skip(f"в {demo_root} нет run с sample/+sample_v2/*_report.json")

    tmp_root = tmp_path_factory.mktemp("f393_compare")
    dst_run = tmp_root / src_run.name
    dst_run.mkdir()
    for sub in ("sample", "sample_v2"):
        shutil.copytree(src_run / sub, dst_run / sub)
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), str(dst_run)],
        capture_output=True, text=True, env=env, encoding="utf-8",
    )
    assert proc.returncode == 0, (
        f"gen_v2_compare_th232.py failed:\n"
        f"stdout={proc.stdout}\nstderr={proc.stderr}"
    )
    out = dst_run / "v2_compare" / "v2_compare_report.html"
    return out.read_text(encoding="utf-8")


def test_compare_tables_have_sortable_attr(compare_html: str):
    """F-393 — все peak tables (primary, secondary, diff) — data-sortable."""
    html = compare_html
    # Минимум 4 sortable таблицы: primary V2, primary prod, matched, only_v2
    n = html.count('data-sortable="true"')
    assert n >= 4, (
        f"compare HTML должно содержать ≥4 sortable таблиц, найдено: {n}"
    )


def test_compare_sortable_th_onclick(compare_html: str):
    """F-393 — th в peak tables имеют onclick="sortTable(this)"."""
    assert "onclick=\"sortTable(this)\"" in compare_html, (
        "compare HTML must have onclick=\"sortTable(this)\" on sortable th"
    )


def test_compare_sort_js_util_present(compare_html: str):
    """F-393 — vanilla JS sortTable() utility внедрена."""
    html = compare_html
    assert "function sortTable" in html, (
        "sortTable() function must be defined inline in compare HTML"
    )
    # Pattern из task spec: rows.sort + cells[colIdx] + parseFloat / localeCompare
    assert "cells[colIdx]" in html, (
        "sort util must index by colIdx into row.cells[]"
    )


def test_compare_default_sort_on_E_column(compare_html: str):
    """F-393 — колонка E помечена data-default-sort="asc"."""
    assert 'data-default-sort="asc"' in compare_html, (
        "хотя бы один th должен нести data-default-sort=\"asc\" "
        "(маркер default sort by E asc)"
    )


def test_compare_primary_peaks_sorted_by_E_asc(compare_html: str):
    """F-393 — initial primary peaks rows ordered by E asc в обеих колонках."""
    html = compare_html
    # Найти все <table class="peaks"> (primary peaks tables) и проверить
    # порядок E_keV (2nd numeric cell). Парсим только первые 2 (V2 + prod
    # primary peaks); diff-tables идут позже и имеют другой контракт.
    tables = re.findall(
        r"<table[^>]*class=['\"]peaks['\"][^>]*>(.*?)</table>",
        html, flags=re.DOTALL,
    )
    if len(tables) < 2:
        pytest.skip("less than 2 peak tables in compare HTML")
    n_checked = 0
    for tbl in tables[:2]:
        rows_html = re.findall(r"<tr>(.*?)</tr>", tbl, flags=re.DOTALL)
        # skip header row
        if not rows_html or "<th" in rows_html[0]:
            rows_html = rows_html[1:]
        energies = []
        for r in rows_html:
            cells = re.findall(r"<td[^>]*>(.*?)</td>", r, flags=re.DOTALL)
            if len(cells) >= 2:
                m = re.search(r"\d+(?:\.\d+)?", cells[1])
                if m:
                    energies.append(float(m.group(0)))
        if len(energies) >= 3:
            n_checked += 1
            for j in range(1, len(energies)):
                assert energies[j] >= energies[j - 1] - 1e-6, (
                    f"primary peaks table not sorted by E asc: "
                    f"{energies[j - 1]} > {energies[j]}"
                )
    assert n_checked >= 1, (
        "no primary peaks table with >=3 rows found for sort check"
    )


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(pytest.main([__file__, "-v"]))
