# -*- coding: utf-8 -*-
"""F-371 / v1.18.24.6 — setView() пересчитывает Y-axis range при смене view.

Регрессия: пользователь сообщил «При первом нажатии на кнопку "фон" спектр
фона не верно масштабируется. При дальнейших переключениях ошибка пропадает»
(2026-06-01).

Root cause: `setView(view)` в `interactive_v1_17_2.html` обновлял datasets,
annotations, legend — но НЕ пересчитывал `chart.options.scales.y.min/max`.
Y-bounds сохранялись от sample-view (cps logMax=5e2, counts logMax=2e5 или
auto-fit от sample-max). Когда первый клик переключает на bg (масштабированный
к k=t_sample/t_bg ≈ 0.24), bg-данные оказываются «прижаты к низу» log-шкалы
до самого верхнего sample-bound. Любое второе действие (units toggle,
log/linear, второй клик view) triggered setYScale → bounds пересчитывались.

Fix: в конце `setView()` добавлен явный recompute Y-min/max по той же
canonical логике что setYScale — log → yAxisDefaults().logMin/logMax,
linear → Math.max(maxDataPoint(), threshold)*1.1.
"""
from __future__ import annotations
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]


def _read_template():
    tmpl = REPO / "scripts" / "gamma" / "reporting" / "templates" / "interactive_v1_17_2.html"
    if not tmpl.is_file():
        pytest.skip(f"template missing: {tmpl}")
    return tmpl.read_text(encoding="utf-8")


def test_F371_setview_contains_yaxis_recompute_block():
    """Template должен содержать F-371 recompute block в setView()."""
    tmpl = _read_template()
    assert "F-371" in tmpl, "F-371 marker не найден в template"
    # Recompute должен ставить Y min/max в зависимости от scale type
    assert "scales.y.type === 'logarithmic' ? 'log' : 'linear'" in tmpl, (
        "F-371 fix должен определять текущий Y-scale type"
    )
    assert "chart.options.scales.y.min = _yd.logMin" in tmpl, (
        "F-371: log mode → logMin"
    )
    assert "chart.options.scales.y.max = _yd.logMax" in tmpl, (
        "F-371: log mode → logMax"
    )
    assert "maxDataPoint()" in tmpl, "F-371: linear mode → maxDataPoint() recompute"


def test_F371_recompute_inside_setview_body():
    """Recompute должен находиться ВНУТРИ setView() — иначе он не
    срабатывает при view-switch. Проверяем что блок появляется между
    `setView = function(view){` и закрывающей `};`."""
    tmpl = _read_template()
    # Match function body
    m = re.search(
        r"setView\s*=\s*function\s*\(\s*view\s*\)\s*\{(.*?)\n\s{4}\};",
        tmpl, re.DOTALL,
    )
    assert m, "не найден setView function body в template"
    body = m.group(1)
    assert "F-371" in body, (
        "F-371 recompute должен быть ВНУТРИ setView body, не снаружи"
    )
    # Recompute должен идти ДО chart.update — иначе update берёт старые bounds
    f371_idx = body.find("F-371")
    update_idx = body.rfind("chart.update")
    assert f371_idx < update_idx, (
        "F-371 recompute должен предшествовать chart.update('none')"
    )


def test_F371_demo_html_has_recompute_in_setview():
    """End-to-end: сгенерированный HTML demo содержит F-371 recompute."""
    demo = (
        REPO / "demo_reports" / "v1_18_24_th232_full" / "sample"
        / "Th232_Маринелли_0cm_report.html"
    )
    if not demo.is_file():
        pytest.skip(f"demo HTML missing: {demo}")
    html = demo.read_text(encoding="utf-8")
    # Ищем setView function body
    m = re.search(
        r"setView\s*=\s*function\s*\(\s*view\s*\)\s*\{(.*?)\n\s{4}\};",
        html, re.DOTALL,
    )
    if not m:
        pytest.skip("setView function not found in demo (no bg)")
    body = m.group(1)
    assert "F-371" in body, (
        "F-371 recompute marker отсутствует в demo HTML setView body"
    )
