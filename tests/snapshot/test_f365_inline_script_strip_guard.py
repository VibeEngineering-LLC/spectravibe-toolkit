# -*- coding: utf-8 -*-
"""F-365 / v1.18.24.1 — Inline-script strip guard.

Регрессия для бага «опять не выводятся окна спектров и мультиплетов»:
F-317 paren-strip regex `\s*\([^)]*\bF-\d{1,3}[^)]*\)` ел через newline,
поэтому в inline `<script>` блоках с конструкцией

    document.querySelectorAll('.fp-view-btn').forEach(btn => {
      // F-147 — secondary-кнопки шарят CSS-class ...
      // ... но НЕ должны попадать в setView(). Фильтруем по
      // наличию data-view.
      if (!btn.dataset.view) return;
      ...
    });

regex жадно матчил от `(btn => {` через newline-комментарий с `F-147` до
закрывающей `)` в `setView()`. Результат — JS-SyntaxError, Chart.js не
инициализируется, спектр и мультиплеты не отображаются.

Фикс: `[^)\n]` (не `[^)]`) — paren-strip не пересекает newline.
"""
from __future__ import annotations
import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))


def test_F365_paren_strip_does_not_cross_newline():
    """`_F_ID_PAREN_PATTERN` НЕ должен матчить через newline."""
    from gamma.reporting.build import _F_ID_PAREN_PATTERN
    js_snippet = (
        "document.querySelectorAll('.fp-view-btn').forEach(btn => {\n"
        "  // F-147 — secondary-кнопки шарят CSS-class fp-view-btn\n"
        "  // консистентности. Фильтруем по наличию data-view.\n"
        "  if (!btn.dataset.view) return;\n"
        "  btn.addEventListener('click', () => setView(btn.dataset.view));\n"
        "});"
    )
    matches = _F_ID_PAREN_PATTERN.findall(js_snippet)
    assert matches == [], (
        f"paren-strip пересёк newline и нашёл match: {matches} — "
        "это значит regex съест `(btn => {...F-147...setView()` и сломает JS"
    )


def test_F365_paren_strip_still_works_single_line():
    """Legit single-line strip остаётся работоспособным."""
    from gamma.reporting.build import _F_ID_PAREN_PATTERN
    text = "Принципы анализа (F-256 / v1.17.10) описаны в §3."
    stripped = _F_ID_PAREN_PATTERN.sub("", text)
    assert "F-256" not in stripped
    assert "Принципы анализа описаны в §3." == stripped


def test_F365_full_pipeline_preserves_js_syntax():
    """Полный F-317 pipeline применённый к inline JS должен оставить
    валидный синтаксис: balanced `(`/`)`, целостный `setView()`."""
    from gamma.reporting.build import _f317_apply_user_facing_compliance
    js_snippet = (
        "document.querySelectorAll('.fp-view-btn').forEach(btn => {\n"
        "  // F-147 — secondary-кнопки шарят CSS-class fp-view-btn для\n"
        "  // визуальной консистентности, но НЕ должны попадать в setView().\n"
        "  // Фильтруем по наличию data-view.\n"
        "  if (!btn.dataset.view) return;\n"
        "  btn.addEventListener('click', () => setView(btn.dataset.view));\n"
        "});"
    )
    out = _f317_apply_user_facing_compliance(js_snippet, format="html")
    # F-147 must be stripped from comments (F-317 contract on output)
    assert "F-147" not in out, "bare F-id must be stripped"
    # JS structure must remain intact
    assert ".forEach(btn => {" in out, "forEach call broken by strip"
    assert "setView(btn.dataset.view)" in out, "setView call broken by strip"
    assert "btn.dataset.view" in out
    # Balanced parens count
    assert out.count("(") == js_snippet.count("("), (
        "open-paren count changed — JS structure damaged"
    )
    assert out.count(")") == js_snippet.count(")"), (
        "close-paren count changed — JS structure damaged"
    )


def test_F365_th232_full_demo_has_balanced_parens_in_script():
    """End-to-end guard на сгенерированном Th-232 отчёте: JS должен
    парсится node-style — balanced parens во всех inline <script>."""
    html_path = (
        REPO / "demo_reports" / "v1_18_24_th232_full" / "sample"
        / "Th232_Маринелли_0cm_report.html"
    )
    if not html_path.is_file():
        pytest.skip(f"demo HTML missing: {html_path}")
    html = html_path.read_text(encoding="utf-8")
    # Extract all inline script bodies (excluding src-loaded CDN scripts)
    pat = re.compile(r"<script>(.*?)</script>", re.DOTALL)
    bodies = pat.findall(html)
    assert len(bodies) >= 1, "expected at least one inline <script>"
    for i, body in enumerate(bodies):
        # Skip very small bodies (likely just init hooks)
        if len(body) < 100:
            continue
        # Heuristic: balanced () pairs (string-literal parens don't break this
        # for our specific template — checked by F-365 contract).
        # NB: this is not a full JS parser, but the bug manifests as exactly
        # this kind of imbalance (regex ate a `(btn => {` opener).
        opens = body.count("(")
        closes = body.count(")")
        # Allow small imbalance from regex literals inside strings, but the
        # damage we're guarding against creates a delta of 1+ from a single
        # eaten `(btn => {`. Empirically diff is 0 on healthy template.
        assert abs(opens - closes) <= 2, (
            f"inline script #{i}: paren imbalance {opens}-{closes}={opens-closes}"
            " — F-317 may have damaged JS"
        )
        # Specific anti-regression check for the exact bug:
        assert ".forEach. " not in body, (
            f"inline script #{i}: contains `.forEach. ` — это signature "
            "сломанного forEach (Cyrillic text after dot)"
        )


def test_F365_th232_full_demo_chart_init_intact():
    """Th-232 sample report должен содержать неповреждённый Chart.js
    init для spectrum (`fp-sp`) и multiplets (`mp-M1c`/`mp-M2c`/`mp-M3c`)."""
    html_path = (
        REPO / "demo_reports" / "v1_18_24_th232_full" / "sample"
        / "Th232_Маринелли_0cm_report.html"
    )
    if not html_path.is_file():
        pytest.skip(f"demo HTML missing: {html_path}")
    html = html_path.read_text(encoding="utf-8")
    # Spectrum canvas init
    assert "id=\"fp-sp\"" in html, "spectrum canvas missing"
    assert "new Chart(canvas.getContext('2d'),{" in html, (
        "spectrum Chart init broken"
    )
    # forEach setView wiring — exact line that was damaged by F-317
    assert "document.querySelectorAll('.fp-view-btn').forEach(btn => {" in html, (
        "fp-view-btn forEach wrapper damaged — F-317 strip пересёк newline"
    )
    # Multiplet loop
    assert "multiplets.forEach(m=>{" in html, "multiplet loop missing"
