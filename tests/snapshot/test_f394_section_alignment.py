"""F-394 / v1.18.27 — section title horizontal alignment (CSS only).

Контракт: все h2/h3 в reporting секциях имеют согласованный CSS:
- `text-align: left` (explicit baseline, не наследуется случайно)
- `padding-left: 0` (consistent left margin)
- `vertical-align: middle` (для inline icons/badges)

Проверяем два места:
1. `scripts/gamma/reporting/templates/interactive_v1_17_2.html` — все
   `.fp-mp-block h3`, `.fp-notes h3`, `.gost-references h2`,
   `.passport-comparison h2` имеют unified text-align/padding-left.
2. `scripts/gen_v2_compare_th232.py` — base h2/h3 имеют unified rules
   и нет inconsistent inline `margin-top` overrides у h3-в-секции.
"""
from __future__ import annotations

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = REPO_ROOT / "scripts" / "gamma" / "reporting" / "templates" / "interactive_v1_17_2.html"
COMPARE_SCRIPT = REPO_ROOT / "scripts" / "gen_v2_compare_th232.py"


def _read(path: Path) -> str:
    assert path.exists(), f"missing: {path}"
    return path.read_text(encoding="utf-8")


def test_interactive_template_headers_have_unified_alignment():
    """Все четыре rule-set для h2/h3 в reporting секциях содержат
    `text-align:left` и `padding-left:0` и `vertical-align:middle`."""
    css = _read(TEMPLATE)

    required = ("text-align:left", "padding-left:0", "vertical-align:middle")
    selectors = (
        ".fp-mp-block h3",
        ".fp-notes h3",
        ".gost-references h2",
        ".passport-comparison h2",
    )

    for sel in selectors:
        # Найти строку правила (rule starts with selector, ends with `}`)
        m = re.search(re.escape(sel) + r"\{[^}]*\}", css)
        assert m, f"selector {sel!r} not found in template"
        rule = m.group(0)
        for prop in required:
            assert prop in rule, (
                f"selector {sel!r} missing required prop {prop!r}; rule={rule!r}"
            )


def test_v2_compare_h2_h3_have_unified_alignment():
    """В gen_v2_compare_th232.py base h2 и h3 CSS rules содержат
    text-align:left + padding-left:0 + vertical-align:middle."""
    src = _read(COMPARE_SCRIPT)

    # Rules are inside f-string with `{{` doubles. Extract from
    # "h2{{...}}" and "h3{{...}}" patterns.
    for tag in ("h2", "h3"):
        m = re.search(tag + r"\{\{[^}]*\}\}", src)
        assert m, f"base {tag} CSS rule not found in v2_compare script"
        rule = m.group(0)
        for prop in ("text-align:left", "padding-left:0", "vertical-align:middle"):
            assert prop in rule, (
                f"v2_compare base {tag} CSS missing {prop!r}; rule={rule!r}"
            )


def test_v2_compare_h3_no_inconsistent_margin_top_overrides():
    """Inline `style=\"margin-top:8px;\"` / `margin-top:18px;` на h3
    создавали visual imbalance — после F-394 убраны (только базовый
    margin из CSS rule остаётся)."""
    src = _read(COMPARE_SCRIPT)
    # На h3-тегах (открывающий) не должно остаться inline margin-top
    h3_with_margin = re.findall(r'<h3[^>]*style="[^"]*margin-top:[^"]*"', src)
    assert not h3_with_margin, (
        f"unexpected inline margin-top overrides on <h3>: {h3_with_margin}"
    )
