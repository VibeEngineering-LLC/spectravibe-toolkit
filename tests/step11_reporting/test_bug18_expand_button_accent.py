"""BUG-18 / v1.18.30+ (Agent B) — «Развернуть» expand button must use accent color.

Symptom: кнопка "⛶ Развернуть" в header'е графика была outline-only
(`border:0.5px solid var(--border-secondary)`, background transparent),
из-за чего теряется среди других header-элементов. Primary CTA должна
визуально выделяться.

Fix: CSS обновлён на accent (#6566D7 light / #9D97E5 dark) с rgba 10%
fill, hover/focus заполняются полностью.

Тест проверяет HTML-template (статический шаблон, без рендеринга
analyze_and_report — быстрый contract check).
"""
from __future__ import annotations

import re
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = (REPO_ROOT / "scripts" / "gamma" / "reporting" / "templates"
            / "interactive_v1_17_2.html")


def _read_template() -> str:
    assert TEMPLATE.exists(), f"missing template: {TEMPLATE}"
    return TEMPLATE.read_text(encoding="utf-8")


def test_bug18_expand_button_has_accent_color():
    """`.fp-expand-btn` CSS должен содержать accent border (#6566D7)
    и accent color (#6566D7) — не outline-only secondary."""
    css = _read_template()
    # Match `.fp-expand-btn{ ... }` selector block (greedy until closing brace).
    m = re.search(r"\.fp-expand-btn\s*\{([^}]+)\}", css)
    assert m, "BUG-18: `.fp-expand-btn` CSS selector not found in template"
    block = m.group(1)
    # Accent color #6566D7 must appear at least twice (border + color).
    n_accent = block.count("#6566D7")
    assert n_accent >= 2, (
        "BUG-18: `.fp-expand-btn` must use accent color #6566D7 for both "
        "border AND text/color; found {} occurrences in block:\n{}"
        .format(n_accent, block)
    )
    # Background should be filled (not transparent) — rgba(101,102,215,0.10)
    assert "rgba(101,102,215" in block, (
        "BUG-18: `.fp-expand-btn` should have rgba accent fill, not "
        "transparent background; got block:\n{}".format(block)
    )


def test_bug18_expand_button_has_hover_invert():
    """Hover state must invert: background=accent, color=white."""
    css = _read_template()
    m = re.search(r"\.fp-expand-btn:hover\s*\{([^}]+)\}", css)
    assert m, "BUG-18: `.fp-expand-btn:hover` rule missing"
    block = m.group(1)
    assert "#6566D7" in block, (
        "BUG-18: hover state should fill with accent background; got:\n{}"
        .format(block)
    )
    assert "#ffffff" in block or "#fff" in block.lower(), (
        "BUG-18: hover state should set color to white for contrast; got:\n{}"
        .format(block)
    )


def test_bug18_dark_mode_variant_exists():
    """`@media (prefers-color-scheme: dark)` should override accent
    to a lighter shade (#9D97E5 / #CECBF6) for dark UI contrast."""
    css = _read_template()
    # Find the dark-mode block that mentions .fp-expand-btn
    m = re.search(
        r"@media\s*\(\s*prefers-color-scheme:\s*dark\s*\)\s*\{[^{}]*"
        r"\.fp-expand-btn\s*\{[^}]+\}",
        css,
    )
    assert m, (
        "BUG-18: dark-mode media-query override for `.fp-expand-btn` "
        "is missing — без него на dark UI accent #6566D7 теряет контраст."
    )


if __name__ == "__main__":
    test_bug18_expand_button_has_accent_color()
    test_bug18_expand_button_has_hover_invert()
    test_bug18_dark_mode_variant_exists()
    print("OK")
