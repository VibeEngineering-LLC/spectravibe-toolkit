# -*- coding: utf-8 -*-
"""F-317 issue #36 — context-aware `_F_ID_PATTERN.sub` refactor tests.

Pre-refactor `scripts/gamma/reporting/build.py:121` did
`text = _F_ID_PATTERN.sub("", text)` blindly across the entire body.
This produced two failure modes:
  (a) legitimate prose «по правилу F-89.» → «по правилу .» (orphan
      punctuation residue). Evidence: BUG-6 workaround comments at
      interactive_html.py:2244-2252 and :2336-2341 — Agent B rewrote
      prose to *avoid* bare F-IDs to dodge this damage.
  (b) F-IDs intentionally embedded in `` `inline code` ``, ```fenced
      blocks```, `<code>`, `<pre>`, `<script>` (developer-facing
      surfaces — e.g. an example showing how F-317 itself is called)
      got nuked.

Fix: span-aware substitution (`_sub_outside_protected`) +
post-strip residue cleanup (`_cleanup_strip_residue`).

Contract from §1.3 brief: "only strip user-facing surfaces, preserve
developer-facing surfaces".
"""
from __future__ import annotations
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))


# ──────────────────────────────────────────────────────────────────
# Span protection — F-IDs inside dev-facing surfaces are preserved
# ──────────────────────────────────────────────────────────────────

def test_F317_36_backtick_inline_code_preserves_fid():
    """`F-317` inside markdown inline-code backticks must survive strip."""
    from gamma.reporting.build import _f317_apply_user_facing_compliance
    out = _f317_apply_user_facing_compliance(
        "Strip rule `F-317` applies to user-facing body.", format="md",
    )
    assert "`F-317`" in out, (
        f"F-317 inside backticks must be preserved, got: {out!r}"
    )


def test_F317_36_html_code_tag_preserves_fid():
    """F-id inside <code>…</code> survives strip (dev-facing surface)."""
    from gamma.reporting.build import _f317_apply_user_facing_compliance
    out = _f317_apply_user_facing_compliance(
        "<code>build.py implements F-317</code> strip.", format="html",
    )
    assert "<code>build.py implements F-317</code>" in out


def test_F317_36_html_pre_tag_preserves_fid():
    """F-id inside <pre>…</pre> survives strip."""
    from gamma.reporting.build import _f317_apply_user_facing_compliance
    src = "<pre>F-317 example\n  see comment</pre>"
    out = _f317_apply_user_facing_compliance(src, format="html")
    assert "<pre>F-317 example" in out


def test_F317_36_html_script_tag_preserves_fid_comments():
    """F-id inside <script>…</script> JS comments must survive (the
    F-365 fix prevents JS *syntax* damage; this fix prevents the F-IDs
    in `// F-NN —` comment-line markers from being silently stripped).
    """
    from gamma.reporting.build import _f317_apply_user_facing_compliance
    src = (
        "<script>\n"
        "// F-397 — bg block data\n"
        "// F-393 — sortable peak table\n"
        "doStuff();\n"
        "</script>"
    )
    out = _f317_apply_user_facing_compliance(src, format="html")
    assert "// F-397" in out
    assert "// F-393" in out
    assert "doStuff();" in out


def test_F317_36_fenced_codeblock_preserves_fid():
    """Triple-backtick fenced block — F-IDs preserved."""
    from gamma.reporting.build import _f317_apply_user_facing_compliance
    src = (
        "Header\n\n"
        "```python\n"
        "# F-317 example\n"
        "_f317_apply_user_facing_compliance(text)\n"
        "```\n\n"
        "Tail."
    )
    out = _f317_apply_user_facing_compliance(src, format="md")
    assert "# F-317 example" in out, "fenced code-block F-id stripped"


# ──────────────────────────────────────────────────────────────────
# Residue cleanup — prose strip leaves no orphan punctuation
# (closes BUG-6 root cause; workarounds at interactive_html.py
#  :2244-2252 and :2336-2341 are no longer required)
# ──────────────────────────────────────────────────────────────────

def test_F317_36_prose_strip_cleans_orphan_period():
    """«Метод σ по правилу F-89.» → «Метод σ по правилу.» (no orphan)."""
    from gamma.reporting.build import _f317_apply_user_facing_compliance
    out = _f317_apply_user_facing_compliance(
        "Метод σ по правилу F-89.", format="md",
    )
    assert out == "Метод σ по правилу.", f"got: {out!r}"


def test_F317_36_heading_strip_cleans_empty_prefix():
    """«### F-145: title» → «### title» (no empty `### :` residue)."""
    from gamma.reporting.build import _f317_apply_user_facing_compliance
    out = _f317_apply_user_facing_compliance(
        "### F-145: двухфазная самокалибровка", format="md",
    )
    assert out == "### двухфазная самокалибровка", f"got: {out!r}"


def test_F317_36_paired_strip_cleans_orphan_slash():
    """«see F-317/F-365 for context» → «see for context» (no orphan «/»)."""
    from gamma.reporting.build import _f317_apply_user_facing_compliance
    out = _f317_apply_user_facing_compliance(
        "see F-317/F-365 for context", format="md",
    )
    assert out == "see for context", f"got: {out!r}"


# ──────────────────────────────────────────────────────────────────
# Regression baseline — pre-refactor good behaviour still holds
# ──────────────────────────────────────────────────────────────────

def test_F317_36_parenthetical_strip_still_works():
    """Parenthetical F-id strip — legacy v1.18.15 contract still holds."""
    from gamma.reporting.build import _f317_apply_user_facing_compliance
    out = _f317_apply_user_facing_compliance(
        "Принципы анализа (F-256 / v1.17.10) описаны в §3.", format="md",
    )
    assert out == "Принципы анализа описаны в §3.", f"got: {out!r}"


def test_F317_36_parenthetical_strip_does_not_cross_newline():
    """F-365 fix — paren-strip class `[^)\\n]` does NOT match through
    newlines. Re-asserted here so the issue-#36 refactor cannot regress
    F-365 even by accident."""
    from gamma.reporting.build import _F_ID_PAREN_PATTERN
    js = (
        "document.querySelectorAll('.fp-view-btn').forEach(btn => {\n"
        "  // F-147 — secondary buttons\n"
        "  setView(btn.dataset.view);\n"
        "});"
    )
    matches = _F_ID_PAREN_PATTERN.findall(js)
    assert matches == [], f"paren-strip crossed newline: {matches}"


def test_F317_36_external_api_symbols_preserved():
    """Refactor must not rename public-ish helpers — F-365 tests import
    `_F_ID_PAREN_PATTERN` and `_f317_apply_user_facing_compliance` by
    name."""
    from gamma.reporting import build as B
    assert hasattr(B, "_F_ID_PATTERN")
    assert hasattr(B, "_F_ID_PAREN_PATTERN")
    assert hasattr(B, "_KT_ID_PAREN_PATTERN")
    assert hasattr(B, "_f317_apply_user_facing_compliance")


# ──────────────────────────────────────────────────────────────────
# Edge cases — span tokeniser robustness
# ──────────────────────────────────────────────────────────────────

def test_F317_36_unclosed_code_span_protects_to_eot():
    """Defensive: unclosed `` ` `` opener protects to end-of-text rather
    than crashing or stripping arbitrarily."""
    from gamma.reporting.build import _f317_apply_user_facing_compliance
    out = _f317_apply_user_facing_compliance(
        "Mention `F-317 unclosed", format="md",
    )
    assert "F-317" in out  # protected by unclosed-opener fallback


def test_F317_36_mixed_protected_and_prose():
    """Mixed input — F-id in prose stripped, F-id in backticks preserved,
    both in same string."""
    from gamma.reporting.build import _f317_apply_user_facing_compliance
    src = "Принципы (F-256) описаны в §3. Implementation: `F-317` helper."
    out = _f317_apply_user_facing_compliance(src, format="md")
    assert "F-256" not in out  # paren-strip
    assert "`F-317`" in out    # backtick-protected
