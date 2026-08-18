"""SEC-02 — IAEA fetcher must reject non-whitelisted nuclide names.

`scripts/gamma/data/iaea_fetcher._normalize_nuclide_name` at parent SHA
e2bddbb falls through with `return s` on regex non-match (line 78). `s`
is just the lowercased input with hyphens stripped — NO whitelist, NO
sanitization. The unsanitized string is then interpolated into the
IAEA API URL at `_build_url` (line 103):

    f"{IAEA_API_URL}?fields=decay_rads&nuclides={nuc_norm}&rad_types=g"

An adversary supplying e.g. `"Cs-137?evil=1"` lowercases to
`"cs137?evil=1"`; neither regex matches; the function returns
`"cs137?evil=1"`; `_build_url` then interpolates an attacker-
controlled `?evil=` query parameter into the URL. Worse, a name like
`"../etc/passwd"` (lowercases to `"../etc/passwd"`, hyphens stripped)
could escape the cache directory under `_cache_path`.

Fix per censor envelope 7827a815 (DECISION b): at `_normalize_nuclide_name`
line 78, replace `return s` with `raise ValueError(...)`. Fail-loud,
no sanitised fallback. Both upstream callers `_cache_path` and
`_build_url` then fail-loud on adversary input instead of silently
passing unsanitised characters downstream.

Tests:

  1. `test_valid_nuclide_names_accepted` — known-good names (Cs-137,
     Tc-99m, Th-232, Pb-212, K-40) must pass through and return
     canonical IAEA form (137cs, 99mtc, 232th, 212pb, 40k).

  2. `test_invalid_chars_raise` — URL-injection characters and path
     traversal (`Cs-137?evil=1`, `Cs/137`, `Cs%2D137`, `../etc/passwd`)
     must raise `ValueError`.

  3. `test_empty_and_oversized_raise` — empty string and 100-char
     oversized nuclide names must raise `ValueError`. The 100-char
     bound is a reasonable upper bound (longest real nuclide name
     is e.g. `Md-258m`, ~7 chars; mass-symbol ≤ 7 chars realistically).

  4. `test_url_does_not_contain_injection` — defence-in-depth: calling
     `_build_url("Cs-137?evil=1")` must raise (not return a URL with
     `?evil=` in the nuclide segment).

Red-without-fix evidence: `_tmp/red_sec02_parent_20260606.txt`
captures parent SHA e2bddbb pytest output where invalid chars
silently fall through (return `s` instead of raising).
Post-fix green: `_tmp/green_sec02_post_20260606.txt`.
"""

from __future__ import annotations

import pytest


def test_valid_nuclide_names_accepted():
    """Known-good nuclide names normalise to canonical IAEA form."""
    from gamma.data.iaea_fetcher import _normalize_nuclide_name as norm

    # Mass-number-first canonical IAEA form (per module docstring:
    # "234TH" → "234th", "Th-234" → "234th").
    assert norm("Cs-137") == "137cs"
    assert norm("cs137") == "137cs"
    assert norm("137CS") == "137cs"
    assert norm("Th-232") == "232th"
    assert norm("Pb-212") == "212pb"
    assert norm("K-40") == "40k"

    # Metastable suffix preserved. The current regex
    # `^([a-z]+)(\d+)([a-z]*)$` against the hyphen-stripped lowercase
    # "tc99m" captures sym="tc", num="99", meta="m" and joins as
    # `{num}{sym}{meta}` → "99tcm". That is the canonical form this
    # codebase uses; SEC-02 does not change normalisation, only adds
    # fail-loud on non-match.
    assert norm("Tc-99m") == "99tcm"


def test_invalid_chars_raise():
    """URL-injection chars and path traversal must raise ValueError."""
    from gamma.data.iaea_fetcher import _normalize_nuclide_name as norm

    # URL-parameter injection.
    with pytest.raises(ValueError, match=r"[Ii]nvalid nuclide"):
        norm("Cs-137?evil=1")

    # Embedded path separator — could otherwise corrupt `_cache_path`.
    with pytest.raises(ValueError, match=r"[Ii]nvalid nuclide"):
        norm("Cs/137")

    # Percent-encoded hyphen — passes `replace('-', '')` no-op then fails
    # regex; pre-fix would silently return "cs%2d137".
    with pytest.raises(ValueError, match=r"[Ii]nvalid nuclide"):
        norm("Cs%2D137")

    # Path traversal.
    with pytest.raises(ValueError, match=r"[Ii]nvalid nuclide"):
        norm("../etc/passwd")

    # Pure letters (no number) — not a valid nuclide spec.
    with pytest.raises(ValueError, match=r"[Ii]nvalid nuclide"):
        norm("Cesium")

    # Pure number (no symbol) — not a valid nuclide spec.
    with pytest.raises(ValueError, match=r"[Ii]nvalid nuclide"):
        norm("137")


def test_empty_and_oversized_raise():
    """Empty and absurdly long names must raise ValueError."""
    from gamma.data.iaea_fetcher import _normalize_nuclide_name as norm

    with pytest.raises(ValueError, match=r"[Ii]nvalid nuclide"):
        norm("")

    with pytest.raises(ValueError, match=r"[Ii]nvalid nuclide"):
        norm("    ")  # whitespace-only after strip

    # Oversized input — far beyond any real nuclide name length.
    with pytest.raises(ValueError, match=r"[Ii]nvalid nuclide"):
        norm("Cs-" + "1" * 100)


def test_url_does_not_contain_injection():
    """Defence-in-depth: `_build_url` must raise on adversary input.

    `_build_url("Cs-137?evil=1")` must NOT return a URL containing
    `?evil=` in the nuclide segment; it must propagate the underlying
    `_normalize_nuclide_name` ValueError.
    """
    from gamma.data.iaea_fetcher import _build_url

    with pytest.raises(ValueError, match=r"[Ii]nvalid nuclide"):
        _build_url("Cs-137?evil=1")

    with pytest.raises(ValueError, match=r"[Ii]nvalid nuclide"):
        _build_url("../etc/passwd")
