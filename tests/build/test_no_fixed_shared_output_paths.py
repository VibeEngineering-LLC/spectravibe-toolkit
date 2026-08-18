"""P1-3b / P1-3c meta-guard — xdist shared-path race antipatterns.

Reds if any of the three fit-overlay / sortable-table test modules
reintroduces a fixed, on-disk output path written by a ``scope="module"``
(or ``scope="session"``) fixture — that's the **P1-3b** module-fixture
half of the race.

Reds if any ``tests/**/*.py`` reintroduces the broader **P1-3c** pattern:
**any** test function or fixture that hardcodes a scratch output dir of
the form ``demo_reports/_test_…`` or ``_tmp/…`` as a string assignment.
Under ``pytest -n auto`` (``--dist load``) test items distribute across
workers; two workers writing the **same** fixed dir create a write race:
worker B can read ``report.json`` while worker A is mid-write → partial
JSON missing top-level keys (e.g. ``primary_feps``) → intermittent
``KeyError`` / ``assert "primary_feps" in report`` failure. The empirical
trigger was 9 ``test_f396_fep_ru_translate.py`` test functions sharing
``demo_reports/_test_f396_fep_ru`` (full-suite ``-n auto ×3`` flaked on
run 2/3, 2026-06-05).

Fix uniformly: every test function gets ``tmp_path`` and writes
``out = str(tmp_path)``; every shared helper threads ``tmp_path`` /
``tmp_path_factory`` through; every module/session fixture uses
``tmp_path_factory.mktemp(...)``. pytest gives each worker a unique
``basetemp`` (``--basetemp=…/popen-gw0/...``) → no cross-worker collision.

This guard asserts the antipatterns are absent from source, so a
regression reds deterministically (a behavioral race repro is only
probabilistic).

Cite (file:line, pre-fix 369bd76 / ee7ffd1, 2026-06-05):

P1-3b (module-fixture, already fixed at ee7ffd1):
  - tests/step11_reporting/test_f_fit_view.py:51 — ``_OUT_DIR = "_tmp/test_f_fit_view"``.
  - tests/snapshot/test_f_fit_view_v2.py:45 — ``_OUT_DIR = "_tmp/test_f_fit_view_v2"``.
  - tests/snapshot/test_f_fit_view_v2.py:56 — ``cached_dir = "_tmp/test_f_fit_view/th232"`` (cross-file read).
  - tests/snapshot/test_f393_sortable_tables.py:43 — ``out = "demo_reports/_test_f393_sortable"``.

P1-3c (test-function / helper, fixed in this commit) — historical fixed dirs:
  - tests/step01_io_and_metadata/test_bg_only_environment.py:23 — ``demo_reports/_test_bg_only``.
  - tests/step11_reporting/test_anonymization.py:28 — ``demo_reports/_test_anonymize``.
  - tests/step11_reporting/test_bug14_bg_multiplet_section_hidden.py:48,98 — ``_test_bug14_bg_only``, ``_test_bug14_sample``.
  - tests/step11_reporting/test_interactive_report.py:27 — ``demo_reports/_test_interactive``.
  - tests/step11_reporting/test_bug43_bg_multiplet_511_fix.py:144 — ``_test_bug43_sample_multiplets``.
  - tests/step11_reporting/test_rows_sorted_ascending.py:20 — ``demo_reports/_test_sorted``.
  - tests/step11_reporting/test_no_en_leak.py:190 — ``demo_reports/_test_no_en``.
  - tests/step11_reporting/test_pdf_artefact.py:20 — ``demo_reports/_test_pdf_artefact``.
  - tests/step11_reporting/test_cost_footer.py:13 — ``demo_reports/_test_cost_footer``.
  - tests/step07_identification/test_chain_completeness.py:69 — ``_test_chain_completeness``.
  - tests/snapshot/test_f386_terminology_no_old_ru.py:65 — ``_test_f386_terminology``.
  - tests/step07_identification/test_th_composite_present.py:17 — ``_test_th_composite``.
  - tests/snapshot/test_f396_fep_ru_translate.py:41,203 — ``_test_f396_fep_ru`` (the empirical f396 flake).
"""
from __future__ import annotations

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

FIT_VIEW = REPO_ROOT / "tests" / "step11_reporting" / "test_f_fit_view.py"
FIT_VIEW_V2 = REPO_ROOT / "tests" / "snapshot" / "test_f_fit_view_v2.py"
F393 = REPO_ROOT / "tests" / "snapshot" / "test_f393_sortable_tables.py"

# ─────────────────────────────────────────────────────────────────────
# P1-3b (legacy, kept verbatim) — fixed-shared-write antipatterns for
# the three fit-overlay / sortable-table test modules. Each MUST be
# absent post-fix.
# ─────────────────────────────────────────────────────────────────────
_RE_OUT_DIR_TMP = re.compile(r"""_OUT_DIR\s*=\s*["']_tmp/""")
_RE_CACHED_DIR = re.compile(r"""cached_dir\s*=\s*["']_tmp/test_f_fit_view""")
_RE_DEMO_FIXED = re.compile(r"""["']demo_reports/_test_f393_sortable["']""")

# ─────────────────────────────────────────────────────────────────────
# P1-3c (broadened) — fixed scratch-dir WRITE assignment scanner that
# walks every ``tests/**/*.py``. Flags ANY of:
#   <indent> <var> = "demo_reports/_test…"
#   <indent> <var> = "_tmp/…"
# where <var> ∈ {out, output_dir, _OUT_DIR, cached_dir, cached}.
#
# Anchored to ^\s* so the regex only fires on assignment statements at
# line start, NOT on legitimate READS like:
#   ROOT / "demo_reports/v1_18_20/foo"  (no `<var> = ` LHS — not flagged)
#   pytest.skip("demo_reports/ не найден")  (string literal in call — no LHS)
#   os.environ["GAMMA_DEMO_REPORTS_DIR"]  (dict subscript — no LHS)
# All three patterns lack ``^\s*<varname>\s*=`` immediately before the
# scratch literal, so the predicate does not fire on them.
# ─────────────────────────────────────────────────────────────────────
_SCRATCH_RE = re.compile(
    r"""^\s*(out|output_dir|_OUT_DIR|cached_dir|cached)\s*=\s*["'](demo_reports/_test|_tmp/)"""
)

# Files this guard itself quotes the antipattern in (as test data /
# docstrings). They are EXEMPT from the scan; otherwise the scanner
# would always red on its own evidence strings. The exemption is a
# small whitelist of explicit absolute filenames — no glob.
_GUARD_FILE = Path(__file__).resolve()
_EXEMPT_FILES = {_GUARD_FILE}


def _iter_test_py_files():
    """Yield every ``tests/**/*.py`` Path under repo root, excluding
    the guard file itself (which quotes the historical antipattern as
    test data)."""
    tests_root = REPO_ROOT / "tests"
    for p in tests_root.rglob("*.py"):
        if p.resolve() in _EXEMPT_FILES:
            continue
        yield p


def test_fit_view_has_no_fixed_tmp_out_dir():
    """test_f_fit_view.py must not pin ``_OUT_DIR = "_tmp/..."``.

    A module-scoped fixture writing that fixed dir collides across xdist
    workers (P1-3b). Use ``tmp_path_factory.mktemp(...)`` instead.
    """
    assert FIT_VIEW.is_file(), f"missing test file: {FIT_VIEW}"
    src = FIT_VIEW.read_text(encoding="utf-8")
    assert not _RE_OUT_DIR_TMP.search(src), (
        "tests/step11_reporting/test_f_fit_view.py reintroduced a fixed "
        "`_OUT_DIR = \"_tmp/...\"` shared write path. Under `pytest -n auto` "
        "the module-scoped th232_result fixture writes this once per worker "
        "→ concurrent writes → reader hits a partial report.json (ERROR at "
        "fixture setup). Use tmp_path_factory.mktemp(...) (P1-3b)."
    )


def test_fit_view_v2_has_no_fixed_tmp_out_dir_or_cross_file_cache():
    """test_f_fit_view_v2.py must not pin ``_OUT_DIR`` nor read v1's fixed dir.

    The cross-file ``cached_dir = "_tmp/test_f_fit_view/th232"`` read is the
    P1-3b race amplifier: v2 (worker B) reads v1's report.json while v1
    (worker A) is still writing it.
    """
    assert FIT_VIEW_V2.is_file(), f"missing test file: {FIT_VIEW_V2}"
    src = FIT_VIEW_V2.read_text(encoding="utf-8")
    assert not _RE_OUT_DIR_TMP.search(src), (
        "tests/snapshot/test_f_fit_view_v2.py reintroduced a fixed "
        "`_OUT_DIR = \"_tmp/...\"` shared write path (P1-3b xdist race). "
        "Use tmp_path_factory.mktemp(...)."
    )
    assert not _RE_CACHED_DIR.search(src), (
        "tests/snapshot/test_f_fit_view_v2.py reintroduced the cross-file "
        "cache read `cached_dir = \"_tmp/test_f_fit_view/th232\"`. This let "
        "v2 read test_f_fit_view.py's report.json mid-write under "
        "`pytest -n auto` (P1-3b). v2 must run its own analyze_and_report "
        "into tmp_path_factory.mktemp(...)."
    )


def test_f393_sortable_has_no_fixed_demo_reports_out_dir():
    """test_f393_sortable_tables.py must not write fixed demo_reports/ dir.

    The ``interactive_html`` module-scoped fixture formerly wrote
    ``demo_reports/_test_f393_sortable`` — a fixed shared path (same xdist
    write race, plus it pollutes the repo's demo_reports/). Use
    tmp_path_factory.mktemp(...).
    """
    assert F393.is_file(), f"missing test file: {F393}"
    src = F393.read_text(encoding="utf-8")
    assert not _RE_DEMO_FIXED.search(src), (
        "tests/snapshot/test_f393_sortable_tables.py reintroduced the fixed "
        "write path `demo_reports/_test_f393_sortable` in the interactive_html "
        "module fixture (P1-3b xdist race + demo_reports/ pollution). Use "
        "tmp_path_factory.mktemp(...)."
    )


def test_no_fixed_scratch_dir_assignment_anywhere_under_tests():
    """P1-3c GENERAL scanner — no test under ``tests/**`` may assign a
    fixed scratch output dir.

    Catches both module/session fixtures (P1-3b regression class) AND
    per-test-function writes (P1-3c — the f396 flake class: 9 test funcs
    sharing ``demo_reports/_test_f396_fep_ru`` → ``primary_feps``
    KeyError on partial read).

    Anchored to ``^\\s*<var> = "scratch_literal"`` so legitimate READS
    of ``demo_reports/`` (path-construction via ``ROOT / "..."``, skip
    messages, ``os.environ`` lookups) are NOT flagged — only the
    assignment-to-scratch-string antipattern.
    """
    offenders: list[tuple[str, int, str]] = []
    for path in _iter_test_py_files():
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            # Skip unreadable files (binary garbage, encoding gremlins);
            # nothing we can scan, nothing to flag.
            continue
        for lineno, line in enumerate(text.splitlines(), start=1):
            if _SCRATCH_RE.search(line):
                offenders.append((str(path.relative_to(REPO_ROOT)), lineno, line.strip()))

    if offenders:
        msg_lines = [
            "P1-3c: the following test sources reintroduced a fixed scratch "
            "output dir assignment (`<var> = \"demo_reports/_test…\"` or "
            "`<var> = \"_tmp/…\"`). Under `pytest -n auto` two test items "
            "writing the SAME dir race on report.json → intermittent "
            "KeyError / partial-JSON failure. Replace with `tmp_path` "
            "fixture param (test functions) or `tmp_path_factory.mktemp(...)` "
            "(module/session fixtures). Offending lines:",
        ]
        for relpath, lineno, snippet in offenders:
            msg_lines.append(f"  {relpath}:{lineno}: {snippet}")
        raise AssertionError("\n".join(msg_lines))


def test_meta_guard_self_check_would_fire():
    """Self-check: the regex predicates actually detect the regressions.

    Confirms the assertions above are not vacuous — i.e. they WOULD red if
    the fixed-shared-write antipatterns were present. Exercised against
    captured pre-fix source snippets WITHOUT mutating any file
    (the 'reds without fix' half of the DoD, statically).

    Cite (pre-fix snippets captured 2026-06-05):
      - P1-3b @ 369bd76: _OUT_DIR v1/v2 + cached_dir + f393 demo_reports.
      - P1-3c @ ee7ffd1: 13 test files with ``out = "demo_reports/_test…"``
        or ``out = "_tmp/…"`` assignments (f396 was the empirical flake).
    """
    # (a) P1-3b legacy predicates fire on historical fixture lines.
    historical_out_dir_v1 = '_OUT_DIR = "_tmp/test_f_fit_view"'
    historical_out_dir_v2 = '_OUT_DIR = "_tmp/test_f_fit_view_v2"'
    assert _RE_OUT_DIR_TMP.search(historical_out_dir_v1)
    assert _RE_OUT_DIR_TMP.search(historical_out_dir_v2)

    # (b) cross-file cache predicate fires on the historical v2 line.
    historical_cached = '    cached_dir = "_tmp/test_f_fit_view/th232"'
    assert _RE_CACHED_DIR.search(historical_cached)

    # (c) demo_reports predicate fires on the historical f393 fixture line.
    historical_demo = '    out = "demo_reports/_test_f393_sortable"'
    assert _RE_DEMO_FIXED.search(historical_demo)

    # (d) P1-3c general scratch-assignment regex fires on the historical
    # f396 line — the empirical flake source (full-suite -n auto ×3 run 2
    # on 2026-06-05, ee7ffd1 main tree: missing `primary_feps`).
    historical_f396 = '    out = "demo_reports/_test_f396_fep_ru"'
    assert _SCRATCH_RE.search(historical_f396), (
        "self-check: P1-3c general scanner regex must fire on the captured "
        "f396 pre-fix line"
    )
    # And on the _tmp/ variant.
    historical_tmp = '    out = "_tmp/test_anything"'
    assert _SCRATCH_RE.search(historical_tmp)

    # (e) The regex must NOT fire on legitimate READ patterns — these are
    # the false-positive classes the docstring promises to skip. Each
    # appears verbatim in production tests (e.g. test_f390_compare_2col.py,
    # test_f393_sortable_tables.py, test_v1_18_18_5_kits_and_toggle.py).
    legitimate_reads = [
        '    p = ROOT / "demo_reports/v1_18_20/M_th_легкий_report.html"',
        '    pytest.skip("demo_reports/ не найден (set GAMMA_DEMO_REPORTS_DIR)")',
        '    base = os.environ["GAMMA_DEMO_REPORTS_DIR"]',
        # path joining (no `out =` LHS — different var name with /-join, or
        # f-string interpolation — not a literal scratch assignment).
        '    report_path = demo_root / "demo_reports/v1_18_20/foo.html"',
    ]
    for r in legitimate_reads:
        assert not _SCRATCH_RE.search(r), (
            f"self-check: P1-3c general scanner regex must NOT fire on "
            f"legitimate read pattern: {r!r}"
        )
