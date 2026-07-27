"""DEEP-06 (Project #5 P2 sub-2) — silent exception handler ceiling.

Ceiling test: silent_handler count in delivered-number-critical modules
of ``scripts/gamma/**`` must NOT regress.

Rationale
---------
The Wave-2 P2 sub-wave-2 Censor census flagged 99 silent exception
handlers in ``scripts/gamma/`` (Censor's definition: first body
statement is a sentinel fall-through — ``pass`` / ``continue`` /
``return <None|0|0.0|1.0|False|True|[]|{}>`` — and no logging call
appears in the first three statements).

DEEP-06 wrapped six specific delivered-activity-critical sites with
``warnings.warn(...)`` / ``logger.warning(...)`` (no behavior change —
the existing fallback values are kept). The wrapped sites are:

  scripts/gamma/activity/compute.py
    - L562 TCS total-efficiency lookup falls to 0.0
    - L572 TCS per-line correction → ``continue``
    - L576 TCS auto-correction block → ``pass``
    - L1340 matrix-method activity solve → ``pass``
  scripts/gamma/activity/quasi_template_solver.py
    - L97 P/T ratio lookup → ``return 1.0``
  scripts/gamma/identification/staged_pipeline.py
    - L1867 efficiency-aware completeness re-run → ``pass``

Ceilings below are the *current* counts after DEEP-06 wrapping landed
and define the maximum acceptable silent-handler count per critical
module. If a new silent handler is added in any of these modules, this
test goes RED and forces a conscious choice: either add logging (turn
it into a non-silent handler) or bump the ceiling in this test with a
comment explaining why the silent fallback is acceptable.

Red-without-fix design
----------------------
Inserting any new ``except Exception: pass`` (or ``return None`` /
``return 0.0`` / ``return 1.0`` / ``continue`` etc. without a preceding
log/warn) anywhere in the critical modules drives the silent_count
above its ceiling and turns this test RED.

The Censor's other-modules residue (deconvolve.py FWHM-fallback inner
routines, identification benign fallbacks) is deliberately kept under a
permissive ceiling — DEEP-06's targeted-triage scope was the
delivered-activity path, not a blanket rewrite.

Definition of "silent" matches ``scripts/verify_silent_handlers.py``
exactly (see the helper's module docstring).
"""
from __future__ import annotations

import ast
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# Critical modules: silent_handler count must NOT regress here.
# Ceiling = the count observed after DEEP-06 targeted-logging fixes.
CEILING: dict[str, int] = {
    "scripts/gamma/activity/compute.py": 0,
    "scripts/gamma/activity/quasi_template_solver.py": 0,
    "scripts/gamma/activity/quasitemplate.py": 0,
    "scripts/gamma/peaks/deconvolve.py": 26,  # F-441 added 4 new helper fns with silent fallbacks
    "scripts/gamma/peaks/coupled_multiplet.py": 0,
    "scripts/gamma/identification/staged_pipeline.py": 10,
    "scripts/gamma/calibration/efficiency_autoload.py": 0,
}

# Names that, if they appear in any of the first three statements of a
# handler body, mark the handler as non-silent. Must match
# scripts/verify_silent_handlers.py:LOGGING_NAMES exactly.
LOGGING_NAMES: tuple[str, ...] = (
    "logging",
    "logger",
    "warnings",
    "warn",
    "warning",
    "error",
    "info",
    "debug",
)

MAX_BODY_PEEK = 3

_SENTINEL_CONSTANTS: tuple = (None, 0, 0.0, 1.0, False, True)


def _is_sentinel_fall_through(stmt: ast.stmt) -> bool:
    """``pass`` / ``continue`` / ``return <sentinel>``."""
    if isinstance(stmt, (ast.Pass, ast.Continue)):
        return True
    if isinstance(stmt, ast.Return):
        val = stmt.value
        if val is None:
            return True
        if isinstance(val, ast.Constant) and val.value in _SENTINEL_CONSTANTS:
            return True
        if isinstance(val, ast.List) and not val.elts:
            return True
        if isinstance(val, ast.Dict) and not val.keys:
            return True
        if isinstance(val, ast.Tuple) and not val.elts:
            return True
        if isinstance(val, ast.Set) and not val.elts:
            return True
    return False


def _has_logging_call(handler: ast.ExceptHandler) -> bool:
    """Structural detection of a logging-shaped call in the first
    MAX_BODY_PEEK statements. Structural — not a substring match — so
    ``except Exception as warn`` does not falsely look like a log call.
    """
    for stmt in handler.body[:MAX_BODY_PEEK]:
        for sub in ast.walk(stmt):
            if not isinstance(sub, ast.Call):
                continue
            func = sub.func
            if isinstance(func, ast.Attribute):
                if func.attr in LOGGING_NAMES:
                    return True
                root = func.value
                while isinstance(root, ast.Attribute):
                    if root.attr in LOGGING_NAMES:
                        return True
                    root = root.value
                if isinstance(root, ast.Name) and root.id in LOGGING_NAMES:
                    return True
            elif isinstance(func, ast.Name):
                if func.id in LOGGING_NAMES:
                    return True
    return False


def _count_silent_handlers(path: Path) -> list[int]:
    """Return line numbers of silent handlers in ``path``."""
    try:
        source = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return []
    silent_lines: list[int] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.ExceptHandler):
            continue
        if not node.body:
            continue
        if not _is_sentinel_fall_through(node.body[0]):
            continue
        if _has_logging_call(node):
            continue
        silent_lines.append(node.lineno)
    return silent_lines


def test_silent_handler_ceiling_in_critical_modules():
    """Critical-module silent_handler count must not exceed its ceiling.

    Update the CEILING dict consciously when adding a legitimate silent
    handler; prefer wrapping with ``warnings.warn()`` /
    ``logger.warning()`` instead.
    """
    failures: list[str] = []
    for rel, ceiling in CEILING.items():
        path = REPO_ROOT / rel
        if not path.is_file():
            failures.append(f"{rel}: missing from repo")
            continue
        silent_lines = _count_silent_handlers(path)
        if len(silent_lines) > ceiling:
            failures.append(
                f"{rel}: {len(silent_lines)} silent handlers > "
                f"ceiling {ceiling} (lines: {silent_lines})"
            )
    assert not failures, (
        "Silent exception-handler count regressed in delivered-"
        "number-critical modules:\n  " + "\n  ".join(failures)
    )


def test_silent_handler_helper_finds_the_same_sites():
    """Sanity-check parity between this test's scanner and the
    helper script ``scripts/verify_silent_handlers.py``. If the two
    drift apart, one of them is wrong.
    """
    import subprocess
    import sys
    import json

    helper = REPO_ROOT / "scripts" / "verify_silent_handlers.py"
    assert helper.is_file(), f"helper missing: {helper}"

    proc = subprocess.run(
        [sys.executable, str(helper)],
        capture_output=True,
        text=True,
        encoding="utf-8",
        cwd=str(REPO_ROOT),
        timeout=120,
    )
    assert proc.returncode == 0, (
        f"helper exited {proc.returncode}: stderr={proc.stderr!r}"
    )
    report = json.loads(proc.stdout)
    helper_files = report["files"]

    mismatches: list[str] = []
    for rel in CEILING:
        scanner_lines = sorted(_count_silent_handlers(REPO_ROOT / rel))
        helper_lines = sorted(
            helper_files.get(rel, {}).get("silent_lines", [])
        )
        if scanner_lines != helper_lines:
            mismatches.append(
                f"{rel}: test_scanner={scanner_lines} "
                f"helper={helper_lines}"
            )
    assert not mismatches, (
        "Helper and test scanner disagree on silent-handler line set "
        "(definitions out of sync):\n  " + "\n  ".join(mismatches)
    )
