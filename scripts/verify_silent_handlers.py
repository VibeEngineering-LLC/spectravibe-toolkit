from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""Scan scripts/gamma/**/*.py for silent exception handlers.

DEEP-06 (Project #5 critique remediation, Wave 2 P2 sub-wave 2):
maintainability / observability survey of exception handlers in the
delivered-number-critical code paths.

Definition of "silent" used here (matches the Censor's S5-deep-P2P3
definition, see brief DEEP-06 §"Count definition"):

  An ``ast.ExceptHandler`` is silent if **both** hold:

  1. The **first** body statement is one of the sentinel-fall-through
     shapes — ``pass``, ``continue``, or ``return <sentinel>`` where
     ``<sentinel>`` is ``None``, ``0``, ``0.0``, ``1.0``, ``[]``,
     ``{}``, ``False``, or ``True``.

  2. No logging-shaped call (``logging.*``, ``logger.*``,
     ``warnings.warn``, ``.warning(`` / ``.error(`` / ``.info(`` /
     ``.debug(``) appears in the first ``MAX_BODY_PEEK`` statements.

That is the Censor's 99-silent definition: handlers that visibly fall
through to a sentinel without raising or logging.  Handlers that
``warnings.warn(...)`` and then continue to a sentinel are NOT counted
as silent — adding a warn turns a silent site into an observable one.

Handlers that re-raise are not silent either (rule covered by 1: a
``raise`` statement is not a sentinel fall-through).

Output: JSON document on stdout.  Schema::

    {
      "total_excepts": int,
      "silent_count": int,
      "files": {
        "<relative path>": {
          "total_excepts": int,
          "silent_count": int,
          "silent_lines": [int, ...]
        },
        ...
      },
      "sites_critical": [
        {"file": "...", "line": int, "module": "compute.py"}, ...
      ],
      "sites_other": [
        {"file": "...", "line": int}, ...
      ]
    }

Exit code:
  0 — normal, even when silent handlers are present.
  1 — only when ``--strict`` is given and ``silent_count`` exceeds
      ``--ceiling`` (integer, default 9999).

The script never writes files and never modifies the tree.

F-155: lives in ``scripts/`` (Python may live only there).
F-115: emits relative paths only — never absolute operator paths.
"""


import argparse
import ast
import json
import sys
from pathlib import Path
from typing import Iterable

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Names that, if they appear in ``ast.dump(stmt)``, mark the handler as
# non-silent.  We accept both attribute-call style (``logger.warning(``,
# ``logging.error(``) and the stdlib ``warnings.warn`` flow.
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

# Inspect at most this many statements at the top of the handler body
# when classifying.  Three is enough to cover the common pattern
# ``log()  # then  sentinel = ...; return sentinel``.
MAX_BODY_PEEK = 3

# Modules whose silent handlers affect delivered activity numbers.
# These are the "critical" sites for DEEP-06 triage.
CRITICAL_MODULES: frozenset[str] = frozenset({
    "compute.py",
    "quasi_template_solver.py",
    "quasitemplate.py",
    "deconvolve.py",
    "coupled_multiplet.py",
    "staged_pipeline.py",
    "efficiency_autoload.py",
})


# ---------------------------------------------------------------------------
# AST helpers
# ---------------------------------------------------------------------------


_SENTINEL_CONSTANTS: tuple = (None, 0, 0.0, 1.0, False, True)


def _is_sentinel_fall_through(stmt: ast.stmt) -> bool:
    """Return True if ``stmt`` is a pass / continue / return <sentinel>."""
    if isinstance(stmt, (ast.Pass, ast.Continue)):
        return True
    if isinstance(stmt, ast.Return):
        val = stmt.value
        if val is None:
            return True  # bare `return`
        if isinstance(val, ast.Constant) and val.value in _SENTINEL_CONSTANTS:
            return True
        # Empty collection literals.
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
    """Return True if any of the first MAX_BODY_PEEK statements is a
    logging-shaped call (``logger.warning(...)``, ``warnings.warn(...)``,
    ``logging.error(...)``, ...).

    Uses a structural walk rather than a substring match on
    ``ast.dump`` so the exception variable name itself never produces
    a false positive (e.g. ``except Exception as warn`` would otherwise
    match the substring ``warn``).
    """
    for stmt in handler.body[:MAX_BODY_PEEK]:
        for sub in ast.walk(stmt):
            if not isinstance(sub, ast.Call):
                continue
            func = sub.func
            # Attribute call: foo.warning(...), logging.error(...),
            # warnings.warn(...).
            if isinstance(func, ast.Attribute):
                if func.attr in LOGGING_NAMES:
                    return True
                # x.logger.warning(...) — match the receiver chain too.
                root = func.value
                while isinstance(root, ast.Attribute):
                    if root.attr in LOGGING_NAMES:
                        return True
                    root = root.value
                if isinstance(root, ast.Name) and root.id in LOGGING_NAMES:
                    return True
            elif isinstance(func, ast.Name):
                # Bare names: warn(...), warning(...).
                if func.id in LOGGING_NAMES:
                    return True
    return False


def _is_silent(handler: ast.ExceptHandler) -> bool:
    """Return True if the handler swallows the exception silently.

    Censor-style rule: first body statement is a sentinel fall-through
    (``pass`` / ``continue`` / ``return <sentinel>``) AND no logging
    call is present in the first MAX_BODY_PEEK statements.
    """
    if not handler.body:
        return False
    if not _is_sentinel_fall_through(handler.body[0]):
        return False
    if _has_logging_call(handler):
        return False
    return True


def _scan_file(path: Path) -> tuple[int, list[int]]:
    """Return ``(total_excepts, silent_lineno_list)`` for ``path``.

    On SyntaxError the file is reported as ``(0, [])`` — the AST scan
    never aborts the run because of a single bad source file.
    """
    try:
        source = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return 0, []
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return 0, []

    total = 0
    silent_lines: list[int] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ExceptHandler):
            total += 1
            if _is_silent(node):
                silent_lines.append(node.lineno)
    return total, silent_lines


def _iter_targets(root: Path) -> Iterable[Path]:
    """Yield every ``*.py`` file under ``root/scripts/gamma`` sorted by path."""
    base = root / "scripts" / "gamma"
    if not base.is_dir():
        return iter(())
    return sorted(base.rglob("*.py"))


# ---------------------------------------------------------------------------
# Report construction
# ---------------------------------------------------------------------------


def build_report(repo_root: Path) -> dict:
    """Aggregate the AST scan into the JSON-serialisable report dict."""
    report: dict = {
        "total_excepts": 0,
        "silent_count": 0,
        "files": {},
        "sites_critical": [],
        "sites_other": [],
    }

    for path in _iter_targets(repo_root):
        total, silent_lines = _scan_file(path)
        if total == 0 and not silent_lines:
            continue
        rel = path.relative_to(repo_root).as_posix()
        report["files"][rel] = {
            "total_excepts": total,
            "silent_count": len(silent_lines),
            "silent_lines": silent_lines,
        }
        report["total_excepts"] += total
        report["silent_count"] += len(silent_lines)

        bucket = (
            "sites_critical"
            if path.name in CRITICAL_MODULES
            else "sites_other"
        )
        for ln in silent_lines:
            entry = {"file": rel, "line": ln}
            if bucket == "sites_critical":
                entry["module"] = path.name
            report[bucket].append(entry)

    return report


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "DEEP-06 census of silent exception handlers in "
            "scripts/gamma/**/*.py."
        ),
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help=(
            "Project root (default: parent directory of this script's "
            "containing folder)."
        ),
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help=(
            "Exit with code 1 when silent_count exceeds --ceiling. "
            "Without --strict the script always exits 0."
        ),
    )
    parser.add_argument(
        "--ceiling",
        type=int,
        default=9999,
        help=(
            "Maximum acceptable silent_count when --strict is given. "
            "Default 9999 (effectively no ceiling)."
        ),
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(list(sys.argv[1:] if argv is None else argv))
    report = build_report(args.repo_root)
    json.dump(report, sys.stdout, indent=2, sort_keys=False)
    sys.stdout.write("\n")
    if args.strict and report["silent_count"] > args.ceiling:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
