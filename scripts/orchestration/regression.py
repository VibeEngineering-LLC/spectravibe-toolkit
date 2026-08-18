#!/usr/bin/env python3
"""Concise regression runner — Agent D orchestration tool.

Wraps pytest with agent-friendly output: counts only for full suite,
short trace for --fast. Avoids dumping 1525 lines of test names.

Usage:
    python scripts/orchestration/regression.py              # full suite
    python scripts/orchestration/regression.py --fast PAT   # -k PAT with trace
    python scripts/orchestration/regression.py --counts     # only PASS/FAIL/XFAIL/XPASS counts
    python scripts/orchestration/regression.py --failed     # only show FAILED tests

Exit codes:
    0 — all green
    1 — failures present
    2 — wrong invocation
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

# Force UTF-8 on stdout/stderr to survive pytest output on Windows cp1251 consoles.
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

ROOT = Path(__file__).resolve().parents[2]
PYTEST_SUMMARY_RE = re.compile(
    r"=+\s*(?:(\d+)\s+passed)?(?:,\s*(\d+)\s+failed)?(?:,\s*(\d+)\s+skipped)?"
    r"(?:,\s*(\d+)\s+xfailed)?(?:,\s*(\d+)\s+xpassed)?[^=]*=+",
    re.IGNORECASE,
)


def _run(cmd: list[str], capture: bool = True) -> tuple[int, str, str]:
    p = subprocess.run(
        cmd, cwd=ROOT, capture_output=capture, text=True, encoding="utf-8", errors="replace"
    )
    return p.returncode, (p.stdout or ""), (p.stderr or "")


def _summarize(out: str) -> dict:
    """Parse pytest summary line."""
    last = ""
    for line in reversed(out.splitlines()):
        if "passed" in line or "failed" in line or "error" in line.lower():
            if "=" in line or "in " in line:
                last = line
                break
    m = re.search(r"(\d+)\s+passed", out)
    f = re.search(r"(\d+)\s+failed", out)
    s = re.search(r"(\d+)\s+skipped", out)
    xf = re.search(r"(\d+)\s+xfailed", out)
    xp = re.search(r"(\d+)\s+xpassed", out)
    return {
        "passed": int(m.group(1)) if m else 0,
        "failed": int(f.group(1)) if f else 0,
        "skipped": int(s.group(1)) if s else 0,
        "xfailed": int(xf.group(1)) if xf else 0,
        "xpassed": int(xp.group(1)) if xp else 0,
        "summary_line": last,
    }


def _failures(out: str) -> list[str]:
    return [ln.strip() for ln in out.splitlines() if "FAILED" in ln and "::" in ln]


def run_full(counts_only: bool = False, failed_only: bool = False) -> int:
    cmd = [sys.executable, "-m", "pytest", "-q", "--tb=line", "--no-header"]
    rc, out, err = _run(cmd)
    sm = _summarize(out)
    fails = _failures(out)

    if counts_only:
        print(f"P={sm['passed']} F={sm['failed']} S={sm['skipped']} XF={sm['xfailed']} XP={sm['xpassed']}")
        return rc

    print(f"REGRESSION: {sm['summary_line'] or 'unparsed'}")
    print(f"  pass={sm['passed']}  fail={sm['failed']}  skip={sm['skipped']}  xfail={sm['xfailed']}  xpass={sm['xpassed']}")
    if fails:
        print(f"\nFAILED ({len(fails)}):")
        for line in fails[:20]:
            print(f"  {line}")
        if len(fails) > 20:
            print(f"  ... and {len(fails) - 20} more")
    if failed_only and not fails:
        print("(no failures)")
    if rc != 0 and not fails and err:
        print(f"\nSTDERR tail:\n{err[-500:]}", file=sys.stderr)
    return rc


def run_fast(pattern: str) -> int:
    cmd = [sys.executable, "-m", "pytest", "-v", "--tb=short", "-k", pattern]
    rc, out, err = _run(cmd, capture=True)
    sys.stdout.write(out)
    if rc != 0 and err:
        sys.stderr.write(err)
    return rc


def main(argv: list[str]) -> int:
    if len(argv) > 1 and argv[1] == "--fast":
        if len(argv) < 3:
            print("usage: regression.py --fast PATTERN", file=sys.stderr)
            return 2
        return run_fast(argv[2])
    if len(argv) > 1 and argv[1] == "--counts":
        return run_full(counts_only=True)
    if len(argv) > 1 and argv[1] == "--failed":
        return run_full(failed_only=True)
    return run_full()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
