from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""
version_check.py — CI gate: README badge == SKILL_VERSION == git tag.

Part of Project #5 P0/P1 remediation (CRITIQUE_AND_PLAN.md QUAL-PROC-01).

Usage:
    python scripts/version_check.py              # full check (README + SKILL_VERSION + git tag)
    python scripts/version_check.py --allow-no-tag  # skip tag check (in-progress dev clone)

Exit codes:
    0 — all markers match
    1 — mismatch found (prints diff)
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent


def read_skill_version() -> str | None:
    """
    Read SKILL_VERSION from scripts/gamma/reporting/json_report.py via regex.
    Does NOT import the module — keeps this script dependency-free.
    """
    fpath = ROOT / "scripts" / "gamma" / "reporting" / "json_report.py"
    text = fpath.read_text(encoding="utf-8")
    m = re.search(r'^SKILL_VERSION\s*=\s*"(v[\d\.]+)"', text, re.MULTILINE)
    if m:
        return m.group(1)
    return None


def read_readme_version() -> str | None:
    """
    Read version badge from README.md.
    Looks for first vX.Y.Z on a line starting with '> **'.
    """
    fpath = ROOT / "README.md"
    text = fpath.read_text(encoding="utf-8")
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("> **"):
            m = re.search(r"v(\d+\.\d+\.\d+)", stripped)
            if m:
                return "v" + m.group(1)
    return None


def read_git_tag() -> str | None:
    """Return a tag that points exactly at HEAD, or None.

    Uses `git tag --points-at HEAD` instead of `git describe --tags --abbrev=0`
    so a stale ancestor tag (e.g. v1.30.0 reachable from HEAD when HEAD is at
    v1.30.1 but the v1.30.1 tag is not yet visible) returns None rather than
    the stale value. This makes the gate race-immune in CI when the branch
    push and the tag push are not atomic, and lets --allow-no-tag behave
    correctly in dev clones with un-tagged HEAD.
    """
    try:
        result = subprocess.run(
            ["git", "tag", "--points-at", "HEAD"],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            timeout=10,
        )
        if result.returncode == 0:
            tags = [
                t.strip()
                for t in result.stdout.splitlines()
                if t.strip().startswith("v")
            ]
            if tags:
                return tags[0]
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description="Check version markers consistency.")
    parser.add_argument(
        "--allow-no-tag",
        action="store_true",
        help="Skip git-tag check if no tag is found (useful for in-progress dev clones).",
    )
    args = parser.parse_args()

    skill_version = read_skill_version()
    readme_version = read_readme_version()
    git_tag = read_git_tag()

    ok = True
    messages: list[str] = []

    if skill_version is None:
        messages.append("FAIL: Could not read SKILL_VERSION from scripts/gamma/reporting/json_report.py")
        ok = False

    if readme_version is None:
        messages.append("FAIL: Could not find version badge in README.md (expected '> **vX.Y.Z' line)")
        ok = False

    if ok and skill_version != readme_version:
        messages.append(
            f"MISMATCH: README badge={readme_version!r} != SKILL_VERSION={skill_version!r}"
        )
        ok = False

    if args.allow_no_tag:
        if git_tag is None:
            messages.append("NOTE: no tag points at HEAD — skipping tag check (--allow-no-tag)")
        else:
            if ok and git_tag != skill_version:
                messages.append(
                    f"MISMATCH: git tag={git_tag!r} != SKILL_VERSION={skill_version!r}"
                )
                ok = False
    else:
        if git_tag is None:
            messages.append(
                "FAIL: no git tag points at HEAD. "
                "Either tag the release commit (push branch+tag together with --follow-tags) "
                "or pass --allow-no-tag to skip this check."
            )
            ok = False
        else:
            if ok and git_tag != skill_version:
                messages.append(
                    f"MISMATCH: git tag={git_tag!r} != SKILL_VERSION={skill_version!r}"
                )
                ok = False

    if ok:
        version = skill_version or "(unknown)"
        print(f"OK: README={version} == SKILL_VERSION={version} == tag={git_tag or '(skipped)'}")
        return 0
    else:
        print("VERSION CHECK FAILED:")
        for msg in messages:
            print(f"  {msg}")
        print(f"  README={readme_version!r}  SKILL_VERSION={skill_version!r}  git_tag={git_tag!r}")
        return 1


if __name__ == "__main__":
    sys.exit(main())