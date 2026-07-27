from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""verify_worktree_rerun — independent read-only worktree harness.

Draft per censor envelope 3e2a72bc (item #3). Implements the (a)(b)(c) guards
that close the fail-open hole exposed by re-run bw6hcfm5v:

  (a) ABORT if the worktree path already exists (no reuse — defends against
      orphan dirs with stale content).
  (b) git rev-parse HEAD must equal expected SHA before pytest runs.
  (c) ABORT on any worktree-add / rev-parse failure. Never fall through to
      pytest on a stale checkout.

Additional:
  (+) Captures subprocess.returncode after each pytest run -> explicit
      per-run rc list.
  (+) Emits ``all_rc_zero`` flag and exits non-zero if any run failed,
      so wrapper rc faithfully reflects pytest rc.

This is P3-tooling, separate commit from DEEP-06.

Status: live — landed in P3-tooling commit (Wave-2 closeout).
"""

import argparse
import datetime
import json
import os
import pathlib
import shlex
import subprocess
import sys
from typing import List, Tuple


def _run(cmd: List[str], cwd: pathlib.Path | None = None,
         capture: bool = True, check: bool = False) -> subprocess.CompletedProcess:
    """Run a subprocess; return CompletedProcess. ``check=True`` raises CalledProcessError."""
    return subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        capture_output=capture,
        text=True,
        check=check,
    )


def _now_stamp() -> str:
    """Stamp safe for use in a path / filename (no colons, no dots)."""
    now = datetime.datetime.now(datetime.UTC)
    return now.strftime("%Y%m%d%H%M%S")


def _unique_worktree_path(base: pathlib.Path, sha_short: str) -> pathlib.Path:
    """Compose a unique worktree path: <base>-<sha_short>-<UTC_timestamp>.

    Censor requirement (3e2a72bc #4): never reuse a path while .git/worktrees
    metadata for it is still locked (Google Drive lock). Timestamp suffix
    guarantees uniqueness across re-runs.
    """
    return base.with_name(f"{base.name}-{sha_short}-{_now_stamp()}")


def _assert_path_absent(path: pathlib.Path) -> None:
    """Guard (a): hard-exit if path already exists. No silent fall-through."""
    if path.exists():
        sys.stderr.write(
            f"FATAL (guard a): worktree path already exists: {path}\n"
            f"This script never reuses existing worktree paths. "
            f"Censor doctrine: pre-existing dir + worktree add 'already exists' "
            f"failure was the root cause of fail-open run bw6hcfm5v.\n"
        )
        sys.exit(2)


def _git_worktree_add(repo: pathlib.Path, wt: pathlib.Path, sha: str) -> None:
    """Guard (c): worktree add must succeed. Any failure -> hard exit."""
    r = _run(["git", "-C", str(repo), "worktree", "add", "--detach", str(wt), sha])
    if r.returncode != 0:
        sys.stderr.write(
            f"FATAL (guard c): git worktree add failed (rc={r.returncode}).\n"
            f"stdout: {r.stdout}\n"
            f"stderr: {r.stderr}\n"
            f"NEVER falling through to pytest on a stale or missing checkout.\n"
        )
        sys.exit(3)


def _assert_head(wt: pathlib.Path, expected_sha: str) -> str:
    """Guard (b): rev-parse HEAD must equal expected SHA (prefix or full).

    Returns full HEAD sha. Hard-exits on mismatch.
    """
    r = _run(["git", "-C", str(wt), "rev-parse", "HEAD"])
    if r.returncode != 0:
        sys.stderr.write(
            f"FATAL (guard b): git rev-parse HEAD failed in worktree {wt} "
            f"(rc={r.returncode}). stderr: {r.stderr}\n"
            f"This is the bw6hcfm5v failure mode: orphan dir without a valid "
            f".git pointer. ABORTING before pytest.\n"
        )
        sys.exit(4)
    head = r.stdout.strip()
    if not head.startswith(expected_sha):
        sys.stderr.write(
            f"FATAL (guard b): HEAD mismatch. expected prefix: {expected_sha}, "
            f"got: {head}. ABORTING before pytest.\n"
        )
        sys.exit(5)
    return head


def _pytest_run(wt: pathlib.Path, pytest_args: List[str]) -> Tuple[int, str]:
    """One pytest invocation in worktree. Returns (returncode, tail_of_output)."""
    cmd = [sys.executable, "-m", "pytest"] + pytest_args
    r = subprocess.run(
        cmd, cwd=str(wt), capture_output=True, text=True, check=False
    )
    # Keep only last 6 lines for log compactness; rc remains accurate.
    tail = "\n".join((r.stdout + r.stderr).splitlines()[-6:])
    return r.returncode, tail


def _cleanup(repo: pathlib.Path, wt: pathlib.Path) -> str:
    """Best-effort cleanup. Returns human-readable status. Never raises."""
    r = _run(["git", "-C", str(repo), "worktree", "remove", "--force", str(wt)])
    if r.returncode == 0:
        return "cleanup: OK"
    # Permission denied on Google Drive lock is an accepted cleanup artifact
    # per censor 3e2a72bc — pytest rc is decisive, wrapper rc is not.
    return (
        f"cleanup: failed (rc={r.returncode}) "
        f"— typically Google Drive lock on .git/worktrees/<name> metadata; "
        f"accepted as cleanup-only artifact when all pytest rc==0. "
        f"stderr-tail: {r.stderr.strip().splitlines()[-1] if r.stderr else ''}"
    )


def main(argv: List[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Independent read-only pytest re-run in a fresh, unique "
                    "git worktree. Honours censor (a)(b)(c) doctrine."
    )
    p.add_argument("--repo", required=True, type=pathlib.Path,
                   help="path to main git repository")
    p.add_argument("--sha", required=True,
                   help="full or unambiguous-prefix SHA to verify against")
    p.add_argument("--worktree-base", required=True, type=pathlib.Path,
                   help="base path for worktree; unique suffix is appended")
    p.add_argument("--n-runs", type=int, default=10,
                   help="how many pytest invocations (default: 10)")
    p.add_argument("--pytest-args", default="-n auto -q -ra",
                   help="quoted pytest CLI args (default: '-n auto -q -ra')")
    p.add_argument("--log", type=pathlib.Path, default=None,
                   help="optional path to write a full structured log JSON")
    args = p.parse_args(argv)

    repo: pathlib.Path = args.repo.resolve()
    if not (repo / ".git").exists():
        sys.stderr.write(f"FATAL: --repo {repo} is not a git repo\n")
        return 1

    sha_short = args.sha[:7]
    wt = _unique_worktree_path(args.worktree_base.resolve(), sha_short)

    print(f"=== verify_worktree_rerun ===")
    print(f"repo : {repo}")
    print(f"sha  : {args.sha}")
    print(f"wt   : {wt}")
    print(f"n    : {args.n_runs}")

    # ── Guard (a): path absent ──
    _assert_path_absent(wt)
    print(f"guard (a) PASS: path is fresh")

    # ── Guard (c): worktree add succeeds ──
    _git_worktree_add(repo, wt, args.sha)
    print(f"guard (c) PASS: worktree add rc=0")

    # ── Guard (b): HEAD matches ──
    head_full = _assert_head(wt, args.sha)
    print(f"guard (b) PASS: HEAD == {head_full}")

    # ── pytest runs ──
    rcs: List[int] = []
    summaries: List[str] = []
    for i in range(1, args.n_runs + 1):
        print(f"=== run {i}/{args.n_runs} pytest ===")
        rc, tail = _pytest_run(wt, shlex.split(args.pytest_args))
        rcs.append(rc)
        summaries.append(tail)
        print(tail)
        print(f"run {i} pytest rc={rc}")

    all_zero = all(r == 0 for r in rcs)
    print("=== SUMMARY: per-run pytest exit codes ===")
    for i, rc in enumerate(rcs, 1):
        print(f"  run {i:2d}: rc={rc}")
    print(f"all_rc_zero: {all_zero}")

    # ── Cleanup (best effort) ──
    cleanup_status = _cleanup(repo, wt)
    print(f"=== {cleanup_status} ===")

    # ── Optional structured log ──
    if args.log:
        args.log.parent.mkdir(parents=True, exist_ok=True)
        args.log.write_text(
            json.dumps(
                {
                    "repo": str(repo),
                    "sha_requested": args.sha,
                    "head_actual": head_full,
                    "worktree": str(wt),
                    "n_runs": args.n_runs,
                    "pytest_args": args.pytest_args,
                    "per_run_rc": rcs,
                    "per_run_tail": summaries,
                    "all_rc_zero": all_zero,
                    "cleanup": cleanup_status,
                },
                ensure_ascii=False,
                indent=2,
            ),
            encoding="utf-8",
        )
        print(f"=== log -> {args.log} ===")

    # ── Wrapper rc faithfully reflects pytest rc ──
    # Censor doctrine: wrapper rc 0 iff all pytest rc 0.
    return 0 if all_zero else 1


if __name__ == "__main__":
    raise SystemExit(main())
