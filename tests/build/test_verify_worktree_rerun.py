"""Tests for scripts/verify_worktree_rerun.py.

Status: live — landed alongside scripts/verify_worktree_rerun.py.

Test inventory (7 tests; docstring synced to actual function names per
censor cosmetic finding 2026-06-06):

  Guard coverage:
    - test_guard_a_aborts_on_existing_path                 ↔ guard (a) red-without-fix
    - test_guard_b_aborts_on_head_mismatch                 ↔ guard (b) *via guard (c) path*
                                                             (literal non-existent SHA
                                                              triggers (c) exit 3 first;
                                                              kept as guard-chain
                                                              coverage, not isolated
                                                              guard (b) test)
    - test_guard_b_aborts_when_head_prefix_does_not_match  ↔ guard (b) ISOLATED
                                                             (monkeypatched rev-parse
                                                              → exit 5)
    - test_guard_c_aborts_on_worktree_add_failure          ↔ guard (c) red-without-fix

  Wrapper rc contract:
    - test_wrapper_rc_zero_when_all_runs_pass              ↔ closes bw6hcfm5v hole
                                                             (cleanup error MUST NOT
                                                              mask pytest rc when all
                                                              runs are GREEN)
    - test_wrapper_rc_nonzero_when_any_pytest_fails        ↔ any per-run rc != 0
                                                             propagates to wrapper rc

  Structured log:
    - test_structured_log_written_when_log_path_given      ↔ --log emits JSON with
                                                             per_run_rc array +
                                                             all_rc_zero flag +
                                                             sha_requested/head_actual

Run from repo root:
    python -m pytest tests/build/test_verify_worktree_rerun.py -v
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest


HELPER_PATH = Path(__file__).resolve().parents[2] / "scripts" / "verify_worktree_rerun.py"


def _import_helper():
    """Import the helper as a module — works regardless of CWD."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("verify_worktree_rerun", HELPER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[attr-defined]
    return mod


@pytest.fixture()
def helper():
    if not HELPER_PATH.exists():
        pytest.skip(f"helper not landed yet: {HELPER_PATH}")
    return _import_helper()


def _mk_repo(tmp_path: Path) -> Path:
    """Initialise a throwaway git repo with one commit. Returns repo path."""
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init", "-q"], cwd=repo, check=True)
    subprocess.run(["git", "config", "user.email", "t@t"], cwd=repo, check=True)
    subprocess.run(["git", "config", "user.name", "t"], cwd=repo, check=True)
    (repo / "README.md").write_text("seed", encoding="utf-8")
    subprocess.run(["git", "add", "."], cwd=repo, check=True)
    subprocess.run(
        ["git", "commit", "-q", "-m", "seed"], cwd=repo, check=True
    )
    return repo


def _head(repo: Path) -> str:
    return subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "HEAD"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()


# ───────────────────────── Guard (a) ─────────────────────────

def test_guard_a_aborts_on_existing_path(helper, tmp_path):
    """Pre-existing worktree path -> hard exit before any git op.

    Red-without-fix: if _assert_path_absent is removed from helper, this test
    would let helper attempt `worktree add` on existing dir, fail with
    'already exists', and (without guard (c)) silently fall through to pytest
    on the stale content. The (a) check stops it earliest possible.
    """
    repo = _mk_repo(tmp_path)
    sha = _head(repo)
    # Pre-create the would-be worktree path:
    existing = tmp_path / "wt"
    (existing.parent / f"wt-{sha[:7]}").mkdir(parents=True, exist_ok=False)
    # Patch _unique_worktree_path so we know the exact path the helper picks:
    target = existing.parent / f"wt-{sha[:7]}"
    helper._unique_worktree_path = lambda base, sha_short: target  # type: ignore[assignment]

    with pytest.raises(SystemExit) as exc:
        helper.main([
            "--repo", str(repo),
            "--sha", sha,
            "--worktree-base", str(existing),
            "--n-runs", "1",
        ])
    assert exc.value.code == 2, "guard (a) must exit code 2 on existing path"


# ───────────────────────── Guard (b) ─────────────────────────

def test_guard_b_aborts_on_head_mismatch(helper, tmp_path, monkeypatch):
    """Guard-chain coverage: non-existent SHA literal -> guard (c) exit 3 (NOT isolated guard (b)).

    Censor cosmetic finding (2026-06-06): this test feeds ``deadbeef…`` which fails at
    ``git worktree add`` (no such object) — guard (c) triggers FIRST with exit 3 before
    guard (b)'s rev-parse-HEAD check is reached. The assertion at :128 accepts both
    exit codes ``(3, 5)`` because the actual path taken is (c), not (b).

    Kept as guard-CHAIN coverage (pytest must not run when ANY guard fails) — NOT a
    red-without-fix for guard (b) in isolation. Guard (b) isolated coverage lives in
    ``test_guard_b_aborts_when_head_prefix_does_not_match`` (:134, exit 5 via
    monkeypatched rev-parse).
    """
    repo = _mk_repo(tmp_path)
    actual_head = _head(repo)
    wrong_sha = "deadbeefdeadbeefdeadbeefdeadbeefdeadbeef"

    # Use real path + real worktree add (cheap on tiny seed repo):
    wt_base = tmp_path / "wt"

    pytest_called = {"n": 0}

    def _fake_pytest_run(wt_path, args):
        pytest_called["n"] += 1
        return (0, "OK")

    monkeypatch.setattr(helper, "_pytest_run", _fake_pytest_run)

    with pytest.raises(SystemExit) as exc:
        helper.main([
            "--repo", str(repo),
            "--sha", wrong_sha,
            "--worktree-base", str(wt_base),
            "--n-runs", "1",
        ])

    # NOTE: worktree add WILL fail because wrong_sha doesn't exist.
    # That's guard (c) catching it first — exit code 3.
    # If we want to test guard (b) in isolation, we need a sha that exists
    # but mismatches expected. Helper already accepts prefix match, so the
    # mismatch test is more naturally done by feeding a wrong prefix:
    assert exc.value.code in (3, 5), (
        f"expected guard (b) or (c) abort, got {exc.value.code}"
    )
    assert pytest_called["n"] == 0, "pytest must NOT run on guard failure"


def test_guard_b_aborts_when_head_prefix_does_not_match(helper, tmp_path, monkeypatch):
    """Variant: simulate rev-parse returning unexpected sha via monkeypatch."""
    repo = _mk_repo(tmp_path)
    real_head = _head(repo)
    expected = real_head  # we pass the real SHA

    # But monkeypatch _assert_head's underlying _run to return a *different*
    # head when the helper calls rev-parse. We do this by wrapping subprocess.run.
    original_run = subprocess.run

    def fake_run(cmd, *a, **kw):
        if "rev-parse" in cmd and "HEAD" in cmd:
            return SimpleNamespace(
                returncode=0,
                stdout="cafebabecafebabecafebabecafebabecafebabe\n",
                stderr="",
            )
        return original_run(cmd, *a, **kw)

    monkeypatch.setattr(subprocess, "run", fake_run)

    pytest_called = {"n": 0}
    def _fake_pytest_run(wt_path, args):
        pytest_called["n"] += 1
        return (0, "OK")
    monkeypatch.setattr(helper, "_pytest_run", _fake_pytest_run)

    with pytest.raises(SystemExit) as exc:
        helper.main([
            "--repo", str(repo),
            "--sha", expected,
            "--worktree-base", str(tmp_path / "wt"),
            "--n-runs", "1",
        ])
    assert exc.value.code == 5, "guard (b) must exit code 5 on HEAD prefix mismatch"
    assert pytest_called["n"] == 0, "pytest must NOT run on guard (b) failure"


# ───────────────────────── Guard (c) ─────────────────────────

def test_guard_c_aborts_on_worktree_add_failure(helper, tmp_path, monkeypatch):
    """worktree add returns non-zero -> hard exit before pytest."""
    repo = _mk_repo(tmp_path)
    sha = _head(repo)

    original_run = subprocess.run

    def fake_run(cmd, *a, **kw):
        if "worktree" in cmd and "add" in cmd:
            return SimpleNamespace(
                returncode=128,
                stdout="",
                stderr="fatal: simulated add failure",
            )
        return original_run(cmd, *a, **kw)

    monkeypatch.setattr(subprocess, "run", fake_run)

    pytest_called = {"n": 0}
    def _fake_pytest_run(wt_path, args):
        pytest_called["n"] += 1
        return (0, "OK")
    monkeypatch.setattr(helper, "_pytest_run", _fake_pytest_run)

    with pytest.raises(SystemExit) as exc:
        helper.main([
            "--repo", str(repo),
            "--sha", sha,
            "--worktree-base", str(tmp_path / "wt"),
            "--n-runs", "1",
        ])
    assert exc.value.code == 3, "guard (c) must exit code 3 on worktree add failure"
    assert pytest_called["n"] == 0, "pytest must NOT run on guard (c) failure"


# ───────────────────────── Wrapper rc ──────────────────────

def test_wrapper_rc_zero_when_all_runs_pass(helper, tmp_path, monkeypatch):
    """Happy path: all guards pass, all pytest rc=0, cleanup OK or not — wrapper exits 0."""
    repo = _mk_repo(tmp_path)
    sha = _head(repo)

    monkeypatch.setattr(helper, "_pytest_run", lambda wt, args: (0, "passed"))
    # Simulate cleanup failure on Google Drive lock — must NOT affect wrapper rc:
    def fake_cleanup(repo, wt):
        return "cleanup: failed (rc=128) — Permission denied"
    monkeypatch.setattr(helper, "_cleanup", fake_cleanup)

    rc = helper.main([
        "--repo", str(repo),
        "--sha", sha,
        "--worktree-base", str(tmp_path / "wt"),
        "--n-runs", "3",
    ])
    assert rc == 0, (
        "wrapper rc must be 0 when ALL pytest runs are 0, even if cleanup fails. "
        "Closes bw6hcfm5v hole (wrapper rc != pytest rc)."
    )


def test_wrapper_rc_nonzero_when_any_pytest_fails(helper, tmp_path, monkeypatch):
    """Any pytest rc != 0 -> wrapper rc != 0."""
    repo = _mk_repo(tmp_path)
    sha = _head(repo)

    counter = {"n": 0}
    def fake_pytest(wt, args):
        counter["n"] += 1
        return (1 if counter["n"] == 2 else 0, "failed")
    monkeypatch.setattr(helper, "_pytest_run", fake_pytest)
    monkeypatch.setattr(helper, "_cleanup", lambda r, w: "cleanup: OK")

    rc = helper.main([
        "--repo", str(repo),
        "--sha", sha,
        "--worktree-base", str(tmp_path / "wt"),
        "--n-runs", "3",
    ])
    assert rc == 1, "wrapper rc must be 1 if any per-run pytest rc != 0"


def test_structured_log_written_when_log_path_given(helper, tmp_path, monkeypatch):
    """--log path -> JSON file with per_run_rc array and all_rc_zero flag."""
    repo = _mk_repo(tmp_path)
    sha = _head(repo)

    monkeypatch.setattr(helper, "_pytest_run", lambda wt, args: (0, "OK"))
    monkeypatch.setattr(helper, "_cleanup", lambda r, w: "cleanup: OK")
    log_path = tmp_path / "out" / "rerun.json"

    helper.main([
        "--repo", str(repo),
        "--sha", sha,
        "--worktree-base", str(tmp_path / "wt"),
        "--n-runs", "2",
        "--log", str(log_path),
    ])

    assert log_path.exists(), "structured log must be written when --log given"
    data = json.loads(log_path.read_text(encoding="utf-8"))
    assert data["per_run_rc"] == [0, 0]
    assert data["all_rc_zero"] is True
    assert data["sha_requested"] == sha
    assert data["head_actual"].startswith(sha[:7])
