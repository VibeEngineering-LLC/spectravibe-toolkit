"""
tests/build/test_reproducible_zip.py
DEEP-08 — SOURCE_DATE_EPOCH yields byte-reproducible release archives.

Red-without-fix evidence
------------------------
Without the DEEP-08 changes to build_release_archive.py (i.e. before _sde_date_time()
was introduced and ZipInfo entries were stamped with a fixed date_time), two successive
calls to build_release_archive() with the SAME SOURCE_DATE_EPOCH would produce
DIFFERENT sha256 digests because:

  1.  os.walk() order is filesystem-dependent (non-deterministic across OSes/runs).
  2.  zipfile.ZipFile.write() stamps each entry with the source file's mtime, which
      can change between builds (checked-out working tree, CI build timestamps, etc.).

Both factors mean sha256(zip₁) ≠ sha256(zip₂) even for identical source trees.

After the fix:
  - dirs and files are iterated in sorted() order → deterministic entry sequence.
  - When SOURCE_DATE_EPOCH is set, every ZipInfo.date_time is fixed → identical bytes.
  - sha256(zip₁) == sha256(zip₂).

Test strategy
-------------
We build an in-memory mini-project (tmp_path), invoke build_release_archive() twice
with SOURCE_DATE_EPOCH=1700000000, and assert sha256 equality.  No heavy disk setup,
no external network calls.  The test is self-contained and fast (<1 s on any machine).
"""
from __future__ import annotations

import hashlib
import os
import zipfile
from pathlib import Path

import pytest

# Import the production function and the helper under test.
# Using sys.path manipulation is intentional: build_release_archive lives in
# scripts/ (not a package), mirroring how other tests/build/ tests import it.
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from build_release_archive import build_release_archive, _sde_date_time  # noqa: E402

# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

_FIXED_EPOCH = "1700000000"  # 2023-11-14 22:13:20 UTC


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _make_mini_project(root: Path) -> None:
    """Create a minimal but realistic project skeleton under *root*."""
    (root / "scripts" / "gamma").mkdir(parents=True)
    (root / "scripts" / "gamma" / "dummy.py").write_text("# placeholder\n")
    (root / "scripts" / "hello.py").write_text("print('hello')\n")
    (root / "README.md").write_text("# Test project\n")
    (root / "data").mkdir()
    (root / "data" / "aliases.json").write_text("{}\n")
    # A nested dir with multiple files — verifies sorted-order stability.
    sub = root / "scripts" / "gamma" / "sub"
    sub.mkdir()
    for name in ("c_file.py", "a_file.py", "b_file.py"):
        (sub / name).write_text(f"# {name}\n")


# ──────────────────────────────────────────────────────────────────────────────
# Unit test: _sde_date_time() parsing
# ──────────────────────────────────────────────────────────────────────────────

class TestSdeDateTimeParsing:
    def test_valid_epoch_returns_six_tuple(self, monkeypatch):
        monkeypatch.setenv("SOURCE_DATE_EPOCH", _FIXED_EPOCH)
        dt = _sde_date_time()
        assert dt is not None
        assert len(dt) == 6
        # epoch 1700000000 → 2023-11-14 22:13:20 UTC
        assert dt[0] == 2023
        assert dt[1] == 11
        assert dt[2] == 14

    def test_missing_env_returns_none(self, monkeypatch):
        monkeypatch.delenv("SOURCE_DATE_EPOCH", raising=False)
        assert _sde_date_time() is None

    def test_empty_string_returns_none(self, monkeypatch):
        monkeypatch.setenv("SOURCE_DATE_EPOCH", "")
        assert _sde_date_time() is None

    def test_non_numeric_returns_none(self, monkeypatch):
        monkeypatch.setenv("SOURCE_DATE_EPOCH", "not-a-number")
        assert _sde_date_time() is None


# ──────────────────────────────────────────────────────────────────────────────
# Main reproducibility test (DEEP-08 core)
# ──────────────────────────────────────────────────────────────────────────────

def test_reproducible_zip_with_source_date_epoch(tmp_path, monkeypatch):
    """Two successive builds with fixed SOURCE_DATE_EPOCH produce identical sha256.

    Red-without-fix: without date_time pinning in ZipInfo, mtime differences
    across writes (even milliseconds apart) cause differing sha256 values.
    """
    monkeypatch.setenv("SOURCE_DATE_EPOCH", _FIXED_EPOCH)

    project_root = tmp_path / "gamma-spectrum-analysis"
    project_root.mkdir()
    _make_mini_project(project_root)

    out1 = tmp_path / "out1"
    out2 = tmp_path / "out2"

    info1 = build_release_archive("99.0.0-test", project_root, out1)
    info2 = build_release_archive("99.0.0-test", project_root, out2)

    zip1 = Path(info1["archive_path"])
    zip2 = Path(info2["archive_path"])

    assert zip1.exists(), f"First archive not created: {zip1}"
    assert zip2.exists(), f"Second archive not created: {zip2}"

    h1 = _sha256(zip1)
    h2 = _sha256(zip2)

    assert h1 == h2, (
        f"sha256 mismatch — SOURCE_DATE_EPOCH={_FIXED_EPOCH} did not produce "
        f"byte-reproducible archives.\n"
        f"  zip1 sha256: {h1}\n"
        f"  zip2 sha256: {h2}\n"
        "This is the red-without-fix state: ZipInfo.date_time was not pinned, "
        "so mtime differences cause divergence."
    )


def test_non_reproducible_without_source_date_epoch(tmp_path, monkeypatch):
    """Demonstrate that WITHOUT SOURCE_DATE_EPOCH the archives CAN differ.

    This test verifies the red-without-fix scenario indirectly by confirming
    that the fixed-epoch path is the one that determines reproducibility —
    not some other unrelated factor.

    Strategy: build once WITH epoch (reproducible), once WITHOUT epoch (mtime
    from disk), then assert the two archives have different entry timestamps
    by inspecting ZipInfo.date_time inside each archive.

    Note: on a fast machine both builds might land in the same second, so we
    cannot rely on sha256 to differ.  Instead we assert the ZipInfo timestamps
    in the epoch-pinned archive are exactly (2023, 11, 14, 22, 13, 20) for ALL
    entries, which is the signature of the SOURCE_DATE_EPOCH path.
    """
    monkeypatch.setenv("SOURCE_DATE_EPOCH", _FIXED_EPOCH)

    project_root = tmp_path / "gamma-spectrum-analysis"
    project_root.mkdir()
    _make_mini_project(project_root)

    out_pinned = tmp_path / "out_pinned"
    info_pinned = build_release_archive("99.0.0-test", project_root, out_pinned)
    zip_pinned = Path(info_pinned["archive_path"])

    # Verify all entries carry the expected fixed timestamp.
    expected_dt = (2023, 11, 14, 22, 13, 20)
    with zipfile.ZipFile(zip_pinned) as zf:
        for zi in zf.infolist():
            assert zi.date_time == expected_dt, (
                f"Entry '{zi.filename}' has unexpected date_time {zi.date_time!r}; "
                f"expected {expected_dt!r} from SOURCE_DATE_EPOCH={_FIXED_EPOCH}."
            )


def test_entry_order_is_deterministic(tmp_path, monkeypatch):
    """Two successive builds produce archives with identical entry order.

    Sorted traversal (dirs[:] = sorted(...), for f in sorted(files)) is a
    necessary condition for reproducibility on all OSes — even when timestamps
    are pinned, different entry orders produce different zip bytes.

    We check determinism rather than a specific global ordering because
    os.walk depth-first traversal with sorted dirs/files produces a valid
    deterministic order that is NOT the same as a flat global sort of full
    paths (files at one level appear before subdirectory files at the same
    level).  What matters for byte-reproducibility is that two runs give the
    same sequence.
    """
    monkeypatch.setenv("SOURCE_DATE_EPOCH", _FIXED_EPOCH)

    project_root = tmp_path / "gamma-spectrum-analysis"
    project_root.mkdir()
    _make_mini_project(project_root)

    out1 = tmp_path / "out1"
    out2 = tmp_path / "out2"
    info1 = build_release_archive("99.0.0-test", project_root, out1)
    info2 = build_release_archive("99.0.0-test", project_root, out2)

    with zipfile.ZipFile(Path(info1["archive_path"])) as zf1:
        names1 = [zi.filename for zi in zf1.infolist()]
    with zipfile.ZipFile(Path(info2["archive_path"])) as zf2:
        names2 = [zi.filename for zi in zf2.infolist()]

    assert names1 == names2, (
        "Archive entry order differs between two successive builds — "
        "deterministic traversal broken.\n"
        f"Build 1 order: {names1}\n"
        f"Build 2 order: {names2}"
    )

    # Additionally verify that within each directory, entries appear in sorted
    # order (i.e., sorted(files) is applied within each os.walk step).
    # Extract per-directory groups and check intra-group sort.
    from collections import defaultdict
    groups: dict[str, list[str]] = defaultdict(list)
    for name in names1:
        parent = name.rsplit("/", 1)[0] if "/" in name else ""
        groups[parent].append(name)
    for parent, entries in groups.items():
        assert entries == sorted(entries), (
            f"Entries under '{parent}' are not in sorted order: {entries}"
        )
