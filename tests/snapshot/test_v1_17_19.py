# -*- coding: utf-8 -*-
"""v1.17.19 delivery tests — F-293 books reorganization."""
from __future__ import annotations
import os, sys, tempfile, hashlib
from pathlib import Path
import pytest

SCRIPTS = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts")
)
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)


# ──────────────────────────────────────────────────────────────────
# F-293 — verify_books_inventory
# ──────────────────────────────────────────────────────────────────

def _make_fake_project(tmp_path: Path, index_filenames: list[str],
                       fs_filenames: list[str]) -> tuple[Path, Path]:
    """Создать минимальный project skeleton с INDEX.md и books_library/."""
    root = tmp_path / "proj"
    books = root / "books_library"
    refbooks = root / "references" / "books"
    books.mkdir(parents=True)
    refbooks.mkdir(parents=True)
    # INDEX.md
    lines = ["# Test INDEX", ""]
    for i, fn in enumerate(index_filenames, start=1):
        lines.append(f"### {i}. `{fn}`")
        lines.append("- Описание тестового файла.")
        lines.append("")
    (refbooks / "INDEX.md").write_text("\n".join(lines), encoding="utf-8")
    # books_library/ files
    for fn in fs_filenames:
        (books / fn).write_bytes(b"FAKE CONTENT " + fn.encode("utf-8"))
    return root, books


def test_F293_verifier_pass_when_match(tmp_path):
    from verify_books_inventory import verify_books_inventory
    root, _ = _make_fake_project(
        tmp_path,
        index_filenames=["a.pdf", "b.pdf"],
        fs_filenames=["a.pdf", "b.pdf"],
    )
    res = verify_books_inventory(root)
    assert res["ok"]
    assert res["n_files"] == 2
    assert not res["missing"]
    assert not res["untracked"]


def test_F293_verifier_detects_missing(tmp_path):
    from verify_books_inventory import verify_books_inventory
    root, _ = _make_fake_project(
        tmp_path,
        index_filenames=["a.pdf", "b.pdf", "c.pdf"],
        fs_filenames=["a.pdf"],   # b, c — отсутствуют
    )
    res = verify_books_inventory(root)
    assert not res["ok"]
    assert "b.pdf" in res["missing"]
    assert "c.pdf" in res["missing"]


def test_F293_verifier_detects_untracked(tmp_path):
    from verify_books_inventory import verify_books_inventory
    root, _ = _make_fake_project(
        tmp_path,
        index_filenames=["a.pdf"],
        fs_filenames=["a.pdf", "x.pdf", "y.pptx"],
    )
    res = verify_books_inventory(root)
    assert not res["ok"]
    assert "x.pdf" in res["untracked"]
    assert "y.pptx" in res["untracked"]


def test_F293_verifier_system_files_ignored(tmp_path):
    """desktop.ini / Thumbs.db / .DS_Store не должны считаться untracked."""
    from verify_books_inventory import verify_books_inventory
    root, books = _make_fake_project(
        tmp_path, index_filenames=["a.pdf"], fs_filenames=["a.pdf"],
    )
    (books / "desktop.ini").write_text("fake", encoding="utf-8")
    (books / "Thumbs.db").write_bytes(b"thumb-fake")
    (books / ".DS_Store").write_bytes(b"ds-fake")
    res = verify_books_inventory(root)
    assert res["ok"]
    assert "desktop.ini" not in res["untracked"]
    assert "Thumbs.db" not in res["untracked"]


def test_F293_verifier_env_override(tmp_path, monkeypatch):
    """env GAMMA_BOOKS_LIBRARY_DIR должен переопределять default путь."""
    from verify_books_inventory import (
        resolve_books_library_dir, verify_books_inventory,
    )
    # Создаём project с пустым books_library, а реальные книги в другой папке
    root = tmp_path / "proj"
    (root / "books_library").mkdir(parents=True)
    (root / "references" / "books").mkdir(parents=True)
    (root / "references" / "books" / "INDEX.md").write_text(
        "### 1. `external.pdf`\n", encoding="utf-8",
    )
    external = tmp_path / "external_library"
    external.mkdir()
    (external / "external.pdf").write_bytes(b"FAKE")

    monkeypatch.setenv("GAMMA_BOOKS_LIBRARY_DIR", str(external))
    resolved = resolve_books_library_dir(root)
    assert resolved == external.resolve()

    res = verify_books_inventory(root)
    assert res["ok"]
    assert res["n_files"] == 1


# ──────────────────────────────────────────────────────────────────
# F-293 — build_books_archive
# ──────────────────────────────────────────────────────────────────

def test_F293_books_archive_contains_manifest_and_files(tmp_path):
    import zipfile
    from build_books_archive import build_books_archive
    root, books = _make_fake_project(
        tmp_path,
        index_filenames=["x.pdf", "y.pdf"],
        fs_filenames=["x.pdf", "y.pdf"],
    )
    out_dir = tmp_path / "out"
    info = build_books_archive(
        root=root, out_dir=out_dir, books_dir=books, tag="2026-05-30",
    )
    assert info["n_files"] >= 4   # INDEX + 2 files + MANIFEST
    archive = Path(info["archive_path"])
    assert archive.exists()
    with zipfile.ZipFile(archive) as zf:
        names = set(zf.namelist())
    assert "INDEX.md" in names
    assert "MANIFEST.sha256" in names
    # Имена внутри books_library/ префиксованы
    assert any("books_library" in n and n.endswith("x.pdf") for n in names)


def test_F293_books_archive_rotation_preserves_prev(tmp_path):
    from build_books_archive import build_books_archive
    root, books = _make_fake_project(
        tmp_path, index_filenames=["x.pdf"], fs_filenames=["x.pdf"],
    )
    out_dir = tmp_path / "out"
    # Первый build
    info1 = build_books_archive(root, out_dir, books, "2026-05-30")
    # Подмена контента
    (books / "x.pdf").write_bytes(b"NEW CONTENT")
    # Второй build с тем же тегом → должен сделать .prev.zip
    info2 = build_books_archive(root, out_dir, books, "2026-05-30")
    assert Path(info2["archive_path"]).exists()
    prev = Path(info2["archive_path"]).with_name(
        Path(info2["archive_path"]).stem + ".prev.zip"
    )
    assert prev.exists(), "Rotation не создал .prev.zip"
