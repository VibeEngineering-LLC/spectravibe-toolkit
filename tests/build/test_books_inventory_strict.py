"""
tests/build/test_books_inventory_strict.py
V126-03 — Guard that `--strict` is a superset of default mode.

Red-without-fix: before the V126-03 one-line fix to main(), the `--strict` branch
returned 0 even when the inventory diff (references/books/INDEX.md vs books_library/)
indicated a mismatch, as long as all INDEX.md files were present and no archives
existed.  This silently disarmed the CI gate.

After the fix the `--strict` exit code is:
  0   iff   strict_res["ok"]  AND  res["ok"]
  1   otherwise  (either strict checks OR inventory mismatch, or both)

Test: test_strict_mode_fails_on_inventory_mismatch
  Arrange: all required INDEX.md files present + no archives in tree
           (so strict_res["ok"] is True), but one untracked PDF in
           books_library/ that is NOT listed in references/books/INDEX.md
           (so res["ok"] is False).
  Assert: main() returns 1.
  Without fix: main() returned 0.  With fix: 1.
"""

from pathlib import Path

import verify_books_inventory as vbi


def test_strict_mode_fails_on_inventory_mismatch(tmp_path, monkeypatch):
    """--strict must exit 1 when inventory diff exists, even if INDEX.md files are present."""
    # Arrange: create a books_library/ with all required INDEX.md files present
    # but one extra untracked file (inventory mismatch).
    books_dir = tmp_path / "books_library"
    books_dir.mkdir()

    # Master INDEX.md — lists only existing_book.pdf (not untracked_book.pdf).
    (books_dir / "INDEX.md").write_text("### 1. `existing_book.pdf`\n", encoding="utf-8")

    # Required subtrees with their INDEX.md files.
    for sub in ["01_methodology_pdf", "_corpus_pages"]:
        (books_dir / sub).mkdir()
        (books_dir / sub / "INDEX.md").write_text("", encoding="utf-8")

    lsrm = books_dir / "Документация ЛСРМ"
    lsrm.mkdir()
    (lsrm / "INDEX.md").write_text("", encoding="utf-8")

    # _corpus_pages/README.md must be present.
    (books_dir / "_corpus_pages" / "README.md").write_text("", encoding="utf-8")

    # Actual files on disk: one tracked, one untracked.
    (books_dir / "existing_book.pdf").write_bytes(b"pdf-content")
    (books_dir / "untracked_book.pdf").write_bytes(b"pdf-content")

    # references/books/INDEX.md references only existing_book.pdf.
    ref_books = tmp_path / "references" / "books"
    ref_books.mkdir(parents=True)
    (ref_books / "INDEX.md").write_text("### 1. `existing_book.pdf`\n", encoding="utf-8")

    # Act: call main() with --strict and explicit paths so the tmp_path is used.
    result = vbi.main(["--root", str(tmp_path), "--books-dir", str(books_dir), "--strict"])

    # Assert:
    # Without the V126-03 fix: result == 0 (WRONG — inventory mismatch ignored).
    # After the V126-03 fix:   result == 1 (CORRECT — mismatch propagates under --strict).
    assert result == 1, (
        "--strict must exit 1 when inventory mismatch exists "
        "(untracked_book.pdf not in references/books/INDEX.md), "
        "even when all INDEX.md files are present and no archives exist. "
        "Failing assertion: `result == 1` (got 0 before V126-03 fix)."
    )


def test_strict_mode_ok_when_fully_consistent(tmp_path):
    """--strict must exit 0 when inventory is OK and all INDEX.md files are present."""
    books_dir = tmp_path / "books_library"
    books_dir.mkdir()

    # Master INDEX.md — lists exactly the file present on disk.
    (books_dir / "INDEX.md").write_text("### 1. `existing_book.pdf`\n", encoding="utf-8")

    for sub in ["01_methodology_pdf", "_corpus_pages"]:
        (books_dir / sub).mkdir()
        (books_dir / sub / "INDEX.md").write_text("", encoding="utf-8")

    lsrm = books_dir / "Документация ЛСРМ"
    lsrm.mkdir()
    (lsrm / "INDEX.md").write_text("", encoding="utf-8")

    (books_dir / "_corpus_pages" / "README.md").write_text("", encoding="utf-8")
    (books_dir / "existing_book.pdf").write_bytes(b"pdf-content")

    ref_books = tmp_path / "references" / "books"
    ref_books.mkdir(parents=True)
    (ref_books / "INDEX.md").write_text("### 1. `existing_book.pdf`\n", encoding="utf-8")

    result = vbi.main(["--root", str(tmp_path), "--books-dir", str(books_dir), "--strict"])
    assert result == 0, (
        "--strict must exit 0 when inventory is fully consistent "
        "and all INDEX.md files are present."
    )


def test_default_mode_still_fails_on_inventory_mismatch(tmp_path):
    """Default (non-strict) mode must still exit 1 on inventory mismatch (regression guard)."""
    books_dir = tmp_path / "books_library"
    books_dir.mkdir()
    (books_dir / "untracked_book.pdf").write_bytes(b"pdf-content")

    ref_books = tmp_path / "references" / "books"
    ref_books.mkdir(parents=True)
    # INDEX.md lists a file that is not on disk (missing) — triggers mismatch.
    (ref_books / "INDEX.md").write_text("### 1. `expected_book.pdf`\n", encoding="utf-8")

    result = vbi.main(["--root", str(tmp_path), "--books-dir", str(books_dir)])
    assert result == 1, (
        "Default (non-strict) mode must exit 1 when inventory mismatch exists."
    )


def test_index_md_total_count_matches_disk_aggregation(tmp_path):
    """DEEP-10: INDEX.md parsed entry count must equal top-level fs file count.

    Red-without-fix scenario: INDEX.md uses ``### N. `filename``` headings for
    files that are NOT in books_library/ top-level (e.g. nested in Документация
    ЛСРМ/).  The HEADING_RE parser counts those entries, but scan_books_library()
    does not find them — leading to len(index_names) != len(fs_names) and
    verify_books_inventory returning ok=False.

    After the DEEP-10 fix those headings use ``### N-nested. filename`` (no
    backtick-wrapped filename), so HEADING_RE skips them.  Now the parsed count
    equals the top-level fs count and ok=True.
    """
    books_dir = tmp_path / "books_library"
    books_dir.mkdir()

    ref_books = tmp_path / "references" / "books"
    ref_books.mkdir(parents=True)

    # --- Scenario that is RED without the DEEP-10 fix ---
    # INDEX.md has a nested entry in ### N. `filename` format (old style).
    # The file only exists inside a nested subfolder, NOT at top-level.
    nested_subdir = books_dir / "Документация ЛСРМ"
    nested_subdir.mkdir()
    (nested_subdir / "nested_doc.pdf").write_bytes(b"pdf-content")

    # Top-level file that IS tracked.
    (books_dir / "top_level_book.pdf").write_bytes(b"pdf-content")

    # INDEX.md (pre-fix style): both files under ### N. `...` — the nested
    # one should NOT be here, causing mismatch.
    bad_index = (
        "### 1. `top_level_book.pdf`\n"
        "### 2. `nested_doc.pdf`\n"  # file only in subfolder — WRONG
    )
    (ref_books / "INDEX.md").write_text(bad_index, encoding="utf-8")

    result_bad = vbi.main(["--root", str(tmp_path), "--books-dir", str(books_dir)])
    assert result_bad == 1, (
        "Pre-fix INDEX.md with nested file listed as top-level ### N. `filename` "
        "entry must produce exit code 1 (mismatch: file exists only in subfolder, "
        "not at books_library/ root). "
        "This is the RED-without-fix check for DEEP-10."
    )

    # --- Scenario that is GREEN after the DEEP-10 fix ---
    # Nested entry uses ### N-nested. filename (no backtick-filename) — HEADING_RE
    # does not match it, so parsed count == 1 == top-level fs count.
    good_index = (
        "### 1. `top_level_book.pdf`\n"
        "### 2-nested. nested_doc.pdf (in Документация ЛСРМ/)\n"
    )
    (ref_books / "INDEX.md").write_text(good_index, encoding="utf-8")

    result_good = vbi.main(["--root", str(tmp_path), "--books-dir", str(books_dir)])
    assert result_good == 0, (
        "Post-fix INDEX.md with nested entry as ### N-nested. (no backtick filename) "
        "must produce exit code 0: parsed count (1) == top-level fs count (1). "
        "Failing assertion: post-fix result must be 0."
    )
