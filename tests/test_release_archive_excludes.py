"""Regression tests for scripts/build_release_archive.py exclusion contract.

Covers two recently-fixed leaks (synthetic tmp tree, never the real repo):

  - nested-zip doubling (#1ii): any in-tree .zip/.7z must be excluded so the
    asset cannot embed prior release artifacts (e.g. 1_Version/.../*.zip).
  - copyright (#2): verbatim text extracts must NOT ship —
    references/_extracted_corpus/*.md and references/Lsrm_algorithmic_foundations.txt.
    Runtime RAG works off the curated references/knowledge_*.json indices
    (see gamma/knowledge/rag_search.py:154-174, corpus_path=None), which MUST
    still ship.

conftest.py adds scripts/ to sys.path, so the top-level module imports directly.
"""
from __future__ import annotations

import zipfile
from pathlib import Path

from build_release_archive import build_release_archive


def _make_fake_root(tmp_path: Path) -> Path:
    """Build a synthetic repo tree under tmp_path/repo and return its path."""
    root = tmp_path / "repo"

    def write(rel: str, text: str = "x") -> None:
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")

    # normal code — must ship
    write("main.py", "print('hi')\n")
    # in-tree nested release zip — must be EXCLUDED (#1ii)
    write("1_Version/v9.9.9/SpectraVibe_v9.9.9.zip", "PK\x03\x04 fake zip\n")
    # verbatim text extract — must be EXCLUDED (#2, _extracted_corpus dir)
    write("references/_extracted_corpus/sample.md", "verbatim book text\n")
    # already-excluded subdir (defense-in-depth) — EXCLUDED
    write("references/_extracted_corpus/_converted_legacy/x.md", "legacy\n")
    # verbatim LSRM book text — EXCLUDED (#2, EXCLUDE_FILES)
    write("references/Lsrm_algorithmic_foundations.txt", "lsrm verbatim\n")
    # curated RAG indices — must SHIP
    write("references/knowledge_bm25.json", '{"idx": "bm25"}\n')
    write("references/knowledge_index.json", '{"idx": "curated"}\n')

    return root


def _archive_members(tmp_path: Path) -> list[str]:
    root = _make_fake_root(tmp_path)
    out_dir = tmp_path / "out"
    info = build_release_archive("9.9.9", root, out_dir)
    assert info["n_files"] >= 1
    with zipfile.ZipFile(info["archive_path"]) as zf:
        # normalise separators so suffix checks are OS-independent
        return [name.replace("\\", "/") for name in zf.namelist()]


def test_nested_zip_and_7z_excluded(tmp_path):
    members = _archive_members(tmp_path)
    offenders = [m for m in members if m.lower().endswith((".zip", ".7z"))]
    assert offenders == [], f"nested archive(s) leaked into asset: {offenders}"


def test_extracted_corpus_excluded(tmp_path):
    members = _archive_members(tmp_path)
    offenders = [m for m in members if "_extracted_corpus/" in m]
    assert offenders == [], f"verbatim corpus leaked into asset: {offenders}"


def test_lsrm_verbatim_text_excluded(tmp_path):
    members = _archive_members(tmp_path)
    offenders = [m for m in members if Path(m).name == "Lsrm_algorithmic_foundations.txt"]
    assert offenders == [], f"verbatim LSRM text leaked into asset: {offenders}"


def test_curated_indices_and_code_still_ship(tmp_path):
    members = _archive_members(tmp_path)
    basenames = {Path(m).name for m in members}
    for required in ("knowledge_bm25.json", "knowledge_index.json", "main.py"):
        assert required in basenames, f"{required} missing from asset (must ship)"
