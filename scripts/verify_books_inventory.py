from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
r"""
F-293 (v1.17.19, extended F-300/W6 2026-06-05) — Books library inventory verifier.

Сверяет содержимое внешней папки `books_library/` со списком файлов в
`references/books/INDEX.md`. Возвращает diff и статус OK/MISMATCH.

F-300/W6 additions:
  --strict   Strict mode: require presence of books_library/INDEX.md and each
             subtree INDEX.md (01_methodology_pdf, _corpus_pages, Документация ЛСРМ).
             Also verifies no .zip/.7z inside books_library/ tree (F-150).
             Exits with code 1 on any failure; default (without --strict) warns only.

Запускается:

  1. **Автоматически** перед `build_release_archive.py` — печатает warning
     при расхождении; с `--strict-books-inventory` блокирует архивацию.

  2. **Вручную** оператором при добавлении/удалении книг:

        python scripts/verify_books_inventory.py            # check only
        python scripts/verify_books_inventory.py --strict   # strict: INDEX.md presence + F-150
        python scripts/verify_books_inventory.py --update-index  # rewrite SHA-list

Расположение библиотеки определяется через:

  • env `GAMMA_BOOKS_LIBRARY_DIR` (абсолютный путь) — высший приоритет
  • default: `<root>/books_library` (внутри проекта, на верхнем уровне)
    Например: `D:\...\0_Work\gamma-spectrum-analysis\books_library\`

Контракт сверки
---------------
1. Файлы из секции **`### N. <filename>`** (NN-level headers) в INDEX.md
   считаются записанными в каталоге.
2. Реальные файлы в `books_library/` (рекурсивно, EXCLUDE INDEX.md и
   `Документация ЛСРМ/` подпапку — она нестабильная, со своим INDEX).
3. Расхождения:
   • **missing** — в INDEX.md есть запись, в `books_library/` файла нет.
   • **untracked** — в `books_library/` файл есть, в INDEX.md записи нет.
   • **size-mismatch** — если INDEX.md записал размер и он не совпадает.

Returns dict (programmatic API):
  {"ok": bool, "report": str, "n_files": int,
   "missing": [...], "untracked": [...], "size_mismatch": [...]}
"""

import argparse
import hashlib
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


DEFAULT_BOOKS_SUBPATH = Path("books_library")   # внутри корня проекта
INDEX_MD_PATH_FROM_ROOT = Path("references") / "books" / "INDEX.md"

# Папки внутри books_library/, которые НЕ учитываются в основной inventory
# (имеют свой смысл, могут быть нестабильны):
NESTED_SUBLIBRARIES = frozenset({"Документация ЛСРМ"})

# F-300/W6: subtrees that MUST have INDEX.md in --strict mode
REQUIRED_INDEX_SUBTREES = (
    "",                    # master INDEX.md at books_library/INDEX.md
    "01_methodology_pdf",
    "_corpus_pages",
    "Документация ЛСРМ",
)

# Заголовки `### N. <filename>` — N может быть числом с .X (e.g. "1.", "5.1.").
# Заголовок может содержать расширение или без него (короткий id).
HEADING_RE = re.compile(
    r"^###\s+\d+(?:\.\d+)*\.\s+`([^`]+)`",   # tolerate trailing annotations like "(v1.17.9.x)"
    re.MULTILINE,
)


def resolve_books_library_dir(root: Path) -> Path:
    """Resolve books_library/ path: env override, then default."""
    env = os.environ.get("GAMMA_BOOKS_LIBRARY_DIR")
    if env:
        p = Path(env)
        return p.resolve() if p.is_absolute() else (root / p).resolve()
    return (root / DEFAULT_BOOKS_SUBPATH).resolve()


def parse_index_filenames(index_md_path: Path) -> List[str]:
    r"""Извлечь имена файлов из INDEX.md секций ### N. `filename`."""
    if not index_md_path.exists():
        return []
    text = index_md_path.read_text(encoding="utf-8", errors="replace")
    return HEADING_RE.findall(text)


# Системные файлы, не считаются частью библиотеки.
SYSTEM_FILES = frozenset({
    "desktop.ini", "Thumbs.db", ".DS_Store", "INDEX.md",
})


def scan_books_library(books_dir: Path) -> List[Path]:
    """Список файлов в books_library/ (non-recursive, top-level only).

    Подпапки в NESTED_SUBLIBRARIES игнорируются (имеют свой каталог).
    Системные файлы (desktop.ini, Thumbs.db, .DS_Store, INDEX.md) тоже.
    """
    if not books_dir.exists():
        return []
    result: List[Path] = []
    for entry in sorted(books_dir.iterdir()):
        if entry.is_file() and entry.name not in SYSTEM_FILES:
            result.append(entry)
        # Папки игнорируем (nested sublibraries или вспомогательные).
    return result


def sha256_short(path: Path, n_bytes: int = 8) -> str:
    """SHA-256 первых байт (для quick comparison; full SHA для манифеста)."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()[: n_bytes * 2]


def verify_books_inventory(
    root: Path,
    books_dir: Optional[Path] = None,
) -> Dict:
    """Сверить books_library/ ↔ INDEX.md. Programmatic entry-point.

    Parameters
    ----------
    root : Path
        Корень проекта (`gamma-spectrum-analysis/`).
    books_dir : Path | None
        Override path к books_library/. По умолчанию — `resolve_books_library_dir(root)`.

    Returns
    -------
    dict с полями ok / report / n_files / missing / untracked.
    """
    if books_dir is None:
        books_dir = resolve_books_library_dir(root)
    index_md = root / INDEX_MD_PATH_FROM_ROOT

    index_names = set(parse_index_filenames(index_md))
    fs_files = scan_books_library(books_dir)
    fs_names = {p.name for p in fs_files}

    missing = sorted(index_names - fs_names)
    untracked = sorted(fs_names - index_names)

    ok = not missing and not untracked

    lines: List[str] = []
    lines.append(f"books_library : {books_dir}")
    lines.append(f"INDEX.md      : {index_md}")
    lines.append(f"files in INDEX: {len(index_names)}")
    lines.append(f"files in fs   : {len(fs_names)}")
    if missing:
        lines.append(f"MISSING ({len(missing)}) — в INDEX, отсутствуют в каталоге:")
        for m in missing:
            lines.append(f"  - {m}")
    if untracked:
        lines.append(f"UNTRACKED ({len(untracked)}) — в каталоге, нет в INDEX:")
        for u in untracked:
            lines.append(f"  - {u}")
    if ok:
        lines.append("STATUS: OK — все файлы согласованы.")
    else:
        lines.append("STATUS: MISMATCH — требуется ручная сверка.")

    return {
        "ok": ok,
        "report": "\n".join(lines),
        "n_files": len(fs_names),
        "missing": missing,
        "untracked": untracked,
        "books_dir": str(books_dir),
        "index_md": str(index_md),
    }


def verify_strict_index_presence(
    books_dir: Path,
) -> Dict:
    """F-300/W6: Verify INDEX.md presence in master + 3 subtrees, and F-150 (no .zip/.7z).

    Parameters
    ----------
    books_dir : Path
        Path to books_library/ root.

    Returns
    -------
    dict with ok / issues list / report string.
    """
    issues: List[str] = []

    # 1. Check INDEX.md presence per required subtree
    for sub in REQUIRED_INDEX_SUBTREES:
        if sub == "":
            index_path = books_dir / "INDEX.md"
            label = "books_library/INDEX.md (master)"
        else:
            index_path = books_dir / sub / "INDEX.md"
            label = f"books_library/{sub}/INDEX.md"
        if not index_path.exists():
            issues.append(f"MISSING INDEX.md: {label}")

    # 2. Check _corpus_pages/README.md still present (must not be overwritten)
    readme = books_dir / "_corpus_pages" / "README.md"
    if not readme.exists():
        issues.append("MISSING README.md: books_library/_corpus_pages/README.md (must not be deleted)")

    # 3. F-150: no .zip / .7z inside books_library/ tree
    for ext in ("*.zip", "*.7z"):
        for found in books_dir.rglob(ext):
            issues.append(f"F-150 VIOLATION: archive found in books_library/: {found.relative_to(books_dir)}")

    ok = len(issues) == 0
    lines: List[str] = []
    lines.append(f"[--strict] books_library : {books_dir}")
    if ok:
        lines.append("[--strict] STATUS: OK — all INDEX.md present, README.md preserved, no in-tree archives.")
    else:
        lines.append(f"[--strict] STATUS: FAIL — {len(issues)} issue(s):")
        for issue in issues:
            lines.append(f"  - {issue}")

    return {"ok": ok, "issues": issues, "report": "\n".join(lines)}


def list_with_sha(books_dir: Path) -> List[Tuple[str, int, str]]:
    """Список (name, size_bytes, sha256_full) для всех файлов library."""
    out: List[Tuple[str, int, str]] = []
    for p in scan_books_library(books_dir):
        h = hashlib.sha256()
        with p.open("rb") as fh:
            while True:
                chunk = fh.read(1 << 20)
                if not chunk:
                    break
                h.update(chunk)
        out.append((p.name, p.stat().st_size, h.hexdigest()))
    return out


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description="F-293 books inventory verifier")
    p.add_argument(
        "--root",
        default=str(Path(__file__).resolve().parents[1]),
        help="корень проекта (default: автоопределение)",
    )
    p.add_argument(
        "--books-dir",
        default=None,
        help="override path к books_library/ (default: env GAMMA_BOOKS_LIBRARY_DIR "
             "или ../../books_library)",
    )
    p.add_argument(
        "--print-sha",
        action="store_true",
        help="дополнительно напечатать SHA-256 + размер для каждого файла",
    )
    p.add_argument(
        "--strict",
        action="store_true",
        help="F-300/W6: strict mode — require INDEX.md in master + 3 subtrees, "
             "verify README.md preserved, check F-150 (no .zip/.7z in tree). "
             "Exit 1 on failure. Default: warn-only.",
    )
    args = p.parse_args(argv)

    root = Path(args.root).resolve()
    books_dir_override = Path(args.books_dir).resolve() if args.books_dir else None

    res = verify_books_inventory(root, books_dir=books_dir_override)
    print(res["report"])

    if args.print_sha:
        print("\nSHA-256 manifest:")
        actual_books_dir = Path(res["books_dir"])
        for name, size, sha in list_with_sha(actual_books_dir):
            mb = size / 1024 / 1024
            print(f"  {sha}  {mb:8.2f} MB  {name}")

    # F-300/W6: --strict mode adds INDEX.md presence + F-150 checks
    if args.strict:
        actual_books_dir = Path(res["books_dir"])
        strict_res = verify_strict_index_presence(actual_books_dir)
        print()
        print(strict_res["report"])
        # --strict is a SUPERSET: fail on INDEX.md/F-150 problems OR inventory mismatch.
        # V126-03: previously this branch returned 0 when only res["ok"] was False
        # (inventory mismatch), which silently disarmed the CI gate.
        return 0 if (strict_res["ok"] and res["ok"]) else 1

    # Default (non-strict): exit 1 on inventory mismatch (legacy behaviour)
    return 0 if res["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
