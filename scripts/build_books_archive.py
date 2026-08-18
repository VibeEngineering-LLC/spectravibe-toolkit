from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""
F-293 (v1.17.19) — Books library archive builder.

Standalone packager для внешней библиотеки книг (`books_library/`).
**НЕ** запускается автоматически — только по явному команде оператора:

    python scripts/build_books_archive.py
    python scripts/build_books_archive.py --tag 2026-05-30
    python scripts/build_books_archive.py --books-dir /path/to/library

По умолчанию архив попадает в `1_Version/books_library/`:

    1_Version/books_library/gamma-books_vYYYY-MM-DD.zip

Где `YYYY-MM-DD` — текущая дата (или `--tag` override). Это отдельный
подкаталог `1_Version/books_library/`, чтобы не смешивать с релизными
zip'ами проекта (`SpectraVibe_vX.Y.Z.zip`).

Контракт упаковки:
  - Архивируется весь `books_library/` (включая `Документация ЛСРМ/`
    sublibrary, если присутствует), вместе с `INDEX.md` из проекта
    как `INDEX.md` в корне zip.
  - SHA-256 манифест добавляется как `MANIFEST.sha256` в корне zip
    (для повторной проверки целостности).
  - Compresslevel = 6 (баланс времени/размера; PDF плохо жмутся).

Логика обновления:
  - Архив пишется атомарно через `.tmp` → rename.
  - Если архив с таким именем уже существует — переименовываем
    предыдущий в `gamma-books_vYYYY-MM-DD.prev.zip` (один rotation).
"""

import argparse
import hashlib
import os
import shutil
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import List


# Импорт verify_books_inventory для resolve_books_library_dir / list_with_sha
SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from verify_books_inventory import (
    resolve_books_library_dir,
    list_with_sha,
    INDEX_MD_PATH_FROM_ROOT,
    NESTED_SUBLIBRARIES,
)


# F-333 / v1.18.18.6 — папки внутри books_library/, которые НЕ архивируются
# (regenerable артефакты, не часть документной библиотеки):
#   - `_corpus_pages/` — rendered PNG страниц для multimodal-чтения,
#     генерируются `_render_scanned.py` / `_ocr_*.py` за секунды.
ARCHIVE_EXCLUDE_DIRS = frozenset({"_corpus_pages"})


def _today_tag() -> str:
    """YYYY-MM-DD в UTC (стабильный во всех timezone)."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%d")


def build_books_archive(
    root: Path,
    out_dir: Path,
    books_dir: Path,
    tag: str,
) -> dict:
    """Собрать books_library/ в один zip + INDEX.md + manifest.

    Returns dict с архивной статистикой.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    archive_name = f"gamma-books_v{tag}.zip"
    out_path = out_dir / archive_name
    tmp_path = out_path.with_suffix(".zip.tmp")

    if not books_dir.exists():
        raise FileNotFoundError(f"books_library/ не найдена: {books_dir}")

    n_files = 0
    total_raw = 0

    with zipfile.ZipFile(
        tmp_path, "w", zipfile.ZIP_DEFLATED, compresslevel=6
    ) as zf:
        # 1. INDEX.md из проекта (если есть) — как INDEX.md в корне zip
        index_md = root / INDEX_MD_PATH_FROM_ROOT
        if index_md.exists():
            zf.write(index_md, arcname="INDEX.md")
            n_files += 1
            total_raw += index_md.stat().st_size

        # 2. Все файлы books_library/ (recursive), плюс nested sublibraries.
        # F-333: ARCHIVE_EXCLUDE_DIRS отсекает regenerable артефакты
        # (например `_corpus_pages/` — rendered PNG для multimodal-чтения).
        for current, dirs, files in os.walk(books_dir):
            # In-place prune excluded subdirs so os.walk не входит в них.
            dirs[:] = [d for d in dirs if d not in ARCHIVE_EXCLUDE_DIRS]
            # Архивируем nested sublibraries как отдельные подпапки
            # (Документация ЛСРМ/ итд) — это OK, размер всё равно
            # большой, но это и есть назначение этого архива.
            cur_path = Path(current)
            for f in files:
                fp = cur_path / f
                # Внутри-папочного INDEX.md тоже сохраняем (если есть)
                rel = fp.relative_to(books_dir)
                arcname = Path("books_library") / rel
                zf.write(fp, arcname=str(arcname))
                n_files += 1
                total_raw += fp.stat().st_size

        # 3. MANIFEST.sha256 — только для top-level files (не nested sublibs)
        manifest_lines: List[str] = []
        manifest_lines.append(f"# gamma-books v{tag} — MANIFEST.sha256")
        manifest_lines.append(f"# Generated: {datetime.now(timezone.utc).isoformat()}")
        manifest_lines.append(f"# Source: {books_dir}")
        manifest_lines.append("")
        for name, size, sha in list_with_sha(books_dir):
            mb = size / 1024 / 1024
            manifest_lines.append(f"{sha}  {size:>14d}  {mb:>10.2f}MB  {name}")
        zf.writestr("MANIFEST.sha256", "\n".join(manifest_lines) + "\n")
        n_files += 1

    # Rotation: если архив с тем же именем уже существует — backup
    if out_path.exists():
        backup_path = out_path.with_name(out_path.stem + ".prev.zip")
        if backup_path.exists():
            backup_path.unlink()
        out_path.rename(backup_path)

    tmp_path.rename(out_path)

    size_comp = out_path.stat().st_size
    return {
        "archive_path": str(out_path),
        "n_files": n_files,
        "size_compressed_mb": size_comp / 1024 / 1024,
        "size_uncompressed_mb": total_raw / 1024 / 1024,
        "ratio_pct": 100 * size_comp / max(total_raw, 1),
    }


def main(argv: List[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="F-293 — упаковка внешней books_library/ в отдельный zip",
    )
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
        "--out-dir",
        default=None,
        help="каталог для архива (default: ../1_Version/books_library)",
    )
    p.add_argument(
        "--tag",
        default=None,
        help="версия архива (default: сегодняшняя дата YYYY-MM-DD UTC)",
    )
    args = p.parse_args(argv)

    root = Path(args.root).resolve()
    books_dir = (
        Path(args.books_dir).resolve()
        if args.books_dir
        else resolve_books_library_dir(root)
    )
    out_dir = (
        Path(args.out_dir).resolve()
        if args.out_dir
        else (root.parent.parent / "1_Version" / "books_library")
    )
    tag = args.tag or _today_tag()

    print(f"F-293 books library archive build")
    print(f"  root      : {root}")
    print(f"  books_dir : {books_dir}")
    print(f"  out_dir   : {out_dir}")
    print(f"  tag       : {tag}")
    print()

    if not books_dir.exists():
        print(f"ERROR: books_library/ не найдена: {books_dir}", file=sys.stderr)
        return 2

    info = build_books_archive(root, out_dir, books_dir, tag)
    print(f"archive    : {info['archive_path']}")
    print(f"files      : {info['n_files']}")
    print(f"compressed : {info['size_compressed_mb']:.2f} MiB")
    print(f"uncomp     : {info['size_uncompressed_mb']:.2f} MiB")
    print(f"ratio      : {info['ratio_pct']:.1f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
