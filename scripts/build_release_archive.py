from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""
F-150 / v1.17.8 — Release archive builder (ЗАКРЕПЛЕНО НАВСЕГДА).
F-293 / v1.17.19 — extended exclusions (books moved out, build-cache pruned).

Канонический скрипт для упаковки релиза в `1_Version/`.

Контракт «упаковать» (фиксированный, не менять без явной директивы):

  ВКЛЮЧАЕТСЯ в архив:
    - весь рабочий код (`scripts/`, `tests/`)
    - вся документация (`*.md`, `SKILL.md`, `handoff.md`, `handoff_ru.md`,
      `KNOWN_AND_FIXED_ISSUES.md`, `ROADMAP_v1_17_8_plus.md`, `README.md`)
    - все стабилизированные demo артефакты `references/demo_contract_v1_17_2/`,
      `demo_reports/v1_17_6/`, `demo_reports/v1_17_7/`, `demo_reports/v1_17_8/`
    - `references/books/INDEX.md` (каталог библиотеки знаний — указывает
      на внешнюю папку `../../books_library/`)
    - `references/iaea_cache/`, `references/nuclide_libraries/`,
      `references/*.md` (методологические разделы 01-09)
    - `references/knowledge_bm25.json`, `references/knowledge_index.json`
      (curated RAG-индексы — единственное, что runtime читает для поиска;
      см. gamma/knowledge/rag_search.py:154-174, corpus_path=None)
    - детекторные ассеты `detectors/Gamma-1S/` (efficiency, certificates,
      reference_spectra, averaged_backgrounds, lsrm-libraries)
    - `data/{aliases,anchor_patterns,nuclides,xrf_lines}.json`

  ИСКЛЮЧАЕТСЯ из архива (расширено в v1.17.19):
    - **F-293**: вся папка `books_library/` (в корне проекта) — содержит
      рабочие PDF/PPT/PPTX/DOC/DOCX библиотеки знаний. Архивируется
      ОТДЕЛЬНО через `scripts/build_books_archive.py` →
      `1_Version/books_library/gamma-books_vYYYY-MM-DD.zip`. В `references/books/`
      остаётся только `INDEX.md` как каталог-описатель.
    - **F-293**: `references/_converted_tmp/` (build-cache PPT→PPTX, ≈40 МБ)
    - **F-293/F-333**: `references/_extracted_corpus/_pages/` (legacy локация
      rendered PNG страниц, ≈63 МБ). F-333 / v1.18.18.6 — PNG-оригиналы
      перенесены в `books_library/_corpus_pages/{lsrm_activity_2014,vartanov}/`
      (внешняя библиотека, уже исключена целиком). Текстовые `.md` экстракты
      в `_extracted_corpus/` СОХРАНЯЮТСЯ для RAG. Оба пути (`_pages` и
      `_corpus_pages`) в exclude-листе для defense-in-depth.
    - **F-293**: `references/_extracted_corpus/_converted_legacy/` и
      `_converted_odt/` (intermediate file-conversion outputs).
    - **copyright**: ВСЯ папка `references/_extracted_corpus/` (95 `.md`,
      ≈6.5 MiB) — VERBATIM текстовые экстракты копирайтных книг/ГОСТ.
      Runtime RAG читает `references/knowledge_*.json` (curated индексы,
      которые ВКЛЮЧАЮТСЯ выше), а НЕ эти `.md` — см. rag_search.py:154-174
      (corpus_path=None). Сырые экстракты в архив не идут (copyright),
      остаются только в репозитории. (Прежний docstring ошибочно числил
      их ВКЛЮЧЁННЫМИ — исправлено.)
    - **copyright**: `references/Lsrm_algorithmic_foundations.txt` —
      verbatim текст книги LSRM. Runtime использует `knowledge_*.json`,
      не этот файл. (Прежде ошибочно числился ВКЛЮЧЁННЫМ — исправлено.)
    - На-всякий-случай legacy guard: PDF/PPT/PPTX/DOC/DOCX в
      `references/books/` (если кто-то снова положит) — exclude.
    - `__pycache__/`, `.pytest_cache/`, `.pyc`/.pyo`/.pyd`
    - переходные demo дирректории `demo_reports/v1_15*/`, `_test_*/`,
      `_smoke/`, `_f145_smoke/`, `_tmp/`
    - системные `.DS_Store`, `desktop.ini`, `Thumbs.db`
    - VCS/IDE: `.git/`, `.vscode/`, `.idea/`

  АВТО-ПРОВЕРКА:
    Перед упаковкой автоматически запускается `verify_books_inventory.py`
    с режимом `--check` против `INDEX.md`. При несоответствии пишет
    предупреждение в stderr но НЕ блокирует архивацию (книги внешние).
    Для блокирующей проверки: `--strict-books-inventory`.

Использование:
    python scripts/build_release_archive.py 1.17.19
    python scripts/build_release_archive.py 1.17.19 --strict-books-inventory

  Результат: `../1_Version/SpectraVibe_v1.17.19.zip`
"""

import argparse
import os
import sys
import time
import zipfile
from pathlib import Path


# ──────────────────────────────────────────────────────────────────
# F-150 контракт исключений — ЗАКРЕПЛЕНО НАВСЕГДА
# ──────────────────────────────────────────────────────────────────

EXCLUDE_DIRS = frozenset({
    "__pycache__", ".pytest_cache",
    "_tmp",
    # переходные test-демо
    "_test_anonymize", "_test_bg_only", "_test_chain_completeness",
    "_test_cost_footer", "_test_interactive", "_test_no_en",
    "_test_pdf_artefact", "_test_sorted", "_test_th_composite",
    "_smoke", "_f145_smoke",
    # переходные релизы (v1.15.x) — обновлено через ОИСН-16 рефакторинг
    "v1_15", "v1_15_1", "v1_15_2",
    # VCS / IDE
    ".git", ".vscode", ".idea",
    # F-293 (v1.17.19) — build-cache корпусной экстракции
    "_converted_tmp",     # references/_converted_tmp/ (PPT→PPTX, ≈40 МБ)
    "_converted_legacy",  # references/_extracted_corpus/_converted_legacy/
    "_converted_odt",     # references/_extracted_corpus/_converted_odt/
    "_pages",             # references/_extracted_corpus/_pages/ (legacy; F-333 переехало в books_library/_corpus_pages/)
    "_corpus_pages",      # books_library/_corpus_pages/ — F-333 (defense-in-depth; books_library уже исключена ниже)
    # copyright — references/_extracted_corpus/*.md are VERBATIM text extracts
    # of copyrighted books/ГОСТ (95 files, ~6.5 MiB). Runtime uses
    # references/knowledge_*.json (curated), NOT these .md — verified
    # rag_search.py:154-174. Excluded from asset; raw corpus stays in repo only.
    "_extracted_corpus",
    # F-293 (v1.17.19) — внешняя библиотека книг (внутри проекта на root level)
    # Архивируется отдельно через scripts/build_books_archive.py →
    # 1_Version/books_library/gamma-books_vYYYY-MM-DD.zip
    "books_library",
    # F-115 (2026-06-05) — Gamma-1S LSRM raw spectra working copy.
    # Operator-private LSRM data (paths + cert-S/N). Visual templates
    # built from this are anonymised (provenance.constituent_raw_ingest_paths
    # → basename) before commit; raw .spe / passports never ship.
    # Defence-in-depth duplicate of .gitignore entry. See
    # detectors/Gamma-1S/README.md §7. NOTE: canonical folder = Gamma-1S
    # (НЕ Gamma-1S; кириллическая «Гамма-1с» === латинский омоглиф «Gamma-1S»,
    # user lock 2026-06-05).
    "raw_lsrm",
    # Agent mailboxes / dev-state — не входит в релиз (F-155 + briefing v1.20.0)
    "_state",
    # F-384 (v1.18.25.3) — папка демо-отчётов. Перенесена в
    # ../demo_reports/ (или путь, заданный через GAMMA_DEMO_REPORTS_DIR),
    # генерируется во время analyze --full-report, в архив релиза не
    # входит (15+ MiB на каждый референтный набор; быстро растёт).
    # Контракт: первый запуск CLI создаёт папку через
    # gamma.data.demo_reports_root.ensure_demo_reports_root().
    "demo_reports",
})

EXCLUDE_EXTS = frozenset({
    ".pyc", ".pyo", ".pyd",
    # F-150: any .zip/.7z in-tree is build output / nested release artifact —
    # never ship it (prevents the nested-zip doubling cascade). Zero legitimate
    # tracked .zip/.7z exist (verified git ls-files).
    ".zip", ".7z",
})

EXCLUDE_FILES = frozenset({
    ".DS_Store", "desktop.ini", "Thumbs.db",
    # verbatim LSRM book text; runtime uses knowledge_*.json, not this.
    "Lsrm_algorithmic_foundations.txt",
})

# Microsoft Office создаёт временные lock-файлы `~$*.docx` / `~$*.pptx`
# при открытом документе; они не имеют постоянного содержимого и часто
# locked. F-150 расширение v1.17.10 — пропускаем по префиксу.
EXCLUDE_FILENAME_PREFIXES = ("~$",)

# F-293 (v1.17.19): books physically moved to external books_library/.
# Legacy guard: если кто-то снова положит binary в references/books/ —
# исключаем (kept INDEX.md only).
BOOKS_RELATIVE_DIR = ("references", "books")
BOOKS_BINARY_EXTS = frozenset({
    ".pdf", ".djvu", ".ppt", ".pptx", ".doc", ".docx", ".odt",
    ".spe", ".efr", ".lib", ".cpt", ".cfw", ".efa", ".tc", ".tmp",
    ".cen", ".efd", ".fr3", ".src", ".rar", ".mtx", ".mdb",
    ".jpg", ".jpeg", ".png", ".bmp", ".gif", ".cdr",
    ".ini", ".cnf", ".zip", ".7z",
})


def _is_in_books_binary(file_path: Path, root: Path) -> bool:
    """True если файл — binary внутри references/books/ (любая глубина).

    F-293 v1.17.19: книги перенесены в внешний `books_library/`. В
    `references/books/` должен остаться ТОЛЬКО `INDEX.md`. Этот guard
    защищает от случайного коммита binary обратно в проект.
    """
    if file_path.suffix.lower() not in BOOKS_BINARY_EXTS:
        return False
    try:
        rel = file_path.resolve().relative_to(root.resolve())
    except ValueError:
        return False
    parts = rel.parts
    return (
        len(parts) >= 3
        and parts[0] == BOOKS_RELATIVE_DIR[0]
        and parts[1] == BOOKS_RELATIVE_DIR[1]
    )


def _sde_date_time() -> tuple[int, int, int, int, int, int] | None:
    """Return a fixed ZipInfo.date_time tuple from SOURCE_DATE_EPOCH, or None.

    SOURCE_DATE_EPOCH is the standard env var used by reproducible-builds.org
    to pin timestamps inside build artefacts.  When set to a valid UNIX epoch
    integer, all ZipInfo entries receive that timestamp so two identical builds
    produce byte-for-byte identical archives (DEEP-08).

    Returns a 6-tuple (year, month, day, hour, minute, second) in local time,
    or None if SOURCE_DATE_EPOCH is absent / invalid.
    """
    raw = os.environ.get("SOURCE_DATE_EPOCH", "").strip()
    if not raw:
        return None
    try:
        epoch = int(raw)
    except ValueError:
        return None
    t = time.gmtime(epoch)
    return (t.tm_year, t.tm_mon, t.tm_mday, t.tm_hour, t.tm_min, t.tm_sec)


def build_release_archive(version: str, root: Path, out_dir: Path) -> dict:
    """Собрать архив релиза согласно F-150 контракту.

    Parameters
    ----------
    version : str
        Номер версии без префикса "v" (например "1.17.8").
    root : Path
        Путь к рабочему каталогу `gamma-spectrum-analysis/`.
    out_dir : Path
        Каталог для архива (обычно `1_Version/`).

    Returns
    -------
    dict
        Сводка: archive_path, n_files, size_compressed_mb,
        size_uncompressed_mb, skipped_books, skipped_other_count.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    archive_name = f"SpectraVibe_v{version}.zip"
    out_path = out_dir / archive_name

    # DEEP-08: SOURCE_DATE_EPOCH pins all entry timestamps for byte-reproducible
    # archives.  None → fall back to each file's mtime (original behaviour).
    fixed_dt = _sde_date_time()

    skipped_books: list[str] = []
    skipped_other = 0
    n_files = 0
    total_raw = 0

    # Collect all qualifying paths first so we can sort them.
    # Sorted order is required for deterministic zip entry sequence (DEEP-08).
    qualifying: list[Path] = []
    for current, dirs, files in os.walk(root):
        # prune dirs in place (os.walk topdown=True default)
        dirs[:] = sorted(d for d in dirs if d not in EXCLUDE_DIRS)
        for f in sorted(files):
            if f in EXCLUDE_FILES:
                skipped_other += 1
                continue
            if any(f.startswith(pfx) for pfx in EXCLUDE_FILENAME_PREFIXES):
                skipped_other += 1
                continue
            ext = os.path.splitext(f)[1].lower()
            if ext in EXCLUDE_EXTS:
                skipped_other += 1
                continue
            fp = Path(current) / f
            if _is_in_books_binary(fp, root):
                skipped_books.append(fp.name)
                continue
            qualifying.append(fp)

    with zipfile.ZipFile(out_path, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for fp in qualifying:
            arcname = str(fp.relative_to(root.parent))
            if fixed_dt is not None:
                # Build a ZipInfo with pinned timestamp for reproducibility.
                zi = zipfile.ZipInfo(filename=arcname, date_time=fixed_dt)
                zi.compress_type = zipfile.ZIP_DEFLATED
                zf.writestr(zi, fp.read_bytes())
            else:
                zf.write(fp, arcname=arcname)
            n_files += 1
            total_raw += fp.stat().st_size

    size_comp = out_path.stat().st_size
    return {
        "archive_path": str(out_path),
        "n_files": n_files,
        "size_compressed_mb": size_comp / 1024 / 1024,
        "size_uncompressed_mb": total_raw / 1024 / 1024,
        "ratio_pct": 100 * size_comp / max(total_raw, 1),
        "skipped_books": skipped_books,
        "skipped_other_count": skipped_other,
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="F-150 / F-293 — релизный архиватор (книги/PNG/build-cache исключаются)",
    )
    p.add_argument("version", help="версия без 'v' (например '1.17.19')")
    p.add_argument(
        "--root",
        default=str(Path(__file__).resolve().parents[1]),
        help="путь к gamma-spectrum-analysis/ (default: автоопределение)",
    )
    p.add_argument(
        "--out-dir",
        default=None,
        help="каталог для архива (default: ../1_Version)",
    )
    p.add_argument(
        "--strict-books-inventory",
        action="store_true",
        help="блокировать сборку при расхождении books_library/ vs INDEX.md "
             "(default: warning без блокировки)",
    )
    p.add_argument(
        "--skip-books-inventory",
        action="store_true",
        help="пропустить проверку books_library/ (offline / CI без библиотеки)",
    )
    args = p.parse_args(argv)

    root = Path(args.root).resolve()
    if not (root / "scripts" / "gamma").exists():
        print(f"ERROR: не похоже на корень skill: {root}", file=sys.stderr)
        return 2

    # Default архивный каталог: соседний с 0_Work/, на уровне
    # «1 Скилы/». root = 0_Work/gamma-spectrum-analysis/, поэтому
    # root.parent.parent = «1 Скилы/».
    out_dir = Path(args.out_dir).resolve() if args.out_dir \
        else (root.parent.parent / "1_Version")

    # F-293 pre-check: books_library inventory vs INDEX.md
    if not args.skip_books_inventory:
        try:
            sys.path.insert(0, str(root / "scripts"))
            from verify_books_inventory import verify_books_inventory  # type: ignore

            res = verify_books_inventory(root)
            if not res["ok"]:
                msg = "WARNING (F-293): books inventory mismatch:\n" + res["report"]
                print(msg, file=sys.stderr)
                if args.strict_books_inventory:
                    print("FAIL: --strict-books-inventory задан, остановка.",
                          file=sys.stderr)
                    return 3
            else:
                print(f"F-293 books inventory OK: {res['n_files']} файлов "
                      f"в books_library/ согласованы с INDEX.md")
        except Exception as e:
            print(f"WARNING: books inventory verifier не отработал: {e}",
                  file=sys.stderr)

    print(f"F-150 / F-293 release archive build")
    print(f"  root    : {root}")
    print(f"  out_dir : {out_dir}")
    print(f"  version : {args.version}")
    print()
    info = build_release_archive(args.version, root, out_dir)
    print(f"archive    : {info['archive_path']}")
    print(f"files      : {info['n_files']}")
    print(f"compressed : {info['size_compressed_mb']:.2f} MiB")
    print(f"uncomp     : {info['size_uncompressed_mb']:.2f} MiB")
    print(f"ratio      : {info['ratio_pct']:.1f}%")
    print()
    print(f"PDFs excluded ({len(info['skipped_books'])}):")
    for b in info["skipped_books"]:
        print(f"  - {b}")
    print(f"other excluded: {info['skipped_other_count']} files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
