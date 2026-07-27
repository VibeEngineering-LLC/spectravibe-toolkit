# -*- coding: utf-8 -*-
from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""F-256 + F-154 — Sync `gost_ref_num` field in references/knowledge_index.json.

Каноничный (F-154 — навсегда сохраняемый) инструмент для поддержания
двухслойной схемы ссылок (см. references/REFERENCES.md §0):
- **Layer 1** (RAG-ID) — `entries[*].id` и `entries[*].book` в knowledge_index.json
- **Layer 2** (ГОСТ № источника) — добавляется поле `gost_ref_num` в каждую
  запись (entry) и в каждую книгу (book), чтобы reporting/citation_translator.py
  мог быстро резолвить.

Запуск:
    PYTHONIOENCODING=utf-8 python scripts/sync_knowledge_index_gost_refs.py
    # с --dry-run — только показать diff без записи
    PYTHONIOENCODING=utf-8 python scripts/sync_knowledge_index_gost_refs.py --dry-run

При каждом расширении RAG (новые книги/entries) запускать ПОВТОРНО — скрипт
идемпотентен.

Источник истины маппинга prefix → GOST-№: REFERENCES.md §0 (продублирован ниже).
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, Tuple

# Импорт мапы из citation_translator — единый источник истины
sys.path.insert(0, str(Path(__file__).resolve().parent / "gamma" / "reporting"))
from citation_translator import RAG_PREFIX_TO_GOST  # noqa: E402

# Маппинг book-key (в knowledge_index.json -> books) → GOST-№
BOOK_KEY_TO_GOST: Dict[str, int] = {
    "lsrm_algorithmic_foundations": 7,
    "lsrm_format_specification": 8,
    "pgs_gilmore_2008": 19,
    "shendrik_scintillators_pt1": 14,
    "shendrik_scintillators_pt2": 15,
    "experiment_results_analysis": 18,
    "budyka_glossary": 13,
    "budyka_textbook": 12,
    "mda_basics_ru": 17,
    "dose_rate_lsrm_2000": 6,
    "precision_measurements_ru": 10,
    "precision_measurements_ru_v2": 10,
    "lsrm_activity_counting_samples_2024": 5,
    "vartanov_practical_scint_djvu": 16,
    "lsrm_precision_methods_kuvykin_2023": 10,
    "lsrm_nuclear_materials_kuvykin_2023": 11,
    "ortec_gammavision_v9_a66": 22,
    "sklearn_mixture": 23,
}


def resolve_gost_num_for_entry(entry_id: str) -> int | None:
    """Резолвит GOST-№ для конкретного entry-id (например 'LSRM-ACT-11')."""
    best_prefix = ""
    for prefix in RAG_PREFIX_TO_GOST:
        if entry_id.startswith(prefix) and len(prefix) > len(best_prefix):
            best_prefix = prefix
    if not best_prefix:
        return None
    return RAG_PREFIX_TO_GOST[best_prefix][0]


def sync(index_path: Path, *, dry_run: bool = False) -> Tuple[int, int, list]:
    """Добавить gost_ref_num во все entries+books.

    Returns: (entries_updated, books_updated, list_of_unresolved_ids)
    """
    idx = json.loads(index_path.read_text(encoding="utf-8"))
    entries_updated = 0
    books_updated = 0
    unresolved: list[str] = []

    # Books
    for key, meta in idx.get("books", {}).items():
        target = BOOK_KEY_TO_GOST.get(key)
        if target is None:
            unresolved.append(f"BOOK:{key}")
            continue
        if meta.get("gost_ref_num") != target:
            meta["gost_ref_num"] = target
            books_updated += 1

    # Entries
    for entry in idx.get("entries", []):
        eid = entry.get("id", "")
        # Сначала пробуем через book-key (точнее), включая hybrid-keys типа
        # "lsrm_algorithmic_foundations + internal binding"
        book_key = entry.get("book", "")
        target = BOOK_KEY_TO_GOST.get(book_key)
        if target is None:
            # Fallback по book-key с prefix-match (для hybrid-source entries)
            for canonical_key, num in BOOK_KEY_TO_GOST.items():
                if book_key.startswith(canonical_key):
                    target = num
                    break
        if target is None:
            # Fallback — по entry-id префиксу
            target = resolve_gost_num_for_entry(eid)
        if target is None:
            unresolved.append(f"ENTRY:{eid}")
            continue
        if entry.get("gost_ref_num") != target:
            entry["gost_ref_num"] = target
            entries_updated += 1

    # Метаданные апдейта
    idx["gost_refs_synced_at"] = "2026-05-30"
    idx["gost_refs_sync_tool"] = "scripts/sync_knowledge_index_gost_refs.py"

    if not dry_run:
        index_path.write_text(
            json.dumps(idx, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

    return entries_updated, books_updated, unresolved


def _main() -> int:
    parser = argparse.ArgumentParser(
        description="F-256+F-154 Sync gost_ref_num in knowledge_index.json (Layer 1↔Layer 2 mapping)."
    )
    parser.add_argument(
        "--index",
        default="references/knowledge_index.json",
        help="Путь к knowledge_index.json (по умолчанию относительно cwd проекта).",
    )
    parser.add_argument("--dry-run", action="store_true", help="Не записывать, только diff-print.")
    args = parser.parse_args()

    path = Path(args.index)
    if not path.exists():
        print(f"[sync_gost_refs] FAIL: not found: {path}", file=sys.stderr)
        return 2

    n_entries, n_books, unresolved = sync(path, dry_run=args.dry_run)
    print(f"[sync_gost_refs] Books updated: {n_books}", flush=True)
    print(f"[sync_gost_refs] Entries updated: {n_entries}", flush=True)
    if unresolved:
        print(f"[sync_gost_refs] UNRESOLVED ({len(unresolved)}):", flush=True)
        for u in unresolved[:20]:
            print(f"  - {u}", flush=True)
        if len(unresolved) > 20:
            print(f"  ... +{len(unresolved) - 20} more", flush=True)
    if args.dry_run:
        print("[sync_gost_refs] (dry-run; no changes written)", flush=True)
    return 0 if not unresolved else 1


if __name__ == "__main__":
    raise SystemExit(_main())
