# -*- coding: utf-8 -*-
from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""F-256 + F-154 — Retrofit ad-hoc citation patterns → Layer 1 [RAG-ID].

Каноничный инструмент конверсии legacy-форм ссылок в стандартный Layer 1
формат проекта (см. references/REFERENCES.md §0).

Конвертирует (в Markdown/text файлах):
- `LSRM §<X>`      → `[LSRM-Algo-<X>]`
- `Gilmore §<X>`   → `[GILMORE-<X>]`
- `Будыка §<X>`    → `[BUDYKA-<X>]`
- `Shendrik §<X>`  → `[SHENDRIK-<X>]`
- `Shendrik pt.<n> гл.<m>` → `[SHENDRIK-<n>-<m>]`

Идемпотентен — повторный запуск не двоит.

Запуск:
    PYTHONIOENCODING=utf-8 python scripts/retrofit_citations_to_layer1.py \\
        KNOWN_AND_FIXED_ISSUES.md AUDIT_v2_MERGED.md ROADMAP_v1_17_8_plus.md \\
        [--dry-run]

После retrofit проверяется регулярка по citation_translator: все Layer 1
ссылки должны быть резолвуемыми (см. report-unresolved флаг).
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable, List, Tuple

# Pattern -> replacement; порядок важен (более специфичные сверху)
RETROFIT_RULES: List[Tuple[re.Pattern, str]] = [
    # Shendrik pt.<N> гл.<M> → [SHENDRIK-<N>-<M>]
    (re.compile(r"\bShendrik\s+pt\.\s*(\d+)\s+гл\.\s*(\d+)"), r"[SHENDRIK-\1-\2]"),
    # Будыка §X.Y(.Z) — Cyrillic
    (re.compile(r"\bБудыка\s+§(\d+(?:\.\d+)*)"), r"[BUDYKA-\1]"),
    # Budyka §X.Y(.Z) — Latin (used in v1.17.9.3 audit text)
    (re.compile(r"\bBudyka\s+§(\d+(?:\.\d+)*)"), r"[BUDYKA-\1]"),
    # Gilmore §X.Y(.Z) — поддерживает Ch. N тоже
    (re.compile(r"\bGilmore\s+§(\d+(?:\.\d+)*)"), r"[GILMORE-\1]"),
    (re.compile(r"\bGilmore\s+Ch\.\s*(\d+)"), r"[GILMORE-\1]"),
    # LSRM §X.Y(.Z) → LSRM-Algo (главный документ — Algorithmic Foundations)
    # ВАЖНО: не цеплять "LSRM-PREC", "LSRM-ACT" и пр. (уже в Layer 1)
    (re.compile(r"\bLSRM\s+§(\d+(?:\.\d+)*)"), r"[LSRM-Algo-\1]"),
    # Shendrik §X (без pt) → ч.2 (FWHM используется в pt.2)
    (re.compile(r"\bShendrik\s+§(\d+(?:\.\d+)*)"), r"[SHENDRIK-2-\1]"),
    # Experiment §X
    (re.compile(r"\bExperiment\s+§(\d+(?:\.\d+)*)"), r"[EXPERIMENT-\1]"),
]


def retrofit_text(text: str) -> Tuple[str, int]:
    """Применить все правила. Возвращает (новый_текст, число_замен)."""
    total = 0
    for pattern, repl in RETROFIT_RULES:
        text, n = pattern.subn(repl, text)
        total += n
    return text, total


def retrofit_file(path: Path, *, dry_run: bool = False) -> int:
    src = path.read_text(encoding="utf-8")
    out, n = retrofit_text(src)
    if n and not dry_run:
        path.write_text(out, encoding="utf-8")
    return n


def _main() -> int:
    parser = argparse.ArgumentParser(
        description="F-256+F-154 Retrofit ad-hoc citation patterns → Layer 1 [RAG-ID]."
    )
    parser.add_argument("files", nargs="+", help="Markdown/text файлы для retrofit.")
    parser.add_argument("--dry-run", action="store_true", help="Не записывать, только показать счёт.")
    args = parser.parse_args()

    total = 0
    for f in args.files:
        p = Path(f)
        if not p.exists():
            print(f"[retrofit] SKIP (not found): {p}", file=sys.stderr)
            continue
        n = retrofit_file(p, dry_run=args.dry_run)
        total += n
        print(f"[retrofit] {p}: {n} substitutions", flush=True)

    print(f"[retrofit] TOTAL: {total} substitutions"
          f"{' (dry-run, not written)' if args.dry_run else ''}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
