"""
F-151 / v1.17.9 — RAG corpus extractor (PDF → JSON).

ЗАКРЕПЛЁННЫЙ ИНСТРУМЕНТ. Не запускается в обычном рабочем потоке;
вызывается по явной команде когда нужно построить полнотекстовый
корпус из PDF библиотеки знаний:

    python -m gamma.knowledge.rag_extract \\
        --books-dir references/books \\
        --out references/knowledge_corpus.json

Поведение:
    - Перебирает все *.pdf в --books-dir.
    - Для каждой страницы извлекает текст через pypdf.PdfReader.
    - Чанкит постранично + по абзацам (1500..3000 символов).
    - Распознаёт заголовки разделов (RU «§N.N» / «Глава N» /
      EN «Chapter N» / «N.N.N» в начале абзаца).
    - Сохраняет в JSON:
        {
          "version": "1.0",
          "generated_at": "<передаётся через --timestamp>",
          "books": {
              "<short_id>": {
                  "file": "Lsrm_algorithmic_foundations.pdf",
                  "n_pages": 55,
                  "chunks": [
                      {"page": 7, "section": "§8.4", "text": "…"},
                      ...
                  ]
              }
          }
        }

Совместимость с RAG search:
    rag_search.py подхватит этот файл автоматически, если он есть.
    Без него поиск работает по manually-curated knowledge_index.json
    (F-138 INDEX.md + ссылки из кода).

Замечания:
    - pypdf не идеально работает с двухколоночной вёрсткой
      (Gilmore 2008): порядок строк может смешиваться. Это
      ОК для BM25/keyword поиска, но не для семантического анализа.
    - При запуске на больших PDF (Gilmore 389 страниц)
      первичная экстракция занимает 30-60 секунд.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional

try:
    import pypdf
except ImportError:
    pypdf = None  # type: ignore


# ──────────────────────────────────────────────────────────────────
# Section heading detection
# ──────────────────────────────────────────────────────────────────

# RU/EN markers — strict patterns, applied only at line start
_SECTION_PATTERNS: List[re.Pattern[str]] = [
    re.compile(r"^\s*§\s*(\d+(?:\.\d+){0,3})\b"),                  # §8.4.4
    re.compile(r"^\s*(\d+\.\d+(?:\.\d+){0,2})\s+[A-ZА-Я]"),        # 8.4.4 Compton
    re.compile(r"^\s*(Глава|Раздел)\s+(\d+)", re.IGNORECASE),
    re.compile(r"^\s*(Chapter|Section)\s+(\d+)", re.IGNORECASE),
    re.compile(r"^\s*(Приложение|Appendix)\s+([A-ZА-Я])", re.IGNORECASE),
]

# Minimum chunk size (chars) — below this we keep merging with the next
_CHUNK_MIN = 600
# Soft target — try to break around this size
_CHUNK_SOFT = 1800
# Hard ceiling — even mid-paragraph, force-split here
_CHUNK_MAX = 3500


def _detect_section(line: str) -> Optional[str]:
    """Return canonical section marker if the line starts a new section."""
    for pat in _SECTION_PATTERNS:
        m = pat.match(line)
        if m:
            # Use the first captured group as marker label
            groups = [g for g in m.groups() if g]
            return " ".join(groups[:2])
    return None


def _chunk_page_text(text: str, page_num: int) -> List[Dict]:
    """Split a single page's text into chunks with section metadata.

    Falls back to fixed-size windows if no sections detected.
    """
    if not text or not text.strip():
        return []

    lines = text.splitlines()
    chunks: List[Dict] = []
    cur_section: Optional[str] = None
    cur_text: List[str] = []

    def _flush(section: Optional[str], buf: List[str]) -> None:
        s = "\n".join(buf).strip()
        if not s:
            return
        # Split overly large buffers on paragraph boundaries
        if len(s) <= _CHUNK_MAX:
            chunks.append({"page": page_num, "section": section, "text": s})
            return
        # Try to split on double-newline first
        parts = re.split(r"\n\s*\n", s)
        if len(parts) == 1:
            # Fall back to size-based windows
            for i in range(0, len(s), _CHUNK_SOFT):
                chunks.append({
                    "page": page_num,
                    "section": section,
                    "text": s[i:i + _CHUNK_SOFT].strip(),
                })
            return
        buf2: List[str] = []
        for part in parts:
            buf2.append(part)
            if sum(len(p) for p in buf2) >= _CHUNK_SOFT:
                chunks.append({
                    "page": page_num,
                    "section": section,
                    "text": "\n\n".join(buf2).strip(),
                })
                buf2 = []
        if buf2:
            chunks.append({
                "page": page_num,
                "section": section,
                "text": "\n\n".join(buf2).strip(),
            })

    for line in lines:
        sec = _detect_section(line)
        if sec and cur_text:
            _flush(cur_section, cur_text)
            cur_text = []
            cur_section = sec
        elif sec:
            cur_section = sec
        cur_text.append(line)

    _flush(cur_section, cur_text)

    # Merge undersized chunks with neighbours
    merged: List[Dict] = []
    for ch in chunks:
        if merged and len(merged[-1]["text"]) < _CHUNK_MIN and \
                merged[-1]["section"] == ch["section"]:
            merged[-1]["text"] += "\n" + ch["text"]
        else:
            merged.append(ch)
    return merged


def _short_id(filename: str) -> str:
    """File stem, lowercased, used as the corpus key."""
    return Path(filename).stem.lower()


def extract_book(pdf_path: Path) -> Dict:
    """Extract one PDF to a corpus dict (page-chunks)."""
    if pypdf is None:
        raise RuntimeError(
            "pypdf is not installed; cannot extract PDF corpus. "
            "Run: pip install pypdf"
        )
    reader = pypdf.PdfReader(str(pdf_path))
    chunks: List[Dict] = []
    for i, page in enumerate(reader.pages, start=1):
        try:
            text = page.extract_text() or ""
        except Exception:
            text = ""
        chunks.extend(_chunk_page_text(text, page_num=i))
    return {
        "file": pdf_path.name,
        "n_pages": len(reader.pages),
        "n_chunks": len(chunks),
        "chunks": chunks,
    }


def extract_library(books_dir: Path, timestamp: Optional[str] = None) -> Dict:
    """Extract every PDF under books_dir into a corpus JSON."""
    out: Dict = {
        "version": "1.0",
        "generated_at": timestamp or "unstamped",
        "books": {},
    }
    pdfs = sorted(books_dir.glob("*.pdf"))
    for pdf in pdfs:
        sid = _short_id(pdf.name)
        out["books"][sid] = extract_book(pdf)
    return out


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        prog="gamma.knowledge.rag_extract",
        description=(
            "F-151 — извлечение корпуса из библиотеки PDF. Запускается "
            "только по явной команде; обычная работа RAG-поиска идёт "
            "через manually-curated knowledge_index.json."
        ),
    )
    p.add_argument(
        "--books-dir",
        default="references/books",
        help="Каталог с PDF (default: references/books)",
    )
    p.add_argument(
        "--out",
        default="references/knowledge_corpus.json",
        help="Куда сохранить корпус (default: references/knowledge_corpus.json)",
    )
    p.add_argument(
        "--timestamp",
        default=None,
        help="ISO-timestamp для маркера generated_at (опционально)",
    )
    args = p.parse_args(argv)

    books_dir = Path(args.books_dir).resolve()
    if not books_dir.is_dir():
        print(f"ERROR: books-dir не найден: {books_dir}", file=sys.stderr)
        return 2

    corpus = extract_library(books_dir, timestamp=args.timestamp)
    out_path = Path(args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        json.dumps(corpus, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    total = sum(b["n_chunks"] for b in corpus["books"].values())
    print(f"OK: {len(corpus['books'])} книг, {total} chunks → {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
