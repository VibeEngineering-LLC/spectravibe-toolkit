"""
F-151 / v1.17.9 — RAG index builder (BM25 over knowledge_index.json
+ optional knowledge_corpus.json).

Архитектура двухслойная:

1. **Manually-curated layer** (`references/knowledge_index.json`)
   - ~50-100 структурированных записей из INDEX.md + ссылок в коде.
   - Каждая запись: book / section / topic / keywords / summary /
     formula / code_citations. Высокая точность, нулевая шумность.
   - Это основной слой для RAG-поиска при обычной работе.

2. **Full-text corpus layer** (`references/knowledge_corpus.json`)
   - Опциональный, собирается на запрос через `rag_extract.py`.
   - Постраничные чанки из всех 6 PDF, ~10000-20000 chunks суммарно.
   - Подключается автоматически, если файл существует.

BM25 формула (Okapi BM25):
    score(q, d) = Σ_t∈q IDF(t) · (tf(t,d)·(k1+1)) / (tf(t,d) + k1·(1 - b + b·|d|/avgdl))

Параметры по умолчанию: k1=1.5, b=0.75 (стандартные).
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ──────────────────────────────────────────────────────────────────
# Default paths (relative to gamma-spectrum-analysis root)
# ──────────────────────────────────────────────────────────────────

DEFAULT_INDEX_JSON = "references/knowledge_index.json"
DEFAULT_CORPUS_JSON = "references/knowledge_corpus.json"
DEFAULT_BUILT_INDEX = "references/knowledge_bm25.json"

# BM25 params
_BM25_K1 = 1.5
_BM25_B = 0.75

# Tokenization: lowercase, split on non-word, keep RU + EN + digits
_TOKEN_RE = re.compile(r"[A-Za-zА-Яа-яЁё][A-Za-zА-Яа-яЁё0-9\-]*")


# ──────────────────────────────────────────────────────────────────
# Data structures
# ──────────────────────────────────────────────────────────────────

@dataclass
class IndexedDoc:
    """One BM25 document — a curated entry or a corpus chunk."""
    doc_id: str
    source_layer: str  # "curated" | "corpus"
    book: str
    section: str
    title: str
    text: str
    tokens: List[str] = field(default_factory=list)
    # Tracking metadata
    page_from: Optional[int] = None
    page_to: Optional[int] = None
    keywords: List[str] = field(default_factory=list)
    formula: Optional[str] = None
    code_citations: List[str] = field(default_factory=list)


@dataclass
class BM25Index:
    """Serializable BM25 index."""
    docs: List[IndexedDoc]
    df: Dict[str, int]            # term -> document frequency
    avgdl: float                  # average document length
    n_docs: int
    books_meta: Dict[str, Dict]   # book_short_id -> metadata

    def to_dict(self) -> Dict:
        return {
            "version": "1.0",
            "n_docs": self.n_docs,
            "avgdl": self.avgdl,
            "df": self.df,
            "docs": [
                {
                    "doc_id": d.doc_id,
                    "source_layer": d.source_layer,
                    "book": d.book,
                    "section": d.section,
                    "title": d.title,
                    "text": d.text,
                    "tokens": d.tokens,
                    "page_from": d.page_from,
                    "page_to": d.page_to,
                    "keywords": d.keywords,
                    "formula": d.formula,
                    "code_citations": d.code_citations,
                }
                for d in self.docs
            ],
            "books_meta": self.books_meta,
        }

    @classmethod
    def from_dict(cls, d: Dict) -> "BM25Index":
        docs = [
            IndexedDoc(
                doc_id=x["doc_id"],
                source_layer=x["source_layer"],
                book=x["book"],
                section=x["section"],
                title=x.get("title", ""),
                text=x["text"],
                tokens=x.get("tokens", []),
                page_from=x.get("page_from"),
                page_to=x.get("page_to"),
                keywords=x.get("keywords", []),
                formula=x.get("formula"),
                code_citations=x.get("code_citations", []),
            )
            for x in d["docs"]
        ]
        return cls(
            docs=docs,
            df=d["df"],
            avgdl=d["avgdl"],
            n_docs=d["n_docs"],
            books_meta=d.get("books_meta", {}),
        )


# ──────────────────────────────────────────────────────────────────
# Tokenization
# ──────────────────────────────────────────────────────────────────

def tokenize(text: str) -> List[str]:
    """Lowercase + extract RU/EN words (alpha + digits + hyphen)."""
    if not text:
        return []
    return [t.lower() for t in _TOKEN_RE.findall(text)]


# ──────────────────────────────────────────────────────────────────
# Building from curated knowledge_index.json
# ──────────────────────────────────────────────────────────────────

def _curated_to_doc(entry: Dict, books_meta: Dict) -> IndexedDoc:
    """Convert one curated entry into an IndexedDoc."""
    book_id = entry["book"]
    book_title = books_meta.get(book_id, {}).get("title", book_id)
    text_parts: List[str] = []
    if entry.get("topic_ru"):
        text_parts.append(entry["topic_ru"])
    if entry.get("topic_en"):
        text_parts.append(entry["topic_en"])
    if entry.get("summary_ru"):
        text_parts.append(entry["summary_ru"])
    if entry.get("formula"):
        text_parts.append("Формула: " + entry["formula"])
    if entry.get("keywords"):
        text_parts.append("Ключевые слова: " + ", ".join(entry["keywords"]))
    if entry.get("code_citations"):
        text_parts.append("Код: " + "; ".join(entry["code_citations"]))
    full_text = "\n".join(text_parts)

    title_parts = [entry.get("topic_ru") or entry.get("topic_en") or ""]
    title_parts.append(f"[{book_title}, {entry.get('section','')}]")
    title = " ".join(title_parts).strip()

    tokens = tokenize(full_text)
    # Boost keywords: count each twice
    for kw in entry.get("keywords", []):
        tokens.extend(tokenize(kw))

    return IndexedDoc(
        doc_id=entry["id"],
        source_layer="curated",
        book=book_id,
        section=entry.get("section", ""),
        title=title,
        text=full_text,
        tokens=tokens,
        page_from=entry.get("page_from"),
        page_to=entry.get("page_to"),
        keywords=entry.get("keywords", []),
        formula=entry.get("formula"),
        code_citations=entry.get("code_citations", []),
    )


def _corpus_to_docs(corpus: Dict) -> List[IndexedDoc]:
    """Flatten corpus JSON into IndexedDoc list."""
    docs: List[IndexedDoc] = []
    for book_id, book in corpus.get("books", {}).items():
        for i, chunk in enumerate(book.get("chunks", [])):
            text = chunk.get("text", "")
            section = chunk.get("section") or ""
            page = chunk.get("page", 0)
            doc_id = f"CORPUS-{book_id}-p{page}-c{i:04d}"
            docs.append(IndexedDoc(
                doc_id=doc_id,
                source_layer="corpus",
                book=book_id,
                section=section,
                title=f"[{book_id}, p.{page}, {section}]",
                text=text,
                tokens=tokenize(text),
                page_from=page,
                page_to=page,
                keywords=[],
                formula=None,
                code_citations=[],
            ))
    return docs


def build_bm25_index(
    knowledge_index_path: Path,
    corpus_path: Optional[Path] = None,
) -> BM25Index:
    """Build a BM25 index over curated + (optional) corpus layers."""
    if not knowledge_index_path.exists():
        raise FileNotFoundError(
            f"knowledge_index.json не найден: {knowledge_index_path}"
        )

    with knowledge_index_path.open("r", encoding="utf-8") as f:
        ki = json.load(f)
    books_meta = ki.get("books", {})
    entries = ki.get("entries", [])
    docs = [_curated_to_doc(e, books_meta) for e in entries]

    if corpus_path and corpus_path.exists():
        with corpus_path.open("r", encoding="utf-8") as f:
            corpus = json.load(f)
        docs.extend(_corpus_to_docs(corpus))

    # Compute df + avgdl
    df: Dict[str, int] = {}
    for d in docs:
        for t in set(d.tokens):
            df[t] = df.get(t, 0) + 1
    total_len = sum(len(d.tokens) for d in docs)
    avgdl = (total_len / len(docs)) if docs else 0.0

    return BM25Index(
        docs=docs,
        df=df,
        avgdl=avgdl,
        n_docs=len(docs),
        books_meta=books_meta,
    )


# ──────────────────────────────────────────────────────────────────
# BM25 scoring
# ──────────────────────────────────────────────────────────────────

def bm25_score(
    index: BM25Index,
    query_tokens: List[str],
    doc: IndexedDoc,
    k1: float = _BM25_K1,
    b: float = _BM25_B,
) -> float:
    """Standard Okapi BM25 score for query against one document."""
    score = 0.0
    dl = len(doc.tokens)
    if dl == 0 or index.avgdl == 0:
        return 0.0
    # Count tf in doc
    tf: Dict[str, int] = {}
    for t in doc.tokens:
        tf[t] = tf.get(t, 0) + 1
    for q in query_tokens:
        if q not in tf:
            continue
        df_q = index.df.get(q, 0)
        if df_q == 0:
            continue
        # Robertson-Spärck Jones IDF (with +1 smoothing)
        idf = math.log(1 + (index.n_docs - df_q + 0.5) / (df_q + 0.5))
        f_qd = tf[q]
        denom = f_qd + k1 * (1 - b + b * dl / index.avgdl)
        score += idf * f_qd * (k1 + 1) / denom
    return score


# ──────────────────────────────────────────────────────────────────
# Persistence
# ──────────────────────────────────────────────────────────────────

def save_index(index: BM25Index, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(index.to_dict(), ensure_ascii=False),
        encoding="utf-8",
    )


def load_index_from_file(path: Path) -> BM25Index:
    with path.open("r", encoding="utf-8") as f:
        return BM25Index.from_dict(json.load(f))


# ──────────────────────────────────────────────────────────────────
# Public rebuild API (called from gamma.knowledge.__init__)
# ──────────────────────────────────────────────────────────────────

def rebuild_index(
    root: Optional[Path] = None,
    include_corpus: bool = True,
) -> Tuple[BM25Index, Path]:
    """Rebuild the BM25 index from knowledge_index.json (+ optional corpus)
    and save to knowledge_bm25.json.

    Returns (index, output_path).
    """
    if root is None:
        # Auto-detect: walk up from this file to find references/
        here = Path(__file__).resolve()
        # scripts/gamma/knowledge/rag_index.py → repo root is parents[3]
        root = here.parents[3]
    ki_path = root / DEFAULT_INDEX_JSON
    corpus_path = root / DEFAULT_CORPUS_JSON if include_corpus else None
    out_path = root / DEFAULT_BUILT_INDEX
    index = build_bm25_index(ki_path, corpus_path)
    save_index(index, out_path)
    return index, out_path


# ──────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────

def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        prog="gamma.knowledge.rag_index",
        description="F-151 — построить BM25 индекс из knowledge_index.json",
    )
    p.add_argument("--rebuild", action="store_true",
                   help="Пересобрать индекс (BM25 over curated [+corpus])")
    p.add_argument("--no-corpus", action="store_true",
                   help="Не подключать knowledge_corpus.json даже если есть")
    p.add_argument("--root", default=None,
                   help="Корень skill (default: автоопределение)")
    p.add_argument("--stats", action="store_true",
                   help="Показать статистику индекса")
    args = p.parse_args(argv)

    root = Path(args.root).resolve() if args.root else None

    if args.rebuild:
        index, out_path = rebuild_index(
            root=root,
            include_corpus=not args.no_corpus,
        )
        n_curated = sum(1 for d in index.docs if d.source_layer == "curated")
        n_corpus = sum(1 for d in index.docs if d.source_layer == "corpus")
        print(f"OK: BM25 index → {out_path}")
        print(f"  total docs:    {index.n_docs}")
        print(f"  curated:       {n_curated}")
        print(f"  corpus chunks: {n_corpus}")
        print(f"  avgdl:         {index.avgdl:.1f} tokens")
        print(f"  vocab size:    {len(index.df)} terms")
        return 0

    if args.stats:
        from gamma.knowledge.rag_search import load_index as _load
        try:
            index = _load()
            n_curated = sum(1 for d in index.docs if d.source_layer == "curated")
            n_corpus = sum(1 for d in index.docs if d.source_layer == "corpus")
            print(f"docs:    {index.n_docs}")
            print(f"curated: {n_curated}")
            print(f"corpus:  {n_corpus}")
            print(f"avgdl:   {index.avgdl:.1f}")
            print(f"vocab:   {len(index.df)}")
            return 0
        except FileNotFoundError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 2

    p.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
