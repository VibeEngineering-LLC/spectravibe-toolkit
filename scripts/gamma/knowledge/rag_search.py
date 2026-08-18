"""
F-152 / v1.17.9 — RAG search API.

Публичные функции (4 паттерна):

    rag_query(query, k=5)   — top-k чанков, ранжированных BM25
    rag_explain(topic)      — связный объяснительный ответ (summary +
                              top-1 хит + цитата)
    rag_cite(topic)         — каноническая цитата для использования
                              в коде/отчётах (book §section, page-range)
    rag_verify(claim)       — проверка утверждения: есть ли в
                              библиотеке подтверждающий источник

Все 4 функции работают поверх manually-curated layer (всегда есть)
+ corpus layer (если собран). Кеширование загруженного индекса —
через `functools.lru_cache`.

Контракт F-153 (doc-driven decisions):
    Новые методологические решения (F-rules) в коде должны
    цитировать источник через `rag_cite()` или явную ссылку
    на запись `knowledge_index.json`. Иначе F-rule отклоняется
    как «не методологически обоснованный».

Пример использования в коде:

    >>> from gamma.knowledge import rag_query, rag_cite
    >>> hits = rag_query("Compton step erfc NaI h_step")
    >>> for h in hits[:3]:
    ...     print(f"{h.score:.2f}  {h.book} {h.section}")
    LSRM-8.4.4  lsrm_algorithmic_foundations §8.4.4
    ...
    >>> cite = rag_cite("FWHM calibration quadratic")
    >>> print(cite.formatted())
    [Shendrik pt.2, Часть 2, p.10-40] FWHM²(E) = a + b·E + c·E²
"""
from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from gamma.knowledge.rag_index import (
    BM25Index,
    DEFAULT_BUILT_INDEX,
    DEFAULT_INDEX_JSON,
    bm25_score,
    build_bm25_index,
    load_index_from_file,
    tokenize,
)


# ──────────────────────────────────────────────────────────────────
# Result types
# ──────────────────────────────────────────────────────────────────

@dataclass
class RagHit:
    """One ranked hit returned by rag_query()."""
    doc_id: str
    score: float
    book: str
    section: str
    title: str
    text: str
    page_from: Optional[int] = None
    page_to: Optional[int] = None
    keywords: List[str] = None  # type: ignore[assignment]
    formula: Optional[str] = None
    code_citations: List[str] = None  # type: ignore[assignment]
    source_layer: str = "curated"

    def __post_init__(self):
        if self.keywords is None:
            self.keywords = []
        if self.code_citations is None:
            self.code_citations = []

    def formatted_citation(self, book_title: Optional[str] = None) -> str:
        """Human-readable citation: [Book, §section, p.X-Y]."""
        bk = book_title or self.book
        pages = ""
        if self.page_from is not None:
            if self.page_to and self.page_to != self.page_from:
                pages = f", p.{self.page_from}-{self.page_to}"
            else:
                pages = f", p.{self.page_from}"
        sec = f", {self.section}" if self.section else ""
        return f"[{bk}{sec}{pages}]"


@dataclass
class RagExplanation:
    """Structured explanation: short answer + top hits + citation."""
    topic: str
    short_answer: str
    top_hits: List[RagHit]
    primary_citation: str
    formula: Optional[str] = None


@dataclass
class RagCitation:
    """Canonical citation for use in code/reports."""
    doc_id: str
    book: str
    book_title: str
    section: str
    page_from: Optional[int]
    page_to: Optional[int]
    topic: str
    formula: Optional[str] = None

    def formatted(self) -> str:
        pages = ""
        if self.page_from is not None:
            if self.page_to and self.page_to != self.page_from:
                pages = f", p.{self.page_from}-{self.page_to}"
            else:
                pages = f", p.{self.page_from}"
        sec = f", {self.section}" if self.section else ""
        return f"[{self.book_title}{sec}{pages}] {self.topic}"


@dataclass
class RagVerdict:
    """Result of rag_verify(claim) — is the claim supported?"""
    claim: str
    supported: bool
    confidence: float       # 0..1
    supporting_hits: List[RagHit]
    reason: str


# ──────────────────────────────────────────────────────────────────
# Index loading (with caching)
# ──────────────────────────────────────────────────────────────────

@lru_cache(maxsize=4)
def _load_cached(path_str: str) -> BM25Index:
    return load_index_from_file(Path(path_str))


def _autoresolve_root() -> Path:
    """Find the gamma-spectrum-analysis repo root from this file."""
    here = Path(__file__).resolve()
    # scripts/gamma/knowledge/rag_search.py → repo root is parents[3]
    return here.parents[3]


def load_index(root: Optional[Path] = None) -> BM25Index:
    """Load the BM25 index, building it on-the-fly if not yet saved.

    Strategy:
        1. If `references/knowledge_bm25.json` exists, load it.
        2. Else, build from `references/knowledge_index.json` (no corpus)
           and return without writing — pure read-only fallback.
        3. If neither exists, raise FileNotFoundError.
    """
    if root is None:
        root = _autoresolve_root()
    built = root / DEFAULT_BUILT_INDEX
    if built.exists():
        return _load_cached(str(built))
    curated = root / DEFAULT_INDEX_JSON
    if curated.exists():
        return build_bm25_index(curated, corpus_path=None)
    raise FileNotFoundError(
        f"Ни {built}, ни {curated} не найдены. "
        "Запустите: python -m gamma.knowledge.rag_index --rebuild"
    )


# ──────────────────────────────────────────────────────────────────
# Pattern 1 — query (low-level top-k retrieval)
# ──────────────────────────────────────────────────────────────────

def rag_query(
    query: str,
    k: int = 5,
    book_filter: Optional[List[str]] = None,
    min_score: float = 0.0,
    index: Optional[BM25Index] = None,
) -> List[RagHit]:
    """Return top-k BM25-ranked hits for query.

    Parameters
    ----------
    query : str
        Free-form query (RU / EN keywords work).
    k : int
        Number of hits to return.
    book_filter : list of book short_id, optional
        Restrict search to these books.
    min_score : float
        Drop hits with BM25 score below this threshold.
    index : BM25Index, optional
        Pre-loaded index (test injection).
    """
    if index is None:
        index = load_index()
    q_tokens = tokenize(query)
    if not q_tokens:
        return []

    scored: List[Tuple[float, "IndexedDoc"]] = []  # type: ignore[name-defined]
    for d in index.docs:
        if book_filter and d.book not in book_filter:
            continue
        s = bm25_score(index, q_tokens, d)
        if s > min_score:
            scored.append((s, d))
    scored.sort(key=lambda t: t[0], reverse=True)
    scored = scored[:k]

    hits: List[RagHit] = []
    for s, d in scored:
        hits.append(RagHit(
            doc_id=d.doc_id,
            score=s,
            book=d.book,
            section=d.section,
            title=d.title,
            text=d.text,
            page_from=d.page_from,
            page_to=d.page_to,
            keywords=list(d.keywords),
            formula=d.formula,
            code_citations=list(d.code_citations),
            source_layer=d.source_layer,
        ))
    return hits


# ──────────────────────────────────────────────────────────────────
# Pattern 2 — explain (compose a short answer from top hits)
# ──────────────────────────────────────────────────────────────────

def rag_explain(
    topic: str,
    k: int = 3,
    index: Optional[BM25Index] = None,
) -> RagExplanation:
    """Return a structured explanation of a topic with citations.

    Algorithm:
        1. rag_query(topic, k)
        2. short_answer = top-1 hit's summary (first 400 chars)
        3. primary_citation = top-1 hit's canonical citation
    """
    if index is None:
        index = load_index()
    hits = rag_query(topic, k=k, index=index)
    if not hits:
        return RagExplanation(
            topic=topic,
            short_answer="(не найдено в библиотеке знаний)",
            top_hits=[],
            primary_citation="",
            formula=None,
        )
    top = hits[0]
    book_title = index.books_meta.get(top.book, {}).get("title", top.book)
    short = top.text.split("\n", 1)[-1][:500] if "\n" in top.text else top.text[:500]
    return RagExplanation(
        topic=topic,
        short_answer=short,
        top_hits=hits,
        primary_citation=top.formatted_citation(book_title=book_title),
        formula=top.formula,
    )


# ──────────────────────────────────────────────────────────────────
# Pattern 3 — cite (canonical citation for code/reports)
# ──────────────────────────────────────────────────────────────────

def rag_cite(
    topic: str,
    index: Optional[BM25Index] = None,
) -> Optional[RagCitation]:
    """Return the canonical citation for a topic (top-1 hit), or None.

    Used in code comments / F-rules / report footers to demonstrate
    methodological backing.
    """
    if index is None:
        index = load_index()
    hits = rag_query(topic, k=1, index=index)
    if not hits:
        return None
    h = hits[0]
    book_title = index.books_meta.get(h.book, {}).get("title", h.book)
    return RagCitation(
        doc_id=h.doc_id,
        book=h.book,
        book_title=book_title,
        section=h.section,
        page_from=h.page_from,
        page_to=h.page_to,
        topic=h.title,
        formula=h.formula,
    )


# ──────────────────────────────────────────────────────────────────
# Pattern 4 — verify (is this claim supported by the library?)
# ──────────────────────────────────────────────────────────────────

def rag_verify(
    claim: str,
    min_score: float = 1.5,
    k: int = 5,
    index: Optional[BM25Index] = None,
) -> RagVerdict:
    """Check whether the library contains evidence for `claim`.

    Heuristic:
        - supported = True if top-1 BM25 score ≥ min_score
        - confidence = top-1 score / (top-1 score + 2*median_of_topk)
        - supporting_hits = all hits above min_score

    Used as a guard before making strong methodological claims.
    """
    if index is None:
        index = load_index()
    hits = rag_query(claim, k=k, index=index)
    if not hits:
        return RagVerdict(
            claim=claim,
            supported=False,
            confidence=0.0,
            supporting_hits=[],
            reason="нет ни одного хита в библиотеке знаний",
        )
    supporting = [h for h in hits if h.score >= min_score]
    top_s = hits[0].score
    supported = top_s >= min_score
    if supported and len(hits) >= 2:
        median_rest = hits[len(hits) // 2].score
        conf = top_s / (top_s + 2 * max(median_rest, 0.01))
    else:
        conf = min(1.0, top_s / max(min_score, 0.01))
    if supported:
        reason = (
            f"{len(supporting)} хит(ов) выше порога {min_score}; "
            f"top-1 score={top_s:.2f}"
        )
    else:
        reason = (
            f"top-1 score={top_s:.2f} ниже порога {min_score} — "
            "слабое подтверждение"
        )
    return RagVerdict(
        claim=claim,
        supported=supported,
        confidence=conf,
        supporting_hits=supporting,
        reason=reason,
    )


# ──────────────────────────────────────────────────────────────────
# CLI passthrough — for `gamma rag query …`
# ──────────────────────────────────────────────────────────────────

def cli_query(args) -> int:
    try:
        hits = rag_query(args.query, k=args.k)
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    if args.json:
        out = [
            {
                "doc_id": h.doc_id,
                "score": h.score,
                "book": h.book,
                "section": h.section,
                "title": h.title,
                "page_from": h.page_from,
                "page_to": h.page_to,
                "formula": h.formula,
                "code_citations": h.code_citations,
                "text": h.text,
            }
            for h in hits
        ]
        print(json.dumps(out, ensure_ascii=False, indent=2))
    else:
        if not hits:
            print("Нет хитов.")
            return 0
        print(f"Top-{len(hits)} по запросу: «{args.query}»")
        print("=" * 70)
        for i, h in enumerate(hits, 1):
            print(f"\n[{i}] {h.doc_id}  score={h.score:.2f}  ({h.source_layer})")
            print(f"    Book   : {h.book}")
            print(f"    Section: {h.section}")
            if h.page_from is not None:
                p = (f"{h.page_from}-{h.page_to}"
                     if h.page_to and h.page_to != h.page_from
                     else str(h.page_from))
                print(f"    Page   : {p}")
            if h.formula:
                print(f"    Formula: {h.formula}")
            # First 300 chars
            preview = h.text[:300].replace("\n", " ")
            print(f"    {preview}…")
            if h.code_citations:
                print(f"    Code   : {h.code_citations[0]}")
    return 0


def cli_explain(args) -> int:
    try:
        exp = rag_explain(args.topic, k=args.k)
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    print(f"Тема: {exp.topic}")
    print(f"Цитата: {exp.primary_citation}")
    print()
    print(exp.short_answer)
    if exp.formula:
        print()
        print(f"Формула: {exp.formula}")
    return 0


def cli_cite(args) -> int:
    try:
        cite = rag_cite(args.topic)
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    if cite is None:
        print("(не найдено)")
        return 1
    print(cite.formatted())
    if cite.formula:
        print(f"Формула: {cite.formula}")
    return 0


def cli_verify(args) -> int:
    try:
        verdict = rag_verify(args.claim, min_score=args.min_score)
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    status = "ПОДТВЕРЖДЕНО" if verdict.supported else "не подтверждено"
    print(f"Утверждение: {verdict.claim}")
    print(f"Статус     : {status} (confidence={verdict.confidence:.2f})")
    print(f"Причина    : {verdict.reason}")
    if verdict.supporting_hits:
        print(f"\nИсточники ({len(verdict.supporting_hits)}):")
        for h in verdict.supporting_hits[:3]:
            print(f"  · {h.doc_id} ({h.book} {h.section}) score={h.score:.2f}")
    return 0 if verdict.supported else 1


def main(argv: Optional[List[str]] = None) -> int:
    import argparse
    p = argparse.ArgumentParser(
        prog="gamma.knowledge.rag_search",
        description="F-152 — RAG search API (query/explain/cite/verify)",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    pq = sub.add_parser("query", help="top-k chunks по запросу")
    pq.add_argument("query")
    pq.add_argument("-k", type=int, default=5)
    pq.add_argument("--json", action="store_true")
    pq.set_defaults(func=cli_query)

    pe = sub.add_parser("explain", help="связное объяснение темы")
    pe.add_argument("topic")
    pe.add_argument("-k", type=int, default=3)
    pe.set_defaults(func=cli_explain)

    pc = sub.add_parser("cite", help="каноническая цитата top-1")
    pc.add_argument("topic")
    pc.set_defaults(func=cli_cite)

    pv = sub.add_parser("verify", help="проверить утверждение")
    pv.add_argument("claim")
    pv.add_argument("--min-score", type=float, default=1.5)
    pv.set_defaults(func=cli_verify)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
