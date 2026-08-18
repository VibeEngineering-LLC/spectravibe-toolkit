"""
gamma.knowledge — F-151..F-153 / v1.17.9 RAG library.

Retrieval-Augmented Generation поверх библиотеки `references/books/`.

Слои:
    rag_extract.py   — извлечение текста из PDF (pypdf), постранично,
                       с распознаванием секций (RU/EN заголовки)
    rag_index.py     — chunking + BM25 индекс, JSON-сериализация
    rag_search.py    — query API + helper-паттерны (ask/explain/cite/verify)
    cli.py           — subcommand `gamma rag {query,rebuild,stats}`

Контракты:
    - F-151 (RAG index): индекс хранится в
      `references/knowledge_index.json` + `..._corpus.json`,
      пересобирается командой `gamma rag rebuild`. Все знания из
      `references/books/*.pdf` доступны для текстового поиска.

    - F-152 (RAG search API): публичный API
        `from gamma.knowledge import rag_query, rag_explain, rag_cite`
      возвращает top-k чанков с цитатой (book, page, section)
      для использования в коде/отчётах.

    - F-153 (doc-driven decisions): новые методологические решения
      (F-rules) обязаны цитировать источник из библиотеки знаний
      через `rag_cite()`. Это закрепляется в SKILL.md.
"""

__all__ = [
    "rag_query",
    "rag_explain",
    "rag_cite",
    "rag_verify",
    "rebuild_index",
    "load_index",
]

from gamma.knowledge.rag_search import (
    rag_query,
    rag_explain,
    rag_cite,
    rag_verify,
    load_index,
)
from gamma.knowledge.rag_index import rebuild_index
