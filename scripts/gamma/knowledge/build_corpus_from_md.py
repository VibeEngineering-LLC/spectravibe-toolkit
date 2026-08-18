"""
Строит knowledge_corpus.json из OCR .md файлов (формат ## PAGE N).

Используется вместо rag_extract.py для документов, распознанных
через Unlimited-OCR (страницы сохранены как page_NN/result.md и
затем объединены merge_ocr_pages.py в единый .md с разметкой PAGE N).

Запуск:
    python -m gamma.knowledge.build_corpus_from_md --rebuild

Результат объединяется с knowledge_corpus.json (merge-режим):
если файл уже существует, обновляются только новые/изменённые книги.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ──────────────────────────────────────────────────────────────────
# Chunk params (same as rag_extract.py)
# ──────────────────────────────────────────────────────────────────
_CHUNK_MIN = 600
_CHUNK_SOFT = 1800
_CHUNK_MAX = 3500

_SECTION_PATTERNS: List[re.Pattern] = [
    re.compile(r"^\s*§\s*(\d+(?:\.\d+){0,3})\b"),
    re.compile(r"^\s*(\d+\.\d+(?:\.\d+){0,2})\s+[A-ZА-Я]"),
    re.compile(r"^\s*(Глава|Раздел)\s+(\d+)", re.IGNORECASE),
    re.compile(r"^\s*(Chapter|Section)\s+(\d+)", re.IGNORECASE),
    re.compile(r"^\s*(Приложение|Appendix)\s+([A-ZА-Я])", re.IGNORECASE),
]

_PAGE_HEADER = re.compile(r"^##\s+PAGE\s+(\d+)\s*$", re.MULTILINE)

ROOT = Path(__file__).resolve().parents[3]
CORPUS_PATH = ROOT / "references" / "knowledge_corpus.json"

# ──────────────────────────────────────────────────────────────────
# Документы для индексирования
# ──────────────────────────────────────────────────────────────────
CORP = ROOT / "references" / "_extracted_corpus"
REF  = ROOT / "references"

SOURCES: List[Tuple[str, Path, str]] = [
    # (short_id, md_path, human_title)
    (
        "lsrm_algorithmic_foundations_ocr",
        REF / "Lsrm_algorithmic_foundations_ocr.md",
        "Алгоритмические основы SpectraLine (UL-OCR 2026-07-04)",
    ),
    (
        "aktivnost_v_schetnyh_obrazcah",
        CORP / "Документация ЛСРМ" / "01_methodology_pdf"
              / "Активность в счетных образцах. Методика измерений на гамма-спектрометрах с использоваонием ПО СпектраЛайн.pdf.md",
        "Активность радионуклидов в счётных образцах. Методика измерений (ЛСРМ/Аспект, 2014)",
    ),
    (
        "pretsizionnye_izmerenia",
        CORP / "Документация ЛСРМ" / "02_topical_pdf" / "Прецизионные измерения.pdf.md",
        "Прецизионные измерения (ЛСРМ)",
    ),
    (
        "5_2_yadernye_materialy",
        CORP / "Документация ЛСРМ" / "02_topical_pdf"
              / "5_2_Практическая спектрометрия-ядерные материалы.pdf.md",
        "Гамма-спектрометрические методы анализа ядерных материалов (Кувыкин)",
    ),
    (
        "rukovodstvo_spektroskopista_solovev",
        CORP / "Руководство_спектроскописта_V1.05_Соловьев_2024.pdf.md",
        "Руководство спектроскописта. Практическая гамма-спектроскопия с RadiaCode-103. Соловьёв, 2024",
    ),
    (
        "spectralinexx_2_0_basic_functions",
        CORP / "spectralinexx_2.0_basic_functions_rus.pdf.md",
        "Семейство программ SpectraLine. Прецизионная обработка спектров. Описание основных функций (ООО ЛСРМ)",
    ),
]


# ──────────────────────────────────────────────────────────────────
# Section detection
# ──────────────────────────────────────────────────────────────────

def _detect_section(line: str) -> Optional[str]:
    for pat in _SECTION_PATTERNS:
        m = pat.match(line)
        if m:
            groups = [g for g in m.groups() if g]
            return " ".join(groups[:2])
    return None


def _chunk_text(text: str, page_num: int) -> List[Dict]:
    if not text or not text.strip():
        return []
    lines = text.splitlines()
    chunks: List[Dict] = []
    cur_section: Optional[str] = None
    cur_text: List[str] = []

    def _flush(section, buf):
        s = "\n".join(buf).strip()
        if not s:
            return
        if len(s) <= _CHUNK_MAX:
            chunks.append({"page": page_num, "section": section, "text": s})
            return
        parts = re.split(r"\n\s*\n", s)
        if len(parts) == 1:
            for i in range(0, len(s), _CHUNK_SOFT):
                chunks.append({"page": page_num, "section": section,
                                "text": s[i:i + _CHUNK_SOFT].strip()})
            return
        buf2: List[str] = []
        for part in parts:
            buf2.append(part)
            if sum(len(p) for p in buf2) >= _CHUNK_SOFT:
                chunks.append({"page": page_num, "section": section,
                                "text": "\n\n".join(buf2).strip()})
                buf2 = []
        if buf2:
            chunks.append({"page": page_num, "section": section,
                            "text": "\n\n".join(buf2).strip()})

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

    merged: List[Dict] = []
    for ch in chunks:
        if merged and len(merged[-1]["text"]) < _CHUNK_MIN and \
                merged[-1]["section"] == ch["section"]:
            merged[-1]["text"] += "\n" + ch["text"]
        else:
            merged.append(ch)
    return merged


# ──────────────────────────────────────────────────────────────────
# MD → corpus entry
# ──────────────────────────────────────────────────────────────────

def extract_md(md_path: Path) -> Dict:
    text = md_path.read_text(encoding="utf-8")
    # Split on ## PAGE N headers
    splits = _PAGE_HEADER.split(text)
    # splits = [pre, page_num_str, page_body, page_num_str, page_body, ...]
    chunks: List[Dict] = []
    i = 1
    while i < len(splits) - 1:
        page_num = int(splits[i])
        body = splits[i + 1]
        chunks.extend(_chunk_text(body, page_num))
        i += 2
    n_pages = max((c["page"] for c in chunks), default=0)
    return {
        "file": md_path.name,
        "n_pages": n_pages,
        "n_chunks": len(chunks),
        "source": "unlimited-ocr-md",
        "chunks": chunks,
    }


# ──────────────────────────────────────────────────────────────────
# Build / merge corpus
# ──────────────────────────────────────────────────────────────────

def build_corpus(sources: List[Tuple[str, Path, str]]) -> Dict:
    corpus: Dict = {"version": "1.1", "generated_at": "2026-07-05", "books": {}}
    for sid, md_path, title in sources:
        if not md_path.exists():
            print(f"  SKIP (нет файла): {md_path.name}", file=sys.stderr)
            continue
        entry = extract_md(md_path)
        entry["title"] = title
        corpus["books"][sid] = entry
        print(f"  {sid}: {entry['n_pages']} стр., {entry['n_chunks']} chunks")
    return corpus


def merge_corpus(existing: Dict, new: Dict) -> Dict:
    merged = dict(existing)
    merged["books"] = dict(existing.get("books", {}))
    for sid, book in new.get("books", {}).items():
        merged["books"][sid] = book
    merged["generated_at"] = new.get("generated_at", "")
    return merged


# ──────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────

def main(argv=None):
    p = argparse.ArgumentParser(description="Строит corpus из OCR .md файлов")
    p.add_argument("--rebuild", action="store_true",
                   help="Собрать/обновить knowledge_corpus.json")
    p.add_argument("--out", default=str(CORPUS_PATH),
                   help=f"Выходной файл (default: {CORPUS_PATH})")
    args = p.parse_args(argv)

    if not args.rebuild:
        p.print_help()
        return 1

    out_path = Path(args.out).resolve()
    print(f"Строю корпус из {len(SOURCES)} MD источников...")
    new_corpus = build_corpus(SOURCES)

    # Merge с существующим корпусом если есть
    if out_path.exists():
        with out_path.open(encoding="utf-8") as f:
            existing = json.load(f)
        result = merge_corpus(existing, new_corpus)
        print(f"Merge с существующим: {len(existing.get('books', {}))} → {len(result['books'])} книг")
    else:
        result = new_corpus

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    total = sum(b["n_chunks"] for b in result["books"].values())
    print(f"OK: {len(result['books'])} книг, {total} chunks → {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
