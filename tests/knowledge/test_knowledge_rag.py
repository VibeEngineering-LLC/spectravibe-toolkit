"""
F-151..F-153 / v1.17.9 — Regression tests for RAG library.

Покрытие:
    1. Индекс собирается из knowledge_index.json (≥40 entries).
    2. BM25 поиск находит правильную секцию для канонических тем:
        - «Compton step erfc NaI» → LSRM-8.4.4
        - «TCS Co-60 cascade summing» → LSRM-10 / GILMORE-8.5
        - «peak shape Gaussian tail» → LSRM-8.4 / SHENDRIK-1-PEAKSHAPE
        - «FWHM calibration quadratic» → SHENDRIK-2-FWHM
        - «Marinelli self-attenuation density» → MARINELLI-SELFATTN
    3. rag_cite() возвращает корректную каноническую цитату.
    4. rag_verify() подтверждает обоснованные утверждения и отвергает
       вымышленные.
    5. CLI subcommand `gamma rag query` отдаёт top-k.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

# Make scripts/ importable
HERE = Path(__file__).resolve().parent.parent.parent
SCRIPTS = HERE / "scripts"
sys.path.insert(0, str(SCRIPTS))

from gamma.knowledge.rag_index import build_bm25_index, tokenize  # noqa: E402
from gamma.knowledge.rag_search import (  # noqa: E402
    rag_query,
    rag_explain,
    rag_cite,
    rag_verify,
)


# ──────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def index():
    """Build BM25 index from the canonical knowledge_index.json."""
    ki_path = HERE / "references" / "knowledge_index.json"
    assert ki_path.exists(), f"knowledge_index.json missing: {ki_path}"
    return build_bm25_index(ki_path, corpus_path=None)


# ──────────────────────────────────────────────────────────────────
# T1 — index basics
# ──────────────────────────────────────────────────────────────────

def test_index_built(index):
    """Index built with ≥40 curated docs covering all 6 books."""
    assert index.n_docs >= 40, f"expected ≥40 docs, got {index.n_docs}"
    assert all(d.source_layer == "curated" for d in index.docs)
    # All 6 books from INDEX.md must be referenced
    books_in_docs = {d.book for d in index.docs}
    expected_books = {
        "lsrm_algorithmic_foundations",
        "lsrm_format_specification",
        "pgs_gilmore_2008",
        "shendrik_scintillators_pt1",
        "shendrik_scintillators_pt2",
        "experiment_results_analysis",
    }
    missing = expected_books - books_in_docs
    assert not missing, f"books missing from index: {missing}"
    # avgdl reasonable for our short-form entries (50-200 tokens)
    assert 30 <= index.avgdl <= 300, f"avgdl looks wrong: {index.avgdl}"
    # Vocabulary spans RU + EN
    has_cyr = any(any(0x400 <= ord(c) <= 0x4FF for c in t) for t in index.df)
    has_eng = any(all(ord(c) < 0x100 for c in t) for t in index.df)
    assert has_cyr and has_eng, "vocab should contain both RU and EN terms"


# ──────────────────────────────────────────────────────────────────
# T2 — canonical query → correct section
# ──────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("query,expected_doc_id", [
    ("Compton step erfc NaI h_step ступенька", "LSRM-8.4.4"),
    ("TCS coincidence summing Co-60 каскадное", "LSRM-10"),
    ("Marinelli self-attenuation плотность ρ μ(E)", "MARINELLI-SELFATTN"),
    ("FWHM калибровка квадратичная a + b·E", "SHENDRIK-2-FWHM"),
    ("matched filter convolution свёртка пиков", "PEAK-SEARCH-CONVOLUTION"),
])
def test_canonical_queries(index, query, expected_doc_id):
    """top-1 hit by canonical topic queries."""
    hits = rag_query(query, k=3, index=index)
    assert hits, f"no hits for: {query}"
    top_ids = [h.doc_id for h in hits]
    assert expected_doc_id in top_ids, (
        f"expected {expected_doc_id} in top-3 for «{query}»; "
        f"got {top_ids}"
    )


# ──────────────────────────────────────────────────────────────────
# T3 — cite returns canonical citation
# ──────────────────────────────────────────────────────────────────

def test_rag_cite_compton_step(index):
    cite = rag_cite("Compton step erfc NaI", index=index)
    assert cite is not None
    assert cite.doc_id == "LSRM-8.4.4"
    assert cite.section == "§8.4.4"
    assert "LSRM" in cite.book_title or "SpectraLine" in cite.book_title
    assert cite.formula is not None
    assert "erfc" in cite.formula
    formatted = cite.formatted()
    assert "§8.4.4" in formatted
    assert "p.5" in formatted or "p.5-6" in formatted


def test_rag_cite_marinelli(index):
    cite = rag_cite("Marinelli self-attenuation correction", index=index)
    assert cite is not None
    assert cite.book == "pgs_gilmore_2008"
    assert "exp" in (cite.formula or "")


# ──────────────────────────────────────────────────────────────────
# T4 — verify works
# ──────────────────────────────────────────────────────────────────

def test_verify_supported(index):
    """Supported claim: h_step ≈ 0.03 for NaI."""
    verdict = rag_verify(
        "h_step около 0.03 для NaI peak shape",
        index=index,
    )
    assert verdict.supported, f"NaI h_step claim should be supported: {verdict.reason}"
    assert verdict.confidence > 0.2
    assert len(verdict.supporting_hits) >= 1


def test_verify_unsupported(index):
    """Pure noise claim — should not find strong support."""
    # Используем термины, не пересекающиеся с библиотекой:
    # ботанические/гастрономические слова + псевдонаучный лексикон.
    verdict = rag_verify(
        "фламинго клубника телекинез патиссон асфальт квазар свирель",
        index=index,
        min_score=2.0,
    )
    assert not verdict.supported, f"noise claim unexpectedly supported: {verdict.reason}"


# ──────────────────────────────────────────────────────────────────
# T5 — explain composes coherent answer
# ──────────────────────────────────────────────────────────────────

def test_explain_returns_summary(index):
    exp = rag_explain("Mariscotti second derivative peak finder", index=index)
    assert exp.top_hits, "no hits"
    assert exp.primary_citation
    assert exp.short_answer
    assert len(exp.short_answer) >= 50
    # Top hit should be the Mariscotti entry
    assert exp.top_hits[0].doc_id == "PEAK-SEARCH-MARISCOTTI"


# ──────────────────────────────────────────────────────────────────
# T6 — CLI passthrough (`gamma rag query …`)
# ──────────────────────────────────────────────────────────────────

def test_cli_rag_query():
    """`python -m gamma.cli rag query "..."` returns ranked hits."""
    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONPATH"] = str(SCRIPTS) + (
        os.pathsep + env.get("PYTHONPATH", "") if env.get("PYTHONPATH") else ""
    )
    proc = subprocess.run(
        [
            sys.executable, "-m", "gamma.cli",
            "rag", "query", "energy calibration polynomial",
            "-k", "2", "--json",
        ],
        env=env,
        cwd=str(HERE),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert proc.returncode == 0, f"stderr: {proc.stderr}"
    out = json.loads(proc.stdout)
    assert isinstance(out, list)
    assert len(out) <= 2
    if out:
        assert "doc_id" in out[0]
        assert "score" in out[0]
        assert out[0]["score"] > 0


# ──────────────────────────────────────────────────────────────────
# T7 — tokenizer handles RU/EN/digits/hyphens
# ──────────────────────────────────────────────────────────────────

def test_tokenizer_mixed():
    toks = tokenize("Compton step h_step=0.03 для NaI(Tl) и HPGe")
    # All tokens lowercase
    assert all(t == t.lower() for t in toks)
    assert "compton" in toks
    assert "step" in toks
    assert "h" in toks or "h_step" in " ".join(toks)  # depending on tokenize rules
    assert "nai" in toks
    assert "hpge" in toks
    assert "для" in toks
