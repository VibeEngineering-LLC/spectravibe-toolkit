# -*- coding: utf-8 -*-
"""F-070 W3 — Visual similarity scoring for canonical gamma-spectrum templates.

Purpose
-------
Score an unknown query spectrum against the W2 library of canonical visual
templates (24 entries, schema 0.1, see
`audit/_rag/visual_templates/SCHEMA.md`) and return the top-K matches with
a 3-level verdict label (`match` / `ambiguous` / `mismatch`).

The module is **deterministic numeric code**: no Ollama, no ML inference.

API surface (load-bearing — Agent B's W3 JSON-report wiring depends on these
exact names/shapes, do not change without an outbox heads-up):

    compute_query_vector(counts, energy_calib, dim=128) -> np.ndarray
    load_templates(geometry_class=None, index_path=...) -> list[dict]
    score_against_templates(query_vector, templates, top_k=3) -> list[dict]

Public threshold constants:

    THRESHOLD_MATCH            = 0.93   # adjusted cosine >= → "match"
    THRESHOLD_AMBIGUOUS_LOWER  = 0.85   # [0.85, 0.93) → "ambiguous"
                                        #     < 0.85   → "mismatch"
    TIER_C_DOWNWEIGHT          = 0.70   # multiplied into raw cosine for tier C
    STALE_REFERENCE_AGE_YEARS  = 15.0   # decay_age threshold for stale badge

Encoder spec (replicates W2 `feature_vector_64` bit-for-bit)
------------------------------------------------------------
- Energy grid: 0–3000 keV, 128 equal bins (Δ = 23.4375 keV/bin).
- Channel→energy: polynomial sum c_k * ch^k with the LSRM ENERGY coefficients.
- Continuum suppression: each bin v_i is replaced by max(0, v_i - 1.2 * m_i)
  where m_i is the minimum over a ±7-bin window (15-bin total).
- Log scale: log10(suppressed + 1).
- L2-normalized; the resulting unit vector is the dot-product partner for
  cosine similarity against `template["feature_vector"]["values"]`.

Threshold rationale: see `audit/_rag/visual_templates/SIMILARITY_POLICY.md`.

Anti-hallucination
------------------
- Every emitted similarity score is a deterministic function of the query
  counts, the published canonical feature vector (in the template JSON),
  and the fixed downweight/threshold constants in this module.
- `cert_reference_dates` are read verbatim from the canonical template
  provenance + (for merged templates) the constituent raw-ingest JSONs;
  no date is inferred.
- `decay_age_years` is whatever the template recorded; if absent we surface
  `None` rather than guessing.

Out of scope
------------
- Nuclide identification verdict. This module returns a *similarity*
  judgement only. Final ID combines this with F-157 isolation / ratio scores
  inside Agent B's HTML card.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np

# ---------------------------------------------------------------------------
# Public constants (also recorded in SIMILARITY_POLICY.md — keep in sync)
# ---------------------------------------------------------------------------

THRESHOLD_MATCH: float = 0.93
THRESHOLD_AMBIGUOUS_LOWER: float = 0.85
TIER_C_DOWNWEIGHT: float = 0.70
STALE_REFERENCE_AGE_YEARS: float = 15.0

FEATURE_VECTOR_DIM: int = 128
FEATURE_VECTOR_E_MAX_KEV: float = 3000.0
FEATURE_VECTOR_ENCODING: str = "peak_emphasis_log_continuum_suppressed"
CONTINUUM_SUPPRESSION_HALF_WINDOW: int = 7
CONTINUUM_SUPPRESSION_FACTOR: float = 1.2

REPO_ROOT: Path = Path(__file__).resolve().parents[2]
DEFAULT_INDEX_PATH: Path = REPO_ROOT / "audit" / "_rag" / "VISUAL_TEMPLATES_INDEX.json"


# ---------------------------------------------------------------------------
# 1. Query-vector encoder — VERBATIM port of W2 feature_vector_64
# ---------------------------------------------------------------------------

def _channel_to_energy_func(coefs: list[float]):
    """Reproduce W2 channel→energy polynomial: sum c_k * ch^k."""
    coefs = list(coefs)

    def f(ch: float) -> float:
        return sum(c * (ch ** k) for k, c in enumerate(coefs))

    return f


def compute_query_vector(
    counts: np.ndarray,
    energy_calib: dict,
    dim: int = FEATURE_VECTOR_DIM,
) -> np.ndarray:
    """Encode a query spectrum into the W2 128-bin feature vector.

    Parameters
    ----------
    counts : 1-D array-like of channel counts. Treated as `len(counts) == N`.
    energy_calib : dict with at least:
        - ``slope_keV_per_ch`` (a1) and ``offset_keV`` (a0); OR
        - ``coefficients`` — full polynomial list ``[c0, c1, c2, ...]``.
        The full-polynomial form is the W2 LSRM-native shape; the
        2-parameter form is convenience for callers that only have a
        linear calibration. If both are supplied, ``coefficients`` wins.
    dim : feature-vector length. Default 128; W2 templates are dim=128 and
        passing a different value will break cosine compatibility.

    Returns
    -------
    np.ndarray of shape ``(dim,)``, dtype float64, L2-normalized
    (``np.linalg.norm == 1`` modulo ~1e-12).

    Implementation note
    -------------------
    Internally the routine uses **pure-Python list arithmetic** mirroring
    the W2 helper (`audit/_drafts/_ollama_helpers/_session_2026-06-04/
    F-070_W2_extract_templates.py::feature_vector_64`) so that re-encoding
    the canonical raw-ingest record yields cosine ≥ 0.999 against the
    stored canonical `feature_vector.values` (regression covered by
    `tests/rag/test_visual_similarity.py::
    test_query_vector_matches_canonical_encoding_cs137`).
    """
    counts_list = list(np.asarray(counts).astype(int).tolist())
    n_bins = int(dim)
    e_max = float(FEATURE_VECTOR_E_MAX_KEV)

    coefs = energy_calib.get("coefficients")
    if coefs is None:
        offset = float(energy_calib.get("offset_keV", 0.0))
        slope = float(energy_calib["slope_keV_per_ch"])
        coefs = [offset, slope]
    coefs = list(coefs)
    ch_to_E = _channel_to_energy_func(coefs)

    bin_width = e_max / n_bins
    bins = [0.0] * n_bins
    for c, val in enumerate(counts_list):
        E = ch_to_E(c)
        if E < 0 or E >= e_max:
            continue
        b = int(E / bin_width)
        if 0 <= b < n_bins:
            bins[b] += val

    # Continuum suppression: subtract 1.2× local minimum over ±7-bin window.
    suppressed = []
    half_window = CONTINUUM_SUPPRESSION_HALF_WINDOW
    factor = CONTINUUM_SUPPRESSION_FACTOR
    for i, v in enumerate(bins):
        win_lo = max(0, i - half_window)
        win_hi = min(n_bins, i + half_window + 1)
        local_min = min(bins[win_lo:win_hi])
        suppressed.append(max(0.0, v - factor * local_min))

    log_vec = [math.log10(v + 1.0) for v in suppressed]
    norm = math.sqrt(sum(v * v for v in log_vec))
    if norm < 1e-12:
        return np.zeros(n_bins, dtype=np.float64)
    return np.asarray([v / norm for v in log_vec], dtype=np.float64)


# ---------------------------------------------------------------------------
# 2. Template loading + cert-date resolution
# ---------------------------------------------------------------------------

def _read_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def _resolve_cert_reference_dates(template: dict, repo_root: Path) -> list[str]:
    """Build the (unique, sorted) list of certificate reference dates for a template.

    Canonical templates store one `provenance.certificate_reference_date`,
    but multi-constituent templates (e.g. K-40 marinelli, Cs-137 marinelli)
    may pull from raw-ingest records with DIFFERENT passport reference
    dates. The HTML card needs the FULL set, so we walk the constituent
    raw-ingest paths and union their passport `reference_date` values for
    the template's nuclide.
    """
    nuclide = template.get("nuclide")
    prov = template.get("provenance", {}) or {}
    dates: list[str] = []
    primary = prov.get("certificate_reference_date")
    if primary:
        dates.append(primary)

    raw_paths = prov.get("constituent_raw_ingest_paths") or []
    for rel in raw_paths:
        rp = repo_root / rel
        if not rp.exists():
            continue
        try:
            raw = _read_json(rp)
        except Exception:
            # Anti-hallucination: skip silently, surface only verifiable
            # dates from on-disk files (don't synthesize).
            continue
        for entry in raw.get("passport_entries") or []:
            if entry.get("nuclide") == nuclide:
                rd = entry.get("reference_date")
                if rd:
                    dates.append(rd)

    # dedupe preserving sortable order
    return sorted(set(d for d in dates if d))


def load_templates(
    geometry_class: str | None = None,
    index_path: Path | str = DEFAULT_INDEX_PATH,
) -> list[dict]:
    """Load canonical templates from VISUAL_TEMPLATES_INDEX.

    Parameters
    ----------
    geometry_class : optional filter — one of
        {"pointlike_5cm", "marinelli_0cm", "petri_60ml",
         "denta_100ml", "denta_120ml"}. If None, all templates load.
    index_path : path to `VISUAL_TEMPLATES_INDEX.json` (default = repo
        canonical).

    Returns
    -------
    List of fully-resolved template dicts (the canonical JSON content).
    Two extra keys are injected per record for downstream convenience:

    - ``_index_entry`` : dict, the matching entry from the index registry.
    - ``cert_reference_dates`` : list[str], deduped passport reference
      dates (top-level for cheap access by the scorer).
    """
    idx_path = Path(index_path)
    idx = _read_json(idx_path)
    repo_root = idx_path.resolve().parent.parent.parent  # …/audit/_rag → repo root
    out: list[dict] = []
    for entry in idx.get("entries", []):
        if geometry_class and entry.get("geometry_class") != geometry_class:
            continue
        # F-070 W4 S0 (2026-06-05): skip templates flagged pending operator
        # review (currently denta_100ml `other_denta_100ml` fallback). Defensive
        # for future `_pending_review/` folder population in S2+.
        if entry.get("pending_review"):
            continue
        rel = entry["path"]
        tp = repo_root / rel
        if not tp.exists():
            # Try resolving as already-absolute path
            tp = Path(rel)
            if not tp.exists():
                continue
        tmpl = _read_json(tp)
        tmpl["_index_entry"] = entry
        tmpl["cert_reference_dates"] = _resolve_cert_reference_dates(tmpl, repo_root)
        out.append(tmpl)
    return out


# ---------------------------------------------------------------------------
# 3. Cosine scoring + verdict
# ---------------------------------------------------------------------------

def _verdict(adjusted_cosine: float) -> str:
    if adjusted_cosine >= THRESHOLD_MATCH:
        return "match"
    if adjusted_cosine >= THRESHOLD_AMBIGUOUS_LOWER:
        return "ambiguous"
    return "mismatch"


def _cosine(a: np.ndarray, b: np.ndarray) -> float:
    """Cosine between two real vectors. Inputs may or may not be unit-normed.

    Returns 0.0 if either has zero L2 norm (guards against degenerate
    all-zero feature vectors).
    """
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if na < 1e-12 or nb < 1e-12:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


def score_against_templates(
    query_vector: np.ndarray,
    templates: list[dict],
    top_k: int = 3,
) -> list[dict]:
    """Score `query_vector` against `templates` and return top-K matches.

    Sort key is the RAW cosine (descending). Tier-C downweight is applied
    AFTER sort, in the per-record verdict step — so a perfect Tier-C
    match can still appear at #1 even though its adjusted score will
    drop the verdict from "match" to "ambiguous".

    (Design choice: surfacing the raw ordering keeps the operator-facing
    table physically interpretable — a Cs-137 self-match should be at #1
    even if Tier-C tags it stale; the badge tells them why the verdict
    label degraded.)

    Each output entry:

        {
            "template_id":          str,
            "nuclide":              str,
            "geometry_class":       str,
            "tier":                 "A" | "B" | "C",
            "cosine_raw":           float,
            "cosine_adjusted":      float,
            "verdict":              "match" | "ambiguous" | "mismatch",
            "decay_age_years":      float | None,
            "stale_reference":      bool,
            "cert_reference_dates": list[str],
        }
    """
    if query_vector is None or len(query_vector) == 0:
        return []
    q = np.asarray(query_vector, dtype=np.float64)

    scored: list[dict[str, Any]] = []
    for tmpl in templates:
        fv = tmpl.get("feature_vector") or {}
        values = fv.get("values") or []
        if not values:
            continue
        v = np.asarray(values, dtype=np.float64)
        if v.shape[0] != q.shape[0]:
            # Dim mismatch is a programmer error, not a runtime fallback.
            # Skip silently here; tests pin the dim invariant.
            continue
        raw = _cosine(q, v)
        tier = (tmpl.get("tier") or "").upper()
        adjusted = raw * (TIER_C_DOWNWEIGHT if tier == "C" else 1.0)

        prov = tmpl.get("provenance") or {}
        decay_age = prov.get("decay_age_years")
        try:
            decay_age_f = float(decay_age) if decay_age is not None else None
        except (TypeError, ValueError):
            decay_age_f = None
        stale = bool(decay_age_f is not None and decay_age_f >= STALE_REFERENCE_AGE_YEARS)

        scored.append({
            "template_id": tmpl.get("template_id"),
            "nuclide": tmpl.get("nuclide"),
            "geometry_class": tmpl.get("geometry_class"),
            "tier": tier or None,
            "cosine_raw": raw,
            "cosine_adjusted": adjusted,
            "verdict": _verdict(adjusted),
            "decay_age_years": decay_age_f,
            "stale_reference": stale,
            "cert_reference_dates": list(tmpl.get("cert_reference_dates") or []),
        })

    scored.sort(key=lambda r: r["cosine_raw"], reverse=True)
    if top_k is None or top_k <= 0:
        return scored
    return scored[:top_k]


# ---------------------------------------------------------------------------
# 4. Convenience: end-to-end wrapper (not part of the load-bearing surface,
#    but useful for the JSON-report wiring and for tests).
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SimilarityVerdict:
    template_id: str
    nuclide: str
    geometry_class: str
    tier: str | None
    cosine_raw: float
    cosine_adjusted: float
    verdict: str
    decay_age_years: float | None
    stale_reference: bool
    cert_reference_dates: tuple[str, ...]


def score_query(
    counts: np.ndarray,
    energy_calib: dict,
    geometry_class: str | None = None,
    top_k: int = 3,
    index_path: Path | str = DEFAULT_INDEX_PATH,
) -> list[dict]:
    """One-shot: encode → load → score. Returns the same shape as
    `score_against_templates`. Provided as a thin convenience for B's
    JSON-report wiring; the three primitives above remain the load-bearing
    surface."""
    q = compute_query_vector(counts, energy_calib)
    templates = load_templates(geometry_class=geometry_class, index_path=index_path)
    return score_against_templates(q, templates, top_k=top_k)


__all__ = [
    "THRESHOLD_MATCH",
    "THRESHOLD_AMBIGUOUS_LOWER",
    "TIER_C_DOWNWEIGHT",
    "STALE_REFERENCE_AGE_YEARS",
    "FEATURE_VECTOR_DIM",
    "FEATURE_VECTOR_E_MAX_KEV",
    "FEATURE_VECTOR_ENCODING",
    "compute_query_vector",
    "load_templates",
    "score_against_templates",
    "score_query",
    "SimilarityVerdict",
]
