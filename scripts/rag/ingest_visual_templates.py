# -*- coding: utf-8 -*-
"""Ingest .spe files from the Gamma-1S archive into raw-ingest JSON records.

This module is the Wave-1 deliverable for F-070 (RAG visual spectrum templates).
It reads candidate .spe files from an operator-supplied corpus root, classifies
their geometry, and emits one JSON record per file into the directory structure
that Agent A will classify in Wave 2.

Invocation
----------
    python -m scripts.rag.ingest_visual_templates \\
        --corpus-root <path-to-archive> \\
        --output-root audit/_rag/visual_templates \\
        [--dry-run]

    # Or as a direct script:
    python scripts/rag/ingest_visual_templates.py --corpus-root <path> ...

Output structure
----------------
Canonical geometries::

    <output-root>/_raw_ingest/<geometry>/<template_id>.json

Drift-study spectra (Поверка-2016 Маринелли + Точка 25см)::

    <output-root>/_drift_study/_raw_ingest_poverka2016/<geometry>/<template_id>.json

Поверка-2016 Петри/Дента are canonical candidates (A's QC decides retention)::

    <output-root>/_raw_ingest/<geometry>/<template_id>.json

Block-list (skipped, no record emitted)
-----------------------------------------
- ``M_cs_тяж_2001-2005.spe``                  — heavy matrix, not canonical
- ``Фон_*.spe`` / ``фон_*.spe`` / ``bg_*.spe`` — backgrounds
- Any file under a ``Временная нестабильность`` directory — stability runs

F-115 compliance
-----------------
``source_file`` uses the ``<CORPUS>`` placeholder; operator-absolute paths
never leak into committed JSON files. ``absolute_source_path`` stores the
real path for the current session (not committed to git).
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import traceback
from pathlib import Path
from typing import Optional

# ─── Project path bootstrap (mirrors build_spectra_index.py) ─────────────────
PROJ_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJ_ROOT / "scripts"))

from gamma.io.lsrm_spe import read_lsrm_spe  # noqa: E402  (after sys.path fix)

# ─── Constants ────────────────────────────────────────────────────────────────

CORPUS_PLACEHOLDER = "<CORPUS>"
SCHEMA_VERSION = "0.1-raw"

# Canonical nuclide list (normalised form).
# Short forms found in filenames mapped to canonical.
_NUCLIDE_ALIASES: dict[str, str] = {
    "Cs137": "Cs-137",
    "K40": "K-40",
    "Ra226": "Ra-226",
    "Th232": "Th-232",
    "Th228": "Th-228",
    "Am241": "Am-241",
    "Co60": "Co-60",
    "Ba133": "Ba-133",
    "Cd109": "Cd-109",
    "Eu152": "Eu-152",
    "Y88": "Y-88",
    "Na22": "Na-22",
}

_VALID_NUCLIDES: frozenset[str] = frozenset({
    "Am-241", "Cs-137", "K-40", "Th-232", "Ra-226", "Th-228",
    "Y-88", "Eu-152", "Ba-133", "Cd-109", "Co-60", "Na-22",
    "Смесь-AmTiCsEu",
} | set(_NUCLIDE_ALIASES.keys()))

# Regex patterns for geometry classification (applied to the stringified path).
# Order matters — more-specific patterns first.
_RE_DRIFT_MARINELLI = re.compile(r"Поверка-2016.+Маринелли", re.IGNORECASE)
_RE_DRIFT_POINT25   = re.compile(r"Поверка-2016.+Точка 25см", re.IGNORECASE)
_RE_POVERKA_PETRI   = re.compile(r"Поверка-2016.+Чашка Петри 60мл", re.IGNORECASE)
_RE_POVERKA_DENTA   = re.compile(r"Поверка-2016.+Дента-100мл", re.IGNORECASE)
_RE_POINT5CM        = re.compile(r"(Точечная-5см|5cm)", re.IGNORECASE)
_RE_MARINELLI       = re.compile(r"(Маринелли|Marinelli)", re.IGNORECASE)

# Block-list patterns applied to the filename only (not the full path).
_BLOCKLIST_EXACT = frozenset({"M_cs_тяж_2001-2005.spe"})
_RE_BLOCKLIST_NAME = re.compile(r"^(Фон_|фон_|bg_)", re.IGNORECASE)

# Year extraction from filename (e.g. ``_2017_``, ``_2019_``, or trailing ``_2023``).
# Guard: only accept 4-digit tokens in the plausible range 1990–2035 (certificate
# numbers like SRC-05 must not be interpreted as years).
_RE_YEAR = re.compile(r"[_\s](\d{4})(?:[_\s]|\.spe$|$)")
_YEAR_MIN = 1990
_YEAR_MAX = 2035

# Geometry-code for template_id
_GEOMETRY_CODE: dict[str, str] = {
    "pointlike_5cm":  "POINT5CM",
    "marinelli_0cm":  "MARINELLI0CM",
    "petri_60ml":     "PETRI60ML",
    "denta_100ml":    "DENTA100ML",
    "pointlike_25cm": "POINT25CM",
}


# ─── Public classification helpers ────────────────────────────────────────────

def classify_geometry(path: Path) -> tuple[str, bool]:
    """Return ``(geometry_class, is_drift_study)`` for the given .spe path.

    ``is_drift_study=True`` means the file belongs to the Поверка-2016
    drift-study and should be routed to ``_drift_study/_raw_ingest_poverka2016/``
    rather than ``_raw_ingest/``.

    Args:
        path: Absolute or relative filesystem path to the .spe file.

    Returns:
        Tuple of (geometry_class string, is_drift_study bool).
        geometry_class is one of: pointlike_5cm, marinelli_0cm, petri_60ml,
        denta_100ml, pointlike_25cm, or "unknown".
    """
    path_str = str(path)

    # Drift-study routes (Поверка-2016 Маринелли and Точка 25см) — check BEFORE
    # the generic Маринелли classifier so they don't fall through to canonical.
    if _RE_DRIFT_MARINELLI.search(path_str):
        return "marinelli_0cm", True
    if _RE_DRIFT_POINT25.search(path_str):
        return "pointlike_25cm", True

    # Поверка-2016 Петри / Дента — canonical candidates, NOT drift study.
    if _RE_POVERKA_PETRI.search(path_str):
        return "petri_60ml", False
    if _RE_POVERKA_DENTA.search(path_str):
        return "denta_100ml", False

    # Generic Точечная-5см / 5cm
    if _RE_POINT5CM.search(path_str):
        return "pointlike_5cm", False

    # Generic Маринелли
    if _RE_MARINELLI.search(path_str):
        return "marinelli_0cm", False

    return "unknown", False


def extract_nuclide(filename: str) -> Optional[str]:
    """Extract and normalise the nuclide from a .spe filename.

    Extraction strategy:
    1. Strip the ``.spe`` suffix.
    2. Special-case: ``Смесь_AmTiCsEu`` prefix → ``Смесь-AmTiCsEu``.
    3. Split on ``_``; the **first token** is the nuclide candidate.
    4. Normalise via ``_NUCLIDE_ALIASES``; validate against ``_VALID_NUCLIDES``.
    5. Return ``None`` if no match.

    Cyrillic-prefixed files (e.g. ``РИСН №379_Am-Ti-Eu-Cs_...``) that have
    a non-nuclide first token return ``None`` and are skipped by the ingestor.

    Args:
        filename: Bare filename (e.g. ``Am-241_045_02_2019_Точечная-5см_5cm.spe``).

    Returns:
        Canonical nuclide string or ``None``.
    """
    stem = filename
    if stem.lower().endswith(".spe"):
        stem = stem[:-4]

    # Special: Смесь_AmTiCsEu
    if stem.startswith("Смесь_"):
        suffix = stem[len("Смесь_"):].split("_")[0]
        return f"Смесь-{suffix}" if suffix else "Смесь-AmTiCsEu"

    tokens = stem.split("_")
    if not tokens:
        return None

    first = tokens[0]
    # Normalise alias (e.g. Cs137 → Cs-137)
    canonical = _NUCLIDE_ALIASES.get(first, first)
    if canonical in _VALID_NUCLIDES:
        return canonical

    return None


def is_blocklisted(path: Path) -> bool:
    """Return ``True`` if the file must be skipped (block-list enforced).

    Three block rules (F-070 §12 Risks):
    1. Exact filename ``M_cs_тяж_2001-2005.spe``.
    2. Filename starts with ``Фон_``, ``фон_``, or ``bg_``.
    3. Any ancestor directory named ``Временная нестабильность``.

    Args:
        path: Path to the .spe file.

    Returns:
        True if the file should be skipped.
    """
    name = path.name

    if name in _BLOCKLIST_EXACT:
        return True

    if _RE_BLOCKLIST_NAME.match(name):
        return True

    # Check every ancestor component for the stability-run directory name.
    for part in path.parts[:-1]:  # exclude filename itself
        if "Временная нестабильность" in part:
            return True

    return False


def extract_year(filename: str) -> Optional[str]:
    """Extract a 4-digit year from the filename.

    Looks for patterns like ``_2017_``, ``_2019.spe``, or ``_2023`` in the
    filename. Returns the first match as a string, or ``None``.

    Args:
        filename: Bare filename.

    Returns:
        Year string (e.g. ``"2017"``) or ``None``.
    """
    for m in _RE_YEAR.finditer(filename):
        candidate = m.group(1)
        if _YEAR_MIN <= int(candidate) <= _YEAR_MAX:
            return candidate
    return None


def make_template_id(nuclide: str, geometry: str, year: Optional[str]) -> str:
    """Generate the ``VT-{NUCLIDE}-{GEOMETRY}-{YEAR}`` template ID.

    Normalisation:
    - Nuclide: remove hyphens, uppercase → ``Cs-137`` → ``CS137``.
    - Geometry: map via ``_GEOMETRY_CODE`` → ``marinelli_0cm`` → ``MARINELLI0CM``.
    - Year: ``"2017"`` or ``"NONE"`` when not extracted.

    Args:
        nuclide: Canonical nuclide string.
        geometry: Geometry class string.
        year: 4-digit year string or ``None``.

    Returns:
        Template ID string, e.g. ``VT-CS137-MARINELLI0CM-2017``.
    """
    nuc_code = re.sub(r"[-–]", "", nuclide).upper()
    geom_code = _GEOMETRY_CODE.get(geometry, geometry.upper())
    yr = year if year else "NONE"
    return f"VT-{nuc_code}-{geom_code}-{yr}"


# ─── Core ingest ──────────────────────────────────────────────────────────────

def ingest_one(path: Path, corpus_root: Path) -> Optional[dict]:
    """Read a single .spe file and return a raw-ingest JSON record.

    Returns ``None`` if the file is block-listed, the nuclide cannot be
    extracted, or the geometry is unknown. Parse errors are caught and
    re-raised so the caller can log them.

    The returned dict matches the ``schema_version: "0.1-raw"`` shape::

        {
          "template_id":          str,
          "source_file":          str,   # <CORPUS>/... placeholder
          "absolute_source_path": str,   # verbatim for anti-hallucination
          "geometry_class":       str,
          "nuclide":              str,
          "source_epoch":         str | null,
          "raw_channels":         list[int],
          "n_channels":           int,
          "live_time_s":          float | null,
          "real_time_s":          float | null,
          "dead_time_pct":        float | null,
          "energy_calibration":   dict | null,
          "fwhm_calibration":     dict | null,
          "spe_description":      str,
          "schema_version":       "0.1-raw"
        }

    Args:
        path: Absolute path to the .spe file.
        corpus_root: Root used to compute the ``<CORPUS>``-prefixed relative path.

    Returns:
        Record dict or ``None`` if skipped.
    """
    if is_blocklisted(path):
        return None

    geometry, is_drift_study = classify_geometry(path)
    if geometry == "unknown":
        return None

    nuclide = extract_nuclide(path.name)
    if nuclide is None:
        return None

    year = extract_year(path.name)
    template_id = make_template_id(nuclide, geometry, year)

    # Read the .spe file (delegate entirely to existing parser — do NOT rewrite).
    spec = read_lsrm_spe(str(path))

    # Build relative source_file with placeholder (F-115).
    try:
        rel = path.relative_to(corpus_root)
        source_file = f"{CORPUS_PLACEHOLDER}/{rel.as_posix()}"
    except ValueError:
        # Path not under corpus_root — use filename only.
        source_file = f"{CORPUS_PLACEHOLDER}/{path.name}"

    # Energy calibration dict (a0, a1, a2 from the stored polynomial).
    energy_cal: Optional[dict] = None
    if spec.energy_cal:
        coefs = list(spec.energy_cal)
        energy_cal = {
            f"a{i}": float(c) for i, c in enumerate(coefs)
        }

    # FWHM calibration (coefficients + model label).
    fwhm_cal: Optional[dict] = None
    if spec.stored_fwhm_calibration is not None:
        fwhm_cal = {
            "coefficients": [float(c) for c in spec.stored_fwhm_calibration.coefficients],
            "model": spec.stored_fwhm_calibration.model,
        }

    # Dead-time percentage.
    live = float(spec.live_time) if spec.live_time is not None else None
    real = float(spec.real_time) if spec.real_time is not None else None
    dead_pct: Optional[float] = None
    if live is not None and real is not None and real > 0:
        dead_pct = round((real - live) / real * 100.0, 4)

    record: dict = {
        "template_id":          template_id,
        "source_file":          source_file,
        "absolute_source_path": str(path.resolve()),
        "geometry_class":       geometry,
        "is_drift_study":       is_drift_study,
        "nuclide":              nuclide,
        "source_epoch":         year,
        "raw_channels":         spec.counts.tolist(),
        "n_channels":           int(spec.n_channels),
        "live_time_s":          live,
        "real_time_s":          real,
        "dead_time_pct":        dead_pct,
        "energy_calibration":   energy_cal,
        "fwhm_calibration":     fwhm_cal,
        "spe_description":      (spec.comments or "").strip(),
        "schema_version":       SCHEMA_VERSION,
    }
    return record


def _output_path_for_record(record: dict, output_root: Path) -> Path:
    """Compute the output JSON path for a record.

    Canonical geometries::

        <output_root>/_raw_ingest/<geometry>/<template_id>.json

    Drift-study spectra::

        <output_root>/_drift_study/_raw_ingest_poverka2016/<geometry>/<template_id>.json

    Args:
        record: Ingested record dict (must have ``is_drift_study``,
                ``geometry_class``, ``template_id``).
        output_root: Root of the visual_templates output tree.

    Returns:
        Absolute output Path (not yet created).
    """
    tid = record["template_id"]
    geom = record["geometry_class"]
    if record.get("is_drift_study"):
        return output_root / "_drift_study" / "_raw_ingest_poverka2016" / geom / f"{tid}.json"
    return output_root / "_raw_ingest" / geom / f"{tid}.json"


def ingest_corpus(
    corpus_root: Path,
    output_root: Path,
    *,
    dry_run: bool = False,
    verbose: bool = True,
) -> dict:
    """Scan *corpus_root* recursively for .spe files and ingest each one.

    Args:
        corpus_root: Root directory to scan (the archive directory).
        output_root: Root of the ``audit/_rag/visual_templates`` tree.
        dry_run: When ``True``, log what would be written without creating files.
        verbose: Emit progress lines to stderr.

    Returns:
        Summary dict with keys ``total``, ``ingested``, ``skipped``,
        ``errors``, ``records`` (list of dicts).
    """
    spe_files = sorted(corpus_root.rglob("*.spe"))
    if verbose:
        print(f"[ingest_visual_templates] Found {len(spe_files)} .spe files under {corpus_root}",
              file=sys.stderr)

    ingested: list[dict] = []
    skipped = 0
    errors: list[str] = []

    for path in spe_files:
        try:
            rec = ingest_one(path, corpus_root)
        except Exception as exc:
            errors.append(f"{path}: {type(exc).__name__}: {exc}")
            if verbose:
                print(f"[ingest_visual_templates] ERROR {path.name}: {exc}", file=sys.stderr)
            continue

        if rec is None:
            skipped += 1
            if verbose:
                print(f"[ingest_visual_templates] SKIP {path.name}", file=sys.stderr)
            continue

        out_path = _output_path_for_record(rec, output_root)
        if dry_run:
            if verbose:
                print(f"[ingest_visual_templates] DRY-RUN would write: {out_path}", file=sys.stderr)
        else:
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(
                json.dumps(rec, ensure_ascii=False, indent=2),
                encoding="utf-8",
            )
            if verbose:
                print(f"[ingest_visual_templates] WROTE {out_path.name}", file=sys.stderr)

        ingested.append(rec)

    summary = {
        "total":    len(spe_files),
        "ingested": len(ingested),
        "skipped":  skipped,
        "errors":   len(errors),
        "error_details": errors,
        "records":  ingested,
    }
    if verbose:
        print(
            f"[ingest_visual_templates] DONE — "
            f"ingested={len(ingested)} skipped={skipped} errors={len(errors)}",
            file=sys.stderr,
        )
    return summary


# ─── CLI ──────────────────────────────────────────────────────────────────────

def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Ingest Gamma-1S archive .spe files into raw-ingest JSON records "
            "for RAG visual-template classification (F-070 Wave 1)."
        )
    )
    p.add_argument(
        "--corpus-root", required=True,
        help="Root directory to scan (e.g. detectors/Gamma-1S/reference_spectra/archive). "
             "NOT stored in the output — replaced by <CORPUS> placeholder per F-115.",
    )
    p.add_argument(
        "--output-root",
        default=str(PROJ_ROOT / "audit" / "_rag" / "visual_templates"),
        help="Root of the visual_templates output tree. "
             "Default: <repo>/audit/_rag/visual_templates.",
    )
    p.add_argument(
        "--dry-run", action="store_true",
        help="Log intended writes without creating any files.",
    )
    p.add_argument(
        "--quiet", action="store_true",
        help="Suppress per-file progress lines on stderr.",
    )
    return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = _parse_args(argv)
    corpus_root = Path(args.corpus_root).resolve()
    output_root = Path(args.output_root).resolve()

    if not corpus_root.exists():
        print(f"[ingest_visual_templates] Corpus root does not exist: {corpus_root}",
              file=sys.stderr)
        return 1

    summary = ingest_corpus(
        corpus_root=corpus_root,
        output_root=output_root,
        dry_run=args.dry_run,
        verbose=not args.quiet,
    )

    print(
        f"[ingest_visual_templates] Summary: "
        f"total={summary['total']} ingested={summary['ingested']} "
        f"skipped={summary['skipped']} errors={summary['errors']}",
        file=sys.stderr,
    )
    return 0 if summary["errors"] == 0 else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        traceback.print_exc()
        sys.exit(1)
