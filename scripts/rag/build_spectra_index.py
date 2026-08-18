# -*- coding: utf-8 -*-
"""Build (or refresh) ``audit/_rag/SPECTRA_INDEX.json``.

A2 design proposal: ``_state/agent_a/outbox/2026-06-04_spectra_rag_design.md``.
Prototype: ``audit/_drafts/_ollama_helpers/2026-06-04_lsrm_batch_eval.py``.

This module productionises the prototype with:

1. **Stable IDs**: ``RAG-SPEC-NNNN`` assigned in deterministic sha256-ASC
   sort order. Re-runs preserve existing IDs (looked up by sha256) and
   only assign new IDs to genuinely-new content.

2. **Byte-identical re-runs**: the schema payload (everything under the
   top-level keys excluding ``_meta.generated_at_utc``) is deterministic
   across invocations. Tests assert this.

3. **F-115 compliance**: all paths use the ``<LSRM>`` placeholder. The
   real operator root is passed via ``--lsrm-root`` (no default, no
   hardcoded path leaks into git).

4. **31-worker multiprocessing.Pool** for the ``.spe`` re-scan (Phase 1
   MAXIMUM-parallelism mandate).

CLI::

    python scripts/rag/build_spectra_index.py \\
        --lsrm-root C:\\LSRM \\
        --output audit/_rag/SPECTRA_INDEX.json

    python scripts/rag/build_spectra_index.py \\
        --lsrm-root C:\\LSRM --check        # dry run, reports only

Exit code 0 on success, 1 on hard failure.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
import traceback
from datetime import datetime
from multiprocessing import Pool, freeze_support
from pathlib import Path
from typing import Any, Optional

PROJ_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJ_ROOT / "scripts"))

LSRM_PLACEHOLDER = "<LSRM>"

# ─── Scope sets (mirror the prototype's verified P0-P3 windows) ──────────────
# Each tuple: (subdir_relative_to_lsrm_root, glob_pattern, detector_tag, priority)
P0_GLOBS = [
    (r"Work\BG\Gamma-1S", "**/*.spe", "Gamma-1S", "P0"),
]
P1_GLOBS = [
    (r"Work\GP\HPGe(20%)", "**/*.spe", "GP_HPGe20", "P1"),
    (r"Work\Handy\Handy(HPGe)", "**/*.spe", "Handy_HPGe", "P1"),
    (r"Work\Handy\Handy(NaI)", "**/*.spe", "Handy_NaI", "P1"),
    (r"Work\Handy\Handy(LaBr)", "**/*.spe", "Handy_LaBr", "P1"),
]
P2_GLOBS = [
    (r"Work\NM\HPGe(20%)", "**/*.spe", "NM_HPGe20", "P2"),
]
P3_GLOBS = [
    (r"Work\Simple\Alpha(Demo)", "**/*.spe", "Simple_Alpha", "P3"),
    (r"Work\Simple\HPGe(Demo)", "**/*.spe", "Simple_HPGe", "P3"),
    (r"Work\Simple\NaI(Demo)", "**/*.spe", "Simple_NaI", "P3"),
    (r"Work\Simple\SiLi(Demo)", "**/*.spe", "Simple_SiLi", "P3"),
    (r"Work\Simple\TeCd(Demo)", "**/*.spe", "Simple_TeCd", "P3"),
]
P2_SAMPLE_CAP = 10
P3_SAMPLE_CAP = 5

SCHEMA_VERSION = "0.2"
GENERATOR_PATH = "scripts/rag/build_spectra_index.py"

# F-300 W3: SCHEMA_VERSION bumped 0.1 → 0.2 (2026-06-05).
# Changelog (also surfaced in _meta.schema_changelog):
#   0.2 — Added per-record ``detector_folder`` field tying each spectrum to
#         its canonical ``detectors/<class>/`` folder. Enables per-class
#         SPECTRA_MANIFEST.json emission (F-300/W3 deliverable).
#   0.1 — Initial schema (v1.25.x); A2 design doc 2026-06-04.
SCHEMA_CHANGELOG: dict[str, str] = {
    "0.1": "Initial schema (v1.25.x).",
    "0.2": "Added detector_folder field per F-300 W3 (2026-06-05).",
}

# F-300 W3: detector_tag → detector_folder canonical mapping.
# Each new-class (W2 created folders) tag from ``_enumerate_targets`` ⇒
# the matching ``detectors/<class>/`` folder. Legacy Gamma-1S synonym maps
# to Gamma-1S (post-2026-06-05 taxonomy lock).
_DETECTOR_TAG_TO_FOLDER: dict[str, str] = {
    "Gamma-1S": "detectors/Gamma-1S/",
    "Gamma-1S": "detectors/Gamma-1S/",  # legacy synonym (canonicalized → C)
    "Handy_LaBr": "detectors/Handy_LaBr/",
    "Handy_HPGe": "detectors/Handy_HPGe/",
    "Handy_NaI": "detectors/Handy_NaI/",
    "GP_HPGe20": "detectors/GP_HPGe20/",
    "Simple_HPGe": "detectors/Simple_HPGe/",
    "Simple_NaI": "detectors/Simple_NaI/",
    "Simple_TeCd": "detectors/Simple_TeCd/",
    "Simple_SiLi": "detectors/Simple_SiLi/",
    "Simple_Alpha": "detectors/Simple_Alpha/",
}


def _detector_folder_for(detector_tag: str) -> str:
    """Return the canonical ``detectors/<class>/`` folder for a detector tag.

    Pure function. Empty / whitespace-only inputs return ``""`` silently
    (legitimate "no detector tag yet" case). Non-empty *unmapped* tags
    return ``""`` AND emit a :class:`UserWarning` so the caller / test
    suite has a visible signal that ``_DETECTOR_TAG_TO_FOLDER`` is out of
    sync with the data feed (V126-04, 2026-06-05). The mapping table is
    exhaustive over all 11 detector tags currently emitted by
    ``_enumerate_targets`` (10 new-classes plus legacy ``Gamma-1S``
    synonym).
    """
    if not detector_tag:
        return ""
    key = detector_tag.strip()
    if not key:
        return ""
    folder = _DETECTOR_TAG_TO_FOLDER.get(key, "")
    if not folder:
        import warnings
        warnings.warn(
            f"_detector_folder_for: unmapped detector tag {key!r}; "
            f"detector_folder will be empty string in SPECTRA_INDEX. "
            f"Add it to _DETECTOR_TAG_TO_FOLDER.",
            UserWarning,
            stacklevel=2,
        )
    return folder

# ─── Canonical geometry remaps (operator-typo collapse) ──────────────────────
# User lock 2026-06-05 (operator typo confirmation):
#   «Дента-120 (2 records) — это опечатка, должно быть Дента-120мл»
# Applied at .spe-header ingestion time so the canonical form propagates to
# both the per-record ``geometry`` field AND the ``indexes.by_geometry``
# secondary lookup. Pre-patch: ``Дента-120`` = 2 records, ``Дента-120мл`` =
# 23 records, separate keys; post-patch: ``Дента-120мл`` = 25, ``Дента-120``
# absent.
#
# F-300 W4 additions (2026-06-05): canonicalize 4 latin geometry strings
# emitted by the Handy_* / GP_HPGe20 / Simple_* .spe headers into their
# cyrillic equivalents, so by_geometry shows ONE key per logical geometry
# rather than two (latin-token and cyrillic-token siblings).
#
# Brief lock: «Do NOT merge Маринелли ≡ Маринелли 1л — volume distinction
# is meaningful». Only the latin↔cyrillic synonym pairs collapse here.
#
# Out of scope for this remap (see brief 2026-06-05): Дента-100 (3 records
# — PENDING source reconciliation per
# detectors/Gamma-1S/README.md §"Дента-100 — pending reconciliation").
_CANONICAL_GEOMETRY_REMAPS: dict[str, str] = {
    "Дента-120": "Дента-120мл",
    # F-300 W4 additions (2026-06-05) — Handy_*/GP_HPGe20 latin↔cyrillic:
    "Point24": "Точечная",       # Handy_LaBr (logical equivalence, ≈23 records)
    "Point-15cm": "Точечный",    # Handy_HPGe / Handy_NaI (≈70 records)
    "MARINELLI": "Маринелли",    # GP_HPGe20 case variant (≈11 records)
    "Marinelli": "Маринелли",    # GP_HPGe20 case variant (1 record)
}


def _canonicalize_geometry(raw: str) -> str:
    """Apply operator-typo remap to a raw .spe geometry string.

    Pure function; idempotent; keys not in the remap table pass through
    unchanged. Empty/whitespace-only inputs return ``""``.
    """
    if not raw:
        return ""
    key = raw.strip()
    return _CANONICAL_GEOMETRY_REMAPS.get(key, key)


# ─── Canonical detector remaps (taxonomy lock: cyrillic↔latin homoglyph) ─────
# User lock 2026-06-05:
#   «Гамма-1с (кириллица) = Gamma-1S (латиница, омоглиф)». В проекте одна
#   физическая NaI-63×63 станция, её canonical name = Gamma-1S, её папка =
#   detectors/Gamma-1S/.
#
# The legacy ``Gamma-1S`` token is emitted by the LSRM-header parser
# (cyrillic «Гамма-1С» → latin transliteration of cyrillic «С» = «S» by
# sound). v1.11.1 released ``data/aliases.json:detector`` containing
# ``Gamma-1S`` as a legacy synonym for backwards compat with that parser
# (and BUG-40 raises a defensive warning when this aliasing fires). The
# remap below applies the post-2026-06-05 canon at SPECTRA_INDEX-build
# time so:
#   - per-record ``detector_tag`` carries the canonical ``Gamma-1S``;
#   - ``indexes.by_detector`` keys on ``Gamma-1S`` (matches the canonical
#     folder name ``detectors/Gamma-1S/`` and the 24 fitted VT-*.json
#     ``station_observed_on: "Gamma-1S"`` produced by F-070 W4 S0).
#
# Out of scope for this remap: ``data/aliases.json`` (released invariant
# v1.11.1) and ``tests/step04_detector_type/test_bug40_cyrillic_latin_warning.py``
# (released BUG-40 defensive warning) — both retain ``Gamma-1S`` as the
# legacy ingestion-side synonym.
_CANONICAL_DETECTOR_REMAPS: dict[str, str] = {
    "Gamma-1S": "Gamma-1S",
}


def _canonicalize_detector(raw: str) -> str:
    """Apply taxonomy-lock remap to a raw detector tag.

    Pure function; idempotent; keys not in the remap table pass through
    unchanged. Empty/whitespace-only inputs return ``""``.
    """
    if not raw:
        return ""
    key = raw.strip()
    return _CANONICAL_DETECTOR_REMAPS.get(key, key)

TIER_DEFINITIONS = {
    "A": "passport_with_ref_date AND energy_cal AND counts>1000 AND dead<50% AND live>0",
    "B": "energy_cal AND counts>1000 AND live>0 (no passport, not 2012-10-08 sentinel)",
    "C": "counts>0 AND live>0, failed A/B (sentinel date or no calibration)",
    "D": "live=0 OR counts<100 OR parse error",
}


# ─── Module-level helpers ────────────────────────────────────────────────────

def _rel_placeholder(abs_path: Path, lsrm_root: Path) -> str:
    """F-115 — replace operator's LSRM root with ``<LSRM>`` placeholder.

    Falls back to absolute path stringification only if the file is not
    under ``lsrm_root`` (this should not happen for files yielded by the
    scope-glob enumerator, but is defensive).
    """
    try:
        rel = abs_path.relative_to(lsrm_root)
        return f"{LSRM_PLACEHOLDER}\\{rel}"
    except ValueError:
        return str(abs_path)


def _sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _energy_range(energy_cal, n_channels) -> tuple[Optional[float], Optional[float]]:
    if not energy_cal or n_channels <= 0:
        return (None, None)
    def _eval(ch):
        return sum(a * (ch ** i) for i, a in enumerate(energy_cal))
    return (float(_eval(0)), float(_eval(n_channels - 1)))


def _enumerate_targets(lsrm_root: Path) -> list[tuple[str, str, str]]:
    """Return list of ``(abs_path_str, detector_tag, priority)`` tuples.

    Applies sample caps to P2/P3 (deterministic by sorted file order).
    """
    out: list[tuple[str, str, str]] = []
    for (subdir, pattern, tag, prio) in (
        P0_GLOBS + P1_GLOBS + P2_GLOBS + P3_GLOBS
    ):
        base = lsrm_root / subdir
        if not base.exists():
            continue
        matches = sorted(base.glob(pattern))
        if prio == "P2":
            matches = matches[:P2_SAMPLE_CAP]
        elif prio == "P3":
            matches = matches[:P3_SAMPLE_CAP]
        for p in matches:
            if p.is_file():
                out.append((str(p), tag, prio))
    return out


# ─── Worker (top-level so multiprocessing.Pool can pickle it on Windows) ─────

def _parse_one(args: tuple[str, str, str, str]) -> dict:
    """Parse a single ``.spe`` and return the per-record dict.

    The dict shape here is the *intermediate* shape (matches the prototype
    eval output 1:1 except for schema-version-ish keys). The intermediate
    shape is later projected to the SPECTRA_INDEX.json schema by
    :func:`_to_index_record`.

    Args tuple: ``(abs_path_str, detector_tag, priority, lsrm_root_str)``.
    """
    abs_path_str, detector_tag, priority, lsrm_root_str = args
    abs_path = Path(abs_path_str)
    lsrm_root = Path(lsrm_root_str)
    # Canonicalize detector tag at ingestion (see ``_CANONICAL_DETECTOR_REMAPS``
    # above) so the canonical form propagates to the per-record
    # ``detector_tag`` field AND the ``indexes.by_detector`` secondary
    # lookup. Scope-glob tuples retain the legacy ``Gamma-1S`` literal
    # (matches the LSRM folder name on disk: ``Work\BG\Gamma-1S``).
    detector_tag = _canonicalize_detector(detector_tag)
    rec: dict = {
        "rel_path": _rel_placeholder(abs_path, lsrm_root),
        "detector_tag": detector_tag,
        "priority": priority,
    }
    try:
        raw = abs_path.read_bytes()
        rec["file_size_bytes"] = len(raw)
        rec["sha256"] = _sha256_bytes(raw)

        from gamma.io.lsrm_spe import read_lsrm_spe
        spec = read_lsrm_spe(str(abs_path))

        live_t = float(spec.live_time)
        real_t = float(spec.real_time)
        n_ch = int(spec.n_channels)

        rec["live_time_s"] = live_t
        rec["real_time_s"] = real_t
        rec["dead_time_frac"] = (
            float((real_t - live_t) / real_t) if real_t > 0 else None
        )
        rec["channels"] = n_ch
        try:
            rec["total_counts"] = int(spec.counts.sum())
        except Exception:
            rec["total_counts"] = None

        rec["calib_energy"] = (
            [float(c) for c in spec.energy_cal] if spec.energy_cal else None
        )
        rec["energy_cal_source"] = getattr(spec, "energy_cal_source", None)
        if spec.stored_fwhm_calibration:
            rec["calib_fwhm"] = [
                float(c) for c in spec.stored_fwhm_calibration.coefficients
            ]
            rec["calib_fwhm_model"] = spec.stored_fwhm_calibration.model
        else:
            rec["calib_fwhm"] = None
            rec["calib_fwhm_model"] = None

        e_min, e_max = _energy_range(spec.energy_cal, n_ch)
        rec["energy_min_keV"] = e_min
        rec["energy_max_keV"] = e_max

        rec["comment_raw"] = (spec.comments or "").strip()
        # Canonicalize operator-typo geometry strings at ingestion (see
        # ``_CANONICAL_GEOMETRY_REMAPS`` above). This ensures both the
        # per-record ``geometry`` field and the ``indexes.by_geometry``
        # secondary lookup show the canonical form.
        rec["geometry"] = _canonicalize_geometry(spec.geometry or "")
        rec["sample_id"] = spec.sample_id or ""
        rec["detector_id_header"] = spec.detector_id or ""

        passport_entries: list[dict] = []
        ex = spec.extras.get("lsrm_passport") if hasattr(spec, "extras") else None
        if ex:
            passport_entries = [
                {
                    "nuclide": e.get("nuclide"),
                    "value": e.get("value"),
                    "unit": e.get("unit"),
                    "uncertainty_pct": e.get("uncertainty_pct"),
                    "reference_date": e.get("reference_date"),
                    "is_specific_activity": e.get("is_specific_activity"),
                    "source": "f330_canonical",
                }
                for e in ex
            ]

        # Поверка 2016 fallback parser (regex on COMMENT lines like
        # "Nuc A=v unit dA=u% date") — duplicate of A's pre-BUG-49 helper.
        # Kept inline so we do not depend on the draft probe script.
        if not passport_entries and rec["comment_raw"]:
            for entry in _parse_poverka_2016_comment(rec["comment_raw"]):
                passport_entries.append(entry)

        rec["passport_entries"] = passport_entries
        rec["passport_count"] = len(passport_entries)

        if spec.start_datetime is not None:
            rec["acq_started_at"] = spec.start_datetime.isoformat()
        else:
            rec["acq_started_at"] = None

        # ─── Quality tier scoring (deterministic, per A2 design §1) ──
        flags: list[str] = []
        if live_t <= 0:
            flags.append("zero_live_time")
        if rec.get("total_counts") is not None and rec["total_counts"] < 100:
            flags.append("low_counts_lt100")
        if rec.get("dead_time_frac") is not None and rec["dead_time_frac"] > 0.5:
            flags.append("high_dead_time")
        if not rec["calib_energy"]:
            flags.append("no_energy_calibration")
        if (rec.get("acq_started_at") or "").startswith("2012-10-08"):
            flags.append("sentinel_demo_date_2012_10_08")

        has_passport = (
            rec["passport_count"] > 0
            and any(e.get("reference_date") for e in passport_entries)
        )
        has_calib = bool(rec["calib_energy"])
        counts_ok = (
            rec.get("total_counts") is not None and rec["total_counts"] > 1000
        )
        dt_ok = (
            rec.get("dead_time_frac") is None
            or rec["dead_time_frac"] < 0.5
        )

        if (
            live_t > 0
            and counts_ok
            and dt_ok
            and has_calib
            and has_passport
        ):
            tier = "A"
        elif live_t > 0 and counts_ok and has_calib:
            if "sentinel_demo_date_2012_10_08" in flags:
                tier = "C"
            else:
                tier = "B"
        elif (
            live_t > 0
            and rec.get("total_counts") is not None
            and rec["total_counts"] > 0
        ):
            tier = "C"
        else:
            tier = "D"

        rec["quality_tier"] = tier
        rec["quality_flags"] = flags
        rec["error"] = None
    except Exception as ex:
        rec["quality_tier"] = "D"
        rec["quality_flags"] = ["parse_error"]
        rec["error"] = f"{type(ex).__name__}: {ex}"
    return rec


# Поверка 2016 COMMENT format regex.
# Examples:
#   "Cs-137 A=10500 Бк dA=2.5% 30.05.2016"
#   "K-40   A=2540 Бк/кг dA=10% 17.05.2016"
_POVERKA_2016_RE = re.compile(
    r"(?P<nuc>[A-Za-zА-Яа-я]+-\d+[mМ]?)\s+"
    r"A\s*=\s*(?P<val>[\d.,eE+\-]+)\s*"
    r"(?P<unit>Бк/?кг|Бк·кг⁻¹|Бк/л|Bq/?kg|Бк)?\s*"
    r"dA\s*=\s*(?P<unc>[\d.,]+)\s*%\s*"
    r"(?P<date>\d{2}\.\d{2}\.(?:\d{4}|\d{2}))",
    re.UNICODE,
)


def _parse_poverka_2016_comment(comment: str) -> list[dict]:
    """Best-effort fallback regex for Поверка 2016 layout.

    Mirror of the helper in the prototype (kept inline to avoid a draft-
    script dependency). BUG-49 wave will fold this into
    ``gamma.io.lsrm_passport``; until then this provides supplementary
    A-tier extraction.
    """
    out: list[dict] = []
    for m in _POVERKA_2016_RE.finditer(comment):
        date_raw = m.group("date")
        # Normalise DD.MM.YY → DD.MM.20YY (LSRM tree starts in late 1990s,
        # but the Поверка 2016 layout post-dates 2000 — safe heuristic).
        parts = date_raw.split(".")
        if len(parts[-1]) == 2:
            year = int(parts[-1])
            yyyy = 1900 + year if year >= 90 else 2000 + year
            date_iso = f"{yyyy:04d}-{int(parts[1]):02d}-{int(parts[0]):02d}"
        else:
            date_iso = f"{int(parts[2]):04d}-{int(parts[1]):02d}-{int(parts[0]):02d}"
        try:
            value = float(m.group("val").replace(",", "."))
        except ValueError:
            continue
        try:
            unc = float(m.group("unc").replace(",", "."))
        except ValueError:
            unc = None
        unit = m.group("unit") or "Бк"
        is_specific = "/" in unit or "·" in unit
        out.append({
            "nuclide": m.group("nuc"),
            "value": value,
            "unit": unit,
            "uncertainty_pct": unc,
            "reference_date": date_iso,
            "is_specific_activity": is_specific,
            "source": "poverka_2016_re",
        })
    return out


# ─── Index assembly ──────────────────────────────────────────────────────────

def _to_index_record(intermediate: dict, spectrum_id: str) -> dict:
    """Project the prototype's intermediate record to the SPECTRA_INDEX schema.

    Schema fields (per design §2):

    - ``spectrum_id`` — stable ``RAG-SPEC-NNNN``
    - ``sha256`` — dedup primary key
    - ``rel_path`` — F-115 ``<LSRM>\\...`` placeholder path
    - ``quality_tier`` — A/B/C/D
    - ``detector_tag``, ``detector_id_header``, ``geometry``
    - ``live_time_s``, ``real_time_s``, ``dead_time_frac``, ``channels``
    - ``total_counts``, ``acq_started_at``
    - ``energy_cal``, ``fwhm_cal``, ``fwhm_model``, ``energy_range_keV``
    - ``passport`` (list of ``{nuclide, value_Bq, uncertainty_pct,
      reference_date, is_specific_activity}``)
    - ``tags`` (computed: tier, geometry hints, drift-cohort tag)
    - ``use_cases`` (computed from tier × passport × calib)
    - ``linked_rag_methodology`` (default empty list; manual augmentation
      via RAG_INDEX.json side, not here)
    - ``quality_flags`` (preserved for diagnostics)
    """
    passport: list[dict] = []
    for e in intermediate.get("passport_entries") or []:
        try:
            value_bq = float(e["value"]) if e.get("value") is not None else None
        except (TypeError, ValueError):
            value_bq = None
        try:
            unc = float(e["uncertainty_pct"]) if e.get("uncertainty_pct") is not None else None
        except (TypeError, ValueError):
            unc = None
        passport.append({
            "nuclide": e.get("nuclide"),
            "value_Bq": value_bq,
            "unit": e.get("unit"),
            "uncertainty_pct": unc,
            "reference_date": e.get("reference_date"),
            "is_specific_activity": bool(e.get("is_specific_activity")),
            "parser_source": e.get("source"),
        })

    tier = intermediate.get("quality_tier", "D")
    geometry = intermediate.get("geometry") or ""
    detector_tag = intermediate.get("detector_tag", "?")

    # Tag synthesis (small, intentional set per design §2 illustration).
    tags: list[str] = [
        f"tier-{tier}",
        f"priority-{intermediate.get('priority', '?')}",
    ]
    if tier == "A":
        tags.append("production")
    geom_lower = geometry.lower()
    if "маринел" in geom_lower or "marinelli" in geom_lower:
        tags.append("marinelli=true")
    else:
        tags.append("marinelli=false")
    if "точечн" in geom_lower or "point" in geom_lower:
        tags.append("point-source")
    if intermediate.get("comment_raw") and any(
        s.startswith("poverka_2016") for s in
        (e.get("parser_source") for e in passport)
    ):
        tags.append("poverka-2016")

    # Drift-cohort tag: only if A-tier with passport AND single-nuclide.
    nuclides_in_passport = sorted({
        p["nuclide"] for p in passport if p.get("nuclide")
    })
    if tier == "A" and len(nuclides_in_passport) == 1:
        tags.append(f"drift-cohort:{detector_tag}/{nuclides_in_passport[0]}")

    # Use cases derived from tier + capabilities.
    use_cases: list[str] = []
    if tier == "A":
        use_cases.append("F-QC-01-validation")
        if intermediate.get("calib_fwhm"):
            use_cases.append("FWHM-anchor")
        if passport and intermediate.get("calib_energy"):
            use_cases.append("efficiency-anchor")
    elif tier == "B":
        use_cases.append("background-survey")

    return {
        "spectrum_id": spectrum_id,
        "sha256": intermediate["sha256"],
        "rel_path": intermediate["rel_path"],
        "quality_tier": tier,
        "detector_tag": detector_tag,
        "detector_folder": _detector_folder_for(detector_tag),
        "detector_id_header": intermediate.get("detector_id_header") or "",
        "geometry": geometry,
        "sample_id": intermediate.get("sample_id") or "",
        "priority": intermediate.get("priority"),
        "live_time_s": intermediate.get("live_time_s"),
        "real_time_s": intermediate.get("real_time_s"),
        "dead_time_frac": intermediate.get("dead_time_frac"),
        "channels": intermediate.get("channels"),
        "total_counts": intermediate.get("total_counts"),
        "energy_cal": intermediate.get("calib_energy"),
        "energy_cal_source": intermediate.get("energy_cal_source"),
        "fwhm_cal": intermediate.get("calib_fwhm"),
        "fwhm_model": intermediate.get("calib_fwhm_model"),
        "energy_range_keV": [
            intermediate.get("energy_min_keV"),
            intermediate.get("energy_max_keV"),
        ],
        "acq_started_at": intermediate.get("acq_started_at"),
        "passport": passport,
        "tags": tags,
        "use_cases": use_cases,
        "quality_flags": intermediate.get("quality_flags") or [],
        "parse_error": intermediate.get("error"),
        "linked_rag_methodology": [],
    }


def _build_indexes(records: list[dict]) -> dict:
    """Build the secondary lookup tables (`by_nuclide`, `by_detector`, ...).

    All output lists are sorted ASC by ``spectrum_id`` for deterministic
    re-runs.
    """
    by_nuclide: dict[str, list[str]] = {}
    by_detector: dict[str, list[str]] = {}
    by_geometry: dict[str, list[str]] = {}
    by_quality_tier: dict[str, list[str]] = {"A": [], "B": [], "C": [], "D": []}
    by_passport_year: dict[str, list[str]] = {}
    drift_cohort_members: dict[tuple[str, str], dict] = {}

    for r in records:
        sid = r["spectrum_id"]
        det = r.get("detector_tag") or "?"
        geom = r.get("geometry") or ""
        tier = r.get("quality_tier", "D")

        by_detector.setdefault(det, []).append(sid)
        if geom:
            by_geometry.setdefault(geom, []).append(sid)
        by_quality_tier.setdefault(tier, []).append(sid)

        for p in r.get("passport") or []:
            nuc = p.get("nuclide")
            if nuc:
                by_nuclide.setdefault(nuc, []).append(sid)
            ref_date = p.get("reference_date") or ""
            if len(ref_date) >= 4:
                yr = ref_date[:4]
                by_passport_year.setdefault(yr, []).append(sid)

        # Drift cohort accumulation: A-tier with single passport nuclide.
        if tier == "A":
            unique_nucs = sorted({
                p.get("nuclide") for p in (r.get("passport") or [])
                if p.get("nuclide")
            })
            if len(unique_nucs) == 1:
                key = (det, unique_nucs[0])
                cohort = drift_cohort_members.setdefault(key, {
                    "cohort_id": f"DRIFT-{unique_nucs[0]}-{det}",
                    "detector_tag": det,
                    "nuclide": unique_nucs[0],
                    "members": [],
                    "years": set(),
                })
                cohort["members"].append(sid)
                ts = r.get("acq_started_at") or ""
                if len(ts) >= 4:
                    cohort["years"].add(ts[:4])

    # Deduplicate + sort all lists.
    def _sorted_unique(lst: list[str]) -> list[str]:
        return sorted(set(lst))

    by_nuclide = {k: _sorted_unique(v) for k, v in sorted(by_nuclide.items())}
    by_detector = {k: _sorted_unique(v) for k, v in sorted(by_detector.items())}
    by_geometry = {k: _sorted_unique(v) for k, v in sorted(by_geometry.items())}
    by_quality_tier = {
        k: _sorted_unique(by_quality_tier.get(k, []))
        for k in ("A", "B", "C", "D")
    }
    by_passport_year = {
        k: _sorted_unique(v) for k, v in sorted(by_passport_year.items())
    }

    # Promote drift cohorts: only keep cohorts with ≥3 distinct year stamps
    # (matches A2 design §1 "drift-study quadruples" criterion).
    drift_cohorts: list[dict] = []
    for (det, nuc), info in sorted(drift_cohort_members.items()):
        years = sorted(info["years"])
        if len(years) >= 3:
            drift_cohorts.append({
                "cohort_id": info["cohort_id"],
                "detector_tag": det,
                "nuclide": nuc,
                "members": _sorted_unique(info["members"]),
                "year_span": years,
            })
    # Stable sort: detector ASC, nuclide ASC (already sorted above).

    return {
        "by_nuclide": by_nuclide,
        "by_detector": by_detector,
        "by_geometry": by_geometry,
        "by_quality_tier": by_quality_tier,
        "by_passport_year": by_passport_year,
        "drift_cohorts": drift_cohorts,
    }


def _build_duplicates(records: list[dict]) -> list[dict]:
    """Compute SHA-256 duplicate groups across the index.

    Each entry: ``{sha256, spectrum_ids: [...sorted...]}``.
    """
    by_sha: dict[str, list[str]] = {}
    for r in records:
        by_sha.setdefault(r["sha256"], []).append(r["spectrum_id"])
    dups = []
    for sha, ids in sorted(by_sha.items()):
        if len(ids) > 1:
            dups.append({"sha256": sha, "spectrum_ids": sorted(set(ids))})
    return dups


# ─── Stable-ID assignment ────────────────────────────────────────────────────

def _load_existing_id_map(output_path: Path) -> dict[str, str]:
    """Return ``{sha256: spectrum_id}`` from an existing payload.

    If the file does not exist or cannot be parsed, returns an empty map
    (genesis case — all IDs newly assigned).
    """
    if not output_path.exists():
        return {}
    try:
        data = json.loads(output_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    out: dict[str, str] = {}
    for rec in data.get("spectra", []):
        sha = rec.get("sha256")
        sid = rec.get("spectrum_id")
        if sha and sid:
            out[sha] = sid
    return out


def _assign_stable_ids(
    intermediates: list[dict],
    existing_map: dict[str, str],
) -> list[tuple[dict, str]]:
    """Map each intermediate record to a stable ``RAG-SPEC-NNNN`` id.

    Strategy:
    1. Collapse by sha256 (first occurrence wins for the chosen
       intermediate; duplicates surface in the dedup table).
    2. Reuse IDs from ``existing_map`` where sha256 matches.
    3. Assign new IDs in sha256-ASC order, starting from the next free
       slot after the highest existing ID.
    """
    # Deduplicate by sha256, keeping the intermediate record with the
    # ASC-smallest rel_path (deterministic choice; multiprocessing.Pool
    # may yield records in non-deterministic order so we must impose a
    # fixed tie-breaker on the canonical winner).
    by_sha: dict[str, dict] = {}
    for r in sorted(
        intermediates,
        key=lambda x: (x.get("sha256") or "", x.get("rel_path") or ""),
    ):
        sha = r.get("sha256")
        if not sha:
            continue
        if sha not in by_sha:
            by_sha[sha] = r

    # Determine the high-water mark from existing IDs.
    max_id = 0
    for sid in existing_map.values():
        m = re.match(r"^RAG-SPEC-(\d+)$", sid)
        if m:
            n = int(m.group(1))
            if n > max_id:
                max_id = n

    out: list[tuple[dict, str]] = []
    next_free = max_id + 1
    for sha in sorted(by_sha):
        rec = by_sha[sha]
        if sha in existing_map:
            sid = existing_map[sha]
        else:
            sid = f"RAG-SPEC-{next_free:04d}"
            next_free += 1
        out.append((rec, sid))
    return out


# ─── Top-level pipeline ──────────────────────────────────────────────────────

def build_payload(
    lsrm_root: Path,
    existing_output: Optional[Path] = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """End-to-end build of the SPECTRA_INDEX.json payload.

    Returns a fully-assembled dict ready to write. Does not write to disk.

    Args:
        lsrm_root: operator's LSRM directory (must exist).
        existing_output: optional path to a previous SPECTRA_INDEX.json
            for ID-preservation lookup. ``None`` ⇒ genesis run.
        verbose: emit progress lines to stderr.
    """
    targets = _enumerate_targets(lsrm_root)
    if verbose:
        print(f"[build_spectra_index] Enumerated {len(targets)} target .spe files.",
              file=sys.stderr)
    if not targets:
        raise RuntimeError(
            f"No .spe targets enumerated under {lsrm_root}; check scope globs."
        )

    work = [(p, det, prio, str(lsrm_root)) for (p, det, prio) in targets]

    n_workers = min(31, max(1, len(work) // 8))
    if verbose:
        print(f"[build_spectra_index] Spawning Pool({n_workers}) ...",
              file=sys.stderr)

    intermediates: list[dict] = []
    with Pool(processes=n_workers) as pool:
        for i, rec in enumerate(
            pool.imap_unordered(_parse_one, work, chunksize=8), 1
        ):
            intermediates.append(rec)
            if verbose and i % 100 == 0:
                print(f"[build_spectra_index]   {i}/{len(work)} done",
                      file=sys.stderr)

    existing_map = _load_existing_id_map(existing_output) if existing_output else {}
    if verbose:
        print(f"[build_spectra_index] Loaded {len(existing_map)} existing IDs "
              f"from previous output.", file=sys.stderr)

    assigned = _assign_stable_ids(intermediates, existing_map)
    if verbose:
        new_ids = sum(1 for inter, sid in assigned if inter["sha256"] not in existing_map)
        print(f"[build_spectra_index] Assigned IDs: {len(assigned)} total "
              f"({new_ids} new, {len(assigned) - new_ids} preserved).",
              file=sys.stderr)

    # Project intermediates → index records.
    records = [_to_index_record(inter, sid) for inter, sid in assigned]
    # Records are already in sha256-ASC order via _assign_stable_ids().
    # Re-sort by spectrum_id for output stability (equivalent ordering).
    records.sort(key=lambda r: r["spectrum_id"])

    indexes = _build_indexes(records)
    duplicates = _build_duplicates_with_dup_paths(
        records=records, all_intermediates=intermediates,
    )

    # Per-record extraction log (provenance, mirrors the Ollama helper
    # convention even though no Ollama is used here).
    extraction_log = {
        "total_intermediate_records": len(intermediates),
        "total_unique_sha256": len(records),
        "n_parse_errors": sum(1 for r in intermediates if r.get("error")),
        "n_tier_A": sum(1 for r in records if r["quality_tier"] == "A"),
        "n_tier_B": sum(1 for r in records if r["quality_tier"] == "B"),
        "n_tier_C": sum(1 for r in records if r["quality_tier"] == "C"),
        "n_tier_D": sum(1 for r in records if r["quality_tier"] == "D"),
        "n_drift_cohorts": len(indexes["drift_cohorts"]),
        "n_dedup_groups": len(duplicates),
        "scope_globs": {
            "P0": [g[0] for g in P0_GLOBS],
            "P1": [g[0] for g in P1_GLOBS],
            "P2_cap": P2_SAMPLE_CAP,
            "P2": [g[0] for g in P2_GLOBS],
            "P3_cap": P3_SAMPLE_CAP,
            "P3": [g[0] for g in P3_GLOBS],
            "skipped_oos": {
                "ADA_alpha": "Out of F-150 gamma scope",
                "BG_Beta-1S": "Out of F-150 gamma scope",
                "NM_Pu": "Separate NM methodology contract",
                "NM_U": "Separate NM methodology contract",
            },
        },
    }

    payload = {
        "schema_version": SCHEMA_VERSION,
        "generator": GENERATOR_PATH,
        "lsrm_root_placeholder": LSRM_PLACEHOLDER,
        "tier_definitions": TIER_DEFINITIONS,
        "spectra": records,
        "indexes": indexes,
        "duplicates": duplicates,
        "_extraction_log": extraction_log,
        "_meta": {
            "generated_at_utc": datetime.utcnow().isoformat() + "Z",
            "parser": "gamma.io.lsrm_spe.read_lsrm_spe",
            "passport_parsers": ["f330_canonical", "poverka_2016_re_fallback"],
            "id_scheme": "RAG-SPEC-NNNN, deterministic sha256-ASC assignment; preserved across re-runs.",
            "schema_version": SCHEMA_VERSION,
            "schema_changelog": dict(SCHEMA_CHANGELOG),
            "deterministic_keys_note": (
                "Re-runs over the same LSRM tree produce byte-identical "
                "payloads EXCLUDING _meta.generated_at_utc. Tests assert "
                "this invariant."
            ),
        },
    }
    return payload


def _build_duplicates_with_dup_paths(
    records: list[dict],
    all_intermediates: list[dict],
) -> list[dict]:
    """Dedup table that also surfaces ALL distinct rel_paths per sha.

    The schema says ``duplicates[*].spectrum_ids`` (per design §5). We
    augment with ``paths`` for human inspection — multiple file system
    copies of the same content under one spectrum_id.
    """
    # Map sha → list of intermediate-record rel_paths
    paths_by_sha: dict[str, list[str]] = {}
    for r in all_intermediates:
        sha = r.get("sha256")
        if sha:
            paths_by_sha.setdefault(sha, []).append(r["rel_path"])
    id_by_sha = {r["sha256"]: r["spectrum_id"] for r in records}

    out = []
    for sha in sorted(paths_by_sha):
        paths = sorted(set(paths_by_sha[sha]))
        if len(paths) > 1:
            sid = id_by_sha.get(sha)
            if not sid:
                continue
            out.append({
                "sha256": sha,
                "spectrum_ids": [sid],  # one canonical id per content
                "paths": paths,
            })
    return out


def write_payload(payload: dict[str, Any], output_path: Path) -> None:
    """Serialise payload with sorted keys + UTF-8 (deterministic on-disk)."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
    )


# ─── CLI ─────────────────────────────────────────────────────────────────────

def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build/refresh audit/_rag/SPECTRA_INDEX.json from "
                    "an operator's LSRM .spe tree."
    )
    p.add_argument(
        "--lsrm-root", required=True,
        help="Absolute path to the LSRM root (e.g. C:\\LSRM). NOT stored "
             "in the output — replaced by <LSRM> placeholder per F-115.",
    )
    p.add_argument(
        "--output",
        default=str(PROJ_ROOT / "audit" / "_rag" / "SPECTRA_INDEX.json"),
        help="Output JSON path. Default: <repo>/audit/_rag/SPECTRA_INDEX.json.",
    )
    p.add_argument(
        "--check", action="store_true",
        help="Dry run — build payload in memory and report tier counts to "
             "stderr without writing the output file.",
    )
    p.add_argument(
        "--quiet", action="store_true",
        help="Suppress per-batch progress lines on stderr.",
    )
    return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    t0 = time.time()
    args = _parse_args(argv)
    lsrm_root = Path(args.lsrm_root).resolve()
    output_path = Path(args.output).resolve()

    if not lsrm_root.exists():
        print(f"[build_spectra_index] LSRM root does not exist: {lsrm_root}",
              file=sys.stderr)
        return 1

    payload = build_payload(
        lsrm_root=lsrm_root,
        existing_output=output_path if output_path.exists() else None,
        verbose=not args.quiet,
    )

    log = payload["_extraction_log"]
    print(f"[build_spectra_index] tiers: "
          f"A={log['n_tier_A']} B={log['n_tier_B']} "
          f"C={log['n_tier_C']} D={log['n_tier_D']} "
          f"(unique sha256 = {log['total_unique_sha256']}, "
          f"intermediates = {log['total_intermediate_records']})",
          file=sys.stderr)

    if args.check:
        print(f"[build_spectra_index] --check requested; NOT writing output.",
              file=sys.stderr)
    else:
        write_payload(payload, output_path)
        size_kb = output_path.stat().st_size / 1024
        print(f"[build_spectra_index] Wrote {output_path} ({size_kb:.1f} KB)",
              file=sys.stderr)

    dt = time.time() - t0
    print(f"[build_spectra_index] DONE in {dt:.1f}s.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    freeze_support()
    try:
        sys.exit(main())
    except Exception:
        traceback.print_exc()
        sys.exit(1)
