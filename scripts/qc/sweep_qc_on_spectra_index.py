"""
F-QC-01 sweep over SPECTRA_INDEX A-tier (or other tier) fixtures.

Reads `audit/_rag/SPECTRA_INDEX.json`, resolves fixture paths against an
LSRM tree root, reads each `.spe` via `gamma.io.lsrm_spe.read_lsrm_spe`,
applies a lightweight proxy object into `build_spectrum_qc`, and emits a
per-fixture pass/fail verdict + aggregate counts.

Design notes
------------
* The sweep does NOT run the full staged pipeline (no peak-search, no
  identification).  Instead it builds a minimal `_QcProxy` object that
  populates only the attributes `build_spectrum_qc` actually reads:
    - fwhm_at_661        (from build_fwhm_model + fwhm_keV_at_energy)
    - fwhm_model_source  (from build_fwhm_model)
    - efficiency_curve   (attempted autoload via efficiency_autoload)
    - efficiency_source  (same)
    - spec.live_time     (from read_lsrm_spe)
    - seven_line_check   = None  (requires full calibration chain → skip)
    - background_subtraction = None  (no paired background → skip)
    - bg_quality_check   = None  (no peaks table → skip)
  This gives reliable verdicts for criteria 2 (FWHM) and 3 (efficiency)
  and conservative-pass for 1 (energy drift), 4 (bg drift), 5 (per-peak).

* Idempotent: if output file already exists, re-running overwrites it with
  a fresh sweep.

* 31-worker multiprocessing (Phase 1 MAXIMUM mandate, spectravibe-dev SKILL.md §«Compute policy»).

RAG-ID: [F-QC-01], [RAG-041]
Cite: spectrum_qc_methodology_v2_2026-06-03.md; KNOWN_AND_FIXED_ISSUES.md:1292
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import os
import sys
import traceback
from datetime import date
from multiprocessing import Pool
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# ─── Project path bootstrap ─────────────────────────────────────────────────
_HERE = Path(__file__).resolve()
_PROJECT_ROOT = _HERE.parent.parent.parent   # scripts/qc → scripts → project root
_SCRIPTS_DIR = _PROJECT_ROOT / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

# ─── Project imports ─────────────────────────────────────────────────────────
from gamma.io.lsrm_spe import read_lsrm_spe
from gamma.identification.staged_pipeline import build_fwhm_model, fwhm_keV_at_energy
from gamma.reporting.spectrum_qc_aggregator import build_spectrum_qc

import math as _math


def _fwhm_at_661_from_stored(spec) -> tuple:
    """
    Compute FWHM at 661.657 keV from the LSRM stored calibration if available.

    The LSRM header stores FWHM as a polynomial in sqrt(E) per BUG-22
    (confirmed in lsrm_spe.py comments):
        FWHM_keV(E) = sum(c_k * sqrt(E)^k)

    Falls back to build_fwhm_model() (internal FWHM²(E) model) when the
    stored calibration is absent.

    Returns (fwhm_keV: float, source: str).
    """
    E = 661.657
    sfwc = getattr(spec, "stored_fwhm_calibration", None)
    if sfwc is not None and sfwc.coefficients:
        coefs = sfwc.coefficients
        z = _math.sqrt(E)
        fwhm = sum(c * z ** i for i, c in enumerate(coefs))
        if fwhm > 0.1:
            return fwhm, "lsrm_stored_fwhm_polynomial_sqrt_E"
    # Fallback: pipeline FWHM² model from peaks table or NaI default
    model, src = build_fwhm_model(spec)
    return fwhm_keV_at_energy(model, E), src

# Efficiency autoload (may be unavailable for non-Gamma-1S detectors)
try:
    from gamma.calibration.efficiency_autoload import find_efr_file
    from gamma.calibration.efficiency import fit_efficiency_from_efr_file
    _EFF_AUTOLOAD_AVAILABLE = True
except ImportError:
    _EFF_AUTOLOAD_AVAILABLE = False

_SPECTRA_INDEX_PATH = _PROJECT_ROOT / "audit" / "_rag" / "SPECTRA_INDEX.json"
_WORKERS = 31
_LSRM_PLACEHOLDER = "<LSRM>"

# Per-detector-class FWHM reference at 661 keV (keV).
# Used by criterion 2 so that NaI and HPGe are evaluated against their own
# expected resolution — a single reference=8.5 keV is wrong for NaI.
#
# IMPORTANT: The default in spectrum_qc_aggregator.py (_FWHM_REFERENCE_KEV=8.5)
# is physically incorrect for NaI 63x63.  NaI FWHM at 661 keV is ~47 keV
# (7.1% resolution), NOT 8.5 keV.  8.5 keV corresponds to ~1.3% resolution
# which is neither NaI nor HPGe.  This sweep uses correct per-detector
# references; the aggregator default is flagged as a CROSS_TERRITORY finding.
#
# Reference values (FWHM at 661.657 keV):
#   NaI 63x63 (Gamma-1S/1S): ~47 keV  (default_NaI_63x63 model: sqrt(2.95*661)=44.2
#                                        + c·E² correction → 46.95 keV)
#                               Use 47 keV ±15% → [40, 54] keV window
#   HPGe (coaxial/BEGe)     :  2.0 keV (~0.3% at 661 keV)
#   LaBr3 (Handy_LaBr)      : ~13 keV  (~2% at 661 keV)
_FWHM_REFERENCE_BY_DETECTOR: Dict[str, float] = {
    "Gamma-1S":   47.0,   # NaI 63x63: 7.1% FWHM at 661 keV
    "Gamma-1S":   47.0,   # NaI 63x63: same detector class
    "Handy_HPGe":  1.5,   # HPGe coaxial (Handy): 1.46-2.1 keV FWHM at 661 keV measured
    "Handy_LaBr": 20.0,   # LaBr3 (Handy): ~3% FWHM at 661 keV (19.59 keV measured)
    "GP_HPGe20":   1.3,   # GP HPGe 20% coaxial: ~0.2% FWHM at 661 keV (1.28-1.30 keV measured)
}
_FWHM_REFERENCE_DEFAULT: float = 47.0   # NaI fallback (most common in LSRM)


# ════════════════════════════════════════════════════════════════════════════
# Lightweight QC proxy (avoids full staged pipeline)
# ════════════════════════════════════════════════════════════════════════════

class _SpecProxy:
    """Minimal spectrum proxy — only live_time needed by criterion 4."""
    def __init__(self, live_time: float) -> None:
        self.live_time = live_time


class _QcProxy:
    """
    Minimal StagedAnalysisResult substitute for build_spectrum_qc.

    Only the fields that build_spectrum_qc reads are populated:
      fwhm_at_661, fwhm_model_source — from build_fwhm_model()
      efficiency_curve, efficiency_source — attempted autoload
      spec — _SpecProxy(live_time=…)
      seven_line_check = None   → criterion 1 conservative pass
      background_subtraction = None → criterion 4 conservative pass
      bg_quality_check = None   → criterion 5 no-peaks pass
    """
    def __init__(
        self,
        fwhm_at_661: float,
        fwhm_model_source: str,
        live_time: float,
        efficiency_curve: Optional[object] = None,
        efficiency_source: str = "",
    ) -> None:
        self.fwhm_at_661 = fwhm_at_661
        self.fwhm_model_source = fwhm_model_source
        self.spec = _SpecProxy(live_time)
        self.efficiency_curve = efficiency_curve
        self.efficiency_source = efficiency_source
        self.seven_line_check = None          # criterion 1: conservative pass
        self.background_subtraction = None    # criterion 4: conservative pass
        self.bg_quality_check = None          # criterion 5: no peaks tested


# ════════════════════════════════════════════════════════════════════════════
# Per-fixture worker (called by Pool.map — must be module-level)
# ════════════════════════════════════════════════════════════════════════════

def _process_fixture(args: Tuple) -> Dict[str, Any]:
    """
    Process one spectrum entry.  Returns a result dict.
    Called in worker pool — exceptions are caught and recorded.
    """
    entry, lsrm_root = args
    spec_id = entry.get("spectrum_id", "?")
    sha256 = entry.get("sha256", "")
    rel_path = entry.get("rel_path", "")
    detector_tag = entry.get("detector_tag", "")
    geometry = entry.get("geometry", "")
    nuclides = [p.get("nuclide", "") for p in entry.get("passport", [])]

    # --- Resolve path ---
    abs_path = rel_path.replace(_LSRM_PLACEHOLDER, str(lsrm_root))

    result_base = {
        "spectrum_id": spec_id,
        "sha256": sha256,
        "rel_path": rel_path,        # F-115: original placeholder path, not absolute
        "detector_tag": detector_tag,
        "geometry": geometry,
        "nuclides": nuclides,
    }

    # --- Read .spe ---
    try:
        spec = read_lsrm_spe(abs_path)
    except Exception as exc:
        return {
            **result_base,
            "read_error": str(exc),
            "spectrum_qc": None,
            "overall_passed": None,
            "criteria_verdicts": {},
        }

    # --- Build FWHM at 661 keV ---
    # Uses stored LSRM calibration (poly_sqrt_E) when available (BUG-22).
    # Falls back to pipeline FWHM² model from peaks table or NaI default.
    try:
        fwhm_at_661, src = _fwhm_at_661_from_stored(spec)
    except Exception as exc:
        return {
            **result_base,
            "read_error": f"fwhm_model: {exc}",
            "spectrum_qc": None,
            "overall_passed": None,
            "criteria_verdicts": {},
        }

    # --- Try efficiency autoload ---
    eff_curve: Optional[object] = None
    eff_source: str = ""
    if _EFF_AUTOLOAD_AVAILABLE:
        # DEEP-01 (Project #5 wave 2 P1-1): distinguish "no .efr" from
        # "found .efr, fit failed".  The former is silent (criterion 3
        # fails as before with loaded=False).  The latter MUST emit a
        # warning so an operator running the sweep sees the corrupted
        # calibration file by basename — silent pass here previously
        # hid broken .efr files from the sweep report entirely.
        efr_path = None
        try:
            efr_path = find_efr_file(geometry, detector_tag)
        except Exception as exc:
            logger.warning(
                "sweep_qc: find_efr_file failed for geometry=%r detector=%r "
                "(%s: %s) — treating as 'no .efr available'.",
                geometry, detector_tag, type(exc).__name__, exc,
            )
        if efr_path:
            try:
                eff_curve = fit_efficiency_from_efr_file(efr_path)
                eff_source = os.path.basename(efr_path)   # F-115: leaf name only
            except Exception as exc:
                logger.warning(
                    "sweep_qc: efficiency fit failed for %r (%s: %s) — "
                    "criterion 3 will report loaded=False; investigate the .efr.",
                    os.path.basename(efr_path),
                    type(exc).__name__, exc,
                )
                eff_curve = None
                eff_source = "fit_failed"  # operator-visible marker
                # criterion 3 will still fail (loaded=False) but the
                # source string surfaces the cause distinct from "no file".

    # --- Build QC proxy ---
    proxy = _QcProxy(
        fwhm_at_661=fwhm_at_661,
        fwhm_model_source=src,
        live_time=float(spec.live_time),
        efficiency_curve=eff_curve,
        efficiency_source=eff_source,
    )

    # --- Run F-QC-01 aggregator ---
    # Use per-detector FWHM reference so criterion 2 is meaningful for both
    # NaI (8.5 keV) and HPGe (~2 keV). Single reference=8.5 only correct for NaI.
    fwhm_ref = _FWHM_REFERENCE_BY_DETECTOR.get(detector_tag, _FWHM_REFERENCE_DEFAULT)
    try:
        qc = build_spectrum_qc(proxy, fwhm_reference_keV=fwhm_ref)
    except Exception as exc:
        return {
            **result_base,
            "read_error": f"build_spectrum_qc: {exc}",
            "spectrum_qc": None,
            "overall_passed": None,
            "criteria_verdicts": {},
        }

    # --- Extract per-criterion verdicts (compact) ---
    criteria_verdicts = {
        "c1_energy_drift": {
            "passed": qc["energy_drift"]["passed"],
            "available": qc["energy_drift"]["available"],
            "max_residual_keV": qc["energy_drift"].get("max_residual_keV"),
            "note": qc["energy_drift"].get("note", ""),
        },
        "c2_fwhm_stability": {
            "passed": qc["fwhm_stability"]["passed"],
            "available": qc["fwhm_stability"]["available"],
            "fwhm_at_661_keV": qc["fwhm_stability"].get("fwhm_at_661_keV"),
            "rel_deviation": qc["fwhm_stability"].get("rel_deviation"),
            "reference_keV_used": fwhm_ref,
            "fwhm_model_source": src,
            "note": qc["fwhm_stability"].get("note", ""),
        },
        "c3_efficiency_qa": {
            "passed": qc["efficiency_qa"]["passed"],
            "efficiency_loaded": qc["efficiency_qa"].get("efficiency_loaded"),
            "note": qc["efficiency_qa"].get("note", ""),
        },
        "c4_bg_drift": {
            "passed": qc["bg_drift"]["passed"],
            "available": qc["bg_drift"]["available"],
            "note": qc["bg_drift"].get("note", ""),
        },
        "c5_peak_z_roi": {
            "passed": (qc["n_failed"] == 0),
            "n_peaks_tested": qc["n_peaks_tested"],
            "n_failed": qc["n_failed"],
        },
        "c6_sensitivity": {
            "passed": True,   # placeholder per v1.21.0
            "available": False,
        },
    }

    return {
        **result_base,
        "read_error": None,
        "spectrum_qc": {
            "overall_passed": qc["overall_passed"],
            "n_peaks_tested": qc["n_peaks_tested"],
            "n_passed": qc["n_passed"],
            "n_failed": qc["n_failed"],
        },
        "overall_passed": qc["overall_passed"],
        "criteria_verdicts": criteria_verdicts,
    }


# ════════════════════════════════════════════════════════════════════════════
# Aggregation helpers
# ════════════════════════════════════════════════════════════════════════════

def _aggregate(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Compute aggregate stats from per-fixture results."""
    n_total = len(results)
    n_error = sum(1 for r in results if r.get("read_error"))
    n_pass = sum(1 for r in results if r.get("overall_passed") is True)
    n_fail = sum(1 for r in results if r.get("overall_passed") is False)
    n_null = n_total - n_pass - n_fail  # errors / None

    # Per-criterion fail counts
    criterion_keys = [
        "c1_energy_drift", "c2_fwhm_stability", "c3_efficiency_qa",
        "c4_bg_drift", "c5_peak_z_roi", "c6_sensitivity",
    ]
    per_criterion_fail = {k: 0 for k in criterion_keys}
    for r in results:
        for k in criterion_keys:
            v = r.get("criteria_verdicts", {}).get(k, {})
            if v.get("passed") is False:
                per_criterion_fail[k] += 1

    # Per-detector pass/fail
    from collections import defaultdict, Counter
    by_detector: Dict[str, Counter] = defaultdict(Counter)
    by_geometry: Dict[str, Counter] = defaultdict(Counter)
    for r in results:
        det = r.get("detector_tag", "unknown")
        geo = r.get("geometry", "unknown")
        verdict = "pass" if r.get("overall_passed") is True else (
            "fail" if r.get("overall_passed") is False else "error"
        )
        by_detector[det][verdict] += 1
        by_geometry[geo][verdict] += 1

    # Top FWHM failures — extract fwhm_at_661_keV values
    fwhm_failures = []
    for r in results:
        c2 = r.get("criteria_verdicts", {}).get("c2_fwhm_stability", {})
        if c2.get("passed") is False:
            fwhm_failures.append({
                "spectrum_id": r["spectrum_id"],
                "detector_tag": r.get("detector_tag", ""),
                "geometry": r.get("geometry", ""),
                "fwhm_at_661_keV": c2.get("fwhm_at_661_keV"),
                "rel_deviation": c2.get("rel_deviation"),
                "fwhm_model_source": c2.get("fwhm_model_source", ""),
            })
    # sort by rel_deviation descending
    fwhm_failures.sort(key=lambda x: (x.get("rel_deviation") or 0), reverse=True)

    return {
        "n_total": n_total,
        "n_pass": n_pass,
        "n_fail": n_fail,
        "n_error": n_error,
        "n_null": n_null,
        "pass_rate": round(n_pass / n_total, 4) if n_total else None,
        "per_criterion_fail": per_criterion_fail,
        "by_detector": {k: dict(v) for k, v in by_detector.items()},
        "by_geometry": {k: dict(v) for k, v in by_geometry.items()},
        "top_fwhm_failures_sample": fwhm_failures[:20],
    }


# ════════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════════

def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="F-QC-01 sweep over SPECTRA_INDEX fixtures",
    )
    parser.add_argument(
        "--lsrm-root",
        default="C:\\LSRM",
        help="Root of the LSRM spectrum tree (default: C:\\LSRM)",
    )
    parser.add_argument(
        "--tier",
        default="A",
        help="Quality tier to sweep (default: A)",
    )
    parser.add_argument(
        "--output",
        default=str(
            _PROJECT_ROOT / "audit" / "_drafts"
            / f"f_qc_01_sweep_{date.today().isoformat()}.json"
        ),
        help="Output JSON path (default: audit/_drafts/f_qc_01_sweep_<DATE>.json)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=_WORKERS,
        help=f"Multiprocessing workers (default: {_WORKERS})",
    )
    parser.add_argument(
        "--smoke",
        type=int,
        default=0,
        help="If > 0, only process first N fixtures (smoke test)",
    )
    args = parser.parse_args(argv)

    lsrm_root = Path(args.lsrm_root)
    if not lsrm_root.exists():
        print(f"ERROR: LSRM root not found: {lsrm_root}", file=sys.stderr)
        return 1

    # --- Load SPECTRA_INDEX ---
    if not _SPECTRA_INDEX_PATH.exists():
        print(f"ERROR: SPECTRA_INDEX not found: {_SPECTRA_INDEX_PATH}", file=sys.stderr)
        return 1

    with open(_SPECTRA_INDEX_PATH, encoding="utf-8") as f:
        index_data = json.load(f)

    all_spectra: List[Dict[str, Any]] = index_data.get("spectra", [])
    tier_ids = set(index_data.get("indexes", {}).get("by_quality_tier", {}).get(args.tier, []))
    tier_entries = [s for s in all_spectra if s.get("spectrum_id") in tier_ids]

    if not tier_entries:
        print(f"ERROR: no entries for tier={args.tier}", file=sys.stderr)
        return 1

    if args.smoke > 0:
        tier_entries = tier_entries[: args.smoke]
        print(f"[SMOKE] limiting to {len(tier_entries)} fixtures", file=sys.stderr)

    n = len(tier_entries)
    print(f"Sweeping {n} tier-{args.tier} fixtures via {args.workers} workers ...",
          file=sys.stderr)

    # --- Worker pool ---
    worker_args = [(entry, str(lsrm_root)) for entry in tier_entries]

    if args.workers > 1 and n > 1:
        with Pool(processes=min(args.workers, n)) as pool:
            results = pool.map(_process_fixture, worker_args)
    else:
        results = [_process_fixture(a) for a in worker_args]

    errors = [r for r in results if r.get("read_error")]
    print(f"Done. {len(results)} processed, {len(errors)} errors.", file=sys.stderr)

    # --- Aggregate ---
    aggregate = _aggregate(results)

    # --- Build output ---
    output = {
        "_meta": {
            "schema_version": "1.0.0",
            "generated_at": date.today().isoformat(),
            "tier": args.tier,
            "n_fixtures": n,
            "lsrm_root": str(lsrm_root),   # intentionally NOT F-115 scrubbed here
            # (this file lives in audit/_drafts, not in a released report)
            "spectra_index_path": str(_SPECTRA_INDEX_PATH),
            "f_rule": "F-QC-01",
            "rag_id": "RAG-041",
            "cite": "spectrum_qc_methodology_v2_2026-06-03.md; KNOWN_AND_FIXED_ISSUES.md:1292",
        },
        "aggregate": aggregate,
        "fixtures": results,
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(output, f, ensure_ascii=False, indent=2)

    print(f"Output: {out_path}", file=sys.stderr)
    print(
        f"Pass={aggregate['n_pass']} Fail={aggregate['n_fail']} "
        f"Error={aggregate['n_error']} Rate={aggregate['pass_rate']}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
