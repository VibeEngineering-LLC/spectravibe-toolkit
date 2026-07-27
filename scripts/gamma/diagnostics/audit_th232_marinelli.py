"""
Audit Th-232 chain identification accuracy on Marinelli A-tier fixtures.

Filters SPECTRA_INDEX for `quality_tier='A' AND nuclide='Th-232' AND
geometry in {Маринелли, Дента-120мл, Петри-60, Дента-100}` (A-tier subset).

For each fixture runs ``analyze_lsrm_spe`` and builds a confusion matrix
per nuclide of the Th-232 chain (Ac-228, Pb-212, Bi-212, Tl-208, Th-228).

Output: JSON to stdout (or ``--out`` path). Optional ``--filter
geometry=Маринелли`` to subset.

Usage (PowerShell)::

    python scripts\\gamma\\diagnostics\\audit_th232_marinelli.py `
        --lsrm-root C:\\LSRM --out audit/_drafts/th232_audit.json

Anti-hallucination: each per-fixture verdict cites
``RAG-SPEC-NNNN`` + the actual list of identified nuclides from the
pipeline output. No expected/synthesised data.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import traceback
from pathlib import Path
from typing import Any, Iterable

# Make scripts/ importable when run as script
_ROOT = Path(__file__).resolve().parents[3]
_SCRIPTS = _ROOT / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))


# Expected Th-232-chain progeny that SHOULD be present for an A-tier
# Th-232 reference source (≥1 in the F-372 branching-corrected library
# AND in the canonical ENSDF ownership tables).
TH232_EXPECTED_PROGENY = {
    "Th-232",      # parent
    "Ac-228",      # 911/969/338 keV — strongest direct β decay
    "Pb-212",      # 238.6 keV — main
    "Bi-212",      # 727/1620 keV — weak but characteristic
    "Tl-208",      # 583/2614/860/510 keV — hallmark
    "Th-228",      # 84.4 keV (low-E, often unresolved)
}

TH232_CORE = {"Ac-228", "Pb-212", "Tl-208"}  # must-have minimum
TH232_FULL = {"Ac-228", "Pb-212", "Tl-208", "Bi-212"}


def _resolve_lsrm_path(rel_path: str, lsrm_root: Path) -> Path:
    return Path(str(rel_path).replace("<LSRM>", str(lsrm_root)))


def _filter_th232_a_tier(
    spectra_index: dict,
    geometry_substring: str | None = None,
) -> list[dict]:
    spectra = spectra_index["spectra"]
    by_id = {s["spectrum_id"]: s for s in spectra}
    a_tier = set(spectra_index["indexes"]["by_quality_tier"].get("A", []))
    th232 = set(spectra_index["indexes"]["by_nuclide"].get("Th-232", []))
    hits = [by_id[sid] for sid in (a_tier & th232) if sid in by_id]
    if geometry_substring:
        hits = [h for h in hits if geometry_substring.lower() in
                (h.get("geometry_normalized") or h.get("geometry") or "").lower()]
    return sorted(hits, key=lambda h: h["spectrum_id"])


def _run_pipeline(spe_path: Path) -> dict:
    """Run the staged identification pipeline and reduce its output to a
    flat dict of {spectrum_id, identified_nuclides, expected, hits, misses}."""
    from gamma.identification.staged_pipeline import analyze_lsrm_spe

    res = analyze_lsrm_spe(
        str(spe_path),
        detector_type="NaI",
        sigma_threshold=2.5,
        fwhm_window_multiple=0.5,
        # Stay default (no escalation, no activity calc) — pure ID audit
    )
    # final_detected is a list of NuclideIdentification with .name
    names = sorted({getattr(n, "name", None) or getattr(n, "nuclide", None)
                    for n in (res.final_detected or [])})
    names = [n for n in names if n]
    return {
        "identified_nuclides": names,
        "fwhm_at_661_keV": getattr(res, "fwhm_at_661", None),
        "seven_line_anchors": (
            getattr(res.seven_line_check, "n_anchors_matched", None)
            if res.seven_line_check is not None else None
        ),
        "is_background": getattr(res, "is_background", False),
        "detector_canonical": getattr(res, "detector_canonical", ""),
        "geometry_canonical": getattr(res, "geometry_canonical", ""),
    }


def _classify_th232_chain(identified: Iterable[str]) -> dict:
    """Return a per-progeny confusion-matrix-style record."""
    ident_set = set(identified)
    tp = sorted(TH232_EXPECTED_PROGENY & ident_set)
    fn = sorted(TH232_EXPECTED_PROGENY - ident_set)
    fp_set = ident_set - TH232_EXPECTED_PROGENY
    # Filter out chain non-progeny that are EXPECTED background (K-40,
    # Bi-214 backround, Pb-214) — keep them in FP list but mark them.
    background_natural = {"K-40", "Bi-214", "Pb-214", "Pb-210", "Ra-226"}
    fp_bg = sorted(fp_set & background_natural)
    fp_other = sorted(fp_set - background_natural)
    return {
        "true_positive_chain": tp,
        "false_negative_chain": fn,
        "false_positive_background": fp_bg,  # natural BG, not necessarily wrong
        "false_positive_unexpected": fp_other,  # real mis-IDs
        "n_core_hit": len(set(tp) & TH232_CORE),
        "n_core_expected": len(TH232_CORE),
        "core_match": set(tp) >= TH232_CORE,
        "full_chain_match": set(tp) >= TH232_FULL,
    }


def audit(lsrm_root: Path, spectra_index_path: Path,
          geometry_filter: str | None = None,
          max_fixtures: int | None = None) -> dict:
    with spectra_index_path.open("r", encoding="utf-8") as fh:
        idx = json.load(fh)

    fixtures = _filter_th232_a_tier(idx, geometry_filter)
    if max_fixtures:
        fixtures = fixtures[:max_fixtures]

    per_fixture: list[dict] = []
    for s in fixtures:
        sid = s["spectrum_id"]
        rel = s["rel_path"]
        spe = _resolve_lsrm_path(rel, lsrm_root)
        record: dict[str, Any] = {
            "spectrum_id": sid,
            "rel_path": rel,
            "geometry": s.get("geometry_normalized") or s.get("geometry"),
            "detector": s.get("detector_tag"),
            "sample_id": s.get("sample_id"),
            "passport_Bq": (s.get("passport") or [{}])[0].get("value_Bq"),
        }
        if not spe.exists():
            record["error"] = f"spe missing: {spe}"
            per_fixture.append(record)
            continue
        try:
            out = _run_pipeline(spe)
            record.update(out)
            record["chain_classification"] = _classify_th232_chain(
                out["identified_nuclides"]
            )
        except Exception as ex:
            record["error"] = f"{type(ex).__name__}: {ex}"
            record["traceback"] = traceback.format_exc(limit=4)
        per_fixture.append(record)

    # Aggregate
    n = len(per_fixture)
    ok = [r for r in per_fixture if "error" not in r]
    n_core = sum(1 for r in ok if r["chain_classification"]["core_match"])
    n_full = sum(1 for r in ok if r["chain_classification"]["full_chain_match"])
    # Per-nuclide hit rates
    nuc_hits: dict[str, int] = {n: 0 for n in TH232_EXPECTED_PROGENY}
    for r in ok:
        for tp in r["chain_classification"]["true_positive_chain"]:
            nuc_hits[tp] = nuc_hits.get(tp, 0) + 1
    summary = {
        "n_fixtures": n,
        "n_runs_ok": len(ok),
        "n_runs_errored": n - len(ok),
        "n_core_match (Ac-228+Pb-212+Tl-208)": n_core,
        "n_full_chain_match (+Bi-212)": n_full,
        "per_progeny_hit_rate": {
            nuc: f"{nuc_hits.get(nuc,0)}/{len(ok)}"
            for nuc in sorted(TH232_EXPECTED_PROGENY)
        },
        "geometry_filter": geometry_filter,
    }
    return {
        "_meta": {
            "lsrm_root": str(lsrm_root),
            "spectra_index": str(spectra_index_path),
            "audit_module": "scripts/gamma/diagnostics/audit_th232_marinelli.py",
        },
        "summary": summary,
        "per_fixture": per_fixture,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Audit Th-232 chain identification accuracy on A-tier fixtures."
    )
    parser.add_argument("--lsrm-root", type=Path, default=Path(r"C:\LSRM"))
    parser.add_argument(
        "--spectra-index",
        type=Path,
        default=_ROOT / "audit" / "_rag" / "SPECTRA_INDEX.json",
    )
    parser.add_argument("--geometry-filter", default=None,
                        help="Substring filter, e.g. 'Маринелли' or 'Дента'")
    parser.add_argument("--max-fixtures", type=int, default=None)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args(argv)

    if not args.lsrm_root.exists():
        print(f"ERROR: LSRM root not found: {args.lsrm_root}", file=sys.stderr)
        return 2
    if not args.spectra_index.exists():
        print(f"ERROR: SPECTRA_INDEX not found: {args.spectra_index}",
              file=sys.stderr)
        return 2

    result = audit(
        lsrm_root=args.lsrm_root,
        spectra_index_path=args.spectra_index,
        geometry_filter=args.geometry_filter,
        max_fixtures=args.max_fixtures,
    )

    text = json.dumps(result, indent=2, ensure_ascii=False, default=str)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text, encoding="utf-8")
        print(f"Wrote {args.out}")
    else:
        try:
            sys.stdout.reconfigure(encoding="utf-8")  # type: ignore[attr-defined]
        except Exception:
            pass
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
