#!/usr/bin/env python
from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""lsrm_distribute.py — multi-detector LSRM ingest pipeline (v1.0).

Distributes LSRM source tree across 11+ detector folders per
detectors/Gamma-1S/README.md §9 crystal-class map source-pins.

Modes:
    --dry-run     : print plan, do NOT copy/create anything
    --apply       : execute copies with per-file checksum verify

Version: v1.0 (promoted 2026-06-06 from DRAFT after verified apply run).
Baseline manifest: _tmp/lsrm_source_sha256_pre_20260606.json
Cross-ref: detectors/Gamma-1S/README.md §9 (crystal-class source-pin map)

Tier-2-LOCAL only. No commit. No push. Idempotent (re-runs skip already-present).
"""
import argparse
import hashlib
import json
import shutil
import sys
import time
from pathlib import Path


# === Source-pin → target detector mapping (per README §9 + 2026-06-06 inventory) ===
# Each rule: (source-pin tuple of path-parts, target detector_id).
# Order matters: more specific pins first.
MAPPING = [
    # Gamma-1S (BG/Gamma-1S → omoглыф lock to Gamma-1S)
    (("BG", "Gamma-1S"), "Gamma-1S"),
    # GP HPGe20%
    (("GP", "HPGe(20%)"), "GP_HPGe20"),
    # NM (Nuclear Materials) — all three sub-streams on NM_HPGe20 station
    (("NM", "HPGe(20%)"), "NM_HPGe20"),
    (("NM", "Pu"), "NM_HPGe20"),
    (("NM", "U"), "NM_HPGe20"),
    # Handy
    (("Handy", "Handy(HPGe)"), "Handy_HPGe"),
    (("Handy", "Handy(LaBr)"), "Handy_LaBr"),
    (("Handy", "Handy(NaI)"), "Handy_NaI"),
    # Simple demos
    (("Simple", "Alpha(Demo)"), "Simple_Alpha"),
    (("Simple", "HPGe(Demo)"), "Simple_HPGe"),
    (("Simple", "NaI(Demo)"), "Simple_NaI"),
    (("Simple", "SiLi(Demo)"), "Simple_SiLi"),
    (("Simple", "TeCd(Demo)"), "Simple_TeCd"),
    # ADA — NEW skeleton (B1)
    (("ADA", "AlphaDuoSmall"), "ADA_AlphaDuoSmall"),
    # Beta-1S — NEW skeleton (B2, pending operator decision)
    (("BG", "Beta-1S"), "Beta-1S"),
]

# Root-level NM/*.lib files belong to NM_HPGe20 (per inventory: NM = Nuclear Materials campaign)
NM_LIB_FILES = {"Rn222 в угле.lib", "Осколки деления.lib", "Уран.lib"}

# Root meta-files (uncategorized — placed in _meta_unsorted/ for later sorting)
ROOT_META = {"Th-232.enx"}

# Detectors that exist (F-300 W2 + Gamma-1S). Others need skeleton creation.
EXISTING_DETECTORS = {
    "Gamma-1S", "GP_HPGe20", "Handy_HPGe", "Handy_LaBr", "Handy_NaI",
    "Simple_Alpha", "Simple_HPGe", "Simple_NaI", "Simple_SiLi", "Simple_TeCd",
    "AtomSpectra",
}

NEW_SKELETONS = {"ADA_AlphaDuoSmall", "NM_HPGe20", "Beta-1S"}


def classify(rel_path: Path) -> tuple[str, str]:
    """Return (target_detector_id, sub_path_within_raw_lsrm_Work_or_None).

    Returns ("_meta_unsorted", rel_path_str) for files that don't match any pin.
    """
    parts = rel_path.parts
    # Root-level NM/*.lib files
    if len(parts) == 2 and parts[0] == "NM" and parts[1] in NM_LIB_FILES:
        return ("NM_HPGe20", str(rel_path))
    # Root meta files
    if len(parts) == 1 and parts[0] in ROOT_META:
        return ("_meta_unsorted", parts[0])
    # Standard 2-part pin matching
    if len(parts) >= 2:
        pin = (parts[0], parts[1])
        for rule_pin, det_id in MAPPING:
            if pin == rule_pin:
                return (det_id, str(rel_path))
    return ("_meta_unsorted", str(rel_path))


def sha256_of(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def ensure_detector_skeleton(repo_root: Path, det_id: str, dry_run: bool):
    """Create minimal skeleton for new detector (mirror Gamma-1S pattern lite)."""
    det_dir = repo_root / "detectors" / det_id
    raw_lsrm = det_dir / "raw_lsrm" / "Work"
    if dry_run:
        if not det_dir.exists():
            print(f"  [DRY] would create skeleton: detectors/{det_id}/ (raw_lsrm/Work/, README.md placeholder)")
        return
    raw_lsrm.mkdir(parents=True, exist_ok=True)
    readme = det_dir / "README.md"
    if not readme.exists():
        readme.write_text(
            f"# detectors/{det_id}/\n\n"
            f"Detector station skeleton. Created automatically by "
            f"distribute_lsrm_tree.py on 2026-06-06 LSRM ingest. "
            f"Mirror full detectors/Gamma-1S/README.md pattern when populated.\n\n"
            f"raw_lsrm/Work/ contains operator's LSRM SpectraLine source spectra.\n"
            f"This folder is gitignored — local-only working copy.\n",
            encoding="utf-8",
        )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True, help="LSRM Work tree root")
    ap.add_argument("--repo", required=True, help="gamma-spectrum-analysis repo root")
    ap.add_argument("--baseline-manifest", required=True, help="Pre-computed source sha256 manifest")
    ap.add_argument("--mode", choices=["dry-run", "apply"], default="dry-run")
    ap.add_argument("--include-beta", action="store_true", help="Include BG/Beta-1S (default skip pending operator decision)")
    ap.add_argument("--out", required=True, help="Output post-distribute manifest")
    args = ap.parse_args()

    source = Path(args.source)
    repo = Path(args.repo)
    manifest_in = json.loads(Path(args.baseline_manifest).read_text(encoding="utf-8"))

    plan = []  # list of (src_rel, dst_rel, det_id, sha256)
    skipped = []
    new_skel_needed = set()

    for rel_str, info in manifest_in["files"].items():
        rel = Path(rel_str)
        det_id, sub_path = classify(rel)
        if det_id == "Beta-1S" and not args.include_beta:
            skipped.append((rel_str, "beta-skipped-by-default"))
            continue
        if det_id == "_meta_unsorted":
            dst_rel = Path("_inbox_raw_lsrm") / "_meta_unsorted" / sub_path
        else:
            dst_rel = Path("detectors") / det_id / "raw_lsrm" / "Work" / sub_path
            if det_id not in EXISTING_DETECTORS:
                new_skel_needed.add(det_id)
        plan.append({
            "source_rel": rel_str,
            "target_rel": str(dst_rel).replace("\\", "/"),
            "detector_id": det_id,
            "sha256_source": info["sha256"],
            "size": info["size"],
        })

    print(f"mode: {args.mode}")
    print(f"plan size: {len(plan)}")
    print(f"skipped: {len(skipped)}")
    print(f"new skeletons needed: {sorted(new_skel_needed)}")
    print()
    per_det = {}
    for entry in plan:
        per_det.setdefault(entry["detector_id"], 0)
        per_det[entry["detector_id"]] += 1
    for det, n in sorted(per_det.items(), key=lambda x: -x[1]):
        print(f"  {det}: {n} files")
    print()

    if args.mode == "dry-run":
        # Save plan only
        out = {
            "mode": "dry-run",
            "source_root": str(source),
            "plan_size": len(plan),
            "skipped": skipped,
            "new_skeletons_needed": sorted(new_skel_needed),
            "per_detector": per_det,
            "plan": plan,
        }
        Path(args.out).write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"dry-run plan saved: {args.out}")
        return 0

    # === apply mode ===
    for det_id in new_skel_needed:
        ensure_detector_skeleton(repo, det_id, dry_run=False)

    copied = []
    verified = []
    failed = []
    t0 = time.time()

    for entry in plan:
        src = source / entry["source_rel"]
        dst = repo / entry["target_rel"]
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.exists():
            # Idempotency: compare sha256 — if matches, skip
            existing_sha = sha256_of(dst)
            if existing_sha == entry["sha256_source"]:
                verified.append({"target": entry["target_rel"], "status": "already-present-match"})
                continue
            else:
                failed.append({"target": entry["target_rel"], "reason": "exists-with-different-content"})
                continue
        try:
            shutil.copy2(src, dst)
            actual_sha = sha256_of(dst)
            if actual_sha != entry["sha256_source"]:
                failed.append({
                    "target": entry["target_rel"],
                    "reason": "checksum-mismatch-after-copy",
                    "expected": entry["sha256_source"],
                    "actual": actual_sha,
                })
            else:
                copied.append({"target": entry["target_rel"], "sha256": actual_sha})
                verified.append({"target": entry["target_rel"], "status": "fresh-copy-verified"})
        except (OSError, PermissionError) as e:
            failed.append({"target": entry["target_rel"], "reason": f"copy-error: {e}"})

    out = {
        "mode": "apply",
        "source_root": str(source),
        "repo_root": str(repo),
        "elapsed_s": time.time() - t0,
        "copied_count": len(copied),
        "verified_count": len(verified),
        "failed_count": len(failed),
        "skipped": skipped,
        "new_skeletons_created": sorted(new_skel_needed),
        "per_detector": per_det,
        "failures": failed,
    }
    Path(args.out).write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"copied: {len(copied)}, verified: {len(verified)}, failed: {len(failed)}, elapsed: {time.time()-t0:.1f}s")
    print(f"manifest saved: {args.out}")
    return 0 if not failed else 2


if __name__ == "__main__":
    sys.exit(main())
