#!/usr/bin/env python3
"""Project state snapshot — Agent D orchestration tool.

Writes _state/project_state.json with machine-readable status:
- SKILL_VERSION, latest HANDOFF
- regression baseline (from HANDOFF, not live)
- open-indicator checks (C-T4-* housekeeping, Tier 1 plan coverage)
- file inventory stats

Purpose: agents A/B/C read THIS file instead of grep'ing the codebase.
Run once per phase. Cost: <1 second, no network, stdlib only.

Usage:
    python scripts/orchestration/snapshot.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from datetime import datetime

ROOT = Path(__file__).resolve().parents[2]
STATE_DIR = ROOT / "_state"
STATE_DIR.mkdir(exist_ok=True)


def skill_version() -> str | None:
    p = ROOT / "scripts" / "gamma" / "reporting" / "json_report.py"
    if not p.exists():
        return None
    m = re.search(r'SKILL_VERSION\s*=\s*"([^"]+)"', p.read_text(encoding="utf-8"))
    return m.group(1) if m else None


def latest_handoff() -> dict | None:
    files = sorted(ROOT.glob("HANDOFF_v*.md"), reverse=True)
    if not files:
        return None
    f = files[0]
    return {
        "name": f.name,
        "size_kb": f.stat().st_size // 1024,
        "mtime": datetime.fromtimestamp(f.stat().st_mtime).isoformat(timespec="seconds"),
    }


def known_issues() -> dict:
    p = ROOT / "KNOWN_AND_FIXED_ISSUES.md"
    if not p.exists():
        return {"exists": False}
    text = p.read_text(encoding="utf-8")
    return {
        "exists": True,
        "size_kb": p.stat().st_size // 1024,
        "lines": text.count("\n"),
        "rotation_needed": p.stat().st_size > 100 * 1024,  # >100 KB
    }


def desktop_ini_in_git() -> int:
    """C-T4-01 indicator: GDrive pollution count in .git/. -1 = no git repo."""
    git = ROOT / ".git"
    if not git.exists():
        return -1
    return sum(1 for _ in git.rglob("desktop.ini"))


def stale_handoffs() -> list[str]:
    """C-T4-12 indicator: all HANDOFFs except the latest."""
    files = sorted(ROOT.glob("HANDOFF_v*.md"))
    return [f.name for f in files[:-1]] if len(files) > 1 else []


def stale_v1_17_zips() -> dict:
    """C-T4-03 indicator: stale v1.17.x archives in ../../1_Version/."""
    versions = ROOT.parent / "1_Version"
    if not versions.exists():
        # try ../../1_Version relative to project (typical layout)
        versions = ROOT.parent.parent / "1_Version"
    if not versions.exists():
        return {"exists": False}
    zips = list(versions.glob("SpectraVibe_v1.17*.zip"))
    return {
        "exists": True,
        "count": len(zips),
        "total_mb": sum(z.stat().st_size for z in zips) // (1024 * 1024),
        "names": [z.name for z in zips[:6]],  # first 6 for brevity
    }


def gitignore_status() -> dict:
    p = ROOT / ".gitignore"
    expected = ("desktop.ini", "Thumbs.db", ".DS_Store", "*.swp", ".idea/", ".vscode/")
    if not p.exists():
        return {"exists": False, "missing": list(expected)}
    text = p.read_text(encoding="utf-8")
    missing = [pat for pat in expected if pat not in text]
    return {"exists": True, "missing": missing}


def tier1_plans() -> dict:
    """Tier 1 audit plan coverage (15 HIGH F-rules from v1.17.9.5 ORTEC delta)."""
    plans_dir = ROOT / "audit" / "_plans"
    tier1 = [
        "F-167", "F-168", "F-169-REV", "F-170", "F-171", "F-172",
        "F-173", "F-174", "F-176", "F-178", "F-179", "F-180",
        "F-241", "F-242", "F-243",
    ]
    if not plans_dir.exists():
        return {"dir_exists": False, "total": len(tier1), "found": 0, "files": []}
    files = sorted(plans_dir.glob("F-*.md"))
    found = []
    for fp in files:
        for tid in tier1:
            stem_check = tid.replace("-REV", "_REV")
            if fp.stem.upper().startswith(tid.upper()) or fp.stem.upper().startswith(stem_check.upper()):
                found.append(fp.name)
                break
    return {
        "dir_exists": True,
        "total": len(tier1),
        "found_count": len(found),
        "found_files": found,
        "missing_ids": [t for t in tier1 if not any(t.upper() in f.upper() for f in found)],
    }


def file_stats() -> dict:
    g = ROOT / "scripts" / "gamma"
    t = ROOT / "tests"
    return {
        "scripts_gamma_py": sum(1 for _ in g.rglob("*.py")) if g.exists() else None,
        "tests_py": sum(1 for _ in t.rglob("*.py")) if t.exists() else None,
        "audit_plans_md": sum(1 for _ in (ROOT / "audit" / "_plans").rglob("*.md"))
            if (ROOT / "audit" / "_plans").exists() else 0,
    }


def main() -> int:
    state = {
        "schema": "agent_d/state/v1",
        "generated_by": "scripts/orchestration/snapshot.py",
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "baseline": {
            "skill_version": skill_version(),
            "handoff_latest": latest_handoff(),
        },
        "regression_baseline_from_handoff": {
            "_note": "Static from HANDOFF_v1.18.28; run regression.py for live counts.",
            "expected_pass": 1525,
            "expected_fail": 1,
            "expected_xfail": 1,
            "expected_xpass": 2,
            "known_fail_test": "test_th232_demo_chain_ratio_improves_with_bi212_tcs (TD-2)",
        },
        "open_indicators": {
            "TD-2_status": "open (pre-existing, root: nuclide_library._CACHE leak per commit 6d65d6b)",
            "C-T4-01_desktop_ini_in_git_count": desktop_ini_in_git(),
            "C-T4-11_known_issues": known_issues(),
            "C-T4-12_stale_handoffs": stale_handoffs(),
            "C-T4-03_stale_v1_17_zips": stale_v1_17_zips(),
            "C-T4-02_gitignore": gitignore_status(),
            "F-392.1_json_export": "pending — add step_anchor_energies_keV + step_intensity_pct to _build_continuum_block",
            "F-389.1_xfail_cleanup": "pending — remove 2 XPASS @xfail decorators in test_f389_v2_activity_parity.py",
        },
        "tier1_audit_plans": tier1_plans(),
        "stats": file_stats(),
    }
    out = STATE_DIR / "project_state.json"
    out.write_text(json.dumps(state, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"OK wrote {out.relative_to(ROOT)} ({out.stat().st_size} bytes)")
    print(f"  baseline: SKILL_VERSION={state['baseline']['skill_version']}")
    print(f"  HANDOFF:  {state['baseline']['handoff_latest']['name'] if state['baseline']['handoff_latest'] else 'NONE'}")
    print(f"  .git/desktop.ini count: {state['open_indicators']['C-T4-01_desktop_ini_in_git_count']}")
    ki = state["open_indicators"]["C-T4-11_known_issues"]
    if ki.get("exists"):
        print(f"  KNOWN_AND_FIXED_ISSUES.md: {ki['size_kb']} KB (rotation_needed={ki['rotation_needed']})")
    t1 = state["tier1_audit_plans"]
    print(f"  Tier 1 plans: {t1.get('found_count', 0)}/{t1.get('total', 15)} (missing: {len(t1.get('missing_ids', []))})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
