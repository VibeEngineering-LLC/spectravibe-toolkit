#!/usr/bin/env python3
from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""F-331 / v1.18.18.5 — Build canonical reference-kit folder structure.

Creates `detectors/Gamma-1S/reference_spectra/reference_kits/` with a
{geometry}/{nuclide}/ hierarchy. Each terminal folder holds:
- `sample.spe`     — copy of the chosen exemplar spectrum
- `background.spe` — copy of the matching averaged background

Leftover files from the legacy detector subfolder tree are moved to
`reference_spectra/archive/` preserving relative paths. The legacy root
folder is then removed if empty.

NOTE: this script is a one-shot migration — its target legacy folder
(`reference_spectra/Gamma-1S_NaI_63x63_USB_SN-01`) has already been
retired in v1.18.18.5. Kept for provenance + idempotent re-execution
(no-op when folder absent).

Usage
-----
    python scripts/build_reference_kits.py
    python scripts/build_reference_kits.py --dry-run

Idempotent: re-running detects existing kits and skips them.
"""

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

REPO = Path(__file__).resolve().parent.parent
DET_ROOT = REPO / "detectors" / "Gamma-1S"
# Original legacy folder name — pre-v1.18.18.5. Migration is one-shot;
# kept here for provenance / idempotent re-run.
LEGACY_ROOT = (
    DET_ROOT / "reference_spectra" / "Gamma-1S_NaI_63x63_USB_SN-01"
)
KITS_ROOT = DET_ROOT / "reference_spectra" / "reference_kits"
ARCHIVE_ROOT = DET_ROOT / "reference_spectra" / "archive"
BG_ROOT = DET_ROOT / "data" / "averaged_backgrounds"


# ─── Kit specification ──────────────────────────────────────────────
# (geometry, nuclide, src_rel_to_LEGACY_ROOT, bg_rel_to_BG_ROOT)
KITS: List[Tuple[str, str, str, str]] = [
    # ── Marinelli 1 L (4 кэВ Бк/кг калибраторов с passport-данными
    #     в COMMENT поля .spe — F-330 auto-routing)
    ("Marinelli_1L", "Cs-137",
     "M_cs_легкий_2001-2005.spe",
     "bg_2016_marinelli_water_marinelli.spe"),
    ("Marinelli_1L", "K-40",
     "M_k_легкий_2001-2005.spe",
     "bg_2016_marinelli_water_marinelli.spe"),
    ("Marinelli_1L", "Ra-226",
     "M_ra_легкий_2001-2007.spe",
     "bg_2016_marinelli_water_marinelli.spe"),
    ("Marinelli_1L", "Th-232",
     "M_th_легкий_2001-2005.spe",
     "bg_2016_marinelli_water_marinelli.spe"),

    # ── Точечный источник 5 cm (комплект Поверка-2016 #SRC-05 + 2017/
    #     2019/2023 расширения для нуклидов отсутствующих в #SRC-05)
    ("Point_5cm", "Am-241",
     "Поверка-2016/Точка 5см/Am-241 42.13_Точечная-5см_5cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
    ("Point_5cm", "Ba-133",
     "Поверка-2016/Точка 5см/Ba-133 #SRC-05_Точечная-5см_5cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
    ("Point_5cm", "Bi-207",
     "Bi-207__176_04_2017_Точечная-5см_5cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
    ("Point_5cm", "Cd-109",
     "Поверка-2016/Точка 5см/Cd-109 #SRC-05_Точечная-5см_5cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
    ("Point_5cm", "Ce-139",
     "Поверка-2016/Точка 5см/Ce-139_591_Точечная-5см_5cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
    ("Point_5cm", "Co-57",
     "Поверка-2016/Точка 5см/Co-57 #SRC-05_Точечная-5см_5cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
    ("Point_5cm", "Co-60",
     "Поверка-2016/Точка 5см/Co-60 #SRC-05_Точечная-5см_5cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
    ("Point_5cm", "Cs-137",
     "Поверка-2016/Точка 5см/Cs-137 #SRC-05_Точечная-5см_5cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
    ("Point_5cm", "Eu-152",
     "Поверка-2016/Точка 5см/Eu-152 #SRC-05_Точечная-5см_5cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
    ("Point_5cm", "Mn-54",
     "Поверка-2016/Точка 5см/Mn-54_587_Точечная-5см_5cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
    ("Point_5cm", "Na-22",
     "Поверка-2016/Точка 5см/Na-22_585_Точечная-5см_5cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
    ("Point_5cm", "Th-228",
     "Поверка-2016/Точка 5см/Th-228 #SRC-05_Точечная-5см_5cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
    ("Point_5cm", "Y-88",
     "Поверка-2016/Точка 5см/Y-88_589_Точечная-5см_5cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
    ("Point_5cm", "Zn-65",
     "Zn-65__342_2019_Точечная-5см_5cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),

    # ── Точечный источник 25 cm
    ("Point_25cm", "Am-241",
     "Поверка-2016/Точка 25см/Am-241 42.13_Точечная-25см_25cm.spe",
     "bg_2016_open_lid_point25cm.spe"),
    ("Point_25cm", "Ba-133",
     "Поверка-2016/Точка 25см/Ba-133 #SRC-05_Точечная-25см_25cm.spe",
     "bg_2016_open_lid_point25cm.spe"),
    ("Point_25cm", "Cd-109",
     "Поверка-2016/Точка 25см/Cd-109 #SRC-05_Точечная-25см_25cm.spe",
     "bg_2016_open_lid_point25cm.spe"),
    ("Point_25cm", "Ce-139",
     "Поверка-2016/Точка 25см/Ce-139_591_Точечная-25см_25cm.spe",
     "bg_2016_open_lid_point25cm.spe"),
    ("Point_25cm", "Co-60",
     "Поверка-2016/Точка 25см/Co-60 #SRC-05_Точечная-25см_25cm.spe",
     "bg_2016_open_lid_point25cm.spe"),
    ("Point_25cm", "Cs-137",
     "Точечная-25см/Cs-137 №SRC-02_Точечная-25см_25cm.spe",
     "bg_2016_open_lid_point25cm.spe"),
    ("Point_25cm", "Eu-152",
     "Поверка-2016/Точка 25см/Eu-152 #SRC-05_Точечная-25см_25cm.spe",
     "bg_2016_open_lid_point25cm.spe"),
    ("Point_25cm", "Mn-54",
     "Поверка-2016/Точка 25см/Mn-54_587_Точечная-25см_25cm.spe",
     "bg_2016_open_lid_point25cm.spe"),
    ("Point_25cm", "Na-22",
     "Точечная-25см/Na-22 #01.22_Точечная-25см_25cm.spe",
     "bg_2016_open_lid_point25cm.spe"),
    ("Point_25cm", "Th-228",
     "Точечная-25см/Th-228 №309_Точечная-25см_25cm.spe",
     "bg_2016_open_lid_point25cm.spe"),
    ("Point_25cm", "Y-88",
     "Точечная-25см/Y-88 №SRC-06_Точечная-25см_25cm.spe",
     "bg_2016_open_lid_point25cm.spe"),

    # ── Петри 60 мл (специфический empty-shield bg не подобран;
    #     используем generic empty-shield-point5cm как ближайший)
    ("Petri_60mL", "Cs-137",
     "Петри-60мл/Cs137_420-7-14_Петри-60мл_0cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
    ("Petri_60mL", "K-40",
     "Петри-60мл/K40_420-7-20_Петри-60мл_0cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
    ("Petri_60mL", "Ra-226",
     "Петри-60мл/Ra226_420-7-18_Петри-60мл_0cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
    ("Petri_60mL", "Th-232",
     "Петри-60мл/Th232_420-7-17_Петри-60мл_0cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),

    # ── Дента 120 мл
    ("Denta_120mL", "Cs-137",
     "Дента-120мл/Cs137_420-7-14_Дента-120мл_0cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
    ("Denta_120mL", "K-40",
     "Дента-120мл/K40_420-7-20_Дента-120мл_0cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
    ("Denta_120mL", "Ra-226",
     "Дента-120мл/Ra226_420-7-18_Дента-120мл_0cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
    ("Denta_120mL", "Th-232",
     "Дента-120мл/Th232_420-7-17_Дента-120мл_0cm.spe",
     "bg_2016_empty_shield_point5cm.spe"),
]


def _short_hash(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()[:10]


def build_kits(dry_run: bool = False) -> Dict[str, object]:
    """Copy chosen spectra into kits/{geometry}/{nuclide}/ structure."""
    KITS_ROOT.mkdir(parents=True, exist_ok=True)
    used_sources: set[Path] = set()
    rows: List[Dict] = []
    skipped: List[Dict] = []
    for geometry, nuc, src_rel, bg_name in KITS:
        src = LEGACY_ROOT / src_rel
        bg = BG_ROOT / bg_name
        if not src.is_file():
            skipped.append({
                "geometry": geometry, "nuclide": nuc,
                "reason": f"missing source: {src_rel}",
            })
            continue
        if not bg.is_file():
            skipped.append({
                "geometry": geometry, "nuclide": nuc,
                "reason": f"missing bg: {bg_name}",
            })
            continue
        used_sources.add(src.resolve())

        leaf = KITS_ROOT / geometry / nuc
        sample_target = leaf / f"sample_{src.name}"
        bg_target = leaf / f"background_{bg.name}"

        action = "create"
        if leaf.exists() and sample_target.exists() and bg_target.exists():
            action = "skip-exists"
        else:
            if not dry_run:
                leaf.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, sample_target)
                shutil.copy2(bg, bg_target)

        rows.append({
            "geometry": geometry,
            "nuclide": nuc,
            "source_path": str(src_rel),
            "sample_kit_path": str(sample_target.relative_to(REPO)),
            "background_kit_path": str(bg_target.relative_to(REPO)),
            "source_md5_short": _short_hash(src),
            "bg_md5_short": _short_hash(bg),
            "action": action,
        })

    return {"kits": rows, "skipped": skipped, "used_sources": used_sources}


def archive_leftovers(used_sources: set[Path], dry_run: bool = False) -> Dict:
    """Move everything in LEGACY_ROOT not chosen for a kit to ARCHIVE_ROOT,
    preserving relative paths."""
    ARCHIVE_ROOT.mkdir(parents=True, exist_ok=True)
    moved: List[Dict] = []
    kept: List[Dict] = []
    for f in sorted(LEGACY_ROOT.rglob("*")):
        if f.is_dir():
            continue
        # Always preserve desktop.ini at its sibling location
        if f.name.lower() == "desktop.ini":
            continue
        rel = f.relative_to(LEGACY_ROOT)
        target = ARCHIVE_ROOT / rel
        if f.resolve() in used_sources:
            kept.append({"path": str(rel), "reason": "selected for kit"})
            continue
        if target.exists():
            # Already moved on prior run; remove duplicate from legacy
            if not dry_run:
                f.unlink()
            moved.append({"src": str(rel), "dst": str(target.relative_to(REPO)),
                          "action": "removed-duplicate-of-archive"})
            continue
        if not dry_run:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(f), str(target))
        moved.append({"src": str(rel), "dst": str(target.relative_to(REPO)),
                      "action": "moved"})
    return {"moved": moved, "kept": kept}


def remove_used_sources_from_legacy(used_sources: set[Path], dry_run: bool = False) -> int:
    """After kits are copied, drop the originals from LEGACY_ROOT so the
    legacy tree can be retired (kits hold the canonical копии)."""
    removed = 0
    for src in used_sources:
        if src.is_file():
            if not dry_run:
                src.unlink()
            removed += 1
    return removed


def prune_empty_dirs(root: Path, dry_run: bool = False) -> int:
    """Remove empty directories bottom-up. Ignores desktop.ini files
    (Windows artefact) which are treated as virtually empty."""
    removed = 0
    for d in sorted(
        (p for p in root.rglob("*") if p.is_dir()),
        key=lambda p: -len(p.parts),
    ):
        children = [c for c in d.iterdir() if c.name.lower() != "desktop.ini"]
        if not children:
            for c in d.iterdir():
                if not dry_run:
                    c.unlink()
            if not dry_run:
                try:
                    d.rmdir()
                except OSError:
                    pass
            removed += 1
    return removed


def write_readme(rows: List[Dict]) -> None:
    """Generate README.md inside reference_kits/ with the kit manifest."""
    lines = [
        "# Reference Kits — Gamma-1S NaI 63×63 USB № SN-01",
        "",
        "> F-331 / v1.18.18.5 — Канонические комплекты «образец + фон»",
        "> для использования в регрессионных и приёмочных тестах.",
        "",
        "## Структура",
        "",
        "```",
        "reference_kits/",
        "├── {geometry}/",
        "│   └── {nuclide}/",
        "│       ├── sample_{...}.spe        — спектр эталонного источника",
        "│       └── background_{...}.spe    — релевантный усреднённый фон",
        "```",
        "",
        "Все остальные исторические спектры из легаси-папки",
        "`Gamma-1S_NaI_63x63_USB_SN-01/` перенесены в",
        "`detectors/Gamma-1S/reference_spectra/archive/` с сохранением",
        "относительной структуры.",
        "",
        "## Манифест комплектов",
        "",
        "| Geometry | Nuclide | Sample (md5) | Background (md5) |",
        "|---|---|---|---|",
    ]
    for r in rows:
        sname = Path(r["sample_kit_path"]).name
        bname = Path(r["background_kit_path"]).name
        lines.append(
            f"| {r['geometry']} | {r['nuclide']} | "
            f"`{sname}` ({r['source_md5_short']}) | "
            f"`{bname}` ({r['bg_md5_short']}) |"
        )
    lines.extend([
        "",
        "## Использование в тестах",
        "",
        "```python",
        "from pathlib import Path",
        "KIT = Path('detectors/Gamma-1S/reference_spectra/reference_kits')",
        "from gamma.reporting import analyze_and_report",
        "",
        "# Marinelli 1L Cs-137 kit",
        "kit = KIT / 'Marinelli_1L' / 'Cs-137'",
        "sample = next(kit.glob('sample_*.spe'))",
        "bg = next(kit.glob('background_*.spe'))",
        "",
        "artefacts = analyze_and_report(",
        "    str(sample),",
        "    output_dir='./out',",
        "    background_path=str(bg),",
        "    sample_mass_kg=0.570,",
        ")",
        "```",
        "",
        "Passport activities для Marinelli автоматически читаются из",
        "COMMENT-поля .spe (F-330 v1.18.18.4) — `passport_activity_Bq`",
        "передавать не нужно.",
    ])
    (KITS_ROOT / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    if not LEGACY_ROOT.exists():
        print(f"ERROR: legacy root not found: {LEGACY_ROOT}", file=sys.stderr)
        return 1
    if not BG_ROOT.exists():
        print(f"ERROR: bg root not found: {BG_ROOT}", file=sys.stderr)
        return 1

    print(f"[F-331] Building reference kits …")
    res = build_kits(dry_run=args.dry_run)
    rows = res["kits"]
    used_sources = res["used_sources"]
    print(f"  {len(rows)} kits processed, {len(res['skipped'])} skipped")
    for s in res["skipped"]:
        print(f"  SKIPPED: {s['geometry']} {s['nuclide']} — {s['reason']}")

    if not args.dry_run:
        write_readme(rows)

    print(f"[F-331] Archiving leftovers …")
    arch = archive_leftovers(used_sources, dry_run=args.dry_run)
    print(f"  moved {len(arch['moved'])} non-kit files to archive/")

    print(f"[F-331] Removing kit-source copies from legacy tree …")
    rm = remove_used_sources_from_legacy(used_sources, dry_run=args.dry_run)
    print(f"  removed {rm} source files (now live in kits/)")

    print(f"[F-331] Pruning empty legacy dirs …")
    pruned = prune_empty_dirs(LEGACY_ROOT, dry_run=args.dry_run)
    print(f"  pruned {pruned} empty dirs")
    # Try to remove the legacy root itself if empty
    if LEGACY_ROOT.exists() and not args.dry_run:
        remaining = [c for c in LEGACY_ROOT.iterdir()
                     if c.name.lower() != "desktop.ini"]
        if not remaining:
            for c in LEGACY_ROOT.iterdir():
                c.unlink()
            try:
                LEGACY_ROOT.rmdir()
                print(f"  removed empty legacy root: {LEGACY_ROOT.name}")
            except OSError as e:
                print(f"  legacy root not removed: {e}")

    # Manifest dump
    manifest = {
        "version": "v1.18.18.5",
        "detector": "Gamma-1S_NaI_63x63_USB_SN-01",
        "kits": rows,
        "archived_files": len(arch["moved"]),
        "removed_kit_sources": rm,
    }
    if not args.dry_run:
        (KITS_ROOT / "MANIFEST.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
    print(f"\n[F-331] Done. Kits root: {KITS_ROOT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
