from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""
F-43 (v1.7.21) / F-44 (v1.7.22) — Generate canonical long-exposure
background spectra from the 2016 + 2024 Поверка archives.

Walks the registered background-context folders, groups by geometry,
and writes one canonical .spe per (background_context, geometry) into
`data/averaged_backgrounds/`.

**F-44 correction**: The original F-43 generator naively summed every
file in each context folder, but the LSRM Spectraline acquisition
software emits CUMULATIVE checkpoint files (`..._01.spe` = 1h, ..._02
= 2h cumulative, ..._N = N h cumulative). Summing them inflates both
counts and live_time by ~N/2 while preserving the rate; σ-claim was
therefore false (σ_red ∝ √N not realised). F-44 routes through
`average_lsrm_spectra` which auto-detects the cumulative pattern via
`detect_cumulative_pattern` and switches to "longest-file" semantics
(σ-reduction = 1.0, but counts/live-time correctly represent ONE
independent long-exposure measurement).

This script is run-once at packaging time; the resulting .spe files
are checked into `data/averaged_backgrounds/` so downstream code can
just read them via the standard `read_spectrum` API. Re-run after the
archive is updated.

Usage:
    PYTHONPATH=scripts python build_averaged_backgrounds.py [--dry-run]

Output goes to `data/averaged_backgrounds/` relative to the project
root. The script prints a summary table at the end.
"""


import argparse
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from gamma.io.average_lsrm import (  # noqa: E402
    average_lsrm_spectra, write_lsrm_spe,
    CalibrationMismatchError, IdentityMismatchError,
)
from gamma.io.readers import read_spectrum  # noqa: E402


# ============================================================================
# Discovery: groups in the 2016 archive
# ============================================================================

# F-83 (v1.12.0): both inputs and outputs are Gamma-1S-specific assets.
from gamma.detectors.gamma1s import (
    DEFAULT_REFERENCE_DIR as PROJECT_REF_BASE,
    AVERAGED_BACKGROUNDS_DIR as OUT_DIR,
)

# Background contexts: each entry is a folder under PROJECT_REF_BASE.
# `key` becomes part of the output filename. `pair_with_geometries`
# documents which sample geometries the resulting bg pairs with
# (informational — used in the MANIFEST and provenance).
#
# Per LSRM operator convention:
#   • Marinelli sample bg → Marinelli vessel filled with water
#     (matrix-attenuated bg matching sample attenuation).
#   • Дента-120мл, Чашка-Петри 60мл, Точечная-5см sample bg → empty
#     shielding chamber, closed lid (low ambient γ).
#   • Точечная-25см sample bg → open lid (no shielding above 25cm
#     sample position — direct ambient γ).
#
# Per F-43/F-44, each context folder is a CUMULATIVE LSRM checkpoint
# set; `average_lsrm_spectra` auto-detects this and returns the
# longest single file as the canonical bg.

BG_CONTEXTS = [
    # ---- 2016 archive ----
    {
        "key": "2016_marinelli_water",
        "subdir": "Поверка-2016/Фон вода",
        "description": ("2016 Marinelli + water bg "
                        "(matrix-matched for Marinelli sample geometry)"),
        "pair_with_geometries": ["Маринелли"],
        "min_files_per_geometry": 5,
    },
    {
        "key": "2016_empty_shield",
        "subdir": "Поверка-2016/фон пустая защита",
        "description": ("2016 empty shielding bg "
                        "(for Дента / Чашка-60 / Точечная-5см sample geometries)"),
        "pair_with_geometries": ["Дента-120мл", "Чашка Петри 60мл",
                                 "Точечная-5см"],
        "min_files_per_geometry": 5,
    },
    {
        "key": "2016_open_lid",
        "subdir": "Поверка-2016/Фон с открытыми крышками",
        "description": ("2016 open lid bg "
                        "(for Точечная-25см sample geometry)"),
        "pair_with_geometries": ["Точечная-25см"],
        "min_files_per_geometry": 5,
    },
    # ---- 2024 archive (F-44 additions) ----
    {
        "key": "2024_marinelli_water_closed_lid",
        "subdir": "Поверка-2024/Фон закр кр",
        "description": ("2024 Marinelli + water + closed lid bg "
                        "(matrix-matched for Marinelli samples; "
                        "newer epoch closer to 2023+ certs)"),
        "pair_with_geometries": ["Маринелли"],
        "min_files_per_geometry": 5,
    },
    {
        "key": "2024_open_lid",
        "subdir": "Поверка-2024/Фон откр кр",
        "description": ("2024 open lid bg "
                        "(for Точечная-25см samples; newer epoch)"),
        "pair_with_geometries": ["Точечная-25см"],
        "min_files_per_geometry": 5,
    },
]


def _short_geo_name(geometry: str) -> str:
    """Map Cyrillic geometry name to ASCII filename token."""
    g = (geometry or "unknown").strip().lower()
    table = {
        "маринелли": "marinelli",
        "точечная-5см": "point5cm",
        "точечная-25см": "point25cm",
        "дента-100": "denta100",
        "дента-120мл": "denta120",
        "петри-60": "petri60",
        "чашка-60": "petri60",
    }
    for k, v in table.items():
        if k in g:
            return v
    # Fallback: strip non-ASCII and collapse spaces
    safe = "".join(c if c.isalnum() else "_" for c in g)
    return safe.strip("_") or "unknown"


def discover_groups() -> list[dict]:
    """Walk each context dir, group .spe files by geometry."""
    groups = []
    for ctx in BG_CONTEXTS:
        d = PROJECT_REF_BASE / ctx["subdir"]
        if not d.is_dir():
            print(f"  ⚠ {d} does not exist, skipping context "
                  f"{ctx['key']!r}", file=sys.stderr)
            continue
        files = sorted(f for f in os.listdir(d) if f.endswith(".spe"))
        by_geo: dict[str, list[Path]] = defaultdict(list)
        for f in files:
            full = d / f
            try:
                s = read_spectrum(str(full))
            except Exception as e:
                print(f"  ⚠ failed to read {f}: {e}", file=sys.stderr)
                continue
            by_geo[s.geometry or "unknown"].append(full)

        for geo, paths in by_geo.items():
            if len(paths) < ctx["min_files_per_geometry"]:
                print(f"  · {ctx['key']} / {geo}: only {len(paths)} files, "
                      f"skipping (need ≥{ctx['min_files_per_geometry']})",
                      file=sys.stderr)
                continue
            groups.append({
                "context_key": ctx["key"],
                "context_description": ctx["description"],
                "pair_with_geometries": ctx.get("pair_with_geometries", []),
                "geometry": geo,
                "geometry_short": _short_geo_name(geo),
                "paths": paths,
                "n": len(paths),
                "out_filename":
                    f"bg_{ctx['key']}_{_short_geo_name(geo)}.spe",
            })
    return groups


# ============================================================================
# Main
# ============================================================================

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true",
                    help="Show what would be generated, write nothing.")
    args = ap.parse_args(argv)

    if not args.dry_run:
        OUT_DIR.mkdir(parents=True, exist_ok=True)

    groups = discover_groups()
    if not groups:
        print("  No qualifying groups found.", file=sys.stderr)
        return 1

    print("=" * 72)
    print(f"F-43/F-44 background-archive generator (dry-run={args.dry_run})")
    print(f"Archive: {PROJECT_REF_BASE}")
    print(f"Output : {OUT_DIR}")
    print("=" * 72)
    print()

    summary_rows = []
    manifest = []

    for g in groups:
        ctx = g["context_key"]
        geo = g["geometry"]
        n = g["n"]
        out_name = g["out_filename"]

        print(f"-- {ctx} / {geo} (N={n})")
        try:
            avg = average_lsrm_spectra(
                [str(p) for p in g["paths"]],
                sample_id=f"bg_{ctx}_{g['geometry_short']}",
                comment=(f"Background archive entry, context={ctx}, "
                         f"geometry={geo}, N={n} "
                         f"({g['context_description']})"),
            )
        except (CalibrationMismatchError, IdentityMismatchError) as e:
            print(f"   FAIL: {type(e).__name__}: {e}")
            continue

        prov = avg.extras["averaging_provenance"]
        mode = avg.extras["averaging_mode"]
        live_h = avg.live_time / 3600.0
        sigma_red = avg.extras["averaging_sigma_reduction"]
        sum_counts = int(sum(int(c) for c in avg.counts))
        print(f"   mode={mode}, live={live_h:.1f} h, σ-red={sigma_red:.2f}×, "
              f"sum_counts={sum_counts}")

        out_path = OUT_DIR / out_name
        if not args.dry_run:
            write_lsrm_spe(avg, str(out_path), type_label="Фон")
            sidecar = out_path.with_suffix(".provenance.json")
            full_prov = dict(prov)
            full_prov["pair_with_geometries"] = g["pair_with_geometries"]
            full_prov["context_description"] = g["context_description"]
            with open(sidecar, "w", encoding="utf-8") as f:
                json.dump(full_prov, f, ensure_ascii=False, indent=2,
                          default=str)
            print(f"   wrote {out_path.name}")
            print(f"   wrote {sidecar.name}")

        summary_rows.append((ctx, geo, n, mode, live_h, sigma_red,
                             sum_counts, out_name))
        manifest.append({
            "out_filename": out_name,
            "context_key": ctx,
            "geometry": geo,
            "pair_with_geometries": g["pair_with_geometries"],
            "n_inputs": n,
            "aggregation_mode": mode,
            "total_live_time_h": round(live_h, 2),
            "sigma_reduction": round(sigma_red, 3),
            "total_counts": sum_counts,
        })

    if not args.dry_run:
        manifest_path = OUT_DIR / "MANIFEST.json"
        with open(manifest_path, "w", encoding="utf-8") as f:
            json.dump({
                "module": "build_averaged_backgrounds.py",
                "version": "v1.7.22 (F-43 + F-44)",
                "items": manifest,
            }, f, ensure_ascii=False, indent=2)
        print(f"\nManifest: {manifest_path}")

    print()
    print("=" * 120)
    print(f"{'Context':<40} {'Geometry':<16} {'N':>3} "
          f"{'mode':<18} {'live (h)':>9} {'σ-red':>6} {'Σcounts':>11}")
    print("-" * 120)
    for row in summary_rows:
        ctx, geo, n, mode, live_h, sr, sc, _ = row
        print(f"{ctx:<40} {geo:<16} {n:>3} "
              f"{mode:<18} {live_h:>9.1f} {sr:>6.2f}× {sc:>11d}")
    print("=" * 120)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
