#!/usr/bin/env python3
from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""F-311 / v1.18.11 — Regenerate demo reports with v1.18.x opt-in integrations.

Запускает analyze_and_report на 4 Gamma-1S Marinelli .spe фикстурах с включёнными
TCS / Cutshall / matrix-method поправками. Output → demo_reports/v1_18_11/.

Usage
-----
    python scripts/regen_demo_reports.py [--output-dir demo_reports/v1_18_11]

Output artefacts per fixture:
- {stem}_report.json
- {stem}_report.md
- {stem}_spectrum.png  (если matplotlib доступен)
- {stem}.html
"""

import argparse
import json
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

FIXTURES_DIR = REPO / "evals" / "fixtures"

FIXTURES = [
    ("Cs-137",  FIXTURES_DIR / "M_cs_легкий_2001-2005.spe", 0.570),
    ("K-40",    FIXTURES_DIR / "M_k_легкий_2001-2005.spe",  0.665),
    ("Ra-226",  FIXTURES_DIR / "M_ra_легкий_2001-2007.spe", 0.622),
    ("Th-232",  FIXTURES_DIR / "M_th_легкий_2001-2005.spe", 0.550),
]

# F-313 / v1.18.12 — explicit background path для archive фикстур.
# F-131 эвристика максимум 90 days_apart по умолчанию; archive sample 2001-2005
# vs bg 2016 = ~11 лет → reject. Передаём фон явно (это разрешено F-135
# контрактом: explicit background_path всегда применяется).
DEFAULT_BG_FOR_MARINELLI = (
    REPO / "detectors" / "Gamma-1S" / "data" / "averaged_backgrounds"
    / "bg_2016_marinelli_water_marinelli.spe"
)


def run_demo(spe_path: Path, sample_mass_kg: float,
             output_dir: Path, *, with_corrections: bool = True,
             ) -> dict:
    """Run full-report pipeline on one fixture."""
    from gamma.reporting import analyze_and_report

    kwargs = dict(
        sample_mass_kg=sample_mass_kg,
        write_json=True,
        write_markdown=True,
        write_plots=False,    # требует matplotlib + большой размер
        write_html=True,
        allow_stage2=True,
        complete_workflow=True,
        # F-313 / v1.18.12 (P0): explicit bg path для archive фикстур.
        # Auto-search F-131 для них не работает (date filter).
        background_path=str(DEFAULT_BG_FOR_MARINELLI),
    )
    if with_corrections:
        # F-309 / v1.18.8 opt-in flags на all 3 corrections
        kwargs.update(
            enable_tcs_correction=True,
            tcs_detector_id="Gamma-1S",
            # Для Marinelli legacy легкий ρ~1.0-1.2 → Cutshall ок
            sample_density_g_cm3=1.2,
            enable_cutshall_self_abs=True,
            cutshall_calib_density_g_cm3=1.0,
            # Matrix method для Ra-226 / Th-232 multi-nuclide
            enable_matrix_method=True,
            matrix_method_energy_tolerance_keV=2.0,
            # F-322 / v1.18.16 — F-96 bg-anchors для multiplet deconv
            # (закрывает M_th M3 chi2_red=12.94 — добавляет 511 annihilation
            # как constraint в clusters содержащие линии в ±60 keV окне 511).
            enable_f96_bg_anchors=True,
        )

    return analyze_and_report(
        str(spe_path),
        output_dir=str(output_dir),
        **kwargs,
    )


def summarize(art: dict) -> dict:
    """Extract key numbers for diff."""
    res = art.get("result")
    if res is None:
        return {"error": "no result"}
    activities = getattr(res, "activities", None) or []
    return {
        "n_nuclides": len(activities),
        "activities_Bq": {
            getattr(a, "nuclide", "?"): float(getattr(a, "A_Bq", 0.0))
            for a in activities
        },
        "geometry": getattr(res, "geometry_canonical", ""),
        "detector": getattr(res, "detector_canonical", ""),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", default=str(REPO / "demo_reports" / "v1_18_11"),
        help="Where to write the demo artefacts",
    )
    parser.add_argument(
        "--no-corrections", action="store_true",
        help="Build baseline reports without v1.18.x corrections (для diff)",
    )
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {out}")
    print(f"With v1.18.x corrections: {not args.no_corrections}")
    print()

    summary_all = {}
    for nuc, path, mass in FIXTURES:
        print(f"=== {nuc} ({path.name}, {mass} kg) ===")
        try:
            art = run_demo(
                path, mass, out,
                with_corrections=not args.no_corrections,
            )
            summary = summarize(art)
            summary_all[nuc] = summary
            print(f"  detected {summary['n_nuclides']} nuclides, "
                  f"geometry={summary.get('geometry', '?')}")
            for n, A in summary.get("activities_Bq", {}).items():
                print(f"    {n:10s}  {A:.3e} Bq")
        except Exception as e:
            summary_all[nuc] = {"error": f"{type(e).__name__}: {e}"}
            print(f"  ERROR: {e}", file=sys.stderr)
        print()

    # Write consolidated summary
    summary_file = out / "DEMO_SUMMARY.json"
    summary_file.write_text(
        json.dumps(summary_all, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"Summary: {summary_file}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
