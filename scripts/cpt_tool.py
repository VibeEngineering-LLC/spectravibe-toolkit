#!/usr/bin/env python3
from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""F-301 / v1.18.3 — LSRM .cpt CLI utility.

CLI обёртка вокруг `gamma.io.cpt_io` (F-301) + `gamma.peaks.peak_image_tabulated`
(F-299) для работы с LSRM Calibrated Peak Template файлами.

Usage
-----
Чтение .cpt в JSON:
    python scripts/cpt_tool.py read input.cpt [-o output.json]

Создание .cpt из calibration pairs:
    python scripts/cpt_tool.py build \\
        --detector-id Gamma-1S --detector-class NaI --diameter-mm 63 \\
        --anchor 122.0:12.5 --anchor 662.0:46.5 --anchor 1332.0:78.0 \\
        --out template.cpt

Inspect (dump-info без записи):
    python scripts/cpt_tool.py inspect input.cpt

References
----------
- F-301 (v1.17.21) — XML I/O
- F-299 (v1.17.21) — tabulated peak image
- LSRM SpectraLine .cpt format
"""

import argparse
import json
import os
import sys
from pathlib import Path

# Make repo `scripts/` available on path when run as standalone script.
HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))


def _cmd_read(args: argparse.Namespace) -> int:
    from gamma.io.cpt_io import read_cpt_file
    src = Path(args.input)
    if not src.exists():
        print(f"ERROR: file not found: {src}", file=sys.stderr)
        return 2
    img = read_cpt_file(src, strict=args.strict)
    payload = {
        "detector_id": img.detector_id,
        "detector_class": img.detector_class,
        "crystal_diameter_mm": img.crystal_diameter_mm,
        "source_metadata": img.source_metadata,
        "notes": img.notes,
        "anchors": [
            {
                "E_keV": a.E_keV,
                "fwhm_keV": a.fwhm_keV,
                "tail_fraction": a.tail_fraction,
                "tail_slope_inv_keV": a.tail_slope_inv_keV,
                "step_height_frac": a.step_height_frac,
                "asymmetry": a.asymmetry,
                "weight": a.weight,
            }
            for a in img.anchors
        ],
    }
    text = json.dumps(payload, indent=2, ensure_ascii=False)
    if args.output:
        Path(args.output).write_text(text, encoding="utf-8")
        print(f"OK: wrote {args.output} ({len(img.anchors)} anchors)")
    else:
        print(text)
    return 0


def _parse_anchor_spec(spec: str) -> tuple:
    """Parse 'E_keV:fwhm_keV' anchor spec."""
    parts = spec.split(":")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            f"--anchor expects 'E_keV:fwhm_keV', got {spec!r}"
        )
    try:
        E = float(parts[0])
        f = float(parts[1])
    except ValueError as e:
        raise argparse.ArgumentTypeError(
            f"--anchor values not numeric: {spec!r}: {e}"
        )
    if E <= 0 or f <= 0:
        raise argparse.ArgumentTypeError(
            f"--anchor E and FWHM must be > 0, got {spec!r}"
        )
    return (E, f)


def _cmd_build(args: argparse.Namespace) -> int:
    from gamma.peaks.peak_image_tabulated import (
        build_anchors_from_calibration,
    )
    from gamma.io.cpt_io import write_cpt_file
    if not args.anchor:
        print(
            "ERROR: at least one --anchor E_keV:fwhm_keV required",
            file=sys.stderr,
        )
        return 2
    img = build_anchors_from_calibration(
        detector_id=args.detector_id,
        detector_class=args.detector_class,
        crystal_diameter_mm=args.diameter_mm,
        calibration_pairs=args.anchor,
    )
    issues = img.validate()
    if issues:
        print("WARNINGS:", file=sys.stderr)
        for i in issues:
            print(f"  - {i}", file=sys.stderr)
    write_cpt_file(img, args.out)
    print(
        f"OK: wrote {args.out} "
        f"({len(img.anchors)} anchors, "
        f"detector={img.detector_id}/{img.detector_class})"
    )
    return 0


def _cmd_inspect(args: argparse.Namespace) -> int:
    from gamma.io.cpt_io import read_cpt_file
    from gamma.peaks.peak_image_logspline import interpolate_peak_shape
    src = Path(args.input)
    if not src.exists():
        print(f"ERROR: file not found: {src}", file=sys.stderr)
        return 2
    img = read_cpt_file(src, strict=False)
    print(f"=== .cpt: {src.name} ===")
    print(f"  detector_id      : {img.detector_id}")
    print(f"  detector_class   : {img.detector_class}")
    print(f"  crystal_diameter : {img.crystal_diameter_mm:.2f} mm")
    print(f"  source_metadata  : {img.source_metadata}")
    print(f"  notes            : {img.notes}")
    print(f"  anchors          : {len(img.anchors)}")
    for a in img.anchors:
        print(
            f"    E={a.E_keV:8.2f}  FWHM={a.fwhm_keV:6.2f} keV  "
            f"tail={a.tail_fraction:.4f}  step={a.step_height_frac:.4f}"
        )
    issues = img.validate()
    if issues:
        print()
        print("VALIDATION WARNINGS:")
        for i in issues:
            print(f"  - {i}")
    # Estimate FWHM @ 662 keV (NaI cal reference)
    try:
        fwhm_at_662 = interpolate_peak_shape(img, 662.0)
        pct = (fwhm_at_662.fwhm_keV / 662.0) * 100.0
        print(
            f"  FWHM%@662keV     : {pct:.2f}% "
            f"({'extrap' if fwhm_at_662.was_extrapolated else 'interp'})"
        )
    except Exception as e:
        print(f"  FWHM%@662keV     : N/A ({e})")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="cpt_tool",
        description="LSRM .cpt (Calibrated Peak Template) CLI utility.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_read = sub.add_parser("read", help="Read .cpt → JSON")
    p_read.add_argument("input", help="Input .cpt file")
    p_read.add_argument("-o", "--output", help="Output .json (default stdout)")
    p_read.add_argument(
        "--strict", action="store_true",
        help="Strict mode: raise on unknown tags",
    )
    p_read.set_defaults(func=_cmd_read)

    p_build = sub.add_parser(
        "build", help="Build .cpt from calibration anchors",
    )
    p_build.add_argument(
        "--detector-id", required=True,
        help="Detector identifier (e.g. 'Gamma-1S')",
    )
    p_build.add_argument(
        "--detector-class", default="NaI",
        choices=["NaI", "HPGe", "CsI", "LaBr", "CeBr"],
        help="Detector type (default NaI)",
    )
    p_build.add_argument(
        "--diameter-mm", type=float, default=63.0,
        help="Crystal diameter in mm (default 63)",
    )
    p_build.add_argument(
        "--anchor", action="append", type=_parse_anchor_spec,
        metavar="E_keV:fwhm_keV", default=[],
        help="Add calibration anchor (can be repeated)",
    )
    p_build.add_argument(
        "--out", required=True,
        help="Output .cpt file path",
    )
    p_build.set_defaults(func=_cmd_build)

    p_inspect = sub.add_parser(
        "inspect", help="Inspect .cpt content (human-readable dump)",
    )
    p_inspect.add_argument("input", help="Input .cpt file")
    p_inspect.set_defaults(func=_cmd_inspect)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
