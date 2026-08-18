#!/usr/bin/env python3
from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""
CLI for the gamma-spectrum file format converter.

Usage:
    python convert_spectrum.py IN.spe OUT.n42
    python convert_spectrum.py IN.xml OUT.spe --out-format lsrm_spe
    python convert_spectrum.py --list-formats

Format ids (use with --in-format / --out-format):
    lsrm_spe        LSRM SpectraLine binary .spe (CP-1251 + uint32-LE)
    lsrm_spe_text   LSRM SpectraLine ASCII .spe ($-section export)
    becqmoni_xml    BecqMoni / AtomSpectra ResultDataFile (.xml)
    n42_2012        ANSI/IEEE N42.42-2012 (.n42, .xml)

The format of both input and output is autodetected from the file's
content (input) or extension (output). Pass --in-format / --out-format
to override.
"""


import argparse
import sys
from pathlib import Path

# Ensure scripts/ is importable when run directly
_THIS = Path(__file__).resolve().parent
if str(_THIS) not in sys.path:
    sys.path.insert(0, str(_THIS))

from gamma.io import format_registry as fr
from gamma.io.convert import convert_spectrum


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="convert_spectrum",
        description="Convert gamma spectra between LSRM .spe / BecqMoni XML "
                    "/ ANSI N42-2012 / IAEA SPE formats.",
    )
    parser.add_argument("input", nargs="?", help="input spectrum file")
    parser.add_argument("output", nargs="?", help="output spectrum file")
    parser.add_argument(
        "--in-format", default=None,
        help="override input format detection (one of: "
             "lsrm_spe, iaea_spe, becqmoni_xml, n42_2012)",
    )
    parser.add_argument(
        "--out-format", default=None,
        help="override output format selection",
    )
    parser.add_argument(
        "--apply-energy-ceiling", action="store_true",
        help="apply the project's 3 MeV energy ceiling on read",
    )
    parser.add_argument(
        "--ceiling-keV", type=float, default=None,
        help="per-call ceiling override (keV); implies --apply-energy-ceiling",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="print conversion summary",
    )
    parser.add_argument(
        "--list-formats", action="store_true",
        help="print supported formats and exit",
    )
    args = parser.parse_args(argv)

    if args.list_formats:
        _print_format_table()
        return 0

    if not args.input or not args.output:
        parser.error("INPUT and OUTPUT paths are required "
                     "(or pass --list-formats).")
        return 2

    apply_ceil = args.apply_energy_ceiling or (args.ceiling_keV is not None)
    convert_spectrum(
        args.input,
        args.output,
        in_format=args.in_format,
        out_format=args.out_format,
        apply_energy_ceiling=apply_ceil,
        ceiling_keV=args.ceiling_keV,
        verbose=args.verbose,
    )
    return 0


def _print_format_table() -> None:
    print("ID            R W  EXT            LABEL")
    print("-" * 70)
    for spec in fr.list_formats():
        r = "R" if spec.reader else "-"
        w = "W" if spec.writer else "-"
        exts = ",".join(spec.extensions)
        print(f"{spec.id:<13} {r} {w}  {exts:<14} {spec.label}")


if __name__ == "__main__":
    raise SystemExit(main())
