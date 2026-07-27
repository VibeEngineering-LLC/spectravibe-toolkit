# -*- coding: utf-8 -*-
"""Convert the Gamma-1S reference kits (LSRM .spe) to BecqMoni XML.

Each kit leaf directory holds a matched pair — one sample spectrum and the
background measured in the same geometry. BecqMoni carries the background
inside the same file as `<BackgroundEnergySpectrum>`, so one XML is emitted
per pair rather than one per .spe.

Layout in  -> detectors/<detector>/reference_spectra/reference_kits/<geometry>/<nuclide>/*.spe
Layout out -> detectors/<detector>/reference_spectra/reference_kits_becqmoni/<geometry>/<nuclide>/<sample>.xml

Usage:
    python scripts/convert/reference_kits_to_becqmoni.py
    python scripts/convert/reference_kits_to_becqmoni.py --dry-run
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "scripts" / "convert"))

from device_guid import device_guid  # noqa: E402

from gamma.io.lsrm_spe import read_lsrm_spe          # noqa: E402
from gamma.io.becqmoni_xml import write_becqmoni_xml  # noqa: E402

DEFAULT_KITS = REPO / "detectors" / "Gamma-1S" / "reference_spectra" / "reference_kits"
DEFAULT_OUT = REPO / "detectors" / "Gamma-1S" / "reference_spectra" / "reference_kits_becqmoni"

# Background files are named either `background_*.spe` (canonical, written by
# build_reference_kits.py) or carry the raw LSRM operator name «Фон …» for the
# pairs that were pinned by hand.
BG_PREFIXES = ("background_", "bg_")
BG_TOKENS = ("фон",)


# BecqMoni stores SampleInfo.Weight in kilograms and SampleInfo.Volume in
# litres, then scales for display: `DCSampleInfoView.LoadFormContents()`
# assigns `Weight * 1000` / `Volume * 1000` to a NumericUpDown whose range is
# 0.001…100000 in gram/millilitre mode and 0.001…100 in kg/litre mode.
# LSRM .spe carries grams and millilitres, so passing the raw numbers through
# blows the control up: a 120 ml Denta vial becomes 120 l -> 120000 ml ->
# ArgumentOutOfRangeException before the document is drawn. A point source
# recorded as 0 fails the same check from the other side (< Minimum).
BQ_MIN = 0.001
BQ_MAX = 100.0
BQ_NEUTRAL = 1.0  # what BecqMoni itself writes when the field is unknown


def to_becqmoni_unit(raw) -> str:
    """LSRM grams/millilitres -> BecqMoni kilograms/litres, range-checked."""
    try:
        value = float(str(raw).split(";")[0].replace(",", "."))
    except (TypeError, ValueError):
        return f"{BQ_NEUTRAL:g}"
    if value <= 0:
        return f"{BQ_NEUTRAL:g}"
    value = min(max(value / 1000.0, BQ_MIN), BQ_MAX)
    return f"{value:g}"


def is_background(path: Path) -> bool:
    name = path.name.lower()
    return name.startswith(BG_PREFIXES) or any(tok in name for tok in BG_TOKENS)


def find_pairs(kits_dir: Path):
    """Yield (leaf_dir, sample_path, background_path) for every kit entry."""
    for leaf in sorted(p for p in kits_dir.rglob("*") if p.is_dir()):
        spes = sorted(leaf.glob("*.spe"))
        if not spes:
            continue
        bgs = [p for p in spes if is_background(p)]
        samples = [p for p in spes if not is_background(p)]
        if len(samples) != 1 or len(bgs) != 1:
            yield leaf, samples, bgs  # reported as a problem by the caller
            continue
        yield leaf, samples[0], bgs[0]


def convert(kits_dir: Path, out_dir: Path, dry_run: bool = False) -> int:
    ok, problems = 0, []
    for leaf, sample, bg in find_pairs(kits_dir):
        rel = leaf.relative_to(kits_dir)
        if isinstance(sample, list):
            problems.append(
                f"{rel}: expected 1 sample + 1 background, got "
                f"{len(sample)} sample(s) / {len(bg)} background(s)")
            continue

        spec = read_lsrm_spe(str(sample))
        spec.background_embedded = read_lsrm_spe(str(bg))
        spec.background_link = bg.name
        if not spec.sample_id:
            spec.sample_id = sample.stem
        if not spec.detector_id:
            spec.detector_id = "Gamma-1S"
        # The writer falls back to uuid4() when device_guid is empty, so every
        # run would emit fresh GUIDs — BecqMoni would register each rebuild as
        # a new device, and git would see all 40 files change. device_guid()
        # keeps it stable per detector; the salt it uses lives outside the
        # repository, otherwise the GUID would give the detector name away.
        if not spec.device_guid:
            spec.device_guid = device_guid(spec.detector_id)
        # The writer copies these straight into <Weight>/<Volume>; normalise
        # the units here so BecqMoni can open the result (see to_becqmoni_unit).
        spec.extras["lsrm_samplemass"] = to_becqmoni_unit(
            spec.extras.get("lsrm_samplemass", ""))
        spec.extras["lsrm_samplevolume"] = to_becqmoni_unit(
            spec.extras.get("lsrm_samplevolume", ""))

        dst = out_dir / rel / (sample.stem + ".xml")
        print(f"{rel}/{sample.name}"
              f"  ->  {dst.relative_to(out_dir.parent)}"
              f"   [{spec.n_channels or len(spec.counts)} ch, "
              f"live {spec.live_time:.0f} s]")
        if not dry_run:
            dst.parent.mkdir(parents=True, exist_ok=True)
            write_becqmoni_xml(spec, str(dst))
        ok += 1

    print(f"\nconverted: {ok}")
    if problems:
        print(f"problems: {len(problems)}")
        for p in problems:
            print("  ! " + p)
        return 1
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--kits-dir", type=Path, default=DEFAULT_KITS)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    if not args.kits_dir.is_dir():
        print(f"kits dir not found: {args.kits_dir}", file=sys.stderr)
        return 2
    return convert(args.kits_dir, args.out_dir, args.dry_run)


if __name__ == "__main__":
    raise SystemExit(main())
