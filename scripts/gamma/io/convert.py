"""
High-level format conversion entry point.

`convert_spectrum(in_path, out_path)` autodetects both input and output
formats via the registry, reads the input into the canonical Spectrum
intermediate, and writes it out in the target format. Lossy conversions
report dropped fields when `verbose=True`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from gamma.spectrum import Spectrum
from gamma.io import format_registry as fr


def convert_spectrum(
    in_path: str,
    out_path: str,
    *,
    in_format: Optional[str] = None,
    out_format: Optional[str] = None,
    apply_energy_ceiling: bool = False,
    ceiling_keV: Optional[float] = None,
    verbose: bool = False,
) -> Spectrum:
    """
    Convert a spectrum file between any two supported formats.

    Args:
        in_path: source file path.
        out_path: destination file path.
        in_format: registry id of the input format; auto-detected if None.
        out_format: registry id of the output format; inferred from
            `out_path` extension if None.
        apply_energy_ceiling: passed to the input reader. Default False
            for conversion (we want the full channel range preserved).
        ceiling_keV: per-call override of `ENERGY_CEILING_KEV`. Used only
            when `apply_energy_ceiling=True`.
        verbose: print a one-line summary of dropped fields, if any.

    Returns:
        The intermediate Spectrum (for chaining or inspection).
    """
    in_p = Path(in_path)
    out_p = Path(out_path)

    in_fmt = in_format or fr.detect_format(str(in_p))
    out_fmt = out_format
    if out_fmt is None:
        # Infer from out extension
        ext = out_p.suffix.lower()
        matches = [
            s for s in fr.list_formats()
            if ext in s.extensions and s.writer is not None
        ]
        if not matches:
            raise ValueError(
                f"Cannot infer output format from extension {ext!r}. "
                f"Pass out_format explicitly."
            )
        # When extension is ambiguous (.spe → lsrm vs iaea, .xml → becqmoni
        # vs n42), prefer the same format as input if applicable, else
        # take the first registered.
        if any(m.id == in_fmt for m in matches):
            out_fmt = in_fmt
        else:
            out_fmt = matches[0].id

    # ---- read ----
    read_kwargs = {
        "apply_energy_ceiling": apply_energy_ceiling,
    }
    if apply_energy_ceiling and ceiling_keV is not None:
        read_kwargs["ceiling_keV"] = ceiling_keV
    spec = fr.read(str(in_p), fmt_id=in_fmt, **read_kwargs)

    # ---- write ----
    fr.write(spec, str(out_p), fmt_id=out_fmt)

    if verbose:
        dropped = _summarise_lossy_fields(spec, in_fmt, out_fmt)
        if dropped:
            print(f"[convert] {in_fmt} -> {out_fmt}: dropped {', '.join(dropped)}")
        else:
            print(f"[convert] {in_fmt} -> {out_fmt}: clean round-trip "
                  f"(n_channels={spec.n_channels}, "
                  f"live={spec.live_time:.2f}s, real={spec.real_time:.2f}s)")

    return spec


def _summarise_lossy_fields(
    spec: Spectrum, in_fmt: str, out_fmt: str
) -> list[str]:
    """List Spectrum fields that will not survive the chosen out format."""
    losses: list[str] = []
    if out_fmt == "lsrm_spe_text":
        # LSRM ASCII SPE has no SampleInfo/DeviceConfigReference/full FWHM
        # block beyond $SHAPE_CAL; richer metadata is dropped.
        if spec.detector_id and in_fmt != "lsrm_spe_text":
            losses.append("device_metadata")
        if spec.extras:
            losses.append(f"extras({len(spec.extras)})")
    if out_fmt == "n42_2012":
        if spec.extras.get("lsrm_peaks_table") and in_fmt == "lsrm_spe":
            losses.append("lsrm_peaks_table")
        if spec.stored_fwhm_calibration and in_fmt != "n42_2012":
            losses.append("fwhm_calibration")
    return losses
