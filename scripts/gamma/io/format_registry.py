"""
Format registry for the gamma-spectrum converter.

Maps file extensions and content sniffers to format identifiers, and
exposes the bidirectional reader/writer dispatch used by
`gamma.io.convert.convert_spectrum`.

A "format id" is a short stable string (e.g. `"lsrm_spe"`, `"n42_2012"`)
that identifies one of our supported formats. Adding a new format is a
single dict-update in `_FORMATS`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

from gamma.spectrum import Spectrum

from gamma.io.atomspectra_txt import (
    looks_like_atomspectra_txt,
    read_atomspectra_txt,
)
from gamma.io.atomspectra_xml import read_atomspectra_xml
from gamma.io.becqmoni_xml import write_becqmoni_xml
from gamma.io.lsrm_spe import read_lsrm_spe, write_lsrm_spe
from gamma.io.lsrm_spe_text import (
    read_lsrm_spe_text,
    write_lsrm_spe_text,
    looks_like_lsrm_spe_text,
)
from gamma.io.n42_2012 import read_n42_2012, write_n42_2012
from gamma.io.progress_spc import (
    looks_like_progress_spc,
    read_progress_spc,
)
from gamma.io.aspect_spc import (
    looks_like_aspect_spc,
    read_aspect_spc,
)


# ============================================================================
# Format records
# ============================================================================

@dataclass(frozen=True)
class FormatSpec:
    id: str
    label: str
    extensions: tuple
    reader: Optional[Callable[..., Spectrum]] = None
    writer: Optional[Callable[..., None]] = None
    sniffer: Optional[Callable[[bytes], bool]] = None


_FORMATS: dict[str, FormatSpec] = {}


def register(spec: FormatSpec) -> None:
    _FORMATS[spec.id] = spec


# ----- bootstrapping: builtin formats -----

register(FormatSpec(
    id="lsrm_spe",
    label="LSRM SpectraLine (binary .spe)",
    extensions=(".spe",),
    reader=read_lsrm_spe,
    writer=write_lsrm_spe,
    sniffer=lambda head: b"SPECTR=" in head[:4096] and not looks_like_lsrm_spe_text(head),
))

register(FormatSpec(
    id="lsrm_spe_text",
    label="LSRM SpectraLine (ASCII $-section .spe)",
    extensions=(".spe",),
    reader=read_lsrm_spe_text,
    writer=write_lsrm_spe_text,
    sniffer=looks_like_lsrm_spe_text,
))

register(FormatSpec(
    id="becqmoni_xml",
    label="BecqMoni / AtomSpectra ResultDataFile XML",
    extensions=(".xml",),
    reader=read_atomspectra_xml,
    writer=write_becqmoni_xml,
    sniffer=lambda head: b"<ResultDataFile" in head[:4096],
))

register(FormatSpec(
    id="n42_2012",
    label="ANSI/IEEE N42.42-2012 XML",
    extensions=(".n42", ".xml"),
    reader=read_n42_2012,
    writer=write_n42_2012,
    sniffer=lambda head: b"<RadInstrumentData" in head[:4096],
))

register(FormatSpec(
    id="atomspectra_txt",
    label="AtomSpectra mobile FORMAT: 3 (text)",
    extensions=(".txt",),
    reader=read_atomspectra_txt,
    writer=None,
    sniffer=looks_like_atomspectra_txt,
))

register(FormatSpec(
    id="progress_spc",
    label="Progress / Amplituda binary .spc (SF9I)",
    extensions=(".spc",),
    reader=read_progress_spc,
    writer=None,
    sniffer=looks_like_progress_spc,
))

register(FormatSpec(
    id="aspect_spc",
    label="ASPECT (Dubna) binary .spc (512-byte tSpectrHeader)",
    extensions=(".spc",),
    reader=read_aspect_spc,
    writer=None,
    sniffer=looks_like_aspect_spc,
))


# ============================================================================
# Public API
# ============================================================================

def list_formats() -> list[FormatSpec]:
    return list(_FORMATS.values())


def get_format(fmt_id: str) -> FormatSpec:
    if fmt_id not in _FORMATS:
        raise KeyError(
            f"Unknown format {fmt_id!r}. Known: {sorted(_FORMATS)}"
        )
    return _FORMATS[fmt_id]


def detect_format(path: str) -> str:
    """
    Identify the format of an existing file. Strategy:
      1. Sniff the first 8 KiB through all registered sniffers; if exactly
         one matches, use it. If multiple match, pick the one whose
         extension also matches (tie-break).
      2. Otherwise fall back to extension lookup.

    Returns the format id. Raises ValueError if no match.
    """
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"File not found: {path}")

    head = p.read_bytes()[:8192]
    matches = [
        spec for spec in _FORMATS.values()
        if spec.sniffer is not None and spec.sniffer(head)
    ]
    if len(matches) == 1:
        return matches[0].id
    if len(matches) > 1:
        ext = p.suffix.lower()
        for spec in matches:
            if ext in spec.extensions:
                return spec.id
        # No extension-confirmed tie-break — return the first match
        return matches[0].id

    # Fall back to extension only
    ext = p.suffix.lower()
    candidates = [s for s in _FORMATS.values() if ext in s.extensions]
    if len(candidates) == 1:
        return candidates[0].id
    if len(candidates) > 1:
        # Ambiguous extension with no sniff hit — bias toward the first
        # registered (LSRM .spe before IAEA .spe because most local files
        # are LSRM in this lab).
        return candidates[0].id

    raise ValueError(
        f"Cannot identify format for {path!r}. "
        f"Extension {ext!r} not in known formats and no sniffer matched."
    )


def read(path: str, fmt_id: Optional[str] = None, **kwargs) -> Spectrum:
    """Read a spectrum file via the registry."""
    fmt_id = fmt_id or detect_format(path)
    spec = get_format(fmt_id)
    if spec.reader is None:
        raise NotImplementedError(f"Format {fmt_id!r} has no reader")
    return spec.reader(path, **kwargs)


def write(spectrum: Spectrum, path: str, fmt_id: Optional[str] = None,
          **kwargs) -> None:
    """Write a spectrum via the registry. fmt_id required for writers."""
    if fmt_id is None:
        # Try to infer from the output extension
        ext = Path(path).suffix.lower()
        candidates = [
            s for s in _FORMATS.values()
            if ext in s.extensions and s.writer is not None
        ]
        if not candidates:
            raise ValueError(
                f"No writer registered for extension {ext!r}. "
                f"Pass fmt_id explicitly."
            )
        # Prefer formats that are unambiguously identified by ext.
        spec = candidates[0]
        fmt_id = spec.id
    else:
        spec = get_format(fmt_id)

    if spec.writer is None:
        raise NotImplementedError(f"Format {fmt_id!r} has no writer")
    spec.writer(spectrum, path, **kwargs)
