"""
Unified spectrum reader.

read_spectrum(path) dispatches to the right format reader based on the
file extension AND optional content sniff (see `gamma.io.format_registry`).

Supported formats (read+write unless noted):
  .spe  → LSRM SpectraLine (binary header + uint32-LE counts)
  .spe  → IAEA SPE (ASCII $-section format)  ← disambiguated by content
  .xml  → BecqMoni / AtomSpectra ResultDataFile
  .xml  → ANSI N42-42-2012 RadInstrumentData ← disambiguated by content
  .n42  → ANSI N42-42-2012

Add new formats by registering them in `gamma.io.format_registry`.
"""

from __future__ import annotations

from pathlib import Path

from gamma.spectrum import Spectrum
from gamma.io import format_registry as _fr


def read_spectrum(path: str, **kwargs) -> Spectrum:
    """
    Dispatch to a format-specific reader.

    The format is selected by sniffing file content first, then falling
    back to extension. All keyword arguments are forwarded to the
    underlying reader. All readers accept the energy-ceiling contract:

      - `apply_energy_ceiling: bool = False` — drop channels whose
        calibrated energy exceeds the ceiling. **Default flipped from
        True to False in v1.18.32 (BUG-9, 2026-06-03)** to stop silently
        dropping the high-energy tail of low-channel-count NaI files
        (e.g. a 1024-ch Gamma-1S spectrum with a0≈-8 keV lost 21 trailing
        channels under the old default). Use `True` explicitly when the
        3 MeV trim is desired, or apply
        `gamma.spectrum.trim_to_working_energy` after reading.
      - `ceiling_keV: float | None = None` — per-call override of the
        `ENERGY_CEILING_KEV` constant (3000 keV by project scope).

    `read_atomspectra_xml` additionally accepts `parse_background: bool = True`.
    """
    return _fr.read(str(path), **kwargs)
