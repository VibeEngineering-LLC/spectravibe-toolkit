"""
gamma — gamma-ray spectrum analysis toolkit.

Public surface (Phase 0):
    gamma.spectrum.Spectrum            — the parsed-spectrum dataclass
    gamma.spectrum.ENERGY_CEILING_KEV  — 3000 keV upper-energy bound
    gamma.io.readers.read_spectrum     — dispatch reader by extension
    gamma.io.atomspectra_xml.read_atomspectra_xml — AtomSpectra reader
    gamma.io.background.resolve_external_background
    gamma.cli                          — `python -m gamma.cli`

More modules (calibration, peaks, identification, physics, reporting)
will appear in Phases 1–3.
"""

from gamma.spectrum import Spectrum, ENERGY_CEILING_KEV  # noqa: F401

__version__ = "1.27.9"
