"""
BUG-48 — per-detector FWHM reference lookup.

Tests for ``_fwhm_reference_keV(detector_class)`` in
``scripts/gamma/reporting/spectrum_qc_aggregator.py``.

Previously ``_FWHM_REFERENCE_KEV = 8.5`` was the single hard-coded default,
which caused criterion 2 to always FAIL for NaI spectra
(rel_dev = |47 - 8.5| / 8.5 ≈ 4.52, far above the 0.15 threshold).

Fix (BUG-48): replaced with ``_FWHM_REFERENCE_BY_DETECTOR`` dict +
``_fwhm_reference_keV(detector_class)`` helper with ``"default"`` fallback
of 47.0 keV (NaI-grade).

Sources:
    NaI 47.0 keV — staged_pipeline.py:485 ``_DEFAULT_NAI_FWHM_MODEL``
        → ``fwhm_keV_at_energy(model, 661.66) ≈ 46.95 keV``; LSRM-9.4 §3.2;
        RAG-043 (Gilmore & Joss §6.4).
    HPGe 1.5 keV — RAG-043 Gilmore & Joss §6.4 (~0.23% at 661 keV).
    LaBr 20.0 keV — manufacturer spec ~3.0% at 661 keV; RAG-043.
"""
from gamma.reporting.spectrum_qc_aggregator import _fwhm_reference_keV


def test_fwhm_reference_nai():
    """NaI reference must be 47.0 keV (7.1% resolution at 661 keV)."""
    assert _fwhm_reference_keV("NaI") == 47.0


def test_fwhm_reference_hpge():
    """HPGe reference must be 1.5 keV (~0.23% resolution at 661 keV)."""
    assert _fwhm_reference_keV("HPGe") == 1.5


def test_fwhm_reference_labr():
    """LaBr reference must be 20.0 keV (~3% resolution at 661 keV)."""
    assert _fwhm_reference_keV("LaBr") == 20.0


def test_fwhm_reference_unknown_uses_default():
    """Unknown and empty detector class must fall back to default (47.0 keV, NaI-grade)."""
    assert _fwhm_reference_keV("unknown") == 47.0
    assert _fwhm_reference_keV("") == 47.0
    assert _fwhm_reference_keV("CsI") == 47.0
