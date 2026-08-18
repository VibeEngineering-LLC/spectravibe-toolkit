"""
BUG-1 / 2026-06-02 — Spectrum.sample_mass_kg / sample_mass_uncertainty_kg
from LSRM SAMPLEMASS field.

Source spectrum: detectors/Gamma-1S/reference_spectra/archive/
                 Th232_420-7-17_Маринелли_0cm.spe (line 21):
    SAMPLEMASS=1600.0;16.0
    SAMPLEVOLUME=1000.0;10.0

LSRM convention (gamma/io/lsrm_spe.py docstring lines 35-40):
    SAMPLEMASS=value;uncertainty       (grams; grams)
    SAMPLEVOLUME=value;uncertainty     (millilitres; millilitres)

The reader divides grams by 1000 to populate `Spectrum.sample_mass_kg`
(typed); SAMPLEVOLUME flows through verbatim (millilitres).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

import pytest

from gamma.io.readers import read_spectrum
from gamma.io.lsrm_spe import _parse_value_err_pair_full


TH232_REF = (
    "detectors/Gamma-1S/reference_spectra/archive/"
    "Th232_420-7-17_Маринелли_0cm.spe"
)


def test_th232_sample_mass_kg_from_samplemass():
    """SAMPLEMASS=1600.0;16.0 (grams) → 1.6 kg ± 0.016 kg."""
    spec = read_spectrum(TH232_REF)
    assert spec.sample_mass_kg == pytest.approx(1.6, abs=1e-9), (
        f"Expected 1.6 kg from SAMPLEMASS=1600.0, got {spec.sample_mass_kg}"
    )
    assert spec.sample_mass_uncertainty_kg == pytest.approx(0.016, abs=1e-9), (
        f"Expected 0.016 kg uncertainty from SAMPLEMASS=...;16.0, got "
        f"{spec.sample_mass_uncertainty_kg}"
    )


def test_th232_sample_volume_ml_from_samplevolume():
    """SAMPLEVOLUME=1000.0;10.0 (mL) → 1000.0 mL ± 10.0 mL."""
    spec = read_spectrum(TH232_REF)
    assert spec.sample_volume_ml == pytest.approx(1000.0, abs=1e-9)
    assert spec.sample_volume_uncertainty_ml == pytest.approx(10.0, abs=1e-9)


def test_th232_extras_back_compat():
    """Legacy `extras["lsrm_sample_mass_kg"]` (F-140 contract) is preserved
    so that staged_pipeline's F-378 mismatch check and F-140 auto-extract
    continue to work without modification."""
    spec = read_spectrum(TH232_REF)
    assert spec.extras.get("lsrm_sample_mass_kg") == pytest.approx(1.6)
    assert spec.extras.get("lsrm_mass_source") == "sample_mass_field"


def test_parse_value_err_pair_full_basic():
    """Unit-test the helper on the canonical LSRM forms."""
    # Standard `value;uncertainty`
    assert _parse_value_err_pair_full("1600.0;16.0") == (1600.0, 16.0)
    # Value only (no semicolon) — uncertainty becomes None
    assert _parse_value_err_pair_full("1600.0") == (1600.0, None)
    # Empty / None
    assert _parse_value_err_pair_full("") is None
    assert _parse_value_err_pair_full(None) is None
    # Comma-decimal tolerated by _safe_float
    assert _parse_value_err_pair_full("1600,0;16,0") == (1600.0, 16.0)


def test_parse_value_err_pair_full_unparseable_uncertainty():
    """Unparseable uncertainty must not break value extraction."""
    # Garbage uncertainty → (value, None)
    assert _parse_value_err_pair_full("1600.0;xyz") == (1600.0, None)


def test_mass_field_absent_when_no_samplemass():
    """When the source format lacks SAMPLEMASS, `sample_mass_kg` is None.

    Cs-137 reference fixture is a calibration-source `.spe` that DOES carry
    SAMPLEMASS; for a negative case we synthesize a minimal Spectrum and
    check the default. (Reading every other format would belong in the
    respective reader test files.)
    """
    from gamma.spectrum import Spectrum
    import numpy as np
    spec = Spectrum(
        counts=np.zeros(10, dtype=np.int64),
        live_time=1.0,
        real_time=1.0,
    )
    assert spec.sample_mass_kg is None
    assert spec.sample_mass_uncertainty_kg is None
    assert spec.sample_volume_ml is None
    assert spec.sample_volume_uncertainty_ml is None


if __name__ == "__main__":
    test_th232_sample_mass_kg_from_samplemass()
    test_th232_sample_volume_ml_from_samplevolume()
    test_th232_extras_back_compat()
    test_parse_value_err_pair_full_basic()
    test_parse_value_err_pair_full_unparseable_uncertainty()
    test_mass_field_absent_when_no_samplemass()
    print("✓ All BUG-1 LSRM SAMPLEMASS tests passed.")
