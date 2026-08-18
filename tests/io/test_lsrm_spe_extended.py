"""
Tests for extended Lsrm SpectraLine SPE reader features per official
format specification (PDF supplied by user):
  - ENERGY_QUALITY (chi2, integral_nonlinearity, calibration_nonlinearity)
  - PEAKS table parsing (per Lsrm spec §7.5.2.1)
  - ZONES table parsing (per Lsrm spec §7.5.2.1)
  - FWHM_ORT JSON detection
  - Orthogonal polynomial fallback flag
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.io.readers import read_spectrum
from gamma.io.lsrm_spe import _parse_peaks_table, _parse_zones_table


# Cs-137 reference file is known to have all extended fields
CS137_PATH = ("detectors/Gamma-1S/reference_spectra/archive/"
              "Cs-137__163_2017.spe")


def test_lsrm_spe_basic_calibration():
    """Verify the polynomial energy calibration is read correctly."""
    spec = read_spectrum(CS137_PATH)
    assert spec.energy_cal is not None
    assert len(spec.energy_cal) >= 2
    # ENERGY=4,... → 5 coefficients (degree-4 polynomial)
    # First coefficient should be near -10.05
    assert -11 < spec.energy_cal[0] < -9, f"a0={spec.energy_cal[0]}"
    # Should give ~661 keV at the Cs-137 photopeak (around ch=232)
    E = spec.channel_to_energy(232)
    assert 655 < E < 665, f"Cs-137 photopeak at ch=232 gave E={E:.2f}"
    print(f"  ✓ test_lsrm_spe_basic_calibration (E@ch=232 = {E:.2f} keV)")


def test_lsrm_spe_energy_quality_parsed():
    """ENERGY_QUALITY field should be parsed into a 3-value list."""
    spec = read_spectrum(CS137_PATH)
    qual = spec.extras.get("lsrm_energy_quality")
    assert qual is not None, "ENERGY_QUALITY not parsed"
    assert len(qual) == 3, f"Expected 3 values, got {len(qual)}"
    # First value is chi2 (typically >10)
    assert qual[0] > 0
    print(f"  ✓ test_lsrm_spe_energy_quality_parsed: chi2={qual[0]:.2f}, "
          f"int_nonlin={qual[1]:.4f}, cal_nonlin={qual[2]:.4f}")


def test_lsrm_spe_peaks_table_parsed():
    """PEAKS table should be parsed into structured rows."""
    spec = read_spectrum(CS137_PATH)
    peaks = spec.extras.get("lsrm_peaks_table")
    assert peaks is not None, "PEAKS table not parsed"
    assert len(peaks) >= 1
    p0 = peaks[0]
    # Required fields
    assert "position_ch" in p0
    assert "energy_keV" in p0
    assert "fwhm_keV" in p0
    assert "area" in p0
    # Cs-137 661 keV should be in there
    cs_peak = next((p for p in peaks if 655 < p["energy_keV"] < 670), None)
    assert cs_peak is not None, "Cs-137 661 keV not in PEAKS table"
    print(f"  ✓ test_lsrm_spe_peaks_table_parsed "
          f"({len(peaks)} peaks; Cs-137 at E={cs_peak['energy_keV']:.2f})")


def test_lsrm_spe_zones_table_parsed():
    """ZONES table should be parsed into structured rows."""
    spec = read_spectrum(CS137_PATH)
    zones = spec.extras.get("lsrm_zones_table")
    assert zones is not None, "ZONES table not parsed"
    assert len(zones) >= 1
    z0 = zones[0]
    assert "left_bound" in z0
    assert "right_bound" in z0
    assert "n_peaks_in_zone" in z0
    assert "minimize" in z0
    # Bounds should be sensible (in keV — typically 100-3000 range)
    assert 0 < z0["left_bound"] < z0["right_bound"]
    print(f"  ✓ test_lsrm_spe_zones_table_parsed ({len(zones)} zones, "
          f"first=[{z0['left_bound']:.1f}, {z0['right_bound']:.1f}], "
          f"n_peaks={z0['n_peaks_in_zone']})")


def test_parse_peaks_table_isolated():
    """Direct unit test of _parse_peaks_table on synthetic input."""
    # From Lsrm spec §7.5.2.1 example:
    table = ("114.898 0.039 39.680 0.014 0.868 0.028 19773 140 17.455 0 0\n"
             "130.071 0.071 45.154 0.026 0.887 0.049 6600 96 17.455 0 0\n"
             "133.620 0.170 46.435 0.061 0.891 0.094 1787 55 17.455 0 0")
    parsed = _parse_peaks_table(table)
    assert len(parsed) == 3
    assert parsed[0]["position_ch"] == 114.898
    assert parsed[0]["energy_keV"] == 39.680
    assert parsed[0]["area"] == 19773
    assert parsed[0]["d_area"] == 140
    print(f"  ✓ test_parse_peaks_table_isolated (3 peaks parsed)")


def test_parse_zones_table_isolated():
    """Direct unit test of _parse_zones_table."""
    # From Lsrm spec §7.5.2.1 example
    table = ("33.944 52.345 1 FWHM, Position, Step, Linear -1\n"
             "116.568 127.032 1 FWHM, Position, Step, Linear -1\n"
             "238.526 250.795 1 FWHM, Position, Step, Linear -1")
    parsed = _parse_zones_table(table)
    assert len(parsed) == 3
    assert parsed[0]["left_bound"] == 33.944
    assert parsed[0]["right_bound"] == 52.345
    assert parsed[0]["n_peaks_in_zone"] == 1
    assert parsed[0]["bg_polynomial_degree"] == -1
    assert "FWHM" in parsed[0]["minimize"]
    print(f"  ✓ test_parse_zones_table_isolated (3 zones parsed)")


def test_lsrm_spe_metadata_preserved():
    """All metadata fields should be carried through to spec.extras."""
    spec = read_spectrum(CS137_PATH)
    # SHIFR and other identity fields
    assert spec.sample_id != "", "SHIFR not parsed"
    assert spec.detector_id != "", "DETECTOR not parsed"
    assert spec.geometry != "", "GEOMETRY not parsed"
    # Distance for Точечная-5см should be 5
    distance_str = spec.extras.get("lsrm_distance")
    assert distance_str is not None, "DISTANCE not preserved"
    print(f"  ✓ test_lsrm_spe_metadata_preserved "
          f"(SHIFR={spec.sample_id!r}, GEOMETRY={spec.geometry!r}, "
          f"DISTANCE={distance_str})")


def test_lsrm_spe_no_unknown_fields():
    """Verify that recently-added fields no longer appear as 'unknown'
    when iterating over our reference set."""
    # This is more of a sanity check than a test of any specific feature;
    # we just verify the readers don't crash on any of the reference
    # spectra.
    ref_dir = Path("detectors/Gamma-1S/reference_spectra/archive")
    n_files = 0
    for path in ref_dir.glob("*.spe"):
        try:
            spec = read_spectrum(str(path))
            assert spec.counts is not None and len(spec.counts) > 0
            assert spec.live_time >= 0
            n_files += 1
        except Exception as e:
            assert False, f"Failed to read {path.name}: {e}"
    print(f"  ✓ test_lsrm_spe_no_unknown_fields ({n_files} files read OK)")


if __name__ == "__main__":
    print("Running Lsrm SPE reader extended tests...\n")
    test_lsrm_spe_basic_calibration()
    test_lsrm_spe_energy_quality_parsed()
    test_lsrm_spe_peaks_table_parsed()
    test_lsrm_spe_zones_table_parsed()
    test_parse_peaks_table_isolated()
    test_parse_zones_table_isolated()
    test_lsrm_spe_metadata_preserved()
    test_lsrm_spe_no_unknown_fields()
    print("\n✓ All Lsrm SPE reader extended tests passed.")
