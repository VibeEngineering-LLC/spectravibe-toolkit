"""
Regression test for T41 (BUG-40 (b) hardening).

Validates the standalone validator
``gamma.calibration.efficiency_provenance.check_efr_detector_match``,
which catches the silent CONTENT-fallback class: the .efr file was
loaded successfully, but its ``[detector;…]`` header / ``Detector=``
metadata records a DIFFERENT physical instance (serial-year) than the
spectrum's CONFIGNAME/DETECTOR fields.

Real BUG-40 (b) incident (2026-06-23):
``detectors/Gamma-1S/efficiency/Gamma-1S_NaI_63x63_USB_SN-01/
…Marinelli.efr`` is named after Gamma-1S but its header is
``УДС-ГЦ-63х63-USB №SN-01`` (Поверка-2024) while the spectrum's
CONFIGNAME is ``Гамма-1С №SN-02`` (Поверка-2016). Serial 0086 ≠
0221 → wrong physical instrument's efficiency curve → activity bias
of −96% to −97% on Am-241/Ti-44.

The existing ``cyrillic_to_latin_collision`` predicate only catches
PATH-level homoglyph; CONTENT-level mismatch was silent until this
validator was added.
"""
from __future__ import annotations

from gamma.calibration.efficiency_provenance import (
    extract_serial_year,
    check_efr_detector_match,
)


# ---------------------------------------------------------------------------
# extract_serial_year — regex-level coverage
# ---------------------------------------------------------------------------

def test_extract_serial_year_cyrillic_no_sign():
    assert extract_serial_year("Гамма-1С №SN-02") == ("0221", "16")


def test_extract_serial_year_latin_ascii_no():
    assert extract_serial_year("UDS-GTs-63x63-USB No SN-01") == ("0086", "16")


def test_extract_serial_year_latin_ascii_n_dot():
    assert extract_serial_year("Detector N-SN-01") == ("0086", "16")


def test_extract_serial_year_no_match_returns_none():
    assert extract_serial_year("Гамма-1С") is None
    assert extract_serial_year("Unknown detector") is None
    assert extract_serial_year("") is None
    assert extract_serial_year(None) is None


def test_extract_serial_year_three_digit_serial():
    assert extract_serial_year("Detector №123-16") == ("123", "16")


# ---------------------------------------------------------------------------
# check_efr_detector_match — integration with .efr parser
# ---------------------------------------------------------------------------

def _write_efr(tmp_path, detector_header: str, detector_metadata: str) -> str:
    efr = tmp_path / "test.efr"
    text = (
        f"[{detector_header};Marinelli-1L;Cs-137]\n"
        f"Detector={detector_metadata}\n"
        f"Volume,ml=1000\n"
        f"Distance,cm=0\n"
        f"Density,g/cm3=1.0\n"
        f"100=1.0e-02,5.0,Cs-137,1000,10,85.0\n"
    )
    # .efr is CP-1251 per format spec
    efr.write_bytes(text.encode("cp1251"))
    return str(efr)


def test_check_efr_detector_match_real_mismatch_case(tmp_path):
    """Real BUG-40 (b) reproduction: 0221 vs 0086 different physical units."""
    efr_path = _write_efr(
        tmp_path,
        detector_header="УДС-ГЦ-63х63-USB №SN-01",
        detector_metadata="УДС-ГЦ-63х63-USB №SN-01",
    )
    result = check_efr_detector_match(efr_path, "Гамма-1С №SN-02")
    assert result is not None
    assert result["code"] == "EFFICIENCY_DETECTOR_SERIAL_MISMATCH"
    assert result["expected_serial_year"] == ["0221", "16"]
    assert result["actual_serial_year"] == ["0086", "16"]
    assert result["efr_file_basename"] == "test.efr"
    assert "Гамма-1С" in result["expected_detector"]
    assert "SN-01" in result["actual_detector"]


def test_check_efr_detector_match_same_serial_no_warning(tmp_path):
    """When both sides agree → silent None (don't cry wolf)."""
    efr_path = _write_efr(
        tmp_path,
        detector_header="Гамма-1С №SN-02",
        detector_metadata="Гамма-1С №SN-02",
    )
    result = check_efr_detector_match(efr_path, "Гамма-1С №SN-02")
    assert result is None


def test_check_efr_detector_match_expected_no_serial_silent(tmp_path):
    """If expected detector has no extractable serial — silent None."""
    efr_path = _write_efr(
        tmp_path,
        detector_header="УДС-ГЦ-63х63-USB №SN-01",
        detector_metadata="УДС-ГЦ-63х63-USB №SN-01",
    )
    result = check_efr_detector_match(efr_path, "Gamma-1S")
    assert result is None


def test_check_efr_detector_match_actual_no_serial_silent(tmp_path):
    """If .efr Detector= has no extractable serial — silent None."""
    efr_path = _write_efr(
        tmp_path,
        detector_header="Some detector",
        detector_metadata="Some detector",
    )
    result = check_efr_detector_match(efr_path, "Гамма-1С №SN-02")
    assert result is None


def test_check_efr_detector_match_empty_inputs_silent():
    assert check_efr_detector_match("", "Гамма-1С №SN-02") is None
    assert check_efr_detector_match("nonexistent.efr", "") is None
    assert check_efr_detector_match("", "") is None


def test_check_efr_detector_match_unreadable_efr_silent(tmp_path):
    """Parser failure → silent None (validator never blocks loading)."""
    bad = tmp_path / "broken.efr"
    bad.write_bytes(b"\xff\xfe garbage no sections \x00\x00")
    result = check_efr_detector_match(str(bad), "Гамма-1С №SN-02")
    assert result is None


def test_check_efr_detector_match_serial_only_in_block_detector_field(tmp_path):
    """If Detector= metadata is missing but block detector is set, prefer block."""
    efr = tmp_path / "test.efr"
    text = (
        "[УДС-ГЦ-63х63-USB №SN-01;Marinelli-1L;Cs-137]\n"
        "Volume,ml=1000\n"
        "100=1.0e-02,5.0,Cs-137,1000,10,85.0\n"
    )
    efr.write_bytes(text.encode("cp1251"))
    result = check_efr_detector_match(str(efr), "Гамма-1С №SN-02")
    assert result is not None
    assert result["actual_serial_year"] == ["0086", "16"]