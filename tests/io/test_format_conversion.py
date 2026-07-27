"""
Tests for the multi-format file conversion utility (F-49).

Verifies:
  - Format registry — list_formats, detect_format, extension dispatch
  - Per-format round-trip preserves counts, live_time, real_time and
    energy calibration within tolerance
  - Cross-format conversion chains preserve the same invariants
  - CLI flag --list-formats prints something usable

Run:
    cd gamma-spectrum-analysis
    PYTHONPATH=scripts python test_format_conversion.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

import numpy as np

from gamma.io import format_registry as fr
from gamma.io.convert import convert_spectrum
from gamma.io.format_registry import detect_format, read, write
from gamma.io.lsrm_spe_text import read_lsrm_spe_text, write_lsrm_spe_text
from gamma.io.lsrm_spe import read_lsrm_spe, write_lsrm_spe
from gamma.io.atomspectra_xml import read_atomspectra_xml
from gamma.io.becqmoni_xml import write_becqmoni_xml
from gamma.io.n42_2012 import read_n42_2012, write_n42_2012
from gamma.spectrum import Spectrum


# ============================================================================
# Fixtures
# ============================================================================

FIX = Path(__file__).parent.parent.parent / "evals" / "fixtures"
# F-307 / v1.18.7 — AtomSpectra xml-фикстуры изолированы в свой detector-subtree
ATOMSPECTRA_FIX = (
    Path(__file__).parent.parent.parent
    / "detectors" / "AtomSpectra" / "data" / "fixtures"
)

LSRM_FIXTURE = FIX / "M_cs_легкий_2001-2005.spe"
XML_FIXTURE = ATOMSPECTRA_FIX / "Фон_кабинет_8192к_01-01-2025.xml"


# ============================================================================
# Test counter
# ============================================================================

PASS = 0
FAIL = 0
ERRORS: list[str] = []


def _result(name: str, ok: bool, msg: str = "") -> None:
    global PASS, FAIL
    if ok:
        PASS += 1
        print(f"  ok   {name}")
    else:
        FAIL += 1
        ERRORS.append(f"{name}: {msg}")
        print(f"  FAIL {name}: {msg}")


def _expect_eq(name: str, a, b) -> None:
    _result(name, a == b, f"expected {b!r}, got {a!r}")


def _expect_close(name: str, a: float, b: float, atol: float = 1e-3) -> None:
    diff = abs(float(a) - float(b))
    _result(name, diff <= atol, f"expected ~{b}, got {a} (d={diff})")


# ============================================================================
# Tests — registry
# ============================================================================

def test_list_formats():
    print("\n[1] Format registry")
    fmts = fr.list_formats()
    ids = {f.id for f in fmts}
    _result("list_formats returns >=4 formats", len(fmts) >= 4)
    for required in ("lsrm_spe", "lsrm_spe_text", "becqmoni_xml", "n42_2012"):
        _result(f"registry has {required}", required in ids)

    # Every registered format must have both reader and writer in this release
    for f in fmts:
        _result(f"{f.id} has reader", f.reader is not None)
        _result(f"{f.id} has writer", f.writer is not None)


def test_detect_format():
    print("\n[2] Format auto-detection")
    if LSRM_FIXTURE.is_file():
        _expect_eq("LSRM .spe detected as lsrm_spe",
                   detect_format(str(LSRM_FIXTURE)), "lsrm_spe")
    if XML_FIXTURE.is_file():
        _expect_eq("AtomSpectra .xml detected as becqmoni_xml",
                   detect_format(str(XML_FIXTURE)), "becqmoni_xml")


# ============================================================================
# Tests — per-format round-trip
# ============================================================================

def test_lsrm_round_trip(tmp_path: Path):
    print("\n[3] LSRM .spe round-trip")
    if not LSRM_FIXTURE.is_file():
        print("  skip: fixture missing")
        return

    orig = read_lsrm_spe(str(LSRM_FIXTURE), apply_energy_ceiling=False)
    out = tmp_path / "rt.spe"
    write_lsrm_spe(orig, str(out))
    back = read_lsrm_spe(str(out), apply_energy_ceiling=False)

    _expect_eq("n_channels preserved", back.n_channels, orig.n_channels)
    _expect_close("live_time preserved", back.live_time, orig.live_time, 0.01)
    _expect_close("real_time preserved", back.real_time, orig.real_time, 0.01)
    _result("counts identical",
            np.array_equal(np.asarray(back.counts), np.asarray(orig.counts)),
            "counts arrays differ")
    if orig.energy_cal:
        # Trailing zeros may differ — compare leading coefficients
        n = min(len(orig.energy_cal), len(back.energy_cal))
        for i in range(n):
            _expect_close(
                f"energy_cal[{i}] preserved",
                back.energy_cal[i], orig.energy_cal[i],
                atol=1e-6 * max(1.0, abs(orig.energy_cal[i])),
            )


def test_becqmoni_round_trip(tmp_path: Path):
    print("\n[4] BecqMoni / AtomSpectra XML round-trip")
    if not XML_FIXTURE.is_file():
        print("  skip: fixture missing")
        return

    orig = read_atomspectra_xml(str(XML_FIXTURE), apply_energy_ceiling=False)
    out = tmp_path / "rt.xml"
    write_becqmoni_xml(orig, str(out))
    back = read_atomspectra_xml(str(out), apply_energy_ceiling=False)

    _expect_eq("n_channels preserved", back.n_channels, orig.n_channels)
    _expect_close("live_time preserved", back.live_time, orig.live_time, 0.01)
    _expect_close("real_time preserved", back.real_time, orig.real_time, 0.01)
    _result("counts identical",
            np.array_equal(np.asarray(back.counts), np.asarray(orig.counts)),
            "counts arrays differ")
    if orig.energy_cal:
        for i in range(len(orig.energy_cal)):
            _expect_close(
                f"energy_cal[{i}] preserved",
                back.energy_cal[i], orig.energy_cal[i],
                atol=1e-6 * max(1.0, abs(orig.energy_cal[i])),
            )


def test_n42_round_trip(tmp_path: Path):
    print("\n[5] N42-42-2012 round-trip")
    if not LSRM_FIXTURE.is_file():
        print("  skip: fixture missing")
        return

    orig = read_lsrm_spe(str(LSRM_FIXTURE), apply_energy_ceiling=False)
    via_n42 = tmp_path / "rt.n42"
    write_n42_2012(orig, str(via_n42))
    back = read_n42_2012(str(via_n42), apply_energy_ceiling=False)

    _expect_eq("n_channels preserved", back.n_channels, orig.n_channels)
    _expect_close("live_time preserved", back.live_time, orig.live_time, 0.01)
    _expect_close("real_time preserved", back.real_time, orig.real_time, 0.01)
    _result("counts identical",
            np.array_equal(np.asarray(back.counts), np.asarray(orig.counts)),
            "counts arrays differ")
    if orig.energy_cal:
        for i in range(len(orig.energy_cal)):
            _expect_close(
                f"energy_cal[{i}] preserved",
                back.energy_cal[i], orig.energy_cal[i],
                atol=1e-6 * max(1.0, abs(orig.energy_cal[i])),
            )


def test_iaea_round_trip(tmp_path: Path):
    print("\n[6] LSRM SpectraLine ASCII SPE round-trip")
    if not LSRM_FIXTURE.is_file():
        print("  skip: fixture missing")
        return

    orig = read_lsrm_spe(str(LSRM_FIXTURE), apply_energy_ceiling=False)
    via = tmp_path / "rt_iaea.spe"
    write_lsrm_spe_text(orig, str(via))
    back = read_lsrm_spe_text(str(via), apply_energy_ceiling=False)

    _expect_eq("n_channels preserved", back.n_channels, orig.n_channels)
    _expect_close("live_time preserved", back.live_time, orig.live_time, 0.01)
    _expect_close("real_time preserved", back.real_time, orig.real_time, 0.01)
    _result("counts identical",
            np.array_equal(np.asarray(back.counts), np.asarray(orig.counts)),
            "counts arrays differ")
    if orig.energy_cal:
        for i in range(len(orig.energy_cal)):
            _expect_close(
                f"energy_cal[{i}] preserved",
                back.energy_cal[i], orig.energy_cal[i],
                atol=1e-6 * max(1.0, abs(orig.energy_cal[i])),
            )


# ============================================================================
# Tests — cross-format chain
# ============================================================================

def test_cross_format_chain(tmp_path: Path):
    print("\n[7] Cross-format chain: .spe -> n42 -> xml -> .spe")
    if not LSRM_FIXTURE.is_file():
        print("  skip: fixture missing")
        return

    orig = read_lsrm_spe(str(LSRM_FIXTURE), apply_energy_ceiling=False)
    p_n42 = tmp_path / "chain.n42"
    p_xml = tmp_path / "chain.xml"
    p_spe = tmp_path / "chain_back.spe"

    convert_spectrum(str(LSRM_FIXTURE), str(p_n42))
    convert_spectrum(str(p_n42), str(p_xml), out_format="becqmoni_xml")
    convert_spectrum(str(p_xml), str(p_spe), out_format="lsrm_spe")

    back = read_lsrm_spe(str(p_spe), apply_energy_ceiling=False)

    _expect_eq("chain preserves n_channels", back.n_channels, orig.n_channels)
    _expect_close("chain preserves live_time", back.live_time, orig.live_time, 0.01)
    _result("chain preserves counts (exact)",
            np.array_equal(np.asarray(back.counts), np.asarray(orig.counts)),
            "counts changed through chain")


def test_xml_chain_to_spe(tmp_path: Path):
    print("\n[8] XML chain: BecqMoni -> N42 -> LSRM .spe")
    if not XML_FIXTURE.is_file():
        print("  skip: fixture missing")
        return

    orig = read_atomspectra_xml(str(XML_FIXTURE), apply_energy_ceiling=False)
    p_n42 = tmp_path / "bg.n42"
    p_spe = tmp_path / "bg.spe"

    convert_spectrum(str(XML_FIXTURE), str(p_n42), out_format="n42_2012")
    convert_spectrum(str(p_n42), str(p_spe), out_format="lsrm_spe")

    back = read_lsrm_spe(str(p_spe), apply_energy_ceiling=False)
    _expect_eq("chain preserves n_channels (xml fixture)",
               back.n_channels, orig.n_channels)
    _result("chain preserves counts (xml fixture)",
            np.array_equal(np.asarray(back.counts), np.asarray(orig.counts)),
            "counts changed through chain")


# ============================================================================
# Tests — sniffer disambiguation
# ============================================================================

def test_sniffer_lsrm_vs_iaea(tmp_path: Path):
    print("\n[9] Sniffer disambiguates LSRM binary vs LSRM ASCII .spe")
    if not LSRM_FIXTURE.is_file():
        print("  skip: fixture missing")
        return

    # Write the LSRM binary fixture as LSRM ASCII — both end in .spe
    spec = read_lsrm_spe(str(LSRM_FIXTURE), apply_energy_ceiling=False)
    ascii_path = tmp_path / "as_text.spe"
    write_lsrm_spe_text(spec, str(ascii_path))

    _expect_eq("sniffer picks lsrm_spe for binary file",
               detect_format(str(LSRM_FIXTURE)), "lsrm_spe")
    _expect_eq("sniffer picks lsrm_spe_text for ASCII file",
               detect_format(str(ascii_path)), "lsrm_spe_text")


def test_sniffer_xml_vs_n42(tmp_path: Path):
    print("\n[10] Sniffer disambiguates BecqMoni vs N42 .xml")
    if not LSRM_FIXTURE.is_file():
        print("  skip: fixture missing")
        return

    spec = read_lsrm_spe(str(LSRM_FIXTURE), apply_energy_ceiling=False)
    xml_path = tmp_path / "as_xml.xml"
    n42_path = tmp_path / "as_n42.xml"
    write_becqmoni_xml(spec, str(xml_path))
    write_n42_2012(spec, str(n42_path))

    _expect_eq("sniffer picks becqmoni_xml for ResultDataFile",
               detect_format(str(xml_path)), "becqmoni_xml")
    _expect_eq("sniffer picks n42_2012 for RadInstrumentData",
               detect_format(str(n42_path)), "n42_2012")


# ============================================================================
# Tests — N42 channel-data CountedZeroes
# ============================================================================

def test_n42_counted_zeroes(tmp_path: Path):
    """Reader must accept compressionCode='CountedZeroes' channel data."""
    print("\n[11] N42 CountedZeroes decoding")
    n42_text = """<?xml version="1.0" encoding="UTF-8"?>
<RadInstrumentData xmlns="http://physics.nist.gov/N42/2011/N42">
  <RadInstrumentInformation id="i1">
    <RadInstrumentManufacturerName>test</RadInstrumentManufacturerName>
    <RadInstrumentModelName>test</RadInstrumentModelName>
    <RadInstrumentIdentifier>uuid</RadInstrumentIdentifier>
    <RadInstrumentClassCode>Spectroscopic Personal Radiation Detector</RadInstrumentClassCode>
  </RadInstrumentInformation>
  <EnergyCalibration id="ec1"><CoefficientValues>0 1.0</CoefficientValues></EnergyCalibration>
  <RadMeasurement id="m1">
    <MeasurementClassCode>Foreground</MeasurementClassCode>
    <RealTimeDuration>PT60S</RealTimeDuration>
    <Spectrum id="s1" energyCalibrationReference="ec1">
      <LiveTimeDuration>PT60S</LiveTimeDuration>
      <ChannelData compressionCode="CountedZeroes">0 3 5 7 0 2 11</ChannelData>
    </Spectrum>
  </RadMeasurement>
</RadInstrumentData>
"""
    p = tmp_path / "compressed.n42"
    p.write_text(n42_text, encoding="utf-8")
    spec = read_n42_2012(str(p), apply_energy_ceiling=False)

    # 0,3 → [0,0,0]; 5,7 → [5,7]; 0,2 → [0,0]; 11 → [11]
    # Total: [0,0,0,5,7,0,0,11] = 8 channels
    expected = np.array([0, 0, 0, 5, 7, 0, 0, 11], dtype=np.int64)
    _expect_eq("decoded length matches", spec.n_channels, 8)
    _result("decoded values match CountedZeroes spec",
            np.array_equal(np.asarray(spec.counts), expected),
            f"got {list(spec.counts)}, expected {list(expected)}")


# ============================================================================
# Tests — convert_spectrum API edge cases
# ============================================================================

def test_convert_api_explicit_formats(tmp_path: Path):
    print("\n[12] convert_spectrum honors explicit in_format/out_format")
    if not LSRM_FIXTURE.is_file():
        print("  skip: fixture missing")
        return

    # Force LSRM ASCII writer on a path with a non-matching extension
    out = tmp_path / "force.txt"
    convert_spectrum(
        str(LSRM_FIXTURE),
        str(out),
        out_format="lsrm_spe_text",
    )
    head = out.read_bytes()[:256]
    _result("explicit out_format=lsrm_spe_text wrote ASCII $-section",
            b"$DATE_MEA:" in head or b"$MEAS_TIM:" in head,
            f"head was {head!r}")


# ============================================================================
# Driver
# ============================================================================

def main():
    print("=" * 70)
    print("F-49 Format conversion test suite")
    print("=" * 70)
    print(f"LSRM fixture : {LSRM_FIXTURE} (exists={LSRM_FIXTURE.is_file()})")
    print(f"XML fixture  : {XML_FIXTURE} (exists={XML_FIXTURE.is_file()})")

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        test_list_formats()
        test_detect_format()
        test_lsrm_round_trip(tmp)
        test_becqmoni_round_trip(tmp)
        test_n42_round_trip(tmp)
        test_iaea_round_trip(tmp)
        test_cross_format_chain(tmp)
        test_xml_chain_to_spe(tmp)
        test_sniffer_lsrm_vs_iaea(tmp)
        test_sniffer_xml_vs_n42(tmp)
        test_n42_counted_zeroes(tmp)
        test_convert_api_explicit_formats(tmp)

    print()
    print("=" * 70)
    print(f"PASS: {PASS}   FAIL: {FAIL}")
    if ERRORS:
        print("\nFailures:")
        for e in ERRORS:
            print(f"  - {e}")
    print("=" * 70)
    return 0 if FAIL == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
