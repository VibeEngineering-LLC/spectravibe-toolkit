"""
Golden-fixture tests for the AtomSpectra mobile FORMAT: 3 (.txt) reader.

Reuses the project's golden-fixture pattern (see test_reader_api.py).
Fixture: detectors/AtomSpectra/data/fixtures/Spectrum-2024-12-13_20-50-58-ДПР_радона.txt
Origin:  operator capture 2026-06-27 (#IO-1 / Task #55), 8192 ch, 5739.30 s live.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

import numpy as np

from gamma.io import format_registry as _fr
from gamma.io.atomspectra_txt import (
    looks_like_atomspectra_txt,
    read_atomspectra_txt,
)
from gamma.io.readers import read_spectrum
from gamma.spectrum import ENERGY_CEILING_KEV


TXT_FIXTURE = "detectors/AtomSpectra/data/fixtures/Spectrum-2024-12-13_20-50-58-ДПР_радона.txt"


# ---------------------------------------------------------------------------
# Sniffer
# ---------------------------------------------------------------------------

def test_sniffer_recognises_format3_banner():
    assert looks_like_atomspectra_txt(b"FORMAT: 3\nrest...") is True
    assert looks_like_atomspectra_txt(b"# not the magic") is False
    assert looks_like_atomspectra_txt(b"") is False
    print("  ✓ test_sniffer_recognises_format3_banner")


def test_registry_detects_format():
    fmt = _fr.detect_format(TXT_FIXTURE)
    assert fmt == "atomspectra_txt", f"expected atomspectra_txt, got {fmt!r}"
    print(f"  ✓ test_registry_detects_format ({fmt})")


# ---------------------------------------------------------------------------
# Direct reader
# ---------------------------------------------------------------------------

def test_reader_basic_fields():
    s = read_atomspectra_txt(TXT_FIXTURE)
    assert s.source_format == "atomspectra_txt"
    assert s.n_channels_raw == 8192
    assert len(s.counts) == s.n_channels
    assert s.counts.dtype == np.int64
    assert math.isclose(s.live_time, 5739.30, abs_tol=1e-6), (
        f"live_time={s.live_time} != 5739.30"
    )
    assert math.isclose(s.real_time, s.live_time, abs_tol=1e-9), (
        f"real_time fallback broken: real={s.real_time} live={s.live_time}"
    )
    assert s.extras["real_time_source"] == "fallback_from_live_time"
    print(f"  ✓ test_reader_basic_fields "
          f"(n={s.n_channels}, live={s.live_time}, real={s.real_time})")


def test_reader_counts_sanity_against_declared_sum():
    s = read_atomspectra_txt(TXT_FIXTURE)
    # Line 2 of the file declares "Counts: 86934". Verify the parsed array
    # sums to that value (plus the overflow channel if dropped).
    declared = 86934
    sum_kept = int(s.counts.sum())
    sum_overflow_dropped = sum_kept  # if dropped_overflow_count > 0 we lost a marker bucket
    assert sum_kept == declared, (
        f"counts.sum()={sum_kept} != declared {declared} on line 2 "
        f"(dropped_overflow_count={s.dropped_overflow_count})"
    )
    print(f"  ✓ test_reader_counts_sanity_against_declared_sum "
          f"(sum={sum_kept})")


def test_reader_calibration():
    s = read_atomspectra_txt(TXT_FIXTURE)
    assert s.energy_cal_degree == 3
    assert s.energy_cal is not None
    assert len(s.energy_cal) == 4
    # Low-to-high. Values verbatim from lines 12-15 of the fixture.
    expected = (-4.51604127884, 0.639000058174, 2.88695264317e-05, -4.09163902759e-09)
    for i, (got, exp) in enumerate(zip(s.energy_cal, expected)):
        assert math.isclose(got, exp, rel_tol=1e-10, abs_tol=1e-15), (
            f"coeff[{i}] = {got} != expected {exp}"
        )
    assert s.energy_cal_source == "stored"
    # Sanity: E(channel ~1023) → ~Cs-137 661.7 keV neighbourhood given a1≈0.639.
    e_at_1000 = s.channel_to_energy(1000)
    assert 600 < e_at_1000 < 700, f"E(1000)={e_at_1000} outside sanity window"
    print(f"  ✓ test_reader_calibration (deg=3, E(1000)≈{e_at_1000:.1f} keV)")


def test_reader_label_and_timestamps():
    s = read_atomspectra_txt(TXT_FIXTURE)
    assert s.sample_id == "ДПР радона", f"sample_id={s.sample_id!r}"
    assert s.comments == "ДПР радона", f"comments={s.comments!r}"
    assert s.start_datetime is not None, "start_datetime must be parsed"
    assert s.start_datetime.year == 2024
    assert s.start_datetime.month == 12
    assert s.start_datetime.day == 13
    assert s.start_datetime.hour == 20  # localised to +03:00 from summary
    assert s.start_datetime.minute == 50
    assert s.start_datetime.second == 58
    assert s.extras["timezone_offset_seconds"] == 3 * 3600
    print(f"  ✓ test_reader_label_and_timestamps "
          f"(start={s.start_datetime.isoformat()})")


def test_reader_header_field_passthrough():
    s = read_atomspectra_txt(TXT_FIXTURE)
    assert s.extras["format_version"] == 3
    assert s.extras["header_field_4"] == "0"
    assert s.extras["header_field_5"] == "0.0"
    assert s.extras["header_field_6"] == "0.0"
    assert "Counts: 86934" in s.extras["header_summary"]
    assert s.extras["trailing_extra_lines"] == 0
    print("  ✓ test_reader_header_field_passthrough")


# ---------------------------------------------------------------------------
# Energy ceiling parity with atomspectra_xml reader
# ---------------------------------------------------------------------------

def test_default_keeps_full_range():
    s = read_atomspectra_txt(TXT_FIXTURE)
    expected_n = s.n_channels_raw - s.dropped_overflow_count
    assert s.n_channels == expected_n, (
        f"default call must keep full decoded range: n_channels={s.n_channels} "
        f"!= n_channels_raw - dropped_overflow_count = {expected_n}"
    )
    assert s.extras["dropped_high_energy_count"] == 0
    assert s.energy_max_keV_kept is not None
    # With a1≈0.639 keV/ch and 8192 channels the upper edge sits well above 3 MeV.
    assert s.energy_max_keV_kept > ENERGY_CEILING_KEV, (
        f"fixture should reach above {ENERGY_CEILING_KEV} keV "
        f"(got e_max={s.energy_max_keV_kept})"
    )
    print(f"  ✓ test_default_keeps_full_range "
          f"(n={s.n_channels}, e_max={s.energy_max_keV_kept:.1f})")


def test_explicit_trim_at_3000():
    s = read_atomspectra_txt(TXT_FIXTURE, apply_energy_ceiling=True)
    assert s.energy_max_keV_kept is not None
    assert s.energy_max_keV_kept <= ENERGY_CEILING_KEV, (
        f"explicit-trim e_max {s.energy_max_keV_kept} > {ENERGY_CEILING_KEV}"
    )
    assert s.extras["dropped_high_energy_count"] > 0
    print(f"  ✓ test_explicit_trim_at_3000 "
          f"(n={s.n_channels}, e_max={s.energy_max_keV_kept:.1f})")


def test_custom_ceiling_shrinks_array():
    default = read_atomspectra_txt(TXT_FIXTURE)
    custom = read_atomspectra_txt(
        TXT_FIXTURE, apply_energy_ceiling=True, ceiling_keV=1500.0
    )
    assert custom.energy_max_keV_kept is not None
    assert custom.energy_max_keV_kept <= 1500.0
    assert custom.n_channels < default.n_channels
    print(f"  ✓ test_custom_ceiling_shrinks_array "
          f"(custom n={custom.n_channels} vs default n={default.n_channels})")


# ---------------------------------------------------------------------------
# Registry-level dispatch (the conversion entry point)
# ---------------------------------------------------------------------------

def test_read_spectrum_via_registry():
    s = read_spectrum(TXT_FIXTURE)
    assert s.source_format == "atomspectra_txt"
    assert s.n_channels_raw == 8192
    assert int(s.counts.sum()) == 86934
    print(f"  ✓ test_read_spectrum_via_registry (n={s.n_channels})")


# ---------------------------------------------------------------------------
# Negative paths
# ---------------------------------------------------------------------------

def test_rejects_non_format3(tmp_path=None):
    import tempfile
    with tempfile.NamedTemporaryFile(
        "w", suffix=".txt", encoding="utf-8", delete=False
    ) as fp:
        fp.write("# not a FORMAT: 3 file\n0\n1\n2\n")
        bad_path = fp.name
    try:
        try:
            read_atomspectra_txt(bad_path)
        except ValueError as exc:
            assert "FORMAT: 3" in str(exc)
            print("  ✓ test_rejects_non_format3 (ValueError as expected)")
            return
        raise AssertionError("expected ValueError for non-FORMAT: 3 file")
    finally:
        import os
        os.unlink(bad_path)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("Running AtomSpectra FORMAT: 3 reader tests...\n")
    test_sniffer_recognises_format3_banner()
    test_registry_detects_format()
    test_reader_basic_fields()
    test_reader_counts_sanity_against_declared_sum()
    test_reader_calibration()
    test_reader_label_and_timestamps()
    test_reader_header_field_passthrough()
    test_default_keeps_full_range()
    test_explicit_trim_at_3000()
    test_custom_ceiling_shrinks_array()
    test_read_spectrum_via_registry()
    test_rejects_non_format3()
    print("\nAll AtomSpectra FORMAT: 3 reader tests passed.")
