"""
Tests for the ASPECT (АСПЕКТ / Dubna) `.spc` binary reader.

Sample corpus: 4 unique real spectra from SpecUtils issue #47
(https://github.com/sandialabs/SpecUtils/issues/47), NaI(Tl) 25 mm² Ø25 mm
handheld probes measuring Cs-137 and Co-60 sources in January 2016.

All fixtures are uncalibrated (empty `energy[4][16]` / `fwhm[4][16]` blocks) —
downstream calibration is a detector-specific decision, not the reader's
business.

Ground truth was cross-checked byte-by-byte against the vendor header
`references/vendor/aspect/spheader_en.h` (translation of the Russian original
`spheader_ru_cp1251.h`).
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

import numpy as np
import pytest

from gamma.io.aspect_spc import (
    looks_like_aspect_spc,
    read_aspect_spc,
)
from gamma.io import format_registry as _fr


FIXTURE_DIR = Path(__file__).resolve().parent.parent / "data" / "aspect_spc"


# ---------------------------------------------------------------------------
# Reference sample data — see references/ASPECT_SPC_FORMAT_SPEC.md §11
# ---------------------------------------------------------------------------
SAMPLES = {
    "cs137_test2_25mm2.spc": {
        "n_channels": 8192,
        "sum_counts": 51065,
        "live_time_s": 7434.877,
        "real_time_s": 7435.634,
        "peak_ch": 74,
        "last_nonzero_ch": 6999,
        "start_datetime": datetime(2016, 1, 13, 21, 0, 13),
    },
    "cs137_test1.spc": {
        "n_channels": 8192,
        "sum_counts": 281957,
        "live_time_s": 5618.064,
        "real_time_s": 5623.009,
        "peak_ch": 77,
        "last_nonzero_ch": 3189,
        "start_datetime": datetime(2016, 1, 10, 23, 56, 44),
    },
    "co60_test2.spc": {
        "n_channels": 8192,
        "sum_counts": 244126,
        "live_time_s": 4872.964,
        "real_time_s": 4877.246,
        "peak_ch": 74,
        "last_nonzero_ch": 3189,
        "start_datetime": datetime(2016, 1, 10, 23, 44, 19),
    },
    "co60_cs137_25m2.spc": {
        "n_channels": 8192,
        "sum_counts": 270090,
        "live_time_s": 4969.813,
        "real_time_s": 4974.207,
        "peak_ch": 74,
        "last_nonzero_ch": 3192,
        "start_datetime": datetime(2016, 1, 14, 19, 55, 4),
    },
}


@pytest.mark.parametrize("fixture_name", sorted(SAMPLES.keys()))
def test_sniffer_recognizes_all_fixtures(fixture_name):
    path = FIXTURE_DIR / fixture_name
    head = path.read_bytes()[:1024]
    assert looks_like_aspect_spc(head), f"Sniffer missed {fixture_name}"


@pytest.mark.parametrize("fixture_name", sorted(SAMPLES.keys()))
def test_detect_format_returns_aspect_spc(fixture_name):
    path = FIXTURE_DIR / fixture_name
    assert _fr.detect_format(str(path)) == "aspect_spc"


@pytest.mark.parametrize("fixture_name", sorted(SAMPLES.keys()))
def test_reader_basic_fields(fixture_name):
    expected = SAMPLES[fixture_name]
    path = FIXTURE_DIR / fixture_name
    spec = read_aspect_spc(str(path))

    assert spec.source_format == "aspect_spc"
    assert spec.n_channels == expected["n_channels"]
    assert spec.n_channels_raw == expected["n_channels"]
    assert len(spec.counts) == expected["n_channels"]
    assert isinstance(spec.counts, np.ndarray)
    assert spec.counts.dtype == np.int64

    assert int(spec.counts.sum()) == expected["sum_counts"]
    assert spec.live_time == pytest.approx(expected["live_time_s"], abs=1e-3)
    assert spec.real_time == pytest.approx(expected["real_time_s"], abs=1e-3)
    assert int(spec.counts.argmax()) == expected["peak_ch"]

    nz = int((spec.counts > 0).nonzero()[0].max())
    assert nz == expected["last_nonzero_ch"]

    assert spec.start_datetime == expected["start_datetime"]


@pytest.mark.parametrize("fixture_name", sorted(SAMPLES.keys()))
def test_all_fixtures_are_uncalibrated(fixture_name):
    """All 4 SpecUtils issue #47 fixtures ship with empty energy/fwhm blocks."""
    path = FIXTURE_DIR / fixture_name
    spec = read_aspect_spc(str(path))
    assert spec.energy_cal is None
    assert spec.energy_cal_source == ""
    assert spec.extras.get("aspect_uncalibrated") is True
    assert "aspect_fwhm_coeffs" not in spec.extras


@pytest.mark.parametrize("fixture_name", sorted(SAMPLES.keys()))
def test_reader_extras_provenance(fixture_name):
    """Raw ms values preserved in extras for anti-hallucination cross-check."""
    path = FIXTURE_DIR / fixture_name
    spec = read_aspect_spc(str(path))

    live_ms = spec.extras["aspect_live_time_ms"]
    real_ms = spec.extras["aspect_real_time_ms"]
    assert isinstance(live_ms, int)
    assert isinstance(real_ms, int)
    assert live_ms / 1000.0 == pytest.approx(spec.live_time, abs=1e-6)
    assert real_ms / 1000.0 == pytest.approx(spec.real_time, abs=1e-6)

    assert spec.extras["aspect_buffer"] == spec.n_channels
    assert spec.extras["aspect_first"] == 0
    last_val = spec.extras["aspect_last"]
    assert last_val in (spec.n_channels - 1, 0)


def test_size_mismatch_rejected(tmp_path):
    """Reader must refuse a file whose size does not match 512 + N × 4."""
    src = FIXTURE_DIR / "cs137_test2_25mm2.spc"
    bogus = tmp_path / "bogus.spc"
    bogus.write_bytes(src.read_bytes()[:-8])
    with pytest.raises(ValueError, match="size mismatch"):
        read_aspect_spc(str(bogus))


def test_too_short_rejected(tmp_path):
    """Files shorter than 516 bytes are rejected outright."""
    bogus = tmp_path / "tiny.spc"
    bogus.write_bytes(b"\x00" * 100)
    with pytest.raises(ValueError, match="too short"):
        read_aspect_spc(str(bogus))


def test_sniffer_rejects_progress_magic():
    """Progress .spc (SF9I magic) must not trigger the ASPECT sniffer."""
    fake_progress = b"SF9I" + b"KEY=value " * 60
    fake_progress = fake_progress.ljust(1024, b" ")
    assert not looks_like_aspect_spc(fake_progress)


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_aspect_spc(str(tmp_path / "does_not_exist.spc"))
