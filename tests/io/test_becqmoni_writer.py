"""F-RPT-05 / v1.18.29 — BecqMoni writer `n_channels` invariant.

Контракт: writer always emits `len(counts)` в `<NumberOfChannels>`,
никогда не использует pre-truncation `spec.n_channels_raw`.

Иначе output XML заявляет N каналов в заголовке, но в `<Spectrum>` их меньше —
BecqMoni reader либо падает с EOF, либо дочитывает нулями (data corruption).

Reproduction baseline (наблюдалось в v1.18.27):
  Th232_..._Маринелли_0cm.spe       — 1024 каналов сырых
  → analyze pipeline → spec.counts trimmed до 1003 (energy ceiling 3000 keV)
  → write_becqmoni_xml             → declared 1024 / actual 1003 (mismatch)
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path
from xml.etree import ElementTree as ET

REPO = Path(__file__).resolve().parent.parent.parent
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


TH232_FIXTURE = (
    REPO / "detectors" / "Gamma-1S" / "reference_spectra" / "reference_kits"
    / "Marinelli_1L" / "Th-232" / "Th232_420-7-17_Маринелли_0cm.spe"
)


def _read_1024_channel_spe():
    """Read the Th-232 1024-channel .spe with apply_energy_ceiling=False
    so the writer test sees the full 1024 channels (not the 1003 trim).
    """
    from gamma.io.lsrm_spe import read_lsrm_spe
    return read_lsrm_spe(str(TH232_FIXTURE), apply_energy_ceiling=False)


def test_writer_emits_len_counts_not_n_channels_raw():
    """Read 1024-channel .spe with ceiling=False → write BecqMoni XML →
    re-read via AtomSpectra parser → assert round-trip preserves 1024.
    """
    if not TH232_FIXTURE.exists():
        import pytest
        pytest.skip(f"fixture not present: {TH232_FIXTURE}")

    from gamma.io.becqmoni_xml import write_becqmoni_xml
    from gamma.io.atomspectra_xml import read_atomspectra_xml

    spec_in = _read_1024_channel_spe()
    assert spec_in.n_channels == 1024, (
        f"sanity: 1024-channel fixture expected, got n_channels={spec_in.n_channels}"
    )
    assert len(spec_in.counts) == 1024
    assert spec_in.n_channels_raw == 1024

    with tempfile.TemporaryDirectory() as td:
        out_path = Path(td) / "th232_test.bq.xml"
        write_becqmoni_xml(spec_in, str(out_path))

        # Re-read via the matching reader; disable ceiling on the reader
        # too so we test writer fidelity (1024 → 1024), not reader trim.
        spec_out = read_atomspectra_xml(
            str(out_path), apply_energy_ceiling=False,
        )

        # Round-trip invariant (F-RPT-05)
        assert spec_out.n_channels == 1024, (
            f"round-trip n_channels mismatch: in=1024 vs out={spec_out.n_channels}"
        )
        assert len(spec_out.counts) == 1024, (
            f"round-trip len(counts) mismatch: in=1024 vs out={len(spec_out.counts)}"
        )

        # Declared `<NumberOfChannels>` must equal len(counts)
        tree = ET.parse(str(out_path))
        root = tree.getroot()
        # Find first EnergySpectrum/NumberOfChannels (depth-first OK because the
        # main block is the first energy block).
        noc_elem = root.find(".//EnergySpectrum/NumberOfChannels")
        assert noc_elem is not None, "writer did not emit <NumberOfChannels>"
        declared = int(noc_elem.text or 0)
        assert declared == len(spec_in.counts), (
            f"declared NumberOfChannels={declared} != "
            f"len(counts)={len(spec_in.counts)} — invariant violation"
        )
        assert declared == 1024


def test_writer_invariant_after_truncation():
    """Если pipeline truncated `counts` (например 1024 → 1003 от energy
    ceiling), writer must emit 1003, не 1024 (n_channels_raw).
    """
    if not TH232_FIXTURE.exists():
        import pytest
        pytest.skip(f"fixture not present: {TH232_FIXTURE}")

    from gamma.io.lsrm_spe import read_lsrm_spe
    from gamma.io.becqmoni_xml import write_becqmoni_xml

    # Read WITH energy ceiling — counts get trimmed from 1024 to ~1003.
    spec = read_lsrm_spe(str(TH232_FIXTURE), apply_energy_ceiling=True)
    assert spec.n_channels_raw == 1024, "sanity: raw should remember 1024"
    assert spec.n_channels < 1024, "sanity: counts should have been trimmed"
    n_actual = len(spec.counts)

    with tempfile.TemporaryDirectory() as td:
        out_path = Path(td) / "th232_trim.bq.xml"
        write_becqmoni_xml(spec, str(out_path))
        tree = ET.parse(str(out_path))
        noc_elem = tree.getroot().find(".//EnergySpectrum/NumberOfChannels")
        assert noc_elem is not None
        declared = int(noc_elem.text or 0)
        # F-RPT-05: must equal len(counts), NOT n_channels_raw=1024.
        assert declared == n_actual, (
            f"writer wrote {declared} but counts has {n_actual} entries "
            f"(n_channels_raw={spec.n_channels_raw}); F-RPT-05 violation"
        )
        assert declared != spec.n_channels_raw


# ============================================================================
# BUG-8 / 2026-06-02 — BecqMoni writer correctness against the .NET reference.
#
# Three fixes:
#   • A2 — counts must be a verbatim pass-through (no live/real scaling, no
#     dead-time correction). Cross-reference: `BecquerelMonitor/EnergySpectrum.cs`
#     stores `Spectrum/DataPoint` as raw `int` channel content.
#   • C3 — StartTime / EndTime must carry a local TZ offset. Reference BecqMoni
#     files always emit `…+HH:MM`. `BecquerelMonitor/ResultData.cs` declares
#     them as `DateTime`, which the .NET XML serialiser stamps with the local
#     offset.
#   • B6 — ValidPulseCount and TotalPulseCount must be emitted unconditionally.
#     `BecquerelMonitor/EnergySpectrum.cs:218-219` declares them as `long`
#     (default 0); the BecqMoni UI relies on them to display count rate.
# ============================================================================


def test_writer_counts_passthrough_no_scaling():
    """A2: writer must NOT apply live/real scaling or any other normalisation.

    Loads the Th-232 fixture, writes a .bq.xml, asserts:
      sum(spec.counts)  ==  sum(<DataPoint> values in the XML)
      len(spec.counts)  ==  count of <DataPoint> elements
    """
    if not TH232_FIXTURE.exists():
        import pytest
        pytest.skip(f"fixture not present: {TH232_FIXTURE}")

    from gamma.io.becqmoni_xml import write_becqmoni_xml

    spec = _read_1024_channel_spe()
    in_sum = int(sum(int(c) for c in spec.counts))
    in_len = len(spec.counts)

    with tempfile.TemporaryDirectory() as td:
        out_path = Path(td) / "th232_counts.bq.xml"
        write_becqmoni_xml(spec, str(out_path))
        tree = ET.parse(str(out_path))
        root = tree.getroot()
        dps = root.findall(".//EnergySpectrum/Spectrum/DataPoint")
        assert len(dps) == in_len, (
            f"DataPoint count mismatch: in={in_len} out={len(dps)}"
        )
        out_sum = sum(int(d.text or 0) for d in dps)
        # A2 invariant: BIT-FOR-BIT equality, not approximate. If this fails,
        # something in the writer is scaling counts (live/real ratio, dead-time,
        # normalisation by mass, etc.) — that's the bug to find.
        assert out_sum == in_sum, (
            f"counts sum mismatch (A2 violation): in={in_sum} out={out_sum} "
            f"delta={out_sum - in_sum} ratio={out_sum / in_sum if in_sum else 'inf':.6f}"
        )


def test_writer_starttime_has_timezone_offset():
    """C3: StartTime must carry an ISO-8601 timezone suffix.

    LSRM .spe MEASBEGIN is a naive datetime; the writer must attach the local
    machine TZ offset before serialising. Reference BecqMoni files all show
    `±HH:MM` after the seconds (and fractional seconds) field.
    """
    if not TH232_FIXTURE.exists():
        import pytest
        pytest.skip(f"fixture not present: {TH232_FIXTURE}")

    import re
    from gamma.io.becqmoni_xml import write_becqmoni_xml

    spec = _read_1024_channel_spe()
    assert spec.start_datetime is not None, "sanity: Th-232 fixture has MEASBEGIN"
    # Sanity: LSRM datetime is naive
    assert spec.start_datetime.tzinfo is None, (
        "sanity check: LSRM-parsed datetime is expected to be tz-naive"
    )

    with tempfile.TemporaryDirectory() as td:
        out_path = Path(td) / "th232_tz.bq.xml"
        write_becqmoni_xml(spec, str(out_path))
        tree = ET.parse(str(out_path))
        root = tree.getroot()

        st_elem = root.find(".//StartTime")
        assert st_elem is not None, "writer did not emit <StartTime>"
        st = st_elem.text or ""

        # ISO-8601 with TZ: YYYY-MM-DDThh:mm:ss[.fff]±hh:mm
        tz_re = re.compile(
            r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(\.\d+)?[+-]\d{2}:\d{2}$"
        )
        assert tz_re.match(st), (
            f"StartTime {st!r} does not carry a TZ offset (C3 violation)"
        )


def test_writer_emits_pulse_counts_unconditionally():
    """B6: ValidPulseCount and TotalPulseCount must be present in the XML,
    even when the upstream LSRM .spe did not record explicit pulse counters.

    Invariants enforced:
      • both elements present
      • both positive
      • ValidPulseCount <= TotalPulseCount
    """
    if not TH232_FIXTURE.exists():
        import pytest
        pytest.skip(f"fixture not present: {TH232_FIXTURE}")

    from gamma.io.becqmoni_xml import write_becqmoni_xml

    spec = _read_1024_channel_spe()
    # Sanity: LSRM .spe does not populate these
    assert spec.valid_pulse_count is None, (
        "sanity: LSRM .spe should not provide explicit ValidPulseCount"
    )
    assert spec.total_pulse_count is None, (
        "sanity: LSRM .spe should not provide explicit TotalPulseCount"
    )

    with tempfile.TemporaryDirectory() as td:
        out_path = Path(td) / "th232_pulses.bq.xml"
        write_becqmoni_xml(spec, str(out_path))
        tree = ET.parse(str(out_path))
        root = tree.getroot()

        vpc_elem = root.find(".//EnergySpectrum/ValidPulseCount")
        tpc_elem = root.find(".//EnergySpectrum/TotalPulseCount")
        assert vpc_elem is not None, (
            "B6 violation: <ValidPulseCount> not emitted"
        )
        assert tpc_elem is not None, (
            "B6 violation: <TotalPulseCount> not emitted"
        )
        vpc = int(vpc_elem.text or 0)
        tpc = int(tpc_elem.text or 0)
        assert vpc > 0, f"ValidPulseCount should be positive, got {vpc}"
        assert tpc > 0, f"TotalPulseCount should be positive, got {tpc}"
        assert vpc <= tpc, (
            f"BecqMoni invariant violated: ValidPulseCount={vpc} > "
            f"TotalPulseCount={tpc}"
        )
