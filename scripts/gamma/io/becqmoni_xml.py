"""
BecqMoni / AtomSpectra XML writer.

Both BecqMoni (Am6er/BecqMoni, Nuclear edition) and AtomSpectra PRO emit
the same `ResultDataFile` → `ResultDataList` → `ResultData` structure.
SpecUtils calls this parser `RadiaCode`. Our existing reader for this
format lives in `gamma.io.atomspectra_xml.read_atomspectra_xml`; this
module adds the matching writer.

Structure produced (see `atomspectra_xml.py` docstring for full notes):

    <ResultDataFile>
      <FormatVersion>120920</FormatVersion>
      <ResultDataList>
        <ResultData>
          <SampleInfo>
            <Name/> <Location/> <Time/>
            <Weight>1</Weight> <Volume>1</Volume> <Note/>
          </SampleInfo>
          <DeviceConfigReference>
            <Name>...</Name> <Guid>...</Guid>
          </DeviceConfigReference>
          <StartTime>...</StartTime>
          <EndTime>...</EndTime>
          <EnergySpectrum>
            <NumberOfChannels>N</NumberOfChannels>
            <ChannelPitch>1</ChannelPitch>
            <EnergyCalibration>
              <PolynomialOrder>order</PolynomialOrder>
              <Coefficients>
                <Coefficient>c0</Coefficient>
                ...
              </Coefficients>
            </EnergyCalibration>
            <ValidPulseCount>...</ValidPulseCount>
            <TotalPulseCount>...</TotalPulseCount>
            <MeasurementTime>real_time_s</MeasurementTime>
            <LiveTime>live_time_s</LiveTime>
            <NumberOfSamples>0</NumberOfSamples>
            <Spectrum>
              <DataPoint>0</DataPoint> ...
            </Spectrum>
          </EnergySpectrum>
          [<BackgroundEnergySpectrum>...</BackgroundEnergySpectrum>]
          <Visible>true</Visible>
        </ResultData>
      </ResultDataList>
    </ResultDataFile>

Conventions matching the reader:
  - PolynomialOrder is **degree** (NOT degree+1). Number of coefficients
    = PolynomialOrder + 1.
  - Coefficients are low-to-high.
  - MeasurementTime is real time in seconds (float).
  - LiveTime is live time in seconds (float).
  - DataPoint values are integers.
"""

from __future__ import annotations

import uuid
from datetime import datetime
from typing import Optional
# SEC-01 (P1): This module is WRITER-ONLY — it accepts a Spectrum dataclass
# and emits XML via Element/SubElement/ElementTree.write. There is NO parser
# entry point (no ET.parse / ET.fromstring on user data) and therefore no
# attack surface for XXE / billion-laughs / DOCTYPE entity expansion.
# Stdlib xml.etree.ElementTree is correct here; defusedxml is parse-only and
# intentionally not imported. If a reader is ever added to this module, swap
# the parse call to `defusedxml.ElementTree.parse / fromstring` (mirroring
# atomspectra_xml.py / n42_2012.py / cpt_io.py / lsrm_library.py).
from xml.etree import ElementTree as ET

import numpy as np

from gamma.spectrum import Spectrum


FORMAT_VERSION = "120920"  # Same value AtomSpectra PRO writes.


def write_becqmoni_xml(
    spec: Spectrum,
    path: str,
    *,
    pretty: bool = True,
) -> None:
    """
    Write a Spectrum to BecqMoni/AtomSpectra ResultDataFile XML.

    Round-trip pairing: this format is identical to what
    `gamma.io.atomspectra_xml.read_atomspectra_xml` reads.

    If `spec.background_embedded` is set, it is written as a
    <BackgroundEnergySpectrum> sibling of <EnergySpectrum>.
    """
    root = ET.Element("ResultDataFile")

    fv = ET.SubElement(root, "FormatVersion")
    fv.text = FORMAT_VERSION

    rd_list = ET.SubElement(root, "ResultDataList")
    rd = ET.SubElement(rd_list, "ResultData")

    # --- SampleInfo ---
    sample = ET.SubElement(rd, "SampleInfo")
    ET.SubElement(sample, "Name").text = spec.sample_id or ""
    ET.SubElement(sample, "Location").text = spec.extras.get("location", "")
    ET.SubElement(sample, "Time").text = (
        _format_iso(spec.start_datetime) if spec.start_datetime else ""
    )
    ET.SubElement(sample, "Weight").text = str(
        spec.extras.get("lsrm_samplemass", "").split(";")[0] or "1"
    )
    ET.SubElement(sample, "Volume").text = str(
        spec.extras.get("lsrm_samplevolume", "").split(";")[0] or "1"
    )
    ET.SubElement(sample, "Note").text = spec.comments or ""

    # --- DeviceConfigReference ---
    dev = ET.SubElement(rd, "DeviceConfigReference")
    ET.SubElement(dev, "Name").text = spec.detector_id or "Unknown"
    ET.SubElement(dev, "Guid").text = spec.device_guid or str(uuid.uuid4())

    # --- Background link (filename hint only — actual data goes embedded) ---
    if spec.background_link:
        ET.SubElement(rd, "BackgroundSpectrumFile").text = spec.background_link

    # --- StartTime / EndTime ---
    # BUG-8 / Fix C3 (2026-06-02): BecqMoni reference files all carry a local
    # timezone offset (e.g. `+03:00`). `BecquerelMonitor/ResultData.cs` declares
    # these as `DateTime`, and the .NET deserialiser parses both with and
    # without TZ — but display semantics differ (Local vs UTC interpretation).
    # If the upstream datetime is naive (LSRM .spe MEASBEGIN provides no TZ),
    # we attach the local machine TZ so BecqMoni renders it correctly instead
    # of falling back to UTC.
    if spec.start_datetime is not None:
        ET.SubElement(rd, "StartTime").text = _format_iso_with_tz(spec.start_datetime)
    if spec.end_datetime is not None:
        ET.SubElement(rd, "EndTime").text = _format_iso_with_tz(spec.end_datetime)

    # --- EnergySpectrum (the main block) ---
    _append_energy_block(rd, spec, tag="EnergySpectrum")

    # --- BackgroundEnergySpectrum (optional) ---
    if spec.background_embedded is not None:
        _append_energy_block(rd, spec.background_embedded, tag="BackgroundEnergySpectrum")

    ET.SubElement(rd, "Visible").text = "true"

    if pretty:
        _indent(root)
    tree = ET.ElementTree(root)
    tree.write(str(path), encoding="utf-8", xml_declaration=True)


# BUG-BQ3 (2026-07-27) — BecqMoni evaluates degree 2, 3 and 4 only.
# `BecquerelMonitor/PolynomialEnergyCalibration.cs:ChannelToEnergy()` branches
# on polynomialOrder == 4 / 3 / 2; anything else falls through to
# `coefficients[1] * n + coefficients[0]`, a silent drop to a straight line.
# LSRM emits degree-5 calibrations for a handful of files and the top term is
# not negligible there — dropping it costs ~53 keV at channel 8192 on the
# HPGe backgrounds and far more on a curved NaI calibration.
#
# The writer still emits the source coefficients unchanged. Re-fitting to
# degree 4 was tried and rejected: on the well-behaved HPGe curves it lands
# within 0.21 keV, but on a strongly curved NaI calibration the best degree-4
# approximation is 57 keV off, and silently altering a stored calibration is
# worse than reporting it. Affected spectra are flagged in the per-class
# INDEX.json under `becqmoni_reads_as_linear`.
BQ_MAX_POLY_ORDER = 4


def _append_energy_block(
    parent: ET.Element,
    spec: Spectrum,
    *,
    tag: str,
) -> None:
    """Add an <EnergySpectrum> or <BackgroundEnergySpectrum> block."""
    block = ET.SubElement(parent, tag)
    counts = np.asarray(spec.counts, dtype=np.int64)
    n = int(len(counts))

    # F-RPT-05 / v1.18.29 — writer ВСЕГДА эмитит `len(counts)`, никогда
    # не использует pre-truncation `spec.n_channels_raw`. Иначе output
    # XML заявляет N каналов в заголовке, но в <Spectrum> их меньше —
    # BecqMoni reader падает или дочитывает нулями.
    ET.SubElement(block, "NumberOfChannels").text = str(n)
    ET.SubElement(block, "ChannelPitch").text = str(spec.channel_pitch or 1)

    # EnergyCalibration: PolynomialOrder = degree (= len(coefs) - 1).
    cal = ET.SubElement(block, "EnergyCalibration")
    if spec.energy_cal:
        # Emitted verbatim, including degree 5+ — see BQ_MAX_POLY_ORDER above
        # for what BecqMoni does with it. Operator decision 2026-07-27: keep
        # the source calibration intact and flag the affected files in the
        # per-class INDEX.json instead of re-fitting.
        coefs = list(spec.energy_cal)
        ET.SubElement(cal, "PolynomialOrder").text = str(len(coefs) - 1)
        parent_c = ET.SubElement(cal, "Coefficients")
        for c in coefs:
            ET.SubElement(parent_c, "Coefficient").text = f"{c:.10g}"
    else:
        # BUG-BQ1 (2026-07-03): BecqMoni PolynomialEnergyCalibration.ChannelToEnergy
        # reads coeffs[0] unconditionally → IndexOutOfRangeException on empty block.
        # Emit identity calibration (E=ch, 1 keV/channel) so BecqMoni opens the file.
        ET.SubElement(cal, "PolynomialOrder").text = "1"
        parent_c = ET.SubElement(cal, "Coefficients")
        ET.SubElement(parent_c, "Coefficient").text = "0"
        ET.SubElement(parent_c, "Coefficient").text = "1"

    # BUG-8 / Fix B6 (2026-06-02): ValidPulseCount and TotalPulseCount must be
    # emitted unconditionally. `BecquerelMonitor/EnergySpectrum.cs:218-219`
    # declares both as `long` with default 0; the reader tolerates absence,
    # but the BecqMoni UI relies on them to render the count-rate / total-pulses
    # widget. When the upstream format does not carry explicit pulse counters
    # (LSRM .spe in particular), we fall back to `sum(counts)` — for an
    # un-rejected spectrum this is the same number BecqMoni would compute
    # from `Spectrum/DataPoint` integration. We use the same value for both
    # to preserve the invariant `ValidPulseCount <= TotalPulseCount`.
    counts_sum_fallback = int(counts.sum())
    if spec.valid_pulse_count is not None:
        valid_pc = int(spec.valid_pulse_count)
    else:
        valid_pc = counts_sum_fallback
    if spec.total_pulse_count is not None:
        total_pc = int(spec.total_pulse_count)
    else:
        total_pc = counts_sum_fallback
    # Defensive: if the upstream provided both but somehow valid > total
    # (corrupted source), clamp to keep the BecqMoni invariant intact.
    if valid_pc > total_pc:
        valid_pc = total_pc
    ET.SubElement(block, "ValidPulseCount").text = str(valid_pc)
    ET.SubElement(block, "TotalPulseCount").text = str(total_pc)

    ET.SubElement(block, "MeasurementTime").text = f"{float(spec.real_time):.4f}"
    ET.SubElement(block, "LiveTime").text = f"{float(spec.live_time):.6f}"
    ET.SubElement(block, "NumberOfSamples").text = "0"

    # BUG-8 / Fix A2 (2026-06-02): DataPoint values are the RAW counts from
    # the source spectrum (e.g. LSRM `SPECTR=` block). No live/real scaling,
    # no dead-time correction, no normalisation. Any such correction is the
    # responsibility of the downstream analyser. The BecqMoni-spec audit
    # raised a ~0.99x discrepancy vs. a reference BecqMoni-saved .bq.xml of
    # the same .spe; our writer is a verified pass-through (see
    # `test_writer_counts_passthrough_no_scaling`), so any residual delta is
    # in the BecqMoni-side import path, not here.
    sp = ET.SubElement(block, "Spectrum")
    for v in counts:
        ET.SubElement(sp, "DataPoint").text = str(int(v))


# ============================================================================
# Helpers
# ============================================================================

def _format_iso(dt: datetime) -> str:
    if dt is None:
        return ""
    return dt.isoformat()


def _format_iso_with_tz(dt: datetime) -> str:
    """ISO-8601 timestamp with a guaranteed timezone offset.

    BUG-8 / Fix C3: BecqMoni reference .bq.xml files always carry a local
    timezone suffix (e.g. `+03:00`) in StartTime/EndTime. Naive datetimes
    are interpreted as local time and the machine offset is attached.
    Aware datetimes pass through their existing offset.
    """
    if dt is None:
        return ""
    if dt.tzinfo is None or dt.tzinfo.utcoffset(dt) is None:
        # Naive datetime — interpret as local wall-clock and attach the
        # machine's current TZ offset. `astimezone()` with no args uses the
        # local zone in Python 3.6+.
        dt = dt.astimezone()
    return dt.isoformat()


def _indent(elem: ET.Element, level: int = 0, step: str = "  ") -> None:
    """Pretty-print backport (mirrors n42_2012._indent)."""
    pad = "\n" + step * level
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = pad + step
        for i, child in enumerate(elem):
            _indent(child, level + 1, step)
            if i == len(elem) - 1:
                if not child.tail or not child.tail.strip():
                    child.tail = pad
            else:
                if not child.tail or not child.tail.strip():
                    child.tail = pad + step
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = pad
