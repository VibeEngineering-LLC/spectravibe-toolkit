"""
AtomSpectra XML reader.

Format: ResultDataFile produced by AtomSpectra PRO and Lsrm SpectraLine.

Container structure (single ResultData per file in practice):

    <ResultDataFile>
      <FormatVersion>120920</FormatVersion>
      <ResultDataList>
        <ResultData>
          <SampleInfo>
            <Name/> <Location/> <Time>...</Time>
            <Weight>1</Weight> <Volume>1</Volume> <Note/>
          </SampleInfo>
          <DeviceConfigReference>
            <Name>AtomSpectra PRO 8192к</Name>
            <Guid>...</Guid>
          </DeviceConfigReference>
          <ROIConfigReference/>
          <BackgroundSpectrumFile>name.xml</BackgroundSpectrumFile>
          <StartTime>...</StartTime>
          <EndTime>...</EndTime>
          <PresetTime>360000000</PresetTime>
          <EnergySpectrum>
            <NumberOfChannels>8192</NumberOfChannels>
            <ChannelPitch>1</ChannelPitch>
            <EnergyCalibration>
              <PolynomialOrder>3</PolynomialOrder>
              <Coefficients>
                <Coefficient>...</Coefficient> × (order + 1)
              </Coefficients>
            </EnergyCalibration>
            <ValidPulseCount>...</ValidPulseCount>
            <TotalPulseCount>...</TotalPulseCount>
            <MeasurementTime>252678</MeasurementTime>   ← real time, seconds
            <LiveTime>252671.7348...</LiveTime>          ← live time, seconds
            <NumberOfSamples>0</NumberOfSamples>
            <Spectrum>
              <DataPoint>0</DataPoint> ... × NumberOfChannels
            </Spectrum>
          </EnergySpectrum>
          [<BackgroundEnergySpectrum>... (sample files only) ...</BackgroundEnergySpectrum>]
          <Visible>true</Visible>
          <PulseCollection>...</PulseCollection>
          [<SimpleSqrtFwhmCalibration>...</SimpleSqrtFwhmCalibration>  (sample files only)]
        </ResultData>
      </ResultDataList>
    </ResultDataFile>

Key facts established from real files:

  1. **PolynomialOrder = degree** (NOT degree − 1). Number of coefficients
     = PolynomialOrder + 1. Verified: degree 4 → 5 coeffs (background file),
     degree 3 → 4 coeffs (sample file).

  2. **Last channel is an overflow marker**. Real value orders of magnitude
     above the local tail (e.g. 128939 with tail ~1-7, or 48154 with tail ~1-2).
     We drop it. This is per the plan and verified on both files.

  3. **BackgroundSpectrumFile is a filename hint, not a payload**. The
     actual background counts are EMBEDDED in <BackgroundEnergySpectrum>
     (sample files), which has its own EnergyCalibration that may differ
     from the sample's (per SKILL.md §1.4: background may have been
     measured with a different gain settling).

  4. **Coefficients are low-to-high**: E(N) = c[0] + c[1]·N + c[2]·N² + ...
     Verified by comparing with peaks at known channels.

  5. **MeasurementTime is real time in seconds; LiveTime is live time in
     seconds**. Both as in their respective tags. PresetTime is in some
     internal ticks and is ignored.

  6. **Energies above ENERGY_CEILING_KEV are dropped at read time.** This
     trims the counts array to the channels whose calibrated energy lies
     below the ceiling (3000 keV per project scope).

  7. **SimpleSqrtFwhmCalibration**, when present, carries the vendor's
     FWHM calibration: list of calibration peaks with (channel, energy,
     fwhm_channels), polynomial coefficients, peak-type (1 = Hypermet),
     left/right tail parameters, χ²/N_dof. We expose it for diagnostic
     and for downstream FWHM(E) verification, NOT trust it blindly
     (per SKILL.md §6 stored-calibration policy).
"""

from __future__ import annotations

import os
import re
from datetime import datetime
from pathlib import Path
from typing import Optional
from xml.etree import ElementTree as ET

# SEC-01 (P1) hardening: parse untrusted user .xml via defusedxml to block
# billion-laughs / quadratic-blowup / DOCTYPE entity DoS. Stdlib ET is kept
# for Element/SubElement construction and type hints (defusedxml is parse-only).
from defusedxml.ElementTree import parse as _safe_parse  # noqa: E402

import numpy as np

from gamma.spectrum import (
    ENERGY_CEILING_KEV,
    FwhmCalPeak,
    Spectrum,
    StoredFwhmCalibration,
)
from gamma.io.filename_hints import parse_filename


# ============================================================================
# AtomSpectra XML
# ============================================================================

def _parse_iso_datetime(text: str) -> Optional[datetime]:
    """Parse ISO 8601 with timezone; falls back to naive datetime."""
    if not text:
        return None
    try:
        # Python <3.11 fromisoformat doesn't accept trailing 'Z' or offsets like
        # +03:00 with subseconds longer than 6 digits. AtomSpectra writes
        # microseconds with 7 digits. Strip surplus fractional digits.
        # Example: "2025-12-23T22:45:14.2427895+03:00"
        m = re.match(
            r"(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})"
            r"(\.\d{1,7})?"
            r"([+-]\d{2}:\d{2})?$",
            text.strip(),
        )
        if m:
            base = m.group(1)
            frac = m.group(2) or ""
            tz = m.group(3) or ""
            # Truncate fractional to 6 digits (microseconds)
            if frac:
                frac = frac[:7]  # ".XXXXXX"
            iso = base + frac + tz
            try:
                return datetime.fromisoformat(iso)
            except ValueError:
                # Drop fractional altogether if Python still chokes
                return datetime.fromisoformat(base + tz)
        return datetime.fromisoformat(text.strip())
    except (ValueError, TypeError):
        return None


def _text(elem: Optional[ET.Element], default: str = "") -> str:
    if elem is None or elem.text is None:
        return default
    return elem.text.strip()


def _int(elem: Optional[ET.Element], default: Optional[int] = None) -> Optional[int]:
    t = _text(elem, "")
    if t == "":
        return default
    try:
        return int(t)
    except ValueError:
        try:
            return int(float(t))
        except ValueError:
            return default


def _float(elem: Optional[ET.Element], default: Optional[float] = None) -> Optional[float]:
    t = _text(elem, "")
    if t == "":
        return default
    try:
        return float(t)
    except ValueError:
        return default


def _parse_energy_calibration(node: ET.Element) -> tuple:
    """
    Parse <EnergyCalibration> block.

    Returns: (coefficients_low_to_high, degree).
    Empty tuple and None if absent.
    """
    if node is None:
        return (), None

    order_elem = node.find("PolynomialOrder")
    degree = _int(order_elem)

    coefs_parent = node.find("Coefficients")
    if coefs_parent is None:
        return (), degree

    coefs = []
    for c in coefs_parent.findall("Coefficient"):
        v = _float(c)
        if v is not None:
            coefs.append(v)

    return tuple(coefs), degree


def _parse_spectrum_array(node: ET.Element) -> np.ndarray:
    """
    Parse <Spectrum><DataPoint>...</DataPoint>* into a numpy array.

    Uses iterparse-equivalent direct walk; ElementTree is fast enough for
    our scale (8192 channels × ~10 bytes per element ≈ 80 KB DOM, in the
    600-KB files we tested).
    """
    if node is None:
        return np.array([], dtype=np.int64)
    points = node.findall("DataPoint")
    arr = np.empty(len(points), dtype=np.int64)
    for i, dp in enumerate(points):
        try:
            arr[i] = int(dp.text)
        except (TypeError, ValueError):
            arr[i] = 0
    return arr


def _detect_and_drop_overflow(counts: np.ndarray) -> tuple:
    """
    The last channel of MCA output usually accumulates an overflow marker
    that swamps every event above the full-scale energy. Heuristic:
    if the last channel is at least 100× the 95th percentile of the upper
    tail (last 1% of channels, excluding itself), it's overflow — drop it.

    Returns (trimmed_counts, dropped_count).
    """
    if len(counts) < 100:
        return counts, 0

    tail_start = max(1, int(0.99 * len(counts)))
    tail = counts[tail_start:-1]
    if len(tail) == 0:
        return counts, 0

    p95 = np.percentile(tail, 95) if tail.size else 0
    last = counts[-1]
    # Robust check: overflow if last >= 100x the upper-tail 95th percentile
    # AND >= 1000 in absolute count.
    if last >= max(1000, 100 * max(p95, 1)):
        return counts[:-1].copy(), 1
    return counts, 0


def _channel_energies(coeffs: tuple, n: int) -> Optional[np.ndarray]:
    """
    Evaluate the calibrated energy of every channel 0..n-1 via Horner's
    method on the low-to-high coefficient list. Returns None if no usable
    calibration (need at least a0 and a1).
    """
    if not coeffs or len(coeffs) < 2 or n == 0:
        return None
    channels = np.arange(n, dtype=np.float64)
    energies = np.zeros(n, dtype=np.float64)
    for c in reversed(coeffs):
        energies = energies * channels + c
    return energies


def _apply_energy_ceiling(counts: np.ndarray, coeffs: tuple,
                          ceiling_keV: float) -> tuple:
    """
    Drop channels whose calibrated energy exceeds the ceiling.

    Returns (trimmed_counts, n_dropped_high_energy, energy_max_kept).
    If no calibration, no trimming and energy_max_kept is None.
    """
    energies = _channel_energies(coeffs, len(counts))
    if energies is None:
        return counts, 0, None

    n = len(counts)
    mask = energies <= ceiling_keV
    if mask.all():
        return counts, 0, float(energies[-1])

    # First channel where energy exceeds ceiling
    first_above = int(np.argmax(~mask)) if (~mask).any() else n
    kept = counts[:first_above].copy()
    return (
        kept,
        n - first_above,
        float(energies[first_above - 1]) if first_above > 0 else None,
    )


def _parse_calibration_peaks(node: ET.Element) -> list:
    """Parse <CalibrationPeaks><Peak>...</Peak></CalibrationPeaks>."""
    out = []
    if node is None:
        return out
    for p in node.findall("Peak"):
        ch = _int(p.find("Channel"))
        e = _float(p.find("Energy"))
        fw = _float(p.find("FWHM"))
        if ch is None or e is None or fw is None:
            continue
        out.append(FwhmCalPeak(channel=ch, energy_keV=e, fwhm_channels=fw))
    return out


def _parse_fwhm_calibration(node: ET.Element) -> Optional[StoredFwhmCalibration]:
    """
    Parse a SimpleSqrtFwhmCalibration block.

    Returns None if the node is absent or empty. Otherwise a populated
    StoredFwhmCalibration.
    """
    if node is None:
        return None

    peaks = _parse_calibration_peaks(node.find("CalibrationPeaks"))

    coefs_parent = node.find("Coefficients")
    coefs = ()
    if coefs_parent is not None:
        cs = []
        for c in coefs_parent.findall("Coefficient"):
            v = _float(c)
            if v is not None:
                cs.append(v)
        coefs = tuple(cs)

    return StoredFwhmCalibration(
        calibration_peaks=peaks,
        coefficients=coefs,
        peak_type=_int(node.find("PeakType")),
        left_tail=_float(node.find("ExpGaussExpLeftTail")),
        right_tail=_float(node.find("ExpGaussExpRightTail")),
        chi2_per_dof=_float(node.find("Chi2pNdp")),
        model="SimpleSqrtFwhm",
    )


def _parse_energy_spectrum_block(
    block: ET.Element,
    source_path: str,
    is_background: bool,
    *,
    apply_energy_ceiling: bool = False,
    ceiling_keV: Optional[float] = None,
) -> Optional[Spectrum]:
    """
    Parse one <EnergySpectrum> or <BackgroundEnergySpectrum> element into
    a Spectrum dataclass. Handles overflow drop and energy ceiling.

    The keyword-only `apply_energy_ceiling` / `ceiling_keV` mirror the
    public reader signature; see `read_atomspectra_xml` for semantics.
    """
    if block is None:
        return None

    n_channels_raw_elem = block.find("NumberOfChannels")
    n_channels_raw = _int(n_channels_raw_elem) or 0

    channel_pitch = _int(block.find("ChannelPitch"), default=1) or 1

    coeffs, deg = _parse_energy_calibration(block.find("EnergyCalibration"))

    valid_count = _int(block.find("ValidPulseCount"))
    total_count = _int(block.find("TotalPulseCount"))

    real_time = _float(block.find("MeasurementTime")) or 0.0
    live_time = _float(block.find("LiveTime")) or 0.0

    # Format compatibility:
    #   - Newer AtomSpectra firmware writes LiveTime explicitly (in seconds,
    #     fractional) — use it directly.
    #   - Older AtomSpectra firmware writes <LiveTime>0</LiveTime> as a
    #     placeholder; live_time is not recorded and must be approximated
    #     by MeasurementTime (i.e., assume dead time is negligible).
    # Both formats are accepted. When fallback is used, we flag the
    # spectrum so downstream code knows dead-time correction is
    # unavailable on this file (any dead-time correction step should
    # either skip it or rely on ValidPulseCount/TotalPulseCount ratio).
    live_time_from_fallback = False
    if live_time == 0.0 and real_time > 0.0:
        live_time = real_time
        live_time_from_fallback = True

    counts = _parse_spectrum_array(block.find("Spectrum"))

    # Overflow trim
    counts, dropped_overflow = _detect_and_drop_overflow(counts)

    # Energy ceiling trim (3 MeV per project scope, overridable per call)
    if apply_energy_ceiling:
        ceiling = ceiling_keV if ceiling_keV is not None else ENERGY_CEILING_KEV
        counts, dropped_high_e, e_max = _apply_energy_ceiling(
            counts, coeffs, ceiling
        )
    else:
        dropped_high_e = 0
        energies = _channel_energies(coeffs, len(counts))
        e_max = float(energies[-1]) if energies is not None else None

    spec = Spectrum(
        counts=counts,
        live_time=live_time,
        real_time=real_time,
        source_path=source_path,
        source_format="atomspectra_xml",
        is_background=is_background,
        valid_pulse_count=valid_count,
        total_pulse_count=total_count,
        dropped_overflow_count=dropped_overflow,
        n_channels_raw=n_channels_raw,
        n_channels=len(counts),
        channel_pitch=channel_pitch,
        energy_cal=coeffs if coeffs else None,
        energy_cal_degree=deg,
        energy_cal_source="stored" if coeffs else "",
        energy_max_keV_kept=e_max,
    )
    spec.extras["dropped_high_energy_count"] = dropped_high_e
    spec.extras["live_time_from_fallback"] = live_time_from_fallback
    return spec


def read_atomspectra_xml(
    path: str,
    *,
    parse_background: bool = True,
    apply_energy_ceiling: bool = False,
    ceiling_keV: Optional[float] = None,
) -> Spectrum:
    """
    Read an AtomSpectra XML file. Returns the primary Spectrum.

    Args:
        path: filesystem path to the .xml file.
        parse_background: if True (default), the embedded
            BackgroundEnergySpectrum (when present) is parsed and
            attached as `spec.background_embedded` with is_background=True.
        apply_energy_ceiling: if True, channels whose calibrated energy
            exceeds the ceiling are dropped from `counts`. The same
            policy is applied to the embedded background spectrum so the
            two arrays stay in register. **Default False since BUG-9
            (v1.18.32, 2026-06-03):** the reader keeps every decoded
            channel by default to avoid silent channel loss; opt in
            explicitly when the trim is desired.
        ceiling_keV: per-call override of `ENERGY_CEILING_KEV` (3000 keV
            by project scope). When None, the module constant is used.
            Ignored if `apply_energy_ceiling` is False.

    Raises:
      FileNotFoundError if path does not exist.
      ValueError if the file is not a valid ResultDataFile.
    """
    path = str(path)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"AtomSpectra XML not found: {path}")

    # ElementTree handles UTF-8 with BOM transparently (it strips the BOM
    # before parsing). lxml would do the same.
    # SEC-01: _safe_parse blocks DOCTYPE / billion-laughs entity expansion.
    tree = _safe_parse(path)
    root = tree.getroot()

    if root.tag != "ResultDataFile":
        raise ValueError(
            f"Not an AtomSpectra ResultDataFile (root is {root.tag!r}): {path}"
        )

    rd_list = root.find("ResultDataList")
    if rd_list is None:
        raise ValueError(f"ResultDataList missing in {path}")
    rd_nodes = rd_list.findall("ResultData")
    if not rd_nodes:
        raise ValueError(f"No ResultData entries in {path}")

    # Use the first ResultData. Per the plan, multi-ResultData is rare
    # and we flag it.
    rd = rd_nodes[0]
    multi_results = len(rd_nodes) > 1

    # ----- Sample info -----
    si = rd.find("SampleInfo")
    sample_name = _text(si.find("Name")) if si is not None else ""
    sample_location = _text(si.find("Location")) if si is not None else ""
    sample_note = _text(si.find("Note")) if si is not None else ""
    sample_time = _parse_iso_datetime(_text(si.find("Time"))) if si is not None else None

    # ----- Device -----
    dc = rd.find("DeviceConfigReference")
    device_name = _text(dc.find("Name")) if dc is not None else ""
    device_guid = _text(dc.find("Guid")) if dc is not None else ""

    # ----- Background link (string only) -----
    bsf = rd.find("BackgroundSpectrumFile")
    background_link = _text(bsf) if bsf is not None else ""
    if background_link == "":
        background_link = None

    # ----- Times -----
    start_dt = _parse_iso_datetime(_text(rd.find("StartTime")))
    end_dt = _parse_iso_datetime(_text(rd.find("EndTime")))

    # ----- Primary energy spectrum -----
    es = rd.find("EnergySpectrum")
    spec = _parse_energy_spectrum_block(
        es, source_path=path, is_background=False,
        apply_energy_ceiling=apply_energy_ceiling,
        ceiling_keV=ceiling_keV,
    )
    if spec is None:
        raise ValueError(f"No EnergySpectrum block in {path}")

    # Fill in identity / provenance
    spec.sample_id = sample_name
    spec.geometry = sample_location  # AtomSpectra uses Location loosely
    spec.comments = sample_note
    spec.detector_id = device_name
    spec.device_guid = device_guid
    spec.start_datetime = start_dt
    spec.end_datetime = end_dt
    spec.file_created_datetime = sample_time
    spec.background_link = background_link
    spec.filename_tokens = parse_filename(path)
    if spec.filename_tokens.get("is_background_hint"):
        spec.is_background = True
    spec.extras["multi_resultdata_in_file"] = multi_results
    spec.extras["weight_from_file"] = (
        _float(si.find("Weight")) if si is not None else None
    )
    spec.extras["volume_from_file"] = (
        _float(si.find("Volume")) if si is not None else None
    )

    # ----- Embedded background -----
    if parse_background:
        bes = rd.find("BackgroundEnergySpectrum")
        if bes is not None:
            bg_spec = _parse_energy_spectrum_block(
                bes, source_path=path, is_background=True,
                apply_energy_ceiling=apply_energy_ceiling,
                ceiling_keV=ceiling_keV,
            )
            if bg_spec is not None:
                # Inherit some identity from the host
                bg_spec.detector_id = device_name
                bg_spec.device_guid = device_guid
                bg_spec.sample_id = (
                    f"{sample_name} (embedded background)" if sample_name
                    else "embedded background"
                )
                bg_spec.extras["source_section"] = "BackgroundEnergySpectrum"
                spec.background_embedded = bg_spec

    # ----- Stored FWHM calibration (sample files only, usually) -----
    fwhm_node = rd.find("SimpleSqrtFwhmCalibration")
    spec.stored_fwhm_calibration = _parse_fwhm_calibration(fwhm_node)

    return spec
