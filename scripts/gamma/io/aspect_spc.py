"""
ASPECT (АСПЕКТ / Dubna) binary `.spc` reader.

ASPECT — Russian gamma-spectroscopy vendor (JINR Dubna) that ships handheld and
lab scintillator instruments with proprietary Windows-side software. Binary
`.spc` layout follows the `tSpectrHeader` C struct in the vendor's `spheader.h`
(SSOT: Russian original `Формат файла спектра.h`, cp1251-encoded; English
translation kept alongside as ChatGPT interpretation of the Russian header).

    bytes   0..512:  fixed-size `tSpectrHeader` (`#pragma pack(1)`)
    bytes 512..EOF:  buffer × uint32 LE       - channel counts

The 512-byte header layout (verbatim from `spheader.h`):

    offset  size  field           notes
     0      16    device[16]      device name (ASCII, NUL-padded)
    16       1    unit            ADC number
    17       1    section         section number
    18       6    sampleTime[6]   sample collection D,M,Y,H,M,S (year = 2000+Y)
    24      16    weight[16]      mass in kg, ASCII float
    40      16    volume[16]      volume in L, ASCII float
    56       6    startTime[6]    spectrum acquisition D,M,Y,H,M,S
    62       4    liveTime        DWORD LE, MILLISECONDS
    66       4    realTime        DWORD LE, MILLISECONDS
    70       2    buffer          WORD LE, ADC buffer size (== channel count N)
    72       2    gain            WORD LE, ADC resolution
    74       2    offset          WORD LE, ADC offset
    76       2    lowerLevel      WORD LE, LLD
    78       2    upperLevel      WORD LE, ULD
    80       2    first           WORD LE, first channel
    82       2    last            WORD LE, last channel (== N-1 or 0)
    84       2    chWidth         WORD LE, channel width
    86      20    adcCont[10]     10 × WORD LE, extra ADC params
   106       1    prepar          preparation type (0 = none, 1 = ashing)
   107       6    preparTime[6]   preparation D,M,Y,H,M,S
   113      15    reserve[15]     reserved
   128      64    energy[4][16]   4 × 16-byte ASCII floats (a0,a1,a2,a3)
   192      64    fwhm[4][16]     4 × 16-byte ASCII floats (FWHM coefficients)
   256     256    comment[4][64]  4 × 64-byte ASCII lines

Time units — liveTime/realTime are **milliseconds** (empirical: 5 samples give
plausible 22-124 min ranges for scintillator lab measurements; if seconds, the
values would translate to years). Reader divides by 1000 to normalize to
seconds, mirroring `Spectrum.live_time` / `real_time` contract.

Energy calibration — ASCII floats in bytes 128..191 (four 16-byte slots). When
all four slots are empty/whitespace/zero, calibration is missing and the reader
emits a warning in `extras["aspect_uncalibrated"]`. Cs-137 anchor at ch 74/77 in
the test corpus gives a rough 8.6-8.9 keV/ch bootstrap for 25 mm² Ø25 mm probes,
but that is a downstream detector-specific decision, not the reader's business.

`source_format` = "aspect_spc".

Provenance:
- Vendor header: `references/vendor/aspect/spheader_ru_cp1251.h` (SSOT, Russian).
- English translation: `references/vendor/aspect/spheader_en.h` (ChatGPT).
- Sample corpus: `references/vendor/aspect/aspect_samples_issue47.zip` (SpecUtils
  issue #47 attachment, public).
- SpecUtils integration reference: https://github.com/sandialabs/SpecUtils/issues/47
  (their `src/SpecFile_aspect.cpp` implements the same struct).
"""

from __future__ import annotations

import os
import struct
from datetime import datetime
from typing import Optional

import numpy as np

from gamma.io.filename_hints import parse_filename
from gamma.spectrum import Spectrum


_HEADER_SIZE = 512
_VALID_N_CHANNELS = frozenset({1024, 2048, 4096, 8192, 16384})


def _decode_ascii(raw: bytes) -> str:
    """Decode NUL/space-padded ASCII slot, tolerating cp1251 in comments."""
    stripped = raw.split(b"\x00", 1)[0].rstrip(b" \t\r\n\x00")
    if not stripped:
        return ""
    try:
        return stripped.decode("ascii").strip()
    except UnicodeDecodeError:
        return stripped.decode("cp1251", errors="replace").strip()


def _parse_ascii_float(raw: bytes) -> Optional[float]:
    """Parse a 16-byte ASCII float slot. Returns None on empty/invalid."""
    s = _decode_ascii(raw)
    if not s:
        return None
    try:
        return float(s.replace(",", "."))
    except ValueError:
        return None


def _parse_dmyhms(raw: bytes) -> Optional[datetime]:
    """Parse the 6-byte DMYHMS timestamp. Year is Y + 2000."""
    if len(raw) < 6:
        return None
    d, m, y, hh, mm, ss = raw[0], raw[1], raw[2], raw[3], raw[4], raw[5]
    if d == 0 and m == 0 and y == 0 and hh == 0 and mm == 0 and ss == 0:
        return None
    try:
        return datetime(2000 + y, m, d, hh, mm, ss)
    except ValueError:
        return None


def _dmyhms_looks_valid(raw: bytes) -> bool:
    """Range-check DMYHMS without raising."""
    if len(raw) < 6:
        return False
    d, m, y, hh, mm, ss = raw[0], raw[1], raw[2], raw[3], raw[4], raw[5]
    return (
        1 <= d <= 31
        and 1 <= m <= 12
        and 0 <= y <= 99
        and 0 <= hh <= 23
        and 0 <= mm <= 59
        and 0 <= ss <= 59
    )


def looks_like_aspect_spc(head: bytes) -> bool:
    """
    Sniffer for ASPECT `.spc`.

    The format has no magic bytes, so we validate the 512-byte header
    structurally:
      - buffer (u16 LE @ offset 70) ∈ {1024, 2048, 4096, 8192, 16384}
      - last (u16 LE @ offset 82) == buffer - 1  OR  last == 0
      - startTime[6] (bytes 56..62) is a valid DMYHMS tuple
      - first (u16 LE @ offset 80) == 0
      - gain (u16 LE @ offset 72) plausible (≤ 65535, i.e. any u16)

    Progress `.spc` (SF9I magic at bytes 0..3) is rejected via the `buffer`
    check — SF9I magic decodes to buffer = 0x3946 = 14662, which is not in
    the valid set.
    """
    if len(head) < _HEADER_SIZE:
        return False
    buffer_val = int.from_bytes(head[70:72], "little")
    if buffer_val not in _VALID_N_CHANNELS:
        return False
    last_val = int.from_bytes(head[82:84], "little")
    if last_val != buffer_val - 1 and last_val != 0:
        return False
    first_val = int.from_bytes(head[80:82], "little")
    if first_val != 0:
        return False
    if not _dmyhms_looks_valid(head[56:62]):
        return False
    return True


def _read_energy_coeffs(header: bytes) -> Optional[tuple]:
    """
    Parse energy[4][16] block at offset 128. Returns LOW-to-HIGH tuple, or
    None if all slots are empty (uncalibrated file).
    """
    coeffs = []
    all_zero = True
    for i in range(4):
        raw = header[128 + i * 16 : 128 + (i + 1) * 16]
        val = _parse_ascii_float(raw)
        if val is None:
            coeffs.append(0.0)
        else:
            coeffs.append(val)
            if val != 0.0:
                all_zero = False
    if all_zero:
        return None
    while len(coeffs) > 1 and coeffs[-1] == 0.0:
        coeffs.pop()
    return tuple(coeffs)


def _read_fwhm_coeffs(header: bytes) -> Optional[tuple]:
    """Parse fwhm[4][16] block at offset 192. Same layout as energy block."""
    coeffs = []
    all_zero = True
    for i in range(4):
        raw = header[192 + i * 16 : 192 + (i + 1) * 16]
        val = _parse_ascii_float(raw)
        if val is None:
            coeffs.append(0.0)
        else:
            coeffs.append(val)
            if val != 0.0:
                all_zero = False
    if all_zero:
        return None
    while len(coeffs) > 1 and coeffs[-1] == 0.0:
        coeffs.pop()
    return tuple(coeffs)


def _read_comments(header: bytes) -> list[str]:
    """Parse comment[4][64] block at offset 256."""
    lines = []
    for i in range(4):
        raw = header[256 + i * 64 : 256 + (i + 1) * 64]
        s = _decode_ascii(raw)
        if s:
            lines.append(s)
    return lines


def read_aspect_spc(path: str) -> Spectrum:
    """
    Read an ASPECT `.spc` file. Returns a Spectrum.

    Raises:
        FileNotFoundError: path does not exist.
        ValueError: file too short, buffer field implausible, or size
            inconsistent with `512 + buffer * 4`.
    """
    path = str(path)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"ASPECT .spc not found: {path}")

    with open(path, "rb") as f:
        raw = f.read()

    if len(raw) < _HEADER_SIZE + 4:
        raise ValueError(
            f"ASPECT .spc too short: {len(raw)} bytes (need ≥ 516): {path}"
        )

    header = raw[:_HEADER_SIZE]

    # Fixed-layout fields via struct
    (
        device_raw,
        unit,
        section,
        sample_time_raw,
        weight_raw,
        volume_raw,
        start_time_raw,
        live_time_ms,
        real_time_ms,
        buffer_val,
        gain,
        offset_adc,
        lower_level,
        upper_level,
        first_ch,
        last_ch,
        ch_width,
    ) = struct.unpack("<16sBB6s16s16s6sIIHHHHHHHH", header[:86])

    if buffer_val not in _VALID_N_CHANNELS:
        raise ValueError(
            f"ASPECT .spc: implausible buffer={buffer_val} (expected one of "
            f"{sorted(_VALID_N_CHANNELS)}): {path}"
        )

    expected_size = _HEADER_SIZE + buffer_val * 4
    if len(raw) != expected_size:
        raise ValueError(
            f"ASPECT .spc size mismatch: buffer={buffer_val} implies "
            f"{expected_size} bytes, got {len(raw)}: {path}"
        )

    if last_ch != buffer_val - 1 and last_ch != 0:
        raise ValueError(
            f"ASPECT .spc: last={last_ch} inconsistent with buffer={buffer_val} "
            f"(expected {buffer_val - 1} or 0): {path}"
        )

    counts = np.frombuffer(raw[_HEADER_SIZE:], dtype="<u4").astype(np.int64)
    assert len(counts) == buffer_val, "size check passed but frombuffer length mismatch"

    live_time_s = live_time_ms / 1000.0
    real_time_s = real_time_ms / 1000.0

    energy_coeffs = _read_energy_coeffs(header)
    fwhm_coeffs = _read_fwhm_coeffs(header)
    energy_cal_source = "stored" if energy_coeffs is not None else ""

    weight_kg = _parse_ascii_float(weight_raw)
    volume_l = _parse_ascii_float(volume_raw)

    device_name = _decode_ascii(device_raw)
    comments_list = _read_comments(header)
    comments_text = " | ".join(comments_list) if comments_list else ""

    start_dt = _parse_dmyhms(start_time_raw)
    sample_dt = _parse_dmyhms(sample_time_raw)

    end_dt = None
    if start_dt is not None and real_time_s > 0:
        from datetime import timedelta
        end_dt = start_dt + timedelta(seconds=real_time_s)

    energy_max_keV_kept: Optional[float] = None
    if energy_coeffs is not None:
        channels = np.arange(len(counts), dtype=np.float64)
        energies = np.zeros_like(channels)
        for c in reversed(energy_coeffs):
            energies = energies * channels + c
        energy_max_keV_kept = float(energies[-1])

    spec = Spectrum(
        counts=counts,
        live_time=live_time_s,
        real_time=real_time_s,
        source_path=path,
        source_format="aspect_spc",
        sample_id="",
        geometry="",
        detector_id=device_name,
        comments=comments_text,
        is_background=False,
        sample_mass_kg=weight_kg,
        sample_volume_ml=(volume_l * 1000.0) if volume_l is not None else None,
        n_channels_raw=buffer_val,
        n_channels=len(counts),
        channel_pitch=1,
        energy_cal=energy_coeffs,
        energy_cal_degree=(len(energy_coeffs) - 1) if energy_coeffs else None,
        energy_cal_source=energy_cal_source,
        energy_max_keV_kept=energy_max_keV_kept,
        start_datetime=start_dt,
    )
    if end_dt is not None:
        spec.end_datetime = end_dt

    spec.filename_tokens = parse_filename(path)
    if spec.filename_tokens.get("is_background_hint"):
        spec.is_background = True

    spec.extras["aspect_unit"] = int(unit)
    spec.extras["aspect_section"] = int(section)
    spec.extras["aspect_buffer"] = int(buffer_val)
    spec.extras["aspect_gain"] = int(gain)
    spec.extras["aspect_offset"] = int(offset_adc)
    spec.extras["aspect_lower_level"] = int(lower_level)
    spec.extras["aspect_upper_level"] = int(upper_level)
    spec.extras["aspect_first"] = int(first_ch)
    spec.extras["aspect_last"] = int(last_ch)
    spec.extras["aspect_ch_width"] = int(ch_width)
    spec.extras["aspect_live_time_ms"] = int(live_time_ms)
    spec.extras["aspect_real_time_ms"] = int(real_time_ms)
    if sample_dt is not None:
        spec.extras["aspect_sample_collection_datetime"] = sample_dt.isoformat()
    if comments_list:
        spec.extras["aspect_comments_lines"] = comments_list
    if energy_coeffs is None:
        spec.extras["aspect_uncalibrated"] = True
    if fwhm_coeffs is not None:
        spec.extras["aspect_fwhm_coeffs"] = list(fwhm_coeffs)

    return spec
