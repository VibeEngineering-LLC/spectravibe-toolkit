"""
AtomSpectra mobile FORMAT: 3 (text) reader.

A line-oriented plain-text spectrum format produced by the AtomSpectra
Android application as an alternative to the desktop ResultDataFile XML
(handled by `atomspectra_xml.py`). One file = one ResultData; 8192
channels typical on the mobile detector.

Layout decoded from real captures (see #IO-1, 2026-06-27):

    line  1:  FORMAT: 3                            — magic / version banner
    line  2:  <date> <time> <tz> Counts: N, ~cps: X, Time: T s
              ^ human-readable summary; not parsed for primary fields,
                preserved verbatim in `extras["header_summary"]` for audit
    line  3:  1734112258532                        — start_datetime, Unix epoch ms
    line  4:  0                                    — service field, semantics unknown
    line  5:  0.0                                  — service field, semantics unknown
    line  6:  0.0                                  — service field, semantics unknown
              ^ lines 4-6 verbatim in extras["header_field_{4,5,6}"]
                pending vendor documentation
    line  7:  <user label>                         — sample_id / comments (Cyrillic OK)
    line  8:  <blank>
    line  9:  5739.300000                          — live time, seconds (float)
    line 10:  8192                                 — number of channels (int)
    line 11:  3                                    — polynomial degree (int)
    lines 12..(11 + degree + 1):
              calibration coefficients LOW-TO-HIGH — same convention as XML reader
              E(N) = c[0] + c[1]·N + c[2]·N² + ...
    lines (12 + degree)..(11 + degree + n_channels):
              channel counts (non-negative int)
    optional trailing blank line.

Field-level decisions:
  - real_time is not written by FORMAT: 3 → fallback `real_time = live_time`,
    flagged via `extras["real_time_source"] = "fallback_from_live_time"`.
    Mobile detector pipeline does not separately report dead-time; this
    matches the device behaviour.
  - Overflow-marker heuristic is reused verbatim from `atomspectra_xml`.
  - Energy-ceiling contract mirrors `read_atomspectra_xml`
    (apply_energy_ceiling=False default since BUG-9 / v1.18.32).
  - `source_format` = "atomspectra_txt".
"""

from __future__ import annotations

import os
import re
from datetime import datetime, timedelta, timezone
from typing import Optional

import numpy as np

from gamma.io.filename_hints import parse_filename
from gamma.spectrum import ENERGY_CEILING_KEV, Spectrum


_FORMAT_MAGIC = b"FORMAT: 3"


def looks_like_atomspectra_txt(head: bytes) -> bool:
    """Sniffer: True iff `FORMAT: 3` banner appears in the first 512 bytes."""
    return _FORMAT_MAGIC in head[:512]


def _parse_unix_epoch_ms(text: str) -> Optional[datetime]:
    """Parse Unix epoch milliseconds → timezone-aware UTC datetime."""
    text = text.strip()
    if not text:
        return None
    try:
        ms = int(text)
    except ValueError:
        return None
    try:
        return datetime.fromtimestamp(ms / 1000.0, tz=timezone.utc)
    except (OverflowError, OSError, ValueError):
        return None


def _parse_summary_tz_offset(line: str) -> Optional[timezone]:
    """Extract `±HHMM` offset from the line-2 summary → timezone object."""
    m = re.search(r"([+-])(\d{2})(\d{2})\b", line)
    if not m:
        return None
    sign = 1 if m.group(1) == "+" else -1
    hh = int(m.group(2))
    mm = int(m.group(3))
    return timezone(sign * timedelta(hours=hh, minutes=mm))


def _detect_and_drop_overflow(counts: np.ndarray) -> tuple:
    """Overflow heuristic mirrored from atomspectra_xml._detect_and_drop_overflow."""
    if len(counts) < 100:
        return counts, 0
    tail_start = max(1, int(0.99 * len(counts)))
    tail = counts[tail_start:-1]
    if len(tail) == 0:
        return counts, 0
    p95 = np.percentile(tail, 95) if tail.size else 0
    last = counts[-1]
    if last >= max(1000, 100 * max(p95, 1)):
        return counts[:-1].copy(), 1
    return counts, 0


def _channel_energies(coeffs: tuple, n: int) -> Optional[np.ndarray]:
    """Horner evaluation of E(channel) over [0, n). None if no usable cal."""
    if not coeffs or len(coeffs) < 2 or n == 0:
        return None
    channels = np.arange(n, dtype=np.float64)
    energies = np.zeros(n, dtype=np.float64)
    for c in reversed(coeffs):
        energies = energies * channels + c
    return energies


def _apply_energy_ceiling(counts: np.ndarray, coeffs: tuple,
                          ceiling_keV: float) -> tuple:
    """Drop channels above the ceiling. Mirrors atomspectra_xml helper."""
    energies = _channel_energies(coeffs, len(counts))
    if energies is None:
        return counts, 0, None
    n = len(counts)
    mask = energies <= ceiling_keV
    if mask.all():
        return counts, 0, float(energies[-1])
    first_above = int(np.argmax(~mask)) if (~mask).any() else n
    kept = counts[:first_above].copy()
    return (
        kept,
        n - first_above,
        float(energies[first_above - 1]) if first_above > 0 else None,
    )


def read_atomspectra_txt(
    path: str,
    *,
    apply_energy_ceiling: bool = False,
    ceiling_keV: Optional[float] = None,
) -> Spectrum:
    """
    Read an AtomSpectra mobile FORMAT: 3 text file. Returns a Spectrum.

    Args:
        path: filesystem path to the .txt file.
        apply_energy_ceiling: drop channels above the ceiling if True.
            Default False — mirrors `read_atomspectra_xml` post BUG-9
            (v1.18.32) so raw decoded range is kept by default.
        ceiling_keV: per-call override of `ENERGY_CEILING_KEV`.
            Ignored when `apply_energy_ceiling` is False.

    Raises:
        FileNotFoundError: path does not exist.
        ValueError: missing FORMAT: 3 banner, malformed metrology block,
            implausible polynomial degree, or truncated counts.
    """
    path = str(path)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"AtomSpectra FORMAT: 3 not found: {path}")

    with open(path, "r", encoding="utf-8-sig") as f:
        lines = [ln.rstrip("\r\n") for ln in f.readlines()]

    while lines and lines[-1] == "":
        lines.pop()

    if not lines or _FORMAT_MAGIC.decode("ascii") not in lines[0]:
        head = lines[0] if lines else ""
        raise ValueError(
            f"Not an AtomSpectra FORMAT: 3 file (first line {head!r}): {path}"
        )

    if len(lines) < 12:
        raise ValueError(
            f"FORMAT: 3 file too short ({len(lines)} lines): {path}"
        )

    summary_line = lines[1]
    start_dt = _parse_unix_epoch_ms(lines[2])
    tz_from_summary = _parse_summary_tz_offset(summary_line)
    if start_dt is not None and tz_from_summary is not None:
        start_dt = start_dt.astimezone(tz_from_summary)

    header_field_4 = lines[3]
    header_field_5 = lines[4]
    header_field_6 = lines[5]
    user_label = lines[6]
    # lines[7] is the blank separator; skipped intentionally.

    cursor = 8
    while cursor < len(lines) and lines[cursor] == "":
        cursor += 1

    try:
        live_time = float(lines[cursor]); cursor += 1
        n_channels_raw = int(lines[cursor]); cursor += 1
        degree = int(lines[cursor]); cursor += 1
    except (IndexError, ValueError) as exc:
        raise ValueError(
            f"Malformed FORMAT: 3 metrology block near line {cursor + 1} "
            f"of {path}: {exc}"
        ) from exc

    if degree < 0 or degree > 10:
        raise ValueError(
            f"Implausible polynomial degree {degree} in {path} (expected 0..10)"
        )

    n_coeffs = degree + 1
    if cursor + n_coeffs > len(lines):
        raise ValueError(
            f"FORMAT: 3 truncated in calibration block ({path}): "
            f"need {n_coeffs} coefficients, only {len(lines) - cursor} lines left"
        )

    coeffs = []
    for i in range(n_coeffs):
        try:
            coeffs.append(float(lines[cursor + i]))
        except ValueError as exc:
            raise ValueError(
                f"Bad calibration coefficient at line {cursor + i + 1} "
                f"of {path}: {exc}"
            ) from exc
    cursor += n_coeffs

    remaining = len(lines) - cursor
    if remaining < n_channels_raw:
        raise ValueError(
            f"FORMAT: 3 declares {n_channels_raw} channels but only "
            f"{remaining} count lines remain in {path}"
        )

    counts = np.empty(n_channels_raw, dtype=np.int64)
    for i in range(n_channels_raw):
        try:
            counts[i] = int(lines[cursor + i])
        except ValueError as exc:
            raise ValueError(
                f"Bad count at line {cursor + i + 1} of {path}: {exc}"
            ) from exc
    trailing_extra_lines = remaining - n_channels_raw

    counts, dropped_overflow = _detect_and_drop_overflow(counts)

    coeffs_tuple = tuple(coeffs)
    if apply_energy_ceiling:
        ceiling = ceiling_keV if ceiling_keV is not None else ENERGY_CEILING_KEV
        counts, dropped_high_e, e_max = _apply_energy_ceiling(
            counts, coeffs_tuple, ceiling
        )
    else:
        dropped_high_e = 0
        energies = _channel_energies(coeffs_tuple, len(counts))
        e_max = float(energies[-1]) if energies is not None else None

    real_time = live_time

    spec = Spectrum(
        counts=counts,
        live_time=live_time,
        real_time=real_time,
        source_path=path,
        source_format="atomspectra_txt",
        sample_id=user_label,
        comments=user_label,
        is_background=False,
        dropped_overflow_count=dropped_overflow,
        n_channels_raw=n_channels_raw,
        n_channels=len(counts),
        channel_pitch=1,
        energy_cal=coeffs_tuple if coeffs_tuple else None,
        energy_cal_degree=degree,
        energy_cal_source="stored" if coeffs_tuple else "",
        energy_max_keV_kept=e_max,
        start_datetime=start_dt,
    )

    spec.filename_tokens = parse_filename(path)
    if spec.filename_tokens.get("is_background_hint"):
        spec.is_background = True

    spec.extras["dropped_high_energy_count"] = dropped_high_e
    spec.extras["real_time_source"] = "fallback_from_live_time"
    spec.extras["format_version"] = 3
    spec.extras["header_summary"] = summary_line
    spec.extras["header_field_4"] = header_field_4
    spec.extras["header_field_5"] = header_field_5
    spec.extras["header_field_6"] = header_field_6
    spec.extras["trailing_extra_lines"] = trailing_extra_lines
    if tz_from_summary is not None:
        spec.extras["timezone_offset_seconds"] = int(
            tz_from_summary.utcoffset(None).total_seconds()
        )

    return spec
