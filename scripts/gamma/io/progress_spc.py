"""
Progress / Amplituda binary .spc reader.

Progress-family NaI/CsI gamma-spectrometer software (Amplituda, device types
`ScGammaSp` / `ScBetaSp` / `ScAlphaR`). Binary spectrum layout:

    bytes 0..3:    b"SF9I"                          - magic
    bytes 4..M:    cp1251 text header,
                   space-separated `KEY=VALUE` pairs; nested `{...}` sub-maps
                   for `CfgMap`, `GeomMap`, `ResultsMap`, `TaskMap`.
                   Padded with trailing spaces up to `M`.
    bytes M..M+4:  little-endian uint32 = binary_block_size (bytes) = SpLen*4
    bytes M+4..:   SpLen x uint32 LE                - channel counts

Header keys used (top-level):
    DeviceName, DeviceNick, DeviceType, Geometry,
    SpLen, LTime, RTime,
    KevPerCh, NullEnCh,
    StartTime, SavedToJournal, SourceFile,
    TaskName, ResList,
    Res137Cs, Err137Cs, Res226Ra, Err226Ra, Res232Th, Err232Th, Res40K, Err40K,
    Unit137Cs / Unit226Ra / Unit232Th / Unit40K,
    TaskMap (nested)                                - ProbeInfo, ProbeID, MassA, SepMethod, ...

Energy calibration is LINEAR in the source: `E(ch) = (ch - NullEnCh) * KevPerCh`.
Converted to LOW-to-HIGH polynomial form of `Spectrum.energy_cal`:
    a0 = -NullEnCh * KevPerCh
    a1 = KevPerCh

`StartTime` format is `DD.MM.YYYY HH:MM:SS` local wall clock (no timezone).
Local machine TZ is attached via `datetime.astimezone()` on read, mirroring the
BecqMoni writer TZ-attach policy (BUG-8 / Fix C3).

Sample mass in `TaskMap.MassA` is grams; converted to kilograms.

Device-computed activities and their sigmas are preserved verbatim in
`extras["progress_stored_results"]`.

`source_format` = "progress_spc".
"""

from __future__ import annotations

import os
import re
from datetime import datetime, timedelta
from typing import Optional

import numpy as np

from gamma.io.filename_hints import parse_filename
from gamma.spectrum import ENERGY_CEILING_KEV, Spectrum


_MAGIC = b"SF9I"
_KEY_RE = re.compile(r"([A-Za-z_][A-Za-z0-9_/()]*)=")


def looks_like_progress_spc(head: bytes) -> bool:
    """Sniffer: True iff file starts with `SF9I` magic bytes."""
    return head[:4] == _MAGIC


def _find_brace_blocks(text: str):
    """Return list of (start, end_exclusive) covering top-level `{...}`."""
    blocks = []
    i = 0
    n = len(text)
    while i < n:
        if text[i] == "{":
            depth = 1
            j = i + 1
            while j < n and depth > 0:
                if text[j] == "{":
                    depth += 1
                elif text[j] == "}":
                    depth -= 1
                j += 1
            blocks.append((i, j))
            i = j
        else:
            i += 1
    return blocks


def _mask_braces(text: str):
    """Replace top-level `{...}` blocks with placeholder tokens."""
    blocks = _find_brace_blocks(text)
    parts = []
    saved = []
    last = 0
    for k, (s, e) in enumerate(blocks):
        parts.append(text[last:s])
        parts.append(f"\x01{k}\x02")
        saved.append(text[s + 1 : e - 1])
        last = e
    parts.append(text[last:])
    return "".join(parts), saved


def _restore(s: str, saved) -> str:
    return re.sub(r"\x01(\d+)\x02", lambda m: "{" + saved[int(m.group(1))] + "}", s)


def _parse_header(text: str) -> dict:
    """
    Parse Progress header text into a dict. Nested `{...}` maps recurse.

    Values may contain spaces (e.g. `DeviceName=... NaI`, `Color=0 0 192`,
    `DeviceERegion=300 3000`); each value spans from `KEY=` to the next
    top-level `KEY=` marker or end-of-text. Nested `{...}` are masked before
    the top-level scan so their internal `=` do not confuse it.
    """
    masked, saved = _mask_braces(text)
    d = {}
    matches = list(_KEY_RE.finditer(masked))
    for idx, m in enumerate(matches):
        key = m.group(1)
        val_start = m.end()
        val_end = matches[idx + 1].start() if idx + 1 < len(matches) else len(masked)
        val = masked[val_start:val_end].rstrip()
        val = _restore(val, saved)
        val = val.strip()
        if val.startswith("{") and val.endswith("}"):
            d[key] = _parse_header(val[1:-1])
        else:
            d[key] = val
    return d


def _parse_progress_datetime(s: str) -> Optional[datetime]:
    if not s:
        return None
    try:
        return datetime.strptime(s.strip(), "%d.%m.%Y %H:%M:%S")
    except ValueError:
        return None


def _channel_energies(coeffs, n: int):
    if not coeffs or len(coeffs) < 2 or n == 0:
        return None
    channels = np.arange(n, dtype=np.float64)
    energies = np.zeros(n, dtype=np.float64)
    for c in reversed(coeffs):
        energies = energies * channels + c
    return energies


def _apply_energy_ceiling(counts, coeffs, ceiling_keV: float):
    energies = _channel_energies(coeffs, len(counts))
    if energies is None:
        return counts, 0, None
    n = len(counts)
    mask = energies <= ceiling_keV
    if mask.all():
        return counts, 0, float(energies[-1])
    first_above = int(np.argmax(~mask))
    kept = counts[:first_above].copy()
    return (
        kept,
        n - first_above,
        float(energies[first_above - 1]) if first_above > 0 else None,
    )


def read_progress_spc(
    path: str,
    *,
    apply_energy_ceiling: bool = False,
    ceiling_keV: Optional[float] = None,
) -> Spectrum:
    """
    Read a Progress / Amplituda binary .spc file. Returns a Spectrum.

    Args:
        path: filesystem path to the .spc file.
        apply_energy_ceiling: drop channels above the ceiling if True.
            Default False mirrors `read_atomspectra_xml` post BUG-9 (v1.18.32)
            so the raw decoded range is preserved by default.
        ceiling_keV: per-call override of `ENERGY_CEILING_KEV`. Ignored when
            `apply_energy_ceiling=False`.

    Raises:
        FileNotFoundError: path does not exist.
        ValueError: missing `SF9I` magic, missing `SpLen`, implausible size,
            or binary-block marker mismatch.
    """
    path = str(path)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Progress .spc not found: {path}")

    with open(path, "rb") as f:
        raw = f.read()

    if raw[:4] != _MAGIC:
        raise ValueError(f"Not a Progress SF9I .spc file (magic {raw[:4]!r}): {path}")

    text_all = raw.decode("cp1251", errors="replace")
    m = re.search(r"(?<![A-Za-z_])SpLen=(\d+)", text_all)
    if not m:
        raise ValueError(f"Progress .spc missing SpLen key in header: {path}")
    sp_len = int(m.group(1))
    if sp_len <= 0 or sp_len > (1 << 20):
        raise ValueError(f"Implausible SpLen={sp_len} in {path}")

    bin_size = sp_len * 4
    marker_pos = len(raw) - bin_size - 4
    if marker_pos < 4:
        raise ValueError(
            f"Progress .spc too short: SpLen={sp_len} implies bin_size={bin_size}, "
            f"file size {len(raw)}"
        )
    marker = int.from_bytes(raw[marker_pos : marker_pos + 4], "little")
    if marker != bin_size:
        raise ValueError(
            f"Progress .spc binary-block marker mismatch: header SpLen*4={bin_size}, "
            f"marker={marker} at offset {marker_pos} in {path}"
        )

    header_text = raw[4:marker_pos].decode("cp1251", errors="replace")
    fields = _parse_header(header_text)

    counts = np.frombuffer(raw[marker_pos + 4 :], dtype="<u4").astype(np.int64)
    if len(counts) != sp_len:
        raise ValueError(
            f"Progress .spc channel-count mismatch: SpLen={sp_len}, "
            f"binary carries {len(counts)} channels ({path})"
        )

    live_time = float(fields.get("LTime", 0.0) or 0.0)
    real_time = float(fields.get("RTime", live_time) or live_time)

    coeffs_tuple: Optional[tuple] = None
    energy_cal_source = ""
    try:
        kev_per_ch = float(fields["KevPerCh"])
        null_en_ch = float(fields["NullEnCh"])
        coeffs_tuple = (-null_en_ch * kev_per_ch, kev_per_ch)
        energy_cal_source = "stored"
    except (KeyError, ValueError):
        pass

    device_name = (fields.get("DeviceName") or "").strip()
    device_nick = (fields.get("DeviceNick") or "").strip()
    detector_id = device_name or device_nick
    geometry = (fields.get("Geometry") or "").strip()

    task_map = fields.get("TaskMap") or {}
    if not isinstance(task_map, dict):
        task_map = {}
    sample_id = (task_map.get("ProbeInfo") or "").strip()
    task_name = (fields.get("TaskName") or "").strip()

    sample_mass_kg = None
    mass_a = task_map.get("MassA")
    if mass_a:
        try:
            sample_mass_kg = float(mass_a) / 1000.0
        except ValueError:
            pass

    start_dt = _parse_progress_datetime(fields.get("StartTime", ""))
    if start_dt is not None:
        start_dt = start_dt.astimezone()
    saved_dt = _parse_progress_datetime(fields.get("SavedToJournal", ""))
    if saved_dt is not None:
        saved_dt = saved_dt.astimezone()

    end_dt = None
    if start_dt is not None and real_time > 0:
        end_dt = start_dt + timedelta(seconds=real_time)

    dropped_high_e = 0
    e_max: Optional[float] = None
    if apply_energy_ceiling and coeffs_tuple is not None:
        ceiling = ceiling_keV if ceiling_keV is not None else ENERGY_CEILING_KEV
        counts, dropped_high_e, e_max = _apply_energy_ceiling(counts, coeffs_tuple, ceiling)
    else:
        energies = _channel_energies(coeffs_tuple, len(counts)) if coeffs_tuple else None
        e_max = float(energies[-1]) if energies is not None else None

    spec = Spectrum(
        counts=counts,
        live_time=live_time,
        real_time=real_time,
        source_path=path,
        source_format="progress_spc",
        sample_id=sample_id,
        geometry=geometry,
        detector_id=detector_id,
        comments=task_name,
        is_background=False,
        sample_mass_kg=sample_mass_kg,
        n_channels_raw=sp_len,
        n_channels=len(counts),
        channel_pitch=1,
        energy_cal=coeffs_tuple,
        energy_cal_degree=(len(coeffs_tuple) - 1) if coeffs_tuple else None,
        energy_cal_source=energy_cal_source,
        energy_max_keV_kept=e_max,
        start_datetime=start_dt,
        end_datetime=end_dt,
        file_created_datetime=saved_dt,
    )

    spec.filename_tokens = parse_filename(path)
    if spec.filename_tokens.get("is_background_hint"):
        spec.is_background = True

    stored_results = {}
    for nuc in ("137Cs", "226Ra", "232Th", "40K"):
        val = fields.get(f"Res{nuc}")
        err = fields.get(f"Err{nuc}")
        unit = fields.get(f"Unit{nuc}")
        if val is None:
            continue
        try:
            stored_results[nuc] = {
                "activity": float(val),
                "sigma": float(err) if err else None,
                "unit": unit,
            }
        except ValueError:
            pass
    if stored_results:
        spec.extras["progress_stored_results"] = stored_results

    spec.extras["dropped_high_energy_count"] = dropped_high_e
    spec.extras["dropped_overflow_count"] = 0
    spec.extras["progress_device_nick"] = device_nick
    spec.extras["progress_device_type"] = (fields.get("DeviceType") or "").strip()
    spec.extras["progress_task_name"] = task_name
    if fields.get("SourceFile"):
        spec.extras["progress_source_file"] = fields["SourceFile"].strip()
    for tm_key in ("ProbeID", "ProbeType", "ProbeMethod", "SepMethod", "SepComponents"):
        val = task_map.get(tm_key)
        if val:
            spec.extras[f"progress_task_{tm_key.lower()}"] = val

    return spec