"""
Lsrm SpectraLine `.spe` file reader.

Format: a CP-1251 text header of `KEY=VALUE\\r\\n` lines, ending with a
literal `SPECTR=` marker that is followed immediately by a raw binary
block of channel counts as little-endian unsigned 32-bit integers. The
total file consists of one continuous header–data stream with no
delimiters between the marker and the binary block.

Verified empirically against:
  - M_cs_легкий_2001-2005.spe  (Cs-137 control source, 1023 ch, NaI Гамма-1С)
  - M_cs_тяж_2001-2005.spe     (Cs-137, heavy density geometry)
  - M_k_легкий_2001-2005.spe   (K-40 control source)
  - M_ra_легкий_2001-2007.spe  (Ra-226 control source)
  - M_th_легкий_2001-2005.spe  (Th-232 control source)

Header keys observed (all five fixtures):
  SHIFR=…                — sample identifier
  NOMER=…                — sample number
  TYPE=…                 — measurement type (Калибровка = "Calibration")
  CONFIGNAME=…           — instrument configuration (e.g. "Gamma-1S")
  MEASBEGIN=DD-MM-YY HH:MM:SS  — measurement start date/time
  PREPBEGIN=DD-MM-YY     — sample preparation start
  PREPEND=DD-MM-YY       — sample preparation end
  TLIVE=NNNN.NN          — live time, seconds (float)
  TREAL=NNNN.NN          — real time, seconds (float)
  OPERATOR=…
  GEOMETRY=…             — e.g. "Маринелли"
  DETECTOR=…             — e.g. "Гамма-1С" (NaI 63×63)
  SETTYPE=…              — set type
  CONTTYPE=…             — container type
  MATERIAL=…             — container material
  DISTANCE=NN.N          — source-detector distance, cm
  DETRADIUS=NN.N         — detector radius, cm
  RAWMASS=val;err
  PROBEMASS=val;err
  SAMPLEMASS=val;err
  RAWVOLUME=val;err
  PROBEVOLUME=val;err
  SAMPLEVOLUME=val;err
  ENBOUNDS=lo,hi         — energy range in use
  ENERGYLOW=...          — low-energy region settings (often zeros)
  ENERGY=N,a0,a1,a2,a3,…  — N is degree-marker (e.g. 3), then 7
                           polynomial coefficients low-to-high
  FWHM=N,c0,c1,c2,…      — N is degree-marker (typically 2-3), then up
                           to 7 polynomial coefficients. BUG-22 /
                           2026-06-02: the polynomial argument is
                           z = √E_keV, NOT E directly, despite the
                           internal model label `lsrm_fwhm_polynomial_in_E`
                           and prior code comments here. Confirmed by
                           LSRM «Алгоритмические основы» §8.3
                           («Калибровка по полуширине»):
                             FWHM_keV(E) = Σ_k c_k · z^k,  z = √E_keV
                           Verified against the Th-232 Marinelli
                           archive fixture (coefs ≈ (8.97, −0.898,
                           0.143, −0.00169)): evaluated as a
                           polynomial in √E it yields physically
                           reasonable values (≈52 keV at 661 keV →
                           7.8% NaI, ≈111 keV at 2614 keV → 4.3%); a
                           naïve polynomial-in-E evaluation produces
                           negative FWHM (−14 888 keV at 238 keV)
                           and is the BUG-22 root cause. The stored
                           coefficients are returned as-is on
                           `StoredFwhmCalibration.coefficients`;
                           downstream evaluation must use the √E
                           argument. (Convention also differs from
                           AtomSpectra `SimpleSqrtFwhm`, where the
                           model is FWHM²(N)=c₀+c₁·N in channels.)
  FWHMCALIBRATIONFILE=…  — provenance file path(s)
  ENERGYCALIBRATIONFILE=… — provenance file path
  COMMENT=…              — free-text comment
  PEAKS=N                — number of peak entries in the embedded peak
                           table (one line per peak, tab-separated)
  …peak rows…
  SPECTR=<binary>        — binary block immediately follows the `=`

Token economy: returns a `Spectrum` dataclass; counts array on it but
nothing else heavy. Filename hints are pushed through the existing
`gamma.io.filename_hints.parse_filename`.
"""

from __future__ import annotations

import struct
from datetime import datetime
from pathlib import Path
from typing import Iterable, Optional

import numpy as np

from gamma.spectrum import (
    Spectrum,
    StoredFwhmCalibration,
    FwhmCalPeak,
    ENERGY_CEILING_KEV,
)
from gamma.io.filename_hints import parse_filename


# ============================================================================
# Public entry point
# ============================================================================

def read_lsrm_spe(
    path: str,
    *,
    apply_energy_ceiling: bool = False,
    ceiling_keV: Optional[float] = None,
) -> Spectrum:
    """
    Read an Lsrm SpectraLine .spe file.

    Args:
        path: filesystem path to the .spe file.
        apply_energy_ceiling: if True, drop trailing channels whose
            energy under the stored calibration exceeds the ceiling.
            **Default False since BUG-9 (v1.18.32, 2026-06-03):** silently
            trimming a 1024-channel LSRM file to 1003 channels when a0≈-8
            and gain≈3 keV/ch surprised users; the reader now keeps every
            decoded channel by default. Set True explicitly at the call
            site if you need the project-scope 3 MeV trim. See
            `gamma.spectrum.trim_to_working_energy` for a post-read trim
            that is easier to audit than this reader-stage knob.
        ceiling_keV: per-call override of `ENERGY_CEILING_KEV` (3000 keV
            by project scope). When None, the module constant is used.
            Ignored if `apply_energy_ceiling` is False.

    Returns:
        Populated `Spectrum`.
    """
    p = Path(path)
    raw = p.read_bytes()

    # ----- locate SPECTR= and split header from binary -----
    marker = b"SPECTR="
    mi = raw.find(marker)
    if mi < 0:
        raise ValueError(
            f"Lsrm .spe missing required SPECTR= marker: {path!r}"
        )
    header_bytes = raw[:mi]
    binary_bytes = raw[mi + len(marker):]

    # ----- parse header (CP-1251 KEY=VALUE\r\n) -----
    header = header_bytes.decode("cp1251", errors="replace")
    fields = _parse_kv_lines(header)

    # ----- parse binary counts (uint32 LE) -----
    n_full = len(binary_bytes) // 4
    if n_full == 0:
        raise ValueError(
            f"Lsrm .spe data block empty or unaligned: {path!r}"
        )
    counts = np.frombuffer(
        binary_bytes[: n_full * 4], dtype="<u4"
    ).astype(np.int64)

    # Build Spectrum dataclass
    spec = Spectrum(
        counts=counts,
        live_time=_get_float(fields, "TLIVE", default=0.0),
        real_time=_get_float(fields, "TREAL", default=0.0),
        source_path=str(p),
        source_format="lsrm_spe",
    )

    # --- identity / metadata ---
    spec.sample_id = fields.get("SHIFR", "")
    spec.operator = fields.get("OPERATOR", "")
    spec.geometry = fields.get("GEOMETRY", "")
    spec.detector_id = fields.get("DETECTOR", "")
    spec.comments = fields.get("COMMENT", "")

    # F-330 / v1.18.18.4 — auto-extract passport activities из COMMENT.
    # LSRM конвенция: «Nuclide - X Бк/кг (Y%) [от DD.MM.YY]». Stashing
    # parsed entries в extras["lsrm_passport"] позволяет wrapper.py
    # auto-передать их в `passport_activity_Bq` (F-326) без ручного
    # ввода пользователем. Decay correction применяется в wrapper при
    # наличии reference_date.
    if spec.comments:
        try:
            from gamma.io.lsrm_passport import parse_lsrm_passport_comment
            _entries = parse_lsrm_passport_comment(spec.comments)
            if _entries:
                spec.extras["lsrm_passport"] = [
                    {
                        "nuclide": e.nuclide,
                        "value": float(e.value),
                        "unit": e.unit,
                        "uncertainty_pct": float(e.uncertainty_pct),
                        "reference_date": (
                            e.reference_date.isoformat()
                            if e.reference_date else None
                        ),
                        "is_specific_activity": e.is_specific_activity,
                        "raw_match": e.raw_match,
                    }
                    for e in _entries
                ]
        except Exception:
            # Defensive: never let passport parsing break .spe reading.
            pass

    # TYPE=Калибровка → mark as calibration source via extras; not is_background
    type_val = fields.get("TYPE", "")
    if type_val:
        spec.extras["lsrm_type"] = type_val
    if "CONFIGNAME" in fields:
        spec.extras["lsrm_config"] = fields["CONFIGNAME"]

    # --- datetime ---
    if "MEASBEGIN" in fields:
        spec.start_datetime = _parse_lsrm_datetime(fields["MEASBEGIN"])

    # --- masses / volumes (carry through as extras; skill scope excludes
    #     activity computations so we don't act on these) ---
    for key in ("RAWMASS", "PROBEMASS", "SAMPLEMASS",
                "RAWVOLUME", "PROBEVOLUME", "SAMPLEVOLUME",
                "DISTANCE", "DETRADIUS", "MATERIAL"):
        if key in fields:
            spec.extras[f"lsrm_{key.lower()}"] = fields[key]

    # --- F-130 / v1.17.7: автоматическое определение плотности образца ---
    # Источник приоритета (в порядке убывания достоверности):
    #   1. MATERIAL JSON → поле "Ro" (г/см³) — прямая запись плотности
    #   2. SAMPLEMASS / SAMPLEVOLUME → ρ = m/V
    #   3. PROBEMASS  / PROBEVOLUME  → fallback
    # Результат: extras["lsrm_sample_density_g_cm3"] (float).
    # Используется downstream (F-122 self-attenuation), если CLI флаг
    # `--sample-density-g-cm3` не задан явно.
    auto_density = _auto_extract_density(fields)
    if auto_density is not None:
        spec.extras["lsrm_sample_density_g_cm3"] = float(auto_density[0])
        spec.extras["lsrm_density_source"] = auto_density[1]

    # --- F-140 / v1.17.7: автоматическое определение массы образца (кг).
    # SAMPLEMASS / PROBEMASS поля LSRM хранят массу в граммах. Pipeline
    # принимает массу в килограммах через kwarg `sample_mass_kg`. F-140
    # автоматически извлекает её, если CLI флаг `--sample-mass-kg` не
    # задан явно. Источник приоритета:
    #   1. SAMPLEMASS (г) → /1000 → кг + source="sample_mass_field"
    #   2. PROBEMASS  (г) → /1000 → кг + source="probe_mass_field"
    auto_mass = _auto_extract_mass_kg(fields)
    if auto_mass is not None:
        spec.extras["lsrm_sample_mass_kg"] = float(auto_mass[0])
        spec.extras["lsrm_mass_source"] = auto_mass[1]

    # --- BUG-1 / 2026-06-02: typed surface for SAMPLEMASS / SAMPLEVOLUME ---
    # LSRM stores SAMPLEMASS=<grams>;<grams_uncertainty>. Pipeline (run_skill)
    # historically reached into spec.extras["lsrm_sample_mass_kg"] which was
    # set only for SAMPLEMASS/PROBEMASS in (1 g .. 100 kg) range. A typed
    # Spectrum.sample_mass_kg / sample_mass_uncertainty_kg pair removes the
    # need for downstream callers to know the LSRM-specific extras key, and
    # carries the uncertainty (previously dropped). SAMPLEVOLUME is exposed
    # the same way for completeness (LSRM stores it in millilitres).
    pair = _parse_value_err_pair_full(fields.get("SAMPLEMASS", ""))
    if pair is not None:
        m_g, u_g = pair
        # Sanitary range mirrors `_auto_extract_mass_kg`: 1 g .. 100 kg.
        if 1.0 <= m_g <= 100_000.0:
            spec.sample_mass_kg = float(m_g) / 1000.0
            if u_g is not None and u_g >= 0.0:
                spec.sample_mass_uncertainty_kg = float(u_g) / 1000.0
    vpair = _parse_value_err_pair_full(fields.get("SAMPLEVOLUME", ""))
    if vpair is not None:
        v_ml, vu_ml = vpair
        # LSRM SAMPLEVOLUME is in millilitres. Sanity: 0.1 ml .. 100 L.
        if 0.1 <= v_ml <= 100_000.0:
            spec.sample_volume_ml = float(v_ml)
            if vu_ml is not None and vu_ml >= 0.0:
                spec.sample_volume_uncertainty_ml = float(vu_ml)

    # --- energy calibration ---
    # Lsrm files can carry energy calibration in TWO forms:
    #   1. ENERGY=N,a0,a1,a2,... — ordinary polynomial (preferred when present)
    #   2. ENERGY_ZONESCOUNT/ENERGY_ZONE_<z>/ENERGY_CURVE_<z>_<k> — orthogonal
    #      polynomial decomposition (used for piecewise calibration over zones)
    # In all observed files where both are present, the ordinary polynomial is
    # the evaluated form of the orthogonal expansion. We prefer the ordinary
    # polynomial when it is present, and fall back to orthopoly only if it is
    # entirely absent (currently never observed in our fixtures).
    e_line = fields.get("ENERGY")
    if e_line:
        coefs = _parse_polynomial_line(e_line)
        if coefs:
            # Trim trailing zeros for a clean degree value
            while len(coefs) > 1 and coefs[-1] == 0.0:
                coefs.pop()
            spec.energy_cal = tuple(coefs)
            spec.energy_cal_degree = max(0, len(coefs) - 1)
            spec.energy_cal_source = "stored"
    elif "ENERGY_ZONESCOUNT" in fields:
        # Fallback: try to extract a polynomial from the orthopoly expansion.
        # For single-zone calibration (the only case observed), we can read
        # the highest-order polynomial coefficients from ENERGY_CURVE_<z>_<k>
        # — these are the coefficients of P_k(channel), the k-th orthogonal
        # polynomial. The final E(ch) requires combining them with the
        # weights in ENERGY_CURVE_<z>. We don't implement this fully here
        # because the ordinary polynomial form is always provided in practice;
        # for diagnostic purposes we report what we found.
        spec.extras["lsrm_energy_orthopoly_only"] = True
        spec.extras["lsrm_energy_zones"] = fields.get("ENERGY_ZONESCOUNT")

    # Record energy calibration quality if present
    if "ENERGY_QUALITY" in fields:
        # Format: chi2, integral_nonlinearity, calibration_nonlinearity
        qual_parts = [_safe_float(p) for p in fields["ENERGY_QUALITY"].split(",")]
        qual_parts = [p for p in qual_parts if p is not None]
        if qual_parts:
            spec.extras["lsrm_energy_quality"] = qual_parts

    # --- FWHM calibration ---
    fwhm_line = fields.get("FWHM")
    if fwhm_line:
        fw_coefs = _parse_polynomial_line(fwhm_line)
        if fw_coefs:
            while len(fw_coefs) > 1 and fw_coefs[-1] == 0.0:
                fw_coefs.pop()
            spec.stored_fwhm_calibration = StoredFwhmCalibration(
                calibration_peaks=[],
                coefficients=tuple(fw_coefs),
                model="lsrm_fwhm_polynomial_in_E",
            )
    elif "FWHM_ORT" in fields:
        # Orthogonal polynomial FWHM model — JSON-encoded.
        # Format: {"Polynomials": [{"Power":N, "Vector":[...], "Matrix":[[...]],
        #                            "LeftBound":..., "RightBound":..., ...}],
        #          "SewType": int}
        # Vector + Matrix together encode an orthogonal polynomial expansion.
        # We store the raw JSON for forensic use; full evaluation requires
        # reconstructing the orthogonal basis (deferred).
        spec.extras["lsrm_fwhm_orthopoly_json"] = fields["FWHM_ORT"]

    # --- SPECTRSIZE field (used in .spm/.spex multi-spectrum files; here
    #     informational since the binary block size already tells us) ---
    if "SPECTRSIZE" in fields:
        try:
            spec.extras["lsrm_spectr_size_field"] = int(fields["SPECTRSIZE"])
        except ValueError:
            pass

    # --- embedded peak table (informational; we re-find peaks ourselves) ---
    # Per Lsrm spec §7.5.2.1: each row is space-separated:
    #   position d_pos energy d_energy fwhm d_fwhm area d_area chi2 \
    #     left right (then optional ';'-separated nuclide list)
    # We store both the raw text for forensics and a parsed list for
    # downstream use (Phase 2.1b multiplet deconvolution can seed from
    # these pre-found peaks).
    peaks_n = _get_int(fields, "PEAKS", default=0)
    if peaks_n and "PEAKS_TABLE" in fields:
        spec.extras["lsrm_peaks_table_raw"] = fields["PEAKS_TABLE"]
        parsed_peaks = _parse_peaks_table(fields["PEAKS_TABLE"])
        if parsed_peaks:
            spec.extras["lsrm_peaks_table"] = parsed_peaks

    # --- ZONES section (pre-marked multiplet zones — useful for Phase 2.1b) ---
    # Per Lsrm spec §7.5.2.1: each row gives left, right, n_peaks_in_zone,
    # minimize parameters, polynomial degree
    zones_n = _get_int(fields, "ZONES", default=0)
    if zones_n and "ZONES_TABLE" in fields:
        spec.extras["lsrm_zones_table_raw"] = fields["ZONES_TABLE"]
        parsed_zones = _parse_zones_table(fields["ZONES_TABLE"])
        if parsed_zones:
            spec.extras["lsrm_zones_table"] = parsed_zones

    # --- record raw channel count BEFORE trimming ---
    spec.n_channels_raw = int(n_full)
    spec.channel_pitch = 1

    # --- BUG-9 / 2026-06-03: calibration drift diagnostic ---------------------
    # When the stored energy calibration has a0 < 0, the first few channels
    # map to E < 0. Per user clarification, this is **not a data bug** — it is
    # a diagnostic signature of a spectrometer that has drifted left of its
    # original calibration anchor (zero-channel pedestal shifted into negative
    # energy). The channels and their counts MUST be preserved: dropping them
    # would erase the very evidence of the drift.
    #
    # This block:
    #   1. Detects the condition (a0 < 0).
    #   2. Counts how many leading channels currently map to E < 0.
    #   3. Stashes both as diagnostic metadata in `extras`.
    #   4. DOES NOT trim or alter `spec.counts`. The reader already preserves
    #      every channel from the binary block (`n_full = len(binary)//4`);
    #      this comment block formalises that contract for BUG-9.
    if spec.energy_cal is not None and len(spec.energy_cal) >= 1:
        a0 = float(spec.energy_cal[0])
        if a0 < 0.0:
            spec.extras["calibration_drift_left"] = True
            spec.extras["calibration_drift_a0_keV"] = a0
            # Count leading channels with E < 0 (under the stored polynomial).
            n_neg = 0
            for ch in range(n_full):
                if _energy_at(ch, spec.energy_cal) < 0.0:
                    n_neg += 1
                else:
                    break
            spec.extras["calibration_drift_neg_energy_channels"] = int(n_neg)
        else:
            spec.extras["calibration_drift_left"] = False

    # --- apply energy ceiling if requested and a calibration is available ---
    if apply_energy_ceiling and spec.energy_cal is not None:
        ceiling = ceiling_keV if ceiling_keV is not None else ENERGY_CEILING_KEV
        keep_n, e_max_kept = _apply_ceiling(counts, spec.energy_cal, ceiling)
        if keep_n < n_full:
            spec.counts = counts[:keep_n].copy()
            spec.dropped_overflow_count = 0  # nothing flagged as overflow here
            spec.energy_max_keV_kept = float(e_max_kept)
        else:
            spec.energy_max_keV_kept = float(
                _energy_at(n_full - 1, spec.energy_cal)
            )
    else:
        if spec.energy_cal is not None:
            spec.energy_max_keV_kept = float(
                _energy_at(n_full - 1, spec.energy_cal)
            )

    spec.n_channels = int(len(spec.counts))

    # --- filename hints (uses existing token parser, RU-aware) ---
    spec.filename_tokens = parse_filename(p.name)
    # If filename token parser flags background, propagate (rare for .spe
    # which usually contains TYPE=Калибровка, but check anyway)
    if spec.filename_tokens.get("is_background_hint"):
        spec.is_background = True

    return spec


# ============================================================================
# Header parsing
# ============================================================================

def _parse_kv_lines(header: str) -> dict:
    """
    Split CP-1251 decoded header into KEY=VALUE pairs.

    Lines without a `=` that follow `PEAKS=N` are accumulated under the
    synthetic key `PEAKS_TABLE` until N rows have been collected. The
    same logic applies to `ZONES=N` → `ZONES_TABLE`. This matches the
    Lsrm convention where rows come immediately after the count line.
    """
    fields: dict = {}
    peaks_table: list = []
    zones_table: list = []
    accumulating: Optional[str] = None  # "peaks" or "zones" or None
    expected_rows_remaining = 0

    for raw_line in header.split("\r\n"):
        line = raw_line.rstrip("\r")
        if not line:
            # blank line ends an accumulation block
            accumulating = None
            expected_rows_remaining = 0
            continue

        # If we're accumulating peaks/zones, eat lines until we hit the
        # promised count or a new KEY=VALUE
        if accumulating is not None and expected_rows_remaining > 0:
            # Detect if this line is actually a new KEY= field — some
            # files write fewer rows than declared
            if "=" in line and not line[0].isdigit() and not line[0] in "-+. ":
                accumulating = None
                expected_rows_remaining = 0
                # fall through to normal KEY=VALUE handling below
            else:
                if accumulating == "peaks":
                    peaks_table.append(line)
                elif accumulating == "zones":
                    zones_table.append(line)
                expected_rows_remaining -= 1
                if expected_rows_remaining == 0:
                    accumulating = None
                continue

        if "=" in line:
            key, _, value = line.partition("=")
            key = key.strip()
            value = value.strip()
            fields[key] = value
            if key == "PEAKS":
                try:
                    expected_rows_remaining = int(value)
                    accumulating = "peaks" if expected_rows_remaining > 0 else None
                except ValueError:
                    accumulating = None
                    expected_rows_remaining = 0
            elif key == "ZONES":
                try:
                    expected_rows_remaining = int(value)
                    accumulating = "zones" if expected_rows_remaining > 0 else None
                except ValueError:
                    accumulating = None
                    expected_rows_remaining = 0

    if peaks_table:
        fields["PEAKS_TABLE"] = "\n".join(peaks_table)
    if zones_table:
        fields["ZONES_TABLE"] = "\n".join(zones_table)
    return fields


def _parse_peaks_table(table_text: str) -> list:
    """
    Parse the PEAKS rows table per Lsrm spec §7.5.2.1.

    Each row is whitespace-separated:
      position d_pos energy d_energy fwhm d_fwhm area d_area chi2 \\
        left_bound right_bound [;-separated nuclide list]

    Some files include only the first 10 numeric fields; nuclide list is
    optional and follows in trailing tokens.

    Returns: list of dicts with parsed fields. Returns empty list on any
    parse failure for robustness.
    """
    result = []
    for line in table_text.splitlines():
        line = line.strip()
        if not line:
            continue
        # Split off nuclide list if any
        nuclides = []
        if ";" in line:
            parts = line.split(";")
            line = parts[0].strip()
            nuclides = [p.strip() for p in parts[1:] if p.strip()]
        tokens = line.split()
        if len(tokens) < 7:
            continue
        try:
            row = {
                "position_ch": float(tokens[0]),
                "d_position_ch": float(tokens[1]),
                "energy_keV": float(tokens[2]),
                "d_energy_keV": float(tokens[3]),
                "fwhm_keV": float(tokens[4]),
                "d_fwhm_keV": float(tokens[5]),
                "area": float(tokens[6]),
            }
            if len(tokens) > 7:
                row["d_area"] = float(tokens[7])
            if len(tokens) > 8:
                row["chi2"] = float(tokens[8])
            if len(tokens) > 9:
                row["left_bound"] = float(tokens[9])
            if len(tokens) > 10:
                row["right_bound"] = float(tokens[10])
            if nuclides:
                row["nuclides"] = nuclides
            result.append(row)
        except ValueError:
            continue
    return result


def _parse_zones_table(table_text: str) -> list:
    """
    Parse the ZONES rows table per Lsrm spec §7.5.2.1.

    Each row format (whitespace-separated):
      left_bound right_bound n_peaks_in_zone minimize_params bg_poly_degree

    Where minimize_params is a comma-separated list (e.g.
    "FWHM, Position, Step, Linear") embedded in the row.

    Useful for Phase 2.1b multiplet deconvolution: each ZONE with
    n_peaks_in_zone > 1 is a candidate for joint Gaussian fit.

    Returns: list of dicts.
    """
    result = []
    for line in table_text.splitlines():
        line = line.strip()
        if not line:
            continue
        # The minimize_params contains commas/spaces; we need to find
        # boundary by recognising the trailing numeric field (poly degree).
        # Robust approach: parse leading 3 floats, trailing 1 int, middle as
        # text.
        tokens = line.split()
        if len(tokens) < 5:
            continue
        try:
            left = float(tokens[0])
            right = float(tokens[1])
            n_peaks = int(tokens[2])
            # last token is the bg polynomial degree
            bg_degree = int(tokens[-1])
            # middle tokens form the minimize specification
            minimize_str = " ".join(tokens[3:-1])
            result.append({
                "left_bound": left,
                "right_bound": right,
                "n_peaks_in_zone": n_peaks,
                "minimize": minimize_str,
                "bg_polynomial_degree": bg_degree,
            })
        except (ValueError, IndexError):
            continue
    return result


def _parse_polynomial_line(value: str) -> list:
    """
    Parse Lsrm `ENERGY=` / `FWHM=` polynomial line.

    Format: `N, a0, a1, a2, …` where N is a degree-marker (small int)
    followed by 7 coefficients (the slot count is fixed at 7 in all
    observed files; degree-marker tells how many of them are
    meaningful, but the unused tail is filled with zeros). We ignore
    the degree-marker and return the full coefficient list low-to-
    high; the caller trims trailing zeros.

    `D` exponent form (`-4.313128E-07`) and `,` decimals (`3,008911`)
    have both been observed across the fixtures and are tolerated.
    """
    parts = [p.strip() for p in value.split(",") if p.strip()]
    if not parts:
        return []

    # First element is the degree marker
    try:
        _degree_marker = int(parts[0])
    except ValueError:
        # Some files might omit it — fall through and treat all as coefs
        return [_safe_float(s) for s in parts]

    coefs = []
    for s in parts[1:]:
        v = _safe_float(s)
        if v is None:
            continue
        coefs.append(v)
    return coefs


def _safe_float(s: str) -> Optional[float]:
    """Parse a float, tolerating `D` exponents and stray commas."""
    if s is None:
        return None
    s = s.strip()
    if not s:
        return None
    # Lsrm uses 'E' for exponents but be defensive
    s2 = s.replace("D", "E").replace("d", "E")
    try:
        return float(s2)
    except ValueError:
        # Last resort: replace comma-as-decimal
        if "," in s2 and "." not in s2:
            try:
                return float(s2.replace(",", "."))
            except ValueError:
                return None
        return None


def _get_float(fields: dict, key: str, default: float = 0.0) -> float:
    v = _safe_float(fields.get(key, ""))
    return float(v) if v is not None else default


def _get_int(fields: dict, key: str, default: int = 0) -> int:
    raw = fields.get(key, "")
    if not raw:
        return default
    try:
        return int(raw)
    except ValueError:
        return default


def _parse_value_err_pair(raw: str) -> Optional[float]:
    """F-130 / v1.17.7 — распарсить пару «значение;погрешность» из
    LSRM-поля типа SAMPLEMASS=1600.0;16.0. Возвращает только value
    (для расчёта плотности погрешность не требуется в первом приближении).
    """
    if not raw:
        return None
    if ";" in raw:
        head, _ = raw.split(";", 1)
        return _safe_float(head)
    return _safe_float(raw)


def _parse_value_err_pair_full(raw: str) -> Optional[tuple]:
    """BUG-1 / 2026-06-02 — распарсить пару «значение;погрешность» из
    LSRM-поля (SAMPLEMASS=1600.0;16.0 → (1600.0, 16.0)).

    Возвращает кортеж (value, uncertainty). Когда uncertainty не задан
    или не парсится — возвращает (value, None). Возвращает None если и
    value не распарсился. Единицы — те же, что в исходном поле (граммы
    для SAMPLEMASS / SAMPLEMASS, миллилитры для SAMPLEVOLUME).
    """
    if not raw:
        return None
    if ";" in raw:
        head, tail = raw.split(";", 1)
        v = _safe_float(head)
        if v is None:
            return None
        u = _safe_float(tail)
        return (v, u)
    v = _safe_float(raw)
    if v is None:
        return None
    return (v, None)


def _auto_extract_density(fields: dict) -> Optional[tuple]:
    """F-130 / v1.17.7 — автоматическое определение ρ_образца (г/см³).

    Логика приоритета:
      1. MATERIAL JSON содержит поле "Ro" (LSRM-конвенция) — это
         прямая запись плотности → ("material_ro", ρ).
      2. SAMPLEMASS (г) и SAMPLEVOLUME (мл = см³) оба заданы → ρ = m/V.
      3. PROBEMASS (г) и PROBEVOLUME (мл) оба заданы → fallback.

    Возвращает (ρ, source_label) либо None, если не удалось вывести.
    Все три источника валидны только при ρ ∈ (0.1, 10.0) г/см³ —
    это санитарный диапазон для гамма-спектрометрии матриц
    (вода=1.0, ОИСН-16=1.6, песок≈2.0, типовые верхние пределы).
    """
    import json as _json

    # ─── 1. MATERIAL.Ro прямой ───────────────────────────────────────
    mat_raw = fields.get("MATERIAL", "")
    if mat_raw:
        try:
            mat = _json.loads(mat_raw)
            ro = mat.get("Ro")
            if ro is not None:
                ro_f = float(ro)
                if 0.1 <= ro_f <= 10.0:
                    return (ro_f, "material_ro")
        except (ValueError, _json.JSONDecodeError, TypeError):
            pass

    # ─── 2. SAMPLEMASS / SAMPLEVOLUME ────────────────────────────────
    m_g = _parse_value_err_pair(fields.get("SAMPLEMASS", ""))
    v_ml = _parse_value_err_pair(fields.get("SAMPLEVOLUME", ""))
    if m_g and v_ml and v_ml > 0:
        rho = m_g / v_ml
        if 0.1 <= rho <= 10.0:
            return (float(rho), "sample_mass_over_volume")

    # ─── 3. PROBEMASS / PROBEVOLUME fallback ─────────────────────────
    m_g = _parse_value_err_pair(fields.get("PROBEMASS", ""))
    v_ml = _parse_value_err_pair(fields.get("PROBEVOLUME", ""))
    if m_g and v_ml and v_ml > 0:
        rho = m_g / v_ml
        if 0.1 <= rho <= 10.0:
            return (float(rho), "probe_mass_over_volume")

    return None


def _auto_extract_mass_kg(fields: dict) -> Optional[tuple]:
    """F-140 / v1.17.7 — авто-извлечение массы образца в КИЛОГРАММАХ.

    LSRM поля SAMPLEMASS / PROBEMASS хранят массу в граммах. Pipeline
    работает в кг через `sample_mass_kg`. Возвращает (mass_kg, source)
    либо None при отсутствии данных. Санитарный диапазон 0.001..100 кг.
    """
    for key, label in (("SAMPLEMASS", "sample_mass_field"),
                       ("PROBEMASS", "probe_mass_field")):
        m_g = _parse_value_err_pair(fields.get(key, ""))
        if m_g and 1.0 <= m_g <= 100_000.0:    # 1 г .. 100 кг
            return (float(m_g) / 1000.0, label)
    return None


def _parse_lsrm_datetime(s: str) -> Optional[datetime]:
    """
    Parse Lsrm datetime in `DD-MM-YY HH:MM:SS[.fff]` or `DD-MM-YY` format.

    The two-digit year is mapped via the classical pivot: 50–99 → 19YY,
    00–49 → 20YY. Sources with this format date from 1999 onward; the
    pivot covers their realistic span (1999 → 2049).

    Lsrm spectrometers since ~2020 emit timestamps with sub-second
    precision (e.g. "21-10-24 14:16:05.80"), so the parser tries the
    fractional-seconds variant first before falling back to the
    integer-seconds and date-only formats. Discovered while validating
    F-29 activity against .src certificate dates (F-30).
    """
    s = (s or "").strip()
    if not s:
        return None
    for fmt in (
        "%d-%m-%y %H:%M:%S.%f",     # 21-10-24 14:16:05.80
        "%d-%m-%y %H:%M:%S",        # 21-10-24 14:16:05
        "%d-%m-%y",                 # 21-10-24
    ):
        try:
            dt = datetime.strptime(s, fmt)
            # %y already does pivot at 1969/2070; not what we want for these
            # files but close enough for practical lab use. Override if needed.
            return dt
        except ValueError:
            pass
    return None


# ============================================================================
# Energy ceiling
# ============================================================================

def _energy_at(ch: int, coefs) -> float:
    """Evaluate polynomial low-to-high at channel ch."""
    return sum(a * (ch ** i) for i, a in enumerate(coefs))


def _apply_ceiling(counts, coefs, ceiling_keV: float):
    """Find the largest channel index whose energy is ≤ ceiling."""
    n = len(counts)
    # The polynomial is monotonically increasing within the channel
    # range for well-formed calibrations — locate the cutoff by simple
    # search rather than algebra to stay robust to degree-3+ shapes.
    keep = n
    for ch in range(n - 1, -1, -1):
        if _energy_at(ch, coefs) <= ceiling_keV:
            keep = ch + 1
            break
        keep = ch
    e_max_kept = _energy_at(max(0, keep - 1), coefs)
    return keep, e_max_kept


# ============================================================================
# Writer
# ============================================================================

def write_lsrm_spe(
    spec: "Spectrum",
    path: str,
    *,
    extra_header: Optional[dict] = None,
) -> None:
    """
    Write a Spectrum to LSRM SpectraLine `.spe` format.

    The file consists of a CP-1251 KEY=VALUE\\r\\n text header that ends
    with the literal marker `SPECTR=`, immediately followed by the
    binary counts block as little-endian uint32 channel values.

    Args:
        spec: Spectrum to serialize.
        path: output filesystem path.
        extra_header: optional dict[str, str] of additional KEY=VALUE
            lines to emit before SPECTR=. Caller may use it to preserve
            arbitrary metadata not covered by Spectrum fields.

    Notes:
        - SHIFR ← spec.sample_id
        - TLIVE/TREAL ← live_time / real_time (4 decimal places, dot)
        - MEASBEGIN ← spec.start_datetime in DD-MM-YY HH:MM:SS
        - GEOMETRY/DETECTOR/OPERATOR/COMMENT mapped from Spectrum fields
        - ENERGY=3, c0, c1, c2, c3, 0, 0, 0 — degree marker + 7 slots
          (matches the convention discovered in the reader)
        - FWHM emitted when spec.stored_fwhm_calibration is set
        - Counts written as uint32-LE; values clipped to uint32 max
          (4,294,967,295) for safety.
    """
    from gamma.spectrum import Spectrum  # local import: writer-only dep

    if not isinstance(spec, Spectrum):
        raise TypeError(f"write_lsrm_spe expects a Spectrum, got {type(spec)!r}")

    p = Path(path)
    counts = np.asarray(spec.counts, dtype=np.int64)
    # Clip safely to uint32 range
    counts_u32 = np.clip(counts, 0, 0xFFFFFFFF).astype("<u4")

    # ----- assemble header dict (ordered) -----
    h: list[tuple[str, str]] = []

    def add(k: str, v) -> None:
        if v is None:
            return
        s = v if isinstance(v, str) else str(v)
        if s != "":
            h.append((k, s))

    # Identity
    add("SHIFR", spec.sample_id)
    add("NOMER", spec.extras.get("lsrm_nomer", ""))
    add("TYPE", spec.extras.get("lsrm_type", ""))
    add("CONFIGNAME", spec.extras.get("lsrm_config", ""))

    # Times
    if spec.start_datetime is not None:
        add("MEASBEGIN", _format_lsrm_datetime(spec.start_datetime))
    add("PREPBEGIN", spec.extras.get("lsrm_prepbegin", ""))
    add("PREPEND", spec.extras.get("lsrm_prepend", ""))
    add("TLIVE", f"{float(spec.live_time):.2f}")
    add("TREAL", f"{float(spec.real_time):.2f}")

    # Operator / geometry / detector
    add("OPERATOR", spec.operator)
    add("GEOMETRY", spec.geometry)
    add("DETECTOR", spec.detector_id)
    add("SETTYPE", spec.extras.get("lsrm_settype", ""))
    add("CONTTYPE", spec.extras.get("lsrm_conttype", ""))

    for k in ("material", "distance", "detradius",
              "rawmass", "probemass", "samplemass",
              "rawvolume", "probevolume", "samplevolume"):
        v = spec.extras.get(f"lsrm_{k}")
        if v:
            add(k.upper(), v)

    # ENBOUNDS
    if spec.energy_max_keV_kept:
        add("ENBOUNDS", f"0,{int(round(spec.energy_max_keV_kept))}")

    # ENERGY polynomial — emit as "3,a0,a1,a2,a3,0,0,0" (degree marker = degree
    # of polynomial, followed by 7 coefficient slots, low-to-high, trailing 0s).
    if spec.energy_cal:
        coefs = list(spec.energy_cal)
        degree = len(coefs) - 1
        # Pad/truncate to 7 slots
        slots = coefs + [0.0] * max(0, 7 - len(coefs))
        slots = slots[:7]
        coef_strs = [_format_lsrm_float(c) for c in slots]
        add("ENERGY", ",".join([str(degree)] + coef_strs))

    # FWHM polynomial (if available)
    if spec.stored_fwhm_calibration and spec.stored_fwhm_calibration.coefficients:
        fw = list(spec.stored_fwhm_calibration.coefficients)
        fw_deg = len(fw) - 1
        slots = fw + [0.0] * max(0, 7 - len(fw))
        slots = slots[:7]
        coef_strs = [_format_lsrm_float(c) for c in slots]
        add("FWHM", ",".join([str(fw_deg)] + coef_strs))

    add("COMMENT", spec.comments)

    # SPECTRSIZE (informational; reader does not require it)
    add("SPECTRSIZE", str(int(len(counts_u32))))

    # Caller-provided extra header lines (last, so they can override)
    if extra_header:
        for k, v in extra_header.items():
            add(str(k), str(v))

    # ----- assemble bytes -----
    header_lines = "".join(f"{k}={v}\r\n" for k, v in h)
    header_bytes = header_lines.encode("cp1251", errors="replace")
    marker = b"SPECTR="
    binary_bytes = counts_u32.tobytes()

    p.write_bytes(header_bytes + marker + binary_bytes)


def _format_lsrm_datetime(dt: datetime) -> str:
    """Emit `DD-MM-YY HH:MM:SS.ff` (sub-second precision, two digits)."""
    return dt.strftime("%d-%m-%y %H:%M:%S") + f".{dt.microsecond // 10000:02d}"


def _format_lsrm_float(x: float) -> str:
    """
    Mimic the float format observed in real LSRM files.

    Coefficients use either `0,0` style or scientific `-4.313128E-07`.
    For round-trip safety we use a generic %.10g format which the
    reader parses with `_safe_float`.
    """
    return f"{float(x):.10g}"
