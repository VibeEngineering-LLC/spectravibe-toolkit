"""
LSRM SpectraLine ASCII `.spe` reader and writer.

This is the **text-export** format produced by ЛСРМ SpectraLine — distinct
from the **binary** LSRM `.spe` (CP-1251 KEY=VALUE header + `SPECTR=`
marker + uint32-LE counts; handled in `gamma.io.lsrm_spe`). Both share
the `.spe` extension and are disambiguated by content sniff (binary has
`SPECTR=` marker; text has `$`-section headers like `$DATA:`).

Despite superficial similarity to the IAEA / ORTEC GammaVision SPE
format, the LSRM ASCII variant follows ЛСРМ-specific conventions which
this module implements strictly. The two are NOT identical — see
"Differences from IAEA/ORTEC" below.

Canonical structure (each section opened by `$TAG:` on its own line,
content lines follow until the next `$` or EOF):

    $SPEC_ID:                       (optional in LSRM, mandatory in IAEA)
    Sample identifier
    $SPEC_REM:                      (optional in LSRM)
    Free-form comments
    $DATE_MEA:
    09/28/2019 14:31:31             (MM/DD/YYYY HH:MM:SS) — or DD-MM-YY
    $MEAS_TIM:
    97829.567706599992 97829.567706599992    (live, real — high-precision)
    $DATA:
    0 4999                          (first_ch  last_ch)
    0                               (counts, one per line, exactly
    389                              last_ch - first_ch + 1 entries)
    509
    ...
    $ROI:                           (optional in LSRM)
    0
    $ENER_FIT:                      (typically the polynomial degree-1
    -3.92435 0.780871                fit for backward compat)
    $MCA_CAL:
    3                               (LSRM: N = NUMBER OF COEFFICIENTS,
    -3.92435 0.780871 2.01905e-06    NOT polynomial degree)
                                    (LSRM: NO unit suffix; IAEA appends "keV")
    $SHAPE_CAL:                     (optional in LSRM, mandatory in IAEA)
    3
    1.2 0.05 0.0
    $ENDRECORD:

Differences from IAEA/ORTEC GammaVision SPE (the SpecUtils `SpeIaea`
parser): although both formats use `$`-section ASCII, they differ in
several places:

  - **$MCA_CAL first line**: in LSRM it is N = number of coefficients
    (3 means 3 floats follow); in IAEA the same number is the
    polynomial degree (3 means cubic, 4 floats follow). Often the same
    number serves both purposes because real cal polynomials are usually
    quadratic (3 coefs = degree 2), but the spec is different. We
    always trust the LSRM convention here.
  - **$MCA_CAL coefficient line**: LSRM emits raw floats only; IAEA
    appends a unit token (`keV` or `MeV`).
  - **$MEAS_TIM**: LSRM writes high-precision floats (sub-millisecond);
    IAEA / GammaVision writes integers (whole seconds) or single decimals.
  - **$SPEC_ID / $SPEC_REM / $ROI / $SHAPE_CAL**: all optional in LSRM;
    IAEA emits them even when empty.
  - **Float notation**: LSRM uses lowercase `e` (`2.01905e-06`); IAEA
    GammaVision uses uppercase `E` with three exponent digits
    (`2.01905E-006`).

References:
  - SpectraLine User Manual (Russian), ЛСРМ — section "Экспорт спектра в
    текстовом формате".
  - SpecUtils `LsrmSpe` parser (https://github.com/sandialabs/SpecUtils);
    treated as a distinct format from `SpeIaea`.
"""

from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np

from gamma.spectrum import (
    ENERGY_CEILING_KEV,
    Spectrum,
    StoredFwhmCalibration,
)


# ============================================================================
# Reader
# ============================================================================

def read_lsrm_spe_text(
    path: str,
    *,
    apply_energy_ceiling: bool = False,
    ceiling_keV: Optional[float] = None,
) -> Spectrum:
    """Read an LSRM SpectraLine ASCII `.spe` file. See module docstring.

    `apply_energy_ceiling` defaults to False since BUG-9 (v1.18.32,
    2026-06-03) — same change as for the binary `read_lsrm_spe`. Opt in
    explicitly when the project-scope 3 MeV trim is needed.
    """
    p = Path(path)
    text = p.read_text(encoding="utf-8", errors="replace")

    sections = _parse_sections(text)

    sample_id = _first_nonempty_line(sections.get("SPEC_ID", ""))
    comments = sections.get("SPEC_REM", "").strip()

    start_dt = _parse_lsrm_datetime(sections.get("DATE_MEA", ""))

    live_time, real_time = 0.0, 0.0
    meas_line = _first_nonempty_line(sections.get("MEAS_TIM", ""))
    if meas_line:
        toks = meas_line.split()
        if len(toks) >= 2:
            try:
                live_time = float(toks[0])
                real_time = float(toks[1])
            except ValueError:
                pass

    counts = _parse_data_section(sections.get("DATA", ""))

    # Energy calibration: prefer $MCA_CAL (LSRM canonical), fall back to
    # $ENER_FIT (legacy linear-only fit kept for backward compat).
    coeffs = _parse_mca_cal_lsrm(sections.get("MCA_CAL", ""))
    if not coeffs:
        line = _first_nonempty_line(sections.get("ENER_FIT", ""))
        if line:
            try:
                vs = [float(x) for x in line.split()]
                if vs:
                    coeffs = tuple(vs)
            except ValueError:
                coeffs = ()

    fwhm_coefs = _parse_shape_cal_lsrm(sections.get("SHAPE_CAL", ""))

    n_full = int(len(counts))
    e_max_kept: Optional[float] = None
    if apply_energy_ceiling and coeffs and n_full:
        ceiling = ceiling_keV if ceiling_keV is not None else ENERGY_CEILING_KEV
        keep_n, e_max_kept = _apply_ceiling(counts, coeffs, ceiling)
        if keep_n < n_full:
            counts = counts[:keep_n].copy()
    elif coeffs and n_full:
        e_max_kept = float(sum(c * (n_full - 1) ** i for i, c in enumerate(coeffs)))

    spec = Spectrum(
        counts=counts,
        live_time=live_time,
        real_time=real_time,
        source_path=str(p),
        source_format="lsrm_spe_text",
        sample_id=sample_id,
        comments=comments,
        start_datetime=start_dt,
        n_channels_raw=n_full,
        n_channels=int(len(counts)),
        channel_pitch=1,
        energy_cal=coeffs if coeffs else None,
        energy_cal_degree=(len(coeffs) - 1) if coeffs else None,
        energy_cal_source="stored" if coeffs else "",
        energy_max_keV_kept=e_max_kept,
    )

    # --- BUG-9 / 2026-06-03: calibration drift diagnostic ---------------------
    # When a0 < 0 the first channels map to E < 0 — a diagnostic signature of
    # spectrometer drift, not a data bug. Channels are PRESERVED; we just
    # surface the condition in extras for downstream consumers. See the
    # matching block in `lsrm_spe.py` for the full rationale.
    if coeffs and len(coeffs) >= 1:
        a0 = float(coeffs[0])
        if a0 < 0.0:
            spec.extras["calibration_drift_left"] = True
            spec.extras["calibration_drift_a0_keV"] = a0
            # Count leading channels with E < 0 in the kept array (so the
            # number reflects what downstream code will actually see). For
            # apply_energy_ceiling=True the leading channels are unaffected
            # (ceiling only trims the high-E tail).
            n_neg = 0
            for ch in range(int(len(counts))):
                e = sum(c * ch ** i for i, c in enumerate(coeffs))
                if e < 0.0:
                    n_neg += 1
                else:
                    break
            spec.extras["calibration_drift_neg_energy_channels"] = int(n_neg)
        else:
            spec.extras["calibration_drift_left"] = False

    if fwhm_coefs:
        spec.stored_fwhm_calibration = StoredFwhmCalibration(
            calibration_peaks=[],
            coefficients=tuple(fwhm_coefs),
            model="lsrm_shape_cal_polynomial",
        )
    return spec


# ============================================================================
# Writer
# ============================================================================

def write_lsrm_spe_text(spec: Spectrum, path: str) -> None:
    """
    Write a Spectrum в формате LSRM SpectraLine ASCII `.spe`.

    Strict conformance to LSRM documentation:
      - `$MCA_CAL` first line is **N = number of coefficients**
        (NOT polynomial degree as in IAEA).
      - `$MCA_CAL` coefficient line carries NO unit suffix.
      - `$MEAS_TIM` emits high-precision floats (12 significant digits).
      - `$SPEC_ID`, `$SPEC_REM`, `$ROI`, `$SHAPE_CAL` blocks are emitted
        only when the corresponding Spectrum field is populated
        (matching the minimal-section style observed in real LSRM exports).
      - Float exponents use lowercase `e` (LSRM convention).
      - Date format: MM/DD/YYYY HH:MM:SS.

    Output text uses LF line endings (the LSRM SpectraLine reader is
    tolerant; the format spec doesn't mandate CRLF).
    """
    counts = np.asarray(spec.counts, dtype=np.int64)
    n = int(len(counts))

    lines: list[str] = []

    def section(name: str, body_lines: list[str]) -> None:
        lines.append(f"${name}:")
        lines.extend(body_lines)

    # Optional sections — only when content is meaningful
    if spec.sample_id:
        section("SPEC_ID", [spec.sample_id])
    if spec.comments:
        section("SPEC_REM", spec.comments.splitlines() or [""])

    # Date (always emitted when present)
    if spec.start_datetime is not None:
        section("DATE_MEA", [spec.start_datetime.strftime("%m/%d/%Y %H:%M:%S")])

    # Times: high-precision floats per LSRM convention
    section(
        "MEAS_TIM",
        [f"{float(spec.live_time):.12g} {float(spec.real_time):.12g}"],
    )

    # Counts
    data_lines = [f"0 {n - 1 if n > 0 else 0}"]
    data_lines.extend(str(int(v)) for v in counts)
    section("DATA", data_lines)

    # Energy calibration
    if spec.energy_cal:
        coefs = list(spec.energy_cal)
        # $ENER_FIT is conventionally the linear (a0 a1) — emit if we have ≥ 2
        if len(coefs) >= 2:
            section(
                "ENER_FIT",
                [_lsrm_format_coefs(coefs)],
            )
        # $MCA_CAL: first line N (number of coefficients), then coefficient
        # line WITHOUT unit suffix.
        section(
            "MCA_CAL",
            [
                str(len(coefs)),
                _lsrm_format_coefs(coefs),
            ],
        )

    # FWHM (shape) calibration
    if spec.stored_fwhm_calibration and spec.stored_fwhm_calibration.coefficients:
        fw = list(spec.stored_fwhm_calibration.coefficients)
        section(
            "SHAPE_CAL",
            [str(len(fw)), _lsrm_format_coefs(fw)],
        )

    section("ENDRECORD", [])

    Path(path).write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


# ============================================================================
# Format sniff helper (used by format_registry to disambiguate `.spe`)
# ============================================================================

def looks_like_lsrm_spe_text(head_bytes: bytes) -> bool:
    """
    Return True if the leading bytes resemble an LSRM ASCII `.spe`.

    Heuristic: file starts (after optional BOM/whitespace) with a `$` line,
    and contains `$DATA:` and/or `$MEAS_TIM:` in the first few KB. Binary
    LSRM `.spe` never starts with `$` (it starts with `SHIFR=` or similar
    CP-1251 KEY=VALUE line).
    """
    head = head_bytes[:512].lstrip(b"\xef\xbb\xbf").lstrip()
    if not head.startswith(b"$"):
        return False
    needles = (b"$DATA:", b"$MEAS_TIM:", b"$DATE_MEA:")
    return any(n in head_bytes[:4096] for n in needles)


# ============================================================================
# Helpers — section parser
# ============================================================================

def _parse_sections(text: str) -> dict:
    """Split text into {SECTION_NAME: body_text} (section name is uppercase)."""
    sections: dict = {}
    current_name: Optional[str] = None
    current_body: list[str] = []
    for raw_line in text.splitlines():
        m = re.match(r"^\$([A-Z_0-9]+):\s*$", raw_line.strip())
        if m:
            if current_name is not None:
                sections[current_name] = "\n".join(current_body)
            current_name = m.group(1)
            current_body = []
        else:
            current_body.append(raw_line)
    if current_name is not None:
        sections[current_name] = "\n".join(current_body)
    return sections


def _first_nonempty_line(body: str) -> str:
    for line in body.splitlines():
        if line.strip():
            return line.strip()
    return ""


def _parse_data_section(body: str) -> np.ndarray:
    """Parse $DATA: section. First non-empty line is `start end`; rest are counts."""
    lines = [ln.strip() for ln in body.splitlines() if ln.strip()]
    if not lines:
        return np.array([], dtype=np.int64)
    counts_text = lines[1:]
    flat: list[int] = []
    for ln in counts_text:
        for tok in ln.split():
            try:
                flat.append(int(float(tok)))
            except ValueError:
                continue
    if not flat:
        for tok in (" ".join(counts_text)).split():
            try:
                flat.append(int(float(tok)))
            except ValueError:
                continue
    return np.array(flat, dtype=np.int64)


def _parse_mca_cal_lsrm(body: str) -> tuple:
    """
    Parse $MCA_CAL section per LSRM convention:
        N           <- number of coefficients
        c0 c1 ...   <- N floats, NO unit suffix

    Tolerant of stray unit tokens for cross-IAEA compatibility (e.g. when
    a file labels itself with `$MCA_CAL: 3 / -3.9 0.78 2e-6 keV`).
    """
    lines = [ln.strip() for ln in body.splitlines() if ln.strip()]
    if len(lines) < 2:
        return ()
    try:
        n_coefs = int(lines[0])
    except ValueError:
        n_coefs = 0
    tokens = lines[1].split()
    numeric: list = []
    for tok in tokens:
        try:
            numeric.append(float(tok))
        except ValueError:
            break  # hit a non-numeric token (likely a unit) — stop
    if n_coefs > 0:
        numeric = numeric[:n_coefs]
    return tuple(numeric)


def _parse_shape_cal_lsrm(body: str) -> tuple:
    lines = [ln.strip() for ln in body.splitlines() if ln.strip()]
    if len(lines) < 2:
        return ()
    try:
        n_coefs = int(lines[0])
    except ValueError:
        n_coefs = 0
    tokens = lines[1].split()
    numeric: list = []
    for tok in tokens:
        try:
            numeric.append(float(tok))
        except ValueError:
            break
    if n_coefs > 0:
        numeric = numeric[:n_coefs]
    return tuple(numeric)


def _parse_lsrm_datetime(body: str) -> Optional[datetime]:
    s = _first_nonempty_line(body)
    if not s:
        return None
    for fmt in (
        "%m/%d/%Y %H:%M:%S",
        "%m/%d/%y %H:%M:%S",
        "%d-%m-%Y %H:%M:%S",
        "%d-%m-%y %H:%M:%S",
        "%d/%m/%Y %H:%M:%S",
        "%Y-%m-%dT%H:%M:%S",
        "%Y-%m-%d %H:%M:%S",
    ):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            continue
    return None


def _lsrm_format_coefs(coefs) -> str:
    """
    Format a coefficient list in LSRM convention:
      - lowercase `e` for exponents
      - 6 significant digits (sufficient for stored calibrations)
      - space-separated
    """
    parts = []
    for c in coefs:
        s = f"{float(c):g}"
        # `%g` already uses lowercase 'e'; safety net:
        parts.append(s.replace("E", "e"))
    return " ".join(parts)


# ============================================================================
# Helpers — energy ceiling (mirrors n42_2012._apply_ceiling)
# ============================================================================

def _apply_ceiling(counts, coefs, ceiling_keV: float):
    """Return (keep_n, e_max_kept)."""
    n = int(len(counts))
    if n == 0 or not coefs:
        return n, None
    energies = np.zeros(n, dtype=np.float64)
    channels = np.arange(n, dtype=np.float64)
    for c in reversed(coefs):
        energies = energies * channels + c
    mask = energies <= ceiling_keV
    if mask.all():
        return n, float(energies[-1])
    first_above = int(np.argmax(~mask))
    if first_above == 0:
        return 0, None
    return first_above, float(energies[first_above - 1])
