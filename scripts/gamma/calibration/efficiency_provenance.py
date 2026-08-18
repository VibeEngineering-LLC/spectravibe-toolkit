"""
T41 (BUG-40 (b) hardening) — efficiency-file detector content fingerprint gate.

Anti-hallucination class: silent **content** fallback when the file-path lookup
succeeds (efficiency_autoload found an .efr) but the .efr file's `Detector=`
field encodes a DIFFERENT physical instance than the spectrum was acquired on.

Real incident (BUG-40 (b), 2026-06-23): `detectors/Gamma-1S/efficiency/
Gamma-1S_NaI_63x63_USB_SN-01/...Marinelli.efr` is named after Gamma-1S but
the .efr `[detector;geometry;source]` header records `Detector=
УДС-ГЦ-63х63-USB №SN-01` (Поверка-2024) while the .spe spectrum's
CONFIGNAME is `Гамма-1С №SN-02` (Поверка-2016). Serial 0086 != 0221 ->
wrong physical instrument's efficiency curve -> activity bias of -96% to -97%
on Am-241/Ti-44, +9.5% on Cs-137. The path-level cyrillic_to_latin_collision
predicate does NOT catch this (path is fine, content is wrong).

This module provides the missing content-side check:

  check_efr_detector_match(efr_path, expected_lsrm_config) -> dict | None

Returns a mismatch descriptor (suitable for embedding into
`StagedResult.detector_fallback`) when serial-year extracted from the .efr's
`Detector=` field differs from serial-year extracted from the .spe's
CONFIGNAME. Returns None on match, on missing serials (cannot decide), or on
read failure (no false-positive on broken files; the existing
EFFICIENCY_FIT_FAILED sentinel covers that path).
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional


_SERIAL_YEAR_RE = re.compile(r"N[oO]?\s*[-_]?\s*(\d{3,4})[-_](\d{2})", re.UNICODE)
_SERIAL_YEAR_RE_CYR = re.compile(r"№\s*(\d{3,4})[-_](\d{2})")


def extract_serial_year(detector_str: str) -> Optional[tuple]:
    """Pull (serial, year) out of an LSRM-format detector string.

    Accepts both Cyrillic numero sign U+2116 (`Гамма-1С №SN-02`) and ASCII
    `No SN-02` / `N SN-02` variants. Returns None if nothing parses.
    """
    if not detector_str:
        return None
    m = _SERIAL_YEAR_RE_CYR.search(detector_str)
    if m:
        return (m.group(1), m.group(2))
    m = _SERIAL_YEAR_RE.search(detector_str)
    if m:
        return (m.group(1), m.group(2))
    return None


def check_efr_detector_match(
    efr_path: str,
    expected_detector_str: str,
) -> Optional[dict]:
    """Compare serial-year from .efr `Detector=` against expected source.

    Args:
        efr_path: absolute path to a successfully-loaded .efr/.efa file
            (caller has already validated parseability through
            ``fit_efficiency_from_efr_file``).
        expected_detector_str: source-of-truth detector string from the
            spectrum being analyzed. Prefer ``spec.extras["lsrm_config"]``
            (LSRM CONFIGNAME); fall back to ``spec.detector_id``.

    Returns:
        ``None`` if (a) serials match, (b) either side has no extractable
        serial (insufficient evidence — silent, do not cry wolf), or
        (c) the file cannot be read (a separate failure mode already
        surfaced by ``EFFICIENCY_FIT_FAILED``).

        Otherwise a dict suitable for embedding into
        ``detector_fallback["efficiency_detector_mismatch"]``::

            {
                "code": "EFFICIENCY_DETECTOR_SERIAL_MISMATCH",
                "expected_detector": str,
                "actual_detector": str,
                "expected_serial_year": [serial, year],
                "actual_serial_year": [serial, year],
                "efr_file_basename": str,
            }

        F-115: only the basename is included; absolute operator paths are
        not embedded into operator-facing report fields.
    """
    if not expected_detector_str or not efr_path:
        return None

    expected_sy = extract_serial_year(expected_detector_str)
    if not expected_sy:
        return None

    try:
        from gamma.io.lsrm_efficiency import read_efficiency_file
        eff_file = read_efficiency_file(efr_path)
    except Exception:
        return None

    if not eff_file.blocks:
        return None

    block = eff_file.blocks[0]
    actual_str = block.metadata.get("Detector", "") or block.detector or ""
    actual_sy = extract_serial_year(actual_str)
    if not actual_sy:
        return None

    if actual_sy == expected_sy:
        return None

    return {
        "code": "EFFICIENCY_DETECTOR_SERIAL_MISMATCH",
        "expected_detector": expected_detector_str,
        "actual_detector": actual_str,
        "expected_serial_year": list(expected_sy),
        "actual_serial_year": list(actual_sy),
        "efr_file_basename": Path(efr_path).name,
    }


__all__ = [
    "extract_serial_year",
    "check_efr_detector_match",
]