"""
Reader for Lsrm SpectraLine nuclide library XML files (.lib).

Format: SpectraLine Nuclear Data Information File (LSRM), schema
`library_type="gamma"`, library_version 2.0, database_version 26+.

Encoding: Windows-1251 (CP-1251), declared in the XML prolog.

Top-level structure:

    <?xml version="1.0" encoding="windows-1251"?>
    <Library library_type="gamma" library_version="2.0" database_version="26">
      <Comment></Comment>
      <Nuclide name="Cs-137" half_life_value="30,05" half_life_unit="year"
               gamma_constant="0,0789" atomic_mass="137">
        <Line energy="661,657" d_energy="0,003" intensity="85,1"
              d_intensity="0,2" [line_type="X"] [used="false"]/>
        ...
      </Nuclide>
      ...
    </Library>

Notes:
  • Decimal separator is COMMA (Russian locale convention)
  • `line_type="X"` marks X-ray (typically excluded from γ-only analyses)
  • `used="false"` marks lines the calibration didn't use (low intensity,
    interferences, etc.) — we preserve them but mark as `used=False`
  • Half-life units encountered: year, day, hour, second
  • This is a DETECTOR-SPECIFIC library that has been pre-curated for
    Gamma-1S NaI 63×63; line selection reflects what's reliably
    resolvable on that hardware
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
# SEC-01 (P1) hardening: parse untrusted .lib content via defusedxml to block
# billion-laughs / DOCTYPE entity DoS. Stdlib ET kept for ParseError type
# compatibility in caller error handling (defusedxml re-uses ET.ParseError).
from defusedxml.ElementTree import fromstring as _safe_fromstring
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass(frozen=True)
class LibraryLine:
    """One γ-line in a nuclide library entry."""
    energy_keV: float
    energy_uncertainty_keV: float
    intensity_pct: float
    intensity_uncertainty_pct: float
    line_type: Optional[str] = None       # "X" for X-ray, None for γ
    used: bool = True                     # False if pre-marked as unused


@dataclass(frozen=True)
class LibraryNuclide:
    """One nuclide entry from a library."""
    name: str
    half_life_value: Optional[float]
    half_life_unit: Optional[str]         # "year"/"day"/"hour"/"second"
    gamma_constant: Optional[float]
    atomic_mass: Optional[float]
    lines: tuple                          # tuple of LibraryLine

    @property
    def half_life_seconds(self) -> Optional[float]:
        """Convert half-life to seconds for consistency."""
        if self.half_life_value is None or self.half_life_unit is None:
            return None
        unit = self.half_life_unit.lower()
        factor = {
            "second": 1.0, "s": 1.0, "sec": 1.0,
            "minute": 60.0, "min": 60.0,
            "hour": 3600.0, "h": 3600.0, "hr": 3600.0,
            "day": 86400.0, "d": 86400.0,
            "year": 365.25 * 86400.0, "y": 365.25 * 86400.0, "yr": 365.25 * 86400.0,
        }.get(unit, None)
        if factor is None:
            return None
        return self.half_life_value * factor


@dataclass(frozen=True)
class LsrmLibrary:
    """Parsed Lsrm nuclide library."""
    path: str
    library_type: str                    # "gamma" / "alpha" / ...
    library_version: str
    database_version: str
    comment: str
    nuclides: tuple                      # tuple of LibraryNuclide

    def get(self, name: str) -> Optional[LibraryNuclide]:
        """Lookup by name (case-insensitive)."""
        nm = name.lower()
        for n in self.nuclides:
            if n.name.lower() == nm:
                return n
        return None

    def names(self) -> list:
        return [n.name for n in self.nuclides]


def _parse_decimal(value: Optional[str]) -> Optional[float]:
    """Parse Russian-locale float (comma decimal) or English float."""
    if value is None or value == "":
        return None
    s = value.strip().replace(",", ".")
    try:
        return float(s)
    except ValueError:
        return None


def read_lsrm_library(path) -> LsrmLibrary:
    """
    Parse an Lsrm SpectraLine nuclide library file.

    Args:
        path: file path (string or Path). CP-1251 XML.

    Returns:
        LsrmLibrary.

    Raises:
        FileNotFoundError, xml.etree.ElementTree.ParseError on bad XML.
    """
    path = Path(path)
    with open(path, "rb") as f:
        data = f.read()
    # XML declares encoding="windows-1251"; let ET handle decoding.
    # SEC-01: _safe_fromstring blocks DOCTYPE / entity expansion.
    root = _safe_fromstring(data)

    library_type = root.get("library_type", "")
    library_version = root.get("library_version", "")
    database_version = root.get("database_version", "")
    comment_el = root.find("Comment")
    comment = comment_el.text or "" if comment_el is not None else ""

    nuclides = []
    for nuc_el in root.findall("Nuclide"):
        name = nuc_el.get("name", "").strip()
        hl_value = _parse_decimal(nuc_el.get("half_life_value"))
        hl_unit = nuc_el.get("half_life_unit")
        gamma_const = _parse_decimal(nuc_el.get("gamma_constant"))
        atomic_mass = _parse_decimal(nuc_el.get("atomic_mass"))

        lines = []
        for line_el in nuc_el.findall("Line"):
            E = _parse_decimal(line_el.get("energy"))
            dE = _parse_decimal(line_el.get("d_energy")) or 0.0
            I = _parse_decimal(line_el.get("intensity")) or 0.0
            dI = _parse_decimal(line_el.get("d_intensity")) or 0.0
            line_type = line_el.get("line_type")  # "X" or None
            used_raw = line_el.get("used", "true")
            used = used_raw.lower() != "false"
            if E is None:
                continue
            lines.append(LibraryLine(
                energy_keV=E,
                energy_uncertainty_keV=dE,
                intensity_pct=I,
                intensity_uncertainty_pct=dI,
                line_type=line_type,
                used=used,
            ))

        nuclides.append(LibraryNuclide(
            name=name,
            half_life_value=hl_value,
            half_life_unit=hl_unit,
            gamma_constant=gamma_const,
            atomic_mass=atomic_mass,
            lines=tuple(lines),
        ))

    return LsrmLibrary(
        path=str(path),
        library_type=library_type,
        library_version=library_version,
        database_version=database_version,
        comment=comment,
        nuclides=tuple(nuclides),
    )


def merge_lsrm_library_into_internal(
    lsrm_lib: LsrmLibrary,
    include_xrays: bool = False,
    include_unused: bool = False,
    exclude_nuclides: tuple = ("Th-228",),
) -> dict:
    """
    Convert an LsrmLibrary into the internal nuclide_library dict shape
    used by gamma.data.nuclide_library.

    Internal shape (per nuclide):
        {
            "T_half_s": float,
            "lines": [[E_keV, I_pct, dI_pct], ...]
        }

    Args:
        lsrm_lib: parsed LsrmLibrary
        include_xrays: include X-ray lines (line_type=="X")? Default False.
        include_unused: include lines marked used="false"? Default False.
        exclude_nuclides: tuple of nuclide names to skip. Default:
            ("Th-228",) — Th-228 is a Th-232-chain daughter and in
            secular equilibrium its lines overlap Pb-212/Ac-228/Tl-208
            (the same chain). Including both Th-228 and Th-232 creates
            duplicate identification candidates for the same physical
            decay chain. We use Th-232 as the parent identifier and
            its daughters (Pb-212, Bi-212, Tl-208, Ac-228) as chain
            indicators; Th-228 by itself is ambiguous and excluded.

    Returns:
        dict {nuclide_name: {"T_half_s": ..., "lines": ...}}
    """
    exclude_set = set(exclude_nuclides)
    out = {}
    for nuc in lsrm_lib.nuclides:
        if nuc.name in exclude_set:
            continue
        lines_filtered = []
        for ln in nuc.lines:
            if not include_xrays and ln.line_type == "X":
                continue
            if not include_unused and not ln.used:
                continue
            lines_filtered.append([
                ln.energy_keV, ln.intensity_pct, ln.intensity_uncertainty_pct,
            ])
        if not lines_filtered:
            continue
        lines_filtered.sort(key=lambda x: x[0])
        entry = {"lines": lines_filtered}
        T = nuc.half_life_seconds
        if T is not None:
            entry["T_half_s"] = T
        out[nuc.name] = entry
    return out


__all__ = [
    "LibraryLine", "LibraryNuclide", "LsrmLibrary",
    "read_lsrm_library", "merge_lsrm_library_into_internal",
]
