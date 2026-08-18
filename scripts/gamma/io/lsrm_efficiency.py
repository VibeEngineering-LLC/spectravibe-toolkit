"""
Reader for Lsrm SpectraLine efficiency files (.efa / .efr).

Two related formats, both CP-1251 text with CRLF line endings:

  • `.efa` — Aggregated efficiency for one geometry. Single header
    block followed by energy points combining all reference sources.

  • `.efr` — Raw per-source measurements. Each block is one reference
    source measurement; multiple blocks per file.

Both formats use a section-based syntax similar to INI files:

    [detector_name;geometry_name[;source_id]]
    Key=Value
    ...
    E_keV=epsilon,dEpsilon_pct,nuclide,S_counts,dS_counts,I_pct
    ...

Energy lines (those whose key parses as a float) are the efficiency
data points. Other lines are metadata.

Reference: Lsrm SpectraLine documentation, format version 1.7.11918+
(observed on Gamma-1S NaI 63×63 supplied as test fixtures).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass(frozen=True)
class EfficiencyPoint:
    """One reference (energy, efficiency) calibration point."""

    energy_keV: float
    efficiency: float
    efficiency_uncertainty_pct: float
    source_nuclide: str
    peak_area_counts: float
    peak_area_uncertainty_counts: float
    intensity_pct: float


@dataclass(frozen=True)
class EfficiencyBlock:
    """One block from an .efa or .efr file — header + metadata + points."""

    header: str                          # e.g. "[detector;geometry;source]"
    detector: str
    geometry: str
    source_id: Optional[str]             # set only in .efr per-source blocks
    metadata: dict                       # all Key=Value pairs from header
    points: tuple                        # tuple of EfficiencyPoint

    @property
    def volume_ml(self) -> Optional[float]:
        try:
            return float(self.metadata.get("Volume,ml", "0"))
        except (ValueError, TypeError):
            return None

    @property
    def distance_cm(self) -> Optional[float]:
        try:
            return float(self.metadata.get("Distance,cm", "0"))
        except (ValueError, TypeError):
            return None

    @property
    def density_g_per_cm3(self) -> Optional[float]:
        raw = self.metadata.get("Density,g/cm3", "")
        try:
            return float(raw)
        except (ValueError, TypeError):
            return None  # may be "not essential"


@dataclass(frozen=True)
class EfficiencyFile:
    """Parsed contents of an .efa or .efr file."""

    path: str
    file_format: str              # "efa" or "efr"
    blocks: tuple                 # tuple of EfficiencyBlock

    def all_points(self) -> list:
        """Flatten all efficiency points from all blocks."""
        result = []
        for b in self.blocks:
            result.extend(b.points)
        return result

    def by_geometry(self) -> dict:
        """Group blocks by geometry name."""
        result = {}
        for b in self.blocks:
            result.setdefault(b.geometry, []).append(b)
        return result


# Regex for energy lines (key is a float)
_ENERGY_LINE = re.compile(
    r"^\s*([+-]?\d+(?:\.\d*)?(?:[eE][+-]?\d+)?)\s*=\s*(.*?)\s*$"
)
_HEADER_LINE = re.compile(r"^\s*\[(.+)\]\s*$")
_KV_LINE = re.compile(r"^\s*([^=]+?)\s*=\s*(.*?)\s*$")


def _parse_energy_value(value_str: str) -> Optional[tuple]:
    """
    Parse the right-hand side of an energy line.

    Expected format (after the '='):
        epsilon, dEps_pct, nuclide, S_counts, dS_counts, I_pct[, extra]

    Some .efr files have trailing extra columns (e.g. " No" flag).
    Returns 6-tuple or None on parse failure.
    """
    parts = [p.strip() for p in value_str.split(",")]
    if len(parts) < 6:
        return None
    try:
        eps = float(parts[0])
        d_eps_pct = float(parts[1])
        nuclide = parts[2]
        s_counts = float(parts[3])
        ds_counts = float(parts[4])
        i_pct = float(parts[5])
    except (ValueError, IndexError):
        return None
    return (eps, d_eps_pct, nuclide, s_counts, ds_counts, i_pct)


def read_efficiency_file(path) -> EfficiencyFile:
    """
    Parse an .efa or .efr Lsrm efficiency file.

    Args:
        path: file path (string or Path). Must exist and be CP-1251
              encoded text.

    Returns:
        EfficiencyFile with all parsed blocks.

    Raises:
        FileNotFoundError, UnicodeDecodeError on I/O failure.
    """
    path = Path(path)
    file_format = path.suffix.lstrip(".").lower()
    if file_format not in ("efa", "efr"):
        raise ValueError(f"Unknown file format: {path.suffix} "
                         f"(expected .efa or .efr)")

    with open(path, "rb") as f:
        data = f.read()
    try:
        text = data.decode("cp1251")
    except UnicodeDecodeError:
        text = data.decode("latin1")

    blocks = []
    current_header = None
    current_metadata = {}
    current_points = []

    def _finalise_block():
        if current_header is None:
            return None
        # Parse header components
        parts = current_header.split(";")
        detector = parts[0].strip() if len(parts) > 0 else ""
        geometry = parts[1].strip() if len(parts) > 1 else ""
        source_id = parts[2].strip() if len(parts) > 2 else None
        return EfficiencyBlock(
            header=current_header,
            detector=detector,
            geometry=geometry,
            source_id=source_id,
            metadata=dict(current_metadata),
            points=tuple(current_points),
        )

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        # Section header?
        m = _HEADER_LINE.match(line)
        if m:
            # Finalise previous block, start new one
            prev = _finalise_block()
            if prev is not None:
                blocks.append(prev)
            current_header = m.group(1)
            current_metadata = {}
            current_points = []
            continue

        # Energy line? (key parses as a number)
        m_energy = _ENERGY_LINE.match(line)
        if m_energy:
            try:
                E_keV = float(m_energy.group(1))
                if E_keV > 0:
                    parsed = _parse_energy_value(m_energy.group(2))
                    if parsed is not None:
                        eps, d_eps_pct, nuc, s_c, ds_c, i_pct = parsed
                        current_points.append(EfficiencyPoint(
                            energy_keV=E_keV,
                            efficiency=eps,
                            efficiency_uncertainty_pct=d_eps_pct,
                            source_nuclide=nuc,
                            peak_area_counts=s_c,
                            peak_area_uncertainty_counts=ds_c,
                            intensity_pct=i_pct,
                        ))
                        continue
            except ValueError:
                pass  # fall through to KV handling

        # Generic Key=Value
        m_kv = _KV_LINE.match(line)
        if m_kv:
            key = m_kv.group(1).strip()
            value = m_kv.group(2).strip()
            current_metadata[key] = value

    # Final block
    last = _finalise_block()
    if last is not None:
        blocks.append(last)

    return EfficiencyFile(
        path=str(path),
        file_format=file_format,
        blocks=tuple(blocks),
    )


__all__ = [
    "EfficiencyPoint", "EfficiencyBlock", "EfficiencyFile",
    "read_efficiency_file",
]
