"""
XRF line catalog lookups.

Loads `data/xrf_lines.json` once on first access. The catalog covers the
frequently-used elements only (variant b from the Phase 1.1 design); for
rarely-encountered Z, refer to references/08b_xrf_lines_catalog.md.

Key API:
    get_element(symbol) -> dict | None
    lookup_by_energy(E_keV, tolerance_keV) -> list[XrfLine]
    elements_in_catalog() -> list[str]
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Optional

from gamma.data import DATA_DIR

_XRF_PATH = DATA_DIR / "xrf_lines.json"

_CACHE: Optional[dict] = None


@dataclass(frozen=True)
class XrfLine:
    """A single XRF line, with the element and line label."""
    element: str
    Z: int
    line: str            # "Ka1", "Ka2", "Kb1", "Kb2", "La1", "Lb1", "Lb2", "Lg1"
    E_keV: float
    delta_keV: float = 0.0

    def __repr__(self) -> str:
        return (f"XrfLine({self.element} {self.line} @ {self.E_keV:.2f} keV, "
                f"Z={self.Z}, Δ={self.delta_keV:+.2f})")


_K_LABELS = ("Ka1", "Ka2", "Kb1", "Kb2")
_L_LABELS = ("La1", "Lb1", "Lb2", "Lg1")


def _load() -> dict:
    global _CACHE
    if _CACHE is not None:
        return _CACHE
    if not _XRF_PATH.is_file():
        raise FileNotFoundError(
            f"xrf_lines.json not found at {_XRF_PATH}. "
            f"Expected: <skill_root>/data/xrf_lines.json"
        )
    with _XRF_PATH.open("r", encoding="utf-8") as f:
        data = json.load(f)
    _CACHE = {k: v for k, v in data.items() if not k.startswith("_")}
    return _CACHE


def reset_cache() -> None:
    global _CACHE
    _CACHE = None


def elements_in_catalog() -> list:
    """Return all element symbols in the catalog."""
    return sorted(_load().keys())


def get_element(symbol: str) -> Optional[dict]:
    """Return the full record for an element, or None."""
    return _load().get(symbol)


def lookup_by_energy(E_keV: float, tolerance_keV: float) -> list:
    """
    Find all XRF lines within ±tolerance_keV of E_keV.

    Returns a list of XrfLine, sorted by |delta_keV| ascending.

    Only catalog elements are searched. Use references/08b for rare elements.
    """
    hits = []
    for symbol, rec in _load().items():
        Z = rec.get("Z", 0)

        # K lines
        K = rec.get("K") or []
        for label, energy in zip(_K_LABELS, K):
            if energy is None:
                continue
            delta = energy - E_keV
            if abs(delta) <= tolerance_keV:
                hits.append(XrfLine(
                    element=symbol, Z=Z, line=label,
                    E_keV=energy, delta_keV=delta,
                ))

        # L lines (heavy elements only)
        L = rec.get("L") or []
        for label, energy in zip(_L_LABELS, L):
            if energy is None:
                continue
            delta = energy - E_keV
            if abs(delta) <= tolerance_keV:
                hits.append(XrfLine(
                    element=symbol, Z=Z, line=label,
                    E_keV=energy, delta_keV=delta,
                ))

    hits.sort(key=lambda h: abs(h.delta_keV))
    return hits


def expected_partner(line: XrfLine) -> Optional[float]:
    """
    Given an observed XRF line, return the expected energy of its partner
    (for doublet/triplet check).

      Ka1 -> Ka2
      Ka2 -> Ka1
      Kb1 -> Kb2 if available, else Ka1
      Kb2 -> Kb1

    Returns None if no useful partner is defined.
    """
    rec = get_element(line.element)
    if rec is None:
        return None
    K = rec.get("K") or []
    if line.line == "Ka1":
        return K[1] if len(K) > 1 and K[1] is not None else None
    if line.line == "Ka2":
        return K[0] if len(K) > 0 and K[0] is not None else None
    if line.line == "Kb1":
        if len(K) > 3 and K[3] is not None:
            return K[3]
        return K[0] if len(K) > 0 and K[0] is not None else None
    if line.line == "Kb2":
        return K[2] if len(K) > 2 and K[2] is not None else None
    return None
