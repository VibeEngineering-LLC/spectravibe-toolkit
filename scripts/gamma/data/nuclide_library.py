"""
Nuclide library lookups.

Loads `data/nuclides.json` once on first access; subsequent calls hit the
in-memory cache. Module-level state is acceptable here because the data
is small (~5 KB) and read-only.

Key API:
    get_nuclide(name) -> dict | None
    lookup_by_energy(E_keV, tolerance_keV) -> list[NuclideLine]
    nuclides_in_chain(chain_name) -> list[str]
    list_nuclides() -> list[str]

Each "NuclideLine" returned by lookup_by_energy is a small dataclass
holding (nuclide_name, E_keV, I_gamma_pct, sigma_I_pct, delta_keV).
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from gamma.data import DATA_DIR

_NUCLIDES_PATH = DATA_DIR / "nuclides.json"

# Module-level cache. Loaded on first access.
_CACHE: Optional[dict] = None


@dataclass(frozen=True)
class NuclideLine:
    """A single library line, with the parent nuclide name."""
    nuclide: str
    E_keV: float
    I_gamma_pct: float
    sigma_I_pct: float
    delta_keV: float = 0.0  # set by lookup_by_energy

    def __repr__(self) -> str:
        return (f"NuclideLine({self.nuclide} @ {self.E_keV:.2f} keV, "
                f"I={self.I_gamma_pct:.2f}%, Δ={self.delta_keV:+.2f})")


def _load() -> dict:
    """Load the JSON file lazily. Cached after first call."""
    global _CACHE
    if _CACHE is not None:
        return _CACHE
    if not _NUCLIDES_PATH.is_file():
        raise FileNotFoundError(
            f"nuclides.json not found at {_NUCLIDES_PATH}. "
            f"Expected location: <skill_root>/data/nuclides.json"
        )
    with _NUCLIDES_PATH.open("r", encoding="utf-8") as f:
        data = json.load(f)
    # Strip metadata key (underscore-prefixed entries are schema/docs)
    _CACHE = {k: v for k, v in data.items() if not k.startswith("_")}
    return _CACHE


def reset_cache() -> None:
    """Clear the in-memory cache. Useful for tests."""
    global _CACHE
    _CACHE = None


def load_external_library(
    library_path: str,
    *,
    merge_mode: str = "override",
    include_xrays: bool = False,
    include_unused: bool = False,
    split_chains: bool = True,
) -> int:
    """
    Load an external Lsrm SpectraLine library (.lib XML) and merge it
    into the active in-memory library.

    Use this when you have a detector-specific Lsrm library that is
    more comprehensive or more accurate than the built-in
    `data/nuclides.json`. The detector-specific Lsrm libraries
    typically:

      • Include more nuclides (e.g. Bi-207, Cd-109, Y-88, Ti-44,
        Mn-54, Nb-94, Ce-139, Co-57, Eu-154) than the built-in set
      • Use only "true" γ-lines (excluding 511 keV annihilation and
        X-ray escape peaks that are detector-dependent)
      • Carry detector-specific intensity uncertainties from past
        calibrations
      • Bundle natural-chain daughters (Pb-214, Bi-214) under the
        parent name (Ra-226). This is valid for sealed equilibrium
        calibration sources but WRONG for environmental samples where
        Rn-222 escape breaks the equilibrium.

    Args:
        library_path: path to the .lib XML file
        merge_mode: "override" (default) replaces any existing entries
            with the same name; "supplement" only adds nuclides not
            already present.
        include_xrays: include lines marked line_type="X"
        include_unused: include lines marked used="false"
        split_chains: when True (default), decompose combined chain
            entries (Ra-226, Th-232, U-238, U-235) into separate
            entries for each daughter (Pb-214, Bi-214, Tl-208,
            Pb-212, Ac-228, Bi-212, Th-234, Pa-234m, Th-231, etc.)
            using ENSDF-based line ownership in
            `gamma.data.chain_decomposer`. Set to False to preserve
            the Lsrm combined-chain bundling (only valid for sealed
            equilibrium samples).

    Returns:
        Number of nuclide entries loaded (counting each split daughter
        separately when split_chains=True).
    """
    from gamma.io.lsrm_library import (
        read_lsrm_library, merge_lsrm_library_into_internal,
    )

    lib = read_lsrm_library(library_path)
    external = merge_lsrm_library_into_internal(
        lib, include_xrays=include_xrays, include_unused=include_unused,
    )

    # Optionally split natural-chain parent entries into individual daughters
    if split_chains:
        from gamma.data.chain_decomposer import (
            split_chain_entry, LSRM_CHAIN_PARENT_NAMES,
        )
        new_external = dict(external)
        for parent_name in list(external.keys()):
            if parent_name not in LSRM_CHAIN_PARENT_NAMES:
                continue
            combined_lines = external[parent_name].get("lines", [])
            parent_T = external[parent_name].get("T_half_s")
            split_result = split_chain_entry(combined_lines, parent_name)
            # Remove the old combined entry
            del new_external[parent_name]
            # Add per-daughter entries; preserve parent T½ on the
            # parent entry if it still has lines
            for nuc_name, lines in split_result.items():
                entry = {"lines": lines}
                if nuc_name == parent_name and parent_T is not None:
                    entry["T_half_s"] = parent_T
                # Daughters inherit no T½ here; downstream code should
                # consult ENSDF separately if T½ is needed for daughters
                if nuc_name in new_external and merge_mode == "supplement":
                    continue  # don't overwrite existing entry
                new_external[nuc_name] = entry
        external = new_external

    # REL-01 (AUDIT_v2 §3 / §6 P0-1): rebuild `_CACHE` as a new dict
    # rather than mutating the existing dict in place. Pre-fix code did
    # `current = _load(); current.update(external)` — but `_load()`
    # returns a reference to the module-global `_CACHE`, so every prior
    # caller that still held a `_load()` reference (e.g. another test
    # function on the same xdist worker) saw the mutation, producing
    # documented order-dependent regression Th-232 chain ratio
    # 2.30× → 7.63× under `pytest -n auto`. The fixed contract is
    # snapshot-semantic: new callers see the merged library, prior
    # references remain an immutable snapshot of what they observed.
    global _CACHE
    # Ensure the cache is populated before we read it.
    _load()
    if merge_mode == "override":
        _CACHE = {**_CACHE, **external}
        return len(external)
    elif merge_mode == "supplement":
        new_cache = dict(_CACHE)
        added = 0
        for name, entry in external.items():
            if name not in new_cache:
                new_cache[name] = entry
                added += 1
        _CACHE = new_cache
        return added
    else:
        raise ValueError(f"Unknown merge_mode {merge_mode!r}; "
                         f"use 'override' or 'supplement'")


def list_nuclides() -> list[str]:
    """Return all nuclide names in the library."""
    return sorted(_load().keys())


def get_nuclide(name: str) -> Optional[dict]:
    """
    Return the full library record for one nuclide, or None if absent.
    The record has: T_half_s, lines, and optionally parent, daughters,
    chain, is_cascade, ic_xrays.
    """
    return _load().get(name)


def lookup_by_energy(
    E_keV: float,
    tolerance_keV: float,
    *,
    min_I_pct: float = 0.0,
    nuclides: Optional[list] = None,
) -> list:
    """
    Find all library lines within ±tolerance_keV of E_keV.

    Args:
        E_keV: target energy
        tolerance_keV: half-width of the search window
        min_I_pct: skip lines with intensity below this threshold (default 0)
        nuclides: restrict search to this subset (default: search all)

    Returns:
        List of NuclideLine, sorted by absolute |delta| ascending.
    """
    lib = _load()
    if nuclides is not None:
        items = [(n, lib[n]) for n in nuclides if n in lib]
    else:
        items = list(lib.items())

    hits = []
    for name, rec in items:
        for line in rec.get("lines", []):
            # line format: [E, I, sigma_I]
            line_E, line_I = line[0], line[1]
            line_sigma = line[2] if len(line) > 2 else 0.0
            if line_I < min_I_pct:
                continue
            delta = line_E - E_keV
            if abs(delta) <= tolerance_keV:
                hits.append(NuclideLine(
                    nuclide=name,
                    E_keV=line_E,
                    I_gamma_pct=line_I,
                    sigma_I_pct=line_sigma,
                    delta_keV=delta,
                ))

    hits.sort(key=lambda h: abs(h.delta_keV))
    return hits


def nuclides_in_chain(chain_name: str) -> list:
    """
    Return all nuclide names belonging to a decay chain.
    Chain names follow the parent: "U-238", "U-235", "Th-232".
    """
    lib = _load()
    return sorted(n for n, r in lib.items() if r.get("chain") == chain_name)


def parent_of(nuclide: str) -> Optional[str]:
    """Return parent nuclide name, or None if root / not in library."""
    rec = get_nuclide(nuclide)
    return rec.get("parent") if rec else None


def daughters_of(nuclide: str) -> list:
    """Return immediate-daughter list, or empty list."""
    rec = get_nuclide(nuclide)
    return rec.get("daughters", []) if rec else []


def is_cascade(nuclide: str) -> bool:
    """True if the nuclide is a known cascade emitter (TCS-affected)."""
    rec = get_nuclide(nuclide)
    return bool(rec.get("is_cascade", False)) if rec else False


def characteristic_line(nuclide: str) -> Optional[tuple]:
    """
    Return the line with highest I_gamma as a proxy for the characteristic
    line. (True 'characteristic line' = lowest MDA, which depends on
    background and detector — but for a small candidate list this proxy
    is usually adequate. The full MDA-based choice is made by the
    identification module at runtime.)

    Returns (E_keV, I_gamma_pct, sigma_I_pct) or None.
    """
    rec = get_nuclide(nuclide)
    if rec is None or not rec.get("lines"):
        return None
    best = max(rec["lines"], key=lambda l: l[1])
    return (best[0], best[1], best[2] if len(best) > 2 else 0.0)


def decay_correct(
    rate: float, nuclide: str, delta_t_s: float
) -> Optional[float]:
    """
    Decay-correct a count rate from time t back to t_0 = t - delta_t_s.

    rate_0 = rate(t) * exp(ln2 * delta_t / T_half)

    Returns None if T_half is unknown for this nuclide.
    """
    rec = get_nuclide(nuclide)
    if rec is None:
        return None
    T = rec.get("T_half_s")
    if T is None or T == 0:
        return None
    if math.isinf(T):
        return rate  # stable: no decay
    return rate * math.exp(math.log(2) * delta_t_s / T)


# ---------------------------------------------------------------------------
# Bundled Lsrm chain libraries (opt-in)
# ---------------------------------------------------------------------------
#
# F-39 (v1.7.17): Two Lsrm-native libraries shipped in
# `detectors/Gamma-1S/lsrm-libraries/` extend the built-in 27-nuclide JSON to
# 47+ nuclides on demand:
#
#   NaI-Etl+Esc.lib  - NaI-tuned representation of Th-232, Cs-137, K-40,
#                      Ra-226 chains. With split_chains=True it decomposes
#                      Th-232 into Ac-228, Pb-212, Bi-212, Tl-208,
#                      Ra-224, Th-228 -- exactly the daughter set needed
#                      to measure Th-228 cert sources via secular-
#                      equilibrium γ-emission.
#
#   ОСГИ.lib         - 33 ОСГИ certified source nuclides incl. Eu-154,
#                      Eu-155, Ce-144, Sn-113, Hg-203, Rh-106, Sb-125,
#                      Ir-192, Co-56, Ag-110m, Ta-182, Cs-134, Ru-103,
#                      Zr-95+, Ho-166m, plus chain members (Th-231,
#                      Th-234, U-232, TI-44).
#
# `load_lsrm_chain_libs()` is **opt-in** -- the default `_load()` path
# uses only `data/nuclides.json` (deterministic behaviour, K-03). Call
# this once at session/process startup to supplement the library:
#
#     from gamma.data.nuclide_library import load_lsrm_chain_libs
#     load_lsrm_chain_libs()
#     # ... now Tl-208 has full 5-line record, etc.

# F-83 (v1.12.0): the LSRM nuclide-library bundle (.lib files shipped by Lsrm
# SpectraLine) is a Gamma-1S-specific asset — see detectors/Gamma-1S/README.md.
# When other detector subtrees are added (deferred), each ships its own
# .lib bundle and points the loader at its own path resolver.
from gamma.detectors.gamma1s import LSRM_LIBRARIES_DIR as _LSRM_LIB_DIR
_NAI_CHAIN_LIB = _LSRM_LIB_DIR / "NaI-Etl+Esc.lib"
_OSGI_LIB = _LSRM_LIB_DIR / "ОСГИ.lib"


def load_lsrm_chain_libs(
    *,
    include_nai_chain: bool = True,
    include_osgi: bool = True,
    merge_mode: str = "supplement",
    split_chains: bool = True,
) -> dict:
    """
    Opt-in supplemental loader for the two bundled Lsrm libraries.

    Args:
        include_nai_chain: load NaI-Etl+Esc.lib (Th-232, Cs-137, K-40,
            Ra-226 chain decomposition tuned for NaI detectors).
        include_osgi: load ОСГИ.lib (33 ОСГИ certified-source nuclides
            with their full intensity + d_intensity records).
        merge_mode: "supplement" (default, keep existing entries) or
            "override" (Lsrm values win).
        split_chains: when True, chain-parent entries (Th-232, Ra-226,
            etc.) are decomposed into daughter-named records via
            `gamma.data.chain_decomposer`. Set False to keep the Lsrm
            sealed-source bundling.

    Returns:
        Dict {library_name: n_nuclides_added}.

    Effect:
        Mutates the in-memory cache. Subsequent `get_nuclide()`,
        `list_nuclides()`, `lookup_by_energy()` calls see the merged
        library. Call `reset_cache()` to revert to the built-in JSON
        only.
    """
    added: dict = {}
    if include_nai_chain and _NAI_CHAIN_LIB.is_file():
        n = load_external_library(
            str(_NAI_CHAIN_LIB),
            merge_mode=merge_mode,
            include_xrays=False,
            split_chains=split_chains,
        )
        added["NaI-Etl+Esc.lib"] = n
    if include_osgi and _OSGI_LIB.is_file():
        n = load_external_library(
            str(_OSGI_LIB),
            merge_mode=merge_mode,
            include_xrays=False,
            split_chains=split_chains,
        )
        added["ОСГИ.lib"] = n
    return added
