"""BUG-36 (Wave 3): Ti-44 + Sc-44 must be present in the DEFAULT nuclide
library (`data/nuclides.json`), not only via opt-in
`load_lsrm_chain_libs(ОСГИ)`.

Without these defaults, AmTiCsEu fixture cannot identify Ti-44 by its
characteristic 67.87 keV / 78.34 keV lines, nor reach the daughter
Sc-44 1157.02 keV line that is used in equilibrium-based activity
measurement.

Source provenance per number (NNDC ENSDF, retrieved 2026-06-04 via
https://www.nndc.bnl.gov/ensdf/):

* Ti-44 half-life: 60.0 ± 1.1 years = 60.0 * 365.25 * 86400 s
  = 1.893456e9 s (NNDC ENSDF, A=44 dataset, parent Ti-44 EC).
* Ti-44 γ lines: 67.8679 keV (I_g = 92.95 ± 0.06 %), 78.337 keV
  (I_g = 96.4 ± 0.07 %), 146.222 keV (I_g = 0.091 ± 0.005 %)
  — NNDC ENSDF Ti-44 dataset, γ-emission table.
* Sc-44 (ground state) half-life: 3.97 ± 0.04 h = 3.97 * 3600
  = 1.4292e4 s.  NB: in radioactive equilibrium with Ti-44 parent
  (T½_parent >> T½_daughter ⇒ secular equilibrium A_Sc44 ≈ A_Ti44).
* Sc-44 γ line: 1157.02 keV (I_g = 99.9 ± 0.4 %) — NNDC ENSDF Sc-44
  dataset.  511 keV β+ annihilation is diagnostic only (not added as
  it is detector-dependent, mirroring exclusion of 511 keV in
  load_external_library/include_xrays semantics).

Equilibrium link is encoded with the existing canonical fields used
throughout the library (cf. Ra-226 → daughters: ["Pb-214", "Bi-214"]
and Pb-214 → parent: "Rn-222"):

* Ti-44.daughters = ["Sc-44"], Ti-44.chain = "Ti-44"
* Sc-44.parent    = "Ti-44",   Sc-44.chain = "Ti-44"

No new schema fields invented (per brief — only canonical
parent/daughters/chain used).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
NUCLIDES_PATH = REPO_ROOT / "data" / "nuclides.json"


@pytest.fixture(scope="module")
def nuclides() -> dict:
    return json.loads(NUCLIDES_PATH.read_text(encoding="utf-8"))


def test_ti44_in_default_library(nuclides: dict) -> None:
    """BUG-36: Ti-44 must appear in the default library, not only via
    the opt-in load_lsrm_chain_libs(ОСГИ) supplemental loader."""
    assert "Ti-44" in nuclides, (
        "Ti-44 missing from default data/nuclides.json — AmTiCsEu fixture "
        "identification will fail (Ti-44 unreachable without "
        "load_lsrm_chain_libs)"
    )
    entry = nuclides["Ti-44"]
    assert "T_half_s" in entry
    assert "lines" in entry and len(entry["lines"]) >= 1


def test_sc44_in_default_library(nuclides: dict) -> None:
    """Sc-44 daughter must appear so equilibrium-based activity
    measurement through 1157.02 keV is available."""
    assert "Sc-44" in nuclides, (
        "Sc-44 missing from default data/nuclides.json — Ti-44 activity "
        "via Sc-44 1157 keV daughter γ unreachable"
    )
    entry = nuclides["Sc-44"]
    assert "T_half_s" in entry
    assert "lines" in entry and len(entry["lines"]) >= 1


def test_ti44_sc44_equilibrium_link(nuclides: dict) -> None:
    """Parent/daughter cross-reference must be intact (canonical
    fields parent/daughters — see Ra-226 → Pb-214/Bi-214 pattern).
    This is the schema-compatible way to encode the secular-equilibrium
    relationship (T½_Ti44=60 y >> T½_Sc44=3.97 h).

    Note: the `chain` field was deliberately cleared in F-LIB-EXTENSION
    (commit 3fa32ea, 2026-06-07) — Ti-44 is not a NORM chain, so the
    canonical parent/daughters link is used instead of the chain shortcut.
    """
    ti44 = nuclides["Ti-44"]
    sc44 = nuclides["Sc-44"]
    assert "Sc-44" in ti44.get("daughters", []), (
        f"Ti-44.daughters must include Sc-44; got {ti44.get('daughters')!r}"
    )
    assert sc44.get("parent") == "Ti-44", (
        f"Sc-44.parent must be 'Ti-44'; got {sc44.get('parent')!r}"
    )


def test_ti44_char_line_67_keV(nuclides: dict) -> None:
    """Ti-44 characteristic γ-line at 67.87 keV (I≈93%) must be in
    lines list within 0.1 keV (NNDC ENSDF Ti-44 dataset)."""
    lines = nuclides["Ti-44"]["lines"]
    energies = [line[0] for line in lines]
    assert any(abs(E - 67.87) <= 0.1 for E in energies), (
        f"Ti-44 67.87 keV characteristic line missing; got energies={energies}"
    )


def test_sc44_char_line_1157_keV(nuclides: dict) -> None:
    """Sc-44 characteristic γ-line at 1157.02 keV (I≈99.9%) must be in
    lines list within 0.1 keV (NNDC ENSDF Sc-44 dataset). This is the
    line used in secular-equilibrium activity measurement of Ti-44."""
    lines = nuclides["Sc-44"]["lines"]
    energies = [line[0] for line in lines]
    assert any(abs(E - 1157.02) <= 0.1 for E in energies), (
        f"Sc-44 1157.02 keV characteristic line missing; "
        f"got energies={energies}"
    )


def test_ti44_half_life_60_years(nuclides: dict) -> None:
    """Sanity: Ti-44 T½ = 60.0 y = 1.893e9 s (NNDC ENSDF)."""
    T = nuclides["Ti-44"]["T_half_s"]
    expected = 60.0 * 365.25 * 86400.0  # = 1.893456e9
    rel_err = abs(T - expected) / expected
    assert rel_err < 0.05, (
        f"Ti-44 T½={T} s, expected ≈{expected:.3e} s (60 y); "
        f"rel_err={rel_err:.3f}"
    )


def test_sc44_half_life_3p97_hours(nuclides: dict) -> None:
    """Sanity: Sc-44 T½ = 3.97 h = 14292 s (NNDC ENSDF)."""
    T = nuclides["Sc-44"]["T_half_s"]
    expected = 3.97 * 3600.0  # = 14292
    rel_err = abs(T - expected) / expected
    assert rel_err < 0.05, (
        f"Sc-44 T½={T} s, expected ≈{expected:.0f} s (3.97 h); "
        f"rel_err={rel_err:.3f}"
    )
