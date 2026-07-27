"""F-441 - isolated-peak classifier (Rayleigh isolation) unit tests.

Per brief _state/agent_a/inbox/2026-06-13_F-441_isolated_peak_classifier.md §4.1.

5 unit tests covering the _is_isolated_peak classifier directly:
  1. Tl-208 583 keV with Ac-228 weak neighbour (0.93%) -> isolated
  2. Tl-208 2614 keV with no neighbours -> isolated
  3. Strong neighbour in window (>= 3%) -> NOT isolated
  4. Weak neighbour in window (< 3%) -> isolated
  5. Strong neighbour OUTSIDE window -> isolated (window respect)

Defaults under test:
  window_fwhm = 1.0 (Rayleigh criterion)
  min_neighbor_I_pct = 3.0
  self_match_keV = 0.1
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

from gamma.peaks.deconvolve import (
    _f441_flatten_library_lines,
    _f441_fwhm_keV_at_channel,
    _is_isolated_peak,
)


# ---------------------------------------------------------------------------
# Minimal fakes (no real spectra / identification - unit level).
# ---------------------------------------------------------------------------

@dataclass
class _FakeLineMatch:
    nuclide: str
    library_E_keV: float
    peak_channel: float
    library_I_pct: float = 0.0
    peak_area: float = 100.0
    peak_area_uncertainty: float = 10.0
    peak_area_source: str = "cowell"


class _LinearSpec:
    """Channel<->energy at 1.0 keV/channel for predictable FWHM_keV math."""

    def channel_to_energy(self, ch):
        return float(ch)

    def energy_to_channel(self, E):
        return float(E)


def _fwhm_naI_at(ch):
    """NaI 63x63 - approximate FWHM(channels) at channel ch (1 keV/ch).

    Empirical (project default v1.31.3): ~47 keV at 583, ~64 keV at 1461,
    ~95 keV at 2614. Use a simple linear fit through (583, 47) and
    (2614, 95): slope = 48/2031 ~= 0.02364, intercept ~= 33.2.
    """
    return 33.2 + 0.02364 * float(ch)


# ---------------------------------------------------------------------------
# 1. Tl-208 583.19 keV - Ac-228 weak neighbour 562 at 0.93% -> isolated
# ---------------------------------------------------------------------------

def test_tl208_583_isolated_against_weak_ac228_neighbour():
    """Tl-208 583 keV: Ac-228 562 (0.93%) and 572 (0.18%) within ~47 keV
    FWHM window are both BELOW 3% threshold -> classifier returns True.

    Reference: brief §4.1 case 1.
    """
    m = _FakeLineMatch(
        nuclide="Tl-208", library_E_keV=583.19, peak_channel=583.0,
        library_I_pct=84.5,
    )
    spec = _LinearSpec()
    all_lib_lines = [
        ("Tl-208", 583.19, 84.5),
        ("Ac-228", 562.50, 0.93),
        ("Ac-228", 572.30, 0.18),
        ("Ac-228", 463.00, 4.40),  # outside window
        ("Ac-228", 794.95, 4.34),  # outside window
    ]
    assert _is_isolated_peak(m, _fwhm_naI_at, spec, all_lib_lines) is True


# ---------------------------------------------------------------------------
# 2. Tl-208 2614.51 keV - no library neighbour in window -> isolated
# ---------------------------------------------------------------------------

def test_tl208_2614_isolated_no_neighbour():
    """Tl-208 2614 keV: empty +/-FWHM window -> isolated.

    Reference: brief §4.1 case 2.
    """
    m = _FakeLineMatch(
        nuclide="Tl-208", library_E_keV=2614.51, peak_channel=2614.0,
        library_I_pct=99.0,
    )
    spec = _LinearSpec()
    all_lib_lines = [
        ("Tl-208", 2614.51, 99.0),
        ("Tl-208", 583.19, 84.5),    # far outside
        ("Ac-228", 911.20, 25.8),    # far outside
    ]
    assert _is_isolated_peak(m, _fwhm_naI_at, spec, all_lib_lines) is True


# ---------------------------------------------------------------------------
# 3. Strong neighbour in window -> NOT isolated
# ---------------------------------------------------------------------------

def test_strong_neighbour_in_window_blocks_isolation():
    """A hypothetical 5%-intensity neighbour within +/-FWHM window
    blocks isolation - classifier returns False.

    Reference: brief §4.1 case 3 (synthetic strong-neighbour stress test).
    """
    m = _FakeLineMatch(
        nuclide="Tl-208", library_E_keV=583.19, peak_channel=583.0,
        library_I_pct=84.5,
    )
    spec = _LinearSpec()
    # FWHM at 583 ~= 47 keV - place neighbour at 510 (delta = 73 keV)
    # is OUTSIDE 1*FWHM. Place at 605 (delta = 22 keV) is INSIDE window.
    all_lib_lines = [
        ("Tl-208", 583.19, 84.5),
        ("FAKE", 605.00, 5.0),  # strong (>=3%) AND inside ~47 keV window
    ]
    assert _is_isolated_peak(m, _fwhm_naI_at, spec, all_lib_lines) is False


# ---------------------------------------------------------------------------
# 4. Weak neighbour in window -> isolated (window allows weak neighbours)
# ---------------------------------------------------------------------------

def test_weak_neighbour_in_window_allows_isolation():
    """A 0.5%-intensity neighbour inside the window is BELOW the 3%
    threshold and does not block isolation -> classifier returns True.

    Reference: brief §4.1 case 4.
    """
    m = _FakeLineMatch(
        nuclide="Tl-208", library_E_keV=583.19, peak_channel=583.0,
        library_I_pct=84.5,
    )
    spec = _LinearSpec()
    all_lib_lines = [
        ("Tl-208", 583.19, 84.5),
        ("FAKE", 605.00, 0.5),  # weak (<3%), inside window
    ]
    assert _is_isolated_peak(m, _fwhm_naI_at, spec, all_lib_lines) is True


# ---------------------------------------------------------------------------
# 5. Strong neighbour OUTSIDE window -> isolated (window respect)
# ---------------------------------------------------------------------------

def test_strong_neighbour_outside_window_does_not_block():
    """A 50% neighbour FAR outside +/-1*FWHM window must NOT block
    isolation. Confirms the window bound is actually respected.

    Reference: brief §4.1 case 5.
    """
    m = _FakeLineMatch(
        nuclide="Tl-208", library_E_keV=583.19, peak_channel=583.0,
        library_I_pct=84.5,
    )
    spec = _LinearSpec()
    # FWHM ~= 47 keV at 583. Place strong neighbour at 700 (delta = 117 keV)
    # comfortably outside +/-47 keV window.
    all_lib_lines = [
        ("Tl-208", 583.19, 84.5),
        ("FAKE", 700.00, 50.0),  # very strong but OUTSIDE window
    ]
    assert _is_isolated_peak(m, _fwhm_naI_at, spec, all_lib_lines) is True


# ---------------------------------------------------------------------------
# Defensive coverage: empty pool, degenerate spec, self-match exclusion
# ---------------------------------------------------------------------------

def test_empty_library_pool_returns_isolated():
    """No detected nuclides -> empty pool -> isolated (vacuously true)."""
    m = _FakeLineMatch(
        nuclide="Tl-208", library_E_keV=583.19, peak_channel=583.0,
        library_I_pct=84.5,
    )
    spec = _LinearSpec()
    assert _is_isolated_peak(m, _fwhm_naI_at, spec, []) is True


def test_degenerate_fwhm_returns_not_isolated_conservatively():
    """If FWHM cannot be evaluated (returns 0), classifier returns False
    so pre-F-441 routing kicks in. Defensive contract.
    """
    def _zero_fwhm(_ch):
        return 0.0

    m = _FakeLineMatch(
        nuclide="Tl-208", library_E_keV=583.19, peak_channel=583.0,
        library_I_pct=84.5,
    )
    spec = _LinearSpec()
    all_lib_lines = [("Tl-208", 583.19, 84.5)]
    assert _is_isolated_peak(m, _zero_fwhm, spec, all_lib_lines) is False


def test_self_match_excluded_from_neighbour_test():
    """The line itself (dE < self_match_keV) must not block isolation,
    even with very high library intensity. Otherwise every line would
    fail the isolation test against itself.
    """
    m = _FakeLineMatch(
        nuclide="Tl-208", library_E_keV=583.19, peak_channel=583.0,
        library_I_pct=84.5,
    )
    spec = _LinearSpec()
    # Pool contains the same line at 583.19 with very high I_gamma -
    # must be skipped by self-match guard.
    all_lib_lines = [("Tl-208", 583.19, 84.5)]
    assert _is_isolated_peak(m, _fwhm_naI_at, spec, all_lib_lines) is True


# ---------------------------------------------------------------------------
# _f441_flatten_library_lines coverage
# ---------------------------------------------------------------------------

def test_flatten_library_lines_basic():
    """Verify flatten produces (nuclide, E, I) tuples for declared nuclides."""
    lib = {
        "Tl-208": {"lines": [(583.19, 84.5), (2614.51, 99.0)]},
        "Ac-228": {"lines": [(911.20, 25.8)]},
    }
    flat = _f441_flatten_library_lines(lib, ["Tl-208", "Ac-228"])
    assert len(flat) == 3
    e_keVs = sorted(t[1] for t in flat)
    assert e_keVs == [583.19, 911.20, 2614.51]


def test_flatten_library_lines_only_detected():
    """Lines of NON-detected nuclides must be excluded from the pool."""
    lib = {
        "Tl-208": {"lines": [(583.19, 84.5)]},
        "Cs-137": {"lines": [(661.66, 85.1)]},
    }
    flat = _f441_flatten_library_lines(lib, ["Tl-208"])  # Cs-137 NOT detected
    assert len(flat) == 1
    assert flat[0][0] == "Tl-208"


def test_flatten_library_lines_handles_dict_records():
    """Both list/tuple and dict line records are accepted by the flattener."""
    lib = {
        "Tl-208": {"lines": [
            {"E_keV": 583.19, "I_pct": 84.5},
            {"E_keV": 2614.51, "I_pct": 99.0},
        ]},
    }
    flat = _f441_flatten_library_lines(lib, ["Tl-208"])
    assert len(flat) == 2
    e_keVs = sorted(t[1] for t in flat)
    assert e_keVs == [583.19, 2614.51]


def test_fwhm_keV_at_channel_basic():
    """Linear spec at 1 keV/channel: FWHM_keV equals FWHM_channels exactly."""
    spec = _LinearSpec()

    def _const_fwhm_ch(_ch):
        return 47.0

    keV = _f441_fwhm_keV_at_channel(583.0, _const_fwhm_ch, spec)
    assert abs(keV - 47.0) < 1e-9


def test_fwhm_keV_at_channel_zero_fwhm_returns_zero():
    """Degenerate FWHM(channels) = 0 -> returns 0.0."""
    spec = _LinearSpec()
    assert _f441_fwhm_keV_at_channel(500.0, lambda c: 0.0, spec) == 0.0


def test_fwhm_keV_at_channel_none_spec_returns_zero():
    """None spec -> returns 0.0 (defensive)."""
    assert _f441_fwhm_keV_at_channel(500.0, lambda c: 47.0, None) == 0.0