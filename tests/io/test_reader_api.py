"""
Symmetric energy-ceiling API for both file readers.

Both `read_atomspectra_xml` and `read_lsrm_spe` (and the shared
`read_spectrum` dispatcher) expose:
  - `apply_energy_ceiling: bool = False` — set True to drop channels
    above the ceiling. Default flipped to False in v1.18.32 (BUG-9):
    silently trimming a 1024-ch LSRM file to 1003 channels surprised
    users; the reader now keeps every decoded channel by default.
  - `ceiling_keV: float | None = None` — per-call override of the
    `ENERGY_CEILING_KEV` constant (3000 keV by project scope)

The embedded BackgroundEnergySpectrum (when present in a .xml file)
inherits the same trim policy as the primary spectrum so that the two
arrays stay channel-aligned.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

from gamma.io.readers import read_spectrum
from gamma.spectrum import ENERGY_CEILING_KEV


XML_BG = "detectors/AtomSpectra/data/fixtures/Фон_кабинет_8192к_01-01-2025.xml"
XML_WITH_EMBEDDED_BG = "detectors/AtomSpectra/data/fixtures/Cs137_0_см.xml"
SPE_CS137 = "evals/fixtures/M_cs_легкий_2001-2005.spe"


# ---------------------------------------------------------------------------
# AtomSpectra XML reader
# ---------------------------------------------------------------------------

def test_atomspectra_default_keeps_full_range():
    """
    Default call keeps every decoded channel (BUG-9 / v1.18.32: the
    reader no longer trims at 3 MeV by default; trim is opt-in).
    """
    s = read_spectrum(XML_BG)
    assert s.energy_max_keV_kept is not None
    expected_n = s.n_channels_raw - s.dropped_overflow_count
    assert s.n_channels == expected_n, (
        f"default call must keep full decoded range: n_channels {s.n_channels} "
        f"should equal n_channels_raw - dropped_overflow_count = {expected_n}"
    )
    assert s.extras["dropped_high_energy_count"] == 0, (
        "default call must not drop any high-energy channels"
    )
    assert s.energy_max_keV_kept > ENERGY_CEILING_KEV, (
        f"fixture should reach above 3000 keV when not trimmed "
        f"(got e_max={s.energy_max_keV_kept})"
    )
    print(f"  ✓ test_atomspectra_default_keeps_full_range "
          f"(n_ch={s.n_channels}, e_max={s.energy_max_keV_kept:.1f})")


def test_atomspectra_explicit_trim_at_3000():
    """`apply_energy_ceiling=True` opts back in to the 3 MeV trim."""
    s = read_spectrum(XML_BG, apply_energy_ceiling=True)
    assert s.energy_max_keV_kept is not None
    assert s.energy_max_keV_kept <= ENERGY_CEILING_KEV, (
        f"explicit-trim e_max {s.energy_max_keV_kept} should be <= "
        f"{ENERGY_CEILING_KEV}"
    )
    assert s.extras["dropped_high_energy_count"] > 0, (
        "fixture is expected to have channels above 3000 keV"
    )
    print(f"  ✓ test_atomspectra_explicit_trim_at_3000 "
          f"(n_ch={s.n_channels}, e_max={s.energy_max_keV_kept:.1f})")


def test_atomspectra_apply_false_keeps_full_range():
    """apply_energy_ceiling=False keeps every decoded channel (minus overflow)."""
    s = read_spectrum(XML_BG, apply_energy_ceiling=False)
    expected_n = s.n_channels_raw - s.dropped_overflow_count
    assert s.n_channels == expected_n, (
        f"with ceiling disabled, n_channels {s.n_channels} should equal "
        f"n_channels_raw - dropped_overflow_count = {expected_n}"
    )
    assert s.extras["dropped_high_energy_count"] == 0
    assert s.energy_max_keV_kept is not None
    assert s.energy_max_keV_kept > ENERGY_CEILING_KEV, (
        f"fixture should reach above 3000 keV when not trimmed "
        f"(got e_max={s.energy_max_keV_kept})"
    )
    print(f"  ✓ test_atomspectra_apply_false_keeps_full_range "
          f"(n_ch={s.n_channels}, e_max={s.energy_max_keV_kept:.1f})")


def test_atomspectra_custom_ceiling():
    """ceiling_keV=1500 produces a strictly shorter array than default."""
    default = read_spectrum(XML_BG)
    # `ceiling_keV` alone has no effect now that the default toggle is
    # False; the caller must opt in to trimming via apply_energy_ceiling.
    custom = read_spectrum(
        XML_BG, apply_energy_ceiling=True, ceiling_keV=1500.0
    )
    assert custom.energy_max_keV_kept is not None
    assert custom.energy_max_keV_kept <= 1500.0, (
        f"custom ceiling e_max {custom.energy_max_keV_kept} should be <= 1500"
    )
    assert custom.n_channels < default.n_channels, (
        f"custom ceiling should keep fewer channels: "
        f"{custom.n_channels} vs default {default.n_channels}"
    )
    print(f"  ✓ test_atomspectra_custom_ceiling "
          f"(custom n_ch={custom.n_channels} e_max={custom.energy_max_keV_kept:.1f} "
          f"vs default n_ch={default.n_channels})")


def test_atomspectra_embedded_bg_inherits():
    """Embedded BackgroundEnergySpectrum inherits the trim policy."""
    # apply_energy_ceiling=False — both primary and embedded keep full range
    s = read_spectrum(XML_WITH_EMBEDDED_BG, apply_energy_ceiling=False)
    assert s.background_embedded is not None, (
        f"fixture {XML_WITH_EMBEDDED_BG} expected to have embedded background"
    )
    p_full = s.n_channels_raw - s.dropped_overflow_count
    b_full = (s.background_embedded.n_channels_raw
              - s.background_embedded.dropped_overflow_count)
    assert s.n_channels == p_full
    assert s.background_embedded.n_channels == b_full
    # ceiling_keV=1500 — both shrink to the same ceiling (channel counts may
    # differ because primary and background have independent calibrations,
    # but each must satisfy e_max <= 1500). Opt-in toggle required since
    # v1.18.32 (BUG-9).
    s2 = read_spectrum(
        XML_WITH_EMBEDDED_BG, apply_energy_ceiling=True, ceiling_keV=1500.0
    )
    assert s2.background_embedded is not None
    assert s2.energy_max_keV_kept is not None and s2.energy_max_keV_kept <= 1500.0
    assert (s2.background_embedded.energy_max_keV_kept is not None
            and s2.background_embedded.energy_max_keV_kept <= 1500.0)
    print(f"  ✓ test_atomspectra_embedded_bg_inherits "
          f"(apply=False: primary={s.n_channels}, bg={s.background_embedded.n_channels}; "
          f"ceiling=1500: primary={s2.n_channels}, bg={s2.background_embedded.n_channels})")


# ---------------------------------------------------------------------------
# Lsrm SpectraLine .spe reader
# ---------------------------------------------------------------------------

def test_lsrm_spe_default_keeps_full_range():
    """
    Default .spe call keeps every decoded channel — no implicit 3 MeV
    trim (BUG-9 / v1.18.32). The reproducer for the user-reported
    «В файлах лсрм 1024 канала, в процессе конвертации идет потеря»
    failure: a 1024-ch LSRM file with a0≈-8 keV used to silently shrink
    to 1003 channels under the old default.
    """
    s = read_spectrum(SPE_CS137)
    assert s.energy_max_keV_kept is not None
    assert s.n_channels == s.n_channels_raw, (
        f"default LSRM read must keep full range: "
        f"n_channels={s.n_channels} vs raw={s.n_channels_raw}"
    )
    assert s.energy_max_keV_kept > ENERGY_CEILING_KEV, (
        f"fixture should reach above 3000 keV when not trimmed "
        f"(got e_max={s.energy_max_keV_kept})"
    )
    print(f"  ✓ test_lsrm_spe_default_keeps_full_range "
          f"(n_ch={s.n_channels}, raw={s.n_channels_raw}, "
          f"e_max={s.energy_max_keV_kept:.1f})")


def test_lsrm_spe_explicit_trim_at_3000():
    """`apply_energy_ceiling=True` opts back in to the 3 MeV trim."""
    s = read_spectrum(SPE_CS137, apply_energy_ceiling=True)
    assert s.energy_max_keV_kept is not None
    assert s.energy_max_keV_kept <= ENERGY_CEILING_KEV
    assert s.n_channels < s.n_channels_raw, (
        f"fixture is expected to have channels above 3000 keV "
        f"(n_ch={s.n_channels}, raw={s.n_channels_raw})"
    )
    print(f"  ✓ test_lsrm_spe_explicit_trim_at_3000 "
          f"(n_ch={s.n_channels}, raw={s.n_channels_raw}, "
          f"e_max={s.energy_max_keV_kept:.1f})")


def test_lsrm_spe_apply_false_keeps_full_range():
    """apply_energy_ceiling=False keeps every decoded channel."""
    s = read_spectrum(SPE_CS137, apply_energy_ceiling=False)
    assert s.n_channels == s.n_channels_raw, (
        f"with ceiling disabled, n_channels {s.n_channels} should equal "
        f"raw {s.n_channels_raw}"
    )
    assert s.energy_max_keV_kept is not None
    assert s.energy_max_keV_kept > ENERGY_CEILING_KEV
    print(f"  ✓ test_lsrm_spe_apply_false_keeps_full_range "
          f"(n_ch={s.n_channels}, e_max={s.energy_max_keV_kept:.1f})")


def test_lsrm_spe_custom_ceiling():
    """ceiling_keV=400 produces a strictly shorter .spe array than default."""
    default = read_spectrum(SPE_CS137)
    custom = read_spectrum(
        SPE_CS137, apply_energy_ceiling=True, ceiling_keV=400.0
    )
    assert custom.energy_max_keV_kept is not None
    assert custom.energy_max_keV_kept <= 400.0
    assert custom.n_channels < default.n_channels
    print(f"  ✓ test_lsrm_spe_custom_ceiling "
          f"(custom n_ch={custom.n_channels} e_max={custom.energy_max_keV_kept:.1f} "
          f"vs default n_ch={default.n_channels})")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("Running reader API symmetry tests...\n")
    test_atomspectra_default_keeps_full_range()
    test_atomspectra_explicit_trim_at_3000()
    test_atomspectra_apply_false_keeps_full_range()
    test_atomspectra_custom_ceiling()
    test_atomspectra_embedded_bg_inherits()
    test_lsrm_spe_default_keeps_full_range()
    test_lsrm_spe_explicit_trim_at_3000()
    test_lsrm_spe_apply_false_keeps_full_range()
    test_lsrm_spe_custom_ceiling()
    print("\nAll 9 reader API tests passed.")
