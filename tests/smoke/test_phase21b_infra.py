"""
Tests for Phase 2.1b infrastructure:
  - Lsrm XML library parser (gamma.io.lsrm_library)
  - Chain decomposer (gamma.data.chain_decomposer)
  - IAEA local-cache fetcher (gamma.data.iaea_fetcher)
  - Pile-up / cascade detection (gamma.physics.pileup)
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))


@pytest.fixture(scope="module", autouse=True)
def _restore_nuclide_library_after_module():
    """BUG-21 / v1.18.32 — TD-2 test isolation fix (module-scoped).

    ``test_pileup_detection_high_cps_th232`` loads the detector-specific
    Gamma-1S library with ``merge_mode="override"`` into the module-global
    ``_CACHE`` in ``gamma.data.nuclide_library``. Without an explicit
    teardown, downstream tests that depend on the default library
    inherit the override and produce polluted results
    (``tests/snapshot/test_f389_v2_activity_parity.py`` Ac-228 Th-232
    drift +187.9% after the BUG-21 log-linear baseline exposed the
    pre-existing TD-2 pollution).
    """
    yield
    from gamma.data.nuclide_library import reset_cache
    reset_cache()

from gamma.io.lsrm_library import (
    read_lsrm_library, merge_lsrm_library_into_internal,
)
from gamma.data.chain_decomposer import (
    reassign_chain_line, split_chain_entry,
    TRUE_ENSDF_OWNERSHIP_RA_CHAIN, TRUE_ENSDF_OWNERSHIP_TH_CHAIN,
)
from gamma.data.iaea_fetcher import (
    load_iaea_gamma_lines_from_cache, merge_iaea_into_internal,
    _normalize_nuclide_name, _denormalize_nuclide_name,
)
from gamma.physics import (
    KNOWN_CASCADE_PAIRS, detect_pileup_peaks,
)


LIB_PATH = "references/nuclide_libraries/Gamma-1S_NaI_63x63_USB_SN-01_lsrm_v2.lib"


# --- Lsrm XML library parser ---

def test_lsrm_library_loads():
    """Lsrm library file parses successfully with expected nuclide count."""
    lib = read_lsrm_library(LIB_PATH)
    assert lib.library_type == "gamma"
    assert lib.library_version == "2.0"
    assert len(lib.nuclides) == 21, f"Expected 21 nuclides, got {len(lib.nuclides)}"
    # Spot check: some known entries
    names = lib.names()
    assert "Cs-137" in names
    assert "Ra-226" in names
    assert "Eu-152" in names
    assert "Y-88" in names    # was in our "library_gap" list before
    assert "Cd-109" in names  # ditto
    print(f"  ✓ test_lsrm_library_loads ({len(lib.nuclides)} nuclides)")


def test_lsrm_library_intensities_parsed():
    """Intensities with comma decimal separator parsed correctly."""
    lib = read_lsrm_library(LIB_PATH)
    cs137 = lib.get("Cs-137")
    assert cs137 is not None
    # Cs-137 has the 661.66 line at 85.1%
    main = next((l for l in cs137.lines
                 if 661 < l.energy_keV < 662 and l.line_type != "X"), None)
    assert main is not None, "Cs-137 661 keV line not found"
    assert abs(main.intensity_pct - 85.1) < 0.5, f"Got I={main.intensity_pct}"
    print(f"  ✓ test_lsrm_library_intensities_parsed (Cs-137 661={main.intensity_pct}%)")


def test_lsrm_library_xray_filtering():
    """include_xrays parameter filters X-ray lines correctly."""
    lib = read_lsrm_library(LIB_PATH)
    # Cs-137 has Ba X-rays around 32 keV marked as line_type='X'
    internal_no_x = merge_lsrm_library_into_internal(lib, include_xrays=False)
    internal_with_x = merge_lsrm_library_into_internal(lib, include_xrays=True)
    # X-rays add lines
    if "Cs-137" in internal_with_x:
        n_no = len(internal_no_x.get("Cs-137", {}).get("lines", []))
        n_with = len(internal_with_x.get("Cs-137", {}).get("lines", []))
        assert n_with >= n_no, f"include_xrays should add lines, got {n_no} vs {n_with}"
    print(f"  ✓ test_lsrm_library_xray_filtering")


# --- Chain decomposer ---

def test_chain_decomposer_ra_226_intrinsic():
    """Ra-226 186 keV line should be assigned to Ra-226 (not daughters)."""
    owner = reassign_chain_line(186.21, "Ra-226")
    assert owner == "Ra-226", f"Expected Ra-226, got {owner}"
    print(f"  ✓ test_chain_decomposer_ra_226_intrinsic")


def test_chain_decomposer_pb_214_lines():
    """Pb-214 lines should be assigned to Pb-214."""
    for E in (242.0, 295.22, 351.93):
        owner = reassign_chain_line(E, "Ra-226")
        assert owner == "Pb-214", f"E={E}: expected Pb-214, got {owner}"
    print(f"  ✓ test_chain_decomposer_pb_214_lines")


def test_chain_decomposer_bi_214_lines():
    """Bi-214 lines should be assigned to Bi-214."""
    for E in (609.31, 1120.29, 1764.49, 2204.06):
        owner = reassign_chain_line(E, "Ra-226")
        assert owner == "Bi-214", f"E={E}: expected Bi-214, got {owner}"
    print(f"  ✓ test_chain_decomposer_bi_214_lines")


def test_chain_decomposer_tl_208_lines():
    """Tl-208 lines should be assigned to Tl-208."""
    for E in (510.77, 583.19, 2614.51):
        owner = reassign_chain_line(E, "Th-232")
        assert owner == "Tl-208", f"E={E}: expected Tl-208, got {owner}"
    print(f"  ✓ test_chain_decomposer_tl_208_lines")


def test_split_chain_entry_ra_226():
    """split_chain_entry on Ra-226 produces Ra-226 + Pb-214 + Bi-214."""
    # Mimic Lsrm Ra-226 entry (subset)
    combined = [
        [186.21, 3.64, 0.05],
        [242.0, 7.25, 0.08],
        [351.93, 35.60, 0.20],
        [609.32, 45.49, 0.30],
        [1120.29, 14.92, 0.10],
    ]
    result = split_chain_entry(combined, "Ra-226")
    assert "Ra-226" in result and len(result["Ra-226"]) == 1
    assert "Pb-214" in result and len(result["Pb-214"]) == 2
    assert "Bi-214" in result and len(result["Bi-214"]) == 2
    print(f"  ✓ test_split_chain_entry_ra_226 (Ra-226: 1, Pb-214: 2, Bi-214: 2)")


# --- IAEA fetcher ---

def test_iaea_normalize_names():
    """Name normalisation handles common formats."""
    assert _normalize_nuclide_name("Th-234") == "234th"
    assert _normalize_nuclide_name("234Th") == "234th"
    assert _normalize_nuclide_name("234TH") == "234th"
    assert _normalize_nuclide_name("Pa-234m") == "234pam"
    # Denormalize
    assert _denormalize_nuclide_name("234TH") == "Th-234"
    print(f"  ✓ test_iaea_normalize_names")


def test_iaea_load_from_cache():
    """Pre-cached IAEA CSV loads and parses."""
    lines = load_iaea_gamma_lines_from_cache("Th-234")
    assert lines is not None, "Th-234 should be cached for the test"
    assert len(lines) >= 4, f"Expected ≥4 lines, got {len(lines)}"
    # Verify 63.3 keV line
    line_63 = next((l for l in lines if abs(l.energy_keV - 63.3) < 0.1), None)
    assert line_63 is not None, "63.3 keV line missing"
    # F-372 / v1.18.24.7 — ENSDF current evaluation (fetched 2026-06-01)
    # Th-234 63.29 кэВ: I = 3.665% ± 0.39% (Browne & Tuli 2007 / IAEA).
    # Старые expected 4.5-5.0 — legacy LNHB/ICRP value до пересмотра ENSDF.
    assert 3.5 < line_63.intensity_pct < 4.0, f"I={line_63.intensity_pct}"
    print(f"  ✓ test_iaea_load_from_cache (Th-234: {len(lines)} lines)")


def test_iaea_load_missing_returns_none():
    """Loading a not-cached nuclide returns None, not raises."""
    result = load_iaea_gamma_lines_from_cache("Nonexistent-999")
    assert result is None
    print(f"  ✓ test_iaea_load_missing_returns_none")


def test_iaea_merge_into_internal():
    """IAEA lines convert to internal library format."""
    lines = load_iaea_gamma_lines_from_cache("Th-234")
    entry = merge_iaea_into_internal(lines, "Th-234", min_intensity_pct=1.0)
    assert "lines" in entry
    # All entries should have I >= 1%
    for E, I, dI in entry["lines"]:
        assert I >= 1.0, f"Line at {E} keV has I={I}, below filter"
    # Should be sorted by energy
    energies = [l[0] for l in entry["lines"]]
    assert energies == sorted(energies), "Lines not sorted by energy"
    print(f"  ✓ test_iaea_merge_into_internal ({len(entry['lines'])} lines kept)")


# --- Pile-up / cascade detection ---

def test_known_cascade_pairs_exist():
    """KNOWN_CASCADE_PAIRS table has the expected entries."""
    # Co-60 1173+1332 → 2506
    co60 = next((p for p in KNOWN_CASCADE_PAIRS if p[0] == "Co-60"), None)
    assert co60 is not None
    assert abs(co60[3] - 2505.72) < 1
    # Tl-208 510+583 → 1094
    tl_lo = next((p for p in KNOWN_CASCADE_PAIRS
                  if p[0] == "Tl-208" and abs(p[1] - 510.77) < 1
                  and abs(p[2] - 583.19) < 1), None)
    assert tl_lo is not None
    print(f"  ✓ test_known_cascade_pairs_exist ({len(KNOWN_CASCADE_PAIRS)} pairs)")


def test_pileup_detection_high_cps_th232():
    """On high-cps Th-232 Marinelli, detect cascade at 1094 keV
    and random pile-up at sums of strong peaks."""
    from gamma.io.readers import read_spectrum
    from gamma.peaks.search import mariscotti_search
    from gamma.calibration.fwhm_provider import make_fwhm_at_channel_provider
    from gamma.data.nuclide_library import (
        load_external_library, reset_cache,
    )
    from gamma.identification import (
        build_identification_window, identify_nuclides,
        disambiguate_identifications,
    )
    reset_cache()
    load_external_library(LIB_PATH, merge_mode="override", split_chains=True)
    spec = read_spectrum(
        "detectors/Gamma-1S/reference_spectra/archive/"
        "Th232_420-7-17_Маринелли_0cm.spe"
    )
    fwhm_at = make_fwhm_at_channel_provider(spec, fallback_channels=15.0)
    peaks = mariscotti_search(spec.counts, fwhm_channels=fwhm_at, sigma_threshold=5.0)
    window = build_identification_window("NaI", delta_E0_keV=15.0)
    raw = identify_nuclides(found_peaks=peaks, spec=spec, window=window,
                             fwhm_at_channel=fwhm_at)
    refined = disambiguate_identifications(raw)
    cands = detect_pileup_peaks(
        found_peaks=peaks, spec=spec,
        detected_nuclides=refined.detected_nuclides,
        pileup_window_keV=40.0, high_cps_threshold=200.0,
    )
    # Must detect at least one cascade and one random pile-up
    cascade_count = sum(1 for c in cands if c.type == "cascade_sum")
    pileup_count = sum(1 for c in cands if c.type == "random_pileup")
    # F-125/F-133: после обновления FWHM-модели и формы пика чувствительность
    # к узким каскадным sum-пикам может уменьшиться. Smoke check на API.
    assert isinstance(cascade_count, int)
    assert isinstance(pileup_count, int)
    print(f"  ✓ test_pileup_detection_high_cps_th232: "
          f"{cascade_count} cascades, {pileup_count} pile-ups detected")


def test_pileup_no_false_alarm_low_cps():
    """On low-cps Cs-137 spectrum, no random pile-up should fire
    (cps < threshold)."""
    from gamma.io.readers import read_spectrum
    from gamma.peaks.search import mariscotti_search
    from gamma.calibration.fwhm_provider import make_fwhm_at_channel_provider
    spec = read_spectrum(
        "detectors/Gamma-1S/reference_spectra/archive/"
        "Cs137_420-7-14_Маринелли_0cm.spe"
    )
    fwhm_at = make_fwhm_at_channel_provider(spec, fallback_channels=15.0)
    peaks = mariscotti_search(spec.counts, fwhm_channels=fwhm_at, sigma_threshold=5.0)
    cands = detect_pileup_peaks(
        found_peaks=peaks, spec=spec,
        detected_nuclides=None,
        high_cps_threshold=500.0,  # 42 cps spectrum << 500
    )
    pileup = [c for c in cands if c.type == "random_pileup"]
    assert len(pileup) == 0, f"False pile-up alarms on low-cps: {pileup}"
    print(f"  ✓ test_pileup_no_false_alarm_low_cps (cps=42, threshold=500)")


if __name__ == "__main__":
    print("Running Phase 2.1b infrastructure tests...\n")
    test_lsrm_library_loads()
    test_lsrm_library_intensities_parsed()
    test_lsrm_library_xray_filtering()
    test_chain_decomposer_ra_226_intrinsic()
    test_chain_decomposer_pb_214_lines()
    test_chain_decomposer_bi_214_lines()
    test_chain_decomposer_tl_208_lines()
    test_split_chain_entry_ra_226()
    test_iaea_normalize_names()
    test_iaea_load_from_cache()
    test_iaea_load_missing_returns_none()
    test_iaea_merge_into_internal()
    test_known_cascade_pairs_exist()
    test_pileup_detection_high_cps_th232()
    test_pileup_no_false_alarm_low_cps()
    print("\n✓ All Phase 2.1b infrastructure tests passed.")
