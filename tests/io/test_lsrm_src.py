"""
Tests for Phase 2.1f .src certificate parser (F-30, v1.7.8):
  • All 7 reference .src files parse without error
  • Decimal-comma values and DD.MM.YYYY dates parse correctly
  • Activity entries split as <value>,<sigma_pct> at the LAST comma
  • Sigma=N confidence level honoured in conversion to 1σ
  • Source lookup: exact, fuzzy substring, candidates
  • Multi-nuclide and multi-sub-source configurations parse correctly
  • End-to-end certificate validation against F-29 activity result
"""

from __future__ import annotations

import math
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.io.lsrm_src import (
    Certificate,
    CertificateActivity,
    CertificateSource,
    CertificateSubSource,
    read_certificate_file,
    read_certificate_files,
    find_certificate_for_nuclide,
)


CERT_DIR = Path("detectors/Gamma-1S/certificates")
OSGI_2024 = CERT_DIR / "АСПЕКТ_ОСГИ_2024.src"
OSGI_2024_2 = CERT_DIR / "Аспектр_ОСГИ_2024_2.src"
MARINELLI = CERT_DIR / "Эталон_Маринелли__Аспект_2017_.src"
PETRI = CERT_DIR / "Эталон_Петри-60__Аспект_2017_.src"
DENTA_100 = CERT_DIR / "Эталон_Дента100мл__Аспект2017_.src"
DENTA_120 = CERT_DIR / "Эталон_Дента120мл__Аспект2017_.src"
FILTER_AFA = CERT_DIR / "Фильтр_АФА__Аспект_2018_.src"


# ═════════════════════════════════════════════════════════════════════
# Group 1 — basic parsing
# ═════════════════════════════════════════════════════════════════════

def test_all_seven_files_parse():
    """All 7 reference .src files must parse without exception."""
    expected = [
        (OSGI_2024,    34, 40),     # 34 sources, 40 activities
        (OSGI_2024_2,  34, 40),
        (FILTER_AFA,    3,  5),
        (DENTA_100,     4, 13),
        (DENTA_120,    12, 27),
        (MARINELLI,     9, 24),
        (PETRI,         5, 17),
    ]
    for path, n_src, n_act in expected:
        cert = read_certificate_file(path)
        assert isinstance(cert, Certificate)
        n_src_actual = len(cert.sources)
        n_act_actual = sum(len(s.all_activities()) for s in cert.sources.values())
        assert n_src_actual == n_src, \
            f"{path.name}: {n_src_actual} sources vs expected {n_src}"
        assert n_act_actual == n_act, \
            f"{path.name}: {n_act_actual} activities vs expected {n_act}"
    print(f"  ✓ test_all_seven_files_parse (101 sources, 166 activities total)")


def test_cs137_163_2017_certificate():
    """The Cs-137 №SRC-02 source — known values verified externally
    against the spectrum's COMMENT field 'A=106000 Бк dA=3% 19-05-2017'.
    """
    cert = read_certificate_file(OSGI_2024)
    src = cert.find_source("Cs-137 №SRC-02")
    assert src is not None
    assert src.geometry == "Точечная"
    assert src.reference_datetime == datetime(2017, 5, 19)
    act = src.get_activity("Cs-137")
    assert act is not None
    assert math.isclose(act.A_Bq, 106000.0, rel_tol=1e-9)
    assert math.isclose(act.sigma_pct, 3.0, rel_tol=1e-9)
    assert act.confidence_sigma == 2
    assert act.unit == "Bq"
    # Conversion to 1σ
    assert math.isclose(act.sigma_1pct(), 1.5, rel_tol=1e-9)
    assert math.isclose(act.sigma_1_Bq(), 1590.0, rel_tol=1e-6)
    print(f"  ✓ test_cs137_163_2017_certificate "
          f"(A={act.A_Bq:.0f} Bq, σ_1={act.sigma_1pct()}%, "
          f"date={src.reference_datetime.date()})")


def test_general_sigma_default_when_missing():
    """Confidence_sigma defaults to 1 if [General] Sigma is missing."""
    # All real files have Sigma=2 — synthesise an in-memory test instead.
    import tempfile, os
    content = (
        "[Sets]\r\n"
        "Foo=\r\n"
        "[Foo]\r\n"
        "Geometry=Точечная\r\n"
        "Date=01.01.2020\r\n"
        "Time=0:00:00\r\n"
        "Units=Bq\r\n"
        "[Foo,structure]\r\n"
        "FooSub=\r\n"
        "[Foo,FooSub,Act]\r\n"
        "Cs-137=1000,5\r\n"
    )
    with tempfile.NamedTemporaryFile(
        suffix=".src", delete=False, mode="wb",
    ) as fp:
        fp.write(content.encode("cp1251"))
        tmp = fp.name
    try:
        cert = read_certificate_file(tmp)
        assert cert.confidence_sigma == 1
        act = cert.find_source("Foo").get_activity("Cs-137")
        # sigma_1pct == sigma_pct when confidence_sigma=1
        assert act.sigma_1pct() == act.sigma_pct
    finally:
        os.unlink(tmp)
    print(f"  ✓ test_general_sigma_default_when_missing (defaults to 1σ)")


def test_decimal_comma_in_metadata():
    """`Thick,mm=10,1` must parse as 10.1, not crash on the comma."""
    cert = read_certificate_file(PETRI)
    # Every Petri source has Thick,mm=10,1
    src = cert.find_source("Петри-60 ОИСН №420/7_р16")
    assert src is not None
    assert math.isclose(src.thickness_mm, 10.1, rel_tol=1e-9)
    print(f"  ✓ test_decimal_comma_in_metadata "
          f"(Thick,mm 10,1 → {src.thickness_mm})")


def test_activity_value_sigma_split():
    """`Cs-137=106000,3` → A=106000, σ_pct=3 (last comma is the splitter)."""
    cert = read_certificate_file(OSGI_2024)
    src = cert.find_source("Cs-137 №SRC-02")
    act = src.get_activity("Cs-137")
    assert act.A_Bq == 106000.0
    assert act.sigma_pct == 3.0
    print(f"  ✓ test_activity_value_sigma_split")


def test_units_bq_vs_bq_per_kg():
    """Point sources use 'Bq', volumetric standards use 'Bq/kg'."""
    cert_p = read_certificate_file(OSGI_2024)
    cs137 = cert_p.find_source("Cs-137 №SRC-02").get_activity("Cs-137")
    assert cs137.unit == "Bq"

    cert_v = read_certificate_file(PETRI)
    # Petri-60мл sources have Activity unit=Bq/kg (older "Activity unit" key)
    petri = cert_v.find_source("Петри-60 ОИСН №420/7_р16")
    cs137_petri = petri.get_activity("Cs-137")
    assert cs137_petri.unit == "Bq/kg", f"unit={cs137_petri.unit!r}"
    print(f"  ✓ test_units_bq_vs_bq_per_kg "
          f"(point={cs137.unit!r}, Petri={cs137_petri.unit!r})")


# ═════════════════════════════════════════════════════════════════════
# Group 2 — source lookup
# ═════════════════════════════════════════════════════════════════════

def test_find_source_exact_and_case_insensitive():
    """find_source uses normalised lower-case + whitespace strip."""
    cert = read_certificate_file(OSGI_2024)
    # Exact match
    assert cert.find_source("Cs-137 №SRC-02") is not None
    # Different casing
    assert cert.find_source("cs-137 №SRC-02") is not None
    # Surrounding whitespace
    assert cert.find_source("  Cs-137 №SRC-02  ") is not None
    # Non-existent
    assert cert.find_source("NoSuchSource") is None
    print(f"  ✓ test_find_source_exact_and_case_insensitive")


def test_find_source_fuzzy_strips_punctuation():
    """Fuzzy lookup ignores #, №, /, _ when matching."""
    cert = read_certificate_file(OSGI_2024)
    # Cert has "Cs-137 №SRC-02"; query without №
    src = cert.find_source_fuzzy("cs-137 SRC-02")
    assert src is not None
    assert "163" in src.name
    # Underscore / space variation: cert "Co-60 №043 02.2019"
    src2 = cert.find_source_fuzzy("Co-60 043 02 2019")
    assert src2 is not None
    print(f"  ✓ test_find_source_fuzzy_strips_punctuation "
          f"(matched {src.name!r} and {src2.name!r})")


def test_find_source_candidates_returns_all_matches():
    """For ambiguous queries (e.g. 'Co-60'), return all candidates."""
    cert = read_certificate_file(OSGI_2024)
    candidates = cert.find_source_candidates("Co-60")
    # OSGI 2024 contains: Co-60 №043 02.2019, Co-60 №189.2019,
    # Co-60 №479 09.2016, Co-60 #SRC-05 → ≥3 plausible matches
    assert len(candidates) >= 3, \
        f"Expected ≥3 Co-60 sources, got {len(candidates)}"
    all_names = [c.name for c in candidates]
    assert any("043" in n for n in all_names)
    print(f"  ✓ test_find_source_candidates_returns_all_matches "
          f"({len(candidates)} Co-60 candidates)")


def test_find_certificate_for_nuclide_across_files():
    """Convenience helper finds nuclide across multiple certificates."""
    paths = [OSGI_2024, MARINELLI, PETRI]
    certs = list(read_certificate_files(paths).values())
    # Cs-137 with hint "SRC-02" → must hit OSGI 2024's point source
    hit = find_certificate_for_nuclide(certs, "Cs-137",
                                        source_hint="SRC-02")
    assert hit is not None
    src, act = hit
    assert "163" in src.name
    assert act.nuclide == "Cs-137"
    # K-40 in any → must hit a Marinelli/Petri (point sources don't have K-40)
    hit2 = find_certificate_for_nuclide(certs, "K-40")
    assert hit2 is not None
    src2, act2 = hit2
    assert act2.nuclide == "K-40"
    print(f"  ✓ test_find_certificate_for_nuclide_across_files")


def test_get_activity_returns_none_for_missing_nuclide():
    """Querying a nuclide not in the source returns None."""
    cert = read_certificate_file(OSGI_2024)
    src = cert.find_source("Cs-137 №SRC-02")
    assert src.get_activity("Sr-90") is None
    assert src.get_activity("Pu-239") is None
    print(f"  ✓ test_get_activity_returns_none_for_missing_nuclide")


# ═════════════════════════════════════════════════════════════════════
# Group 3 — multi-nuclide and multi-sub-source structures
# ═════════════════════════════════════════════════════════════════════

def test_marinelli_multi_nuclide_in_single_act_block():
    """Маринелли05 ОИСН-42280 has 4 nuclides in one Act block."""
    cert = read_certificate_file(MARINELLI)
    src = cert.find_source("Маринелли05 ОИСН-42280.73585")
    assert src is not None
    activities = src.all_activities()
    nuclides = {a.nuclide for a in activities}
    assert nuclides == {"Cs-137", "Cd-109", "Am-241", "Eu-152"}, \
        f"got {nuclides}"
    cs = src.get_activity("Cs-137")
    assert cs.A_Bq == 8800.0 and cs.sigma_pct == 5.0
    print(f"  ✓ test_marinelli_multi_nuclide_in_single_act_block "
          f"({len(activities)} activities)")


def test_petri_multi_sub_source_structure():
    """Петри-60 ОИСН has 4 sub-sources (Cs/Th/Ra/K each in own Act block)."""
    cert = read_certificate_file(PETRI)
    src = cert.find_source("Петри-60 ОИСН №420/7_р16")
    assert src is not None
    assert len(src.sub_sources) == 4
    sub_names = [s.name for s in src.sub_sources]
    assert any("Cs137" in n for n in sub_names)
    assert any("K40" in n for n in sub_names)
    # Each sub-source carries exactly one nuclide
    for sub in src.sub_sources:
        assert len(sub.activities) == 1
    # Flattened all_activities should yield 4 entries
    assert len(src.all_activities()) == 4
    nuclides = {a.nuclide for a in src.all_activities()}
    assert nuclides == {"Cs-137", "Th-232", "Ra-226", "K-40"}
    print(f"  ✓ test_petri_multi_sub_source_structure "
          f"({len(src.sub_sources)} sub-sources, 4 nuclides)")


def test_geometry_diversity():
    """All seven distinct geometries appear in the reference set."""
    geometries = set()
    for path in (OSGI_2024, MARINELLI, PETRI, DENTA_100,
                 DENTA_120, FILTER_AFA):
        cert = read_certificate_file(path)
        for src in cert.sources.values():
            if src.geometry:
                geometries.add(src.geometry)
    # Expect at least the major variants
    expected = {"Точечная", "Маринелли", "Петри-60мл",
                "Дента-100мл", "Дента-120мл"}
    missing = expected - geometries
    assert not missing, f"Missing geometries: {missing}"
    print(f"  ✓ test_geometry_diversity ({len(geometries)} geometries)")


def test_repr_smoke():
    """All repr methods return strings containing key identifiers."""
    cert = read_certificate_file(OSGI_2024)
    src = cert.find_source("Cs-137 №SRC-02")
    act = src.get_activity("Cs-137")
    assert "Cs-137" in repr(act)
    assert "Bq" in repr(act)
    assert "163" in repr(src)
    assert "Точечная" in repr(src)
    assert "ОСГИ" in repr(cert) or "АСПЕКТ" in repr(cert)
    print(f"  ✓ test_repr_smoke")


# ═════════════════════════════════════════════════════════════════════
# Group 4 — end-to-end validation against F-29 activity calculation
# ═════════════════════════════════════════════════════════════════════

def test_e2e_cs137_certificate_validation():
    """End-to-end: spectrum + certificate → decay-corrected activity
    must agree with certificate within combined uncertainty.

    Pipeline: read .spe → calibrate → identify → ε(E) → activity
    with decay correction → compare to certificate.

    Cs-137 №SRC-02 source:
      • Certificate: A_ref = 106000 ± 3% (2σ) Bq @ 19.05.2017
      • Spectrum:    measured 21.10.2024 (Δt = 2712 days)
      • T½(Cs-137) = 30.05 yr → decay factor ≈ 1.187 (ref/meas)

    Acceptance criterion: |A_measured - A_cert| ≤ 2 × combined 1σ.
    """
    from gamma.io.readers import read_spectrum
    from gamma.calibration.efficiency import fit_efficiency_from_efr_file
    from gamma.calibration.fwhm_provider import make_fwhm_at_channel_provider
    from gamma.peaks.search import mariscotti_search
    from gamma.identification.window import build_identification_window
    from gamma.identification.identify import identify_nuclides
    from gamma.identification.disambiguate import disambiguate_identifications
    from gamma.activity import compute_activity

    spe_path = ("detectors/Gamma-1S/reference_spectra/"
                "archive/Cs-137__163_2017.spe")
    efr_path = ("detectors/Gamma-1S/efficiency/Gamma-1S_NaI_63x63_USB_SN-01/"
                "УДС-ГЦ-63х63-USB__SN-01_-_Точечная-5см.efr")

    spec = read_spectrum(spe_path)
    assert spec.start_datetime is not None, \
        "spectrum date must parse (F-30 fix to _parse_lsrm_datetime)"

    cert = read_certificate_file(OSGI_2024)
    src = cert.find_source("Cs-137 №SRC-02")
    act_cert = src.get_activity("Cs-137")

    eff = fit_efficiency_from_efr_file(efr_path, degree=3)
    fwhm_at = make_fwhm_at_channel_provider(spec, fallback_channels=15.0)
    peaks = mariscotti_search(
        spec.counts, fwhm_channels=fwhm_at, sigma_threshold=5.0,
    )
    window = build_identification_window("NaI", delta_E0_keV=15.0)
    raw_id = identify_nuclides(
        found_peaks=peaks, spec=spec, window=window,
        fwhm_at_channel=fwhm_at,
    )
    refined = disambiguate_identifications(raw_id)
    cs = next(ni for ni in refined.detected_nuclides if ni.nuclide == "Cs-137")

    result = compute_activity(
        cs, efficiency_curve=eff,
        live_time_s=spec.live_time,
        from_bg_subtracted=False, bg_available=False,
        decay_correction=True,
        measurement_datetime=spec.start_datetime,
        reference_datetime=src.reference_datetime,
    )
    assert result.decay_corrected
    # Δt ≈ 7.43 yr → factor = exp(ln2·7.43/30.05) ≈ 1.187
    assert 1.15 < result.decay_factor < 1.22, \
        f"decay factor {result.decay_factor:.3f} outside expected range"

    deviation_pct = (result.A_Bq - act_cert.A_Bq) / act_cert.A_Bq * 100.0
    combined_sigma_pct = math.sqrt(
        (result.sigma_A_Bq / result.A_Bq * 100.0) ** 2
        + act_cert.sigma_1pct() ** 2
    )
    # Within 2σ combined — this also implies within 1σ for the
    # observed case but the looser bound makes the test robust
    # against ε(E) curve refits etc.
    assert abs(deviation_pct) < 2.0 * combined_sigma_pct, (
        f"|Δ|={abs(deviation_pct):.2f}% exceeds 2σ_combined="
        f"{2*combined_sigma_pct:.2f}%; A_meas={result.A_Bq:.0f}, "
        f"A_cert={act_cert.A_Bq:.0f}"
    )
    print(f"  ✓ test_e2e_cs137_certificate_validation "
          f"(Δ={deviation_pct:+.2f}%, combined 1σ={combined_sigma_pct:.2f}%, "
          f"A_meas={result.A_Bq:.0f}, A_cert={act_cert.A_Bq:.0f} Bq)")


# ═════════════════════════════════════════════════════════════════════
# Runner
# ═════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Running Phase 2.1f .src certificate parser tests...\n")

    # Group 1 — basic parsing (6)
    test_all_seven_files_parse()
    test_cs137_163_2017_certificate()
    test_general_sigma_default_when_missing()
    test_decimal_comma_in_metadata()
    test_activity_value_sigma_split()
    test_units_bq_vs_bq_per_kg()

    # Group 2 — source lookup (5)
    test_find_source_exact_and_case_insensitive()
    test_find_source_fuzzy_strips_punctuation()
    test_find_source_candidates_returns_all_matches()
    test_find_certificate_for_nuclide_across_files()
    test_get_activity_returns_none_for_missing_nuclide()

    # Group 3 — multi-nuclide / multi-sub-source (4)
    test_marinelli_multi_nuclide_in_single_act_block()
    test_petri_multi_sub_source_structure()
    test_geometry_diversity()
    test_repr_smoke()

    # Group 4 — end-to-end certificate validation (1)
    test_e2e_cs137_certificate_validation()

    print("\n✓ All Phase 2.1f .src certificate parser tests passed.")
