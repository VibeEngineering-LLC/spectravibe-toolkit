# -*- coding: utf-8 -*-
"""v1.17.21 — Peak-image .cpt layer (F-299..F-301)."""
from __future__ import annotations
import math, os, sys
from pathlib import Path
import pytest

SCRIPTS = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts")
)
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)


# ──────────────────────────────────────────────────────────────────
# F-299 — Tabulated peak image
# ──────────────────────────────────────────────────────────────────

def test_F299_build_from_calibration_pairs_nai():
    from gamma.peaks.peak_image_tabulated import (
        build_anchors_from_calibration,
    )
    img = build_anchors_from_calibration(
        detector_id="Gamma-1S",
        detector_class="NaI",
        crystal_diameter_mm=63.0,
        calibration_pairs=[
            (122.0, 12.5),
            (662.0, 46.5),
            (1332.0, 78.0),
        ],
    )
    assert len(img.anchors) == 3
    assert img.detector_class == "NaI"
    # NaI defaults applied
    assert img.anchors[0].tail_fraction == pytest.approx(0.03)
    assert img.anchors[0].step_height_frac == pytest.approx(0.05)


def test_F299_validate_catches_negative_fwhm():
    from gamma.peaks.peak_image_tabulated import (
        TabulatedPeakImage, PeakShapeAnchor,
    )
    img = TabulatedPeakImage(
        detector_id="x", detector_class="NaI", crystal_diameter_mm=63.0,
        anchors=[PeakShapeAnchor(E_keV=662.0, fwhm_keV=-1.0)],
    )
    issues = img.validate()
    assert any("FWHM" in i and "≤ 0" in i for i in issues)


def test_F299_validate_nai_resolution_check():
    from gamma.peaks.peak_image_tabulated import (
        TabulatedPeakImage, PeakShapeAnchor,
    )
    # 30% FWHM @ 662 keV для NaI — нефизично
    img = TabulatedPeakImage(
        detector_id="x", detector_class="NaI", crystal_diameter_mm=63.0,
        anchors=[PeakShapeAnchor(E_keV=662.0, fwhm_keV=200.0)],
    )
    issues = img.validate()
    assert any("FWHM%" in i for i in issues)


def test_F299_anchor_at_E_finds_within_tolerance():
    from gamma.peaks.peak_image_tabulated import (
        build_anchors_from_calibration,
    )
    img = build_anchors_from_calibration(
        "Gamma-1S", "NaI", 63.0,
        [(122.0, 12.5), (662.0, 46.5)],
    )
    found = img.anchor_at_E(660.0, tolerance_keV=5.0)
    assert found is not None
    assert found.E_keV == 662.0


def test_F299_estimate_fwhm_pct_at_662():
    from gamma.peaks.peak_image_tabulated import (
        build_anchors_from_calibration,
    )
    img = build_anchors_from_calibration(
        "Gamma-1S", "NaI", 63.0,
        [(662.0, 46.34)],
    )
    pct = img.estimate_fwhm_pct_at_662()
    assert pct == pytest.approx(7.0, abs=0.1)


def test_F299_hpge_defaults_different_from_nai():
    from gamma.peaks.peak_image_tabulated import (
        build_anchors_from_calibration,
    )
    nai = build_anchors_from_calibration(
        "n", "NaI", 63.0, [(662.0, 46.5)],
    )
    hpge = build_anchors_from_calibration(
        "h", "HPGe", 76.0, [(662.0, 1.3)],
    )
    assert nai.anchors[0].tail_fraction > hpge.anchors[0].tail_fraction
    assert nai.anchors[0].step_height_frac > hpge.anchors[0].step_height_frac


# ──────────────────────────────────────────────────────────────────
# F-300 — Log-spline interpolation
# ──────────────────────────────────────────────────────────────────

def test_F300_exact_anchor_returns_exact():
    from gamma.peaks.peak_image_tabulated import (
        build_anchors_from_calibration,
    )
    from gamma.peaks.peak_image_logspline import interpolate_peak_shape
    img = build_anchors_from_calibration(
        "Gamma-1S", "NaI", 63.0,
        [(122.0, 12.5), (662.0, 46.5), (1332.0, 78.0)],
    )
    res = interpolate_peak_shape(img, 662.0, exact_match_tolerance_keV=1.0)
    assert res.fwhm_keV == pytest.approx(46.5, abs=1e-6)
    assert not res.was_extrapolated


def test_F300_interpolate_between_anchors_returns_intermediate():
    from gamma.peaks.peak_image_tabulated import (
        build_anchors_from_calibration,
    )
    from gamma.peaks.peak_image_logspline import interpolate_peak_shape
    img = build_anchors_from_calibration(
        "Gamma-1S", "NaI", 63.0,
        [(100.0, 10.0), (1000.0, 60.0)],
    )
    res = interpolate_peak_shape(img, 316.23)   # log-mid: √(100·1000)
    # In log-log space, log(316.23)≈log(100)+0.5*(log(1000)-log(100))
    # So fwhm should be exp(log(10)+0.5*(log(60)-log(10))) = √(10·60) = √600 ≈ 24.5
    expected = math.sqrt(10.0 * 60.0)
    assert res.fwhm_keV == pytest.approx(expected, rel=0.02)
    assert not res.was_extrapolated


def test_F300_extrapolation_below_min_flagged():
    from gamma.peaks.peak_image_tabulated import (
        build_anchors_from_calibration,
    )
    from gamma.peaks.peak_image_logspline import interpolate_peak_shape
    img = build_anchors_from_calibration(
        "Gamma-1S", "NaI", 63.0, [(500.0, 30.0), (1000.0, 50.0)],
    )
    res = interpolate_peak_shape(img, 100.0)   # ниже min anchor
    assert res.was_extrapolated
    assert res.fwhm_keV > 0


def test_F300_batch_interpolate_returns_list():
    from gamma.peaks.peak_image_tabulated import (
        build_anchors_from_calibration,
    )
    from gamma.peaks.peak_image_logspline import batch_interpolate
    img = build_anchors_from_calibration(
        "Gamma-1S", "NaI", 63.0,
        [(122.0, 12.5), (662.0, 46.5), (1332.0, 78.0)],
    )
    shapes = batch_interpolate(img, [200.0, 500.0, 1000.0])
    assert len(shapes) == 3
    # FWHM монотонно растёт с E
    assert shapes[0].fwhm_keV < shapes[1].fwhm_keV < shapes[2].fwhm_keV


def test_F300_empty_anchors_raises():
    from gamma.peaks.peak_image_tabulated import TabulatedPeakImage
    from gamma.peaks.peak_image_logspline import interpolate_peak_shape
    img = TabulatedPeakImage(
        detector_id="x", detector_class="NaI",
        crystal_diameter_mm=63.0, anchors=[],
    )
    with pytest.raises(ValueError):
        interpolate_peak_shape(img, 500.0)


def test_F300_fwhm_at_E_convenience():
    from gamma.peaks.peak_image_tabulated import (
        build_anchors_from_calibration,
    )
    from gamma.peaks.peak_image_logspline import fwhm_at_E
    img = build_anchors_from_calibration(
        "Gamma-1S", "NaI", 63.0,
        [(122.0, 12.5), (662.0, 46.5)],
    )
    v = fwhm_at_E(img, 662.0)
    assert v == pytest.approx(46.5, abs=1e-6)


# ──────────────────────────────────────────────────────────────────
# F-301 — .cpt XML I/O
# ──────────────────────────────────────────────────────────────────

def test_F301_build_xml_contains_detector_and_anchors():
    from gamma.peaks.peak_image_tabulated import (
        build_anchors_from_calibration,
    )
    from gamma.io.cpt_io import build_cpt_xml
    img = build_anchors_from_calibration(
        "Gamma-1S", "NaI", 63.0,
        [(122.0, 12.5), (662.0, 46.5)],
    )
    xml = build_cpt_xml(img)
    assert '<?xml' in xml
    assert "Gamma-1S" in xml
    assert "NaI" in xml
    assert "662" in xml
    assert "<anchor" in xml
    assert "</peak_template>" in xml


def test_F301_round_trip_preserves_anchors(tmp_path):
    from gamma.peaks.peak_image_tabulated import (
        build_anchors_from_calibration,
    )
    from gamma.io.cpt_io import write_cpt_file, read_cpt_file
    img = build_anchors_from_calibration(
        "Gamma-1S", "NaI", 63.0,
        [(122.0, 12.5), (662.0, 46.5), (1332.0, 78.0)],
    )
    cpt_path = tmp_path / "test.cpt"
    write_cpt_file(img, cpt_path)

    loaded = read_cpt_file(cpt_path)
    assert loaded.detector_id == "Gamma-1S"
    assert loaded.detector_class == "NaI"
    assert loaded.crystal_diameter_mm == pytest.approx(63.0)
    assert len(loaded.anchors) == 3
    for orig, restored in zip(img.anchors, loaded.anchors):
        assert orig.E_keV == pytest.approx(restored.E_keV, abs=0.01)
        assert orig.fwhm_keV == pytest.approx(restored.fwhm_keV, abs=0.01)
        assert orig.tail_fraction == pytest.approx(
            restored.tail_fraction, abs=1e-5,
        )


def test_F301_parse_rejects_non_root_template():
    from gamma.io.cpt_io import parse_cpt_xml
    bad_xml = '<?xml version="1.0"?><wrong_root/>'
    with pytest.raises(ValueError):
        parse_cpt_xml(bad_xml)


def test_F301_parse_rejects_malformed_xml():
    from gamma.io.cpt_io import parse_cpt_xml
    bad_xml = '<peak_template><detector></peak_template>'  # неверная вложенность
    with pytest.raises(ValueError):
        parse_cpt_xml(bad_xml)


def test_F301_parse_warns_on_unknown_tag_non_strict(capsys):
    from gamma.io.cpt_io import parse_cpt_xml
    xml = '''<?xml version="1.0"?>
    <peak_template version="1.17.21">
        <detector id="x" class="NaI" diameter_mm="63"/>
        <anchors>
            <anchor E_keV="662" fwhm_keV="46.5"/>
        </anchors>
        <some_future_tag/>
    </peak_template>'''
    res = parse_cpt_xml(xml, strict=False)
    assert len(res.anchors) == 1
    captured = capsys.readouterr()
    assert "unknown" in captured.err.lower()


def test_F301_strict_mode_rejects_unknown_tag():
    from gamma.io.cpt_io import parse_cpt_xml
    xml = '''<?xml version="1.0"?>
    <peak_template version="1.17.21">
        <detector id="x" class="NaI" diameter_mm="63"/>
        <anchors><anchor E_keV="662" fwhm_keV="46.5"/></anchors>
        <bad_tag/>
    </peak_template>'''
    with pytest.raises(ValueError):
        parse_cpt_xml(xml, strict=True)


def test_F301_preserves_source_metadata_and_notes(tmp_path):
    from gamma.peaks.peak_image_tabulated import (
        TabulatedPeakImage, PeakShapeAnchor,
    )
    from gamma.io.cpt_io import write_cpt_file, read_cpt_file
    img = TabulatedPeakImage(
        detector_id="Gamma-1S", detector_class="NaI", crystal_diameter_mm=63.0,
        anchors=[PeakShapeAnchor(E_keV=662.0, fwhm_keV=46.5, tail_fraction=0.03)],
        source_metadata="Cs-137_2026-05-30",
        notes="Test calibration template",
    )
    p = tmp_path / "meta.cpt"
    write_cpt_file(img, p)
    loaded = read_cpt_file(p)
    assert loaded.source_metadata == "Cs-137_2026-05-30"
    assert loaded.notes == "Test calibration template"
