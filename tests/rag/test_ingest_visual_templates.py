# -*- coding: utf-8 -*-
"""Tests for scripts/rag/ingest_visual_templates.py  (F-070 Wave 1).

Eight required cases (brief §Required unit tests):
  1. Pointlike-5cm Am-241 ingest (real fixture .spe)
  2. Marinelli Cs-137 ingest (mocked)
  3. Block-list — heavy matrix M_cs_тяж_2001-2005.spe
  4. Block-list — background Фон_*.spe
  5. Block-list — stability run under Временная нестабильность/
  6. Поверка-2016 routing — Маринелли goes to _drift_study, NOT canonical
  7. Петри/Дента routing — canonical _raw_ingest/petri_60ml/
  8. Provenance integrity — absolute_source_path matches input verbatim
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

# conftest.py in tests/ adds scripts/ to sys.path automatically.
from rag.ingest_visual_templates import (
    classify_geometry,
    extract_nuclide,
    extract_year,
    ingest_one,
    is_blocklisted,
    make_template_id,
    _output_path_for_record,
)

# ─── Fixtures ─────────────────────────────────────────────────────────────────

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"
REAL_AM241_SPE = FIXTURES_DIR / "Am-241_045_02_2019_Точечная-5см_5cm.spe"

CORPUS_ROOT_MOCK = Path("/mock/corpus/archive")


def _make_mock_spectrum(
    n_channels: int = 1024,
    live_time: float = 1800.0,
    real_time: float = 1802.0,
) -> MagicMock:
    """Return a MagicMock that quacks like a Spectrum object."""
    spec = MagicMock()
    counts = np.zeros(n_channels, dtype=np.int64)
    counts[100:150] = 500  # non-zero values for realism
    spec.counts = counts
    spec.live_time = live_time
    spec.real_time = real_time
    spec.energy_cal = (-8.5, 3.04, 0.0)
    spec.stored_fwhm_calibration = None
    spec.comments = "Test spectrum"
    spec.n_channels = n_channels
    return spec


# ─── Test 1: Pointlike-5cm Am-241 ingest (real .spe fixture) ─────────────────

@pytest.mark.skipif(
    not REAL_AM241_SPE.exists(),
    reason="Am-241 fixture .spe not present in tests/rag/fixtures/",
)
def test_pointlike_5cm_am241_ingest_real_fixture():
    """Test 1 — real Am-241 pointlike-5cm .spe: geometry, nuclide, channels."""
    rec = ingest_one(REAL_AM241_SPE, FIXTURES_DIR)

    assert rec is not None, "ingest_one must not return None for a valid .spe"
    assert rec["geometry_class"] == "pointlike_5cm"
    assert rec["nuclide"] == "Am-241"
    assert isinstance(rec["raw_channels"], list)
    assert len(rec["raw_channels"]) > 0, "raw_channels must be non-empty"
    assert rec["n_channels"] == len(rec["raw_channels"])
    assert rec["schema_version"] == "0.1-raw"


# ─── Test 2: Marinelli Cs-137 ingest (mocked) ────────────────────────────────

def test_marinelli_cs137_ingest():
    """Test 2 — Cs137 Маринелли path emits geometry marinelli_0cm."""
    path = CORPUS_ROOT_MOCK / "Cs137_420-7-14_Маринелли_0cm.spe"
    mock_spec = _make_mock_spectrum()

    with patch("rag.ingest_visual_templates.read_lsrm_spe", return_value=mock_spec):
        rec = ingest_one(path, CORPUS_ROOT_MOCK)

    assert rec is not None
    assert rec["geometry_class"] == "marinelli_0cm"
    assert rec["nuclide"] == "Cs-137"


# ─── Test 3: Block-list — heavy matrix ───────────────────────────────────────

def test_blocklist_heavy_matrix():
    """Test 3 — M_cs_тяж_2001-2005.spe is block-listed and skipped."""
    path = CORPUS_ROOT_MOCK / "M_cs_тяж_2001-2005.spe"

    assert is_blocklisted(path) is True, "is_blocklisted must return True"

    # ingest_one must return None without calling read_lsrm_spe
    with patch("rag.ingest_visual_templates.read_lsrm_spe") as mock_reader:
        rec = ingest_one(path, CORPUS_ROOT_MOCK)
        mock_reader.assert_not_called()

    assert rec is None, "ingest_one must return None for block-listed file"


# ─── Test 4: Block-list — background ─────────────────────────────────────────

def test_blocklist_background_fon():
    """Test 4 — Фон_*.spe is block-listed and skipped."""
    path = CORPUS_ROOT_MOCK / "Фон_закр_кр_вода_01.spe"

    assert is_blocklisted(path) is True

    with patch("rag.ingest_visual_templates.read_lsrm_spe") as mock_reader:
        rec = ingest_one(path, CORPUS_ROOT_MOCK)
        mock_reader.assert_not_called()

    assert rec is None


def test_blocklist_background_bg():
    """Test 4b — bg_*.spe is block-listed."""
    path = CORPUS_ROOT_MOCK / "bg_background_run.spe"
    assert is_blocklisted(path) is True


# ─── Test 5: Block-list — stability run ──────────────────────────────────────

def test_blocklist_stability_run():
    """Test 5 — file under Временная нестабильность/ is block-listed."""
    path = Path("/corpus/Поверка-2016/Временная нестабильность/Cs137_run1.spe")

    assert is_blocklisted(path) is True

    with patch("rag.ingest_visual_templates.read_lsrm_spe") as mock_reader:
        rec = ingest_one(path, Path("/corpus"))
        mock_reader.assert_not_called()

    assert rec is None


# ─── Test 6: Поверка-2016 Маринелли → drift_study, NOT canonical ─────────────

def test_poverka2016_marinelli_routes_to_drift_study():
    """Test 6 — Поверка-2016/Маринелли path is marked drift study."""
    path = Path("/corpus/Поверка-2016/Маринелли/Cs137_420-7-14_Маринелли_0cm.spe")
    mock_spec = _make_mock_spectrum()

    geom, is_drift = classify_geometry(path)
    assert is_drift is True, "Поверка-2016/Маринелли must be drift_study=True"
    assert geom == "marinelli_0cm"

    with patch("rag.ingest_visual_templates.read_lsrm_spe", return_value=mock_spec):
        rec = ingest_one(path, Path("/corpus"))

    assert rec is not None
    assert rec["is_drift_study"] is True

    # Output must NOT go to canonical _raw_ingest/marinelli_0cm/
    out_path = _output_path_for_record(rec, Path("/output"))
    assert "_drift_study" in str(out_path), (
        f"Drift-study record must route to _drift_study/, got: {out_path}"
    )
    assert "_raw_ingest_poverka2016" in str(out_path)
    assert "marinelli_0cm" not in str(
        out_path.parent.parent
    ) or "_drift_study" in str(out_path), (
        "Must not appear in canonical _raw_ingest/marinelli_0cm/"
    )


# ─── Test 7: Поверка-2016 Петри → canonical _raw_ingest/petri_60ml/ ──────────

def test_petri_routes_to_canonical():
    """Test 7 — Поверка-2016/Чашка Петри 60мл/ is petri_60ml, NOT drift study."""
    path = Path("/corpus/Поверка-2016/Чашка Петри 60мл/Th-232 #420_SRC-03 Петри-60_Петри-60.spe")

    geom, is_drift = classify_geometry(path)
    assert geom == "petri_60ml", f"Expected petri_60ml, got {geom!r}"
    assert is_drift is False, "Петри 60мл must NOT be drift study"

    mock_spec = _make_mock_spectrum()
    # Patch nuclide extraction too — filename is atypical; test routing only.
    with patch("rag.ingest_visual_templates.read_lsrm_spe", return_value=mock_spec), \
         patch("rag.ingest_visual_templates.extract_nuclide", return_value="Th-232"):
        rec = ingest_one(path, Path("/corpus"))

    assert rec is not None
    assert rec["geometry_class"] == "petri_60ml"
    assert rec["is_drift_study"] is False

    out_path = _output_path_for_record(rec, Path("/output"))
    assert "_raw_ingest" in str(out_path)
    assert "_drift_study" not in str(out_path)
    assert "petri_60ml" in str(out_path)


# ─── Test 8: Provenance integrity (anti-hallucination) ────────────────────────

def test_provenance_absolute_path_matches_verbatim():
    """Test 8 — absolute_source_path must match str(path.resolve()) verbatim."""
    path = CORPUS_ROOT_MOCK / "Am-241_045_02_2019_Точечная-5см_5cm.spe"
    mock_spec = _make_mock_spectrum()

    with patch("rag.ingest_visual_templates.read_lsrm_spe", return_value=mock_spec):
        rec = ingest_one(path, CORPUS_ROOT_MOCK)

    assert rec is not None
    expected = str(path.resolve())
    assert rec["absolute_source_path"] == expected, (
        f"absolute_source_path mismatch:\n"
        f"  expected: {expected!r}\n"
        f"  got:      {rec['absolute_source_path']!r}"
    )


# ─── Additional unit tests for helper functions ───────────────────────────────

def test_extract_nuclide_aliases():
    """Nuclide alias normalisation (Cs137 → Cs-137, K40 → K-40)."""
    assert extract_nuclide("Cs137_420-7-14_Маринелли_0cm.spe") == "Cs-137"
    assert extract_nuclide("K40_420-7-20_Маринелли_0cm.spe") == "K-40"
    assert extract_nuclide("Ra226_420-7-18_Маринелли_0cm.spe") == "Ra-226"
    assert extract_nuclide("Th232_420-7-17_Маринелли.spe") == "Th-232"


def test_extract_nuclide_canonical_form():
    """Canonical nuclide forms are returned unchanged."""
    assert extract_nuclide("Am-241_045_02_2019_Точечная-5см_5cm.spe") == "Am-241"
    assert extract_nuclide("Eu-152__04_21_Точечная-5см_5cm.spe") == "Eu-152"
    assert extract_nuclide("Y-88__260_2023_Точечная-5см_5cm.spe") == "Y-88"


def test_extract_nuclide_unknown_returns_none():
    """Unrecognised first token returns None."""
    assert extract_nuclide("РИСН_something.spe") is None
    assert extract_nuclide("Unknown_X_Y.spe") is None


def test_extract_year():
    """Year extraction from filename patterns."""
    assert extract_year("Am-241_045_02_2019_Точечная-5см_5cm.spe") == "2019"
    assert extract_year("Cs-137__163_2017.spe") == "2017"
    assert extract_year("Th-228__264_2023_Точечная-5см_5cm.spe") == "2023"
    assert extract_year("Cs137_420-7-14_Маринелли_0cm.spe") is None
    # Certificate number (SRC-05) must NOT be treated as a year.
    assert extract_year("Ba-133__SRC-05_Точечная-5см_5cm.spe") is None


def test_make_template_id():
    """Template ID construction."""
    assert make_template_id("Cs-137", "marinelli_0cm", "2017") == "VT-CS137-MARINELLI0CM-2017"
    assert make_template_id("Am-241", "pointlike_5cm", "2019") == "VT-AM241-POINT5CM-2019"
    assert make_template_id("K-40", "marinelli_0cm", None) == "VT-K40-MARINELLI0CM-NONE"
    assert make_template_id("Th-232", "petri_60ml", "2016") == "VT-TH232-PETRI60ML-2016"


def test_blocklist_fon_lowercase():
    """фон_*.spe (lowercase) is also block-listed."""
    path = CORPUS_ROOT_MOCK / "фон_пустая_защита.spe"
    assert is_blocklisted(path) is True


def test_classify_geometry_denta():
    """Поверка-2016/Дента-100мл → denta_100ml, NOT drift study."""
    path = Path("/corpus/Поверка-2016/Дента-100мл/Th-232_test.spe")
    geom, is_drift = classify_geometry(path)
    assert geom == "denta_100ml"
    assert is_drift is False
