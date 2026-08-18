"""G5 / v1.31.1 — stored-vs-rebuilt E-cal source label.

Unit-level test for `_build_calibration` (json_report) и `_build_subtitle`
(interactive_html): убедиться, что:

  1. Поле `calibration.energy_cal.source` берётся из `spec.energy_cal_source`,
     а не хардкодится «stored».
  2. Поле `reused` = True только когда source == "stored".
  3. Поле `source_label` — человекочитаемая RU-строка для известных кодов
     и `перестроена ({code})` fallback для неизвестных.
  4. `_build_subtitle` добавляет `· калибровка: {label}` chunk, когда
     передан `energy_cal` с непустым label, и не добавляет ничего, когда
     label == "источник не указан".
"""
from __future__ import annotations

import pathlib
import sys
from types import SimpleNamespace

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from gamma.reporting.json_report import _build_calibration
from gamma.reporting.interactive_html import _build_subtitle


def _make_result(source: str):
    spec = SimpleNamespace(
        energy_cal=(0.0, 3.0),
        energy_cal_source=source,
    )
    return SimpleNamespace(
        spec=spec,
        fwhm_model=(1.0, 0.01, 0.0),
        fwhm_model_source="bootstrap",
        fwhm_at_661=7.0,
        seven_line_check=None,
    )


def test_stored_source_marks_reused():
    cal = _build_calibration(_make_result("stored"))["energy_cal"]
    assert cal["source"] == "stored"
    assert cal["reused"] is True
    assert cal["source_label"] == "использована сохранённая E-cal"


def test_f145_self_cal_marks_not_reused():
    cal = _build_calibration(
        _make_result("F-145_multiplet_self_calibration")
    )["energy_cal"]
    assert cal["source"] == "F-145_multiplet_self_calibration"
    assert cal["reused"] is False
    assert cal["source_label"] == "перестроена (F-145 multiplet self-cal)"


def test_bootstrap_source_marks_not_reused():
    cal = _build_calibration(_make_result("bootstrap"))["energy_cal"]
    assert cal["source"] == "bootstrap"
    assert cal["reused"] is False
    assert cal["source_label"] == "перестроена (bootstrap)"


def test_unknown_source_falls_back_to_generic_rebuilt():
    cal = _build_calibration(_make_result("g7_experimental"))["energy_cal"]
    assert cal["source"] == "g7_experimental"
    assert cal["reused"] is False
    assert cal["source_label"] == "перестроена (g7_experimental)"


def test_empty_source_marks_unknown_and_not_reused():
    cal = _build_calibration(_make_result(""))["energy_cal"]
    assert cal["source"] == "unknown"
    assert cal["reused"] is False
    assert cal["source_label"] == "источник не указан"


def _minimal_header():
    return {
        "detector_canonical": "Gamma-1S",
        "n_channels": 1024,
        "live_time_s": 3600.0,
        "background_status": "",
        "background_subtracted": True,
        "sample_filename": "",
        "background_filename": "",
    }


def test_subtitle_appends_cal_chunk_for_stored():
    subtitle = _build_subtitle(
        _minimal_header(),
        diag={},
        energy_cal={
            "source": "stored",
            "reused": True,
            "source_label": "использована сохранённая E-cal",
        },
    )
    assert "· калибровка: использована сохранённая E-cal" in subtitle


def test_subtitle_skips_cal_chunk_when_source_unknown():
    subtitle = _build_subtitle(
        _minimal_header(),
        diag={},
        energy_cal={
            "source": "unknown",
            "reused": False,
            "source_label": "источник не указан",
        },
    )
    assert "калибровка:" not in subtitle


def test_subtitle_without_energy_cal_kwarg_is_backward_compatible():
    subtitle = _build_subtitle(_minimal_header(), diag={})
    assert "калибровка:" not in subtitle
    assert "Gamma-1S" in subtitle