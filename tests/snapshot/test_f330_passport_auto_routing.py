# -*- coding: utf-8 -*-
"""F-330 / v1.18.18.4 — Auto-routing passport activities из LSRM .spe
COMMENT в `passport_activity_Bq` без явной передачи пользователем.

Контракт:
1. spec.extras["lsrm_passport"] populated reader-ом при загрузке .spe.
2. wrapper._auto_passport_from_spec() конвертирует Бк/кг → Бк по массе
   образца, применяет decay correction до даты измерения, возвращает
   passport_dict + meta-dict с provenance.
3. build_report принимает passport_meta и render-ит provenance preamble
   в HTML/MD passport comparison section.
4. Когда user passes passport_activity_Bq явно, auto-routing skipped
   (explicit > auto).
"""
from __future__ import annotations
import sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))


# ─── spec.extras populated by reader ────────────────────────────────

def test_F330_lsrm_reader_populates_passport_extras():
    """После read_lsrm_spe spec.extras["lsrm_passport"] непуст для
    калибровочного эталона с passport-данными в COMMENT."""
    from gamma.io.lsrm_spe import read_lsrm_spe
    fixture = ROOT / "evals" / "fixtures" / "M_cs_легкий_2001-2005.spe"
    spec = read_lsrm_spe(str(fixture))
    assert "lsrm_passport" in spec.extras
    p = spec.extras["lsrm_passport"]
    assert isinstance(p, list)
    assert len(p) == 1
    e = p[0]
    assert e["nuclide"] == "Cs-137"
    assert e["value"] == 1890.0
    assert e["unit"] == "Бк/кг"
    assert e["uncertainty_pct"] == 5.0
    assert e["reference_date"] == "1997-05-30"
    assert e["is_specific_activity"] is True


def test_F330_lsrm_reader_M_k_cyrillic_K_normalized():
    """К-40 в COMMENT → K-40 в spec.extras."""
    from gamma.io.lsrm_spe import read_lsrm_spe
    fixture = ROOT / "evals" / "fixtures" / "M_k_легкий_2001-2005.spe"
    spec = read_lsrm_spe(str(fixture))
    p = spec.extras.get("lsrm_passport") or []
    assert len(p) == 1
    assert p[0]["nuclide"] == "K-40"
    assert p[0]["reference_date"] is None  # M_k не указывает дату


# ─── wrapper auto-routing ───────────────────────────────────────────

def _build_stub_result(extras: dict, mass_kg: float | None, meas_date):
    """Lightweight stub чтобы не запускать full pipeline."""
    from dataclasses import dataclass
    from datetime import datetime

    @dataclass
    class _Spec:
        extras: dict
        start_datetime: object | None = None

    if mass_kg is not None:
        extras = {**extras, "lsrm_sample_mass_kg": mass_kg}

    sd = datetime.combine(meas_date, datetime.min.time()) if meas_date else None

    @dataclass
    class _R:
        spec: _Spec

    return _R(spec=_Spec(extras=extras, start_datetime=sd))


def test_F330_auto_passport_specific_activity_decayed():
    """Specific activity Бк/кг × mass × decay correction → Бк."""
    from gamma.reporting.wrapper import _auto_passport_from_spec
    from datetime import date

    extras = {
        "lsrm_passport": [{
            "nuclide": "Cs-137",
            "value": 1890.0,
            "unit": "Бк/кг",
            "uncertainty_pct": 5.0,
            "reference_date": "1997-05-30",
            "is_specific_activity": True,
        }],
    }
    result = _build_stub_result(extras, mass_kg=0.570,
                                meas_date=date(1999, 8, 4))
    passport, meta = _auto_passport_from_spec(result, explicit=None)
    assert passport is not None
    assert "Cs-137" in passport
    # 1890 Бк/кг × 0.570 кг = 1077.3 Bq; decay 1997→1999 → ≈1024.5 Bq
    assert passport["Cs-137"] == pytest.approx(1024.5, rel=2e-3)
    assert meta["source"] == "spe_comment"
    assert meta["mass_kg"] == 0.570
    assert meta["meas_date"] == "1999-08-04"
    assert meta["decay_corrected"]["Cs-137"] is True
    assert meta["ref_dates"]["Cs-137"] == "1997-05-30"


def test_F330_auto_passport_no_ref_date_no_decay():
    """Без reference_date — decay-correction опущена, A0 × mass."""
    from gamma.reporting.wrapper import _auto_passport_from_spec
    from datetime import date

    extras = {
        "lsrm_passport": [{
            "nuclide": "K-40",
            "value": 2540.0,
            "unit": "Бк/кг",
            "uncertainty_pct": 10.0,
            "reference_date": None,
            "is_specific_activity": True,
        }],
    }
    result = _build_stub_result(extras, mass_kg=0.665,
                                meas_date=date(1999, 8, 4))
    passport, meta = _auto_passport_from_spec(result, explicit=None)
    # 2540 × 0.665 = 1689.1 Bq
    assert passport["K-40"] == pytest.approx(1689.1, rel=1e-3)
    assert meta["decay_corrected"]["K-40"] is False


def test_F330_auto_passport_no_mass_skips_specific():
    """Specific activity без mass_kg → entry пропущен с notes."""
    from gamma.reporting.wrapper import _auto_passport_from_spec
    from datetime import date

    extras = {
        "lsrm_passport": [{
            "nuclide": "Cs-137",
            "value": 1890.0,
            "unit": "Бк/кг",
            "uncertainty_pct": 5.0,
            "reference_date": "1997-05-30",
            "is_specific_activity": True,
        }],
    }
    result = _build_stub_result(extras, mass_kg=None,
                                meas_date=date(1999, 8, 4))
    passport, meta = _auto_passport_from_spec(result, explicit=None)
    # Conversion skipped — passport dict empty → None
    assert passport is None
    assert "масса образца не известна" in meta["notes"]


def test_F330_auto_passport_absolute_Bq_no_mass_needed():
    """Абсолютная Бк (не specific) — mass не нужна."""
    from gamma.reporting.wrapper import _auto_passport_from_spec
    from datetime import date

    extras = {
        "lsrm_passport": [{
            "nuclide": "Co-60",
            "value": 1050.0,
            "unit": "Бк",
            "uncertainty_pct": 5.0,
            "reference_date": "2020-01-01",
            "is_specific_activity": False,
        }],
    }
    result = _build_stub_result(extras, mass_kg=None,
                                meas_date=date(2025, 1, 1))
    passport, meta = _auto_passport_from_spec(result, explicit=None)
    assert passport is not None
    # Co-60 t½ = 5.27 y; Δt = 5y → exp(-ln(2)·5/5.27) = 0.5179
    # 1050 × 0.5179 ≈ 543.8 Bq
    assert passport["Co-60"] == pytest.approx(543.8, rel=1e-2)


def test_F330_explicit_passport_overrides_auto():
    """Когда user passes explicit, auto skipped, meta.source='explicit'."""
    from gamma.reporting.wrapper import _auto_passport_from_spec
    from datetime import date

    extras = {
        "lsrm_passport": [{
            "nuclide": "Cs-137",
            "value": 1890.0,
            "unit": "Бк/кг",
            "uncertainty_pct": 5.0,
            "reference_date": "1997-05-30",
            "is_specific_activity": True,
        }],
    }
    result = _build_stub_result(extras, mass_kg=0.570,
                                meas_date=date(1999, 8, 4))
    passport, meta = _auto_passport_from_spec(
        result, explicit={"Cs-137": 2000.0},
    )
    # Explicit case — function returns (None, meta-source=explicit) and
    # wrapper uses the explicit value directly.
    assert passport is None
    assert meta["source"] == "explicit"


def test_F330_no_passport_extras_returns_none():
    """Spec без lsrm_passport extras → no auto."""
    from gamma.reporting.wrapper import _auto_passport_from_spec
    from datetime import date

    result = _build_stub_result(extras={}, mass_kg=0.5,
                                meas_date=date(2025, 1, 1))
    passport, meta = _auto_passport_from_spec(result, explicit=None)
    assert passport is None
    assert meta["source"] == "none"


# ─── HTML render of provenance ───────────────────────────────────────

def test_F330_html_render_provenance_preamble():
    """_f326_append_passport_comparison HTML branch emits provenance."""
    from gamma.reporting.build import _f326_append_passport_comparison

    html = '<html><body><div class="page">main</div></body></html>'
    passport = {"Cs-137": 1024.5}
    activities = {"Cs-137": 1038.0}
    meta = {
        "source": "spe_comment",
        "mass_kg": 0.570,
        "meas_date": "1999-08-04",
        "ref_dates": {"Cs-137": "1997-05-30"},
        "decay_corrected": {"Cs-137": True},
        "notes": "",
    }
    out = _f326_append_passport_comparison(
        html, format="html", passport=passport,
        activities=activities, passport_meta=meta,
    )
    # F-337.6 / v1.18.19.1 — фраза «Источник паспорта: автоматически
    # извлечено из поля COMMENT файла .spe» удалена по запросу пользователя
    # как избыточная. Provenance теперь видна через mass/meas_date + decay
    # blocks (которые присутствуют тогда и только тогда, когда meta заполнена).
    assert "масса образца 0.570 кг" in out
    assert "дата измерения 1999-08-04" in out
    assert "Cs-137: ref=1997-05-30" in out
    assert "decay-correction применён" in out
    # Comparison table present
    assert "<table class=\"passport-tbl\">" in out
    assert "Cs-137" in out


def test_F330_md_render_provenance_preamble():
    from gamma.reporting.build import _f326_append_passport_comparison

    md = "# Report\n\nSome content.\n"
    passport = {"K-40": 1689.1}
    activities = {"K-40": 1700.0}
    meta = {
        "source": "spe_comment",
        "mass_kg": 0.665,
        "meas_date": "1999-08-04",
        "ref_dates": {"K-40": None},
        "decay_corrected": {"K-40": False},
        "notes": "",
    }
    out = _f326_append_passport_comparison(
        md, format="md", passport=passport,
        activities=activities, passport_meta=meta,
    )
    assert "Сравнение с паспортной удельной активностью" in out
    # F-337.6 / v1.18.19.1 — фраза «Источник паспорта: …» удалена;
    # provenance видна через параметры пересчёта (масса, дата, decay).
    assert "масса образца 0.665 кг" in out
    assert "ref-дата не указана" in out
    assert "| K-40 |" in out


def test_F330_explicit_source_label():
    """F-337.6: explicit source больше не имеет отдельной фразы —
    проверяем что блок рендерится без падения с meta.source='explicit'."""
    from gamma.reporting.build import _f326_append_passport_comparison

    html = '<html><body><div class="page">x</div></body></html>'
    passport = {"Cs-137": 1000.0}
    activities = {"Cs-137": 950.0}
    meta = {
        "source": "explicit",
        "nuclides": ["Cs-137"],
        "mass_kg": None,
        "meas_date": None,
        "ref_dates": {},
        "decay_corrected": {},
    }
    out = _f326_append_passport_comparison(
        html, format="html", passport=passport,
        activities=activities, passport_meta=meta,
    )
    # Comparison table должна быть, deferred message — отсутствовать
    assert 'class="passport-tbl"' in out
    assert "Сравнение не выполнено" not in out
    assert "Cs-137" in out


def test_F330_deferred_message_unchanged_when_no_data():
    """passport=None путь не ломается, выводится deferred message."""
    from gamma.reporting.build import _f326_append_passport_comparison

    html = '<html><body><div class="page">x</div></body></html>'
    out = _f326_append_passport_comparison(
        html, format="html", passport=None, activities={},
        passport_meta={"source": "none"},
    )
    assert "Сравнение не выполнено" in out
    assert "passport_activity_Bq" in out
