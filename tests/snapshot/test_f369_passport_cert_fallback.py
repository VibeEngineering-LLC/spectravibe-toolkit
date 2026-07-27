# -*- coding: utf-8 -*-
"""F-369 / v1.18.24.4 — Passport auto-extract из .src сертификатов LSRM.

Регрессия: пользователь сообщил «опять паспорт источника не читается для
обоих версий» (2026-06-01). Production-pipeline предусматривает
auto-extract паспорта из LSRM .spe COMMENT (F-330), но многие реальные
.spe файлы (включая Th232_420-7-17_Маринелли_0cm.spe) НЕ содержат
inline-паспортных entries — реальные сертификаты лежат в
`detectors/Gamma-1S/certificates/*.src` (LSRM Аспект INI-формат).

Fix: добавлен fallback `_passport_from_certificate()` — после проверки
`spec.extras["lsrm_passport"]` ищет sub-source по serial из filename
(например `420-7-17`) в сертификатах, конвертирует Бк/кг → Бк по mass
из .src, корректирует на распад от reference_datetime до даты измерения.
"""
from __future__ import annotations
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))


def test_F369_extract_serial_from_filename():
    """Helper извлечения serial из имени .spe."""
    from gamma.reporting.wrapper import _extract_serial_from_filename
    assert _extract_serial_from_filename("Th232_420-7-17_Маринелли_0cm.spe") == "420-7-17"
    assert _extract_serial_from_filename("Cs137_420-7-15_Дента-120мл_5cm.spe") == "420-7-15"
    assert _extract_serial_from_filename("Ra226_420-7-18_Петри-60мл_0cm.spe") == "420-7-18"
    # No serial — должно вернуть None
    assert _extract_serial_from_filename("M_cs_легкий_2001-2005.spe") is None
    assert _extract_serial_from_filename("bg_2016.spe") is None


def test_F369_filename_geometry_hint():
    """Извлечение canonical geometry из filename."""
    from gamma.reporting.wrapper import _filename_geometry_hint
    assert _filename_geometry_hint("Th232_420-7-17_Маринелли_0cm.spe") == "маринелли"
    assert _filename_geometry_hint("Cs137_420-7-14_Дента-120мл_0cm.spe") == "дента"
    assert _filename_geometry_hint("Ra226_420-7-19_Петри-60мл_0cm.spe") == "петри"
    # Without recognised hint
    assert _filename_geometry_hint("noname.spe") == ""


def test_F369_cert_passport_th232_marinelli():
    """End-to-end: Th-232 / Marinelli kit находит паспорт через
    fallback и возвращает Th-232 ≈ 3104 Бк (1940 Бк/кг × 1.6 кг)."""
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    from gamma.reporting.wrapper import _auto_passport_from_spec

    sample = (
        REPO / "detectors" / "Gamma-1S" / "reference_spectra"
        / "reference_kits" / "Marinelli_1L" / "Th-232"
        / "Th232_420-7-17_Маринелли_0cm.spe"
    )
    if not sample.is_file():
        pytest.skip(f"demo .spe missing: {sample}")

    result = analyze_lsrm_spe(str(sample), sample_mass_kg=0.5,
                              complete_workflow=True)
    passport, meta = _auto_passport_from_spec(result, explicit=None)

    assert passport is not None, (
        "F-369 fallback должен найти паспорт через .src certificate, "
        f"meta: {meta}"
    )
    assert "Th-232" in passport, f"Th-232 missing in passport: {passport}"
    # 1940 Бк/кг × 1.6 кг = 3104 Бк, после decay correction для Th-232
    # (T½ = 1.4·10¹⁰ years) — практически тот же 3104 Бк
    A_Bq = passport["Th-232"]
    assert 3000 < A_Bq < 3200, (
        f"Th-232 passport activity expected ≈ 3104 Бк, got {A_Bq}"
    )
    # Meta sanity.
    # BUG-49 (2026-06-04): после расширения F-330 parser на Поверка-2016
    # COMMENT format (line «Th-232 A=1940 Бк/кг dA=6% 17-09-2007» в .spe),
    # источник паспорта для этого файла теперь — «spe_comment», а не
    # «cert_src». Оба пути дают идентичную числовую активность
    # (3104 Бк) с одинаковой ref_date / decay-correction / mass — что и
    # является меaningful invariant теста. Принимаем оба source-tag'а,
    # чтобы тест оставался валидным как для .spe с inline-паспортом,
    # так и для .spe без него (где fallback на .src сертификат
    # сохраняется).
    assert meta["source"] in ("cert_src", "spe_comment"), (
        f"meta.source: {meta['source']}"
    )
    if meta["source"] == "cert_src":
        assert meta["cert_file"].endswith(".src"), (
            f"cert_file: {meta['cert_file']}"
        )
        assert "Th232_420-7-17" in (meta["cert_subsource"] or ""), (
            f"cert_subsource: {meta['cert_subsource']}"
        )
    assert meta["ref_dates"]["Th-232"] == "2007-09-17"
    assert meta["decay_corrected"]["Th-232"] is True
    assert meta["mass_kg"] == 1.6


def test_F369_passport_block_appears_in_th232_demo_html():
    """End-to-end: production demo HTML должен содержать ‹section
    passport-comparison› с реальными цифрами Th-232, а НЕ «Сравнение не
    выполнено — данные паспорта источника не переданы»."""
    demo = (
        REPO / "demo_reports" / "v1_18_24_th232_full" / "sample"
        / "Th232_Маринелли_0cm_report.html"
    )
    if not demo.is_file():
        pytest.skip(f"demo HTML missing: {demo}")
    html = demo.read_text(encoding="utf-8")
    # passport block присутствует
    assert "passport-comparison" in html, (
        "passport-comparison секции нет в demo HTML"
    )
    # Регрессионная фраза НЕ должна появляться, т.к. fallback нашёл паспорт
    assert "данные паспорта источника не переданы" not in html, (
        "регрессия F-369: HTML говорит «паспорт не передан» при "
        "наличии .src сертификата для 420-7-17"
    )
    # Должна быть строка с decay-correction и cert-derived activity
    assert "ref=2007-09-17" in html or "2007-09-17" in html, (
        "ref date паспорта (2007-09-17) не отображается в HTML"
    )
