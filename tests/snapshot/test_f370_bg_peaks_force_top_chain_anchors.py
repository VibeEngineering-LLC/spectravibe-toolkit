# -*- coding: utf-8 -*-
"""F-370 / v1.18.24.5 — Forced top labels для диагностических цепочко-анкеров
в bg-peaks (Chart.js view «Фон»).

Регрессия: пользователь сообщил «2614 не указан в фоновом спектре в обоих
версиях». В bg-peaks payload пик Tl-208 2614 кэВ ПРИСУТСТВУЕТ
(intensity≈6.7 native counts/s), но не попадает в top-5 по чистому
intensity-ранкингу (top-5: Pb K-X / Pb-212 238 / K-40 1461 / 28 / 145).
Без `is_top=True` annotation рисуется как маленький dot БЕЗ полной
подписи — пользователь видит «голую» область 2614 на графике фона.

Tl-208 2614 — top Th-232-анкер (99.75% I_γ, чистая ROI), наличие в фоне
— главное доказательство Th-232 цепочки в lab-shielding. Аналогично
K-40 1461 (одиночник), Bi-214 609 (Ra-226 анкер), Pb-212 238 (low-E Th).
"""
from __future__ import annotations
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))


def _run_th232_bg_detection():
    """Helper: возвращает bg_peaks от Th-232 demo."""
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    from gamma.reporting.interactive_html import _detect_bg_peaks

    sample = (
        REPO / "detectors" / "Gamma-1S" / "reference_spectra"
        / "reference_kits" / "Marinelli_1L" / "Th-232"
        / "Th232_420-7-17_Маринелли_0cm.spe"
    )
    bg = (
        REPO / "detectors" / "Gamma-1S" / "reference_spectra"
        / "reference_kits" / "Marinelli_1L" / "Th-232"
        / "Фон закр кр вода_13.spe"
    )
    if not sample.is_file() or not bg.is_file():
        pytest.skip("Th-232 kit files missing")
    res = analyze_lsrm_spe(
        str(sample), background_path=str(bg),
        sample_mass_kg=0.5, complete_workflow=True,
    )
    return _detect_bg_peaks(res, top_n=20)


def test_F370_bg_peaks_tl208_2614_is_top():
    """Tl-208 2614 кэВ должен быть is_top=True даже при слабой intensity."""
    bg = _run_th232_bg_detection()
    matches = [
        p for p in bg
        if (p.get("isotope") == "Tl-208"
            and 2600 <= (p.get("E_lib") or 0) <= 2625)
    ]
    assert matches, (
        f"Tl-208 2614 не найден в bg_peaks (всего {len(bg)} пиков): "
        f"{[(p['e'], p.get('isotope')) for p in bg]}"
    )
    assert any(p.get("is_top") for p in matches), (
        f"Tl-208 2614 не получил is_top=True — F-370 boost регрессия. "
        f"Matches: {matches}"
    )


def test_F370_bg_peaks_k40_1461_is_top():
    """K-40 1461 кэВ — одиночник, всегда top."""
    bg = _run_th232_bg_detection()
    matches = [
        p for p in bg
        if (p.get("isotope") == "K-40"
            and 1450 <= (p.get("E_lib") or 0) <= 1475)
    ]
    assert matches, "K-40 1461 не найден в bg_peaks"
    assert any(p.get("is_top") for p in matches), (
        "K-40 1461 должен быть is_top=True"
    )


def test_F370_bg_peaks_bi214_609_is_top():
    """Bi-214 609 — Ra-226 анкер, всегда top."""
    bg = _run_th232_bg_detection()
    matches = [
        p for p in bg
        if (p.get("isotope") == "Bi-214"
            and 605 <= (p.get("E_lib") or 0) <= 615)
    ]
    if not matches:
        pytest.skip("Bi-214 609 не виден в этом фоне")
    assert any(p.get("is_top") for p in matches), (
        "Bi-214 609 должен быть is_top=True (Ra-226 анкер)"
    )


def test_F370_th232_demo_html_has_tl208_2614_annotated():
    """End-to-end: demo HTML содержит Tl-208 2614 как top-annotated в
    bg_peaks JSON payload (а не только в sample-peaks)."""
    import re
    import json
    demo = (
        REPO / "demo_reports" / "v1_18_24_th232_full" / "sample"
        / "Th232_Маринелли_0cm_report.html"
    )
    if not demo.is_file():
        pytest.skip(f"demo HTML missing: {demo}")
    html = demo.read_text(encoding="utf-8")
    m = re.search(r"const CHART=(\{.*?\});", html, re.DOTALL)
    assert m, "CHART payload not injected"
    payload = json.loads(m.group(1))
    bg_peaks = payload.get("bg_peaks", [])
    tl208_2614 = [
        p for p in bg_peaks
        if (p.get("isotope") == "Tl-208"
            and 2600 <= (p.get("E_lib") or 0) <= 2625)
    ]
    assert tl208_2614, (
        f"Tl-208 2614 не в demo bg_peaks: "
        f"{[(p.get('e'), p.get('isotope')) for p in bg_peaks]}"
    )
    assert any(p.get("is_top") for p in tl208_2614), (
        "Tl-208 2614 не is_top в demo HTML — F-370 regression"
    )
