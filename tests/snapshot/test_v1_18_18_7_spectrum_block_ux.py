"""F-334 / v1.18.18.7 — Spectrum block UX overhaul regression suite.

Подзадачи:
  • F-334.1 — bg-режим показывает реальные пики фона (peak-search + match)
  • F-334.2 — Y-units cps/counts toggle, default cps, disabled в overlay
  • F-334.3 — modal popup с zoom/pan/reset/details
  • F-334.4 — top-5 пиков labelled, rest tooltip-only
"""
from __future__ import annotations

import json
import os
import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, "scripts")

from gamma.reporting import analyze_and_report
from gamma.reporting.interactive_html import (
    _BG_LINES_DICT,
    _detect_bg_peaks,
    _mark_top_peaks,
    _match_bg_isotope,
)


ROOT = Path(__file__).resolve().parents[2]
KIT = ROOT / "detectors/Gamma-1S/reference_spectra/reference_kits"
CS_SAMPLE = KIT / "Marinelli_1L/Cs-137/sample_M_cs_легкий_2001-2005.spe"
CS_BG = KIT / "Marinelli_1L/Cs-137/background_bg_2016_marinelli_water_marinelli.spe"
K_SAMPLE = KIT / "Marinelli_1L/K-40/sample_M_k_легкий_2001-2005.spe"
K_BG = KIT / "Marinelli_1L/K-40/background_bg_2016_marinelli_water_marinelli.spe"


@pytest.fixture(scope="module")
def cs137_html(tmp_path_factory):
    """Generate Cs-137 Marinelli report with bg — used by multiple tests."""
    if not CS_SAMPLE.exists():
        pytest.skip(f"Kit sample missing: {CS_SAMPLE}")
    out_dir = tmp_path_factory.mktemp("f334_cs137")
    res = analyze_and_report(
        str(CS_SAMPLE),
        background_path=str(CS_BG) if CS_BG.exists() else None,
        sample_mass_kg=0.570,
        output_dir=str(out_dir),
        write_html=True, write_plots=False, write_markdown=False,
    )
    html_path = Path(res["html"])
    assert html_path.exists()
    return html_path.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def k40_html(tmp_path_factory):
    """Generate K-40 Marinelli report with bg — used for top-N variation."""
    if not K_SAMPLE.exists():
        pytest.skip(f"Kit sample missing: {K_SAMPLE}")
    out_dir = tmp_path_factory.mktemp("f334_k40")
    res = analyze_and_report(
        str(K_SAMPLE),
        background_path=str(K_BG) if K_BG.exists() else None,
        sample_mass_kg=0.860,
        output_dir=str(out_dir),
        write_html=True, write_plots=False, write_markdown=False,
    )
    return Path(res["html"]).read_text(encoding="utf-8")


# ─────────────────────────────────────────────────────────────────
# F-334.1 — Real background peak markers in bg-view
# ─────────────────────────────────────────────────────────────────

def test_f334_1_bg_lines_dict_has_norm_lines():
    """Static dict должен содержать ключевые NORM-линии."""
    energies = [E for E, _ in _BG_LINES_DICT]
    isotopes = [iso for _, iso in _BG_LINES_DICT]
    # K-40 1460.8 — обязательно
    assert any(abs(E - 1460.8) < 0.2 for E in energies), "K-40 1460.8 missing"
    # Tl-208 2614.5 (Th-цепочка верх)
    assert any(abs(E - 2614.5) < 0.2 for E in energies), "Tl-208 2614.5 missing"
    # Bi-214 609 (U-цепочка ключевая)
    assert any(abs(E - 609.3) < 0.2 for E in energies), "Bi-214 609.3 missing"
    # Аннигиляция 511
    assert any(abs(E - 511.0) < 0.2 for E in energies), "annihilation 511 missing"


def test_f334_1_match_isotope_exact_K40():
    """K-40 1460.8 → matches itself."""
    m = _match_bg_isotope(1460.8)
    assert m is not None
    iso, E_lib = m
    assert "K-40" in iso
    assert abs(E_lib - 1460.8) < 0.5


def test_f334_1_match_isotope_fwhm_window_tolerates_gain_drift():
    """1443 кэВ (gain-drift ~17 кэВ) должен matchится в K-40 при FWHM-window."""
    m = _match_bg_isotope(1443.7)  # default FWHM-based tolerance
    assert m is not None, "FWHM-based tolerance must accept 17 keV gain drift at 1460"
    assert "K-40" in m[0]


def test_f334_1_match_isotope_no_match_far_energy():
    """Энергия далеко от любой библиотечной линии → None.

    Выбираем 2000 keV — провал между Bi-214 1764 (~50 keV FWHM-window)
    и Bi-214 2204 (~55 keV FWHM-window). 2000 находится на 236 от 1764
    и 204 от 2204 — ВНЕ обоих окон.
    """
    m = _match_bg_isotope(2000.0)
    assert m is None, f"Expected no match at 2000 keV; got {m}"


def test_f334_1_detect_bg_peaks_returns_top5_marked(cs137_html):
    """В сгенерированном payload должны быть bg_peaks (Cs-137 с фоном)."""
    # Extract CHART payload JSON
    m = re.search(r"const CHART=(\{.*?\});", cs137_html, re.DOTALL)
    assert m, "CHART payload missing in HTML"
    chart = json.loads(m.group(1))
    assert chart.get("has_background") is True
    bg = chart.get("bg_peaks") or []
    assert len(bg) >= 3, f"Expected ≥3 bg peaks; got {len(bg)}"
    # Top marker count.
    # F-370 / v1.18.24.5 — диагностические chain-anchors (Tl-208 2614,
    # K-40 1461, Bi-214 609, Ac-228 911/969, Pb-212 238, Cs-137 661)
    # форсятся в top независимо от intensity. Поэтому upper-bound
    # расширен с 5 до 12 (top-5 intensity + до 7 forced anchors).
    n_top = sum(1 for p in bg if p.get("is_top"))
    assert 1 <= n_top <= 12, f"Expected 1-12 top peaks; got {n_top}"
    # Каждый peak имеет обязательные поля
    for p in bg:
        assert "e" in p and "label" in p and "is_top" in p
        assert isinstance(p["e"], (int, float))


def test_f334_1_bg_peaks_empty_when_no_background(tmp_path):
    """Без фона bg_peaks=[]."""
    res = analyze_and_report(
        str(CS_SAMPLE),
        background_path=None,
        sample_mass_kg=0.570,
        output_dir=str(tmp_path),
        write_html=True, write_plots=False, write_markdown=False,
    )
    html = Path(res["html"]).read_text(encoding="utf-8")
    m = re.search(r"const CHART=(\{.*?\});", html, re.DOTALL)
    chart = json.loads(m.group(1))
    assert chart.get("has_background") is False
    assert chart.get("bg_peaks") == []


def test_f334_1_setview_bg_uses_bg_peaks(cs137_html):
    """JS setView('bg') должен подменять annotations на CHART.bg_peaks."""
    assert "CHART.bg_peaks" in cs137_html, "JS не ссылается на CHART.bg_peaks"
    # Должна быть условная логика choice peakList
    assert "view === 'bg'" in cs137_html


# ─────────────────────────────────────────────────────────────────
# F-334.2 — cps/counts Y-units toggle
# ─────────────────────────────────────────────────────────────────

def test_f334_2_cps_counts_buttons_present(cs137_html):
    """В HTML присутствуют обе кнопки cps + counts."""
    assert 'data-units="cps"' in cs137_html
    assert 'data-units="counts"' in cs137_html
    # CSS класс fp-units-btn
    assert "fp-units-btn" in cs137_html
    # Handler setYUnits
    assert "function setYUnits" in cs137_html or "setYUnits = function" in cs137_html


def test_f334_2_default_cps_in_payload(cs137_html):
    """Payload default_y_units='cps' (per user Q4)."""
    m = re.search(r"const CHART=(\{.*?\});", cs137_html, re.DOTALL)
    chart = json.loads(m.group(1))
    assert chart.get("default_y_units") == "cps"


def test_f334_2_t_sample_in_payload(cs137_html):
    """t_sample > 0 для cps conversion (sample/clean view)."""
    m = re.search(r"const CHART=(\{.*?\});", cs137_html, re.DOTALL)
    chart = json.loads(m.group(1))
    assert chart.get("t_sample", 0) > 0


def test_f334_2_overlay_disables_cps(cs137_html):
    """В overlay режиме cps button должен быть disabled (нельзя mix)."""
    # Должно быть условие disabled = view === 'overlay'
    assert "view === 'overlay'" in cs137_html
    assert "cpsBtn.disabled" in cs137_html or "disabled = (view === 'overlay')" in cs137_html


def test_f334_2_cps_label_in_y_axis(cs137_html):
    """В JS должен быть 'имп/с' label для cps mode."""
    assert "имп/с" in cs137_html


# ─────────────────────────────────────────────────────────────────
# F-334.3 — Modal popup with zoom/pan/reset/details
# ─────────────────────────────────────────────────────────────────

def test_f334_3_modal_html_present(cs137_html):
    """Modal markup присутствует."""
    assert 'id="fp-modal"' in cs137_html
    assert 'id="fp-sp-modal"' in cs137_html  # modal chart canvas
    assert 'id="fp-modal-details"' in cs137_html
    # 3 кнопки управления
    assert 'data-action="reset-zoom"' in cs137_html
    assert 'data-action="close"' in cs137_html
    assert 'data-action="toggle-details"' in cs137_html


def test_f334_3_expand_button_present(cs137_html):
    """Кнопка «⛶ Развернуть» рядом с chart controls."""
    assert 'id="fp-expand"' in cs137_html
    assert "Развернуть" in cs137_html


def test_f334_3_zoom_cdn_loaded(cs137_html):
    """chartjs-plugin-zoom CDN tag присутствует (graceful fallback if fails)."""
    assert "chartjs-plugin-zoom" in cs137_html
    assert "hammerjs" in cs137_html  # touch support dependency


def test_f334_3_rmb_pan_handler_present(cs137_html):
    """JS handler для right mouse button pan."""
    # contextmenu preventDefault + mousedown button=2 check
    assert "contextmenu" in cs137_html
    assert "button !== 2" in cs137_html or "e.button === 2" in cs137_html or "button === 2" in cs137_html
    assert "panDragging" in cs137_html


def test_f334_3_reset_zoom_function_present(cs137_html):
    """resetModalZoom function — fits Y axis to data range."""
    assert "function resetModalZoom" in cs137_html
    # Должна вызывать resetZoom + autofit
    assert "resetZoom" in cs137_html


def test_f334_3_esc_closes_modal(cs137_html):
    """ESC keydown handler закрывает modal."""
    assert "Escape" in cs137_html
    assert "closeModal" in cs137_html


def test_f334_3_details_panel_clonable(cs137_html):
    """Details panel клонируется внутрь modal."""
    assert "syncModalDetails" in cs137_html
    # Source div ↔ modal details mirror
    assert "fp-modal-details" in cs137_html
    assert "fp-detail" in cs137_html


# ─────────────────────────────────────────────────────────────────
# F-334.4 — Top-5 labels + tooltip-only markers
# ─────────────────────────────────────────────────────────────────

def test_f334_4_peaks_have_is_top_field(cs137_html):
    """Каждый sample peak имеет is_top boolean field."""
    m = re.search(r"const peaks=(\[.*?\]);", cs137_html, re.DOTALL)
    assert m, "peaks array missing"
    peaks = json.loads(m.group(1))
    assert len(peaks) >= 1
    for p in peaks:
        assert "is_top" in p, f"peak {p['id']} missing is_top"
        assert isinstance(p["is_top"], bool)


def test_f334_4_top_count_max_5(k40_html):
    """Top-N count ≤ 5 (K-40 demo может иметь >5 peaks с chain-completeness)."""
    m = re.search(r"const peaks=(\[.*?\]);", k40_html, re.DOTALL)
    assert m
    peaks = json.loads(m.group(1))
    n_top = sum(1 for p in peaks if p.get("is_top"))
    assert n_top <= 5, f"Expected ≤5 top peaks; got {n_top}"


def test_f334_4_mark_top_peaks_unit():
    """_mark_top_peaks unit-test: 7 peaks, 2 с area, остальные 0 → 2 top."""
    peaks = [
        {"id": "p100", "e": 100.0, "label": "A 100", "color": "#000"},
        {"id": "p200", "e": 200.0, "label": "B 200", "color": "#000"},
        {"id": "p300", "e": 300.0, "label": "C 300", "color": "#000"},
        {"id": "p400", "e": 400.0, "label": "D 400", "color": "#000"},
    ]
    report = {
        "primary_feps": [
            {"peak_E_keV": 100.0, "peak_area_counts": 5000.0},
            {"peak_E_keV": 200.0, "peak_area_counts": 100.0},
        ],
        "secondary_peaks": [],
    }
    _mark_top_peaks(peaks, report, top_n=5)
    # Только 2 пика имеют ненулевую area → 2 top
    n_top = sum(1 for p in peaks if p["is_top"])
    assert n_top == 2
    # Высшие area → top=True
    assert next(p for p in peaks if p["id"] == "p100")["is_top"] is True
    assert next(p for p in peaks if p["id"] == "p200")["is_top"] is True
    # Без area → False
    assert next(p for p in peaks if p["id"] == "p300")["is_top"] is False
    assert next(p for p in peaks if p["id"] == "p400")["is_top"] is False


def test_f334_4_tooltip_handler_present(cs137_html):
    """JS tooltip-on-hover handler для non-top пиков."""
    assert "fp-peak-tooltip" in cs137_html
    assert "showPeakTip" in cs137_html or "peakTip" in cs137_html
    # mousemove handler
    assert "mousemove" in cs137_html


def test_f334_4_buildAnnots_respects_is_top(cs137_html):
    """buildAnnots function реализует is_top differentiation."""
    assert "function buildAnnots" in cs137_html
    assert "is_top" in cs137_html
    assert "isTop" in cs137_html  # local variable in buildAnnots


# ─────────────────────────────────────────────────────────────────
# Cross-cutting: regression — no broken existing functionality
# ─────────────────────────────────────────────────────────────────

def test_f334_existing_chart_anchors_still_present(cs137_html):
    """F-114 form anchors не сломаны новыми правками."""
    must_have = [
        '<canvas id="fp-sp"',
        '<table class="fp-tbl">',
        "__initReport",
        "Chart.js не загрузился",
        'data-view="sample"',  # F-332 4-way toggle still there
        'data-view="bg"',
        'data-view="overlay"',
        'data-view="net"',
    ]
    for anchor in must_have:
        assert anchor in cs137_html, f"missing anchor: {anchor!r}"


def test_f334_version_bump():
    """SKILL_VERSION должен быть ≥ v1.18.18.7 (F-334 introduced bump).

    Note: SKILL_VERSION формат `vMAJOR.MINOR.PATCH[.HOTFIX]` — версия
    может уйти в v1.18.19+/v1.19+, не привязывайся к точному номеру.
    """
    from gamma.reporting.json_report import SKILL_VERSION
    import re
    m = re.match(r"^v(\d+)\.(\d+)\.(\d+)(?:\.(\d+))?$", SKILL_VERSION)
    assert m, f"Unexpected SKILL_VERSION format: {SKILL_VERSION}"
    major, minor, patch, hotfix = (int(g) if g else 0 for g in m.groups())
    # F-334 ≥ v1.18.18.7 (in (major, minor, patch, hotfix) order)
    assert (major, minor, patch, hotfix) >= (1, 18, 18, 7), \
        f"Expected ≥ v1.18.18.7; got {SKILL_VERSION}"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
