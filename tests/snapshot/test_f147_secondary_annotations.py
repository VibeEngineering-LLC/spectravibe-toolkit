"""F-147 / v1.18.22.0 — Secondary annotations (Compton/SE/DE/...) toggle.

Validates:
  * _build_peaks помечает secondary feature_kinds полем is_secondary=True
  * primary_fep / chain_completeness получают is_secondary=False
  * Шаблон HTML включает fp-secondary-toggle div, initSecondaryToggle JS,
    let secondaryVisible, fp-sec-btn class
  * Существующий F-332 4-way toggle НЕ сломан (setView не получает
    secondary-кнопки благодаря data-view фильтру)
  * При генерации Th-232 demo ≥1 peak помечен is_secondary=True
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, "scripts")

ROOT = Path(__file__).resolve().parents[2]
KIT = ROOT / "detectors/Gamma-1S/reference_spectra/reference_kits"
TH_SAMPLE = KIT / "Marinelli_1L/Th-232/sample_M_th_легкий_2001-2005.spe"
TH_BG = KIT / "Marinelli_1L/Th-232/background_bg_2016_marinelli_water_marinelli.spe"


# ──────────────────────────────────────────────────────────────────
# Unit tests — _build_peaks classification of is_secondary
# ──────────────────────────────────────────────────────────────────

def test_f147_build_peaks_marks_compton_edge_secondary():
    """_build_peaks → secondary_peaks с feature_kind=compton_edge получают
    is_secondary=True."""
    from gamma.reporting.interactive_html import _build_peaks

    report = {
        "primary_feps": [
            {"nuclide": "Cs-137", "peak_E_keV": 661.7, "library_E_keV": 661.7,
             "peak_area_counts": 50000.0},
        ],
        "secondary_peaks": [
            {"energy_keV": 477.0, "feature_kind": "compton_edge",
             "parent_nuclide": "Cs-137", "parent_line_keV": 661.7,
             "type": "computed_feature"},
        ],
    }
    peaks = _build_peaks(report)
    assert len(peaks) >= 2
    by_kind = {p["feature_kind"]: p for p in peaks}
    assert by_kind["primary_fep"]["is_secondary"] is False
    assert by_kind["compton_edge"]["is_secondary"] is True


def test_f147_build_peaks_marks_se_de_secondary():
    """SE/DE/sum_peak/annihilation/backscatter — is_secondary=True."""
    from gamma.reporting.interactive_html import _build_peaks

    report = {
        "primary_feps": [],
        "secondary_peaks": [
            {"energy_keV": 2103.0, "feature_kind": "single_escape",
             "parent_nuclide": "Tl-208", "parent_line_keV": 2614.5},
            {"energy_keV": 1593.0, "feature_kind": "double_escape",
             "parent_nuclide": "Tl-208", "parent_line_keV": 2614.5},
            {"energy_keV": 511.0, "feature_kind": "annihilation_511"},
            {"energy_keV": 200.0, "feature_kind": "backscatter",
             "parent_nuclide": "Tl-208"},
            {"energy_keV": 1330.0, "feature_kind": "sum_peak",
             "parent_nuclide": "Co-60"},
        ],
    }
    peaks = _build_peaks(report)
    assert len(peaks) == 5
    for p in peaks:
        assert p["is_secondary"] is True, (
            f"{p['feature_kind']} should be is_secondary=True"
        )


def test_f147_build_peaks_cluster_NOT_secondary():
    """composite_cluster — это группа реальных γ-линий (Pb K-РИ + ВК),
    НЕ артефакт детектора → НЕ должен попадать в secondary toggle."""
    from gamma.reporting.interactive_html import _build_peaks

    report = {
        "primary_feps": [],
        "secondary_peaks": [
            {"energy_keV": 80.0, "feature_kind": "composite_cluster",
             "parent_nuclide": "Th-232"},
        ],
    }
    peaks = _build_peaks(report)
    assert len(peaks) == 1
    assert peaks[0]["feature_kind"] == "composite_cluster"
    assert peaks[0]["is_secondary"] is False


# ──────────────────────────────────────────────────────────────────
# Template tests — HTML / JS wiring
# ──────────────────────────────────────────────────────────────────

def test_f147_template_contains_toggle_markup():
    """Шаблон содержит fp-secondary-toggle div + fp-sec-btn class + 2 кнопки."""
    template_path = (
        ROOT / "scripts/gamma/reporting/templates/interactive_v1_17_2.html"
    )
    text = template_path.read_text(encoding="utf-8")
    assert 'id="fp-secondary-toggle"' in text
    assert 'class="fp-sec-btn fp-view-btn active"' in text
    assert 'data-secondary="on"' in text
    assert 'data-secondary="off"' in text
    assert "Вторичные процессы" in text


def test_f147_template_contains_init_function():
    """JS initSecondaryToggle IIFE + secondaryVisible let-binding."""
    template_path = (
        ROOT / "scripts/gamma/reporting/templates/interactive_v1_17_2.html"
    )
    text = template_path.read_text(encoding="utf-8")
    assert "initSecondaryToggle" in text
    assert "let secondaryVisible = true" in text
    assert "fp-sec-btn" in text


def test_f147_template_setView_handler_filters_data_view():
    """F-147 защита: setView() click-handler должен иметь guard
    `if (!btn.dataset.view) return;` чтобы fp-sec-btn (без data-view) не
    попадали туда. Селектор `.fp-view-btn` встречается несколько раз в
    шаблоне (active-state toggle, etc.); фильтр должен быть в handler,
    привязанном к setView(btn.dataset.view)."""
    template_path = (
        ROOT / "scripts/gamma/reporting/templates/interactive_v1_17_2.html"
    )
    text = template_path.read_text(encoding="utf-8")
    # Найти блок где привязан click → setView; рядом должен быть guard.
    idx = text.find("setView(btn.dataset.view)")
    assert idx >= 0, "не нашёл setView(btn.dataset.view) handler"
    # В пределах ±300 символов от него должен быть guard
    region = text[max(0, idx - 300):idx + 300]
    assert "if (!btn.dataset.view) return" in region, (
        "setView click-handler не имеет guard — fp-sec-btn попадёт в setView"
    )


# ──────────────────────────────────────────────────────────────────
# Integration — рендер реального Th-232 demo
# ──────────────────────────────────────────────────────────────────

def test_f147_html_includes_secondary_flag_in_payload(tmp_path):
    """Generated HTML содержит ≥1 is_secondary=true в DATA_PEAKS payload
    при анализе Th-232 (есть SE 2103 + DE 1593 + composite_cluster + ...).
    Composite_cluster — is_secondary=false; SE/DE — true."""
    if not TH_SAMPLE.exists():
        pytest.skip("Th-232 kit sample missing")
    from gamma.reporting import analyze_and_report

    res = analyze_and_report(
        str(TH_SAMPLE),
        background_path=str(TH_BG),
        sample_mass_kg=0.5,
        output_dir=str(tmp_path / "f147"),
        write_html=True, write_plots=False, write_markdown=False,
        write_json=True, write_technical_pdf=False,
    )
    html_path = Path(res["html"])
    assert html_path.exists()
    html = html_path.read_text(encoding="utf-8")

    # Payload содержит флаги
    assert '"is_secondary"' in html
    n_true = len(re.findall(r'"is_secondary"\s*:\s*true', html))
    n_false = len(re.findall(r'"is_secondary"\s*:\s*false', html))
    assert n_true >= 1, "Th-232 demo должен иметь ≥1 secondary peak (SE/DE)"
    assert n_false >= 1, "primary FEPs должны быть is_secondary=false"

    # Toggle UI вкомпилирован
    assert 'id="fp-secondary-toggle"' in html
    assert "Вторичные процессы" in html


# ──────────────────────────────────────────────────────────────────
# Version-bump assertion
# ──────────────────────────────────────────────────────────────────

def test_f147_version_bump():
    """Skill version >= (1, 18, 22, 0) — F-147 landed."""
    from gamma.reporting.json_report import SKILL_VERSION

    m = re.match(r"v(\d+)\.(\d+)\.(\d+)(?:\.(\d+))?", SKILL_VERSION)
    assert m, f"Unparseable SKILL_VERSION: {SKILL_VERSION}"
    parts = tuple(int(p or 0) for p in m.groups())
    assert parts >= (1, 18, 22, 0), (
        f"SKILL_VERSION {SKILL_VERSION} below F-147 baseline v1.18.22.0"
    )
