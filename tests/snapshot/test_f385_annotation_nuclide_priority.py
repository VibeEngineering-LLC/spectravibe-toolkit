# -*- coding: utf-8 -*-
"""F-385 / v1.18.26 — annotation tiebreaker (nuclide-FEP > secondary).

Контракт `_mark_top_peaks`: при коллизии в пределах ±1·FWHM по энергии
(NaI rough FWHM ≈ 6% E + 5 кэВ) secondary peak (SE/DE/sum/annihilation/
backscatter, `is_secondary=True`) уступает primary FEP'у, даже если он
был forced через BOOST_KINDS. Edge case — только secondary в slot →
secondary остаётся top.

Кейс из практики: общий спектр Th-232. DE 1593 (Tl-208 → 2614−1022)
forced as BOOST_KIND. Ac-228 1588 кэВ — primary_fep. Без F-385 label DE
перекрывал бы label Ac-228 на холсте; после F-385 — Ac-228 получает
top, DE остаётся в data но без top label.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))


# ──────────────────────────────────────────────────────────────────
# Core: DE 1593 vs Ac-228 1588 collision
# ──────────────────────────────────────────────────────────────────

def test_f385_de_loses_to_primary_fep_in_collision_window():
    """DE 1593 (BOOST forced) vs Ac-228 1588 (primary_fep) — primary wins.

    Энергии разнесены на 5 кэВ; rough NaI FWHM @ 1590 кэВ ≈ 100 кэВ,
    окно гарантированно покрывает обе линии → DE теряет is_top.
    """
    from gamma.reporting.interactive_html import _build_peaks, _mark_top_peaks

    report = {
        "primary_feps": [
            {
                "nuclide": "Ac-228",
                "peak_E_keV": 1588.0,
                "library_E_keV": 1588.2,
                "peak_area_counts": 12000.0,
            },
        ],
        "secondary_peaks": [
            {
                "energy_keV": 1593.0,
                "feature_kind": "double_escape",
                "parent_nuclide": "Tl-208",
                "parent_line_keV": 2614.5,
                "significance": 3.5,
            },
        ],
    }
    peaks = _build_peaks(report)
    _mark_top_peaks(peaks, report, top_n=5)

    by_kind = {p["feature_kind"]: p for p in peaks}
    assert "primary_fep" in by_kind, "Ac-228 1588 primary_fep отсутствует"
    assert "double_escape" in by_kind, "DE 1593 secondary отсутствует"

    # F-385 contract: primary wins, secondary is in data but not top
    assert by_kind["primary_fep"]["is_top"] is True, (
        "Ac-228 1588 primary_fep должен получить is_top=True"
    )
    assert by_kind["double_escape"]["is_top"] is False, (
        "DE 1593 secondary должен потерять is_top при коллизии с primary"
    )
    assert by_kind["double_escape"]["is_secondary"] is True, (
        "DE остаётся помеченным as is_secondary — это не меняется F-385"
    )


# ──────────────────────────────────────────────────────────────────
# Edge: только secondary в окне → остаётся top (нет жертвы)
# ──────────────────────────────────────────────────────────────────

def test_f385_isolated_secondary_keeps_top():
    """SE 2103 без соседнего primary_fep в окне — остаётся top."""
    from gamma.reporting.interactive_html import _build_peaks, _mark_top_peaks

    report = {
        "primary_feps": [
            # Дальний primary, > 200 кэВ от SE 2103 → не коллидирует
            {
                "nuclide": "Tl-208",
                "peak_E_keV": 2614.5,
                "library_E_keV": 2614.5,
                "peak_area_counts": 50000.0,
            },
        ],
        "secondary_peaks": [
            {
                "energy_keV": 2103.0,
                "feature_kind": "single_escape",
                "parent_nuclide": "Tl-208",
                "parent_line_keV": 2614.5,
                "significance": 4.0,
            },
        ],
    }
    peaks = _build_peaks(report)
    _mark_top_peaks(peaks, report, top_n=5)

    by_kind = {p["feature_kind"]: p for p in peaks}
    # SE 2103 — отдельный slot (FWHM @ 2103 ≈ 130 кэВ; Tl-208 2614 — 511 кэВ
    # дальше → нет коллизии)
    assert by_kind["single_escape"]["is_top"] is True, (
        "Изолированный SE 2103 должен остаться top (нет primary в slot)"
    )
    assert by_kind["primary_fep"]["is_top"] is True, (
        "Tl-208 2614 primary_fep — top по area"
    )


# ──────────────────────────────────────────────────────────────────
# Edge: secondary без primary вообще → top
# ──────────────────────────────────────────────────────────────────

def test_f385_secondary_alone_keeps_top():
    """Без single primary_fep в peaks list — F-385 не должен ничего ломать.

    SE/DE остаются forced top through BOOST_KINDS.
    """
    from gamma.reporting.interactive_html import _build_peaks, _mark_top_peaks

    report = {
        "primary_feps": [],
        "secondary_peaks": [
            {
                "energy_keV": 1593.0,
                "feature_kind": "double_escape",
                "parent_nuclide": "Tl-208",
                "parent_line_keV": 2614.5,
                "significance": 3.5,
            },
        ],
    }
    peaks = _build_peaks(report)
    _mark_top_peaks(peaks, report, top_n=5)
    assert len(peaks) == 1
    assert peaks[0]["feature_kind"] == "double_escape"
    assert peaks[0]["is_top"] is True, (
        "Без primary_fep в peaks — DE остаётся top (forced BOOST)"
    )


# ──────────────────────────────────────────────────────────────────
# Edge: non-BOOST secondary (compton_edge) vs primary тоже работает
# ──────────────────────────────────────────────────────────────────

def test_f385_compton_edge_yields_to_primary():
    """compton_edge не в BOOST_KINDS, но is_secondary=True. Если он попадает
    в top by significance score AND в окне primary — должен тоже уступить."""
    from gamma.reporting.interactive_html import _build_peaks, _mark_top_peaks

    report = {
        "primary_feps": [
            {
                "nuclide": "Cs-137",
                "peak_E_keV": 661.7,
                "library_E_keV": 661.7,
                "peak_area_counts": 100000.0,
            },
            # Близкий primary к compton_edge
            {
                "nuclide": "Bi-214",
                "peak_E_keV": 480.0,
                "library_E_keV": 480.4,
                "peak_area_counts": 8000.0,
            },
        ],
        "secondary_peaks": [
            {
                "energy_keV": 477.0,
                "feature_kind": "compton_edge",
                "parent_nuclide": "Cs-137",
                "parent_line_keV": 661.7,
                "significance": 5.0,
                "type": "computed_feature",
            },
        ],
    }
    peaks = _build_peaks(report)
    _mark_top_peaks(peaks, report, top_n=5)

    by_e = {round(p["e"]): p for p in peaks}
    # Bi-214 480 — primary в окне ±FWHM от Compton 477 (FWHM ≈ 34 кэВ)
    assert by_e[480]["is_top"] is True, "Bi-214 480 primary должен победить"
    # Compton 477 теряет top (collide с Bi-214 480)
    assert by_e[477]["is_top"] is False, (
        "Compton_edge 477 должен уступить primary Bi-214 480 в slot"
    )


# ──────────────────────────────────────────────────────────────────
# Regression: при отсутствии коллизий top-ranking не меняется
# ──────────────────────────────────────────────────────────────────

def test_f385_no_collisions_no_change():
    """Если все peaks разнесены по энергии вне FWHM окон — F-385 ничего
    не должен трогать, существующее ранкирование работает как раньше."""
    from gamma.reporting.interactive_html import _build_peaks, _mark_top_peaks

    report = {
        "primary_feps": [
            {"nuclide": "K-40", "peak_E_keV": 1460.8, "library_E_keV": 1460.8,
             "peak_area_counts": 30000.0},
            {"nuclide": "Cs-137", "peak_E_keV": 661.7, "library_E_keV": 661.7,
             "peak_area_counts": 200000.0},
        ],
        "secondary_peaks": [
            {"energy_keV": 511.0, "feature_kind": "annihilation_511",
             "significance": 5.0},
        ],
    }
    peaks = _build_peaks(report)
    _mark_top_peaks(peaks, report, top_n=5)

    by_e = {round(p["e"]): p for p in peaks}
    # 511 vs Cs-137 661.7 — разница 150 кэВ; FWHM @ 511 ≈ 35 кэВ → нет коллизии
    assert by_e[511]["is_top"] is True, "annihilation_511 forced остаётся top"
    assert by_e[662]["is_top"] is True, "Cs-137 661.7 primary top по area"
    assert by_e[1461]["is_top"] is True, "K-40 1460.8 primary top по area"


# ──────────────────────────────────────────────────────────────────
# Skill version guard
# ──────────────────────────────────────────────────────────────────

def test_f385_module_imports_clean():
    """Smoke: модуль импортируется после F-385 правки без ошибок."""
    from gamma.reporting.interactive_html import _mark_top_peaks  # noqa: F401
    assert _mark_top_peaks is not None
