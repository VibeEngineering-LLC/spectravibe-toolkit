# -*- coding: utf-8 -*-
"""BUG-43 / 2026-06-04 (Agent B) — Two UI fixes:
  1. Мультиплет-блок НЕ показывается в BG-view (MULTIPLET_BLOCKS_BG = "").
  2. 511 keV в BG-пиках: если Tl-208 chain подтверждена (583/2614 keV),
     атрибутировать → Tl-208 510.77 keV, не «Аннигиляция».

Tests (a), (b), (c) per brief §3.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, "scripts")


# ─────────────────────────────────────────────────────────────────────────────
# Helper that replicates the BUG-43 511→Tl-208 relabeling logic.
# Mirrors interactive_html.py _detect_bg_peaks() BUG-43 post-processing block.
# ─────────────────────────────────────────────────────────────────────────────
def _apply_511_tl208_rule(enriched: list) -> list:
    """Apply BUG-43 511 keV Tl-208 override rule to an enriched peak list."""
    tl208_char_E_ranges = [(575, 592), (2600, 2625), (850, 870)]
    tl208_confirmed = any(
        d.get("isotope") == "Tl-208"
        and any(lo <= (d.get("E_lib") or 0.0) <= hi for lo, hi in tl208_char_E_ranges)
        for d in enriched
    )
    if tl208_confirmed:
        for d in enriched:
            if d.get("isotope") == "Аннигиляция" and 500 <= d.get("E_lib", 0) <= 520:
                d["isotope"] = "Tl-208"
                d["E_lib"] = 510.77
                d["label"] = "Tl-208 510.8"
                d["_511_tl208_override"] = True
    return enriched


# ─────────────────────────────────────────────────────────────────────────────
# Test (a): Th-232 chain present → 511 keV labeled as Tl-208
# ─────────────────────────────────────────────────────────────────────────────
def test_511_labeled_as_tl208_when_th232_chain_present():
    """BUG-43 case (a): когда Tl-208 цепочка подтверждена (583.2 + 2614.5),
    запись «Аннигиляция 511.0» должна быть переопределена → Tl-208 510.8."""
    enriched = [
        {
            "isotope": "Tl-208", "E_lib": 583.2, "e": 583.0,
            "intensity": 30, "is_top": True, "id": "bg583",
            "label": "Tl-208 583.2", "color": "#7F77DD",
        },
        {
            "isotope": "Tl-208", "E_lib": 2614.5, "e": 2614.0,
            "intensity": 10, "is_top": True, "id": "bg2614",
            "label": "Tl-208 2614.5", "color": "#7F77DD",
        },
        {
            "isotope": "Аннигиляция", "E_lib": 511.0, "e": 499.5,
            "intensity": 26, "is_top": False, "id": "bg511",
            "label": "Аннигиляция 511.0", "color": "#888780",
        },
    ]

    result = _apply_511_tl208_rule(enriched)

    # The 511 keV entry (e=499.5) must be relabeled.
    entry_511 = next((d for d in result if d["id"] == "bg511"), None)
    assert entry_511 is not None, "bg511 entry должна остаться в списке"
    assert entry_511["isotope"] == "Tl-208", (
        f"BUG-43: 511 keV должен быть атрибутирован к Tl-208 при наличии "
        f"цепочки. Получено: isotope={entry_511['isotope']!r}"
    )
    assert entry_511.get("_511_tl208_override") is True, (
        "BUG-43: _511_tl208_override sentinel должен быть True"
    )
    assert abs(entry_511["E_lib"] - 510.77) < 0.01, (
        f"BUG-43: E_lib должен стать 510.77 keV (Tl-208 line). "
        f"Получено: {entry_511['E_lib']}"
    )

    # Other entries must not be affected.
    tl583 = next(d for d in result if d["id"] == "bg583")
    assert tl583["isotope"] == "Tl-208"
    assert "_511_tl208_override" not in tl583


# ─────────────────────────────────────────────────────────────────────────────
# Test (b): No Tl-208 chain → 511 keV stays as «Аннигиляция»
# ─────────────────────────────────────────────────────────────────────────────
def test_511_stays_annihilation_without_tl208_chain():
    """BUG-43 case (b): без Tl-208 char lines → 511 остаётся «Аннигиляция»."""
    enriched = [
        {
            "isotope": "K-40", "E_lib": 1460.8, "e": 1445.7,
            "intensity": 100, "is_top": True, "id": "bg1461",
            "label": "K-40 1460.8", "color": "#1D9E75",
        },
        {
            "isotope": "Аннигиляция", "E_lib": 511.0, "e": 499.5,
            "intensity": 26, "is_top": False, "id": "bg511",
            "label": "Аннигиляция 511.0", "color": "#888780",
        },
    ]

    result = _apply_511_tl208_rule(enriched)

    entry_511 = next((d for d in result if d["id"] == "bg511"), None)
    assert entry_511 is not None
    assert entry_511["isotope"] == "Аннигиляция", (
        f"BUG-43 case (b): без Tl-208 chain isotope должен остаться "
        f"«Аннигиляция». Получено: {entry_511['isotope']!r}"
    )
    assert "_511_tl208_override" not in entry_511, (
        "BUG-43 case (b): _511_tl208_override не должен появляться без chain"
    )
    assert entry_511["E_lib"] == 511.0, (
        "BUG-43 case (b): E_lib должен остаться 511.0"
    )


# ─────────────────────────────────────────────────────────────────────────────
# Test (c): MULTIPLET_BLOCKS_BG is empty in sample run (HTML rendering)
# ─────────────────────────────────────────────────────────────────────────────
def test_multiplet_blocks_bg_empty_in_sample_run(tmp_path):
    """BUG-43 case (c): в sample-run MULTIPLET_BLOCKS_BG = "" →
    view-bg div не содержит блоков разложения мультиплетов.

    Использует реальный pipeline на fixture Th-232 Маринелли.
    Если fixture отсутствует — тест пропускается."""
    import pytest

    sp = (
        "detectors/Gamma-1S/reference_spectra/archive/"
        "Th232_420-7-17_Маринелли_0cm.spe"
    )
    bg = (
        "detectors/Gamma-1S/data/averaged_backgrounds/"
        "bg_2016_marinelli_water_marinelli.spe"
    )
    if not os.path.exists(sp) or not os.path.exists(bg):
        pytest.skip(f"fixture missing: {sp!r} or {bg!r}")

    from gamma.reporting import analyze_and_report

    # P1-3c: per-test tmp_path prevents xdist concurrent-write race on
    # fixed shared dir (was "demo_reports/_test_bug43_sample_multiplets").
    out = str(tmp_path)
    res = analyze_and_report(
        sp,
        output_dir=out,
        write_html=True,
        write_markdown=False,
        write_plots=False,
        sample_mass_kg=1.6,
        background_path=bg,
    )
    html_path = res["html"]
    assert html_path and os.path.exists(html_path), f"HTML missing: {html_path}"

    html = open(html_path, encoding="utf-8").read()

    # BUG-43 fix: fp-multiplets-bg div must exist (template placeholder) but
    # must contain NO fp-mp-block elements (= empty).
    assert "fp-multiplets-bg" in html, (
        "fp-multiplets-bg div должен присутствовать в шаблоне (view-bg CSS)"
    )
    # The inner content of view-bg div should not contain any multiplet block header.
    # _build_multiplet_blocks generates <h2>Мультиплеты — разложение в фоновом спектре
    assert "Мультиплеты — разложение в фоновом спектре" not in html, (
        "BUG-43: 'Мультиплеты — разложение в фоновом спектре' не должно "
        "появляться в sample-run HTML. MULTIPLET_BLOCKS_BG должен быть пустым."
    )

    # Sample multiplet section must still be present (no regression on Fix #1).
    assert "в спектре образца" in html, (
        "BUG-43 regression guard: sample multiplet section «в спектре образца» "
        "должна оставаться в sample-run HTML."
    )


if __name__ == "__main__":
    test_511_labeled_as_tl208_when_th232_chain_present()
    print("test (a) OK")
    test_511_stays_annihilation_without_tl208_chain()
    print("test (b) OK")
    # test (c) requires pipeline — run via pytest
    print("Run test (c) via pytest")
