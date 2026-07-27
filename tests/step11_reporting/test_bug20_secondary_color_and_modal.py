"""BUG-20 / v1.18.30+ (Agent B) — secondary processes must use distinct
color + «Вторичные процессы» toggle must work in expand/fullscreen modal.

Symptom #1: SE/DE/sum/511/backscatter/compton_edge получали parent's chain
color (Th-232 SE/DE были синие, как primary FEPs Th-цепочки) → оператор не
отличает артефакт детектора от реальной γ-линии.

Symptom #2: «Вторичные процессы» toggle (показать/скрыть) был только в
normal view header. При открытии fullscreen modal (z-index:9999, fixed
inset:0) основной toggle скрывался за overlay → невозможно скрыть SE/DE в
expand-режиме без выхода из fullscreen.

Fix #1: `_COL_SECONDARY = "#E8884F"` (orange) + legend swatch «Вторичные
процессы». Применяется регardless of parent nuclide для всех feature_kinds
в F147_SECONDARY_PHYS_KINDS.

Fix #2: дублирующие кнопки в modal header (id="fp-modal-sec-grp",
.fp-modal-sec-btn) + общий handler `setSecondaryVisibility(on)` сохраняет
sync двусторонний (normal click → modal active state + vice versa).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, "scripts")

from gamma.reporting.interactive_html import (  # noqa: E402
    _build_peaks,
    _build_legend_items,
    _COL_SECONDARY,
    _COL_TH,
    _COL_NAT,
    _COL_PHYS,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = (REPO_ROOT / "scripts" / "gamma" / "reporting" / "templates"
            / "interactive_v1_17_2.html")


# ──── Fix #1 tests ─────────────────────────────────────────────────


def test_bug20_secondary_color_constant_exists():
    """`_COL_SECONDARY` should be defined and distinct from chain/nat/phys."""
    assert _COL_SECONDARY, "BUG-20: _COL_SECONDARY constant is required"
    assert _COL_SECONDARY != _COL_TH, (
        "BUG-20: secondary color must differ from chain color _COL_TH"
    )
    assert _COL_SECONDARY != _COL_NAT, (
        "BUG-20: secondary color must differ from natural color _COL_NAT"
    )
    assert _COL_SECONDARY != _COL_PHYS, (
        "BUG-20: secondary color must differ from generic phys color"
    )
    # Sanity: hex format
    assert re.match(r"^#[0-9A-Fa-f]{6}$", _COL_SECONDARY), (
        "BUG-20: _COL_SECONDARY should be 6-digit hex"
    )


def test_bug20_secondary_peaks_use_orange_regardless_of_parent():
    """SE/DE/sum/backscatter/compton/511 ВСЕГДА получают _COL_SECONDARY,
    even если parent_nuclide=Th-232 (раньше брал _COL_TH через _chain_color)."""
    report = {
        "primary_feps": [],
        "secondary_peaks": [
            {"energy_keV": 2103.0, "type": "single_escape",
             "feature_kind": "single_escape",
             "parent_nuclide": "Tl-208", "significance": 4.0},
            {"energy_keV": 1593.0, "type": "double_escape",
             "feature_kind": "double_escape",
             "parent_nuclide": "Tl-208", "significance": 4.0},
            {"energy_keV": 511.0, "type": "annihilation_511",
             "feature_kind": "annihilation_511",
             "parent_nuclide": "", "significance": 5.0},
            {"energy_keV": 200.0, "type": "backscatter",
             "feature_kind": "backscatter",
             "parent_nuclide": "Cs-137", "significance": 3.5},
        ],
        "diagnostics": {},
    }
    peaks = _build_peaks(report)
    sec_peaks = [p for p in peaks if p.get("is_secondary")]
    assert len(sec_peaks) == 4, (
        "expected 4 secondary peaks; got {}: {!r}".format(len(sec_peaks),
                                                          sec_peaks)
    )
    for p in sec_peaks:
        assert p["color"] == _COL_SECONDARY, (
            "BUG-20: secondary peak {} must use _COL_SECONDARY ({}), got {!r}"
            .format(p["label"], _COL_SECONDARY, p["color"])
        )


def test_bug20_primary_feps_keep_chain_color():
    """Primary FEPs (Tl-208 583, Pb-212 238) НЕ должны менять цвет —
    BUG-20 fix #1 не должен регрессировать primary path."""
    report = {
        "primary_feps": [
            {"peak_E_keV": 238.6, "peak_area_counts": 5000,
             "nuclide": "Pb-212", "library_E_keV": 238.6},
            {"peak_E_keV": 583.2, "peak_area_counts": 8000,
             "nuclide": "Tl-208", "library_E_keV": 583.2},
        ],
        "secondary_peaks": [],
        "diagnostics": {},
    }
    peaks = _build_peaks(report)
    for p in peaks:
        if not p.get("is_secondary"):
            assert p["color"] == _COL_TH, (
                "BUG-20 regression: primary FEP {} of Th chain must keep "
                "_COL_TH ({}), got {!r}".format(p["label"], _COL_TH, p["color"])
            )


def test_bug20_legend_has_secondary_swatch():
    """`_build_legend_items` должен включать swatch для «Вторичные процессы»."""
    legend = _build_legend_items()
    assert _COL_SECONDARY in legend, (
        "BUG-20: legend HTML must include _COL_SECONDARY swatch; got:\n{}"
        .format(legend)
    )
    assert "Вторичные процессы" in legend, (
        "BUG-20: legend must label the secondary-processes swatch in Russian"
    )


# ──── Fix #2 tests (modal toggle wiring) ───────────────────────────


def _read_template() -> str:
    assert TEMPLATE.exists(), f"missing: {TEMPLATE}"
    return TEMPLATE.read_text(encoding="utf-8")


def test_bug20_modal_has_secondary_toggle_buttons():
    """Modal header должен содержать дублирующие кнопки toggle секцию
    `data-modal-secondary="on"` и `="off"` (для fullscreen-режима)."""
    tmpl = _read_template()
    assert 'data-modal-secondary="on"' in tmpl, (
        "BUG-20 fix #2: modal must have toggle button data-modal-secondary='on'"
    )
    assert 'data-modal-secondary="off"' in tmpl, (
        "BUG-20 fix #2: modal must have toggle button data-modal-secondary='off'"
    )
    assert 'id="fp-modal-sec-grp"' in tmpl, (
        "BUG-20 fix #2: modal must have `fp-modal-sec-grp` label container"
    )


def test_bug20_modal_click_handler_wired():
    """Modal click delegate должен вызвать setSecondaryVisibility() при
    клике на дублирующие кнопки."""
    tmpl = _read_template()
    assert "setSecondaryVisibility" in tmpl, (
        "BUG-20 fix #2: shared `setSecondaryVisibility` must be defined"
    )
    # delegate handler dispatches modalSecondary dataset attr
    assert "modalSecondary" in tmpl, (
        "BUG-20 fix #2: modal click delegate must read e.target.dataset"
        ".modalSecondary"
    )


if __name__ == "__main__":
    test_bug20_secondary_color_constant_exists()
    test_bug20_secondary_peaks_use_orange_regardless_of_parent()
    test_bug20_primary_feps_keep_chain_color()
    test_bug20_legend_has_secondary_swatch()
    test_bug20_modal_has_secondary_toggle_buttons()
    test_bug20_modal_click_handler_wired()
    print("OK")
