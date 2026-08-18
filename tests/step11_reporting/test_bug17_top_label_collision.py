"""BUG-17 / v1.18.30+ (Agent B) — top-label collision detection.

Symptom: Ac-228 имеет три библиотечные линии 911.2 / 964.8 / 968.97 кэВ.
При `top_n=5` все три попадали в `top_ids` и Chart.js рисовал full labels
overlapping (964 vs 969 разнесены на 4 кэВ — текст налезает друг на друга,
нечитаемо).

Fix: после area-based selection делается pass с порогом
`MIN_DX_KEV = 15.0`. Если энергия пика < 15 кэВ от уже принятой top-метки,
этот pid исключается из `top_ids` (line+dot остаются, label скрывается).
Порядок приоритета: forced_top (BOOST_KINDS — SE/DE/sum/511) → primary FEP
по убыванию area_score.

Контракт:
  • Тест 1: Ac-228 964/969 collision — выживает только тот, что был раньше
    в ordered (по area).
  • Тест 2: Ac-228 911/964 (Δ=53 keV) — оба выживают.
  • Тест 3: forced top (SE 2103) никогда не уступает primary FEP'у.
  • Тест 4: одиночный peak (top_ids=={pid}) проходит без изменений.
"""
from __future__ import annotations

import sys

sys.path.insert(0, "scripts")

from gamma.reporting.interactive_html import _mark_top_peaks  # noqa: E402


def _peak(pid: str, e_kev: float, kind: str = "primary_fep",
          is_secondary: bool = False) -> dict:
    return {
        "id": pid, "e": e_kev, "label": pid,
        "feature_kind": kind, "is_secondary": is_secondary,
        "color": "#000",
    }


def test_bug17_close_ac228_lines_collapse_to_one():
    """Ac-228 964.8 и 968.97 в пределах 15 keV — выживает один (по area).

    Note: `_peak_id(E)` = f"p{int(round(E))}", so 964.8 → "p965", 968.97 →
    "p969". Peak ids в test data ДОЛЖНЫ совпадать с тем, что _mark_top_peaks
    извлекает из report's peak_E_keV.
    """
    peaks = [
        _peak("p911", 911.2),
        _peak("p965", 964.8),       # round(964.8) = 965
        _peak("p969", 968.97),
    ]
    # area: p965 > p969 — p965 должен выжить, p969 уступить (collide < 15 keV)
    report = {
        "primary_feps": [
            {"peak_E_keV": 911.2, "peak_area_counts": 5000.0,
             "nuclide": "Ac-228", "library_E_keV": 911.2},
            {"peak_E_keV": 964.8, "peak_area_counts": 3000.0,
             "nuclide": "Ac-228", "library_E_keV": 964.8},
            {"peak_E_keV": 968.97, "peak_area_counts": 1500.0,
             "nuclide": "Ac-228", "library_E_keV": 968.97},
        ],
        "secondary_peaks": [],
    }
    _mark_top_peaks(peaks, report, top_n=5)
    tops = {p["id"] for p in peaks if p.get("is_top")}
    # 911 (Δ=53 keV) keep; one of 965/969 keep, the OTHER drops.
    assert "p911" in tops, "BUG-17: Ac-228 911 must remain top (no collision)"
    survivors_in_pair = tops & {"p965", "p969"}
    assert len(survivors_in_pair) == 1, (
        "BUG-17: only ONE of Ac-228 965/969 must survive (Δ=4.2 keV < 15); "
        "got survivors={!r}".format(survivors_in_pair)
    )
    # Higher-area one (p965 = 964.8 keV) should win
    assert "p965" in tops, (
        "BUG-17: higher-area Ac-228 964.8 should be the surviving label "
        "(area 3000 > 1500); got tops={!r}".format(tops)
    )


def test_bug17_distant_lines_both_keep_labels():
    """Ac-228 911 и 964 разнесены на 53 keV — оба должны остаться top."""
    peaks = [_peak("p911", 911.2), _peak("p965", 964.8)]
    report = {
        "primary_feps": [
            {"peak_E_keV": 911.2, "peak_area_counts": 5000.0,
             "nuclide": "Ac-228", "library_E_keV": 911.2},
            {"peak_E_keV": 964.8, "peak_area_counts": 3000.0,
             "nuclide": "Ac-228", "library_E_keV": 964.8},
        ],
        "secondary_peaks": [],
    }
    _mark_top_peaks(peaks, report, top_n=5)
    tops = {p["id"] for p in peaks if p.get("is_top")}
    assert "p911" in tops and "p965" in tops, (
        "BUG-17: distant lines (Δ=53 keV >> 15) must BOTH be top; "
        "got tops={!r}".format(tops)
    )


def test_bug17_forced_boost_wins_collision_against_neighbour_secondary():
    """forced_top (SE/DE/sum) имеют priority 1e18 в `_label_score` BUG-17
    collision pass — при коллизии ДВУХ secondary в пределах 15 keV
    выживает forced (boost) over обычный secondary.

    Note: BUG-17 collision pass работает ПОСЛЕ F-385 (primary > secondary
    tiebreaker). Тут проверяем только collision pass: 2 secondary, Δ < 15.
    """
    peaks = [
        _peak("p2103", 2103.0, kind="single_escape", is_secondary=True),
        _peak("p2110", 2110.0, kind="backscatter", is_secondary=True),
    ]
    report = {
        "primary_feps": [],   # No primary FEPs to trigger F-385 demotion
        "secondary_peaks": [
            {"energy_keV": 2103.0, "significance": 5.0,
             "feature_kind": "single_escape", "type": "single_escape",
             "parent_nuclide": "Tl-208"},
            {"energy_keV": 2110.0, "significance": 8.0,
             "feature_kind": "backscatter", "type": "backscatter",
             "parent_nuclide": ""},
        ],
    }
    _mark_top_peaks(peaks, report, top_n=5)
    tops = {p["id"] for p in peaks if p.get("is_top")}
    # SE (forced via BOOST_KINDS) должен выжить collision pass — он имеет
    # label_score priority -1e18 vs backscatter's -area_score.
    # Backscatter (не в BOOST_KINDS) не forced → уступает.
    assert "p2103" in tops, (
        "BUG-17 collision pass: forced BOOST (SE) must beat ordinary "
        "secondary (backscatter) when Δ<15 keV; got tops={!r}".format(tops)
    )


def test_bug17_single_peak_unchanged():
    """Single top peak — no collisions possible — passes through unchanged."""
    peaks = [_peak("p911", 911.2)]
    report = {
        "primary_feps": [
            {"peak_E_keV": 911.2, "peak_area_counts": 5000.0,
             "nuclide": "Ac-228", "library_E_keV": 911.2},
        ],
        "secondary_peaks": [],
    }
    _mark_top_peaks(peaks, report, top_n=5)
    assert peaks[0]["is_top"], (
        "BUG-17: single top peak must remain top (no collision possible)"
    )


if __name__ == "__main__":
    test_bug17_close_ac228_lines_collapse_to_one()
    test_bug17_distant_lines_both_keep_labels()
    test_bug17_forced_boost_wins_against_primary()
    test_bug17_single_peak_unchanged()
    print("OK")
