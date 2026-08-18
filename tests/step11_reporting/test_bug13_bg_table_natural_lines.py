"""BUG-13 / v1.18.30+ (Agent B) — bg-peak table must include natural-BG lines
that the plot already annotates.

Symptom: на view «Фон» plot подписывает K-40 1460.8 keV (`_detect_bg_peaks`
matches against `_BG_LINES_DICT`), но в таблице «Пики фона» этой строки
нет — `_build_rows(bg_report_view, ...)` берёт `background_primary_feps`
которые проходят chain-filtering и K-40 туда не попадает. Plot и таблица
расходятся → оператор недоумевает «где K-40?».

Fix: `_augment_bg_rows_with_natural_lines(rows_bg, analysis_result)` зовёт
ТУ ЖЕ функцию `_detect_bg_peaks(analysis_result)`, что и plot, и добавляет
строки для каждого matched isotope (K-40, Cs-137, Pb K-XRF, U-235, 511…),
если не покрыт через `primary_feps`.

Test strategy: monkey-patch `_detect_bg_peaks` to return a known list of
matched isotopes — это изолирует augmenter от тяжёлой peak-search infra.
"""
from __future__ import annotations

import sys

sys.path.insert(0, "scripts")

from gamma.reporting import interactive_html as ih  # noqa: E402
from gamma.reporting.interactive_html import (  # noqa: E402
    _augment_bg_rows_with_natural_lines,
    _NATURAL_BG_T_HALF_RU,
)


def _stub_detect(*, isos):
    """Build a stub `_detect_bg_peaks` return value (list-of-dicts) for
    the given (isotope, E_lib, E_meas, intensity) tuples."""
    return [
        {
            "id": "bg{}".format(int(round(em))),
            "e": em,
            "label": "{} {:.1f}".format(iso, el),
            "isotope": iso,
            "E_lib": el,
            "intensity": inten,
            "color": "#888780",
            "is_top": True,
        }
        for (iso, el, em, inten) in isos
    ]


def test_bug13_natural_lines_added_when_plot_finds_them(monkeypatch):
    """When _detect_bg_peaks finds K-40 via library matching, augment must
    add a corresponding row even if `primary_feps` did not include K-40."""
    monkeypatch.setattr(
        ih, "_detect_bg_peaks",
        lambda *_a, **_k: _stub_detect(
            isos=[("K-40", 1460.8, 1461.2, 1200.0)]
        ),
    )
    rows_bg = []
    augmented = _augment_bg_rows_with_natural_lines(rows_bg, analysis_result=None)
    k40_rows = [r for r in augmented if r.get("iso") == "K-40"]
    assert k40_rows, (
        "BUG-13: K-40 must appear in bg-rows after augmentation; "
        "rows={!r}".format(augmented)
    )
    assert k40_rows[0]["t"] == _NATURAL_BG_T_HALF_RU["K-40"]
    assert k40_rows[0]["kls"] == "fp-nat"
    assert k40_rows[0]["section"] == "primary_detected"


def test_bug13_does_not_duplicate_existing_rows(monkeypatch):
    """If K-40 row already in rows_bg (from primary_feps), augment must NOT
    add a duplicate."""
    monkeypatch.setattr(
        ih, "_detect_bg_peaks",
        lambda *_a, **_k: _stub_detect(
            isos=[("K-40", 1460.8, 1461.2, 1200.0)]
        ),
    )
    rows_bg = [
        {
            "peak": "p1461", "iso": "K-40", "kls": "fp-nat",
            "line": "1460.8 / 1460.8", "t": "1.25×10⁹ лет",
            "a": "—", "cmt": "from primary_feps",
            "section": "primary_detected",
        }
    ]
    augmented = _augment_bg_rows_with_natural_lines(rows_bg, analysis_result=None)
    k40_rows = [r for r in augmented if r.get("iso") == "K-40"]
    assert len(k40_rows) == 1, (
        "BUG-13: K-40 must NOT be duplicated; found {} K-40 rows; rows={!r}"
        .format(len(k40_rows), augmented)
    )


def test_bug13_skips_chain_nuclides(monkeypatch):
    """Pb-212/Tl-208 (Th-chain) уже идут через primary_feps path —
    augmenter не должен дублировать их."""
    monkeypatch.setattr(
        ih, "_detect_bg_peaks",
        lambda *_a, **_k: _stub_detect(
            isos=[
                ("Pb-212", 238.6, 238.5, 5000.0),
                ("Tl-208", 2614.5, 2614.0, 800.0),
                ("K-40",   1460.8, 1461.2, 1200.0),
            ]
        ),
    )
    augmented = _augment_bg_rows_with_natural_lines([], analysis_result=None)
    pb212 = [r for r in augmented if r.get("iso") == "Pb-212"]
    tl208 = [r for r in augmented if r.get("iso") == "Tl-208"]
    k40   = [r for r in augmented if r.get("iso") == "K-40"]
    assert not pb212 and not tl208, (
        "BUG-13: chain nuclides (Pb-212/Tl-208) must be skipped by augment — "
        "they belong to primary_feps path. Got: {!r}".format(augmented)
    )
    # But K-40 (natural, non-chain) MUST be added.
    assert k40, "BUG-13: non-chain K-40 must still be added by augment"


def test_bug13_empty_when_detect_returns_nothing(monkeypatch):
    """Defensive: when `_detect_bg_peaks` returns [] (no bg_grid / no peaks),
    augment returns rows_bg unchanged."""
    monkeypatch.setattr(ih, "_detect_bg_peaks", lambda *_a, **_k: [])
    original = [{"peak": "x", "iso": "Cs-137", "kls": "fp-nat",
                 "line": "—", "t": "—", "a": "—", "cmt": "",
                 "section": "primary_detected"}]
    result = _augment_bg_rows_with_natural_lines(original, analysis_result=None)
    assert result == original, (
        "BUG-13: when detect returns nothing, augment is a no-op."
    )


def test_bug13_adds_pb_xrf_and_cs137(monkeypatch):
    """Both Pb K-XRF (Pb-shield fluorescence) и Cs-137 (anthropogenic) могут
    появляться в фоне — оба должны попасть в таблицу."""
    monkeypatch.setattr(
        ih, "_detect_bg_peaks",
        lambda *_a, **_k: _stub_detect(
            isos=[
                ("Pb K-XRF", 40.0, 40.5, 600.0),
                ("Cs-137",   661.7, 661.0, 200.0),
                ("U-235",    185.7, 185.3, 80.0),
            ]
        ),
    )
    augmented = _augment_bg_rows_with_natural_lines([], analysis_result=None)
    isos = {r.get("iso") for r in augmented}
    assert "Pb K-XRF" in isos, (
        "BUG-13: Pb K-XRF (defensive shield fluorescence) must be added"
    )
    assert "Cs-137" in isos, (
        "BUG-13: Cs-137 anthropogenic contamination must be added"
    )
    assert "U-235" in isos, (
        "BUG-13: U-235 (Pb-shield impurity) must be added"
    )
    # Pb K-XRF and U-235 → fp-phys class (not natural-decay isotope)
    pb_row = next(r for r in augmented if r["iso"] == "Pb K-XRF")
    assert pb_row["kls"] == "fp-phys", "Pb K-XRF should be classed as fp-phys"


if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-v"])
