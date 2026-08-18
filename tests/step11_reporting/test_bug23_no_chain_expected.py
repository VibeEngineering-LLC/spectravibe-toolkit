"""BUG-23 / v1.18.31+ (Agent B) — peak tables must NEVER contain a
'chain_expected' section.

Symptom: после BUG-19 sample peak table получила 4 секции, включая «Линии
из цепочки (по равновесию, F-111)» (`section == "chain_expected"`).
Эти строки имеют «—» в observed-E (линия не реально детектирована), но
показывают library_E + «присутствует по цепочке». Пользователь:
«перечислять не найденные линии не нужно».

Fix: убрать chain_expected на уровне data-pipeline:
  1. `_build_rows` больше НЕ инжектит F-111 chain-equilibrium placeholder
     rows (раньше делал при th232_dominant / u238_dominant).
  2. `_sync_peaks_rows_detail` пропускает orphan-row для peaks с
     `feature_kind == "chain_completeness"` — chart markers остаются,
     таблица не загрязняется.
  3. JS template (`interactive_v1_17_2.html`) больше не имеет
     'chain_expected' в `SECTION_ORDER` / `SECTION_LABELS`.

Acceptance: rowsSample / rowsBg payload от `_build_rows` НЕ содержит
ни одной строки с `section == "chain_expected"`, а также не содержит
старого комментария «присутствует по цепочке» / «присутствует по
равновесию» в поле `a`/`cmt`.
"""
from __future__ import annotations

import sys

sys.path.insert(0, "scripts")

from gamma.reporting.interactive_html import (  # noqa: E402
    _build_rows,
    _sync_peaks_rows_detail,
)


# ──── Stubs ────────────────────────────────────────────────────────


class _FakeAnalysisResult:
    def __init__(self, activities=None, sample_mass_kg=1.6):
        self.activities = activities or []
        self.sample_mass_kg = sample_mass_kg


# ──── Tests ────────────────────────────────────────────────────────


def test_bug23_th232_dominant_no_chain_expected_in_rows():
    """При th232_dominant — sample table НЕ должна содержать chain_expected
    placeholder-строк. Раньше Th-232 chain library (Ac-228 209.3, 277.4, 463,
    763 и т.д.) добавлялась в `_build_rows` с «— / E_lib» и «присутствует по
    цепочке»."""
    report = {
        "primary_feps": [
            {"peak_E_keV": 238.6, "peak_area_counts": 5000,
             "nuclide": "Pb-212", "library_E_keV": 238.6,
             "library_I_pct": 43.6},
        ],
        "secondary_peaks": [],
        "diagnostics": {"chain_dominance": {"th232_dominant": True}},
        "identified_nuclides": [],
    }
    rows = _build_rows(report, _FakeAnalysisResult())
    bad = [r for r in rows if r.get("section") == "chain_expected"]
    assert not bad, (
        "BUG-23: th232_dominant case produced {} chain_expected row(s); "
        "expected 0. Rows: {!r}".format(len(bad), bad)
    )
    # Defensive: also check no row has the legacy comment text
    legacy_markers = ("присутствует по цепочке", "присутствует по равновесию")
    for r in rows:
        for fld in ("a", "cmt"):
            v = r.get(fld, "") or ""
            for m in legacy_markers:
                assert m not in v, (
                    "BUG-23: legacy text {!r} leaked into row.{} = {!r}"
                    .format(m, fld, v)
                )


def test_bug23_u238_dominant_no_chain_expected_in_rows():
    """При u238_dominant — sample table НЕ должна содержать chain_expected.
    Раньше добавлялись Pb-214 295.2 / 351.9 и Bi-214 609.3 / 1764.5."""
    report = {
        "primary_feps": [],
        "secondary_peaks": [],
        "diagnostics": {"chain_dominance": {"u238_dominant": True}},
        "identified_nuclides": [],
    }
    rows = _build_rows(report, _FakeAnalysisResult())
    bad = [r for r in rows if r.get("section") == "chain_expected"]
    assert not bad, (
        "BUG-23: u238_dominant produced {} chain_expected row(s); expected "
        "0. Rows: {!r}".format(len(bad), bad)
    )


def test_bug23_bg_view_no_chain_expected():
    """Bg-view rowsBg payload (background_primary_feps) тоже не должен
    содержать chain_expected (фон передаётся с пустыми diagnostics, см.
    _generate_html bg_report_view, но проверим явно)."""
    bg_report_view = {
        "primary_feps": [
            {"peak_E_keV": 1461.0, "peak_area_counts": 8000,
             "nuclide": "K-40", "library_E_keV": 1460.8,
             "library_I_pct": 10.7},
        ],
        "secondary_peaks": [],
        "diagnostics": {},
        "identified_nuclides": [],
    }
    rows_bg = _build_rows(bg_report_view, _FakeAnalysisResult())
    bad = [r for r in rows_bg if r.get("section") == "chain_expected"]
    assert not bad, (
        "BUG-23: bg-view produced {} chain_expected row(s); expected 0"
        .format(len(bad))
    )


def test_bug23_sync_skips_chain_completeness_orphan_peaks():
    """`_sync_peaks_rows_detail` не должен создавать табличную строку для
    peak с feature_kind == 'chain_completeness'. Chart marker / detail
    остаются, но строки в таблице нет."""
    # Orphan chain-completeness peak: in peaks/detail, but not in rows.
    pid = "p209"
    peaks = [{
        "id": pid, "e": 209.3, "label": "Ac-228 209",
        "color": "#5b7fff",
        "feature_kind": "chain_completeness",
        "is_secondary": False,
    }]
    rows = []
    detail = {pid: {"title": "Ac-228 209", "color": "#5b7fff",
                    "desc": "Подтверждена по равновесию цепочки Th-232"}}
    _sync_peaks_rows_detail(peaks, rows, detail)
    # No row should have been created for the chain_completeness peak.
    chain_rows = [r for r in rows if r.get("peak") == pid]
    assert not chain_rows, (
        "BUG-23: sync created {} row(s) for chain_completeness peak; "
        "expected 0. Rows: {!r}".format(len(chain_rows), chain_rows)
    )
    # And of course no row should have section == chain_expected.
    bad = [r for r in rows if r.get("section") == "chain_expected"]
    assert not bad, (
        "BUG-23: sync emitted chain_expected section; got {!r}".format(bad)
    )


def test_bug23_sync_still_creates_rows_for_non_chain_orphan_peaks():
    """Regression guard: для non-chain orphan peaks (например secondary с
    feature_kind != chain_completeness) sync ВСЁ ЕЩЁ должен дополнять
    таблицу — иначе peaks↔rows invariant сломается."""
    pid = "p1593"
    peaks = [{
        "id": pid, "e": 1593.0, "label": "DE Tl-208 1593",
        "color": "#ff9900",
        "feature_kind": "double_escape",
        "is_secondary": True,
    }]
    rows = []
    detail = {pid: {"title": "DE Tl-208 1593", "color": "#ff9900",
                    "desc": "Двойное вылет"}}
    _sync_peaks_rows_detail(peaks, rows, detail)
    # A row SHOULD have been created (non-chain-completeness path).
    created = [r for r in rows if r.get("peak") == pid]
    assert len(created) == 1, (
        "BUG-23: sync should still create row for double_escape orphan "
        "peak (invariant peaks↔rows). Rows: {!r}".format(rows)
    )
    # And it must not be tagged chain_expected.
    assert created[0].get("section") != "chain_expected", (
        "BUG-23: orphan-row from sync must not use chain_expected section"
    )


if __name__ == "__main__":
    test_bug23_th232_dominant_no_chain_expected_in_rows()
    test_bug23_u238_dominant_no_chain_expected_in_rows()
    test_bug23_bg_view_no_chain_expected()
    test_bug23_sync_skips_chain_completeness_orphan_peaks()
    test_bug23_sync_still_creates_rows_for_non_chain_orphan_peaks()
    print("OK")
