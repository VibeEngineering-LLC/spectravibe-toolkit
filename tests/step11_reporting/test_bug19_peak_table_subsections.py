"""BUG-19 / v1.18.30+ (Agent B) — sample peak table must split into distinct
subsections, and per-line activity replaces nuclide weighted mean.

BUG-23 / v1.18.31+ (Agent B): subsection «chain_expected» (F-111 placeholders
из цепочки по равновесию) удалена — пользователь не хочет видеть в таблице
не-детектированные линии. Тест обновлён с 4 секций до 3:
  • primary_detected
  • weak_candidate
  • secondary

Symptom (BUG-19): таблица «Пики образца» сводила в одну колонку разных
entry-type:
  • реально детектированные FEP с подгонкой,
  • placeholder'ы по равновесию цепочки (F-111) — удалены в BUG-23,
  • secondary processes (Comp/SE/DE/sum/backscatter),
  • слабые библиотечные кандидаты.

И каждая строка Ac-228 (911 / 964 / 969 кэВ etc.) показывала ОДИН и тот же
weighted-mean `A = 6167 ± 1026 Бк/кг`, хотя per-line A_i физически разные
(они вычисляются для каждой линии и потом усредняются).

Fix (BUG-19 + BUG-23):
  1. `_build_rows` помечает каждую строку полем `section ∈ {primary_detected,
     weak_candidate, secondary}`. Section "chain_expected" больше не
     производится (BUG-23).
  2. Для chain-нуклидов с известными `analysis_result.activities[*].lines_used`
     поле `a` подставляется per-line A_i ± σ_A_i.
  3. JS-renderer (renderRows) рендерит каждую подсекцию отдельным
     заголовком (colspan=5, class="fp-section-head").
"""
from __future__ import annotations

import sys

sys.path.insert(0, "scripts")

from gamma.reporting.interactive_html import _build_rows  # noqa: E402


# ──── Minimal stubs ────────────────────────────────────────────────


class _FakeLineActivity:
    """Duck-type for gamma.activity.compute.LineActivity."""

    def __init__(self, E_keV: float, A_Bq: float, sigma_A_Bq: float):
        self.E_keV = E_keV
        self.A_Bq = A_Bq
        self.sigma_A_Bq = sigma_A_Bq


class _FakeActivityResult:
    def __init__(self, nuclide: str, lines_used):
        self.nuclide = nuclide
        self.lines_used = tuple(lines_used)


class _FakeAnalysisResult:
    def __init__(self, activities, sample_mass_kg=1.6):
        self.activities = activities
        self.sample_mass_kg = sample_mass_kg


# ──── Tests ────────────────────────────────────────────────────────


def test_bug19_every_row_has_section_field():
    """Каждая строка должна получить `section` — нет «голых» rows
    (старый contract сломан, JS использует section для group-render)."""
    report = {
        "primary_feps": [
            {"peak_E_keV": 238.6, "peak_area_counts": 5000,
             "nuclide": "Pb-212", "library_E_keV": 238.6,
             "library_I_pct": 43.6},
        ],
        "secondary_peaks": [
            {"energy_keV": 2103.0, "type": "single_escape",
             "feature_kind": "single_escape",
             "parent_nuclide": "Tl-208", "significance": 4.0},
        ],
        "diagnostics": {"chain_dominance": {"th232_dominant": True}},
        "identified_nuclides": [],
    }
    ar = _FakeAnalysisResult(activities=[])
    rows = _build_rows(report, ar)
    for r in rows:
        assert "section" in r, (
            "BUG-19: every row must have `section` field; got row={!r}"
            .format(r)
        )


def test_bug19_section_categorization():
    """primary_feps → 'primary_detected'; secondary → 'secondary'.

    BUG-23 / v1.18.31+: section 'chain_expected' больше НЕ производится —
    F-111 chain-equilibrium placeholders удалены из таблицы (по фидбэку
    пользователя «перечислять не найденные линии не нужно»). Тест явно
    проверяет, что при `th232_dominant: True` никакая строка не имеет
    section == 'chain_expected'."""
    report = {
        "primary_feps": [
            {"peak_E_keV": 238.6, "peak_area_counts": 5000,
             "nuclide": "Pb-212", "library_E_keV": 238.6,
             "library_I_pct": 43.6},
        ],
        "secondary_peaks": [
            {"energy_keV": 1593.0, "type": "double_escape",
             "feature_kind": "double_escape",
             "parent_nuclide": "Tl-208", "significance": 4.0},
        ],
        "diagnostics": {"chain_dominance": {"th232_dominant": True}},
        "identified_nuclides": [],
    }
    rows = _build_rows(report, _FakeAnalysisResult(activities=[]))
    by_section = {}
    for r in rows:
        by_section.setdefault(r["section"], []).append(r)
    assert "primary_detected" in by_section, "missing primary_detected section"
    assert any(r["iso"] == "Pb-212" for r in by_section["primary_detected"])
    assert "secondary" in by_section, "missing secondary section"
    assert any(
        "Tl-208" in r.get("iso", "") or r.get("kls") == "fp-phys"
        for r in by_section["secondary"]
    )
    # BUG-23: chain_expected must NEVER appear, even when th232_dominant=True.
    assert "chain_expected" not in by_section, (
        "BUG-23: section 'chain_expected' must not be produced; got rows "
        "{!r}".format(by_section.get("chain_expected"))
    )
    # Allowed sections: only the 3 documented ones.
    assert set(by_section.keys()) <= {"primary_detected",
                                       "weak_candidate",
                                       "secondary"}, (
        "BUG-23: unexpected sections present; got {!r}"
        .format(set(by_section.keys()))
    )


def test_bug19_per_line_activity_overrides_weighted_mean():
    """Для Ac-228 911/964/969 line_activities дают разные A_i — каждая
    строка получает СВОЁ значение, не общий weighted-mean."""
    line_acts = [
        _FakeLineActivity(E_keV=911.2, A_Bq=10000.0, sigma_A_Bq=500.0),
        _FakeLineActivity(E_keV=964.77, A_Bq=8500.0, sigma_A_Bq=900.0),
        _FakeLineActivity(E_keV=968.97, A_Bq=12000.0, sigma_A_Bq=700.0),
    ]
    ar_ac228 = _FakeActivityResult(nuclide="Ac-228", lines_used=line_acts)
    ar = _FakeAnalysisResult(activities=[ar_ac228], sample_mass_kg=1.6)

    report = {
        "primary_feps": [
            {"peak_E_keV": 911.0, "peak_area_counts": 5000,
             "nuclide": "Ac-228", "library_E_keV": 911.2,
             "library_I_pct": 25.8},
            {"peak_E_keV": 964.5, "peak_area_counts": 3000,
             "nuclide": "Ac-228", "library_E_keV": 964.77,
             "library_I_pct": 4.99},
            {"peak_E_keV": 968.8, "peak_area_counts": 4000,
             "nuclide": "Ac-228", "library_E_keV": 968.97,
             "library_I_pct": 15.8},
        ],
        "secondary_peaks": [],
        "diagnostics": {},
        # Existing `identified_nuclides` would give weighted-mean as fallback,
        # but per-line override must take priority.
        "identified_nuclides": [{
            "nuclide": "Ac-228",
            "specific_activity_Bq_per_kg": 6167.0,
            "specific_activity_sigma_Bq_per_kg": 1026.0,
        }],
    }
    rows = _build_rows(report, ar)
    by_lib_e = {}
    for r in rows:
        if r.get("iso") != "Ac-228":
            continue
        # extract library E (last number in "911.0 / 911.2" → 911.2)
        import re
        m = re.search(r"(\d+(?:\.\d+)?)\s*$", r["line"])
        if m:
            by_lib_e[round(float(m.group(1)), 1)] = r["a"]
    # 911.2 → A=10000/1.6=6250, σ=312 — should be different from 964/969
    a911 = by_lib_e.get(911.2, "")
    a964 = by_lib_e.get(964.8, "")  # rounded to 1 decimal
    a969 = by_lib_e.get(969.0, "")
    # Anti-flake: at least 2 of the 3 should match the per-line lookup
    # (rounding 964.77→964.8 happens via _peak_id → row line text)
    assert "6250" in a911, (
        "BUG-19: Ac-228 911 must show per-line A_i (10000/1.6=6250), not "
        "weighted-mean 6167; got `a`={!r}".format(a911)
    )
    # Not all three strings must equal — that's the core fix.
    distinct = {a911, a964, a969} - {""}
    assert len(distinct) >= 2, (
        "BUG-19: Ac-228 per-line activities must differ across lines "
        "(not all = weighted-mean 6167±1026); got distinct values={!r}"
        .format(distinct)
    )


def test_bug19_per_line_falls_back_when_no_match():
    """Если для (nuclide, lib_E) нет line_activity match, остаётся
    weighted-mean specific_activity (backward-compat)."""
    ar = _FakeAnalysisResult(activities=[])   # no line_activities at all
    report = {
        "primary_feps": [
            {"peak_E_keV": 238.6, "peak_area_counts": 5000,
             "nuclide": "Pb-212", "library_E_keV": 238.6,
             "library_I_pct": 43.6},
        ],
        "secondary_peaks": [],
        "diagnostics": {},
        "identified_nuclides": [{
            "nuclide": "Pb-212",
            "specific_activity_Bq_per_kg": 250.0,
            "specific_activity_sigma_Bq_per_kg": 30.0,
        }],
    }
    rows = _build_rows(report, ar)
    pb_rows = [r for r in rows if r.get("iso") == "Pb-212"]
    assert pb_rows, "expected Pb-212 row"
    # weighted-mean fallback: "250 ± 30"
    assert "250" in pb_rows[0]["a"] and "30" in pb_rows[0]["a"], (
        "BUG-19: backward-compat — when line_activities is empty, A column "
        "must fall back to specific_activity_Bq_per_kg ± σ; got `a`={!r}"
        .format(pb_rows[0]["a"])
    )


if __name__ == "__main__":
    test_bug19_every_row_has_section_field()
    test_bug19_section_categorization()
    test_bug19_per_line_activity_overrides_weighted_mean()
    test_bug19_per_line_falls_back_when_no_match()
    print("OK")
