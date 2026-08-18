"""BUG-15 — Shared-peak deduplication in compute_activity.

Контракты:

  Within-nuclide (compute_activity внутренний dedup):
      Если две библиотечные линии ОДНОГО нуклида попадают в один и тот
      же peak_channel с одинаковой peak_area > 0, в weighted-mean вкладывает
      только линия с максимальным library_I_pct. Остальные уходят в
      ``lines_skipped`` с reason="shared_peak_dedupe …".

      Физика: связанная подгонка multiplet'а присваивает одинаковую
      площадь S всем компонентам, чья ширина не позволяет разделить
      интенсивности независимо. Засчитывать каждой линии полную S — это
      double-counting; A_i ∝ 1/I_pct, поэтому слабая (low-I_pct) линия
      даёт раздутую активность.

  Cross-nuclide (compute_activities_for_all внешний dedup):
      Если на одном peak_channel есть линии разных нуклидов И хотя бы
      один из них имеет ``is_characteristic=True`` на этом канале, его
      нуклид становится «владельцем» канала. Линии остальных нуклидов
      на том же канале (false-match через FWHM window) skipped с
      reason="cross_nuclide_peak_owned …".

Bug history: на Th-232 demo Ac-228 ранее давал weighted-mean ~6167 Bq/kg
(+218% от LSRM 1940 Bq/kg). Корни:
  1. Within-nuclide: 904.20 keV (I=0.77%) получал ту же area, что 911.20
     keV (I=25.8%), и давал A_i ≈ 73 000 Bq.
  2. Cross-nuclide: Ac-228 562.50 keV (I=0.87%) матчилось в окно Tl-208
     характеристического 583 keV — A_i ≈ 68 000 Bq.
После фикса: Ac-228 ≈ 2240 Bq/kg, отношение к Pb-212 (single-line, 3156)
≈ 0.71 — внутри допустимого ±30%.
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Tuple

import pytest

sys.path.insert(
    0, str(Path(__file__).resolve().parent.parent.parent / "scripts")
)

from gamma.activity.compute import (
    compute_activity, compute_activities_for_all,
)
from gamma.calibration.efficiency import EfficiencyCurve


# ─── Helpers ────────────────────────────────────────────────────────────


def _flat_efficiency_curve(eps: float = 1.0e-3) -> EfficiencyCurve:
    """ε(E) = eps everywhere (log-log poly degree 0, coef = log eps).
    Конструктивно гарантирует, что A_i зависит ТОЛЬКО от (S, I, t_live),
    не от энергетической зависимости эффективности.
    """
    return EfficiencyCurve(
        coefficients=(math.log(eps),),
        E_min_keV=1.0,
        E_max_keV=3000.0,
        chi2_per_dof=1.0,
        n_points_used=10,
        n_dof=9,
        detector_id="test",
        geometry="test",
    )


@dataclass
class _FakeMatch:
    """Duck-typed stand-in for LineMatch (writable so tests can mutate)."""
    nuclide: str
    library_E_keV: float
    library_I_pct: float
    peak_channel: int
    peak_E_keV: float
    peak_sigma: float = 1.0
    residual_keV: float = 0.0
    is_characteristic: bool = False
    peak_area: float = 0.0
    peak_area_uncertainty: float = 0.0
    peak_area_source: str = "deconvolved_coupled"


@dataclass
class _FakeNuclideId:
    """Duck-typed stand-in for NuclideIdentification."""
    nuclide: str
    detected: bool = True
    reason: str = "test"
    characteristic_line_keV: float = 0.0
    matched_lines: tuple = ()
    confidence: object = None


@dataclass
class _FakeIdResult:
    """Duck-typed stand-in for IdentificationResult."""
    detector_type: str = "test"
    window: object = None
    candidates_considered: int = 0
    detected_nuclides: tuple = ()
    rejected_nuclides: tuple = ()
    unmatched_peaks: tuple = ()
    notes: str = ""


# ─── Within-nuclide dedup ────────────────────────────────────────────────


def test_within_nuclide_shared_peak_dedup_keeps_highest_I():
    """Две линии Cs-137 с одинаковым peak_channel и одной peak_area:
    оставляется только I_pct-доминантная; вторая уходит в skipped."""
    # Cs-137 661.66 keV — характеристическая линия, I=85.1%.
    # Фиктивная вторая «линия» 660.00 keV с I=2.0% попадает в тот же
    # peak_channel (NaI не разделяет → coupled fitter присваивает одинаковую S).
    S = 50_000.0
    m_dom = _FakeMatch(
        nuclide="Cs-137",
        library_E_keV=661.66, library_I_pct=85.1,
        peak_channel=300, peak_E_keV=661.5,
        peak_area=S, peak_area_uncertainty=math.sqrt(S),
        is_characteristic=True,
    )
    m_weak = _FakeMatch(
        nuclide="Cs-137",
        library_E_keV=660.00, library_I_pct=2.0,
        peak_channel=300, peak_E_keV=661.5,
        peak_area=S, peak_area_uncertainty=math.sqrt(S),
        is_characteristic=False,
    )
    nid = _FakeNuclideId(
        nuclide="Cs-137",
        characteristic_line_keV=661.66,
        matched_lines=(m_dom, m_weak),
    )

    res = compute_activity(
        nid,
        efficiency_curve=_flat_efficiency_curve(1.0e-3),
        live_time_s=1000.0,
        from_bg_subtracted=True,
        decay_correction=False,
    )

    # Только одна линия использована (доминантная по I_pct).
    assert res.n_lines_used() == 1, (
        f"expected 1 line after within-nuclide dedup, got "
        f"{res.n_lines_used()}; lines_used="
        f"{[la.E_keV for la in res.lines_used]}"
    )
    used = res.lines_used[0]
    assert math.isclose(used.E_keV, 661.66, abs_tol=0.01), (
        f"wrong line kept: {used.E_keV}, expected 661.66"
    )

    # Слабая линия в lines_skipped с правильной причиной.
    skipped_E = [E for E, _ in res.lines_skipped]
    assert any(math.isclose(E, 660.00, abs_tol=0.01) for E in skipped_E), (
        f"weak line 660.00 should be in skipped, got {res.lines_skipped}"
    )
    skipped_reasons = [reason for _, reason in res.lines_skipped]
    assert any("shared_peak_dedupe" in r for r in skipped_reasons), (
        f"expected 'shared_peak_dedupe' reason, got {skipped_reasons}"
    )


def test_within_nuclide_dedup_prevents_double_counting_A_i():
    """Контракт BUG-15: если weak-I_pct линия не была бы отброшена,
    weighted-mean уплыл бы наверх на десятки. После dedup A_Bq ровно
    совпадает с A_i из одной dominant-I_pct линии."""
    S = 100_000.0
    t_live = 1000.0
    eps = 1.0e-3
    I_dom = 25.8     # ~Ac-228 911.20 keV
    I_weak = 0.77    # ~Ac-228 904.20 keV

    m_dom = _FakeMatch(
        nuclide="Ac-228",
        library_E_keV=911.20, library_I_pct=I_dom,
        peak_channel=350, peak_E_keV=911.0,
        peak_area=S, peak_area_uncertainty=math.sqrt(S),
        is_characteristic=True,
    )
    m_weak = _FakeMatch(
        nuclide="Ac-228",
        library_E_keV=904.20, library_I_pct=I_weak,
        peak_channel=350, peak_E_keV=911.0,
        peak_area=S, peak_area_uncertainty=math.sqrt(S),
    )
    nid = _FakeNuclideId(
        nuclide="Ac-228",
        characteristic_line_keV=911.20,
        matched_lines=(m_dom, m_weak),
    )

    res = compute_activity(
        nid, efficiency_curve=_flat_efficiency_curve(eps),
        live_time_s=t_live, from_bg_subtracted=True,
        decay_correction=False,
    )

    # BUG-15 (v1.18.31): keeps only the dominant-I line, throws away weak.
    # BUG-27 (v1.18.32+): the dominant-I survivor's effective I in the
    # A_i denominator is the SUM I_dom + I_weak (in-blob ΣI), because S
    # physically contains photons from both unresolved library lines.
    #   A = S / (ε · (I_dom + I_weak)/100 · t_live)
    # Numerically: 100000 / (1e-3 · 0.26570 · 1000) = 376 364 Bq
    # (a 2.9% drop from the bare-I_dom value of 387 597 Bq).
    A_expected_bug27 = S / (eps * ((I_dom + I_weak) / 100.0) * t_live)
    A_dom_only = S / (eps * (I_dom / 100.0) * t_live)
    assert math.isclose(res.A_Bq, A_expected_bug27, rel_tol=1e-3), (
        f"A_Bq={res.A_Bq:.1f}, BUG-27 expected (ΣI denom) "
        f"~{A_expected_bug27:.1f}, bare-I_dom would be ~{A_dom_only:.1f}"
    )

    # Контр-проверка: без dedup'а weak-линия даёт A_i ~ 33x больше.
    # Если бы оба A_i учитывались независимо (БЕЗ BUG-15), weighted-mean
    # смещался бы вверх.
    A_weak_would_be = S / (eps * (I_weak / 100.0) * t_live)
    assert A_weak_would_be > 20.0 * A_dom_only, (
        "test sanity check failed: weak/dom A_i ratio is small, "
        "BUG-15 not exercised"
    )
    # Финальный A_Bq не должен взлететь до средневзвешенного с weak-линией.
    assert res.A_Bq < 2.0 * A_dom_only, (
        f"weighted mean drifted upward — within-nuclide dedup not applied: "
        f"A_Bq={res.A_Bq:.0f}, A_dom_only={A_dom_only:.0f}"
    )


# ─── Cross-nuclide ownership ─────────────────────────────────────────────


def test_cross_nuclide_characteristic_owner_skips_non_owners():
    """Если на peak_channel C нуклид A имеет characteristic line, а
    нуклид B имеет ту же peak_channel C через FWHM window matching,
    то линия B на C отбрасывается в compute_activities_for_all
    через cross_nuclide_skip_energies_keV."""
    S = 200_000.0
    eps_curve = _flat_efficiency_curve(1.0e-3)

    # Tl-208 583 keV — characteristic peak (ch=203), area S.
    tl_char = _FakeMatch(
        nuclide="Tl-208",
        library_E_keV=583.19, library_I_pct=84.5,
        peak_channel=203, peak_E_keV=583.0,
        peak_area=S, peak_area_uncertainty=math.sqrt(S),
        is_characteristic=True,
    )
    tl_id = _FakeNuclideId(
        nuclide="Tl-208",
        characteristic_line_keV=583.19,
        matched_lines=(tl_char,),
    )

    # Ac-228 — characteristic: 911 keV (ch=350). Кроме того, library line
    # 562.50 keV (I=0.87%) случайно матчится в FWHM-окно Tl-208 583 (ch=203),
    # получая ту же area S (false cross-match).
    ac_char = _FakeMatch(
        nuclide="Ac-228",
        library_E_keV=911.20, library_I_pct=25.8,
        peak_channel=350, peak_E_keV=911.0,
        peak_area=S, peak_area_uncertainty=math.sqrt(S),
        is_characteristic=True,
    )
    ac_false = _FakeMatch(
        nuclide="Ac-228",
        library_E_keV=562.50, library_I_pct=0.87,
        peak_channel=203, peak_E_keV=583.0,  # упало на Tl-208 583
        peak_area=S, peak_area_uncertainty=math.sqrt(S),
        is_characteristic=False,
    )
    ac_id = _FakeNuclideId(
        nuclide="Ac-228",
        characteristic_line_keV=911.20,
        matched_lines=(ac_char, ac_false),
    )

    id_result = _FakeIdResult(detected_nuclides=(tl_id, ac_id))

    results = compute_activities_for_all(
        id_result,
        efficiency_curve=eps_curve,
        live_time_s=1000.0,
        from_bg_subtracted=True,
        decay_correction=False,
    )

    by_nuc = {r.nuclide: r for r in results}
    assert "Ac-228" in by_nuc and "Tl-208" in by_nuc

    ac_res = by_nuc["Ac-228"]
    # 562.50 keV должна быть skipped с cross_nuclide_peak_owned reason.
    skipped_562 = [
        (E, r) for (E, r) in ac_res.lines_skipped
        if math.isclose(E, 562.50, abs_tol=0.01)
    ]
    assert skipped_562, (
        f"Ac-228 562.50 keV должна быть в skipped (cross-nuclide), "
        f"got skipped={ac_res.lines_skipped}, "
        f"used={[la.E_keV for la in ac_res.lines_used]}"
    )
    assert "cross_nuclide_peak_owned" in skipped_562[0][1], (
        f"wrong skip reason for 562.50: {skipped_562[0][1]}"
    )

    # Tl-208 (владелец канала) сохраняет 583.19.
    tl_res = by_nuc["Tl-208"]
    used_E = [la.E_keV for la in tl_res.lines_used]
    assert any(math.isclose(E, 583.19, abs_tol=0.01) for E in used_E), (
        f"Tl-208 owner должен сохранить 583.19; lines_used={used_E}"
    )


def test_cross_nuclide_no_action_when_only_one_nuclide_on_channel():
    """Если на peak_channel матчится только один нуклид, никакого
    cross-nuclide skip не должно происходить (даже если ничего не
    is_characteristic)."""
    S = 100_000.0
    m = _FakeMatch(
        nuclide="Cs-137",
        library_E_keV=661.66, library_I_pct=85.1,
        peak_channel=300, peak_E_keV=661.5,
        peak_area=S, peak_area_uncertainty=math.sqrt(S),
        is_characteristic=False,
    )
    nid = _FakeNuclideId(
        nuclide="Cs-137",
        characteristic_line_keV=661.66,
        matched_lines=(m,),
    )
    id_result = _FakeIdResult(detected_nuclides=(nid,))

    results = compute_activities_for_all(
        id_result,
        efficiency_curve=_flat_efficiency_curve(1.0e-3),
        live_time_s=1000.0,
        from_bg_subtracted=True,
        decay_correction=False,
    )
    assert len(results) == 1
    r = results[0]
    assert r.n_lines_used() == 1, (
        f"single-nuclide / single-channel must not be skipped, "
        f"got skipped={r.lines_skipped}"
    )


if __name__ == "__main__":
    test_within_nuclide_shared_peak_dedup_keeps_highest_I()
    test_within_nuclide_dedup_prevents_double_counting_A_i()
    test_cross_nuclide_characteristic_owner_skips_non_owners()
    test_cross_nuclide_no_action_when_only_one_nuclide_on_channel()
    print("BUG-15 shared-peak dedup tests PASS.")
