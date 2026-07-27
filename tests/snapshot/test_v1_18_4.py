# -*- coding: utf-8 -*-
"""v1.18.4 — Quasi-template activation: high-level solver wrapper."""
from __future__ import annotations
import math, os, sys
import pytest

SCRIPTS = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts")
)
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)


def _ch_to_keV(ch):
    return float(ch) * 1.0   # 1 keV/channel


def _fwhm_at(E):
    return 0.07 * math.sqrt(max(E, 1.0) * 662.0)


def _eps_at(E):
    if E <= 0:
        return None
    return 0.05 * (662.0 / E) ** 0.8


def _build_synthetic_spectrum(activities, t_live=600.0, n_ch=2048):
    """Construct synthetic spectrum from known activities of given nuclides."""
    from gamma.activity.quasi_template_ppp import (
        NuclideDef, NuclideLine, build_nuclide_template,
    )
    from gamma.data.nuclide_library import get_nuclide
    counts = [0.0] * n_ch
    for nid, A_Bq in activities.items():
        rec = get_nuclide(nid)
        if not rec:
            continue
        lines = []
        for ll in rec.get("lines", []):
            E = float(ll[0])
            I_pct = float(ll[1])
            eps = _eps_at(E)
            if eps and eps > 0 and I_pct > 0:
                lines.append(NuclideLine(
                    E_keV=E, intensity=I_pct / 100.0, efficiency=eps,
                ))
        if not lines:
            continue
        ndef = NuclideDef(nuclide_id=nid, lines=lines)
        tpl = build_nuclide_template(ndef, n_ch, _ch_to_keV, _fwhm_at)
        for ch in range(n_ch):
            counts[ch] += A_Bq * t_live * tpl.counts[ch]
    return counts


# ──────────────────────────────────────────────────────────────────
# Smoke / API
# ──────────────────────────────────────────────────────────────────

def test_F18_4_solver_returns_list_for_known_nuclides():
    from gamma.activity.quasi_template_solver import (
        solve_quasi_template_activities,
    )
    obs = [0.0] * 2048
    res = solve_quasi_template_activities(
        spectrum_counts=obs,
        channel_to_keV=_ch_to_keV,
        fwhm_at_E_func=_fwhm_at,
        efficiency_at_E_func=_eps_at,
        nuclide_ids=["Cs-137"],
        live_time_s=600.0,
    )
    # Если Cs-137 есть в library — returns list. Если нет — пустой list.
    assert isinstance(res, list)


def test_F18_4_zero_live_time_raises():
    from gamma.activity.quasi_template_solver import (
        solve_quasi_template_activities,
    )
    with pytest.raises(ValueError):
        solve_quasi_template_activities(
            spectrum_counts=[0.0] * 100,
            channel_to_keV=_ch_to_keV,
            fwhm_at_E_func=_fwhm_at,
            efficiency_at_E_func=_eps_at,
            nuclide_ids=["Cs-137"],
            live_time_s=0.0,
        )


def test_F18_4_unknown_nuclides_skipped():
    """Unknown nuclide IDs не должны падать — skipped silently."""
    from gamma.activity.quasi_template_solver import (
        solve_quasi_template_activities,
    )
    obs = [0.0] * 100
    res = solve_quasi_template_activities(
        spectrum_counts=obs,
        channel_to_keV=_ch_to_keV,
        fwhm_at_E_func=_fwhm_at,
        efficiency_at_E_func=_eps_at,
        nuclide_ids=["NOT-A-REAL-NUCLIDE-9999"],
        live_time_s=600.0,
    )
    assert res == []


# ──────────────────────────────────────────────────────────────────
# Activity recovery — synthetic round-trip
# ──────────────────────────────────────────────────────────────────

def test_F18_4_recovers_cs137_activity_from_synthetic():
    """Известная активность Cs-137 → должна быть восстановлена."""
    from gamma.activity.quasi_template_solver import (
        solve_quasi_template_activities,
    )
    from gamma.data.nuclide_library import get_nuclide
    if not get_nuclide("Cs-137"):
        pytest.skip("Cs-137 not in library")
    A_true = 1000.0
    t = 600.0
    obs = _build_synthetic_spectrum({"Cs-137": A_true}, t_live=t, n_ch=2048)
    # Disable continuum чтобы избежать шаблон-расхождения
    res = solve_quasi_template_activities(
        spectrum_counts=obs, channel_to_keV=_ch_to_keV,
        fwhm_at_E_func=_fwhm_at, efficiency_at_E_func=_eps_at,
        nuclide_ids=["Cs-137"], live_time_s=t,
        enable_compton_continuum=False,
    )
    assert len(res) == 1
    assert res[0].sigma_method == "quasi_template"
    assert res[0].A_Bq == pytest.approx(A_true, rel=0.05)
    assert "F-302..F-304" in (res[0].notes or "")
    assert "χ²_red=" in (res[0].notes or "")


def test_F18_4_recovers_two_nuclides():
    """Cs-137 + Co-60 synthetic → both restored."""
    from gamma.activity.quasi_template_solver import (
        solve_quasi_template_activities,
    )
    from gamma.data.nuclide_library import get_nuclide
    if not get_nuclide("Cs-137") or not get_nuclide("Co-60"):
        pytest.skip("Cs-137 or Co-60 not in library")
    activities = {"Cs-137": 500.0, "Co-60": 200.0}
    t = 1000.0
    obs = _build_synthetic_spectrum(activities, t_live=t, n_ch=2048)
    res = solve_quasi_template_activities(
        spectrum_counts=obs, channel_to_keV=_ch_to_keV,
        fwhm_at_E_func=_fwhm_at, efficiency_at_E_func=_eps_at,
        nuclide_ids=["Cs-137", "Co-60"], live_time_s=t,
        enable_compton_continuum=False,
    )
    by_name = {r.nuclide: r for r in res}
    assert by_name["Cs-137"].A_Bq == pytest.approx(500.0, rel=0.05)
    assert by_name["Co-60"].A_Bq == pytest.approx(200.0, rel=0.05)


def test_F18_4_notes_contain_chi2_and_acceptance():
    from gamma.activity.quasi_template_solver import (
        solve_quasi_template_activities,
    )
    from gamma.data.nuclide_library import get_nuclide
    if not get_nuclide("Cs-137"):
        pytest.skip("Cs-137 not in library")
    obs = _build_synthetic_spectrum(
        {"Cs-137": 500.0}, t_live=600.0, n_ch=2048,
    )
    res = solve_quasi_template_activities(
        spectrum_counts=obs, channel_to_keV=_ch_to_keV,
        fwhm_at_E_func=_fwhm_at, efficiency_at_E_func=_eps_at,
        nuclide_ids=["Cs-137"], live_time_s=600.0,
        enable_compton_continuum=False,
    )
    assert len(res) == 1
    notes = res[0].notes or ""
    assert "χ²_red=" in notes
    assert "is_accepted=" in notes
    assert "detector=" in notes


def test_F18_4_compton_continuum_toggle():
    """С enable_compton_continuum=True vs False должно влиять на fit (хотя бы notes)."""
    from gamma.activity.quasi_template_solver import (
        solve_quasi_template_activities,
    )
    from gamma.data.nuclide_library import get_nuclide
    if not get_nuclide("Cs-137"):
        pytest.skip("Cs-137 not in library")
    obs = _build_synthetic_spectrum(
        {"Cs-137": 500.0}, t_live=600.0, n_ch=2048,
    )
    res_no = solve_quasi_template_activities(
        spectrum_counts=obs, channel_to_keV=_ch_to_keV,
        fwhm_at_E_func=_fwhm_at, efficiency_at_E_func=_eps_at,
        nuclide_ids=["Cs-137"], live_time_s=600.0,
        enable_compton_continuum=False,
    )
    res_yes = solve_quasi_template_activities(
        spectrum_counts=obs, channel_to_keV=_ch_to_keV,
        fwhm_at_E_func=_fwhm_at, efficiency_at_E_func=_eps_at,
        nuclide_ids=["Cs-137"], live_time_s=600.0,
        enable_compton_continuum=True,
    )
    # Оба возвращают валидный fit
    assert len(res_no) == 1 and len(res_yes) == 1
    assert res_no[0].A_Bq > 0
    assert res_yes[0].A_Bq > 0


def test_F18_4_energy_window_restriction():
    """С energy_window — fit не должен падать; result list корректной длины."""
    from gamma.activity.quasi_template_solver import (
        solve_quasi_template_activities,
    )
    from gamma.data.nuclide_library import get_nuclide
    if not get_nuclide("Cs-137"):
        pytest.skip("Cs-137 not in library")
    obs = _build_synthetic_spectrum(
        {"Cs-137": 500.0}, t_live=600.0, n_ch=2048,
    )
    res = solve_quasi_template_activities(
        spectrum_counts=obs, channel_to_keV=_ch_to_keV,
        fwhm_at_E_func=_fwhm_at, efficiency_at_E_func=_eps_at,
        nuclide_ids=["Cs-137"], live_time_s=600.0,
        energy_window_keV=(500.0, 900.0),
        enable_compton_continuum=False,
    )
    assert len(res) == 1
    assert res[0].A_Bq > 0


def test_F18_4_result_is_activity_result_compatible():
    """Return type должен быть ActivityResult — для downstream совместимости."""
    from gamma.activity.quasi_template_solver import (
        solve_quasi_template_activities,
    )
    from gamma.activity.compute import ActivityResult
    from gamma.data.nuclide_library import get_nuclide
    if not get_nuclide("Cs-137"):
        pytest.skip("Cs-137 not in library")
    obs = _build_synthetic_spectrum(
        {"Cs-137": 500.0}, t_live=600.0, n_ch=2048,
    )
    res = solve_quasi_template_activities(
        spectrum_counts=obs, channel_to_keV=_ch_to_keV,
        fwhm_at_E_func=_fwhm_at, efficiency_at_E_func=_eps_at,
        nuclide_ids=["Cs-137"], live_time_s=600.0,
        enable_compton_continuum=False,
    )
    assert all(isinstance(r, ActivityResult) for r in res)
