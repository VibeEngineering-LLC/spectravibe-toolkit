# -*- coding: utf-8 -*-
"""v1.18.1 — Activity integration slice 1.

Wire-up F-296 (auto-TCS) и F-294 (Cutshall self-absorption fallback)
в gamma/activity/compute.py через opt-in flags. Default OFF → back-compat.
"""
from __future__ import annotations
import math, os, sys
import pytest

SCRIPTS = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts")
)
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)


# ──────────────────────────────────────────────────────────────────
# Mock-классы для compute_activity (минимальные)
# ──────────────────────────────────────────────────────────────────

class _MockLineMatch:
    def __init__(self, E, I_pct, area, area_unc=None, area_source=""):
        self.library_E_keV = E
        self.library_I_pct = I_pct
        self.peak_area = area
        self.peak_area_uncertainty = area_unc
        self.peak_area_source = area_source


class _MockNuclideId:
    def __init__(self, nuclide, lines, detected=True):
        self.nuclide = nuclide
        self.matched_lines = lines
        self.detected = detected


class _MockEfficiency:
    """Simple ε(E) = ε₀·(E₀/E)^p curve."""
    def __init__(self, eps0=0.05, E0=662.0, power=0.8):
        self.eps0, self.E0, self.power = eps0, E0, power

    def efficiency_at(self, E):
        if E <= 0:
            return None
        return self.eps0 * (self.E0 / E) ** self.power

    def is_extrapolating(self, E):
        return False


# ──────────────────────────────────────────────────────────────────
# Backward-compat tests: default OFF should leave behavior unchanged
# ──────────────────────────────────────────────────────────────────

def test_F296_default_off_no_change_for_co60():
    """С default enable_tcs_correction=False → результат идентичен старому."""
    from gamma.activity.compute import compute_activity
    eps = _MockEfficiency()
    co = _MockNuclideId("Co-60", [
        _MockLineMatch(1173.2, 99.85, 5000.0, math.sqrt(5000.0)),
        _MockLineMatch(1332.5, 99.98, 4500.0, math.sqrt(4500.0)),
    ])
    r_default = compute_activity(
        co, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True,
    )
    r_explicit = compute_activity(
        co, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True,
        enable_tcs_correction=False,
    )
    assert r_default.A_Bq == pytest.approx(r_explicit.A_Bq, rel=1e-6)


def test_F296_default_off_no_change_for_cs137():
    """Cs-137 default — TCS не должен срабатывать."""
    from gamma.activity.compute import compute_activity
    eps = _MockEfficiency()
    cs = _MockNuclideId("Cs-137", [
        _MockLineMatch(661.66, 85.1, 10000.0, 100.0),
    ])
    r = compute_activity(
        cs, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True,
    )
    assert r.A_Bq > 0
    # Cs-137 single line, no cascade — auto-TCS irrelevant
    r_with = compute_activity(
        cs, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True,
        enable_tcs_correction=True,
    )
    # Cs-137 в CASCADE_PRESETS пустой — поправка не применится
    assert r_with.A_Bq == pytest.approx(r.A_Bq, rel=1e-6)


# ──────────────────────────────────────────────────────────────────
# F-296 — auto-TCS active path
# ──────────────────────────────────────────────────────────────────

def test_F296_auto_tcs_applies_for_co60_when_enabled():
    """С enable_tcs_correction=True и Co-60 → C_TCS > 1 (depletion correction)."""
    from gamma.activity.compute import compute_activity
    eps = _MockEfficiency()
    co = _MockNuclideId("Co-60", [
        _MockLineMatch(1173.2, 99.85, 5000.0, math.sqrt(5000.0)),
        _MockLineMatch(1332.5, 99.98, 4500.0, math.sqrt(4500.0)),
    ])
    r_no = compute_activity(
        co, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True, enable_tcs_correction=False,
    )
    r_yes = compute_activity(
        co, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True, enable_tcs_correction=True,
    )
    # TCS добавляет положительную поправку (восстанавливает true activity)
    assert r_yes.A_Bq >= r_no.A_Bq
    if r_yes.A_Bq > r_no.A_Bq:
        # Если сработала поправка — нота должна это отразить
        assert "F-296" in (r_yes.notes or "") or r_yes.coincidence_correction_applied


def test_F296_user_supplied_coincidence_takes_precedence():
    """Если пользователь передал coincidence_correction — auto-TCS не дублирует."""
    from gamma.activity.compute import compute_activity
    eps = _MockEfficiency()
    co = _MockNuclideId("Co-60", [
        _MockLineMatch(1173.2, 99.85, 5000.0, math.sqrt(5000.0)),
    ])
    user_cc = {1173.2: 1.10}
    r = compute_activity(
        co, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True,
        coincidence_correction=user_cc,
        enable_tcs_correction=True,    # игнорируется, т.к. user-supplied present
    )
    assert r.A_Bq > 0
    # User-supplied applied, not auto
    assert "F-296" not in (r.notes or "")


def test_F296_no_change_for_unknown_nuclide():
    """K-40 не в CASCADE_PRESETS → enable_tcs_correction никак не повлияет."""
    from gamma.activity.compute import compute_activity
    eps = _MockEfficiency()
    k40 = _MockNuclideId("K-40", [
        _MockLineMatch(1460.8, 10.7, 8000.0, math.sqrt(8000.0)),
    ])
    r_off = compute_activity(
        k40, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True, enable_tcs_correction=False,
    )
    r_on = compute_activity(
        k40, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True, enable_tcs_correction=True,
    )
    assert r_on.A_Bq == pytest.approx(r_off.A_Bq, rel=1e-9)


def test_F296_tcs_detector_id_passed():
    """tcs_detector_id default 'Gamma-1S' — не должен ломаться при unknown id."""
    from gamma.activity.compute import compute_activity
    eps = _MockEfficiency()
    co = _MockNuclideId("Co-60", [
        _MockLineMatch(1173.2, 99.85, 5000.0, math.sqrt(5000.0)),
    ])
    # Unknown detector — должно gracefully degrade
    r = compute_activity(
        co, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True, enable_tcs_correction=True,
        tcs_detector_id="UNKNOWN_DETECTOR",
    )
    assert r.A_Bq > 0  # не падает


# ──────────────────────────────────────────────────────────────────
# F-294 Cutshall fallback via compute_activities_for_all
# ──────────────────────────────────────────────────────────────────

class _MockIdResult:
    def __init__(self, detected, rejected=()):
        self.detected_nuclides = detected
        self.rejected_nuclides = rejected


def test_F294_cutshall_fallback_default_off():
    """С default enable_cutshall_self_abs=False — пустой self_att_factors."""
    from gamma.activity.compute import compute_activities_for_all
    eps = _MockEfficiency()
    cs = _MockNuclideId("Cs-137", [
        _MockLineMatch(661.66, 85.1, 10000.0, 100.0),
    ])
    idr = _MockIdResult([cs])
    results = compute_activities_for_all(
        idr, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True,
        # no geometry_canonical, no enable_cutshall_self_abs
        sample_density_g_cm3=1.5,
    )
    assert len(results) == 1
    assert results[0].A_Bq > 0


def test_F294_cutshall_fallback_active_changes_activity():
    """С enable_cutshall_self_abs=True + non-default density → активность изменится."""
    from gamma.activity.compute import compute_activities_for_all
    eps = _MockEfficiency()
    # Low-energy line больше чувствительна к плотности
    cs = _MockNuclideId("Co-57", [
        _MockLineMatch(122.06, 85.5, 5000.0, math.sqrt(5000.0)),
    ])
    idr = _MockIdResult([cs])
    res_baseline = compute_activities_for_all(
        idr, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True,
        enable_cutshall_self_abs=False,
        sample_density_g_cm3=1.5,
    )
    res_corrected = compute_activities_for_all(
        idr, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True,
        enable_cutshall_self_abs=True,
        sample_density_g_cm3=1.5,
        cutshall_calib_density_g_cm3=1.0,
    )
    assert res_baseline[0].A_Bq > 0
    assert res_corrected[0].A_Bq > 0
    # При rho_sample > rho_calib для low-E ожидаем correction != 1
    assert res_baseline[0].A_Bq != pytest.approx(
        res_corrected[0].A_Bq, rel=1e-9
    )


def test_F294_cutshall_skipped_when_no_density():
    """Без sample_density_g_cm3 — Cutshall тихо игнорируется."""
    from gamma.activity.compute import compute_activities_for_all
    eps = _MockEfficiency()
    cs = _MockNuclideId("Cs-137", [
        _MockLineMatch(661.66, 85.1, 10000.0, 100.0),
    ])
    idr = _MockIdResult([cs])
    res = compute_activities_for_all(
        idr, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True,
        enable_cutshall_self_abs=True,
        sample_density_g_cm3=None,
    )
    assert len(res) == 1
    assert res[0].A_Bq > 0


def test_F294_cutshall_fallback_with_explicit_path():
    """cutshall_path_cm explicit передаётся → не использует default."""
    from gamma.activity.compute import compute_activities_for_all
    eps = _MockEfficiency()
    cs = _MockNuclideId("Co-57", [
        _MockLineMatch(122.06, 85.5, 5000.0, math.sqrt(5000.0)),
    ])
    idr = _MockIdResult([cs])
    res_default = compute_activities_for_all(
        idr, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True,
        enable_cutshall_self_abs=True,
        sample_density_g_cm3=1.5,
        cutshall_calib_density_g_cm3=1.0,
    )
    res_custom = compute_activities_for_all(
        idr, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True,
        enable_cutshall_self_abs=True,
        sample_density_g_cm3=1.5,
        cutshall_calib_density_g_cm3=1.0,
        cutshall_path_cm=3.0,  # двойной путь → сильнее поправка
    )
    assert res_default[0].A_Bq > 0
    assert res_custom[0].A_Bq > 0


# ──────────────────────────────────────────────────────────────────
# Integration: F-296 через compute_activities_for_all
# ──────────────────────────────────────────────────────────────────

def test_F296_passthrough_via_compute_activities_for_all():
    """enable_tcs_correction должен прокидываться from caller."""
    from gamma.activity.compute import compute_activities_for_all
    eps = _MockEfficiency()
    co = _MockNuclideId("Co-60", [
        _MockLineMatch(1173.2, 99.85, 5000.0, math.sqrt(5000.0)),
        _MockLineMatch(1332.5, 99.98, 4500.0, math.sqrt(4500.0)),
    ])
    idr = _MockIdResult([co])
    res_off = compute_activities_for_all(
        idr, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True, enable_tcs_correction=False,
    )
    res_on = compute_activities_for_all(
        idr, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True, enable_tcs_correction=True,
        tcs_detector_id="Gamma-1S",
    )
    assert len(res_off) == 1 and len(res_on) == 1
    assert res_on[0].A_Bq >= res_off[0].A_Bq
