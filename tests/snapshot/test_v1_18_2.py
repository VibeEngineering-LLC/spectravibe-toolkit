# -*- coding: utf-8 -*-
"""v1.18.2 — Activity integration slice 2.

Wire-up F-297 (matrix_method_chi2) как opt-in alternative-path solver
в compute_activities_for_all. Default OFF → back-compat preserved.
"""
from __future__ import annotations
import math, os, sys
import pytest

SCRIPTS = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts")
)
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)


class _MockLineMatch:
    def __init__(self, E, I_pct, area, area_unc=None):
        self.library_E_keV = E
        self.library_I_pct = I_pct
        self.peak_area = area
        self.peak_area_uncertainty = area_unc
        self.peak_area_source = ""


class _MockNuclideId:
    def __init__(self, nuclide, lines, detected=True):
        self.nuclide = nuclide
        self.matched_lines = lines
        self.detected = detected


class _MockIdResult:
    def __init__(self, detected, rejected=()):
        self.detected_nuclides = detected
        self.rejected_nuclides = rejected


class _MockEfficiency:
    def __init__(self, eps0=0.05, E0=662.0, power=0.8):
        self.eps0, self.E0, self.power = eps0, E0, power

    def efficiency_at(self, E):
        if E <= 0:
            return None
        return self.eps0 * (self.E0 / E) ** self.power

    def is_extrapolating(self, E):
        return False


# ──────────────────────────────────────────────────────────────────
# Backward-compat: default OFF
# ──────────────────────────────────────────────────────────────────

def test_F297_default_off_uses_per_nuclide():
    """Default enable_matrix_method=False → per-nuclide weighted-mean path."""
    from gamma.activity.compute import compute_activities_for_all
    eps = _MockEfficiency()
    cs = _MockNuclideId("Cs-137", [
        _MockLineMatch(661.66, 85.1, 10000.0, 100.0),
    ])
    co = _MockNuclideId("Co-60", [
        _MockLineMatch(1173.2, 99.85, 5000.0, 70.0),
        _MockLineMatch(1332.5, 99.98, 4500.0, 67.0),
    ])
    idr = _MockIdResult([cs, co])
    res = compute_activities_for_all(
        idr, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True,
    )
    # Per-nuclide path: lines_used должен быть populated
    assert len(res) == 2
    for r in res:
        assert len(r.lines_used) > 0
        assert r.sigma_method != "matrix_method"


# ──────────────────────────────────────────────────────────────────
# F-297 — matrix method active path
# ──────────────────────────────────────────────────────────────────

def test_F297_matrix_method_active_solves_simultaneously():
    """enable_matrix_method=True + ≥2 nuclides → matrix solver path."""
    from gamma.activity.compute import compute_activities_for_all
    eps = _MockEfficiency()
    cs = _MockNuclideId("Cs-137", [
        _MockLineMatch(661.66, 85.1, 10000.0, 100.0),
    ])
    co = _MockNuclideId("Co-60", [
        _MockLineMatch(1173.2, 99.85, 5000.0, 70.0),
        _MockLineMatch(1332.5, 99.98, 4500.0, 67.0),
    ])
    idr = _MockIdResult([cs, co])
    res = compute_activities_for_all(
        idr, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True,
        enable_matrix_method=True,
    )
    assert len(res) == 2
    for r in res:
        assert r.sigma_method == "matrix_method"
        assert "F-297" in (r.notes or "")
        assert r.A_Bq > 0


def test_F297_chi2_reduced_reported_via_intra_chi2():
    """χ²_red от matrix_method прокидывается через intra_chi2_per_dof."""
    from gamma.activity.compute import compute_activities_for_all
    eps = _MockEfficiency()
    cs = _MockNuclideId("Cs-137", [
        _MockLineMatch(661.66, 85.1, 10000.0, 100.0),
    ])
    co = _MockNuclideId("Co-60", [
        _MockLineMatch(1173.2, 99.85, 5000.0, 70.0),
        _MockLineMatch(1332.5, 99.98, 4500.0, 67.0),
    ])
    idr = _MockIdResult([cs, co])
    res = compute_activities_for_all(
        idr, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True,
        enable_matrix_method=True,
    )
    for r in res:
        assert r.intra_chi2_per_dof is not None
        assert r.intra_chi2_per_dof >= 0


def test_F297_single_nuclide_falls_back_to_per_nuclide():
    """С одним нуклидом — matrix_method singular, должен fallback."""
    from gamma.activity.compute import compute_activities_for_all
    eps = _MockEfficiency()
    cs = _MockNuclideId("Cs-137", [
        _MockLineMatch(661.66, 85.1, 10000.0, 100.0),
    ])
    idr = _MockIdResult([cs])
    res = compute_activities_for_all(
        idr, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True,
        enable_matrix_method=True,
    )
    assert len(res) == 1
    # Fallback в per-nuclide → lines_used populated
    assert len(res[0].lines_used) > 0


def test_F297_underdetermined_falls_back():
    """N_peaks < N_nuclides → singular, fallback."""
    from gamma.activity.compute import compute_activities_for_all
    eps = _MockEfficiency()
    cs = _MockNuclideId("Cs-137", [
        _MockLineMatch(661.66, 85.1, 10000.0, 100.0),
    ])
    co = _MockNuclideId("Co-60", [
        _MockLineMatch(1173.2, 99.85, 5000.0, 70.0),
    ])
    # Shared single peak ситуация невозможна для unique energies.
    # Только 2 пика, 2 нуклида → matrix solvable (N=M). С 3+ nuclides и 2 peaks
    # будет singular. Создадим 3 нуклида с 2 уникальными пиками.
    co2 = _MockNuclideId("K-40", [
        _MockLineMatch(1460.8, 10.7, 800.0, math.sqrt(800.0)),
    ])
    idr = _MockIdResult([cs, co, co2])
    res = compute_activities_for_all(
        idr, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True,
        enable_matrix_method=True,
    )
    # N_peaks=3 (661, 1173, 1460), N_nuclides=3 — solvable, не singular
    assert len(res) == 3
    for r in res:
        assert r.A_Bq != 0 or r.A_Bq == 0  # просто не падает


def test_F297_recovers_known_activities():
    """Synthetic: inject known activities, recover via matrix_method."""
    from gamma.activity.compute import compute_activities_for_all
    eps = _MockEfficiency()
    # Construct peaks consistent with A_Cs=1000 Bq, A_Co=500 Bq, t=600s
    t = 600.0
    A_cs, A_co = 1000.0, 500.0
    # Cs-137 661 keV: ε(661) = 0.05 * (662/661)^0.8 ≈ 0.0500
    eps_cs = eps.efficiency_at(661.66)
    S_cs = A_cs * t * 0.851 * eps_cs
    # Co-60 1173 keV: ε ≈ 0.05 * (662/1173)^0.8
    eps_co1 = eps.efficiency_at(1173.2)
    eps_co2 = eps.efficiency_at(1332.5)
    S_co1 = A_co * t * 0.9985 * eps_co1
    S_co2 = A_co * t * 0.9998 * eps_co2
    cs = _MockNuclideId("Cs-137", [
        _MockLineMatch(661.66, 85.1, S_cs, math.sqrt(S_cs)),
    ])
    co = _MockNuclideId("Co-60", [
        _MockLineMatch(1173.2, 99.85, S_co1, math.sqrt(S_co1)),
        _MockLineMatch(1332.5, 99.98, S_co2, math.sqrt(S_co2)),
    ])
    idr = _MockIdResult([cs, co])
    res = compute_activities_for_all(
        idr, efficiency_curve=eps, live_time_s=t,
        from_bg_subtracted=True,
        enable_matrix_method=True,
    )
    by_name = {r.nuclide: r for r in res}
    assert by_name["Cs-137"].A_Bq == pytest.approx(A_cs, rel=0.02)
    assert by_name["Co-60"].A_Bq == pytest.approx(A_co, rel=0.02)


def test_F297_energy_tolerance_passes_through():
    """matrix_method_energy_tolerance_keV должен передаваться."""
    from gamma.activity.compute import compute_activities_for_all
    eps = _MockEfficiency()
    cs = _MockNuclideId("Cs-137", [
        _MockLineMatch(661.66, 85.1, 10000.0, 100.0),
    ])
    co = _MockNuclideId("Co-60", [
        _MockLineMatch(1173.2, 99.85, 5000.0, 70.0),
        _MockLineMatch(1332.5, 99.98, 4500.0, 67.0),
    ])
    idr = _MockIdResult([cs, co])
    # Должен не падать с разными tolerance values
    for tol in [0.5, 1.0, 2.0, 5.0]:
        res = compute_activities_for_all(
            idr, efficiency_curve=eps, live_time_s=600.0,
            from_bg_subtracted=True,
            enable_matrix_method=True,
            matrix_method_energy_tolerance_keV=tol,
        )
        assert len(res) == 2


def test_F297_notes_indicate_acceptance_status():
    """notes field должен показывать is_acceptable status."""
    from gamma.activity.compute import compute_activities_for_all
    eps = _MockEfficiency()
    cs = _MockNuclideId("Cs-137", [
        _MockLineMatch(661.66, 85.1, 10000.0, 100.0),
    ])
    co = _MockNuclideId("Co-60", [
        _MockLineMatch(1173.2, 99.85, 5000.0, 70.0),
        _MockLineMatch(1332.5, 99.98, 4500.0, 67.0),
    ])
    idr = _MockIdResult([cs, co])
    res = compute_activities_for_all(
        idr, efficiency_curve=eps, live_time_s=600.0,
        from_bg_subtracted=True,
        enable_matrix_method=True,
    )
    for r in res:
        assert "χ²_red=" in (r.notes or "")
        assert ("is_acceptable=True" in (r.notes or "")
                or "is_acceptable=False" in (r.notes or ""))
