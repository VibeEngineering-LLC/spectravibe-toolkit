# -*- coding: utf-8 -*-
"""v1.18.0 — Quasi-template (F-302..F-304)."""
from __future__ import annotations
import math, os, sys
import pytest

SCRIPTS = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts")
)
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)


# ──────────────────────────────────────────────────────────────────
# Common fixtures
# ──────────────────────────────────────────────────────────────────

def _channel_to_keV(ch):
    # Linear: 1 keV/ch с offset 0
    return float(ch) * 1.0


def _fwhm_at(E):
    # NaI ~7% @ 662 keV → const FWHM% scaling
    return 0.07 * math.sqrt(max(E, 1.0) * 662.0)


# ──────────────────────────────────────────────────────────────────
# F-302 — PPP sum-spectra builder
# ──────────────────────────────────────────────────────────────────

def test_F302_cs137_single_line_template():
    from gamma.activity.quasi_template_ppp import (
        NuclideLine, NuclideDef, build_nuclide_template,
    )
    cs137 = NuclideDef(
        nuclide_id="Cs-137",
        lines=[NuclideLine(E_keV=661.66, intensity=0.851, efficiency=0.05)],
    )
    t = build_nuclide_template(
        cs137, n_channels=1024,
        channel_to_keV=_channel_to_keV, fwhm_keV_at=_fwhm_at,
    )
    assert t.nuclide_id == "Cs-137"
    assert t.n_channels == 1024
    assert t.integral() == pytest.approx(0.851 * 0.05, rel=0.01)
    assert not t.has_continuum


def test_F302_co60_two_line_template():
    from gamma.activity.quasi_template_ppp import (
        NuclideLine, NuclideDef, build_nuclide_template,
    )
    co60 = NuclideDef(
        nuclide_id="Co-60",
        lines=[
            NuclideLine(E_keV=1173.2, intensity=0.999, efficiency=0.03),
            NuclideLine(E_keV=1332.5, intensity=0.9998, efficiency=0.028),
        ],
    )
    t = build_nuclide_template(
        co60, n_channels=2048,
        channel_to_keV=_channel_to_keV, fwhm_keV_at=_fwhm_at,
    )
    expected_area = 0.999 * 0.03 + 0.9998 * 0.028
    assert t.integral() == pytest.approx(expected_area, rel=0.01)


def test_F302_empty_nuclide_returns_zero_template():
    from gamma.activity.quasi_template_ppp import (
        NuclideDef, build_nuclide_template,
    )
    empty = NuclideDef(nuclide_id="EMPTY", lines=[])
    t = build_nuclide_template(
        empty, n_channels=512,
        channel_to_keV=_channel_to_keV, fwhm_keV_at=_fwhm_at,
    )
    assert t.integral() == 0.0
    assert "No lines" in (t.notes or "")


def test_F302_zero_n_channels_raises():
    from gamma.activity.quasi_template_ppp import (
        NuclideLine, NuclideDef, build_nuclide_template,
    )
    cs = NuclideDef("X", [NuclideLine(662.0, 1.0, 0.05)])
    with pytest.raises(ValueError):
        build_nuclide_template(
            cs, n_channels=0,
            channel_to_keV=_channel_to_keV, fwhm_keV_at=_fwhm_at,
        )


def test_F302_validate_detects_duplicate_ids():
    from gamma.activity.quasi_template_ppp import (
        PPPTemplate, validate_template_collection,
    )
    templates = [
        PPPTemplate("Cs-137", 512, [1.0] * 512),
        PPPTemplate("Cs-137", 512, [1.0] * 512),
    ]
    issues = validate_template_collection(templates)
    assert any("Duplicate" in i for i in issues)


def test_F302_validate_detects_channel_mismatch():
    from gamma.activity.quasi_template_ppp import (
        PPPTemplate, validate_template_collection,
    )
    templates = [
        PPPTemplate("Cs-137", 512, [1.0] * 512),
        PPPTemplate("K-40", 1024, [1.0] * 1024),
    ]
    issues = validate_template_collection(templates)
    assert any("n_channels" in i for i in issues)


def test_F302_validate_empty_collection():
    from gamma.activity.quasi_template_ppp import validate_template_collection
    assert validate_template_collection([]) == ["Empty template collection"]


def test_F302_template_centroid_alignment():
    """FEP centroid должен попасть в правильный канал."""
    from gamma.activity.quasi_template_ppp import (
        NuclideLine, NuclideDef, build_nuclide_template,
    )
    cs = NuclideDef(
        nuclide_id="Cs-137",
        lines=[NuclideLine(E_keV=662.0, intensity=1.0, efficiency=1.0)],
    )
    t = build_nuclide_template(
        cs, n_channels=1024,
        channel_to_keV=_channel_to_keV, fwhm_keV_at=_fwhm_at,
    )
    # Найти channel с макс counts
    max_ch = max(range(1024), key=lambda i: t.counts[i])
    assert abs(max_ch - 662) <= 1


def test_F302_with_continuum_function():
    from gamma.activity.quasi_template_ppp import (
        NuclideLine, NuclideDef, build_nuclide_template,
    )
    from gamma.activity.compton_continuum import make_continuum_func
    cs = NuclideDef(
        nuclide_id="Cs-137",
        lines=[NuclideLine(E_keV=662.0, intensity=1.0, efficiency=1.0)],
    )
    cont_func = make_continuum_func()
    t_no = build_nuclide_template(
        cs, 1024, _channel_to_keV, _fwhm_at,
    )
    t_with = build_nuclide_template(
        cs, 1024, _channel_to_keV, _fwhm_at,
        continuum_func=cont_func,
    )
    assert t_with.has_continuum
    assert t_with.integral() > t_no.integral()


# ──────────────────────────────────────────────────────────────────
# F-303 — Simultaneous fit solver
# ──────────────────────────────────────────────────────────────────

def test_F303_single_nuclide_recovery():
    """Известная активность Cs-137 → должна быть восстановлена."""
    from gamma.activity.quasi_template_ppp import (
        NuclideLine, NuclideDef, build_nuclide_template,
    )
    from gamma.activity.quasi_template_fit import solve_quasi_template_fit
    cs = NuclideDef(
        nuclide_id="Cs-137",
        lines=[NuclideLine(E_keV=662.0, intensity=0.851, efficiency=0.05)],
    )
    t = build_nuclide_template(cs, 1024, _channel_to_keV, _fwhm_at)
    A_true = 1000.0
    t_live = 600.0
    observed = [A_true * t_live * c for c in t.counts]

    result = solve_quasi_template_fit(
        observed=observed, templates=[t], live_time_s=t_live,
    )
    assert result.converged
    assert result.activities["Cs-137"] == pytest.approx(A_true, rel=0.01)


def test_F303_two_nuclides_separated():
    """Cs-137 + Co-60 с очень разными энергиями — should be separable."""
    from gamma.activity.quasi_template_ppp import (
        NuclideLine, NuclideDef, build_nuclide_template,
    )
    from gamma.activity.quasi_template_fit import solve_quasi_template_fit
    cs = NuclideDef("Cs-137", [NuclideLine(662.0, 0.851, 0.05)])
    co = NuclideDef("Co-60", [
        NuclideLine(1173.2, 0.999, 0.03),
        NuclideLine(1332.5, 0.9998, 0.028),
    ])
    t_cs = build_nuclide_template(cs, 2048, _channel_to_keV, _fwhm_at)
    t_co = build_nuclide_template(co, 2048, _channel_to_keV, _fwhm_at)
    A_cs, A_co, t_live = 500.0, 200.0, 1000.0
    observed = [
        A_cs * t_live * t_cs.counts[ch] + A_co * t_live * t_co.counts[ch]
        for ch in range(2048)
    ]
    result = solve_quasi_template_fit(observed, [t_cs, t_co], t_live)
    assert result.converged
    assert result.activities["Cs-137"] == pytest.approx(A_cs, rel=0.02)
    assert result.activities["Co-60"] == pytest.approx(A_co, rel=0.02)


def test_F303_chi2_acceptance():
    """Synthetic clean data → χ²_red должен быть very small."""
    from gamma.activity.quasi_template_ppp import (
        NuclideLine, NuclideDef, build_nuclide_template,
    )
    from gamma.activity.quasi_template_fit import solve_quasi_template_fit
    cs = NuclideDef("Cs-137", [NuclideLine(662.0, 0.851, 0.05)])
    t = build_nuclide_template(cs, 1024, _channel_to_keV, _fwhm_at)
    observed = [1000.0 * 600.0 * c for c in t.counts]
    result = solve_quasi_template_fit(observed, [t], 600.0)
    # Idealized data: χ²_red << 1
    assert result.is_accepted(chi2_red_max=1.5)


def test_F303_underdetermined_returns_failed():
    from gamma.activity.quasi_template_ppp import (
        NuclideLine, NuclideDef, build_nuclide_template,
    )
    from gamma.activity.quasi_template_fit import solve_quasi_template_fit
    # 1 channel в окне vs 2 nuclides → underdetermined
    cs = NuclideDef("Cs-137", [NuclideLine(662.0, 0.851, 0.05)])
    co = NuclideDef("Co-60", [NuclideLine(1173.2, 0.999, 0.03)])
    t1 = build_nuclide_template(cs, 1024, _channel_to_keV, _fwhm_at)
    t2 = build_nuclide_template(co, 1024, _channel_to_keV, _fwhm_at)
    obs = [0.0] * 1024
    result = solve_quasi_template_fit(
        observed=obs, templates=[t1, t2], live_time_s=600.0,
        energy_window=(500, 501),   # 1 channel < 2 nuclides
    )
    assert not result.converged
    assert any("Underdetermined" in n for n in result.notes)


def test_F303_zero_width_window_raises():
    from gamma.activity.quasi_template_ppp import (
        NuclideLine, NuclideDef, build_nuclide_template,
    )
    from gamma.activity.quasi_template_fit import solve_quasi_template_fit
    cs = NuclideDef("Cs-137", [NuclideLine(662.0, 0.851, 0.05)])
    t = build_nuclide_template(cs, 1024, _channel_to_keV, _fwhm_at)
    with pytest.raises(ValueError, match="zero width"):
        solve_quasi_template_fit(
            observed=[0.0] * 1024, templates=[t],
            live_time_s=600.0, energy_window=(500, 500),
        )


def test_F303_zero_live_time_raises():
    from gamma.activity.quasi_template_ppp import (
        NuclideLine, NuclideDef, build_nuclide_template,
    )
    from gamma.activity.quasi_template_fit import solve_quasi_template_fit
    cs = NuclideDef("Cs-137", [NuclideLine(662.0, 0.851, 0.05)])
    t = build_nuclide_template(cs, 1024, _channel_to_keV, _fwhm_at)
    with pytest.raises(ValueError):
        solve_quasi_template_fit([0.0] * 1024, [t], 0.0)


def test_F303_template_length_mismatch_raises():
    from gamma.activity.quasi_template_ppp import (
        PPPTemplate,
    )
    from gamma.activity.quasi_template_fit import solve_quasi_template_fit
    t = PPPTemplate("X", 512, [1.0] * 512)
    with pytest.raises(ValueError):
        solve_quasi_template_fit([0.0] * 1024, [t], 100.0)


def test_F303_negative_activity_noted():
    """Background-only spectrum + template → may yield A≈0 or slightly negative.
    Should report без crash."""
    from gamma.activity.quasi_template_ppp import (
        NuclideLine, NuclideDef, build_nuclide_template,
    )
    from gamma.activity.quasi_template_fit import solve_quasi_template_fit
    cs = NuclideDef("Cs-137", [NuclideLine(662.0, 0.851, 0.05)])
    t = build_nuclide_template(cs, 1024, _channel_to_keV, _fwhm_at)
    # Flat noise spectrum, no Cs-137 signal в районе 662
    obs = [5.0] * 1024
    result = solve_quasi_template_fit(obs, [t], 100.0)
    assert result.converged
    # Не падает; activity может быть около нуля или negative
    assert isinstance(result.activities["Cs-137"], float)


# ──────────────────────────────────────────────────────────────────
# F-304 — Compton continuum
# ──────────────────────────────────────────────────────────────────

def test_F304_compton_edge_cs137():
    """E=662 keV → E_edge ≈ 477.3 keV (canonical)."""
    from gamma.activity.compton_continuum import compton_edge_keV
    assert compton_edge_keV(662.0) == pytest.approx(477.3, abs=0.5)


def test_F304_compton_edge_co60_high():
    """E=1332 keV → E_edge ≈ 1118 keV."""
    from gamma.activity.compton_continuum import compton_edge_keV
    assert compton_edge_keV(1332.5) == pytest.approx(1118.0, abs=2.0)


def test_F304_compton_edge_zero_returns_zero():
    from gamma.activity.compton_continuum import compton_edge_keV
    assert compton_edge_keV(0.0) == 0.0
    assert compton_edge_keV(-100.0) == 0.0


def test_F304_backscatter_complement():
    """E_back = E - E_edge должно работать для arbitrary E."""
    from gamma.activity.compton_continuum import (
        compton_edge_keV, backscatter_peak_keV,
    )
    E = 662.0
    assert backscatter_peak_keV(E) == pytest.approx(E - compton_edge_keV(E))


def test_F304_continuum_zero_area_returns_zeros():
    from gamma.activity.compton_continuum import compton_continuum_for_line
    counts = compton_continuum_for_line(
        E_line_keV=662.0, compton_area=0.0,
        n_channels=512, channel_to_keV=_channel_to_keV,
    )
    assert all(c == 0.0 for c in counts)


def test_F304_continuum_preserves_area():
    """Integral counts должен ≈ compton_area (после re-normalization)."""
    from gamma.activity.compton_continuum import compton_continuum_for_line
    counts = compton_continuum_for_line(
        E_line_keV=662.0, compton_area=100.0,
        n_channels=1024, channel_to_keV=_channel_to_keV,
    )
    assert sum(counts) == pytest.approx(100.0, rel=0.05)


def test_F304_continuum_below_edge_dominant():
    """Большая часть площади ниже Compton-edge."""
    from gamma.activity.compton_continuum import (
        compton_continuum_for_line, compton_edge_keV,
    )
    E = 662.0
    edge = compton_edge_keV(E)
    counts = compton_continuum_for_line(
        E_line_keV=E, compton_area=1000.0,
        n_channels=1024, channel_to_keV=_channel_to_keV,
    )
    below = sum(counts[ch] for ch in range(int(edge)))
    above = sum(counts[ch] for ch in range(int(edge), 1024))
    assert below > above


def test_F304_make_continuum_func_factory():
    from gamma.activity.compton_continuum import make_continuum_func
    fn = make_continuum_func(fwhm_keV_at=_fwhm_at)
    counts = fn(662.0, 100.0, 1024, _channel_to_keV)
    assert len(counts) == 1024
    assert sum(counts) == pytest.approx(100.0, rel=0.05)


def test_F304_electron_rest_mass_constant():
    from gamma.activity.compton_continuum import ELECTRON_REST_MASS_KEV
    assert ELECTRON_REST_MASS_KEV == pytest.approx(511.0, abs=0.1)
