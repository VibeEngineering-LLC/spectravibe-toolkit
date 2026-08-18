# -*- coding: utf-8 -*-
"""v1.18.5 — Methodology cleanup.

- T-020 verification: seven_line_check.py (F-81) exists and produces correct shape
- T-014 verification: OISN_16_COMPOSITION sums to 1.0, has expected elements
- F-148 (partial): Bi-214 + Bi-212 cascades added to CASCADE_PRESETS
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
# T-020 — Phase D 7-line ЕРН gate verification
# ──────────────────────────────────────────────────────────────────

def test_T020_seven_lines_constant_has_seven_entries():
    from gamma.calibration.seven_line_check import SEVEN_LINES
    assert len(SEVEN_LINES) == 7


def test_T020_seven_lines_canonical_energies():
    """Per LSRM methodology — точные энергии ЕРН реперов."""
    from gamma.calibration.seven_line_check import SEVEN_LINES
    energies = sorted(e for e, _, _ in SEVEN_LINES)
    expected = sorted([
        240.00,    # Pb-212/Pb-214 superposition
        351.93,    # Pb-214
        511.00,    # Tl-208 + annihilation
        1120.29,   # Bi-214
        1460.82,   # K-40
        1764.49,   # Bi-214
        2614.51,   # Tl-208
    ])
    for a, b in zip(energies, expected):
        assert abs(a - b) < 0.01


def test_T020_run_seven_line_check_returns_result():
    """Smoke: function callable + returns proper dataclass with fields."""
    from gamma.calibration.seven_line_check import (
        SevenLineCheck, run_seven_line_check,
    )

    # Mock peak (channel-based) + spec with channel_to_energy
    class _MockPeak:
        def __init__(self, channel, sig=10.0):
            self.channel = channel
            self.significance = sig

    class _MockSpec:
        @staticmethod
        def channel_to_energy(ch):
            return float(ch) * 1.0   # 1 keV/channel

    # 5 of 7 peaks present (channel ≈ keV since gain=1)
    peaks = [
        _MockPeak(240),
        _MockPeak(352),
        _MockPeak(1461),
        _MockPeak(1764),
        _MockPeak(2614),
    ]
    result = run_seven_line_check(
        peaks, _MockSpec(),
        fwhm_provider_keV=lambda E: 0.07 * math.sqrt(max(E, 1.0) * 662.0),
    )
    assert isinstance(result, SevenLineCheck)
    assert result.lines_total == 7
    assert result.lines_present == 5
    assert result.is_reliable     # ≥4
    assert result.quality in ("ok", "drift", "broken", "n/a")


def test_T020_no_peaks_yields_broken_or_na():
    """Без peaks → quality broken/n/a, не падает."""
    from gamma.calibration.seven_line_check import run_seven_line_check

    class _MockSpec:
        @staticmethod
        def channel_to_energy(ch):
            return float(ch) * 1.0

    result = run_seven_line_check(
        [], _MockSpec(),
        fwhm_provider_keV=lambda E: 50.0,
    )
    assert result.lines_present == 0
    assert not result.is_reliable


# ──────────────────────────────────────────────────────────────────
# T-014 — ОИСН-16 composition verification
# ──────────────────────────────────────────────────────────────────

def test_T014_oisn16_composition_sums_to_one():
    """Mass fractions ОИСН-16 должны быть сумма 1.0 ± 1e-3."""
    from gamma.physics.self_attenuation import OISN_16_COMPOSITION
    total = sum(OISN_16_COMPOSITION.values())
    assert total == pytest.approx(1.0, abs=1e-3)


def test_T014_oisn16_has_required_elements():
    """Документация: H, C, N, O, Fe должны быть в составе."""
    from gamma.physics.self_attenuation import OISN_16_COMPOSITION
    for elem in ("H", "C", "N", "O", "Fe"):
        assert elem in OISN_16_COMPOSITION
        assert OISN_16_COMPOSITION[elem] > 0


def test_T014_oisn16_fe_dominant():
    """Fe — главный compound в ОИСН-16 (~71.4%)."""
    from gamma.physics.self_attenuation import OISN_16_COMPOSITION
    fractions = OISN_16_COMPOSITION
    fe = fractions["Fe"]
    other = sum(v for k, v in fractions.items() if k != "Fe")
    assert fe > other


def test_T014_ref_geometry_has_marinelli():
    """REF_GEOMETRY должен содержать Marinelli."""
    from gamma.physics.self_attenuation import REF_GEOMETRY
    # Допустимые ключи: проверяем что Marinelli присутствует под одним
    # из канонических имён.
    has_marinelli = any(
        "marinelli" in k.lower() or "маринелли" in k.lower()
        for k in REF_GEOMETRY.keys()
    )
    assert has_marinelli


# ──────────────────────────────────────────────────────────────────
# F-148 — Bi-214 / Bi-212 cascade addition
# ──────────────────────────────────────────────────────────────────

def test_F148_bi214_in_cascade_presets():
    from gamma.activity.tcs_close_geometry import CASCADE_PRESETS
    assert "Bi-214" in CASCADE_PRESETS
    pairs = CASCADE_PRESETS["Bi-214"]
    assert len(pairs) >= 4
    # 609.31 / 1120.29 / 1764.49 — главные ЕРН-линии Bi-214
    energies = {round(p.E_i_keV, 0) for p in pairs}
    assert 609.0 in energies
    assert 1120.0 in energies
    assert 1764.0 in energies


def test_F148_bi212_in_cascade_presets():
    from gamma.activity.tcs_close_geometry import CASCADE_PRESETS
    assert "Bi-212" in CASCADE_PRESETS
    pairs = CASCADE_PRESETS["Bi-212"]
    assert len(pairs) >= 2


def test_F148_bi214_pairs_have_positive_probability():
    from gamma.activity.tcs_close_geometry import BI214_PAIRS
    for p in BI214_PAIRS:
        assert 0.0 < p.coincidence_probability < 1.0
        assert p.E_i_keV > 0
        assert p.E_j_keV > 0


def test_F148_bi214_compute_tcs_smoke():
    """compute_tcs_correction для Bi-214 не падает с simple eps_T."""
    from gamma.activity.tcs_close_geometry import (
        CASCADE_PRESETS, compute_tcs_correction,
    )
    # Близкая geometry для Marinelli → eps_T ~ 0.05 при close geometry
    def _eps_T(E):
        return 0.05 * (662.0 / max(E, 1.0)) ** 0.4
    res = compute_tcs_correction(
        E_i_keV=609.31,
        nuclide_pairs=CASCADE_PRESETS["Bi-214"],
        total_efficiency_func=_eps_T,
    )
    assert res.correction_factor >= 1.0
    assert res.n_pairs_used >= 1


def test_F148_existing_presets_preserved():
    """Старые presets (Co-60, Eu-152, Cs-137, Ba-133) должны остаться."""
    from gamma.activity.tcs_close_geometry import CASCADE_PRESETS
    for nuc in ("Co-60", "Eu-152", "Cs-137", "Ba-133"):
        assert nuc in CASCADE_PRESETS
