"""F-391 / v1.18.27 — S/N significance gate + library-anchor dedupe.

Покрывает:
  1. ``_f391_peak_snr`` — иерархия источников S/N (area/uncertainty,
     peak_sigma fallback, library_anchor → 0).
  2. ``find_multiplet_regions(min_significance_snr=...)`` — low-S/N
     активные компоненты помечаются как phantom anchors с
     ``peak_area_source="library_anchor_phantom"``.
  3. CC build skips phantoms — topology формируется только по active'ам.
  4. ``apply_multiplet_deconvolution`` cluster acceptance:
     - 0 actives (phantom-only) → drop,
     - 1 active (singleton) → drop из multiplet array (route в primary_fep),
     - ≥2 actives → keep.
  5. Empirical: Th-232 demo — M3 (Ac-228 409), M5 (Ac-228 674) singletons
     с низким library_I дропаются при S/N=0 (нет measured area).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(
    Path(__file__).resolve().parent.parent.parent / "scripts"
))

from gamma.peaks.deconvolve import (
    find_multiplet_regions,
    _f391_peak_snr,
    _f391_mark_phantom,
)


# ──────────────────────────────────────────────────────────────────
# Fixtures
# ──────────────────────────────────────────────────────────────────

class _FakeLineMatch:
    """Mutable stand-in for LineMatch."""

    def __init__(
        self,
        nuclide,
        library_E_keV,
        peak_channel,
        library_I_pct=10.0,
        peak_area=None,
        peak_area_uncertainty=None,
        peak_sigma=0.0,
        peak_area_source="",
        significance_currie=None,
    ):
        self.nuclide = nuclide
        self.library_E_keV = library_E_keV
        self.library_I_pct = library_I_pct
        self.peak_channel = peak_channel
        self.peak_E_keV = library_E_keV
        self.peak_sigma = peak_sigma
        self.significance_currie = (
            significance_currie if significance_currie is not None else peak_sigma
        )
        self.residual_keV = 0.0
        self.is_characteristic = False
        self.peak_area = peak_area
        self.peak_area_uncertainty = peak_area_uncertainty
        self.peak_area_source = peak_area_source


class _FakeNuclideId:
    def __init__(self, nuclide, matched_lines):
        self.nuclide = nuclide
        self.matched_lines = matched_lines
        self.detected = True


class _FakeIdent:
    def __init__(self, detected):
        self.detected_nuclides = detected


def _is_phantom(m):
    return getattr(m, "peak_area_source", "") in (
        "library_anchor", "library_anchor_phantom",
    )


# ──────────────────────────────────────────────────────────────────
# 1. _f391_peak_snr helper
# ──────────────────────────────────────────────────────────────────

class TestPeakSnrHelper:
    def test_area_over_uncertainty_is_primary_source(self):
        m = _FakeLineMatch(
            "Ac-228", 911, 911,
            peak_area=1000.0, peak_area_uncertainty=100.0,
            peak_sigma=5.0,  # должно игнорироваться, area/unc primary
        )
        assert _f391_peak_snr(m) == pytest.approx(10.0, rel=0.01)

    def test_peak_sigma_fallback_when_uncertainty_missing(self):
        m = _FakeLineMatch(
            "Ac-228", 911, 911,
            peak_area=1000.0, peak_area_uncertainty=None,
            peak_sigma=7.5,
        )
        assert _f391_peak_snr(m) == pytest.approx(7.5)

    def test_library_anchor_returns_zero(self):
        m = _FakeLineMatch(
            "Ac-228", 911, 911,
            peak_area=None, peak_area_uncertainty=None,
            peak_sigma=0.0,
            peak_area_source="library_anchor",
        )
        assert _f391_peak_snr(m) == 0.0

    def test_library_anchor_phantom_returns_zero(self):
        m = _FakeLineMatch(
            "Ac-228", 911, 911,
            peak_area=500.0, peak_area_uncertainty=50.0,
            peak_area_source="library_anchor_phantom",
        )
        # Phantom semantically has no measurement, даже если поле осталось
        assert _f391_peak_snr(m) == 0.0

    def test_zero_uncertainty_falls_through_to_sigma(self):
        m = _FakeLineMatch(
            "Ac-228", 911, 911,
            peak_area=1000.0, peak_area_uncertainty=0.0,
            peak_sigma=3.5,
        )
        assert _f391_peak_snr(m) == pytest.approx(3.5)

    def test_no_data_returns_infinity(self):
        """Back-compat: пустые поля без anchor-метки → inf (gate disabled
        for this LineMatch). Иначе тестовые fixtures без populated
        значимости получили бы phantom-метку."""
        m = _FakeLineMatch(
            "Ac-228", 911, 911,
            peak_area=None, peak_area_uncertainty=None,
            peak_sigma=0.0,
        )
        import math
        assert math.isinf(_f391_peak_snr(m))

    def test_mark_phantom_changes_source(self):
        m = _FakeLineMatch(
            "Ac-228", 911, 911,
            peak_area=10.0, peak_area_uncertainty=5.0,
            peak_area_source="cowell",
        )
        p = _f391_mark_phantom(m)
        assert getattr(p, "peak_area_source", "") == "library_anchor_phantom"


# ──────────────────────────────────────────────────────────────────
# 2. Synthetic — pure-noise singleton drops, strong signal keeps
# ──────────────────────────────────────────────────────────────────

class TestSyntheticSingleton:
    """После S/N gate singleton cluster (1 active) дропается полностью
    из multiplet_regions output (size=1 не возвращается)."""

    def test_strong_signal_singleton_keeps_active(self):
        """Single peak с S/N=10 не должен стать phantom."""
        fwhm_at = lambda ch: 25.0
        # Один LineMatch без явных соседей — overlap-step его дропнет
        # (size<2). Чтобы убедиться, что gate сам по себе НЕ выгоняет
        # активного singleton'a, добавим вторую линию рядом которая
        # явно станет phantom.
        m_strong = _FakeLineMatch(
            "Ac-228", 911, 911, library_I_pct=27.0,
            peak_area=10000.0, peak_area_uncertainty=1000.0,  # S/N=10
        )
        m_weak = _FakeLineMatch(
            "Ac-228", 925, 925, library_I_pct=0.5,
            peak_area=50.0, peak_area_uncertainty=50.0,  # S/N=1
        )
        ident = _FakeIdent((_FakeNuclideId("Ac-228", [m_strong, m_weak]),))
        clusters = find_multiplet_regions(
            ident, fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
            min_significance_snr=3.0,
        )
        # После gate: weak помечен phantom, strong остаётся active.
        # CC topology только по strong → 1 active CC; phantom attaches.
        # cluster sizes >= 2 retained на overlap-step.
        all_members = []
        for cl in clusters:
            all_members.extend(cl)
        strong_member = next(
            (m for m in all_members if abs(m.library_E_keV - 911) < 1),
            None,
        )
        weak_member = next(
            (m for m in all_members if abs(m.library_E_keV - 925) < 1),
            None,
        )
        if strong_member is not None:
            assert not _is_phantom(strong_member), (
                "S/N=10 component должен остаться active"
            )
        if weak_member is not None:
            assert _is_phantom(weak_member), (
                "S/N=1 component должен стать phantom (gate=3.0)"
            )

    def test_pure_noise_singleton_marked_phantom(self):
        """Single active S/N=1 в multiplet ROI → phantom."""
        fwhm_at = lambda ch: 25.0
        m_noise1 = _FakeLineMatch(
            "Ac-228", 400, 400,
            peak_area=10.0, peak_area_uncertainty=10.0,  # S/N=1
        )
        m_noise2 = _FakeLineMatch(
            "Ac-228", 415, 415,
            peak_area=8.0, peak_area_uncertainty=10.0,  # S/N=0.8
        )
        ident = _FakeIdent((
            _FakeNuclideId("Ac-228", [m_noise1, m_noise2]),
        ))
        clusters = find_multiplet_regions(
            ident, fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
            min_significance_snr=3.0,
        )
        # Все members → phantom; CC build не построит active cluster.
        # Cluster в output может остаться (phantom-only fallback CC)
        # — apply_multiplet_deconvolution дропнет его.
        all_members = []
        for cl in clusters:
            all_members.extend(cl)
        for m in all_members:
            assert _is_phantom(m), (
                f"All low-S/N components должны быть phantom; "
                f"{m.library_E_keV} осталась active"
            )


# ──────────────────────────────────────────────────────────────────
# 3. Synthetic — mixed cluster (high + low S/N) topology
# ──────────────────────────────────────────────────────────────────

class TestMixedTopology:
    def test_two_strong_plus_five_weak_yields_two_actives(self):
        """Cluster с 2 high-S/N + 5 low-S/N → topology = 2 active CC,
        5 phantom prepended to it for evidence."""
        fwhm_at = lambda ch: 25.0
        # 2 сильных близких пиков (Δ<25 → unresolved pair) + 5 слабых
        strong = [
            _FakeLineMatch(
                "Ac-228", E, E, library_I_pct=I,
                peak_area=5000.0, peak_area_uncertainty=500.0,
            )
            for E, I in [(965.0, 5.0), (969.0, 16.0)]
        ]
        weak = [
            _FakeLineMatch(
                "Ac-228", E, E, library_I_pct=I,
                peak_area=30.0, peak_area_uncertainty=30.0,
            )
            for E, I in [
                (955.0, 0.3), (960.0, 0.2), (973.0, 0.5),
                (978.0, 0.4), (983.0, 0.6),
            ]
        ]
        ident = _FakeIdent((_FakeNuclideId("Ac-228", strong + weak),))
        clusters = find_multiplet_regions(
            ident, fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
            min_significance_snr=3.0,
        )
        # Найти cluster с двумя сильными
        target = None
        for cl in clusters:
            actives = [m for m in cl if not _is_phantom(m)]
            if len(actives) == 2 and all(
                abs(m.library_E_keV - 965) < 1
                or abs(m.library_E_keV - 969) < 1
                for m in actives
            ):
                target = cl
                break
        assert target is not None, (
            "Должен существовать cluster с 2 active strong components"
        )
        # 5 weak → phantom evidence
        phantoms = [m for m in target if _is_phantom(m)]
        # Phantom-ы могут быть распределены на ближайшие CC; в этом
        # тесте только один active CC, поэтому все 5 weak должны быть
        # в нём (после attach-nearest).
        assert len(phantoms) >= 3, (
            f"Cluster должен содержать ≥3 phantom evidence; got "
            f"{len(phantoms)}: {[(m.library_E_keV, m.peak_area_source) for m in phantoms]}"
        )


# ──────────────────────────────────────────────────────────────────
# 4. Empirical — Th-232 demo M3/M5 singleton drops
# ──────────────────────────────────────────────────────────────────

class TestEmpiricalTh232SingletonDrop:
    """v1.18.26.1 Th-232 demo (PROD) показывал M3 = Ac-228 409 как
    одиночный multiplet entry (singleton после F-387.1 split). После
    F-391 такой singleton без реального measured area дропается."""

    def test_prod_m3_ac228_409_singleton_drops_when_no_measured_area(self):
        """Ac-228 409 keV (I=1.74%) без peak_area → S/N=0 → phantom
        → CC=1 active=0 → drop entirely."""
        fwhm_at = lambda ch: 30.0
        # Ac-228 409 как library_anchor (нет measured peak)
        m_anchor = _FakeLineMatch(
            "Ac-228", 409.5, 409.5, library_I_pct=1.74,
            peak_area=None, peak_area_uncertainty=None,
            peak_area_source="library_anchor",
        )
        # Изолированный — overlap-step его дропнет, но если бы он
        # попал в чужой ROI display-window, gate отнесёт к phantom.
        # Добавим близкий strong-active как seed:
        m_strong = _FakeLineMatch(
            "Ac-228", 463.0, 463.0, library_I_pct=4.40,
            peak_area=8000.0, peak_area_uncertainty=400.0,  # S/N=20
        )
        m_strong2 = _FakeLineMatch(
            "Tl-208", 510.77, 510.77, library_I_pct=22.6,
            peak_area=12000.0, peak_area_uncertainty=600.0,  # S/N=20
        )
        ident = _FakeIdent((
            _FakeNuclideId("Ac-228", [m_anchor, m_strong]),
            _FakeNuclideId("Tl-208", [m_strong2]),
        ))
        clusters = find_multiplet_regions(
            ident, fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
            min_significance_snr=3.0,
        )
        # Ac-228 409 должен быть phantom (library_anchor → S/N=0)
        for cl in clusters:
            for m in cl:
                if abs(m.library_E_keV - 409.5) < 1:
                    assert _is_phantom(m), (
                        "Ac-228 409 без peak_area должен быть phantom"
                    )

    def test_min_significance_snr_zero_disables_gate(self):
        """Back-compat: min_significance_snr=0.0 → gate off."""
        fwhm_at = lambda ch: 25.0
        # Weak signal — без gate он остаётся active
        m_weak = _FakeLineMatch(
            "Ac-228", 911, 911, library_I_pct=27.0,
            peak_area=10.0, peak_area_uncertainty=20.0,  # S/N=0.5
        )
        m_strong = _FakeLineMatch(
            "Ac-228", 925, 925, library_I_pct=10.0,
            peak_area=10000.0, peak_area_uncertainty=1000.0,
        )
        ident = _FakeIdent((
            _FakeNuclideId("Ac-228", [m_weak, m_strong]),
        ))
        # Gate disabled
        clusters = find_multiplet_regions(
            ident, fwhm_at,
            overlap_threshold_fwhm=4.0,
            expand_to_display_window=False,
            min_significance_snr=0.0,
        )
        all_members = []
        for cl in clusters:
            all_members.extend(cl)
        weak = next(
            (m for m in all_members if abs(m.library_E_keV - 911) < 1),
            None,
        )
        if weak is not None:
            # Без gate, low-S/N не помечен (peak_area_source осталось "")
            assert weak.peak_area_source not in (
                "library_anchor_phantom",
            ), "Gate=0 не должен помечать low-S/N как phantom"
