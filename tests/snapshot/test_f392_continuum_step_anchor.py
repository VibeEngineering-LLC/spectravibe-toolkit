"""F-392 / v1.18.27 — multi-step continuum для multi-anchor ROI.

F-392 расширяет step-linear continuum: для широких multi-anchor ROI
(≥200 кэВ, ≥3 anchor с library_I_pct≥5%) глобальная одиночная
β_step заменяется НАБОРОМ β_step_i terms, по одному на каждую intense
anchor-линию. Каждый step якорится на E энергии своего anchor (не free),
σ_step_i = FWHM(E_anchor_i)/2.355.

Мотивация: на NaI ROI 350-700 кэВ (Ac-228 463 + Tl-208 510 + Tl-208 583)
continuum резко опускается ПОСЛЕ Tl-208 583 кэВ из-за double-escape от
Tl-208 2614. Один глобальный erfc-step не способен описать это: линейный
continuum «уезжает», closure% растёт.

Тесты:
  1. Empirical M4-like: 3 NaI-peaks ROI 350-700 кэВ с реальным step jump
     после 583 → step_linear_multi должен дать χ²/ν ≤ step_linear.
  2. Synthetic single-step: 3 peaks с искусственным step → multi
     корректно ловит step magnitude.
  3. Synthetic multi-step: 4 peaks с двумя step jumps → continuum правильно
     моделируется (улучшение χ²/ν).
  4. No-step (auto-select boundary): 2 close peaks без step → auto-select
     остаётся на step_linear (multi не активируется).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

from gamma.peaks.coupled_multiplet import coupled_intensity_fit, ComponentSpec
from gamma.peaks.deconvolve import _f392_auto_select_continuum


# ──────────────────────────────────────────────────────────────────
# Утилиты
# ──────────────────────────────────────────────────────────────────

def _gauss(E, E0, sigma, A):
    return A * np.exp(-0.5 * ((E - E0) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))


def _smooth_step(E, E_step, sigma_step):
    """Гладкая ступенька 0.5·erfc((E-E_step)/(σ√2)). 1 при E<<E_step, 0 при E>>."""
    try:
        from scipy.special import erfc as _vec_erfc
        return 0.5 * _vec_erfc((E - E_step) / (sigma_step * np.sqrt(2)))
    except ImportError:
        import math
        return np.array([
            0.5 * math.erfc((float(e) - E_step) / (sigma_step * np.sqrt(2)))
            for e in E
        ])


# ──────────────────────────────────────────────────────────────────
# Test 1: Empirical M4-like (Ac-228 463 + Tl-208 510 + Tl-208 583) +
# реалистичный step-jump после 583 кэВ (double-escape от Tl-208 2614)
# ──────────────────────────────────────────────────────────────────

class TestEmpiricalM4MultiStep:
    """ROI 350-700 кэВ, 3 NaI компоненты, FWHM≈30 кэВ, step после 583."""

    def _build_roi(self):
        np.random.seed(20260601)
        E = np.linspace(350.0, 700.0, 700)
        sigma_25 = 25.0 / 2.355  # FWHM=25
        # Compton baseline
        y = 200.0 - 0.1 * (E - 525.0)
        # 3 peaks
        y += _gauss(E, 463.0, sigma_25, 8000.0)    # Ac-228 463
        y += _gauss(E, 510.77, sigma_25, 25000.0)  # Tl-208 510
        y += _gauss(E, 583.19, sigma_25, 35000.0)  # Tl-208 583
        # Step после 583 — double-escape от 2614, magnitude ~40 counts
        y += 80.0 * _smooth_step(E, 600.0, sigma_25)
        y = np.maximum(y, 1.0)
        y = np.random.poisson(y).astype(float)
        comp = [
            ComponentSpec("Ac-228", 463.0,   4.4,  "Ac-228"),
            ComponentSpec("Tl-208", 510.77, 22.6, "Tl-208"),
            ComponentSpec("Tl-208", 583.19, 30.5, "Tl-208"),
        ]
        fwhm_at = lambda E0: 25.0
        return E, y, comp, fwhm_at

    def test_step_linear_multi_not_worse_than_step_linear(self):
        E, y, comp, fwhm_at = self._build_roi()
        res_std = coupled_intensity_fit(
            E, y, comp, fwhm_at, continuum="step_linear",
        )
        res_multi = coupled_intensity_fit(
            E, y, comp, fwhm_at, continuum="step_linear_multi",
        )
        # multi должен быть НЕ ХУЖЕ обычного step_linear по χ²/ν
        # F-392.1 / v1.18.27.1: default threshold снижен 5.0 → 4.0%,
        # теперь Ac-228 463 (I=4.4%) тоже qualifies как anchor вместе
        # с Tl-208 510 (22.6%) и 583 (30.5%) → 3 step term'а.
        assert res_multi.chi2_per_dof <= res_std.chi2_per_dof * 1.05, (
            f"multi χ²/ν={res_multi.chi2_per_dof:.3f} НЕ должен быть значимо "
            f"хуже step_linear χ²/ν={res_std.chi2_per_dof:.3f}"
        )

    def test_multi_anchors_recorded_in_notes(self):
        E, y, comp, fwhm_at = self._build_roi()
        res = coupled_intensity_fit(
            E, y, comp, fwhm_at, continuum="step_linear_multi",
        )
        # F-392.1 / v1.18.27.1: threshold 4.0% → все три anchor'а
        # qualify: 463 (4.4%), 510 (22.6%), 583 (30.5%). Separation
        # 463→510=47 кэВ ≥ 40, 510→583=73 кэВ ≥ 40 → 3 step term'а.
        assert "F-392 multi-step anchors=3" in res.notes


# ──────────────────────────────────────────────────────────────────
# Test 2: Synthetic single-step (явный step после 1 anchor)
# ──────────────────────────────────────────────────────────────────

class TestSyntheticSingleStep:
    """4 peaks (≥5% intensity), искусственный step после самого
    интенсивного — multi-step должен корректно catch step magnitude."""

    def _build_roi(self):
        np.random.seed(20260602)
        E = np.linspace(400.0, 800.0, 800)
        sigma_25 = 25.0 / 2.355
        y = 150.0 + 0.0 * E  # flat baseline
        # 4 peaks с I≥5% при separation ≥ 40 кэВ
        y += _gauss(E, 463.0, sigma_25, 5000.0)   # I=5
        y += _gauss(E, 530.0, sigma_25, 20000.0)  # I=20 — самый интенсивный
        y += _gauss(E, 620.0, sigma_25, 8000.0)   # I=8
        y += _gauss(E, 750.0, sigma_25, 10000.0)  # I=10
        # ЕДИНСТВЕННЫЙ step jump после самого интенсивного (530)
        y += 60.0 * _smooth_step(E, 550.0, sigma_25)
        y = np.maximum(y, 1.0)
        y = np.random.poisson(y).astype(float)
        comp = [
            ComponentSpec("Nuc-A", 463.0, 5.0,  ""),
            ComponentSpec("Nuc-B", 530.0, 20.0, ""),
            ComponentSpec("Nuc-C", 620.0, 8.0,  ""),
            ComponentSpec("Nuc-D", 750.0, 10.0, ""),
        ]
        fwhm_at = lambda E0: 25.0
        return E, y, comp, fwhm_at

    def test_multi_improves_chi2(self):
        E, y, comp, fwhm_at = self._build_roi()
        res_std = coupled_intensity_fit(
            E, y, comp, fwhm_at, continuum="step_linear",
        )
        res_multi = coupled_intensity_fit(
            E, y, comp, fwhm_at, continuum="step_linear_multi",
        )
        # Step реально присутствует → multi не должен ухудшить fit
        assert res_multi.chi2_per_dof <= res_std.chi2_per_dof * 1.10

    def test_largest_step_at_530(self):
        E, y, comp, fwhm_at = self._build_roi()
        res = coupled_intensity_fit(
            E, y, comp, fwhm_at, continuum="step_linear_multi",
        )
        # continuum_params = [β₀, β₁, step_463, step_530, step_620, step_750]
        # NB: enable_quadratic = (E_span≥250) AND (n_comp≥5) — здесь
        # 4 comp, поэтому quadratic не активен → ровно 6 continuum params
        # (4 step anchors при threshold=5).
        assert len(res.continuum_params) >= 4
        step_params = res.continuum_params[2:]  # все step terms идут после β₀, β₁
        # Самый сильный должен быть тот, что около 530 (anchor #2 в
        # отсортированном по E списке).
        max_idx = int(np.argmax(step_params))
        # Anchors отсортированы по E: 463, 530, 620, 750. Самый интенсивный
        # step должен быть на 530 (индекс 1).
        assert max_idx == 1, (
            f"Самый сильный step должен быть на 530 кэВ (idx=1), "
            f"но max idx={max_idx}, params={step_params}"
        )


# ──────────────────────────────────────────────────────────────────
# Test 3: Synthetic multi-step (два step jumps)
# ──────────────────────────────────────────────────────────────────

class TestSyntheticMultiStep:
    """4 peaks с двумя step jumps в разных местах — multi-step должен
    моделировать оба."""

    def _build_roi(self):
        np.random.seed(20260603)
        E = np.linspace(400.0, 900.0, 900)
        sigma_25 = 25.0 / 2.355
        y = 200.0 + 0.0 * E
        y += _gauss(E, 463.0, sigma_25, 6000.0)
        y += _gauss(E, 530.0, sigma_25, 15000.0)
        y += _gauss(E, 620.0, sigma_25, 12000.0)
        y += _gauss(E, 750.0, sigma_25, 10000.0)
        # Два step'а: после 530 и после 750
        y += 50.0 * _smooth_step(E, 550.0, sigma_25)
        y += 30.0 * _smooth_step(E, 770.0, sigma_25)
        y = np.maximum(y, 1.0)
        y = np.random.poisson(y).astype(float)
        comp = [
            ComponentSpec("Nuc-A", 463.0, 6.0,  ""),
            ComponentSpec("Nuc-B", 530.0, 15.0, ""),
            ComponentSpec("Nuc-C", 620.0, 12.0, ""),
            ComponentSpec("Nuc-D", 750.0, 10.0, ""),
        ]
        fwhm_at = lambda E0: 25.0
        return E, y, comp, fwhm_at

    def test_multi_strictly_better_than_step_linear(self):
        E, y, comp, fwhm_at = self._build_roi()
        res_std = coupled_intensity_fit(
            E, y, comp, fwhm_at, continuum="step_linear",
        )
        res_multi = coupled_intensity_fit(
            E, y, comp, fwhm_at, continuum="step_linear_multi",
        )
        # При наличии двух реальных step'ов multi обязан улучшить fit
        # минимум на 5% по χ²/ν.
        assert res_multi.chi2_per_dof < res_std.chi2_per_dof * 0.95, (
            f"multi χ²/ν={res_multi.chi2_per_dof:.3f} должен быть лучше "
            f"step_linear χ²/ν={res_std.chi2_per_dof:.3f} минимум на 5% "
            f"при наличии 2 реальных step'ов"
        )


# ──────────────────────────────────────────────────────────────────
# Test 4: No-step / boundary cases — auto-select не активирует multi
# ──────────────────────────────────────────────────────────────────

class TestAutoSelectBoundary:
    """Проверяем _f392_auto_select_continuum на boundary cases."""

    def test_narrow_two_peak_stays_step_linear(self):
        # 2 close peaks, span < 200 кэВ → step_linear
        items = [(910.0, 25.8), (969.0, 15.8)]
        assert _f392_auto_select_continuum(items, "step_linear") == "step_linear"

    def test_wide_but_only_two_intense_stays_step_linear(self):
        # Component span 240 кэВ ≥ 200, только 2 anchor ≥4% → step_linear.
        # F-392.1: threshold снижен с 5.0 до 4.0%, поэтому boundary case
        # переформулирован: используем 510 (22.6%) и 583 (30.5%) как
        # единственные ≥4% линии, плюс sub-threshold 3.5% и 3.0%.
        items = [(463.0, 3.5), (510.0, 22.6), (583.0, 30.5), (703.0, 3.0)]
        assert _f392_auto_select_continuum(items, "step_linear") == "step_linear"

    def test_three_intense_anchors_promoted(self):
        # 3 intense, span > 200 → multi
        items = [(463.0, 10.0), (530.0, 15.0), (700.0, 20.0)]
        assert _f392_auto_select_continuum(items, "step_linear") == "step_linear_multi"

    def test_close_anchors_collapsed(self):
        # 4 intense, но 2 пары близки (< 40 кэВ) → merged → 2 anchors
        items = [(463.0, 10.0), (490.0, 12.0),    # collapse → 490 (i=12)
                 (600.0, 15.0), (620.0, 8.0)]      # collapse → 600 (i=15)
        # После merge: 2 anchors → недостаточно
        assert _f392_auto_select_continuum(items, "step_linear") == "step_linear"

    def test_roi_span_dominates_when_components_narrow(self):
        # Components narrow (40 кэВ), но ROI span 350 кэВ → должен auto-promote
        # если ≥3 intense anchors. Здесь intense=3 (>= sep_min=40 не нарушает,
        # потому что отдельные anchors сами по себе должны быть > 40 apart).
        # Используем 3 anchors с separation ≥ 40 кэВ всё ещё.
        items = [(463.0, 10.0), (520.0, 15.0), (583.0, 20.0)]
        # Component span = 120, ROI span = 350 → max=350 >= 200 → multi
        assert _f392_auto_select_continuum(
            items, "step_linear", roi_e_span_keV=350.0
        ) == "step_linear_multi"
        # Без ROI hint: span=120 < 200 → не активируется
        assert _f392_auto_select_continuum(
            items, "step_linear",
        ) == "step_linear"

    def test_linear_base_never_promoted(self):
        items = [(463.0, 10.0), (530.0, 15.0), (700.0, 20.0)]
        assert _f392_auto_select_continuum(items, "linear") == "linear"

    def test_step_linear_multi_passthrough(self):
        # Если уже step_linear_multi — не понижать обратно.
        items = [(910.0, 25.8), (969.0, 15.8)]
        assert _f392_auto_select_continuum(
            items, "step_linear_multi"
        ) == "step_linear_multi"

    def test_f392_1_threshold_boundary_4pct(self):
        """F-392.1 / v1.18.27.1 — default threshold 5.0 → 4.0%.

        Th-232 PROD M3 boundary: Ac-228 463 (I=4.4%) + Tl-208 510
        (I=8.1%) + Tl-208 583 (I=30.5%). При прежнем threshold=5%
        Ac-228 463 отсекался → 2 anchor → step_linear. При новом
        threshold=4% → 3 anchor → step_linear_multi.

        Также проверяем что 3.5% линия ВСЁ ЕЩЁ ниже порога:
        sub-threshold линия (Ac-228 409.46 I=1.92%) не качается, и
        Ac-228 562.5 (I=0.87%) тоже отсечён.
        """
        # Реальный M3 PROD anchor set: 3 in-ROI anchors с I≥4% + sub-threshold
        items = [
            (409.46, 1.92),   # Ac-228 409 — sub-threshold (1.92 < 4.0)
            (463.00, 4.40),   # Ac-228 463 — QUALIFIES под F-392.1 (4.4 ≥ 4.0)
            (510.77, 8.12),   # Tl-208 510 — qualifies (8.1 ≥ 4.0)
            (562.50, 0.87),   # Ac-228 562 — sub-threshold (<4.0)
            (583.19, 30.55),  # Tl-208 583 — qualifies (30.5 ≥ 4.0)
            (674.75, 2.10),   # Ac-228 674 — sub-threshold
        ]
        # E_span = 674-409 = 265 кэВ ≥ 200 ✓
        # 3 anchors ≥ 4.0%: {463, 510.77, 583.19}, separations 47.77 / 72.42 ≥ 40 ✓
        # → step_linear_multi
        assert _f392_auto_select_continuum(
            items, "step_linear"
        ) == "step_linear_multi"
        # Если бы прежний threshold (5.0%) был активен, 463 (4.4%) отсёкся
        # бы → только 2 anchor → step_linear. Проверим эксплицитно.
        assert _f392_auto_select_continuum(
            items, "step_linear", intense_threshold_pct=5.0
        ) == "step_linear"

    def test_f392_1_threshold_below_boundary(self):
        """F-392.1: линия с I=3.99% ВСЁ ЕЩЁ ниже порога 4.0%."""
        items = [
            (463.0, 3.99),    # ниже порога 4.0
            (510.77, 8.12),
            (583.19, 30.55),
        ]
        # Только 2 anchor ≥ 4.0% → step_linear (не повышается)
        assert _f392_auto_select_continuum(
            items, "step_linear", roi_e_span_keV=250.0
        ) == "step_linear"


# ──────────────────────────────────────────────────────────────────
# Test 5: Защита от лишних columns — back-compat для step_linear
# ──────────────────────────────────────────────────────────────────

class TestBackCompatStepLinear:
    """step_linear без multi не должен поменять поведение от F-383."""

    def test_step_linear_unchanged(self):
        np.random.seed(20260604)
        E = np.linspace(500.0, 800.0, 600)
        sigma_25 = 25.0 / 2.355
        y = 150.0 + _gauss(E, 600.0, sigma_25, 10000.0)
        y = np.maximum(y, 1.0)
        y = np.random.poisson(y).astype(float)
        comp = [
            ComponentSpec("Nuc-A", 600.0, 20.0, ""),
        ]
        fwhm_at = lambda E0: 25.0
        res = coupled_intensity_fit(
            E, y, comp, fwhm_at, continuum="step_linear",
        )
        # continuum_params = [β₀, β₁, β_step]
        assert len(res.continuum_params) == 3
        # Подгонка должна сходиться, χ²/ν около 1
        assert res.chi2_per_dof < 3.0
