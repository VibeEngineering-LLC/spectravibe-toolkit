"""F-127 / v1.17.7 — per-line T(E) tail calibration для NaI 63×63."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.peaks.coupled_multiplet import (
    nai_tail_T_at, NAI_T_E_T_REF, NAI_T_E_REF_KEV,
    NAI_T_E_T_MIN, NAI_T_E_T_MAX,
    ComponentSpec, coupled_intensity_fit,
)


def test_nai_T_at_reference_point():
    """T(662 кэВ) = T_REF = 0.7 точно (определение калибровки)."""
    assert abs(nai_tail_T_at(662.0) - 0.7) < 1e-6


def test_nai_T_at_monotonic_with_energy():
    """T(E) монотонно растёт с энергией (хвост становится короче)."""
    T_60 = nai_tail_T_at(60.0)
    T_300 = nai_tail_T_at(300.0)
    T_662 = nai_tail_T_at(662.0)
    T_1461 = nai_tail_T_at(1461.0)
    T_2614 = nai_tail_T_at(2614.0)
    assert T_60 <= T_300 <= T_662 <= T_1461 <= T_2614


def test_nai_T_at_clamps_low_E():
    """T(<10 кэВ) clamp на T_MIN."""
    T = nai_tail_T_at(5.0)
    assert T == NAI_T_E_T_MIN


def test_nai_T_at_clamps_high_E():
    """T(>10 МэВ) clamp на T_MAX."""
    T = nai_tail_T_at(20000.0)
    assert T == NAI_T_E_T_MAX


def test_nai_T_at_zero_energy_returns_ref():
    """T(E≤0) → T_ref (защита от деления на ноль и log(0))."""
    assert nai_tail_T_at(0.0) == NAI_T_E_T_REF
    assert nai_tail_T_at(-100.0) == NAI_T_E_T_REF


def test_coupled_fit_uses_T_at_callable():
    """coupled_intensity_fit принимает tail_T_at callable и применяет."""
    E = np.linspace(900.0, 1000.0, 200)
    counts = np.full_like(E, 50.0)
    # 1 компонента, фит должен пройти без ошибок и записать
    # basis_label = peak_image_T(E)
    res = coupled_intensity_fit(
        E, counts, [ComponentSpec("Test", 950.0, 100.0)],
        lambda x: 10.0,
        use_peak_image=True,
        tail_T_at=nai_tail_T_at,
    )
    assert "peak_image_T(E)" in res.notes


def test_coupled_fit_hardcoded_T_when_no_callable():
    """Без tail_T_at callable — старое поведение с T=tail_param хардкодом."""
    E = np.linspace(900.0, 1000.0, 200)
    counts = np.full_like(E, 50.0)
    res = coupled_intensity_fit(
        E, counts, [ComponentSpec("Test", 950.0, 100.0)],
        lambda x: 10.0,
        use_peak_image=True, tail_param=0.5,
    )
    assert "T=0.50" in res.notes or "peak_image_T=0.50" in res.notes
