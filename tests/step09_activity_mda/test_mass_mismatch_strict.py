# -*- coding: utf-8 -*-
"""F-378 / v1.18.25 — strict mass mismatch detection.

Регрессия: пользователь запускал gamma.cli с --sample-mass-kg=0.5,
у .spe в SAMPLEMASS лежало 1.6 кг. Скил молча использовал CLI value
→ удельная активность Th-232 = 6237 Бк/кг вместо паспортных 1940
(× 3.2). Паспортная сверка дала +221.5 % отклонение — внешне
выглядит как промах эталона, на деле — конфликт массы.

F-378 эмитирует ⚠ WARN в notes + stderr, когда CLI mass и .spe
SAMPLEMASS расходятся > 1 %.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from gamma.identification.staged_pipeline import check_sample_mass_mismatch


# ═════════════════════════════════════════════════════════════════════
# Конфликт обнаружен
# ═════════════════════════════════════════════════════════════════════


def test_F378_detects_user_bug_0_5_vs_1_6():
    """Реальный случай из bug-репорта: CLI=0.5 кг, .spe=1.6 кг."""
    note = check_sample_mass_mismatch(
        cli_mass_kg=0.5,
        spec_extras={"lsrm_sample_mass_kg": 1.6},
    )
    assert note is not None, "должно сработать на пользовательском кейсе"
    assert "F-378" in note
    assert "0.500" in note  # CLI mass в форматированном виде
    assert "1.600" in note  # .spe mass в форматированном виде
    # Разница (1.6 - 0.5) / 1.6 = 0.6875 → "68.8 %"
    assert "68.8" in note or "68.7" in note  # rel_diff в процентах


def test_F378_factor_is_cli_over_spec():
    """factor ≡ cli/spec — показывает «во сколько раз искажена SA»."""
    note = check_sample_mass_mismatch(
        cli_mass_kg=0.5,
        spec_extras={"lsrm_sample_mass_kg": 1.6},
    )
    # SA = A/m. Если cli<spec, SA(cli) > SA(spec) → завышена.
    # Factor = cli/spec показывает направление и величину искажения.
    expected = 0.5 / 1.6  # 0.3125
    assert f"{expected:.2f}×" in note


def test_F378_emits_when_difference_above_1pct():
    """Минимальный порог — 1 % относительной разницы."""
    # 1.005 кг vs 1.00 кг = 0.5 % разницы → ниже порога, не триггерится
    assert check_sample_mass_mismatch(
        cli_mass_kg=1.005,
        spec_extras={"lsrm_sample_mass_kg": 1.0},
    ) is None
    # 1.02 кг vs 1.00 кг = 2 % → выше порога, триггерится
    note = check_sample_mass_mismatch(
        cli_mass_kg=1.02,
        spec_extras={"lsrm_sample_mass_kg": 1.0},
    )
    assert note is not None and "F-378" in note


def test_F378_huge_factor_call_out():
    """Огромная разница (10×) — фактор виден ярко."""
    note = check_sample_mass_mismatch(
        cli_mass_kg=10.0,
        spec_extras={"lsrm_sample_mass_kg": 1.0},
    )
    assert note is not None
    assert "10.00×" in note


# ═════════════════════════════════════════════════════════════════════
# Конфликта нет (или невозможен)
# ═════════════════════════════════════════════════════════════════════


def test_F378_no_warn_when_cli_mass_absent():
    """Без CLI флага — нет конфликта (включается auto-F-140)."""
    assert check_sample_mass_mismatch(
        cli_mass_kg=None,
        spec_extras={"lsrm_sample_mass_kg": 1.6},
    ) is None


def test_F378_no_warn_when_spec_mass_absent():
    """В .spe нет SAMPLEMASS → не с чем сравнивать."""
    assert check_sample_mass_mismatch(
        cli_mass_kg=0.5,
        spec_extras={},
    ) is None
    assert check_sample_mass_mismatch(
        cli_mass_kg=0.5,
        spec_extras=None,
    ) is None


def test_F378_no_warn_when_values_match():
    """Точное совпадение → нет warn."""
    assert check_sample_mass_mismatch(
        cli_mass_kg=1.6,
        spec_extras={"lsrm_sample_mass_kg": 1.6},
    ) is None


def test_F378_no_warn_within_tolerance():
    """Расхождение ≤ 1 % (float rounding) — не warn."""
    # 1.6 vs 1.605 = 0.31 % разницы
    assert check_sample_mass_mismatch(
        cli_mass_kg=1.605,
        spec_extras={"lsrm_sample_mass_kg": 1.6},
    ) is None


def test_F378_tolerance_configurable():
    """rel_tol можно ужесточить (используется в snapshot-тестах)."""
    # 0.5 % разница: проходит при default 1 %, но при tol=0.1 % — warn
    assert check_sample_mass_mismatch(
        cli_mass_kg=1.005,
        spec_extras={"lsrm_sample_mass_kg": 1.0},
        rel_tol=0.01,  # 1 %
    ) is None
    note = check_sample_mass_mismatch(
        cli_mass_kg=1.005,
        spec_extras={"lsrm_sample_mass_kg": 1.0},
        rel_tol=0.001,  # 0.1 %
    )
    assert note is not None


# ═════════════════════════════════════════════════════════════════════
# Защита от мусорных данных
# ═════════════════════════════════════════════════════════════════════


def test_F378_handles_invalid_spec_mass():
    """Невалидный тип в extras — silently None, без падения."""
    assert check_sample_mass_mismatch(
        cli_mass_kg=0.5,
        spec_extras={"lsrm_sample_mass_kg": "не-число"},
    ) is None
    assert check_sample_mass_mismatch(
        cli_mass_kg=0.5,
        spec_extras={"lsrm_sample_mass_kg": None},
    ) is None


def test_F378_handles_zero_or_negative_spec_mass():
    """Нулевая или отрицательная масса в .spe — невалидно, не warn."""
    assert check_sample_mass_mismatch(
        cli_mass_kg=0.5,
        spec_extras={"lsrm_sample_mass_kg": 0.0},
    ) is None
    assert check_sample_mass_mismatch(
        cli_mass_kg=0.5,
        spec_extras={"lsrm_sample_mass_kg": -1.0},
    ) is None
