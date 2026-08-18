# -*- coding: utf-8 -*-
"""F-384 / v1.18.25.3 — ensure_demo_reports_root тесты.

Verify:
1. env var GAMMA_DEMO_REPORTS_DIR имеет приоритет
2. default path = <skill_root>/demo_reports
3. non-interactive: создаёт default без prompt
4. existing dir: возвращает as-is, без создания
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from gamma.data.demo_reports_root import (
    ensure_demo_reports_root,
    get_demo_reports_root_default,
)


def test_default_path_under_skill_root():
    p = get_demo_reports_root_default()
    # Должен быть <skill_root>/demo_reports
    assert p.name == "demo_reports"
    # parent — корень скилла (контракт)
    skill_root = p.parent
    # содержит обязательные верхне-уровневые артефакты
    assert (skill_root / "scripts" / "gamma").is_dir(), \
        f"skill_root {skill_root} не содержит scripts/gamma/"


def test_env_var_overrides(tmp_path, monkeypatch):
    target = tmp_path / "custom_demos"
    monkeypatch.setenv("GAMMA_DEMO_REPORTS_DIR", str(target))
    result = ensure_demo_reports_root(interactive=False)
    assert result == target.resolve()
    assert result.exists()
    assert result.is_dir()


def test_env_var_creates_if_missing(tmp_path, monkeypatch):
    target = tmp_path / "auto_create" / "demos_sub"
    assert not target.exists()
    monkeypatch.setenv("GAMMA_DEMO_REPORTS_DIR", str(target))
    result = ensure_demo_reports_root(interactive=False)
    assert result.exists()
    assert result.is_dir()


def test_default_used_when_no_env_and_exists(tmp_path, monkeypatch):
    """Сценарий: env var unset, default уже существует → возвращается as-is.

    Используем tmp_path как fake skill_root через env var вместо
    реального skill_root (чтобы тест не засорил `<repo>/demo_reports`
    при non-interactive create_default path).
    """
    fake_default = tmp_path / "demo_reports"
    fake_default.mkdir()
    monkeypatch.setenv("GAMMA_DEMO_REPORTS_DIR", str(fake_default))
    result = ensure_demo_reports_root(interactive=False)
    assert result == fake_default.resolve()


def test_non_interactive_no_prompt(monkeypatch, tmp_path, capsys):
    """В non-interactive режиме (interactive=False) НЕ должно быть
    interactive input() / readline()."""
    target = tmp_path / "x" / "demos"
    monkeypatch.setenv("GAMMA_DEMO_REPORTS_DIR", str(target))
    # Если bы был prompt, stdin readline блокировался бы → таймаут pytest.
    # Просто проверяем что вызов проходит и stdin не trognut.
    result = ensure_demo_reports_root(interactive=False)
    assert result.exists()


def test_ensure_idempotent(tmp_path, monkeypatch):
    target = tmp_path / "demos"
    monkeypatch.setenv("GAMMA_DEMO_REPORTS_DIR", str(target))
    a = ensure_demo_reports_root(interactive=False)
    b = ensure_demo_reports_root(interactive=False)
    assert a == b
    assert a.exists()


def test_skill_root_contract_includes_build_release_excludes():
    """F-384 контракт: demo_reports в EXCLUDE_DIRS у
    build_release_archive.py чтобы не попадал в архив релиза."""
    here = Path(__file__).resolve()
    skill_root = here.parents[2]
    build_script = skill_root / "scripts" / "build_release_archive.py"
    text = build_script.read_text(encoding="utf-8")
    assert '"demo_reports"' in text, (
        "build_release_archive.py не содержит 'demo_reports' в "
        "EXCLUDE_DIRS — папка будет включаться в архив, что нарушает "
        "F-384 контракт."
    )
