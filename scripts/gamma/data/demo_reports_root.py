# -*- coding: utf-8 -*-
"""F-384 / v1.18.25.3 — first-run check для папки demo_reports.

Контракт навсегда:
* Демо-отчёты НЕ входят в релизный архив (см. build_release_archive.py
  EXCLUDE_DIRS), но скилл хочет знать куда писать отчёты при ``--full-report``.
* По умолчанию папка ``demo_reports/`` создаётся в **корне скилла**.
* Переопределение через переменную окружения ``GAMMA_DEMO_REPORTS_DIR``.
* При первом запуске CLI (`gamma analyze --full-report`):
    1. Если папка существует — используется молча.
    2. Если переменная ``GAMMA_DEMO_REPORTS_DIR`` задана — берётся оттуда,
       создаётся при отсутствии (молча).
    3. Если stdin это TTY (interactive) — пользователь видит prompt с
       default-предложением, может ввести альтернативный путь.
    4. В non-interactive режиме (CI, pipeline, scripts) — создаётся
       default ``<skill_root>/demo_reports/`` без вопросов.

API:
    ensure_demo_reports_root() -> Path
        Возвращает гарантированно-существующий Path к demo_reports.
        В non-interactive и при заданном env var НЕ выводит ничего лишнего.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path


def _skill_root() -> Path:
    """Корень скилла. cli.py живёт в scripts/gamma/, root = parents[2]."""
    return Path(__file__).resolve().parents[3]


_DEFAULT_NAME = "demo_reports"
_ENV_VAR = "GAMMA_DEMO_REPORTS_DIR"


def get_demo_reports_root_default() -> Path:
    """Default путь: <skill_root>/demo_reports."""
    return _skill_root() / _DEFAULT_NAME


def ensure_demo_reports_root(interactive: bool = True) -> Path:
    """First-run check. Возвращает Path с гарантией существования.

    Порядок:
    1. env var ``GAMMA_DEMO_REPORTS_DIR`` — приоритет, создаётся при
       отсутствии, без prompt.
    2. ``<skill_root>/demo_reports/`` — default. Если уже существует —
       молча возвращает.
    3. Default не существует И stdin это TTY И ``interactive=True`` —
       prompt пользователю с default-предложением.
    4. Default не существует И НЕ TTY (или ``interactive=False``) —
       создаёт default молча.

    Параметры:
        interactive: если False — НЕ запускать prompt даже в TTY
            (используется в тестах/CI).
    """
    # 1. env var
    env_path = os.environ.get(_ENV_VAR)
    if env_path:
        p = Path(env_path).expanduser().resolve()
        p.mkdir(parents=True, exist_ok=True)
        return p

    # 2. default exists
    default = get_demo_reports_root_default()
    if default.exists() and default.is_dir():
        return default

    # 3. interactive prompt
    if interactive and sys.stdin.isatty() and sys.stderr.isatty():
        try:
            print(
                f"\n[gamma] Папка для отчётов не найдена.\n"
                f"        По умолчанию будет создана: {default}",
                file=sys.stderr,
            )
            print(
                "        Введите альтернативный путь или Enter для default: ",
                end="", file=sys.stderr, flush=True,
            )
            user = sys.stdin.readline().strip()
        except (EOFError, KeyboardInterrupt):
            user = ""
        if user:
            chosen = Path(user).expanduser().resolve()
        else:
            chosen = default
    else:
        # 4. non-interactive
        chosen = default

    chosen.mkdir(parents=True, exist_ok=True)
    print(f"[gamma] demo_reports → {chosen}", file=sys.stderr)
    return chosen


__all__ = [
    "ensure_demo_reports_root",
    "get_demo_reports_root_default",
]
