# -*- coding: utf-8 -*-
"""Стабильный идентификатор прибора для BecqMoni.

BecqMoni сопоставляет спектр с записью прибора по `<Guid>`. Если выдавать
случайный идентификатор, при каждой пересборке экспорта появляется новый
прибор, а git видит изменившимися все файлы. Поэтому GUID детерминированный:
один и тот же прибор всегда получает один и тот же идентификатор.

Значение выводится из имени прибора и локального ключа:

    GUID = uuid5(NAMESPACE_URL, key + "/" + detector_id)

Ключ лежит в `spectravibe_guid_salt.txt` рядом с каталогом репозитория либо
задаётся переменной окружения `SPECTRAVIBE_GUID_SALT`, и в репозиторий не
входит. При первом запуске создаётся автоматически.
"""
from __future__ import annotations

import os
import secrets
import uuid
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
KEY_FILE = REPO.parent / "spectravibe_guid_salt.txt"
ENV_VAR = "SPECTRAVIBE_GUID_SALT"


def load_salt(create: bool = True) -> str:
    """Взять ключ из окружения или локального файла; при отсутствии — создать."""
    env = os.environ.get(ENV_VAR)
    if env:
        return env.strip()
    if KEY_FILE.exists():
        key = KEY_FILE.read_text(encoding="utf-8").strip()
        if key:
            return key
    if not create:
        raise FileNotFoundError(
            f"нет ключа: задайте {ENV_VAR} или создайте {KEY_FILE}")
    key = secrets.token_hex(32)
    KEY_FILE.write_text(key + "\n", encoding="utf-8")
    print(f"создан локальный ключ: {KEY_FILE}")
    print("  Без него идентификаторы приборов пересоберутся другими.")
    return key


def device_guid(detector_id: str, salt: str | None = None) -> str:
    """Детерминированный GUID прибора."""
    salt = salt if salt is not None else load_salt()
    return str(uuid.uuid5(uuid.NAMESPACE_URL, f"{salt}/{detector_id}"))
