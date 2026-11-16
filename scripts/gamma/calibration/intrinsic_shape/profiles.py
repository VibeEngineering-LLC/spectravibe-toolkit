"""Профили детектора: пресеты, метаданные, автоопределение.
Профиль = набор пресетов метода (опоры эталона, register, ecal, lowE, окна структур, selftest, ПШПВ, допуски) + метаданные прибора;
поле method -- метод калибровки профиля (La-138 -- "la138_intrinsic"; другие кристаллы -- свой метод, команды метода La-138 их отвергают).
Сгенерировано gen_code.py (qwen3.8:27b) по спецификации audit/specs/profiles_spec.md, 2026-11-15; правка: global в using()."""
from __future__ import annotations

import json
import re
from contextlib import contextmanager
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRESETS_DIR = HERE / "presets"
PROFILES_DIR = PRESETS_DIR / "profiles"

# id профиля начинается с кристалла: "labr_<прибор>" (LaBr3), "nai_<прибор>" (NaI, другой метод); в "title" -- та же
# приставка «LaBr »/«NaI ». Прежние id -- в "aliases" профиля: принимаются с предупреждением «устаревшее имя» (resolve).
DEFAULT_ID = "labr_atomnano5pro"
KINDS = ("reference", "ecal", "lowE", "struct_windows", "selftest")

_ACTIVE: str | None = None  # id активного профиля; None — по умолчанию


def _nai_profiles() -> list[dict]:
    """Профили метода NaI (nai_shape/profiles/*.json) в общем виде: id, title, aliases, method "nai", status; пресетов La-138 нет.
    Числа метода остаются в файле профиля nai_shape и читаются nai_shape.api.load_profile."""
    from gamma.calibration.nai_shape.api import PROFILE_ALIASES, PROFILE_DIR
    out: list[dict] = []
    for f in sorted(PROFILE_DIR.glob("nai_*.json")):
        p = json.loads(f.read_text("utf-8"))
        st = p.get("status", "")
        out.append({"id": p["id"], "title": p.get("device", p["id"]) + (" [экспериментально]" if st == "experimental" else ""),
                    "aliases": sorted(a for a, v in PROFILE_ALIASES.items() if v == p["id"]),
                    "crystal": p.get("crystal"), "size": None, "n_ch": p.get("n_channels"), "method": p.get("method"),
                    "status": st, "status_note": p.get("status_note", ""), "default": False, "match": [], "presets": {}})
    return out


def list_profiles(all_methods: bool = False) -> list[dict]:
    """Профили метода La-138 из PROFILES_DIR (по имени файла); all_methods=True -- плюс профили NaI (метод "nai")."""
    out: list[dict] = []
    for f in sorted(PROFILES_DIR.glob("*.json")):
        prof = json.loads(f.read_text("utf-8"))
        if "id" not in prof or prof["id"] != f.stem:
            raise ValueError(f"некорректный профиль '{f.name}': id не совпадает с именем файла")
        out.append(prof)
    return out + (_nai_profiles() if all_methods else [])


def ids(aliases: bool = False, all_methods: bool = False) -> list[str]:
    """Список id профилей; aliases=True -- плюс устаревшие имена; all_methods=True -- и профили NaI (для выбора в CLI)."""
    ps = list_profiles(all_methods)
    return [p["id"] for p in ps] + ([a for p in ps for a in p.get("aliases", [])] if aliases else [])


def get(pid: str | None = None) -> dict:
    """Профиль по id или устаревшему имени (глубокая копия); ищет среди профилей всех методов (La-138 и NaI)."""
    if pid is None:
        pid = active_id()
    for p in list_profiles(all_methods=True):
        if p["id"] == pid or pid in p.get("aliases", []):  # устаревшее имя -- тот же профиль
            return json.loads(json.dumps(p))
    raise KeyError(f"неизвестный профиль детектора '{pid}'; есть: {', '.join(ids(all_methods=True))}")


def active_id() -> str:
    """Текущий активный id."""
    return _ACTIVE if _ACTIVE is not None else DEFAULT_ID


def activate(pid: str | None) -> str:
    """Установить активный профиль; None — сброс к умолчанию."""
    global _ACTIVE
    if pid is None:
        _ACTIVE = None
        return DEFAULT_ID
    _ACTIVE = get(pid)["id"]   # устаревшее имя -> канонический id
    return _ACTIVE


@contextmanager
def using(pid: str | None):
    """Временно переключает активный профиль."""
    global _ACTIVE
    prev = _ACTIVE
    activate(pid)
    try:
        yield pid
    finally:
        _ACTIVE = prev


def preset_path(kind: str, pid: str | None = None) -> Path:
    """Путь к пресету."""
    if kind not in KINDS:
        raise ValueError(f"неизвестный вид пресета '{kind}'; допустимы: {', '.join(KINDS)}")
    prof = get(pid)
    if kind not in prof.get("presets", {}):
        raise KeyError(f"в профиле '{prof['id']}' нет пресета '{kind}'")
    return PRESETS_DIR / prof["presets"][kind]


def load_preset(kind: str, pid: str | None = None) -> dict:
    """Загрузить пресет из JSON."""
    return json.loads(preset_path(kind, pid).read_text("utf-8"))


def detect(device_name: str | None) -> tuple[str | None, str]:
    """Определить профиль по имени прибора."""
    if not device_name:
        return None, "имя прибора в файле не указано"
    hits: list[str] = []
    for p in list_profiles():
        for pat in p.get("match", []):
            if re.search(pat, device_name, re.IGNORECASE):
                hits.append(p["id"])
                break
    if len(hits) == 1:
        return hits[0], f"прибор '{device_name}' -> профиль '{hits[0]}'"
    if len(hits) > 1:
        return None, f"прибор '{device_name}' подходит под несколько профилей: {', '.join(hits)}"
    return None, f"прибор '{device_name}' не распознан ни одним профилем ({', '.join(ids())})"


def resolve(requested: str | None = None, device_name: str | None = None,
            n_ch: int | None = None) -> tuple[str, list[str]]:
    """Выбор профиля для прогона; возвращает (pid, warnings)."""
    warnings: list[str] = []
    if requested is None:
        requested = "auto"

    if requested != "auto":
        pid = get(requested)["id"]
        if pid != requested:
            warnings.append(f"устаревшее имя профиля '{requested}' -- используйте '{pid}'")
        det, _ = detect(device_name)
        if det is not None and det != pid:
            warnings.append(
                f"выбран профиль '{requested}', а прибор в файле '{device_name}' распознан как '{det}'"
            )
    else:
        det, msg = detect(device_name)
        if det is not None:
            pid = det
        else:
            pid = DEFAULT_ID
            warnings.append(f"{msg}; взят профиль по умолчанию '{DEFAULT_ID}' -- проверьте выбор (--detector)")

    # n_ch -- справочно: число каналов профиля -- каналы его эталона, автопресет переводит пресеты на любое число
    # каналов (AtomNano 5 PRO -- 2048 и 8192), поэтому несовпадение числа каналов не признак чужого профиля
    return pid, warnings


def require_method(method: str, pid: str | None = None) -> None:
    """Проверить совместимость метода с профилем."""
    prof = get(pid)
    if prof.get("method") != method:
        raise ValueError(
            f"профиль '{prof['id']}' использует метод '{prof.get('method')}', "
            f"а эта команда реализует '{method}'"
        )


def describe(pid: str | None = None) -> str:
    """Однострочное описание профиля."""
    p = get(pid)

    def _f(key: str) -> str:
        v = p.get(key)
        return str(v) if v is not None else "?"

    return f"{p['id']}: {_f('title')} ({_f('crystal')}, {_f('size')}, {_f('n_ch')} кан, метод {_f('method')})"
