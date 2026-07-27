# -*- coding: utf-8 -*-
"""Собрать SPECTRA_MAP.md из индекса git.

Источник истины — `git ls-files`: в карту попадает только то, что
действительно опубликовано, поэтому цифры не могут разойтись с репозиторием.
Ссылки на файлы разрешаются по факту, а не по записям индексов классов.

Usage:
    python scripts/convert/build_spectra_map.py
"""
from __future__ import annotations

import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "SPECTRA_MAP.md"

DETECTOR_TYPE = {
    "Gamma-1S": "NaI(Tl) 63×63, УДС-ГЦ",
    "NM_HPGe20": "HPGe коаксиальный 20 %",
    "GP_HPGe20": "HPGe 20 % общего назначения",
    "Handy_HPGe": "HPGe GMX",
    "Handy_LaBr": "LaBr₃",
    "Handy_NaI": "NaI",
    "Simple_HPGe": "HPGe планарный, демо",
    "Simple_NaI": "NaI",
    "Simple_TeCd": "CdTe / CZT",
    "AtomSpectra": "сцинтилляционный, фикстуры",
    "RadiaCode_103": "CsI(Tl), апатит и камень",
}


def tracked() -> list[str]:
    out = subprocess.run(["git", "-C", str(REPO), "ls-files", "-z"],
                         capture_output=True, check=True)
    return [n for n in out.stdout.decode("utf-8").split("\0") if n]


def size_of(rel: str) -> int:
    p = REPO / rel
    try:
        return p.stat().st_size
    except OSError:
        return 0


def collect():
    per = defaultdict(lambda: {"spe": 0, "xml": 0, "bytes": 0})
    for rel in tracked():
        parts = rel.split("/")
        if len(parts) < 2 or parts[0] != "detectors":
            continue
        cls = parts[1]
        per[cls]["bytes"] += size_of(rel)
        if rel.endswith(".spe"):
            per[cls]["spe"] += 1
        elif rel.endswith(".xml"):
            per[cls]["xml"] += 1
    return per


def degree5():
    """Файлы с полиномом 5-й степени — по флагу из INDEX.json классов."""
    rows = []
    tracked_set = set(tracked())
    for rel in sorted(tracked_set):
        if not rel.endswith("INDEX.json") or not rel.startswith("detectors/"):
            continue
        cls = rel.split("/")[1]
        try:
            data = json.loads((REPO / rel).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        base = rel.rsplit("/", 1)[0]
        for e in data.get("entries", []):
            if not e.get("becqmoni_reads_as_linear"):
                continue
            name = e.get("xml") or e.get("spe")
            if not name:
                continue
            full = f"{base}/becqmoni/{name}"
            rows.append((cls, *resolve(full, tracked_set),
                         e.get("energy_cal_degree")))
    return rows


def resolve(full: str, tracked_set: set) -> tuple[str, str]:
    """Найти файл по записи индекса и сказать, опубликован он или локален.

    Индексы классов могут отставать от текущих имён, поэтому прямой путь не
    всегда существует. Ищем по имени с подстановкой на месте метки прибора или
    источника, а файл, лежащий на диске но не в индексе git, помечаем как
    локальный, а не как потерянный.
    """
    import re as _re
    if full in tracked_set:
        return full, "в git"
    directory, _, name = full.rpartition("/")
    pattern = _re.escape(name).replace(r"SRC\-01", r"SRC\-\d+") \
                              .replace(r"SN\-01", r"SN\-\d+")
    pattern = _re.sub(r"(?:SRC|SN)\\-\d+", r"(?:SRC|SN)\\-\\d+", _re.escape(name))
    rx = _re.compile("^" + pattern + "$")
    for cand in tracked_set:
        if cand.startswith(directory + "/") and rx.match(cand.rpartition("/")[2]):
            return cand, "в git"
    disk = REPO / full
    if disk.exists():
        return full, "локально"
    parent = REPO / directory
    if parent.is_dir():
        for cand in parent.iterdir():
            if rx.match(cand.name):
                return f"{directory}/{cand.name}", "локально"
    return full, "не найден"


def local_only() -> dict:
    """Что лежит на диске, но исключено .gitignore."""
    counts = {}
    for sub in ("lsrm", "becqmoni"):
        d = REPO / "detectors" / "Gamma-1S" / "reference_spectra" / sub
        counts[sub] = sum(1 for p in d.rglob("*") if p.is_file()) if d.is_dir() else 0
    return counts


def main() -> int:
    per = collect()
    d5 = degree5()
    loc = local_only()
    tot_spe = sum(v["spe"] for v in per.values())
    tot_xml = sum(v["xml"] for v in per.values())
    tot_mb = sum(v["bytes"] for v in per.values()) / 1024 / 1024

    L = []
    A = L.append
    A("# Карта спектров\n")
    A("Собрано генератором `scripts/convert/build_spectra_map.py` по индексу git —")
    A("здесь только то, что действительно опубликовано. Подробные описи с паспортными")
    A("активностями лежат в `detectors/<класс>/reference_spectra/README.md` и `INDEX.json`.\n")
    A(f"Всего в репозитории: **{tot_xml} файлов BecqMoni XML** и {tot_spe} исходных `.spe`")
    A(f"по {len(per)} классам детекторов, {tot_mb:.0f} МБ.\n")

    A("## Структура\n")
    A("```")
    A("spectravibe-toolkit/")
    A("├── SKILL.md, ARCH.md, FORMAT_REGISTRY.md, NOTES_v1.7_methodology.md")
    A("├── scripts/          библиотека gamma/, CLI-утилиты, конвертеры")
    A("├── tests/            тесты по шагам конвейера")
    A("├── data/             алиасы, нуклиды, XRF-линии, tsc_lookup")
    A("├── references/       контракт демо-отчёта, кэш IAEA, ground-truth ЛСРМ")
    A(f"└── detectors/        {tot_xml} XML · {tot_spe} .spe")
    A("```\n")

    A("## Спектры по классам\n")
    A("| Класс | Тип детектора | `.spe` | **XML** | МБ |")
    A("|---|---|---:|---:|---:|")
    for cls in sorted(per, key=lambda c: -per[c]["xml"]):
        v = per[cls]
        A(f"| `{cls}` | {DETECTOR_TYPE.get(cls, '—')} | {v['spe'] or '—'} "
          f"| **{v['xml']}** | {v['bytes']/1024/1024:.1f} |")
    A(f"| **Итого** | | **{tot_spe}** | **{tot_xml}** | **{tot_mb:.1f}** |\n")

    if any(loc.values()):
        A("### Что лежит локально и в git не идёт\n")
        A("В репозитории публикуются только **курированные** наборы Gamma-1S —")
        A("`reference_kits/` и `reference_kits_becqmoni/`: пары «образец + фон» по")
        A("геометриям, с описью, паспортными активностями и кривыми эффективности.")
        A("Массовая выгрузка рабочего дерева ЛСРМ остаётся на машине оператора.\n")
        for k, v in loc.items():
            A(f"- `detectors/Gamma-1S/reference_spectra/{k}/` — {v} файлов")
        A("")
        A("Пересобирается из исходного дерева одной командой:\n")
        A("```bash")
        A("python scripts/convert/lsrm_tree_publish.py --class Gamma-1S")
        A("```\n")

    A("### Почему у NM_HPGe20 разрыв между `.spe` и XML\n")
    A("Не сконвертированы объёмные исследовательские коллекции — Pu/U LNHB, INR,")
    A("IPPPE, kpti, около 400 спектров. У HPGe 16 384 канала, один XML весит ~590 КБ")
    A("против 33 КБ исходника; это добавило бы ≈250 МБ производных файлов. Исходники")
    A("опубликованы, экспорт собирается по требованию:\n")
    A("```bash")
    A("python scripts/convert/lsrm_tree_publish.py --class NM_HPGe20 --xml-all")
    A("```\n")

    A("---\n")
    A(f"## ⚠️ ВАЖНО: {len(d5)} спектров BecqMoni прочитает с неверной шкалой\n")
    A("`BecquerelMonitor/PolynomialEnergyCalibration.cs` в методе `ChannelToEnergy()`")
    A("разбирает только степени полинома **4, 3 и 2**. Всё, что выше, молча уходит")
    A("в линейную ветку:\n")
    A("```csharp")
    A("if (this.polynomialOrder == 4) { ... }")
    A("if (this.polynomialOrder == 3) { ... }")
    A("if (this.polynomialOrder == 2) { ... }")
    A("return this.coefficients[1] * n + this.coefficients[0];   // сюда попадает 5-я степень")
    A("```\n")
    A("У этих спектров ЛСРМ записал калибровку **5-й степени**. Коэффициенты сохранены")
    A("как в источнике — но при открытии в BecqMoni энергетическая шкала будет неверной,")
    A("причём без единого сообщения об ошибке.\n")
    A("| Класс | Файл | Где | Степень |")
    A("|---|---|---|---:|")
    for cls, full, status, deg in d5:
        A(f"| `{cls}` | `{full.split('/')[-1]}` | {status} | {deg} |")
    A("")
    if any(s == "локально" for _, _, s, _ in d5):
        A("«локально» — файл лежит на машине оператора, в репозиторий не публикуется.\n")
    A("**Что делать.** Берите калибровку из `INDEX.json` соответствующего класса —")
    A("поле `energy_cal`, коэффициенты от младшей степени к старшей:\n")
    A("```")
    A("E(n) = a₀ + a₁·n + a₂·n² + a₃·n³ + a₄·n⁴ + a₅·n⁵")
    A("```\n")
    A("Найти такие файлы можно по флагу `becqmoni_reads_as_linear` в том же индексе.\n")
    A("**Почему не исправлено пересчётом.** Приведение к 4-й степени было реализовано")
    A("и отклонено: на HPGe-кривых оно ложится в пределах 0,21 кэВ, но на сильно")
    A("изогнутой NaI-калибровке лучшее приближение 4-й степени промахивается на 57 кэВ.")
    A("Молча менять сохранённую калибровку хуже, чем сообщить о проблеме — решение")
    A("оператора от 2026-07-27. Простое отбрасывание старшего члена ещё хуже: на канале")
    A("8192 оно даёт ошибку 53 кэВ, а на NaI-случае — более 14 000 кэВ.\n")

    A("---\n")
    A("## Прочие особенности форматов\n")
    A("Конвертеры обходят ещё два дефекта BecqMoni, оба задокументированы в коде:\n")
    A("- **Единицы `SampleInfo`** — `Weight` в килограммах, `Volume` в литрах, диапазон")
    A("  контрола 0,001…100. ЛСРМ даёт граммы и миллилитры; без деления на 1000")
    A("  приложение падает с `ArgumentOutOfRangeException` до отрисовки документа")
    A("  (`gamma/io/becqmoni_xml.py`, BUG-BQ2).")
    A("- **GUID прибора** выводится из имени детектора через `uuid5`, иначе каждая")
    A("  перегенерация создавала бы в BecqMoni новый прибор.\n")

    A("## Обозначения\n")
    A("Приборы и источники обозначены внутренними метками вида `SN-nn` и `SRC-nn` —")
    A("одна метка на один прибор или источник, устойчиво между сборками.\n")

    OUT.write_text("\n".join(L), encoding="utf-8")
    print(f"{OUT.name}: {tot_xml} XML, {tot_spe} .spe, {len(per)} классов, "
          f"{len(d5)} спектров со степенью 5")
    missing = [f for _, f, ok, _ in d5 if not ok]
    if missing:
        print(f"  ВНИМАНИЕ: {len(missing)} ссылок ведут в несуществующие файлы")
        for f in missing:
            print(f"    {f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
