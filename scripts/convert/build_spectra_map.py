# -*- coding: utf-8 -*-
"""Собрать SPECTRA_MAP.md из индекса git.

Источник истины — `git ls-files`: в карту попадает только то, что
действительно опубликовано, поэтому цифры не могут разойтись с репозиторием.
Ссылки на файлы разрешаются по факту, а не по записям индексов классов.

Usage:
    python scripts/convert/build_spectra_map.py
"""
from __future__ import annotations

import re
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
        # >= 3 частей: файлы прямо в detectors/ (CHANGELOG.md) — не класс.
        if len(parts) < 3 or parts[0] != "detectors":
            continue
        cls = parts[1]
        per[cls]["bytes"] += size_of(rel)
        if rel.endswith(".spe"):
            per[cls]["spe"] += 1
        elif rel.endswith(".xml"):
            per[cls]["xml"] += 1
    return per


BLOCK = re.compile(r"<(EnergySpectrum|BackgroundEnergySpectrum)>(.*?)</\1>",
                   re.DOTALL)
ORDER = re.compile(r"<PolynomialOrder>(\d+)</PolynomialOrder>")
MAX_ORDER = 4


def high_order():
    """Файлы, где BecqMoni отвергнет калибровку — по самим XML.

    Читаем оба блока, а не флаг из INDEX.json: BecqMoni проверяет спектр и
    вшитый фон как две независимые калибровки, и опись, построенная только по
    спектру, вторую половину пропускала.
    """
    tracked_set = set(tracked())
    rows = []
    for p in sorted(REPO.joinpath("detectors").rglob("*.xml")):
        rel = p.relative_to(REPO).as_posix()
        try:
            text = p.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        worst = {}
        for m in BLOCK.finditer(text):
            o = ORDER.search(m.group(2))
            if o:
                worst[m.group(1)] = int(o.group(1))
        sample = worst.get("EnergySpectrum")
        bg = worst.get("BackgroundEnergySpectrum")
        if (sample or 0) <= MAX_ORDER and (bg or 0) <= MAX_ORDER:
            continue
        rows.append((rel.split("/")[1], rel,
                     "в git" if rel in tracked_set else "локально",
                     sample, bg))
    return rows


def local_only() -> dict:
    """Что лежит на диске, но исключено .gitignore."""
    counts = {}
    for sub in ("lsrm", "becqmoni", "background"):
        d = REPO / "detectors" / "Gamma-1S" / "reference_spectra" / sub
        counts[sub] = sum(1 for p in d.rglob("*") if p.is_file()) if d.is_dir() else 0
    return counts


def main() -> int:
    per = collect()
    d5 = high_order()
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
    A("`PolynomialEnergyCalibration.CheckCalibration()` отвергает любой порядок выше")
    A("четвёртого:\n")
    A("```csharp")
    A("if (this.polynomialOrder > 4 || this.polynomialOrder < 1) { return false; }")
    A("```\n")
    A("Дальше `DocumentManager.CheckDocument()` пытается спасти положение, но понижает")
    A("степень только за счёт нулевых старших коэффициентов. У настоящего полинома 5-й")
    A("степени их нет, поэтому калибровка заменяется на **дефолтную** — файл открывается,")
    A("шкала теряется, сообщения об ошибке нет.\n")
    A("У этих спектров ЛСРМ записал калибровку 5-й степени. Коэффициенты сохранены как")
    A("в источнике.\n")
    A("| Класс | Файл | Где | Спектр | Фон |")
    A("|---|---|---|---:|---:|")
    for cls, full, status, sdeg, bdeg in d5:
        A(f"| `{cls}` | `{full.split('/')[-1]}` | {status} | {sdeg or '—'} "
          f"| {bdeg or '—'} |")
    A("")
    if any(s == "локально" for _, _, s, _, _ in d5):
        A("«локально» — файл лежит на машине оператора, в репозиторий не публикуется.\n")
    A("**Что делать.** Берите калибровку из `INDEX.json` соответствующего класса —")
    A("поле `energy_cal`, коэффициенты от младшей степени к старшей:\n")
    A("```")
    A("E(n) = a₀ + a₁·n + a₂·n² + a₃·n³ + a₄·n⁴ + a₅·n⁵")
    A("```\n")
    A("Найти такие файлы можно по флагу `becqmoni_reads_as_linear` в том же индексе.\n")
    A("**Фон проверяется отдельно.** BecqMoni валидирует `<EnergySpectrum>` и")
    A("`<BackgroundEnergySpectrum>` как две независимые калибровки, поэтому подшитый фон")
    A("со степенью 5 губил шкалу даже у спектра с линейной калибровкой. Такой фон больше")
    A("не вшивается — остаётся только ссылка `<BackgroundSpectrumFile>`, а признак стоит")
    A("в описи полем `background_dropped_high_order`. Подробности — в")
    A("[`detectors/CHANGELOG.md`](detectors/CHANGELOG.md).\n")
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

    A("## Где лежит фон\n")
    A("Фоновые измерения каждого детектора собраны в отдельную ветку")
    A("`reference_spectra/background/` — `lsrm/` для исходников и `becqmoni/` для XML,")
    A("с сохранением исходной иерархии. Измерения на них ссылаются по имени файла.\n")

    A("## Обозначения\n")
    A("Приборы и источники обозначены внутренними метками вида `SN-nn` и `SRC-nn` —")
    A("одна метка на один прибор или источник, устойчиво между сборками.\n")

    OUT.write_text("\n".join(L), encoding="utf-8")
    print(f"{OUT.name}: {tot_xml} XML, {tot_spe} .spe, {len(per)} классов, "
          f"{len(d5)} спектров со степенью 5")
    bad_bg = [f for _, f, _, _, bdeg in d5 if bdeg and bdeg > MAX_ORDER]
    if bad_bg:
        print(f"  ВНИМАНИЕ: {len(bad_bg)} файлов всё ещё несут вшитый фон степени >4")
        for f in bad_bg[:10]:
            print(f"    {f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
