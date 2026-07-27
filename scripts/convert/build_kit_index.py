# -*- coding: utf-8 -*-
"""Build the catalogue that ships next to the BecqMoni export.

Reads every kit pair, pulls the source passport (activity, uncertainty,
reference date), the acquisition parameters and the geometry data out of the
LSRM headers, and writes:

  reference_kits_becqmoni/README.md   — human-readable catalogue
  reference_kits_becqmoni/INDEX.json  — the same data, machine-readable
  reference_kits_becqmoni/efficiency/ — the SpectraLine efficiency curves for
                                        the geometries covered by the kit

Usage:
    python scripts/convert/build_kit_index.py
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from gamma.io.lsrm_spe import read_lsrm_spe  # noqa: E402

DETECTOR_DIR = REPO / "detectors" / "Gamma-1S"
KITS = DETECTOR_DIR / "reference_spectra" / "reference_kits"
OUT = DETECTOR_DIR / "reference_spectra" / "reference_kits_becqmoni"
EFF_SRC = DETECTOR_DIR / "efficiency" / "Gamma-1S_NaI_63x63_USB_SN-01"

# Kit geometry -> the LSRM geometry token used in the efficiency file names.
# Kit geometry -> the LSRM geometry token used in the efficiency file names.
GEOMETRY_EFFICIENCY = {
    "Marinelli_1L": "Маринелли",
    "Denta_120mL": "Дента",
    "Petri_60mL": "Петри",
    "Point_25cm": "Точечная-25см",
    "Point_5cm": "Точечная-5см",
}

GEOMETRY_NOTE = {
    "Marinelli_1L": "сосуд Маринелли 1 л, вплотную к детектору",
    "Denta_120mL": ("флакон Дента 120 мл, 0 см. У смеси РИСН-379 в заголовке стоит "
                    "`GEOMETRY=Дента-100` — опечатка оператора, сосуд Дента всегда 120 мл; "
                    "100 г в `RAWMASS` — масса наполнения, а не объём"),
    "Petri_60mL": "чашка Петри 60 мл, 0 см",
    "Point_25cm": "точечный источник, 25 см от торца",
    "Point_5cm": "точечный источник, 5 см от торца",
}

BG_PREFIXES = ("background_", "bg_")
BG_TOKENS = ("фон",)


def is_background(path: Path) -> bool:
    name = path.name.lower()
    return name.startswith(BG_PREFIXES) or any(t in name for t in BG_TOKENS)


def read_efa_field(path: Path, key: str):
    """Pull a single `Key=value` line out of a cp1251 SpectraLine file."""
    try:
        for line in path.read_text(encoding="cp1251", errors="replace").splitlines():
            if line.startswith(key + "="):
                return line.split("=", 1)[1].strip()
    except OSError:
        pass
    return None


def xml_sample_info(xml_path: Path) -> dict:
    root = ET.parse(str(xml_path)).getroot()
    def txt(p):
        el = root.find(p)
        return el.text if el is not None else None
    return {
        "weight_kg": float(txt(".//SampleInfo/Weight") or 0),
        "volume_l": float(txt(".//SampleInfo/Volume") or 0),
        "device_guid": txt(".//DeviceConfigReference/Guid"),
    }


def collect() -> list:
    rows = []
    for leaf in sorted(p for p in KITS.rglob("*") if p.is_dir()):
        spes = sorted(leaf.glob("*.spe"))
        if not spes:
            continue
        samples = [p for p in spes if not is_background(p)]
        bgs = [p for p in spes if is_background(p)]
        if len(samples) != 1 or len(bgs) != 1:
            print(f"  ! пропуск {leaf.relative_to(KITS)}: "
                  f"{len(samples)} образцов / {len(bgs)} фонов")
            continue
        sample, bg = samples[0], bgs[0]
        geometry, nuclide = leaf.relative_to(KITS).parts[0], leaf.relative_to(KITS).parts[1]

        s = read_lsrm_spe(str(sample))
        b = read_lsrm_spe(str(bg))
        xml_path = OUT / geometry / nuclide / (sample.stem + ".xml")
        xi = xml_sample_info(xml_path) if xml_path.exists() else {}

        def num(raw):
            try:
                return float(str(raw).split(";")[0].replace(",", "."))
            except (TypeError, ValueError):
                return None

        rows.append({
            "geometry": geometry,
            "nuclide": nuclide,
            "xml": str(xml_path.relative_to(OUT)).replace("\\", "/"),
            "source_spe": sample.name,
            "background_spe": bg.name,
            "source_id": s.sample_id,
            "passport_raw": s.comments,
            "passport": s.extras.get("lsrm_passport"),
            "measurement_type": s.extras.get("lsrm_type"),
            "detector": s.detector_id,
            "daq_config": s.extras.get("lsrm_config"),
            "lsrm_geometry": read_efa_field(sample, "GEOMETRY"),
            "distance_cm": num(s.extras.get("lsrm_distance")),
            "sample_mass_g": num(s.extras.get("lsrm_samplemass")),
            "sample_volume_ml": num(s.extras.get("lsrm_samplevolume")),
            "measured_at": s.start_datetime.isoformat() if s.start_datetime else None,
            "live_time_s": round(s.live_time, 3),
            "real_time_s": round(s.real_time, 3),
            "dead_time_pct": (round(100.0 * (s.real_time - s.live_time) / s.real_time, 2)
                              if s.real_time else None),
            "channels": int(s.n_channels or len(s.counts)),
            "energy_cal": list(s.energy_cal) if s.energy_cal else None,
            "energy_cal_degree": s.energy_cal_degree,
            # BecqMoni evaluates degree 2..4 only (PolynomialEnergyCalibration.cs);
            # above that it silently reads the calibration as a straight line.
            "becqmoni_reads_as_linear": bool(
                s.energy_cal and len(s.energy_cal) - 1 > 4),
            "lsrm_peaks_found": len(s.extras.get("lsrm_peaks_table") or []),
            "background_live_time_s": round(b.live_time, 3),
            "background_channels": int(b.n_channels or len(b.counts)),
            "becqmoni_weight_kg": xi.get("weight_kg"),
            "becqmoni_volume_l": xi.get("volume_l"),
            "device_guid": xi.get("device_guid"),
        })
    return rows


def copy_efficiency() -> list:
    """Copy the curves for the geometries the kit actually covers."""
    dst = OUT / "efficiency"
    dst.mkdir(parents=True, exist_ok=True)
    wanted = sorted(t for t in set(GEOMETRY_EFFICIENCY.values()) if t)
    copied = []
    for token in wanted:
        for src in sorted(EFF_SRC.glob(f"*{token}.ef*")):
            shutil.copy2(src, dst / src.name)
            copied.append({
                "geometry_token": token,
                "file": src.name,
                "kind": "точки эффективности по опорным источникам" if src.suffix == ".efr"
                        else "описание геометрии и матрицы",
                "volume_ml": read_efa_field(src, "Volume,ml"),
                "density_g_cm3": read_efa_field(src, "Density,g/cm3"),
                "distance_cm": read_efa_field(src, "Distance,cm"),
                "built_at": read_efa_field(src, "Date"),
                "spectraline_version": read_efa_field(src, "EfficiencyVersion"),
            })
    readme = EFF_SRC / "README.md"
    if readme.exists():
        shutil.copy2(readme, dst / "README_source.md")
    return copied


def fmt_activity(row) -> str:
    pp = row.get("passport")
    if not pp:
        return "—"
    parts = []
    for e in pp:
        val, unit = e.get("value"), e.get("unit") or "Бк"
        unc, ref = e.get("uncertainty_pct"), e.get("reference_date")
        chunk = f"{e.get('nuclide','?')} {val:,.0f} {unit}".replace(",", " ")
        if unc is not None:
            chunk += f" ±{unc:g} %"
        if ref:
            chunk += f" на {ref}"
        parts.append(chunk)
    return "; ".join(parts)


def write_readme(rows, eff):
    lines = []
    A = lines.append
    A("# Эталонные спектры Gamma-1S в формате BecqMoni\n")
    A("Набор аттестованных источников, измеренных на спектрометрическом комплексе")
    A("**Гамма-1С** (детектор УДС-ГЦ-63×63-USB № SN-01, NaI(Tl) 63×63 мм, ПО Lsrm SpectraLine).")
    A("Каждый файл самодостаточен: спектр образца и фон той же геометрии лежат в одном XML —")
    A("фон записан как `<BackgroundEnergySpectrum>`, отдельный файл подгружать не нужно.\n")
    A(f"Всего записей: **{len(rows)}** в {len(set(r['geometry'] for r in rows))} геометриях.")
    A("Исходные `.spe` — в соседней папке `../reference_kits/`, конвертация —")
    A("`scripts/convert/reference_kits_to_becqmoni.py` (идемпотентна).\n")

    A("## Как читать таблицы\n")
    A("**Активность** взята из поля `COMMENT` заголовка LSRM — это паспорт источника")
    A("на опорную дату, а не результат измерения. Где указано «Бк/кг», это удельная")
    A("активность наполнителя, где «Бк» — активность точечного источника. Распад к дате")
    A("измерения не пересчитан.\n")
    A("Отдельно про смесь РИСН-379: в `.spe` поле `COMMENT` хранит **одну** строку, поэтому")
    A("в таблице стоит только один нуклид из четырёх (Am-241 для Маринелли, Cs-137 для")
    A("Дента и Петри-60). Полный состав источника в спектрах не записан.\n")
    A("**Единицы.** BecqMoni хранит `Weight` в килограммах и `Volume` в литрах, LSRM `.spe` —")
    A("в граммах и миллилитрах. При конвертации выполняется деление на 1000; колонка")
    A("«масса» ниже — исходная, в граммах.\n")
    A("**Мёртвое время** посчитано как `(real − live) / real`, по временам из заголовка.\n")

    for geom in sorted(set(r["geometry"] for r in rows)):
        grp = [r for r in rows if r["geometry"] == geom]
        A(f"## {geom}\n")
        A(f"{GEOMETRY_NOTE.get(geom, '')}\n")
        A("| Нуклид | Источник | Активность паспорта | Масса, г | Дата измерения | Живое, с | Мёртвое, % | Каналов | Файл |")
        A("|---|---|---|---:|---|---:|---:|---:|---|")
        for r in sorted(grp, key=lambda x: x["nuclide"]):
            mass = f"{r['sample_mass_g']:g}" if r["sample_mass_g"] else "—"
            date = (r["measured_at"] or "")[:10]
            dt = f"{r['dead_time_pct']:g}" if r["dead_time_pct"] is not None else "—"
            A(f"| {r['nuclide']} | {r['source_id']} | {fmt_activity(r)} | {mass} | {date} "
              f"| {r['live_time_s']:.0f} | {dt} | {r['channels']} | `{Path(r['xml']).name}` |")
        A("")
        bgs = sorted({(r["background_spe"], r["background_live_time_s"]) for r in grp})
        for name, lt in bgs:
            A(f"Фон: `{name}`, живое время {lt:.0f} с.")
        A("")

    linear = [r for r in rows if r.get("becqmoni_reads_as_linear")]
    if linear:
        A("## Калибровка, которую BecqMoni прочитает неверно\n")
        A("`PolynomialEnergyCalibration.ChannelToEnergy()` разбирает степени 2, 3 и 4;")
        A("всё выше молча уходит в линейную ветку. У этих файлов LSRM записал полином")
        A("5-й степени — коэффициенты сохранены как в источнике, но шкала в BecqMoni")
        A("будет неверной; для расчётов берите `energy_cal` из `INDEX.json`.\n")
        for r in sorted(linear, key=lambda x: x["xml"]):
            A(f"- `{r['xml']}` — степень {r.get('energy_cal_degree')}")
        A("")

    A("## Кривые эффективности\n")
    A("Лежат отдельно, в подпапке `efficiency/`. Формат Lsrm SpectraLine:")
    A("`.efr` — точки эффективности по опорным источникам, `.efa` — описание геометрии")
    A("и матрицы образца.\n")
    A("**Кривые привязаны к конкретному экземпляру детектора** (№ SN-01) и на другой")
    A("прибор не переносятся даже при одинаковой модели — различаются усиление ФЭУ,")
    A("световыход кристалла и тракт электроники.\n")
    A("| Геометрия набора | Файл | Что внутри | Объём, мл | Плотность, г/см³ | Расстояние, см | Собрана |")
    A("|---|---|---|---|---|---|---|")
    for geom in sorted(GEOMETRY_EFFICIENCY):
        if geom not in {r["geometry"] for r in rows}:
            continue
        token = GEOMETRY_EFFICIENCY[geom]
        for e in [x for x in eff if x["geometry_token"] == token]:
            A(f"| {geom} | `{e['file']}` | {e['kind']} | {e['volume_ml'] or '—'} "
              f"| {e['density_g_cm3'] or '—'} | {e['distance_cm'] or '—'} | {e['built_at'] or '—'} |")
    A("")
    A("## Машиночитаемая опись\n")
    A("`INDEX.json` — те же данные плюс коэффициенты энергетической калибровки,")
    A("реальное время, число пиков, найденных SpectraLine, и GUID прибора, записанный в XML.\n")

    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    if not KITS.is_dir():
        print(f"нет набора: {KITS}", file=sys.stderr)
        return 2
    rows = collect()
    eff = copy_efficiency()
    write_readme(rows, eff)
    (OUT / "INDEX.json").write_text(
        json.dumps({"detector": "Gamma-1S / УДС-ГЦ-63х63-USB №SN-01",
                    "entries": rows, "efficiency": eff},
                   ensure_ascii=False, indent=2),
        encoding="utf-8")
    print(f"записей: {len(rows)}, файлов эффективности: {len(eff)}")
    print(f"  {OUT / 'README.md'}")
    print(f"  {OUT / 'INDEX.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
