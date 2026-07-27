# -*- coding: utf-8 -*-
"""Publish an LSRM spectra tree: strip personal data, convert to BecqMoni.

The LSRM working trees carry an `OPERATOR=` line naming the person who ran the
measurement. Nothing else in the header holds personal data — a scan of all
507 .spe files across the HPGe classes found the names in that field only, and
`COMMENT=` holds source passport activities. So publishing means blanking one
line and leaving the rest byte-for-byte intact.

For every detector class the script:
  1. copies each .spe to reference_spectra/lsrm/<relative path>, with
     OPERATOR blanked (cp1251 preserved);
  2. writes a BecqMoni XML next to it under reference_spectra/becqmoni/,
     embedding the class background when one can be identified;
  3. emits README.md and INDEX.json describing what landed.

Usage:
    python scripts/convert/lsrm_tree_publish.py --list
    python scripts/convert/lsrm_tree_publish.py --class Handy_HPGe
    python scripts/convert/lsrm_tree_publish.py --all
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import uuid
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

sys.path.insert(0, str(REPO / "scripts" / "convert"))

from gamma.io.lsrm_spe import read_lsrm_spe          # noqa: E402
from gamma.io.becqmoni_xml import write_becqmoni_xml  # noqa: E402
from device_guid import device_guid                   # noqa: E402

# Where the untracked LSRM trees live. Kept as a parameter rather than a
# constant so the desktop can point at its own checkout.
DEFAULT_SOURCE = Path(
    r"C:\Users\Дмитрий\Мой диск\Дозиметрия\ИИ\1 Скилы\0_Work\gamma-spectrum-analysis")

# Gamma spectrometers only. The alpha and Si(Li) X-ray classes of the LSRM
# tree (ADA_AlphaDuoSmall, Simple_Alpha, Simple_SiLi) are out of scope here.
CLASSES = [
    "Gamma-1S",                                             # NaI 63x63
    "GP_HPGe20", "NM_HPGe20", "Handy_HPGe", "Simple_HPGe",  # HPGe
    "Handy_LaBr", "Handy_NaI", "Simple_NaI", "Simple_TeCd",  # scintillation / CZT
]

BG_PATTERN = re.compile(r"(?i)bckg|background|фон")

# Классы с замороженной раскладкой: фон остаётся там, где лежит в дереве
# измерений, и в отдельную папку не выносится. Решение оператора.
FROZEN_LAYOUT = {"Gamma-1S"}
OPERATOR_LINE = re.compile(rb"^OPERATOR=.*$", re.MULTILINE)

# Already published as a curated kit — skip so the tree is not duplicated.
SKIP_DIRS = {"reference_kits", "reference_kits_becqmoni", "lsrm", "becqmoni"}

# HPGe spectra are 16384-channel, so one BecqMoni XML is ~590 KB against a
# 33 KB .spe — an 18x blow-up. The bulk research collections below would add
# roughly 250 MB of derived files to a 13 MB repository, so their .spe are
# published but the XML is not: `lsrm_tree_publish.py` regenerates them on
# demand. Pass --xml-all to override.
BULK_DIRS = re.compile(r"(?i)Pu_LNHB|U_LNHB|_LNHB|INR|IPPPE|kpti|SIL_|InfiniteU|Model")

# BecqMoni SampleInfo units — see gamma/io/becqmoni_xml.py BUG-BQ2.
BQ_MIN, BQ_MAX, BQ_DEFAULT = 0.001, 100.0, 1.0

# BecqMoni validates the sample and the embedded background as two independent
# calibrations (`DocumentManager.CheckDocument`, foreground at :132, background
# at :204). `PolynomialEnergyCalibration.CheckCalibration` rejects any order
# above 4, and the recovery branch only downgrades when the top coefficients
# are zero — a genuine 5th-degree fit has none, so the calibration is replaced
# by the default one. A background carrying such a fit therefore destroys the
# energy axis of the whole document even when the sample itself is linear.
# Rather than rewrite the stored calibration, we leave the background out and
# keep only the `<BackgroundSpectrumFile>` reference.
BQ_MAX_POLY_ORDER = 4


def poly_degree(spec) -> int | None:
    return (len(spec.energy_cal) - 1) if spec.energy_cal else None


def to_becqmoni_unit(raw) -> str:
    try:
        value = float(str(raw).split(";")[0].replace(",", "."))
    except (TypeError, ValueError):
        return f"{BQ_DEFAULT:g}"
    if value <= 0:
        return f"{BQ_DEFAULT:g}"
    return f"{min(max(value / 1000.0, BQ_MIN), BQ_MAX):g}"


def scrub(src: Path, dst: Path) -> bool:
    """Copy the .spe with OPERATOR blanked. Returns True if a name was removed."""
    raw = src.read_bytes()
    had_name = False
    for m in OPERATOR_LINE.finditer(raw):
        if m.group(0)[len(b"OPERATOR="):].strip():
            had_name = True
    cleaned = OPERATOR_LINE.sub(b"OPERATOR=", raw)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_bytes(cleaned)
    return had_name


def find_backgrounds(spes: list) -> list:
    return [p for p in spes if BG_PATTERN.search(p.name)]


def pick_background(sample: Path, backgrounds: list):
    """Nearest background: same folder, then nearest ancestor, then the class one."""
    if not backgrounds:
        return None
    same = [b for b in backgrounds if b.parent == sample.parent]
    if same:
        return same[0]
    best, best_depth = None, -1
    for b in backgrounds:
        try:
            common = len(set(b.parts) & set(sample.parts))
        except TypeError:
            common = 0
        if common > best_depth:
            best, best_depth = b, common
    return best


def publish_class(name: str, source_root: Path, xml_all: bool = False) -> dict:
    src_dir = source_root / "detectors" / name
    if not src_dir.is_dir():
        return {"class": name, "error": f"нет папки {src_dir}"}

    spes = sorted(p for p in src_dir.rglob("*.spe")
                  if not SKIP_DIRS & set(p.relative_to(src_dir).parts))
    if not spes:
        return {"class": name, "error": "нет .spe"}

    ref = REPO / "detectors" / name / "reference_spectra"
    out_spe = ref / "lsrm"
    out_xml = ref / "becqmoni"
    # Фоны каждого детектора собираются в свою папку, а не рассыпаны по дереву
    # измерений: их переиспользуют разные геометрии, и искать их по вложенным
    # `Spe/Background` неудобно.
    split_bg = name not in FROZEN_LAYOUT
    out_bg_spe = (ref / "background" / "lsrm") if split_bg else out_spe
    out_bg_xml = (ref / "background" / "becqmoni") if split_bg else out_xml
    backgrounds = find_backgrounds(spes)

    rows, scrubbed, failed, xml_skipped = [], 0, [], 0
    for src in spes:
        # Path relative to the class folder, minus the raw_lsrm/archive prefix.
        rel = src.relative_to(src_dir)
        parts = [p for p in rel.parts if p not in ("raw_lsrm", "reference_spectra")]
        rel = Path(*parts)

        is_bg = bool(BG_PATTERN.search(src.name))
        dst = (out_bg_spe if is_bg else out_spe) / rel
        if scrub(src, dst):
            scrubbed += 1

        # Пути в описи — от reference_spectra, чтобы было видно ветку background/.
        spe_rel = (Path("background") / rel) if (is_bg and split_bg) else rel

        make_xml = xml_all or not BULK_DIRS.search(str(rel))
        if not make_xml:
            xml_skipped += 1
            rows.append({"spe": str(spe_rel).replace("\\", "/"), "xml": None,
                         "is_background": is_bg, "xml_skipped_bulk": True})
            continue

        try:
            spec = read_lsrm_spe(str(dst))
        except Exception as exc:                       # noqa: BLE001
            failed.append({"file": str(rel), "error": f"{type(exc).__name__}: {exc}"})
            continue

        bg_path = None if is_bg else pick_background(src, backgrounds)
        bg_degree = None
        bg_dropped = False
        if bg_path is not None:
            bg_rel = bg_path.relative_to(src_dir)
            bg_rel = Path(*[p for p in bg_rel.parts
                            if p not in ("raw_lsrm", "reference_spectra")])
            bg_dst = out_bg_spe / bg_rel
            if not bg_dst.exists():
                scrub(bg_path, bg_dst)
            try:
                bg_spec = read_lsrm_spe(str(bg_dst))
                bg_degree = poly_degree(bg_spec)
                spec.background_link = bg_path.name
                if bg_degree is not None and bg_degree > BQ_MAX_POLY_ORDER:
                    # See BQ_MAX_POLY_ORDER: embedding it would cost the whole
                    # document its energy axis. Keep the link, drop the data.
                    bg_dropped = True
                else:
                    spec.background_embedded = bg_spec
            except Exception:                          # noqa: BLE001
                spec.background_embedded, bg_path = None, None

        if not spec.detector_id:
            spec.detector_id = name
        if not spec.device_guid:
            # Соль живёт вне репозитория — иначе GUID выдаёт имя прибора,
            # см. scripts/convert/device_guid.py
            spec.device_guid = device_guid(spec.detector_id)
        spec.extras["lsrm_samplemass"] = to_becqmoni_unit(
            spec.extras.get("lsrm_samplemass", ""))
        spec.extras["lsrm_samplevolume"] = to_becqmoni_unit(
            spec.extras.get("lsrm_samplevolume", ""))

        xml_path = (out_bg_xml if is_bg else out_xml) / rel.with_suffix(".xml")
        xml_path.parent.mkdir(parents=True, exist_ok=True)
        write_becqmoni_xml(spec, str(xml_path))

        rows.append({
            "spe": str(spe_rel).replace("\\", "/"),
            "xml": str(spe_rel.with_suffix(".xml")).replace("\\", "/"),
            "is_background": is_bg,
            "sample_id": spec.sample_id,
            "passport_raw": spec.comments,
            "passport": spec.extras.get("lsrm_passport"),
            "detector": spec.detector_id,
            "daq_config": spec.extras.get("lsrm_config"),
            "geometry": read_field(dst, "GEOMETRY"),
            "distance_cm": spec.extras.get("lsrm_distance"),
            "measured_at": spec.start_datetime.isoformat() if spec.start_datetime else None,
            "live_time_s": round(spec.live_time, 3),
            "real_time_s": round(spec.real_time, 3),
            "channels": int(spec.n_channels or len(spec.counts)),
            "energy_cal": list(spec.energy_cal) if spec.energy_cal else None,
            "energy_cal_degree": poly_degree(spec),
            # BecqMoni rejects any order above 4 and falls back to its default
            # calibration — see BQ_MAX_POLY_ORDER.
            "becqmoni_reads_as_linear": bool(
                spec.energy_cal and poly_degree(spec) > BQ_MAX_POLY_ORDER),
            "background_spe": (str(bg_path.name) if bg_path is not None else None),
            "background_cal_degree": bg_degree,
            "background_embedded": bg_path is not None and not bg_dropped,
            "background_dropped_high_order": bg_dropped,
        })

    return {"class": name, "entries": rows, "scrubbed": scrubbed,
            "failed": failed, "xml_skipped_bulk": xml_skipped,
            "backgrounds": [b.name for b in backgrounds]}


def read_field(path: Path, key: str):
    try:
        for line in path.read_text(encoding="cp1251", errors="replace").splitlines():
            if line.startswith(key + "="):
                return line.split("=", 1)[1].strip() or None
    except OSError:
        pass
    return None


def write_class_docs(result: dict):
    name = result["class"]
    base = REPO / "detectors" / name / "reference_spectra"
    rows = result["entries"]
    bulk = [r for r in rows if r.get("xml_skipped_bulk")]
    conv = [r for r in rows if not r.get("xml_skipped_bulk")]
    samples = [r for r in conv if not r["is_background"]]
    bgs = [r for r in conv if r["is_background"]]

    L = []
    A = L.append
    A(f"# {name} — спектры LSRM и экспорт в BecqMoni\n")
    A(f"Спектров всего: **{len(rows)}**. Исходные `.spe` — в `lsrm/`, структура папок")
    A("сохранена. Экспорт в BecqMoni — в `becqmoni/`.\n")
    if bulk:
        A(f"Из них {len(bulk)} спектров опубликованы **только как `.spe`**: это объёмные")
        A("исследовательские коллекции (Pu/U LNHB, INR, IPPPE, kpti). У HPGe 16 384 канала,")
        A("поэтому один BecqMoni-XML весит около 590 КБ против 33 КБ исходника — держать их")
        A("в репозитории неоправданно. Собираются по требованию:\n")
        A("```bash")
        A(f"python scripts/convert/lsrm_tree_publish.py --class {name} --xml-all")
        A("```\n")
    A(f"Сконвертировано: {len(conv)} ({len(samples)} образцов, {len(bgs)} фоновых).\n")
    A("## Очистка персональных данных\n")
    A(f"Из заголовков удалено поле `OPERATOR` — оно было заполнено в {result['scrubbed']} файлах.")
    A("Остальные поля не изменялись: `COMMENT` содержит только паспортные активности")
    A("источников, персональных данных в нём нет.\n")
    linear = [r for r in conv if r.get("becqmoni_reads_as_linear")]
    if linear:
        A("## Калибровка, которую BecqMoni прочитает неверно\n")
        A("`PolynomialEnergyCalibration.ChannelToEnergy()` разбирает только степени 2, 3 и 4;")
        A("всё, что выше, молча уходит в линейную ветку `c[1]·n + c[0]`. У перечисленных")
        A("файлов LSRM записал полином **5-й степени** — коэффициенты сохранены как в")
        A("источнике, но энергетическая шкала в BecqMoni будет неверной. Для расчётов")
        A("берите калибровку из `INDEX.json` (поле `energy_cal`).\n")
        for r in sorted(linear, key=lambda x: x["spe"]):
            A(f"- `{r['spe']}` — степень {r.get('energy_cal_degree')}")
        A("")
    if result["failed"]:
        A("## Файлы, которые не читаются\n")
        for f in result["failed"]:
            A(f"- `{f['file']}` — {f['error']}")
        A("")
    A("## Перечень\n")
    A("| Файл | Источник | Паспорт | Геометрия | Детектор | Дата | Живое, с | Каналов | Фон |")
    A("|---|---|---|---|---|---|---:|---:|---|")
    for r in sorted(samples, key=lambda x: x["spe"]):
        A(f"| `{r['spe']}` | {r['sample_id'] or '—'} | {r['passport_raw'] or '—'} "
          f"| {r['geometry'] or '—'} | {r['detector'] or '—'} "
          f"| {(r['measured_at'] or '')[:10] or '—'} | {r['live_time_s']:.0f} "
          f"| {r['channels']} | {r['background_spe'] or '—'} |")
    A("")
    if bgs:
        A("## Фоновые спектры\n")
        A("| Файл | Живое, с | Каналов | Дата |")
        A("|---|---:|---:|---|")
        for r in sorted(bgs, key=lambda x: x["spe"]):
            A(f"| `{r['spe']}` | {r['live_time_s']:.0f} | {r['channels']} "
              f"| {(r['measured_at'] or '')[:10] or '—'} |")
        A("")
    A("`INDEX.json` рядом — те же данные плюс калибровка, реальное время и разобранный паспорт.\n")

    (base / "README.md").write_text("\n".join(L), encoding="utf-8")
    (base / "INDEX.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    ap.add_argument("--class", dest="cls", choices=CLASSES)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--xml-all", action="store_true",
                    help="конвертировать и объёмные исследовательские коллекции")
    args = ap.parse_args(argv)

    if args.list:
        for c in CLASSES:
            n = len(list((args.source / "detectors" / c).rglob("*.spe")))
            print(f"  {c:<12} {n:4d} .spe")
        return 0

    targets = CLASSES if args.all else ([args.cls] if args.cls else [])
    if not targets:
        ap.error("укажите --class, --all или --list")

    for c in targets:
        res = publish_class(c, args.source, xml_all=args.xml_all)
        if "error" in res:
            print(f"{c}: {res['error']}")
            continue
        write_class_docs(res)
        print(f"{c}: {len(res['entries'])} спектров, "
              f"OPERATOR очищен в {res['scrubbed']}, "
              f"не прочитано {len(res['failed'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
