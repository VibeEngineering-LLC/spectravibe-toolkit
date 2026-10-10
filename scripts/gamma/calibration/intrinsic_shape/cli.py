"""intrinsic-shape -- командная строка калибровки по форме собственных структур (python -m gamma.calibration.intrinsic_shape).
Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-10-10; правки Claude -- только тексты справки.

Подкоманды:
  register      — совмещение шкалы образца со шкалой фона по форме: h, ковариация, список минимумов
  subtract      — разностный спектр, два XML для BecqMoni (обнулённый и со знаком) с блоком ПШПВ
  fwhm          — модель FWHM, измерение, XML, конфигурация устройства
  ecal          — энергетическая шкала по реперам: невязки и LOO
  bg-inventory  — инвентаризация файлов фона
  decompose     — собственная и внешняя части фона, доля K-40 под сумм-пиком
  selftest      — мутационный тест метода на паре «фон + образец»
  sample-cal    — калибровка по пробе: разность на шкале пробы, --use-reper-scale -- шкала по реперам (если рекомендована)
"""

import sys
import json
import pathlib
import datetime
import argparse
import numpy as np


def _coeffs(values):
    if len(values) == 1:
        p = pathlib.Path(values[0])
        if p.is_file():
            data = json.loads(p.read_text("utf-8"))
            if "coeffs" in data:
                return [float(v) for v in data["coeffs"]]
            return [float(v) for v in data["c"]]
    return [float(v) for v in values]


def _points(values):
    res = []
    for v in values:
        parts = v.split(":")
        if len(parts) != 2:
            raise SystemExit(f"Ошибка формата точки: '{v}'. Ожидается 'канал:fwhm'")
        try:
            ch = float(parts[0])
            fw = float(parts[1])
        except ValueError:
            raise SystemExit(f"Ошибка формата точки: '{v}'. Ожидается 'канал:fwhm'")
        res.append((ch, fw))
    return res


def _monotone(value):
    if value == "c2>=0":
        return "c2>=0"
    return float(value)


def _manual(args):
    """Ручная поправка опор автопресета: --anchors-ch A37 SUM [EDGE] (каналы фона)."""
    v = getattr(args, "anchors_ch", None)
    return None if not v else dict(zip(("a37", "sum", "edge"), [float(x) for x in v]))


def _cmd_register(args):
    # правка Claude: автопресет (#LABR-8) + устойчивая регистрация (линии образца); явные ключи перекрывают автопресет
    from . import pipeline as pl
    ov = {"mask": args.mask, "anchor_bg": args.anchor_bg, "anchor_sample": args.anchor_sample, "anchor_fwhm": args.anchor_fwhm,
          "anchor_hw": args.anchor_hw, "coarse_window": args.coarse_window}
    pl.run_register(args.bg, args.sample, auto=not args.no_auto, manual=_manual(args), robust=args.robust, overrides=ov,
                    n_ms=args.n_ms, profiles=args.profiles, out_json=args.out, z_peak=args.z_peak)
    return 0


def _cmd_subtract(args):
    # правка Claude: шкала/ПШПВ из JSON шага ecal (составная шкала ниже 60 кэВ, точки ПШПВ пары), окна -- автопресет
    from . import ops, pipeline as pl
    h_arg = args.h[0] if len(args.h) == 1 else args.h
    h = ops.load_h(h_arg)
    reg_json = h_arg if isinstance(h_arg, str) and pathlib.Path(h_arg).is_file() else None
    e_path = pathlib.Path(args.ecal[0]) if len(args.ecal) == 1 else None
    ecal = pl._ecal_from(str(e_path)) if e_path is not None and e_path.is_file() else pl._ecal_from(_coeffs(args.ecal))
    fw_points = _points(args.fwhm_points) if args.fwhm_points else [tuple(p) for p in ecal.get("fwhm_points", [])]
    if not fw_points:
        print("Ошибка: нет точек ПШПВ: укажите --fwhm-points или JSON шага ecal с полем fwhm_points", file=sys.stderr)
        return 2
    pl.run_subtract(args.bg, args.sample, h, ecal, fw_points, args.out_zeroed, args.out_signed, name=args.name,
                    date=args.date, monotone_to=_monotone(args.monotone_to),
                    sum_window=tuple(args.sum_window) if args.sum_window else None, auto=not args.no_auto,
                    manual=_manual(args), reg_json=reg_json, xml_scale=args.xml_scale, lt_out=args.lt_out)
    print(f"[OUT] {args.out_zeroed}")
    print(f"[OUT] {args.out_signed}")
    return 0


def _cmd_sample_cal(args):
    # правка Claude (VERIFY7): разность на шкале пробы; шкала XML -- калибровка файла пробы, шкала по реперам -- только с --use-reper-scale
    from . import samplecal as sc, samplecal_net as sn
    repers = sc.load_json(args.repers)["repers"] if args.repers else []
    try:
        r = sn.run_sample_cal(args.bg, args.sample, repers, mode=args.mode, use_reper_scale=args.use_reper_scale)
    except ValueError as exc:
        print("Ошибка: %s" % exc, file=sys.stderr)
        return 2
    sn.run_sample_subtract(args.bg, args.sample, r, r["fwhm_points"], args.out_zeroed, args.out_signed, name=args.name,
                           date=args.date, lt_out=args.lt_out)
    for w in r["warnings"]:
        print("[WARN] " + w)
    print(f"[OUT] {args.out_zeroed}")
    print(f"[OUT] {args.out_signed}")
    return 0


def _cmd_fwhm(args):
    from . import ops
    from .core import load_spectrum

    points = _points(args.points or [])

    if args.measure:
        if not args.channels:
            print("Ошибка: для --measure требуется --channels", file=sys.stderr)
            return 2
        counts = load_spectrum(args.measure)["counts"]
        res = ops.measure_fwhm(counts, args.channels, args.fwhm0)
        points += [(c, f) for c, f, e in res]

    if not points:
        print("Ошибка: нет точек FWHM. Укажите --points или --measure.", file=sys.stderr)
        return 2

    try:   # правка Claude: понятная ошибка модели ПШПВ вместо трассировки
        c = ops.fwhm_model(points, args.n_ch, _monotone(args.monotone_to))
    except ValueError as exc:
        print("Ошибка: %s" % exc, file=sys.stderr)
        return 2

    if args.xml_in and args.xml_out:
        if not args.ecal:
            print("Ошибка: для --xml-in/--xml-out требуется --ecal", file=sys.stderr)
            return 2
        ops.fwhm_to_xml(
            args.xml_in, args.xml_out, points,
            _coeffs(args.ecal), args.n_ch, _monotone(args.monotone_to)
        )
        print(f"[OUT] {args.xml_out}")

    if args.device:
        ops.device_config(points, args.ch_ref, args.n_ch)

    return 0


def _cmd_ecal(args):
    # правка Claude: реперы пары (линии образца по библиотеке) + ветвь шкалы ниже 100 кэВ (#LABR-9)
    from . import pipeline as pl
    h = args.h[0] if len(args.h) == 1 else [float(v) for v in args.h]
    pl.run_ecal(args.bg, args.sample, h, preset=args.config, auto=not args.no_auto, manual=_manual(args), low=not args.no_low,
                use37=not args.no_37, e37=args.e37, source=args.source, source_nuclide=args.source_nuclide, out_json=args.out)
    return 0


def _cmd_bg_inventory(args):
    from . import ops

    folder = pathlib.Path(args.folder)
    if args.recursive:
        files = sorted(folder.rglob(args.glob))
    else:
        files = sorted(folder.glob(args.glob))

    if not files:
        print(f"В папке {folder} не найдено файлов по шаблону '{args.glob}'", file=sys.stderr)
        return 1

    ops.inventory(files)
    return 0


def _cmd_decompose(args):
    from . import background

    e_ref = _coeffs(args.ecal) if args.ecal else None
    try:   # правка Claude (#LABR-3): вырожденная шкала -- понятная ошибка
        background.check_ecal(e_ref, background.load_any(args.ref)["counts"], "--ecal") if e_ref is not None else None
    except ValueError as exc:
        print("Ошибка: %s" % exc, file=sys.stderr)
        return 2
    background.decompose(
        args.ref, args.unshielded,
        e_ref=e_ref,
        out_json=args.out
    )
    return 0


def _cmd_selftest(args):
    from . import selftest

    return selftest.main(args.bg, args.sample, args.config, alt_bg=args.alt_bg, auto=not args.no_auto, manual=_manual(args))


def main(argv=None):
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass

    parser = argparse.ArgumentParser(
        prog="intrinsic-shape",
        description="Калибровка по форме собственных структур спектра (LaBr3, CeBr3, LYSO/LuAG)"
    )
    subparsers = parser.add_subparsers(dest="cmd", required=True)

    # register
    p_reg = subparsers.add_parser("register", help="Регистрация образца")
    p_reg.add_argument("--bg", required=True, help="Файл фона")
    p_reg.add_argument("--sample", required=True, help="Файл образца")
    p_reg.add_argument("--mask", nargs=2, type=int, default=None, help="каналы функционала регистрации LO HI (шкала фона); по умолчанию -- автопресет (эталон 430 1700)")
    p_reg.add_argument("--anchor-bg", type=float, default=None, help="канал якоря в фоне (LaBr3: пик ~37 кэВ)")
    p_reg.add_argument("--anchor-sample", type=float, default=None, help="канал якоря в образце")
    p_reg.add_argument("--anchor-fwhm", type=float, default=None, help="стартовая ПШПВ якоря, кан")
    p_reg.add_argument("--anchor-hw", type=int, default=None, help="полуокно подгонки якоря, кан")
    p_reg.add_argument("--coarse-window", nargs=2, type=int, default=None, help="окно сильной структуры для грубого наклона LO HI")
    p_reg.add_argument("--n-ms", type=int, default=20, help="число стартов локального спуска")
    p_reg.add_argument("--profiles", action="store_true", help="печатать профили chi2 по c0, c1, c2 (медленно)")
    p_reg.add_argument("--out", help="Файл вывода JSON")
    p_reg.add_argument("--robust", choices=["auto", "on", "off"], default="auto", help="устойчивая регистрация при линиях образца (auto -- по chi2/ndf и z разности)")
    p_reg.add_argument("--z-peak", type=float, default=4.5, help="порог значимости пика образца в разности для исключения окна")
    p_reg.set_defaults(func=_cmd_register)

    # subtract
    p_sub = subparsers.add_parser("subtract", help="Вычитание фона")
    p_sub.add_argument("--bg", required=True, help="Файл фона")
    p_sub.add_argument("--sample", required=True, help="Файл образца")
    p_sub.add_argument("--h", nargs="+", required=True, help="Коэффициенты h (JSON или 3 числа)")
    p_sub.add_argument("--ecal", nargs="+", required=True, help="Коэффициенты ecal (JSON или числа)")
    p_sub.add_argument("--fwhm-points", nargs="+", default=None, help="Точки FWHM 'канал:fwhm'")
    p_sub.add_argument("--out-zeroed", required=True, help="Файл zeroed")
    p_sub.add_argument("--out-signed", required=True, help="Файл signed")
    p_sub.add_argument("--name", help="Имя")
    p_sub.add_argument("--date", default=datetime.date.today().isoformat(), help="Дата")
    p_sub.add_argument("--monotone-to", default="1.5", help="ширина не убывает до канала X*N (число) или c2>=0")
    p_sub.add_argument("--sum-window", nargs=2, type=int, default=None, help="Окно суммирования LO HI (справочно); по умолчанию -- автопресет")
    p_sub.add_argument("--xml-scale", choices=["composite", "main"], default="composite", help="шкала в XML: составная (ветвь ниже 100 кэВ, полином до 4-й степени) или основная")
    p_sub.add_argument("--lt-out", choices=["sample", "bg"], default="sample", help="живое время выходного файла: образца (по умолчанию) или фона; cps разности от выбора не зависит")
    p_sub.set_defaults(func=_cmd_subtract)

    # sample-cal (правка Claude)
    p_sc = subparsers.add_parser("sample-cal", help="Калибровка по пробе: разность на шкале пробы")
    p_sc.add_argument("--bg", required=True, help="Файл фона")
    p_sc.add_argument("--sample", required=True, help="Файл пробы")
    p_sc.add_argument("--repers", default=None, help="JSON набора реперов пробы (сохранённый на вкладке «Шкала по пробе»)")
    p_sc.add_argument("--mode", choices=["visible", "energy"], default=None, help="перевод фона на пробу (по умолчанию -- автоматически)")
    p_sc.add_argument("--use-reper-scale", action="store_true", help="записать шкалу по реперам (только если рекомендована)")
    p_sc.add_argument("--out-zeroed", required=True, help="Файл zeroed")
    p_sc.add_argument("--out-signed", required=True, help="Файл signed")
    p_sc.add_argument("--name", help="Имя")
    p_sc.add_argument("--date", default=datetime.date.today().isoformat(), help="Дата")
    p_sc.add_argument("--lt-out", choices=["sample", "bg"], default="sample", help="живое время выходного файла")
    p_sc.set_defaults(func=_cmd_sample_cal)

    # fwhm
    p_fw = subparsers.add_parser("fwhm", help="Модель FWHM")
    p_fw.add_argument("--points", nargs="*", help="Точки FWHM 'канал:fwhm'")
    p_fw.add_argument("--measure", help="Спектр для измерения FWHM")
    p_fw.add_argument("--channels", nargs="+", type=float, help="Каналы для измерения")
    p_fw.add_argument("--fwhm0", type=float, default=8.0, help="Начальный FWHM")
    p_fw.add_argument("--n-ch", type=int, default=2048, help="Число каналов")
    p_fw.add_argument("--monotone-to", default="1.5", help="ширина не убывает до канала X*N (число) или c2>=0")
    p_fw.add_argument("--xml-in", help="Входной XML")
    p_fw.add_argument("--xml-out", help="Выходной XML")
    p_fw.add_argument("--ecal", nargs="+", help="Коэффициенты ecal")
    p_fw.add_argument("--device", action="store_true", help="Конфигурация устройства")
    p_fw.add_argument("--ch-ref", type=int, default=1309, help="Референсный канал")
    p_fw.set_defaults(func=_cmd_fwhm)

    # ecal
    p_ec = subparsers.add_parser("ecal", help="Калибровка энергии")
    p_ec.add_argument("--bg", required=True, help="Файл фона")
    p_ec.add_argument("--sample", required=True, help="Файл образца")
    p_ec.add_argument("--h", nargs="+", required=True, help="Коэффициенты h")
    p_ec.add_argument("--config", default=None, help="пресет реперов (по умолчанию presets/labr3_ecal.json, пересчитывается автопресетом)")
    p_ec.add_argument("--out", help="Файл вывода JSON")
    p_ec.add_argument("--no-low", action="store_true", help="не строить ветвь шкалы ниже 100 кэВ")
    p_ec.add_argument("--no-37", action="store_true", help="пик ~37 кэВ не входит в подгонку ветви")
    p_ec.add_argument("--e37", type=float, default=None, help="гамма-эквивалент собственного пика ~37 кэВ (по умолчанию 35.75; литература 35.5-36)")
    p_ec.add_argument("--source", help="спектр внешнего источника того же прибора (реперы ниже 60 кэВ)")
    p_ec.add_argument("--source-nuclide", choices=["Am-241", "Ba-133"], default="Am-241", help="нуклид внешнего источника")
    p_ec.set_defaults(func=_cmd_ecal)

    # bg-inventory
    p_inv = subparsers.add_parser("bg-inventory", help="Инвентаризация фона")
    p_inv.add_argument("folder", help="Папка")
    p_inv.add_argument("--glob", default="*.xml", help="Шаблон файлов")
    p_inv.add_argument("--recursive", action="store_true", help="Рекурсивный поиск")
    p_inv.set_defaults(func=_cmd_bg_inventory)

    # decompose
    p_dec = subparsers.add_parser("decompose", help="Декомпозиция фона")
    p_dec.add_argument("--ref", required=True, help="Экранированный фон")
    p_dec.add_argument("--unshielded", nargs="+", required=True, help="Неэкранированные фоны")
    p_dec.add_argument("--ecal", nargs="+", help="Коэффициенты ecal")
    p_dec.add_argument("--out", help="Файл вывода JSON")
    p_dec.set_defaults(func=_cmd_decompose)

    # selftest
    p_st = subparsers.add_parser("selftest", help="Самопроверка")
    p_st.add_argument("--bg", required=True, help="Файл фона")
    p_st.add_argument("--sample", required=True, help="Файл образца")
    p_st.add_argument("--config", default=None, help="Конфигурация")
    p_st.add_argument("--alt-bg", nargs="*", default=[], help="Дополнительные фоны")
    p_st.set_defaults(func=_cmd_selftest)

    for sp in (p_reg, p_sub, p_ec, p_st):   # правка Claude: автопресет (#LABR-8)
        sp.add_argument("--no-auto", action="store_true", help="не пересчитывать пресеты под прибор (каналы эталона 2048)")
        sp.add_argument("--anchors-ch", nargs="+", type=float, default=None, help="ручная поправка опор фона: A37 SUM [EDGE], каналы")
    args = parser.parse_args(argv)

    if args.cmd == "selftest" and args.config is None:
        from . import selftest
        args.config = str(selftest.DEFAULT_CFG)

    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
