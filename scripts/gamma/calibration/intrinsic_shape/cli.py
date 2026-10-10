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


def _cmd_register(args):
    from . import register as reg

    S, B, k, live_s, live_b = reg.load_inputs(args.bg, args.sample)
    print(f"[IN] k={k:.6f} live_s={live_s:.3f} live_b={live_b:.3f} каналов={len(B)}")

    anchor = reg.anchor37(
        B, S,
        x_bg=args.anchor_bg,
        x_sm=args.anchor_sample,
        w0=args.anchor_fwhm,
        hw=args.anchor_hw,
        verbose=True
    )

    mask = np.zeros(len(B), dtype=bool)
    lo, hi = args.mask
    mask[lo:hi + 1] = True

    reg.report(
        S, B, k, anchor, mask,
        out_json=args.out,
        profiles=args.profiles,
        n_ms=args.n_ms,
        coarse_window=tuple(args.coarse_window)
    )

    if args.out:
        print(f"[OUT] {args.out}")
    return 0


def _cmd_subtract(args):
    from . import ops

    h = ops.load_h(args.h[0] if len(args.h) == 1 else args.h)
    ecal = _coeffs(args.ecal)
    fw_points = _points(args.fwhm_points)

    ops.make_net(
        args.bg, args.sample, h, ecal, fw_points,
        args.out_zeroed, args.out_signed,
        name=args.name,
        date=args.date,
        monotone_to=_monotone(args.monotone_to),
        sum_window=tuple(args.sum_window)
    )

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
        ops.device_config(points, args.ch_ref)

    return 0


def _cmd_ecal(args):
    from . import register as reg
    from . import ops

    S, B, k, live_s, live_b = reg.load_inputs(args.bg, args.sample)
    cfg_path = pathlib.Path(args.config)
    cfg = json.loads(cfg_path.read_text("utf-8"))

    h = ops.load_h(args.h[0] if len(args.h) == 1 else args.h)
    res = ops.calibrate(cfg, B, S, h)

    if args.out:
        out_path = pathlib.Path(args.out)
        out_path.write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
        print(f"[OUT] {args.out}")

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
    background.decompose(
        args.ref, args.unshielded,
        e_ref=e_ref,
        out_json=args.out
    )
    return 0


def _cmd_selftest(args):
    from . import selftest

    return selftest.main(args.bg, args.sample, args.config, alt_bg=args.alt_bg)


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
    p_reg.add_argument("--mask", nargs=2, type=int, default=[430, 1700], help="каналы функционала регистрации LO HI (шкала фона)")
    p_reg.add_argument("--anchor-bg", type=float, default=28.9, help="канал якоря в фоне (LaBr3: пик ~37 кэВ)")
    p_reg.add_argument("--anchor-sample", type=float, default=28.6, help="канал якоря в образце")
    p_reg.add_argument("--anchor-fwhm", type=float, default=3.3, help="стартовая ПШПВ якоря, кан")
    p_reg.add_argument("--anchor-hw", type=int, default=6, help="полуокно подгонки якоря, кан")
    p_reg.add_argument("--coarse-window", nargs=2, type=int, default=[780, 900], help="окно сильной структуры для грубого наклона LO HI")
    p_reg.add_argument("--n-ms", type=int, default=20, help="число стартов локального спуска")
    p_reg.add_argument("--profiles", action="store_true", help="печатать профили chi2 по c0, c1, c2 (медленно)")
    p_reg.add_argument("--out", help="Файл вывода JSON")
    p_reg.set_defaults(func=_cmd_register)

    # subtract
    p_sub = subparsers.add_parser("subtract", help="Вычитание фона")
    p_sub.add_argument("--bg", required=True, help="Файл фона")
    p_sub.add_argument("--sample", required=True, help="Файл образца")
    p_sub.add_argument("--h", nargs="+", required=True, help="Коэффициенты h (JSON или 3 числа)")
    p_sub.add_argument("--ecal", nargs="+", required=True, help="Коэффициенты ecal (JSON или числа)")
    p_sub.add_argument("--fwhm-points", nargs="+", required=True, help="Точки FWHM 'канал:fwhm'")
    p_sub.add_argument("--out-zeroed", required=True, help="Файл zeroed")
    p_sub.add_argument("--out-signed", required=True, help="Файл signed")
    p_sub.add_argument("--name", help="Имя")
    p_sub.add_argument("--date", default=datetime.date.today().isoformat(), help="Дата")
    p_sub.add_argument("--monotone-to", default="1.5", help="ширина не убывает до канала X*N (число) или c2>=0")
    p_sub.add_argument("--sum-window", nargs=2, type=int, default=[352, 2046], help="Окно суммирования LO HI")
    p_sub.set_defaults(func=_cmd_subtract)

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
    p_ec.add_argument("--config", default=str(pathlib.Path(__file__).resolve().parent / "presets" / "labr3_ecal.json"), help="Конфигурация")
    p_ec.add_argument("--out", help="Файл вывода JSON")
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

    args = parser.parse_args(argv)

    if args.cmd == "selftest" and args.config is None:
        from . import selftest
        args.config = str(selftest.DEFAULT_CFG)

    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
