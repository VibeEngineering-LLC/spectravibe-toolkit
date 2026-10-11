"""Обработчики без интерфейса: кнопки окна вызывают эти функции в фоновом потоке (isgui.widgets.Job).
Сгенерировано gen_code.py (qwen3.8:27b) по specs/backend.md, 2026-10-10; правки Claude помечены."""
import io
import sys
import re
import json
import pathlib
import datetime
import contextlib
import numpy as np


class UserError(Exception):
    """Ошибка, понятная пользователю (показывается текстом, без трассировки)."""


class _Tee(io.TextIOBase):
    def __init__(self, sink):
        self.sink = sink
        self.parts = []

    def write(self, s):
        self.parts.append(s)
        if self.sink is not None:
            self.sink(s)
        return len(s)

    def flush(self):
        pass

    def text(self):
        return "".join(self.parts)


def capture(fn, *args, log=None, **kwargs):
    t = _Tee(log)
    with contextlib.redirect_stdout(t):
        result = fn(*args, **kwargs)
    return (result, t.text())


def friendly(exc) -> str:
    if isinstance(exc, UserError):
        return str(exc)
    if type(exc).__name__ == "Cancelled":   # правка Claude: отмена из окна
        return "Отменено: " + str(exc)
    if isinstance(exc, FileNotFoundError):
        return "Файл не найден: " + str(getattr(exc, "filename", "") or exc)
    if isinstance(exc, PermissionError):
        return "Нет доступа к файлу (открыт в другой программе?): " + str(getattr(exc, "filename", "") or exc)
    if "ParseError" in type(exc).__name__:
        return "Файл не является корректным XML: " + str(exc)
    if isinstance(exc, UnicodeDecodeError):
        return "Файл не в кодировке UTF-8 (не XML BecqMoni?)"
    if isinstance(exc, ValueError):
        return "Ошибка данных: " + str(exc)
    if isinstance(exc, KeyError):
        return "В файле нет ожидаемого поля: " + str(exc)
    return "Внутренняя ошибка: %s: %s" % (type(exc).__name__, exc)


def check_input(path) -> pathlib.Path:
    if not str(path).strip():
        raise UserError("Файл не выбран")
    p = pathlib.Path(path)
    if not p.is_file():
        raise UserError("Файл не найден: %s" % p)
    if p.suffix.lower() != ".xml":
        raise UserError("Нужен файл спектра BecqMoni в формате XML (.xml), выбран: %s" % p.name)
    with open(p, "rb") as f:
        head = f.read(4096)
    if b"<" not in head:
        raise UserError("Файл %s не похож на XML" % p.name)
    return p


def spectrum_info(path) -> dict:
    p = check_input(path)
    try:
        from gamma.calibration.intrinsic_shape.core import load_spectrum
        s = load_spectrum(p)
    except UserError:
        raise
    except Exception as e:
        raise UserError("Не удалось прочитать спектр %s: %s" % (p.name, friendly(e)))
    
    counts = s["counts"]
    if len(counts) == 0 or counts.sum() <= 0:
        raise UserError("В файле %s нет отсчётов спектра" % p.name)
    if s["live"] <= 0:
        raise UserError("В файле %s нет живого времени (live time)" % p.name)
        
    return {
        "path": str(p),
        "name": p.name,
        "live": float(s["live"]),
        "real": float(s["real"]),
        "n_ch": int(len(counts)),
        "total": float(counts.sum()),
        "cal": list(s["cal"]),
        "start": s["start"],
        "counts": counts
    }


def pair_info(bg_path, sample_path) -> dict:
    b = spectrum_info(bg_path)
    s = spectrum_info(sample_path)
    
    if b["n_ch"] != s["n_ch"]:
        raise UserError("Разное число каналов: фон %d, образец %d" % (b["n_ch"], s["n_ch"]))
        
    # правка Claude (#LABR-8): окна пресетов пересчитываются по собственным структурам (автопресет) -- предупреждение
    # только для «образца», который является продолжением того же измерения, что и фон
    from gamma.calibration.intrinsic_shape.pipeline import nested_pair
    warning = nested_pair(bg_path, sample_path)
        
    return {
        "bg": b,
        "sample": s,
        "k": s["live"] / b["live"],
        "warning": warning
    }


INTERP_HELP = {
    "несколько минимумов": "Найдено несколько почти одинаково хороших совмещений. Это нормально; погрешность совмещения взята как их разброс.",
    "на краю диапазона": "Лучшее совмещение лежит на границе перебора: дрейф больше ожидаемого. Результату доверять нельзя без проверки.",
    "chi2/ndf": "Качество совпадения форм заметно отличается от ожидаемого: спектры могут быть от разных приборов или условий.",
    "якорь": "Пик ~37 кэВ тянет решение: проверить окно якоря (каналы пика ~37 кэВ в фоне и образце).",
    "разброс h > 1": "Неопределённость совмещения больше 1 канала: разность в этой области ненадёжна.",
    "профиль": "Поиск минимума неустойчив: профиль нашёл лучшую точку или не параболичен. Нужен разбор специалистом."
}


def explain(item) -> str:
    for key, val in INTERP_HELP.items():
        if key in item:
            return val
    return "Пункт требует разбора специалистом."


def parse_interp(text) -> list:
    lines = text.splitlines()
    idx = None
    for i, line in enumerate(lines):
        if line.strip() == "[ТРЕБУЕТ ТОЛКОВАНИЯ]":
            idx = i
            break
    if idx is None:
        return []
    
    items = []
    for line in lines[idx + 1:]:
        if line.startswith("  "):
            stripped = line.strip()
            if stripped != "(нет пунктов)":
                items.append(stripped)
        else:
            break
    return items


def run_register(bg_path, sample_path, mask=(430, 1700), anchor_bg=28.9, anchor_sample=28.6, anchor_fwhm=3.3,
                 anchor_hw=6, coarse_window=(780, 900), n_ms=20, out_json=None, log=None) -> dict:
    pair_info(bg_path, sample_path)
    
    def _go():
        from gamma.calibration.intrinsic_shape import register as reg
        from gamma.calibration.intrinsic_shape.core import poly_eval
        
        S, B, k, live_s, live_b = reg.load_inputs(bg_path, sample_path)
        print("[IN] k=%.6f live_s=%.3f live_b=%.3f каналов=%d" % (k, live_s, live_b, len(B)))
        
        anchor = reg.anchor37(B, S, x_bg=anchor_bg, x_sm=anchor_sample, w0=anchor_fwhm, hw=anchor_hw, verbose=True)
        
        m = np.zeros(len(B), dtype=bool)
        m[mask[0]:mask[1] + 1] = True
        
        r, out = reg.report(S, B, k, anchor, m, out_json=out_json, profiles=False, n_ms=n_ms, coarse_window=tuple(coarse_window))
        return (anchor, out)

    (anchor, out), text = capture(_go, log=log)
    
    items = parse_interp(text)
    c = out["c"]
    from gamma.calibration.intrinsic_shape.core import poly_eval
    z = (poly_eval(c, anchor[0]) - anchor[1]) / anchor[2]
    
    return {
        "h": [float(v) for v in c],
        "perr": out["perr"],
        "chi2": out["chi2"],
        "ndf": out["ndf"],
        "chi2ndf": out["chi2ndf"],
        "k": out["k"],
        "n_minima": len(out["minima"]),
        "n_close": out["spread"]["n"],
        "grid_on_edge": out["grid_on_edge"],
        "anchor_z": float(z),
        "interp": items,
        "interp_explained": [(it, explain(it)) for it in items],
        "h_json": out,
        "log": text
    }


def default_ecal_preset() -> pathlib.Path:
    from gamma.calibration.intrinsic_shape import ops
    return pathlib.Path(ops.__file__).resolve().parent / "presets" / "labr3_ecal.json"


def run_ecal(bg_path, sample_path, h, preset=None, out_json=None, log=None) -> dict:
    if not isinstance(h, list) or len(h) != 3:
        raise UserError("Сначала выполните «Совместить»: нет коэффициентов h")
        
    path = pathlib.Path(preset) if preset else default_ecal_preset()   # правка Claude: строка -> Path
    cfg = json.loads(path.read_text("utf-8"))
    
    def _go():
        from gamma.calibration.intrinsic_shape import register as reg
        from gamma.calibration.intrinsic_shape import ops
        S, B, k, _, _ = reg.load_inputs(bg_path, sample_path)
        return ops.calibrate(cfg, B, S, list(h))

    res, text = capture(_go, log=log)
    
    checks = [line.strip() for line in text.splitlines() if line.startswith("[CHK]")]
    
    if out_json:
        with open(out_json, "w", encoding="utf-8") as f:
            f.write(json.dumps(res, ensure_ascii=False, indent=1))
            
    result = dict(res)
    result["checks"] = checks
    result["log"] = text
    result["preset"] = str(path)
    return result


def default_fwhm_points() -> list:
    return [(29.20, 2.47), (142.87, 6.14), (211.97, 7.29), (338.14, 8.62), (847.10, 11.63)]


def parse_points(text) -> list:
    tokens = re.split(r"[\s,;]+", text.strip())
    points = []
    for tok in tokens:
        if not tok:
            continue
        parts = tok.split(":")
        if len(parts) != 2:
            raise UserError("Неверная точка ПШПВ «%s»: нужно канал:ширина, например 142.87:6.14" % tok)
        try:
            ch = float(parts[0])
            fwhm = float(parts[1])
            points.append((ch, fwhm))
        except ValueError:
            raise UserError("Неверная точка ПШПВ «%s»: нужно канал:ширина, например 142.87:6.14" % tok)
            
    if len(points) < 3:
        raise UserError("Нужно не меньше 3 точек ПШПВ")
        
    return points


def file_check(path, intended, signed) -> dict:
    text = pathlib.Path(path).read_text("utf-8")
    raw = re.findall(r"<DataPoint>([^<]*)</DataPoint>", text)
    all_int = len(raw) > 0 and all(re.fullmatch(r"-?\d+", v.strip()) for v in raw)
    
    if all_int:
        vals = np.array([int(v) for v in raw])
    else:
        vals = np.zeros(0)
        
    intended_arr = np.asarray(intended)
    n = min(len(vals), len(intended_arr)) - 1
    match = all_int and n > 0 and np.array_equal(vals[:n], intended_arr[:n].astype(np.int64))
    
    has_fwhm = "<SqrtFwhmCalibration>" in text
    min_ok = signed or (all_int and vals.min() >= 0)
    
    return {
        "path": str(path),
        "name": pathlib.Path(path).name,
        "n_points": len(raw),
        "all_int": bool(all_int),
        "match": bool(match),
        "has_fwhm": has_fwhm,
        "sum": int(vals.sum()) if all_int else 0,
        "min": int(vals.min()) if all_int and len(vals) else 0,
        "signed": signed,
        "min_ok": bool(min_ok)
    }


def run_subtract(bg_path, sample_path, h, ecal, fw_points, out_zeroed, out_signed, name=None, date=None, monotone_to=1.5,
                 sum_window=(352, 2046), z_limit=2.0, log=None) -> dict:
    pair = pair_info(bg_path, sample_path)
    
    if not isinstance(h, list) or len(h) != 3:
        raise UserError("Сначала выполните «Совместить»: нет коэффициентов h")
        
    if not ecal:
        raise UserError("Сначала выполните «Шкала»: нет энергетической шкалы")
        
    if not out_zeroed or not out_signed:
        raise UserError("Укажите файлы для записи результата")
        
    p_bg = pathlib.Path(bg_path).resolve()
    p_sm = pathlib.Path(sample_path).resolve()
    p_z = pathlib.Path(out_zeroed).resolve()
    p_s = pathlib.Path(out_signed).resolve()
    
    if p_z == p_s:
        raise UserError("Файлы результата не должны совпадать")
    if p_z == p_bg or p_z == p_sm or p_s == p_bg or p_s == p_sm:
        raise UserError("Нельзя записывать результат поверх исходных файлов")
        
    if date is None:
        date = datetime.date.today().isoformat()
        
    def _go():
        from gamma.calibration.intrinsic_shape import ops
        return ops.make_net(bg_path, sample_path, list(h), list(ecal), list(fw_points), out_zeroed, out_signed, name=name, date=date, monotone_to=monotone_to, sum_window=tuple(sum_window))

    res, text = capture(_go, log=log)
    
    resid = np.asarray(res["resid"], float)
    var = np.asarray(res["var"], float)
    lo, hi = sum_window
    tot = float(resid[lo:hi + 1].sum())
    sig = float(np.sqrt(var[lo:hi + 1].sum()))
    z = tot / sig if sig > 0 else 0.0
    z_ok = abs(z) <= z_limit
    
    files = []
    files.append(file_check(out_zeroed, res["zeroed"], False))
    files.append(file_check(out_signed, res["signed"], True))
    
    from gamma.calibration.intrinsic_shape import fwhm as fwm
    fw_ok, fw_vals, fw_vertex, fw_margin = fwm.check(res["fwhm"], len(resid))
    
    from gamma.calibration.intrinsic_shape import register as reg
    from gamma.calibration.intrinsic_shape.core import poly_eval, rebin_conserve_var
    
    S, B, k, _, _ = reg.load_inputs(bg_path, sample_path)
    edges = poly_eval(list(h), np.arange(len(S) + 1) - 0.5)
    s_on_bg, _ = rebin_conserve_var(np.asarray(S, float), edges, len(B))
    kb = s_on_bg - resid
    energy = poly_eval(list(ecal), np.arange(len(B), dtype=float))
    
    all_ok = z_ok and fw_ok and all(f["match"] and f["all_int"] and f["has_fwhm"] and f["min_ok"] for f in files)
    
    return {
        "k": float(res["k"]),
        "sum": tot,
        "sigma": sig,
        "z": float(z),
        "z_ok": z_ok,
        "z_limit": z_limit,
        "files": files,
        "fwhm": [float(v) for v in res["fwhm"]],
        "fwhm_ok": bool(fw_ok),
        "fwhm_values": fw_vals,
        "fwhm_margin": float(fw_margin),
        "all_ok": bool(all_ok),
        "channels": np.arange(len(B)),
        "energy": energy,
        "sample": s_on_bg,
        "bg_k": kb,
        "resid": resid,
        "zeroed": np.asarray(res["zeroed"]),
        "signed": np.asarray(res["signed"]),
        "n_negative": int((resid < 0).sum()),
        "log": text
    }


def run_fwhm_model(points, n_ch=2048, monotone_to=1.5, log=None) -> dict:
    def _go():
        from gamma.calibration.intrinsic_shape import ops
        return ops.fwhm_model(list(points), n_ch, monotone_to)

    c, text = capture(_go, log=log)
    
    from gamma.calibration.intrinsic_shape import fwhm as fwm
    ok, vals, vertex, margin = fwm.check(c, n_ch)
    
    rows = []
    for ch, f in points:
        model = float(np.sqrt(max(c[0] + c[1] * ch + c[2] * ch * ch, 0.0)))
        rows.append((ch, f, model, (model - f) / f * 100.0))
        
    x = np.arange(n_ch, dtype=float)
    y = np.sqrt(np.maximum(c[0] + c[1] * x + c[2] * x * x, 0.0))
    
    return {
        "coeffs": [float(v) for v in c],
        "ok": bool(ok),
        "values": vals,
        "margin": float(margin),
        "rows": rows,
        "x": x,
        "y": y,
        "log": text
    }


DEVICE_HINT = (
    "Куда вписать: файл %APPDATA%\\BecqMoni\\config\\device\\<имя прибора>.xml, раздел PeakDetectionMethodConfig, поля FWHM_AT_0, Ch_Fwhm, Width_Fwhm.\n"
    "Перед правкой ЗАКРОЙТЕ BecqMoni (иначе он перезапишет файл при выходе). Сделайте копию файла перед изменением.\n"
    "Эти значения действуют для всех спектров прибора без собственного блока ПШПВ."
)


def run_device(points, ch_ref=1309, n_ch=2048, log=None) -> dict:
    def _go():
        from gamma.calibration.intrinsic_shape import ops
        return ops.device_config(list(points), ch_ref, n_ch)

    (F0, ch, W), text = capture(_go, log=log)
    
    text_for_clipboard = "FWHM_AT_0=%s\nCh_Fwhm=%s\nWidth_Fwhm=%s" % (F0, ch, W)
    
    return {
        "FWHM_AT_0": F0,
        "Ch_Fwhm": ch,
        "Width_Fwhm": W,
        "clipboard": text_for_clipboard,
        "hint": DEVICE_HINT,
        "log": text
    }


def run_inventory(folder, recursive=False, pattern="*.xml", log=None) -> dict:
    p = pathlib.Path(folder)
    if not p.is_dir():
        raise UserError("Папка не найдена: %s" % p)
        
    files = sorted(p.rglob(pattern) if recursive else p.glob(pattern))
    if not files:
        raise UserError("В папке нет файлов %s" % pattern)
        
    def _go():
        from gamma.calibration.intrinsic_shape import ops
        return ops.inventory(files)

    (rows, groups), text = capture(_go, log=log)
    
    return {
        "rows": rows,
        "dup_groups": groups,
        "n_files": len(files),
        "n_unique": len({r["md5"] for r in rows}),
        "log": text
    }


def run_decompose(ref_path, unshielded_paths, ecal=None, out_json=None, log=None) -> dict:
    check_input(ref_path)
    for path in unshielded_paths:
        check_input(path)
    if not unshielded_paths:
        raise UserError("Выберите хотя бы один неэкранированный фон")
        
    def _go():
        from gamma.calibration.intrinsic_shape import background as bgm
        return bgm.decompose(ref_path, list(unshielded_paths), e_ref=list(ecal) if ecal else None, out_json=out_json)

    data, text = capture(_go, log=log)
    
    return {
        "summaries": data.get("summaries", {}),
        "log": text
    }


def run_selftest(bg_path, sample_path, alt_bg=(), manual=None, log=None) -> dict:
    pair_info(bg_path, sample_path)
    
    def _go():
        from gamma.calibration.intrinsic_shape import selftest as st
        return st.main(bg_path, sample_path, None, alt_bg=list(alt_bg), manual=manual)

    rc, text = capture(_go, log=log)
    
    itog = [l.strip() for l in text.splitlines() if l.startswith("ИТОГ")]
    
    return {
        "rc": int(rc),
        "green": rc == 0,
        "itog": itog,
        "interp": parse_interp(text),
        "log": text
    }
