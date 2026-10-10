"""Высокоуровневые операции командной строки intrinsic-shape: шкала по реперам, инвентарь фонов, модель ПШПВ,
разностные спектры. Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-10-10; правки Claude помечены."""

import json
import pathlib
import numpy as np

from .core import load_spectrum, load_pair, poly_eval, poly_deriv, window_quality
from .ecal import fit_poly, leave_one_out, scan
from .register import apply as reg_apply
from .peakfit import fit_peak
from .fwhm import fit_fwhm, fwhm, check, fwhm_block, insert_block, device_defaults
from .netxml import write_xml
from .background import load_any


def load_h(path_or_list):
    """Загружает коэффициенты отображения h из JSON-файла или принимает список из 3 чисел."""
    if isinstance(path_or_list, (str, pathlib.Path)):
        p = pathlib.Path(path_or_list)
        if p.exists():
            with open(p, "r", encoding="utf-8") as f:
                data = json.load(f)
            return [float(v) for v in data["c"]]
    return [float(v) for v in path_or_list]


def calibrate(cfg, bg_counts, sm_counts, h):
    """Калибровка энергетической шкалы по опорным точкам."""
    deg = cfg.get("degree", 2)
    h_err_ch = cfg.get("h_err_ch", 0.15)
    points = cfg["points"]
    
    ch_list = []
    che_list = []
    E_list = []
    Ee_list = []
    point_details = []
    
    for p in points:
        spectrum = p["spectrum"]
        name = p["name"]
        kind = p["kind"]
        x0 = p["x0"]
        w0 = p["w0"]
        hw = p["hw"]
        E = p["E"]
        dE = p["dE"]
        
        counts = bg_counts if spectrum == "bg" else sm_counts
        kind_s = "peak" if kind == "peak" else "edge"
        
        if kind == "given":   # правка Claude: репер, измеренный заранее на шкале фона (линии образца в разности)
            c, e, cs = float(p["ch"]), float(p["ch_err"]), np.array([p["ch"]])
        else:
            c, e, cs = scan(counts, {"x0": x0, "w0": w0, "hw": hw}, kind_s)
        
        if spectrum == "sample":
            # перенос на шкалу фона: c_bg = h(c), ошибка = hypot(h'(c) * e, h_err_ch)
            # Ensure c is float for scalar operations
            c_float = float(c)
            e_float = float(e)
            
            dc_dc = poly_deriv(h, c_float)
            e_bg = np.hypot(dc_dc * e_float, h_err_ch)
            c_bg = poly_eval(h, c_float)
            
            c = c_bg
            e = e_bg
            
        # Print line
        print(f"[A] {spectrum} {name} ch_bg={c:.3f} +-{e:.3f} окна {np.round(cs, 2)}")
        
        ch_list.append(c)
        che_list.append(e)
        E_list.append(E)
        Ee_list.append(dE)
        
        point_details.append({
            "name": name,
            "ch": c,
            "ch_err": e,
            "E": E,
            "E_err": dE
        })
        
    ch = np.array(ch_list)
    che = np.array(che_list)
    E = np.array(E_list)
    Ee = np.array(Ee_list)
    
    fit = fit_poly(ch, che, E, Ee, deg)
    loo = leave_one_out(ch, che, E, Ee, deg)
    
    coeffs = list(fit["coeffs"])
    chi2ndf = fit["chi2ndf"]
    rms = fit["rms"]
    
    print(f"[CAL] deg{deg} coeffs={coeffs} chi2/ndf={chi2ndf:.3f}")
    
    loo_rms = float(np.sqrt(np.nanmean(loo ** 2))) if not np.all(np.isnan(loo)) else float("nan")   # правка Claude: rms, не среднее |LOO|
    
    for i, pd in enumerate(point_details):
        E_fit = poly_eval(coeffs, pd["ch"])
        residual = pd["E"] - E_fit   # правка Claude: невязка = E_ref - E_fit (как в прототипе)
        loo_val = float(loo[i])   # правка Claude: NaN не подменяется нулём
        print(f"[CAL] {pd['name']:15} ch_bg={pd['ch']:.3f} E_ref={pd['E']:.2f} E_fit={E_fit:.3f} невязка={residual:+.3f} LOO={loo_val:+.3f}")
        
        # Update point_details with residual and loo for return
        pd["residual"] = residual
        pd["loo"] = loo_val
        
    print(f"[CAL] rms невязок {rms:.3f} кэВ; LOO rms {loo_rms:.3f} кэВ")
    
    # Checks
    checks = cfg.get("checks", [])
    for item in checks:
        spectrum = item["spectrum"]
        name = item["name"]
        x0 = item["x0"]
        w0 = item["w0"]
        hw = item["hw"]
        
        counts = bg_counts if spectrum == "bg" else sm_counts
        c, e, _ = scan(counts, {"x0": x0, "w0": w0, "hw": hw}, "peak")
        
        if spectrum == "sample":
            c = poly_eval(h, float(c))
            
        E_val = poly_eval(coeffs, float(c))
        print(f"[CHK] {name} ch_bg={c:.3f} -> E={E_val:.2f} кэВ (в шкалу не входил)")
        
    return {
        "coeffs": coeffs,
        "chi2ndf": float(chi2ndf),
        "points": point_details,
        "rms": float(rms),
        "loo_rms": float(loo_rms)
    }


def inventory(files, windows=None):
    """Инвентаризация файлов фона: проверка дублей и статистика по окнам."""
    if windows is None:
        windows = {
            "BaK": (22, 37),
            "lo": (40, 420),
            "b430": (430, 1700),
            "sum": (830, 860),
            "a3": (1585, 1635)
        }
        
    rows = []
    md5_map = {}
    
    for f in files:
        data = load_any(f)
        name = data["name"]
        n_raw = data["n_raw"]
        live = data["live"]
        start = data["start"]
        md5 = data["md5"]
        counts = data["counts"]
        
        total_rate = counts.sum() / live if live > 0 else 0.0
        
        rates = {}
        for wname, (lo, hi) in windows.items():
            w_rate = counts[lo:hi].sum() / live if live > 0 else 0.0
            rates[wname] = w_rate
            
        row = {
            "path": str(f),
            "n_raw": n_raw,
            "live": live,
            "start": start,
            "md5": md5,
            "rate": total_rate,
            "rates": rates
        }
        rows.append(row)
        
        if md5 not in md5_map:
            md5_map[md5] = []
        md5_map[md5].append(str(f))
        
        # Print line
        rate_str = " ".join([f"{wname}={rates[wname]:.5f}" for wname in windows])
        print(f"{name} n_raw={n_raw} live={live:.3f} start={start} rate={total_rate:.5f} {rate_str} md5={md5}")
        
    dup_groups = [v for v in md5_map.values() if len(v) > 1]
    n_files = len(files)
    n_unique = len(md5_map)
    n_dups = len(dup_groups)
    
    print(f"[EQ] файлов {n_files}, различных измерений {n_unique}, групп дублей {n_dups}")
    
    for group in dup_groups:
        eq_str = " = ".join(group)
        print(f"[EQ] одинаковый счёт: {eq_str}")
        
    return rows, dup_groups


def measure_fwhm(counts, channels, fwhm0, max_rel_err=0.2, width_range=(0.4, 2.5)):
    """Измерение ПШПВ в указанных каналах. Правка Claude: точка отбраковывается (с причиной), если
    σ(ПШПВ)/ПШПВ > max_rel_err или ПШПВ вне width_range * стартовой ширины (в разностном спектре -- нет пика)."""
    results = []
    for x in channels:
        w0 = fwhm0(x) if callable(fwhm0) else fwhm0
        half_window = max(6, int(round(2 * w0)))
        
        res = fit_peak(counts, x, w0, half_window=half_window)
        
        if res is None:
            print(f"[FWHM] ch {x}: подгонка не сошлась")
        else:
            centroid = res["centroid"]
            fwhm_val = res["fwhm"]
            fwhm_err = res["fwhm_err"]
            rel = fwhm_err / fwhm_val if fwhm_val > 0 else float("inf")
            lo_w, hi_w = width_range[0] * w0, width_range[1] * w0
            if not (rel <= max_rel_err and lo_w <= fwhm_val <= hi_w):
                print(f"[FWHM] ch {x}: ОТБРОШЕНА -- центроид {centroid:.2f} ПШПВ {fwhm_val:.2f} +- {fwhm_err:.2f} "
                      f"(отн. σ {rel:.2f}, порог {max_rel_err}; допустимая ширина {lo_w:.2f}..{hi_w:.2f})")
                continue
            results.append((centroid, fwhm_val, fwhm_err))
            print(f"[FWHM] ch {x}: центроид {centroid:.2f} ПШПВ {fwhm_val:.2f} +- {fwhm_err:.2f}")
            
    return results


def fwhm_model(points, n_ch, monotone_to=1.5):
    """Построение модели ПШПВ и проверка."""
    x = [p[0] for p in points]
    f = [p[1] for p in points]
    if len(points) < 3:   # правка Claude (#LABR-3): модель с 3 коэффициентами по 1-2 точкам данными не определена
        print("[FWHM] ПРЕДУПРЕЖДЕНИЕ: точек ПШПВ %d < 3 -- коэффициенты c0, c1, c2 данными не определены "
              "(подобраны ограничениями модели); ширина вдали от точек ненадёжна, добавьте точки" % len(points))
    
    c = fit_fwhm(n_ch, x, f, monotone_to=monotone_to)
    
    ok, vals, vertex, margin = check(c, n_ch)
    
    print(f"[FWHM] модель c0,c1,c2 = {c}")
    
    for ch, fwhm_meas in points:
        fwhm_model_val = fwhm(c, ch)
        residual = fwhm_model_val - fwhm_meas
        pct = (residual / fwhm_meas * 100) if fwhm_meas != 0 else 0.0
        print(f"[FWHM] ch {ch:.2f}: измерено {fwhm_meas:.2f}, модель {fwhm_model_val:.2f}, невязка {residual:+.2f} кан ({pct:+.1f} %)")
        
    print(f"[FWHM] CheckCalibration: годен={ok} ПШПВ(0,500,1309,N-1)={vals} вершина параболы={vertex:.1f} запас за N-1={margin:.1f} кан")
    
    if not ok:
        raise ValueError("Модель ПШПВ не прошла проверку")
        
    return c


def fwhm_to_xml(xml_in, xml_out, points, ecal, n_ch, monotone_to=1.5):
    """Вставка блока калибровки ПШПВ в XML-файл."""
    c = fwhm_model(points, n_ch, monotone_to)
    
    with open(xml_in, "r", encoding="utf-8", newline="") as f:
        text = f.read()
        
    nl = "\r\n" if "\r\n" in text else "\n"
    
    block = fwhm_block(nl, c, points, ecal)
    new_text = insert_block(text, block)
    
    with open(xml_out, "w", encoding="utf-8", newline="") as f:
        f.write(new_text)
        
    return c


def device_config(points, ch_ref=1309, n_ch=2048):
    """Расчёт параметров устройства для конфигурации ПШПВ."""
    x = [p[0] for p in points]
    f = [p[1] for p in points]
    
    F0, ch, W, coefs = device_defaults(x, f, ch_ref)
    
    print(f"FWHM_AT_0={F0} Ch_Fwhm={ch} Width_Fwhm={W}")
    
    for ch_val, fwhm_meas in points:
        fwhm_model_val = np.sqrt(coefs[0] + coefs[1] * ch_val)
        residual = fwhm_model_val - fwhm_meas
        pct = (residual / fwhm_meas * 100) if fwhm_meas != 0 else 0.0
        print(f"[FWHM] ch {ch_val:.2f}: измерено {fwhm_meas:.2f}, модель {fwhm_model_val:.2f}, невязка {residual:+.2f} кан ({pct:+.1f} %)")
        
    ok, vals, vertex, margin = check(coefs, n_ch)   # правка Claude: число каналов прибора
    print(f"[FWHM] CheckCalibration: годен={ok} ПШПВ(0,500,1309,N-1)={vals} вершина параболы={vertex:.1f} запас за N-1={margin:.1f} кан")
    
    return F0, ch, W


def make_net(bg_path, sample_path, h, ecal, fw_points, out_zeroed, out_signed, name=None, date="", monotone_to=1.5, sum_window=(352, 2046)):
    """Построение чистого спектра (вычитание фона) и запись в XML."""
    bg, sm = load_pair(bg_path, sample_path)
    B = bg["counts"]
    S = sm["counts"]
    k = sm["live"] / bg["live"]
    
    resid, var = reg_apply(h, S, B, k)
    
    net_zeroed = np.maximum(0.0, np.round(resid))
    net_signed = np.round(resid)
    
    fwc = fwhm_model(fw_points, len(B), monotone_to)
    
    if name is None:
        name = pathlib.Path(sm["name"]).stem + " - фон (рег. по форме)"
        
    write_xml(sample_path, out_zeroed, net_zeroed, False, ecal, h, k, fwc, fw_points, name, date)
    write_xml(sample_path, out_signed, net_signed, True, ecal, h, k, fwc, fw_points, name, date)
    
    min_val = np.min(resid)
    neg_channels = np.sum(resid < 0)
    clamped_count = int(np.sum(net_signed < 0))   # правка Claude: каналы, которые кламп обнулил (целое < 0)
    
    print(f"[NET] k={k:.6f} min={min_val:.2f} отриц. каналов={int(neg_channels)} изменено клампом={clamped_count}")
    
    lo, hi = sum_window
    S_ = resid[lo:hi+1].sum()
    sg = np.sqrt(var[lo:hi+1].sum())
    z = S_ / sg if sg > 0 else 0.0
    print(f"[ACC] сумма net кан {lo}-{hi} = {S_:.1f} +- {sg:.1f} (z={z:+.2f})")
    
    # Acceptance checks
    for out_path, intended_data in [(out_zeroed, net_zeroed), (out_signed, net_signed)]:
        re_read = load_spectrum(out_path, n_ch=len(B))
        counts_re = re_read["counts"]
        
        # Ignore last channel
        n_check = len(intended_data) - 1
        match = np.array_equal(counts_re[:n_check], intended_data[:n_check])
        
        n_channels = len(counts_re)
        sum_val = counts_re.sum()
        live_val = re_read["live"]
        
        print(f"[ACC] {pathlib.Path(out_path).name}: каналов={n_channels} сумма={sum_val:.0f} live={live_val:.3f} совпадает={match}")
        
    return {
        "k": k,
        "resid": resid,
        "var": var,
        "fwhm": fwc,
        "zeroed": net_zeroed,
        "signed": net_signed
    }
