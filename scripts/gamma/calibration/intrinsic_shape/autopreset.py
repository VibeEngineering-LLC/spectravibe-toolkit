"""
Модуль автоматического пересчёта пресетов калибровки для нового спектрометра.
Определяет опорные структуры (пик 37 кэВ, сумма, край бета-спектра) в фоне,
строит полиномиальное отображение каналов эталона на каналы прибора и
пересчитывает все канальные параметры пресетов.
Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-10-10
"""

import copy
import json
import pathlib
import numpy as np

from .core import poly_eval, poly_deriv, poly_inverse, smooth_gauss
from .peakfit import fit_peak
from .ecal import scan

HERE = pathlib.Path(__file__).resolve().parent
REF_PATH = HERE / "presets" / "labr3_reference.json"


def load_ref(path=None) -> dict:
    """Загружает эталонный файл пресетов."""
    p = path or REF_PATH
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)


def _peak_ok(r, fwhm_exp) -> (bool, str):
    """Проверяет результат подгонки пика."""
    if r is None:
        return False, "подгонка не сошлась"
    
    centroid = r.get("centroid")
    centroid_err = r.get("centroid_err")
    fwhm = r.get("fwhm")
    signif = r.get("signif")
    
    if not all(np.isfinite(v) for v in [centroid, centroid_err, fwhm]):
        return False, "нечисловой результат"
    
    if signif < 10:
        return False, "значимость %.1f < 10" % signif
    
    if centroid_err > 0.1 * fwhm:
        return False, "ошибка центроида велика"
    
    if not (0.5 * fwhm_exp <= fwhm <= 2.0 * fwhm_exp):
        return False, "ширина %.2f вне 0.5..2 от ожидаемой %.2f" % (fwhm, fwhm_exp)
    
    return True, ""


def find_a37(counts, ref) -> dict:
    """Ищет пик 37 кэВ в начале шкалы."""
    n = len(counts)
    r0 = n / ref["n_ch"]
    
    nz = np.where(counts > 0)[0]
    if len(nz) == 0:
        return {"ok": False, "reason": "пустой спектр"}
    
    i0 = nz[0]
    lo = int(i0 + max(2, round(2 * r0)))
    hi = max(lo + 5, n // 16)
    
    sm = smooth_gauss(counts.astype(float), max(1.0, r0))
    xg = lo + int(np.argmax(sm[lo:hi]))
    
    r_est = xg / ref["anchors"]["a37"]["x"]
    fwhm_exp = ref["anchors"]["a37"]["fwhm"] * r_est
    hw = max(6, int(round(ref["anchors"]["a37"]["hw"] * r_est)))
    
    res = fit_peak(counts, float(xg), fwhm_exp, half_window=hw)
    ok, why = _peak_ok(res, fwhm_exp)
    
    if res is not None:
        x = float(res["centroid"])
        err = float(res["centroid_err"])
        fwhm = float(res["fwhm"])
        signif = float(res["signif"])
    else:
        x = float(xg)
        err = float("nan")
        fwhm = float(fwhm_exp)
        signif = 0.0
    
    return {
        "ok": ok,
        "reason": why,
        "x": x,
        "err": err,
        "fwhm": fwhm,
        "signif": signif,
        "x_guess": float(xg),
        "manual": False
    }


def find_sum(counts, a37_x, ref, cal=None) -> dict:
    """Ищет суммарный пик La-138 (~1471 кэВ)."""
    A = ref["anchors"]
    n = len(counts)
    
    x_ratio = a37_x * A["sum"]["x"] / A["a37"]["x"]
    
    x_cal = None
    if cal is not None and len(cal) >= 2:
        try:
            x_cal = float(poly_inverse(list(cal), 1471.0, hi=n + 200))
        except Exception:
            x_cal = None
    
    if x_cal is not None and (not np.isfinite(x_cal) or not (0 < x_cal < n)):
        x_cal = None
    
    note = ""
    if x_cal is not None and abs(x_cal / x_ratio - 1) < 0.10:
        x_pred = x_cal
    else:
        x_pred = x_ratio
        if x_cal is not None:
            note = "калибровка файла (%.1f) расходится с отношением к пику 37 кэВ (%.1f) больше 10 %%; взято отношение" % (x_cal, x_ratio)
    
    r_est = x_pred / A["sum"]["x"]
    lo = max(0, int(x_pred * 0.94))
    hi = min(n - 1, int(x_pred * 1.06) + 1)
    
    sm = smooth_gauss(counts.astype(float), max(1.0, r_est))
    xg = lo + int(np.argmax(sm[lo:hi + 1]))
    
    fwhm_exp = A["sum"]["fwhm"] * r_est
    hw = max(7, int(round(A["sum"]["hw"] * r_est)))
    
    res = fit_peak(counts, float(xg), fwhm_exp, half_window=hw)
    ok, why = _peak_ok(res, fwhm_exp)
    
    if res is not None:
        x = float(res["centroid"])
        err = float(res["centroid_err"])
        fwhm = float(res["fwhm"])
        signif = float(res["signif"])
    else:
        x = float(xg)
        err = float("nan")
        fwhm = float(fwhm_exp)
        signif = 0.0
    
    return {
        "ok": ok,
        "reason": why,
        "x": x,
        "err": err,
        "fwhm": fwhm,
        "signif": signif,
        "x_guess": float(xg),
        "manual": False,
        "x_pred": float(x_pred),
        "x_cal": float(x_cal) if x_cal is not None else None,
        "note": note
    }


def find_edge(counts, m, ref) -> dict:
    """Ищет край бета-спектра La-138 (788.7 кэВ)."""
    E = ref["anchors"]["edge"]
    x0 = float(poly_eval(m, E["x"]))
    s = float(poly_deriv(m, E["x"]))
    
    w0 = E["sigma"] * s
    hw = max(8, int(round(E["hw"] * s)))
    
    try:
        c, e, cs = scan(counts, {"x0": x0, "w0": w0, "hw": hw}, "edge")
    except Exception:
        c, e = float("nan"), float("nan")
    
    c = float(c) if np.isfinite(c) else float("nan")
    e = float(e) if np.isfinite(e) else float("nan")
    
    ok = False
    reason = ""
    
    if not np.isfinite(c) or not np.isfinite(e):
        reason = "подгонка края не сошлась"
    elif abs(c - x0) >= 0.25 * hw:
        reason = "край далеко от ожидаемого (%.1f против %.1f)" % (c, x0)
    elif e >= 0.2 * hw:
        reason = "ошибка края велика"
    else:
        ok = True
    
    return {
        "ok": ok,
        "reason": reason,
        "x": c if np.isfinite(c) else x0,
        "err": e if np.isfinite(e) else float("nan"),
        "sigma": float(w0),
        "x_pred": float(x0),
        "manual": False
    }


def build_map(pairs, ref, r_hint) -> (list, dict):
    """Строит полиномиальное отображение по опорным точкам."""
    if len(pairs) < 2:
        raise ValueError("для отображения нужно не меньше двух опор")
    
    x_ref = np.array([p[0] for p in pairs], dtype=float)
    x_dev = np.array([p[1] for p in pairs], dtype=float)
    errs = np.array([p[2] for p in pairs], dtype=float)
    
    # Заменяем NaN ошибки на 1.0 * r_hint
    errs = np.where(np.isfinite(errs), errs, 1.0 * r_hint)
    # Ограничиваем ошибку снизу
    errs = np.maximum(errs, 0.05 * r_hint)
    
    w = 1.0 / (errs ** 2)
    
    # Линейная подгонка: x_dev = a + b * x_ref
    # Матрица дизайна [1, x_ref]
    A_lin = np.column_stack([np.ones_like(x_ref), x_ref])
    # Взвешенные нормальные уравнения: (A^T W A) theta = A^T W y
    W = np.diag(w)
    ATA = A_lin.T @ W @ A_lin
    ATy = A_lin.T @ W @ x_dev
    
    try:
        theta_lin = np.linalg.solve(ATA, ATy)
        a, b = theta_lin[0], theta_lin[1]
    except np.linalg.LinAlgError:
        # Fallback: обычная подгонка
        a, b = np.polyfit(x_ref, x_dev, 1)[::-1]
    
    resid_lin = x_dev - (a + b * x_ref)
    
    kind = "линейное"
    coeffs = [float(a), float(b), 0.0]
    
    if len(pairs) >= 3 and np.max(np.abs(resid_lin)) > ref["quad_tol_ch"] * r_hint:
        # Квадратичная подгонка: x_dev = a + b*x + c*x^2
        A_quad = np.column_stack([np.ones_like(x_ref), x_ref, x_ref ** 2])
        ATA_q = A_quad.T @ W @ A_quad
        ATy_q = A_quad.T @ W @ x_dev
        
        try:
            theta_quad = np.linalg.solve(ATA_q, ATy_q)
            a_q, b_q, c_q = theta_quad
            coeffs = [float(a_q), float(b_q), float(c_q)]
            kind = "квадратичное"
        except np.linalg.LinAlgError:
            pass
    
    resid = x_dev - np.array([poly_eval(coeffs, x) for x in x_ref])
    scale = float(poly_deriv(coeffs, ref["n_ch"] / 2))
    
    info = {
        "kind": kind,
        "resid_lin": [float(v) for v in resid_lin],
        "resid": [float(v) for v in resid],
        "scale": scale
    }
    
    return coeffs, info


def map_x(m, x):
    """Отображает координату x через полином m."""
    return float(poly_eval(m, x))


def map_scale(m, x):
    """Возвращает производную отображения в точке x."""
    return float(poly_deriv(m, x))


def map_win(m, lo, hi, n):
    """Отображает окно [lo, hi] в координаты прибора."""
    lo_new = max(0, int(round(map_x(m, lo))))
    hi_new = min(n - 2, int(round(map_x(m, hi))))
    return [lo_new, hi_new]


def is_identity(m, n, ref):
    """Проверяет, является ли отображение тождественным."""
    if n != ref["n_ch"]:
        return False
    
    x = np.linspace(0, n - 1, 200)
    y = np.array([poly_eval(m, xi) for xi in x])
    max_diff = np.max(np.abs(y - x))
    
    return max_diff < ref["identity_tol_ch"]


def verify(counts, m, ref, found) -> list:
    """Проверяет опорные точки после отображения."""
    results = []
    
    for name in ("a37", "sum", "edge"):
        if name not in found:
            continue
        
        anchor = ref["anchors"][name]
        x_map = map_x(m, anchor["x"])
        s = map_scale(m, anchor["x"])
        
        if name == "edge":
            w0 = anchor["sigma"] * s
            hw = max(8, int(round(anchor["hw"] * s)))
            try:
                c, e, cs = scan(counts, {"x0": x_map, "w0": w0, "hw": hw}, "edge")
                x_meas = float(c) if np.isfinite(c) else float("nan")
            except Exception:
                x_meas = float("nan")
            tol = ref["verify_tol_edge_sigma"] * anchor["sigma"] * s
        else:
            fw = anchor["fwhm"] * s
            hw = max(6, int(round(anchor["hw"] * s)))
            res = fit_peak(counts, x_map, fw, half_window=hw)
            x_meas = float(res["centroid"]) if res is not None else float("nan")
            tol = ref["verify_tol_fwhm"] * fw
        
        diff = x_meas - x_map
        ok = np.isfinite(diff) and abs(diff) <= tol
        
        results.append({
            "name": name,
            "x_map": float(x_map),
            "x_meas": float(x_meas),
            "diff": float(diff),
            "tol": float(tol),
            "ok": bool(ok)
        })
    
    return results


def map_register_kw(m, ref, a37_bg, a37_sm, n):
    """Пересчитывает параметры регистра."""
    R = ref["register"]
    s37 = map_scale(m, R["anchor_bg"])
    b = map_scale(m, ref["n_ch"] / 2)
    
    return {
        "mask": map_win(m, *R["mask"], n),
        "anchor_bg": float(a37_bg),
        "anchor_sample": float(a37_sm),
        "anchor_fwhm": float(R["anchor_fwhm"] * s37),
        "anchor_hw": max(6, int(round(R["anchor_hw"] * s37))),
        "coarse_window": map_win(m, *R["coarse_window"], n),
        "c2_lo": float(R["c2_lo"] / b),
        "c2_hi": float(R["c2_hi"] / b),
        "c2_step": float(R["c2_step"] / b),
        "test_ch": [round(map_x(m, v), 1) for v in R["test_ch"]]
    }


def map_sum_window(m, ref, n):
    """Пересчитывает окно суммы."""
    lo, hi = ref["sum_window"]
    if hi >= ref["n_ch"] - 2:
        return [max(0, int(round(map_x(m, lo)))), n - 2]
    else:
        return map_win(m, lo, hi, n)


def map_fwhm_points(m, pts):
    """Пересчитывает точки FWHM."""
    return [(round(map_x(m, x), 2), round(f * map_scale(m, x), 3)) for x, f in pts]


def map_ecal_cfg(m, cfg):
    """Пересчитывает конфигурацию ecal."""
    cfg_copy = copy.deepcopy(cfg)
    
    for p in cfg_copy["points"]:
        s = map_scale(m, p["x0"])
        p["w0"] = p["w0"] * s
        p["hw"] = max(5, int(round(p["hw"] * s)))
        p["x0"] = map_x(m, p["x0"])
    
    if "checks" in cfg_copy:
        for p in cfg_copy["checks"]:
            s = map_scale(m, p["x0"])
            p["w0"] = p["w0"] * s
            p["hw"] = max(5, int(round(p["hw"] * s)))
            p["x0"] = map_x(m, p["x0"])
    
    cfg_copy["_mapped"] = True
    cfg_copy["_map"] = list(m)
    
    return cfg_copy


def map_selftest_cfg(m, cfg, ref, a37_bg, a37_sm, n):
    """Пересчитывает конфигурацию selftest."""
    cfg_copy = copy.deepcopy(cfg)
    b = map_scale(m, ref["n_ch"] / 2)
    
    # Окна
    if "windows" in cfg_copy:
        for name, win in cfg_copy["windows"].items():
            cfg_copy["windows"][name] = map_win(m, win[0], win[1], n)
    
    # Полная маска
    if "full_mask" in cfg_copy:
        cfg_copy["full_mask"] = map_win(m, cfg_copy["full_mask"][0], cfg_copy["full_mask"][1], n)
    
    # Интегралы
    if "integrals" in cfg_copy:
        new_integrals = []
        for pair in cfg_copy["integrals"]:
            lo, hi = pair
            if hi >= ref["n_ch"] - 2:
                new_integrals.append([max(0, int(round(map_x(m, lo)))), n - 2])   # правка Claude: целый канал
            else:
                new_integrals.append(map_win(m, lo, hi, n))
        cfg_copy["integrals"] = new_integrals
    
    # Пивот
    if "pivot_joint" in cfg_copy:
        cfg_copy["pivot_joint"] = map_x(m, cfg_copy["pivot_joint"])
    
    # Сдвиги
    for key in ("shift_scan", "shift_scan_wide"):
        if key in cfg_copy:
            cfg_copy[key] = {kk: vv * b for kk, vv in cfg_copy[key].items()}   # правка Claude: скан -- словарь lo/hi/step
    
    # Сетка
    if "joint_grid" in cfg_copy:
        jg = cfg_copy["joint_grid"]
        jg["d_lo"] = jg["d_lo"] * b
        jg["d_hi"] = jg["d_hi"] * b
        jg["d_step"] = jg["d_step"] * b
    
    # Якорь
    if "anchor" in cfg_copy:
        anc = cfg_copy["anchor"]
        s_anc = map_scale(m, anc["x_bg"])
        anc["x_bg"] = float(a37_bg)
        anc["x_sm"] = float(a37_sm)
        anc["w0"] = anc["w0"] * s_anc
        anc["hw"] = max(6, int(round(anc["hw"] * s_anc)))
    
    # Мутации
    if "mutations" in cfg_copy:
        for mut in cfg_copy["mutations"]:
            if mut.get("kind") == "shift":
                mut["value"] = mut["value"] * b
            elif mut.get("kind") == "stretch":
                mut["pivot"] = map_x(m, mut["pivot"])
    
    # Регистрация
    reg_kw = map_register_kw(m, ref, a37_bg, a37_sm, n)
    cfg_copy["registration"] = dict(cfg_copy.get("registration", {}))
    cfg_copy["registration"].update({
        "coarse_window": reg_kw["coarse_window"],
        "c2_lo": reg_kw["c2_lo"],
        "c2_hi": reg_kw["c2_hi"],
        "c2_step": reg_kw["c2_step"],
        "test_ch": reg_kw["test_ch"]
    })
    
    cfg_copy["_mapped"] = True
    cfg_copy["_map"] = list(m)
    
    return cfg_copy


def autopreset(B, S, cal_bg=None, manual=None, ref=None, verbose=True) -> dict:
    """Основная функция автоматического пересчёта пресетов."""
    ref = ref or load_ref()
    B = np.asarray(B, dtype=float)
    S = np.asarray(S, dtype=float)
    n = len(B)
    r0 = n / ref["n_ch"]
    warnings = []
    
    # Ручные значения
    if manual is None:
        manual = {}
    
    # Поиск опор
    a37 = find_a37(B, ref)
    a37s = find_a37(S, ref)
    
    # Ручные значения для a37 и a37_sample
    if "a37" in manual:
        a37 = {
            "ok": True,
            "reason": "",
            "x": float(manual["a37"]),
            "err": 0.05 * r0,
            "fwhm": ref["anchors"]["a37"]["fwhm"] * r0,
            "signif": float("nan"),
            "manual": True
        }
    
    if "a37_sample" in manual:
        a37s = {
            "ok": True,
            "reason": "",
            "x": float(manual["a37_sample"]),
            "err": 0.05 * r0,
            "fwhm": ref["anchors"]["a37"]["fwhm"] * r0,
            "signif": float("nan"),
            "manual": True
        }
    
    # Поиск суммы
    sm_ = find_sum(B, a37["x"], ref, cal_bg)
    
    if "sum" in manual:
        sm_ = {
            "ok": True,
            "reason": "",
            "x": float(manual["sum"]),
            "err": 0.05 * r0,
            "fwhm": ref["anchors"]["sum"]["fwhm"] * r0,
            "signif": float("nan"),
            "manual": True,
            "x_pred": float(manual["sum"]),
            "x_cal": None,
            "note": ""
        }
    
    # Временное отображение для поиска края
    A = ref["anchors"]
    pairs_m1 = [
        (A["a37"]["x"], a37["x"], a37["err"]),
        (A["sum"]["x"], sm_["x"], sm_["err"])
    ]
    m1, _ = build_map(pairs_m1, ref, r0)
    
    # Поиск края
    if "edge" in manual:
        edge = {
            "ok": True,
            "reason": "",
            "x": float(manual["edge"]),
            "err": 0.5 * r0,
            "sigma": ref["anchors"]["edge"]["sigma"] * r0,
            "x_pred": float(manual["edge"]),
            "manual": True
        }
    else:
        edge = find_edge(B, m1, ref)
    
    # Сбор пар для основного отображения
    pairs = []
    anchor_names = ["a37", "sum", "edge"]
    anchor_data = {"a37": a37, "sum": sm_, "edge": edge}
    
    for name in anchor_names:
        d = anchor_data[name]
        if d["ok"]:
            pairs.append((A[name]["x"], d["x"], d["err"]))
        else:
            warnings.append("опора %s не найдена уверенно: %s" % (name, d["reason"]))
    
    if len(pairs) < 2:
        # Используем все опоры, даже неуверенные
        pairs = []
        for name in anchor_names:
            d = anchor_data[name]
            pairs.append((A[name]["x"], d["x"], d["err"]))
        warnings.append("меньше двух уверенных опор -- отображение по неуверенным, проверьте вручную")
    
    if sm_.get("note"):
        warnings.append(sm_["note"])
    
    # Основное отображение
    m, info = build_map(pairs, ref, r0)
    
    # Проверка
    checks = verify(B, m, ref, anchor_data)
    for c in checks:
        if not c["ok"]:
            warnings.append(
                "проверка опоры %s: повторная подгонка %.2f, по отображению %.2f, расхождение %.2f > допуска %.2f" %
                (c["name"], c["x_meas"], c["x_map"], c["diff"], c["tol"])
            )
    
    # Проверка тождественности
    ident = bool(is_identity(m, n, ref))   # правка Claude: bool, не np.bool_
    
    if ident:
        register_kw = dict(ref["register"])
        sum_window = list(ref["sum_window"])
        fwhm_points = [tuple(p) for p in ref["fwhm_points"]]
        map_used = [0.0, 1.0, 0.0]
    else:
        register_kw = map_register_kw(m, ref, a37["x"], a37s["x"], n)
        sum_window = map_sum_window(m, ref, n)
        fwhm_points = map_fwhm_points(m, ref["fwhm_points"])
        map_used = list(m)
    
    ok = bool(all(anchor_data[name]["ok"] for name in anchor_names) and all(c["ok"] for c in checks))
    
    result = {
        "n_ch": n,
        "r": info["scale"],
        "map": [float(v) for v in m],
        "map_used": map_used,
        "kind": info["kind"],
        "resid": info["resid"],
        "identity": ident,
        "anchors": anchor_data,
        "sample_a37": a37s,
        "checks": checks,
        "warnings": warnings,
        "ok": ok,
        "register_kw": register_kw,
        "sum_window": sum_window,
        "fwhm_points": fwhm_points
    }
    
    if verbose:
        print("[AUTO] каналов %d, эталон %d; опоры фона:" % (n, ref["n_ch"]))
        for name in anchor_names:
            d = anchor_data[name]
            suffix = ""
            if d.get("manual"):
                suffix += ", вручную"
            if not d["ok"]:
                suffix += ", НЕУВЕРЕННО: %s" % d["reason"]
            print("[AUTO]   %-5s канал %.3f ± %.3f%s" % (name, d["x"], d["err"], suffix))
        
        suffix_s = ""
        if a37s.get("manual"):
            suffix_s += ", вручную"
        print("[AUTO]   пик 37 кэВ образца канал %.3f%s" % (a37s["x"], suffix_s))
        
        m0, m1_c, m2_c = m[0], m[1], m[2] if len(m) > 2 else 0.0
        print("[AUTO] отображение %s: x_прибор = %.6g + %.8g*x + %.4g*x^2 (по возрастанию степени), масштаб %.5f кан/кан эталона" %
              (info["kind"], m0, m1_c, m2_c, info["scale"]))
        
        resid_str = ", ".join("%+.3f" % r for r in info["resid"])
        print("[AUTO] невязки опор: " + resid_str)
        
        for c in checks:
            status = "OK" if c["ok"] else "НЕ ПРОЙДЕНА"
            print("[AUTO] проверка %s: подгонка %.3f, по отображению %.3f, расхождение %+.3f (допуск %.3f) %s" %
                  (c["name"], c["x_meas"], c["x_map"], c["diff"], c["tol"], status))
        
        if ident:
            print("[AUTO] отображение ≈ тождественное -- пресеты без пересчёта")
        else:
            print("[AUTO] пресеты пересчитаны: маска %s, грубое окно %s, якорь фона %.3f, образца %.3f, окно суммы %s" %
                  (register_kw["mask"], register_kw["coarse_window"], register_kw["anchor_bg"], register_kw["anchor_sample"], sum_window))
        
        for w in warnings:
            print("[AUTO] ПРЕДУПРЕЖДЕНИЕ: " + w)
    
    return result


def ecal_cfg_for(result, cfg):
    """Возвращает конфигурацию ecal, пересчитанную при необходимости."""
    if result["identity"]:
        return copy.deepcopy(cfg)
    else:
        return map_ecal_cfg(result["map"], cfg)


def selftest_cfg_for(result, cfg, ref=None):
    """Возвращает конфигурацию selftest, пересчитанную при необходимости."""
    ref = ref or load_ref()
    if result["identity"]:
        return copy.deepcopy(cfg)
    else:
        return map_selftest_cfg(
            result["map"],
            cfg,
            ref,
            result["anchors"]["a37"]["x"],
            result["sample_a37"]["x"],
            result["n_ch"]
        )


def fwhm_points_for(result, pts):
    """Возвращает точки FWHM, пересчитанные при необходимости."""
    if result["identity"]:
        return [tuple(p) for p in pts]
    else:
        return map_fwhm_points(result["map"], pts)
