"""Вердикт по окнам собственных структур на остатке
без модели, подогнанной на тех же каналах: базовая линия окна -- по боковым полосам вне окна.

Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-10-10"""

import json
import math
import pathlib

import numpy as np

from .core import poly_eval, poly_inverse
from .robust import step_col, gauss_col

HERE = pathlib.Path(__file__).resolve().parent
HIDE_MAX = 1.0 / 3.0      # доля скрытых каналов в любой половине окна, выше -- «не оценено»
SIDE_FWHM = 2.0           # ширина боковой полосы, ПШПВ
SIDE_MIN = 3              # минимум каналов в каждой боковой полосе
NFREE_MIN = 5             # минимум свободных каналов окна (после подгонки K-40)
Z_THR = 3.0
INFL_MAX = 2.0            # рост шума суммы/диполя от неопределённости базовой линии, выше -- «не оценено»


def load_windows_kev(path=None):
    """Читает JSON с окнами и возвращает список словарей."""
    if path is None:
        path = HERE / "presets" / "labr3_struct_windows.json"
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    windows = data.get("windows", [])
    result = []
    for w in windows:
        result.append({
            "name": w.get("name", ""),
            "e_lo": float(w.get("e_lo", 0.0)),
            "e_hi": float(w.get("e_hi", 0.0)),
            "source": w.get("source", "")
        })
    return result


def windows_to_channels(wins_kev, scale, n):
    """Преобразует окна в кэВ в каналы."""
    result = {}
    for w in wins_kev:
        name = w["name"]
        e_lo = w["e_lo"]
        e_hi = w["e_hi"]
        source = w.get("source", "")
        try:
            lo = int(math.ceil(poly_inverse(scale, e_lo, hi=n + 200)))
            hi = int(math.floor(poly_inverse(scale, e_hi, hi=n + 200)))
        except Exception:
            continue
        lo = max(0, lo)
        hi = min(n - 2, hi)
        if hi - lo < 4:
            continue
        result[name] = {
            "lo": lo,
            "hi": hi,
            "e_lo": e_lo,
            "e_hi": e_hi,
            "source": source
        }
    return result


def _hidden(lo, hi, excl):
    """Возвращает булев массив: True, если канал входит в исключённые пары."""
    length = hi - lo + 1
    hid = np.zeros(length, dtype=bool)
    for pair in excl:
        a, b = pair
        # Каналы в диапазоне [lo, hi], которые попадают в [a, b]
        start = max(lo, a)
        end = min(hi, b)
        if start <= end:
            idx_start = start - lo
            idx_end = end - lo
            hid[idx_start:idx_end + 1] = True
    return hid


def window_verdict(net, var, name, w, excl, fwhm_fn, steps_ch=(), k40=None, z_thr=Z_THR, hide_max=HIDE_MAX, side_fwhm=SIDE_FWHM, side_min=SIDE_MIN):
    """Вычисляет вердикт для одного окна."""
    lo, hi = w["lo"], w["hi"]
    n = len(net)
    centre = 0.5 * (lo + hi)
    win = np.arange(lo, hi + 1)
    
    hid = _hidden(lo, hi, excl)
    n_win = len(win)
    hidden_frac = float(hid.mean())
    
    # Половины
    left_mask = win < centre
    right_mask = win > centre
    
    frac_l = 0.0
    frac_r = 0.0
    if np.any(left_mask):
        frac_l = float(hid[left_mask].mean())
    if np.any(right_mask):
        frac_r = float(hid[right_mask].mean())
    
    # Боковые полосы
    wb = max(side_min, int(math.ceil(side_fwhm * fwhm_fn(centre))))
    
    # Левая полоса: max(0, lo - wb) .. lo - 1
    side_l_lo = max(0, lo - wb)
    side_l_hi = lo - 1
    # Правая полоса: hi + 1 .. min(n - 2, hi + wb)
    side_r_lo = hi + 1
    side_r_hi = min(n - 2, hi + wb)
    
    # Формируем списки каналов боковых полос, исключая excl и var <= 0
    def filter_side(lo_s, hi_s):
        if lo_s > hi_s:
            return []
        channels = []
        for ch in range(lo_s, hi_s + 1):
            # Проверка на исключение
            is_excl = False
            for pair in excl:
                a, b = pair
                if a <= ch <= b:
                    is_excl = True
                    break
            if is_excl:
                continue
            if var[ch] <= 0:
                continue
            channels.append(ch)
        return channels
    
    side_l_ch = filter_side(side_l_lo, side_l_hi)
    side_r_ch = filter_side(side_r_lo, side_r_hi)
    
    n_side_l = len(side_l_ch)
    n_side_r = len(side_r_ch)
    
    # Начальная запись результата
    res = {
        "name": name,
        "lo": lo,
        "hi": hi,
        "e_lo": w["e_lo"],
        "e_hi": w["e_hi"],
        "n_win": n_win,
        "hidden_frac": hidden_frac,
        "hidden_frac_l": frac_l,
        "hidden_frac_r": frac_r,
        "n_side_l": n_side_l,
        "n_side_r": n_side_r,
        "n_used": int(np.sum(~hid)),
        "z_sum": 0.0,
        "z_dip": 0.0,
        "z_chi2": 0.0,
        "chi2": 0.0,
        "ndf": 0,
        "evaluated": False,
        "ok": False,
        "status": "не оценено",
        "reason": "",
        "dependent": False,
        "k40": None,
        "n_steps": 0,
        "noise_infl": 0.0
    }
    
    # Проверка скрытых долей
    if max(frac_l, frac_r) > hide_max:
        res["reason"] = "скрыто линиями образца %.0f %% окна (левая половина %.0f %%, правая %.0f %%), порог %.0f %% на половину" % (
            hidden_frac * 100, frac_l * 100, frac_r * 100, hide_max * 100
        )
        return res
    
    # Проверка боковых полос
    if n_side_l < side_min or n_side_r < side_min:
        res["reason"] = "нет боковых полос вне окна (каналов слева %d, справа %d)" % (n_side_l, n_side_r)
        return res
    
    # Каналы подгонки базовой линии
    xs = np.array(side_l_ch + side_r_ch, dtype=int)
    
    if len(xs) == 0:
        res["reason"] = "нет каналов для подгонки базовой линии"
        return res
    
    # Матрица базовой линии
    # Столбцы: 1, (x - centre) / max(1.0, hi - lo)
    denom = max(1.0, float(hi - lo))
    
    # Определяем, какие ступени попадают в диапазон [lo - wb, hi + wb]
    steps_in_range = []
    for s in steps_ch:
        if lo - wb <= s <= hi + wb:
            steps_in_range.append(s)
    
    def build_matrix(x_arr):
        """Строит матрицу M для каналов x_arr."""
        n_ch = len(x_arr)
        # Два базовых столбца
        M = np.zeros((n_ch, 2 + len(steps_in_range)))
        M[:, 0] = 1.0
        M[:, 1] = (x_arr - centre) / denom
        
        # Столбцы ступеней
        for i, s in enumerate(steps_in_range):
            # step_col(n, s, fwhm_fn(s) / 2.3548)
            # fwhm_fn(s) возвращает ПШПВ в каналах
            fwhm_s = fwhm_fn(s)
            sigma_s = fwhm_s / 2.3548
            col = step_col(n, s, sigma_s)
            M[:, 2 + i] = col[x_arr]
        
        return M
    
    M_full = build_matrix(xs)
    
    # Взвешенный МНК
    Wt = 1.0 / np.sqrt(var[xs])
    A = M_full * Wt[:, None]
    y = net[xs] * Wt
    
    # Проверка ранга
    rank = np.linalg.matrix_rank(A)
    n_cols = A.shape[1]
    
    if rank < n_cols:
        # Отбрасываем столбцы ступеней (оставляем два первых)
        A = A[:, :2]
        M_full = M_full[:, :2]
        steps_in_range = []
    
    beta, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    cov_b = np.linalg.pinv(A.T @ A)
    
    res["n_steps"] = len(steps_in_range)
    
    # Каналы окна (нескрытые)
    u_idx = win[~hid]
    
    if len(u_idx) == 0:
        res["reason"] = "нет свободных каналов в окне"
        return res
    
    Mw = build_matrix(u_idx)
    u = net[u_idx] - Mw @ beta
    C = np.diag(var[u_idx]) + Mw @ cov_b @ Mw.T
    
    # K-40
    nfit = 0
    if k40 is not None and lo <= k40["ch"] <= hi:
        g = gauss_col(n, k40["ch"], k40["sigma"])[u_idx]
        Ci = np.linalg.inv(C)
        a = (g @ Ci @ u) / (g @ Ci @ g)
        sa = 1.0 / np.sqrt(g @ Ci @ g)
        # площадь K-40 образца не бывает отрицательной; отрицательный «K-40» -- это пересчёт фона
        # (k или сдвиг у сумм-пика), его нельзя поглощать гауссианой. На границе a = 0 параметр не расходуется.
        nfit = 1
        if a < 0.0:
            a, nfit = 0.0, 0
        u = u - a * g
        res["k40"] = {
            "ch": float(k40["ch"]),
            "sigma": float(k40["sigma"]),
            "amp": float(a),
            "amp_err": float(sa),
            "area": float(a * k40["sigma"] * math.sqrt(2 * math.pi)),
            "area_err": float(sa * k40["sigma"] * math.sqrt(2 * math.pi))
        }
    
    ndf = len(u_idx) - nfit
    
    if ndf < NFREE_MIN:
        res["reason"] = "мало свободных каналов окна (%d)" % ndf
        return res
    
    one = np.ones(len(u_idx))
    sgn = np.sign(u_idx - centre)
    
    denom_sum = np.sqrt(one @ C @ one)
    if denom_sum == 0:
        z_sum = 0.0
    else:
        z_sum = (one @ u) / denom_sum
    
    denom_dip = np.sqrt(sgn @ C @ sgn)
    if denom_dip == 0:
        z_dip = 0.0
    else:
        z_dip = (sgn @ u) / denom_dip
    
    chi2 = u @ np.linalg.solve(C, u)
    z_chi2 = (chi2 - ndf) / np.sqrt(2 * ndf)

    # базовая линия по полосам плохо определена (ступень внутри окна против наклона) -- шум суммы и
    # диполя растёт, проверка теряет чувствительность. Рост > INFL_MAX раз -- «не оценено», а не зелёное.
    v0 = float(np.sum(var[u_idx]))
    infl = max(math.sqrt((one @ C @ one) / v0), math.sqrt((sgn @ C @ sgn) / max(float(sgn @ (var[u_idx] * sgn)), 1e-12)))
    res["noise_infl"] = float(infl)
    if infl > INFL_MAX:
        res["reason"] = "базовая линия по боковым полосам плохо определена: шум проверки вырос в %.1f раза (порог %.1f)" % (
            infl, INFL_MAX)
        return res

    res["evaluated"] = True
    res["z_sum"] = float(z_sum)
    res["z_dip"] = float(z_dip)
    res["z_chi2"] = float(z_chi2)
    res["chi2"] = float(chi2)
    res["ndf"] = int(ndf)
    res["n_used"] = int(len(u_idx))
    
    # Причины
    reasons = []
    if abs(z_sum) > z_thr:
        reasons.append("сумма")
    if abs(z_dip) > z_thr:
        reasons.append("диполь")
    if z_chi2 > z_thr:
        reasons.append("избыток χ²")
    if z_chi2 < -z_thr:
        reasons.append("χ² ниже статистики: подозрение на вложенную или зависимую пару")
        res["dependent"] = True
    
    res["ok"] = len(reasons) == 0
    res["status"] = "OK" if res["ok"] else "ВНЕ НОРМЫ"
    res["reason"] = "; ".join(reasons)
    
    return res


def structure_verdict(net, var, windows, excl, fwhm_fn, steps_ch=(), k40=None, **kw):
    """Вычисляет вердикты для всех окон."""
    net = np.asarray(net, dtype=float)
    var = np.asarray(var, dtype=float)
    
    results = []
    for name, w in windows.items():
        res = window_verdict(net, var, name, w, excl, fwhm_fn, steps_ch=steps_ch, k40=k40, **kw)
        results.append(res)
    
    return results


def summarize(checks, min_eval_frac=0.5):
    """Сводит результаты в итоговую статистику."""
    n = len(checks)
    ev = [c for c in checks if c["evaluated"]]
    n_ok = sum(1 for c in ev if c["ok"])
    
    struct_ok = (len(ev) >= math.ceil(min_eval_frac * n)) and len(ev) > 0 and all(c["ok"] for c in ev)
    dependent = any(c["dependent"] for c in checks)
    
    return {
        "struct_ok": bool(struct_ok),
        "n": n,
        "n_eval": len(ev),
        "n_ok": n_ok,
        "n_not_eval": n - len(ev),
        "dependent": bool(dependent),
        "not_eval": [c["name"] for c in checks if not c["evaluated"]]
    }


def format_table(checks, title):
    """Форматирует результаты в таблицу."""
    lines = [title + " СПРАВОЧНО: самопроверка неполная, см. ограничения"]
    
    for c in checks:
        name = c["name"]
        e_lo = c["e_lo"]
        e_hi = c["e_hi"]
        lo = c["lo"]
        hi = c["hi"]
        hidden_frac = c["hidden_frac"]
        n_used = c["n_used"]
        z_sum = c["z_sum"]
        z_dip = c["z_dip"]
        z_chi2 = c["z_chi2"]
        status = c["status"]
        reason = c["reason"]
        
        line = "[STRUCT] %-13s %6.0f-%-6.0f кэВ кан %d-%d скрыто %3.0f %% каналов %d: z_sum %+6.2f z_dip %+6.2f z_chi2 %+6.2f %s" % (
            name, e_lo, e_hi, lo, hi, hidden_frac * 100, n_used, z_sum, z_dip, z_chi2, status
        )
        
        if reason:
            line += " (" + reason + ")"
        
        lines.append(line)
        
        # Дополнительная строка для K-40
        if c["k40"] is not None:
            k40 = c["k40"]
            line_k40 = "[STRUCT]   K-40 в окне: кан %.1f, σ %.2f кан, площадь %.0f ± %.0f отсчётов" % (
                k40["ch"], k40["sigma"], k40["area"], k40["area_err"]
            )
            lines.append(line_k40)
    
    # Итоговая строка
    summ = summarize(checks)
    not_eval_str = ", ".join(summ["not_eval"]) if summ["not_eval"] else "нет"
    
    line_sum = "[STRUCT] итог (справочно, в all_ok не входит): оценено %d из %d, в норме %d; не оценено: %s" % (
        summ["n_eval"], summ["n"], summ["n_ok"], not_eval_str
    )
    
    if summ["dependent"]:
        line_sum += "; ПОДОЗРЕНИЕ НА ВЛОЖЕННУЮ/ЗАВИСИМУЮ ПАРУ"
    
    lines.append(line_sum)
    
    return "\n".join(lines)
