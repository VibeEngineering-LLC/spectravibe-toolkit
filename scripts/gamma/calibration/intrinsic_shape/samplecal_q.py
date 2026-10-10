"""Качество калибровки по пробе: степень шкалы по LOO-СКО, сравнение с калибровкой файла и абсолютный порог, двойной фит мультиплета, ранжирование линий (близость × интенсивность × ряды, видимые в пробе; дублеты одного нуклида — эффективная энергия), проверка физичности усиления при переносе. Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-10-10."""

import math
import numpy as np
from scipy.optimize import curve_fit
from .core import poly_eval, poly_deriv, poly_inverse
from . import pairref as pr

DEG_MAX = 4                # BecqMoni допускает полином степени <= 4
CENTROID_ERR_MAX = 0.25    # доля ПШПВ: ошибка центроида больше -- репер не создаётся
SIGNIF_REFUSE = 3.0        # значимость ниже -- репер не создаётся
MULTIPLET_RATIO = 1.5      # ПШПВ больше ожидаемой во столько раз -- «вероятно мультиплет»
GAIN_RANGE = (0.95, 1.05)  # ±5 % от усиления по калибровкам файлов: наблюдаемый разброс 0.993-1.013 (±1.3 %), запас ~4x
SERIES = {"Ra-226": "U", "Pb-214": "U", "Bi-214": "U", "Pb-210": "U",
          "Ac-228": "Th", "Pb-212": "Th", "Bi-212": "Th", "Tl-208": "Th",
          "K-40": "K", "Cs-137": "Cs"}


def select_degree(cfg, ctx, verbose=True):
    from . import pipeline as pl
    n = len(cfg["points"])
    top = max(1, min(DEG_MAX, n - 2))
    all_results = []
    for deg in range(1, top + 1):
        res = pl._guarded_calibrate(dict(cfg, degree=deg), ctx, [0.0, 1.0, 0.0], False)
        loo = res.get("loo_rms", float("nan"))   # правка Claude: LOO по ВСЕМ реперам (крайние -- честная цена охвата)
        all_results.append((deg, loo, res))
    
    valid = [(d, l, r) for d, l, r in all_results if not r.get("rejected") and np.isfinite(l)]
    
    if not valid:
        fallback_deg = min(2, top)
        res = pl._guarded_calibrate(dict(cfg, degree=fallback_deg), ctx, [0.0, 1.0, 0.0], verbose)
        res["deg_table"] = []
        return res
    
    best_deg, best_loo, best_res = min(valid, key=lambda x: (round(x[1], 3), x[0]))
    
    if verbose:
        parts = []
        for d, l, r in all_results:
            suffix = "" if not r.get("rejected") else " (отвергнута)"
            parts.append("%d: %.3f кэВ%s" % (d, l, suffix))
        print("[SCAL] степень шкалы по LOO-СКО: " + ", ".join(parts) + " -> выбрана %d" % best_deg)
    
    result = dict(best_res)
    result["deg"] = best_deg
    result["deg_table"] = [[d, float(l)] for d, l, r in all_results]
    return result


def fit_doublet(counts, x_click, fwhm_ch):
    y = np.asarray(counts, float)
    hw = int(math.ceil(2.5 * fwhm_ch))
    lo = max(0, int(round(x_click)) - hw)
    hi = min(len(y) - 1, int(round(x_click)) + hw)
    
    if hi - lo < 8:
        return None
    
    x = np.arange(lo, hi + 1, dtype=float)
    yy = y[lo:hi + 1]
    sig = np.sqrt(np.maximum(yy, 1.0))
    s0 = fwhm_ch / 2.3548
    
    def f(x, a1, c1, a2, c2, s, b0, b1):
        return a1 * np.exp(-0.5 * ((x - c1) / s) ** 2) + a2 * np.exp(-0.5 * ((x - c2) / s) ** 2) + b0 + b1 * (x - x_click)
    
    b0_init = min(yy[:3].mean(), yy[-3:].mean())
    b1_init = 0.0
    c1_init = x_click - 0.35 * fwhm_ch
    c2_init = x_click + 0.35 * fwhm_ch
    a_init = max(yy.max() - b0_init, 1.0) / 2.0
    s_init = s0
    
    p0 = [a_init, c1_init, a_init, c2_init, s_init, b0_init, b1_init]
    
    bounds = ([0, lo, 0, lo, 0.5 * s0, -np.inf, -np.inf],
              [np.inf, hi, np.inf, hi, 1.5 * s0, np.inf, np.inf])
    
    try:
        popt, pcov = curve_fit(f, x, yy, p0=p0, bounds=bounds, sigma=sig, absolute_sigma=True, maxfev=20000)
    except Exception:
        return None
    
    perr = np.sqrt(np.diag(pcov))
    if not np.all(np.isfinite(perr)):
        return None
    
    a1, c1, a2, c2, s, b0, b1 = popt
    da1, dc1, da2, dc2, ds, db0, db1 = perr
    
    results = []
    for a, c, da, dc in [(a1, c1, da1, dc1), (a2, c2, da2, dc2)]:
        A = a * s * math.sqrt(2 * math.pi)
        if a > 0:
            dA = A * math.hypot(da / a, ds / s)
        else:
            dA = float('inf')
        signif = A / dA if dA > 0 else 0.0
        results.append({
            "ch": float(c),
            "ch_err": float(dc),
            "fwhm": float(2.3548 * s),
            "signif": float(signif),
            "area": float(A),
            "warn": []
        })
    
    results.sort(key=lambda r: r["ch"])
    return results


def merge_lines(lines, fwhm_keV):
    if not lines:
        return []
    
    # Группируем по нуклиду
    from collections import defaultdict
    by_nuclide = defaultdict(list)
    for l in lines:
        by_nuclide[l["nuclide"]].append(l)
    
    merged = []
    for nuclide, nuclide_lines in by_nuclide.items():
        # Сортируем по энергии
        nuclide_lines.sort(key=lambda x: x["E"])
        
        # Группируем в цепочки
        groups = []
        current_group = [nuclide_lines[0]]
        
        for i in range(1, len(nuclide_lines)):
            prev = nuclide_lines[i - 1]
            curr = nuclide_lines[i]
            if abs(curr["E"] - prev["E"]) < 0.5 * fwhm_keV:
                current_group.append(curr)
            else:
                groups.append(current_group)
                current_group = [curr]
        groups.append(current_group)
        
        # Обрабатываем каждую группу
        for group in groups:
            if len(group) == 1:
                l = group[0]
                merged.append({
                    "nuclide": l["nuclide"],
                    "E": l["E"],
                    "I": l["I"],
                    "components": []
                })
            else:
                I_total = sum(l["I"] for l in group)
                if I_total == 0:
                    E_eff = sum(l["E"] for l in group) / len(group)
                else:
                    E_eff = sum(l["E"] * l["I"] for l in group) / I_total
                
                components = [(l["E"], l["I"]) for l in group]
                dE = (max(l["E"] for l in group) - min(l["E"] for l in group)) / 2.0
                
                merged.append({
                    "nuclide": group[0]["nuclide"],
                    "E": E_eff,
                    "I": I_total,
                    "components": components,
                    "dE": dE
                })
    
    return merged


def rank_candidates(lines, E_est, fwhm_keV, visible_series=(), n_max=15):
    m = merge_lines(lines, fwhm_keV)
    candidates = []
    
    for l in m:
        d = abs(l["E"] - E_est)
        I = l["I"]
        w = math.sqrt(max(I, 0.01) / 100.0)
        
        nuclide = l["nuclide"]
        ser = SERIES.get(nuclide)
        if ser in visible_series:
            w *= 2.0
        
        score = d / w if w > 0 else float('inf')
        foreign = nuclide in pr.EXTRA_NUCLIDES
        dE_est = l["E"] - E_est
        
        note_parts = []
        if l.get("components"):
            comp_str = " + ".join("%.2f" % e for e, i in l["components"])
            note_parts.append("дублет: " + comp_str)
        if foreign:
            note_parts.append("чужой нуклид (не из рядов пробы)")
        if ser in visible_series:
            note_parts.append("ряд %s виден в пробе" % ser)
        
        note = "; ".join(note_parts) if note_parts else ""
        
        candidates.append({
            "nuclide": nuclide,
            "E": l["E"],
            "I": I,
            "dE_est": dE_est,
            "score": score,
            "foreign": foreign,
            "note": note
        })
    
    candidates.sort(key=lambda c: (c["foreign"], c["score"]))
    return candidates[:n_max]


def library_energy(nuclide, E, fwhm_keV, tol=0.05):
    lines = pr.library_lines((nuclide,), min_int=0.0, e_min=0.0, e_max=1e5)
    
    # Сначала одиночные линии
    for l in lines:
        if abs(l["E"] - E) <= tol:
            return (l["E"], 0.05)
    
    # Затем дублеты
    merged = merge_lines(lines, fwhm_keV)
    for g in merged:
        if g.get("components") and abs(g["E"] - E) <= tol:
            return (g["E"], max(g.get("dE", 0.0), 0.05))
    
    return None


def series_seen(repers):
    seen = set()
    for r in repers:
        nuclide = r.get("nuclide")
        if nuclide in SERIES:
            seen.add(SERIES[nuclide])
    return seen


def gain_check(c1, cal_src, cal_dst, x_mid):
    E = poly_eval(cal_src, x_mid)
    x2 = float(poly_inverse(list(cal_dst), E, hi=4.0 * x_mid + 1000.0))   # правка Claude: по умолчанию сетка обращения -- до 2248 кан
    g_file = float(poly_deriv(cal_src, x_mid)) / float(poly_deriv(cal_dst, x2))
    ratio = c1 / g_file
    ok = GAIN_RANGE[0] <= ratio <= GAIN_RANGE[1]
    return (ok, ratio, g_file)
