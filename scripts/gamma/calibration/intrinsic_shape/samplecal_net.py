"""Калибровка по пробе,
архитектура «разность на шкале пробы»: шкала XML по умолчанию — калибровка файла пробы; шкала по реперам рекомендуется
при охвате, LOO ≤ 3 кэВ и СКО калибровки файла на реперах больше LOO с запасом (эталон — энергии реперов), независимые линии
(по шкале реперов) — доп. проверка; записывается только решением пользователя; фон переводится на сетку пробы (по видимым
собственным структурам LaBr3 линейно либо точно через энергии) с сохранением дисперсии; XML — на сетке пробы.
Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-10-11."""

import pathlib
import numpy as np
from .core import load_spectrum, load_pair, poly_eval, poly_deriv, rebin_conserve_var
from .peakfit import find_peaks
from . import samplecal as sc
from . import pairref as pr
from . import robust as rb
from . import fwhm as F
from .netxml import write_xml
from .ops import fwhm_model

E_LO, E_HI = 100.0, 2700.0          # кэВ: область интереса (выше Tl-208 2614.5 -- только альфа-структуры LaBr3)
COVER_MIN = 0.60                    # доля [E_LO, E_HI], которую должны охватывать реперы
N_REPERS_MIN = 4
IND_SIGNIF_MIN = 5.0                # значимость независимой линии
IND_MIN = 3                         # независимых линий для проверки замены
LOO_MAX = 3.0                       # кэВ
DEG_MAX = 4
AMBIG_FRAC = 0.30                   # другая линия библиотеки с интенсивностью >= 30 % ближе 1 ПШПВ -- линия неоднозначна
INTRINSIC_E = (35.75, 788.74, 1435.8, 1471.0)   # собственные структуры LaBr3 -- не независимые линии


def load_own(bg_path, sample_path):
    """Правка Claude: фон и проба -- каждый на своём числе каналов (load_pair приводит пробу к сетке фона);
    последний канал-переполнение обнуляется."""
    bg, sm = load_spectrum(bg_path), load_spectrum(sample_path)
    bg["counts"][-1] = 0.0
    sm["counts"][-1] = 0.0
    return bg, sm


def inverse_extrap(C, E, lo, hi):
    """Обращение возрастающего полинома с линейной экстраполяцией за краями."""
    g = np.linspace(lo, hi, 40001)
    y = poly_eval(C, g)
    
    # Оставить монотонный префикс
    diffs = np.diff(y)
    # Найти первое место, где diff <= 0
    non_mono_idx = np.where(diffs <= 0)[0]
    if len(non_mono_idx) > 0:
        cut_idx = non_mono_idx[0] + 1
        g = g[:cut_idx]
        y = y[:cut_idx]
    
    if len(g) < 10:
        raise ValueError("шкала не монотонна")
    
    E = np.asarray(E, float)
    scalar_input = E.ndim == 0
    if scalar_input:
        E = E.reshape(1)
    
    x = np.interp(E, y, g)
    
    # Экстраполяция ниже y[0]
    mask_low = E < y[0]
    if np.any(mask_low):
        slope_low = (y[1] - y[0]) / (g[1] - g[0])
        x[mask_low] = g[0] + (E[mask_low] - y[0]) / slope_low
    
    # Экстраполяция выше y[-1]
    mask_high = E > y[-1]
    if np.any(mask_high):
        slope_high = (y[-1] - y[-2]) / (g[-1] - g[-2])
        x[mask_high] = g[-1] + (E[mask_high] - y[-1]) / slope_high
    
    if scalar_input:
        return float(x[0])
    return x


def select_degree_full(points_cfg, ctx, verbose=True):
    """Правка Claude: степень по LOO-СКО по ВСЕМ реперам, не выше n-2 и 4 -- общая реализация samplecal_q.select_degree."""
    if len(points_cfg["points"]) - 2 < 1:
        return {"rejected": "реперов %d -- степень не определена" % len(points_cfg["points"]), "loo_rms": np.nan, "deg": 0, "deg_table": []}
    return sc.q.select_degree(points_cfg, ctx, verbose)


def coverage_ok(Es):
    """Проверка охвата реперами области интереса."""
    Es = np.asarray(Es, float)
    if len(Es) < N_REPERS_MIN:
        return False, "реперов %d < %d" % (len(Es), N_REPERS_MIN)
    
    span = np.max(Es) - np.min(Es)
    if span < COVER_MIN * (E_HI - E_LO):
        return False, "охват %.0f кэВ < %.0f %% от [%g, %g] кэВ (%.0f кэВ)" % (
            span, COVER_MIN * 100, E_LO, E_HI, COVER_MIN * (E_HI - E_LO))
    
    third = (E_HI - E_LO) / 3
    lower_bound = E_LO + third
    upper_bound = E_HI - third
    
    if not np.any(Es <= lower_bound):
        return False, "нет репера в нижней трети [%g, %.1f] кэВ" % (E_LO, lower_bound)
    
    if not np.any(Es >= upper_bound):
        return False, "нет репера в верхней трети [%.1f, %g] кэВ" % (upper_bound, E_HI)
    
    return True, ""


# правка Claude (VERIFY7, «эталон -- энергии реперов»): гейт -- по энергиям реперов; независимые линии -- по шкале реперов;
# замена шкалы -- только явным применением. Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-10-11.
REC_RATIO = 1.5          # СКО калибровки файла на реперах должно превышать LOO шкалы по реперам не меньше чем в 1.5 раза
REC_MARGIN = 1.0         # ... и не меньше чем на 1 кэВ (шум на медальоне: LOO 1.79 кэВ, файл на реперах 1.43 кэВ)
IND_WORSE_RATIO = 1.5
IND_WORSE_MARGIN = 1.0   # кэВ: «заметно хуже» на независимых линиях -- предупреждение, не запрет

def _find_lines(S, scale, repers):
    """Пики пробы со значимостью >= 5, отождествлённые по шкале scale (шкала по реперам, не калибровка файла)."""
    n = len(S)
    fw = lambda x: sc.fwhm_guess(x, n)
    pk = find_peaks(S, fw, threshold=4.0, ch_min=max(15, int(0.01 * n)), ch_max=n - 20, min_signif=IND_SIGNIF_MIN)
    lib = pr.library_lines(pr.SAMPLE_NUCLIDES, min_int=1.0, e_min=E_LO, e_max=E_HI)
    strong = [l["E"] for l in lib if l["I"] >= 20.0 and l["E"] > 1500.0]
    excl = list(INTRINSIC_E) + [511.0] + [e - 511.0 for e in strong] + [e - 1022.0 for e in strong]
    
    lines = []
    for r in pk:
        f = sc.fit_click(S, r["centroid"], fw(r["centroid"]))
        if f is None or f["signif"] < IND_SIGNIF_MIN:
            continue
        
        E = float(poly_eval(scale, f["ch"]))
        fk = f["fwhm"] * abs(float(poly_deriv(scale, f["ch"])))
        
        if E < E_LO or E > E_HI:
            continue
        if any(abs(E - e) <= fk for e in excl):
            continue
        if any(abs(q["ch"] - f["ch"]) <= f["fwhm"] or abs(q["E"] - E) <= fk for q in repers):
            continue
        
        near = [l for l in lib if abs(l["E"] - E) <= 0.5 * fk]
        if not near:
            continue
        
        top = max(near, key=lambda l: l["I"])
        if any(l is not top and abs(l["E"] - top["E"]) <= fk and l["I"] >= AMBIG_FRAC * top["I"] for l in lib):   # правка Claude: неоднозначность -- по всей библиотеке, не только по near
            continue
        
        lines.append({
            "E": top["E"],
            "nuclide": top["nuclide"],
            "ch": f["ch"],
            "ch_err": f["ch_err"],
            "signif": f["signif"],
            "fwhm": f["fwhm"]
        })
    
    # дубли по E -- остаётся линия с наибольшей значимостью
    if lines:
        lines.sort(key=lambda x: x["E"])
        unique_lines = []
        i = 0
        while i < len(lines):
            best = lines[i]
            j = i + 1
            while j < len(lines) and lines[j]["E"] == best["E"]:
                if lines[j]["signif"] > best["signif"]:
                    best = lines[j]
                j += 1
            unique_lines.append(best)
            i = j
        lines = unique_lines
    
    return lines

def independent_lines(counts, scale, repers, verbose=True):
    """Независимые линии -- отождествление по шкале реперов, два прохода: второй -- по шкале, уточнённой реперами и линиями
    первого прохода (только для меток); калибровка файла не используется."""
    S = np.asarray(counts, float)
    lines = _find_lines(S, scale, repers)
    
    # второй проход -- по шкале, уточнённой реперами и линиями первого прохода (только метки)
    deg = len(scale) - 1
    xs = [q["ch"] for q in repers] + [l["ch"] for l in lines]
    ys = [q["E"] for q in repers] + [l["E"] for l in lines]
    
    if len(xs) >= deg + 3:
        ref = list(np.polynomial.polynomial.polyfit(xs, ys, deg))
        lines = _find_lines(S, ref, repers)
    
    if verbose:
        energies_str = ", ".join("%.1f" % l["E"] for l in lines)
        print("[SCAL] независимых линий (по шкале реперов): %d (%s)" % (len(lines), energies_str))
    
    return lines

def rms_on(C, lines):
    if not lines:
        return np.nan
    return float(np.sqrt(np.mean(np.square([poly_eval(C, l["ch"]) - l["E"] for l in lines]))))

def decide_scale(counts, cal_s, repers, verbose=True):
    """Рекомендация шкалы XML. Эталон -- энергии, присвоенные реперам: СКО калибровки файла на реперах против LOO шкалы
    по реперам. По умолчанию шкала XML -- калибровка файла пробы (replaced=False); замена -- только apply_reper_scale()."""
    S = np.asarray(counts, float)
    good = sc.valid_repers([r for r in repers if r["E"] >= E_LO])[0]
    
    out = {
        "coeffs": list(cal_s),
        "coeffs_file": list(cal_s),
        "coeffs_new": [],
        "recommend": False,
        "replaced": False,
        "why": "",
        "warnings": [],
        "deg": 0,
        "loo_rms": np.nan,
        "deg_table": [],
        "points": [],
        "independent": [],
        "rms_new": np.nan,
        "rms_file": np.nan,
        "rms_file_repers": np.nan,
        "chi2ndf": np.nan,
        "rms": np.nan,
        "rows_repers": [],
        "rows_ind": []
    }
    
    def stop(why):
        out["why"] = why
        if verbose:
            print("[SCAL] шкала XML: калибровка файла пробы -- " + why)
        return out
    
    if good:
        out["rms_file_repers"] = rms_on(cal_s, good)
    
    ok, why = coverage_ok([r["E"] for r in good])
    if not ok:
        return stop(why)
    
    points = [
        {
            "name": r.get("name", ""),
            "spectrum": "bg",
            "kind": "given",
            "ch": r["ch"],
            "ch_err": r.get("ch_err", 0.0),
            "E": r["E"],
            "dE": r.get("dE", 0.0),
            "x0": r["ch"],
            "w0": 3.0,
            "hw": 0
        } for r in good
    ]
    cfg = {"degree": 2, "points": points}
    ctx = {"B": S, "S": S, "cal": list(cal_s), "n": len(S), "auto": {"anchors": {}}}
    res = select_degree_full(cfg, ctx, verbose)
    
    out["deg"] = res.get("deg", 0)
    out["deg_table"] = res.get("deg_table", [])
    out["loo_rms"] = res.get("loo_rms", np.nan)
    out["points"] = res.get("points", [])
    out["chi2ndf"] = res.get("chi2ndf", np.nan)
    out["rms"] = res.get("rms", np.nan)
    out["coeffs_new"] = list(res.get("coeffs", []))
    
    if res.get("rejected"):
        return stop("шкала по реперам отвергнута: " + res["rejected"])
    
    loo = out["loo_rms"]
    if not np.isfinite(loo) or loo > LOO_MAX:
        return stop("LOO-СКО по всем реперам %.2f кэВ > %.1f кэВ" % (loo, LOO_MAX))
    
    loo_by = {p.get("name"): p.get("loo", np.nan) for p in out["points"]}
    out["rows_repers"] = [
        (
            r.get("name", ""),
            float(r["E"]),
            float(r["ch"]),
            float(poly_eval(cal_s, r["ch"]) - r["E"]),
            float(poly_eval(out["coeffs_new"], r["ch"]) - r["E"]),
            float(loo_by.get(r.get("name", ""), np.nan))
        ) for r in good
    ]
    
    ind = independent_lines(S, out["coeffs_new"], good, verbose)
    out["independent"] = ind
    out["rms_new"] = rms_on(out["coeffs_new"], ind)
    out["rms_file"] = rms_on(cal_s, ind)
    out["rows_ind"] = [
        (
            l["nuclide"],
            float(l["E"]),
            float(l["ch"]),
            float(poly_eval(cal_s, l["ch"]) - l["E"]),
            float(poly_eval(out["coeffs_new"], l["ch"]) - l["E"])
        ) for l in ind
    ]
    
    if len(ind) < IND_MIN:
        out["warnings"].append("доп. проверки нет: независимых линий %d < %d" % (len(ind), IND_MIN))
    else:
        if out["rms_new"] > max(IND_WORSE_RATIO * out["rms_file"], out["rms_file"] + IND_WORSE_MARGIN):
            out["warnings"].append("на %d независимых линиях шкала по реперам заметно хуже калибровки файла: %.2f против %.2f кэВ" % (len(ind), out["rms_new"], out["rms_file"]))
    
    rf = out["rms_file_repers"]
    if not (rf >= REC_RATIO * loo and rf - loo >= REC_MARGIN):
        return stop("СКО калибровки файла на реперах %.2f кэВ не превышает LOO шкалы по реперам %.2f кэВ с запасом (x%.1f и +%.1f кэВ)" % (rf, loo, REC_RATIO, REC_MARGIN))
    
    out["recommend"] = True
    if verbose:
        print("[SCAL] рекомендована шкала по реперам: степень %d, LOO %.2f кэВ, калибровка файла на реперах %.2f кэВ" % (out["deg"], loo, rf))
    
    return out

def apply_reper_scale(dec):
    """Применение шкалы по реперам (решение пользователя): только при рекомендации, иначе ValueError с причиной."""
    if not dec.get("recommend"):
        raise ValueError("шкала по реперам не рекомендована: " + (dec.get("why") or "нет шкалы"))
    return dict(dec, coeffs=list(dec["coeffs_new"]), replaced=True)

def log_decision(dec):
    """Журнал: решение и обе шкалы."""
    print("[SCAL] калибровка файла пробы: " + ", ".join("%.10g" % c for c in dec["coeffs_file"]))
    
    if dec["coeffs_new"]:
        print("[SCAL] шкала по реперам (степень %d, LOO %.2f кэВ): " % (dec["deg"], dec["loo_rms"]) + ", ".join("%.10g" % c for c in dec["coeffs_new"]))
    
    print("[SCAL] СКО калибровки файла на реперах %.2f кэВ; на независимых линиях: по реперам %.2f, файл %.2f кэВ (%d линий)" % (
        dec["rms_file_repers"], dec["rms_new"], dec["rms_file"], len(dec["independent"])
    ))
    
    if dec["replaced"]:
        decision_str = "шкала по реперам (применена)"
    elif dec["recommend"]:
        decision_str = "калибровка файла пробы; шкала по реперам рекомендована, применяет пользователь"
    else:
        decision_str = "калибровка файла пробы: " + dec["why"]
    
    print("[SCAL] решение: шкала XML -- " + decision_str)
    
    for w in dec["warnings"]:
        print("[SCAL] предупреждение: " + w)


def bg_edges_on_sample(n_b, transfer):
    """Границы каналов фона в координатах каналов пробы."""
    xb = np.arange(n_b + 1) - 0.5
    
    if transfer["mode"] == "visible":
        c0, c1 = transfer["h"]
        edges = (xb - c0) / c1
    elif transfer["mode"] == "energy":
        E = poly_eval(transfer["cal_b"], xb)
        edges = inverse_extrap(transfer["scale_s"], E, -0.5, transfer["n_s"] - 0.5 + 200)
    else:
        raise ValueError("неизвестный режим перевода фона: %s" % transfer["mode"])
    
    # Проверка строгого возрастания
    if np.any(np.diff(edges) <= 0):
        raise ValueError("перевод фона на пробу: границы каналов не возрастают")
    
    return edges


def bg_on_sample(B, n_s, transfer):
    """Перевод фона на сетку пробы с сохранением дисперсии."""
    edges = bg_edges_on_sample(len(B), transfer)
    vals, var = rebin_conserve_var(np.asarray(B, float), edges, n_s)
    return vals, var


def run_sample_cal(bg_path, sample_path, repers, use_intrinsic=True, mode=None, verbose=True, use_reper_scale=False):
    """Основная функция калибровки по пробе."""
    rb.check_cancel("калибровка по пробе")
    
    bg, sm = load_own(bg_path, sample_path)
    S = sm["counts"].astype(float)
    B = bg["counts"].astype(float)
    cal_s = list(sm["cal"])
    cal_b = list(bg["cal"])
    n_s = len(S)
    warns = []
    
    # Опоры пробы
    intr, S_anc = sc.intrinsic_repers(S, cal_s)
    intr = sc.check_intrinsic(S_anc, intr, repers, S, cal_s, verbose)
    warns += ["опора %s пробы отброшена: %s" % (k, v["why"]) for k, v in S_anc.items() if v.get("why")]
    
    # Все реперы
    allr = list(repers)
    if use_intrinsic:
        existing_names = {r.get("name") for r in allr}
        for r in intr:   # правка Claude: intr -- список реперов
            if r["name"] not in existing_names:
                allr.append(r)
    
    # Решение о шкале
    dec = decide_scale(S, cal_s, allr, verbose)
    if use_reper_scale:   # правка Claude: замену применяет пользователь (кнопка GUI / --use-reper-scale)
        dec = apply_reper_scale(dec)
    warns += dec["warnings"]
    if dec["recommend"] and not dec["replaced"]:
        warns.append("шкала по реперам рекомендована -- применить: «Записать шкалу по реперам» (по умолчанию -- калибровка файла пробы)")
    elif not dec["replaced"]:
        warns.append("шкала XML -- калибровка файла пробы: " + dec["why"])
    if verbose:
        log_decision(dec)
    
    # Перевод фона на пробу
    rb.check_cancel("перевод фона на пробу")
    B_anc = sc.anchors_of(B, cal_b)
    reg = sc.register_by_visible(S_anc, B_anc)
    
    if mode is None:
        mode = "visible" if reg is not None else "energy"
    
    if mode == "visible":
        if reg is None:
            raise ValueError("видимых общих структур меньше двух")
        transfer = {
            "mode": "visible",
            "h": [reg["h"][0], reg["h"][1]],
            "resid": reg["resid"],
            "keys": reg["keys"]
        }
    else:  # mode == "energy"
        from . import pipeline as pl
        own_b = pl.run_ecal(bg_path, bg_path, [0.0, 1.0, 0.0], low=False, verbose=False)
        cal_b_own, how = sc.bg_own_scale(own_b, cal_b, B)
        transfer = {
            "mode": "energy",
            "cal_b": cal_b_own,
            "scale_s": dec["coeffs"],
            "n_s": n_s,
            "bg_scale_how": how
        }
        warns.append("видимых собственных структур в пробе меньше двух -- фон переведён на шкалу пробы через энергии (фон: %s)" % how)
    
    if verbose:
        if mode == "visible":
            print("[SCAL] фон переведён на пробу по видимым структурам")
        else:
            print("[SCAL] фон переведён на пробу через энергии")
    
    # ПШПВ на сетке пробы
    pts = []
    for r in allr:
        fwhm_val = r.get("fwhm", np.nan)
        if not np.isfinite(fwhm_val):
            continue
        warn_list = r.get("warn", [])
        if any("ПШПВ" in w for w in warn_list):
            continue
        source = r.get("source")
        signif = r.get("signif", np.nan)
        if source == sc.SRC_LABR or (np.isfinite(signif) and signif >= sc.SIGNIF_MIN):
            pts.append((r["ch"], fwhm_val))
    
    pts.sort(key=lambda x: x[0])
    
    # правка Claude: SLSQP в fit_fwhm сходится на 8192 кан только с точкой у начала шкалы -- сначала точки пробы,
    # затем они же с опорной точкой прибора на 1.5 % шкалы, затем опорные точки прибора
    low = (0.015 * n_s, sc.fwhm_guess(0.015 * n_s, n_s))
    ref = [(x * n_s / 2048.0, sc.fwhm_guess(x * n_s / 2048.0, n_s)) for x in (29.2, 142.87, 211.97, 338.14, 847.1)]   # каналы точек ПШПВ пресета
    for cand, note in ((pts, ""), (sorted(pts + [low]), "к точкам пробы добавлена опорная точка ПШПВ прибора у начала шкалы"),
                       (ref, "точки ПШПВ пробы не дают годной модели -- взяты опорные точки ПШПВ прибора")):
        if len(cand) < 3:
            continue
        try:
            F.fit_fwhm(n_s, [p[0] for p in cand], [p[1] for p in cand])
        except ValueError:
            continue
        pts = cand
        if note:
            warns.append(note)
        break
    else:
        pts = ref
        warns.append("модель ПШПВ не подобрана ни по одному набору точек -- взяты опорные точки ПШПВ прибора без проверки")

    return {
        "mode": mode,
        "transfer": transfer,
        "decision": dec,
        "coeffs": dec["coeffs"],
        "coeffs_xml": dec["coeffs"],
        "rejected": "" if dec["replaced"] else dec["why"],
        "points": dec["points"],
        "chi2ndf": dec["chi2ndf"],
        "rms": dec["rms"],
        "loo_rms": dec["loo_rms"],
        "deg": dec["deg"],
        "fwhm_points": [[float(a), float(b)] for a, b in pts],
        "warnings": warns,
        "repers": allr,
        "anchors_sample": S_anc,
        "anchors_bg": B_anc,
        "visible_sample": [k for k in S_anc if S_anc[k]["visible"]],
        "k": float(sm["live"] / bg["live"]),
        "n_ch": n_s,
        "lu": False,
        "lines": []
    }


def run_sample_subtract(bg_path, sample_path, cal_res, fw_points, out_zeroed, out_signed, name=None, date="", monotone_to=1.5, lt_out="sample", verbose=True):
    """Вычитание фона на сетке пробы."""
    rb.check_cancel("перед вычитанием")
    
    bg, sm = load_own(bg_path, sample_path)
    S = sm["counts"].astype(float)
    B = bg["counts"].astype(float)
    k = sm["live"] / bg["live"]
    
    Bs, varB = bg_on_sample(B, len(S), cal_res["transfer"])
    resid = S - k * Bs
    var = np.maximum(S + k * k * varB, 1.0)
    
    if lt_out not in ("sample", "bg"):
        raise ValueError("lt_out должен быть 'sample' или 'bg'")
    
    ref = sm if lt_out == "sample" else bg
    f_out = ref["live"] / sm["live"]
    times = None if lt_out == "sample" else (ref["live"], ref["real"])
    
    zer = np.maximum(0.0, np.round(resid * f_out))
    sig = np.round(resid * f_out)
    
    fwc = fwhm_model(fw_points, len(S), monotone_to)
    C = list(cal_res["coeffs_xml"])
    
    if cal_res["transfer"]["mode"] == "visible":
        H = [cal_res["transfer"]["h"][0], cal_res["transfer"]["h"][1], 0.0]
    else:
        H = [0.0, 1.0, 0.0]
    
    if name is None:
        name = pathlib.Path(sm["name"]).stem + " - фон (на шкале пробы)"
    
    write_xml(sample_path, out_zeroed, zer, False, C, H, k, fwc, fw_points, name, date, times=times)
    write_xml(sample_path, out_signed, sig, True, C, H, k, fwc, fw_points, name, date, times=times)
    
    if verbose:
        print("[NET] разность на сетке пробы: net_cps = S/LT_s - B'/LT_b (B' -- фон на сетке пробы, %s), каналов %d" % (
            cal_res["transfer"]["mode"], len(S)))
    
    return {
        "k": k,
        "resid": resid,
        "var": var,
        "fwhm": fwc,
        "zeroed": zer,
        "signed": sig,
        "lt_out": lt_out,
        "live_out": float(ref["live"]),
        "live_s": float(sm["live"]),
        "live_b": float(bg["live"]),
        "coeffs_xml": C,
        "bg_on_sample": Bs,
        "n_ch": len(S)
    }
