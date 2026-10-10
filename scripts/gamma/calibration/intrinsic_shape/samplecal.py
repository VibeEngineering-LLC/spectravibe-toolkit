"""Калибровка по пробе: реперы, выбранные на спектре пробы (пики образца с энергиями из ядерных данных проекта или
ручными), плюс видимые в пробе собственные структуры LaBr3; шкала пробы, совмещение фона по видимым структурам,
набор реперов в JSON.

Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-10-10"""

import json, math, pathlib
import numpy as np
from .core import load_spectrum, poly_eval, poly_deriv, poly_inverse
from .peakfit import fit_peak
from . import autopreset as ap
from . import pairref as pr
from . import register as reg
from . import robust as rb

HERE = pathlib.Path(__file__).resolve().parent
FORMAT = "labr3_samplecal/1"
SIGNIF_MIN = 5.0          # ниже -- предупреждение «незначимый пик»
FWHM_RATIO_MAX = 1.5      # ПШПВ пика отличается от ожидаемой больше чем во столько раз -- предупреждение
VISIBLE_SIGNIF = 10.0     # собственная структура «видна» в спектре при значимости опоры >= 10 (край -- если ok)
SRC_LIB, SRC_MANUAL, SRC_LABR = "библиотека", "ручная", "LaBr3"


def intrinsic_energies():
    """Читает пресеты LaBr3 и возвращает словарь энергий собственных структур."""
    with open(HERE / "presets" / "labr3_ecal.json", "r", encoding="utf-8") as f:
        ecal = json.load(f)
    with open(HERE / "presets" / "labr3_lowE.json", "r", encoding="utf-8") as f:
        lowE = json.load(f)

    points = ecal["points"]
    edge_pt = next(p for p in points if p["name"] == "La138_789_edge")
    sum_pt = next(p for p in points if p["name"] == "La138_sum1473")

    return {
        "a37": (lowE["e37_gamma"], lowE["de37"], "La138_37_intr"),
        "edge": (edge_pt["E"], edge_pt["dE"], "La138_789_edge"),
        "sum": (sum_pt["E"], sum_pt["dE"], "La138_sum1473"),
    }


def fit_click(counts, x_click, fwhm_ch, fwhm_expected=None):
    """Подгонка пика в окрестности клика. Возвращает словарь с параметрами или None."""
    hw = max(3, int(math.ceil(2.0 * fwhm_ch)))
    r = fit_peak(np.asarray(counts, float), float(x_click), float(fwhm_ch), half_window=hw)

    if r is None:
        return None
    if not np.isfinite(r["centroid"]) or not np.isfinite(r["centroid_err"]) or r["centroid_err"] <= 0:
        return None
    if not (r.get("area", 1.0) > 0 and np.isfinite(r["fwhm"]) and r["centroid_err"] < r["fwhm"]
            and abs(r["centroid"] - x_click) <= hw and 0.0 <= r["centroid"] <= len(counts) - 1 and r["signif"] > 0.0):
        return None   # правка Claude: подгонка выродилась -- площадь <= 0, центроид не определён (ошибка >= ПШПВ) или вне окна

    warn = []
    if r["signif"] < SIGNIF_MIN:
        warn.append("незначимый пик (значимость %.1f < %.0f)" % (r["signif"], SIGNIF_MIN))
    if fwhm_expected is not None:
        ratio = r["fwhm"] / fwhm_expected
        if ratio < 1.0 / FWHM_RATIO_MAX or ratio > FWHM_RATIO_MAX:
            warn.append("ПШПВ %.2f кан отличается от ожидаемой %.2f больше чем в %.1f раза" % (r["fwhm"], fwhm_expected, FWHM_RATIO_MAX))

    return {
        "ch": float(r["centroid"]),
        "ch_err": float(r["centroid_err"]),
        "fwhm": float(r["fwhm"]),
        "fwhm_err": float(r.get("fwhm_err", float("nan"))),
        "signif": float(r["signif"]),
        "warn": warn,
    }


def candidate_lines(E_est, nuclides=None, min_int=0.5, n_max=15, e_min=10.0, e_max=4000.0):
    """Возвращает список кандидатов на линии из библиотеки, отсортированный по близости к E_est."""
    if nuclides is None:
        nuclides = tuple(pr.SAMPLE_NUCLIDES) + tuple(pr.EXTRA_NUCLIDES)
    lines = pr.library_lines(nuclides, min_int=min_int, e_min=e_min, e_max=e_max)
    lines = sorted(lines, key=lambda l: abs(l["E"] - E_est))
    result = []
    for l in lines[:n_max]:
        d = dict(l)
        d["dE_est"] = l["E"] - E_est
        result.append(d)
    return result


def make_reper(fit, E, source, nuclide=None, dE=None):
    """Создаёт словарь репера из результатов подгонки."""
    if source not in (SRC_LIB, SRC_MANUAL, SRC_LABR):
        raise ValueError("неизвестный источник репера: %s" % source)
    if source == SRC_LIB and not nuclide:
        raise ValueError("для линии из библиотеки нужен нуклид")

    if dE is None:
        if source == SRC_LIB:
            dE = 0.05
        else:
            dE = 0.5

    if nuclide:
        name = nuclide.replace("-", "") + "_" + str(int(round(E)))
    else:
        name = "manual_%d" % round(E)

    return {
        "name": name,
        "ch": fit["ch"],
        "ch_err": fit["ch_err"],
        "fwhm": fit.get("fwhm", float("nan")),
        "signif": fit.get("signif", float("nan")),
        "E": float(E),
        "dE": float(dE),
        "source": source,
        "nuclide": nuclide or "",
        "warn": list(fit.get("warn", [])),
    }


def anchors_of(counts, cal):
    """Возвращает словарь опор LaBr3 с их видимостью и параметрами."""
    a = _autopreset_anchors(counts, cal)
    result = {}
    for key in ("a37", "edge", "sum"):
        A = a.get(key, {})
        x = A.get("x", float("nan"))
        err = A.get("err", float("nan"))
        fwhm = A.get("fwhm", float("nan"))
        signif = A.get("signif", float("nan"))
        
        visible = (
            A.get("ok", False) and
            np.isfinite(x) and
            np.isfinite(err) and
            err > 0 and
            (key == "edge" or (np.isfinite(signif) and signif >= VISIBLE_SIGNIF))
        )
        
        result[key] = {
            "x": float(x),
            "err": float(err),
            "fwhm": float(fwhm),
            "signif": float(signif),
            "visible": bool(visible),
        }
    return result


def intrinsic_repers(counts, cal):
    """Возвращает список реперов для видимых собственных структур LaBr3 и словарь опор."""
    anc = anchors_of(counts, cal)
    en = intrinsic_energies()
    repers = []
    for key in ("a37", "edge", "sum"):
        if anc[key]["visible"]:
            E, dE, name = en[key]
            fit_dict = {
                "ch": anc[key]["x"],
                "ch_err": anc[key]["err"],
                "fwhm": anc[key]["fwhm"],
                "signif": anc[key]["signif"],
            }
            r = make_reper(fit_dict, E, SRC_LABR, nuclide="La-138", dE=dE)
            r["name"] = name
            repers.append(r)
    return repers, anc


def valid_repers(repers):
    """Разделяет реперы на хорошие и отброшенные."""
    good = []
    dropped = []
    for r in repers:
        if (np.isfinite(r["ch"]) and np.isfinite(r["ch_err"]) and np.isfinite(r["E"]) and r["ch_err"] > 0):
            good.append(r)
        else:
            dropped.append(r["name"])
    return good, dropped


def scale_from_repers(repers, counts, cal, verbose=True):
    """Строит шкалу по реперам."""
    from . import pipeline as pl

    good, dropped = valid_repers(repers)
    if verbose:
        for name in dropped:
            print("[SCAL] репер %s отброшен: положение не определено" % name)

    if not good:
        raise ValueError("Шкала по пробе не построена: реперов нет")

    points = []
    for r in good:
        points.append({
            "name": r["name"],
            "spectrum": "bg",
            "kind": "given",
            "ch": r["ch"],
            "ch_err": r["ch_err"],
            "E": r["E"],
            "dE": r["dE"],
            "x0": r["ch"],
            "w0": 3.0,
            "hw": 0,
        })

    cfg = {"degree": 2, "points": points}
    ctx = {
        "B": np.asarray(counts, float),
        "S": np.asarray(counts, float),
        "cal": list(cal),
        "n": len(counts),
        "auto": {"anchors": {}},
    }

    res = pl._guarded_calibrate(cfg, ctx, [0.0, 1.0, 0.0], verbose)
    return dict(res, dropped=dropped, n_repers=len(good))


def save_json(path, repers, device, n_ch, anchors):
    """Сохраняет набор реперов в JSON."""
    data = {
        "format": FORMAT,
        "device": str(device),
        "n_ch": int(n_ch),
        "anchors": {
            k: {"x": float(anchors[k]["x"]), "visible": bool(anchors[k]["visible"])}
            for k in anchors
        },
        "repers": [
            {
                "name": r["name"],
                "ch": float(r["ch"]),
                "ch_err": float(r["ch_err"]),
                "fwhm": float(r["fwhm"]),
                "E": float(r["E"]),
                "dE": float(r["dE"]),
                "source": r["source"],
                "nuclide": r["nuclide"],
            }
            for r in repers
        ],
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=1)


def load_json(path):
    """Загружает набор реперов из JSON."""
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if data.get("format") != FORMAT:
        raise ValueError("не набор реперов калибровки по пробе: %s" % path)
    return data


def apply_to(data, counts, cal, refit=True, verbose=True):
    """Переносит набор реперов на другой спектр того же прибора."""
    dst = anchors_of(counts, cal)
    src_anchors = data["anchors"]
    
    common_keys = [k for k in ("a37", "edge", "sum") if src_anchors[k]["visible"] and dst[k]["visible"]]
    if len(common_keys) < 2:
        raise ValueError("перенос набора реперов: общих видимых опор меньше двух")

    src_x = np.array([src_anchors[k]["x"] for k in common_keys])
    dst_x = np.array([dst[k]["x"] for k in common_keys])
    m = np.polyfit(src_x, dst_x, 1)
    c0, c1 = float(m[1]), float(m[0])

    if verbose:
        print("[SCAL] перенос набора: отображение опор x' = %.4f + %.6f*x по %d опорам" % (c0, c1, len(common_keys)))

    # Получаем новые собственные структуры LaBr3
    intr_repers, _ = intrinsic_repers(counts, cal)
    intr_names = set(r["name"] for r in intr_repers)

    new_repers = []
    for r in data["repers"]:
        if r["source"] == SRC_LABR:
            # Берем заново из intrinsic_repers
            match = next((ir for ir in intr_repers if ir["name"] == r["name"]), None)
            if match:
                new_repers.append(match)
            else:
                # Если не нашли, пропускаем или оставляем как есть? Спецификация говорит "берутся заново", значит только видимые.
                # Если не видим, то не включаем.
                pass
        else:
            x = c0 + c1 * r["ch"]
            fwhm_orig = r.get("fwhm", float("nan"))
            fwhm_new = fwhm_orig * c1 if np.isfinite(fwhm_orig) else 3.0 * c1
            
            if refit:
                f = fit_click(counts, x, fwhm_new)
                if f is not None and abs(f["ch"] - x) <= 2.0 * fwhm_new:
                    new_r = dict(r)
                    new_r["ch"] = f["ch"]
                    new_r["ch_err"] = f["ch_err"]
                    new_r["fwhm"] = f["fwhm"]
                    new_r["signif"] = f["signif"]
                    new_r["warn"] = f["warn"]
                    new_repers.append(new_r)
                else:
                    new_r = dict(r)
                    new_r["ch"] = x
                    new_r["ch_err"] = math.hypot(r["ch_err"] * c1, 0.5)
                    new_r["warn"] = list(r.get("warn", [])) + ["пик на новом спектре не подогнан, перенос по опорам"]
                    new_repers.append(new_r)
            else:
                new_r = dict(r)
                new_r["ch"] = x
                new_r["ch_err"] = r["ch_err"] * abs(c1)
                if np.isfinite(fwhm_orig):
                    new_r["fwhm"] = fwhm_orig * abs(c1)
                new_repers.append(new_r)

    return new_repers, [c0, c1]


def register_by_visible(S_anc, B_anc):
    """Строит линейное отображение между видимыми опорами пробы и фона."""
    keys = [k for k in ("a37", "edge", "sum") if S_anc[k]["visible"] and B_anc[k]["visible"]]
    if len(keys) < 2:
        return None
    
    src_x = np.array([S_anc[k]["x"] for k in keys])
    dst_x = np.array([B_anc[k]["x"] for k in keys])
    m = np.polyfit(src_x, dst_x, 1)
    c1, c0 = float(m[0]), float(m[1])
    
    resid = {}
    for k in keys:
        pred = c0 + c1 * S_anc[k]["x"]   # правка Claude: в генерации c0 и c1 были переставлены
        resid[k] = float(B_anc[k]["x"] - pred)
        
    return {
        "h": [c0, c1, 0.0],
        "keys": keys,
        "resid": resid,
    }


def h_by_energy(scale_s, cal_b, n, e_lo=20.0, e_hi=3000.0, info=None):
    """Строит отображение каналов пробы на каналы фона через энергии."""
    x = np.arange(n, dtype=float)
    E = poly_eval(scale_s, x)
    xb = np.asarray(poly_inverse(list(cal_b), E, hi=n + 200), float)
    ok = np.isfinite(xb) & (xb > 0.5) & (xb < n - 1.5) & (E >= e_lo) & (E <= e_hi)   # правка Claude: только область линий; вне шкалы фона обращение не годно
    if ok.sum() < 3:
        raise ValueError("перевод фона на шкалу пробы через энергии: шкалы пробы и фона не перекрываются")
    m = np.polyfit(x[ok], xb[ok], 2)
    if info is not None:   # правка Claude: погрешность квадратичного h против точного перевода через энергии
        info["h_err_max"] = float(np.max(np.abs(np.polyval(m, x[ok]) - xb[ok])))
    return [float(m[2]), float(m[1]), float(m[0])]


def run_sample_cal(bg_path, sample_path, repers, use_intrinsic=True, verbose=True):
    """Выполняет полную калибровку по пробе."""
    rb.check_cancel("калибровка по пробе")
    S, B, k, live_s, live_b = reg.load_inputs(bg_path, sample_path)
    cal_s = list(load_spectrum(sample_path)["cal"])
    cal_b = list(load_spectrum(bg_path)["cal"])
    n = len(B)

    # 2. Собственные структуры
    intr, S_anc = intrinsic_repers(S, cal_s)
    intr = check_intrinsic(S_anc, intr, repers, S, cal_s, verbose)   # правка Claude: ложные опоры в сильной пробе
    allr = list(repers)
    if use_intrinsic:
        existing_names = set(r["name"] for r in allr)
        for r in intr:
            if r["name"] not in existing_names:
                allr.append(r)
    
    if verbose:
        print("[SCAL] реперов всего %d, из них видимых в пробе собственных структур LaBr3 %d" % (len(allr), len(intr)))

    # 3. Шкала пробы
    # правка Claude: основная шкала -- реперы от 100 кэВ; 37 кэВ и прочие ниже 100 -- только в ветвь ниже 100 кэВ (как в паре)
    sc = scale_from_repers([r for r in allr if r["E"] >= 100.0], S, cal_s, verbose)

    # 4. Совмещение фона
    rb.check_cancel("совмещение фона")
    B_anc = anchors_of(B, cal_b)
    reg_v = register_by_visible(S_anc, B_anc)
    
    if reg_v is not None:
        h = reg_v["h"]
        mode = "visible"
        if verbose:
            resid_str = ", ".join("%s=%.2f" % (k, v) for k, v in reg_v["resid"].items())
            print("[SCAL] фон совмещён по видимым собственным структурам %s: h = %s, невязки %s кан" % (
                reg_v["keys"], h, resid_str))
    else:
        scale_for_h = cal_s if sc.get("rejected") else sc["coeffs"]
        from . import pipeline as pl   # правка Claude: «своя шкала» фона -- по его собственным структурам LaBr3
        own_b = pl.run_ecal(bg_path, bg_path, [0.0, 1.0, 0.0], low=False, verbose=False)
        cal_b_own, how = bg_own_scale(own_b, cal_b, B)
        if verbose:
            print("[SCAL] шкала фона для перевода через энергии: %s %s" % (how, cal_b_own))
        hinfo = {}
        h = h_by_energy(scale_for_h, cal_b_own, n, info=hinfo)
        if verbose:
            print("[SCAL] квадратичный h против точного перевода через энергии (20-3000 кэВ): макс. отклонение %.2f кан%s" % (
                hinfo["h_err_max"], " -- ПРЕДУПРЕЖДЕНИЕ: шкалы файлов плохо согласуются квадратичным h" if hinfo["h_err_max"] > 2.0 else ""))
        mode = "energy"
        if verbose:
            print("[SCAL] видимых собственных структур в пробе меньше двух -- фон переводится на шкалу пробы через энергии (каждый файл по своей шкале)")

    # 5. Реперы на каналах фона
    good_repers, _ = valid_repers(allr)
    repers_bg = []
    for r in good_repers:
        ch_bg = poly_eval(h, r["ch"])
        deriv = poly_deriv(h, r["ch"])
        ch_err_bg = r["ch_err"] * abs(deriv)
        fwhm_bg = r["fwhm"] * abs(deriv) if np.isfinite(r["fwhm"]) else float("nan")
        
        new_r = dict(r)
        new_r["ch"] = ch_bg
        new_r["ch_err"] = ch_err_bg
        new_r["fwhm"] = fwhm_bg
        repers_bg.append(new_r)

    scb = scale_from_repers([r for r in repers_bg if r["E"] >= 100.0], B, cal_b, verbose)
    low = None if scb.get("rejected") else low_branch(scb["coeffs"], repers_bg, B, verbose)

    # 6. Точки ПШПВ
    fwhm_points = []
    for r in repers_bg:
        if np.isfinite(r["fwhm"]) and not any("ПШПВ" in w for w in r.get("warn", [])):   # правка Claude: без подозрительных ПШПВ
            if r["source"] == SRC_LABR:
                fwhm_points.append((r["ch"], r["fwhm"]))
            else:
                if np.isfinite(r["signif"]) and r["signif"] >= SIGNIF_MIN:
                    fwhm_points.append((r["ch"], r["fwhm"]))
    
    fwhm_points.sort(key=lambda p: p[0])
    
    # Убрать дубли по каналу ближе 1 кан
    if fwhm_points:
        unique_points = [fwhm_points[0]]
        for p in fwhm_points[1:]:
            if abs(p[0] - unique_points[-1][0]) >= 1.0:
                unique_points.append(p)
        fwhm_points = unique_points

    if len(fwhm_points) < 3:
        ap_res = ap.autopreset(B, S, cal_b, verbose=False)
        ref = ap.load_ref()
        extra_points = ap.fwhm_points_for(ap_res, ref["fwhm_points"])
        fwhm_points.extend([(p[0], p[1]) for p in extra_points])
        fwhm_points.sort(key=lambda p: p[0])
        # Повторно убрать дубли
        if fwhm_points:
            unique_points = [fwhm_points[0]]
            for p in fwhm_points[1:]:
                if abs(p[0] - unique_points[-1][0]) >= 1.0:
                    unique_points.append(p)
            fwhm_points = unique_points

    from . import fwhm as F   # правка Claude: модель ПШПВ должна строиться (CheckCalibration) -- иначе точки фона
    try:
        F.fit_fwhm(n, [p[0] for p in fwhm_points], [p[1] for p in fwhm_points])
    except ValueError as exc:
        ref_pts = [(p[0], p[1]) for p in ap.fwhm_points_for(ap.autopreset(B, B, cal_b, verbose=False), ap.load_ref()["fwhm_points"])]
        if verbose:
            print("[SCAL] ПРЕДУПРЕЖДЕНИЕ: точки ПШПВ пробы не дают годной модели (%s) -- взяты опорные точки ПШПВ фона %s"
                  % (str(exc).split(";")[0], ["%.1f:%.2f" % q for q in ref_pts]))
        fwhm_points = ref_pts

    # 7. Результат
    return {
        "mode": mode,
        "h": h,
        "visible_sample": [k for k in S_anc if S_anc[k]["visible"]],
        "anchors_sample": S_anc,
        "anchors_bg": B_anc,
        "reg_visible": reg_v,
        "scale_sample": sc,
        "coeffs": scb["coeffs"],
        "coeffs_xml": low["xml"]["coeffs"] if low else scb["coeffs"],
        **({"low": low} if low else {}),
        "points": scb["points"],
        "chi2ndf": scb["chi2ndf"],
        "rms": scb.get("rms"),
        "loo_rms": scb.get("loo_rms"),
        "rejected": scb.get("rejected", ""),
        "fwhm_points": [[float(a), float(b)] for a, b in fwhm_points],
        "repers": allr,
        "k": k,
    }


def reper_from_choice(fit, E, nuclide=None, tol=0.05):
    """Правка Claude: источник энергии по выбору в диалоге -- линия библиотеки (нуклид и E совпали с таблицей) или «ручная»."""
    if nuclide:
        lines = pr.library_lines((nuclide,), min_int=0.0, e_min=0.0, e_max=1e5)
        hit = [l for l in lines if abs(l["E"] - float(E)) <= tol]
        if hit:
            return make_reper(fit, hit[0]["E"], SRC_LIB, nuclide=nuclide)
    return make_reper(fit, E, SRC_MANUAL, nuclide=nuclide)


def fwhm_guess(x, n_ch):
    """Правка Claude: ожидаемая ПШПВ (кан) в канале x по опорным точкам ПШПВ пресета, пересчёт на число каналов."""
    ref = ap.load_ref()
    r0 = n_ch / ref["n_ch"]
    pts = sorted((p[0] * r0, p[1] * r0) for p in ref["fwhm_points"])
    xs, fs = [p[0] for p in pts], [p[1] for p in pts]
    f2 = np.polyfit(xs, np.square(fs), 1)   # F^2 = a + b*x
    return float(math.sqrt(max(np.polyval(f2, float(x)), 1.0)))


def _autopreset_anchors(counts, cal):
    """Правка Claude: опоры автопресета; если сумм-пик не найден или β-край с ним не согласован (прогноз по отношению к пику 37 кэВ не учитывает
    смещение нуля), окно поиска сдвигается по сетке отношений 0.75..1.25 и берётся самый значимый годный пик."""
    counts = np.asarray(counts, float)
    a = ap.autopreset(counts, counts, list(cal), verbose=False)["anchors"]
    if (a.get("sum", {}).get("ok") and a.get("edge", {}).get("ok")) or not a.get("a37", {}).get("ok"):
        return a
    ref = ap.load_ref()
    cands = [ap.find_sum(counts, a["a37"]["x"] * r, ref, None) for r in np.arange(0.75, 1.2501, 0.08)]
    # кандидат принимается, только если по паре (37, сумм-пик) β-край найден там, где его предсказывает отображение
    # эталона; из согласованных -- ближайший к прогнозу чистого усиления (смещение нуля АЦП мало), а не самый
    # значимый: сильная линия пробы значимее сумм-пика
    x_gain = a["a37"]["x"] * ref["anchors"]["sum"]["x"] / ref["anchors"]["a37"]["x"]
    for s in sorted((c for c in cands if c.get("ok")), key=lambda c: abs(c["x"] - x_gain)):
        b = ap.autopreset(counts, counts, list(cal), manual={"sum": s["x"]}, verbose=False)["anchors"]
        if b.get("edge", {}).get("ok"):
            b["sum"] = dict(s)
            return b
    return a


def low_branch(C, repers_bg, B, verbose=True):
    """Правка Claude: ветвь ниже 100 кэВ по правилам пары (lowscale.report) -- пик 37 кэВ пробы и реперы < 100 кэВ."""
    from . import lowscale as lw
    lowe = json.loads((HERE / "presets" / "labr3_lowE.json").read_text(encoding="utf-8"))
    pts = [{"name": r["name"], "E": r["E"], "dE": r["dE"], "ch": r["ch"], "ch_err": max(r["ch_err"], 0.01),
            "kind": "intrinsic37" if r["name"] == "La138_37_intr" else "photon"} for r in repers_bg if r["E"] < 100.0]
    n = len(B)
    thr = next((i + 1 for i in range(n) if B[i] > 0), 1)
    try:
        if not pts:
            raise ValueError("реперов ниже 100 кэВ нет")
        lr = lw.report(list(C), pts, n, thr, e37=lowe["e37_gamma"], de37=lowe["de37"], verbose=verbose)
    except Exception as exc:
        if verbose:
            print("[LOWE] ветвь ниже 100 кэВ не построена: %s" % exc)
        return None
    if not np.isfinite(lr["fit"]["sigma_s"]):
        if verbose:
            print("[LOWE] ветвь ниже 100 кэВ не принята: s без погрешности (реперов мало)")
        return None
    return {"s": lr["fit"]["s"], "sigma_s": lr["fit"]["sigma_s"], "shape": lr["fit"]["shape"], "dK": lr["fit"]["dK"],
            "e37": lowe["e37_gamma"], "points": pts, "residuals": lr["residuals"], "loo": lr["loo"],
            "e_threshold": lr["e_threshold"], "xml": lr["xml"], "dev_main": lr["dev_main"],
            "warn_extrapolation": lr["warn_extrapolation"],
            "band": {k: list(np.asarray(lr["band"][k])[::4]) for k in ("ch", "E", "sigma")}}


INTR_TOL = {"a37": 15.0, "edge": 0.03 * 788.74, "sum": 0.02 * 1471.0}   # кэВ: допуск опоры против шкалы по реперам пробы


def check_intrinsic(S_anc, intr, user_repers, counts, cal, verbose=True):
    """Правка Claude: «видимая» в сильной пробе опора LaBr3 может быть линией пробы (напр. рентген 75-93 кэВ вместо 37);
    при двух и более реперах пробы от 100 кэВ опора, чья энергия по их шкале дальше допуска INTR_TOL, отбрасывается."""
    user = valid_repers([r for r in user_repers if r["E"] >= 100.0])[0]
    if len(user) < 2:
        return intr
    pre = scale_from_repers(user, counts, cal, verbose=False)
    if pre.get("rejected"):
        return intr
    en = intrinsic_energies()
    for k, (E, dE, name) in en.items():
        if S_anc[k]["visible"]:
            Ex = float(poly_eval(pre["coeffs"], S_anc[k]["x"]))
            if abs(Ex - E) > INTR_TOL[k]:
                S_anc[k]["visible"] = False
                S_anc[k]["why"] = "по шкале реперов пробы %.1f кэВ вместо %.1f" % (Ex, E)
                if verbose:
                    print("[SCAL] опора %s в пробе (канал %.2f) отброшена: %s -- это линия пробы" % (k, S_anc[k]["x"], S_anc[k]["why"]))
    keep = {en[k][2] for k in en if S_anc[k]["visible"]}
    return [r for r in intr if r["name"] in keep]


BG_SPAN_MIN = 3.0   # отношение крайних энергий точек фона, при котором поправка калибровки файла применяется


def bg_own_scale(own_b, cal_b, B):
    """Правка Claude: «своя шкала» фона -- по его собственным структурам; если полином по ним отвергнут, а точек >= 2 --
    калибровка файла с линейной поправкой по энергии E' = a + b*E_файла (кривизна файла сохраняется)."""
    from . import pipeline as pl
    if not own_b.get("rejected"):
        return list(own_b["coeffs"]), "по собственным структурам фона"
    pts_b = [(q["ch"], q["E"]) for q in own_b.get("points", []) if np.isfinite(q["ch"])]
    Es = [q[1] for q in pts_b]
    if len(pts_b) >= 2 and max(Es) >= BG_SPAN_MIN * min(Es):   # узкий охват (789 и 1471) -- поправка экстраполирует вниз хуже файла
        ab = np.polyfit(poly_eval(list(cal_b), np.array([q[0] for q in pts_b])), [q[1] for q in pts_b], 1)
        C2 = [float(ab[1] + ab[0] * cal_b[0])] + [float(ab[0] * c) for c in cal_b[1:]]
        if not pl._scale_defect({"coeffs": C2, "chi2ndf": 0.0}, {"n": len(B), "B": B}):
            return C2, "калибровка файла с поправкой E' = %.2f + %.5f*E по собственным структурам фона (%s кэВ)" % (
                ab[1], ab[0], ", ".join("%.1f" % q[1] for q in pts_b))
    return list(cal_b), "калибровка файла"
