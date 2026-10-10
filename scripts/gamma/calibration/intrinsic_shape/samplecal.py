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
from . import samplecal_q as q   # правка Claude (круг 1): качество шкалы, ранжирование линий, двойной фит

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
    # правка Claude (круг 1): без ступеньки -- у слабых пиков на наклонном континууме её свободная высота сдвигает центроид
    # на 4-6 кан (медальон, 1764 кэВ: 3686.6 со ступенькой, 3681.1 без, независимая подгонка 3680.6)
    r = fit_peak(np.asarray(counts, float), float(x_click), float(fwhm_ch), half_window=hw, step=False)

    if r is None:
        return None
    if not np.isfinite(r["centroid"]) or not np.isfinite(r["centroid_err"]) or r["centroid_err"] <= 0:
        return None
    if not (r.get("area", 1.0) > 0 and np.isfinite(r["fwhm"]) and r["centroid_err"] <= q.CENTROID_ERR_MAX * r["fwhm"]
            and abs(r["centroid"] - x_click) <= hw and 0.0 <= r["centroid"] <= len(counts) - 1 and r["signif"] >= q.SIGNIF_REFUSE):
        return None   # правка Claude: подгонка выродилась -- площадь <= 0, центроид не определён (ошибка >= ПШПВ) или вне окна

    warn = []
    if r["signif"] < SIGNIF_MIN:
        warn.append("незначимый пик (значимость %.1f < %.0f)" % (r["signif"], SIGNIF_MIN))
    if fwhm_expected is not None:
        ratio = r["fwhm"] / fwhm_expected
        if ratio > q.MULTIPLET_RATIO:   # правка Claude (круг 1): широкий пик -- вероятно мультиплет
            warn.append("ПШПВ %.2f кан больше ожидаемой %.2f в %.1f раза -- вероятно мультиплет, предложен двойной фит" % (r["fwhm"], fwhm_expected, ratio))
        elif ratio < 1.0 / FWHM_RATIO_MAX:
            warn.append("ПШПВ %.2f кан отличается от ожидаемой %.2f больше чем в %.1f раза" % (r["fwhm"], fwhm_expected, FWHM_RATIO_MAX))

    return {
        "ch": float(r["centroid"]),
        "ch_err": float(r["centroid_err"]),
        "fwhm": float(r["fwhm"]),
        "fwhm_err": float(r.get("fwhm_err", float("nan"))),
        "signif": float(r["signif"]),
        "warn": warn,
    }


def candidate_lines(E_est, nuclides=None, min_int=0.5, n_max=15, e_min=10.0, e_max=4000.0, fwhm_keV=None, visible_series=()):
    """Кандидаты линий: близость x интенсивность x ряды, видимые в пробе; чужие нуклиды ниже; дублеты одного нуклида --
    эффективная энергия (правка Claude, круг 1)."""
    if nuclides is None:
        nuclides = tuple(pr.SAMPLE_NUCLIDES) + tuple(pr.EXTRA_NUCLIDES)
    lines = pr.library_lines(nuclides, min_int=min_int, e_min=e_min, e_max=e_max)
    return q.rank_candidates(lines, E_est, fwhm_keV or max(3.0, 0.03 * E_est), visible_series, n_max)   # правка Claude (круг 1)


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
    try:   # правка Claude (круг 1): пустое окно поиска («argmax of an empty sequence») -- опоры невидимы, без падения
        a = _autopreset_anchors(counts, cal)
    except (ValueError, IndexError, KeyError) as exc:
        a = {k: {"ok": False, "reason": str(exc)} for k in ("a37", "edge", "sum")}
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

    res = q.select_degree(cfg, ctx, verbose)   # правка Claude (архитектура «на шкале пробы»): степень по LOO, без гейтов -- решение в samplecal_net
    return dict(res, dropped=dropped, n_repers=len(good))


def save_json(path, repers, device, n_ch, anchors, cal=None):
    """Сохраняет набор реперов в JSON (правка Claude, круг 1: калибровка файла, значимость и предупреждения)."""
    data = {
        "cal": [float(c) for c in cal] if cal is not None else None,
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
                "signif": float(r.get("signif", float("nan"))),
                "warn": list(r.get("warn", [])),
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
    for r in data.get("repers", []):   # правка Claude (круг 1): наборы без значимости/предупреждений -- без KeyError
        r.setdefault("signif", float("nan"))
        r.setdefault("warn", [])
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

    if data.get("cal"):   # правка Claude (круг 1): физичность усиления против калибровок файлов
        ok, ratio, g_file = q.gain_check(c1, data["cal"], cal, float(np.mean(src_x)))
    else:
        ok, ratio, g_file = q.GAIN_RANGE[0] <= c1 <= q.GAIN_RANGE[1], c1, 1.0
    if not ok:   # правка Claude (VERIFY7 п.4): согласованные опоры -- виновата калибровка файла приёмника, не опоры
        dev = np.abs(c0 + c1 * src_x - dst_x) / np.array([fwhm_guess(x, len(counts)) for x in dst_x])
        agree = (len(common_keys) >= 3 and float(dev.max()) <= 0.5) or (len(common_keys) == 2 and q.GAIN_RANGE[0] <= c1 <= q.GAIN_RANGE[1])
        tail = ("опоры согласованы между собой -- вероятно, неверна калибровка файла приёмника" if agree
                else "опоры найдены неверно")
        raise ValueError("перенос набора реперов: нефизичное усиление x' = %.2f + %.4f*x (отношение %.3f к усилению по "
                         "калибровкам файлов %.4f, допустимо %.2f-%.2f) -- %s" % (c0, c1, ratio, g_file, *q.GAIN_RANGE, tail))
    dropped = []
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
                if f is not None and abs(f["ch"] - x) <= 0.5 * fwhm_new and f["signif"] >= SIGNIF_MIN:   # правка Claude: промах > 0.5 ПШПВ или значимость < 5 -- отбросить
                    new_r = dict(r)
                    new_r["ch"] = f["ch"]
                    new_r["ch_err"] = f["ch_err"]
                    new_r["fwhm"] = f["fwhm"]
                    new_r["signif"] = f["signif"]
                    new_r["warn"] = f["warn"]
                    new_repers.append(new_r)
                else:
                    dropped.append("репер %s отброшен: на новом спектре нет пика значимостью >= 5 в 0.5 ПШПВ от прогноза (канал %.1f)%s" % (
                        r["name"], x, "" if f is None else ", подогнан в %.1f" % f["ch"]))
            else:
                new_r = dict(r)
                new_r["ch"] = x
                new_r["ch_err"] = r["ch_err"] * abs(c1)
                if np.isfinite(fwhm_orig):
                    new_r["fwhm"] = fwhm_orig * abs(c1)
                new_repers.append(new_r)

    if verbose:
        for msg in dropped:
            print("[SCAL] ПРЕДУПРЕЖДЕНИЕ: " + msg)
    return new_repers, [c0, c1], dropped


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


def reper_from_choice(fit, E, nuclide=None, tol=0.05, fwhm_keV=None):
    """Правка Claude: источник энергии по выбору в диалоге -- линия библиотеки или дублет одного нуклида (эффективная
    энергия по интенсивностям), иначе «ручная»."""
    if nuclide:
        hit = q.library_energy(nuclide, float(E), fwhm_keV or max(3.0, 0.03 * float(E)), tol)
        if hit:
            return make_reper(fit, hit[0], SRC_LIB, nuclide=nuclide, dE=hit[1])
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


INTR_TOL = {"a37": 15.0, "edge": 0.03 * 788.74, "sum": 0.02 * 1471.0}   # кэВ: допуск опоры против шкалы по реперам пробы


def check_intrinsic(S_anc, intr, user_repers, counts, cal, verbose=True):
    """Правка Claude: «видимая» в сильной пробе опора LaBr3 может быть линией пробы (напр. рентген 75-93 кэВ вместо 37);
    при двух и более реперах пробы от 100 кэВ опора, чья энергия по их шкале дальше допуска INTR_TOL, отбрасывается."""
    user = valid_repers([r for r in user_repers if r["E"] >= 100.0])[0]
    if len(user) < 2:
        return intr
    pre = scale_from_repers(user, counts, cal, verbose=False)
    # правка Claude (круг 1): отвергнутая шкала = калибровка файла -- лучшая доступная, проверка опор всё равно по ней
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


BG_EXTRA_E = (609.312, 1120.294, 1764.491, 2614.511)   # линии фона (Bi-214, Tl-208) -- доп. опоры шкалы фона; K-40 1460.8 не берётся (рядом сумм-структура LaBr3)
# правка Claude (VERIFY7 п.3): поправка шкалы фона -- больше двух опор, если они есть. Сгенерировано gen_code.py (qwen3.8:27b), 2026-10-11.

def bg_own_scale(own_b, cal_b, B):
    """Правка Claude: «своя шкала» фона -- по его собственным структурам; если полином по ним отвергнут, а точек >= 2 -- калибровка
файла с поправкой по энергии E' = g(E_файла), опоры: собственные структуры + линии фона BG_EXTRA_E, найденные со значимостью >= 5
в пределах 0.5 ПШПВ от прогноза; поправка квадратичная при >= 4 опорах, иначе линейная. Последняя опора ниже 2614.5 кэВ --
в описании предупреждение об экстраполяции. Возвращает (коэффициенты, описание)."""
    from . import pipeline as pl
    n = len(B)
    if not own_b.get("rejected"):
        C = list(own_b["coeffs"])
        Es0 = [q["E"] for q in own_b.get("points", []) if np.isfinite(q["ch"])]
        how = "по собственным структурам фона"
        if Es0 and max(Es0) < 2614.5 - 1.0:
            how += "; экстраполяция выше последней опоры %.0f кэВ" % max(Es0)
        return C, how
    pts = [(float(q["ch"]), float(q["E"])) for q in own_b.get("points", []) if np.isfinite(q["ch"])]
    if len(pts) < 2 or max(e for _, e in pts) < BG_SPAN_MIN * min(e for _, e in pts):
        return list(cal_b), "калибровка файла"
    ab = np.polyfit(poly_eval(list(cal_b), np.array([p[0] for p in pts])), [p[1] for p in pts], 1)
    C1 = [float(ab[1] + ab[0] * cal_b[0])] + [float(ab[0] * c) for c in cal_b[1:]]
    extra = []
    for E in BG_EXTRA_E:
        if any(abs(E - e) < 1.0 for _, e in pts):
            continue
        x0 = float(poly_inverse(C1, E, hi=n + 200))
        if not (20 < x0 < n - 20):
            continue
        fw = fwhm_guess(x0, n)
        f = fit_click(B, x0, fw)
        if f is not None and f["signif"] >= SIGNIF_MIN and abs(f["ch"] - x0) <= 0.5 * f["fwhm"]:
            extra.append((float(f["ch"]), E))
    allp = pts + extra
    deg = 2 if len(allp) >= 4 else 1
    g = np.polynomial.polynomial.polyfit(poly_eval(list(cal_b), np.array([p[0] for p in allp])), [p[1] for p in allp], deg)
    out = np.zeros(1)
    pw = np.ones(1)
    for a in g:
        out = np.polynomial.polynomial.polyadd(out, a * pw)
        pw = np.polynomial.polynomial.polymul(pw, np.asarray(cal_b, float))
    C2 = [float(v) for v in out]
    if pl._scale_defect({"coeffs": C2, "chi2ndf": 0.0}, {"n": n, "B": B}):
        return list(cal_b), "калибровка файла"
    Es = sorted(e for _, e in allp)
    how = "калибровка файла с поправкой степени %d по %d опорам фона (%s кэВ)" % (deg, len(allp), ", ".join("%.1f" % e for e in Es))
    if Es[-1] < 2614.5 - 1.0:
        how += "; экстраполяция выше последней опоры %.0f кэВ" % Es[-1]
    return C2, how


def anchors_for_save(counts, cal, repers, verbose=True):
    """Правка Claude (круг 1): опоры для JSON -- после проверки против шкалы по реперам пробы (ложная a37 не сохраняется)."""
    anc = anchors_of(counts, cal)
    check_intrinsic(anc, [], repers, counts, cal, verbose)
    return anc


