"""
Модуль парной калибровки: реперы и точки ПШПВ, извлечённые из самой пары
(фон + образец), а не из пресета, сделанного под один конкретный образец.
Фон: внутренние структуры. Образец: гамма-линии в спектре разности net,
идентифицированные по энергиям из библиотечных данных проекта.
Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-10-10
"""

import csv
import copy
import numpy as np
from scipy.optimize import curve_fit
from gamma.data import DATA_DIR
from gamma.data.nuclide_library import get_nuclide
from .core import poly_eval, poly_inverse, poly_deriv
from .peakfit import fit_peak

# правка Claude: по умолчанию -- природные ряды, K-40 и Cs-137 (типичный образец среды). Калибровочные нуклиды
# (Am-241, Ba-133, Co-60, Cs-134, Eu-152) в общем списке блокировали «изолированность» природных линий
# (Cs-134 604.7 рядом с Bi-214 609) и давали ложные отождествления; их можно передать явно.
SAMPLE_NUCLIDES = ("Cs-137", "K-40", "Pb-214", "Bi-214", "Tl-208", "Ac-228", "Pb-212")
NET_Z_NEG = 3.0   # правка Claude: z суммы разности <= -3 -- «значимо отрицательна», линии образца в ней не ищутся
EXTRA_NUCLIDES = ("Ra-226", "Am-241", "Ba-133", "Co-60", "Cs-134", "Eu-152")
AVOID_KEV = ((30.0, 45.0), (505.0, 517.0), (1425.0, 1485.0))   # пик ~37 кэВ, аннигиляция/Tl-208 510.77, сумм-пик La-138 (K-40 там моделируется отдельно)
ME_KEV = 510.999


def xray_blend(element, lines):
    """
    Читает CSV с X-лучами и возвращает (E, dE) для заданного элемента и набора линий.
    E - интенсивностно-взвешенное среднее энергий, dE - максимальное из 0.05 и
    интенсивностно-взвешенного стандартного отклонения энергий.
    """
    path = DATA_DIR / "xray_lines_lsrm.csv"
    selected = []
    with open(path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["element"] == element and row["line"] in lines:
                selected.append(row)

    if not selected:
        raise ValueError("нет линий %s %s в %s" % (element, lines, path))

    energies = []
    intensities = []
    for row in selected:
        e = float(row["energy_keV"])
        i = float(row["intensity"])
        energies.append(e)
        intensities.append(i)

    energies = np.array(energies)
    intensities = np.array(intensities)

    total_int = np.sum(intensities)
    if total_int <= 0:
        raise ValueError("нет линий %s %s в %s" % (element, lines, path))

    weights = intensities / total_int
    E = float(np.sum(weights * energies))
    var = float(np.sum(weights * (energies - E) ** 2))
    dE = max(0.05, float(np.sqrt(var)))

    return E, dE


def library_lines(nuclides=SAMPLE_NUCLIDES, min_int=3.0, e_min=40.0, e_max=3000.0):
    """
    Возвращает список линий для всех изотопов в nuclides с интенсивностью >= min_int
    и энергией в диапазоне [e_min, e_max]. Отсортировано по энергии.
    """
    result = []
    for name in nuclides:
        nuclide = get_nuclide(name)
        if nuclide is None:
            continue
        lines = nuclide.get("lines", [])
        for line in lines:
            E = float(line[0])
            I = float(line[1])
            if I >= min_int and e_min <= E <= e_max:
                result.append({"nuclide": name, "E": E, "I": I})

    result.sort(key=lambda x: x["E"])
    return result


def all_lines(nuclides=SAMPLE_NUCLIDES):
    """
    Возвращает ВСЕ линии для всех изотопов (без отсечения по интенсивности),
    с энергией в диапазоне 20..3000 кэВ. Отсортировано по энергии.
    """
    return library_lines(nuclides, min_int=0.0, e_min=20.0, e_max=3000.0)


def isolated(lines, fwhm_kev, others, min_sep=1.0, frac=0.1, avoid=AVOID_KEV):
    """
    Отбирает изолированные линии: энергия не в запрещённых окнах и нет других линий
    с близкой энергией и достаточной интенсивностью.
    """
    result = []
    for L in lines:
        E_L = L["E"]
        I_L = L["I"]
        nuclide_L = L["nuclide"]

        # Проверка на запрещённые окна
        in_avoid = False
        for lo, hi in avoid:
            if lo <= E_L <= hi:
                in_avoid = True
                break
        if in_avoid:
            continue

        # Проверка на близость других линий
        fw = fwhm_kev(E_L)
        sep_threshold = min_sep * fw
        blocked = False
        for O in others:
            if O["nuclide"] == nuclide_L and abs(O["E"] - E_L) < 1e-9:
                continue  # та же линия
            if abs(O["E"] - E_L) < sep_threshold and O["I"] >= frac * I_L:
                blocked = True
                break
        if blocked:
            continue

        result.append(L)

    return result


def fit_net_peak(y, var, x0, fwhm0, hw):
    """
    Взвешенная подгонка гауссиана + линейный фон в окне вокруг x0.
    Возвращает словарь с параметрами или None при неудаче.
    """
    n = len(y)
    lo = int(np.clip(round(x0) - hw, 0, n - 1))
    hi = int(np.clip(round(x0) + hw, 0, n - 1))

    if hi - lo + 1 < 7:
        return None

    x = np.arange(lo, hi + 1)
    y_win = y[lo:hi + 1]
    var_win = var[lo:hi + 1]
    sigma = np.sqrt(np.maximum(var_win, 1.0))

    # Начальные параметры
    edge_vals = np.concatenate([y_win[:3], y_win[-3:]])
    H0 = max(float(np.max(y_win) - np.median(edge_vals)), 1.0)
    c0 = float(x0)
    s0 = fwhm0 / 2.3548
    a0_0 = float(np.mean([np.mean(y_win[:3]), np.mean(y_win[-3:])]))
    a1_0 = 0.0

    def model(x, H, c, s, a0, a1):
        return H * np.exp(-0.5 * ((x - c) / s) ** 2) + a0 + a1 * (x - x0)

    p0 = [H0, c0, s0, a0_0, a1_0]
    bounds = ([0, x0 - 0.5 * fwhm0 - 2, 0.3, -np.inf, -np.inf],
              [np.inf, x0 + 0.5 * fwhm0 + 2, 3 * fwhm0 / 2.3548, np.inf, np.inf])

    try:
        popt, pcov = curve_fit(model, x, y_win, p0=p0, sigma=sigma,
                               bounds=bounds, absolute_sigma=True, maxfev=20000)
    except Exception:
        return None

    H, c, s, a0, a1 = popt
    H = float(H)
    c = float(c)
    s = float(s)

    # Расчёт chi2 и ndf
    y_pred = model(x, *popt)
    residuals = (y_win - y_pred) / sigma
    chi2 = float(np.sum(residuals ** 2))
    ndf = len(x) - len(popt)
    chi2ndf = chi2 / max(ndf, 1)
    scale = float(np.sqrt(max(1.0, chi2ndf)))

    # Ошибки параметров
    diag_pcov = np.diag(pcov)
    if np.any(np.isnan(diag_pcov)):
        return None
    errors = np.sqrt(diag_pcov) * scale

    H_err = float(errors[0])
    c_err = float(errors[1])
    s_err = float(errors[2])

    # Площадь и её ошибка
    sqrt2pi = float(np.sqrt(2 * np.pi))
    area = H * s * sqrt2pi
    var_area = (s * sqrt2pi) ** 2 * H_err ** 2 + (H * sqrt2pi) ** 2 * s_err ** 2
    # Ковариация H и s из pcov
    cov_Hs = float(pcov[0, 2]) * scale ** 2
    var_area += 2 * (s * sqrt2pi) * (H * sqrt2pi) * cov_Hs
    area_err = float(np.sqrt(max(0.0, var_area)))

    fwhm = 2.3548 * s
    fwhm_err = 2.3548 * s_err

    signif = area / area_err if area_err > 0 else 0.0

    return {
        "centroid": c,
        "centroid_err": c_err,
        "fwhm": fwhm,
        "fwhm_err": fwhm_err,
        "area": area,
        "area_err": area_err,
        "signif": signif,
        "chi2ndf": chi2ndf
    }


def prov_scale(cal, ch_sum, e_sum=1471.0):
    """
    Временная энергетическая шкала: калибровка cal умножается на f = e_sum / poly_eval(cal, ch_sum),
    чтобы сумм-пик La-138 попал в e_sum.
    """
    E_at_ch_sum = poly_eval(cal, ch_sum)
    if E_at_ch_sum == 0:
        raise ValueError("poly_eval(cal, ch_sum) == 0")
    f = e_sum / E_at_ch_sum
    return [float(c * f) for c in cal]


def find_net_lines(net, var, scale, fwhm_ch, nuclides=SAMPLE_NUCLIDES, z_min=8.0,
                   min_int=3.0, e_min=40.0, e_max=3000.0, verbose=True, tol_rel=0.03, min_sep=1.0):
    """
    Ищет гамма-линии в спектре разности net, используя библиотеку изотопов.
    Возвращает список принятых линий.
    """
    n = len(net)

    def ch_of(E):
        return poly_inverse(scale, E, hi=n + 200)

    def fwhm_kev(E):
        ch = ch_of(E)
        return fwhm_ch(ch) * poly_deriv(scale, ch)

    cand = isolated(library_lines(nuclides, min_int, e_min, e_max), fwhm_kev, all_lines(nuclides), min_sep=min_sep)

    accepted = []
    for L in cand:
        E = L["E"]
        I = L["I"]
        nuclide = L["nuclide"]

        x0 = ch_of(E)
        if x0 < 5 or x0 > n - 10:
            if verbose:
                print("[LINE] %-8s %8.2f кэВ: пропущена (канал вне диапазона)" % (nuclide, E))
            continue

        fw = fwhm_ch(x0)
        hw = max(5, int(round(1.6 * fw)))
        # правка Claude: предварительная шкала может ошибаться на проценты (калибровка файла) -- сначала максимум
        # разности в окне +-max(1.5 ПШПВ, tol_rel*x0), подгонка -- от него (иначе центроид упирается в границу)
        dx = max(1.5 * fw, tol_rel * x0)
        a_, b_ = int(max(0, x0 - dx)), int(min(n - 1, x0 + dx))
        xs = a_ + int(np.argmax(np.convolve(net[a_:b_ + 1], np.ones(3) / 3.0, mode="same"))) if b_ > a_ + 2 else x0
        r = fit_net_peak(net, var, float(xs), fw, hw)

        reason = None
        if r is None:
            reason = "подгонка не сошлась"
        elif r["signif"] < z_min:
            reason = "значимость < %g" % z_min
        elif abs(r["centroid"] - x0) > 0.5 * fw + 2 + tol_rel * x0:
            reason = "центроид далеко"
        elif not (0.6 * fw <= r["fwhm"] <= 1.7 * fw):
            reason = "ширина вне 0.6..1.7"
        elif r["centroid_err"] >= 0.25 * fw:
            reason = "ошибка центроида велика"

        if reason is None:
            record = {
                "nuclide": nuclide,
                "E": float(E),
                "I": float(I),
                "ch": float(r["centroid"]),
                "ch_err": float(r["centroid_err"]),
                "fwhm": float(r["fwhm"]),
                "fwhm_err": float(r["fwhm_err"]),
                "signif": float(r["signif"]),
                "x_pred": float(x0),
                "E_prov": float(poly_eval(scale, r["centroid"]))
            }
            accepted.append(record)
            if verbose:
                print("[LINE] %-8s %8.2f кэВ: кан %.2f ± %.2f (ожид. %.1f), ПШПВ %.2f, значимость %.1f %s" %
                      (nuclide, E, r["centroid"], r["centroid_err"], x0, r["fwhm"], r["signif"], "принята"))
        else:
            if verbose:
                print("[LINE] %-8s %8.2f кэВ: кан %.2f ± %.2f (ожид. %.1f), ПШПВ %.2f, значимость %.1f %s" %
                      (nuclide, E, r["centroid"] if r else x0, r["centroid_err"] if r else 0.0, x0,
                       r["fwhm"] if r else fw, r["signif"] if r else 0.0, "отклонена: " + reason))

    accepted.sort(key=lambda x: x["E"])
    return accepted


def is_lu_sample(net, var, ecal_cfg, h, fwhm_ch, z_min=10.0):
    """
    Проверяет, является ли образец Lu-176, ищя линии 202 и 307 кэВ.
    """
    lu_points = []
    for p in ecal_cfg.get("points", []):   # правка Claude: конфиг -- словарь points/checks
        if p.get("name") in ("Lu176_202", "Lu176_307") and p.get("spectrum") == "sample":
            lu_points.append(p)

    if len(lu_points) < 2:
        return False

    for p in lu_points:
        x0 = p["x0"]
        x_bg = poly_eval(h, x0)
        fw = fwhm_ch(x_bg)
        hw = max(5, int(round(1.6 * fw)))
        r = fit_net_peak(net, var, x_bg, fw, hw)
        if r is None or r["signif"] < z_min:
            return False

    return True


def pair_ecal_cfg(preset_cfg, B, net, var, h, scale, fwhm_ch, lu=None, verbose=True, use_lines=True, **line_kw):
    """
    Формирует конфигурацию калибровки для пары (фон + образец).
    Возвращает (cfg, lines).
    """
    if lu is None:
        lu = is_lu_sample(net, var, preset_cfg, h, fwhm_ch)

    if lu:
        cfg = copy.deepcopy(preset_cfg)
        if verbose:
            print("[PAIR] образец с Lu-176: реперы образца -- из пресета")
        return cfg, []

    cfg = copy.deepcopy(preset_cfg)

    # Оставляем только точки с spectrum == "bg"
    bg_points = [p for p in cfg["points"] if p.get("spectrum") == "bg"]   # правка Claude

    # Для точек типа "peak" проверяем наличие пика
    kept_bg = []
    for p in bg_points:
        if p.get("kind") == "peak":
            r = fit_peak(B, p["x0"], p["w0"], half_window=int(p["hw"]))
            if r is None or r["signif"] < 5:
                if verbose:
                    print("[PAIR] репер фона %s отброшен: нет пика (значимость %.1f)" % (p.get("name", "?"), r["signif"] if r else 0.0))
                continue
        kept_bg.append(p)

    # Оставляем только проверки с spectrum == "bg"
    checks = [c for c in cfg.get("checks", []) if c.get("spectrum") == "bg"]

    # Ищем линии в образце
    # правка Claude: линии образца ищутся в разности, только если она не значимо отрицательна (фон из другой
    # обстановки / перепутанные спектры: образец ниже фона·k -- «линии» разности ложные); собственные реперы -- из фона
    nn = len(net)
    a, b = int(60 * nn / 2048), int(1700 * nn / 2048) + 1
    z_net = float(np.sum(net[a:b]) / np.sqrt(max(np.sum(var[a:b]), 1.0)))
    if not use_lines:   # правка Claude (LABR-11 п.0): регистрация не сошлась -- положения линий в разности смещены
        lines = []
        if verbose:
            print("[PAIR] устойчивая регистрация не сошлась: линии образца из разности в реперы не берутся (их положения на шкале фона смещены)")
    elif z_net <= -NET_Z_NEG:
        lines = []
        if verbose:
            print("[PAIR] разность значимо отрицательна (z суммы %.1f в каналах %d-%d): образец ниже фона·k, линии "
                  "образца не ищутся -- реперы только собственные структуры фона" % (z_net, a, b - 1))
    else:
        lines = find_net_lines(net, var, scale, fwhm_ch, verbose=verbose, **line_kw)

    # Добавляем линии как реперы
    new_points = []
    for L in lines:
        E = L["E"]
        # Проверяем, нет ли уже репера с близкой энергией
        already = False
        for p in kept_bg:
            if abs(p.get("E", 0) - E) < 2.0:
                already = True
                break
        if already:
            if verbose:
                print("[PAIR] линия %s %.2f уже есть среди реперов фона" % (L["nuclide"], E))
            continue

        name = "%s_%d" % (L["nuclide"].replace("-", ""), round(E))
        point = {
            "name": name,
            "spectrum": "net",
            "kind": "given",
            "ch": L["ch"],
            "ch_err": max(L["ch_err"], 0.02),
            "E": E,
            "dE": 0.05,
            "x0": L["ch"],
            "w0": L["fwhm"],
            "hw": 0
        }
        new_points.append(point)

    # Формируем итоговую конфигурацию
    cfg["points"] = kept_bg + new_points
    cfg["checks"] = checks
    cfg["_pair"] = True

    if verbose:
        print("[PAIR] реперов шкалы: фон %d, линии образца %d" % (len(kept_bg), len(new_points)))

    return cfg, lines


def pair_fwhm_points(anchors, lines, preset_points, lu, min_points=3, verbose=True):
    """
    Формирует список точек (ch, fwhm) для калибровки ПШПВ.
    """
    if lu:
        return [tuple(p) for p in preset_points]

    points = [
        (anchors["a37"]["x"], anchors["a37"]["fwhm"]),
        (anchors["sum"]["x"], anchors["sum"]["fwhm"])
    ]
    for L in lines:
        if L["signif"] >= 10 and L["fwhm_err"] <= 0.2 * L["fwhm"]:
            points.append((L["ch"], L["fwhm"]))

    points.sort(key=lambda x: x[0])
    points = [(round(float(ch), 2), round(float(fw), 3)) for ch, fw in points]

    if len(points) < min_points:
        if verbose:
            print("[PAIR] ПРЕДУПРЕЖДЕНИЕ: точек ПШПВ пары %d < %d -- взяты точки пресета" % (len(points), min_points))
        return [tuple(p) for p in preset_points]

    if verbose:
        print("[PAIR] точки ПШПВ пары: " + str(points))

    return points


def compton_edge_kev(E):
    """
    Энергия комптоновского края для гамма-кванта энергии E (кэВ).
    """
    return E - E / (1 + 2 * E / ME_KEV)


def chain_centres(peak_chs, scale, fwhm_ch, n, nuclides=SAMPLE_NUCLIDES, min_int=2.0, id_frac=0.5, tol_fwhm=0.5):
    """Правка Claude: нуклид отождествлён, если найденный пик совпал (|E - E_lib| < tol_fwhm*ПШПВ) с одной из его
    сильнейших линий (I >= id_frac * I_max); при нескольких кандидатах берётся первый по порядку списка nuclides.
    Возвращает (каналы ВСЕХ линий отождествлённых нуклидов с I >= min_int %, список нуклидов)."""
    def ch_of(E):
        return float(poly_inverse(scale, E, hi=n + 200))
    strong = {}
    for nm in nuclides:
        ls = library_lines([nm], 0.0, 20.0, 3000.0)
        if ls:
            imax = max(L["I"] for L in ls)
            strong[nm] = [L for L in ls if L["I"] >= id_frac * imax]
    found = []
    for pc in peak_chs:
        E = float(poly_eval(scale, pc))
        fk = fwhm_ch(pc) * float(poly_deriv(scale, pc))
        for nm in nuclides:
            if any(abs(L["E"] - E) < tol_fwhm * fk for L in strong.get(nm, [])):
                if nm not in found:
                    found.append(nm)
                break
    chs = [ch_of(L["E"]) for L in library_lines(found, min_int, 20.0, 3000.0)]
    return [c for c in chs if 0 < c < n - 1], sorted(found)
