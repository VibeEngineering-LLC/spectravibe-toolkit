"""Модуль самокалибровки спектра детектора NaI(Tl) по линиям естественного фона."""

import numpy as np
from scipy.ndimage import gaussian_filter1d

from gamma.calibration.intrinsic_shape.core import poly_eval, poly_deriv, poly_inverse
from gamma.calibration.intrinsic_shape.peakfit import fit_peak
from gamma.calibration.intrinsic_shape.fwhm import fit_fwhm, fwhm as fwhm_eval


def find_strong(counts, file_cal, e_strong, fwhm_rel662, search_rel):
    """Найти сильную опорную линию и вернуть (отношение каналов, результат подгонки).

    Parameters
    ----------
    counts : array_like
        Массив счётов (float).
    file_cal : sequence
        Коэффициенты шкалы файла (по возрастанию степени).
    e_strong : float
        Энергия опорной сильной линии, кэВ.
    fwhm_rel662 : float
        Относительная ПШПВ при 661.7 кэВ.
    search_rel : float
        Относительная полуширина окна поиска.

    Returns
    -------
    tuple
        (ratio, result_dict), где ratio = измеренный канал / канал по шкале файла.
    """
    n = len(counts)
    x_pred = poly_inverse(file_cal, e_strong, hi=n + 200)
    slope = poly_deriv(file_cal, x_pred)  # кэВ/кан
    fw = fwhm_rel662 * np.sqrt(661.7 * e_strong) / slope  # ПШПВ в каналах
    sg = fw / 2.3548

    # Фильтр разности гауссиан (сглаживание «шаром»)
    d = (gaussian_filter1d(counts, sg, mode="nearest")
         - gaussian_filter1d(counts, 3.0 * sg, mode="nearest"))

    lo = max(0, int(x_pred * (1.0 - search_rel)))
    hi = min(n - 1, int(x_pred * (1.0 + search_rel)))
    xm = lo + int(np.argmax(d[lo:hi + 1]))

    r = fit_peak(counts, float(xm), fw, step=False)
    if r is None:
        raise ValueError(f"опорная линия {e_strong:.1f} кэВ не найдена")

    return (r["centroid"] / x_pred, r)


def calibrate(counts, file_cal, lines, e_strong=1460.8, fwhm_rel662=0.07,
              search_rel=0.05, min_z=5.0, deg=2):
    """Самокалибровка спектра по известным линиям естественного фона.

    Parameters
    ----------
    counts : array_like
        Массив счётов (float).
    file_cal : sequence
        Коэффициенты шкалы файла (по возрастанию степени).
    lines : list of dict
        Список опорных линий: ``{"E": кэВ, "sys_kev": число}``.
    e_strong : float, optional
        Энергия сильной опорной линии (по умолчанию 1460.8 кэВ, K-40).
    fwhm_rel662 : float, optional
        Относительная ПШПВ при 661.7 кэВ.
    search_rel : float, optional
        Относительная полуширина окна поиска.
    min_z : float, optional
        Минимальная значимость пика.
    deg : int, optional
        Степень полинома шкалы E(канал).

    Returns
    -------
    dict
        Результат калибровки (ecal, ratio, lines, n_accepted, resid_kev, fwc, R662, n).
    """
    c = np.array(counts, dtype=np.float64, copy=True)
    c[-1] = 0.0
    n = len(c)

    # Определяем отношение «фактический канал / канал по шкале файла»
    ratio, _ = find_strong(c, file_cal, e_strong, fwhm_rel662, search_rel)

    # Обход всех заданных линий
    records = []
    for ln in lines:
        E = ln["E"]
        sys_kev = ln.get("sys_kev", 0.0)

        xf = poly_inverse(file_cal, E, hi=n + 200)
        x0 = xf * ratio
        slope = poly_deriv(file_cal, xf) / ratio  # кэВ/кан в окрестности линии
        fw0 = fwhm_rel662 * np.sqrt(661.7 * E) / slope  # ожидаемая ПШПВ, кан

        # Проверка попадания в спектр
        if x0 < 3.0 or x0 > n - 3.0:
            records.append({
                "E": float(E), "sys_kev": float(sys_kev),
                "x0": float(x0), "ch": None, "ch_err": None,
                "fwhm_ch": None, "fwhm_err": None,
                "signif": None, "slope": float(slope),
                "accepted": False, "reason": "вне шкалы",
            })
            continue

        r = fit_peak(c, float(x0), float(fw0), step=False)

        # Каскад критериев принятия
        if r is None:
            accepted, reason = False, "нет подгонки"
        elif r["signif"] < min_z:
            accepted, reason = False, "z<min"
        elif abs(r["centroid"] - x0) > 0.5 * fw0:
            accepted, reason = False, "сдвиг>0.5ПШПВ"
        elif not (0.5 * fw0 <= r["fwhm"] <= 2.0 * fw0):
            accepted, reason = False, "ширина"
        elif r["centroid_err"] >= 0.2 * fw0:
            accepted, reason = False, "ошибка центроида"
        else:
            accepted, reason = True, ""

        if accepted:
            records.append({
                "E": float(E), "sys_kev": float(sys_kev),
                "x0": float(x0), "ch": float(r["centroid"]),
                "ch_err": float(r["centroid_err"]),
                "fwhm_ch": float(r["fwhm"]), "fwhm_err": float(r["fwhm_err"]),
                "signif": float(r["signif"]), "slope": float(slope),
                "accepted": True, "reason": reason,
            })
        else:
            records.append({
                "E": float(E), "sys_kev": float(sys_kev),
                "x0": float(x0), "ch": None, "ch_err": None,
                "fwhm_ch": None, "fwhm_err": None,
                "signif": None, "slope": float(slope),
                "accepted": False, "reason": reason,
            })

    # Отбираем принятые линии
    acc = [rec for rec in records if rec["accepted"]]
    n_acc = len(acc)
    if n_acc < 2:
        raise ValueError(f"мало опорных линий фона: {n_acc}")

    # --- Шкала E(канал) ---
    dg = min(deg, n_acc - 1)
    ch_arr = np.array([rec["ch"] for rec in acc], dtype=np.float64)
    E_arr = np.array([rec["E"] for rec in acc], dtype=np.float64)

    # Стандартная неопределённость энергии линии
    sig_E = np.array([
        np.sqrt((rec["slope"] * rec["ch_err"]) ** 2 + rec["sys_kev"] ** 2 + 0.01)
        for rec in acc
    ], dtype=np.float64)

    # Взвешенная полиномиальная подгонка (numpy.polyfit: коэффициенты по убыванию)
    p = np.polyfit(ch_arr, E_arr, dg, w=1.0 / sig_E)
    ecal = list(reversed(p))  # по возрастанию степени
    while len(ecal) < 3:
        ecal.append(0.0)

    # Невязки в кэВ
    resid_kev = [float(E_arr[i] - poly_eval(ecal, ch_arr[i])) for i in range(n_acc)]

    # --- Модель ПШПВ ---
    fwc = None
    if n_acc >= 3:
        fwhm_arr = np.array([rec["fwhm_ch"] for rec in acc], dtype=np.float64)
        try:
            fwc = fit_fwhm(n, x=ch_arr, f=fwhm_arr)
        except ValueError:
            fwc = None

    # --- Относительная ПШПВ при 661.7 кэВ ---
    if fwc is not None:
        x662 = poly_inverse(ecal, 661.7, hi=n + 200)
        R662 = float(fwhm_eval(fwc, x662) * poly_deriv(ecal, x662) / 661.7)
    else:
        vals = np.array([
            rec["fwhm_ch"] * rec["slope"] / np.sqrt(661.7 * rec["E"])
            for rec in acc
        ], dtype=np.float64)
        R662 = float(np.median(vals))

    return {
        "ecal": ecal,
        "ratio_vs_file": float(ratio),
        "lines": records,
        "n_accepted": n_acc,
        "resid_kev": resid_kev,
        "fwc": fwc,
        "R662": R662,
        "n": n,
    }


def fwhm_ch_fn(cal):
    """Возвращает функцию расчёта ПШПВ в каналах по результату :func:`calibrate`.

    Parameters
    ----------
    cal : dict
        Словарь, возвращаемый :func:`calibrate`.

    Returns
    -------
    callable
        Функция ``f(x)``: x — канал (скаляр или массив) → ПШПВ в каналах,
        не менее 1.0.
    """
    fwc = cal["fwc"]
    ecal = cal["ecal"]
    R662 = cal["R662"]

    def f(x):
        if fwc is not None:
            y = fwhm_eval(fwc, x)
        else:
            E = np.maximum(poly_eval(ecal, x), 1.0)
            slope = np.maximum(poly_deriv(ecal, x), 1e-6)
            y = R662 * np.sqrt(661.7 * E) / slope
        y = np.maximum(y, 1.0)
        if np.ndim(x) == 0:
            return float(y)
        return y

    return f


def energy_to_ch(cal, E):
    """Преобразовать энергию (кэВ) в номер канала по шкале ``cal``.

    Parameters
    ----------
    cal : dict
        Словарь, возвращаемый :func:`calibrate`.
    E : float or array_like
        Энергия в кэВ.

    Returns
    -------
    float or ndarray
        Номер канала.
    """
    return poly_inverse(cal["ecal"], E, hi=cal["n"] + 200)
