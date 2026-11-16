"""Проверки разности в окнах."""

import numpy as np

__all__ = ["window_tests", "verdict", "masked_chi2"]


def window_tests(net, var, e2ch, fwhm_kev, lines_kev,
                 half_fwhm=1.0, exclude_kev=(), exclude_fwhm=1.5,
                 sideband=False, tags=None):
    """Посчитать z-статистики и chi²/ndf для каждого окна.

    Параметры
    ----------
    net : array_like
        Разность спектров (сигнал − фон), по каналам.
    var : array_like
        Дисперсия разности по каналам.
    e2ch : callable
        Функция энергия (кэВ) → канал (float).
    fwhm_kev : callable
        Функция энергия (кэВ) → ПШПВ (кэВ).
    lines_kev : sequence of float
        Энергии (кэВ), для которых проверяются окна.
    half_fwhm : float, optional
        Полуширина окна в единицах ПШПВ (по умолчанию 1.0).
    exclude_kev : sequence of float, optional
        Энергии линий образца для исключения соседних окон.
    exclude_fwhm : float, optional
        Порог исключения в единицах ПШПВ (по умолчанию 1.5).
    sideband : bool, optional
        Включить оценку и вычитание базовой линии по боковым полосам.
    tags : dict или None, optional
        Словарь {энергия: метка} для подписей.

    Возврат
    -------
    list of dict
        Записи по каждому окну.
    """
    net = np.asarray(net, dtype=np.float64)
    var = np.asarray(var, dtype=np.float64)
    nch = len(net)
    records = []

    for E in lines_kev:
        xc = float(e2ch(E))
        fw = float(e2ch(E + fwhm_kev(E) / 2.0) - e2ch(E - fwhm_kev(E) / 2.0))
        lo = int(round(xc - half_fwhm * fw))
        hi = int(round(xc + half_fwhm * fw))
        lo = max(0, lo)
        hi = min(nch - 1, hi)

        # Проверка исключения из-за близости к линии образца
        excluded = any(abs(E - Ex) < exclude_fwhm * fwhm_kev(E) for Ex in exclude_kev)

        x = np.arange(lo, hi + 1, dtype=np.float64)
        r = net[lo:hi + 1].copy()
        v = var[lo:hi + 1].copy()
        var_base = 0.0
        sideband_ok = False

        # Оценка и вычитание базовой линии по боковым полосам
        if sideband:
            w = max(3, int(round(fw)))
            L_lo = max(0, lo - w)
            L_hi = min(nch - 1, lo - 1)
            R_lo = max(0, hi + 1)
            R_hi = min(nch - 1, hi + w)

            sb_x = np.concatenate([
                np.arange(L_lo, L_hi + 1, dtype=np.float64),
                np.arange(R_lo, R_hi + 1, dtype=np.float64),
            ])
            sb_i = sb_x.astype(np.int64)   # ручная правка: индексы целые (в генерации был float -> IndexError)
            sb_y = net[sb_i].copy()
            sb_v = var[sb_i].copy()

            if len(sb_x) >= 4 and np.all(sb_v > 0):
                sideband_ok = True
                # Нормированные переменные: модель a + b*(x - xc)
                t = sb_x - xc
                W = 1.0 / sb_v
                # Матрица A: строки [1, t_i]
                # N = A^T W A (2×2), rhs = A^T W y
                S = np.array([
                    [np.sum(W), np.sum(W * t)],
                    [np.sum(W * t), np.sum(W * t * t)],
                ])
                rhs = np.array([np.sum(W * sb_y), np.sum(W * t * sb_y)])
                S_inv = np.linalg.inv(S)
                ab = S_inv @ rhs
                a_fit, b_fit = ab[0], ab[1]
                # Ковариация параметров = S_inv
                # Вычитание базовой линии из окна
                r -= a_fit + b_fit * (x - xc)
                # Дисперсия суммы базовой линии по окну: u^T C u
                u = np.array([len(x), np.sum(x - xc)])
                var_base = float(u @ S_inv @ u)

        # Статистики
        v_total = v.sum() + var_base
        z_sum = float(r.sum() / np.sqrt(v_total)) if v_total > 0 else 0.0

        sgn = np.sign(x - xc)
        denom_dip = np.sqrt(np.sum(sgn ** 2 * v))
        z_dip = float(np.sum(sgn * r) / denom_dip) if denom_dip > 0 else 0.0

        chi2ndf = float(np.sum(r * r / v) / len(r)) if len(r) > 0 else 0.0

        rec = {
            "E": float(E),
            "tag": (tags.get(E, "") if tags is not None else ""),
            "lo": int(lo),
            "hi": int(hi),
            "excluded": bool(excluded),
            "sum": float(r.sum()),
            "sigma": float(np.sqrt(v_total)),
            "z_sum": z_sum,
            "z_dip": z_dip,
            "chi2ndf": chi2ndf,
            "sideband_ok": bool(sideband_ok),
        }
        records.append(rec)

    return records


def verdict(tests, z_thr=3.0, use_sum_for=None):
    """Свести результаты window_tests в итоговый вердикт.

    Параметры
    ----------
    tests : list of dict
        Записи из :func:`window_tests`.
    z_thr : float, optional
        Порог по |z| (по умолчанию 3.0).
    use_sum_for : set или None, optional
        Множество энергий, для которых дополнительно проверяется z_sum.
        Если None — z_sum проверяется для всех.

    Возврат
    -------
    dict
        Сводная информация: число окон, максимумы, список «плохих» E.
    """
    considered = [t for t in tests if not t["excluded"]]
    n = len(considered)
    n_ok = 0
    bad = []
    max_abs_z_sum = 0.0
    max_abs_z_dip = 0.0

    for t in considered:
        abs_zs = abs(t["z_sum"])
        abs_zd = abs(t["z_dip"])
        max_abs_z_sum = max(max_abs_z_sum, abs_zs)
        max_abs_z_dip = max(max_abs_z_dip, abs_zd)

        ok = abs_zd <= z_thr
        if use_sum_for is None:
            ok = ok and (abs_zs <= z_thr)
        else:
            if t["E"] in use_sum_for:
                ok = ok and (abs_zs <= z_thr)

        if ok:
            n_ok += 1
        else:
            bad.append(t["E"])

    return {
        "n": int(n),
        "n_ok": int(n_ok),
        "max_abs_z_sum": float(max_abs_z_sum),
        "max_abs_z_dip": float(max_abs_z_dip),
        "bad": bad,
        "ok": bool(n_ok == n and n > 0),
    }


def masked_chi2(z, mask):
    """Chi² по подмножеству каналов, заданному маской.

    Параметры
    ----------
    z : array_like
        Значения z-статистик (r / sigma) по каналам.
    mask : array_like of bool
        Маска выбранных каналов.

    Возврат
    -------
    dict
        {"chi2": float, "n": int, "chi2_per_ch": float}
    """
    z = np.asarray(z, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool)
    z_sel = z[mask]
    chi2 = float(np.sum(z_sel ** 2))
    n = int(mask.sum())
    chi2_per_ch = chi2 / n if n > 0 else 0.0
    return {"chi2": chi2, "n": n, "chi2_per_ch": chi2_per_ch}
