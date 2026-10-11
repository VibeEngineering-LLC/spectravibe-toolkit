"""
Устойчивая регистрация спектра образца на фоне по форме внутренних структур детектора.
Модуль добавляет линейные возмущающие члены (гладкая базовая линия и гауссиана K-40),
которые профилируются внутри хи-квадрат, а также итеративно исключает окна вокруг
значимых пиков образца, найденных в спектре разности.
Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-10-10
"""

import time
import numpy as np
from .core import poly_eval, rebin_conserve_var
from . import register as reg

CANCEL = {"flag": False}   # флаг отмены, ставит окно (UI-поток), проверяется между итерациями


class Cancelled(Exception):
    """Операция отменена пользователем."""


def check_cancel(where=""):
    """точка отмены между этапами (Шкала, Вычесть, Совместить)."""
    if CANCEL["flag"]:
        raise Cancelled("операция отменена пользователем" + (" (" + where + ")" if where else ""))


def hat_basis(n, lo, hi, spacing):
    """
    Построение базиса кусочно-линейных "шапок" (hat functions).
    Узлы равномерно распределены от lo до hi с шагом spacing.
    Возвращает матрицу (n, K), где K - количество узлов.
    """
    K = max(2, int(round((hi - lo) / spacing)) + 1)
    knots = np.linspace(lo, hi, K)
    x = np.arange(n, dtype=float)
    basis = np.zeros((n, K), dtype=float)

    for i in range(K):
        if i == 0:
            # Левый край: константа 1 для x <= knot0, линейное падение до 0 на [knot0, knot1]
            if K > 1:
                mask_left = x <= knots[0]
                mask_slope = (x > knots[0]) & (x <= knots[1])
                basis[mask_left, i] = 1.0
                if np.any(mask_slope):
                    basis[mask_slope, i] = (knots[1] - x[mask_slope]) / (knots[1] - knots[0])
            else:
                basis[:, i] = 1.0
        elif i == K - 1:
            # Правый край: линейный рост от 0 на [knot_{K-2}, knot_{K-1}], константа 1 для x >= knot_{K-1}
            mask_slope = (x >= knots[i - 1]) & (x < knots[i])
            mask_right = x >= knots[i]
            if np.any(mask_slope):
                basis[mask_slope, i] = (x[mask_slope] - knots[i - 1]) / (knots[i] - knots[i - 1])
            basis[mask_right, i] = 1.0
        else:
            # Внутренние узлы: треугольная форма
            mask_left = (x >= knots[i - 1]) & (x < knots[i])
            mask_right = (x >= knots[i]) & (x <= knots[i + 1])
            if np.any(mask_left):
                basis[mask_left, i] = (x[mask_left] - knots[i - 1]) / (knots[i] - knots[i - 1])
            if np.any(mask_right):
                basis[mask_right, i] = (knots[i + 1] - x[mask_right]) / (knots[i + 1] - knots[i])

    return basis


def step_col(n, x0, sigma):
    """ступень (комптоновский край линии образца) 0.5*erfc((x-x0)/(sqrt2*sigma)) -- 1 слева, 0 справа."""
    from scipy.special import erfc
    x = np.arange(n, dtype=float)
    return 0.5 * erfc((x - x0) / (np.sqrt(2.0) * sigma))


def gauss_col(n, x0, sigma):
    """
    Построение гауссовой колонки.
    Возвращает массив (n,) с значениями exp(-0.5*((x - x0)/sigma)**2).
    """
    x = np.arange(n, dtype=float)
    return np.exp(-0.5 * ((x - x0) / sigma) ** 2)


class RobustProblem(reg.Problem):
    """
    Расширение класса Problem для учета линейных возмущающих членов (базиса).
    """

    def __init__(self, S, B, k, mask, anchor, basis, bg_var=None):
        super().__init__(S, B, k, mask, anchor, bg_var)
        self.basis = None if basis is None else np.asarray(basis, dtype=float)
        self.n_nuis = 0 if basis is None else self.basis.shape[1]

    def _parts(self, c):
        """
        Вычисление разности и весов.
        Возвращает (d, w) или None, если отображение не монотонно.
        """
        e = poly_eval(c, self.x_edges)
        if np.any(np.diff(e) <= 0):
            return None

        sv, vs = rebin_conserve_var(self.S, e, len(self.B))
        bb, vb = reg.blur_bg(lambda x: poly_eval(c, x), self.B, self.bg_var)
        w = 1.0 / np.sqrt(np.maximum(vs + self.k ** 2 * vb, 1.0))
        return (sv - self.k * bb), w

    def nuisance(self, c, parts=None, need_cov=True):
        """
        Профилирование линейных возмущающих членов.
        Возвращает (coef, cov) или (None, None).
        """
        if parts is None:
            parts = self._parts(c)
        if parts is None or self.basis is None:
            return None, None

        d, w = parts
        m = self.mask
        if getattr(self, "_Bm", None) is None:   # базис на маске -- один раз, не на каждом шаге
            self._Bm = self.basis[m]
        wm = w[m]
        Aw = self._Bm * wm[:, None]
        G = Aw.T @ Aw
        rhs = Aw.T @ (d[m] * wm)
        try:   # нормальные уравнения K x K вместо SVD (n_маски x K) на каждой оценке χ²
            coef = np.linalg.solve(G, rhs)
        except np.linalg.LinAlgError:
            coef = np.linalg.lstsq(Aw, d[m] * wm, rcond=None)[0]
        cov = np.linalg.pinv(G) if need_cov else None   # ковариация нужна только итогу, не сетке
        return coef, cov

    def resid(self, c):
        """
        Вычисление взвешенного остатка с учетом возмущающих членов.
        """
        parts = self._parts(c)
        if parts is None:
            return None

        d, w = parts
        if self.basis is not None:
            coef, _ = self.nuisance(c, parts=parts, need_cov=False)
            d = d - self.basis @ coef

        return d * w

    def ndf(self, npar=3):
        """
        Количество степеней свободы с учетом возмущающих членов.
        """
        return self.nres - npar - self.n_nuis


def find_peaks_mf(resid, var, fwhm_fn, lo, hi, z_min=6.0):
    """
    Поиск пиков методом согласованного фильтра (matched filter).
    Возвращает список кортежей (channel, z).
    """
    n = len(resid)
    z_scores = np.zeros(hi - lo + 1)
    channels = np.arange(lo, hi + 1)

    for i, c in enumerate(channels):
        fw = fwhm_fn(c)
        hw = max(3, int(round(1.5 * fw)))
        a = max(0, c - hw)
        b = min(n - 1, c + hw)

        # Гауссов шаблон с нулевым средним
        x_win = np.arange(a, b + 1, dtype=float)
        sigma = fw / 2.3548
        g = np.exp(-0.5 * ((x_win - c) / sigma) ** 2)
        g -= np.mean(g)

        var_win = var[a:b + 1]
        resid_win = resid[a:b + 1]

        # Избегаем деления на ноль
        valid = var_win > 0
        if not np.any(valid):
            z_scores[i] = 0.0
            continue

        wv = (g * g) / var_win
        sum_wv = np.sum(wv)
        if sum_wv <= 0:
            z_scores[i] = 0.0
            continue

        amp = np.sum(g * resid_win / var_win) / sum_wv
        z = amp * np.sqrt(sum_wv)
        z_scores[i] = z

    # Поиск локальных максимумов
    peaks = []
    for i, c in enumerate(channels):
        if z_scores[i] <= z_min:
            continue

        fw = fwhm_fn(c)
        local_radius = max(1, int(0.5 * fw))
        start = max(0, i - local_radius)
        end = min(len(z_scores) - 1, i + local_radius)

        if z_scores[i] == np.max(z_scores[start:end + 1]):
            # диполь рассогласования (рядом провал отклика сравнимой величины) -- не линия образца
            w2 = max(2, int(2.0 * fw))
            if np.min(z_scores[max(0, i - w2):min(len(z_scores), i + w2 + 1)]) < -0.8 * z_scores[i]:
                continue
            peaks.append((int(c), float(z_scores[i])))

    return sorted(peaks, key=lambda x: x[0])


def merge_windows(centres, fwhm_fn, half_fwhm, n):
    """
    Объединение перекрывающихся окон вокруг центров пиков.
    Возвращает список окон [lo, hi].
    """
    if not centres:
        return []

    windows = []
    for c in centres:
        fw = fwhm_fn(c)
        lo = int(np.floor(c - half_fwhm * fw))
        hi = int(np.ceil(c + half_fwhm * fw))
        lo = max(0, lo)
        hi = min(n - 1, hi)
        windows.append([lo, hi])

    # Сортировка по нижнему краю
    windows.sort(key=lambda x: x[0])

    # Объединение перекрывающихся или касающихся окон
    merged = []
    for w in windows:
        if not merged:
            merged.append(w)
        else:
            last = merged[-1]
            if w[0] <= last[1] + 1:
                last[1] = max(last[1], w[1])
            else:
                merged.append(w)

    return merged


def mask_minus(mask, windows):
    """
    Создание копии маски с установленными в False каналами из окон.
    """
    new_mask = mask.copy()
    for w in windows:
        lo, hi = w
        new_mask[lo:hi + 1] = False
    return new_mask


def structure_checks(net, var, model, windows, excl=(), z_thr=3.0):
    """
    Проверка структур в окнах.
    Возвращает список словарей с результатами проверок.
    """
    if model is None:
        model = np.zeros_like(net)

    r = net - model
    results = []

    # Преобразование excl в список окон для проверки
    excl_windows = list(excl)

    for name, w in windows.items():
        lo, hi = w
        # Индексы каналов в окне, не входящие в исключенные окна
        idx = []
        for ch in range(lo, hi + 1):
            in_excl = False
            for ew in excl_windows:
                if ew[0] <= ch <= ew[1]:
                    in_excl = True
                    break
            if not in_excl:
                idx.append(ch)

        n_used = len(idx)
        if n_used < 5:
            results.append({
                "name": name,
                "lo": lo,
                "hi": hi,
                "n_used": n_used,
                "z_sum": 0.0,
                "z_dip": 0.0,
                "z_chi2": 0.0,
                "ok": False,
                "reason": "окно закрыто линиями образца"
            })
            continue

        idx_arr = np.array(idx)
        rr = r[idx_arr]
        vv = var[idx_arr]

        # Проверка суммы
        z_sum = rr.sum() / np.sqrt(vv.sum())

        # Проверка "впадины" (dip)
        centre = (lo + hi) / 2.0
        sgn = np.sign(idx_arr - centre)
        # Если все знаки одинаковы или нулевые, z_dip = 0
        if np.all(sgn == 0):
            z_dip = 0.0
        else:
            denom = np.sqrt(np.sum(sgn * sgn * vv))
            if denom > 0:
                z_dip = np.sum(sgn * rr) / denom
            else:
                z_dip = 0.0

        # Проверка хи-квадрат
        chi2 = np.sum(rr * rr / vv)
        z_chi2 = (chi2 - n_used) / np.sqrt(2 * n_used)

        ok = (abs(z_sum) <= z_thr) and (abs(z_dip) <= z_thr) and (z_chi2 <= z_thr)

        results.append({
            "name": name,
            "lo": lo,
            "hi": hi,
            "n_used": n_used,
            "z_sum": float(z_sum),
            "z_dip": float(z_dip),
            "z_chi2": float(z_chi2),
            "ok": bool(ok),
            "reason": ""
        })

    return results


def register_robust(S, B, k, mask0, anchor, fwhm_fn, k40_ch=None, mode="auto",
                    knot_ch=150.0, z_peak=6.0, excl_fwhm=2.5, protect_fwhm=1.0,
                    max_iter=6, chi2_auto=1.3, z_auto=5.0, verbose=True, chi2_stop=20.0, compton_fn=None, exclude=True, extra_fn=None, **reg_kw):
    """
    Устойчивая регистрация спектра образца на фоне.
    """
    S = np.asarray(S, dtype=float)
    B = np.asarray(B, dtype=float)
    n = len(B)

    # Определение границ маски
    true_indices = np.where(mask0)[0]
    if len(true_indices) == 0:
        raise ValueError("Маска не содержит истинных значений")
    lo0 = true_indices[0]
    hi0 = true_indices[-1]

    # Простая регистрация
    t_pl = time.time()
    plain = reg.register(S, B, k, mask0, anchor, verbose=False, **reg_kw)
    if verbose:
        print("[TIME] простая регистрация: %.1f с" % (time.time() - t_pl))
    if CANCEL["flag"]:
        raise Cancelled("операция отменена пользователем")
    net0, var0 = reg.apply(plain['c'], S, B, k)
    z_plain = net0[mask0].sum() / np.sqrt(var0[mask0].sum())

    # Принятие решения о режиме
    if mode == "off":
        robust = False
    elif mode == "on":
        robust = True
    else:  # auto
        robust = (plain['chi2ndf'] > chi2_auto) or (abs(z_plain) > z_auto)

    if verbose:
        if robust:
            msg = "устойчивая регистрация (линии образца)"
        else:
            msg = "простая (образец не добавляет структур)"
        print(f"[ROBUST] простая регистрация: chi2/ndf={plain['chi2ndf']:.4f}, z суммы разности в маске={z_plain:+.1f} -> {msg}")

    if not robust:
        return {
            "mode": "plain",
            "reg": plain,
            "c": plain['c'],
            "chi2ndf_plain": plain['chi2ndf'],
            "chi2ndf": plain['chi2ndf'],
            "z_plain": float(z_plain),
            "excluded": [],
            "peaks": [],
            "skipped": [],
            "k40": None,
            "mask": mask0.copy(),
            "basis": None,
            "coef": None,
            "model": None,
            "iterations": 0,
            "converged": True,
            "steps_ch": [],
            "same_basis": None,
            "stop_reason": ""
        }

    # Построение базиса
    # комптоновские края линий образца (K-40 и найденных пиков, compton_fn: канал пика -> канал края)
    # -- ступени в базисе (структуры образца, не собственные); K-40 -- последний столбец.
    s40 = fwhm_fn(k40_ch) / 2.3548 if k40_ch is not None else None
    steps_ch = []   # каналы ступеней последнего базиса -- для вердикта по окнам
    def make_basis(peak_chs):
        cols = [hat_basis(n, lo0, hi0, knot_ch)]
        steps_ch.clear()
        for pc in (list(peak_chs) + ([k40_ch] if k40_ch is not None else [])) if compton_fn else []:
            ce = compton_fn(pc)
            if ce is not None and lo0 < ce < hi0:
                cols.append(step_col(n, ce, fwhm_fn(ce) / 2.3548)[:, None])
                steps_ch.append(float(ce))
        if k40_ch is not None:
            cols.append(gauss_col(n, k40_ch, s40)[:, None])
        return np.column_stack(cols)
    basis = make_basis([])

    excl = []
    peaks_all = []
    skipped = []
    converged = False
    excl_used = None
    rr = None
    c = None
    prob = None
    net = None
    var = None
    coef = None
    cov = None
    model = None
    it = 0

    hist, stop_reason, chi_prev = [], "", float("inf")
    for it in range(max_iter):
        if CANCEL["flag"]:   # отмена из окна -- между итерациями
            raise Cancelled("операция отменена пользователем (итерация %d)" % it)
        t_it = time.time()
        skipped = []   # только пики последней итерации
        mask = mask_minus(mask0, excl)
        prob = RobustProblem(S, B, k, mask, anchor, basis)
        rr = reg.register(S, B, k, mask, anchor, verbose=False, problem=prob, **reg_kw)
        c = rr['c']

        net, var = reg.apply(c, S, B, k)
        coef, cov = prob.nuisance(c)
        model = basis @ coef

        # Поиск пиков в разности
        pk = find_peaks_mf(net - model, var, fwhm_fn, lo0, hi0, z_peak)

        keep = []
        for (ch, z) in pk:
            if k40_ch is not None and abs(ch - k40_ch) <= protect_fwhm * fwhm_fn(k40_ch):
                # Проверка, не добавлен ли уже
                if (ch, z) not in skipped:
                    skipped.append((ch, z))
            else:
                keep.append((ch, z))

        # extra_fn -- каналы прочих линий нуклидов, отождествлённых по найденным пикам (библиотека линий)
        extra = list(extra_fn([p[0] for p in keep])) if extra_fn else []
        if k40_ch is not None:   # зона сумм-пика La-138 не исключается и для линий библиотеки
            extra = [x for x in extra if abs(x - k40_ch) > protect_fwhm * fwhm_fn(k40_ch)]
        centres = [p[0] for p in keep] + extra
        new = merge_windows(centres, fwhm_fn, excl_fwhm, n) if exclude else []
        basis = make_basis([p[0] for p in keep])

        if verbose:
            print(f"[ROBUST] итерация {it}: chi2/ndf={rr['chi2ndf']:.4f}, найдено пиков образца {len(keep)}, исключено окон {len(new)}")
            for (ch, z) in keep:
                print(f"[ROBUST]   пик кан {ch} z={z:.1f}")

        if verbose:
            print("[TIME] ROBUST итерация %d: %.1f с" % (it, time.time() - t_it))
        if new == excl and basis.shape[1] == prob.basis.shape[1]:
            peaks_all = keep
            converged = True
            excl_used = excl
            break
        elif it >= 1 and (new in hist or (rr['chi2ndf'] > chi2_stop and rr['chi2ndf'] >= chi_prev)):   # досрочная остановка
            stop_reason = ("набор исключённых окон повторяется (колебание)" if new in hist else
                           "χ²/ndf %.1f > %.0f и не падает (было %.1f) на итерации %d" % (rr['chi2ndf'], chi2_stop, chi_prev, it))
            basis, peaks_all, excl_used = prob.basis, keep, excl   # итог -- последняя выполненная регистрация
            if verbose:
                print("[ROBUST] ПРЕДУПРЕЖДЕНИЕ: не сошлась: образец с множеством линий по всей шкале (%s) -- "
                      "остановлено досрочно на итерации %d из %d" % (stop_reason, it, max_iter))
            break
        else:
            hist.append(excl)
            chi_prev = rr['chi2ndf']
            excl = new
            peaks_all = keep
            excl_used = excl

    # Финальная регистрация, если не была выполнена в последней итерации с тем же excl
    if excl_used != excl or basis.shape[1] != prob.basis.shape[1]:
        mask = mask_minus(mask0, excl)
        prob = RobustProblem(S, B, k, mask, anchor, basis)
        rr = reg.register(S, B, k, mask, anchor, verbose=False, problem=prob, **reg_kw)
        c = rr['c']
        net, var = reg.apply(c, S, B, k)
        coef, cov = prob.nuisance(c)
        model = basis @ coef
        it += 1  # Учет дополнительной итерации

    # Информация о K-40
    k40info = None
    if k40_ch is not None:
        amp = float(coef[-1])
        amp_err = float(np.sqrt(cov[-1, -1] * max(1, rr['chi2ndf'])))
        area = amp * s40 * np.sqrt(2 * np.pi)
        area_err = amp_err * s40 * np.sqrt(2 * np.pi)
        k40info = {
            "ch": k40_ch,
            "sigma": float(s40),
            "amp": amp,
            "amp_err": amp_err,
            "area": float(area),
            "area_err": float(area_err)
        }

    # сравнение простой и устойчивой на ОДНОМ наборе каналов (итоговая маска) и ОДНОЙ модели
    m_fin = mask_minus(mask0, excl)
    p_nom = reg.Problem(S, B, k, m_fin, anchor)
    p_mod = RobustProblem(S, B, k, m_fin, anchor, basis)
    same = {"n_ch": int(m_fin.sum()),
            "no_model": [p_nom.chi2(plain['c']) / p_nom.ndf(), p_nom.chi2(c) / p_nom.ndf()],
            "with_model": [p_mod.chi2(plain['c']) / p_mod.ndf(), p_mod.chi2(c) / p_mod.ndf()]}
    if verbose:
        print("[ROBUST] на одних каналах (%d, итоговая маска): без модели образца chi2/ndf простая %.3f -> устойчивая %.3f; "
              "с моделью образца %.3f -> %.3f" % (same["n_ch"], *same["no_model"], *same["with_model"]))
        if not converged:
            print("[ROBUST] ПРЕДУПРЕЖДЕНИЕ: исключение линий не сошлось за %d итераций -- итог не устойчив" % max_iter)
    if verbose:
        print(f"[ROBUST] итог: chi2/ndf {rr['chi2ndf']:.4f} (простая {plain['chi2ndf']:.4f}), окон исключено {len(excl)}, каналов {np.sum(mask_minus(mask0, excl))} из {n}, сходимость {converged}")
        for w in excl:
            print(f"[ROBUST] исключено кан {w[0]}-{w[1]}")
        if k40info:
            print(f"[ROBUST] K-40 под сумм-пиком La-138: кан {k40_ch:.1f}, площадь {k40info['area']:.0f} ± {k40info['area_err']:.0f} отсчётов")
        for (ch, z) in skipped:
            print(f"[ROBUST] пик кан {ch} z={z:.1f} у K-40 -- не исключён (моделируется гауссианой)")

    return {
        "mode": "robust",
        "reg": rr,
        "c": c,
        "chi2ndf_plain": plain['chi2ndf'],
        "chi2ndf": rr['chi2ndf'],
        "z_plain": float(z_plain),
        "excluded": excl,
        "peaks": peaks_all,
        "skipped": skipped,
        "k40": k40info,
        "mask": mask_minus(mask0, excl),
        "basis": basis,
        "coef": coef,
        "model": model,
        "iterations": it + 1,
        "converged": converged,
        "plain": plain,
        "steps_ch": list(steps_ch),
        "same_basis": same,
        "stop_reason": stop_reason
    }
