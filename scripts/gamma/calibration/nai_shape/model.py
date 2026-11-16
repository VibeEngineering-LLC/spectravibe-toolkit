"""Модуль совмещения спектров NaI с учётом дрейфа усиления/смещения и уширения."""

import numpy as np
from scipy.optimize import least_squares
from scipy.special import ndtr

from gamma.calibration.intrinsic_shape.core import poly_eval, rebin_conserve_var
from gamma.calibration.intrinsic_shape.register import blur_bg


def blur_var(y, var, sigma):
    """Сохраняющее счёт размытие с переменной шириной по каналу.

    Parameters
    ----------
    y : ndarray, shape (n,)
        Исходные значения.
    var : ndarray, shape (n,)
        Вариации.
    sigma : ndarray, shape (n,)
        σ (в каналах) для каждого исходного канала j.

    Returns
    -------
    (out, vout) : tuple of ndarray
    """
    n = len(y)
    sigma = np.asarray(sigma, dtype=float)
    if np.max(sigma) < 1e-3:
        return y.copy(), var.copy()

    K = int(np.ceil(4.0 * np.max(sigma)))
    sg = np.maximum(sigma, 1e-3)
    offsets = np.arange(-K, K + 1)  # (2K+1,)

    # Веса: (2K+1, n)
    w = np.exp(-0.5 * (offsets[:, None] / sg[None, :]) ** 2)
    w = w / w.sum(axis=0, keepdims=True)

    out = np.zeros(n)
    vout = np.zeros(n)

    # Целевые индексы: (2K+1, n)
    target = np.arange(n)[None, :] + offsets[:, None]
    contrib_out = w * y[None, :]
    contrib_var = w ** 2 * var[None, :]

    valid = (target >= 0) & (target < n)
    np.add.at(out, target[valid], contrib_out[valid])
    np.add.at(vout, target[valid], contrib_var[valid])
    return out, vout


def coefs(p, h0, n_b):
    """Вычисляет 3 коэффициента отображения h по параметрам p = [d, gN, q, s]."""
    g = 1.0 + p[1] / n_b
    h0 = np.asarray(h0, dtype=float)
    if h0.size < 3:
        h0 = np.pad(h0, (0, 3 - h0.size))
    c = g * h0[:3].copy()
    c[0] += p[0]
    c[2] += p[2] / n_b ** 2   # ручная правка: q -- добавка в каналах НА канале n_b (было p[2]/n_b -> q*n_b)
    return c


class Problem:
    """Задача совмещения спектра образца S со спектром фона B."""

    def __init__(self, S, B, k, mask, h0, basis=None, ref_ch=None, blur_eq=True):
        self.S = np.asarray(S, dtype=float)
        self.B = np.asarray(B, dtype=float)
        self.k = float(k)
        self.mask = np.asarray(mask, dtype=bool)
        self.h0 = np.asarray(h0, dtype=float)
        self.basis = basis
        self.n_s = len(self.S)
        self.n_b = len(self.B)
        self.m = 0 if basis is None else basis.shape[1]
        self.nmask = int(self.mask.sum())
        self.ref_ch = self.n_b / 2.0 if ref_ch is None else float(ref_ch)
        self.blur_eq = blur_eq

    def evaluate(self, p):
        """Вычисляет все промежуточные величины для параметров p.

        Returns
        -------
        dict or None
        """
        c = coefs(p, self.h0, self.n_b)
        edges = poly_eval(c, np.arange(self.n_s + 1) - 0.5)
        if not np.all(np.diff(edges) > 0.0):
            return None

        sh, vs = rebin_conserve_var(self.S, edges, self.n_b)

        if self.blur_eq:
            bb, vb = blur_bg(lambda x: poly_eval(c, x), self.B, self.B)
        else:
            bb, vb = self.B.copy(), self.B.copy()

        s = float(p[3])
        sig = abs(s) * np.sqrt(np.maximum(np.arange(self.n_b), 0) / self.ref_ch)
        if s > 0:
            bb, vb = blur_var(bb, vb, sig)
        elif s < 0:
            sh, vs = blur_var(sh, vs, sig)

        r0 = sh - self.k * bb
        var = np.maximum(vs + self.k ** 2 * vb, 1.0)
        w = 1.0 / np.sqrt(var)

        if self.basis is not None:
            b_m = self.basis[self.mask]
            r0_m = r0[self.mask]
            w_m = w[self.mask]
            amps = np.linalg.lstsq(b_m * w_m[:, None], r0_m * w_m, rcond=None)[0]
            fit_lin = self.basis @ amps
        else:
            amps = np.empty(0)
            fit_lin = np.zeros(self.n_b)

        z = (r0 - fit_lin) * w
        return {"c": c, "sh": sh, "vs": vs, "bb": bb, "vb": vb, "net": r0,
                "var": var, "amps": amps, "fit_lin": fit_lin, "z": z}

    def resid(self, p):
        """Вектор взвешенных остатков на каналах маски."""
        ev = self.evaluate(p)
        if ev is None:
            return np.full(self.nmask, 1e3)
        return ev["z"][self.mask]

    def chi2(self, p):
        """χ² = сумма квадратов resid."""
        r = self.resid(p)
        return float(np.sum(r ** 2))

    def ndf(self, n_free):
        """Число степеней свободы."""
        return self.nmask - n_free - self.m


def fit(prob, free, p_init=None, d_grid=None, g_grid=None, n_starts=3):
    """Многофазовая оптимизация свободных параметров free.

    Parameters
    ----------
    prob : Problem
    free : list of int
        Индексы свободных параметров (0=d, 1=gN, 2=q, 3=s).
    p_init : ndarray or None
        Начальный вектор параметров (длина 4).
    d_grid, g_grid : ndarray or None
        Сетки для глобального поиска.
    n_starts : int
        Число стартовых точек.

    Returns
    -------
    dict
    """
    if p_init is None:
        p_init = np.zeros(4)
    p_init = np.asarray(p_init, dtype=float).copy()

    # --- Поиск стартовых точек ---
    if d_grid is not None and g_grid is not None:
        d_grid = np.asarray(d_grid, dtype=float)
        g_grid = np.asarray(g_grid, dtype=float)
        nd, ng = len(d_grid), len(g_grid)
        G = np.zeros((nd, ng))
        for i in range(nd):
            for j in range(ng):
                p_t = p_init.copy()
                p_t[0] = d_grid[i]
                p_t[1] = g_grid[j]
                G[i, j] = prob.chi2(p_t)

        # Локальные минимумы в окрестности 3×3
        candidates = []
        for i in range(nd):
            for j in range(ng):
                i0, i1 = max(0, i - 1), min(nd - 1, i + 1)
                j0, j1 = max(0, j - 1), min(ng - 1, j + 1)
                window = G[i0:i1 + 1, j0:j1 + 1]
                if G[i, j] <= window.min():
                    candidates.append((G[i, j], i, j))
        candidates.sort(key=lambda x: x[0])
        candidates = candidates[:n_starts]

        start_ps = []
        for _, i, j in candidates:
            p_s = p_init.copy()
            p_s[0] = d_grid[i]
            p_s[1] = g_grid[j]
            start_ps.append(p_s)
    else:
        G = None
        start_ps = [p_init.copy()]

    # --- Локальная оптимизация из каждой стартовой точки ---
    results = []
    for p_s in start_ps:
        def fun(x):
            p_local = p_s.copy()
            p_local[list(free)] = x
            return prob.resid(p_local)

        x0 = p_s[list(free)]
        res = least_squares(fun, x0, method="trf", diff_step=1e-3,
                            x_scale=1.0, max_nfev=400)
        p_best = p_s.copy()
        p_best[list(free)] = res.x
        chi2_val = prob.chi2(p_best)
        results.append((chi2_val, p_best))

    results.sort(key=lambda x: x[0])
    best_chi2, best_p = results[0]

    return {
        "p": best_p,
        "chi2": float(best_chi2),
        "ndf": prob.ndf(len(free)),
        "minima": [{"p": list(pp), "chi2": float(c)} for c, pp in results],
        "grid": G,
    }


def covariance(prob, p, free, step=0.01):
    """Ковариационная матрица свободных параметров (центральные разности)."""
    p = np.asarray(p, dtype=float)
    n_free = len(free)
    n_res = prob.nmask

    J = np.zeros((n_res, n_free))
    for i, idx in enumerate(free):
        p_plus = p.copy()
        p_minus = p.copy()
        p_plus[idx] += step
        p_minus[idx] -= step
        J[:, i] = (prob.resid(p_plus) - prob.resid(p_minus)) / (2.0 * step)

    chi2 = prob.chi2(p)
    ndf = prob.ndf(n_free)
    scale = max(1.0, chi2 / ndf)
    cov = np.linalg.pinv(J.T @ J) * scale

    return {
        "cov": cov,
        "err": np.sqrt(np.diag(cov)),
        "chi2ndf": chi2 / ndf,
        "scale": scale,
        "J": J,
    }


def h_sigma(p, cov, free, h0, n_b, chans):
    """Неопределённость отображения h (в каналах) в точках chans.

    sigma_h = sqrt(diag(D @ cov @ D.T)), где D — производные h по свободным
    параметрам (численно, центральная разность, шаг 0.01).
    """
    p = np.asarray(p, dtype=float)
    chans = np.asarray(chans, dtype=float)
    step = 0.01
    n_free = len(free)

    D = np.zeros((len(chans), n_free))
    for i, idx in enumerate(free):
        p_plus = p.copy()
        p_minus = p.copy()
        p_plus[idx] += step
        p_minus[idx] -= step
        h_plus = poly_eval(coefs(p_plus, h0, n_b), chans)
        h_minus = poly_eval(coefs(p_minus, h0, n_b), chans)
        D[:, i] = (h_plus - h_minus) / (2.0 * step)

    # diag(D @ cov @ D.T)
    sigma_h = np.sqrt(np.einsum('ij,jk,ik->i', D, cov, D))
    return sigma_h


def hats(n, knots):
    """Базис кусочно-линейных «шапок» по узловым точкам knots.

    Returns
    -------
    ndarray, shape (n, len(knots))
    """
    knots = np.asarray(knots, dtype=float)
    m = len(knots)
    x = np.arange(n)
    out = np.zeros((n, m))
    for i in range(m):
        e_i = np.zeros(m)
        e_i[i] = 1.0
        out[:, i] = np.interp(x, knots, e_i)
    return out


def gauss_cols(n, centers, sigmas):
    """Столбцы гауссиан единичной площади на сетке каналов 0..n-1.

    Интеграл по каналу j: Φ((j+0.5−c)/σ) − Φ((j−0.5−c)/σ).

    Returns
    -------
    ndarray, shape (n, len(centers))
    """
    centers = np.asarray(centers, dtype=float)
    sigmas = np.asarray(sigmas, dtype=float)
    x = np.arange(n)
    lo = x - 0.5
    hi = x + 0.5
    out = (ndtr((hi[:, None] - centers[None, :]) / sigmas[None, :])
           - ndtr((lo[:, None] - centers[None, :]) / sigmas[None, :]))
    return out
