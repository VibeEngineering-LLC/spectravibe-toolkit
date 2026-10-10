"""Общие функции метода калибровки по форме собственных структур (пространство КАНАЛОВ).

Чтение -- читателем toolkit (gamma.io.readers.read_spectrum), без своей копии.
Здесь: полиномы (коэффициенты по возрастанию степени), перебинировка с сохранением
счёта и честной дисперсией, качество вычитания в окне.
Перенесено из прототипа labr3_selfcal (common.py, subtract.py), 2026-10-10.
"""
import pathlib
import numpy as np
from scipy.ndimage import gaussian_filter1d

from gamma.io.readers import read_spectrum

N_CH = 2048   # число каналов, на котором метод проверен; функции берут длину из данных


def load_spectrum(path, n_ch=None) -> dict:
    """counts дополняется нулями до числа каналов файла (читатель может отбросить канал переполнения)."""
    s = read_spectrum(str(path))
    raw = np.asarray(s.counts, dtype=float)
    n = int(n_ch or getattr(s, "n_channels_raw", 0) or len(raw))
    counts = np.zeros(n, dtype=np.float64)
    n_copy = min(len(raw), n)
    if n_copy > 0:
        counts[:n_copy] = raw[:n_copy]
    return {
        "counts": counts,
        "live": float(s.live_time),
        "real": float(s.real_time),
        "cal": tuple(float(c) for c in s.energy_cal),
        "start": str(getattr(s, "start_datetime", "")),
        "name": pathlib.Path(path).name,
        "path": str(path),
        "n_read": int(len(raw)),
    }


def load_pair(bg_path, sample_path, zero_last=True) -> tuple:
    """Фон и образец (образец -- на число каналов фона); последний канал-переполнение обнуляется."""
    bg = load_spectrum(bg_path)
    sample = load_spectrum(sample_path, n_ch=len(bg["counts"]))
    if zero_last:
        bg["counts"][-1] = 0.0
        sample["counts"][-1] = 0.0
    return (bg, sample)


def poly_eval(coeffs, x):
    x_arr = np.asarray(x, dtype=float)
    scalar = x_arr.ndim == 0
    result = np.zeros_like(x_arr, dtype=float)
    for c in reversed(coeffs):
        result = result * x_arr + c
    if scalar:
        return float(result)
    return result


def poly_deriv(coeffs, x):
    if len(coeffs) <= 1:
        x_arr = np.asarray(x, dtype=float)
        scalar = x_arr.ndim == 0
        if scalar:
            return 0.0
        return np.zeros_like(x_arr, dtype=float)
    deriv_coeffs = [i * coeffs[i] for i in range(1, len(coeffs))]
    return poly_eval(deriv_coeffs, x)


def poly_inverse(coeffs, E, lo=-0.5, hi=N_CH - 0.5 + 200, n_grid=40001):
    g = np.linspace(lo, hi, n_grid)
    y = poly_eval(coeffs, g)
    # Правка вручную (Claude): полином файла фона 4-й степени перестаёт расти
    # около канала ~2240 — берём монотонный префикс сетки, за ним линейная экстраполяция.
    dy = np.diff(y)
    bad = np.nonzero(dy <= 0)[0]
    if bad.size:
        if bad[0] < 10:
            raise ValueError("calibration not monotonic on grid")
        g, y = g[:bad[0]], y[:bad[0]]
    E_arr = np.asarray(E, dtype=float)
    scalar = E_arr.ndim == 0
    x = np.interp(E_arr, y, g)
    s_hi = (y[-1] - y[-2]) / (g[-1] - g[-2])
    s_lo = (y[1] - y[0]) / (g[1] - g[0])
    x = np.where(E_arr > y[-1], g[-1] + (E_arr - y[-1]) / s_hi, x)
    x = np.where(E_arr < y[0], g[0] + (E_arr - y[0]) / s_lo, x)
    inside = (E_arr >= y[0]) & (E_arr <= y[-1])
    for _ in range(3):
        dx = poly_deriv(coeffs, x)
        x = np.where(inside, x - (poly_eval(coeffs, x) - E_arr) / dx, x)
    if scalar:
        return float(x)
    return x


def rebin_conserve(counts, src_edges_in_dst, n_dst=N_CH) -> np.ndarray:
    src_edges_in_dst = np.asarray(src_edges_in_dst, dtype=float)
    if len(src_edges_in_dst) != len(counts) + 1:
        raise ValueError("src_edges_in_dst must have length len(counts) + 1")
    if not np.all(np.diff(src_edges_in_dst) > 0):
        raise ValueError("src_edges_in_dst must be strictly increasing")
    cum = np.concatenate([[0.0], np.cumsum(counts)])
    dst_edges = np.arange(n_dst + 1) - 0.5
    F = np.interp(dst_edges, src_edges_in_dst, cum, left=0.0, right=cum[-1])
    return np.diff(F)


def rebin_conserve_var(counts, src_edges_in_dst, n_dst=N_CH):
    """Правка Claude (после первого прогона stage5): та же перебинировка, плюс честная
    дисперсия var_j = sum_i f_ij^2 n_i. С дисперсией = счёт хи-квадрат получал
    модуляцию с периодом 1 канал (дробный сдвиг усредняет соседние каналы) -> ложные минимумы."""
    e = np.asarray(src_edges_in_dst, dtype=float)
    n = np.asarray(counts, dtype=float)
    a, b = e[:-1], e[1:]
    w = b - a
    vals = np.zeros(n_dst); var = np.zeros(n_dst)
    j0 = np.floor(a + 0.5).astype(int)
    for o in range(0, int(np.ceil(w.max())) + 2):
        j = j0 + o
        f = np.clip(np.minimum(b, j + 0.5) - np.maximum(a, j - 0.5), 0, None) / w
        ok = (j >= 0) & (j < n_dst) & (f > 0)
        np.add.at(vals, j[ok], f[ok] * n[ok]); np.add.at(var, j[ok], f[ok] ** 2 * n[ok])
    return vals, var


def edges_from_mapping(mapping, n_src=N_CH) -> np.ndarray:
    x = np.arange(n_src + 1) - 0.5
    return np.asarray(mapping(x), dtype=float)


def smooth_gauss(y, sigma_ch) -> np.ndarray:
    return gaussian_filter1d(y, sigma_ch, mode="nearest")


def tag_print(tag, *items):
    print("[" + tag + "] " + " ".join(str(i) for i in items), flush=True)


# --- перенесено из прототипа subtract.py:19-65 (net, window_quality) ---
def net(s_on_b, bg_counts, k):
    """Compute residual and variance after subtracting live-time-scaled background."""
    residual = np.asarray(s_on_b) - k * bg_counts
    s_var = getattr(s_on_b, "var_s", None)
    variance = (np.asarray(s_var) if s_var is not None else np.asarray(s_on_b)) + k * k * bg_counts
    variance = np.maximum(variance, 1.0)
    return residual, variance


def window_quality(resid, var, lo, hi, center=None, detrend=False):
    """Quantify subtraction quality in a channel window [lo, hi]."""
    r = resid[lo:hi + 1]
    v = var[lo:hi + 1]
    x = np.arange(lo, hi + 1, dtype=float)
    n = len(r)

    if detrend:
        w = 1.0 / np.sqrt(v)
        coeffs = np.polyfit(x, r, 1, w=w)
        r = r - np.polyval(coeffs, x)

    sum_res = r.sum()
    sum_sigma = np.sqrt(v.sum())
    z_sum = sum_res / sum_sigma

    chi2 = np.sum(r * r / v)
    ndf = n - (2 if detrend else 0)
    chi2ndf = chi2 / ndf
    p_equiv_z = (chi2 - ndf) / np.sqrt(2 * ndf)

    if center is None:
        center = (lo + hi) / 2.0
    w_dip = np.sign(x - center)
    z_dip = np.sum(w_dip * r) / np.sqrt(np.sum(w_dip * w_dip * v))

    return {
        'lo': lo,
        'hi': hi,
        'n': n,
        'sum_res': sum_res,
        'sum_sigma': sum_sigma,
        'z_sum': z_sum,
        'chi2ndf': chi2ndf,
        'z_chi2': p_equiv_z,
        'z_dip': z_dip,
        'detrended': detrend
    }
