"""Подгонка пика (гаусс + линейный фон + ступень) и края (erfc) в каналах. Перенесено из прототипа labr3_selfcal/peaks.py.

Отличие от gamma.peaks.area.gaussian_fit_area: ступень под пиком и оценка ошибки центроида с масштабом chi2/ndf;
find_peaks использует поиск пиков toolkit (gamma.peaks.search.mariscotti_search)."""
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erfc
from scipy.ndimage import gaussian_filter1d

from gamma.peaks.search import mariscotti_search


def gauss_lin_step(x, H, c, s, a0, a1, h):
    return H * np.exp(-0.5 * ((x - c) / s) ** 2) + a0 + a1 * (x - c) + h * 0.5 * erfc((x - c) / (np.sqrt(2) * s))


def fit_peak(counts, x0, fwhm0, half_window=None, step=True):
    if half_window is None:
        half_window = max(4, int(round(1.6 * fwhm0)))
    lo = max(0, int(round(x0)) - half_window)
    hi = min(len(counts) - 1, int(round(x0)) + half_window)
    if hi - lo + 1 < 7:
        return None
    y = counts[lo:hi + 1]
    x = np.arange(lo, hi + 1, dtype=float)
    sigma = np.sqrt(np.maximum(y, 1.0))
    first3 = y[:3]
    last3 = y[-3:]
    a0_init = (np.mean(first3) + np.mean(last3)) / 2.0
    a1_init = (np.mean(last3) - np.mean(first3)) / (x[-1] - x[0]) if x[-1] != x[0] else 0.0
    H_init = max(np.max(y) - a0_init, 1.0)
    c_init = x0
    s_init = fwhm0 / 2.3548
    h_init = 0.0
    p0 = [H_init, c_init, s_init, a0_init, a1_init, h_init]
    bounds_lo = [0, x0 - 0.5 * fwhm0 - 2, 0.3, -np.inf, -np.inf, -np.inf]
    bounds_hi = [np.inf, x0 + 0.5 * fwhm0 + 2, 3 * fwhm0, np.inf, np.inf, np.inf]
    try:
        if step:
            popt, pcov = curve_fit(gauss_lin_step, x, y, p0=p0, sigma=sigma, absolute_sigma=True,
                                   bounds=(bounds_lo, bounds_hi), maxfev=20000)
            nparams = 6
        else:
            def model5(x, H, c, s, a0, a1):
                return H * np.exp(-0.5 * ((x - c) / s) ** 2) + a0 + a1 * (x - c)
            popt, pcov = curve_fit(model5, x, y, p0=p0[:5], sigma=sigma, absolute_sigma=True,
                                   bounds=(bounds_lo[:5], bounds_hi[:5]), maxfev=20000)
            nparams = 5
            popt = np.append(popt, 0.0)
            pcov = np.pad(pcov, ((0, 1), (0, 1)))
    except Exception:
        return None
    H, c, s, a0, a1, h = popt
    y_model = gauss_lin_step(x, H, c, s, a0, a1, h)
    chi2 = np.sum(((y - y_model) / sigma) ** 2)
    ndf = len(y) - nparams
    chi2ndf = chi2 / max(ndf, 1)
    scale = np.sqrt(max(1.0, chi2ndf))
    perr = np.sqrt(np.diag(pcov)) * scale
    H_err, c_err, s_err, a0_err, a1_err, h_err = perr
    area = H * s * np.sqrt(2 * np.pi)
    # Error propagation for area = H*s*sqrt(2*pi)
    # d(area)/dH = s*sqrt(2*pi), d(area)/ds = H*sqrt(2*pi)
    # var_area = (s*sqrt(2*pi))^2 * var_H + (H*sqrt(2*pi))^2 * var_s + 2*(s*sqrt(2*pi))*(H*sqrt(2*pi))*cov_Hs
    sqrt_2pi = np.sqrt(2 * np.pi)
    var_H = pcov[0, 0]
    var_s = pcov[2, 2]
    cov_Hs = pcov[0, 2]
    var_area = (s * sqrt_2pi) ** 2 * var_H + (H * sqrt_2pi) ** 2 * var_s + 2 * (s * sqrt_2pi) * (H * sqrt_2pi) * cov_Hs
    area_err = np.sqrt(max(var_area, 0)) * scale
    fwhm = 2.3548 * s
    fwhm_err = 2.3548 * s_err
    signif = area / area_err if area_err > 0 else 0.0
    return {
        'centroid': c,
        'centroid_err': c_err,
        'fwhm': fwhm,
        'fwhm_err': fwhm_err,
        'area': area,
        'area_err': area_err,
        'height': H,
        'signif': signif,
        'chi2ndf': chi2ndf,
        'lo': lo,
        'hi': hi,
        'step': h
    }


def find_peaks(counts, fwhm_fn, threshold=4.0, ch_min=15, ch_max=2040, min_signif=3.0):
    found = mariscotti_search(counts, fwhm_fn, sigma_threshold=threshold)
    candidates = []
    for p in found:
        if ch_min <= p.channel <= ch_max:
            r = fit_peak(counts, p.channel, p.fwhm_channels)
            if r is not None and r['signif'] >= min_signif and r['area'] > 0 and r['centroid_err'] < r['fwhm']:
                r['search_channel'] = p.channel
                r['search_signif'] = float(p.significance)
                candidates.append(r)
    candidates.sort(key=lambda x: x['centroid'])
    kept = []
    for r in candidates:
        duplicate = False
        for k in kept:
            if abs(r['centroid'] - k['centroid']) < 0.5 * max(r['fwhm'], k['fwhm']):
                if r['signif'] > k['signif']:
                    kept.remove(k)
                    break
                else:
                    duplicate = True
                    break
        if not duplicate:
            kept.append(r)
    kept.sort(key=lambda x: x['centroid'])
    return kept


def edge_model(x, xe, s, Hl, Hr, al, ar, rising):
    t = 0.5 * erfc(-(x - xe) / (np.sqrt(2) * s)) if rising else 0.5 * erfc((x - xe) / (np.sqrt(2) * s))
    L = Hl + al * (x - xe)
    R = Hr + ar * (x - xe)
    if rising:
        return L * (1 - t) + R * t
    else:
        return L * t + R * (1 - t)


def fit_edge(counts, xe0, s0, half_window, rising, fix_sigma=False):
    lo = max(0, int(xe0) - half_window)
    hi = min(len(counts) - 1, int(xe0) + half_window)
    if hi - lo + 1 < 7:
        return None
    y = counts[lo:hi + 1]
    x = np.arange(lo, hi + 1, dtype=float)
    sigma = np.sqrt(np.maximum(y, 1.0))
    Hl_init = np.mean(y[:5])
    Hr_init = np.mean(y[-5:])
    al_init = 0.0
    ar_init = 0.0
    xe_init = xe0
    s_init = s0
    try:
        if fix_sigma:
            def model_fixed(x, xe, Hl, Hr, al, ar):
                return edge_model(x, xe, s0, Hl, Hr, al, ar, rising)
            p0 = [xe_init, Hl_init, Hr_init, al_init, ar_init]
            bounds_lo = [xe0 - half_window / 2, -np.inf, -np.inf, -np.inf, -np.inf]
            bounds_hi = [xe0 + half_window / 2, np.inf, np.inf, np.inf, np.inf]
            popt, pcov = curve_fit(model_fixed, x, y, p0=p0, sigma=sigma, absolute_sigma=True,
                                   bounds=(bounds_lo, bounds_hi), maxfev=20000)
            nparams = 5
            xe, Hl, Hr, al, ar = popt
            s = s0
            s_err = 0.0
        else:
            def model_full(x, xe, s, Hl, Hr, al, ar):
                return edge_model(x, xe, s, Hl, Hr, al, ar, rising)
            p0 = [xe_init, s_init, Hl_init, Hr_init, al_init, ar_init]
            bounds_lo = [xe0 - half_window / 2, 0.3, -np.inf, -np.inf, -np.inf, -np.inf]
            bounds_hi = [xe0 + half_window / 2, 4 * s0 + 1, np.inf, np.inf, np.inf, np.inf]
            popt, pcov = curve_fit(model_full, x, y, p0=p0, sigma=sigma, absolute_sigma=True,
                                   bounds=(bounds_lo, bounds_hi), maxfev=20000)
            nparams = 6
            xe, s, Hl, Hr, al, ar = popt
            s_err = np.sqrt(np.diag(pcov)[1])
    except Exception:
        return None
    y_model = edge_model(x, xe, s, Hl, Hr, al, ar, rising)
    chi2 = np.sum(((y - y_model) / sigma) ** 2)
    ndf = len(y) - nparams
    chi2ndf = chi2 / max(ndf, 1)
    scale = np.sqrt(max(1.0, chi2ndf))
    perr = np.sqrt(np.diag(pcov)) * scale
    xe_err = perr[0]
    return {
        'edge': xe,
        'edge_err': xe_err,
        'sigma': s,
        'sigma_err': s_err,
        'Hl': Hl,
        'Hr': Hr,
        'chi2ndf': chi2ndf,
        'lo': lo,
        'hi': hi,
        'rising': rising
    }


def find_edges(counts, sigma_fn, min_z=5.0, ch_min=15, ch_max=2040, step=4):
    results = []
    for c in range(ch_min, ch_max, step):
        s = max(1.0, sigma_fn(c))
        win = int(4 * s)
        lo = max(0, c - win)
        hi = min(len(counts) - 1, c + win)
        if hi - lo + 1 < 3:
            continue
        x = np.arange(lo, hi + 1, dtype=float)
        y = counts[lo:hi + 1]
        # Ручная правка: знак в спецификации был неверен — d/dc сглаженного
        # спектра = sum y_i*(i-c)/s^2*g; z>0 — счёт растёт с каналом.
        w = (x - c) / (s ** 2) * np.exp(-0.5 * ((x - c) / s) ** 2)
        d = np.sum(w * y)
        var = np.sum(w ** 2 * np.maximum(y, 1.0))
        if var <= 0:
            continue
        z = d / np.sqrt(var)
        results.append((c, z, s))
    edges = []
    for i in range(1, len(results) - 1):
        c, z, s = results[i]
        if abs(z) >= min_z and abs(z) > abs(results[i - 1][1]) and abs(z) > abs(results[i + 1][1]):
            edges.append({
                'channel': c,
                'z': z,
                'kind': 'rising' if z > 0 else 'falling',
                'sigma': s
            })
    edges.sort(key=lambda x: x['channel'])
    return edges


def fwhm_fn_from_points(ch, fwhm):
    ch = np.asarray(ch, dtype=float)
    fwhm = np.asarray(fwhm, dtype=float)
    p1, p0 = np.polyfit(ch, fwhm ** 2, 1)
    def fn(x):
        return float(np.sqrt(max(p0 + p1 * x, 1.0)))
    fn.coef = (p0, p1)
    return fn
