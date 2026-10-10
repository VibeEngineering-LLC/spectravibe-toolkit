"""Энергетическая шкала по реперам: взвешенный полином с эффективной дисперсией, невязки, LOO.
Перенесено из прототипа labr3_selfcal/calib.py и make_net_xml.py:47-60 (scan -- разброс центроида по 5 окнам подгонки).
Отличие от gamma.calibration.energy_fit.polynomial_energy_fit: веса по ошибкам канала и энергии, LOO, ковариация."""
import numpy as np
from .core import poly_eval, poly_deriv, poly_inverse
from .peakfit import fit_peak, fit_edge


def fit_poly(ch, ch_err, E, E_err, deg):
    """Fit a weighted polynomial E(ch) of given degree using effective variance iteration."""
    ch = np.asarray(ch, dtype=float)
    ch_err = np.asarray(ch_err, dtype=float)
    E = np.asarray(E, dtype=float)
    E_err = np.asarray(E_err, dtype=float)
    n = len(ch)
    if n != len(ch_err) or n != len(E) or n != len(E_err):
        raise ValueError("All inputs must have equal length")
    if deg < 1:
        raise ValueError("deg must be >= 1")
    if n <= deg:
        raise ValueError("n must be > deg")

    # Initial slope estimate
    slope0 = (np.max(E) - np.min(E)) / (np.max(ch) - np.min(ch))

    def _wls(var):
        A = np.column_stack([ch**i for i in range(deg + 1)])
        W = 1.0 / var
        ATA = A.T @ (W[:, None] * A)
        ATE = A.T @ (W * E)
        coeffs = np.linalg.solve(ATA, ATE)
        cov = np.linalg.inv(ATA)
        return coeffs, cov

    # Pass 1
    var = E_err**2 + (ch_err * slope0)**2
    coeffs, cov = _wls(var)

    # Pass 2
    slope = poly_deriv(coeffs, ch)
    var = E_err**2 + (ch_err * slope)**2
    coeffs, cov = _wls(var)

    residuals = E - poly_eval(coeffs, ch)
    pulls = residuals / np.sqrt(var)
    chi2 = np.sum(pulls**2)
    ndf = n - (deg + 1)
    chi2ndf = chi2 / ndf if ndf > 0 else np.nan
    rms = np.sqrt(np.mean(residuals**2))
    max_abs = np.max(np.abs(residuals))

    return {
        'deg': deg,
        'coeffs': tuple(float(c) for c in coeffs),
        'cov': cov,
        'residuals': residuals,
        'pulls': pulls,
        'chi2': chi2,
        'ndf': ndf,
        'chi2ndf': chi2ndf,
        'rms': rms,
        'max_abs': max_abs,
    }


def energy_err(fit, ch):
    """Standard uncertainty of E(ch) from coefficient covariance."""
    ch = np.asarray(ch, dtype=float)
    deg = fit['deg']
    cov = fit['cov']
    if ch.ndim == 0:
        g = np.array([ch**i for i in range(deg + 1)])
        return float(np.sqrt(g @ cov @ g))
    else:
        out = np.empty_like(ch, dtype=float)
        for idx, c in enumerate(ch):
            g = np.array([c**i for i in range(deg + 1)])
            out[idx] = np.sqrt(g @ cov @ g)
        return out


def compare_degrees(ch, ch_err, E, E_err, degs=(1, 2, 3)):
    """Fit multiple degrees and return list of fit dicts for feasible degrees."""
    n = len(ch)
    fits = []
    for d in degs:
        if n > d:
            fits.append(fit_poly(ch, ch_err, E, E_err, d))
    return fits


def choose_degree(fits):
    """Choose the best degree based on chi2ndf criterion."""
    if not fits:
        raise ValueError("No fits provided")

    valid = [f for f in fits if f['ndf'] > 0]
    if not valid:
        best = max(fits, key=lambda f: f['deg'])
        return {**best, 'reason': 'No fits with ndf>0; chose highest degree'}

    # Check criterion: chi2ndf <= 1.0 + 2.0*sqrt(2.0/ndf)
    for f in sorted(valid, key=lambda x: x['deg']):
        threshold = 1.0 + 2.0 * np.sqrt(2.0 / f['ndf'])
        if f['chi2ndf'] <= threshold:
            return {**f, 'reason': f"Lowest degree with chi2ndf <= {threshold:.3f}"}

    # Fallback: minimal chi2ndf among valid
    best = min(valid, key=lambda f: f['chi2ndf'])
    return {**best, 'reason': 'No degree met criterion; chose minimal chi2ndf'}


def leave_one_out(ch, ch_err, E, E_err, deg):
    """Cross-validation residuals: fit without each point, predict residual."""
    ch = np.asarray(ch, dtype=float)
    ch_err = np.asarray(ch_err, dtype=float)
    E = np.asarray(E, dtype=float)
    E_err = np.asarray(E_err, dtype=float)
    n = len(ch)
    out = np.full(n, np.nan)
    for i in range(n):
        if n - 1 <= deg:
            continue
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        try:
            fit = fit_poly(ch[mask], ch_err[mask], E[mask], E_err[mask], deg)
            e_pred = poly_eval(fit['coeffs'], ch[i])
            out[i] = E[i] - e_pred
        except Exception:
            out[i] = np.nan
    return out


def channel_of_energy(fit, E):
    """Invert the polynomial to find channel for given energy."""
    return poly_inverse(fit['coeffs'], E)


# --- перенесено из прототипа make_net_xml.py:47-60 (scan: разброс центроида по 5 окнам) ---
def scan(counts, a, kind):
    cs, es = [], []
    for f in (0.7, 0.85, 1.0, 1.15, 1.3):
        hw = max(5, int(round(a["hw"] * f)))
        r = fit_peak(counts, a["x0"], a["w0"], half_window=hw) if kind == "peak" else fit_edge(counts, a["x0"], a["w0"], hw, True)
        if r is None:
            continue
        cs.append(r["centroid"] if kind == "peak" else r["edge"])
        es.append(r["centroid_err"] if kind == "peak" else r["edge_err"])
    cs, es = np.array(cs), np.array(es)
    ok = (np.abs(cs - np.median(cs)) < 3.0) & (es < 3.0)      # выбросы подгонки (расходящееся окно) отбрасываем
    print("      окон принято %d из %d" % (ok.sum(), len(cs)))
    cs, es = cs[ok], es[ok]
    return cs.mean(), float(np.hypot(cs.std(ddof=1), np.mean(es))), cs
