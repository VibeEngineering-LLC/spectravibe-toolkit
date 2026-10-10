"""Разложение фона LaBr3 на собственную и внешнюю части; инвентарь фонов (дубли, группы).

Метод (перенос отношений площадей внешних линий): каждый НЕэкранированный фон калибруется по своим внешним
линиям; отношения K-40/прокси и Tl-208/прокси умножаются на площадь прокси в эталонном (экранированном) фоне;
внешняя часть эталона = сумма аналитических пиков, собственная = эталон минус внешняя. Доля K-40 под сумм-пиком La-138.
Перенесено из прототипа labr3_selfcal/templates (extlines.py, decompose.py:28-182, tplcore.py:12-46,213-235,
inventory.py, dupmatrix.py), 2026-10-10.
"""
import json, pathlib, hashlib
import numpy as np
from scipy.special import erfc

from gamma.io.readers import read_spectrum
from .core import poly_eval, poly_inverse


# --- перенесено из tplcore.py:12-46 (load_any: 8192 -> 2048 суммированием по 4) ---
def load_any(path) -> dict:
    s = read_spectrum(str(path))
    raw = np.asarray(s.counts, dtype=float)
    n_raw = len(raw)
    md5 = hashlib.md5(raw.astype(np.int64).tobytes()).hexdigest()[:8]
    cal = [float(c) for c in s.energy_cal]
    if n_raw > 4096:
        factor = 4
        if n_raw < 8192:
            raw = np.pad(raw, (0, 8192 - n_raw), mode='constant')
        else:
            raw = raw[:8192]
        counts = raw.reshape(2048, 4).sum(axis=1)
        c = np.arange(2048.0)
        E = poly_eval(cal, 4.0 * c + 1.5)
        deg = min(max(len(cal) - 1, 1), 3)
        cal2048 = list(np.polynomial.polynomial.polyfit(c, E, deg))
    else:
        if n_raw < 2048:
            raw = np.pad(raw, (0, 2048 - n_raw), mode='constant')
        else:
            raw = raw[:2048]
        counts = raw
        cal2048 = cal
    return {
        "counts": counts,
        "live": float(s.live_time),
        "cal": tuple(cal2048),
        "start": str(s.start_datetime),
        "md5": md5,
        "n_raw": n_raw,
        "name": pathlib.Path(path).name,
        "path": str(path)
    }


# --- перенесено из tplcore.py:213-235 (аналитический пик: гаусс с левым экспоненциальным хвостом) ---
def analytic_peak_counts(y_lo, y_hi, mu, sigma, T=0.75, step=0.0):
    if sigma <= 0:
        return 0.0
    yy = np.arange(mu - 40 * sigma - 5, mu + 12 * sigma + 5, sigma / 20)
    dy = yy[1] - yy[0] if len(yy) > 1 else 1.0
    
    z = (yy - mu) / sigma
    density = np.zeros_like(yy)
    mask = z >= -T
    density[mask] = np.exp(-0.5 * z[mask] ** 2)
    density[~mask] = np.exp(T * z[~mask] + T * T / 2)
    
    integral = np.sum(density) * dy
    if integral > 0:
        density = density / integral
        
    step_density = step * 0.5 * erfc((yy - mu) / (sigma * np.sqrt(2)))
    total_density = density + step_density
    
    cum = np.cumsum(total_density) * dy
    # fix Claude: float() broke array input (bins); return array
    return np.interp(y_hi, yy, cum, left=0.0, right=cum[-1]) - np.interp(y_lo, yy, cum, left=0.0, right=cum[-1])



# --- перенесено из decompose.py:28-182 (fit_peak_diff, amp_of_structure) ---
def fit_peak_diff(d, v, mu0, sig0, mu_range, bg_deg=2, half_width_sig=3.5):
    # Fixed window based on initial guess and mu_range
    lo = int(np.floor(mu0 - mu_range - half_width_sig*sig0 - 3))
    hi = int(np.ceil(mu0 + mu_range + half_width_sig*sig0 + 3))
    lo = max(1, lo)
    hi = min(2046, hi)
    
    x = np.arange(lo, hi + 1)
    n = len(x)
    xc = (lo + hi) / 2.0
    
    # Variance flooring
    v_win = v[lo:hi+1]
    v_min = np.min(v_win[v_win > 0]) if np.any(v_win > 0) else 1e-30
    v_win = np.maximum(v_win, v_min)
    w = 1.0 / v_win
    
    # Grid
    mus = np.arange(mu0 - mu_range, mu0 + mu_range + 1e-9, 0.1)
    sigs = sig0 * np.arange(0.6, 1.6001, 0.05)
    
    best_chi2 = np.inf
    best_mu = mu0
    best_sig = sig0
    best_coef = None
    best_C = None
    
    # Precompute background columns
    bg_cols = []
    for p in range(bg_deg + 1):
        bg_cols.append(((x - xc) / 10.0) ** p)
    bg_matrix = np.column_stack(bg_cols)
    
    # Iterate grid
    chi2_grid = np.zeros((len(sigs), len(mus)))
    
    for i_sig, sig in enumerate(sigs):
        for i_mu, mu in enumerate(mus):
            # Peak column
            # analytic_peak_counts(x_lo, x_hi, mu, sigma, fwhm_ratio, offset)
            # Assuming signature: analytic_peak_counts(x_left, x_right, mu, sigma, fwhm_ratio, offset)
            # The prompt says: analytic_peak_counts(x - 0.5, x + 0.5, mu, sigma, 0.75, 0.0)
            # This implies x is the channel index.
            peak_col = analytic_peak_counts(x - 0.5, x + 0.5, mu, sig, 0.75, 0.0)
            
            M = np.column_stack([peak_col] + bg_cols)
            
            # Weighted LSQ
            # Minimize sum w * (d - M @ coef)^2
            # Normal equations: (M^T W M) coef = M^T W d
            W = np.diag(w) # For small n, explicit is fine, but let's use matrix ops
            # M is (n, k), W is (n, n)
            # A = M.T @ W @ M
            # b = M.T @ W @ d
            # To avoid creating large W, use:
            # A = (M * w[:, None]).T @ M
            # b = (M * w[:, None]).T @ d
            
            A = (M * w[:, None]).T @ M
            b = (M * w[:, None]).T @ d[lo:hi+1]  # fix Claude: window slice (generated code used the full array)
            
            try:
                coef = np.linalg.solve(A, b)
            except np.linalg.LinAlgError:
                chi2_grid[i_sig, i_mu] = np.inf
                continue
                
            resid = d[lo:hi+1] - M @ coef
            chi2 = np.sum(w * resid**2)
            chi2_grid[i_sig, i_mu] = chi2
            
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_mu = mu
                best_sig = sig
                best_coef = coef
                best_C = np.linalg.inv(A)
                
    # Count local minima in mu at best sigma
    i_sig_best = np.argmin(chi2_grid, axis=0) # This is wrong, we need the row corresponding to best_sig
    # Find index of best_sig in sigs
    i_sig_best_idx = np.argmin(np.abs(sigs - best_sig))
    row_chi2 = chi2_grid[i_sig_best_idx, :]
    
    n_min = 0
    for i in range(1, len(row_chi2) - 1):
        if row_chi2[i] < row_chi2[i-1] and row_chi2[i] < row_chi2[i+1]:
            n_min += 1
    # Check endpoints? Usually local minima are interior.
    
    # Final stats
    ndf = n - (bg_deg + 1) - 1 # n - (1 peak + bg_deg+1 bg)
    if ndf <= 0:
        ndf = 1
    s = max(1.0, best_chi2 / ndf)
    
    area = best_coef[0]
    area_err = np.sqrt(best_C[0, 0] * s)
    if area_err == 0:
        z = 0.0
    else:
        z = area / area_err
        
    return {
        "mu": best_mu,
        "sigma": best_sig,
        "fwhm": 2.3548 * best_sig,
        "area": area,
        "area_err": area_err,
        "z": z,
        "chi2ndf": best_chi2 / ndf,
        "n_min": n_min,
        "lo": lo,
        "hi": hi
    }

def amp_of_structure(d, v, ref, lo, hi, bg_deg=1):
    x = np.arange(lo, hi + 1)
    n = len(x)
    xc = (lo + hi) / 2.0
    
    v_win = v[lo:hi+1]
    v_min = np.min(v_win[v_win > 0]) if np.any(v_win > 0) else 1e-30
    v_win = np.maximum(v_win, v_min)
    w = 1.0 / v_win
    
    d_win = d[lo:hi+1]
    ref_win = ref[lo:hi+1]
    
    bg_cols = []
    for p in range(bg_deg + 1):
        bg_cols.append(((x - xc) / 10.0) ** p)
        
    M = np.column_stack([ref_win] + bg_cols)
    
    A = (M * w[:, None]).T @ M
    b = (M * w[:, None]).T @ d_win
    
    try:
        coef = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        return 0.0, np.inf
        
    resid = d_win - M @ coef
    chi2 = np.sum(w * resid**2)
    
    ndf = n - M.shape[1]
    if ndf <= 0:
        ndf = 1
    s = max(1.0, chi2 / ndf)
    
    a = coef[0]
    a_err = np.sqrt(np.linalg.inv(A)[0, 0] * s)
    
    return a, a_err

# --- перенесено из extlines.py:15-343 ---
LINES = {"Pb214_295": 295.22, "Pb214_352": 351.93, "Tl208_583": 583.19, "Bi214_609": 609.32, "Ac228_911": 911.2,
         "Bi214_1120": 1120.29, "K40": 1460.82, "Bi214_1764": 1764.49, "Tl208_2614": 2614.51}
CAL_LINES = ["Pb214_352", "Bi214_609", "Bi214_1120", "K40", "Bi214_1764", "Tl208_2614"]
PROXIES = ["Pb214_352", "Bi214_609", "Tl208_583", "Ac228_911", "Bi214_1120", "Bi214_1764"]
TARGETS = ["K40", "Tl208_2614"]

def fwhm_ch(E, cal):
    # FWHM in channels for energy E keV
    fwhm_kev = 0.033 * np.sqrt(662.0 * E)
    ch = poly_inverse(cal, E)
    # Local dispersion dE/dch
    e_plus = poly_eval(cal, ch + 0.5)
    e_minus = poly_eval(cal, ch - 0.5)
    dE_dch = e_plus - e_minus
    if dE_dch == 0:
        return 1.0
    return fwhm_kev / dE_dch

def selfcal(sp, tag):
    # Two-pass self-calibration
    cal = sp["cal"]
    counts = sp["counts"]
    var = np.maximum(counts, 1.0)
    
    for p in [1, 2]:
        kept = []
        for name in CAL_LINES:
            E = LINES[name]
            mu0 = poly_inverse(cal, E)
            sig0 = fwhm_ch(E, cal) / 2.3548
            mu_range = 12.0 if p == 1 else 4.0
            r = fit_peak_diff(counts, var, mu0, sig0, mu_range)
            if r["z"] > 5:
                kept.append((name, E, r))
                print(f"[SELF] {tag} pass={p} {name} E={E} mu={r['mu']:.2f} fwhm={r['fwhm']:.2f} area={r['area']:.0f}±{r['area_err']:.0f} z={r['z']:.1f}")
            else:
                print(f"[SELF] {tag} pass={p} {name} E={E} mu={r['mu']:.2f} fwhm={r['fwhm']:.2f} area={r['area']:.0f}±{r['area_err']:.0f} z={r['z']:.1f} (rejected)")
        
        if len(kept) >= 3:
            # Weighted least squares fit for quadratic calibration
            mus = np.array([k[2]["mu"] for k in kept])
            Es = np.array([k[1] for k in kept])
            ws = np.array([k[2]["z"]**2 for k in kept])
            
            # Design matrix for c0 + c1*ch + c2*ch^2
            A = np.vstack([np.ones_like(mus), mus, mus**2]).T
            W = np.diag(ws)
            # Solve (A^T W A) c = A^T W E
            ATA = A.T @ W @ A
            ATE = A.T @ W @ Es
            try:
                c = np.linalg.solve(ATA, ATE)
                cal = [c[0], c[1], c[2]]
            except np.linalg.LinAlgError:
                pass # keep previous cal
        elif len(kept) >= 2:
            # Linear fit
            mus = np.array([k[2]["mu"] for k in kept])
            Es = np.array([k[1] for k in kept])
            ws = np.array([k[2]["z"]**2 for k in kept])
            
            A = np.vstack([np.ones_like(mus), mus]).T
            W = np.diag(ws)
            ATA = A.T @ W @ A
            ATE = A.T @ W @ Es
            try:
                c = np.linalg.solve(ATA, ATE)
                cal = [c[0], c[1], 0.0]
            except np.linalg.LinAlgError:
                pass
        else:
            # Keep file calibration
            pass
            
        # Print calibration and residuals
        if len(kept) > 0:
            resids = []
            for name, E, r in kept:
                mu = r["mu"]
                E_calc = poly_eval(cal, mu)
                resids.append(E - E_calc)
            resid_str = ",".join([f"{r:.2f}" for r in resids])
            print(f"[SELFCAL] {tag} c={cal[0]:.4f},{cal[1]:.6f},{cal[2]:.4e} resid_keV={resid_str}")
        else:
            print(f"[SELFCAL] {tag} c={cal[0]:.4f},{cal[1]:.6f},{cal[2]:.4e} resid_keV=none")
            
    return cal

def line_areas(counts, var, cal, live, tag):
    # Fit all lines and return areas in cps
    results = {}
    for name in LINES:
        E = LINES[name]
        mu0 = poly_inverse(cal, E)
        sig0 = fwhm_ch(E, cal) / 2.3548
        mu_range = 4.0
        if name == "Tl208_2614":
            mu_range = 8.0
            
        r = fit_peak_diff(counts, var, mu0, sig0, mu_range)
        
        # Convert area to cps
        rate = r["area"] / live
        err = r["area_err"] / live
        
        print(f"[AREA] {tag} {name} E={E} mu={r['mu']:.2f} fwhm={r['fwhm']:.2f} rate={rate:.3e}±{err:.1e} cps z={r['z']:+.1f} chi2/ndf={r['chi2ndf']:.2f}")
        
        results[name] = {
            "rate": rate,
            "err": err,
            "z": r["z"],
            "mu": r["mu"],
            "sigma": r["sigma"]
        }
    return results

def ref_windows(counts):
    """Правка Claude (#LABR-8): окна сумм-пика (810-880, ядро 832-862) и α1 (1290-1400), заданные в каналах эталона
    2048, -- в каналы этого прибора (load_any приводит к 2048 суммированием, но усиление у приборов разное)."""
    from . import autopreset as ap
    base = {"sum": (810, 880), "core": (832, 862), "a1": (1290, 1400)}
    a = ap.autopreset(counts, counts, None, verbose=False)
    if a["identity"]:
        return base
    print("[WIN] окна эталона пересчитаны по собственным структурам: масштаб %.4f" % a["r"])
    return {k: tuple(ap.map_win(a["map"], lo, hi, len(counts))) for k, (lo, hi) in base.items()}


def check_ecal(cal, counts, what):
    """Правка Claude (#LABR-3): понятная ошибка вместо трассировки/мусора при вырожденной шкале (например --ecal 0 1 0).
    Шкала должна быть конечной, возрастать на всех каналах и ставить собственный сумм-пик La-138 около 1471 кэВ."""
    from . import autopreset as ap
    c = [float(v) for v in cal]
    if len(c) < 2 or not all(np.isfinite(c)):
        raise ValueError("%s: нужны конечные коэффициенты c0 c1 [c2 ...] по возрастанию степени, получено %s" % (what, c))
    e = poly_eval(c, np.arange(len(counts) + 1, dtype=float))
    if np.any(np.diff(e) <= 0):
        raise ValueError("%s: шкала %s не возрастает на каналах 0..%d" % (what, c, len(counts)))
    ref = ap.load_ref()
    s = ap.find_sum(counts, ap.find_a37(counts, ref)["x"], ref)
    E = float(poly_eval(c, s["x"]))
    if not 1300.0 <= E <= 1650.0:
        raise ValueError("%s: шкала %s ставит собственный сумм-пик La-138 (канал %.1f) на %.0f кэВ, ожидается ~1471 кэВ; "
                         "проверьте коэффициенты (порядок по возрастанию степени: c0 c1 c2)" % (what, c, s["x"], E))


def decompose(ref_path, unshielded_paths, e_ref=None, out_json=None):
    """ref_path -- экранированный (эталонный) фон; unshielded_paths -- неэкранированные фоны (внешние линии видны);
    e_ref -- приближённая шкала эталона кэВ(канал) по возрастанию (по умолчанию -- калибровка файла)."""
    ref = load_any(ref_path)
    E_REF = list(e_ref) if e_ref is not None else list(ref["cal"])
    check_ecal(E_REF, ref["counts"], "--ecal" if e_ref is not None else "калибровка файла эталона")   # правка Claude (#LABR-3)
    ar = {}
    for u in unshielded_paths:
        path = pathlib.Path(u)
        tag = path.stem
        sp = load_any(path)
        cal = selfcal(sp, tag)
        ar[tag] = line_areas(sp["counts"], np.maximum(sp["counts"], 1.0), cal, sp["live"], tag)
        
    # Reference areas
    arA = line_areas(ref["counts"], np.maximum(ref["counts"], 1.0), E_REF, ref["live"], "A")
    print("[NOTE] K40 and Tl208_2614 in A overlap the intrinsic sum peak / alpha band 1: their direct fits in A are NOT used.")
    
    # Ratio transfer
    transfers = {t: [] for t in TARGETS}
    
    for u_tag, u_areas in ar.items():
        for p in PROXIES:
            if p not in u_areas or p not in arA:
                continue
            if u_areas[p]["z"] > 3 and arA[p]["z"] > 2:
                for t in TARGETS:
                    if t not in u_areas:
                        continue
                    if u_areas[t]["z"] > 3:
                        R = u_areas[t]["rate"] / u_areas[p]["rate"]
                        est = arA[p]["rate"] * R
                        
                        # Relative errors
                        rel_err_u_t = u_areas[t]["err"] / u_areas[t]["rate"] if u_areas[t]["rate"] > 0 else 0
                        rel_err_u_p = u_areas[p]["err"] / u_areas[p]["rate"] if u_areas[p]["rate"] > 0 else 0
                        rel_err_A_p = arA[p]["err"] / arA[p]["rate"] if arA[p]["rate"] > 0 else 0
                        
                        # Quadrature of relative errors
                        rel_err_est = np.sqrt(rel_err_u_t**2 + rel_err_u_p**2 + rel_err_A_p**2)
                        err_est = est * rel_err_est
                        
                        print(f"[XFER] {u_tag} proxy={p} target={t} R={R:.3f} est_A={est:.3e}±{err_est:.1e} cps")
                        
                        transfers[t].append({
                            "est": est,
                            "err": err_est,
                            "u": u_tag,
                            "p": p
                        })
                        
    # Summaries
    summaries = {}
    for t in TARGETS:
        if transfers[t]:
            ests = np.array([x["est"] for x in transfers[t]])
            errs = np.array([x["err"] for x in transfers[t]])
            
            n = len(ests)
            median = np.median(ests)
            p16 = np.percentile(ests, 16)
            p84 = np.percentile(ests, 84)
            
            # Weighted mean
            weights = 1.0 / (errs**2)
            wmean = np.sum(weights * ests) / np.sum(weights)
            werr = 1.0 / np.sqrt(np.sum(weights))
            
            # Chi2/ndf
            chi2 = np.sum(weights * (ests - wmean)**2)
            ndf = n - 1
            chi2ndf = chi2 / ndf if ndf > 0 else 0.0
            
            print(f"[XSUM] {t} n={n} median={median:.3e} p16={p16:.3e} p84={p84:.3e} wmean={wmean:.3e}±{werr:.1e} chi2ndf={chi2ndf:.2f}")
            
            summaries[t] = {
                "n": n,
                "median": float(median),
                "p16": float(p16),
                "p84": float(p84),
                "wmean": float(wmean),
                "werr": float(werr),
                "chi2ndf": float(chi2ndf)
            }
        else:
            print(f"[XSUM] {t} n=0")
            summaries[t] = {"n": 0}

    # Fractions
    W = ref_windows(ref["counts"])   # правка Claude (#LABR-8): окна в каналах этого прибора
    S_sum = np.sum(ref["counts"][W["sum"][0]:W["sum"][1] + 1]) / ref["live"]
    S_core = np.sum(ref["counts"][W["core"][0]:W["core"][1] + 1]) / ref["live"]
    
    # K40 fractions
    if "K40" in summaries and summaries["K40"]["n"] > 0:
        f_median = summaries["K40"]["median"] / S_sum
        f_p16 = summaries["K40"]["p16"] / S_sum
        f_p84 = summaries["K40"]["p84"] / S_sum
        f_wmean = summaries["K40"]["wmean"] / S_sum
        f_werr = summaries["K40"]["werr"] / S_sum
        print(f"[KFRAC] window=810-880 f_median={f_median:.3f} f_p16={f_p16:.3f} f_p84={f_p84:.3f} f_wmean={f_wmean:.3f}±{f_werr:.3f}")
        
        f_median_core = summaries["K40"]["median"] / S_core
        f_p16_core = summaries["K40"]["p16"] / S_core
        f_p84_core = summaries["K40"]["p84"] / S_core
        f_wmean_core = summaries["K40"]["wmean"] / S_core
        f_werr_core = summaries["K40"]["werr"] / S_core
        print(f"[KFRAC] window=832-862 f_median={f_median_core:.3f} f_p16={f_p16_core:.3f} f_p84={f_p84_core:.3f} f_wmean={f_wmean_core:.3f}±{f_werr_core:.3f}")
    else:
        print("[KFRAC] window=810-880 no data")
        print("[KFRAC] window=832-862 no data")
        
    # Tl208_2614 fractions
    # Alpha1 window 1290-1400
    S_alpha1 = np.sum(ref["counts"][W["a1"][0]:W["a1"][1] + 1]) / ref["live"]
    
    if "Tl208_2614" in summaries and summaries["Tl208_2614"]["n"] > 0:
        f_median = summaries["Tl208_2614"]["median"] / S_alpha1
        f_p16 = summaries["Tl208_2614"]["p16"] / S_alpha1
        f_p84 = summaries["Tl208_2614"]["p84"] / S_alpha1
        f_wmean = summaries["Tl208_2614"]["wmean"] / S_alpha1
        f_werr = summaries["Tl208_2614"]["werr"] / S_alpha1
        print(f"[TLFRAC] window=1290-1400 f_median={f_median:.3f} f_p16={f_p16:.3f} f_p84={f_p84:.3f} f_wmean={f_wmean:.3f}±{f_werr:.3f}")
    else:
        print("[TLFRAC] window=1290-1400 no data")

    # External peak component of A
    x = np.arange(2048.0)
    E_A = np.zeros(2048)
    E_A_lo = np.zeros(2048)
    E_A_hi = np.zeros(2048)
    
    for name in LINES:
        E = LINES[name]
        
        # Determine amplitude
        if name not in TARGETS:
            if name in arA and arA[name]["z"] > 2:
                amp = arA[name]["rate"]
                mu = arA[name]["mu"]
            else:
                continue
            amp_lo = amp
            amp_hi = amp
        else:
            # Use median transfer estimate
            if name in summaries and summaries[name]["n"] > 0:
                amp = summaries[name]["median"]
                amp_lo = summaries[name]["p16"]
                amp_hi = summaries[name]["p84"]
            else:
                continue
            mu = poly_inverse(E_REF, E)
            
        sigma = fwhm_ch(E, E_REF) / 2.3548
        
        # Add peaks
        peak = analytic_peak_counts(x - 0.5, x + 0.5, mu, sigma, 0.75, 0.0)
        E_A += amp * peak
        E_A_lo += amp_lo * peak
        E_A_hi += amp_hi * peak
        
        print(f"[EA] {name} amp={amp:.3e} cps mu={mu:.2f} sigma={sigma:.2f}")

    # Intrinsic component
    rates_A = ref["counts"] / ref["live"]
    var_A = np.maximum(ref["counts"], 1.0) / (ref["live"] ** 2)
    I = rates_A - E_A
    
    # Check negative channels
    neg_mask = I < -3 * np.sqrt(var_A)
    neg_count = np.sum(neg_mask[20:2040])
    print(f"[ICHK] negative_channels={neg_count}")
    
    # Prepare JSON data
    json_data = {
        "areas": {
            "ref": {k: {kk: float(vv) for kk, vv in v.items()} for k, v in arA.items()},
            "high": {u: {k: {kk: float(vv) for kk, vv in v.items()} for k, v in u_areas.items()} for u, u_areas in ar.items()}
        },
        "transfers": {t: transfers[t] for t in TARGETS},
        "summaries": summaries,
        "intrinsic_rate": [float(v) for v in I],
        "external_rate": [float(v) for v in E_A],
    }
    if out_json:
        with open(out_json, "w", encoding="utf-8") as f:
            json.dump(json_data, f, ensure_ascii=False, indent=1)
    return json_data
