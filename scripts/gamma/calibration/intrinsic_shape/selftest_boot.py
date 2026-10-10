"""Проверки сдвига/масштаба/амплитуды для selftest: глобальный минимум на сетке + параметрический бутстрап.
Перенесено из прототипа labr3_selfcal/robust_shift.py (сгенерирован gen_code.py, qwen3.8:27b), 2026-10-10."""
import numpy as np
from .core import rebin_conserve_var, poly_eval, N_CH, window_quality
from . import register as sr


def resample(S, B, bg_var, rng):
    S_r = rng.poisson(np.maximum(S, 0)).astype(float)
    if bg_var is None:
        B_r = rng.poisson(np.maximum(B, 0)).astype(float)
    else:
        B_r = B + rng.normal(0.0, np.sqrt(np.maximum(bg_var, 0.0)))
    return S_r, B_r

def refine_min(ds, X):
    i = int(np.argmin(X))
    if i == 0 or i == len(ds) - 1:
        return (float(ds[i]), True)
    den = X[i+1] - 2*X[i] + X[i-1]
    if den <= 0:
        return (float(ds[i]), False)
    off = 0.5 * (X[i-1] - X[i+1]) / den
    off = np.clip(off, -0.5, 0.5)
    step = ds[1] - ds[0]
    return (float(ds[i] + off * step), False)

def window_shifts(c, S, B, k, bg_var, wins, ds):
    c = np.asarray(c, float)
    X = {n: np.zeros(len(ds)) for n, lo, hi in wins}
    for a, d in enumerate(ds):
        cd = c.copy()
        cd[0] += d
        r, v = sr.apply(cd, S, B, k, bg_var)
        q = r * r / v
        for n, lo, hi in wins:
            X[n][a] = q[lo:hi+1].sum()
    return {n: refine_min(ds, X[n]) for n, lo, hi in wins}

def joint_grid(c, S, B, k, bg_var, mask, pivot, g):
    ds = np.arange(g["d_lo"], g["d_hi"] + 1e-9, g["d_step"])
    ss = np.arange(g["s_lo"], g["s_hi"] + 1e-12, g["s_step"])
    c = np.asarray(c, float)
    xe = np.arange(len(S) + 1) - 0.5
    best_x = np.inf
    best_d = 0.0
    best_s = 0.0
    best_a = 0.0
    best_i_d = 0
    best_i_s = 0
    for i_s, s in enumerate(ss):
        for i_d, d in enumerate(ds):
            cc = (1 + s) * c
            cc[0] += d - s * pivot
            e = poly_eval(cc, xe)
            if np.any(np.diff(e) <= 0):
                continue
            sv, vs = rebin_conserve_var(S, e, len(B))
            bb, vb = sr.blur_bg(lambda x, cc=cc: poly_eval(cc, x), B, bg_var)
            w = 1.0 / np.maximum(vs + k * k * vb, 1.0)
            m = mask
            kb = k * bb
            denom = np.sum(kb[m] ** 2 * w[m])
            if denom == 0:
                continue
            a = np.sum(sv[m] * kb[m] * w[m]) / denom
            X = np.sum((sv[m] - a * kb[m]) ** 2 * w[m])
            if X < best_x:
                best_x = X
                best_d = d
                best_s = s
                best_a = a
                best_i_d = i_d
                best_i_s = i_s
    at_edge = best_i_d in (0, len(ds) - 1) or best_i_s in (0, len(ss) - 1)
    return (float(best_d), float(best_s), float(best_a), bool(at_edge))

def shifts_cfg(c, S, B, k, bg_var, wins, ds, cfg):
    """Правка Claude (v3): α-окна из cfg['wide_scan_windows'] сканируются по cfg['shift_scan_wide'] (+-4 кан)."""
    wide = set(cfg.get("wide_scan_windows", [])) if "shift_scan_wide" in cfg else set()
    narrow = [w for w in wins if w[0] not in wide]
    out = window_shifts(c, S, B, k, bg_var, narrow, ds) if narrow else {}
    wl = [w for w in wins if w[0] in wide]
    if wl:
        sw = cfg["shift_scan_wide"]
        out.update(window_shifts(c, S, B, k, bg_var, wl, np.arange(sw["lo"], sw["hi"] + 1e-9, sw["step"])))
    return out

def check_set_boot(c, S, B, k, bg_var, cfg, names, prefix):
    th = cfg["thresholds"]
    ordered = sorted(names, key=lambda n: cfg["windows"][n][0])
    wins = [(n, cfg["windows"][n][0], cfg["windows"][n][1]) for n in ordered]
    scan = cfg["shift_scan"]
    ds = np.arange(scan["lo"], scan["hi"] + 1e-9, scan["step"])
    g = cfg["joint_grid"]
    mask = np.zeros(len(B), dtype=bool)
    for n, lo, hi in wins:
        mask[lo:hi+1] = True
    
    results = []
    
    # (a) window_quality checks
    r, v = sr.apply(c, S, B, k, bg_var)
    for n, lo, hi in wins:
        q = window_quality(r, v, lo, hi)
        z_dip = float(q['z_dip'])
        z_chi2 = float(q['z_chi2'])
        results.append((prefix + n + ":z_dip", z_dip, float(th["z_dip"]), abs(z_dip) <= th["z_dip"]))
        results.append((prefix + n + ":z_chi2", z_chi2, float(th["z_chi2"]), z_chi2 <= th["z_chi2"]))
    
    # (b) Point estimates
    sh0 = shifts_cfg(c, S, B, k, bg_var, wins, ds, cfg)
    j0 = joint_grid(c, S, B, k, bg_var, mask, cfg["pivot_joint"], g)
    
    # (c) Bootstrap
    rng = np.random.default_rng(cfg["seed"])
    boot_deltas = {n: [] for n in ordered}
    boot_dists = []
    boot_jd = []
    boot_js = []
    boot_ja = []
    
    for b in range(cfg["n_boot"]):
        S_r, B_r = resample(S, B, bg_var, rng)
        shb = shifts_cfg(c, S_r, B_r, k, bg_var, wins, ds, cfg)
        jb = joint_grid(c, S_r, B_r, k, bg_var, mask, cfg["pivot_joint"], g)
        
        for n in ordered:
            boot_deltas[n].append(shb[n][0])
        
        if len(ordered) >= 2:
            first_n = ordered[0]
            last_n = ordered[-1]
            boot_dists.append(shb[last_n][0] - shb[first_n][0])
        
        boot_jd.append(jb[0])
        boot_js.append(jb[1])
        boot_ja.append(jb[2])
    
    # (d) Sigmas
    sd_win = {}
    for n in ordered:
        if len(boot_deltas[n]) > 1:
            std_val = np.std(boot_deltas[n], ddof=1)
        else:
            std_val = 0.0
        sd_win[n] = float(np.sqrt(std_val**2 + scan["step"]**2 / 12))
    
    if len(boot_dists) > 1:
        std_dist = np.std(boot_dists, ddof=1)
    else:
        std_dist = 0.0
    sd_dist = float(np.sqrt(std_dist**2 + 2 * scan["step"]**2 / 12))
    
    if len(boot_jd) > 1:
        std_jd = np.std(boot_jd, ddof=1)
    else:
        std_jd = 0.0
    sd_jd = float(np.sqrt(std_jd**2 + g["d_step"]**2 / 12))
    
    if len(boot_js) > 1:
        std_js = np.std(boot_js, ddof=1)
    else:
        std_js = 0.0
    sd_js = float(np.sqrt(std_js**2 + g["s_step"]**2 / 12))
    
    if len(boot_ja) > 1:
        std_ja = np.std(boot_ja, ddof=1)
    else:
        std_ja = 0.0
    sd_ja = float(std_ja)
    
    # (e) Checks
    for n in ordered:
        delta0 = sh0[n][0]
        at_edge = sh0[n][1]
        if at_edge:
            z_val = float('inf')
            passed = False
        else:
            z_val = float(delta0 / sd_win[n])
            passed = abs(z_val) <= th["z_shift_win"]
        results.append((prefix + n + ":shift_z", z_val, float(th["z_shift_win"]), passed))
    
    if len(ordered) >= 2:
        first_n = ordered[0]
        last_n = ordered[-1]
        delta0_first = sh0[first_n][0]
        delta0_last = sh0[last_n][0]
        edge_first = sh0[first_n][1]
        edge_last = sh0[last_n][1]
        if edge_first or edge_last:
            z_val = float('inf')
            passed = False
        else:
            z_val = float((delta0_last - delta0_first) / sd_dist)
            passed = abs(z_val) <= th["z_dist"]
        results.append((prefix + "dist:" + first_n + "->" + last_n, z_val, float(th["z_dist"]), passed))
    
    j0_d, j0_s, j0_a, j0_edge = j0
    if j0_edge:
        z_jd = float('inf')
        z_js = float('inf')
        passed_jd = False
        passed_js = False
    else:
        z_jd = float(j0_d / sd_jd)
        z_js = float(j0_s / sd_js)
        passed_jd = abs(z_jd) <= th["z_joint_shift"]
        passed_js = abs(z_js) <= th["z_joint_scale"]
    results.append((prefix + "joint:shift_z", z_jd, float(th["z_joint_shift"]), passed_jd))
    results.append((prefix + "joint:scale_z", z_js, float(th["z_joint_scale"]), passed_js))
    
    if sd_ja == 0:
        z_ja = float('inf')
        passed_ja = False
    else:
        z_ja = float((j0_a - 1) / sd_ja)
        passed_ja = abs(z_ja) <= th["z_joint_amp"]
    results.append((prefix + "joint:amp_z", z_ja, float(th["z_joint_amp"]), passed_ja))
    
    # (f) Diagnostic
    diag_parts = ["[BOOT] %s n_boot=%d" % (prefix, cfg["n_boot"])]
    for n in ordered:
        diag_parts.append("%s: d=%+.3f sd=%.3f" % (n, sh0[n][0], sd_win[n]))
    diag_parts.append("joint: d=%+.3f sd=%.3f s=%+.2e sd=%.2e a=%.4f sd=%.4f" % (j0_d, sd_jd, j0_s, sd_js, j0_a, sd_ja))
    print(" ".join(diag_parts))
    
    return results
