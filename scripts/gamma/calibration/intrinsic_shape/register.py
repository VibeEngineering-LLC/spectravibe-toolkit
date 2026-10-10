"""register -- регистрация шкалы образца на шкалу фона глобальным полиномом h(x)=c0+c1*x+c2*x^2.
Глобальный поиск: сетка (c1, c2) с жёстким якорем + многостарт локального спуска; фон размывается эквивалентным
ядром перебинировки ВНУТРИ chi2; якорь (для LaBr3 -- пик ~37 кэВ) с сигмой из разброса по окнам подгонки.
Перенесено из прототипа labr3_selfcal/shape_reg.py (сгенерирован gen_code.py, qwen3.8:27b), 2026-10-10.
"""
import json, pathlib
import numpy as np
from scipy.optimize import least_squares
from .core import load_pair, rebin_conserve_var, poly_eval, N_CH
from .peakfit import fit_peak


U = 1000.0
SC = np.array([1.0, U, U * U])
TEST_CH = (29.0, 337.0, 845.0, 1340.0, 1610.0)


def blur_bg(mapping, b, var_b=None):
    """
    Source: make_net_xml.py:35-44 (generalised: any length, optional input variance)
    """
    b = np.asarray(b, float)
    x = np.arange(len(b), dtype=float)
    d = np.asarray(mapping(x), float) - x
    f = d - np.floor(d)
    v = f * (1 - f)
    w0 = 1 - v
    w1 = v / 2
    vb = b if var_b is None else np.asarray(var_b, float)

    out = w0 * b
    out[1:] += w1[:-1] * b[:-1]
    out[:-1] += w1[1:] * b[1:]

    var = w0**2 * vb
    var[1:] += w1[:-1]**2 * vb[:-1]
    var[:-1] += w1[1:]**2 * vb[1:]

    return out, var


def scan_centroid(counts, x0, w0, hw):
    """
    Source: make_net_xml.py:58-71 (peak branch only)
    """
    cs = []
    es = []
    # правка Claude: при hw=6 множители 0.7 и 0.85 дают одно и то же окно (5) -> дубликат занижал разброс; берём различные окна
    for h in sorted({max(5, int(round(hw * f))) for f in (0.7, 0.85, 1.0, 1.15, 1.3)}):
        r = fit_peak(counts, x0, w0, half_window=h)
        if r is not None:
            cs.append(r['centroid'])
            es.append(r['centroid_err'])
    
    cs = np.array(cs)
    es = np.array(es)
    ok = (np.abs(cs - np.median(cs)) < 3.0) & (es < 3.0)
    if ok.sum() < 2:
        raise RuntimeError("Not enough valid centroids in scan_centroid")
    
    cs = cs[ok]
    es = es[ok]
    
    return {
        'mean': float(cs.mean()),
        'spread': float(cs.std(ddof=1)),
        'fit_err': float(es.mean()),
        'sigma': float(np.hypot(cs.std(ddof=1), es.mean())),
        'centroids': [float(v) for v in cs],
        'n_ok': int(ok.sum()),
        'n_tot': len(ok)
    }


def anchor37(bg_counts, sm_counts, x_bg=28.9, x_sm=28.6, w0=3.3, hw=6, verbose=True):
    a = scan_centroid(bg_counts, x_bg, w0, hw)
    s = scan_centroid(sm_counts, x_sm, w0, hw)
    
    xs = s['mean']
    xb = a['mean']
    sigma = float(np.hypot(a['sigma'], s['sigma']))
    
    if verbose:
        print(f"[ANCHOR] bg: mean={a['mean']:.4f} spread={a['spread']:.4f} fit_err={a['fit_err']:.4f} n={a['n_ok']}/{a['n_tot']} centroids={a['centroids']}")
        print(f"[ANCHOR] sm: mean={s['mean']:.4f} spread={s['spread']:.4f} fit_err={s['fit_err']:.4f} n={s['n_ok']}/{s['n_tot']} centroids={s['centroids']}")
        print(f"[ANCHOR] h({xs:.4f}) = {xb:.4f} +- {sigma:.4f} (sigma from window spread, not fit error)")
        
    return xs, xb, sigma


class Problem:
    def __init__(self, S, B, k, mask, anchor, bg_var=None, blur=True):
        self.S = np.asarray(S, float)
        self.B = np.asarray(B, float)
        self.k = float(k)
        self.mask = np.asarray(mask, bool)
        self.anchor = tuple(anchor)
        self.bg_var = None if bg_var is None else np.asarray(bg_var, float)
        self.blur = blur
        self.x_edges = np.arange(len(self.S) + 1) - 0.5
        self.nres = int(self.mask.sum()) + 1

    def resid(self, c):
        e = poly_eval(c, self.x_edges)
        if np.any(np.diff(e) <= 0):
            return None
        s, vs = rebin_conserve_var(self.S, e, len(self.B))
        if self.blur:
            bb, vb = blur_bg(lambda x: poly_eval(c, x), self.B, self.bg_var)
        else:
            bb = self.B
            vb = self.B if self.bg_var is None else self.bg_var
        return (s - self.k * bb) / np.sqrt(np.maximum(vs + self.k ** 2 * vb, 1.0))

    def resid_vec(self, c):
        r = self.resid(c)
        if r is None:
            return np.full(self.nres, 1e3)
        a = (poly_eval(c, self.anchor[0]) - self.anchor[1]) / self.anchor[2]
        return np.concatenate([r[self.mask], [a]])

    def chi2(self, c):
        return float(np.sum(self.resid_vec(c) ** 2))

    def ndf(self, npar=3):
        return self.nres - npar


def grid_search(prob, c1s, c2s):
    """
    Generalised from verify/v5_grid.py:9-13 (hard anchor c0), :18-23 (grid), :27-31 (3x3 local minima)
    """
    xs, xb = prob.anchor[0], prob.anchor[1]
    G = np.full((len(c1s), len(c2s)), np.inf)
    
    for i, c1 in enumerate(c1s):
        for j, c2 in enumerate(c2s):
            c0 = xb - c1 * xs - c2 * xs ** 2
            G[i, j] = prob.chi2([c0, c1, c2])
            
    starts = []
    for i in range(len(c1s)):
        for j in range(len(c2s)):
            if not np.isfinite(G[i, j]):
                continue
            # 3x3 neighborhood
            i_min = max(0, i - 1)
            i_max = min(len(c1s) - 1, i + 1)
            j_min = max(0, j - 1)
            j_max = min(len(c2s) - 1, j + 1)
            
            neigh = G[i_min:i_max+1, j_min:j_max+1]
            if G[i, j] == np.min(neigh):
                c0 = xb - c1s[i] * xs - c2s[j] * xs ** 2
                starts.append((G[i, j], [c0, c1s[i], c2s[j]]))
                
    starts.sort(key=lambda x: x[0])
    return G, starts


def local_fit(prob, c_start, fixed=None, max_nfev=300):
    q0 = np.asarray(c_start, float) * SC
    if fixed is None:
        free = [0, 1, 2]
    else:
        idx, value = fixed
        free = [i for i in range(3) if i != idx]
        q0[idx] = value * SC[idx]
        
    # правка Claude (#DBG-1): относительный diff_step=1e-4 на q1~1000 давал шаг 0.1 кан-экв. (h(1610) +0.16 кан),
    # сравнимый с рябью chi2 по c1 -> спуск останавливался не в минимуме (точка verify_v6 была ниже на 2.0).
    # Оптимизируем смещение p = q - q0 (каналы при x=1000), абсолютный шаг ~0.005 кан.
    def fun(p):
        q = q0.copy()
        q[free] = q0[free] + p
        return prob.resid_vec(q / SC)

    res = least_squares(fun, np.zeros(len(free)), method="trf", diff_step=5e-3, max_nfev=max_nfev)
    q = q0.copy()
    q[free] = q0[free] + res.x
    c = q / SC
    return c, float(np.sum(res.fun ** 2))


def hvals(c, chans=TEST_CH):
    return np.array([poly_eval(c, x) for x in chans])


def multistart(prob, starts, n_max=20, tol_ch=0.02):
    results = []
    for g, c in starts[:n_max]:
        c_fit, chi2 = local_fit(prob, c)
        results.append({
            'c': c_fit,
            'chi2': chi2,
            'h': hvals(c_fit),
            'start_chi2': g
        })
        
    results.sort(key=lambda x: x['chi2'])
    kept = []
    for r in results:
        is_dup = False
        for k in kept:
            if np.max(np.abs(r['h'] - k['h'])) < tol_ch:
                is_dup = True
                break
        if not is_dup:
            kept.append(r)
            
    return kept


def covariance(prob, c, step=0.01):
    q = c * SC
    J = np.zeros((prob.nres, 3))
    for i in range(3):
        dq = np.zeros(3)
        dq[i] = step
        J[:, i] = (prob.resid_vec((q + dq) / SC) - prob.resid_vec((q - dq) / SC)) / (2 * step)
        
    chi2 = prob.chi2(c)
    chi2ndf = chi2 / prob.ndf()
    cov_q = np.linalg.pinv(J.T @ J) * max(1.0, chi2ndf)
    cov = cov_q / np.outer(SC, SC)
    return cov, np.sqrt(np.diag(cov)), chi2ndf


def profile_grid(prob, idx, values, c1s, c2s, n_ref=3):
    """Правка Claude: профиль chi2 по параметру idx через сетку по двум другим (c0 из жёсткого якоря, если он свободен)
    + локальное уточнение с фиксированным idx из n_ref лучших ячеек; спуск из минимумов (profile) застревал."""
    xs, xb = prob.anchor[0], prob.anchor[1]
    out = np.zeros(len(values))
    for n, v in enumerate(values):
        cells = []
        for c1 in ([v] if idx == 1 else c1s):
            for c2 in ([v] if idx == 2 else c2s):
                c0 = v if idx == 0 else xb - c1 * xs - c2 * xs ** 2
                cells.append((prob.chi2([c0, c1, c2]), [c0, c1, c2]))
        cells.sort(key=lambda t: t[0])
        out[n] = min([cells[0][0]] + [local_fit(prob, cc, fixed=(idx, v), max_nfev=150)[1] for _, cc in cells[:n_ref]])
    return out


def profile(prob, minima, idx, values, n_starts=5):
    out = np.zeros(len(values))
    for i, v in enumerate(values):
        best = np.inf
        for m in minima[:n_starts]:
            c_s = m['c'].copy()
            c_s[idx] = v
            c_f, chi2 = local_fit(prob, c_s, fixed=(idx, v), max_nfev=150)
            best = min(best, chi2)
        out[i] = best
    return out


def mapping_spread(minima, dchi2=4.0, chans=TEST_CH):
    c2min = minima[0]['chi2']
    sel = [m for m in minima if m['chi2'] - c2min <= dchi2]
    H = np.array([hvals(m['c'], chans) for m in sel])
    return {
        'n': len(sel),
        'chans': list(chans),
        'lo': H.min(0),
        'hi': H.max(0),
        'range': H.max(0) - H.min(0),
        'std': (H.std(0, ddof=1) if len(sel) > 1 else np.zeros(len(chans)))
    }


def coarse_c1(S, B, anchor, window=(780, 900)):
    """Грубый наклон по максимуму сильной структуры в окне (для LaBr3 -- сумм-пик La-138) и якорю."""
    lo, hi = window
    kernel = np.ones(7) / 7
    sS = np.convolve(S, kernel, mode="same")
    sB = np.convolve(B, kernel, mode="same")
    pS = lo + np.argmax(sS[lo:hi])
    pB = lo + np.argmax(sB[lo:hi])
    return (pB - anchor[1]) / (pS - anchor[0])


def register(S, B, k, mask, anchor, bg_var=None, c1_center=None, c1_half=0.008, c1_step=0.0002, c2_lo=-1.0e-5, c2_hi=2.0e-6, c2_step=2.0e-7, n_ms=20, verbose=True, coarse_window=(780, 900)):
    prob = Problem(S, B, k, mask, anchor, bg_var)
    if c1_center is None:
        c1_center = coarse_c1(S, B, anchor, coarse_window)
        
    c1s = np.arange(c1_center - c1_half, c1_center + c1_half + 1e-12, c1_step)
    c2s = np.arange(c2_lo, c2_hi + 1e-15, c2_step)
    
    G, starts = grid_search(prob, c1s, c2s)
    minima = multistart(prob, starts, n_max=n_ms)
    # правка Claude: дополнительная мелкая сетка (шаг /10) вокруг лучших 3 минимумов + спуск -> добавляется в список минимумов
    xs_, xb_ = anchor[0], anchor[1]
    extra = []
    for m in minima[:3]:
        fc1 = m['c'][1] + np.arange(-10, 11) * c1_step / 10
        fc2 = m['c'][2] + np.arange(-10, 11) * c2_step / 10
        _, st = grid_search(prob, fc1, fc2)
        extra += st[:3]
    minima = multistart(prob, [(m['chi2'], m['c']) for m in minima] + extra, n_max=len(minima) + len(extra))
    best = minima[0]
    
    cov, perr, chi2ndf = covariance(prob, best['c'])
    spread = mapping_spread(minima)
    
    ib, jb = np.unravel_index(np.argmin(G), G.shape)
    on_edge = ib in (0, len(c1s) - 1) or jb in (0, len(c2s) - 1)
    
    if verbose:
        print(f"[GRID] {len(c1s)} x {len(c2s)} cells, grid local minima: {len(starts)}, best grid chi2={G[ib, jb]:.1f} at c1={c1s[ib]:.5f} c2={c2s[jb]:.3e}, on_edge={on_edge}")
        print(f"[MIN] distinct local minima after multistart: {len(minima)}")
        for i, m in enumerate(minima):
            h_str = ", ".join([f"{v:.3f}" for v in m['h']])
            print(f"[MIN] #{i} chi2={m['chi2']:.2f} dchi2={m['chi2'] - best['chi2']:.2f} c=[{m['c'][0]:.6f}, {m['c'][1]:.7f}, {m['c'][2]:.4e}] h(29,337,845,1340,1610)={h_str}")
            
    return {
        'c': best['c'],
        'chi2': best['chi2'],
        'ndf': prob.ndf(),
        'chi2ndf': chi2ndf,
        'cov': cov,
        'perr': perr,
        'minima': minima,
        'spread': spread,
        'n_grid_minima': len(starts),
        'grid_on_edge': bool(on_edge),
        'prob': prob,
        'c1s': c1s,
        'c2s': c2s,
        'G': G
    }


def apply(c, S, B, k, bg_var=None, blur=True):
    """
    Generalised from make_net_xml.py:47-54 (net_for)
    """
    e = poly_eval(c, np.arange(len(S) + 1) - 0.5)
    s, vs = rebin_conserve_var(np.asarray(S, float), e, len(B))
    if blur:
        bb, vb = blur_bg(lambda x: poly_eval(c, x), np.asarray(B, float), bg_var)
    else:
        bb = B
        vb = B if bg_var is None else bg_var
    return s - k * bb, np.maximum(vs + k ** 2 * vb, 1.0)


def load_inputs(bg_path, sample_path):
    """S, B, k = live_S/live_B, live_S, live_B; последний канал (переполнение) обнулён."""
    bg, sm = load_pair(bg_path, sample_path)
    B = bg['counts'].copy()
    S = sm['counts'].copy()
    k = sm['live'] / bg['live']
    return S, B, k, sm['live'], bg['live']


def report(S, B, k, anchor, mask, out_json=None, profiles=True, **reg_kw):
    """Регистрация + печать h, ковариации, минимумов, профилей, «требует толкования»; JSON -- в out_json."""
    r = register(S, B, k, mask, anchor, **reg_kw)
    
    c = r['c']
    perr = r['perr']
    cov = r['cov']
    
    print(f"[H] c = [{c[0]:.7g}, {c[1]:.9g}, {c[2]:.6g}]")
    print(f"[H] perr = [{perr[0]:.3g}, {perr[1]:.3g}, {perr[2]:.3g}]")
    print(f"[H] chi2={r['chi2']:.2f} ndf={r['ndf']} chi2/ndf={r['chi2ndf']:.4f}")
    
    for i in range(3):
        row = " ".join([f"{cov[i, j]:.4e}" for j in range(3)])
        print(f"[COV] {row}")
        
    corr = cov / np.outer(perr, perr)
    for i in range(3):
        row = " ".join([f"{corr[i, j]:.4f}" for j in range(3)])
        print(f"[CORR] {row}")
        
    z = (poly_eval(c, anchor[0]) - anchor[1]) / anchor[2]
    print(f"[H] anchor pull z={z:+.2f}")
    
    prof_lower = []
    for idx in ((0, 1, 2) if profiles else ()):
        within = [m for m in r['minima'] if m['chi2'] - r['chi2'] <= 25.0]
        half = max(6.0 * r['perr'][idx], 1.2 * max(abs(m['c'][idx] - r['c'][idx]) for m in within))
        values = r['c'][idx] + np.linspace(-half, half, 25)
        P = profile_grid(r['prob'], idx, values, r['c1s'], r['c2s'])
        print(f"[PROF] c{idx} (dchi2 к глобальному минимуму {r['chi2']:.2f}; min профиля {P.min() - r['chi2']:+.2f})")
        for v, p in zip(values, P):
            print(f"[PROF] c{idx}={v:.7g} dchi2={p - r['chi2']:.2f}")
        if P.min() < r['chi2'] - 0.5:
            prof_lower.append(f"профиль c{idx} нашёл точку ниже глобального минимума на {r['chi2'] - P.min():.2f} — минимум не достигнут")
        # правка Claude: широкий профиль (шаг ~ разброс минимумов / 10) грубее 1 sigma -> добавлен мелкий профиль +-4 perr
        vf = r['c'][idx] + np.linspace(-4.0, 4.0, 17) * r['perr'][idx]
        Pf = profile_grid(r['prob'], idx, vf, r['c1s'], r['c2s'])
        print(f"[PROF] c{idx} мелкий (+-4 perr): " + " ".join(f"{p - r['chi2']:.2f}" for p in Pf))
        if Pf.min() < r['chi2'] - 0.5:
            prof_lower.append(f"мелкий профиль c{idx} нашёл точку ниже глобального минимума на {r['chi2'] - Pf.min():.2f}")
        d = Pf - r['chi2'] - 1.0
        lo_x = hi_x = None
        for a in range(8, 0, -1):
            if d[a] <= 0 < d[a - 1]:
                lo_x = vf[a - 1] + (vf[a] - vf[a - 1]) * d[a - 1] / (d[a - 1] - d[a]); break
        for a in range(8, 16):
            if d[a] <= 0 < d[a + 1]:
                hi_x = vf[a] + (vf[a + 1] - vf[a]) * d[a] / (d[a] - d[a + 1]); break
        if lo_x is not None and hi_x is not None:
            print(f"[PROF] c{idx} мелкий: dchi2=1 пересечения [{lo_x:.7g}, {hi_x:.7g}], полуширина {(hi_x - lo_x) / 2:.3g} (Jacobian perr {r['perr'][idx]:.3g})")
        else:
            print(f"[PROF] c{idx} мелкий: пересечение dchi2=1 не найдено в +-4 perr (профиль не параболический)")
            prof_lower.append(f"профиль c{idx}: нет пересечения dchi2=1 в +-4 perr")
        mask_1sig = P - r['chi2'] <= 1.0
        if np.any(mask_1sig):
            v_lo = values[mask_1sig].min()
            v_hi = values[mask_1sig].max()
            hw = (v_hi - v_lo) / 2
            print(f"[PROF] c{idx} dchi2<=1 interval: [{v_lo:.7g}, {v_hi:.7g}], half-width {hw:.3g} (Jacobian perr {r['perr'][idx]:.3g})")
            
    sp = r['spread']
    print(f"[SYST] minima within dchi2<=4: {sp['n']}")
    for i, ch in enumerate(sp['chans']):
        g = np.array([1.0, ch, ch**2])
        stat = np.sqrt(g @ cov @ g)
        print(f"[SYST] h({ch:g}): min={sp['lo'][i]:.3f} max={sp['hi'][i]:.3f} range={sp['range'][i]:.3f} std={sp['std'][i]:.3f} stat={stat:.3f}")
        
    print("[ТРЕБУЕТ ТОЛКОВАНИЯ]")
    items = list(prof_lower)
    if sp['n'] > 1:
        items.append(f"несколько минимумов в пределах dchi2<=4: {sp['n']}; систематика отображения = их разброс")
    if r['grid_on_edge']:
        items.append("минимум сетки на краю диапазона — расширить сетку")
    if abs(r['chi2ndf'] - 1) > 0.1:
        items.append(f"chi2/ndf={r['chi2ndf']:.3f} далеко от 1")
    if abs(z) > 2:
        items.append(f"якорь ~37 кэВ тянет решение: z={z:+.2f}")
    for i, ch in enumerate(sp['chans']):
        if sp['range'][i] > 1.0:
            items.append(f"разброс h > 1 канала на канале {ch:g}")
            
    if items:
        for it in items:
            print(f"  {it}")
    else:
        print("  (нет пунктов)")
        
    out_data = {
        'c': [float(x) for x in c],
        'perr': [float(x) for x in perr],
        'cov': [[float(x) for x in row] for row in cov],
        'chi2': float(r['chi2']),
        'ndf': int(r['ndf']),
        'chi2ndf': float(r['chi2ndf']),
        'anchor': [float(x) for x in anchor],
        'k': float(k),
        'minima': [{'c': [float(x) for x in m['c']], 'chi2': float(m['chi2'])} for m in r['minima']],
        'spread': {
            'n': int(sp['n']),
            'chans': [float(x) for x in sp['chans']],
            'lo': [float(x) for x in sp['lo']],
            'hi': [float(x) for x in sp['hi']],
            'range': [float(x) for x in sp['range']],
            'std': [float(x) for x in sp['std']]
        },
        'n_grid_minima': int(r['n_grid_minima']),
        'grid_on_edge': bool(r['grid_on_edge'])
    }
    
    if out_json:
        pathlib.Path(out_json).parent.mkdir(parents=True, exist_ok=True)
        with open(out_json, "w", encoding="utf-8") as f:
            json.dump(out_data, f, ensure_ascii=False, indent=1)
    return r, out_data
