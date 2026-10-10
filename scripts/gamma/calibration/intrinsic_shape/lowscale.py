"""
Низкоэнергетическая ветвь (< E_join = 100 кэВ) энергетической шкалы сцинтиллятора LaBr3(Ce),
стыкованная непрерывно с основной полиномиальной шкалой.
Модель: однопараметрическое падение отклика (non-proportionality) + ступень на K-крае La.
Источники:
- K-край La 38.925 кэВ -- X-ray Data Booklet (LBNL/CXRO);
- величина ступени отклика на K-крае La 1.7 % -- синхротронные измерения LaBr3:Ce, arXiv:1010.5887 (группа TU Delft); падение отклика у K-краёв -- Feng et al., arXiv:2312.16658;
- общее падение фотонного отклика 15 % в диапазоне 9-100 кэВ -- arXiv:1010.5887 (используется только как априорный масштаб параметра s);
- форма кривой между точками (квадратичная по ln E) -- минимальная гладкая параметризация автора модуля, НЕ из источника; альтернативная форма (линейная по ln E) даёт оценку модельной неопределённости;
- гамма-эквивалент собственного пика ~37 кэВ -- параметр (литература 35.5-36 кэВ на гамма-шкалах: Cheng 2020 arXiv:2007.07552, Zeng 2016 arXiv:1608.01109; 37.45 кэВ при фотонной калибровке 32-40 кэВ: Li 2011 arXiv:1112.5815).
Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-10-10
"""

import numpy as np
from .core import poly_eval, poly_deriv

E_K_LA = 38.925      # кэВ, K-край La
DK_STEP = 0.017      # относительная ступень отклика на K-крае La
S_PRIOR = 0.15       # падение отклика 9-100 кэВ
E_JOIN = 100.0       # кэВ, стык с основной шкалой
E37_GAMMA = 35.75    # кэВ, гамма-эквивалент собственного пика (литература 35.5-36)
DE37 = 0.5           # кэВ
RANGES = ((17.0, 37.0), (37.0, 60.0), (60.0, 100.0), (100.0, 300.0), (300.0, 1500.0), (1500.0, 3000.0))


def response(E, s, shape="quad", dK=DK_STEP, e_join=E_JOIN, e_k=E_K_LA):
    """
    Вычисляет относительный отклик R(E) для энергии E.
    R(E) = Rs * (1 + dK) для E < e_k, иначе Rs.
    Rs = 1 - s * u^2 (quad) или 1 - s * u (lin), где u = ln(e_join / E) для E < e_join, иначе 0.
    """
    E = np.asarray(E, dtype=float)
    scalar_input = E.ndim == 0
    if scalar_input:
        E = E.reshape(1)
    
    # u(E) = ln(E_join / E) for E < E_join, else 0
    u = np.zeros_like(E)
    mask_low = E < e_join
    u[mask_low] = np.log(e_join / E[mask_low])
    
    # Rs
    if shape == "quad":
        Rs = 1.0 - s * u**2
    elif shape == "lin":
        Rs = 1.0 - s * u
    else:
        raise ValueError("Неизвестная форма: %s" % shape)
    
    # Staircase at K-edge
    R = np.where(E < e_k, Rs * (1.0 + dK), Rs)
    
    if scalar_input:
        return float(R[0])
    return R


def g_of_e(E, s, **kw):
    """
    g(E) = E * response(E, s, **kw).
    """
    E = np.asarray(E, dtype=float)
    scalar_input = E.ndim == 0
    if scalar_input:
        E = E.reshape(1)
    
    R = response(E, s, **kw)
    g = E * R
    
    if scalar_input:
        return float(g[0])
    return g


def e_of_l(L, s, shape="quad", dK=DK_STEP, e_join=E_JOIN, e_k=E_K_LA, e_min=5.0):   # правка Claude: сетка от 5 кэВ (порог прибора ~15 кэВ); от 1 кэВ модель немонотонна уже при s>0.085
    """
    Обратная функция e_of_l(L) = g^{-1}(L).
    Использует плотную сетку для интерполяции.
    """
    L = np.asarray(L, dtype=float)
    scalar_input = L.ndim == 0
    if scalar_input:
        L = L.reshape(1)
    
    # Build grids
    EA = np.geomspace(e_min, e_k, 4000, endpoint=False)
    EB = np.geomspace(e_k, e_join, 4000)
    
    gA = g_of_e(EA, s, shape=shape, dK=dK, e_join=e_join, e_k=e_k)
    gB = g_of_e(EB, s, shape=shape, dK=dK, e_join=e_join, e_k=e_k)
    
    # Check monotonicity
    if np.any(np.diff(gA) <= 0) or np.any(np.diff(gB) <= 0):
        raise ValueError("модель отклика немонотонна при s=%g" % s)
    
    E_out = np.zeros_like(L)
    
    for i in range(len(L)):
        Li = L[i]
        if Li >= e_join:
            E_out[i] = Li
        elif Li >= gB[0]:
            E_out[i] = np.interp(Li, gB, EB)
        else:
            # Interpolate on A
            if Li < gA[0]:
                # Linear extrapolation: E = L / response(e_min)
                R_min = response(e_min, s, shape=shape, dK=dK, e_join=e_join, e_k=e_k)
                E_out[i] = Li / R_min
            else:
                E_out[i] = np.interp(Li, gA, EA)
    
    if scalar_input:
        return float(E_out[0])
    return E_out


def branch_energy(C, ch, s, **kw):
    """
    Энергия ветви: E(ch) = e_of_l(poly_eval(C, ch), s, **kw).
    """
    L = poly_eval(C, ch)
    return e_of_l(L, s, **kw)


def fit_s(C, points, shape="quad", dK=DK_STEP, e37=E37_GAMMA, de37=DE37, use37=True, s_grid=None):
    """
    Подгонка параметра s по реперам.
    """
    if s_grid is None:
        s_grid = np.linspace(-0.3, 0.3, 1201)   # правка Claude: s -- поправка ОТНОСИТЕЛЬНО основного полинома, знак данными
    
    # Prepare points
    used_points = []
    for p in points:
        if p["kind"] == "intrinsic37":
            if use37:
                used_points.append({
                    "name": p["name"],
                    "E": e37,
                    "dE": de37,
                    "ch": p["ch"],
                    "ch_err": p["ch_err"]
                })
        else:
            if p["E"] < E_JOIN:
                used_points.append(p)
    
    n_used = len(used_points)
    if n_used == 0:
        raise ValueError("нет реперов ниже %g кэВ" % E_JOIN)
    
    # Precompute ch and ch_err
    chs = np.array([p["ch"] for p in used_points])
    ch_errs = np.array([p["ch_err"] for p in used_points])
    Es = np.array([p["E"] for p in used_points])
    dEs = np.array([p["dE"] for p in used_points])
    
    # Slopes
    slopes = poly_deriv(C, chs)
    
    # Compute chi2 for each s in grid
    chi2_vals = np.full_like(s_grid, np.nan)
    valid_mask = np.zeros_like(s_grid, dtype=bool)
    
    for i, s_val in enumerate(s_grid):
        try:
            E_branch = branch_energy(C, chs, s_val, shape=shape, dK=dK)
            res = Es - E_branch
            sigma_sq = dEs**2 + (slopes * ch_errs)**2
            # Avoid division by zero
            sigma_sq = np.maximum(sigma_sq, 1e-12)
            chi2 = np.sum((res**2) / sigma_sq)
            chi2_vals[i] = chi2
            valid_mask[i] = True
        except Exception:
            pass
    
    if not np.any(valid_mask):
        raise ValueError("нет валидных значений s")
    
    # Find minimum
    valid_indices = np.where(valid_mask)[0]
    chi2_valid = chi2_vals[valid_indices]
    min_idx_in_valid = np.argmin(chi2_valid)
    min_idx = valid_indices[min_idx_in_valid]
    
    s_best = s_grid[min_idx]
    chi2_min = chi2_vals[min_idx]
    
    # Refine by parabola
    if 0 < min_idx < len(s_grid) - 1:
        x0 = s_grid[min_idx - 1]
        x1 = s_grid[min_idx]
        x2 = s_grid[min_idx + 1]
        y0 = chi2_vals[min_idx - 1]
        y1 = chi2_vals[min_idx]
        y2 = chi2_vals[min_idx + 1]
        
        # Parabola through (x0,y0), (x1,y1), (x2,y2)
        # y = a*x^2 + b*x + c
        # Solve for a, b, c
        A = np.array([
            [x0**2, x0, 1],
            [x1**2, x1, 1],
            [x2**2, x2, 1]
        ])
        b_vec = np.array([y0, y1, y2])
        try:
            coeffs = np.linalg.solve(A, b_vec)
            a, b, c = coeffs
            if a > 0:
                s_refined = -b / (2 * a)
                if s_grid[min_idx - 1] < s_refined < s_grid[min_idx + 1]:
                    s_best = s_refined
                    # Recompute chi2 at refined s
                    try:
                        E_branch = branch_energy(C, chs, s_best, shape=shape, dK=dK)
                        res = Es - E_branch
                        sigma_sq = dEs**2 + (slopes * ch_errs)**2
                        sigma_sq = np.maximum(sigma_sq, 1e-12)
                        chi2_min = np.sum((res**2) / sigma_sq)
                    except Exception:
                        pass
        except np.linalg.LinAlgError:
            pass
    
    # Sigma_s
    sigma_s = np.nan
    # Search for chi2 = chi2_min + 1
    target = chi2_min + 1.0
    # Look left
    left_cross = None
    for i in range(min_idx, -1, -1):
        if not valid_mask[i]:
            continue
        if chi2_vals[i] >= target:
            left_cross = i
            break
    # Look right
    right_cross = None
    for i in range(min_idx, len(s_grid)):
        if not valid_mask[i]:
            continue
        if chi2_vals[i] >= target:
            right_cross = i
            break
    
    if left_cross is not None and right_cross is not None:
        # Interpolate
        s_left = s_grid[left_cross]
        s_right = s_grid[right_cross]
        # Linear interpolation of s where chi2 crosses target
        # Actually, we want the s values where chi2 = target.
        # We have points (s_left, chi2_left) and (s_right, chi2_right)
        # We want s such that chi2(s) = target.
        # Since chi2 is convex, we can interpolate.
        # But we have discrete points. Let's just take the average of the crossing points?
        # Better: interpolate s as a function of chi2.
        # s = s_left + (target - chi2_left) * (s_right - s_left) / (chi2_right - chi2_left)
        # This is not quite right because we have two separate crossings.
        # Let's find the s values more accurately.
        # For left side: between s_grid[left_cross-1] and s_grid[left_cross]
        # For right side: between s_grid[right_cross] and s_grid[right_cross+1]
        
        # Simplified: use the grid points where it crosses
        # Actually, the spec says "search on the grid". So we can just take the grid points.
        # But we need sigma_s. Let's assume the crossing is between the grid points.
        # We'll interpolate s for chi2 = target.
        
        # Left crossing: between s_grid[left_cross-1] and s_grid[left_cross]
        if left_cross > 0 and valid_mask[left_cross-1]:
            s_l1 = s_grid[left_cross-1]
            s_l2 = s_grid[left_cross]
            c_l1 = chi2_vals[left_cross-1]
            c_l2 = chi2_vals[left_cross]
            if c_l2 != c_l1:
                s_cross_left = s_l1 + (target - c_l1) * (s_l2 - s_l1) / (c_l2 - c_l1)
            else:
                s_cross_left = s_l1
        else:
            s_cross_left = s_grid[left_cross]
            
        # Right crossing: between s_grid[right_cross] and s_grid[right_cross+1]
        if right_cross < len(s_grid) - 1 and valid_mask[right_cross+1]:
            s_r1 = s_grid[right_cross]
            s_r2 = s_grid[right_cross+1]
            c_r1 = chi2_vals[right_cross]
            c_r2 = chi2_vals[right_cross+1]
            if c_r2 != c_r1:
                s_cross_right = s_r1 + (target - c_r1) * (s_r2 - s_r1) / (c_r2 - c_r1)
            else:
                s_cross_right = s_r1
        else:
            s_cross_right = s_grid[right_cross]
            
        sigma_s = (s_cross_right - s_cross_left) / 2.0
    
    # Scale by sqrt(max(1, chi2_min / max(ndf, 1)))
    ndf = n_used - 1
    scale = np.sqrt(max(1.0, chi2_min / max(ndf, 1.0)))
    sigma_s = sigma_s * scale
    
    return {
        "s": float(s_best),
        "sigma_s": float(sigma_s),
        "chi2": float(chi2_min),
        "ndf": int(ndf),
        "shape": shape,
        "dK": float(dK),
        "e37": float(e37),
        "used": [p["name"] for p in used_points],
        "n": int(n_used)
    }


def residuals(C, points, fit, e37=E37_GAMMA):
    """
    Вычисляет невязки до и после подгонки.
    """
    res_list = []
    for p in points:
        if p["kind"] == "intrinsic37":
            E = e37
        else:
            if p["E"] >= E_JOIN:
                continue
            E = p["E"]
        
        ch = p["ch"]
        before = E - poly_eval(C, ch)
        after = E - branch_energy(C, ch, fit["s"], shape=fit["shape"], dK=fit["dK"])
        
        res_list.append({
            "name": p["name"],
            "E": float(E),
            "ch": float(ch),
            "before": float(before),
            "after": float(after)
        })
    return res_list


def loo(C, points, shape="quad", dK=DK_STEP, e37=E37_GAMMA, de37=DE37, use37=True):
    """
    Leave-one-out кросс-валидация.
    """
    loo_dict = {}
    
    # Prepare used points
    used_points = []
    for p in points:
        if p["kind"] == "intrinsic37":
            if use37:
                used_points.append({
                    "name": p["name"],
                    "E": e37,
                    "dE": de37,
                    "ch": p["ch"],
                    "ch_err": p["ch_err"],
                    "kind": "photon"   # правка Claude: энергия уже заменена на e37; без kind fit_s падал -> LOO = nan
                })
        else:
            if p["E"] < E_JOIN:
                used_points.append(p)
    
    for i, p in enumerate(used_points):
        # Remove point i
        points_without = [used_points[j] for j in range(len(used_points)) if j != i]
        
        if len(points_without) == 0:
            loo_dict[p["name"]] = np.nan
            continue
        
        try:
            fit_without = fit_s(C, points_without, shape=shape, dK=dK, e37=e37, de37=de37, use37=use37)
            E_branch = branch_energy(C, p["ch"], fit_without["s"], shape=shape, dK=dK)
            loo_res = p["E"] - E_branch
            loo_dict[p["name"]] = float(loo_res)
        except Exception:
            loo_dict[p["name"]] = np.nan
    
    return loo_dict


def band(C, fit, ch_grid, alt_fit=None):
    """
    Полоса неопределённости шкалы.
    """
    ch_grid = np.asarray(ch_grid, dtype=float)
    
    E = branch_energy(C, ch_grid, fit["s"], shape=fit["shape"], dK=fit["dK"])
    
    # Stat
    if np.isnan(fit["sigma_s"]):
        stat = np.zeros_like(ch_grid)
    else:
        try:
            E_plus = branch_energy(C, ch_grid, fit["s"] + fit["sigma_s"], shape=fit["shape"], dK=fit["dK"])
            E_minus = branch_energy(C, ch_grid, fit["s"] - fit["sigma_s"], shape=fit["shape"], dK=fit["dK"])
            stat = np.abs(E_plus - E_minus) / 2.0
        except ValueError:   # правка Claude: s +- sigma вне области монотонности -- полоса по одному краю
            E_one = branch_energy(C, ch_grid, fit["s"] - np.sign(fit["s"]) * fit["sigma_s"], shape=fit["shape"], dK=fit["dK"])
            stat = np.abs(E_one - E)
    
    # Shape
    if alt_fit is not None:
        E_alt = branch_energy(C, ch_grid, alt_fit["s"], shape=alt_fit["shape"], dK=alt_fit["dK"])
        shape_band = np.abs(E - E_alt)
    else:
        shape_band = np.zeros_like(ch_grid)
    
    # K-edge
    kedge = np.where(E < E_K_LA, 0.5 * fit["dK"] * E, 0.0)
    
    sigma = np.sqrt(stat**2 + shape_band**2 + kedge**2)
    
    return {
        "ch": ch_grid,
        "E": E,
        "sigma": sigma,
        "stat": stat,
        "shape": shape_band,
        "kedge": kedge
    }


def composite_energy(C, ch, fit):
    """
    Композитная энергия: ветвь ниже E_JOIN, основная шкала выше.
    """
    return branch_energy(C, ch, fit["s"], shape=fit["shape"], dK=fit["dK"])


def fit_poly_xml(C, fit, n_ch, e_lo=17.0, e_hi=3000.0, max_deg=4, tol_fn=None, max_damage=1.0):
    """
    Приближение композитной шкалы полиномом степени <= 4.
    """
    ch_grid = np.arange(n_ch, dtype=float)
    E_model = composite_energy(C, ch_grid, fit)
    
    # Use channels where e_lo <= E_model <= e_hi
    mask = (E_model >= e_lo) & (E_model <= e_hi)
    ch_used = ch_grid[mask]
    E_used = E_model[mask]
    
    if len(ch_used) == 0:
        return {
            "coeffs": list(C),
            "deg": len(C) - 1,
            "max_dev": {},
            "fallback": True
        }
    
    # Weights
    # правка Claude: допуск -- полоса неопределённости ветви (tol_fn), не меньше 0.3 кэВ; выше 100 кэВ основная шкала точна
    tol = np.maximum(0.3, tol_fn(ch_used)) if tol_fn is not None else np.maximum(0.15, 0.002 * E_used)
    w = 1.0 / tol**2
    
    # Scaled variable t = ch / n_ch
    t = ch_used / n_ch
    
    # Try degrees from max_deg down to 2
    for d in range(max_deg, 1, -1):
        # Design matrix [t^0, t^1, ..., t^d]
        A = np.vander(t, N=d + 1, increasing=True)
        
        # Weighted least squares
        Aw = A * np.sqrt(w)[:, np.newaxis]
        bw = E_used * np.sqrt(w)
        
        try:
            b, residuals, rank, sv = np.linalg.lstsq(Aw, bw, rcond=None)
        except np.linalg.LinAlgError:
            continue
        
        # Convert to raw channel coefficients
        # E = sum_{j=0}^d b_j * (ch/n_ch)^j = sum_{j=0}^d (b_j / n_ch^j) * ch^j
        a = [b[j] / (n_ch ** j) for j in range(d + 1)]
        
        # Check becqmoni_ok
        ok, reason = becqmoni_ok(a, n_ch)
        if ok:
            # Compute max deviation per range
            dev = poly_eval(a, ch_used) - E_used
            max_dev = {}
            for r in RANGES:
                mask_r = (E_used >= r[0]) & (E_used < r[1])
                if np.any(mask_r):
                    max_dev["%g-%g" % r] = float(np.max(np.abs(dev[mask_r])))
                else:
                    max_dev["%g-%g" % r] = np.nan
            
            # правка Claude: полином, портящий шкалу выше 100 кэВ больше max_damage кэВ, не пишется (запасной -- основная)
            damage = max([v for k, v in max_dev.items() if float(k.split('-')[0]) >= E_JOIN and np.isfinite(v)] or [0.0])
            if damage > max_damage:
                return {"coeffs": list(C), "deg": len(C) - 1, "max_dev": max_dev, "fallback": True, "damage": float(damage), "deg_tried": d}
            return {
                "coeffs": a,
                "deg": d,
                "max_dev": max_dev,
                "fallback": False,
                "damage": float(damage)
            }
    
    # Fallback
    return {
        "coeffs": list(C),
        "deg": len(C) - 1,
        "max_dev": {},
        "fallback": True
    }


def becqmoni_ok(coeffs, n_ch):
    """
    Проверка полинома для целевой программы.
    """
    d = len(coeffs) - 1
    if d < 1 or d > 4:
        return (False, "степень вне диапазона 1..4")
    
    # Last coefficient non-zero
    if coeffs[-1] == 0:
        return (False, "старший коэффициент равен нулю")
    
    # All coefficients finite
    for c in coeffs:
        if not np.isfinite(c):
            return (False, "коэффициент не конечен")
    
    # Discriminant for d == 2
    if d == 2:
        c0, c1, c2 = coeffs
        disc = c1**2 - 4 * c2 * c0
        if disc < 0:
            return (False, "дискриминант отрицателен")
    
    # Check monotonicity and range
    ch = np.arange(n_ch + 1, dtype=float)
    E = poly_eval(coeffs, ch)
    
    if np.any(E <= -1000) or np.any(E >= 100000):
        return (False, "энергия вне диапазона")
    
    if np.any(np.diff(E) < 0):
        return (False, "шкала не монотонна")
    
    return (True, "")


def dev_main(C, fit, n_ch, e_lo=17.0, e_hi=3000.0):
    """
    Отклонение основной шкалы от композитной модели.
    """
    ch_grid = np.arange(n_ch, dtype=float)
    E_main = poly_eval(C, ch_grid)
    E_comp = composite_energy(C, ch_grid, fit)
    
    mask = (E_comp >= e_lo) & (E_comp <= e_hi)
    ch_used = ch_grid[mask]
    E_main_used = E_main[mask]
    E_comp_used = E_comp[mask]
    
    dev = E_main_used - E_comp_used
    
    max_dev = {}
    for r in RANGES:
        mask_r = (E_comp_used >= r[0]) & (E_comp_used < r[1])
        if np.any(mask_r):
            max_dev["%g-%g" % r] = float(np.max(np.abs(dev[mask_r])))
        else:
            max_dev["%g-%g" % r] = np.nan
    
    return max_dev


def report(C, points, n_ch, threshold_ch, shape="quad", dK=DK_STEP, e37=E37_GAMMA, de37=DE37, use37=True, verbose=True):
    """
    Полный отчёт по низкоэнергетической ветви.
    """
    # Fit main shape
    fit = fit_s(C, points, shape=shape, dK=dK, e37=e37, de37=de37, use37=use37)
    
    # Fit alternative shape
    alt_shape = "lin" if shape == "quad" else "quad"
    alt = fit_s(C, points, shape=alt_shape, dK=dK, e37=e37, de37=de37, use37=use37)
    
    # Residuals
    res = residuals(C, points, fit, e37)
    
    # LOO
    lo = loo(C, points, shape=shape, dK=dK, e37=e37, de37=de37, use37=use37)
    
    # Band
    grid = np.arange(int(threshold_ch), n_ch, dtype=float)
    bd = band(C, fit, grid, alt)
    
    # Threshold energy
    e_thr = composite_energy(C, float(threshold_ch), fit)
    
    # XML polynomial
    xml = fit_poly_xml(C, fit, n_ch, tol_fn=lambda ch: np.interp(ch, bd["ch"], bd["sigma"]))
    
    # Deviation of main scale
    dm = dev_main(C, fit, n_ch)
    
    # Check for extrapolation warning
    has_low_point = False
    for p in points:
        if p["kind"] != "intrinsic37" and p["E"] < 37.0:
            has_low_point = True
            break
    warn_extrapolation = not has_low_point
    
    if verbose:
        # Print fit info
        suffix = " (в подгонке)" if use37 else " (вне подгонки)"
        print("[LOWE] ветвь ниже %g кэВ: форма %s, s = %.4f ± %.4f (априори порядок %.2f), ступень на K-крае La %.3f, гамма-эквивалент пика 37 кэВ %.2f ± %.2f кэВ%s" % (
            E_JOIN, shape, fit["s"], fit["sigma_s"], S_PRIOR, dK, e37, de37, suffix
        ))
        
        print("[LOWE] chi2/ndf = %.2f / %d; альтернативная форма %s: s = %.4f, chi2 = %.2f" % (
            fit["chi2"], fit["ndf"], alt_shape, alt["s"], alt["chi2"]
        ))
        
        # Print residuals
        for r in res:
            loo_val = lo.get(r["name"], np.nan)
            print("[LOWE] %-14s E=%7.2f кан %8.2f: невязка до %+7.2f кэВ, после %+7.2f кэВ, LOO %+7.2f" % (
                r["name"], r["E"], r["ch"], r["before"], r["after"], loo_val
            ))
        
        print("[LOWE] порог прибора: канал %d -> %.1f кэВ (по ветви)" % (int(threshold_ch), e_thr))
        
        # Band at specific energies
        target_energies = [17, 20, 25, 30, 37, 45, 55, 75, 100]
        band_strs = []
        for e_target in target_energies:
            if e_target < e_thr:
                continue
            # Find closest channel
            idx = np.argmin(np.abs(bd["E"] - e_target))
            e_actual = bd["E"][idx]
            sigma_actual = bd["sigma"][idx]
            band_strs.append("%.0f кэВ ± %.2f" % (e_actual, sigma_actual))
        
        if band_strs:
            print("[LOWE] полоса неопределённости шкалы: " + ", ".join(band_strs))
        
        if warn_extrapolation:
            print("[LOWE] ПРЕДУПРЕЖДЕНИЕ: ниже 37 кэВ -- экстраполяция, нужен репер Am-241/Ba-133")
        
        # XML polynomial
        fallback_str = (", НЕ ЗАПИСЫВАЕТСЯ: портит шкалу выше 100 кэВ на %.2f кэВ (> 1 кэВ), в XML -- основная шкала" % xml.get("damage", float("nan"))) if xml["fallback"] else ""
        max_dev_strs = []
        for k, v in xml["max_dev"].items():
            if np.isnan(v):
                max_dev_strs.append("%s: nan" % k)
            else:
                max_dev_strs.append("%s: %.2f кэВ" % (k, v))
        print("[LOWE] полином для XML: степень %d%s; макс. отклонение от модели по диапазонам: %s" % (
            xml["deg"], fallback_str, ", ".join(max_dev_strs)
        ))
        
        # Deviation of main scale
        dm_strs = []
        for k, v in dm.items():
            if np.isnan(v):
                dm_strs.append("%s: nan" % k)
            else:
                dm_strs.append("%s: %.2f кэВ" % (k, v))
        print("[LOWE] основная шкала (без ветви) против модели: %s" % ", ".join(dm_strs))
    
    return {
        "fit": fit,
        "alt": alt,
        "residuals": res,
        "loo": lo,
        "band": bd,
        "e_threshold": float(e_thr),
        "xml": xml,
        "dev_main": dm,
        "warn_extrapolation": warn_extrapolation
    }
