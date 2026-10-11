"""selftest -- мутационный тест регистрации по форме на паре «фон + образец».
Перенесено из прототипа labr3_selfcal/test_shape_reg.py (сгенерирован gen_code.py, qwen3.8:27b), 2026-10-10.
Итог считается по НОВЫМ провалам сверх провалов исправной версии; известные провалы исправной -- в конфиге
(known_base_failures), reg:dchi2 -- INFO (согласованность с аргминимумом, не правильность)."""
import sys, json, pathlib, hashlib
import numpy as np
from scipy.optimize import least_squares
from .core import load_spectrum, rebin_conserve_var, poly_eval, N_CH, window_quality
from . import register as sr

HERE = pathlib.Path(__file__).resolve().parent
DEFAULT_CFG = HERE / "presets" / "labr3_selftest.json"


def mask_windows(cfg, names, n=N_CH):
    m = np.zeros(n, dtype=bool)   # длина -- из данных (любое число каналов)
    for n in names:
        lo, hi = cfg["windows"][n]
        m[lo:hi + 1] = True
    return m


def mask_range(lo, hi, n=N_CH):
    m = np.zeros(n, dtype=bool)
    m[lo:hi + 1] = True
    return m


def mutate_c(c, mut):
    c = np.asarray(c, float).copy()
    if mut is None:
        return c
    kind = mut.get("kind")
    if kind == "shift":
        c[0] += mut["value"]
    elif kind == "stretch":
        v = mut["value"]
        p = mut["pivot"]
        c = (1 + v) * c
        c[0] -= v * p
    return c


def win_shift(c, S, B, k, bg_var, lo, hi, scan):
    ds = np.arange(scan["lo"], scan["hi"] + 1e-9, scan["step"])
    X = np.empty(len(ds))
    for i, d in enumerate(ds):
        cd = np.asarray(c, float).copy()
        cd[0] += d
        r, v = sr.apply(cd, S, B, k, bg_var)
        X[i] = np.sum(r[lo:hi + 1] ** 2 / v[lo:hi + 1])
    i = int(np.argmin(X))
    if i == 0 or i == len(ds) - 1:
        return (float(ds[i]), float("nan"), True)
    den = X[i + 1] - 2 * X[i] + X[i - 1]
    if den <= 0:
        return (float(ds[i]), float("nan"), True)
    step = scan["step"]
    delta = ds[i] - step * (X[i + 1] - X[i - 1]) / (2 * den)
    A = den / step ** 2
    n = hi - lo + 1
    sigma = np.sqrt(2 / A) * np.sqrt(max(1.0, X[i] / n))
    return (float(delta), float(sigma), False)


def joint_fit(c, S, B, k, bg_var, mask, pivot):
    c = np.asarray(c, float)
    xe = np.arange(len(S) + 1) - 0.5

    def cmap(p):
        s = p[1] / 1000.0
        cc = (1 + s) * c
        cc[0] += p[0] - s * pivot
        return cc

    def fun(p):
        cc = cmap(p)
        e = poly_eval(cc, xe)
        if np.any(np.diff(e) <= 0):
            return np.full(int(mask.sum()), 1e3)
        sv, vs = rebin_conserve_var(S, e, len(B))
        bb, vb = sr.blur_bg(lambda x: poly_eval(cc, x), B, bg_var)
        ak = p[2] * k
        r = (sv - ak * bb) / np.sqrt(np.maximum(vs + ak ** 2 * vb, 1.0))
        return r[mask]

    best = None
    for d0 in (-1.0, -0.5, 0.0, 0.5, 1.0):
        res = least_squares(fun, [d0, 0.0, 1.0], method="trf", diff_step=1e-4, max_nfev=200)
        if best is None or res.cost < best.cost:
            best = res

    chi2 = 2 * best.cost
    ndf = mask.sum() - 3
    cov = np.linalg.pinv(best.jac.T @ best.jac) * max(1.0, chi2 / ndf)
    err = np.sqrt(np.diag(cov))
    return {
        "delta": float(best.x[0]),
        "s": float(best.x[1] / 1000),
        "a": float(best.x[2]),
        "sd": float(err[0]),
        "ss": float(err[1] / 1000),
        "sa": float(err[2]),
        "chi2ndf": float(chi2 / ndf),
    }


def check_set(c, S, B, k, bg_var, cfg, names, prefix):
    # при sigma_method=bootstrap -- глобальный минимум на сетке + бутстрап-сигма (robust_shift.py)
    if cfg.get("sigma_method") == "bootstrap":
        from . import selftest_boot as robust_shift
        return robust_shift.check_set_boot(c, S, B, k, bg_var, cfg, names, prefix)
    th = cfg["thresholds"]
    r, v = sr.apply(c, S, B, k, bg_var)
    ordered = sorted(names, key=lambda n: cfg["windows"][n][0])
    checks = []
    shifts = {}
    for n in ordered:
        lo, hi = cfg["windows"][n]
        q = window_quality(r, v, lo, hi)
        checks.append((prefix + n + ":z_dip", float(q["z_dip"]), float(th["z_dip"]), abs(q["z_dip"]) <= th["z_dip"]))
        checks.append((prefix + n + ":z_chi2", float(q["z_chi2"]), float(th["z_chi2"]), q["z_chi2"] <= th["z_chi2"]))
        d, sg, edge = win_shift(c, S, B, k, bg_var, lo, hi, cfg["shift_scan"])
        shifts[n] = (d, sg, edge)
        if edge or not np.isfinite(sg):
            z = float("inf")
        else:
            z = d / sg
        checks.append((prefix + n + ":shift_z", float(z), float(th["z_shift_win"]), abs(z) <= th["z_shift_win"]))

    if len(ordered) >= 2:
        first = ordered[0]
        last = ordered[-1]
        d1, sg1, e1 = shifts[first]
        d2, sg2, e2 = shifts[last]
        if e1 or e2 or not np.isfinite(sg1) or not np.isfinite(sg2):
            z = float("inf")
        else:
            z = (d2 - d1) / np.hypot(sg2, sg1)
        checks.append((prefix + "dist:" + first + "->" + last, float(z), float(th["z_dist"]), abs(z) <= th["z_dist"]))

    union_mask = np.zeros(len(B), dtype=bool)
    for n in ordered:
        lo, hi = cfg["windows"][n]
        union_mask[lo:hi + 1] = True
    j = joint_fit(c, S, B, k, bg_var, union_mask, cfg["pivot_joint"])
    checks.append((prefix + "joint:shift_z", float(j["delta"] / j["sd"]), float(th["z_joint_shift"]), abs(j["delta"] / j["sd"]) <= th["z_joint_shift"]))
    checks.append((prefix + "joint:scale_z", float(j["s"] / j["ss"]), float(th["z_joint_scale"]), abs(j["s"] / j["ss"]) <= th["z_joint_scale"]))
    checks.append((prefix + "joint:amp_z", float((j["a"] - 1) / j["sa"]), float(th["z_joint_amp"]), abs((j["a"] - 1) / j["sa"]) <= th["z_joint_amp"]))
    return checks


def integral_checks(c, S, B, k, bg_var, cfg, prefix):
    th = cfg["thresholds"]
    r, v = sr.apply(c, S, B, k, bg_var)
    checks = []
    for lo, hi in cfg["integrals"]:
        z = r[lo:hi + 1].sum() / np.sqrt(v[lo:hi + 1].sum())
        checks.append((prefix + "int:%d-%d" % (lo, hi), float(z), float(th["z_int"]), abs(z) <= th["z_int"]))
    return checks


def register_all(S, B, k, bg_var, cfg):
    anchor = sr.anchor37(B, S, verbose=False, **cfg.get("anchor", {}))
    masks = {
        "full": mask_range(*cfg["full_mask"], n=len(B)),
        "A": mask_windows(cfg, cfg["half_A"], len(B)),
        "B": mask_windows(cfg, cfg["half_B"], len(B)),
    }
    if "info_holdout" in cfg:   # информационная регистрация вне зачёта
        masks["info"] = mask_windows(cfg, cfg["info_holdout"]["fit"], len(B))
    regs = {}
    for key in masks:
        regs[key] = sr.register(S, B, k, masks[key], anchor, bg_var=bg_var, verbose=False, **cfg["registration"])   # n_ms и окна регистрации (автопресет)
        c = regs[key]["c"]
        hv = np.round(sr.hvals(c, cfg["registration"].get("test_ch")), 3)
        print("[REG] fit=%s c=[%.6f, %.7f, %.4e] chi2/ndf=%.4f h(%s)=%s" % (
            key, c[0], c[1], c[2], regs[key]["chi2ndf"], ",".join("%g" % v for v in cfg["registration"].get("test_ch", sr.TEST_CH)), hv))
    return regs


def run_checks(regs, S, B, k, bg_var, cfg, mut=None):
    checks = []
    c_full = mutate_c(regs["full"]["c"], mut)
    checks.extend(check_set(c_full, S, B, k, bg_var, cfg, list(cfg["windows"].keys()), "full|"))
    checks.extend(integral_checks(c_full, S, B, k, bg_var, cfg, "full|"))

    c_a = mutate_c(regs["A"]["c"], mut)
    checks.extend(check_set(c_a, S, B, k, bg_var, cfg, cfg["half_B"], "fitA->B|"))

    c_b = mutate_c(regs["B"]["c"], mut)
    checks.extend(check_set(c_b, S, B, k, bg_var, cfg, cfg["half_A"], "fitB->A|"))
    # dchi2 функционала регистрации при проверяемом h против h регистрации (full) -- ТОЛЬКО INFO, вне зачёта:
    # сравнивает h с аргминимумом той же регистрации (у исправной всегда +0.000), т.е. проверяет
    # согласованность с аргминимумом, а не правильность (вторая проверка 2026-10-10, замечание 2).
    if "dchi2_reg" in cfg["thresholds"]:
        dc = regs["full"]["prob"].chi2(c_full) - regs["full"]["chi2"]
        thr = cfg["thresholds"]["dchi2_reg"]
        print("[INFO вне зачёта: согласованность с аргминимумом, не правильность] full|reg:dchi2 %+10.3f thr %.1f %s"
              % (dc, thr, "PASS" if dc <= thr else "FAIL"))
    if "info" in regs:   # вне зачёта: ограничение метода (без сумм-пика и β-края)
        info = check_set(mutate_c(regs["info"]["c"], mut), S, B, k, bg_var, cfg,
                         cfg["info_holdout"]["check"], "INFO fitBold->Aold|")
        for ch in info:
            print("[INFO вне зачёта, ограничение метода] %-45s %+8.3f thr %.1f %s" % (ch[0], ch[1], ch[2], "PASS" if ch[3] else "FAIL"))
    return checks


def load_alt_bg(fname, B, live_b, cfg=None):
    """Другой фон, зарегистрированный на шкалу эталонного (мутация «чужой фон»)."""
    alt = load_spectrum(fname, n_ch=len(B))
    A = alt["counts"].copy()
    A[-1] = 0.0
    print("[ALTBG] %s live=%.1f identical_to_reference=%s" % (fname, alt["live"], np.array_equal(A, B)))
    if np.array_equal(A, B):
        raise RuntimeError("background identical; mutation is vacuous")
    kk = alt["live"] / live_b
    cfg = cfg or {}   # якорь, маска и окна регистрации -- из конфига (автопресет)
    a = cfg.get("anchor", {})
    anchor = sr.anchor37(B, A, verbose=False, **({"x_bg": a["x_bg"], "x_sm": a["x_sm"], "w0": a["w0"], "hw": a["hw"]} if a else {}))
    rk = dict(cfg.get("registration", {}), n_ms=10)
    r = sr.register(A, B, kk, mask_range(*cfg.get("full_mask", (430, 1700)), n=len(B)), anchor, verbose=False, **rk)
    e = poly_eval(r["c"], np.arange(len(A) + 1) - 0.5)
    B2, var2 = rebin_conserve_var(A, e, len(B))
    print("[ALTBG] registered onto reference grid: c=%s chi2/ndf=%.4f" % (r["c"], r["chi2ndf"]))
    return (B2, var2, alt["live"])


def summarize(label, checks):
    failed = [c for c in checks if not c[3]]
    print("[%s] checks=%d failed=%d" % (label, len(checks), len(failed)))
    for c in failed:
        print("  FAIL %s = %+.3f (thr %.1f)" % (c[0], c[1], c[2]))
    return len(failed)


def main(bg_path, sample_path, cfg_path=None, alt_bg=(), auto=True, manual=None):
    """Код возврата 0: у исправной нет провалов сверх known_base_failures и каждая мутация даёт >=1 НОВЫЙ провал."""
    cfg_path = pathlib.Path(cfg_path or DEFAULT_CFG)
    raw = cfg_path.read_bytes()
    cfg = json.loads(raw.decode("utf-8"))
    print("[CFG] %s sha256=%s" % (cfg_path.name, hashlib.sha256(raw).hexdigest()))

    S, B, k, live_s, live_b = sr.load_inputs(bg_path, sample_path)
    if auto:   # окна конфига -- в каналах этого прибора (автопресет по собственным структурам фона)
        from . import autopreset as ap
        cfg = ap.selftest_cfg_for(ap.autopreset(B, S, load_spectrum(bg_path)["cal"], manual=manual), cfg)

    regs0 = register_all(S, B, k, None, cfg)
    base = run_checks(regs0, S, B, k, None, cfg)
    for c in base:
        print("[BASE] %-40s %+8.3f thr %.1f %s" % (c[0], c[1], c[2], "PASS" if c[3] else "FAIL"))
    nb = summarize("BASE", base)

    mut_results = []
    alt = [{"name": "bg=" + pathlib.Path(f).stem, "kind": "bg", "file": str(f)} for f in alt_bg]
    for mut in list(cfg["mutations"]) + alt:
        kind = mut["kind"]
        if kind in ("shift", "stretch"):
            checks = run_checks(regs0, S, B, k, None, cfg, mut)
        elif kind == "k":
            k2 = k * (1 + mut["value"])
            regs = register_all(S, B, k2, None, cfg)
            checks = run_checks(regs, S, B, k2, None, cfg)
        elif kind == "bg":
            f = pathlib.Path(mut["file"])
            if not f.is_file():
                print("[SKIP] мутация %s: файл фона не найден (%s)" % (mut["name"], mut["file"]))
                continue
            B2, var2, live2 = load_alt_bg(f, B, live_b, cfg)
            k2 = live_s / live2
            regs = register_all(S, B2, k2, var2, cfg)
            checks = run_checks(regs, S, B2, k2, var2, cfg)
        else:
            checks = run_checks(regs0, S, B, k, None, cfg, mut)
        nf = summarize("MUT " + mut["name"], checks)
        differs = [a[3] for a in checks] != [b[3] for b in base]
        base_failed = {c[0] for c in base if not c[3]}
        new = [c[0] for c in checks if not c[3] and c[0] not in base_failed]
        n_inf = sum(1 for c in checks if not c[3] and c[0] not in base_failed and not np.isfinite(c[1]))
        print("[MUT %s] новых провалов сверх базовых: %d (по конечному z: %d; только «минимум на краю скана» = inf: %d) %s"
              % (mut["name"], len(new), len(new) - n_inf, n_inf, new))
        mut_results.append((mut["name"], nf, len(checks), differs, len(new), n_inf))

    print("[SUMMARY] %-22s failed %3d / %3d" % ("BASE", nb, len(base)))
    for name, nf, ntot, differs, nnew, ninf in mut_results:
        print("[SUMMARY] %-22s failed %3d / %3d  новых %3d  %s" % (name, nf, ntot, nnew, "ПОЙМАНА" if nnew > 0 else "НЕ ПОЙМАНА"))

    caught = sum(1 for r in mut_results if r[3])
    print("[SA-3a] mutants whose verdict vector differs from the correct pipeline: %d of %d" % (caught, len(mut_results)))

    print("[ТРЕБУЕТ ТОЛКОВАНИЯ]")
    interp = []
    for name, nf, ntot, differs, nnew, ninf in mut_results:
        if nnew == 0:
            interp.append("  не поймано (новых провалов нет): %s" % name)
        elif nnew == ninf:
            interp.append("  поймано ТОЛЬКО «минимумом на краю скана» (inf), конечным z не ловится: %s" % name)
    for c in base:
        if abs(c[1]) > 0.8 * c[2]:
            interp.append("  близко к порогу: %s = %+.3f" % (c[0], c[1]))
    if not interp:
        print("  (нет пунктов)")
    else:
        for line in interp:
            print(line)

    known = set(cfg.get("known_base_failures", []))
    base_failed = sorted(c[0] for c in base if not c[3])
    unexpected = [n for n in base_failed if n not in known]
    caught_new = sum(1 for r in mut_results if r[4] > 0)
    print("ИТОГ (по НОВЫМ провалам): исправная -- провалов %d, из них известных (задокументированных) %d %s, "
          "неожиданных %d %s; мутаций с >=1 новым провалом %d из %d" % (
              len(base_failed), len(base_failed) - len(unexpected), [n for n in base_failed if n in known],
              len(unexpected), unexpected, caught_new, len(mut_results)))
    ok = not unexpected and caught_new == len(mut_results)
    print("ИТОГ: %s" % ("ЗЕЛЁНЫЙ" if ok else "КРАСНЫЙ"))
    return 0 if ok else 1
