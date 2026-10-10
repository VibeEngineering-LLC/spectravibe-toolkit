"""
Модуль pipeline.py: оркестрация калибровки и вычитания фона.
Объединяет автоматический пресет, робастную регистрацию,
энергетическую шкалу (пара реперных точек + низкочастотная ветвь)
и проверку собственных структур детектора.
Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-10-10
"""

import json
import pathlib
import numpy as np
from .core import load_spectrum, poly_eval, poly_deriv, poly_inverse
from .peakfit import fit_peak, fwhm_fn_from_points
from . import register as reg
from . import autopreset as ap
from . import robust as rb
from . import pairref as pr
from . import lowscale as lw
from . import ops
from . import verdict as vd

STRUCT_NAMES = ("La_beta_edge", "La_compton", "La_sum_K40", "Tl208_alpha1", "alpha2", "alpha3")
E_SUM = 1471.0
E_K40 = 1460.82


def _load_json(path):
    return json.loads(pathlib.Path(path).read_text("utf-8"))


def nested_pair(bg_path, sample_path, S=None, B=None):
    """Правка Claude: «образец» -- продолжение того же набора, что и «фон» (одно время начала, счёт образца >= фона
    во всех каналах). Тогда разность -- не «образец минус фон», а chi2/ndf регистрации ~ 0.2 (счёт коррелирован)."""
    b, s = load_spectrum(bg_path), load_spectrum(sample_path)
    B = b["counts"] if B is None else B
    S = s["counts"] if S is None else S
    if b["start"] and b["start"] == s["start"] and len(B) == len(S) and np.all(S[:-1] >= B[:-1]):
        return ("образец -- продолжение того же измерения, что и фон (время начала %s у обоих, счёт образца не меньше "
                "фона во всех каналах): разность не является «образец минус фон»; нужен независимый фон" % b["start"])
    return ""


def setup(bg_path, sample_path, auto=True, manual=None, verbose=True):
    S, B, k, live_s, live_b = reg.load_inputs(bg_path, sample_path)
    cal = list(load_spectrum(bg_path)["cal"])
    nested = nested_pair(bg_path, sample_path, S, B)
    if nested and verbose:
        print("[PAIR] ПРЕДУПРЕЖДЕНИЕ: " + nested)
    
    a = ap.autopreset(B, S, cal, manual=manual, verbose=verbose)
    
    if not auto:
        ref = ap.load_ref()
        a["register_kw"] = dict(ref["register"])
        a["sum_window"] = list(ref["sum_window"])
        a["identity"] = True
        if verbose:
            print("[AUTO] автопресет отключён (--no-auto): пресеты в каналах эталона")

    A = a["anchors"]
    fw = fwhm_fn_from_points([A["a37"]["x"], A["sum"]["x"]], [A["a37"]["fwhm"], A["sum"]["fwhm"]])
    
    try:
        scale = pr.prov_scale(cal, A["sum"]["x"], E_SUM)
    except Exception:
        scale = [0.0, E_SUM / A["sum"]["x"], 0.0]
        if verbose:
            print("[WARN] не удалось построить шкалу по эталону, использована линейная аппроксимация")

    k40 = A["sum"]["x"] - (E_SUM - E_K40) / poly_deriv(scale, A["sum"]["x"])
    
    # правка Claude (раунд 4, п.2): окна структур -- в энергиях (ядерные данные), каналы -- через шкалу пары
    windows = vd.windows_to_channels(vd.load_windows_kev(), scale, len(B))

    return {
        "S": S, "B": B, "k": k, "live_s": live_s, "live_b": live_b,
        "cal": cal, "auto": a, "fwhm_fn": fw, "scale": scale,
        "k40_ch": float(k40), "windows": windows, "n": len(B),
        "bg_path": str(bg_path), "sample_path": str(sample_path), "nested": nested
    }


def verdict_kw(ctx, excluded, steps_ch):
    """Правка Claude (раунд 4): общие аргументы вердикта по окнам -- исключения, ступени, K-40 (канал, сигма)."""
    k40 = {"ch": float(ctx["k40_ch"]), "sigma": float(ctx["fwhm_fn"](ctx["k40_ch"]) / 2.3548)}
    return {"excl": [[int(a), int(b)] for a, b in excluded], "fwhm_fn": ctx["fwhm_fn"],
            "steps_ch": [float(s) for s in steps_ch], "k40": k40}


def compton_fn_for(ctx):
    scale = ctx["scale"]
    n = ctx["n"]
    
    def _compton(ch):
        E_peak = poly_eval(scale, ch)
        E_edge = pr.compton_edge_kev(E_peak)
        if E_edge < 50.0:
            return None
        try:
            return poly_inverse(scale, E_edge, hi=n + 200)
        except Exception:
            return None
    return _compton


def extra_fn_for(ctx, verbose=True):
    calls = []
    def _extra(peak_chs):
        calls.append(1)
        if len(calls) == 1:   # правка Claude: на итерации 0 пики -- ещё и диполи рассогласования; отождествлять с итерации 1
            return []
        chs, nuclides = pr.chain_centres(peak_chs, ctx["scale"], ctx["fwhm_fn"], ctx["n"])
        ctx["nuclides"] = nuclides
        if nuclides and verbose:
            print("[ROBUST] отождествлены по найденным пикам: %s; исключаются окна всех их линий >= 1 %%" % ", ".join(nuclides))
        return chs
    return _extra


def run_register(bg_path, sample_path, auto=True, manual=None, robust="auto", overrides=None, n_ms=20, profiles=False, out_json=None, z_peak=4.5, verbose=True):
    ctx = setup(bg_path, sample_path, auto, manual, verbose)
    kw = dict(ctx["auto"]["register_kw"])
    if overrides:
        for key, v in overrides.items():
            if v is not None:
                kw[key] = v

    if verbose:
        print("[IN] k=%.6f live_s=%.3f live_b=%.3f каналов=%d" % (ctx["k"], ctx["live_s"], ctx["live_b"], ctx["n"]))

    anchor = reg.anchor37(ctx["B"], ctx["S"], x_bg=kw["anchor_bg"], x_sm=kw["anchor_sample"], w0=kw["anchor_fwhm"], hw=int(kw["anchor_hw"]), verbose=verbose)
    
    mask = np.zeros(ctx["n"], dtype=bool)
    lo, hi = kw["mask"]
    mask[lo:hi+1] = True

    rkw = {
        "n_ms": n_ms,
        "coarse_window": tuple(kw["coarse_window"]),
        "c2_lo": kw["c2_lo"],
        "c2_hi": kw["c2_hi"],
        "c2_step": kw["c2_step"],
        "test_ch": kw["test_ch"]
    }

    rr = rb.register_robust(
        ctx["S"], ctx["B"], ctx["k"], mask, anchor, ctx["fwhm_fn"],
        k40_ch=ctx["k40_ch"], mode=robust, verbose=verbose,
        compton_fn=compton_fn_for(ctx), extra_fn=extra_fn_for(ctx, verbose),
        z_peak=z_peak, protect_fwhm=3.0, **rkw
    )   # правка Claude: окно сумм-пика La-138 (+-3 ПШПВ от K-40) не исключается -- K-40 моделируется

    r, out = reg.report(ctx["S"], ctx["B"], ctx["k"], anchor, rr["mask"], out_json=None, profiles=profiles, r=rr["reg"], **rkw)
    
    net, var = reg.apply(rr["c"], ctx["S"], ctx["B"], ctx["k"])
    
    # правка Claude (раунд 4, п.1): вердикт -- на остатке БЕЗ модели образца; простая и устойчивая -- на одной основе
    # (те же окна, исключения, ступени, K-40 только в окне сумм-пика)
    vkw = verdict_kw(ctx, rr["excluded"], rr.get("steps_ch", []))
    checks = vd.structure_verdict(net, var, ctx["windows"], **vkw)
    c_plain = (rr.get("plain") or {}).get("c", rr["c"])
    net_p, var_p = reg.apply(c_plain, ctx["S"], ctx["B"], ctx["k"])
    checks_plain = vd.structure_verdict(net_p, var_p, ctx["windows"], **vkw)

    scale = ctx["scale"]
    excluded_kev = [[int(lo), int(hi), float(poly_eval(scale, lo)), float(poly_eval(scale, hi))] for lo, hi in rr["excluded"]]

    if verbose:
        print("[ROBUST] режим %s: chi2/ndf простой %.4f -> итог %.4f; исключено окон %d" % (rr["mode"], rr["chi2ndf_plain"], rr["chi2ndf"], len(rr["excluded"])))
        for e in excluded_kev:
            print("[ROBUST]   кан %d-%d (%.0f-%.0f кэВ)" % (e[0], e[1], e[2], e[3]))
        if rr["k40"]:
            print("[ROBUST] K-40 образца под сумм-пиком: площадь %.0f ± %.0f отсчётов (кан %.1f)" % (rr["k40"]["area"], rr["k40"]["area_err"], rr["k40"]["ch"]))
        if rr["mode"] == "robust" and not rr.get("converged", True):
            print("[ROBUST] ПРЕДУПРЕЖДЕНИЕ: устойчивая регистрация не сошлась -- итог не принимается (all_ok=False)")
        print(vd.format_table(checks, "[STRUCT] окна собственных структур, регистрация %s (остаток без модели образца, "
                              "базовая линия по боковым полосам):" % rr["mode"]))
        if rr["mode"] == "robust":
            print(vd.format_table(checks_plain, "[STRUCT] те же окна, простая регистрация (та же основа):"))

    a = ctx["auto"]
    A = a["anchors"]
    
    out["auto"] = {
        "map": [float(v) for v in a["map"]],
        "identity": bool(a["identity"]),
        "kind": a.get("kind", "unknown"),
        "r": float(a["r"]),
        "anchors": {
            nm: {
                "x": float(A[nm]["x"]),
                "err": float(A[nm]["err"]),
                "ok": bool(A[nm]["ok"]),
                "manual": bool(A[nm].get("manual", False))
            } for nm in A
        },
        "sample_a37": float(a["sample_a37"]["x"]),
        "warnings": list(a["warnings"]),
        "sum_window": [int(v) for v in a["sum_window"]],
        "fwhm_points": [[float(v) for v in p] for p in a["fwhm_points"]],
        "register_kw": {k: (list(v) if isinstance(v, (list, tuple)) else v) for k, v in kw.items()}
    }
    
    out["robust"] = {
        "mode": rr["mode"],
        "chi2ndf_plain": float(rr["chi2ndf_plain"]),
        "chi2ndf": float(rr["chi2ndf"]),
        "z_plain": float(rr["z_plain"]),
        "excluded": excluded_kev,
        "k40": rr["k40"],
        "nuclides": ctx.get("nuclides", []),
        "converged": bool(rr.get("converged", True)),
        "iterations": int(rr.get("iterations", 0)),
        "steps_ch": [float(v) for v in rr.get("steps_ch", [])],
        "excluded_ch": [[int(a), int(b)] for a, b in rr["excluded"]],
        "same_basis": rr.get("same_basis"),
        "model": (None if rr["model"] is None else [round(float(v), 3) for v in rr["model"]])
    }

    out["structure"] = checks
    out["structure_plain"] = checks_plain
    out["structure_summary"] = vd.summarize(checks)
    out["structure_plain_summary"] = vd.summarize(checks_plain)
    out["k40_ch"] = float(ctx["k40_ch"])
    out["scale_prov"] = [float(v) for v in ctx["scale"]]

    if out_json:
        p = pathlib.Path(out_json)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
        if verbose:
            print("[OUT] " + str(p))

    return {
        "ctx": ctx,
        "rr": rr,
        "out": out,
        "net": net,
        "var": var,
        "checks": checks,
        "anchor": anchor,
        "h": [float(v) for v in rr["c"]]
    }


def _h_from(h):
    if isinstance(h, dict) and "c" in h:
        return list(h["c"])
    if isinstance(h, (str, pathlib.Path)):
        p = pathlib.Path(h)
        if p.exists():
            return list(_load_json(p)["c"])
    return list(h)


def low_points(ctx, cal_res, h, net, var, lu, lowe, C, verbose=True):
    points = []
    A = ctx["auto"]["anchors"]
    n = ctx["n"]
    
    points.append({
        "name": "La138_37_intr",
        "E": lowe["e37_gamma"],
        "dE": lowe["de37"],
        "ch": A["a37"]["x"],
        "ch_err": max(A["a37"]["err"], 0.01),
        "kind": "intrinsic37"
    })
    
    for p in cal_res["points"]:
        if p["E"] < 100.0:
            points.append({
                "name": p["name"],
                "E": p["E"],
                "dE": p["E_err"],
                "ch": p["ch"],
                "ch_err": p["ch_err"],
                "kind": "photon"
            })
    
    def _add_point(pt):
        for existing in points:
            if abs(pt["E"] - existing["E"]) < 1.0:
                return False
        points.append(pt)
        return True

    for item in lowe.get("bg_xray", []):
        E, dE = pr.xray_blend(item["element"], item["lines"])
        if item.get("mix_element"):   # правка Claude: в фоне K-рентген Pb экрана смешан с Bi K от Pb-214 (доля неизвестна)
            E2, dE2 = pr.xray_blend(item["mix_element"], item["lines"])
            E, dE = (E + E2) / 2.0, max(dE, dE2, abs(E2 - E) / 2.0)
        try:
            x0 = poly_inverse(C, E, hi=n + 200)
        except Exception:
            continue
        fw = ctx["fwhm_fn"](x0)
        r = fit_peak(ctx["B"], x0, fw, half_window=max(5, int(round(1.6 * fw))))
        if r and r["signif"] >= lowe["min_signif"] and abs(r["centroid"] - x0) <= 1.5 * fw:
            _add_point({
                "name": item.get("name", "bg_%s" % item["element"]),
                "E": E,
                "dE": dE,
                "ch": r["centroid"],
                "ch_err": r["centroid_err"],
                "kind": "photon"
            })
            if verbose:
                print("[LOWE] репер %s %.2f кэВ: кан %.2f ± %.2f, значимость %.1f" % (item["element"], E, r["centroid"], r["centroid_err"], r["signif"]))
        else:
            if verbose:
                sig = r["signif"] if r else 0.0
                print("[LOWE] репер %s %.2f кэВ: не найден (значимость %.1f / центроид далеко)" % (item["element"], E, sig))

    if lu:
        for item in lowe.get("lu_xray", []):
            E, dE = pr.xray_blend(item["element"], item["lines"])
            try:
                x0 = poly_inverse(C, E, hi=n + 200)
            except Exception:
                continue
            fw = ctx["fwhm_fn"](x0)
            r = pr.fit_net_peak(net, var, x0, fw, max(5, int(round(1.6 * fw))))
            if r and r["signif"] >= lowe["min_signif"] and abs(r["centroid"] - x0) <= 1.5 * fw:
                _add_point({
                    "name": item.get("name", "lu_%s" % item["element"]),
                    "E": E,
                    "dE": dE,
                    "ch": r["centroid"],
                    "ch_err": r["centroid_err"],
                    "kind": "photon"
                })
        
        for item in lowe.get("lu_gamma", []):
            E = item["E"]
            dE = item.get("dE", 0.05)
            try:
                x0 = poly_inverse(C, E, hi=n + 200)
            except Exception:
                continue
            fw = ctx["fwhm_fn"](x0)
            r = pr.fit_net_peak(net, var, x0, fw, max(5, int(round(1.6 * fw))))
            if r and r["signif"] >= lowe["min_signif"] and abs(r["centroid"] - x0) <= 1.5 * fw:
                _add_point({
                    "name": item.get("name", "lu_gamma_%.1f" % E),
                    "E": E,
                    "dE": dE,
                    "ch": r["centroid"],
                    "ch_err": r["centroid_err"],
                    "kind": "photon"
                })

    return points


def source_points(src_path, ctx, lowe, nuclide, C, verbose=True):
    rs = run_register(ctx["bg_path"], src_path, auto=True, robust="on", n_ms=10, verbose=False)
    hs = rs["h"]
    src = reg.load_inputs(ctx["bg_path"], src_path)[0]
    n = ctx["n"]
    
    lines = []
    src_cfg = lowe["sources"][nuclide]
    
    for line in pr.library_lines([src_cfg["gamma"]], lowe["source_min_int"], 15.0, lowe["source_e_max"]):
        lines.append({
            "name": "%s_%.1f" % (nuclide, line["E"]),
            "E": line["E"],
            "dE": 0.05,
            "weight": line["I"]
        })
    
    for item in src_cfg.get("xray", []):
        E, dE = pr.xray_blend(item["element"], item["lines"])
        lines.append({
            "name": "%s_xray_%s" % (nuclide, item["element"]),
            "E": E,
            "dE": dE,
            "weight": 1.0
        })
    
    merged = []
    for line in lines:
        x = poly_inverse(C, line["E"], hi=n + 200)
        fwhm_kev = ctx["fwhm_fn"](x) * poly_deriv(C, x)
        
        found = False
        for m in merged:
            if abs(line["E"] - m["E"]) < 1.0 * fwhm_kev:
                w_total = m["weight"] + line["weight"]
                m["E"] = (m["E"] * m["weight"] + line["E"] * line["weight"]) / w_total
                m["weight"] = w_total
                m["dE"] = max(m["dE"], line["dE"], abs(line["E"] - m["E"]) / 2.0)
                found = True
                break
        if not found:
            merged.append(dict(line))
    
    points = []
    grid = np.arange(n, dtype=float)
    eval_hs = poly_eval(hs, grid)
    
    for line in merged:
        E = line["E"]
        dE = line["dE"]
        try:
            x_bg = poly_inverse(C, E, hi=n + 200)
        except Exception:
            continue
        
        x_src = float(np.interp(x_bg, eval_hs, grid))
        fw_src = ctx["fwhm_fn"](x_bg) / max(poly_deriv(hs, x_src), 1e-6)
        
        r = fit_peak(src, x_src, fw_src, half_window=max(5, int(round(1.6 * fw_src))))
        
        if r and r["signif"] >= lowe["min_signif"] and abs(r["centroid"] - x_src) <= 1.5 * fw_src:
            ch_bg = poly_eval(hs, r["centroid"])
            ch_err = float(np.hypot(poly_deriv(hs, r["centroid"]) * r["centroid_err"], 0.15))
            points.append({
                "name": line["name"],
                "E": E,
                "dE": dE,
                "ch": ch_bg,
                "ch_err": ch_err,
                "kind": "photon"
            })
            if verbose:
                print("[SRC] %s %.2f кэВ: кан источника %.2f -> кан фона %.2f, значимость %.1f" % (line["name"], E, r["centroid"], ch_bg, r["signif"]))
        else:
            if verbose:
                print("[SRC] %s %.2f кэВ: не найден" % (line["name"], E))
    
    if verbose:
        print("[SRC] реперов внешнего источника %s: %d" % (nuclide, len(points)))
    
    return points


def run_ecal(bg_path, sample_path, h, preset=None, auto=True, manual=None, low=True, use37=True, e37=None, source=None, source_nuclide=None, robust_json=None, out_json=None, verbose=True):
    ctx = setup(bg_path, sample_path, auto, manual, verbose=False)
    hh = _h_from(h)
    n = ctx["n"]
    
    preset_cfg = _load_json(preset or (ap.HERE / "presets" / "labr3_ecal.json"))
    cfg0 = ap.ecal_cfg_for(ctx["auto"], preset_cfg) if auto else preset_cfg
    
    net, var = reg.apply(hh, ctx["S"], ctx["B"], ctx["k"])
    
    # правка Claude: два прохода -- (1) по предварительной шкале с допуском 3 %, (2) по шкале прохода 1 без допуска
    cfg, lines = pr.pair_ecal_cfg(cfg0, ctx["B"], net, var, hh, ctx["scale"], ctx["fwhm_fn"], verbose=False,
                                  tol_rel=0.03, min_sep=1.5)
    if lines and len(cfg["points"]) >= 4:
        sc2 = list(ops.calibrate(cfg, ctx["B"], ctx["S"], hh)["coeffs"]) if verbose is not None else None
        if verbose:
            print("[PAIR] проход 2: шкала прохода 1 %s" % ["%.6g" % v for v in sc2])
        cfg, lines = pr.pair_ecal_cfg(cfg0, ctx["B"], net, var, hh, sc2, ctx["fwhm_fn"], verbose=verbose, tol_rel=0.0)
    elif verbose:
        cfg, lines = pr.pair_ecal_cfg(cfg0, ctx["B"], net, var, hh, ctx["scale"], ctx["fwhm_fn"], verbose=True,
                                      tol_rel=0.03, min_sep=1.5)
    lu = len(lines) == 0 and pr.is_lu_sample(net, var, cfg0, hh, ctx["fwhm_fn"])
    
    res = ops.calibrate(cfg, ctx["B"], ctx["S"], hh)
    
    C = list(res["coeffs"])
    if verbose and np.isfinite(res["chi2ndf"]) and res["chi2ndf"] > 3.0:   # правка Claude
        print("[CAL] ПРЕДУПРЕЖДЕНИЕ: chi2/ndf = %.1f -- реперы не описываются полиномом степени %d в пределах погрешностей "
              "(макс. |невязка| %.1f кэВ); собственные структуры (сумм-пик, β-край) имеют свой гамма-эквивалент"
              % (res["chi2ndf"], len(C) - 1, max(abs(p["residual"]) for p in res["points"])))
    result = dict(res)
    result["lu"] = bool(lu)
    result["lines"] = lines
    result["ecal_points_cfg"] = cfg["points"]
    
    fw_lines = list(lines)   # правка Claude: + пики фона (кроме аннигиляции -- доплеровское уширение) как точки ПШПВ
    for p in ([] if lu else cfg["points"]):
        if p.get("spectrum") == "bg" and p.get("kind") == "peak" and not p["name"].startswith(("Ann", "La138")):
            r = fit_peak(ctx["B"], p["x0"], p["w0"], half_window=int(p["hw"]))
            if r is not None and r["signif"] >= 10:
                fw_lines.append({"ch": r["centroid"], "fwhm": r["fwhm"], "fwhm_err": r["fwhm_err"], "signif": r["signif"]})
    fw_pts = pr.pair_fwhm_points(ctx["auto"]["anchors"], fw_lines, ap.fwhm_points_for(ctx["auto"], ap.load_ref()["fwhm_points"]), lu, verbose=verbose)
    result["fwhm_points"] = [list(p) for p in fw_pts]
    
    lr = None
    if low:
        lowe = _load_json(ap.HERE / "presets" / "labr3_lowE.json")
        if e37 is not None:
            lowe["e37_gamma"] = float(e37)
        
        pts = low_points(ctx, res, hh, net, var, lu, lowe, C, verbose)
        if source:
            pts += source_points(source, ctx, lowe, source_nuclide or "Am-241", C, verbose)
        
        thr = 1
        for i in range(n):
            if ctx["B"][i] > 0:
                thr = i + 1
                break
        
        try:
            lr = lw.report(C, pts, n, thr, e37=lowe["e37_gamma"], de37=lowe["de37"], use37=use37, verbose=verbose)
        except ValueError as exc:
            if verbose:
                print("[LOWE] ветвь не построена: %s" % exc)
            lr = None
        
        if lr:
            band_ch = lr["band"]["ch"]
            band_E = lr["band"]["E"]
            band_sigma = lr["band"]["sigma"]
            
            result["low"] = {
                "s": lr["fit"]["s"],
                "sigma_s": lr["fit"]["sigma_s"],
                "shape": lr["fit"]["shape"],
                "dK": lr["fit"]["dK"],
                "e37": lowe["e37_gamma"],
                "points": pts,
                "residuals": lr["residuals"],
                "loo": lr["loo"],
                "e_threshold": lr["e_threshold"],
                "xml": lr["xml"],
                "dev_main": lr["dev_main"],
                "warn_extrapolation": lr["warn_extrapolation"],
                "band": {
                    "ch": list(band_ch[::4]),
                    "E": list(band_E[::4]),
                    "sigma": list(band_sigma[::4])
                }
            }
            result["coeffs_xml"] = lr["xml"]["coeffs"]
    
    if not low or lr is None:
        result["coeffs_xml"] = C
    
    def _json_safe(obj):
        if isinstance(obj, dict):
            return {k: _json_safe(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [_json_safe(v) for v in obj]
        elif isinstance(obj, (np.integer,)):
            return int(obj)
        elif isinstance(obj, (np.floating,)):
            return float(obj)
        elif isinstance(obj, (np.bool_,)):
            return bool(obj)
        elif isinstance(obj, np.ndarray):
            return _json_safe(obj.tolist())
        else:
            return obj
    
    result = _json_safe(result)
    
    if out_json:
        p = pathlib.Path(out_json)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(result, ensure_ascii=False, indent=1), encoding="utf-8")
        if verbose:
            print("[OUT] " + str(p))
    
    return result


def _ecal_from(ecal):
    if isinstance(ecal, dict):
        return ecal
    if isinstance(ecal, (str, pathlib.Path)):
        return _load_json(ecal)
    if isinstance(ecal, list):
        return {"coeffs": list(ecal), "coeffs_xml": list(ecal)}
    raise ValueError("Неподдерживаемый тип ecal")


def run_subtract(bg_path, sample_path, h, ecal, fw_points, out_zeroed, out_signed, name=None, date="", monotone_to=1.5, sum_window=None, auto=True, manual=None, reg_json=None, xml_scale="composite", verbose=True, lt_out="sample"):
    ctx = setup(bg_path, sample_path, auto, manual, verbose=False)
    hh = _h_from(h)
    E = _ecal_from(ecal)
    
    C_xml = E.get("coeffs_xml", E["coeffs"]) if xml_scale == "composite" else E["coeffs"]
    
    sw = sum_window or ctx["auto"]["sum_window"]
    
    res = ops.make_net(
        bg_path, sample_path, hh, C_xml, fw_points,
        out_zeroed, out_signed, name=name, date=date,
        monotone_to=monotone_to, sum_window=tuple(sw), lt_out=lt_out
    )
    
    # правка Claude (раунд 4): вердикт без модели образца; исключения и ступени -- из регистрации, окна -- через шкалу пары
    excluded, steps, converged, model = [], [], True, None
    if reg_json is not None:
        if isinstance(reg_json, (str, pathlib.Path)):
            reg_json = _load_json(reg_json)
        rj = reg_json.get("robust", {})
        excluded = rj.get("excluded_ch") or [[int(e[0]), int(e[1])] for e in rj.get("excluded", [])]
        steps, converged = rj.get("steps_ch", []), bool(rj.get("converged", True))
    windows = vd.windows_to_channels(vd.load_windows_kev(), E["coeffs"], ctx["n"])
    checks = vd.structure_verdict(res["resid"], res["var"], windows, **verdict_kw(ctx, excluded, steps))
    if verbose:
        if not converged:
            print("[ROBUST] ПРЕДУПРЕЖДЕНИЕ: устойчивая регистрация не сошлась -- итог не принимается")
        print(vd.format_table(checks, "[STRUCT] окна собственных структур (остаток без модели образца, базовая линия "
                              "по боковым полосам, окна по шкале пары):"))
        print("[ACC] сумма разности по всей области -- справочно (у образца со своими линиями не равна нулю)")
    
    result = dict(res)
    result["checks"] = checks
    result["checks_summary"] = vd.summarize(checks)
    result["converged"] = converged
    result["model"] = model
    result["excluded"] = excluded
    result["coeffs_xml"] = C_xml
    result["sum_window"] = list(sw)
    
    return result
