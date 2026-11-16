"""Совмещение по линиям для сильных образцов (метод nai): только положения линий образца и фона, без модели
формы отклика (форма зависит от геометрии, положения — нет). Режимы: peaks (≥ 2 линии в образце), peaks1 (одна линия,
усиление — из шкал файлов, систематика one_line_gain_sys). Выбор режима и флаг шкалы — api.run (run_auto — обёртка).
ВНЕ поддерживаемой области метода (итерация 4): результат помечается «шкала не проверена» (api.run, ключи
scale_verified/flag). peaks1 НЕ восстанавливает дрейф усиления: усиление берётся из шкал файлов, по линии — только сдвиг."""

import numpy as np
from gamma.calibration.intrinsic_shape.core import poly_eval, poly_deriv, poly_inverse
from . import bgcal, binned, checks, pipeline, surrogate


def _fit_ecal(recs, deg=2, shift=None):
    """Подгонка полиномиальной калибровки по принятым линиям (взвешенная по неопределённости энергии).

    recs — список записей линий с ключами ch, E, slope, ch_err, sys_kev;
    deg — степень полинома; shift — None или массив поправок к каналам (Монте-Карло).
    Возвращает список коэффициентов по возрастанию, дополненный нулями до 3.
    """
    ch = np.array([r["ch"] for r in recs], dtype=float)
    if shift is not None:
        ch = ch + np.asarray(shift, dtype=float)
    E = np.array([r["E"] for r in recs], dtype=float)
    slope = np.array([r["slope"] for r in recs], dtype=float)
    ch_err = np.array([r["ch_err"] for r in recs], dtype=float)
    sys_kev = np.array([r["sys_kev"] for r in recs], dtype=float)
    sigE = np.sqrt((slope * ch_err) ** 2 + sys_kev ** 2 + 0.01)
    dg = min(deg, len(recs) - 1)
    if dg < 0:
        return [0.0, 0.0, 0.0]
    p = np.polyfit(ch, E, dg, w=1.0 / sigE)
    # np.polyfit возвращает коэффициенты по убыванию степени; переворачиваем
    result = [float(v) for v in p[::-1]]
    while len(result) < 3:
        result.append(0.0)
    return result


def _ch_sigma(recs):
    """Стандартное отклонение канала линии: sqrt(ch_err^2 + (sys_kev/slope)^2)."""
    ch_err = np.array([r["ch_err"] for r in recs], dtype=float)
    sys_kev = np.array([r["sys_kev"] for r in recs], dtype=float)
    slope = np.array([r["slope"] for r in recs], dtype=float)
    return np.sqrt(ch_err ** 2 + (sys_kev / slope) ** 2)


def _map(ecal_s, ecal_b, xs, n_b, deg):
    """Построение отображения канал_образца → канал_фона по линиям.

    ecal_s — калибровка образца (коэфф. по возрастанию);
    ecal_b — калибровка фона;
    xs — каналы образца в рабочем диапазоне;
    n_b — число каналов фона;
    deg — степень полинома отображения.
    """
    Es = poly_eval(ecal_s, xs)
    xb = poly_inverse(ecal_b, Es, hi=n_b + 200)
    c = [float(v) for v in np.polyfit(xs, xb, deg)[::-1]]
    while len(c) < 3:
        c.append(0.0)
    return c


def sample_anchor_lines(profile, nuclides, fwhm_kev):
    """Формирует набор линий для калибровки образца.

    Возвращает объединение:
      1) якорных линий всех нуклидов (profile["anchor_lines"]);
      2) линий фона (profile["calib_lines"]), удалённых более чем на
         1.5·fwhm_kev(E) от всех энергий линий образца и от всех якорей образца.
    Дубликаты (|ΔE| < 0.5 кэВ) удаляются.
    """
    # Якорные линии всех нуклидов
    anchors = []
    for nuc in nuclides:
        for rec in profile["anchor_lines"][nuc]:
            anchors.append({"E": rec["E"], "sys_kev": rec["sys_kev"]})

    # Энергии линий образца
    sample_energies = pipeline.sample_lines(profile, nuclides)
    anchor_energies = [a["E"] for a in anchors]

    # Добавляем линии фона, не перекрытые линиями образца и якорями
    result = list(anchors)
    for line in profile["calib_lines"]:
        E = line["E"]
        margin = 1.5 * fwhm_kev(E)
        # Удалённость от всех линий образца
        if any(abs(E - E_s) <= margin for E_s in sample_energies):
            continue
        # Удалённость от всех якорей образца
        if any(abs(E - E_a) <= margin for E_a in anchor_energies):
            continue
        result.append({"E": E, "sys_kev": line.get("sys_kev", 0.1)})

    # Удаление дубликатов по энергии (|ΔE| < 0.5 кэВ)
    result.sort(key=lambda r: r["E"])
    dedup = []
    for r in result:
        if dedup and abs(r["E"] - dedup[-1]["E"]) < 0.5:
            continue
        dedup.append(r)
    return dedup


def _register_by_lines(bg, sm, profile, nuclides, n_mc=200, seed=0, verbose=False):
    """Совмещение по линиям (режим peaks / peaks1).

    Использует только положения принятых линий образца и фона.
    Возвращает dict результата в формате, совместимом с pipeline._run_pair.
    """
    # --- Шаг 1: подготовка спектров ---
    B = bg["counts"].astype(float).copy()
    B[-1] = 0.0
    S = sm["counts"].astype(float).copy()
    S[-1] = 0.0
    n_b = len(B)
    n_s = len(S)
    k = sm["live"] / bg["live"]

    # --- Шаг 2: калибровка фона ---
    cal_b = bgcal.calibrate(
        B, bg["cal"], profile["calib_lines"],
        profile["e_strong"], profile["fwhm_rel662"],
        profile["search_rel"], profile.get("min_z", 5.0),
    )
    fw_ch = bgcal.fwhm_ch_fn(cal_b)
    e2ch = lambda E: float(bgcal.energy_to_ch(cal_b, E))
    fwhm_kev = lambda E: float(fw_ch(e2ch(E)) * poly_deriv(cal_b["ecal"], e2ch(E)))

    # --- Шаг 3: калибровка образца ---
    lines_s = sample_anchor_lines(profile, nuclides, fwhm_kev)
    e_strong_s = profile["anchor_strong"][nuclides[0]]

    try:
        cal_s = bgcal.calibrate(
            S, sm["cal"], lines_s, e_strong_s,
            profile["fwhm_rel662"], profile["search_rel"],
            profile.get("min_z", 5.0),
        )
    except ValueError:
        cal_s = None

    lo, hi = profile["fit_range_kev"]

    # Принятые записи линий
    rb = [r for r in cal_b["lines"] if r.get("accepted", False)]

    # --- Определение режима и построение c ---
    if cal_s is None or cal_s["n_accepted"] < 2:
        # ===== Режим peaks1: одна линия =====
        ratio, rp = bgcal.find_strong(
            S, sm["cal"], e_strong_s,
            profile["fwhm_rel662"], profile["search_rel"],
        )
        xl = float(rp["centroid"])
        xl_err = float(rp["centroid_err"])
        h0 = pipeline.initial_map(bg, sm, n_b, n_s, profile["fit_range_kev"])
        h0 = list(h0)
        yl = e2ch(e_strong_s)
        c = list(h0)
        c[0] += yl - poly_eval(h0, xl)
        # Шкала образца — шкала файла
        ecal_sample = None
        sample_cal_lines = []
        one_line = {"E": float(e_strong_s), "x": xl, "x_err": xl_err}
        n_acc_s = 1
    else:
        # ===== Режим peaks: ≥ 2 линии =====
        xs = np.arange(n_s, dtype=float)
        Es = poly_eval(cal_s["ecal"], xs)
        xs = xs[(Es >= lo) & (Es <= hi)]
        deg = profile.get("peaks_deg", 1)
        c = _map(cal_s["ecal"], cal_b["ecal"], xs, n_b, deg)
        ecal_sample = cal_s["ecal"]
        sample_cal_lines = cal_s["lines"]
        one_line = None
        n_acc_s = cal_s["n_accepted"]
        rs = [r for r in cal_s["lines"] if r.get("accepted", False)]

    # --- Шаг 5: контрольные точки ---
    test_kev = profile["test_kev"]
    yb = np.array([e2ch(E) for E in test_kev], dtype=float)
    x_test = surrogate.inverse_h(c, yb)
    x_test = np.asarray(x_test, dtype=float)

    # --- Шаг 6: Монте-Карло ---
    rng = np.random.default_rng(seed)
    H = np.zeros((n_mc, len(x_test)), dtype=float)

    if one_line is not None:
        # MC для peaks1
        sigma_rb = _ch_sigma(rb)
        for i in range(n_mc):
            sh_b = rng.normal(0.0, sigma_rb, size=len(rb))
            eb = _fit_ecal(rb, shift=sh_b)
            yl_i = float(poly_inverse(eb, e_strong_s, hi=n_b + 200))
            xl_i = xl + rng.normal(0.0, xl_err)
            g_i = 1.0 + rng.normal(0.0, profile.get("one_line_gain_sys", 0.005))
            ci = [h0[0] * g_i, h0[1] * g_i, h0[2] * g_i]
            ci[0] += yl_i - poly_eval(ci, xl_i)
            H[i] = poly_eval(ci, x_test)
    else:
        # MC для peaks
        sigma_rb = _ch_sigma(rb)
        sigma_rs = _ch_sigma(rs)
        deg = profile.get("peaks_deg", 1)
        for i in range(n_mc):
            sh_b = rng.normal(0.0, sigma_rb, size=len(rb))
            sh_s = rng.normal(0.0, sigma_rs, size=len(rs))
            eb = _fit_ecal(rb, shift=sh_b)
            es = _fit_ecal(rs, shift=sh_s)
            ci = _map(es, eb, xs, n_b, deg)
            H[i] = poly_eval(ci, x_test)

    sigma_h = H.std(axis=0, ddof=1) if n_mc > 1 else np.zeros(len(x_test))

    # --- Шаг 7: контрольные точки на сетке фона ---
    Eb = poly_eval(cal_b["ecal"], np.arange(n_b, dtype=float))
    mask = (Eb >= lo) & (Eb <= hi)
    prob = binned.BinnedProblem(
        S, B, k, mask, c, None,
        ref_ch=n_b / 2.0, fw_ch=fw_ch,
        bin_fwhm=profile.get("bin_fwhm", 0.5),
    )
    ev = prob.evaluate_channels(np.zeros(4))

    # --- Шаг 8: оконные тесты ---
    lines = [E for E in pipeline.sample_lines(profile, nuclides) if lo - 50 <= E <= hi + 50]
    wl = [w["E"] for w in profile["null_windows"]]
    tags = {w["E"]: w["tag"] for w in profile["null_windows"]}
    tests = checks.window_tests(
        ev["net"], ev["var"], e2ch, fwhm_kev,
        wl, 1.0, lines, 1.5, True, tags,
    )
    vd = checks.verdict(tests, profile.get("z_thr", 3.0), None)

    # --- Шаг 9: результат ---
    mode = "peaks1" if one_line is not None else "peaks"
    out = {
        "mode": mode,
        "bg": bg["name"],
        "sample": sm["name"],
        "k": float(k),
        "h0": c,
        "c": c,
        "p": {"d": 0.0, "gN": 0.0, "q": 0.0, "s": 0.0},
        "err": {},
        "chi2ndf": None,
        "test_kev": [float(E) for E in test_kev],
        "h_test": [float(v) for v in yb],
        "x_test": [float(v) for v in x_test],
        "sigma_h_test": [float(v) for v in sigma_h],
        "ecal_bg": cal_b["ecal"],
        "ecal_sample": ecal_sample,
        "R662": cal_b["R662"],
        "fwc": cal_b["fwc"],  # список или None — не приводим к float
        "bg_lines": cal_b["lines"],
        "sample_cal_lines": sample_cal_lines,
        "lines_sample": lines,
        "tests": tests,
        "verdict": vd,
        # ревизия 6: внутренний результат по умолчанию не проверен; флаг снимает только api.run
        "scale_verified": False,
        "flag": "шкала не проверена",
        "status": "experimental",
        "resid_out": None,
        "arrays": {
            "net": ev["net"],
            "var": ev["var"],
            "fit_lin": ev["fit_lin"],
            "mask": mask,
            "Eb": Eb,
        },
        "lt": [sm["live"], sm["real"]],
    }
    if one_line is not None:
        out["one_line"] = one_line

    if verbose:
        print(f"[NAI-PEAKS] режим={mode}")
        print(f"[NAI-PEAKS] c = {c}")
        print(f"[NAI-PEAKS] sigma_h_test = {list(sigma_h)}")
        print(f"[NAI-PEAKS] линии фона: {len(rb)}, образца: {n_acc_s}")

    return out


def bg_fraction(bg, sm, profile):
    """Доля фона в спектре образца в рабочем диапазоне.

    f = (live_s / live_b) · ΣB[маска_b] / max(ΣS[маска_s], 1)
    Маски строятся по fit_range_kev на шкале файла каждого спектра.
    Последний канал каждого спектра не учитывается.
    """
    lo, hi = profile["fit_range_kev"]

    # Шкала файла образца
    Es = poly_eval(sm["cal"], np.arange(len(sm["counts"]), dtype=float))
    mask_s = (Es >= lo) & (Es <= hi)

    # Шкала файла фона
    Eb = poly_eval(bg["cal"], np.arange(len(bg["counts"]), dtype=float))
    mask_b = (Eb >= lo) & (Eb <= hi)

    S_arr = sm["counts"].astype(float).copy()
    B_arr = bg["counts"].astype(float).copy()
    S_arr[-1] = 0.0
    B_arr[-1] = 0.0

    num = (sm["live"] / bg["live"]) * float(B_arr[mask_b].sum())
    den = max(float(S_arr[mask_s].sum()), 1.0)
    return float(num / den)


def run_auto(bg, sm, profile, nuclides, verbose=False, force=None, kind="sample"):
    """Совместимость: делегирует api.run (выбор режима, флаг шкалы, вердикт). Ревизия 5."""
    from . import api
    return api.run(bg, sm, profile, nuclides, kind=kind, force=force, verbose=verbose)
