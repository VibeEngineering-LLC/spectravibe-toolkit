"""Модуль pipeline пакета gamma.calibration.nai_shape.

Реализует полный цикл калибровки формы пиков NAI-детектора:
подготовка спектров, начальное отображение, фит нелинейной модели,
контрольные окна и запись результатов (CSV, XML, JSON).
"""

import json
import os
import pathlib
import tempfile

import numpy as np

from gamma.calibration.intrinsic_shape.core import load_spectrum, poly_eval, poly_deriv, poly_inverse
from gamma.calibration.intrinsic_shape.fwhm import fit_fwhm
from gamma.calibration.intrinsic_shape import netxml
from gamma.io.becqmoni_xml import write_becqmoni_xml
from gamma.spectrum import Spectrum
from . import bgcal, model, checks, surrogate, binned, search


# Индексы параметров модели p = [d, gN, q, s]
PARAM = {"d": 0, "gN": 1, "q": 2, "s": 3}


# ревизия 7: токен _API_TOKEN убран (обходился импортом); синонимы старых имён профилей
_PROFILE_ALIASES = {"gamma1s.json": "nai_gamma1s.json", "gs2020.json": "nai_gs2020.json"}


def _alias_profile(path):
    """Старое имя профиля → новое, с предупреждением DeprecationWarning."""
    p = str(path)
    base = os.path.basename(p)
    if base in _PROFILE_ALIASES and not os.path.exists(p):
        import warnings
        new = os.path.join(os.path.dirname(p), _PROFILE_ALIASES[base])
        warnings.warn("профиль %s переименован в %s" % (base, _PROFILE_ALIASES[base]), DeprecationWarning, stacklevel=3)
        return new
    return p
def load_profile(path):
    """Загружает JSON-профиль конфигурации (utf-8)."""
    path = _alias_profile(path)  # ревизия 6: старые имена профилей — синонимы
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def sample_lines(profile, nuclides=None):
    """Возвращает отсортированный список уникальных энергий (кэВ) линий образца.

    Parameters
    ----------
    profile : dict
        Профиль, содержащий ключ ``sample_lines`` — словарь
        {нуклид: [энергии кэВ]}.
    nuclides : list | None
        Если None — берутся все нуклиды словаря.

    Returns
    -------
    list[float]

    Raises
    ------
    KeyError
        Если запрошен нуклид отсутствует в ``profile["sample_lines"]``.
    """
    lines_dict = profile["sample_lines"]
    if nuclides is None:
        nuclides = list(lines_dict.keys())

    energies = set()
    for nu in nuclides:
        if nu not in lines_dict:
            raise KeyError(
                f"Нуклид {nu!r} не найден в sample_lines. "
                f"Доступные: {sorted(lines_dict.keys())}"
            )
        energies.update(float(e) for e in lines_dict[nu])

    return sorted(energies)


def initial_map(bg, sm, n_b, n_s, fit_range, deg=1):
    """Начальное квадратичное отображение «канал образца → канал фона».

    Строится по шкалам ФАЙЛОВ на пересечении энергетических диапазонов.

    Returns
    -------
    list[float]
        Коэффициенты [a0, a1, a2] (обратный порядок полинома).
    """
    xs = np.arange(n_s, dtype=float)
    E = poly_eval(sm["cal"], xs)

    lo, hi = fit_range[0], fit_range[1]
    sel = (E >= lo) & (E <= hi)
    xs_sel = xs[sel]
    E_sel = E[sel]

    xb = poly_inverse(bg["cal"], E_sel, hi=n_b + 200)

    # Решение по данным (2026-11-15): только ЛИНЕЙНАЯ часть. Квадратичная часть отношения шкал файлов -- артефакт
    # перекалибровки файла (Гамма-1С, вода_08 против вода_16: h0 2-й степени давал смещение 1 кан у K-40 при
    # σ_h 0.03 кан); физический дрейф усиления -- масштаб + смещение, нелинейность -- свободным q.
    coeffs = np.polyfit(xs_sel, xb, deg)
    h0 = [float(v) for v in reversed(coeffs)] + [0.0] * (2 - deg)
    return h0


def _run_pair(bg, sm, profile, mode="sample", nuclides=None, free=None, verbose=True):
    """Основная функция: фит пары спектров (фон / образец или фон / фон).

    Parameters
    ----------
    bg, sm : dict
        Результаты ``load_spectrum`` (фон и образец соответственно).
    profile : dict
        Конфигурация из ``load_profile``.
    mode : str
        ``"sample"`` — образец со своими линиями;
        ``"background"`` — фон против фона.
    nuclides : list | None
        Ограничение нуклидов (только mode="sample").
    free : list | None
        Свободные параметры (перекрывает профиль).
    verbose : bool
        Печать промежуточных результатов.

    Returns
    -------
    dict
        Сериализуемый результат (см. спецификацию шага 13).
    """
    # ------------------------------------------------------------------
    # 1. Подготовка массивов
    # ------------------------------------------------------------------
    B = bg["counts"].astype(float).copy()
    B[-1] = 0.0
    S = sm["counts"].astype(float).copy()
    S[-1] = 0.0
    n_b = len(B)
    n_s = len(S)
    k = sm["live"] / bg["live"]

    # ------------------------------------------------------------------
    # 2. Калибровка фона
    # ------------------------------------------------------------------
    cal_b = bgcal.calibrate(
        B,
        bg["cal"],
        profile["calib_lines"],
        profile["e_strong"],
        profile["fwhm_rel662"],
        profile["search_rel"],
        profile.get("min_z", 5.0),
    )
    fw_ch = bgcal.fwhm_ch_fn(cal_b)
    ecal = cal_b["ecal"]
    e2ch = lambda E: float(bgcal.energy_to_ch(cal_b, E))
    fwhm_kev = lambda E: float(fw_ch(e2ch(E)) * poly_deriv(ecal, e2ch(E)))

    # ------------------------------------------------------------------
    # 3. Начальное отображение
    # ------------------------------------------------------------------
    h0 = initial_map(bg, sm, n_b, n_s, profile["fit_range_kev"])

    # ------------------------------------------------------------------
    # 4. Маска по энергетическому диапазону
    # ------------------------------------------------------------------
    lo, hi = profile["fit_range_kev"]
    Eb = poly_eval(ecal, np.arange(n_b, dtype=float))
    mask = (Eb >= lo) & (Eb <= hi)
    idxs = np.where(mask)[0]
    ch_lo = int(idxs[0])
    ch_hi = int(idxs[-1])

    # ------------------------------------------------------------------
    # 5. Базис (только для образца)
    # ------------------------------------------------------------------
    if mode == "sample":
        all_lines = sample_lines(profile, nuclides)
        lines = [E for E in all_lines if lo - 50.0 <= E <= hi + 50.0]
        centers = np.array([e2ch(E) for E in lines], dtype=float)
        sigmas = np.array([fw_ch(c) / 2.3548 for c in centers], dtype=float)
        G = model.gauss_cols(n_b, centers, sigmas)

        # Узлы шапок (hat functions)
        knots = [ch_lo]
        x = float(ch_lo)
        while True:   # ручная правка: в генерации узлы внутри цикла не добавлялись (было 2 шапки)
            x += profile["hat_spacing_fwhm"] * fw_ch(x)
            if x >= ch_hi:
                break
            knots.append(x)
        knots.append(ch_hi)
        H = model.hats(n_b, np.array(knots, dtype=float))
        basis = np.hstack([G, H])
    elif profile.get("hats_background", False):
        # фон против фона: гладкая подложка поглощает различие континуума (радон, загрузка); совмещение держится на пиках и краях масштаба меньше шага шапок
        spacing = profile.get("hat_spacing_bg_fwhm", 4.0)
        knots = [float(ch_lo)]
        x = float(ch_lo)
        while True:
            x += spacing * fw_ch(x)
            if x >= ch_hi:
                break
            knots.append(x)
        knots.append(float(ch_hi))
        basis = model.hats(n_b, np.array(knots, dtype=float))
        lines = []
    else:
        basis = None
        lines = []

    # ------------------------------------------------------------------
    # 6. Задача оптимизации
    # ------------------------------------------------------------------
    prob = binned.BinnedProblem(S, B, k, mask, h0, basis, ref_ch=n_b / 2.0, fw_ch=fw_ch, bin_fwhm=profile.get("bin_fwhm", 0.125))

    # ------------------------------------------------------------------
    # 7. Свободные параметры
    # ------------------------------------------------------------------
    names = free or profile["free_" + mode]
    idx = [PARAM[n] for n in names]

    # ------------------------------------------------------------------
    # 8. Сетка начальных значений и фит
    # ------------------------------------------------------------------
    step_g = profile["grid_step_fwhm"] * fw_ch(ch_hi) * n_b / ch_hi
    span_g = profile["gain_span"] * n_b
    g_grid = np.arange(-span_g, span_g + step_g / 2.0, step_g)

    step_d = profile["grid_step_fwhm"] * fw_ch(ch_lo)
    d_grid = np.arange(
        -profile["offset_span_ch"],
        profile["offset_span_ch"] + step_d / 2.0,
        step_d,
    )

    # ручная правка: сетка по развязанным сдвигам в двух опорных точках (20 % и 85 % диапазона) вместо (d, gN) -- search.py
    res = search.anchor_grid_fit(prob, idx, n_b, (ya := ch_lo + 0.2 * (ch_hi - ch_lo)), (yb := ch_lo + 0.85 * (ch_hi - ch_lo)), fw_ch(ya), fw_ch(yb), profile["offset_span_ch"] + profile["gain_span"] * ya, profile["offset_span_ch"] + profile["gain_span"] * yb, profile["grid_step_fwhm"], profile.get("n_starts", 5))
    p = res["p"]

    # ------------------------------------------------------------------
    # 9. Ковариация и оценка
    # ------------------------------------------------------------------
    cv = model.covariance(prob, p, idx)
    ev = prob.evaluate_channels(p)
    c = ev["c"]

    # ------------------------------------------------------------------
    # 10. Контрольные каналы
    # ------------------------------------------------------------------
    test_kev = profile["test_kev"]
    yb = np.array([e2ch(E) for E in test_kev], dtype=float)
    xs_test = surrogate.inverse_h(c, yb)
    sh = model.h_sigma(p, cv["cov"], idx, h0, n_b, xs_test)

    # ------------------------------------------------------------------
    # 11. Окна контроля и вердикт
    # ------------------------------------------------------------------
    tags = {w["E"]: w["tag"] for w in profile["null_windows"]}
    wl = [w["E"] for w in profile["null_windows"]]
    tests = checks.window_tests(
        ev["net"],
        ev["var"],
        e2ch,
        fwhm_kev,
        wl,
        half_fwhm=1.0,
        exclude_kev=lines,
        exclude_fwhm=1.5,
        sideband=(mode == "sample"),
        tags=tags,
    )

    if mode == "background":
        use_sum = [E for E in wl if E not in profile.get("radon_variable_kev", [])]
    else:
        use_sum = None

    vd = checks.verdict(tests, profile.get("z_thr", 3.0), use_sum)

    # Остаток модели вне линий образца
    out = mask.copy()
    for E in lines:
        half = 1.5 * fwhm_kev(E)

        out &= ~(np.abs(Eb - E) < half)   # ручная правка: сравнение энергии с энергией (было с каналом)

    rm = checks.masked_chi2(ev["z"], out)
    rz = (ev["net"] - ev["fit_lin"])[out].sum() / np.sqrt(ev["var"][out].sum())

    # ------------------------------------------------------------------
    # 12. Печать (verbose)
    # ------------------------------------------------------------------
    if verbose:
        _print_report(p, cv, idx, names, res, ev, sh, xs_test, tests, vd, rm, rz, test_kev, yb)

    # ------------------------------------------------------------------
    # 13. Формирование результата
    # ------------------------------------------------------------------
    p_dict = {}
    for i, n in enumerate(["d", "gN", "q", "s"]):
        p_dict[n] = float(p[i])

    err_dict = {}
    for i, n in enumerate(["d", "gN", "q", "s"]):
        if n in names:
            err_dict[n] = float(cv["err"][names.index(n)])   # ручная правка: индекс в списке свободных

    result = {
        "mode": mode,
        "bg": bg["name"],
        "sample": sm["name"],
        "k": float(k),
        "h0": h0,
        "c": list(c),
        "free": names,
        "p": p_dict,
        "err": err_dict,
        "chi2": float(res["chi2"]),
        "ndf": int(res["ndf"]),
        "chi2ndf": float(cv["chi2ndf"]),
        "minima": res["minima"][:5],
        "test_kev": test_kev,
        "h_test": list(yb),
        "x_test": list(xs_test),
        "sigma_h_test": list(sh),
        "ecal_bg": ecal,
        "R662": cal_b["R662"],
        "fwc": cal_b["fwc"],
        "bg_lines": cal_b["lines"],
        "lines_sample": lines,
        "tests": tests,
        "verdict": vd,
        # ревизия 6: внутренний результат по умолчанию не проверен; флаг снимает только api.run
        "scale_verified": False,
        "flag": "шкала не проверена",
        "status": "experimental",
        "resid_out": {
            "chi2_per_ch": float(rm["chi2_per_ch"]),
            "n": int(rm["n"]),
            "z_sum": float(rz),
        },
        "arrays": {
            "net": ev["net"],
            "var": ev["var"],
            "fit_lin": ev["fit_lin"],
            "mask": mask,
            "Eb": Eb,
        },
        "lt": [sm["live"], sm["real"]],
    }
    return result


def _print_report(p, cv, idx, names, res, ev, sh, xs_test, tests, vd, rm, rz, test_kev, yb):
    """Выводит читаемый отчёт в stdout с префиксом ``[NAI]``."""
    all_names = ["d", "gN", "q", "s"]
    print("[NAI] === Параметры ===")
    for i, n in enumerate(all_names):
        val = float(p[i])
        if n in names:
            err = float(cv["err"][names.index(n)])
            print(f"[NAI]   {n:4s} = {val:12.6g} ± {err:.4g}")
        else:
            print(f"[NAI]   {n:4s} = {val:12.6g}  (фикс.)")

    print(f"[NAI] χ² = {res['chi2']:.4f},  ndf = {res['ndf']},  χ²/ndf = {cv['chi2ndf']:.4f}")

    print("[NAI] === Контрольные каналы (h ± σ) ===")
    for i, (E, y, x, s) in enumerate(zip(test_kev, yb, xs_test, sh)):
        print(f"[NAI]   E={E:8.1f} кэВ → h={y:8.2f}, x={x:8.2f} ± {s:.3f}")

    print("[NAI] === Окна контроля ===")
    print("[NAI]   {:>8s}  {:<12s}  {:>8s}  {:>8s}  {:>8s}".format(
        "E", "tag", "z_sum", "z_dip", "excluded"))
    for t in tests:
        print("[NAI]   {:8.1f}  {:<12s}  {:8.2f}  {:8.2f}  {:8s}".format(
            t["E"], t.get("tag", ""), t.get("z_sum", 0.0),
            t.get("z_dip", 0.0), "да" if t.get("excluded") else "нет"))

    print(f"[NAI] Вердикт: {vd}")
    print(f"[NAI] Остаток вне линий: χ²/ndf={rm['chi2_per_ch']:.3f} (n={rm['n']}), z_sum={rz:.3f}")


def _jsonable(o):
    """Преобразует numpy-объекты в сериализуемые типы Python.

    ndarray → list, numpy число → float/int, numpy bool → bool.
    """
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def _write_outputs(res, sm, out_dir, stem, date=""):
    """Записывает CSV, два XML и JSON-сводку в каталог ``out_dir``.

    Returns
    -------
    dict
        Путь к каждому файлу: ``{"csv", "xml_zeroed", "xml_signed", "json"}``.
    """
    # ревизия 7: та же проверка согласованности, что в api.write_outputs; результат без статуса api.run не пишется
    from . import api as _api
    _api._check_consistent(res)
    os.makedirs(out_dir, exist_ok=True)

    net = res["arrays"]["net"]
    var = res["arrays"]["var"]
    Eb = res["arrays"]["Eb"]
    n_b = len(net)
    lt, rt = res["lt"]

    # ------------------------------------------------------------------
    # 2. CSV
    # ------------------------------------------------------------------
    csv_path = os.path.join(out_dir, f"{stem}_net_cps.csv")
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("channel,energy_kev,net_cps,sigma_cps\n")
        for ch in range(n_b):
            e_kev = Eb[ch]
            net_cps = net[ch] / lt
            sig_cps = np.sqrt(var[ch]) / lt
            f.write("%d,%.3f,%.6g,%.6g\n" % (ch, e_kev, net_cps, sig_cps))

    # ------------------------------------------------------------------
    # 3. ПШПВ
    # ------------------------------------------------------------------
    fw_peaks = [
        (L["ch"], L["fwhm_ch"])
        for L in res["bg_lines"]
        if L["accepted"]
    ]
    fwc = res["fwc"]
    if fwc is None:
        xs = np.linspace(0.1 * n_b, 0.9 * n_b, 5)
        _fw_fn = bgcal.fwhm_ch_fn({
            "fwc": None,
            "R662": res["R662"],
            "ecal": res["ecal_bg"],
            "n": n_b,
        })
        f_vals = [_fw_fn(int(x)) for x in xs]
        fwc = fit_fwhm(n_b, x=xs, f=f_vals)

    # ------------------------------------------------------------------
    # 4. Шаблон XML
    # ------------------------------------------------------------------
    src = None
    sm_path = sm["path"]
    if sm_path.endswith(".xml") and len(sm["counts"]) == n_b:  # ревизия 6: читатель отбрасывает канал переполнения (n_read = n-1), шаблон — исходный XML образца
        src = sm_path
    else:
        spec = Spectrum(
            counts=np.zeros(n_b, dtype=np.int64),
            live_time=lt,
            real_time=rt,
            sample_id=stem,
            energy_cal=tuple(res["ecal_bg"]),
            valid_pulse_count=0,
            total_pulse_count=0,
        )
        tmp_path = os.path.join(tempfile.gettempdir(), f"{stem}_tmp_template.xml")
        write_becqmoni_xml(spec, tmp_path)
        src = tmp_path

    # ------------------------------------------------------------------
    # 5. Два XML-файла
    # ------------------------------------------------------------------
    xml_zeroed_path = os.path.join(out_dir, f"{stem}_net_zeroed.xml")
    xml_signed_path = os.path.join(out_dir, f"{stem}_net_signed.xml")

    data_zeroed = np.clip(net, 0.0, None)
    netxml.write_xml(
        src, xml_zeroed_path, data_zeroed,
        signed=False,
        C=res["ecal_bg"],
        H=res["c"],
        k=res["k"],
        fwc=fwc,
        fw_peaks=fw_peaks,
        name=stem + " (обнулённый)",
        date=date,
        times=(lt, rt),
    )

    netxml.write_xml(
        src, xml_signed_path, net,
        signed=True,
        C=res["ecal_bg"],
        H=res["c"],
        k=res["k"],
        fwc=fwc,
        fw_peaks=fw_peaks,
        name=stem + " (со знаком)",
        date=date,
        times=(lt, rt),
    )

    # ------------------------------------------------------------------
    # 6. JSON-сводка
    # ------------------------------------------------------------------
    json_path = os.path.join(out_dir, f"{stem}_summary.json")
    summary = {k: v for k, v in res.items() if k != "arrays"}
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=1, default=_jsonable)

    # ------------------------------------------------------------------
    # 7. Возврат путей
    # ------------------------------------------------------------------
    return {
        "csv": csv_path,
        "xml_zeroed": xml_zeroed_path,
        "xml_signed": xml_signed_path,
        "json": json_path,
    }
