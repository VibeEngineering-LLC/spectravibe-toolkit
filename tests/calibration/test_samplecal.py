"""Тесты калибровки по пробе (samplecal): синтетика с известной шкалой, перенос набора реперов,
источник энергии, отбраковка NaN, степень по числу реперов, XML для BecqMoni. Сгенерировано gen_code.py
(qwen3.8:27b) по спецификации, 2026-10-10."""

import re
import math
import pathlib
import numpy as np
import pytest
from gamma.calibration.intrinsic_shape.core import load_spectrum, poly_eval, poly_deriv, poly_inverse, rebin_conserve
from gamma.calibration.intrinsic_shape import samplecal as sc
from gamma.calibration.intrinsic_shape import pipeline as pl
from gamma.calibration.intrinsic_shape import lowscale as lw
from gamma.calibration.intrinsic_shape import fwhm as F

FIX = pathlib.Path(__file__).resolve().parents[1] / "fixtures" / "intrinsic_shape"
BG = FIX / "labr3_background.xml"
ECAL_REF = [-0.48496698, 1.35007262, 4.58971597e-04]
ECAL_LIN = [-0.48, 1.62]   # линейная шкала для проверки степени по числу реперов
THR = 16   # первый ненулевой канал фона фикстуры (порог прибора)
LINES = [("Pb-214", 351.932), ("Bi-214", 609.312), ("Cs-137", 661.657), ("Bi-214", 1120.287), ("Bi-214", 1764.494), ("Tl-208", 2614.511)]


def _stretch(counts, a, b, n):
    """Растянуть спектр: сдвиг a, масштаб b, новая длина n."""
    edges = a + b * (np.arange(len(counts) + 1) - 0.5)
    return rebin_conserve(np.asarray(counts, float), edges, n)


def _peaks(n, lines, area=30000.0, scale=ECAL_REF, a=0.0, b=1.0):
    """Сгенерировать синтетический спектр с гауссовыми пиками.
    
    Возвращает (y, centers), где centers — список центров в порядке lines.
    """
    x = np.arange(n, dtype=float)
    y = np.zeros(n, dtype=float)
    centers = []
    for nuc, E in lines:
        # Инверсия шкалы: энергия -> канал
        ch_ref = float(poly_inverse(scale, E))
        # Центр в новых координатах
        c = a + b * ch_ref
        # ПШПВ в новых координатах
        w = b * sc.fwhm_guess(c / b, n)
        s = w / 2.355
        # Гауссиана
        gauss = area / (math.sqrt(2 * math.pi) * s) * np.exp(-0.5 * ((x - c) / s) ** 2)
        y += gauss
        centers.append(c)
    return y, centers


def _sigma_E(scale, r):
    """Оценка погрешности энергии для репера r."""
    return math.hypot(r["dE"], float(poly_deriv(scale, r["ch"])) * r["ch_err"])


def _repers_by_click(counts, centers, lines, shift=2.0):
    """Построить список реперов, используя fit_click для каждого пика."""
    repers = []
    for c, (nuc, E) in zip(centers, lines):
        f = sc.fit_click(counts, c + shift, sc.fwhm_guess(c, len(counts)))
        assert f is not None, f"fit_click failed for {nuc} at {E}"
        repers.append(sc.reper_from_choice(f, E, nuc))
    return repers


def _write_like_bg(tmp_path, name, counts):
    """Записать спектр в XML-файл, имитирующий структуру BG."""
    t = BG.read_text("utf-8")
    vals = iter(np.round(np.maximum(counts, 0)).astype(int))
    result = re.sub(
        r"<DataPoint>[^<]*</DataPoint>",
        lambda m: "<DataPoint>%d</DataPoint>" % next(vals),
        t
    )
    path = tmp_path / name
    path.write_text(result, encoding="utf-8")
    return path


def test_samplecal_synthetic_known_scale_recovered():
    """Синтетика с известной шкалой: шкала должна быть восстановлена."""
    n = 2048
    x = np.arange(n, dtype=float)
    rng = np.random.default_rng(7)
    cont = 400.0 * np.exp(-x / 400.0) + 20.0
    y, centers = _peaks(n, LINES)
    counts = rng.poisson(cont + y).astype(float)
    counts[:THR] = 0.0   # правка Claude: порог прибора, как у боевого спектра

    rep = _repers_by_click(counts, centers, LINES)
    res = sc.scale_from_repers(rep, counts, ECAL_REF, verbose=False)

    assert not res.get("rejected"), f"Результат отклонен: {res}"
    assert len(res["coeffs"]) == 3, f"Ожидалось 3 коэффициента, получено {len(res['coeffs'])}"

    # Проверка точности восстановления энергии
    for c, (nuc, E), r in zip(centers, LINES, rep):
        E_rec = float(poly_eval(res["coeffs"], c))
        sigma = _sigma_E(res["coeffs"], r)
        assert abs(E_rec - E) <= 3.0 * sigma, (
            f"Энергия {nuc} {E} не восстановлена: получено {E_rec}, σ={sigma}"
        )

    # Все реперы должны быть из библиотеки
    assert all(r["source"] == sc.SRC_LIB for r in rep), "Не все реперы имеют источник SRC_LIB"


def test_samplecal_json_transfer_same_scale(tmp_path):
    """Перенос набора реперов через JSON: шкала должна сохраниться."""
    d = load_spectrum(BG)
    B = d["counts"].astype(float)
    cal = list(d["cal"])
    n = len(B)

    y, centers = _peaks(n, LINES[:5])
    s1 = B + y
    rep1 = _repers_by_click(s1, centers, LINES[:5])
    intr, anc1 = sc.intrinsic_repers(s1, cal)
    rep1 = rep1 + intr

    # Сохранение и загрузка JSON
    p = tmp_path / "repers.json"
    sc.save_json(p, rep1, "LaBr3 test", n, anc1)
    data = sc.load_json(p)

    # Растяжение спектра
    s2 = _stretch(s1, 3.0, 1.05, n)
    rep2, m = sc.apply_to(data, s2, cal, verbose=False)

    # Построение шкал
    sc1 = sc.scale_from_repers(rep1, s1, cal, verbose=False)
    sc2 = sc.scale_from_repers(rep2, s2, cal, verbose=False)

    assert not sc1.get("rejected"), f"Шкала 1 отклонена: {sc1}"
    assert not sc2.get("rejected"), f"Шкала 2 отклонена: {sc2}"

    # Проверка параметров растяжения
    assert abs(m[1] - 1.05) < 0.01, f"Масштаб: ожидалось 1.05, получено {m[1]}"
    assert abs(m[0] - 3.0) < 1.5, f"Сдвиг: ожидалось 3.0, получено {m[0]}"

    # Проверка согласованности энергий для реперов, не из LaBr
    count_non_labr = 0
    for r in rep2:
        if r["source"] != sc.SRC_LABR:
            count_non_labr += 1
            y2 = r["ch"]
            E2 = poly_eval(sc2["coeffs"], y2)
            # Обратное преобразование координаты
            y1 = (y2 - 3.0) / 1.05
            E1 = poly_eval(sc1["coeffs"], y1)
            sigma = _sigma_E(sc2["coeffs"], r)
            assert abs(E2 - E1) <= 3.0 * sigma, (
                f"Несогласованность энергий: E2={E2}, E1={E1}, σ={sigma}"
            )

    assert count_non_labr == 5, f"Ожидалось 5 реперов не из LaBr, получено {count_non_labr}"


def test_samplecal_manual_energy_marked_and_kept_in_json(tmp_path):
    """Ручная энергия должна быть помечена и сохранена в JSON."""
    fit = {"ch": 500.0, "ch_err": 0.1, "fwhm": 9.0, "signif": 50.0}
    
    r_lib = sc.reper_from_choice(fit, 661.657, "Cs-137")
    r_man = sc.reper_from_choice(fit, 700.0, "Cs-137")
    r_none = sc.reper_from_choice(fit, 1234.5)

    assert r_lib["source"] == sc.SRC_LIB, f"Ожидалось SRC_LIB, получено {r_lib['source']}"
    assert r_man["source"] == sc.SRC_MANUAL, f"Ожидалось SRC_MANUAL, получено {r_man['source']}"
    assert r_none["source"] == sc.SRC_MANUAL, f"Ожидалось SRC_MANUAL, получено {r_none['source']}"
    assert sc.SRC_MANUAL == "ручная", f"Ожидалось 'ручная', получено {sc.SRC_MANUAL}"

    # Сохранение и загрузка
    anc = {k: {"x": 1.0, "visible": False} for k in ("a37", "edge", "sum")}
    sc.save_json(tmp_path / "r.json", [r_lib, r_man, r_none], "dev", 2048, anc)
    data = sc.load_json(tmp_path / "r.json")

    sources = [r["source"] for r in data["repers"]]
    assert sources == [sc.SRC_LIB, sc.SRC_MANUAL, sc.SRC_MANUAL], (
        f"Источники не сохранены: {sources}"
    )


def test_samplecal_nan_centroid_dropped():
    """Репер с NaN центроидом должен быть отброшен."""
    rep = [
        {
            "name": "p%d" % i,
            "ch": float(poly_inverse(ECAL_REF, E)),
            "ch_err": 0.05,
            "E": E,
            "dE": 0.05,
            "source": sc.SRC_LIB,
            "nuclide": "X",
            "fwhm": 8.0,
            "signif": 50.0
        }
        for i, E in enumerate((351.932, 661.657, 1120.287, 1764.494))
    ]
    # Добавляем плохой репер с NaN
    rep.append({
        "name": "bad",
        "ch": float("nan"),
        "ch_err": float("nan"),
        "E": 609.312,
        "dE": 0.05,
        "source": sc.SRC_LIB,
        "nuclide": "Bi-214",
        "fwhm": float("nan"),
        "signif": float("nan")
    })

    counts = np.full(2048, 100.0)
    counts[:THR] = 0.0   # правка Claude: порог прибора, как у боевого спектра
    res = sc.scale_from_repers(rep, counts, ECAL_REF, verbose=False)

    assert res["dropped"] == ["bad"], f"Ожидалось ['bad'], получено {res['dropped']}"
    assert res["n_repers"] == 4, f"Ожидалось 4 репера, получено {res['n_repers']}"
    assert np.all(np.isfinite(res["coeffs"])), "Коэффициенты содержат NaN/Inf"
    assert not res.get("rejected"), f"Результат отклонен: {res}"


def test_samplecal_three_repers_not_quadratic():
    """С тремя реперами должна строиться линейная шкала (2 коэффициента)."""
    rep = [
        {
            "name": "p%d" % i,
            "ch": (E - ECAL_LIN[0]) / ECAL_LIN[1],   # правка Claude: линейная истинная шкала
            "ch_err": 0.05,
            "E": E,
            "dE": 0.05,
            "source": sc.SRC_LIB,
            "nuclide": "X",
            "fwhm": 8.0,
            "signif": 50.0
        }
        for i, E in enumerate((351.932, 661.657, 1764.494))
    ]

    counts = np.full(2048, 100.0)
    counts[:THR] = 0.0   # правка Claude: порог прибора, как у боевого спектра
    res = sc.scale_from_repers(rep, counts, ECAL_LIN, verbose=False)

    assert len(res["coeffs"]) == 2, f"Ожидалось 2 коэффициента, получено {len(res['coeffs'])}"
    assert not res.get("rejected"), f"Результат отклонен: {res}"


def test_samplecal_xml_passes_becqmoni_check(tmp_path):
    """XML-вывод должен проходить проверку BecqMoni."""
    d = load_spectrum(BG)
    B = d["counts"].astype(float)
    n = len(B)

    y, centers = _peaks(n, LINES[:5], area=60000.0, a=2.0, b=1.03)
    S = _stretch(B, 2.0, 1.03, n) + y
    smp = _write_like_bg(tmp_path, "sample.xml", S)
    rep = _repers_by_click(S, centers, LINES[:5], shift=1.0)

    cal_res = sc.run_sample_cal(BG, smp, rep, verbose=False)
    out = pl.run_subtract(
        BG, smp, cal_res["h"], cal_res, cal_res["fwhm_points"],
        tmp_path / "z.xml", tmp_path / "s.xml",
        verbose=False
    )

    text = (tmp_path / "z.xml").read_text("utf-8")

    # Извлечение энергетических коэффициентов
    m = re.search(r"<EnergyCalibration>.*?<Coefficients>(.*?)</Coefficients>", text, re.S)
    assert m is not None, "Не найдена секция EnergyCalibration"
    ce = [float(v) for v in re.findall(r"<Coefficient>([^<]+)</Coefficient>", m.group(1))]

    # Извлечение коэффициентов ПШПВ
    m_fwhm = re.search(r"<SqrtFwhmCalibration>.*?<Coefficients>(.*?)</Coefficients>", text, re.S)
    assert m_fwhm is not None, "Не найдена секция SqrtFwhmCalibration"
    cf = [float(v) for v in re.findall(r"<Coefficient>([^<]+)</Coefficient>", m_fwhm.group(1))]

    # Проверки
    assert cal_res["mode"] == "visible", f"Ожидалось mode='visible', получено {cal_res['mode']}"
    assert not cal_res["rejected"], f"Результат отклонен: {cal_res}"
    
    ok_energy, msg_energy = lw.becqmoni_ok(ce, n)[:2]
    assert ok_energy is True, f"Проверка BecqMoni для энергии не пройдена: {msg_energy}"
    
    ok_fwhm, msg_fwhm = F.check(cf, n)[:2]   # правка Claude: check возвращает 4 значения
    assert ok_fwhm is True, f"Проверка ПШПВ не пройдена: {msg_fwhm}"
    
    # Проверка точности для линии Cs-137
    c_bg = float(poly_inverse(ECAL_REF, 661.657))
    E_rec = float(poly_eval(ce, c_bg))
    assert abs(E_rec - 661.657) <= 1.5, (
        f"Энергия Cs-137 не восстановлена: получено {E_rec}, ожидалось 661.657"
    )


def test_cancel_flag_stops_ecal_subtract_and_samplecal(tmp_path):
    """LABR-11 п.7: флаг отмены останавливает «Шкалу», «Вычесть» и «Шкалу по пробе» между этапами (правка Claude)."""
    from gamma.calibration.intrinsic_shape import robust as rb
    rb.CANCEL["flag"] = True
    try:
        with pytest.raises(rb.Cancelled):
            pl.run_ecal(BG, BG, [0.0, 1.0, 0.0], verbose=False)
        with pytest.raises(rb.Cancelled):
            pl.run_subtract(BG, BG, [0.0, 1.0, 0.0], ECAL_REF, [(29.2, 2.47), (847.1, 11.63)],
                            tmp_path / "z.xml", tmp_path / "s.xml", verbose=False)
        with pytest.raises(rb.Cancelled):
            sc.run_sample_cal(BG, BG, [], verbose=False)
    finally:
        rb.CANCEL["flag"] = False
    assert not (tmp_path / "z.xml").exists()


def test_samplecal_37_only_low_branch_and_wide_peak_not_in_fwhm(tmp_path):
    """Найдено живым прогоном (правка Claude): 37 кэВ -- только в ветвь ниже 100 кэВ, не в основной полином; пик с
    подозрительной ПШПВ -- не в точки ПШПВ; модель ПШПВ из точек строится (CheckCalibration)."""
    d = load_spectrum(BG)
    B = d["counts"].astype(float)
    n = len(B)
    y, centers = _peaks(n, LINES[:5], area=60000.0, a=2.0, b=1.03)
    smp = _write_like_bg(tmp_path, "sample.xml", _stretch(B, 2.0, 1.03, n) + y)
    rep = _repers_by_click(_stretch(B, 2.0, 1.03, n) + y, centers, LINES[:5], shift=1.0)
    rep[0] = dict(rep[0], warn=["ПШПВ 40.00 кан отличается от ожидаемой 27.00 больше чем в 1.5 раза"])
    r = sc.run_sample_cal(BG, smp, rep, verbose=False)
    assert "a37" in r["visible_sample"] and all(p["E"] >= 100.0 for p in r["points"])
    assert all(p["E"] >= 100.0 for p in r["scale_sample"]["points"])   # и в шкале пробы
    assert "low" in r and r["coeffs_xml"] == r["low"]["xml"]["coeffs"]
    ch0 = float(poly_eval(r["h"], rep[0]["ch"]))
    assert all(abs(p[0] - ch0) > 1.0 for p in r["fwhm_points"])
    assert F.check(F.fit_fwhm(n, [p[0] for p in r["fwhm_points"]], [p[1] for p in r["fwhm_points"]]), n)[0] is True


def test_samplecal_false_intrinsic_anchor_dropped():
    """Найдено живым прогоном (правка Claude): «пик 37 кэВ» сильной пробы на 82 кэВ по шкале её реперов отбрасывается,
    верная опора остаётся."""
    rep = [{"name": "p%d" % i, "ch": float(poly_inverse(ECAL_REF, E)), "ch_err": 0.05, "E": E, "dE": 0.05,
            "source": sc.SRC_LIB, "nuclide": "X", "fwhm": 8.0, "signif": 50.0} for i, E in enumerate((351.932, 661.657, 1120.287, 1764.494))]
    counts = np.full(2048, 100.0)
    counts[:THR] = 0.0
    anc = {"a37": {"x": float(poly_inverse(ECAL_REF, 82.0)), "visible": True},
           "edge": {"x": float(poly_inverse(ECAL_REF, 788.74)), "visible": True}, "sum": {"x": 0.0, "visible": False}}
    intr = [{"name": "La138_37_intr"}, {"name": "La138_789_edge"}]
    out = sc.check_intrinsic(anc, intr, rep, counts, ECAL_REF, verbose=False)
    assert [r["name"] for r in out] == ["La138_789_edge"] and anc["a37"]["visible"] is False and anc["edge"]["visible"] is True


def test_samplecal_bg_own_scale_corrects_file_calibration():
    """Правка Claude: шкала фона по своим структурам (2 точки, полином отвергнут) -- калибровка файла с поправкой
    E' = a + b*E_файла проходит через обе точки."""
    cal_b = [-5.0, 1.40, 1.0e-4]
    pts = [{"ch": 300.0, "E": 450.0}, {"ch": 900.0, "E": 1400.0}]
    B = np.full(2048, 100.0)
    B[:THR] = 0.0
    C, how = sc.bg_own_scale({"rejected": "порог", "points": pts}, cal_b, B)
    assert "поправкой" in how
    for q in pts:
        assert abs(float(poly_eval(C, q["ch"])) - q["E"]) < 1e-6
    narrow = [{"ch": 500.0, "E": 725.0}, {"ch": 850.0, "E": 1262.0}]   # охват меньше 3 раз (поправка порог прошла бы) -- файл
    C2, how2 = sc.bg_own_scale({"rejected": "порог", "points": narrow}, cal_b, B)
    assert how2 == "калибровка файла" and C2 == cal_b


def test_samplecal_energy_mode_sample_line_lands_at_its_energy(tmp_path):
    """Правка Claude: проба без видимых структур LaBr3 -- режим «через энергии»; линия пробы на сетке фона имеет по шкале
    XML свою энергию, отклонение квадратичного h ограничено."""
    n = 2048
    x = np.arange(n, dtype=float)
    y, centers = _peaks(n, LINES[:5], area=60000.0, a=2.0, b=1.03)
    S = np.random.default_rng(3).poisson(3000.0 * np.exp(-x / 300.0) + 50.0 + y).astype(float)
    S[:THR] = 0.0
    smp = _write_like_bg(tmp_path, "sample.xml", S)
    r = sc.run_sample_cal(BG, smp, _repers_by_click(S, centers, LINES[:5], shift=1.0), verbose=False)
    assert r["mode"] == "energy" and not r["rejected"]
    for c, (nuc, E) in zip(centers, LINES[:5]):
        assert abs(float(poly_eval(r["coeffs_xml"], poly_eval(r["h"], c))) - E) <= 1.5, (nuc, E)


def test_fit_click_degenerate_fit_rejected_weak_peak_kept():
    """Правка Claude: щелчок без пика (нули, ниже порога) -- центроид не определён, репер не создаётся; слабый, но
    подогнанный пик (значимость < 5) -- с предупреждением."""
    B = load_spectrum(BG)["counts"].astype(float)
    assert sc.fit_click(np.zeros(2048), 500.0, 9.0, 9.0) is None
    assert sc.fit_click(B, 5.0, sc.fwhm_guess(5.0, 2048)) is None
    assert sc.fit_click(B, 2.0, sc.fwhm_guess(2.0, 2048)) is None
    f = sc.fit_click(B, 300.0, sc.fwhm_guess(300.0, 2048))
    assert f is not None and f["signif"] < sc.SIGNIF_MIN and any("незначимый" in w for w in f["warn"])
