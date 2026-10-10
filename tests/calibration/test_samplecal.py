"""Тесты калибровки по пробе (samplecal, samplecal_q): нелинейная истинная шкала и независимые ширины, пулы
центроидов, перенос набора реперов, источник энергии и дублеты, отказы fit_click, двойной фит, ранжирование линий,
собственная шкала фона, отмена. Решение о шкале XML и разность на сетке пробы -- test_samplecal_net*.py (живые
проверки на обезличенной фикстуре, без пропуска). Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-10-10;
правка Claude 2026-10-11 -- тесты удалённого кода сняты."""

import os
import re
import math
import pathlib
import numpy as np
import pytest
from gamma.calibration.intrinsic_shape.core import load_spectrum, poly_eval, poly_deriv, poly_inverse, rebin_conserve
from gamma.calibration.intrinsic_shape import samplecal as sc
from gamma.calibration.intrinsic_shape import samplecal_q as q
from gamma.calibration.intrinsic_shape import pipeline as pl
from gamma.calibration.intrinsic_shape import lowscale as lw
from gamma.calibration.intrinsic_shape import fwhm as F
from gamma.calibration.intrinsic_shape import robust as rb
from gamma.calibration.intrinsic_shape import samplecal_net as sn

FIX = pathlib.Path(__file__).resolve().parents[1] / "fixtures" / "intrinsic_shape"
BG = FIX / "labr3_background.xml"
ECAL_REF = [-0.48496698, 1.35007262, 4.58971597e-04]
TRUE3 = [2.0, 1.30, 4.0e-4, 6.0e-8]          # нелинейная (кубическая) истинная шкала синтетики
THR = 16                                       # первый ненулевой канал фона фикстуры (порог прибора)
LINES = [("Pb-214", 351.932), ("Bi-214", 609.312), ("Cs-137", 661.657), ("Bi-214", 1120.294), ("Bi-214", 1764.491), ("Tl-208", 2614.511)]
LINES8 = [("Pb-212", 238.632), ("Pb-214", 351.932), ("Tl-208", 583.187), ("Cs-137", 661.657), ("Ac-228", 911.204),
          ("Bi-214", 1120.294), ("Bi-214", 1764.491), ("Tl-208", 2614.511)]
MED_CLICKS = [(628.3, "Pb-212", 238.632), (864.9, "Ac-228", 338.32), (898.9, "Pb-214", 351.932), (1431.4, "Tl-208", 583.187),
              (2121.8, "Ac-228", 911.0), (2236.5, "Ac-228", 967.96), (2530.2, "Bi-214", 1120.294), (3680.6, "Bi-214", 1764.491),
              (4978.3, "Tl-208", 2614.511)]


def fwhm_true_keV(E):
    return math.sqrt(4.0 + 0.9 * E + 2.0e-4 * E * E)


def _stretch(counts, a, b, n):
    edges = a + b * (np.arange(len(counts) + 1) - 0.5)
    return rebin_conserve(np.asarray(counts, float), edges, n)


def _peaks_true(n, lines, scale, area=30000.0):
    x = np.arange(n, dtype=float)
    y = np.zeros(n)
    centers = []
    for nuc, E in lines:
        c = float(poly_inverse(scale, E, hi=n + 200))
        w = fwhm_true_keV(E) / float(poly_deriv(scale, c))
        s = w / 2.3548
        y += area / (math.sqrt(2 * math.pi) * s) * np.exp(-0.5 * ((x - c) / s) ** 2)
        centers.append(c)
    return y, centers


def _synth(seed, scale=TRUE3, lines=LINES8, n=2048):
    x = np.arange(n, dtype=float)
    y, cs = _peaks_true(n, lines, scale)
    S = np.random.default_rng(seed).poisson(400.0 * np.exp(-x / 400.0) + 20.0 + y).astype(float)
    S[:THR] = 0.0
    return S, cs


def _peaks(n, lines, area=30000.0, a=0.0, b=1.0):
    x = np.arange(n, dtype=float)
    y = np.zeros(n)
    centers = []
    for nuc, E in lines:
        c = a + b * poly_inverse(ECAL_REF, E)
        w = b * sc.fwhm_guess(c / b, n)
        s = w / 2.3548
        y += area / (math.sqrt(2 * math.pi) * s) * np.exp(-0.5 * ((x - c) / s) ** 2)
        centers.append(c)
    return y, centers


def _sigma_E(scale, r):
    return math.hypot(r["dE"], float(poly_deriv(scale, r["ch"])) * r["ch_err"])


def _repers_by_click(counts, centers, lines, shift=2.0):
    reps = []
    for c, (nuc, E) in zip(centers, lines):
        f = sc.fit_click(counts, c + shift, sc.fwhm_guess(c, len(counts)))
        assert f is not None, (nuc, E)
        reps.append(sc.reper_from_choice(f, E, nuc))
    return reps


def _write_like(tmp_path, name, counts, cal=None):
    """Спектр в XML по образцу фона фикстуры; при cal -- своя энергетическая калибровка (правка Claude: ответ модели
    не заменял DataPoint)."""
    text = BG.read_text("utf-8")
    vals = iter(np.round(np.maximum(counts, 0)).astype(int))
    text = re.sub(r"<DataPoint>[^<]*</DataPoint>", lambda m: "<DataPoint>%d</DataPoint>" % next(vals), text)
    if cal is not None:
        coeff = "".join("<Coefficient>%r</Coefficient>" % float(c) for c in cal)
        text = re.sub(r"(<EnergyCalibration>.*?<PolynomialOrder>)\d+(</PolynomialOrder>\s*<Coefficients>).*?(</Coefficients>)",
                      lambda m: m.group(1) + str(len(cal) - 1) + m.group(2) + coeff + m.group(3), text, count=1, flags=re.S)
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


def _flat_repers(Es, scale):
    reps = []
    for i, E in enumerate(Es):
        reps.append({
            "name": "p%d" % i,
            "ch": float(poly_inverse(scale, E, hi=4000)),
            "ch_err": 0.05,
            "E": E,
            "dE": 0.05,
            "source": sc.SRC_LIB,
            "nuclide": "X",
            "fwhm": 8.0,
            "signif": 50.0,
            "warn": []
        })
    return reps


def _flat(n=2048):
    c = np.full(n, 100.0)
    c[:THR] = 0.0
    return c


def test_cubic_true_scale_recovered_with_degree_selection():
    S, cs = _synth(7)
    rep = _repers_by_click(S, cs, LINES8)
    res = sc.scale_from_repers(rep, S, [0.0, 1.6, 0.0], verbose=False)
    assert not res.get("rejected")
    assert res["deg"] >= 3
    for r in rep:
        E_rec = poly_eval(res["coeffs"], r["ch"])
        assert abs(E_rec - r["E"]) <= 3 * _sigma_E(res["coeffs"], r)
    g = np.linspace(poly_inverse(TRUE3, 150.0, hi=2300), poly_inverse(TRUE3, 2200.0, hi=2300), 300)
    assert np.max(np.abs(poly_eval(res["coeffs"], g) - poly_eval(TRUE3, g))) <= 1.0


def test_centroid_pulls_unit_width():
    z = []
    for seed in range(100, 106):
        S, cs = _synth(seed)
        for c in cs:
            f = sc.fit_click(S, c + 1.0, sc.fwhm_guess(c, len(S)))
            z.append((f["ch"] - c) / f["ch_err"])
    std_z = np.std(z)
    assert 0.6 <= std_z <= 1.6, f"std(z) = {std_z}"


def test_json_transfer_same_scale_and_physical_gain(tmp_path):
    d = load_spectrum(BG)
    B = d["counts"].astype(float)
    cal = d["cal"]
    n = len(B)
    y, centers = _peaks(n, LINES[:5])
    s1 = B + y
    rep1 = _repers_by_click(s1, centers, LINES[:5])
    intr, anc1 = sc.intrinsic_repers(s1, cal)
    rep1 = rep1 + intr
    p = tmp_path / "r.json"
    sc.save_json(p, rep1, "LaBr3 test", n, anc1, cal=cal)
    data = sc.load_json(p)
    assert data["cal"] == [float(c) for c in cal]
    for r in data["repers"]:
        assert "signif" in r
        assert "warn" in r

    s2 = _stretch(s1, 3.0, 1.03, n)
    rep2, m, dropped = sc.apply_to(data, s2, cal, verbose=False)
    sc1 = sc.scale_from_repers(rep1, s1, [0.0, 1.6, 0.0], verbose=False)
    sc2 = sc.scale_from_repers(rep2, s2, [0.0, 1.6, 0.0], verbose=False)
    assert not sc1.get("rejected") and not sc2.get("rejected")
    assert abs(m[1] - 1.03) < 0.01 and abs(m[0] - 3.0) < 1.5
    non_labr = [r for r in rep2 if r["source"] != sc.SRC_LABR]
    for r in non_labr:
        y_val = r["ch"]
        x_val = (y_val - 3.0) / 1.03
        assert abs(poly_eval(sc2["coeffs"], y_val) - poly_eval(sc1["coeffs"], x_val)) <= 3 * _sigma_E(sc2["coeffs"], r)
    assert len(non_labr) == 5
    assert dropped == []

    s3 = _stretch(s1, 0.0, 1.5, n)
    with pytest.raises(ValueError, match="нефизичное усиление|опор"):
        sc.apply_to(data, s3, cal, verbose=False)


def test_manual_library_and_doublet_energy():
    fit = {"ch": 500.0, "ch_err": 0.1, "fwhm": 9.0, "signif": 50.0}
    assert sc.reper_from_choice(fit, 661.657, "Cs-137")["source"] == sc.SRC_LIB
    assert sc.reper_from_choice(fit, 700.0, "Cs-137")["source"] == sc.SRC_MANUAL
    assert sc.reper_from_choice(fit, 662.5, "Cs-137")["source"] == sc.SRC_MANUAL
    assert sc.reper_from_choice(fit, 1234.5)["source"] == sc.SRC_MANUAL
    c = sc.candidate_lines(968.0)[0]
    assert c["nuclide"] == "Ac-228"
    assert abs(c["E"] - (964.77 * 4.99 + 968.97 * 15.8) / (4.99 + 15.8)) < 0.05
    assert "дублет" in c["note"]
    r = sc.reper_from_choice(fit, c["E"], "Ac-228")
    assert r["source"] == sc.SRC_LIB
    assert r["dE"] >= 2.0


def test_nan_dropped_and_three_repers_linear():
    bad = dict(_flat_repers([609.312], ECAL_REF)[0], name="bad", ch=float("nan"), ch_err=float("nan"))
    res = sc.scale_from_repers(_flat_repers([351.932, 661.657, 1120.287, 1764.494], ECAL_REF) + [bad], _flat(), [6.0, 1.35007262, 4.58971597e-04], verbose=False)
    assert res["dropped"] == ["bad"]
    assert res["n_repers"] == 4
    assert all(np.isfinite(c) for c in res["coeffs"])
    assert not res.get("rejected")

    lin = [-0.48, 1.62]
    res3 = sc.scale_from_repers(_flat_repers([351.932, 661.657, 1764.494], lin), _flat(), [6.0, 1.62], verbose=False)
    assert len(res3["coeffs"]) == 2
    assert not res3.get("rejected")


def test_intrinsic_energies_are_nuclear_data():
    en = sc.intrinsic_energies()
    assert abs(en["a37"][0] - 35.75) < 1e-6
    assert abs(en["edge"][0] - 788.74) < 1e-6
    assert abs(en["sum"][0] - 1471.006) < 0.01


def test_fit_click_refusals_multiplet_and_doublet():
    B = load_spectrum(BG)["counts"].astype(float)
    assert sc.fit_click(np.zeros(2048), 500.0, 9.0, 9.0) is None
    assert sc.fit_click(B, 5.0, sc.fwhm_guess(5.0, 2048)) is None
    
    x = np.arange(2048, dtype=float)
    base = np.full(2048, 400.0)
    def pk(A):
        return A / (math.sqrt(2 * math.pi) * 4.0) * np.exp(-0.5 * ((x - 1000.0) / 4.0) ** 2)
    
    # значимость ~2.5 (< 3) -- отказ; при площади 500 -- ~3.6: принят с предупреждением (замер, круг 1)
    f2 = sc.fit_click(base + pk(340.0), 1000.0, 9.42, 9.42)
    assert f2 is None
    
    f4 = sc.fit_click(base + pk(500.0), 1000.0, 9.42, 9.42)
    assert f4 is not None
    assert q.SIGNIF_REFUSE <= f4["signif"] < sc.SIGNIF_MIN
    assert any("незначимый" in t for t in f4["warn"])
    
    w = 20.0
    s = w / 2.3548
    yd = 50.0 + 20000.0 / (math.sqrt(2 * math.pi) * s) * (np.exp(-0.5 * ((x - 990.0) / s) ** 2) + np.exp(-0.5 * ((x - 1010.0) / s) ** 2))
    Sd = np.random.default_rng(1).poisson(yd).astype(float)
    f = sc.fit_click(Sd, 1000.0, w, fwhm_expected=w)
    assert any("мультиплет" in t for t in f["warn"])
    d = q.fit_doublet(Sd, 1000.0, w)
    assert len(d) == 2
    assert abs(d[0]["ch"] - 990.0) <= 0.15 * w
    assert abs(d[1]["ch"] - 1010.0) <= 0.15 * w


def test_candidate_ranking_foreign_last():
    c = sc.candidate_lines(346.0, visible_series={"U"})
    assert c[0]["nuclide"] == "Pb-214"
    assert abs(c[0]["E"] - 351.932) < 0.01
    foreign_indices = [i for i, x in enumerate(c) if x["foreign"]]
    non_foreign_indices = [i for i, x in enumerate(c) if not x["foreign"]]
    if foreign_indices and non_foreign_indices:
        assert min(foreign_indices) > max(non_foreign_indices)
    for i in foreign_indices:
        assert "чужой" in c[i]["note"]
    c30 = sc.candidate_lines(346.0, visible_series={"U"}, n_max=200)
    # интенсивность решает: ближайшая к 664 кэВ -- слабая Bi-214 665.45 (1.5 %), первой должна быть Cs-137 661.66 (85 %)
    assert sc.candidate_lines(664.0)[0]["nuclide"] == "Cs-137"
    assert any(x["nuclide"] == "Eu-152" for x in c30)


def test_anchors_of_degenerate_spectrum_no_crash():
    a = sc.anchors_of(np.zeros(8192), [0.0, 0.36, 3.0e-5])
    assert all(a[k]["visible"] is False for k in ("a37", "edge", "sum"))


def test_bg_own_scale_corrects_file_calibration():
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



def test_false_intrinsic_anchor_dropped():
    """Найдено живым прогоном (правка Claude): «пик 37 кэВ» сильной пробы на 82 кэВ по шкале её реперов отбрасывается,
    верная опора остаётся."""
    rep = [{"name": "p%d" % i, "ch": float(poly_inverse(ECAL_REF, E)), "ch_err": 0.05, "E": E, "dE": 0.05,
            "source": sc.SRC_LIB, "nuclide": "X", "fwhm": 8.0, "signif": 50.0} for i, E in enumerate((351.932, 661.657, 1120.287, 1764.494))]
    counts = np.full(2048, 100.0)
    counts[:THR] = 0.0
    anc = {"a37": {"x": float(poly_inverse(ECAL_REF, 82.0)), "visible": True},
           "edge": {"x": float(poly_inverse(ECAL_REF, 788.74)), "visible": True}, "sum": {"x": 0.0, "visible": False}}
    intr = [{"name": "La138_37_intr"}, {"name": "La138_789_edge"}]
    out = sc.check_intrinsic(anc, intr, rep, counts, [6.0, 1.35007262, 4.58971597e-04], verbose=False)
    assert [r["name"] for r in out] == ["La138_789_edge"] and anc["a37"]["visible"] is False and anc["edge"]["visible"] is True



def test_cancel_flag_stops_stages(tmp_path):
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
            sn.run_sample_cal(BG, BG, [], verbose=False)
    finally:
        rb.CANCEL["flag"] = False
    assert not (tmp_path / "z.xml").exists()



