# tests/calibration/test_samplecal_net.py
# Тесты архитектуры «разность на шкале пробы» (samplecal_net): решение о шкале XML (охват, независимые
# линии, LOO), живые проверки на обезличенной фикстуре образца 8192 (Th/U) и фона 8192 — без пропуска, перевод фона на
# сетку пробы (видимые структуры и точный перевод через энергии), дисперсия разности, XML на сетке пробы, пороги
# переноса JSON. Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-10-11.

import re
import math
import pathlib
import numpy as np
import pytest
from scipy.ndimage import gaussian_filter1d
from gamma.calibration.intrinsic_shape.core import load_spectrum, poly_eval, poly_deriv, poly_inverse, rebin_conserve
from gamma.calibration.intrinsic_shape import samplecal as sc
from gamma.calibration.intrinsic_shape import samplecal_q as q
from gamma.calibration.intrinsic_shape import samplecal_net as sn
from gamma.calibration.intrinsic_shape import lowscale as lw

FIX = pathlib.Path(__file__).resolve().parents[1] / "fixtures" / "intrinsic_shape"
BG = FIX / "labr3_background.xml"
MED = FIX / "labr3_8192_sample_th_u.xml"
MED_BG = FIX / "labr3_8192_background.xml"
ECAL_REF = [-0.48496698, 1.35007262, 4.58971597e-04]
TRUE3 = [2.0, 1.30, 4.0e-4, 6.0e-8]
THR = 16
MED_CLICKS = [(628.3, "Pb-212", 238.632), (864.9, "Ac-228", 338.32), (898.9, "Pb-214", 351.932), (1431.4, "Tl-208", 583.187),
              (2121.8, "Ac-228", 911.0), (2236.5, "Ac-228", 967.96), (2530.2, "Bi-214", 1120.294), (3680.6, "Bi-214", 1764.491),
              (4978.3, "Tl-208", 2614.511)]
REPER_LINES = [("Pb-212", 238.632), ("Pb-214", 351.932), ("Tl-208", 583.187), ("Cs-137", 661.657), ("Ac-228", 911.204),
               ("Bi-214", 1120.294), ("Bi-214", 1764.491), ("Tl-208", 2614.511)]
IND_LINES = [("Pb-214", 295.224), ("Ac-228", 463.004), ("Bi-214", 1238.11), ("Bi-214", 2204.06)]


def fwhm_true_keV(E):
    return math.sqrt(4.0 + 0.9 * E + 2.0e-4 * E * E)


def _gauss(n, c, w, area):
    x = np.arange(n, dtype=float)
    s = w / 2.3548
    return area / (math.sqrt(2 * math.pi) * s) * np.exp(-0.5 * ((x - c) / s) ** 2)


def _lines_on(scale, lines, n, area=30000.0):
    y = np.zeros(n)
    cs = []
    for nuc, E in lines:
        c = float(poly_inverse(scale, E, hi=n + 200))
        cs.append(c)
        y += _gauss(n, c, fwhm_true_keV(E) / poly_deriv(scale, c), area)
    return y, cs


def _write_like(tmp_path, name, counts, cal=None):
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


def _med_repers(clicks=MED_CLICKS, shift_E=None):
    """Реперы фикстуры образца по кликам (как PeakPicker; при «мультиплет» -- двойной фит). Правка Claude: генерация
    подгоняла по массиву каналов вместо клика и искала «мультиплет» как элемент списка."""
    S = load_spectrum(MED)["counts"].astype(float)
    reps = []
    for x, nuc, E in clicks:
        fw = sc.fwhm_guess(x, len(S))
        f = sc.fit_click(S, x, fw, fwhm_expected=fw)
        if f and any("мультиплет" in w for w in f["warn"]):
            dd = q.fit_doublet(S, x, fw)
            f = min(dd, key=lambda c: abs(c["ch"] - x)) if dd else None
        if f:
            r = sc.reper_from_choice(f, E, nuc)
            if shift_E and abs(E - shift_E[0]) < 0.01:
                r = dict(r, E=r["E"] + shift_E[1])
            reps.append(r)
    return reps


def _flat_repers(Es, scale):
    return [
        {"name": "p%d" % i, "ch": float(poly_inverse(scale, E, hi=4000)), "ch_err": 0.05, "E": E, "dE": 0.05,
         "source": sc.SRC_LIB, "nuclide": "X", "fwhm": 8.0, "signif": 50.0, "warn": []}
        for i, E in enumerate(Es)
    ]


def test_medallion_decision_not_worse_than_file():
    d = load_spectrum(MED)
    S = d["counts"].astype(float)
    cal = list(d["cal"])
    dec = sn.decide_scale(S, cal, _med_repers(), verbose=False)
    assert len(dec["independent"]) >= sn.IND_MIN
    if dec["replaced"]:
        assert dec["rms_new"] < dec["rms_file"]
    else:
        assert dec["coeffs"] == cal
        assert dec["why"]
    assert sn.rms_on(dec["coeffs"], dec["independent"]) <= sn.rms_on(cal, dec["independent"]) + 1e-9
    # правка Claude: пересчёт СКО своей формулой (не sn.rms_on) и исключения ложных «независимых» линий
    own = lambda C: math.sqrt(np.mean([(np.polyval(C[::-1], l["ch"]) - l["E"]) ** 2 for l in dec["independent"]]))
    assert abs(dec["rms_file"] - own(cal)) < 1e-6 and abs(dec["rms_new"] - own(dec["coeffs_new"])) < 1e-6
    bad = (1460.82, 2118.51, 1588.2, 510.77)   # K-40 на сумм-пике La-138, вылеты Tl-208, аннигиляция
    assert not any(abs(l["E"] - e) < 1.0 for l in dec["independent"] for e in bad)


def test_medallion_narrow_set_forbidden():
    d = load_spectrum(MED)
    S = d["counts"].astype(float)
    cal = list(d["cal"])
    clicks = [c for c in MED_CLICKS if c[2] in (351.932, 583.187, 911.0, 967.96)]
    dec = sn.decide_scale(S, cal, _med_repers(clicks), verbose=False)
    assert not dec["replaced"] and not dec["recommend"]
    assert "охват" in dec["why"] or "трети" in dec["why"]
    assert dec["coeffs"] == cal


def test_medallion_extreme_reper_error_rejected():
    d = load_spectrum(MED)
    S = d["counts"].astype(float)
    cal = list(d["cal"])
    rep = _med_repers(shift_E=(2614.511, 25.0))
    dec = sn.decide_scale(S, cal, rep, verbose=False)
    assert not dec["replaced"] and not dec["recommend"] and "LOO" in dec["why"]   # правка Claude (VERIFY7)
    assert dec["coeffs"] == cal


def test_synthetic_cubic_replaced_and_checked():
    n = 2048
    x = np.arange(n, dtype=float)
    y, cs = _lines_on(TRUE3, REPER_LINES, n)
    yi, ci = _lines_on(TRUE3, IND_LINES, n)
    S = np.random.default_rng(11).poisson(400.0 * np.exp(-x / 400.0) + 20.0 + y + yi).astype(float)
    S[:THR] = 0
    file_cal = [TRUE3[0] + 3.0, TRUE3[1] * 1.003, TRUE3[2], TRUE3[3]]
    rep = [sc.reper_from_choice(sc.fit_click(S, c + 1.0, sc.fwhm_guess(c, n)), E, nuc) for c, (nuc, E) in zip(cs, REPER_LINES)]
    dec = sn.decide_scale(S, file_cal, rep, verbose=False)
    assert dec["recommend"] and not dec["replaced"] and dec["coeffs"] == file_cal   # правка Claude: применяет пользователь
    dec = sn.apply_reper_scale(dec)
    assert dec["replaced"] and dec["coeffs"] == dec["coeffs_new"]
    assert len(dec["independent"]) >= 3
    assert dec["rms_new"] < 1.0
    assert dec["rms_new"] < dec["rms_file"]
    own = lambda C: math.sqrt(np.mean([(np.polyval(C[::-1], l["ch"]) - l["E"]) ** 2 for l in dec["independent"]]))
    assert abs(dec["rms_file"] - own(file_cal)) < 1e-6 and abs(dec["rms_new"] - own(dec["coeffs"])) < 1e-6   # правка Claude



def test_synthetic_extreme_reper_error_rejected_by_loo():
    """Правка Claude: крайний репер (2614.5) +25 кэВ, независимые линии только до 1238 кэВ -- кросс-проверка его не
    видит, отвергнуть обязана LOO по всем реперам (крайние включены)."""
    n = 2048
    x = np.arange(n, dtype=float)
    y, cs = _lines_on(TRUE3, REPER_LINES, n)
    yi, ci = _lines_on(TRUE3, IND_LINES[:3], n)
    S = np.random.default_rng(11).poisson(400.0 * np.exp(-x / 400.0) + 20.0 + y + yi).astype(float)
    S[:THR] = 0
    file_cal = [TRUE3[0] + 3.0, TRUE3[1] * 1.003, TRUE3[2], TRUE3[3]]
    rep = [sc.reper_from_choice(sc.fit_click(S, c + 1.0, sc.fwhm_guess(c, n)), E, nuc) for c, (nuc, E) in zip(cs, REPER_LINES)]
    rep[-1] = dict(rep[-1], E=rep[-1]["E"] + 25.0)
    dec = sn.decide_scale(S, file_cal, rep, verbose=False)
    assert not dec["replaced"]
    assert "LOO" in dec["why"]
    assert dec["coeffs"] == file_cal

def test_no_independent_lines_warning_not_forbidden():
    n = 2048
    x = np.arange(n, dtype=float)
    y, cs = _lines_on(TRUE3, REPER_LINES, n)
    S = np.random.default_rng(11).poisson(400.0 * np.exp(-x / 400.0) + 20.0 + y).astype(float)
    S[:THR] = 0
    file_cal = [TRUE3[0] + 3.0, TRUE3[1] * 1.003, TRUE3[2], TRUE3[3]]
    rep = [sc.reper_from_choice(sc.fit_click(S, c + 1.0, sc.fwhm_guess(c, n)), E, nuc) for c, (nuc, E) in zip(cs, REPER_LINES)]
    dec = sn.decide_scale(S, file_cal, rep, verbose=False)
    assert not dec["replaced"]   # правка Claude (VERIFY7): нет независимых линий -- предупреждение, не запрет
    assert any("доп. проверки нет" in w for w in dec["warnings"])
    assert dec["recommend"]


def test_bg_to_sample_visible_cancels_structures(tmp_path):
    B = load_spectrum(BG)["counts"].astype(float)
    n = len(B)
    S0 = rebin_conserve(B, 2.0 + 1.03 * (np.arange(n + 1) - 0.5), n)
    c = 2.0 + 1.03 * float(poly_inverse(ECAL_REF, 661.657))
    pk = np.random.default_rng(2).poisson(_gauss(n, c, 1.03 * sc.fwhm_guess(c / 1.03, n), 30000.0)).astype(float)
    smp = _write_like(tmp_path, "s.xml", S0 + pk)
    r = sn.run_sample_cal(BG, smp, [], verbose=False)
    assert r["mode"] == "visible"
    x = sn.run_sample_subtract(BG, smp, r, r["fwhm_points"], tmp_path / "z.xml", tmp_path / "s2.xml", verbose=False)
    res, var = x["resid"], x["var"]
    
    # Окно вокруг сумм-пика
    c_sum = 2.0 + 1.03 * 846.9
    a = max(0, int(c_sum) - 25)
    b = min(n, int(c_sum) + 25)
    assert abs(res[a:b].sum()) <= 3 * math.sqrt(var[a:b].sum())
    
    # Окно вокруг края
    c_edge = 2.0 + 1.03 * 497.6
    a = max(0, int(c_edge) - 25)
    b = min(n, int(c_edge) + 25)
    assert abs(res[a:b].sum()) <= 3 * math.sqrt(var[a:b].sum())
    
    # Пик: окно c ± 2.5 * ПШПВ
    fwhm_pts = 1.03 * sc.fwhm_guess(c / 1.03, n)
    half_w = 2.5 * fwhm_pts
    a = max(0, int(c - half_w))
    b = min(n, int(c + half_w))
    assert abs(res[a:b].sum() - pk.sum()) <= 3 * math.sqrt(var[a:b].sum())
    # правка Claude: вне пика фон, переведённый на пробу, совпадает с пробой поканально (точный перевод;
    # проба построена тем же переводом без шума; порог прибора до кан. 40 исключён) -- chi2/ndf << 1; ошибка формы перевода его поднимает
    m = np.ones(n, bool); m[a:b] = False; m[:40] = False; m[n - 300:] = False
    assert np.sum(res[m] ** 2 / var[m]) / m.sum() < 0.05


def test_bg_to_sample_energy_exact(tmp_path):
    TB = ECAL_REF
    bgp = _write_like(tmp_path, "b.xml", load_spectrum(BG)["counts"].astype(float), cal=TB)
    TS = [1.0, 1.28, 4.2e-4, 5.0e-8]
    B = load_spectrum(bgp)["counts"].astype(float)
    n = len(B)
    # правка Claude: проба строится независимо от проверяемого кода -- края каналов фона в каналах пробы точным
    # обращением шкалы пробы (poly_inverse по каждому краю), не через sn.bg_on_sample
    xb = np.arange(n + 1) - 0.5
    S0 = rebin_conserve(B, np.array([float(poly_inverse(TS, poly_eval(TB, v), hi=n + 400)) for v in xb]), n)
    c = float(poly_inverse(TS, 661.657, hi=n + 200))
    pk = np.random.default_rng(2).poisson(_gauss(n, c, fwhm_true_keV(661.657) / poly_deriv(TS, c), 30000.0)).astype(float)
    smp = _write_like(tmp_path, "s.xml", S0 + pk, cal=TS)
    r = sn.run_sample_cal(bgp, smp, [], mode="energy", verbose=False)
    assert r["mode"] == "energy"
    x = sn.run_sample_subtract(bgp, smp, r, r["fwhm_points"], tmp_path / "z.xml", tmp_path / "s2.xml", verbose=False)
    res, var = x["resid"], x["var"]
    
    # Каналы пробы для сумм-пика и края
    c_sum_s = float(poly_inverse(TS, poly_eval(TB, 846.9), hi=n + 200))
    c_edge_s = float(poly_inverse(TS, poly_eval(TB, 497.6), hi=n + 200))
    
    # Окно вокруг сумм-пика
    a = max(0, int(c_sum_s) - 25)
    b = min(n, int(c_sum_s) + 25)
    assert abs(res[a:b].sum()) <= 3 * math.sqrt(var[a:b].sum())
    
    # Окно вокруг края
    a = max(0, int(c_edge_s) - 25)
    b = min(n, int(c_edge_s) + 25)
    assert abs(res[a:b].sum()) <= 3 * math.sqrt(var[a:b].sum())
    
    # Пик: окно c ± 2.5 * ПШПВ
    fwhm_pts = fwhm_true_keV(661.657) / poly_deriv(TS, c)
    half_w = 2.5 * fwhm_pts
    a = max(0, int(c - half_w))
    b = min(n, int(c + half_w))
    assert abs(res[a:b].sum() - pk.sum()) <= 3 * math.sqrt(var[a:b].sum())
    # правка Claude: вне пика фон, переведённый на пробу, совпадает с пробой поканально (точный перевод;
    # проба построена тем же переводом без шума; порог прибора до кан. 40 исключён) -- chi2/ndf << 1; ошибка формы перевода его поднимает
    m = np.ones(n, bool); m[a:b] = False; m[:40] = False; m[n - 300:] = False
    assert np.sum(res[m] ** 2 / var[m]) / m.sum() < 0.05
