# tests/calibration/test_samplecal_net2.py
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
    centers = []
    for nuc, E in lines:
        c = float(poly_inverse(scale, E, hi=n + 200))
        centers.append(c)
        y += _gauss(n, c, fwhm_true_keV(E) / poly_deriv(scale, c), area)
    return y, centers


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


def test_difference_variance_unit_pulls(tmp_path):
    mu_b = gaussian_filter1d(load_spectrum(BG)["counts"].astype(float), 3.0)
    n = len(mu_b)
    rng = np.random.default_rng(4)
    B = rng.poisson(mu_b).astype(float)
    bgp = _write_like(tmp_path, "b.xml", B, cal=ECAL_REF)
    TS = [1.0, 1.28, 4.2e-4, 5.0e-8]
    mu_s, _ = sn.bg_on_sample(mu_b, n, {"mode": "energy", "cal_b": ECAL_REF, "scale_s": TS, "n_s": n})
    smp = _write_like(tmp_path, "s.xml", rng.poisson(mu_s).astype(float), cal=TS)
    r = sn.run_sample_cal(bgp, smp, [], mode="energy", verbose=False)
    x = sn.run_sample_subtract(bgp, smp, r, r["fwhm_points"], tmp_path / "z.xml", tmp_path / "s2.xml", verbose=False)
    m = mu_s > 20
    pull = x["resid"][m] / np.sqrt(x["var"][m])
    assert 0.8 <= np.std(pull) <= 1.2


def test_xml_on_sample_grid(tmp_path):
    r = sn.run_sample_cal(MED_BG, MED, [], verbose=False)
    sn.run_sample_subtract(MED_BG, MED, r, r["fwhm_points"], tmp_path / "z.xml", tmp_path / "s.xml", verbose=False)
    z = load_spectrum(tmp_path / "z.xml")
    assert len(z["counts"]) == 8192
    assert list(z["cal"]) == list(load_spectrum(MED)["cal"])
    assert lw.becqmoni_ok(list(z["cal"]), 8192)[0] is True
    assert z["live"] == load_spectrum(MED)["live"]
    sn.run_sample_subtract(MED_BG, MED, r, r["fwhm_points"], tmp_path / "z2.xml", tmp_path / "s2.xml", verbose=False, lt_out="bg")
    assert load_spectrum(tmp_path / "z2.xml")["live"] == load_spectrum(MED_BG)["live"]


def test_transfer_thresholds(tmp_path):
    d = load_spectrum(BG)
    B = d["counts"].astype(float)
    cal = list(d["cal"])
    n = len(B)
    
    strong_lines = [("Pb-214", 351.932), ("Cs-137", 661.657), ("Bi-214", 1120.294)]
    weak_line = ("Bi-214", 1764.491)
    
    y = np.zeros(n)
    repers = []
    
    for nuc, E in strong_lines:
        c = float(poly_inverse(cal, E, hi=n + 200))
        fw = sc.fwhm_guess(c, n)
        y += _gauss(n, c, fw, 40000.0)
        f = sc.fit_click(B + y, c, fw)
        if f:
            repers.append(sc.reper_from_choice(f, E, nuc))
            
    c_w = float(poly_inverse(cal, weak_line[1], hi=n + 200))
    fw_w = sc.fwhm_guess(c_w, n)
    y += _gauss(n, c_w, fw_w, 1500.0)
    f_w = sc.fit_click(B + y, c_w, fw_w)
    if f_w:
        repers.append(sc.reper_from_choice(f_w, weak_line[1], weak_line[0]))
        
    # правка Claude: репер 2000 кэВ сильный в исходном спектре и слабый (значимость < 5) в новом -- порог переноса
    c_q = float(poly_inverse(cal, 2000.0, hi=n + 200))
    fw_q = sc.fwhm_guess(c_q, n)
    f_q = sc.fit_click(B + y + _gauss(n, c_q, fw_q, 3000.0), c_q, fw_q)
    repers.append(sc.reper_from_choice(f_q, 2000.0, "Tl-208"))
    s1 = B + y + _gauss(n, c_q, fw_q, 3000.0)
    intr, anc = sc.intrinsic_repers(s1, cal)
    sc.save_json(tmp_path / "r.json", repers + intr, "t", n, anc, cal=cal)
    data = sc.load_json(tmp_path / "r.json")
    
    # Shift Cs-137 reper
    for r in data["repers"]:
        if "Cs137" in r["name"] or (r.get("nuclide") == "Cs-137" and abs(r["E"] - 661.657) < 1.0):
            # Find fwhm for this reper to shift
            # We need to estimate fwhm or use stored one if available
            # Spec says: "сдвинуть ch на +0.8 * fwhm"
            # Assuming r has 'fwhm' or we can estimate it. Let's assume r['fwhm'] exists or use a default.
            # In _flat_repers we set fwhm=8.0. In reper_from_choice it should be set.
            fwhm_val = r.get("fwhm", 8.0)
            r["ch"] += 0.8 * fwhm_val
            break
            
    # правка Claude: смещение 10 кан (перенос обязан учитывать c0), линия 2000 кэВ в новом спектре -- площадь 160 (значимость подгонки ~3.9, ниже порога 5)
    s2 = rebin_conserve(B + y, 10.0 + 1.03 * (np.arange(n + 1) - 0.5), n) + _gauss(n, 10.0 + 1.03 * c_q, 1.03 * fw_q, 160.0)
    new, m, dropped = sc.apply_to(data, s2, cal, verbose=False)
    assert abs(m[1] - 1.03) < 0.01
    names = [r["name"] for r in new]
    assert all(r["name"] in names for r in (repers[0], repers[2]))   # Cs-137 сдвинут выше и отброшен по замыслу
    assert any(abs(r["E"] - 2000.0) < 0.5 for r in data["repers"])
    assert not any(abs(r["E"] - 2000.0) < 0.5 for r in new)
    assert "Cs137_662" not in [r["name"] for r in new]
    assert any("Cs137_662" in t for t in dropped)
    assert all(r.get("signif", 99.0) >= sc.SIGNIF_MIN for r in new if r["source"] != sc.SRC_LABR)
    
    s3 = rebin_conserve(s1, 1.08 * (np.arange(n + 1) - 0.5), n)
    with pytest.raises(ValueError, match="нефизичное усиление|опор"):
        sc.apply_to(data, s3, cal, verbose=False)


def test_fit_click_1764_matches_independent_fit():
    S = load_spectrum(MED)["counts"].astype(float)
    f = sc.fit_click(S, 3680.6, sc.fwhm_guess(3680.6, len(S)))
    assert abs(f["ch"] - 3680.55) <= 2.0


def test_medallion_sum_anchor_not_visible():
    d = load_spectrum(MED)
    a = sc.anchors_of(d["counts"].astype(float), list(d["cal"]))
    assert a["sum"]["visible"] is False


def test_fwhm_points_exclude_weak_reper(tmp_path):
    """Правка Claude: по спецификации -- синтетика как в тесте 4, реперы подгонкой клика (генерация задавала ПШПВ 8 кан)."""
    n = 2048
    x = np.arange(n, dtype=float)
    y, cs = _lines_on(TRUE3, REPER_LINES, n)
    S = np.random.default_rng(11).poisson(400.0 * np.exp(-x / 400.0) + 20.0 + y).astype(float)
    S[:THR] = 0
    rep = [sc.reper_from_choice(sc.fit_click(S, c + 1.0, sc.fwhm_guess(c, n)), E, nuc) for c, (nuc, E) in zip(cs, REPER_LINES)]
    rep[2] = dict(rep[2], signif=4.0)
    smp = _write_like(tmp_path, "s.xml", S, cal=TRUE3)
    r = sn.run_sample_cal(BG, smp, rep, verbose=False)
    assert all(abs(p[0] - rep[2]["ch"]) > 1.0 for p in r["fwhm_points"])
    assert any(abs(p[0] - rep[3]["ch"]) < 1.0 for p in r["fwhm_points"])



def test_degree_tie_lower():
    lin = [-0.48, 1.62]
    rep = _flat_repers([238.6, 351.9, 583.2, 911.2, 1120.3, 1764.5], lin)
    c = np.full(2048, 100.0)
    c[:THR] = 0
    res = sc.scale_from_repers(rep, c, lin, verbose=False)
    assert res["deg"] == 1
