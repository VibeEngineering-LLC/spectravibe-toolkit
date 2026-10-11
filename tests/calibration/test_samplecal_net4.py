"""Тесты раунда 2 гейта «реперы пользователя» (VERIFY8): уникальность реперов, охват по третям и размаху, порог LOO,
вето независимых линий, предупреждение о систематической ошибке, шкала фона (экстраполяция в обеих ветвях, квадратичная
поправка, порог значимости доп. опор), коды CLI. Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-10-11."""

import numpy as np
import pytest
from gamma.calibration.intrinsic_shape.core import load_spectrum, poly_inverse
from gamma.calibration.intrinsic_shape import samplecal as sc
from gamma.calibration.intrinsic_shape import samplecal_net as sn
from gamma.calibration.intrinsic_shape import cli
from test_samplecal_net import _med_repers, _flat_repers, MED, MED_BG, BG

def test_duplicate_repers_dropped():
    d = load_spectrum(MED)
    S = d["counts"].astype(float)
    cal = list(d["cal"])
    three = [r for r in _med_repers() if round(r["E"]) in (239, 911, 2615)]
    assert len(three) == 3
    dups = three + [dict(r, name=r["name"] + "_2") for r in three]
    good, dropped = sc.valid_repers(dups)
    assert len(good) == 3 and len(dropped) == 3 and all("(повтор" in x for x in dropped)
    same_E = three + [dict(three[1], name="other", ch=three[1]["ch"] + 300.0)]
    assert len(sc.valid_repers(same_E)[0]) == 3
    dec = sn.decide_scale(S, [cal[0] + 20.0] + list(cal[1:]), dups, verbose=False)
    assert not dec["recommend"]
    assert "реперов 3" in dec["why"]

def test_coverage_thirds_and_span():
    ok, why = sn.coverage_ok([1000.0, 1400.0, 1800.0, 2200.0, 2614.5])
    assert not ok and "нижней трети" in why
    ok, why = sn.coverage_ok([238.6, 500.0, 700.0, 900.0])
    assert not ok and "охват" in why
    ok, why = sn.coverage_ok([900.0, 1000.0, 1500.0, 1900.0])
    assert not ok and "охват" in why
    assert sn.coverage_ok([238.6, 911.2, 1764.5, 2614.5])[0]

def test_loo_gate_blocks(monkeypatch):
    lin = [-0.48, 1.62]
    rep = _flat_repers([238.6, 351.9, 583.2, 911.2, 1120.3, 1764.5, 2614.5], lin)
    fake = {"coeffs": list(lin), "deg": 1, "deg_table": [], "loo_rms": 5.0, "points": [], "chi2ndf": 0.0, "rms": 0.0}
    monkeypatch.setattr(sn, "select_degree_full", lambda *a, **k: dict(fake))
    dec = sn.decide_scale(np.full(2048, 100.0), [lin[0] + 40.0, lin[1]], rep, verbose=False)
    assert not dec["recommend"]
    assert "LOO" in dec["why"]
    assert "для замены нужно" in dec["need"] and "LOO 5.00" in dec["need"]

def test_independent_veto_and_main_case_not_vetoed():
    d = load_spectrum(MED)
    S = d["counts"].astype(float)
    cal = list(d["cal"])
    rep = _med_repers()
    dec = sn.decide_scale(S, cal, rep, verbose=False)
    assert len(dec["independent"]) >= 3
    if dec["rms_new"] > dec["rms_file"]:
        assert not dec["recommend"] and "хуже калибровки файла" in dec["why"]
    bad = [cal[0] + 15.0] + list(cal[1:])
    dec2 = sn.decide_scale(S, bad, rep, verbose=False)
    assert dec2["recommend"]
    assert len(dec2["independent"]) >= 3 and dec2["rms_new"] < dec2["rms_file"]
    assert not any("хуже калибровки файла" in w for w in dec2["warnings"])

def test_systematic_shift_warning():
    d = load_spectrum(MED)
    S = d["counts"].astype(float)
    cal = list(d["cal"])
    rep = [dict(r, E=r["E"] + 25.0) for r in _med_repers()]
    dec = sn.decide_scale(S, cal, rep, verbose=False)
    assert any("систематическая ошибка" in w for w in dec["warnings"])
    dec0 = sn.decide_scale(S, cal, _med_repers(), verbose=False)
    assert not any("систематическая ошибка" in w for w in dec0["warnings"])

def test_bg_scale_extrapolation_both_branches():
    B = np.full(2048, 100.0)
    C, how = sc.bg_own_scale({"coeffs": [0.0, 1.0], "points": [{"ch": 100.0, "E": 351.9}, {"ch": 900.0, "E": 1471.0}]}, [0.0, 1.0], B)
    assert "экстраполяция выше последней опоры 1471" in how
    C2, how2 = sc.bg_own_scale({"coeffs": [0.0, 1.0], "points": [{"ch": 100.0, "E": 351.9}, {"ch": 1900.0, "E": 2614.5}]}, [0.0, 1.0], B)
    assert "экстраполяция" not in how2

def test_bg_scale_quadratic_and_significance():
    d = load_spectrum(MED_BG)
    B = d["counts"].astype(float)
    cal = list(d["cal"])
    n = len(B)
    own = {"rejected": "тест", "points": [{"ch": float(poly_inverse(cal, E, hi=n + 200)), "E": E} for E in (351.932, 1120.294)]}
    C, how = sc.bg_own_scale(own, cal, B)
    assert "степени 2" in how
    d2 = load_spectrum(BG)
    B2 = d2["counts"].astype(float)
    cal2 = list(d2["cal"])
    n2 = len(B2)
    own2 = {"rejected": "тест", "points": [{"ch": float(poly_inverse(cal2, E, hi=n2 + 200)), "E": E} for E in (351.932, 1120.294)]}
    C2, how2 = sc.bg_own_scale(own2, cal2, B2)
    assert "по 2 опорам" in how2

def test_cli_error_codes(tmp_path):
    rc = cli.main(["sample-cal", "--bg", str(MED_BG), "--sample", str(tmp_path / "нет.xml"), "--out-zeroed", str(tmp_path / "z.xml"), "--out-signed", str(tmp_path / "s.xml")])
    assert rc == 3


def test_xml_degree_rule(monkeypatch):
    """Правка Claude (VERIFY8): решение о шкале XML -- степень не выше n-3, степень 3+ только от 8 реперов."""
    from gamma.calibration.intrinsic_shape import pipeline as pl
    monkeypatch.setattr(pl, "_guarded_calibrate", lambda cfg, ctx, h, v: {"loo_rms": 5.0 - cfg["degree"], "coeffs": [0.0] * (cfg["degree"] + 1), "points": []})
    for n, top in ((4, 1), (5, 2), (7, 2), (8, 4)):
        assert sn.select_degree_full({"points": [{}] * n}, {}, verbose=False)["deg"] == top


def test_independent_labels_do_not_depend_on_file_calibration():
    """Правка Claude (VERIFY8, N4): метки независимых линий -- по шкале реперов; при любой порче файла набор тот же."""
    from test_samplecal_net3 import _corrupt, DIST
    d = load_spectrum(MED); S = d["counts"].astype(float); cal = list(d["cal"]); rep = _med_repers()
    ref = [l["E"] for l in sn.decide_scale(S, cal, rep, verbose=False)["independent"]]
    assert len(ref) >= 3
    for f in DIST.values():
        assert [l["E"] for l in sn.decide_scale(S, _corrupt(cal, f), rep, verbose=False)["independent"]] == ref


def test_bg_extra_anchor_needs_significance():
    """Правка Claude (VERIFY8, N13): доп. опора фона со значимостью 3.5-4.6 (< 5) не берётся, с 6.6 -- берётся."""
    n, cal = 2048, [0.0, 1.4]; x = np.arange(n, dtype=float)
    own = {"rejected": "тест", "points": [{"ch": 351.932 / 1.4, "E": 351.932}, {"ch": 1120.294 / 1.4, "E": 1120.294}]}
    for area, k in ((450.0, 2), (800.0, 3)):
        B = np.random.default_rng(3).poisson(np.full(n, 200.0)).astype(float)
        c = 1764.491 / 1.4; s = sc.fwhm_guess(c, n) / 2.3548
        B += area / (np.sqrt(2 * np.pi) * s) * np.exp(-0.5 * ((x - c) / s) ** 2)
        assert "по %d опорам" % k in sc.bg_own_scale(own, cal, B)[1]


def test_second_pass_adds_line():
    """Правка Claude (VERIFY8, N17): реперы 338/911/968/1120/2615 -- первый проход по шкале реперов не узнаёт Pb-214 295.2,
    второй (шкала, уточнённая линиями первого прохода) -- узнаёт."""
    d = load_spectrum(MED); S = d["counts"].astype(float)
    rep = [r for r in _med_repers() if round(r["E"]) in (338, 911, 968, 1120, 2615)]
    dec = sn.decide_scale(S, list(d["cal"]), rep, verbose=False)
    first = [round(l["E"], 1) for l in sn._find_lines(S, dec["coeffs_new"], rep)]
    assert 295.2 not in first and 295.2 in [round(l["E"], 1) for l in dec["independent"]]
