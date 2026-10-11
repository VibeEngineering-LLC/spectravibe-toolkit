"""Тесты доработки по VERIFY7 «эталон -- энергии реперов»: рекомендация шкалы по энергиям реперов (главный случай -- испорченная калибровка
файла медальона), запас порога, применение только по рекомендации, CLI --use-reper-scale, вычитание на сетке пробы (LiveTime, дисперсия,
канал переполнения, округление и отсечение, h в XML), собственная шкала фона, сообщение переноса, шкала диалога. Сгенерировано gen_code.py
(qwen3.8:27b) по спецификации, 2026-10-11."""

import re
import math
import functools
import numpy as np
import pytest
from numpy.polynomial import polynomial as P
from gamma.calibration.intrinsic_shape.core import load_spectrum, poly_eval, poly_inverse, rebin_conserve
from gamma.calibration.intrinsic_shape import samplecal as sc
from gamma.calibration.intrinsic_shape import samplecal_q as q
from gamma.calibration.intrinsic_shape import samplecal_net as sn
from gamma.calibration.intrinsic_shape import pipeline as pl
from gamma.calibration.intrinsic_shape import cli
from test_samplecal_net import _med_repers, _write_like, _gauss, _flat_repers, BG, MED, MED_BG, ECAL_REF

EXTRA = (("Pb-214", 295.224), ("Bi-214", 609.312), ("Bi-214", 1238.11), ("Bi-214", 1377.669))
DIST = {"shift15": [15.0, 1.0], "gain103": [0.0, 1.03], "quad20": [0.0, 1.0, 20.0 / 2000.0 ** 2]}


@functools.lru_cache(maxsize=1)
def _truth():
    """Истина шкалы медальона: кубика по 9 реперам и 4 дополнительным линиям,
    подогнанным по каналам (калибровка файла не участвует)."""
    d = load_spectrum(MED)
    cal = list(d["cal"])
    ext = _med_repers([(float(poly_inverse(cal, E, hi=8400)), nuc, E) for nuc, E in EXTRA])
    allr = _med_repers() + ext
    return tuple(float(c) for c in P.polyfit([r["ch"] for r in allr], [r["E"] for r in allr], 3))


def _corrupt(cal, f):
    """Порча калибровки файла: E' = f(E_файла) как полином по каналу."""
    return [float(v) for v in functools.reduce(P.polyadd, [a * P.polypow(np.asarray(cal, float), k) for k, a in enumerate(f)])]


def _write_med(tmp_path, name, cal):
    text = MED.read_text("utf-8")
    coeff = "".join("<Coefficient>%r</Coefficient>" % float(c) for c in cal)
    text = re.sub(r"(<EnergyCalibration>.*?<PolynomialOrder>)\d+(</PolynomialOrder>\s*<Coefficients>).*?(</Coefficients>)",
                  lambda m: m.group(1) + str(len(cal) - 1) + m.group(2) + coeff + m.group(3), text, count=1, flags=re.S)
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


def _datapoints(path):
    return np.array([float(v) for v in re.findall(r"<DataPoint>([^<]*)</DataPoint>", path.read_text("utf-8"))])


def _set_live(path, live):
    text = path.read_text("utf-8")
    text = re.sub(r"<LiveTime>.*?</LiveTime>", "<LiveTime>%r</LiveTime>" % float(live), text, count=1)
    path.write_text(text, encoding="utf-8")


def _vis_pair(tmp_path, lt=3.0, overflow=False):
    """Пара для вычитания: фон -- новая пуассоновская реализация фикстуры фона,
    проба -- математическое ожидание фона, растянутое x' = 2 + 1.03x, умноженное на lt (LiveTime пробы = lt*LiveTime фона), плюс пик 662 кэВ;
    шум фона и пробы независимы."""
    d = load_spectrum(BG)
    B = d["counts"].astype(float)
    n = len(B)
    live_b = float(d["live"])
    rng = np.random.default_rng(7)
    bgp = _write_like(tmp_path, "b.xml", rng.poisson(B).astype(float))
    S0 = rebin_conserve(B, 2.0 + 1.03 * (np.arange(n + 1) - 0.5), n) * lt
    c = 2.0 + 1.03 * float(poly_inverse(ECAL_REF, 661.657))
    w = 1.03 * sc.fwhm_guess(c / 1.03, n)
    pk = _gauss(n, c, w, 30000.0)
    S = rng.poisson(S0 + pk).astype(float)
    if overflow:
        S[-1] = 1.0e6
    smp = _write_like(tmp_path, "s.xml", S)
    _set_live(smp, live_b * lt)
    return bgp, smp, c, w, float(pk.sum()), live_b


def _run(tmp_path, bgp, smp, lt_out="sample", mode=None):
    r = sn.run_sample_cal(bgp, smp, [], mode=mode, verbose=False)
    x = sn.run_sample_subtract(bgp, smp, r, r["fwhm_points"], tmp_path / "z.xml", tmp_path / "s2.xml", lt_out=lt_out, verbose=False)
    return r, x


@pytest.mark.parametrize("dist", sorted(DIST))
def test_main_case_corrupted_file_calibration(dist):
    d = load_spectrum(MED)
    S = d["counts"].astype(float)
    bad = _corrupt(list(d["cal"]), DIST[dist])
    dec = sn.decide_scale(S, bad, _med_repers(), verbose=False)
    assert dec["recommend"]
    assert not dec["replaced"] and dec["coeffs"] == bad
    assert len(dec["independent"]) >= 3
    assert dec["rms_new"] < dec["rms_file"]
    app = sn.apply_reper_scale(dec)
    assert app["replaced"]
    truth = list(_truth())
    for E in (300.0, 600.0, 1000.0, 1500.0, 2000.0):
        x = float(poly_inverse(truth, E, hi=8400))
        assert abs(float(poly_eval(app["coeffs"], x)) - E) <= 2.0


def test_control_not_recommended():
    d = load_spectrum(MED)
    S = d["counts"].astype(float)
    cal = list(d["cal"])
    rep = _med_repers()
    dec = sn.decide_scale(S, cal, rep, verbose=False)
    assert not dec["recommend"] and not dec["replaced"]
    assert dec["coeffs"] == cal
    assert "не превышает" in dec["why"] or "хуже калибровки файла" in dec["why"]   # правка Claude (VERIFY8): вето независимых линий
    own = math.sqrt(np.mean([(np.polyval(cal[::-1], r["ch"]) - r["E"]) ** 2 for r in rep if r["E"] >= 100.0]))
    assert abs(own - dec["rms_file_repers"]) < 1e-6
    with pytest.raises(ValueError, match="не рекомендована"):
        sn.apply_reper_scale(dec)


def test_margin_threshold():
    lin = [-0.48, 1.62]
    rep = _flat_repers([238.6, 351.9, 583.2, 911.2, 1120.3, 1764.5, 2614.5], lin)
    c = np.full(2048, 100.0)
    c[:16] = 0.0
    for shift, expected in ((0.9, False), (1.2, True)):
        dec = sn.decide_scale(c, [lin[0] + shift, lin[1]], rep, verbose=False)
        assert dec["loo_rms"] < 0.2
        assert dec["recommend"] is expected
        if expected is False:
            assert "не превышает" in dec["why"]


def test_repers_below_100_keV_ignored():
    lin = [-0.48, 1.62]
    rep = _flat_repers([50.0, 400.0, 1900.0, 2600.0], lin)
    c = np.full(2048, 100.0)
    dec = sn.decide_scale(c, [lin[0] + 5.0, lin[1]], rep, verbose=False)
    assert not dec["recommend"]
    assert "реперов 3" in dec["why"]


def test_cli_use_reper_scale(tmp_path):
    d = load_spectrum(MED)
    bad = _corrupt(list(d["cal"]), DIST["shift15"])
    smp = _write_med(tmp_path, "s.xml", bad)
    rep = _med_repers()
    sc.save_json(tmp_path / "r.json", rep, "t", len(d["counts"]), sc.anchors_of(d["counts"].astype(float), bad), cal=bad)
    base = ["sample-cal", "--bg", str(MED_BG), "--sample", str(smp), "--repers", str(tmp_path / "r.json")]
    assert cli.main(base + ["--out-zeroed", str(tmp_path / "z1.xml"), "--out-signed", str(tmp_path / "s1.xml")]) == 0
    assert np.allclose(load_spectrum(tmp_path / "z1.xml")["cal"], bad, rtol=1e-9, atol=1e-12)
    assert cli.main(base + ["--use-reper-scale", "--out-zeroed", str(tmp_path / "z2.xml"), "--out-signed", str(tmp_path / "s2.xml")]) == 0
    r = sn.run_sample_cal(MED_BG, smp, rep, verbose=False)
    assert r["decision"]["recommend"]
    got = list(load_spectrum(tmp_path / "z2.xml")["cal"])
    assert len(got) == len(r["decision"]["coeffs_new"])
    assert np.allclose(got, r["decision"]["coeffs_new"], rtol=1e-9, atol=1e-12)
    ok = _write_med(tmp_path, "ok.xml", list(d["cal"]))
    sc.save_json(tmp_path / "r0.json", rep, "t", len(d["counts"]), sc.anchors_of(d["counts"].astype(float), list(d["cal"])), cal=list(d["cal"]))
    assert cli.main(["sample-cal", "--bg", str(MED_BG), "--sample", str(ok), "--repers", str(tmp_path / "r0.json"), "--use-reper-scale", "--out-zeroed", str(tmp_path / "z3.xml"), "--out-signed", str(tmp_path / "s3.xml")]) == 2


def test_subtract_livetime_ratio_and_pulls(tmp_path):
    """Правка Claude: перевод через энергии с точными шкалами (подогнанное отображение по видимым структурам даёт
    систематику сдвига ~0.2 кан, она расширяет пулы до 1.35-1.55); шум фона и пробы независимы, LiveTime пробы = 3*фона."""
    d = load_spectrum(BG)
    B = d["counts"].astype(float)
    n = len(B)
    TS = [1.0, 1.28, 4.2e-4, 5.0e-8]
    rng = np.random.default_rng(7)
    bgp = _write_like(tmp_path, "b.xml", rng.poisson(B).astype(float), cal=ECAL_REF)
    xb = np.arange(n + 1) - 0.5
    S0 = rebin_conserve(B, np.array([float(poly_inverse(TS, poly_eval(ECAL_REF, v), hi=n + 400)) for v in xb]), n) * 3.0
    smp = _write_like(tmp_path, "s.xml", rng.poisson(S0).astype(float), cal=TS)
    _set_live(smp, float(d["live"]) * 3.0)
    r, x = _run(tmp_path, bgp, smp, mode="energy")
    assert abs(x["k"] - 3.0) < 1e-9
    res, var = x["resid"], x["var"]
    m = np.zeros(n, bool)
    m[40:n - 300] = True
    pulls = res[m] / np.sqrt(var[m])
    assert 0.85 <= float(np.std(pulls)) <= 1.15
    assert abs(float(np.mean(pulls))) < 0.1


def test_overflow_channel_zeroed(tmp_path):
    bgp, smp, c, w, area, live_b = _vis_pair(tmp_path, lt=1.0, overflow=True)
    r, x = _run(tmp_path, bgp, smp)
    assert abs(x["resid"][-1]) < 1000.0
    assert _datapoints(tmp_path / "s2.xml")[-1] > -1000.0 and _datapoints(tmp_path / "z.xml")[-1] < 1000.0


def test_xml_zeroed_signed_rounding(tmp_path):
    bgp, smp, c, w, area, live_b = _vis_pair(tmp_path, lt=3.0)
    r, x = _run(tmp_path, bgp, smp)
    dz, ds = _datapoints(tmp_path / "z.xml"), _datapoints(tmp_path / "s2.xml")
    assert (ds < 0).any()
    assert (dz >= 0).all()
    assert np.array_equal(ds, np.round(x["resid"]))
    assert np.array_equal(dz, np.maximum(0.0, np.round(x["resid"])))


def test_lt_out_background(tmp_path):
    bgp, smp, c, w, area, live_b = _vis_pair(tmp_path, lt=3.0)
    r, x = _run(tmp_path, bgp, smp, lt_out="bg")
    f = x["live_b"] / x["live_s"]
    assert abs(f - 1.0 / 3.0) < 1e-9
    ds = _datapoints(tmp_path / "s2.xml")
    assert np.array_equal(ds, np.round(x["resid"] * f))
    lt = float(re.search(r"<LiveTime>([^<]*)</LiveTime>", (tmp_path / "s2.xml").read_text("utf-8")).group(1))
    assert abs(lt - x["live_b"]) < 1e-6


def test_visible_h_written_to_xml(tmp_path):
    bgp, smp, c, w, area, live_b = _vis_pair(tmp_path, lt=1.0)
    r, x = _run(tmp_path, bgp, smp)
    h = r["transfer"]["h"]
    assert abs(h[1] - 1.0) > 0.01
    text = (tmp_path / "z.xml").read_text("utf-8")
    assert ("Отображение h(x)=%.7g+%.7g*x" % (h[0], h[1])) in text


def test_energy_mode_uses_own_background_scale(tmp_path, monkeypatch):
    d = load_spectrum(BG)
    B = d["counts"].astype(float)
    n = len(B)
    cal_file = [c * 1.02 for c in ECAL_REF]
    bgp = _write_like(tmp_path, "b.xml", B, cal=cal_file)
    smp = _write_like(tmp_path, "s.xml", B, cal=ECAL_REF)
    pts = [{"ch": float(poly_inverse(ECAL_REF, E, hi=n + 200)), "E": E} for E in (351.932, 1120.294)]
    monkeypatch.setattr(pl, "run_ecal", lambda *a, **k: {"rejected": "тест", "points": pts, "coeffs": []})
    r = sn.run_sample_cal(bgp, smp, [], mode="energy", verbose=False)
    cb = r["transfer"]["cal_b"]
    assert list(cb) != cal_file
    for p in pts:
        assert abs(float(poly_eval(cb, p["ch"])) - p["E"]) < 0.5
    assert "экстраполяция выше последней опоры" in r["transfer"]["bg_scale_how"]
    assert any("экстраполирована" in w for w in r["warnings"])   # правка Claude (VERIFY8): отдельное предупреждение в окне


def test_background_scale_extra_anchors():
    d = load_spectrum(MED_BG)
    B = d["counts"].astype(float)
    cal = list(d["cal"])
    n = len(B)
    own = {"rejected": "тест", "points": [{"ch": float(poly_inverse(cal, E, hi=n + 200)), "E": E} for E in (351.932, 1120.294)]}
    C, how = sc.bg_own_scale(own, cal, B)
    m = re.search(r"по (\d+) опорам", how)
    assert m and int(m.group(1)) >= 3
    assert "экстраполяция" not in how


def test_transfer_message_receiver_calibration(tmp_path):
    d = load_spectrum(BG)
    B = d["counts"].astype(float)
    cal = list(d["cal"])
    intr, anc = sc.intrinsic_repers(B, cal)
    sc.save_json(tmp_path / "r.json", intr, "t", len(B), anc, cal=cal)
    data = sc.load_json(tmp_path / "r.json")
    with pytest.raises(ValueError, match="неверна калибровка файла приёмника"):
        sc.apply_to(data, B, [c * 1.2 for c in cal], verbose=False)


def test_dialog_scale_from_repers():
    assert q.current_scale([{"ch": 100.0, "E": 150.0}]) is None
    two = [{"ch": 100.0, "E": 150.0}, {"ch": 500.0, "E": 750.0}]
    C2 = q.current_scale(two)
    assert len(C2) == 2
    assert all(abs(float(poly_eval(C2, r["ch"])) - r["E"]) < 1e-9 for r in two)
    four = [{"ch": x, "E": 1.0 + 1.5 * x + 1e-4 * x * x} for x in (100.0, 400.0, 900.0, 1500.0)]
    C4 = q.current_scale(four)
    assert len(C4) == 3 and abs(C4[2] - 1e-4) < 1e-9


def test_margin_ratio_and_difference(monkeypatch):
    """Правка Claude: оба условия запаса по отдельности -- LOO шкалы по реперам задаётся подменой select_degree_full,
    СКО файла на точных реперах равно сдвигу калибровки файла."""
    lin = [-0.48, 1.62]
    rep = _flat_repers([238.6, 351.9, 583.2, 911.2, 1120.3, 1764.5, 2614.5], lin)
    c = np.full(2048, 100.0)
    for loo, shift, expected in ((2.5, 3.6, False), (2.5, 3.9, True), (0.5, 1.2, False), (0.5, 1.6, True)):
        fake = {"coeffs": list(lin), "deg": 1, "deg_table": [], "loo_rms": loo, "points": [], "chi2ndf": 0.0, "rms": 0.0}
        monkeypatch.setattr(sn, "select_degree_full", lambda *a, **k: dict(fake))
        dec = sn.decide_scale(c, [lin[0] + shift, lin[1]], rep, verbose=False)
        assert abs(dec["rms_file_repers"] - shift) < 1e-6
        assert dec["recommend"] is expected, (loo, shift)


def test_degree4_needs_seven_repers(monkeypatch):
    """Правка Claude (VERIFY7 п.4): LOO убывает со степенью -- степень 4 выбирается только при >= 7 реперах."""
    calls = []

    def fake(cfg, ctx, h, verbose):
        calls.append(cfg["degree"])
        return {"loo_rms": 5.0 - cfg["degree"], "coeffs": [0.0] * (cfg["degree"] + 1), "points": []}

    monkeypatch.setattr(pl, "_guarded_calibrate", fake)
    for n, top in ((6, 3), (7, 4)):
        calls.clear()
        res = q.select_degree({"points": [{}] * n}, {}, verbose=False)
        assert max(calls) == top and res["deg"] == top
