import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np, pytest
from scipy.special import ndtr
from numpy.random import default_rng

from gamma.calibration.nai_shape import checks, peakreg, strongcheck, surrogate
from gamma.calibration.intrinsic_shape.core import poly_eval
from test_nai_shape_pipeline import synth, spec, PROFILE


def one_line_pair(seed, c_true):
    """Генерирует пару (фон, образец) с одной линией Cs-137 на слабом континууме."""
    rng = default_rng(seed)
    # Фон: только континуум и линии фона, масштаб 20
    bg = spec(rng.poisson(synth(scale=20.0)), 1000.0, "bg")
    # Образец: плоский континуум 5.0 + одна линия Cs-137
    truth = 5.0 + synth(scale=0.0, extra=((661.7, 2.0e6),))
    sm = spec(surrogate.drift_events(rng.poisson(truth), c_true, 1024, rng), 1000.0, "sm")
    return bg, sm


# --- Тесты verdict ---

def test_verdict_flags_negative_dipole():
    """Отрицательный z_dip должен сделать вердикт неудачным."""
    t = [{"E": 1.0, "excluded": False, "z_sum": 0.0, "z_dip": -5.0}]
    assert checks.verdict(t, 3.0, None)["ok"] is False


def test_verdict_flags_z_sum():
    """Отрицательный z_sum флажит; при use_sum_for=[] z_sum игнорируется."""
    t = [{"E": 1.0, "excluded": False, "z_sum": 5.0, "z_dip": 0.0}]
    # z_sum = 5.0 > 3.0 → ok = False
    assert checks.verdict(t, 3.0, None)["ok"] is False
    # use_sum_for пустой список → z_sum не проверяется → ok = True
    assert checks.verdict(t, 3.0, [])["ok"] is True


# --- Тесты peakreg (режим peaks1) ---

def test_peaks1_recovers_offset():
    """В режиме peaks1 сдвиг (c[0]) восстанавливается в пределах неопределённости."""
    c_true = [3.0, 1.0, 0.0]
    bg, sm = one_line_pair(31, c_true)
    r = peakreg._register_by_lines(bg, sm, PROFILE, ["Cs-137"])
    assert r["mode"] == "peaks1"
    xs = np.array(r["x_test"])
    err = abs(poly_eval(r["c"], xs) - poly_eval(c_true, xs))
    # Ошибка меньше 4σ + допуск
    assert np.all(err < 4 * np.array(r["sigma_h_test"]) + 0.1)


def test_peaks1_sigma_includes_gain_sys():
    """Систематическая неопределённость gain_sys входит в σ_h_test."""
    bg, sm = one_line_pair(32, [0.0, 1.0, 0.0])
    r = peakreg._register_by_lines(bg, sm, PROFILE, ["Cs-137"])
    xl = r["one_line"]["x"]
    x = np.array(r["x_test"])
    s = np.array(r["sigma_h_test"])
    # Точка, самая далёкая от линии
    i = int(np.argmax(np.abs(x - xl)))
    gain_sys = PROFILE.get("one_line_gain_sys", 0.005)
    assert s[i] > 0.5 * gain_sys * abs(x[i] - xl)


def test_peaks_sigma_matches_scatter():
    """σ_h по Монте-Карло согласована с разбросом o(900) между реализациями."""
    h = []
    sig = []
    for seed in range(100, 108):  # 8 независимых реализаций
        rng = default_rng(seed)
        bg = spec(rng.poisson(synth(scale=20.0)), 1000.0, "bg")
        truth_s = synth(scale=10.0, extra=((661.7, 3.0e6),))
        sm = spec(surrogate.drift_events(rng.poisson(truth_s), [3.0, 1.015, 0.0], 1024, rng), 1000.0, "sm")
        r = peakreg._register_by_lines(bg, sm, PROFILE, ["Cs-137"])
        # Значение полинома и σ в точке E=900 кэВ (линейная интерполяция)
        h.append(poly_eval(r["c"], 900.0))
        sig.append(float(np.interp(900.0, r["x_test"], r["sigma_h_test"])))
    ratio = np.std(h, ddof=1) / np.mean(sig)
    assert 0.4 < ratio < 2.5


# --- Тесты peakreg.run_auto ---

def test_run_auto_flags_unverified_scale():
    """При длительном фоне (20000 с) шкала считается непроверенной."""
    rng = default_rng(40)
    # Фон снят за 20000 с — доля фона в образце мала → шкала не подтверждена
    bg = spec(rng.poisson(synth(scale=20.0)), 20000.0, "bg")
    truth_s = synth(scale=10.0, extra=((661.7, 3.0e6),))
    sm = spec(surrogate.drift_events(rng.poisson(truth_s), [3.0, 1.015, 0.0], 1024, rng), 1000.0, "sm")
    r = peakreg.run_auto(bg, sm, PROFILE, ["Cs-137"])
    assert r["scale_verified"] is False
    assert r["flag"] == "шкала не проверена"
    assert r["verdict"]["counted"] is False


def test_run_auto_verified_for_background_like_sample():
    """Образец того же состава, что и фон → шкала подтверждена."""
    rng = default_rng(41)
    bg = spec(rng.poisson(synth(scale=20.0)), 1000.0, "bg")
    sm = spec(rng.poisson(synth(scale=10.0)), 500.0, "sm")  # доля фона ≈ 1
    r = peakreg.run_auto(bg, sm, PROFILE, ["K-40"])
    assert r["scale_verified"] is True
    assert r["flag"] == ""


# --- Тесты strongcheck ---

def test_zones_compton_backscatter():
    """При I >= icomp добавляются энергии Комптона и обратного комптона."""
    z = strongcheck.zones([(661.7, 85.0)])
    Ec = 2 * 661.7**2 / (511 + 2 * 661.7)
    Eb = 661.7 * 511 / (511 + 2 * 661.7)
    assert len(z) == 3
    assert np.allclose(sorted(z), sorted([661.7, Ec, Eb]))
    # При I < icomp — только сама линия
    assert strongcheck.zones([(661.7, 2.0)]) == [661.7]


def test_select_windows_drops_near_zone():
    """Окно, близкое к зоне, отбрасывается; далёкое — сохраняется."""
    kept, dropped = strongcheck.select_windows([600.0, 1500.0], [620.0], lambda E: 30.0)
    assert kept == [1500.0]
    assert dropped == [600.0]
