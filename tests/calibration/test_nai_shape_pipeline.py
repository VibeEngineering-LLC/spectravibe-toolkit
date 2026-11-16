import numpy as np
import pytest
from scipy.special import ndtr
from numpy.random import default_rng

from gamma.calibration.nai_shape import pipeline, search, checks, model, validate, peakreg, surrogate
from gamma.calibration.intrinsic_shape.core import poly_eval


# ─────────────────────────────────────────────
# Вспомогательные функции
# ─────────────────────────────────────────────

def synth(n=1024, kpc=3.0, scale=1.0, extra=()):
    """Синтетический спектр: континуум + гауссовы пики (интеграл через ndtr)."""
    ch = np.arange(n)
    E = kpc * ch
    # Континуум: exp-затухание + постоянный фон
    result = scale * (4000.0 * np.exp(-E / 400.0) + 30.0)

    # Основные пики: (E0_keV, площадь)
    base_peaks = [(238.6, 6000.0), (609.3, 3000.0), (1460.8, 8000.0),
                  (1764.5, 1500.0), (2614.5, 2500.0)]
    for E0, A in base_peaks:
        fwhm = 0.07 * np.sqrt(661.7 * E0)
        sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        lo = E - kpc / 2.0
        hi = E + kpc / 2.0
        result += scale * A * (ndtr((hi - E0) / sigma) - ndtr((lo - E0) / sigma))

    # Дополнительные пики (без умножения на scale)
    for E0, A in extra:
        fwhm = 0.07 * np.sqrt(661.7 * E0)
        sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        lo = E - kpc / 2.0
        hi = E + kpc / 2.0
        result += A * (ndtr((hi - E0) / sigma) - ndtr((lo - E0) / sigma))

    return result


def spec(counts, live, name="s"):
    """Формирует словарь-спектр, ожидаемый интерфейсами."""
    return {
        "counts": np.array(counts, dtype=float),
        "live": live,
        "real": live,
        "cal": (0.0, 3.0, 0.0),
        "name": name,
        "path": name + ".spe",
        "n_read": len(counts),
        "start": "",
    }


# Профиль калибровки — фиксированный словарь параметров
PROFILE = {
    "fit_range_kev": [100.0, 2800.0],
    "e_strong": 1460.8,
    "fwhm_rel662": 0.07,
    "search_rel": 0.05,
    "min_z": 5.0,
    "calib_lines": [
        {"E": 238.6, "sys_kev": 0.0},
        {"E": 609.3, "sys_kev": 0.0},
        {"E": 1460.8, "sys_kev": 0.0},
        {"E": 1764.5, "sys_kev": 0.0},
        {"E": 2614.5, "sys_kev": 0.0},
    ],
    "gain_span": 0.04,
    "offset_span_ch": 8.0,
    "grid_step_fwhm": 0.25,
    "n_starts": 5,
    "bin_fwhm": 0.5,
    "hats_background": True,
    "hat_spacing_bg_fwhm": 4.0,
    "free_background": ["d", "gN"],
    "free_sample": ["d", "gN"],
    "hat_spacing_fwhm": 3.0,
    "test_kev": [238.6, 609.3, 1460.8, 2614.5],
    "null_windows": [
        {"E": 238.6, "tag": ""},
        {"E": 609.3, "tag": ""},
        {"E": 1460.8, "tag": ""},
        {"E": 1764.5, "tag": ""},
        {"E": 2614.5, "tag": ""},
    ],
    "radon_variable_kev": [],
    "z_thr": 3.0,
    "sample_lines": {"K-40": [1460.8], "Cs-137": [661.7]},
    "shape_min_bg_fraction": 0.5,
    "peaks_deg": 1,
    "anchor_strong": {"Cs-137": 661.7, "K-40": 1460.8},
    "anchor_lines": {
        "Cs-137": [{"E": 661.7, "sys_kev": 0.0}],
        "K-40": [{"E": 1460.8, "sys_kev": 0.0}],
    },
}


def bg_pair(c_true, seed):
    """Генерирует пару (фон, образец) с известным дрейфом c_true."""
    rng = default_rng(seed)
    truth = synth(scale=20.0)
    B = rng.poisson(truth)
    base = rng.poisson(truth)
    S = surrogate.drift_events(base, c_true, 1024, rng, thin=0.5)
    return (spec(B, 1000.0, "bg"), spec(S, 500.0, "sm"))


# ─────────────────────────────────────────────
# Тесты
# ─────────────────────────────────────────────

def test_run_pair_k_is_live_ratio():
    """Коэффициент k равен отношению live-времени образца к фону."""
    bg, sm = bg_pair([0.0, 1.0, 0.0], 11)
    r = pipeline._run_pair(bg, sm, PROFILE, mode="background")
    assert abs(r["k"] - 0.5) < 1e-12


def test_run_pair_net_zero_for_background():
    """При нулевом дрейфе нетто-сигнал статистически нullevel."""
    bg, sm = bg_pair([0.0, 1.0, 0.0], 12)
    r = pipeline._run_pair(bg, sm, PROFILE, mode="background")
    a = r["arrays"]
    m = a["mask"]
    z = a["net"][m].sum() / np.sqrt(a["var"][m].sum())
    assert abs(z) < 4


def test_run_pair_recovers_drift():
    """Пайплайн восстанавливает дрейф в пределах 4σ при σ_h < 0.3 канала."""
    c_true = [2.0, 1.012, 0.0]
    bg, sm = bg_pair(c_true, 13)
    r = pipeline._run_pair(bg, sm, PROFILE, mode="background")
    xs = np.array(r["x_test"])
    err = np.abs(poly_eval(r["c"], xs) - poly_eval(c_true, xs))
    sig = np.array(r["sigma_h_test"])
    assert np.all(err < 4 * sig)
    assert np.all(sig < 0.3)


def test_anchors_to_p_roundtrip():
    """Вектор параметров p точно воспроизводит заданные сдвиги на якорях."""
    p = search.anchors_to_p(1.5, -2.0, 200.0, 800.0, 1024)
    f = lambda y: p[0] + p[1] * y / 1024
    assert abs(f(200.0) - 1.5) < 1e-9
    assert abs(f(800.0) + 2.0) < 1e-9


def test_verdict_flags_dipole():
    """Вердикт помечает линию как bad при z_dip > порога, иначе ok."""
    t = [{"E": 1.0, "excluded": False, "z_sum": 0.0, "z_dip": 5.0}]
    v = checks.verdict(t, 3.0, None)
    assert v["ok"] is False
    assert v["bad"] == [1.0]
    # При z_dip ниже порога — всё хорошо
    t2 = [{"E": 1.0, "excluded": False, "z_sum": 0.0, "z_dip": 2.0}]
    v2 = checks.verdict(t2, 3.0, None)
    assert v2["ok"] is True


def test_h_sigma_analytic():
    """Аналитическая формула σ_h для p=0, free=[0,1], h0=[0,1,0]."""
    cov = np.diag([0.04, 0.09])
    chans = np.array([0.0, 500.0, 1000.0])
    s = model.h_sigma(np.zeros(4), cov, [0, 1], [0.0, 1.0, 0.0], 1000, chans)
    expected = np.sqrt(0.04 + 0.09 * (chans / 1000.0) ** 2)
    assert np.allclose(s, expected, rtol=1e-6)


def test_surrogate_cases_independent_streams(monkeypatch):
    """Каждое (случай, seed) генерирует независимую реализацию фона."""
    saved = []

    def fake_run_pair(bg, sm, profile, mode="sample", nuclides=None, free=None, verbose=False):
        saved.append(bg["counts"].copy())
        return {
            "c": [0.0, 1.0, 0.0],
            "x_test": [100.0, 200.0, 300.0, 400.0],
            "sigma_h_test": [1.0] * 4,
            "chi2ndf": 1.0,
        }

    monkeypatch.setattr(validate.pipeline, "_run_pair", fake_run_pair)
    bg = spec(np.rint(synth(scale=1.0)), 1000.0)
    prof = dict(PROFILE)
    validate.surrogate_cases(bg, prof, [{"d": 0.0, "g": 1.0}, {"d": 1.0, "g": 1.0}], seeds=(1, 2))
    # Порядок: (случай0, seed1), (случай0, seed2), (случай1, seed1), (случай1, seed2)
    assert len(saved) == 4
    # Разные случаи при одном seed — разные реализации
    assert not np.array_equal(saved[0], saved[2])
    # Один случай при разных seeds — тоже разные реализации
    assert not np.array_equal(saved[0], saved[1])


def test_segment_rejects_non_nested(monkeypatch):
    """segment() бросает ValueError при невыполненном условии вложенности; иначе вычитает."""
    def fake_load(path):
        if path == "b":
            return spec(np.array([5.0, 5.0, 5.0]), 10.0)
        elif path == "a":
            return spec(np.array([1.0, 6.0, 1.0]), 5.0)
        elif path == "a2":
            return spec(np.array([1.0, 1.0, 1.0]), 5.0)
        raise FileNotFoundError(path)

    monkeypatch.setattr(validate, "load_spectrum", fake_load)
    # a[1]=6 > b[1]=5 → не вложен → ошибка
    with pytest.raises(ValueError):
        validate.segment("b", "a")
    # a2 везде <= b → корректное вычитание
    r = validate.segment("b", "a2")
    assert np.allclose(r["counts"], [4.0, 4.0, 4.0])
    assert r["live"] == 5.0


def test_register_by_lines_recovers_drift():
    """peakreg._register_by_lines восстанавливает дрейф по линиям Cs-137."""
    rng = default_rng(21)
    bg = spec(rng.poisson(synth(scale=20.0)), 1000.0, "bg")
    truth_s = synth(scale=10.0, extra=((661.7, 3.0e6),))
    c_true = [3.0, 1.015, 0.0]
    S = surrogate.drift_events(rng.poisson(truth_s), c_true, 1024, rng)
    sm = spec(S, 1000.0, "sm")
    r = peakreg._register_by_lines(bg, sm, PROFILE, ["Cs-137"])
    xs = np.array(r["x_test"])
    err = np.abs(poly_eval(r["c"], xs) - poly_eval(c_true, xs))
    assert np.all(err < 4 * np.array(r["sigma_h_test"]) + 0.05)


def test_run_auto_selects_peaks_for_strong_sample():
    """run_auto выбирает режим 'peaks' при сильном сигнале образца (bg_fraction < 0.5)."""
    rng = default_rng(22)
    # ручная правка: фон снят за 20000 с (k = 0.05) -- иначе доля фона k·ΣB/ΣS > 1 и образец не «сильный»
    bg = spec(rng.poisson(synth(scale=20.0)), 20000.0, "bg")
    truth_s = synth(scale=10.0, extra=((661.7, 3.0e6),))
    c_true = [3.0, 1.015, 0.0]
    S = surrogate.drift_events(rng.poisson(truth_s), c_true, 1024, rng)
    sm = spec(S, 1000.0, "sm")
    r = peakreg.run_auto(bg, sm, PROFILE, ["Cs-137"])
    assert r["mode_auto"] == "peaks"
    assert r["bg_fraction"] < 0.5
