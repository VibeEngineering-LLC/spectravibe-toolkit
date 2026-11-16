import numpy as np
import pytest
from scipy.special import ndtr

from gamma.calibration.nai_shape import surrogate, checks, model, binned, bgcal
from gamma.calibration.intrinsic_shape.core import poly_eval


def synth_bg(n=1024, keV_per_ch=3.0, scale=1.0):
    """Ожидаемый спектр NaI-фона: континуум + гауссовы пики (мат. ожидание, без шума)."""
    ch = np.arange(n, dtype=float)
    E = keV_per_ch * ch  # энергия канала, кэВ

    # Континуальный фон
    spectrum = scale * 4000.0 * np.exp(-E / 400.0) + scale * 30.0

    # Гауссовы пики (только те, что попадают в шкалу)
    peaks = [(238.6, 6000), (609.3, 3000), (1460.8, 8000), (1764.5, 1500), (2614.5, 2500)]
    E_max = keV_per_ch * (n - 1)

    for E0, A in peaks:
        if E0 > E_max:
            continue
        fwhm_kev = 0.07 * np.sqrt(661.7 * E0)
        sigma_kev = fwhm_kev / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        mu_ch = E0 / keV_per_ch
        sigma_ch = sigma_kev / keV_per_ch
        # Интеграл гауссианы по границам канала через ndtr
        lo = (ch - mu_ch) / sigma_ch
        hi = (ch + 1.0 - mu_ch) / sigma_ch
        spectrum += scale * A * (ndtr(hi) - ndtr(lo))

    return spectrum


# ─────────────────────────────────────────────────────────────────────────────
# Тесты
# ─────────────────────────────────────────────────────────────────────────────

def test_split_poisson_conserves_counts():
    """Разбиение Пуассона сохраняет общее число событий поканально."""
    counts = np.rint(synth_bg())
    rng = np.random.default_rng(1)
    a, b = surrogate.split_poisson(counts, 0.3, rng)
    assert np.array_equal(a + b, counts)
    assert 0.25 < a.sum() / counts.sum() < 0.35


def test_inverse_h_roundtrip():
    """Обратная функция: poly_eval(c, inverse_h(c, y)) ≈ y."""
    c = [5.0, 1.03, 2e-5]
    y = np.linspace(0, 1000, 101)
    x = surrogate.inverse_h(c, y)
    assert np.max(np.abs(poly_eval(c, x) - y)) < 1e-8


def test_drift_events_moves_peak():
    """Пик после дрейфа смещён на величину h⁻¹(y) − x."""
    spec = np.zeros(2000)
    spec[1000] = 200000.0
    c = [10.0, 1.05, 0.0]
    out = surrogate.drift_events(spec, c, 2000, np.random.default_rng(2))
    centroid = np.sum(np.arange(2000) * out) / np.sum(out)
    expected = (1000 - 10) / 1.05
    assert abs(centroid - expected) < 0.05


def test_blur_var_conserves_sum():
    """Размытие с дисперсией: сумма сохранена, максимум не вырос, дисперсии ≥ 0."""
    y = synth_bg()
    sigma = np.full(len(y), 3.0)
    out, vout = model.blur_var(y, y, sigma)
    rel = abs(out[50:-50].sum() - y[50:-50].sum()) / y[50:-50].sum()
    assert rel < 1e-3
    assert (vout >= 0).all()
    assert out.max() < y.max()


def test_sum_bins_conserves():
    """Сумма по бинам равна сумме событий в исходном окне каналов."""
    y = np.arange(100, dtype=float)
    edges = binned.make_edges(10, 89, lambda x: 16.0, 0.25)
    assert edges[0] == 9.5
    assert edges[-1] == 89.5
    assert binned.sum_bins(y, edges).sum() == y[10:90].sum()


def test_binned_recovers_drift():
    """Фит бинированной модели восстанавливает истинные коэффициенты дрейфа."""
    n = 1024
    truth = synth_bg(n, scale=40.0)
    rng = np.random.default_rng(3)
    B = rng.poisson(truth).astype(float)
    base = rng.poisson(truth).astype(float)
    c_true = [3.0, 1.012, 0.0]
    S = surrogate.drift_events(base, c_true, n, rng)

    mask = np.zeros(n, dtype=bool)
    mask[34:931] = True  # активная область каналов

    def fw(x):
        """Полная ширина на полумесяце в каналах."""
        return np.maximum(1.0, 0.07 * np.sqrt(661.7 * np.maximum(3.0 * x, 1.0)) / 3.0)  # векторизуемо (binned v2)

    prob = binned.BinnedProblem(S, B, 1.0, mask, [0.0, 1.0, 0.0], fw_ch=fw)
    res = model.fit(prob, [0, 1], np.zeros(4), np.arange(-8, 8.1, 1.0), np.arange(-40, 40.1, 4.0), 3)
    cv = model.covariance(prob, res["p"], [0, 1])
    c_fit = model.coefs(res["p"], [0, 1, 0], n)

    # Сравниваем предсказание модели с истинным дрейфом в нескольких каналах
    # допуск 4σ_h (σ из ковариации, ~0.01-0.06 кан при этом счёте) вместо фиксированных 0.5 кан (VERIFY.md S4: сдвиг 0.4 выживал)
    xs = np.array([100.0, 300.0, 500.0, 800.0])
    sig = model.h_sigma(res["p"], cv["cov"], [0, 1], [0.0, 1.0, 0.0], n, xs)
    assert np.all(np.abs(poly_eval(c_fit, xs) - poly_eval(c_true, xs)) < 4 * sig) and np.all(sig < 0.2)


def test_binned_chi2_no_ripple():
    """Chi² монотонно не убывает при увеличении сдвига (отсутствие «ряби»)."""
    n = 1024
    truth = synth_bg(n, scale=40.0)
    B = np.rint(truth)
    S = np.rint(truth)

    mask = np.zeros(n, dtype=bool)
    mask[34:931] = True

    def fw(x):
        return np.maximum(1.0, 0.07 * np.sqrt(661.7 * np.maximum(3.0 * x, 1.0)) / 3.0)  # векторизуемо (binned v2)

    prob = binned.BinnedProblem(S, B, 1.0, mask, [0.0, 1.0, 0.0], fw_ch=fw)

    shifts = [0.0, 0.25, 0.5, 0.75, 1.0]
    chi2_vals = [prob.chi2([d, 0.0, 0.0, 0.0]) for d in shifts]
    for i in range(len(chi2_vals) - 1):
        # строго возрастает: χ² чувствителен к дробному сдвигу (нет ступенек/плато от привязки к каналам)
        assert chi2_vals[i + 1] > chi2_vals[i] + 1e-3


def test_window_tests_dipole():
    """Дипольная аномалия детектируется; окно с excluded помечается корректно."""
    net = np.zeros(1000)
    var = np.ones(1000) * 100.0
    e2ch = lambda E: E / 3.0
    fwhm = lambda E: 60.0

    # Без сигнала: обе z-статистики нулевые
    res = checks.window_tests(net, var, e2ch, fwhm, [1500.0])
    assert abs(res[0]["z_sum"]) < 1e-9
    assert abs(res[0]["z_dip"]) < 1e-9

    # Диполь: +50 ниже пика, −50 выше (вокруг канала 500 ≈ E=1500 кэВ)
    net[480:500] += 50.0
    net[501:521] -= 50.0
    res = checks.window_tests(net, var, e2ch, fwhm, [1500.0])
    assert res[0]["z_dip"] < -5.0

    # Исключённая линия → excluded = True
    res = checks.window_tests(net, var, e2ch, fwhm, [1500.0], exclude_kev=[1500.0])
    assert res[0]["excluded"] is True


def test_bgcal_recovers_scale():
    """Калибровка восстанавливает энергетическую шкалу с ошибкой < 3 кэВ."""
    n = 1024
    truth = synth_bg(n, scale=20.0)
    counts = np.random.default_rng(4).poisson(truth).astype(float)
    # Файл-шкала с намеренной ошибкой 2 %
    file_cal = [0.0, 3.0 * 1.02, 0.0]
    lines = [{"E": 238.6}, {"E": 609.3}, {"E": 1460.8}, {"E": 1764.5}, {"E": 2614.5}]
    r = bgcal.calibrate(counts, file_cal, lines, min_z=5.0)

    accepted = [ln for ln in r["lines"] if ln["accepted"]]
    assert len(accepted) >= 4
    for ln in accepted:
        E = ln["E"]
        assert abs(poly_eval(r["ecal"], E / 3.0) - E) < 3.0


def test_coefs_quadratic_scaling():
    """Квадратичный член p[2] масштабируется как p[2]/n_b²; линейный — p[1]/n_b."""
    # p[2]=1, n_b=1000 → c2 = 1e-6, h(1000)−1000 = 1e-6·10⁶ = 1
    c = model.coefs([0.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0], 1000)
    assert abs(poly_eval(c, 1000) - 1000 - 1.0) < 1e-9

    # p=[2, 10, 0, 0] → c = [2, 1.01, 0]
    c2 = model.coefs([2.0, 10.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1000)
    assert np.allclose(c2, [2.0, 1.01, 0.0])
