"""Суррогатные спектры с известным дрейфом."""

import numpy as np

__all__ = ["split_poisson", "drift_events", "inverse_h"]


def split_poisson(counts, frac, rng):
    """Биномиальное прореживание спектра на две независимые пуассоновские части.

    Параметры
    ----------
    counts : array_like
        Числа отсчётов (могут быть нецелыми).
    frac : float
        Доля, приходящаяся на первую часть (0 ≤ frac ≤ 1).
    rng : numpy.random.Generator
        Генератор случайных чисел.

    Возврат
    -------
    tuple of ndarray
        (a, b) — две float-массива, сумма которых равна rint(counts).
    """
    n = np.rint(counts).astype(np.int64)
    n = np.maximum(n, 0)
    a = rng.binomial(n, frac)
    b = n - a
    return a.astype(np.float64), b.astype(np.float64)


def inverse_h(c, y):
    """Обратное отображение х = h⁻¹(у) для c = [c0, c1, c2].

    h(x) = c0 + c1*x + c2*x**2.

    Параметры
    ----------
    c : sequence of float, длина 3
        Коэффициенты [c0, c1, c2].
    y : float или ndarray
        Значение(я) на шкале фона.

    Возврат
    -------
    ndarray или float
        Соответствующее(е) значение(я) на шкале образца.
    """
    c0, c1, c2 = c
    y = np.asarray(y, dtype=np.float64)
    if abs(c2) < 1e-15:
        return (y - c0) / c1
    # Устойчивая формула квадратного корня (избегает потери точности)
    disc = c1 * c1 + 4.0 * c2 * (y - c0)
    disc = np.maximum(disc, 0.0)
    return 2.0 * (y - c0) / (c1 + np.sqrt(disc))


def drift_events(counts, c, n_out, rng, thin=1.0, sigma_extra=0.0, ref_ch=None):
    """Покомпонентное моделирование событий с дрейфом шкалы.

    Параметры
    ----------
    counts : array_like
        Спектр на сетке фона (истинная шкала), длина N.
    c : sequence of float, длина 3
        Коэффициенты отображения h(x) = c0 + c1*x + c2*x**2
        (канал образца x → канал фона y).
    n_out : int
        Число каналов выходного спектра образца.
    rng : numpy.random.Generator
        Генератор случайных чисел.
    thin : float, optional
        Коэффициент прореживания (по умолчанию 1.0).
    sigma_extra : float, optional
        Дополнительное уширение (по умолчанию 0.0).
    ref_ch : float или None, optional
        Справочный канал для нормировки уширения.

    Возврат
    -------
    ndarray
        Float-массив длины n_out — спектр образца.
    """
    # Шаг 1: целое число событий, опциональное прореживание
    n = np.rint(counts).astype(np.int64)
    n = np.maximum(n, 0)
    if thin < 1.0:
        n = rng.binomial(n, thin)

    # Шаг 2: позиции событий на шкале фона
    idx = np.repeat(np.arange(len(n), dtype=np.float64), n)
    y = idx - 0.5 + rng.random(len(idx))

    # Шаг 3: дополнительное уширение
    if sigma_extra > 0.0:
        ref = ref_ch if ref_ch is not None else len(n) / 2.0
        y = y + rng.normal(0.0, 1.0, len(y)) * sigma_extra * np.sqrt(np.maximum(y, 0.0) / ref)

    # Шаг 4: обратное отображение y → x
    x = inverse_h(c, y)

    # Шаг 5: бинаризация и гистограмма
    ix = np.floor(x + 0.5).astype(np.int64)
    mask = (ix >= 0) & (ix < n_out)
    ix = ix[mask]
    hist = np.bincount(ix, minlength=n_out)[:n_out]
    return hist.astype(np.float64)
