"""binned.py — калибровка NAI: бинный подход, версия 2.

Почему версия 2
---------------
Версия 1 интерполировала накопленный счёт ОБРАЗЦА в дробных краях бинов фона.
Измерено (суррогат «тождество», Гамма-1С):
χ² в истинной точке 454.7, в точке d = 0.43 кан — 383.1, при d = ±1 — 531/557:
интерполяция шумного спектра сглаживает его в дробных краях и занижает χ²
вне целых сдвигов (рябь), смещение больше статистической ошибки.

Версия 2:
  • Бины задаются на сетке ОБРАЗЦА с целыми краями — образец никогда не
    интерполируется.
  • Интерполируется только шаблон ФОНА, предварительно сглаженный гауссианой
    σ_s ≥ 1 кан (шум на масштабе канала подавлен, интерполяция не меняет
    структуру).
  • Образец сглаживается тем же ядром на своей сетке (операция не зависит от
    h, разрешение обеих сторон одинаково).
  • Дисперсия суммы по бину вычисляется по НЕсглаженным счётам (сумма
    сглаженного по бину шире ядра ≈ сумма исходного).
  • Параметр s < 0 в версии 2 не поддерживается: трактуется как s = 0.
"""

import numpy as np
from scipy.ndimage import gaussian_filter1d
from gamma.calibration.intrinsic_shape.core import poly_eval, rebin_conserve_var
from .model import coefs, blur_var
from .surrogate import inverse_h


# ---------------------------------------------------------------------------
# Вспомогательные функции
# ---------------------------------------------------------------------------

def make_edges(ch_lo, ch_hi, fw_ch, bin_fwhm, min_w=1):
    """Построить края бинов по границам каналов образца.

    Края лежат на сетке образца (целые каналы ± 0.5). Ширина каждого бина
    пропорциональна ПШПВ в середине бина: w = bin_fwhm * FWHM(channel).

    Параметры
    ---------
    ch_lo, ch_hi : int
        Первый и последний канал образца (включительно).
    fw_ch : callable
        Функция канал → ПШПВ (в каналах).
    bin_fwhm : float
        Коэффициент: bin_width = bin_fwhm * FWHM.
    min_w : int
        Минимальная ширина бина (каналы).

    Возвращает
    ----------
    numpy.ndarray, shape (nbins + 1,)
    """
    e = ch_lo - 0.5
    edges = [e]
    while e < ch_hi + 0.5:
        w = max(min_w, int(round(bin_fwhm * fw_ch(e + 0.5))))
        e = min(e + w, ch_hi + 0.5)
        edges.append(e)
    return np.asarray(edges, dtype=np.float64)


def sum_bins(y, edges):
    """Суммировать каналы y по бинам с краями edges (на сетке образца)."""
    starts = (edges[:-1] + 0.5).astype(np.intp)
    y = np.asarray(y)[:int(edges[-1] + 0.5)]
    return np.add.reduceat(y, starts, axis=0)


def smooth_var_sigma(n, fw, kappa, s_min):
    """Массив σ_s по n каналам: max(s_min, kappa * FWHM(ch))."""
    return np.maximum(s_min, kappa * fw(np.arange(n, dtype=np.float64)))


def smooth(y, sig):
    """Сглаживание спектра с (возможно) переменным σ.

    Если все σ равны — гауссианский фильтр; иначе — blur_var (ядро
    переменной ширины).
    """
    if np.allclose(sig, sig[0]):
        return gaussian_filter1d(y, sig[0], mode="nearest")
    return blur_var(y, y, sig)[0]


# ---------------------------------------------------------------------------
# Основной класс
# ---------------------------------------------------------------------------

class BinnedProblem:
    """Бинная задача калибровки NAI (версия 2).

    Бины определены на сетке образца; шаблон фона интерполируется
    (после сглаживания) в краях бинов. Параметр s < 0 трактуется как s = 0.
    """

    def __init__(self, S, B, k, mask, h0, basis=None, ref_ch=None,
                 fw_ch=None, bin_fwhm=0.125, kappa=0.1, s_min=1.0,
                 min_w=2, margin=0.03):
        """
        Параметры
        ---------
        S : array, shape (n_s,)
            Спектр образца (float).
        B : array, shape (n_b,)
            Спектр фона (float).
        k : float
            Отношение времён наработки LT_s / LT_b.
        mask : array bool, shape (n_b,)
            Рабочий диапазон по каналам ФОНА.
        h0 : array, shape (3,)
            Начальные коэффициенты h (образец → фон).
        basis : array, shape (n_b, m), or None
            Линейный базис фона (на сетке фона).
        ref_ch : float, optional
            Опорный канал для нормировки s; по умолчанию n_b / 2.
        fw_ch : callable
            Канал фона → ПШПВ (каналы).
        bin_fwhm : float
        kappa : float
        s_min : float
        min_w : int
        margin : float
            Доля длины, на которую сузивается диапазон образца.
        """
        self.S = np.asarray(S, dtype=np.float64)
        self.B = np.asarray(B, dtype=np.float64)
        self.k = k
        self.mask = np.asarray(mask, dtype=bool)
        self.h0 = np.asarray(h0, dtype=np.float64)
        self.basis = basis
        self.fw_ch = fw_ch
        self.n_s = len(self.S)
        self.n_b = len(self.B)

        if ref_ch is None:
            ref_ch = self.n_b / 2.0
        self.ref_ch = ref_ch

        # --- Диапазон образца (по маске фона, обратное преобразование h0) ---
        idx = np.flatnonzero(self.mask)
        xa, xb = inverse_h(self.h0, [idx[0], idx[-1]])
        L = xb - xa
        x_lo = int(np.ceil(xa + margin * L))
        x_hi = int(np.floor(xb - margin * L))

        # --- ПШПВ в каналах образца ---
        def fw_s(x):
            return self.fw_ch(poly_eval(self.h0, x)) / self.h0[1]

        # --- Бины на сетке образца ---
        self.edges = make_edges(x_lo, x_hi, fw_s, bin_fwhm, min_w)
        self.nbins = len(self.edges) - 1
        self.nmask = self.nbins

        # --- Образец: сглаживание + биннинг ---
        sig_s = smooth_var_sigma(self.n_s, fw_s, kappa, s_min)
        self.Ss = smooth(self.S, sig_s)
        self.Sb = sum_bins(self.Ss, self.edges)
        # Дисперсия бинов — по НЕсглаженным счётам
        self.Sv = np.maximum(sum_bins(self.S, self.edges), 1.0)

        # --- Фон: сглаживание + cumulative ---
        self.sig_b = smooth_var_sigma(self.n_b, self.fw_ch, kappa, s_min)
        self.Bs = smooth(self.B, self.sig_b)
        self.cumBs = np.concatenate([[0.0], np.cumsum(self.Bs)])
        self.cumB = np.concatenate([[0.0], np.cumsum(self.B)])

        # --- Базис ---
        self.m = 0
        if basis is not None:
            self.m = basis.shape[1]
            self.cumBasis = np.vstack([np.zeros((1, self.m)),
                                       np.cumsum(basis, axis=0)])

    # ------------------------------------------------------------------
    # Внутренние методы
    # ------------------------------------------------------------------

    def _interp_cum(self, cum, ye):
        """Интерполяция cumulative-массива в краях ye, возврат разностей.

        cum : 1D (n_b+1,) или 2D (n_b+1, m)
        ye  : 1D (nbins+1,) — края бинов в координатах ФОНА
        """
        x = np.arange(self.n_b + 1, dtype=np.float64) - 0.5
        if cum.ndim == 1:
            val = np.interp(ye, x, cum)
        else:
            val = np.empty((len(ye), cum.shape[1]), dtype=np.float64)
            for j in range(cum.shape[1]):
                val[:, j] = np.interp(ye, x, cum[:, j])
        return np.diff(val, axis=0)

    # ------------------------------------------------------------------
    # Оценка
    # ------------------------------------------------------------------

    def evaluate(self, p):
        """Вычислить бинную оценку для параметров p = [d, gN, q, s].

        Возвращает dict или None (если края не строго возрастают).
        Параметр s < 0 трактуется как s = 0.
        """
        # 1. Коэффициенты и края в координатах фона
        c = coefs(p, self.h0, self.n_b)
        ye = poly_eval(c, self.edges)
        if not np.all(np.diff(ye) > 0.0):
            return None

        # 2. Шаблон фона (сглаженный + интерполированный)
        s = p[3]
        if s > 0.0:
            sig = abs(s) * np.sqrt(np.arange(self.n_b) / self.ref_ch)
            bb, _ = blur_var(self.B, self.B, sig)
            cumT = np.concatenate([[0.0], np.cumsum(smooth(bb, self.sig_b))])
        else:
            # s < 0 не поддерживается в v2 — трактуем как s = 0
            cumT = self.cumBs

        # 3. Интерполяция шаблона и его дисперсии в бины
        Tb = self._interp_cum(cumT, ye)
        Tv = self._interp_cum(self.cumB, ye)

        # 4. Нейтральная разность и веса
        r0 = self.Sb - self.k * Tb
        var = np.maximum(self.Sv + self.k * self.k * Tv, 1.0)
        w = 1.0 / np.sqrt(var)

        # 5. Линейный базис
        if self.m > 0:
            Ab = self._interp_cum(self.cumBasis, ye)  # (nbins, m)
            amps, *_ = np.linalg.lstsq(Ab * w[:, None], r0 * w, rcond=None)
            fit_lin = Ab @ amps
        else:
            amps = np.zeros(0)
            fit_lin = np.zeros(self.nbins)

        # 6. Остатки
        z = (r0 - fit_lin) * w

        return {
            "c": c,
            "Sb": self.Sb,
            "Tb": Tb,
            "net": r0,
            "var": var,
            "amps": amps,
            "fit_lin": fit_lin,
            "z": z,
        }

    def resid(self, p):
        """Вектор стандартизированных остатков (или placeholder при ошибке)."""
        ev = self.evaluate(p)
        if ev is None:
            return np.full(self.nbins, 1e3)
        return ev["z"]

    def chi2(self, p):
        """χ² = сумма квадратов остатков."""
        z = self.resid(p)
        return float(np.sum(z * z))

    def ndf(self, n_free):
        """Число степеней свободы: nbins − n_free − m."""
        return self.nbins - n_free - self.m

    # ------------------------------------------------------------------
    # Поканальный вывод (на сетке фона, НЕсглаженные спектры)
    # ------------------------------------------------------------------

    def evaluate_channels(self, p, amps=None):
        """Поканальная разность на сетке ФОНА для визуализации.

        Использует НЕСглаженные спектры. Параметр s < 0 → s = 0.
        """
        c = coefs(p, self.h0, self.n_b)
        # Образец → сетка фона (ребиннинг с сохранением дисперсии)
        sh, vs = rebin_conserve_var(
            self.S, poly_eval(c, np.arange(self.n_s + 1, dtype=np.float64) - 0.5),
            self.n_b,
        )

        # Фон
        bb = self.B.copy()
        vb = self.B.copy()
        s = p[3]
        if s > 0.0:
            sig = abs(s) * np.sqrt(np.arange(self.n_b) / self.ref_ch)
            bb, vb = blur_var(bb, vb, sig)

        net = sh - self.k * bb
        var = np.maximum(vs + self.k * self.k * vb, 1.0)

        # Базис
        if self.m > 0:
            if amps is None:
                ev = self.evaluate(p)
                amps = ev["amps"] if ev is not None else np.zeros(self.m)
            fit_lin = self.basis @ amps
        else:
            amps = np.zeros(0)
            fit_lin = np.zeros(self.n_b)

        z = (net - fit_lin) / np.sqrt(var)

        return {
            "c": c,
            "sh": sh,
            "vs": vs,
            "bb": bb,
            "vb": vb,
            "net": net,
            "var": var,
            "amps": amps,
            "fit_lin": fit_lin,
            "z": z,
        }
