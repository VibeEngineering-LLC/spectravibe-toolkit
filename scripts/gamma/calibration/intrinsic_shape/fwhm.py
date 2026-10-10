"""Модель ПШПВ для BecqMoni: ПШПВ[кан] = sqrt(c0 + c1*ch + c2*ch^2) (блок SqrtFwhmCalibration).
Условия BecqMoni (SqrtFwhmCalibration.CheckCalibration): F^2 не убывает на 0..N-1, F(0)>0; здесь ещё F(0)>=1 кан.
Перенесено из прототипа labr3_selfcal (fwhm_model.py, fwhm_xml.py, check_fwhm_xml.py, fit_device_defaults.py), 2026-10-10.
Исправление по второй проверке: вершина параболы стояла ровно на N-1 (нулевой запас) -- теперь ширина не убывает
до канала monotone_to*N (по умолчанию 1.5*N), либо c2 >= 0 (режим "c2>=0")."""
import re
import numpy as np
from scipy.optimize import minimize, least_squares

from .core import poly_eval

# Точки тестовой пары LaBr3 (канал фона, ПШПВ в каналах): Lu-176 202, 307, сумма 202+307, ядро сумм-пика La-138, пик ~37 кэВ
FW_X = np.array([29.20, 142.87, 211.97, 338.14, 847.10])
FW_F = np.array([2.47, 6.14, 7.29, 8.62, 11.63])


def fwhm(c, x):
    x = np.asarray(x, float)
    return np.sqrt(np.maximum(c[0] + c[1] * x + c[2] * x * x, 0.0))


def fit_fwhm(n_ch, x=FW_X, f=FW_F, monotone_to=1.5):
    """monotone_to: число -- F^2 не убывает до канала monotone_to*n_ch; "c2>=0" -- выпуклая парабола (c2 >= 0)."""
    x, f = np.asarray(x, float), np.asarray(f, float)
    if monotone_to == "c2>=0":
        mono = {"type": "ineq", "fun": lambda c: c[2]}
    else:
        x_top = float(monotone_to) * n_ch - 1
        mono = {"type": "ineq", "fun": lambda c: c[1] + 2 * c[2] * x_top}
    cons = [mono, {"type": "ineq", "fun": lambda c: c[0] - 1.0}]
    best = None
    for s0 in ([1, .2, -5e-5], [3, .2, -3e-5], [1.5, .3, -1e-4], [1, .1, 0.0]):
        r = minimize(lambda c: np.sum(((fwhm(c, x) - f) / f) ** 2), s0, constraints=cons, method="SLSQP", options={"maxiter": 500, "ftol": 1e-14})
        if r.success and (best is None or r.fun < best.fun):
            best = r
    pts = ", ".join("%.1f:%.2f" % (a, b) for a, b in zip(x, f))
    hint = ("проверьте точки: ПШПВ должна расти с каналом; уберите точки, где подгонка не сошлась "
            "(в разностном спектре -- нет пика), или смягчите ограничение --monotone-to c2>=0")
    if best is None:
        raise ValueError("модель ПШПВ не подобрана ни из одного старта; точки (канал:ПШПВ) %s; %s" % (pts, hint))
    c = [float(v) for v in best.x]
    g = np.arange(n_ch)
    d = np.diff(c[0] + c[1] * g + c[2] * g * g)
    if not (np.all(d >= 0) and c[0] > 0):
        raise ValueError("модель ПШПВ c=%s нарушает CheckCalibration: F^2 убывает в %d каналах из %d, F(0)^2=%.4g; "
                         "точки (канал:ПШПВ) %s; %s" % (["%.6g" % v for v in c], int((d < 0).sum()), n_ch - 1, c[0], pts, hint))
    return c


# --- перенесено из прототипа fwhm_xml.py (блок SqrtFwhmCalibration) ---


def fwhm_block(nl, coefs, peaks, ecal):
    """peaks: список (канал, ПШПВ в каналах); coefs: c0, c1, c2; ecal: коэффициенты шкалы E(ch) по возрастанию.
    Channel в формате -- целое (все образцы BecqMoni/AtomSpectra; читатель toolkit берёт int), поэтому Energy
    пересчитывается из ЦЕЛОГО канала (вторая проверка: целый Channel был не согласован с дробной Energy)."""
    def tag(ind, name, val):
        return " " * ind + "<%s>%s</%s>" % (name, val, name) + nl
    rows = [(int(round(x)), float(f)) for x, f in peaks]
    pk = "".join(" " * 10 + "<Peak>" + nl + tag(12, "Channel", ch) + tag(12, "Energy", "%.17g" % poly_eval(ecal, ch))
                 + tag(12, "FWHM", "%.17g" % f) + " " * 10 + "</Peak>" + nl for ch, f in rows)
    cf = "".join(tag(10, "Coefficient", "%.17g" % c) for c in coefs)
    return (" " * 6 + "<SqrtFwhmCalibration>" + nl + " " * 8 + "<CalibrationPeaks>" + nl + pk + " " * 8 + "</CalibrationPeaks>" + nl
            + " " * 8 + "<Coefficients>" + nl + cf + " " * 8 + "</Coefficients>" + nl + tag(8, "PeakType", 0)
            + tag(8, "ExpGaussExpLeftTail", 1) + tag(8, "ExpGaussExpRightTail", 1) + tag(8, "Chi2pNdp", -1)
            + " " * 6 + "</SqrtFwhmCalibration>")


# --- перенесено из прототипа check_fwhm_xml.py:3-7 (повтор SqrtFwhmCalibration.CheckCalibration) ---
def check(c, n=2048):
    """(годен, ПШПВ в каналах 0/500/1309/n-1, вершина параболы -c1/(2 c2) или inf, запас вершины в каналах за n-1)."""
    f = lambda x: max(c[0] + c[1] * x + c[2] * x * x, 0.0) ** 0.5
    mono = all(f(i - 1) <= f(i) for i in range(1, n))
    vertex = -c[1] / (2 * c[2]) if c[2] < 0 else float("inf")
    return mono and f(0) > 0, [round(f(k), 3) for k in (0, 500, 1309, n - 1)], vertex, vertex - (n - 1)


# --- перенесено из прототипа fit_device_defaults.py (умолчание прибора BecqMoni, FwhmCalibration.cs:92-121) ---
def device_defaults(x=FW_X, f=FW_F, ch_ref=1309):
    """FWHM_AT_0 / Ch_Fwhm / Width_Fwhm: прямая ПШПВ^2 через (0, F0^2) и (Ch, W^2), модель SimpleSqrt."""
    x, f = np.asarray(x, float), np.asarray(f, float)
    r = least_squares(lambda c: (np.sqrt(np.maximum(c[0] + c[1] * x, 0)) - f) / f, [1, .2], bounds=([1e-3, 0], [np.inf, np.inf]))
    F0 = round(float(np.sqrt(r.x[0])), 2)
    W = round(float(np.sqrt(r.x[0] + r.x[1] * ch_ref)), 2)
    return F0, int(ch_ref), W, [F0 ** 2, (W ** 2 - F0 ** 2) / int(ch_ref), 0.0]


def insert_block(xml_text, block):
    """Вставить/заменить SqrtFwhmCalibration в первом ResultData (перед </ResultData>)."""
    t = re.sub(r"[ \t]*<SqrtFwhmCalibration>.*?</SqrtFwhmCalibration>\r?\n", "", xml_text, count=1, flags=re.S)
    nl = "\r\n" if "\r\n" in t else "\n"
    return re.sub(r"(\r?\n)([ \t]*</ResultData>)", lambda m: m.group(1) + block + nl + m.group(2), t, count=1)
