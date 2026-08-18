# -*- coding: utf-8 -*-
"""ICRP-74 Table A.21 — фотонные конверсионные коэффициенты (операционные величины).

ПЕРВОИСТОЧНИК: ICRP Publication 74 (Annals of the ICRP vol.26 No.3-4, 1996),
  "Conversion Coefficients for use in Radiological Protection against External Radiation".
  Table A.21, внутр. стр.179 = PDF-стр.194 файла  ANIB_26_3-4.pdf
  Заголовок: "Conversion coefficients for the ambient dose equivalent, H*(10), and
  directional dose equivalent, H'(0.07,0 deg), from photon fluence and air kerma free-in-air."
  Footnote a: "Data compiled from ICRU Report 47 (1992a) using Hubbell and Seltzer (1995)."

ИЗВЛЕЧЕНИЕ: baidu/Unlimited-OCR (VLM, "китайский OCR"), скан стр.194 @200 dpi.
  Провенанс OCR: dose_curve/icrp74_data/ocr_A21/page_194/result.md.
  §23-СВЕРКА: все 25 строк колонки H*(10)/Phi сверены вручную против
  dose_curve/icrp74_data/verify_p194.jpg (@175 dpi) — 25/25 совпали, 0 переворотов. 2026-07-06.

ЕДИНИЦЫ:
  E                     — энергия фотона, MeV
  hstar10_per_Ka        — H*(10)/Ka, Sv/Gy
  hprime007_per_Ka      — H'(0.07,0deg)/Ka, Sv/Gy
  Ka_per_fluence        — Ka/Phi, pGy·cm²
  hstar10_per_fluence   — H*(10)/Phi, pSv·cm²  <-- ОПЕРАЦИОННАЯ h_Phi(E), ЦЕЛЬ
  hprime007_per_fluence — H'(0.07,0deg)/Phi, pSv·cm²
"""
import numpy as np

# E, H*(10)/Ka, H'(0.07)/Ka, Ka/Phi[pGy·cm²], H*(10)/Phi[pSv·cm²], H'(0.07)/Phi[pSv·cm²]
_A21 = np.array([
    [0.010, 0.008, 0.95, 7.60,  0.061, 7.20],
    [0.015, 0.26,  0.99, 3.21,  0.83,  3.19],
    [0.020, 0.61,  1.05, 1.73,  1.05,  1.81],
    [0.030, 1.10,  1.22, 0.739, 0.81,  0.90],
    [0.040, 1.47,  1.41, 0.438, 0.64,  0.62],
    [0.050, 1.67,  1.53, 0.328, 0.55,  0.50],
    [0.060, 1.74,  1.59, 0.292, 0.51,  0.47],
    [0.080, 1.72,  1.61, 0.308, 0.53,  0.49],
    [0.100, 1.65,  1.55, 0.372, 0.61,  0.58],
    [0.150, 1.49,  1.42, 0.600, 0.89,  0.85],
    [0.200, 1.40,  1.34, 0.856, 1.20,  1.15],
    [0.300, 1.31,  1.31, 1.38,  1.80,  1.80],
    [0.400, 1.26,  1.26, 1.89,  2.38,  2.38],
    [0.500, 1.23,  1.23, 2.38,  2.93,  2.93],
    [0.600, 1.21,  1.21, 2.84,  3.44,  3.44],
    [0.800, 1.19,  1.19, 3.69,  4.38,  4.38],
    [1.000, 1.17,  1.17, 4.47,  5.20,  5.20],
    [1.500, 1.15,  1.15, 6.12,  6.90,  6.90],
    [2.000, 1.14,  1.14, 7.51,  8.60,  8.60],
    [3.000, 1.13,  1.13, 9.89,  11.1,  11.1],
    [4.000, 1.12,  1.12, 12.0,  13.4,  13.4],
    [5.000, 1.11,  1.11, 13.9,  15.5,  15.5],
    [6.000, 1.11,  1.11, 15.8,  17.6,  17.6],
    [8.000, 1.11,  1.11, 19.5,  21.6,  21.6],
    [10.00, 1.10,  1.10, 23.2,  25.6,  25.6],
])

E_MEV                 = _A21[:, 0]
HSTAR10_PER_KA        = _A21[:, 1]
HPRIME007_PER_KA      = _A21[:, 2]
KA_PER_FLUENCE        = _A21[:, 3]   # pGy·cm²
HSTAR10_PER_FLUENCE   = _A21[:, 4]   # pSv·cm²  <-- h_Phi(E)
HPRIME007_PER_FLUENCE = _A21[:, 5]   # pSv·cm²


def h_phi_interp(E_keV):
    """H*(10)/Phi [pSv·cm²] интерполяцией по log(E). Вне 10 keV–10 MeV — клип к краю."""
    E_mev = np.asarray(E_keV, float) / 1000.0
    E_mev = np.clip(E_mev, E_MEV[0], E_MEV[-1])
    return np.interp(np.log(E_mev), np.log(E_MEV), HSTAR10_PER_FLUENCE)


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    print("ICRP-74 Table A.21 — H*(10)/Phi [pSv·cm²], %d rows, %.3f-%.1f MeV"
          % (len(E_MEV), E_MEV[0], E_MEV[-1]))
    for e, h in zip(E_MEV, HSTAR10_PER_FLUENCE):
        print(f"  {e:7.3f} MeV : {h:6.3f}")
