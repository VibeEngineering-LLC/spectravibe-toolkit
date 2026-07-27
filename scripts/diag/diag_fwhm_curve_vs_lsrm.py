# -*- coding: utf-8 -*-
"""Proof that my pipeline's FWHM(E) curve (source='default_NaI_63x63') does NOT
match this spectrum. Compares the default-NaI model curve at each LSRM line vs
the LSRM 'Параметры пиков' measured FWHM, then fits a corrected curve to the
measured points. Operator 2026-06-20: 'калибровку по FWHM нужно проверять,
слепо на кривую нельзя полагаться'."""
import numpy as np

# My pipeline's current FWHM model (from report calibration.fwhm_cal):
A, B, C = 0.0, 2.950048, 0.0005764  # FWHM^2 = A + B*E + C*E^2  (default_NaI_63x63)
def fwhm_model(E):
    return float(np.sqrt(A + B*E + C*E*E))

# LSRM Gamma-1S ground truth, 'Параметры пиков' (operator image 2026-06-20):
# (energy_keV, measured_FWHM_keV, measured_FWHM_pct)
LSRM = [
    (208.129, 21.644, 10.399), (237.508, 24.029, 10.117),
    (327.625, 29.292, 8.941), (337.944, 29.881, 8.842),
    (459.984, 34.666, 7.536), (507.612, 36.835, 7.255),
    (580.077, 43.116, 7.433), (900.228, 55.823, 6.201),
    (953.433, 58.215, 6.104), (957.978, 58.402, 6.096),
    (1581.207, 90.683, 5.735), (1613.500, 92.099, 5.702),
    (1623.624, 92.419, 5.692), (2612.857, 112.796, 4.317),
]

print('E_keV   LSRM_FWHM  MY_model  Δ(my-lsrm)  Δ%')
es = np.array([e for e,_,_ in LSRM]); fw = np.array([f for _,f,_ in LSRM])
for e, f, _ in LSRM:
    m = fwhm_model(e); d = m - f
    print('%8.1f  %8.2f  %8.2f  %+8.2f  %+6.1f%%' % (e, f, m, d, 100*d/f))

# Fit corrected curve FWHM^2 = a + b*E + c*E^2 to the LSRM measured points
y = fw**2
M = np.vstack([np.ones_like(es), es, es*es]).T
coef, *_ = np.linalg.lstsq(M, y, rcond=None)
print('\nCorrected FWHM^2 = a + b*E + c*E^2  fitted to LSRM measured:')
print('  a=%.6g  b=%.6g  c=%.6g' % (coef[0], coef[1], coef[2]))
fit = np.sqrt(np.clip(M @ coef, 0, None))
print('  max |residual| = %.3f keV   FWHM@661 = %.3f keV'
      % (np.max(np.abs(fit - fw)), np.sqrt(coef[0]+coef[1]*661+coef[2]*661**2)))
print('  (my default curve FWHM@661 = %.3f keV)' % fwhm_model(661))
