# -*- coding: utf-8 -*-
"""Подскилл расчёта мощности дозы (МД) из гамма-спектра — gamma.dose.

Три метода (форма кривой Σ=1, привязка к внешнему якорю MD_TOTAL [мкЗв/ч]):
  (3)   ЛСРМ:      w = f10(E)·(μen/ρ)_air·E   — без ε детектора
  (2)   ICRP-74:   w = h_Φ(E)                  — counts≈флюенс, без ε
  (2+ε) ICRP-74:   w = h_Φ(E)/ε_tot(E)         — counts→флюенс через ε NaI

Провенанс данных/формул и границы точности — README.md.
"""
from .icrp74_hphi import h_phi_interp
from .nai_efficiency import eps_tot, load_nai_mu_rho, RHO_NAI, D_CM
from .lsrm_coeffs import (
    parse_spectrum, muen_air_interp, f10_interp, E_LO, E_HI, DEFAULT_SPECTRUM,
)

__all__ = [
    "h_phi_interp", "eps_tot", "load_nai_mu_rho", "RHO_NAI", "D_CM",
    "parse_spectrum", "muen_air_interp", "f10_interp", "E_LO", "E_HI",
    "DEFAULT_SPECTRUM",
]