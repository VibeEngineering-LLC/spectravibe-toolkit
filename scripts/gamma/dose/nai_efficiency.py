# -*- coding: utf-8 -*-
"""Модельная кривая эффективности ε_tot(E) кристалла NaI(Tl) 2×2″.

Полная эффективность взаимодействия (доля фотонов, испытавших ЛЮБОЕ взаимодействие
на пути через кристалл толщиной d):
    ε_tot(E) = 1 − exp(−(μ/ρ)·ρ·d)
μ/ρ NaI — mixture-rule (Bragg) из NIST z11(Na)+z53(I) (data/nai_mu_rho_nist.csv,
w_Na=0.153373, w_I=0.846627; регенерация — data/build_nai_mu_rho.py).

ГРАНИЦА ТОЧНОСТИ: ε_tot — полная эфф. взаимодействия, НЕ фотопиковая (ε_ph). Точна
для «сколько фотонов провзаимодействовало», приблизительна как счётная эффективность
на комптон-континууме (частичный энерго-депозит). Строгий counts→флюенс требует полной
матрицы отклика (Monte-Carlo). Метод (2+ε) применяет ε_tot как первопорядковую поправку.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np

RHO_NAI = 3.667   # г/см³ — плотность NaI(Tl)
D_CM = 5.08       # см — толщина кристалла по оси (2 дюйма)
DEFAULT_CSV = Path(__file__).parent / "data" / "nai_mu_rho_nist.csv"


def load_nai_mu_rho(csv_path=DEFAULT_CSV):
    """Читает CSV μ/ρ NaI → (E_MeV[], mu_rho[]). Пропускает # и строку-заголовок."""
    E, mu = [], []
    for ln in open(csv_path, encoding="utf-8"):
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        p = ln.split(",")
        try:
            e, m = float(p[0]), float(p[1])
        except (ValueError, IndexError):
            continue
        E.append(e); mu.append(m)
    return np.array(E), np.array(mu)


_NAI_E, _NAI_MU = load_nai_mu_rho()


def eps_tot(E_keV, rho=RHO_NAI, d_cm=D_CM):
    """ε_tot = 1−exp(−(μ/ρ)·ρ·d). log-log интерполяция μ/ρ, клип к краям таблицы."""
    e = np.clip(np.asarray(E_keV, float) / 1000.0, _NAI_E[0], _NAI_E[-1])
    mu = np.exp(np.interp(np.log(e), np.log(_NAI_E), np.log(_NAI_MU)))
    return 1.0 - np.exp(-mu * rho * d_cm)


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    print("ε_tot NaI 2×2\" (ρ=%.3f г/см³, d=%.2f см) на реперах:" % (RHO_NAI, D_CM))
    for ek in (352, 609, 1120, 1461, 1764, 2614):
        print(f"  {ek:5.0f} кэВ : ε_tot = {float(eps_tot(ek)):.3f}")