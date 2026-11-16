# -*- coding: utf-8 -*-
"""Критерий 3б v3: выбор окон проверки фона по составу образца
(audit/nai_runs/criteria_v3.md) и проверка чувствительности к известной ошибке отображения.
ИТЕРАЦИЯ 4: вне пути приёмки — только диагностика сильных образцов; результат несёт counted=False и не идёт в зачёт."""

import json
import pathlib

import numpy as np

from gamma.calibration.intrinsic_shape.core import poly_deriv
from . import bgcal, binned, checks

# ── Константы ──────────────────────────────────────────────────────────────────

GROUPS = {
    "K-40":   ["K-40"],
    "Cs-137": ["Cs-137", "Ba-137m"],
    "Ra-226": ["Ra-226", "Pb-214", "Bi-214"],
    "Th-232": ["Ac-228", "Pb-212", "Bi-212", "Tl-208"],
    "Th-228": ["Pb-212", "Bi-212", "Tl-208"],
}

POSITRON = {"Na-22", "Zn-65"}

SUMS = {
    "Co-60": [2505.7],
    "Na-22": [1785.5],
    "Y-88":  [2734.1],
}

M_E = 511.0

DATA = (
    pathlib.Path(__file__).resolve().parents[4] / "data" / "nuclides.json"
)

# ── Вспомогательные функции ───────────────────────────────────────────────────


def sample_lines(nuclides, imin=1.0, emin=60.0, emax=3000.0, data_path=DATA):
    """Собрать список линий (E, I) для данного набора нуклидов."""
    with open(data_path, encoding="utf-8") as f:
        D = json.load(f)

    raw = []
    for nuc in nuclides:
        members = GROUPS.get(nuc, [nuc])
        for m in members:
            if m not in D:
                continue
            for ln in D[m]["lines"]:
                E, I = ln[0], ln[1]
                if I >= imin and emin <= E <= emax:
                    raw.append((E, I))
        # аннигиляционное излучение
        if nuc in POSITRON:
            raw.append((M_E, 100.0))
        # линии суммирования
        if nuc in SUMS:
            for E_sum in SUMS[nuc]:
                raw.append((E_sum, 100.0))

    # сортировка и удаление повторов (|ΔE| < 0.3 кэВ)
    raw.sort(key=lambda p: p[0])
    out = []
    for E, I in raw:
        if out and abs(E - out[-1][0]) < 0.3:
            continue
        out.append((E, I))
    return out


def zones(lines, icomp=5.0):
    """Построить список энергий зон (основные + комптон-пик/край для интенсивных линий)."""
    result = []
    for E, I in lines:
        result.append(float(E))
        if I >= icomp:
            Ec = 2.0 * E * E / (M_E + 2.0 * E)
            Eb = E * M_E / (M_E + 2.0 * E)
            result.append(float(Ec))
            result.append(float(Eb))
    return sorted(result)


def select_windows(windows_kev, zone_kev, fwhm_kev):
    """Исключить окна, попадающие в зоны линий образца.

    Возвращает (kept, dropped) — списки энергий окон (кэВ).
    """
    kept = []
    dropped = []
    for W in windows_kev:
        excluded = False
        fwhm_W = fwhm_kev(W)
        for X in zone_kev:
            if abs(W - X) < fwhm_W + 2.0 * fwhm_kev(X):
                excluded = True
                break
        (dropped if excluded else kept).append(W)
    return kept, dropped


def cal_from_result(r):
    """Извлечь словарь калибровки из результата run_auto."""
    return {
        "ecal":  r["ecal_bg"],
        "fwc":   r["fwc"],
        "R662":  r["R662"],
        "n":     len(r["arrays"]["net"]),
    }


def check(r, bg, sm, profile, kept, c_override=None):
    """Прогнать checks.window_tests по сохранным окнам.

    Если c_override задан — пересчитать сетку фона с модифицированными
    коэффициентами отображения.
    """
    cal = cal_from_result(r)
    fw = bgcal.fwhm_ch_fn(cal)
    e2ch = lambda E: float(bgcal.energy_to_ch(cal, E))
    fwhm_kev = lambda E: float(fw(e2ch(E)) * poly_deriv(cal["ecal"], e2ch(E)))

    if c_override is None:
        net = r["arrays"]["net"]
        var = r["arrays"]["var"]
    else:
        B = bg["counts"].astype(float).copy()
        B[-1] = 0
        S = sm["counts"].astype(float).copy()
        S[-1] = 0
        prob = binned.BinnedProblem(
            S, B, r["k"], r["arrays"]["mask"],
            list(c_override), None,
            ref_ch=len(B) / 2,
            fw_ch=fw,
            bin_fwhm=profile.get("bin_fwhm", 0.5),
        )
        ev = prob.evaluate_channels(np.zeros(4))
        net = ev["net"]
        var = ev["var"]

    tests = checks.window_tests(
        net, var, e2ch, fwhm_kev, kept,
        1.0, (), 1.5, True, None,
    )
    return tests


def evaluate_pair(r, bg, sm, profile, nuclides,
                  eps_list=(0.005, 0.01), z_thr=3.0, min_windows=2):
    """Полная оценка пары (r, bg/sm) по критерию 3б v3.

    1. Выбор окон по составу образца.
    2. Проверка z_dip в сохранных окнах.
    3. Проверка чувствительности к известной ошибке (масштабирование c[1]).
    """
    # --- подготовка калибровки и линий ---
    cal = cal_from_result(r)
    fw = bgcal.fwhm_ch_fn(cal)
    e2ch = lambda E: float(bgcal.energy_to_ch(cal, E))
    fwhm_kev = lambda E: float(fw(e2ch(E)) * poly_deriv(cal["ecal"], e2ch(E)))

    lines = sample_lines(nuclides)
    zk = zones(lines)

    all_windows = [w["E"] for w in profile["null_windows"]]
    kept, dropped = select_windows(all_windows, zk, fwhm_kev)

    # --- ранний выход, если окон мало ---
    if len(kept) < min_windows:
        return {
            "counted": False,
            "flag":   "шкала не проверена",
            "status": "нечем проверить",
            "kept":   kept,
            "dropped": dropped,
            "tests":  [],
            "sens":   {},
        }

    # --- проверка фона в сохранных окнах ---
    tests = check(r, bg, sm, profile, kept)
    ok = all(abs(t["z_dip"]) <= z_thr for t in tests)

    # --- чувствительность к ошибке отображения ---
    sens = {}
    for eps in eps_list:
        c = list(r["c"])
        c[1] = c[1] * (1.0 + eps)
        tb = check(r, bg, sm, profile, kept, c)
        sens[eps] = {
            "caught":        any(abs(t["z_dip"]) > z_thr for t in tb),
            "max_abs_z_dip": max(abs(t["z_dip"]) for t in tb),
        }

    return {
        "counted":       False,   # итерация 4: диагностика вне зачёта
        "flag":          "шкала не проверена",
        "status":        "прошла" if ok else "не прошла",
        "kept":          kept,
        "dropped":       dropped,
        "tests":         tests,
        "sens":          sens,
        "max_abs_z_dip": max(abs(t["z_dip"]) for t in tests),
        "z_dip":         [t["z_dip"] for t in tests],
    }
