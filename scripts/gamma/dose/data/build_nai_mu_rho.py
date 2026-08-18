# -*- coding: utf-8 -*-
r"""
build_nai_mu_rho.py
===================
Строит массовый коэффициент ослабления mu/rho и массовый коэффициент поглощения
энергии mu_en/rho для NaI (натрий йодистый, sodium iodide) по ПРАВИЛУ СМЕСИ из
ЭЛЕМЕНТНЫХ таблиц NIST (X-Ray Mass Attenuation Coefficients), т.к. составная
страница NIST ComTab/sodiumiodide.html стабильно отдаёт HTTP 404.

Источник (дословно, извлечено WebFetch 2026-07-06):
  Na (Z=11): https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z11.html
  I  (Z=53): https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z53.html

Правило смеси (Bragg additivity, как у самого NIST для соединений):
  (mu/rho)_NaI    = w_Na*(mu/rho)_Na       + w_I*(mu/rho)_I
  (mu_en/rho)_NaI = w_Na*(mu_en/rho)_Na    + w_I*(mu_en/rho)_I

Массовые доли (атомные массы NIST/IUPAC):
  M(Na) = 22.98977, M(I) = 126.90447, M(NaI) = 149.89424
  w_Na = 0.153372, w_I = 0.846628

Сетка энергий = объединение узлов обеих таблиц. Основная сетка Na и I совпадает
узел-в-узел (1 keV ... 20 MeV) -> там правило смеси применяется ДОСЛОВНО. Йод
добавляет узлы своих краёв (L: 4.55710e-3, 4.70229e-3, 4.85210e-3, 5.18810e-3;
K: 3.31694e-2 ~33.17 keV, дубль до/после сохранён), которых нет у Na. В этих узлах
значение I берётся ДОСЛОВНО, а Na (гладкое там) получается ЛОГ-ЛОГ интерполяцией.
"""
import math
import os
import sys

# Na (Z=11) -- z11.html, извлечено WebFetch 2026-07-06:
NA = [
    (1.00000e-03, 6.542e+02, 6.522e+02),
    (1.03542e-03, 5.960e+02, 5.941e+02),
    (1.07210e-03, 5.429e+02, 5.410e+02),
    (1.07210e-03, 6.435e+03, 6.320e+03),
    (1.50000e-03, 3.194e+03, 3.151e+03),
    (2.00000e-03, 1.521e+03, 1.504e+03),
    (3.00000e-03, 5.070e+02, 5.023e+02),
    (4.00000e-03, 2.261e+02, 2.238e+02),
    (5.00000e-03, 1.194e+02, 1.178e+02),
    (6.00000e-03, 7.030e+01, 6.915e+01),
    (8.00000e-03, 3.018e+01, 2.941e+01),
    (1.00000e-02, 1.557e+01, 1.499e+01),
    (1.50000e-02, 4.694e+00, 4.313e+00),
    (2.00000e-02, 2.057e+00, 1.759e+00),
    (3.00000e-02, 7.197e-01, 4.928e-01),
    (4.00000e-02, 3.969e-01, 2.031e-01),
    (5.00000e-02, 2.804e-01, 1.063e-01),
    (6.00000e-02, 2.268e-01, 6.625e-02),
    (8.00000e-02, 1.796e-01, 3.761e-02),
    (1.00000e-01, 1.585e-01, 2.931e-02),
    (1.50000e-01, 1.335e-01, 2.579e-02),
    (2.00000e-01, 1.199e-01, 2.635e-02),
    (3.00000e-01, 1.029e-01, 2.771e-02),
    (4.00000e-01, 9.185e-02, 2.833e-02),
    (5.00000e-01, 8.372e-02, 2.845e-02),
    (6.00000e-01, 7.736e-02, 2.830e-02),
    (8.00000e-01, 6.788e-02, 2.760e-02),
    (1.00000e+00, 6.100e-02, 2.669e-02),
    (1.25000e+00, 5.454e-02, 2.549e-02),
    (1.50000e+00, 4.968e-02, 2.437e-02),
    (2.00000e+00, 4.282e-02, 2.249e-02),
    (3.00000e+00, 3.487e-02, 1.997e-02),
    (4.00000e+00, 3.037e-02, 1.842e-02),
    (5.00000e+00, 2.753e-02, 1.743e-02),
    (6.00000e+00, 2.559e-02, 1.675e-02),
    (8.00000e+00, 2.319e-02, 1.595e-02),
    (1.00000e+01, 2.181e-02, 1.552e-02),
    (1.50000e+01, 2.023e-02, 1.508e-02),
    (2.00000e+01, 1.970e-02, 1.496e-02),
]

# I (Z=53) -- z53.html, извлечено WebFetch 2026-07-06:
I = [
    (1.00000e-03, 9.096e+03, 9.078e+03),
    (1.03542e-03, 8.465e+03, 8.448e+03),
    (1.07210e-03, 7.863e+03, 7.847e+03),
    (1.07210e-03, 8.198e+03, 8.181e+03),
    (1.50000e-03, 3.919e+03, 3.908e+03),
    (2.00000e-03, 1.997e+03, 1.988e+03),
    (3.00000e-03, 7.420e+02, 7.351e+02),
    (4.00000e-03, 3.607e+02, 3.548e+02),
    (4.55710e-03, 2.592e+02, 2.537e+02),
    (4.55710e-03, 7.550e+02, 7.121e+02),
    (4.70229e-03, 7.123e+02, 6.724e+02),
    (4.85210e-03, 6.636e+02, 6.270e+02),
    (4.85210e-03, 8.943e+02, 8.375e+02),
    (5.00000e-03, 8.430e+02, 7.903e+02),
    (5.18810e-03, 7.665e+02, 7.198e+02),
    (5.18810e-03, 8.837e+02, 8.283e+02),
    (6.00000e-03, 6.173e+02, 5.822e+02),
    (8.00000e-03, 2.922e+02, 2.777e+02),
    (1.00000e-02, 1.626e+02, 1.548e+02),
    (1.50000e-02, 5.512e+01, 5.208e+01),
    (2.00000e-02, 2.543e+01, 2.363e+01),
    (3.00000e-02, 8.561e+00, 7.622e+00),
    (3.31694e-02, 6.553e+00, 5.744e+00),
    (3.31694e-02, 3.582e+01, 1.188e+01),
    (4.00000e-02, 2.210e+01, 9.616e+00),
    (5.00000e-02, 1.232e+01, 6.573e+00),
    (6.00000e-02, 7.579e+00, 4.518e+00),
    (8.00000e-02, 3.510e+00, 2.331e+00),
    (1.00000e-01, 1.942e+00, 1.342e+00),
    (1.50000e-01, 6.978e-01, 4.742e-01),
    (2.00000e-01, 3.663e-01, 2.295e-01),
    (3.00000e-01, 1.771e-01, 9.257e-02),
    (4.00000e-01, 1.217e-01, 5.650e-02),
    (5.00000e-01, 9.701e-02, 4.267e-02),
    (6.00000e-01, 8.313e-02, 3.598e-02),
    (8.00000e-01, 6.749e-02, 2.962e-02),
    (1.00000e+00, 5.841e-02, 2.646e-02),
    (1.25000e+00, 5.111e-02, 2.399e-02),
    (1.50000e+00, 4.647e-02, 2.243e-02),
    (2.00000e+00, 4.124e-02, 2.092e-02),
    (3.00000e+00, 3.716e-02, 2.059e-02),
    (4.00000e+00, 3.607e-02, 2.142e-02),
    (5.00000e+00, 3.608e-02, 2.250e-02),
    (6.00000e+00, 3.655e-02, 2.357e-02),
    (8.00000e+00, 3.815e-02, 2.553e-02),
    (1.00000e+01, 4.002e-02, 2.714e-02),
    (1.50000e+01, 4.455e-02, 2.980e-02),
    (2.00000e+01, 4.823e-02, 3.101e-02),
]

M_NA = 22.98977
M_I = 126.90447
M_NAI = M_NA + M_I
W_NA = M_NA / M_NAI
W_I = M_I / M_NAI

EPS = 1e-9


def loglog_interp(e, table, col):
    xs = [row[0] for row in table]
    ys = [row[col] for row in table]
    lo = None
    hi = None
    for k in range(len(xs) - 1):
        x0, x1 = xs[k], xs[k + 1]
        if x1 <= x0:
            continue
        if x0 <= e <= x1:
            lo, hi = k, k + 1
            break
    if lo is None:
        raise ValueError("energy out of range")
    x0, y0 = xs[lo], ys[lo]
    x1, y1 = xs[hi], ys[hi]
    lx = math.log(e)
    lx0, lx1 = math.log(x0), math.log(x1)
    ly0, ly1 = math.log(y0), math.log(y1)
    ly = ly0 + (ly1 - ly0) * (lx - lx0) / (lx1 - lx0)
    return math.exp(ly)


def na_lookup(e, dup_after):
    matches = [row for row in NA if abs(row[0] - e) < EPS]
    if matches:
        if len(matches) == 2:
            row = matches[1] if dup_after else matches[0]
        else:
            row = matches[0]
        return row[1], row[2], True
    return (loglog_interp(e, NA, 1), loglog_interp(e, NA, 2), False)


def i_lookup(e, dup_idx):
    matches = [row for row in I if abs(row[0] - e) < EPS]
    if len(matches) >= 2:
        row = matches[min(dup_idx, len(matches) - 1)]
    else:
        row = matches[0]
    return row[1], row[2]


def build_grid():
    grid = []
    seen = {}
    for (e, _mr, _me) in I:
        c = seen.get(round(e, 12), 0)
        grid.append((e, c))
        seen[round(e, 12)] = c + 1
    i_energies = set(seen.keys())
    for (e, _mr, _me) in NA:
        if round(e, 12) not in i_energies:
            grid.append((e, 0))
    grid.sort(key=lambda t: (t[0], t[1]))
    return grid


def main():
    grid = build_grid()
    rows = []
    interp_energies = []
    for (e, dup_idx) in grid:
        i_mr, i_me = i_lookup(e, dup_idx)
        na_mr, na_me, exact = na_lookup(e, dup_after=(dup_idx >= 1))
        mu_rho = W_NA * na_mr + W_I * i_mr
        mu_en_rho = W_NA * na_me + W_I * i_me
        rows.append((e, mu_rho, mu_en_rho))
        if not exact:
            interp_energies.append(e)

    interp_list = ", ".join("%.5e" % x for x in sorted(set(interp_energies)))

    header_lines = [
        "# NaI (sodium iodide, натрий йодистый) mass attenuation coefficients",
        "# Source: MIXTURE-RULE (Bragg additivity) из ЭЛЕМЕНТНЫХ таблиц NIST X-Ray Mass Attenuation Coefficients,",
        "#   Na (Z=11) https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z11.html",
        "#   I  (Z=53) https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z53.html",
        "#   (составная страница ComTab/sodiumiodide.html отдаёт HTTP 404 -> рассчитано по правилу смеси).",
        "# Extracted/computed 2026-07-06. w_Na=%.6f, w_I=%.6f (M_Na=%.5f, M_I=%.5f, M_NaI=%.5f)." % (W_NA, W_I, M_NA, M_I, M_NAI),
        "# (mu/rho)_NaI = w_Na*(mu/rho)_Na + w_I*(mu/rho)_I; аналогично mu_en/rho. Единицы cm2/g, энергия MeV.",
        "# Сетка = узлы NIST (объединение Na+I). На совпадающих узлах — дословно из NIST.",
        "# В узлах L/K-краёв йода без узла Na (%s MeV) значение Na — лог-лог интерполяция по соседним узлам Na (Na там гладок); значение I — дословно NIST." % interp_list,
        "# K-край йода 0.0331694 MeV сохранён двумя строками (до/после края), не схлопнут.",
        "energy_MeV,mu_rho_cm2_g,mu_en_rho_cm2_g",
    ]

    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "nai_mu_rho_nist.csv")
    lines = list(header_lines)
    for (e, mr, me) in rows:
        lines.append("%.6E,%.6E,%.6E" % (e, mr, me))
    text = "\n".join(lines) + "\n"

    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)

    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    print("OK: %s" % out_path)
    print("rows_data=%d" % len(rows))
    print("energy_range_MeV=%.5E .. %.5E" % (rows[0][0], rows[-1][0]))
    print("interp_Na_nodes=%s" % interp_list)
    print("--- first 3 data rows ---")
    for r in rows[:3]:
        print("%.6E,%.6E,%.6E" % r)
    print("--- last 2 data rows ---")
    for r in rows[-2:]:
        print("%.6E,%.6E,%.6E" % r)


if __name__ == "__main__":
    main()