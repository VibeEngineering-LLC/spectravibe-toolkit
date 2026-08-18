# -*- coding: utf-8 -*-
"""Драйвер: сравнение трёх методов мощности дозы из гамма-спектра (форма кривой).

  (3)   ЛСРМ:      w = f10(E)·(μen/ρ)_air·E     — без ε детектора
  (2)   ICRP-74:   w = h_Φ(E)                    — counts≈флюенс, без ε
  (2+ε) ICRP-74:   w = h_Φ(E)/ε_tot(E)           — counts→флюенс через ε NaI

Все нормированы Σ=1 (доля вклада в дозу), привязаны к внешнему якорю MD_TOTAL [мкЗв/ч].
ε входит ТОЛЬКО в метод (2+ε) (делитель); в (3) и (2) — НЕ участвует.

CLI:
  python compare_methods.py [--spectrum PATH] [--md-anchor 0.09] [--out-dir DIR]
  Дефолт --spectrum = «Фон подвал радон.xml» (первый прогон, воспроизводимость).
  --md-anchor — якорь МД прибора [мкЗв/ч], СПЕЦИФИЧЕН для спектра (0.09 = подвал радон).
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).parent))

import icrp74_hphi as icrp
from nai_efficiency import eps_tot, RHO_NAI, D_CM
from lsrm_coeffs import (
    parse_spectrum, muen_air_interp, f10_interp, E_LO, E_HI, DEFAULT_SPECTRUM,
)

_ANCHOR_LINES = [(352, "Pb-214"), (609, "Bi-214"), (1120, "Bi-214"),
                 (1461, "K-40"), (1764, "Bi-214"), (2614, "Tl-208")]
_BANDS = [(50, 200), (200, 500), (500, 800), (800, 1200),
          (1200, 1600), (1600, 2000), (2000, 3000)]


def compute(spectrum_path, md_anchor):
    """Возвращает dict с E, N, тремя нормированными кривыми c3/c2/c2e и весами."""
    s = parse_spectrum(spectrum_path)
    E = s["energy_keV"]
    N = s["counts"] / s["live"]
    inr = (E >= E_LO) & (E <= E_HI)

    muen = muen_air_interp(E)
    f10 = f10_interp(E)
    hphi = icrp.h_phi_interp(E)     # pSv·cm²
    eps = eps_tot(E)

    w3 = np.where(inr, f10 * muen * E, 0.0)      # (3) ЛСРМ — без ε
    w2 = np.where(inr, hphi, 0.0)                 # (2) ICRP-74 — без ε
    w2e = np.where(inr, hphi / eps, 0.0)         # (2+ε) — ε здесь (делитель)

    def curve(w):
        c = w * N
        tot = c.sum()
        return c / tot if tot > 0 else c

    return {"s": s, "E": E, "N": N, "inr": inr, "muen": muen, "f10": f10,
            "hphi": hphi, "eps": eps, "md": md_anchor,
            "c3": curve(w3), "c2": curve(w2), "c2e": curve(w2e)}


def report(r):
    s, E, N = r["s"], r["E"], r["N"]
    c3, c2, c2e = r["c3"], r["c2"], r["c2e"]
    print("=" * 74)
    print("СРАВНЕНИЕ ТРЁХ МЕТОДОВ МОЩНОСТИ ДОЗЫ — форма кривой (Σ=1)")
    print("Спектр: " + s["path"].name)
    print("=" * 74)
    print(f"каналов={len(N)}, live={s['live']:.0f} с, rate={s['rate_cps']:.2f} имп/с")
    print(f"диапазон {E_LO:.0f}-{E_HI:.0f} кэВ, каналов в диапазоне={int(r['inr'].sum())}")
    print(f"якорь MD_TOTAL={r['md']:.3f} мкЗв/ч\n")
    print("  (3)   ЛСРМ:      w = f10·(μen/ρ)air·E       [без ε детектора]")
    print("  (2)   ICRP-74:   w = h_Φ(E)                  [counts≈флюенс, без ε]")
    print("  (2+ε) ICRP-74:   w = h_Φ(E)/ε_tot(E)         [counts→флюенс через ε NaI]\n")

    print(f"  {'полоса,кэВ':>12} | {'(3)ЛСРМ':>9} | {'(2)ICRP':>9} | {'(2+ε)ICRP':>10}")
    print("  " + "-" * 52)
    for lo, hi in _BANDS:
        m = (E >= lo) & (E < hi)
        print(f"  {lo:5.0f}-{hi:<6.0f}| {c3[m].sum()*100:8.2f}% | "
              f"{c2[m].sum()*100:8.2f}% | {c2e[m].sum()*100:9.2f}%")
    print()
    for name, c in [("(3)ЛСРМ", c3), ("(2)ICRP", c2), ("(2+ε)ICRP", c2e)]:
        cum = np.cumsum(c)
        e50 = E[min(int(np.searchsorted(cum, 0.5)), len(E) - 1)]
        e90 = E[min(int(np.searchsorted(cum, 0.9)), len(E) - 1)]
        print(f"  {name:>10}: 50% дозы ниже {e50:6.0f} кэВ, 90% ниже {e90:6.0f} кэВ")


def plot(r, out_dir):
    s, E, N = r["s"], r["E"], r["N"]
    c3, c2, c2e = r["c3"], r["c2"], r["c2e"]
    fig, ax = plt.subplots(3, 1, figsize=(12, 12), sharex=True)
    ax[0].step(E, N, where="mid", lw=0.7, color="0.35")
    ax[0].set_yscale("log"); ax[0].set_ylabel("имп/с·канал")
    ax[0].set_title("Спектр: " + s["path"].name)
    ax[0].grid(True, alpha=0.3); ax[0].set_xlim(0, 3000)

    ax[1].step(E, c3, where="mid", lw=1.0, color="tab:red", label="(3) ЛСРМ  f10·μen·E")
    ax[1].step(E, c2, where="mid", lw=1.0, color="tab:blue", label="(2) ICRP-74  h_Φ")
    ax[1].step(E, c2e, where="mid", lw=1.2, color="tab:green", label="(2+ε) ICRP-74  h_Φ/ε")
    ax[1].set_ylabel("Доля дозы (Σ=1)")
    ax[1].set_title("Кривая мощности дозы — 3 метода (форма, привязка "
                    + str(r["md"]) + " мкЗв/ч)")
    ax[1].legend(); ax[1].grid(True, alpha=0.3); ax[1].set_xlim(0, 3000)
    for e, lbl in _ANCHOR_LINES:
        ax[1].axvline(e, color="0.6", ls=":", lw=0.5)
        ax[1].text(e, ax[1].get_ylim()[1] * 0.95, lbl, rotation=90,
                   fontsize=7, va="top", ha="right", color="0.4")

    Egrid = np.linspace(50, 3000, 600)
    ax[2].plot(Egrid, eps_tot(Egrid), color="tab:green", lw=1.6)
    ax[2].set_xlabel("Энергия, кэВ"); ax[2].set_ylabel("ε_tot NaI 2×2\"")
    ax[2].set_title("Кривая эффективности ε_tot(E)=1−exp(−(μ/ρ)·ρ·d), NaI ρ=%.3f d=%.2f см "
                    "— делитель в (2+ε); в (3)/(2) НЕ участвует" % (RHO_NAI, D_CM))
    ax[2].grid(True, alpha=0.3); ax[2].set_xlim(0, 3000); ax[2].set_ylim(0, 1.03)
    for e in (352, 609, 1120, 1461, 1764, 2614):
        ev = float(eps_tot(e))
        ax[2].axvline(e, color="0.6", ls=":", lw=0.5)
        ax[2].plot(e, ev, "o", color="tab:red", ms=4)
        ax[2].annotate(f"{ev:.2f}", (e, ev), textcoords="offset points",
                       xytext=(4, 4), fontsize=7, color="0.3")
    fig.tight_layout()
    png = Path(out_dir) / "compare_dose_methods.png"
    fig.savefig(png, dpi=130)
    return png


def save_csv(r, out_dir):
    s, E, N = r["s"], r["E"], r["N"]
    c3, c2, c2e, md = r["c3"], r["c2"], r["c2e"], r["md"]
    csv = Path(out_dir) / "compare_dose_methods.csv"
    np.savetxt(csv, np.column_stack([s["ch"], E, N, r["muen"], r["f10"], r["hphi"],
                                     r["eps"], c3, c2, c2e, c3 * md, c2 * md, c2e * md]),
               delimiter=",",
               header="channel,energy_keV,rate_cps,muen_air,f10,h_phi_pSvcm2,eps_tot,"
                      "c3_lsrm,c2_icrp,c2e_icrp,d3_uSvh,d2_uSvh,d2e_uSvh",
               comments="", fmt="%.6g")
    return csv


def main(argv=None):
    ap = argparse.ArgumentParser(description="Сравнение 3 методов МД из гамма-спектра")
    ap.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM,
                    help="путь к спектру (XML AtomSpectra/BecqMoni); дефолт — подвал радон")
    ap.add_argument("--md-anchor", type=float, default=0.09,
                    help="якорь МД прибора [мкЗв/ч], специфичен для спектра (дефолт 0.09)")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="каталог вывода PNG/CSV (дефолт — рядом со спектром)")
    a = ap.parse_args(argv)
    out_dir = a.out_dir or a.spectrum.parent
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    r = compute(a.spectrum, a.md_anchor)
    report(r)
    png = plot(r, out_dir)
    csv = save_csv(r, out_dir)
    print(f"\nГрафик: {png}")
    print(f"CSV: {csv}")


if __name__ == "__main__":
    main()