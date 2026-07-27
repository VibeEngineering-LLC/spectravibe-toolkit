"""Диагностика формы подложки (continuum) под ключевыми пиками Th-232.

Для каждого пика/мультиплета строит:
  - серая ступенька = sample (counts из .spe)
  - фиолетовая пунктирная = НАША текущая подложка (из json fit_overlay)
  - оранжевая = идеальная erf-step модель (Compton step от главного пика)
  - голубая = total fit envelope

Сохраняет один PNG. Чисто диагностика, ничего в pipeline не правит.
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import erfc

sys.stdout.reconfigure(encoding="utf-8")

PROJ = Path(__file__).resolve().parents[2]
REPORT_DIR = PROJ / "demo_reports" / "2026-06-20_22-20_th232_marinelli_full_with_comparison" / "sample_v2"
JSON_PATH = REPORT_DIR / "Th232_Маринелли_0cm_report.json"
OUT_PATH = REPORT_DIR / "DIAG_continuum_shape.png"
SPE_PATH = PROJ / "detectors" / "Gamma-1S" / "reference_spectra" / "reference_kits" / "Marinelli_1L" / "Th-232" / "Th232_420-7-17_Маринелли_0cm.spe"

sys.path.insert(0, str(PROJ / "scripts"))
from gamma.io.readers import read_spectrum

spec = read_spectrum(str(SPE_PATH))
counts = np.asarray(spec.counts, dtype=float)
live = float(getattr(spec, "live_time", None) or 1.0)
cps = counts / live

cal = getattr(spec, "energy_cal_keV", None)
if cal is None or len(cal) < 2:
    a = b = 0
    c2 = 0
    E_ch = np.arange(len(counts), dtype=float) * 3.0
else:
    a = cal[0]; b = cal[1]; c2 = cal[2] if len(cal) > 2 else 0.0
    ch = np.arange(len(counts), dtype=float)
    E_ch = a + b * ch + c2 * ch * ch

data = json.load(open(JSON_PATH, encoding="utf-8"))
fo = data["fit_overlay"]


def erf_step_model(E, E0, sigma, B_low, B_high):
    """erf-step: B(E) = B_low + (B_high - B_low) * 0.5 * erfc((E - E0)/(sqrt(2)*sigma))."""
    return B_low + (B_high - B_low) * 0.5 * erfc((E - E0) / (np.sqrt(2.0) * sigma))


def fit_erf_step_simple(E, baseline_their):
    """Грубая erf-step через 4 параметра: E0 = середина диапазона,
    sigma = 1/8 диапазона, B_low = последняя точка, B_high = первая точка.

    Это не fit, а наглядная демонстрация: 'вот как ДОЛЖНА выглядеть монотонная
    подложка с erf-step через те же концевые значения'.
    """
    E0 = 0.5 * (E[0] + E[-1])
    sigma = (E[-1] - E[0]) / 8.0
    B_low = float(baseline_their[-1])
    B_high = float(baseline_their[0])
    return erf_step_model(np.asarray(E), E0, sigma, B_low, B_high)


fig = plt.figure(figsize=(15, 14))
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1])
axes = [
    fig.add_subplot(gs[0, 0]),
    fig.add_subplot(gs[0, 1]),
    fig.add_subplot(gs[1, 0]),
    fig.add_subplot(gs[1, 1]),
    fig.add_subplot(gs[2, :]),
]
fig.suptitle(
    "Диагностика формы подложки под пиками Th-232 (Marinelli, Gamma-1С)\n"
    "фиолет: наша подложка · оранж: идеальная erf-step · голуб: fit envelope · сер: данные",
    fontsize=12,
)

cases = []

p2614 = next(
    p for p in fo["peaks"]
    if p.get("nuclide") == "Tl-208" and abs(p.get("energy_keV", 0) - 2614) < 5
)
cg = p2614["continuum_grid"]
cases.append({
    "title": "Tl-208 2614 (singlet) — χ²/dof по 2614 fit'у",
    "E": np.asarray(cg["energies"]),
    "their_baseline": np.asarray(cg["values"]),
    "total": None,
    "comps": [],
})

mc_by_id = {m.get("cluster_id"): m for m in fo["multiplet_continua"]}

for ck, title in [
    ("M1", "Cluster 1 (751-1111) — Ac-228 911/965/969 + Tl-208 861"),
    ("M2", "Cluster 2 (1431-1787) — Ac-228 1588/Bi-212 1620/Ac-228 1631"),
]:
    m = mc_by_id.get(ck)
    if m is None:
        continue
    cases.append({
        "title": title,
        "E": np.asarray(m["E_keV"]),
        "their_baseline": np.asarray(m["continuum"]),
        "total": np.asarray(m["total"]),
        "comps": m.get("components", []),
    })

multiplet_338 = next(
    (m for m in fo["multiplet_continua"]
     if min(m["E_keV"]) > 100 and max(m["E_keV"]) < 400),
    None,
)
if multiplet_338:
    cases.append({
        "title": f"Cluster 122-390 — Pb-212 239/Ra-224 241/... (χ²=986!)",
        "E": np.asarray(multiplet_338["E_keV"]),
        "their_baseline": np.asarray(multiplet_338["continuum"]),
        "total": np.asarray(multiplet_338["total"]),
        "comps": multiplet_338.get("components", []),
    })

for ax, case in zip(axes, cases):
    E = case["E"]
    e0, e1 = float(E[0]), float(E[-1])
    msk = (E_ch >= e0) & (E_ch <= e1)
    E_data = E_ch[msk]
    cnt_data = counts[msk]

    if case["total"] is None:
        bin_width_keV = float(np.mean(np.diff(E)))
        ch_width_keV = float(np.mean(np.diff(E_data))) if len(E_data) >= 2 else 1.0
        bins_per_chunk = max(1, int(round(bin_width_keV / ch_width_keV)))
        cnt_in_chunks = []
        for ei in E:
            j0 = int(np.searchsorted(E_data, ei - bin_width_keV/2))
            j1 = int(np.searchsorted(E_data, ei + bin_width_keV/2))
            if j1 > j0:
                cnt_in_chunks.append(float(cnt_data[j0:j1].sum()))
            else:
                cnt_in_chunks.append(0.0)
        cnt_in_chunks = np.asarray(cnt_in_chunks)
        ax.step(E, cnt_in_chunks, color="gray", lw=0.9, where="mid",
                label=f"образец (counts на {bin_width_keV:.0f} кэВ)")
        their = case["their_baseline"]
        ideal = fit_erf_step_simple(E, case["their_baseline"])
    else:
        ax.plot(E_data, cnt_data, color="gray", lw=0.8, label="образец (counts/канал)")
        their = case["their_baseline"]
        ideal = fit_erf_step_simple(E, case["their_baseline"])
        ax.plot(E, case["total"], color="#1976d2", lw=1.6, label="total fit envelope")

    ax.plot(E, their, color="#7b1fa2", lw=2.0, ls="--",
            label="НАША подложка", marker="o", markersize=4)
    ax.plot(E, ideal, color="#e65100", lw=1.8,
            label="идеал erf-step (через крайние точки)")

    is_mono_their = bool(np.all(np.diff(their) <= 1e-6))
    n_extrema = int(np.sum(np.diff(np.sign(np.diff(their))) != 0))
    monobadge = "монотонна ✓" if is_mono_their else f"НЕ монотонна (экстремумов: {n_extrema}) ⚠"
    ax.set_title(f"{case['title']}\n→ {monobadge}", fontsize=10)
    ax.set_xlabel("Энергия, кэВ")
    ax.set_ylabel("Счёт")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3)

ax5 = axes[4]
case0 = cases[0]
E0 = case0["E"]
their0 = case0["their_baseline"]
ideal0 = fit_erf_step_simple(E0, their0)
ax5.plot(E0, their0, color="#7b1fa2", lw=2.5, ls="--",
         label="НАША подложка (counts/канал)", marker="o", markersize=7)
ax5.plot(E0, ideal0, color="#e65100", lw=2.0,
         label="идеал erf-step (Compton-шаг от 2614)")
for ei, vi in zip(E0, their0):
    ax5.annotate(f"{vi:.0f}", (ei, vi), xytext=(0, 8),
                 textcoords="offset points", ha="center", fontsize=8, color="#7b1fa2")
ax5.axvline(2614.5, color="#1976d2", lw=1, ls=":", label="центр Tl-208 2614")
ax5.set_title(
    "ZOOM ① Tl-208 2614 — ТОЛЬКО подложка (без пика)\n"
    "→ НАША подложка: 543→497→326→359→489→316 (2 экстремума, синусоида) | "
    "идеал: монотонный спад от 543 до 79",
    fontsize=10,
)
ax5.set_xlabel("Энергия, кэВ")
ax5.set_ylabel("Счёт/канал в подложке")
ax5.legend(fontsize=9, loc="upper right")
ax5.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUT_PATH, dpi=80, bbox_inches="tight")
print(f"PNG сохранён: {OUT_PATH}")
print(f"absolute: {OUT_PATH.resolve()}")

OUT_ZOOM = REPORT_DIR / "DIAG_continuum_2614_zoom.png"
fig2, ax = plt.subplots(figsize=(12, 7))
ax.plot(E0, their0, color="#7b1fa2", lw=2.8, ls="--",
        label="НАША подложка (counts/канал)", marker="o", markersize=9)
ax.plot(E0, ideal0, color="#e65100", lw=2.5,
        label="идеал erf-step (Compton-шаг)")
for ei, vi in zip(E0, their0):
    ax.annotate(f"{vi:.0f}", (ei, vi), xytext=(0, 10),
                textcoords="offset points", ha="center", fontsize=10, color="#4a148c")
ax.axvline(2614.5, color="#1976d2", lw=1.2, ls=":", label="центр Tl-208 2614")
ax.set_title(
    "Tl-208 2614 — форма подложки (zoom без пика)\n"
    "НАША: 543→497→326→359→489→316 — НЕ монотонна (2 экстремума, синусоида)  ⚠\n"
    "Идеал erf-step: монотонный спад через те же концевые точки (543→79)",
    fontsize=11,
)
ax.set_xlabel("Энергия, кэВ", fontsize=11)
ax.set_ylabel("Счёт/канал в подложке", fontsize=11)
ax.legend(fontsize=10, loc="upper right")
ax.grid(True, alpha=0.4)
plt.tight_layout()
plt.savefig(OUT_ZOOM, dpi=85, bbox_inches="tight")
print(f"ZOOM PNG: {OUT_ZOOM}")
print(f"absolute: {OUT_ZOOM.resolve()}")