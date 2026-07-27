# -*- coding: utf-8 -*-
"""F-449 / agent-a-math-3 (2026-06-20) — debug which anchors bootstrap
collects on Th-232 vs Ra-226 vs ANY .spe. Mirrors
_bootstrap_fwhm_from_significant_peaks internals so the rejection reason
for each candidate is visible.

Usage:
    PYTHONIOENCODING=utf-8 py -3.11 scripts/diag/diag_bootstrap_anchors_v2.py [<spe>]

Default targets are both Th-232 and Ra-226 Marinelli archives.
"""
from __future__ import annotations
import sys, math, os
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

import numpy as np
from gamma.io.readers import read_spectrum
from gamma.identification.staged_pipeline import (
    _DEFAULT_NAI_FWHM_MODEL, _eval_fwhm2_quadratic,
)
from gamma.peaks.search import mariscotti_search
from gamma.calibration.fwhm_provider import _measure_fwhm_channels


def trace(spe_path: str, *, sn_min: float = 5.0,
          energy_min_keV: float = 200.0, energy_max_keV: float = 2700.0,
          isolation_factor: float = 2.0, height_min: float = 30.0,
          min_anchors: int = 4):
    print(f"\n========== {spe_path} ==========")
    spec = read_spectrum(spe_path)
    counts = np.asarray(spec.counts, dtype=np.float64)
    n_ch = len(counts)
    print(f"channels={n_ch}, e_cal={spec.energy_cal}")
    if not spec.energy_cal or len(spec.energy_cal) < 2:
        print("  no e-cal — abort")
        return
    gain = abs(float(spec.energy_cal[1]))
    seed_fwhm_keV = math.sqrt(_eval_fwhm2_quadratic(_DEFAULT_NAI_FWHM_MODEL, 661.0))
    seed_fwhm_ch = max(3.0, seed_fwhm_keV / gain)
    iso_dist_ch = isolation_factor * seed_fwhm_ch
    print(f"  seed_fwhm@661={seed_fwhm_keV:.2f} keV → seed_fwhm_ch={seed_fwhm_ch:.2f}")
    print(f"  isolation_dist_ch={iso_dist_ch:.2f} ({isolation_factor}·seed)")
    print(f"  filter: sn≥{sn_min}, height≥{height_min}, E∈[{energy_min_keV},{energy_max_keV}]")

    found = mariscotti_search(
        counts=counts,
        fwhm_channels=float(seed_fwhm_ch),
        sigma_threshold=float(sn_min),
        min_separation_factor=1.0,
        edge_margin=10,
    )
    if not found:
        print("  Mariscotti: 0 candidates")
        return
    chs_sorted = sorted(p.channel for p in found)
    print(f"\n  Mariscotti candidates: {len(found)}")
    print(f"  {'ch':>5}  {'E_keV':>8}  {'sig':>6}  {'h':>8}  {'reject':>40}")
    print(f"  {'-'*5}  {'-'*8}  {'-'*6}  {'-'*8}  {'-'*40}")
    anchors = []
    for p in sorted(found, key=lambda q: q.channel):
        reasons = []
        E_keV = spec.channel_to_energy(int(p.channel))
        sig = getattr(p, "significance", 0.0)
        h = getattr(p, "height", 0.0)
        if sig < sn_min:
            reasons.append(f"sn<{sn_min}")
        if h < height_min:
            reasons.append(f"h<{height_min:.0f}")
        if E_keV is None or not (energy_min_keV <= E_keV <= energy_max_keV):
            reasons.append("E_out_of_band")
        # isolation
        nearest = min((abs(q - p.channel) for q in chs_sorted if q != p.channel),
                      default=float("inf"))
        if nearest < iso_dist_ch:
            reasons.append(f"iso(nearest={nearest:.1f}<{iso_dist_ch:.1f})")
        # measure FWHM
        fwhm_ch = _measure_fwhm_channels(counts, int(p.channel), seed_fwhm_ch)
        if fwhm_ch is None or fwhm_ch <= 1.0:
            reasons.append("fwhm_ch_bad")
        ch_f = float(p.channel)
        dE_dN = sum(
            i * float(a) * (ch_f ** (i - 1))
            for i, a in enumerate(spec.energy_cal) if i > 0
        )
        if dE_dN <= 0:
            reasons.append("dE/dch≤0")
        E_str = f"{E_keV:.1f}" if E_keV is not None else "?"
        if reasons:
            print(f"  {int(p.channel):5d}  {E_str:>8}  {sig:6.2f}  {h:8.1f}  {','.join(reasons)}")
        else:
            fwhm_keV = fwhm_ch * dE_dN
            anchors.append((E_keV, fwhm_keV, int(p.channel), sig, h))
            print(f"  {int(p.channel):5d}  {E_str:>8}  {sig:6.2f}  {h:8.1f}  KEEP fwhm={fwhm_keV:.2f}")
    print(f"\n  KEPT anchors: {len(anchors)} (need ≥{min_anchors})")
    if len(anchors) >= 3:
        Es = np.array([a[0] for a in anchors])
        Fs = np.array([a[1] for a in anchors])
        print(f"\n  Range coverage: E_min={Es.min():.1f}  E_max={Es.max():.1f}")
        bands = [(200, 800), (800, 1500), (1500, 2700)]
        for lo, hi in bands:
            n_b = int(((Es >= lo) & (Es < hi)).sum())
            print(f"    [{lo},{hi}): {n_b}")
        # try sqrt-form fit FWHM^2 = a + b*E
        A2 = np.vstack([np.ones_like(Es), Es]).T
        sol, *_ = np.linalg.lstsq(A2, Fs * Fs, rcond=None)
        a, b = float(sol[0]), float(sol[1])
        print(f"\n  sqrt-fit: a={a:+.4g}  b={b:+.4g}")
        print(f"  sqrt-FWHM samples:")
        for E in (60, 100, 200, 295, 352, 500, 661, 1461, 2614):
            v = a + b * E
            fw = math.sqrt(max(v, 0.0))
            print(f"    E={E:>4}  FWHM={fw:6.2f}  (default={math.sqrt(_eval_fwhm2_quadratic(_DEFAULT_NAI_FWHM_MODEL, E)):6.2f})")


def main():
    if len(sys.argv) >= 2:
        trace(sys.argv[1])
        return
    targets = []
    th232 = REPO / "detectors" / "Gamma-1S" / "reference_spectra" / "archive" / "Th232_420-7-17_Маринелли_0cm.spe"
    ra226 = REPO / "detectors" / "Gamma-1S" / "reference_spectra" / "archive" / "Ra226_420-7-18_Маринелли_0cm.spe"
    if th232.is_file():
        targets.append(str(th232))
    if ra226.is_file():
        targets.append(str(ra226))
    # external archive on operator's box
    extra_th = r"C:\LSRM\Work\BG\Gamma-1S\Spe - поверки\Поверка 2024\Маринелли\Th232_420-7-17_Маринелли_0cm.spe"
    extra_ra = r"C:\LSRM\Work\BG\Gamma-1S\Spe - поверки\Поверка 2024\Маринелли\Ra226_420-7-18_Маринелли_0cm.spe"
    for p in (extra_th, extra_ra):
        if os.path.isfile(p):
            targets.append(p)
    if not targets:
        print("no targets found")
        return
    for t in targets:
        trace(t)


if __name__ == "__main__":
    main()
