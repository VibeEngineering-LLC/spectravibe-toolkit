# -*- coding: utf-8 -*-
"""F-449 / agent-a-math-2 (2026-06-20) — verify which FWHM branch
build_fwhm_model picks for Ra-226 Marinelli, and what the resulting curve
looks like vs the default NaI 63x63.

Background: after Правка B (bootstrap-from-significant-peaks), the M1
Pb-214 multiplet self-calibration test went from χ²/ν=61.69 baseline to
χ²/ν=242.52. We need to know:
  (1) Which build_fwhm_model branch fired for Ra-226 — bootstrap or
      default?
  (2) If bootstrap, what coefficients and how many anchors?
  (3) FWHM(295.2 keV) and FWHM(351.9 keV) — these are the Pb-214 doublet
      energies; if bootstrap produced a much wider σ at these energies
      the deconvolution would explode in χ².

Run: PYTHONIOENCODING=utf-8 py -3.11 scripts/diag/diag_bootstrap_ra226_fwhm.py
"""
from __future__ import annotations
import sys, math
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from gamma.io.readers import read_spectrum
from gamma.identification.staged_pipeline import (
    build_fwhm_model, _DEFAULT_NAI_FWHM_MODEL, _DEFAULT_NAI_FWHM_MODEL_OBJ,
    _bootstrap_fwhm_from_significant_peaks,
    FwhmModel,
)


def fwhm_at(model, E):
    """F-452 polymorphic: FwhmModel callable | legacy 3-tuple."""
    if isinstance(model, FwhmModel):
        return model(E)
    a, b, c = model
    v = a + b * E + c * E * E
    return math.sqrt(max(v, 0.0))


SPE = (
    REPO / "detectors" / "Gamma-1S" / "reference_spectra"
    / "archive" / "Ra226_420-7-18_Маринелли_0cm.spe"
)
print(f"=== Ra-226 .spe = {SPE}")
spec = read_spectrum(str(SPE))
print(f"channels={len(spec.counts)} energy_cal={spec.energy_cal}")
print(f"has lsrm_peaks_table = {bool(spec.extras and spec.extras.get('lsrm_peaks_table'))}")

print("\n=== Try bootstrap directly ===")
bs = _bootstrap_fwhm_from_significant_peaks(spec)
if bs is None:
    print("bootstrap returned None (insufficient anchors)")
else:
    bm, lbl, n = bs
    print(f"bootstrap OK: source={lbl} n_anchors={n}")
    print(f"  model=(a={bm[0]:.6g}, b={bm[1]:.6g}, c={bm[2]:.6g})")

print("\n=== Full build_fwhm_model() pipeline ===")
model, src = build_fwhm_model(spec)
print(f"final: source={src}")
assert isinstance(model, FwhmModel), f"build_fwhm_model должен вернуть FwhmModel, got {type(model)}"
print(f"  kind={model.kind}; coefficients(len={len(model.coefficients)})={tuple(f'{c:.6g}' for c in model.coefficients)}")
if model.kind == "quad_fwhm2_in_E":
    a, b, c = model.coefficients
    print(f"  legacy quad form: a={a:.6g}, b={b:.6g}, c={c:.6g}")
elif model.kind == "lsrm_poly_sqrt_E":
    print("  F-452 LSRM poly-4 in sqrt(E).")

print("\n=== FWHM curve comparison at Pb-214 doublet ===")
print(f"{'E_keV':>8s}  {'default':>9s}  {'chosen':>9s}  {'Δ(keV)':>8s}  {'Δ(%)':>7s}")
for E in (60.0, 100.0, 200.0, 295.21, 351.92, 661.7, 1460.8, 2614.5):
    f_def = fwhm_at(_DEFAULT_NAI_FWHM_MODEL, E)
    f_new = fwhm_at(model, E)
    d = f_new - f_def
    pct = 100 * d / f_def if f_def > 0 else float("nan")
    print(f"{E:8.2f}  {f_def:9.3f}  {f_new:9.3f}  {d:+8.3f}  {pct:+6.1f}%")

print("\n=== Σ at Pb-214 (M1 deconvolution input) ===")
for E in (295.21, 351.92):
    s_def = fwhm_at(_DEFAULT_NAI_FWHM_MODEL, E) / 2.354820045
    s_new = fwhm_at(model, E) / 2.354820045
    print(f"  E={E:.2f}  σ_default={s_def:.3f} keV  σ_chosen={s_new:.3f} keV  ratio={s_new/s_def:.3f}")
