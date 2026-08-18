# -*- coding: utf-8 -*-
"""Why did Th-232 report show source=default_NaI_63x63 instead of
bootstrap_from_significant_peaks? Inspect lsrm_peaks_table on the
Th-232 .spe and trace which build_fwhm_model branch fires."""
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
    build_fwhm_model, _DEFAULT_NAI_FWHM_MODEL,
    _bootstrap_fwhm_from_significant_peaks,
    _model_is_pathological,
    FwhmModel,
)

SPE = r"C:\LSRM\Work\BG\Gamma-1S\Spe - поверки\Поверка 2024\Маринелли\Th232_420-7-17_Маринелли_0cm.spe"
if not Path(SPE).is_file():
    # fallback to project archive
    SPE = str(REPO / "detectors" / "Gamma-1S" / "reference_spectra" / "archive" / "Th232_420-7-17_Маринелли_0cm.spe")
print(f"=== Th-232 .spe = {SPE}")
spec = read_spectrum(SPE)
print(f"channels={len(spec.counts)}")
pks = spec.extras.get("lsrm_peaks_table") if spec.extras else None
print(f"lsrm_peaks_table present: {bool(pks)}; len={len(pks) if pks else 0}")
if pks:
    print(f"first 3 rows: {pks[:3]}")
    print(f"...rows with E>0 and FWHM>0:")
    valid = [(p.get('energy_keV'), p.get('fwhm_keV')) for p in pks
             if p.get('energy_keV', 0) > 0 and p.get('fwhm_keV', 0) > 0]
    print(f"valid count = {len(valid)}; sample = {valid[:5]}")

print("\n=== Try bootstrap directly (would it find anchors?) ===")
bs = _bootstrap_fwhm_from_significant_peaks(spec)
if bs is None:
    print("bootstrap returned None")
else:
    bm, lbl, n = bs
    print(f"bootstrap OK: source={lbl} n_anchors={n}")
    print(f"  model=(a={bm[0]:.6g}, b={bm[1]:.6g}, c={bm[2]:.6g})")
    def fwhm(model, E):
        a, b, c = model
        v = a + b * E + c * E * E
        return math.sqrt(max(v, 0.0))
    for E in (208.13, 237.51, 327.63, 459.98, 580.08, 900.23, 953.43, 2612.86):
        print(f"  FWHM({E}) = {fwhm(bm, E):.3f}")

print("\n=== Full build_fwhm_model() pipeline ===")
model, src = build_fwhm_model(spec)
print(f"final: source={src}")
assert isinstance(model, FwhmModel), f"build_fwhm_model должен вернуть FwhmModel, got {type(model)}"
print(f"  kind={model.kind}; coefficients(len={len(model.coefficients)})={tuple(f'{c:.6g}' for c in model.coefficients)}")
if model.kind == "quad_fwhm2_in_E":
    a, b, c = model.coefficients
    print(f"  legacy quad form: a={a:.6g}, b={b:.6g}, c={c:.6g}")
    print(f"  pathological={_model_is_pathological(model.coefficients)}")
elif model.kind == "lsrm_poly_sqrt_E":
    print("  F-452 LSRM poly-4 in sqrt(E); _model_is_pathological skipped (quad-only guard).")
print("  FWHM(E) sample:")
for E in (208.13, 237.51, 327.63, 459.98, 580.08, 900.23, 953.43, 2612.86):
    print(f"    FWHM({E}) = {model(E):.3f}")
