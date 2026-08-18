# -*- coding: utf-8 -*-
"""F-449 forensic: trace the EXACT origin of the 2614 singlet's reported FWHM.

The 20:07 report shows source='default_NaI_63x63' (curve→107.9 keV @2614) yet the
singlet fit_overlay entry carries FWHM=116.63 (sigma=49.528). Static analysis could
not reconcile this. This diag runs the REAL v2 pipeline in-process, captures the
in-memory result, and prints — for every matched_line near 2614 — the gauss_sigma_keV
plus what the calibration curve WOULD give, so we see whether the singlet sigma is
the calibration value or an independent free fit.

Run:  PYTHONIOENCODING=utf-8 py -3.11 scripts/diag/diag_trace_singlet_sigma.py
"""
from __future__ import annotations
import sys, os, math
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

SAMPLE = r"C:\LSRM\Work\BG\Gamma-1S\Spe - поверки\Поверка 2024\Маринелли\Th232_420-7-17_Маринелли_0cm.spe"
BG = r"C:\LSRM\Work\BG\Gamma-1S\Spe - поверки\Поверка 2024\Фон закр кр\Фон закр кр вода_13.spe"

# ── instrument identify.py singlet sigma cache ────────────────────────────────
import gamma.identification.identify as _idmod  # noqa: E402

_orig_identify = _idmod.identify_nuclides
_cache_log = []


def _patched_identify(*args, **kwargs):
    res = _orig_identify(*args, **kwargs)
    # res holds matched_lines; dump every near-2614 line's sigma
    for ni in (getattr(res, "detected_nuclides", None) or ()):
        for m in getattr(ni, "matched_lines", ()) or ():
            e = getattr(m, "peak_E_keV", None)
            if e is not None and 2580 < float(e) < 2650:
                _cache_log.append((
                    ni.nuclide,
                    float(e),
                    getattr(m, "gauss_sigma_keV", None),
                    getattr(m, "peak_sigma", None),
                    getattr(m, "peak_area_source", None),
                    getattr(m, "is_characteristic", None),
                ))
    return res


_idmod.identify_nuclides = _patched_identify
# staged_pipeline imported the symbol by-name; patch there too
import gamma.identification.staged_pipeline as _sp  # noqa: E402
if hasattr(_sp, "identify_nuclides"):
    _sp.identify_nuclides = _patched_identify

# ── run the real v2 pipeline (same kwargs run_skill builds) ───────────────────
from gamma.experimental.v2_integration import analyze_and_report_v2  # noqa: E402

out_dir = REPO / "demo_reports" / "_diag_trace_singlet"
out_dir.mkdir(parents=True, exist_ok=True)

print("Running analyze_and_report_v2 ... (this takes ~1-2 min)")
art = analyze_and_report_v2(
    SAMPLE,
    output_dir=str(out_dir),
    background_path=BG,
    sample_mass_kg=1.6,
    background_auto="apply",
    peak_search_method="mariscotti",
    allow_stage2=True,
    allow_stage3=False,
    write_json=True, write_markdown=False, write_html=False,
    write_plots=False, write_technical_pdf=False,
    complete_workflow=True,
)

result = art.get("result")
fwhm_model = getattr(result, "fwhm_model", None)
fwhm_src = getattr(result, "fwhm_model_source", None)
print("\n=== result.fwhm_model =", fwhm_model, " source =", fwhm_src)


def _curve(E):
    if not fwhm_model:
        return None
    a, b, c = fwhm_model
    v = a + b * E + c * E * E
    return math.sqrt(v) if v > 0 else None


print("\n=== identify.py singlet-cache near 2614 (in-memory matched_lines) ===")
if not _cache_log:
    print("  (no near-2614 matched_lines captured)")
for nuc, e, gsk, psig, src, ischar in _cache_log:
    cv = _curve(e)
    fwhm_from_gsk = gsk * 2.35482 if gsk else None
    print("  %-7s E=%.2f  gauss_sigma_keV=%s  FWHM_from_gsk=%s  curve_FWHM=%s  area_src=%s char=%s"
          % (nuc, e,
             ("%.4f" % gsk) if gsk else "None",
             ("%.3f" % fwhm_from_gsk) if fwhm_from_gsk else "None",
             ("%.3f" % cv) if cv else "None",
             src, ischar))

# ── now read what landed in the JSON fit_overlay ──────────────────────────────
import json, glob  # noqa: E402
cands = glob.glob(str(out_dir / "**" / "*_report.json"), recursive=True)
if cands:
    d = json.load(open(cands[0], encoding="utf-8"))
    print("\n=== JSON fit_overlay peaks near 2614 ===")
    for p in d.get("fit_overlay", {}).get("peaks", []):
        e = p.get("energy_keV") or 0
        if e and 2580 < e < 2650:
            s = p.get("sigma_keV")
            print("  E=%.2f sigma_keV=%s FWHM=%.3f source=%s"
                  % (e, s, (s or 0) * 2.35482, p.get("source")))
    fc = d.get("calibration", {}).get("fwhm_cal", {})
    print("\n=== JSON fwhm_cal source =", fc.get("source"), " coef =", fc.get("coefficients"))
else:
    print("\n(no JSON report produced)")
