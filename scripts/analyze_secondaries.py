from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""
Secondary-peak characterisation on clean Cs-137 and K-40 reference fixtures
(F-37 / Variant "reference samples library", v1.7.15).

For every Cs-137 (E_gamma = 661.66 keV) and K-40 (E_gamma = 1460.82 keV)
spectrum available across all measured geometries on the Gamma-1S
NaI 63x63 detector, identify each detected Mariscotti peak as one of:

  photopeak       - the primary gamma photopeak
  compton_edge    - 2 E^2 / (m_e + 2 E)
  backscatter     - E / (1 + 2 E / m_e); photon scattering off shielding
  single_escape   - E - 511.0  (pair production, one annihilation gamma escapes)
  double_escape   - E - 1022.0 (both annihilation gammas escape)
  xray_escape     - E - 28.0   (iodine K X-ray escape from NaI)
  ic_xray         - Ba K X-rays at 32 keV (Cs-137 only; from library `ic_xrays`)
  k40_natural     - the natural-background K-40 1460.82 keV (in Cs-137 spectra)
  unknown         - anything else (typically faint pile-up / sample contamination)

For each (geometry, source ID) row we report:
  - photopeak area S_p
  - per-secondary: position residual (measured - theoretical), area S_s,
    and the geometry-independent intensity ratio R = S_s / S_p
  - FWHM (measured) of the secondary

The output is a JSON catalogue `data/secondary_peaks.json` (per-nuclide
expected features with statistics across geometries) plus a stdout
diagnostic table.

Methodology references:
  - Knoll, "Radiation Detection and Measurement" 4th Ed., chap. 10 (Compton
    edge geometry), chap. 11.A.5 (backscatter), chap. 12.B.2 (escape peaks).
  - Gilmore & Joss, "Practical Gamma-ray Spectrometry" 3rd Ed., chap. 2
    (interactions) and chap. 6 (NaI(Tl) spectral artefacts).
"""


import json
import math
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).parent))

import numpy as np

from gamma.io.readers import read_spectrum
from gamma.peaks.search import mariscotti_search
from validate_certs import make_lsrm_fwhm_provider


# ---------------------------------------------------------------------------
# Physics
# ---------------------------------------------------------------------------

M_E_KEV = 510.9989  # electron rest mass energy


def compton_edge(E_keV: float) -> float:
    return 2.0 * E_keV * E_keV / (M_E_KEV + 2.0 * E_keV)


def backscatter(E_keV: float) -> float:
    return E_keV / (1.0 + 2.0 * E_keV / M_E_KEV)


@dataclass(frozen=True)
class ExpectedFeature:
    name: str
    E_keV: float
    tolerance_keV: float = 30.0   # NaI is broad; allow ~half FWHM
    note: str = ""


def expected_features_for(nuclide: str, E_gamma: float) -> list:
    """Theoretical secondary-peak positions for a primary gamma line."""
    feats = [
        ExpectedFeature("photopeak",     E_gamma, 30.0, "primary photopeak"),
        ExpectedFeature("compton_edge",  compton_edge(E_gamma), 60.0,
                        "max energy electron from 180-deg Compton scattering"),
        ExpectedFeature("backscatter",   backscatter(E_gamma), 40.0,
                        "photon at 180-deg scattering off shielding"),
    ]
    # Pair production products only relevant if E_gamma >= ~1100 keV
    if E_gamma > 1022.0 + 100.0:
        feats.append(ExpectedFeature("single_escape", E_gamma - 511.0, 40.0,
                                     "one annihilation gamma escapes detector"))
        feats.append(ExpectedFeature("double_escape", E_gamma - 1022.0, 40.0,
                                     "both annihilation gammas escape detector"))
    # Iodine K X-ray escape from NaI (only resolvable for low-E primary
    # gammas where 28 keV is a significant fraction)
    if E_gamma < 200.0:
        feats.append(ExpectedFeature("xray_escape", E_gamma - 28.0, 15.0,
                                     "iodine K X-ray escape from NaI crystal"))
    # Cs-137 only: Ba K X-rays at ~32 keV from internal conversion
    if nuclide == "Cs-137":
        feats.append(ExpectedFeature("ic_xray_Ba_Ka", 32.0, 10.0,
                                     "Ba K-alpha X-rays from 8% IC decay branch"))
    # Background K-40 contaminates every long-integration spectrum
    if nuclide == "Cs-137":
        feats.append(ExpectedFeature("k40_natural", 1460.82, 30.0,
                                     "natural-background K-40 photopeak"))
    # Backscatter from a Compton-scatter off the source itself produces an
    # annihilation peak only if E > 1022 - skip
    return feats


# ---------------------------------------------------------------------------
# Fixture inventory
# ---------------------------------------------------------------------------

ROOT = Path(__file__).parent.parent
# F-83 (v1.12.0): Gamma-1S reference spectra live under detectors/Gamma-1S/.
from gamma.detectors.gamma1s import (
    DEFAULT_REFERENCE_DIR as REF,
    SECONDARY_PEAKS_PATH as _OUT_SECONDARY_PEAKS_PATH,
)


@dataclass
class Fixture:
    nuclide: str       # "Cs-137" or "K-40"
    E_primary: float   # primary photopeak energy in keV
    geometry: str
    path: Path
    label: str         # short tag for the table


def gather_fixtures() -> list:
    out = []
    # ---- Cs-137 ----
    cs_E = 661.66
    cs_files = [
        ("5cm",        REF / "Cs-137__163_2017.spe", "5cm/SRC-02"),
        ("25cm",       REF / "Точечная-25см" / "Cs-137 №SRC-02_Точечная-25см_25cm.spe", "25cm/SRC-02"),
        ("Денnта-120мл", REF / "Дента-120мл" / "Cs137_420-7-14_Дента-120мл_0cm.spe", "Дента/7-14"),
        ("Денnта-120мл", REF / "Дента-120мл" / "Cs137_420-7-15_Дента-120мл_0cm.spe", "Дента/7-15"),
        ("Петри-60мл", REF / "Петри-60мл" / "Cs137_420-7-14_Петри-60мл_0cm.spe", "Петри/7-14"),
        ("Петри-60мл", REF / "Петри-60мл" / "Cs137_420-7-15_Петри-60мл_0cm.spe", "Петри/7-15"),
        ("Маринелли",  REF / "Cs137_420-7-14_Маринелли_0cm.spe", "Мар/7-14"),
        ("Маринелли",  REF / "Cs137_420-7-15_Маринелли_0cm.spe", "Мар/7-15"),
        ("M-source",   REF / "M_cs_легкий_2001-2005.spe", "M_cs/легкий"),
        ("M-source",   REF / "M_cs_тяж_2001-2005.spe", "M_cs/тяж"),
    ]
    for g, p, tag in cs_files:
        if p.is_file():
            out.append(Fixture("Cs-137", cs_E, g, p, tag))
    # ---- K-40 ----
    k_E = 1460.82
    k_files = [
        ("Дента-120мл", REF / "Дента-120мл" / "K40_420-7-20_Дента-120мл_0cm.spe", "Дента/7-20"),
        ("Дента-120мл", REF / "Дента-120мл" / "K40_420-7-21_Дента-120мл_0cm.spe", "Дента/7-21"),
        ("Петри-60мл", REF / "Петри-60мл"  / "K40_420-7-20_Петри-60мл_0cm.spe", "Петри/7-20"),
        ("Петри-60мл", REF / "Петри-60мл"  / "K40_420-7-21_Петри-60мл_0cm.spe", "Петри/7-21"),
        ("Маринелли",  REF / "K40_420-7-20_Маринелли_0cm.spe", "Мар/7-20"),
        ("Маринелли",  REF / "K40_420-7-21_Маринелли_0cm.spe", "Мар/7-21"),
        ("M-source",   REF / "M_k_легкий_2001-2005.spe", "M_k/легкий"),
    ]
    for g, p, tag in k_files:
        if p.is_file():
            out.append(Fixture("K-40", k_E, g, p, tag))
    return out


# ---------------------------------------------------------------------------
# Per-fixture analysis
# ---------------------------------------------------------------------------

@dataclass
class MatchedFeature:
    feature: str
    expected_E: float
    measured_E: float
    residual_keV: float
    area: float
    fwhm_keV: float
    sigma_significance: float
    ratio_to_photopeak: Optional[float] = None


def analyze_one(fx: Fixture) -> list:
    spec = read_spectrum(str(fx.path))
    fwhm_at = make_lsrm_fwhm_provider(spec)
    counts = np.asarray(spec.counts, dtype=np.float64)
    found = mariscotti_search(counts, fwhm_channels=fwhm_at, sigma_threshold=3.0)

    # Compute area for each found peak via simple ROI sum minus linear baseline
    # over +/- 2 sigma (Cowell-style). For secondaries we don't need Lsrm-table
    # accuracy -- areas are diagnostic, not certificate-level quantitative.
    def cowell_area(peak_ch: int, fwhm_ch: float) -> tuple:
        sigma = fwhm_ch / 2.355
        roi = int(round(2.5 * sigma))
        lo = max(0, peak_ch - roi)
        hi = min(counts.size, peak_ch + roi + 1)
        # baseline = mean of 1 sigma outside each side
        bg_l_lo = max(0, peak_ch - int(round(4 * sigma)))
        bg_l_hi = max(bg_l_lo + 1, peak_ch - int(round(2.5 * sigma)))
        bg_r_lo = min(counts.size - 1, peak_ch + int(round(2.5 * sigma)))
        bg_r_hi = min(counts.size, peak_ch + int(round(4 * sigma)))
        if bg_l_hi > bg_l_lo and bg_r_hi > bg_r_lo:
            bg_per_ch = 0.5 * (counts[bg_l_lo:bg_l_hi].mean()
                               + counts[bg_r_lo:bg_r_hi].mean())
        else:
            bg_per_ch = 0.0
        gross = float(counts[lo:hi].sum())
        net = gross - bg_per_ch * (hi - lo)
        return max(0.0, net), abs(net)**0.5 + bg_per_ch ** 0.5

    feats = expected_features_for(fx.nuclide, fx.E_primary)
    matched: list = []
    photopeak_area: Optional[float] = None
    photopeak_E: Optional[float] = None

    # Resolve each expected feature -> nearest found peak within tolerance
    used_peaks = set()
    for exp in feats:
        if exp.E_keV <= 0:
            continue
        best = None
        best_dE = float("inf")
        for p in found:
            if p.channel in used_peaks:
                continue
            E_meas = spec.channel_to_energy(p.channel)
            dE = abs(E_meas - exp.E_keV)
            if dE < exp.tolerance_keV and dE < best_dE:
                best = (p, E_meas, dE)
                best_dE = dE
        if best is None:
            continue
        peak, E_meas, dE = best
        used_peaks.add(peak.channel)
        # FWHM in keV from local FWHM provider
        fw_keV = fwhm_at(peak.channel) * float(spec.energy_cal[1]
                                                if len(spec.energy_cal) > 1 else 1.0)
        area, _ = cowell_area(peak.channel, peak.fwhm_channels)
        mf = MatchedFeature(
            feature=exp.name,
            expected_E=exp.E_keV,
            measured_E=E_meas,
            residual_keV=E_meas - exp.E_keV,
            area=area,
            fwhm_keV=fw_keV,
            sigma_significance=peak.significance,
        )
        matched.append(mf)
        if exp.name == "photopeak":
            photopeak_area = area
            photopeak_E = E_meas

    # Fill in ratios once photopeak is known
    if photopeak_area and photopeak_area > 0:
        for m in matched:
            if m.feature != "photopeak":
                m.ratio_to_photopeak = m.area / photopeak_area
    return matched, photopeak_E, photopeak_area


# ---------------------------------------------------------------------------
# Pretty printing
# ---------------------------------------------------------------------------

def fmt_pct(v: Optional[float]) -> str:
    if v is None:
        return "    -    "
    return f"{100.0*v:>7.3f}%"


def main() -> int:
    fixtures = gather_fixtures()
    print(f"Loaded {len(fixtures)} fixtures "
          f"({sum(1 for f in fixtures if f.nuclide=='Cs-137')} Cs-137 + "
          f"{sum(1 for f in fixtures if f.nuclide=='K-40')} K-40)\n")

    # Per-nuclide aggregator -- collect ratios across geometries
    aggregator: dict = {}     # {(nuclide, feature) -> [(geometry, ratio, residual_keV), ...]}
    fixture_rows: list = []   # for stdout table

    for fx in fixtures:
        try:
            matched, pp_E, pp_area = analyze_one(fx)
        except Exception as e:
            print(f"  ERR {fx.label}: {e}")
            continue
        if pp_area is None or pp_area <= 0:
            print(f"  ! {fx.nuclide} {fx.geometry:>14} {fx.label:>14}: photopeak not found")
            continue
        fixture_rows.append((fx, pp_E, pp_area, matched))
        for m in matched:
            if m.feature == "photopeak":
                continue
            key = (fx.nuclide, m.feature)
            aggregator.setdefault(key, []).append({
                "geometry": fx.geometry,
                "label": fx.label,
                "measured_E_keV": m.measured_E,
                "expected_E_keV": m.expected_E,
                "residual_keV": m.residual_keV,
                "ratio_to_photopeak": m.ratio_to_photopeak,
                "fwhm_keV": m.fwhm_keV,
            })

    # === stdout: per-fixture annotated peaks ===
    print("=" * 110)
    print("Per-fixture matched secondary features (E_meas - E_theory, S_secondary / S_photopeak):")
    print("=" * 110)
    for fx, pp_E, pp_area, matched in fixture_rows:
        print(f"\n  [{fx.nuclide}] {fx.geometry:>14} | {fx.label:>14} | "
              f"E_pp={pp_E:.1f} keV  S_pp={pp_area:.4e}")
        for m in matched:
            r_str = fmt_pct(m.ratio_to_photopeak) if m.feature != "photopeak" else "  (ref)  "
            print(f"     {m.feature:>14}  E={m.measured_E:>7.1f} (Δ={m.residual_keV:>+6.1f})  "
                  f"S={m.area:>9.0f}  FWHM={m.fwhm_keV:>5.1f}  R={r_str}")

    # === aggregator: mean / std of ratios per nuclide-feature ===
    print("\n" + "=" * 110)
    print("Закономерности: mean R = S_secondary / S_photopeak across geometries")
    print("=" * 110)
    print(f"{'Nuclide':>8}  {'Feature':>16}  {'E_theory':>9}  {'n':>3}  "
          f"{'mean R':>9}  {'std R':>9}  {'min R':>9}  {'max R':>9}  {'<dE>':>7}")
    print("-" * 110)
    cat: dict = {}   # JSON-serialisable catalog
    for (nuc, feat), entries in sorted(aggregator.items()):
        ratios = [e["ratio_to_photopeak"] for e in entries
                  if e["ratio_to_photopeak"] is not None]
        if not ratios:
            continue
        E_theory = entries[0]["expected_E_keV"]
        residuals = [e["residual_keV"] for e in entries]
        mean_r = sum(ratios) / len(ratios)
        var_r = sum((r - mean_r)**2 for r in ratios) / max(1, len(ratios)-1)
        std_r = math.sqrt(var_r)
        print(f"{nuc:>8}  {feat:>16}  {E_theory:>9.2f}  {len(ratios):>3}  "
              f"{mean_r:>9.4f}  {std_r:>9.4f}  "
              f"{min(ratios):>9.4f}  {max(ratios):>9.4f}  "
              f"{sum(residuals)/len(residuals):>+7.2f}")
        cat.setdefault(nuc, {"primary_E_keV": None, "features": []})
        cat[nuc]["features"].append({
            "name": feat,
            "expected_E_keV": E_theory,
            "n_observations": len(ratios),
            "mean_intensity_ratio": mean_r,
            "std_intensity_ratio": std_r,
            "min_intensity_ratio": min(ratios),
            "max_intensity_ratio": max(ratios),
            "mean_position_residual_keV": sum(residuals) / len(residuals),
            "observations": entries,
        })

    # Set primary E_gamma per nuclide for JSON header
    for fx, *_ in fixture_rows:
        cat.setdefault(fx.nuclide, {"primary_E_keV": fx.E_primary, "features": []})
        cat[fx.nuclide]["primary_E_keV"] = fx.E_primary

    # === write the catalog ===
    out_json = _OUT_SECONDARY_PEAKS_PATH
    out_json.parent.mkdir(parents=True, exist_ok=True)
    catalog = {
        "_schema": {
            "version": "0.1",
            "primary_E_keV": "energy of the primary photopeak gamma",
            "features": "list of secondary features; each has expected_E_keV, "
                        "mean/std/min/max intensity ratio S_sec / S_photopeak across "
                        "measured geometries, mean position residual (E_meas - E_theory)",
            "observations": "per-fixture entries that fed the statistics",
        },
        "_sources": {
            "physics": "Knoll 'Radiation Detection and Measurement' 4th Ed., "
                       "chap. 10 (Compton); chap. 11.A.5 (backscatter); chap. 12.B.2 "
                       "(escape peaks). Gilmore & Joss 3rd Ed., chap. 2 + chap. 6 (NaI artefacts).",
            "fixtures": "Gamma-1S NaI 63x63 USB #SN-01 reference measurements (2024 поверка)",
        },
        "nuclides": cat,
    }
    with out_json.open("w", encoding="utf-8") as f:
        json.dump(catalog, f, indent=2, ensure_ascii=False)
    print(f"\nCatalog saved: {out_json}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
