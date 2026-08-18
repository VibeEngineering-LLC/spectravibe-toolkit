from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""
Range/shape characterisation of problem isotopes for inference of their
presence in unknown samples (F-38 prep, v1.7.16).

This is the methodological refinement of v1.7.15's `analyze_secondaries.py`.
Where v1.7.15 reported mean +/- std per (nuclide, feature), v1.7.16
reports the full empirical RANGE — {min, p10, median, p90, max} of
position residual and intensity ratio — plus per-feature SHAPE
descriptors (FWHM, low-tail asymmetry index, edge-width for Compton
shoulders) so identification can ask not just "is there a peak near
X keV?" but "does the peak's position/shape match the range observed
across all known measurement conditions of this parent?".

Key methodological points the user emphasised:

  Position floats with isotope activity, geometry, and other
  conditions; characterise the range, not a single value.

  Identify the expected SHAPE of the parent's own photopeak (so we
  can confirm a candidate is a genuine photopeak vs a broader
  secondary), AND the shape of each secondary.

  Focus on problem isotopes whose secondaries can be confused with
  genuine photopeaks of other nuclides — those whose presence the
  catalog should help infer.

Problem-isotope set (rationale = secondary collides with another
nuclide's real gamma line):

  Cs-137    Compton edge ~478 keV  ↔ Bi-214 503 keV; Be-7 478 keV
            backscatter ~184 keV   ↔ U-235 185 keV
            Ba Ka IC X-rays ~32 keV (CONFIRMS Cs-137 vs Bi-214 609)

  K-40      Compton edge ~1243 keV ↔ Co-60 1173 keV !!! (close)
            backscatter ~217 keV   ↔ Pb-212 238.6 keV, U-235 185
            single escape ~950 keV ↔ Bi-214 934 keV

  Co-60     two photopeaks 1173, 1332 + sum coincidence 2506
            two Compton edges 963, 1118
            two backscatters 210, 214 (close together, often merge)

  Na-22     511 annihilation + 1274 photopeak
            sum coincidence 1786
            Compton edges 341, 1062

  Tl-208    multi-line cascade: 583, 2614, 510, 860
            backscatter 178, 232 (overlaps Pb-212 238.6)
            DE peak 1593 (from 2614)

  Bi-214    main 609 + many others (1764, 1120, 1238 ...)
            Compton edge 415  ↔ Ba-133 383

For each (parent, feature) the catalog records:

  position_keV : {min, p10, median, p90, max, n_obs, n_geometries}
  position_dependence : {by_geometry: {geo: [E_obs values]}}
  intensity_ratio : {min, p10, median, p90, max} of S_secondary / S_photopeak
  shape :
    photopeak : { FWHM_keV_range, asymmetry_range, gaussian_chi2 }
    compton_edge : { effective_FWHM_range, edge_drop_amplitude_range }
    backscatter : { FWHM_keV_range, high_E_tail_index }
  conflict_lines : [ {nuclide, library_E_keV, distance_keV} ... ]
                   real gamma lines from other nuclides within p10..p90
                   of this feature's observed position range

Methodology:
  - Each fixture is run through mariscotti_search at its own per-channel
    FWHM model.
  - Each expected feature gets the best Mariscotti hit within an
    energy-dependent tolerance (FWHM-scaled).
  - For photopeak: Gaussian fit on +/- 2.5*sigma window, report
    FWHM_keV (full width half maximum), and asymmetry =
    (area_left_half / area_total) − 0.5 (positive = right-skewed).
  - For Compton edge: width of the falling shoulder estimated from
    distance between half-plateau and quarter-plateau intensities.
  - For backscatter: Gaussian + high-E tail (multi-scatter), report
    FWHM and tail_index = (sum above peak+1*sigma) / (peak area).
  - Conflicts: scan `gamma.data.nuclide_library` for γ-lines from
    OTHER nuclides whose energy falls inside p10..p90 of our feature's
    observed range.

Output:
  - `data/secondary_peaks_v2.json` — rich catalog (consumed by future
    identification logic).
  - stdout summary tables per problem isotope.
"""


import json
import math
import re
import statistics
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Callable

sys.path.insert(0, str(Path(__file__).parent))

import numpy as np

from gamma.io.readers import read_spectrum
from gamma.peaks.search import mariscotti_search
from gamma.data.nuclide_library import (
    get_nuclide, list_nuclides, load_lsrm_chain_libs,
)
from gamma.physics.secondary_peaks import (
    compton_edge_keV, backscatter_keV, expected_features_for,
    M_E_KEV,
)
from validate_certs import make_lsrm_fwhm_provider

# F-39 (v1.7.17): supplement the in-memory library so Tl-208 / Pb-212 /
# Ac-228 records (Th-232 chain daughters) are visible to the conflict
# detector below.
load_lsrm_chain_libs()


# ---------------------------------------------------------------------------
# Inventory: which fixtures belong to which parent isotope
# ---------------------------------------------------------------------------

ROOT = Path(__file__).parent.parent
# F-83 (v1.12.0): Gamma-1S reference spectra live under detectors/Gamma-1S/.
from gamma.detectors.gamma1s import (
    DEFAULT_REFERENCE_DIR as REF,
    SECONDARY_PEAKS_V2_PATH as _OUT_SECONDARY_PEAKS_V2_PATH,
)
FIX = ROOT / "evals" / "fixtures"


def classify_geometry(p: Path) -> str:
    s = str(p)
    if "25см" in s or "25cm" in s.lower(): return "point25cm"
    if "марин" in s.lower() or "marinelli" in s.lower(): return "marinelli"
    if "дента" in s.lower() or "denta" in s.lower(): return "denta"
    if "петри" in s.lower() or "petri" in s.lower(): return "petri"
    if "5см" in s or "5cm" in s.lower(): return "point5cm"
    if "_0_см" in s or "_0см" in s: return "point0cm_xml"
    if p.name.startswith("M_"): return "M-source-legacy"
    return "unknown"


_NUCLIDE_RE = re.compile(
    r"(Am-?241|Ba-?133|Cd-?109|Ce-?139|Co-?57|Co-?60|Cs-?137|Cs137|"
    r"Eu-?152|Mn-?54|Na-?22|Th-?228|Th-?232|Th232|Y-?88|K40|K-?40|"
    r"Ra226|Ra-?226|Bi-?207|M_cs|M_k|M_ra|M_th)",
    re.IGNORECASE,
)
_NUC_NORM = {
    "AM241": "Am-241", "BA133": "Ba-133", "CD109": "Cd-109", "CE139": "Ce-139",
    "CO57": "Co-57", "CO60": "Co-60", "CS137": "Cs-137", "EU152": "Eu-152",
    "MN54": "Mn-54", "NA22": "Na-22", "TH228": "Th-228", "TH232": "Th-232",
    "Y88": "Y-88", "K40": "K-40", "RA226": "Ra-226", "BI207": "Bi-207",
    "M_CS": "Cs-137", "M_K": "K-40", "M_RA": "Ra-226", "M_TH": "Th-232",
}


def classify_nuclide(name: str) -> Optional[str]:
    m = _NUCLIDE_RE.search(name)
    if not m: return None
    key = m.group(1).upper().replace("-", "")
    return _NUC_NORM.get(key, key)


def gather_inventory() -> dict:
    """Map nuclide name -> list of (geometry, Path) for all fixtures."""
    inv: dict = {}
    for D in [REF, FIX]:
        if not D.is_dir(): continue
        for f in sorted(list(D.rglob("*.spe")) + list(D.rglob("*.xml"))):
            s = str(f).lower()
            if "фон" in s or "нестабильность" in s or "subtract" in s:
                continue
            nuc = classify_nuclide(f.name)
            if not nuc: continue
            inv.setdefault(nuc, []).append((classify_geometry(f), f))
    return inv


# ---------------------------------------------------------------------------
# Problem isotope catalog: which features to characterise per nuclide
# ---------------------------------------------------------------------------
#
# `primary_lines` lists every γ-line whose photopeak and secondaries we
# track. `conflict_keys` is the descriptive set of confusions this
# isotope generates downstream (informational; the actual conflicts are
# computed from the nuclide library scanning).

PROBLEM_ISOTOPES = {
    "Cs-137": {
        "primary_lines": [661.66],
        "conflict_note": "Compton edge collides with Bi-214 503; Be-7 478. "
                         "Backscatter collides with U-235 185.",
    },
    "K-40": {
        "primary_lines": [1460.82],
        "conflict_note": "Compton edge at ~1243 keV is dangerously close to Co-60 1173 keV. "
                         "Backscatter at ~217 keV collides with Pb-212 238.6 keV.",
    },
    "Co-60": {
        "primary_lines": [1173.23, 1332.49],
        "conflict_note": "Both photopeaks close (159 keV apart) so each "
                         "sits in the other's backscatter region.",
    },
    "Na-22": {
        "primary_lines": [511.00, 1274.54],
        "conflict_note": "511 keV annihilation is multi-source. "
                         "1274 keV Compton edge at ~1062 keV overlaps Bi-214 1120.",
    },
    "Y-88": {
        "primary_lines": [898.04, 1836.06],
        "conflict_note": "898 keV is well separated; 1836 keV escapes contaminate "
                         "Compton continuum from K-40.",
    },
    "Th-228": {
        # We characterise the daughter Tl-208 secondaries which dominate
        # the spectrum in secular equilibrium.
        "primary_lines": [583.19, 2614.51],
        "conflict_note": "Tl-208 583 backscatter (~178 keV) collides with "
                         "Pb-212 238.6 region. 2614 single escape (2104) is "
                         "the highest-energy clean line on NaI.",
    },
    # F-39: explicit chain daughters. Th-228 fixtures are aliased to
    # Tl-208 / Pb-212 / Ac-228 below in `_PARENT_ALIASES` so their
    # γ-lines (which dominate the Th-228 spectrum in secular
    # equilibrium) are characterised separately.
    "Tl-208": {
        "primary_lines": [583.19, 2614.51, 510.77, 860.56],
        "conflict_note": "510.77 ↔ Na-22 511 annihilation (always present "
                         "in positron-emitter spectra). 860 keV close to "
                         "Bi-214 768/806 doublet on NaI.",
    },
    "Pb-212": {
        "primary_lines": [238.63, 300.09],
        "conflict_note": "238.63 ↔ U-235 238 (Th-234 chain) and Cs-137 "
                         "backscatter region overlap.",
    },
    "Ac-228": {
        "primary_lines": [911.20, 968.97, 338.32, 463.00, 794.95, 209.25],
        "conflict_note": "911 ↔ Co-60 1173 Compton edge tail. 209 ↔ K-40 "
                         "backscatter and U-235 205 keV region.",
    },
}

# Some "problem isotopes" don't have their own dedicated fixtures: they
# appear inside chain-equilibrium parent sources. Map them to the
# physical fixture-providing parent so the inventory scan picks them up.
_PARENT_ALIASES = {
    "Tl-208": "Th-228",
    "Pb-212": "Th-228",
    "Ac-228": "Th-228",
}


# ---------------------------------------------------------------------------
# Per-spectrum analysis
# ---------------------------------------------------------------------------

@dataclass
class FeatureObservation:
    feature_name: str
    parent_nuclide: str
    primary_E_keV: float              # WHICH photopeak this secondary belongs to
    expected_E_keV: float             # theoretical position of this feature
    measured_E_keV: float
    residual_keV: float
    area: float
    fwhm_keV: float
    fwhm_theoretical_keV: float       # from the FWHM(E) provider at this E
    fwhm_ratio: float                 # measured / theoretical
    asymmetry: float                  # left-half-area/total - 0.5
    geometry: str
    fixture_label: str


def gauss_fit_local(counts: np.ndarray, peak_ch: int, fwhm_ch: float,
                    energy_at: Callable[[int], float]):
    """Local Gaussian + linear baseline fit on +/- 2.5*sigma window.

    Returns (E_centroid, FWHM_keV, area, asymmetry, residual_rms).
    On failure returns Nones.
    """
    sigma_ch = max(1.0, fwhm_ch / 2.355)
    roi = int(round(2.5 * sigma_ch))
    lo = max(0, peak_ch - roi)
    hi = min(counts.size, peak_ch + roi + 1)
    if hi - lo < 5: return None
    x = np.arange(lo, hi, dtype=np.float64)
    y = counts[lo:hi].astype(np.float64)
    # Linear baseline from a wider window (1*sigma outside the peak ROI)
    bg_l_lo = max(0, peak_ch - int(round(4 * sigma_ch)))
    bg_l_hi = max(bg_l_lo + 1, peak_ch - int(round(2.5 * sigma_ch)))
    bg_r_lo = min(counts.size - 1, peak_ch + int(round(2.5 * sigma_ch)))
    bg_r_hi = min(counts.size, peak_ch + int(round(4 * sigma_ch)))
    if bg_l_hi <= bg_l_lo or bg_r_hi <= bg_r_lo: return None
    bg_l = counts[bg_l_lo:bg_l_hi].mean()
    bg_r = counts[bg_r_lo:bg_r_hi].mean()
    # local linear baseline slope from the two sides
    x_l = 0.5 * (bg_l_lo + bg_l_hi - 1)
    x_r = 0.5 * (bg_r_lo + bg_r_hi - 1)
    if x_r != x_l:
        slope = (bg_r - bg_l) / (x_r - x_l)
    else:
        slope = 0.0
    intercept = bg_l - slope * x_l
    bg = slope * x + intercept
    net = np.maximum(0.0, y - bg)
    if net.sum() <= 0: return None
    # Gaussian via weighted moments
    total = float(net.sum())
    centroid_ch = float((net * x).sum() / total)
    var = float((net * (x - centroid_ch) ** 2).sum() / total)
    sigma_ch_fit = max(0.5, math.sqrt(max(var, 1.0)))
    fwhm_ch_fit = sigma_ch_fit * 2.355
    E_centroid = float(energy_at(centroid_ch))
    E_lo = float(energy_at(lo))
    E_hi = float(energy_at(hi - 1))
    if E_hi == E_lo: return None
    fwhm_keV = fwhm_ch_fit * (E_hi - E_lo) / (hi - lo - 1) * (hi - lo - 1)  # in keV
    # Simpler: fwhm_keV = fwhm_ch_fit * gain at centroid
    fwhm_keV = fwhm_ch_fit * abs(energy_at(int(round(centroid_ch + 0.5)))
                                  - energy_at(int(round(centroid_ch - 0.5))))
    # Asymmetry: fraction of area to the left of centroid - 0.5
    left_total = float(net[x < centroid_ch].sum())
    asym = left_total / total - 0.5
    return E_centroid, fwhm_keV, total, asym


def analyze_fixture(parent: str, E_lines: list, fx_path: Path, geometry: str) -> list:
    """Find every expected feature for `parent` in `fx_path` and return
    a list of FeatureObservation. One observation per (E_line, feature)."""
    try:
        spec = read_spectrum(str(fx_path))
    except Exception:
        return []
    fwhm_at = make_lsrm_fwhm_provider(spec)
    counts = np.asarray(spec.counts, dtype=np.float64)
    found = mariscotti_search(counts, fwhm_channels=fwhm_at, sigma_threshold=3.0)
    if not found: return []

    # Pre-compute energy for each found peak
    peak_E = [(p, spec.channel_to_energy(p.channel)) for p in found]

    observations: list = []
    # Photopeak area per primary line for ratio normalisation
    photopeak_area: dict = {}
    for E_p in E_lines:
        # Find best photopeak match
        best = None
        best_dE = 100.0
        for p, E in peak_E:
            dE = abs(E - E_p)
            if dE < best_dE:
                best = (p, E); best_dE = dE
        if best is None or best_dE > 30.0:
            continue
        p, E_meas = best
        fit = gauss_fit_local(counts, p.channel, p.fwhm_channels,
                              spec.channel_to_energy)
        if fit is None: continue
        E_c, fwhm_keV, area, asym = fit
        if area <= 0: continue
        fw_theo = fwhm_at(p.channel) * abs(spec.energy_cal[1]
                                            if len(spec.energy_cal) > 1 else 1.0)
        observations.append(FeatureObservation(
            feature_name="photopeak", parent_nuclide=parent,
            primary_E_keV=E_p,
            expected_E_keV=E_p, measured_E_keV=E_c,
            residual_keV=E_c - E_p, area=area,
            fwhm_keV=fwhm_keV, fwhm_theoretical_keV=fw_theo,
            fwhm_ratio=fwhm_keV / max(1.0, fw_theo), asymmetry=asym,
            geometry=geometry, fixture_label=fx_path.name,
        ))
        photopeak_area[E_p] = area

    # For each E_line: scan its expected secondaries
    used = {o.measured_E_keV for o in observations}
    for E_p in E_lines:
        if E_p not in photopeak_area: continue
        S_pp = photopeak_area[E_p]
        feats = expected_features_for(parent, E_p)
        for ef in feats:
            if ef.name == "photopeak": continue
            best = None; best_dE = ef.tolerance_keV
            for p, E in peak_E:
                if E in used: continue
                dE = abs(E - ef.E_keV)
                if dE < best_dE:
                    best = (p, E); best_dE = dE
            if best is None: continue
            p, E_meas = best
            used.add(E_meas)
            fit = gauss_fit_local(counts, p.channel, p.fwhm_channels,
                                  spec.channel_to_energy)
            if fit is None: continue
            E_c, fwhm_keV, area, asym = fit
            fw_theo = fwhm_at(p.channel) * abs(spec.energy_cal[1]
                                                if len(spec.energy_cal) > 1 else 1.0)
            observations.append(FeatureObservation(
                feature_name=ef.name, parent_nuclide=parent,
                primary_E_keV=E_p,
                expected_E_keV=ef.E_keV, measured_E_keV=E_c,
                residual_keV=E_c - ef.E_keV,
                area=area * (1.0 / max(1e-9, S_pp)),  # store as RATIO directly
                fwhm_keV=fwhm_keV, fwhm_theoretical_keV=fw_theo,
                fwhm_ratio=fwhm_keV / max(1.0, fw_theo), asymmetry=asym,
                geometry=geometry, fixture_label=fx_path.name,
            ))

    return observations


# ---------------------------------------------------------------------------
# Aggregation: per (parent, feature) quantile statistics
# ---------------------------------------------------------------------------

def quantile(values, q: float) -> float:
    """Linear-interp quantile, q in [0, 1]."""
    if not values: return float("nan")
    a = sorted(values)
    pos = q * (len(a) - 1)
    lo = int(math.floor(pos)); hi = int(math.ceil(pos))
    if lo == hi: return a[lo]
    return a[lo] + (a[hi] - a[lo]) * (pos - lo)


def summarise(values: list) -> dict:
    if not values:
        return {"n": 0, "min": None, "p10": None, "median": None,
                "p90": None, "max": None, "mean": None, "std": None}
    return {
        "n": len(values),
        "min": min(values),
        "p10": quantile(values, 0.10),
        "median": quantile(values, 0.50),
        "p90": quantile(values, 0.90),
        "max": max(values),
        "mean": statistics.mean(values),
        "std": statistics.pstdev(values) if len(values) > 1 else 0.0,
    }


# ---------------------------------------------------------------------------
# Conflict detection
# ---------------------------------------------------------------------------

def find_conflicts(parent: str, E_range: tuple, min_I_pct: float = 1.0) -> list:
    """Return real γ-lines from OTHER nuclides falling inside `E_range`."""
    lo, hi = E_range
    conflicts = []
    for nuc in list_nuclides():
        if nuc == parent: continue
        rec = get_nuclide(nuc) or {}
        for line in rec.get("lines", []):
            E = float(line[0]); I = float(line[1]) if len(line) > 1 else 0.0
            if I < min_I_pct: continue
            if lo <= E <= hi:
                conflicts.append({
                    "nuclide": nuc, "library_E_keV": E,
                    "library_I_pct": I,
                    "distance_to_range_lo_keV": E - lo,
                    "distance_to_range_hi_keV": hi - E,
                })
    conflicts.sort(key=lambda c: -c["library_I_pct"])
    return conflicts


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    inv = gather_inventory()
    print(f"=== Inventory ===")
    for nuc, items in sorted(inv.items()):
        if nuc not in PROBLEM_ISOTOPES: continue
        geo_count: dict = {}
        for g, _ in items:
            geo_count[g] = geo_count.get(g, 0) + 1
        print(f"  {nuc:>7}: {len(items):>2} fixtures   "
              + "  ".join(f"{g}:{c}" for g, c in sorted(geo_count.items())))

    # Run analysis per problem isotope
    # Key by (parent, primary_E_keV, feature) so multi-line nuclides
    # (Co-60 1173/1332, Th-228 583/2614, Na-22 511/1274) don't conflate
    # different secondaries into one bucket.
    all_obs: dict = {}
    for parent, cfg in PROBLEM_ISOTOPES.items():
        # F-39: chain daughters (Tl-208, Pb-212, Ac-228) inherit fixtures
        # from their physical-source parent (Th-228 sealed source).
        fixture_source = _PARENT_ALIASES.get(parent, parent)
        fixtures = inv.get(fixture_source, [])
        if not fixtures:
            print(f"\n!  {parent}: no fixtures available")
            continue
        for geo, p in fixtures:
            obs = analyze_fixture(parent, cfg["primary_lines"], p, geo)
            for o in obs:
                key = (o.parent_nuclide, round(o.primary_E_keV, 2),
                       o.feature_name)
                all_obs.setdefault(key, []).append(o)

    # Build catalog
    catalog: dict = {
        "_schema": {
            "version": "0.2",
            "description": "Range/shape characterisation of problem isotopes",
            "fields": {
                "primary_lines": "gamma-line photopeaks tracked for this parent",
                "features": {
                    "<feature_name>": {
                        "expected_E_keV": "theoretical/library position",
                        "position_keV": "quantile summary of observed centroid energy across fixtures",
                        "position_residual_keV": "observed - expected, quantile summary",
                        "intensity_ratio": "S_secondary / S_photopeak quantile summary",
                        "fwhm_keV": "measured FWHM quantile summary",
                        "fwhm_ratio": "measured / theoretical FWHM (>>1 = wider than photopeak = backscatter, etc.)",
                        "asymmetry": "left-half-area/total - 0.5; positive = right-skewed",
                        "by_geometry": "per-geometry observations",
                        "conflict_lines": "real gamma lines from other nuclides within p10..p90 position range",
                    }
                }
            }
        },
        "nuclides": {},
    }

    print("\n" + "=" * 110)
    print("RANGE/SHAPE catalog per problem isotope (per primary line)")
    print("=" * 110)
    for (parent, primary_E, feature), obs_list in sorted(all_obs.items()):
        positions = [o.measured_E_keV for o in obs_list]
        residuals = [o.residual_keV for o in obs_list]
        ratios = [o.area for o in obs_list]  # already stored as ratio for secondaries
        fwhms = [o.fwhm_keV for o in obs_list]
        fwhm_ratios = [o.fwhm_ratio for o in obs_list]
        asyms = [o.asymmetry for o in obs_list]
        expected_E = obs_list[0].expected_E_keV
        n_geometries = len({o.geometry for o in obs_list})
        position_p10 = quantile(positions, 0.10)
        position_p90 = quantile(positions, 0.90)

        cat_nuc = catalog["nuclides"].setdefault(parent, {
            "primary_lines": PROBLEM_ISOTOPES[parent]["primary_lines"],
            "conflict_note": PROBLEM_ISOTOPES[parent]["conflict_note"],
            "by_primary_line": {},
        })
        primary_key = f"{primary_E:.2f}"
        primary_entry = cat_nuc["by_primary_line"].setdefault(primary_key, {
            "features": {},
        })

        # Conflicts: real γ-lines within p10..p90 observed range
        conflicts = []
        if feature != "photopeak":
            conflicts = find_conflicts(parent, (position_p10, position_p90),
                                       min_I_pct=1.0)

        entry = {
            "expected_E_keV": expected_E,
            "position_keV": summarise(positions),
            "position_residual_keV": summarise(residuals),
            "intensity_ratio": (summarise(ratios)
                                if feature != "photopeak" else None),
            "fwhm_keV": summarise(fwhms),
            "fwhm_ratio_to_theory": summarise(fwhm_ratios),
            "asymmetry": summarise(asyms),
            "n_geometries": n_geometries,
            "by_geometry": {},
            "conflict_lines": conflicts,
        }
        for o in obs_list:
            entry["by_geometry"].setdefault(o.geometry, []).append({
                "fixture": o.fixture_label,
                "measured_E_keV": o.measured_E_keV,
                "ratio_to_photopeak": (o.area if feature != "photopeak" else None),
                "fwhm_keV": o.fwhm_keV,
                "asymmetry": o.asymmetry,
            })
        primary_entry["features"][feature] = entry

        # Stdout summary
        print(f"\n  [{parent} from {primary_E:.1f}] {feature:>14}  E_exp={expected_E:.2f}  "
              f"n={len(obs_list)}/{n_geometries}geo")
        print(f"     position    : "
              f"[{summarise(positions)['min']:>7.1f} .. p10={position_p10:>7.1f} .. "
              f"med={summarise(positions)['median']:>7.1f} .. p90={position_p90:>7.1f} .. "
              f"{summarise(positions)['max']:>7.1f}]  std={summarise(positions)['std']:.2f}")
        print(f"     residual    : "
              f"min={summarise(residuals)['min']:>+6.1f}  med={summarise(residuals)['median']:>+6.1f}  "
              f"max={summarise(residuals)['max']:>+6.1f}")
        if feature != "photopeak":
            print(f"     ratio S/S_pp: "
                  f"min={summarise(ratios)['min']:.4f}  med={summarise(ratios)['median']:.4f}  "
                  f"max={summarise(ratios)['max']:.4f}")
        print(f"     FWHM_keV    : "
              f"[{summarise(fwhms)['min']:>5.1f} .. med={summarise(fwhms)['median']:>5.1f} .. "
              f"{summarise(fwhms)['max']:>5.1f}]  ratio_to_theory med="
              f"{summarise(fwhm_ratios)['median']:.2f}")
        if conflicts[:3]:
            for c in conflicts[:3]:
                print(f"     CONFLICT    : {c['nuclide']:>8} @ {c['library_E_keV']:>7.2f} "
                      f"keV (I={c['library_I_pct']:.1f}%)")

    # Save catalog
    out_json = _OUT_SECONDARY_PEAKS_V2_PATH
    with out_json.open("w", encoding="utf-8") as f:
        json.dump(catalog, f, indent=2, ensure_ascii=False)
    print(f"\nCatalog saved: {out_json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
