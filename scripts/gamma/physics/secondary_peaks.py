"""
Reference-sample-derived catalog of secondary spectral features.

A "secondary peak" is anything that appears in a γ-spectrum but is NOT
a photopeak of a real nuclide γ-emission: Compton-edge maxima,
backscatter peaks (single 180° scatter off shielding), single/double
escape peaks from pair production, iodine K X-ray escape from the NaI
crystal, IC X-rays from the decaying nuclide itself, and the always-
present natural K-40 background line.

Identifying secondary features explicitly serves two purposes:

  1. **Anti-misidentification**. A peak detected at ~480 keV on a
     Cs-137-dominated spectrum is overwhelmingly the Cs-137 Compton
     edge, not a Bi-214 503 keV line. The identification module can
     ask `expected_secondary_position(parent="Cs-137", feature=
     "compton_edge")` to know what to NOT credit as a new nuclide.

  2. **Cross-validation**. A claimed Cs-137 photopeak that has NO
     detectable backscatter peak at ~190 keV and NO IC X-ray at ~28
     keV is more likely a Bi-214 609 keV misidentification (especially
     on NaI). The presence of expected secondaries CONFIRMS the
     parent.

Two access paths:

    from gamma.physics.secondary_peaks import (
        compton_edge_keV, backscatter_keV,
        expected_features_for,  # theoretical positions
        load_catalog,            # measured statistics across geometries
    )

The theoretical functions are pure physics and have no dependencies.
The empirical catalog reads `detectors/Gamma-1S/data/secondary_peaks.json`
(built by `analyze_secondaries.py` from real Gamma-1S NaI 63×63 measurements).
The path is resolved via `gamma.detectors.gamma1s.SECONDARY_PEAKS_PATH` —
no other detector subtree currently supplies a calibrated catalog.

═══════════════════════════════════════════════════════════════════════
Empirical observations on the Gamma-1S NaI 63x63 fixture set (v1.7.15)
═══════════════════════════════════════════════════════════════════════

Across 17 reference fixtures (10 Cs-137 + 7 K-40 at 4 distinct
geometries), the patterns from measured spectra are:

| Feature                 | E_theory  | <E_meas − E_theory> | mean R = S/S_pp |
|-------------------------|----------:|--------------------:|----------------:|
| Cs-137 backscatter      | 184.3 keV |       +8.1 keV      |     7.3 % ±2.9% |
| Cs-137 Compton edge     | 477.3 keV |      −37.0 keV      |     3.0 % ±1.0% |
| Cs-137 Ba Kα IC X-ray   |  32.0 keV |       −5.9 keV      |     8.4 % ±4.1% |
| K-40   backscatter      | 217.5 keV |      +14.9 keV      |    12.2 % ±3.9% |
| K-40   Compton edge     |1243.4 keV |      −53.2 keV      |     7.7 % ±2.6% |
| K-40   single escape    | 949.8 keV |       +5.2 keV      |     2.9 % ±0.8% |

Three robust patterns:

(A) The **Compton edge centroid sits BELOW** its analytical position
    by ~35–55 keV on NaI 63×63. This is not a calibration drift — it
    affects both Cs-137 and K-40 identically. Root cause is the
    Mariscotti algorithm: the analytical Compton edge is a step
    discontinuity, broadened by detector resolution into a falling
    shoulder. The peak-finder picks the local maximum of the second-
    derivative response, which sits below the analytical step by
    approximately 0.7·FWHM. For Cs-137 at 477 keV, FWHM ≈ 50 keV →
    expected shift ≈ −35 keV (observed: −37). For K-40 at 1243 keV,
    FWHM ≈ 75 keV → expected shift ≈ −52 keV (observed: −53). The
    rule of thumb: **E_Compton_observed ≈ E_Compton_theory − 0.7·FWHM**.

(B) The **backscatter peak sits ABOVE** its analytical position by
    ~+8 to +15 keV. The shift correlates with extended-source
    geometry (worst on Дента/Петри/Маринелли, smallest on point 25 cm),
    consistent with multi-path / multi-scatter contributions adding
    to the single-180° base. The intensity ratio S_bs / S_pp is the
    cleanest geometry diagnostic — **>10% at close point geometry,
    ~3% at distant point geometry, ~7% at extended-source containers**.

(C) **Natural K-40 background contaminates every long-integration
    Cs-137 spectrum** at 0.3–10 % of the Cs-137 photopeak area
    (depending on container material, sample mass, and integration
    time). The peak appears at 1444–1472 keV with calibration drift;
    identification should EXPECT it and never credit it to a
    spuriously-added K-40 source unless the magnitude vastly exceeds
    expected background level.

Methodology references:
  - Knoll, "Radiation Detection and Measurement" 4th Ed., chap. 10
    (Compton scattering kinematics), chap. 11.A.5 (backscatter peak
    formation), chap. 12.B.2 (single/double escape).
  - Gilmore & Joss, "Practical Gamma-ray Spectrometry" 3rd Ed.,
    chap. 2 (interactions) and chap. 6 (NaI(Tl) artefacts).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from gamma.data import DATA_DIR


# Electron rest mass energy [keV]. Codata-22 to 6 sig fig.
M_E_KEV = 510.999


def compton_edge_keV(E_keV: float) -> float:
    """Analytical Compton edge (max electron KE at 180° photon scatter).

    On NaI the observed peak sits ~0.7·FWHM BELOW this value. See
    `compton_edge_observed_keV` for the empirically-corrected position.
    """
    return 2.0 * E_keV * E_keV / (M_E_KEV + 2.0 * E_keV)


def backscatter_keV(E_keV: float) -> float:
    """Analytical 180° backscatter photon energy.

    On NaI the observed peak sits ~+8 to +15 keV ABOVE this value due
    to multi-path contributions in extended sources. See
    `backscatter_observed_keV` for the empirically-corrected position.
    """
    return E_keV / (1.0 + 2.0 * E_keV / M_E_KEV)


def compton_edge_observed_keV(E_keV: float, fwhm_at_edge_keV: float) -> float:
    """Empirical: Compton edge observed centroid on NaI scintillators.

    Rule of thumb from the Gamma-1S catalogue:
        E_obs ≈ E_theory − 0.7·FWHM(E_theory)
    """
    return compton_edge_keV(E_keV) - 0.7 * fwhm_at_edge_keV


def backscatter_observed_keV(E_keV: float, geometry: str = "point_far") -> float:
    """Empirical: backscatter peak observed centroid on NaI scintillators.

    Geometry-conditional shift from the Gamma-1S catalogue:
        point_5cm / point_close   → +14 keV
        extended_source           → +10 keV
        point_25cm / point_far    → +5 keV
    """
    shifts = {
        "point_5cm":        14.0,
        "point_close":      14.0,
        "extended_source":  10.0,
        "marinelli":        10.0,
        "denta":            10.0,
        "petri":            10.0,
        "point_25cm":        5.0,
        "point_far":         5.0,
    }
    return backscatter_keV(E_keV) + shifts.get(geometry, 10.0)


@dataclass(frozen=True)
class ExpectedFeature:
    name: str
    E_keV: float
    tolerance_keV: float = 30.0
    note: str = ""


def expected_features_for(nuclide: str, E_gamma: float) -> list:
    """Build the theoretical secondary-feature list for one γ-line.

    Used by `analyze_secondaries.py` AND by future identification-side
    "anti-misidentification" logic. The same function emits BOTH so the
    cataloguing and the consuming logic stay consistent.
    """
    feats = [
        ExpectedFeature("photopeak",    E_gamma, 30.0, "primary photopeak"),
        ExpectedFeature("compton_edge", compton_edge_keV(E_gamma), 60.0,
                        "max electron KE from 180° Compton scatter"),
        ExpectedFeature("backscatter",  backscatter_keV(E_gamma), 40.0,
                        "180° backscatter photon off shielding"),
    ]
    if E_gamma > 1022.0 + 100.0:
        feats.append(ExpectedFeature("single_escape", E_gamma - 511.0, 40.0,
                                     "one annihilation γ escapes"))
        feats.append(ExpectedFeature("double_escape", E_gamma - 1022.0, 40.0,
                                     "both annihilation γs escape"))
    if E_gamma < 200.0:
        feats.append(ExpectedFeature("xray_escape", E_gamma - 28.0, 15.0,
                                     "iodine K X-ray escape from NaI"))
    if nuclide == "Cs-137":
        feats.append(ExpectedFeature("ic_xray_Ba_Ka", 32.0, 10.0,
                                     "Ba K-α X-rays from 8% IC branch"))
        feats.append(ExpectedFeature("k40_natural", 1460.82, 30.0,
                                     "natural-background K-40"))
    return feats


# ---------------------------------------------------------------------------
# Empirical catalog (built by analyze_secondaries.py from real spectra)
# ---------------------------------------------------------------------------

# F-83 (v1.12.0): secondary-peak catalogues were calibrated on Gamma-1S
# reference spectra and live under detectors/Gamma-1S/data/. The resolver
# module is the single source of truth for these paths.
from gamma.detectors.gamma1s import (
    SECONDARY_PEAKS_PATH as _CATALOG_PATH,
    SECONDARY_PEAKS_V2_PATH as _CATALOG_V2_PATH,
)
_CATALOG_CACHE: Optional[dict] = None
_CATALOG_V2_CACHE: Optional[dict] = None


def load_catalog() -> dict:
    """Load the empirical secondary-peaks catalogue, cached.

    Returns:
        dict with shape:
        ```
        {
            "_schema": {...},
            "_sources": {...},
            "nuclides": {
                "Cs-137": {
                    "primary_E_keV": 661.66,
                    "features": [
                        {
                            "name": "backscatter",
                            "expected_E_keV": 184.3,
                            "n_observations": 10,
                            "mean_intensity_ratio": 0.0726,
                            "std_intensity_ratio": 0.0297,
                            "mean_position_residual_keV": +8.12,
                            "observations": [...]
                        },
                        ...
                    ]
                },
                "K-40": {...}
            }
        }
        ```
    """
    global _CATALOG_CACHE
    if _CATALOG_CACHE is not None:
        return _CATALOG_CACHE
    if not _CATALOG_PATH.is_file():
        # Empty catalog -- consumer code falls back to theoretical formulas
        _CATALOG_CACHE = {"_schema": {}, "_sources": {}, "nuclides": {}}
        return _CATALOG_CACHE
    with _CATALOG_PATH.open("r", encoding="utf-8") as fh:
        _CATALOG_CACHE = json.load(fh)
    return _CATALOG_CACHE


def reset_catalog_cache() -> None:
    """Clear in-memory cache (useful for tests)."""
    global _CATALOG_CACHE
    _CATALOG_CACHE = None


def empirical_ratio(nuclide: str, feature: str) -> Optional[dict]:
    """Look up the measured intensity ratio for a (nuclide, feature) pair.

    Returns:
        Dict with keys `mean`, `std`, `min`, `max`, `n_observations`,
        `mean_residual_keV`, or None if no data.
    """
    cat = load_catalog()
    nucs = cat.get("nuclides", {})
    rec = nucs.get(nuclide, {})
    for f in rec.get("features", []):
        if f.get("name") == feature:
            return {
                "mean": f.get("mean_intensity_ratio"),
                "std":  f.get("std_intensity_ratio"),
                "min":  f.get("min_intensity_ratio"),
                "max":  f.get("max_intensity_ratio"),
                "n_observations": f.get("n_observations"),
                "mean_residual_keV": f.get("mean_position_residual_keV"),
            }
    return None


def load_catalog_v2() -> dict:
    """Load the v0.2 range/shape catalogue built by `analyze_problem_isotopes.py`.

    The v2 schema characterises **range** (min, p10, median, p90, max)
    instead of mean +/- std. It is keyed by `(nuclide, primary_E_keV,
    feature)` so multi-line parents (Co-60 1173/1332, Th-228 583/2614,
    Na-22 511/1274) keep their secondaries separated.

    Per (parent, primary_E, feature) record:
      - position_keV        : observed centroid quantile summary
      - position_residual_keV : observed - theoretical quantile summary
      - intensity_ratio     : S_secondary / S_photopeak quantile summary (None for photopeak)
      - fwhm_keV            : measured FWHM quantile summary
      - fwhm_ratio_to_theory: measured / theoretical FWHM (>1.2 typical for backscatter)
      - asymmetry           : left-half-area/total - 0.5 quantile summary
      - n_geometries        : how many geometries contributed
      - by_geometry         : per-geometry raw observations
      - conflict_lines      : real gamma-lines from OTHER nuclides falling in
                              p10..p90 of this feature's observed position range

    Returns:
        dict with keys `_schema`, `nuclides` — same shape as the JSON.
    """
    global _CATALOG_V2_CACHE
    if _CATALOG_V2_CACHE is not None:
        return _CATALOG_V2_CACHE
    if not _CATALOG_V2_PATH.is_file():
        _CATALOG_V2_CACHE = {"_schema": {}, "nuclides": {}}
        return _CATALOG_V2_CACHE
    with _CATALOG_V2_PATH.open("r", encoding="utf-8") as fh:
        _CATALOG_V2_CACHE = json.load(fh)
    return _CATALOG_V2_CACHE


def reset_catalog_v2_cache() -> None:
    global _CATALOG_V2_CACHE
    _CATALOG_V2_CACHE = None


def position_range(nuclide: str, primary_E_keV: float, feature: str,
                   *, span: str = "p10p90") -> Optional[tuple]:
    """Empirical position range for a (nuclide, primary line, feature).

    Args:
        span: which range to return.
            "minmax"   -> (min, max) of observed positions
            "p10p90"   -> (p10, p90) — recommended for inference (90% CI)
            "iqr"      -> (p25, p75) — robust narrow band

    Returns:
        (E_low_keV, E_high_keV) or None if the catalog has no entry.
    """
    cat = load_catalog_v2()
    rec = cat.get("nuclides", {}).get(nuclide, {})
    by_pl = rec.get("by_primary_line", {})
    key = f"{primary_E_keV:.2f}"
    if key not in by_pl:
        # Try wider key-match tolerance
        for k in by_pl:
            try:
                if abs(float(k) - primary_E_keV) < 0.5:
                    key = k; break
            except ValueError:
                continue
        else:
            return None
    feat = by_pl[key].get("features", {}).get(feature)
    if not feat:
        return None
    pos = feat.get("position_keV", {})
    if span == "minmax":
        return (pos.get("min"), pos.get("max"))
    if span == "p10p90":
        return (pos.get("p10"), pos.get("p90"))
    if span == "iqr":
        # Not pre-computed; fall back to median +/- 0.5*(p90-p10)
        med = pos.get("median"); p10 = pos.get("p10"); p90 = pos.get("p90")
        if med is None or p10 is None or p90 is None:
            return None
        half = 0.25 * (p90 - p10)
        return (med - half, med + half)
    return None


def matches_secondary(parent: str, observed_E_keV: float,
                      *, feature: Optional[str] = None,
                      span: str = "p10p90") -> list:
    """Is `observed_E_keV` consistent with one of `parent`'s known secondaries?

    Used by identification logic to demote a candidate nuclide if its
    only matched line falls inside the observed range of a parent's
    known secondary (e.g. Bi-214 503 keV vs Cs-137 Compton edge).

    Args:
        parent: the already-identified nuclide whose secondaries we check.
        observed_E_keV: energy of the candidate peak to test.
        feature: limit the check to one specific feature (default: any).
        span: "p10p90" (recommended) or "minmax".

    Returns:
        List of `{primary_E_keV, feature, range, distance_keV}` dicts —
        one entry per (primary_line, feature) of `parent` that contains
        `observed_E_keV` in its position range. Empty list = no match.
    """
    cat = load_catalog_v2()
    rec = cat.get("nuclides", {}).get(parent, {})
    out = []
    for primary_key, pdata in rec.get("by_primary_line", {}).items():
        for feat_name, fdata in pdata.get("features", {}).items():
            if feature is not None and feat_name != feature:
                continue
            pos = fdata.get("position_keV", {})
            if span == "minmax":
                lo, hi = pos.get("min"), pos.get("max")
            else:
                lo, hi = pos.get("p10"), pos.get("p90")
            if lo is None or hi is None:
                continue
            if lo <= observed_E_keV <= hi:
                median = pos.get("median")
                out.append({
                    "primary_E_keV": float(primary_key),
                    "feature": feat_name,
                    "range": (lo, hi),
                    "median_E_keV": median,
                    "distance_to_median_keV": (observed_E_keV - median
                                                if median else None),
                })
    return out


__all__ = [
    "M_E_KEV",
    "compton_edge_keV", "backscatter_keV",
    "compton_edge_observed_keV", "backscatter_observed_keV",
    "ExpectedFeature", "expected_features_for",
    "load_catalog", "reset_catalog_cache", "empirical_ratio",
    "load_catalog_v2", "reset_catalog_v2_cache",
    "position_range", "matches_secondary",
]
