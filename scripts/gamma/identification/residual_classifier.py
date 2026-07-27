"""
Residual peak classifier (F-74 / v1.11.1).

Once Stage-1 identification is done, the list of `unmatched_peaks`
contains everything Mariscotti found that wasn't claimed by any
nuclide's library line. But on NaI 63×63 the *majority* of such peaks
are well-understood physics, not new nuclides:

  • elemental K-shell XRF from Pb shielding, Th matrix, U matrix
    (the Th K-XRF at ~93/106 keV is dominant in any Th-rich sample)
  • backscatter / Compton edges / single+double escapes of the
    strongest detected primary photopeaks — catalogued empirically in
    `detectors/Gamma-1S/data/secondary_peaks_v2.json` (Lsrm methodology, F-40)
  • sum peaks E_a + E_b where both parents are already in the spectrum
  • Lsrm-software's own ROI residuals near the energy ceiling

Only AFTER subtracting these explainable categories does the residual
list represent peaks that genuinely cannot be accounted for by Stage 1
ЕРН. That count is what should drive the recommendation to escalate
to Stage 2 (technogenic) — not the raw σ≥4 count.

Per user methodology (15.11.2025): "Лучше не фантазировать, а
переспросить пользователя." A spectrum where every residual is XRF or
chain-secondary should NOT trigger Stage 2 escalation.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from gamma.physics.secondary_peaks import matches_secondary
from gamma.data.nuclide_library import get_nuclide


# ──────────────────────────────────────────────────────────────────
# XRF catalog loading
# ──────────────────────────────────────────────────────────────────

_XRF_PATH = Path(__file__).resolve().parents[3] / "data" / "xrf_lines.json"


@lru_cache(maxsize=1)
def _load_xrf_catalog() -> Dict[str, Dict]:
    if not _XRF_PATH.is_file():
        return {}
    with _XRF_PATH.open("r", encoding="utf-8") as f:
        d = json.load(f)
    return {k: v for k, v in d.items() if not k.startswith("_")}


def _match_xrf(peak_E: float, tolerance_keV: float) -> List[Tuple[str, str, float, float]]:
    """
    For each element in the catalog, check if peak_E matches any of its
    K-shell (Kα1, Kα2, Kβ1, Kβ2) or L-shell lines within tolerance.

    Returns list of (element, shell, library_E_keV, residual_keV).
    """
    cat = _load_xrf_catalog()
    out = []
    for el, info in cat.items():
        if not isinstance(info, dict):
            continue
        for shell in ("K", "L"):
            lines = info.get(shell) or []
            for lib_E in lines:
                if lib_E is None:
                    continue
                try:
                    lib_E_f = float(lib_E)
                except (TypeError, ValueError):
                    continue
                d = abs(peak_E - lib_E_f)
                if d <= tolerance_keV:
                    out.append((el, shell, lib_E_f, d))
    return out


# ──────────────────────────────────────────────────────────────────
# Result dataclass
# ──────────────────────────────────────────────────────────────────

# Classification labels — explicit set for downstream switch logic.
LBL_XRF = "xrf"
LBL_CHAIN_SECONDARY = "chain_secondary"
LBL_SINGLE_ESCAPE = "single_escape"
LBL_DOUBLE_ESCAPE = "double_escape"
LBL_SUM_PEAK = "sum_peak"
LBL_ANNIHILATION = "annihilation_511"
LBL_EDGE_OF_RANGE = "edge_of_range"
LBL_TRUE_UNMATCHED = "true_unmatched"


@dataclass
class ResidualClassification:
    peak_E_keV: float
    sigma: float
    label: str                                  # one of LBL_*
    note: str                                   # human-readable description
    parent_nuclide: str = ""                    # for chain_secondary / escape / sum
    parent_line_keV: Optional[float] = None
    feature_kind: str = ""                      # backscatter / compton_edge / etc.
    element: str = ""                           # for xrf
    delta_keV: Optional[float] = None


# ──────────────────────────────────────────────────────────────────
# Classification logic
# ──────────────────────────────────────────────────────────────────

def _strong_lines_of(nuclide_name: str, intensity_floor: float = 5.0) -> List[float]:
    """Library lines of `nuclide_name` with I ≥ intensity_floor (in %)."""
    nuc = get_nuclide(nuclide_name)
    if not nuc:
        return []
    return [
        float(L[0]) for L in nuc.get("lines", [])
        if len(L) > 1 and float(L[1]) >= intensity_floor
    ]


def classify_residual(
    peak_E_keV: float,
    sigma: float,
    *,
    detected_nuclide_names: List[str],
    fwhm_at_keV: float,
    energy_max_keV: float = 3000.0,
    annihilation_tolerance_keV: Optional[float] = None,
) -> ResidualClassification:
    """
    Classify a single unmatched peak.

    Order of checks (first hit wins):
      1. Annihilation 511 ± FWHM
      2. Single escape (E_parent − 511) for each strong detected line
      3. Double escape (E_parent − 1022) for each strong detected line
      4. Sum peak  E_a + E_b for any two strong detected lines
      5. Chain-secondary (backscatter, compton_edge, sum, ce_band) per
         `secondary_peaks_v2` catalog
      6. XRF — Pb / Th / U / W / Sn / Cd / Ba K-shell lines
      7. Near energy ceiling (> energy_max_keV − 1·FWHM) → likely
         truncation / Lsrm-software artefact
      8. else → true_unmatched

    Parameters
    ----------
    peak_E_keV : float
        Centroid energy of the residual peak.
    sigma : float
        Mariscotti significance of the peak.
    detected_nuclide_names : list[str]
        Names of nuclides already detected in Stage 1 (post-disambiguate).
    fwhm_at_keV : float
        FWHM at `peak_E_keV` used as matching tolerance.
    energy_max_keV : float
        Upper edge of the kept channel range — peaks within ~1 FWHM of
        the ceiling are flagged as edge_of_range.
    annihilation_tolerance_keV : float, optional
        Annihilation tolerance override. Default = max(FWHM/2, 8 keV)
        — broad to absorb both the Doppler-broadened 511 and the
        adjacent Tl-208 510.77 photopeak.

    Returns
    -------
    ResidualClassification
    """
    tol = max(fwhm_at_keV * 0.5, 8.0)
    ann_tol = annihilation_tolerance_keV \
        if annihilation_tolerance_keV is not None else max(fwhm_at_keV * 0.6, 10.0)

    # 1. Annihilation 511
    if abs(peak_E_keV - 511.0) <= ann_tol:
        return ResidualClassification(
            peak_E_keV=peak_E_keV, sigma=sigma,
            label=LBL_ANNIHILATION,
            note=(f"511 кэВ — аннигиляция e⁺e⁻ от μ-мезонов космики "
                  f"(всегда присутствует на NaI вне свинцовой защиты)"),
            parent_line_keV=511.0,
            delta_keV=abs(peak_E_keV - 511.0),
        )

    # Build the list of strong detected library lines once
    parent_strong_lines: List[Tuple[str, float]] = []
    for n in detected_nuclide_names:
        for E in _strong_lines_of(n, intensity_floor=5.0):
            parent_strong_lines.append((n, E))

    # 2. Single escape
    for parent, E_p in parent_strong_lines:
        E_se = E_p - 511.0
        if E_se > 100 and abs(peak_E_keV - E_se) <= tol:
            return ResidualClassification(
                peak_E_keV=peak_E_keV, sigma=sigma,
                label=LBL_SINGLE_ESCAPE,
                note=f"single escape от {parent} {E_p:.0f} кэВ (E−511={E_se:.0f})",
                parent_nuclide=parent, parent_line_keV=E_p,
                feature_kind="single_escape",
                delta_keV=abs(peak_E_keV - E_se),
            )

    # 3. Double escape
    for parent, E_p in parent_strong_lines:
        E_de = E_p - 1022.0
        if E_de > 100 and abs(peak_E_keV - E_de) <= tol:
            return ResidualClassification(
                peak_E_keV=peak_E_keV, sigma=sigma,
                label=LBL_DOUBLE_ESCAPE,
                note=f"double escape от {parent} {E_p:.0f} кэВ (E−1022={E_de:.0f})",
                parent_nuclide=parent, parent_line_keV=E_p,
                feature_kind="double_escape",
                delta_keV=abs(peak_E_keV - E_de),
            )

    # 4. Sum peak (E_a + E_b, both detected, both ≥ 5% intensity)
    if peak_E_keV > 500:
        for i, (na, Ea) in enumerate(parent_strong_lines):
            for nb, Eb in parent_strong_lines[i:]:
                E_sum = Ea + Eb
                if abs(peak_E_keV - E_sum) <= tol and E_sum <= energy_max_keV + 50:
                    same = na == nb
                    return ResidualClassification(
                        peak_E_keV=peak_E_keV, sigma=sigma,
                        label=LBL_SUM_PEAK,
                        note=(f"sum-peak {na} {Ea:.0f} + "
                              f"{'(self)' if same else nb} {Eb:.0f}"
                              f" = {E_sum:.0f} кэВ"),
                        parent_nuclide=na,
                        parent_line_keV=E_sum,
                        feature_kind="sum_peak",
                        delta_keV=abs(peak_E_keV - E_sum),
                    )

    # 4b. TD-3 sum-peak fallback (cascade-coincidence) at I ≥ 3%.
    # Chain co-emissions with intensity ≥ 3% can produce visible sum
    # peaks via cascade coincidence even when one partner falls below
    # the 5% threshold of the primary pass. Classic case: Tl-208 583 +
    # Ac-228 1588 ≈ 2171 keV in Th-232 chain. Relaxed only for the
    # 4b pass (sum peak), not for escapes.
    if peak_E_keV > 500:
        parent_medium_lines: List[Tuple[str, float]] = []
        for n in detected_nuclide_names:
            for E in _strong_lines_of(n, intensity_floor=3.0):
                parent_medium_lines.append((n, E))
        strong_set = {(n, E) for n, E in parent_strong_lines}
        for i, (na, Ea) in enumerate(parent_medium_lines):
            for nb, Eb in parent_medium_lines[i:]:
                # Skip pairs already exhausted by the 5% pass above.
                if ((na, Ea) in strong_set and (nb, Eb) in strong_set):
                    continue
                E_sum = Ea + Eb
                if abs(peak_E_keV - E_sum) <= tol and E_sum <= energy_max_keV + 50:
                    same = na == nb
                    return ResidualClassification(
                        peak_E_keV=peak_E_keV, sigma=sigma,
                        label=LBL_SUM_PEAK,
                        note=(f"sum-peak (cascade, I≥3%) {na} {Ea:.0f} + "
                              f"{'(self)' if same else nb} {Eb:.0f}"
                              f" = {E_sum:.0f} кэВ"),
                        parent_nuclide=na,
                        parent_line_keV=E_sum,
                        feature_kind="sum_peak",
                        delta_keV=abs(peak_E_keV - E_sum),
                    )

    # 5. Chain-secondary via v2 catalog
    for parent in detected_nuclide_names:
        hits = matches_secondary(parent, peak_E_keV, span="p10p90")
        non_pp = [h for h in hits if h.get("feature") != "photopeak"]
        if non_pp:
            h = non_pp[0]
            lo, hi = h.get("range", (peak_E_keV, peak_E_keV))
            return ResidualClassification(
                peak_E_keV=peak_E_keV, sigma=sigma,
                label=LBL_CHAIN_SECONDARY,
                note=(f"{h.get('feature','?')} of {parent} "
                      f"({h.get('primary_line_keV','?')} кэВ) "
                      f"empirical range [{lo:.0f}..{hi:.0f}]"),
                parent_nuclide=parent,
                parent_line_keV=h.get("primary_line_keV"),
                feature_kind=h.get("feature", ""),
            )

    # 6. XRF
    xrf_hits = _match_xrf(peak_E_keV, tolerance_keV=tol)
    if xrf_hits:
        # Prefer the closest hit
        el, shell, lib_E, dE = sorted(xrf_hits, key=lambda x: x[3])[0]
        return ResidualClassification(
            peak_E_keV=peak_E_keV, sigma=sigma,
            label=LBL_XRF,
            note=(f"{el} {shell}-XRF (lib {lib_E:.1f} кэВ, Δ={dE:.1f})"),
            element=el,
            parent_line_keV=lib_E,
            delta_keV=dE,
        )

    # 7. Near energy ceiling
    if peak_E_keV > energy_max_keV - fwhm_at_keV:
        return ResidualClassification(
            peak_E_keV=peak_E_keV, sigma=sigma,
            label=LBL_EDGE_OF_RANGE,
            note=(f"в пределах 1·FWHM от потолка ({energy_max_keV:.0f} кэВ) — "
                  "возможен truncation-артефакт"),
        )

    # F-143 / v1.17.7 — для сильных моно-нуклидов (Cs-137 / K-40) region
    # 100-400 keV содержит выраженный Compton-continuum от photopeak.
    # Library-window matches здесь часто ложные (U-235 185.7 / Ra-226 186 /
    # Bi-211 / etc.). Если binding к single-isotope известен, классифицируем
    # эти peaks как `chain_secondary` со feature_kind="compton_residual_of_<nuc>".
    for parent in detected_nuclide_names:
        if parent == "Cs-137" and 100.0 <= peak_E_keV <= 400.0:
            return ResidualClassification(
                peak_E_keV=peak_E_keV, sigma=sigma,
                label=LBL_CHAIN_SECONDARY,
                note=(f"Compton-continuum от Cs-137 662 кэВ "
                      f"(E={peak_E_keV:.1f} в диапазоне 100-400 кэВ — "
                      f"вторичный процесс, не U-235/Ra-226)."),
                parent_nuclide=parent,
                parent_line_keV=661.66,
                feature_kind="compton_residual",
            )
        if parent == "K-40" and 200.0 <= peak_E_keV <= 1300.0:
            return ResidualClassification(
                peak_E_keV=peak_E_keV, sigma=sigma,
                label=LBL_CHAIN_SECONDARY,
                note=(f"Compton-continuum от K-40 1461 кэВ "
                      f"(E={peak_E_keV:.1f} в диапазоне 200-1300 кэВ)."),
                parent_nuclide=parent,
                parent_line_keV=1460.82,
                feature_kind="compton_residual",
            )
        # F-329 / v1.18.18.3 (ROADMAP TD-2) — Co-60 имеет ДВА FEP
        # 1173.23 + 1332.50 keV; Compton-края соответственно ≈ 963 / 1118 кэВ.
        # Residual continuum 200-1200 кэВ часто даёт ложные библиотечные
        # совпадения (Bi-214 609 / Cs-137 662 — particularly insidious
        # на calibration-source спектрах).
        if parent == "Co-60" and 200.0 <= peak_E_keV <= 1200.0:
            return ResidualClassification(
                peak_E_keV=peak_E_keV, sigma=sigma,
                label=LBL_CHAIN_SECONDARY,
                note=(f"Compton-continuum от Co-60 1173/1332 кэВ "
                      f"(E={peak_E_keV:.1f} в диапазоне 200-1200 кэВ). "
                      "Compton-края 963/1118 кэВ дают плотный continuum."),
                parent_nuclide=parent,
                parent_line_keV=1173.23,
                feature_kind="compton_residual",
            )
        # F-329 — Na-22: 1274.5 keV FEP + 511 back-to-back annihilation.
        # Compton-edge от 1275 ≈ 1062 кэВ; континуум 150-400 кэВ
        # часто содержит back-to-back артефакты + residual continuum
        # от обеих линий.
        if parent == "Na-22" and 150.0 <= peak_E_keV <= 400.0:
            return ResidualClassification(
                peak_E_keV=peak_E_keV, sigma=sigma,
                label=LBL_CHAIN_SECONDARY,
                note=(f"Compton-continuum от Na-22 1274.5 кэВ + back-to-back "
                      f"511 кэВ (E={peak_E_keV:.1f} в диапазоне 150-400 кэВ)."),
                parent_nuclide=parent,
                parent_line_keV=1274.5,
                feature_kind="compton_residual",
            )

    # 8. Truly unmatched
    return ResidualClassification(
        peak_E_keV=peak_E_keV, sigma=sigma,
        label=LBL_TRUE_UNMATCHED,
        note="не классифицируется как вторичная ЕРН, рентген. флуор., ускользание, сумм. или край",
    )


def classify_residuals(
    unmatched_peaks: List,
    spec,
    detected_nuclides: List,
    fwhm_provider_keV,
    *,
    energy_max_keV: Optional[float] = None,
    sigma_floor: float = 2.5,
) -> List[ResidualClassification]:
    """
    Classify every unmatched peak from a staged-identification result.

    Parameters
    ----------
    unmatched_peaks : list of FoundPeak
        Peaks not claimed by any nuclide after disambiguate.
    spec : Spectrum
    detected_nuclides : list of NuclideIdentification
        Detected nuclide list (post-disambiguate).
    fwhm_provider_keV : Callable[[float], float]
        Function returning FWHM in keV at any energy.
    energy_max_keV : float, optional
        Defaults to `spec.energy_max_keV_kept`.
    sigma_floor : float
        Peaks below this significance are skipped (returned with label
        empty in the list element's place).

    Returns
    -------
    list[ResidualClassification]
        Same length as `unmatched_peaks` for the σ ≥ floor entries,
        skipping the rest.
    """
    if energy_max_keV is None:
        energy_max_keV = float(getattr(spec, "energy_max_keV_kept", 3000.0))
    detected_names = [n.nuclide for n in detected_nuclides]
    out: List[ResidualClassification] = []
    for p in unmatched_peaks:
        if p.significance < sigma_floor:
            continue
        e = spec.channel_to_energy(p.channel)
        fwhm = max(2.0, float(fwhm_provider_keV(e)))
        out.append(classify_residual(
            peak_E_keV=e, sigma=float(p.significance),
            detected_nuclide_names=detected_names,
            fwhm_at_keV=fwhm,
            energy_max_keV=energy_max_keV,
        ))
    return out


__all__ = [
    "classify_residual",
    "classify_residuals",
    "ResidualClassification",
    "LBL_XRF", "LBL_CHAIN_SECONDARY", "LBL_SINGLE_ESCAPE",
    "LBL_DOUBLE_ESCAPE", "LBL_SUM_PEAK", "LBL_ANNIHILATION",
    "LBL_EDGE_OF_RANGE", "LBL_TRUE_UNMATCHED",
]
