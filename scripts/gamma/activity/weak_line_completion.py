"""F-440 / v1.30.0 -- Two-phase weak-line completion (Phase 2).

Formula:
    S_expected(E_w) = A * I_gamma * eps * t_live * f_self_abs * f_TCS
    sigma_S = S_expected * sqrt((sigmaA/A)**2 + (sigmaI/I)**2 + (sigmaEps/Eps)**2)

Caller injects all inputs. Pure physics module.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple


DEFAULT_MIN_GROUPING_SNR = 5.0
DEFAULT_MIN_GROUPING_INTENSITY_PCT = 3.0


@dataclass(frozen=True)
class CompletedLine:
    E_keV: float
    I_pct: float
    sigma_I_pct_relative: float
    epsilon: float
    epsilon_extrapolated: bool
    f_self_abs: float
    f_TCS: float
    S_expected: float
    sigma_S_expected: float
    in_ROI_of_multiplet: Optional[str] = None


@dataclass(frozen=True)
class FittedLineSummary:
    E_keV: float
    I_pct: float
    S_measured: Optional[float]
    peak_area_source: str


@dataclass(frozen=True)
class WeakContamination:
    strong_nuclide: str
    from_nuclide: str
    from_lines_keV: Tuple[float, ...]
    S_contamination: float
    fraction_of_strong_pct: float
    multiplet_label: str = ""


@dataclass(frozen=True)
class NuclideCompletion:
    nuclide: str
    phase1_activity_Bq: Optional[float]
    phase1_activity_sigma_Bq: Optional[float]
    phase1_specific_activity_Bq_per_kg: Optional[float]
    phase1_specific_activity_sigma_Bq_per_kg: Optional[float]
    phase1_fitted_lines: Tuple[FittedLineSummary, ...]
    phase2_completed_lines: Tuple[CompletedLine, ...]
    fitted_area_total: float
    completed_area_total: float
    completeness_pct: float
    weak_contamination_into_strong_peaks: Tuple[WeakContamination, ...]



def _safe_float(x) -> Optional[float]:
    try:
        v = float(x)
        if math.isnan(v) or math.isinf(v):
            return None
        return v
    except (TypeError, ValueError):
        return None


def _line_records(lib_entry) -> List[Tuple[float, float, float]]:
    out: List[Tuple[float, float, float]] = []
    if not lib_entry:
        return out
    lines = lib_entry.get("lines") or []
    for line in lines:
        try:
            if isinstance(line, (list, tuple)) and len(line) >= 2:
                E = float(line[0])
                I_val = float(line[1])
                sigma_I = float(line[2]) if len(line) >= 3 else 0.0
            elif isinstance(line, dict):
                E = float(line.get("E_keV") or 0.0)
                I_val = float(line.get("I_pct") or 0.0)
                sigma_I = float(line.get("sigma_I_pct") or 0.0)
            else:
                continue
            if E <= 0 or I_val < 0:
                continue
            out.append((E, I_val, sigma_I))
        except (TypeError, ValueError):
            continue
    return out


def _lookup_factor(
    factors: Optional[Dict[float, float]], E: float, tol_keV: float = 0.5
) -> float:
    if not factors:
        return 1.0
    for key_E, key_f in factors.items():
        try:
            if abs(float(key_E) - float(E)) <= tol_keV:
                return float(key_f)
        except (TypeError, ValueError):
            continue
    return 1.0


def _is_weak_match(line_match, min_snr: float, min_intensity_pct: float) -> bool:
    """True iff line is below Phase 1 strong-criteria.

    F-440 rule: a line is STRONG (Phase 1) when it has a real measured area
    (peak_area > 0 AND source != 'library_anchor_phantom') AND I_gamma passes
    the intensity gate. Lines below intensity gate, or with phantom/missing
    area, are WEAK (Phase 2).
    """
    I_pct = float(getattr(line_match, "library_I_pct", 0.0) or 0.0)
    src = str(getattr(line_match, "peak_area_source", "") or "")
    if "phantom" in src.lower():
        return True
    if min_intensity_pct > 0.0 and I_pct < min_intensity_pct:
        return True
    area = getattr(line_match, "peak_area", None)
    area_unc = getattr(line_match, "peak_area_uncertainty", None)
    has_area = False
    try:
        if area is not None and float(area) > 0:
            has_area = True
    except (TypeError, ValueError):
        has_area = False
    if not has_area:
        sig = getattr(line_match, "significance_currie", None)
        if sig is not None:
            try:
                return float(sig) < min_snr
            except (TypeError, ValueError):
                pass
        return True
    if area_unc in (None, 0, 0.0):
        # Measured area present, no uncertainty: trust as strong (passed
        # identification on real data; common for cowell / library_anchor_strong).
        return False
    try:
        snr = float(area) / float(area_unc)
        return snr < min_snr
    except (TypeError, ValueError, ZeroDivisionError):
        return False



def complete_weak_lines(
    activities,
    matched_lines_by_nuclide: Dict[str, List[Any]],
    nuclide_library: Dict[str, dict],
    efficiency_curve,
    t_live: float,
    *,
    sample_mass_kg: Optional[float] = None,
    specific_activities_Bq_per_kg: Optional[Dict[str, Tuple[float, float]]] = None,
    self_attenuation_factors: Optional[Dict[float, float]] = None,
    tcs_factors: Optional[Dict[Tuple[str, float], float]] = None,
    min_grouping_snr: float = DEFAULT_MIN_GROUPING_SNR,
    min_grouping_intensity_pct: float = DEFAULT_MIN_GROUPING_INTENSITY_PCT,
    multiplet_rois: Optional[List[dict]] = None,
    contamination_threshold_pct: float = 1.0,
    min_library_intensity_for_completeness: float = 0.1,
) -> Dict[str, NuclideCompletion]:
    if activities is None:
        return {}
    activities_by_nuc: Dict[str, Any] = {}
    for ar in activities:
        nuc = getattr(ar, "nuclide", None)
        if nuc:
            activities_by_nuc[nuc] = ar
    spec_acts = specific_activities_Bq_per_kg or {}
    nuclide_strong_S: Dict[str, Dict[float, float]] = {}
    nuclide_completed: Dict[str, List[CompletedLine]] = {}
    nuclide_fitted: Dict[str, List[FittedLineSummary]] = {}

    for nuc, matches in matched_lines_by_nuclide.items():
        ar = activities_by_nuc.get(nuc)
        A_Bq = _safe_float(getattr(ar, "A_Bq", None)) if ar else None
        sigma_A_Bq = _safe_float(getattr(ar, "sigma_A_Bq", None)) if ar else None
        rel_sigma_A = None
        if A_Bq is not None and A_Bq > 0 and sigma_A_Bq is not None:
            rel_sigma_A = sigma_A_Bq / A_Bq
        strong_S_by_E: Dict[float, float] = {}
        fitted_summaries: List[FittedLineSummary] = []
        for m in matches:
            E_lib = _safe_float(getattr(m, "library_E_keV", None))
            if E_lib is None:
                continue
            I_pct = float(getattr(m, "library_I_pct", 0.0) or 0.0)
            area = getattr(m, "peak_area", None)
            src = str(getattr(m, "peak_area_source", "") or "")
            if not _is_weak_match(m, min_grouping_snr, min_grouping_intensity_pct):
                fitted_summaries.append(FittedLineSummary(
                    E_keV=E_lib, I_pct=I_pct,
                    S_measured=_safe_float(area),
                    peak_area_source=src,
                ))
                if area is not None:
                    try:
                        strong_S_by_E[round(E_lib, 2)] = float(area)
                    except (TypeError, ValueError):
                        pass
        nuclide_strong_S[nuc] = strong_S_by_E
        nuclide_fitted[nuc] = fitted_summaries

        if A_Bq is None or A_Bq <= 0 or efficiency_curve is None or t_live <= 0:
            nuclide_completed[nuc] = []
            continue
        lib_entry = nuclide_library.get(nuc) if nuclide_library else None
        lib_lines = _line_records(lib_entry)
        if not lib_lines:
            nuclide_completed[nuc] = []
            continue
        strong_E_set = {round(s.E_keV, 2) for s in fitted_summaries}
        completed: List[CompletedLine] = []
        for E_lib, I_pct_lib, sigma_I_lib_abs in lib_lines:
            if round(E_lib, 2) in strong_E_set:
                continue
            # F-440: skip ultra-weak library lines (noise tier) so completeness
            # reflects the physically observable portion only.
            if (min_library_intensity_for_completeness > 0.0
                    and I_pct_lib < min_library_intensity_for_completeness):
                continue
            try:
                eps = float(efficiency_curve.efficiency_at(E_lib))
            except Exception:
                continue
            if eps is None or eps <= 0:
                continue
            try:
                eps_extrap = bool(efficiency_curve.is_extrapolating(E_lib))
            except Exception:
                eps_extrap = False
            f_self = _lookup_factor(self_attenuation_factors, E_lib)
            f_tcs = 1.0
            if tcs_factors:
                key = (nuc, round(E_lib, 2))
                if key in tcs_factors:
                    try:
                        f_tcs = float(tcs_factors[key])
                    except (TypeError, ValueError):
                        f_tcs = 1.0
            S_exp = A_Bq * (I_pct_lib / 100.0) * eps * t_live * f_self * f_tcs
            sigma_A_rel = rel_sigma_A or 0.0
            sigma_I_rel = (sigma_I_lib_abs / I_pct_lib) if I_pct_lib > 0 else 0.0
            sigma_eps_rel = 0.10 if eps_extrap else 0.03
            rel_sigma_S = math.sqrt(sigma_A_rel ** 2 + sigma_I_rel ** 2 + sigma_eps_rel ** 2)
            sigma_S = S_exp * rel_sigma_S
            in_roi: Optional[str] = None
            if multiplet_rois:
                for roi in multiplet_rois:
                    try:
                        E_lo = float(roi.get("E_lo_keV") or 0)
                        E_hi = float(roi.get("E_hi_keV") or 0)
                        if E_lo <= E_lib <= E_hi:
                            in_roi = str(roi.get("label") or "")
                            break
                    except (TypeError, ValueError):
                        continue
            completed.append(CompletedLine(
                E_keV=E_lib, I_pct=I_pct_lib,
                sigma_I_pct_relative=sigma_I_rel,
                epsilon=eps, epsilon_extrapolated=eps_extrap,
                f_self_abs=f_self, f_TCS=f_tcs,
                S_expected=S_exp, sigma_S_expected=sigma_S,
                in_ROI_of_multiplet=in_roi,
            ))
        nuclide_completed[nuc] = completed

    roi_strong_per_nuc: Dict[str, Dict[str, float]] = {}
    if multiplet_rois:
        for roi in multiplet_rois:
            label = str(roi.get("label") or "")
            E_lo = float(roi.get("E_lo_keV") or 0)
            E_hi = float(roi.get("E_hi_keV") or 0)
            roi_strong_per_nuc.setdefault(label, {})
            for nuc, strong_map in nuclide_strong_S.items():
                S_in_roi = 0.0
                for E_key, S_val in strong_map.items():
                    if E_lo <= E_key <= E_hi:
                        S_in_roi += S_val
                if S_in_roi > 0:
                    roi_strong_per_nuc[label][nuc] = S_in_roi

    out: Dict[str, NuclideCompletion] = {}
    for nuc in matched_lines_by_nuclide.keys():
        ar = activities_by_nuc.get(nuc)
        A_Bq = _safe_float(getattr(ar, "A_Bq", None)) if ar else None
        sigma_A_Bq = _safe_float(getattr(ar, "sigma_A_Bq", None)) if ar else None
        spec_pair = spec_acts.get(nuc) if isinstance(spec_acts, dict) else None
        if spec_pair is not None:
            try:
                A_Bq_kg = _safe_float(spec_pair[0])
                sigma_A_Bq_kg = _safe_float(spec_pair[1])
            except (TypeError, IndexError):
                A_Bq_kg = None
                sigma_A_Bq_kg = None
        else:
            A_Bq_kg = None
            sigma_A_Bq_kg = None
            if A_Bq is not None and sample_mass_kg and sample_mass_kg > 0:
                A_Bq_kg = A_Bq / sample_mass_kg
                if sigma_A_Bq is not None:
                    sigma_A_Bq_kg = sigma_A_Bq / sample_mass_kg

        fitted = tuple(nuclide_fitted.get(nuc, []))
        completed = tuple(nuclide_completed.get(nuc, []))
        S_fitted_total = sum((f.S_measured or 0.0) for f in fitted if f.S_measured is not None)
        S_completed_total = sum(c.S_expected for c in completed)
        denom = S_fitted_total + S_completed_total
        completeness_pct = 100.0 * S_fitted_total / denom if denom > 0 else 100.0

        contam_records: List[WeakContamination] = []
        for roi_label, strong_map in roi_strong_per_nuc.items():
            completed_in_roi = [c for c in completed if c.in_ROI_of_multiplet == roi_label]
            if not completed_in_roi:
                continue
            S_contam_total = sum(c.S_expected for c in completed_in_roi)
            for other_nuc, S_strong in strong_map.items():
                if other_nuc == nuc:
                    continue
                if S_strong <= 0:
                    continue
                frac_pct = 100.0 * S_contam_total / S_strong
                if frac_pct >= contamination_threshold_pct:
                    contam_records.append(WeakContamination(
                        strong_nuclide=other_nuc,
                        from_nuclide=nuc,
                        from_lines_keV=tuple(round(c.E_keV, 2) for c in completed_in_roi),
                        S_contamination=S_contam_total,
                        fraction_of_strong_pct=frac_pct,
                        multiplet_label=roi_label,
                    ))

        out[nuc] = NuclideCompletion(
            nuclide=nuc,
            phase1_activity_Bq=A_Bq,
            phase1_activity_sigma_Bq=sigma_A_Bq,
            phase1_specific_activity_Bq_per_kg=A_Bq_kg,
            phase1_specific_activity_sigma_Bq_per_kg=sigma_A_Bq_kg,
            phase1_fitted_lines=fitted,
            phase2_completed_lines=completed,
            fitted_area_total=S_fitted_total,
            completed_area_total=S_completed_total,
            completeness_pct=completeness_pct,
            weak_contamination_into_strong_peaks=tuple(contam_records),
        )

    return out



def to_json_block(completions: Dict[str, NuclideCompletion]) -> Dict[str, dict]:
    out: Dict[str, dict] = {}
    for nuc, c in completions.items():
        out[nuc] = {
            "phase1_fitted_lines": [
                {
                    "E": round(f.E_keV, 3),
                    "I_pct": round(f.I_pct, 4),
                    "S_measured": (round(f.S_measured, 3) if f.S_measured is not None else None),
                    "peak_area_source": f.peak_area_source,
                }
                for f in c.phase1_fitted_lines
            ],
            "phase1_activity_Bq": (round(c.phase1_activity_Bq, 4) if c.phase1_activity_Bq is not None else None),
            "phase1_activity_sigma_Bq": (round(c.phase1_activity_sigma_Bq, 4) if c.phase1_activity_sigma_Bq is not None else None),
            "phase1_activity_Bq_kg": (round(c.phase1_specific_activity_Bq_per_kg, 4) if c.phase1_specific_activity_Bq_per_kg is not None else None),
            "phase1_activity_sigma_Bq_kg": (round(c.phase1_specific_activity_sigma_Bq_per_kg, 4) if c.phase1_specific_activity_sigma_Bq_per_kg is not None else None),
            "phase2_completed_lines": [
                {
                    "E": round(cl.E_keV, 3),
                    "I_pct": round(cl.I_pct, 4),
                    "epsilon": cl.epsilon,
                    "epsilon_extrapolated": cl.epsilon_extrapolated,
                    "f_self_abs": round(cl.f_self_abs, 4),
                    "f_TCS": round(cl.f_TCS, 4),
                    "S_expected": round(cl.S_expected, 3),
                    "sigma": round(cl.sigma_S_expected, 3),
                    "in_ROI_of_M": cl.in_ROI_of_multiplet,
                }
                for cl in c.phase2_completed_lines
            ],
            "fitted_area_total": round(c.fitted_area_total, 3),
            "completed_area_total": round(c.completed_area_total, 3),
            "completeness_pct": round(c.completeness_pct, 3),
            "weak_contamination_into_strong_peaks": [
                {
                    "strong_nuclide": w.strong_nuclide,
                    "from": w.from_nuclide,
                    "from_lines_keV": list(w.from_lines_keV),
                    "S_contamination": round(w.S_contamination, 3),
                    "fraction_of_strong_pct": round(w.fraction_of_strong_pct, 4),
                    "multiplet_label": w.multiplet_label,
                }
                for w in c.weak_contamination_into_strong_peaks
            ],
        }
    return out


__all__ = [
    "CompletedLine", "FittedLineSummary", "WeakContamination",
    "NuclideCompletion",
    "complete_weak_lines", "to_json_block",
    "DEFAULT_MIN_GROUPING_SNR", "DEFAULT_MIN_GROUPING_INTENSITY_PCT",
]
