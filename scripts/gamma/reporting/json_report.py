"""
JSON report assembler — primary output of Step 11.

Per `references/06_report_format.md`, the JSON report is the machine-
readable artefact written to disk; the chat layer only shows a 3–8
line summary. This module accepts a `StagedAnalysisResult` (from
`gamma.identification.staged_pipeline.analyze_lsrm_spe`) and returns a
dict ready for `json.dump`.

Schema version: **0.8** (see SCHEMA_VERSION constant; was 0.1 at Phase 3.1
inception, bumped through 0.2/0.3/0.4 during F-QC and detector-folder waves,
0.5 as of v1.21.0+ F-QC-01 unified spectrum_qc, 0.6 as of v1.31.0 G4
bg_carryover annotation on sample-side peaks, 0.7 as of v1.31.1 G5
calibration.energy_cal.source label / reused flag from spec.energy_cal_source,
0.8 as of v1.31.2 G2/G3 — intensity-ratio chi² gate top-level block and
continuum_diagnostic block inside header).
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional

from gamma.reporting.environment import classify_environment, continuum_diagnostic
# F-QC-01 / v1.21.0 — unified 6-criterion spectrum_qc aggregator (RAG-041)
from gamma.reporting.spectrum_qc_aggregator import build_spectrum_qc as _build_spectrum_qc
# F-452 / v1.33.0 — FwhmModel polymorphic API (callable + legacy 3-tuple shim)
from gamma.identification.staged_pipeline import (
    fwhm_keV_at_energy as _fwhm_keV_at_energy,
    fwhm_model_legacy_abc as _fwhm_model_legacy_abc,
)
# G4 / v1.31.0 — annotate sample-side peaks that mirror a bg-side peak.
# Pure annotation pass, reuses bg_staged_result.final_detected.
from gamma.identification.bg_carryover import (
    build_bg_peak_catalog as _build_bg_peak_catalog,
    mark_bg_carryover as _mark_bg_carryover,
)
# F-RPT-08 / v1.19.1 — ISO 11929 §5.4.3/§5.4.4 decision threshold + detection limit
from gamma.math.iso_11929_thresholds import (
    decision_threshold as _iso_decision_threshold,
    detection_limit as _iso_detection_limit,
)


SCHEMA_VERSION = "0.8"
SKILL_VERSION = "v1.32.0"


# v1.2.18 — operator-locked render «приклеить сумму к спектру»: когда
# Phase A free-centroid side-fit нашла ненулевые сдвиги центроидов,
# multiplet_continua overlay перерисовывается на fitted positions
# (E_lib + shift_k). Площади и σ остаются из основного locked-fit
# (shift'ы малы → area integral Gaussian'a инвариант с точностью до o(shift²)).
# v1.30.3 / F-445 — `_compute_cluster_global_shift` PROMOTED из render-слоя в
# calibration: now lives as public `compute_cluster_global_shift` in
# `gamma.calibration.cluster_shift_anchors`. F-445 Phase B/C использует
# его same function для построения cluster-level anchor'ов перед E-cal refit.
# Render-time override остаётся как fallback (срабатывает, если F-445
# Phase C откатился и данные так и не сдвинулись).
from gamma.calibration.cluster_shift_anchors import (
    compute_cluster_global_shift as _compute_cluster_global_shift,
)


def _rebuild_overlay_on_fitted_centroids(
    E_arr_raw, cont_arr_raw, shifts_kev, rebuild_specs,
    comp_list, fallback_total,
):
    """Return (total_arr_rounded, comp_list_with_updated_g_curve)."""
    import math as _m
    has_shift = any(abs(float(s)) > 1e-9 for s in (shifts_kev or ()))
    if not (has_shift and rebuild_specs and E_arr_raw):
        return ([round(float(v), 3) for v in (fallback_total or ())], comp_list)
    total_v2 = list(cont_arr_raw)
    for (E_fit_c, amp_c, sigma_c) in rebuild_specs:
        if sigma_c <= 0:
            continue
        inv = 1.0 / (2.0 * sigma_c * sigma_c)
        for i, Ex in enumerate(E_arr_raw):
            dE = Ex - E_fit_c
            total_v2[i] += amp_c * _m.exp(-dE * dE * inv)
    for ci, (E_fit_c, amp_c, sigma_c) in enumerate(rebuild_specs):
        if ci >= len(comp_list) or sigma_c <= 0:
            continue
        inv = 1.0 / (2.0 * sigma_c * sigma_c)
        comp_list[ci]["g_curve"] = [
            round(amp_c * _m.exp(-(Ex - E_fit_c) ** 2 * inv), 3)
            for Ex in E_arr_raw
        ]
    return ([round(v, 3) for v in total_v2], comp_list)




def _safe_float(x) -> Optional[float]:
    """JSON-friendly float: None for NaN / inf."""
    if x is None:
        return None
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    if math.isnan(v) or math.isinf(v):
        return None
    return v


def _filename_leaf(path: str) -> str:
    if not path:
        return ""
    sep_pos = max(path.rfind("/"), path.rfind("\\"))
    return path[sep_pos + 1:] if sep_pos >= 0 else path


def _build_header(result) -> Dict[str, Any]:
    """Section 1 — file/measurement header."""
    spec = result.spec
    sp = getattr(spec, "source_path", "") or ""
    live = float(getattr(spec, "live_time", 0.0) or 0.0)
    real = float(getattr(spec, "real_time", 0.0) or 0.0)
    dead_pct: Optional[float]
    if real > 0 and live > 0:
        dead_pct = max(0.0, 100.0 * (1.0 - live / real))
    else:
        dead_pct = None

    start_dt = getattr(spec, "start_datetime", None)
    start_iso = start_dt.isoformat() if start_dt is not None else None

    bg_sub = result.background_subtraction
    bg_present = bool(result.is_background) or bg_sub is not None

    env = classify_environment(result)

    # F-89a / v1.15.2 — background status MUST always be surfaced.
    bg_status = getattr(result, "background_status", "") or (
        "subtracted_from_external_file" if bg_sub is not None
        else ("embedded_present_not_subtracted"
              if getattr(spec, "background_embedded", None) is not None
              else "absent_no_subtraction")
    )

    # F-144 / v1.17.7 — имя фон-файла должно быть в отчёте.
    # Источник приоритета:
    #   1. auto_background_applied_path (F-131/F-135 apply mode)
    #   2. background_subtraction.bg_source_path (явный --background-path)
    #   3. embedded_path / link
    bg_filename = ""
    bg_path = ""
    auto_bg_path = getattr(result, "auto_background_applied_path", None)
    if auto_bg_path:
        bg_path = str(auto_bg_path)
        bg_filename = _filename_leaf(bg_path)
    elif bg_sub is not None:
        for attr in ("bg_source_path", "source_path", "bg_path"):
            v = getattr(bg_sub, attr, None)
            if v:
                bg_path = str(v); bg_filename = _filename_leaf(bg_path)
                break

    return {
        "filename": _filename_leaf(sp),
        "sample_filename": _filename_leaf(sp),
        "source_path": sp,
        # F-144
        "background_filename": bg_filename,
        "background_path": bg_path,
        "sample_id": getattr(spec, "sample_id", "") or "",
        "operator": getattr(spec, "operator", "") or "",
        "start_datetime": start_iso,
        "live_time_s": _safe_float(live),
        "real_time_s": _safe_float(real),
        "dead_time_pct": _safe_float(dead_pct),
        "geometry": getattr(spec, "geometry", "") or "",
        "geometry_canonical": result.geometry_canonical or "",
        "detector_type": result.detector_type,
        "detector_id": getattr(spec, "detector_id", "") or "",
        "detector_canonical": result.detector_canonical or "",
        # BUG-39 / BUG-40 — silent-fallback record for the detector profile.
        # Schema: {"requested": str, "actual": str, "reason": str, "human": str,
        #          "cyrillic_to_latin_collision": bool}.
        # `reason == "profile_loaded_no_fallback"` means no substitution happened.
        # F-115 (2026-06-21): strip `original_raw` here — it carries the
        # verbatim header (incl. serial №NNNN-NN) and would leak under the
        # anonymization regression. The homoglyph signal survives via the
        # `cyrillic_to_latin_collision` boolean; raw header stays in stderr
        # logs for operator debugging but never reaches the JSON artifact.
        "detector_fallback": {
            k: v for k, v in (getattr(result, "detector_fallback", None) or {}).items()
            if k != "original_raw"
        },
        "device_guid": getattr(spec, "device_guid", "") or "",
        "background_present": bg_present,
        "background_subtracted": bg_sub is not None,
        "background_status": bg_status,
        "environment": env,
        # G3 / v1.31.2 — continuum-level diagnostic (annotation, no decision).
        # Three integral metrics surfaced next to the categorical environment
        # label: total_cps, bg_line_dominance_pct, environment_hint_by_cps.
        # Divergence between env (Pb-XRF/ERN logic) and hint_by_cps is the
        # signal the operator wants to see.
        "continuum_diagnostic": continuum_diagnostic(result),
        "sample_type_hint": result.sample_type_hint or "",
        "sample_type_canonical": result.sample_type_canonical or "",
        # F-89b — surface filename binding hypothesis at the top.
        "filename_isotope_hints": list(getattr(result, "filename_isotope_hints", []) or []),
        "filename_chains_claimed": list(getattr(result, "filename_chains_claimed", []) or []),
        "analysis_mode": result.analysis_mode or "",
        "energy_ceiling_keV": _safe_float(getattr(spec, "ENERGY_CEILING_KEV", 3000.0)),
        "energy_max_keV_kept": _safe_float(getattr(spec, "energy_max_keV_kept", None)),
        "n_channels": int(getattr(spec, "n_channels", 0) or 0),
        "n_channels_raw": int(getattr(spec, "n_channels_raw", 0) or 0),
        "dropped_high_energy_count": int(
            getattr(spec, "dropped_overflow_count", 0) or 0
        ),
    }


_ENERGY_CAL_SOURCE_LABELS = {
    "stored": "использована сохранённая E-cal",
    "F-145_multiplet_self_calibration": "перестроена (F-145 multiplet self-cal)",
    "bootstrap": "перестроена (bootstrap)",
    "manual": "перестроена (ручная)",
    "": "источник не указан",
}


def _energy_cal_source_label(source: str) -> str:
    return _ENERGY_CAL_SOURCE_LABELS.get(
        source, f"перестроена ({source})" if source else "источник не указан"
    )


def _build_calibration(result) -> Dict[str, Any]:
    """Section 2 — calibration."""
    spec = result.spec
    e_cal = getattr(spec, "energy_cal", ()) or ()
    # F-452: FwhmModel (callable) → legacy 3-tuple для schema поля
    # `fwhm_cal.coefficients` (FWHM²(E)=a+b·E+c·E²). Для lsrm_poly_sqrt_E
    # ветки это least-squares refit (±5-7 keV documented); честный poly-4
    # остаётся в callable `result.fwhm_model(E)` и используется всеми
    # downstream-консьюмерами через `fwhm_keV_at_energy`.
    a, b, c = _fwhm_model_legacy_abc(result.fwhm_model)
    raw_source = str(getattr(spec, "energy_cal_source", "") or "")
    return {
        "energy_cal": {
            "degree": max(0, len(e_cal) - 1),
            "coefficients": [_safe_float(x) for x in e_cal],
            "source": raw_source or "unknown",
            "reused": raw_source == "stored",
            "source_label": _energy_cal_source_label(raw_source),
        },
        "fwhm_cal": {
            "model": "FWHM^2 = a + b*E + c*E^2",
            "coefficients": [_safe_float(a), _safe_float(b), _safe_float(c)],
            "source": result.fwhm_model_source,
            "fwhm_at_661_keV": _safe_float(result.fwhm_at_661),
        },
        "seven_line_check": _serialize_seven_line(result.seven_line_check),
    }


def _serialize_seven_line(slc) -> Optional[Dict[str, Any]]:
    if slc is None:
        return None
    per_line = []
    for r in getattr(slc, "per_line", []) or []:
        per_line.append({
            "expected_keV": _safe_float(getattr(r, "expected_keV", None)),
            "description": getattr(r, "description", ""),
            "line_type": getattr(r, "line_type", ""),
            "found": bool(getattr(r, "found", False)),
            "observed_keV": _safe_float(getattr(r, "observed_keV", None)),
            "residual_keV": _safe_float(getattr(r, "residual_keV", None)),
            "fwhm_at_keV": _safe_float(getattr(r, "fwhm_at_keV", None)),
            "residual_fwhm_fraction": _safe_float(getattr(r, "residual_fwhm_fraction", None)),
        })
    return {
        "lines_present": int(getattr(slc, "lines_present", 0) or 0),
        "lines_total": int(getattr(slc, "lines_total", 7) or 7),
        "max_residual_keV": _safe_float(getattr(slc, "max_residual_keV", None)),
        "mean_residual_keV": _safe_float(getattr(slc, "mean_residual_keV", None)),
        "max_residual_fwhm_fraction": _safe_float(getattr(slc, "max_residual_fwhm_fraction", None)),
        "quality": getattr(slc, "quality", "n/a"),
        "quality_note": getattr(slc, "quality_note", ""),
        "per_line": per_line,
    }


def _build_primary_feps(result) -> List[Dict[str, Any]]:
    """Section 4 — primary FEPs (one row per matched line of detected nuclides).

    F-391 / v1.18.27 — phantom anchors (peak_area_source ∈ {library_anchor,
    library_anchor_phantom}) НЕ попадают в primary_feps. Они представляют
    library-evidence для identified нуклидов без реального измеренного
    signal — их место в multiplet_deconvolutions.components (phantom flag)
    для evidence, но НЕ в primary FEP table. Это устраняет дубли
    «Identified nuclides» + peak list (П.9, П.10 v1.18.27 user feedback).
    """
    spec = result.spec
    live = float(getattr(spec, "live_time", 0.0) or 1.0)
    PHANTOM_SOURCES = {"library_anchor", "library_anchor_phantom"}
    out = []
    for ni in result.final_detected:
        nuc = ni.nuclide
        for m in ni.matched_lines:
            # F-391 — drop phantom anchors из primary_feps
            if str(m.peak_area_source or "") in PHANTOM_SOURCES:
                continue
            area = _safe_float(m.peak_area)
            area_unc = _safe_float(m.peak_area_uncertainty)
            rate_cps: Optional[float] = None
            rate_sigma_cps: Optional[float] = None
            if area is not None and live > 0:
                rate_cps = area / live
                if area_unc is not None:
                    rate_sigma_cps = area_unc / live
            # BUG-34 Phase 3b R2: prefer Semantic-B gauss_sigma_keV (W3
            # writer), fall back to legacy peak_sigma for callers/fixtures
            # that have not migrated yet. Numerically identical when only
            # peak_sigma is set. See PLAN_v1_18_32_to_v1_19_0 §3 Phase 2.
            _sigma = getattr(m, "gauss_sigma_keV", None)
            if _sigma is None:
                _sigma = m.peak_sigma if m.peak_sigma else None
            # F3 (Censor 2026-06-21): NLLS pinned-to-zero artifact detector.
            # При S < max(1.0, |σ_S|) пик статистически неотличим от нуля
            # (ISO 11929: значимая площадь должна превышать σ). Сохраняем
            # численные значения для аудита, но помечаем флагом — рендереры
            # покажут upper-limit нотацию вместо unphysical 1e-20.
            is_upper_limit_artifact = False
            if area is not None and area_unc is not None:
                threshold = max(1.0, abs(area_unc))
                if abs(area) < threshold:
                    is_upper_limit_artifact = True
            out.append({
                "nuclide": nuc,
                "library_E_keV": _safe_float(m.library_E_keV),
                "library_I_pct": _safe_float(m.library_I_pct),
                "peak_channel": int(m.peak_channel),
                "peak_E_keV": _safe_float(m.peak_E_keV),
                "fwhm_keV": _safe_float(_sigma * 2.355) if _sigma is not None else None,
                "residual_keV": _safe_float(m.residual_keV),
                "peak_area_counts": area,
                "peak_area_sigma": area_unc,
                "rate_cps": _safe_float(rate_cps),
                "rate_sigma_cps": _safe_float(rate_sigma_cps),
                "peak_area_source": m.peak_area_source,
                "is_characteristic": m.is_characteristic,
                "is_upper_limit_artifact": is_upper_limit_artifact,
            })
    return out


def _build_secondary_peaks(result) -> List[Dict[str, Any]]:
    """Section 5 — secondary peak table (residuals NOT classified as xrf/true)."""
    out = []
    for rc in result.residual_classifications:
        if rc.label in ("xrf", "true_unmatched", "edge_of_range"):
            # XRF goes to its own table; true_unmatched goes to "unidentified".
            continue
        out.append({
            "channel": None,
            "energy_keV": _safe_float(rc.peak_E_keV),
            "significance": _safe_float(rc.sigma),
            "type": rc.label,
            "feature_kind": rc.feature_kind,
            "parent_nuclide": rc.parent_nuclide,
            "parent_line_keV": _safe_float(rc.parent_line_keV),
            "note": rc.note,
        })
    return out


def _build_elemental_xrf(result) -> List[Dict[str, Any]]:
    """Section 6 — elemental XRF table.

    For Gamma-1S the typical entries are:
    * Pb K-XRF (shield fluorescence) — element="Pb"
    * Ba K-XRF from ¹³⁷ᵐBa IC daughter of Cs-137 — element="Ba"
    * Th / U K-XRF for natural-radon-chain matrices
    """
    by_element: Dict[str, Dict[str, Any]] = {}
    for rc in result.residual_classifications:
        if rc.label != "xrf":
            continue
        el = rc.element or "?"
        entry = by_element.setdefault(el, {
            "element": el,
            "observed_lines": [],
            "n_observed": 0,
            "mechanism": "",
            "note": rc.note,
        })
        entry["observed_lines"].append({
            "energy_keV": _safe_float(rc.peak_E_keV),
            "significance": _safe_float(rc.sigma),
            "library_E_keV": _safe_float(rc.parent_line_keV),
            "delta_keV": _safe_float(rc.delta_keV),
        })
        entry["n_observed"] += 1
    # Assign mechanism heuristically per Gamma-1S catalogue.
    # Pb → fluorescence_shield; Ba (with Cs-137 detected) → IC daughter
    # of Cs-137; Th / U → matrix fluorescence; W → shield variant.
    detected_names = {n.nuclide for n in result.final_detected}
    for el, entry in by_element.items():
        if el == "Pb":
            entry["mechanism"] = "fluorescence_shield"
        elif el == "Ba" and "Cs-137" in detected_names:
            entry["mechanism"] = "IC_Cs-137"
        elif el in {"Th", "U"}:
            entry["mechanism"] = "fluorescence_matrix"
        elif el == "W":
            entry["mechanism"] = "fluorescence_shield_collimator"
        else:
            entry["mechanism"] = "unknown"
    return list(by_element.values())


def _build_identified_nuclides(result) -> List[Dict[str, Any]]:
    """Section 7 — identified nuclides with CI, activity, Bq/kg."""
    spec = result.spec
    live = float(getattr(spec, "live_time", 0.0) or 1.0)

    activities_by_nuclide = {}
    if result.activities:
        for ar in result.activities:
            activities_by_nuclide[ar.nuclide] = ar
    specific_by_nuclide = result.specific_activities_Bq_per_kg or {}

    # CI gating dict (confirmed/tentative/noise)
    ci_tier = {}
    if result.ci_gating is not None:
        for tier in ("confirmed", "tentative", "noise"):
            for n in getattr(result.ci_gating, tier, []) or []:
                ci_tier[getattr(n, "nuclide", "")] = tier

    # F-93: 50% σ_A/A → upper-limit gate per LSRM §11
    # Pre-compute an L_D per nuclide using the smallest MDA_Bq among the
    # nuclide's lines (best detection threshold available).
    ld_by_nuclide: Dict[str, float] = {}
    if result.mda_per_line is not None:
        for (nuclide, _E), mda in result.mda_per_line.items():
            try:
                if mda.MDA_Bq is None or _safe_float(mda.MDA_Bq) is None:
                    continue
                cur = ld_by_nuclide.get(nuclide)
                if cur is None or mda.MDA_Bq < cur:
                    ld_by_nuclide[nuclide] = mda.MDA_Bq
            except Exception:
                continue

    UPPER_LIMIT_SIGMA_THRESHOLD = 0.50  # LSRM §11

    # F-RPT-08 / v1.19.1 — Pre-build per-(nuclide, E_keV) lookup for ISO 11929
    # thresholds.  Key: (nuclide, round(E_keV, 2)).  Value: MdaResult.
    # Populated from result.mda_per_line (same source as ld_by_nuclide above).
    mda_by_nuclide_energy: Dict = {}
    if result.mda_per_line is not None:
        for (nuclide, e_key), mda in result.mda_per_line.items():
            mda_by_nuclide_energy[(nuclide, round(float(e_key or 0), 2))] = mda

    mass_kg = getattr(result, "sample_mass_kg", None)

    out = []
    for ni in result.final_detected:
        ar = activities_by_nuclide.get(ni.nuclide)
        spec_act = specific_by_nuclide.get(ni.nuclide)
        # Sum peak rate over matched lines (cps, gross)
        peak_rate = 0.0
        for m in ni.matched_lines:
            if m.peak_area is not None and live > 0:
                peak_rate += m.peak_area / live

        # F-93: decide upper-limit reporting
        is_upper_limit = False
        upper_limit_Bq = None
        sigma_rel = None
        if ar is not None and ar.A_Bq and ar.sigma_A_Bq is not None and ar.A_Bq > 0:
            sigma_rel = ar.sigma_A_Bq / ar.A_Bq
            if sigma_rel > UPPER_LIMIT_SIGMA_THRESHOLD:
                is_upper_limit = True
                upper_limit_Bq = ld_by_nuclide.get(ni.nuclide)

        # F-RPT-08 / v1.19.1 — ISO 11929 §5.4.3/§5.4.4 per-nuclide thresholds
        # using the characteristic (highest-intensity) matched line.
        # Inputs derived from mda_per_line when available; else None.
        iso_dt: Optional[float] = None
        iso_dl: Optional[float] = None
        if mass_kg is not None and mass_kg > 0 and live > 0 and result.mda_per_line is not None:
            # Find MDA record for characteristic line (or first matched line)
            char_mda = None
            for m in ni.matched_lines:
                e_rounded = round(float(m.library_E_keV or 0), 2)
                key = (ni.nuclide, e_rounded)
                if key in mda_by_nuclide_energy:
                    char_mda = mda_by_nuclide_energy[key]
                    if m.is_characteristic:
                        break  # prefer the characteristic line
            if char_mda is not None:
                try:
                    bg_cnt = float(char_mda.background_counts_in_ROI or 0.0)
                    eff    = float(char_mda.efficiency or 0.0)
                    intens = float(char_mda.intensity or 0.0)  # decimal (not pct)
                    # gross_counts = net_peak_area + bg_counts.
                    # Find the matched line whose library_E_keV matches char_mda.
                    char_e_rounded = round(float(char_mda.line_energy_keV or 0), 2)
                    net_area: Optional[float] = None
                    for m in ni.matched_lines:
                        if round(float(m.library_E_keV or 0), 2) == char_e_rounded:
                            net_area = m.peak_area
                            break
                    if net_area is None:
                        # Fallback: first matched line with a measured area
                        for m in ni.matched_lines:
                            if m.peak_area is not None:
                                net_area = m.peak_area
                                break
                    gross_cnt = max(0.0, (net_area or 0.0)) + max(0.0, bg_cnt)
                    iso_dt = _iso_decision_threshold(
                        gross_cnt, bg_cnt, eff, intens, mass_kg, live
                    )
                    iso_dl = _iso_detection_limit(
                        gross_cnt, bg_cnt, eff, intens, mass_kg, live
                    )
                except Exception:
                    iso_dt = None
                    iso_dl = None

        out.append({
            "nuclide": ni.nuclide,
            "tier": ci_tier.get(ni.nuclide, "detected"),
            "n_matched_lines": len(ni.matched_lines),
            "characteristic_line_keV": _safe_float(ni.characteristic_line_keV),
            "confidence_index": _safe_float(
                ni.confidence.CI if ni.confidence is not None else None
            ),
            "confidence_level": ni.confidence_level()
                if ni.confidence is not None else "n/a",
            "peak_rate_cps": _safe_float(peak_rate),
            "activity_Bq": _safe_float(ar.A_Bq) if ar is not None else None,
            "activity_sigma_Bq": _safe_float(ar.sigma_A_Bq) if ar is not None else None,
            "activity_sigma_method": (ar.sigma_method if ar is not None else None),
            "activity_relative_sigma": _safe_float(sigma_rel),
            "is_upper_limit": is_upper_limit,
            "upper_limit_Bq": _safe_float(upper_limit_Bq),
            "specific_activity_Bq_per_kg":
                _safe_float(spec_act[0]) if spec_act is not None else None,
            "specific_activity_sigma_Bq_per_kg":
                _safe_float(spec_act[1]) if spec_act is not None else None,
            "decay_corrected": bool(ar.decay_corrected) if ar is not None else False,
            "decay_factor": _safe_float(
                ar.decay_factor if ar is not None else None
            ),
            "cascade_warning": (ar.cascade_warning if ar is not None else None) or "",
            "matched_lines": [
                {
                    "library_E_keV": _safe_float(m.library_E_keV),
                    "library_I_pct": _safe_float(m.library_I_pct),
                    "peak_E_keV": _safe_float(m.peak_E_keV),
                    "is_characteristic": m.is_characteristic,
                }
                for m in ni.matched_lines
            ],
            # F-RPT-08 / v1.19.1 — ISO 11929-1:2019 §5.4.3/§5.4.4 characteristic limits
            # in Bq/kg.  None when mass, efficiency, or mda data not available.
            # Cross-refs: RAG-005, RAG-008, RAG-009, RAG-022.
            "decision_threshold_Bq_per_kg": _safe_float(iso_dt),
            "detection_limit_Bq_per_kg": _safe_float(iso_dl),
            "reason": ni.reason,
        })
    return out


def _build_unidentified_peaks(result) -> List[Dict[str, Any]]:
    """Section 8 — true-unmatched residuals after classification."""
    out = []
    for rc in result.residual_classifications:
        if rc.label != "true_unmatched":
            continue
        out.append({
            "energy_keV": _safe_float(rc.peak_E_keV),
            "significance": _safe_float(rc.sigma),
            "label": rc.label,
            "note": rc.note,
        })
    return out


def _build_mda(result) -> List[Dict[str, Any]]:
    """Section 11 — ISO 11929 MDA table."""
    if result.mda_per_line is None:
        return []
    rows = []
    for (nuclide, E_keV_key), mda in result.mda_per_line.items():
        rows.append({
            "nuclide": nuclide,
            "line_E_keV": _safe_float(mda.line_energy_keV),
            "decision_threshold_counts": _safe_float(mda.decision_threshold_counts),
            "detection_limit_counts": _safe_float(mda.detection_limit_counts),
            "decision_threshold_cps": _safe_float(mda.decision_threshold_cps),
            "detection_limit_cps": _safe_float(mda.detection_limit_cps),
            "MDA_Bq": _safe_float(mda.MDA_Bq),
            "efficiency": _safe_float(mda.efficiency),
            "intensity_pct": _safe_float(mda.intensity * 100.0),
            "k_alpha": _safe_float(mda.k_alpha),
            "live_time_s": _safe_float(mda.live_time_s),
            "notes": mda.notes,
        })
    # Sort by (nuclide, line_E_keV)
    rows.sort(key=lambda r: (r["nuclide"] or "~", r["line_E_keV"] or 0.0))
    return rows


def _build_deconvolutions(result) -> List[Dict[str, Any]]:
    """Section 10 — multiplet deconvolution outputs."""
    if not result.deconvolution_results:
        return []
    out = []
    for d in result.deconvolution_results:
        components = []
        # F-145 / v1.17.8 — per-component centroid shift from Phase A.
        # Выровнено по components: shifts[k] = E_fitted - E_passport для k-ой
        # компоненты. Пустой tuple = Phase A не запускалась для кластера.
        shifts = list(getattr(d, "centroid_shifts_keV", ()) or ())
        for k, (comp, area, unc) in enumerate(zip(
            getattr(d, "components", []) or [],
            getattr(d, "areas", []) or [],
            getattr(d, "area_uncertainties", []) or [],
        )):
            comp_entry = {
                "nuclide": getattr(comp, "nuclide", ""),
                "line_E_keV": _safe_float(getattr(comp, "line_E_keV", None)),
                "center_channel": _safe_float(getattr(comp, "center_channel", None)),
                "fwhm_channels": _safe_float(getattr(comp, "fwhm_channels", None)),
                "deconvolved_area": _safe_float(area),
                "deconvolved_area_sigma": _safe_float(unc),
                # F-387.1 — provenance: active fit components отмечены
                # как "deconvolved" (или "deconvolved_coupled" — узнаётся
                # downstream HTML). Используем literal "deconvolved" для
                # consistency с identification.peak_area_source.
                "peak_area_source": "deconvolved",
            }
            # F-145: добавить centroid_shift_keV если Phase A была запущена
            if k < len(shifts):
                comp_entry["F145_centroid_shift_keV"] = _safe_float(shifts[k])
            components.append(comp_entry)
        # F-387.1 / v1.18.26.1 — подмешать phantom anchors (отрезанные
        # top-K cap после Rayleigh-CC split). Они НЕ участвовали в fit
        # (area=0, нет separate Gaussian), но сохраняются для evidence /
        # identification visibility.
        for comp in (getattr(d, "phantom_components", ()) or ()):
            components.append({
                "nuclide": getattr(comp, "nuclide", ""),
                "line_E_keV": _safe_float(getattr(comp, "line_E_keV", None)),
                "center_channel": _safe_float(
                    getattr(comp, "center_channel", None)
                ),
                "fwhm_channels": _safe_float(
                    getattr(comp, "fwhm_channels", None)
                ),
                "deconvolved_area": None,
                "deconvolved_area_sigma": None,
                "peak_area_source": "library_anchor_phantom",
                "is_phantom_anchor": True,
            })
        out.append({
            "cluster_id": str(getattr(d, "cluster_id", "") or ""),
            "converged": bool(getattr(d, "converged", False)),
            "chi2_per_dof": _safe_float(getattr(d, "chi2_per_dof", None)),
            # F-145 / v1.17.8 — Phase A χ²/ν (free-centroid pre-fit), даёт
            # opportunity diagnostic: если phase_A_chi2 ≪ chi2_per_dof, то
            # locked-passport теряет точность из-за дрейфа E(N).
            "F145_phase_A_chi2_per_dof": _safe_float(
                getattr(d, "phase_A_chi2_per_dof", None)
            ),
            "F145_phase_A_converged": bool(
                getattr(d, "phase_A_converged", False)
            ),
            "n_components": len(components),
            "components": components,
            "degenerate_pairs": list(getattr(d, "degenerate_pairs", []) or []),
            # F-392 / v1.18.27 + F-392.1 / v1.18.28 — surface chosen continuum
            # model so downstream consumers (HTML, run_skill.py sanity,
            # external dashboards) can see whether auto-select activated
            # step_linear_multi. Values: "linear", "step_linear",
            # "step_linear_multi", "quadratic" (or "" if absent).
            "continuum_model": str(getattr(d, "continuum_model", "") or ""),
            # F-392.1 / v1.18.29 — multi-step continuum diagnostics. Energies
            # (keV) of intense-anchor компонент, на которых ставится отдельный
            # β_step_i term в "step_linear_multi" модели; пустой список для
            # других continuum моделей. step_intensity_pct — порог library_I_pct
            # (typ. 4.0%); null когда multi-step не применялся.
            "step_anchor_energies_keV": [
                _safe_float(e)
                for e, _sigma in (getattr(d, "multi_step_anchors", ()) or ())
            ],
            "step_intensity_pct": _safe_float(
                getattr(d, "multi_step_intensity_threshold_pct", None)
            ),
        })
    return out


def _erfc_step_continuum(
    energies: List[float],
    e0: float,
    sigma: float,
    cont_left: float,
    cont_right: float,
) -> List[float]:
    """ГОСТ / Gilmore & Joss §6.5 — ступенчатая подстилающая под фотопиком.

    Подстилающая (continuum) под пиком — НЕ линейная хорда, а гладкая
    ступенька: высокая на низкоэнергетической (левой, комптоновской) стороне,
    плавно спадает ПОД пиком к низкому уровню справа. Форма = дополнительная
    функция ошибок erfc, центрированная на энергии пика e0 с масштабом перехода
    σ (= ширина пика):

        B(E) = cont_left·frac + cont_right·(1 − frac),
        frac = ½·erfc((E − e0)/(√2·σ))   — 1 на левом плече → 0 на правом.

    На E ≪ e0 → cont_left (касается чёрного спектра на ЛЕВОМ плече),
    на E ≫ e0 → cont_right (касается чёрного на ПРАВОМ плече),
    в e0 → среднее уровней плеч (continuum под пиком, ниже самого пика).

    cont_left / cont_right берутся из плеч ROI, привязанных к измеренному
    (gross) спектру → ступенька касается чёрной кривой на обоих плечах.
    DISPLAY-ONLY — квантование/активность/сертификат не затрагиваются.
    """
    if not energies:
        return []
    if sigma <= 0.0:
        # дегенерат (нет ширины): линейная хорда как fallback
        e_lo, e_hi = energies[0], energies[-1]
        span = (e_hi - e_lo) or 1.0
        return [cont_left + (cont_right - cont_left) * (E - e_lo) / span
                for E in energies]
    inv = 1.0 / (math.sqrt(2.0) * sigma)
    out: List[float] = []
    for E in energies:
        frac = 0.5 * math.erfc((E - e0) * inv)
        out.append(cont_left * frac + cont_right * (1.0 - frac))
    return out


def _build_fit_overlay(result) -> Dict[str, Any]:
    """F-FIT-VIEW v3 / v1.22.6 — build fit_overlay section for interactive HTML.

    Emits per-peak Gaussian parameters (singlets) and per-cluster continuum +
    total-fit arrays (multiplets) so the HTML toggle can reconstruct the
    full spectrum decomposition overlay.

    v2 additions (Task #66): covers secondary_peaks (source="secondary",
    orange), background_primary_feps (source="background", gray dashed),
    and unidentified_peaks (source="unidentified", yellow dashed).

    v3 additions (Task #67): each singlet/secondary/background/unidentified
    peak gains an optional ``continuum_grid`` field — a small (N≈11) array of
    {energies, values} sampled at ±3σ around the peak energy from the raw
    spectrum continuum shoulders.  The JS uses this to lift each Gaussian
    onto the local continuum baseline so the displayed photopeak matches the
    operator's visual reading of the raw spectrum.

    Returns a dict with:
      ``peaks``              — list of peaks with Gaussian params, source in:
                               singlet | multiplet_component | secondary |
                               background | unidentified
      ``multiplet_continua`` — list of per-cluster {E_keV, continuum, total,
                               components}
    """
    import math as _math
    SQRT_2PI = _math.sqrt(2.0 * _math.pi)

    # ── Helper: FWHM in keV from calibration model (a + b·E + c·E²) ──────────
    fwhm_model = getattr(result, "fwhm_model", None)
    def _fwhm_kev(E: float) -> float:
        if fwhm_model is None:
            return max(E * 0.07, 1.0)
        # F-452: polymorphic — FwhmModel callable | legacy 3-tuple.
        return _fwhm_keV_at_energy(fwhm_model, float(E))

    # ── Helper: rough amplitude from raw spectrum counts at given energy ───────
    spec = getattr(result, "spec", None)
    counts_arr = None
    if spec is not None:
        raw = getattr(spec, "counts", None)
        if raw is not None:
            try:
                import numpy as _np
                counts_arr = _np.asarray(raw, dtype=_np.float64)
            except Exception:
                counts_arr = None

    # Live time for activity-normalization of library_coverage peaks (v1.2.26).
    t_live = 1.0
    if spec is not None:
        _lt = getattr(spec, 'live_time', None)
        if _lt:
            try:
                _lt_f = float(_lt)
                if _lt_f > 0:
                    t_live = _lt_f
            except Exception:
                pass

    def _amp_from_spectrum(e: float, sigma: float) -> Optional[float]:
        """Estimate peak amplitude from raw spectrum counts near energy e.

        Uses the raw count rate at the nearest channel as a proxy for amplitude.
        Returns None when spectrum data is unavailable.
        """
        if spec is None or counts_arr is None:
            return None
        try:
            ch = spec.energy_to_channel(e)
            if ch is None:
                return None
            ch_idx = int(ch)
            n = len(counts_arr)
            if not (0 <= ch_idx < n):
                return None
            # Take max in a 1-sigma window around the centroid
            sigma_ch = 1.0
            try:
                ch_lo = spec.energy_to_channel(e - sigma)
                ch_hi = spec.energy_to_channel(e + sigma)
                if ch_lo is not None and ch_hi is not None:
                    lo = max(0, int(ch_lo))
                    hi = min(n, int(ch_hi) + 1)
                    if hi > lo:
                        sigma_ch = max(1.0, (hi - lo) / 2.0)
            except Exception:
                pass
            # Use the channel count as amplitude proxy (not area-normalised)
            val = float(counts_arr[ch_idx])
            return round(val, 2) if val > 0 else None
        except Exception:
            return None

    def _amp_net_from_spectrum(e: float, sigma: float) -> Optional[float]:
        """Net peak amplitude at energy e above local continuum estimate.

        Peak max in [E-sigma, E+sigma] minus mean of left/right shoulders
        sampled at the integration-ROI baseline wings ([1.75*FWHM, 2.5*FWHM]
        each side).  Returns None when unavailable or net <= 0.
        """
        if spec is None or counts_arr is None:
            return None
        try:
            import numpy as _npn
            n = len(counts_arr)
            ch_p_lo = spec.energy_to_channel(e - sigma)
            ch_p_hi = spec.energy_to_channel(e + sigma)
            if ch_p_lo is None or ch_p_hi is None:
                return None
            p_lo = max(0, int(ch_p_lo))
            p_hi = min(n, max(p_lo + 1, int(ch_p_hi) + 1))
            if p_hi <= p_lo:
                return None
            peak_count = float(_npn.max(counts_arr[p_lo:p_hi]))
            cont_vals: list = []
            # Sync continuum window to integration ROI (area.py window_factor=2.5):
            # ROI half-width = 2.5*FWHM; baseline wings = outer 30% (Cowell §5.2.5)
            # → shoulders at [1.75*FWHM, 2.5*FWHM] each side.  FWHM = 2.355*sigma.
            fwhm = sigma * 2.355
            roi_half = 2.5 * fwhm
            wing_inner = 1.75 * fwhm
            for (ca, cb) in [
                (e - roi_half, e - wing_inner),
                (e + wing_inner, e + roi_half),
            ]:
                ch_a = spec.energy_to_channel(ca)
                ch_b = spec.energy_to_channel(cb)
                if ch_a is None or ch_b is None:
                    continue
                a = max(0, int(ch_a))
                b = min(n, max(a + 1, int(ch_b)))
                if b > a:
                    cont_vals.append(float(_npn.mean(counts_arr[a:b])))
            # Use mean of two shoulder estimates for continuum — matches the linear
            # interpolation used in continuum_grid, so that g + continuum_at_peak
            # equals the measured spectrum value at the peak channel.
            cont_est = float(_npn.mean(cont_vals)) if cont_vals else 0.0
            net = peak_count - cont_est
            return round(max(0.0, net), 2) if net > 0 else None
        except Exception:
            return None

    def _dE_per_channel(e: float) -> Optional[float]:
        # F-FIT-VIEW units fix (2026-06-20): local keV/ch. Singlet/bg Gaussian is
        # drawn in keV but area is counts over CHANNELS → height needs dE/dch:
        # area·(dE/dch)/(σ_keV·√2π) = area/(σ_ch·√2π) (else ~3× short at 2614).
        if spec is None:
            return None
        try:
            ch = spec.energy_to_channel(e)
            if ch is None:
                return None
            ch_i = int(round(float(ch)))
            e0 = spec.channel_to_energy(ch_i)
            e1 = spec.channel_to_energy(ch_i + 1)
            if e0 is None or e1 is None:
                return None
            d = abs(float(e1) - float(e0))
            return d if d > 0 else None
        except Exception:
            return None

    def _continuum_grid_for_peak(e: float, sigma: float) -> Optional[Dict[str, Any]]:
        """F-FIT-VIEW v3 — build continuum_grid for a singlet/non-multiplet peak.

        Samples the raw spectrum at the integration-ROI baseline wings
        ([1.75*FWHM, 2.5*FWHM] each side) and produces N≈11 continuum sample
        points across [E-2.5*FWHM, E+2.5*FWHM] shaped as a ГОСТ / Gilmore §6.5
        erfc STEP (high on the left Compton side → drops under the peak → low on
        the right), not a straight chord.  Window synced to area.py
        window_factor=2.5.  Used by the JS to lift the Gaussian display onto the
        local Compton/scatter floor.

        Returns None when spectrum data is unavailable (JS falls back to y=0).
        """
        if spec is None or counts_arr is None:
            return None
        try:
            import numpy as _np
            n = len(counts_arr)
            # Continuum window synced to integration ROI (area.py window_factor=2.5):
            # half-width 2.5*FWHM, baseline wings = outer 30% (Cowell §5.2.5).
            fwhm = sigma * 2.355
            roi_half = 2.5 * fwhm
            wing_inner = 1.75 * fwhm
            # Left shoulder: mean of counts in [E-2.5*FWHM, E-1.75*FWHM]
            ch_lo = spec.energy_to_channel(e - roi_half)
            ch_hi = spec.energy_to_channel(e - wing_inner)
            if ch_lo is None or ch_hi is None:
                return None
            lo_lo = max(0, int(ch_lo))
            lo_hi = max(lo_lo + 1, int(ch_hi))
            lo_hi = min(n, lo_hi)
            if lo_hi <= lo_lo:
                return None
            cont_left = float(_np.mean(counts_arr[lo_lo:lo_hi]))
            # Right shoulder: mean of counts in [E+1.75*FWHM, E+2.5*FWHM]
            ch_rl = spec.energy_to_channel(e + wing_inner)
            ch_rh = spec.energy_to_channel(e + roi_half)
            if ch_rl is None or ch_rh is None:
                return None
            ri_lo = min(n - 1, max(0, int(ch_rl)))
            ri_hi = min(n, max(ri_lo + 1, int(ch_rh)))
            if ri_hi <= ri_lo:
                return None
            cont_right = float(_np.mean(counts_arr[ri_lo:ri_hi]))
            # Build 11 sample points linearly spaced across [E-2.5*FWHM, E+2.5*FWHM]
            e_lo = e - roi_half
            e_hi = e + roi_half
            n_pts = 11
            energies = [round(e_lo + i * (e_hi - e_lo) / (n_pts - 1), 2)
                        for i in range(n_pts)]
            # ГОСТ / Gilmore §6.5 erfc STEP between the two shoulder levels
            # (high on the left Compton side → drops under the peak → low on the
            # right). Чистая erfc-step, без clamp к data: clamp давал «синусоиду»
            # в провалах между пиками (display-artifact, операторская диагностика
            # 2026-06-21 на Tl-208 2614: 543→497→326→359→489→316 — 2 экстремума
            # вместо монотонного спада). DISPLAY-ONLY — quantification untouched.
            step = _erfc_step_continuum(energies, e, sigma, cont_left, cont_right)
            values = [round(max(0.0, float(step[i])), 3) for i in range(n_pts)]
            return {"energies": energies, "values": values}
        except Exception:
            return None

    def _grouped_continuum_grids(peak_list):
        """Gilmore §9.7 / ORTEC §6.5.1 — shared continuum for groups of close peaks.

        Groups peaks within 2·FWHM of each other (multiplet criterion).
        Baseline anchors are sampled OUTSIDE the full group zone so overlapping
        peaks share one linear baseline instead of each sampling its own shoulder
        (which would land on a neighbour tail and inflate the continuum estimate).

        peak_list: list of (energy_keV, sigma_keV).
        Returns dict {energy_keV: continuum_grid}.
        """
        if not peak_list or spec is None or counts_arr is None:
            return {}
        try:
            import numpy as _np
            _nc = len(counts_arr)
            _srt = sorted(peak_list, key=lambda p: p[0])
            groups: list = []
            _cur: list = [_srt[0]]
            for _gi in range(1, len(_srt)):
                _pe, _ps = _srt[_gi - 1]
                _ce, _cs = _srt[_gi]
                if _ce - _pe <= 2.0 * 2.355 * (_ps + _cs) / 2.0:
                    _cur.append(_srt[_gi])
                else:
                    groups.append(_cur)
                    _cur = [_srt[_gi]]
            groups.append(_cur)

            try:
                _y = _np.maximum(_np.asarray(counts_arr, dtype=float), 0.0)
                _v = _np.log(_np.log(_np.sqrt(_y + 1.0) + 1.0) + 1.0)
                for _m in range(1, 25):
                    _vp = _v.copy()
                    if _nc > 2 * _m:
                        _v[_m:_nc - _m] = _np.minimum(
                            _vp[_m:_nc - _m],
                            0.5 * (_vp[:_nc - 2 * _m] + _vp[2 * _m:]))
                _bl_arr = (_np.exp(_np.exp(_v) - 1.0) - 1.0)
                _bl_arr = _np.maximum(_bl_arr * _bl_arr - 1.0, 0.0)
            except Exception:
                _bl_arr = _np.zeros(_nc, dtype=float)

            def _ch_floor(ea, eb):
                ca = spec.energy_to_channel(ea)
                cb = spec.energy_to_channel(eb)
                if ca is None or cb is None:
                    return 0.0
                a = max(0, int(ca))
                b = min(_nc, max(a + 1, int(cb) + 1))
                if b <= a:
                    return 0.0
                try:
                    return float(_np.mean(_bl_arr[a:b]))
                except Exception:
                    return 0.0

            def _snip_at(eg):
                _ch = spec.energy_to_channel(eg)
                if _ch is None:
                    return 0.0
                _ci = max(0, min(_nc - 1, int(_ch)))
                try:
                    return float(max(0.0, _bl_arr[_ci]))
                except Exception:
                    return 0.0

            _res: Dict[float, Any] = {}
            for _grp in groups:
                for (_ep, _sp_pk) in _grp:
                    _gl = _ep - 4.0 * _sp_pk
                    _gh = _ep + 4.0 * _sp_pk
                    _npt = 17
                    _en = [_gl + _j * (_gh - _gl) / (_npt - 1) for _j in range(_npt)]
                    _vv = [round(_snip_at(_eg), 3) for _eg in _en]
                    _res[_ep] = {"energies": [round(_eg, 2) for _eg in _en],
                                 "values": _vv}
            return _res
        except Exception:
            return {}

    # Activity-normalization helpers (v1.2.26) ---------------------------------
    # K_eff(E) = area/(I_gamma*t_live): implicit Activity*eps(E)*m product.
    # Scales lib_coverage lines to I_gamma-consistent amplitudes.

    def _I_gamma_for_line(nuc_name, energy_keV, tol=3.0):
        """Return tabulated gamma intensity fraction for closest nuclide line."""
        try:
            from gamma.data.nuclide_library import get_nuclide as _gnl
            nrec = _gnl(nuc_name)
            if not nrec:
                return 0.0
            best_i, best_de = 0.0, 1e9
            lines_key = 'lines'
            for ll in (nrec.get(lines_key) or []):
                le = float(ll[0]) if ll and len(ll) > 0 else None
                li = float(ll[1]) if ll and len(ll) > 1 else 0.0
                if le is None:
                    continue
                de = abs(le - energy_keV)
                if de < best_de and de <= tol:
                    best_de, best_i = de, li
            return best_i / 100.0
        except Exception:
            return 0.0


    def _build_nuc_K_eff(pk_list, tl):
        """Per-nuclide K_eff(E) from fitted singlet/multiplet_component peaks.

        K_eff = area / (I_gamma * t_live)  -- implicit Activity*eps(E)*m_kg.
        Returns {nuclide: sorted [(E, K_eff)]}.
        """
        import math as _mh
        _r = {}
        tl = tl if tl and tl > 0 else 1.0
        for pk in pk_list:
            if pk.get('source') not in ('singlet', 'multiplet_component'):
                continue
            nuc = pk.get('nuclide', '')
            e_pk = float(pk.get('energy_keV', 0.0))
            amp_pk = float(pk.get('amp_counts', 0.0) or 0.0)
            sig_pk = float(pk.get('sigma_keV', 0.0) or 0.0)
            if not nuc or amp_pk <= 0 or sig_pk <= 0:
                continue
            I_g = _I_gamma_for_line(nuc, e_pk)
            # Require I_g >= 1% to prevent very weak lines contaminating K_eff (v1.2.27).
            # Tiny I_g denominators produce K_eff 10-300x too high, distorting lib_coverage amps.
            if I_g < 0.01:
                continue
            area_pk = amp_pk * sig_pk * _mh.sqrt(2.0 * _mh.pi)
            K = area_pk / (I_g * tl)
            _r.setdefault(nuc, []).append((e_pk, K))
        for _nk in _r:
            _r[_nk].sort(key=lambda x: x[0])
        return _r

    def _interp_K_eff(k_pts, E):
        """Linearly interpolate K_eff at energy E from [(E, K)] list."""
        if not k_pts:
            return None
        if len(k_pts) == 1:
            return k_pts[0][1]
        k_s = sorted(k_pts, key=lambda x: x[0])
        e_v = [p[0] for p in k_s]
        k_v = [p[1] for p in k_s]
        if E <= e_v[0]:
            if len(k_s) >= 2:
                sl = (k_v[1] - k_v[0]) / max(e_v[1] - e_v[0], 0.1)
                return max(0.1, k_v[0] + sl * (E - e_v[0]))
            return k_v[0]
        if E >= e_v[-1]:
            if len(k_s) >= 2:
                sl = (k_v[-1] - k_v[-2]) / max(e_v[-1] - e_v[-2], 0.1)
                return max(0.1, k_v[-1] + sl * (E - e_v[-1]))
            return k_v[-1]
        for idx in range(len(k_s) - 1):
            if e_v[idx] <= E <= e_v[idx + 1]:
                f = (E - e_v[idx]) / max(e_v[idx + 1] - e_v[idx], 0.1)
                return k_v[idx] + f * (k_v[idx + 1] - k_v[idx])
        return k_v[-1]

    def _build_global_composite(pk_list, mc_list, spec_obj, cnts_arr):
        """Global composite: sum of all fitted Gaussians + continuum.

        Returns {energies:[...], values:[...]} in counts (400 points). None on failure.
        Continuum: from cluster continua where available, else smoothed spectrum.
        """
        if spec_obj is None or cnts_arr is None:
            return None
        try:
            import numpy as _np_gc
            n_ch = len(cnts_arr)
            try:
                e_start = float(spec_obj.channel_to_energy(0))
                e_end = float(spec_obj.channel_to_energy(n_ch - 1))
            except Exception:
                e_start, e_end = 25.0, 3000.0
            N = 400
            E_grid = _np_gc.linspace(max(e_start, 25.0), min(e_end, 3100.0), N)
            gauss_sum = _np_gc.zeros(N)
            for pk in pk_list:
                amp = float(pk.get('amp_counts', 0.0) or 0.0)
                sig = float(pk.get('sigma_keV', 0.0) or 0.0)
                E0 = float(pk.get('energy_keV', 0.0) or 0.0)
                if amp <= 0 or sig <= 0:
                    continue
                dE = E_grid - E0
                gauss_sum += amp * _np_gc.exp(-dE * dE / (2.0 * sig * sig))
            continuum = _np_gc.zeros(N)
            cont_covered = _np_gc.zeros(N, dtype=bool)
            for mc in mc_list:
                E_cl = mc.get('E_keV') or []
                cont_cl = mc.get('continuum') or []
                if not E_cl or not cont_cl or len(E_cl) != len(cont_cl):
                    continue
                E_cl_a = _np_gc.array(E_cl, dtype=float)
                ct_cl_a = _np_gc.array(cont_cl, dtype=float)
                e_lo_c = float(E_cl_a[0])
                e_hi_c = float(E_cl_a[-1])
                mask = (E_grid >= e_lo_c) & (E_grid <= e_hi_c)
                if mask.any():
                    continuum[mask] = _np_gc.interp(E_grid[mask], E_cl_a, ct_cl_a)
                    cont_covered[mask] = True
            uncov = ~cont_covered
            if uncov.any():
                sw = max(20, n_ch // 50)
                kernel = _np_gc.ones(2 * sw + 1) / (2 * sw + 1)
                cs = _np_gc.convolve(cnts_arr, kernel, mode='same')
                for ii in _np_gc.where(uncov)[0]:
                    try:
                        ch = spec_obj.energy_to_channel(float(E_grid[ii]))
                        if ch is not None:
                            ci = max(0, min(n_ch - 1, int(ch)))
                            continuum[ii] = float(cs[ci])
                    except Exception:
                        pass
            total = gauss_sum + continuum
            return {
                'energies': [round(float(e), 2) for e in E_grid],
                'values': [round(float(v), 3) for v in total],
            }
        except Exception:
            return None




    # ── 1. Singlet peaks (primary_feps not in multiplet clusters) ────────────
    # Collect energy ranges covered by multiplet ROIs so we don't double-emit.
    decons = getattr(result, "deconvolution_results", None) or []
    multiplet_E_ranges = []
    for d in decons:
        roi_lo = int(getattr(d, "roi_low_ch", 0))
        roi_hi = int(getattr(d, "roi_high_ch", 0))
        if spec is not None and roi_lo < roi_hi:
            try:
                e_lo = float(spec.channel_to_energy(roi_lo))
                e_hi = float(spec.channel_to_energy(roi_hi))
                multiplet_E_ranges.append((min(e_lo, e_hi), max(e_lo, e_hi)))
            except Exception:
                pass

    def _in_multiplet(e: float) -> bool:
        for lo, hi in multiplet_E_ranges:
            if lo - 5.0 <= e <= hi + 5.0:
                return True
        return False

    peaks_out: List[Dict[str, Any]] = []
    PHANTOM_SOURCES = {"library_anchor", "library_anchor_phantom"}
    seen_peaks: set = set()  # deduplicate by (nuclide, round(peak_E_keV, 0))
    for ni in (getattr(result, "final_detected", None) or []):
        nuc = ni.nuclide
        for m in ni.matched_lines:
            if str(m.peak_area_source or "") in PHANTOM_SOURCES:
                continue
            e = _safe_float(m.peak_E_keV)
            if e is None:
                continue
            area = _safe_float(m.peak_area) or 0.0
            if area <= 0:
                continue
            # Deduplicate: same nuclide + same energy (within 1 keV)
            dedup_key = (nuc, round(e))
            if dedup_key in seen_peaks:
                continue
            seen_peaks.add(dedup_key)
            # fwhm: prefer gauss_sigma_keV * 2.355, else legacy peak_sigma
            _sigma_attr = getattr(m, "gauss_sigma_keV", None)
            if _sigma_attr is None:
                _sigma_attr = getattr(m, "peak_sigma", None)
            if _sigma_attr and float(_sigma_attr) > 0:
                sigma = float(_sigma_attr)
            else:
                sigma = _fwhm_kev(e) / 2.355
            # amp = counts/channel peak HEIGHT; area is counts over channels, sigma
            # is keV → multiply by local dE/dch (F-FIT-VIEW units fix). DISPLAY-ONLY:
            # quantification uses peak_area_counts directly, cert untouched.
            amp = area / (sigma * SQRT_2PI)
            _dE = _dE_per_channel(e)
            if _dE and _dE > 0:
                amp *= _dE
            source = "singlet" if not _in_multiplet(e) else "multiplet_component"
            # Build compact label: «Cs-137 662»
            lib_e = _safe_float(m.library_E_keV)
            label = f"{nuc} {round(lib_e)}" if lib_e else f"{nuc} {round(e)}"
            # v3: attach continuum_grid for non-multiplet peaks so JS can lift
            # the Gaussian onto the local Compton/scatter floor.
            # multiplet_component peaks are rendered via g_curve (continuum handled
            # by multiplet_continua cluster array) — skip continuum_grid for them.
            cg = None if source == "multiplet_component" else _continuum_grid_for_peak(e, sigma)
            entry: Dict[str, Any] = {
                "peak_id": f"p{round(e)}",
                "nuclide": nuc,
                "energy_keV": round(e, 2),
                "amp_counts": round(amp, 2),
                "sigma_keV": round(sigma, 3),
                "source": source,
                "label": label,
            }
            if cg is not None:
                entry["continuum_grid"] = cg
            peaks_out.append(entry)

    # G4 / v1.31.0 — annotate sample-side peaks that also appear in the
    # background identification. We do this **before** secondary/bg/unid
    # entries are appended so the catalog match runs against a clean
    # sample-side slice (sources singlet | multiplet_component) and operator
    # tooling can rely on the annotation regardless of later peak emissions.
    _bg_staged_for_carryover = getattr(result, "background_staged_result", None)
    _bg_catalog = _build_bg_peak_catalog(_bg_staged_for_carryover)
    if _bg_catalog:
        _mark_bg_carryover(peaks_out, _bg_catalog, _fwhm_kev)

    # Per-nuclide K_eff normalization table from fitted peaks (v1.2.26).
    _nuc_K_table = _build_nuc_K_eff(peaks_out, t_live)

    # ── 2. Secondary peaks (residual classifications — not xrf / unmatched) ──
    # source="secondary" — orange Gaussians in frontend overlay.
    # RC objects have peak_E_keV and significance (S/N), but no direct area.
    # Use _amp_from_spectrum() to get a visual height proxy.
    SECONDARY_SKIP_LABELS = {"xrf", "true_unmatched", "edge_of_range"}
    seen_secondary: set = set()
    for rc in (getattr(result, "residual_classifications", None) or []):
        if getattr(rc, "label", "") in SECONDARY_SKIP_LABELS:
            continue
        e = _safe_float(getattr(rc, "peak_E_keV", None))
        if e is None:
            continue
        dedup_key = round(e)
        if dedup_key in seen_secondary:
            continue
        seen_secondary.add(dedup_key)
        sigma = _fwhm_kev(e) / 2.355
        amp = _amp_from_spectrum(e, sigma)
        feat = getattr(rc, "feature_kind", "") or getattr(rc, "label", "") or ""
        parent = getattr(rc, "parent_nuclide", "") or ""
        label_parts = [feat or "secondary"]
        if parent:
            label_parts.append(parent)
        label_parts.append(str(round(e)))
        label = " ".join(label_parts)
        sec_entry: Dict[str, Any] = {
            "peak_id": f"ps{round(e)}",
            "nuclide": parent or feat or "?",
            "energy_keV": round(e, 2),
            "amp_counts": amp,
            "sigma_keV": round(sigma, 3),
            "source": "secondary",
            "label": label,
        }
        # v3: lift secondary Gaussians onto local continuum
        sec_cg = _continuum_grid_for_peak(e, sigma)
        if sec_cg is not None:
            sec_entry["continuum_grid"] = sec_cg
        peaks_out.append(sec_entry)

    # ── 3. Background primary FEPs (source="background", gray dashed) ─────────
    # From bg_result.final_detected — same structure as primary FEPs.
    bg_result = getattr(result, "background_staged_result", None)
    if bg_result is not None:
        bg_fwhm_model = getattr(bg_result, "fwhm_model", None)
        def _bg_fwhm_kev(E: float) -> float:
            if bg_fwhm_model is None:
                return max(E * 0.07, 1.0)
            # F-452: polymorphic — FwhmModel callable | legacy 3-tuple.
            return _fwhm_keV_at_energy(bg_fwhm_model, float(E))

        # v3: background spectrum continuum helper (uses bg_result.spec)
        bg_spec = getattr(bg_result, "spec", None)
        bg_counts_arr = None
        if bg_spec is not None:
            _bg_raw = getattr(bg_spec, "counts", None)
            if _bg_raw is not None:
                try:
                    import numpy as _np2
                    bg_counts_arr = _np2.asarray(_bg_raw, dtype=_np2.float64)
                except Exception:
                    bg_counts_arr = None

        def _bg_continuum_grid(e: float, sigma: float) -> Optional[Dict[str, Any]]:
            if bg_spec is None or bg_counts_arr is None:
                return None
            try:
                import numpy as _np3
                n_bg = len(bg_counts_arr)
                # Continuum window synced to integration ROI (area.py window_factor=2.5):
                # half-width 2.5*FWHM, baseline wings = outer 30% (Cowell §5.2.5).
                fwhm = sigma * 2.355
                roi_half = 2.5 * fwhm
                wing_inner = 1.75 * fwhm
                ch_lo = bg_spec.energy_to_channel(e - roi_half)
                ch_hi = bg_spec.energy_to_channel(e - wing_inner)
                if ch_lo is None or ch_hi is None:
                    return None
                lo_lo = max(0, int(ch_lo))
                lo_hi = max(lo_lo + 1, min(n_bg, int(ch_hi)))
                if lo_hi <= lo_lo:
                    return None
                cont_left = float(_np3.mean(bg_counts_arr[lo_lo:lo_hi]))
                ch_rl = bg_spec.energy_to_channel(e + wing_inner)
                ch_rh = bg_spec.energy_to_channel(e + roi_half)
                if ch_rl is None or ch_rh is None:
                    return None
                ri_lo = max(0, int(ch_rl))
                ri_hi = min(n_bg, max(ri_lo + 1, int(ch_rh)))
                if ri_hi <= ri_lo:
                    return None
                cont_right = float(_np3.mean(bg_counts_arr[ri_lo:ri_hi]))
                e_lo = e - roi_half
                e_hi = e + roi_half
                n_pts = 11
                energies = [round(e_lo + i * (e_hi - e_lo) / (n_pts - 1), 2)
                            for i in range(n_pts)]
                # ГОСТ / Gilmore §6.5 erfc STEP без clamp: clamp давал «синусоиду»
                # в провалах между пиками (display-artifact). DISPLAY-ONLY
                # (зеркало sample-helper'а _continuum_grid_for_peak, 2026-06-21).
                step = _erfc_step_continuum(energies, e, sigma, cont_left, cont_right)
                values = [round(max(0.0, float(step[i])), 3) for i in range(n_pts)]
                return {"energies": energies, "values": values}
            except Exception:
                return None

        seen_bg: set = set()
        for ni in (getattr(bg_result, "final_detected", None) or []):
            nuc = ni.nuclide
            for m in ni.matched_lines:
                if str(getattr(m, "peak_area_source", "") or "") in PHANTOM_SOURCES:
                    continue
                e = _safe_float(getattr(m, "peak_E_keV", None))
                if e is None:
                    continue
                area = _safe_float(getattr(m, "peak_area", None)) or 0.0
                if area <= 0:
                    continue
                dedup_key = (nuc, round(e))
                if dedup_key in seen_bg:
                    continue
                seen_bg.add(dedup_key)
                _sigma_attr = getattr(m, "gauss_sigma_keV", None)
                if _sigma_attr is None:
                    _sigma_attr = getattr(m, "peak_sigma", None)
                if _sigma_attr and float(_sigma_attr) > 0:
                    sigma = float(_sigma_attr)
                else:
                    sigma = _bg_fwhm_kev(e) / 2.355
                # Same dE/dch height correction as the singlet path — background
                # peaks render via gaussianPoints too (F-FIT-VIEW units fix).
                amp = area / (sigma * SQRT_2PI)
                _dE = _dE_per_channel(e)
                if _dE and _dE > 0:
                    amp *= _dE
                lib_e = _safe_float(getattr(m, "library_E_keV", None))
                label = (f"{nuc} (bg) {round(lib_e)}"
                         if lib_e else f"{nuc} (bg) {round(e)}")
                bg_entry: Dict[str, Any] = {
                    "peak_id": f"pb{round(e)}",
                    "nuclide": f"{nuc} (bg)",
                    "energy_keV": round(e, 2),
                    "amp_counts": round(amp, 2),
                    "sigma_keV": round(sigma, 3),
                    "source": "background",
                    "label": label,
                }
                # v3: lift background Gaussians onto bg-spectrum continuum
                bg_cg = _bg_continuum_grid(e, sigma)
                if bg_cg is not None:
                    bg_entry["continuum_grid"] = bg_cg
                peaks_out.append(bg_entry)

    # ── 4. Unidentified peaks (true_unmatched residuals) ──────────────────────
    # source="unidentified" — yellow dashed Gaussians.
    seen_unident: set = set()
    for rc in (getattr(result, "residual_classifications", None) or []):
        if getattr(rc, "label", "") != "true_unmatched":
            continue
        e = _safe_float(getattr(rc, "peak_E_keV", None))
        if e is None:
            continue
        dedup_key = round(e)
        if dedup_key in seen_unident:
            continue
        seen_unident.add(dedup_key)
        sigma = _fwhm_kev(e) / 2.355
        amp = _amp_from_spectrum(e, sigma)
        label = f"? {round(e)}"
        unid_entry: Dict[str, Any] = {
            "peak_id": f"pu{round(e)}",
            "nuclide": "?",
            "energy_keV": round(e, 2),
            "amp_counts": amp,
            "sigma_keV": round(sigma, 3),
            "source": "unidentified",
            "label": label,
        }
        # v3: lift unidentified Gaussians onto local continuum
        unid_cg = _continuum_grid_for_peak(e, sigma)
        if unid_cg is not None:
            unid_entry["continuum_grid"] = unid_cg
        peaks_out.append(unid_entry)

    # ── 4.5 Library-coverage Gaussians (confirmed nuclide lines absent from overlay) ──
    # source="library_coverage" -- green dashed Gaussians for significant (I >= 0.5%)
    # confirmed-nuclide lines that have zero/missing amplitude in the overlay due to:
    # phantom-zeroing, being outside TH232_FORCED_CLUSTERS component list, or
    # simply unmatched in the identification step.  Amplitude = net peak counts
    # above local continuum, so the Gaussian height matches the actual spectrum peak.
    # Examples: Ac-228 129 keV (in bg cluster only), Ac-228 726 keV (phantom-zeroed),
    # Ac-228 794.95 keV (missing from TH232_FORCED_CLUSTERS M1 components entirely).
    LIB_COV_I_MIN = 0.5

    def _in_overlay_already(lib_e: float) -> bool:
        """True if any peak with positive amplitude is within FWHM/2 of lib_e.

        Tolerance-based check instead of integer rounding so that a fitted peak
        at a calibration-shifted position (e.g. Pb-212 at 234.64 vs library
        238.63 keV) correctly suppresses a duplicate library_coverage entry.
        Also catches within-section duplicates since we append to peaks_out.
        """
        _tol = _fwhm_kev(lib_e) * 0.5
        for _pk in peaks_out:
            # background-source peaks represent the background spectrum Gaussian,
            # not a fitted sample line; they should not suppress library_coverage.
            if _pk.get("source") == "background":
                continue
            if (_pk.get("amp_counts") or 0.0) > 0.0:
                if abs(_pk.get("energy_keV", -99999.0) - lib_e) <= _tol:
                    return True
        return False

    # Pass 1: collect all lib_coverage candidates for grouped continuum.
    _lc_cands: list = []
    for _ni in (getattr(result, "final_detected", None) or []):
        _nuc_p1 = _ni.nuclide
        try:
            from gamma.data.nuclide_library import get_nuclide as _gn_p1
            _nrec_p1 = _gn_p1(_nuc_p1)
        except Exception:
            _nrec_p1 = None
        if _nrec_p1 is None:
            continue
        for _ll in (_nrec_p1.get("lines") or []):
            _le = float(_ll[0]) if len(_ll) > 0 else None
            _li = float(_ll[1]) if len(_ll) > 1 else 0.0
            if _le is None or _li < LIB_COV_I_MIN:
                continue
            if _in_overlay_already(_le):
                continue
            _sg = _fwhm_kev(_le) / 2.355
            if (_amp_net_from_spectrum(_le, _sg) or 0) <= 0:
                continue
            _lc_cands.append((_le, _sg, _nuc_p1))
    # Pass 2: build grouped continua (Gilmore §9.7), then create entries.
    _lc_cont = _grouped_continuum_grids([(_e, _s) for _e, _s, _ in _lc_cands])
    _lc_seen: set = set()
    for (_lib_E, _sigma_cov, _nuc) in _lc_cands:
        _key = round(_lib_E, 2)
        if _key in _lc_seen:
            continue
        _lc_seen.add(_key)
        # Activity-normalized amplitude (v1.2.26): K_eff * I_gamma * t_live / (sigma*sqrt2pi)
        _I_g_cov = _I_gamma_for_line(_nuc, _lib_E)
        _K_cov = _interp_K_eff(_nuc_K_table.get(_nuc, []), _lib_E)
        if _K_cov is not None and _I_g_cov > 0 and _sigma_cov > 0:
            import math as _mc26
            _amp_cov = round(
                max(0.1, _K_cov * _I_g_cov * t_live / (_sigma_cov * _mc26.sqrt(2.0 * _mc26.pi))),
                2,
            )
        else:
            _amp_cov = _amp_net_from_spectrum(_lib_E, _sigma_cov)
            if _amp_cov is None or _amp_cov <= 0:
                continue
        _cg_cov = _lc_cont.get(_lib_E) or _continuum_grid_for_peak(_lib_E, _sigma_cov)
        _cov_entry: Dict[str, Any] = {
            "peak_id": f"plc{round(_lib_E)}",
            "nuclide": _nuc,
            "energy_keV": _key,
            "amp_counts": round(_amp_cov, 2),
            "sigma_keV": round(_sigma_cov, 3),
            "source": "library_coverage",
            "label": f"{_nuc} {round(_lib_E)}",
        }
        if _cg_cov is not None:
            _cov_entry["continuum_grid"] = _cg_cov
        peaks_out.append(_cov_entry)

    # ── 5. Multiplet continua (per-cluster overlay arrays) ────────────────────
    multiplet_continua_out: List[Dict[str, Any]] = []
    for d in decons:
        has_overlay = bool(
            getattr(d, "overlay_E_keV", None)
            and getattr(d, "overlay_continuum", None)
            and getattr(d, "overlay_total", None)
        )
        if not has_overlay:
            continue
        E_arr_raw = [float(v) for v in d.overlay_E_keV]
        cont_arr_raw = [float(v) for v in d.overlay_continuum]
        E_arr = [round(v, 2) for v in E_arr_raw]
        cont_arr = [round(v, 3) for v in cont_arr_raw]

        comp_list: List[Dict[str, Any]] = []
        overlays = getattr(d, "overlay_components", ()) or ()
        # F-145: per-component fitted centroid shift (Phase A free-centroid fit).
        shifts_kev = list(getattr(d, "centroid_shifts_keV", ()) or ())
        # v1.2.18 cache for shifted-render rebuild (E_fit, amp, sigma)
        _rebuild_specs: List[Tuple[float, float, float]] = []
        for k, (comp, area) in enumerate(
            zip(getattr(d, "components", []) or [],
                getattr(d, "areas", []) or [])
        ):
            nuc = getattr(comp, "nuclide", "") or ""
            E_line = float(getattr(comp, "line_E_keV", 0.0) or 0.0)
            shift_k = float(shifts_kev[k]) if k < len(shifts_kev) else 0.0
            E_fit = E_line + shift_k
            area_f = float(area or 0.0)
            if area_f <= 0:
                continue
            # sigma from fwhm_channels (in keV domain via spec channel_to_energy)
            fwhm_ch = _safe_float(getattr(comp, "fwhm_channels", None))
            if fwhm_ch and fwhm_ch > 0 and spec is not None:
                try:
                    ch_center = int(getattr(comp, "center_channel", None) or
                                    spec.energy_to_channel(E_line))
                    e0 = float(spec.channel_to_energy(ch_center))
                    e1 = float(spec.channel_to_energy(ch_center + int(fwhm_ch)))
                    sigma = abs(e1 - e0) / 2.355
                except Exception:
                    sigma = _fwhm_kev(E_line) / 2.355
            else:
                sigma = _fwhm_kev(E_line) / 2.355
            amp = area_f / (max(sigma, 0.01) * SQRT_2PI)
            # Per-component g-curve (g_plus_cont - continuum = Gaussian only)
            if k < len(overlays):
                g_pure = [
                    round(max(float(g) - float(b), 0.0), 3)
                    for g, b in zip(overlays[k], d.overlay_continuum)
                ]
            else:
                g_pure = []
            label = f"{nuc} {round(E_line)}"
            # v1.2.18 — cache spec for shifted-render rebuild below
            _rebuild_specs.append((E_fit, amp, sigma))
            comp_list.append({
                "nuclide": nuc,
                "energy_keV": round(E_line, 2),
                "fit_centroid_keV": round(E_fit, 3),
                "centroid_shift_keV": round(shift_k, 3),
                "sigma_keV": round(sigma, 3),
                "amp_counts": round(amp, 2),
                "label": label,
                "g_curve": g_pure,
            })

        # v1.2.19 — single Δ_cluster (preserves library spacing).
        # Skip override on phantom clusters (all amps ≈ 0).
        _max_amp_sigma = max(
            (float(s[1]) * float(s[2] or 1.0) for s in (_rebuild_specs or ())),
            default=0.0,
        )
        if (_rebuild_specs and comp_list
                and len(_rebuild_specs) == len(comp_list)
                and _max_amp_sigma > 1.0):
            lib_pos = [float(c["energy_keV"]) for c in comp_list]
            strongest = max(
                range(len(_rebuild_specs)),
                key=lambda k: float(_rebuild_specs[k][1]) * float(_rebuild_specs[k][2] or 1.0),
            )
            delta_cluster = _compute_cluster_global_shift(
                lib_pos[strongest], float(_rebuild_specs[strongest][2] or 1.0),
                E_arr_raw, cont_arr_raw, counts_arr, spec,
            )
            # F-445 / v1.30.3 — when E-cal moved the data, Δ_cluster
            # measured here is residual (typ. ≤ 0.1 keV). Skip rebuild
            # in that case so render and data stay aligned. Render
            # override only fires when F-445 Phase C did NOT converge
            # (rare). |Δ| ≥ 0.1 keV → real residual → still adjust.
            if delta_cluster is not None and abs(delta_cluster) >= 0.1:
                _rebuild_specs = [
                    (lib_pos[k] + delta_cluster,
                     _rebuild_specs[k][1], _rebuild_specs[k][2])
                    for k in range(len(_rebuild_specs))
                ]
                for ci in range(len(comp_list)):
                    comp_list[ci]["fit_centroid_keV"] = round(lib_pos[ci] + delta_cluster, 3)
                    comp_list[ci]["centroid_shift_keV"] = round(delta_cluster, 3)
                shifts_kev = [delta_cluster] * len(_rebuild_specs)

        # v1.2.18 helper: redraw total + per-comp g_curve at shifted centroids.
        total_arr, comp_list = _rebuild_overlay_on_fitted_centroids(
            E_arr_raw, cont_arr_raw, shifts_kev, _rebuild_specs,

            comp_list, fallback_total=getattr(d, "overlay_total", ()) or (),
        )
        multiplet_continua_out.append({
            "cluster_id": str(getattr(d, "cluster_id", "") or ""),
            "E_keV": E_arr,
            "continuum": cont_arr,
            "total": total_arr,
            "components": comp_list,
        })

    # Global composite (v1.2.26): activity-scaled Gaussians + continuum.
    _gc = _build_global_composite(peaks_out, multiplet_continua_out, spec, counts_arr)
    return {
        "peaks": peaks_out,
        "multiplet_continua": multiplet_continua_out,
        "global_composite": _gc,
    }


def _build_spectrum_qc_block(result) -> Optional[Dict[str, Any]]:
    """F-QC-01 / v1.21.0 — build unified spectrum_qc block.

    Returns None if aggregator fails (never blocks report assembly).
    Backward compat: result.bg_quality_check remains as-is in pipeline;
    its data is re-exposed inside the spectrum_qc.peak_z_roi sub-block.
    Cite: KNOWN_AND_FIXED_ISSUES.md:1292, RAG-041, PLAN_v1_20_to_v1_21.md §P0-5.
    """
    try:
        return _build_spectrum_qc(result)
    except Exception:
        # Never fail the report assembly on QC aggregator errors.
        # Fallback: return legacy bg_quality_check shape when available,
        # wrapped in a minimal spectrum_qc envelope for backward compat.
        bqc = getattr(result, "bg_quality_check", None)
        return bqc


def _build_diagnostics(result) -> Dict[str, Any]:
    """Section 12 — diagnostics block."""
    cmp = result.completeness
    spec = result.spec
    real = float(getattr(spec, "real_time", 0.0) or 0.0)
    live = float(getattr(spec, "live_time", 0.0) or 0.0)
    dead_pct = (100.0 * (1.0 - live / real)) if real > 0 and live > 0 else 0.0
    n_xrf = sum(1 for rc in result.residual_classifications if rc.label == "xrf")
    n_ann = sum(1 for rc in result.residual_classifications if rc.label == "annihilation_511")
    n_esc = sum(1 for rc in result.residual_classifications
                if rc.label in ("single_escape", "double_escape"))
    n_sum = sum(1 for rc in result.residual_classifications if rc.label == "sum_peak")

    cascade_nuclides = []
    if result.activities:
        for ar in result.activities:
            if ar.cascade_warning:
                cascade_nuclides.append(ar.nuclide)

    # F-88 chain dominance + K-40 overlap warning surfacing.
    # F-89d adds suppressed_chains / suppression_reason / chain_filtered_out.
    cd = result.chain_dominance
    chain_dominance_block = None
    if cd is not None:
        chain_dominance_block = {
            "th232_dominant": bool(cd.th232),
            "u238_dominant": bool(cd.u238),
            "th232_evidence": list(cd.th232_evidence),
            "u238_evidence": list(cd.u238_evidence),
            "th232_strength_sigma": _safe_float(cd.th232_strength_sigma),
            "u238_strength_sigma": _safe_float(cd.u238_strength_sigma),
            "reason": cd.reason,
            # F-89d
            "suppressed_chains": list(cd.suppressed_chains),
            "suppression_reason": cd.suppression_reason or "",
            "chain_filtered_out_nuclides": list(
                getattr(result, "chain_filtered_out", []) or []
            ),
        }

    # F-119 / v1.17.5 — Chain-equilibrium guard
    chain_equilibrium_block = {}
    try:
        from gamma.activity.compute import chain_equilibrium_guard
        if result.activities:
            chain_equilibrium_block = chain_equilibrium_guard(result.activities)
    except Exception:
        chain_equilibrium_block = {}

    return {
        "measurement_environment": classify_environment(result),
        "dead_time_pct": _safe_float(dead_pct),
        "dead_time_correction_applied": False,  # K-NN: see roadmap
        "pile_up_indicator": dead_pct > 5.0,
        "tcs_correction_applied": bool(result.activities),
        "cascade_warning_nuclides": cascade_nuclides,
        "annihilation_511_observed": n_ann > 0,
        "n_escape_peaks": n_esc,
        "n_sum_peaks": n_sum,
        "n_xrf_residuals": n_xrf,
        "background_subtracted": result.background_subtraction is not None,
        "intrinsic_activity_signature": _intrinsic_signature(result),
        "calibration_quality":
            getattr(result.seven_line_check, "quality", "")
            if result.seven_line_check is not None else "",
        "completeness_dc_pct":
            _safe_float(getattr(cmp, "dc_percent", None)) if cmp is not None else None,
        "completeness_flag":
            getattr(cmp, "flag", None) if cmp is not None else None,
        "fwhm_model_source": result.fwhm_model_source,
        "efficiency_source": result.efficiency_source,
        "efficiency_loaded": result.efficiency_curve is not None,
        # F-88 v1.15.1 — user-priority anchor verdicts + chain dominance
        "chain_dominance": chain_dominance_block,
        "k40_ac228_overlap_warning": bool(result.k40_ac228_overlap_warning),
        # F-119 v1.17.5 — chain-equilibrium diagnostic
        "chain_equilibrium": chain_equilibrium_block,
        # F-129 / v1.17.7 — peak search method dispatch + сравнение
        "peak_search_method": getattr(result, "peak_search_method", "mariscotti"),
        "peak_search_method_comparison":
            _serialize_peak_search_comparison(
                getattr(result, "peak_search_method_comparison", None)
            ),
        # F-131 / v1.17.7 — auto-background search diagnostic block
        "auto_background_search": {
            "mode": getattr(result, "auto_background_mode", "off"),
            "candidates": getattr(result, "auto_background_candidates", None),
            "applied_path": getattr(result, "auto_background_applied_path", None),
        },
        # F-QC-01 / v1.21.0 — unified 6-criterion spectrum_qc block (RAG-041).
        # Aggregates: energy_drift, fwhm_stability, efficiency_qa, bg_drift,
        # per-peak ROI z (BUG-35 pass-through), sensitivity (placeholder).
        # overall_passed = AND of all available criteria.
        # Backward compat: result.bg_quality_check field still populated by
        # staged_pipeline.py:1170-1179 and re-exposed here in peak_z_roi.
        # Cite: spectrum_qc_methodology_v2_2026-06-03.md, KNOWN_AND_FIXED_ISSUES.md:1292.
        "spectrum_qc": _build_spectrum_qc_block(result),
    }


def _serialize_peak_search_comparison(cmp: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """F-129 / v1.17.7 — отфильтровать «сырые» FoundPeak объекты из
    результата compare_peak_methods так, чтобы JSON-сериализатор не падал
    на нестандартных типах. Сохраняем только метрики и счётчики.
    """
    if not cmp:
        return None
    return {
        "primary_method": cmp.get("primary_method", "mariscotti"),
        "secondary_method": cmp.get("secondary_method", "convolution"),
        "n_mariscotti": int(cmp.get("n_mariscotti", cmp.get("n_a", 0))),
        "n_convolution": int(cmp.get("n_convolution", cmp.get("n_b", 0))),
        "n_agreed": int(len(cmp.get("agreed", []) or [])),
        "n_only_mariscotti": int(len(cmp.get("a_only", []) or [])),
        "n_only_convolution": int(len(cmp.get("b_only", []) or [])),
        "agreement_fraction": _safe_float(cmp.get("agreement_fraction", 0.0)),
        "mean_residual_channels": _safe_float(cmp.get("mean_residual_channels", 0.0)),
        "tolerance_channels": _safe_float(cmp.get("tolerance_channels", 1.5)),
    }


def _build_priority_findings(result) -> List[Dict[str, Any]]:
    """F-88 / v1.15.1 — surface the user-priority anchor verdicts."""
    out = []
    for pf in (result.priority_findings or []):
        sig = pf.signal
        out.append({
            "order": sig.order,
            "label": sig.label,
            "nuclide_or_chain": sig.nuclide_or_chain,
            "chain": sig.chain,
            "required_lines_keV": list(sig.required_lines_keV),
            "minimum_required": sig.minimum_required,
            "matched": bool(pf.matched),
            "matched_lines_keV": list(pf.matched_lines_keV),
            "missing_lines_keV": list(pf.missing_lines_keV),
            "max_significance_sigma": _safe_float(pf.significance),
            "rationale": sig.note,
            "note": pf.note,
        })
    return out


def _intrinsic_signature(result) -> Dict[str, Any]:
    """Detector intrinsic-activity signature for the Gamma-1S complex.

    Per `detectors/Gamma-1S/references/05_intrinsic_detector_activity.md`:
    NaI(Tl) has negligible intrinsic activity (trace ⁴⁰K from Na handling
    < 0.01 Bq/cm³). The only routinely-observable artefact is iodine
    K-escape peaks at E_γ − 28.6 keV (Kα) and E_γ − 32.3 keV (Kβ) for
    every strong γ above the I K-binding (33.17 keV). No ¹³⁸La, no
    ²²⁷Ac (those are LaBr₃/CeBr₃ signatures).
    """
    detector = (result.detector_canonical or result.detector_type or "").lower()
    is_gamma1s = "gamma-1c" in detector or detector == "nai"
    if not is_gamma1s:
        return {
            "detector_canonical": result.detector_canonical or result.detector_type,
            "signature": "not characterised for this detector",
        }
    return {
        "detector_canonical": result.detector_canonical or "Gamma-1S",
        "Bq_per_cm3": None,   # negligible per reference 05
        "expected_artefacts": [
            {
                "kind": "I_K_escape_Ka",
                "rule": "E_gamma - 28.6 keV",
                "note": "low-E daughter peak for every γ > ~50 keV",
            },
            {
                "kind": "I_K_escape_Kb",
                "rule": "E_gamma - 32.3 keV",
                "note": "weaker than Kα escape; same population mechanism",
            },
        ],
        "absent_signatures": ["La-138", "Ac-227 series", "Ge K-escape",
                              "Cd K-escape", "Te K-escape"],
    }


def _build_warnings(result) -> List[Any]:
    """Build the operator-facing ``warnings`` list.

    Returns a heterogeneous list of either:
      * ``str`` — bilingual EN/RU operator-text warnings (legacy form);
      * ``dict`` — structured warnings with a machine-readable ``code``
        field, used by anti-hallucination consumers (BUG-40) and by
        F-386-compliant RU renderers.

    Dict entries currently emitted:
      * ``DETECTOR_CYRILLIC_LATIN_FALLBACK`` (BUG-40) — Cyrillic detector
        header (e.g. ``"Гамма-1С"``) silently mapped to a Latin canonical
        whose profile triggers a fallback. Fields: ``code``, ``message``
        (EN), ``original_detector``, ``mapped_to``, ``severity``.

    Downstream consumers (markdown_report.py, html_report.py,
    interactive_html.py, chat_summary.py, cli.py) MUST handle both
    ``str`` and ``dict`` shapes. RU-only renderers MUST NOT render the
    English ``message`` field verbatim — render a localised equivalent
    from the structured fields instead.
    """
    warnings: List[Any] = []
    if result.efficiency_curve is None and (
        result.activities is not None or result.mda_per_line is not None
    ):
        warnings.append(
            "Efficiency curve not loaded — Bq / MDA suite may be incomplete."
        )
    if result.activities:
        for ar in result.activities:
            if ar.cascade_warning:
                warnings.append(
                    f"Cascade summing warning for {ar.nuclide}: {ar.cascade_warning}"
                )
    if result.next_stage_recommended and result.next_stage_reason:
        warnings.append(
            f"Stage-{result.next_stage_recommended} escalation suggested: "
            f"{result.next_stage_reason}"
        )
    # BUG-39 / BUG-40 — surface silent detector-profile fallback.
    fb = getattr(result, "detector_fallback", None)
    if fb and isinstance(fb, dict):
        reason = fb.get("reason", "")
        if reason and reason != "profile_loaded_no_fallback":
            human = fb.get("human") or (
                f"Detector profile fallback: requested "
                f"{fb.get('requested', '?')!r}, using {fb.get('actual', '?')!r} "
                f"(reason: {reason})."
            )
            warnings.append(human)
            # BUG-40 — additionally emit a structured dict warning when the
            # winning canonicalization came from a Cyrillic raw string and
            # the resolved canonical is pure ASCII (homoglyph substitution).
            # Anti-hallucination consumers key off ``code`` rather than
            # English text. Per KFI:1401-1422 brief: severity = MEDIUM.
            if fb.get("cyrillic_to_latin_collision"):
                original_name = fb.get("original_raw") or fb.get("requested", "?")
                canonical = fb.get("requested") or fb.get("actual", "?")
                warnings.append({
                    "code": "DETECTOR_CYRILLIC_LATIN_FALLBACK",
                    "message": (
                        f"Detector name '{original_name}' (Cyrillic) mapped to "
                        f"'{canonical}' (Latin); profile '{canonical}' loaded. "
                        f"Original Cyrillic detector profile not registered. "
                        f"Activity results may reflect efficiency mismatch."
                    ),
                    "original_detector": original_name,
                    "mapped_to": canonical,
                    "severity": "MEDIUM",
                })
        # T41 (BUG-40 (b) hardening) — silent content-fallback class.
        # The path-level cyrillic_to_latin_collision predicate (above)
        # catches name-level homoglyph fallback only. The detector
        # directory may be NAMED correctly while its .efr's
        # `[detector;…]` header records a different physical instance
        # (serial-year mismatch). Surface that as a separate structured
        # warning regardless of the surrounding profile_fallback reason.
        em = fb.get("efficiency_detector_mismatch") if isinstance(fb, dict) else None
        if em and isinstance(em, dict):
            # F-115 anonymization: surface only the structured
            # serial/year integers, NOT the full detector header
            # strings (they encode the certified-source S/N and the
            # exact instrument model — both are PII per F-115). The
            # source-of-truth detector strings stay in the staged_pipeline
            # logger.warning line, which is operator-only and never
            # written into the report artefacts.
            exp_sy = em.get("expected_serial_year") or []
            act_sy = em.get("actual_serial_year") or []
            exp_sn = exp_sy[0] if len(exp_sy) >= 1 else "?"
            exp_yy = exp_sy[1] if len(exp_sy) >= 2 else "?"
            act_sn = act_sy[0] if len(act_sy) >= 1 else "?"
            act_yy = act_sy[1] if len(act_sy) >= 2 else "?"
            warnings.append({
                "code": "EFFICIENCY_DETECTOR_SERIAL_MISMATCH",
                "message": (
                    f"Efficiency curve detector serial mismatch: spectrum "
                    f"detector serial={exp_sn} year=20{exp_yy} vs loaded "
                    f".efr detector serial={act_sn} year=20{act_yy}. "
                    f"Efficiency curve belongs to a different physical "
                    f"instrument; activity results may be biased."
                ),
                "expected_serial_year": list(exp_sy) if exp_sy else None,
                "actual_serial_year": list(act_sy) if act_sy else None,
                "severity": "HIGH",
            })
    return warnings


# ---------------------------------------------------------------------------
# F-070 W3 — Visual similarity block (geometry inference + scoring)
# ---------------------------------------------------------------------------

# Geometry inference table: LSRM COMMENT GEOMETRY= value → canonical class.
# Matched by substring (case-insensitive).  Order matters: more specific first.
_GEOMETRY_SUBSTR_MAP: List[tuple] = [
    ("точ.5см",      "pointlike_5cm"),
    ("точечная-5см", "pointlike_5cm"),
    ("маринелли",    "marinelli_0cm"),
    ("marinelli",    "marinelli_0cm"),
    ("дента-120",    "denta_120ml"),
    ("дента-100",    "denta_100ml"),
    ("петри-60",     "petri_60ml"),
    ("чашка петри",  "petri_60ml"),
    ("petri",        "petri_60ml"),
]


def _infer_geometry_from_comment(comment: str) -> Optional[str]:
    """Return canonical geometry_class from LSRM COMMENT string, or None.

    Scans for 'GEOMETRY=<value>' token; if found maps to canonical class
    using _GEOMETRY_SUBSTR_MAP (substring, case-insensitive).  Returns None
    when GEOMETRY token is absent or value is unrecognised.

    Cite: F-070 W3 brief §Geometry detection (5 classes + None fallback).
    """
    if not comment:
        return None
    import re as _re
    m = _re.search(r"GEOMETRY\s*=\s*([^\r\n;,]+)", comment, _re.IGNORECASE)
    if not m:
        return None
    geo_val = m.group(1).strip().lower()
    for substr, geom_class in _GEOMETRY_SUBSTR_MAP:
        if substr in geo_val:
            return geom_class
    return None


def _build_visual_similarity(result) -> Dict[str, Any]:
    """F-070 W3 / v1.24.0 — visual_similarity JSON block.

    Calls A's score_against_templates() and formats the top-3 matches into
    the JSON shape specified in the W3 brief.

    Returns `enabled=False` with a `reason` key when:
    - background-only spectrum (consistent with bg-aware UI)
    - templates index unavailable (test/CI without audit/_rag/)
    - no geometry match AND best cosine_adjusted < THRESHOLD_AMBIGUOUS_LOWER

    Cite: _state/agent_b/inbox/2026-06-04_F-070-W3_html_card_json_wiring.md §JSON report block shape
    """
    # Guard: background-only spectra don't get visual similarity.
    from gamma.reporting.environment import classify_environment
    env = classify_environment(result)
    if env == "background_only":
        return {"enabled": False, "reason": "background_only_spectrum"}

    # Lazy import A's API — fails gracefully if module absent.
    # Module is importable as `rag.visual_similarity` when scripts/ is on sys.path
    # (conftest.py adds scripts/ to sys.path for tests; run_skill.py does the same).
    try:
        from rag.visual_similarity import (
            compute_query_vector,
            load_templates,
            score_against_templates,
            THRESHOLD_MATCH,
            THRESHOLD_AMBIGUOUS_LOWER,
            TIER_C_DOWNWEIGHT,
            STALE_REFERENCE_AGE_YEARS,
            FEATURE_VECTOR_DIM,
        )
    except ImportError:
        return {"enabled": False, "reason": "module_unavailable"}

    # Geometry from LSRM COMMENT block.
    spec = result.spec
    comment = getattr(spec, "comment", "") or ""
    query_geometry = _infer_geometry_from_comment(comment)

    # Energy calibration for compute_query_vector.
    e_cal = list(getattr(spec, "energy_cal", ()) or ())
    if len(e_cal) >= 2:
        energy_calib: Dict[str, Any] = {"coefficients": e_cal}
    else:
        slope = e_cal[1] if len(e_cal) > 1 else 1.0
        offset = e_cal[0] if len(e_cal) > 0 else 0.0
        energy_calib = {"slope_keV_per_ch": slope, "offset_keV": offset}

    # Load templates (filtered by geometry or all 24 if unknown).
    try:
        templates = load_templates(geometry_class=query_geometry)
    except Exception:
        return {"enabled": False, "reason": "templates_unavailable"}

    if not templates:
        return {"enabled": False, "reason": "templates_unavailable"}

    # Build query vector from spectrum counts.
    counts = getattr(spec, "counts", None)
    if counts is None:
        return {"enabled": False, "reason": "no_counts"}

    try:
        import numpy as _np_vs
        q = compute_query_vector(_np_vs.asarray(counts), energy_calib)
    except Exception:
        return {"enabled": False, "reason": "vector_encoding_failed"}

    # Score and return top-3.
    top_k = 3
    matches_raw = score_against_templates(q, templates, top_k=top_k)

    # Disabled fallback: unknown geometry AND best adjusted cosine < mismatch threshold.
    if (query_geometry is None
            and matches_raw
            and matches_raw[0]["cosine_adjusted"] < THRESHOLD_AMBIGUOUS_LOWER):
        return {
            "enabled": False,
            "reason": "no_match_above_mismatch_threshold",
        }

    matches_out = []
    for m in matches_raw:
        cosine_raw = m["cosine_raw"]
        cosine_adj = m["cosine_adjusted"]
        decay_age = m["decay_age_years"]
        matches_out.append({
            "template_id": m["template_id"],
            "nuclide": m["nuclide"],
            "geometry_class": m["geometry_class"],
            "tier": m["tier"],
            "cosine_raw": round(float(cosine_raw), 4) if cosine_raw is not None else None,
            "cosine_adjusted": round(float(cosine_adj), 4) if cosine_adj is not None else None,
            "verdict": m["verdict"],
            "decay_age_years": round(float(decay_age), 1) if decay_age is not None else None,
            "stale_reference": m["stale_reference"],
            "cert_reference_dates": list(m.get("cert_reference_dates") or []),
        })

    # verdict_summary: best match verdict and nuclide.
    best_verdict = matches_out[0]["verdict"] if matches_out else "mismatch"
    best_nuclide = (
        matches_out[0]["nuclide"]
        if matches_out and best_verdict != "mismatch"
        else None
    )

    return {
        "enabled": True,
        "policy": {
            "threshold_match": THRESHOLD_MATCH,
            "threshold_ambiguous_lower": THRESHOLD_AMBIGUOUS_LOWER,
            "tier_c_downweight": TIER_C_DOWNWEIGHT,
            "stale_reference_age_years": STALE_REFERENCE_AGE_YEARS,
        },
        "query_geometry": query_geometry,
        "query_vector_dim": FEATURE_VECTOR_DIM,
        "top_k": top_k,
        "matches": matches_out,
        "verdict_summary": best_verdict,
        "verdict_summary_nuclide": best_nuclide,
    }


def _inject_compton_edge_se_for_confirmed(
    secondary_peaks: List[Dict[str, Any]],
    result,
) -> List[Dict[str, Any]]:
    """F-141 / v1.17.7 — инжекция Compton edge / Backscatter / SE / DE
    ТОЛЬКО для выбранных интенсивных photopeak'ов с выраженным
    Compton-континуумом.

    Подтверждение даётся не для каждой confirmed линии — это создавало бы
    шум, — а только для **strong-line** нуклидов с физически выраженным
    плечом континуума:

        Cs-137  662   — выраженное плато, обратное рассеяние ~184 кэВ
        K-40    1461  — широкий Compton continuum, edge ~1244
        Tl-208  2614  — самый интенсивный плечо на NaI; SE=2103, DE=1592
        Co-60   1173, 1332 — двойной плечо; backscatter ~210
        Bi-214  1764  — заметный Compton continuum

    Дополнительное требование: peak_area соответствующей линии должна
    превышать порог `PEAK_AREA_THRESHOLD` — слабые линии не дают
    наблюдаемого Compton continuum.

    Формулы (Knoll 4ed §10, Gilmore §7):
      • Compton edge: E_CE = 2·E_γ² / (511 + 2·E_γ)
      • Backscatter:  E_BS = E_γ / (1 + 2·E_γ/511)
      • Single escape: E_SE = E_γ − 511 (только при E_γ > 1022)
      • Double escape: E_DE = E_γ − 1022 (только при E_γ > 2044)
    """
    if not secondary_peaks:
        secondary_peaks = []

    # Каталог «strong-line» нуклидов с обязательным Compton-anchor.
    # Key: (nuclide, E_γ_keV). Value: набор features, которые
    # инжектируются при превышении peak_area порога.
    STRONG_LINE_ANCHORS = {
        ("Cs-137", 661.66): ("compton_edge", "backscatter"),
        ("K-40",   1460.82): ("compton_edge", "backscatter"),
        ("Tl-208", 2614.51): ("compton_edge", "backscatter",
                              "single_escape", "double_escape"),
        ("Co-60",  1173.23): ("compton_edge", "backscatter"),
        ("Co-60",  1332.49): ("compton_edge", "backscatter",
                              "single_escape"),
        ("Bi-214", 1764.49): ("compton_edge", "backscatter",
                              "single_escape"),
        ("Bi-214", 1120.29): ("compton_edge",),
        ("Bi-214", 609.31):  ("compton_edge", "backscatter"),
        ("Ac-228", 911.20):  ("compton_edge", "backscatter"),
        ("Ac-228", 968.97):  ("compton_edge",),
        ("Pb-214", 351.93):  ("compton_edge",),
    }
    PEAK_AREA_THRESHOLD = 1000.0  # counts — порог «выраженного» photopeak

    # Существующие entries (избежать дублирования по feature+parent)
    existing = set()
    for sp in secondary_peaks:
        fk = (sp.get("feature_kind") or "").lower()
        parent = sp.get("parent_nuclide") or ""
        E_sp = sp.get("energy_keV") or 0.0
        # Дубли по (feature_kind, parent_nuclide) или (feature_kind, E)
        existing.add((fk, parent))
        existing.add((fk, round(E_sp, 0)))

    # Собрать confirmed peak_area по (nuc, library_E_keV)
    confirmed_areas: Dict[Tuple[str, float], float] = {}
    for ni in (result.final_detected or []):
        for m in (ni.matched_lines or ()):
            E_lib = float(getattr(m, "library_E_keV", 0.0))
            pa = float(getattr(m, "peak_area", 0.0) or 0.0)
            for (nuc, E_anchor), _ in STRONG_LINE_ANCHORS.items():
                if ni.nuclide == nuc and abs(E_lib - E_anchor) < 1.0:
                    confirmed_areas[(nuc, E_anchor)] = max(
                        confirmed_areas.get((nuc, E_anchor), 0.0), pa
                    )

    additions = []
    for (nuc, E_g), features in STRONG_LINE_ANCHORS.items():
        pa = confirmed_areas.get((nuc, E_g), 0.0)
        if pa < PEAK_AREA_THRESHOLD:
            continue   # слабый photopeak — Compton plateau не виден
        for feat in features:
            if feat == "compton_edge":
                E_x = 2.0 * E_g * E_g / (511.0 + 2.0 * E_g)
                note = (f"Комптоновский край фотопика {nuc} {E_g:.1f} кэВ "
                        f"(Knoll §10: E_CE={E_x:.1f} кэВ); "
                        f"plateau физически присутствует "
                        f"(площадь ФЭП = {pa:.0f} > {PEAK_AREA_THRESHOLD:.0f}).")
            elif feat == "backscatter":
                E_x = E_g / (1.0 + 2.0 * E_g / 511.0)
                note = (f"Обратнорассеянное излучение от {nuc} {E_g:.1f} кэВ "
                        f"(Gilmore §7: E_BS={E_x:.1f} кэВ).")
            elif feat == "single_escape":
                if E_g <= 1022.0:
                    continue
                E_x = E_g - 511.0
                note = (f"Single escape от {nuc} {E_g:.1f} кэВ "
                        f"(E_SE = E_γ − 511 = {E_x:.1f} кэВ).")
            elif feat == "double_escape":
                if E_g <= 2044.0:
                    continue
                E_x = E_g - 1022.0
                note = (f"Double escape от {nuc} {E_g:.1f} кэВ "
                        f"(E_DE = E_γ − 1022 = {E_x:.1f} кэВ).")
            else:
                continue
            # Дубль если уже есть (feat, parent_nuclide) ИЛИ (feat, E_x)
            if (feat, nuc) in existing or (feat, round(E_x, 0)) in existing:
                continue
            additions.append({
                "channel": None,
                "energy_keV": round(E_x, 1),
                "significance": None,
                "type": "расчётная_сигнатура",
                "feature_kind": feat,
                "parent_nuclide": nuc,
                "parent_line_keV": E_g,
                "note": note,
            })
            existing.add((feat, nuc))
            existing.add((feat, round(E_x, 0)))
    return secondary_peaks + additions


def _filter_and_augment_secondary_peaks(
    secondary_peaks: List[Dict[str, Any]],
    diagnostics: Dict[str, Any],
) -> List[Dict[str, Any]]:
    """F-110 (D-08, D-09).

    * Drop diffuse zones (``backscatter_region`` / ``broad_compton_plateau``)
      that the user explicitly objected to (Зоны 200-280, ~510).
    * If Th-232 chain is dominant, prepend the mandatory 73-90 кэВ
      composite cluster entry.
    """
    # D-09 — drop diffuse zone artefacts
    DROP_TYPES = {"backscatter_region", "broad_compton_plateau"}
    filtered = [sp for sp in (secondary_peaks or [])
                if (sp.get("type") or "") not in DROP_TYPES]

    cd = (diagnostics or {}).get("chain_dominance") or {}
    if cd.get("th232_dominant"):
        # D-08 — mandatory 73-90 кэВ composite when Th-232 chain dominant
        composite_entry = {
            "channel": None,
            "energy_keV": 80.0,
            "significance": None,
            "type": "cluster",
            "feature_kind": "composite_cluster",
            "parent_nuclide": "Th-232",
            "parent_line_keV": None,
            "note": (
                "Композит 73-90 кэВ: Pb K-РИ от внутренней конверсии "
                "Pb-212 239 (α_K≈1.4) + ВК Tl-208/Bi-212 + Th-228 84.37 "
                "+ ≪1% флуоресценция Pb-50 экрана. "
                "F-110: не разлагается."
            ),
        }
        # Insert at top, but only if not already present
        already = any(
            (sp.get("type") == "cluster"
             and 73.0 <= (sp.get("energy_keV") or 0.0) <= 90.0)
            for sp in filtered
        )
        if not already:
            filtered = [composite_entry] + filtered
    return filtered


# F-440 / v1.30.0 -- Two-phase weak-line completion (Phase 2).
# Wire complete_weak_lines() into the JSON report: rebuild a per-nuclide
# matched_lines dict from result.final_detected, gather multiplet ROIs
# from result.deconvolution_results, and serialize via to_json_block().
def _build_weak_line_completion(result):
    """Phase 2 weak-line completion block. None when prerequisites absent."""
    try:
        from gamma.activity.weak_line_completion import (
            complete_weak_lines,
            to_json_block,
            DEFAULT_MIN_GROUPING_SNR,
            DEFAULT_MIN_GROUPING_INTENSITY_PCT,
        )
    except ImportError:
        return None
    activities = getattr(result, 'activities', None)
    eff_curve = getattr(result, 'efficiency_curve', None)
    spec = getattr(result, 'spec', None)
    t_live = float(getattr(spec, 'live_time', 0.0) or 0.0)
    if activities is None or eff_curve is None or t_live <= 0:
        return None
    final_detected = getattr(result, 'final_detected', []) or []
    matched_by_nuc = {}
    for ni in final_detected:
        nuc = getattr(ni, 'nuclide', None)
        if not nuc:
            continue
        ml = getattr(ni, 'matched_lines', ())
        matched_by_nuc[nuc] = list(ml)
    if not matched_by_nuc:
        return None
    # DEBUG F-440: dump live matched_lines for diagnostics
    import os
    if os.environ.get('F440_DEBUG'):
        import sys
        for nuc, ml in matched_by_nuc.items():
            if nuc not in ('Pb-212','Ac-228'):
                continue
            print(f'[F440-DEBUG] {nuc}: {len(ml)} matched lines', file=sys.stderr)
            for m in ml[:20]:
                print(f'  E={getattr(m,"library_E_keV",None)} I={getattr(m,"library_I_pct",None)} area={getattr(m,"peak_area",None)} src={getattr(m,"peak_area_source",None)}', file=sys.stderr)
    library = {}
    try:
        from gamma.data.nuclide_library import get_nuclide
        for nuc in matched_by_nuc.keys():
            try:
                entry = get_nuclide(nuc)
                if entry is None:
                    continue
                if hasattr(entry, 'to_dict'):
                    library[nuc] = entry.to_dict()
                elif isinstance(entry, dict):
                    library[nuc] = entry
                else:
                    lines = getattr(entry, 'lines', None)
                    if lines is not None:
                        library[nuc] = {'lines': list(lines)}
            except Exception:
                continue
    except ImportError:
        return None
    deconv = getattr(result, 'deconvolution_results', None) or []
    multiplet_rois = []
    for idx, d in enumerate(deconv):
        try:
            ch_lo = int(getattr(d, 'roi_low_ch', 0) or 0)
            ch_hi = int(getattr(d, 'roi_high_ch', 0) or 0)
            if ch_hi <= ch_lo:
                continue
            E_lo = float(spec.channel_to_energy(ch_lo))
            E_hi = float(spec.channel_to_energy(ch_hi))
            label = str(getattr(d, 'cluster_id', None) or f'M{idx + 1}')
            if E_lo > 0 and E_hi > E_lo:
                multiplet_rois.append({'label': label, 'E_lo_keV': E_lo, 'E_hi_keV': E_hi})
        except (TypeError, ValueError, AttributeError):
            continue
    spec_acts = getattr(result, 'specific_activities_Bq_per_kg', None)
    sample_mass_kg = getattr(result, 'sample_mass_kg', None)
    completions = complete_weak_lines(
        activities=activities,
        matched_lines_by_nuclide=matched_by_nuc,
        nuclide_library=library,
        efficiency_curve=eff_curve,
        t_live=t_live,
        sample_mass_kg=sample_mass_kg,
        specific_activities_Bq_per_kg=spec_acts,
        self_attenuation_factors=None,
        tcs_factors=None,
        # F-440 Phase 2 classification thresholds are FIXED at 5.0/3.0
        # regardless of grouping defaults (which can be 0.0 in pipeline).
        min_grouping_snr=5.0,
        min_grouping_intensity_pct=3.0,
        multiplet_rois=multiplet_rois if multiplet_rois else None,
        contamination_threshold_pct=1.0,
    )
    if not completions:
        return None
    return {
        'min_grouping_snr': DEFAULT_MIN_GROUPING_SNR,
        'min_grouping_intensity_pct': DEFAULT_MIN_GROUPING_INTENSITY_PCT,
        'per_nuclide': to_json_block(completions),
    }


# G2 / v1.31.2 -- intensity-ratio chi^2 gate emit builder.
def _build_intensity_ratio_chi2_gate(result):
    try:
        from gamma.identification.cross_check import (
            intensity_ratio_chi2_gate,
            INTENSITY_RATIO_CHI2_STRICT_THRESHOLD as STRICT_T,
            INTENSITY_RATIO_CHI2_LENIENT_THRESHOLD as LENIENT_T,
        )
    except ImportError:
        return None
    ident = getattr(result, "identification_result", None)
    if ident is None:
        return None
    per_nuc = intensity_ratio_chi2_gate(ident)
    if not per_nuc:
        return None
    return {
        "strict_threshold": STRICT_T,
        "lenient_threshold": LENIENT_T,
        "per_nuclide": per_nuc,
    }

def build_json_report(result) -> Dict[str, Any]:
    """Assemble the full JSON report dict from a StagedAnalysisResult.

    The returned dict is fully JSON-serializable (no NaN / inf;
    datetimes as ISO-8601 strings; tuples flattened to lists).

    Schema version: SCHEMA_VERSION (currently "0.6"). Bump on
    breaking field changes.
    """
    # F-110 / D-08, D-09 — adjust secondary peaks before assembly.
    # F-141 / v1.17.7 — обязательная инжекция Compton edge / SE / DE
    # для каждой confirmed photopeak.
    raw_secondary = _build_secondary_peaks(result)
    diagnostics = _build_diagnostics(result)
    secondary = _filter_and_augment_secondary_peaks(raw_secondary, diagnostics)
    secondary = _inject_compton_edge_se_for_confirmed(secondary, result)

    # F-397 / v1.18.27 — pre-build bg peak block (если detection прогонялась
    # на фоне). Используем те же builder helpers, что и для sample, чтобы
    # схема bg-полей повторяла sample-схему byte-for-byte. Пустые списки
    # когда фон не анализировался.
    bg_result = getattr(result, "background_staged_result", None)
    if bg_result is not None:
        bg_primary_feps = _build_primary_feps(bg_result)
        bg_secondary_raw = _build_secondary_peaks(bg_result)
        bg_secondary = _filter_and_augment_secondary_peaks(bg_secondary_raw, {})
        bg_deconvolutions = _build_deconvolutions(bg_result)
    else:
        bg_primary_feps = []
        bg_secondary = []
        bg_deconvolutions = []

    report = {
        "schema_version": SCHEMA_VERSION,
        "skill_version": SKILL_VERSION,
        "header": _build_header(result),
        "calibration": _build_calibration(result),
        # F-88 v1.15.1 — express priority verdicts (user methodology)
        "priority_express_findings": _build_priority_findings(result),
        "primary_feps": _build_primary_feps(result),
        "secondary_peaks": secondary,
        # F-397 / v1.18.27 — фон-only peak block (для HTML toggle «Фон»).
        # Пустые списки когда фон не анализировался.
        "background_primary_feps": bg_primary_feps,
        "background_secondary_peaks": bg_secondary,
        "background_multiplet_deconvolutions": bg_deconvolutions,
        "elemental_xrf": _build_elemental_xrf(result),
        "identified_nuclides": _build_identified_nuclides(result),
        "unidentified_peaks": _build_unidentified_peaks(result),
        "completeness": {
            "dc_pct": _safe_float(
                getattr(result.completeness, "dc_percent", None)
                if result.completeness is not None else None
            ),
            "flag": (getattr(result.completeness, "flag", None)
                     if result.completeness is not None else None),
        },
        "mda": _build_mda(result),
        "multiplet_deconvolutions": _build_deconvolutions(result),
        # F-FIT-VIEW / v1.22.1 — fit overlay data for interactive HTML toggle.
        # Contains per-peak Gaussian params (singlets) + per-cluster continuum +
        # total-fit arrays (multiplets). Used by the «Подгонка» toggle button.
        "fit_overlay": _build_fit_overlay(result),
        # F-070 W3 / v1.24.0 — visual similarity scoring against canonical
        # spectrum templates (24 entries, geometry-filtered). Emitted after
        # fit_overlay, before decision summary. enabled=False when bg-only,
        # templates unavailable, or no match above mismatch threshold.
        "visual_similarity": _build_visual_similarity(result),
        # F-145 / v1.17.8 — two-phase multiplet self-calibration diagnostic.
        # Pass-through из StagedAnalysisResult.multiplet_self_calibration_diag.
        # None — F-145 не запускалась (нет мультиплетов / не NaI).
        "multiplet_self_calibration": getattr(
            result, "multiplet_self_calibration_diag", None
        ),
        "diagnostics": diagnostics,
        "warnings": _build_warnings(result),
        "pipeline_notes": list(result.notes or []),
        # F-440 / v1.30.0 — two-phase weak-line completion. Per-nuclide
        # projection of weak (S/N<5 OR I<3%) library lines from Phase 1
        # activity; per-nuclide completeness; cross-nuclide contamination
        # detection inside multiplet ROIs. None when no efficiency curve.
        "weak_line_completion": _build_weak_line_completion(result),
        # G2 / v1.31.2 — intensity-ratio chi² gate per detected nuclide.
        # Annotation only; thresholds 1.5 (strict) / 3.0 (lenient).
        "intensity_ratio_chi2_gate": _build_intensity_ratio_chi2_gate(result),
    }

    # G4 / v1.31.0 — annotate primary_feps with bg_carryover. fit_overlay.peaks
    # already get the annotation inside _build_fit_overlay; primary_feps is the
    # source the HTML peaks table renders from, so we mark it here using the
    # same catalog.
    _bg_cat_feps = _build_bg_peak_catalog(bg_result)
    if _bg_cat_feps:
        _fm_feps = getattr(result, "fwhm_model", None)
        def _fwhm_feps(E: float) -> float:
            if _fm_feps is None:
                return max(E * 0.07, 1.0)
            # F-452: polymorphic — FwhmModel callable | legacy 3-tuple.
            return _fwhm_keV_at_energy(_fm_feps, float(E))
        _mark_bg_carryover(
            report["primary_feps"], _bg_cat_feps, _fwhm_feps,
            sample_sources=None, energy_field="peak_E_keV",
        )

    # F-115 (D-10) — strip personal / confidential identifiers from the
    # whole report dict before any artefact is written.
    from gamma.reporting.anonymize import anonymize_report_inplace
    anonymize_report_inplace(report)

    return report


__all__ = [
    "SCHEMA_VERSION", "SKILL_VERSION",
    "build_json_report",
]
