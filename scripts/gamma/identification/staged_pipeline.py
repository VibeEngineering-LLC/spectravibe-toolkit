"""
Staged identification pipeline (F-69 / v1.11.1).

Orchestrator that ties together the existing pieces:

    read .spe (LSRM)
       │
       ▼
   build empirical FWHM(E)  ← from `lsrm_peaks_table` when stored
       │                       polynomial form is unusable
       ▼
   Mariscotti peak search
       │
       ▼
   Stage 1 identify_nuclides  ← only ЕРН candidates
       │
       ▼
   disambiguate_identifications  ← Lsrm Rule 1..5
       │
       ▼
   compute unmatched-residual diagnostic
       │
       ▼
  → Stage 2/3 only if user_allowed OR residuals justify it
       │
       ▼
   return StagedAnalysisResult

Per user methodology (15.11.2025):
  • Default is Stage 1 only ("don't fantasize").
  • Stage 2 (Cs-137/Cs-134/Co-60/I-131) needs explicit `allow_stage2=True`
    or an `auto_escalate=True` policy + significant unmatched residuals.
  • Stage 3 (Na-22, Be-7, Am-241, calibration, medical) is always
    opt-in.

The 511 keV peak in Th-232-rich samples is dominated by Tl-208 510.77;
the disambiguate Rule 2 removes Na-22 claims when Tl-208 is confirmed,
so Na-22 false-positives are handled automatically.
"""

from __future__ import annotations

import json
import logging
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple, Callable, Any, Sequence

import numpy as np

logger = logging.getLogger(__name__)

from gamma.io.readers import read_spectrum
from gamma.io.bg_control import (
    validate_background as _f243_validate_background,
    check_bg_quality as _f243_check_bg_quality,
    bg_z_test as _bg_z_test,
)
from gamma.peaks.search import mariscotti_search, FoundPeak
from gamma.identification.identify import (
    identify_nuclides, IdentificationResult, NuclideIdentification,
)
from gamma.identification.disambiguate import disambiguate_identifications
from gamma.identification.window import (
    identification_window_from_fwhm,
    build_id_window_k_fwhm,
)
from gamma.identification.id_window import normalize_detector_class
from gamma.identification.ern_set import (
    candidates_for_stage, stage_of_nuclide,
    ERN_7_REFERENCE_LINES_keV,
)
from gamma.identification.residual_classifier import (
    classify_residuals, ResidualClassification,
    LBL_TRUE_UNMATCHED, LBL_XRF, LBL_CHAIN_SECONDARY,
    LBL_SINGLE_ESCAPE, LBL_DOUBLE_ESCAPE, LBL_SUM_PEAK,
    LBL_ANNIHILATION, LBL_EDGE_OF_RANGE,
)
from gamma.identification.anchor_ranks import (
    find_anchor_matches, confirm_express_patterns,
    seed_calibration_anchors,
    AnchorMatch, PatternConfirmation, AnchorSeedResult,
    # F-88 — user-prioritized express order + chain dominance
    derive_priority_findings, derive_chain_dominance,
    PriorityFinding, ChainDominance,
    TH232_PROXY_NUCLIDES, U238_PROXY_NUCLIDES,
    CALIBRATION_RANK_START,
)
from gamma.calibration.seven_line_check import (
    run_seven_line_check, SevenLineCheck,
)
from gamma.calibration.bg_subtract_energy import (
    subtract_background, BackgroundSubtractionResult,
)
from gamma.calibration.efficiency_autoload import (
    load_efficiency_for_geometry, EFFICIENCY_FIT_FAILED,
)
from gamma.identification.ci_gating import (
    gate_identifications, CIGating,
)
from gamma.identification.completeness import (
    compute_completeness, CompletenessResult,
)
# F-84 / v1.13.0 — Round 5: activities + MDA + multiplet deconvolution
from gamma.activity.compute import (
    compute_activities_for_all, ActivityResult,
)
from gamma.identification.mda import (
    mda_for_peak, MdaResult,
)
from gamma.peaks.deconvolve import (
    apply_multiplet_deconvolution, DeconvolutionResult,
    run_chain_forced_multiplets, TH232_FORCED_CLUSTERS,
)
# F-145 / v1.17.8 — Two-phase multiplet self-calibration
from gamma.calibration.multiplet_self_calibration import (
    recalibrate_from_multiplet_centroids,
    SelfCalibrationDiag,
    PHASE_D_CENTROID_TOLERANCE_FRAC,
)


def _f445_refresh_walk(detected_nuclides, ch2e, dc_replace):
    """Walk detected_nuclides and rebuild peak_E_keV per matched line."""
    new_detected = []
    for ni in detected_nuclides:
        new_ml = []
        for m in ni.matched_lines:
            try:
                ch = float(m.peak_channel)
                E_new = float(ch2e(ch))
                resid_new = abs(E_new - float(m.library_E_keV))
                new_ml.append(dc_replace(
                    m, peak_E_keV=E_new, residual_keV=resid_new,
                ))
            except Exception:
                new_ml.append(m)
        new_detected.append(dc_replace(ni, matched_lines=tuple(new_ml)))
    return new_detected


def _f445_refresh_match_energies(id_result, spec, IR_cls, dc_replace):
    """F-445: recompute LineMatch.peak_E_keV with current spec.energy_cal."""
    if id_result is None or spec is None:
        return id_result
    try:
        ch2e = spec.channel_to_energy
    except Exception:
        return id_result
    new_detected = _f445_refresh_walk(
        id_result.detected_nuclides, ch2e, dc_replace,
    )
    return IR_cls(
        detector_type=id_result.detector_type,
        window=id_result.window,
        candidates_considered=id_result.candidates_considered,
        detected_nuclides=tuple(new_detected),
        rejected_nuclides=id_result.rejected_nuclides,
        unmatched_peaks=id_result.unmatched_peaks,
        notes=id_result.notes,
    )


def _f445_note_cal_accepted(diag, chi2_A, chi2_D) -> None:
    """F-445: log carve-out (spec.cal updated despite Phase D rollback)."""
    if diag is None:
        return
    diag.reason += (
        " | F-445 cluster-Δ: spec.energy_cal обновлён несмотря на "
        "Phase D rollback (chi2_A=" + format(chi2_A, ".2f")
        + " → chi2_D=" + format(chi2_D, ".2f") + ")"
    )


def _f445_note_phase_d_rollback(diag, chi2_A, chi2_D) -> None:
    """Legacy F-145 path: Phase D rollback → keep old cal."""
    if diag is None:
        return
    diag.phase_C_applied = False
    diag.reason += (
        " | фаза Д откат: χ²_sum " + format(chi2_A, ".2f")
        + " → " + format(chi2_D, ".2f") + " (хуже)"
    )


def _f445_is_cluster_delta_used(diag) -> bool:
    """F-445: True if at least one anchor came from cluster-Δ collector."""
    if diag is None:
        return False
    try:
        return any(
            "cluster_delta_" in str(a.get("source", ""))
            for a in (diag.anchors_used or [])
        )
    except Exception:
        return False


def _f445_build_continuum_arrays(forced_clusters):
    """F-445: {cluster_id: (E_arr_list, cont_arr_list)} from overlay arrays."""
    out = {}
    for _fc in forced_clusters or ():
        _cid = str(getattr(_fc, "cluster_id", "") or "")
        _e = getattr(_fc, "overlay_E_keV", None)
        _c = getattr(_fc, "overlay_continuum", None)
        if _cid and _e and _c and len(_e) == len(_c):
            out[_cid] = (list(_e), list(_c))
    return out


def _f453_build_singleton_extras(
    anchor_matches, fwhm_provider_keV, forced_clusters,
):
    """F-453 (BUG-38 follow-up, 2026-06-23) — singleton-anchor fallback для F-145.

    На short NaI fixtures (AmTiCsEu, Cs-Co) `n_multiplets_seen=0` потому что
    нет Th/U `chain_dominance` forced_clusters; F-145 Phase B/C тогда
    замолкает с `reason='ни один мультиплет... не прошёл фильтр'`. Этот
    helper собирает confirmed singleton anchors (Am-241/Cs-137/Eu-152/Ti-44)
    из `anchor_matches` и возвращает их как `extra_anchors` tuple-list для
    `recalibrate_from_multiplet_centroids(extra_anchors=...)`. F-145 Phase C
    тогда видит ≥3 anchor'а через ``extra_anchors`` branch и refit'ит δ(N).

    Фильтр:
      • anchor.nuclide non-empty
      • partner_required_but_missing == False
      • (nuclide, round(E_passport, 1)) NOT в активных forced_clusters
        (по components — чтобы не дублировать линию, которую multiplet
        машина уже учитывает)

    Возвращает: list of (E_passport_keV, channel_obs, source_label).
    Source label: ``F-453_singleton_<nuclide>_<E>kev``.
    """
    in_cluster = set()
    for fc in forced_clusters or ():
        for comp in getattr(fc, "components", []) or ():
            nuc = str(getattr(comp, "nuclide", "") or "")
            E_lib_raw = (
                getattr(comp, "line_E_keV", None)
                if getattr(comp, "line_E_keV", None) is not None
                else getattr(comp, "E_keV", 0.0)
            )
            E_lib = float(E_lib_raw or 0.0)
            if nuc and E_lib > 0:
                in_cluster.add((nuc, round(E_lib, 1)))
    extras = []
    for am in anchor_matches or ():
        anchor = getattr(am, "anchor", None)
        if anchor is None or not getattr(anchor, "nuclide", ""):
            continue
        if getattr(am, "partner_required_but_missing", False):
            continue
        nuc = str(anchor.nuclide)
        E_lib = float(getattr(anchor, "energy_keV", 0.0) or 0.0)
        if E_lib <= 0:
            continue
        if (nuc, round(E_lib, 1)) in in_cluster:
            continue
        ch_obs_raw = getattr(am, "peak_channel", None)
        if ch_obs_raw is None:
            continue
        ch_obs = float(ch_obs_raw)
        extras.append((
            E_lib, ch_obs,
            f"F-453_singleton_{nuc}_{int(round(E_lib))}keV",
        ))
    return extras or None


def _serialize_f145_diag(diag: SelfCalibrationDiag) -> dict:
    """F-145: SelfCalibrationDiag → JSON-safe dict для StagedAnalysisResult."""
    return {
        "attempted": bool(diag.attempted),
        "phase_A_run": bool(diag.phase_A_run),
        "phase_B_passed": bool(diag.phase_B_passed),
        "phase_C_applied": bool(diag.phase_C_applied),
        "n_multiplets_seen": int(diag.n_multiplets_seen),
        "n_multiplets_phase_A_converged": int(
            diag.n_multiplets_phase_A_converged),
        "n_anchors_collected": int(diag.n_anchors_collected),
        "n_anchors_after_filter": int(diag.n_anchors_after_filter),
        "anchors_used": list(diag.anchors_used),
        "old_energy_cal": diag.old_energy_cal,
        "new_energy_cal": diag.new_energy_cal,
        "old_residual_max_keV": diag.old_residual_max_keV,
        "new_residual_max_keV": diag.new_residual_max_keV,
        "degree_used": diag.degree_used,
        "reason": str(diag.reason),
        "phase_A_chi2_per_mult": dict(diag.phase_A_chi2_per_mult),
        # F-446: adaptive Phase C degree policy diagnostics
        "delta_degree_used": getattr(diag, "delta_degree_used", None),
        "delta_const_keV": getattr(diag, "delta_const_keV", None),
        "degree_choice_reason": str(getattr(diag, "degree_choice_reason", "") or ""),
        "accepted_cluster_ids": list(getattr(diag, "accepted_cluster_ids", []) or []),
    }
from gamma.physics.cascade_summing import (
    compute_tcs_corrections, peak_to_total_NaI_for_geometry,
)
from gamma.data.nuclide_library import get_nuclide


# ──────────────────────────────────────────────────────────────────
# F-116 / v1.17.5 — Th-chain and U-chain member sets (out-of-chain
# Stage-3 suppression). Mirrors TH232_PROXY_NUCLIDES / U238_PROXY_NUCLIDES
# from anchor_ranks but is kept here so the suppression logic is local
# and self-contained.
# ──────────────────────────────────────────────────────────────────
_TH_CHAIN_NUCLIDES = frozenset({
    "Th-232", "Tl-208", "Pb-212", "Ac-228", "Bi-212", "Th-228", "Ra-224",
})
_U_CHAIN_NUCLIDES = frozenset({
    "U-238", "Bi-214", "Pb-214", "Pb-210", "Ra-226", "Po-214", "Th-234",
})
# F-130 / v1.19.1 — Primordial nuclides that are ubiquitous in natural
# background but belong to neither the Th-232 nor U-238 decay chains.
# Always kept by F-116 suppression regardless of chain dominance.
_BACKGROUND_NATURAL_NUCLIDES = frozenset({
    "K-40",  # 1460.82 keV; primordial, present in all natural-background spectra
})


def _apply_f116_out_of_chain_suppression(
    final_detected: list,
    chain_dominance_out,
    filename_isotope_hints,
) -> Tuple[list, List[str]]:
    """F-116 / v1.19.1 — Out-of-chain Stage-3 suppression (pure helper).

    Extracted in v1.19.1 (#124) to make the union-of-chains fix unit-testable
    in isolation (see tests/test_f116_dual_chain.py). The inline body inside
    run_staged() calls this helper; behaviour is preserved.

    Behaviour:
      - both chains DOMINANT → chain_set = union(Th, U)
        (FIX #124: previously only Th was used → spurious U-cycle suppression)
      - only Th DOMINANT → chain_set = Th set
      - only U  DOMINANT → chain_set = U  set
      - neither DOMINANT → returns (final_detected unchanged, [])

    Filtering rules (preserved from v1.17.5):
      hints_in_chain_only ⇒ keep ni if nuc ∈ chain_set ∪ hint_set,
      else require σ ≥ 10 AND ≥2 matched lines AND CI ≥ 5 to keep.

    Returns (kept_list, out_of_chain_suppressed_messages).
    """
    out_of_chain_suppressed: List[str] = []
    if chain_dominance_out is None or not (
        chain_dominance_out.th232 or chain_dominance_out.u238
    ):
        return final_detected, out_of_chain_suppressed

    chain_set_parts: List[frozenset] = []
    chain_labels: List[str] = []
    if chain_dominance_out.th232:
        chain_set_parts.append(_TH_CHAIN_NUCLIDES)
        chain_labels.append("Th-232")
    if chain_dominance_out.u238:
        chain_set_parts.append(_U_CHAIN_NUCLIDES)
        chain_labels.append("U-238")
    chain_set = frozenset().union(*chain_set_parts)
    chain_label = "+".join(chain_labels)
    hint_set = set(filename_isotope_hints or ())

    hints_in_chain_only = (not hint_set) or hint_set.issubset(chain_set)
    if not hints_in_chain_only:
        return final_detected, out_of_chain_suppressed

    kept = []
    for ni in final_detected:
        nuc = ni.nuclide
        if nuc in chain_set or nuc in hint_set or nuc in _BACKGROUND_NATURAL_NUCLIDES:
            kept.append(ni)
            continue
        n_matched = len(ni.matched_lines)
        best_sigma = 0.0
        for m in ni.matched_lines:
            area = getattr(m, "peak_area", None)
            area_unc = getattr(m, "peak_area_uncertainty", None)
            if area and area_unc and area_unc > 0:
                s = float(area) / float(area_unc)
                if s > best_sigma:
                    best_sigma = s
        ci = (ni.confidence.CI if ni.confidence is not None else 0.0)
        if best_sigma >= 10.0 and n_matched >= 2 and (ci or 0.0) >= 5.0:
            kept.append(ni)
            continue
        out_of_chain_suppressed.append(
            f"F-116: подавлен {nuc} (вне цепочки {chain_label}, "
            f"σ={best_sigma:.1f}, матчей={n_matched})"
        )
    return kept, out_of_chain_suppressed


# ──────────────────────────────────────────────────────────────────
# Result dataclasses
# ──────────────────────────────────────────────────────────────────

@dataclass
class StageResult:
    """Outcome of one identification stage."""
    stage: int                                       # 1, 2, or 3
    candidates_considered: List[str]
    detected: List[NuclideIdentification]            # after disambiguate
    rejected: List[NuclideIdentification]
    unmatched_peaks: List[FoundPeak]                 # after this stage
    notes: List[str] = field(default_factory=list)


@dataclass
class StagedAnalysisResult:
    """Full output of `analyze_lsrm_spe`."""
    spec: object                                     # Spectrum
    peaks: List[FoundPeak]
    detector_type: str
    fwhm_at_661: float                               # keV
    fwhm_model: Tuple[float, float, float]           # (a, b, c) in FWHM²=a+bE+cE²
    fwhm_model_source: str                           # "lsrm_peaks_table" / "default_NaI"
    stages: List[StageResult]
    final_detected: List[NuclideIdentification]      # union after all stages
    final_unmatched: List[FoundPeak]
    next_stage_recommended: Optional[int]            # None = stop, else 2 or 3
    next_stage_reason: str
    # Hints (raw matched text from filename / header)
    sample_type_hint: str
    geometry_hint: str
    detector_hint: str = ""
    # Canonical names (F-78) — preferred for downstream code, file lookup, docs
    sample_type_canonical: str = ""
    geometry_canonical: str = ""
    detector_canonical: str = ""
    is_background: bool = False                      # from filename hint
    # F-74 residual classification
    residual_classifications: List[ResidualClassification] = field(default_factory=list)
    # F-79/F-80 anchor-rank + express pattern results
    anchor_matches: List[AnchorMatch] = field(default_factory=list)
    pattern_confirmations: List[PatternConfirmation] = field(default_factory=list)
    # F-81 7-line ЕРН calibration verification (Lsrm methodology §9)
    seven_line_check: Optional[SevenLineCheck] = None
    # F-60 CI-gating (confirmed / tentative / noise tiers)
    ci_gating: Optional[CIGating] = None
    # F-61 Dose Contribution completeness metric
    completeness: Optional[CompletenessResult] = None
    # Analysis mode tag (per user methodology, mode is derived from
    # is_background): "background_7line" or "sample_anchor_rank"
    analysis_mode: str = ""
    # F-58 background subtraction (when a paired background path is given)
    background_subtraction: Optional[BackgroundSubtractionResult] = None
    # F-57 efficiency curve (auto-loaded by canonical geometry name)
    efficiency_curve: Optional[object] = None    # EfficiencyCurve | None
    efficiency_source: str = ""                  # path to .efr or ""
    # BUG-39 / BUG-40 (v1.22.0 Wave 6; F2-A renormalisation 2026-06-21) —
    # silent-fallback record for the detector profile. Populated
    # unconditionally; reason can be "profile_loaded_no_fallback"
    # (no warning) or "profile_not_on_disk" (canonical has no JSON
    # profile). The legacy "efficiency_tbd_using_fallback_profile" reason
    # was retired in F2-A together with the bogus Gamma-1S stub profile
    # that drove it. Surfaced in report.json `warnings` and in HTML/MD
    # reports when `reason != "profile_loaded_no_fallback"`.
    # Schema: {"requested": str, "actual": str, "reason": str, "human": str}.
    detector_fallback: Optional[dict] = None
    # F-84 / v1.13.0 — Round 5: activities, MDA, multiplet deconvolution
    # Activities per detected nuclide (Bq). Populated only when
    # `compute_activities=True` is requested AND an efficiency curve is
    # loaded. None when the call did not request activities.
    activities: Optional[List[ActivityResult]] = None
    # Specific activities (Bq/kg) derived by dividing each ActivityResult
    # by `sample_mass_kg`. Keyed by nuclide name. None when no mass given.
    specific_activities_Bq_per_kg: Optional[dict] = None
    sample_mass_kg: Optional[float] = None
    # Per-line ISO 11929 detection limits. Key: (nuclide_or_None, E_keV).
    # Populated for the standard MDA suite (Cs-137, Co-60, K-40, ²¹⁴Bi,
    # ²⁰⁸Tl, ²²⁸Ac) plus every detected nuclide's lines when the
    # efficiency curve is loaded. None when MDA was not requested.
    mda_per_line: Optional[dict] = None
    # Multiplet deconvolution outputs (one per cluster), when
    # `apply_deconvolution=True` is requested. Empty list when no
    # cluster was found; None when deconvolution was not requested.
    deconvolution_results: Optional[List[DeconvolutionResult]] = None
    # F-87 / v1.15.0 — Step 5β opt-in calibration refit diagnostic.
    # Populated only when `recalibrate_on_anchor_disagreement=True`.
    # Keys: attempted, applied, old_residual_max_keV, new_residual_max_keV,
    # old_energy_cal, new_energy_cal, n_anchors_used.
    recalibration_diag: dict = field(default_factory=dict)
    # F-88 / v1.15.1 — User-prioritized express anchor order +
    # chain-dominance hard-prior for Step 7 identification.
    # `priority_findings` is ordered per USER_PRIORITY_ORDER (1..6),
    # `chain_dominance` carries `th232`/`u238` flags and evidence.
    # `k40_ac228_overlap_warning` fires when Th-dominant AND K-40
    # priority signal matched — flags the Ac-228 1459.20 contamination
    # of the 1460.82 keV K-40 peak on NaI 63×63.
    priority_findings: List[PriorityFinding] = field(default_factory=list)
    chain_dominance: Optional[ChainDominance] = None
    k40_ac228_overlap_warning: bool = False
    # F-89 / v1.15.2 — filename binding hypothesis (SKILL.md §7A.1):
    # canonical isotope hints extracted from the filename drive the
    # Stage-1 candidate list and (via chain suppression) protect
    # against false-positive chain identification on single-isotope
    # sources. `chain_filtered_out` lists the nuclides dropped by the
    # F-89d suppression rule for transparency.
    filename_isotope_hints: List[str] = field(default_factory=list)
    filename_chains_claimed: List[str] = field(default_factory=list)
    chain_filtered_out: List[str] = field(default_factory=list)
    # F-89a / v1.15.2 — explicit background-subtraction status surfaced
    # so the report header can never silently omit it. One of:
    #   "subtracted_from_external_file"   — background_path was used
    #   "embedded_present_not_subtracted" — bg in same file but unused
    #   "absent_no_subtraction"           — neither; cps include nat. bg
    background_status: str = ""
    # F-129 / v1.17.7 — выбранный метод поиска пиков ("mariscotti" /
    # "convolution" / "compare") и результат сравнения двух методов
    # (только в режиме "compare"; иначе None).
    peak_search_method: str = "mariscotti"
    peak_search_method_comparison: Optional[dict] = None
    # F-131 / v1.17.7 — авто-поиск фонового спектра. ``auto_background_mode``
    # = "off" / "suggest" / "apply". ``auto_background_candidates`` — список
    # сериализованных кандидатов (см. BackgroundCandidate.to_dict). Может
    # быть пустым (кандидатов нет) либо None (auto-search не запускался).
    auto_background_mode: str = "off"
    auto_background_candidates: Optional[List[dict]] = None
    auto_background_applied_path: Optional[str] = None
    # F-332 / v1.18.18.5 — chart-toggle support. После background
    # subtraction `spec.counts` хранит NET; чтобы интерактивный отчёт
    # мог переключаться между {Образец / Фон / Наложение / Чистый}, мы
    # отдельно фиксируем gross sample counts + scaled bg counts на
    # сетке образца. None если фон не применялся.
    gross_counts: Optional[Any] = None              # ndarray sample's original counts
    bg_counts_on_sample_grid: Optional[Any] = None  # ndarray bg aligned + scaled to sample live-time
    bg_live_time: Optional[float] = None
    bg_scale_factor: Optional[float] = None
    # F-145 / v1.17.8 — Two-phase multiplet self-calibration diagnostic.
    # Сериализованная SelfCalibrationDiag из gamma.calibration.
    # multiplet_self_calibration. None — F-145 не запускалась (нет
    # мультиплетов, не NaI и т.п.).
    multiplet_self_calibration_diag: Optional[dict] = None
    # F-397 / v1.18.27 — детекция пиков в фоновом спектре (когда задействован
    # auto-bg или explicit background_path). Эти три поля содержат результат
    # независимого peak-search прогона по фону: используются HTML toggle
    # «Фон» для обновления peak block. None — фон не анализировался
    # (background_path не задан / complete_workflow=False / fallback fail).
    background_staged_result: Optional["StagedAnalysisResult"] = None
    # F-QC-01 / v1.19.1 — per-peak Poisson |z|-test results (BUG-35 / RAG-022).
    # Populated only when a background was subtracted AND peaks were found.
    # Structure: {"peak_z_roi": [{"e_lo_keV": float, "e_hi_keV": float,
    #   "peak_energy_keV": float, "z": float, "abs_z": float, "tier": str,
    #   "passed": bool, "B1": int, "B2": int, "note": str}, ...]}
    # None when background not subtracted or no peaks found.
    bg_quality_check: Optional[dict] = None
    notes: List[str] = field(default_factory=list)


# ──────────────────────────────────────────────────────────────────
# F-378 / v1.18.25 — sample mass mismatch detector
# ──────────────────────────────────────────────────────────────────

# Tolerance (relative) below which CLI vs .spe mass values are treated
# as the same. 1 % covers float rounding in the SAMPLEMASS LSRM field
# (stored in grams, converted to kg) without catching real conflicts.
_MASS_MISMATCH_REL_TOL = 0.01


def check_sample_mass_mismatch(
    cli_mass_kg: Optional[float],
    spec_extras: Optional[dict],
    *,
    rel_tol: float = _MASS_MISMATCH_REL_TOL,
) -> Optional[str]:
    """F-378: detect a conflict between user-supplied --sample-mass-kg
    and the SAMPLEMASS field stored in the .spe metadata.

    The .spe LSRM format carries the true sample mass (used by the
    instrument operator at the measurement step). When the CLI flag
    is also passed and the two disagree by more than ``rel_tol``
    (default 1 %), the specific activity computed downstream will be
    silently wrong by ``cli / spec``× — and the passport comparison
    will diverge by the same factor. This helper returns a warning
    string in such a case; the caller appends it to report notes and
    mirrors it to stderr.

    Args:
        cli_mass_kg: value of the CLI/kwarg --sample-mass-kg
            (None if the user did not pass it — no conflict possible).
        spec_extras: ``spec.extras`` from the LSRM reader; we read
            ``lsrm_sample_mass_kg`` (already converted from grams).
        rel_tol: relative tolerance below which values are treated as
            equal. Default 1 % (covers float rounding).

    Returns:
        Warning string (Russian, prefixed with ⚠ F-378) when both
        values are present and disagree beyond ``rel_tol``; None
        otherwise.

    Side effects: none — pure function for testability.
    """
    if cli_mass_kg is None:
        return None
    if not spec_extras:
        return None
    spec_mass_raw = spec_extras.get("lsrm_sample_mass_kg")
    if spec_mass_raw is None:
        return None
    try:
        spec_mass_kg = float(spec_mass_raw)
    except (TypeError, ValueError):
        return None
    if spec_mass_kg <= 0:
        return None
    rel_diff = abs(float(cli_mass_kg) - spec_mass_kg) / spec_mass_kg
    if rel_diff <= rel_tol:
        return None
    factor = float(cli_mass_kg) / spec_mass_kg
    return (
        f"⚠ F-378: расхождение массы образца — передано "
        f"--sample-mass-kg = {float(cli_mass_kg):.3f} кг, в .spe "
        f"SAMPLEMASS = {spec_mass_kg:.3f} кг (разница "
        f"{rel_diff * 100:.1f} %). Используется значение из CLI; "
        f"удельная активность может быть искажена в {factor:.2f}× "
        f"раз. Перепроверьте --sample-mass-kg или уберите флаг, "
        f"чтобы скил взял массу из .spe."
    )


# ──────────────────────────────────────────────────────────────────
# Empirical FWHM model from stored LSRM peaks table
# ──────────────────────────────────────────────────────────────────

# Default NaI 63×63 (Gamma-1S USB SN-01) fallback.
#
# F-125 / v1.17.6 — REFIT по 26 анкер-точкам из калибровочных
# спектров Gamma-1S (Cs-137 / Co-60 / Eu-152 / Th-228 / Ra-226 /
# K-40 в Marinelli и точечной геометрии). Использован constrained
# fit FWHM²(E) = a + b·E + c·E² с границами a ≥ 0, b ≥ 0 (физически
# обоснованная неотрицательность).
#
# Метрики аппроксимации:
#   • rms ≈ 2.24 кэВ на 26 точках (84–2615 кэВ)
#   • E=661.66 → FWHM=46.95 кэВ (R=7.1 %; измерено 45.74)
#   • E=2614.5 → FWHM=107.95 кэВ (измерено 106.99)
#
# Предыдущая модель v1.17.5 = (-980.7, 6.846, -0.000797) была
# откалибрована только на 3 точках (84, 238, 1460, 2614 кэВ из
# `Фон Вода 19-11-2025`) и завышала FWHM в области 600-1500 кэВ
# на 8-15 % (M1 χ²/ν → 31.72 вместо контракта 17.02 v1.17.2).
_DEFAULT_NAI_FWHM_MODEL = (0.0, 2.950048, 0.000576400)


# F-452 (2026-06-21) — единый объект FWHM-модели для build_fwhm_model.
#
# До F-452 build_fwhm_model возвращал плоский 3-tuple (a,b,c) формы
# FWHM²(E) = a + b·E + c·E² (квадратичная-в-E). LSRM-полином 4-й степени
# в sqrt(E) — реальная физическая модель FWHM(E) для калибровочных
# спектров Гамма-1С (LSRM "Алгоритмические основы" §8.3, RAG-043) — при
# попытке fit-ить его как 3-параметрическую квадратичную (F-160
# ground-truth ветка, NNLS) теряет ~5-7 кэВ точности на anchor-точках
# (max|ΔFWHM|, документировано в F-160 ALERT-сообщениях).
#
# FwhmModel позволяет вернуть LSRM-полином ЧЕСТНО как poly-4 sqrt(E),
# сохранив обратную совместимость для квадратичной ветки (bootstrap,
# lsrm_peaks_table, default_NaI).
@dataclass(frozen=True)
class FwhmModel:
    """FWHM(E) модель из build_fwhm_model.

    Две формы:
      - ``quad_fwhm2_in_E``: FWHM(E) = √(max(a + b·E + c·E², 0.01)),
        ``coefficients = (a, b, c)``. Используется bootstrap (F-449),
        lsrm_peaks_table (3+ rows quadratic / 1-2 rows α·√E),
        default_NaI_63x63 fallback.
      - ``lsrm_poly_sqrt_E``: FWHM(E) = Σ_{k=0..n} c_k · √E^k,
        ``coefficients = (c_0, c_1, ..., c_n)``, типично n=4 (LSRM
        канонический полином). Используется F-160 ground-truth ветка
        (lstsq poly-4 на anchor-точках из
        ``references/lsrm_ground_truth/<base>/fwhm_calibration_lsrm.json``).

    Класс CALLABLE: ``model(E_keV)`` возвращает FWHM в кэВ. Это
    позволяет всем потребителям (``_make_fwhm_at_channel``,
    ``fwhm_provider_keV=lambda E: fwhm_keV_at_energy(...)``)
    работать model-agnostic.

    Числовой floor (max(val, 0.01)) сохранён в quad-форме как
    last-resort safety net (BUG-43 / RAG-043), не имеет физического
    смысла. В LSRM-форме floor = 0.1 кэВ (полином by-construction
    положителен на калибровочном диапазоне).

    F-452-FU (2026-06-22) — линейная экстраполяция выше E_max_anchor:
    poly-4 в √E неограниченно экстраполируется за anchor-диапазон
    (LSRM anchor-точки типично до ~1.6 MeV; Tl-208 главная линия 2614 keV
    лежит в чистой extrapolation). Без clamp коэффициент c_4 (часто
    отрицательный) приводит к runaway-FWHM, что разрушает peak-fit на
    high-E пиках (Tl-208 = 0 в Th-232 fixture — root cause 3 failing tests).
    Когда заданы все три ``linear_extrap_*``-поля, для E > linear_extrap_start_E
    возвращается ``intercept + slope·(E - start_E)`` (slope >= 0,
    physical sanity: FWHM монотонно растёт с E).
    """

    kind: str  # "quad_fwhm2_in_E" | "lsrm_poly_sqrt_E"
    coefficients: Tuple[float, ...]
    # F-452-FU (2026-06-22) — clamp выше E_max_anchor (линейная extrap).
    # Все три поля задаются вместе или все остаются None.
    linear_extrap_start_E: Optional[float] = None
    linear_extrap_intercept: Optional[float] = None
    linear_extrap_slope: Optional[float] = None

    def __call__(self, E_keV: float) -> float:
        E = max(float(E_keV), 5.0)
        if (
            self.linear_extrap_start_E is not None
            and self.linear_extrap_intercept is not None
            and self.linear_extrap_slope is not None
            and E > float(self.linear_extrap_start_E)
        ):
            val = (
                float(self.linear_extrap_intercept)
                + float(self.linear_extrap_slope)
                * (E - float(self.linear_extrap_start_E))
            )
            return max(val, 0.1)
        if self.kind == "quad_fwhm2_in_E":
            a, b, c = self.coefficients
            val = a + b * E + c * E * E
            return math.sqrt(max(val, 0.01))
        if self.kind == "lsrm_poly_sqrt_E":
            z = math.sqrt(E)
            out = 0.0
            for k, c in enumerate(self.coefficients):
                out += float(c) * (z ** k)
            return max(out, 0.1)
        raise ValueError(
            f"Unknown FwhmModel kind: {self.kind!r} "
            f"(expected 'quad_fwhm2_in_E' or 'lsrm_poly_sqrt_E')"
        )


def _wrap_quad(coefs: Tuple[float, float, float]) -> "FwhmModel":
    """Helper: упаковать legacy (a,b,c) в FwhmModel(kind='quad_fwhm2_in_E')."""
    return FwhmModel(kind="quad_fwhm2_in_E", coefficients=tuple(float(c) for c in coefs))


_DEFAULT_NAI_FWHM_MODEL_OBJ = _wrap_quad(_DEFAULT_NAI_FWHM_MODEL)


def _fit_alpha_sqrt_E_model(
    Es: list, Fs: list
) -> Tuple[float, float, float]:
    """
    1-parameter model FWHM(E) = α·√E for scintillators with poor counting
    at high E. Returned in the same (a, b, c) form for FWHM²(E):
        FWHM²(E) = α²·E  →  (a, b, c) = (0, α², 0)
    """
    arr_E = np.array(Es); arr_F = np.array(Fs)
    # Weighted by 1/F (typical Gilmore convention) — give equal weight here
    alpha2 = np.mean((arr_F ** 2) / arr_E)
    return (0.0, float(alpha2), 0.0)


# BUG-41 / Wave 7 / 2026-06-05 — pathology-detection thresholds for the
# `lsrm_peaks_table_quadratic` fit. The LSRM PEAKS= table is sometimes too
# sparse or too clustered to constrain `FWHM^2(E) = a + b*E + c*E^2` over the
# full spectrum range. When the resulting discriminant goes negative at low
# E (BUG-37/41 root cause for the AmTiCsEu Marinelli fixture: val = -174.5
# + 2.67*E + 0.0004*E^2 -> val(59.5) = -15.4), the 0.01 floor in
# `fwhm_keV_at_energy` clamps FWHM to 0.1 keV — far below real NaI 63x63
# FWHM at 60 keV (~12 keV). Match window 0.3 keV vs required 36 keV ->
# Am-241 59.54 keV never matched.
#
# Sanity-check the fitted quadratic at the **identification-critical**
# low/mid energies. Test energies are picked to cover Am-241 59.54 keV
# (the BUG-41 motivating case), the Ba-133 81 keV / 122 keV cluster,
# Cs-137 661 keV, K-40 1460 keV, Co-60 1332 keV / 1173 keV. Below 60 keV
# the quadratic-in-E approximation of a true sqrt(E)-shaped FWHM curve
# routinely produces small (or even slightly negative) values without
# being identification-blocking — so we do NOT test E<60 keV, to avoid
# false-positive rejections of physically reasonable fits. (Verified by
# sampling the operator-certified LSRM polynomial
# ``(-6.2914, 1.7443, 0.004948)`` and refitting as FWHM^2 quadratic:
# val(20 keV) = -59, val(40 keV) = -9 — both **below threshold**, yet
# the underlying polynomial is the canonical Gamma-1S calibration.)
#
# Threshold 1.0 keV^2 -> FWHM=1.0 keV. Real NaI 63x63 FWHM at the test
# energies is far above this floor (>=12 keV at 60 keV per default-NaI
# model), so any fit returning val < 1 keV^2 at any of these E is
# unambiguously broken (typical: negative-discriminant from a sparse or
# clustered ``lsrm_peaks_table`` -> matches window collapses to <1 keV).
_FWHM_PATHOLOGY_TEST_ENERGIES_keV = (60.0, 100.0, 200.0, 500.0, 1000.0)
_FWHM_PATHOLOGY_VAL_THRESHOLD_keV2 = 1.0


def _eval_fwhm2_quadratic(model: Tuple[float, float, float], E: float) -> float:
    """Evaluate FWHM^2(E) = a + b*E + c*E^2 for a quadratic-in-E model."""
    a, b, c = model
    return a + b * E + c * E * E


def _eval_lsrm_sqrt_E_polynomial(coefs: tuple, E: float) -> float:
    """Evaluate the LSRM stored FWHM polynomial.

    Per LSRM "Algorithmic Foundations" section 8.3 (RAG-043,
    ``lsrm_act_2014``) and verified in
    ``scripts/gamma/io/lsrm_spe.py:44-68``, the stored polynomial is
        FWHM_keV(E) = sum_k c_k * z^k,   z = sqrt(E_keV)
    NOT a polynomial in E directly despite the legacy label
    ``lsrm_fwhm_polynomial_in_E``. Used as the preferred BUG-41
    fallback when the empirical `lsrm_peaks_table` quadratic fit is
    pathological at low E.
    """
    if not coefs or E <= 0:
        return float("nan")
    z = math.sqrt(float(E))
    out = 0.0
    for k, c in enumerate(coefs):
        out += float(c) * (z ** k)
    return out


def _convert_lsrm_sqrt_E_to_fwhm2_quadratic(
    coefs: tuple,
    sample_energies_keV: Sequence[float] = (
        50.0, 80.0, 122.0, 200.0, 500.0, 1000.0, 1500.0, 2000.0,
    ),
) -> Optional[Tuple[float, float, float]]:
    """Convert the stored LSRM sqrt(E) polynomial into the ``(a, b, c)``
    form of ``FWHM^2(E) = a + b*E + c*E^2`` used by the rest of the
    pipeline.

    Samples the LSRM polynomial at a spread of energies covering the
    physical valid range (default: 50 .. 2000 keV) and least-squares
    fits ``FWHM^2(E_i) ~ a + b*E_i + c*E_i^2``. Sampling starts at 50 keV
    because the LSRM polynomial extrapolates to unphysical values below
    its calibration range (typically ~84 keV — the lowest PEAKS= row in
    the Gamma-1S corpus).

    Returns None if the LSRM polynomial itself is broken (NaN or all
    non-positive over the sample grid). Reference for the polynomial
    form: ``scripts/gamma/io/lsrm_spe.py:44-68`` (BUG-22 sqrt(E)
    argument convention).
    """
    Es: list = []
    F2s: list = []
    for E in sample_energies_keV:
        fwhm = _eval_lsrm_sqrt_E_polynomial(coefs, E)
        if not math.isfinite(fwhm) or fwhm <= 0:
            continue
        Es.append(float(E))
        F2s.append(float(fwhm * fwhm))
    if len(Es) < 3:
        return None
    Es_a = np.asarray(Es, dtype=np.float64)
    F2_a = np.asarray(F2s, dtype=np.float64)
    A = np.vstack([np.ones_like(Es_a), Es_a, Es_a ** 2]).T
    sol, *_ = np.linalg.lstsq(A, F2_a, rcond=None)
    return (float(sol[0]), float(sol[1]), float(sol[2]))


def _model_is_pathological(
    model: Tuple[float, float, float],
    test_energies_keV: Sequence[float] = _FWHM_PATHOLOGY_TEST_ENERGIES_keV,
    val_threshold_keV2: float = _FWHM_PATHOLOGY_VAL_THRESHOLD_keV2,
) -> bool:
    """True iff ``FWHM^2(E)`` drops below ``val_threshold_keV2`` for any
    ``E`` in ``test_energies_keV``. A pathological quadratic produces a
    spuriously narrow FWHM at low E that downstream identification
    windows cannot recover from (BUG-41).
    """
    for E in test_energies_keV:
        if _eval_fwhm2_quadratic(model, E) < val_threshold_keV2:
            return True
    return False


def _record_fwhm_fallback_warning(spec, message: str) -> None:
    """Push a fallback warning into ``spec.extras["fwhm_model_warnings"]``
    so downstream reporting can surface it on the
    ``report.json:warnings`` channel. The list is created lazily;
    preserves any pre-existing entries from other calls in the same
    pipeline run.
    """
    if spec.extras is None:
        spec.extras = {}
    bucket = spec.extras.setdefault("fwhm_model_warnings", [])
    bucket.append(str(message))


def _resolve_pathology_fallback(
    spec,
    *,
    original_source: str,
    original_model: Tuple[float, float, float],
) -> Tuple[Tuple[float, float, float], str]:
    """Pick the BUG-41 fallback model when the empirical fit went
    pathological.

    Strategy:
      1. **Prefer LSRM stored sqrt(E) polynomial** (operator-certified
         for the calibrated detector). Convert it to ``FWHM^2(E)`` form
         by sampling at energies >= 50 keV (above the LSRM lowest stored
         anchor in our fixtures).
      2. **Cross-check** the converted LSRM model against
         ``_DEFAULT_NAI_FWHM_MODEL`` at a low-E checkpoint (E=60 keV).
         If the LSRM-converted FWHM at 60 keV is at least 70 % of the
         default-NaI FWHM at 60 keV (~12 keV), keep LSRM stored;
         otherwise the polynomial is extrapolating below its calibration
         range and we use ``_DEFAULT_NAI_FWHM_MODEL`` (fit on 26 anchors
         across 84-2614 keV).
      3. If the LSRM polynomial is itself unusable (NaN / all
         non-positive over the sample grid, or no
         ``StoredFwhmCalibration`` present), fall back to
         ``_DEFAULT_NAI_FWHM_MODEL`` directly.

    The brief preference order (revalidation outbox section 5 BUG-41)
    is LSRM-stored > default; we honour it except where step 2's
    consistency check fires.

    Always records a ``fwhm_model_warnings`` entry on ``spec.extras``
    so the report layer can surface it on ``report.json:warnings``.

    Source labels use ASCII suffixes (``_rejected_to_<fallback>``) to
    avoid console-encoding issues on cp1251 hosts.
    """
    sf = getattr(spec, "stored_fwhm_calibration", None)
    if (
        sf is not None
        and sf.coefficients
        and sf.model == "lsrm_fwhm_polynomial_in_E"
    ):
        converted = _convert_lsrm_sqrt_E_to_fwhm2_quadratic(sf.coefficients)
        if converted is not None:
            E_check = 60.0
            lsrm_check = math.sqrt(
                max(_eval_fwhm2_quadratic(converted, E_check), 0.01)
            )
            default_check = math.sqrt(
                max(
                    _eval_fwhm2_quadratic(_DEFAULT_NAI_FWHM_MODEL, E_check),
                    0.01,
                )
            )
            converted_ok = not _model_is_pathological(
                converted,
                test_energies_keV=(60.0, 100.0, 200.0, 500.0),
            )
            if converted_ok and lsrm_check >= 0.70 * default_check:
                _record_fwhm_fallback_warning(
                    spec,
                    f"BUG-41: FWHM model '{original_source}' fit "
                    f"a={original_model[0]:.3g} b={original_model[1]:.3g} "
                    f"c={original_model[2]:.3g} is pathological at low E "
                    f"(FWHM2<{_FWHM_PATHOLOGY_VAL_THRESHOLD_keV2} keV2 at "
                    f"one of {tuple(_FWHM_PATHOLOGY_TEST_ENERGIES_keV)} keV); "
                    "falling back to LSRM stored sqrt(E) polynomial "
                    f"coefs={tuple(round(c, 6) for c in sf.coefficients)} "
                    f"-> (a,b,c)=({converted[0]:.4g}, {converted[1]:.4g}, "
                    f"{converted[2]:.4g}). FWHM(60 keV) via converted="
                    f"{lsrm_check:.2f} keV (vs default NaI 63x63 "
                    f"= {default_check:.2f} keV).",
                )
                return converted, (
                    f"{original_source}_rejected_to_lsrm_stored_sqrtE"
                )
            _record_fwhm_fallback_warning(
                spec,
                f"BUG-41: FWHM model '{original_source}' pathological; "
                f"LSRM stored sqrt(E) polynomial extrapolation at E=60 keV "
                f"gives FWHM={lsrm_check:.2f} keV (< 70% of default NaI "
                f"63x63 {default_check:.2f} keV) -> extrapolation regime; "
                "using _DEFAULT_NAI_FWHM_MODEL (fit on 26 anchors across "
                "84-2614 keV).",
            )
            return _DEFAULT_NAI_FWHM_MODEL, (
                f"{original_source}_rejected_to_default_NaI_63x63"
            )
        _record_fwhm_fallback_warning(
            spec,
            f"BUG-41: FWHM model '{original_source}' pathological AND "
            f"LSRM stored sqrt(E) polynomial conversion failed "
            f"(coefficients={tuple(sf.coefficients)}); "
            "falling back to _DEFAULT_NAI_FWHM_MODEL.",
        )
        return _DEFAULT_NAI_FWHM_MODEL, (
            f"{original_source}_rejected_to_default_NaI_63x63"
        )
    _record_fwhm_fallback_warning(
        spec,
        f"BUG-41: FWHM model '{original_source}' pathological and no "
        "StoredFwhmCalibration available; falling back to "
        "_DEFAULT_NAI_FWHM_MODEL.",
    )
    return _DEFAULT_NAI_FWHM_MODEL, (
        f"{original_source}_rejected_to_default_NaI_63x63"
    )


# F-160 (2026-06-20) — auto-loader LSRM ground-truth FWHM-anchors.
# Оператор: «всегда нужно делать калибровку по fwhm полагаться на
# расчетную кривую нельзя». Эталонные спектры поверочного набора
# Гамма-1С имеют зафиксированную LSRM-кривую FWHM(E) в
# references/lsrm_ground_truth/<base>/fwhm_calibration_lsrm.json
# (14 anchor-точек, polynomial deg=4, χ²=2.0332 для Th232 Marinelli).
# Маппинг .spe-basename → папка определяется в _index.json (явный,
# без транслитерации; «Маринелли» в .spe vs «Marinelli» в references/).
_GROUND_TRUTH_INDEX_PATH = Path("references/lsrm_ground_truth/_index.json")


def _load_ground_truth_fwhm_anchors(spec):
    """Если для spec.source_path есть mapping в _index.json → загрузить
    anchor-точки (E_keV, fwhm_keV_measured) из соответствующего
    fwhm_calibration_lsrm.json. Иначе None.

    Возвращает: (list[float] Es, list[float] Fs, str gt_folder_name) или None.
    """
    sp = getattr(spec, "source_path", None)
    if not sp:
        return None
    base = Path(str(sp)).name
    try:
        idx_data = json.loads(_GROUND_TRUTH_INDEX_PATH.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError, OSError) as _gt_idx_exc:
        logger.debug(
            "F-160: _index.json не прочитан (%s); auto-loader пропущен.",
            _gt_idx_exc,
        )
        return None
    mapping = idx_data.get("mapping", {}) if isinstance(idx_data, dict) else {}
    gt_folder = mapping.get(base)
    if not gt_folder:
        return None
    gt_dir = _GROUND_TRUTH_INDEX_PATH.parent / gt_folder
    fwhm_json = gt_dir / "fwhm_calibration_lsrm.json"
    if not fwhm_json.exists():
        logger.warning(
            "F-160: mapping для %s указывает на %s, но %s не существует.",
            base, gt_folder, fwhm_json,
        )
        return None
    try:
        data = json.loads(fwhm_json.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as _gt_json_exc:
        logger.warning(
            "F-160: %s не читается (%s); auto-loader пропущен.",
            fwhm_json, _gt_json_exc,
        )
        return None
    anchors = data.get("anchors", []) if isinstance(data, dict) else []
    Es, Fs = [], []
    for a in anchors:
        e = a.get("energy_keV")
        f = a.get("fwhm_keV_measured")
        if e and f and e > 0 and f > 0:
            Es.append(float(e))
            Fs.append(float(f))
    if len(Es) < 3:
        logger.warning(
            "F-160: %s содержит только %d валидных anchor-точек (нужно ≥3); "
            "auto-loader пропущен.", fwhm_json, len(Es),
        )
        return None
    return Es, Fs, gt_folder


# F-449 / agent-a-math-3 (2026-06-20) — hardening Правки B.
#
# Bootstrap guard constants — see _bootstrap_fwhm_from_significant_peaks.
# - Sqrt-form FWHM(E) = sqrt(a + b·E)  (NaI scintillator standard, 2 params).
# - Bounds on (a, b):
#     a ∈ [0, 100]   (electronic noise, typ. 0..20 keV²)
#     b ∈ [0.3, 15]  (physical NaI dispersion, typ. 1..6)
# - Monotonicity grid (FWHM must grow monotonically on this energy grid).
# - Distribution-guard: ≥6 anchors total, ≥1 anchor in each of 3 sub-bands.
# - Sanity-vs-default: bootstrap must not diverge from _DEFAULT_NAI_FWHM_MODEL
#   by more than 20 keV (absolute) or 30 % (relative) on the comparison grid.
_BOOTSTRAP_A_MIN_keV2 = 0.0
_BOOTSTRAP_A_MAX_keV2 = 100.0
_BOOTSTRAP_B_MIN_keV = 0.3
_BOOTSTRAP_B_MAX_keV = 15.0
_BOOTSTRAP_MONOTONICITY_GRID_keV = (60.0, 100.0, 200.0, 500.0, 1000.0,
                                    1500.0, 2000.0, 2614.0, 3000.0)
_BOOTSTRAP_MIN_ANCHORS = 6
_BOOTSTRAP_DIST_BANDS_keV = ((200.0, 800.0), (800.0, 1500.0),
                             (1500.0, 2700.0))
_BOOTSTRAP_SANITY_GRID_keV = (208.0, 238.0, 338.0, 460.0, 580.0,
                              900.0, 953.0, 1581.0, 2613.0)
_BOOTSTRAP_SANITY_ABS_MAX_keV = 20.0
_BOOTSTRAP_SANITY_REL_MAX = 0.30


def _bootstrap_fwhm_from_significant_peaks(
    spec,
    *,
    min_anchors: int = _BOOTSTRAP_MIN_ANCHORS,
    energy_min_keV: float = 200.0,
    energy_max_keV: float = 2700.0,
    sn_min: float = 5.0,
    height_min: float = 10.0,
    isolation_factor: float = 1.5,
) -> Optional[Tuple[Tuple[float, float, float], str, int]]:
    """F-449 / agent-a-math-3 (2026-06-20) — bootstrap FWHM(E) from
    significant, isolated peaks of the spectrum itself when no
    trustworthy stored curve is available. **Hardened from
    agent-a-math-2 prototype** (Правка B v2).

    Methodology (CLAUDE.md «FWHM любого пика», method 2 of the F-rule):
      1. Run a seeded Mariscotti search with σ_th=sn_min.
      2. Filter for high-S/N, isolated, in-band anchors:
           - significance ≥ sn_min
           - net height ≥ height_min  (10 counts; Mariscotti S/N gate
             already filters faint peaks more aggressively than the
             agent-a-math-2 prototype's 30, which dropped real Th-232
             ³⁴⁴ daughter lines on this corpus).
           - E ∈ [energy_min_keV, energy_max_keV]
           - isolation: nearest other Mariscotti candidate ≥
             ``isolation_factor·seed_fwhm_ch`` (default 1.5 instead of
             the prototype's 2.0 — at 1.5·FWHM two equal Gaussians dip
             ~50 % between centres, which is the textbook FWHM-criterion
             for "resolved peaks"; 2.0 was over-strict and dropped the
             Th-232 583/635 + 911/969 pairs that ARE physically usable).
      3. Measure FWHM in channels by half-height crossing on raw counts
         (fwhm_provider._measure_fwhm_channels).
      4. Convert to keV via local gain dE/dch from spec.energy_cal.
      5. Distribution-guard (a-priori, BEFORE the fit):
           - total ≥ ``_BOOTSTRAP_MIN_ANCHORS`` (6)
           - at least 1 anchor in each band of
             ``_BOOTSTRAP_DIST_BANDS_keV`` (200..800, 800..1500,
             1500..2700) — sqrt-fit on a single band is unconstrained
             and gives unphysical extrapolation (Ra-226 prototype: 4
             anchors all 600..1800 keV → b<0 → FWHM(60)=66 keV).
      6. Fit the **sqrt-form** ``FWHM²(E) = a + b·E`` by unconstrained
         lstsq on Y=FWHM² vs X=[1, E]. 2-parameter form (NaI standard
         √E for scintillator dispersion) is **statistically more
         stable** on small samples than the 3-parameter polynomial that
         the prototype used (which over-fitted to 4 high-E points and
         flipped b negative).
      7. Apply guards in order:
           (a) bounds: a ∈ [0, 100], b ∈ [0.3, 15] — outside ⇒ reject
           (b) monotonicity on ``_BOOTSTRAP_MONOTONICITY_GRID_keV``:
               FWHM(E_{i+1}) > FWHM(E_i) for every neighbouring pair
           (c) sanity vs ``_DEFAULT_NAI_FWHM_MODEL`` on
               ``_BOOTSTRAP_SANITY_GRID_keV``: ``max |ΔFWHM|`` must be
               ≤ ``_BOOTSTRAP_SANITY_ABS_MAX_keV`` (20 keV) AND
               ≤ ``_BOOTSTRAP_SANITY_REL_MAX`` (30 %) of the default
               curve at the worst point. Default-NaI is generic but
               physically reasonable; bootstrap must lie in its
               neighbourhood, not diverge in 2× as the prototype did.
         Any failure ⇒ reject (caller falls back to
         ``_DEFAULT_NAI_FWHM_MODEL``).

    Returns ``(model, source_label, n_anchors)`` on success, else None.

    The anchor energy window [energy_min_keV, energy_max_keV] is
    deliberately conservative: below 200 keV NaI FWHM curve is steep
    and a single mis-isolated XRF line can pull the fit; above 2700 keV
    only Tl-208 2614 anchors usefully exist on Th-232 spectra.
    """
    counts = getattr(spec, "counts", None)
    if counts is None or len(counts) < 100:
        return None
    if not spec.energy_cal or len(spec.energy_cal) < 2:
        return None
    try:
        from gamma.peaks.search import mariscotti_search
        from gamma.calibration.fwhm_provider import _measure_fwhm_channels
    except Exception as _bs_imp_exc:
        logger.warning(
            "F-449 bootstrap FWHM: helper import failed (%s); "
            "falling back to default NaI curve.", _bs_imp_exc,
        )
        return None
    counts_arr = np.asarray(counts, dtype=np.float64)
    # Seed FWHM in channels for the search. Use the default NaI 63x63
    # curve evaluated at 661 keV as a reasonable scintillator seed:
    # ~47 keV / gain. This is only a SEED for kernel sizing; the actual
    # FWHM is measured from counts (half-height crossing).
    gain = abs(float(spec.energy_cal[1]))
    if gain <= 0:
        return None
    seed_fwhm_keV = math.sqrt(_eval_fwhm2_quadratic(
        _DEFAULT_NAI_FWHM_MODEL, 661.0
    ))
    seed_fwhm_ch = max(3.0, seed_fwhm_keV / gain)
    try:
        found = mariscotti_search(
            counts=counts_arr,
            fwhm_channels=float(seed_fwhm_ch),
            sigma_threshold=float(sn_min),
            min_separation_factor=1.0,
            edge_margin=10,
        )
    except Exception as _bs_search_exc:
        logger.warning(
            "F-449 bootstrap FWHM: Mariscotti seed search failed (%s); "
            "no anchors collected.", _bs_search_exc,
        )
        return None
    if not found:
        return None

    iso_dist_ch = isolation_factor * seed_fwhm_ch
    _DOMINANT_RATIO = 2.0  # F-rule «всегда калибровать»: dominant-of-doublet rescue
    anchors: List[Tuple[float, float]] = []  # (E_keV, FWHM_keV)
    for p in found:
        if getattr(p, "significance", 0.0) < sn_min:
            continue
        if getattr(p, "height", 0.0) < height_min:
            continue
        E_keV = spec.channel_to_energy(int(p.channel))
        if E_keV is None:
            continue
        E_keV = float(E_keV)
        if not (energy_min_keV <= E_keV <= energy_max_keV):
            continue
        # Isolation: nearest other Mariscotti candidate ≥ isolation_factor *
        # seed_fwhm_ch away. At 1.5·FWHM equal Gaussians dip to ~50 %
        # between centres (FWHM-criterion). Dominant-of-doublet rescue
        # (F-160 «всегда калибровать»): если данный кандидат имеет
        # значимость ≥ _DOMINANT_RATIO × значимость каждого «слишком
        # близкого» соседа, он эффективно действует как singlet
        # (Mariscotti FWHM на доминирующем пике слабо искажается слабым
        # соседом). Без этого Ac-228 911 (σ=71) бракуется парой к
        # Ac-228 968 (σ=33) → пустая bootstrap-полоса [800,1500) на
        # Th-232 → весь bootstrap rejected → fallback на default NaI.
        too_close = False
        p_sig = float(getattr(p, "significance", 0.0))
        for q in found:
            if q.channel == p.channel:
                continue
            if abs(q.channel - p.channel) < iso_dist_ch:
                q_sig = float(getattr(q, "significance", 0.0))
                if p_sig < _DOMINANT_RATIO * q_sig:
                    too_close = True
                    break
        if too_close:
            continue
        # Measure FWHM in channels from raw counts.
        fwhm_ch = _measure_fwhm_channels(counts_arr, int(p.channel), seed_fwhm_ch)
        if fwhm_ch is None or fwhm_ch <= 1.0:
            continue
        # Local gain dE/dch from polynomial e-cal.
        ch_f = float(p.channel)
        dE_dN = sum(
            i * float(a) * (ch_f ** (i - 1))
            for i, a in enumerate(spec.energy_cal) if i > 0
        )
        if dE_dN <= 0:
            continue
        fwhm_keV = fwhm_ch * dE_dN
        if fwhm_keV <= 0:
            continue
        anchors.append((E_keV, fwhm_keV))

    # Guard (d) — distribution: total + per-band coverage.
    if len(anchors) < min_anchors:
        logger.info(
            "F-449 bootstrap FWHM: only %d anchors (need ≥%d); "
            "falling back to default NaI curve.",
            len(anchors), min_anchors,
        )
        return None
    Es = np.array([a[0] for a in anchors], dtype=np.float64)
    Fs = np.array([a[1] for a in anchors], dtype=np.float64)
    band_counts = []
    for lo, hi in _BOOTSTRAP_DIST_BANDS_keV:
        n_in_band = int(((Es >= lo) & (Es < hi)).sum())
        band_counts.append(n_in_band)
        if n_in_band < 1:
            logger.info(
                "F-449 bootstrap FWHM: no anchors in band [%g,%g) "
                "(band coverage %s); falling back to default NaI curve.",
                lo, hi, band_counts,
            )
            return None

    # Fit sqrt-form FWHM² = a + b·E (2 params, lstsq). NaI standard.
    A = np.vstack([np.ones_like(Es), Es]).T
    try:
        sol, *_ = np.linalg.lstsq(A, Fs * Fs, rcond=None)
    except Exception as _bs_lstsq_exc:
        logger.warning(
            "F-449 bootstrap FWHM: sqrt-form lstsq failed (%s); "
            "rejecting bootstrap.", _bs_lstsq_exc,
        )
        return None
    a_fit, b_fit = float(sol[0]), float(sol[1])
    # Express in the (a, b, c) tuple so the rest of the pipeline keeps
    # using a single quadratic-in-E representation. c = 0 for sqrt-form.
    model: Tuple[float, float, float] = (a_fit, b_fit, 0.0)

    # Guard (b1) — coefficient bounds.
    if not (_BOOTSTRAP_A_MIN_keV2 <= a_fit <= _BOOTSTRAP_A_MAX_keV2):
        logger.info(
            "F-449 bootstrap FWHM: a=%.4g outside bounds [%.1f, %.1f]; "
            "falling back to default NaI curve.",
            a_fit, _BOOTSTRAP_A_MIN_keV2, _BOOTSTRAP_A_MAX_keV2,
        )
        return None
    if not (_BOOTSTRAP_B_MIN_keV <= b_fit <= _BOOTSTRAP_B_MAX_keV):
        logger.info(
            "F-449 bootstrap FWHM: b=%.4g outside bounds [%.2f, %.1f]; "
            "falling back to default NaI curve.",
            b_fit, _BOOTSTRAP_B_MIN_keV, _BOOTSTRAP_B_MAX_keV,
        )
        return None

    # Guard (c) — monotonicity on grid (defensive; b≥0.3 already
    # implies monotone √(a+bE), but check explicitly to catch any
    # future model-form change).
    grid = _BOOTSTRAP_MONOTONICITY_GRID_keV
    prev = math.sqrt(max(_eval_fwhm2_quadratic(model, grid[0]), 0.0))
    for E_next in grid[1:]:
        nxt = math.sqrt(max(_eval_fwhm2_quadratic(model, E_next), 0.0))
        if not (nxt > prev):
            logger.info(
                "F-449 bootstrap FWHM: monotonicity violated at %g→%g keV "
                "(%.2f → %.2f); falling back to default NaI curve.",
                grid[0], E_next, prev, nxt,
            )
            return None
        prev = nxt

    # Guard (e) — sanity vs default NaI on the LSRM-anchor grid.
    worst_abs = 0.0
    worst_rel = 0.0
    for E_chk in _BOOTSTRAP_SANITY_GRID_keV:
        f_bs = math.sqrt(max(_eval_fwhm2_quadratic(model, E_chk), 0.0))
        f_def = math.sqrt(max(_eval_fwhm2_quadratic(
            _DEFAULT_NAI_FWHM_MODEL, E_chk), 0.0))
        delta = abs(f_bs - f_def)
        rel = delta / f_def if f_def > 0 else float("inf")
        worst_abs = max(worst_abs, delta)
        worst_rel = max(worst_rel, rel)
    if (worst_abs > _BOOTSTRAP_SANITY_ABS_MAX_keV
            or worst_rel > _BOOTSTRAP_SANITY_REL_MAX):
        logger.info(
            "F-449 bootstrap FWHM: diverges from default NaI by "
            "max %.2f keV (limit %.1f) / %.1f%% (limit %.0f%%); "
            "falling back to default NaI curve.",
            worst_abs, _BOOTSTRAP_SANITY_ABS_MAX_keV,
            100.0 * worst_rel, 100.0 * _BOOTSTRAP_SANITY_REL_MAX,
        )
        return None

    return model, "bootstrap_from_significant_peaks_sqrt", len(anchors)


def build_fwhm_model(spec) -> Tuple["FwhmModel", str]:
    """
    Build an empirical FWHM(E) model from the spectrum, returned as
    ``FwhmModel`` (callable, model-agnostic).

    Priority order:
      1. Embedded ``lsrm_peaks_table`` with >=3 PEAKS= rows -> quadratic
         fit. If the resulting fit is **pathological at low E**
         (FWHM^2(E) < ``_FWHM_PATHOLOGY_VAL_THRESHOLD_keV2`` keV^2 for any
         test energy in ``_FWHM_PATHOLOGY_TEST_ENERGIES_keV``), it is
         rejected and we fall through to step 1b / 1c via
         ``_resolve_pathology_fallback``.
      1b. Pathological quadratic -> convert the **stored LSRM sqrt(E)
          polynomial** (``StoredFwhmCalibration.coefficients``, model
          ``lsrm_fwhm_polynomial_in_E``) into FWHM^2(E) form by
          sampling-then-refitting. Preferred over the generic default
          because the operator-certified .spe polynomial reflects the
          actual detector. A ``warnings`` entry is pushed to
          ``spec.extras["fwhm_model_warnings"]``.
      1c. If the stored sqrt(E) polynomial is also unavailable or
          extrapolates below 70 % of default NaI 63x63 FWHM at E=60 keV,
          fall back to ``_DEFAULT_NAI_FWHM_MODEL``.
      2. Embedded ``lsrm_peaks_table`` with 1-2 rows -> 1-parameter
         alpha*sqrt(E) model. Same pathology guard applied.
      3. **F-449 / agent-a-math-2 (2026-06-20)** — bootstrap FWHM(E)
         from significant, isolated peaks of the spectrum itself (method
         2 of the operator-locked F-rule «FWHM любого пика»). This is the
         intended path when no trustworthy stored FWHM curve is
         available: free-σ on isolated anchors → fit FWHM²(E)=a+b·E+c·E²
         → curve becomes the lock for all downstream consumers. The
         hard-coded LSRM stored sqrt(E) polynomial is NOT used as a
         fallback here — verified empirically to overshoot the default
         NaI 63x63 curve by +17 % at 900 keV on this corpus, so a
         spectrum-specific bootstrap is methodically preferred.
      4. Default NaI 63x63 (Gamma-1S) model — last-resort fallback when
         bootstrap could not collect >= 4 isolated anchors.

    The numerical floor ``max(val, 0.01)`` in ``fwhm_keV_at_energy``
    is retained as an absolute last-resort safety net (no semantic
    meaning).

    BUG-41 fix (Wave 7, 2026-06-05): step 1b added; the old code path
    ``lsrm_peaks_table_quadratic`` -> return-as-is meant any negative-
    discriminant fit (observed on the AmTiCsEu Marinelli fixture
    val(59.5) = -15.4) propagated to identification and collapsed the
    match window to 0.3 keV. See KFI BUG-37 / BUG-41 and revalidation
    outbox ``_state/agent_a/outbox/2026-06-04_v1_22_0_AmTiCsEu_revalidation.md``
    section 5.

    Returns ``(FwhmModel, source_label)``. ``source_label`` reflects
    the branch finally chosen, so reports can show e.g.
    ``"lsrm_peaks_table_quadratic_rejected_to_default_NaI_63x63"``.

    F-452 (2026-06-21): F-160 ground-truth ветка теперь возвращает
    ``FwhmModel(kind='lsrm_poly_sqrt_E', coefficients=(c_0..c_4))`` —
    честный полином 4-й степени в √E (как в LSRM «Алгоритмических
    основах» §8.3), вместо 3-параметрической NNLS-аппроксимации в виде
    квадратичной FWHM²(E)=a+b·E+c·E². Это убирает документированный
    систематический сдвиг ±5-7 кэВ на anchor-точках. Остальные ветви
    (peaks_table, bootstrap, default) — без изменений в физической
    модели, обёрнуты в FwhmModel(kind='quad_fwhm2_in_E') для единого API.
    """
    # F-160 (2026-06-20, оператор Дмитрий: «всегда нужно делать
    # калибровку по fwhm полагаться на расчетную кривую нельзя»):
    # для эталонных спектров поверочного набора Гамма-1С — LSRM-кривая
    # FWHM(E) из references/lsrm_ground_truth/<base>/ имеет
    # АБСОЛЮТНЫЙ приоритет над bootstrap и default. Маппинг .spe→папка
    # определяется в references/lsrm_ground_truth/_index.json (явный,
    # без транслитерации). Это первый источник правды для калибровки.
    #
    # F-452 (2026-06-21): fit honest poly-4 in √E directly — это та же
    # физическая модель, что LSRM использует для своего эталонного
    # полинома (RAG-043, lsrm_act_2014 §8.3, model
    # `lsrm_fwhm_polynomial_in_E`). Unconstrained lstsq на FWHM(E) (не
    # FWHM²) даёт оптимум в L²-метрике на anchor-точках. Anchor-sanity
    # check (pred ≥ 0.5·obs И pred > 0) защищает от вырожденного fit.
    gt = _load_ground_truth_fwhm_anchors(spec)
    if gt is not None:
        Es_gt, Fs_gt, gt_base = gt
        Es_a = np.array(Es_gt, dtype=np.float64)
        Fs_a = np.array(Fs_gt, dtype=np.float64)
        # Design matrix for FWHM(E) = Σ_{k=0..4} c_k · √E^k
        z = np.sqrt(Es_a)
        A = np.vstack([z ** k for k in range(5)]).T  # shape (N, 5)
        coefs, *_ = np.linalg.lstsq(A, Fs_a, rcond=None)
        gt_model = FwhmModel(
            kind="lsrm_poly_sqrt_E",
            coefficients=tuple(float(c) for c in coefs),
        )
        # Anchor-sanity: ни в одной anchor-точке fit не должен
        # предсказывать FWHM ниже 50 % от observed (грубый sanity на
        # случай bug-а в anchor-файле) И должен быть строго положителен.
        pred_F = np.array([gt_model(float(E)) for E in Es_gt])
        rel_err_ok = bool(np.all(pred_F >= 0.5 * Fs_a) and np.all(pred_F > 0))
        if rel_err_ok:
            max_abs = float(np.max(np.abs(pred_F - Fs_a)))
            _record_fwhm_fallback_warning(
                spec,
                f"F-160/F-452 FWHM-калибровка из LSRM ground-truth: "
                f"{len(Es_gt)} anchor-точек из references/lsrm_ground_truth/"
                f"{gt_base}/fwhm_calibration_lsrm.json; lstsq-fit "
                f"FWHM(E)=Σ c_k·√E^k (deg=4) с (c0..c4)=("
                f"{coefs[0]:.4g}, {coefs[1]:.4g}, {coefs[2]:.4g}, "
                f"{coefs[3]:.4g}, {coefs[4]:.4g}); "
                f"max|ΔFWHM|={max_abs:.2f} кэВ на anchor-точках; "
                f"source='lsrm_ground_truth_reference_poly4_sqrtE'.",
            )
            return gt_model, "lsrm_ground_truth_reference_poly4_sqrtE"
        else:
            _record_fwhm_fallback_warning(
                spec,
                f"F-160/F-452 ALERT: LSRM ground-truth anchor-точки "
                f"загружены ({len(Es_gt)} шт из {gt_base}), но "
                f"poly-4 lstsq-fit отвергнут: rel_err_ok={rel_err_ok}, "
                f"min(pred_F)={float(np.min(pred_F)):.3f} кэВ. Откат на "
                f"bootstrap/lsrm_peaks_table/default.",
            )

    pks = spec.extras.get("lsrm_peaks_table") if spec.extras else None
    if pks:
        Es, Fs = [], []
        for p in pks:
            e = p.get("energy_keV"); f = p.get("fwhm_keV")
            if e and f and e > 0 and f > 0:
                Es.append(float(e)); Fs.append(float(f))
        if len(Es) >= 3:
            Es_a = np.array(Es); Fs_a = np.array(Fs)
            A = np.vstack([np.ones_like(Es_a), Es_a, Es_a**2]).T
            coefs, *_ = np.linalg.lstsq(A, Fs_a**2, rcond=None)
            quad_coefs = (
                float(coefs[0]), float(coefs[1]), float(coefs[2]),
            )
            if not _model_is_pathological(quad_coefs):
                return _wrap_quad(quad_coefs), "lsrm_peaks_table_quadratic"
            fb_coefs, fb_label = _resolve_pathology_fallback(
                spec,
                original_source="lsrm_peaks_table_quadratic",
                original_model=quad_coefs,
            )
            return _wrap_quad(fb_coefs), fb_label
        if 1 <= len(Es) <= 2:
            alpha_coefs = _fit_alpha_sqrt_E_model(Es, Fs)
            if not _model_is_pathological(alpha_coefs):
                return _wrap_quad(alpha_coefs), "lsrm_peaks_table_alpha_sqrt_E"
            fb_coefs, fb_label = _resolve_pathology_fallback(
                spec,
                original_source="lsrm_peaks_table_alpha_sqrt_E",
                original_model=alpha_coefs,
            )
            return _wrap_quad(fb_coefs), fb_label
    # F-449 / agent-a-math-2 (2026-06-20): bootstrap from significant
    # peaks of THIS spectrum. Preferred over generic default NaI 63x63
    # when no lsrm_peaks_table is available — generic-curve mismatch
    # caused the Th-232 Marinelli FWHM report skew (default_NaI_63x63
    # underestimates at 200..600 keV and overestimates at 2614 keV vs
    # LSRM ground-truth on this spectrum).
    bs = _bootstrap_fwhm_from_significant_peaks(spec)
    if bs is not None:
        bs_coefs, label, n_anchors = bs
        _record_fwhm_fallback_warning(
            spec,
            f"F-449 bootstrap: collected {n_anchors} isolated significant "
            f"peaks from spectrum; fitted FWHM^2(E)=a+b*E+c*E^2 with "
            f"(a,b,c)=({bs_coefs[0]:.4g}, {bs_coefs[1]:.4g}, "
            f"{bs_coefs[2]:.4g}); source='{label}'.",
        )
        return _wrap_quad(bs_coefs), label
    # F-160 «всегда нужно делать калибровку по FWHM, полагаться на
    # расчётную кривую нельзя» — fallback на default_NaI_63x63 здесь
    # должен СРАЗУ быть виден оператору как HIGH-severity warning, а не
    # утонуть в pipeline_notes. Содержательно ничего страшного не
    # произошло (FWHM-кривая есть), но это сигнал что bootstrap-гарды
    # отвергли все anchor-варианты — в этом случае FWHM-numbers и все
    # derived (площади, активности) на высоких энергиях занижены.
    _record_fwhm_fallback_warning(
        spec,
        "F-160 ALERT: FWHM bootstrap отвергнут всеми гардами — используется "
        "generic-кривая default_NaI_63x63. Это нарушает F-rule «всегда "
        "калибровать FWHM по спектру». Возможные причины: <6 isolated "
        "anchors, пустая bootstrap-полоса (200..800/800..1500/1500..2700), "
        "a/b out-of-bounds, monotonicity/sanity guard. Проверь WARNING-логи "
        "F-449 в этом прогоне.",
    )
    return _DEFAULT_NAI_FWHM_MODEL_OBJ, "default_NaI_63x63"


def fwhm_keV_at_energy(model, e_keV: float) -> float:
    """Evaluate FWHM(E) at energy ``e_keV``.

    F-452 (2026-06-21): accepts either a ``FwhmModel`` instance (new
    default API returned by ``build_fwhm_model``) or a legacy 3-tuple
    ``(a, b, c)`` interpreted as ``FWHM²(E) = a + b·E + c·E²`` for
    backward compatibility with internal helpers that still pass tuples.
    """
    if isinstance(model, FwhmModel):
        return model(e_keV)
    a, b, c = model
    e = max(float(e_keV), 5.0)
    val = a + b * e + c * e * e
    # BUG-43 / Wave 7 / 2026-06-05 — DOC-ONLY CORRECTION.
    #
    # The previous comment claimed: "any physical scintillator model
    # gives val >> 1 above E=5 keV". This is FALSE for an
    # ``lsrm_peaks_table_quadratic`` fit when the embedded PEAKS= table
    # is sparse or clustered. Observed pathology on the AmTiCsEu
    # Marinelli fixture (revalidation outbox section 5 BUG-41):
    #     val(E) = -174.5 + 2.67*E + 0.0004*E^2
    # which gives val(59.5) = -15.4 (negative discriminant), val(20) =
    # -121, val(40) = -67. The 0.01 floor then yields FWHM = 0.1 keV —
    # far below real NaI 63x63 FWHM at 60 keV (~12 keV) — and the
    # downstream match window collapses to 0.3 keV vs the required
    # ~36 keV. Am-241 59.54 keV characteristic line was never matched.
    #
    # **The real BUG-41 fix lives in ``build_fwhm_model``** above
    # (pathology detection + LSRM stored sqrt(E) polynomial fallback
    # + _DEFAULT_NAI_FWHM_MODEL last resort + ``fwhm_model_warnings``
    # surfaced on ``spec.extras``). By the time a model reaches this
    # function it has already been sanity-checked at the test energies
    # ``_FWHM_PATHOLOGY_TEST_ENERGIES_keV`` and either passed or been
    # replaced with a safe fallback.
    #
    # The ``max(val, 0.01)`` floor below is retained as an absolute
    # last-resort numerical safety net (e.g. ``e < 5 keV`` after clamp
    # combined with adversarial fits). It has no semantic meaning.
    #
    # See ``audit/_rag/RAG_INDEX.json`` RAG-043, Gilmore section 6.4,
    # and KFI BUG-37 / BUG-41 / BUG-43.
    return math.sqrt(max(val, 0.01))


def fwhm_model_legacy_abc(model) -> Tuple[float, float, float]:
    """Backward-compat: представить FWHM-модель как legacy 3-tuple
    ``(a, b, c)`` для формы ``FWHM²(E) = a + b·E + c·E²``.

    F-452 (2026-06-21): после введения ``FwhmModel`` (callable, kind +
    coefficients) часть reporting-слоя (``json_report._build_calibration``,
    ``plots._fwhm_keV_at``, ``interactive_html`` overlays) и legacy schema
    ``report.json:fwhm_cal.coefficients`` опираются на жёсткий
    3-параметрический контракт. Этот helper мостит legacy схему:

      * ``FwhmModel(kind='quad_fwhm2_in_E')`` — coefficients passthrough.
      * ``FwhmModel(kind='lsrm_poly_sqrt_E')`` — least-squares refit
        FWHM²(E) = a+b·E+c·E² на репрезентативной сетке (60..2614 keV).
        Вернёт лучшее 3-параметрическое приближение с документированной
        ошибкой ±5-7 keV на anchor-точках (для legacy полей; F-452
        честная poly-4 sqrt(E) хранится в самом FwhmModel и используется
        всеми callable-консьюмерами через ``fwhm_keV_at_energy``).
      * legacy 3-tuple — passthrough.

    Использовать ТОЛЬКО для legacy schema/displays; для evaluation FWHM(E)
    использовать ``fwhm_keV_at_energy(model, E)`` (model-agnostic).
    """
    if isinstance(model, FwhmModel):
        if model.kind == "quad_fwhm2_in_E":
            a, b, c = model.coefficients
            return float(a), float(b), float(c)
        # lsrm_poly_sqrt_E → refit на сетке для legacy 3-параметрической
        # schema; ошибка ±5-7 keV здесь — намеренная цена legacy contract,
        # callable-путь model(E) остаётся честным poly-4.
        Es = np.array(
            [60.0, 100.0, 200.0, 500.0, 800.0, 1200.0, 1800.0, 2614.0],
            dtype=np.float64,
        )
        Fs = np.array([model(float(E)) for E in Es], dtype=np.float64)
        A = np.vstack([np.ones_like(Es), Es, Es ** 2]).T
        coefs, *_ = np.linalg.lstsq(A, Fs ** 2, rcond=None)
        return float(coefs[0]), float(coefs[1]), float(coefs[2])
    a, b, c = model
    return float(a), float(b), float(c)


# ──────────────────────────────────────────────────────────────────
# Peak search
# ──────────────────────────────────────────────────────────────────

def _make_fwhm_at_channel(spec, fwhm_model) -> Callable[[int], float]:
    def fwhm_at_ch(ch: int) -> float:
        e = spec.channel_to_energy(ch)
        gain = abs(spec.energy_cal[1])
        return max(2.0, fwhm_keV_at_energy(fwhm_model, e) / gain)
    return fwhm_at_ch


def _run_peak_search(
    spec, fwhm_at_ch, sigma_threshold=2.5,
    *, method: str = "mariscotti",
    filter_narrow_peaks: bool = False,
    narrow_peak_fwhm_ratio: float = 0.3,
) -> Tuple[List[FoundPeak], Optional[dict]]:
    """F-129 / v1.17.7 — выбор метода поиска пиков.

    method:
      • "mariscotti"  — Mariscotti second-derivative (default, v1.0+);
      • "convolution" — matched-filter свёртка (F-124 / v1.17.6);
      • "compare"     — оба метода + сравнение в diagnostics.

    Возвращает (peaks, comparison_dict). comparison_dict не None только
    в режиме "compare" — это словарь от ``compare_peak_methods``.
    """
    counts = spec.counts.astype(float)
    method = (method or "mariscotti").lower().strip()
    if method not in ("mariscotti", "convolution", "compare"):
        method = "mariscotti"
    mar_kwargs = dict(
        counts=counts, fwhm_channels=fwhm_at_ch,
        sigma_threshold=sigma_threshold,
        min_separation_factor=0.6, edge_margin=10,
        # F-139 / v1.17.7 — отбраковка узких шумовых пиков. Прокинуто
        # из analyze_lsrm_spe(). Pipeline default OFF для back-compat с
        # synthetic-тестами; реальные демо включают через CLI флаг.
        filter_narrow_peaks=bool(filter_narrow_peaks),
        min_fwhm_ratio=float(narrow_peak_fwhm_ratio),
    )
    if method == "mariscotti":
        return mariscotti_search(**mar_kwargs), None
    # Импорт ленивый — F-124 модуль может быть исключён в кастомных сборках
    from gamma.peaks.convolution_search import (
        convolution_peak_search, compare_peak_methods,
    )
    conv_kwargs = dict(
        counts=counts, fwhm_channels=fwhm_at_ch,
        sigma_threshold=sigma_threshold,
        min_separation_factor=0.6, edge_margin=10,
    )
    if method == "convolution":
        return convolution_peak_search(**conv_kwargs), None
    # compare: запускаем оба, primary = mariscotti (для совместимости
    # downstream), comparison прикладываем в diagnostics.
    peaks_mar = mariscotti_search(**mar_kwargs)
    peaks_conv = convolution_peak_search(**conv_kwargs)
    cmp = compare_peak_methods(
        peaks_mar, peaks_conv, tolerance_channels=1.5,
    )
    cmp["primary_method"] = "mariscotti"
    cmp["secondary_method"] = "convolution"
    cmp["n_mariscotti"] = len(peaks_mar)
    cmp["n_convolution"] = len(peaks_conv)
    return peaks_mar, cmp


# ──────────────────────────────────────────────────────────────────
# Stage runner
# ──────────────────────────────────────────────────────────────────

def _run_stage(
    spec, peaks, stage_num, fwhm_at_ch, window,
    extra_candidates: Optional[List[str]] = None,
    already_detected_names: Optional[List[str]] = None,
    line_window_overrides_keV: Optional[dict] = None,
    # F-449 (agent-a-math-2 2026-06-20): energy-domain FWHM curve so
    # identify_nuclides can populate LineMatch.gauss_sigma_keV directly
    # from FWHM(E) instead of round-tripping channels→keV (avoids the
    # local-dE/dch artifact on non-linear e-cal; see identify.py).
    fwhm_keV_at_energy_fn: Optional[Callable[[float], float]] = None,
) -> StageResult:
    candidates = candidates_for_stage(stage_num)
    if extra_candidates:
        for c in extra_candidates:
            if c not in candidates:
                candidates.append(c)
    raw = identify_nuclides(
        found_peaks=peaks,
        spec=spec,
        candidate_nuclides=candidates,
        window=window,
        compute_peak_areas=True,
        fwhm_at_channel=fwhm_at_ch,
        fwhm_keV_at_energy_fn=fwhm_keV_at_energy_fn,
        line_window_overrides_keV=line_window_overrides_keV,
    )
    # Disambiguation (Lsrm Rule 1..5) — removes Na-22@511 when Tl-208 wins,
    # Cs-134@604.7 when Bi-214@609 wins, Ra-224@241 when Pb-212@239 wins,
    # etc.
    try:
        disambig = disambiguate_identifications(raw)
    except Exception as e:
        disambig = raw
        notes = [f"disambiguate failed: {e!r} — using raw identifications"]
    else:
        notes = []
    detected = list(disambig.detected_nuclides)
    rejected = list(disambig.rejected_nuclides)
    return StageResult(
        stage=stage_num,
        candidates_considered=candidates,
        detected=detected,
        rejected=rejected,
        unmatched_peaks=list(disambig.unmatched_peaks),
        notes=notes,
    )


# ──────────────────────────────────────────────────────────────────
# Residual diagnostic — decide whether Stage 2/3 is warranted
# ──────────────────────────────────────────────────────────────────

def _diagnose_residuals(
    stage1: StageResult,
    fwhm_model,
    spec,
    *,
    sigma_floor: float = 4.0,
) -> Tuple[Optional[int], str, List[ResidualClassification]]:
    """
    Decide whether to recommend escalation to Stage 2.

    F-74: each unmatched peak is first classified (XRF, chain-secondary,
    escape, sum, annihilation, edge_of_range, or true_unmatched). Only
    `true_unmatched` peaks count toward the Stage-2 recommendation. This
    prevents the "fantasizing" failure mode where escalation is suggested
    even though every residual is well-explained physics (Pb-XRF, Tl-208
    backscatter at 235, etc.).

    Returns (recommended_stage, reason, classifications).
    """
    fwhm_at = lambda e: fwhm_keV_at_energy(fwhm_model, e)
    classifications = classify_residuals(
        stage1.unmatched_peaks, spec, stage1.detected,
        fwhm_provider_keV=fwhm_at,
        sigma_floor=sigma_floor,
    )
    true_unmatched = [c for c in classifications
                      if c.label == LBL_TRUE_UNMATCHED]
    explained_n = len(classifications) - len(true_unmatched)

    if not true_unmatched:
        return None, (
            f"Все {len(classifications)} значимых пика (σ≥{sigma_floor}) "
            f"объяснены: {explained_n} как рентген.флуор./вторичные по цепи/ускользания/сумм. — "
            "эскалация в Stage 2 не требуется."
        ), classifications

    # Check if any true-unmatched peaks fall in a technogenic window
    tech_lines = [
        (661.66, "Cs-137"), (604.72, "Cs-134"), (795.86, "Cs-134"),
        (1173.23, "Co-60"), (1332.49, "Co-60"), (364.49, "I-131"),
    ]
    hits = []
    for c in true_unmatched:
        for el, name in tech_lines:
            if abs(c.peak_E_keV - el) < 0.5 * fwhm_at(el):
                hits.append((c, el, name))
                break

    if hits:
        names = sorted({h[2] for h in hits})
        peaks_desc = ", ".join(f"{c.peak_E_keV:.0f} кэВ (σ={c.sigma:.0f})"
                               for c, _, _ in hits[:5])
        return 2, (
            f"Найдено {len(true_unmatched)} истинно-неидентифицированных пиков "
            f"σ≥{sigma_floor} (объяснено {explained_n} как вторичные/рентген.флуор.), "
            f"из них {len(hits)} в окне техногенных линий "
            f"({', '.join(names)}): {peaks_desc}. "
            "Запросить подтверждение пользователя для Stage 2."
        ), classifications

    return 2, (
        f"Найдено {len(true_unmatched)} истинно-неидентифицированных пиков "
        f"σ≥{sigma_floor} (объяснено {explained_n} как вторичные/рентген.флуор.), не "
        "попадающих в окно техногенных линий. Запросить подтверждение "
        "пользователя — возможна редкая экзотика или калибровочный артефакт."
    ), classifications


# ──────────────────────────────────────────────────────────────────
# Public entry point
# ──────────────────────────────────────────────────────────────────

def analyze_lsrm_spe(
    path: str,
    *,
    detector_type: str = "NaI",
    sigma_threshold: float = 2.5,
    fwhm_window_multiple: float = 0.5,
    allow_stage2: bool = False,
    allow_stage3: bool = False,
    auto_escalate: bool = False,
    user_confirmed_stage2_nuclides: Optional[List[str]] = None,
    user_confirmed_stage3_nuclides: Optional[List[str]] = None,
    background_path: Optional[str] = None,
    # F-84 / v1.13.0 — Round 5 keywords (all default OFF for
    # backward compatibility with v1.12.0 callers)
    apply_deconvolution: bool = False,
    deconvolution_overlap_fwhm: float = 1.0,
    compute_activities: bool = False,
    sample_mass_kg: Optional[float] = None,
    # F-122 / v1.17.6 — self-attenuation для Marinelli / Дента / Петри
    sample_density_g_cm3: Optional[float] = None,
    matrix_composition: Optional[dict] = None,
    # F-129 / v1.17.7 — выбор метода поиска пиков
    peak_search_method: str = "mariscotti",
    # F-139 / v1.17.7 — отбраковка узких пиков (FWHM_meas < ratio·FWHM(E))
    # F-316 / v1.18.14: автоматически включаем default ON для NaI detector
    # (filter_narrow_peaks=None → resolved через detector_class в коде).
    filter_narrow_peaks: Optional[bool] = None,
    narrow_peak_fwhm_ratio: float = 0.3,
    # F-131 / F-135 / v1.17.7 — авто-поиск + авто-применение
    # подходящего фонового спектра. CLI default = "apply" (фон
    # вычитается ВСЕГДА при наличии подходящего кандидата) — это
    # закреплённое правило безопасности. Pipeline kwarg default
    # = "suggest" для test isolation (synthetic-тесты ожидают v1.17.6
    # behavior без auto-bg). reporting/wrapper.py подменяет default
    # на "apply" при вызове через analyze_and_report().
    # "off" — отключить; "suggest" — только записать в notes; "apply" —
    # вычесть лучшего кандидата.
    background_auto: str = "suggest",
    background_auto_max_days: int = 90,
    compute_mda: bool = False,
    mda_suite_extra_lines_keV: Optional[List[Tuple[str, float]]] = None,
    reference_datetime=None,
    # F-85 / v1.14.0 — Step 11 convenience: complete_workflow turns on
    # Steps 8/9/10 autonomous defaults per SKILL.md methodology
    complete_workflow: bool = False,
    # F-87 / v1.15.0 — Step 5β opt-in calibration refit hook. When True
    # and the seeded anchors disagree with the stored energy calibration
    # by more than `recalibration_threshold_fwhm` × FWHM at the anchor
    # energy, the bootstrap refits E(N) on the anchor channels (deg ≤ 4)
    # and re-seeds once. Default is False to preserve the v1.14.0
    # contract — the stored .spe calibration is trusted by default.
    recalibrate_on_anchor_disagreement: bool = False,
    recalibration_threshold_fwhm: float = 0.3,
    # F-167 / v1.17.10 — каноничное ID-окно ±k·FWHM(E) per [LSRM-Algo-9].
    # Когда True, Phase D matching использует `build_id_window_k_fwhm`
    # с k=1.5 для NaI / k=1.0 для HPGe вместо legacy эвристики
    # `δE₀·√(E/E_ref)`. Это даёт корректное масштабирование окна на
    # широких пиках NaI выше 1500 кэВ.
    #
    # **Default = False (opt-in) для v1.17.10** — широкое окно (3× legacy)
    # переразмечает matched-peaks: близкие линии разных нуклидов
    # (например, Ac-228 1459 и K-40 1461) перехватывают друг у друга
    # пики, что меняет выходные nuclide lists и ломает 7 snapshot-тестов.
    # В v1.17.11 default переключится на True после adjust пакетов
    # `test_filename_binding`, `test_priority_express`,
    # `test_staged_identification` под канонические ЛСРМ-результаты.
    use_lsrm_id_window: bool = False,
    # F-167 опциональный override k (для нестандартных детекторов).
    id_window_k_override: Optional[float] = None,
    # F-167 / v1.17.10 — явное указание канонического класса детектора
    # ("NaI" | "CsI" | "LaBr" | "CeBr" | "HPGe" | "CdZnTe").
    # Когда None (default), класс выводится из `detector_type` через
    # `normalize_detector_class()`. Используй, когда `detector_type`
    # имеет нестандартный формат (например, "Колибри-1М NaI 63x63" из
    # filename hints) или когда хочешь явно задать HPGe для k=1.0.
    detector_class: Optional[str] = None,
    # F-307 / v1.18.7 — pass-through opt-in флагов из v1.18.1..v1.18.4
    # в compute_activities_for_all. Активируются только при compute_activities=True.
    # Default OFF для back-compat. См. docs/archive/KNOWN_AND_FIXED_ISSUES_v1.16-v1.17.md F-294..F-297 (F-336 archive split).
    enable_tcs_correction: bool = False,
    tcs_detector_id: str = "Gamma-1S",
    enable_cutshall_self_abs: bool = False,
    cutshall_path_cm: Optional[float] = None,
    cutshall_calib_density_g_cm3: float = 1.0,
    enable_matrix_method: bool = False,
    matrix_method_energy_tolerance_keV: float = 1.0,
    # F-322 / v1.18.16 — opt-in F-96 bg-anchors injection в multiplet deconv
    # (закрывает M_th M3 chi2_red=12.94 issue — добавляет 511 annihilation
    # как constraint в clusters содержащие линии в ±60 keV окне 511).
    enable_f96_bg_anchors: bool = False,
    # F-397 / v1.18.27 — internal guard: skip recursive bg-analysis hook.
    # Когда True, не запускать вторичный analyze_lsrm_spe(bg_path) для
    # peak detection в фоне (используется когда сам вызов уже относится к
    # фону — иначе бесконечная рекурсия). Default False — top-level calls
    # анализируют фон при наличии bg path и `complete_workflow=True`.
    _skip_background_analysis: bool = False,
) -> StagedAnalysisResult:
    """
    Run the staged identification pipeline on a single LSRM .spe file.

    Parameters
    ----------
    path : str
        Path to the .spe file (LSRM SpectraLine binary header + counts).
    detector_type : str, default 'NaI'
        Detector type for the identification-window builder. Project scope
        is LSRM NaI; other types are out of scope per current methodology.
    sigma_threshold : float, default 2.5
        Mariscotti peak-search significance threshold.
    fwhm_window_multiple : float, default 0.5
        Identification-window half-width as a multiple of FWHM(E).
    allow_stage2 : bool, default False
        If True, run Stage 2 (Cs-137/Cs-134/Co-60/I-131) even if Stage 1
        explains the spectrum. Use with care — per methodology, default
        is False and Stage 2 should be opt-in by the user.
    allow_stage3 : bool, default False
        Same for Stage 3 (Na-22, Be-7, Am-241, calibration sources).
    auto_escalate : bool, default False
        If True, *automatically* escalate to Stage 2 when residual
        diagnostic suggests it. Disabled by default — the user prefers
        explicit confirmation.
    user_confirmed_stage2_nuclides : list[str], optional
        If provided, Stage 2 runs with only these candidates (subset of
        TECHNOGENIC_STAGE2). Used after the user replies to a Stage 2
        consultation prompt.
    user_confirmed_stage3_nuclides : list[str], optional
        Same for Stage 3.
    background_path : str, optional
        Path to a paired background `.spe`. When given, energy-rebinned
        background subtraction is performed (F-58) and downstream
        activity / MDA results are flagged as `from_bg_subtracted=True`.
    apply_deconvolution : bool, default False
        **Round 5 (F-84)**: when True, the multiplet deconvolution
        post-pass `gamma.peaks.deconvolve.apply_multiplet_deconvolution`
        runs after disambiguation. For each cluster of overlapping
        identified lines (separation ≤ `deconvolution_overlap_fwhm`·FWHM),
        component areas are recovered by constrained NNLS fit at fixed
        positions/widths. Updated `LineMatch.peak_area` carries the
        deconvolved value with `peak_area_source="deconvolved"`. The
        canonical use case is the **600–680 keV multiplet** on Th-rich
        NaI spectra (Bi-214 609.31 / Cs-137 661.66 / Cs-134 604.72/795.86).
    deconvolution_overlap_fwhm : float, default 1.0
        Overlap threshold (in units of FWHM at the cluster centre) used
        by `find_multiplet_regions` to decide what counts as a multiplet.
        Default 1.0 (FWHM-touching pairs).
    compute_activities : bool, default False
        **Round 5 (F-84)**: when True AND an efficiency curve is loaded,
        per-nuclide activity is computed via `compute_activities_for_all`
        (Lsrm §8.4 weighted-mean across matched lines). Returned in
        `result.activities`. Decay correction is applied automatically
        when `reference_datetime` is supplied.
    sample_mass_kg : float, optional
        **Round 5 (F-84)**: when given AND `compute_activities=True`,
        per-nuclide specific activity is reported in `Bq/kg`
        (`result.specific_activities_Bq_per_kg`). Use 0.0 or None to
        report only absolute Bq.
    compute_mda : bool, default False
        **Round 5 (F-84)**: when True AND an efficiency curve is loaded,
        per-line MDA per ISO 11929 / Lsrm §6.3 is computed for (i)
        every detected nuclide's matched lines and (ii) the standard
        ЕРН/technogenic suite. Returned in `result.mda_per_line` as
        `{(nuclide, E_keV): MdaResult}`.
    mda_suite_extra_lines_keV : list[(nuclide, E_keV)], optional
        Extra lines to include in the MDA suite beyond the default
        Cs-137 661.66 / Co-60 1173+1332 / K-40 1460.82 / Bi-214 1764 /
        Tl-208 2614 / Ac-228 911 list.
    reference_datetime : datetime, optional
        Source certificate / reference epoch for decay correction.
        Forwarded to `compute_activity` when activities are requested.
    complete_workflow : bool, default False
        **F-85 / v1.14.0 (Step 11)**: convenience umbrella that
        autonomously enables Steps 8 (`apply_deconvolution=True`,
        `deconvolution_overlap_fwhm=3.0` to capture the 600–680 keV
        Bi-214/Cs-137/Cs-134 cluster on Th-rich NaI), 9
        (`compute_activities=True`, `compute_mda=True`) and 10
        (already always-on via residual_classifier and the v1.7.16
        secondary_peaks_v2 catalog). The flag respects per-argument
        overrides — if a caller passes `compute_mda=False` together
        with `complete_workflow=True`, MDA stays off. Default OFF for
        backward compatibility.

    Returns
    -------
    StagedAnalysisResult
        Complete analysis report including per-stage breakdown plus
        — when requested — multiplet deconvolution, per-nuclide
        activities (Bq and Bq/kg), and per-line MDA.
    """
    # F-85 / v1.14.0 — Step 11 umbrella. When `complete_workflow=True`,
    # apply Round-5 hooks autonomously (Steps 8/9/10 per SKILL.md). The
    # individual flags can be overridden by passing them explicitly to
    # True before this point (they'd already be True), but cannot be
    # turned OFF by the umbrella — that's the point of the umbrella.
    if complete_workflow:
        apply_deconvolution = True
        # Co-60 1173/1332 doublet, K-40 1460 + Compton, and the Th-chain
        # 600-680 keV cluster all need overlap_threshold > 1·FWHM to be
        # treated as multiplets on NaI 63×63 — push to 3.0 only if the
        # caller did not override.
        if deconvolution_overlap_fwhm == 1.0:
            deconvolution_overlap_fwhm = 3.0
        compute_activities = True
        compute_mda = True

    spec = read_spectrum(path)

    # BUG-48 / v1.22.x — hard calibration self-consistency gate.
    # Catches non-monotone energy axes, broken FWHM models and out-of-band
    # resolution **before** any peak search runs. Result is stored in
    # ``spec.extras['calibration_gate']`` for downstream consumers
    # (report.json warnings channel, CLI summary). Does not raise — the
    # gate ALWAYS surfaces; ``passed=False`` means downstream identification
    # SHOULD treat the spectrum as suspect, but the legacy contract is
    # preserved: ID still runs so users see what the broken cal produced.
    # See scripts/gamma/calibration/calibration_gate.py for criteria.
    try:
        from gamma.calibration.calibration_gate import evaluate_calibration_gate
        _cal_gate_result = evaluate_calibration_gate(spec)
        if spec.extras is None:
            spec.extras = {}
        spec.extras["calibration_gate"] = _cal_gate_result.as_dict()
    except Exception as _e:  # pragma: no cover - defensive
        if spec.extras is None:
            spec.extras = {}
        spec.extras["calibration_gate"] = {
            "passed": True,  # do not block on gate internals failure
            "hard_failures": [],
            "soft_warnings": [],
            "criteria_evaluated": [],
            "detector_class_hint": "unknown",
            "E_low_keV": None,
            "E_high_keV": None,
            "reason": f"calibration_gate evaluation error (skipped): {_e!r}",
        }

    # F-130 / v1.17.7 — auto-detect sample density from .spe metadata.
    # Если CLI/wrapper не передали явное значение sample_density_g_cm3,
    # пробуем взять из spec.extras["lsrm_sample_density_g_cm3"]
    # (заполняется reader-ом из MATERIAL.Ro / SAMPLEMASS÷SAMPLEVOLUME).
    auto_density_source = None
    if sample_density_g_cm3 is None:
        try:
            auto_rho = spec.extras.get("lsrm_sample_density_g_cm3") if spec.extras else None
        except AttributeError:
            auto_rho = None
        if auto_rho is not None:
            try:
                sample_density_g_cm3 = float(auto_rho)
                auto_density_source = (
                    spec.extras.get("lsrm_density_source", "auto")
                    if spec.extras else "auto"
                )
            except (TypeError, ValueError):
                sample_density_g_cm3 = None

    # F-140 / v1.17.7 — auto-detect sample mass (kg) from .spe metadata.
    # SAMPLEMASS поле LSRM хранит в граммах; reader делит на 1000.
    # Активируется только когда CLI --sample-mass-kg не задан явно.
    auto_mass_source = None
    if sample_mass_kg is None:
        try:
            auto_m = spec.extras.get("lsrm_sample_mass_kg") if spec.extras else None
        except AttributeError:
            auto_m = None
        if auto_m is not None:
            try:
                sample_mass_kg = float(auto_m)
                auto_mass_source = (
                    spec.extras.get("lsrm_mass_source", "auto")
                    if spec.extras else "auto"
                )
            except (TypeError, ValueError):
                sample_mass_kg = None

    # F-378 / v1.18.25 — strict mass mismatch check.
    # Triggers only when CLI mass won (auto_mass_source is None) — see
    # `check_sample_mass_mismatch` for the rule. Mirror в stderr для
    # видимости при интерактивном запуске.
    mass_mismatch_note = None
    if auto_mass_source is None:
        mass_mismatch_note = check_sample_mass_mismatch(
            cli_mass_kg=sample_mass_kg,
            spec_extras=spec.extras if hasattr(spec, "extras") else None,
        )
        if mass_mismatch_note:
            import sys as _sys
            _sys.stderr.write(mass_mismatch_note + "\n")

    # ─── F-131 / v1.17.7 — авто-поиск подходящего фонового спектра ───
    # Если пользователь не задал --background-path явно, у spec нет
    # embedded BG и режим background_auto != "off", запускаем
    # эвристический поиск в той же папке + типовых местах хранения
    # фонов (data/averaged_backgrounds, *Фон*/ подпапки детектора).
    auto_bg_candidates_list = []
    auto_bg_applied_note = None
    auto_bg_suggest_note = None
    if (background_auto in ("suggest", "apply")
            and not background_path
            and getattr(spec, "background_embedded", None) is None):
        try:
            from gamma.io.background_search import (
                find_background_candidates,
                render_suggestion_note,
                render_applied_note,
            )
            auto_bg_candidates_list = find_background_candidates(
                spec, path, max_days_apart=int(background_auto_max_days),
            )
            if auto_bg_candidates_list:
                best = auto_bg_candidates_list[0]
                if background_auto == "apply":
                    # Подставляем background_path так, чтобы F-58 ниже
                    # вычел этот файл штатным путём.
                    background_path = str(best.path)
                    auto_bg_applied_note = render_applied_note(best)
                else:
                    # suggest: только нарративная заметка, gross-расчёт
                    # продолжится.
                    auto_bg_suggest_note = render_suggestion_note(best)
        except Exception:
            # F-131 — никогда не падает: при любой ошибке продолжаем
            # на gross-спектре.
            auto_bg_candidates_list = []
            auto_bg_applied_note = None
            auto_bg_suggest_note = None

    # ─── F-58: optional background subtraction ───
    # If a paired background path is provided, subtract it (re-binned to
    # the sample's energy grid, scaled by live-time ratio) and continue
    # the pipeline on the NET counts. The original `spec.counts` is
    # replaced in-place; the subtraction result is preserved in the
    # output for diagnostics.
    bg_sub_result: Optional[BackgroundSubtractionResult] = None
    # F-332 / v1.18.18.5 — captured for 4-way chart toggle in HTML report.
    f332_gross_counts: Optional[np.ndarray] = None
    f332_bg_on_grid: Optional[np.ndarray] = None
    f332_bg_live_time: Optional[float] = None
    f332_bg_scale: Optional[float] = None
    # F-243 / v1.18.29 — BG control pre-check (sum_Y, t_live, rate ratio).
    # Warning-only: failures are surfaced via pipeline notes; pipeline does
    # NOT abort.
    f243_bg_warning: Optional[str] = None
    bg_spec = None  # F-QC-01: keep in scope for post-peak z-test
    if background_path:
        bg_spec = read_spectrum(background_path)
        # F-332: stash gross counts BEFORE replacement.
        f332_gross_counts = np.array(spec.counts, dtype=np.float64).copy()
        f332_bg_live_time = float(getattr(bg_spec, "live_time", 0.0) or 0.0)
        # ─── F-243 BG control gate (non-blocking) ───
        # Hybrid layer per AGENT_A_BRIEF_F-243.md (Agent D decision 2026-06-02):
        # use the high-level sample-relative `check_bg_quality(sample, bg)` —
        # it evaluates flux_drift / sum_y_stat / live_time_min with thresholds
        # tied to *this* sample (Gate 3: bg.live_time ≥ 0.5 × sample.live_time).
        # The low-level numeric kernel `validate_background` remains imported
        # above for unit-tests and future absolute-threshold consumers.
        try:
            _bg_check = _f243_check_bg_quality(spec, bg_spec)
            if not _bg_check.passed:
                f243_bg_warning = (
                    "F-243 BG control: " + "; ".join(_bg_check.notes)
                )
        except Exception:
            # Never block pipeline on a sanity-check internal failure.
            f243_bg_warning = None
        bg_sub_result = subtract_background(spec, bg_spec)
        f332_bg_on_grid = np.array(
            bg_sub_result.background_counts_on_sample_grid, dtype=np.float64
        ).copy()
        # F-451 B1 (operator-locked 2026-06-22): scale_factor is the
        # LEGACY t_s/t_bg ratio. Downstream spec keeps **sample-scale**
        # counts (net = sample − bg·(t_s/t_bg), legacy "bg up to sample")
        # and the original `spec.live_time` / `spec.real_time` — so that
        # peak detection, FWHM bootstrap, NNLS multiplet deconvolution
        # and activity formulas continue to operate at the full sample
        # count statistics. The F-451 "к меньшему" direction
        # (`applied_scale`, `scale_direction`, `effective_live_time`)
        # lives entirely inside `BackgroundSubtractionResult` for σ
        # propagation and audit/diagnostics — NOT inside `spec`. The
        # cps-invariant net_cps = sample_rate − bg_rate holds identically
        # in either scale.
        f332_bg_scale = float(bg_sub_result.scale_factor)
        _bg_on_grid_arr = np.asarray(
            bg_sub_result.background_counts_on_sample_grid, dtype=np.float64
        )
        _net_in_sample_scale = (
            np.asarray(f332_gross_counts, dtype=np.float64)
            - _bg_on_grid_arr * f332_bg_scale
        )
        spec.counts = np.maximum(_net_in_sample_scale, 0.0)
        # spec.live_time / spec.real_time intentionally untouched (B1).
        if hasattr(spec, "extras") and isinstance(spec.extras, dict):
            spec.extras["background_subtracted"] = True
            # F-451 audit-only fields: direction + applied_scale +
            # effective_live_time of the σ-propagation. Downstream
            # computations do NOT depend on these; they exist for
            # reporting/diagnostics only.
            spec.extras["background_subtraction_applied_scale"] = float(
                bg_sub_result.applied_scale
            )
            spec.extras["background_subtraction_scale_direction"] = str(
                bg_sub_result.scale_direction
            )
            spec.extras["background_subtraction_effective_live_time"] = float(
                bg_sub_result.effective_live_time
            )

    # Build FWHM model
    fwhm_model, fwhm_src = build_fwhm_model(spec)
    fwhm_at_ch = _make_fwhm_at_channel(spec, fwhm_model)
    fwhm_661 = fwhm_keV_at_energy(fwhm_model, 661.66)

    # Identification window
    # F-167 / v1.17.10 — canonical ±k·FWHM(E) per [LSRM-Algo-9]
    # вместо legacy эвристики δE₀·√(E/E_ref). Переключается флагом
    # `use_lsrm_id_window` (default True). Когда False — старый путь
    # для обратной совместимости.
    if use_lsrm_id_window:
        # F-167: явный detector_class имеет приоритет над эвристической
        # нормализацией из detector_type. Порядок:
        #   detector_class > normalize_detector_class(detector_type) > "NaI"
        det_class = detector_class or normalize_detector_class(detector_type)
        window = build_id_window_k_fwhm(
            det_class,
            fwhm_provider_keV=lambda E: fwhm_keV_at_energy(fwhm_model, E),
            k_override=id_window_k_override,
        )
        try:
            from gamma.utils.logger import get_logger
            get_logger(__name__).info(
                "F-167 ID-window: scaling=k_fwhm | detector=%s | k=%.2f | "
                "FWHM@661=%.2f keV | window@661=%.2f keV",
                det_class, window.k_fwhm or 0.0, fwhm_661,
                window.window_keV(REFERENCE_ENERGY_KEV),
            )
        except Exception:
            pass
    else:
        window = identification_window_from_fwhm(
            fwhm_661, detector_type=detector_type,
            fwhm_multiple=fwhm_window_multiple,
        )

    # F-129 / v1.17.7 — peak search (Mariscotti default, convolution opt-in)
    # F-139 / v1.17.7 — filter_narrow_peaks прокидывается из CLI/kwarg
    # F-316 / v1.18.14 (REVERTED): попытка default-ON для NaI отсекала
    # legitimate Bi-214 609/1764 keV pair (regression test caught).
    # Filter остаётся **opt-in** (default None → False). Operators
    # могут включить через CLI флаг --filter-narrow-peaks. Для default-ON
    # перехода необходимо tuning narrow_peak_fwhm_ratio для NaI с учётом
    # Mariscotti FWHM bias на широких NaI peaks (Bi-214 1764 имеет
    # measured FWHM ~ 50 keV vs expected 60 keV — ratio 0.83, ниже
    # default threshold 0.3 не должен отсекать, но текущая логика
    # _run_peak_search видимо более агрессивна. См. v1.18.15+).
    _filter_resolved = (
        False if filter_narrow_peaks is None
        else bool(filter_narrow_peaks)
    )
    peaks, peak_method_compare = _run_peak_search(
        spec, fwhm_at_ch, sigma_threshold,
        method=peak_search_method,
        filter_narrow_peaks=_filter_resolved,
        narrow_peak_fwhm_ratio=narrow_peak_fwhm_ratio,
    )

    # ── F-QC-01 / v1.19.1 — per-peak Poisson |z|-test (BUG-35 / RAG-022) ──
    # Run AFTER peak search so we know which energies to test.
    # Uses f332_gross_counts (sample counts BEFORE bg subtraction) and
    # bg_spec.counts (reference background), with a ±1 FWHM ROI per peak.
    # Warning-only: never blocks pipeline. Results stored in bg_quality_check.
    _bg_z_check_result: Optional[dict] = None
    if (
        bg_spec is not None
        and f332_gross_counts is not None
        and peaks
        and getattr(spec, "energy_cal", None)
        and getattr(bg_spec, "energy_cal", None)
    ):
        try:
            _z_roi_entries = []
            for _pk in peaks:
                try:
                    _e_center = float(spec.channel_to_energy(_pk.channel))
                    # ROI width = ±1 FWHM (keV) around peak center.
                    _fwhm_keV = fwhm_keV_at_energy(fwhm_model, _e_center)
                    _e_lo = max(0.0, _e_center - _fwhm_keV)
                    _e_hi = _e_center + _fwhm_keV
                    # Sum counts within ROI from gross sample and bg.
                    from gamma.io.bg_control import _roi_sum_counts as _z_roi_sum
                    _B1 = _z_roi_sum(f332_gross_counts, spec.energy_cal, _e_lo, _e_hi)
                    _B2 = _z_roi_sum(bg_spec.counts, bg_spec.energy_cal, _e_lo, _e_hi)
                    _zr = _bg_z_test(_B1, _B2)
                    _z_roi_entries.append({
                        "peak_energy_keV": round(_e_center, 2),
                        "e_lo_keV": round(_e_lo, 2),
                        "e_hi_keV": round(_e_hi, 2),
                        "z": float("nan") if _zr.z != _zr.z else round(_zr.z, 3),
                        "abs_z": float("nan") if _zr.abs_z != _zr.abs_z else round(_zr.abs_z, 3),
                        "tier": _zr.tier,
                        "passed": bool(_zr.passed),
                        "B1": int(_zr.B1),
                        "B2": int(_zr.B2),
                        "note": _zr.note,
                    })
                except Exception:
                    continue
            if _z_roi_entries:
                _n_pass = sum(1 for e in _z_roi_entries if e["passed"])
                _n_fail = len(_z_roi_entries) - _n_pass
                _bg_z_check_result = {
                    "n_peaks_tested": len(_z_roi_entries),
                    "n_passed": _n_pass,
                    "n_failed": _n_fail,
                    "overall_passed": _n_fail == 0,
                    "peak_z_roi": _z_roi_entries,
                }
        except Exception:
            # Never block pipeline on QC diagnostic failure.
            _bg_z_check_result = None

    # Filename hints + canonical resolution (F-78). Header values from
    # the .spe file (`spec.geometry`, `spec.detector_id`) are also
    # passed through the alias registry — these are the most reliable
    # source since they come from the laboratory's instrument config,
    # not from arbitrary filenames.
    ft = getattr(spec, "filename_tokens", None) or {}
    sample_type = ft.get("sample_type", "")
    geometry_hint = ft.get("geometry", "") or (spec.geometry or "")
    detector_hint = (spec.detector_id or "")
    is_background = bool(ft.get("is_background_hint", False))
    # F-89b / v1.15.2 — canonical isotope labels from filename, used as
    # the binding hypothesis for Step-7 identification per SKILL.md
    # §7A.1 ("Nuclides suggested by filename/metadata — highest
    # priority"). chains_claimed_from_filename is used by F-89d to
    # decide whether competing chains (e.g. U-238 on a Th-232 source)
    # should be suppressed.
    filename_isotope_hints: List[str] = list(ft.get("isotope_hints", []) or [])
    filename_chains_claimed: List[str] = list(ft.get("chains_claimed", []) or [])

    # BUG-40 — track which raw string drove the winning canonicalization,
    # so the warning emitter can flag Cyrillic → Latin homoglyph collisions
    # (e.g. CONFIGNAME "Гамма-1С" → canonical "Gamma-1S") on profile
    # fallback. ``detector_canon_source_raw`` is the original (possibly
    # Cyrillic) text whose canonicalize() round-trip produced the winning
    # ``detector_canon``. Empty when canonicalization fell through to
    # ``ft.get("detector_canonical")`` (pre-canonicalised by filename
    # hints) or to an exception.
    detector_canon_source_raw: str = ""
    try:
        from gamma.data.aliases import canonicalize
        sample_type_canon = (ft.get("sample_type_canonical") or "") \
            or (canonicalize("sample_type", sample_type) or "")
        geometry_canon = (ft.get("geometry_canonical") or "") \
            or (canonicalize("geometry", spec.geometry or "") or "") \
            or (canonicalize("geometry", geometry_hint) or "")
        # Track which raw string wins for detector canonicalization.
        _ft_det = ft.get("detector_canonical") or ""
        _from_hint = canonicalize("detector", detector_hint) or "" \
            if not _ft_det else ""
        _from_type = canonicalize("detector", detector_type) or "" \
            if not _ft_det and not _from_hint else ""
        detector_canon = _ft_det or _from_hint or _from_type
        if _ft_det:
            # Pre-canonicalised — original raw not tracked through this path.
            detector_canon_source_raw = ""
        elif _from_hint:
            detector_canon_source_raw = str(detector_hint or "")
        elif _from_type:
            detector_canon_source_raw = str(detector_type or "")
        # BUG-39 / BUG-40 — also probe the LSRM CONFIGNAME header.
        # The CONFIGNAME field carries the COMPLEX identifier (e.g.
        # "Гамма-1С №SN-02"), whereas DETECTOR carries only the
        # head-model serial (often shared between physically distinct
        # complexes). When CONFIGNAME yields a more-specific canonical
        # than the head-model match (e.g. Gamma-1S vs the generic
        # Gamma-1S head), prefer CONFIGNAME.
        lsrm_config = ""
        try:
            lsrm_config = str(spec.extras.get("lsrm_config", "")) if spec.extras else ""
        except Exception:
            lsrm_config = ""
        if lsrm_config:
            config_canon = canonicalize("detector", lsrm_config) or ""
            if config_canon and config_canon != detector_canon:
                # CONFIGNAME wins — that's the complex-level identity.
                detector_canon = config_canon
                detector_canon_source_raw = lsrm_config
    except Exception:
        sample_type_canon = geometry_canon = detector_canon = ""
        detector_canon_source_raw = ""

    # BUG-39 — detect silent fallback to a substitute detector profile.
    # BUG-40 — augment the fallback record with the original raw string
    # (when available) and a flag indicating Cyrillic → Latin homoglyph
    # collision. Downstream report builders read these fields to emit
    # the operator-facing structured warning.
    try:
        from gamma.detectors.profile import detect_silent_fallback
        from gamma.data.aliases import cyrillic_to_latin_collision
        _fallback = detect_silent_fallback(detector_canon)
        detector_fallback_dict = _fallback.as_dict()
        detector_fallback_dict["original_raw"] = detector_canon_source_raw
        detector_fallback_dict["cyrillic_to_latin_collision"] = bool(
            cyrillic_to_latin_collision(detector_canon_source_raw, detector_canon)
        )
    except Exception:
        detector_fallback_dict = None

    def _fwhm_provider_keV(e: float) -> float:
        return fwhm_keV_at_energy(fwhm_model, e)

    # ══════════════════════════════════════════════════════════════
    # Step 5α (F-87 / v1.15.0) — anchor seeding for calibration
    # ──────────────────────────────────────────────────────────────
    # Per SKILL.md Step 5 ("use invariant peak patterns: Co-60
    # doublet, K-40 1460.82, Bi-214 series, Tl-208 2614.5, Cs-137
    # 661.66, Pb XRF triplet, LaBr3 self-pattern"), the express
    # heuristic — F-79 anchor-rank + F-80 express-patterns — is the
    # canonical way to obtain anchor identifications before any
    # calibration refit. Until v1.14.0 these passes lived inside
    # Stage-1 identification; v1.15.0 reorders them to Step 5α to
    # match the methodology. F-81 (7-line ЕРН check, Lsrm §9) remains
    # a separate Pass C — that is the *final calibration
    # verification*, not a seeding step.
    seed_mode = "background" if is_background else "sample"
    seed = seed_calibration_anchors(
        peaks, spec,
        mode=seed_mode,
        fwhm_provider_keV=_fwhm_provider_keV,
        window_fwhm_multiple=fwhm_window_multiple,
    )
    anchor_matches = seed.anchor_matches
    pattern_confirmations = seed.pattern_confirmations

    # ──────────────────────────────────────────────────────────────
    # Step 5β (F-87c/d / v1.15.0) — opt-in calibration refit
    # ──────────────────────────────────────────────────────────────
    # When `recalibrate_on_anchor_disagreement=True` and the seeded
    # anchors disagree with the stored energy calibration by more
    # than `recalibration_threshold_fwhm` × FWHM, refit E(N) on the
    # anchor channels (degree ≤ 4) and re-seed once with the new
    # calibration. Default is False to preserve the v1.14.0 contract.
    recalibration_diag: dict = {
        "attempted": False,
        "applied": False,
        "old_residual_max_keV": None,
        "new_residual_max_keV": None,
        "old_energy_cal": None,
        "new_energy_cal": None,
        "n_anchors_used": 0,
    }
    # F-453 (BUG-38 follow-up, 2026-06-23): auto-trigger when residuals
    # are clearly out-of-tolerance even without explicit opt-in. Threshold
    # 0.5·FWHM (вдвое выше standard 0.3·FWHM) — "однозначный" drift, выше
    # шума singleton-fitting. Closes BUG-38 on AmTiCsEu Marinelli where
    # F-145 multiplet self-cal is silent (n_multiplets_seen=0, no Th/U
    # forced_clusters chains). Explicit kwarg=True still wins.
    f453_auto_diag: dict = {"fired": False, "reason": "kwarg=True override"}
    if not recalibrate_on_anchor_disagreement:
        from gamma.calibration.anchor_recalibration import (
            should_auto_recalibrate as _f453_should_auto,
        )
        f453_fire, f453_auto_diag = _f453_should_auto(
            anchor_matches,
            fwhm_provider_keV=_fwhm_provider_keV,
            drift_frac_threshold=0.5,
            min_anchors=3,
        )
    else:
        f453_fire = True
    recalibration_diag["f453_auto_trigger"] = f453_auto_diag
    if recalibrate_on_anchor_disagreement or f453_fire:
        try:
            from gamma.calibration.anchor_recalibration import (
                recalibrate_energy_if_anchors_disagree,
            )
            new_cal, recalibration_diag = recalibrate_energy_if_anchors_disagree(
                spec, anchor_matches,
                threshold_fraction_of_fwhm=recalibration_threshold_fwhm,
                fwhm_provider_keV=_fwhm_provider_keV,
            )
            recalibration_diag["f453_auto_trigger"] = f453_auto_diag
            if new_cal is not None:
                spec.energy_cal = tuple(new_cal)
                # Re-seed once with the refitted calibration.
                seed = seed_calibration_anchors(
                    peaks, spec,
                    mode=seed_mode,
                    fwhm_provider_keV=_fwhm_provider_keV,
                    window_fwhm_multiple=fwhm_window_multiple,
                )
                anchor_matches = seed.anchor_matches
                pattern_confirmations = seed.pattern_confirmations
        except ImportError:
            # Module not yet present — leave stored calibration in place.
            pass

    # ══════════════════════════════════════════════════════════════
    # Step 5α′ (F-88 / v1.15.1) — user-priority order + chain dominance
    # ──────────────────────────────────────────────────────────────
    # Per user methodology (2026-05-29), the express anchors must be
    # evaluated in a specific priority order:
    #   1. 2615 Tl-208      — trump card; locks Th-232 chain
    #   2. 1461 K-40        — overlap warning if Th-dominant (Ac-228 1459)
    #   3.  662 Cs-137
    #   4. 1173+1332 Co-60  (paired)
    #   5.  609+1764 Bi-214 (paired, U-238 chain)
    #   6.   59.5 Am-241
    #
    # Chain dominance is computed from the seeded anchors + patterns
    # and HARD-PASSED to Step 7 identification: when th232/u238 is
    # True, the corresponding chain proxies (Tl-208/Pb-212/Ac-228 or
    # Bi-214/Pb-214/Pb-210/Ra-226) become strong-prior candidates.
    priority_findings_out = derive_priority_findings(
        anchor_matches, pattern_confirmations,
    )
    # F-89d — pass filename chain-claim set so derive_chain_dominance
    # can apply the suppression rule (Th-only filename ⇒ U-238 needs
    # the Bi-214 quartet, not just the Ra-pair).
    chain_dominance_out = derive_chain_dominance(
        anchor_matches, pattern_confirmations,
        filename_chains_claimed=set(filename_chains_claimed) or None,
    )

    # F-88 K-40 / Ac-228 overlap warning. On NaI 63×63 (FWHM ~85 keV
    # at 1460), the K-40 1460.82 line and the Ac-228 1459.20 line
    # (I = 0.85%) are unresolvable as a doublet. When Th-232 is
    # dominant AND K-40 is seen via the priority signal, the K-40
    # peak area is contaminated by Ac-228 — flag it.
    k40_overlap_warning_out = False
    for pf in priority_findings_out:
        if pf.signal.order == 2 and pf.matched and chain_dominance_out.th232:
            k40_overlap_warning_out = True
            break

    # ══════════════════════════════════════════════════════════════
    # Step 5γ (F-81 / canonical) — 7-line ЕРН calibration verification
    # ──────────────────────────────────────────────────────────────
    # Per Lsrm methodology §9, this is the canonical final calibration
    # check on background and ЕРН-rich spectra. Runs always but
    # mode-tagged in analysis_mode below.
    seven_line = run_seven_line_check(
        peaks, spec,
        fwhm_provider_keV=_fwhm_provider_keV,
        window_fwhm_multiple=fwhm_window_multiple,
    )

    # ─── F-123 / v1.17.6 — построить line-specific window overrides ───
    # При доминантной цепочке Th-232 расширяем окно для Pb-212 238 кэВ
    # до ±2.5·FWHM (вместо стандартных ±0.5·FWHM_at_661 ≈ ±28 кэВ при
    # FWHM_661 ≈ 56 кэВ). На NaI 63×63 FWHM(238) ≈ 25 кэВ → окно ≈
    # ±62 кэВ. Это уверенно ловит Pb-212 в Th-цепи, где он часто
    # промахивается из-за наложения с Pb-XR 73-90 кэВ и Th-228 84.4.
    # Аналогично для U-238: Pb-214 295 / 352 — расширяем окно при
    # u238-dominance (в дополнение к F-121 forced multiplets).
    line_window_overrides_keV: dict = {}
    if chain_dominance_out is not None:
        if chain_dominance_out.th232:
            f238 = _fwhm_provider_keV(238.63)
            line_window_overrides_keV[("Pb-212", 238.63)] = float(2.5 * f238)
        if chain_dominance_out.u238:
            for E_pb in (295.22, 351.93):
                f = _fwhm_provider_keV(E_pb)
                line_window_overrides_keV[("Pb-214", round(E_pb, 2))] = float(2.5 * f)

    # ─── Stage 1: ЕРН + filename-hinted isotopes ───
    # F-89e / v1.15.2 — SKILL.md §7A.1 mandates that nuclides suggested
    # by filename/metadata are the highest-priority candidates. Adding
    # them to Stage 1 means a `Cs137_*.spe` fixture gets Cs-137 into
    # the candidate list without requiring an explicit allow_stage2=True
    # from the caller. The cross-stage candidate de-duplication in
    # _run_stage handles overlaps with the ЕРН list (Tl-208 / K-40 /
    # Bi-214 / Ac-228 / Pb-212 / Pb-214 / Cs-137 already in Stage 1
    # for some lists; for others — like Cs-137 in technogenic Stage 2 —
    # the filename hint promotes them into Stage 1.)
    # F-449 (agent-a-math-2 2026-06-20): energy-domain FWHM curve so
    # identify_nuclides populates LineMatch.gauss_sigma_keV directly from
    # the FWHM(E) curve, bypassing the legacy round-trip channels→keV
    # (avoids the local-dE/dch artifact on non-linear e-cal — Th-232
    # Marinelli 2614 singlet inflated 107.84 → 116.63 keV; identify.py).
    _fwhm_keV_fn_for_id = lambda E: fwhm_keV_at_energy(fwhm_model, E)
    stage1 = _run_stage(
        spec, peaks, 1, fwhm_at_ch, window,
        extra_candidates=filename_isotope_hints or None,
        line_window_overrides_keV=(line_window_overrides_keV or None),
        fwhm_keV_at_energy_fn=_fwhm_keV_fn_for_id,
    )

    stages = [stage1]
    final_detected = list(stage1.detected)
    final_unmatched = stage1.unmatched_peaks
    already_names = [n.nuclide for n in stage1.detected]

    # ─── Residual diagnostic (F-74) ───
    rec_stage, rec_reason, residual_cls = _diagnose_residuals(
        stage1, fwhm_model, spec
    )

    # ─── F-60 CI-gating with cross-promotion from anchors/patterns ───
    anchor_nuc = {m.anchor.nuclide for m in anchor_matches
                  if m.anchor.nuclide
                  and not m.partner_required_but_missing
                  and m.anchor.rank <= 7}
    pattern_nuc = set()
    for pc in pattern_confirmations:
        if pc.confirmed and pc.pattern.nuclide:
            # Express patterns name either a specific nuclide or "Th-232 chain"
            if pc.pattern.nuclide in {"Co-60", "Bi-214", "Pb-214",
                                     "Cs-137", "Eu-152"}:
                pattern_nuc.add(pc.pattern.nuclide)
            # Chain patterns promote multiple proxies
            elif pc.pattern.nuclide == "Th-232 chain":
                pattern_nuc.update({"Tl-208", "Ac-228", "Pb-212"})

    # F-88 hard-pass — chain dominance overrides per-line filtering.
    # User methodology: "Данные о наличии тория должны жёстко
    # передаваться на этап идентификации пиков." When the trump-card
    # rule fires (Tl-208 2614 alone at σ≥5, even without partner
    # anchors), all Th-chain proxies become confirmed candidates for
    # CI-gating. Same for U-238 chain via Bi-214 Ra-pair.
    if chain_dominance_out.th232:
        pattern_nuc.update(TH232_PROXY_NUCLIDES)
    if chain_dominance_out.u238:
        pattern_nuc.update(U238_PROXY_NUCLIDES)

    ci_gating = gate_identifications(
        final_detected,
        anchor_confirmed_nuclides=anchor_nuc,
        pattern_confirmed_nuclides=pattern_nuc,
    )

    # ─── F-61 Completeness (DC %) ───
    # Build (E, area) tuples for residuals that ended up true_unmatched
    # — these are the only ones that count as "unidentified dose".
    unmatched_for_dc: List[Tuple[float, float]] = []
    for rc in residual_cls:
        if rc.label == LBL_TRUE_UNMATCHED:
            # find the matching peak for area
            for p in stage1.unmatched_peaks:
                if abs(spec.channel_to_energy(p.channel) - rc.peak_E_keV) < 0.5:
                    unmatched_for_dc.append(
                        (rc.peak_E_keV, float(p.area_estimate or 0))
                    )
                    break
    completeness_result = compute_completeness(
        detected_nuclides=final_detected,
        unmatched_peaks_E_area=unmatched_for_dc,
        efficiency_at=None,    # F-57 will wire this when avail
    )

    # ─── F-57 auto-load efficiency curve by canonical geometry ───
    eff_curve = None
    eff_source = ""
    if geometry_canon:
        from gamma.calibration.efficiency_autoload import find_efr_file
        loaded = load_efficiency_for_geometry(
            spec.geometry or geometry_hint or "",
            spec.detector_id or "",
        )
        # DEEP-01 (Project #5 wave 2 P1-1): the loader now distinguishes
        # "no file found" (None) from "file found, fit failed"
        # (EFFICIENCY_FIT_FAILED sentinel).  Surface the failure case
        # explicitly so the operator sees that the report is
        # efficiency-uncorrected by accident rather than by absence
        # of calibration data.
        if loaded is EFFICIENCY_FIT_FAILED:
            eff_curve = None
            eff_source = "fit_failed"   # operator-visible reporting field
            logger.warning(
                "staged_pipeline: efficiency fit failed for geometry=%r "
                "detector=%r — proceeding without efficiency correction; "
                "see efficiency_autoload warning above for the offending file.",
                spec.geometry or geometry_hint or "",
                spec.detector_id or "",
            )
        else:
            eff_curve = loaded
            eff_source = find_efr_file(spec.geometry or geometry_hint or "",
                                       spec.detector_id or "") or ""
            # T41 (BUG-40 (b) hardening) — efficiency-file detector content
            # fingerprint gate. The path-level cyrillic_to_latin_collision
            # predicate cannot catch the case where the .efr was found
            # under the right detector directory but its `[detector;…]`
            # header records a DIFFERENT physical instance (serial-year
            # mismatch). Real incident: detectors/Gamma-1S/efficiency/
            # …_SN-01/…Marinelli.efr loaded for a spectrum whose
            # CONFIGNAME is Гамма-1С №SN-02 → activity bias −96% to
            # −97% on Am-241/Ti-44. Surface the mismatch into the
            # existing detector_fallback record.
            if eff_source and detector_fallback_dict is not None:
                try:
                    from gamma.calibration.efficiency_provenance import (
                        check_efr_detector_match,
                    )
                    _expected_str = (
                        spec.extras.get("lsrm_config", "")
                        or spec.detector_id
                        or ""
                    )
                    _mismatch = check_efr_detector_match(
                        eff_source, _expected_str
                    )
                    if _mismatch:
                        # F-115: detector strings + basename carry the
                        # certified-source / instrument S/N — keep them
                        # in the operator log only, never in the
                        # serialised detector_fallback dict (which lands
                        # in the JSON report verbatim).
                        detector_fallback_dict[
                            "efficiency_detector_mismatch"
                        ] = {
                            "code": _mismatch["code"],
                            "expected_serial_year": list(
                                _mismatch.get("expected_serial_year") or []
                            ),
                            "actual_serial_year": list(
                                _mismatch.get("actual_serial_year") or []
                            ),
                        }
                        logger.warning(
                            "T41 efficiency-file detector serial mismatch: "
                            "spectrum=%r vs .efr=%r (file: %s) — efficiency "
                            "curve belongs to a different physical "
                            "instrument; activity values may be biased.",
                            _mismatch["expected_detector"],
                            _mismatch["actual_detector"],
                            _mismatch["efr_file_basename"],
                        )
                except Exception as _t41_exc:
                    # Validator must never block efficiency loading.
                    # Log at DEBUG so silent-handler ceiling (DEEP-06)
                    # is not regressed — this branch IS observed when
                    # diagnostics fail (corrupt .efr, etc).
                    logger.debug(
                        "T41 validator skipped on %s: %s",
                        eff_source, _t41_exc,
                    )

        # Re-run completeness with efficiency if curve was loaded
        if eff_curve is not None:
            try:
                completeness_result = compute_completeness(
                    detected_nuclides=final_detected,
                    unmatched_peaks_E_area=unmatched_for_dc,
                    efficiency_at=eff_curve.efficiency_at,
                )
            except Exception as exc:  # DEEP-06
                logger.warning(
                    "efficiency-aware completeness re-run failed (%r); "
                    "keeping efficiency-free result — completeness score "
                    "will not reflect ε(E) weighting.",
                    exc,
                )

    # Mode tagging (background vs sample)
    analysis_mode = "background_7line" if is_background else "sample_anchor_rank"

    # ─── Stage 2 if explicitly allowed or auto-escalation enabled ───
    run_stage2 = False
    if allow_stage2:
        run_stage2 = True
    elif auto_escalate and rec_stage == 2:
        run_stage2 = True
    elif user_confirmed_stage2_nuclides is not None:
        run_stage2 = True

    if run_stage2:
        # Combine ЕРН + technogenic: identify_nuclides will pick the best
        # match overall, then disambiguate resolves overlaps.
        extra_cand = user_confirmed_stage2_nuclides or candidates_for_stage(2)
        # Re-run identification with cumulative ЕРН+S2 candidate list so
        # disambiguate has full context.
        combined = list(stage1.candidates_considered)
        for c in extra_cand:
            if c not in combined:
                combined.append(c)
        stage2 = _run_stage(spec, peaks, 2, fwhm_at_ch, window,
                            extra_candidates=combined,
                            line_window_overrides_keV=(line_window_overrides_keV or None),
                            fwhm_keV_at_energy_fn=_fwhm_keV_fn_for_id)
        stage2.stage = 2
        stages.append(stage2)
        final_detected = list(stage2.detected)
        final_unmatched = stage2.unmatched_peaks

    # ─── Stage 3 (opt-in only) ───
    if allow_stage3 or user_confirmed_stage3_nuclides is not None:
        extra_cand = user_confirmed_stage3_nuclides or candidates_for_stage(3)
        combined = [n.nuclide for n in final_detected] + extra_cand
        # Pull ALL candidates (stage 1+2+3) for disambiguate-aware run.
        all_cand = list(stage1.candidates_considered)
        for c in extra_cand:
            if c not in all_cand:
                all_cand.append(c)
        if run_stage2:
            for c in candidates_for_stage(2):
                if c not in all_cand:
                    all_cand.append(c)
        # F-456 / BUG-Y (2026-06-23): early singleton pre-calibration for Stage 3.
        # Stage 3 runs with stored SPE calibration, which for short NaI fixtures
        # (AmTiCsEu) has 10+ keV drift at high E (Ti-44 1157, Eu-152 1408).
        # This causes: (a) false positives Tc-99m/In-111/Co-57 matching mis-shifted
        # peaks, (b) Ti-44 missed because 10.9 keV shift exceeds matching window.
        # Fix: apply F-453-style singleton recalibration before Stage 3 runs.
        # IMPORTANT: the pipeline's anchor_matches (from Step 5α) uses max_rank=10,
        # which excludes calibration-tier ranks 15-17 (Sc-44/Ti-44/Eu-152). These
        # high-E anchors are exactly what's needed to correct the 10+ keV high-E drift.
        # Solution: build a separate full-rank anchor set (max_rank=99) so the
        # fixture-fingerprint gate (Cs-137 + Am-241 both visible) can include Ti-44
        # and Eu-152. The full F-453 recal still runs later and refines.
        _f456_am_full = find_anchor_matches(
            peaks, spec,
            fwhm_provider_keV=_fwhm_provider_keV,
            max_rank=99,
            window_fwhm_multiple=fwhm_window_multiple,
        )
        _f456_extras = _f453_build_singleton_extras(
            _f456_am_full, _fwhm_provider_keV, []
        )
        # F-456 guard: only pre-calibrate when fixture-fingerprint fired
        # (calibration-tier rank ≥ CALIBRATION_RANK_START present, i.e.
        # AmTiCsEu with Am-241 + Cs-137 both visible). For natural spectra
        # (Ra-226, Th-232) the fingerprint gate in find_anchor_matches
        # suppresses ranks 15-17 → _f456_am_full has only standard naturals
        # → skip F-456 to avoid corrupting energy_cal before Phase D.
        _f456_calib_tier_present = any(
            getattr(getattr(_am, "anchor", None), "rank", 0)
            >= CALIBRATION_RANK_START
            for _am in (_f456_am_full or ())
        )
        if _f456_extras and _f456_calib_tier_present:
            try:
                _f456_cal, _f456_diag = recalibrate_from_multiplet_centroids(
                    spec, [],
                    fwhm_provider_keV=_fwhm_provider_keV,
                    extra_anchors=_f456_extras,
                    use_cluster_global=False,
                )
                if _f456_cal is not None and _f456_diag.phase_C_applied:
                    spec.energy_cal = _f456_cal
                    if hasattr(spec, "energy_cal_degree"):
                        spec.energy_cal_degree = _f456_diag.degree_used or 2
                    if hasattr(spec, "energy_cal_source"):
                        spec.energy_cal_source = "F-456_precal_for_stage3"
            except Exception as _exc:
                logger.debug("F-456 pre-cal for Stage 3 skipped: %s", _exc)
        stage3 = _run_stage(spec, peaks, 3, fwhm_at_ch, window,
                            extra_candidates=all_cand,
                            line_window_overrides_keV=(line_window_overrides_keV or None),
                            fwhm_keV_at_energy_fn=_fwhm_keV_fn_for_id)
        stage3.stage = 3
        stages.append(stage3)
        final_detected = list(stage3.detected)
        final_unmatched = stage3.unmatched_peaks

    # ─── F-89a background-subtraction status — always surfaced ───
    # Per user feedback 2026-05-29: the report must NEVER silently
    # omit whether the spectrum was background-subtracted. Three
    # (now four, F-131) possible states:
    if bg_sub_result is not None:
        # F-131: если фон был авто-подобран И применён, помечаем
        # отдельно, чтобы отчёт отличал ручной --background-path
        # от автоматического подбора.
        if auto_bg_applied_note is not None:
            background_status_value = "auto_resolved_from_directory"
        else:
            background_status_value = "subtracted_from_external_file"
    elif getattr(spec, "background_embedded", None) is not None:
        # An embedded background was present but we did NOT subtract
        # it automatically — note that explicitly so the user knows.
        background_status_value = "embedded_present_not_subtracted"
    else:
        background_status_value = "absent_no_subtraction"

    # ─── F-89d post-filter: drop nuclides from chains suppressed by
    # the filename-binding rule. On a `Th232_*.spe` source the user
    # has bound the spectrum to the Th-232 chain only; any U-chain
    # identification (Bi-214, Pb-214, Pb-210, Ra-226) is, by user
    # methodology ("Радия точно нет"), a false positive from NaI
    # 63×63 Compton overlap (609 keV peak ≡ Tl-208 583 shifted by
    # fit imprecision).
    chain_filtered_out_names: List[str] = []
    if chain_dominance_out.suppressed_chains:
        suppressed_proxy_sets: List[set] = []
        for ch in chain_dominance_out.suppressed_chains:
            if ch == "U-238":
                suppressed_proxy_sets.append(set(U238_PROXY_NUCLIDES))
            elif ch == "Th-232":
                suppressed_proxy_sets.append(set(TH232_PROXY_NUCLIDES))
        if suppressed_proxy_sets:
            forbidden = set.union(*suppressed_proxy_sets)
            # Don't drop nuclides that were explicitly hinted by
            # filename — those are user-asserted.
            forbidden -= set(filename_isotope_hints or ())
            kept = []
            for ni in final_detected:
                if ni.nuclide in forbidden:
                    chain_filtered_out_names.append(ni.nuclide)
                else:
                    kept.append(ni)
            final_detected = kept

    # ──────────────────────────────────────────────────────────────
    # F-116 / v1.17.5 — Out-of-chain Stage-3 suppression.
    #
    # Когда цепочка Th-232 жёстко доминирует И filename-подсказки
    # содержат ТОЛЬКО члены Th-цепи (Th-232, Tl-208, Pb-212, Ac-228,
    # Bi-212, Th-228), мы отбрасываем все нуклиды, найденные на
    # Stage 3, которые: (а) не принадлежат Th-цепи и (б) не присутствуют
    # в filename-подсказках, ЕСЛИ они не подтверждены жёстко:
    #   σ ≥ 10 AND ≥ 2 matched lines AND CI ≥ 5.
    # Цель — устранить ложноположительный Co-57 (122 keV) и подобные
    # артефакты из «попутных» Stage-3 кандидатов, видных на Th-232
    # фоне (D-22 / Co-57 false-positive в SESSION_DIRECTIVES).
    #
    # Зеркальное правило применяется к доминантной U-238 цепи.
    # ──────────────────────────────────────────────────────────────
    # F-116-v2 (#124, 2026-06-03): union обеих цепей если обе DOMINANT.
    # Раньше тут стоял `if/else` — при одновременной dominance Th-232 и
    # U-238 (типичный natural-background сценарий) выбиралась только
    # Th-цепь, и члены U-цепи (Bi-214, Pb-214, Ra-226, Pb-210) отбрасывались
    # как «вне цепочки Th-232» — spurious suppression реальных background-
    # компонентов. Реализация вынесена в _apply_f116_out_of_chain_suppression()
    # для unit-test изоляции (tests/test_f116_dual_chain.py).
    final_detected, out_of_chain_suppressed = _apply_f116_out_of_chain_suppression(
        final_detected, chain_dominance_out, filename_isotope_hints,
    )
    # REVIEW NOTE 2026-06-03 (#124): pre-v1.19.1 behaviour was `if/elif` —
    # chose Th-only chain_set when both chains DOMINANT, suppressing U-cycle
    # members (Bi-214, Pb-214 etc.) despite U-238 being independently flagged
    # DOMINANT. Fixed by union-of-flagged-chains. Backward compat preserved
    # for single-chain cases. См. PATCHED.json line 996 review_note.

    # ──────────────────────────────────────────────────────────────
    # F-136 / v1.17.7 — Filename single-isotope binding suppression.
    #
    # Контракт навсегда: если имя файла связывает источник с одним
    # моно-нуклидом (Cs-137 / K-40 / Co-60 / Am-241 / Na-22 ...) И
    # ни одна цепочка (Th-232 / U-238) не доминирует, то ВСЕ остальные
    # detected нуклиды должны иметь ВЕСЬ комплект подтверждений:
    #   (1) ≥ 2 matched_lines с peak_area > 0 (реальные интегрированные
    #       пики, не library window-match без peak'а)
    #   (2) ХОТЯ БЫ ОДНА линия совпадает с библиотечной в пределах
    #       0.5·FWHM(E) — это «характеристическая» линия. Без этого
    #       library window-match при FWHM=25 кэВ может тащить любой
    #       Compton plateau как «подтверждение».
    #
    # Иначе Stage-2/3 нуклиды типа U-235 (185 keV → Compton residual)
    # с большими ΔE-сдвигами систематически проникают как false positives.
    #
    # F-116 уже покрывает случай С доминантной цепочкой; F-136
    # симметрично покрывает источники БЕЗ цепочечной dominance.
    # ──────────────────────────────────────────────────────────────
    if (filename_isotope_hints
            and chain_dominance_out is not None
            and not (chain_dominance_out.th232 or chain_dominance_out.u238)):
        hint_set = set(filename_isotope_hints)
        kept_f136: List[NuclideIdentification] = []
        for ni in final_detected:
            nuc = ni.nuclide
            if nuc in hint_set:
                kept_f136.append(ni)
                continue
            n_real_peaks = 0
            best_dE_ratio = float("inf")    # min |ΔE|/FWHM по матчам
            for m in ni.matched_lines:
                pa = getattr(m, "peak_area", None)
                if pa is None or float(pa or 0.0) <= 0.0:
                    continue
                n_real_peaks += 1
                # |ΔE|/FWHM
                ch = getattr(m, "peak_channel", None)
                E_lib = getattr(m, "library_E_keV", None)
                if ch is None or E_lib is None:
                    continue
                try:
                    E_meas = float(spec.channel_to_energy(int(ch)))
                except Exception:
                    continue
                fwhm = _fwhm_provider_keV(float(E_lib))
                if fwhm <= 0:
                    continue
                ratio = abs(E_meas - float(E_lib)) / fwhm
                if ratio < best_dE_ratio:
                    best_dE_ratio = ratio
            # F-136 + F-143 / v1.17.7: 0.25·FWHM — ужесточённый критерий
            # (LSRM Алгоритмические основы §9). На NaI 63×63 это ~6-7 кэВ
            # при E=200 keV, отсекает library-window matches от Compton-
            # residual (Cs-137 ROI 100-300 keV где регистрируются ложные
            # U-235 185.7 keV / Ra-226 186.2 keV / Bi-211 / ... matches).
            confirmed_by_dE = best_dE_ratio <= 0.25
            if n_real_peaks >= 2 and confirmed_by_dE:
                kept_f136.append(ni)
            else:
                out_of_chain_suppressed.append(
                    f"F-136: подавлен {nuc} (источник связан с "
                    f"{', '.join(sorted(hint_set))}; "
                    f"real_peaks={n_real_peaks}/2 нужно, "
                    f"|ΔE|/FWHM_min={best_dE_ratio:.2f}/0.25 нужно)"
                )
        final_detected = kept_f136

    # ──────────────────────────────────────────────────────────────
    # F-335.4 / v1.18.18.8 — Chain-proxy single-line guard.
    #
    # Контракт навсегда: dochterный нуклид Th-232 / U-238 цепочки
    # (Pb-212, Tl-208, Ac-228, Bi-212, Th-228, Pb-214, Bi-214, Pb-210)
    # принимается ТОЛЬКО при выполнении ХОТЯ БЫ ОДНОГО условия:
    #   (a) Соответствующая цепочка (Th-232 / U-238) объявлена
    #       доминантной по chain_dominance (≥2 anchor линий или
    #       trump-card 2614);
    #   (b) Нуклид имеет ≥2 собственных matched_lines в этом спектре;
    #   (c) Имя файла explicitly содержит этот нуклид как hint
    #       (например `Pb212_calib.spe`).
    #
    # Иначе — drop (false ID).  Примеры:
    #   • K-40 sample: Phase D matched Pb-212 238.6 (одиночный peak
    #     рядом с K-40 backscatter 217). Нет других Th-линий. → drop.
    #   • Cs-137 sample: residual в 200 кэВ matched как U-235/Pb-214.
    #     Одиночная линия, цепочка не dominant. → drop.
    #   • Th-232 demo: Tl-208 583+2614 + Ac-228 911+338 → th232_dominant
    #     = True → Pb-212 проходит (пункт (a)).
    #
    # F-136 покрывает только filename-binding случай; F-335.4
    # симметрично покрывает безымённые/общие спектры.
    # ──────────────────────────────────────────────────────────────
    f335_4_suppressed: List[str] = []
    th_dom = bool(chain_dominance_out and chain_dominance_out.th232)
    u_dom = bool(chain_dominance_out and chain_dominance_out.u238)
    hint_set_335 = set(filename_isotope_hints or ())
    kept_f335_4: List[NuclideIdentification] = []
    for ni in final_detected:
        nuc = ni.nuclide
        is_th_proxy = nuc in _TH_CHAIN_NUCLIDES and nuc != "Th-232"
        is_u_proxy = nuc in _U_CHAIN_NUCLIDES and nuc != "U-238"
        if not (is_th_proxy or is_u_proxy):
            kept_f335_4.append(ni)
            continue
        # Pass (a): chain dominant
        if (is_th_proxy and th_dom) or (is_u_proxy and u_dom):
            kept_f335_4.append(ni)
            continue
        # Pass (b): ≥2 matched_lines
        if len(ni.matched_lines) >= 2:
            kept_f335_4.append(ni)
            continue
        # Pass (c): filename hint explicit for this nuclide
        if nuc in hint_set_335:
            kept_f335_4.append(ni)
            continue
        # Drop with diagnostic
        chain_name = "Th-232" if is_th_proxy else "U-238"
        f335_4_suppressed.append(
            f"F-335.4: подавлен {nuc} (одиночная линия, цепочка "
            f"{chain_name} не доминирует, имя файла не содержит {nuc})"
        )
    final_detected = kept_f335_4
    # Merge into existing out_of_chain_suppressed log so диагностики
    # видны в pipeline_notes.
    if f335_4_suppressed:
        out_of_chain_suppressed.extend(f335_4_suppressed)

    notes = []
    notes.append(f"FWHM model source: {fwhm_src}")
    # F-243 / v1.18.29 — surface bg-control pre-check failures as a non-
    # blocking warning. Computed at the Step 2 bg-subtraction branch.
    if f243_bg_warning:
        notes.append(f243_bg_warning)
    # F-131 / v1.17.7 — нарративные заметки про авто-поиск фона.
    if auto_bg_applied_note:
        notes.append(auto_bg_applied_note)
    elif auto_bg_suggest_note:
        notes.append(auto_bg_suggest_note)
    # F-140 / v1.17.7 — нарративная заметка про авто-извлечение массы.
    if auto_mass_source and sample_mass_kg is not None:
        notes.append(
            f"F-140: масса образца авто-определена = {sample_mass_kg*1000:.1f} г "
            f"(источник: {auto_mass_source})"
        )
    # F-378 / v1.18.25 — warn о расхождении CLI mass vs .spe mass.
    if mass_mismatch_note:
        notes.append(mass_mismatch_note)

    # F-142 / v1.17.7 — Cs-Kα 32 кэВ self-calibration check.
    # При binding=Cs-137: ожидаем Ba Kα ~32 кэВ как характеристический
    # secondary peak от 661.66 IC. Если найденный пик в окне 20-50 кэВ
    # отличается от 32 на > 0.5·FWHM(32), энергетическая калибровка
    # вероятно смещена. Эмитируем warning для пользователя.
    if (filename_isotope_hints and "Cs-137" in filename_isotope_hints):
        BA_KA_E = 32.06
        fwhm_at_32 = _fwhm_provider_keV(BA_KA_E)
        nearest_ba = None
        nearest_dE = float("inf")
        for pk in peaks:
            try:
                E_pk = float(spec.channel_to_energy(int(pk.channel)))
            except Exception:
                continue
            if 20.0 <= E_pk <= 50.0:
                dE = abs(E_pk - BA_KA_E)
                if dE < nearest_dE:
                    nearest_dE = dE
                    nearest_ba = E_pk
        if nearest_ba is not None and fwhm_at_32 > 0:
            ratio = nearest_dE / fwhm_at_32
            if ratio > 0.5:
                notes.append(
                    f"F-142: WARNING — возможен сдвиг энергетической "
                    f"калибровки на низких E. Ba Kα ожидается на "
                    f"{BA_KA_E:.1f} кэВ (характеристический secondary от "
                    f"Cs-137 661.66 IC), но ближайший пик в окне 20-50 кэВ "
                    f"найден на {nearest_ba:.2f} кэВ "
                    f"(|ΔE|={nearest_dE:.2f}, ΔE/FWHM={ratio:.2f}). "
                    f"Рекомендуется проверить калибровку E(N) на низких "
                    f"энергиях или пересчитать через --recalibrate-on-"
                    f"anchor-disagreement."
                )
            else:
                notes.append(
                    f"F-142: Ba Kα калибровка OK — пик найден на "
                    f"{nearest_ba:.2f} кэВ (ожидается {BA_KA_E:.1f}, "
                    f"|ΔE|={nearest_dE:.2f}, ΔE/FWHM={ratio:.2f})."
                )
    # F-129 / v1.17.7 — нарративная заметка о методе поиска пиков, если
    # отличается от Mariscotti (default).
    if peak_search_method != "mariscotti":
        notes.append(
            f"F-129: метод поиска пиков — '{peak_search_method}' "
            f"(альтернатива Mariscotti)."
        )
        if peak_method_compare is not None:
            notes.append(
                f"F-129: сравнение Mariscotti vs convolution: "
                f"совпадение {peak_method_compare.get('agreement_fraction', 0)*100:.0f}%, "
                f"найдено пиков {peak_method_compare.get('n_mariscotti', 0)} / "
                f"{peak_method_compare.get('n_convolution', 0)} "
                f"(средняя |Δch|={peak_method_compare.get('mean_residual_channels', 0):.2f})."
            )
    # F-123 — раскрытие применённых per-line оверрайдов окна
    if line_window_overrides_keV:
        items = ", ".join(
            f"{nuc}@{E:.1f}=±{w:.1f}кэВ"
            for (nuc, E), w in sorted(line_window_overrides_keV.items())
        )
        notes.append(
            f"F-123: расширенные окна идентификации (доминантная цепочка): {items}"
        )
    # F-116 — surface suppression notes (RU narrative)
    for note_line in out_of_chain_suppressed:
        notes.append(note_line)
    if not run_stage2 and rec_stage == 2:
        notes.append("RECOMMENDATION: " + rec_reason)
    elif not run_stage2:
        notes.append("Stage 1 ЕРН достаточен — Stage 2/3 не требуются.")
    # F-89b/d notes
    if filename_isotope_hints:
        notes.append(
            f"Step 7A.1: filename isotope hints driving candidate list: "
            f"{', '.join(filename_isotope_hints)}"
        )
    if chain_filtered_out_names:
        notes.append(
            f"F-89d chain suppression: dropped {', '.join(chain_filtered_out_names)} "
            f"from identifications. "
            + chain_dominance_out.suppression_reason
        )

    # F-88 chain dominance + K-40 overlap warnings (user methodology)
    if chain_dominance_out.th232:
        notes.append(
            f"Step 5α: Th-232 chain DOMINANT — "
            f"{chain_dominance_out.reason}. "
            f"Tl-208 / Pb-212 / Ac-228 / Bi-212 confirmed as "
            f"strong-prior candidates."
        )
    if chain_dominance_out.u238:
        notes.append(
            f"Step 5α: U-238 chain DOMINANT — "
            f"{chain_dominance_out.reason}. "
            f"Bi-214 / Pb-214 / Pb-210 confirmed as strong-prior candidates."
        )
    if k40_overlap_warning_out:
        notes.append(
            "WARNING (F-88): K-40 1460.82 keV peak overlaps Ac-228 "
            "1459.20 keV (I=0.85%) — on NaI 63×63 the doublet is "
            "unresolvable. K-40 area is contaminated; confirmation "
            "requires either deconvolution against confirmed Tl-208 "
            "anchor or a separate Ac-228 reference measurement."
        )

    # ══════════════════════════════════════════════════════════════
    # F-84 / Round 5 (v1.13.0): post-identification wiring of
    # multiplet deconvolution → activities (Bq, Bq/kg) → MDA
    # ══════════════════════════════════════════════════════════════
    deconvolutions_out: Optional[List[DeconvolutionResult]] = None
    activities_out: Optional[List[ActivityResult]] = None
    specific_activities_out: Optional[dict] = None
    mda_per_line_out: Optional[dict] = None

    # F-145 / v1.17.8 — multiplet self-calibration diagnostic; populated
    # внутри блока forced_clusters если запускается. None — F-145 не пытался.
    f145_diag: Optional[SelfCalibrationDiag] = None

    # ─── F-84b: targeted multiplet deconvolution ───
    # Default OFF. When ON, this rebinds peak_area / peak_area_source on
    # any LineMatch whose channel lies in a multiplet cluster — most
    # notably the 600–680 keV Bi-214 609 / Cs-137 661 / Cs-134 604+795
    # cluster on Th-rich samples. Isolated lines are untouched.
    #
    # F-118 / v1.17.5 — Chain-aware multiplet auto-discovery:
    # когда Th-232 цепочка ДОМИНАНТНА, всегда эмиттируются жёстко
    # закреплённые ROI M1 (754-1114 keV) и M2 (1430-1786 keV) с
    # коннектами по библиотечным интенсивностям (см. coupled_multiplet).
    # Эти результаты идут ПЕРВЫМИ в списке кластеров, легаси
    # free-parameter подгонка — после, по возрастанию E (D-03/F-108).
    id_result_for_activities = None
    if apply_deconvolution and final_detected:
        # Reconstruct a minimal IdentificationResult container so that
        # `apply_multiplet_deconvolution` can iterate detected_nuclides
        # consistently regardless of which stage produced them.
        from gamma.identification.identify import IdentificationResult as _IR
        synthetic_id = _IR(
            detector_type=detector_type,
            window=window,
            candidates_considered=sum(
                len(s.candidates_considered) for s in stages
            ),
            detected_nuclides=tuple(final_detected),
            rejected_nuclides=(),
            unmatched_peaks=tuple(final_unmatched),
            notes="",
        )
        try:
            # F-118 / F-121: chain-forced clusters first
            #   • Th-232 dominant → M1 + M2 (F-118)
            #   • Ra-226 / U-238 dominant → U1 + U2 + U3 (F-121, v1.17.6)
            # F-120 (v1.17.6): peak-image (Гаусс + хвост) для NaI.
            # F-145 (v1.17.8): two-phase self-calibration —
            #   Phase A — forced_clusters_pA с free_centroids=True (выяв. drift)
            #   Phase B — фильтр Phase A результатов в gamma.calibration.
            #               multiplet_self_calibration.recalibrate_*
            #   Phase C — polyfit E(N) на anchor'ах фитированных центроидов
            #   фаза Д — повторный forced_clusters на новой калибровке,
            #               free_centroids=False (locked-passport)
            # Если Phase C не применил калибровку → forced_clusters = Phase A
            # (free-centroid side-fit отбрасывается; используются locked
            # площади из основного fit'a Phase A).
            forced_clusters = run_chain_forced_multiplets(
                spec, fwhm_at_ch, _fwhm_provider_keV,
                chain_dominance_out, filename_isotope_hints,
                use_peak_image=True, detector_type=detector_type,
                # F-127 / v1.17.7: per-line T(E) калибровка для NaI 63×63
                use_T_E_model=True,
                # F-126 / v1.17.7: нелинейный refinement (σ_scale + dE)
                # после линейного NNLS старта. Принимается только если
                # χ²/ν улучшается ≥5 % и амплитуды остаются физически
                # валидны (см. coupled_multiplet.py).
                nonlinear_refine=True,
                # F-145 / v1.17.8: Phase A free-centroid side-fit
                free_centroids=True,
            )

            # F-145 Phase B + C — попытаться refit'нуть E(N) по центроидам.
            # F-445 / v1.30.3: pass counts + continuum arrays for cluster-Δ.
            f145_new_cal = None
            f445_counts_arr = getattr(spec, "counts", None)
            f445_cont_arrays = _f445_build_continuum_arrays(forced_clusters)
            # F-453 (BUG-38 follow-up, 2026-06-23): singleton-anchor
            # fallback для short NaI fixtures без chain forced_clusters
            # (AmTiCsEu, Cs-Co). Если мультиплетов мало, F-145 машина
            # подберёт singleton anchors из anchor_matches через extra_anchors.
            f453_singleton_extras = _f453_build_singleton_extras(
                anchor_matches, _fwhm_provider_keV, forced_clusters,
            )
            try:
                f145_new_cal, f145_diag = recalibrate_from_multiplet_centroids(
                    spec, forced_clusters,
                    fwhm_provider_keV=_fwhm_provider_keV,
                    counts_arr=f445_counts_arr,
                    continuum_arrays=f445_cont_arrays or None,
                    use_cluster_global=True,
                    extra_anchors=f453_singleton_extras,
                )
            except Exception as exc_f145:
                f145_diag = SelfCalibrationDiag(
                    attempted=True, reason=f"F-145 refit exception: {exc_f145!r}",
                )

            # F-145 фаза Д — если калибровка применена, пересчитать
            # forced clusters на новой шкале (locked-passport).
            if f145_new_cal is not None:
                from dataclasses import replace as _dc_replace_spec
                try:
                    spec_recal = _dc_replace_spec(
                        spec,
                        energy_cal=tuple(f145_new_cal),
                        energy_cal_degree=(
                            f145_diag.degree_used if f145_diag else None
                        ),
                        energy_cal_source="F-145_multiplet_self_calibration",
                    )
                    # F-145 фаза Д — СМЯГЧЁННЫЙ locked-passport.
                    # Центроиды свободны в малом окне ±tolerance·FWHM, что
                    # компенсирует невозможность точного покрытия всего
                    # спектра 2-3 anchor'ами polyfit'a. Дефолтный tolerance
                    # 0.15·FWHM — компромисс между ригидностью паспортных
                    # энергий и адаптацией к остаточному нелинейному
                    # дрейфу шкалы. См. PHASE_D_CENTROID_TOLERANCE_FRAC.
                    # F-453 carve-out (BUG-38 follow-up, 2026-06-23):
                    # короткие NaI fixtures (AmTiCsEu, Cs-Co) → пустой
                    # forced_clusters (chain_dominance не активен). F-145
                    # refit построен на singleton extras; Phase D refit
                    # тоже вернёт пустой список, и стандартное условие
                    # `chi2_D <= chi2_A AND forced_clusters_D` отвергнет
                    # cal. При пустом forced_clusters χ² regression
                    # risk = 0 by definition — принять cal безусловно.
                    f453_singleton_only = (
                        not forced_clusters
                        and f453_singleton_extras
                    )
                    if f453_singleton_only:
                        # Skip Phase D refit entirely (нечего пересчитывать).
                        forced_clusters_D = []
                        chi2_A_sum = 0.0
                        chi2_D_sum = 0.0
                    else:
                        forced_clusters_D = run_chain_forced_multiplets(
                            spec_recal, fwhm_at_ch, _fwhm_provider_keV,
                            chain_dominance_out, filename_isotope_hints,
                            use_peak_image=True, detector_type=detector_type,
                            use_T_E_model=True,
                            nonlinear_refine=True,
                            free_centroids=True,
                            centroid_window_frac=PHASE_D_CENTROID_TOLERANCE_FRAC,
                        )
                        # Проверка: χ² фаза Д действительно лучше Phase A locked
                        chi2_A_sum = sum(
                            getattr(d, "chi2_per_dof", 0.0) for d in forced_clusters
                        )
                        chi2_D_sum = sum(
                            getattr(d, "chi2_per_dof", 0.0) for d in forced_clusters_D
                        )
                    if f453_singleton_only:
                        # Принять cal без Phase D refit.
                        spec = spec_recal
                        if f145_diag is not None:
                            f145_diag.phase_C_applied = True
                            f145_diag.reason += (
                                " | F-453 carve-out: forced_clusters пуст,"
                                " singleton-only refit принят безусловно"
                            )
                    elif chi2_D_sum <= chi2_A_sum and forced_clusters_D:
                        # Принять фаза Д — обновить spec и forced_clusters
                        spec = spec_recal
                        forced_clusters = forced_clusters_D
                        if f145_diag is not None:
                            f145_diag.reason += (
                                f" | фаза Д принят: χ²_sum {chi2_A_sum:.2f} → "
                                f"{chi2_D_sum:.2f}"
                            )
                    else:
                        # Phase D rollback: keep old cal regardless of
                        # whether anchors came from per-component or
                        # cluster-Δ path. The carve-out experiment showed
                        # χ² regressions on Ra-226 (M1 61.69→66.56) when
                        # we accepted cal-only without Phase D refit, so
                        # we revert to the original "all-or-nothing" rule.
                        # F-445 still emits cluster-Δ anchors when per-
                        # component anchors fail filter — but downstream
                        # acceptance requires Phase D χ²_D ≤ χ²_A.
                        _f445_note_phase_d_rollback(
                            f145_diag, chi2_A_sum, chi2_D_sum,
                        )
                except Exception as exc_pD:
                    if f145_diag is not None:
                        f145_diag.phase_C_applied = False
                        f145_diag.reason += f" | фаза Д exception: {exc_pD!r}"

            # Legacy free-parameter multiplet auto-discovery
            # F-322 / v1.18.16 — pass through F-96 bg-anchors flag для
            # включения через analyze_lsrm_spe.
            # F-447 wire-in rolled back 2026-06-15: V1 (drop guest) и V2
            # (phantom-mark guest) обе разрушают Tl-208 на Th-232 demo —
            # см. 1_Version/v1.30.x/F-447_guest_only/NEGATIVE_RESULT.md.
            # Helpers _f447_* и kwarg f440_guest_only_filter оставлены в
            # deconvolve.py как dormant API для будущего исследования.
            new_id, free_deconv_list = apply_multiplet_deconvolution(
                synthetic_id, spec, fwhm_at_ch,
                overlap_threshold_fwhm=deconvolution_overlap_fwhm,
                continuum="step_linear",
                enable_f96_bg_anchors=enable_f96_bg_anchors,
            )

            # Defect 3 / F-118: пропустить пересечения легаси-кластеров
            # с уже учтёнными forced-кластерами.
            # F-381 / v1.18.25.2 — substantial-overlap criterion: auto
            # cluster дропается только если ≥50% его длины пересекается
            # с forced ROI (или forced полностью внутри auto). Раньше
            # любое касание сбрасывало auto — после F-374+F-381 expansion
            # auto-M3 на Ac-228 503-583 кэВ дотягивался до forced M1
            # 750+ кэВ кончиком и отбрасывался целиком.
            forced_ranges = [
                (d.roi_low_ch, d.roi_high_ch) for d in forced_clusters
            ]
            def _overlaps_forced(d) -> bool:
                d_len = max(1, d.roi_high_ch - d.roi_low_ch)
                for lo, hi in forced_ranges:
                    # Intersection length
                    inter_lo = max(d.roi_low_ch, lo)
                    inter_hi = min(d.roi_high_ch, hi)
                    if inter_hi <= inter_lo:
                        continue
                    inter_len = inter_hi - inter_lo
                    # Drop if ≥50% of auto inside forced, OR forced fully
                    # inside auto (substantial overlap either way).
                    if inter_len / d_len >= 0.5:
                        return True
                    f_len = max(1, hi - lo)
                    if inter_len / f_len >= 0.9:
                        return True
                return False
            free_filtered = [
                d for d in free_deconv_list if not _overlaps_forced(d)
            ]
            # Forced first, free after; both sorted by ROI low channel.
            forced_sorted = sorted(
                forced_clusters, key=lambda d: d.roi_low_ch
            )
            free_sorted = sorted(
                free_filtered, key=lambda d: d.roi_low_ch
            )
            combined_deconv = list(forced_sorted) + list(free_sorted)

            # Defect 3 / Дефект 3 — заменить peak_area на связанной площади
            # ВО ВСЕХ LineMatch, попадающих в любой ROI деконволюции
            # (источник = deconvolved_coupled для F-118, deconvolved
            # для легаси). Маппинг по ближайшему E_keV в окне 1.5·FWHM.
            from dataclasses import replace as _dc_replace
            from gamma.identification.identify import (
                IdentificationResult as _IR2,
                NuclideIdentification as _NI, LineMatch as _LM,
            )
            # Сначала собрать словарь (E_lib_keV, nuclide) → (area, sigma, source)
            comp_replacements: dict = {}
            for d in combined_deconv:
                if not getattr(d, "converged", False):
                    continue
                src_label = (
                    "deconvolved_coupled"
                    if str(getattr(d, "method", "")).startswith("coupled_")
                    else "deconvolved"
                )
                for comp, area, unc in zip(
                    d.components, d.areas, d.area_uncertainties,
                ):
                    key = (str(comp.nuclide), round(float(comp.line_E_keV), 2))
                    # Если несколько кластеров покрывают одну линию,
                    # forced (первый в списке) перекрывает легаси.
                    if key in comp_replacements:
                        continue
                    comp_replacements[key] = (
                        float(area), float(unc), src_label
                    )

            # Применить замену к каждому matched_lines
            # Defect 3: для F-118 forced-кластеров добавляем НОВЫЕ
            # LineMatch для компонент, которых не было в matched_lines
            # (например, Ac-228 911/964.77/969 могли не пройти Stage 1
            # пик-поиск — связанная подгонка их «открывает»).
            new_detected2 = []
            covered_by_nuclide: dict = {}  # nuclide → set of rounded E
            for ni in new_id.detected_nuclides:
                new_matches = []
                changed = False
                covered_by_nuclide.setdefault(ni.nuclide, set())
                for m in ni.matched_lines:
                    key = (str(m.nuclide), round(float(m.library_E_keV), 2))
                    repl = comp_replacements.get(key)
                    if repl is None:
                        # Также пробуем «приблизительный» матч в пределах 1.5·FWHM
                        # на случай минорных расхождений в E_lib.
                        E_lib = float(m.library_E_keV)
                        fwhm_E = float(_fwhm_provider_keV(E_lib))
                        best_key = None
                        best_d = 1e9
                        for k, v in comp_replacements.items():
                            if k[0] != str(m.nuclide):
                                continue
                            d_keV = abs(k[1] - E_lib)
                            if d_keV < 1.5 * fwhm_E and d_keV < best_d:
                                best_d = d_keV
                                best_key = k
                        if best_key is not None:
                            repl = comp_replacements.get(best_key)
                            if repl is not None:
                                covered_by_nuclide[ni.nuclide].add(best_key[1])
                    else:
                        covered_by_nuclide[ni.nuclide].add(key[1])
                    if repl is not None:
                        area, unc, src = repl
                        new_matches.append(_dc_replace(
                            m,
                            peak_area=area,
                            peak_area_uncertainty=unc,
                            peak_area_source=src,
                        ))
                        changed = True
                    else:
                        new_matches.append(m)

                # Дополняем новые LineMatch для компонент мультиплета,
                # которые ещё не присутствуют (только из forced/coupled
                # кластеров — чтобы не плодить ложные совпадения).
                for (rk_nuc, rk_E), (area, unc, src) in comp_replacements.items():
                    if rk_nuc != ni.nuclide:
                        continue
                    if rk_E in covered_by_nuclide[ni.nuclide]:
                        continue
                    if src != "deconvolved_coupled":
                        continue
                    if area <= 0:
                        continue
                    # Конвертация E_lib → канал, центр-канал по средней
                    # калибровке. Для активности достаточно area/source.
                    try:
                        ch_center = int(round(spec.energy_to_channel(rk_E)))
                    except Exception:
                        ch_center = 0
                    fwhm_keV_here = float(_fwhm_provider_keV(rk_E))
                    new_lm = _LM(
                        nuclide=ni.nuclide,
                        library_E_keV=float(rk_E),
                        library_I_pct=0.0,  # будет заменено ниже
                        peak_channel=ch_center,
                        peak_E_keV=float(rk_E),
                        peak_sigma=fwhm_keV_here / 2.355,
                        residual_keV=0.0,
                        is_characteristic=False,
                        peak_area=float(area),
                        peak_area_uncertainty=float(unc),
                        peak_area_source=src,
                        # BUG-34 Phase 1+2: explicit Gaussian sigma alias
                        gauss_sigma_keV=fwhm_keV_here / 2.355,
                    )
                    # Получить I_pct из библиотеки
                    rec = get_nuclide(ni.nuclide) or {}
                    lib_lines = rec.get("lines", [])
                    for ll in lib_lines:
                        if abs(float(ll[0]) - rk_E) < 0.5:
                            new_lm = _dc_replace(
                                new_lm, library_I_pct=float(ll[1]),
                            )
                            break
                    new_matches.append(new_lm)
                    covered_by_nuclide[ni.nuclide].add(rk_E)
                    changed = True

                if changed:
                    new_detected2.append(_dc_replace(
                        ni, matched_lines=tuple(new_matches),
                    ))
                else:
                    new_detected2.append(ni)

            id_result_for_activities = _IR2(
                detector_type=new_id.detector_type,
                window=new_id.window,
                candidates_considered=new_id.candidates_considered,
                detected_nuclides=tuple(new_detected2),
                rejected_nuclides=new_id.rejected_nuclides,
                unmatched_peaks=new_id.unmatched_peaks,
                notes=new_id.notes,
            )
            # F-445: refresh peak_E_keV on every LineMatch using the
            # post-F-145 spec.energy_cal so that primary_feps residuals
            # reflect the data-domain shift (not just the render overlay).
            # No-op when spec was not recalibrated.
            id_result_for_activities = _f445_refresh_match_energies(
                id_result_for_activities, spec, _IR2, _dc_replace,
            )
            final_detected = list(id_result_for_activities.detected_nuclides)
            deconvolutions_out = list(combined_deconv)
            n_coupled_replaced = sum(
                1 for ni in id_result_for_activities.detected_nuclides
                for m in ni.matched_lines
                if m.peak_area_source == "deconvolved_coupled"
            )
            n_legacy_replaced = sum(
                1 for ni in id_result_for_activities.detected_nuclides
                for m in ni.matched_lines
                if m.peak_area_source == "deconvolved"
            )
            notes.append(
                f"F-118 связанная подгонка цепочки Th-232: "
                f"{len(forced_clusters)} жёстко-закреплённых кластер(а); "
                f"всего деконволюций {len(deconvolutions_out)}, "
                f"заменено площадей: coupled={n_coupled_replaced}, "
                f"legacy={n_legacy_replaced}."
            )
        except Exception as exc:
            deconvolutions_out = []
            notes.append(
                f"F-118: связанная подгонка пропущена из-за ошибки: {exc!r}"
            )

    # ─── F-84c: per-nuclide activity in Bq, optional Bq/kg ───
    bg_available = bg_sub_result is not None
    from_bg_subtracted = bg_available

    if compute_activities and eff_curve is not None and final_detected:
        if id_result_for_activities is None:
            from gamma.identification.identify import IdentificationResult as _IR
            id_result_for_activities = _IR(
                detector_type=detector_type,
                window=window,
                candidates_considered=sum(
                    len(s.candidates_considered) for s in stages
                ),
                detected_nuclides=tuple(final_detected),
                rejected_nuclides=(),
                unmatched_peaks=tuple(final_unmatched),
                notes="",
            )
        # Pre-compute TCS factor dicts for cascade nuclides. The
        # close-geometry P/T scaling is picked from geometry_canon when
        # registered; unknown geometries fall back to the 5 cm reference.
        pt_func = peak_to_total_NaI_for_geometry(spec.geometry or "")
        coincidence_corrections = {}
        for ni in final_detected:
            try:
                cc = compute_tcs_corrections(
                    ni.nuclide, eff_curve, p_t_func=pt_func,
                )
                if cc:
                    coincidence_corrections[ni.nuclide] = cc
            except Exception:
                # Non-cascade nuclide or insufficient data — skip TCS.
                pass

        # Convert reference_datetime / measurement_datetime to datetime
        # objects acceptable by compute_activity (decay correction).
        meas_dt = getattr(spec, "start_datetime", None)
        try:
            activities_out = list(compute_activities_for_all(
                id_result_for_activities,
                efficiency_curve=eff_curve,
                live_time_s=float(spec.live_time),
                from_bg_subtracted=from_bg_subtracted,
                bg_available=bg_available,
                force_gross=not bg_available,
                coincidence_corrections=coincidence_corrections,
                decay_correction=(reference_datetime is not None
                                  and meas_dt is not None),
                reference_datetime=reference_datetime,
                measurement_datetime=meas_dt,
                # F-122 / v1.17.6 — self-attenuation wiring
                geometry_canonical=geometry_canon,
                sample_density_g_cm3=sample_density_g_cm3,
                matrix_composition=matrix_composition,
                # F-307 / v1.18.7 — pass-through opt-in flags
                enable_tcs_correction=enable_tcs_correction,
                tcs_detector_id=tcs_detector_id,
                enable_cutshall_self_abs=enable_cutshall_self_abs,
                cutshall_path_cm=cutshall_path_cm,
                cutshall_calib_density_g_cm3=cutshall_calib_density_g_cm3,
                enable_matrix_method=enable_matrix_method,
                matrix_method_energy_tolerance_keV=matrix_method_energy_tolerance_keV,
            ))
            # F-122: добавить нарративную заметку про коррекцию
            _geom_low = (geometry_canon or "").lower()
            if (sample_density_g_cm3 is not None
                    and geometry_canon
                    and ("marinelli" in _geom_low
                         or "маринелли" in _geom_low)):
                src_label = (
                    f" (источник: {auto_density_source})"
                    if auto_density_source else " (источник: CLI флаг)"
                )
                notes.append(
                    f"F-122: применена коррекция самопоглощения для "
                    f"геометрии {geometry_canon} (ρ_образец="
                    f"{sample_density_g_cm3:.3f} г/см³){src_label}. "
                    f"Per-line факторы F_ref/F_sample учтены в активности."
                )
        except Exception as exc:
            activities_out = []
            notes.append(
                f"compute_activities_for_all failed: {exc!r}"
            )

        # Specific activity in Bq/kg when sample mass is provided.
        if (sample_mass_kg is not None
                and sample_mass_kg > 0
                and activities_out):
            specific_activities_out = {}
            for ar in activities_out:
                if ar.is_valid():
                    spec_act = ar.A_Bq / sample_mass_kg
                    spec_unc = ar.sigma_A_Bq / sample_mass_kg
                    specific_activities_out[ar.nuclide] = (
                        float(spec_act), float(spec_unc),
                    )

    elif compute_activities and eff_curve is None:
        notes.append(
            "compute_activities requested but efficiency curve is "
            "unavailable for this geometry — activities skipped."
        )

    # ─── F-84d: per-line ISO 11929 MDA ───
    if compute_mda and eff_curve is not None:
        mda_per_line_out = {}
        # Standard MDA suite — always evaluated (regardless of detection)
        # so the report can state "Cs-137 not present; MDA = X Bq".
        default_suite: List[Tuple[str, float]] = [
            ("Cs-137",  661.66),
            ("Co-60",  1173.23),
            ("Co-60",  1332.49),
            ("K-40",   1460.82),
            ("Bi-214",  609.31),
            ("Bi-214", 1764.49),
            ("Tl-208", 2614.51),
            ("Ac-228",  911.20),
            ("Cs-134",  604.72),
            ("Cs-134",  795.86),
        ]
        if mda_suite_extra_lines_keV:
            for nucl, E in mda_suite_extra_lines_keV:
                default_suite.append((str(nucl), float(E)))

        # Augment the suite with every detected line so identification
        # carries its own MDA.
        for ni in final_detected:
            for m in ni.matched_lines:
                default_suite.append((ni.nuclide, float(m.library_E_keV)))

        # De-duplicate at the (nuclide, rounded_E) granularity.
        seen = set()
        for entry in default_suite:
            key_round = (entry[0], round(entry[1], 2))
            if key_round in seen:
                continue
            seen.add(key_round)
            nuclide, E = entry[0], float(entry[1])
            # Look up library intensity. If the nuclide line is unknown
            # we record a placeholder MDA (inf) so the caller can see
            # the gap rather than silently dropping the row.
            rec = get_nuclide(nuclide) or {}
            lib_lines = rec.get("lines", [])
            I_pct = 0.0
            for ll in lib_lines:
                if abs(float(ll[0]) - E) < 0.5:
                    I_pct = float(ll[1])
                    break
            if I_pct <= 0.0:
                continue   # Unknown line — skip (caller will not see it).
            # Estimate ROI background counts as 2·FWHM(E) window worth
            # of mean baseline. We use the smoothed counts around the
            # expected channel.
            try:
                ch_center = int(round(spec.energy_to_channel(E)))
            except Exception:
                continue
            fwhm_ch = max(1.0, float(fwhm_at_ch(ch_center)))
            roi_half = int(round(1.5 * fwhm_ch))
            lo = max(0, ch_center - roi_half)
            hi = min(spec.n_channels, ch_center + roi_half + 1)
            counts_arr = spec.counts
            try:
                roi = counts_arr[lo:hi]
                # Build a continuum estimate from the wings: 2 wings of
                # `roi_half` channels each. Take the mean of the two wings
                # as the per-channel baseline.
                left_lo = max(0, lo - roi_half)
                right_hi = min(spec.n_channels, hi + roi_half)
                left = counts_arr[left_lo:lo]
                right = counts_arr[hi:right_hi]
                wing_concat = list(left) + list(right)
                if wing_concat:
                    bg_per_ch = float(np.mean(wing_concat))
                else:
                    bg_per_ch = float(np.mean(roi)) if len(roi) else 0.0
                bg_in_roi = max(0.0, bg_per_ch * float(len(roi)))
            except Exception:
                bg_in_roi = 0.0
            try:
                eps_E = float(eff_curve.efficiency_at(E))
            except Exception:
                eps_E = 0.0
            if eps_E <= 0.0:
                continue
            try:
                res = mda_for_peak(
                    line_energy_keV=E,
                    background_counts_in_ROI=bg_in_roi,
                    live_time_s=float(spec.live_time),
                    efficiency=eps_E,
                    intensity_pct=I_pct,
                )
                mda_per_line_out[(nuclide, round(E, 2))] = res
            except Exception:
                pass

    # ──────────────────────────────────────────────────────────────
    # F-397 / v1.18.27 — bg peak detection (secondary pipeline run on
    # фоновый спектр). Активируется когда ВКЛ complete_workflow,
    # есть валидный bg path (explicit или auto-applied), и текущий
    # вызов сам не является bg-проходом (`_skip_background_analysis=False`).
    # Результат — отдельный StagedAnalysisResult, который HTML/JSON
    # reporters используют для рендера peak block в режиме «Фон».
    # ──────────────────────────────────────────────────────────────
    background_staged_result_out: Optional[StagedAnalysisResult] = None
    if (complete_workflow
            and not _skip_background_analysis
            and background_path
            and bg_sub_result is not None):
        try:
            background_staged_result_out = analyze_lsrm_spe(
                background_path,
                detector_type=detector_type,
                sigma_threshold=sigma_threshold,
                fwhm_window_multiple=fwhm_window_multiple,
                # Полный workflow на фоне — peaks + multiplets +
                # secondary peaks (нужен для рендера peak block).
                complete_workflow=True,
                # Защита от рекурсии: фон не должен искать свой собственный фон.
                background_path=None,
                background_auto="off",
                _skip_background_analysis=True,
                # Сохраняем выбор метода и фильтра пользователя.
                peak_search_method=peak_search_method,
                filter_narrow_peaks=filter_narrow_peaks,
                narrow_peak_fwhm_ratio=narrow_peak_fwhm_ratio,
                use_lsrm_id_window=use_lsrm_id_window,
                id_window_k_override=id_window_k_override,
                detector_class=detector_class,
                # Активность/MDA на фоне обычно бессмысленны (нет mass-context)
                # — отключаем явно, чтобы не падать на отсутствии efficiency.
                compute_activities=False,
                compute_mda=False,
                # Без mass/density — это фон, не образец.
                sample_mass_kg=None,
                sample_density_g_cm3=None,
            )
        except Exception:
            # F-397 — никогда не падает: bg-detection это диагностика,
            # её отсутствие не должно валить основной анализ.
            background_staged_result_out = None

    return StagedAnalysisResult(
        spec=spec,
        peaks=peaks,
        detector_type=detector_type,
        fwhm_at_661=fwhm_661,
        fwhm_model=fwhm_model,
        fwhm_model_source=fwhm_src,
        stages=stages,
        final_detected=final_detected,
        final_unmatched=final_unmatched,
        next_stage_recommended=rec_stage if not run_stage2 else None,
        next_stage_reason=rec_reason,
        sample_type_hint=sample_type,
        geometry_hint=geometry_hint,
        detector_hint=detector_hint,
        sample_type_canonical=sample_type_canon,
        geometry_canonical=geometry_canon,
        detector_canonical=detector_canon,
        is_background=is_background,
        residual_classifications=residual_cls,
        anchor_matches=anchor_matches,
        pattern_confirmations=pattern_confirmations,
        seven_line_check=seven_line,
        ci_gating=ci_gating,
        completeness=completeness_result,
        analysis_mode=analysis_mode,
        background_subtraction=bg_sub_result,
        efficiency_curve=eff_curve,
        efficiency_source=eff_source,
        detector_fallback=detector_fallback_dict,
        activities=activities_out,
        specific_activities_Bq_per_kg=specific_activities_out,
        sample_mass_kg=sample_mass_kg,
        mda_per_line=mda_per_line_out,
        deconvolution_results=deconvolutions_out,
        recalibration_diag=recalibration_diag,
        priority_findings=priority_findings_out,
        chain_dominance=chain_dominance_out,
        k40_ac228_overlap_warning=k40_overlap_warning_out,
        filename_isotope_hints=filename_isotope_hints,
        filename_chains_claimed=filename_chains_claimed,
        chain_filtered_out=chain_filtered_out_names,
        background_status=background_status_value,
        peak_search_method=peak_search_method,
        peak_search_method_comparison=peak_method_compare,
        # F-131 / v1.17.7 — авто-поиск фона
        auto_background_mode=str(background_auto or "off"),
        auto_background_candidates=(
            [c.to_dict() for c in auto_bg_candidates_list]
            if auto_bg_candidates_list else None
        ),
        # F-325 / v1.18.18.1 — record bg path for any subtraction (auto OR
        # explicit). Without this, when caller passed `background_path=...`
        # explicitly, the JSON output had `background_filename: ""` even
        # though subtraction succeeded. Reports then showed
        # "фон вычтен" without naming the file, leaving the operator
        # uncertain about reproducibility.
        auto_background_applied_path=(
            background_path
            if (auto_bg_applied_note or bg_sub_result is not None)
            else None
        ),
        # F-145 / v1.17.8 — multiplet self-calibration diagnostic
        multiplet_self_calibration_diag=(
            _serialize_f145_diag(f145_diag)
            if f145_diag is not None else None
        ),
        # F-397 / v1.18.27 — фон-как-StagedAnalysisResult для HTML toggle
        background_staged_result=background_staged_result_out,
        # F-332 / v1.18.18.5 — chart toggle (preserve gross + scaled-bg)
        gross_counts=f332_gross_counts,
        bg_counts_on_sample_grid=f332_bg_on_grid,
        bg_live_time=f332_bg_live_time,
        bg_scale_factor=f332_bg_scale,
        # F-QC-01 / v1.19.1 — per-peak Poisson |z|-test (BUG-35 / RAG-022)
        bg_quality_check=_bg_z_check_result,
        notes=notes,
    )


__all__ = [
    "StagedAnalysisResult", "StageResult",
    "analyze_lsrm_spe",
    "build_fwhm_model", "fwhm_keV_at_energy", "fwhm_model_legacy_abc",
    "FwhmModel",
]
