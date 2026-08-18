"""
F-QC-01 / v1.21.0 — Unified 6-criterion spectrum_qc aggregator.

Assembles the top-level ``diagnostics.spectrum_qc`` block in report.json
from all 6 quality-control criteria defined in
``audit/_drafts/spectrum_qc_methodology_v2_2026-06-03.md``.

Criteria (methodology v2, §Contents table):
    1. Energy calibration drift   — seven_line_check.max_residual_keV ≤ threshold
    2. Peak FWHM stability        — fwhm_at_661 within ±15% of reference
    3. Efficiency curve QA        — efficiency_curve loaded (binary gate v1.21)
    4. Background drift (bg_drift)— rate-normalised z-test (bg_z_test_rates)
                                    or legacy integer z when live-times equal
    5. Per-peak ROI z-test        — already wired by BUG-35; pass-through from
                                    bg_quality_check in StagedAnalysisResult
    6. Sensitivity quarterly      — placeholder (None) in v1.21.0;
                                    Phase 2 RC follow-up with live data

``overall_passed`` logic: AND of all criteria whose data is available.
Unavailable / null criteria count as PASS (conservative: don't false-reject
when data is missing).

Backward compat: existing ``bg_quality_check`` field in ``StagedAnalysisResult``
remains untouched — its data is read here and re-exposed inside spectrum_qc.

RAG-ID: [F-QC-01], [RAG-041]
Methodology cite: spectrum_qc_methodology_v2_2026-06-03.md §criteria 1-6
Rate z-test cite: Gilmore & Joss 2nd Ed. §5.5 (eq. 5.21), F-157
BUG-35 / per-peak z-test: bg_control.py:351-438, RAG-022
"""
from __future__ import annotations

import math
import warnings
from typing import Any, Dict, List, Optional

from gamma.io.bg_control import (
    Z_TIER_BORDERLINE_MAX,
    Z_TIER_STABLE_MAX,
    bg_z_test_rates,
)

# ─── Criterion thresholds (hard defaults; can be overridden by caller) ───
_ENERGY_DRIFT_THRESHOLD_KEV: float = 1.0   # methodology v2 §crit 7 (HPGe ≤0.1%)
_FWHM_REL_THRESHOLD: float = 0.15         # ±15% relative FWHM deviation

# Per-detector-class FWHM reference at 661 keV.
# Sources:
#   NaI  47.0 keV — fwhm_keV_at_energy(_DEFAULT_NAI_FWHM_MODEL, 661.66) = 46.95 keV
#                   staged_pipeline.py:485 _DEFAULT_NAI_FWHM_MODEL=(0.0,2.950048,0.000576400)
#                   ≈ 7.1% resolution; LSRM-9.4 §3.2; RAG-043 Gilmore §6.4.
#   HPGe  1.5 keV — ≈ 0.23% resolution at 661 keV (Gilmore §6.4, RAG-043).
#   LaBr 20.0 keV — ≈ 3.0% resolution at 661 keV (manufacturer spec; RAG-043).
#   default 47.0 — NaI-grade fallback for unknown detector classes.
_FWHM_REFERENCE_BY_DETECTOR: dict = {
    "NaI":     47.0,
    "HPGe":     1.5,
    "LaBr":    20.0,
    "default": 47.0,
}


def _fwhm_reference_keV(detector_class: str) -> float:
    """Return reference FWHM (keV) at 661 keV for the given detector class.

    Parameters
    ----------
    detector_class : str
        One of ``"NaI"``, ``"HPGe"``, ``"LaBr"``.  Unknown values fall back
        to ``"default"`` (47.0 keV, NaI-grade).

    Returns
    -------
    float
        Reference FWHM in keV.

    Sources
    -------
    NaI 47.0 keV: staged_pipeline.py:485 ``_DEFAULT_NAI_FWHM_MODEL``
        → ``fwhm_keV_at_energy(model, 661.66) ≈ 46.95 keV``; LSRM-9.4 §3.2;
        RAG-043 (Gilmore & Joss §6.4).
    HPGe 1.5 keV: RAG-043 Gilmore & Joss §6.4 (~0.23% at 661 keV).
    LaBr 20.0 keV: manufacturer spec ~3.0% at 661 keV; RAG-043.
    """
    return _FWHM_REFERENCE_BY_DETECTOR.get(
        str(detector_class), _FWHM_REFERENCE_BY_DETECTOR["default"]
    )


def _safe_float(x) -> Optional[float]:
    """None-safe float; None for NaN/inf."""
    if x is None:
        return None
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    if math.isnan(v) or math.isinf(v):
        return None
    return v


# ════════════════════════════════════════════════════════════════════
# Criterion 1 — Energy calibration drift
# Methodology v2 §criterion 7: ΔE/E check via seven-line anchor residuals.
# Threshold: max residual ≤ 1.0 keV (conservative for NaI; HPGe would be 0.1%).
# Source: IAEA TRS-295, EURATOM gamma spectrometry QA guidance (cited in
# methodology v2 §criterion 7 extended cite-list).
# ════════════════════════════════════════════════════════════════════

def _criterion_energy_drift(
    result,
    threshold_keV: float = _ENERGY_DRIFT_THRESHOLD_KEV,
) -> Dict[str, Any]:
    """Criterion 1: energy calibration drift via seven_line_check residuals.

    Returns dict with keys: available (bool), max_residual_keV (float|None),
    threshold_keV (float), passed (bool), note (str).

    Cite: methodology v2 §criterion 7 (ΔE/E); IAEA TRS-295.
    """
    slc = getattr(result, "seven_line_check", None)
    if slc is None:
        return {
            "available": False,
            "max_residual_keV": None,
            "threshold_keV": threshold_keV,
            "passed": True,   # unavailable → conservative pass
            "note": "seven_line_check not available; criterion skipped (pass)",
        }
    max_res = _safe_float(getattr(slc, "max_residual_keV", None))
    if max_res is None:
        return {
            "available": False,
            "max_residual_keV": None,
            "threshold_keV": threshold_keV,
            "passed": True,
            "note": "seven_line_check.max_residual_keV is None; criterion skipped (pass)",
        }
    passed = max_res <= threshold_keV
    return {
        "available": True,
        "max_residual_keV": max_res,
        "threshold_keV": threshold_keV,
        "passed": passed,
        "note": (
            f"max_residual={max_res:.3f} keV ≤ {threshold_keV:.1f} keV → stable"
            if passed else
            f"max_residual={max_res:.3f} keV > {threshold_keV:.1f} keV → drift detected"
        ),
    }


# ════════════════════════════════════════════════════════════════════
# Criterion 2 — Peak FWHM stability
# Methodology v2 §criterion 1 (σ_N = √N counting stat) + §criterion 2
# (relative uncertainty). FWHM stability is operational QA gate.
# Source: IAEA AQ-48 (practical QA/QC reference, criteria 2 cite).
# ════════════════════════════════════════════════════════════════════

def _criterion_fwhm_stability(
    result,
    rel_threshold: float = _FWHM_REL_THRESHOLD,
    reference_keV: Optional[float] = None,
    detector_class: str = "default",
) -> Dict[str, Any]:
    """Criterion 2: FWHM at 661 keV within ±15% of reference.

    Parameters
    ----------
    result : StagedAnalysisResult
        Must have ``fwhm_at_661`` and optionally ``fwhm_model_source``.
    rel_threshold : float
        Max allowed relative FWHM deviation (default 0.15 = 15%).
    reference_keV : float, optional
        Explicit reference FWHM in keV.  When *None* (default), the
        reference is resolved via ``_fwhm_reference_keV(detector_class)``.
        Passing an explicit value preserves backward compatibility.
    detector_class : str
        Detector class for automatic reference lookup when ``reference_keV``
        is None.  One of ``"NaI"``, ``"HPGe"``, ``"LaBr"``; unknown values
        use ``"default"`` (47.0 keV, NaI-grade).

    Returns dict with keys: available, fwhm_at_661_keV, reference_keV,
    rel_deviation, rel_threshold, passed, note.

    Cite: methodology v2 §criterion 2; IAEA AQ-48; RAG-043.
    """
    if reference_keV is None:
        reference_keV = _fwhm_reference_keV(detector_class)
    fwhm = _safe_float(getattr(result, "fwhm_at_661", None))
    fwhm_source = getattr(result, "fwhm_model_source", "") or ""

    if fwhm is None:
        return {
            "available": False,
            "fwhm_at_661_keV": None,
            "reference_keV": reference_keV,
            "rel_deviation": None,
            "rel_threshold": rel_threshold,
            "fwhm_model_source": fwhm_source,
            "passed": True,   # unavailable → conservative pass
            "note": "fwhm_at_661 not available; criterion skipped (pass)",
        }
    if reference_keV <= 0:
        return {
            "available": False,
            "fwhm_at_661_keV": fwhm,
            "reference_keV": reference_keV,
            "rel_deviation": None,
            "rel_threshold": rel_threshold,
            "fwhm_model_source": fwhm_source,
            "passed": True,
            "note": "reference_keV=0; criterion skipped (pass)",
        }
    rel_dev = abs(fwhm - reference_keV) / reference_keV
    passed = rel_dev <= rel_threshold
    return {
        "available": True,
        "fwhm_at_661_keV": fwhm,
        "reference_keV": reference_keV,
        "rel_deviation": _safe_float(rel_dev),
        "rel_threshold": rel_threshold,
        "fwhm_model_source": fwhm_source,
        "passed": passed,
        "note": (
            f"FWHM={fwhm:.2f} keV, ref={reference_keV:.2f} keV, "
            f"rel_dev={rel_dev:.3f} ≤ {rel_threshold:.2f} → stable"
            if passed else
            f"FWHM={fwhm:.2f} keV, ref={reference_keV:.2f} keV, "
            f"rel_dev={rel_dev:.3f} > {rel_threshold:.2f} → FWHM drift"
        ),
    }


# ════════════════════════════════════════════════════════════════════
# Criterion 3 — Efficiency curve QA (binary gate v1.21)
# Methodology v2 implied by criteria 3 (detection limit consistency)
# and F-EFF-02 roadmap. Full residual-based QA deferred to F-EFF-02.
# Source: IAEA AQ-48 criteria 3 (efficiency QA).
# ════════════════════════════════════════════════════════════════════

def _criterion_efficiency_qa(result) -> Dict[str, Any]:
    """Criterion 3: efficiency curve loaded (binary gate for v1.21.0).

    Full residual-based QA is deferred to F-EFF-02 (Phase 2 RC).
    Returns dict with keys: available, efficiency_loaded, efficiency_source,
    passed, note.

    Note: the ``note`` field uses only the leaf filename (not the full path)
    to avoid leaking absolute paths or serial numbers into the report text —
    F-115 anonymization compliance.

    Cite: methodology v2 (criteria 3/8), IAEA AQ-48 criteria 3.
    """
    eff_curve = getattr(result, "efficiency_curve", None)
    eff_source = getattr(result, "efficiency_source", "") or ""
    loaded = eff_curve is not None
    return {
        "available": True,
        "efficiency_loaded": loaded,
        "efficiency_source": eff_source,   # full path — _scrub_dict_paths anonymizes it
        "passed": loaded,
        # F-115: note must NOT include path/filename tokens (not scrubbed by anonymize).
        "note": (
            "efficiency_curve loaded"
            if loaded else
            "efficiency_curve not loaded → activities/MDA unreliable"
        ),
    }


# ════════════════════════════════════════════════════════════════════
# Criterion 4 — Background drift F-test / z-test
# Methodology v2 §criterion 5 (ПОЛНОСТЬЮ ПЕРЕПИСАН section):
#   B₁, B₂ — интегральные счёты фона в ROI за ОДИНАКОВОЕ live time
#   Rate form: when t1 ≠ t2 → bg_z_test_rates(c1,t1,c2,t2)
#   Cite: Gilmore & Joss §5.5, F-157, methodology v2 §5
#   Tiers: |z|<2 stable, 2-3 borderline, ≥3 reject (methodology v2 §crit 5)
#
# P0-6 wiring: bg_z_test_rates applies when sample_live_time ≠ bg_live_time.
# When bg is not available, criterion is skipped (conservative pass).
# ════════════════════════════════════════════════════════════════════

def _criterion_bg_drift(result) -> Dict[str, Any]:
    """Criterion 4: background drift rate-normalised z-test.

    Uses bg_z_test_rates(c1,t1,c2,t2) from bg_control.py:351-438 for
    unequal live-times. Falls back to bg_quality_check overall pass when
    equal live-times (BUG-35 path already handles that).

    Returns dict with keys: available, method, z, is_significant,
    sample_live_time_s, bg_live_time_s, sample_sum_counts, bg_sum_counts,
    passed, note.

    Cite: Gilmore & Joss 2nd Ed. §5.5, methodology v2 §criterion 5.
    """
    bg_sub = getattr(result, "background_subtraction", None)
    if bg_sub is None:
        return {
            "available": False,
            "method": None,
            "z": None,
            "is_significant": None,
            "sample_live_time_s": None,
            "bg_live_time_s": None,
            "sample_sum_counts": None,
            "bg_sum_counts": None,
            "passed": True,   # no background → skip criterion, conservative pass
            "note": "no background_subtraction; criterion skipped (pass)",
        }

    spec = getattr(result, "spec", None)
    sample_t = _safe_float(getattr(spec, "live_time", None) if spec else None)
    bg_t = _safe_float(getattr(bg_sub, "bg_live_time", None))
    sample_c = _safe_float(getattr(bg_sub, "sample_sum_counts", None))
    bg_c = _safe_float(getattr(bg_sub, "bg_sum_counts", None))

    # If counts or live-times not available, fall back to legacy gate result
    if None in (sample_t, bg_t, sample_c, bg_c) or sample_t <= 0 or bg_t <= 0:
        # Try to recover from existing bg_quality_check
        bqc = getattr(result, "bg_quality_check", None)
        if bqc is not None:
            legacy_passed = bool(bqc.get("overall_passed", True))
            return {
                "available": True,
                "method": "legacy_bg_quality_check",
                "z": None,
                "is_significant": None,
                "sample_live_time_s": sample_t,
                "bg_live_time_s": bg_t,
                "sample_sum_counts": sample_c,
                "bg_sum_counts": bg_c,
                "passed": legacy_passed,
                "note": "live-time or count data incomplete; using legacy bg_quality_check",
            }
        return {
            "available": False,
            "method": None,
            "z": None,
            "is_significant": None,
            "sample_live_time_s": sample_t,
            "bg_live_time_s": bg_t,
            "sample_sum_counts": sample_c,
            "bg_sum_counts": bg_c,
            "passed": True,
            "note": "live-time or count data unavailable; criterion skipped (pass)",
        }

    # Rate-normalised z-test (Gilmore & Joss §5.5)
    z, is_sig = bg_z_test_rates(sample_c, sample_t, bg_c, bg_t)

    # Tier assignment per methodology v2 §criterion 5
    if math.isnan(z):
        tier = "undefined"
        passed = True   # degenerate input — skip
        note = "bg_z_test_rates returned NaN (degenerate input); criterion skipped (pass)"
    elif abs(z) < Z_TIER_STABLE_MAX:
        tier = "stable"
        passed = True
        note = f"|z|={abs(z):.2f} < {Z_TIER_STABLE_MAX:.1f} → bg drift stable"
    elif abs(z) < Z_TIER_BORDERLINE_MAX:
        tier = "borderline"
        passed = True   # borderline = monitor, not reject
        note = f"|z|={abs(z):.2f} in [{Z_TIER_STABLE_MAX:.1f},{Z_TIER_BORDERLINE_MAX:.1f}) → monitor"
    else:
        tier = "reject"
        passed = False
        note = f"|z|={abs(z):.2f} ≥ {Z_TIER_BORDERLINE_MAX:.1f} → systematic bg drift"

    return {
        "available": True,
        "method": "rate_normalised",
        "z": _safe_float(z),
        "is_significant": is_sig,
        "tier": tier,
        "sample_live_time_s": sample_t,
        "bg_live_time_s": bg_t,
        "sample_sum_counts": sample_c,
        "bg_sum_counts": bg_c,
        "passed": passed,
        "note": note,
    }


# ════════════════════════════════════════════════════════════════════
# Criterion 5 — Per-peak ROI z-test (BUG-35 pass-through)
# Already wired in staged_pipeline.py:1126-1179 via bg_quality_check.
# Re-exposed here in the unified spectrum_qc block for F-QC-01.
# ════════════════════════════════════════════════════════════════════

def _criterion_peak_z_roi(result) -> tuple[Optional[List], int, int, int]:
    """Criterion 5: extract per-peak ROI z-test from existing bg_quality_check.

    Returns (peak_z_roi_list, n_peaks_tested, n_passed, n_failed).
    When bg_quality_check is None: returns ([], 0, 0, 0).

    Cite: BUG-35, RAG-022, methodology v2 §criterion 5.
    """
    bqc = getattr(result, "bg_quality_check", None)
    if bqc is None:
        return [], 0, 0, 0
    peak_z_roi = bqc.get("peak_z_roi", []) or []
    n_tested = int(bqc.get("n_peaks_tested", len(peak_z_roi)))
    n_passed = int(bqc.get("n_passed", sum(1 for e in peak_z_roi if e.get("passed"))))
    n_failed = int(bqc.get("n_failed", n_tested - n_passed))
    return peak_z_roi, n_tested, n_passed, n_failed


# ════════════════════════════════════════════════════════════════════
# Criterion 6 — Sensitivity quarterly (placeholder v1.21.0)
# Phase 2 RC follow-up: requires live quarterly detector efficiency data.
# Methodology v2 §criterion 8 (Detection limit consistency).
# ════════════════════════════════════════════════════════════════════

def _criterion_sensitivity(_result) -> Optional[Dict[str, Any]]:
    """Criterion 6: sensitivity quarterly (placeholder in v1.21.0).

    Returns None in v1.21.0. Phase 2 RC follow-up requires live
    detector efficiency trend data for the quarterly sensitivity check.

    Cite: methodology v2 §criterion 8 (Detection limit consistency),
    ISO 11929-4 (Repeated measurements / uncertainty consistency).
    """
    # Phase 2 RC TODO: query detector efficiency trend DB, compute
    # relative sensitivity drift vs rolling 90-day baseline.
    return None


# ════════════════════════════════════════════════════════════════════
# Top-level aggregator
# ════════════════════════════════════════════════════════════════════

def build_spectrum_qc(
    result,
    *,
    energy_drift_threshold_keV: float = _ENERGY_DRIFT_THRESHOLD_KEV,
    fwhm_rel_threshold: float = _FWHM_REL_THRESHOLD,
    fwhm_reference_keV: Optional[float] = None,
    detector_class: str = "default",
) -> Dict[str, Any]:
    """Build the unified ``spectrum_qc`` block for report.json.

    Aggregates all 6 QC criteria from methodology v2 into a single dict.
    ``overall_passed`` is the AND of all available (non-skipped) criteria.

    Parameters
    ----------
    result : StagedAnalysisResult
        Pipeline result object. Must have attributes:
        ``seven_line_check``, ``fwhm_at_661``, ``fwhm_model_source``,
        ``efficiency_curve``, ``efficiency_source``,
        ``background_subtraction`` (optional), ``bg_quality_check``
        (optional, from BUG-35), ``spec``.
    energy_drift_threshold_keV : float
        Max allowed energy residual in keV (criterion 1). Default 1.0.
    fwhm_rel_threshold : float
        Max allowed relative FWHM deviation (criterion 2). Default 0.15.
    fwhm_reference_keV : float, optional
        Explicit reference FWHM at 661 keV for criterion 2.  When *None*
        (default), the reference is resolved via
        ``_fwhm_reference_keV(detector_class)``.  Pass an explicit value
        for backward-compatible call-sites that supply a fixed reference.
    detector_class : str
        Detector class for automatic FWHM reference lookup.  One of
        ``"NaI"``, ``"HPGe"``, ``"LaBr"``; unknown → ``"default"``
        (47.0 keV, NaI-grade).  Used only when ``fwhm_reference_keV``
        is *None*.  BUG-48 fix: was hard-coded to 8.5 keV (wrong for NaI).

    Returns
    -------
    dict
        Keys: ``n_peaks_tested``, ``n_passed``, ``n_failed``,
        ``overall_passed``, ``peak_z_roi``, ``energy_drift``,
        ``fwhm_stability``, ``efficiency_qa``, ``bg_drift``,
        ``sensitivity``.

    Notes
    -----
    Backward compat: ``result.bg_quality_check`` is not removed or
    modified; its data is re-exposed in ``peak_z_roi`` + criterion counts.

    F-QC-01 cite: KNOWN_AND_FIXED_ISSUES.md:1292 (original definition),
    PLAN_v1_20_to_v1_21.md §P0-5 lines 124-148.
    """
    # Criterion 1 — energy drift
    energy_drift = _criterion_energy_drift(
        result, threshold_keV=energy_drift_threshold_keV
    )

    # Criterion 2 — FWHM stability
    # BUG-48: fwhm_reference_keV=None triggers per-detector lookup via detector_class.
    # Explicit fwhm_reference_keV still accepted for backward compat.
    fwhm_stability = _criterion_fwhm_stability(
        result,
        rel_threshold=fwhm_rel_threshold,
        reference_keV=fwhm_reference_keV,
        detector_class=detector_class,
    )

    # Criterion 3 — efficiency QA
    efficiency_qa = _criterion_efficiency_qa(result)

    # Criterion 4 — background drift
    bg_drift = _criterion_bg_drift(result)

    # Criterion 5 — per-peak ROI z (pass-through from BUG-35)
    peak_z_roi, n_tested, n_passed_peaks, n_failed_peaks = _criterion_peak_z_roi(result)

    # Criterion 6 — sensitivity (placeholder)
    sensitivity = _criterion_sensitivity(result)

    # overall_passed: AND of available criteria
    # Unavailable criteria (.passed=True by convention) do not false-reject.
    criteria_passed = [
        energy_drift["passed"],
        fwhm_stability["passed"],
        efficiency_qa["passed"],
        bg_drift["passed"],
        # Criterion 5: peak-level overall_passed from bg_quality_check
        (n_failed_peaks == 0),   # zero failures when no peaks tested = pass
        # Criterion 6: placeholder → always pass
        True,
    ]
    overall_passed = all(criteria_passed)

    return {
        "n_peaks_tested": n_tested,
        "n_passed": n_passed_peaks,
        "n_failed": n_failed_peaks,
        "overall_passed": overall_passed,
        "peak_z_roi": peak_z_roi,
        "energy_drift": energy_drift,
        "fwhm_stability": fwhm_stability,
        "efficiency_qa": efficiency_qa,
        "bg_drift": bg_drift,
        "sensitivity": sensitivity,
    }


__all__ = [
    "build_spectrum_qc",
    "_fwhm_reference_keV",
    "_FWHM_REFERENCE_BY_DETECTOR",
    "_criterion_energy_drift",
    "_criterion_fwhm_stability",
    "_criterion_efficiency_qa",
    "_criterion_bg_drift",
    "_criterion_peak_z_roi",
    "_criterion_sensitivity",
]
