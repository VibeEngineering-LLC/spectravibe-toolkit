"""
bg_carryover — annotate sample-side peaks that originate from the background.

After `subtract_background()` we run identification on **both** the sample
(net) spectrum and the background spectrum separately. Lines like K-40 1461,
Tl-208 2614 or Bi-214 609 routinely sit in both: they live in the bg spectrum
(walls / Marinelli vessel material) and they survive in the net after
subtraction as residuals. Without an explicit marker the operator sees them
in the "detected in sample" table and has to mentally cross-reference the
parallel "background peaks" block to know which lines are really sample-side
vs which are bg residuals.

This module wires that cross-reference into the data: for each sample-side
peak it looks for a matching bg-side peak (same nuclide stem, energy within
±N·FWHM(E)) and attaches a compact `bg_carryover` dict to the sample entry.
JSON / HTML / Markdown renderers consume that field and emit a «фон» pill or
note.

Design contract:
- Pure annotation. Does NOT change activities, does NOT subtract anything
  twice. The bg counts are already gone — we just label what's left.
- Reuses already-computed `bg_staged_result.final_detected`. No new peak
  search runs.
- Energy-only match. Intensity / residual_ratio is intentionally not gated
  on — the goal is to surface origin, not to second-guess subtraction
  quality. (`S_bg_counts` is still emitted as informational so the operator
  can eyeball magnitude.)

API:
    catalog = build_bg_peak_catalog(bg_staged_result)
    # fit_overlay.peaks  (energy_keV, source filter on singlets/multiplet):
    mark_bg_carryover(peaks_out, catalog, fwhm_e_keV)
    # primary_feps rows  (peak_E_keV, no source filter):
    mark_bg_carryover(primary_feps, catalog, fwhm_e_keV,
                      sample_sources=None, energy_field="peak_E_keV")

`peaks_out` is mutated in place: matching sample-side entries gain a
`"bg_carryover": {...}` key. Non-matching entries are untouched.
"""
from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional


_BG_NUC_SUFFIX = " (bg)"
"""Suffix appended by json_report._build_fit_overlay when emitting bg-side
entries (label format `f"{nuc} (bg)"`)."""


def _strip_bg_suffix(nuclide: str) -> str:
    """Return the nuclide stem without the trailing ' (bg)' marker."""
    if nuclide.endswith(_BG_NUC_SUFFIX):
        return nuclide[: -len(_BG_NUC_SUFFIX)]
    return nuclide


def build_bg_peak_catalog(bg_staged_result: Any) -> List[Dict[str, Any]]:
    """Flatten a staged IdentifyResult into a list of bg peak descriptors.

    Each entry: ``{"nuclide": "K-40", "E_keV": 1461.13, "S_counts": 12.5,
    "peak_id": "pb1461"}``.

    Returns ``[]`` if the input is falsy / has no ``final_detected``.
    """
    if bg_staged_result is None:
        return []
    detected = getattr(bg_staged_result, "final_detected", None) or []
    catalog: List[Dict[str, Any]] = []
    for ni in detected:
        nuc = getattr(ni, "nuclide", None)
        if not nuc:
            continue
        for m in (getattr(ni, "matched_lines", None) or []):
            e = getattr(m, "peak_E_keV", None)
            if e is None:
                continue
            try:
                e_f = float(e)
            except (TypeError, ValueError):
                continue
            area = getattr(m, "peak_area", None)
            try:
                area_f = float(area) if area is not None else 0.0
            except (TypeError, ValueError):
                area_f = 0.0
            if area_f <= 0:
                continue
            catalog.append({
                "nuclide": str(nuc),
                "E_keV": e_f,
                "S_counts": area_f,
                "peak_id": f"pb{round(e_f)}",
            })
    return catalog


def mark_bg_carryover(
    peaks_out: List[Dict[str, Any]],
    bg_catalog: List[Dict[str, Any]],
    fwhm_e_keV: Optional[Callable[[float], float]],
    *,
    window_n_fwhm: float = 1.5,
    sample_sources: Optional[tuple] = ("singlet", "multiplet_component"),
    energy_field: str = "energy_keV",
    nuclide_field: str = "nuclide",
) -> int:
    """Annotate sample-side peaks that match a bg-side peak.

    Match criteria: identical nuclide stem (after stripping the `(bg)`
    suffix from bg entries) AND ``|E_sample - E_bg| ≤ window_n_fwhm ·
    FWHM(E_sample)``. When ``fwhm_e_keV`` is None we fall back to a flat
    3 keV window (NaI-class default).

    Mutates ``peaks_out`` in place. Returns the number of annotated peaks.
    """
    if not peaks_out or not bg_catalog:
        return 0

    by_nuclide: Dict[str, List[Dict[str, Any]]] = {}
    for bg in bg_catalog:
        by_nuclide.setdefault(_strip_bg_suffix(bg["nuclide"]), []).append(bg)

    annotated = 0
    for entry in peaks_out:
        if sample_sources is not None and entry.get("source") not in sample_sources:
            continue
        nuc = entry.get(nuclide_field)
        if not nuc:
            continue
        candidates = by_nuclide.get(_strip_bg_suffix(str(nuc)))
        if not candidates:
            continue
        e_sample = entry.get(energy_field)
        if e_sample is None:
            continue
        try:
            e_sample_f = float(e_sample)
        except (TypeError, ValueError):
            continue
        if fwhm_e_keV is not None:
            try:
                fwhm = float(fwhm_e_keV(e_sample_f))
            except Exception:
                fwhm = 3.0
            if fwhm <= 0:
                fwhm = 3.0
        else:
            fwhm = 3.0
        window = window_n_fwhm * fwhm

        best = None
        best_delta = float("inf")
        for bg in candidates:
            delta = abs(bg["E_keV"] - e_sample_f)
            if delta <= window and delta < best_delta:
                best = bg
                best_delta = delta
        if best is None:
            continue
        entry["bg_carryover"] = {
            "matched": True,
            "bg_peak_id": best["peak_id"],
            "E_bg_keV": round(best["E_keV"], 2),
            "delta_E_keV": round(best_delta, 2),
            "S_bg_counts": round(best["S_counts"], 2),
        }
        annotated += 1
    return annotated


__all__ = [
    "build_bg_peak_catalog",
    "mark_bg_carryover",
]