"""
F-445 / v1.30.3 - Cluster-level Delta_cluster anchors for E-calibration.

Operator-locked contract (2026-06-11, HARD):
    "Do calibration peaks may sit off their library positions, but spacing
    between centroids must match the library spacing modulo calibration
    nonlinearity" (operator wording, paraphrased to ASCII).

F-445 idea:
  * Old F-145 Phase B/C used per-component I_pct-weighted shifts from
    Phase A free_centroids. On close lines this can swap centroids
    (Ac-228 964.77 / 968.97 in M1 on Th-232, observed 2026-06-11).
  * Delta_cluster = one shift for the whole cluster, chosen so that the
    strongest component sits on local max of (counts - continuum) in the
    +-0.7 * FWHM_lib window. Spacing is preserved, operator contract holds.
  * v1.29.1 (`json_report._compute_cluster_global_shift`) lived only in
    the overlay layer; data did not move, sum-curve did.
  * F-445 promotes Delta_cluster to the *real* E-cal: anchors
    (E_strongest_lib, channel(E_strongest_lib + Delta)) go into Phase C
    polyfit, shifting the spectrum scale.

This module:
  * `compute_cluster_global_shift` - public port from render layer.
  * `collect_cluster_global_anchors` - per-cluster anchor builder for
    F-145 Phase B/C.

See also:
  * LSRM Algorithmic Foundations sec.8.4.4 (multiplet self-calibration).
  * Gilmore & Joss 3rd Ed. sec.6.4 (energy calibration drift).
  * F-145 - `gamma.calibration.multiplet_self_calibration`.
"""
from __future__ import annotations

from typing import Callable, Iterable, List, Optional, Any


SQRT_2PI = 2.5066282746310002


def compute_cluster_global_shift(
    E_lib_strongest,
    sigma_strongest,
    E_arr_raw,
    cont_arr_raw,
    counts_arr,
    spec,
):
    """Return Delta (keV) such that the strongest peak sits on local max of
    (counts - continuum) in window +-0.7 * FWHM_lib.

    Returns None when:
      * inputs are missing / inconsistent;
      * no positive net signal in the window (best_net <= 0);
      * |Delta| > 1.5 * FWHM_lib (runaway clamp).

    Verbatim port from `gamma.reporting.json_report._compute_cluster_global_shift`
    v1.29.1 (operator-locked render). Preserves all guards and numerical
    semantics so downstream tests/results stay bit-identical when invoked
    with the same inputs.
    """
    E_arr_list = list(E_arr_raw) if E_arr_raw is not None else []
    cont_arr_list = list(cont_arr_raw) if cont_arr_raw is not None else []
    if counts_arr is None or spec is None or not E_arr_list:
        return None
    try:
        fwhm_lib = max(float(sigma_strongest) * 2.355, 1.0)
    except Exception:
        return None
    win = max(0.7 * fwhm_lib, 8.0)
    n = len(E_arr_list)
    if n == 0 or len(cont_arr_list) != n:
        return None

    def _cont_at(E):
        if E <= E_arr_list[0]:
            return float(cont_arr_list[0])
        if E >= E_arr_list[-1]:
            return float(cont_arr_list[-1])
        lo, hi = 0, n - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if E_arr_list[mid] <= E:
                lo = mid
            else:
                hi = mid
        x0, x1 = E_arr_list[lo], E_arr_list[hi]
        y0, y1 = cont_arr_list[lo], cont_arr_list[hi]
        if x1 <= x0:
            return float(y0)
        t = (E - x0) / (x1 - x0)
        return float(y0 + (y1 - y0) * t)

    try:
        ch_lo_raw = spec.energy_to_channel(E_lib_strongest - win)
        ch_hi_raw = spec.energy_to_channel(E_lib_strongest + win)
        if ch_lo_raw is None or ch_hi_raw is None:
            return None
        ch_lo = max(0, int(ch_lo_raw))
        ch_hi = min(len(counts_arr), int(ch_hi_raw) + 1)
        if ch_hi - ch_lo < 5:
            return None
        best_E = float(E_lib_strongest)
        best_net = -1e30
        for ch in range(ch_lo, ch_hi):
            try:
                E_ch = float(spec.channel_to_energy(ch))
            except Exception:
                continue
            net = float(counts_arr[ch]) - _cont_at(E_ch)
            if net > best_net:
                best_net = net
                best_E = E_ch
        if best_net <= 0:
            return None
        delta = best_E - float(E_lib_strongest)
        if abs(delta) > 1.5 * fwhm_lib:
            return None
        return float(delta)
    except Exception:
        return None


def collect_cluster_global_anchors(
    spec,
    forced_clusters: Iterable,
    counts_arr,
    continuum_arrays: Optional[dict],
    E_arr,
    fwhm_provider_keV: Callable[[float], float],
    *,
    phantom_amp_sigma_floor: float = 1.0,
) -> List[Any]:
    """Build per-cluster Delta_cluster anchors for F-445 Phase B/C.

    Parameters
    ----------
    spec : Spectrum
        Sample spectrum (used for energy_to_channel / channel_to_energy).
    forced_clusters : iterable of DeconvolutionResult
        Phase A coupled-fit results. Each must carry `components`, `areas`,
        `overlay_E_keV`, `overlay_continuum` (rendered overlay arrays).
    counts_arr : array-like
        Per-channel counts (preferably netto, see staged_pipeline wiring).
    continuum_arrays : dict[cluster_id -> (E_arr_list, cont_arr_list)] or None
        Optional override. If None, falls back to cluster.overlay_E_keV /
        cluster.overlay_continuum.
    E_arr : array-like
        Per-channel energy axis. Currently unused (continuum interpolation
        uses overlay_E_keV in keV domain). Kept for API symmetry with
        caller-side wiring.
    fwhm_provider_keV : Callable[[float], float]
        FWHM(E) in keV.

    Returns
    -------
    List[CentroidAnchor] - one per cluster where Delta_cluster computed.

    Phantom guard: if max area / sqrt(2*pi) <= phantom_amp_sigma_floor,
    the cluster is skipped (matches render-layer guard).
    """
    # Lazy import to avoid module-load cycle
    from gamma.calibration.multiplet_self_calibration import CentroidAnchor

    anchors: List[Any] = []
    for ci, cluster in enumerate(forced_clusters):
        comps = list(getattr(cluster, "components", []) or [])
        areas = list(getattr(cluster, "areas", []) or [])
        if not comps or len(areas) != len(comps):
            continue
        # Strongest = max area (sigma cancels in amp*sigma = area/sqrt(2*pi))
        strongest_idx = -1
        max_area = -1.0
        for k, a in enumerate(areas):
            try:
                af = float(a or 0.0)
            except Exception:
                af = 0.0
            if af > max_area:
                max_area = af
                strongest_idx = k
        if strongest_idx < 0 or max_area <= 0:
            continue
        comp_strong = comps[strongest_idx]
        E_lib_strong_raw = (
            getattr(comp_strong, "line_E_keV", None)
            if getattr(comp_strong, "line_E_keV", None) is not None
            else getattr(comp_strong, "E_keV", 0.0)
        )
        try:
            E_lib_strong = float(E_lib_strong_raw or 0.0)
        except Exception:
            continue
        if E_lib_strong <= 0:
            continue
        # sigma_strongest in keV - defines +-0.7 * FWHM_lib search window.
        fwhm_lib = max(float(fwhm_provider_keV(E_lib_strong) or 0.0), 1e-6)
        sigma_strong_keV = fwhm_lib / 2.355
        # Phantom guard: amp_strong * sigma_strong = area_strong / sqrt(2*pi).
        amp_sigma_strong = max_area / SQRT_2PI
        if amp_sigma_strong <= phantom_amp_sigma_floor:
            continue
        # Continuum arrays
        cid_str = str(getattr(cluster, "cluster_id", "") or "")
        if continuum_arrays is not None and cid_str in continuum_arrays:
            E_arr_raw, cont_arr_raw = continuum_arrays[cid_str]
        else:
            E_arr_raw = list(getattr(cluster, "overlay_E_keV", ()) or ())
            cont_arr_raw = list(getattr(cluster, "overlay_continuum", ()) or ())
        if not E_arr_raw or not cont_arr_raw or len(E_arr_raw) != len(cont_arr_raw):
            continue
        delta = compute_cluster_global_shift(
            E_lib_strong, sigma_strong_keV,
            E_arr_raw, cont_arr_raw, counts_arr, spec,
        )
        if delta is None:
            continue
        try:
            ch_fit = float(spec.energy_to_channel(E_lib_strong + delta))
        except Exception:
            continue
        nuc = str(getattr(comp_strong, "nuclide", "") or "")
        cid = cid_str or ("M" + str(ci))
        source_label = (
            "cluster_delta_" + cid + "_" + nuc + "_"
            + str(int(round(E_lib_strong)))
        )
        anchor = CentroidAnchor(
            nuclide=nuc,
            E_passport_keV=E_lib_strong,
            E_fitted_keV=E_lib_strong + delta,
            channel_fitted=ch_fit,
            fwhm_keV=fwhm_lib,
            drift_fraction_of_fwhm=abs(delta) / max(fwhm_lib, 1e-6),
            source=source_label,
        )
        anchors.append(anchor)
    return anchors


__all__ = [
    "compute_cluster_global_shift",
    "collect_cluster_global_anchors",
]
