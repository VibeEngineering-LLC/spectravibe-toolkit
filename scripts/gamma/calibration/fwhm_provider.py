"""
FWHM-at-channel provider builder.

Builds a `Callable[[int], float]` that returns the expected peak FWHM
(in channel units) at any channel of a parsed spectrum. The provider
is the single source of truth used by:
  - `gamma.peaks.search.mariscotti_search` (adaptive filter sizing)
  - `gamma.calibration.stored_check.check_stored_calibration`
    (per-anchor match window)
  - future `gamma.calibration.bootstrap` improvements (anchor-pattern
    tolerance scaling)

Source hierarchy (first hit wins):
  1. Vendor-stored FWHM model coefficients
     (e.g. AtomSpectra `SimpleSqrtFwhmCalibration`: FWHM(N) = c0 + c1*N
      - note: NOT `c0 + c1*sqrt(N)` as a naive reading of the name might
      suggest; verified against the embedded `<CalibrationPeaks>` values
      from real fixtures, both peaks reproduce to <0.01 channel)
  2. Linear interpolation across `<CalibrationPeaks>` measured pairs
     (channel, fwhm_channels)
  3. F-449 bootstrap: FWHM(E) fitted from significant, isolated peaks of
     the spectrum itself (only when `bootstrap_from_peaks=True` and the
     stored sources above are unavailable). source-label
     "bootstrap_significant_peaks_*".
  4. Caller-supplied fallback in channels (default 10.0) - last resort.

The returned callable is safe to call at any channel index; out-of-
range channels are clamped to the nearest endpoint of the supported
range. Negative or zero outputs of any underlying model are bounded
to a minimum of 1.0 channel so that downstream kernel-builders cannot
construct a degenerate filter.

F-449 (operator-locked 2026-06-16): the returned callable carries a
`.fwhm_source` attribute naming the source actually used
("stored_model" | "calibration_peaks" | "bootstrap_significant_peaks_*"
 | "fallback") for report provenance.
"""

from __future__ import annotations

import math
from typing import Callable, List, Optional, Tuple


# F-449 (operator-locked 2026-06-16): bootstrap FWHM(E) from significant
# peaks of the spectrum when no vendor FWHM model / calibration peaks are
# usable, instead of returning a generic constant. The generic constant
# remains the last-resort fallback. Provenance of the chosen source is
# exposed on the returned callable via the `.fwhm_source` attribute.
#
# Selection of bootstrap anchor peaks (all must hold):
#   - significance (Currie L_C) >= `_BOOTSTRAP_MIN_SIGMA`
#   - net height >= `_BOOTSTRAP_MIN_HEIGHT` counts (statistics floor for
#     a reliable half-height width read)
#   - isolated: nearest neighbouring found peak is >= `_BOOTSTRAP_ISO_FAC`
#     seed-FWHM away on both sides (otherwise the half-height edge runs
#     into a neighbour and inflates the measured width)
# At least `_BOOTSTRAP_MIN_ANCHORS` clean anchors are required for a
# quadratic FWHM^2 = a + b*E + c*E^2 fit; with 2 anchors a linear
# FWHM^2 = a + b*E fit is used; with fewer, bootstrap is declined and the
# generic constant fallback is used.
_BOOTSTRAP_MIN_SIGMA = 8.0
_BOOTSTRAP_MIN_HEIGHT = 30.0
_BOOTSTRAP_ISO_FAC = 1.5
_BOOTSTRAP_MIN_ANCHORS = 3
_BOOTSTRAP_MIN_ANCHORS_LINEAR = 2


def _measure_fwhm_channels(counts, ch: int, seed_fwhm: float) -> Optional[float]:
    """Measure a peak's FWHM (channels) by half-height crossing.

    Estimates the local continuum from the ROI wings (mean of the outer
    fifth on each side, the same convention as the area-fit wings) and
    finds the half-maximum crossings around `ch`. Returns None when the
    statistics are too poor to read a width (net height <= 0, or the
    half-level is never crossed within the search range).
    """
    n = len(counts)
    if ch < 1 or ch >= n - 1:
        return None
    search_range = max(3, int(round(2.5 * seed_fwhm)))
    lo = max(0, ch - search_range)
    hi = min(n, ch + search_range + 1)
    if hi - lo < 5:
        return None
    wing = max(1, (hi - lo) // 5)
    left_bg = sum(float(counts[i]) for i in range(lo, lo + wing)) / wing
    right_bg = sum(float(counts[i]) for i in range(hi - wing, hi)) / wing
    baseline = 0.5 * (left_bg + right_bg)
    peak_val = float(counts[ch])
    net_height = peak_val - baseline
    if net_height <= 0.0:
        return None
    half_level = baseline + 0.5 * net_height
    # Walk right until below half-level, then linear-interpolate the crossing.
    rr = ch
    while rr < hi - 1 and float(counts[rr]) >= half_level:
        rr += 1
    if float(counts[rr]) >= half_level:
        return None  # never dropped below - width unbounded on right
    y1, y0 = float(counts[rr]), float(counts[rr - 1])
    x_right = (rr - 1) + (half_level - y0) / (y1 - y0) if y1 != y0 else float(rr)
    ll = ch
    while ll > lo and float(counts[ll]) >= half_level:
        ll -= 1
    if float(counts[ll]) >= half_level:
        return None
    y0, y1 = float(counts[ll]), float(counts[ll + 1])
    x_left = ll + (half_level - y0) / (y1 - y0) if y1 != y0 else float(ll)
    fwhm = x_right - x_left
    if fwhm <= 0.0:
        return None
    return float(fwhm)


def _bootstrap_fwhm_anchors(
    spec, *, seed_fwhm: float, min_channels: float,
) -> List[Tuple[float, float]]:
    """Detect significant, isolated peaks and measure their FWHM.

    Returns a list of `(E_keV, fwhm_keV)` anchor pairs suitable for a
    FWHM^2(E) fit. Requires `spec.counts` and a 2+ term `spec.energy_cal`.
    Empty list => bootstrap not possible (caller uses constant fallback).
    """
    counts = getattr(spec, "counts", None)
    if counts is None or len(counts) < 50:
        return []
    if not spec.energy_cal or len(spec.energy_cal) < 2:
        return []
    try:
        import numpy as _np

        from gamma.peaks.search import mariscotti_search
    except Exception:
        return []

    counts_arr = _np.asarray(counts, dtype=float)
    # Scalar seed FWHM => the internal search does NOT call back into this
    # provider, so there is no recursion.
    try:
        found = mariscotti_search(
            counts_arr,
            fwhm_channels=float(seed_fwhm),
            sigma_threshold=_BOOTSTRAP_MIN_SIGMA,
            min_separation_factor=1.0,
        )
    except Exception:
        return []
    if not found:
        return []

    channels_sorted = sorted(p.channel for p in found)
    anchors: List[Tuple[float, float]] = []
    for p in found:
        if p.significance < _BOOTSTRAP_MIN_SIGMA:
            continue
        if p.height < _BOOTSTRAP_MIN_HEIGHT:
            continue
        # Isolation: nearest other found peak must be >= ISO_FAC*seed away.
        iso_dist = _BOOTSTRAP_ISO_FAC * seed_fwhm
        too_close = False
        for q_ch in channels_sorted:
            if q_ch == p.channel:
                continue
            if abs(q_ch - p.channel) < iso_dist:
                too_close = True
                break
        if too_close:
            continue
        fwhm_ch = _measure_fwhm_channels(counts_arr, int(p.channel), seed_fwhm)
        if fwhm_ch is None or fwhm_ch < min_channels:
            continue
        E_keV = spec.channel_to_energy(int(p.channel))
        if E_keV is None or E_keV <= 0:
            continue
        # local gain dE/dN for channel->keV width conversion
        dE_dN = sum(
            i * a * (float(p.channel) ** (i - 1))
            for i, a in enumerate(spec.energy_cal) if i > 0
        )
        if dE_dN <= 0:
            continue
        fwhm_keV = fwhm_ch * dE_dN
        if fwhm_keV <= 0:
            continue
        anchors.append((float(E_keV), float(fwhm_keV)))
    return anchors


def _fit_fwhm_sq_model(
    anchors: List[Tuple[float, float]],
) -> Optional[Tuple[Tuple[float, ...], str]]:
    """Fit FWHM^2(E) = a + b*E (+ c*E^2) to bootstrap anchors.

    Returns `(coeffs, source_label)` or None when fewer than
    `_BOOTSTRAP_MIN_ANCHORS_LINEAR` anchors are available or the linear
    algebra fails. `coeffs` is `(a, b)` for linear or `(a, b, c)` for
    quadratic. The fit is rejected (None) if it produces a non-positive
    FWHM^2 at the lowest anchor energy (pathological / non-physical).
    """
    if len(anchors) < _BOOTSTRAP_MIN_ANCHORS_LINEAR:
        return None
    try:
        import numpy as _np
    except Exception:
        return None
    Es = _np.array([a[0] for a in anchors], dtype=float)
    Fsq = _np.array([a[1] ** 2 for a in anchors], dtype=float)
    if len(anchors) >= _BOOTSTRAP_MIN_ANCHORS:
        A = _np.vstack([_np.ones_like(Es), Es, Es ** 2]).T
        try:
            coefs, *_ = _np.linalg.lstsq(A, Fsq, rcond=None)
        except Exception:
            return None
        coeffs = (float(coefs[0]), float(coefs[1]), float(coefs[2]))
        label = "bootstrap_significant_peaks_quadratic"
    else:
        A = _np.vstack([_np.ones_like(Es), Es]).T
        try:
            coefs, *_ = _np.linalg.lstsq(A, Fsq, rcond=None)
        except Exception:
            return None
        coeffs = (float(coefs[0]), float(coefs[1]))
        label = "bootstrap_significant_peaks_linear"
    # Reject non-physical fit (negative FWHM^2 at the anchor span).
    E_min = float(Es.min())
    val_min = sum(c * (E_min ** i) for i, c in enumerate(coeffs))
    if val_min <= 0:
        return None
    return coeffs, label


def make_fwhm_at_channel_provider(
    spec,
    *,
    fallback_channels: float = 10.0,
    min_channels: float = 1.0,
    bootstrap_from_peaks: bool = False,
) -> Callable[[int], float]:
    """
    Build a callable `fwhm_at_channel(ch) -> fwhm_channels` for `spec`.

    Resolves the FWHM source in this order:

    1. `spec.stored_fwhm_calibration.model == "SimpleSqrtFwhm"` with
       at least two coefficients -> use FWHM(N) = sqrt(c0 + c1*N).
       This is the AtomSpectra `SimpleSqrtFwhmCalibration` model. The
       name is misleading: the "Sqrt" refers to the outer square root
       of an affine combination of N, not to sqrt(N) inside.
    2. `spec.stored_fwhm_calibration.calibration_peaks` with >=1 entry
       carrying a measured `fwhm_channels` -> use piecewise-linear
       interpolation in channel space, with constant extrapolation
       beyond the endpoints. With a single point, returns that point's
       FWHM at every channel.
    3. F-449 (only if `bootstrap_from_peaks=True`): fit FWHM^2(E) from
       significant, isolated peaks measured directly off `spec.counts`,
       then evaluate that curve at any channel. Used when there is no
       usable vendor model / calibration peaks. source-label
       "bootstrap_significant_peaks_quadratic|linear".
    4. Otherwise -> return `fallback_channels` everywhere (last resort).

    The returned callable carries a `.fwhm_source` attribute naming the
    source actually used (F-449 provenance).

    Args:
        spec: a `gamma.spectrum.Spectrum` instance
        fallback_channels: FWHM in channels to return when no stored
            information is usable.
        min_channels: lower clamp on the returned value (the Mariscotti
            kernel becomes degenerate at FWHM < 1).
        bootstrap_from_peaks: F-449 opt-in. When True and stored sources
            are unavailable, attempt the significant-peaks FWHM(E)
            bootstrap before falling back to the constant. Default False
            preserves the legacy constant-fallback behaviour exactly.

    Returns:
        A `Callable[[int], float]` mapping channel index -> FWHM in
        channels, with a `.fwhm_source` attribute.
    """
    sf = getattr(spec, "stored_fwhm_calibration", None)

    # --- Optional: physical scintillator floor ---
    # For scintillator detectors (NaI, CeBr, LaBr, CsI), FWHM(E) ~ sqrt(E)
    # because counting statistics in the photomultiplier dominate. When
    # the vendor FWHM model breaks down at low channels (e.g.
    # SimpleSqrtFwhm with c0+c1*N <= 0 for N below |c0|/c1), naive
    # fallback to min_channels=1 makes Mariscotti miss real low-energy
    # peaks (32 keV NaI K X-ray escape; 46 keV Pb-210 from shielding;
    # 75 keV Pb XRF). To handle this, we compute an alpha coefficient
    # FWHM_keV = alpha*sqrt(E_keV) from a known calibration point and use
    # it as a floor.
    #
    # Selection of the reference point for alpha:
    # 1. Prefer a measured cal_peak with both channel AND fwhm_channels
    #    populated (the most reliable empirical anchor).
    # 2. Else use the geometric midpoint of valid SimpleSqrtFwhm range
    #    (N = 2*|c0|/c1, well into the regime where the model works).
    # 3. Else skip the floor entirely.
    alpha_floor_keV_sqrt = None
    if spec.energy_cal:
        a0_cal = float(spec.energy_cal[0])
        a1_cal = float(spec.energy_cal[1])
        if a1_cal > 0:
            ref_N = None
            ref_FWHM_ch = None
            # Option 1: measured cal peak
            if sf is not None and sf.calibration_peaks:
                for cp in sf.calibration_peaks:
                    if cp.fwhm_channels and cp.fwhm_channels > 0:
                        ref_N = float(cp.channel)
                        ref_FWHM_ch = float(cp.fwhm_channels)
                        break
            # Option 2: midpoint of valid SimpleSqrtFwhm range
            if ref_N is None and sf is not None and sf.model == "SimpleSqrtFwhm" \
                    and len(sf.coefficients) >= 2:
                c0 = float(sf.coefficients[0])
                c1 = float(sf.coefficients[1])
                if c1 > 0:
                    N_low_valid = max(1.0, -c0 / c1) if c0 < 0 else 1.0
                    ref_N = N_low_valid * 2.0  # well into valid regime
                    arg = c0 + c1 * ref_N
                    if arg > 0:
                        ref_FWHM_ch = math.sqrt(arg)
            if ref_N is not None and ref_FWHM_ch is not None:
                # Convert FWHM_ch to FWHM_keV (local gain at ref_N)
                local_gain = a1_cal  # linear approximation; good enough for alpha scaling
                ref_E_keV = a0_cal + a1_cal * ref_N
                ref_FWHM_keV = ref_FWHM_ch * local_gain
                if ref_E_keV > 0 and ref_FWHM_keV > 0:
                    alpha_floor_keV_sqrt = ref_FWHM_keV / math.sqrt(ref_E_keV)

    def _physical_floor_channels(ch: int) -> float:
        """Scintillator physical model FWHM floor at given channel."""
        if alpha_floor_keV_sqrt is None or not spec.energy_cal:
            return float(min_channels)
        a0_cal = float(spec.energy_cal[0])
        a1_cal = float(spec.energy_cal[1])
        E_keV = a0_cal + a1_cal * float(ch)
        if E_keV <= 0 or a1_cal <= 0:
            return float(min_channels)
        FWHM_keV = alpha_floor_keV_sqrt * math.sqrt(E_keV)
        FWHM_ch = FWHM_keV / a1_cal
        return max(float(min_channels), FWHM_ch)

    # --- Source 1: stored SimpleSqrtFwhm coefficients ---
    if sf is not None and sf.model == "SimpleSqrtFwhm" and len(sf.coefficients) >= 2:
        c0 = float(sf.coefficients[0])
        c1 = float(sf.coefficients[1])

        def from_model(ch: int) -> float:
            arg = c0 + c1 * float(ch)
            if arg <= 0:
                # Model breaks down - use physical floor
                return _physical_floor_channels(ch)
            v = math.sqrt(arg)
            # Take max of model and physical floor - physical model is
            # a hard lower bound for any scintillator FWHM.
            floor = _physical_floor_channels(ch)
            return max(float(min_channels), v, floor)

        # Sanity-check against calibration peaks if both are available.
        # If the model produces values that disagree with the recorded
        # measurements by more than 30% at the calibration peaks, we
        # fall through to the interpolation path - the model has been
        # mis-parsed or is on a different convention.
        cal_peaks = list(sf.calibration_peaks or [])
        if cal_peaks:
            disagreement = False
            for cp in cal_peaks:
                if cp.fwhm_channels and cp.fwhm_channels > 0:
                    predicted = from_model(cp.channel)
                    if predicted <= 0:
                        disagreement = True
                        break
                    rel = abs(predicted - cp.fwhm_channels) / cp.fwhm_channels
                    if rel > 0.30:
                        disagreement = True
                        break
            if not disagreement:
                from_model.fwhm_source = "stored_model"
                return from_model

    # --- Source 2: piecewise-linear interpolation on calibration peaks ---
    cal_peaks = []
    if sf is not None:
        cal_peaks = [
            (int(cp.channel), float(cp.fwhm_channels))
            for cp in (sf.calibration_peaks or [])
            if cp.fwhm_channels and cp.fwhm_channels > 0
        ]
    cal_peaks.sort()
    if cal_peaks:
        channels = [p[0] for p in cal_peaks]
        fwhms = [p[1] for p in cal_peaks]

        if len(cal_peaks) == 1:
            const = max(float(min_channels), fwhms[0])

            def single_const(_ch: int) -> float:
                return const

            single_const.fwhm_source = "calibration_peaks"
            return single_const

        def from_interp(ch: int) -> float:
            x = float(ch)
            if x <= channels[0]:
                val = fwhms[0]
            elif x >= channels[-1]:
                val = fwhms[-1]
            else:
                # Locate the segment containing x via linear scan
                # (small number of cal points - fine).
                val = fwhms[-1]
                for i in range(len(channels) - 1):
                    if channels[i] <= x <= channels[i + 1]:
                        x0, x1 = channels[i], channels[i + 1]
                        y0, y1 = fwhms[i], fwhms[i + 1]
                        if x1 == x0:
                            val = y0
                        else:
                            t = (x - x0) / (x1 - x0)
                            val = y0 + t * (y1 - y0)
                        break
            floor = _physical_floor_channels(ch)
            return max(float(min_channels), val, floor)

        from_interp.fwhm_source = "calibration_peaks"
        return from_interp

    # --- Source 3 (F-449): bootstrap FWHM(E) from significant peaks ---
    # Only when explicitly opted in and the stored sources above failed.
    # Seed FWHM for the internal peak search is the constant fallback
    # (scalar => no recursion into this provider).
    if bootstrap_from_peaks:
        seed = max(float(min_channels), float(fallback_channels))
        anchors = _bootstrap_fwhm_anchors(
            spec, seed_fwhm=seed, min_channels=min_channels,
        )
        fit = _fit_fwhm_sq_model(anchors)
        if fit is not None:
            coeffs, source_label = fit

            def from_bootstrap(ch: int) -> float:
                if not spec.energy_cal or len(spec.energy_cal) < 2:
                    return max(float(min_channels), float(fallback_channels))
                E_keV = spec.channel_to_energy(int(ch))
                if E_keV is None:
                    return max(float(min_channels), float(fallback_channels))
                E = max(float(E_keV), 5.0)
                val = sum(c * (E ** i) for i, c in enumerate(coeffs))
                if val <= 0:
                    return _physical_floor_channels(ch)
                fwhm_keV = math.sqrt(val)
                dE_dN = sum(
                    i * a * (float(ch) ** (i - 1))
                    for i, a in enumerate(spec.energy_cal) if i > 0
                )
                if dE_dN <= 0:
                    return max(float(min_channels), float(fallback_channels))
                fwhm_ch = fwhm_keV / dE_dN
                floor = _physical_floor_channels(ch)
                return max(float(min_channels), fwhm_ch, floor)

            from_bootstrap.fwhm_source = source_label
            from_bootstrap.bootstrap_n_anchors = len(anchors)
            from_bootstrap.bootstrap_coeffs = tuple(coeffs)
            return from_bootstrap

    # --- Source 4: caller-supplied constant fallback (last resort) ---
    fb = max(float(min_channels), float(fallback_channels))

    def constant_fb(ch: int) -> float:
        floor = _physical_floor_channels(ch)
        return max(fb, floor)

    constant_fb.fwhm_source = "fallback"
    return constant_fb


__all__ = ["make_fwhm_at_channel_provider"]