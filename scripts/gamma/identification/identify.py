"""
Main identification algorithm per Lsrm Algorithmic Foundations §6.

The Lsrm methodology inverts the naive approach:

   ❌ Naive: for each peak, search the library for matching lines
   ✅ Lsrm:  for each candidate nuclide, search the spectrum for its
            characteristic line; if found, look for all other lines

The Lsrm approach is more robust because:
  - It uses prior knowledge (nuclide line patterns) to disambiguate
    multiple-candidate situations.
  - It explicitly requires the **characteristic line** (lowest MDA)
    before accepting a nuclide as "present".
  - The order of operations gives natural confidence: nuclides found
    via their characteristic line are stronger identifications than
    nuclides matched only by secondary lines.

Three-step algorithm (Lsrm §6):

  Step 1 — Detection:
    For each candidate nuclide in the library:
      a. Compute MDA for each of its library lines (skipping any
         already-matched as fixed nuclides' lines).
      b. Find the characteristic line (lowest MDA).
      c. Check if a found peak lies within δE(E_characteristic) of
         that line's energy.
      d. If yes → nuclide is "detected", proceed to Step 2.
         If no  → reject, move to next candidate.

  Step 2 — Line matching:
    For each detected nuclide:
      a. For each of its library lines, look for a found peak within
         δE of the line's energy.
      b. Build (nuclide, line, peak) match list.

  Step 3 — Promotion / Confidence assessment:
    For each detected nuclide:
      a. Compute Confidence Index from matched lines.
      b. Compute Dose Contribution (DC = % of spectrum's dose
         explained by identified lines — completeness metric per
         Lsrm §14.2).
      c. Flag low-CI identifications for operator review.

This algorithm does NOT compute activities (that's a separate
calculation in Phase 2.1 once all peak areas are known). Identification
focuses on the qualitative question: "which nuclides are present?"

Reference: Lsrm Algorithmic Foundations 2022 §6, §14.2, §14.3.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from gamma.data.nuclide_library import get_nuclide, list_nuclides
from gamma.identification.window import (
    IdentificationWindow, build_identification_window,
)
from gamma.identification.mda import (
    MdaResult, mda_for_peak, characteristic_line_of_nuclide,
)
from gamma.identification.confidence import (
    ConfidenceIndexResult, confidence_index,
)


@dataclass(frozen=True)
class LineMatch:
    """One found peak matched to one library line of a nuclide."""

    nuclide: str
    library_E_keV: float
    library_I_pct: float
    peak_channel: int
    peak_E_keV: float
    peak_sigma: float
    residual_keV: float
    is_characteristic: bool  # True if this is the lowest-MDA line

    # Phase 2.1a: net peak area from integration (optional).
    # When None, downstream code falls back to peak_sigma as proxy.
    # When populated, proportionality and activity calculations use
    # this value directly (more accurate than σ proxy).
    peak_area: Optional[float] = None
    peak_area_uncertainty: Optional[float] = None

    # F-34 / v1.7.12: provenance of `peak_area`. One of
    #   "" / "cowell" / "lsrm_peaks_table" / "deconvolved" / "failed"
    # Set by `gamma.peaks.area.get_peak_area` (returns the source as a
    # third tuple element) and by the multiplet-deconvolution post-pass
    # `gamma.peaks.deconvolve.apply_multiplet_deconvolution`.
    peak_area_source: str = ""

    # BUG-34 Phase 1+2 (NON-BREAKING) → Phase 3 W1 (v1.21.0): explicit
    # successors to the polysemic `peak_sigma` field. Both populated
    # by all writers; readers (compute.py:798, json_report.py:209)
    # prefer `gauss_sigma_keV` for FWHM math, prefer
    # `significance_currie` for significance gates.
    #
    # NOTE on `peak_sigma` (legacy, deprecated semantic):
    # historically polysemic — W1 writers (identify.py characteristic
    # & secondary line) store Currie significance there; W2 writer
    # (staged_pipeline.py phantom-multiplet, line 2130) and W3 phantom
    # anchor writers (deconvolve.py 1060,1247) store FWHM/2.355 in keV.
    # NEW code must read `significance_currie` or `gauss_sigma_keV`
    # explicitly; `peak_sigma` retained only for backward compat.
    #
    # Field semantics (post-v1.21.0):
    #   - significance_currie : Currie L_C significance (dimensionless),
    #     set by identify.py W1 writers (characteristic & secondary).
    #   - gauss_sigma_keV     : Gaussian σ in keV (= FWHM/2.355), set
    #     by identify.py W1 (from `fwhm_at_channel` × keV/ch conv),
    #     staged_pipeline.py W2 (multiplet phantoms), and
    #     deconvolve.py W3 (library-anchor phantoms).
    #
    # See audit/_drafts/BUG-34_peak_sigma_forensics_2026-06-03.md §7
    # Option A, §9 Steps 1-2 and audit/_plans/PLAN_v1_20_to_v1_21.md
    # § P0-2 (BUG-34 W1/W2 writer normalisation).
    significance_currie: Optional[float] = None
    gauss_sigma_keV: Optional[float] = None


@dataclass(frozen=True)
class NuclideIdentification:
    """Identification result for one candidate nuclide."""

    nuclide: str
    detected: bool          # True if characteristic line was found
    reason: str             # Why detected (or not)
    characteristic_line_keV: float
    matched_lines: tuple    # tuple of LineMatch
    confidence: Optional[ConfidenceIndexResult] = None

    def confidence_level(self) -> str:
        if not self.detected or self.confidence is None:
            return "n/a"
        return self.confidence.confidence_level()

    def __repr__(self) -> str:
        if not self.detected:
            return f"NuclideIdentification({self.nuclide}: NOT DETECTED — {self.reason})"
        ci_str = f", CI={self.confidence.CI:.2f}" if self.confidence else ""
        return (f"NuclideIdentification({self.nuclide}: detected, "
                f"{len(self.matched_lines)} lines{ci_str})")


@dataclass(frozen=True)
class IdentificationResult:
    """Result of running identification across all candidate nuclides."""

    detector_type: str
    window: IdentificationWindow
    candidates_considered: int
    detected_nuclides: tuple   # tuple of NuclideIdentification (detected only)
    rejected_nuclides: tuple   # tuple of NuclideIdentification (not detected)
    unmatched_peaks: tuple     # peaks not assigned to any identified line
    notes: str = ""

    def n_detected(self) -> int:
        return len(self.detected_nuclides)

    def summary(self) -> str:
        lines = [
            f"Identification on {self.detector_type} ({self.window.delta_E0_keV:.1f} keV window @ 661):",
            f"  candidates considered: {self.candidates_considered}",
            f"  detected: {self.n_detected()}",
        ]
        for ni in self.detected_nuclides:
            ci = ni.confidence
            ci_str = f"CI={ci.CI:.2f} ({ci.confidence_level()})" if ci else "no CI"
            lines.append(f"    ✓ {ni.nuclide}: {len(ni.matched_lines)} lines, {ci_str}")
        if self.unmatched_peaks:
            lines.append(f"  unmatched peaks: {len(self.unmatched_peaks)}")
        return "\n".join(lines)


def identify_nuclides(
    *,
    found_peaks: list,
    spec,
    candidate_nuclides: Optional[list] = None,
    window: Optional[IdentificationWindow] = None,
    efficiency_model: Optional[callable] = None,
    compute_peak_areas: bool = True,
    fwhm_at_channel: Optional[callable] = None,
    # F-449 (operator-locked 2026-06-16): when provided, the per-peak
    # Gaussian σ in keV is computed DIRECTLY from this energy-domain
    # FWHM curve (σ = FWHM_keV(E_peak)/2.3548), bypassing the legacy
    # round-trip keV→channels→keV via `fwhm_at_channel` + two-point
    # `channel_to_energy` difference. The round-trip introduced a
    # local-dE/dch artifact under non-linear energy calibration
    # (Th-232 Marinelli, deg-4 e-cal at 2614 keV: local dE/dch=3.032 vs
    # nominal 2.804 → ×1.081 inflation of the singlet FWHM from
    # curve-true 107.84 keV to reported 116.63 keV).
    #
    # Pass `fwhm_keV_at_energy_fn=lambda E: fwhm_keV_at_energy(fwhm_model, E)`
    # to use the canonical FWHM(E) curve. When None, the legacy
    # `fwhm_at_channel`-based path is used (preserved for tests /
    # fixtures that do not have access to an energy-domain FWHM curve).
    fwhm_keV_at_energy_fn: Optional[callable] = None,
    # F-123 / v1.17.6 — per-(nuclide, library_E_keV) window overrides.
    # Ключ (nuclide_str, round(E_keV, 2)) → custom half-window в keV.
    # Используется для расширения окна 238 кэВ (Pb-212) до ±2.5·FWHM
    # при доминантной цепочке Th-232.
    line_window_overrides_keV: Optional[dict] = None,
) -> IdentificationResult:
    """
    Run Lsrm-style identification on a calibrated spectrum.

    Args:
        found_peaks: peaks from mariscotti_search() with channel info
        spec: Spectrum with energy_cal — peaks must map to keV
        candidate_nuclides: list of nuclide names to try. If None, uses
                          all nuclides in the library.
        window: identification window. If None, built from spec.detector_id.
        efficiency_model: optional callable(E_keV) → efficiency for MDA
                        computation. If None, MDA is not computed and
                        characteristic line is the highest-intensity
                        library line (simplification).
        compute_peak_areas: if True (default), integrate each matched
                        peak with Cowell method and populate
                        LineMatch.peak_area. This enables accurate
                        intensity-ratio proportionality checks
                        downstream (Phase 2.1a). Set to False to skip
                        for performance.
        fwhm_at_channel: optional callable(ch) → FWHM in channels.
                        Used for peak area ROI sizing. Falls back to
                        FoundPeak.fwhm_channels if not provided.

    Returns:
        IdentificationResult. LineMatch objects include peak_area
        (counts) and peak_area_uncertainty when compute_peak_areas=True.
    """
    if spec.energy_cal is None or len(spec.energy_cal) < 2:
        raise ValueError("Spectrum must have energy_cal of at least degree 1")

    # Build window from detector type if not provided
    if window is None:
        det_type = getattr(spec, "detector_type", None) or \
                   getattr(spec, "detector_id", "") or "NaI"
        det_tokens = spec.filename_tokens.get("detector", []) if hasattr(spec, "filename_tokens") else []
        if det_tokens:
            det_type = det_tokens[0]
        window = build_identification_window(det_type)

    # Candidate list
    if candidate_nuclides is None:
        candidate_nuclides = list_nuclides()

    # Pre-compute peak energies once
    peak_E_keV = [(p, spec.channel_to_energy(p.channel)) for p in found_peaks]

    # BUG-34 W1 writer normalisation (v1.21.0): pre-compute per-peak
    # Gaussian σ in keV so that LineMatch.gauss_sigma_keV is populated
    # with the correct semantic (FWHM/2.355 in keV) rather than relying
    # on downstream callers to fall back to the polysemic `peak_sigma`
    # field (which W1 stores Currie significance in — wrong unit for
    # readers like compute.py:798 chain-sibling FWHM gate).
    #
    # F-449 (operator-locked 2026-06-16, agent-a-math-2 2026-06-20):
    # PREFERRED path is the direct energy-domain FWHM curve
    # `fwhm_keV_at_energy_fn(E)` — σ_keV = FWHM_keV(E_peak) / 2.354820045.
    # This avoids the legacy keV→ch→keV round-trip artifact under
    # non-linear energy calibration. On Th-232 Marinelli (deg-4 e-cal),
    # the singlet 2614 keV reported FWHM=116.63 keV because local
    # dE/dch=3.032 ≠ nominal gain 2.804 → ×1.081 inflation of the
    # curve-true 107.84 keV. Locking σ directly to the curve removes
    # the artifact.
    #
    # Fallback (legacy round-trip via two-point `channel_to_energy`
    # difference) preserved when no `fwhm_keV_at_energy_fn` is passed
    # (tests, fixtures, callers that have only a channel-domain FWHM
    # provider). Same as before — robust to non-linear e-cal in the
    # sense that gain is not assumed constant, but still affected by
    # local dE/dch when the FWHM is wide and the e-cal is curved.
    #
    # Source: BUG-34 carry-forward, audit/_plans/PLAN_v1_20_to_v1_21.md
    # § P0-2 acceptance #1; F-449 σ-lock from FWHM(E) curve per
    # CLAUDE.md «FWHM любого пика = калибровка FWHM(E)» operator rule.
    peak_gauss_sigma_keV_cache: dict = {}  # channel → σ in keV (or None)
    for p in found_peaks:
        try:
            # F-449 preferred path: direct from energy-domain FWHM curve.
            if fwhm_keV_at_energy_fn is not None:
                _, e_peak_kev = next(
                    (pp for pp in peak_E_keV if pp[0] is p),
                    (None, None),
                )
                if e_peak_kev is None:
                    e_peak_kev = float(spec.channel_to_energy(p.channel))
                fwhm_keV_here = float(fwhm_keV_at_energy_fn(float(e_peak_kev)))
                sigma_keV = (
                    fwhm_keV_here / 2.354820045 if fwhm_keV_here > 0 else None
                )
                peak_gauss_sigma_keV_cache[p.channel] = sigma_keV
                continue
            # Legacy fallback: two-point keV difference around the peak.
            if fwhm_at_channel is not None:
                fwhm_ch = float(fwhm_at_channel(p.channel))
            else:
                fwhm_ch = float(getattr(p, "fwhm_channels", 0.0) or 0.0)
            if fwhm_ch <= 0:
                peak_gauss_sigma_keV_cache[p.channel] = None
                continue
            half_ch = fwhm_ch / 2.0
            ch_lo = float(p.channel) - half_ch
            ch_hi = float(p.channel) + half_ch
            E_lo = float(spec.channel_to_energy(ch_lo))
            E_hi = float(spec.channel_to_energy(ch_hi))
            fwhm_keV_here = abs(E_hi - E_lo)
            sigma_keV = fwhm_keV_here / 2.354820045 if fwhm_keV_here > 0 else None
            peak_gauss_sigma_keV_cache[p.channel] = sigma_keV
        except Exception:
            peak_gauss_sigma_keV_cache[p.channel] = None

    # Pre-compute peak areas if requested. We do this once per unique
    # peak (not per LineMatch) since multiple library lines from
    # different nuclides may match the same spectrum peak.
    #
    # F-31 (v1.7.9): use the `get_peak_area` helper which prefers the
    # Lsrm-software-fitted area from `spec.extras["lsrm_peaks_table"]`
    # when available (gaussian + step baseline, much more accurate on
    # closely-spaced peaks like Co-60 1173/1332), and falls back to
    # Cowell when the spectrum has no such table (e.g. AtomSpectra XML).
    # AUDIT-F4 (2026-06-25): кэш площадей поднят на уровень спектра.
    # Ранее `peak_area_cache = {}` пересоздавался на каждом вызове
    # identify_nuclides; staged_pipeline.py вызывает функцию до 3 раз
    # (stage1/2/3) → одни и те же площади считались повторно.
    # Ключ включает fwhm_channels: при уточнении FWHM-модели между
    # стадиями ключ изменится → произойдёт корректный пересчёт, без
    # риска stale-значения. Кэш хранится атрибутом spec, НЕ через
    # spec.extras — иначе cost_estimator._step2_environment посчитал
    # бы его как «поле-расширение» и сдвинул tokens-bookkeeping.
    peak_area_cache: dict = {}  # channel → (area, uncertainty, source)
    spec_area_cache = getattr(spec, "_audit_f4_peak_area_cache", None)
    if spec_area_cache is None:
        spec_area_cache = {}
        try:
            object.__setattr__(spec, "_audit_f4_peak_area_cache", spec_area_cache)
        except (AttributeError, TypeError):
            pass  # spec frozen / __slots__ → degrade to call-local cache
    if compute_peak_areas and hasattr(spec, "counts"):
        from gamma.peaks.area import get_peak_area
        for p in found_peaks:
            try:
                if fwhm_at_channel is not None:
                    fwhm_ch = float(fwhm_at_channel(p.channel))
                else:
                    fwhm_ch = float(getattr(p, "fwhm_channels", 5.0))
                if fwhm_ch <= 0:
                    fwhm_ch = 5.0
                key = (int(p.channel), round(fwhm_ch, 6))
                if key in spec_area_cache:
                    area, dA, src = spec_area_cache[key]
                else:
                    area, dA, src = get_peak_area(
                        spec,
                        peak_channel=int(p.channel),
                        fwhm_channels=fwhm_ch,
                    )
                    spec_area_cache[key] = (area, dA, src)
                peak_area_cache[p.channel] = (area, dA, src)
            except Exception:
                peak_area_cache[p.channel] = (None, None, "failed")

    detected_list = []
    rejected_list = []
    matched_peaks_set = set()

    for nuc_name in candidate_nuclides:
        nuc = get_nuclide(nuc_name)
        if not nuc:
            continue

        lines = nuc.get("lines", [])
        if not lines:
            rejected_list.append(NuclideIdentification(
                nuclide=nuc_name, detected=False,
                reason="No library lines",
                characteristic_line_keV=0.0,
                matched_lines=(),
            ))
            continue

        # Step 1: find the characteristic line.
        # Simplification when no efficiency model: characteristic line
        # = highest-intensity line. With efficiency model (future):
        # characteristic line = lowest-MDA line, requires fitting a
        # baseline under each line first.
        #
        # F-459 (BUG-Y, 2026-06-23): cascade char-line search.
        # For nuclides whose primary characteristic line (highest I%) is
        # absent from the peak list due to known physical suppression (not
        # because the nuclide is absent), fall back to the next-highest-
        # intensity line. Currently only Eu-152 qualifies: its 121.78 keV
        # primary (I=28.53%) is buried under Am-241 Compton continuum in
        # mixed calibration sources (AmTiCsEu Marinelli), while its second
        # line 344.28 keV (I=26.59%) is clearly visible.
        #
        # Design intent: this is NOT a general "try multiple char lines"
        # mechanism. Stage 1/2 nuclides whose primary line is absent are
        # correctly rejected (nuclide not in sample). Stage 3 nuclides may
        # appear in known source mixtures where physics suppresses a line.
        # The whitelist prevents Co-60 from matching Sc-44's 1157 keV peak
        # via its cascade-secondary 1173 keV line, Bi-212 from spuriously
        # detecting via 1620.50 keV (I=1.47%), etc.
        _CASCADE_WHITELIST: frozenset = frozenset({
            "Eu-152",   # 121.78 keV (I=28.53%) hidden under Am-241 Compton;
                        # cascade to 344.28 keV (I=26.59%) in AmTiCsEu mix
        })
        _CHAR_CANDIDATES = 3  # try at most this many lines as char line
        _MIN_CASCADE_I_PCT = 20.0  # cascade fallback line must have ≥ this I%

        lines_sorted_by_I = sorted(
            lines, key=lambda L: -(L[1] if len(L) > 1 else 0)
        )
        characteristic_match = None
        char_E = float(lines_sorted_by_I[0][0])   # primary char line energy
        char_I = float(lines_sorted_by_I[0][1]) if len(lines_sorted_by_I[0]) > 1 else 0.0
        best_residual = float("inf")
        _tried_char: list = []  # for rejection reason

        for _cand_line in lines_sorted_by_I[:_CHAR_CANDIDATES]:
            _c_E = float(_cand_line[0])
            _c_I = float(_cand_line[1]) if len(_cand_line) > 1 else 0.0
            if _cand_line is not lines_sorted_by_I[0]:
                # Cascade fallback: only for whitelisted nuclides with
                # sufficiently intense secondary lines.
                if nuc_name not in _CASCADE_WHITELIST:
                    break
                if _c_I < _MIN_CASCADE_I_PCT:
                    break
            _c_window = window.window_keV(_c_E)
            if line_window_overrides_keV:
                _key_ov = (nuc_name, round(_c_E, 2))
                if _key_ov in line_window_overrides_keV:
                    _c_window = float(line_window_overrides_keV[_key_ov])
            _tried_char.append(f"{_c_E:.1f}±{_c_window:.1f}")
            _best_here = float("inf")
            _match_here = None
            for peak, E_keV in peak_E_keV:
                residual = abs(E_keV - _c_E)
                if residual <= _c_window and residual < _best_here:
                    _match_here = (peak, E_keV)
                    _best_here = residual
            if _match_here is not None:
                characteristic_match = _match_here
                best_residual = _best_here
                char_E = _c_E
                char_I = _c_I
                break  # use first (highest-I) char line that has a peak

        # Step 1b: apply F-123 override to the chosen char_E.
        # (Override already applied in cascade above for each candidate;
        # the chosen char_E already used the correct window.)

        if characteristic_match is None:
            rejected_list.append(NuclideIdentification(
                nuclide=nuc_name, detected=False,
                reason=f"Characteristic line not found within window "
                       f"(tried top-{min(_CHAR_CANDIDATES, len(lines_sorted_by_I))} "
                       f"by intensity: {', '.join(_tried_char)} keV)",
                characteristic_line_keV=char_E,
                matched_lines=(),
            ))
            continue

        # Step 2: now find all other lines for this nuclide.
        matches = []
        char_area, char_area_unc, char_src = peak_area_cache.get(
            characteristic_match[0].channel, (None, None, "")
        )
        matches.append(LineMatch(
            nuclide=nuc_name,
            library_E_keV=char_E,
            library_I_pct=char_I,
            peak_channel=characteristic_match[0].channel,
            peak_E_keV=characteristic_match[1],
            peak_sigma=characteristic_match[0].significance,
            residual_keV=best_residual,
            is_characteristic=True,
            peak_area=char_area,
            peak_area_uncertainty=char_area_unc,
            peak_area_source=char_src,
            # BUG-34 Phase 1+2: explicit Currie significance alias
            significance_currie=characteristic_match[0].significance,
            # BUG-34 W1 writer normalisation (v1.21.0): populate
            # Gaussian σ in keV from peak FWHM so downstream readers
            # (compute.py:798 chain-sibling gate, json_report.py:209)
            # do not fall back to legacy `peak_sigma` which here holds
            # Currie significance (dimensionless, NOT keV).
            gauss_sigma_keV=peak_gauss_sigma_keV_cache.get(
                characteristic_match[0].channel
            ),
        ))
        matched_peaks_set.add(characteristic_match[0].channel)

        for line in lines:
            E = float(line[0])
            I = float(line[1]) if len(line) > 1 else 0.0
            if abs(E - char_E) < 0.01:
                continue
            line_window = window.window_keV(E)
            # F-123: override for specific (nuclide, library_E)
            if line_window_overrides_keV:
                key_ov = (nuc_name, round(E, 2))
                if key_ov in line_window_overrides_keV:
                    line_window = float(line_window_overrides_keV[key_ov])
            best_resid = float("inf")
            best_peak = None
            for peak, E_keV in peak_E_keV:
                residual = abs(E_keV - E)
                if residual <= line_window and residual < best_resid:
                    best_peak = (peak, E_keV)
                    best_resid = residual
            if best_peak is not None:
                peak_area, peak_area_unc, peak_src = peak_area_cache.get(
                    best_peak[0].channel, (None, None, "")
                )
                matches.append(LineMatch(
                    nuclide=nuc_name,
                    library_E_keV=E,
                    library_I_pct=I,
                    peak_channel=best_peak[0].channel,
                    peak_E_keV=best_peak[1],
                    peak_sigma=best_peak[0].significance,
                    residual_keV=best_resid,
                    is_characteristic=False,
                    peak_area=peak_area,
                    peak_area_uncertainty=peak_area_unc,
                    peak_area_source=peak_src,
                    # BUG-34 Phase 1+2: explicit Currie significance alias
                    significance_currie=best_peak[0].significance,
                    # BUG-34 W1 writer normalisation (v1.21.0): see
                    # characteristic-line LineMatch above for rationale.
                    gauss_sigma_keV=peak_gauss_sigma_keV_cache.get(
                        best_peak[0].channel
                    ),
                ))
                matched_peaks_set.add(best_peak[0].channel)

        # Step 3: compute Confidence Index
        matched_lines_for_CI = [
            {"E_keV": m.library_E_keV, "I_pct": m.library_I_pct}
            for m in matches
        ]
        ci_result = confidence_index(
            nuc_name, matched_lines_for_CI, window.window_keV,
        )

        detected_list.append(NuclideIdentification(
            nuclide=nuc_name, detected=True,
            reason=f"Characteristic line {char_E:.1f} keV matched at "
                   f"ch={characteristic_match[0].channel}, "
                   f"E={characteristic_match[1]:.2f}, "
                   f"Δ={best_residual:.2f} keV",
            characteristic_line_keV=char_E,
            matched_lines=tuple(matches),
            confidence=ci_result,
        ))

    # Find unmatched peaks
    unmatched = tuple(
        p for p, _ in peak_E_keV
        if p.channel not in matched_peaks_set
    )

    return IdentificationResult(
        detector_type=window.detector_type,
        window=window,
        candidates_considered=len(candidate_nuclides),
        detected_nuclides=tuple(detected_list),
        rejected_nuclides=tuple(rejected_list),
        unmatched_peaks=unmatched,
    )


__all__ = [
    "LineMatch", "NuclideIdentification", "IdentificationResult",
    "identify_nuclides",
]
