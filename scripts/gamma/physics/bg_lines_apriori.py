"""
Background-line a-priori injection — LSRM §9.4
("Разметка спектра с выставлением фоновых линий").

Why this exists (vs the existing channel-wise background subtraction)
-------------------------------------------------------------------
Channel-wise background subtraction (`physics/background_subtraction.py`)
requires that the *energy calibrations* of the sample and the background
spectrum agree to within fractions of a resolution element. When they
don't agree — even by a fraction of a keV — subtracting one from the
other introduces a synthetic "dipole" (positive + negative bin) at
every background peak. That dipole then propagates into the fit
residual and the area uncertainty of every nearby sample peak.

For Gamma-1S running with Aspect electronics, calibration drifts of
~0.5–1 keV between successive Marinelli measurements (≥ 12 h apart)
are common — small enough to be invisible by eye but large enough to
corrupt subtraction at intense background lines (²⁰⁸Tl 2614,
⁴⁰K 1461).

LSRM §9.4 recommends a different approach:

  * Keep the background's PEAK LIST (positions + rates), not the
    channel-wise spectrum.
  * In the SAMPLE-spectrum fit, **inject each background peak as an
    extra component with a known a-priori rate** (= net background
    rate, in cps).
  * The peak position is allowed to float within the identification
    window so calibration drift between sample and background is
    absorbed by the fit.
  * The amplitude is constrained to (rate · t_live_sample), with
    uncertainty (rate σ) inherited from the background measurement.

This module provides:

  • `BackgroundPeak` — one a-priori prior on a background peak.
  • `BackgroundLineLibrary` — a collection of priors with helpful
    builders from typical bg spectra (the project's averaged
    backgrounds under `data/averaged_bg/` are the canonical input).
  • `inject_priors_into_roi()` — given an ROI fit definition (peak
    positions + initial widths), augment with background priors so
    the downstream multiplet fit can include them.

This module is **additive**. The existing channel-wise subtraction is
retained as the default fallback. Callers can opt-in to a-priori bg
treatment via the new orchestrator kwarg `bg_strategy="a_priori"`
(default remains "channel_subtract") once wiring is complete.

Reference
---------
LSRM Algorithmic Foundations 2022, §9.4 «Разметка спектра с
выставлением фоновых линий», стр. 9-3.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, List, Optional, Sequence, Tuple


@dataclass(frozen=True)
class BackgroundPeak:
    """
    One background prior — an LSRM-§9.4-style a-priori peak.

    The fit will add this as an extra component whose AMPLITUDE is
    constrained to (rate · t_live) ± σ. The position is floated within
    the identification window; the width is fixed from calibration.
    """
    E_keV: float                  # energy of the bg peak
    rate_cps: float               # net background rate (already
                                  # bg-from-bg-spectrum, in cps)
    sigma_rate_cps: float = 0.0   # 1σ uncertainty of rate_cps
    nuclide_hint: str = ""        # e.g. "K-40", "Tl-208", "Bi-214"
    notes: str = ""

    def expected_amplitude(self, t_live_s: float) -> float:
        """Net counts expected in the sample's live time."""
        return self.rate_cps * t_live_s

    def amplitude_sigma(self, t_live_s: float) -> float:
        return self.sigma_rate_cps * t_live_s


@dataclass(frozen=True)
class BackgroundLineLibrary:
    """A collection of background priors (the project's bg peak list)."""
    peaks: tuple = ()                  # tuple[BackgroundPeak, ...]
    source_label: str = ""             # e.g. averaged-empty-Marinelli-2016
    t_live_source_s: Optional[float] = None
    notes: str = ""

    def peaks_in_range(self, E_low_keV: float, E_high_keV: float) -> List[BackgroundPeak]:
        return [p for p in self.peaks if E_low_keV <= p.E_keV <= E_high_keV]

    def peak_near(self, E_keV: float, window_keV: float) -> Optional[BackgroundPeak]:
        """Return the closest bg peak within window, else None."""
        candidates = [p for p in self.peaks if abs(p.E_keV - E_keV) <= window_keV]
        if not candidates:
            return None
        return min(candidates, key=lambda p: abs(p.E_keV - E_keV))

    def __len__(self) -> int:
        return len(self.peaks)


def build_library_from_peak_list(
    peaks_data: Iterable,
    source_label: str = "",
    t_live_source_s: Optional[float] = None,
) -> BackgroundLineLibrary:
    """
    Build a BackgroundLineLibrary from any iterable of objects with
    `.E_keV`, `.rate_cps`, optional `.sigma_rate_cps` and `.nuclide_hint`
    attributes, OR from dicts with those same keys.
    """
    out: List[BackgroundPeak] = []
    for p in peaks_data:
        if isinstance(p, dict):
            out.append(BackgroundPeak(
                E_keV=float(p["E_keV"]),
                rate_cps=float(p["rate_cps"]),
                sigma_rate_cps=float(p.get("sigma_rate_cps", 0.0)),
                nuclide_hint=str(p.get("nuclide_hint", "")),
                notes=str(p.get("notes", "")),
            ))
        else:
            out.append(BackgroundPeak(
                E_keV=float(p.E_keV),
                rate_cps=float(p.rate_cps),
                sigma_rate_cps=float(getattr(p, "sigma_rate_cps", 0.0)),
                nuclide_hint=str(getattr(p, "nuclide_hint", "")),
                notes=str(getattr(p, "notes", "")),
            ))
    # sort by energy
    out.sort(key=lambda p: p.E_keV)
    return BackgroundLineLibrary(
        peaks=tuple(out),
        source_label=source_label,
        t_live_source_s=t_live_source_s,
        notes=f"{len(out)} bg priors",
    )


@dataclass(frozen=True)
class BackgroundPriorInjection:
    """
    Output of `inject_priors_into_roi()` — the extra components that
    the downstream fit should include.
    """
    extra_peak_positions_keV: tuple = ()    # positions to fit
    expected_amplitudes_counts: tuple = ()  # prior on net counts
    expected_amplitudes_sigma_counts: tuple = ()
    peak_labels: tuple = ()                 # for diagnostics

    def n_extra(self) -> int:
        return len(self.extra_peak_positions_keV)


def inject_priors_into_roi(
    library: BackgroundLineLibrary,
    *,
    E_low_keV: float,
    E_high_keV: float,
    t_live_sample_s: float,
    existing_peak_positions_keV: Sequence[float] = (),
    min_separation_keV: float = 1.0,
) -> BackgroundPriorInjection:
    """
    Pick background priors in the ROI [E_low_keV, E_high_keV] that are
    not already covered by an existing peak in the sample-side fit
    list.

    Parameters
    ----------
    library : BackgroundLineLibrary
        The bg peak list to inject from.
    E_low_keV, E_high_keV : ROI bounds for this fit.
    t_live_sample_s : live time of the sample spectrum (for prior
        amplitude conversion).
    existing_peak_positions_keV : peaks the fit is already aware of
        (sample-side) — bg peaks within `min_separation_keV` of these
        are dropped to avoid degenerate coupling.
    min_separation_keV : merge threshold (default 1 keV — finer than
        a NaI FWHM is unproductive).

    Returns
    -------
    BackgroundPriorInjection
    """
    candidates = library.peaks_in_range(E_low_keV, E_high_keV)
    if not candidates:
        return BackgroundPriorInjection()

    existing = list(existing_peak_positions_keV)

    pos_list: List[float] = []
    amp_list: List[float] = []
    sig_list: List[float] = []
    label_list: List[str] = []

    for p in candidates:
        # skip if too close to an existing peak
        if existing and any(abs(p.E_keV - e) < min_separation_keV for e in existing):
            continue
        pos_list.append(p.E_keV)
        amp_list.append(p.expected_amplitude(t_live_sample_s))
        sig_list.append(p.amplitude_sigma(t_live_sample_s))
        label_list.append(p.nuclide_hint or f"bg@{p.E_keV:.1f}")

    return BackgroundPriorInjection(
        extra_peak_positions_keV=tuple(pos_list),
        expected_amplitudes_counts=tuple(amp_list),
        expected_amplitudes_sigma_counts=tuple(sig_list),
        peak_labels=tuple(label_list),
    )


__all__ = [
    "BackgroundPeak",
    "BackgroundLineLibrary",
    "build_library_from_peak_list",
    "BackgroundPriorInjection",
    "inject_priors_into_roi",
]
