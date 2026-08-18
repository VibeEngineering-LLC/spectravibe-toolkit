"""
Core data structures for spectrum analysis.

Single point of truth: a parsed spectrum, regardless of source format,
is a Spectrum dataclass with these fields. All downstream code consumes
this dataclass without caring about the original file format.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional


# ============================================================================
# Energy ceiling (per skill scope — see SKILL.md and references/01)
# ============================================================================

# Energy axis is truncated at this value at the reader stage. Any channel
# whose calibrated energy exceeds ENERGY_CEILING_KEV is dropped from the
# counts array and not visible to downstream steps. Rationale: most
# laboratory-relevant gamma lines lie below 3 MeV (²⁰⁸Tl 2614.5 is the
# practical upper anchor); higher energies typically host pile-up,
# cosmic-ray secondaries, and overflow markers — not signal of interest
# for the workflows this skill implements.
ENERGY_CEILING_KEV = 3000.0


@dataclass
class FwhmCalPeak:
    """A single peak used by the instrument's FWHM calibration."""
    channel: int
    energy_keV: float
    fwhm_channels: float


@dataclass
class StoredFwhmCalibration:
    """
    FWHM calibration as stored in the file. Format varies by vendor; we
    keep a vendor-neutral representation: the calibration peaks, the
    fitted coefficients (low-to-high order), and shape parameters.
    """
    calibration_peaks: list = field(default_factory=list)  # list[FwhmCalPeak]
    coefficients: tuple = ()
    peak_type: Optional[int] = None
    left_tail: Optional[float] = None
    right_tail: Optional[float] = None
    chi2_per_dof: Optional[float] = None
    model: str = ""  # e.g., "SimpleSqrtFwhm"


@dataclass
class Spectrum:
    """
    Parsed gamma-ray spectrum. Format-independent intermediate.

    Conventions:
      - counts[i] is the count in channel i (channel index 0 is the first
        recorded channel, no offset).
      - Channels above ENERGY_CEILING_KEV are dropped during reading; the
        kept length is exposed via n_channels.
      - Overflow markers (typically the very last channel of MCA output
        with a value orders of magnitude above the local tail) are
        dropped by readers and reported in dropped_overflow_count.
      - energy_cal coefficients are stored low-to-high so that:
            E(N) = sum(a_i * N**i for i, a_i in enumerate(energy_cal))
      - All times are in seconds. start_datetime is timezone-aware where
        available.
    """
    counts: object  # np.ndarray, but kept untyped to avoid hard numpy dep here
    live_time: float
    real_time: float

    # Identity / provenance
    source_path: str = ""
    source_format: str = ""              # "atomspectra_xml", "spe", "csv", ...
    sample_id: str = ""
    operator: str = ""
    geometry: str = ""
    detector_id: str = ""
    device_guid: str = ""
    comments: str = ""
    is_background: bool = False

    # Timing
    start_datetime: Optional[datetime] = None
    end_datetime: Optional[datetime] = None
    file_created_datetime: Optional[datetime] = None

    # Counts/pulses (vendor-reported)
    valid_pulse_count: Optional[int] = None
    total_pulse_count: Optional[int] = None
    dropped_overflow_count: int = 0       # how many channels were trimmed as overflow

    # Sample mass / volume (vendor-stored; LSRM SAMPLEMASS / SAMPLEVOLUME).
    # BUG-1 / 2026-06-02: typed surface so downstream code does not have to
    # dig into `extras["lsrm_sample_mass_kg"]`. Units:
    #   - sample_mass_kg, sample_mass_uncertainty_kg — kilograms
    #     (LSRM SAMPLEMASS stored in grams; reader divides by 1000)
    #   - sample_volume_ml, sample_volume_uncertainty_ml — millilitres
    #     (LSRM SAMPLEVOLUME stored in millilitres; passed through verbatim)
    # All four are None when the reader could not extract the field.
    sample_mass_kg: Optional[float] = None
    sample_mass_uncertainty_kg: Optional[float] = None
    sample_volume_ml: Optional[float] = None
    sample_volume_uncertainty_ml: Optional[float] = None

    # Energy axis
    n_channels_raw: int = 0               # before energy ceiling / overflow trim
    n_channels: int = 0                   # after trim; == len(counts)
    channel_pitch: int = 1
    energy_cal: Optional[tuple] = None    # (a0, a1, a2, ...) low-to-high
    energy_cal_degree: Optional[int] = None
    energy_cal_source: str = ""           # "stored" | "bootstrap" | "manual"
    energy_max_keV_kept: Optional[float] = None  # max energy after trimming

    # FWHM calibration (vendor-stored, if any)
    stored_fwhm_calibration: Optional[StoredFwhmCalibration] = None

    # Linked / embedded background
    # AtomSpectra: BackgroundSpectrumFile carries the filename hint, and
    # the actual background spectrum may be embedded in
    # <BackgroundEnergySpectrum> in the same file. When that is the case,
    # background_embedded holds the parsed sub-spectrum and
    # background_link is set for traceability.
    background_link: Optional[str] = None
    background_embedded: Optional["Spectrum"] = None

    # Filename-derived hints (token parsing — see io/filename_hints.py)
    filename_tokens: dict = field(default_factory=dict)

    # Anything we want to keep but don't have a typed slot for
    extras: dict = field(default_factory=dict)

    # ------------------------------------------------------------------
    # Derived properties
    # ------------------------------------------------------------------

    @property
    def dead_time_fraction(self) -> float:
        if self.real_time <= 0:
            return 0.0
        return max(0.0, 1.0 - self.live_time / self.real_time)

    @property
    def dead_time_pct(self) -> float:
        return 100.0 * self.dead_time_fraction

    def channel_to_energy(self, ch) -> Optional[float]:
        """Apply the energy calibration polynomial. Returns None if no cal."""
        if self.energy_cal is None:
            return None
        if hasattr(ch, "__iter__"):
            return [self.channel_to_energy(c) for c in ch]
        return sum(a * ch ** i for i, a in enumerate(self.energy_cal))

    def energy_to_channel(self, e: float) -> Optional[float]:
        """Inverse of channel_to_energy. Linear closed-form; Newton otherwise."""
        if self.energy_cal is None:
            return None
        if len(self.energy_cal) == 2:
            a0, a1 = self.energy_cal
            if a1 == 0:
                return None
            return (e - a0) / a1
        # Newton iteration
        N = (e - self.energy_cal[0]) / max(self.energy_cal[1], 1e-9)
        for _ in range(20):
            E_N = self.channel_to_energy(N)
            dE_dN = sum(i * a * N ** (i - 1)
                        for i, a in enumerate(self.energy_cal) if i > 0)
            if dE_dN == 0:
                break
            step = (E_N - e) / dE_dN
            N = N - step
            if abs(step) < 1e-6:
                break
        return N

    # ------------------------------------------------------------------
    # Serialization for JSON CLI summary and parser_evals
    # ------------------------------------------------------------------

    def to_summary_dict(self) -> dict:
        """JSON-serializable summary. Excludes large arrays (counts)."""
        d = {
            "source_path": self.source_path,
            "source_format": self.source_format,
            "sample_id": self.sample_id,
            "operator": self.operator,
            "detector_id": self.detector_id,
            "device_guid": self.device_guid,
            "is_background": self.is_background,
            "live_time_s": self.live_time,
            "real_time_s": self.real_time,
            "dead_time_pct": round(self.dead_time_pct, 4),
            "valid_pulse_count": self.valid_pulse_count,
            "total_pulse_count": self.total_pulse_count,
            "dropped_overflow_count": self.dropped_overflow_count,
            "n_channels_raw": self.n_channels_raw,
            "n_channels": self.n_channels,
            "channel_pitch": self.channel_pitch,
            "energy_cal": list(self.energy_cal) if self.energy_cal else None,
            "energy_cal_degree": self.energy_cal_degree,
            "energy_cal_source": self.energy_cal_source,
            "energy_max_keV_kept": self.energy_max_keV_kept,
            "energy_ceiling_keV": ENERGY_CEILING_KEV,
            "start_datetime": (self.start_datetime.isoformat()
                               if self.start_datetime else None),
            "end_datetime": (self.end_datetime.isoformat()
                             if self.end_datetime else None),
            "file_created_datetime": (self.file_created_datetime.isoformat()
                                      if self.file_created_datetime else None),
            "background_link": self.background_link,
            "has_background_embedded": self.background_embedded is not None,
            "filename_tokens": self.filename_tokens,
        }

        if self.stored_fwhm_calibration is not None:
            sf = self.stored_fwhm_calibration
            d["stored_fwhm_calibration"] = {
                "model": sf.model,
                "coefficients": list(sf.coefficients),
                "peak_type": sf.peak_type,
                "left_tail": sf.left_tail,
                "right_tail": sf.right_tail,
                "chi2_per_dof": sf.chi2_per_dof,
                "n_calibration_peaks": len(sf.calibration_peaks),
                "calibration_peaks": [
                    {"channel": p.channel,
                     "energy_keV": p.energy_keV,
                     "fwhm_channels": p.fwhm_channels}
                    for p in sf.calibration_peaks
                ],
            }
        else:
            d["stored_fwhm_calibration"] = None

        # Top-level counts summary, not the array itself
        try:
            import numpy as np
            if self.counts is not None and len(self.counts) > 0:
                d["total_counts"] = int(np.sum(self.counts))
                d["peak_channel_value"] = int(np.max(self.counts))
                d["peak_channel_index"] = int(np.argmax(self.counts))
            else:
                d["total_counts"] = 0
                d["peak_channel_value"] = 0
                d["peak_channel_index"] = None
        except Exception:
            pass

        if self.background_embedded is not None:
            d["background_embedded_summary"] = (
                self.background_embedded.to_summary_dict()
            )

        return d


# ============================================================================
# Post-read trim helper (BUG-9 / v1.18.32, 2026-06-03)
# ============================================================================

def trim_to_working_energy(
    spec: Spectrum,
    max_keV: float = ENERGY_CEILING_KEV,
) -> Spectrum:
    """
    Drop trailing channels whose calibrated energy exceeds `max_keV`.

    Introduced in v1.18.32 as the explicit, auditable counterpart to the
    reader-stage ``apply_energy_ceiling`` knob (whose default was flipped
    to ``False`` in the same release — see BUG-9). Use this when a
    pipeline truly needs the 3 MeV scope cut and wants the action to be
    visible in the call site rather than hidden in the reader.

    Returns the same ``spec`` object, mutated in place:
      * ``spec.counts`` is replaced with the trimmed array;
      * ``spec.n_channels`` is updated to the new length;
      * ``spec.energy_max_keV_kept`` is set to the largest kept energy;
      * ``spec.n_channels_raw`` is left untouched so the original channel
        count remains visible for audit.

    If the spectrum has no energy calibration the call is a no-op.

    Args:
        spec: parsed :class:`Spectrum`.
        max_keV: upper energy ceiling. Defaults to
            :data:`ENERGY_CEILING_KEV` (3000 keV by project scope).

    Returns:
        The (mutated) ``spec`` for fluent chaining.
    """
    if spec is None or spec.energy_cal is None or spec.counts is None:
        return spec
    try:
        import numpy as _np
        n = int(len(spec.counts))
    except Exception:
        return spec
    if n == 0:
        return spec
    coeffs = tuple(spec.energy_cal)
    keep = n
    for ch in range(n - 1, -1, -1):
        e = sum(a * (ch ** i) for i, a in enumerate(coeffs))
        if e <= max_keV:
            keep = ch + 1
            break
        keep = ch
    if keep < n:
        spec.counts = _np.asarray(spec.counts)[:keep].copy()
        spec.n_channels = int(keep)
    if keep > 0:
        e_max = sum(a * ((keep - 1) ** i) for i, a in enumerate(coeffs))
        spec.energy_max_keV_kept = float(e_max)
    if spec.background_embedded is not None:
        trim_to_working_energy(spec.background_embedded, max_keV=max_keV)
    return spec
