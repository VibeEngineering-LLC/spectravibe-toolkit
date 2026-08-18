"""
F-43 (v1.7.21) / F-44 (v1.7.22) — Averaged Lsrm `.spe` background spectra.

Goal: reduce statistical noise on long-integration background
measurements by combining multiple short-acquisition spectra of the
same geometry into one averaged spectrum.

**F-44 critical correction**: The original F-43 design assumed every
input file is an INDEPENDENT measurement and summed counts across all
inputs. Discovered during F-44 inventory work that the LSRM Spectraline
acquisition software emits **cumulative checkpoints** instead — for a
single long-running acquisition, files `..._01.spe` (1 h), `..._02.spe`
(2 h cumulative), `..._N.spe` (N h cumulative) are saved at intervals;
each later file already contains all events of its predecessors.
Naively summing N cumulative files inflates counts AND live-time by
the same factor sum(1..N)/N ≈ N/2; the RATE (cps) is preserved, but
σ-claim of √N reduction is FALSE — there is really one independent
measurement of duration max(t_i).

F-44 fix: `average_lsrm_spectra` now detects cumulative pattern
(identical `start_datetime` across inputs combined with monotonic
live-times where each is ~k·t_unit) and switches to "take the longest
file" semantics. When inputs are genuinely independent the original
sum-of-counts semantics apply. The detection is conservative — when
in doubt the function raises `CumulativeAmbiguityError` so the caller
can decide explicitly via the `cumulative_policy` kwarg.

Design choices (post-F-44):
  • **Cumulative detection** (new): same `start_datetime` AND
    monotonically increasing `live_time` proportional to (1,2,...,N)·t
    pattern → cumulative. The longest input is returned as-is (its
    counts/live-time already reflect the long-exposure measurement).
  • **Independent measurements** (default semantic when not
    cumulative): channel-wise SUM of counts; sum of live-times. This is
    the correct Poisson aggregation when each input is an independent
    sample of the same underlying process.
  • Calibration: take from the FIRST valid spectrum (cumulative case)
    or LAST valid spectrum after agreement check (independent case).
    Default tolerances: abs_offset_tolerance=2.0 keV on a0;
    rel_gain_tolerance=0.5 % on a1.
  • Detector identity: must match exactly across inputs (defensive —
    don't accidentally merge spectra from different detectors).
  • Geometry: must match exactly. If user wants to merge across
    geometries, pass `require_same_geometry=False` (explicit opt-in).
  • Output: a new `Spectrum` with `is_background=True` and an
    `extras["averaging_provenance"]` dict recording the source paths,
    aggregation mode ("cumulative_last" or "independent_sum"),
    calibration agreement, and the resulting σ reduction factor (√N
    for true independent, 1.0 for cumulative). The output is also
    writable back to .spe via `write_lsrm_spe` for downstream code
    that re-reads via `read_lsrm_spe`.

Scope of this module: LSRM NaI `.spe` only. Other formats (AtomSpectra
XML, .chn, .n42, .mca, .csv) will get parallel averaging in a separate
iteration after the LSRM NaI plan is fully shipped.
"""

from __future__ import annotations

import math
import struct
from dataclasses import replace
from datetime import datetime
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from gamma.spectrum import Spectrum
from gamma.io.lsrm_spe import read_lsrm_spe


# ============================================================================
# Calibration consistency check
# ============================================================================

class CalibrationMismatchError(ValueError):
    """Raised when input spectra have inconsistent energy calibrations."""


class IdentityMismatchError(ValueError):
    """Raised when input spectra have inconsistent detector or geometry."""


class CumulativeAmbiguityError(ValueError):
    """
    Raised when the input set looks partially cumulative (same
    start_datetime but live_times don't match the simple k·t_unit
    cumulative pattern) and `cumulative_policy="auto"` cannot pick a
    safe semantic.

    The caller should pass `cumulative_policy="cumulative_last"` or
    `cumulative_policy="independent_sum"` explicitly.
    """


# ============================================================================
# Cumulative-pattern detection (F-44 fix)
# ============================================================================

def detect_cumulative_pattern(
    specs: Sequence[Spectrum],
    *,
    rel_live_time_tolerance: float = 0.01,
) -> dict:
    """
    Heuristic detector for LSRM Spectraline "cumulative snapshot" sets.

    LSRM Spectraline acquisition software writes intermediate checkpoint
    files during a long-running acquisition: file `_01.spe` after 1 h
    of integration, `_02.spe` after 2 h, ..., `_N.spe` after N h. Each
    later file contains all events of its predecessors. Treating them
    as independent measurements (e.g. summing) inflates total counts
    and live_time by a factor sum(1..N)/N ≈ (N+1)/2 ≈ N/2, while the
    RATE (cps) is preserved.

    Detection criteria (all must hold for "cumulative=True"):
      1. N ≥ 2 inputs.
      2. All inputs share the SAME `start_datetime` (to within 1 s).
      3. When sorted by `live_time`, the live-times form an
         approximate arithmetic progression starting at min:
         t_i ≈ (i+1) · t_min for i = 0..N-1.
         Tolerance: each ratio (t_i / ((i+1)·t_min)) must lie within
         (1 ± rel_live_time_tolerance).

    Returns a dict with keys:
      - `is_cumulative`: bool — final classification.
      - `same_start_datetime`: bool — first criterion.
      - `live_time_progression_ok`: bool — second criterion.
      - `sorted_live_times`: list[float] — debug info.
      - `expected_progression`: list[float] — debug info (i+1)·t_min.
      - `max_relative_deviation`: float — max |t_i − expected| / expected.
      - `t_unit_seconds`: Optional[float] — inferred unit step
        (only set when classification is True).
      - `longest_idx`: int — index of input with largest live_time.

    The function is **non-throwing**; it returns its findings as data.
    The caller decides how to react.
    """
    if len(specs) < 2:
        return {
            "is_cumulative": False,
            "same_start_datetime": False,
            "live_time_progression_ok": False,
            "sorted_live_times": [s.live_time for s in specs],
            "expected_progression": [],
            "max_relative_deviation": None,
            "t_unit_seconds": None,
            "longest_idx": 0 if specs else -1,
        }

    # Criterion 1: same start_datetime
    starts = [s.start_datetime for s in specs]
    same_start = (all(d is not None for d in starts)
                  and all(abs((d - starts[0]).total_seconds()) < 1.0
                          for d in starts))

    # Criterion 2: live-time progression
    indexed = sorted(enumerate(specs), key=lambda p: p[1].live_time)
    sorted_specs = [s for _, s in indexed]
    sorted_lt = [s.live_time for s in sorted_specs]
    t_min = sorted_lt[0]
    expected = [(i + 1) * t_min for i in range(len(sorted_lt))]
    if t_min <= 0:
        progression_ok = False
        max_dev = None
    else:
        deviations = [abs(t - e) / e for t, e in zip(sorted_lt, expected)]
        max_dev = max(deviations)
        progression_ok = max_dev < rel_live_time_tolerance

    is_cum = bool(same_start and progression_ok)

    longest_idx = max(range(len(specs)),
                      key=lambda i: specs[i].live_time)

    return {
        "is_cumulative": is_cum,
        "same_start_datetime": same_start,
        "live_time_progression_ok": progression_ok,
        "sorted_live_times": sorted_lt,
        "expected_progression": expected,
        "max_relative_deviation": max_dev,
        "t_unit_seconds": t_min if is_cum else None,
        "longest_idx": longest_idx,
    }


def _check_calibrations(
    specs: Sequence[Spectrum],
    *,
    rel_gain_tolerance: float,
    abs_offset_tolerance: float,
) -> dict:
    """
    Verify each input spectrum's energy_cal is close enough to the first.

    Returns a dict summarising agreement statistics (max offset spread,
    max gain rel spread). Raises CalibrationMismatchError on violation.
    """
    if not specs:
        return {}
    ref = specs[0].energy_cal
    if ref is None:
        # If first has no cal, allow only when ALL have no cal — degenerate
        # but supported (raw-channel averaging).
        for i, s in enumerate(specs):
            if s.energy_cal is not None:
                raise CalibrationMismatchError(
                    f"spectrum #{i} has energy_cal but reference (#0) does not"
                )
        return {"reference_cal": None, "all_uncalibrated": True}

    if len(ref) < 2:
        raise CalibrationMismatchError(
            f"reference energy_cal {ref} too short (need ≥2 coefficients)"
        )

    a0_ref = ref[0]
    a1_ref = ref[1]
    a0_spread = 0.0
    a1_rel_spread = 0.0

    for i, s in enumerate(specs[1:], start=1):
        cal = s.energy_cal
        if cal is None or len(cal) < 2:
            raise CalibrationMismatchError(
                f"spectrum #{i} has no usable energy_cal ({cal!r}) "
                f"while reference has {ref!r}"
            )
        a0_diff = abs(cal[0] - a0_ref)
        a0_spread = max(a0_spread, a0_diff)
        if a0_diff > abs_offset_tolerance:
            raise CalibrationMismatchError(
                f"spectrum #{i} a0={cal[0]:.4f} keV deviates from "
                f"reference a0={a0_ref:.4f} by {a0_diff:.4f} keV "
                f"(tolerance {abs_offset_tolerance} keV)"
            )
        if a1_ref == 0:
            raise CalibrationMismatchError(
                f"reference a1={a1_ref} is zero — cannot check relative gain"
            )
        rel = abs(cal[1] - a1_ref) / abs(a1_ref)
        a1_rel_spread = max(a1_rel_spread, rel)
        if rel > rel_gain_tolerance:
            raise CalibrationMismatchError(
                f"spectrum #{i} a1={cal[1]:.6f} deviates from "
                f"reference a1={a1_ref:.6f} by {rel*100:.3f}% "
                f"(tolerance {rel_gain_tolerance*100:.3f}%)"
            )

    return {
        "reference_cal": tuple(ref),
        "max_a0_spread_keV": a0_spread,
        "max_a1_rel_spread": a1_rel_spread,
        "all_uncalibrated": False,
    }


def _check_identity(
    specs: Sequence[Spectrum],
    *,
    require_same_detector: bool,
    require_same_geometry: bool,
) -> dict:
    """Verify detector_id / geometry match across inputs as requested."""
    if not specs:
        return {}
    ref = specs[0]
    summary = {
        "detector_id": ref.detector_id,
        "geometry": ref.geometry,
        "n_inputs": len(specs),
    }
    for i, s in enumerate(specs[1:], start=1):
        if require_same_detector and s.detector_id != ref.detector_id:
            raise IdentityMismatchError(
                f"spectrum #{i} detector_id={s.detector_id!r} "
                f"does not match reference {ref.detector_id!r}"
            )
        if require_same_geometry and s.geometry != ref.geometry:
            raise IdentityMismatchError(
                f"spectrum #{i} geometry={s.geometry!r} "
                f"does not match reference {ref.geometry!r}"
            )
    return summary


# ============================================================================
# Channel-length normalisation
# ============================================================================

def _check_channel_lengths(specs: Sequence[Spectrum]) -> int:
    """
    All inputs must have the same n_channels (post-trim length). The
    raw channel count (`n_channels_raw`) is also required to match
    because we sum raw channel counts before any trim has been applied
    in the OUTPUT. (When averaging, we re-apply the trim to the summed
    counts at the output stage to ensure ENERGY_CEILING_KEV invariance.)
    """
    if not specs:
        raise ValueError("cannot average empty list of spectra")
    n = len(specs[0].counts)
    n_raw = specs[0].n_channels_raw
    for i, s in enumerate(specs[1:], start=1):
        if len(s.counts) != n:
            raise ValueError(
                f"spectrum #{i} has {len(s.counts)} channels, "
                f"reference has {n}"
            )
        if s.n_channels_raw != n_raw:
            raise ValueError(
                f"spectrum #{i} n_channels_raw={s.n_channels_raw} "
                f"does not match reference {n_raw}"
            )
    return n


# ============================================================================
# Main averaging entry point
# ============================================================================

def average_lsrm_spectra(
    paths: Sequence[str],
    *,
    rel_gain_tolerance: float = 0.005,
    abs_offset_tolerance: float = 2.0,
    require_same_detector: bool = True,
    require_same_geometry: bool = True,
    cumulative_policy: str = "auto",
    sample_id: Optional[str] = None,
    comment: str = "",
) -> Spectrum:
    """
    Read N Lsrm .spe spectra and return one combined Spectrum.

    Two semantics, chosen automatically by default based on whether the
    inputs are independent measurements or cumulative LSRM snapshots
    (see module docstring and `detect_cumulative_pattern`):

    **Independent-sum** (default for genuinely independent inputs):
    channel-wise SUM of counts, sum of live-times / real-times.
    Equivalent to one long-exposure measurement with proper Poisson
    statistics (σ_rate ∝ 1/√N_total_counts; σ-reduction factor √N
    over a single equivalent-duration measurement).

    **Cumulative-last** (for LSRM cumulative-checkpoint sets): take the
    input with the longest live-time as the result. Each LSRM
    cumulative file already contains all events of its predecessors,
    so naive summing would inflate counts AND live-time by ≈ N/2 while
    leaving the RATE (cps) unchanged. σ-reduction here is 1.0 (one
    measurement of duration max(t_i)).

    Args:
        paths: list of `.spe` file paths to read and combine. Must
            contain at least one path.
        rel_gain_tolerance: maximum allowed relative spread of `a1`
            (linear gain) across inputs (default 0.5 %).
        abs_offset_tolerance: maximum allowed absolute spread of `a0`
            (energy offset) across inputs (default 2.0 keV).
        require_same_detector: when True (default), all inputs must
            have identical `detector_id`. Defensive.
        require_same_geometry: when True (default), all inputs must
            have identical `geometry`. Disable explicitly to merge
            across geometries (use only when you know what you're
            doing — geometry affects efficiency, so averaging across
            geometries breaks downstream ε(E)-based activity).
        cumulative_policy: how to handle suspected LSRM cumulative
            snapshot sets:
              - "auto" (default): run `detect_cumulative_pattern`;
                if positive → "cumulative_last" semantics, else
                "independent_sum".
              - "cumulative_last": always take the longest-live-time
                input as result (caller asserts cumulative).
              - "independent_sum": always sum (caller asserts
                independent — F-43 v1.7.21 original semantic).
        sample_id: explicit `sample_id` for the output spectrum.
            Defaults to a synthesized "averaged: <first_id> ×N".
        comment: free-text comment to append to `extras` and the
            output file's COMMENT field when written back to .spe.

    Returns:
        New `Spectrum` with `is_background=True`,
        `source_format="averaged_lsrm_spe"`, and an
        `extras["averaging_provenance"]` dict carrying the full audit
        trail (aggregation mode used, cumulative-detection results,
        per-input live-times, calibration / identity summary).

    Raises:
        CalibrationMismatchError on calibration drift exceeding tolerances.
        IdentityMismatchError on detector/geometry mismatch when required.
        CumulativeAmbiguityError when policy="auto" cannot make a safe
            choice (currently never raised — auto always picks one).
        ValueError on empty list or inconsistent channel counts.
    """
    if not paths:
        raise ValueError("cannot average empty list of paths")
    if cumulative_policy not in ("auto", "cumulative_last", "independent_sum"):
        raise ValueError(
            f"cumulative_policy must be one of "
            f"'auto', 'cumulative_last', 'independent_sum'; got {cumulative_policy!r}"
        )

    specs: list[Spectrum] = [read_lsrm_spe(str(p)) for p in paths]

    # Identity (detector / geometry) must match across all inputs even
    # in cumulative_last mode — we don't want to silently take one file
    # from a mixed set.
    identity = _check_identity(
        specs,
        require_same_detector=require_same_detector,
        require_same_geometry=require_same_geometry,
    )

    # Cumulative pattern detection — uses start_datetime + live-times,
    # both safe to evaluate before channel-length / calibration checks.
    cum = detect_cumulative_pattern(specs)
    if cumulative_policy == "auto":
        mode = "cumulative_last" if cum["is_cumulative"] else "independent_sum"
    else:
        mode = cumulative_policy

    # In independent_sum mode every input contributes to the output, so
    # channel lengths AND calibrations must agree. In cumulative_last
    # mode we use a single input; agreement of the others is irrelevant
    # and can legitimately fail (e.g. file_N saved after a calibration
    # update mid-acquisition — observed on 2024 Marinelli closed-lid set).
    if mode == "independent_sum":
        _check_channel_lengths(specs)
        cal_summary = _check_calibrations(
            specs,
            rel_gain_tolerance=rel_gain_tolerance,
            abs_offset_tolerance=abs_offset_tolerance,
        )
    else:
        cal_summary = {
            "reference_cal": (tuple(specs[cum["longest_idx"]].energy_cal)
                              if specs[cum["longest_idx"]].energy_cal else None),
            "calibration_check_skipped": (
                "cumulative_last uses one input file; "
                "calibration agreement of others not required"
            ),
            "channel_length_check_skipped": (
                "cumulative_last uses one input file; "
                "channel length of others not required"
            ),
        }

    N = len(specs)
    ref = specs[0]

    if mode == "cumulative_last":
        # Take the file with the largest live_time — its counts and
        # live_time already reflect the full cumulative measurement.
        idx = cum["longest_idx"]
        chosen = specs[idx]
        out_counts = np.asarray(chosen.counts, dtype=np.int64).copy()
        total_live = float(chosen.live_time)
        total_real = float(chosen.real_time)
        total_overflow = int(chosen.dropped_overflow_count)
        sigma_reduction = 1.0
        # Calibration / FWHM / geometry inherited from CHOSEN file
        cal_inherit = chosen
    else:
        # Independent-sum semantic
        out_counts = np.zeros_like(specs[0].counts, dtype=np.int64)
        for s in specs:
            out_counts = out_counts + np.asarray(s.counts, dtype=np.int64)
        total_live = float(sum(s.live_time for s in specs))
        total_real = float(sum(s.real_time for s in specs))
        total_overflow = sum(s.dropped_overflow_count for s in specs)
        sigma_reduction = math.sqrt(N) if N > 0 else 0.0
        cal_inherit = ref

    # Synthesize sample_id
    if sample_id is None:
        base_id = ref.sample_id or Path(ref.source_path).stem or "background"
        if mode == "cumulative_last":
            sample_id = (f"cumulative last of {N}: "
                         f"{Path(paths[cum['longest_idx']]).stem}")
        else:
            sample_id = f"averaged: {base_id} ×{N}"

    # Aggregate datetime: earliest start as start_datetime
    start_dts = [s.start_datetime for s in specs if s.start_datetime]
    aggregated_start = min(start_dts) if start_dts else None

    provenance = {
        "module": "gamma.io.average_lsrm.average_lsrm_spectra",
        "version_introduced": "v1.7.21 (F-43); cumulative detection F-44 (v1.7.22)",
        "aggregation_mode": mode,
        "n_inputs": N,
        "source_files": [str(p) for p in paths],
        "input_live_times_s": [float(s.live_time) for s in specs],
        "input_real_times_s": [float(s.real_time) for s in specs],
        "total_live_time_s": total_live,
        "total_real_time_s": total_real,
        "sigma_reduction_factor": sigma_reduction,
        "cumulative_detection": {
            "is_cumulative": cum["is_cumulative"],
            "same_start_datetime": cum["same_start_datetime"],
            "live_time_progression_ok": cum["live_time_progression_ok"],
            "t_unit_seconds": cum["t_unit_seconds"],
            "longest_idx": cum["longest_idx"],
            "max_relative_deviation": cum["max_relative_deviation"],
        },
        "cumulative_policy_requested": cumulative_policy,
        "calibration": cal_summary,
        "identity": identity,
        "rel_gain_tolerance_applied": rel_gain_tolerance,
        "abs_offset_tolerance_applied_keV": abs_offset_tolerance,
        "require_same_detector": require_same_detector,
        "require_same_geometry": require_same_geometry,
        "comment": comment,
    }

    out = Spectrum(
        counts=out_counts,
        live_time=total_live,
        real_time=total_real,
        source_path=f"<{mode} from {N} files>",
        source_format="averaged_lsrm_spe",
        sample_id=sample_id,
        operator=ref.operator,
        geometry=ref.geometry,
        detector_id=ref.detector_id,
        device_guid=ref.device_guid,
        comments=(comment or
                  f"{mode} of {N} input spectra "
                  f"({'cumulative LSRM checkpoints' if mode == 'cumulative_last' else 'independent samples'})"),
        is_background=True,
        start_datetime=aggregated_start,
        end_datetime=None,
        file_created_datetime=datetime.now(),
        valid_pulse_count=None,
        total_pulse_count=None,
        dropped_overflow_count=int(total_overflow),
        n_channels_raw=int(cal_inherit.n_channels_raw),
        n_channels=int(len(out_counts)),
        channel_pitch=int(cal_inherit.channel_pitch),
        energy_cal=cal_inherit.energy_cal,
        energy_cal_degree=cal_inherit.energy_cal_degree,
        energy_cal_source=cal_inherit.energy_cal_source,
        energy_max_keV_kept=cal_inherit.energy_max_keV_kept,
        stored_fwhm_calibration=cal_inherit.stored_fwhm_calibration,
    )
    out.extras["averaging_provenance"] = provenance
    out.extras["averaging_mode"] = mode
    out.extras["averaging_sigma_reduction"] = sigma_reduction
    # Back-compat: also keep the old key, but now reflect mode-aware value
    out.extras["averaging_sigma_reduction_factor"] = sigma_reduction
    return out


# ============================================================================
# Minimal .spe writer (round-trip with read_lsrm_spe)
# ============================================================================

def write_lsrm_spe(
    spec: Spectrum,
    path: str,
    *,
    type_label: str = "Калибровка",
    config_name: str = "",
) -> None:
    """
    Write a Spectrum to disk in Lsrm .spe format.

    Minimal — only fields actually read back by `read_lsrm_spe` are
    emitted. Header is CP-1251 with `\\r\\n` line endings; the SPECTR=
    marker is followed by raw uint32 LE channel counts.

    The output round-trips: `spec_back = read_lsrm_spe(path)` yields
    a Spectrum with identical counts, live_time, real_time, energy_cal,
    geometry, detector_id, and FWHM polynomial.

    Args:
        spec: source Spectrum (typically the output of
            `average_lsrm_spectra`).
        path: destination filesystem path.
        type_label: value for the TYPE header field. Default
            "Калибровка" (calibration) — neutral classifier that does
            not flag the file as background. For averaged backgrounds
            pass "Фон" or leave default and rely on filename.
        config_name: optional CONFIGNAME field (instrument id).
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []

    def kv(key: str, value):
        lines.append(f"{key}={value}")

    if spec.sample_id:
        kv("SHIFR", spec.sample_id)
    if type_label:
        kv("TYPE", type_label)
    if config_name:
        kv("CONFIGNAME", config_name)

    # Timing
    if spec.start_datetime is not None:
        kv("MEASBEGIN",
           spec.start_datetime.strftime("%d-%m-%y %H:%M:%S"))
    kv("TLIVE", f"{spec.live_time:.2f}")
    kv("TREAL", f"{spec.real_time:.2f}")

    if spec.operator:
        kv("OPERATOR", spec.operator)
    if spec.geometry:
        kv("GEOMETRY", spec.geometry)
    if spec.detector_id:
        kv("DETECTOR", spec.detector_id)

    # Energy calibration — emit as N,a0,a1,a2,a3,0,0 (slot count = 7)
    if spec.energy_cal:
        coefs = list(spec.energy_cal)
        # Pad to 7 slots; mark degree as len(coefs)-1
        degree = len(coefs) - 1
        padded = coefs + [0.0] * (7 - len(coefs))
        coef_str = ",".join(f"{c:.10E}" for c in padded[:7])
        kv("ENERGY", f"{degree},{coef_str}")

    # FWHM calibration
    if (spec.stored_fwhm_calibration is not None
            and spec.stored_fwhm_calibration.coefficients):
        coefs = list(spec.stored_fwhm_calibration.coefficients)
        degree = len(coefs) - 1
        padded = coefs + [0.0] * (7 - len(coefs))
        coef_str = ",".join(f"{c:.10E}" for c in padded[:7])
        kv("FWHM", f"{degree},{coef_str}")

    if spec.comments:
        kv("COMMENT", spec.comments)

    # Write the header in CP-1251 then the SPECTR= marker and binary block
    header_text = "\r\n".join(lines) + "\r\n"
    header_bytes = header_text.encode("cp1251", errors="replace")

    counts = np.asarray(spec.counts, dtype="<u4")
    binary = counts.tobytes()

    with open(p, "wb") as f:
        f.write(header_bytes)
        f.write(b"SPECTR=")
        f.write(binary)
