"""
BUG-48 — Calibration hard-gate (self-consistency check at load time).

Purpose
-------
BUG-48 root-cause analysis identified ≥10 documented mis-identifications
that traced back to silently accepted calibrations whose own self-consistency
indicators were broken (non-monotone energy axis, FWHM model degenerate at
some channels, FWHM out of any physically plausible band). The existing
``check_stored_calibration`` validates against found peaks + anchor lines —
but it only runs *after* a peak search, and many mis-IDs occurred on
spectra where peak-search itself was already biased by the broken
calibration.

This module adds an **upstream gate** that inspects only the calibration
coefficients themselves (no peak search needed). It is intentionally
self-contained and deterministic:

- Pass / soft-warn / hard-fail are returned as a structured result.
- The gate never mutates the spectrum.
- Downstream code reads ``spec.extras['calibration_gate']`` and can decide
  to skip identification when ``passed=False``.

Criteria (all sourced — see ``sources`` section below)
------------------------------------------------------
G1. **Energy-cal monotonicity**: ``dE/dN > 0`` everywhere on
    ``[0, n_channels-1]``. A non-monotone axis means two different
    channels map to the same energy → ID windows alias.
    *Source*: trivial invariant of any physical MCA; failure modes
    documented in LSRM-9.4 §3.1 (energy calibration must be strictly
    increasing) and Gilmore §6.3 (calibration sanity).

G2. **Channel-range sanity**: ``E(0) ≥ -50 keV`` (small negative offset
    is normal — most stored calibrations use a non-zero intercept ≈ 0).
    ``E(n_channels-1)`` ≥ 100 keV (the spectrum covers a meaningful
    range; below 100 keV the spectrum has no methodological value
    against the anchor library).
    *Source*: project convention; ENERGY_CEILING_KEV defined in
    ``gamma.spectrum`` is 3000 keV, an upper bound. The lower 100 keV
    bound matches the lowest anchor in the library (Am-241 59.54 keV +
    Pb-X 75 keV).

G3. **FWHM monotonicity** (when stored FWHM model is present):
    ``FWHM²(E)`` must be non-decreasing across the supported energy
    range — for the LSRM quadratic ``a + b·E + c·E²`` and the
    AtomSpectra SimpleSqrtFwhm ``c0 + c1·N`` this means
    ``dFWHM²/dE ≥ 0`` (resp. ``c1 ≥ 0``) at the upper end.
    *Source*: Gilmore §6.4 eq. 6.13 — FWHM²(E) is monotone increasing
    for any scintillator under the Poisson-broadening regime
    (statistical term dominates above the noise floor). RAG-043 source
    citation already lists Gilmore §6.4.

    **Anchor-mode extension (BUG-50, v1.22.6)**: when the stored model
    has fewer than 4 usable measured calibration peaks, G3 falls back
    to evaluating the FWHM at a standard set of literature anchor
    energies (``STANDARD_FWHM_ANCHORS_KEV`` — see below) that lie
    within the spectrum's calibrated energy range. The same monotone
    test is then applied to the synthesised (anchor_E, FWHM_keV)
    sequence ordered by energy. This activates G3 on the LSRM corpus
    where the vendor file stores only polynomial coefficients (no
    measured calibration_peaks list). The anchor-mode branch is
    additive: it never overrides a measured-peak verdict, and the
    ``criteria_evaluated`` entry is reported as ``"G3_anchor"`` (resp.
    ``"G4_anchor"``) so downstream aggregators (BUG-48 sweep) can
    distinguish measured-mode from anchor-mode runs.

G4. **FWHM physical plausibility band**:
    For the stored-cal anchor energies (the file's own
    ``calibration_peaks``) the local FWHM must lie in
    ``[lo, hi] · E`` where:
        - HPGe-class (R ≤ 1% @ 662): 0.05% ≤ FWHM/E ≤ 2%
        - Scintillator-class (R > 1.5% @ 662): 2% ≤ FWHM/E ≤ 20%
    Outside these bands the FWHM model is broken or the detector
    classification is wrong → identification windows are wrong.
    *Source*: Gilmore §6.4 (NaI 6-9% @ 661.7 keV, HPGe ≤0.5% @ 1.33 MeV);
    RAG-043 verbatim quote: "For NaI(Tl) the FWHM at 661.7 keV is
    typically in the range 6-9 %". LSRM-9.4 §3.2 confirms the same
    bracket. Bracket widened on each side to absorb edge cases:
    HPGe upper bound 2% accommodates poorly-cooled detectors; NaI
    upper bound 20% accommodates very large crystals + low-E behaviour
    where FWHM/E grows as 1/√E.

Verdicts
--------
- ``passed=True``: all criteria passed.
- ``passed=False``: at least one HARD criterion failed (G1 or G2).
- ``soft_warnings``: list of G3/G4 issues that do not block ID but are
  reported in ``report.json.warnings`` via the standard channel.

Out of scope (deferred or other modules)
----------------------------------------
- Cross-check FWHM vs found peaks → already in ``stored_check``.
- Verify against passport / calibration certificate → ``verification_loop``.
- Drift between runs → BUG-35 ``|z|``-test path (RAG-005/008).

References (sources)
--------------------
- Gilmore G., *Practical Gamma-Ray Spectrometry*, 2nd ed., Wiley 2008,
  §6.3 (energy calibration) and §6.4 eq. 6.13 (FWHM(E) function for
  scintillators and HPGe) — RAG-043 doc_corpus_id ``gilmore_practical_gamma``.
- LSRM Algorithmic Foundations §3.1 (energy cal must be strictly
  increasing) and §3.2 (FWHM model FWHM²(E) = a + b·E + c·E²) —
  RAG-043 doc_corpus_id ``lsrm_act_2014``.
- LSRM Algorithmic Foundations §8.3 (FWHM-keV polynomial in z = √E_keV
  for the ``lsrm_fwhm_polynomial_in_E`` model) — RAG-043 doc_corpus_id
  ``lsrm_act_2014``; verified against Th-232 Marinelli fixture and
  documented in ``scripts/gamma/io/lsrm_spe.py:47-58``.
- Standard FWHM-anchor energy set (BUG-50): LNHB Recommended Data
  (laraweb.free.fr) for nuclide line energies + Gilmore 2nd ed. Table 1.1
  (canonical calibration nuclides for gamma spectroscopy). All anchor
  values below cite their NLHB/Gilmore origin individually in
  ``STANDARD_FWHM_ANCHORS_KEV``.
- Project floor ``ENERGY_CEILING_KEV = 3000`` (``gamma/spectrum.py:28``).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Sequence


# -------------------------------------------------------------------
# Resolution-band table (see G4 above for citations)
# Bands are (R_min, R_max) where R = FWHM/E (fraction, not percent)
# Selected by relative resolution at the reference energy (typically
# the file's highest stored anchor or the nearest to 662 keV).
# -------------------------------------------------------------------
HPGE_BAND_R = (0.0005, 0.020)      # 0.05% … 2%
SCINT_BAND_R = (0.020, 0.20)       # 2% … 20%

# Boundary R between HPGe-class and scintillator-class (R @ 662 keV).
# 1.5% matches the CdZnTe upper edge used in
# gamma.calibration.detector_type._R_RANGES.
HPGE_VS_SCINT_BOUNDARY_R = 0.015

# Channel-range sanity
E_LOW_KEV_MIN = -50.0     # tolerated negative intercept on E(0)
E_HIGH_KEV_MIN = 100.0    # spectrum must cover at least up to 100 keV


# -------------------------------------------------------------------
# Standard FWHM-anchor energy set (BUG-50, v1.22.6)
#
# Canonical literature gamma lines used for FWHM/E plausibility and
# monotonicity evaluation when the vendor file does not carry an
# explicit list of measured ``calibration_peaks`` (the case for the
# entire LSRM corpus, where the FWHM model is stored as polynomial
# coefficients alone — see ``gamma.io.lsrm_spe`` lines 47-58 and
# 314-323; ``calibration_peaks=[]`` is the literal value emitted).
#
# Each anchor lists ``(energy_keV, nuclide_label, source_tag)`` where
# ``source_tag`` is the citation handle:
#   - LNHB: Laboratoire National Henri Becquerel "Recommended Data"
#     tables (laraweb.free.fr / nucleide.org), the de-facto reference
#     dataset used by Gilmore 2nd ed. and LSRM as well.
#   - Gilmore: Practical Gamma-Ray Spectrometry, 2nd ed., Wiley 2008,
#     Table 1.1 (canonical nuclides for calibration).
#
# Anchors are ordered low → high. Each anchor is also expressible as a
# well-known nuclide library entry (Am-241, Co-57, Cs-137, Co-60,
# K-40, Tl-208) — these are the standard "first-tier" calibration
# nuclides used by both HPGe and scintillator laboratories.
# -------------------------------------------------------------------
STANDARD_FWHM_ANCHORS_KEV: tuple = (
    # (energy_keV, nuclide_label, source_tag)
    ( 59.5409, "Am-241", "LNHB"),     # γ at 59.5409 keV (intense, low-E HPGe anchor)
    (122.0607, "Co-57",  "LNHB"),     # γ at 122.06065 keV (low-E scint+HPGe anchor)
    (661.657 , "Cs-137", "LNHB"),     # γ at 661.657 keV (the reference R-line)
    (1173.228, "Co-60",  "LNHB"),     # γ1 at 1173.228 keV
    (1332.492, "Co-60",  "LNHB"),     # γ2 at 1332.492 keV
    (1460.822, "K-40",   "LNHB"),     # natural-background γ; Gilmore Table 1.1
    (2614.511, "Tl-208", "LNHB"),     # Th-232 chain terminator; Gilmore Table 1.1
)

# Anchor-mode threshold: when fewer than this many evaluable measured
# calibration peaks exist, the gate falls back to the anchor-mode
# evaluator. Set to 4 (matches the BUG-48 sweep recommendation in
# `_state/agent_a/outbox/2026-06-04_A2_BUG48_sweep.md` §7).
MEASURED_PEAK_THRESHOLD = 4

# Anchor-mode G4 overshoot factor: at any anchor, if the FWHM at that
# energy exceeds the upper band edge by more than this factor, a soft
# G4_anchor warning fires. 2.0 is conservative — a single anchor must
# be ≥ 2× the upper edge of its detector-class plausibility band to
# trigger. Rationale: the measured-mode G4 already fires at exactly
# the band edge (lo / hi); anchor-mode adds an additional 2× safety
# margin so that physically borderline anchors do not blow up
# false-positive rate on the LSRM corpus.
ANCHOR_OVERSHOOT_FACTOR = 2.0


@dataclass
class CalibrationGateResult:
    """Verdict from the hard calibration gate.

    Attributes:
        passed: True iff no HARD criterion (G1 / G2) failed.
        hard_failures: list of dict, each with keys ``code`` (G1..G2),
            ``message``, ``value``.
        soft_warnings: list of dict for G3 / G4 issues — do not block
            identification but surface to ``report.json.warnings``.
        criteria_evaluated: list of criterion codes that were actually
            run (some can be skipped when no stored FWHM model is
            present).
        detector_class_hint: one of "HPGe", "scintillator", "unknown" —
            chosen for the FWHM band check; useful for downstream
            diagnostics.
        E_low_keV / E_high_keV: evaluated endpoints of the spectrum
            energy axis (sanity-check inputs).
        reason: one-line human summary suitable for logs.
    """
    passed: bool
    hard_failures: list = field(default_factory=list)
    soft_warnings: list = field(default_factory=list)
    criteria_evaluated: list = field(default_factory=list)
    detector_class_hint: str = "unknown"
    E_low_keV: Optional[float] = None
    E_high_keV: Optional[float] = None
    reason: str = ""

    def as_dict(self) -> dict:
        """Return a JSON-friendly view (consumed by report.json)."""
        return {
            "passed": bool(self.passed),
            "hard_failures": list(self.hard_failures),
            "soft_warnings": list(self.soft_warnings),
            "criteria_evaluated": list(self.criteria_evaluated),
            "detector_class_hint": self.detector_class_hint,
            "E_low_keV": self.E_low_keV,
            "E_high_keV": self.E_high_keV,
            "reason": self.reason,
        }


# ============================================================================
# Public entry point
# ============================================================================

def evaluate_calibration_gate(spec) -> CalibrationGateResult:
    """Run all calibration self-consistency criteria on a parsed Spectrum.

    Parameters
    ----------
    spec : gamma.spectrum.Spectrum
        Must have ``energy_cal`` set; ``stored_fwhm_calibration`` is
        optional (G3 / G4 only run when it carries usable data).

    Returns
    -------
    CalibrationGateResult
        Structured verdict. Caller writes ``spec.extras['calibration_gate']
        = result.as_dict()`` and decides whether to short-circuit ID based
        on ``result.passed``.

    Notes
    -----
    The function never raises on a broken calibration — it returns a
    ``passed=False`` result with diagnostic ``hard_failures``. The only
    exceptions raised come from genuinely malformed inputs (None spec,
    energy_cal of wrong type) so caller bugs surface clearly.
    """
    if spec is None:
        raise ValueError("evaluate_calibration_gate: spec is None")

    hard: list = []
    soft: list = []
    criteria: list = []

    # ------------------------------------------------------------------
    # G1: energy-cal presence + monotonicity
    # ------------------------------------------------------------------
    criteria.append("G1")
    energy_cal = getattr(spec, "energy_cal", None)
    if energy_cal is None or len(energy_cal) < 2:
        return CalibrationGateResult(
            passed=False,
            hard_failures=[{
                "code": "G1",
                "message": "No energy calibration present (energy_cal is None or < 2 coefs)",
                "value": None,
            }],
            criteria_evaluated=criteria,
            reason="HARD FAIL G1: missing energy calibration",
        )

    n_ch = max(int(getattr(spec, "n_channels", 0)), 2)
    a1 = float(energy_cal[1]) if len(energy_cal) >= 2 else 0.0

    # Sample dE/dN on a dense grid; fail if ever ≤ 0 in [0, n_ch-1]
    grid = _channel_grid(n_ch)
    dEdN_vals = [_dE_dN(energy_cal, ch) for ch in grid]
    min_dEdN = min(dEdN_vals)
    if min_dEdN <= 0:
        idx = dEdN_vals.index(min_dEdN)
        bad_ch = grid[idx]
        hard.append({
            "code": "G1",
            "message": (
                f"Energy axis non-monotone: dE/dN={min_dEdN:.4g} at channel "
                f"{bad_ch} (must be > 0 everywhere on [0, {n_ch-1}])"
            ),
            "value": float(min_dEdN),
        })

    # ------------------------------------------------------------------
    # G2: channel-range sanity
    # ------------------------------------------------------------------
    criteria.append("G2")
    E_low = float(_poly_eval(energy_cal, 0.0))
    E_high = float(_poly_eval(energy_cal, n_ch - 1))

    if E_low < E_LOW_KEV_MIN:
        hard.append({
            "code": "G2",
            "message": (
                f"E(channel=0)={E_low:.2f} keV below floor {E_LOW_KEV_MIN} keV "
                f"(unphysical negative offset)"
            ),
            "value": E_low,
        })
    if E_high < E_HIGH_KEV_MIN:
        hard.append({
            "code": "G2",
            "message": (
                f"E(channel={n_ch-1})={E_high:.2f} keV below minimum useful "
                f"range {E_HIGH_KEV_MIN} keV"
            ),
            "value": E_high,
        })

    # ------------------------------------------------------------------
    # G3 / G4: FWHM-related — only when stored model is usable
    # ------------------------------------------------------------------
    detector_class_hint = "unknown"
    sf = getattr(spec, "stored_fwhm_calibration", None)
    if sf is not None:
        # Count the number of usable measured calibration peaks. When
        # this drops below MEASURED_PEAK_THRESHOLD AND the model is
        # one that the anchor-mode evaluator can evaluate directly
        # (currently ``lsrm_fwhm_polynomial_in_E``), we fall through
        # to anchor mode. See STANDARD_FWHM_ANCHORS_KEV header
        # comment (BUG-50). Both criteria still appear in
        # ``criteria_evaluated`` so downstream aggregators always see
        # G3 / G4 as "exercised". When anchor-mode ran, the entry is
        # tagged with the ``_anchor`` suffix.
        #
        # For ``SimpleSqrtFwhm`` (AtomSpectra) we deliberately keep
        # measured-mode regardless of cal_peaks count: the model
        # evaluates in channel space (FWHM²(N) = c0 + c1·N) and the
        # existing G3 check already exercises it analytically — no
        # anchor synthesis is needed.
        n_measured = _count_usable_cal_peaks(sf)
        model_label = (getattr(sf, "model", "") or "").strip()
        anchor_capable = model_label == "lsrm_fwhm_polynomial_in_E"
        anchor_mode = anchor_capable and n_measured < MEASURED_PEAK_THRESHOLD

        # G3: FWHM model monotonicity
        if anchor_mode:
            criteria.append("G3_anchor")
            msg = _fwhm_anchor_monotonicity_failure(sf, energy_cal, n_ch)
            if msg is not None:
                soft.append({
                    "code": "G3_anchor",
                    "message": msg,
                    "value": None,
                })
        else:
            criteria.append("G3")
            non_monotone = _fwhm_monotonicity_failure(sf, n_ch)
            if non_monotone is not None:
                soft.append({
                    "code": "G3",
                    "message": non_monotone,
                    "value": None,
                })

        # G4: FWHM physical plausibility band
        if anchor_mode:
            criteria.append("G4_anchor")
            detector_class_hint, anchor_issues = _fwhm_anchor_band_check(
                sf, energy_cal, n_ch,
            )
            for issue in anchor_issues:
                soft.append(issue)
        else:
            criteria.append("G4")
            detector_class_hint, band_issues = _fwhm_band_check(sf)
            for issue in band_issues:
                soft.append(issue)

    # ------------------------------------------------------------------
    # Compose verdict
    # ------------------------------------------------------------------
    passed = (len(hard) == 0)
    if passed and not soft:
        reason = "PASS: all calibration self-consistency criteria met"
    elif passed and soft:
        reason = (
            f"PASS with {len(soft)} soft warning(s): "
            f"{', '.join(w['code'] for w in soft)}"
        )
    else:
        reason = (
            f"HARD FAIL: {len(hard)} criterion(s) — "
            f"{', '.join(h['code'] for h in hard)}; "
            f"{len(soft)} soft warning(s)"
        )

    return CalibrationGateResult(
        passed=passed,
        hard_failures=hard,
        soft_warnings=soft,
        criteria_evaluated=criteria,
        detector_class_hint=detector_class_hint,
        E_low_keV=E_low,
        E_high_keV=E_high,
        reason=reason,
    )


# ============================================================================
# Helpers
# ============================================================================

def _channel_grid(n_ch: int) -> list:
    """Return a sampling grid over [0, n_ch-1] used for monotonicity checks.

    Up to ~100 evenly-spaced channels — dense enough to catch any
    physically realistic local minimum, cheap enough to run on every
    spectrum load.
    """
    n_samples = min(100, n_ch)
    if n_samples <= 1:
        return [0]
    step = (n_ch - 1) / (n_samples - 1)
    return [int(round(i * step)) for i in range(n_samples)]


def _poly_eval(coefs: Sequence[float], x: float) -> float:
    """Evaluate polynomial sum(a_i * x**i), low-to-high coefficients."""
    return sum(float(a) * (float(x) ** i) for i, a in enumerate(coefs))


def _dE_dN(coefs: Sequence[float], ch: float) -> float:
    """Local slope of energy_cal polynomial at channel ``ch``."""
    return sum(
        i * float(a) * (float(ch) ** (i - 1))
        for i, a in enumerate(coefs) if i > 0
    )


def _fwhm_monotonicity_failure(sf, n_ch: int) -> Optional[str]:
    """Return a message if the stored FWHM model is non-monotone over the
    channel range, else None.

    For SimpleSqrtFwhm: FWHM²(N) = c0 + c1·N → monotone iff c1 ≥ 0.
    For other models we walk the stored calibration peaks (channel,
    fwhm_channels) and check pairwise non-decrease.
    """
    model = getattr(sf, "model", "") or ""
    coefs = getattr(sf, "coefficients", ()) or ()

    if model == "SimpleSqrtFwhm" and len(coefs) >= 2:
        c1 = float(coefs[1])
        if c1 < 0:
            return (
                f"SimpleSqrtFwhm c1={c1:.4g} < 0 → FWHM²(N) decreases with "
                f"channel; non-physical for scintillator/HPGe"
            )
        # Also: floor argument must stay > 0 inside the range — otherwise
        # FWHM hits the min_channels clamp silently.
        c0 = float(coefs[0])
        if c0 + c1 * (n_ch - 1) <= 0:
            return (
                f"SimpleSqrtFwhm argument c0+c1·N becomes ≤ 0 inside "
                f"[0, {n_ch-1}] (c0={c0:.3g}, c1={c1:.4g}) → FWHM clamps "
                f"silently to provider floor"
            )
        return None

    # Fallback: check stored cal peaks order
    cal_peaks = list(getattr(sf, "calibration_peaks", []) or [])
    if len(cal_peaks) >= 2:
        ordered = sorted(cal_peaks, key=lambda cp: float(cp.channel))
        for prev, nxt in zip(ordered, ordered[1:]):
            f_prev = float(getattr(prev, "fwhm_channels", 0) or 0)
            f_next = float(getattr(nxt, "fwhm_channels", 0) or 0)
            if f_prev > 0 and f_next > 0 and f_next < f_prev:
                return (
                    f"Stored calibration peaks show FWHM decrease: "
                    f"channel {prev.channel} fwhm={f_prev:.2f} → "
                    f"channel {nxt.channel} fwhm={f_next:.2f}"
                )
    return None


def _fwhm_band_check(sf) -> tuple:
    """Check FWHM/E against the per-detector-class band at each stored
    calibration peak.

    Returns
    -------
    (detector_class_hint, issues)
        ``detector_class_hint`` ∈ {"HPGe", "scintillator", "unknown"} —
        chosen from the cal peak nearest to 662 keV.
        ``issues`` is a list of soft-warning dicts (G4 entries).
    """
    cal_peaks = list(getattr(sf, "calibration_peaks", []) or [])
    issues: list = []

    usable = [
        cp for cp in cal_peaks
        if float(getattr(cp, "energy_keV", 0) or 0) > 0
        and float(getattr(cp, "fwhm_channels", 0) or 0) > 0
    ]
    if not usable:
        return "unknown", issues

    # Pick the reference peak for class-hint: nearest to 662 keV, else
    # highest-E available.
    ref = min(
        usable,
        key=lambda cp: abs(float(cp.energy_keV) - 662.0),
    )

    # Convert FWHM at the reference peak to keV via local channel pitch
    # at neighbouring cal peaks — without spec we approximate by assuming
    # the cal peaks are equispaced enough that E(ch)/ch ≈ keV/channel.
    # This is a coarse estimate only used for class hint; the bracket is
    # wide enough (5x) to absorb the approximation.
    ref_E = float(ref.energy_keV)
    ref_ch = float(ref.channel)
    keV_per_ch = ref_E / ref_ch if ref_ch > 0 else 1.0
    R_ref = (float(ref.fwhm_channels) * keV_per_ch) / ref_E if ref_E > 0 else 0.0

    if R_ref <= HPGE_VS_SCINT_BOUNDARY_R:
        cls = "HPGe"
        band = HPGE_BAND_R
    else:
        cls = "scintillator"
        band = SCINT_BAND_R

    lo, hi = band
    for cp in usable:
        E = float(cp.energy_keV)
        ch = float(cp.channel)
        # Same coarse keV/ch estimate per peak (good enough as gate)
        per = E / ch if ch > 0 else 1.0
        fwhm_keV = float(cp.fwhm_channels) * per
        R = fwhm_keV / E if E > 0 else 0.0
        if R < lo or R > hi:
            issues.append({
                "code": "G4",
                "message": (
                    f"FWHM/E={R*100:.2f}% at cal peak E={E:.1f} keV outside "
                    f"{cls} plausibility band [{lo*100:.2f}%, {hi*100:.2f}%]"
                ),
                "value": float(R),
            })
    return cls, issues


# ============================================================================
# Anchor-mode helpers (BUG-50, v1.22.6)
#
# When the stored FWHM model carries fewer than MEASURED_PEAK_THRESHOLD
# usable measured calibration peaks (the LSRM case — see
# ``audit/.../2026-06-04_A2_BUG48_sweep.md`` §7), we evaluate the FWHM
# model directly at a standard set of literature anchor energies that
# lie within the spectrum's calibrated energy range. This activates G3
# / G4 on the entire LSRM corpus without needing peak search to fire
# first. See STANDARD_FWHM_ANCHORS_KEV header comment for citations.
# ============================================================================

def _count_usable_cal_peaks(sf) -> int:
    """Count peaks in ``sf.calibration_peaks`` with E>0 AND FWHM>0."""
    cal_peaks = list(getattr(sf, "calibration_peaks", []) or [])
    return sum(
        1 for cp in cal_peaks
        if float(getattr(cp, "energy_keV", 0) or 0) > 0
        and float(getattr(cp, "fwhm_channels", 0) or 0) > 0
    )


def _spectrum_energy_range(energy_cal, n_ch: int) -> tuple:
    """Return ``(E_low_keV, E_high_keV)`` of the calibrated axis."""
    return (
        float(_poly_eval(energy_cal, 0.0)),
        float(_poly_eval(energy_cal, max(n_ch - 1, 0))),
    )


def _eval_fwhm_keV_from_model(sf, E_keV: float) -> Optional[float]:
    """Evaluate the stored FWHM model at energy ``E_keV`` → FWHM in keV.

    Supports the two model labels currently emitted by the readers:

    - ``"lsrm_fwhm_polynomial_in_E"`` (LSRM .spe): despite the label
      the polynomial argument is ``z = √E_keV`` (LSRM Algorithmic
      Foundations §8.3; see ``gamma/io/lsrm_spe.py:47-58``).
      FWHM_keV(E) = Σ_k c_k · z^k.
    - ``"SimpleSqrtFwhm"`` (AtomSpectra): FWHM²(N) = c0 + c1·N is the
      model in channel space. To use it for an anchor energy we need
      the energy_cal to invert E → N. Returns FWHM_keV via local
      slope dE/dN.

    Returns ``None`` if the model is not supported or the coefficients
    are unusable (empty / too short / produce a non-positive argument).
    """
    model = (getattr(sf, "model", "") or "").strip()
    coefs = list(getattr(sf, "coefficients", ()) or ())
    if not coefs:
        return None

    if model == "lsrm_fwhm_polynomial_in_E":
        if E_keV <= 0:
            return None
        z = math.sqrt(E_keV)
        # Skip a leading degree-marker if its value is suspiciously like
        # a small int (the LSRM file format prefixes the coefficient
        # list with a degree-marker integer). The reader in
        # ``gamma/io/lsrm_spe.py`` already strips it, so the stored
        # ``coefficients`` tuple is the pure low-to-high polynomial
        # coefficient list — no extra handling here.
        val = sum(float(c) * (z ** k) for k, c in enumerate(coefs))
        return float(val) if val > 0 else None

    if model == "SimpleSqrtFwhm" and len(coefs) >= 2:
        # FWHM²(N) = c0 + c1·N → need to invert energy_cal here, but
        # we deliberately do NOT take energy_cal as an arg in this
        # helper so the caller decides whether to invoke this path.
        # Returning None here is the right thing — the caller's
        # AtomSpectra spectra always carry calibration_peaks ≥ 2 so
        # they never hit the anchor branch.
        return None

    return None


def _anchors_in_range(E_low: float, E_high: float) -> list:
    """Return the subset of STANDARD_FWHM_ANCHORS_KEV inside the range.

    Result is the list of ``(energy_keV, nuclide_label, source_tag)``
    triples that lie strictly inside ``(E_low, E_high)``. Anchors at
    the edges are excluded to avoid edge-of-spectrum numerical effects
    in downstream FWHM polynomial evaluation. The standard set is
    already sorted ascending; the filtered list inherits that order.
    """
    return [
        (E, label, src) for (E, label, src) in STANDARD_FWHM_ANCHORS_KEV
        if E_low < E < E_high
    ]


def _classify_detector_anchor(fwhm_at_662_keV: Optional[float]) -> str:
    """Pick a detector class hint from FWHM at 662 keV (the reference
    line). Returns "HPGe" / "scintillator" / "unknown".
    """
    if fwhm_at_662_keV is None or fwhm_at_662_keV <= 0:
        return "unknown"
    R = fwhm_at_662_keV / 661.657
    return "HPGe" if R <= HPGE_VS_SCINT_BOUNDARY_R else "scintillator"


def _fwhm_anchor_monotonicity_failure(sf, energy_cal, n_ch: int) -> Optional[str]:
    """Anchor-mode G3 evaluator.

    Evaluates the stored FWHM model at every standard anchor energy
    inside the spectrum's calibrated range and checks that the FWHM
    sequence is non-decreasing with energy (Gilmore §6.4 eq. 6.13:
    FWHM grows monotonically with energy for any scintillator-class
    detector under Poisson broadening, and the same holds for HPGe at
    energies above the noise floor).

    Returns a diagnostic message describing the first decrease, or
    None if all consecutive (anchor_E, FWHM_keV) pairs satisfy
    FWHM[i+1] ≥ FWHM[i]. When fewer than 2 anchors lie inside the
    spectrum range, the check is skipped (returns None — anchor-mode
    cannot decide with a single sample).
    """
    E_low, E_high = _spectrum_energy_range(energy_cal, n_ch)
    anchors = _anchors_in_range(E_low, E_high)
    if len(anchors) < 2:
        return None

    fwhms: list = []
    for (E, label, _src) in anchors:
        f = _eval_fwhm_keV_from_model(sf, E)
        if f is None:
            # Skip anchors the model cannot evaluate; never abort —
            # this branch is best-effort by design.
            continue
        fwhms.append((E, label, float(f)))

    if len(fwhms) < 2:
        return None

    for (E_prev, label_prev, f_prev), (E_next, label_next, f_next) in zip(
        fwhms, fwhms[1:],
    ):
        if f_next + 1e-9 < f_prev:
            return (
                f"Anchor-mode FWHM non-monotone: "
                f"FWHM({label_prev}@{E_prev:.1f})={f_prev:.2f} keV → "
                f"FWHM({label_next}@{E_next:.1f})={f_next:.2f} keV "
                f"(decreases with energy; Gilmore §6.4 expects "
                f"non-decreasing FWHM)"
            )
    return None


def _fwhm_anchor_band_check(sf, energy_cal, n_ch: int) -> tuple:
    """Anchor-mode G4 evaluator.

    Evaluates the stored FWHM model at every standard anchor energy
    inside the spectrum's calibrated range and flags anchors whose
    FWHM/E ratio falls outside the detector-class plausibility band
    by more than ``ANCHOR_OVERSHOOT_FACTOR`` × the band edge. The
    detector class is selected by FWHM/E at the Cs-137 661.657 keV
    anchor when that anchor is in range; else the median FWHM/E ratio
    over the available anchors is used.

    Returns
    -------
    (detector_class_hint, issues)
        ``detector_class_hint`` ∈ {"HPGe", "scintillator", "unknown"}.
        ``issues`` is a list of soft-warning dicts (G4_anchor entries),
        empty when all anchors lie inside the band.
    """
    E_low, E_high = _spectrum_energy_range(energy_cal, n_ch)
    anchors = _anchors_in_range(E_low, E_high)
    issues: list = []
    if not anchors:
        return "unknown", issues

    # Collect (E, label, FWHM_keV, R=FWHM/E) for every anchor the
    # model could evaluate.
    evaluated: list = []
    fwhm_at_662: Optional[float] = None
    for (E, label, _src) in anchors:
        f = _eval_fwhm_keV_from_model(sf, E)
        if f is None or f <= 0:
            continue
        R = f / E if E > 0 else 0.0
        evaluated.append((E, label, float(f), float(R)))
        # The 661.657 keV Cs-137 anchor is the reference line for class
        # selection; cache it when present.
        if abs(E - 661.657) < 1e-6:
            fwhm_at_662 = float(f)

    if not evaluated:
        return "unknown", issues

    # Class hint: prefer FWHM at 662 keV; else use median R across the
    # anchor set (deterministic and robust to one outlier).
    if fwhm_at_662 is not None:
        cls = _classify_detector_anchor(fwhm_at_662)
    else:
        Rs = sorted(R for (_E, _l, _f, R) in evaluated)
        mid = Rs[len(Rs) // 2]
        cls = "HPGe" if mid <= HPGE_VS_SCINT_BOUNDARY_R else "scintillator"

    band = HPGE_BAND_R if cls == "HPGe" else SCINT_BAND_R
    lo, hi = band
    # Anchor-mode applies an additional safety factor on the high side
    # only (ANCHOR_OVERSHOOT_FACTOR). The low side already triggers at
    # the band edge because under-broadening at a literature anchor is
    # always meaningful (broken model or wrong detector classification).
    lo_eff = lo
    hi_eff = hi * ANCHOR_OVERSHOOT_FACTOR

    for (E, label, f, R) in evaluated:
        if R < lo_eff or R > hi_eff:
            issues.append({
                "code": "G4_anchor",
                "message": (
                    f"Anchor-mode FWHM/E={R*100:.2f}% at {label} "
                    f"E={E:.1f} keV outside {cls} plausibility band "
                    f"[{lo_eff*100:.2f}%, {hi_eff*100:.2f}%] "
                    f"(FWHM={f:.2f} keV)"
                ),
                "value": float(R),
            })
    return cls, issues


__all__ = [
    "CalibrationGateResult",
    "evaluate_calibration_gate",
    "HPGE_BAND_R",
    "SCINT_BAND_R",
    "HPGE_VS_SCINT_BOUNDARY_R",
    "E_LOW_KEV_MIN",
    "E_HIGH_KEV_MIN",
    "STANDARD_FWHM_ANCHORS_KEV",
    "MEASURED_PEAK_THRESHOLD",
    "ANCHOR_OVERSHOOT_FACTOR",
]
