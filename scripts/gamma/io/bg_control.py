"""
F-243 BG control: pre-analysis check of background spectrum quality.

Three independent gates per ISO 11929:2019 §9.2 and Gilmore & Joss 3rd Ed. §5.6:

    1.  |F - F_ref| / F_ref < rate_tolerance     (default 10 %)
    2.  sum(Y) >= min_sum_counts                  (default 1000 counts)
    3.  t_live >= min_live_time_s                 (default 600 s)

The function is purely synchronous, dependency-free, and returns a frozen
dataclass. It is intentionally NON-BLOCKING — the caller decides whether
to warn-only or to abort. The default pipeline integration (Phase 3 of
v1.18.29) wires it as a warning-only gate that appends a line to
``pipeline_notes`` when any of the three gates fail.

RAG-ID:
    [LSRM-Algo-BG-Stability], [ISO-11929-§9.2], [Gilmore-§5.6]
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field


@dataclass(frozen=True)
class BgControlResult:
    """Outcome of a single background quality check.

    Attributes
    ----------
    ok : bool
        True iff all three gates passed.
    rate_ratio : float
        ``F / F_ref`` (NaN when ``F_ref <= 0`` and the ratio check is skipped).
    sum_counts : float
        Total counts in the background spectrum (``sum(Y)``).
    live_time_s : float
        Background live time in seconds.
    failures : tuple
        Tuple of human-readable strings, one per failed gate. Empty when
        ``ok=True``.
    """

    ok: bool
    rate_ratio: float
    sum_counts: float
    live_time_s: float
    failures: tuple


def validate_background(
    F: float,
    F_ref: float,
    sum_Y: float,
    t_live: float,
    *,
    rate_tolerance: float = 0.1,
    min_sum_counts: float = 1000.0,
    min_live_time_s: float = 600.0,
) -> BgControlResult:
    """Validate a background spectrum against three independent quality gates.

    Parameters
    ----------
    F : float
        Sample count rate (cps).
    F_ref : float
        Background reference count rate (cps). When ``<= 0`` the rate
        consistency check is skipped (``rate_ratio = NaN``) and only the
        absolute checks remain.
    sum_Y : float
        Total counts in the background spectrum.
    t_live : float
        Background live time in seconds.
    rate_tolerance : float, optional
        Maximum permissible ``|F - F_ref| / F_ref`` (default 0.10).
    min_sum_counts : float, optional
        Minimum ``sum(Y)`` (default 1000).
    min_live_time_s : float, optional
        Minimum live time in seconds (default 600).

    Returns
    -------
    BgControlResult
        ``ok=True`` iff all checked gates pass. The ``failures`` tuple
        lists the human-readable reason for every gate that failed.
    """
    F = float(F)
    F_ref = float(F_ref)
    sum_Y = float(sum_Y)
    t_live = float(t_live)
    rate_tolerance = float(rate_tolerance)
    min_sum_counts = float(min_sum_counts)
    min_live_time_s = float(min_live_time_s)

    failures: list[str] = []

    if F_ref > 0:
        ratio = F / F_ref
        diff = abs(F - F_ref) / F_ref
        if diff > rate_tolerance:
            failures.append(
                f"rate mismatch: |F-F0|/F0={diff:.3f} > {rate_tolerance:.3f}"
            )
    else:
        ratio = float("nan")

    if sum_Y < min_sum_counts:
        failures.append(
            f"sum_counts={sum_Y:.0f} < {min_sum_counts:.0f}"
        )
    if t_live < min_live_time_s:
        failures.append(
            f"t_live={t_live:.0f}s < {min_live_time_s:.0f}s"
        )

    return BgControlResult(
        ok=(len(failures) == 0),
        rate_ratio=ratio,
        sum_counts=sum_Y,
        live_time_s=t_live,
        failures=tuple(failures),
    )


def validate_background_roi(
    F_roi: float,
    F_ref_roi: float,
    sum_Y_roi: float,
    t_live: float,
    *,
    rate_tolerance: float = 0.1,
    min_sum_counts: float = 1000.0,
    min_live_time_s: float = 600.0,
) -> BgControlResult:
    """Validate background quality for a single ROI energy window.

    Scoped variant of ``validate_background()`` — applies the same three-gate
    check but to per-ROI inputs (count rate and counts within [E_lo, E_hi]).
    Delegates directly; no logic duplication. RAG-ID: RAG-008.
    """
    return validate_background(
        F=float(F_roi),
        F_ref=float(F_ref_roi),
        sum_Y=float(sum_Y_roi),
        t_live=float(t_live),
        rate_tolerance=rate_tolerance,
        min_sum_counts=min_sum_counts,
        min_live_time_s=min_live_time_s,
    )


# ════════════════════════════════════════════════════════════════════
# BUG-35 — Statistical |z|-test gate (per-peak ROI bg stability)
# Per ISO 11929-2:2019 §6 + Currie 1968 Eq 17 (RAG-009 cite-list).
# Self-correcting successor to the RAG-005/008 engineering 10% rule.
# ════════════════════════════════════════════════════════════════════
#
# Formula (single-count Poisson form, ISO 11929-2:2019 §6 — see
# audit/_rag/RAG_INDEX.json RAG-005 source iso_11929_2_2019 verbatim
# quote, and RAG-009 cite-list line "Criterion 4 — Poisson z-test"):
#
#     z = (B1 − B2) / √(B1 + B2)
#
# where B1, B2 are integer Poisson counts in two ROIs (sample bg-window
# vs reference bg). Tier system per methodology v2 §criterion 5
# (audit/_drafts/spectrum_qc_methodology_v2_2026-06-03.md, transcribed
# 2026-06-03):
#
#   * |z| < 2.0  → "stable"     (PASS)
#   * 2.0 ≤ |z| < 3.0 → "borderline" (recount recommended; soft PASS)
#   * |z| ≥ 3.0  → "reject"     (FAIL — non-stationarity per ISO §6)
#
# Self-correcting w.r.t. count level: at strong bg the absolute tolerable
# difference grows as √(B1+B2); at weak bg the same z-threshold loosens
# the absolute tolerance. This is the property the 10% engineering rule
# lacks (it is too tight for weak bg, too loose for strong bg).
#
# Knoll 'Radiation Detection and Measurement' 4th ed. §3.IV.B documents
# the Gaussian-normal-deviate interpretation that grounds the 2/3 tier
# boundaries (0.13 % one-sided tail at |z|=3.0, 2.3 % at |z|=2.0).


# Tier boundaries — pinned literals per methodology v2 §criterion 5.
# DO NOT change without RAG-022 update.
Z_TIER_STABLE_MAX = 2.0
Z_TIER_BORDERLINE_MAX = 3.0


@dataclass(frozen=True)
class BgZTestResult:
    """Outcome of a single Poisson |z|-test for bg stability.

    Implements ISO 11929-2:2019 §6 propagation:
        z = (B1 − B2) / √(B1 + B2)
    where B1, B2 are integer Poisson counts.

    Attributes
    ----------
    z : float
        Signed z-statistic ``(B1 − B2) / √(B1 + B2)``. NaN when
        ``B1 + B2 <= 0`` (no counts; test undefined).
    abs_z : float
        ``abs(z)`` or NaN.
    tier : str
        One of ``"stable"`` (|z| < 2.0), ``"borderline"`` (2.0 ≤ |z| < 3.0),
        ``"reject"`` (|z| ≥ 3.0), or ``"undefined"`` (B1+B2 ≤ 0).
    passed : bool
        ``True`` iff ``tier`` ∈ ``{"stable", "borderline"}``. Methodology
        v2 §5 treats borderline as a soft PASS — operator-visible warning
        but not a hard reject. ``"reject"`` and ``"undefined"`` are FAIL.
    B1 : int
        Sample bg-window counts (rounded if input was float).
    B2 : int
        Reference bg counts.
    note : str
        Human-readable diagnostic (empty when ``tier == "stable"``).
    """

    z: float
    abs_z: float
    tier: str
    passed: bool
    B1: int
    B2: int
    note: str


def bg_z_test(B1: float, B2: float) -> BgZTestResult:
    """Two-count Poisson |z|-test for bg stability per ISO 11929-2:2019 §6.

    Implements the single-count Poisson form of the §6 stationarity
    test cited verbatim in RAG-005 / RAG-008 / RAG-009::

        z = (B1 − B2) / √(B1 + B2)

    Tier boundaries from methodology v2 §criterion 5
    (`audit/_drafts/spectrum_qc_methodology_v2_2026-06-03.md`)::

        |z| < 2.0       → "stable"      (PASS, normal-deviate < 2.3 % tail)
        2.0 ≤ |z| < 3.0 → "borderline"  (soft PASS, recount recommended)
        |z| ≥ 3.0       → "reject"      (FAIL, non-stationarity per §6)

    Parameters
    ----------
    B1, B2 : float
        Integer Poisson counts in two background regions / measurements.
        Floats are accepted and rounded to ``int`` (counts must be whole
        for Poisson semantics — fractional inputs likely indicate a
        rate-vs-count confusion at the caller).

    Returns
    -------
    BgZTestResult
        Frozen result with signed ``z``, ``abs_z``, tier label, pass-flag,
        and human-readable note.

    Notes
    -----
    *Edge case — both counts zero or sum non-positive*: ``z`` is
    mathematically undefined; we return ``tier="undefined"`` with
    ``passed=False`` to signal the caller (do NOT silently treat zero
    counts as stable — that hides empty-spectrum bugs upstream).

    *Self-correction*: at strong bg (B1+B2=10⁶), |z|=3 allows ≈ 3 000
    count drift; at weak bg (B1+B2=100), it allows only ≈ 30 counts.
    The 10% engineering rule cannot do this (RAG-005/008 supersession).

    See Also
    --------
    validate_background : 10% rule kernel (kept for backwards compat).
    check_bg_quality : wrapper supporting both 10% and z-test gates.
    """
    B1 = float(B1)
    B2 = float(B2)
    n1 = int(round(B1))
    n2 = int(round(B2))

    total = B1 + B2
    if total <= 0.0:
        return BgZTestResult(
            z=float("nan"),
            abs_z=float("nan"),
            tier="undefined",
            passed=False,
            B1=n1,
            B2=n2,
            note="z-test undefined: B1+B2 ≤ 0 (empty spectrum or zeroed bg)",
        )

    z = (B1 - B2) / math.sqrt(total)
    abs_z = abs(z)

    if abs_z < Z_TIER_STABLE_MAX:
        tier = "stable"
        note = ""
    elif abs_z < Z_TIER_BORDERLINE_MAX:
        tier = "borderline"
        note = (
            f"BG |z|={abs_z:.2f} ∈ [2,3): borderline non-stationarity "
            f"(B1={n1}, B2={n2}). Recount recommended per ISO 11929-2:2019 §6."
        )
    else:
        tier = "reject"
        note = (
            f"BG |z|={abs_z:.2f} ≥ 3: non-stationary background "
            f"(B1={n1}, B2={n2}). Reject per ISO 11929-2:2019 §6."
        )

    return BgZTestResult(
        z=z,
        abs_z=abs_z,
        tier=tier,
        passed=tier in ("stable", "borderline"),
        B1=n1,
        B2=n2,
        note=note,
    )


# ════════════════════════════════════════════════════════════════════
# Wave 4 (2026-06-04) — rate-normalised z-test for unequal live-times
# Gilmore & Joss "Practical Gamma-Ray Spectrometry" 2nd Ed. §5.5
# (counting statistics; comparison of two count rates).
# ════════════════════════════════════════════════════════════════════
#
# `bg_z_test(B1, B2)` above assumes EQUAL live-time: it pretends the
# two integer count totals are directly comparable. In practice
# sample/bg spectra often have very unequal live-times (e.g. 3600 s vs
# 86400 s), and the raw-count z-test then biases toward "reject" simply
# because the longer measurement collected more counts.
#
# The rate-form propagates Poisson variance through the division by t:
#
#     R_i = c_i / t_i,   Var(R_i) = c_i / t_i²
#     z   = (R1 − R2) / √( Var(R1) + Var(R2) )
#         = (c1/t1 − c2/t2) / √( c1/t1² + c2/t2² )
#
# When t1 == t2 == t the formula reduces to
#     z = (c1 − c2) / √(c1 + c2)
# i.e. exactly the integer-count form of `bg_z_test`. The new function
# therefore generalises, not replaces, the existing API — the old
# signature is preserved for the existing pipeline gates.
#
# Cite: Gilmore & Joss §5.5 (eq. 5.21 — comparison of two count rates),
# F-157 (LSRM > Будыка > Gilmore — Gilmore is canonical here, no
# conflict with LSRM scope which does not specify this variant),
# F-243 (bg subtraction safety family).


def bg_z_test_rates(
    c1: float, t1: float, c2: float, t2: float
) -> tuple[float, bool]:
    """Rate-normalised z-test for two Poisson count rates with unequal live-times.

    For each spectrum we estimate the rate ``R_i = c_i / t_i`` (cps) and
    its Poisson variance ``Var(R_i) = c_i / t_i²``. The two-sample
    z-statistic for the difference of rates is then::

        z = (c1/t1 − c2/t2) / √(c1/t1² + c2/t2²)

    Returns a 2-tuple ``(z, is_significant)`` where ``is_significant``
    is ``|z| > 3.0`` (3-sigma threshold per the same tier convention
    as ``bg_z_test`` — see methodology v2 §criterion 5).

    Parameters
    ----------
    c1, c2 : float
        Integer Poisson counts in the two regions. Floats are accepted
        (e.g. ROI sums after numpy slicing) and are NOT rounded here —
        the rate formula is variance-correct for any non-negative real
        ``c_i``.
    t1, t2 : float
        Live-times in seconds. Both must be strictly positive.

    Returns
    -------
    (z, is_significant) : tuple[float, bool]
        ``z`` is the signed rate-difference z-statistic. ``NaN`` when
        ``c1 + c2 ≤ 0`` (no counts; variance undefined) or when either
        live-time is non-positive (degenerate). ``is_significant`` is
        ``True`` iff ``|z| > 3.0``; for ``NaN`` z it is ``False`` —
        callers wanting "fail-closed" semantics on undefined inputs
        should check ``math.isnan(z)`` separately.

    Notes
    -----
    *Equivalence with the integer-count form*: when ``t1 == t2 == t``,
    ``z`` here equals ``(c1 − c2) / √(c1 + c2)`` exactly — see
    `test_bg_control_z_test_rates.py::test_equal_time_matches_integer_form`.

    *Why not round c_i to int*: the rate-form variance is
    ``c/t²``, which is variance-correct for any non-negative real
    estimator of mean counts. Rounding would lose information from
    ROI sums that span partial channels.

    *Tier system*: this function returns the boolean 3-sigma gate only
    (matching the deferred Followup #3 spec in
    ``_state/agent_a/outbox/2026-06-04_backlog_top1_bug35_z_test.md``).
    For the full {stable / borderline / reject} tier the caller can
    compare ``abs(z)`` against ``Z_TIER_STABLE_MAX`` /
    ``Z_TIER_BORDERLINE_MAX``.

    See Also
    --------
    bg_z_test : integer-count form (assumes equal live-time).
    Z_TIER_STABLE_MAX, Z_TIER_BORDERLINE_MAX : tier boundaries.

    References
    ----------
    Gilmore & Joss, "Practical Gamma-Ray Spectrometry" 2nd Ed., §5.5.
    F-157 LSRM > Будыка > Gilmore precedence (no conflict here).
    """
    c1 = float(c1)
    c2 = float(c2)
    t1 = float(t1)
    t2 = float(t2)

    # Degenerate inputs → undefined z, not significant.
    if t1 <= 0.0 or t2 <= 0.0:
        return (float("nan"), False)
    if c1 + c2 <= 0.0:
        return (float("nan"), False)
    # Strict-negative counts are nonsensical for Poisson and would make
    # the variance term negative under the sqrt — fail closed.
    if c1 < 0.0 or c2 < 0.0:
        return (float("nan"), False)

    r1 = c1 / t1
    r2 = c2 / t2
    var = c1 / (t1 * t1) + c2 / (t2 * t2)
    if var <= 0.0:
        # Only reachable if c1 == c2 == 0 (handled above) — defensive.
        return (float("nan"), False)

    z = (r1 - r2) / math.sqrt(var)
    is_significant = abs(z) > Z_TIER_BORDERLINE_MAX  # 3.0
    return (z, is_significant)


# ════════════════════════════════════════════════════════════════════
# F-243 Hybrid layer (per AGENT_A_BRIEF_F-243.md, Agent D decision 2026-06-02)
# ════════════════════════════════════════════════════════════════════
#
# `validate_background` above is a low-level *numeric* kernel: scalars in,
# absolute thresholds, no domain knowledge. It stays for unit-tests and
# future composability.
#
# `check_bg_quality` below is the high-level *sample-relative* wrapper that
# the F-243 brief specifies — it accepts two Spectrum-like objects and
# evaluates three gates where the thresholds are tied to the sample being
# analysed (flux drift sample-vs-bg; Poisson σ of bg counts; bg live time
# as a fraction of sample live time). This is the API the pipeline uses.
#
# The wrapper does NOT call the kernel internally — the underlying physics
# differs (sample-relative vs absolute), and conflating them would mask
# either contract. Both APIs are exported.


def _sum_counts(x) -> float:
    """Lazy total-counts that works for list / tuple / numpy ndarray.

    Uses ``x.sum()`` if available (numpy fast path), otherwise falls back to
    ``sum(x)``. Returns ``0.0`` for empty / falsy inputs.
    """
    if x is None:
        return 0.0
    if hasattr(x, "sum") and hasattr(x, "size"):
        # numpy ndarray — also guards against empty
        return float(x.sum()) if int(x.size) else 0.0
    try:
        return float(sum(x))
    except TypeError:
        return 0.0


def _roi_sum_counts(counts, energy_cal: tuple, e_lo: float, e_hi: float) -> float:
    """Sum counts in channels whose calibrated energy is within [e_lo, e_hi].

    Horner evaluation: E(i) = sum(a_k * i**k for k, a_k in enumerate(energy_cal)).
    Uses numpy fast path when available; pure-Python fallback otherwise.
    """
    if not energy_cal or counts is None:
        return 0.0
    e_lo, e_hi = float(e_lo), float(e_hi)
    _ec = [float(a) for a in energy_cal]
    try:
        import numpy as np
        arr = np.asarray(counts, dtype=np.float64)
        n = int(arr.size)
        if n == 0:
            return 0.0
        idx = np.arange(n, dtype=np.float64)
        energies = np.zeros(n, dtype=np.float64)
        for k, a in enumerate(_ec):
            energies += a * idx ** k
        mask = (energies >= e_lo) & (energies <= e_hi)
        return float(arr[mask].sum())
    except ImportError:
        total = 0.0
        try:
            items = list(counts)
        except TypeError:
            return 0.0
        for i, c in enumerate(items):
            e = sum(a * (i ** k) for k, a in enumerate(_ec))
            if e_lo <= e <= e_hi:
                total += float(c)
        return total


@dataclass
class BgQualityReport:
    """High-level F-243 BG quality report (sample-relative).

    Mirrors the brief contract: ``passed`` is the AND of all gates,
    ``gates`` maps each gate name to ``{value, threshold, passed}``,
    ``notes`` lists one human-readable string per failing gate (empty
    when ``passed=True``).

    Note: not frozen — the brief specifies ``dict``/``list`` fields that
    are conceptually mutable. Treat as immutable downstream regardless.
    """

    passed: bool
    gates: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)


def check_bg_quality(
    sample_spec,
    bg_spec,
    *,
    flux_drift_tolerance: float = 0.10,
    poisson_sigma_tolerance: float = 0.10,
    live_time_min_ratio: float = 0.50,
    roi_keV_list=None,
    peak_z_roi_keV_list=None,
) -> BgQualityReport:
    """Sample-relative BG quality check (F-243 spec, three gates).

    Implements the high-level contract from `AGENT_A_BRIEF_F-243.md`:

      * **Gate 1 — flux_drift**
        ``|rate_sample − rate_bg| / rate_bg ≤ flux_drift_tolerance``
        where ``rate_x = sum(x.counts) / x.live_time``. Catches gross
        flux instability between sample and bg measurements.

      * **Gate 2 — sum_y_stat (Poisson σ of bg)**
        ``1 / sqrt(N_total_bg) ≤ poisson_sigma_tolerance``
        i.e. bg integral has at most a `tolerance` relative Poisson
        uncertainty. Default ``0.10`` → N ≥ 100 counts.

      * **Gate 3 — live_time_min**
        ``bg.live_time ≥ live_time_min_ratio × sample.live_time``
        adaptive floor — short sample ⇒ short bg is acceptable.

    Parameters
    ----------
    sample_spec, bg_spec : object
        Anything that exposes ``counts`` (iterable / ndarray) and
        ``live_time`` (float, seconds). Duck-typed — no isinstance check.
    flux_drift_tolerance, poisson_sigma_tolerance, live_time_min_ratio
        Per-gate thresholds, overridable for tests / tier-2 tuning.
    roi_keV_list : list of (float, float) or None, optional
        Per-ROI flux drift check (F-243.1 / RAG-008). Each ``(E_lo, E_hi)``
        in keV adds a ``"roi_{E_lo:.0f}_{E_hi:.0f}"`` gate that applies
        Gate 1 within that energy window. Requires ``energy_cal`` on both
        specs. ``None`` or ``[]`` → bit-identical to v1.18.29.
    peak_z_roi_keV_list : list of (float, float) or None, optional
        Per-peak-ROI ISO 11929-2:2019 §6 Poisson |z|-test gate (BUG-35 /
        RAG-022). Each ``(E_lo, E_hi)`` adds a
        ``"z_roi_{E_lo:.0f}_{E_hi:.0f}"`` gate where
        ``z = (B1 − B2) / √(B1 + B2)`` with B1 = sample-window counts,
        B2 = bg-window counts in the same energy window. Tier:
        |z| < 2 stable PASS, [2, 3) borderline soft-PASS,
        ≥ 3 reject FAIL. Self-correcting w.r.t. count level — replaces
        engineering 10% rule of RAG-005/008. Requires ``energy_cal`` on
        both specs. ``None`` or ``[]`` → gate not added (backwards compat).

    Returns
    -------
    BgQualityReport
        ``passed = all(g["passed"] for g in gates.values())``.
        ``gates`` keys: ``"flux_drift"``, ``"sum_y_stat"``,
        ``"live_time_min"`` (always); ``"roi_<lo>_<hi>"`` (when
        roi_keV_list is non-empty). ``notes`` lists one short message per
        failing gate.

    RAG-ID
    ------
    [LSRM-Algo-BG-Stability], [ISO-11929-§9.2], [Gilmore-§5.6], [RAG-008]
    """
    sample_counts = getattr(sample_spec, "counts", None)
    bg_counts     = getattr(bg_spec,     "counts", None)
    sample_tlive  = float(getattr(sample_spec, "live_time", 0.0) or 0.0)
    bg_tlive      = float(getattr(bg_spec,     "live_time", 0.0) or 0.0)

    n_total_sample = _sum_counts(sample_counts)
    n_total_bg     = _sum_counts(bg_counts)

    # ── Gate 1 — flux drift (sample-relative) ────────────────────────
    if sample_tlive > 0 and bg_tlive > 0 and n_total_bg > 0:
        rate_sample = n_total_sample / sample_tlive
        rate_bg     = n_total_bg     / bg_tlive
        flux_drift  = abs(rate_sample - rate_bg) / rate_bg
    else:
        flux_drift = float("inf")
    g1 = {
        "value": flux_drift,
        "threshold": float(flux_drift_tolerance),
        "passed": flux_drift <= flux_drift_tolerance,
    }

    # ── Gate 2 — sum_y_stat (Poisson σ of BG integral) ───────────────
    if n_total_bg > 0:
        sigma = 1.0 / math.sqrt(n_total_bg)
    else:
        sigma = float("inf")
    g2 = {
        "value": sigma,
        "threshold": float(poisson_sigma_tolerance),
        "passed": sigma <= poisson_sigma_tolerance,
    }

    # ── Gate 3 — live_time_min (sample-relative) ─────────────────────
    if sample_tlive > 0:
        ratio = bg_tlive / sample_tlive
    else:
        ratio = float("inf") if bg_tlive > 0 else 0.0
    g3 = {
        "value": ratio,
        "threshold": float(live_time_min_ratio),
        "passed": ratio >= live_time_min_ratio,
    }

    gates = {"flux_drift": g1, "sum_y_stat": g2, "live_time_min": g3}
    notes: list = []
    if not g1["passed"]:
        notes.append(
            f"BG flux drift {flux_drift*100:.1f}% > "
            f"{flux_drift_tolerance*100:.0f}%"
        )
    if not g2["passed"]:
        notes.append(
            f"BG Poisson σ={sigma*100:.1f}% > "
            f"{poisson_sigma_tolerance*100:.0f}% "
            f"(N_total_bg={n_total_bg:.0f})"
        )
    if not g3["passed"]:
        notes.append(
            f"BG live_time={ratio*100:.0f}% of sample "
            f"({bg_tlive:.0f}s / {sample_tlive:.0f}s) < "
            f"{live_time_min_ratio*100:.0f}%"
        )

    # ── ROI-windowed flux drift (F-243.1 / RAG-008) ──────────────────
    if roi_keV_list:
        sample_ecal = getattr(sample_spec, "energy_cal", None)
        bg_ecal     = getattr(bg_spec,     "energy_cal", None)
        if not sample_ecal or not bg_ecal:
            raise ValueError(
                "check_bg_quality: roi_keV_list requires energy_cal on "
                "both sample_spec and bg_spec"
            )
        for e_lo, e_hi in roi_keV_list:
            sum_roi_s  = _roi_sum_counts(sample_counts, sample_ecal, e_lo, e_hi)
            sum_roi_bg = _roi_sum_counts(bg_counts,     bg_ecal,     e_lo, e_hi)
            f_roi     = sum_roi_s  / sample_tlive if sample_tlive > 0 else 0.0
            f_ref_roi = sum_roi_bg / bg_tlive     if bg_tlive     > 0 else 0.0
            roi_res = validate_background_roi(
                F_roi=f_roi,
                F_ref_roi=f_ref_roi,
                sum_Y_roi=sum_roi_bg,
                t_live=bg_tlive,
                rate_tolerance=flux_drift_tolerance,
            )
            key = f"roi_{e_lo:.0f}_{e_hi:.0f}"
            gates[key] = {
                "value": roi_res.rate_ratio,
                "threshold": float(flux_drift_tolerance),
                "passed": roi_res.ok,
                "failures": roi_res.failures,
            }
            if not roi_res.ok:
                for msg in roi_res.failures:
                    notes.append(f"ROI [{e_lo:.0f},{e_hi:.0f} keV]: {msg}")

    # ── Per-peak ROI Poisson |z|-test (BUG-35 / RAG-022) ─────────────
    # ISO 11929-2:2019 §6: z = (B1 − B2) / √(B1 + B2).
    # B1 = sample-window counts, B2 = bg-window counts. Self-correcting
    # successor to the 10% rule in RAG-005/008.
    if peak_z_roi_keV_list:
        sample_ecal = getattr(sample_spec, "energy_cal", None)
        bg_ecal     = getattr(bg_spec,     "energy_cal", None)
        if not sample_ecal or not bg_ecal:
            raise ValueError(
                "check_bg_quality: peak_z_roi_keV_list requires energy_cal "
                "on both sample_spec and bg_spec"
            )
        for e_lo, e_hi in peak_z_roi_keV_list:
            B1 = _roi_sum_counts(sample_counts, sample_ecal, e_lo, e_hi)
            B2 = _roi_sum_counts(bg_counts,     bg_ecal,     e_lo, e_hi)
            z_res = bg_z_test(B1, B2)
            key = f"z_roi_{e_lo:.0f}_{e_hi:.0f}"
            gates[key] = {
                "value": z_res.z,
                "abs_z": z_res.abs_z,
                "tier": z_res.tier,
                "B1": z_res.B1,
                "B2": z_res.B2,
                "threshold_stable": Z_TIER_STABLE_MAX,
                "threshold_reject": Z_TIER_BORDERLINE_MAX,
                "passed": z_res.passed,
            }
            if z_res.note:
                notes.append(
                    f"ROI [{e_lo:.0f},{e_hi:.0f} keV] z-test: {z_res.note}"
                )

    return BgQualityReport(
        passed=all(g["passed"] for g in gates.values()),
        gates=gates,
        notes=notes,
    )


__all__ = [
    "BgControlResult",
    "validate_background",
    "validate_background_roi",
    "BgQualityReport",
    "check_bg_quality",
    "BgZTestResult",
    "bg_z_test",
    "bg_z_test_rates",
    "Z_TIER_STABLE_MAX",
    "Z_TIER_BORDERLINE_MAX",
]
