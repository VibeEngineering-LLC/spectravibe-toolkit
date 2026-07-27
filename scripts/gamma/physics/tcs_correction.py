"""True Coincidence Summing (TCS / TSC) corrections — physics module skeleton.

This module is **Phase 5 of v1.19.0** TCS integration. It exposes a callable
correction API but is **NOT wired into the activity-computation pipeline** —
pipeline integration is deferred to v1.19.1 per
``audit/_drafts/PLAN_v1_18_32_to_v1_19_0_TCS_INTEGRATION.md`` §9 (Risk R5.1).

Public surface (5 functions + 1 inner physics helper)::

    fc_lookup(...)                       # detector-class-gated table lookup
    apply_tcs_correction(...)            # A_corr = A_obs × fc + GUM σ
    apply_tcs_correction_gated(...)      # top-level wrapper enforcing D1/D3
    eta_eff_marinelli_scaling(...)       # STUB → 1.0 + DeprecationWarning
    angular_correlation_factor(...)      # STUB → 1.0 + DeprecationWarning
    _compute_tscf(...)                   # inner: summing-out + summing-in

Module-level constants::

    ALLOWED_DETECTOR_CLASSES   # full enum recognized by the skill
    SUPPORTED_DETECTOR_CLASSES # subset with calibrated TSCF data in v1.19.0
    MIN_DETECTOR_CLASS_CI      # Step-4 detector-class CI threshold (D1)

Detector-class gating (D1 / D3 — PLAN §5)
-----------------------------------------

* **D3 fail-closed**: ``detector_class`` not in
  ``ALLOWED_DETECTOR_CLASSES`` → ``ValueError``. Garbage input is a
  programming bug; we refuse to silently produce a report.
* **D1 fail-soft**: low confidence in detector-class identification
  (``detector_class_CI < MIN_DETECTOR_CLASS_CI``) → no correction
  applied, ``tcs_status = "low-detector-CI"``.
* **D3 fail-soft**: known but not in ``SUPPORTED_DETECTOR_CLASSES``
  (e.g. NaI-Tl, LaBr3-Ce) → no correction applied,
  ``tcs_status = "out-of-scope-for-detector"``.

References
----------
* Sima O. (2018) ICRM GSWG, Paris, June 2018 — coincidence summing
  corrections in γ-spectrometry (65-slide deck). Authoritative source
  for the closed-form formulas implemented in ``_compute_tscf``.
  See ``references/12_tcs_overview.md`` for citation breakdown.
* Andreev D.S. et al. (1972) Instrum. Exp. Tech. 15:1358 — original
  recursion for summing-out.
* Semkow T.M. et al. (1990) NIM A290:437 — matrix formalism (basis
  of EFFTRAN), separates summing-out and summing-in.
* Giubrone G. et al. (2016) JER 158-159:114 — TSCF values for
  HPGe-coaxial in petri/Marinelli geometries (data source for
  ``data/tsc_lookup.json``).
* Ordonez-Ródenas J. et al. (2019) RPC 155:244 — natural-chain TSCF
  extension (Bi-214, Bi-212, Tl-208) on the same ORTEC GMX40
  detector.
* PLAN: ``audit/_drafts/PLAN_v1_18_32_to_v1_19_0_TCS_INTEGRATION.md``
  §5 Phase 5 — module spec, Refinements A/B/C, gating decisions D1-D3.
"""
from __future__ import annotations

import json
import math
import warnings
from pathlib import Path
from typing import Any, Callable, Iterable, Optional, Tuple, Union

__all__ = [
    "ALLOWED_DETECTOR_CLASSES",
    "SUPPORTED_DETECTOR_CLASSES",
    "MIN_DETECTOR_CLASS_CI",
    "fc_lookup",
    "apply_tcs_correction",
    "apply_tcs_correction_gated",
    "eta_eff_marinelli_scaling",
    "angular_correlation_factor",
]


# ---------------------------------------------------------------------------
# Module-level constants — DO NOT rename; tests reference these by name.
# ---------------------------------------------------------------------------

#: Full set of detector classes the skill knows about (Step 4 enum).
#: A value outside this set is treated as a programming error (D3
#: fail-closed in ``apply_tcs_correction_gated``).
ALLOWED_DETECTOR_CLASSES = frozenset({
    "HPGe-coaxial",
    "HPGe-BEGE",
    "NaI-Tl",
    "LaBr3-Ce",
    "CeBr3",
    "CdZnTe",
})

#: Subset of ``ALLOWED_DETECTOR_CLASSES`` for which calibrated TSCF
#: data is shipped in v1.19.0. Future detector-class expansion requires
#: BOTH adding entries to ``data/tsc_lookup.json`` AND extending this
#: set. Keeping this conservative is intentional — uncalibrated values
#: would silently bias activity by 10-30% on coincidence cascades
#: (Co-60, Eu-152, Cs-134, natural-chain progeny).
SUPPORTED_DETECTOR_CLASSES = frozenset({"HPGe-coaxial"})

#: Minimum acceptable Step-4 detector-class CI before applying any
#: TSCF correction. Below this threshold the wrapper falls back to
#: uncorrected output with a ``"low-detector-CI"`` status (D1).
MIN_DETECTOR_CLASS_CI: float = 0.7


# ---------------------------------------------------------------------------
# Lookup table — loaded lazily from data/tsc_lookup.json.
# ---------------------------------------------------------------------------

# Repo root: <root>/scripts/gamma/physics/tcs_correction.py → parents[3].
_REPO_ROOT = Path(__file__).resolve().parents[3]
_LOOKUP_PATH = _REPO_ROOT / "data" / "tsc_lookup.json"

_LOOKUP_CACHE: Optional[list] = None


def _load_lookup() -> list:
    """Return list of TSCF entries from ``data/tsc_lookup.json``.

    Cached after first call. Per Phase 4 schema (PLAN §3), each entry
    has at least: ``nuclide``, ``energy_keV``, ``geometry``,
    ``detector_class``, ``tscf``, ``tscf_uncertainty``.
    """
    global _LOOKUP_CACHE
    if _LOOKUP_CACHE is None:
        if not _LOOKUP_PATH.exists():
            _LOOKUP_CACHE = []
        else:
            payload = json.loads(_LOOKUP_PATH.read_text(encoding="utf-8"))
            _LOOKUP_CACHE = list(payload.get("entries", []))
    return _LOOKUP_CACHE


def reset_lookup_cache() -> None:
    """Clear the in-memory TSCF lookup cache.

    REL-02 (AUDIT_v2 §6 P0-1): called by the autouse cache-reset fixture
    in ``tests/conftest.py`` to guarantee per-test isolation under
    ``pytest -n auto``. Production code should not need this; ``_load_lookup``
    is idempotent and the JSON payload is small (~few hundred entries).
    """
    global _LOOKUP_CACHE
    _LOOKUP_CACHE = None


def _coerce_geometry(geometry_class: Any) -> str:
    """Return the string form of a GeometryClass enum or pass through str."""
    if hasattr(geometry_class, "value"):
        return str(geometry_class.value)
    return str(geometry_class)


# ---------------------------------------------------------------------------
# fc_lookup — detector-class-gated table lookup.
# ---------------------------------------------------------------------------

def fc_lookup(
    nuclide: str,
    line_E_keV: float,
    geometry_class: Any,
    detector_class: str,
    tolerance_keV: float = 0.5,
) -> Optional[Tuple[float, float]]:
    """Lookup ``(fc, sigma_fc)`` for one line in ``data/tsc_lookup.json``.

    The ``detector_class`` filter is applied **before** the energy
    tolerance match (per PLAN §5 step 1 fc_lookup description). This
    means a call with ``detector_class="NaI-Tl"`` returns ``None``
    even if the same ``(nuclide, line_E_keV, geometry_class)`` triple
    has a matching HPGe-coaxial entry — the table is class-specific,
    cross-class extrapolation is forbidden (R5.4).

    Per Giubrone 2016 §5 finding that TSCFs are weakly dependent on
    matrix absorber within a single geometry-class, we lookup by
    geometry alone (matrix density is recorded per-entry as metadata
    but not used as a primary filter in v1.19.0).

    Parameters
    ----------
    nuclide : str
        Canonical nuclide identifier, e.g. ``"Co-60"``.
    line_E_keV : float
        Gamma-line energy in keV.
    geometry_class : str | GeometryClass
        Geometry identifier. Either a ``GeometryClass`` enum member
        or its string value (e.g. ``"petri_100ml_water"``).
    detector_class : str
        Detector-class string. MUST be one of
        ``ALLOWED_DETECTOR_CLASSES``; mismatch causes a ``None``
        return (the gated wrapper above translates this to a
        ``"no-cascade-data"`` status). The caller is responsible for
        ensuring ``detector_class`` is a known enum member;
        ``apply_tcs_correction_gated`` does this check upstream.
    tolerance_keV : float, optional
        Energy match tolerance (half-window). Default 0.5 keV.

    Returns
    -------
    Optional[Tuple[float, float]]
        ``(fc, sigma_fc)`` on hit, ``None`` on miss.
    """
    entries = _load_lookup()
    if not entries:
        return None
    geom = _coerce_geometry(geometry_class)

    # Filter: detector_class AND nuclide AND geometry, BEFORE E-tolerance.
    candidates = [
        e
        for e in entries
        if e.get("detector_class") == detector_class
        and e.get("nuclide") == nuclide
        and e.get("geometry") == geom
    ]
    if not candidates:
        return None

    # Closest-E match within tolerance.
    best: Optional[dict] = None
    best_dE = math.inf
    for e in candidates:
        dE = abs(float(e["energy_keV"]) - float(line_E_keV))
        if dE <= tolerance_keV and dE < best_dE:
            best = e
            best_dE = dE
    if best is None:
        return None

    return (float(best["tscf"]), float(best["tscf_uncertainty"]))


# ---------------------------------------------------------------------------
# apply_tcs_correction — A_corr = A_obs × fc + GUM σ propagation.
# ---------------------------------------------------------------------------

def apply_tcs_correction(
    A_observed: float,
    A_uncertainty: float,
    fc: float,
    sigma_fc: float,
) -> Tuple[float, float]:
    """Apply scalar TCS correction with GUM uncertainty propagation.

    ``FC = R_observed / R_no-coinc`` per Sima 2018 §1.1 definition, so
    the **true** activity recovered from observed counts is
    ``A_true = A_obs × FC``.

    Uncertainty propagation (independent ``A_obs`` and ``fc``):

    .. math::

        \\sigma_{A_\\text{corr}}^{2}
            = (fc \\cdot \\sigma_{A_\\text{obs}})^2
              + (A_\\text{obs} \\cdot \\sigma_{fc})^2

    See ``references/12_tcs_overview.md`` §9 (uncertainty propagation)
    for the derivation and references to Sima 2018 §7.1-§7.2.

    Parameters
    ----------
    A_observed : float
        Observed activity (any consistent unit, e.g. Bq).
    A_uncertainty : float
        ``1σ`` uncertainty on ``A_observed``.
    fc : float
        Correction factor from ``fc_lookup`` (TSCF, ≥ 1 for cascade
        loss-dominated lines).
    sigma_fc : float
        ``1σ`` uncertainty on ``fc``.

    Returns
    -------
    Tuple[float, float]
        ``(A_corrected, sigma_corrected)``.
    """
    A_corr = A_observed * fc
    sigma_corr = math.sqrt(
        (fc * A_uncertainty) ** 2 + (A_observed * sigma_fc) ** 2
    )
    return (A_corr, sigma_corr)


# ---------------------------------------------------------------------------
# apply_tcs_correction_gated — top-level wrapper enforcing D1/D3 gating.
# ---------------------------------------------------------------------------

def apply_tcs_correction_gated(
    A_observed: float,
    A_uncertainty: float,
    nuclide: str,
    line_E_keV: float,
    geometry_class: Any,
    detector_class: Optional[str],
    detector_class_CI: float,
    report_diagnostics: dict,
) -> Tuple[float, float, str]:
    """Top-level TCS-correction wrapper enforcing detector-class gating.

    Behavior per PLAN §5 Phase 5 D1/D3 lockdown:

    * **D3 fail-closed**: ``detector_class`` not in
      ``ALLOWED_DETECTOR_CLASSES`` (typo, ``None``, garbage) →
      ``ValueError``. No silent recovery.
    * **D1 fail-soft**: ``detector_class_CI < MIN_DETECTOR_CLASS_CI``
      → ``tcs_status = "low-detector-CI"``, return uncorrected
      ``(A_observed, A_uncertainty)``, warning appended.
    * **D3 fail-soft**: ``detector_class`` in ``ALLOWED`` but not in
      ``SUPPORTED_DETECTOR_CLASSES`` (e.g. NaI-Tl, LaBr3-Ce,
      HPGe-BEGE) → ``tcs_status = "out-of-scope-for-detector"``,
      return uncorrected, warning appended naming the specific class.
    * No matching lookup entry → ``tcs_status = "no-cascade-data"``.
    * Otherwise → ``tcs_status = "applied"`` and corrected values
      returned.

    Every non-``"applied"`` status appends to
    ``report_diagnostics["tcs_warnings"]`` with a structured entry
    ``{nuclide, line_E_keV, detector_class, status, message, ...}``.
    These are rendered in v1.19.1+ in both the Step-11 diagnostic
    block AND the primary-FEP table ``tcs_status`` column (D2 parity).

    .. note::
        v1.19.0 ships this wrapper **callable but not wired** —
        compute_activity / json_report integration is v1.19.1
        (PLAN Risk R5.1 explicitly forbids leaking into pipeline
        here).
    """
    # D3 fail-closed: ValueError for unknown enum members. Crucially,
    # ``None`` and typos like "NaI" land here. We refuse to silently
    # produce a report on garbage input.
    if detector_class not in ALLOWED_DETECTOR_CLASSES:
        raise ValueError(
            f"Unknown detector_class={detector_class!r}; allowed: "
            f"{sorted(ALLOWED_DETECTOR_CLASSES)}"
        )

    # D1 fail-soft: detector class identified, but Step-4 CI too low.
    if detector_class_CI < MIN_DETECTOR_CLASS_CI:
        tcs_status = "low-detector-CI"
        message = (
            f"TSC correction not applied: detector class identification "
            f"below confidence threshold (CI={detector_class_CI:.2f} < "
            f"{MIN_DETECTOR_CLASS_CI}); cannot select TSCF table."
        )
        report_diagnostics.setdefault("tcs_warnings", []).append({
            "nuclide": nuclide,
            "line_E_keV": line_E_keV,
            "detector_class": detector_class,
            "detector_class_CI": detector_class_CI,
            "status": tcs_status,
            "message": message,
        })
        return (A_observed, A_uncertainty, tcs_status)

    # D3 fail-soft: known enum but no calibrated TSCF data in v1.19.0.
    if detector_class not in SUPPORTED_DETECTOR_CLASSES:
        tcs_status = "out-of-scope-for-detector"
        message = (
            f"TSC correction not applied for {nuclide}@{line_E_keV} keV: "
            f"no calibrated TSCF data for detector_class={detector_class}. "
            f"Reported activity may be biased by 10-30% on coincidence "
            f"cascades (Co-60, Eu-152, Cs-134, natural-chain progeny). "
            f"Calibration deferred to v1.20.x."
        )
        report_diagnostics.setdefault("tcs_warnings", []).append({
            "nuclide": nuclide,
            "line_E_keV": line_E_keV,
            "detector_class": detector_class,
            "status": tcs_status,
            "message": message,
        })
        return (A_observed, A_uncertainty, tcs_status)

    # Supported detector class with adequate CI → attempt lookup.
    result = fc_lookup(nuclide, line_E_keV, geometry_class, detector_class)
    if result is None:
        tcs_status = "no-cascade-data"
        message = (
            f"No TSCF lookup entry for {nuclide}@{line_E_keV} keV in "
            f"{_coerce_geometry(geometry_class)}."
        )
        report_diagnostics.setdefault("tcs_warnings", []).append({
            "nuclide": nuclide,
            "line_E_keV": line_E_keV,
            "detector_class": detector_class,
            "status": tcs_status,
            "message": message,
        })
        return (A_observed, A_uncertainty, tcs_status)

    fc, sigma_fc = result
    A_corr, sigma_corr = apply_tcs_correction(A_observed, A_uncertainty, fc, sigma_fc)
    return (A_corr, sigma_corr, "applied")


# ---------------------------------------------------------------------------
# eta_eff_marinelli_scaling — STUB (v1.19.0).
# ---------------------------------------------------------------------------

def eta_eff_marinelli_scaling(
    E_i_keV: float,
    E_j_keV: float,
    geometry: Any,
) -> float:
    """STUB: η_eff scaling for volume sources — returns ``1.0`` in v1.19.0.

    Full implementation requires the per-detector LS-curve ``l(E)`` of
    Vidmar-Korun NIMA 556 (2006) 543 + the per-geometry numerical
    table of Arnold-Sima JRNC 248 (2001) 365 (Marinelli +16% to +44%
    enhancement of η_eff over plain η). Neither data set is shipped
    in v1.19.0 — see ``references/12_tcs_overview.md`` Gaps.

    Emits ``DeprecationWarning`` to prevent silent reliance on the
    placeholder value.
    """
    warnings.warn(
        "eta_eff scaling not implemented; requires LS-curve l(E) per "
        "detector — see references/12_tcs_overview.md",
        DeprecationWarning,
        stacklevel=2,
    )
    return 1.0


# ---------------------------------------------------------------------------
# angular_correlation_factor — STUB (v1.19.0).
# ---------------------------------------------------------------------------

def angular_correlation_factor(
    geometry: Any,
    delta_omega_over_4pi: float,
) -> float:
    """STUB: W(θ) angular correlation factor — returns ``1.0`` in v1.19.0.

    Full implementation requires the per-cascade a₂/a₄ coefficients
    from NuDat 3.0 level schemes + the geometry-dependent
    averaging-over-solid-angle integral of Sima ARI 47 (1996) 919 /
    Sima 2018 §5 page 30 (~+1% to +10% depending on ΔΩ/4π).

    Deferred to v1.20.x; the slot exists in v1.19.0 only to lock the
    API surface so v1.20.x can wire it without breaking callers
    (Refinement C, PLAN §5 step 1).

    Emits ``DeprecationWarning`` to prevent silent reliance on the
    placeholder value.
    """
    warnings.warn(
        "angular correlation factor not implemented; requires NuDat "
        "a2/a4 coefficients + per-geometry solid-angle integral — "
        "see references/12_tcs_overview.md",
        DeprecationWarning,
        stacklevel=2,
    )
    return 1.0


# ---------------------------------------------------------------------------
# _compute_tscf — inner physics helper (summing-out + summing-in).
# ---------------------------------------------------------------------------

def _compute_tscf(
    line: Any,
    detector_response: Any,
    angular_correlation_provider: Optional[Callable[[Any, Any], float]] = None,
) -> dict:
    """Inner closed-form TSCF for one gamma line, exposing summing-out
    and summing-in as **separate multiplicative components**.

    Andreev/Semkow formalism (Andreev IET 15 (1972) 1358; Semkow NIM
    A290 (1990) 437; Sima 2018 §1.2-§1.3 slides 23-34, Refinement A
    in PLAN §5):

    **Summing-out** (FEP loss at energy ``E``)::

        TSCF_out(E) = 1 / (1 − Σ_j ε_T(E_j) · P_cascade(E → E_j) · W(θ))

    Each coincident gamma ``j`` with total efficiency ``ε_T(E_j)`` and
    cascade probability ``P_cascade`` pulls events out of the FEP at
    ``E`` into the sum bin at ``E + E_j``. ``TSCF_out ≥ 1``;
    multiplicative correction makes recovered activity larger.

    **Summing-in** (sum-peak creation at ``E_sum = E_i + E_j``)::

        TSCF_in(E_sum) = 1 + Σ ε_pp(E_i) · ε_pp(E_j) · P_cascade · W(θ)
                            / ε_pp(E_sum)

    Coincident events ADD to a sum-peak that doesn't correspond to
    any real cascade transition. This is the physically opposite
    contribution to TSCF_out (sum-peak gain rather than FEP loss),
    and the two terms scale differently across detector classes
    (R5.4 NaI case: TSCF_out grows 2-4× while TSCF_in shrinks 2-3×).
    Hence keeping them as distinct keys in the return dict.

    The final correction combines both effects multiplicatively when
    both channels are active::

        TSCF_total = TSCF_out × TSCF_in

    Parameters
    ----------
    line : object
        Must expose:

        * ``E_keV`` (float) — line energy
        * ``coincident_partners`` (iterable) — each with ``E_keV``,
          ``P_cascade``, and ``level_j`` for the angular-correlation
          provider
        * ``is_sum_peak_candidate`` (bool)
        * ``contributing_pairs`` (iterable of ``(i, j)`` pairs where
          each member exposes ``E_keV``, ``P_cascade``, ``level``)
          — only consulted if ``is_sum_peak_candidate`` is True
    detector_response : object
        Exposes ``eps_T(E)`` (total efficiency) and ``eps_pp(E)``
        (FEP / "photopeak" efficiency), both callable with a single
        keV argument.
    angular_correlation_provider : callable | None
        If provided, called as ``W = provider(level_i, level_j)`` for
        each summed pair. ``None`` (default) means **isotropic
        emission**, ``W ≡ 1`` — the v1.19.0 assumption locked by
        Refinement C.

    Returns
    -------
    dict
        ``{"tscf_out": float, "tscf_in": float, "tscf_total": float,
        "sigma": float}``.

        The ``sigma`` field is a placeholder for full per-line
        uncertainty propagation (v1.20.x — see references doc §9);
        in v1.19.0 it is computed as the relative-uncertainty
        envelope of ``tscf_out`` only (assuming 2% on each ε_T term,
        independent — same scale as Sima-Lépy 2016 typical ~few %).
        Wired-pipeline integration will replace this with the
        per-entry ``tscf_uncertainty`` from ``data/tsc_lookup.json``.
    """
    # --- Summing-OUT (FEP loss at line.E_keV).
    loss = 0.0
    for partner in getattr(line, "coincident_partners", ()) or ():
        if angular_correlation_provider is not None:
            W = float(
                angular_correlation_provider(
                    getattr(line, "level_i", None),
                    getattr(partner, "level_j", None),
                )
            )
        else:
            W = 1.0
        eps_T_j = float(detector_response.eps_T(partner.E_keV))
        loss += eps_T_j * float(partner.P_cascade) * W

    if loss >= 1.0:
        # Pathological input — would invert the correction. Cap below
        # 1.0 to keep the closed-form expression well-defined. In a
        # realistic detector ε_T stays well below the cascade-resummed
        # bound (1 − Σ ε_T < 1 always for HPGe-coaxial at expected
        # cascade multiplicities).
        loss = 1.0 - 1e-12
    tscf_out = 1.0 / (1.0 - loss)

    # --- Summing-IN (sum-peak creation at line.E_keV if it IS a sum peak).
    gain = 0.0
    if getattr(line, "is_sum_peak_candidate", False):
        eps_pp_sum = float(detector_response.eps_pp(line.E_keV))
        if eps_pp_sum > 0.0:
            for pair in getattr(line, "contributing_pairs", ()) or ():
                i_member, j_member = pair
                if angular_correlation_provider is not None:
                    W = float(
                        angular_correlation_provider(
                            getattr(i_member, "level", None),
                            getattr(j_member, "level", None),
                        )
                    )
                else:
                    W = 1.0
                eps_pp_i = float(detector_response.eps_pp(i_member.E_keV))
                eps_pp_j = float(detector_response.eps_pp(j_member.E_keV))
                P_cascade = float(
                    getattr(i_member, "P_cascade", 0.0)
                )
                gain += eps_pp_i * eps_pp_j * P_cascade * W / eps_pp_sum
    tscf_in = 1.0 + gain

    tscf_total = tscf_out * tscf_in

    # Placeholder σ — see docstring. Uses 2% relative-uncertainty
    # propagation on the loss term (independent ε_T's). Will be
    # replaced by per-entry tscf_uncertainty from the lookup table
    # in v1.19.1 pipeline wiring.
    rel_sigma_out = (
        0.02 * math.sqrt(max(loss, 0.0)) / max(1.0 - loss, 1e-12)
        if loss > 0
        else 0.0
    )
    sigma = abs(tscf_total) * rel_sigma_out

    return {
        "tscf_out": float(tscf_out),
        "tscf_in": float(tscf_in),
        "tscf_total": float(tscf_total),
        "sigma": float(sigma),
    }
