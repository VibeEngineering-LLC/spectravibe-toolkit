"""
Anchor-based energy calibration refit — Step 5β (F-87c / v1.15.0).

Per SKILL.md §5 "Bootstrap energy calibration (conditional)":
    Decision rule: skip this step if the file's stored calibration
    produces residuals < 0.3·FWHM at all confirmed library lines.

The express anchor heuristic (F-79 anchor-rank + F-80 express-patterns)
run at Step 5α gives us a set of AnchorMatch records — practical-
visibility nuclide lines confirmed against the found peaks. This
module turns those matches into a residual check and, if the stored
calibration disagrees, refits E(N) on the anchor channels.

The refit reuses :func:`gamma.calibration.energy_fit.polynomial_energy_fit`
with degree ≤ 4 (the SKILL.md cap). Start at degree 1; raise only if
needed. Co-60 partner-required-but-missing matches are excluded —
they are unsafe anchors. Single matches (n < 3) cannot give a useful
deg-1 fit and are rejected outright (`applied=False`).

Returns a tuple ``(new_energy_cal_or_None, diagnostic_dict)``:

    new_energy_cal : tuple | None
        Refitted low-to-high polynomial coefficients (same convention
        as ``Spectrum.energy_cal``). ``None`` when no refit was
        performed (either disagreement below threshold or too few
        anchors).
    diagnostic_dict : dict
        Always populated:
        - attempted: bool — whether the check ran end-to-end
        - applied: bool — whether a refit was actually returned
        - old_residual_max_keV: float | None
        - new_residual_max_keV: float | None
        - old_residual_max_fraction_of_fwhm: float | None
        - n_anchors_used: int
        - old_energy_cal: list | None
        - new_energy_cal: list | None
        - degree_used: int | None
        - reason: str
"""
from __future__ import annotations

from typing import Iterable, Optional, Tuple, List, Dict, Any

from gamma.calibration.energy_fit import (
    polynomial_energy_fit, MAX_POLYNOMIAL_DEGREE,
)


def _empty_diag(reason: str = "") -> Dict[str, Any]:
    return {
        "attempted": False,
        "applied": False,
        "old_residual_max_keV": None,
        "new_residual_max_keV": None,
        "old_residual_max_fraction_of_fwhm": None,
        "n_anchors_used": 0,
        "old_energy_cal": None,
        "new_energy_cal": None,
        "degree_used": None,
        "reason": reason,
    }


def recalibrate_energy_if_anchors_disagree(
    spec,
    anchor_matches: Iterable,
    *,
    threshold_fraction_of_fwhm: float = 0.3,
    fwhm_provider_keV,
    min_anchors: int = 3,
    max_degree: int = MAX_POLYNOMIAL_DEGREE,
) -> Tuple[Optional[Tuple[float, ...]], Dict[str, Any]]:
    """Check anchor residuals against stored E(N) and refit if needed.

    Parameters
    ----------
    spec : Spectrum
        Source of the stored ``energy_cal`` and ``channel_to_energy``.
    anchor_matches : iterable of AnchorMatch
        From :func:`gamma.identification.anchor_ranks.seed_calibration_anchors`.
    threshold_fraction_of_fwhm : float, default 0.3
        SKILL.md §5 threshold: residual exceeding this fraction of
        FWHM at any anchor triggers a refit. The same fraction is
        used as the target residual for the new fit.
    fwhm_provider_keV : callable
        ``f(E_keV) -> FWHM_keV``.
    min_anchors : int, default 3
        Minimum anchor count to attempt a refit. Below this the
        bootstrap is unsafe (degree-1 needs ≥2 points; we require ≥3
        so degree-1 has at least one residual to check).
    max_degree : int, default 4
        SKILL.md hard cap. ``polynomial_energy_fit`` will start at
        degree 1 and walk up only as needed.

    Returns
    -------
    (new_cal, diag) : tuple
        ``new_cal`` is the refitted coefficient tuple (low-to-high) or
        ``None`` when no refit is needed / possible. ``diag`` is
        always populated; see module docstring for keys.
    """
    diag = _empty_diag("init")
    diag["attempted"] = True

    # Stored calibration tuple (may be empty/identity)
    stored_cal = tuple(getattr(spec, "energy_cal", ()) or ())
    diag["old_energy_cal"] = list(stored_cal) if stored_cal else None

    # Filter usable anchors: nuclide non-empty, partner not missing,
    # finite delta. Exclude 511 (no nuclide assigned) and U-235/Ra-226
    # ambiguity by the anchor.nuclide check.
    usable: List[Tuple[int, float, float]] = []   # (ch, E_lib, delta)
    for am in anchor_matches:
        anchor = getattr(am, "anchor", None)
        if anchor is None or not anchor.nuclide:
            continue
        if getattr(am, "partner_required_but_missing", False):
            continue
        ch = int(getattr(am, "peak_channel"))
        E_lib = float(anchor.energy_keV)
        delta = abs(float(getattr(am, "delta_keV")))
        usable.append((ch, E_lib, delta))

    diag["n_anchors_used"] = len(usable)

    if not usable:
        diag["reason"] = "no usable anchors (all empty / partner-missing)"
        return None, diag

    # Compute max residual fraction of FWHM
    max_frac = 0.0
    max_abs = 0.0
    for ch, E_lib, delta in usable:
        f = max(1e-6, float(fwhm_provider_keV(E_lib)))
        frac = delta / f
        if frac > max_frac:
            max_frac = frac
        if delta > max_abs:
            max_abs = delta

    diag["old_residual_max_keV"] = max_abs
    diag["old_residual_max_fraction_of_fwhm"] = max_frac

    if max_frac < threshold_fraction_of_fwhm:
        diag["reason"] = (
            f"stored calibration within tolerance "
            f"({max_frac:.3f} < {threshold_fraction_of_fwhm})"
        )
        return None, diag

    if len(usable) < min_anchors:
        diag["reason"] = (
            f"insufficient anchors for safe refit "
            f"({len(usable)} < {min_anchors})"
        )
        return None, diag

    # Refit E(N) on the anchor channels.
    channels = [u[0] for u in usable]
    energies = [u[1] for u in usable]

    # Target the same threshold for the refit's residual goal.
    target_keV = threshold_fraction_of_fwhm * min(
        max(1e-6, fwhm_provider_keV(E)) for _, E, _ in usable
    )

    fit = polynomial_energy_fit(
        channels=channels,
        energies=energies,
        max_degree=max_degree,
        target_residual_keV=target_keV,
        min_degree=1,
    )

    # Compute new residual max in keV vs library energies
    predicted = fit.predict(channels)
    new_resids = [abs(float(p) - float(E)) for p, E in zip(predicted, energies)]
    new_max = max(new_resids) if new_resids else 0.0
    diag["new_residual_max_keV"] = new_max
    diag["degree_used"] = int(fit.degree)
    new_cal = tuple(fit.coefficients)
    diag["new_energy_cal"] = list(new_cal)

    # Only "apply" if the refit actually improved residuals
    if new_max < max_abs:
        diag["applied"] = True
        diag["reason"] = (
            f"refit improved residuals: "
            f"{max_abs:.2f} → {new_max:.2f} keV "
            f"(deg {fit.degree}, {len(usable)} anchors)"
        )
        return new_cal, diag
    else:
        diag["reason"] = (
            f"refit did not improve residuals "
            f"({max_abs:.2f} → {new_max:.2f} keV), keeping stored cal"
        )
        return None, diag


def should_auto_recalibrate(
    anchor_matches: Iterable,
    *,
    fwhm_provider_keV,
    drift_frac_threshold: float = 0.5,
    min_anchors: int = 3,
) -> Tuple[bool, Dict[str, Any]]:
    """F-453 (BUG-38 follow-up) — auto-trigger Step 5β без kwarg-opt-in.

    Standard `recalibrate_on_anchor_disagreement=False` default preserves
    the v1.14.0 contract, but on short NaI fixtures (AmTiCsEu, Cs-Co) the
    stored ADC→keV calibration can drift well past 0.3·FWHM on low-energy
    anchors with NO F-145 multiplet self-cal recovery (n_multiplets_seen=0
    when no Th/U chains are present). This helper returns True iff:

      • ≥ ``min_anchors`` usable anchors in ``anchor_matches`` AND
      • max |Δ|/FWHM over usable anchors > ``drift_frac_threshold``
        (вдвое выше стандартного 0.3·FWHM — "однозначный" drift, выше
        шума singleton-fitting).

    When True, caller proceeds with the standard
    :func:`recalibrate_energy_if_anchors_disagree` (0.3·FWHM target).
    """
    diag: Dict[str, Any] = {
        "fired": False,
        "max_frac_of_fwhm": 0.0,
        "n_usable_anchors": 0,
        "threshold": drift_frac_threshold,
        "reason": "",
    }
    usable_count = 0
    max_frac = 0.0
    for am in anchor_matches:
        anchor = getattr(am, "anchor", None)
        if anchor is None or not getattr(anchor, "nuclide", ""):
            continue
        if getattr(am, "partner_required_but_missing", False):
            continue
        E_lib = float(getattr(anchor, "energy_keV", 0.0) or 0.0)
        if E_lib <= 0:
            continue
        delta = abs(float(getattr(am, "delta_keV", 0.0) or 0.0))
        fwhm = max(1e-6, float(fwhm_provider_keV(E_lib)))
        frac = delta / fwhm
        if frac > max_frac:
            max_frac = frac
        usable_count += 1
    diag["max_frac_of_fwhm"] = max_frac
    diag["n_usable_anchors"] = usable_count
    if usable_count < min_anchors:
        diag["reason"] = (
            f"insufficient usable anchors "
            f"({usable_count} < {min_anchors})"
        )
        return False, diag
    if max_frac <= drift_frac_threshold:
        diag["reason"] = (
            f"max drift within tolerance "
            f"({max_frac:.3f} <= {drift_frac_threshold})"
        )
        return False, diag
    diag["fired"] = True
    diag["reason"] = (
        f"F-453 auto-trigger: max |Δ|/FWHM = {max_frac:.3f} > "
        f"{drift_frac_threshold} on {usable_count} anchors"
    )
    return True, diag


__all__ = [
    "recalibrate_energy_if_anchors_disagree",
    "should_auto_recalibrate",
]
