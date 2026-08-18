"""
F-283 (v1.17.13, T-003) — NaI pure-Gaussian peak-shape policy guard.

LSRM Algorithmic Foundations 2022 §8.4.2.1 explicitly states that NaI
peak shape (Gamma-1S 63×63 specifically) is **pure Gaussian** —
without a measurable low-energy exponential tail. This contradicts the
generic ORTEC / Genie HPGe-style peak_image model (Gauss + tail + step).

**Project state** (as of v1.17.12):
- `coupled_intensity_fit(..., use_peak_image=True, tail_param=0.7)`
  uses Gauss+tail+step. F-127 calibrated `T(E)` improves χ²/ν on
  Th-232 demo from 37.68 to ≤25, suggesting **measurable** tail on
  the demo fixtures (possibly artifact of background subtraction).

**T-003 policy decision** (this module):
- Provide opt-in **pure-Gaussian guard**: `force_pure_gaussian_for(detector_class)`
  returns `(tail_param, h_step) = (0.0, 0.0)` for NaI/CsI, preserving
  current values for other classes.
- Pipeline callers (staged_pipeline, demo scripts) can override the
  default `tail_param` / `h_step` via this helper to satisfy strict
  LSRM §8.4.2.1 compliance.
- Default `coupled_intensity_fit` behaviour is UNCHANGED (current
  F-120/F-127/F-133 retain calibrated tail+step).

Reference
---------
- LSRM Algorithmic Foundations 2022 §8.4.2.1 "Калибровка по форме линии"
  ("Для NaI(Tl) форма ФЭП — гауссиан")
- Project F-127 (calibrated T(E) on NaI 63×63) — opt-in performance
  improvement that nonetheless violates strict §8.4.2.1 reading
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple


# Class-set requiring pure-Gaussian per LSRM §8.4.2.1.
PURE_GAUSSIAN_DETECTOR_CLASSES = frozenset({"NaI", "CsI"})


@dataclass(frozen=True)
class PeakShapePolicy:
    """Bundle of peak-shape parameters per LSRM-strict policy."""
    tail_param: float        # T-tail parameter; 0 = pure Gaussian
    h_step: float            # Compton-step fraction; 0 = no step
    use_T_E_model: bool      # whether to apply F-127 T(E) calibration
    note: str = ""


def force_pure_gaussian_for(
    detector_class: str,
) -> Optional[PeakShapePolicy]:
    """Return a strict-LSRM pure-Gaussian policy for NaI/CsI; None otherwise.

    Callers use this when running staged_pipeline with strict §8.4.2.1
    compliance::

        policy = force_pure_gaussian_for("NaI")
        if policy:
            tail_param = policy.tail_param
            h_step = policy.h_step
        # else: use the function's default calibrated tail/step

    Parameters
    ----------
    detector_class : str
        Detector class string ("NaI", "HPGe", etc).

    Returns
    -------
    PeakShapePolicy if detector_class in {"NaI", "CsI"}; else None.
    """
    dc = str(detector_class).strip()
    if dc in PURE_GAUSSIAN_DETECTOR_CLASSES:
        return PeakShapePolicy(
            tail_param=0.0,
            h_step=0.0,
            use_T_E_model=False,
            note=(
                "F-283/T-003: strict LSRM §8.4.2.1 для NaI/CsI — pure Gaussian "
                "без tail. Disable F-127 T(E) model. Используйте, если "
                "приоритет — соответствие методике ЛСРМ, а не минимизация χ²/ν."
            ),
        )
    return None


def recommend_peak_shape(
    detector_class: str,
    *,
    strict_lsrm: bool = False,
) -> Tuple[float, float, bool]:
    """Convenience: return (tail_param, h_step, use_T_E_model) tuple.

    When `strict_lsrm=True` AND detector_class ∈ {NaI, CsI} —
    returns pure-Gaussian policy (0, 0, False).

    Otherwise returns project defaults: (T_TAIL_DEFAULT_NAI=0.7,
    H_STEP_DEFAULT_NAI=0.03, True) for NaI/CsI, и (0, 0, False) для всего
    остального.
    """
    if strict_lsrm:
        policy = force_pure_gaussian_for(detector_class)
        if policy is not None:
            return policy.tail_param, policy.h_step, policy.use_T_E_model
    # Project defaults
    dc = str(detector_class).strip()
    if dc in PURE_GAUSSIAN_DETECTOR_CLASSES:
        # NaI default with calibrated tail+step (F-120/F-127/F-133)
        return 0.7, 0.03, True
    # HPGe / LaBr / CeBr / other → no tail, no step (already sharp)
    return 0.0, 0.0, False


__all__ = [
    "PURE_GAUSSIAN_DETECTOR_CLASSES",
    "PeakShapePolicy",
    "force_pure_gaussian_for",
    "recommend_peak_shape",
]
