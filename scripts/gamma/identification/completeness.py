"""
Dose Contribution (DC) completeness metric (F-61 / v1.11.1).

Per Lsrm Algorithmic Foundations §14.2:

  DC_unidentified = D_unidentified / D_total · 100 %

where the dose-proxy D is computed as

  D = Σ (S_i · E_i) / ε(E_i)

over a set of peaks; S_i is the net peak area in counts and E_i is the
energy in keV. ε(E_i) is the photopeak efficiency at E_i.

Interpretation:
  • DC close to 0 %    — all dose-relevant peaks are identified;
                         identification is essentially complete.
  • DC > 10 %          — significant unidentified contribution; the
                         report should flag this and (per user
                         methodology) consider asking the user about
                         Stage 2/3 candidates.
  • DC > 30 %          — incomplete identification; analysis quality
                         is questionable.

Efficiency-free fallback:
  If no efficiency model is available (the common case for ad-hoc
  analysis of arbitrary `.spe` files), we use S_i · E_i directly as a
  dose proxy. The ratio DC = D_unident / D_total is still meaningful
  because efficiency cancels approximately when peaks are at similar
  energies; the metric is approximate when energies span the full
  detector range, but the gross "explained vs unexplained" partition
  remains informative.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Optional


@dataclass
class CompletenessResult:
    """Output of `compute_completeness`."""
    dose_identified: float          # Σ over confirmed+tentative nuclide lines
    dose_unidentified: float        # Σ over residual true_unmatched + un-explained
    dose_total: float
    dc_percent: float               # 100 · dose_unidentified / dose_total
    n_identified_lines: int
    n_unidentified_peaks: int
    flag: str                       # "complete" / "marginal" / "incomplete"
    note: str = ""

    @property
    def is_complete(self) -> bool:
        return self.dc_percent < 10.0


# ──────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────

def _line_dose(
    energy_keV: float,
    area: float,
    efficiency_at: Optional[Callable[[float], float]],
) -> float:
    """Dose proxy for a single matched line."""
    if efficiency_at is not None:
        try:
            e = efficiency_at(energy_keV)
            if e and e > 0:
                return float(area) * float(energy_keV) / float(e)
        except Exception:
            pass
    # Fallback: just area·E
    return float(area) * float(energy_keV)


def compute_completeness(
    detected_nuclides: List,
    unmatched_peaks_E_area: List,
    *,
    efficiency_at: Optional[Callable[[float], float]] = None,
) -> CompletenessResult:
    """
    Compute the Dose Contribution metric.

    Args:
        detected_nuclides: list of NuclideIdentification — only their
            `matched_lines` are summed.
        unmatched_peaks_E_area: list of (energy_keV, area) tuples for
            peaks that ended up as `true_unmatched` after the residual
            classifier (F-74). Pre-classified explanations (XRF, sum,
            escape, chain-secondary) are EXCLUDED from "unidentified" —
            they count as identified non-FEP features.
        efficiency_at: optional callable(E_keV) → efficiency. When
            absent, the dose proxy is area·E (no 1/ε factor).

    Returns:
        CompletenessResult.
    """
    dose_id = 0.0
    n_lines = 0
    for nid in detected_nuclides:
        for m in getattr(nid, "matched_lines", ()):
            area = getattr(m, "peak_area", None)
            E = getattr(m, "library_E_keV", None)
            if area is not None and E is not None and area > 0:
                dose_id += _line_dose(float(E), float(area), efficiency_at)
                n_lines += 1

    dose_unid = 0.0
    n_unid = 0
    n_unid_no_area = 0
    for E_a in unmatched_peaks_E_area:
        if len(E_a) != 2:
            continue
        E, area = E_a
        if E is None:
            continue
        # Count the residual even if area is unknown / zero — it still
        # represents an unidentified peak that the operator should review.
        n_unid += 1
        if area is None or area <= 0:
            n_unid_no_area += 1
            continue
        dose_unid += _line_dose(float(E), float(area), efficiency_at)

    dose_total = dose_id + dose_unid
    if dose_total <= 0:
        return CompletenessResult(
            dose_identified=0.0, dose_unidentified=0.0, dose_total=0.0,
            dc_percent=0.0, n_identified_lines=n_lines,
            n_unidentified_peaks=n_unid,
            flag="n/a",
            note="Невозможно вычислить DC — нет площадей пиков.",
        )
    dc_pct = 100.0 * dose_unid / dose_total

    if n_unid_no_area > 0:
        # Promote to marginal if many residuals lack area data — the metric
        # underestimates true DC and the analyst should be aware.
        if n_unid_no_area > n_unid // 2:
            note_suffix = (f" (предупреждение: {n_unid_no_area}/{n_unid} "
                           "остаточных пиков без оценки площади — DC занижен).")
        else:
            note_suffix = ""
    else:
        note_suffix = ""

    if dc_pct < 10.0 and n_unid == 0:
        flag = "complete"
        note = (f"DC = {dc_pct:.1f} % — идентификация полная "
                f"({n_lines} линий идентифицировано, 0 остаточных).")
    elif dc_pct < 10.0:
        flag = "complete"
        note = (f"DC = {dc_pct:.1f} % — идентификация в основном полная "
                f"({n_lines} линий идентифицировано, {n_unid} остаточных пиков "
                f"с малой долей дозы).{note_suffix}")
    elif dc_pct < 30.0:
        flag = "marginal"
        note = (f"DC = {dc_pct:.1f} % — заметная неидентифицированная "
                "доля дозы; рассмотреть Stage 2.")
    else:
        flag = "incomplete"
        note = (f"DC = {dc_pct:.1f} % — идентификация неполная; "
                "требуется Stage 2/3 или ручной разбор.")

    return CompletenessResult(
        dose_identified=dose_id,
        dose_unidentified=dose_unid,
        dose_total=dose_total,
        dc_percent=dc_pct,
        n_identified_lines=n_lines,
        n_unidentified_peaks=n_unid,
        flag=flag,
        note=note,
    )


__all__ = [
    "CompletenessResult",
    "compute_completeness",
]
