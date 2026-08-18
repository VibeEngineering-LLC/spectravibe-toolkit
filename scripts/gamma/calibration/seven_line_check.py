"""
Seven-line ЕРН calibration verification (F-81 / v1.11.1).

Implementation of the Lsrm methodology described in
NOTES_v1.7_methodology.md §9 — "NaI калибровка на ЕРН/ФОН-спектрах":

| E (keV) | Состав                                          | Тип            |
|---------|-------------------------------------------------|----------------|
| 240     | Pb-212 (238.63) + Pb-214 (241.98) superposition | always present |
| 351.93  | Pb-214 (Ra-226 chain)                           | clean          |
| 511     | Tl-208 (510.77) + annihilation                  | superposition  |
| 1120.29 | Bi-214 (Ra-226 chain)                           | clean          |
| 1460.82 | K-40                                            | clean          |
| 1764.49 | Bi-214 (Ra-226 chain)                           | clean          |
| 2614.51 | Tl-208 (Th-232 chain)                           | clean          |

The check:
  1. For each of the 7 reference lines, look for a Mariscotti peak
     within `match_window` of the expected energy.
  2. Record the residual (E_observed − E_expected) per line.
  3. Quality metric: number of lines present out of 7. Per methodology,
     ≥4-5 of 7 is sufficient for reliable calibration.
  4. Compute the mean and max absolute residual; flag calibration as
     "ok" if max residual < 0.3·FWHM at the corresponding line, "drift"
     if > 0.3·FWHM, "broken" if many lines missing or max residual is
     comparable to FWHM.

Applies to: ANY spectrum, but is primarily a **background-mode** check
per user methodology (15.11.2025): "для поиска пиков в фоне следуем
методике набор 7 реперных линий". Sample spectra typically have their
own anchor lines (Tl-208 2615, K-40 1461) which overlap with the
seven; the check still works.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, List, Optional, Tuple


# Per methodology §9 — the canonical 7 anchor lines.
SEVEN_LINES: Tuple[Tuple[float, str, str], ...] = (
    (240.00,  "Pb-212 (238.63) + Pb-214 (241.98)", "superposition"),
    (351.93,  "Pb-214",                            "clean"),
    (511.00,  "Tl-208 (510.77) + annihilation",    "superposition"),
    (1120.29, "Bi-214",                            "clean"),
    (1460.82, "K-40",                              "clean"),
    (1764.49, "Bi-214",                            "clean"),
    (2614.51, "Tl-208",                            "clean"),
)


@dataclass
class LineCheckResult:
    expected_keV: float
    description: str
    line_type: str
    found: bool
    observed_keV: Optional[float] = None
    residual_keV: Optional[float] = None
    sigma: Optional[float] = None
    fwhm_at_keV: Optional[float] = None
    residual_fwhm_fraction: Optional[float] = None  # |residual| / FWHM


@dataclass
class SevenLineCheck:
    """Outcome of the 7-line ЕРН calibration check."""
    per_line: List[LineCheckResult]
    lines_present: int                  # how many of 7 were found
    lines_total: int = 7
    max_residual_keV: Optional[float] = None
    mean_residual_keV: Optional[float] = None
    max_residual_fwhm_fraction: Optional[float] = None
    quality: str = "n/a"               # "ok" / "drift" / "broken" / "n/a"
    quality_note: str = ""

    @property
    def is_reliable(self) -> bool:
        """Per methodology: ≥4 of 7 lines present = reliable."""
        return self.lines_present >= 4


# ──────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────

def run_seven_line_check(
    found_peaks,
    spec,
    *,
    fwhm_provider_keV: Callable[[float], float],
    window_fwhm_multiple: float = 0.5,
    drift_threshold_fwhm_fraction: float = 0.3,
) -> SevenLineCheck:
    """
    Run the 7-line calibration check.

    Args:
        found_peaks: peaks from Mariscotti
        spec: Spectrum (used for channel-to-energy)
        fwhm_provider_keV: callable(E) -> FWHM in keV
        window_fwhm_multiple: half-width of the line-match window as a
            multiple of FWHM at the target energy. Default 0.5.
        drift_threshold_fwhm_fraction: a per-line residual exceeding
            this fraction of FWHM is flagged as drifted. Default 0.3
            per the methodology rule "residuals < 0.3·FWHM at anchor
            lines" for accepting stored calibration.

    Returns:
        SevenLineCheck with per-line results and an overall quality
        verdict.
    """
    peak_E = [(p, spec.channel_to_energy(p.channel)) for p in found_peaks]
    per_line: List[LineCheckResult] = []
    residuals: List[float] = []
    residual_fractions: List[float] = []

    for E_expect, desc, ltype in SEVEN_LINES:
        fwhm = max(2.0, fwhm_provider_keV(E_expect))
        tol = fwhm * window_fwhm_multiple
        # Find the closest Mariscotti peak inside the tolerance window
        best: Optional[Tuple[float, float, float]] = None  # (E_obs, |dE|, sigma)
        for p, e in peak_E:
            d = abs(e - E_expect)
            if d <= tol and (best is None or d < best[1]):
                best = (e, d, float(p.significance))
        if best is None:
            per_line.append(LineCheckResult(
                expected_keV=E_expect, description=desc, line_type=ltype,
                found=False, fwhm_at_keV=fwhm,
            ))
        else:
            E_obs, d, sig = best
            frac = d / fwhm if fwhm > 0 else None
            per_line.append(LineCheckResult(
                expected_keV=E_expect, description=desc, line_type=ltype,
                found=True,
                observed_keV=E_obs,
                residual_keV=E_obs - E_expect,
                sigma=sig,
                fwhm_at_keV=fwhm,
                residual_fwhm_fraction=frac,
            ))
            residuals.append(abs(E_obs - E_expect))
            if frac is not None:
                residual_fractions.append(frac)

    lines_present = sum(1 for r in per_line if r.found)

    max_res = max(residuals) if residuals else None
    mean_res = sum(residuals) / len(residuals) if residuals else None
    max_frac = max(residual_fractions) if residual_fractions else None

    # Quality verdict
    if lines_present == 0:
        quality = "broken"
        note = ("Ни одна из 7 ЕРН-реперных линий не найдена — "
                "энергетическая калибровка не верифицируется или ЕРН-фона нет.")
    elif lines_present < 4:
        quality = "broken"
        note = (f"Только {lines_present}/7 ЕРН-линий найдено — "
                f"для надёжной калибровки нужно ≥4. Возможен сильный дрейф "
                f"или искусственно-чистый спектр.")
    elif max_frac is not None and max_frac > 1.0:
        quality = "broken"
        note = (f"{lines_present}/7 ЕРН-линий найдено, но max|Δ| = "
                f"{max_res:.1f} кэВ ({max_frac:.0%} от FWHM) — "
                "калибровка не валидна.")
    elif max_frac is not None and max_frac > drift_threshold_fwhm_fraction:
        quality = "drift"
        note = (f"{lines_present}/7 ЕРН-линий найдено, max|Δ| = "
                f"{max_res:.1f} кэВ ({max_frac:.0%} от FWHM) — "
                f"допустимый дрейф (>{drift_threshold_fwhm_fraction:.0%} FWHM). "
                "Сохранённая калибровка приемлема для идентификации, но не "
                "для активностей ниже 1 σ_E.")
    else:
        quality = "ok"
        note = (f"{lines_present}/7 ЕРН-линий найдено, max|Δ| = "
                f"{max_res:.1f} кэВ ({max_frac:.0%} от FWHM при пороге "
                f"{drift_threshold_fwhm_fraction:.0%}) — калибровка ok.")

    return SevenLineCheck(
        per_line=per_line,
        lines_present=lines_present,
        max_residual_keV=max_res,
        mean_residual_keV=mean_res,
        max_residual_fwhm_fraction=max_frac,
        quality=quality,
        quality_note=note,
    )


__all__ = [
    "SEVEN_LINES",
    "LineCheckResult", "SevenLineCheck",
    "run_seven_line_check",
]
