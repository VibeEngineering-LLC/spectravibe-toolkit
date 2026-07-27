"""
F-302..F-304 / v1.18.4 — High-level quasi-template activity solver.

Объединяет F-302 (PPP templates), F-303 (WLS fit), F-304 (Compton continuum),
F-295 (P/T ratio) и F-300 (FWHM-at-E) в один production-ready API:

    solve_quasi_template_activities(
        spectrum_counts, channel_to_keV, fwhm_at_E_func,
        nuclide_ids, live_time_s, ...
    ) -> list[ActivityResult]

Result совместим с downstream-кодом, который ожидает `ActivityResult`
(json_report, markdown_report, etc.) — `sigma_method="quasi_template"`,
`intra_chi2_per_dof=χ²_red`, `notes` содержит метаданные fit-а.

Это **opt-in алтернативный путь** для production pipeline — не замена
существующего `compute_activities_for_all` (peak-area-based), а
полностью независимый full-spectrum solver.

References
----------
- ЛСРМ Algorithmic Foundations 2022 § 13 «Квазишаблонный метод»
- F-302/F-303/F-304 docstrings
"""
from __future__ import annotations

import warnings
from typing import Callable, Iterable, List, Optional, Sequence

from gamma.activity.compute import ActivityResult


def solve_quasi_template_activities(
    spectrum_counts: Sequence[float],
    channel_to_keV: Callable[[float], float],
    fwhm_at_E_func: Callable[[float], float],
    efficiency_at_E_func: Callable[[float], float],
    nuclide_ids: Iterable[str],
    live_time_s: float,
    *,
    background_counts: Optional[Sequence[float]] = None,
    energy_window_keV: Optional[tuple] = None,
    detector_id: str = "Gamma-1S",
    enable_compton_continuum: bool = True,
    chi2_red_acceptance_max: float = 1.5,
) -> List[ActivityResult]:
    """Production-ready quasi-template activity solver.

    Parameters
    ----------
    spectrum_counts : Sequence[float]
        Observed gross spectrum (counts per channel).
    channel_to_keV : Callable[[float], float]
        Energy calibration channel → keV.
    fwhm_at_E_func : Callable[[float], float]
        Resolution calibration energy keV → FWHM keV.
    efficiency_at_E_func : Callable[[float], float]
        Photopeak efficiency curve ε_FEP(E).
    nuclide_ids : Iterable[str]
        Library nuclides для template construction (must exist в
        gamma.data.nuclide_library).
    live_time_s : float
    background_counts : Optional[Sequence[float]]
        External background spectrum (same length). None → zero baseline.
    energy_window_keV : Optional[(E_lo, E_hi)]
        Restrict fit to this energy range (excludes LLD/ULD).
    detector_id : str
        Used for F-295 P/T lookup (default 'Gamma-1S').
    enable_compton_continuum : bool
        If True — inject F-304 Compton pedestal in templates.
    chi2_red_acceptance_max : float
        Threshold для `is_acceptable` flag в notes.

    Returns
    -------
    list[ActivityResult] — один per nuclide_id, в порядке передачи.
    """
    if live_time_s <= 0:
        raise ValueError(f"live_time_s must be > 0, got {live_time_s}")

    from gamma.activity.quasi_template_ppp import (
        NuclideDef, NuclideLine,
        build_templates_for_library,
    )
    from gamma.activity.quasi_template_fit import solve_quasi_template_fit
    from gamma.data.nuclide_library import get_nuclide

    cont_func = None
    pt_func = None
    if enable_compton_continuum:
        from gamma.activity.compton_continuum import make_continuum_func
        from gamma.activity.pt_ratio_nai import pt_ratio_for_detector
        cont_func = make_continuum_func(fwhm_keV_at=fwhm_at_E_func)

        def _pt(E):
            try:
                return pt_ratio_for_detector(E, detector_id)
            except Exception as exc:  # DEEP-06
                warnings.warn(
                    f"P/T ratio lookup failed for detector="
                    f"{detector_id!r} at E={float(E):.2f} keV ({exc!r}); "
                    f"falling back to 1.0 — no Compton-pedestal scaling "
                    f"will be applied at this energy.",
                    stacklevel=2,
                )
                return 1.0
        pt_func = _pt

    nuclide_defs: List[NuclideDef] = []
    skipped_ids: List[str] = []
    for nid in nuclide_ids:
        rec = get_nuclide(nid)
        if not rec:
            skipped_ids.append(nid)
            continue
        lib_lines = rec.get("lines", [])
        if not lib_lines:
            skipped_ids.append(nid)
            continue
        ppp_lines: List[NuclideLine] = []
        for ll in lib_lines:
            E = float(ll[0])
            I_pct = float(ll[1])
            if E <= 0 or I_pct <= 0:
                continue
            eps = efficiency_at_E_func(E)
            if eps is None or eps <= 0:
                continue
            ppp_lines.append(NuclideLine(
                E_keV=E, intensity=I_pct / 100.0, efficiency=eps,
            ))
        if ppp_lines:
            nuclide_defs.append(NuclideDef(nuclide_id=nid, lines=ppp_lines))

    if not nuclide_defs:
        return []

    n_channels = len(spectrum_counts)

    def _ch_to_keV(ch):
        return channel_to_keV(float(ch))

    templates = build_templates_for_library(
        nuclide_defs, n_channels, _ch_to_keV,
        fwhm_at_E_func, cont_func, pt_func,
    )

    energy_window_channels = None
    if energy_window_keV is not None:
        E_lo, E_hi = energy_window_keV
        ch_lo = max(0, int(round((E_lo - channel_to_keV(0.0))
                                 / max(channel_to_keV(1.0) - channel_to_keV(0.0), 1e-9))))
        ch_hi = min(
            n_channels,
            int(round((E_hi - channel_to_keV(0.0))
                      / max(channel_to_keV(1.0) - channel_to_keV(0.0), 1e-9))) + 1,
        )
        if ch_hi > ch_lo:
            energy_window_channels = (ch_lo, ch_hi)

    fit_res = solve_quasi_template_fit(
        observed=spectrum_counts,
        templates=templates,
        live_time_s=live_time_s,
        background=background_counts,
        energy_window=energy_window_channels,
    )

    results: List[ActivityResult] = []
    is_accepted = fit_res.is_accepted(chi2_red_acceptance_max)
    note_base = (
        f"F-302..F-304 / v1.18.4 quasi-template "
        f"χ²_red={fit_res.chi2_red:.3f}, dof={fit_res.dof}, "
        f"converged={fit_res.converged}, "
        f"is_accepted={is_accepted}, "
        f"detector={detector_id}"
    )
    if fit_res.notes:
        note_base += " | " + "; ".join(fit_res.notes)
    if skipped_ids:
        note_base += (
            f" | skipped_unknown_nuclides: {','.join(skipped_ids[:5])}"
            + ("…" if len(skipped_ids) > 5 else "")
        )

    for nuc_def in nuclide_defs:
        nid = nuc_def.nuclide_id
        A = float(fit_res.activities.get(nid, 0.0))
        sigma = float(fit_res.sigma_activities.get(nid, 0.0))
        results.append(ActivityResult(
            nuclide=nid,
            A_Bq=A, sigma_A_Bq=sigma,
            lines_used=(),
            lines_skipped=(),
            intra_chi2_per_dof=float(fit_res.chi2_red)
            if fit_res.converged else None,
            sigma_method="quasi_template",
            from_bg_subtracted=(background_counts is not None),
            force_gross_override=False,
            notes=note_base,
        ))
    return results


__all__ = [
    "solve_quasi_template_activities",
]
