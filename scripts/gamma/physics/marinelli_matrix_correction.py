"""
Wave-3 / 2026-06-04 — Marinelli matrix calibration diagnostics (Agent A).

ОПИСАНИЕ
========

Этот модуль — **diagnostic-and-research** компаньон к
`gamma.physics.self_attenuation` (F-122 / K-20 thin-slab self-attenuation
correction). Он НЕ подключается в production-пайплайн compute_activities
по умолчанию.

Контекст
--------

Wave 2 на 4 эталонных Marinelli 1 L фикстурах (Cs-137 / K-40 / Ra-226 / Th-232,
паспорт цитирован в tests/snapshot/test_marinelli_certificate_residual.py)
показал residuals:

  • Ra-226 → Pb-214 (295/352 keV)  =  −1.2 %   ✓ (Phase-1 exit-criterion <5%)
  • Th-232 → Pb-212 (238 keV)      =  −11.5 %
  • K-40 direct (1460.83 keV)      =  −12.8 %
  • Cs-137 direct (661.66 keV)     =  not detected on the fixture

Wave-3 phase-1 диагностика (см. outbox `2026-06-04_wave3_…`) обнаружила:

  1. **F-122 silently disabled by canonical-name lookup bug.**
     `compute_activities_for_all` передаёт `geometry_canonical="marinelli_1L"`,
     но `REF_GEOMETRY` в `gamma.physics.self_attenuation` зарегистрирован
     по Cyrillic key `"Маринелли"`. Lookup промахивается, F-122 не
     активируется. Этот промах *эмпирически защищает* current результаты
     (см. ниже п. 2).

  2. **Применение F-122 (canonical-alias patched) ухудшает residuals.**
     Empirical run на тех же фикстурах с `REF_GEOMETRY["marinelli_1L"]`
     зарегистрированным:

       | Kit / line          | base resid | F-122 resid | direction |
       |---------------------|------------|-------------|-----------|
       | Ra-226 / Pb-214 295 | −1.2 %     | −14.5 %     | WORSE     |
       | K-40 1460           | −12.8 %    | −18.3 %     | WORSE     |
       | Th-232 / Pb-212 238 | −11.5 %    | −11.5 %     | unchanged |
       | Cs-137 661          | n/d        | n/d         | n/a       |

     Причина: K-20 thin-slab `F = (1 − exp(−μρt))/(μρt)` с
     ρ_ref=1.6 g/cm³ / ρ_sample≈0.6 g/cm³ даёт corr ≈ 0.86 (Pb-214 295)
     и 0.94 (K-40 1460). A_corrected = A_meas × corr < A_meas. Но
     emperical A_meas УЖЕ под-оценивает паспорт (residual отрицательный).
     Дополнительное умножение на <1 удаляет последний оставшийся резерв.

  3. **Per-block .efr точки** имеют каждый свой ρ_ref (Th=1.6, K=1.55,
     Cs=1.66, Ra=1.67 g/cm³ — почти одинаковые). Polynomial degree-3 fit
     по 15 anchors из 4 разных блоков сглаживает per-block структуру.
     Худший residual fit↔anchor: Pb-214 295 keV = +13.4 % (fit eps
     ВЫШЕ anchor eps на 13 %, и это *компенсирует* что-то ещё, видимо
     библиотечную intensity или background subtraction).

Что делает этот модуль
----------------------

Предоставляет **opt-in research-grade** утилиты:

  • `nearest_block_anchor(E)` — найти ближайшую (по log-E) калибровочную
    anchor-точку из .efr и вернуть её eps вместо polynomial fit.
    Применимо для диагностики «сколько процентов residual'а идёт от
    polynomial smoothing vs от реальной physics».

  • `per_block_efficiency_lookup(efr_path, energy_keV, nuclide)` — найти
    anchor, лежащий в блоке именно для запрошенного nuclide. Это
    *самый точный* способ (использует калибровочную точку для той же
    линии), но только для энергий, фактически измеренных в калибровке.

  • `diagnose_residual_decomposition(...)` — структурный отчёт:
    decompose observed residual into (polynomial-fit deviation,
    per-block density mismatch, residual unknown).

Эти утилиты НЕ изменяют production compute path. Они — инструмент для
дальнейшего расследования; wave-4 может на их основе решить, какой
из подходов даёт лучший balance precision/robustness.

References
----------
* ЛСРМ Algorithmic Foundations 2022, §8.5 «efficiency calibration» —
  per-source vs aggregated (.efr vs .efa) discussion.
* Gilmore & Joss 3rd ed §8.7 «Sample-related effects» — self-attenuation
  derivation, limit of thin-slab approximation.
* Knoll 4th ed §10.III.5 — volume source self-absorption integral form.
* Currie L.A., Anal Chem 40 (1968) 586 — efficiency uncertainty
  propagation (referenced via RAG-005).
* NIST XCOM Photon Cross Section Database (Berger et al., 2010), used
  by `gamma.physics.self_attenuation` for μ/ρ tables.

References to specific files / lines
-----------------------------------

* Wave-2 outbox: `_state/agent_a/outbox/2026-06-04_wave2_marinelli_cert_residual.md`
  (4-kit residual table, lines 28-39).
* F-122 implementation: `scripts/gamma/physics/self_attenuation.py` lines
  168-182 (`REF_GEOMETRY` Cyrillic-only key — the lookup bug).
* F-122 wiring: `scripts/gamma/activity/compute.py` lines 1060-1092
  (`if geometry_canonical and sample_density_g_cm3 is not None`).
* Per-nuclide .efr blocks: `detectors/Gamma-1S/efficiency/Gamma-1S_NaI_63x63_USB_SN-01/
  УДС-ГЦ-63х63-USB__SN-01_-_Маринелли.efr` lines 1, 22, 43, 64 (4
  header blocks for Th-232/K-40/Cs-137/Ra-226 respectively).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Sequence, Tuple


# ──────────────────────────────────────────────────────────────────────
# Convenience: canonical-name table for F-122 lookup (research-only).
# ──────────────────────────────────────────────────────────────────────
#
# NOTE: deliberately *separate* from `self_attenuation.REF_GEOMETRY`. We
# do NOT extend `REF_GEOMETRY` itself because empirical wave-3 evidence
# (outbox 2026-06-04_wave3_…) shows that activating F-122 on the four
# Marinelli reference fixtures DEGRADES residuals (Ra-226 / Pb-214
# anchor jumps from −1.2 % to −14.5 %; K-40 from −12.8 % to −18.3 %).
# Until wave-4 understands the over-correction direction, the F-122
# silent disable via canonical-name mismatch acts as an inadvertent
# safety net. Modifying `REF_GEOMETRY` directly would activate the
# regression on the entire fleet.
CANONICAL_REF_GEOMETRY_ALIASES: dict = {
    "marinelli_1L":  "Маринелли",
    "marinelli_05L": "Маринелли-0.5",   # not in REF_GEOMETRY yet, reserved
}


@dataclass(frozen=True)
class AnchorLookupResult:
    """Result of nearest-anchor efficiency lookup."""

    energy_keV: float
    nearest_anchor_energy_keV: float
    nearest_anchor_efficiency: float
    nearest_anchor_nuclide: str
    log_distance: float       # |log(E_q) − log(E_anchor)|


def nearest_block_anchor(
    energy_keV: float,
    anchors: Sequence[Tuple[float, float, str]],
) -> Optional[AnchorLookupResult]:
    """
    Find the calibration anchor closest to a query energy in log-E space.

    Parameters
    ----------
    energy_keV : float
        Query energy (keV).
    anchors : sequence of (E_keV, efficiency, nuclide) triples
        Calibration anchors, typically harvested from an .efr file via
        `gamma.io.lsrm_efficiency.read_efficiency_file`.

    Returns
    -------
    AnchorLookupResult or None if `anchors` is empty.

    Notes
    -----
    Log-space distance is used because efficiency varies smoothly in
    log-E (cf. Gilmore §8.5: «efficiency follows a near-power-law in E»).
    """
    if not anchors:
        return None
    if energy_keV <= 0:
        return None
    log_q = math.log(energy_keV)
    best = min(anchors, key=lambda a: abs(math.log(a[0]) - log_q))
    return AnchorLookupResult(
        energy_keV=energy_keV,
        nearest_anchor_energy_keV=float(best[0]),
        nearest_anchor_efficiency=float(best[1]),
        nearest_anchor_nuclide=str(best[2]),
        log_distance=abs(math.log(best[0]) - log_q),
    )


def per_block_efficiency_lookup(
    energy_keV: float,
    nuclide: str,
    anchors: Sequence[Tuple[float, float, str]],
    tolerance_keV: float = 1.0,
) -> Optional[float]:
    """
    Return ε from the calibration block matching `nuclide`, if an anchor
    falls within `tolerance_keV` of `energy_keV`.

    Parameters
    ----------
    energy_keV : float
        Query energy (keV).
    nuclide : str
        Nuclide name (case-insensitive); typically the same name used in
        the .efr Material/source line (e.g. "K-40", "Th-232", "Cs-137").
    anchors : sequence of (E_keV, efficiency, nuclide) triples
    tolerance_keV : float, default 1.0
        Maximum |E_q − E_anchor| (keV) to accept the anchor.

    Returns
    -------
    Anchor efficiency ε, or None if no anchor for that nuclide is within
    tolerance.

    Examples
    --------
    >>> anchors = [(1460.822, 9.743e-3, "K-40"),
    ...            (661.657, 1.871e-2, "Cs-137")]
    >>> per_block_efficiency_lookup(1460.8, "K-40", anchors)
    0.009743
    >>> per_block_efficiency_lookup(1460.8, "K-40", anchors,
    ...                             tolerance_keV=0.001) is None
    True
    """
    if not anchors:
        return None
    q_nuc = nuclide.strip().lower()
    for E_a, eps_a, nuc_a in anchors:
        if nuc_a.strip().lower() != q_nuc:
            continue
        if abs(float(E_a) - float(energy_keV)) <= tolerance_keV:
            return float(eps_a)
    return None


def harvest_anchors_from_efr(efr_path: str) -> list:
    """
    Read all (E_keV, ε, nuclide) anchor triples from an .efr file.

    Convenience wrapper around `gamma.io.lsrm_efficiency.read_efficiency_file`.

    Parameters
    ----------
    efr_path : str
        Path to a Lsrm SpectraLine .efr or .efa file.

    Returns
    -------
    list of (E_keV, efficiency, nuclide) triples, sorted by ascending E_keV.

    Notes
    -----
    Used by `diagnose_residual_decomposition` and the wave-3 sanity test
    `tests/snapshot/test_marinelli_matrix_correction.py`.
    """
    # Late import — keeps this module free of mandatory I/O dependency
    # for the simple lookup utilities above.
    from gamma.io.lsrm_efficiency import read_efficiency_file
    ef = read_efficiency_file(efr_path)
    triples = []
    for block in ef.blocks:
        for p in block.points:
            triples.append(
                (float(p.energy_keV),
                 float(p.efficiency),
                 str(p.source_nuclide))
            )
    triples.sort(key=lambda t: t[0])
    return triples


@dataclass(frozen=True)
class ResidualDecomposition:
    """Decomposition of an observed certificate residual into contributors."""

    energy_keV: float
    nuclide: str
    fit_efficiency: float            # ε from polynomial fit
    anchor_efficiency: Optional[float]   # ε at nearest matching anchor (None if absent)
    fit_minus_anchor_pct: Optional[float]   # (fit − anchor) / anchor × 100
    f122_correction: Optional[float]        # F_ref / F_sample (None if cannot compute)
    f122_residual_shift_pct: Optional[float]   # 100 × (corr − 1.0)
    notes: str = ""


def diagnose_residual_decomposition(
    *,
    energy_keV: float,
    nuclide: str,
    fit_efficiency: float,
    anchors: Sequence[Tuple[float, float, str]] = (),
    rho_sample_g_cm3: Optional[float] = None,
    rho_ref_g_cm3: float = 1.6,
    thickness_cm: float = 3.1,
) -> ResidualDecomposition:
    """
    Decompose the contribution of efficiency-fit smoothing and K-20
    self-attenuation correction at one energy.

    Useful for wave-3 diagnostic outbox tables and the sanity test.

    Returns
    -------
    ResidualDecomposition (frozen dataclass).

    Notes
    -----
    Anchor diff direction:  fit/anchor > 1 → A_calc UNDER-estimated by
    that factor (fit ε is too large, A = S/(ε I t) → too small).

    F-122 direction: corr < 1 → if applied, A_calc multiplied by corr →
    A becomes smaller (more under-estimation).
    """
    anchor_eps = per_block_efficiency_lookup(
        energy_keV, nuclide, anchors, tolerance_keV=2.0,
    )
    fit_minus_anchor_pct: Optional[float] = None
    if anchor_eps is not None and anchor_eps > 0:
        fit_minus_anchor_pct = (fit_efficiency - anchor_eps) / anchor_eps * 100.0

    f122_corr: Optional[float] = None
    if rho_sample_g_cm3 is not None and rho_sample_g_cm3 > 0:
        try:
            from gamma.physics.self_attenuation import (
                correction_factor as _cf,
                OISN_16_COMPOSITION as _oisn,
            )
            f122_corr = _cf(
                energy_keV,
                rho_sample_g_cm3=rho_sample_g_cm3,
                rho_ref_g_cm3=rho_ref_g_cm3,
                thickness_cm=thickness_cm,
                composition=_oisn,
            )
        except Exception:  # pragma: no cover — defensive
            f122_corr = None

    f122_shift_pct: Optional[float] = None
    if f122_corr is not None:
        f122_shift_pct = (f122_corr - 1.0) * 100.0

    return ResidualDecomposition(
        energy_keV=float(energy_keV),
        nuclide=nuclide,
        fit_efficiency=float(fit_efficiency),
        anchor_efficiency=anchor_eps,
        fit_minus_anchor_pct=fit_minus_anchor_pct,
        f122_correction=f122_corr,
        f122_residual_shift_pct=f122_shift_pct,
        notes=(
            "fit_minus_anchor_pct > 0 ⇒ polynomial fit overpredicts ε at this "
            "energy; A_calc under-estimated by that fraction. "
            "f122_residual_shift_pct < 0 ⇒ applying F-122 would further "
            "reduce A_calc (sample less dense than ref)."
        ),
    )


__all__ = [
    "CANONICAL_REF_GEOMETRY_ALIASES",
    "AnchorLookupResult",
    "ResidualDecomposition",
    "nearest_block_anchor",
    "per_block_efficiency_lookup",
    "harvest_anchors_from_efr",
    "diagnose_residual_decomposition",
]
