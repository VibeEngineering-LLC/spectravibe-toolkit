"""
F-296 (v1.17.20, T-007 + T-082) — True Coincidence Summing correction
                                    for close-geometry NaI.

Когда несколько γ-квантов каскадного распада регистрируются в пределах
**resolving time** детектора, они суммируются в одно событие — обычно
**ВНЕ** FEP-пика своих исходных энергий. Это называется True Coincidence
Summing (TCS), и в close geometry (Marinelli на NaI) эффект может
достигать 5-30 % для каскадных нуклидов (Co-60, Eu-152, Ba-133).

Формула (Gilmore §8.6, упрощённая для FEP-пика E_i)
---------------------------------------------------
Доля «потерянных» событий для линии E_i из-за TCS с другой линией E_j:

    L_ij = p_ij · ε_T(E_j) · W(θ_ij)

где:
  • p_ij — angular-corrected coincidence probability (берётся из
           branching ratios + angular correlation коэффициентов)
  • ε_T(E_j) = ε_FEP(E_j) / P/T(E_j) — total efficiency на E_j
  • W(θ_ij) — angular correlation factor; для NaI close geometry ≈ 1
              (большой solid angle усредняет анизотропию)

Correction factor (multiplicative, > 1 means measured area
underestimates true rate):

    C_TCS,i = 1 / (1 - Σ_j L_ij)

True activity:
    A_true = A_apparent · C_TCS

«Sum-IN» effect (когда два кванта попадают вместе в один пик-сумма)
учитывается отдельным additive term — он обычно мал для NaI
(меньше resolving time, шире energy-window resolution).

Cascade probability data
------------------------
В этом модуле определены preset-каскады для часто используемых
эталонов (Co-60, Eu-152, Cs-137 одиночный, Ba-133, Y-88), с
табличными branching ratios. Полные decay-схемы для произвольных
нуклидов требуют parser ENSDF — out of scope.

Допущения
---------
1. Angular correlation W(θ) = 1 (close geometry NaI Marinelli).
2. Sum-IN effect ignored для NaI low-resolution (peaks overlap anyway).
3. Random coincidences (chance pile-up) обрабатываются отдельным модулем.

References
----------
- Gilmore & Joss «Practical Gamma-ray Spectrometry» 3rd Ed. § 8.6
- ЛСРМ Algorithmic Foundations 2022 § 10 «Каскадное суммирование»
- Schima FJ, Hoppes DD «Tables for cascade summing corrections» (1983)
- Debertin K, Helmer RG «Gamma- and X-ray spectrometry» (1988) §4.5
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence


@dataclass(frozen=True)
class CascadeLine:
    """Одна γ-линия в каскадной decay-схеме."""
    E_keV: float
    branching_ratio: float    # P(emission per disintegration), decimal


@dataclass(frozen=True)
class CascadePair:
    """Пара коинцидентных линий E_i и E_j (для расчёта L_ij)."""
    E_i_keV: float
    E_j_keV: float
    coincidence_probability: float    # p_ij


# Preset-каскады для часто используемых эталонов.
# Source: ENSDF (NNDC), упрощено до главных пар.

CO60_CASCADE = [
    CascadeLine(1173.2, 0.9985),
    CascadeLine(1332.5, 0.9998),
]
# Co-60 — практически 100 % cascade between two photons.
CO60_PAIRS = [
    CascadePair(1173.2, 1332.5, coincidence_probability=0.998),
    CascadePair(1332.5, 1173.2, coincidence_probability=0.998),
]

EU152_CASCADE = [
    CascadeLine(121.78, 0.2858),
    CascadeLine(244.70, 0.0755),
    CascadeLine(344.28, 0.2658),
    CascadeLine(778.90, 0.1297),
    CascadeLine(964.06, 0.1462),
    CascadeLine(1085.84, 0.1015),
    CascadeLine(1112.07, 0.1356),
    CascadeLine(1408.01, 0.2085),
]
# Для Eu-152 — сложный каскад; preset покрывает доминантные пары.
EU152_PAIRS = [
    CascadePair(121.78, 1408.01, coincidence_probability=0.135),
    CascadePair(344.28, 1112.07, coincidence_probability=0.102),
    CascadePair(778.90, 344.28,  coincidence_probability=0.080),
    CascadePair(964.06, 411.12,  coincidence_probability=0.055),
    CascadePair(1085.84, 344.28, coincidence_probability=0.060),
    CascadePair(1112.07, 344.28, coincidence_probability=0.102),
    CascadePair(1408.01, 121.78, coincidence_probability=0.135),
]

CS137_CASCADE = [
    CascadeLine(661.66, 0.851),    # одиночная линия, no significant cascades
]
CS137_PAIRS: List[CascadePair] = []

BA133_CASCADE = [
    CascadeLine(80.998, 0.3239),
    CascadeLine(276.40, 0.0716),
    CascadeLine(302.85, 0.1834),
    CascadeLine(356.01, 0.6205),
    CascadeLine(383.85, 0.0894),
]
BA133_PAIRS = [
    CascadePair(80.998, 356.01, coincidence_probability=0.18),
    CascadePair(276.40, 80.998, coincidence_probability=0.05),
    CascadePair(302.85, 80.998, coincidence_probability=0.10),
    CascadePair(356.01, 80.998, coincidence_probability=0.18),
    CascadePair(383.85, 80.998, coincidence_probability=0.06),
]

# F-148 (v1.18.5) — Bi-214 cascade pairs (Ra-226 chain, U-238 series).
# Bi-214 → Po-214 (β⁻, 19.9 min) — dominant cascades между основными
# γ-линиями ЕРН-маркеров: 609.31, 1120.29, 1764.49 keV.
# Branching ratios per ENSDF (NNDC); coincidence_probability scaled
# с учётом angular correlation (NaI close geometry → W≈1).
BI214_CASCADE = [
    CascadeLine(609.31, 0.4549),
    CascadeLine(768.36, 0.0489),
    CascadeLine(1120.29, 0.1491),
    CascadeLine(1238.11, 0.0589),
    CascadeLine(1377.67, 0.0400),
    CascadeLine(1764.49, 0.1531),
    CascadeLine(2204.21, 0.0497),
]
# Pairs: главные совпадения, отвечающие за >1% TCS effect на Marinelli NaI.
BI214_PAIRS = [
    CascadePair(609.31, 1120.29, coincidence_probability=0.069),
    CascadePair(609.31, 1764.49, coincidence_probability=0.071),
    CascadePair(1120.29, 609.31, coincidence_probability=0.069),
    CascadePair(1764.49, 609.31, coincidence_probability=0.071),
    CascadePair(1238.11, 609.31, coincidence_probability=0.027),
]

# F-148 (v1.18.5) — Bi-212 cascade pairs (Th-232 chain).
# Bi-212 → Tl-208/Po-212; основные линии 727.33, 1620.50 keV.
BI212_CASCADE = [
    CascadeLine(727.33, 0.0658),
    CascadeLine(785.37, 0.0110),
    CascadeLine(1620.50, 0.0149),
]
BI212_PAIRS = [
    CascadePair(727.33, 1620.50, coincidence_probability=0.010),
    CascadePair(1620.50, 727.33, coincidence_probability=0.010),
]

# Lookup: (nuclide_id, list_of_pairs).
CASCADE_PRESETS: Dict[str, List[CascadePair]] = {
    "Co-60": CO60_PAIRS,
    "Eu-152": EU152_PAIRS,
    "Cs-137": CS137_PAIRS,
    "Ba-133": BA133_PAIRS,
    "Bi-214": BI214_PAIRS,
    "Bi-212": BI212_PAIRS,
}


@dataclass(frozen=True)
class TCSCorrectionResult:
    E_i_keV: float
    sum_L_ij: float              # Σ_j p_ij · ε_T(E_j)
    correction_factor: float     # C_TCS = 1/(1 - sum_L_ij)
    n_pairs_used: int
    is_significant: bool         # True если correction > 5 %


def compute_tcs_correction(
    E_i_keV: float,
    nuclide_pairs: Sequence[CascadePair],
    total_efficiency_func: Callable[[float], float],
    energy_match_tolerance_keV: float = 1.0,
    significant_threshold_pct: float = 5.0,
) -> TCSCorrectionResult:
    """Compute TCS correction factor для линии E_i нуклида.

    Parameters
    ----------
    E_i_keV : float
        Энергия линии, для которой считаем correction.
    nuclide_pairs : sequence of CascadePair
        Список коинцидентных пар (берётся из CASCADE_PRESETS[nuclide]).
    total_efficiency_func : callable
        Функция ε_T(E_j) — total efficiency на коинцидентной энергии.
        Используйте `pt_ratio_nai.total_efficiency_from_fep`.
    energy_match_tolerance_keV : float
        Допуск для парирования E_i с E_pair.E_i_keV (default 1 keV).
    significant_threshold_pct : float
        Порог для флага is_significant (default 5 %).

    Returns
    -------
    TCSCorrectionResult с C_TCS и метаданными.
    """
    sum_L = 0.0
    n_used = 0
    for pair in nuclide_pairs:
        if abs(pair.E_i_keV - E_i_keV) > energy_match_tolerance_keV:
            continue
        eps_T_j = total_efficiency_func(pair.E_j_keV)
        if eps_T_j <= 0:
            continue
        sum_L += pair.coincidence_probability * eps_T_j
        n_used += 1

    if sum_L >= 1.0:
        # Безопасный clamp — иначе деление на ≤0
        sum_L = 0.999

    C_TCS = 1.0 / (1.0 - sum_L)
    pct = 100.0 * (C_TCS - 1.0)
    return TCSCorrectionResult(
        E_i_keV=E_i_keV,
        sum_L_ij=sum_L,
        correction_factor=C_TCS,
        n_pairs_used=n_used,
        is_significant=pct >= significant_threshold_pct,
    )


def compute_tcs_correction_for_nuclide(
    E_i_keV: float,
    nuclide_id: str,
    total_efficiency_func: Callable[[float], float],
) -> TCSCorrectionResult:
    """Convenience: TCS correction по nuclide_id из CASCADE_PRESETS."""
    if nuclide_id not in CASCADE_PRESETS:
        raise KeyError(
            f"Unknown nuclide '{nuclide_id}'. "
            f"Known: {sorted(CASCADE_PRESETS.keys())}"
        )
    return compute_tcs_correction(
        E_i_keV=E_i_keV,
        nuclide_pairs=CASCADE_PRESETS[nuclide_id],
        total_efficiency_func=total_efficiency_func,
    )


def is_tcs_significant_for_geometry(
    sample_distance_cm: float, threshold_cm: float = 10.0,
) -> bool:
    """Эвристика: TCS значим при close geometry (≤ 10 cm от детектора).

    Для Marinelli (≈ 0 cm от поверхности) → ВСЕГДА True.
    Для point source 25 cm от детектора → False.
    """
    return sample_distance_cm <= threshold_cm


# ---------------------------------------------------------------------------
# #PTB-3 (2026-07-02) — Volume-source averaging, PTB-2018 Annex D Eq. (D17)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class VolumeElement:
    """Дискретный элемент объёма источника для D17-усреднения.

    Позиционная зависимость эффективностей факторизуется как
    ε(E; r_i) = s_i · ε_ref(E) — solid-angle-доминированное приближение
    (энергонезависимый позиционный масштаб, та же логика, что в
    efficiency-transfer по Ω̄-отношениям, PTB-2018 Eq. 33).

    Attributes
    ----------
    weight : float
        ΔV_i — объём (или доля объёма) элемента; нормируется внутри
        `compute_tcs_correction_volume`, абсолютная шкала не важна.
    fep_efficiency_scale : float
        s_p,i = ε_p(r_i) / ε_p,ref — масштаб FEP-эффективности в точке r_i
        относительно референсной (той, для которой калибрована
        total_efficiency_func).
    total_efficiency_scale : float
        s_t,i = ε_t(r_i) / ε_t,ref — масштаб total efficiency в точке r_i.
    """
    weight: float
    fep_efficiency_scale: float
    total_efficiency_scale: float


def compute_tcs_correction_volume(
    E_i_keV: float,
    nuclide_pairs: Sequence[CascadePair],
    total_efficiency_func: Callable[[float], float],
    volume_elements: Sequence[VolumeElement],
    energy_match_tolerance_keV: float = 1.0,
    significant_threshold_pct: float = 5.0,
) -> TCSCorrectionResult:
    """TCS correction для объёмного источника — PTB-2018 Annex D Eq. (D17).

    Для объёмного источника (Marinelli, сосуд Дента) вероятность
    регистрации в FEP-пике усредняется по объёму (Eq. D17):

        P̄_D(E_i) = (1/V) ∫_V ε_p(E_i; r) · [1 − Σ_j p_ij · ε_t(E_j; r)] dV

    Кажущаяся активность извлекается через объёмно-усреднённую
    FEP-эффективность ε̄_p, поэтому correction factor:

        C_TCS = ε̄_p / P̄_D
              = 1 / (1 − Σ_j p_ij · ε_T(E_j) · ⟨s_p·s_t⟩ / ⟨s_p⟩)

    где ⟨·⟩ — взвешенное по ΔV_i среднее, s_p/s_t — позиционные масштабы
    из `VolumeElement`. Эффективная total efficiency в loss-term — это
    **ε_p-взвешенное** объёмное среднее ε_t: элементы, дающие больший
    вклад в пик, сильнее страдают от summing-out.

    Полный расчёт Eq. D17 требует Monte Carlo (PTB-2018, Annex D.1);
    данная дискретизация — детерминированное приближение при
    факторизации ε(E; r) = s(r)·ε_ref(E). При всех масштабах = 1
    вырождается в точечную `compute_tcs_correction`.

    Parameters
    ----------
    E_i_keV, nuclide_pairs, total_efficiency_func,
    energy_match_tolerance_keV, significant_threshold_pct
        Как в `compute_tcs_correction`; total_efficiency_func — ε_T(E) в
        референсной точке (обычно центр/эффективный центр объёма).
    volume_elements : sequence of VolumeElement
        Дискретизация объёма. Масштабы s_p/s_t поставляет вызывающий
        (например, из efficiency-transfer solid-angle модели).

    Raises
    ------
    ValueError
        Пустой volume_elements; отрицательный weight/scale;
        нулевая суммарная FEP-взвешенная масса (Σ w_i·s_p,i ≤ 0).
    """
    if not volume_elements:
        raise ValueError("volume_elements must be non-empty")

    w_sp = 0.0      # Σ w_i · s_p,i
    w_sp_st = 0.0   # Σ w_i · s_p,i · s_t,i
    for el in volume_elements:
        if el.weight < 0 or el.fep_efficiency_scale < 0 or el.total_efficiency_scale < 0:
            raise ValueError(
                "VolumeElement weight/scales must be non-negative, got "
                f"(w={el.weight}, s_p={el.fep_efficiency_scale}, "
                f"s_t={el.total_efficiency_scale})"
            )
        w_sp += el.weight * el.fep_efficiency_scale
        w_sp_st += el.weight * el.fep_efficiency_scale * el.total_efficiency_scale

    if w_sp <= 0:
        raise ValueError("Σ weight·fep_efficiency_scale must be > 0")

    eff_total_factor = w_sp_st / w_sp   # ⟨s_p·s_t⟩ / ⟨s_p⟩

    sum_L = 0.0
    n_used = 0
    for pair in nuclide_pairs:
        if abs(pair.E_i_keV - E_i_keV) > energy_match_tolerance_keV:
            continue
        eps_T_j = total_efficiency_func(pair.E_j_keV)
        if eps_T_j <= 0:
            continue
        sum_L += pair.coincidence_probability * eps_T_j * eff_total_factor
        n_used += 1

    if sum_L >= 1.0:
        sum_L = 0.999   # тот же безопасный clamp, что в точечной версии

    C_TCS = 1.0 / (1.0 - sum_L)
    pct = 100.0 * (C_TCS - 1.0)
    return TCSCorrectionResult(
        E_i_keV=E_i_keV,
        sum_L_ij=sum_L,
        correction_factor=C_TCS,
        n_pairs_used=n_used,
        is_significant=pct >= significant_threshold_pct,
    )


__all__ = [
    "CascadeLine",
    "CascadePair",
    "CO60_CASCADE", "CO60_PAIRS",
    "EU152_CASCADE", "EU152_PAIRS",
    "CS137_CASCADE", "CS137_PAIRS",
    "BA133_CASCADE", "BA133_PAIRS",
    "BI214_CASCADE", "BI214_PAIRS",
    "BI212_CASCADE", "BI212_PAIRS",
    "CASCADE_PRESETS",
    "TCSCorrectionResult",
    "VolumeElement",
    "compute_tcs_correction",
    "compute_tcs_correction_for_nuclide",
    "compute_tcs_correction_volume",
    "is_tcs_significant_for_geometry",
]
