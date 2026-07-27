"""
F-288 (v1.17.17, T-015 + T-016 + T-036 + T-037 + T-057 + T-080) —
Statistical thresholds per Currie 1968 / ISO 11929:2019 / ЛСРМ §6.3.

Сборник unified detection / decision / quantification thresholds для
γ-спектрометрии, в одном месте, со ссылками на каноничные источники.

  • **T-015 Currie paired-blank L_D** — `L_D = 2.71 + 4.65·√bg` (Currie
    1968) для случая, когда фон оценивается из ОТДЕЛЬНОГО фонового
    измерения той же длительности (paired-blank), а не из крыльев
    target ROI. Это базовый Lochamy/Currie протокол.

  • **T-016 ISO 11929 upper-limit semantics** — когда detection failed
    (A < L_D), отчёт активности должен быть в формате
    `< L_U = A + 1.645·σ_A`, а НЕ просто `< L_D`. Различие важно для
    регуляторов: upper-limit — physical bound на возможную активность,
    L_D — характеристика прибора.

  • **T-036 L_Q (limit of quantitation)** — `L_Q = 10·σ_0`, или
    эквивалентно `L_Q ≈ 14.1 + 14.1·√bg` для paired-blank. Activity
    может быть отчётна с относительной точностью ≤10% только если
    A ≥ L_Q. Между L_D и L_Q — detectable but not quantifiable.

  • **T-037 MIA / LLD / MDA / L_l-t labelling** — terminology gloss.
    Не алгоритмическое; помогает регулятору быстро понять, какой
    threshold используется. См. `THRESHOLD_GLOSSARY`.

  • **T-057 ISO 11929 quadratic Eq.6.3-6** — полная форма L_D с
    type-B uncertainty: `L_D² − 2·L_D·(σ_0² + L_C²) + L_C⁴ = 0`,
    учитывает u_rel(g) (efficiency, intensity, background-time
    uncertainty). Для non-trivial uncertainty даёт большее значение
    чем простая `L_D = 2·L_C + k²`.

  • **T-080 ONE-best-line MDA** — для нуклида с N линиями отчёт MDA
    делается по ОДНОЙ best line (highest yield·ε / √bg), а НЕ как
    среднее по линиям. Currie/Lochamy: averaging нескольких линий
    смешивает разные характеристики прибора и даёт misleading MDA.

References
----------
- Currie L.A., Anal. Chem. 40 (1968) 586 — оригинальные L_C/L_D/L_Q
- ISO 11929:2019 — Determination of characteristic limits
- ЛСРМ Algorithmic Foundations 2022 §6.3 (formulas 6.3-6, 6.3-7)
- Lochamy J.C., "The Minimum-Detectable-Activity Concept" 1981
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Sequence


# k-quantiles (one-sided 95% / 99%).
K_ALPHA_95 = 1.645
K_ALPHA_99 = 2.326


# ──────────────────────────────────────────────────────────────────
# T-015 — Currie paired-blank L_D
# ──────────────────────────────────────────────────────────────────

def currie_paired_blank_LD(background_counts: float, k: float = K_ALPHA_95) -> float:
    """Currie 1968 paired-blank detection limit.

        L_D = k² + 2·k·√bg   (≈ 2.71 + 4.65·√bg для k=1.645)

    Применять если фон оценивается из ОТДЕЛЬНОГО bg измерения той же
    длительности что и sample.
    """
    bg = max(0.0, float(background_counts))
    return k * k + 2.0 * k * math.sqrt(bg)


def currie_paired_blank_LC(background_counts: float, k: float = K_ALPHA_95) -> float:
    """Currie 1968 paired-blank decision threshold.

        L_C = k·√(2·bg)
    """
    bg = max(0.0, float(background_counts))
    return k * math.sqrt(2.0 * bg)


# ──────────────────────────────────────────────────────────────────
# T-016 — ISO 11929 upper-limit semantics
# ──────────────────────────────────────────────────────────────────

def iso_11929_upper_limit_Bq(
    A_Bq: float, sigma_A_Bq: float, k: float = K_ALPHA_95,
) -> float:
    """ISO 11929:2019 upper-limit L_U для случая failed detection.

        L_U = A + k·σ_A   (one-sided 95% CI upper bound)

    Это **physical** upper bound на activity, не characteristic
    instrument threshold. Использовать в reports как `< L_U`, а НЕ
    `< L_D`.
    """
    return float(A_Bq) + k * float(sigma_A_Bq)


def report_as_upper_limit(
    A_Bq: float, sigma_A_Bq: float, L_C_Bq: float,
    k: float = K_ALPHA_95,
) -> bool:
    """Решить, должен ли результат отчётываться как `< L_U` (upper limit).

    Правило: A < L_C → результат не значим → upper limit reporting.
    """
    return float(A_Bq) < float(L_C_Bq)


# ──────────────────────────────────────────────────────────────────
# T-036 — L_Q (limit of quantitation)
# ──────────────────────────────────────────────────────────────────

K_QUANTIFICATION = 10.0   # default для ≤10% относительной точности


def limit_of_quantitation_LQ(
    background_counts: float, k_Q: float = K_QUANTIFICATION,
) -> float:
    """L_Q — порог количественного определения (Currie 1968).

        L_Q = k_Q · σ_0 = k_Q · √bg  (paired-blank)

    Activity может быть отчётна с relative uncertainty ≤ 1/k_Q
    только если A ≥ L_Q. k_Q=10 → ≤10%; k_Q=3 → ≤33% (loose).
    """
    bg = max(0.0, float(background_counts))
    return k_Q * math.sqrt(bg)


# ──────────────────────────────────────────────────────────────────
# T-037 — MIA / LLD / MDA / L_l-t glossary
# ──────────────────────────────────────────────────────────────────

THRESHOLD_GLOSSARY = {
    "L_C": (
        "Decision threshold (Currie): if N_observed > L_C → presence "
        "declared with probability (1-α) of false positive. Used "
        "ALONE for sample-vs-bg discrimination."
    ),
    "L_D": (
        "Detection limit (Currie): smallest true activity that would "
        "be detected with probability (1-β). Reports as `MDA = L_D/(ε·I·t)`. "
        "ISO 11929 calls this 'detection limit'."
    ),
    "L_Q": (
        "Quantification limit: smallest activity at which result is "
        "quantifiable with target relative uncertainty (default ≤10% "
        "at k_Q=10). Between L_D and L_Q — detectable but not "
        "quantifiable."
    ),
    "L_U": (
        "Upper limit (ISO 11929 §6): used WHEN A < L_C — report "
        "`A_Bq < L_U = A + 1.645·σ_A` instead of `< L_D` or `not detected`. "
        "Provides physical bound on possible activity."
    ),
    "MIA": (
        "Minimum Identifiable Activity (alt. Minimum Detectable "
        "Identifiable Activity): obsolete synonym for L_D in MDA-context "
        "reports. Kept for compatibility with legacy ОСГИ certificates."
    ),
    "LLD": (
        "Lower Limit of Detection (older terminology): in modern usage = "
        "L_D in counts. Different from MDA which is L_D converted to Bq."
    ),
    "MDA": (
        "Minimum Detectable Activity (in Bq or Bq/kg): canonical term for "
        "L_D / (ε · I · t). Use this in all customer-facing reports."
    ),
    "L_l-t": (
        "Long-term detection limit (МИ 1916-88): L_D averaged over many "
        "repeat measurements; used for routine-monitoring threshold."
    ),
}


def explain_threshold(label: str) -> str:
    return THRESHOLD_GLOSSARY.get(label.upper().replace("-", "_"),
                                   THRESHOLD_GLOSSARY.get(label, ""))


# ──────────────────────────────────────────────────────────────────
# T-057 — ISO 11929 quadratic L_D (Eq. 6.3-6)
# ──────────────────────────────────────────────────────────────────

def iso_11929_LD_quadratic(
    *,
    L_C_counts: float,
    sigma_0_counts: float,
    u_rel_g: float,
    k: float = K_ALPHA_95,
) -> float:
    """ISO 11929:2019 §6.3 Eq. 6.3-6 — quadratic L_D.

    Учитывает type-B uncertainty `u_rel(g)` (efficiency, intensity,
    measurement-time uncertainty) — в дополнение к counting statistics.

    Solve for L_D:
        L_D² − 2·k·L_D·(σ_0² + u_rel_g²·L_C²)^(1/2) − L_C² = 0
        L_D = L_C + k·√(σ_0² + u_rel_g² · L_D²)

    Iterative solution (Newton-like):
    """
    L_D = 2.0 * float(L_C_counts) + k * k    # initial guess (simple formula)
    sigma_0_sq = float(sigma_0_counts) ** 2
    u_rel_sq = float(u_rel_g) ** 2
    for _ in range(50):
        rhs = float(L_C_counts) + k * math.sqrt(
            max(0.0, sigma_0_sq + u_rel_sq * L_D * L_D)
        )
        if abs(rhs - L_D) < 1e-6 * max(1.0, L_D):
            return rhs
        L_D = 0.5 * (L_D + rhs)
    return L_D


# ──────────────────────────────────────────────────────────────────
# T-080 — ONE-best-line MDA
# ──────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class _LineCandidate:
    line_E_keV: float
    intensity_decimal: float    # I (decimal, not percent)
    efficiency: float           # ε at this energy
    bg_counts_in_ROI: float


def best_mda_line(
    line_candidates: Sequence[_LineCandidate],
) -> Optional[_LineCandidate]:
    """T-080 — найти line с минимальной MDA для нуклида.

    Метрика: best_line = argmin_i MDA_i, где MDA_i ∝ √bg_i / (I_i · ε_i).
    Equivalently — argmax_i (I_i · ε_i) / √bg_i.

    Returns
    -------
    The candidate с наименьшей MDA, или None если список пуст.

    Контракт T-080 / Currie: НЕ усреднять MDA по линиям нуклида;
    отчёт по ОДНОЙ best line, потому что averaging смешивает разные
    instrument characteristics (efficiency, bg) и даёт misleading MDA.
    """
    candidates = list(line_candidates)
    if not candidates:
        return None

    def _score(c: _LineCandidate) -> float:
        bg_sqrt = math.sqrt(max(c.bg_counts_in_ROI, 1.0))
        denom = c.intensity_decimal * c.efficiency
        if denom <= 0:
            return 0.0
        return denom / bg_sqrt   # больше — лучше (меньше MDA)

    best = max(candidates, key=_score)
    return best


__all__ = [
    "K_ALPHA_95", "K_ALPHA_99", "K_QUANTIFICATION",
    "currie_paired_blank_LD", "currie_paired_blank_LC",
    "iso_11929_upper_limit_Bq", "report_as_upper_limit",
    "limit_of_quantitation_LQ",
    "THRESHOLD_GLOSSARY", "explain_threshold",
    "iso_11929_LD_quadratic",
    "best_mda_line",
]
