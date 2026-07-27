"""F-132 / v1.17.7 — Обязательная оценка стоимости анализа в токенах.

Зафиксировано навсегда в HTML / Markdown отчётах:

  * **HTML**: ВСЕГДА выводит итоговый footer с количеством токенов и
    процентом от 5-часовой бесплатной сессии. Раньше (до F-132) footer
    показывался только если пользователь передавал CLI флаги
    `--cost-tokens`/`--cost-session-pct`/`--cost-detail`; теперь
    отсутствие флагов → выводится авто-оценка из этого модуля.

  * **Markdown**: новый обязательный раздел «Оценка стоимости анализа»
    с таблицей по этапам Step 1-11 (парсинг → калибровка → поиск пиков →
    идентификация → деконволюция → активности → MDA → классификация →
    отчёт) + итог + % от сессии.

Модель оценки — эвристическая, основана на типичной активности
агента-аналитика на каждом этапе. Чем сложнее спектр (больше пиков,
нуклидов, мультиплетов), тем выше оценка. По умолчанию полный
анализ простого спектра (Cs-137) стоит ~6-10k токенов; сложного
многонуклидного спектра (Th-232 chain + ОИСН-16 + 4 мультиплета) —
~25-40k токенов.

Базовый бюджет сессии: 200 000 токенов (default, настраивается).
Это приближение «5-часовой бесплатной сессии Claude»; точное число
зависит от тарифа. CLI флаг `--cost-session-token-budget` позволяет
переопределить.

Если пользователь передал `--cost-tokens N` явно, эта цифра
используется как override (и в Markdown, и в HTML); поэтапная
таблица в Markdown всё равно выводится из авто-оценки.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any


# ─── базовый бюджет 5-часовой сессии ──────────────────────────────────
#
# Анcropic Claude Free tier: типично ~190-220k токенов на 5h сессию
# в зависимости от модели. Берём 200_000 как разумную середину;
# пользователь может переопределить через CLI флаг.
DEFAULT_SESSION_TOKEN_BUDGET = 200_000


# ─── каталог этапов pipeline'а ────────────────────────────────────────
#
# Каждый этап описан фиксированной базовой стоимостью + множителем
# сложности (читается из StagedAnalysisResult полей). Числа эмпирические;
# при значительном изменении pipeline их можно подкрутить.

@dataclass(frozen=True)
class StageCostEstimate:
    """Оценка стоимости одного этапа pipeline."""
    stage_id: str           # "step_1", "step_5α", ...
    stage_name_ru: str      # «Шаг 1: парсинг файла + метаданные»
    tokens_baseline: int    # фиксированная стоимость этапа
    tokens_complexity: int  # надбавка за сложность (n_peaks, n_nuclides, ...)
    tokens_total: int       # baseline + complexity
    why: str                # объяснение надбавки (для прозрачности)

    def to_dict(self) -> dict:
        return {
            "stage_id": self.stage_id,
            "stage_name_ru": self.stage_name_ru,
            "tokens_baseline": int(self.tokens_baseline),
            "tokens_complexity": int(self.tokens_complexity),
            "tokens_total": int(self.tokens_total),
            "why": self.why,
        }


@dataclass(frozen=True)
class CostEstimate:
    """Полная оценка стоимости анализа."""
    tokens_total: int
    session_token_budget: int
    session_pct: float
    by_stage: List[StageCostEstimate] = field(default_factory=list)
    override_used: bool = False   # True если использовался CLI --cost-tokens
    detail: str = ""
    # P0-8: actual Claude output-token count for this analysis run.
    # Sourced from CLI --cost-tokens output_tokens arg when available;
    # defaults to 0 (unknown). Triggers COST_HIGH_OUTPUT_TOKENS alarm at >=20k.
    claude_output_tokens: int = 0

    def to_dict(self) -> dict:
        return {
            "tokens_total": int(self.tokens_total),
            "session_token_budget": int(self.session_token_budget),
            "session_pct": float(self.session_pct),
            "session_pct_formatted": f"{self.session_pct:.1f}%",
            "by_stage": [s.to_dict() for s in self.by_stage],
            "override_used": bool(self.override_used),
            "detail": self.detail,
            "claude_output_tokens": int(self.claude_output_tokens),
        }


# ─── per-stage оценщики ───────────────────────────────────────────────

def _step1_parsing(result) -> StageCostEstimate:
    """Шаг 1: парсинг файла + чтение метаданных."""
    n_ch = int(getattr(result.spec, "n_channels", 0) or 0)
    extra = 0
    why = "базовая стоимость чтения .spe заголовка + двоичных данных"
    if n_ch > 8192:
        extra = 200
        why += f"; +200 на {n_ch} каналов > 8192"
    return StageCostEstimate(
        stage_id="step_1",
        stage_name_ru="Шаг 1: чтение файла + метаданные",
        tokens_baseline=400,
        tokens_complexity=extra,
        tokens_total=400 + extra,
        why=why,
    )


def _step2_environment(result) -> StageCostEstimate:
    """Шаг 2: классификация условий измерения (фон / источник)."""
    n_extras = len(getattr(result.spec, "extras", {}) or {})
    extra = min(200, n_extras * 5)
    return StageCostEstimate(
        stage_id="step_2",
        stage_name_ru="Шаг 2: классификация условий измерения (фон / источник)",
        tokens_baseline=600,
        tokens_complexity=extra,
        tokens_total=600 + extra,
        why=f"+{extra} на чтение {n_extras} полей-расширений",
    )


def _step3_calibration(result) -> StageCostEstimate:
    """Шаги 3–5: калибровка энергии + модель ПШПВ + подбор якорей."""
    recal_diag = getattr(result, "recalibration_diag", {}) or {}
    recal_extra = 600 if recal_diag.get("applied") else 0
    why = "оценка модели ПШПВ + сидинг калибровки E(N)"
    if recal_extra:
        why += "; +600 на пересчёт калибровки (F-87 шаг 5β)"
    return StageCostEstimate(
        stage_id="step_3_4_5",
        stage_name_ru="Шаги 3–5: калибровка E(N), модель ПШПВ, подбор якорей",
        tokens_baseline=1200,
        tokens_complexity=recal_extra,
        tokens_total=1200 + recal_extra,
        why=why,
    )


def _step5_priority(result) -> StageCostEstimate:
    """Шаг 5α/γ — экспресс-якоря + доминантная цепочка + 7-линейная проверка."""
    pf = getattr(result, "priority_findings", []) or []
    n_pf = len(pf)
    cd = getattr(result, "chain_dominance", None)
    extra = n_pf * 80
    if cd is not None and (cd.th232 or cd.u238):
        extra += 200
    return StageCostEstimate(
        stage_id="step_5_express",
        stage_name_ru="Шаг 5α/γ: экспресс-якоря + доминантная цепочка (F-88)",
        tokens_baseline=900,
        tokens_complexity=extra,
        tokens_total=900 + extra,
        why=(f"+{n_pf*80} на {n_pf} приоритетных совпадений"
             + ("; +200 на анализ доминантной цепочки"
                if (cd and (cd.th232 or cd.u238)) else "")),
    )


def _step6_peak_search(result) -> StageCostEstimate:
    """Шаг 6: поиск пиков (Марискотти / свёртка / сравнение)."""
    n_peaks = len(getattr(result, "peaks", []) or [])
    method = getattr(result, "peak_search_method", "mariscotti")
    method_ru = {
        "mariscotti": "Марискотти",
        "convolution": "свёртка",
        "compare": "сравнение",
    }.get(method, "Марискотти")
    extra = n_peaks * 40
    method_extra = 0
    why = f"+{n_peaks*40} на {n_peaks} найденных пиков"
    if method in ("convolution", "compare"):
        method_extra = 400 if method == "convolution" else 800
        why += f"; +{method_extra} на метод «{method_ru}» (F-129)"
    return StageCostEstimate(
        stage_id="step_6_peak_search",
        stage_name_ru="Шаг 6: поиск пиков (Марискотти / свёртка)",
        tokens_baseline=600,
        tokens_complexity=extra + method_extra,
        tokens_total=600 + extra + method_extra,
        why=why,
    )


def _step7_identification(result) -> StageCostEstimate:
    """Шаг 7: идентификация нуклидов (этапы 1/2/3)."""
    stages = getattr(result, "stages", []) or []
    n_stages = len(stages)
    n_cand = sum(len(s.candidates_considered) for s in stages)
    n_det = sum(len(getattr(s, "detected", [])) for s in stages)
    extra = n_stages * 300 + n_cand * 25 + n_det * 80
    return StageCostEstimate(
        stage_id="step_7_identification",
        stage_name_ru="Шаг 7: идентификация нуклидов (этапы 1..3)",
        tokens_baseline=1500,
        tokens_complexity=extra,
        tokens_total=1500 + extra,
        why=(f"+{n_stages*300} на {n_stages} этап(ов); "
             f"+{n_cand*25} на {n_cand} кандидатов; "
             f"+{n_det*80} на {n_det} подтверждений"),
    )


def _step8_deconvolution(result) -> StageCostEstimate:
    """Шаг 8: разложение мультиплетов (F-117/F-118/F-121 + F-126)."""
    dr = getattr(result, "deconvolution_results", None) or []
    n_clusters = len(dr)
    nl_refines = sum(1 for d in dr
                     if "nl_refine" in (getattr(d, "method", "") or ""))
    extra = n_clusters * 500 + nl_refines * 400
    return StageCostEstimate(
        stage_id="step_8_deconvolution",
        stage_name_ru="Шаг 8: разложение мультиплетов (F-117 / F-126)",
        tokens_baseline=800,
        tokens_complexity=extra,
        tokens_total=800 + extra,
        why=(f"+{n_clusters*500} на {n_clusters} кластеров"
             + (f"; +{nl_refines*400} на {nl_refines} нелинейных уточнений"
                if nl_refines else "")),
    )


def _step9_activities(result) -> StageCostEstimate:
    """Шаг 9: активности + ISO 11929 МДА + поправки."""
    acts = getattr(result, "activities", None) or []
    mda = getattr(result, "mda_per_line", None) or {}
    n_acts = len(acts)
    n_mda = len(mda)
    tcs_n = sum(1 for a in acts
                if getattr(a, "coincidence_correction_applied", False))
    sa = getattr(result, "spec", None)
    has_self_att = bool(
        sa and (sa.extras or {}).get("lsrm_sample_density_g_cm3")
    )
    extra = n_acts * 150 + n_mda * 20 + tcs_n * 80
    if has_self_att:
        extra += 200
    why = (f"+{n_acts*150} на {n_acts} нуклидов; "
           f"+{n_mda*20} на {n_mda} линий МДА")
    if tcs_n:
        why += f"; +{tcs_n*80} на каскадные поправки (F-128 Bi-212)"
    if has_self_att:
        why += "; +200 на самопоглощение F-122"
    return StageCostEstimate(
        stage_id="step_9_activities",
        stage_name_ru="Шаг 9: активности + МДА (ISO 11929) + поправки",
        tokens_baseline=1200,
        tokens_complexity=extra,
        tokens_total=1200 + extra,
        why=why,
    )


def _step10_residuals(result) -> StageCostEstimate:
    """Шаг 10: классификация остаточных пиков."""
    rc = getattr(result, "residual_classifications", []) or []
    unmatched = getattr(result, "final_unmatched", []) or []
    n_rc = len(rc)
    n_un = len(unmatched)
    extra = n_rc * 40 + n_un * 30
    return StageCostEstimate(
        stage_id="step_10_residuals",
        # F-386.1 / v1.18.28 (Agent B) — «вылет», не «ускользание» (F-386).
        stage_name_ru=("Шаг 10: классификация остаточных пиков "
                       "(рентген.флуор. / вылет / сумма / край)"),
        tokens_baseline=500,
        tokens_complexity=extra,
        tokens_total=500 + extra,
        why=(f"+{n_rc*40} на {n_rc} классифицированных; "
             f"+{n_un*30} на {n_un} несопоставленных"),
    )


def _step11_reporting(result) -> StageCostEstimate:
    """Шаг 11: сборка отчётов и графиков."""
    n_det = len(getattr(result, "final_detected", []) or [])
    n_peaks = len(getattr(result, "peaks", []) or [])
    extra = n_det * 200 + n_peaks * 30
    return StageCostEstimate(
        stage_id="step_11_reporting",
        stage_name_ru="Шаг 11: сборка отчётов (машино- + текстовый + интерактивный) + графики",
        tokens_baseline=2200,
        tokens_complexity=extra,
        tokens_total=2200 + extra,
        why=(f"+{n_det*200} на {n_det} нуклидов в таблицах; "
             f"+{n_peaks*30} на {n_peaks} пиков на графиках"),
    )


# ─── публичное API ────────────────────────────────────────────────────

def estimate_cost_per_stage(result) -> List[StageCostEstimate]:
    """Вернуть список оценок по 9 этапам Step 1..11 pipeline."""
    return [
        _step1_parsing(result),
        _step2_environment(result),
        _step3_calibration(result),
        _step5_priority(result),
        _step6_peak_search(result),
        _step7_identification(result),
        _step8_deconvolution(result),
        _step9_activities(result),
        _step10_residuals(result),
        _step11_reporting(result),
    ]


def estimate_total_cost(
    result,
    *,
    session_token_budget: int = DEFAULT_SESSION_TOKEN_BUDGET,
    cost_tokens_override: Optional[int] = None,
    detail_override: Optional[str] = None,
) -> CostEstimate:
    """Полная оценка стоимости анализа: per-stage + итог + % сессии.

    Parameters
    ----------
    result : StagedAnalysisResult
        Готовый результат анализа.
    session_token_budget : int
        Бюджет 5-часовой сессии в токенах (default 200_000).
    cost_tokens_override : int or None
        Если задано, используется как итог вместо авто-суммы по этапам.
        Per-stage таблица всё равно строится из авто-оценки — это
        прозрачная разбивка, на случай если пользователь хочет
        свериться с реальностью.
    detail_override : str or None
        Свободная строка описания (для CLI flag `--cost-detail`).

    Returns
    -------
    CostEstimate
        Полный объект с per-stage списком, итогом, % сессии.
    """
    per_stage = estimate_cost_per_stage(result)
    auto_total = sum(s.tokens_total for s in per_stage)
    if cost_tokens_override is not None and cost_tokens_override > 0:
        total = int(cost_tokens_override)
        override = True
    else:
        total = auto_total
        override = False
    pct = (100.0 * total / max(1, session_token_budget))
    detail = detail_override or (
        "Авто-оценка (F-132): итог = сумма по этапам. "
        "Поэтапная разбивка приведена в текстовом отчёте, "
        "раздел «Оценка стоимости анализа»."
    )
    return CostEstimate(
        tokens_total=total,
        session_token_budget=int(session_token_budget),
        session_pct=float(pct),
        by_stage=per_stage,
        override_used=override,
        detail=detail,
    )


__all__ = [
    "StageCostEstimate",
    "CostEstimate",
    "estimate_cost_per_stage",
    "estimate_total_cost",
    "DEFAULT_SESSION_TOKEN_BUDGET",
]
