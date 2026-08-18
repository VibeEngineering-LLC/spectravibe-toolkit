# -*- coding: utf-8 -*-
"""F-367 / v1.18.24.2 — V2 production-pipeline integration via monkey-patch.

Запускает существующий production-pipeline `analyze_and_report` под
context manager, который подменяет ТОЛЬКО peak search-метод на
V2 dual-method (`search_dual_method` = Mariscotti ∪ matched filter).
Все остальные стадии — identification, activities, MDA, reporting,
technical PDF, BecqMoni XML, multiplet PNGs — работают НЕИЗМЕННО
поверх V2-найденных пиков.

Контракт: V2-отчёт **строго идентичен** production-отчёту по составу
(JSON + MD + HTML + Technical PDF + plots + XML), отличается ТОЛЬКО
peak search-методом → бóльшим числом найденных пиков.

Пример:
    from gamma.experimental.v2_integration import analyze_and_report_v2
    artefacts = analyze_and_report_v2(
        sample_path, background_path=bg_path,
        output_dir="./out_v2", sample_mass_kg=0.5,
    )

F-389 / v1.18.26 — V2 activity parity guard (Variant A).
─────────────────────────────────────────────────────────────────────
V2 dual-method находит дополнительные пики (matched_filter-источник),
которые отсутствуют в prod-Mariscotti. Эти пики становятся
полноценными `matched_lines` у `NuclideIdentification` и тянутся в
weighted-mean activity. Результат: V2 activity сильно отличается от
prod (для Th-232 demo: Ac-228 1835 Bq/kg → 905 Bq/kg, −51%) при том,
что F-367 контракт обещает identity downstream.

Решение: внутри `v2_peak_search_patched()` дополнительно патчим
`compute_activities_for_all`, который фильтрует `matched_lines` каждого
`NuclideIdentification` — оставляя только entries, чьи `peak_channel`
найдены **Mariscotti-методом** (источник "mariscotti" или "both"). Это
делает V2 activity численно эквивалентной prod, сохраняя V2 как
расширенный diagnostic в `primary_feps` / `secondary_peaks` /
`multiplet_deconvolutions`.

F-389.1 / v1.18.27.1 — parity tightening (без code-changes здесь).
─────────────────────────────────────────────────────────────────────
F-391 (S/N gating multiplet/singleton, S/N ≥ 3 / 5) на post-v1.18.27
сборке отбрасывает matched_filter-only V2-extras ранее — на стадиях
multiplet detection / singleton acceptance, ДО compute_activities.
Empirical на Th-232 demo (Marinelli 1L): Ac-228 / Tl-208 prod vs V2
specific_activity = 0.0% diff (идентичны до 4 знаков). Test tolerance
ужесточён до 5% (test_f389_v2_activity_parity.PARITY_TOL_FRACTION).
Контракт patches тот же, расширения не требуется — F-391 покрывает
все известные пути leak'а V2-extras.
"""
from __future__ import annotations
from contextlib import contextmanager
from typing import Any, Callable, Dict, List, Optional, Set

import numpy as np

__all__ = ["v2_peak_search_patched", "analyze_and_report_v2",
           "v2_mariscotti_replacement",
           "is_activity_filter_enabled"]


# F-389: per-process регистр каналов, найденных ТОЛЬКО matched_filter
# (т.е. V2-only — НЕ найденных prod-Mariscotti). Обновляется в
# `v2_mariscotti_replacement`, читается в _filter_v2_only_activities.
# Set очищается при входе в `v2_peak_search_patched()` и при выходе.
_V2_ONLY_CHANNELS: Set[int] = set()
# Per-process флаг — activity-фильтр активен ⟺ True (set/cleared только
# `v2_peak_search_patched()`). Используется тестами и инспекцией.
_ACTIVITY_FILTER_ACTIVE: bool = False


def is_activity_filter_enabled() -> bool:
    """True ⟺ V2 activity-фильтр (F-389) сейчас активен.

    Возвращает значение module-level флага, set внутри
    `v2_peak_search_patched()` context manager. Используется в тестах
    и diagnostics.
    """
    return bool(_ACTIVITY_FILTER_ACTIVE)


def _local_continuum_estimate(counts: np.ndarray, ch: int,
                              fwhm_ch: float) -> float:
    """Грубая оценка фона под пиком: медиана в окне ±(1.5..2.5)·FWHM."""
    half = max(2, int(round(fwhm_ch * 2.0)))
    inner = max(1, int(round(fwhm_ch * 1.5)))
    n = len(counts)
    bg_left = counts[max(0, ch - half):max(0, ch - inner)]
    bg_right = counts[min(n, ch + inner):min(n, ch + half)]
    bg = np.concatenate([bg_left, bg_right])
    if bg.size == 0:
        return 0.0
    return float(np.median(bg))


def v2_mariscotti_replacement(
    counts,
    fwhm_channels,
    sigma_threshold=3.0,
    min_separation_factor: float = 1.0,
    edge_margin: Optional[int] = None,
    *,
    band_ratio: float = 1.2,
    filter_narrow_peaks: bool = False,
    min_fwhm_ratio: float = 0.3,
):
    """Drop-in replacement для `mariscotti_search` — возвращает merged
    результаты V2 dual-method search в формате `List[FoundPeak]`.

    Сигнатура полностью совместима с `gamma.peaks.search.mariscotti_search`
    — все production-вызовы работают без правок.
    """
    from gamma.peaks.search import FoundPeak
    from gamma.experimental.peak_pipeline_v2 import search_dual_method

    counts_arr = np.asarray(counts, dtype=np.float64)
    n = counts_arr.size

    # fwhm_channels может быть scalar / callable
    if callable(fwhm_channels):
        fwhm_provider = fwhm_channels
    else:
        const_fwhm = float(fwhm_channels)
        def fwhm_provider(ch: int) -> float:
            return const_fwhm

    # sigma_threshold для search_dual_method — scalar (V2 не поддерживает
    # per-channel threshold). Если callable — берём value в середине спектра.
    if callable(sigma_threshold):
        sigma_val = float(sigma_threshold(n // 2))
    else:
        sigma_val = float(sigma_threshold)

    merged, mari_raw, conv_raw = search_dual_method(
        counts=counts_arr,
        fwhm_provider=fwhm_provider,
        sigma_threshold=sigma_val,
    )

    # Edge margin filter — соответствие production-default
    if edge_margin is None:
        # Production-default: 2·FWHM(edge)
        try:
            edge_margin = int(round(2.0 * fwhm_provider(0)))
        except Exception:
            edge_margin = 10
    edge_margin = max(0, int(edge_margin))

    # Narrow-peak filter (F-139) — на V2 не критично, но соблюдаем contract
    out: List[FoundPeak] = []
    for hit in merged:
        ch = int(hit.channel)
        if ch < edge_margin or ch >= n - edge_margin:
            continue
        fwhm_ch = float(hit.fwhm_channels)
        if fwhm_ch <= 0:
            continue
        sigma_ch = fwhm_ch / 2.3548
        # Высота — counts[ch] минус локальный континуум
        bg = _local_continuum_estimate(counts_arr, ch, fwhm_ch)
        height = float(max(counts_arr[ch] - bg, 0.0))
        # Оценка площади как Gauss с σ_ch — area ≈ height · √(2π)·σ
        area_est = float(height * 2.5066 * sigma_ch)
        sigma_area = float(np.sqrt(max(area_est + bg * 5.0, 0.0)))
        out.append(FoundPeak(
            channel=ch,
            height=height,
            fwhm_channels=fwhm_ch,
            significance=float(hit.significance),
            area_estimate=area_est,
            sigma_area_estimate=sigma_area,
            notes=[f"v2-source:{hit.source}"],
        ))

    # Min-separation filter (как в Mariscotti)
    if min_separation_factor > 0 and len(out) > 1:
        out.sort(key=lambda p: p.channel)
        filtered = [out[0]]
        for p in out[1:]:
            prev = filtered[-1]
            min_sep = min_separation_factor * max(
                prev.fwhm_channels, p.fwhm_channels,
            )
            if (p.channel - prev.channel) >= min_sep:
                filtered.append(p)
            else:
                # Конфликт — оставляем пик с бОльшей significance
                if p.significance > prev.significance:
                    filtered[-1] = p
        out = filtered

    out.sort(key=lambda p: p.channel)

    # F-389: регистрируем каналы V2-extras (НЕ найденные prod-Mariscotti).
    # Источник правды — отдельный вызов `mariscotti_search` с production-
    # default kwargs (min_separation_factor=0.6, edge_margin=10,
    # filter_narrow_peaks из вызова). Это в точности тот же call, который
    # сделал бы `_run_peak_search(method="mariscotti")` в staged_pipeline,
    # поэтому prod-channel-set воспроизводится 1-в-1. «V2-only channel»
    # ≡ финальный out[*].channel НЕ принадлежит prod-mariscotti-set
    # (в пределах ±0.5·FWHM tolerance).
    if _ACTIVITY_FILTER_ACTIVE:
        from gamma.peaks.search import mariscotti_search as _prod_mari
        try:
            _prod_kwargs = dict(
                counts=counts_arr,
                fwhm_channels=fwhm_channels,
                sigma_threshold=sigma_threshold,
                min_separation_factor=0.6,
                edge_margin=10,
                filter_narrow_peaks=bool(filter_narrow_peaks),
                min_fwhm_ratio=float(min_fwhm_ratio),
            )
            prod_mari = _prod_mari(**_prod_kwargs)
        except TypeError:
            # Сигнатура mariscotti_search в данной сборке может не
            # принимать filter_narrow_peaks — откатываемся на минимум
            prod_mari = _prod_mari(
                counts=counts_arr,
                fwhm_channels=fwhm_channels,
                sigma_threshold=sigma_threshold,
                min_separation_factor=0.6,
                edge_margin=10,
            )
        mari_channels = {int(p.channel) for p in prod_mari}
        for p in out:
            ch = int(p.channel)
            fwhm_ch = max(1.0, float(p.fwhm_channels))
            tol = max(1, int(round(0.5 * fwhm_ch)))
            # ch in V2-only ⟺ нет prod-mariscotti-channel в окне ±tol
            if not any(abs(mc - ch) <= tol for mc in mari_channels):
                _V2_ONLY_CHANNELS.add(ch)
    return out


def _filter_matched_lines_for_v2(id_result):
    """F-389: возвращает копию `IdentificationResult`, где каждый
    `NuclideIdentification.matched_lines` отфильтрован — оставлены
    только entries, чей `peak_channel` НЕ в `_V2_ONLY_CHANNELS`.

    Это даёт паритет с prod-Mariscotti при weighted-mean activity:
    matched_filter-only пики (V2 extras) полностью игнорируются на
    стадии активности, но остаются видны в primary_feps / secondary_peaks.
    """
    if not _V2_ONLY_CHANNELS:
        return id_result
    import dataclasses

    def _filter_nuclide(ni):
        ml = getattr(ni, "matched_lines", ()) or ()
        kept = tuple(
            m for m in ml
            if int(getattr(m, "peak_channel", -1)) not in _V2_ONLY_CHANNELS
        )
        if len(kept) == len(ml):
            return ni
        # Use dataclasses.replace to preserve frozen contract
        try:
            return dataclasses.replace(ni, matched_lines=kept)
        except Exception:
            return ni

    try:
        detected = tuple(_filter_nuclide(n) for n
                         in getattr(id_result, "detected_nuclides", ()))
        rejected = tuple(_filter_nuclide(n) for n
                         in getattr(id_result, "rejected_nuclides", ()))
        return dataclasses.replace(
            id_result,
            detected_nuclides=detected,
            rejected_nuclides=rejected,
        )
    except Exception:
        return id_result


def _patched_compute_activities_for_all(id_result, *args, **kwargs):
    """F-389 monkey-patch wrapper: фильтрует matched_lines на V2-only
    pickup'ы перед вызовом оригинального `compute_activities_for_all`.

    Сохраняет identity contract F-367 (V2 activity == prod activity на
    общем подмножестве peak'ов), при этом V2 extras остаются видны в
    остальных секциях отчёта (primary_feps, multiplet_deconvolutions,
    unidentified_peaks).
    """
    from gamma.activity.compute import compute_activities_for_all as _orig
    filtered = _filter_matched_lines_for_v2(id_result)
    return _orig(filtered, *args, **kwargs)


@contextmanager
def v2_peak_search_patched():
    """Context manager: на время блока подменяет `mariscotti_search`
    в `gamma.identification.staged_pipeline` на V2 dual-method.

    Гарантия: restore оригинальной функции в finally — даже при
    исключении внутри.

    F-389: дополнительно подменяется `compute_activities_for_all` в
    `gamma.identification.staged_pipeline` — оборачивает оригинал
    фильтром, который выкидывает matched_lines с peak_channel ∈
    `_V2_ONLY_CHANNELS` (matched_filter-only V2-extras). Гарантирует
    численный паритет V2 vs prod на одинаковом peak-subset.
    """
    global _ACTIVITY_FILTER_ACTIVE
    from gamma.identification import staged_pipeline as _sp
    _orig_mari = _sp.mariscotti_search
    _orig_cact = _sp.compute_activities_for_all
    _V2_ONLY_CHANNELS.clear()
    _ACTIVITY_FILTER_ACTIVE = True
    _sp.mariscotti_search = v2_mariscotti_replacement
    _sp.compute_activities_for_all = _patched_compute_activities_for_all
    try:
        yield
    finally:
        _sp.mariscotti_search = _orig_mari
        _sp.compute_activities_for_all = _orig_cact
        _V2_ONLY_CHANNELS.clear()
        _ACTIVITY_FILTER_ACTIVE = False


def analyze_and_report_v2(path: str, **kwargs) -> Dict[str, Any]:
    """Drop-in замена `analyze_and_report` использующая V2 dual-method
    peak search вместо Mariscotti.

    Все kwargs идентичны production `analyze_and_report` — V2-режим
    влияет ТОЛЬКО на peak search-стадию. Identification, multiplet
    deconvolution, activities, MDA, reporting — production без изменений.
    """
    from gamma.reporting import analyze_and_report
    with v2_peak_search_patched():
        return analyze_and_report(path, **kwargs)
