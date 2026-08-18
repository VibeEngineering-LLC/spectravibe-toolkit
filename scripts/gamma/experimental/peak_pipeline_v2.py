"""
F-354 / v1.18.24.0 — экспериментальный peak-pipeline v2.

Запускается ПАРАЛЛЕЛЬНО с production analyze_lsrm_spe; **не заменяет**
production-логику. Цель — собрать сравнительные данные по двум
направлениям, обсуждённым в сессии:

  1. **Поиск пиков двумя методами одновременно** (Mariscotti d²Gauss +
     matched filter Gauss-конволюция). Объединение результатов даёт
     полный список кандидатов; различия записываются в diff-отчёт.
     См. peaks/search.py + peaks/convolution_search.py.

  2. **Автоматическое обнаружение мультиплетов** (без жёстких таблиц
     TH232_FORCED_CLUSTERS / RA226_FORCED_CLUSTERS из deconvolve.py).
     Для каждого набора близких библиотечных линий (попадающих в одно
     окно ROI ≈ 2·FWHM·n_components) запускается coupled_intensity_fit
     с component-list, собранным из реальной библиотеки нуклидов.
     См. peaks/coupled_multiplet.py.

  3. **Сравнение с production**: ComparisonReport содержит diff по числу
     пиков, разделённым мультиплетам, площадям компонент,
     обнаруженным/пропущенным линиям.

Контракт:
  * НЕ импортируется из gamma.identification.staged_pipeline.
  * Вызывается только из prompt-driven scripts / tests / TBD CLI флага.
  * F-115 anonymization применяется через анализ через analyze_lsrm_spe
    при сравнении — здесь raw spectrum только в памяти.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np


# ──────────────────────────────────────────────────────────────────
# Result dataclasses
# ──────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class PeakHit:
    """Объединённый pick от двух методов.

    Если пик найден обоими — `source = "both"`, иначе `"mariscotti"` или
    `"matched_filter"`. Канал и FWHM берутся от того метода, который
    выдал бОльшую `significance` (S/σ).
    """
    channel: int
    energy_keV: float
    fwhm_channels: float
    significance: float
    source: str               # "mariscotti" | "matched_filter" | "both"
    mari_significance: Optional[float] = None
    conv_significance: Optional[float] = None


@dataclass(frozen=True)
class MultipletClusterCandidate:
    """Подозреваемый мультиплет: набор библиотечных линий, попадающих
    в одно общее окно ROI ≈ 2·FWHM·n_components.

    Issue #37 / v1.18.31+: `phantom_components` хранит библиотечные линии,
    которые попали в CC, но были вытеснены top-K cap'ом. Они НЕ участвуют
    в fit'е (downstream использует только `components`), но сохраняются
    для evidence/диагностики, и пост-merge dedup (см. ниже) использует их
    при разрешении конфликтов «одна линия в двух кластерах».
    """
    cluster_id: str           # авто-генерированный, "auto_M{N}_cc{k}"
    E_lo_keV: float
    E_hi_keV: float
    components: Tuple[Tuple[str, float, float, str], ...]
    # (nuclide, E_keV, I_gamma_pct, group)
    n_components: int
    detection_reason: str     # человеко-читаемое объяснение
    found_peaks_in_roi: Tuple[int, ...] = ()   # channels из search
    # Issue #37 — компоненты, вытесненные top-K cap'ом. По UX отображаются
    # в evidence/HTML note, но НЕ fit'ятся (downstream берёт `components`).
    phantom_components: Tuple[Tuple[str, float, float, str], ...] = ()


@dataclass
class PipelineV2Result:
    """Полный output v2-pipeline."""
    found_peaks: List[PeakHit]
    mariscotti_peaks: List[Any] = field(default_factory=list)     # FoundPeak
    matched_filter_peaks: List[Any] = field(default_factory=list)
    multiplet_candidates: List[MultipletClusterCandidate] = field(default_factory=list)
    coupled_fits: List[Any] = field(default_factory=list)         # CoupledFitResult
    notes: List[str] = field(default_factory=list)
    detector_class: str = "NaI"                                   # из identification


@dataclass
class ComparisonReport:
    """Diff между production и v2."""
    # Production-side
    production_n_peaks: int = 0
    production_n_multiplets: int = 0
    production_multiplet_ids: List[str] = field(default_factory=list)

    # v2-side
    v2_n_peaks: int = 0
    v2_n_multiplets: int = 0
    v2_multiplet_ids: List[str] = field(default_factory=list)

    # Diff
    peaks_only_in_v2: List[Tuple[int, float]] = field(default_factory=list)
    peaks_only_in_production: List[Tuple[int, float]] = field(default_factory=list)
    multiplets_only_in_v2: List[str] = field(default_factory=list)
    multiplets_only_in_production: List[str] = field(default_factory=list)

    # Free-form
    notes: List[str] = field(default_factory=list)


# ──────────────────────────────────────────────────────────────────
# 1. Dual-method search
# ──────────────────────────────────────────────────────────────────

def search_dual_method(
    counts: np.ndarray,
    fwhm_provider: Callable[[int], float],
    *,
    sigma_threshold: float = 3.0,
    merge_tolerance_channels: float = 1.5,
    channel_to_energy: Optional[Callable[[int], float]] = None,
) -> Tuple[List[PeakHit], List[Any], List[Any]]:
    """Запускает Mariscotti + matched filter, объединяет результаты.

    Returns
    -------
    (merged, mariscotti_raw, matched_filter_raw)
        merged — список PeakHit с тегом source ∈ {"mariscotti", "matched_filter", "both"}
        mariscotti_raw, matched_filter_raw — оригинальные FoundPeak списки от каждого метода
    """
    from gamma.peaks.search import mariscotti_search
    from gamma.peaks.convolution_search import (
        convolution_peak_search, compare_peak_methods,
    )

    mari = mariscotti_search(counts, fwhm_channels=fwhm_provider,
                              sigma_threshold=sigma_threshold)
    conv = convolution_peak_search(counts, fwhm_channels=fwhm_provider,
                                    sigma_threshold=sigma_threshold)

    cmp = compare_peak_methods(mari, conv,
                                tolerance_channels=merge_tolerance_channels)

    def _to_E(ch: int) -> float:
        if channel_to_energy is None:
            return float(ch)
        return float(channel_to_energy(int(ch)))

    merged: List[PeakHit] = []
    # «оба» — берём более сильный (по significance)
    for pa, pb in cmp["agreed"]:
        if pa.significance >= pb.significance:
            top = pa
        else:
            top = pb
        merged.append(PeakHit(
            channel=int(top.channel),
            energy_keV=_to_E(top.channel),
            fwhm_channels=float(top.fwhm_channels),
            significance=float(top.significance),
            source="both",
            mari_significance=float(pa.significance),
            conv_significance=float(pb.significance),
        ))
    # только Mariscotti
    for p in cmp["a_only"]:
        merged.append(PeakHit(
            channel=int(p.channel),
            energy_keV=_to_E(p.channel),
            fwhm_channels=float(p.fwhm_channels),
            significance=float(p.significance),
            source="mariscotti",
            mari_significance=float(p.significance),
            conv_significance=None,
        ))
    # только matched filter
    for p in cmp["b_only"]:
        merged.append(PeakHit(
            channel=int(p.channel),
            energy_keV=_to_E(p.channel),
            fwhm_channels=float(p.fwhm_channels),
            significance=float(p.significance),
            source="matched_filter",
            mari_significance=None,
            conv_significance=float(p.significance),
        ))
    merged.sort(key=lambda h: h.channel)
    return merged, mari, conv


# ──────────────────────────────────────────────────────────────────
# 2. Multiplet auto-detection (без FORCED_CLUSTERS)
# ──────────────────────────────────────────────────────────────────

# Минимальный библиотечный набор для популярных цепочек.
# В production источник — `data/nuclides.json`. Здесь — sample.
_DEFAULT_LIB: Dict[str, List[Tuple[float, float]]] = {
    # nuclide → [(E_keV, I_γ_pct), ...] — ключевые линии для Gamma-1S
    "Ac-228":  [(338.32, 11.3), (911.20, 25.8), (964.77, 4.99),
                (968.97, 15.8), (1588.20, 3.22), (1630.6, 1.60)],
    "Tl-208":  [(277.4, 6.6), (510.77, 22.6), (583.19, 30.6),
                (860.6, 12.5), (2614.51, 99.0)],
    "Pb-212":  [(238.63, 43.6), (300.09, 3.30)],
    "Bi-212":  [(727.33, 6.65), (1620.50, 1.49)],
    "Pb-214":  [(295.22, 18.42), (351.93, 35.6)],
    "Bi-214":  [(609.31, 45.49), (665.45, 1.531), (768.36, 4.892),
                (806.17, 1.262), (1120.29, 14.92), (1238.11, 5.834),
                (1764.49, 15.31)],
    "K-40":    [(1460.82, 10.66)],
    "Cs-137":  [(661.66, 85.10)],
}


def detect_multiplet_clusters(
    found_peaks: Sequence[PeakHit],
    fwhm_keV_at: Callable[[float], float],
    *,
    library: Optional[Dict[str, List[Tuple[float, float]]]] = None,
    chain_filter: Optional[Sequence[str]] = None,
    overlap_n_fwhm: float = 1.2,
    roi_extend_fwhm: float = 2.5,
    expand_to_display_window: bool = True,
    display_window_fwhm: float = 3.0,
    # F-387 / v1.18.26 → F-387.1 / v1.18.26.1 — Rayleigh-CC unresolved
    # criterion (зеркально find_multiplet_regions в gamma.peaks.deconvolve).
    # Pair (a, b) unresolved ⟺ |ΔE| < factor · FWHM_avg(a, b), где
    # FWHM_avg = (FWHM_a + FWHM_b)/2. **Семантика изменена в v1.18.26.1**:
    # было factor·FWHM_min (один min глобально, keep-cluster-as-monolith),
    # стало factor·FWHM_avg (среднее пары) с CC-split + top-K cap.
    # Default 1.1 (NaI 63×63 calibration override v1.18.26.1) — Rayleigh
    # base 1.0 расширен +10% для покрытия borderline pairs с валли-фил
    # от интенсивных соседей (Ac-228 911 vs 969 factor=1.01, 463 vs 504
    # factor=1.03). См. `find_multiplet_regions` docstring в deconvolve.py.
    unresolved_separation_fwhm_factor: float = 1.1,
    # F-387.1 — top-K cap по library I_pct в sub-cluster крупнее этого
    # значения. Остальные становятся phantom anchors (в этой v2-функции
    # phantom'ы помечаются как кортеж `(nuc, E, I, group, phantom=True)`
    # — но компоненты v2 — это tuple, поэтому phantom помечается через
    # group="phantom" sentinel? Нет — V2 cluster.components — это
    # `tuple[(nuc, E, I, group), ...]`, без флага. Минимально-инвазивно:
    # `group=""` для phantom — но он уже используется как «independent».
    # Решение: phantom компоненты в V2 **удаляются из cluster.components**
    # (downstream V2 fit принимает все группы как fit'ные). Это упрощает
    # путь и адекватно для V2 (experimental). Для evidence в V2 нет
    # «identification»: cluster.components = только fit'ные.
    max_components_per_cluster: int = 3,
    # F-391 / v1.18.27 — S/N significance gate (зеркально
    # find_multiplet_regions). В V2 семантика проще: фильтруем
    # input found_peaks по PeakHit.significance ≥ threshold ДО seed-pick
    # loop — клiente получит multiplet ТОЛЬКО когда есть реальный signal
    # ≥ 3σ в ROI. F-381 library-anchor enrichment добавляет phantom-линии
    # без peak — они уже не создадут собственный pick_roi. Cluster с
    # одним active pick (singleton) дропается; ≥2 active → multiplet.
    min_significance_snr: float = 3.0,
) -> List[MultipletClusterCandidate]:
    """Автоматическое обнаружение мультиплетов БЕЗ хардкода FORCED_CLUSTERS.

    Алгоритм:
      1. Для каждого pick формируем список библиотечных линий, попадающих
         в окно ±roi_extend_fwhm·FWHM(E_peak).
      2. Если ≥2 библиотечных линий попали в окно → это кандидат на
         мультиплет.
      3. Объединяем пересекающиеся ROI (transitive closure по
         overlap_n_fwhm·FWHM).
      4. Для каждого финального кластера собираем уникальный список
         компонент (nuclide, E, I, group=nuclide).
      5. F-374 — если ``expand_to_display_window=True`` (default),
         каждый кластер расширяется до своего chart-окна
         (±display_window_fwhm·FWHM) и в него втягиваются ВСЕ
         библиотечные линии из ``library``, попадающие внутрь
         расширенного окна. Это гарантирует, что любые
         идентифицированные пики в окне отображения учитываются
         deconvolution-фитом (а не игнорируются как "не overlap").
      6. F-387 — финальный фильтр LSRM Алгоритмические основы §9.4:
         cluster выживает как unresolved multiplet ⟺ ∃ пара компонент
         с |ΔE| < unresolved_separation_fwhm_factor · FWHM_min. Иначе
         все пары разрешимы (выше Вартанов §6 Fig.20 порога) → drop.

    Returns
    -------
    list[MultipletClusterCandidate]
        Кандидаты, отсортированные по E_центр.
    """
    lib = library if library is not None else _DEFAULT_LIB
    if chain_filter is not None:
        chain_filter = set(chain_filter)
        lib = {k: v for k, v in lib.items() if k in chain_filter}

    # F-391 / v1.18.27 — S/N significance gate (LSRM-9.4 / Gilmore §5.5).
    # Фильтруем входной список found_peaks по PeakHit.significance ≥
    # threshold ДО seed-pick loop. Это гарантирует, что multiplet'ы
    # формируются ТОЛЬКО вокруг реально измеренных пиков с S/N ≥ 3 —
    # шумовые pick'и (e.g. локальные max на Compton continuum NaI) не
    # порождают artefactных multiplet'ов с library lines, попавших в
    # их окно.
    if min_significance_snr > 0.0:
        found_peaks = [
            h for h in found_peaks
            if float(getattr(h, "significance", 0.0) or 0.0)
            >= min_significance_snr
        ]

    # Шаг 1: для каждого pick найти библиотечные линии в ROI
    pick_rois: List[Dict[str, Any]] = []
    for hit in found_peaks:
        E_peak = hit.energy_keV
        if E_peak <= 0:
            continue
        fwhm = fwhm_keV_at(E_peak)
        roi_lo = E_peak - roi_extend_fwhm * fwhm
        roi_hi = E_peak + roi_extend_fwhm * fwhm
        components_in_roi: List[Tuple[str, float, float, str]] = []
        for nuc, lines in lib.items():
            for E_lib, I_pct in lines:
                if roi_lo <= E_lib <= roi_hi:
                    components_in_roi.append((nuc, E_lib, I_pct, nuc))
        if len(components_in_roi) >= 2:
            pick_rois.append({
                "peak_ch": hit.channel,
                "peak_E": E_peak,
                "roi_lo": roi_lo,
                "roi_hi": roi_hi,
                "components": components_in_roi,
                "fwhm": fwhm,
            })

    if not pick_rois:
        return []

    # Шаг 2: объединение пересекающихся ROI (union-find)
    n = len(pick_rois)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)

    pick_rois.sort(key=lambda r: r["peak_E"])
    for i in range(n):
        for j in range(i + 1, n):
            # Пересечение ROI?
            if pick_rois[i]["roi_hi"] >= pick_rois[j]["roi_lo"]:
                union(i, j)
            else:
                # Слишком далеко — break, ROIs отсортированы
                if pick_rois[j]["peak_E"] - pick_rois[i]["peak_E"] > \
                        overlap_n_fwhm * (pick_rois[i]["fwhm"] + pick_rois[j]["fwhm"]):
                    break

    groups: Dict[int, List[int]] = {}
    for idx in range(n):
        groups.setdefault(find(idx), []).append(idx)

    # Шаг 3: каждый group → MultipletClusterCandidate
    clusters: List[MultipletClusterCandidate] = []
    for k, (root, idxs) in enumerate(sorted(groups.items(),
                                              key=lambda kv: pick_rois[kv[1][0]]["peak_E"])):
        # Слить ROI
        E_lo = min(pick_rois[i]["roi_lo"] for i in idxs)
        E_hi = max(pick_rois[i]["roi_hi"] for i in idxs)
        # Собрать уникальные компоненты
        seen_components: Dict[Tuple[str, float], Tuple[str, float, float, str]] = {}
        for i in idxs:
            for c in pick_rois[i]["components"]:
                key = (c[0], round(c[1], 2))
                if key not in seen_components:
                    seen_components[key] = c

        # F-374 — расширить окно до display-window и втянуть ВСЕ
        # библиотечные линии в нём (не только seed-overlap pick'и). Это
        # покрывает identified пики, которые формально не вошли в seed
        # clustering, но видны в chart-окне мультиплета.
        if expand_to_display_window:
            # Используем FWHM на краях для расширения
            fwhm_lo = fwhm_keV_at(E_lo)
            fwhm_hi = fwhm_keV_at(E_hi)
            E_lo_disp = E_lo - display_window_fwhm * fwhm_lo
            E_hi_disp = E_hi + display_window_fwhm * fwhm_hi
            for nuc, lines in lib.items():
                for E_lib, I_pct in lines:
                    if E_lo_disp <= E_lib <= E_hi_disp:
                        key = (nuc, round(E_lib, 2))
                        if key not in seen_components:
                            seen_components[key] = (nuc, E_lib, I_pct, nuc)
            # Также расширим ROI к display-window-границам, чтобы fit
            # имел достаточно channels вокруг новых компонент.
            E_lo = E_lo_disp
            E_hi = E_hi_disp

        components = tuple(sorted(seen_components.values(), key=lambda c: c[1]))
        # Только если ≥2 компоненты — иначе это одиночка
        if len(components) < 2:
            continue
        reason = (
            f"auto-detected: {len(components)} библиотечные линии в ROI "
            f"[{E_lo:.1f}, {E_hi:.1f}] кэВ"
        )
        found_channels = tuple(pick_rois[i]["peak_ch"] for i in idxs)
        clusters.append(MultipletClusterCandidate(
            cluster_id=f"auto_M{k+1}",
            E_lo_keV=float(E_lo),
            E_hi_keV=float(E_hi),
            components=components,
            n_components=len(components),
            detection_reason=reason,
            found_peaks_in_roi=found_channels,
        ))

    # F-387 / v1.18.26 → F-387.1 / v1.18.26.1 — Rayleigh-CC split +
    # top-K cap (зеркально gamma.peaks.deconvolve.find_multiplet_regions).
    # Алгоритм:
    #   1. Vertices = cl.components, edge(a, b) ⟺
    #      |ΔE| < factor · FWHM_avg(a, b).
    #   2. BFS → connected components.
    #   3. CC размера 1 → isolated singleton (cluster size=1 на выходе;
    #      downstream v2 fit обработает как trivial 1-component fit).
    #      CC размера ≥ 2 → unresolved multiplet sub-cluster.
    #   4. Sub-cluster крупнее `max_components_per_cluster` → top-K по
    #      I_gamma_pct остаются, остальные **удаляются**
    #      (V2 — experimental, без identification evidence).
    #
    # Default factor 1.0 = Rayleigh. References — см. find_multiplet_regions.
    if clusters and unresolved_separation_fwhm_factor > 0.0:
        from collections import deque as _deque
        new_clusters: List[MultipletClusterCandidate] = []
        for cl in clusters:
            comps = list(cl.components)  # mutable
            n_c = len(comps)
            if n_c < 2:
                # F-387.2: singleton input — drop (route в primary_feps
                # через identification_result). Раньше append'ился как
                # «trivial fit» — теперь UX-консистентность с PROD path.
                continue

            # Step 1: adjacency через Rayleigh edges
            fwhms = [float(fwhm_keV_at(float(c[1]))) for c in comps]
            Es = [float(c[1]) for c in comps]
            adj: List[List[int]] = [[] for _ in range(n_c)]
            for i in range(n_c):
                for j in range(i + 1, n_c):
                    fwhm_avg = 0.5 * (fwhms[i] + fwhms[j])
                    if fwhm_avg <= 0:
                        continue
                    if abs(Es[i] - Es[j]) < (
                        unresolved_separation_fwhm_factor * fwhm_avg
                    ):
                        adj[i].append(j)
                        adj[j].append(i)

            # Step 2: BFS → CCs
            visited = [False] * n_c
            cc_groups: List[List[int]] = []
            for start in range(n_c):
                if visited[start]:
                    continue
                queue = _deque([start])
                visited[start] = True
                cc: List[int] = []
                while queue:
                    u = queue.popleft()
                    cc.append(u)
                    for v in adj[u]:
                        if not visited[v]:
                            visited[v] = True
                            queue.append(v)
                cc_groups.append(cc)

            # Step 3+4: каждый CC → новый MultipletClusterCandidate
            # Issue #37 — `_cc{idx}` теперь per-input-cluster (`sub_idx`),
            # не глобальный — IDs стабильны и читаемы.
            sub_idx = 0
            for cc_idxs in cc_groups:
                sub_comps = [comps[i] for i in cc_idxs]
                sub_comps.sort(key=lambda c: c[1])
                # F-387.2 / v1.18.27.1 — CC=1 (singleton) → drop из
                # multiplet output. «1-component multiplet» физически
                # не имеет смысла — это просто peak, его место в
                # primary_feps (через identification_result.matched_lines),
                # не в multiplet_deconvolutions. Зеркально logiке
                # `apply_multiplet_deconvolution` в deconvolve.py.
                # Раньше (F-391 / v1.18.27): singleton'ы оставались
                # как cluster size=1 для downstream «trivial fit», но
                # это создавало UX-аномалию «multiplet с одной линией».
                if len(sub_comps) < 2:
                    continue
                # Top-K cap (issue #37): below-top-K компоненты теперь
                # сохраняются как `phantom_components` вместо полного
                # удаления. Это (1) даёт evidence для UX, (2) позволяет
                # пост-merge dedup'у разрешать конфликты «одна линия в
                # двух CC» (см. ниже).
                phantom_sub: List[Tuple[str, float, float, str]] = []
                if (len(sub_comps) > max_components_per_cluster
                        and max_components_per_cluster > 0):
                    by_intensity = sorted(
                        sub_comps,
                        key=lambda c: float(c[2] or 0.0),
                        reverse=True,
                    )
                    kept = by_intensity[:max_components_per_cluster]
                    phantom_sub = sorted(
                        by_intensity[max_components_per_cluster:],
                        key=lambda c: c[1],
                    )
                    sub_comps = sorted(kept, key=lambda c: c[1])
                E_lo_sub = sub_comps[0][1]
                E_hi_sub = sub_comps[-1][1]
                # Расширим до chart-window (как в шаге F-374 выше)
                fwhm_lo_sub = float(fwhm_keV_at(E_lo_sub))
                fwhm_hi_sub = float(fwhm_keV_at(E_hi_sub))
                E_lo_disp = E_lo_sub - display_window_fwhm * fwhm_lo_sub
                E_hi_disp = E_hi_sub + display_window_fwhm * fwhm_hi_sub
                reason_sub = (
                    f"auto-detected (F-387.1 CC split): "
                    f"{len(sub_comps)} компоненты в ROI "
                    f"[{E_lo_disp:.1f}, {E_hi_disp:.1f}] кэВ"
                )
                new_clusters.append(MultipletClusterCandidate(
                    cluster_id=cl.cluster_id + f"_cc{sub_idx}",
                    E_lo_keV=float(E_lo_disp),
                    E_hi_keV=float(E_hi_disp),
                    components=tuple(sub_comps),
                    n_components=len(sub_comps),
                    detection_reason=reason_sub,
                    found_peaks_in_roi=cl.found_peaks_in_roi,
                    phantom_components=tuple(phantom_sub),
                ))
                sub_idx += 1
        clusters = new_clusters

    # Issue #37 / v1.18.31+ — POST-MERGE DEDUP PASS.
    # Симптом: F-374 display-window expansion (lines 377-393) расширяет
    # ROI каждого union-find кластера на ±display_window_fwhm·FWHM, что
    # для близких mega-кластеров приводит к ПЕРЕСЕЧЕНИЮ окон. Библиотечные
    # линии в overlap-зоне попадают в ОБА кластера → после F-387.1 CC-split
    # одна и та же (nuc, E_keV) линия может встретиться в двух
    # `new_clusters`. Сочетание с top-K cap'ом (см. выше) даёт:
    #   - одна линия (low I_pct) выпадает из «правильного» CC через top-K
    #     и переотрожается в соседнем CC как «mini-cluster overkill»
    #   - симптом #37: «mini-clusters that should fold + real multiplet
    #     broken into ≥2 sub-clusters».
    #
    # Dedup rule (LSRM Алгоритмические основы §9.4 — multiplet uniqueness):
    # для каждой library-line (nuc, round(E, 2)), встречающейся в N≥2
    # clusters, оставляем её только в том кластере, где она ближе к
    # центроиду cluster.components (взвешенному по I_pct). После dedup'а
    # драпаем кластеры, ставшие <2 компонент.
    if clusters and len(clusters) >= 2:
        # Index: (nuc, round(E,2)) → list of (cluster_idx, distance_to_centroid)
        line_occurrences: Dict[Tuple[str, float], List[Tuple[int, float]]] = {}
        for ci, cl in enumerate(clusters):
            if not cl.components:
                continue
            # Intensity-weighted centroid
            total_I = sum(float(c[2] or 0.0) for c in cl.components)
            if total_I > 0:
                centroid = sum(
                    float(c[1]) * float(c[2] or 0.0) for c in cl.components
                ) / total_I
            else:
                centroid = 0.5 * (cl.E_lo_keV + cl.E_hi_keV)
            for c in cl.components:
                key = (c[0], round(float(c[1]), 2))
                dist = abs(float(c[1]) - centroid)
                line_occurrences.setdefault(key, []).append((ci, dist))

        # For lines with ≥2 occurrences, keep in closest-centroid cluster
        # only — remove from others.
        to_remove: Dict[int, set] = {}  # cluster_idx → set of keys to drop
        for key, occs in line_occurrences.items():
            if len(occs) < 2:
                continue
            occs_sorted = sorted(occs, key=lambda t: t[1])
            keeper_idx = occs_sorted[0][0]
            for ci, _ in occs_sorted[1:]:
                to_remove.setdefault(ci, set()).add(key)

        if to_remove:
            deduped: List[MultipletClusterCandidate] = []
            for ci, cl in enumerate(clusters):
                drop_keys = to_remove.get(ci, set())
                if not drop_keys:
                    deduped.append(cl)
                    continue
                new_comps = tuple(
                    c for c in cl.components
                    if (c[0], round(float(c[1]), 2)) not in drop_keys
                )
                new_phantoms = tuple(
                    c for c in cl.phantom_components
                    if (c[0], round(float(c[1]), 2)) not in drop_keys
                )
                # Drop cluster if <2 active components remain
                if len(new_comps) < 2:
                    continue
                # Recompute display window from remaining active components
                E_lo_sub = new_comps[0][1]
                E_hi_sub = new_comps[-1][1]
                fwhm_lo_sub = float(fwhm_keV_at(float(E_lo_sub)))
                fwhm_hi_sub = float(fwhm_keV_at(float(E_hi_sub)))
                E_lo_disp = float(E_lo_sub) - display_window_fwhm * fwhm_lo_sub
                E_hi_disp = float(E_hi_sub) + display_window_fwhm * fwhm_hi_sub
                deduped.append(MultipletClusterCandidate(
                    cluster_id=cl.cluster_id,
                    E_lo_keV=E_lo_disp,
                    E_hi_keV=E_hi_disp,
                    components=new_comps,
                    n_components=len(new_comps),
                    detection_reason=(
                        cl.detection_reason + " [#37 dedup]"
                    ),
                    found_peaks_in_roi=cl.found_peaks_in_roi,
                    phantom_components=new_phantoms,
                ))
            clusters = deduped

        # Subset-merge: if cluster B's components ⊂ cluster A's components,
        # drop B (it's redundant — A already covers everything B has).
        if len(clusters) >= 2:
            keys_per_cluster = [
                set((c[0], round(float(c[1]), 2)) for c in cl.components)
                for cl in clusters
            ]
            keep_mask = [True] * len(clusters)
            for i in range(len(clusters)):
                if not keep_mask[i]:
                    continue
                for j in range(len(clusters)):
                    if i == j or not keep_mask[j]:
                        continue
                    # If j's components are a strict subset of i's, drop j
                    if (keys_per_cluster[j]
                            and keys_per_cluster[j] < keys_per_cluster[i]):
                        keep_mask[j] = False
            clusters = [c for c, k in zip(clusters, keep_mask) if k]

    return clusters


# ──────────────────────────────────────────────────────────────────
# 3. Decomposition (coupled_fit на каждый кластер)
# ──────────────────────────────────────────────────────────────────

def decompose_multiplets(
    counts: np.ndarray,
    energy_axis_keV: np.ndarray,
    energy_to_channel: Callable[[float], float],
    channel_to_energy: Callable[[int], float],
    fwhm_keV_at: Callable[[float], float],
    clusters: Sequence[MultipletClusterCandidate],
    *,
    use_peak_image: bool = True,
    h_step: float = 0.03,
    continuum: str = "step_linear",
    # BUG-32ζ / task #82 — phantom-inclusion-in-fit.
    # phantom_inclusive=False (default) → текущий v1.18.31+ путь: phantom'ы
    # удаляются из fit, остаются только cl.components. С phantom_inclusive=
    # True phantom'ы передаются в coupled_intensity_fit как явно
    # помеченные independent-компоненты с Tikhonov zero-prior penalty
    # (lambda_phantom_rel). Защищает от «phantom flux absorption» в
    # kept-компоненты, когда top-K cap демотировал реальную линию.
    # lambda_phantom_rel = 1e-3 — empirical default (~0.1% от median(yw)).
    phantom_inclusive: bool = False,
    lambda_phantom_rel: float = 1e-3,
) -> List[Any]:
    """Запускает coupled_intensity_fit на каждый MultipletClusterCandidate.

    Returns
    -------
    list[CoupledFitResult]  (см. gamma.peaks.coupled_multiplet)
    """
    from gamma.peaks.coupled_multiplet import (
        coupled_intensity_fit, ComponentSpec,
    )

    results = []
    for cl in clusters:
        comp_specs = [
            ComponentSpec(nuclide=nuc, E_keV=E, I_gamma_pct=I, group=group)
            for (nuc, E, I, group) in cl.components
        ]
        # BUG-32ζ — build phantom specs только если флаг активен; иначе
        # пустой tuple → coupled_intensity_fit short-circuit'нёт phantom
        # block и сохранит pre-BUG-32ζ control flow.
        if phantom_inclusive and cl.phantom_components:
            phantom_specs = [
                ComponentSpec(nuclide=nuc, E_keV=E, I_gamma_pct=I, group=group)
                for (nuc, E, I, group) in cl.phantom_components
            ]
            lam_rel = float(lambda_phantom_rel)
        else:
            phantom_specs = []
            lam_rel = 0.0
        roi_lo_ch = max(0, int(round(energy_to_channel(cl.E_lo_keV))))
        roi_hi_ch = min(len(counts), int(round(energy_to_channel(cl.E_hi_keV))))
        if roi_hi_ch <= roi_lo_ch + 5:
            continue
        try:
            fit = coupled_intensity_fit(
                energy_keV=energy_axis_keV[roi_lo_ch:roi_hi_ch],
                counts=counts[roi_lo_ch:roi_hi_ch],
                components=comp_specs,
                fwhm_at=fwhm_keV_at,
                continuum=continuum,
                roi_low_ch=roi_lo_ch,
                cluster_id=cl.cluster_id,
                title=cl.detection_reason,
                use_peak_image=use_peak_image,
                h_step=h_step,
                phantom_components=phantom_specs,
                lambda_phantom_rel=lam_rel,
            )
            results.append(fit)
        except Exception as e:
            # graceful — fit может не сойтись на пограничных случаях
            results.append({
                "cluster_id": cl.cluster_id,
                "error": f"{type(e).__name__}: {e}",
                "skipped": True,
            })
    return results


# ──────────────────────────────────────────────────────────────────
# 4. Orchestrator
# ──────────────────────────────────────────────────────────────────

def run_v2_pipeline(
    spec,
    *,
    sigma_threshold: float = 3.0,
    chain_filter: Optional[Sequence[str]] = None,
    fwhm_fallback_R662_pct: float = 7.0,
    # BUG-32ζ / task #82 — forwarded to decompose_multiplets.
    phantom_inclusive: bool = False,
    lambda_phantom_rel: float = 1e-3,
) -> PipelineV2Result:
    """Один-call orchestrator.

    Steps:
      1. dual-method search
      2. cluster detection
      3. decompose

    Parameters
    ----------
    spec : Spectrum
        Любой объект с .counts, .channel_to_energy, .energy_to_channel.
    sigma_threshold : float
        Порог S/σ для обоих search-методов.
    chain_filter : sequence, optional
        Ограничить библиотеку до этих нуклидов (например ("Ac-228", "Tl-208", ...))
        Если None — вся библиотека.
    fwhm_fallback_R662_pct : float
        Если stored FWHM-cal даёт ерунду (R(122)>15%) — используется
        эмпирическая модель R(662)=N% для NaI.
    """
    from gamma.calibration.fwhm_provider import make_fwhm_at_channel_provider
    import math

    counts = np.asarray(spec.counts, dtype=float)
    ch_to_E = spec.channel_to_energy
    E_to_ch = spec.energy_to_channel
    fwhm_ch_provider = make_fwhm_at_channel_provider(spec)

    energy_axis = np.array([ch_to_E(i) for i in range(len(counts))])

    # FWHM(E) wrapper с fallback на R(662)=N% модель если stored cal даёт ерунду
    def fwhm_keV_at(E: float) -> float:
        ch = int(round(E_to_ch(E)))
        if ch < 0 or ch >= len(counts):
            return 0.07 * (E * 662.0) ** 0.5
        fwhm_ch = fwhm_ch_provider(ch)
        keV_per_ch = ch_to_E(min(ch + 1, len(counts) - 1)) - ch_to_E(ch)
        fwhm = fwhm_ch * keV_per_ch
        # Fallback при выходе за разумные пределы
        if fwhm < 5.0 or fwhm > E * 0.3:
            return (fwhm_fallback_R662_pct / 100.0) * math.sqrt(E * 662.0)
        return fwhm

    notes: List[str] = []

    # 1. Search
    merged, mari_raw, conv_raw = search_dual_method(
        counts, fwhm_ch_provider,
        sigma_threshold=sigma_threshold,
        channel_to_energy=ch_to_E,
    )
    notes.append(
        f"search: Mariscotti={len(mari_raw)}, matched_filter={len(conv_raw)}, "
        f"merged={len(merged)} (overlap "
        f"{sum(1 for h in merged if h.source == 'both')}, "
        f"mari_only={sum(1 for h in merged if h.source == 'mariscotti')}, "
        f"conv_only={sum(1 for h in merged if h.source == 'matched_filter')})"
    )

    # 2. Detect clusters
    clusters = detect_multiplet_clusters(
        merged, fwhm_keV_at,
        chain_filter=chain_filter,
    )
    notes.append(f"clusters detected: {len(clusters)}")

    # 3. Decompose
    fits = decompose_multiplets(
        counts, energy_axis, E_to_ch, ch_to_E, fwhm_keV_at, clusters,
        phantom_inclusive=phantom_inclusive,
        lambda_phantom_rel=lambda_phantom_rel,
    )
    n_converged = sum(
        1 for f in fits
        if not isinstance(f, dict) and getattr(f, "converged", False)
    )
    notes.append(f"coupled_fits: {len(fits)} attempted, {n_converged} converged")

    return PipelineV2Result(
        found_peaks=merged,
        mariscotti_peaks=list(mari_raw),
        matched_filter_peaks=list(conv_raw),
        multiplet_candidates=list(clusters),
        coupled_fits=list(fits),
        notes=notes,
    )


# ──────────────────────────────────────────────────────────────────
# 5. Comparison with production
# ──────────────────────────────────────────────────────────────────

def compare_with_production(
    sample_path: str,
    *,
    background_path: Optional[str] = None,
    sample_mass_kg: Optional[float] = None,
    sigma_threshold: float = 3.0,
    chain_filter: Optional[Sequence[str]] = None,
) -> Tuple[Any, PipelineV2Result, ComparisonReport]:
    """Запускает PRODUCTION (analyze_lsrm_spe) И v2-pipeline на одном spec.

    Returns
    -------
    (production_result, v2_result, comparison_report)
    """
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    from gamma.io.readers import read_spectrum

    # Production
    production = analyze_lsrm_spe(
        sample_path,
        background_path=background_path,
        sample_mass_kg=sample_mass_kg,
        complete_workflow=True,
    )

    # v2 — на том же спектре (через тот же read_spectrum)
    spec = read_spectrum(sample_path)
    v2 = run_v2_pipeline(
        spec,
        sigma_threshold=sigma_threshold,
        chain_filter=chain_filter,
    )

    # ──── Сравнение ────
    rpt = ComparisonReport()

    # Production peaks — список FoundPeak в `result.peaks` (после search,
    # ДО идентификации). Энергия выводится через калибровку.
    prod_peaks = []
    try:
        spec = production.spec
        ch_to_E = spec.channel_to_energy
        for fp in (getattr(production, "peaks", None) or []):
            ch = int(getattr(fp, "channel", 0))
            E = float(ch_to_E(ch))
            prod_peaks.append((ch, E))
    except Exception:
        pass
    rpt.production_n_peaks = len(prod_peaks)

    # Production multiplets — `deconvolution_results` (list of CoupledFitResult)
    # либо `multiplet_results` (legacy)
    prod_mult_ids = []
    try:
        for m in (getattr(production, "deconvolution_results", None) or []):
            prod_mult_ids.append(getattr(m, "id", "?"))
        if not prod_mult_ids:
            for m in (getattr(production, "multiplet_results", None) or []):
                prod_mult_ids.append(getattr(m, "id", "?"))
    except Exception:
        pass
    rpt.production_n_multiplets = len(prod_mult_ids)
    rpt.production_multiplet_ids = list(prod_mult_ids)

    # v2
    rpt.v2_n_peaks = len(v2.found_peaks)
    v2_peaks = [(h.channel, h.energy_keV) for h in v2.found_peaks]
    rpt.v2_n_multiplets = len(v2.multiplet_candidates)
    rpt.v2_multiplet_ids = [c.cluster_id for c in v2.multiplet_candidates]

    # Diff peaks (по каналу с допуском ±3)
    def _match(a, b, tol=3):
        return abs(a - b) <= tol
    v2_set = {ch for ch, _ in v2_peaks}
    prod_set = {ch for ch, _ in prod_peaks}
    rpt.peaks_only_in_v2 = sorted([
        (ch, E) for ch, E in v2_peaks
        if not any(_match(ch, p) for p in prod_set)
    ])
    rpt.peaks_only_in_production = sorted([
        (ch, E) for ch, E in prod_peaks
        if not any(_match(ch, v) for v in v2_set)
    ])
    rpt.multiplets_only_in_v2 = [
        c.cluster_id for c in v2.multiplet_candidates
        if c.cluster_id not in prod_mult_ids
    ]
    rpt.multiplets_only_in_production = [
        m for m in prod_mult_ids if m not in rpt.v2_multiplet_ids
    ]

    rpt.notes.extend(v2.notes)
    rpt.notes.append(
        f"production primary_feps={len(prod_peaks)}, multiplets={len(prod_mult_ids)}"
    )
    return production, v2, rpt


__all__ = [
    "PeakHit",
    "MultipletClusterCandidate",
    "PipelineV2Result",
    "ComparisonReport",
    "search_dual_method",
    "detect_multiplet_clusters",
    "decompose_multiplets",
    "run_v2_pipeline",
    "compare_with_production",
]
