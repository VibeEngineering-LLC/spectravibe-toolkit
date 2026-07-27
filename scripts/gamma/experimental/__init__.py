"""F-354 / v1.18.24.0 — экспериментальный pipeline-стек.

Запускается **параллельно** с production analyze_lsrm_spe (не заменяет!)
для сравнения новых алгоритмов на реальных данных.

Сейчас содержит:
  * peak_pipeline_v2 — dual-method search (Mariscotti ∪ matched filter) +
    автоматическое обнаружение мультиплетов БЕЗ FORCED_CLUSTERS таблиц +
    coupled_intensity_fit на каждый детектированный кластер.

Контракт: эти модули НЕ должны импортироваться из production
identification.staged_pipeline. Они вызываются только из:
  - CLI флаг `--experimental-compare` (TBD)
  - tests/snapshot/test_f354_*.py
  - prompt-driven scripts (как в текущем walkthrough)
"""
from gamma.experimental.peak_pipeline_v2 import (
    PipelineV2Result,
    ComparisonReport,
    search_dual_method,
    detect_multiplet_clusters,
    decompose_multiplets,
    run_v2_pipeline,
    compare_with_production,
)

__all__ = [
    "PipelineV2Result",
    "ComparisonReport",
    "search_dual_method",
    "detect_multiplet_clusters",
    "decompose_multiplets",
    "run_v2_pipeline",
    "compare_with_production",
]
