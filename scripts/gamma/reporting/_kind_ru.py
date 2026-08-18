# -*- coding: utf-8 -*-
"""TD-4 / v1.18.30 — Shared feature_kind Russian translations.

Единый источник RU-меток для типов вторичных пиков (`feature_kind` /
`residual_classification.label`). Устраняет дублирование между
`markdown_report._CELL_TRANSLATIONS` и `plots._PLOT_KIND_RU`.

Две таблицы:
* ``FEATURE_KIND_RU`` — полные метки (для таблиц Markdown, Technical PDF,
  interactive HTML).
* ``FEATURE_KIND_RU_SHORT`` — сокращённые метки (для аннотаций графиков
  matplotlib, где место ограничено).

Оба словаря используют одни и те же ключи. Потребители:
* ``gamma.reporting.markdown_report`` — ``_CELL_TRANSLATIONS`` (merge)
* ``gamma.reporting.plots`` — ``_PLOT_KIND_RU``
* ``gamma.reporting.interactive_html`` (через ``_CELL_TRANSLATIONS``)
* при необходимости — ``gamma.reporting.technical_pdf``

Расширение: добавьте новый ключ в оба словаря одновременно.
Если сокращение совпадает с полной формой — просто скопируйте.
"""
from __future__ import annotations

from typing import Dict

# ──────────────────────────────────────────────────────────────────────────────
# Полные метки (для текстовых отчётов)
# ──────────────────────────────────────────────────────────────────────────────
FEATURE_KIND_RU: Dict[str, str] = {
    "single_escape": "пик одиночного вылета",
    "double_escape": "пик двойного вылета",
    "annihilation_511": "аннигиляция 511",
    "sum_peak": "сумм. пик",
    "compton_edge": "Комптон-край",
    "compton_plateau": "Комптон-плато",
    "backscatter": "обратное рассеяние",
    "backscatter_region": "обл. обратного рассеяния",
    "fluorescence_shield": "флуоресценция защиты",
    "fluorescence_matrix": "флуоресценция матрицы",
    "fluorescence_shield_collimator": "флуоресценция коллиматора",
    "broad_compton_plateau": "широкое Комптон-плато",
    "chain_secondary": "вторичная по цепочке",
    "composite_cluster": "композит-кластер",
    # Escape-peak семейство (F-386 hard-lock: «вылет», не «ускользание»)
    "I_K_escape_Ka": "I K-вылет (Kα)",
    "I_K_escape_Kb": "I K-вылет (Kβ)",
}

# ──────────────────────────────────────────────────────────────────────────────
# Сокращённые метки (для аннотаций matplotlib и компактных UI-элементов)
# ──────────────────────────────────────────────────────────────────────────────
FEATURE_KIND_RU_SHORT: Dict[str, str] = {
    "single_escape": "пик одиноч. вылета",
    "double_escape": "пик двойн. вылета",
    "annihilation_511": "511",
    "sum_peak": "сумм. пик",
    "compton_edge": "Комптон-край",
    "compton_plateau": "Комптон-плато",
    "backscatter": "обр. рассеяние",
    "backscatter_region": "обл. обр. рассеяния",
    "fluorescence_shield": "флуор. защиты",
    "fluorescence_matrix": "флуор. матрицы",
    "fluorescence_shield_collimator": "флуор. коллим.",
    "broad_compton_plateau": "шир. Комптон-плато",
    "chain_secondary": "по цепочке",
    "composite_cluster": "композ. кластер",
    "I_K_escape_Ka": "I K-вылет Kα",
    "I_K_escape_Kb": "I K-вылет Kβ",
}


def feature_kind_ru(kind: str, *, short: bool = False) -> str:
    """Translate a ``feature_kind`` / ``label`` string to Russian.

    Args:
        kind: Raw identifier (e.g. ``"single_escape"``).
        short: If True, return abbreviated form (suitable for plot labels).
               If False (default), return full form.

    Returns:
        Russian translation, or *kind* unchanged if not in the mapping.
    """
    table = FEATURE_KIND_RU_SHORT if short else FEATURE_KIND_RU
    return table.get(kind, kind)
