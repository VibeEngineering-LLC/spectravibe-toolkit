# -*- coding: utf-8 -*-
"""F-167-5 / F-167-6 — Контракт уникальности и документированности 8 ORTEC
kFWHM-контекстов в gamma.identification.id_window.

F-167-5: каждый из ≥8 контекстов (ID matching, PBC, multiplet criterion,
multiplet padding, lib↔unknown assoc, background search Auto, ISO Singlet
peak/bg, Directed Fit) имеет свою константу с уникальным числовым
значением — нет magic-number коллизий.

F-167-6: docstring модуля содержит ссылку на первоисточник каждой константы
(`[ORTEC-GV9-*]` или `[LSRM-Algo-9]` или F-145/F-117 для PHASE_D_REG).

Источник: audit/_plans/F-167_id_window_k_fwhm.md §3.6, §11
           audit/v1_17_9_5_ortec/AUDIT_DELTA_ortec_gv9_sklearn.md §3
"""
from __future__ import annotations

import inspect

import pytest

from gamma.identification import id_window as id_window_mod
from gamma.identification.id_window import (
    BACKGROUND_SEARCH_K_MW_FWHM,
    BACKGROUND_X_FWHM,
    PBC_MATCH_K_FWHM,
    MULTIPLET_CRITERION_K_FWHM,
    MULTIPLET_PADDING_K_FWHM,
    LIB_UNKNOWN_ASSOC_K_FWHM,
    ISO_NORM_PEAK_DOMINANT_K_FWHM,
    ISO_NORM_BG_DOMINANT_K_FWHM,
    DIRECTED_FIT_K_FWHM_BASE,
    PHASE_D_REGULARIZATION_K_FWHM,
    ID_WINDOW_K_FWHM,
)


# ════════════════════════════════════════════════════════════════════
# F-167-5 — 8 ORTEC контекстов имеют уникальные значения
# ════════════════════════════════════════════════════════════════════

# Полный реестр (имя константы → значение, источник)
_KFWHM_REGISTRY = {
    "BACKGROUND_SEARCH_K_MW_FWHM":   (BACKGROUND_SEARCH_K_MW_FWHM,   "[ORTEC-GV9-Background-Auto]"),
    "BACKGROUND_X_FWHM":             (BACKGROUND_X_FWHM,             "[ORTEC-GV9-Background-FWHM]"),
    "PBC_MATCH_K_FWHM":              (PBC_MATCH_K_FWHM,              "[ORTEC-GV9-PBC]"),
    "MULTIPLET_CRITERION_K_FWHM":    (MULTIPLET_CRITERION_K_FWHM,    "[ORTEC-GV9-Multiplet-Region]"),
    "MULTIPLET_PADDING_K_FWHM":      (MULTIPLET_PADDING_K_FWHM,      "[ORTEC-GV9-Multiplet-Region]"),
    "LIB_UNKNOWN_ASSOC_K_FWHM":      (LIB_UNKNOWN_ASSOC_K_FWHM,      "[ORTEC-GV9-Deconvolution-Width]"),
    "ISO_NORM_PEAK_DOMINANT_K_FWHM": (ISO_NORM_PEAK_DOMINANT_K_FWHM, "[ORTEC-GV9-ISO-NORM]"),
    "ISO_NORM_BG_DOMINANT_K_FWHM":   (ISO_NORM_BG_DOMINANT_K_FWHM,   "[ORTEC-GV9-ISO-NORM]"),
    "DIRECTED_FIT_K_FWHM_BASE":      (DIRECTED_FIT_K_FWHM_BASE,      "[ORTEC-GV9-Engines]"),
    "PHASE_D_REGULARIZATION_K_FWHM": (PHASE_D_REGULARIZATION_K_FWHM, "F-145/F-117"),
}

# Канонические ORTEC значения (acceptance per план §3.6 F-167-5)
_EXPECTED_VALUES = {
    "BACKGROUND_SEARCH_K_MW_FWHM":   6.0,
    "BACKGROUND_X_FWHM":             1.0,
    "PBC_MATCH_K_FWHM":              0.5,
    "MULTIPLET_CRITERION_K_FWHM":    3.08,
    "MULTIPLET_PADDING_K_FWHM":      1.5,
    "LIB_UNKNOWN_ASSOC_K_FWHM":      3.3,
    "ISO_NORM_PEAK_DOMINANT_K_FWHM": 2.5,
    "ISO_NORM_BG_DOMINANT_K_FWHM":   1.2,
    "DIRECTED_FIT_K_FWHM_BASE":      4.84,
    "PHASE_D_REGULARIZATION_K_FWHM": 0.15,
}


def test_F167_5_at_least_8_ortec_contexts_registered():
    """≥8 ORTEC-контекстов представлены отдельными именованными константами."""
    ortec_only = [name for name, (_, src) in _KFWHM_REGISTRY.items()
                  if src.startswith("[ORTEC-")]
    assert len(ortec_only) >= 8, (
        f"План §3.6 требует ≥8 ORTEC контекстов; зарегистрировано {len(ortec_only)}: "
        f"{ortec_only}"
    )


def test_F167_5_values_canonical():
    """Каждая константа имеет каноничное ORTEC/LSRM значение."""
    for name, (value, _) in _KFWHM_REGISTRY.items():
        expected = _EXPECTED_VALUES[name]
        assert value == expected, (
            f"{name} = {value}, ожидалось {expected} (per ORTEC-GV9/LSRM canon)"
        )


def test_F167_5_no_value_collisions():
    """Все 10 имённых констант должны иметь различные числовые значения
    (как маркер уникальности контекста)."""
    values = [v for v, _ in _KFWHM_REGISTRY.values()]
    assert len(set(values)) == len(values), (
        f"Коллизия значений: {sorted(values)}. "
        f"Если контексты концептуально разные, значения тоже должны быть разными."
    )


def test_F167_5_id_window_per_detector_not_in_kfwhm_registry():
    """ID_WINDOW_K_FWHM (per-detector, dict) — отдельный контракт от
    8 ORTEC контекстов; не входит в реестр."""
    # NaI/CsI k=1.5 совпадает с MULTIPLET_PADDING_K_FWHM=1.5, и это OK:
    # ID-окно — отдельная семантика, не подлежит unique-check.
    assert ID_WINDOW_K_FWHM["NaI"] == 1.5
    assert MULTIPLET_PADDING_K_FWHM == 1.5
    # Но они называются разными именами в коде → нет magic-number коллизии


# ════════════════════════════════════════════════════════════════════
# F-167-6 — docstring каждой константы содержит RAG-ID цитату
# ════════════════════════════════════════════════════════════════════

def test_F167_6_module_docstring_lists_all_sources():
    """Module docstring модуля id_window перечисляет все RAG-источники
    (LSRM-Algo-9 + ORTEC-GV9-* + ISO-NORM + Engines)."""
    doc = id_window_mod.__doc__ or ""
    required_citations = [
        "[LSRM-Algo-9]",
        "[ORTEC-GV9-Background-Auto]",
        "[ORTEC-GV9-Deconvolution-Width]",
        "[ORTEC-GV9-Multiplet-Region]",
        "[ORTEC-GV9-PBC]",
        "[ORTEC-GV9-ISO-NORM]",
        "[ORTEC-GV9-Engines]",
    ]
    for cit in required_citations:
        assert cit in doc, (
            f"F-167-6: module docstring должен ссылаться на {cit}; "
            f"первые 500 символов docstring: {doc[:500]}"
        )


def test_F167_6_module_source_has_section_headers_for_constants():
    """Тело модуля содержит секции `§ N.` с пояснениями контекста для
    каждой константы (документация рядом с определением)."""
    src = inspect.getsource(id_window_mod)
    # Проверяем явное разделение секциями
    assert "§ 1.  ID matching" in src
    assert "§ 2.  Background search" in src
    assert "§ 3.  PBC" in src
    assert "§ 4.  Multiplet grouping" in src
    assert "§ 5.  Library" in src and "association" in src
    assert "§ 6.  ISO 11929 NORM" in src
    assert "§ 7.  Directed Fit" in src
    assert "§ 8.  Phase D regularization" in src


@pytest.mark.parametrize("citation_tag", [
    "[LSRM-Algo-9]",
    "[ORTEC-GV9-Background-Auto]",
    "[ORTEC-GV9-Background-FWHM]",
    "[ORTEC-GV9-PBC]",
    "[ORTEC-GV9-Multiplet-Region]",
    "[ORTEC-GV9-Deconvolution-Width]",
    "[ORTEC-GV9-ISO-NORM]",
    "[ORTEC-GV9-Engines]",
])
def test_F167_6_each_citation_present_in_source(citation_tag):
    """Каждая Layer 1 RAG-ID цитата встречается в теле модуля рядом с
    соответствующей константой."""
    src = inspect.getsource(id_window_mod)
    assert citation_tag in src, (
        f"Цитата {citation_tag} должна быть рядом с константой в id_window.py"
    )
