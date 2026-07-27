# -*- coding: utf-8 -*-
"""F-167 — Полное покрытие тестами для id_window.py + build_id_window_k_fwhm.

Покрывает F-167-1 (NaI K-40 включение/исключение), F-167-2 (HPGe K-40),
F-167-3 (таблица k для NaI/CsI/LaBr/CeBr/HPGe/CdZnTe), F-167-4 (regression:
сравнение легаси `δE₀·√(E/Eref)` vs canonical `k·FWHM(E)`).

F-167-5 / F-167-6 — в `test_kfwhm_contexts.py` (подзадача 8).
"""
from __future__ import annotations

import math

import pytest

from gamma.identification.id_window import (
    ID_WINDOW_K_FWHM,
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
    id_window_keV,
    normalize_detector_class,
)


def test_id_window_k_fwhm_table_canonical_values():
    """LSRM-Algo-9: NaI/CsI k=1.5, HPGe/LaBr/CeBr/CdZnTe k=1.0."""
    assert ID_WINDOW_K_FWHM["NaI"] == 1.5
    assert ID_WINDOW_K_FWHM["CsI"] == 1.5
    assert ID_WINDOW_K_FWHM["LaBr"] == 1.0
    assert ID_WINDOW_K_FWHM["CeBr"] == 1.0
    assert ID_WINDOW_K_FWHM["HPGe"] == 1.0
    assert ID_WINDOW_K_FWHM["CdZnTe"] == 1.0


def test_ortec_kfwhm_constants_unique_and_set():
    """F-167-5 stub: 8 ORTEC-контекстов имеют уникальные значения."""
    contexts = {
        "background_auto":     BACKGROUND_SEARCH_K_MW_FWHM,   # 6.0
        "background_xfwhm":    BACKGROUND_X_FWHM,             # 1.0
        "pbc":                 PBC_MATCH_K_FWHM,              # 0.5
        "mult_criterion":      MULTIPLET_CRITERION_K_FWHM,    # 3.08
        "mult_padding":        MULTIPLET_PADDING_K_FWHM,      # 1.5
        "lib_unknown":         LIB_UNKNOWN_ASSOC_K_FWHM,      # 3.3
        "iso_peak":            ISO_NORM_PEAK_DOMINANT_K_FWHM, # 2.5
        "iso_bg":              ISO_NORM_BG_DOMINANT_K_FWHM,   # 1.2
        "directed_fit":        DIRECTED_FIT_K_FWHM_BASE,      # 4.84
        "phase_d_reg":         PHASE_D_REGULARIZATION_K_FWHM, # 0.15
    }
    # 10 различных контекстов — 10 различных значений (без коллизий)
    assert len(set(contexts.values())) == len(contexts), (
        f"Коллизия значений k·FWHM: {contexts}"
    )


def test_id_window_keV_nai_K40():
    """F-167-1: K-40 1461 кэВ на NaI, FWHM=120 кэВ → окно 360 кэВ полное."""
    w = id_window_keV(1461.0, 120.0, "NaI")
    assert math.isclose(w, 360.0, rel_tol=1e-9)
    # ±-сравнение: half-width = 180 кэВ
    assert (w / 2.0) == 180.0


def test_id_window_keV_hpge_K40():
    """F-167-2: K-40 1461 кэВ на HPGe, FWHM=1.8 кэВ → окно 3.6 кэВ полное."""
    w = id_window_keV(1461.0, 1.8, "HPGe")
    assert math.isclose(w, 3.6, rel_tol=1e-9)


def test_normalize_detector_class():
    """normalize_detector_class покрывает legacy-форматы."""
    assert normalize_detector_class("NaI 63×63") == "NaI"
    assert normalize_detector_class("HPGe Canberra GR2018") == "HPGe"
    assert normalize_detector_class("LaBr3") == "LaBr"
    assert normalize_detector_class("CeBr3") == "CeBr"
    assert normalize_detector_class("CsI(Tl)") == "CsI"
    assert normalize_detector_class("CdZnTe") == "CdZnTe"
    assert normalize_detector_class(None) == "NaI"
    assert normalize_detector_class("") == "NaI"


def test_normalize_detector_class_invalid_falls_back_to_nai():
    """Unknown ввод → NaI (профиль Gamma-1S per F-157)."""
    assert normalize_detector_class("garbage_detector") == "NaI"


@pytest.mark.parametrize("detector,k", [
    ("NaI", 1.5), ("CsI", 1.5),
    ("LaBr", 1.0), ("CeBr", 1.0), ("HPGe", 1.0), ("CdZnTe", 1.0),
])
def test_id_window_keV_scales_linearly_with_fwhm(detector, k):
    """Для любого detector: window = 2·k·FWHM (формула канонична)."""
    fwhm = 5.0
    expected = 2.0 * k * fwhm
    assert math.isclose(id_window_keV(1000.0, fwhm, detector), expected, rel_tol=1e-9)


# ════════════════════════════════════════════════════════════════════
# F-167-1 — NaI K-40 1461 кэВ: матч/не-матч в окне ±k·FWHM
# ════════════════════════════════════════════════════════════════════

def test_F167_1_nai_K40_match_decisions():
    """F-167-1: ±k·FWHM = ±180 кэВ; (1369, 1461) и (1553, 1461) внутри,
    (1280, 1461) и (1642, 1461) снаружи."""
    fwhm_at_1461 = 120.0
    half_window = id_window_keV(1461.0, fwhm_at_1461, "NaI") / 2.0
    assert half_window == 180.0

    # Includes (within ±180 keV)
    assert abs(1369.0 - 1461.0) <= half_window  # Δ=92, in
    assert abs(1553.0 - 1461.0) <= half_window  # Δ=92, in
    # Excludes (outside ±180 keV)
    assert abs(1280.0 - 1461.0) > half_window   # Δ=181, NOT in
    assert abs(1642.0 - 1461.0) > half_window   # Δ=181, NOT in


# ════════════════════════════════════════════════════════════════════
# F-167-4 — Regression: legacy sqrt_E vs canonical k·FWHM
# ════════════════════════════════════════════════════════════════════

def _toy_fwhm_provider(E_keV: float) -> float:
    """Toy FWHM model: ~60 кэВ @ 661, scaling as √E (typical NaI 63×63)."""
    return 60.0 * math.sqrt(E_keV / 661.66)


@pytest.mark.parametrize("E,fwhm_E,expect_canonical_wider", [
    (661.66, 60.0, True),    # at ref: canonical=90, legacy=30 → 3×
    (1461.0, 89.2, True),    # at K-40: canonical=134, legacy=44 → 3×
    (2614.51, 119.4, True),  # at Tl-208: canonical=179, legacy=60 → 3×
])
def test_F167_4_canonical_wider_than_legacy(E, fwhm_E, expect_canonical_wider):
    """F-167-4: canonical k·FWHM окно должно быть **шире** legacy на всех
    типичных энергиях NaI (660-2614 keV). Это и есть фикс false-negatives."""
    canonical_half = id_window_keV(E, fwhm_E, "NaI") / 2.0
    # Legacy сравнение: δE₀ = 0.5·FWHM_661 = 30 keV
    legacy_half = 30.0 * math.sqrt(E / 661.66)
    if expect_canonical_wider:
        assert canonical_half > legacy_half, (
            f"F-167 canonical ({canonical_half:.1f}) must be wider than "
            f"legacy ({legacy_half:.1f}) at E={E}"
        )
        # Ratio must be ~3× (k_canonical=1.5, k_legacy=0.5 → 3:1)
        ratio = canonical_half / legacy_half
        assert 2.8 < ratio < 3.2, f"Expected ~3× wider, got {ratio:.2f}×"


# ════════════════════════════════════════════════════════════════════
# build_id_window_k_fwhm factory tests
# ════════════════════════════════════════════════════════════════════

def test_build_id_window_k_fwhm_nai():
    """Фабрика возвращает scaling='k_fwhm' + k=1.5 для NaI."""
    from gamma.identification.window import build_id_window_k_fwhm
    w = build_id_window_k_fwhm("NaI", fwhm_provider_keV=_toy_fwhm_provider)
    assert w.scaling == "k_fwhm"
    assert w.k_fwhm == 1.5
    assert w.detector_type == "NaI"
    # window_keV(661) = k * FWHM(661) = 1.5 * 60 = 90 keV (half-width)
    assert math.isclose(w.window_keV(661.66), 90.0, rel_tol=1e-6)


def test_build_id_window_k_fwhm_hpge():
    """Фабрика для HPGe: k=1.0."""
    from gamma.identification.window import build_id_window_k_fwhm
    hpge_fwhm = lambda E: 1.5 + 0.0005 * E  # typical HPGe: ~1.8 keV @ 661
    w = build_id_window_k_fwhm("HPGe", fwhm_provider_keV=hpge_fwhm)
    assert w.k_fwhm == 1.0
    # window_keV(1461) = 1.0 * (1.5 + 0.0005·1461) = 1.0 * 2.23 = 2.23 keV
    assert math.isclose(w.window_keV(1461.0), 2.2305, rel_tol=1e-4)


def test_build_id_window_k_fwhm_with_k_override():
    """k_override переопределяет per-detector default."""
    from gamma.identification.window import build_id_window_k_fwhm
    w = build_id_window_k_fwhm(
        "NaI", fwhm_provider_keV=_toy_fwhm_provider, k_override=2.0,
    )
    assert w.k_fwhm == 2.0  # not the NaI default 1.5
    assert math.isclose(w.window_keV(661.66), 120.0, rel_tol=1e-6)  # 2.0·60


def test_build_id_window_k_fwhm_normalizes_legacy_strings():
    """Фабрика принимает legacy detector_type strings и нормализует."""
    from gamma.identification.window import build_id_window_k_fwhm
    w1 = build_id_window_k_fwhm("NaI 63×63", fwhm_provider_keV=_toy_fwhm_provider)
    w2 = build_id_window_k_fwhm("HPGe Canberra GR2018", fwhm_provider_keV=lambda E: 2.0)
    assert w1.detector_type == "NaI" and w1.k_fwhm == 1.5
    assert w2.detector_type == "HPGe" and w2.k_fwhm == 1.0


def test_build_id_window_k_fwhm_safe_fallback_on_bad_fwhm():
    """При FWHM≤0 или NaN — fallback на δE₀ (защитный, не падает)."""
    from gamma.identification.window import build_id_window_k_fwhm

    def bad_provider(E):
        if E < 100.0:
            return float("nan")
        if E > 3000.0:
            return -1.0  # отрицательный
        return 60.0 * math.sqrt(E / 661.66)

    w = build_id_window_k_fwhm("NaI", fwhm_provider_keV=bad_provider)
    # Valid range: returns k·FWHM
    assert math.isclose(w.window_keV(661.66), 90.0, rel_tol=1e-6)
    # NaN at low E: fallback to delta_E0_keV
    result_low = w.window_keV(50.0)
    assert math.isfinite(result_low) and result_low > 0
    # Negative at high E: fallback to delta_E0_keV
    result_hi = w.window_keV(3500.0)
    assert math.isfinite(result_hi) and result_hi > 0


def test_k_fwhm_scaling_raises_without_provider():
    """ValueError если scaling='k_fwhm' но fwhm_provider_keV не задан."""
    from gamma.identification.window import IdentificationWindow
    w = IdentificationWindow(
        detector_type="NaI", delta_E0_keV=15.0, scaling="k_fwhm",
        k_fwhm=1.5, fwhm_provider_keV=None,
    )
    with pytest.raises(ValueError, match="F-167"):
        w.window_keV(661.66)


# ════════════════════════════════════════════════════════════════════
# F-167-3 — расширенная таблица детекторов
# ════════════════════════════════════════════════════════════════════

def test_F167_3_per_detector_window_widths_at_K40():
    """F-167-3: на K-40 1461 кэВ (FWHM=90 keV NaI / 2.5 HPGe) окна различаются
    по типу детектора в правильную сторону."""
    nai_window = id_window_keV(1461.0, 90.0, "NaI")  # k=1.5
    csi_window = id_window_keV(1461.0, 90.0, "CsI")  # k=1.5
    hpge_window = id_window_keV(1461.0, 2.5, "HPGe")  # k=1.0
    labr_window = id_window_keV(1461.0, 30.0, "LaBr")  # k=1.0
    assert nai_window == csi_window  # NaI/CsI share k=1.5
    assert hpge_window < labr_window < nai_window  # rank by FWHM·k
    assert nai_window == 270.0
    assert hpge_window == 5.0
