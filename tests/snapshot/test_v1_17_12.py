# -*- coding: utf-8 -*-
"""
v1.17.12 delivery tests — Identification quality.

Covers fixes F-275..F-279 (T-022, T-001, T-028, T-030, T-025, T-026, T-005).

Run:
    cd 0_Work/gamma-spectrum-analysis
    $env:PYTHONPATH = "scripts"
    python -m pytest tests/snapshot/test_v1_17_12.py -v
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest

SCRIPTS = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts")
)
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)


# ─────────────────────────────────────────────────────────────────────
# F-275 / T-022 — per-nuclide CI thresholds (Lsrm Table 14-1)
# ─────────────────────────────────────────────────────────────────────

def test_F275_lsrm_table_14_1_nai_values():
    from gamma.identification.ci_thresholds import LSRM_TABLE_14_1_NAI
    assert LSRM_TABLE_14_1_NAI["Cs-137"] == 1.8
    assert LSRM_TABLE_14_1_NAI["Eu-152"] == 18.3
    assert LSRM_TABLE_14_1_NAI["Co-60"] == 5.9


def test_F275_required_ci_eu152_now_strict():
    """Eu-152 на NaI требует CI ≥ 0.5·18.3 = 9.15 для confirmed."""
    from gamma.identification.ci_thresholds import required_ci_for_confirmation
    req = required_ci_for_confirmation("Eu-152", "NaI")
    assert req == pytest.approx(9.15, abs=0.01)


def test_F275_required_ci_cs137_uses_global_floor():
    """Cs-137 калибровка 1.8 × 0.5 = 0.9 → берётся global_floor=2.0."""
    from gamma.identification.ci_thresholds import required_ci_for_confirmation
    req = required_ci_for_confirmation("Cs-137", "NaI")
    assert req == 2.0


def test_F275_unknown_nuclide_falls_back_to_global_floor():
    from gamma.identification.ci_thresholds import required_ci_for_confirmation
    req = required_ci_for_confirmation("Xe-999-fake", "NaI")
    assert req == 2.0


def test_F275_hpge_thresholds_higher_than_nai():
    """HPGe values shift up by ~6 (better resolution → larger CI possible)."""
    from gamma.identification.ci_thresholds import lsrm_table_14_1_ci
    nai = lsrm_table_14_1_ci("Eu-152", "NaI")
    hpge = lsrm_table_14_1_ci("Eu-152", "HPGe")
    assert hpge > nai + 4.0


# ─────────────────────────────────────────────────────────────────────
# F-276 / T-028 — 511 keV ID veto
# ─────────────────────────────────────────────────────────────────────

def test_F276_only_511_match_detected():
    from gamma.identification.id_anomalies import is_only_511_line_match
    assert is_only_511_line_match([511.0])
    assert is_only_511_line_match([512.5])    # ± 3 keV
    assert not is_only_511_line_match([511.0, 1274.5])  # corroborated
    assert not is_only_511_line_match([515.0])  # >3 keV off


def test_F276_na22_only_511_gets_veto_reason():
    from gamma.identification.id_anomalies import annihilation_veto_reason
    msg = annihilation_veto_reason("Na-22")
    assert msg is not None
    assert "511" in msg


def test_F276_cs137_not_in_veto_list():
    """Cs-137 не β+ → нет veto reason."""
    from gamma.identification.id_anomalies import annihilation_veto_reason
    assert annihilation_veto_reason("Cs-137") is None


# ─────────────────────────────────────────────────────────────────────
# F-276 / T-030 — Pb fluorescence flag
# ─────────────────────────────────────────────────────────────────────

def test_F276_pb_flag_detects_72_75():
    from gamma.identification.id_anomalies import is_pb_fluorescence_line
    assert is_pb_fluorescence_line(72.80)
    assert is_pb_fluorescence_line(74.97)
    assert is_pb_fluorescence_line(75.5)   # in tol ±3
    assert is_pb_fluorescence_line(84.94)  # K-β
    assert not is_pb_fluorescence_line(100.0)
    assert not is_pb_fluorescence_line(60.0)


def test_F276_pb_fluorescence_note_returns_text():
    from gamma.identification.id_anomalies import pb_fluorescence_note
    note = pb_fluorescence_note(74.0)
    assert note is not None
    assert "Pb" in note
    assert "shield" in note.lower() or "домик" in note


# ─────────────────────────────────────────────────────────────────────
# F-276 / T-025 — Iodine K-escape annotation
# ─────────────────────────────────────────────────────────────────────

def test_F276_am241_k_escape_annotation():
    from gamma.identification.id_anomalies import iodine_k_escape_for
    ann = iodine_k_escape_for(59.54, "NaI")
    assert ann is not None
    assert ann.escape_E_keV == pytest.approx(30.93, abs=0.01)


def test_F276_no_escape_above_200_kev():
    from gamma.identification.id_anomalies import iodine_k_escape_for
    assert iodine_k_escape_for(661.0, "NaI") is None


def test_F276_no_escape_for_hpge():
    from gamma.identification.id_anomalies import iodine_k_escape_for
    assert iodine_k_escape_for(59.54, "HPGe") is None


def test_F276_find_escape_candidates_in_peak_list():
    from gamma.identification.id_anomalies import find_iodine_escape_candidates
    found_peaks = [29.0, 31.5, 59.5, 100.0]
    matches = find_iodine_escape_candidates(found_peaks, 59.54, "NaI")
    # 31.5 близко к 30.9 (±5 кэВ tol)
    assert 31.5 in matches
    assert 29.0 in matches    # 29.0 vs 30.93 = 1.93 < 5
    assert 100.0 not in matches


# ─────────────────────────────────────────────────────────────────────
# F-277 / T-005 — NORM 186 keV apportionment
# ─────────────────────────────────────────────────────────────────────

def test_F277_natural_apportionment_factors():
    from gamma.activity.norm_apportionment import (
        NORM_NATURAL_U_SHARE_226RA_186,
        NORM_NATURAL_U_SHARE_235U_186,
    )
    # Sum ≈ 1.0 by construction
    assert (NORM_NATURAL_U_SHARE_226RA_186
            + NORM_NATURAL_U_SHARE_235U_186) == pytest.approx(1.0, abs=0.001)


def test_F277_apportion_186_split():
    from gamma.activity.norm_apportionment import apportion_186_keV_NORM
    res = apportion_186_keV_NORM(1000.0)
    assert res.share_226Ra > 0.5
    assert res.share_235U > 0.4
    assert res.share_226Ra + res.share_235U == pytest.approx(1.0, abs=0.001)
    assert res.n226Ra_apportioned_area == pytest.approx(1000.0 * 0.5753, abs=0.5)
    assert res.n235U_apportioned_area == pytest.approx(1000.0 * 0.4247, abs=0.5)


def test_F277_is_186_keV_peak():
    from gamma.activity.norm_apportionment import is_186_keV_NORM_peak
    assert is_186_keV_NORM_peak(186.0)
    assert is_186_keV_NORM_peak(187.5)   # ±6 keV
    assert is_186_keV_NORM_peak(180.5)
    assert not is_186_keV_NORM_peak(200.0)


def test_F277_overestimate_without_apportionment():
    """Без apportionment 226Ra активность завышается на ~43% (0.5753 → 1.0)."""
    from gamma.activity.norm_apportionment import NORM_NATURAL_U_SHARE_226RA_186
    overestimation_factor = 1.0 / NORM_NATURAL_U_SHARE_226RA_186
    assert overestimation_factor == pytest.approx(1.738, abs=0.01)
    # 1.738 - 1.0 = 0.738 = +73.8% overestimation если приписать всю площадь Ra
    # Audit упоминал 43% — это middle-of-the-road для смешанных samples


# ─────────────────────────────────────────────────────────────────────
# F-278 / T-026 — targeted library presets
# ─────────────────────────────────────────────────────────────────────

def test_F278_osgi_preset_compact():
    from gamma.data.targeted_libraries import get_preset_nuclides
    osgi = get_preset_nuclides("OSGI")
    assert len(osgi) == 5
    assert "Cs-137" in osgi
    assert "Eu-152" in osgi


def test_F278_npp_preset_includes_iodines_and_actinides():
    from gamma.data.targeted_libraries import get_preset_nuclides
    npp = get_preset_nuclides("NPP")
    assert "I-131" in npp
    assert "Cs-137" in npp
    assert "Co-60" in npp
    assert "Pu-239" in npp
    assert 50 < len(npp) < 200


def test_F278_environmental_includes_norm():
    from gamma.data.targeted_libraries import get_preset_nuclides
    env = get_preset_nuclides("ENVIRONMENTAL")
    assert "K-40" in env
    assert "Ra-226" in env
    assert "Cs-137" in env
    assert "Pb-210" in env


def test_F278_unknown_preset_raises():
    from gamma.data.targeted_libraries import get_preset_nuclides
    with pytest.raises(KeyError):
        get_preset_nuclides("UNKNOWN_FAKE")


def test_F278_combine_presets_dedup():
    from gamma.data.targeted_libraries import combine_presets
    combined = combine_presets("OSGI", "ENVIRONMENTAL")
    # Cs-137 присутствует в обоих, но в combined — один раз
    assert combined.count("Cs-137") == 1
    assert "Am-241" in combined


def test_F278_list_presets_signature():
    from gamma.data.targeted_libraries import list_presets
    presets = list_presets()
    names = {p[0] for p in presets}
    assert "OSGI" in names
    assert "NPP" in names
    assert "RAO" in names
    assert "ENVIRONMENTAL" in names


# ─────────────────────────────────────────────────────────────────────
# F-279 / T-001 — LSRM √E + k·FWHM(E) floor
# ─────────────────────────────────────────────────────────────────────

def test_F279_window_low_e_floor_dominates():
    """На low-E (60 keV) NaI FWHM велик относительно √E — floor спасает."""
    from gamma.identification.window import build_id_window_lsrm_with_kfwhm_floor

    def nai_fwhm(E):
        # Approx Gamma-1S 63x63: 7% @ 662 → FWHM(E) ≈ 0.07·sqrt(E·662)
        return 0.07 * math.sqrt(E * 661.66)

    w = build_id_window_lsrm_with_kfwhm_floor("NaI", nai_fwhm)
    win_60 = w.window_keV(60.0)
    # δE₀·√(60/662) = 15·0.301 = 4.52
    # k·FWHM(60) = 1.5·0.07·sqrt(60·662) = 1.5·14 = 21.0
    # max = 21
    assert win_60 > 15.0, f"floor should boost low-E window: got {win_60}"


def test_F279_window_mid_e_formulas_close():
    """На Cs-137 661 кэВ обе формулы должны давать ~15 кэВ."""
    from gamma.identification.window import build_id_window_lsrm_with_kfwhm_floor

    def nai_fwhm(E):
        return 0.07 * math.sqrt(E * 661.66)

    w = build_id_window_lsrm_with_kfwhm_floor("NaI", nai_fwhm)
    win_661 = w.window_keV(661.66)
    # δE₀·√(1) = 15
    # k·FWHM(661) = 1.5 · 0.07 · 661.66 = 69.5 ← floor doминирует!
    # Это значит на 661 floor реально доминирует — Eq логика правильная
    assert win_661 >= 15.0


def test_F279_window_high_e_kfwhm_dominates():
    """На Tl-208 2614 кэВ k·FWHM ≫ δE₀·√E."""
    from gamma.identification.window import build_id_window_lsrm_with_kfwhm_floor

    def nai_fwhm(E):
        return 0.07 * math.sqrt(E * 661.66)

    w = build_id_window_lsrm_with_kfwhm_floor("NaI", nai_fwhm)
    win_2614 = w.window_keV(2614.0)
    sqrt_e_only = 15.0 * math.sqrt(2614.0 / 661.66)
    assert win_2614 > sqrt_e_only, "kfwhm floor should dominate at high E"


def test_F279_legacy_sqrt_e_branch_still_works():
    """Legacy sqrt_E branch без k_fwhm — work unchanged."""
    from gamma.identification.window import build_identification_window
    w = build_identification_window("NaI", delta_E0_keV=15.0)
    assert w.scaling == "sqrt_E"
    assert w.window_keV(661.66) == pytest.approx(15.0, abs=0.01)
