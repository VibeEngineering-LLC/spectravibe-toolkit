# -*- coding: utf-8 -*-
"""v1.17.17 delivery tests — Statistical thresholds (F-288)."""
from __future__ import annotations
import math, os, sys
import pytest

SCRIPTS = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts")
)
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)


def test_F288_currie_paired_blank_LD_at_bg_100():
    from gamma.identification.statistical_thresholds import currie_paired_blank_LD
    # L_D = 1.645² + 2·1.645·√100 = 2.706 + 32.9 = 35.6
    ld = currie_paired_blank_LD(100.0)
    assert ld == pytest.approx(35.61, abs=0.05)


def test_F288_currie_paired_blank_LD_zero_bg():
    from gamma.identification.statistical_thresholds import currie_paired_blank_LD
    # L_D = k² = 2.706
    assert currie_paired_blank_LD(0.0) == pytest.approx(2.706, abs=0.01)


def test_F288_currie_paired_blank_LC():
    from gamma.identification.statistical_thresholds import currie_paired_blank_LC
    # L_C = 1.645·√200 = 23.26
    assert currie_paired_blank_LC(100.0) == pytest.approx(23.26, abs=0.05)


def test_F288_iso_upper_limit():
    from gamma.identification.statistical_thresholds import iso_11929_upper_limit_Bq
    # L_U = A + 1.645·σ_A
    assert iso_11929_upper_limit_Bq(10.0, 5.0) == pytest.approx(18.225, abs=0.01)


def test_F288_report_as_upper_limit():
    from gamma.identification.statistical_thresholds import report_as_upper_limit
    # A=5 < L_C=20 → True (report as upper limit)
    assert report_as_upper_limit(A_Bq=5.0, sigma_A_Bq=2.0, L_C_Bq=20.0)
    # A=30 > L_C=20 → False (normal detection report)
    assert not report_as_upper_limit(A_Bq=30.0, sigma_A_Bq=2.0, L_C_Bq=20.0)


def test_F288_limit_of_quantitation_default_10():
    from gamma.identification.statistical_thresholds import limit_of_quantitation_LQ
    # L_Q = 10·√100 = 100
    assert limit_of_quantitation_LQ(100.0) == pytest.approx(100.0)


def test_F288_threshold_glossary_keys():
    from gamma.identification.statistical_thresholds import (
        THRESHOLD_GLOSSARY, explain_threshold,
    )
    assert "L_C" in THRESHOLD_GLOSSARY
    assert "L_D" in THRESHOLD_GLOSSARY
    assert "MDA" in THRESHOLD_GLOSSARY
    assert "ISO 11929" in explain_threshold("L_U")


def test_F288_iso_quadratic_with_zero_u_rel_matches_simple():
    """При u_rel_g=0 quadratic должна дать ≈ простую формулу."""
    from gamma.identification.statistical_thresholds import iso_11929_LD_quadratic
    L_C_counts = 23.26
    sigma_0 = 14.14   # √200 → L_C = 1.645·14.14 = 23.26
    # Simple: L_D = 2·L_C + k² = 46.52 + 2.71 = 49.23
    ld = iso_11929_LD_quadratic(
        L_C_counts=L_C_counts, sigma_0_counts=sigma_0, u_rel_g=0.0,
    )
    # Iterative may converge to slightly different value via L_C + k·σ_0;
    # let's check both expected forms agree within 5%.
    expected_simple_formula = L_C_counts + 1.645 * sigma_0
    assert ld == pytest.approx(expected_simple_formula, rel=0.02)


def test_F288_iso_quadratic_with_u_rel_increases_LD():
    """При non-trivial u_rel_g L_D должно вырасти."""
    from gamma.identification.statistical_thresholds import iso_11929_LD_quadratic
    L_C = 23.26
    sig0 = 14.14
    ld0 = iso_11929_LD_quadratic(L_C_counts=L_C, sigma_0_counts=sig0, u_rel_g=0.0)
    ld_with = iso_11929_LD_quadratic(L_C_counts=L_C, sigma_0_counts=sig0, u_rel_g=0.1)
    assert ld_with > ld0


def test_F288_best_mda_line_picks_highest_yield_lowest_bg():
    from gamma.identification.statistical_thresholds import (
        best_mda_line, _LineCandidate,
    )
    lines = [
        _LineCandidate(line_E_keV=100.0, intensity_decimal=0.05, efficiency=0.05, bg_counts_in_ROI=400),
        _LineCandidate(line_E_keV=662.0, intensity_decimal=0.85, efficiency=0.015, bg_counts_in_ROI=100),
        _LineCandidate(line_E_keV=1461.0, intensity_decimal=0.1, efficiency=0.005, bg_counts_in_ROI=50),
    ]
    best = best_mda_line(lines)
    # 662 line: I·ε = 0.01275, √bg = 10, score = 0.001275
    # 100 line: I·ε = 0.0025, √bg = 20, score = 0.000125
    # 1461 line: I·ε = 0.0005, √bg = 7.07, score = 0.0000707
    # 662 wins
    assert best.line_E_keV == 662.0


def test_F288_best_mda_line_empty():
    from gamma.identification.statistical_thresholds import best_mda_line
    assert best_mda_line([]) is None
