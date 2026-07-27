# -*- coding: utf-8 -*-
"""F-330 / v1.18.18.4 — LSRM .spe COMMENT passport parser.

Извлекает паспортные активности из поля COMMENT файлов калибровочных
источников. Закрывает разрыв «оператор должен вручную ввести данные
паспорта», когда они уже физически присутствуют в .spe.

Контракт:
- 4 real Gamma-1S фикстуры (M_cs / M_k / M_ra / M_th) дают по одному
  entry каждая.
- Cyrillic К → Latin K в имени нуклида.
- Reference date парсится в date(YYYY, MM, DD) с эпохой-эвристикой
  (97 → 1997, 25 → 2025).
- Decay correction для Cs-137 (t½ = 30.07 y): 1997-05-30 → 1999-08-04
  даёт фактор ≈0.9510 (расчётно).
- Long-lived nuclides (K-40, Ra-226, Th-232) имеют known half-lives,
  но фактор decay-correction ≈1 на масштабах years.
"""
from __future__ import annotations
import sys
from datetime import date
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from gamma.io.lsrm_passport import (
    parse_lsrm_passport_comment,
    decay_correct,
    half_life_seconds,
    _canonicalize_nuclide,
    _normalize_nuclide_symbol,
)


# ─── Real-fixture parsing ───────────────────────────────────────────

REAL_COMMENTS = {
    "M_cs": "(Sum/T =    70.67 Sum*I/T =  8219.98) Cs-137 - 1890 Бк/кг (5%) от 30.05.97г.",
    "M_k":  "(Sum/T =    21.56 Sum*I/T =  3770.29) К-40 - 2540 Бк/кг (10%)",
    "M_ra": "(Sum/T =   189.30 Sum*I/T = 21016.32)  Ra-226 - 1850 Бк/кг(10%)",
    "M_th": "(Sum/T =   277.48 Sum*I/T = 29649.73) Th-232 - 2200 Бк/кг(5%)",
}


def test_F330_M_cs_full_parse():
    """M_cs: Cs-137 - 1890 Бк/кг (5%) от 30.05.97г."""
    es = parse_lsrm_passport_comment(REAL_COMMENTS["M_cs"])
    assert len(es) == 1
    e = es[0]
    assert e.nuclide == "Cs-137"
    assert e.value == 1890.0
    assert e.unit == "Бк/кг"
    assert e.uncertainty_pct == 5.0
    assert e.reference_date == date(1997, 5, 30)
    assert e.is_specific_activity is True


def test_F330_M_k_cyrillic_K_normalized():
    """M_k: К-40 (Cyrillic К) → K-40 (Latin K)."""
    es = parse_lsrm_passport_comment(REAL_COMMENTS["M_k"])
    assert len(es) == 1
    assert es[0].nuclide == "K-40"
    assert es[0].uncertainty_pct == 10.0
    assert es[0].reference_date is None  # no «от ...»


def test_F330_M_ra_no_space_before_paren():
    """M_ra: Бк/кг(10%) — no space перед скобкой — must still match."""
    es = parse_lsrm_passport_comment(REAL_COMMENTS["M_ra"])
    assert len(es) == 1
    assert es[0].nuclide == "Ra-226"
    assert es[0].value == 1850.0


def test_F330_M_th_full_parse():
    """M_th: Th-232 - 2200 Бк/кг(5%)."""
    es = parse_lsrm_passport_comment(REAL_COMMENTS["M_th"])
    assert len(es) == 1
    e = es[0]
    assert e.nuclide == "Th-232"
    assert e.value == 2200.0
    assert e.uncertainty_pct == 5.0


# ─── Edge cases ─────────────────────────────────────────────────────

def test_F330_empty_comment():
    assert parse_lsrm_passport_comment("") == []


def test_F330_only_accounting_no_passport():
    """COMMENT с только Sum/T бухгалтерией — без passport-данных."""
    es = parse_lsrm_passport_comment("(Sum/T = 12.34 Sum*I/T = 56.78)")
    assert es == []


def test_F330_decimal_comma():
    """Россиская локаль: «1890,5» вместо «1890.5»."""
    es = parse_lsrm_passport_comment("Cs-137 - 1890,5 Бк/кг (5%)")
    assert es[0].value == 1890.5


def test_F330_decimal_dot():
    es = parse_lsrm_passport_comment("Cs-137 - 1890.5 Бк/кг (5%)")
    assert es[0].value == 1890.5


def test_F330_scientific_notation():
    es = parse_lsrm_passport_comment("Co-60 - 1.05e3 Бк (5%)")
    assert es[0].value == pytest.approx(1050.0)
    assert es[0].is_specific_activity is False


def test_F330_kbq_kilo_prefix():
    es = parse_lsrm_passport_comment("Cs-137 - 1.890 kBq/kg (5%)")
    assert len(es) == 1
    # kBq prefix → multiplier 1000
    assert es[0].value_Bq_per_kg() == pytest.approx(1890.0)


def test_F330_year_two_digit_19xx():
    """Двухзначный год 97 → 1997."""
    es = parse_lsrm_passport_comment("Cs-137 - 1000 Бк (5%) от 30.05.97")
    assert es[0].reference_date == date(1997, 5, 30)


def test_F330_year_two_digit_20xx():
    """Двухзначный год 25 → 2025."""
    es = parse_lsrm_passport_comment("Cs-137 - 1000 Бк (5%) от 30.05.25")
    assert es[0].reference_date == date(2025, 5, 30)


def test_F330_year_four_digit():
    es = parse_lsrm_passport_comment("Cs-137 - 1000 Бк (5%) от 30.05.1997")
    assert es[0].reference_date == date(1997, 5, 30)


def test_F330_metastable_nuclide():
    """Tc-99m, Ba-137m — суффикс m."""
    es = parse_lsrm_passport_comment("Tc-99m - 500 Бк (5%)")
    assert es[0].nuclide == "Tc-99m"


def test_F330_multinuclide():
    """Один COMMENT с двумя нуклидами."""
    c = "Cs-137 - 1000 Бк/кг (5%) от 01.01.20; Co-60 - 800 Бк/кг (3%) от 01.01.20"
    es = parse_lsrm_passport_comment(c)
    assert len(es) == 2
    nuclides = sorted(e.nuclide for e in es)
    assert nuclides == ["Co-60", "Cs-137"]


# ─── Normalization helpers ──────────────────────────────────────────

def test_F330_cyrillic_normalize_K():
    assert _normalize_nuclide_symbol("К-40") == "K-40"


def test_F330_cyrillic_normalize_C():
    assert _normalize_nuclide_symbol("С-137") == "C-137"


def test_F330_canonicalize_uppercase():
    """«CS-137» → «Cs-137»."""
    assert _canonicalize_nuclide("CS-137") == "Cs-137"


def test_F330_canonicalize_single_letter():
    assert _canonicalize_nuclide("K-40") == "K-40"


def test_F330_canonicalize_metastable():
    assert _canonicalize_nuclide("TC-99M") == "Tc-99m"


# ─── Decay correction ──────────────────────────────────────────────

def test_F330_half_life_cs137():
    """Cs-137 t½ = 30.07 y = 9.4894e8 s."""
    T = half_life_seconds("Cs-137")
    assert T is not None
    assert T == pytest.approx(9.4894e8, rel=1e-4)


def test_F330_half_life_unknown_returns_none():
    assert half_life_seconds("Xx-999") is None


def test_F330_decay_cs137_1997_to_1999():
    """Cs-137 1890 Bq в 1997-05-30 → 1999-08-04 ≈ 1798 Bq.

    Δt = 2.181 y. Fraction = exp(-ln(2)·2.181/30.07) = 0.9511.
    1890 × 0.9511 ≈ 1797.6 Bq.
    """
    A_corr = decay_correct(1890.0, "Cs-137", date(1997, 5, 30), date(1999, 8, 4))
    assert A_corr == pytest.approx(1797.6, rel=2e-3)


def test_F330_decay_k40_negligible():
    """K-40 t½ = 1.28e9 y → over 50 years correction <0.001%."""
    A_corr = decay_correct(1000.0, "K-40", date(1970, 1, 1), date(2026, 1, 1))
    assert A_corr is not None
    assert A_corr == pytest.approx(1000.0, rel=1e-5)


def test_F330_decay_unknown_nuclide_returns_none():
    """Unknown nuclide → None (no half-life in table)."""
    A_corr = decay_correct(1000.0, "Xx-999", date(2020, 1, 1), date(2025, 1, 1))
    assert A_corr is None


def test_F330_decay_zero_dt_returns_A0():
    A_corr = decay_correct(1000.0, "Cs-137", date(2020, 1, 1), date(2020, 1, 1))
    assert A_corr == 1000.0


# ─── BUG-49: Поверка-2016 lab COMMENT format ─────────────────────────
#
# Operator lab «Поверка 2016» (Gamma-1S, 22 Tier-1 fixtures verified
# via probe 2026-06-04: _state/agent_a/outbox/2026-06-04_lsrm_spe_descriptions.json)
# writes passport data в COMMENT в формате:
#   «<Nuc> A=<v> <Бк|Бк/кг> dA=<u>% DD-MM-YYYY»
# Production F-330 regex returns 0/22 (silent failure); the extended
# parser below matches 22/22.
#
# Each test case cites the source .spe fixture by shifr / filename token.


def test_BUG49_poverka_2016_Am_241_Bq():
    """Source #42.13: «Am-241 A=118000 Бк dA=5% 03-12-2013».

    Probe JSON entry: per_file[…].comment_raw для
    <LSRM>\\Work\\BG\\Gamma-1S\\Spe - поверки\\Поверка 2016\\Точка 5см\\
    Am-241 42.13_<sample>_5cm.spe.
    """
    es = parse_lsrm_passport_comment("Am-241 A=118000 Бк dA=5% 03-12-2013")
    assert len(es) == 1
    e = es[0]
    assert e.nuclide == "Am-241"
    assert e.value == 118000.0
    assert e.unit == "Бк"
    assert e.uncertainty_pct == 5.0
    assert e.reference_date == date(2013, 12, 3)
    assert e.is_specific_activity is False
    assert e.value_Bq() == 118000.0
    assert e.value_Bq_per_kg() is None


def test_BUG49_poverka_2016_Cs_137_Bq_per_kg():
    """Marinelli 420-7-14 (Cs-137_420-7-14_<sample>_0cm.spe):
    «Cs-137 A=1760 Бк/кг dA=5% 24-05-2002» — mass-specific activity for
    filled-volume Marinelli geometry."""
    es = parse_lsrm_passport_comment("Cs-137 A=1760 Бк/кг dA=5% 24-05-2002")
    assert len(es) == 1
    e = es[0]
    assert e.nuclide == "Cs-137"
    assert e.value == 1760.0
    assert e.unit == "Бк/кг"
    assert e.uncertainty_pct == 5.0
    assert e.reference_date == date(2002, 5, 24)
    assert e.is_specific_activity is True
    assert e.value_Bq() is None
    assert e.value_Bq_per_kg() == 1760.0


def test_BUG49_poverka_2016_e_notation():
    """Source #SRC-05 Cd-109 (Cd-109 #SRC-05_<sample>_5cm.spe):
    «Cd-109 A=1.033E6 Бк dA=2% 01-10-2008» — scientific notation in
    value field must be parsed."""
    es = parse_lsrm_passport_comment("Cd-109 A=1.033E6 Бк dA=2% 01-10-2008")
    assert len(es) == 1
    e = es[0]
    assert e.nuclide == "Cd-109"
    assert e.value == pytest.approx(1.033e6)
    assert e.reference_date == date(2008, 10, 1)
    assert e.value_Bq() == pytest.approx(1.033e6)


def test_BUG49_poverka_2016_multiline_source_SRC-05():
    """Source #SRC-05 calibration-source set: 7 nuclides distributed
    across 7 .spe files (one per nuclide), each with one-line COMMENT.

    Concatenating their COMMENT contents emulates a multi-nuclide
    COMMENT, which the parser must split into 7 entries. Each line's
    canonical value is from probe outbox JSON entries:
      Ba-133 A=44100, Cd-109 A=1.033E6, Co-57 A=99500,
      Co-60 A=107800, Cs-137 A=94200, Eu-152 A=46700, Th-228 A=37700.
    All dated 2008-10-01, all dA=2%.
    """
    multi = (
        "Ba-133 A=44100 Бк dA=2% 01-10-2008\n"
        "Cd-109 A=1.033E6 Бк dA=2% 01-10-2008\n"
        "Co-57 A=99500 Бк dA=2% 01-10-2008\n"
        "Co-60 A=107800 Бк dA=2% 01-10-2008\n"
        "Cs-137 A=94200 Бк dA=2% 01-10-2008\n"
        "Eu-152 A=46700 Бк dA=2% 01-10-2008\n"
        "Th-228 A=37700 Бк dA=2% 01-10-2008\n"
    )
    es = parse_lsrm_passport_comment(multi)
    assert len(es) == 7
    nucs = sorted(e.nuclide for e in es)
    assert nucs == [
        "Ba-133", "Cd-109", "Co-57", "Co-60",
        "Cs-137", "Eu-152", "Th-228",
    ]
    # All same reference date
    assert all(e.reference_date == date(2008, 10, 1) for e in es)
    # All Бк (volumetric), uncertainty 2%
    assert all(e.unit == "Бк" for e in es)
    assert all(e.uncertainty_pct == 2.0 for e in es)
    # Spot-check one value
    cd109 = next(e for e in es if e.nuclide == "Cd-109")
    assert cd109.value == pytest.approx(1.033e6)


def test_BUG49_poverka_2016_missing_date_no_match():
    """Поверка-2016 format requires DD-MM-YYYY suffix. If date is
    missing (malformed COMMENT), the new regex does NOT match — the
    fallback parser is strict about the full line shape, mirroring
    F-330 «skip-malformed-silently» convention. Returns []."""
    # No date suffix → fallback regex requires the date group, so no match.
    es = parse_lsrm_passport_comment("Am-241 A=118000 Бк dA=5%")
    # Note: F-330 canonical regex also won't match (no dash, no «-»).
    # Expected behaviour: no entries (skip malformed silently).
    assert es == []


def test_BUG49_poverka_2016_K40_cyrillic_Bq_per_kg():
    """Marinelli K40_420-7-20 fixture:
    «K-40 A=2530 Бк/кг dA=6% 24-05-2002».

    Verifies that the Latin «K-40» nuclide symbol parses correctly with
    Бк/кг unit and date 2002-05-24. Probe JSON entry confirms
    activity_Bq_per_kg = 2530.0."""
    es = parse_lsrm_passport_comment("K-40 A=2530 Бк/кг dA=6% 24-05-2002")
    assert len(es) == 1
    e = es[0]
    assert e.nuclide == "K-40"
    assert e.value_Bq_per_kg() == 2530.0
    assert e.reference_date == date(2002, 5, 24)
    assert e.uncertainty_pct == 6.0


def test_BUG49_poverka_2016_cyrillic_K_normalized():
    """К-40 (Cyrillic К) in Поверка-2016 format → K-40 (Latin K).

    Defensive: probe data shows operator lab uses Latin K, but the
    parser must remain robust to Cyrillic К look-alike (matches behaviour
    of canonical F-330 path)."""
    es = parse_lsrm_passport_comment("К-40 A=2530 Бк/кг dA=6% 24-05-2002")
    assert len(es) == 1
    assert es[0].nuclide == "K-40"


def test_BUG49_mixed_canonical_and_poverka_2016():
    """COMMENT containing both canonical F-330 and Поверка-2016 lines
    must parse both without double-counting. Overlap-detection guard
    in parser ensures each text span is owned by exactly one regex."""
    mixed = (
        "Cs-137 - 1890 Бк/кг (5%) от 30.05.97г. "
        "Co-60 A=107800 Бк dA=2% 01-10-2008"
    )
    es = parse_lsrm_passport_comment(mixed)
    nucs = sorted(e.nuclide for e in es)
    assert nucs == ["Co-60", "Cs-137"]
    cs = next(e for e in es if e.nuclide == "Cs-137")
    co = next(e for e in es if e.nuclide == "Co-60")
    assert cs.reference_date == date(1997, 5, 30)
    assert cs.value == 1890.0
    assert co.reference_date == date(2008, 10, 1)
    assert co.value == 107800.0


def test_BUG49_regression_canonical_F330_still_passes():
    """Regression guard: canonical F-330 «Nuc - val unit (unc%) от
    date» format must continue to parse identically after BUG-49
    extension. Same fixture as test_F330_M_cs_full_parse."""
    es = parse_lsrm_passport_comment(
        "(Sum/T = 70.67 Sum*I/T = 8219.98) "
        "Cs-137 - 1890 Бк/кг (5%) от 30.05.97г."
    )
    assert len(es) == 1
    e = es[0]
    assert e.nuclide == "Cs-137"
    assert e.value == 1890.0
    assert e.unit == "Бк/кг"
    assert e.uncertainty_pct == 5.0
    assert e.reference_date == date(1997, 5, 30)
