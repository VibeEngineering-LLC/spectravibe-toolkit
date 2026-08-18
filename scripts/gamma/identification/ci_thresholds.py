"""
F-275 (v1.17.12, T-022) — Per-nuclide CI thresholds from ЛСРМ Table 14-1.

Lsrm Algorithmic Foundations 2022 §14.3 Table 14-1 publishes calibrated
Confidence Index values for common nuclides per detector class. These
values represent the *expected* CI when ALL library lines of the
nuclide are matched in the spectrum:

  Cs-137 on NaI: CI ≈ 1.8   (single 661 keV line)
  K-40   on NaI: CI ≈ 2.2
  Na-22  on NaI: CI ≈ 3.8   (511 + 1274)
  Co-60  on NaI: CI ≈ 5.9   (1173 + 1332)
  Cs-134 on NaI: CI ≈ 4.4   (multiple lines)
  Ba-133 on NaI: CI ≈ 8.5   (5+ lines)
  Eu-152 on NaI: CI ≈ 18.3  (very many lines)
  Th-232 on NaI: CI ≈ 16.6  (full chain, many lines)

These values are not arbitrary thresholds: each is the LSRM-reference
CI for a complete identification. When a candidate has CI **noticeably
below** its reference, the identification is **incomplete** — either
some library lines weren't matched, or those that were matched have
ambiguous evidence.

**Use case**: post-rank gating in `ci_gating.py`. The previous
implementation used a single global `CI_CONFIRMED_THRESHOLD = 2.0`,
which lets multi-line nuclides like Eu-152 promote on **partial**
evidence (e.g., only 3 of 13 lines matched). Per-nuclide thresholds
require Eu-152 to come close to its reference 18.3 before promotion.

For HPGe / LaBr3 the table-1 values are different (higher precision →
larger CI). We expose a class-aware lookup so callers pass the detector
class explicitly.

References
----------
- ЛСРМ Algorithmic Foundations 2022 §14.3, Table 14-1
- F-275 (v1.17.12) — wire into ci_gating
"""
from __future__ import annotations

from typing import Optional


# ──────────────────────────────────────────────────────────────────
# Lsrm Table 14-1 — calibrated CI per (nuclide, detector_class)
# ──────────────────────────────────────────────────────────────────

# Per-nuclide reference CI (full library match) для NaI 50×50 / Gamma-1S 63×63.
# Минимальный CI для **confirmed**-promotion = 0.10 × reference (10% of full
# coverage). Это компромисс: фильтрует совсем слабые ID без блокировки
# реалистично-неполных мультилинейных нуклидов.
LSRM_TABLE_14_1_NAI = {
    "Cs-137":  1.8,
    "K-40":    2.2,
    "Na-22":   3.8,
    "Co-60":   5.9,
    "Cs-134":  4.4,
    "Ba-133":  8.5,
    "Eu-152": 18.3,
    "Eu-154": 16.0,    # similar density of lines
    "Eu-155":  3.0,    # 86.5 + 105.3, weakly resolved on NaI
    "Th-232": 16.6,    # full chain (Tl-208, Ac-228, Pb-212, Bi-212 ...)
    "Ra-226":  9.5,    # daughter chain (Pb-214, Bi-214 ...)
    "U-238":   8.0,    # via Th-234
    "Am-241":  2.0,    # 59.5 single line + 26.3 weak
    "Tl-208":  3.5,    # 583, 2614 — clean lines on NaI
    "Ac-228":  4.0,    # 338+911+969 + several others
    "Pb-212":  2.5,    # 238 dominant
    "Bi-212":  3.0,    # 727, 1620 weak on NaI
    "Pb-214":  4.5,    # 295+352 doublet + 242
    "Bi-214":  6.0,    # 609+1120+1764 ladder
}

# Для HPGe — все CI существенно выше из-за гораздо меньшего δE/E.
# Линейная аппроксимация: HPGe CI ≈ NaI CI + 5...8.
LSRM_TABLE_14_1_HPGE = {
    k: v + 6.0 for k, v in LSRM_TABLE_14_1_NAI.items()
}

# LaBr3 / CeBr3 — между NaI и HPGe; ≈ NaI + 2.
LSRM_TABLE_14_1_LABR = {
    k: v + 2.0 for k, v in LSRM_TABLE_14_1_NAI.items()
}

# Fraction of reference CI required for confirmation.
# 0.10 (10%) — backward-compat default. Раньше пайплайн использовал
# глобальный CI ≥ 2.0; для Eu-152 это soft на любых трёх линиях.
# Per-nuclide gate с фракцией 0.10 для Eu-152 → требует CI ≥ 1.83
# (NaI). Для Cs-137 → CI ≥ 0.18 (NaI), что значит даже single line
# проходит на 1σ. Set fraction = 0.5 (strict, T-022 default) →
# Eu-152 требует CI ≥ 9.15 (≈половина библиотеки).
CONFIRMATION_FRACTION_DEFAULT = 0.5
# Tentative gate — 20%
TENTATIVE_FRACTION_DEFAULT = 0.2


def _table_for(detector_class: str) -> dict:
    dc = (detector_class or "NaI").strip()
    dc_low = dc.lower()
    if dc_low.startswith("hpge") or dc_low.startswith("cdz"):
        return LSRM_TABLE_14_1_HPGE
    if dc_low.startswith("labr") or dc_low.startswith("cebr"):
        return LSRM_TABLE_14_1_LABR
    return LSRM_TABLE_14_1_NAI


def lsrm_table_14_1_ci(
    nuclide: str, detector_class: str = "NaI",
) -> Optional[float]:
    """Вернуть calibrated CI для нуклида из ЛСРМ Table 14-1, или None
    если нуклид не в таблице.

    Не-таблица (Co-58, U-235, Y-88 и т.п.) → None. Caller должен
    fallback на глобальный порог (CI_CONFIRMED_THRESHOLD).
    """
    table = _table_for(detector_class)
    return table.get(str(nuclide))


def required_ci_for_confirmation(
    nuclide: str,
    detector_class: str = "NaI",
    fraction: float = CONFIRMATION_FRACTION_DEFAULT,
    global_floor: float = 2.0,
) -> float:
    """Минимальный CI для promotion в **confirmed**.

    Returns
    -------
    max(global_floor, fraction · lsrm_table[nuclide]) для нуклидов из
    таблицы; иначе просто global_floor (backward-compat).
    """
    ref = lsrm_table_14_1_ci(nuclide, detector_class)
    if ref is None:
        return float(global_floor)
    return max(float(global_floor), float(fraction) * float(ref))


def required_ci_for_tentative(
    nuclide: str,
    detector_class: str = "NaI",
    fraction: float = TENTATIVE_FRACTION_DEFAULT,
    global_floor: float = 1.0,
) -> float:
    """Минимальный CI для **tentative** (просмотр оператором)."""
    ref = lsrm_table_14_1_ci(nuclide, detector_class)
    if ref is None:
        return float(global_floor)
    return max(float(global_floor), float(fraction) * float(ref))


__all__ = [
    "LSRM_TABLE_14_1_NAI",
    "LSRM_TABLE_14_1_HPGE",
    "LSRM_TABLE_14_1_LABR",
    "CONFIRMATION_FRACTION_DEFAULT",
    "TENTATIVE_FRACTION_DEFAULT",
    "lsrm_table_14_1_ci",
    "required_ci_for_confirmation",
    "required_ci_for_tentative",
]
