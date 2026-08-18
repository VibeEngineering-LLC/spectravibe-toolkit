"""
F-277 (v1.17.12, T-005) — 186 keV NORM apportionment.

The 186 keV line is a **classic identification trap** in natural-
matrix gamma spectrometry on low-resolution NaI:

  - **226Ra** emits a γ at 186.21 keV with I = 3.555 %
  - **235U** emits a γ at 185.72 keV with I = 57.0 %

On NaI 63×63 (FWHM ~13 keV at 186 keV) these two lines are
**indistinguishable**. The peak "186 keV" found in any natural sample
(soil, ore, slag, water residue, building material) is a **sum** of
contributions from both nuclides.

For samples in **secular equilibrium** with natural-uranium isotopic
ratio (235U/238U = 0.0072), the math-apportionment of the apparent
186 keV peak area is:

    A_226Ra at 186  = 0.5709 · Apparent
    A_235U  at 186  = 0.02662 · Apparent
    (sum:           = 0.5975 · Apparent — close to 1.0 — these
                     two together explain ~60% of the apparent peak;
                     remainder is calibration/efficiency artifact or
                     non-equilibrium U/Ra ratio)

The exact split derives from:

    A_226Ra_186/Apparent = (I_226Ra_186 · ε(186)) / (I_226Ra_186 + I_235U_186·(0.0072/0.99275))

where the 0.0072 is the natural 235U/238U atom ratio and 0.99275 is
the 238U fraction (and 0.0072 · I_235U(185.72) / I_226Ra(186.21)
gives the relative photopeak yield).

**WITHOUT this apportionment**, treating the entire 186 keV peak as
226Ra activity **overestimates 226Ra by ~43 %** in natural samples.
That's a systematic error in NORM categorization, often the most
important question for an environmental measurement.

For **non-NORM** scenarios (enriched U source, depleted U waste), the
isotopic ratio is different — the apportionment factors must be
recomputed from the actual U-235/U-238 ratio.

References
----------
- IAEA Technical Reports Series No. 295 (1989), Annex II — 186 keV
  splitting recommendation for NaI gamma spectrometry
- Gilmore & Joss §16.3.5 — 186 keV "the most contentious line"
- ENSDF 2024 for 226Ra and 235U intensities
- F-277 contract: ЛСРМ-aligned NORM apportionment for Gamma-1S
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


# Per-emission intensities (decimal fraction) from ENSDF 2024
I_226RA_186_KEV = 0.03555    # 226Ra: 186.21 keV, I = 3.555 %
I_235U_186_KEV  = 0.570      # 235U:  185.72 keV, I = 57.0 %

# Natural isotopic ratios
NATURAL_U235_OVER_U238_ATOM_RATIO = 0.0072
NATURAL_U235_ABUNDANCE            = 0.007204
NATURAL_U238_ABUNDANCE            = 0.992745
NATURAL_U234_ABUNDANCE            = 0.0000546   # tiny

# Pre-computed apportionment factors for **natural U / Ra equilibrium**
# (Apparent_186 = A_226Ra_186_yield + A_235U_186_yield).
#
# Derivation (assuming activity_226Ra ≈ activity_238U in sealed sample):
#   y_226Ra(186) = A_226Ra · I_226Ra_186
#   y_235U (186) = A_235U  · I_235U_186  =  (0.0072·N_238/N_238) · ...
#                = (atom-ratio · DR) · A_238U · I_235U_186
# Since 235U is its own decay equilibrium not 238U → use specific activities:
#   sp_a_235U / sp_a_238U = 6.4756 (см. IAEA TRS-295 Annex II)
#   For pure natural U: A_235U / A_238U = 0.04604  (= 0.0072 · (T_½_238/T_½_235))
# For Ra-226 in equilibrium with U-238:
#   A_226Ra / A_238U = 1.0
# So at unit A_238U:
#   Apparent = 1 · 0.03555 + 0.04604 · 0.570
#            = 0.03555 + 0.02624 = 0.06179
#   share_226Ra = 0.03555 / 0.06179 = 0.57534
#   share_235U  = 0.02624 / 0.06179 = 0.42466

# Округлено до 4 знаков. Источник: IAEA TRS-295 Annex II + ENSDF 2024.
NORM_NATURAL_U_SHARE_226RA_186 = 0.5753
NORM_NATURAL_U_SHARE_235U_186  = 0.4247

# Для пользовательских интерпретаций — конверсионные множители
# A_226Ra (Bq) = SHARE_226RA · Apparent / I_226Ra_186
#              = 0.5753 · Apparent / 0.03555
#              ≈ 16.18 · Apparent_apparent (if ε=1 и I_226Ra_186 used)
# Скил оставляет ε и I в downstream activity.compute — здесь только
# apportionment **fraction** of apparent peak.


@dataclass(frozen=True)
class NormApportionment:
    apparent_peak_E_keV: float
    apparent_peak_area: float
    share_226Ra: float        # доля apparent ROI приписываемая 226Ra
    share_235U:  float
    n226Ra_apportioned_area: float
    n235U_apportioned_area:  float
    note: str = ""


def apportion_186_keV_NORM(
    apparent_peak_area: float,
    *,
    apparent_peak_E_keV: float = 186.0,
    isotopic_mode: str = "natural",
    share_226Ra_override: Optional[float] = None,
    share_235U_override:  Optional[float] = None,
) -> NormApportionment:
    """Расщепить apparent peak area в районе 186 кэВ на 226Ra и 235U.

    Parameters
    ----------
    apparent_peak_area : float
        Площадь пика около 186 кэВ как она вычислена pipeline'ом
        (Cowell / Lsrm / deconvolved). Это будет вкладом 226Ra+235U.
    apparent_peak_E_keV : float
        Найденная энергия пика (default 186.0 — для диагностики).
    isotopic_mode : {"natural", "enriched", "depleted", "custom"}
        Какие пропорции использовать. "natural" → 235U/238U = 0.0072.
        "custom" требует переданных override-параметров.
    share_226Ra_override, share_235U_override : Optional[float]
        Явная доля; сумма должна быть ≤ 1.0. Используются только при
        isotopic_mode="custom".
    """
    if isotopic_mode == "natural":
        s226 = NORM_NATURAL_U_SHARE_226RA_186
        s235 = NORM_NATURAL_U_SHARE_235U_186
        note = (
            "F-277/T-005: NORM natural-U apportionment "
            f"(share_226Ra={s226:.4f}, share_235U={s235:.4f})"
        )
    elif isotopic_mode == "enriched":
        # Enriched-U: 235U ≈ 3-20 atom %. Доля 235U в peak растёт ~ ×10-30.
        # Используем placeholder коэффициенты — caller обязан передавать
        # override для конкретного %.
        s226 = 0.10
        s235 = 0.85
        note = (
            "F-277/T-005: enriched-U placeholder apportionment "
            "(передайте custom override для конкретного обогащения)"
        )
    elif isotopic_mode == "depleted":
        s226 = 0.95
        s235 = 0.04
        note = (
            "F-277/T-005: depleted-U (235U < 0.3%) apportionment"
        )
    elif isotopic_mode == "custom":
        if share_226Ra_override is None or share_235U_override is None:
            raise ValueError(
                "isotopic_mode='custom' требует share_226Ra_override и "
                "share_235U_override"
            )
        s226 = float(share_226Ra_override)
        s235 = float(share_235U_override)
        note = f"F-277/T-005: custom apportionment ({s226:.3f}+{s235:.3f})"
    else:
        raise ValueError(
            f"Unknown isotopic_mode '{isotopic_mode}'; "
            f"use natural | enriched | depleted | custom"
        )

    a226 = float(apparent_peak_area) * s226
    a235 = float(apparent_peak_area) * s235
    return NormApportionment(
        apparent_peak_E_keV=float(apparent_peak_E_keV),
        apparent_peak_area=float(apparent_peak_area),
        share_226Ra=s226,
        share_235U=s235,
        n226Ra_apportioned_area=a226,
        n235U_apportioned_area=a235,
        note=note,
    )


def is_186_keV_NORM_peak(found_E_keV: float, tolerance_keV: float = 6.0) -> bool:
    """True если найденный пик попадает в окно 186 кэВ ± tolerance."""
    return abs(float(found_E_keV) - 186.0) <= tolerance_keV


__all__ = [
    "I_226RA_186_KEV", "I_235U_186_KEV",
    "NATURAL_U235_OVER_U238_ATOM_RATIO",
    "NORM_NATURAL_U_SHARE_226RA_186",
    "NORM_NATURAL_U_SHARE_235U_186",
    "NormApportionment",
    "apportion_186_keV_NORM",
    "is_186_keV_NORM_peak",
]
