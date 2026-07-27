"""
Staged candidate-nuclide sets for two-pass identification (F-70 / v1.11.1).

Identification is run in three stages, from most-likely to least-likely:

  Stage 1 — ЕРН (естественные радионуклиды):
    • K-40 (always)
    • Th-232 chain proxies: Ac-228, Pb-212, Tl-208, Bi-212
    • Ra-226 chain proxies: Pb-214, Bi-214, Ra-226 (186)
    • U-238 chain proxies: Th-234, Pa-234m
    • Other naturals: U-235, Pb-210, Th-228, Ra-224, Ra-228
    • Plus XRF lines from Pb shielding and Th matrix (handled in
      `gamma.identification.xrf_step` separately)

    Methodology reference: NOTES_v1.7_methodology.md §9 — Lsrm 7-line
    NaI calibration set on ЕРН/Фон-спектрах.

  Stage 2 — Техногенные ("Chernobyl set"):
    • Cs-137, Cs-134 (fallout, Chernobyl, fuel-cycle contamination)
    • Co-60 (industrial / cobalt sources)
    • I-131 (medical / fresh fallout)
    • Sr-90 (β-emitter — only detectable via bremsstrahlung continuum,
      not via FEP; included for diagnostics)

  Stage 3 — Прочее (medical, calibration, exotic):
    • Na-22, Be-7 (calibration sources or rare contamination)
    • Am-241, Eu-152, Ba-133 (laboratory calibration sources)
    • Mn-54, Zn-65, Co-57, Tc-99m, Ga-67, In-111 (medical/industrial)

Per user methodology (15.11.2025):
  - Stage 2 and Stage 3 are only invoked when Stage 1 leaves significant
    unidentified peaks or multiplet broadening that the ЕРН set cannot
    explain.
  - Better to ask the user than fantasize: if Stage 1 explains most of
    the spectrum, do NOT spawn Stage 2/3 without user confirmation.
  - Annihilation 511 keV peak in Th-232-rich samples is dominated by
    Tl-208 510.77 (I≈22.6% after Bi-212 branching); Na-22 should not be
    proposed unless 1274.5 keV is also confirmed with correct ratio.
"""

from __future__ import annotations

from typing import List, Optional


# --- Stage 1: ЕРН (natural radionuclides) ----------------------------
# Order: dominant chain proxies first so identify_nuclides processes them
# before weaker chain members that share photopeaks on NaI resolution.
ERN_STAGE1: List[str] = [
    # K-40 — universal
    "K-40",
    # Th-232 chain dominant proxies
    "Tl-208",  # 2614.51 + 583.19 + 510.77 + 860.56 + 277.36 (I=35.85/30.6/8.1)
    "Ac-228",  # 911.20 + 968.97 + 338.32 + 463.0 + 794.95 + 209.25 (I=25.8/15.8/11.3)
    "Pb-212",  # 238.63 (I=43.6) + 300.09 — dominant 240-keV NaI peak in Th
    "Bi-212",  # 727.33 (I=6.67) + 1620.50 + 785.37
    # Ra-226 / U-238 chain dominant proxies
    "Bi-214",  # 609.31 (I=45.49) + 1120 + 1764 + 2204 + 1238 + 768
    "Pb-214",  # 351.93 (I=35.6) + 295.22 + 241.98 (I=7.43) + 53.23
    "Ra-226",  # 186.21 (I=3.59) — weak, often blended with U-235 185.71
    # U-238 chain — uranium-indicating daughters
    "Th-234",  # 63.3 + 92.4 + 92.8 (weak, low-energy XRF region)
    "Pa-234m", # 1001.03
    # Other natural / quasi-natural
    "Pb-210",  # 46.5 — independent (lead shielding contamination)
    "U-235",   # 143.76 + 185.71 + 163.36 + 205.31 — only on U-enriched samples
    # ── Intentionally EXCLUDED from default Stage 1 (F-73a / v1.11.1) ──
    # Th-228 (84.4 keV, I=1.22%) — dominated by Pb-XRF triplet 72.8/75.0/84.4
    #   on NaI 63×63 (FWHM ~7 keV at 80 keV). Add manually if testing a
    #   pure Th-228 sealed source.
    # Ra-224 (241.0 keV, I=4.1%; 81.0 keV, I=1.27%) — its 241 ROI is
    #   physically 91% Pb-212 (I=43.6 vs 4.1) in any sample where both
    #   are in equilibrium. On NaI, including both in the candidate set
    #   makes disambiguate Rule 3 incorrectly award the peak to Ra-224
    #   (which has 2 weak matches over Pb-212's 1 strong). Add only
    #   when the user specifically needs Ra-224 quantification.
]


# --- Stage 2: Технические/чернобыльские ------------------------------
TECHNOGENIC_STAGE2: List[str] = [
    "Cs-137",  # 661.66 — Chernobyl, weapon fallout, medical
    "Cs-134",  # 604.72 + 795.86 + 569.33 + 801.95 + 1365.19
    "Co-60",   # 1173.23 + 1332.49
    "I-131",   # 364.49 + 636.99
    # Sr-90 is a pure β emitter — not detectable as FEP
    # Y-90 same — only bremsstrahlung
]


# --- Stage 3: Экзотика / медицина / прочее ---------------------------
# Per user policy: should NEVER be auto-proposed; require explicit user
# confirmation. Add 511 keV here so that Tl-208 wins by Lsrm Rule 2 in
# `disambiguate` — but Na-22 is still available as a candidate if the
# user explicitly enables Stage 3.
EXOTIC_STAGE3: List[str] = [
    "Na-22",   # 511 + 1274.54
    "Be-7",    # 477.6 — cosmogenic, air filters / fresh vegetation
    "Am-241",  # 59.54 + 26.34 — calibration / smoke detector / radwaste
    "Eu-152",  # multi-line — calibration source
    "Ba-133",  # multi-line — calibration source
    "Mn-54",   # 834.85
    "Zn-65",   # 1115.55
    "Co-57",   # 122.06 + 136.47
    "Tc-99m",  # 140.51 — medical
    "Ga-67",   # 93.31 + 184.58 + 300.22 + 393.53 — medical
    "In-111",  # 171.28 + 245.40 — medical
    # BUG-36 (Wave 3 C1, 2026-06-05) — Ti-44 + Sc-44 calibration set.
    # Ti-44 is a long-lived (T½ = 60 y) β+ EC calibration source whose
    # short-lived daughter Sc-44 (T½ = 3.97 h) is in secular equilibrium
    # in any aged Ti-44 source. The AmTiCsEu mixed-calibration Marinelli
    # fixture (`detectors/Gamma-1S/reference_spectra/archive/Поверка-2016/`)
    # carries both. Library entries are present in `data/nuclides.json`
    # (lines 1514–1552) but identification missed them prior to this
    # fix because neither nuclide was in any stage candidate set.
    # Both are tagged `tcs_type=calibration_source` / `chain_progeny`
    # in the JSON library — explicit calibration set, Stage 3 territory.
    "Ti-44",   # 67.87 + 78.34 + 146.22 — long-lived β+ EC calib source
    "Sc-44",   # 1157.02 — daughter of Ti-44 (secular eq), 99.9% BR
]


# --- Public helpers --------------------------------------------------

def candidates_for_stage(stage: int) -> List[str]:
    """Return the candidate list for a single stage (1, 2, or 3)."""
    if stage == 1:
        return list(ERN_STAGE1)
    if stage == 2:
        return list(TECHNOGENIC_STAGE2)
    if stage == 3:
        return list(EXOTIC_STAGE3)
    raise ValueError(f"Stage must be 1, 2, or 3 — got {stage}")


def candidates_up_to_stage(stage: int) -> List[str]:
    """Cumulative candidate list for stages 1..stage."""
    out: List[str] = []
    for s in range(1, stage + 1):
        out.extend(candidates_for_stage(s))
    # Deduplicate while preserving order
    seen = set()
    return [c for c in out if not (c in seen or seen.add(c))]


def stage_of_nuclide(nuclide: str) -> Optional[int]:
    """Reverse lookup: which stage does this nuclide belong to?"""
    if nuclide in ERN_STAGE1:
        return 1
    if nuclide in TECHNOGENIC_STAGE2:
        return 2
    if nuclide in EXOTIC_STAGE3:
        return 3
    return None


# --- Lsrm 7-line ЕРН/Фон NaI calibration set --------------------------
# Per NOTES_v1.7_methodology.md §9, the canonical NaI calibration anchors
# on a natural-background spectrum. Used by the staged orchestrator for
# (a) stored-cal verification and (b) empirical FWHM(E) fitting when the
# `lsrm_fwhm_polynomial_in_E` stored form is unusable.
ERN_7_REFERENCE_LINES_keV = [
    240.0,    # Pb-212 (238.63) + Pb-214 (241.98) superposition
    351.93,   # Pb-214 (Ra-226 chain) — clean
    511.0,    # Tl-208 (510.77) + annihilation superposition
    1120.29,  # Bi-214 (Ra-226 chain) — clean
    1460.82,  # K-40 — clean
    1764.49,  # Bi-214 (Ra-226 chain) — clean
    2614.51,  # Tl-208 (Th-232 chain) — clean
]


__all__ = [
    "ERN_STAGE1",
    "TECHNOGENIC_STAGE2",
    "EXOTIC_STAGE3",
    "candidates_for_stage",
    "candidates_up_to_stage",
    "stage_of_nuclide",
    "ERN_7_REFERENCE_LINES_keV",
]
