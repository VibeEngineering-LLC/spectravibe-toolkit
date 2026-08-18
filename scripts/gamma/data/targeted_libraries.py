"""
F-278 (v1.17.12, T-026) — Targeted nuclide library presets.

Per ЛСРМ §6 principle of least action: a candidate library should be
**tailored to the analysis task**. Running identification against the
full ENSDF library (~2000 nuclides) introduces noise and increases
false-positive risk, especially on low-resolution NaI Gamma-1S 63×63.

This module exposes 5 curated presets:

  - `OSGI`         — calibration / QA (5 nuclides: Cs-137, Co-60, K-40, Am-241, Eu-152)
  - `OSGI_K`       — extended ОСГИ-К (~20 nuclides including Ra-226 chain)
  - `NPP`          — nuclear power plant operational monitoring (~80 nuclides)
  - `RAO`          — radioactive waste characterization (~100 nuclides)
  - `ENVIRONMENTAL` — environmental & NORM (~30 nuclides)

Each preset returns a list of nuclide names compatible with
`gamma.data.nuclide_library.get_nuclide`. Use them as the
`candidate_nuclides` filter in `identify_nuclides()` or
`analyze_lsrm_spe(library_preset="NPP")`.

Reference: ЛСРМ Algorithmic Foundations 2022 §6, §13 (template method
benefits from smaller, focused libraries).
"""
from __future__ import annotations

from typing import Iterable


# ───── ОСГИ (наиболее частая калибровочная задача) ──────────────────
OSGI_NUCLIDES = (
    "Cs-137",
    "Co-60",
    "K-40",       # часто как стандартный фоновый эталон
    "Am-241",
    "Eu-152",
)


# ───── ОСГИ-К (расширенный) + типичные ЕРН ─────────────────────────
OSGI_K_NUCLIDES = OSGI_NUCLIDES + (
    "Na-22",
    "Ba-133",
    "Cd-109",
    "Mn-54",
    "Zn-65",
    "Cs-134",
    "Ra-226", "Pb-214", "Bi-214",
    "Th-232", "Tl-208", "Ac-228", "Pb-212", "Bi-212",
)


# ───── NPP (АЭС) — operational monitoring + accident-marker isotopes ─
NPP_NUCLIDES = (
    # Fission products (short-lived, ~days-weeks)
    "I-131", "I-132", "I-133", "I-134", "I-135",
    "Te-132",
    "Cs-134", "Cs-136", "Cs-137",
    "Sr-89", "Sr-90", "Y-90",
    "Ba-140", "La-140",
    "Ru-103", "Ru-106", "Rh-106",
    "Mo-99", "Tc-99m",
    "Zr-95", "Nb-95",
    "Ce-141", "Ce-144", "Pr-144",
    "Nd-147",
    "Pm-147",
    # Activation products
    "Mn-54", "Mn-56",
    "Fe-55", "Fe-59",
    "Co-58", "Co-60",
    "Ni-63", "Ni-65",
    "Cu-64",
    "Zn-65",
    "Cr-51",
    "Ag-110m",
    "Sb-122", "Sb-124", "Sb-125",
    # Actinides (accident / fuel-handling)
    "Np-239",
    "Pu-238", "Pu-239", "Pu-240", "Pu-241",
    "Am-241", "Am-243",
    "Cm-242", "Cm-244",
    "U-235", "U-238", "U-234",
    "Th-232", "Th-228", "Th-230", "Th-234",
    "Pa-233", "Pa-234m",
    # Background / NORM
    "K-40",
    "Ra-226", "Pb-210", "Pb-214", "Bi-214",
    "Tl-208", "Ac-228", "Pb-212", "Bi-212",
    "Be-7",
    # Calibration crosschecks
    "Cs-137", "Co-60", "Eu-152",
)


# ───── RAO (РАО — radioactive waste) ────────────────────────────────
RAO_NUCLIDES = NPP_NUCLIDES + (
    # Additional reactor and waste-process species
    "Sr-85", "Sr-86", "Sr-91", "Sr-92",
    "Y-87", "Y-88", "Y-91",
    "Zr-89", "Zr-93", "Zr-97",
    "Tc-99",
    "Ru-105",
    "I-125", "I-129", "I-130",
    "Xe-131m", "Xe-133", "Xe-135",
    "Cs-138",
    "Ba-137m",
    "La-138",
    "Hf-181",
    "Ta-182",
    "W-187",
    "Re-188",
    "Ir-192", "Ir-194",
    "Au-198",
    "Hg-203",
    "Tl-201",
    "Pb-203",
    "Bi-207", "Bi-210",
    "Po-210",
    "At-217",
    "Fr-221",
    "Ac-225",
    # Long-lived for repository characterization
    "Np-237",
    "Pu-242",
    "Cm-243", "Cm-245", "Cm-246",
    "Cf-249", "Cf-252",
)


# ───── ENVIRONMENTAL / NORM ─────────────────────────────────────────
ENVIRONMENTAL_NUCLIDES = (
    # NORM core
    "K-40",
    "Ra-226", "Pb-210", "Po-210",
    "Pb-214", "Bi-214",
    "Ra-228", "Th-232", "Th-228", "Th-230", "Th-234",
    "Tl-208", "Ac-228", "Pb-212", "Bi-212",
    "U-234", "U-235", "U-238",
    "Pa-231",
    # Cosmogenic
    "Be-7",
    "Na-22",     # cosmogenic
    # Anthropogenic fallout
    "Cs-137",
    "Cs-134",
    "Sr-90",
    "I-131",
    "Am-241",
    "Pu-238", "Pu-239", "Pu-240",
    # Medical / industrial release markers
    "I-125", "I-129",
    "Co-60",
)


# Dispatch table
PRESETS = {
    "OSGI":          OSGI_NUCLIDES,
    "OSGI_K":        OSGI_K_NUCLIDES,
    "NPP":           NPP_NUCLIDES,
    "RAO":           RAO_NUCLIDES,
    "ENVIRONMENTAL": ENVIRONMENTAL_NUCLIDES,
}


def get_preset_nuclides(preset_name: str) -> tuple:
    """Вернуть tuple нуклидов для заданного пресета.

    Parameters
    ----------
    preset_name : str
        Один из "OSGI", "OSGI_K", "NPP", "RAO", "ENVIRONMENTAL"
        (case-insensitive).

    Raises
    ------
    KeyError если preset_name неизвестен.
    """
    key = str(preset_name).strip().upper()
    if key not in PRESETS:
        raise KeyError(
            f"Unknown library preset '{preset_name}'; available: "
            f"{sorted(PRESETS.keys())}"
        )
    return tuple(PRESETS[key])


def list_presets() -> tuple:
    """Список имён доступных пресетов с числом нуклидов."""
    return tuple((k, len(v)) for k, v in PRESETS.items())


def combine_presets(*preset_names: str) -> tuple:
    """Объединить несколько пресетов (с дедупликацией, порядок сохраняется
    по первому появлению)."""
    seen = set()
    result = []
    for name in preset_names:
        for n in get_preset_nuclides(name):
            if n not in seen:
                seen.add(n)
                result.append(n)
    return tuple(result)


__all__ = [
    "OSGI_NUCLIDES",
    "OSGI_K_NUCLIDES",
    "NPP_NUCLIDES",
    "RAO_NUCLIDES",
    "ENVIRONMENTAL_NUCLIDES",
    "PRESETS",
    "get_preset_nuclides",
    "list_presets",
    "combine_presets",
]
