"""
Decay-chain line assignment per ENSDF/IAEA references.

Lsrm libraries typically bundle all lines of a decay chain (Ra-226 → ...
→ Pb-210; Th-232 → ... → Pb-208) under the PARENT name (Ra-226 or
Th-232). This is correct for **secular equilibrium** sources (calibration
KIs, sealed bottles aged > 10 half-lives of the longest-lived chain
member). It is INCORRECT for general samples where:

  • Rn-222 (T½=3.82 d, gas) escapes from Ra-226 sources → Pb-214/Bi-214
    are independent of Ra-226 activity.
  • Pb-210 (T½=22.3 y) is dominated by lead-shielding contamination,
    NOT by sample chain (already covered in K-NN limitations).
  • Geochemical fractionation can separate elements at scales much
    longer than chain half-lives (e.g. Ra/Th separation in
    sedimentary processes).

For correct identification and activity calculation, each daughter
nuclide must be tracked SEPARATELY. This module provides the
ENSDF-based line assignment that lets us split combined chain
entries into proper individual nuclides.

Source: ENSDF (Evaluated Nuclear Structure Data File, IAEA) line
energies and intensities. Tolerances reflect typical NaI 50×50
ambiguity at each energy.

If the project's Lsrm library is incomplete (e.g. missing U-238
daughters Th-234, Pa-234m, etc.), users should add entries from
ENSDF/IAEA references; this module's TRUE_ENSDF_OWNERSHIP table
documents the canonical assignments.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


# Canonical ENSDF/IAEA line ownership for natural decay chains.
# Each tuple: (energy_keV, true_emitting_nuclide, intensity_pct, tolerance_keV)
# Tolerance reflects how close lines must be (in keV) to be considered
# the same as one of these reference lines. NaI 50×50 FWHM at e.g.
# 1000 keV is ~60 keV, so a 3 keV tolerance is plenty.
#
# Per ENSDF retrieved 2024-2025:
TRUE_ENSDF_OWNERSHIP_RA_CHAIN = [
    # Energy   Owner       I_pct  tol_keV   ENSDF reference comment
    (186.21,  "Ra-226",     3.59, 1.5),   # Intrinsic Ra-226
    (242.00,  "Pb-214",     7.27, 1.5),   # Rn-222 daughter
    (258.87,  "Pb-214",     0.53, 1.5),
    (274.80,  "Pb-214",     0.42, 1.5),
    (295.22,  "Pb-214",    18.41, 1.5),   # Rn-222 daughter
    (351.93,  "Pb-214",    35.60, 1.5),   # Rn-222 daughter
    (462.00,  "Pb-214",     0.21, 2.0),
    (480.43,  "Pb-214",     0.34, 2.0),
    (487.09,  "Pb-214",     0.43, 2.0),
    (533.66,  "Pb-214",     0.18, 2.0),
    (580.13,  "Pb-214",     0.36, 2.0),
    (609.31,  "Bi-214",    45.49, 1.5),   # Bi-214 main line
    (665.45,  "Bi-214",     1.53, 1.5),
    (703.11,  "Bi-214",     0.46, 1.5),
    (719.86,  "Bi-214",     0.40, 1.5),
    (768.36,  "Bi-214",     4.89, 1.5),
    (785.96,  "Pb-214",     1.06, 1.5),
    (806.17,  "Bi-214",     1.26, 1.5),
    (821.20,  "Bi-214",     0.13, 2.0),
    (904.30,  "Bi-214",     0.71, 2.0),
    (934.06,  "Bi-214",     3.11, 1.5),
    (964.05,  "Bi-214",     0.36, 2.0),
    (1052.0,  "Bi-214",     0.31, 2.0),
    (1070.0,  "Bi-214",     0.27, 2.0),
    (1120.29, "Bi-214",    14.92, 1.5),
    (1133.66, "Bi-214",     0.27, 2.0),
    (1155.19, "Bi-214",     1.69, 1.5),
    (1207.68, "Bi-214",     0.45, 2.0),
    (1238.12, "Bi-214",     5.83, 1.5),
    (1280.96, "Bi-214",     1.43, 1.5),
    (1303.65, "Bi-214",     0.46, 2.0),
    (1377.67, "Bi-214",     3.99, 1.5),
    (1385.31, "Bi-214",     0.78, 2.0),
    (1401.50, "Bi-214",     1.39, 1.5),
    (1407.99, "Bi-214",     2.39, 1.5),
    (1509.21, "Bi-214",     2.13, 1.5),
    (1538.59, "Bi-214",     0.40, 2.0),
    (1583.22, "Bi-214",     0.71, 2.0),
    (1599.31, "Bi-214",     0.33, 2.0),
    (1661.27, "Bi-214",     1.05, 1.5),
    (1683.99, "Bi-214",     0.13, 2.0),
    (1729.60, "Bi-214",     2.88, 1.5),
    (1764.49, "Bi-214",    15.30, 1.5),
    (1838.36, "Bi-214",     0.34, 2.0),
    (1847.43, "Bi-214",     2.02, 1.5),
    (2010.78, "Bi-214",     0.50, 2.0),
    (2021.83, "Bi-214",     0.20, 2.0),
    (2118.51, "Bi-214",     1.16, 1.5),
    (2204.06, "Bi-214",     4.92, 1.5),
    (2293.36, "Bi-214",     0.31, 2.0),
    (2447.70, "Bi-214",     1.55, 1.5),
    (46.54,   "Pb-210",     4.25, 1.5),   # NOT chain-validating; usually shielding contamination
]

TRUE_ENSDF_OWNERSHIP_TH_CHAIN = [
    # Th-232 chain (Th-228 → Ra-224 → ... → Pb-208)
    (238.63,  "Pb-212",    43.60, 1.5),   # Pb-212 main
    (277.36,  "Tl-208",     6.31, 1.5),
    (300.09,  "Pb-212",     3.30, 1.5),
    (328.00,  "Ac-228",     2.95, 1.5),
    (338.32,  "Ac-228",    11.27, 1.5),
    (409.46,  "Ac-228",     1.92, 1.5),
    (463.00,  "Ac-228",     4.40, 1.5),
    (510.77,  "Tl-208",     8.12, 1.5),   # NB: indistinguishable from annihilation on NaI
    (562.50,  "Ac-228",     0.86, 2.0),
    (583.19,  "Tl-208",    30.50, 1.5),   # Tl-208 main
    (727.33,  "Bi-212",     6.74, 1.5),
    (755.30,  "Ac-228",     1.00, 1.5),
    (763.13,  "Tl-208",     1.81, 1.5),
    (785.40,  "Bi-212",     1.10, 1.5),
    (794.95,  "Ac-228",     4.25, 1.5),
    (835.71,  "Ac-228",     1.61, 2.0),
    (860.56,  "Tl-208",     4.50, 1.5),
    (904.20,  "Ac-228",     0.83, 2.0),
    (911.20,  "Ac-228",    25.80, 1.5),
    (964.77,  "Ac-228",     5.11, 1.5),
    (968.97,  "Ac-228",    15.80, 1.5),
    (1078.6,  "Bi-212",     0.55, 2.0),
    (1153.5,  "Ac-228",     0.13, 2.0),
    (1247.1,  "Ac-228",     0.50, 2.0),
    (1459.13, "Ac-228",     0.83, 2.0),
    (1495.91, "Ac-228",     0.86, 2.0),
    (1512.7,  "Bi-212",     0.31, 2.0),
    (1580.5,  "Ac-228",     0.69, 2.0),
    (1588.20, "Ac-228",     3.22, 1.5),
    (1620.50, "Bi-212",     1.51, 1.5),
    (1630.63, "Ac-228",     1.51, 1.5),
    (1638.28, "Ac-228",     0.49, 2.0),
    (1647.4,  "Tl-208",     0.077, 2.0),
    (1666.5,  "Ac-228",     0.18, 2.0),
    (1685.0,  "Ac-228",     0.13, 2.0),
    (1759.5,  "Ac-228",     0.13, 2.0),
    (2614.51, "Tl-208",    35.85, 1.5),   # Tl-208 hallmark (raw 99.75% × 0.3594 Bi-212→Tl-208 α-branching, per F-372.1; was 99.75 — typo, inconsistent with sibling Tl-208 entries already branching-corrected: 583.19=30.50, 510.77=8.12, 860.56=4.50)
]

TRUE_ENSDF_OWNERSHIP_U238_CHAIN = [
    # U-238 chain (above Ra-226)
    (49.55,   "Th-234",     0.069, 2.0),
    (63.30,   "Th-234",     4.84, 1.5),  # Th-234 daughter of U-238
    (92.40,   "Th-234",     5.58, 1.5),  # Th-234 (sum of 92.38 + 92.79)
    (1001.03, "Pa-234m",    0.84, 1.5),  # Pa-234m daughter of Th-234
    # (Ra-226 onwards already in TRUE_ENSDF_OWNERSHIP_RA_CHAIN)
]

TRUE_ENSDF_OWNERSHIP_U235_CHAIN = [
    # U-235 itself
    (143.76,  "U-235",     10.96, 1.5),   # U-235 strong line
    (163.36,  "U-235",      5.08, 1.5),
    (185.72,  "U-235",     57.00, 1.5),   # NB: unresolvable from Ra-226 186.21 on NaI
    (205.31,  "U-235",      5.01, 1.5),
    # Daughters Th-231 (T½ 25.5 h):
    (84.21,   "Th-231",     6.6, 2.0),
    (89.95,   "Th-231",     0.95, 2.0),
]


# Combined lookup table: all known chain lines with their true ENSDF owner.
ALL_CHAIN_LINES = (
    TRUE_ENSDF_OWNERSHIP_RA_CHAIN
    + TRUE_ENSDF_OWNERSHIP_TH_CHAIN
    + TRUE_ENSDF_OWNERSHIP_U238_CHAIN
    + TRUE_ENSDF_OWNERSHIP_U235_CHAIN
)


# Names of chain-parent nuclides whose Lsrm-library entries should be
# split. The Lsrm convention is to use the chain-parent name as a
# combined entry; physically correct identification requires splitting.
LSRM_CHAIN_PARENT_NAMES = {"Ra-226", "Th-232", "U-238", "U-235"}


def reassign_chain_line(
    energy_keV: float,
    parent_name: str,
    *,
    fallback_tolerance_keV: float = 2.0,
) -> Optional[str]:
    """
    For one line of a Lsrm-style chain-parent entry, return the TRUE
    ENSDF emitting nuclide.

    Args:
        energy_keV: line energy
        parent_name: chain-parent name (Ra-226, Th-232, ...) — used as
            fallback if no specific match found
        fallback_tolerance_keV: default tolerance for matching when the
            ENSDF table doesn't specify one

    Returns:
        True ENSDF nuclide name (Pb-214, Bi-214, Tl-208, etc.) or
        None if no match within tolerance — caller should treat such
        lines as "unknown chain member" and ignore them for nuclide-
        specific proportionality, but still use them as anchors for
        energy calibration.
    """
    best_match = None
    best_dE = float("inf")
    for E_ref, owner, I, tol in ALL_CHAIN_LINES:
        dE = abs(energy_keV - E_ref)
        actual_tol = tol if tol > 0 else fallback_tolerance_keV
        if dE <= actual_tol and dE < best_dE:
            best_dE = dE
            best_match = owner
    return best_match


def split_chain_entry(
    lsrm_lines: list,
    parent_name: str,
) -> dict:
    """
    Split a Lsrm-style combined chain entry into separate nuclide entries
    based on ENSDF line assignment.

    Args:
        lsrm_lines: list of (E_keV, I_pct, dI_pct) triples — the
            combined chain entry's `lines` field from
            merge_lsrm_library_into_internal()
        parent_name: chain-parent name (Ra-226, Th-232, ...)

    Returns:
        dict {nuclide_name: [[E_keV, I_pct, dI_pct], ...]} — multiple
        nuclides, each with their own true ENSDF lines.

    Notes:
        Lines that don't match any known ENSDF reference are assigned
        to the parent (fallback). This handles minor library entries
        without losing information.
    """
    result = {}
    unassigned = []
    for triple in lsrm_lines:
        E, I, dI = triple[0], triple[1], triple[2] if len(triple) > 2 else 0.0
        owner = reassign_chain_line(E, parent_name)
        if owner is None:
            unassigned.append([E, I, dI])
            continue
        result.setdefault(owner, []).append([E, I, dI])
    # Assign unmatched lines back to the parent
    if unassigned:
        result.setdefault(parent_name, []).extend(unassigned)
    # Sort each nuclide's lines by energy
    for k in result:
        result[k].sort(key=lambda x: x[0])
    return result


__all__ = [
    "TRUE_ENSDF_OWNERSHIP_RA_CHAIN",
    "TRUE_ENSDF_OWNERSHIP_TH_CHAIN",
    "TRUE_ENSDF_OWNERSHIP_U238_CHAIN",
    "TRUE_ENSDF_OWNERSHIP_U235_CHAIN",
    "ALL_CHAIN_LINES",
    "LSRM_CHAIN_PARENT_NAMES",
    "reassign_chain_line",
    "split_chain_entry",
]
