from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""
Extended multi-source certificate validation (Variant B', v1.7.14;
F-45 bg upgrade in v1.7.23).

Runs the full F-31a + F-31b + F-35 pipeline on every available point-source
5 cm .spe fixture and compares the measured activity against the
certificate value (decay-corrected to the spectrum's measurement date).

Pipeline per fixture (NaI 63x63 mm, Gamma-1S, 5 cm point geometry):
  1. read_spectrum (Lsrm SpectraLine .spe; embeds peak table + FWHM cal)
  2. subtract bg_2016_empty_shield_point5cm.spe (F-43 averaged, F-44
     cumulative_last semantic; 15-hour live time, matrix-matched
     "empty shield closed lid" geometry for point-5cm samples)
  3. build FWHM(channel) from spec FWHM polynomial in E
  4. mariscotti_search peaks
  5. identify_nuclides (target nuclide only -- single-source ground truth)
  6. disambiguate_identifications
  7. apply_multiplet_deconvolution
  8. compute_tcs_corrections (analytic Knoll/Gilmore TCS, NaI P/T)
  9. compute_activity with tcs_method_scale defaults
       (F-35: Lsrm peak-table area gets scale=0; others get full TCS)
 10. decay-correct from certificate date to measurement date
 11. deviation %  =  (A_measured - A_cert_decayed) / A_cert_decayed * 100

Library-gap nuclides (Th-228, Bi-207, Cd-109, Y-88) are reported as
"no library record" rather than skipped silently.

Output:
  - deviation matrix printed to stdout (Russian text)
  - cert_validation_matrix.csv saved next to this script

Methodology references:
  - Gilmore & Joss, "Practical Gamma-ray Spectrometry" 3rd Ed., chap. 8
  - Knoll, "Radiation Detection and Measurement" 4th Ed., chap. 12.D, 17.6
  - Lsrm Algorithmic Foundations 2022 sec. 8.4, sec. 10
"""


import math
import os
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable, Optional

sys.path.insert(0, str(Path(__file__).parent))

import numpy as np

from gamma.io.readers import read_spectrum
from gamma.io.lsrm_src import read_certificate_file
from gamma.calibration.efficiency import fit_efficiency_from_efr_file
from gamma.calibration.bg_subtract_dual_mode import (
    subtract_background, apply_subtraction_to_spectrum,
)
from gamma.peaks.search import mariscotti_search
from gamma.peaks.deconvolve import apply_multiplet_deconvolution
from gamma.identification.identify import identify_nuclides
from gamma.identification.disambiguate import disambiguate_identifications
from gamma.physics.cascade_summing import (
    compute_tcs_corrections,
    peak_to_total_NaI_for_geometry,
)
from gamma.physics.self_attenuation import (
    correction_factor as k20_correction_factor,
    weighted_mean_correction as k20_weighted_correction,
    REF_GEOMETRY as K20_REF_GEOMETRY,
    OISN_16_COMPOSITION,
)
from gamma.activity import compute_activity
from gamma.data.nuclide_library import get_nuclide, load_lsrm_chain_libs

# F-39 (v1.7.17): supplement the built-in 27-nuclide library with the
# bundled Lsrm chain libraries (NaI-Etl+Esc.lib + ОСГИ.lib). This brings
# 47+ nuclides incl. the full Th-232 daughter chain (Tl-208, Pb-212,
# Bi-212, Ac-228, Ra-224), which is required to measure the Th-228
# cert source via its detectable secular-equilibrium daughters.
load_lsrm_chain_libs()


# ---------------------------------------------------------------------------
# Project paths
# ---------------------------------------------------------------------------

ROOT = Path(__file__).parent.parent

# F-83 (v1.12.0): paths moved into detectors/Gamma-1S/. Resolver lives in
# gamma.detectors.gamma1s — see detectors/Gamma-1S/README.md.
from gamma.detectors.gamma1s import (
    DEFAULT_REFERENCE_DIR as REF_DIR,
    DEFAULT_EFFICIENCY_DIR as EFF_DIR,
    CERTIFICATES_DIR as CERT_DIR,
    AVERAGED_BACKGROUNDS_DIR as BG_DIR,
)

# F-46 (v1.7.24): multi-geometry expansion. Resolvers map each
# fixture's `geometry` token to (a) the appropriate efficiency
# curve and (b) the canonical averaged background per F-44 pairing
# rules. The current set of supported geometries:
#
#   Точечная-5см    -- point source 5 cm:   efr Точечная-5см  / bg empty_shield_point5cm
#   Точечная-25см   -- point source 25 cm:  efr Точечная-25см / bg open_lid_point25cm
#   Дента-120мл     -- Denta-120 cylinder:  efr Дента         / bg empty_shield_point5cm
#   Петри-60мл      -- Petri-60 dish:       efr Петри         / bg empty_shield_point5cm
#   Маринелли       -- Marinelli vessel:    efr Маринелли     / bg marinelli_water_marinelli
#
# Pairing rationale (NOTES_v1.7_methodology.md §v1.7.22):
#   - Marinelli sample (matrix-attenuating vessel) → Marinelli+water bg
#     (matrix-matched attenuation envelope so bg subtraction cancels
#     uniform-room γ-flux correctly).
#   - Closed-geometry samples in the shield (Дента, Петри, Точечная-5см)
#     → "empty shield, closed lid" bg (no source-induced attenuation
#     beyond the shield itself).
#   - Точечная-25см at 25 cm distance is above shield top, lid open
#     during measurement → "open lid" bg.
#
# F-45 fix (v1.7.23): previously hardcoded single Marinelli+water bg
# for all point-5cm fixtures (wrong geometry); F-45 fixed Точ-5см alone.
# F-46 generalizes to all 5 geometries.

EFF_PATHS = {
    "Точечная-5см":   EFF_DIR / "УДС-ГЦ-63х63-USB__SN-01_-_Точечная-5см.efr",
    "Точечная-25см":  EFF_DIR / "УДС-ГЦ-63х63-USB__SN-01_-_Точечная-25см.efr",
    "Дента-120мл":    EFF_DIR / "УДС-ГЦ-63х63-USB__SN-01_-_Дента.efr",
    "Петри-60мл":     EFF_DIR / "УДС-ГЦ-63х63-USB__SN-01_-_Петри.efr",
    "Маринелли":      EFF_DIR / "УДС-ГЦ-63х63-USB__SN-01_-_Маринелли.efr",
}

BG_PATHS = {
    "Точечная-5см":   BG_DIR / "bg_2016_empty_shield_point5cm.spe",
    "Точечная-25см":  BG_DIR / "bg_2016_open_lid_point25cm.spe",
    "Дента-120мл":    BG_DIR / "bg_2016_empty_shield_point5cm.spe",
    "Петри-60мл":     BG_DIR / "bg_2016_empty_shield_point5cm.spe",
    "Маринелли":      BG_DIR / "bg_2016_marinelli_water_marinelli.spe",
}

# Per-geometry certificate paths. F-46 multi-cert: point sources use
# the АСПЕКТ_ОСГИ_2024 cert; volume sources (Marinelli/Дента/Петри)
# use their dedicated 2017 cert files.
CERT_PATHS = {
    "Точечная-5см":   CERT_DIR / "АСПЕКТ_ОСГИ_2024.src",
    "Точечная-25см":  CERT_DIR / "АСПЕКТ_ОСГИ_2024.src",
    "Дента-120мл":    CERT_DIR / "Эталон_Дента120мл__Аспект2017_.src",
    "Петри-60мл":     CERT_DIR / "Эталон_Петри-60__Аспект_2017_.src",
    "Маринелли":      CERT_DIR / "Эталон_Маринелли__Аспект_2017_.src",
}

# Legacy single-bg path for diagnostic comparison (pre-F-45 baseline).
# BG_PATH_LEGACY = REF_DIR / "Фон_закр_кр_вода_01.spe"

# Back-compat module-level constants (existing tests in
# test_chain_proxy.py grep validate_certs.py for `BG_PATH` / `CERT_PATH`
# / `EFF_5CM` symbols; keep them aliased to the Точ-5см defaults so the
# string-search tests stay green).
BG_PATH = BG_PATHS["Точечная-5см"]
EFF_5CM = EFF_PATHS["Точечная-5см"]
CERT_PATH = CERT_PATHS["Точечная-5см"]


# ---------------------------------------------------------------------------
# Fixture <-> certificate map
# ---------------------------------------------------------------------------
#
# Each row maps a 5 cm .spe fixture to its certificate source name and
# the target nuclide whose activity the certificate reports. Source-name
# strings are matched against the parsed Certificate via
# `find_source_fuzzy` (tokenized substring match), so minor punctuation
# variants in the .src file work without exact equality.

@dataclass
class CertFixture:
    nuclide: str               # nuclide whose ACTIVITY we measure ("Co-60")
    spe_filename: str          # filename relative to REF_DIR (use forward-slash
                                # subpath for non-root geometries, e.g.
                                # "Точечная-25см/Cs-137 №SRC-02_...spe")
    cert_source_hint: str      # substring used to find the .src source
    # F-39: chain-proxy fields. When set, the cert reports activity of
    # `cert_nuclide` (e.g. "Th-228" parent) but only the daughter
    # `nuclide` (e.g. "Pb-212") has detectable γ-lines. In secular
    # equilibrium A_daughter = `chain_branching` * A_parent, where the
    # branching factor is 1.0 for direct-chain daughters (Pb-212 from
    # Th-228 via Ra-224 → ... → Pb-212, all 1:1) and < 1 for branched
    # paths (Tl-208 gets only 36% of Bi-212 decays).
    cert_nuclide: Optional[str] = None         # parent name in cert
    chain_branching: float = 1.0               # A_daughter / A_parent
    # F-46: geometry token — keys into EFF_PATHS/BG_PATHS/CERT_PATHS so
    # each fixture resolves to its own efficiency curve, background, and
    # certificate file. Default Точечная-5см preserves all v1.7.23 rows
    # without modification.
    geometry: str = "Точечная-5см"
    # F-51 (v1.11.0) — K-22: Th-232 chain-equilibrium correction.
    #
    # When the cert reports the PARENT activity of a Th-232 source and
    # the chain was reset to zero at cert reference date (e.g. chemical
    # separation during source preparation), the measurable daughter
    # γ-emission rate (e.g. Tl-208 at 2614 keV) is suppressed by the
    # Ra-228 bottleneck (T½ = 5.75 y) until the chain re-grows toward
    # secular equilibrium. The harness applies the in-growth factor
    #
    #     A_daughter@meas = A_parent@meas × (1 − 2^(−Δt / T_Ra-228))
    #
    # where Δt = meas_dt − cert_ref_dt. Default `True` means the cert
    # already reports the equilibrium daughter rate (no correction
    # needed) — this preserves v1.10.1 behavior for all 40 fixtures.
    # Setting `False` is supported for the 3 heavy Marinelli/Дента/
    # Петри Th-232 fixtures (cert 420-7-17 / 2007-09 ref) where the
    # 17-year-old chain has grown to ≈ 87 % of equilibrium and the
    # empirical Δ was −13 to −16 % before this correction (see K-22 in
    # KNOWN_AND_FIXED_ISSUES.md).
    chain_at_cert_equilibrium: bool = True
    # Half-life of the chain bottleneck (seconds). Required when
    # chain_at_cert_equilibrium=False; ignored otherwise. For Th-232
    # the bottleneck is Ra-228 (T½ = 5.75 y → 1.81417e8 s; see
    # RA228_T_HALF_S constant below).
    chain_bottleneck_T_half_s: Optional[float] = None


# F-51 (v1.11.0) — K-22 closure: Ra-228 half-life as the Th-232 → ... →
# Tl-208 chain bottleneck. T½ = 5.75 years (commonly cited 5.75 ± 0.04 y;
# ENSDF: 5.75 y). One Julian year = 365.25 d × 86400 s = 3.15576e7 s.
RA228_T_HALF_S = 5.75 * 365.25 * 86400.0   # ≈ 1.81417e8 s


FIXTURES = [
    CertFixture("Cs-137", "Cs-137__163_2017.spe",
                "Cs-137 SRC-02"),
    CertFixture("Co-60",  "Co-60__043_02_2019_Точечная-5см_5cm.spe",
                "Co-60 043 02.2019"),
    CertFixture("Na-22",  "Na-22__01_22.spe",
                "Na-22 01.22"),
    CertFixture("Eu-152", "Eu-152__04_21_Точечная-5см_5cm.spe",
                "Eu-152 04.21"),
    CertFixture("Ba-133", "Ba-133__SRC-05_Точечная-5см_5cm.spe",
                "SRC-05"),
    CertFixture("Am-241", "Am-241_045_02_2019_Точечная-5см_5cm.spe",
                "Am-241 045 02.2019"),
    CertFixture("Zn-65",  "Zn-65__342_2019_Точечная-5см_5cm.spe",
                "Zn-65 342.2019"),
    # Bi-207 and Cd-109 sit inside the "ОСГИ 5431" multi-nuclide
    # certificate (Ti-44 + Cd-109 + Bi-207 + Th-228 at 25.05.2017);
    # `find_source_fuzzy('5431')` resolves to that compound source and
    # `get_activity(nuclide)` then picks the right sub-source.
    CertFixture("Y-88",   "Y-88__260_2023_Точечная-5см_5cm.spe",
                "Y-88 SRC-06"),
    CertFixture("Bi-207", "Bi-207__176_04_2017_Точечная-5см_5cm.spe",
                "5431"),
    CertFixture("Cd-109", "Cd-109__175_04_2017_Точечная-5см_5cm.spe",
                "5431"),
    # Th-228 cert source: the certificate reports the parent Th-228
    # activity (129000 Bq), but only daughters (Tl-208, Pb-212, Bi-212,
    # Ac-228) emit detectable γ-lines on NaI. F-39 closes this gap via
    # chain-proxy: we measure Pb-212 (direct daughter with 1:1
    # branching: Th-228 → Ra-224 → Rn-220 → Po-216 → Pb-212, all
    # alpha/beta with no branching), then compare A_Pb-212 directly
    # to the cert's A_Th-228 entry. Pb-212 has the cleanest detectable
    # line (238.63 keV at I=43.6%, well separated).
    CertFixture("Pb-212", "Th-228__264_2023_Точечная-5см_5cm.spe",
                "Th-228 SRC-01",
                cert_nuclide="Th-228", chain_branching=1.0),
    # F-41 (v1.7.19): Tl-208 as alternative chain proxy for Th-228.
    # In the Th-232/228 chain, Bi-212 (which is in 1:1 secular
    # equilibrium with Th-228 through Pb-212) decays by two competing
    # modes:
    #   • β⁻ to Po-212 (branching 64.06%)  → quickly α to Pb-208
    #   • α to Tl-208  (branching 35.94%) → β⁻ to Pb-208 with γ
    # Tl-208 therefore sees only 0.3594 of the parent Th-228 disintegration
    # rate. However, the bundled Lsrm chain library (NaI-Etl+Esc.lib via
    # load_lsrm_chain_libs) reports Tl-208 γ-line intensities ALREADY
    # pre-scaled by this β-branching factor: e.g. 583.19 keV at I=30.6%
    # (= 0.3594 × 84.5% raw Tl-208 branching from ENSDF). So when
    # compute_activity inverts the lib intensities it recovers the
    # PARENT Th-228 activity directly, and chain_branching stays 1.0.
    # This gives an independent cross-validation path for Th-228 using
    # different γ-lines (583/860/2614 keV) than the Pb-212 row above
    # (which uses only 238.63 keV).
    CertFixture("Tl-208", "Th-228__264_2023_Точечная-5см_5cm.spe",
                "Th-228 SRC-01",
                cert_nuclide="Th-228", chain_branching=1.0),

    # -----------------------------------------------------------------
    # F-46 (v1.7.24) — Точечная-25см expansion.
    # -----------------------------------------------------------------
    # Same АСПЕКТ_ОСГИ_2024 cert as 5cm rows, but different physical
    # source IDs measured at the 25 cm distance (open-lid bg per F-44).
    # Three single-nuclide direct rows + one Th-228 chain proxy via
    # Tl-208 (mirroring the F-41 design — Pb-212 single-line proxy not
    # included here to avoid repeating the F-41 low-E 238 keV
    # self-absorption finding at a second geometry).
    CertFixture("Cs-137", "Точечная-25см/Cs-137 №SRC-02_Точечная-25см_25cm.spe",
                "Cs-137 SRC-02",
                geometry="Точечная-25см"),
    CertFixture("Na-22",  "Точечная-25см/Na-22 #01.22_Точечная-25см_25cm.spe",
                "Na-22 01.22",
                geometry="Точечная-25см"),
    CertFixture("Y-88",   "Точечная-25см/Y-88 №SRC-06_Точечная-25см_25cm.spe",
                "Y-88 SRC-06",
                geometry="Точечная-25см"),
    # Th-228 №SRC-07 (point at 25 cm) — A_cert = 100 000 Bq at
    # 2021-04-26. Distinct physical source from №SRC-01 used in the
    # 5 cm rows. We measure Tl-208 as the chain proxy (cleaner high-E
    # lines on NaI than Pb-212's single 238 keV line; F-41 showed
    # Δ=−0.07% at 5cm — geometry-cross-validates the methodology).
    CertFixture("Tl-208", "Точечная-25см/Th-228 №309_Точечная-25см_25cm.spe",
                "Th-228 SRC-07",
                cert_nuclide="Th-228", chain_branching=1.0,
                geometry="Точечная-25см"),

    # -----------------------------------------------------------------
    # F-46b (v1.7.25) — Marinelli volume-source expansion.
    # -----------------------------------------------------------------
    # All Marinelli sources are calibration soils/concentrates of
    # known specific activity (Bq/kg) prepared by АСПЕКТ in two
    # variants:
    #   • р06 (cert ref 2002-05-24) — smaller containers (570-660 g)
    #     with higher specific activity. Sub-IDs end in 14/16/18/20.
    #   • р16 (cert ref 2007-09-17) — larger containers (1550-1670 g)
    #     with slightly lower specific activity. Sub-IDs end in
    #     15/17/19/21.
    # Cert reports A in Bq/kg; harness multiplies by sub_source.mass_g
    # to compare against absolute A_meas. Pairing: Marinelli +
    # water-filled bg (matrix-matched per F-44 rules). Decay
    # correction uses parent T½ — negligible for K-40 (1.25×10⁹ y),
    # Ra-226 (1600 y) and Th-232 (1.4×10¹⁰ y); ~50% drop for Cs-137
    # over 17-22 y.
    CertFixture("Cs-137", "Cs137_420-7-14_Маринелли_0cm.spe",
                "420-7_р06", geometry="Маринелли"),
    CertFixture("Cs-137", "Cs137_420-7-15_Маринелли_0cm.spe",
                "420/7_р16", geometry="Маринелли"),
    CertFixture("K-40",   "K40_420-7-20_Маринелли_0cm.spe",
                "420-7_р06", geometry="Маринелли"),
    CertFixture("K-40",   "K40_420-7-21_Маринелли_0cm.spe",
                "420/7_р16", geometry="Маринелли"),
    # Ra-226 chain proxy via Bi-214 (9 lines incl. 609/1120/1764 keV).
    # In sealed Marinelli the Rn-222 buffer is retained → Pb-214 and
    # Bi-214 reach secular equilibrium with Ra-226 in ~25-30 days
    # post-sealing (Rn-222 T½=3.8 d → ≥6 half-lives). 420-series
    # sources are ≥20 years old, deep equilibrium. Lib Bi-214
    # intensities are direct ENSDF (no β-branching pre-scaling),
    # so chain_branching=1.0 yields A_Bi-214 ≈ A_Ra-226.
    CertFixture("Bi-214", "Ra226_420-7-18_Маринелли_0cm.spe",
                "420-7_р06",
                cert_nuclide="Ra-226", chain_branching=1.0,
                geometry="Маринелли"),
    CertFixture("Bi-214", "Ra226_420-7-19_Маринелли_0cm.spe",
                "420/7_р16",
                cert_nuclide="Ra-226", chain_branching=1.0,
                geometry="Маринелли"),
    # Th-232 chain proxy via Tl-208 (5 lines incl. 583/2614 keV). In
    # secular equilibrium the entire Th-232 → Ra-228 → Ac-228 → Th-228
    # → … → Bi-212 → Tl-208 chain shares the same disintegration
    # rate; Tl-208 receives only 0.3594 of Bi-212's decay (Bi-212
    # α-branches to Tl-208 35.94% of the time, β-decays to Po-212
    # 64.06%). The Lsrm chain library already pre-scales Tl-208
    # γ-intensities by 0.3594, so compute_activity inverting the
    # library intensities recovers parent A_Th-232 directly with
    # chain_branching=1.0 (same trick as F-41 for Th-228). Ra-228
    # bottleneck (T½=5.75 y) reaches 93% equilibrium in 22 years
    # since 2002, full equilibrium in 17 years since 2007 — well
    # within the 22-year and 17-year ages of these sources.
    # Marinelli Th-232 sources: the р06 variant has no
    # `Th232_420-7-16_Маринелли_0cm.spe` measurement on file (cert
    # entry exists but no spectrum). Use Th-232 №420-17031 (cert
    # ref 2017-06-05, 640 g @ 860 Bq/kg) as the second Marinelli
    # Th-232 fixture — different physical source from 420-7-17 but
    # same Tl-208 chain-proxy methodology applies.
    # F-51 (v1.11.0) K-22: 420-7-17 cert ref 2007-09 → 17-year-old chain
    # at meas (2024) ≈ 87 % of full equilibrium. Marking
    # chain_at_cert_equilibrium=False makes the harness scale cert@meas
    # by the in-growth factor (1 − 2^(−Δt/T_Ra-228)).
    CertFixture("Tl-208", "Th232_420-7-17_Маринелли_0cm.spe",
                "420/7_р16",
                cert_nuclide="Th-232", chain_branching=1.0,
                geometry="Маринелли",
                chain_at_cert_equilibrium=False,
                chain_bottleneck_T_half_s=RA228_T_HALF_S),
    CertFixture("Tl-208", "Th-232_420-17031_Маринелли_0cm.spe",
                "420-17031",
                cert_nuclide="Th-232", chain_branching=1.0,
                geometry="Маринелли"),

    # -----------------------------------------------------------------
    # F-46c (v1.7.25) — Дента-120мл cylinder expansion.
    # -----------------------------------------------------------------
    # Same 4 nuclides × 2 source variants as Marinelli, but in 120 ml
    # cylindrical container (Дента geometry — 6.7 cm Ø × 4.0 cm
    # height). Masses: р06 variant ≈66-79 g, р16 ≈186-200 g. Pairing:
    # empty shield closed lid bg (per F-44; matches Точечная-5см
    # geometry of bg_2016_empty_shield_point5cm).
    CertFixture("Cs-137", "Дента-120мл/Cs137_420-7-14_Дента-120мл_0cm.spe",
                "420-7_р06", geometry="Дента-120мл"),
    CertFixture("Cs-137", "Дента-120мл/Cs137_420-7-15_Дента-120мл_0cm.spe",
                "420-7_р16", geometry="Дента-120мл"),
    CertFixture("K-40",   "Дента-120мл/K40_420-7-20_Дента-120мл_0cm.spe",
                "420-7_р06", geometry="Дента-120мл"),
    CertFixture("K-40",   "Дента-120мл/K40_420-7-21_Дента-120мл_0cm.spe",
                "420-7_р16", geometry="Дента-120мл"),
    CertFixture("Bi-214", "Дента-120мл/Ra226_420-7-18_Дента-120мл_0cm.spe",
                "420-7_р06",
                cert_nuclide="Ra-226", chain_branching=1.0,
                geometry="Дента-120мл"),
    CertFixture("Bi-214", "Дента-120мл/Ra226_420-7-19_Дента-120мл_0cm.spe",
                "420-7_р16",
                cert_nuclide="Ra-226", chain_branching=1.0,
                geometry="Дента-120мл"),
    # F-51 (v1.11.0) K-22: see Marinelli 420-7-17 note above.
    CertFixture("Tl-208", "Дента-120мл/Th232_420-7-17_Дента-120мл_0cm.spe",
                "420-7_р16",
                cert_nuclide="Th-232", chain_branching=1.0,
                geometry="Дента-120мл",
                chain_at_cert_equilibrium=False,
                chain_bottleneck_T_half_s=RA228_T_HALF_S),
    CertFixture("Tl-208", "Дента-120мл/Th-232_420-17031_Дента-120мл_0cm.spe",
                "420-17031",   # 2017 cert, single Th-232 entry
                cert_nuclide="Th-232", chain_branching=1.0,
                geometry="Дента-120мл"),

    # -----------------------------------------------------------------
    # F-46c (v1.7.25) — Петри-60мл shallow dish expansion.
    # -----------------------------------------------------------------
    # Same 4 nuclides × 2 source variants in 60 ml Petri dish (close-
    # geometry thin-source configuration). Masses: р06 variant
    # ≈33-40 g, р16 ≈93-100 g. Pairing: empty shield closed lid bg
    # (per F-44).
    CertFixture("Cs-137", "Петри-60мл/Cs137_420-7-14_Петри-60мл_0cm.spe",
                "420-7_р06", geometry="Петри-60мл"),
    CertFixture("Cs-137", "Петри-60мл/Cs137_420-7-15_Петри-60мл_0cm.spe",
                "420/7_р16", geometry="Петри-60мл"),
    CertFixture("K-40",   "Петри-60мл/K40_420-7-20_Петри-60мл_0cm.spe",
                "420-7_р06", geometry="Петри-60мл"),
    CertFixture("K-40",   "Петри-60мл/K40_420-7-21_Петри-60мл_0cm.spe",
                "420/7_р16", geometry="Петри-60мл"),
    CertFixture("Bi-214", "Петри-60мл/Ra226_420-7-18_Петри-60мл_0cm.spe",
                "420-7_р06",
                cert_nuclide="Ra-226", chain_branching=1.0,
                geometry="Петри-60мл"),
    CertFixture("Bi-214", "Петри-60мл/Ra226_420-7-19_Петри-60мл_0cm.spe",
                "420/7_р16",
                cert_nuclide="Ra-226", chain_branching=1.0,
                geometry="Петри-60мл"),
    # F-51 (v1.11.0) K-22: see Marinelli 420-7-17 note above.
    CertFixture("Tl-208", "Петри-60мл/Th232_420-7-17_Петри-60мл_0cm.spe",
                "420/7_р16",
                cert_nuclide="Th-232", chain_branching=1.0,
                geometry="Петри-60мл",
                chain_at_cert_equilibrium=False,
                chain_bottleneck_T_half_s=RA228_T_HALF_S),
    CertFixture("Tl-208", "Петри-60мл/Th-232_420-17031_Петри-60мл_0cm.spe",
                "420_17031",   # 2017 cert single Th-232 entry (Петри cert
                               # uses underscore between "420" and "17031";
                               # Дента/Marinelli certs use hyphen — both
                               # supported by find_source_fuzzy tokenization)
                cert_nuclide="Th-232", chain_branching=1.0,
                geometry="Петри-60мл"),
]


# ---------------------------------------------------------------------------
# FWHM(channel) provider for Lsrm .spe ('lsrm_fwhm_polynomial_in_E' model)
# ---------------------------------------------------------------------------
#
# The .spe FWHM polynomial is actually FWHM_keV(sqrt(E_keV)) -- NOT
# FWHM(E) as a naive reading of the field name suggests. Verified on
# Co-60 5cm fixture: at E=662 the formula returns 45.6 keV (NaI ~7%, so
# ~46 keV) and at E=1332 it returns 71.6 keV (NaI ~5.4%, so ~72 keV).
# Both match independent NaI 63x63 typical FWHM specs.
#
# Build the (channel -> FWHM in channels) callable as:
#   E(N)          = sum a_k * N^k          (Horner on energy_cal)
#   FWHM_keV(E)   = sum c_k * sqrt(E)^k    (Horner on FWHM coefs, in sqrt(E))
#   dE/dN at N    = sum k * a_k * N^(k-1)  (energy_cal derivative)
#   FWHM_chans    = FWHM_keV / |dE/dN|

def make_lsrm_fwhm_provider(spec) -> Callable[[float], float]:
    fwhm_coefs = (spec.stored_fwhm_calibration.coefficients
                  if spec.stored_fwhm_calibration else (3.0,))
    e_cal = spec.energy_cal
    e_deriv = tuple(k * e_cal[k] for k in range(1, len(e_cal)))

    def fwhm_at(ch: float) -> float:
        N = float(ch)
        E = 0.0
        for k in range(len(e_cal) - 1, -1, -1):
            E = E * N + e_cal[k]
        if E <= 0:
            return 5.0
        sqE = math.sqrt(E)
        fw_keV = 0.0
        for k in range(len(fwhm_coefs) - 1, -1, -1):
            fw_keV = fw_keV * sqE + fwhm_coefs[k]
        if fw_keV <= 0:
            return 5.0
        dEdN = 0.0
        for k in range(len(e_deriv) - 1, -1, -1):
            dEdN = dEdN * N + e_deriv[k]
        if dEdN <= 0:
            return 5.0
        return max(1.5, fw_keV / dEdN)
    return fwhm_at


# ---------------------------------------------------------------------------
# Per-fixture pipeline
# ---------------------------------------------------------------------------

@dataclass
class RunResult:
    nuclide: str
    spe_filename: str
    measurement_datetime: Optional[datetime]
    cert_reference_datetime: Optional[datetime]
    cert_A_Bq: Optional[float]
    cert_A_decayed_Bq: Optional[float]
    measured_A_Bq: Optional[float]
    measured_sigma_Bq: Optional[float]
    deviation_pct: Optional[float]
    n_lines_used: int
    note: str
    # F-46b (v1.7.25): track geometry on the result so the per-geometry
    # summary table at the end can aggregate without re-matching against
    # FIXTURES. Default "Точечная-5см" preserves backward compat when
    # an exception path constructs a RunResult without geometry.
    geometry: str = "Точечная-5см"


def _decay_to_meas(A_cert: float, T_half_s: Optional[float],
                   t_cert: datetime, t_meas: datetime) -> float:
    if not T_half_s or t_cert is None or t_meas is None:
        return A_cert
    dt_s = (t_meas - t_cert).total_seconds()
    return A_cert * math.exp(-math.log(2.0) * dt_s / T_half_s)


def run_one(fx: CertFixture, eff_curve, bg_spec, cert) -> RunResult:
    spe_path = REF_DIR / fx.spe_filename
    note_parts = []

    if not spe_path.exists():
        return RunResult(fx.nuclide, fx.spe_filename, None, None, None,
                         None, None, None, None, 0, "fixture missing",
                         geometry=fx.geometry)

    # Library check: we identify fx.nuclide (which may be a daughter)
    lib = get_nuclide(fx.nuclide)
    if not lib:
        src = cert.find_source_fuzzy(fx.cert_source_hint)
        cert_lookup_name = fx.cert_nuclide or fx.nuclide
        cert_A = src.get_activity(cert_lookup_name) if src else None
        return RunResult(
            fx.nuclide, fx.spe_filename, None,
            src.reference_datetime if src else None,
            cert_A.A_Bq if cert_A else None,
            None, None, None, None, 0,
            "library gap: no record for this nuclide",
            geometry=fx.geometry,
        )

    # Read spectrum
    spec = read_spectrum(str(spe_path))
    meas_dt = spec.start_datetime
    live = float(spec.live_time)

    # Certificate lookup. cert_nuclide_name may differ from fx.nuclide
    # when we're identifying a daughter (Pb-212) to infer the parent
    # (Th-228) activity from a chain-equilibrium source.
    cert_nuclide_name = fx.cert_nuclide or fx.nuclide
    src = cert.find_source_fuzzy(fx.cert_source_hint)
    if src is None:
        return RunResult(fx.nuclide, fx.spe_filename, meas_dt, None, None,
                         None, None, None, None, 0,
                         f"cert source not found for hint {fx.cert_source_hint!r}",
                         geometry=fx.geometry)

    # F-46b (v1.7.25): for volume sources (Marinelli/Дента/Петри) the
    # cert reports activity per unit mass (Bq/kg) and we must multiply
    # by the sub_source's `mass_g` to recover absolute Bq. Point
    # sources continue to report absolute Bq directly. We need the
    # sub_source object (not just CertificateActivity) to access
    # mass_g. Walk sub_sources to find the one carrying the target
    # cert_nuclide_name.
    target_sub = None
    cert_act = None
    for ss in src.sub_sources:
        for act in ss.activities:
            if act.nuclide == cert_nuclide_name:
                target_sub = ss
                cert_act = act
                break
        if cert_act is not None:
            break
    if cert_act is None:
        # Fall back to legacy compound-source path (e.g., ОСГИ 5431
        # multi-nuclide cert where get_activity walks all sub_sources
        # by nuclide name). This keeps Bi-207 and Cd-109 5cm fixtures
        # working unchanged.
        cert_act = src.get_activity(cert_nuclide_name)
        if cert_act is None:
            return RunResult(fx.nuclide, fx.spe_filename, meas_dt,
                             src.reference_datetime, None, None, None, None, None, 0,
                             f"cert source has no entry for {cert_nuclide_name!r}",
                             geometry=fx.geometry)
        target_sub = None  # mass unavailable for legacy fallback

    # Apply Bq/kg → Bq conversion if cert is per-mass. Point sources
    # use unit "Bq" or "" (absolute); volume sources use "Bq/kg".
    unit = (getattr(cert_act, "unit", "") or "").strip().lower()
    if unit.endswith("bq/kg"):
        if target_sub is None or target_sub.mass_g is None:
            return RunResult(fx.nuclide, fx.spe_filename, meas_dt,
                             src.reference_datetime, cert_act.A_Bq, None, None,
                             None, None, 0,
                             f"cert reports {cert_act.A_Bq:.2e} {cert_act.unit} "
                             "but sub_source.mass_g unavailable — cannot convert",
                             geometry=fx.geometry)
        # convert Bq/kg × g / 1000 g/kg = Bq
        A_cert_absolute = cert_act.A_Bq * target_sub.mass_g / 1000.0
        cert_A_reported_for_csv = A_cert_absolute  # absolute, post-mass-multiplication
    elif unit in ("bq", ""):
        A_cert_absolute = cert_act.A_Bq
        cert_A_reported_for_csv = cert_act.A_Bq
    else:
        return RunResult(fx.nuclide, fx.spe_filename, meas_dt,
                         src.reference_datetime, cert_act.A_Bq, None, None,
                         None, None, 0,
                         f"unsupported cert unit {cert_act.unit!r}",
                         geometry=fx.geometry)

    # For chain-proxy fixtures the cert nuclide is the parent. Decay-
    # correct using the PARENT's half-life (the daughter is in secular
    # equilibrium so its observed activity tracks the parent's decay).
    if fx.cert_nuclide:
        parent_lib = get_nuclide(fx.cert_nuclide) or {}
        T_half_s = float(parent_lib.get("T_half_s") or 0) or None
    else:
        T_half_s = float(lib.get("T_half_s") or 0) or None
    A_cert_at_meas = _decay_to_meas(
        A_cert_absolute, T_half_s, src.reference_datetime, meas_dt,
    )
    # Apply chain branching: expected daughter activity = parent * branching
    A_cert_at_meas *= fx.chain_branching

    # F-51 (v1.11.0) — K-22 chain-equilibrium correction. When the cert
    # reports parent activity and the daughter chain was reset to zero
    # at cert ref date (e.g. chemical separation during preparation),
    # the daughter γ-emission rate at meas time is suppressed by the
    # bottleneck-isotope in-growth factor (1 − 2^(−Δt / T_bottleneck)).
    # For Th-232 sources the bottleneck is Ra-228 (T½ = 5.75 y); 17 y of
    # in-growth → ≈ 0.870 of equilibrium → 13 % daughter under-emission
    # vs naive cert@meas. See K-22 in KNOWN_AND_FIXED_ISSUES.md.
    if (not fx.chain_at_cert_equilibrium
            and fx.chain_bottleneck_T_half_s
            and src.reference_datetime is not None
            and meas_dt is not None):
        dt_s = (meas_dt - src.reference_datetime).total_seconds()
        T_b = float(fx.chain_bottleneck_T_half_s)
        eq_factor = 1.0 - math.exp(-math.log(2.0) * dt_s / T_b)
        A_cert_at_meas *= eq_factor
        note_parts.append(f"K22 eq={eq_factor:.3f}")

    # F-47b (v1.8.0) diagnostic: compute and record per-fixture density
    # ratio ρ_sample / ρ_ref (ρ_ref from .efr Material; ρ_sample from
    # sub_source mass / .efr Volume). This is documented in note_parts
    # but NO correction is applied to A_meas — implementing matrix
    # attenuation correction requires μ(E)/ρ tables for the ОИСН-16
    # matrix material across all measurement energies. See K-18 in
    # KNOWN_AND_FIXED_ISSUES.md for the methodology gap and
    # quantitative bias bounds. The diagnostic helps the reader
    # interpret which fixtures may be affected by matrix attenuation
    # mismatch.
    if (target_sub is not None and target_sub.mass_g is not None
            and unit.endswith("bq/kg")):
        # Look up reference volume/density for this geometry from the
        # .efr metadata cached on the EfficiencyCurve object (if the
        # parser exposes it; otherwise fall back to a small hard-coded
        # table compiled from manual .efr inspection 2024-05).
        REF_DENSITY = {
            "Маринелли":    (1000.0, 1.60),  # vol_ml, ρ_ref g/cm³
            "Дента-120мл":  ( 120.0, 1.66),
            "Петри-60мл":   (  60.0, 1.60),
        }
        if fx.geometry in REF_DENSITY:
            vol_ml, rho_ref = REF_DENSITY[fx.geometry]
            rho_sample = target_sub.mass_g / vol_ml
            rho_ratio = rho_sample / rho_ref
            note_parts.append(
                f"ρ_sample/ρ_ref={rho_ratio:.2f}"
            )

    # Background subtraction
    bg_sub = subtract_background(
        spec, bg_spec, user_confirmed_applicable=True,
    )
    spec_net = apply_subtraction_to_spectrum(spec, bg_sub)

    # Peak search with the .spe-derived FWHM polynomial
    fwhm_at = make_lsrm_fwhm_provider(spec_net)
    counts = np.asarray(spec_net.counts, dtype=np.float64)
    found = mariscotti_search(
        counts,
        fwhm_channels=fwhm_at,
        sigma_threshold=3.0,
    )
    note_parts.append(f"{len(found)} peaks")

    # Identification restricted to the target nuclide. We deliberately
    # SKIP disambiguate_identifications here: it implements
    # mixture-resolution rules (chain-vs-positron, intensity-ratio
    # proportionality for "rare" nuclides). In a known single-source
    # cert spectrum those rules either no-op (Co-60, Cs-137) or actively
    # reject the target (Na-22 fails proportionality vs the Mariscotti
    # significance ratio because Mariscotti sigma scales as
    # height/sqrt(B), not with peak area, and the check does not divide
    # by epsilon(E)). The decision is recorded in NOTES.
    id_res = identify_nuclides(
        found_peaks=found,
        spec=spec_net,
        candidate_nuclides=[fx.nuclide],
        compute_peak_areas=True,
        fwhm_at_channel=fwhm_at,
    )
    id_res, decons = apply_multiplet_deconvolution(
        id_res, spec_net, fwhm_at,
        overlap_threshold_fwhm=1.0,
    )
    if decons:
        note_parts.append(f"{len(decons)} multiplet(s)")

    target_id = next(
        (ni for ni in id_res.detected_nuclides if ni.nuclide == fx.nuclide),
        None,
    )
    if target_id is None or not target_id.matched_lines:
        return RunResult(fx.nuclide, fx.spe_filename, meas_dt,
                         src.reference_datetime, cert_A_reported_for_csv,
                         A_cert_at_meas, None, None, None, 0,
                         "; ".join(note_parts +
                                   ["nuclide not detected"]),
                         geometry=fx.geometry)

    # TCS — K-21 (v1.9.0): use per-geometry effective P/T to model
    # close-geometry cascade-coincidence enhancement. For point-5cm
    # the geometry factor is 1.0 (Gilmore reference), preserving the
    # F-31b behavior. For close-geometry samples (Marinelli, Дента,
    # Петри at 0 cm distance), the factor is < 1, scaling P/T down
    # which scales ε_T up and produces a larger TCS correction. See
    # GEOMETRY_PT_FACTOR in `gamma.physics.cascade_summing` and K-21
    # in KNOWN_AND_FIXED_ISSUES.md.
    pt_for_geom = peak_to_total_NaI_for_geometry(fx.geometry)
    tcs = compute_tcs_corrections(fx.nuclide, eff_curve,
                                  p_t_func=pt_for_geom)

    # Activity. min_intensity_pct=5.0 filters low-intensity outliers
    # (Eu-152's 3.12% I_pct line at 443.96 keV gets ~17x the weight of
    # the well-resolved 121.78 line because its tiny stat sigma
    # dominates Gaussian-on-NaI uncertainty propagation, dragging the
    # weighted mean down by ~60%). We drop to a 0 floor for nuclides
    # whose entire library catalog sits below 5% (Cd-109 has a single
    # 88 keV line at I=3.66%) -- otherwise no lines would survive.
    lib_lines = lib.get("lines", []) or []
    has_strong_line = any(
        (len(L) > 1 and float(L[1]) >= 5.0) for L in lib_lines
    )
    min_I = 5.0 if has_strong_line else 0.0
    act = compute_activity(
        target_id,
        efficiency_curve=eff_curve,
        live_time_s=live,
        from_bg_subtracted=True,
        coincidence_correction=tcs,
        decay_correction=False,    # we already decay-corrected the cert side
        measurement_datetime=meas_dt,
        reference_datetime=meas_dt,
        min_intensity_pct=min_I,
    )

    if not act.is_valid():
        return RunResult(fx.nuclide, fx.spe_filename, meas_dt,
                         src.reference_datetime, cert_A_reported_for_csv,
                         A_cert_at_meas, None, None, None,
                         len(act.lines_used),
                         "; ".join(note_parts +
                                   [f"activity not valid: {act.notes}"]),
                         geometry=fx.geometry)

    # K-20 (v1.9.0): apply self-attenuation correction for volume
    # samples in matrix geometries. Correction is a weighted mean over
    # the lines actually used by compute_activity, weighted by their
    # inverse-variance (which matches compute_activity's own
    # weighting). For point geometries (Точечная-5см/25см),
    # K20_REF_GEOMETRY has no entry and correction = 1.0.
    k20_corr_applied = 1.0
    if (fx.geometry in K20_REF_GEOMETRY
            and target_sub is not None
            and target_sub.mass_g is not None
            and unit.endswith("bq/kg")):
        vol_ml, rho_ref, t_cm = K20_REF_GEOMETRY[fx.geometry]
        rho_sample = target_sub.mass_g / vol_ml
        # weight each line by 1/σ²(A_line) to match compute_activity's
        # inverse-variance-weighted mean.
        E_keVs = [la.E_keV for la in act.lines_used]
        weights = []
        for la in act.lines_used:
            sigma_A = getattr(la, "sigma_A_Bq", None) or 0
            if sigma_A > 0:
                weights.append(1.0 / (sigma_A * sigma_A))
            else:
                weights.append(0.0)
        if any(w > 0 for w in weights) and E_keVs:
            k20_corr_applied = k20_weighted_correction(
                E_keVs, weights,
                rho_sample_g_cm3=rho_sample,
                rho_ref_g_cm3=rho_ref,
                thickness_cm=t_cm,
                composition=OISN_16_COMPOSITION,
            )
        # Annotate the note column
        note_parts.append(f"K20×{k20_corr_applied:.3f}")

    A_meas_corrected = act.A_Bq * k20_corr_applied
    sigma_meas_corrected = act.sigma_A_Bq * k20_corr_applied
    deviation = (A_meas_corrected - A_cert_at_meas) / A_cert_at_meas * 100.0

    if tcs:
        note_parts.append(f"TCS={len(tcs)} line(s)")
    sources = {la.E_keV: None for la in act.lines_used}
    return RunResult(
        fx.nuclide, fx.spe_filename, meas_dt,
        src.reference_datetime, cert_A_reported_for_csv, A_cert_at_meas,
        A_meas_corrected, sigma_meas_corrected, deviation,
        len(act.lines_used),
        "; ".join(note_parts),
        geometry=fx.geometry,
    )


# ---------------------------------------------------------------------------
# Pretty printing
# ---------------------------------------------------------------------------

def fmt_A(v: Optional[float]) -> str:
    if v is None or (isinstance(v, float) and (math.isnan(v) or v <= 0)):
        return "      —    "
    if v >= 1e5:
        return f"{v:>9.3e}"
    return f"{v:>10.1f}"


def fmt_dev(v: Optional[float]) -> str:
    if v is None:
        return "    —   "
    return f"{v:+7.2f}%"


# F-47a (v1.8.0): per-geometry polynomial degree tuning. Default
# degree=3 over-fits Точ-25см's 20-anchor data; degree=5 drops
# chi²/dof from 2.51 to 1.74. Other geometries are best at degree=3.
# Дента .efr has only 13 anchor points and is data-quality-limited
# (chi²/dof≈15 for all degrees 1-5); fit reported as accepted
# limitation K-17, no degree tuning improves it. See
# KNOWN_AND_FIXED_ISSUES.md K-17 for the quantitative analysis.
EFF_DEGREE = {
    "Точечная-5см":  3,   # chi²/dof=6.95
    "Точечная-25см": 5,   # chi²/dof=1.74 (was 2.51 at degree=3)
    "Дента-120мл":   3,   # chi²/dof≈15 invariant of degree (K-17)
    "Петри-60мл":    4,   # chi²/dof=14.28 (was 15.04 at degree=3)
    "Маринелли":     3,   # chi²/dof=3.72
}


def _resolve_geometry_resources(geom: str, *, cache):
    """Return (eff_curve, bg_spec, cert) for the given geometry token,
    lazily loading and caching the three resources per geometry. The
    cache is mutated in place. Raises KeyError if the geometry is not
    registered in EFF_PATHS/BG_PATHS/CERT_PATHS."""
    if geom in cache:
        return cache[geom]
    eff_path = EFF_PATHS[geom]
    bg_path = BG_PATHS[geom]
    cert_path = CERT_PATHS[geom]
    degree = EFF_DEGREE.get(geom, 3)
    print(f"[{geom}] loading resources …")
    print(f"  eff = {eff_path.name}  (degree={degree})")
    eff = fit_efficiency_from_efr_file(str(eff_path), degree=degree)
    print(f"    -> {eff.E_min_keV:.1f}–{eff.E_max_keV:.1f} keV, "
          f"n_pts={eff.n_points_used}, chi²/dof={eff.chi2_per_dof:.2f}")
    print(f"  bg  = {bg_path.name}")
    bg = read_spectrum(str(bg_path))
    print(f"    -> live_time={bg.live_time:.0f}s, n_channels={bg.n_channels}")
    print(f"  cert = {cert_path.name}")
    cert = read_certificate_file(cert_path)
    print(f"    -> {len(cert.sources)} sources, sigma={cert.confidence_sigma}")
    cache[geom] = (eff, bg, cert)
    return cache[geom]


def main() -> int:
    # Lazy per-geometry resource cache. Each geometry's (eff, bg, cert)
    # is loaded once on first fixture-hit and reused for all subsequent
    # fixtures of that geometry — saves ~3-5 seconds per geometry over
    # naive reload.
    geom_cache: dict = {}

    rows = []
    for fx in FIXTURES:
        try:
            eff, bg, cert = _resolve_geometry_resources(fx.geometry, cache=geom_cache)
        except KeyError as e:
            print(f"-- {fx.nuclide:>7} ({fx.spe_filename}) --")
            print(f"    SKIP: unknown geometry {fx.geometry!r}: {e}")
            rows.append(RunResult(fx.nuclide, fx.spe_filename, None, None,
                                  None, None, None, None, None, 0,
                                  f"unknown geometry {fx.geometry!r}",
                                  geometry=fx.geometry))
            continue
        print(f"-- {fx.nuclide:>7} ({fx.spe_filename}) [{fx.geometry}] --")
        try:
            r = run_one(fx, eff, bg, cert)
        except Exception as e:
            r = RunResult(fx.nuclide, fx.spe_filename, None, None, None,
                          None, None, None, None, 0, f"EXC: {e}",
                          geometry=fx.geometry)
        rows.append(r)
        print(f"    cert A={fmt_A(r.cert_A_Bq)}  "
              f"cert@meas A={fmt_A(r.cert_A_decayed_Bq)}  "
              f"meas A={fmt_A(r.measured_A_Bq)}  "
              f"dev={fmt_dev(r.deviation_pct)}  "
              f"({r.note})")

    # Pretty table
    print()
    print("=" * 110)
    print("Матрица отклонений: измеренная активность vs сертификат "
          "(decay-corrected to spectrum date)")
    print("=" * 110)
    header = (f"{'Нуклид':>7} | {'A_cert, Bq':>11} | "
              f"{'A_cert@meas, Bq':>15} | {'A_изм, Bq':>11} | "
              f"{'σ_A, Bq':>9} | {'Δ, %':>8} | {'n_лин':>5} | comment")
    # F-46 (v1.7.25): group rows by geometry in the printed table.
    # Within each geometry block, preserve fixture insertion order.
    by_geom: dict = {}
    for r in rows:
        by_geom.setdefault(r.geometry, []).append(r)
    print(header)
    print("-" * 110)
    for geom in EFF_PATHS.keys():  # iterate in canonical order
        geom_rows = by_geom.get(geom, [])
        if not geom_rows:
            continue
        print(f"--- {geom} ({len(geom_rows)} fixtures) ---")
        for r in geom_rows:
            print(f"{r.nuclide:>7} | {fmt_A(r.cert_A_Bq):>11} | "
                  f"{fmt_A(r.cert_A_decayed_Bq):>15} | "
                  f"{fmt_A(r.measured_A_Bq):>11} | "
                  f"{fmt_A(r.measured_sigma_Bq):>9} | "
                  f"{fmt_dev(r.deviation_pct):>8} | "
                  f"{r.n_lines_used:>5} | {r.note}")
    print("=" * 110)

    # F-46 per-geometry summary table. Compute mean / max |Δ| within
    # each geometry block. Single-geometry harnesses (v1.7.23 and
    # earlier) will print a single-row block.
    print()
    print("=" * 90)
    print("Per-geometry summary")
    print("=" * 90)
    print(f"{'Geometry':<16} | {'n_total':>7} | {'n_meas':>6} | "
          f"{'mean|Δ|':>8} | {'max|Δ|':>8} | "
          f"{'efr_chi²/dof':>12} | {'bg_file':>30}")
    print("-" * 90)
    for geom in EFF_PATHS.keys():
        geom_rows = by_geom.get(geom, [])
        if not geom_rows:
            continue
        measurable = [r for r in geom_rows
                      if r.deviation_pct is not None]
        if measurable:
            devs = [abs(r.deviation_pct) for r in measurable]
            mean_d = sum(devs) / len(devs)
            max_d = max(devs)
            mean_str = f"{mean_d:>7.2f}%"
            max_str = f"{max_d:>7.2f}%"
        else:
            mean_str = max_str = "    —   "
        # F-46d: per-geometry efficiency curve chi²/dof from the cache.
        # (geom_cache populated above stays in scope.)
        chi = geom_cache.get(geom, (None,))[0]
        chi_str = f"{chi.chi2_per_dof:>11.2f}" if chi else "      —    "
        bg_name = BG_PATHS[geom].name[:30]
        print(f"{geom:<16} | {len(geom_rows):>7} | {len(measurable):>6} | "
              f"{mean_str} | {max_str} | "
              f"{chi_str} | {bg_name:>30}")
    print("=" * 90)

    # CSV
    csv_path = ROOT / "cert_validation_matrix.csv"
    with open(csv_path, "w", encoding="utf-8") as fh:
        fh.write("nuclide,spe_filename,meas_datetime,cert_ref_datetime,"
                 "cert_A_Bq,cert_A_at_meas_Bq,measured_A_Bq,measured_sigma_Bq,"
                 "deviation_pct,n_lines_used,note\n")
        for r in rows:
            fh.write(",".join([
                r.nuclide,
                r.spe_filename,
                r.measurement_datetime.isoformat() if r.measurement_datetime else "",
                r.cert_reference_datetime.isoformat() if r.cert_reference_datetime else "",
                f"{r.cert_A_Bq:.6g}" if r.cert_A_Bq is not None else "",
                f"{r.cert_A_decayed_Bq:.6g}" if r.cert_A_decayed_Bq is not None else "",
                f"{r.measured_A_Bq:.6g}" if r.measured_A_Bq is not None else "",
                f"{r.measured_sigma_Bq:.6g}" if r.measured_sigma_Bq is not None else "",
                f"{r.deviation_pct:.4f}" if r.deviation_pct is not None else "",
                str(r.n_lines_used),
                r.note.replace(",", ";"),
            ]) + "\n")
    print(f"\nCSV saved: {csv_path}")

    # F-41 (v1.7.19): chain-proxy cross-validation block. When two or
    # more independent daughters of the same parent (e.g. Pb-212 AND
    # Tl-208 both measure Th-228) produce a measurable activity on the
    # SAME physical source (= same .spe file), print them side-by-side
    # with their ratio. Consistency within a few %% is independent
    # evidence that the chain-proxy methodology + efficiency curve +
    # TCS corrections all behave coherently across widely separated
    # energies (238 keV for Pb-212 vs 583/860/2614 keV for Tl-208).
    # Inconsistency >5 %% would flag a systematic error.
    #
    # F-46 (v1.7.25): pair grouping by (cert_nuclide, spe_filename)
    # instead of cert_nuclide alone. Otherwise ratios cross physical
    # sources/geometries and the comparison loses interpretability
    # (different cert ref dates, different masses → ratios reflect
    # cert metadata as well as detector response).
    same_source_proxies: dict = {}
    for r in rows:
        fx = next((f for f in FIXTURES if f.spe_filename == r.spe_filename
                   and f.nuclide == r.nuclide
                   and f.geometry == r.geometry), None)
        if fx and fx.cert_nuclide and r.measured_A_Bq is not None:
            key = (fx.cert_nuclide, r.spe_filename, r.geometry)
            same_source_proxies.setdefault(key, []).append((r, fx))
    if any(len(v) >= 2 for v in same_source_proxies.values()):
        print()
        print("=" * 110)
        print("Cross-validation of chain proxies (F-41 / F-46) — same .spe, ≥2 daughters")
        print("=" * 110)
        for (parent, spe, geom), proxy_rows in same_source_proxies.items():
            if len(proxy_rows) < 2:
                continue
            print(f"\nParent: {parent}   [{geom}]   spe={spe}")
            print(f"  {'daughter':>10}  {'A_meas, Bq':>12}  "
                  f"{'A_cert@meas, Bq':>16}  {'Δ vs cert, %':>13}")
            for r, fx in proxy_rows:
                print(f"  {r.nuclide:>10}  "
                      f"{fmt_A(r.measured_A_Bq):>12}  "
                      f"{fmt_A(r.cert_A_decayed_Bq):>16}  "
                      f"{fmt_dev(r.deviation_pct):>13}")
            # Pairwise ratio
            for i in range(len(proxy_rows)):
                for j in range(i + 1, len(proxy_rows)):
                    a, _ = proxy_rows[i]
                    b, _ = proxy_rows[j]
                    if (a.measured_A_Bq and b.measured_A_Bq
                            and b.measured_A_Bq > 0):
                        ratio = a.measured_A_Bq / b.measured_A_Bq
                        diff_pct = (ratio - 1.0) * 100.0
                        flag = "  OK" if abs(diff_pct) <= 5.0 else "  >5%"
                        print(f"  ratio {a.nuclide}/{b.nuclide} = "
                              f"{ratio:.4f}  ({diff_pct:+.2f}%){flag}")
        print("=" * 110)

    # Exit code 0 if all measured fixtures within ±10%, 1 otherwise
    deviations = [r.deviation_pct for r in rows if r.deviation_pct is not None]
    if not deviations:
        print("WARNING: no fixture produced a measurable activity")
        return 1
    bad = [r for r in rows
           if r.deviation_pct is not None and abs(r.deviation_pct) > 10.0]
    if bad:
        print(f"\n{len(bad)} fixture(s) deviate by >10%: "
              f"{', '.join(r.nuclide for r in bad)}")
    measured = [r for r in rows if r.deviation_pct is not None]
    print(f"\nMeasured {len(measured)}/{len(rows)} fixtures.")
    print(f"Mean |Δ| = {sum(abs(d) for d in deviations)/len(deviations):.2f}%")
    print(f"Max  |Δ| = {max(abs(d) for d in deviations):.2f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
