"""Geometry-class enum for True Coincidence Summing (TCS) lookup.

This module is the canonical source of geometry identifiers used by
``data/tsc_lookup.json`` (the TSCF lookup table populated in v1.19.0
Phase 4 from Giubrone et al. 2016 and Ordoñez et al. 2019).

Per ``PLAN_v1_18_32_to_v1_19_0_TCS_INTEGRATION.md`` §3 Phase 4 step 3
(geometry classes enum), TSCF is primarily a function of detector
solid angle (sample geometry), with sample matrix density a much
weaker secondary axis (Ordoñez 2019 §3, ``ordonez_2025.txt:454-458``).
Therefore the v1.19.0 lookup table is indexed by ``geometry`` first;
matrix density is recorded per-entry as a metadata field
(``matrix_density_g_cm3``) for future per-density interpolation.

Coverage in v1.19.0
-------------------

Populated geometries (have at least one TSCF entry in
``data/tsc_lookup.json`` for detector_class ``HPGe-coaxial``):

* ``PETRI_25ML_WATER`` -- Ordoñez 2019 PS25W (point-like)
* ``PETRI_100ML_WATER`` -- Ordoñez 2019 PS100W = Giubrone 2016 PGAQ
* ``PETRI_100ML_SEA_SAND`` -- Ordoñez 2019 PS100SS = Giubrone 2016 PGSI
* ``PETRI_100ML_ZR_SAND`` -- Ordoñez 2019 PS100ZrS (high-density ZrSiO4)
* ``PETRI_15ML_AQUEOUS`` -- Giubrone 2016 PPAQ (close geometry, calibration sources only)
* ``MARINELLI_500ML_WATER`` -- Ordoñez 2019 MS500W (Marinelli beaker)

Stub geometries (declared for API stability, NO TSCF data in v1.19.0;
flagged for v1.19.1+ Garcia-Talavera / own-MC research per
``PLAN_v1_18_32_to_v1_19_0_TCS_INTEGRATION.md`` §9.2 deferred items):

* ``MARINELLI_1000ML`` -- larger Marinelli geometry
* ``POINT_SOURCE_ENDCAP`` -- point source on endcap
* ``POINT_SOURCE_10CM`` -- point source 10 cm
* ``POINT_SOURCE_25CM`` -- point source 25 cm
* ``WELL_COUNTER`` -- well-type detector

The string values of every enum member are exactly the strings used in
the ``geometry`` field of ``data/tsc_lookup.json`` entries; tests in
``tests/data/test_tsc_lookup.py`` enforce this 1-to-1 correspondence.

References
----------
* Giubrone et al. (2016) J. Environ. Radioact. 158-159:114-118 --
  PGAQ / PGSI / PPAQ definitions
  (``audit/_drafts/giubrone_2016_tsc_extracts.md`` §1)
* Ordoñez-Rodenas et al. (2019) Radiat. Phys. Chem. 155:244-247 --
  PS25W / PS100W / PS100SS / PS100ZrS / MS500W definitions
  (``audit/_drafts/ordonez_2019_tsc_extracts.md`` §3.2)
"""

from __future__ import annotations

from enum import Enum


class GeometryClass(str, Enum):
    """Sample-geometry identifiers for TSCF lookup.

    Inherits from ``str`` so that ``GeometryClass.PETRI_25ML_WATER ==
    "petri_25ml_water"`` evaluates True; this lets JSON-loaded values
    compare directly to enum members without explicit casting.
    """

    # --- Populated (TSCF entries exist in data/tsc_lookup.json) ---
    PETRI_25ML_WATER = "petri_25ml_water"
    PETRI_100ML_WATER = "petri_100ml_water"
    PETRI_100ML_SEA_SAND = "petri_100ml_sea_sand"
    PETRI_100ML_ZR_SAND = "petri_100ml_zr_sand"
    PETRI_15ML_AQUEOUS = "petri_15ml_aqueous"
    MARINELLI_500ML_WATER = "marinelli_500ml_water"

    # --- Stub (no TSCF data in v1.19.0; v1.19.1+ research candidates) ---
    MARINELLI_1000ML = "marinelli_1000ml"
    POINT_SOURCE_ENDCAP = "point_source_endcap"
    POINT_SOURCE_10CM = "point_source_10cm"
    POINT_SOURCE_25CM = "point_source_25cm"
    WELL_COUNTER = "well_counter"


# Public set of populated geometries (have v1.19.0 TSCF data).
# Used by tests and (in v1.19.1+) by the TCS correction pipeline to
# distinguish "no entry" (lookup miss within a populated geometry) from
# "geometry not calibrated" (stub geometry, warning required).
POPULATED_GEOMETRIES = frozenset(
    {
        GeometryClass.PETRI_25ML_WATER,
        GeometryClass.PETRI_100ML_WATER,
        GeometryClass.PETRI_100ML_SEA_SAND,
        GeometryClass.PETRI_100ML_ZR_SAND,
        GeometryClass.PETRI_15ML_AQUEOUS,
        GeometryClass.MARINELLI_500ML_WATER,
    }
)

STUB_GEOMETRIES = frozenset(
    {
        GeometryClass.MARINELLI_1000ML,
        GeometryClass.POINT_SOURCE_ENDCAP,
        GeometryClass.POINT_SOURCE_10CM,
        GeometryClass.POINT_SOURCE_25CM,
        GeometryClass.WELL_COUNTER,
    }
)

__all__ = ["GeometryClass", "POPULATED_GEOMETRIES", "STUB_GEOMETRIES"]
