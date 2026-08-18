"""
Auto-load efficiency curve by canonical geometry name (F-57 / v1.11.1).

Looks up a fitted efficiency curve from the active detector's efficiency
archive based on the canonical geometry name resolved via
`gamma.data.aliases`. Currently the only supported detector is **Gamma-1S**
— its .efr files live under
`detectors/Gamma-1S/efficiency/<unit>/<geometry>.efr` and the root is
exposed by `gamma.detectors.gamma1s.EFFICIENCY_DIR` (F-83 / v1.12.0).

Lookup procedure:
  1. Determine canonical geometry: `gamma.data.aliases.canonicalize("geometry", raw)`
  2. Map canonical name → expected .efr filename pattern
  3. Walk `EFFICIENCY_DIR/<any-unit-dir>/` and try to match
  4. Load + fit via `gamma.calibration.efficiency.fit_efficiency_from_efr_file`
  5. Cache by (detector, geometry) tuple for the process lifetime

Return-value discipline (DEEP-01, v1.26.2 — Project #5 wave 2 P1-1):
``load_efficiency_for_geometry`` distinguishes three outcomes:

  * **No candidate .efr found** → returns ``None`` (silent — operator
    pipeline degrades to qualitative identification only).
  * **Candidate found AND fit succeeded** → returns the
    :class:`~gamma.calibration.efficiency.EfficiencyCurve` instance.
  * **Candidate found BUT fit failed** (broken/truncated .efr, parser
    raised, fit raised) → returns the module-level singleton
    :data:`EFFICIENCY_FIT_FAILED` AND emits ``logger.warning`` naming the
    offending file (basename only per F-115) + exception. The sentinel is
    explicitly *not* ``None`` so call sites can branch on it; treating it
    as truthy ``not None`` without isinstance-checking will fall back to
    the qualitative pipeline but still surface the warning in logs.

Before DEEP-01 the fit-failed branch silently returned ``None`` —
indistinguishable from "no file found" — so an operator running on a
geometry with a corrupted .efr would unknowingly ship an
efficiency-uncorrected report.
"""

from __future__ import annotations

import logging
from functools import lru_cache
from pathlib import Path
from typing import Optional, Tuple, Union

from gamma.calibration.efficiency import (
    fit_efficiency_from_efr_file, EfficiencyCurve,
)
from gamma.data.aliases import canonicalize
from gamma.detectors.gamma1s import EFFICIENCY_DIR as _EFR_ROOT


logger = logging.getLogger(__name__)


class _EfficiencyFitFailedSentinel:
    """
    Singleton marker for "candidate .efr found but fit failed" (DEEP-01).

    Distinct from ``None`` (which means "no candidate found"). Frozen via
    ``__slots__ = ()`` so it cannot accumulate state and ``__bool__``
    is False so legacy ``if curve:`` guards keep degrading gracefully —
    only ``is`` checks (or explicit ``isinstance``) lift the curtain.
    """

    __slots__ = ()

    def __repr__(self) -> str:  # pragma: no cover - debug aid
        return "<EFFICIENCY_FIT_FAILED>"

    def __bool__(self) -> bool:
        return False


#: Module-level singleton returned by :func:`load_efficiency_for_geometry`
#: when a candidate .efr was located but ``fit_efficiency_from_efr_file``
#: raised. Compare with ``is`` (identity), never equality.
EFFICIENCY_FIT_FAILED = _EfficiencyFitFailedSentinel()


# Map canonical geometry → list of substring patterns to search for in
# .efr filenames. The first matching .efr is used.
_GEOMETRY_FILENAME_PATTERNS = {
    "marinelli_1L":  ["Маринелли", "Marinelli", "marinelli"],
    "marinelli_05L": ["Маринелли-0.5", "Маринелли_0.5", "Marinelli 0.5", "Marinelli500"],
    "denta_120mL":   ["Дента", "Denta"],
    "petri_60mL":    ["Петри", "Petri"],
    "point_25cm":    ["Точечная-25", "Точ-25", "Point25", "point25", "25cm"],
    "point_10cm":    ["Точечная-10", "Точ-10", "Point10", "10cm"],
    "point_5cm":     ["Точечная-5", "Точ-5", "Point5", "5cm"],
    "point_0cm":     ["0cm", "0см", "контакт", "contact"],
    "well":          ["Колодец", "well", "Well"],
}


@lru_cache(maxsize=64)
def find_efr_file(
    geometry_raw: str,
    detector_raw: str = "",
) -> Optional[str]:
    """
    Locate the .efr file matching a geometry and (optionally) detector.

    Args:
        geometry_raw: raw geometry text (e.g. "Маринелли 1л", "Marinelli")
        detector_raw: raw detector text (used as a hint to pick a subdir
            when multiple detectors are present)

    Returns:
        Absolute path string to the matching .efr file, or None.
    """
    if not _EFR_ROOT.is_dir():
        return None

    geom_canonical = canonicalize("geometry", geometry_raw)
    if not geom_canonical:
        return None
    patterns = _GEOMETRY_FILENAME_PATTERNS.get(geom_canonical)
    if not patterns:
        return None

    # Determine candidate detector subdirectories
    det_canonical = canonicalize("detector", detector_raw) if detector_raw else None
    candidate_dirs = []
    for d in sorted(_EFR_ROOT.iterdir()):
        if not d.is_dir():
            continue
        if det_canonical:
            # Match detector by canonical name embedded in directory name
            d_canon = canonicalize("detector", d.name)
            if d_canon == det_canonical:
                candidate_dirs.insert(0, d)  # prefer
            else:
                candidate_dirs.append(d)
        else:
            candidate_dirs.append(d)

    for d in candidate_dirs:
        for pat in patterns:
            for f in d.glob("*.efr"):
                if pat.lower() in f.name.lower():
                    return str(f)
    return None


@lru_cache(maxsize=32)
def load_efficiency_for_geometry(
    geometry_raw: str,
    detector_raw: str = "",
    *,
    degree: int = 3,
) -> Union[Optional[EfficiencyCurve], _EfficiencyFitFailedSentinel]:
    """
    Find + fit + return the efficiency curve for a given geometry.

    Three return states (DEEP-01, v1.26.2):

      * ``None``                    — no matching .efr file found in
                                      ``detectors/Gamma-1S/efficiency/``
                                      (silent).
      * :class:`EfficiencyCurve`    — a candidate was located and the
                                      log-log polynomial fit succeeded.
      * :data:`EFFICIENCY_FIT_FAILED` — a candidate was located but the
                                      reader or fit raised. The function
                                      additionally emits a ``logger.warning``
                                      with the file basename (F-115) and
                                      the exception class/message.

    The sentinel evaluates falsy in boolean context, so legacy
    ``if curve is not None`` and ``if curve`` guards both degrade
    gracefully (efficiency-uncorrected). New code SHOULD branch on
    ``curve is EFFICIENCY_FIT_FAILED`` to surface the failure to the
    operator.
    """
    efr_path = find_efr_file(geometry_raw, detector_raw)
    if not efr_path:
        return None
    try:
        return fit_efficiency_from_efr_file(efr_path, degree=degree)
    except Exception as exc:
        # F-115: basename only — never log absolute operator paths.
        # The underlying readers/fitters embed the full path in their
        # exception messages (e.g. efficiency.py:278
        # ``ValueError("No efficiency points in: <abs path>")``); scrub
        # it out before logging so the warning stays portable.
        basename = Path(efr_path).name
        exc_msg = str(exc).replace(efr_path, basename)
        logger.warning(
            "efficiency_autoload: fit failed for %r (%s: %s) — "
            "returning EFFICIENCY_FIT_FAILED sentinel; "
            "downstream activity computation will skip efficiency correction "
            "and surface the failure to the operator.",
            basename,
            type(exc).__name__,
            exc_msg,
        )
        return EFFICIENCY_FIT_FAILED


__all__ = [
    "find_efr_file",
    "load_efficiency_for_geometry",
    "EFFICIENCY_FIT_FAILED",
]
