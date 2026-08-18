"""gamma.detectors.gamma1s — path resolver for the Gamma-1S spectrometric complex.

The Gamma-1S complex (internal label «Гамма-1С» — cyrillic 'С' U+0421;
ASCII canonical 'Gamma-1S' = transliteration of cyrillic «С»; legacy
'Gamma-1C' alias kept for backward compatibility per F2-A 2026-06-21)
consists of:

- **Detector head:** УДС-ГЦ-63×63 (NaI(Tl) 63×63 mm crystal, БДЭГ-63×63-USB) by Aspect.
- **DAQ software:** Lsrm SpectraLine.
- **Canonical alias:** ``Gamma-1S`` (case-insensitive; synonyms in
  ``data/aliases.json`` → ``detector.Gamma-1S``: УДС-ГЦ-63×63, БДЭГ-63×63,
  Колибри-1М, Гамма-1С, ...).

All detector-specific assets live under ``detectors/Gamma-1S/`` at the project
root. The constants exposed by this module are the **canonical paths** to those
assets — calling code MUST go through them rather than hard-coding string
fragments, because future detector subtrees follow the same naming convention.

Layout under ``detectors/Gamma-1S/``::

    certificates/                 # passports of standard sources
    efficiency/                   # .efr efficiency curves per geometry
    reference_spectra/            # .spe reference measurements
    lsrm-libraries/               # LSRM SpectraLine nuclide libraries
    references/
        05_intrinsic_detector_activity.md
        07_dead_time_correction.md
    data/
        averaged_backgrounds/     # averaged-background .spe per geometry
        secondary_peaks.json      # Cs-137 + K-40 catalog
        secondary_peaks_v2.json   # 9-isotope rich catalog

Isolation policy (v1.12.0 / F-83): this module is the **only** place that knows
the on-disk layout for Gamma-1S. Algorithms in ``gamma.peaks``,
``gamma.identification``, ``gamma.calibration``, ``gamma.activity`` and
``gamma.physics`` import constants from here.
"""

from __future__ import annotations

from pathlib import Path

__all__ = [
    "DETECTOR_NAME",
    "DETECTOR_CANONICAL",
    "DEFAULT_REFERENCE_SUBDIR",
    "DETECTOR_ROOT",
    "CERTIFICATES_DIR",
    "EFFICIENCY_DIR",
    "REFERENCE_SPECTRA_DIR",
    "LSRM_LIBRARIES_DIR",
    "REFERENCES_DIR",
    "DATA_DIR",
    "AVERAGED_BACKGROUNDS_DIR",
    "SECONDARY_PEAKS_PATH",
    "SECONDARY_PEAKS_V2_PATH",
    "INTRINSIC_ACTIVITY_REF",
    "DEAD_TIME_REF",
    "DEFAULT_REFERENCE_DIR",
    "DEFAULT_EFFICIENCY_DIR",
]

# ---------------------------------------------------------------------------
# Identity
# ---------------------------------------------------------------------------

#: Canonical detector name used in reports and registries.
DETECTOR_NAME: str = "Gamma-1S"

#: Same as :data:`DETECTOR_NAME`; provided so callers can read it as the
#: canonical alias (matches ``data/aliases.json`` → ``detector.Gamma-1S``).
DETECTOR_CANONICAL: str = "Gamma-1S"

#: Subdirectory name used by both ``efficiency/`` and ``reference_spectra/``
#: for the primary calibrated Gamma-1S unit (serial SN-01). Reused widely in
#: tests and helpers.
DEFAULT_REFERENCE_SUBDIR: str = "Gamma-1S_NaI_63x63_USB_SN-01"

# ---------------------------------------------------------------------------
# Roots
# ---------------------------------------------------------------------------

# scripts/gamma/detectors/gamma1s.py  →  project root is parents[3]
DETECTOR_ROOT: Path = Path(__file__).resolve().parents[3] / "detectors" / "Gamma-1S"

# ---------------------------------------------------------------------------
# Top-level subdirectories
# ---------------------------------------------------------------------------

#: Source-passport spreadsheets / PDFs.
CERTIFICATES_DIR: Path = DETECTOR_ROOT / "certificates"

#: .efr efficiency curve archive (per-geometry subfolders).
EFFICIENCY_DIR: Path = DETECTOR_ROOT / "efficiency"

#: .spe reference spectra archive (per-geometry subfolders).
REFERENCE_SPECTRA_DIR: Path = DETECTOR_ROOT / "reference_spectra"

#: LSRM SpectraLine nuclide library bundle (.lib / .nlb / ...).
LSRM_LIBRARIES_DIR: Path = DETECTOR_ROOT / "lsrm-libraries"

#: Markdown references that are Gamma-1S-specific (intrinsic activity, dead-time
#: A/B coefficients).
REFERENCES_DIR: Path = DETECTOR_ROOT / "references"

#: Data subtree (averaged backgrounds + secondary-peak catalogs).
DATA_DIR: Path = DETECTOR_ROOT / "data"

# ---------------------------------------------------------------------------
# Common artefact paths
# ---------------------------------------------------------------------------

#: Averaged-background .spe files per geometry (built by ``build_averaged_backgrounds.py``).
AVERAGED_BACKGROUNDS_DIR: Path = DATA_DIR / "averaged_backgrounds"

#: Cs-137 + K-40 secondary-peak catalog (built by ``analyze_secondaries.py``).
SECONDARY_PEAKS_PATH: Path = DATA_DIR / "secondary_peaks.json"

#: 9-isotope rich secondary-peak catalog (built by ``analyze_problem_isotopes.py``;
#: consumed by ``gamma.identification.disambiguate`` F-40 rule and the residual
#: classifier).
SECONDARY_PEAKS_V2_PATH: Path = DATA_DIR / "secondary_peaks_v2.json"

#: Intrinsic-activity methodology reference (NaI(Tl) 63×63 specific).
INTRINSIC_ACTIVITY_REF: Path = REFERENCES_DIR / "05_intrinsic_detector_activity.md"

#: Dead-time correction methodology + A,B coefficients (Lsrm §15 specific to the УДС-ГЦ).
DEAD_TIME_REF: Path = REFERENCES_DIR / "07_dead_time_correction.md"

# ---------------------------------------------------------------------------
# Convenience: default per-unit subdirectories (serial SN-01)
# ---------------------------------------------------------------------------

#: Default reference-spectra subdirectory for serial SN-01.
#: F-331 / v1.18.18.5 — после реорганизации `reference_spectra/`
#: исторический flat-каталог `<SUBDIR>/` упразднён; legacy spectra
#: лежат в `reference_spectra/archive/` (canonical kits — в
#: `reference_spectra/reference_kits/`). Path константа сохраняет
#: интерфейс ('исходные раздавленные файлы серии SN-01') и
#: указывает на archive, чтобы старые скрипты + тесты продолжали
#: работать.
DEFAULT_REFERENCE_DIR: Path = REFERENCE_SPECTRA_DIR / "archive"

#: Default efficiency subdirectory for serial SN-01 (matches the .efr archive layout).
DEFAULT_EFFICIENCY_DIR: Path = EFFICIENCY_DIR / DEFAULT_REFERENCE_SUBDIR

#: F-331 / v1.18.18.5 — canonical reference-kit root (per-geometry,
#: per-nuclide leaves with sample + matching background).
DEFAULT_KIT_DIR: Path = REFERENCE_SPECTRA_DIR / "reference_kits"
