"""
F-274 (v1.17.11, T-035) — Source-set validator for energy calibration.

Some calibration source / detector-class combinations are problematic
even though they don't fail outright. This module emits **warnings**
(never errors) so the operator can plan acquisition time or substitute
a more suitable source.

Background
----------
A common silent-failure mode on NaI Gamma-1S 63x63:

  Operator runs a **dedicated source measurement** (Eu-152 or Th-232)
  expecting standard peaks. If acquisition time is too short, weak
  high-energy lines (Tl-208 2614, Eu-152 1408) don't accumulate
  enough counts to be Mariscotti-found, and the polynomial fit
  silently drops to ≤ MIN_ANCHORS or fails the range gate (F-270).

  The user sees only "calibration failed" without knowing it was
  acquisition-time-bound. This validator surfaces the expected
  difficulty BEFORE the measurement.

Background-mode 7-line ЕРН calibration (long acquisitions, hours) is
unaffected — the 2614 line accumulates over time and is reliably
found via ``seven_line_check.run_seven_line_check``.

Returns
-------
``list[SourceSetWarning]`` — never raises. Empty list means OK.

References
----------
- ЛСРМ Алгоритмические основы 2022 § 8 (calibration source choice)
- Gilmore & Joss "Practical Gamma-ray Spectrometry" 3rd Ed., § 6.2
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional


# Минимальная рекомендованная длительность набора для пар (source, line, detector_class)
# в секундах. Используется как WARNING-порог: фактическое время < этого → warning.
# Значения взяты из практики работы с Gamma-1S 63×63 (ЛСРМ tables).
_MIN_ACQUISITION_SECONDS = {
    # (source_label, line_E_keV_min, detector_class) → seconds
    ("Th-232", 2614.0, "NaI"):  1800,   # 2614 keV NaI: 30 мин при ~1 кБк
    ("Th-232", 2614.0, "CsI"):  1800,
    ("Eu-152", 1408.0, "NaI"):  900,    # 1408 keV NaI: 15 мин при ~5 кБк (I_γ=21%)
    ("Eu-152", 1408.0, "CsI"):  900,
    ("Co-60",  1332.5, "NaI"):  600,    # 1332.5 keV NaI: 10 мин при ~1 кБк
    ("K-40",   1460.8, "NaI"):  3600,   # K-40 4 кБк ОИСН — 1 час
}

# Источники, которые на NaI/CsI системно требуют осторожности
# (плохо разрешаемые мультиплеты или очень слабые ключевые линии).
_NAI_PROBLEMATIC_SOURCES = {
    "Eu-152": (
        "Eu-152 на NaI/CsI имеет до 6 близких линий в диапазоне "
        "121–1408 кэВ. Разрешение Gamma-1S 63×63 ~7% при 662 кэВ ⇒ "
        "линии 245+295, 344+367+411, 779+867+964+1086+1112 сливаются. "
        "Рекомендуется ОСГИ-3 (Cs-137 + Co-60 + K-40 + Am-241) для "
        "чистых одиночных anchor'ов."
    ),
}


@dataclass(frozen=True)
class SourceSetWarning:
    """One warning about a (source, line, detector) combination."""
    severity: str          # "info" | "warn" | "advisory"
    code: str              # "F274-acquisition" | "F274-source-mismatch"
    message: str
    source_label: str = ""
    line_keV: Optional[float] = None
    detector_class: str = ""


def validate_source_set_for_detector(
    *,
    source_label: str,
    anchor_energies_keV: Iterable[float],
    detector_class: str,
    acquisition_time_s: Optional[float] = None,
) -> List[SourceSetWarning]:
    """Validate a source / detector / acquisition-time combination.

    Parameters
    ----------
    source_label : str
        Calibration source name (e.g. "Th-232", "Eu-152", "Co-60",
        "K-40", "ОСГИ-3", "ОСГИ-К", "Cs-137"). Used for matching
        against the recommendation table.
    anchor_energies_keV : iterable of float
        The list of line energies you intend to use as calibration
        anchors. Each one is checked against the per-(source,line,det)
        recommendation table.
    detector_class : str
        One of "HPGe", "CdZnTe", "LaBr3", "CeBr3", "NaI", "CsI".
    acquisition_time_s : float or None
        Live time of the planned acquisition in seconds. If None, the
        check is skipped (only source-mismatch warnings emitted).

    Returns
    -------
    list of SourceSetWarning
        Empty list means no warnings.
    """
    warnings: List[SourceSetWarning] = []
    src = str(source_label).strip()
    det = str(detector_class).strip()

    # Source-mismatch advisory (e.g. Eu-152 on NaI)
    if det in ("NaI", "CsI") and src in _NAI_PROBLEMATIC_SOURCES:
        warnings.append(SourceSetWarning(
            severity="advisory",
            code="F274-source-mismatch",
            message=_NAI_PROBLEMATIC_SOURCES[src],
            source_label=src,
            detector_class=det,
        ))

    # Per-line acquisition-time check (tolerance-based lookup ±2 keV)
    if acquisition_time_s is not None and acquisition_time_s > 0:
        for E in anchor_energies_keV:
            try:
                E_f = float(E)
            except Exception:
                continue
            # Найти ближайшую запись в таблице (tolerance 2 кэВ)
            min_t = None
            line_label = None
            for (s_key, e_key, d_key), t_key in _MIN_ACQUISITION_SECONDS.items():
                if s_key == src and d_key == det and abs(e_key - E_f) <= 2.0:
                    min_t = t_key
                    line_label = e_key
                    break
            if min_t is not None and acquisition_time_s < min_t:
                warnings.append(SourceSetWarning(
                    severity="warn",
                    code="F274-acquisition",
                    message=(
                        f"Линия {E_f:.1f} кэВ источника {src} "
                        f"на {det} требует ≥{min_t}s набора "
                        f"для надёжной идентификации (задано "
                        f"{acquisition_time_s:.0f}s). Возможен silent "
                        f"failure калибровки (пик не будет найден "
                        f"Mariscotti)."
                    ),
                    source_label=src,
                    line_keV=E_f,
                    detector_class=det,
                ))

    return warnings


__all__ = [
    "SourceSetWarning",
    "validate_source_set_for_detector",
]
