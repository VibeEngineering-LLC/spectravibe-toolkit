"""gamma.detectors.profile — JSON profile loader + silent-fallback detector.

BUG-40 / BUG-39 (Wave 6, v1.22.0; F2-A renormalisation 2026-06-21):

* **BUG-40**: each spectrometric complex (``Gamma-1S``, future
  ``AtomSpectra``, …) has a metadata profile at
  ``references/detectors/<canonical>.json``. The profile carries FWHM
  polynomial, efficiency-source pointer and provenance.
* **BUG-39**: when a spectrum's canonical detector resolves to a name
  whose profile is missing on disk (e.g. a future hypothetical
  ``AtomSpectra`` complex without ``references/detectors/AtomSpectra.json``),
  the loader returns a :class:`DetectorFallback` record. Pipeline
  propagates it to ``report.json`` ``warnings`` and to the HTML/MD
  reports — the substitution is no longer silent.

The loader is intentionally cheap and side-effect free; algorithms that
need detector paths still go through ``gamma.detectors.gamma1s`` (or
future ``gamma.detectors.<name>``). This module is the **metadata
registry** layer that sits ABOVE those resolvers.

F2-A note (2026-06-21): the legacy Case-2 ``efficiency_fallback_to``
stub branch was removed together with the bogus ``Gamma-1S`` stub
profile that drove it. Only Case 1 (profile missing on disk) and
Case 3 (clean load) remain.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Dict, Optional, Any


__all__ = [
    "DetectorFallback",
    "DetectorProfile",
    "PROFILES_DIR",
    "load_detector_profile",
    "detect_silent_fallback",
]


#: Root of the detector profile registry (BUG-40).
PROFILES_DIR: Path = (
    Path(__file__).resolve().parents[3] / "references" / "detectors"
)


@dataclass
class DetectorProfile:
    """Loaded metadata for one spectrometric complex.

    Attributes:
        canonical: canonical detector name (e.g. ``"Gamma-1S"``).
        kind: free-form string (``"spectrometric_complex"``).
        validation_status: ``"primary"`` (no stub-pending profiles since
            F2-A 2026-06-21; field kept for forward-compat with future
            workflows that introduce additional statuses).
        efficiency_source_kind: ``"directory"`` (sole remaining kind
            since F2-A; the ``"TBD_pending_calibration_data"`` variant
            was removed together with the bogus Gamma-1S stub profile).
        raw: the full parsed JSON dictionary (for callers that need
            fields the dataclass does not surface).
        source_path: absolute path to the JSON profile file.
    """

    canonical: str
    kind: str = "spectrometric_complex"
    validation_status: str = "primary"
    efficiency_source_kind: str = "directory"
    raw: Dict[str, Any] = field(default_factory=dict)
    source_path: Optional[Path] = None


@dataclass
class DetectorFallback:
    """Record describing a silent-fallback event (BUG-39).

    Emitted when the pipeline's resolved canonical detector is not
    backed by a complete profile + physical assets, and the loader
    chose a substitute.

    Attributes:
        requested: canonical detector the pipeline asked for.
        actual: canonical detector actually used (substitute).
        reason: short machine-readable code:
            ``"profile_not_on_disk"`` — no JSON profile for ``requested``;
            ``"profile_loaded_no_fallback"`` — profile loaded cleanly,
            no fallback applied (no warning emitted).
            (The legacy ``"efficiency_tbd_using_fallback_profile"`` code
            was retired in F2-A 2026-06-21 together with the bogus
            Gamma-1S stub profile that drove it.)
        human_ru: operator-facing Russian message suitable for the
            RU markdown report (F-386 EN-leak gate).
        human_en: same message in English for the HTML/JSON ``warnings``
            channel.
        human: bilingual concatenation ``"<EN>. RU: <RU>"`` — preserved
            for backward compatibility with callers that read a single
            ``human`` field; new callers should pick ``human_ru`` or
            ``human_en`` directly.
    """

    requested: str
    actual: str
    reason: str
    human_ru: str
    human_en: str

    @property
    def human(self) -> str:
        if self.human_en and self.human_ru:
            return f"{self.human_en} RU: {self.human_ru}"
        return self.human_en or self.human_ru or ""

    def as_dict(self) -> Dict[str, str]:
        return {
            "requested": self.requested,
            "actual": self.actual,
            "reason": self.reason,
            "human": self.human,
            "human_ru": self.human_ru,
            "human_en": self.human_en,
        }


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------

@lru_cache(maxsize=16)
def _load_profile_file(canonical: str) -> Optional[DetectorProfile]:
    """Read ``references/detectors/<canonical>.json`` if it exists.

    Returns None when the file is absent. Never raises on malformed JSON
    (returns None and the caller emits a fallback record).
    """
    if not canonical:
        return None
    path = PROFILES_DIR / f"{canonical}.json"
    if not path.is_file():
        return None
    try:
        with path.open("r", encoding="utf-8") as f:
            raw = json.load(f)
    except (OSError, json.JSONDecodeError):
        return None

    eff = raw.get("efficiency_source", {}) or {}
    return DetectorProfile(
        canonical=raw.get("canonical", canonical),
        kind=raw.get("kind", "spectrometric_complex"),
        validation_status=raw.get("validation_status", "primary"),
        efficiency_source_kind=eff.get("kind", "directory"),
        raw=raw,
        source_path=path,
    )


def load_detector_profile(canonical: str) -> Optional[DetectorProfile]:
    """Public reader. ``None`` when the profile JSON does not exist."""
    return _load_profile_file(canonical)


# ---------------------------------------------------------------------------
# Silent-fallback detection (BUG-39)
# ---------------------------------------------------------------------------

def detect_silent_fallback(canonical: str) -> DetectorFallback:
    """Return a :class:`DetectorFallback` describing the loader outcome.

    Two cases (since F2-A 2026-06-21 the Case-2 stub-fallback branch
    was removed together with the bogus Gamma-1S stub profile):

    1. ``canonical`` has no JSON profile on disk ⇒ ``reason =
       "profile_not_on_disk"``, ``actual = "Gamma-1S"`` (current
       hard-coded fallback — the only fully-implemented detector
       subtree in the codebase, F-83). Operator-facing warning emitted.

    2. ``canonical`` has a complete profile with usable efficiency
       assets ⇒ ``reason = "profile_loaded_no_fallback"``. **No warning
       should be emitted** by callers; check via
       :func:`should_emit_warning`.
    """
    if not canonical:
        # Defensive: empty canonical → treat as missing profile.
        return DetectorFallback(
            requested="",
            actual="Gamma-1S",
            reason="profile_not_on_disk",
            human_en=(
                "Detector profile not resolved from spectrum header; "
                "pipeline fell back to Gamma-1S defaults."
            ),
            human_ru=(
                "Профиль детектора не распознан по заголовку спектра — "
                "применены параметры Gamma-1S по умолчанию."
            ),
        )

    prof = load_detector_profile(canonical)
    if prof is None:
        # Case 1 — JSON profile missing on disk.
        return DetectorFallback(
            requested=canonical,
            actual="Gamma-1S",
            reason="profile_not_on_disk",
            human_en=(
                f"Detector profile references/detectors/{canonical}.json "
                f"not found on disk — pipeline fell back to Gamma-1S. "
                f"Quantitative results may carry a detector-cert "
                f"mismatch bias."
            ),
            human_ru=(
                f"Профиль детектора {canonical} отсутствует на диске — "
                f"применены параметры Gamma-1S; возможно смещение "
                f"расчётов из-за несоответствия сертификата детектора."
            ),
        )

    # Case 2 — clean load (F2-A 2026-06-21: removed Case-2 stub-fallback
    # branch together with the bogus Gamma-1S stub profile + the
    # efficiency_fallback_to field on DetectorProfile).
    return DetectorFallback(
        requested=canonical,
        actual=canonical,
        reason="profile_loaded_no_fallback",
        human_en="",
        human_ru="",
    )


def should_emit_warning(fallback: DetectorFallback) -> bool:
    """True when the fallback record should surface in `warnings`."""
    return fallback.reason != "profile_loaded_no_fallback"


# ---------------------------------------------------------------------------
# Cache utilities (for tests)
# ---------------------------------------------------------------------------

def clear_cache() -> None:
    """Drop the LRU cache (tests that monkey-patch PROFILES_DIR call this)."""
    _load_profile_file.cache_clear()
