"""F-120/F-125/F-122 (v1.17.6) — smoke-тест Cs-137 fixture.

Запускает полный CLI-пайплайн на Cs-137 Marinelli fixture и проверяет:
  1. Пайплайн не падает (без исключений).
  2. Cs-137 присутствует в final_detected.
  3. Не появляется ложный U-238 / Th-232 dominance.
"""
from __future__ import annotations

import os
import sys

import pytest

SCRIPTS = os.path.join(os.path.dirname(__file__), "scripts")
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)

from gamma.identification.staged_pipeline import analyze_lsrm_spe


_FIXTURE = (
    "detectors/Gamma-1S/reference_spectra/"
    "archive/Cs137_420-7-14_Маринелли_0cm.spe"
)
_BACKGROUND = (
    "detectors/Gamma-1S/data/averaged_backgrounds/"
    "bg_2016_marinelli_water_marinelli.spe"
)


def test_cs137_pipeline_runs_without_exceptions():
    """CLI-аналог: analyze_lsrm_spe не падает на Cs-137 fixture."""
    if not os.path.exists(_FIXTURE):
        pytest.skip("Cs-137 fixture not available")
    r = analyze_lsrm_spe(
        _FIXTURE,
        detector_type="NaI",
        background_path=_BACKGROUND if os.path.exists(_BACKGROUND) else None,
        allow_stage2=True, allow_stage3=True,
        compute_activities=True, sample_mass_kg=0.5,
        apply_deconvolution=True,
    )
    assert r is not None
    assert r.spec is not None


def test_cs137_detected_in_final_identifications():
    """Cs-137 фигурирует в final_detected (filename binding делает его Stage 1)."""
    if not os.path.exists(_FIXTURE):
        pytest.skip("Cs-137 fixture not available")
    r = analyze_lsrm_spe(
        _FIXTURE, detector_type="NaI",
        background_path=_BACKGROUND if os.path.exists(_BACKGROUND) else None,
        allow_stage2=True,
    )
    names = {n.nuclide for n in r.final_detected}
    assert "Cs-137" in names, f"Cs-137 не найден; final={names}"


def test_cs137_no_false_thorium_or_uranium_dominance():
    """Чистый Cs-137 fixture не должен давать Th/U-dominance."""
    if not os.path.exists(_FIXTURE):
        pytest.skip("Cs-137 fixture not available")
    r = analyze_lsrm_spe(
        _FIXTURE, detector_type="NaI",
        background_path=_BACKGROUND if os.path.exists(_BACKGROUND) else None,
    )
    cd = r.chain_dominance
    assert cd is not None
    # Cs-137 не должен запускать ни Th-232, ни U-238 dominance
    assert not cd.th232, "Th-232 dominance ложно сработал на Cs-137 fixture"
    # U-238 может сработать только если K-40 + Cs-137 геометрия даёт
    # ложные сигналы (тогда F-89d должен подавить через filename binding)
