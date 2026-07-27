"""F-120/F-125 (v1.17.6) — smoke-тест K-40 fixture."""
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
    "archive/K40_420-7-20_Маринелли_0cm.spe"
)
_BACKGROUND = (
    "detectors/Gamma-1S/data/averaged_backgrounds/"
    "bg_2016_marinelli_water_marinelli.spe"
)


def test_k40_pipeline_runs():
    """analyze_lsrm_spe не падает на K-40 fixture."""
    if not os.path.exists(_FIXTURE):
        pytest.skip("K-40 fixture not available")
    r = analyze_lsrm_spe(
        _FIXTURE, detector_type="NaI",
        background_path=_BACKGROUND if os.path.exists(_BACKGROUND) else None,
        allow_stage2=True, allow_stage3=True,
        compute_activities=True, sample_mass_kg=0.5,
        apply_deconvolution=True,
    )
    assert r is not None


def test_k40_detected():
    """K-40 в final_detected."""
    if not os.path.exists(_FIXTURE):
        pytest.skip("K-40 fixture not available")
    r = analyze_lsrm_spe(
        _FIXTURE, detector_type="NaI",
        background_path=_BACKGROUND if os.path.exists(_BACKGROUND) else None,
        allow_stage2=True,
    )
    names = {n.nuclide for n in r.final_detected}
    assert "K-40" in names, f"K-40 не найден; final={names}"


def test_k40_no_th_or_u_dominance():
    """K-40 fixture не запускает Th/U-dominance."""
    if not os.path.exists(_FIXTURE):
        pytest.skip("K-40 fixture not available")
    r = analyze_lsrm_spe(
        _FIXTURE, detector_type="NaI",
        background_path=_BACKGROUND if os.path.exists(_BACKGROUND) else None,
    )
    cd = r.chain_dominance
    if cd is None:
        return
    # K-40 не должен запускать Th-232 dominance (1461 — близка к Ac-228 1459
    # но без других Th-проксей не должен сработать)
    assert not cd.th232, "Th-232 dominance ложно сработал на K-40 fixture"
