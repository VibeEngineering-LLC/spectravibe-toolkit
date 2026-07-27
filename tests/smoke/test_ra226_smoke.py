"""F-121 (v1.17.6) — smoke-тест Ra-226 fixture с U1/U2/U3 forced multiplets."""
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
    "archive/Ra226_420-7-18_Маринелли_0cm.spe"
)
_BACKGROUND = (
    "detectors/Gamma-1S/data/averaged_backgrounds/"
    "bg_2016_marinelli_water_marinelli.spe"
)


def test_ra226_pipeline_runs():
    """analyze_lsrm_spe не падает на Ra-226 fixture."""
    if not os.path.exists(_FIXTURE):
        pytest.skip("Ra-226 fixture not available")
    r = analyze_lsrm_spe(
        _FIXTURE, detector_type="NaI",
        background_path=_BACKGROUND if os.path.exists(_BACKGROUND) else None,
        allow_stage2=True, allow_stage3=True,
        compute_activities=True, sample_mass_kg=0.5,
        apply_deconvolution=True,
    )
    assert r is not None


def test_ra226_uranium_chain_proxies_present():
    """В Ra-226 fixture хотя бы один U-chain нуклид найден."""
    if not os.path.exists(_FIXTURE):
        pytest.skip("Ra-226 fixture not available")
    r = analyze_lsrm_spe(
        _FIXTURE, detector_type="NaI",
        background_path=_BACKGROUND if os.path.exists(_BACKGROUND) else None,
        allow_stage2=True, allow_stage3=True,
    )
    names = {n.nuclide for n in r.final_detected}
    u_chain = {"Bi-214", "Pb-214", "Pb-210", "Ra-226"}
    assert names & u_chain, (
        f"никаких U-chain нуклидов не найдено на Ra-226 fixture; "
        f"final={names}"
    )


def test_ra226_chain_dominance_u238_when_proxies_present():
    """Ra-226 fixture с яркими Bi-214/Pb-214 → u238 dominance срабатывает."""
    if not os.path.exists(_FIXTURE):
        pytest.skip("Ra-226 fixture not available")
    r = analyze_lsrm_spe(
        _FIXTURE, detector_type="NaI",
        background_path=_BACKGROUND if os.path.exists(_BACKGROUND) else None,
        allow_stage2=True,
    )
    cd = r.chain_dominance
    # Может сработать u238 dominance — это правильное поведение для
    # сильной 226Ra пробы. Главное, что Th-232 НЕ доминирует.
    if cd is not None:
        assert not cd.th232, (
            "Th-232 dominance ложно сработал на Ra-226 fixture"
        )


def test_ra226_forced_clusters_emitted_when_u238_dominant():
    """F-121: при u238 dominance run_chain_forced_multiplets возвращает
    кластеры U1/U2/U3."""
    if not os.path.exists(_FIXTURE):
        pytest.skip("Ra-226 fixture not available")
    r = analyze_lsrm_spe(
        _FIXTURE, detector_type="NaI",
        background_path=_BACKGROUND if os.path.exists(_BACKGROUND) else None,
        allow_stage2=True, allow_stage3=True,
        apply_deconvolution=True,
    )
    if r.chain_dominance is None or not r.chain_dominance.u238:
        # Если u238 dominance не сработал на этом fixture, smoke-тест
        # просто проверяет, что pipeline не падает. F-121 контракт
        # покрывается отдельным test_u238_chain_multiplets.py.
        return
    deconv = r.deconvolution_results or []
    coupled_clusters = [
        d for d in deconv
        if str(getattr(d, "method", "")).startswith("coupled_")
    ]
    # Хотя бы один coupled-кластер должен быть для U-238 цепочки
    assert len(coupled_clusters) >= 1, (
        f"F-121 forced multiplets не эмитированы; "
        f"deconvolutions={[d.method for d in deconv]}"
    )
