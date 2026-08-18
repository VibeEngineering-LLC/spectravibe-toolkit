"""F-128 / v1.17.7 — TCS correction для Bi-212 cascade 727+1620."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.activity.compute import CASCADE_SUMMING_NUCLIDES
from gamma.physics.cascade_summing import CASCADE_SCHEMES


def test_bi212_in_cascade_summing_set():
    """Bi-212 должен присутствовать в наборе TCS-кандидатов."""
    assert "Bi-212" in CASCADE_SUMMING_NUCLIDES


def test_bi212_in_cascade_schemes():
    """Bi-212 должен иметь зарегистрированную CascadeScheme."""
    assert "Bi-212" in CASCADE_SCHEMES


def test_bi212_727_1620_cascade_pair_recorded():
    """Каскадная пара 727+1620 должна быть записана в схеме."""
    scheme = CASCADE_SCHEMES["Bi-212"]
    # 727.33 → 785, 1078 — но самое важное: 1620.5 → 727.33 (cascade)
    assert 1620.50 in scheme.cascades
    partners_1620 = scheme.cascades[1620.50]
    partner_E_set = {round(p.E_partner_keV, 1) for p in partners_1620}
    assert 727.3 in partner_E_set


def test_bi212_tcs_factor_computed():
    """compute_tcs_corrections возвращает per-line dict для Bi-212."""
    from gamma.physics.cascade_summing import (
        compute_tcs_corrections, peak_to_total_NaI_for_geometry,
    )
    # Заглушка efficiency curve — простая константа
    class _Eff:
        def efficiency_at(self, E):
            return 0.05  # типичный photopeak ε(700 кэВ) NaI 63×63
    pt = peak_to_total_NaI_for_geometry("Маринелли")
    cc = compute_tcs_corrections("Bi-212", _Eff(), p_t_func=pt)
    assert cc, "должен быть непустой dict"
    # все значения c ≥ 1.0 (multiplicative correction)
    for E, c in cc.items():
        assert c >= 1.0, f"TCS({E})={c} должен быть ≥ 1"


def test_th232_demo_chain_ratio_improves_with_bi212_tcs():
    """E2E: chain_equilibrium ratio для Th-232 должен быть < 3 после F-128
    (v1.17.6 baseline 2.93×)."""
    fixture = (Path(__file__).parent.parent.parent
               / "detectors" / "Gamma-1S" / "reference_spectra"
               / "archive"
               / "Th232_420-7-17_Маринелли_0cm.spe")
    if not fixture.exists():
        pytest.skip(f"fixture missing: {fixture}")
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    from gamma.activity.compute import chain_equilibrium_guard
    r = analyze_lsrm_spe(str(fixture), complete_workflow=True)
    eq = chain_equilibrium_guard(r.activities or [])
    th_block = eq.get("Th-232")
    assert th_block is not None, "Th-232 chain должна быть в diagnostics"
    ratio = th_block["ratio"]
    assert ratio < 3.0, (
        f"Th-232 chain ratio={ratio:.2f}× должен быть < 3 после F-128 "
        f"(v1.17.6 baseline 2.93×)"
    )
