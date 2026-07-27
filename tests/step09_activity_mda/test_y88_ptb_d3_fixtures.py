"""
PTB-2018 Annex D.3 worked-example fixtures for Y-88 TCS correction (#PTB-4).

Source: "gamma_spekt_grundl" ISSN 1865-8725, Version March 2018,
Annex D.3 "Worked example for Y-88", pages gamma-SPEKT/GRUNDL-63..66
(PDF pages 67-70; read visually from 150-DPI renders, 2026-07-02).

Reference input data (Tab. D2 / Tab. D3, page -63):
    transition 2->1: E = 898.04 keV,  eps_t(E) = 0.094644
    transition 1->0: E = 1836.07 keV, eps_t(E) = 0.075176
    P_q(898) = 0.937, P_q(1836) = 0.9938

Reference results:
    f_898  = 1 / (1 - 0.075176 * 0.99984)          = 1.0813   (page -64)
    f_1836 = 1 / (1 - (0.937/0.9938) * 0.094644)   = 1.098    (page -66)

Project model (cascade_summing.py, Y-88 scheme):
    C(898)  = 1 / (1 - 1.000 * eps_T(1836))   -> 1.08129 (delta 1.3e-5 vs PTB)
    C(1836) = 1 / (1 - 0.94  * eps_T(898))    -> 1.09765 (delta 3.5e-4 vs PTB)

The residuals come from the project's rounded cascade probabilities
(1.000 vs PTB 0.99984; 0.94 vs PTB 0.937/0.9938 = 0.94285) and stay
well inside the tolerances asserted below.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.physics.cascade_summing import (
    compute_tcs_corrections,
    tcs_correction_factor,
)


# PTB Tab. D3 total efficiencies (page gamma-SPEKT/GRUNDL-63)
EPS_T_898 = 0.094644
EPS_T_1836 = 0.075176

# PTB Tab. D2 emission probabilities per level feed (page -63)
PQ_898 = 0.937
PQ_1836 = 0.9938

# PTB reference correction factors (pages -64 / -66)
F_898_PTB = 1.0813
F_1836_PTB = 1.098


class PTBTotalEfficiency:
    """Stub feeding PTB Tab. D3 eps_t values straight through.

    Used with a unit peak-to-total function so that
    eps_T = eps_p / P = eps_p equals the PTB table value exactly.
    """

    TABLE = {898.04: EPS_T_898, 1836.07: EPS_T_1836}

    def efficiency_at(self, E_keV: float) -> float:
        for E_ref, eps in self.TABLE.items():
            if abs(E_ref - E_keV) < 2.0:
                return eps
        return 0.0


def unit_pt(E_keV: float) -> float:
    return 1.0


def test_f898_matches_ptb_reference():
    """C(898) reproduces PTB f_898 = 1.0813 (Annex D.3.1, page -64)."""
    eff = PTBTotalEfficiency()
    c = tcs_correction_factor("Y-88", 898.04, eff, p_t_func=unit_pt)
    # Exact project-model value: partner (1836.06, p=1.000)
    assert math.isclose(c, 1.0 / (1.0 - 1.000 * EPS_T_1836), rel_tol=1e-12)
    # PTB reference (their p = 1/(1+alpha) = 0.99984; delta 1.3e-5)
    assert abs(c - F_898_PTB) < 5e-4, (
        f"f_898: project {c:.6f} vs PTB {F_898_PTB} "
        f"(delta {abs(c - F_898_PTB):.2e} > 5e-4)"
    )
    print(f"  ok test_f898_matches_ptb_reference (C={c:.5f} vs 1.0813)")


def test_f1836_matches_ptb_reference():
    """C(1836) reproduces PTB f_1836 = 1.098 (Annex D.3.2, page -66)."""
    eff = PTBTotalEfficiency()
    c = tcs_correction_factor("Y-88", 1836.06, eff, p_t_func=unit_pt)
    # Exact project-model value: partner (898.04, p=0.94)
    assert math.isclose(c, 1.0 / (1.0 - 0.94 * EPS_T_898), rel_tol=1e-12)
    # PTB reference uses p = P_q(898)/P_q(1836) = 0.94285; delta 3.5e-4
    assert abs(c - F_1836_PTB) < 1e-3, (
        f"f_1836: project {c:.6f} vs PTB {F_1836_PTB} "
        f"(delta {abs(c - F_1836_PTB):.2e} > 1e-3)"
    )
    print(f"  ok test_f1836_matches_ptb_reference (C={c:.5f} vs 1.098)")


def test_f1836_ptb_exact_probability_formula():
    """Transcription guard: PTB's own expression on page -66 evaluates
    to their printed 1.098 -- 1/(1 - (0.937/0.9938)*0.094644)."""
    f = 1.0 / (1.0 - (PQ_898 / PQ_1836) * EPS_T_898)
    assert abs(f - F_1836_PTB) < 1e-4, f"PTB self-check failed: {f:.6f}"
    print(f"  ok test_f1836_ptb_exact_probability_formula ({f:.5f})")


def test_y88_corrections_use_partner_energy_not_own():
    """C(E_i) must depend on eps_t at the PARTNER energy: zeroing the
    898-keV efficiency leaves C(898) unchanged but drops C(1836) to 1."""
    class Only1836(PTBTotalEfficiency):
        TABLE = {1836.07: EPS_T_1836}

    eff = Only1836()
    c898 = tcs_correction_factor("Y-88", 898.04, eff, p_t_func=unit_pt)
    c1836 = tcs_correction_factor("Y-88", 1836.06, eff, p_t_func=unit_pt)
    assert math.isclose(c898, 1.0 / (1.0 - EPS_T_1836), rel_tol=1e-12)
    assert c1836 == 1.0
    print("  ok test_y88_corrections_use_partner_energy_not_own")


def test_compute_tcs_corrections_y88_bundle():
    """Dispatcher returns both Y-88 lines with the same PTB-anchored
    values as the per-line calls."""
    eff = PTBTotalEfficiency()
    cc = compute_tcs_corrections("Y-88", eff, p_t_func=unit_pt)
    assert set(cc.keys()) == {898.04, 1836.06}
    assert abs(cc[898.04] - F_898_PTB) < 5e-4
    assert abs(cc[1836.06] - F_1836_PTB) < 1e-3
    print("  ok test_compute_tcs_corrections_y88_bundle")


if __name__ == "__main__":
    test_f898_matches_ptb_reference()
    test_f1836_matches_ptb_reference()
    test_f1836_ptb_exact_probability_formula()
    test_y88_corrections_use_partner_energy_not_own()
    test_compute_tcs_corrections_y88_bundle()
    print("All #PTB-4 Y-88 fixture tests passed.")