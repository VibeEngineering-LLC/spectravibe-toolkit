"""Unit tests for ``scripts/gamma/physics/tcs_correction.py`` (Phase 5).

Test plan (v1.19.0 Phase 5 — PLAN §5 step 3):

* 8 base tests — fc_lookup hits/misses/tolerance + apply_tcs_correction
  arithmetic + 2 STUB DeprecationWarning checks.
* 4 detector-class gating tests — fc_lookup detector_class filter,
  apply_tcs_correction_gated out-of-scope behavior (NaI-Tl), the same
  parametrized over remaining unsupported classes, and the happy path
  for HPGe-coaxial.
* 2 inner physics-refinement tests — _compute_tscf returns separate
  tscf_out / tscf_in keys (Refinement A), isotropic-W default
  matches hand-computed reference within 1e-6 rel tol (Refinement C).
* 2 D1/D3 fail-soft tests — low detector_class_CI, and the
  parametrized known-but-unsupported case extending gating test 3.
* 1 D3 fail-closed test — unknown detector_class raises ValueError
  with the expected message (combines "NaI" typo + None into one
  test for coverage; matches PLAN §5 test 16 description).
* 1 D2 status-warnings parity test — SYNTHETIC json_report dict; we
  do NOT touch production reporting code (Risk R5.1 forbids).
* 1 dict-shape sanity test — _compute_tscf returns exactly the four
  expected keys.

Total: 19 named tests; ``test_apply_gated_unsupported_parametrized``
and ``test_known_but_unsupported_fails_soft_parametrized`` are
``@pytest.mark.parametrize``'d so pytest will report >19 individual
test items.

This file uses scaffolds drafted via local Ollama qwen3-coder:30b
(see ``audit/_drafts/_ollama_helpers/gen_tcs_module_scaffolds.py``)
and reviewed/fixed by Claude on 2026-06-03 — bugs corrected:

* Captured ``report_diagnostics`` as a local variable before passing
  (Ollama scaffold inlined the dict literal, losing the reference).
* ``test_summing_in_and_out_signs_opposite`` uses TWO distinct
  ``_compute_tscf`` calls — one for the 1173-keV FEP line (tscf_out
  loss), one for the 2505.7-keV sum-peak entry (tscf_in gain) —
  rather than expecting a single non-sum-peak line to exhibit gain
  (which contradicts the Refinement-A formulas).
* ``test_apply_gated_hpge_coaxial_applied`` renames the unpacked
  value to ``A_corr`` so the assertion ``A_corr > 100.0`` reads
  correctly (Ollama used ``A_obs`` which loses the meaning).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import pytest

from gamma.physics.tcs_correction import (
    ALLOWED_DETECTOR_CLASSES,
    MIN_DETECTOR_CLASS_CI,
    SUPPORTED_DETECTOR_CLASSES,
    _compute_tscf,
    angular_correlation_factor,
    apply_tcs_correction,
    apply_tcs_correction_gated,
    eta_eff_marinelli_scaling,
    fc_lookup,
)
from gamma.data.geometry_classes import GeometryClass


# ---------------------------------------------------------------------------
# Known good baseline: Co-60 1173 keV in petri_100ml_water on HPGe-coaxial.
# Source: data/tsc_lookup.json (Giubrone 2016 Table 4, PGAQ).
# ---------------------------------------------------------------------------
EXPECTED_CO60_FC: Tuple[float, float] = (1.078, 0.012)


# ---------------------------------------------------------------------------
# Test doubles for _compute_tscf (Refinements A & C).
# ---------------------------------------------------------------------------


class _DetectorResponseStub:
    """Detector with hard-coded ε_T / ε_pp lookup dictionaries.

    Realistic-order HPGe-coaxial values used in tests:

    * ε_T(1332 keV) ≈ 0.01  — total efficiency at the 1332-keV partner
      energy (Co-60 cascade).
    * ε_pp(1173 keV) ≈ 0.005, ε_pp(1332 keV) ≈ 0.0045,
      ε_pp(2505.7 keV) ≈ 0.0025 — FEP efficiencies on the cascade
      lines and the 1173+1332 sum-peak.
    """

    def __init__(self, eps_T_dict: dict, eps_pp_dict: dict) -> None:
        self._eps_T = dict(eps_T_dict)
        self._eps_pp = dict(eps_pp_dict)

    def eps_T(self, E: float) -> float:
        return self._eps_T.get(E, 0.0)

    def eps_pp(self, E: float) -> float:
        return self._eps_pp.get(E, 0.0)


@dataclass
class _PartnerStub:
    E_keV: float
    P_cascade: float
    level_j: int = 0
    level: int = 0  # for contributing_pairs use too


@dataclass
class _LineStub:
    E_keV: float
    coincident_partners: List[_PartnerStub] = field(default_factory=list)
    is_sum_peak_candidate: bool = False
    contributing_pairs: List[Tuple[_PartnerStub, _PartnerStub]] = field(default_factory=list)
    level_i: int = 0


# ===========================================================================
# 1. test_fc_lookup_hit_co60_1173
# ===========================================================================
def test_fc_lookup_hit_co60_1173() -> None:
    """Exact-energy hit returns the known Giubrone 2016 Table 4 value."""
    result = fc_lookup(
        nuclide="Co-60",
        line_E_keV=1173.0,
        geometry_class="petri_100ml_water",
        detector_class="HPGe-coaxial",
        tolerance_keV=0.5,
    )
    assert result == EXPECTED_CO60_FC


# ===========================================================================
# 2. test_fc_lookup_miss_unknown_nuclide
# ===========================================================================
def test_fc_lookup_miss_unknown_nuclide() -> None:
    """Nuclide not present in lookup table → None."""
    result = fc_lookup(
        nuclide="Xe-999",
        line_E_keV=1173.0,
        geometry_class="petri_100ml_water",
        detector_class="HPGe-coaxial",
        tolerance_keV=0.5,
    )
    assert result is None


# ===========================================================================
# 3. test_fc_lookup_miss_unknown_geometry
# ===========================================================================
def test_fc_lookup_miss_unknown_geometry() -> None:
    """Geometry not present in lookup table → None."""
    result = fc_lookup(
        nuclide="Co-60",
        line_E_keV=1173.0,
        geometry_class="moon_rock",
        detector_class="HPGe-coaxial",
        tolerance_keV=0.5,
    )
    assert result is None


# ===========================================================================
# 4. test_fc_lookup_within_tolerance
# ===========================================================================
def test_fc_lookup_within_tolerance() -> None:
    """1172.5 keV with tolerance=0.5 keV hits the 1173-keV entry."""
    result = fc_lookup(
        nuclide="Co-60",
        line_E_keV=1172.5,
        geometry_class="petri_100ml_water",
        detector_class="HPGe-coaxial",
        tolerance_keV=0.5,
    )
    assert result == EXPECTED_CO60_FC


# ===========================================================================
# 5. test_fc_lookup_outside_tolerance
# ===========================================================================
def test_fc_lookup_outside_tolerance() -> None:
    """1170 keV is outside tolerance=0.5 keV from 1173 → None."""
    result = fc_lookup(
        nuclide="Co-60",
        line_E_keV=1170.0,
        geometry_class="petri_100ml_water",
        detector_class="HPGe-coaxial",
        tolerance_keV=0.5,
    )
    assert result is None


# ===========================================================================
# 6. test_apply_tcs_correction_arithmetic
# ===========================================================================
def test_apply_tcs_correction_arithmetic() -> None:
    """A_corr = A_obs × fc; σ via GUM with independent terms.

    Hand-computed:
      A_corr = 100.0 × 1.10 = 110.0
      σ_corr = sqrt((1.10 × 5.0)² + (100.0 × 0.01)²)
             = sqrt(30.25 + 1.00) = sqrt(31.25) ≈ 5.5901699...
    """
    A_corr, sigma_corr = apply_tcs_correction(
        A_observed=100.0,
        A_uncertainty=5.0,
        fc=1.10,
        sigma_fc=0.01,
    )
    assert math.isclose(A_corr, 110.0, rel_tol=1e-12)
    assert math.isclose(sigma_corr, math.sqrt(31.25), rel_tol=1e-9)


# ===========================================================================
# 7. test_eta_eff_marinelli_scaling_stub
# ===========================================================================
def test_eta_eff_marinelli_scaling_stub() -> None:
    """STUB returns 1.0 and emits DeprecationWarning."""
    with pytest.warns(DeprecationWarning, match="eta_eff scaling not implemented"):
        result = eta_eff_marinelli_scaling(1000.0, 500.0, "marinelli_500ml_water")
    assert result == 1.0


# ===========================================================================
# 8. test_angular_correlation_factor_stub
# ===========================================================================
def test_angular_correlation_factor_stub() -> None:
    """STUB returns 1.0 and emits DeprecationWarning."""
    with pytest.warns(DeprecationWarning, match="angular correlation factor not implemented"):
        result = angular_correlation_factor("marinelli_500ml_water", 0.30)
    assert result == 1.0


# ===========================================================================
# 9. test_fc_lookup_naitl_returns_none
# ===========================================================================
def test_fc_lookup_naitl_returns_none() -> None:
    """detector_class filter is applied BEFORE energy tolerance.

    Even though Co-60 1173 keV PGAQ has an HPGe-coaxial entry in the
    table, requesting NaI-Tl class returns None — cross-class
    extrapolation is forbidden (R5.4).
    """
    result = fc_lookup(
        nuclide="Co-60",
        line_E_keV=1173.0,
        geometry_class="petri_100ml_water",
        detector_class="NaI-Tl",
        tolerance_keV=0.5,
    )
    assert result is None


# ===========================================================================
# 10. test_apply_gated_naitl_out_of_scope
# ===========================================================================
def test_apply_gated_naitl_out_of_scope() -> None:
    """NaI-Tl with high CI → out-of-scope-for-detector, A unchanged, warning logged."""
    report_diagnostics: dict = {}
    A_out, sigma_out, status = apply_tcs_correction_gated(
        A_observed=100.0,
        A_uncertainty=5.0,
        nuclide="Co-60",
        line_E_keV=1173.0,
        geometry_class="petri_100ml_water",
        detector_class="NaI-Tl",
        detector_class_CI=0.95,
        report_diagnostics=report_diagnostics,
    )
    assert status == "out-of-scope-for-detector"
    assert A_out == 100.0
    assert sigma_out == 5.0
    assert len(report_diagnostics["tcs_warnings"]) == 1
    entry = report_diagnostics["tcs_warnings"][0]
    assert entry["nuclide"] == "Co-60"
    assert entry["detector_class"] == "NaI-Tl"
    assert entry["status"] == "out-of-scope-for-detector"
    assert "NaI-Tl" in entry["message"]


# ===========================================================================
# 11. test_apply_gated_unsupported_parametrized (D3 fail-soft fan-out)
# ===========================================================================
@pytest.mark.parametrize("detector_class", ["LaBr3-Ce", "CeBr3", "CdZnTe", "HPGe-BEGE"])
def test_apply_gated_unsupported_parametrized(detector_class: str) -> None:
    """All known-but-unsupported classes fail soft with out-of-scope status."""
    report_diagnostics: dict = {}
    A_out, sigma_out, status = apply_tcs_correction_gated(
        A_observed=100.0,
        A_uncertainty=5.0,
        nuclide="Co-60",
        line_E_keV=1173.0,
        geometry_class="petri_100ml_water",
        detector_class=detector_class,
        detector_class_CI=0.95,
        report_diagnostics=report_diagnostics,
    )
    assert status == "out-of-scope-for-detector"
    assert A_out == 100.0
    assert sigma_out == 5.0
    assert len(report_diagnostics["tcs_warnings"]) == 1
    assert report_diagnostics["tcs_warnings"][0]["detector_class"] == detector_class


# ===========================================================================
# 12. test_apply_gated_hpge_coaxial_applied
# ===========================================================================
def test_apply_gated_hpge_coaxial_applied() -> None:
    """HPGe-coaxial + high CI + Co-60 1173 in PGAQ → status='applied', A_corr > A_obs."""
    report_diagnostics: dict = {}
    A_corr, sigma_corr, status = apply_tcs_correction_gated(
        A_observed=100.0,
        A_uncertainty=5.0,
        nuclide="Co-60",
        line_E_keV=1173.0,
        geometry_class="petri_100ml_water",
        detector_class="HPGe-coaxial",
        detector_class_CI=0.95,
        report_diagnostics=report_diagnostics,
    )
    assert status == "applied"
    # fc = 1.078 → A_corr = 107.8
    assert math.isclose(A_corr, 100.0 * 1.078, rel_tol=1e-12)
    # σ_corr = sqrt((1.078 × 5)² + (100 × 0.012)²) ≈ sqrt(29.052 + 1.44) ≈ 5.521
    expected_sigma = math.sqrt((1.078 * 5.0) ** 2 + (100.0 * 0.012) ** 2)
    assert math.isclose(sigma_corr, expected_sigma, rel_tol=1e-9)
    assert report_diagnostics.get("tcs_warnings", []) == []


# ===========================================================================
# 13. test_summing_in_and_out_signs_opposite (Refinement A)
# ===========================================================================
def test_summing_in_and_out_signs_opposite() -> None:
    """Refinement A: summing-out (loss) and summing-in (sum-peak gain)
    are separate, physically-opposite contributions.

    Two distinct calls:

    * A 1173-keV Co-60 FEP line with one coincident partner at 1332
      keV → ``tscf_out > 1.0`` (loss correction); no sum-peak
      channel → ``tscf_in == 1.0``.
    * A 2505.7-keV sum-peak entry with contributing pair
      (1173, 1332) and ``is_sum_peak_candidate=True`` → ``tscf_in
      > 1.0`` separately (the sum-peak channel adds counts).

    Both keys are present in the returned dict — not collapsed.
    """
    # Stub detector with realistic HPGe-coaxial ε_T / ε_pp at Co-60 energies.
    eps_T = {1173.0: 0.012, 1332.0: 0.010}
    eps_pp = {1173.0: 0.005, 1332.0: 0.0045, 2505.7: 0.0025}
    detector = _DetectorResponseStub(eps_T, eps_pp)

    # --- FEP line at 1173 keV: partner at 1332 keV, no sum-peak candidate.
    fep_partner = _PartnerStub(E_keV=1332.0, P_cascade=1.0, level_j=1)
    fep_line = _LineStub(
        E_keV=1173.0,
        coincident_partners=[fep_partner],
        is_sum_peak_candidate=False,
        contributing_pairs=[],
    )
    res_fep = _compute_tscf(fep_line, detector)
    assert set(res_fep.keys()) == {"tscf_out", "tscf_in", "tscf_total", "sigma"}
    assert res_fep["tscf_out"] > 1.0  # FEP loss → tscf_out > 1
    assert res_fep["tscf_in"] == 1.0  # No sum-peak channel for this line
    # Sanity: tscf_out = 1 / (1 - eps_T(1332) × P × W) = 1 / (1 - 0.010)
    assert math.isclose(res_fep["tscf_out"], 1.0 / 0.99, rel_tol=1e-9)

    # --- Sum-peak entry at 2505.7 keV: contributing pair (1173, 1332).
    sum_i = _PartnerStub(E_keV=1173.0, P_cascade=1.0, level=1)
    sum_j = _PartnerStub(E_keV=1332.0, P_cascade=1.0, level=0)
    sum_line = _LineStub(
        E_keV=2505.7,
        coincident_partners=[],  # No FEP-loss channel for the sum-peak itself
        is_sum_peak_candidate=True,
        contributing_pairs=[(sum_i, sum_j)],
    )
    res_sum = _compute_tscf(sum_line, detector)
    assert set(res_sum.keys()) == {"tscf_out", "tscf_in", "tscf_total", "sigma"}
    # No coincident partners listed for the sum-peak entry → no loss → tscf_out = 1.
    assert res_sum["tscf_out"] == 1.0
    # Summing-in: gain = eps_pp(1173) × eps_pp(1332) × P × W / eps_pp(2505.7)
    #            = 0.005 × 0.0045 × 1.0 × 1.0 / 0.0025 = 0.009
    expected_gain = 0.005 * 0.0045 * 1.0 * 1.0 / 0.0025
    assert math.isclose(res_sum["tscf_in"], 1.0 + expected_gain, rel_tol=1e-9)
    assert res_sum["tscf_in"] > 1.0


# ===========================================================================
# 14. test_isotropic_default_W_equals_1 (Refinement C)
# ===========================================================================
def test_isotropic_default_W_equals_1() -> None:
    """Default angular_correlation_provider=None → W ≡ 1.

    Hand-computed reference for Co-60 1173 line with one partner at
    1332 keV, ε_T(1332)=0.01, P_cascade=1.0:

        tscf_out = 1 / (1 − 0.01 × 1.0 × 1.0) = 1 / 0.99 = 1.0101010101...
        tscf_in = 1.0 (no sum-peak channel)
        tscf_total = tscf_out × tscf_in = 1.0101010101...

    Must match to within 1e-6 relative tolerance.
    """
    eps_T = {1332.0: 0.010}
    eps_pp = {1173.0: 0.005, 1332.0: 0.0045}
    detector = _DetectorResponseStub(eps_T, eps_pp)

    partner = _PartnerStub(E_keV=1332.0, P_cascade=1.0)
    line = _LineStub(
        E_keV=1173.0,
        coincident_partners=[partner],
        is_sum_peak_candidate=False,
    )
    result = _compute_tscf(line, detector, angular_correlation_provider=None)

    expected_tscf_out = 1.0 / (1.0 - 0.01 * 1.0 * 1.0)
    assert math.isclose(result["tscf_out"], expected_tscf_out, rel_tol=1e-6)
    assert result["tscf_in"] == 1.0
    assert math.isclose(result["tscf_total"], expected_tscf_out, rel_tol=1e-6)


# ===========================================================================
# 15. test_low_detector_CI_fails_soft (D1)
# ===========================================================================
def test_low_detector_CI_fails_soft() -> None:
    """detector_class_CI=0.5 < 0.7 → low-detector-CI, A unchanged,
    warning message contains 'CI=0.50' and the threshold '0.7'.
    """
    report_diagnostics: dict = {}
    A_out, sigma_out, status = apply_tcs_correction_gated(
        A_observed=100.0,
        A_uncertainty=5.0,
        nuclide="Co-60",
        line_E_keV=1173.0,
        geometry_class="petri_100ml_water",
        detector_class="HPGe-coaxial",
        detector_class_CI=0.5,
        report_diagnostics=report_diagnostics,
    )
    assert status == "low-detector-CI"
    assert A_out == 100.0
    assert sigma_out == 5.0
    assert len(report_diagnostics["tcs_warnings"]) == 1
    entry = report_diagnostics["tcs_warnings"][0]
    assert entry["status"] == "low-detector-CI"
    assert entry["detector_class_CI"] == 0.5
    assert "CI=0.50" in entry["message"]
    assert "0.7" in entry["message"]


# ===========================================================================
# 16. test_unknown_detector_class_raises_ValueError (D3 fail-closed)
# ===========================================================================
def test_unknown_detector_class_raises_ValueError() -> None:
    """Unknown enum member (typo, None) → ValueError with
    'Unknown detector_class' and sorted(ALLOWED_DETECTOR_CLASSES).
    """
    # Typo case: "NaI" (missing -Tl).
    with pytest.raises(ValueError, match="Unknown detector_class"):
        apply_tcs_correction_gated(
            A_observed=100.0,
            A_uncertainty=5.0,
            nuclide="Co-60",
            line_E_keV=1173.0,
            geometry_class="petri_100ml_water",
            detector_class="NaI",  # typo
            detector_class_CI=0.95,
            report_diagnostics={},
        )

    # None case.
    with pytest.raises(ValueError) as excinfo:
        apply_tcs_correction_gated(
            A_observed=100.0,
            A_uncertainty=5.0,
            nuclide="Co-60",
            line_E_keV=1173.0,
            geometry_class="petri_100ml_water",
            detector_class=None,
            detector_class_CI=0.95,
            report_diagnostics={},
        )
    msg = str(excinfo.value)
    assert "Unknown detector_class" in msg
    # The full sorted list of allowed classes must be referenced so the
    # caller can self-diagnose the typo. Spot-check a couple of names.
    assert "HPGe-coaxial" in msg
    assert "NaI-Tl" in msg


# ===========================================================================
# 17. test_known_but_unsupported_fails_soft_parametrized (D3 fail-soft)
# ===========================================================================
@pytest.mark.parametrize(
    "detector_class",
    ["NaI-Tl", "LaBr3-Ce", "CeBr3", "CdZnTe", "HPGe-BEGE"],
)
def test_known_but_unsupported_fails_soft_parametrized(detector_class: str) -> None:
    """Every known-but-unsupported class → out-of-scope-for-detector,
    warning appended; high CI (0.95) so D1 branch does NOT fire.

    Extends ``test_apply_gated_unsupported_parametrized`` by also
    including ``NaI-Tl`` in the parametrize list.
    """
    report_diagnostics: dict = {}
    A_out, sigma_out, status = apply_tcs_correction_gated(
        A_observed=100.0,
        A_uncertainty=5.0,
        nuclide="Co-60",
        line_E_keV=1173.0,
        geometry_class="petri_100ml_water",
        detector_class=detector_class,
        detector_class_CI=0.95,
        report_diagnostics=report_diagnostics,
    )
    assert status == "out-of-scope-for-detector"
    assert (A_out, sigma_out) == (100.0, 5.0)
    assert len(report_diagnostics["tcs_warnings"]) == 1
    assert report_diagnostics["tcs_warnings"][0]["status"] == "out-of-scope-for-detector"


# ===========================================================================
# 18. test_tcs_status_column_matches_warnings_list (D2 parity)
# ===========================================================================
def test_tcs_status_column_matches_warnings_list() -> None:
    """D2 parity: every non-'applied' row has a matching entry in
    tcs_warnings.

    SYNTHETIC json_report dict — we do NOT modify production
    reporting code in Phase 5 (Risk R5.1). This test only locks the
    invariant the v1.19.1 wiring must satisfy.
    """
    json_report = {
        "fep_table": [
            {"nuclide": "Co-60", "line_E_keV": 1173.0, "tcs_status": "applied"},
            {"nuclide": "Co-60", "line_E_keV": 1332.0, "tcs_status": "applied"},
            {
                "nuclide": "Cs-134",
                "line_E_keV": 604.7,
                "tcs_status": "out-of-scope-for-detector",
            },
            {
                "nuclide": "Cs-134",
                "line_E_keV": 795.9,
                "tcs_status": "low-detector-CI",
            },
        ],
        "tcs_warnings": [
            {
                "nuclide": "Cs-134",
                "line_E_keV": 604.7,
                "status": "out-of-scope-for-detector",
                "message": "...",
            },
            {
                "nuclide": "Cs-134",
                "line_E_keV": 795.9,
                "status": "low-detector-CI",
                "message": "...",
            },
        ],
    }

    fep_statuses = {row["tcs_status"] for row in json_report["fep_table"]}
    warning_statuses = {w["status"] for w in json_report["tcs_warnings"]}
    assert fep_statuses <= ({"applied"} | warning_statuses), (
        f"FEP status column has values not present in warnings: "
        f"{fep_statuses - {'applied'} - warning_statuses}"
    )

    # Per-row check: every non-'applied' row has at least one matching
    # warnings entry on (nuclide, line_E_keV).
    for row in json_report["fep_table"]:
        if row["tcs_status"] == "applied":
            continue
        matches = [
            w
            for w in json_report["tcs_warnings"]
            if w["nuclide"] == row["nuclide"]
            and w["line_E_keV"] == row["line_E_keV"]
            and w["status"] == row["tcs_status"]
        ]
        assert matches, f"No tcs_warnings entry matches FEP row {row!r}"


# ===========================================================================
# 19. test_compute_tscf_returns_dict_with_required_keys
# ===========================================================================
def test_compute_tscf_returns_dict_with_required_keys() -> None:
    """_compute_tscf returns exactly the four expected keys, all floats."""
    eps_T = {1332.0: 0.010}
    eps_pp = {1173.0: 0.005, 1332.0: 0.0045}
    detector = _DetectorResponseStub(eps_T, eps_pp)

    partner = _PartnerStub(E_keV=1332.0, P_cascade=1.0)
    line = _LineStub(
        E_keV=1173.0,
        coincident_partners=[partner],
        is_sum_peak_candidate=False,
    )
    result = _compute_tscf(line, detector)

    assert set(result.keys()) == {"tscf_out", "tscf_in", "tscf_total", "sigma"}
    for key in ("tscf_out", "tscf_in", "tscf_total", "sigma"):
        assert isinstance(result[key], float)
