"""BUG-41 / Wave 7 — pathology-detection fallback for ``build_fwhm_model``.

Root cause (revalidation outbox section 5 BUG-41 and brief
``_state/agent_a/inbox/2026-06-05_C2_FWHM_BUG37_BUG41.md``):
  When the embedded ``lsrm_peaks_table`` is sparse, clustered, or
  ill-distributed, the least-squares ``FWHM^2(E) = a + b*E + c*E^2``
  fit produces an unphysical negative discriminant at low E. The
  ``max(val, 0.01)`` numerical floor in ``fwhm_keV_at_energy`` then
  emits FWHM = 0.1 keV (BUG-37 lowered the floor from 1.0 -> 0.01;
  BUG-41 observed the floor is **architectually insufficient** because
  even FWHM = 0.1 keV is far below real NaI 63x63 FWHM at 60 keV
  (~12 keV) -> match window collapses to 0.3 keV vs required ~36 keV
  -> Am-241 59.54 keV characteristic line never matches).

Fix (this branch ``agent-a-wave3-c2-fwhm``):
  ``build_fwhm_model`` now sanity-checks the fitted quadratic at the
  identification-critical energies (60, 100, 200, 500, 1000 keV). If
  FWHM^2(E) < 1.0 keV^2 at any test E, the fit is rejected and we fall
  back via ``_resolve_pathology_fallback``:
    1. Convert the operator-stored LSRM sqrt(E) polynomial
       (``StoredFwhmCalibration.coefficients``, model
       ``lsrm_fwhm_polynomial_in_E``) to FWHM^2(E) quadratic form by
       sampling-then-refitting at energies E >= 50 keV.
    2. Cross-check the converted model at E=60 keV against
       ``_DEFAULT_NAI_FWHM_MODEL``. If the converted model gives
       FWHM(60) >= 70 % of default-NaI FWHM(60), use it; otherwise the
       polynomial extrapolates below its calibration range and we use
       ``_DEFAULT_NAI_FWHM_MODEL`` (fit on 26 anchors across
       84-2614 keV) instead.
    3. If no usable LSRM stored polynomial is present, fall back to
       ``_DEFAULT_NAI_FWHM_MODEL`` directly.
  Every fallback pushes an entry into
  ``spec.extras["fwhm_model_warnings"]`` for downstream surfacing.

Verification target (brief section "Verification"):
  ``fwhm_keV_at_energy(model, 59.5)`` MUST return >= 10.0 keV for the
  AmTiCsEu Marinelli fixture after fix.
"""
from __future__ import annotations

import math
import os
import sys
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import pytest

SCRIPTS = os.path.join(
    os.path.dirname(__file__), os.pardir, os.pardir, "scripts",
)
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)

from gamma.identification.staged_pipeline import (  # noqa: E402
    _DEFAULT_NAI_FWHM_MODEL,
    _FWHM_PATHOLOGY_TEST_ENERGIES_keV,
    _FWHM_PATHOLOGY_VAL_THRESHOLD_keV2,
    _convert_lsrm_sqrt_E_to_fwhm2_quadratic,
    _eval_fwhm2_quadratic,
    _eval_lsrm_sqrt_E_polynomial,
    _model_is_pathological,
    build_fwhm_model,
    fwhm_keV_at_energy,
    FwhmModel,
)
from gamma.spectrum import StoredFwhmCalibration  # noqa: E402


# --- Real AmTiCsEu fixture data --------------------------------------------
# Verbatim PEAKS= block from
# ``detectors/Gamma-1S/reference_spectra/archive/Поверка-2016/Маринелли/
# Смесь_AmTiCsEu_Маринелли.spe`` (15 rows). Reproducing this table here
# (instead of opening the .spe at test time) lets the test run without
# touching the operator-confidential fixture path on disk.
# Quadratic fit on (E, FWHM^2) -> coefs (-174.497, 2.668, 0.000401) ->
# val(60 keV) = -12.99 keV^2 -> pathology guard fires.
_AMTICSEU_PEAKS_FIXTURE = [
    {"energy_keV": 57.159, "fwhm_keV": 6.345},
    {"energy_keV": 72.0,    "fwhm_keV": 8.106},
    {"energy_keV": 115.784, "fwhm_keV": 12.36},
    {"energy_keV": 151.063, "fwhm_keV": 15.205},
    {"energy_keV": 241.359, "fwhm_keV": 21.492},
    {"energy_keV": 344.418, "fwhm_keV": 27.331},
    {"energy_keV": 512.017, "fwhm_keV": 35.331},
    {"energy_keV": 660.908, "fwhm_keV": 41.542},
    {"energy_keV": 775.406, "fwhm_keV": 45.958},
    {"energy_keV": 965.161, "fwhm_keV": 52.595},
    {"energy_keV": 1146.484, "fwhm_keV": 58.519},
    {"energy_keV": 1400.093, "fwhm_keV": 66.215},
    {"energy_keV": 1689.038, "fwhm_keV": 74.277},
    {"energy_keV": 2148.181, "fwhm_keV": 86.224},
    {"energy_keV": 2641.234, "fwhm_keV": 98.134},
]
# Verbatim FWHM= block from the same .spe file (sqrt(E) polynomial,
# operator-certified, model ``lsrm_fwhm_polynomial_in_E``).
_AMTICSEU_LSRM_COEFS = (-6.2914323508675, 1.7442782325498, 0.0049483216398)


# --- Mock spec helper -------------------------------------------------------

@dataclass
class _MockSpec:
    """Minimal stand-in for ``Spectrum`` carrying only the fields touched
    by ``build_fwhm_model``: ``extras`` (peak table + warnings) and
    ``stored_fwhm_calibration``.
    """
    extras: Dict[str, Any] = field(default_factory=dict)
    stored_fwhm_calibration: Optional[StoredFwhmCalibration] = None


def _make_spec(
    peaks_table: List[Dict[str, float]],
    lsrm_coefs: Optional[tuple] = None,
) -> _MockSpec:
    sf = (
        StoredFwhmCalibration(
            coefficients=lsrm_coefs,
            model="lsrm_fwhm_polynomial_in_E",
        )
        if lsrm_coefs is not None
        else None
    )
    return _MockSpec(
        extras={"lsrm_peaks_table": peaks_table},
        stored_fwhm_calibration=sf,
    )


# --- Pathology detector unit tests -----------------------------------------

class TestModelIsPathological:
    """Unit tests for ``_model_is_pathological``."""

    def test_amticseu_fixture_is_pathological(self):
        """The AmTiCsEu Marinelli fixture quadratic
        ``val(E) = -174.5 + 2.67*E + 0.0004*E^2`` gives val(60) = -12.86
        -> below 1.0 threshold at the lowest test energy.
        """
        model = (-174.5, 2.67, 0.0004)
        assert _model_is_pathological(model) is True
        # Explicit val at the failing test energy
        a, b, c = model
        E = 60.0
        assert a + b * E + c * E * E == pytest.approx(-12.86, abs=0.01)

    def test_default_NaI_model_not_pathological(self):
        """The canonical ``_DEFAULT_NAI_FWHM_MODEL`` must never trip the
        pathology guard.
        """
        assert _model_is_pathological(_DEFAULT_NAI_FWHM_MODEL) is False

    def test_well_conditioned_4anchor_fit_not_pathological(self):
        """A FWHM^2 quadratic fitted from 4 anchors covering 60 .. 1332 keV
        with realistic NaI 63x63 FWHM values must not be flagged.
        """
        import numpy as np

        peaks = [(60.0, 12.0), (122.0, 15.0), (661.66, 42.0), (1332.5, 60.0)]
        Es = np.array([p[0] for p in peaks])
        Fs = np.array([p[1] for p in peaks])
        A = np.vstack([np.ones_like(Es), Es, Es ** 2]).T
        coefs, *_ = np.linalg.lstsq(A, Fs ** 2, rcond=None)
        assert _model_is_pathological(tuple(coefs)) is False

    def test_threshold_constant(self):
        """The class-level constant must remain 1.0 keV^2."""
        assert _FWHM_PATHOLOGY_VAL_THRESHOLD_keV2 == 1.0

    def test_test_energies_exclude_below_60keV(self):
        """The pathology test grid must NOT probe below 60 keV, otherwise
        legitimate quadratic-in-E approximations of sqrt(E)-shaped FWHM
        curves will false-positive (see comment block in module).
        """
        assert min(_FWHM_PATHOLOGY_TEST_ENERGIES_keV) >= 60.0


# --- LSRM polynomial evaluation tests --------------------------------------

class TestLsrmSqrtEPolynomial:
    """Unit tests for the LSRM stored polynomial form, per BUG-22
    convention (``scripts/gamma/io/lsrm_spe.py:44-68``):
    ``FWHM_keV(E) = sum_k c_k * sqrt(E_keV) ** k``.
    """

    def test_lsrm_polynomial_form_uses_sqrt_E(self):
        """For coefs ``(-6.2914, 1.7443, 0.004948)`` (AmTiCsEu fixture)
        the value at E=60 keV must equal
        ``-6.2914 + 1.7443*sqrt(60) + 0.004948*60``.
        """
        coefs = (-6.2914, 1.7443, 0.004948)
        E = 60.0
        expected = -6.2914 + 1.7443 * math.sqrt(E) + 0.004948 * E
        got = _eval_lsrm_sqrt_E_polynomial(coefs, E)
        assert got == pytest.approx(expected, abs=1e-9)

    def test_lsrm_polynomial_at_661(self):
        """LSRM polynomial at Cs-137 661.66 keV anchor: should be near
        the real NaI 63x63 FWHM of ~42 keV.
        """
        coefs = (-6.2914, 1.7443, 0.004948)
        got = _eval_lsrm_sqrt_E_polynomial(coefs, 661.66)
        assert 40.0 < got < 45.0


class TestConvertLsrmSqrtEToFwhm2Quadratic:
    """Unit tests for the sqrt(E)-polynomial -> FWHM^2 quadratic
    conversion via sample-then-refit."""

    def test_conversion_returns_three_tuple(self):
        coefs = (-6.2914, 1.7443, 0.004948)
        out = _convert_lsrm_sqrt_E_to_fwhm2_quadratic(coefs)
        assert out is not None
        assert len(out) == 3
        assert all(isinstance(c, float) for c in out)

    def test_conversion_matches_lsrm_at_661(self):
        """The converted quadratic at E=661.66 keV must reproduce the LSRM
        polynomial value within sampling residuals (< 1 keV).
        """
        coefs = (-6.2914, 1.7443, 0.004948)
        converted = _convert_lsrm_sqrt_E_to_fwhm2_quadratic(coefs)
        f_lsrm = _eval_lsrm_sqrt_E_polynomial(coefs, 661.66)
        f_conv = math.sqrt(
            max(_eval_fwhm2_quadratic(converted, 661.66), 0.01)
        )
        assert abs(f_conv - f_lsrm) < 1.0

    def test_conversion_returns_none_for_empty_coefs(self):
        assert _convert_lsrm_sqrt_E_to_fwhm2_quadratic(()) is None

    def test_conversion_returns_none_for_all_negative_coefs(self):
        """All-negative coefs -> all sample points <= 0 -> too few valid
        samples -> None.
        """
        coefs = (-100.0, -10.0, -0.1)
        assert _convert_lsrm_sqrt_E_to_fwhm2_quadratic(coefs) is None


# --- Integration: build_fwhm_model end-to-end ------------------------------

class TestBuildFwhmModelPathologyFallback:
    """End-to-end behaviour of ``build_fwhm_model`` on pathological
    fixtures."""

    def test_amticseu_pathological_fit_falls_back_to_default_NaI(self):
        """The AmTiCsEu fixture: real 15-anchor PEAKS= table produces the
        ``(-174.497, 2.668, 0.000401)`` quadratic with val(60)=-12.99 keV^2.
        After fix: pathology detected -> LSRM polynomial converted ->
        cross-check fails (LSRM extrapolates too low at E=60) -> default
        NaI 63x63 used -> FWHM(59.5) >= 10 keV.

        Peaks taken from
        ``detectors/Gamma-1S/reference_spectra/archive/Poverka-2016/
        Marinelli/Smesh_AmTiCsEu_Marinelli.spe`` PEAKS= block (15 rows).
        """
        peaks = _AMTICSEU_PEAKS_FIXTURE
        spec = _make_spec(peaks, lsrm_coefs=_AMTICSEU_LSRM_COEFS)

        model, source = build_fwhm_model(spec)
        # Source label records the rejection.
        assert "rejected_to" in source
        # Returned model must give FWHM(59.5) >= 10 keV (brief target).
        fwhm = fwhm_keV_at_energy(model, 59.5)
        assert fwhm >= 10.0, (
            f"FWHM(59.5) = {fwhm:.3f} keV < 10 keV target after fix"
        )
        # A warning was recorded.
        warnings = spec.extras.get("fwhm_model_warnings", [])
        assert len(warnings) >= 1
        assert "BUG-41" in warnings[0]

    def test_amticseu_no_stored_calibration_falls_back_to_default(self):
        """Same pathological PEAKS= table (real AmTiCsEu fixture), but no
        ``StoredFwhmCalibration``: must still fall back to
        ``_DEFAULT_NAI_FWHM_MODEL`` directly.
        """
        peaks = _AMTICSEU_PEAKS_FIXTURE
        spec = _make_spec(peaks, lsrm_coefs=None)
        model, source = build_fwhm_model(spec)
        # F-452: build_fwhm_model теперь возвращает FwhmModel; для
        # default_NaI fallback это quad-форма с теми же legacy coefs.
        assert isinstance(model, FwhmModel)
        assert model.kind == "quad_fwhm2_in_E"
        assert model.coefficients == _DEFAULT_NAI_FWHM_MODEL
        assert "rejected_to_default_NaI_63x63" in source
        fwhm = fwhm_keV_at_energy(model, 59.5)
        assert fwhm >= 10.0

    def test_well_conditioned_fit_is_kept_unchanged(self):
        """A well-conditioned 4-anchor PEAKS= table must NOT be flagged
        as pathological and must be returned with the original source
        label ``lsrm_peaks_table_quadratic``.
        """
        peaks = [
            {"energy_keV": 60.0, "fwhm_keV": 12.0},
            {"energy_keV": 122.0, "fwhm_keV": 15.0},
            {"energy_keV": 661.66, "fwhm_keV": 42.0},
            {"energy_keV": 1332.5, "fwhm_keV": 60.0},
        ]
        spec = _make_spec(peaks, lsrm_coefs=(-6.2914, 1.7443, 0.004948))
        model, source = build_fwhm_model(spec)
        assert source == "lsrm_peaks_table_quadratic"
        assert "rejected_to" not in source
        # FWHM(59.5) ~ 10.08 keV per arithmetic check
        fwhm = fwhm_keV_at_energy(model, 59.5)
        assert fwhm >= 10.0
        # No warnings written for the good-fit path.
        assert "fwhm_model_warnings" not in spec.extras or not spec.extras[
            "fwhm_model_warnings"
        ]

    def test_no_peaks_table_returns_default(self):
        """Empty / missing ``lsrm_peaks_table`` -> default model.

        F-160 (2026-06-20): bootstrap-фоллбек на default_NaI_63x63
        теперь намеренно visible (F-160 ALERT), потому что оператор
        зафиксировал «всегда нужно делать калибровку по FWHM,
        полагаться на расчётную кривую нельзя». Тест допускает
        warning о fallback, но проверяет что модель/источник
        соответствуют default-ветке.
        """
        spec = _MockSpec(extras={}, stored_fwhm_calibration=None)
        model, source = build_fwhm_model(spec)
        # F-452: FwhmModel wrapper around legacy default tuple.
        assert isinstance(model, FwhmModel)
        assert model.kind == "quad_fwhm2_in_E"
        assert model.coefficients == _DEFAULT_NAI_FWHM_MODEL
        assert source == "default_NaI_63x63"
        # F-160: warning о fallback допустим (и желателен); если есть —
        # должен явно сигнализировать о bootstrap rejection.
        warnings = spec.extras.get("fwhm_model_warnings", [])
        for w in warnings:
            assert "F-160" in w or "bootstrap" in w.lower(), (
                f"Unexpected warning on default-fallback path: {w!r}"
            )

    def test_alpha_sqrt_E_path_pathology_guard(self):
        """1-2-anchor alpha*sqrt(E) path also passes through the
        pathology guard. A single anchor at high E with very low FWHM is
        unrealistic and might trip the guard; with realistic 2-anchor
        input the guard must NOT trip.
        """
        # Realistic 2-anchor: 661 keV with FWHM 42 keV (NaI 63x63-like).
        peaks = [
            {"energy_keV": 661.66, "fwhm_keV": 42.0},
        ]
        spec = _make_spec(peaks, lsrm_coefs=None)
        model, source = build_fwhm_model(spec)
        # alpha = FWHM/sqrt(E) = 42/sqrt(661.66) = 1.6329 -> alpha^2 ~ 2.666
        # val(60) = 2.666 * 60 = 159.96 keV^2 -> FWHM ~ 12.6 keV >> floor
        assert source == "lsrm_peaks_table_alpha_sqrt_E"
        assert "rejected_to" not in source
        fwhm = fwhm_keV_at_energy(model, 59.5)
        assert fwhm >= 10.0


# --- Brief target: FWHM(59.5) >= 10 keV -- top-level sanity ----------------

def test_brief_target_fwhm_at_59p5_keV_geq_10():
    """Highest-level brief verification: for the AmTiCsEu fixture after
    fix, ``fwhm_keV_at_energy(model, 59.5)`` MUST return >= 10.0 keV.

    Uses the same pathological PEAKS= table that produced BUG-37/41 in
    the v1.21.0 baseline + the operator-certified LSRM stored
    polynomial. Documents the manual arithmetic in the assert message.

    Manual arithmetic (default NaI 63x63 fallback path, expected
    behaviour because LSRM-converted FWHM(60) ~ 5.10 keV is < 70 % of
    default NaI FWHM(60) ~ 13.38 keV -> cross-check fails -> default
    NaI):
        FWHM^2(59.5) = 0.0 + 2.950048 * 59.5 + 0.000576400 * 59.5^2
                     = 0.0 + 175.5279 + 2.04047
                     = 177.5684
        FWHM(59.5)   = sqrt(177.5684) = 13.3255 keV
    """
    peaks = _AMTICSEU_PEAKS_FIXTURE
    spec = _make_spec(peaks, lsrm_coefs=_AMTICSEU_LSRM_COEFS)
    model, _src = build_fwhm_model(spec)
    fwhm = fwhm_keV_at_energy(model, 59.5)
    assert fwhm >= 10.0, (
        f"BUG-41 brief target: FWHM(59.5) = {fwhm:.4f} keV < 10 keV. "
        f"Expected ~13.33 keV from default NaI 63x63 fallback. Source "
        f"label was '{_src}'."
    )
    # Tighter expectation: should be very close to 13.33 keV
    assert fwhm == pytest.approx(13.325, abs=0.01)
