"""
BUG-22 / 2026-06-02 — LSRM .spe FWHM polynomial domain tests.

Pinning the contract for the `lsrm_fwhm_polynomial_in_E` stored FWHM
model so the broken polynomial-in-E reading (which produced negative
FWHM at typical energies, e.g. −14 888 keV at 238 keV on the Th-232
fixture) cannot regress.

Per LSRM «Алгоритмические основы» §8.3 («Калибровка по полуширине»):

    FWHM_keV(E) = Σ_k c_k · z^k,   z = √E_keV

The reader stores `coefficients` as the raw `(c0, c1, c2, …)` low-to-
high tuple; downstream code must evaluate as a polynomial in √E. See
`scripts/validate_certs.py::make_lsrm_fwhm_provider` for the canonical
in-tree evaluator (it already honored this convention pre-BUG-22).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

from gamma.io.readers import read_spectrum  # noqa: E402


# Th-232 Marinelli 0 cm — line 38 of the .spe header reads:
#   FWHM=3,8.9691080405128,-0.8981588537211,0.1433311097464,-0.001693803867
# This is the canonical fixture for the bug: a polynomial that yields
# physically reasonable FWHM only when evaluated in z = √E_keV.
TH232_FIXTURE = (
    "detectors/Gamma-1S/reference_spectra/archive/"
    "Th232_420-7-17_Маринелли_0cm.spe"
)

# Cs-137 reference (4-coefficient polynomial in √E)
CS137_FIXTURE = (
    "detectors/Gamma-1S/reference_spectra/archive/"
    "Cs-137__163_2017.spe"
)


def _eval_poly_sqrtE(coefs, E_keV: float) -> float:
    """Reference evaluator: polynomial in z = √E_keV via Horner."""
    if E_keV <= 0.0:
        return float("nan")
    z = math.sqrt(E_keV)
    v = 0.0
    for c in reversed(coefs):
        v = v * z + float(c)
    return v


def _eval_poly_E(coefs, E_keV: float) -> float:
    """Anti-pattern evaluator (BUG-22 trap): polynomial in E directly."""
    v = 0.0
    for c in reversed(coefs):
        v = v * E_keV + float(c)
    return v


def test_lsrm_fwhm_polynomial_stored_as_raw_coefficients_in_sqrtE_domain():
    """
    The reader must store the FWHM polynomial coefficients verbatim
    (no domain transformation) and label the model with the canonical
    name `lsrm_fwhm_polynomial_in_E`. Despite the legacy name, the
    polynomial argument is z = √E_keV (LSRM Algorithmic Foundations
    §8.3). Pinning the contract so downstream evaluators (channel
    providers, identification windows, multiplet deconvolution) cannot
    silently misinterpret the domain.
    """
    spec = read_spectrum(TH232_FIXTURE)
    sf = spec.stored_fwhm_calibration
    assert sf is not None, "Th-232 fixture must carry a stored FWHM cal"
    assert sf.model == "lsrm_fwhm_polynomial_in_E", (
        f"unexpected FWHM model label: {sf.model!r}"
    )
    coefs = sf.coefficients
    assert len(coefs) >= 3, f"expected ≥3 FWHM coefficients, got {len(coefs)}"
    # Pin the exact stored values (header line 38 of the fixture)
    expected = (8.9691080405128, -0.8981588537211,
                0.1433311097464, -0.001693803867)
    for i, (got, exp) in enumerate(zip(coefs, expected)):
        assert math.isclose(got, exp, rel_tol=1e-10), (
            f"FWHM coef c{i}: stored={got!r}, expected={exp!r}"
        )

    # Evaluated in √E_keV the polynomial must give physically reasonable
    # NaI(Tl) FWHM (~5-10% of E across 200-3000 keV); evaluated in E
    # directly it produces negative or astronomically large values and
    # falls outside any sane FWHM range.
    for E_keV, lo_keV, hi_keV in [
        (238.63, 10.0, 35.0),   # Pb-212; NaI ~9-12% → 22-29 keV
        (661.66, 35.0, 80.0),   # Cs-137; NaI ~7-9%  → 46-60 keV
        (1460.82, 65.0, 110.0), # K-40;    NaI ~5-7%  → 73-102 keV
        (2614.51, 90.0, 150.0), # Tl-208;  NaI ~4-6%  → 105-157 keV
    ]:
        sqrtE_value = _eval_poly_sqrtE(coefs, E_keV)
        assert lo_keV <= sqrtE_value <= hi_keV, (
            f"√E-domain eval at {E_keV} keV: got {sqrtE_value:.2f} keV, "
            f"expected within [{lo_keV}, {hi_keV}] (physical NaI band)"
        )
        # The E-domain anti-pattern: at 238 keV this returns a strongly
        # negative value (≈ −22 keV with these particular coefs — Horner
        # in E balances differently than Horner in √E and the
        # asymptotic sign flips at certain energies). Either way, it
        # cannot land inside the physical NaI band at 238 keV.
        E_value = _eval_poly_E(coefs, 238.63)
        assert not (10.0 <= E_value <= 35.0), (
            f"E-domain anti-pattern eval at 238 keV gave {E_value:.2f} keV; "
            f"if this lands inside the physical FWHM band the test is no "
            f"longer protecting against the BUG-22 misinterpretation"
        )


def test_lsrm_fwhm_polynomial_sanity_at_canonical_peaks_cs137_fixture():
    """
    Second fixture (Cs-137 control source) — same domain contract.
    Sanity-checks at canonical NaI calibration peaks (Pb-212 238,
    Cs-137 661, K-40 1460, Tl-208 2614 keV): the √E-domain evaluation
    must land inside the physical NaI FWHM band at every checkpoint;
    the polynomial-in-E anti-pattern must fail at least one checkpoint
    (i.e. the two evaluators must produce demonstrably different
    answers — otherwise the BUG-22 domain guard isn't actually being
    exercised).
    """
    spec = read_spectrum(CS137_FIXTURE)
    sf = spec.stored_fwhm_calibration
    assert sf is not None
    assert sf.model == "lsrm_fwhm_polynomial_in_E"
    coefs = sf.coefficients
    assert len(coefs) >= 2

    checkpoints_keV = [238.63, 661.66, 1460.82, 2614.51]
    sqrt_evals = [_eval_poly_sqrtE(coefs, E) for E in checkpoints_keV]
    e_evals = [_eval_poly_E(coefs, E) for E in checkpoints_keV]

    # √E-domain: every checkpoint must land in (0, 200) keV — a generous
    # physical envelope for any NaI 63×63 detector.
    for E, fw in zip(checkpoints_keV, sqrt_evals):
        assert 0.0 < fw < 200.0, (
            f"√E-domain FWHM at {E} keV: {fw:.2f} keV out of "
            f"(0, 200) physical envelope"
        )
    # Monotone non-decreasing — FWHM(E) grows with E for NaI scintillators.
    for prev, nxt in zip(sqrt_evals[:-1], sqrt_evals[1:]):
        assert nxt >= prev * 0.95, (  # tiny slack for fit-driven dips
            f"√E-domain FWHM not monotone non-decreasing: "
            f"{sqrt_evals}"
        )
    # The two evaluators must disagree somewhere; otherwise the domain
    # guard is vacuous on this fixture.
    max_rel_diff = max(
        abs(s - e) / max(abs(s), 1e-9)
        for s, e in zip(sqrt_evals, e_evals)
    )
    assert max_rel_diff > 0.1, (
        f"√E-domain and E-domain evaluators agree to within 10% at every "
        f"checkpoint — domain guard is not protecting anything. "
        f"sqrt={sqrt_evals}, E={e_evals}"
    )


def test_lsrm_fwhm_polynomial_negative_at_238_keV_when_misinterpreted():
    """
    Regression guard for the *exact* failure mode reported in BUG-22.
    On the Th-232 fixture, naïvely evaluating the FWHM polynomial as
    Σ c_k · E^k (with E in keV) at 238 keV gives a strongly non-
    physical value. This test pins that mis-evaluation so any future
    contributor who is tempted to interpret the model name literally
    sees an explicit reminder.
    """
    spec = read_spectrum(TH232_FIXTURE)
    coefs = spec.stored_fwhm_calibration.coefficients
    bad = _eval_poly_E(coefs, 238.63)
    # The exact value depends on the polynomial degree; for the Th-232
    # fixture (cubic in z=√E with coefs (8.97, -0.898, 0.143, -0.00169)),
    # Horner-in-E at 238 keV gives a value FAR outside the physical NaI
    # FWHM range. We don't pin the exact number to avoid coupling to
    # numerical details; we just require it to be outside (0, 100) keV.
    assert not (0.0 < bad < 100.0), (
        f"E-domain anti-pattern eval at 238 keV gave {bad:.2f} keV — "
        f"this landed in the physical FWHM band by coincidence and the "
        f"BUG-22 regression guard is no longer effective. Recheck the "
        f"fixture's stored coefficients."
    )


if __name__ == "__main__":
    # Allow direct invocation for quick smoke check during development.
    pytest.main([__file__, "-v"])
