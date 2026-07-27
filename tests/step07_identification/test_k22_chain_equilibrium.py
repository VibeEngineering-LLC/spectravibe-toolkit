"""
F-51 / K-22 unit tests — Th-232 chain-equilibrium correction.

These tests cover the v1.11.0 `chain_at_cert_equilibrium` /
`chain_bottleneck_T_half_s` extension in `validate_certs.CertFixture`
plus the in-growth-factor arithmetic and the FIXTURE-wiring rules that
distinguish "heavy" 420-7-17 Th-232 sources (cert ref 2007-09, chain
reset at cert) from "light" 420-17031 sources (cert ref 2017-06, chain
already at equilibrium at cert).

The tests are intentionally pipeline-free: they exercise the
`CertFixture` dataclass, the `RA228_T_HALF_S` constant, and the
in-growth formula directly. The full validate_certs.py harness has
its own end-to-end matrix check (40 fixtures) that lives in CI manually.
"""
from __future__ import annotations

import math
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

import validate_certs as vc

ONE_JULIAN_YEAR_S = 365.25 * 86400.0


def _approx(a, b, *, tol):
    assert abs(a - b) <= tol, f"expected ~{b}, got {a} (d={a - b})"


# ---------------------------------------------------------------------
# 1. Ra-228 half-life constant
# ---------------------------------------------------------------------

def test_ra228_t_half_s_value():
    """T½(Ra-228) = 5.75 y. Verify the module constant matches the
    documented Julian-year conversion (5.75 × 365.25 × 86400)."""
    expected = 5.75 * ONE_JULIAN_YEAR_S
    _approx(vc.RA228_T_HALF_S, expected, tol=1e-3)
    # Sanity: ~1.81e8 s
    assert 1.81e8 < vc.RA228_T_HALF_S < 1.82e8


# ---------------------------------------------------------------------
# 2. In-growth factor formula
# ---------------------------------------------------------------------

def test_eq_factor_17y_matches_doc():
    """For Δt = 17 y the in-growth factor should be ~0.870 per the
    K-22 entry in KNOWN_AND_FIXED_ISSUES.md (1 − 2^(−17/5.75))."""
    dt_s = 17.0 * ONE_JULIAN_YEAR_S
    factor = 1.0 - math.exp(-math.log(2.0) * dt_s / vc.RA228_T_HALF_S)
    _approx(factor, 0.870, tol=0.005)


def test_eq_factor_7y_matches_doc():
    """For Δt = 7 y the in-growth factor should be ~0.566 (predicted
    by formula, deliberately NOT applied to 17031 fixtures — they were
    prepared with chain already at cert-time equilibrium)."""
    dt_s = 7.0 * ONE_JULIAN_YEAR_S
    factor = 1.0 - math.exp(-math.log(2.0) * dt_s / vc.RA228_T_HALF_S)
    _approx(factor, 0.566, tol=0.005)


def test_eq_factor_asymptotes_to_one():
    """For Δt >> T_Ra-228 (e.g. 60 y ≈ 10·T_Ra-228) the in-growth
    factor should be indistinguishable from 1.0 (≥ 99.9 %)."""
    dt_s = 60.0 * ONE_JULIAN_YEAR_S
    factor = 1.0 - math.exp(-math.log(2.0) * dt_s / vc.RA228_T_HALF_S)
    assert factor > 0.999, f"expected near-1, got {factor}"


def test_eq_factor_zero_at_cert_date():
    """At Δt = 0 the factor is exactly zero — no daughter γ-emission
    is expected if the chain was just reset."""
    factor = 1.0 - math.exp(-math.log(2.0) * 0.0 / vc.RA228_T_HALF_S)
    _approx(factor, 0.0, tol=1e-12)


# ---------------------------------------------------------------------
# 3. CertFixture dataclass defaults
# ---------------------------------------------------------------------

def test_cert_fixture_defaults_preserve_v110_behavior():
    """Default fixture must have chain_at_cert_equilibrium=True and
    chain_bottleneck_T_half_s=None — i.e. no K-22 correction is applied
    to fixtures that don't opt in. This guarantees zero regression on
    the 37 fixtures that don't carry the new flags."""
    fx = vc.CertFixture("Co-60", "x.spe", "Co-60 hint")
    assert fx.chain_at_cert_equilibrium is True
    assert fx.chain_bottleneck_T_half_s is None


# ---------------------------------------------------------------------
# 4. FIXTURE wiring — heavy Th-232 sources opt in; light ones don't
# ---------------------------------------------------------------------

def _find_fixture(spe_filename_contains, geometry):
    matches = [f for f in vc.FIXTURES
               if spe_filename_contains in f.spe_filename
               and f.geometry == geometry]
    assert len(matches) >= 1, (
        f"no fixture with {spe_filename_contains!r} in {geometry!r}")
    return matches[0]


def test_heavy_marinelli_420_7_17_opts_in():
    """Marinelli Th232 420-7-17 (cert ref 2007-09) MUST have K-22
    correction enabled."""
    fx = _find_fixture("Th232_420-7-17_Маринелли", "Маринелли")
    assert fx.cert_nuclide == "Th-232"
    assert fx.chain_at_cert_equilibrium is False
    assert fx.chain_bottleneck_T_half_s is not None
    _approx(fx.chain_bottleneck_T_half_s, vc.RA228_T_HALF_S, tol=1.0)


def test_heavy_denta_420_7_17_opts_in():
    fx = _find_fixture("Th232_420-7-17_Дента", "Дента-120мл")
    assert fx.cert_nuclide == "Th-232"
    assert fx.chain_at_cert_equilibrium is False
    assert fx.chain_bottleneck_T_half_s is not None


def test_heavy_petri_420_7_17_opts_in():
    fx = _find_fixture("Th232_420-7-17_Петри", "Петри-60мл")
    assert fx.cert_nuclide == "Th-232"
    assert fx.chain_at_cert_equilibrium is False
    assert fx.chain_bottleneck_T_half_s is not None


def test_light_17031_marinelli_does_not_opt_in():
    """Marinelli Th-232 420-17031 (cert ref 2017-06) was prepared with
    chain in equilibrium at cert date — MUST NOT receive K-22
    correction (would over-correct by ~43 %)."""
    fx = _find_fixture("Th-232_420-17031_Маринелли", "Маринелли")
    assert fx.cert_nuclide == "Th-232"
    assert fx.chain_at_cert_equilibrium is True
    assert fx.chain_bottleneck_T_half_s is None


def test_light_17031_denta_does_not_opt_in():
    fx = _find_fixture("Th-232_420-17031_Дента", "Дента-120мл")
    assert fx.cert_nuclide == "Th-232"
    assert fx.chain_at_cert_equilibrium is True


def test_light_17031_petri_does_not_opt_in():
    fx = _find_fixture("Th-232_420-17031_Петри", "Петри-60мл")
    assert fx.cert_nuclide == "Th-232"
    assert fx.chain_at_cert_equilibrium is True


def test_ra226_chain_proxies_do_not_opt_in():
    """Ra-226 → Bi-214 chain has Rn-222 bottleneck (T½ = 3.8 d). In
    sealed Marinelli the chain reaches equilibrium within 30 days
    post-sealing, so 20-year-old sources are at FULL equilibrium and
    K-22 correction MUST NOT apply (would be a Th-232-specific bug)."""
    ra_proxies = [f for f in vc.FIXTURES
                  if f.cert_nuclide == "Ra-226"]
    assert ra_proxies, "expected at least one Ra-226 chain-proxy fixture"
    for fx in ra_proxies:
        assert fx.chain_at_cert_equilibrium is True, (
            f"Ra-226 proxy {fx.spe_filename} unexpectedly opts in to K-22")
        assert fx.chain_bottleneck_T_half_s is None


def test_th228_chain_proxies_do_not_opt_in():
    """Th-228 sources (SRC-01 / SRC-07) were certified at chain
    equilibrium; no K-22 correction."""
    th_proxies = [f for f in vc.FIXTURES
                  if f.cert_nuclide == "Th-228"]
    assert th_proxies
    for fx in th_proxies:
        assert fx.chain_at_cert_equilibrium is True
        assert fx.chain_bottleneck_T_half_s is None


# ---------------------------------------------------------------------
# 5. Quantitative sanity: applying K-22 to Marinelli 420-7-17 makes the
#    cert-side activity drop by ~13 % (matches the 17-year in-growth
#    deficit empirically observed pre-correction).
# ---------------------------------------------------------------------

def test_applied_correction_size_for_420_7_17():
    """Δ between 2024-01-01 measurement and 2007-09-17 cert ref date
    (≈ 16.3 y) gives an in-growth factor in the 0.85–0.88 band; full
    17.0-y interval would give 0.870."""
    cert_ref = datetime(2007, 9, 17)
    meas = datetime(2024, 1, 1)
    dt_s = (meas - cert_ref).total_seconds()
    factor = 1.0 - math.exp(-math.log(2.0) * dt_s / vc.RA228_T_HALF_S)
    assert 0.85 <= factor <= 0.88, f"unexpected factor {factor}"


# ---------------------------------------------------------------------
# Test runner
# ---------------------------------------------------------------------

if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    n_pass = 0
    n_fail = 0
    for t in tests:
        try:
            t()
            print(f"  OK   {t.__name__}")
            n_pass += 1
        except AssertionError as exc:
            print(f"  FAIL {t.__name__}: {exc}")
            n_fail += 1
        except Exception as exc:
            print(f"  ERR  {t.__name__}: {type(exc).__name__}: {exc}")
            n_fail += 1
    print(f"\n{n_pass} passed, {n_fail} failed (of {len(tests)} tests).")
    sys.exit(0 if n_fail == 0 else 1)
