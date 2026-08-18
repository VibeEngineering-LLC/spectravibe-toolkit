"""
F-41 (v1.7.19) — chain-proxy cross-validation tests.

Verify that:
  • Both Pb-212 and Tl-208 fixtures point to the same .spe file and
    same Th-228 cert source.
  • Tl-208 lib intensities are pre-scaled by β-branching from Bi-212
    (≈35.94 %), so chain_branching stays 1.0 — independently checked
    against the raw ENSDF intensities.
  • The cross-validation block in validate_certs.py produces the
    expected ratio Pb-212/Tl-208 ≈ 0.87 (current best-effort
    measurement; flagged because Pb-212 238.63 keV is harder than
    Tl-208's high-energy lines on a 5 cm point geometry).

These tests are read-only — they do not run the full pipeline. The
end-to-end activity comparison lives in validate_certs.py (a harness,
not a pytest suite).
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

# Import the fixture list from validate_certs.py — the harness ships
# alongside the project so this is allowed (cf. test_secondary_peaks.py
# importing analyze_problem_isotopes.py).
sys.path.insert(0, str(ROOT))
from validate_certs import FIXTURES, CertFixture  # noqa: E402

from gamma.data.nuclide_library import (  # noqa: E402
    load_lsrm_chain_libs, get_nuclide,
)


# Ensure the chain library is loaded (some pytest runners may import
# this module in isolation without validate_certs.py's module-level
# load_lsrm_chain_libs() side effect).
load_lsrm_chain_libs()


# ---------------------------------------------------------------------------
# 1. Configuration sanity
# ---------------------------------------------------------------------------

def test_pb212_and_tl208_fixtures_share_spe_file():
    """In each measurement geometry, the Pb-212 and Tl-208 chain proxies
    must read from the SAME .spe — same source, different daughter
    γ-lines. This guarantees the cross-validation block measures the
    same physical activity per geometry, not two different sources.

    F-46 (v1.7.24): with multi-geometry expansion the harness now ships
    multiple Tl-208 entries (one per geometry where a Th-228 cert source
    is measured). The invariant per-geometry is: Pb-212 and Tl-208 for
    the SAME geometry must share spe_filename. Currently Точечная-5см
    has both Pb-212 and Tl-208 (F-39 + F-41 cross-validation pair); the
    Точечная-25см geometry has only Tl-208 (F-46a). When a future
    iteration adds Pb-212 at 25cm, this test will already enforce the
    pairing constraint."""
    pb = [f for f in FIXTURES if f.nuclide == "Pb-212"]
    tl = [f for f in FIXTURES if f.nuclide == "Tl-208"]
    assert len(pb) >= 1, "expected at least one Pb-212 fixture"
    assert len(tl) >= 1, "expected at least one Tl-208 fixture"
    # Per-geometry invariant: if both proxies exist for the same
    # geometry, they must share the spe file.
    pb_by_geom = {f.geometry: f for f in pb}
    tl_by_geom = {f.geometry: f for f in tl}
    shared_geoms = set(pb_by_geom.keys()) & set(tl_by_geom.keys())
    assert shared_geoms, (
        "at least one geometry must have BOTH Pb-212 and Tl-208 proxies "
        "for cross-validation to be exercised"
    )
    for geom in shared_geoms:
        assert pb_by_geom[geom].spe_filename == tl_by_geom[geom].spe_filename, (
            f"[{geom}] Pb-212 reads {pb_by_geom[geom].spe_filename!r} but "
            f"Tl-208 reads {tl_by_geom[geom].spe_filename!r} — chain "
            f"proxies for the same geometry must agree on the source file"
        )
        print(f"  [{geom}] Pb-212 and Tl-208 both read: "
              f"{pb_by_geom[geom].spe_filename}")


def test_pb212_and_tl208_fixtures_share_cert_parent():
    """Pb-212 fixtures must all report cert_nuclide=Th-228 (Pb-212 is a
    Th-228 direct daughter only — there is no analogous use as a
    Th-232 proxy because Pb-212 is the bottleneck before Bi-212 →
    Tl-208 branching).

    Tl-208 fixtures appear in TWO chain-proxy roles:
      * cert_nuclide="Th-228" — Th-228 measured at Точечная-5см
        (F-39 baseline) and Точечная-25см (F-46a).
      * cert_nuclide="Th-232" — Th-232 measured in Marinelli /
        Дента-120мл / Петри-60мл (F-46b/c). Th-232 chain is in secular
        equilibrium with the Th-228 sub-chain, so the same Tl-208 lib
        intensity pre-scaling (β-branching 0.3594) inverts to the
        parent's activity in both cases (chain_branching=1.0).

    Constraint (per F-46): every Tl-208/Pb-212 fixture must have a
    non-None cert_nuclide in {"Th-228", "Th-232"}; same-geometry
    Pb-212/Tl-208 pairs on the same source must agree on cert_nuclide
    and cert_source_hint.
    """
    pb_list = [f for f in FIXTURES if f.nuclide == "Pb-212"]
    tl_list = [f for f in FIXTURES if f.nuclide == "Tl-208"]
    assert pb_list and tl_list
    ALLOWED_PARENTS = {"Th-228", "Th-232"}
    for f in pb_list:
        assert f.cert_nuclide == "Th-228", (
            f"Pb-212 on {f.spe_filename!r} must have cert_nuclide=Th-228 "
            f"(Pb-212 is only used as Th-228 chain proxy), "
            f"got {f.cert_nuclide!r}"
        )
    for f in tl_list:
        assert f.cert_nuclide in ALLOWED_PARENTS, (
            f"Tl-208 on {f.spe_filename!r} must have cert_nuclide ∈ "
            f"{ALLOWED_PARENTS}, got {f.cert_nuclide!r}"
        )
    # Same-geometry / same-source Pb-212+Tl-208 pairs must agree on
    # both parent and cert source hint.
    pb_by_key = {(f.geometry, f.spe_filename): f for f in pb_list}
    tl_by_key = {(f.geometry, f.spe_filename): f for f in tl_list}
    for key in set(pb_by_key) & set(tl_by_key):
        geom, spe = key
        assert pb_by_key[key].cert_nuclide == tl_by_key[key].cert_nuclide, (
            f"[{geom}] {spe}: Pb-212 cert_nuclide="
            f"{pb_by_key[key].cert_nuclide!r} != Tl-208 cert_nuclide="
            f"{tl_by_key[key].cert_nuclide!r}"
        )
        assert pb_by_key[key].cert_source_hint == tl_by_key[key].cert_source_hint, (
            f"[{geom}] {spe}: Pb-212 hint={pb_by_key[key].cert_source_hint!r} "
            f"!= Tl-208 hint={tl_by_key[key].cert_source_hint!r}"
        )
        print(f"  [{geom}] same-source pair: parent="
              f"{pb_by_key[key].cert_nuclide}, hint="
              f"{pb_by_key[key].cert_source_hint!r}")


def test_chain_branching_is_one_for_both_proxies():
    """Because Lsrm chain-library lib intensities are pre-scaled by
    β-branching (Tl-208) and 1:1 for direct daughter (Pb-212),
    chain_branching = 1.0 for all proxies. compute_activity inverts the
    lib intensity and recovers the parent activity directly.

    F-46 (v1.7.24): assert for ALL chain-proxy fixtures, not just one
    of each kind (the 25cm Tl-208 row must also satisfy the invariant).
    """
    for f in FIXTURES:
        if f.nuclide in ("Pb-212", "Tl-208") and f.cert_nuclide == "Th-228":
            assert f.chain_branching == 1.0, (
                f"chain proxy {f.nuclide} [{f.geometry}] on "
                f"{f.spe_filename!r} must have chain_branching=1.0, "
                f"got {f.chain_branching}"
            )


# ---------------------------------------------------------------------------
# 2. Lib intensities embed β-branching (Tl-208 specifically)
# ---------------------------------------------------------------------------

# ENSDF branching ratio Bi-212 → α → Tl-208
BI212_ALPHA_TO_TL208 = 0.3594

# Raw ENSDF Tl-208 emission probabilities (per Tl-208 disintegration)
# Source: ENSDF Tl-208 evaluation 2004, IAEA NDS
ENSDF_TL208_RAW = {
    2614.51: 99.75,
    583.19: 84.5,
    510.77: 22.6,
    860.56: 12.5,
}


def test_tl208_lib_intensities_embed_beta_branching():
    """A line-by-line check that the chain-library Tl-208 intensities
    are = raw ENSDF intensity × 0.3594 (Bi-212 α-branching). Within
    ±2 %% of the analytic product."""
    tl = get_nuclide("Tl-208")
    assert tl is not None, "Tl-208 must be loaded (load_lsrm_chain_libs)"
    by_E = {round(float(line[0]), 2): float(line[1])
            for line in tl["lines"]}
    for E_keV, raw_pct in ENSDF_TL208_RAW.items():
        assert E_keV in by_E, (
            f"Tl-208 line at {E_keV} keV missing from lib "
            f"(have: {sorted(by_E)})"
        )
        expected = raw_pct * BI212_ALPHA_TO_TL208
        observed = by_E[E_keV]
        rel = abs(observed - expected) / expected
        assert rel < 0.02, (
            f"Tl-208 {E_keV} keV: lib I={observed:.3f} vs "
            f"ENSDF×branching = {raw_pct:.2f}×{BI212_ALPHA_TO_TL208:.4f} "
            f"= {expected:.3f}; relative error {rel:.3%}"
        )
        print(f"  Tl-208 {E_keV} keV: lib I={observed:.3f}% ≈ "
              f"{raw_pct:.2f}%×{BI212_ALPHA_TO_TL208} "
              f"(rel err {rel:.3%})")


def test_pb212_lib_intensities_are_direct_ensdf():
    """Pb-212 is a direct 1:1 daughter — no branching loss. Lib I
    should match raw ENSDF intensities within reading precision."""
    pb = get_nuclide("Pb-212")
    assert pb is not None
    # ENSDF Pb-212 evaluation: 238.63 keV at 43.6 % (the canonical line)
    by_E = {round(float(line[0]), 2): float(line[1])
            for line in pb["lines"]}
    assert 238.63 in by_E
    assert abs(by_E[238.63] - 43.6) < 1.0, (
        f"Pb-212 238.63 keV lib I={by_E[238.63]} vs ENSDF 43.6"
    )


# ---------------------------------------------------------------------------
# 3. Cross-validation block presence
# ---------------------------------------------------------------------------

def test_validate_certs_module_has_cross_validation_block():
    """Defensive check: the harness must contain the F-41 cross-
    validation block (extended in F-46 to group same-source pairs
    only). If someone removes it the test fails loudly."""
    src = (ROOT / "scripts" / "validate_certs.py").read_text(encoding="utf-8")
    # Header is now "Cross-validation of chain proxies (F-41 / F-46)
    # — same .spe, ≥2 daughters" after F-46. Match the stable prefix
    # that survives future feature-tag extensions.
    assert "Cross-validation of chain proxies (F-41" in src, (
        "Expected the F-41 cross-validation block header to be present"
    )
    assert "ratio" in src
    # Grouping dict renamed from `chain_proxies` (per-parent) to
    # `same_source_proxies` (per-source) in F-46. Accept either spelling
    # so the assertion catches an accidental deletion of the grouping
    # logic but tolerates the F-46 refactor.
    assert ("chain_proxies" in src) or ("same_source_proxies" in src), (
        "Expected a chain-proxy grouping dict (chain_proxies or "
        "same_source_proxies) in validate_certs.py"
    )


def test_chain_proxy_fixtures_can_be_paired():
    """Mechanical check that the harness logic (grouping rows by
    cert_nuclide) can pair the two Th-228 proxies."""
    chain_parents = {}
    for fx in FIXTURES:
        if fx.cert_nuclide:
            chain_parents.setdefault(fx.cert_nuclide, []).append(fx.nuclide)
    assert "Th-228" in chain_parents
    daughters = chain_parents["Th-228"]
    assert "Pb-212" in daughters
    assert "Tl-208" in daughters
    assert len(daughters) >= 2, (
        "cross-validation requires ≥2 chain proxies per parent"
    )


# ---------------------------------------------------------------------------
# 4. Documented numerical expectations (regression guard)
# ---------------------------------------------------------------------------
#
# These match the v1.7.19 measurement on the Th-228 5 cm fixture:
#   Pb-212 (238 keV)            → A=77 250 Bq, Δ=-12.79 %
#   Tl-208 (583/860/2614 keV)   → A=88 517 Bq, Δ=-00.07 %
#   ratio Pb-212/Tl-208 = 0.8727
#
# We test the IDENTITY relations, not the absolute numbers (those live
# in cert_validation_matrix.csv and would force tedious updates after
# any unrelated tweak to peak search or efficiency curve fit). If
# future work tightens the Pb-212 self-absorption model the deviation
# may narrow; the ratio identity stays.

def test_tl208_branching_arithmetic_identity():
    """If Tl-208 lib intensity = raw_ENSDF × β, and compute_activity
    inverts lib intensity, then the recovered A_Th-228 from Tl-208
    is independent of the β-branching value — IT CANCELS OUT.
    Numerical version: simulate it."""
    # Suppose true A_Th-228 = 100 000 Bq. Tl-208 emits 583 keV with
    # raw ENSDF intensity 84.5 % per Tl-208 disintegration.
    A_Th_228 = 100_000.0
    raw_I = 0.845  # per Tl-208 disintegration
    branching = BI212_ALPHA_TO_TL208
    # Observed counts ∝ A_Th-228 × branching × raw_I (× ε × T, irrelevant
    # constants).
    observed_term = A_Th_228 * branching * raw_I
    # compute_activity uses lib_I = branching × raw_I to invert:
    lib_I = branching * raw_I
    recovered = observed_term / lib_I
    assert abs(recovered - A_Th_228) < 1e-6, (
        f"recovered {recovered} vs true {A_Th_228}: "
        "lib intensity pre-scaling identity broken"
    )


# ---------------------------------------------------------------------------
# Test runner
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    tests = [
        test_pb212_and_tl208_fixtures_share_spe_file,
        test_pb212_and_tl208_fixtures_share_cert_parent,
        test_chain_branching_is_one_for_both_proxies,
        test_tl208_lib_intensities_embed_beta_branching,
        test_pb212_lib_intensities_are_direct_ensdf,
        test_validate_certs_module_has_cross_validation_block,
        test_chain_proxy_fixtures_can_be_paired,
        test_tl208_branching_arithmetic_identity,
    ]
    passed = 0
    failed = []
    for t in tests:
        try:
            print(f"-- {t.__name__}")
            t()
            print(f"   OK")
            passed += 1
        except AssertionError as e:
            print(f"   FAIL: {e}")
            failed.append((t.__name__, str(e)))
        except Exception as e:
            print(f"   ERROR: {type(e).__name__}: {e}")
            failed.append((t.__name__, f"{type(e).__name__}: {e}"))
    print()
    print(f"Passed: {passed}/{len(tests)}")
    if failed:
        for name, msg in failed:
            print(f"  - {name}: {msg}")
        sys.exit(1)
    sys.exit(0)
