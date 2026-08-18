"""BUG-27 — Ac-228 1588.2 keV per-line A_i outlier (in-blob multiplet
effective-I correction).

Context. v1.18.31 closed BUG-15 (within-nuclide shared-peak dedup) by
picking the highest-I_pct line as the survivor of an unresolved
multiplet. That fix prevented A_i inflation in the *skipped* low-I
lines but left a residual inflation in the survivor itself: the
survivor receives the full S_net of the unresolved blob while the
denominator A_i = S/(ε·I·t) keeps only the survivor's own
library_I_pct. The physical photon population in S includes:

  • All dedup'd same-nuclide siblings on the same peak_channel
    (their I_pct is missing from the denominator).
  • Chain-equilibrium siblings (e.g. Bi-212 1620.50 keV in the
    Ac-228 1588.20 keV blob at NaI 63×63 FWHM=115 keV) whose
    library lines fall within ±FWHM/2 of the observed centroid.

Closure case (v1.18.31 Th-232 demo): A_i(Ac-228 1588.20) = 3881 Bq/kg
vs nuclide-mean = 1802 Bq/kg (ratio 2.15×). After BUG-27 fix:
A_i(1588.20) ≈ 1827 Bq/kg (within +1.5% of nuclide-mean).

Contracts:
  1. Survivor's A_i denominator uses ΣI_pct of the within-nuclide
     dedup group (winner + all skipped same-nuclide partners sharing
     peak_channel).
  2. For wide deconvolved-coupled multiplets (FWHM ≥ 50 keV,
     peak_area_source starts with "deconvolved"), chain-equilibrium
     sibling library lines within ±FWHM/2 of the observed centroid
     are summed into the effective I as well.
  3. Narrow single-line peaks (no dedup partners, no chain siblings
     in window) are unaffected — backward-compatible with all
     existing single-line activity calculations.
  4. Cs-137 661.66 keV (single isolated line, no chain) → no change.
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path

import pytest

sys.path.insert(
    0, str(Path(__file__).resolve().parent.parent.parent / "scripts")
)

from gamma.activity.compute import (
    compute_activity,
    _chain_sibling_I_in_window,
    _NUCLIDE_TO_CHAIN,
)
from gamma.calibration.efficiency import EfficiencyCurve


def _flat_efficiency_curve(eps: float = 1.0e-3) -> EfficiencyCurve:
    return EfficiencyCurve(
        coefficients=(math.log(eps),),
        E_min_keV=1.0,
        E_max_keV=3000.0,
        chi2_per_dof=1.0,
        n_points_used=10,
        n_dof=9,
        detector_id="test",
        geometry="test",
    )


@dataclass
class _FakeMatch:
    nuclide: str
    library_E_keV: float
    library_I_pct: float
    peak_channel: int
    peak_E_keV: float
    peak_sigma: float = 1.0
    residual_keV: float = 0.0
    is_characteristic: bool = False
    peak_area: float = 0.0
    peak_area_uncertainty: float = 0.0
    peak_area_source: str = "deconvolved_coupled"


@dataclass
class _FakeNuclideId:
    nuclide: str
    detected: bool = True
    reason: str = "test"
    characteristic_line_keV: float = 0.0
    matched_lines: tuple = ()
    confidence: object = None


# ─── Helper-unit tests ───────────────────────────────────────────────────


def test_chain_sibling_I_in_window_finds_bi212_for_ac228_1588_blob():
    """At peak_E_obs=1597.35 keV, FWHM=115.5 keV (the v1.18.31 Th-232
    demo Ac-228 1588-blob), Bi-212 1620.50 (I=1.51%) is the only
    chain-equilibrium sibling inside ±FWHM/2."""
    sigma_I = _chain_sibling_I_in_window(
        nuclide="Ac-228",
        peak_E_keV=1597.35,
        fwhm_keV=115.5,
        exclude_E_keV=[1580.53, 1588.20, 1630.63],
    )
    # Bi-212 1620.50 I=1.51% is the only sibling in this window
    # (1597.35 ± 57.75 keV = 1539.6 – 1655.1 keV)
    assert math.isclose(sigma_I, 1.51, abs_tol=0.01), (
        f"expected Bi-212 1620.50 alone (I=1.51%), got {sigma_I}"
    )


def test_chain_sibling_lookup_skips_non_chain_nuclides():
    """Cs-137 is not in any natural decay chain — no sibling lookup."""
    assert "Cs-137" not in _NUCLIDE_TO_CHAIN
    sigma_I = _chain_sibling_I_in_window(
        nuclide="Cs-137",
        peak_E_keV=661.66,
        fwhm_keV=70.0,
    )
    assert sigma_I == 0.0


def test_chain_sibling_lookup_excludes_same_nuclide():
    """Same-nuclide lines are handled by within-nuclide ΣI, not by
    the chain-sibling sum (avoid double-counting)."""
    # At peak_E=1597, FWHM=115, Ac-228 1580.53/1588.20/1630.63 are
    # same-nuclide. The sibling lookup must NOT include them.
    sigma_I = _chain_sibling_I_in_window(
        nuclide="Ac-228",
        peak_E_keV=1597.35,
        fwhm_keV=115.5,
        exclude_E_keV=[],  # don't even need exclude — same-nuclide
                           # is filtered out unconditionally
    )
    # Only Bi-212 1620.50 (1.51%) should be present, NOT
    # Ac-228 1580.53/1588.20/1630.63 totalling 5.33%.
    assert sigma_I < 5.0, (
        f"sibling sum {sigma_I:.3f}% leaked same-nuclide lines"
    )
    assert sigma_I >= 1.5, (
        f"sibling sum {sigma_I:.3f}% missed Bi-212 1620.50"
    )


# ─── End-to-end activity test (the actual BUG-27 closure case) ───────────


def test_ac228_1588_outlier_converges_to_nuclide_mean_within_15pct():
    """BUG-27 closure: per-line A_i for Ac-228 1588.20 keV in the
    Th-232 demo (peak_channel=539, peak_E=1597.35, FWHM=115.5,
    S=20660.92) must converge to within ±15% of the nuclide-mean
    activity computed from the (un-multipleted) 911.20 keV line.

    Setup mirrors v1.18.31 Th-232 demo numbers:
      live_time_s = 11359.164 s
      ε(911.20)   = 1.346e-2  (.efr Marinelli reference point)
      ε(1588.20)  ≈ 9.06e-3   (linear interp K-40 1460.8 → Ra-226 1764.5)
      A_true ≈ 2882 Bq (Ac-228 weighted-mean before fix's 911-only
                       contribution; sample mass 1.6 kg → 1802 Bq/kg)
    """
    t_live = 11359.164

    # Custom non-flat ε: ε at 911 ≈ 1.35e-2, at 1597 ≈ 9.06e-3.
    # log-log poly degree-1 fit through these two points.
    # log10 ε = a + b · log10 E
    # at 911:  log10 1.346e-2 = -1.871
    # at 1597: log10 9.06e-3  = -2.043
    # log10 E:  911 → 2.960, 1597 → 3.203
    # b = (-2.043 - (-1.871)) / (3.203 - 2.960) = -0.708
    # a = -1.871 - (-0.708)·2.960 = 0.225
    # ε(E) = 10^(0.225 - 0.708·log10 E)
    # Convert to natural-log polynomial as EfficiencyCurve expects.
    import math as _m
    a10 = 0.225
    b10 = -0.708
    ln_a = a10 * _m.log(10.0)
    # natural-log polynomial: ln ε = a_ln + b_ln · ln E
    # ln ε = ln(10^a10) + ln(10^(b10·log10 E))
    #      = a10·ln10 + b10·ln E
    # so coefficients in ascending power of ln E are
    # (a10·ln10, b10).
    eps_curve = EfficiencyCurve(
        coefficients=(ln_a, b10),
        E_min_keV=100.0,
        E_max_keV=3000.0,
        chi2_per_dof=1.0,
        n_points_used=10,
        n_dof=9,
        detector_id="test",
        geometry="Marinelli",
    )
    # Sanity: ε(911) ≈ 1.35e-2, ε(1597) ≈ 9.06e-3
    eps_911 = eps_curve.efficiency_at(911.0)
    eps_1597 = eps_curve.efficiency_at(1597.0)
    assert 1.2e-2 < eps_911 < 1.5e-2, f"ε(911)={eps_911:.3e}"
    assert 7e-3 < eps_1597 < 1e-2, f"ε(1597)={eps_1597:.3e}"

    # Build a synthetic Ac-228 source where A_true is set to give
    # S(911.20) ≈ 116 355 counts (matches v1.18.31 demo). Then the
    # 1588.20-blob area follows from the same A_true.
    # S = A · ε · I_decimal · t
    # A = S / (ε · I · t)
    A_true = 116355.0 / (eps_911 * 0.258 * t_live)
    # Synthesize the 1588-blob S as if it contained photons from
    # Ac-228 1580.53 + 1588.20 + 1630.63 + Bi-212 1620.50 (chain
    # sibling at equilibrium). Per the v1.18.31 demo:
    #   S_blob (winner+1580 sub-Gaussian) = 20660.92
    #   S_blob (1630 sub-Gaussian)        = 10266.29
    # We model only the main sub-Gaussian (S=20660) — that's the one
    # the survivor 1588.20 takes after BUG-15 dedup.
    S_1597 = 20660.92

    # 911.20: characteristic, isolated peak, no dedup partners
    m_char = _FakeMatch(
        nuclide="Ac-228", library_E_keV=911.204, library_I_pct=25.8,
        peak_channel=313, peak_E_keV=914.02,
        peak_sigma=274.5 / 2.35,  # convert FWHM → σ
        peak_area=116354.79,
        peak_area_uncertainty=2078.88,
        is_characteristic=True,
        peak_area_source="deconvolved_coupled",
    )
    # 1580.53, 1588.20, 1630.63: all on peak_channel=539, all sharing
    # the same S_blob (BUG-15 group). 1588.20 wins (highest I_pct).
    m_1580 = _FakeMatch(
        nuclide="Ac-228", library_E_keV=1580.53, library_I_pct=0.6,
        peak_channel=539, peak_E_keV=1597.35,
        peak_sigma=115.53 / 2.35,
        peak_area=S_1597,
        peak_area_uncertainty=555.15,
        peak_area_source="deconvolved_coupled",
    )
    m_1588 = _FakeMatch(
        nuclide="Ac-228", library_E_keV=1588.20, library_I_pct=3.22,
        peak_channel=539, peak_E_keV=1597.35,
        peak_sigma=115.53 / 2.35,
        peak_area=S_1597,
        peak_area_uncertainty=555.15,
        peak_area_source="deconvolved_coupled",
    )
    m_1630 = _FakeMatch(
        nuclide="Ac-228", library_E_keV=1630.627, library_I_pct=1.51,
        peak_channel=539, peak_E_keV=1597.35,
        peak_sigma=115.53 / 2.35,
        peak_area=S_1597,
        peak_area_uncertainty=555.15,
        peak_area_source="deconvolved_coupled",
    )

    nid = _FakeNuclideId(
        nuclide="Ac-228",
        characteristic_line_keV=911.204,
        matched_lines=(m_char, m_1580, m_1588, m_1630),
    )

    res = compute_activity(
        nid, efficiency_curve=eps_curve,
        live_time_s=t_live, from_bg_subtracted=True,
        decay_correction=False,
    )

    # Two lines used (911 + 1588 winner of dedup); 1580/1630 skipped.
    assert res.n_lines_used() == 2, (
        f"expected 2 lines (911 + 1588 dedup winner), got "
        f"{res.n_lines_used()}; "
        f"used={[la.E_keV for la in res.lines_used]}"
    )

    # Identify the 1588-line in lines_used and verify its A_i is
    # within ±15% of the 911-line A_i.
    la_911 = next(la for la in res.lines_used
                  if abs(la.E_keV - 911.204) < 0.5)
    la_1588 = next(la for la in res.lines_used
                   if abs(la.E_keV - 1588.20) < 0.5)
    ratio = la_1588.A_Bq / la_911.A_Bq
    assert 0.85 <= ratio <= 1.15, (
        f"BUG-27 not closed: A_i(1588)/A_i(911) = {ratio:.3f}, "
        f"expected within ±15% of 1.0 "
        f"(A_i(911)={la_911.A_Bq:.0f}, A_i(1588)={la_1588.A_Bq:.0f})"
    )

    # Sanity: weighted-mean activity also within ±15% of A_true.
    rel_err = abs(res.A_Bq - A_true) / A_true
    assert rel_err < 0.15, (
        f"weighted-mean activity drift {rel_err*100:.1f}% > 15%; "
        f"A_Bq={res.A_Bq:.0f}, A_true={A_true:.0f}"
    )

    # Cross-check: without the BUG-27 fix, A_i(1588) ≈ 3881 Bq/kg in
    # the v1.18.31 demo (ratio 2.15× to 911). Verify we did NOT just
    # silently drop the line.
    assert la_1588.A_Bq > 0


def test_no_dedup_no_boost_single_line_unchanged():
    """If a line has no dedup partners (no other matched line on
    the same peak_channel) AND its peak_area_source is not in the
    chain-sibling-eligible set, the A_i is unchanged from the bare
    library_I_pct formula. Guarantees backward compat for isolated
    peaks (Cs-137 661.66 keV, K-40 1460.82 keV, etc.)."""
    t_live = 1000.0
    S = 50_000.0
    eps = 1.0e-3
    m = _FakeMatch(
        nuclide="Cs-137", library_E_keV=661.66, library_I_pct=85.1,
        peak_channel=300, peak_E_keV=661.5,
        peak_sigma=30.0 / 2.35,  # 30 keV FWHM, well below 50 threshold
        peak_area=S,
        peak_area_uncertainty=math.sqrt(S),
        is_characteristic=True,
        peak_area_source="cowell",  # not in BUG-27 sibling sources
    )
    nid = _FakeNuclideId(
        nuclide="Cs-137",
        characteristic_line_keV=661.66,
        matched_lines=(m,),
    )
    res = compute_activity(
        nid, efficiency_curve=_flat_efficiency_curve(eps),
        live_time_s=t_live, from_bg_subtracted=True,
        decay_correction=False,
    )
    A_expected = S / (eps * 0.851 * t_live)
    assert math.isclose(res.A_Bq, A_expected, rel_tol=1e-6), (
        f"single-line Cs-137 661.66 drifted by BUG-27: "
        f"A_Bq={res.A_Bq:.3f}, expected {A_expected:.3f}"
    )


def test_narrow_peak_no_chain_sibling_boost():
    """FWHM < 50 keV ⇒ sibling boost does not fire even on chain
    nuclides (narrow peaks already isolate the line)."""
    t_live = 1000.0
    S = 100_000.0
    eps = 1.0e-3
    # Ac-228 911.20 at 36 keV FWHM (typical Marinelli σ ≈ 15) →
    # no sibling boost (narrow), no dedup partners → bare-I formula.
    m = _FakeMatch(
        nuclide="Ac-228", library_E_keV=911.204, library_I_pct=25.8,
        peak_channel=313, peak_E_keV=914.0,
        peak_sigma=36.0 / 2.35,  # 36 keV FWHM
        peak_area=S,
        peak_area_uncertainty=math.sqrt(S),
        is_characteristic=True,
        peak_area_source="deconvolved",
    )
    nid = _FakeNuclideId(
        nuclide="Ac-228",
        characteristic_line_keV=911.204,
        matched_lines=(m,),
    )
    res = compute_activity(
        nid, efficiency_curve=_flat_efficiency_curve(eps),
        live_time_s=t_live, from_bg_subtracted=True,
        decay_correction=False,
    )
    A_bare = S / (eps * 0.258 * t_live)
    # No dedup → I_eff = I_pct = 25.8.
    # No sibling boost: 904.20 (I=0.77%) IS within ±18 keV of 911 but
    # it's same-nuclide (excluded). All other Ac-228 lines outside.
    assert math.isclose(res.A_Bq, A_bare, rel_tol=1e-6), (
        f"narrow Ac-228 911 drifted by BUG-27: "
        f"A_Bq={res.A_Bq:.3f}, expected {A_bare:.3f}"
    )


if __name__ == "__main__":
    test_chain_sibling_I_in_window_finds_bi212_for_ac228_1588_blob()
    test_chain_sibling_lookup_skips_non_chain_nuclides()
    test_chain_sibling_lookup_excludes_same_nuclide()
    test_ac228_1588_outlier_converges_to_nuclide_mean_within_15pct()
    test_no_dedup_no_boost_single_line_unchanged()
    test_narrow_peak_no_chain_sibling_boost()
    print("BUG-27 Ac-228 1588 outlier tests PASS.")
