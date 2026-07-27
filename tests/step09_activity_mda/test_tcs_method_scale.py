"""
Area-method-aware TCS scaling (K-18 fix, F-35 v1.7.13).

`compute_activity` reads `LineMatch.peak_area_source` and scales the
analytic TCS correction accordingly:

  c_effective = 1 + (c_analytic − 1) · scale[area_source]

Default scale: 0.0 for `"lsrm_peaks_table"` (the Lsrm wide-ROI
Gaussian fit already recovers summing-displaced counts), 1.0 for
everything else (Cowell, deconvolved, unknown — full TCS).

Verifies:
  - Default scale: Lsrm-table source ⇒ TCS effectively disabled.
  - Default scale: Cowell / deconvolved / empty / unknown source ⇒
    full TCS applied unchanged (regression-safe for existing tests).
  - Custom `tcs_method_scale` overrides the defaults; partial scaling
    yields the algebraic mid-point.
  - The `notes` field records when scaling was applied.
  - No TCS dict ⇒ `tcs_method_scale` is irrelevant; activity is
    identical to the no-scaling path.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

from gamma.activity import (
    compute_activity, DEFAULT_TCS_METHOD_SCALE,
)
from gamma.identification.identify import (
    NuclideIdentification, LineMatch,
)


# ---------------------------------------------------------------------------
# Stub efficiency curve — borrowed from test_cascade_summing.py
# ---------------------------------------------------------------------------

class StubEfficiency:
    """ε(E) = a·E^b · 10^-2, mimicking a real NaI 5cm curve.

    Defaults: a=12.6, b=-1.05 → ε(1173)=~1.4e-3.
    """

    def __init__(self, a: float = 12.6, b: float = -1.05):
        self.a, self.b = a, b
        self.E_min = 50.0
        self.E_max = 3000.0

    def efficiency_at(self, E_keV: float) -> float:
        return self.a * (E_keV ** self.b) * 1e-2

    def is_extrapolating(self, E_keV: float) -> bool:
        return E_keV < self.E_min or E_keV > self.E_max


def _co60_line(area_source: str, *, E: float = 1173.23,
               area: float = 1.0e5, sigma: float = 1.0e3) -> LineMatch:
    return LineMatch(
        nuclide="Co-60",
        library_E_keV=E, library_I_pct=99.85,
        peak_channel=400, peak_E_keV=E,
        peak_sigma=50.0, residual_keV=0.0,
        is_characteristic=True,
        peak_area=area, peak_area_uncertainty=sigma,
        peak_area_source=area_source,
    )


def _co60_id(lm: LineMatch) -> NuclideIdentification:
    return NuclideIdentification(
        nuclide="Co-60", detected=True,
        reason="synthetic fixture",
        characteristic_line_keV=lm.library_E_keV,
        matched_lines=(lm,),
    )


def _activity_with(lm: LineMatch, *, tcs_dict=None,
                   tcs_method_scale=None) -> float:
    eff = StubEfficiency()
    return compute_activity(
        _co60_id(lm),
        efficiency_curve=eff, live_time_s=1800.0,
        from_bg_subtracted=True,
        coincidence_correction=tcs_dict,
        tcs_method_scale=tcs_method_scale,
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_default_scale_lsrm_table_kills_tcs():
    """Default scale: TCS dict on a Lsrm-table-sourced line ⇒ effective
    c=1.0 ⇒ activity identical to the no-TCS run."""
    tcs = {1173.23: 1.05}
    r_no = _activity_with(_co60_line("lsrm_peaks_table"))
    r_lsrm = _activity_with(_co60_line("lsrm_peaks_table"), tcs_dict=tcs)
    assert math.isclose(r_no.A_Bq, r_lsrm.A_Bq, rel_tol=1e-9), (
        f"Lsrm-table TCS should be zeroed: {r_no.A_Bq} vs {r_lsrm.A_Bq}"
    )
    assert "K-18: TCS scaled by area-method on 1 line(s)" in r_lsrm.notes, (
        f"missing K-18 note in {r_lsrm.notes!r}"
    )
    print(f"  ✓ test_default_scale_lsrm_table_kills_tcs "
          f"(A_no={r_no.A_Bq:.3e}, A_lsrm={r_lsrm.A_Bq:.3e})")


def test_default_scale_cowell_keeps_full_tcs():
    """Default scale: Cowell-sourced area gets full TCS — activity
    scales by exactly the TCS factor."""
    tcs = {1173.23: 1.05}
    r_no = _activity_with(_co60_line("cowell"))
    r_tcs = _activity_with(_co60_line("cowell"), tcs_dict=tcs)
    ratio = r_tcs.A_Bq / r_no.A_Bq
    assert math.isclose(ratio, 1.05, rel_tol=1e-9), (
        f"Cowell TCS ratio {ratio:.6f} ≠ 1.05"
    )
    assert "K-18: TCS scaled" not in r_tcs.notes, (
        f"K-18 note should NOT appear for full-scale path: {r_tcs.notes!r}"
    )
    print(f"  ✓ test_default_scale_cowell_keeps_full_tcs (ratio={ratio:.4f})")


def test_default_scale_deconvolved_keeps_full_tcs():
    """Deconvolved areas also get full TCS by default — they are
    channel-sum-equivalent (linear LSQ-fit of fixed-position Gaussians
    over a step+linear baseline does not recover summing-displaced
    counts, those land outside the photopeak ROI entirely)."""
    tcs = {1173.23: 1.05}
    r_no = _activity_with(_co60_line("deconvolved"))
    r_tcs = _activity_with(_co60_line("deconvolved"), tcs_dict=tcs)
    ratio = r_tcs.A_Bq / r_no.A_Bq
    assert math.isclose(ratio, 1.05, rel_tol=1e-9), (
        f"Deconvolved TCS ratio {ratio:.6f} ≠ 1.05"
    )
    print(f"  ✓ test_default_scale_deconvolved_keeps_full_tcs (ratio={ratio:.4f})")


def test_default_scale_empty_source_keeps_full_tcs():
    """Empty source (legacy LineMatch from old code paths) ⇒ full TCS
    (safe default — never silently disable a correction)."""
    tcs = {1173.23: 1.05}
    r_no = _activity_with(_co60_line(""))
    r_tcs = _activity_with(_co60_line(""), tcs_dict=tcs)
    ratio = r_tcs.A_Bq / r_no.A_Bq
    assert math.isclose(ratio, 1.05, rel_tol=1e-9), (
        f"Empty-source TCS ratio {ratio:.6f} ≠ 1.05"
    )
    print(f"  ✓ test_default_scale_empty_source_keeps_full_tcs (ratio={ratio:.4f})")


def test_unknown_source_label_default_full_tcs():
    """Unknown source labels (e.g. 'foo') fall back to full TCS."""
    tcs = {1173.23: 1.05}
    r_no = _activity_with(_co60_line("foo"))
    r_tcs = _activity_with(_co60_line("foo"), tcs_dict=tcs)
    ratio = r_tcs.A_Bq / r_no.A_Bq
    assert math.isclose(ratio, 1.05, rel_tol=1e-9), (
        f"Unknown-source TCS ratio {ratio:.6f} ≠ 1.05"
    )
    print(f"  ✓ test_unknown_source_label_default_full_tcs (ratio={ratio:.4f})")


def test_custom_scale_partial_lsrm():
    """Caller can override the default with a partial scale (e.g. 0.5
    for Lsrm in a geometry where the wide-ROI fit recovers only half
    the summing-displaced counts)."""
    tcs = {1173.23: 1.10}
    custom = {"lsrm_peaks_table": 0.5}
    r_no = _activity_with(_co60_line("lsrm_peaks_table"))
    r_tcs = _activity_with(
        _co60_line("lsrm_peaks_table"),
        tcs_dict=tcs, tcs_method_scale=custom,
    )
    # Effective c = 1 + (1.10 − 1) · 0.5 = 1.05
    expected_ratio = 1.05
    ratio = r_tcs.A_Bq / r_no.A_Bq
    assert math.isclose(ratio, expected_ratio, rel_tol=1e-9), (
        f"partial Lsrm scale: ratio {ratio:.6f} ≠ expected {expected_ratio}"
    )
    print(f"  ✓ test_custom_scale_partial_lsrm "
          f"(ratio={ratio:.4f}, expected={expected_ratio})")


def test_custom_scale_overrides_only_named_keys():
    """The user dict is merged ON TOP of the defaults, not replacing
    them entirely — an entry only for one source label leaves the
    other defaults intact."""
    tcs = {1173.23: 1.05}
    # Only override 'cowell' to 0.0; lsrm_peaks_table keeps its default 0.0
    custom = {"cowell": 0.0}
    r_cowell = _activity_with(_co60_line("cowell"), tcs_dict=tcs,
                              tcs_method_scale=custom)
    r_no = _activity_with(_co60_line("cowell"))
    assert math.isclose(r_no.A_Bq, r_cowell.A_Bq, rel_tol=1e-9), (
        f"caller override cowell→0 should match no-TCS: "
        f"{r_no.A_Bq} vs {r_cowell.A_Bq}"
    )
    # lsrm_peaks_table default (0.0) should still apply
    r_lsrm = _activity_with(_co60_line("lsrm_peaks_table"), tcs_dict=tcs,
                            tcs_method_scale=custom)
    assert math.isclose(r_no.A_Bq, r_lsrm.A_Bq, rel_tol=1e-9), (
        f"lsrm_peaks_table default 0.0 should survive partial override"
    )
    print(f"  ✓ test_custom_scale_overrides_only_named_keys")


def test_no_tcs_dict_no_scaling_applied():
    """When `coincidence_correction` is None, `tcs_method_scale` is
    inert — activity unchanged and no K-18 note is emitted."""
    r_a = _activity_with(_co60_line("lsrm_peaks_table"))
    r_b = _activity_with(_co60_line("lsrm_peaks_table"),
                         tcs_method_scale={"lsrm_peaks_table": 0.5})
    assert math.isclose(r_a.A_Bq, r_b.A_Bq, rel_tol=1e-9), (
        "scale parameter must be inert without a coincidence_correction"
    )
    assert "K-18: TCS scaled" not in r_a.notes
    assert "K-18: TCS scaled" not in r_b.notes
    print(f"  ✓ test_no_tcs_dict_no_scaling_applied")


def test_default_dict_publishes_canonical_keys():
    """Sanity: `DEFAULT_TCS_METHOD_SCALE` documents every source label
    used elsewhere in the codebase."""
    for key in ("", "cowell", "deconvolved", "failed", "lsrm_peaks_table"):
        assert key in DEFAULT_TCS_METHOD_SCALE, (
            f"missing canonical key {key!r} in DEFAULT_TCS_METHOD_SCALE"
        )
    # Sanity on values
    assert DEFAULT_TCS_METHOD_SCALE["lsrm_peaks_table"] == 0.0
    assert DEFAULT_TCS_METHOD_SCALE["cowell"] == 1.0
    assert DEFAULT_TCS_METHOD_SCALE[""] == 1.0
    print(f"  ✓ test_default_dict_publishes_canonical_keys "
          f"(keys={sorted(DEFAULT_TCS_METHOD_SCALE)})")


def test_mixed_sources_in_one_nuclide():
    """A nuclide with two lines from different methods: each line
    gets its own scaling. Co-60 with 1173 from Lsrm-table and 1332
    from Cowell — only the 1332 contribution receives the full TCS."""
    eff = StubEfficiency()
    lm1 = _co60_line("lsrm_peaks_table", E=1173.23, area=1.0e5)
    lm2 = LineMatch(
        nuclide="Co-60",
        library_E_keV=1332.49, library_I_pct=99.98,
        peak_channel=420, peak_E_keV=1332.49,
        peak_sigma=50.0, residual_keV=0.0,
        is_characteristic=False,
        peak_area=1.0e5, peak_area_uncertainty=1.0e3,
        peak_area_source="cowell",
    )
    ni = NuclideIdentification(
        nuclide="Co-60", detected=True,
        reason="mixed-source test", characteristic_line_keV=1173.23,
        matched_lines=(lm1, lm2),
    )
    tcs = {1173.23: 1.05, 1332.49: 1.05}
    res = compute_activity(
        ni, efficiency_curve=eff, live_time_s=1800.0,
        from_bg_subtracted=True, coincidence_correction=tcs,
    )
    # Per-line: lm1 (Lsrm) has c_eff = 1.0; lm2 (Cowell) has c_eff = 1.05
    line1 = next(la for la in res.lines_used if abs(la.E_keV - 1173.23) < 0.5)
    line2 = next(la for la in res.lines_used if abs(la.E_keV - 1332.49) < 0.5)
    assert math.isclose(line1.correction_factor, 1.0, rel_tol=1e-9), (
        f"Lsrm line c={line1.correction_factor} should be 1.0"
    )
    assert math.isclose(line2.correction_factor, 1.05, rel_tol=1e-9), (
        f"Cowell line c={line2.correction_factor} should be 1.05"
    )
    assert "K-18: TCS scaled by area-method on 1 line(s)" in res.notes, (
        f"note should say exactly 1 line was scaled, got: {res.notes!r}"
    )
    print(f"  ✓ test_mixed_sources_in_one_nuclide "
          f"(line1.c={line1.correction_factor:.3f}, "
          f"line2.c={line2.correction_factor:.3f})")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("Running K-18 TCS-method-scale tests...\n")
    test_default_scale_lsrm_table_kills_tcs()
    test_default_scale_cowell_keeps_full_tcs()
    test_default_scale_deconvolved_keeps_full_tcs()
    test_default_scale_empty_source_keeps_full_tcs()
    test_unknown_source_label_default_full_tcs()
    test_custom_scale_partial_lsrm()
    test_custom_scale_overrides_only_named_keys()
    test_no_tcs_dict_no_scaling_applied()
    test_default_dict_publishes_canonical_keys()
    test_mixed_sources_in_one_nuclide()
    print("\nAll K-18 tests passed.")
