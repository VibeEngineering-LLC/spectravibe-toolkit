"""BUG-TCS#47 — detector_id kwarg drift in compute.py:559-560.

Bug
---
`compute.py:559-560` called `total_efficiency_from_fep(eps_j, E_j,
detector_id=tcs_detector_id, ...)`, but the callee's signature is
`total_efficiency_from_fep(eps_fep, E_keV, crystal_diameter_mm=63.0)` —
no `detector_id` parameter. The kwarg mismatch raised `TypeError`,
which was silently swallowed by the DEEP-06 handler at compute.py:562,
returning 0.0 → TCS correction skipped for every cascade nuclide when
`enable_tcs_correction=True`.

Fix path (Q2 censor plank — 3rd path)
-------------------------------------
The call site now uses `pt_ratio_for_detector(E_keV, tcs_detector_id)`
(public API), which **raises `KeyError`** on unknown detector_id. The
`KeyError` is caught by the generic DEEP-06 silent handler
(`try/except Exception as exc:` at compute.py:558-562). On catch,
`ε_T` is set to `0.0`, so `C_TCS` evaluates to `1.0` (effective
TCS-skip), and a `UserWarning` is emitted with F296 text. The 63 mm
Gamma-1S default is NOT substituted.

Q4 censor plank — stub geometry
-------------------------------
`is_significant=True` requires correction ≥ 5% (per
tcs_close_geometry.py:184/192/236). The stub uses `ε_FEP = 0.05`
constant (Co-60 Marinelli-close geometry — censor-suggested geometry
class). For Gamma-1S 63 mm at 1332.5 keV:
    P/T(1332, 63 mm) ≈ pt_3in3(1332) - 0.02·(76-63)/10
                     ≈ 0.294 - 0.026 = 0.268
    ε_T(1332)         = 0.05 / 0.268 ≈ 0.187
    sum_L (Co-60 1173 → 1332) = 0.998 × 0.187 ≈ 0.186
    C_TCS             = 1/(1 − 0.186) ≈ 1.229 → +22.9% ≥ 5% ✓

Red-then-green design (Gate C)
------------------------------
Before fix: `_eps_T` raises TypeError → DEEP-06 catches → returns 0.0
→ `sum_L = 0` → `C_TCS = 1.0` → `is_significant=False` →
`tcs_auto_applied` stays empty → the note block at compute.py:1094-1098
is never entered → the "F-296 / v1.18.1: auto-TCS correction applied"
substring is absent from `r.notes`. The note assertion fails. RED.

After fix: ε_T returned correctly → `is_significant=True` (correction
≈ 22.9%) → `tcs_auto_applied` populated → note emitted. GREEN.

Three censor-plank Q4 assertions are checked (note string, non-empty
tcs flag via `coincidence_correction_applied`, non-trivial Δactivity
vs TCS-off baseline).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

# Ensure scripts/ is on sys.path (matches sibling tests' convention).
_REPO_ROOT = Path(__file__).resolve().parents[2]
_SCRIPTS = _REPO_ROOT / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from gamma.activity.compute import compute_activity  # noqa: E402


class _StubEff:
    """ε_FEP(E) = 0.05 constant — Co-60 Marinelli-close geometry-class.

    With Gamma-1S 63 mm crystal, this yields correction > 5% (≈ 22.9%
    on the Co-60 1173/1332 cascade) — well above the
    `significant_threshold_pct=5.0` gate at
    tcs_close_geometry.py:192/236.
    """

    def efficiency_at(self, E):
        return 0.05

    def is_extrapolating(self, E, margin_factor=1.1):
        return False


class _MockLine:
    def __init__(self, E, I, area, area_unc, src="cowell"):
        self.library_E_keV = E
        self.library_I_pct = I        # canonical field name in compute.py
        self.I_pct = I                 # legacy alias (some helpers read it)
        self.peak_area = area
        self.peak_area_uncertainty = area_unc   # canonical name
        self.peak_area_unc = area_unc           # legacy alias
        self.peak_area_source = src
        self.peak_channel = None
        self.blob_partner_energies_keV = ()


class _MockNuclideId:
    def __init__(self, name, lines):
        self.nuclide = name
        self.detected = True
        self.matched_lines = lines


def _mk_co60():
    return _MockNuclideId(
        "Co-60",
        [
            _MockLine(1173.23, 99.85, 5000.0, math.sqrt(5000.0)),
            _MockLine(1332.49, 99.98, 4800.0, math.sqrt(4800.0)),
        ],
    )


def test_tcs47_correction_actually_applied():
    """Three-assertion censor-Q4 gate: note string + tcs flag + Δactivity.

    All three must hold simultaneously after fix:
      (i)  "F-296/v1.18.1: auto-TCS correction applied" substring in notes
      (ii) `coincidence_correction_applied` is True (tcs_auto_applied
           was non-empty)
      (iii) post-TCS activity differs non-trivially from TCS-off baseline
           (ratio ≠ 1.0 by > 1% — Co-60 close-geometry expected ~20-25%).

    Before fix: all three FAIL (silent skip → no note, flag stays False,
    activity equal to TCS-off).
    """
    co60_on = _mk_co60()
    co60_off = _mk_co60()

    r_on = compute_activity(
        co60_on,
        efficiency_curve=_StubEff(),
        live_time_s=600.0,
        from_bg_subtracted=True,
        enable_tcs_correction=True,
        tcs_detector_id="Gamma-1S",
    )
    r_off = compute_activity(
        co60_off,
        efficiency_curve=_StubEff(),
        live_time_s=600.0,
        from_bg_subtracted=True,
        enable_tcs_correction=False,
        tcs_detector_id="Gamma-1S",
    )

    notes_on = r_on.notes or ""

    # (i) — note string (censor plank Q4 verbatim substring)
    assert "auto-TCS correction applied" in notes_on, (
        f"Expected 'auto-TCS correction applied' in notes, got: "
        f"{notes_on!r}. TCS was silently skipped — BUG-TCS#47 kwarg "
        f"drift still active in _eps_T closure (compute.py:559-560)."
    )

    # (ii) — tcs_auto_applied was populated (proxy: the
    # coincidence_correction_applied flag is True only when auto-TCS
    # filled coincidence_correction)
    assert r_on.coincidence_correction_applied is True, (
        f"coincidence_correction_applied=False — tcs_auto_applied "
        f"stayed empty (TCS silently skipped). r_on.notes={notes_on!r}."
    )

    # (iii) — non-trivial Δactivity: TCS-on must differ from TCS-off
    # by > 1% (Co-60 close-geometry ε=0.05 stub → ~22.9% bump).
    assert r_off.A_Bq > 0
    assert r_on.A_Bq > 0
    ratio = r_on.A_Bq / r_off.A_Bq
    assert abs(ratio - 1.0) > 0.01, (
        f"TCS-on/TCS-off activity ratio = {ratio:.6f} — Δ < 1%. "
        f"Expected ~1.15-1.25× for Co-60 ε_FEP=0.05 Marinelli stub. "
        f"TCS was applied per the note string but produced no effective "
        f"correction — investigate compute_tcs_correction inputs."
    )


def test_tcs47_unknown_detector_graceful():
    """Backward-compat with test_v1_18_1.py:162-175.

    Unknown detector_id must NOT crash; activity > 0; tcs_detector_id
    falls back to Gamma-1S 63 mm via a WARN-and-default path (not
    silent — censor Q2 plank: the fallback emits warnings.warn).
    """
    co60 = _MockNuclideId(
        "Co-60",
        [_MockLine(1173.23, 99.85, 5000.0, math.sqrt(5000.0))],
    )
    with pytest.warns(UserWarning, match="unknown.*detector|Unknown detector"):
        r = compute_activity(
            co60,
            efficiency_curve=_StubEff(),
            live_time_s=600.0,
            from_bg_subtracted=True,
            enable_tcs_correction=True,
            tcs_detector_id="UNKNOWN_DETECTOR",
        )
    assert r.A_Bq > 0, (
        "Unknown detector_id must not crash; activity must be > 0 "
        "(backward-compat with test_v1_18_1.py:162-175 "
        "test_F296_tcs_detector_id_passed)."
    )
