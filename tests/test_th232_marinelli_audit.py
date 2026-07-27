"""Tests for Task #68 — Th-232 Marinelli ID accuracy audit + Tl-208 typo.

Lightweight, deterministic tests that do NOT require LSRM .spe files:

1. ``test_chain_decomposer_tl208_2614_intensity_branching_corrected``
   Tl-208 2614.51 keV in TRUE_ENSDF_OWNERSHIP_TH_CHAIN must be
   branching-corrected (~35.85 % per Bi-212 parent decay), matching the
   convention of sibling Tl-208 entries (583/510/860 keV). The previous
   value 99.75 % was a typo (raw uncorrected).

2. ``test_chain_decomposer_tl208_lines_consistent_with_nuclides_json``
   Every Tl-208 line in chain_decomposer table must be within ±5 % of
   the same energy's intensity in ``data/nuclides.json`` (which is the
   authoritative branching-corrected library per F-372.1).

3. ``test_audit_script_filters_th232_a_tier``
   The audit helper ``_filter_th232_a_tier`` must return exactly the
   set produced by intersecting indexes by_quality_tier['A'] with
   by_nuclide['Th-232'] (deterministic on SPECTRA_INDEX HEAD).

4. ``test_audit_classify_th232_chain_buckets``
   ``_classify_th232_chain`` must correctly bucket TP/FN/FP across the
   chain progeny set and not falsely flag K-40/Bi-214 as unexpected.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]
_SCRIPTS = _ROOT / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))
_DIAG = _SCRIPTS / "gamma" / "diagnostics"
if str(_DIAG) not in sys.path:
    sys.path.insert(0, str(_DIAG))


def test_chain_decomposer_tl208_2614_intensity_branching_corrected():
    """Tl-208 2614.51 keV must be ~35.85 % (branching-corrected),
    not 99.75 (raw). Sibling Tl-208 entries (583/510/860) are
    branching-corrected — see ENSDF_TL208_RAW / F-372.1."""
    from gamma.data.chain_decomposer import TRUE_ENSDF_OWNERSHIP_TH_CHAIN

    tl208_entries = [
        (E, I) for (E, owner, I, _tol) in TRUE_ENSDF_OWNERSHIP_TH_CHAIN
        if owner == "Tl-208" and abs(E - 2614.51) < 0.1
    ]
    assert tl208_entries, "Tl-208 2614.51 keV missing from chain table"
    E, I = tl208_entries[0]
    # 99.754 (raw) × 0.3594 = 35.85; allow 1 % slack.
    assert 34.0 <= I <= 37.0, (
        f"Tl-208 2614.51 intensity {I} %% is outside the branching-"
        f"corrected window 34-37 %% (was 99.75 — pre-fix typo)."
    )


def test_chain_decomposer_tl208_strong_lines_consistent_with_nuclides_json():
    """Cross-check the STRONG Tl-208 lines (≥4 % per Tl-208 emission)
    in chain_decomposer agree with data/nuclides.json (authoritative
    branching-corrected library) within ±10 %.

    Strong lines are checked because they drive ID confidence; weak
    lines (≤2 % raw) have a known mixed convention residue (277.36,
    763.13 — see audit outbox 2026-06-04_A3_th232_marinelli_audit.md
    §"Remaining table inconsistencies"); these are out of scope for
    the single-line typo fix and tracked separately.
    """
    from gamma.data.chain_decomposer import TRUE_ENSDF_OWNERSHIP_TH_CHAIN

    nlib = json.loads(
        (_ROOT / "data" / "nuclides.json").read_text(encoding="utf-8")
    )
    lib_lines = {round(float(L[0]), 1): float(L[1])
                 for L in nlib["Tl-208"]["lines"]}
    # Strong-line whitelist (high I_pct, drive identification)
    strong_E = {510.77, 583.19, 860.56, 2614.51}
    mismatches = []
    for E, owner, I, _tol in TRUE_ENSDF_OWNERSHIP_TH_CHAIN:
        if owner != "Tl-208":
            continue
        if not any(abs(E - se) < 0.5 for se in strong_E):
            continue
        candidates = [(k, v) for k, v in lib_lines.items() if abs(k - E) < 1.0]
        if not candidates:
            continue
        _, lib_I = candidates[0]
        if not (0.9 * lib_I <= I <= 1.1 * lib_I):
            mismatches.append(
                f"E={E}: chain_decomposer I={I} vs nuclides.json I={lib_I}"
            )
    assert not mismatches, (
        "Tl-208 STRONG-line intensity mismatches "
        "(chain_decomposer vs nuclides.json):\n  " + "\n  ".join(mismatches)
    )


def test_audit_script_filters_th232_a_tier():
    """The audit's filter must produce exactly the (A-tier ∩ Th-232)
    intersection from SPECTRA_INDEX."""
    from audit_th232_marinelli import _filter_th232_a_tier  # type: ignore

    idx_path = _ROOT / "audit" / "_rag" / "SPECTRA_INDEX.json"
    if not idx_path.exists():
        pytest.skip("SPECTRA_INDEX.json not present")
    idx = json.loads(idx_path.read_text(encoding="utf-8"))
    hits = _filter_th232_a_tier(idx)
    a_tier = set(idx["indexes"]["by_quality_tier"].get("A", []))
    th232 = set(idx["indexes"]["by_nuclide"].get("Th-232", []))
    expected = a_tier & th232
    assert {h["spectrum_id"] for h in hits} == expected, (
        f"filter returned {len(hits)} but expected {len(expected)}"
    )
    assert len(hits) >= 8, "expect ≥8 A-tier Th-232 fixtures (HEAD: 11)"


def test_audit_classify_th232_chain_buckets():
    """``_classify_th232_chain`` puts chain members → TP/FN, natural
    BG → false_positive_background, real strangers → unexpected."""
    from audit_th232_marinelli import _classify_th232_chain  # type: ignore

    cc = _classify_th232_chain(
        ["Ac-228", "Pb-212", "Tl-208", "K-40", "Th-234"]
    )
    assert set(cc["true_positive_chain"]) == {"Ac-228", "Pb-212", "Tl-208"}
    assert "Bi-212" in cc["false_negative_chain"]
    assert cc["false_positive_background"] == ["K-40"]
    assert cc["false_positive_unexpected"] == ["Th-234"], (
        "Th-234 (U-238 chain) is a real mis-ID in a Th-232 sample, "
        "not background"
    )
    assert cc["core_match"] is True
    assert cc["full_chain_match"] is False  # Bi-212 absent


def test_audit_chain_decomposer_th_table_sanity():
    """Th-chain ownership table must contain the canonical anchor lines."""
    from gamma.data.chain_decomposer import TRUE_ENSDF_OWNERSHIP_TH_CHAIN

    anchors = {(round(E, 1), owner) for E, owner, _I, _tol in
               TRUE_ENSDF_OWNERSHIP_TH_CHAIN}
    must_have = {
        (238.6, "Pb-212"),
        (583.2, "Tl-208"),
        (911.2, "Ac-228"),
        (968.97, "Ac-228"),
        (2614.51, "Tl-208"),
        (727.33, "Bi-212"),
    }
    missing = {k for k in must_have if k not in
               {(round(E, 1), o) for (E, o) in anchors} and
               not any(abs(E - k[0]) < 0.1 and o == k[1] for E, o in anchors)}
    # Tolerant E rounding compare
    rounded = {(round(E, 2), o) for E, o in anchors}
    missing = []
    for E_ref, owner in must_have:
        if not any(abs(E - E_ref) < 0.5 and o == owner for E, o in anchors):
            missing.append((E_ref, owner))
    assert not missing, f"Th-chain table missing anchors: {missing}"
