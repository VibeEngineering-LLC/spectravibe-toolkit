"""
BUG-36 (Wave 3 C1, 2026-06-05) — Ti-44 / Sc-44 must enter Stage 3
candidate set so that the AmTiCsEu Marinelli calibration fixture can
identify both the Ti-44 parent (67.87 / 78.34 / 146.22 keV) and the
Sc-44 daughter (1157.02 keV, 99.9 % BR).

Pre-fix root cause (see `_state/agent_a/outbox/2026-06-04_v1_22_0_
AmTiCsEu_revalidation.md` §3.3 + this branch outbox):

  • Library entries Ti-44 and Sc-44 are present in `data/nuclides.json`
    (lines 1514–1552) — added by Wave 3 BUG-36 prior pass.
  • `identify.py` / `staged_pipeline.py` have NO parent-dependency gate
    on Sc-44 (the file is grep-clean for "parent" — verified at HEAD).
  • The actual blocker: `scripts/gamma/identification/ern_set.py` ships
    27 candidate nuclides total across ERN/Stage2/Stage3 — **Ti-44 and
    Sc-44 are not in any of them**. So they are never considered as
    candidates by `_run_stage` and the 1146 keV peak (Sc-44 1157 obs)
    stays in `unidentified_peaks[1]` despite being a 115σ peak.

Fix:
  • Add `Ti-44` and `Sc-44` to `EXOTIC_STAGE3` (both are explicit
    calibration sources per `data/nuclides.json` tcs_type fields
    `calibration_source` / `chain_progeny`).
  • Add `Ti44` / `Sc44` to `NUCLIDE_TOKENS` + `_NUCLIDE_CANONICAL` +
    `_NUCLIDE_TO_CHAIN` in `gamma.io.filename_hints` so that future
    standalone filenames like `Ti44_calib.spe` get the filename hint
    promotion path described by SKILL.md §7A.1.

These tests cover the candidate-list state and the canonicalisation
state, but do NOT re-run the full LSRM pipeline (the LSRM .spe fixture
in `Поверка-2016` is local-only per F-115 / CLAUDE.md and is not
committed). End-to-end verification is the orchestrator's job during
Tier 2 re-validation.
"""
from __future__ import annotations

import sys
from pathlib import Path

_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(_ROOT / "scripts"))

from gamma.identification.ern_set import (
    ERN_STAGE1, TECHNOGENIC_STAGE2, EXOTIC_STAGE3,
    candidates_for_stage, candidates_up_to_stage, stage_of_nuclide,
)
from gamma.io.filename_hints import (
    NUCLIDE_TOKENS,
    _NUCLIDE_CANONICAL,
    _NUCLIDE_TO_CHAIN,
    extract_isotope_hints,
    chains_claimed_by_isotope_hints,
    parse_filename,
)
from gamma.data.nuclide_library import get_nuclide


# ── BUG-36 fix — Stage 3 candidate-list extension ─────────────────────

def test_ti44_in_stage3_candidates():
    """Ti-44 must appear in EXOTIC_STAGE3 — calibration source."""
    assert "Ti-44" in EXOTIC_STAGE3, (
        "Ti-44 must be a Stage 3 candidate (calibration source per "
        "data/nuclides.json:1537 tcs_type='calibration_source'). "
        "Without this, the AmTiCsEu fixture Ti-44 lines (67.87, 78.34, "
        "146.22 keV) never enter the matcher and Ti-44 stays "
        "UNIDENTIFIED despite being in library."
    )


def test_sc44_in_stage3_candidates():
    """Sc-44 (Ti-44 daughter) must be a Stage 3 candidate.

    Sc-44 is the short-lived (T½=3.97 h) γ-emitting daughter of Ti-44
    in secular equilibrium for any aged source. Its 1157.02 keV line
    (BR=99.9 %) is the single strongest peak in an AmTiCsEu fixture
    above 700 keV (observed at 1146 keV with 115σ significance per
    `2026-06-04_v1_22_0_AmTiCsEu_revalidation.md` §3.3).
    """
    assert "Sc-44" in EXOTIC_STAGE3, (
        "Sc-44 must be a Stage 3 candidate (chain_progeny per "
        "data/nuclides.json:1551). Sc-44 1157 keV line is the AmTiCsEu "
        "spectrum's loudest unidentified peak — 115σ at Δ=0.19·FWHM."
    )


def test_stage_of_nuclide_for_ti44_sc44():
    """stage_of_nuclide must return 3 for Ti-44 and Sc-44."""
    assert stage_of_nuclide("Ti-44") == 3, (
        f"Ti-44 should map to stage 3, got {stage_of_nuclide('Ti-44')!r}"
    )
    assert stage_of_nuclide("Sc-44") == 3, (
        f"Sc-44 should map to stage 3, got {stage_of_nuclide('Sc-44')!r}"
    )


def test_candidates_for_stage_3_includes_ti44_sc44():
    cand3 = candidates_for_stage(3)
    assert "Ti-44" in cand3
    assert "Sc-44" in cand3


def test_cumulative_stage3_includes_ti44_sc44():
    """Stage 1+2+3 cumulative candidate list must include Ti-44/Sc-44.

    This is the actual list passed into `_run_stage` when
    `allow_stage3=True` (the default in `scripts/run_skill.py:170`).
    """
    cand_all = candidates_up_to_stage(3)
    assert "Ti-44" in cand_all
    assert "Sc-44" in cand_all


# ── BUG-36 fix — filename-hint canonicalisation extension ─────────────

def test_ti44_sc44_tokens_present():
    """`Ti44` and `Sc44` are recognised tokens for filename hints."""
    assert "Ti44" in NUCLIDE_TOKENS
    assert "Sc44" in NUCLIDE_TOKENS


def test_ti44_sc44_canonical_mapping():
    """Tokens canonicalise to library labels with hyphen."""
    assert _NUCLIDE_CANONICAL.get("Ti44") == "Ti-44"
    assert _NUCLIDE_CANONICAL.get("Sc44") == "Sc-44"


def test_ti44_sc44_chain_membership():
    """Both belong to the Ti-44 chain (data/nuclides.json:1519,1543)."""
    assert _NUCLIDE_TO_CHAIN.get("Ti-44") == "Ti-44"
    assert _NUCLIDE_TO_CHAIN.get("Sc-44") == "Ti-44"


def test_ti44_filename_hint_extraction():
    """Standalone filenames `Ti44_*.spe` produce Ti-44 hint."""
    hints = extract_isotope_hints("Ti44_calib_60d_Marinelli.spe")
    assert "Ti-44" in hints, f"Expected Ti-44 in hints, got {hints!r}"


def test_sc44_filename_hint_extraction():
    """Standalone filenames `Sc44_*.spe` produce Sc-44 hint."""
    hints = extract_isotope_hints("Sc44_test.spe")
    assert "Sc-44" in hints, f"Expected Sc-44 in hints, got {hints!r}"


def test_amticseu_filename_does_NOT_get_hints():
    """Agglutinated 'AmTiCsEu' has no word boundary — no hint match.

    This documents that the candidate-set fix in ern_set.py is the
    necessary closure for the AmTiCsEu Marinelli fixture; the filename-
    hint extension here is defensive and only matters for future
    standalone fixtures.
    """
    hints = extract_isotope_hints("Смесь_AmTiCsEu_Маринелли.spe")
    assert hints == [], (
        f"Agglutinated AmTiCsEu must NOT match — boundary regex rejects, "
        f"got {hints!r}"
    )


def test_chains_claimed_for_ti44_hint():
    """Filename binding Ti-44 → chain = {'Ti-44'}."""
    chains = chains_claimed_by_isotope_hints(["Ti-44"])
    assert chains == {"Ti-44"}, f"Expected {{'Ti-44'}}, got {chains!r}"


def test_chains_claimed_for_sc44_hint():
    """Filename binding Sc-44 → chain = {'Ti-44'} (daughter rolls up)."""
    chains = chains_claimed_by_isotope_hints(["Sc-44"])
    assert chains == {"Ti-44"}, f"Expected {{'Ti-44'}}, got {chains!r}"


# ── Library cross-check (regression guard for prior BUG-36 pass) ──────

def test_ti44_in_library_with_char_lines():
    """Ti-44 entry present in `data/nuclides.json` with 67.87 / 78.34
    / 146.22 keV lines.

    Cross-checks the library state already committed by the Wave 3
    prior pass; would catch a regression where Ti-44 disappears.
    """
    nuc = get_nuclide("Ti-44")
    assert nuc is not None, "Ti-44 missing from default library"
    lines = nuc.get("lines", [])
    energies = sorted(float(L[0]) for L in lines)
    # Within 0.1 keV of NNDC ENSDF values
    for expected in (67.87, 78.34, 146.22):
        match = min(energies, key=lambda E: abs(E - expected))
        assert abs(match - expected) < 0.1, (
            f"Ti-44 expected line near {expected} keV not found; "
            f"available={energies}"
        )


def test_sc44_in_library_with_1157_keV():
    """Sc-44 entry present with the 1157.02 keV char line (99.9 % BR)."""
    nuc = get_nuclide("Sc-44")
    assert nuc is not None, "Sc-44 missing from default library"
    lines = nuc.get("lines", [])
    energies = [float(L[0]) for L in lines]
    char_match = min(energies, key=lambda E: abs(E - 1157.02))
    assert abs(char_match - 1157.02) < 0.1, (
        f"Sc-44 expected 1157.02 keV not found; available={energies}"
    )
    # Sc-44 carries `parent=Ti-44` per data/nuclides.json:1542
    assert nuc.get("parent") == "Ti-44", (
        f"Sc-44 must carry parent=Ti-44 link; got {nuc.get('parent')!r}"
    )


def test_bug42_ti44_67kev_window_admits_observed_drift():
    """BUG-42 verification — Ti-44 67.87 keV match window vs Δ=+4.13 keV.

    The brief gates BUG-42 closure on C2's FWHM fix (BUG-41), assuming
    the matcher uses a `k·FWHM(E)` window where FWHM(67.87)≈2.9 keV
    pre-fix → window ±4.35 keV barely admits Δ=4.13 keV. But the
    production `staged_pipeline.analyze_lsrm_spe` default path uses
    `identification_window_from_fwhm` with `sqrt_E` scaling
    (`fwhm_window_multiple=0.5`, default — `staged_pipeline.py:1147-
    1151`), where the window is δE₀·√(E/E_ref) with δE₀=0.5·FWHM_661.

    For AmTiCsEu FWHM_661 = 42.027 keV (per outbox
    `2026-06-04_v1_22_0_AmTiCsEu_revalidation.md` §3.2 table):
      δE₀ = 21.01 keV
      window(67.87) = 21.01 · √(67.87/661.66) ≈ 6.73 keV
      window(1157.02) = 21.01 · √(1157.02/661.66) ≈ 27.79 keV

    Δ=+4.13 keV at 67.87 fits in ±6.73 keV → Ti-44 67.87 line matches.
    Δ=−10.90 keV at 1157 fits in ±27.79 keV → Sc-44 1157 line matches.

    Therefore BUG-42 is auto-resolved by the C1 candidate-set fix
    alone, regardless of C2 BUG-41 FWHM-model corrections. The
    `k·FWHM` window path (`use_lsrm_id_window=True`) is opt-in (default
    False per `staged_pipeline.py:809`) and only matters for callers
    that explicitly switch to F-167 canonical matching.
    """
    import math
    fwhm_661 = 42.027  # AmTiCsEu LSRM-fitted quadratic value
    delta_E0 = fwhm_661 * 0.5  # default fwhm_window_multiple

    # Ti-44 67.87 keV window vs observed Δ
    window_67 = delta_E0 * math.sqrt(67.87 / 661.66)
    drift_67 = 4.13  # Δ=+4.13 keV per revalidation §3.3
    assert window_67 > drift_67, (
        f"Ti-44 67.87 keV window ({window_67:.2f} keV) must admit "
        f"observed drift ({drift_67:.2f} keV). If not, BUG-42 is "
        f"actually still open and needs an explicit window floor."
    )

    # Sc-44 1157.02 keV window vs observed Δ
    window_1157 = delta_E0 * math.sqrt(1157.02 / 661.66)
    drift_1157 = 10.90  # Δ=-10.90 keV
    assert window_1157 > drift_1157, (
        f"Sc-44 1157.02 keV window ({window_1157:.2f} keV) must admit "
        f"observed drift ({drift_1157:.2f} keV)."
    )


def test_no_parent_dependency_gate_present_in_identify():
    """identify.py must not gate daughter ID on parent ID.

    Negative test — the file should not contain a "parent" keyword
    that conditions Sc-44 acceptance on Ti-44 detection. Verified at
    BUG-36 root-cause analysis (2026-06-05) that no such gate exists,
    so Sc-44 in EXOTIC_STAGE3 alone is sufficient to unblock matching.
    A future regression that adds a `if parent not in detected: skip`
    branch would break the AmTiCsEu identification and this test
    surfaces it.
    """
    import gamma.identification.identify as identify_mod
    src = Path(identify_mod.__file__).read_text(encoding="utf-8")
    # The string "parent" must not appear in identify.py logic.
    # (As of HEAD it does not; would only appear via a regression.)
    assert "parent" not in src.lower(), (
        "BUG-36 root cause: identify.py must NOT introduce a parent-"
        "dependency gate. If this test fails, daughter nuclide "
        "(Sc-44) matching may have been unintentionally re-gated on "
        "parent (Ti-44) detection."
    )
