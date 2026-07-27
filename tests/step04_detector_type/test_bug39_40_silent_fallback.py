"""BUG-39 / BUG-40 (Wave 6, v1.22.0; F2-A renormalisation 2026-06-21) —
detector profile loader + silent-fallback warning emission.

* **BUG-40**: ``references/detectors/<canonical>.json`` registry of
  detector profiles. Since F2-A only the primary ``Gamma-1S`` profile
  ships (the bogus stub ``Gamma-1S.json`` that previously co-existed and
  drove the Case-2 stub-fallback branch in :mod:`gamma.detectors.profile`
  was discarded).
* **BUG-39**: when ``detector_canonical`` resolves to a name whose
  profile is missing on disk (e.g. a future hypothetical
  ``AtomSpectra`` complex without a JSON profile), the pipeline emits a
  :class:`gamma.detectors.profile.DetectorFallback` record. ``report.json``
  ``warnings`` surfaces an operator-facing message.

Tests:
    A. ``load_detector_profile`` — JSON profile loading per canonical.
    B. ``detect_silent_fallback`` — fallback case classification
       (only Case 1 ``profile_not_on_disk`` and Case 3
       ``profile_loaded_no_fallback`` remain since F2-A).
    C. ``canonicalize`` — alias resolution maps the Cyrillic CONFIGNAME
       «Гамма-1С» and the legacy ASCII homoglyph ``Gamma-1C`` to the same
       single canonical ``Gamma-1S``.
    D. ``_build_warnings`` — warning emission only when reason is a
       real fallback (backward compat: empty / loaded-clean → no noise).
    E. Regression — Gamma-1S profile carries the LSRM FWHM polynomial
       provenance documented in baseline outbox §1.1 of 2026-06-04.
"""
from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.data.aliases import canonicalize                         # noqa: E402
from gamma.detectors.profile import (                               # noqa: E402
    clear_cache,
    detect_silent_fallback,
    load_detector_profile,
    should_emit_warning,
)
from gamma.reporting.json_report import _build_warnings             # noqa: E402


@pytest.fixture(autouse=True)
def _isolated_profile_cache():
    """Always start tests with a fresh LRU cache so test order is irrelevant."""
    clear_cache()
    yield
    clear_cache()


# ---------------------------------------------------------------------------
# Test A — JSON profile loading (BUG-40)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "canonical, expected_status, expect_loaded",
    [
        ("Gamma-1S", "primary", True),
        ("AtomSpectra", None, False),
        ("Gamma-XYZ", None, False),
        ("", None, False),
    ],
    ids=[
        "Gamma-1S-primary",
        "AtomSpectra-missing-profile",
        "Gamma-XYZ-missing-profile",
        "empty-canonical",
    ],
)
def test_load_detector_profile(canonical, expected_status, expect_loaded):
    profile = load_detector_profile(canonical)
    if not expect_loaded:
        assert profile is None, (
            f"Expected None for missing profile {canonical!r}, got {profile}"
        )
        return
    assert profile is not None, f"Profile for {canonical!r} should load"
    assert profile.canonical == canonical
    assert profile.validation_status == expected_status


# ---------------------------------------------------------------------------
# Test B — silent fallback classification (BUG-39)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "canonical, expected_reason, expected_actual, should_warn",
    [
        ("Gamma-1S", "profile_loaded_no_fallback", "Gamma-1S", False),
        ("Gamma-XYZ", "profile_not_on_disk", "Gamma-1S", True),
        ("AtomSpectra", "profile_not_on_disk", "Gamma-1S", True),
        ("", "profile_not_on_disk", "Gamma-1S", True),
    ],
    ids=[
        "Gamma-1S-clean",
        "Gamma-XYZ-missing",
        "AtomSpectra-missing",
        "empty-canonical",
    ],
)
def test_detect_silent_fallback(
    canonical, expected_reason, expected_actual, should_warn
):
    fb = detect_silent_fallback(canonical)
    assert fb.requested == canonical
    assert fb.actual == expected_actual
    assert fb.reason == expected_reason
    assert should_emit_warning(fb) is should_warn
    if should_warn:
        assert fb.human, "Fallback warning needs a human-readable message"
    else:
        assert fb.human == ""


# ---------------------------------------------------------------------------
# Test C — alias resolution maps every form of «Гамма-1С» to Gamma-1S
# ---------------------------------------------------------------------------
# F2-A canonical normalisation (2026-06-21): cyrillic CONFIGNAME «Гамма-1С»,
# legacy ASCII homoglyph "Gamma-1" + chr(0x43) (built via chr() so bulk
# rename passes never touch this fixture again) and the canonical ASCII
# "Gamma-1S" all resolve to the same canonical. The
# cyrillic_to_latin_collision predicate still fires for cyrillic-raw →
# ASCII-canonical pairs, but no warning is emitted because the profile
# loads cleanly (gate in json_report._build_warnings keys on
# detector_fallback.reason).
@pytest.mark.parametrize(
    "raw_text, expected_canonical",
    [
        ("Гамма-1С №SN-02", "Gamma-1S"),
        ("Гамма-1С", "Gamma-1S"),
        ("гамма-1с", "Gamma-1S"),
        ("Gamma-1S", "Gamma-1S"),
        ("Gamma-1" + chr(0x43), "Gamma-1S"),
        ("УДС-ГЦ-63х63-USB", "Gamma-1S"),
        ("БДЭГ-63×63", "Gamma-1S"),
    ],
    ids=[
        "configname-cyr-with-serial",
        "configname-cyr-bare",
        "configname-cyr-lower",
        "Gamma-1S-canonical",
        "legacy-ASCII-C-alias",
        "head-model-UDS",
        "head-model-BDEG",
    ],
)
def test_alias_resolution_collapses_all_forms_to_Gamma_1S(
    raw_text, expected_canonical
):
    actual = canonicalize("detector", raw_text)
    assert actual == expected_canonical, (
        f"canonicalize('detector', {raw_text!r}) = {actual!r}, "
        f"expected {expected_canonical!r}"
    )


# ---------------------------------------------------------------------------
# Test D — _build_warnings surfaces fallback to operator (BUG-39)
# ---------------------------------------------------------------------------
def _make_minimal_result(detector_fallback):
    """Construct a minimal SimpleNamespace stand-in for StagedAnalysisResult."""
    return SimpleNamespace(
        efficiency_curve=None,
        activities=None,
        mda_per_line=None,
        next_stage_recommended=None,
        next_stage_reason="",
        detector_fallback=detector_fallback,
    )


@pytest.mark.parametrize(
    "detector_fallback, expect_fallback_warning",
    [
        (
            {
                "reason": "profile_not_on_disk",
                "requested": "Gamma-XYZ",
                "actual": "Gamma-1S",
                "human": (
                    "Detector profile Gamma-XYZ.json not found — "
                    "fallback to Gamma-1S."
                ),
            },
            True,
        ),
        (
            {
                "reason": "profile_loaded_no_fallback",
                "requested": "Gamma-1S",
                "actual": "Gamma-1S",
                "human": "",
            },
            False,
        ),
        ({}, False),
        (None, False),
    ],
    ids=[
        "missing-profile-emits-warning",
        "clean-load-no-warning",
        "empty-dict-no-warning",
        "none-no-warning",
    ],
)
def test_build_warnings_detector_fallback(
    detector_fallback, expect_fallback_warning
):
    result = _make_minimal_result(detector_fallback)
    warnings = _build_warnings(result)

    fallback_msgs = [
        w for w in warnings
        if isinstance(w, str) and (
            "Detector profile" in w or "Detector fallback" in w
            or "profile fallback" in w.lower()
        )
    ]
    if expect_fallback_warning:
        assert fallback_msgs, (
            f"Expected a detector-fallback warning for {detector_fallback!r}, "
            f"got warnings = {warnings!r}"
        )
        human = detector_fallback.get("human", "") if detector_fallback else ""
        if human:
            assert any(human in w for w in warnings if isinstance(w, str))
    else:
        assert not fallback_msgs, (
            f"Did NOT expect a detector-fallback warning, but got: {fallback_msgs!r}"
        )


# ---------------------------------------------------------------------------
# Test E (regression guard) — Gamma-1S profile structural fields
# ---------------------------------------------------------------------------
# Per-spectrum FWHM coefficients live in each .spe (FWHM=2 block), not in
# the profile JSON; the profile only declares the LSRM √E convention.
def test_gamma1s_profile_fwhm_convention_from_baseline():
    prof = load_detector_profile("Gamma-1S")
    assert prof is not None
    fwhm = prof.raw.get("fwhm_polynomial", {}) or {}
    assert fwhm.get("convention") == "lsrm_sqrt_E"
    assert "coefficients" in fwhm
    # F2-A 2026-06-21: only "directory" efficiency source kind remains
    # (the legacy "TBD_pending_calibration_data" variant was removed
    # together with the bogus Gamma-1S stub profile that drove it).
    assert prof.efficiency_source_kind == "directory"
