"""
Regression test for DEEP-01 (Project #5 wave 2 P1-1).

Bug: ``scripts/gamma/calibration/efficiency_autoload.py:118-120`` silently
returned ``None`` on BOTH "no candidate .efr file" AND "found a candidate
but the parser/fit raised". An operator running a Marinelli-1L analysis
on a corrupted .efr would unknowingly ship an efficiency-uncorrected
report — the silent ``return None`` was indistinguishable from "no
calibration available, degrade to qualitative".

Fix: the loader now emits ``logger.warning`` and returns a non-None
sentinel ``EFFICIENCY_FIT_FAILED`` on the fit-failed branch, while
"no file found" continues to return ``None`` silently.

Tests:
  * ``test_broken_efr_emits_warning_and_returns_sentinel`` — point the
    loader at a directory containing a deliberately broken
    Marinelli-tagged .efr file. Assert a WARNING was emitted naming the
    file basename and the return value IS the sentinel (NOT None).
    This test would fail (return None, no warning) against the pre-fix
    code at lines 118-120.
  * ``test_no_efr_returns_none_and_silent`` — point the loader at an
    empty directory. Assert return IS None AND no WARNING was emitted.
    Proves the two failure modes are now distinguishable.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pytest

from gamma.calibration import efficiency_autoload
from gamma.calibration.efficiency_autoload import (
    EFFICIENCY_FIT_FAILED, load_efficiency_for_geometry,
    find_efr_file,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#: A deliberately broken .efr payload. The Lsrm reader will decode the
#: CP-1251 bytes without error (latin-1 fallback also OK), parse it as a
#: single header block with ZERO energy points, and
#: ``fit_efficiency_from_efr_file`` will raise
#: ``ValueError("No efficiency points in: ...")`` (efficiency.py:278).
#: This exercises the fit-fail except-branch without us having to fake
#: the fit itself.
_BROKEN_EFR_CONTENT = (
    "[Gamma-1S;Маринелли-1л;BROKEN]\r\n"
    "DetectorType=NaI\r\n"
    "# no energy points at all — fit will raise ValueError\r\n"
)


def _build_broken_efr_tree(tmp_path: Path) -> Path:
    """
    Layout mirroring ``detectors/Gamma-1S/efficiency/<unit>/<geometry>.efr``.

    Returns the synthetic EFR root suitable for monkeypatching
    ``efficiency_autoload._EFR_ROOT``.
    """
    root = tmp_path / "efr_root"
    unit_dir = root / "Gamma-1S_NaI_63x63_USB_TESTUNIT"
    unit_dir.mkdir(parents=True)
    broken = unit_dir / "ТЕСТ-Маринелли-1л.efr"
    broken.write_text(_BROKEN_EFR_CONTENT, encoding="cp1251")
    return root


def _build_empty_efr_tree(tmp_path: Path) -> Path:
    """Same layout, no .efr files anywhere."""
    root = tmp_path / "efr_root_empty"
    (root / "Gamma-1S_NaI_63x63_USB_TESTUNIT").mkdir(parents=True)
    return root


@pytest.fixture(autouse=True)
def _clear_lru_caches():
    """
    Both ``find_efr_file`` and ``load_efficiency_for_geometry`` are
    ``@lru_cache``-decorated at module scope. Each test must start from a
    cold cache or the second test will receive a cached return from the
    first (different ``_EFR_ROOT``).
    """
    find_efr_file.cache_clear()
    load_efficiency_for_geometry.cache_clear()
    yield
    find_efr_file.cache_clear()
    load_efficiency_for_geometry.cache_clear()


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_broken_efr_emits_warning_and_returns_sentinel(
    tmp_path, monkeypatch, caplog,
):
    """
    DEEP-01 red-without-fix anchor.

    A candidate .efr that parses to zero points causes
    ``fit_efficiency_from_efr_file`` to raise. Pre-fix code at
    ``efficiency_autoload.py:118-120`` swallowed the exception and
    returned ``None`` — indistinguishable from "no .efr file found".

    Post-fix: a WARNING is emitted naming the broken file's basename and
    the return value is the module-level singleton
    ``EFFICIENCY_FIT_FAILED`` (NOT None).
    """
    root = _build_broken_efr_tree(tmp_path)
    monkeypatch.setattr(efficiency_autoload, "_EFR_ROOT", root)

    with caplog.at_level(logging.WARNING, logger=efficiency_autoload.__name__):
        result = load_efficiency_for_geometry(
            "Маринелли 1л", "Gamma-1S",
        )

    # Return-value contract: distinct from None.
    assert result is not None, (
        "Pre-fix bug: load_efficiency_for_geometry returned None on a "
        "broken .efr — operator could not distinguish from 'no file'."
    )
    assert result is EFFICIENCY_FIT_FAILED, (
        f"Expected EFFICIENCY_FIT_FAILED sentinel, got {result!r}."
    )

    # Logger contract: WARNING-level record mentioning the broken file.
    warnings = [
        r for r in caplog.records
        if r.levelno == logging.WARNING
        and r.name == efficiency_autoload.__name__
    ]
    assert warnings, (
        "Pre-fix bug: no WARNING was emitted — the broken .efr was "
        "silently dropped from the operator's view."
    )
    msgs = " | ".join(r.getMessage() for r in warnings)
    assert "ТЕСТ-Маринелли-1л.efr" in msgs, (
        f"Warning must name the offending file basename (F-115). "
        f"Got: {msgs!r}"
    )

    # F-115: no absolute operator path leaks into the log message.
    assert str(tmp_path) not in msgs, (
        "F-115 violation: warning message contains the operator's "
        f"absolute path {tmp_path!s}. Got: {msgs!r}"
    )


def test_no_efr_returns_none_and_silent(
    tmp_path, monkeypatch, caplog,
):
    """
    Regression guard for the other branch.

    With NO candidate .efr files in the tree, ``find_efr_file`` returns
    None and ``load_efficiency_for_geometry`` must keep returning None
    silently (no WARNING). This proves the two states are now
    distinguishable: caller sees None vs. EFFICIENCY_FIT_FAILED.
    """
    root = _build_empty_efr_tree(tmp_path)
    monkeypatch.setattr(efficiency_autoload, "_EFR_ROOT", root)

    with caplog.at_level(logging.WARNING, logger=efficiency_autoload.__name__):
        result = load_efficiency_for_geometry(
            "Маринелли 1л", "Gamma-1S",
        )

    assert result is None, (
        f"With no .efr present, loader must return None (not the "
        f"fit-failed sentinel). Got: {result!r}"
    )
    assert result is not EFFICIENCY_FIT_FAILED, (
        "'no file' must NOT collapse into the fit-failed sentinel."
    )

    # No WARNING — 'no file' is a normal degradation path.
    warnings = [
        r for r in caplog.records
        if r.levelno == logging.WARNING
        and r.name == efficiency_autoload.__name__
    ]
    assert not warnings, (
        f"'No .efr found' must be silent; got warnings: "
        f"{[r.getMessage() for r in warnings]}"
    )
