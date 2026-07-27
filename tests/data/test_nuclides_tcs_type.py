"""Backward-compat + sanity tests for the v1.19.0 Phase 4 'tcs_type' field
extension on data/nuclides.json.

The field is OPTIONAL and additive: nuclides without 'tcs_type' must continue
to load normally. Tests verify (a) load-compat, (b) at least the 6
chain_progeny nuclides got the field.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
NUCLIDES_PATH = REPO_ROOT / "data" / "nuclides.json"

ALLOWED_TCS_TYPES = {"chain_progeny", "calibration_source", "single_line"}


@pytest.fixture(scope="module")
def nuclides() -> dict:
    return json.loads(NUCLIDES_PATH.read_text(encoding="utf-8"))


def test_tcs_type_field_optional(nuclides: dict) -> None:
    """Nuclides without 'tcs_type' continue to load: schema is additive.

    Pick a known nuclide that should NOT receive 'tcs_type' in v1.19.0 Phase 4
    (per PLAN §3 Phase 4 only 13 nuclides classified) and verify it still has
    its original required fields.
    """
    # Pick something that should NOT get tcs_type (e.g. Am-241, I-131, Be-7).
    sample_keys = [k for k in ("Be-7", "Am-241", "I-131", "Mn-54") if k in nuclides]
    assert sample_keys, "expected at least one of Be-7/Am-241/I-131/Mn-54 in nuclides.json"
    for key in sample_keys:
        entry = nuclides[key]
        # Must still have the canonical structure
        assert "T_half_s" in entry, f"{key}: T_half_s required field missing"
        assert "lines" in entry, f"{key}: lines required field missing"
        # tcs_type MAY be present but for these picks should NOT be
        # (the test is robust if a future commit adds them: at minimum require
        #  that whatever is there is in the allowed enum)
        if "tcs_type" in entry:
            assert entry["tcs_type"] in ALLOWED_TCS_TYPES, (
                f"{key}: tcs_type={entry['tcs_type']!r} not in {ALLOWED_TCS_TYPES}"
            )


def test_tcs_type_chain_progeny_set(nuclides: dict) -> None:
    """At least 6 nuclides carry tcs_type == 'chain_progeny' (PLAN §3 Phase 4:
    Bi-214, Pb-214, Bi-212, Pb-212, Tl-208, Ac-228).
    """
    chain_progeny = [
        k
        for k, v in nuclides.items()
        if isinstance(v, dict) and v.get("tcs_type") == "chain_progeny"
    ]
    assert len(chain_progeny) >= 6, (
        f"expected >= 6 chain_progeny nuclides, got {len(chain_progeny)}: {chain_progeny}"
    )
    # The 6 named in PLAN must all be present
    required = {"Bi-214", "Pb-214", "Bi-212", "Pb-212", "Tl-208", "Ac-228"}
    missing = required - set(chain_progeny)
    assert not missing, (
        f"required chain_progeny nuclides missing tcs_type: {sorted(missing)}"
    )


def test_tcs_type_values_in_enum(nuclides: dict) -> None:
    """Every nuclide carrying 'tcs_type' uses a value from the allowed enum."""
    bad = [
        (k, v["tcs_type"])
        for k, v in nuclides.items()
        if isinstance(v, dict)
        and "tcs_type" in v
        and v["tcs_type"] not in ALLOWED_TCS_TYPES
    ]
    assert not bad, (
        f"nuclides with tcs_type outside {ALLOWED_TCS_TYPES}: {bad}"
    )
