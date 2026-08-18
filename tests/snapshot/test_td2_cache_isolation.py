"""TD-2 / v1.18.29 — chain-equilibrium ratio must be cache-pollution-immune.

Background
----------
Before v1.18.29 the chain_equilibrium_guard ratio for the Th232 demo
fixture would jump from ≈2.30× (clean state, default JSON library) to
≈7.63× (after `load_lsrm_chain_libs()` was called in the same process)
because:
  • Lsrm's NaI-Etl+Esc.lib bundles all Th-232-chain daughter γ-lines
    under the parent `Th-232` entry.
  • Even with split_chains=True, lines that don't match any ENSDF
    daughter ownership in `chain_decomposer.TRUE_ENSDF_OWNERSHIP_*`
    remain assigned to the parent — so the post-merge `Th-232` record
    still carries 57 low-intensity unassigned lines.
  • When the staged_pipeline runs on a Th-232 spectrum it matches
    these phantom Th-232 lines (e.g. 233 keV @ I=0.11%) against real
    daughter photopeaks and produces a tiny aberrant A_Bq for the
    Th-232 head. That A_Bq then becomes min(A) in the chain, blowing
    the max/min ratio.

Root fix (v1.18.29)
-------------------
`chain_equilibrium_guard` now excludes chain-head nuclides
(`CHAIN_HEADS = {"Th-232", "U-238"}`) from the ratio computation
because Th-232 (T½ = 14 Gyr, α) and U-238 (T½ = 4.47 Gyr, α) emit no
direct γ-rays of practical intensity — any reported activity is a
library artifact. They still appear in the `members` block as
informational (`excluded_from_ratio=True`).

This test pins the cache-pollution-immune behavior as a regression
guard: regardless of whether `load_lsrm_chain_libs()` is called before
the equilibrium guard, the Th-232 chain ratio must remain physical.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Make `scripts/` importable — same convention as tests/conftest.py.
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from gamma.activity.compute import (  # noqa: E402
    CHAIN_HEADS, CHAIN_MEMBERS, chain_equilibrium_guard,
)
from gamma.data.nuclide_library import (  # noqa: E402
    load_lsrm_chain_libs, reset_cache,
)


def test_chain_heads_set_contains_known_heads():
    """CHAIN_HEADS must contain at least Th-232 and U-238 (the two
    natural decay-chain α-parents with no significant direct γ-emission)."""
    assert "Th-232" in CHAIN_HEADS
    assert "U-238" in CHAIN_HEADS


def test_chain_members_still_contain_heads_for_reporting():
    """CHAIN_MEMBERS must still list the chain heads — they appear in
    the `members` diagnostic block with excluded_from_ratio=True."""
    assert "Th-232" in CHAIN_MEMBERS["Th-232"]
    assert "U-238" in CHAIN_MEMBERS["U-238"]


def test_chain_equilibrium_guard_excludes_heads_from_ratio():
    """With a fake activity list including both Th-232 head and three
    daughters in equilibrium, the ratio must be computed from
    daughters ONLY — even if the head has an aberrant A_Bq."""
    from gamma.activity.compute import ActivityResult

    def _ar(name, A_Bq):
        return ActivityResult(
            nuclide=name,
            A_Bq=float(A_Bq),
            sigma_A_Bq=0.05 * float(A_Bq),
            lines_used=(),
            lines_skipped=(),
        )

    activities = [
        _ar("Th-232", 5.0),       # aberrant head (would be min if included)
        _ar("Ac-228", 2000.0),    # daughter
        _ar("Tl-208", 2000.0),    # daughter
        _ar("Pb-212", 2000.0),    # daughter
    ]
    eq = chain_equilibrium_guard(activities)
    th = eq.get("Th-232")
    assert th is not None
    # Ratio is 1.0 (all daughters equal), head is excluded.
    assert th["ratio"] == pytest.approx(1.0, abs=1e-9)
    # Head must still appear in members, but flagged excluded.
    names_in = {m["nuclide"]: m for m in th["members"]}
    assert "Th-232" in names_in
    assert names_in["Th-232"]["excluded_from_ratio"] is True
    # Daughters must NOT be flagged excluded.
    assert names_in["Ac-228"]["excluded_from_ratio"] is False
    assert names_in["Tl-208"]["excluded_from_ratio"] is False
    assert names_in["Pb-212"]["excluded_from_ratio"] is False


def test_th232_demo_chain_ratio_immune_to_lsrm_cache_pollution():
    """E2E: the Th-232 demo chain ratio must stay < 3× whether or not
    `load_lsrm_chain_libs()` was called previously in this process.

    The simulated pollution (calling load_lsrm_chain_libs at module
    import time) historically pushed ratio from ~2.30× to ~7.63× because
    the spurious phantom Th-232 activity became the chain min(A). The
    TD-2 fix isolates chain heads from the ratio.
    """
    fixture = (
        ROOT / "detectors" / "Gamma-1S" / "reference_spectra"
        / "archive" / "Th232_420-7-17_Маринелли_0cm.spe"
    )
    if not fixture.exists():
        pytest.skip(f"fixture missing: {fixture}")

    # Force the polluted state.
    reset_cache()
    load_lsrm_chain_libs()

    from gamma.identification.staged_pipeline import analyze_lsrm_spe

    r = analyze_lsrm_spe(str(fixture), complete_workflow=True)
    eq = chain_equilibrium_guard(r.activities or [])
    th_block = eq.get("Th-232")
    assert th_block is not None
    ratio = th_block["ratio"]
    assert ratio < 3.0, (
        f"Th-232 chain ratio={ratio:.2f}× must be < 3× even after "
        f"load_lsrm_chain_libs() pollutes the library cache "
        f"(TD-2 / v1.18.29 head-exclusion guard)"
    )
