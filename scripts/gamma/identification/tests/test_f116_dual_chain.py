"""Unit tests for F-116-v2 dual-chain union fix (#124, v1.19.1).

Background
----------
Pre-v1.19.1 F-116 used `if/elif`:
    if chain_dominance_out.th232: chain_set = _TH_CHAIN_NUCLIDES
    else:                         chain_set = _U_CHAIN_NUCLIDES

When BOTH Th-232 and U-238 reverse-flagged DOMINANT (typical natural
background), only the Th set was used → U-cycle members (Bi-214, Pb-214,
Pb-210, Ra-226) were classified "вне цепочки Th-232" and dropped as
spurious. v1.19.1 unions all DOMINANT-flagged chain sets, restoring U-cycle
nuclides while preserving single-chain behaviour.

Cases:
  1. Th-only DOMINANT  → Bi-214 (U-cycle, no strong confirm) SUPPRESSED.
  2. U-only  DOMINANT  → Tl-208 (Th-cycle, no strong confirm) SUPPRESSED.
  3. BOTH    DOMINANT  → BOTH Bi-214 AND Tl-208 KEPT (the fix).
  4. Neither DOMINANT  → no F-116 effect; both kept untouched.

Asserts target:
  - membership of `final_detected` (by nuclide name)
  - membership of `out_of_chain_suppressed` log
"""

from __future__ import annotations

import os
import sys
import unittest
from types import SimpleNamespace
from unittest.mock import MagicMock

# Make 'gamma.*' importable when test is run directly via pytest from repo root.
_HERE = os.path.dirname(os.path.abspath(__file__))
_SCRIPTS = os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)

from gamma.identification.staged_pipeline import (  # noqa: E402
    _apply_f116_out_of_chain_suppression,
    _TH_CHAIN_NUCLIDES,
    _U_CHAIN_NUCLIDES,
    _BACKGROUND_NATURAL_NUCLIDES,
)


def _fake_ni(
    nuclide: str,
    *,
    n_matched: int = 1,
    sigma: float = 1.0,
    ci: float = 1.0,
) -> SimpleNamespace:
    """Minimal NuclideIdentification stand-in.

    F-116 inspects only:
      ni.nuclide, ni.matched_lines (len + .peak_area / .peak_area_uncertainty
      via getattr), ni.confidence.CI.
    By default builds a "weak" candidate (1 line, σ≈1, CI≈1) → cannot pass
    the σ≥10 AND ≥2 lines AND CI≥5 strong-confirm gate.
    """
    matched = []
    if n_matched > 0:
        # Encode σ via peak_area / peak_area_uncertainty.
        area_unc = 1.0
        for _ in range(n_matched):
            matched.append(SimpleNamespace(
                peak_area=sigma * area_unc,
                peak_area_uncertainty=area_unc,
            ))
    confidence = SimpleNamespace(CI=ci) if ci is not None else None
    return SimpleNamespace(
        nuclide=nuclide,
        matched_lines=matched,
        confidence=confidence,
    )


def _names(final_detected) -> set:
    return {ni.nuclide for ni in final_detected}


def _suppressed_nucs(messages) -> set:
    """Extract the suppressed-nuclide name from each 'F-116: подавлен X ...'."""
    out = set()
    for m in messages:
        # Format: 'F-116: подавлен <nuc> (вне цепочки ...)'
        parts = m.split()
        if len(parts) >= 3 and parts[0] == "F-116:" and parts[1] == "подавлен":
            out.add(parts[2])
    return out


class TestF116DualChainUnion(unittest.TestCase):

    def setUp(self):
        # Sanity: confirm constants haven't drifted.
        self.assertIn("Bi-214", _U_CHAIN_NUCLIDES)
        self.assertIn("Pb-214", _U_CHAIN_NUCLIDES)
        self.assertIn("Tl-208", _TH_CHAIN_NUCLIDES)
        # K-40 is whitelisted as background-natural (F-130), not in chain sets.
        self.assertIn("K-40", _BACKGROUND_NATURAL_NUCLIDES)
        self.assertNotIn("K-40", _TH_CHAIN_NUCLIDES)
        self.assertNotIn("K-40", _U_CHAIN_NUCLIDES)

    # -- Case 1 -----------------------------------------------------------
    def test_th_only_dominant_suppresses_bi214(self):
        """Th-only DOMINANT: Bi-214 (U-cycle, weak) should be SUPPRESSED."""
        chain_dom = MagicMock(th232=True, u238=False)
        final = [
            _fake_ni("Tl-208"),   # Th-cycle → kept
            _fake_ni("Bi-214"),   # U-cycle, weak → suppressed
        ]
        kept, suppressed = _apply_f116_out_of_chain_suppression(
            final_detected=final,
            chain_dominance_out=chain_dom,
            filename_isotope_hints=None,
        )
        self.assertIn("Tl-208", _names(kept))
        self.assertNotIn("Bi-214", _names(kept))
        self.assertIn("Bi-214", _suppressed_nucs(suppressed))
        # Label should reflect the active chain.
        self.assertTrue(any("Th-232" in m for m in suppressed))

    # -- Case 2 -----------------------------------------------------------
    def test_u_only_dominant_suppresses_tl208(self):
        """U-only DOMINANT: Tl-208 (Th-cycle, weak) should be SUPPRESSED."""
        chain_dom = MagicMock(th232=False, u238=True)
        final = [
            _fake_ni("Bi-214"),   # U-cycle → kept
            _fake_ni("Tl-208"),   # Th-cycle, weak → suppressed
        ]
        kept, suppressed = _apply_f116_out_of_chain_suppression(
            final_detected=final,
            chain_dominance_out=chain_dom,
            filename_isotope_hints=None,
        )
        self.assertIn("Bi-214", _names(kept))
        self.assertNotIn("Tl-208", _names(kept))
        self.assertIn("Tl-208", _suppressed_nucs(suppressed))
        self.assertTrue(any("U-238" in m for m in suppressed))

    # -- Case 3 (the fix) -------------------------------------------------
    def test_both_dominant_keeps_both_chains(self):
        """BOTH DOMINANT: union → Bi-214 AND Tl-208 both KEPT (#124 fix)."""
        chain_dom = MagicMock(th232=True, u238=True)
        final = [
            _fake_ni("Tl-208"),    # Th-cycle → kept (always)
            _fake_ni("Bi-214"),    # U-cycle → kept ONLY with union fix
            _fake_ni("Pb-214"),    # U-cycle → kept ONLY with union fix
            _fake_ni("Ra-226"),    # U-cycle → kept ONLY with union fix
            _fake_ni("Am-241"),    # neither chain, weak → suppressed
        ]
        kept, suppressed = _apply_f116_out_of_chain_suppression(
            final_detected=final,
            chain_dominance_out=chain_dom,
            filename_isotope_hints=None,
        )
        kept_names = _names(kept)
        self.assertIn("Tl-208", kept_names, "Th-cycle member should remain")
        self.assertIn("Bi-214", kept_names, "U-cycle member NOT suppressed under union")
        self.assertIn("Pb-214", kept_names, "U-cycle member NOT suppressed under union")
        self.assertIn("Ra-226", kept_names, "U-cycle member NOT suppressed under union")
        # Am-241 is outside both chains and has no special whitelist.
        self.assertNotIn("Am-241", kept_names)
        self.assertIn("Am-241", _suppressed_nucs(suppressed))
        # Label should reflect both chains being active.
        self.assertTrue(any("Th-232+U-238" in m for m in suppressed))
        # And NO U-cycle / Th-cycle member should appear in suppressed log.
        suppressed_set = _suppressed_nucs(suppressed)
        self.assertFalse(suppressed_set & _U_CHAIN_NUCLIDES)
        self.assertFalse(suppressed_set & _TH_CHAIN_NUCLIDES)

    # -- Case 4 -----------------------------------------------------------
    def test_neither_dominant_no_effect(self):
        """Neither DOMINANT: F-116 inactive; all candidates pass through."""
        chain_dom = MagicMock(th232=False, u238=False)
        final = [
            _fake_ni("Bi-214"),
            _fake_ni("Tl-208"),
            _fake_ni("Cs-137"),
        ]
        kept, suppressed = _apply_f116_out_of_chain_suppression(
            final_detected=final,
            chain_dominance_out=chain_dom,
            filename_isotope_hints=None,
        )
        self.assertEqual(_names(kept), {"Bi-214", "Tl-208", "Cs-137"})
        self.assertEqual(suppressed, [])

    # -- Case 5: F-130 background-natural whitelist -----------------------
    def test_k40_kept_by_background_natural_whitelist(self):
        """K-40 survives F-116 even when BOTH chains are DOMINANT (#130 fix).

        K-40 is primordial and ubiquitous in natural background but belongs
        to neither Th-232 nor U-238 chains. _BACKGROUND_NATURAL_NUCLIDES
        whitelist prevents its spurious suppression in background-only spectra.
        """
        chain_dom = MagicMock(th232=True, u238=True)
        final = [
            _fake_ni("Tl-208"),    # Th-cycle → kept
            _fake_ni("Bi-214"),    # U-cycle → kept
            _fake_ni("K-40"),      # background-natural whitelist → kept (F-130)
            _fake_ni("Am-241"),    # not whitelisted, weak → suppressed
        ]
        kept, suppressed = _apply_f116_out_of_chain_suppression(
            final_detected=final,
            chain_dominance_out=chain_dom,
            filename_isotope_hints=None,
        )
        kept_names = _names(kept)
        self.assertIn("K-40", kept_names, "K-40 must survive due to background-natural whitelist")
        self.assertNotIn("K-40", _suppressed_nucs(suppressed))
        # Confirm Am-241 (not whitelisted) is still suppressed.
        self.assertNotIn("Am-241", kept_names)
        self.assertIn("Am-241", _suppressed_nucs(suppressed))

    # -- Bonus: strong-confirm gate still works under union ---------------
    def test_strong_confirm_bypasses_suppression_under_union(self):
        """Strong outsider (σ≥10, ≥2 lines, CI≥5) is kept even when both chains DOMINANT."""
        chain_dom = MagicMock(th232=True, u238=True)
        strong_outsider = _fake_ni("Cs-137", n_matched=2, sigma=15.0, ci=7.0)
        weak_outsider = _fake_ni("Co-57", n_matched=1, sigma=2.0, ci=1.0)
        final = [strong_outsider, weak_outsider]
        kept, suppressed = _apply_f116_out_of_chain_suppression(
            final_detected=final,
            chain_dominance_out=chain_dom,
            filename_isotope_hints=None,
        )
        kept_names = _names(kept)
        self.assertIn("Cs-137", kept_names, "Strong outsider must survive")
        self.assertNotIn("Co-57", kept_names)
        self.assertIn("Co-57", _suppressed_nucs(suppressed))


if __name__ == "__main__":
    unittest.main()
