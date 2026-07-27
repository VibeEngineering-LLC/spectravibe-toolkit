"""F-391 / v1.18.27 — display dedupe: no phantom-anchor in primary_feps,
no duplicate entries in compare HTML chips/tables.

Покрывает:
  1. ``_build_primary_feps`` фильтрует phantom anchors (library_anchor /
     library_anchor_phantom) из primary_feps JSON.
  2. compare HTML ``_peak_rows`` фильтрует phantom anchors из primary
     peaks table.
  3. compare HTML ``_peak_rows`` дедупит по (nuclide, round(E_keV, 0))
     — оставляет запись с максимальным S/σ.
"""
from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(
    Path(__file__).resolve().parent.parent.parent / "scripts"
))


# ──────────────────────────────────────────────────────────────────
# Fixtures
# ──────────────────────────────────────────────────────────────────

class _FakeLineMatch:
    """Stand-in for LineMatch с minimum fields consumed by _build_primary_feps."""

    def __init__(
        self,
        nuclide,
        library_E_keV,
        peak_channel,
        peak_E_keV=None,
        library_I_pct=10.0,
        peak_area=None,
        peak_area_uncertainty=None,
        peak_sigma=0.0,
        residual_keV=0.0,
        peak_area_source="cowell",
        is_characteristic=False,
    ):
        self.nuclide = nuclide
        self.library_E_keV = library_E_keV
        self.library_I_pct = library_I_pct
        self.peak_channel = peak_channel
        self.peak_E_keV = (
            peak_E_keV if peak_E_keV is not None else library_E_keV
        )
        self.peak_sigma = peak_sigma
        self.residual_keV = residual_keV
        self.is_characteristic = is_characteristic
        self.peak_area = peak_area
        self.peak_area_uncertainty = peak_area_uncertainty
        self.peak_area_source = peak_area_source


class _FakeNi:
    def __init__(self, nuclide, matched_lines):
        self.nuclide = nuclide
        self.matched_lines = matched_lines


# ──────────────────────────────────────────────────────────────────
# 1. JSON _build_primary_feps phantom filter
# ──────────────────────────────────────────────────────────────────

class TestPrimaryFepsPhantomFilter:
    def _make_result(self, matched_lines):
        result = SimpleNamespace()
        result.spec = SimpleNamespace(live_time=1000.0)
        result.final_detected = [_FakeNi("Ac-228", matched_lines)]
        return result

    def test_library_anchor_excluded_from_primary_feps(self):
        from gamma.reporting.json_report import _build_primary_feps

        real = _FakeLineMatch(
            "Ac-228", 911.0, 911,
            peak_area=10000.0, peak_area_uncertainty=500.0,
            peak_area_source="cowell",
        )
        anchor = _FakeLineMatch(
            "Ac-228", 463.0, 463,
            peak_area=None, peak_area_uncertainty=None,
            peak_area_source="library_anchor",
        )
        result = self._make_result([real, anchor])
        feps = _build_primary_feps(result)
        nucl_E = [(f["nuclide"], f["library_E_keV"]) for f in feps]
        assert ("Ac-228", 911.0) in nucl_E
        assert ("Ac-228", 463.0) not in nucl_E

    def test_library_anchor_phantom_excluded(self):
        from gamma.reporting.json_report import _build_primary_feps

        real = _FakeLineMatch(
            "Ac-228", 911.0, 911,
            peak_area=10000.0, peak_area_uncertainty=500.0,
            peak_area_source="deconvolved_coupled",
        )
        phantom = _FakeLineMatch(
            "Ac-228", 969.0, 969,
            peak_area=None, peak_area_uncertainty=None,
            peak_area_source="library_anchor_phantom",
        )
        result = self._make_result([real, phantom])
        feps = _build_primary_feps(result)
        nucl_E = [(f["nuclide"], f["library_E_keV"]) for f in feps]
        assert ("Ac-228", 911.0) in nucl_E
        assert ("Ac-228", 969.0) not in nucl_E

    def test_regular_sources_pass_through(self):
        from gamma.reporting.json_report import _build_primary_feps

        m1 = _FakeLineMatch(
            "Ac-228", 911.0, 911,
            peak_area=10000.0, peak_area_uncertainty=500.0,
            peak_area_source="cowell",
        )
        m2 = _FakeLineMatch(
            "Ac-228", 969.0, 969,
            peak_area=8000.0, peak_area_uncertainty=400.0,
            peak_area_source="deconvolved",
        )
        m3 = _FakeLineMatch(
            "Ac-228", 463.0, 463,
            peak_area=2000.0, peak_area_uncertainty=200.0,
            peak_area_source="lsrm_peaks_table",
        )
        result = self._make_result([m1, m2, m3])
        feps = _build_primary_feps(result)
        assert len(feps) == 3


# ──────────────────────────────────────────────────────────────────
# 2. compare HTML _peak_rows phantom filter + dedupe
# ──────────────────────────────────────────────────────────────────

class TestComparePeakRowsDedupe:
    def test_phantom_anchors_filtered(self):
        from gen_v2_compare_th232 import _peak_rows

        peaks = [
            {
                "nuclide": "Ac-228",
                "peak_channel": 911, "peak_E_keV": 911.0,
                "library_E_keV": 911.0, "library_I_pct": 27.0,
                "peak_area_counts": 10000.0,
                "peak_area_sigma": 500.0,
                "peak_area_source": "cowell",
            },
            {
                "nuclide": "Ac-228",
                "peak_channel": 463, "peak_E_keV": 463.0,
                "library_E_keV": 463.0, "library_I_pct": 4.4,
                "peak_area_counts": None,
                "peak_area_sigma": None,
                "peak_area_source": "library_anchor",
            },
            {
                "nuclide": "Tl-208",
                "peak_channel": 510, "peak_E_keV": 510.77,
                "library_E_keV": 510.77, "library_I_pct": 22.6,
                "peak_area_counts": None,
                "peak_area_sigma": None,
                "peak_area_source": "library_anchor_phantom",
            },
        ]
        rows = _peak_rows(peaks)
        nucl_E = {(r["nuclide"], r["library_E_keV"]) for r in rows}
        assert ("Ac-228", 911.0) in nucl_E
        assert ("Ac-228", 463.0) not in nucl_E
        assert ("Tl-208", 510.77) not in nucl_E

    def test_dedupe_keeps_highest_snr(self):
        """Две записи Ac-228 911 (lsrm + deconvolved_coupled) → одна
        с максимальным S/σ."""
        from gen_v2_compare_th232 import _peak_rows

        peaks = [
            {
                "nuclide": "Ac-228",
                "peak_channel": 911, "peak_E_keV": 911.0,
                "library_E_keV": 911.0, "library_I_pct": 27.0,
                # S/σ = 1000 / 500 = 2
                "peak_area_counts": 1000.0, "peak_area_sigma": 500.0,
                "peak_area_source": "lsrm_peaks_table",
            },
            {
                "nuclide": "Ac-228",
                "peak_channel": 911, "peak_E_keV": 911.4,
                "library_E_keV": 911.0, "library_I_pct": 27.0,
                # S/σ = 10000 / 500 = 20
                "peak_area_counts": 10000.0, "peak_area_sigma": 500.0,
                "peak_area_source": "deconvolved_coupled",
            },
        ]
        rows = _peak_rows(peaks)
        assert len(rows) == 1
        assert rows[0]["sigma"] >= 19.0

    def test_different_nuclide_not_deduped(self):
        """Ac-228 911 + Bi-214 911 (close E) — оба остаются."""
        from gen_v2_compare_th232 import _peak_rows

        peaks = [
            {
                "nuclide": "Ac-228",
                "peak_channel": 911, "peak_E_keV": 911.0,
                "library_E_keV": 911.0, "library_I_pct": 27.0,
                "peak_area_counts": 10000.0, "peak_area_sigma": 500.0,
                "peak_area_source": "cowell",
            },
            {
                "nuclide": "Bi-214",
                "peak_channel": 911, "peak_E_keV": 911.3,
                "library_E_keV": 911.5, "library_I_pct": 0.15,
                "peak_area_counts": 200.0, "peak_area_sigma": 100.0,
                "peak_area_source": "lsrm_peaks_table",
            },
        ]
        rows = _peak_rows(peaks)
        assert len(rows) == 2
        nucs = sorted(r["nuclide"] for r in rows)
        assert nucs == ["Ac-228", "Bi-214"]


# ──────────────────────────────────────────────────────────────────
# 3. compare HTML cluster table phantom rendering
# ──────────────────────────────────────────────────────────────────

class TestClusterTablePhantomRendering:
    def test_phantom_components_styled_separately(self):
        from gen_v2_compare_th232 import _multiplet_clusters_html

        clusters = [{
            "cluster_id": "M_TEST",
            "converged": True,
            "chi2_per_dof": 2.5,
            "closure_pct": 95,
            "E_lo_keV": 900.0, "E_hi_keV": 1000.0,
            "components": [
                {
                    "nuclide": "Ac-228",
                    "line_E_keV": 911.0, "library_I_pct": 27.0,
                    "deconvolved_area": 10000.0,
                    "deconvolved_area_sigma": 500.0,
                    "peak_area_source": "deconvolved_coupled",
                },
                {
                    "nuclide": "Ac-228",
                    "line_E_keV": 969.0, "library_I_pct": 16.0,
                    "deconvolved_area": 8000.0,
                    "deconvolved_area_sigma": 400.0,
                    "peak_area_source": "deconvolved_coupled",
                },
                {
                    "nuclide": "Ac-228",
                    "line_E_keV": 463.0, "library_I_pct": 4.4,
                    "deconvolved_area": 0.0,
                    "deconvolved_area_sigma": 0.0,
                    "peak_area_source": "library_anchor_phantom",
                },
            ],
        }]
        html = _multiplet_clusters_html(clusters, "test", "#abc")
        # Phantom row должен иметь class='phantom'
        assert "phantom" in html
        # Anchor-маркер в label
        assert "якорь" in html
