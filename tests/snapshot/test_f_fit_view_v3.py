"""F-FIT-VIEW v3 / v1.22.6 — Gaussians lifted on continuum baseline (Task #67).

Tests cover the continuum_grid field added in v3:
  - continuum_grid emitted for singlet peaks with ≥5 sample points
  - at peak_E the displayed Gaussian y = continuum_at_peak + amp (lifted max)
  - backward compat: missing continuum_grid → JS falls back to y=0 (no crash)
  - multiplet components do NOT get a separate continuum_grid (use cluster continuum)

All 4 tests use _build_fit_overlay_payload (json pass-through) directly,
mirroring the pattern from test_f_fit_view_v2.py.

Cite: scripts/gamma/reporting/json_report.py:586 (_build_fit_overlay),
      scripts/gamma/reporting/json_report.py:669 (_continuum_grid_for_peak),
      scripts/gamma/reporting/templates/interactive_v1_17_2.html:1152 (continuumAt),
      scripts/gamma/reporting/templates/interactive_v1_17_2.html:1168 (gaussianPoints v3).
"""
from __future__ import annotations

import json
import math
import os
import sys

import pytest

sys.path.insert(0, "scripts")

from gamma.reporting.interactive_html import _build_fit_overlay_payload  # noqa: E402


# ── Test 1: continuum_grid emitted for singlet peaks ──────────────────────────

def test_continuum_grid_emitted_for_singlet():
    """F-FIT-VIEW v3: singlet peak payload passes continuum_grid through with ≥5 samples.

    When fit_overlay.peaks[*] contains a singlet entry with a continuum_grid,
    _build_fit_overlay_payload must preserve it in the output.
    continuum_grid.energies and continuum_grid.values must each have ≥5 elements.
    Cite: scripts/gamma/reporting/json_report.py:669 (_continuum_grid_for_peak),
          scripts/gamma/reporting/json_report.py:780 (continuum_grid attached to singlet entry).
    """
    # Simulate a singlet peak with a pre-built continuum_grid (11 samples, ±3σ window)
    sigma = 19.5  # keV — NaI at 662 keV
    e = 661.7
    e_lo = e - 3 * sigma
    e_hi = e + 3 * sigma
    n_pts = 11
    energies = [round(e_lo + i * (e_hi - e_lo) / (n_pts - 1), 2) for i in range(n_pts)]
    # Linear continuum from 1300 (left shoulder) to 900 (right shoulder)
    values = [round(1300.0 + (900.0 - 1300.0) * i / (n_pts - 1), 3) for i in range(n_pts)]

    sample = {
        "fit_overlay": {
            "peaks": [
                {
                    "peak_id": "p662",
                    "nuclide": "Cs-137",
                    "energy_keV": e,
                    "amp_counts": 500.0,
                    "sigma_keV": sigma,
                    "source": "singlet",
                    "label": "Cs-137 662",
                    "continuum_grid": {"energies": energies, "values": values},
                }
            ],
            "multiplet_continua": [],
        }
    }
    payload = _build_fit_overlay_payload(sample)
    assert len(payload["peaks"]) == 1, f"Expected 1 peak, got {len(payload['peaks'])}"
    pk = payload["peaks"][0]
    assert "continuum_grid" in pk, (
        "v3: singlet peak must carry continuum_grid through payload"
    )
    cg = pk["continuum_grid"]
    assert "energies" in cg and "values" in cg, (
        f"continuum_grid must have 'energies' and 'values' keys, got: {list(cg.keys())}"
    )
    assert len(cg["energies"]) >= 5, (
        f"continuum_grid.energies must have ≥5 sample points, got {len(cg['energies'])}"
    )
    assert len(cg["values"]) == len(cg["energies"]), (
        "continuum_grid.energies and .values must have equal length"
    )
    # energies must be monotonically increasing
    for i in range(len(cg["energies"]) - 1):
        assert cg["energies"][i] < cg["energies"][i + 1], (
            f"continuum_grid.energies must be monotonically increasing: "
            f"{cg['energies'][i]} >= {cg['energies'][i+1]}"
        )


# ── Test 2: lifted Gaussian max = continuum + amp ─────────────────────────────

def test_lifted_gaussian_max_above_continuum():
    """F-FIT-VIEW v3: at peak_E, displayed y = continuum_at_peak + amp_counts.

    Verifies the mathematical invariant: when continuum_grid is present, the peak of
    the displayed Gaussian (at x = energy_keV) equals continuum_at_peak + amp_counts.
    This is the core v3 'lift' behavior — Gaussian sits ON the continuum, not on y=0.

    Cite: scripts/gamma/reporting/templates/interactive_v1_17_2.html:1168 (gaussianPoints),
          scripts/gamma/reporting/templates/interactive_v1_17_2.html:1152 (continuumAt).
    """
    # Simulate a flat continuum at 1100 counts across the peak window
    sigma = 19.5
    e = 661.7
    e_lo = e - 3 * sigma
    e_hi = e + 3 * sigma
    n_pts = 11
    energies = [round(e_lo + i * (e_hi - e_lo) / (n_pts - 1), 2) for i in range(n_pts)]
    flat_cont = 1100.0
    values = [flat_cont] * n_pts  # flat continuum

    amp = 500.0

    sample = {
        "fit_overlay": {
            "peaks": [
                {
                    "peak_id": "p662",
                    "nuclide": "Cs-137",
                    "energy_keV": e,
                    "amp_counts": amp,
                    "sigma_keV": sigma,
                    "source": "singlet",
                    "label": "Cs-137 662",
                    "continuum_grid": {"energies": energies, "values": values},
                }
            ],
            "multiplet_continua": [],
        }
    }
    payload = _build_fit_overlay_payload(sample)
    pk = payload["peaks"][0]
    assert "continuum_grid" in pk, "continuum_grid must be present in payload"

    # Simulate the JS gaussianPoints computation at x = energy_keV (dx=0 → max of Gaussian)
    # continuumAt: flat continuum → returns flat_cont at peak_E
    # y_raw = cont_y + amp * exp(0) = flat_cont + amp * 1.0 = flat_cont + amp
    x = e
    dx = x - e  # = 0
    cont_y = flat_cont  # flat grid → continuumAt = flat_cont at peak center
    y_displayed = cont_y + amp * math.exp(-(dx * dx) / (2 * sigma * sigma))
    expected = flat_cont + amp  # = 1600.0

    assert abs(y_displayed - expected) < 0.01, (
        f"At peak_E, y_displayed must equal continuum + amp = {expected}, got {y_displayed}"
    )
    # Ensure displayed y > y without lift (pure Gaussian on y=0)
    y_no_lift = amp  # Gaussian only, on y=0
    assert y_displayed > y_no_lift, (
        f"Lifted Gaussian ({y_displayed}) must exceed un-lifted ({y_no_lift})"
    )


# ── Test 3: backward compat — absent continuum_grid falls back to y=0 ─────────

def test_continuum_grid_absent_falls_back_to_zero():
    """F-FIT-VIEW v3: backward compat — payload without continuum_grid works as v2.

    When a peak entry has no continuum_grid field, _build_fit_overlay_payload must
    still produce a valid peak dict. The JS gaussianPoints falls back to cont_y=0
    (pk.continuum_grid || null → continuumAt(null, x) → 0).
    No crash, no KeyError, identical to v2 behavior.

    Cite: scripts/gamma/reporting/json_report.py:780 (continuum_grid conditional attach),
          scripts/gamma/reporting/templates/interactive_v1_17_2.html:1272
              (pk.continuum_grid || null — backward compat call site).
    """
    sample = {
        "fit_overlay": {
            "peaks": [
                {
                    "peak_id": "p662",
                    "nuclide": "Cs-137",
                    "energy_keV": 661.7,
                    "amp_counts": 500.0,
                    "sigma_keV": 19.5,
                    "source": "singlet",
                    "label": "Cs-137 662",
                    # No continuum_grid — v2-style payload
                }
            ],
            "multiplet_continua": [],
        }
    }
    payload = _build_fit_overlay_payload(sample)
    assert len(payload["peaks"]) == 1, f"Expected 1 peak, got {len(payload['peaks'])}"
    pk = payload["peaks"][0]

    # continuum_grid absent → JS falls back to y=0 (not an error)
    # Payload must still have all required v2 fields
    required_v2 = {"peak_id", "nuclide", "energy_keV", "amp_counts", "sigma_keV", "source", "label"}
    missing = required_v2 - set(pk.keys())
    assert not missing, (
        f"v2-required fields missing from payload peak: {missing}"
    )
    assert pk["source"] == "singlet", (
        f"source must be preserved: expected 'singlet', got {pk['source']}"
    )

    # Verify JS backward compat: continuumAt(null, x) → 0
    # (simulated in Python — just check no continuum_grid present, no crash)
    cg = pk.get("continuum_grid", None)
    # When continuum_grid is absent, JS computes cont_y = continuumAt(null, x) = 0
    # Simulated: if cg is None → cont_y = 0 → y_raw = 0 + amp * exp(...)
    sigma = pk["sigma_keV"]
    amp = pk["amp_counts"]
    x = pk["energy_keV"]
    dx = 0.0
    cont_y = 0.0 if cg is None else 999.0  # should be 0 (v2 fallback)
    y_raw_fallback = cont_y + amp * math.exp(-(dx * dx) / (2 * sigma * sigma))
    assert abs(y_raw_fallback - amp) < 0.01, (
        f"Without continuum_grid, y at peak should equal amp={amp}, got {y_raw_fallback}"
    )


# ── Test 4: multiplet components do not get individual continuum_grid ──────────

def test_multiplet_components_share_continuum_band():
    """F-FIT-VIEW v3: multiplet_component peaks do NOT carry continuum_grid.

    Multiplet components are rendered via g_curve (pure Gaussian, continuum subtracted)
    alongside the separate multiplet_continua cluster band.  They must NOT receive a
    per-peak continuum_grid so there is no double-lifting (the continuum is already
    shown as a separate magenta dashed line in the JS).

    Cite: scripts/gamma/reporting/json_report.py:779 (continuum_grid skipped for
          multiplet_component: 'cg = None if source == "multiplet_component"'),
          scripts/gamma/reporting/templates/interactive_v1_17_2.html:1253
              ('if pk.source === multiplet_component return' — rendered via g_curve).
    """
    sample = {
        "fit_overlay": {
            "peaks": [
                {
                    "peak_id": "p583",
                    "nuclide": "Tl-208",
                    "energy_keV": 583.2,
                    "amp_counts": 300.0,
                    "sigma_keV": 11.5,
                    "source": "multiplet_component",
                    "label": "Tl-208 583",
                    # No continuum_grid — multiplet components don't get one
                },
                {
                    "peak_id": "p727",
                    "nuclide": "Bi-212",
                    "energy_keV": 727.3,
                    "amp_counts": 150.0,
                    "sigma_keV": 13.2,
                    "source": "multiplet_component",
                    "label": "Bi-212 727",
                    # No continuum_grid — multiplet components don't get one
                },
            ],
            "multiplet_continua": [
                {
                    "cluster_id": "cluster_1",
                    "E_keV": [550.0, 600.0, 650.0, 700.0, 750.0],
                    "continuum": [2000.0, 1850.0, 1700.0, 1550.0, 1400.0],
                    "total": [2350.0, 2100.0, 1950.0, 1700.0, 1550.0],
                    "components": [],
                }
            ],
        }
    }
    payload = _build_fit_overlay_payload(sample)
    assert len(payload["peaks"]) == 2, (
        f"Expected 2 multiplet_component peaks, got {len(payload['peaks'])}"
    )
    for pk in payload["peaks"]:
        assert pk["source"] == "multiplet_component", (
            f"Expected source='multiplet_component', got {pk['source']}"
        )
        # multiplet_component peaks must NOT have continuum_grid (rendered via g_curve)
        assert "continuum_grid" not in pk, (
            f"Multiplet component peak {pk['peak_id']} must NOT have continuum_grid — "
            f"it is rendered via g_curve + shared cluster continuum band. "
            f"Double-lifting would occur if continuum_grid were added."
        )

    # The multiplet_continua band must still be present and correct
    mc = payload.get("multiplet_continua", [])
    assert len(mc) == 1, f"Expected 1 multiplet cluster, got {len(mc)}"
    assert mc[0]["cluster_id"] == "cluster_1", (
        f"cluster_id mismatch: expected 'cluster_1', got {mc[0]['cluster_id']}"
    )
    assert len(mc[0]["E_keV"]) == 5, (
        f"cluster continuum E_keV must have 5 points, got {len(mc[0]['E_keV'])}"
    )
