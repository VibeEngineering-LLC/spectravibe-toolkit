"""F-121 (v1.17.6) — Ra-226/U-238 forced multiplets regression.

Зеркальное правило для F-118 на цепочке U-238. Проверяет, что при
chain_dominance.u238 == True пайплайн всегда эмиттирует три кластера
U1/U2/U3 со ВСЕМИ компонентами по библиотечным интенсивностям.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

SCRIPTS = os.path.join(os.path.dirname(__file__), "scripts")
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)

from gamma.peaks.deconvolve import (
    RA226_FORCED_CLUSTERS, run_chain_forced_multiplets,
)


class _FakeChainDom:
    """Minimal stub of ChainDominance for the forced-multiplet branch."""
    def __init__(self, th232: bool = False, u238: bool = False):
        self.th232 = th232
        self.u238 = u238


def test_ra226_forced_clusters_three_clusters():
    """F-121: ровно 3 кластера U1/U2/U3 c корректными ROI."""
    assert len(RA226_FORCED_CLUSTERS) == 3
    ids = [c["id"] for c in RA226_FORCED_CLUSTERS]
    assert ids == ["U1", "U2", "U3"]
    # ROI: U1 ~ 270-380, U2 ~ 580-830, U3 ~ 1080-1830
    assert RA226_FORCED_CLUSTERS[0]["E_lo_keV"] < 300.0
    assert RA226_FORCED_CLUSTERS[0]["E_hi_keV"] < 400.0
    assert RA226_FORCED_CLUSTERS[1]["E_lo_keV"] > 500.0
    assert RA226_FORCED_CLUSTERS[1]["E_hi_keV"] < 900.0
    assert RA226_FORCED_CLUSTERS[2]["E_lo_keV"] > 1000.0
    assert RA226_FORCED_CLUSTERS[2]["E_hi_keV"] > 1700.0


def test_ra226_forced_clusters_required_lines():
    """F-121: ключевые библиотечные линии Bi-214/Pb-214 присутствуют."""
    all_E = []
    for cluster in RA226_FORCED_CLUSTERS:
        for n, E, I, g in cluster["components"]:
            all_E.append((n, round(float(E), 1)))
    # Pb-214 295, 352
    assert ("Pb-214", 295.2) in all_E or ("Pb-214", 295.2) in all_E
    assert any(n == "Pb-214" and abs(E - 295.2) < 0.5 for n, E in all_E)
    assert any(n == "Pb-214" and abs(E - 351.9) < 0.5 for n, E in all_E)
    # Bi-214 609, 1120, 1238, 1378, 1408, 1764
    for E_expect in (609.3, 1120.3, 1238.1, 1377.7, 1408.0, 1764.5):
        assert any(
            n == "Bi-214" and abs(E - E_expect) < 0.5 for n, E in all_E
        ), f"Bi-214 {E_expect} keV missing from RA226 forced clusters"


def test_ra226_forced_clusters_groups_coupled_by_nuclide():
    """F-121: все компоненты одного нуклида в кластере в одной группе
    (связанная подгонка по библиотечным интенсивностям).
    """
    for cluster in RA226_FORCED_CLUSTERS:
        nuc_to_groups: dict = {}
        for n, E, I, g in cluster["components"]:
            nuc_to_groups.setdefault(n, set()).add(g)
        for n, groups in nuc_to_groups.items():
            # ровно одна группа на нуклид
            assert len(groups) == 1, (
                f"{cluster['id']}: nuclide {n} split across "
                f"groups {groups} — must be single group"
            )
            g = next(iter(groups))
            assert g == n, (
                f"{cluster['id']}: nuclide {n} group label {g!r} "
                f"must equal nuclide name"
            )


def test_run_chain_forced_multiplets_u238_branch_emits_three_clusters():
    """F-121: при chain_dominance.u238=True пайплайн возвращает 3 кластера."""
    # synthetic spectrum: 4096 каналов, калибровка 1 keV/ch + 0
    class _FakeSpec:
        def __init__(self):
            self.counts = np.zeros(4096, dtype=np.float64)
            # Add Gaussians at the key Bi-214 / Pb-214 energies
            for E_kev, A in [
                (295.0, 5000.0), (352.0, 8000.0),
                (609.0, 12000.0), (665.0, 800.0),
                (768.0, 2000.0), (806.0, 500.0),
                (1120.0, 3000.0), (1238.0, 1200.0),
                (1378.0, 800.0), (1408.0, 500.0),
                (1730.0, 600.0), (1764.0, 3500.0),
            ]:
                ch = int(round(E_kev))
                sigma_ch = 12.0  # rough NaI sigma in channels
                ch_arr = np.arange(4096)
                self.counts += A * np.exp(
                    -((ch_arr - ch) / sigma_ch) ** 2 * 0.5
                ) / (sigma_ch * np.sqrt(2 * np.pi))
            # Add background
            self.counts += 50.0

        def energy_to_channel(self, E):
            return float(E)

        def channel_to_energy(self, ch):
            return float(ch)

    spec = _FakeSpec()
    fwhm_keV = lambda E: 12.0 * 2.355 if E > 0 else 5.0  # noqa: E731
    fwhm_ch = lambda ch: 12.0 * 2.355  # noqa: E731

    out = run_chain_forced_multiplets(
        spec, fwhm_ch, fwhm_keV,
        _FakeChainDom(u238=True),
        filename_isotope_hints=None,
        use_peak_image=False, detector_type="NaI",
    )
    # 3 cluster results expected
    assert len(out) == 3, f"expected 3 U-cluster results, got {len(out)}"


def test_run_chain_forced_multiplets_th232_still_returns_two():
    """F-121 не должен ломать F-118: при th232 dominance — 2 кластера M1/M2."""
    class _FakeSpec:
        def __init__(self):
            self.counts = np.zeros(4096, dtype=np.float64)
            # Th-232 Ac-228 lines
            for E_kev, A in [
                (860.6, 800.0), (911.0, 4000.0), (964.8, 800.0),
                (969.0, 2500.0),
                (1588.0, 500.0), (1620.0, 300.0), (1630.0, 300.0),
            ]:
                ch = int(round(E_kev))
                sigma_ch = 12.0
                ch_arr = np.arange(4096)
                self.counts += A * np.exp(
                    -((ch_arr - ch) / sigma_ch) ** 2 * 0.5
                ) / (sigma_ch * np.sqrt(2 * np.pi))
            self.counts += 30.0

        def energy_to_channel(self, E):
            return float(E)

        def channel_to_energy(self, ch):
            return float(ch)

    spec = _FakeSpec()
    fwhm_keV = lambda E: 12.0 * 2.355  # noqa: E731
    fwhm_ch = lambda ch: 12.0 * 2.355  # noqa: E731

    out = run_chain_forced_multiplets(
        spec, fwhm_ch, fwhm_keV,
        _FakeChainDom(th232=True),
        filename_isotope_hints=None,
        use_peak_image=False, detector_type="NaI",
    )
    assert len(out) == 2, f"expected 2 M-cluster results for Th-232, got {len(out)}"
