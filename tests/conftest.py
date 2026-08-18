"""
v1.17.9.1 — pytest conftest для tests/.

Прежде чем тесты импортировали `gamma.*` через явный
`sys.path.insert(0, str(Path(__file__).parent / "scripts"))`,
поскольку тесты лежали в корне. После переезда test_*.py
в каталог tests/ путь до scripts/ изменился — этот conftest
автоматически добавляет ../scripts/ в sys.path для всех тестов
этой папки.

Также:
    - Устанавливает PYTHONIOENCODING=utf-8 (для Cyrillic в stdout).
    - Регистрирует custom marker @pytest.mark.slow.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

# ──────────────────────────────────────────────────────────────────
# Make scripts/ importable for `from gamma...` statements
# ──────────────────────────────────────────────────────────────────
HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

# Some tests still construct `HERE = Path(__file__).resolve().parent`
# expecting that to be the repo root (for accessing references/, demo_reports/,
# detectors/, data/ ...). To minimise rewrites, expose REPO_ROOT via env var:
os.environ.setdefault("GAMMA_SKILL_ROOT", str(REPO_ROOT))

# ──────────────────────────────────────────────────────────────────
# Force UTF-8 for stdout (Windows console fallback)
# ──────────────────────────────────────────────────────────────────
os.environ.setdefault("PYTHONIOENCODING", "utf-8")


def pytest_configure(config):
    """Register custom markers used in this test suite."""
    config.addinivalue_line(
        "markers",
        "slow: marks tests that take >5s (deselect by default)",
    )


# ──────────────────────────────────────────────────────────────────
# REL-02 (AUDIT_v2 §3 / §6 P0-1, 2026-06-06): autouse cache-reset
# fixture for module-global caches.
#
# CI runs `pytest -n auto` (xdist). Within one worker, tests share a
# process and therefore share every module-level cache. Prior to this
# fixture, isolation relied on each test individually remembering to
# call the right `reset_cache()` helpers — a discipline that produced
# the documented Th-232 chain-ratio regression (2.30× ↔ 7.63×,
# order-dependent).
#
# Reset is performed PRE-test (before the `yield`) so that the very
# first call into any module under test sees a clean cache regardless
# of whether the previous test honoured its own cleanup.
# ──────────────────────────────────────────────────────────────────

@pytest.fixture(autouse=True, scope="function")
def _reset_module_caches():
    """Reset module-global caches before every test (REL-02).

    Covered modules (AUDIT_v2 §3, §6 P0-1):

    * ``gamma.data.nuclide_library._CACHE`` — built-in nuclide JSON
      cache (REL-01 mutation point).
    * ``gamma.data.xrf_catalog._CACHE`` — X-ray fluorescence catalogue.
    * ``gamma.data.anchors._CACHE`` — anchor pattern registry.
    * ``gamma.physics.secondary_peaks._CATALOG_CACHE`` and
      ``_CATALOG_V2_CACHE`` — secondary-feature empirical catalogues.
    * ``gamma.physics.tcs_correction._LOOKUP_CACHE`` — TSCF lookup
      table.
    * ``gamma.calibration.efficiency_autoload.find_efr_file`` and
      ``load_efficiency_for_geometry`` ``lru_cache`` decorators.

    All resets are best-effort and tolerant of import failures so that
    a missing optional module never breaks the whole suite (e.g. when
    a fast unit test only exercises ``scripts.utils.*``).
    """
    # gamma.data.nuclide_library — REL-01 mutation point.
    try:
        from gamma.data import nuclide_library as _nl
        _nl.reset_cache()
    except Exception:  # noqa: BLE001 — best-effort cleanup
        pass

    # gamma.data.xrf_catalog
    try:
        from gamma.data import xrf_catalog as _xrf
        _xrf.reset_cache()
    except Exception:  # noqa: BLE001
        pass

    # gamma.data.anchors
    try:
        from gamma.data import anchors as _anchors
        _anchors.reset_cache()
    except Exception:  # noqa: BLE001
        pass

    # gamma.physics.secondary_peaks
    try:
        from gamma.physics import secondary_peaks as _sp
        _sp.reset_catalog_cache()
        _sp.reset_catalog_v2_cache()
    except Exception:  # noqa: BLE001
        pass

    # gamma.physics.tcs_correction (reset helper added in this wave)
    try:
        from gamma.physics import tcs_correction as _tcs
        _tcs.reset_lookup_cache()
    except Exception:  # noqa: BLE001
        pass

    # gamma.calibration.efficiency_autoload — lru_cache decorators.
    try:
        from gamma.calibration import efficiency_autoload as _eff
        _eff.find_efr_file.cache_clear()
        _eff.load_efficiency_for_geometry.cache_clear()
    except Exception:  # noqa: BLE001
        pass

    yield


# ──────────────────────────────────────────────────────────────────
# DEEP-03 (2026-06-05): basename-sort hook REMOVED — it was a no-op fix.
#
# History: v1.17.9.2 added a `pytest_collection_modifyitems` hook that
# sorted collected items by file basename, on the theory that a few
# numerical tests in step07/step09 depended on NOT being preceded by
# step08 multiplets that mutate "scipy fit caches / FWHM provider".
#
# Why the sort never helped:
#   CI runs `pytest -n auto` (pytest-xdist). xdist's default scheduler
#   (`--dist load`) distributes tests to workers BY COLLECTION INDEX,
#   then each worker executes its own slice. Reordering the global
#   collection list does NOT control which tests land on a worker before
#   the previously-quarantined test, nor in what relative order — the
#   basename-sort is therefore a no-op under `-n auto` BY DESIGN.
#
# Real mechanism of the F-410 flake (verifier audit, corrected framing):
#   The previously-quarantined test asserts gauss_sigma_keV to a 1e-6 keV
#   tolerance (test_linematch_writer_normalisation.py). Under parallel
#   workers the residual variation is FP/BLAS non-determinism at that
#   tolerance combined with order-dependence on internal mutated state —
#   NOT a "scipy fit cache" that the basename-sort could have repaired.
#   "scipy fit caches" was an incorrect wave-1 attribution.
#
# Fix: the hook is removed entirely (it had no other purpose). If a
# numerical test ever again proves genuinely order-sensitive under xdist,
# the correct levers are (a) eliminate the shared mutable state in that
# test, (b) tighten/relax the numerical tolerance to the physically
# meaningful bound, or (c) serialise the group via `--dist loadgroup`
# with an `@pytest.mark.xdist_group(...)` mark — NOT a collection-order
# sort, which xdist ignores.
# ──────────────────────────────────────────────────────────────────
