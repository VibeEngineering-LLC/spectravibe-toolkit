"""REL-01 / REL-02 (AUDIT_v2 §3, §6 P0-1): nuclide library cache isolation.

Two related defects:

REL-01 — ``scripts/gamma/data/nuclide_library.load_external_library()``
mutated the module-global ``_CACHE`` in place via the reference returned
by ``_load()``. Any caller holding a prior ``_load()`` reference saw the
mutated state. Under ``pytest -n auto`` (the actual CI command), tests
within one worker contaminated each other, producing a documented
order-dependent regression in the Th-232 chain ratio (2.30x vs 7.63x).

REL-02 — ``tests/conftest.py`` had no autouse fixture resetting
module-global caches between tests. Isolation relied on each individual
test remembering to call the right ``reset_cache()`` helpers. With
``pytest -n auto`` this produced order-dependent pass/fail outcomes.

Tests in this file demonstrate red-without-fix and green-with-fix.

Inventory (5 tests, synced to actual function names per censor LOW-1
finding 2026-06-06 in ACCEPT-WITH-FIX verdict envelope 0c7075a2):

REL-01 (snapshot semantic of ``load_external_library``):

* ``test_load_external_library_does_not_mutate_prior_reference`` —
  override branch (``merge_mode='override'``). Calls the *real*
  ``load_external_library`` against monkey-patched Lsrm reader/merger
  helpers so no .lib fixture is needed on disk.
* ``test_load_external_library_supplement_does_not_mutate_prior_reference`` —
  supplement branch (``merge_mode='supplement'``). Same snapshot contract.

REL-02 (autouse cache-reset fixture in ``tests/conftest.py``):

* ``test_rel02_autouse_fixture_is_enrolled`` — deterministic backstop
  under any xdist scheduler. Pre-fix: fixture doesn't exist ⇒
  ``request.fixturenames`` lacks ``_reset_module_caches`` ⇒ FAIL.
* ``test_rel02_autouse_fixture_resets_nuclide_library_cache`` — asserts
  ``nuclide_library._CACHE is None`` at test entry (relies on fixture
  having reset it before ``yield``).
* ``test_rel02_autouse_fixture_resets_secondary_caches`` — asserts the
  fixture also reset ``xrf_catalog._CACHE``, ``anchors._CACHE``,
  ``secondary_peaks._CATALOG_CACHE`` + ``_CATALOG_V2_CACHE``,
  ``tcs_correction._LOOKUP_CACHE``, and the
  ``gamma.calibration.efficiency_autoload`` ``functools.lru_cache``
  wrappers (``find_efr_file``, ``load_efficiency_for_geometry``).

Rationale for introspection-based REL-02 (instead of paired
``test_a_pollutes`` / ``test_b_detects_leak``): under ``pytest -n auto``
the xdist scheduler (``--dist load``) distributes tests across workers
by collection index — sibling tests in the same file are NOT guaranteed
to land on the same worker, so a pollute/detect pair would be a flaky
leak detector. The enrollment + per-fixture-effect assertion approach
is robust under any scheduler. See lines 167-185 for the inline rationale
block.
"""

from __future__ import annotations


# ----------------------------------------------------------------------
# REL-01: load_external_library must not mutate a prior reference.
# ----------------------------------------------------------------------

def test_load_external_library_does_not_mutate_prior_reference(monkeypatch):
    """REL-01: ``load_external_library`` must not mutate a prior
    ``_load()`` reference.

    Pre-fix code (buggy):
        current = _load()
        current.update(external)              # mutates _CACHE in place

    Post-fix code (correct):
        _CACHE = {**_CACHE, **external}       # new dict, snapshot kept
    """
    from gamma.data import nuclide_library
    from gamma.io import lsrm_library as lsrm_lib_mod

    # Force a clean baseline so we know what to expect.
    nuclide_library.reset_cache()

    # Get a reference to the cache state BEFORE the external load.
    snapshot_ref = nuclide_library._load()
    snapshot_keys_before = set(snapshot_ref.keys())
    # Sanity: snapshot is non-empty (the bundled library has entries).
    assert len(snapshot_keys_before) > 0
    assert "FAKE-RELE01" not in snapshot_keys_before

    # Monkey-patch the Lsrm reader/merger pair inside the
    # nuclide_library module so we don't need a real .lib fixture.
    # load_external_library imports them locally via:
    #     from gamma.io.lsrm_library import (
    #         read_lsrm_library, merge_lsrm_library_into_internal,
    #     )
    # so we patch the source module attributes.
    class _FakeLib:  # placeholder Lsrm object; never inspected
        pass

    fake_lib_obj = _FakeLib()
    fake_external = {
        "FAKE-RELE01": {
            "T_half_s": 1.0,
            "lines": [[100.0, 1.0, 0.0]],
        },
    }

    monkeypatch.setattr(
        lsrm_lib_mod, "read_lsrm_library",
        lambda path: fake_lib_obj,
    )
    monkeypatch.setattr(
        lsrm_lib_mod, "merge_lsrm_library_into_internal",
        lambda lib, include_xrays=False, include_unused=False: dict(fake_external),
    )

    # Invoke the REAL function under test.
    added = nuclide_library.load_external_library(
        "ignored-path.lib",
        merge_mode="override",
        include_xrays=False,
        split_chains=False,  # bypass chain decomposition (irrelevant here)
    )
    assert added == 1, f"expected 1 entry added, got {added}"

    # Verify snapshot_ref is still pristine — this is the core REL-01
    # assertion. Pre-fix the in-place .update() would have leaked the
    # FAKE-RELE01 key into snapshot_ref.
    snapshot_keys_after = set(snapshot_ref.keys())
    leaked = snapshot_keys_after - snapshot_keys_before
    assert snapshot_keys_after == snapshot_keys_before, (
        "REL-01: load_external_library mutated a prior _load() reference. "
        f"Leaked keys: {sorted(leaked)}"
    )
    assert "FAKE-RELE01" not in snapshot_ref, (
        "REL-01: FAKE-RELE01 leaked into the prior cache reference via "
        "in-place mutation."
    )

    # Sanity: a fresh _load() now DOES see the merged entry.
    fresh_view = nuclide_library._load()
    assert "FAKE-RELE01" in fresh_view, (
        "REL-01: post-merge cache must contain the new entry on a fresh "
        f"_load(); got {len(fresh_view)} entries without it."
    )


def test_load_external_library_supplement_does_not_mutate_prior_reference(
        monkeypatch):
    """REL-01: same contract for ``merge_mode='supplement'`` branch."""
    from gamma.data import nuclide_library
    from gamma.io import lsrm_library as lsrm_lib_mod

    nuclide_library.reset_cache()
    snapshot_ref = nuclide_library._load()
    snapshot_keys_before = set(snapshot_ref.keys())
    assert "FAKE-RELE01-SUP" not in snapshot_keys_before

    fake_external = {
        "FAKE-RELE01-SUP": {
            "T_half_s": 1.0,
            "lines": [[200.0, 2.0, 0.0]],
        },
    }
    monkeypatch.setattr(
        lsrm_lib_mod, "read_lsrm_library", lambda path: object(),
    )
    monkeypatch.setattr(
        lsrm_lib_mod, "merge_lsrm_library_into_internal",
        lambda lib, include_xrays=False, include_unused=False: dict(fake_external),
    )

    added = nuclide_library.load_external_library(
        "ignored-path.lib",
        merge_mode="supplement",
        include_xrays=False,
        split_chains=False,
    )
    assert added == 1, f"expected 1 entry added in supplement mode, got {added}"

    snapshot_keys_after = set(snapshot_ref.keys())
    leaked = snapshot_keys_after - snapshot_keys_before
    assert snapshot_keys_after == snapshot_keys_before, (
        "REL-01: supplement-mode load_external_library mutated a prior "
        f"_load() reference. Leaked keys: {sorted(leaked)}"
    )
    assert "FAKE-RELE01-SUP" not in snapshot_ref, (
        "REL-01: FAKE-RELE01-SUP leaked into the prior cache reference "
        "via in-place mutation (supplement mode)."
    )

    fresh_view = nuclide_library._load()
    assert "FAKE-RELE01-SUP" in fresh_view


# ----------------------------------------------------------------------
# REL-02: autouse fixture in tests/conftest.py resets module-global
# caches before each test.
#
# Under ``pytest -n auto`` the xdist scheduler (--dist load, default)
# distributes tests by collection index — tests from sibling functions
# in the same file are NOT guaranteed to land on the same worker, so a
# pure "test_a pollutes / test_b detects leak" pair would be a flaky
# detector under -n auto.
#
# Instead, we exercise the autouse fixture's contract directly via
# pytest's ``request`` introspection: the fixture must (a) be enrolled
# into every test in the suite (``autouse=True``) and (b) actually
# reset ``nuclide_library._CACHE`` to ``None`` before each test runs.
#
# This is robust under any xdist scheduler because it does NOT rely on
# cross-test ordering — each test self-verifies the fixture ran
# correctly for itself.
# ----------------------------------------------------------------------

REL02_AUTOUSE_FIXTURE_NAME = "_reset_module_caches"
"""Name of the autouse fixture in tests/conftest.py that REL-02 adds."""


def test_rel02_autouse_fixture_is_enrolled(request):
    """REL-02: the autouse cache-reset fixture must be enrolled into
    every test.

    Pre-fix: no autouse fixture exists in ``tests/conftest.py`` ⇒
    ``REL02_AUTOUSE_FIXTURE_NAME`` is not in ``request.fixturenames``
    and this test FAILS.

    Post-fix: the fixture is autouse-enrolled and present here.
    """
    assert REL02_AUTOUSE_FIXTURE_NAME in request.fixturenames, (
        f"REL-02: autouse fixture '{REL02_AUTOUSE_FIXTURE_NAME}' is not "
        f"enrolled into this test. Got fixturenames: "
        f"{sorted(request.fixturenames)}. The autouse cache-reset "
        f"fixture in tests/conftest.py is missing."
    )


def test_rel02_autouse_fixture_resets_nuclide_library_cache():
    """REL-02: the autouse fixture must leave ``nuclide_library._CACHE``
    at ``None`` at test start.

    Pre-fix: ``_CACHE`` carries whatever the previous test (or process
    bootstrap) left in it. If a previous test populated the cache and
    no reset fixture exists, ``_CACHE`` is not ``None`` at the start of
    this test and the assertion FAILS.

    Post-fix: the autouse fixture sets ``_CACHE = None`` before each
    test, so this assertion holds.
    """
    from gamma.data import nuclide_library

    # We deliberately do NOT call reset_cache() inside the test — we
    # rely on the autouse fixture to do that. If the fixture is
    # missing, the cache may carry state from an earlier test.
    assert nuclide_library._CACHE is None, (
        "REL-02: nuclide_library._CACHE is not None at test entry. "
        "The autouse cache-reset fixture in tests/conftest.py is "
        "missing or did not reset this module."
    )


def test_rel02_autouse_fixture_resets_secondary_caches():
    """REL-02: the autouse fixture must also reset adjacent
    module-global caches enumerated in AUDIT_v2 §6 P0-1:

    * ``gamma.data.xrf_catalog._CACHE``
    * ``gamma.data.anchors._CACHE``
    * ``gamma.physics.secondary_peaks._CATALOG_CACHE`` and
      ``_CATALOG_V2_CACHE``
    * ``gamma.physics.tcs_correction._LOOKUP_CACHE``
    * ``gamma.calibration.efficiency_autoload.find_efr_file`` and
      ``load_efficiency_for_geometry`` ``functools.lru_cache`` wrappers
      (added 2026-06-06 per censor LOW-2 in envelope 0c7075a2 —
      previously the fixture cleared these but no test asserted it,
      so a silent ``except Exception`` in the fixture could swallow a
      broken/renamed helper unnoticed).
    """
    from gamma.data import xrf_catalog, anchors
    from gamma.physics import secondary_peaks, tcs_correction
    from gamma.calibration import efficiency_autoload

    assert xrf_catalog._CACHE is None, (
        "REL-02: xrf_catalog._CACHE not reset by autouse fixture."
    )
    assert anchors._CACHE is None, (
        "REL-02: anchors._CACHE not reset by autouse fixture."
    )
    assert secondary_peaks._CATALOG_CACHE is None, (
        "REL-02: secondary_peaks._CATALOG_CACHE not reset by autouse "
        "fixture."
    )
    assert secondary_peaks._CATALOG_V2_CACHE is None, (
        "REL-02: secondary_peaks._CATALOG_V2_CACHE not reset by autouse "
        "fixture."
    )
    assert tcs_correction._LOOKUP_CACHE is None, (
        "REL-02: tcs_correction._LOOKUP_CACHE not reset by autouse "
        "fixture."
    )
    # LOW-2 (censor 0c7075a2): explicit assert on efficiency_autoload
    # lru_cache wrappers. If the fixture's
    # `efficiency_autoload.find_efr_file.cache_clear()` call raises
    # AttributeError (e.g. helper renamed) the broad `except Exception`
    # in tests/conftest.py:48-130 swallows it — without this assert the
    # gap is invisible.
    assert efficiency_autoload.find_efr_file.cache_info().currsize == 0, (
        "REL-02: efficiency_autoload.find_efr_file lru_cache not cleared "
        "by autouse fixture (currsize != 0)."
    )
    assert (
        efficiency_autoload.load_efficiency_for_geometry.cache_info().currsize
        == 0
    ), (
        "REL-02: efficiency_autoload.load_efficiency_for_geometry "
        "lru_cache not cleared by autouse fixture (currsize != 0)."
    )
