"""DEEP-03 meta-guard (Project #5 wave 2 P1-3).

Reds if the wave-1 F-410 *wrong-framing* fix tries to come back:

  1. The no-op collection-order basename-sort in ``tests/conftest.py``
     (``pytest_collection_modifyitems`` with ``items.sort(...)``). It was
     removed in DEEP-03 because pytest-xdist ``-n auto`` distributes tests
     by collection INDEX, so reordering the collection list is a no-op by
     design — it only ever masked the real order-dependent test defect.

  2. The F-410 ``@pytest.mark.xfail`` quarantine on
     ``test_bug34_w2_w3_gauss_sigma_through_staged_pipeline_th232`` in
     ``tests/snapshot/test_linematch_writer_normalisation.py``. The flake
     was fixed at source (order-independent test fakes), so the xfail must
     stay reverted.

If either regression lands, the corresponding assertion below fails,
flagging that the wrong "scipy fit caches" framing is being reintroduced
instead of the real order-independence fix.

Cite (file:line, 2026-06-05):
  - tests/conftest.py — basename-sort hook removed (was lines 49-67).
  - tests/snapshot/test_linematch_writer_normalisation.py:291 — xfail
    decorator removed (was 291-302).
  - tests/step08_multiplets/test_deconvolve.py — _FakeNuclideId /
    _FakeIdent converted to dataclasses (the real fix).
"""
from __future__ import annotations

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
CONFTEST = REPO_ROOT / "tests" / "conftest.py"
LINEMATCH = (
    REPO_ROOT / "tests" / "snapshot"
    / "test_linematch_writer_normalisation.py"
)
LINEMATCH_TEST = "test_bug34_w2_w3_gauss_sigma_through_staged_pipeline_th232"


def test_conftest_has_no_basename_collection_sort():
    """conftest must NOT reintroduce the no-op collection-order sort.

    Guards DEEP-03: ``items.sort(`` inside a
    ``pytest_collection_modifyitems`` hook is the exact wrong-framing fix
    that xdist ignores by design.
    """
    assert CONFTEST.is_file(), f"missing conftest: {CONFTEST}"
    src = CONFTEST.read_text(encoding="utf-8")

    assert "items.sort(" not in src, (
        "tests/conftest.py reintroduced `items.sort(` — the no-op "
        "collection-order basename-sort removed in DEEP-03. xdist `-n auto` "
        "distributes by collection INDEX, so this never fixes order "
        "dependence; it only masks real test defects. Fix the order-"
        "dependent test at source instead (see DEEP-03)."
    )
    # Match the hook DEFINITION, not the explanatory comment that documents
    # why it was removed (the comment legitimately names the hook).
    assert not re.search(r"^\s*def\s+pytest_collection_modifyitems\b", src, re.M), (
        "tests/conftest.py reintroduced a `pytest_collection_modifyitems` "
        "hook. DEEP-03 removed it as a no-op masking fix; if you genuinely "
        "need collection reordering, document why xdist does not defeat it."
    )


def test_linematch_writer_test_is_not_xfail_quarantined():
    """The F-410 quarantine on the linematch writer test must stay reverted.

    Guards DEEP-03 Step 4: the test was un-quarantined after the real,
    order-independent fix made ``pytest -n auto`` green x3.
    """
    assert LINEMATCH.is_file(), f"missing test file: {LINEMATCH}"
    src = LINEMATCH.read_text(encoding="utf-8")

    assert LINEMATCH_TEST in src, (
        f"meta-guard stale: expected to find `{LINEMATCH_TEST}` in "
        f"{LINEMATCH.name}. If the test was renamed, update this guard."
    )

    # Locate the test definition and inspect the lines immediately above it
    # for an xfail decorator (the F-410 quarantine form).
    lines = src.splitlines()
    def_idx = next(
        (i for i, ln in enumerate(lines)
         if re.match(rf"\s*def\s+{re.escape(LINEMATCH_TEST)}\s*\(", ln)),
        None,
    )
    assert def_idx is not None, (
        f"could not locate `def {LINEMATCH_TEST}(` in {LINEMATCH.name}"
    )

    # Scan the decorator block directly above the def: walk upward over
    # decorator lines (`@...`), their continuation lines, and comments.
    # An xfail mark anywhere in that contiguous decorator block is the
    # F-410 quarantine regression.
    window = "\n".join(lines[max(0, def_idx - 25):def_idx])
    assert "mark.xfail" not in window, (
        f"F-410 xfail quarantine reappeared above {LINEMATCH_TEST} in "
        f"{LINEMATCH.name}. DEEP-03 reverted it after the real order-"
        f"independence fix. Do NOT re-quarantine; fix the root cause."
    )


def test_meta_guard_self_check_would_fire():
    """Self-check: the guard predicates actually detect the regressions.

    Confirms the assertions above are not vacuous — i.e. they WOULD red if
    the basename-sort or the xfail decorator were present. This is the
    'reds without fix' half of the DoD, exercised without mutating files.
    """
    # (a) basename-sort predicate fires on the historical wrong-framing line.
    historical_sort = (
        "    items.sort(key=lambda it: os.path.basename(str(it.fspath)))"
    )
    assert "items.sort(" in historical_sort

    # (b) xfail predicate fires on the historical F-410 decorator text.
    historical_xfail = "@pytest.mark.xfail(\n    reason=(...),\n    strict=False,\n)"
    assert "mark.xfail" in historical_xfail
