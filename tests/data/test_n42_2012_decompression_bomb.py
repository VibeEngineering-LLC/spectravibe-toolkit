"""SEC-01 — N42-2012 CountedZeroes decompression-bomb hardening.

`scripts/gamma/io/n42_2012.py:_decode_channel_data` previously had an
unbounded `out.extend([0] * max(1, n))` for the `CountedZeroes`
compression branch. An adversary supplying a <100-byte XML snippet
whose token sequence encodes a `0 2000000000` run could force the
process to allocate ~16 GB of zeros (Python ints, sizeof ~28 B), OOM
the analyst's workstation, or kill the staged pipeline silently.

After the SEC-01 fix the cumulative channel count is bounded at
1e7 channels (three orders of magnitude above realistic NaI/HPGe
upper bounds: NaI is typically 1024-16384, HPGe 8192-65536).
A breach raises `ValueError` with the marker substring
"exceeds bound" so the upstream readers can fail-loud rather than
swallowing the OOM.

Tests (3):

  1. `test_counted_zeroes_bounded_normal_case` — happy path. A
     4096-channel decompression with a moderate `0 N` run succeeds
     and returns a numpy array of the expected length.

  2. `test_counted_zeroes_raises_on_bomb` — pathological input. A
     handcrafted text token sequence whose `0` token is followed by
     `2 * 10**7` (just above the cap) must raise
     `ValueError` mentioning "exceeds bound".

  3. `test_counted_zeroes_partial_decode_state` — partial-state
     contract. When the bound is exceeded mid-decode, the function
     must raise rather than return a truncated array masquerading
     as a successful decode. The caller relies on the exception to
     fail the spectrum read; a silently-truncated array would be
     interpreted as a legitimate (just short) spectrum.

Red-without-fix evidence: `_tmp/red_sec01_parent_20260606.txt`
captures parent SHA e2bddbb pytest output with the bomb test
hanging or OOM-ing (Python kills it; reproducibility achieved via
the cap-just-exceeded value of 2e7 — enough to overshoot 1e7 but
not so large that the test itself OOMs the host).
Post-fix green: `_tmp/green_sec01_post_20260606.txt`.
"""

from __future__ import annotations

import numpy as np
import pytest

from gamma.io.n42_2012 import _decode_channel_data


def test_counted_zeroes_bounded_normal_case():
    """Happy path: a moderate CountedZeroes payload decompresses correctly."""
    # 100 explicit values then a run of 500 zeros, encoded as
    # "v1 v2 ... v100 0 500".
    explicit = list(range(1, 101))  # 100 non-zero ints
    run_length = 500
    text = " ".join(str(v) for v in explicit) + f" 0 {run_length}"

    out = _decode_channel_data(text, "CountedZeroes")

    assert isinstance(out, np.ndarray)
    assert out.dtype == np.int64
    # 100 explicit + 500 expanded zeros = 600 total.
    assert len(out) == 100 + run_length
    # First 100 values preserved.
    assert list(out[:100]) == explicit
    # Run of zeros at the tail.
    assert int(out[100:].sum()) == 0


def test_counted_zeroes_raises_on_bomb():
    """Pathological: `0 N` with N just above the 1e7 cap must raise.

    Using N = 2e7 — overshoots the 1e7 cap by 1e7 but is small enough
    that the *raise* path itself does not OOM the test host. The
    pre-fix code WOULD attempt to allocate a 20-million-element list of
    Python ints (~560 MB peak); the post-fix code raises before any
    expansion.
    """
    # Single sentinel value, then `0 20000000` run.
    text = "1 0 20000000"

    with pytest.raises(ValueError, match=r"exceeds bound"):
        _decode_channel_data(text, "CountedZeroes")


def test_counted_zeroes_partial_decode_state():
    """When bound is breached mid-decode, do NOT return a truncated array.

    Contract: the function must `raise`, not silently truncate to 1e7
    channels. A truncated return would masquerade as a legitimate (just
    short) spectrum and propagate downstream.
    """
    # Two valid bins, then a bomb run, then more data that the
    # implementation would never reach.
    text = "5 7 0 20000000 9 11 0 3"

    with pytest.raises(ValueError, match=r"exceeds bound"):
        _decode_channel_data(text, "CountedZeroes")
