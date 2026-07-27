"""
BUG-9 (A) — LSRM `.spe` reader preserves the file's true channel count
regardless of value, with no hardcoded ceiling on N.

User clarification (2026-06-03):
    «уточнение 1024 канала в тех, что присылал. В файлах spe от других
    детекторов может быть другое число кратное 1024.»

This test synthesises minimal LSRM SpectraLine binary `.spe` files (CP-1251
header ending in `SPECTR=` + uint32-LE channel block) at multiple channel
counts that are multiples of 1024, and verifies that
`read_lsrm_spe(apply_energy_ceiling=False)` round-trips every channel.

The energy ceiling is intentionally disabled here so the test isolates the
"channel-count read from binary block" contract from the project's 3-MeV
truncation policy. A separate ceiling-on test asserts the ceiling trims
only the high-E tail (it never affects low channels or the binary read).

If a regression ever reintroduces a hardcoded 1024 / 1003 channel limit
(or any other fixed cap) in the reader, the parametrised cases for
N ∈ {2048, 3072, 4096, 8192} will fail loudly.
"""

from __future__ import annotations

import struct
from pathlib import Path

import numpy as np
import pytest

from gamma.io.lsrm_spe import read_lsrm_spe


# ---------------------------------------------------------------------------
# Helper — build a minimal valid LSRM SpectraLine .spe file in memory.
# Only the fields the reader strictly needs are emitted; everything else is
# optional per the format (see lsrm_spe.py module docstring).
# ---------------------------------------------------------------------------

def _build_minimal_lsrm_spe(
    path: Path,
    counts: np.ndarray,
    *,
    a0: float = 0.0,
    a1: float = 1.0,
    tlive: float = 1000.0,
    treal: float = 1000.0,
) -> None:
    """Write `counts` as a minimal LSRM SpectraLine binary .spe at `path`."""
    n = int(len(counts))
    header_lines = [
        "SHIFR=bug9-synth",
        "TYPE=Калибровка",
        "TLIVE={:.2f}".format(float(tlive)),
        "TREAL={:.2f}".format(float(treal)),
        # ENERGY line: degree marker + 7 coefficient slots. Linear with
        # the requested (a0, a1) is enough to exercise channel↔energy.
        "ENERGY=1,{:.6g},{:.6g},0,0,0,0,0".format(float(a0), float(a1)),
        "SPECTRSIZE={:d}".format(n),
    ]
    header_text = "".join(line + "\r\n" for line in header_lines)
    header_bytes = header_text.encode("cp1251")
    marker = b"SPECTR="
    binary = np.asarray(counts, dtype="<u4").tobytes()
    path.write_bytes(header_bytes + marker + binary)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

# Channel counts: every multiple of 1024 from 1×1024 up to 8×1024. These
# are the realistic spectrometer modes the user cited (1024-channel USB
# NaI, 2048/4096 mid-range, 8192 for AS2 Pro and high-resolution HPGe).
@pytest.mark.parametrize("n_channels", [1024, 2048, 3072, 4096, 8192])
def test_lsrm_spe_preserves_arbitrary_multiple_of_1024(
    tmp_path: Path, n_channels: int
) -> None:
    """Reader must keep all N channels for any N (multiple of 1024)."""
    # Synthesise a distinctive count pattern so a silent truncation would
    # show as a bogus sum / argmax. Position-encoded counts: counts[i] = i+1.
    counts = np.arange(1, n_channels + 1, dtype=np.uint32)
    f = tmp_path / f"synth_{n_channels}.spe"
    _build_minimal_lsrm_spe(f, counts, a0=0.0, a1=3.0)

    spec = read_lsrm_spe(str(f), apply_energy_ceiling=False)

    assert spec.n_channels_raw == n_channels, (
        f"reader trimmed channels at read time: n_channels_raw="
        f"{spec.n_channels_raw} ≠ written {n_channels}"
    )
    assert len(spec.counts) == n_channels, (
        f"counts array length mismatch: {len(spec.counts)} ≠ {n_channels}"
    )
    assert spec.n_channels == n_channels
    # First and last channel values survive intact (catches off-by-one
    # and silent end-of-block clipping).
    assert int(spec.counts[0]) == 1
    assert int(spec.counts[-1]) == n_channels
    # Sum invariant: 1+2+...+N = N(N+1)/2
    expected_sum = n_channels * (n_channels + 1) // 2
    assert int(spec.counts.sum()) == expected_sum


def test_lsrm_spe_ceiling_only_trims_high_E_tail(tmp_path: Path) -> None:
    """
    With apply_energy_ceiling=True (default), only the HIGH-energy tail
    is trimmed — the low-channel block is never touched even if a0<0.

    This guards against any future regression that conflates "trim
    above ceiling" with "trim outside calibrated range".
    """
    n_channels = 4096
    counts = np.arange(1, n_channels + 1, dtype=np.uint32)
    f = tmp_path / "synth_ceiling.spe"
    # a0 = -12 keV (drift left), a1 = 1.0 keV/ch → E(0) = -12,
    # E(4095) ≈ 4083 keV; ceiling at 3 MeV trims a known number of
    # high channels but never the low ones.
    _build_minimal_lsrm_spe(f, counts, a0=-12.0, a1=1.0)

    spec_full = read_lsrm_spe(str(f), apply_energy_ceiling=False)
    spec_ceil = read_lsrm_spe(str(f), apply_energy_ceiling=True)

    # Full read: every channel survives.
    assert spec_full.n_channels_raw == n_channels
    assert len(spec_full.counts) == n_channels

    # Ceiling read: low channels intact, only the high-E tail removed.
    assert spec_ceil.n_channels_raw == n_channels, (
        "n_channels_raw must reflect the binary block, not the trimmed view"
    )
    assert len(spec_ceil.counts) < n_channels, (
        "ceiling at 3 MeV must trim some channels of a 4083-keV spectrum"
    )
    # Low channels preserved bit-for-bit:
    np.testing.assert_array_equal(
        np.asarray(spec_ceil.counts[:50]),
        np.asarray(spec_full.counts[:50]),
    )
    # First channel = 1, second = 2 (the synthetic pattern).
    assert int(spec_ceil.counts[0]) == 1
    assert int(spec_ceil.counts[1]) == 2


def test_lsrm_spe_unaligned_binary_still_reads_full_channels(
    tmp_path: Path,
) -> None:
    """
    Defensive: if a writer mistakenly appends a stray byte after the
    final uint32, the reader takes floor(bytes/4) channels (existing
    contract) — there is no other hardcoded channel cap.
    """
    n_channels = 2048
    counts = np.arange(1, n_channels + 1, dtype=np.uint32)
    f = tmp_path / "synth_unaligned.spe"
    _build_minimal_lsrm_spe(f, counts, a0=0.0, a1=1.5)
    # Append 3 garbage bytes so the binary is no longer aligned.
    with f.open("ab") as fh:
        fh.write(b"\x00\x00\x00")

    spec = read_lsrm_spe(str(f), apply_energy_ceiling=False)
    # The reader uses len(binary)//4 → still 2048 full channels read.
    assert spec.n_channels_raw == n_channels
    assert len(spec.counts) == n_channels
    assert int(spec.counts[-1]) == n_channels
