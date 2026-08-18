"""
BUG-9 (B) — LSRM `.spe` reader preserves channels whose calibrated
energy is < 0 (calibration drift signature).

User clarification (2026-06-03):
    «возможно при калибровке спектр сдвигается влево и нулевые каналы
    уходят в минус по энергии и отбрасываются. Отрицательная энергия в
    первых каналах не баг, свидетельсво о дрейфе калибровки
    спектрометра.»

Tests:
  1. Synthetic .spe with a0 = -12 keV — every channel survives, the
     diagnostic flag `calibration_drift_left=True` is set, and the
     number of negative-E channels is recorded.
  2. Synthetic .spe with a0 = +5 keV (no drift) — flag is False and
     no leading channels are flagged.
  3. Counts of the first few channels match the synthetic pattern
     bit-for-bit (no silent skip of "noise" channels).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from gamma.io.lsrm_spe import read_lsrm_spe


def _build_spe(
    path: Path,
    counts: np.ndarray,
    *,
    a0: float,
    a1: float = 1.0,
) -> None:
    """Minimal LSRM SpectraLine binary .spe writer (test helper)."""
    n = int(len(counts))
    header_lines = [
        "SHIFR=bug9-drift-synth",
        "TYPE=Калибровка",
        "TLIVE=1000.00",
        "TREAL=1000.00",
        "ENERGY=1,{:.6g},{:.6g},0,0,0,0,0".format(float(a0), float(a1)),
        "SPECTRSIZE={:d}".format(n),
    ]
    header_bytes = "".join(line + "\r\n" for line in header_lines).encode("cp1251")
    marker = b"SPECTR="
    binary = np.asarray(counts, dtype="<u4").tobytes()
    path.write_bytes(header_bytes + marker + binary)


def test_negative_a0_preserves_all_channels(tmp_path: Path) -> None:
    """
    With a0 = -12.36 keV (close to the real sopka-mica spectrometer
    drift), the first ~12 channels have E < 0. They MUST be preserved
    in `counts`, and the diagnostic flag MUST be set.
    """
    n_channels = 4096
    # Distinctive pattern: counts[i] = 10 + i, so a silent drop of
    # the first few channels would shift every assertion.
    counts = np.arange(10, 10 + n_channels, dtype=np.uint32)
    f = tmp_path / "drift_left.spe"
    _build_spe(f, counts, a0=-12.36, a1=1.0)

    # Read both with and without ceiling — both must keep low channels.
    for apply_ceiling in (False, True):
        spec = read_lsrm_spe(str(f), apply_energy_ceiling=apply_ceiling)

        # No silent low-end trim:
        assert int(spec.counts[0]) == 10, (
            f"first channel dropped under apply_energy_ceiling={apply_ceiling}: "
            f"counts[0]={int(spec.counts[0])}, expected 10"
        )
        assert int(spec.counts[1]) == 11
        assert int(spec.counts[2]) == 12

        # E(0) is still negative:
        assert spec.channel_to_energy(0) == pytest.approx(-12.36, rel=1e-6)
        assert spec.channel_to_energy(0) < 0.0

        # n_channels_raw reflects the binary block in full, regardless
        # of ceiling behaviour.
        assert spec.n_channels_raw == n_channels

        # Diagnostic flag is set, with the correct a0 stashed and the
        # count of negative-E channels populated.
        assert spec.extras.get("calibration_drift_left") is True
        assert spec.extras.get("calibration_drift_a0_keV") == pytest.approx(
            -12.36, rel=1e-6
        )
        n_neg = spec.extras.get("calibration_drift_neg_energy_channels")
        assert isinstance(n_neg, int)
        # With a0=-12.36, a1=1.0: E(ch)=ch-12.36 → first 13 channels (0..12)
        # have E<0. (ch=12 → -0.36 < 0; ch=13 → 0.64 ≥ 0.)
        assert n_neg == 13, (
            f"expected 13 leading channels with E<0 (a0=-12.36, a1=1.0); "
            f"reader reported {n_neg}"
        )


def test_positive_a0_does_not_flag_drift(tmp_path: Path) -> None:
    """No calibration drift → flag is False, no neg-E channels."""
    n_channels = 1024
    counts = np.arange(1, n_channels + 1, dtype=np.uint32)
    f = tmp_path / "no_drift.spe"
    _build_spe(f, counts, a0=5.0, a1=3.0)

    spec = read_lsrm_spe(str(f), apply_energy_ceiling=False)

    assert spec.extras.get("calibration_drift_left") is False
    assert "calibration_drift_neg_energy_channels" not in spec.extras
    # First channel value preserved:
    assert int(spec.counts[0]) == 1
    # All energies non-negative:
    assert spec.channel_to_energy(0) == pytest.approx(5.0)


def test_n_channels_kept_equals_len_counts_when_no_ceiling(
    tmp_path: Path,
) -> None:
    """
    Invariant: with the ceiling off, n_channels == n_channels_raw ==
    len(counts) regardless of a0 sign. Guards against any future
    "implicit low-channel trim" regression.
    """
    n_channels = 8192
    counts = np.ones(n_channels, dtype=np.uint32)
    f = tmp_path / "invariant.spe"
    _build_spe(f, counts, a0=-30.0, a1=0.4)  # heavy left drift

    spec = read_lsrm_spe(str(f), apply_energy_ceiling=False)
    assert spec.n_channels_raw == n_channels
    assert spec.n_channels == n_channels
    assert len(spec.counts) == n_channels
    assert spec.extras["calibration_drift_left"] is True
