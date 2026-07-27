"""
BUG-9 / v1.18.32 — post-read `trim_to_working_energy` helper.

The reader default was flipped to lossless in v1.18.32 (no automatic
3 MeV trim). For code that genuinely needs the project-scope cut at
analysis stage, the new helper `gamma.spectrum.trim_to_working_energy`
performs the same drop, but at the call site — visible, auditable, and
mutating the Spectrum in place.

This test pins the helper's behaviour against the reader-stage opt-in
trim (`apply_energy_ceiling=True`) on a real Gamma-1S 1024-ch fixture.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from gamma.io.lsrm_spe import read_lsrm_spe  # noqa: E402
from gamma.spectrum import (  # noqa: E402
    ENERGY_CEILING_KEV,
    trim_to_working_energy,
)


FIXTURE = (
    "detectors/Gamma-1S/reference_spectra/reference_kits/"
    "Marinelli_1L/Th-232/Th232_420-7-17_Маринелли_0cm.spe"
)


def _fx() -> Path:
    return (ROOT / FIXTURE).resolve()


def test_trim_helper_matches_reader_stage_trim() -> None:
    """
    `trim_to_working_energy(spec)` on a lossless-read spectrum must
    produce the same channel count and `energy_max_keV_kept` as the
    reader-stage `apply_energy_ceiling=True` path. Different code path,
    identical answer — the new helper is a drop-in replacement that is
    just easier to spot in a diff.
    """
    f = _fx()
    if not f.is_file():
        pytest.skip(f"fixture not present: {f}")

    spec_full = read_lsrm_spe(str(f))                         # lossless
    spec_trim = read_lsrm_spe(str(f), apply_energy_ceiling=True)

    # Trim the lossless copy at analysis stage.
    spec_helper = read_lsrm_spe(str(f))
    trim_to_working_energy(spec_helper)

    assert spec_full.n_channels == 1024
    assert spec_helper.n_channels == spec_trim.n_channels, (
        f"helper trimmed to {spec_helper.n_channels} channels but "
        f"reader-stage trim produced {spec_trim.n_channels}"
    )
    assert spec_helper.energy_max_keV_kept == pytest.approx(
        spec_trim.energy_max_keV_kept, rel=1e-9, abs=1e-6
    )
    # n_channels_raw is preserved by the helper (audit trail).
    assert spec_helper.n_channels_raw == 1024


def test_trim_helper_custom_ceiling() -> None:
    """A custom ceiling truncates more aggressively than the default."""
    f = _fx()
    if not f.is_file():
        pytest.skip(f"fixture not present: {f}")

    spec = read_lsrm_spe(str(f))
    n0 = spec.n_channels
    trim_to_working_energy(spec, max_keV=1500.0)
    assert spec.n_channels < n0, (
        f"custom ceiling=1500 keV must keep fewer channels than default "
        f"(got {spec.n_channels} vs initial {n0})"
    )
    assert spec.energy_max_keV_kept is not None
    assert spec.energy_max_keV_kept <= 1500.0


def test_trim_helper_noop_without_calibration() -> None:
    """If `spec.energy_cal is None` the helper is a no-op."""
    f = _fx()
    if not f.is_file():
        pytest.skip(f"fixture not present: {f}")

    spec = read_lsrm_spe(str(f))
    spec.energy_cal = None
    n_before = int(len(spec.counts))
    trim_to_working_energy(spec, max_keV=ENERGY_CEILING_KEV)
    assert int(len(spec.counts)) == n_before
