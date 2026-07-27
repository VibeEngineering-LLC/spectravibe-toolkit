"""
BUG-9 (real-loss reproducer) — v1.18.32 / 2026-06-03.

The original BUG-9 report was about a 1024-channel LSRM `.spe` file
shrinking to 1003 channels during conversion. Forensic analysis traced
this to the reader-stage 3 MeV energy ceiling (`apply_energy_ceiling=True`
default in `read_lsrm_spe` / `read_atomspectra_xml`): when the stored
calibration has a slightly negative `a0` (e.g. `a0≈-8 keV`, `gain≈3 keV/ch`,
typical for a Gamma-1S NaI head that has drifted left), channel 1023 maps
to E ≈ 3063 keV, just above the ceiling, so the reader silently drops the
last 21 channels.

User quote (verbatim Russian, the canonical trigger for this fix):
    «В файлах лсрм 1024 канала, в процессе конвертации идет потеря.»

Fix in v1.18.32: the default of `apply_energy_ceiling` flipped to False
for all readers. This test locks the new contract by asserting that the
user-facing default path returns exactly the binary block's channel count
on every real-world reproducer fixture cited in the bug report.

If a future refactor reintroduces the lossy default (or any other silent
channel trim at read time), each parametrized case will fail loudly.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Make `scripts/` importable regardless of pytest invocation directory.
ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from gamma.io.lsrm_spe import read_lsrm_spe  # noqa: E402
from gamma.io.readers import read_spectrum  # noqa: E402


# ---------------------------------------------------------------------------
# Real Gamma-1S 1024-channel reproducer fixtures.
# Each tuple: (relative_path, expected_channel_count).
# ---------------------------------------------------------------------------
_FIXTURES: list[tuple[str, int]] = [
    (
        "detectors/Gamma-1S/reference_spectra/reference_kits/"
        "Marinelli_1L/Th-232/Th232_420-7-17_Маринелли_0cm.spe",
        1024,
    ),
    (
        "detectors/Gamma-1S/reference_spectra/archive/"
        "Cs137_420-7-14_Маринелли_0cm.spe",
        1024,
    ),
    (
        "detectors/Gamma-1S/reference_spectra/archive/"
        "Cs137_420-7-15_Маринелли_0cm.spe",
        1024,
    ),
    (
        "detectors/Gamma-1S/reference_spectra/archive/"
        "Поверка-2016/Маринелли/Th232_420-7-16_Маринелли_0cm.spe",
        1024,
    ),
]


def _resolve(rel: str) -> Path:
    """Resolve a fixture path against the project root."""
    return (ROOT / rel).resolve()


@pytest.mark.parametrize(
    "rel_path, expected_n",
    _FIXTURES,
    ids=[
        "Th232_420-7-17_Маринелли (reference_kits)",
        "Cs137_420-7-14_Маринелли (archive)",
        "Cs137_420-7-15_Маринелли (archive)",
        "Th232_420-7-16_Маринелли (Поверка-2016)",
    ],
)
def test_bug9_default_lsrm_read_preserves_all_channels(
    rel_path: str, expected_n: int
) -> None:
    """
    Default `read_lsrm_spe(path)` must return every channel from the
    binary block — no silent loss to the 3 MeV ceiling. Asserts the
    user-facing default path locks the lossless contract.
    """
    f = _resolve(rel_path)
    if not f.is_file():
        pytest.skip(f"fixture not present: {f}")

    # Default call — no kwargs. This is exactly what
    # `gamma.io.readers.read_spectrum(path)` and downstream consumers do.
    spec = read_lsrm_spe(str(f))

    assert spec.n_channels_raw == expected_n, (
        f"binary block decode mismatch: got {spec.n_channels_raw} ≠ "
        f"expected {expected_n} for {rel_path}"
    )
    assert int(len(spec.counts)) == expected_n, (
        f"BUG-9 regression: counts length {len(spec.counts)} ≠ raw "
        f"{expected_n} — reader is silently trimming channels again "
        f"(check apply_energy_ceiling default in lsrm_spe.read_lsrm_spe)"
    )
    assert spec.n_channels == expected_n, (
        f"spec.n_channels {spec.n_channels} ≠ raw {expected_n}"
    )


@pytest.mark.parametrize(
    "rel_path, expected_n",
    _FIXTURES,
    ids=[
        "Th232_420-7-17_Маринелли (reference_kits)",
        "Cs137_420-7-14_Маринелли (archive)",
        "Cs137_420-7-15_Маринелли (archive)",
        "Th232_420-7-16_Маринелли (Поверка-2016)",
    ],
)
def test_bug9_default_read_spectrum_preserves_all_channels(
    rel_path: str, expected_n: int
) -> None:
    """
    Same contract through the public `read_spectrum` dispatcher (what
    `scripts/run_skill.py` actually calls in Phase 1).
    """
    f = _resolve(rel_path)
    if not f.is_file():
        pytest.skip(f"fixture not present: {f}")

    spec = read_spectrum(str(f))

    assert spec.n_channels == expected_n, (
        f"BUG-9 regression in read_spectrum dispatch: "
        f"n_channels={spec.n_channels} ≠ {expected_n} for {rel_path}"
    )
    assert int(len(spec.counts)) == expected_n


def test_bug9_explicit_ceiling_still_works() -> None:
    """
    Opt-in to the legacy 3 MeV trim must still produce a shorter array
    on a known reproducer. This is the symmetric assertion to the
    default-keeps-all test above and pins the new opt-in API.
    """
    f = _resolve(
        "detectors/Gamma-1S/reference_spectra/reference_kits/"
        "Marinelli_1L/Th-232/Th232_420-7-17_Маринелли_0cm.spe"
    )
    if not f.is_file():
        pytest.skip(f"fixture not present: {f}")

    spec_default = read_lsrm_spe(str(f))
    spec_trim = read_lsrm_spe(str(f), apply_energy_ceiling=True)

    assert spec_default.n_channels == 1024
    assert spec_trim.n_channels < spec_default.n_channels, (
        f"explicit apply_energy_ceiling=True must still trim a "
        f"1024-ch fixture with E_max>3000 keV (got "
        f"{spec_trim.n_channels} vs default {spec_default.n_channels})"
    )
    assert spec_trim.energy_max_keV_kept is not None
    assert spec_trim.energy_max_keV_kept <= 3000.0
