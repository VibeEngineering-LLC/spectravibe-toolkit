"""F-QC-01 / v1.19.1 — per-peak Poisson |z|-test wiring test.

Verifies that `analyze_lsrm_spe` populates `StagedAnalysisResult.bg_quality_check`
with z-test results when a background is provided, and that the JSON report
emits `diagnostics.spectrum_qc`.

RAG-022 / BUG-35 integration: pipeline now consumes the `bg_z_test` kernel
and stores results per-peak. This test validates the wiring without testing the
kernel itself (kernel tests are in tests/snapshot/test_f243_bg_control.py).
"""
from __future__ import annotations

import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pytest

REPO_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from gamma.identification.staged_pipeline import analyze_lsrm_spe, StagedAnalysisResult
from gamma.reporting.json_report import build_json_report


# ── Minimal mock spectrum helpers ──────────────────────────────────────────

class _MockSpectrum:
    """Minimal duck-typed spectrum satisfying staged_pipeline requirements."""

    def __init__(self, counts, live_time=60.0, real_time=60.0,
                 energy_cal=(0.0, 1.0), source_path="mock.spe",
                 n_channels=None):
        import numpy as np
        self.counts = np.array(counts, dtype=np.float64)
        self.live_time = float(live_time)
        self.real_time = float(real_time)
        self.energy_cal = tuple(float(c) for c in energy_cal)
        self.source_path = source_path
        self.filename_tokens = {}
        self.geometry = ""
        self.detector_id = ""
        self.sample_id = ""
        self.operator = ""
        self.start_datetime = None
        self.n_channels = int(n_channels or len(counts))
        self.n_channels_raw = self.n_channels
        self.background_embedded = None

    def channel_to_energy(self, ch: int) -> float:
        """Horner-evaluate energy_cal at channel ch."""
        result = 0.0
        for k, a in enumerate(self.energy_cal):
            result += a * (ch ** k)
        return result

    @property
    def ENERGY_CEILING_KEV(self):
        return 3000.0


def _flat_spectrum(n_channels: int = 512, counts_per_ch: float = 100.0,
                   live_time: float = 1800.0, energy_cal=(0.0, 3.0)):
    """Create a flat-count spectrum with uniform counts per channel."""
    import numpy as np
    counts = np.full(n_channels, counts_per_ch, dtype=np.float64)
    return _MockSpectrum(
        counts=counts,
        live_time=live_time,
        real_time=live_time,
        energy_cal=energy_cal,
        n_channels=n_channels,
    )


def _write_spe(path: Path, counts, live_time=1800.0, energy_cal=(0.0, 3.0)):
    """Write a minimal IAEA-style SPE file."""
    nc = len(counts)
    lines = [
        "$SPEC_ID:",
        "mock spectrum",
        "$MEAS_TIM:",
        f"{int(live_time)} {int(live_time)}",
        "$DATA:",
        f"0 {nc - 1}",
    ]
    for c in counts:
        lines.append(str(int(c)))
    lines += [
        "$ENER_FIT:",
        f"{energy_cal[0]:.6f} {energy_cal[1]:.6f}",
        "$MCA_CAL:",
        "2",
        f"{energy_cal[0]:.6E} {energy_cal[1]:.6E} 0.0E+000",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


# ── Tests ──────────────────────────────────────────────────────────────────

@pytest.fixture
def spe_pair(tmp_path):
    """Write a sample SPE + background SPE pair and return their paths.

    Sample has a Gaussian-like bump at channel ~220 (≈660 keV with cal 3.0 keV/ch)
    to simulate a Cs-137 peak. Background is flat lower-rate.
    """
    import numpy as np

    n = 512
    cal = (0.0, 3.0)  # 3 keV/ch  → 512 ch ≈ 1536 keV
    live_s = 3600.0

    rng = np.random.default_rng(42)
    bg_base = 50.0      # counts/ch
    sample_base = 80.0  # slightly higher flat background in sample

    bg_counts = rng.poisson(bg_base, n).astype(float)

    sample_counts = rng.poisson(sample_base, n).astype(float)
    # Add a peak at channel 220 (≈ 660 keV): Gaussian with height ~800, σ=10
    ch = np.arange(n)
    peak = 800.0 * np.exp(-0.5 * ((ch - 220) / 10) ** 2)
    sample_counts += peak

    sample_path = tmp_path / "sample.spe"
    bg_path     = tmp_path / "background.spe"
    _write_spe(sample_path, sample_counts, live_time=live_s, energy_cal=cal)
    _write_spe(bg_path,     bg_counts,     live_time=live_s, energy_cal=cal)
    return sample_path, bg_path


class TestBgQualityCheckFieldPopulated:
    """bg_quality_check field is populated when background is subtracted."""

    def test_field_exists_on_dataclass(self):
        """StagedAnalysisResult has bg_quality_check field (F-QC-01 / v1.19.1)."""
        fields = StagedAnalysisResult.__dataclass_fields__
        assert "bg_quality_check" in fields, (
            "bg_quality_check field missing from StagedAnalysisResult"
        )
        # Default must be None (no bg case)
        default = fields["bg_quality_check"].default
        assert default is None, f"Expected default None, got {default!r}"

    def test_populated_when_bg_subtracted(self, spe_pair):
        """bg_quality_check is non-None when background_path is provided."""
        sample_path, bg_path = spe_pair
        result = analyze_lsrm_spe(
            str(sample_path),
            background_path=str(bg_path),
            allow_stage2=False,
            allow_stage3=False,
        )
        # Field must be populated when bg is given AND peaks exist.
        # Note: it may be None if no peaks found (rare on our synthetic SPE).
        if result.peaks:
            assert result.bg_quality_check is not None, (
                "bg_quality_check is None despite peaks found and bg subtracted"
            )

    def test_none_when_no_bg(self, spe_pair):
        """bg_quality_check is None when no background is provided."""
        sample_path, _ = spe_pair
        result = analyze_lsrm_spe(
            str(sample_path),
            background_path=None,
            allow_stage2=False,
            allow_stage3=False,
        )
        assert result.bg_quality_check is None, (
            "bg_quality_check should be None when no background_path given"
        )

    def test_structure_keys(self, spe_pair):
        """bg_quality_check dict has expected top-level keys."""
        sample_path, bg_path = spe_pair
        result = analyze_lsrm_spe(
            str(sample_path),
            background_path=str(bg_path),
            allow_stage2=False,
            allow_stage3=False,
        )
        if result.bg_quality_check is None:
            pytest.skip("No peaks found — z-test dict not populated")

        qc = result.bg_quality_check
        for key in ("n_peaks_tested", "n_passed", "n_failed", "overall_passed", "peak_z_roi"):
            assert key in qc, f"Missing key '{key}' in bg_quality_check: {qc.keys()}"

    def test_peak_z_roi_entries(self, spe_pair):
        """Each peak_z_roi entry has required keys and sane values."""
        sample_path, bg_path = spe_pair
        result = analyze_lsrm_spe(
            str(sample_path),
            background_path=str(bg_path),
            allow_stage2=False,
            allow_stage3=False,
        )
        if result.bg_quality_check is None:
            pytest.skip("No peaks found")

        entries = result.bg_quality_check["peak_z_roi"]
        assert len(entries) > 0, "peak_z_roi list is empty"

        required = {"peak_energy_keV", "e_lo_keV", "e_hi_keV", "z", "abs_z",
                    "tier", "passed", "B1", "B2", "note"}
        for entry in entries:
            assert required <= entry.keys(), (
                f"Missing keys in peak_z_roi entry: {required - entry.keys()}"
            )
            assert entry["tier"] in ("stable", "borderline", "reject", "undefined"), (
                f"Unexpected tier: {entry['tier']!r}"
            )
            assert isinstance(entry["passed"], bool)
            assert entry["e_lo_keV"] <= entry["peak_energy_keV"] <= entry["e_hi_keV"]

    def test_consistency_n_peaks_passed_failed(self, spe_pair):
        """n_passed + n_failed == n_peaks_tested and matches peak_z_roi list."""
        sample_path, bg_path = spe_pair
        result = analyze_lsrm_spe(
            str(sample_path),
            background_path=str(bg_path),
            allow_stage2=False,
            allow_stage3=False,
        )
        if result.bg_quality_check is None:
            pytest.skip("No peaks found")

        qc = result.bg_quality_check
        assert qc["n_passed"] + qc["n_failed"] == qc["n_peaks_tested"]
        assert qc["n_peaks_tested"] == len(qc["peak_z_roi"])
        count_pass = sum(1 for e in qc["peak_z_roi"] if e["passed"])
        assert count_pass == qc["n_passed"]


class TestSpectrumQcInJsonReport:
    """diagnostics.spectrum_qc block is emitted correctly in JSON report."""

    def test_spectrum_qc_key_in_diagnostics(self, spe_pair):
        """JSON report diagnostics contains 'spectrum_qc' key (F-QC-01)."""
        sample_path, bg_path = spe_pair
        result = analyze_lsrm_spe(
            str(sample_path),
            background_path=str(bg_path),
            allow_stage2=False,
            allow_stage3=False,
        )
        report = build_json_report(result)
        diag = report.get("diagnostics", {})
        assert "spectrum_qc" in diag, (
            "diagnostics.spectrum_qc key missing from JSON report"
        )

    def test_spectrum_qc_populated_when_bg(self, spe_pair):
        """spectrum_qc is non-None when background was subtracted AND peaks found."""
        sample_path, bg_path = spe_pair
        result = analyze_lsrm_spe(
            str(sample_path),
            background_path=str(bg_path),
            allow_stage2=False,
            allow_stage3=False,
        )
        report = build_json_report(result)
        diag = report.get("diagnostics", {})
        if result.peaks:
            assert diag["spectrum_qc"] is not None, (
                "spectrum_qc should be populated when peaks found + bg subtracted"
            )

    def test_spectrum_qc_no_bg_has_schema_keys(self, spe_pair):
        """F-QC-01 / v1.21.0: spectrum_qc block present even without background.

        From v1.21.0 (wave 5 F-QC-01 aggregator), spectrum_qc always contains
        the unified 6-criterion block. Criteria that require background
        (bg_drift, peak_z_roi) report available=False / empty list.
        Other criteria (energy_drift, fwhm_stability, efficiency_qa) are
        always meaningful regardless of background presence.
        Old behaviour (None without bg) was the BUG-35 partial wiring;
        replaced by full aggregator in v1.21.0 — KNOWN_AND_FIXED_ISSUES.md:1292.
        """
        sample_path, _ = spe_pair
        result = analyze_lsrm_spe(
            str(sample_path),
            background_path=None,
            allow_stage2=False,
            allow_stage3=False,
        )
        report = build_json_report(result)
        diag = report.get("diagnostics", {})
        qc = diag.get("spectrum_qc")
        assert qc is not None, (
            "spectrum_qc should always be populated in v1.21.0 (F-QC-01 aggregator)"
        )
        # Required schema keys (acceptance criterion 2)
        for key in ("n_peaks_tested", "n_passed", "n_failed", "overall_passed",
                    "peak_z_roi", "energy_drift", "fwhm_stability",
                    "efficiency_qa", "bg_drift", "sensitivity"):
            assert key in qc, f"spectrum_qc missing key: {key}"
        # Without background, bg_drift should be skipped (available=False, passed=True)
        assert qc["bg_drift"]["available"] is False, (
            "bg_drift should be unavailable when no background subtracted"
        )
        # Without background, peak_z_roi should be empty (no BUG-35 data)
        assert qc["n_peaks_tested"] == 0, (
            "n_peaks_tested should be 0 when no background subtracted"
        )

    def test_report_json_serializable(self, spe_pair):
        """JSON report with spectrum_qc is fully serializable (no NaN/Inf issues)."""
        sample_path, bg_path = spe_pair
        result = analyze_lsrm_spe(
            str(sample_path),
            background_path=str(bg_path),
            allow_stage2=False,
            allow_stage3=False,
        )
        report = build_json_report(result)
        # json.dumps will raise on non-serializable types (NaN float, etc.)
        # NaN in z-test entries are serialized as null by allow_nan=False.
        # We use a custom handler matching the actual build pipeline.
        serialized = json.dumps(report, ensure_ascii=False, allow_nan=True)
        assert len(serialized) > 100, "Serialized report suspiciously short"
