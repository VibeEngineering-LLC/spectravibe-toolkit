"""
Tests for spectrum_qc_aggregator.py — F-QC-01 unified 6-criterion QC block.
Cite: spectrum_qc_methodology_v2_2026-06-03.md, KNOWN_AND_FIXED_ISSUES.md:1292, RAG-041.
"""
import pytest
import math
from unittest.mock import MagicMock
from gamma.reporting.spectrum_qc_aggregator import build_spectrum_qc


def test_energy_drift_pass():
    # seven_line_check.max_residual_keV=0.5, threshold=1.0
    # hand computation: 0.5 <= 1.0 -> passed=True
    slc = MagicMock(); slc.max_residual_keV = 0.5
    result = MagicMock(); result.seven_line_check = slc
    result.fwhm_at_661 = 8.5; result.efficiency_curve = MagicMock()
    result.bg_quality_check = None; result.background_subtraction = None
    qc = build_spectrum_qc(result)
    assert qc["energy_drift"]["passed"] is True
    assert qc["energy_drift"]["max_residual_keV"] == 0.5


def test_energy_drift_fail():
    # max_residual_keV=1.5, threshold=1.0 -> 1.5 > 1.0 -> passed=False
    slc = MagicMock(); slc.max_residual_keV = 1.5
    result = MagicMock(); result.seven_line_check = slc
    result.fwhm_at_661 = 8.5; result.efficiency_curve = MagicMock()
    result.bg_quality_check = None; result.background_subtraction = None
    qc = build_spectrum_qc(result)
    assert qc["energy_drift"]["passed"] is False
    assert qc["overall_passed"] is False


def test_energy_drift_unavailable():
    # seven_line_check=None -> passes by convention
    result = MagicMock(); result.seven_line_check = None
    result.fwhm_at_661 = 8.5; result.efficiency_curve = MagicMock()
    result.bg_quality_check = None; result.background_subtraction = None
    qc = build_spectrum_qc(result)
    assert qc["energy_drift"]["passed"] is True
    assert qc["energy_drift"]["available"] is False


def test_fwhm_stability_pass():
    # fwhm_at_661=8.0, reference=8.5
    # rel_dev = |8.0-8.5|/8.5 = 0.5/8.5 = 0.0588... < 0.15 -> pass
    result = MagicMock(); result.seven_line_check = None
    result.fwhm_at_661 = 8.0; result.fwhm_model_source = "file"
    result.efficiency_curve = MagicMock()
    result.bg_quality_check = None; result.background_subtraction = None
    qc = build_spectrum_qc(result, fwhm_reference_keV=8.5)
    assert qc["fwhm_stability"]["passed"] is True
    assert abs(qc["fwhm_stability"]["rel_deviation"] - 0.0588) < 0.001


def test_fwhm_stability_fail():
    # fwhm_at_661=12.0, reference=8.5
    # rel_dev = |12-8.5|/8.5 = 3.5/8.5 = 0.4117... > 0.15 -> fail
    result = MagicMock(); result.seven_line_check = None
    result.fwhm_at_661 = 12.0; result.fwhm_model_source = "fallback"
    result.efficiency_curve = MagicMock()
    result.bg_quality_check = None; result.background_subtraction = None
    qc = build_spectrum_qc(result, fwhm_reference_keV=8.5)
    assert qc["fwhm_stability"]["passed"] is False
    assert qc["overall_passed"] is False


def test_efficiency_qa_loaded():
    # efficiency_curve is not None -> passed=True
    result = MagicMock(); result.seven_line_check = None; result.fwhm_at_661 = None
    result.efficiency_curve = MagicMock(); result.efficiency_source = "file_loaded"
    result.bg_quality_check = None; result.background_subtraction = None
    qc = build_spectrum_qc(result)
    assert qc["efficiency_qa"]["passed"] is True
    assert qc["efficiency_qa"]["efficiency_loaded"] is True


def test_efficiency_qa_missing():
    # efficiency_curve=None -> passed=False
    result = MagicMock(); result.seven_line_check = None; result.fwhm_at_661 = None
    result.efficiency_curve = None; result.efficiency_source = ""
    result.bg_quality_check = None; result.background_subtraction = None
    qc = build_spectrum_qc(result)
    assert qc["efficiency_qa"]["passed"] is False
    assert qc["overall_passed"] is False


def test_bg_quality_check_passthrough():
    # bg_quality_check with 3 tested, 3 passed -> n_peaks_tested=3 in output
    bqc = {"n_peaks_tested": 3, "n_passed": 3, "n_failed": 0, "overall_passed": True, "peak_z_roi": [{"passed": True}, {"passed": True}, {"passed": True}]}
    result = MagicMock(); result.seven_line_check = None; result.fwhm_at_661 = None
    result.efficiency_curve = MagicMock(); result.bg_quality_check = bqc
    result.background_subtraction = None
    qc = build_spectrum_qc(result)
    assert qc["n_peaks_tested"] == 3
    assert qc["n_passed"] == 3
    assert len(qc["peak_z_roi"]) == 3


def test_bg_drift_unequal_live_times():
    """Criterion 4: rate-normalised z-test applied when sample_t != bg_t.

    Hand computation (Gilmore & Joss §5.5):
      sample_c=50000, sample_t=3600  -> R1=13.888... cps, Var(R1)=50000/3600^2=3.858e-3
      bg_c=900000,    bg_t=86400    -> R2=10.416... cps, Var(R2)=900000/86400^2=1.205e-4
      z = (R1-R2)/sqrt(Var(R1)+Var(R2)) = 3.472/sqrt(3.979e-3) = 3.472/0.06308 ≈ 55.0
      |z| >> 3.0 -> tier=reject, is_significant=True, passed=False
    """
    bg_sub = MagicMock()
    bg_sub.bg_live_time = 86400.0
    bg_sub.sample_sum_counts = 50000.0
    bg_sub.bg_sum_counts = 900000.0
    spec = MagicMock()
    spec.live_time = 3600.0
    result = MagicMock()
    result.seven_line_check = None
    result.fwhm_at_661 = None
    result.efficiency_curve = MagicMock()
    result.bg_quality_check = None
    result.background_subtraction = bg_sub
    result.spec = spec
    qc = build_spectrum_qc(result)
    assert qc["bg_drift"]["method"] == "rate_normalised"
    assert qc["bg_drift"]["available"] is True
    assert qc["bg_drift"]["is_significant"] is True
    assert qc["bg_drift"]["passed"] is False
    assert qc["overall_passed"] is False


def test_bg_drift_stable():
    """Criterion 4: rate-normalised z-test stable when counts match.

    Hand computation:
      sample_c=3600, sample_t=3600  -> R1=1.0 cps
      bg_c=86400,    bg_t=86400    -> R2=1.0 cps
      z = (1.0-1.0)/sqrt(3600/3600^2 + 86400/86400^2) = 0 -> tier=stable
    """
    bg_sub = MagicMock()
    bg_sub.bg_live_time = 86400.0
    bg_sub.sample_sum_counts = 3600.0
    bg_sub.bg_sum_counts = 86400.0
    spec = MagicMock()
    spec.live_time = 3600.0
    result = MagicMock()
    result.seven_line_check = None
    result.fwhm_at_661 = None
    result.efficiency_curve = MagicMock()
    result.bg_quality_check = None
    result.background_subtraction = bg_sub
    result.spec = spec
    qc = build_spectrum_qc(result)
    assert qc["bg_drift"]["method"] == "rate_normalised"
    assert qc["bg_drift"]["z"] == pytest.approx(0.0, abs=1e-9)
    assert qc["bg_drift"]["passed"] is True


def test_sensitivity_placeholder():
    """Criterion 6: sensitivity returns None (Phase 2 RC placeholder)."""
    result = MagicMock()
    result.seven_line_check = None
    result.fwhm_at_661 = None
    result.efficiency_curve = MagicMock()
    result.bg_quality_check = None
    result.background_subtraction = None
    qc = build_spectrum_qc(result)
    assert qc["sensitivity"] is None


def test_output_schema_keys():
    """All required schema keys present in output (acceptance criterion 2)."""
    result = MagicMock()
    result.seven_line_check = None
    result.fwhm_at_661 = None
    result.efficiency_curve = MagicMock()
    result.bg_quality_check = None
    result.background_subtraction = None
    qc = build_spectrum_qc(result)
    required_keys = {
        "n_peaks_tested", "n_passed", "n_failed", "overall_passed",
        "peak_z_roi", "energy_drift", "fwhm_stability", "efficiency_qa",
        "bg_drift", "sensitivity",
    }
    assert required_keys.issubset(set(qc.keys()))


def test_overall_passed_true_all_available():
    """All criteria available and passing -> overall_passed=True.

    Criteria status:
      energy_drift: max_residual=0.3 keV < 1.0 threshold -> pass
      fwhm_stability: 8.5 keV == reference 8.5 keV, rel_dev=0.0 < 0.15 -> pass
      efficiency_qa: curve loaded -> pass
      bg_drift: no background -> skip (pass)
      peak_z_roi: 0 peaks -> pass
      sensitivity: placeholder -> pass
    """
    slc = MagicMock(); slc.max_residual_keV = 0.3
    result = MagicMock()
    result.seven_line_check = slc
    result.fwhm_at_661 = 8.5
    result.fwhm_model_source = "file"
    result.efficiency_curve = MagicMock()
    result.efficiency_source = "file"
    result.bg_quality_check = None
    result.background_subtraction = None
    qc = build_spectrum_qc(result, fwhm_reference_keV=8.5)
    assert qc["overall_passed"] is True
    assert qc["energy_drift"]["passed"] is True
    assert qc["fwhm_stability"]["passed"] is True
    assert qc["efficiency_qa"]["passed"] is True
    assert qc["bg_drift"]["passed"] is True
