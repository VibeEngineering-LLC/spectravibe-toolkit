"""Energy and FWHM calibration."""
from gamma.calibration.energy_fit import (
    EnergyFitResult, polynomial_energy_fit, MAX_POLYNOMIAL_DEGREE,
)
from gamma.calibration.stored_check import (
    StoredCheckResult, check_stored_calibration,
)
from gamma.calibration.bootstrap import (
    BootstrapResult, bootstrap_energy_calibration,
)
from gamma.calibration.subcalibration import (
    SubcalibrationResult, subcalibration_refit,
)
from gamma.calibration.fwhm_fit import (
    FwhmFitResult, fit_fwhm_hpge, fit_fwhm_scintillator,
)
from gamma.calibration.detector_type import (
    DetectorTypeResult, classify_detector,
)
from gamma.calibration.fwhm_provider import (
    make_fwhm_at_channel_provider,
)
from gamma.calibration.calibration_gate import (
    CalibrationGateResult, evaluate_calibration_gate,
)

__all__ = [
    "EnergyFitResult", "polynomial_energy_fit", "MAX_POLYNOMIAL_DEGREE",
    "StoredCheckResult", "check_stored_calibration",
    "BootstrapResult", "bootstrap_energy_calibration",
    "SubcalibrationResult", "subcalibration_refit",
    "FwhmFitResult", "fit_fwhm_hpge", "fit_fwhm_scintillator",
    "DetectorTypeResult", "classify_detector",
    "make_fwhm_at_channel_provider",
    "CalibrationGateResult", "evaluate_calibration_gate",
]
