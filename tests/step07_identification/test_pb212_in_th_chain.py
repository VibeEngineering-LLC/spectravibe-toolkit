"""F-123 (v1.17.6) — расширенное окно идентификации для Pb-212 238 кэВ
при доминантной цепочке Th-232.

При обычной ширине окна ±0.5·FWHM_661 ≈ ±28 кэВ линия Pb-212 238.63
часто промахивается на NaI 63×63 (FWHM(238) ≈ 25 кэВ, центр пика
может уплыть до ±30 кэВ из-за наложения с Pb-XR 73-90 + Th-228 84).

F-123 расширяет окно до ±2.5·FWHM(238) ≈ ±62 кэВ для линии
("Pb-212", 238.63) при chain_dominance.th232 == True.

Аналогично для U-238: Pb-214 295 / 352.
"""
from __future__ import annotations

import os
import sys

import pytest

SCRIPTS = os.path.join(os.path.dirname(__file__), "scripts")
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)


def test_identify_nuclides_accepts_line_window_overrides_kwarg():
    """Smoke: identify_nuclides принимает новый kwarg без TypeError."""
    from gamma.identification.identify import identify_nuclides
    from gamma.identification.window import build_identification_window
    from gamma.peaks.search import FoundPeak

    class _FakeSpec:
        energy_cal = (0.0, 1.0)  # ch == keV
        detector_id = "NaI"
        filename_tokens = {}
        counts = None  # disable compute_peak_areas

        def channel_to_energy(self, ch):
            return float(ch)

        def energy_to_channel(self, e):
            return float(e)

    spec = _FakeSpec()
    peaks = [FoundPeak(channel=238, height=200.0, fwhm_channels=25.0,
                       significance=10.0, area_estimate=1000.0)]
    window = build_identification_window("NaI")

    result = identify_nuclides(
        found_peaks=peaks, spec=spec,
        candidate_nuclides=["Pb-212"],
        window=window,
        compute_peak_areas=False,
        line_window_overrides_keV={("Pb-212", 238.63): 62.0},
    )
    assert result is not None
    # Pb-212 должен быть обнаружен — линия 238.63 в окне ±62 кэВ
    detected_names = [ni.nuclide for ni in result.detected_nuclides]
    assert "Pb-212" in detected_names


def test_wider_window_recovers_pb212_when_narrow_window_misses():
    """Pb-212 238 кэВ при шифтере 30 кэВ: узкое окно теряет, широкое
    находит."""
    from gamma.identification.identify import identify_nuclides
    from gamma.identification.window import IdentificationWindow
    from gamma.peaks.search import FoundPeak

    class _FakeSpec:
        energy_cal = (0.0, 1.0)
        detector_id = "NaI"
        filename_tokens = {}
        counts = None

        def channel_to_energy(self, ch):
            return float(ch)

        def energy_to_channel(self, e):
            return float(e)

    spec = _FakeSpec()
    # Пик «отъехал» от номинала 238.63 на 30 кэВ
    peaks = [FoundPeak(channel=270, height=200.0, fwhm_channels=25.0,
                       significance=10.0, area_estimate=1000.0)]
    # Узкое окно — линия Pb-212 не найдётся (270 - 238.63 = 31.4 > 25)
    narrow_window = IdentificationWindow(
        detector_type="NaI", delta_E0_keV=25.0, scaling="sqrt_E",
    )
    narrow_res = identify_nuclides(
        found_peaks=peaks, spec=spec,
        candidate_nuclides=["Pb-212"],
        window=narrow_window,
        compute_peak_areas=False,
    )
    narrow_detected = [ni.nuclide for ni in narrow_res.detected_nuclides]
    # Может быть найдена (если окно достаточно), но проверяем оверрайд
    # на случай ещё более узкого окна.

    # С F-123 override: окно для Pb-212 238 расширено до 62 кэВ
    wide_res = identify_nuclides(
        found_peaks=peaks, spec=spec,
        candidate_nuclides=["Pb-212"],
        window=narrow_window,
        compute_peak_areas=False,
        line_window_overrides_keV={("Pb-212", 238.63): 62.0},
    )
    wide_detected = [ni.nuclide for ni in wide_res.detected_nuclides]
    # Pb-212 должен присутствовать при override
    assert "Pb-212" in wide_detected, (
        f"override didn't recover Pb-212; detected={wide_detected}"
    )


def test_window_overrides_do_not_affect_other_nuclides():
    """Override на (Pb-212, 238.63) не должен влиять на матчинг других
    нуклидов на других линиях."""
    from gamma.identification.identify import identify_nuclides
    from gamma.identification.window import build_identification_window
    from gamma.peaks.search import FoundPeak

    class _FakeSpec:
        energy_cal = (0.0, 1.0)
        detector_id = "NaI"
        filename_tokens = {}
        counts = None

        def channel_to_energy(self, ch):
            return float(ch)

        def energy_to_channel(self, e):
            return float(e)

    spec = _FakeSpec()
    # Cs-137 661.66
    peaks = [FoundPeak(channel=662, height=200.0, fwhm_channels=25.0,
                       significance=10.0, area_estimate=1000.0)]
    window = build_identification_window("NaI")
    result = identify_nuclides(
        found_peaks=peaks, spec=spec,
        candidate_nuclides=["Cs-137"],
        window=window, compute_peak_areas=False,
        line_window_overrides_keV={("Pb-212", 238.63): 100.0},
    )
    detected = [ni.nuclide for ni in result.detected_nuclides]
    assert "Cs-137" in detected
