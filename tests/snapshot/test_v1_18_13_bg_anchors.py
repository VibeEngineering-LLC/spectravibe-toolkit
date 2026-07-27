# -*- coding: utf-8 -*-
"""v1.18.13 — F-298 / F-314 bg_lines_builder wiring в deconvolve.

Verifies:
- _f298_inject_bg_anchors добавляет bg-линии в components
- Дубликаты с real-нуклидами не вставляются
- Default OFF (enable_f96_bg_anchors=False) → back-compat
- Прокидывается через signature
"""
from __future__ import annotations
import os, sys
from pathlib import Path
import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


def test_F314_helper_imports():
    from gamma.peaks.deconvolve import _f298_inject_bg_anchors
    assert callable(_f298_inject_bg_anchors)


def test_F314_inject_adds_bg_lines_in_window():
    """F-96 bg линии в окне 1400-1500 keV → K-40 1461 должен добавиться."""
    from gamma.peaks.deconvolve import (
        _f298_inject_bg_anchors, MultipletComponent,
    )

    class _MockSpec:
        @staticmethod
        def energy_to_channel(E):
            return float(E)
    # Один real component (например, Eu-152 1408)
    real = [MultipletComponent(
        nuclide="Eu-152", line_E_keV=1408.0,
        library_I_pct=20.85, center_channel=1408.0, fwhm_channels=50.0,
    )]
    augmented = _f298_inject_bg_anchors(
        real, 1400.0, 1500.0, _MockSpec(),
        fwhm_at_channel=lambda ch: 50.0,
    )
    # K-40 1461 должен попасть в окно [1400, 1500]
    bg_nuclides = [c.nuclide for c in augmented if c.nuclide.startswith("bg:")]
    assert any("K-40" in n for n in bg_nuclides), (
        f"K-40 не добавлен как bg-anchor; augmented={[c.nuclide for c in augmented]}"
    )


def test_F314_no_duplicates_with_real_components():
    """Если real component уже на E_keV bg-линии → не дублировать."""
    from gamma.peaks.deconvolve import (
        _f298_inject_bg_anchors, MultipletComponent,
    )

    class _MockSpec:
        @staticmethod
        def energy_to_channel(E):
            return float(E)
    # Real component на 1461 keV (K-40)
    real = [MultipletComponent(
        nuclide="K-40", line_E_keV=1461.0,
        library_I_pct=10.7, center_channel=1461.0, fwhm_channels=50.0,
    )]
    augmented = _f298_inject_bg_anchors(
        real, 1400.0, 1500.0, _MockSpec(),
        fwhm_at_channel=lambda ch: 50.0,
    )
    k40_count = sum(1 for c in augmented if "K-40" in c.nuclide)
    assert k40_count == 1, (
        f"K-40 duplicated: {[c.nuclide for c in augmented]}"
    )


def test_F314_empty_window_returns_originals():
    """Если в окне нет bg-линий → augmented == original."""
    from gamma.peaks.deconvolve import (
        _f298_inject_bg_anchors, MultipletComponent,
    )

    class _MockSpec:
        @staticmethod
        def energy_to_channel(E):
            return float(E)
    real = [MultipletComponent(
        nuclide="X", line_E_keV=2500.0,
        library_I_pct=1.0, center_channel=2500.0, fwhm_channels=50.0,
    )]
    # Окно с очень узким диапазоном где нет канон bg-линий
    augmented = _f298_inject_bg_anchors(
        real, 2495.0, 2505.0, _MockSpec(),
        fwhm_at_channel=lambda ch: 50.0,
    )
    assert len(augmented) == 1


def test_F314_signature_accepts_enable_flag():
    """deconvolve_identified_multiplets должен принимать enable_f96_bg_anchors."""
    from gamma.peaks.deconvolve import deconvolve_identified_multiplets
    import inspect
    sig = inspect.signature(deconvolve_identified_multiplets)
    assert "enable_f96_bg_anchors" in sig.parameters
    assert sig.parameters["enable_f96_bg_anchors"].default is False
    assert "bg_anchor_min_intensity_pct" in sig.parameters


def test_F314_bg_marker_prefix():
    """All bg anchors должны иметь префикс 'bg:' для downstream-разделения."""
    from gamma.peaks.deconvolve import (
        _f298_inject_bg_anchors, MultipletComponent,
    )

    class _MockSpec:
        @staticmethod
        def energy_to_channel(E):
            return float(E)
    real = [MultipletComponent(
        nuclide="Cs-137", line_E_keV=662.0,
        library_I_pct=85.1, center_channel=662.0, fwhm_channels=46.0,
    )]
    augmented = _f298_inject_bg_anchors(
        real, 500.0, 800.0, _MockSpec(),
        fwhm_at_channel=lambda ch: 46.0,
    )
    for c in augmented:
        if c.nuclide != "Cs-137":
            assert c.nuclide.startswith("bg:"), (
                f"BG anchor missing 'bg:' prefix: {c.nuclide!r}"
            )
