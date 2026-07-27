# -*- coding: utf-8 -*-
"""Коэффициенты и парсер спектра для метода ЛСРМ (метод (3)).

Источник формул: "Мощность дозы. Методика расчёта из спектра гамма-излучения"
(ВНИИФТРИ/ЛСРМ/Аспект, Менделеево-Дубна). Извлечено VLM-OCR (baidu/Unlimited-OCR),
сверено с оригиналом (см. README.md, провенанс methodology_ocr/page_04).

Табл.12.2: энергия[МэВ], μen/ρ[м²/кг] в ВОЗДУХЕ, f(10) экспозиц.→эквив.H*(10).
Диапазон методики 50–3000 кэВ (вне — канал в дозу не вносится).
"""
from __future__ import annotations
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

E_LO, E_HI = 50.0, 3000.0  # диапазон методики, кэВ

# спектр первого прогона (подвал радон) — дефолт для воспроизводимости
DEFAULT_SPECTRUM = Path(
    r"<GDRIVE>\Дозиметрия\Спектры\Atom GS5050 PRO MAX\Фон подвал радон.xml"
)

# Табл.12.2: E[МэВ], μen/ρ[м²/кг] воздух, f(10)
_T122 = np.array([
    [0.040, 0.006894, 1.29], [0.050, 0.004031, 1.46], [0.060, 0.003004, 1.52],
    [0.080, 0.002393, 1.51], [0.100, 0.002318, 1.44], [0.150, 0.002494, 1.31],
    [0.200, 0.002672, 1.22], [0.300, 0.002872, 1.15], [0.400, 0.002949, 1.10],
    [0.500, 0.002966, 1.07], [0.600, 0.002953, 1.04], [0.800, 0.002882, 1.02],
    [1.000, 0.002787, 1.01], [1.500, 0.002545, 0.99], [2.000, 0.002342, 0.99],
    [3.000, 0.002054, 0.98],
])
E_MEV, MUEN_AIR, F10 = _T122[:, 0], _T122[:, 1], _T122[:, 2]


def muen_air_interp(E_keV):
    """μen/ρ воздуха [м²/кг] линейной интерполяцией по энергии (кэВ→МэВ)."""
    return np.interp(np.asarray(E_keV, float) / 1000.0, E_MEV, MUEN_AIR)


def f10_interp(E_keV):
    """f(10) экспозиц.→H*(10) линейной интерполяцией по энергии."""
    return np.interp(np.asarray(E_keV, float) / 1000.0, E_MEV, F10)


def parse_spectrum(path) -> dict:
    """AtomSpectra/BecqMoni XML → dict(counts, live, valid, ch, energy_keV, rate_cps).

    Минимальный самодостаточный парсер (полином калибровки .//EnergySpectrum).
    Для других форматов (.spe/.n42/.txt/RadiaCode) — см. README «Точки расширения»:
    gamma.io.format_registry.
    """
    path = Path(path)
    root = ET.parse(path).getroot()
    es = root.find(".//EnergySpectrum")
    coeffs = [float(c.text) for c in es.find("EnergyCalibration/Coefficients")]
    counts = np.array([float(dp.text) for dp in es.find("Spectrum")], dtype=float)
    live = float(es.findtext("LiveTime"))
    valid = float(es.findtext("ValidPulseCount"))
    ch = np.arange(len(counts), dtype=float)
    energy = np.polynomial.polynomial.polyval(ch, coeffs)
    return {"path": path, "counts": counts, "live": live, "valid": valid,
            "ch": ch, "energy_keV": energy, "rate_cps": valid / live}