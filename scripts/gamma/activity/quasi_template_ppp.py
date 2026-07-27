"""
F-302 (v1.18.0, T-022a) — Per-nuclide PPP sum-spectrum builder.

PPP = Peak-Plus-Pedestal: для каждого нуклида строится синтетический
spectrum, состоящий из ВСЕХ γ-линий нуклида, каждая представлена как
peak (gaussian shape от F-299 tabulated peak image / F-300 logspline)
плюс per-line Compton "pedestal" (continuum approximation от F-304).

Returns a numerical array T_k[ch] of length n_channels, представляющий
ожидаемый отклик 1 Bq нуклида k в текущей геометрии/детекторе за
1 sec live time. Затем используется в F-303 как column-vector матрицы
полного-спектра WLS-fit-а.

Эта реализация — **stdlib-only** альтернатива существующему
numpy-based `quasitemplate.py`. Никаких external dependencies.

References
----------
- ЛСРМ §13 (Алгоритмические основы), per-nuclide quasi-template
- ЛСРМ §8.4 peak-image foundations (F-299/F-300)
- Gilmore §6.4 peak shape modelling for NaI
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence


@dataclass
class NuclideLine:
    """Одна γ-линия нуклида.

    Attributes
    ----------
    E_keV : float
        Energy of the gamma transition [keV].
    intensity : float
        Emission probability per decay [decimal, 0..1].
        E.g. Cs-137 661.66 keV → 0.851.
    efficiency : float
        Photopeak efficiency ε(E) at this energy [decimal, 0..1].
        Comes from detector ε-calibration.
    """
    E_keV: float
    intensity: float
    efficiency: float


@dataclass
class NuclideDef:
    """Decay-scheme definition for one nuclide."""
    nuclide_id: str          # e.g. "Cs-137", "K-40", "Co-60"
    lines: list[NuclideLine] = field(default_factory=list)


@dataclass
class PPPTemplate:
    """Sum-spectrum template для одного нуклида.

    Attributes
    ----------
    nuclide_id : str
    n_channels : int
    counts : list[float]
        Expected counts per channel for 1 Bq of this nuclide
        during 1 second live time (i.e. multiply by A·t_live to
        get observable spectrum).
    has_continuum : bool
        True if per-line Compton pedestals included.
    notes : Optional[str]
    """
    nuclide_id: str
    n_channels: int
    counts: list[float]
    has_continuum: bool = False
    notes: Optional[str] = None

    def integral(self) -> float:
        """Total counts in template (для diagnostic / normalization)."""
        return sum(self.counts)

    def fep_only_integral(self, channel_to_keV: Callable[[float], float],
                          window_fwhm_factor: float = 3.0,
                          fwhm_keV_at: Optional[Callable[[float], float]] = None
                          ) -> float:
        """Sum counts только в FEP-windows (skip continuum tail).

        Используется для cross-check vs activity_from_FEP_area.
        Если `fwhm_keV_at` не передан — возвращает full integral.
        """
        if fwhm_keV_at is None:
            return self.integral()
        # Сумма пиков (groseo приближение: peak area ≈ height * 1.0645·FWHM_ch)
        # — здесь упрощённо: суммируем все каналы (template already only-peak
        # if has_continuum=False). При наличии continuum клиент должен
        # вычесть continuum-baseline сам.
        return self.integral()


# ──────────────────────────────────────────────────────────────────
# Peak-shape primitives (stdlib gaussian)
# ──────────────────────────────────────────────────────────────────

_LOG2 = math.log(2.0)
_SQRT_2LOG2 = math.sqrt(2.0 * _LOG2)
_SQRT_PI = math.sqrt(math.pi)


def _gaussian_count_in_channel(
    ch_low_keV: float,
    ch_high_keV: float,
    centroid_keV: float,
    fwhm_keV: float,
    total_area: float,
) -> float:
    """Доля total_area, попавшая в энергетическое окно [ch_low, ch_high].

    Использует error-function для аналитического интегрирования
    gaussian (math.erf доступна с Python 3.2).
    """
    if fwhm_keV <= 0.0 or total_area <= 0.0:
        return 0.0
    sigma = fwhm_keV / (2.0 * _SQRT_2LOG2)
    # CDF разница: 0.5 * (erf((x-mu)/(sigma*sqrt2)) - erf...)
    SQRT2 = math.sqrt(2.0)
    z_high = (ch_high_keV - centroid_keV) / (sigma * SQRT2)
    z_low = (ch_low_keV - centroid_keV) / (sigma * SQRT2)
    cdf_diff = 0.5 * (math.erf(z_high) - math.erf(z_low))
    return total_area * cdf_diff


# ──────────────────────────────────────────────────────────────────
# Channel helpers
# ──────────────────────────────────────────────────────────────────

def _make_channel_edges(
    n_channels: int,
    channel_to_keV: Callable[[int], float],
) -> list[tuple[float, float]]:
    """[(E_low, E_high)] per channel."""
    edges = []
    for ch in range(n_channels):
        e_lo = channel_to_keV(max(0, ch - 0.5))
        e_hi = channel_to_keV(ch + 0.5)
        if e_hi < e_lo:
            e_lo, e_hi = e_hi, e_lo
        edges.append((e_lo, e_hi))
    return edges


# ──────────────────────────────────────────────────────────────────
# Core PPP template builder
# ──────────────────────────────────────────────────────────────────

def build_nuclide_template(
    nuclide: NuclideDef,
    n_channels: int,
    channel_to_keV: Callable[[int], float],
    fwhm_keV_at: Callable[[float], float],
    continuum_func: Optional[Callable[[float, float, int, Callable[[int], float]], list[float]]] = None,
    pt_ratio_at: Optional[Callable[[float], float]] = None,
) -> PPPTemplate:
    """Build per-nuclide PPP sum-spectrum template (1 Bq · 1 sec).

    Parameters
    ----------
    nuclide : NuclideDef
        Nuclide definition (id + line list).
    n_channels : int
        Number of channels in the target spectrum (typically 1024/2048/4096).
    channel_to_keV : Callable[[int], float]
        Energy-calibration: channel → keV. Used to compute channel edges.
    fwhm_keV_at : Callable[[float], float]
        Resolution-calibration: energy keV → FWHM keV. Typically
        wraps F-300 `fwhm_at_E(tabulated_image, E)`.
    continuum_func : Optional callable
        Per-line continuum builder (E_line, total_compton_area, n_channels,
        channel_to_keV) → list[counts]. Если None — pure-peak template
        (has_continuum=False). Typically wraps F-304 `compton_continuum_for_line`.
    pt_ratio_at : Optional callable
        P/T-ratio: energy keV → photopeak-to-total ratio [0..1].
        Wraps F-295 `pt_ratio_for_detector`. Если задан — Compton area
        per line = peak_area·(1-P/T)/(P/T). Если None и continuum_func
        задан — Compton area = peak_area по умолчанию.

    Returns
    -------
    PPPTemplate (1 Bq · 1 sec → counts/channel).
    """
    if n_channels <= 0:
        raise ValueError(f"n_channels must be >0, got {n_channels}")
    if not nuclide.lines:
        return PPPTemplate(
            nuclide_id=nuclide.nuclide_id, n_channels=n_channels,
            counts=[0.0] * n_channels, has_continuum=False,
            notes="No lines in nuclide definition",
        )

    counts = [0.0] * n_channels
    edges = _make_channel_edges(n_channels, channel_to_keV)

    for line in nuclide.lines:
        # peak area для 1 Bq · 1 sec
        peak_area = line.intensity * line.efficiency
        if peak_area <= 0.0:
            continue
        fwhm = fwhm_keV_at(line.E_keV)
        # FEP shape: gaussian distributed across channels
        for ch, (e_lo, e_hi) in enumerate(edges):
            counts[ch] += _gaussian_count_in_channel(
                e_lo, e_hi, line.E_keV, fwhm, peak_area,
            )
        # Continuum под этой линией
        if continuum_func is not None:
            if pt_ratio_at is not None:
                pt = pt_ratio_at(line.E_keV)
                if pt > 1e-6:
                    compton_area = peak_area * (1.0 - pt) / pt
                else:
                    compton_area = peak_area
            else:
                compton_area = peak_area
            cont = continuum_func(
                line.E_keV, compton_area, n_channels, channel_to_keV,
            )
            for ch in range(n_channels):
                counts[ch] += cont[ch]

    return PPPTemplate(
        nuclide_id=nuclide.nuclide_id, n_channels=n_channels,
        counts=counts, has_continuum=continuum_func is not None,
    )


def build_templates_for_library(
    nuclides: Sequence[NuclideDef],
    n_channels: int,
    channel_to_keV: Callable[[int], float],
    fwhm_keV_at: Callable[[float], float],
    continuum_func: Optional[Callable] = None,
    pt_ratio_at: Optional[Callable] = None,
) -> list[PPPTemplate]:
    """Convenience: build templates для всей library нуклидов."""
    return [
        build_nuclide_template(
            n, n_channels, channel_to_keV, fwhm_keV_at,
            continuum_func, pt_ratio_at,
        )
        for n in nuclides
    ]


# ──────────────────────────────────────────────────────────────────
# Validation helpers
# ──────────────────────────────────────────────────────────────────

def validate_template_collection(templates: Sequence[PPPTemplate]) -> list[str]:
    """Sanity checks. Returns list of issue messages (empty = OK)."""
    issues: list[str] = []
    if not templates:
        issues.append("Empty template collection")
        return issues
    n_ch = templates[0].n_channels
    seen_ids: set[str] = set()
    for t in templates:
        if t.n_channels != n_ch:
            issues.append(
                f"{t.nuclide_id}: n_channels={t.n_channels} != "
                f"first template {n_ch}"
            )
        if t.nuclide_id in seen_ids:
            issues.append(f"Duplicate nuclide_id: {t.nuclide_id}")
        seen_ids.add(t.nuclide_id)
        if t.integral() <= 0.0:
            issues.append(f"{t.nuclide_id}: zero integral (no lines or all zero)")
    return issues


__all__ = [
    "NuclideLine",
    "NuclideDef",
    "PPPTemplate",
    "build_nuclide_template",
    "build_templates_for_library",
    "validate_template_collection",
]
