"""
PNG plot rendering for the Step 11 report (F-86a / v1.15.0).

Two artefacts per SKILL.md "Display conventions":

  * build_spectrum_plot(result, output_path)
        Full spectrum overlay — counts-per-second, log Y, no smoothing
        on data points, title with metadata, primary FEPs labelled
        with nuclide+E, secondary peaks in distinct color, optional
        embedded-background overlay at reduced opacity.

  * build_multiplet_plots(result, output_dir) -> List[Path]
        One PNG per resolved multiplet cluster. Raw counts in the ROI,
        component centroids as vertical lines labelled with
        nuclide+E, chi²/dof + convergence annotation in the corner.

Uses the non-interactive ``Agg`` matplotlib backend, so no display
server is required (CI / notebooks / CLI all work).

The module degrades gracefully:
  * If matplotlib is not installed: both functions raise
    ImportError with a clear message — the caller (build_report)
    catches and falls back to the v1.14.0 placeholder.
  * If the result lacks counts/peaks/clusters: returns ``None`` /
    empty list rather than crashing.
"""
from __future__ import annotations

import os
import math
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import numpy as np

from gamma.reporting._kind_ru import FEATURE_KIND_RU_SHORT as _PLOT_KIND_RU_SHARED
# F-452 / v1.33.0 — FwhmModel polymorphic API.
from gamma.identification.staged_pipeline import (
    fwhm_keV_at_energy as _fwhm_keV_at_energy,
)


# ──────────────────────────────────────────────────────────────────
# matplotlib bootstrap (Agg backend, lazy import)
# ──────────────────────────────────────────────────────────────────

def _import_matplotlib():
    """Lazy matplotlib import with Agg backend. Returns plt module."""
    try:
        import matplotlib
        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
        return plt
    except ImportError as e:
        raise ImportError(
            "matplotlib is required for plot generation. "
            "Install with `pip install matplotlib` (>=3.7)."
        ) from e


# ──────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────

# Marker colors — colour-blind friendly per SKILL.md
_COLOR_DATA = "#1f3a5f"          # navy: raw data
_COLOR_BG = "#a8a8a8"            # mid-gray: background overlay
_COLOR_PRIMARY = "#cc3322"       # vermillion: primary FEPs
_COLOR_SECONDARY = "#118866"     # teal: secondary peaks (XRF, escape, sum)
_COLOR_FIT = "#dd9933"           # amber: deconvolution fit envelope
_COLOR_COMPONENT = "#557755"     # olive: individual deconvolution components


def _gaussian(x_ch: np.ndarray, area: float, mu_ch: float,
              fwhm_ch: float) -> np.ndarray:
    """Gaussian normalised so ``area = ∫ g dx``."""
    if fwhm_ch <= 0:
        return np.zeros_like(x_ch)
    sigma = fwhm_ch / 2.3548200450309493
    norm = area / (sigma * math.sqrt(2.0 * math.pi))
    return norm * np.exp(-0.5 * ((x_ch - mu_ch) / sigma) ** 2)


def _safe_filename(text: str, *, max_len: int = 80) -> str:
    """Filesystem-safe slug. Keeps non-ASCII letters (UTF-8 fine)."""
    bad = '<>:"/\\|?*\n\r\t'
    out = "".join("_" if c in bad else c for c in (text or "report"))
    out = out.strip(" ._")
    return (out or "report")[:max_len]


def _energy_axis_for_counts(spec, n: int) -> np.ndarray:
    """Build a per-channel energy array, robust to None calibration."""
    chs = np.arange(n, dtype=np.float64)
    Es = np.empty(n, dtype=np.float64)
    for i in range(n):
        e = spec.channel_to_energy(int(chs[i]))
        Es[i] = float(e) if e is not None else float(i)
    return Es


def _fwhm_keV_at(result, E: float) -> float:
    """Evaluate stored FWHM model at energy E.

    F-452: polymorphic — `result.fwhm_model` теперь FwhmModel (callable)
    или legacy 3-tuple; единая точка вычисления — `fwhm_keV_at_energy`.
    """
    return _fwhm_keV_at_energy(result.fwhm_model, float(E))


def _channel_for_energy(spec, target_E: float, energies: np.ndarray) -> Optional[int]:
    """Find the channel whose calibrated energy is closest to target_E."""
    if energies.size == 0:
        return None
    idx = int(np.argmin(np.abs(energies - target_E)))
    return idx


# ──────────────────────────────────────────────────────────────────
# Public: spectrum plot
# ──────────────────────────────────────────────────────────────────

def build_spectrum_plot(
    result,
    output_path,
    *,
    dpi: int = 120,
    include_background: bool = True,
    max_primary_labels: int = 25,
    max_secondary_labels: int = 15,
) -> Optional[str]:
    """Render the cps-log spectrum overlay PNG. Returns the written path.

    Parameters
    ----------
    result : StagedAnalysisResult
    output_path : str | Path
        Destination ``.png`` file (parent directories are created).
    dpi : int, default 120
    include_background : bool, default True
        If the spec has an embedded background spectrum, overlay it at
        ``alpha=0.35``.
    max_primary_labels : int
        Cap on annotated primary FEP labels to keep the plot legible.
    max_secondary_labels : int
        Same for secondary peaks.

    Returns
    -------
    str | None
        Absolute path of the written PNG, or ``None`` when the spec has
        no counts / no calibration (plot would be meaningless).
    """
    plt = _import_matplotlib()

    spec = result.spec
    counts = np.asarray(spec.counts, dtype=np.float64) if spec.counts is not None else None
    if counts is None or counts.size == 0:
        return None
    live = float(getattr(spec, "live_time", 0.0) or 0.0)
    if live <= 0:
        # cps would divide by zero; fall back to raw counts but warn in title
        rate = counts
        # F-397.4 / v1.18.28.1 (Agent B) — RU labels.
        y_label = "Счёт (нет live time)"
    else:
        rate = counts / live
        y_label = "Скорость счёта, имп/с"

    energies = _energy_axis_for_counts(spec, len(counts))

    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(12, 5.5), dpi=dpi)

    # Raw spectrum — no smoothing
    # F-397.4 / v1.18.28.1 (Agent B) — legend label «Образец», не «Sample».
    ax.plot(energies, rate, drawstyle="steps-mid",
            linewidth=0.6, color=_COLOR_DATA, label="Образец")

    # Background overlay (alpha)
    bg = getattr(spec, "background_embedded", None)
    bg_plotted = False
    if include_background and bg is not None:
        bg_counts = np.asarray(bg.counts, dtype=np.float64) if bg.counts is not None else None
        bg_live = float(getattr(bg, "live_time", 0.0) or 0.0)
        if bg_counts is not None and bg_counts.size > 0 and bg_live > 0:
            bg_E = _energy_axis_for_counts(bg, len(bg_counts))
            bg_rate = bg_counts / bg_live
            ax.plot(bg_E, bg_rate, drawstyle="steps-mid",
                    linewidth=0.5, color=_COLOR_BG, alpha=0.45,
                    label="Background")
            bg_plotted = True

    # Primary FEPs (final_detected nuclides — characteristic line)
    primary_marked: List[Tuple[float, float, str]] = []
    for det in (result.final_detected or [])[:max_primary_labels]:
        # Try to find a single representative line for marker placement
        E_line = getattr(det, "characteristic_line_keV", None)
        if E_line is None:
            lines = getattr(det, "matched_lines_keV", None) or ()
            if lines:
                E_line = float(lines[0])
        if E_line is None or not math.isfinite(float(E_line)):
            continue
        ch = _channel_for_energy(spec, float(E_line), energies)
        if ch is None:
            continue
        y_val = rate[ch] if 0 <= ch < len(rate) else None
        if y_val is None or y_val <= 0:
            continue
        nuclide = getattr(det, "nuclide", "?")
        primary_marked.append((float(E_line), float(y_val), nuclide))

    if primary_marked:
        xs = [p[0] for p in primary_marked]
        ys = [p[1] for p in primary_marked]
        ax.scatter(xs, ys, marker="v", color=_COLOR_PRIMARY, s=45,
                   zorder=5, label="Основные ФЭП")
        for E_line, y_val, nuclide in primary_marked:
            ax.annotate(
                f"{nuclide}\n{E_line:.0f}",
                xy=(E_line, y_val),
                xytext=(0, 12), textcoords="offset points",
                ha="center", fontsize=6.5, color=_COLOR_PRIMARY,
                rotation=0,
            )

    # Secondary peaks (residual classifications — XRF, escape, sum, etc.)
    # TD-4 / v1.18.30 — feature_kind → RU из shared _kind_ru.FEATURE_KIND_RU_SHORT.
    # Расширяйте словарь в _kind_ru.py, а не здесь.
    secondary_marked: List[Tuple[float, float, str]] = []
    n_sec = 0
    for rc in (result.residual_classifications or []):
        if n_sec >= max_secondary_labels:
            break
        label = getattr(rc, "label", "")
        if label in ("", "true_unmatched"):
            continue
        # F-397.4 — translate если есть mapping; иначе оставляем raw label
        # (например для редких feature_kinds которые в map ещё не попали).
        label = _PLOT_KIND_RU_SHARED.get(label, label)
        E_peak = float(getattr(rc, "peak_E_keV", 0.0) or 0.0)
        if E_peak <= 0:
            continue
        ch = _channel_for_energy(spec, E_peak, energies)
        if ch is None:
            continue
        y_val = rate[ch] if 0 <= ch < len(rate) else None
        if y_val is None or y_val <= 0:
            continue
        secondary_marked.append((E_peak, float(y_val), str(label)))
        n_sec += 1

    if secondary_marked:
        xs = [p[0] for p in secondary_marked]
        ys = [p[1] for p in secondary_marked]
        ax.scatter(xs, ys, marker="^", color=_COLOR_SECONDARY, s=30,
                   zorder=4, label="Вторичные")
        for E_peak, y_val, label in secondary_marked:
            ax.annotate(
                f"{label}\n{E_peak:.0f}",
                xy=(E_peak, y_val),
                xytext=(0, -14), textcoords="offset points",
                ha="center", fontsize=5.8, color=_COLOR_SECONDARY,
            )

    # Axes
    ax.set_yscale("log")
    # F-397.4 / v1.18.28.1 (Agent B) — RU axis labels для консистентности
    # с RU report (избегаем смешения keV/кэВ в одном bundle).
    ax.set_xlabel("Энергия, кэВ")
    ax.set_ylabel(y_label)
    ax.grid(True, which="both", linestyle=":", linewidth=0.4, alpha=0.5)

    # Title with metadata
    sp = getattr(spec, "source_path", "") or ""
    leaf = os.path.basename(sp) if sp else "spectrum"
    sd = getattr(spec, "start_datetime", None)
    sd_str = sd.strftime("%Y-%m-%d %H:%M") if sd is not None else "?"
    real = float(getattr(spec, "real_time", 0.0) or 0.0)
    sample_id = getattr(spec, "sample_id", "") or ""
    detector = (result.detector_canonical
                or getattr(spec, "detector_id", None)
                or result.detector_type
                or "?")
    title_lines = [
        leaf,
        f"{sd_str} · live {live:.0f} s / real {real:.0f} s · detector {detector}",
    ]
    if sample_id:
        title_lines.append(f"sample: {sample_id}")
    ax.set_title("\n".join(title_lines), fontsize=9, loc="left")

    ax.legend(loc="upper right", fontsize=7, framealpha=0.85)

    fig.tight_layout()
    fig.savefig(str(out_path), dpi=dpi)
    plt.close(fig)
    return str(out_path)


# ──────────────────────────────────────────────────────────────────
# Public: multiplet cluster plots
# ──────────────────────────────────────────────────────────────────

def build_multiplet_plots(
    result,
    output_dir,
    *,
    dpi: int = 120,
    filename_prefix: str = "multiplet",
) -> List[str]:
    """Render one PNG per resolved multiplet cluster. Returns paths.

    Each plot shows the raw counts in the ROI, the deconvolution fit
    envelope (Gaussian sum + linear continuum) and the individual
    component Gaussians, labelled with nuclide + library energy.
    χ²/dof and convergence are annotated.

    Returns ``[]`` when ``result.deconvolution_results`` is empty / None.
    """
    decons = result.deconvolution_results or []
    if not decons:
        return []

    plt = _import_matplotlib()

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    spec = result.spec
    counts = np.asarray(spec.counts, dtype=np.float64) if spec.counts is not None else None
    if counts is None or counts.size == 0:
        return []
    live = float(getattr(spec, "live_time", 0.0) or 0.0) or 1.0

    paths: List[str] = []
    for i, d in enumerate(decons, start=1):
        roi_lo = int(d.roi_low_ch)
        roi_hi = int(d.roi_high_ch)
        if not (0 <= roi_lo < roi_hi <= len(counts)):
            continue

        # F-134 / v1.17.7 — если CoupledFitResult оставил готовые
        # overlay-массивы (E, data, continuum, total, per-component),
        # используем их напрямую. Это закрывает рассинхрон между точным
        # fit'ом (closure ≈ 0) и канальной реконструкцией PNG-рендером
        # (которая для формы Гаусс+tail+step считалась некорректно).
        if (d.overlay_E_keV and d.overlay_data and
                d.overlay_continuum and d.overlay_total):
            E_at = np.asarray(d.overlay_E_keV, dtype=np.float64)
            roi_counts = np.asarray(d.overlay_data, dtype=np.float64)
            cont = np.asarray(d.overlay_continuum, dtype=np.float64)
            envelope = np.asarray(d.overlay_total, dtype=np.float64)
            # per-component overlays: cont + только эта компонента
            comp_curves: List[Tuple[np.ndarray, str, float]] = []
            overlays = d.overlay_components or ()
            for k, (comp, area) in enumerate(zip(d.components, d.areas)):
                E_line = float(comp.line_E_keV)
                if k < len(overlays):
                    g_plus_cont = np.asarray(overlays[k], dtype=np.float64)
                else:
                    g_plus_cont = cont
                # ((g + cont) - cont) даёт чистый g для совместимости с
                # downstream аннотацией пика
                g = g_plus_cont - cont
                comp_curves.append((g, getattr(comp, "nuclide", "?"), E_line))
        else:
            # Legacy fallback: канальная реконструкция (только для свободных
            # NNLS-fit'ов из apply_multiplet_deconvolution; F-117 coupled
            # сюда не попадает после F-134).
            roi_counts = counts[roi_lo:roi_hi]
            roi_chs = np.arange(roi_lo, roi_hi, dtype=np.float64)
            E_at = np.asarray(
                [spec.channel_to_energy(int(c)) or 0.0 for c in roi_chs],
                dtype=np.float64,
            )
            cp = list(d.continuum_params or ())
            cont = np.zeros_like(roi_chs)
            if len(cp) >= 1:
                cont += cp[0]
            if len(cp) >= 2:
                cont += cp[1] * roi_chs
            comp_curves = []
            for comp, area in zip(d.components, d.areas):
                E_line = float(comp.line_E_keV)
                mu_ch_idx = _channel_for_energy(spec, E_line, E_at)
                mu_ch = (float(roi_chs[mu_ch_idx])
                         if mu_ch_idx is not None
                         else (roi_lo + roi_hi) / 2.0)
                fwhm_keV = _fwhm_keV_at(result, E_line)
                if E_at.size >= 2:
                    dE_per_ch = float(np.mean(np.diff(E_at)))
                    fwhm_ch = (abs(fwhm_keV / dE_per_ch)
                               if dE_per_ch != 0 else 1.0)
                else:
                    fwhm_ch = 1.0
                g = _gaussian(roi_chs, float(area), mu_ch, fwhm_ch)
                comp_curves.append((g, getattr(comp, "nuclide", "?"), E_line))
            envelope = cont.copy()
            for g, _, _ in comp_curves:
                envelope = envelope + g

        # Plot
        fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=dpi)
        ax.plot(E_at, roi_counts, drawstyle="steps-mid",
                linewidth=0.8, color=_COLOR_DATA, label="ROI counts")
        ax.plot(E_at, envelope, linewidth=1.2, color=_COLOR_FIT,
                label="Fit envelope")
        ax.plot(E_at, cont, linewidth=0.7, linestyle="--",
                color=_COLOR_BG, label="Continuum")
        for g, nuclide, E_line in comp_curves:
            ax.plot(E_at, g + cont, linewidth=0.7, alpha=0.7,
                    color=_COLOR_COMPONENT)
            # Annotate at peak
            yv = float(np.max(g + cont))
            ax.annotate(
                f"{nuclide}\n{E_line:.1f}",
                xy=(E_line, yv),
                xytext=(0, 8), textcoords="offset points",
                ha="center", fontsize=6.5, color=_COLOR_COMPONENT,
            )

        ax.set_yscale("linear")
        # F-397.4 / v1.18.28.1 (Agent B) — RU axis labels.
        ax.set_xlabel("Энергия, кэВ")
        ax.set_ylabel("Счёт в окне")
        ax.grid(True, which="both", linestyle=":", linewidth=0.3, alpha=0.5)

        chi = d.chi2_per_dof
        conv = "yes" if d.converged else "no"
        E_lo = float(E_at[0])
        E_hi = float(E_at[-1])
        title = (
            f"Cluster {i}: {E_lo:.0f}–{E_hi:.0f} keV  ·  "
            f"χ²/dof={chi:.2f}  ·  converged: {conv}"
        )
        ax.set_title(title, fontsize=9, loc="left")
        ax.legend(loc="upper right", fontsize=7, framealpha=0.85)

        fig.tight_layout()
        name = f"{filename_prefix}_{i:02d}_{int(E_lo)}-{int(E_hi)}keV.png"
        out_p = out_dir / _safe_filename(name)
        fig.savefig(str(out_p), dpi=dpi)
        plt.close(fig)
        paths.append(str(out_p))

    return paths


# ──────────────────────────────────────────────────────────────────
# Public: bundle (called by build_report)
# ──────────────────────────────────────────────────────────────────

def build_all_plots(
    result,
    output_dir,
    *,
    spectrum_filename: str = "spectrum.png",
    multiplet_subdir: str = "multiplets",
    dpi: int = 120,
) -> dict:
    """Generate both the spectrum overlay and all multiplet plots.

    Returns a dict with keys ``spectrum`` (str | None) and
    ``multiplets`` (List[str]).
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    spec_path = build_spectrum_plot(
        result, out_dir / spectrum_filename, dpi=dpi,
    )
    mp_dir = out_dir / multiplet_subdir
    mp_paths = build_multiplet_plots(result, mp_dir, dpi=dpi)
    return {"spectrum": spec_path, "multiplets": mp_paths}


__all__ = [
    "build_spectrum_plot",
    "build_multiplet_plots",
    "build_all_plots",
]
