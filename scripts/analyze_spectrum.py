#!/usr/bin/env python3
from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""
Reference skeleton for gamma-spectrum analysis following the methodology
described in SKILL.md.

This is a STARTING POINT, not a complete production tool. It demonstrates:
  - Multi-format file reading (.spe, .chn, .txt/.csv); .n42 and .spm stubbed
  - Filename token parsing for metadata-as-prior
  - Mariscotti second-derivative peak search in channel space
  - Bootstrap energy calibration via anchor patterns (60Co doublet, 40K, etc.)
  - FWHM(E) calibration (HPGe quadratic-in-E form)
  - Library-directed search
  - Identification with characteristic-line check and intensity ratios
  - Targeted multiplet deconvolution with fixed positions and library ratios
  - Currie/ISO 11929 MDA
  - Lsrm confidence index CI and Dose Contribution DC

Dependencies: numpy, scipy, lmfit
Optional: matplotlib for plotting

Extend the format readers, peak-shape models, and nuclear-data library for
the actual detector and use case at hand.
"""

import re
import math
import json
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit


# ============================================================================
# Data classes
# ============================================================================

@dataclass
class Spectrum:
    counts: np.ndarray             # 1D array, counts per channel
    live_time: float               # seconds
    real_time: float               # seconds
    start_datetime: Optional[datetime] = None
    sample_id: str = ""
    operator: str = ""
    geometry: str = ""
    detector_id: str = ""
    comments: str = ""
    energy_cal: Optional[tuple] = None      # (a0, a1, a2, ...) such that E = sum(a_i * N^i)
    fwhm_cal: Optional[tuple] = None        # detector-dependent meaning
    filename_tokens: dict = field(default_factory=dict)
    is_background: bool = False
    source_path: str = ""

    @property
    def n_channels(self) -> int:
        return len(self.counts)

    @property
    def dead_time_fraction(self) -> float:
        if self.real_time <= 0:
            return 0.0
        return 1.0 - self.live_time / self.real_time

    def channel_to_energy(self, ch):
        if self.energy_cal is None:
            return None
        return sum(a * ch ** i for i, a in enumerate(self.energy_cal))

    def energy_to_channel(self, e):
        """Inverse: solve E(N) = e numerically."""
        if self.energy_cal is None:
            return None
        # For linear case: N = (E - a0) / a1
        if len(self.energy_cal) == 2:
            a0, a1 = self.energy_cal
            return (e - a0) / a1
        # Otherwise Newton iteration
        N = (e - self.energy_cal[0]) / self.energy_cal[1]
        for _ in range(10):
            E_N = self.channel_to_energy(N)
            dE_dN = sum(i * a * N ** (i - 1) for i, a in enumerate(self.energy_cal) if i > 0)
            N = N - (E_N - e) / dE_dN
        return N

    def fwhm_at(self, E):
        """Return FWHM (keV) at energy E, using fwhm_cal."""
        if self.fwhm_cal is None:
            return None
        # HPGe model: FWHM^2 = a + b*E + c*E^2
        if len(self.fwhm_cal) >= 3:
            a, b, c = self.fwhm_cal[:3]
            val = a + b * E + c * E * E
            return math.sqrt(max(val, 0))
        # Scintillator model: FWHM = k * sqrt(E)
        k = self.fwhm_cal[0]
        return k * math.sqrt(E)


@dataclass
class Peak:
    channel: float
    energy: Optional[float] = None
    fwhm: Optional[float] = None
    area: float = 0.0
    area_sigma: float = 0.0
    significance: float = 0.0
    in_multiplet: bool = False
    assigned_nuclide: str = ""
    type: str = "FEP"   # FEP or one of secondary types


@dataclass
class Nuclide:
    name: str
    lines: list                    # list of (E_keV, I_gamma_percent)
    half_life: float = 0.0         # seconds; 0 = stable
    is_cascade: bool = False       # has true coincidence summing


# ============================================================================
# Filename parsing
# ============================================================================

NUCLIDE_TOKENS = [
    "Cs137", "Co60", "Eu152", "Am241", "Ra226", "Ba133", "Na22",
    "K40", "Th232", "U235", "U238", "Mn54", "Zn65", "Be7", "Mo99",
    "Tc99m", "I131", "F18", "Na24", "Co57",
]
GEOMETRY_TOKENS = ["Marinelli", "Petri", "Dent", "point", "1L", "500mL", "100mL",
                   "25cm", "10cm", "5cm", "20cm"]
DETECTOR_TOKENS = ["HPGe", "NaI", "LaBr", "LaBr3", "CeBr", "CeBr3", "CZT", "CdZnTe"]
SAMPLE_TYPE_TOKENS = ["soil", "water", "air", "filter", "food", "concrete",
                      "metal", "calib", "check_source", "bkg", "background", "fon"]


def parse_filename(filename: str) -> dict:
    """Extract metadata tokens from filename. Returns dict of recognized fields."""
    name = Path(filename).stem
    out = {"nuclides": [], "geometry": "", "date": None, "sample_type": "",
           "detector": "", "is_background_hint": False, "raw": name}
    low = name.lower()

    # Nuclides
    for tok in NUCLIDE_TOKENS:
        if tok.lower() in low or re.search(rf"\b{re.escape(tok)}\b", name, re.IGNORECASE):
            out["nuclides"].append(tok)

    # Geometry
    for tok in GEOMETRY_TOKENS:
        if tok.lower() in low:
            out["geometry"] = tok
            break

    # Detector
    for tok in DETECTOR_TOKENS:
        if tok.lower() in low:
            out["detector"] = tok
            break

    # Sample type / background
    for tok in SAMPLE_TYPE_TOKENS:
        if tok.lower() in low:
            out["sample_type"] = tok
            if tok.lower() in ("bkg", "background", "fon"):
                out["is_background_hint"] = True

    # Date — try YYYYMMDD, YYYY-MM-DD, DD.MM.YYYY
    for pattern, fmt in [
        (r"(\d{4})(\d{2})(\d{2})", lambda m: datetime(int(m[1]), int(m[2]), int(m[3]))),
        (r"(\d{4})-(\d{2})-(\d{2})", lambda m: datetime(int(m[1]), int(m[2]), int(m[3]))),
        (r"(\d{2})\.(\d{2})\.(\d{4})", lambda m: datetime(int(m[3]), int(m[2]), int(m[1]))),
    ]:
        m = re.search(pattern, name)
        if m:
            try:
                out["date"] = fmt(m)
                break
            except (ValueError, KeyError):
                pass

    return out


# ============================================================================
# File readers
# ============================================================================

def read_spe(path: str) -> Spectrum:
    """Read Ortec/Canberra ASCII .spe file."""
    p = Path(path)
    text = p.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    spec = Spectrum(counts=np.array([]), live_time=0, real_time=0)
    spec.source_path = str(p)
    spec.filename_tokens = parse_filename(p.name)

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("$SPEC_ID:"):
            i += 1
            spec.sample_id = lines[i].strip() if i < len(lines) else ""
        elif line.startswith("$SPEC_REM:"):
            i += 1
            rems = []
            while i < len(lines) and not lines[i].startswith("$"):
                rems.append(lines[i].strip())
                i += 1
            spec.comments = "\n".join(rems)
            continue
        elif line.startswith("$DATE_MEA:"):
            i += 1
            try:
                spec.start_datetime = datetime.strptime(lines[i].strip(), "%m/%d/%Y %H:%M:%S")
            except (ValueError, IndexError):
                try:
                    spec.start_datetime = datetime.strptime(lines[i].strip(), "%d-%b-%y %H:%M:%S")
                except (ValueError, IndexError):
                    pass
        elif line.startswith("$MEAS_TIM:"):
            i += 1
            parts = lines[i].split()
            if len(parts) >= 2:
                spec.live_time = float(parts[0])
                spec.real_time = float(parts[1])
        elif line.startswith("$DATA:"):
            i += 1
            ch_range = lines[i].split()
            n0, n1 = int(ch_range[0]), int(ch_range[1])
            i += 1
            counts = []
            while i < len(lines) and not lines[i].startswith("$"):
                for tok in lines[i].split():
                    counts.append(int(tok))
                i += 1
            spec.counts = np.array(counts, dtype=np.float64)
            continue
        elif line.startswith("$ENER_FIT:") or line.startswith("$MCA_CAL:"):
            i += 1
            # Stored calibration — read but treat as a hint, per SKILL.md
            if line.startswith("$MCA_CAL:"):
                i += 1  # MCA_CAL has a count line first
            parts = lines[i].split()
            try:
                spec.energy_cal = tuple(float(x) for x in parts[:3] if x)
            except ValueError:
                pass
        elif line.startswith("$SHAPE_CAL:"):
            i += 1
            i += 1  # count line
            try:
                parts = lines[i].split()
                spec.fwhm_cal = tuple(float(x) for x in parts[:3])
            except (ValueError, IndexError):
                pass
        i += 1

    return spec


def read_csv(path: str) -> Spectrum:
    """Read a simple CSV/text spectrum. Tries to find counts column heuristically."""
    p = Path(path)
    spec = Spectrum(counts=np.array([]), live_time=0, real_time=0)
    spec.source_path = str(p)
    spec.filename_tokens = parse_filename(p.name)

    counts = []
    with p.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith(";"):
                # parse possible metadata in comments
                m = re.search(r"live[_ ]?time[:\s=]+([\d.]+)", line, re.IGNORECASE)
                if m:
                    spec.live_time = float(m.group(1))
                m = re.search(r"real[_ ]?time[:\s=]+([\d.]+)", line, re.IGNORECASE)
                if m:
                    spec.real_time = float(m.group(1))
                continue
            parts = re.split(r"[\s,;]+", line)
            try:
                # Last column = counts; assume integers or floats
                counts.append(float(parts[-1]))
            except ValueError:
                continue

    spec.counts = np.array(counts, dtype=np.float64)
    if spec.live_time == 0:
        # If unknown, fall back to setting live = real = 1 with explicit warning
        spec.live_time = 1.0
        spec.real_time = 1.0
    return spec


def read_spectrum(path: str) -> Spectrum:
    """Dispatch to the right reader by extension."""
    ext = Path(path).suffix.lower()
    if ext == ".spe":
        return read_spe(path)
    if ext in (".csv", ".txt", ".tsv"):
        return read_csv(path)
    raise NotImplementedError(f"Reader for {ext} not implemented in this skeleton. "
                              f"Extend read_spectrum() for .chn / .n42 / .mca / .spm.")


# ============================================================================
# Mariscotti second-derivative peak search
# ============================================================================

def mariscotti_search(counts: np.ndarray, fwhm_channels: float, sigma_threshold: float = 3.0):
    """
    Find peaks via second-derivative-of-Gaussian filter.
    Returns array of channel positions.
    """
    w = max(int(round(fwhm_channels)), 3)
    # half-width
    h = max(w // 2, 2)
    x = np.arange(-h, h + 1)
    sigma = fwhm_channels / 2.355
    g = (1 - (x / sigma) ** 2) * np.exp(-(x / sigma) ** 2 / 2)
    g = g - g.mean()
    g = g / np.sqrt((g * g).sum())

    conv = -np.convolve(counts, g, mode="same")
    # Background variance estimate
    bg = np.maximum(counts, 1.0)
    sigma_conv = np.sqrt(np.convolve(bg, g * g, mode="same"))

    significance = conv / np.maximum(sigma_conv, 1e-9)
    # Find local maxima of significance above threshold
    peaks, props = find_peaks(significance, height=sigma_threshold, distance=int(fwhm_channels))
    return peaks, significance[peaks]


# ============================================================================
# Bootstrap energy calibration
# ============================================================================

# Common anchor lines (keV) for natural background and standard sources
ANCHOR_LINES = {
    "K40": [1460.82],
    "Tl208": [583.19, 860.56, 2614.51],
    "Bi214": [609.31, 1120.29, 1764.49, 2204.21],
    "Pb214": [295.22, 351.93],
    "Pb212": [238.63],
    "Ac228": [338.32, 911.20, 968.97],
    "Cs137": [661.66],
    "Co60": [1173.23, 1332.49],
    "annihilation": [511.0],
    "Pb_XRF_Ka1": [74.97],
    "Pb_XRF_Ka2": [72.80],
    "Pb_XRF_Kb1": [84.94],
    "LaBr_138La": [788.74, 1435.80],
}


def fit_energy_calibration(channels: list, energies: list, degree: int = 2):
    """
    Polynomial fit E(N). Cap degree at 4.
    Returns (coefficients_low_to_high, residuals).
    """
    degree = min(degree, 4)
    n_points = len(channels)
    if n_points <= degree:
        degree = n_points - 1
    coeffs = np.polyfit(channels, energies, degree)
    # numpy returns high-to-low; invert
    coeffs = list(coeffs[::-1])
    predicted = sum(c * np.array(channels) ** i for i, c in enumerate(coeffs))
    residuals = np.array(energies) - predicted
    return tuple(coeffs), residuals


def check_stored_calibration(spec: Spectrum, anchor_peaks: list,
                              tolerance_factor: float = 0.3) -> bool:
    """
    anchor_peaks: list of (channel, true_energy_keV)
    Returns True if stored calibration's residuals < tolerance_factor * FWHM at every anchor.
    """
    if spec.energy_cal is None:
        return False
    for ch, true_E in anchor_peaks:
        predicted_E = spec.channel_to_energy(ch)
        fwhm = spec.fwhm_at(true_E) or 1.0
        if abs(predicted_E - true_E) > tolerance_factor * fwhm:
            return False
    return True


# ============================================================================
# FWHM(E) calibration
# ============================================================================

def fit_fwhm_hpge(energies, fwhms):
    """HPGe: FWHM^2(E) = a + b*E + c*E^2"""
    energies = np.array(energies)
    fwhms2 = np.array(fwhms) ** 2

    def model(E, a, b, c):
        return a + b * E + c * E * E

    popt, _ = curve_fit(model, energies, fwhms2, p0=[1.0, 1e-3, 1e-7])
    return tuple(popt)


def fit_fwhm_scint(energies, fwhms):
    """Scintillator: FWHM(E) = k * sqrt(E + alpha * E^2)"""
    energies = np.array(energies)
    fwhms = np.array(fwhms)

    def model(E, k, alpha):
        return k * np.sqrt(E + alpha * E * E)

    popt, _ = curve_fit(model, energies, fwhms, p0=[0.1, 1e-4])
    return tuple(popt)


# ============================================================================
# Confidence index CI (Lsrm §14.3)
# ============================================================================

def confidence_index(line_energies, energy_sigmas, ratio_sigmas=None):
    """
    CI = log10(1 / prod(delta_E_i * delta_I_j))
    """
    delta_E = [s / E for E, s in zip(line_energies, energy_sigmas) if E > 0]
    product = 1.0
    for d in delta_E:
        product *= max(d, 1e-6)
    if ratio_sigmas:
        for d in ratio_sigmas:
            product *= max(d, 1e-3)
    if product <= 0:
        return 0.0
    return math.log10(1.0 / product)


def dose_contribution(unident_peaks, all_peaks, efficiency_func):
    """
    DC = sum(S_i * E_i / eps(E_i)) over unidentified, divided by sum over all.
    """
    def dose(peaks):
        s = 0.0
        for p in peaks:
            if p.energy is None or p.energy <= 0:
                continue
            eps = max(efficiency_func(p.energy), 1e-6)
            s += p.area * p.energy / eps
        return s

    total = dose(all_peaks)
    if total <= 0:
        return 0.0
    return 100.0 * dose(unident_peaks) / total


# ============================================================================
# Currie / ISO 11929 MDA
# ============================================================================

def mda_iso11929(background_counts, t_live, efficiency, I_gamma_pct,
                 t_background=None, rel_unc_w=0.05, alpha=0.05, beta=0.05):
    """
    Compute MDA in Bq following ISO 11929:2019.
    background_counts = B = continuum count in the ROI
    t_background = background measurement time (default = t_live)
    """
    from scipy.stats import norm
    k_alpha = norm.ppf(1 - alpha)
    k_beta = norm.ppf(1 - beta)

    if t_background is None:
        t_background = t_live

    w = 1.0 / (efficiency * I_gamma_pct / 100.0)
    n0 = background_counts / t_background

    # L_C in counts (Lsrm formula 6.3-8)
    L_C = k_alpha * math.sqrt(background_counts * (1 + (t_live / t_background) ** 2))

    # L_D via quadratic equation (Lsrm formula 6.3-9)
    u_rel2 = rel_unc_w ** 2
    A = 1.0 - u_rel2
    B = -(2 * L_C + k_beta ** 2)
    C = (1 - k_beta ** 2 / k_alpha ** 2) * L_C ** 2
    disc = B * B - 4 * A * C
    if A == 0 or disc < 0:
        L_D = 2 * L_C + 2.71  # Currie fallback
    else:
        L_D = (-B + math.sqrt(disc)) / (2 * A)

    # MDA in Bq
    mda_bq = L_D * w / t_live
    return {"L_C": L_C, "L_D": L_D, "MDA_Bq": mda_bq, "w": w}


# ============================================================================
# Dead-time correction (Lsrm §15)
# ============================================================================

def correct_live_time(spec: Spectrum, A: float, B: float) -> float:
    """
    t_m = A * sum(y_i) + B * sum(y_i * i)
    Returns corrected live time.
    """
    y = spec.counts
    idx = np.arange(len(y))
    t_m = A * y.sum() + B * (y * idx).sum()
    return max(spec.live_time - t_m, 0.001)


# ============================================================================
# Main workflow stub
# ============================================================================

def analyze(spectrum_path: str, background_path: Optional[str] = None,
            user_candidates: Optional[list] = None) -> dict:
    """
    Top-level orchestrator. Follows the 11-step workflow in SKILL.md.
    Returns a structured result dict for report generation.
    """
    # Step 1 — Read file and metadata
    spec = read_spectrum(spectrum_path)
    bg = read_spectrum(background_path) if background_path else None
    if bg is not None:
        bg.is_background = True

    result = {
        "metadata": {
            "filename": spec.source_path,
            "sample_id": spec.sample_id,
            "start_datetime": str(spec.start_datetime) if spec.start_datetime else None,
            "live_time": spec.live_time,
            "real_time": spec.real_time,
            "dead_time_pct": 100 * spec.dead_time_fraction,
            "filename_tokens": spec.filename_tokens,
            "has_background_spectrum": bg is not None,
        },
        "steps_executed": [],
        "steps_skipped": [],
        "warnings": [],
    }

    # Step 2 — environment (placeholder — requires Pb XRF detection)
    # Step 3 — preliminary peak search
    # Default rough FWHM in channels — replace with paspport hint if available
    rough_fwhm = 5 if spec.fwhm_cal else 10
    peak_channels, peak_significances = mariscotti_search(
        spec.counts, fwhm_channels=rough_fwhm, sigma_threshold=3.0)
    peaks = [Peak(channel=float(ch), significance=float(s))
             for ch, s in zip(peak_channels, peak_significances)]
    result["peaks_found"] = len(peaks)
    result["steps_executed"].append("3. Preliminary peak search")

    # Step 4 — detector type (stub: needs FWHM in keV — defer until after step 6)
    # Step 5–6 — calibration (stub: needs anchor identification logic)
    # Step 7 — identification (stub: needs nuclear data library)
    # Step 8 — multiplet deconvolution (stub: needs identified nuclides)
    # Step 9 — MDA per ISO 11929 (use mda_iso11929 once efficiency calibration present)
    # Step 10 — secondary peaks
    # Step 11 — report

    result["status"] = ("skeleton — extend identification, calibration, and "
                        "multiplet deconvolution per SKILL.md")
    return result


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python analyze_spectrum.py <spectrum_file> [background_file]")
        sys.exit(1)
    out = analyze(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
    print(json.dumps(out, indent=2, default=str, ensure_ascii=False))
