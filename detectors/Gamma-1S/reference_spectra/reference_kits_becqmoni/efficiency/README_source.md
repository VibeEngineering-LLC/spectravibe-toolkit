# Lsrm SpectraLine Efficiency Files — Reference Data

This directory holds **example** efficiency files for the
**Gamma-1C NaI 63×63 USB №SN-01** detector, supplied by the
project maintainer on 2026-05-27.

**These files are DETECTOR-SPECIFIC.** They cannot be reused for
any other detector — even nominally identical hardware has different
efficiencies due to PMT gain, crystal manufacturing tolerances,
electronics chain, source-to-window distance, etc. Treat these files
as:

  • Test fixtures for the `.efa` / `.efr` parsers
  • Reference for the file-format schema
  • Source of representative efficiency curves and their structure
    (for documentation and validation)

For real measurements on a different detector, the user must
provide efficiency files calibrated for THAT specific detector.

## File format (Lsrm SpectraLine, v1.7.11918)

Two related formats — both CP-1251 (Windows-1251 Cyrillic) text
with CRLF line endings:

### `.efa` — Aggregated efficiency

One geometry per file. Contains a single header block followed by
energy points. Header is `[detector;geometry]`. Energy lines have
format:

    E_keV = epsilon, dEpsilon_pct, nuclide_source, S_counts, dS_counts, I_pct

where:
  - `E_keV` — gamma-ray energy
  - `epsilon` — efficiency value (decimal, photopeak)
  - `dEpsilon_pct` — relative uncertainty in efficiency (%)
  - `nuclide_source` — reference nuclide producing the line
  - `S_counts`, `dS_counts` — measured peak area and uncertainty
  - `I_pct` — gamma-line emission probability (%)

### `.efr` — Raw / per-source efficiency

Each reference source measurement is a separate block, each with its
own header `[detector;geometry;source-id]` and per-source metadata.
Energy lines have the same format as `.efa`. Used during efficiency
curve fitting to combine multiple source measurements.

## Geometries in the reference set

| Geometry | Volume | Distance | Material | Use case |
|---|---|---|---|---|
| **Маринелли** | 1000 ml | 0 (around detector) | ОИСН-16 ρ=1.6 | Large volume environmental samples |
| **Петри** | 60 ml | 0 (on top) | ОИСН-16 ρ=1.6 | Thin-layer samples (filters, swipes) |
| **Дента** | 120 ml | 0 (close geom.) | ОИСН-16 ρ=1.6 | Small dense samples |
| **Точечная-5см** | 0 | 5 cm | — | Calibration point source close |
| **Точечная-25см** | 0 | 25 cm | — | Calibration point source far |

The "ОИСН-16" matrix in volumetric geometries is the standard
Russian density-matched filler (density 1.6 g/cm³, iron-heavy
composition: 71.4% Fe, 20.6% C, 4.9% O).

## Use in the skill

The `gamma.calibration.efficiency` module (Phase 2.1c, future):
  • Parses `.efa` files to load epsilon(E) reference points
  • Fits efficiency curve ε(E) using log-log polynomial (Lsrm §8.5)
  • Optionally combines multiple source measurements from `.efr`
    for finer point density
  • Exposes ε(E) callable to identification, MDA, activity modules

For Phase 2.1b (multiplet deconvolution), efficiency files are
NOT required — they enter only when computing activities (Phase
2.1c/d).
