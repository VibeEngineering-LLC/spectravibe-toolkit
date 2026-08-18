# Dead-time correction and high-count-rate handling

At high count rates, the spectrometer's effective dead time exceeds what's recorded as analyzer dead time. Random pulse pile-up also distorts areas of strong peaks. Both effects need correction before quantitative activity calculation.

Method follows Lsrm *Algorithmic Foundations* §15.

## Effective dead-time model

Two contributions:
1. Random summing in the analyzer (load-dependent)
2. Pulse-pair resolution loss (also load-dependent)

Lsrm models the combined effective dead time per measurement as:
```
t_m = A · Σy_i + B · Σ(y_i · i)
```
where:
- `y_i` is the count in channel `i` of the spectrum
- `A` and `B` are detector-specific empirical coefficients
- `Σy_i` is the total count (proportional to load)
- `Σ(y_i · i)` weights by channel number (proportional to load × mean energy — captures the energy dependence of pulse-pair resolution)

The corrected live time is:
```
t_live_corrected = t_live − t_m
```
This corrected live time is then used in all activity calculations.

## Calibrating A and B

Three measurements with known sources are needed:

1. **Reference (low-load)**: a single high-energy emitter (typically ⁶⁰Co) at a distance giving < 500 cps total load. Measure long enough that the area of the 1173 keV (or 1332 keV) line has < 1% statistical uncertainty. Compute the **reference rate** `n₀` of the chosen reference line (its count rate at low load, presumed unaffected by pile-up).

2. **Low-energy loading**: keep the ⁶⁰Co source in place; add a low-energy emitter (e.g., ¹³³Ba or ²⁴¹Am) close enough to bring total load to ~5·10⁴ cps. Measure and compute:
 - `n₁` = count rate of the ⁶⁰Co reference line (now suppressed by pile-up)
 - `r₁` = total count rate (Σy_i / t)
 - `p₁` = Σ(y_i · i) / t

3. **High-energy loading**: same procedure but with a high-energy emitter (e.g., ¹⁵²Eu or ¹³⁷Cs) giving ~5·10⁴ cps total. Compute `n₂`, `r₂`, `p₂`.

Relative counting losses:
```
Δ₁ = (n₀ − n₁) / n₀
Δ₂ = (n₀ − n₂) / n₀
```

Solve the 2×2 system for A and B (Lsrm §15.1):
```
A = (Δ₁ · p₂ − Δ₂ · p₁) / (p₂ · r₁ − p₁ · r₂)
B = (Δ₂ · r₁ − Δ₁ · r₂) / (p₂ · r₁ − p₁ · r₂)
```

Store {A, B} per detector in a calibration file. The method is validated up to 5·10⁴ cps with < 5% residual uncertainty.

## Decision rules for applying the correction

- **Dead time < 5%**: skip correction. Default Poisson statistics apply.
- **Dead time 5–30%**: apply correction if A, B are known for this detector. Otherwise flag activities as uncorrected and report the dead time prominently.
- **Dead time > 30%**: correction is essential. If A, B unknown, **ask the user** whether to proceed without correction (results will be biased low) or to defer until A, B are calibrated. Do not silently proceed.

## Pile-up diagnostics in the spectrum itself

Independently of the {A, B} correction, look for these signatures in the spectrum:

1. **Sum-of-self peak** at 2·E for every strong line. Width should be ~√2·FWHM(E). If a peak appears at 2·E with width matching √2·FWHM(E), and no isotope has a line there, it's pile-up.
2. **High-energy continuum tail**: the spectrum has elevated counts above the highest real γ-line. This is the continuous distribution of random sums.
3. **Distorted intensity ratios for cascade nuclides** (e.g., ⁶⁰Co 1173/1332 ratio off from 1.00 after efficiency correction).

If any of these are visible despite the {A, B} correction being applied, the load is exceeding the model's range. Reduce source-to-detector distance, use pile-up rejection electronics, or accept higher uncertainty.

## Distinction from True Coincidence Summing (TCS)

Critical not to confuse with TCS:

| Property | Pile-up (random) | TCS (true coincidence) |
|----------|------------------|------------------------|
| Source | Two unrelated decays in the resolving time | Two γ-rays from the **same** decay cascade |
| Scales with | rate² | (geometry × cascade probability), rate-independent |
| Mitigated by | Lower load, pile-up rejection | Larger source-detector distance |
| Affects which nuclides? | Any, in proportion to their rate | Only cascade emitters (⁶⁰Co, ¹³³Ba, ¹⁵²Eu, ¹³⁴Cs, ²²Na, ⁶⁸Ge/Ga, etc.) |
| Correction model | A, B coefficients (Lsrm §15) | TCS factor per line per geometry (Lsrm §17 + nuclear data) |

Both can be present simultaneously. The order of corrections: first TCS (geometry-dependent, applied per-line), then dead-time / pile-up (load-dependent, applied globally to live time).

## TCS correction (brief — see Lsrm §17)

For cascade nuclides at close geometry, the apparent count rate at each line is altered by coincident summing. The correction factor for line *i* is approximately:
```
K_TCS,i ≈ 1 / Π_j (1 − P_j · ε_total(E_j))
```
where the product runs over all γ-rays coincident with line *i* in the decay scheme, and `ε_total(E_j)` is the **total efficiency** (peak + Compton) at energy E_j.

In practice, K_TCS is computed by Monte Carlo simulation (Lsrm `TccfCalc` module or equivalent: PENELOPE, GEANT4, MCNP) using the ENSDF decay scheme and a model of the detector and source geometry. For routine work, the simplest mitigation is to measure at a far geometry (≥ 25 cm for standard HPGe) where TCS is < 1–2% and corrections become unnecessary.

If TCS is significant and not corrected, the intensity-ratio check in step 7C will fail for cascade nuclides — they will look like they have wrong ratios. This is a diagnostic signal that TCS correction is needed.
