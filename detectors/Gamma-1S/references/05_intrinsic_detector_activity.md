# Intrinsic detector activity

Every detector crystal has its own background — sometimes negligible (HPGe), sometimes substantial enough to dominate parts of the spectrum (LaBr₃). Failure to recognize intrinsic peaks leads to false nuclide identifications. **For sub-Bq sample activities, the intrinsic background is often the limiting factor.**

This reference catalogs what to expect from each common crystal type.

## LaBr₃(Ce) — substantial intrinsic activity (~1.4–1.5 Bq/cm³)

### From ¹³⁸La (natural lanthanum)

¹³⁸La is a primordial isotope with natural abundance 0.0888% in lanthanum and T½ = 1.02 × 10¹¹ years. It decays two ways:

**(a) Electron capture (66.4%) → ¹³⁸Ba***
- Emits 1435.795 keV γ-ray (de-excitation of ¹³⁸Ba 2⁺ state)
- Coincident: Ba K X-rays (from filling the K-shell vacancy left by EC)
  - Ba Kα₂ 31.82 keV
  - Ba Kα₁ 32.19 keV
  - Ba Kβ₁ 36.38 keV
  - Ba Kβ₃ 36.30 keV
  - Ba Kβ₂ 37.26 keV
  - **Weighted-average Kα ~ 32.1 keV, Kβ ~ 36.4 keV** — these are emitted *inside* the LaBr₃ crystal itself, so depending on escape probability:
    - If X-ray fully absorbed in crystal: γ + X-ray summed → peak at 1435.8 + 37.4 = **1473 keV** (Ba K binding energy adds to γ)
    - If X-ray escapes: peak at 1435.8 keV, escape peak at ~1441 keV (1473 − 32 keV escape)
    - Ba K X-ray peaks at 32 and 36 keV from outside-of-cascade events

**(b) β⁻ (33.6%) → ¹³⁸Ce***
- Emits 788.742 keV γ-ray
- Coincident with β⁻ continuum (endpoint 255 keV)
- The β contributes a continuous background 0–255 keV

### From ²²⁷Ac chain contamination

Trace ²²⁷Ac (chemically similar to lanthanum) is co-extracted during crystal growth. T½ = 21.77 y. Daughters (²²⁷Th, ²²³Ra, ²¹⁹Rn, ²¹⁵Po, ²¹¹Bi, ²⁰⁷Tl, stable ²⁰⁷Pb) include 6 α-emitters and 4 β-emitters.

α-particles in LaBr₃ produce light at a different yield than γ — the effective "γ-equivalent" position differs. For a high-resolution analyzer the α peaks land at:
- ²²⁷Th 5.7–6.0 MeV α → ~1.65–1.75 MeV γ-equivalent
- ²²³Ra 5.6–5.7 MeV α → ~1.6 MeV γ-equivalent
- ²¹⁹Rn 6.8 MeV α → ~2.0 MeV γ-equivalent
- ²¹⁵Po 7.4 MeV α → ~2.2 MeV γ-equivalent
- ²¹¹Bi 6.6 MeV α → ~1.9 MeV γ-equivalent

Total ²²⁷Ac chain contribution: ~0.1–0.15 Bq/cm³.

### Summary of LaBr₃ intrinsic peaks (always present)

| Energy (keV) | Origin | Type label |
|--------------|--------|-----------|
| 32.1 | Ba Kα (from ¹³⁸La EC) | intrinsic_La138_BaXray |
| 36.4 | Ba Kβ (from ¹³⁸La EC) | intrinsic_La138_BaXray |
| 788.7 | ¹³⁸La β⁻ branch γ | intrinsic_La138_789 |
| 1435.8 | ¹³⁸La EC branch γ | intrinsic_La138_1436 |
| ~1441 | ¹³⁸La EC, X-ray escape | intrinsic_La138_escape_1441 |
| 1473 | ¹³⁸La EC, X-ray summed | intrinsic_La138_sum_1473 |
| ~1.5–2.2 MeV γ-eq | ²²⁷Ac chain α-peaks | intrinsic_Ac227_alpha |
| 0–255 (continuum) | ¹³⁸La β⁻ | intrinsic_La138_beta |

**Implication for analysis:** any γ-line in the sample at 32, 36, 789, 1436, 1441, or 1473 keV must be checked against intrinsic level before being assigned to a sample nuclide. ¹³⁷Cs is particularly tricky on LaBr₃ — the Ba K X-rays from ¹³⁷ᵐBa internal conversion overlap exactly with the intrinsic Ba X-rays.

## CeBr₃ — modest intrinsic activity

No ¹³⁸La (cerium is stable Ce). The main intrinsic source is ²²⁷Ac contamination from crystal growth, similar to LaBr₃ but lower:
- ²²⁷Ac chain: ~0.02–0.1 Bq/cm³
- α-peaks at the same γ-equivalent energies as LaBr₃ (~1.5–2.2 MeV)
- No γ-lines at 789 or 1436

CeBr₃ is generally preferred over LaBr₃ when intrinsic activity matters (low-count-rate work, high MDA requirements).

## NaI(Tl) — low intrinsic activity but iodine artifacts

The sodium iodide crystal itself has very low intrinsic activity (trace ⁴⁰K from sodium handling is usually < 0.01 Bq/cm³). But:

### Iodine K-escape peaks

When a low-energy γ photoelectrically interacts in NaI and the resulting iodine K-X-ray escapes the crystal, the recorded energy is reduced. For any γ-line of energy E > I K-binding (33.17 keV):

- Iodine Kα-escape: peak at **E − 28.6 keV** (weighted-average Kα)
- Iodine Kβ-escape: peak at **E − 32.3 keV** (weighted-average Kβ)
- Escape intensity decreases rapidly with E above ~50 keV (deeper interaction → X-ray more likely to be absorbed before escaping)

Most visible for low-energy sources (²⁴¹Am 59.5 keV → escapes at 27 and 31 keV).

### Iodine fluorescence

External γ-rays from any source (background, sample) photoelectrically excite iodine in the crystal, producing the I K X-rays themselves as a peak in the spectrum:
- I Kα weighted-average: 28.6 keV
- I Kβ weighted-average: 32.3 keV

These appear as small peaks at low energy in essentially every NaI spectrum when there is significant low-energy γ flux.

### Ba K X-rays from ¹³⁷Cs (not intrinsic, but cooked into the source signature)

For ¹³⁷Cs samples, the Ba K X-rays at 32, 36 keV from internal conversion of ¹³⁷ᵐBa are real source emission, **not intrinsic**. But they look identical to LaBr₃ intrinsic. In NaI they identify the source as a ¹³⁷Cs source.

### ⁴⁰K trace

Sometimes a ⁴⁰K 1460 keV line appears at very low level in NaI spectra, originating from potassium in:
- The crystal package (glass envelope, PMT)
- Light reflector wrapping
- The K used in some PMT photocathodes
- Concrete/walls nearby (more likely environmental than detector)

Discriminating is hard — measure a clean background, subtract.

## HPGe — negligible intrinsic activity, but Ge K-escape

The germanium crystal is highly purified (intrinsic activity is essentially zero — limited by trace ²³²Th, ²³⁸U in the contact materials and cryostat, typically nBq/cm³).

### Ge K-escape (low energy only)

For γ-lines below ~120 keV, photoelectric absorption near the crystal surface can lead to escape of Ge K X-rays:
- Ge Kα weighted: 9.886 keV
- Ge Kβ weighted: 10.98 keV

Escape peaks at E − 9.89 keV (visible) and E − 11.0 keV (usually too close to the main peak to resolve unless E is large).

Most visible for ¹²⁵I 35.5 keV, ²⁴¹Am 59.5 keV, ¹³³Ba 81 keV, ⁵⁷Co 122 keV.

### Other HPGe intrinsic features

- Ge characteristic X-rays from cosmic-ray neutron activation of Ge → ⁷⁵ᵐGe, ⁷¹ᵐGe metastable states (rare, only in deep-underground or long-exposure backgrounds)
- ⁷⁰Ge(n,γ)⁷¹Ge during transport through cosmic flux (very small)

For practical purposes, HPGe is intrinsic-free.

## CdZnTe (CZT)

Cd K-X-ray escape: **E − 23.17 keV** (Cd Kα weighted-average); **E − 26.10 keV** (Cd Kβ).
Te K-X-ray escape: **E − 27.47 keV** (Te Kα weighted-average); **E − 31.0 keV** (Te Kβ).

For low-energy γ (< ~150 keV), each photopeak typically has a Cd Kα-escape and Te Kα-escape companion peak at lower energy, with relative intensities depending on detector size and γ energy.

### ¹¹³Cd β-decay

¹¹³Cd, natural abundance 12.22%, has T½ = 7.7 × 10¹⁵ y (β⁻ to ¹¹³In, E_β,max = 316 keV). Theoretical activity in pure CdZnTe is ~0.4 Bq/cm³, but practical detectors are below this because of isotopically depleted Cd or low-volume detection. Contributes a continuous β background up to ~316 keV.

### ¹¹³ᵐCd

A small fraction of Cd in some CZT batches is ¹¹³ᵐCd (metastable, T½ 14.1 y, β⁻ to ¹¹³In with E_β,max 580 keV). Generally negligible.

## Si and HPGe planar detectors

For X-ray spectrometry detectors (Si(Li), Si-PIN, SDD, HPGe planar):
- Si K-escape: **E − 1.74 keV** (Si Kα). Visible for any line above the Si K-edge.
- Ge K-escape: **E − 9.89 keV**. Same as above.

## Summary table: intrinsic peaks to always check

| Detector | Always check for | Origin | Diagnostic |
|----------|----------------|--------|-----------|
| HPGe (coaxial) | nothing | — | clean detector |
| HPGe (planar/LEGe) | E_γ − 9.89 keV companions | Ge K-escape | low-E γ |
| NaI(Tl) | E_γ − 28.6 / E_γ − 32.3 | I K-escape | low-E γ |
| LaBr₃(Ce) | 32, 36, 789, 1436, 1441, 1473 keV; broad continuum 0–255 keV; α-peaks 1.5–2.2 MeV | ¹³⁸La + ²²⁷Ac | always present, ~1.5 Bq/cm³ |
| CeBr₃ | α-peaks 1.5–2.2 MeV | ²²⁷Ac | low level |
| CdZnTe | E_γ − 23.17, E_γ − 27.47 keV companions | Cd/Te K-escape | low-E γ; β continuum to 316 keV |

When reporting (step 11), the diagnostic block should state the **measured intrinsic activity level** for scintillator detectors based on the ¹³⁸La 1436 (LaBr₃) or ²²⁷Ac α-peak integrals (LaBr₃, CeBr₃), and compare to literature values (~1.4 Bq/cm³ for LaBr₃).
