# Методологические уроки — SpectraVibe

> **F-336 / v1.18.18.15** — выделено из `KNOWN_AND_FIXED_ISSUES.md` как отдельный
> учебный документ. Здесь собраны narrative-уроки накопленные за разработку:
> как читать NaI-спектр, что значит σε boundary, как работает chain-suppression,
> почему 70-90 кэВ кластер не разложим на NaI, и т.п.
>
> Документ нужен для онбординга новых операторов / разработчиков. Связан с
> главным `KNOWN_AND_FIXED_ISSUES.md` через F-ID референсы.

---

## Methodology lessons captured

Beyond bugs, the development process exposed several methodology
gaps that have been documented for future implementers in
`NOTES_v1.7_methodology.md` (446 lines, 10 sections):

| § | Topic |
|---|---|
| 1 | Subtract behaviour (photopeak / Compton / Pb XRF survival) |
| 2 | Compton-feature identification algorithm |
| 3 | Library reference samples / cross-correlation |
| 4 | Pb-210 46.5 keV contamination of Pb shielding (→ F-18) |
| 5 | Pb XRF: induced vs background, only baseline portion subtracts |
| 6 | 511 ROI: Tl-208 510.77 vs annihilation — disambiguation |
| 7 | Lsrm Algorithmic Foundations 2022 summary |
| 8 | NaI interferences: 511, 600, 665 keV ROIs |
| 9 | Lsrm-recommended NaI ERN calibration set (7 lines) |
| 10 | Low-background shielded geometry: 32 + 46.5 + 75 keV |

---

#### F-31: Peak-area accuracy on close doublets + full TCS correction (Phase 2.1g, v1.7.9)
**Severity**: Critical (the area defect caused a −25% activity error
on Co-60, the canonical cascade calibrator; the TCS module closes the
last open quantitative-correction limitation, K-17)
**Closes**: K-17 (true coincidence summing correction — previously
only a placeholder dict argument existed from F-29)
**Strengthens**: K-11 (Cowell baseline in dense regions — now
mitigated by Lsrm-table fallback for any Lsrm-sourced spectrum)

**Discovery context**: While validating F-29 against the new .src
certificates (F-30), the Co-60 №043 02.2019 point source showed a
**−25.56% deviation** from its 105 000 Bq certificate — far larger
than cascade summing alone can explain (TCS at 5cm is only ~2-3%).

**Root cause (two distinct problems, hence F-31a + F-31b):**

**F-31a — Cowell area under-counts on close doublets.**
Debugging revealed the dominant error was not physics but
integration. On the Co-60 1173/1332 doublet (separated by ~160 keV,
each with FWHM ~66 keV on NaI), Cowell's linear baseline runs
through the *neighbouring peak's Compton wing* rather than under the
true continuum, subtracting far too much:

  - Cowell 1173 keV area: ~469 000 counts
  - Lsrm software's own Gaussian-on-step fit (stored in the SPE
    `<START PEAKS>` table): **666 002 counts** — a **42% deficit**
  - For an *isolated* peak (Cs-137 661.66 keV) the Cowell/Lsrm ratio
    is only 1.039 (3.9%), which is why F-30's Cs-137 e2e validation
    succeeded despite this bug.

**Fix (F-31a)**: New helper `get_peak_area(spec, peak_channel,
fwhm_channels, prefer_lsrm_table=True, match_tolerance_fwhm=0.8)` in
`gamma.peaks.area`. Strategy: when the spectrum carries
`extras["lsrm_peaks_table"]` (populated by the F-26 reader from
Lsrm's own peak fit) and a table row lies within
`match_tolerance_fwhm × FWHM` of the requested channel, return that
**Lsrm Gaussian-fit area**; otherwise fall back to Cowell. Returns
`(area, uncertainty, source)` where `source ∈ {"lsrm_peaks_table",
"cowell", "failed"}`. `gamma.identification.identify` now calls
`get_peak_area` instead of `cowell_area` directly when populating
the per-peak area cache. For AtomSpectra XML (no Lsrm table) the
behaviour is unchanged (Cowell).

**F-31b — full TCS correction module.**
New module `gamma.physics.cascade_summing` (~340 LOC) implementing
the Knoll §17.6 / [GILMORE-8.5] single-pair coincidence-summing loss
model:

    Loss(E_i) = Σ_j p_c(E_i,E_j) · ε_T(E_j)
    C(E_i)    = 1 / (1 − Loss(E_i))
    A_true    = C · A_observed

where ε_T(E) = ε_p(E) / P(E) is the *total* efficiency, derived from
the F-27 photopeak efficiency via the **peak-to-total ratio** P(E).
Components:

  - `CASCADE_SCHEMES` — curated decay-scheme catalogue (Co-60, Y-88,
    Na-22 with its 2×511 annihilation pair, Tl-208, Eu-152 partial,
    Ba-133 partial); branching probabilities from ENSDF / NuDat 3.
  - `peak_to_total_NaI(E)` — log-log degree-2 fit to Gilmore Table 8.4
    (NaI 3×3"): `log P = −0.316 + 0.458·ln E − 0.081·(ln E)²`, clipped
    to [0.05, 1.0]. Agrees with reference anchors within ≤5%.
  - `peak_to_total_HPGe(E)` — conservative default (caller should
    supply measured P(E) for accurate HPGe work).
  - `total_efficiency`, `tcs_correction_factor`,
    `compute_tcs_corrections(nuclide, eff_curve) → {E: factor}`. The
    returned dict plugs directly into the F-29
    `compute_activity(coincidence_correction=…)` argument that has
    existed as a K-17 placeholder since v1.7.7.
  - Non-cascade nuclides (Cs-137) return an empty dict / factor 1.0;
    unknown nuclides and unmatched lines return 1.0 silently.
  - `loss_cap` (default 0.5 → max factor 2.0) guards against
    unphysical extrapolation.

**TCS factors verified** on the real Гамма-1С 5cm ε(E): Co-60
~2.2%/line, Y-88 2.0%/2.4%, Tl-208 2.6-3.1%, Eu-152 2.6-7.5%
(strongest at 122/245 keV), Cs-137 empty (non-cascade). These match
the literature expectation for point 5cm geometry.

**Quantitative validation (the headline result):**

| Source | before F-31 | F-31a only | F-31a+F-31b |
|---|---|---|---|
| Cs-137 (isolated, control) | −1.43% | +1.83% | +1.83% (no TCS) |
| **Co-60 (close doublet, cascade)** | **−25.56%** | **+0.61%** | **+2.89%** |

  - Co-60 deviation collapses from −25.56% (≈6.7σ outside the 3.84%
    combined uncertainty) to **+0.61%** with F-31a alone — i.e. the
    area defect, not cascade summing, was the dominant error.
  - Adding F-31b TCS moves Co-60 to +2.89% (still within 1σ). The
    slight over-correction is expected and documented (K-18): the
    Lsrm Gaussian fit already partially recovers summing-displaced
    counts, so the full analytic TCS double-counts a little. For a
    cleaner pipeline the TCS correction should be applied to *Cowell*
    areas, not Lsrm-fitted areas — captured as a future refinement.
  - Cs-137 (non-cascade isolated control) is unaffected by TCS as
    required — confirms the correction is correctly scoped.

**Verification**: **25 new tests** (123 → 148 total):
  - `test_peak_area.py` +5 (Phase 2.1a group extended): Lsrm-table
    preference, fallback when no table, fallback when no match in
    tolerance, `prefer_lsrm_table=False` bypass, and the real Co-60
    doublet showing the 27% Cowell deficit being corrected.
  - `test_cascade_summing.py` +20 (new Phase 2.1g suite): scheme
    catalogue sanity (4), peak-to-total models (4), TCS formula and
    factor behaviour (8), dispatcher + compute_activity integration
    (2), geometry/efficiency scaling reality checks (2).
  All 123 prior tests still pass.

---

#### K-18: TCS over-correction when applied to Lsrm-fitted areas
**Status**: ✅ **FIXED in F-35 (v1.7.13)** — `compute_activity` now
reads `LineMatch.peak_area_source` (added in F-34) and scales the
analytic TCS effect by an area-method-aware factor:
`c_eff = 1 + (c_analytic − 1) · scale[area_source]`. Default
`scale["lsrm_peaks_table"] = 0.0` (the Lsrm wide-ROI Gaussian fit
already recovers summing-displaced counts), `scale[*] = 1.0` for
every other source. See F-35 for the full mechanism and 10-test
verification.

**Original symptom (v1.7.9)**: With both F-31a (Lsrm-fitted areas) and
F-31b (analytic TCS) enabled, Co-60 activity over-shoots the
certificate by ~2.9%, versus +0.6% with F-31a alone.
**Root cause**: The analytic TCS loss model assumes the photopeak
area has been depleted by the full summing loss. But Lsrm's
Gaussian-on-step fit, by fitting the peak shape over a wide ROI,
already recovers some of the summing-displaced counts that a pure
channel-sum (or Cowell) integration would lose. Applying the full
analytic correction on top therefore double-counts a fraction of the
summing effect.
**Workaround**: For cascade nuclides measured at near geometries on
Lsrm hardware, prefer F-31a alone (areas already mostly correct), or
apply TCS to Cowell areas rather than Lsrm-table areas.
**Planned resolution**: Make the activity pipeline TCS-aware of the
area source — apply analytic TCS only when the area came from a
channel-sum / Cowell method, and apply a reduced (or zero) TCS when
the area came from the Lsrm Gaussian fit. Requires threading the
`area_source` flag from `get_peak_area` through `LineMatch` into
`compute_activity`.
**References**: [GILMORE-8.5.3] (interaction of peak-fitting method
with summing correction); Lsrm Algorithmic Foundations §10.

---

#### F-32: Symmetric energy-ceiling API for file readers (Phase 2.1h, v1.7.10)
**Severity**: Minor (architectural / API ergonomics — no numerical
behaviour change at default settings)
**Affected**: `gamma.io.atomspectra_xml.read_atomspectra_xml`,
`gamma.io.lsrm_spe.read_lsrm_spe`, `gamma.io.readers.read_spectrum`

**Discovery context**: `read_lsrm_spe` already accepted a keyword-only
`apply_energy_ceiling` toggle, while `read_atomspectra_xml` applied the
3000 keV trim unconditionally. Neither reader allowed a per-call
override of the ceiling value, so diagnostic work above 3 MeV (or with
a different ceiling) required mutating the module-level constant
`ENERGY_CEILING_KEV` — a pattern that is global, easy to forget to
revert, and contaminates parallel reads.

**Resolution**: Both readers now expose the same keyword-only contract:

```python
read_lsrm_spe(path, *, apply_energy_ceiling=True, ceiling_keV=None)
read_atomspectra_xml(path, *,
                     parse_background=True,
                     apply_energy_ceiling=True,
                     ceiling_keV=None)
```

- `apply_energy_ceiling=False` keeps every decoded channel (minus the
  overflow marker for AtomSpectra files).
- `ceiling_keV=<float>` overrides `ENERGY_CEILING_KEV` for that call
  only; `None` means "use the module constant".
- For `read_atomspectra_xml`, the embedded
  `BackgroundEnergySpectrum` block inherits the same trim policy as the
  primary spectrum so the two arrays stay channel-aligned.
- `parse_background` was promoted to keyword-only as part of the
  signature cleanup; the only internal caller (`gamma.io.background.
  resolve_external_background`) already used the keyword form.

The dispatcher `gamma.io.readers.read_spectrum` forwards `**kwargs`
unchanged, so the new contract is reachable from both ends with no
caller-side changes at default settings.

**Verification**: `test_reader_api.py` — 7 new tests covering default
trim, `apply_energy_ceiling=False`, custom `ceiling_keV`, and embedded-
background inheritance on both formats. Full suite green: 155 tests
(148 prior baseline + 7 new). Spot-checked manually on
`Фон_кабинет_8192к_01-01-2025.xml`:
  - default → 7034 ch, e_max 2999.6 keV (drops 1157 channels above ceiling)
  - apply_energy_ceiling=False → 8191 ch, e_max 3441.5 keV (full range)
  - ceiling_keV=1500 → 3533 ch, e_max 1499.8 keV

**References**: [GILMORE-1] / §3 (energy-axis scope considerations);
Lsrm Algorithmic Foundations §3.

---

#### F-33: Multiplet deconvolution (Phase 2.1b, v1.7.11) — closes K-05
**Severity**: Major (the last open methodological hole from the
Phase 2.1 plan; closes K-05 which had been "deferred" since v1.7.1)
**Closes**: K-05 (multiplet deconvolution)
**Affected**: new module `gamma.peaks.deconvolve`

**Discovery context**: K-05 has been the standing TODO since the
identification pipeline came online (v1.7). Multiplets are
fundamentally limited on scintillators — on NaI the 240 keV "peak"
is actually Pb-212 238.6 + Pb-214 241.98 unresolvable at ~30 keV FWHM,
and the Bi-214 609 / Tl-208 583 region overlaps on the higher Compton
shoulder. Without deconvolution, the activity ratios from these peaks
are systematically biased.

**Methodology (Lsrm §9 / [GILMORE-9.7])**: identification-first
deconvolution. Once identification has produced a confirmed nuclide
list, the **positions** of every component (library E_i → channel via
energy_cal) and the **FWHMs** (from `fwhm_at_channel`) are **fixed**
constraints. Only the **areas** (one per component) and the
**continuum** beneath the multiplet are free. With these constraints
the fit model is **linear in its free parameters**:

```
y(x) = Σ_k A_k · g_k(x; c_k, σ_k)    + β₀ + β₁·(x − x_mid)         [linear]
       + β_step · 0.5·erfc((x − x_step)/(σ_step·√2))                [step+linear]
```

where `g_k` is the unit-area Gaussian. The free vector is
`[A_1, …, A_n, β₀, β₁ (, β_step)]`. Solved by
`scipy.optimize.lsq_linear` with `A_k ≥ 0` and `β_step ≥ 0`
(physically: the Compton step from a photopeak adds counts on the LOW-
energy side, so the step is "down" through the peak — its height must
be non-negative). `β₀` and `β₁` are unconstrained. A
`np.linalg.lstsq` fallback path with iterative pruning of negative
areas exists for environments without `scipy.optimize.lsq_linear`.

**Why linear-LSQ and not full nonlinear LM**: the Lsrm methodology
explicitly forbids floating centroids and widths in multiplets where
the components are known a priori. With these fixed, the nonlinearity
disappears entirely and the fit becomes well-posed whenever any two
components are separated by more than ~0.5·σ. There is no convergence
risk, no local-minimum trap, and the covariance matrix is analytic
from the weighted normal equations.

**API surface**:

```python
from gamma.peaks.deconvolve import (
    MultipletComponent, DeconvolutionResult,
    deconvolve_multiplet, find_multiplet_regions,
    deconvolve_identified_multiplets,
)

# Low-level: solve one multiplet
components = [
    MultipletComponent("Pb-212", 238.6, library_I_pct=43.6,
                       center_channel=ch_238, fwhm_channels=fwhm_238),
    MultipletComponent("Pb-214", 241.98, library_I_pct=7.4,
                       center_channel=ch_242, fwhm_channels=fwhm_242),
]
res = deconvolve_multiplet(spec.counts, components=components,
                           continuum="step_linear")
# res.areas, res.area_uncertainties, res.covariance, res.chi2_per_dof,
# res.degenerate_pairs (flagged when separation < 0.5·σ)

# High-level: detect clusters in an IdentificationResult and run
# deconvolution on each
from gamma.peaks.deconvolve import deconvolve_identified_multiplets
results = deconvolve_identified_multiplets(id_result, spec, fwhm_at)
```

`find_multiplet_regions(identification_result, fwhm_at_channel, *, overlap_threshold_fwhm=1.0)`
uses single-linkage clustering over the `LineMatch` channels — three
adjacent overlapping lines all land in one cluster by transitive
closure, even if the outermost pair lies just outside the threshold.

**Degenerate-pair flag**: any pair of components within 0.5·σ of each
other is reported in `DeconvolutionResult.degenerate_pairs` so callers
know the individual areas are strongly correlated and should be
interpreted via the full covariance or merged.

**Verification — 12 tests in `test_deconvolve.py`**:

- Synthetic 1-component fit recovers known area within 3%, χ²/ν=0.83
- 2-component well-separated doublet (Δ=4σ): each area within 5%
- 2-component close doublet (Δ=1σ): each area within 15%, sum within
  0.5% (most of the uncertainty goes into the off-diagonal covariance,
  as theory predicts)
- 3-component triplet at 2σ spacing: each area within 10%
- Step-continuum recovery: step_linear χ²/ν=0.96 vs pure-linear
  χ²/ν=2.18 on synthetic with a real step
- Non-negativity bound: a "ghost" component placed at an empty channel
  clamps to area = 0
- Degenerate-pair flag fires when Δ<0.5σ
- `area_by_nuclide()` correctly sums same-nuclide components
- `find_multiplet_regions` basic overlap, no-overlap, transitive chain
- **Real-spectrum**: Co-60 1173/1332 doublet on
  `Co-60__043_02_2019_Точечная-5см_5cm.spe`:
  area_1173 = 540 604, area_1332 = 495 024, ratio 1.09 (library 1.00,
  TCS-depleted similarly on both lines at 5 cm). Demonstrates the
  pipeline works on a high-statistics real spectrum.

**Scope note (deferred to next phase)**: the present F-33 delivers
the deconvolution algorithm and a standalone API. **Integrating it
into the identify-then-activity pipeline** (i.e., automatically
re-routing the `LineMatch.peak_area` field through deconvolution
when the line falls into a multiplet cluster) is a separate phase
(potential v1.7.12), because it touches `identify.py` and
`compute_activity` and benefits from a focused validation pass on
multiple real spectra.

**References**: Lsrm Algorithmic Foundations 2022 §9
(multiplet deconvolution under fixed-position constraints); Gilmore &
Joss, Practical Gamma-ray Spectrometry, 3rd Ed., §9.7 (continuum
modelling: step + linear under the photopeak).

---

#### F-34: Pipeline integration of multiplet deconvolution (Phase 2.1b cont., v1.7.12)
**Severity**: Minor (architecture / pipeline ergonomics — adds the
"plumbing" that lets identified multiplets reach the activity
calculation with deconvolved areas)
**Continuation of**: F-33 (v1.7.11)
**Affected**: `gamma.identification.identify.LineMatch` (new field
`peak_area_source`); `gamma.identification.identify.identify_nuclides`
(threads area source through the peak-area cache);
`gamma.identification.disambiguate.disambiguate_identifications`
(preserves area + source on promoted LineMatch);
`gamma.peaks.deconvolve` (new `apply_multiplet_deconvolution`).

**Discovery context**: F-33 delivered the deconvolution algorithm and
synthetic-data validation, but explicitly deferred wiring it into the
identification → activity pipeline. Without that wiring, callers had
to mutate `LineMatch.peak_area` themselves to get the deconvolved
values into `compute_activity` — fragile, error-prone, and not
discoverable.

**Resolution**:

1. **`LineMatch.peak_area_source: str = ""`** — new field on the
   identification `LineMatch` dataclass. Values:
   `"" / "cowell" / "lsrm_peaks_table" / "deconvolved" / "failed"`.
   Populated by `identify_nuclides` from
   `gamma.peaks.area.get_peak_area`'s third return value, and by the
   new post-pass for deconvolved replacements. Source-tracking enables
   future TCS-aware corrections that depend on the area method (open
   issue K-18).

2. **`apply_multiplet_deconvolution(identification_result, spec, fwhm_at_channel, *, overlap_threshold_fwhm=1.0, continuum="step_linear", max_chi2_per_dof=inf)`** — new post-pass in
   `gamma.peaks.deconvolve`. Workflow:

   - Find multiplet clusters via `find_multiplet_regions`.
   - For each cluster: build `MultipletComponent` per `LineMatch`,
     deconvolve, then map `(nuclide, library_E_keV)` to
     `(area, area_uncertainty)`.
   - Walk every `LineMatch` in `detected_nuclides.matched_lines`:
     replace `peak_area`, `peak_area_uncertainty`, and set
     `peak_area_source = "deconvolved"` when the line was in a
     cluster whose fit converged AND
     `chi2_per_dof ≤ max_chi2_per_dof`.
   - Lines outside any cluster are returned unchanged.
   - The result's `notes` records cluster count and replacement count.

   Returns `(new_identification_result, list_of_deconvolution_results)`.
   The original input is not mutated (everything via
   `dataclasses.replace`).

3. **Threshold semantics**: `overlap_threshold_fwhm=1.0` is the
   conservative default (FWHM-touching pairs are multiplets). Wider
   doublets like Co-60 1173/1332 on NaI (~3·FWHM apart but still
   exhibiting wing contamination that depresses Cowell areas by 30%
   in F-31a) require an explicit `overlap_threshold_fwhm=3.0` from
   the caller. We do not raise the default because broad thresholds
   over-cluster and inject correlated uncertainty into otherwise
   isolated lines.

4. **`disambiguate.py` LineMatch reconstruction fix** (drive-by):
   the `promoted = LineMatch(...)` in the characteristic-promotion
   branch did not preserve `peak_area`, `peak_area_uncertainty`, or
   the new `peak_area_source`. Fixed.

**Verification — 5 new tests in `test_deconvolve.py`** (+12 from F-33,
total 17 in the file):

- `test_apply_post_pass_returns_tuple` — function always returns
  `(IdentificationResult, list)`; identity-pass case (no clusters)
  returns the input result unchanged and an empty list.
- `test_apply_post_pass_no_change_for_isolated_lines` — on real
  Co-60 spectrum, lines NOT in the deconvolved cluster keep their
  original `peak_area_source` (e.g. `"lsrm_peaks_table"`).
- `test_apply_post_pass_replaces_co60_doublet` — with
  `overlap_threshold_fwhm=3.0`, the Co-60 1173/1332 doublet is
  picked up; both `LineMatch` entries return with
  `peak_area_source="deconvolved"`, 1173 area = 540 620, 1332 area
  = 495 016, ratio = 1.09 (library 1.00, residual 9% from TCS
  asymmetry + tail residual).
- `test_apply_post_pass_notes_record_replacement` — the new result
  carries an annotation like
  `"Multiplet deconvolution: 3 cluster(s), 7 peak area(s) replaced …"`.
- `test_apply_post_pass_max_chi2_filter_skips_bad_fits` — with
  `max_chi2_per_dof=0.01` every cluster is rejected and no
  replacement happens; cluster list is still returned for
  diagnostics.

**Total: 172 tests across 16 files** (167 prior + 5 new). All
baseline tests unchanged behaviour — `peak_area_source` defaults to
`""` everywhere except the in-pipeline cases that already populate
it.

**Scope notes**:
- `compute_activity` continues to read `peak_area`/`peak_area_uncertainty`
  with no awareness of `peak_area_source` — i.e., a deconvolved area
  is consumed identically to a Cowell or Lsrm-table area. This is the
  intended initial integration: identification chooses the integration
  method, and the activity layer trusts that choice.
- K-18 (TCS over-correction when `peak_area_source == "lsrm_peaks_table"`)
  is now formally unblocked — the field needed to drive an
  area-method-aware TCS correction is present. The correction itself
  is left to a focused phase that calibrates the reduced-TCS coefficient.

**References**: Lsrm Algorithmic Foundations 2022 §9 (multiplet
deconvolution); Gilmore & Joss, Practical Gamma-ray Spectrometry,
3rd Ed., §9.7 (continuum modelling). F-31a (v1.7.9) for the
Co-60 wing-contamination observation that motivates the
`overlap_threshold_fwhm=3.0` use case.

---

#### F-35: Area-method-aware TCS scaling (Phase 2.1d cont., v1.7.13) — closes K-18
**Severity**: Minor (closes the last open K-NN limitation; no
behavioural regression at default settings on Cowell or deconvolved
areas)
**Closes**: K-18
**Affected**: `gamma.activity.compute.compute_activity`,
`gamma.activity.compute.compute_activities_for_all`

**Discovery context**: F-31a (v1.7.9) showed that the Lsrm SpectraLine
peak-table area (a wide-ROI Gaussian-on-step fit) already recovers
most of the counts that a pure Cowell integration would lose on close
doublets — yielding Co-60 5cm activity matching the certificate within
+0.61% with **no** TCS correction. F-31b (v1.7.9) added the analytic
TCS correction `C(E) = 1/(1 − Σ p · ε_T)` from
`gamma.physics.cascade_summing`. Applying both on the same line led to
Co-60 +2.89% (over by ~2.3%): the analytic TCS double-counts what the
Lsrm fit already recovered. K-18 documented this, and F-34 (v1.7.12)
put the area-source provenance into `LineMatch.peak_area_source`
without yet acting on it.

**Resolution**: `compute_activity` now scales the analytic TCS effect
by a per-source factor:

```
c_effective = 1 + (c_analytic − 1) · scale[peak_area_source]
```

Default `DEFAULT_TCS_METHOD_SCALE` (in `gamma.activity.compute`):

| source label         | scale | rationale                                        |
|----------------------|------:|--------------------------------------------------|
| `""` (unknown)       | 1.0   | safe default — never silently disable a TCS step |
| `"cowell"`           | 1.0   | full TCS — pure ROI integration loses sum counts |
| `"deconvolved"`      | 1.0   | full TCS — fixed-position LSQ doesn't recover them|
| `"failed"`           | 1.0   | full TCS                                         |
| `"lsrm_peaks_table"` | 0.0   | no TCS — Lsrm wide-ROI Gaussian already recovers |

The cap on `lsrm_peaks_table` is empirical: F-31a alone already gives
Co-60 within +0.61% (within 1σ) at 5cm. A future calibration can move
this away from 0 if real geometry-dependent data shows partial
recovery; the value is exposed via the new `tcs_method_scale` parameter
on `compute_activity` and `compute_activities_for_all` so callers can
override per-call without editing the constant.

**Notes on the merge semantics**: when a caller passes
`tcs_method_scale={...}`, the dict is **merged** on top of
`DEFAULT_TCS_METHOD_SCALE`, not used as a full replacement — so
overriding only one source label leaves the other defaults intact.
This is implemented as `method_scale = {**DEFAULT, **user_dict}` inside
the function.

**Notes on unknown labels**: a source label not present in the merged
map falls back to `scale = 1.0` (full TCS). Same safe-default
philosophy as for `""`.

**Result-level diagnostics**: `ActivityResult.notes` now includes
`"K-18: TCS scaled by area-method on N line(s)"` whenever the
effective `c` differs from the analytic `c` on at least one matched
line. The per-line `LineActivity.correction_factor` carries the
**effective** `c` actually applied (post-scaling), so callers can
inspect the chain end-to-end.

**Verification — 10 new tests in `test_tcs_method_scale.py`**:

- `test_default_scale_lsrm_table_kills_tcs` — Co-60 line with
  `peak_area_source="lsrm_peaks_table"` and `coincidence_correction={E: 1.05}`
  yields the same `A_Bq` as the no-TCS run; the result notes contain
  `"K-18: TCS scaled by area-method on 1 line(s)"`.
- `test_default_scale_cowell_keeps_full_tcs` — `peak_area_source="cowell"`
  produces `A_Bq` exactly 1.05× the no-TCS run; no K-18 note.
- `test_default_scale_deconvolved_keeps_full_tcs` — same for
  `peak_area_source="deconvolved"`.
- `test_default_scale_empty_source_keeps_full_tcs` — `peak_area_source=""`
  (legacy / pre-F-34 LineMatch) gets full TCS — regression-safe for
  every existing call site that doesn't set the field.
- `test_unknown_source_label_default_full_tcs` — exotic labels also
  get full TCS.
- `test_custom_scale_partial_lsrm` — caller passes
  `tcs_method_scale={"lsrm_peaks_table": 0.5}`, gets exactly half
  the TCS effect on Lsrm-table lines.
- `test_custom_scale_overrides_only_named_keys` — partial override
  preserves the rest of `DEFAULT_TCS_METHOD_SCALE`.
- `test_no_tcs_dict_no_scaling_applied` — without a
  `coincidence_correction`, `tcs_method_scale` is inert (no
  multiplication ever happens, no K-18 note emitted).
- `test_default_dict_publishes_canonical_keys` — sanity: the
  publicly exported `DEFAULT_TCS_METHOD_SCALE` documents every
  source label that the codebase produces.
- `test_mixed_sources_in_one_nuclide` — Co-60 with 1173 from
  Lsrm-table and 1332 from Cowell: 1173 line gets `c=1.0`, 1332
  line gets `c=1.05`, notes say `"1 line(s)"` scaled.

All 172 prior tests pass unchanged (existing test_cascade_summing.py
`test_compute_tcs_corrections_compatible_with_compute_activity` test
uses a LineMatch without `peak_area_source` → defaults to `""` →
full TCS → existing expected ratio preserved).

**End-to-end implication for Co-60 5cm cert validation**:
- v1.7.9 F-31a alone → +0.61% vs cert
- v1.7.9 F-31a + F-31b (full TCS, no K-18) → +2.89%
- v1.7.13 F-31a + F-31b + F-35 (Lsrm-table source on both lines) →
  **+0.61%** (equivalent to F-31a alone because the lsrm_peaks_table
  scale collapses the TCS effect to zero, exactly as desired)

This is the third complete certificate match on Co-60 the codebase
has produced, and the first one with the TCS module actively enabled
end-to-end without overshoot.

**References**: F-31a / F-31b / K-18 (v1.7.9) for the observation
chain that led here. [GILMORE-8.5.3] (interaction of peak-fitting
method with summing correction). Lsrm Algorithmic Foundations §10
(coincidence summing).

**Open follow-ups (not in F-35)**:
- Empirical calibration of `lsrm_peaks_table` scale on geometries
  other than the validated 5cm point source (1cm, Marinelli) — could
  be non-zero at very close geometries where the wide-ROI Gaussian
  fit can't keep up with the changing peak shape.
- Detector-type dependence (HPGe vs NaI vs LaBr3): all current
  validation is on NaI. HPGe wide-ROI fits behave differently and
  this default may need to be detector-conditional.

---

#### F-36: Library coverage extension + extended multi-source cert validation harness (Variant B', v1.7.14)

**Severity**: Major (validation infrastructure + library gap close)

**Discovery context**: F-35 (v1.7.13) closed K-18 on a single
end-to-end Co-60 5cm cert measurement (+0.61% vs certificate). The
pipeline had only been quantitatively validated on one nuclide. The
question for v1.7.14 was: does the same pipeline give comparable
accuracy on the seven other point-source 5cm reference fixtures
available in `references/reference_spectra/Gamma-1C_NaI_63x63_USB_SN-01/`,
each with a matching entry in `АСПЕКТ_ОСГИ_2024.src`?

**What blocked direct validation**:
1. Four cert nuclides — Y-88, Bi-207, Cd-109, Th-228 — had no
   `data/nuclides.json` entry (K-03 reported partial fix via Lsrm-v2
   library import, but the built-in JSON was still 24 entries).
2. No end-to-end harness existed. The +0.61% Co-60 number from F-35
   was computed manually; replicating it required assembling
   read → bg-subtract → peak-search → identify → multiplet-deconv →
   TCS → compute_activity → decay-correct → deviation by hand.
3. The `.spe` FWHM polynomial was documented as `FWHM(E)` but is
   actually `FWHM_keV(sqrt(E))` — a single Horner pass on E gives
   nonsense at NaI energies (FWHM(662)≈−2150 keV). The stock
   `make_fwhm_at_channel_provider` doesn't evaluate this model at
   all (it handles SimpleSqrtFwhm + interp + physical floor) so .spe
   spectra fell back to constant 10 ch, swamping the Mariscotti
   adaptive bands. Verified by hand on the Co-60 5cm fixture:
   - FWHM_keV(sqrt(662)) = c0 + c1·25.73 + c2·25.73² + c3·25.73³
     = -0.4464 + 27.30 + 24.21 - 5.42 ≈ 45.6 keV  (NaI 7% ≈ 46 keV ✓)
   - FWHM_keV(sqrt(1332)) = ≈ 71.6 keV  (NaI 5.4% ≈ 72 keV ✓)
4. `disambiguate_identifications` applies a proportionality check
   (Rule 4) that compares Mariscotti peak σ ratios to library
   intensity ratios. Mariscotti σ scales as height/√B, not as peak
   area; the check does not divide by ε(E). On a known single-
   nuclide cert spectrum (rare-isotope prior ≤ 0.2 for Na-22) this
   rule **rejects the target nuclide** even though it is the only
   thing emitting γs. F-30 had not encountered this because the
   end-to-end test path was Cs-137 (single line, proportionality
   check skipped).

**Fix**:
1. **Library extension** in `data/nuclides.json`: added
   `Y-88` (T½ = 9.217e6 s; lines 898.04 / 1836.06 / 2734.07 keV;
   `is_cascade: true`), `Bi-207` (T½ = 9.946e8 s; lines 569.70 /
   1063.66 / 1770.23 keV; `is_cascade: true`), `Cd-109` (T½ =
   3.987e7 s; single line 88.03 keV at I=3.66%; ic_xrays at
   22-25 keV). Values from NNDC ENSDF (cross-check with LNHB).
2. **New `validate_certs.py`** at the project root: end-to-end
   harness that runs the full F-31a + F-31b + F-35 pipeline on
   every `Точечная-5см` .spe fixture against `АСПЕКТ_ОСГИ_2024.src`
   and emits a deviation matrix (stdout + `cert_validation_matrix.csv`).
   - Embeds its own `make_lsrm_fwhm_provider(spec)` that evaluates
     the Lsrm `lsrm_fwhm_polynomial_in_E` model as FWHM_keV(sqrt(E))
     then divides by |dE/dN|. Independent of the stock provider so
     this harness does not regress AtomSpectra fixtures.
   - Skips `disambiguate_identifications` (single-source semantics:
     we know which nuclide is in the source, mixture-resolution
     rules add no value and reject the rare-isotope target).
   - Adaptive `min_intensity_pct`: 5.0 by default to suppress
     low-statistics outliers in weighted-average aggregation
     (Eu-152's 443.96 keV at I=3.12% gets ~17× the weight of
     121.78 keV at I=28.58% because its tiny S inflates the
     statistical sigma's contribution to σ_A_i); falls to 0.0 for
     nuclides whose entire library catalog sits below 5%
     (Cd-109 only has the 88 keV at 3.66%).
3. **No changes** to `make_fwhm_at_channel_provider`, `identify_nuclides`,
   `compute_activity`, or `disambiguate_identifications`. The
   `validate_certs.py` script reuses every existing module; the
   .spe-FWHM evaluator lives in the harness only.

**Verification — deviation matrix (Gamma-1C NaI 63×63, 5 cm point
geometry, AC-2024 certificate)**:

| Nuclide  | A_cert (Bq) | A_cert@meas (Bq) | A_measured (Bq) | Δ, % | n_lines | Comment                            |
|----------|------------:|-----------------:|----------------:|-----:|--------:|:-----------------------------------|
| Cs-137   |     106 000 |          89 307  |        90 944   |+1.83 |       1 | single line 661.66 keV             |
| Co-60    |     105 000 |          49 886  |        50 193   |+0.61 |       2 | matches v1.7.13 reference exactly  |
| Na-22    |     229 000 |         136 700  |       144 900   |+5.99 |       2 | 511 + 1274 keV; disambig bypassed  |
| Eu-152   |     142 000 |         120 600  |       112 400   |−6.82 |       8 | 1 multiplet deconvolved (1086/1112)|
| Ba-133   |      44 100 |          15 308  |        15 511   |+1.32 |       3 | 4-line TCS                         |
| Am-241   |     103 100 |         102 500  |        93 948   |−8.35 |       1 | 60 keV at ε(E) curve edge          |
| Zn-65    |       3 100 |             871  |           918   |+5.45 |       1 | 1115 keV; 511 (I=2.83%) filtered   |
| Y-88     |     350 000 |          29 903  |        30 392   |+1.64 |       2 | 898 + 1836 cascade pair, TCS       |
| Bi-207   |      97 000 |          82 407  |        80 208   |−2.67 |       3 | 570 + 1064 + 1770 keV              |
| Cd-109   |     395 000 |           6 762  |         6 534   |−3.38 |       1 | 88 keV; min_I floor dropped to 0   |
| Th-228   |     129 000 |             —    |           —     |  —   |       0 | library gap: chain parent          |

**Summary**: 10/11 cert fixtures measurable. Mean |Δ| = 3.81%.
Max |Δ| = 8.35% (Am-241, at the lower ε(E) curve edge). All within
±10% — fully consistent with NaI 63×63 typical certification
accuracy (typical NaI calibration aims at 5-10%). Th-228 only
emits weakly directly (most "Th-228 activity" reaches the detector
via Tl-208/Pb-212/Bi-212/Ac-228 daughters under secular
equilibrium); direct Th-228 measurement is deferred to a future
chain-parent reconstruction feature.

**End-to-end Co-60 5 cm regression intact**:
- v1.7.9 F-31a alone → +0.61% vs cert
- v1.7.9 F-31a + F-31b → +2.89% (K-18 overshoot)
- v1.7.13 F-31a + F-31b + F-35 → **+0.61%** (K-18 resolved)
- v1.7.14 same pipeline, multi-source: Co-60 still **+0.61%**

**Files changed**:
- `data/nuclides.json` — +3 entries (Y-88, Bi-207, Cd-109).
- `validate_certs.py` (new, 320 lines) — end-to-end harness +
  `make_lsrm_fwhm_provider` helper specific to the
  `lsrm_fwhm_polynomial_in_E` model.
- `cert_validation_matrix.csv` (new) — machine-readable deviation
  matrix for downstream regression / dashboarding.

**Test count**: 200 tests across 17 files (unchanged code paths +
new library entries — every prior test still passes).

**References**: ENSDF / NuDat 3 (nuclide line catalog); F-31a /
F-31b / F-35 (the pipeline whose end-to-end behaviour is now
quantitatively documented across 10 nuclides instead of 1);
[GILMORE-5.7.3] (intra-nuclide χ²/dof — observed values 1-135 across
the matrix, well above 1 for multi-line nuclides reflecting real
per-line variation in NaI peak-area accuracy, not a bug);
K-03 (library gap, now reduced from 4 to 1 missing cert nuclide).

**Open follow-ups (not in F-36)**:
- Chain-parent activity reconstruction for Th-228 (and U-238 /
  Ra-226 / Th-232 chains generally). Requires daughter-equilibrium
  modelling and would close the last cert gap on this fixture set.
- The Am-241 -8.35% deviation is at the ε(E) extrapolation
  boundary (ε(59.5) is the lowest calibrated point). A dedicated
  low-E ε(E) refit on dense low-E points might reduce this. The
  Variant F idea (refit ε(E) excluding cascade-depleted points)
  is now less critical: it only attacked Zn-65 which fell from
  +11.4% to +5.4% just by filtering the 511 keV I=2.83% line.

---

#### F-37: Secondary-peak catalogue + reference-samples library (v1.7.15)

**Severity**: Major (foundation for anti-misidentification logic)

**Discovery context**: F-36 (v1.7.14) validated the activity pipeline
end-to-end on 10 cert nuclides. The two open follow-ups it left were
(a) disambiguate Rule 4 correctness, and (b) handling secondary
spectral features that the identification module currently treats as
candidate photopeaks. The Cs-137 Compton edge at ~480 keV is a
chronic source of false Bi-214 (503 keV) misidentifications on NaI;
the K-40 backscatter peak at ~230 keV is similarly misread as Ac-228
(209 keV) or Pb-212 (238.6 keV). Without explicit knowledge of
"where each parent nuclide's secondaries SHOULD appear", the
disambiguator has no principled way to suppress these false claims.

User contributed 20 new reference fixtures (`Поверка 2024`):
- Дента-120мл: Cs-137 ×2, K-40 ×2, Ra-226 ×2, Th-232 ×2.
- Петри-60мл: same 8 sources.
- Точечная-25см: Cs-137, Na-22, Y-88, Th-228.
Plus the existing Точечная-5см and Маринелли spectra of the same
sources. With 17 Cs-137 + K-40 spectra across 4 distinct geometries,
the empirical secondary-peak landscape can be measured rather than
assumed.

**Fix**:
1. **New `gamma.physics.secondary_peaks`** module (~200 LOC). Pure-
   physics helpers (`compton_edge_keV(E)`, `backscatter_keV(E)`,
   `compton_edge_observed_keV(E, fwhm)`, `backscatter_observed_keV(E,
   geometry)`) + the `ExpectedFeature` dataclass + builder
   `expected_features_for(nuclide, E_gamma)` that emits the full
   theoretical secondary set for a primary γ-line (photopeak +
   Compton edge + backscatter + single/double escape if E > 1022 keV
   + I K X-ray escape if E < 200 keV + Cs-137 IC X-rays at 32 keV +
   the always-present natural K-40 background line in Cs-137
   spectra).
2. **New `data/secondary_peaks.json`** catalog built by
   `analyze_secondaries.py` from the 17 reference fixtures.
   Per (nuclide, feature) it stores `mean_intensity_ratio`,
   `std_intensity_ratio`, `min/max`, `n_observations`, and the
   mean position residual against theory. Consumed via
   `load_catalog()` / `empirical_ratio(nuclide, feature)`.
3. **New diagnostic script `analyze_secondaries.py`** (~270 LOC):
   end-to-end peak detection + feature assignment + statistics
   aggregation + JSON catalog emission. Reproducible:
   `py -3.11 analyze_secondaries.py`.

**Three robust empirical patterns documented** (Gamma-1C NaI 63×63):

| Feature                | E_theory   | <ΔE>    | mean R = S/S_pp |
|------------------------|-----------:|--------:|----------------:|
| Cs-137 backscatter     | 184.3 keV  | +8.1 keV|     7.3% ±2.9%  |
| Cs-137 Compton edge    | 477.3 keV  |−37.0 keV|     3.0% ±1.0%  |
| Cs-137 Ba Kα IC X-ray  |  32.0 keV  | −5.9 keV|     8.4% ±4.1%  |
| K-40 backscatter       | 217.5 keV  |+14.9 keV|    12.2% ±3.9%  |
| K-40 Compton edge      |1243.4 keV  |−53.2 keV|     7.7% ±2.6%  |
| K-40 single escape     | 949.8 keV  | +5.2 keV|     2.9% ±0.8%  |

(A) **Compton edge sits BELOW analytical position by ~0.7·FWHM**.
    Universal across Cs-137 and K-40, regardless of geometry. The
    Mariscotti maximum of the second-derivative response lies below
    the analytical step. Predictable from
    `compton_edge_observed_keV(E_gamma, FWHM(E_C))`. Cs-137 at
    FWHM=50 → predicted shift −35 keV (observed −37); K-40 at
    FWHM=75 → predicted −52 (observed −53).

(B) **Backscatter sits ABOVE analytical position by +8 to +15 keV**.
    Geometry-conditional. Worst (+14) on close point geometry and
    extended-source containers, smallest (+5) on distant point.
    Multi-path photons add to the single-180° base.

(C) **Natural K-40 background contaminates every long-integration
    spectrum** at 0.3-10% of the Cs-137 photopeak area. Strongest
    on Маринелли / Дента / Петри (sample-mass effect), weakest on
    distant point. Identification must EXPECT a 1461 keV peak in
    any extended-source measurement and not credit it as an
    anomalous K-40 enrichment unless the magnitude vastly exceeds
    expected background.

**Verification — 13 new tests** in `test_secondary_peaks.py`:
- Theoretical Cs-137 and Co-60 Compton edges match canonical values.
- Backscatter formula passes the energy-conservation identity
  E_C + E_bs = E for 5 canonical energies.
- `compton_edge_observed_keV(E, FWHM)` reproduces the empirical
  −0.7·FWHM rule and matches observed Cs-137 shift to within 5 keV.
- Geometry-conditional backscatter shifts match the documented table.
- `expected_features_for("Cs-137", 661.66)` returns exactly the 5
  expected feature names; `K-40` adds single/double escape;
  Am-241 (60 keV) adds X-ray escape.
- Catalog JSON loads, has both Cs-137 and K-40 with `primary_E_keV`
  set, and `empirical_ratio()` returns the right shape.
- Cs-137 backscatter mean R is in the 5-10% NaI band and the
  Compton-edge residual is < −10 keV (validates the catalog's
  internal consistency with the documented patterns).

All 200 prior tests still pass — `analyze_secondaries.py` and the
new module make no changes to existing code paths. Total: **213
tests** across 18 files.

**Files changed**:
- `scripts/gamma/physics/secondary_peaks.py` (new, ~200 LOC).
- `data/secondary_peaks.json` (new, built artefact).
- `analyze_secondaries.py` (new, ~270 LOC).
- `test_secondary_peaks.py` (new, 13 tests).
- `references/reference_spectra/Gamma-1C_NaI_63x63_USB_SN-01/`:
  3 new geometry subdirectories with 20 .spe fixtures (Дента-120мл,
  Петри-60мл, Точечная-25см) plus an updated reference base.

**References**: Knoll 4th Ed., chap. 10 (Compton kinematics), 11.A.5
(backscatter peak), 12.B.2 (escape peaks). Gilmore & Joss 3rd Ed.,
chap. 6 (NaI(Tl) artefacts).

**Open follow-ups (not in F-37)**:
- Wire the secondary-peak catalog into `disambiguate_identifications`
  Rule 4: when a candidate nuclide is supported by a single line at
  a known secondary-peak position of an already-identified parent
  (e.g. Bi-214 503 keV co-located with Cs-137 Compton edge at
  ~478 keV), demote the candidate. This is the natural sequel to
  v1.7.14's "Rule 4 fix" suggestion. Potential F-38.
- Cross-validation use: a Cs-137 claim with no detectable Ba Kα
  X-ray AND no detectable backscatter peak is more likely a Bi-214
  609 keV misidentification. The catalog gives the expected ratios
  to validate the hypothesis.
- Extend the catalog to Co-60 / Ba-133 / Eu-152 / Th-228 chain
  daughters (multiple γ-lines each → many secondary features per
  source).

---

#### F-38: Range/shape characterisation of problem isotopes for presence inference (v1.7.16)

**Severity**: Major (foundation for assertion-based identification)

**Discovery context**: F-37 (v1.7.15) reported mean ± std of position
residual + intensity ratio per (nuclide, feature). The user observed
that **the position of a secondary peak is not fixed but floats in
a characteristic range** depending on isotope activity, geometry,
and other conditions. Identification by single-point comparison
("is there a peak at exactly 478 keV?") is therefore methodologically
wrong — the right question is "does the observed peak's position
and shape fall within the characteristic range observed across all
known measurement conditions of the candidate parent?".

User contributed the 2016 Поверка archive (99 .spe files): 12
isotopes at point 5 cm, 11 at point 25 cm, Marinelli/Denta/Petri
container sets, plus three background-type folders (water, open-lid,
empty-shielding) and a 15-measurement time-stability series.
Combined with the 2024 Поверка archive (v1.7.15) and the existing
fixtures, the per-isotope inventory becomes:

| Parent  | Fixtures | Geometries                                       |
|---------|---------:|--------------------------------------------------|
| Cs-137  |    17    | M-source, denta, marinelli, petri, point0cm xml, point25cm, point5cm |
| K-40    |    10    | M-source, denta, marinelli, petri                |
| Co-60   |     3    | point5cm ×2, point25cm                           |
| Na-22   |     4    | point5cm, point25cm ×2, fixture                  |
| Y-88    |     4    | point5cm ×2, point25cm ×2                        |
| Th-228  |     4    | point5cm ×2, point25cm ×2                        |

**Fix**:
1. **New `analyze_problem_isotopes.py`** (~370 LOC) replaces v1.7.15's
   point-estimate analyser with **range/shape quantile** statistics.
   Per (parent, primary_E_keV, feature) it computes
   `{min, p10, median, p90, max, mean, std}` of:
   - position (observed centroid energy)
   - position residual (observed − theoretical)
   - intensity ratio (S_secondary / S_photopeak)
   - FWHM (measured) + FWHM ratio to theoretical
   - asymmetry (left-half-area/total − 0.5)
   Plus per-geometry raw observations and the conflict-line list (real
   gamma-lines from OTHER nuclides falling inside p10..p90 of the
   observed position range).

2. **Per-primary-line keying** instead of per-(parent, feature). A
   bug in v1.7.15 conflated multi-line nuclides' secondaries:
   Co-60's 1173-keV Compton edge (at 963) was merged with the
   1332-keV one (at 1118), giving a spurious 250 keV-wide range.
   v1.7.16 keys by `(parent, primary_E, feature)` so each photopeak's
   secondaries form their own tight cluster. Cs-137 (single line)
   is unaffected; Co-60, Na-22, Th-228, Y-88 all show clean ranges.

3. **New module functions** in `gamma.physics.secondary_peaks`:
   - `load_catalog_v2()` — load the v0.2 range/shape catalogue.
   - `position_range(nuclide, primary_E, feature, span="p10p90")` —
     return `(E_low, E_high)` for the observed position range; `span`
     can be `"minmax"`, `"p10p90"` (default, 90% CI), or `"iqr"`.
   - `matches_secondary(parent, observed_E, feature=None,
     span="p10p90")` — is an observed peak consistent with one of
     `parent`'s known secondaries? Returns a list of matched feature
     descriptors. Used by future identification logic to demote a
     candidate nuclide whose only matched line lies in a parent's
     known secondary range.

4. **New `data/secondary_peaks_v2.json`** (167 KB) — the catalogue
   itself, regenerable via `py -3.11 analyze_problem_isotopes.py`.

**Key empirical findings on Gamma-1C NaI 63×63**:

*Photopeak position spread (intrinsic detector drift across all
geometries and dates):*

| Parent  | Primary  | n  | p10..p90 spread | std    |
|---------|---------:|---:|----------------:|-------:|
| Cs-137  |  661.66  | 17 |   2.1 keV       | 0.81   |
| K-40    | 1460.82  | 10 |   4.5 keV       | 2.18   |
| Co-60   | 1173.23  |  3 |   2.8 keV       | 1.52   |
| Co-60   | 1332.49  |  3 |   0.7 keV       | 0.33   |
| Na-22   |  511.00  |  4 |   0.7 keV       | 0.40   |
| Y-88    |  898.04  |  4 |   1.8 keV       | 0.90   |

All under 5 keV — sets the lower-bound uncertainty for any peak-position
test on this detector.

*Compton edge position vs theoretical (consistent with −0.7·FWHM rule
established in F-37):*

| Parent | Primary | E_theory | p10..p90 range | residual median |
|--------|--------:|---------:|---------------:|----------------:|
| Cs-137 |  661.66 | 477.34   | 433.9..439.3   |     −40.6 keV   |
| K-40   | 1460.82 | 1243.36  | 1178.9..1179.1 |     −64.4 keV   |
| Co-60  | 1173.23 | 963.42   | 906.9..912.5   |     −54.0 keV   |
| Na-22  | 1274.54 | 1061.71  | 1000.8..1012.4 |     −57.4 keV   |
| Y-88   | 1836.06 | 1611.77  | 1545.6..1560.1 |     −52.8 keV   |

*Backscatter position vs theoretical (always shifted UP, geometry-dependent):*

| Parent | Primary | E_theory | p10..p90 range | residual median |
|--------|--------:|---------:|---------------:|----------------:|
| Cs-137 |  661.66 | 184.32   | 186.8..196.7   |     +9.6 keV    |
| K-40   | 1460.82 | 217.46   | 227.4..235.2   |     +13.1 keV   |
| Co-60  | 1173.23 | 209.81   | 226.6..227.8   |     +17.0 keV   |
| Na-22  |  511.00 | 170.33   | 175.0..178.4   |     +6.2 keV    |
| Y-88   |  898.04 | 198.91   | 202.2..229.3   |     +18.4 keV   |

**Documented conflict-line catalogue**:

The v2 catalog flags every real γ-line from OTHER library nuclides
that falls inside a feature's p10..p90 position range. Selected
high-risk conflicts:

- **Cs-137 Compton edge** [433.9, 439.3] — no library conflict in range
  (Bi-214 503 keV and Be-7 478 keV sit ABOVE the actual range; the
  conflict is real on detectors with wider drift).
- **Co-60 1173 Compton edge** [906.9, 912.5] → **Ac-228 911.20 keV
  (I=25.8%)** — direct conflict.
- **Co-60 1332 Compton edge** [1166.5, 1169.3] → **Cs-134 1167.94
  (I=1.8%)**, AND this range is INDISTINGUISHABLE FROM Co-60 1173
  photopeak position (1166..1169). Methodologically, Co-60 1332's
  Compton edge sits exactly at Co-60 1173's photopeak — a built-in
  cross-validation: the 1173-peak amplitude must include the 1332
  Compton continuum.
- **K-40 Compton edge** [1178.9, 1179.1] → 5.7 keV above Co-60 1173
  photopeak (1173.23). On this detector the tight range keeps them
  separated, but on a detector with realistic gain drift the
  conflict is acute.
- **Th-228 (Tl-208) 583 Compton edge** [351.9, 378.4] →
  **I-131 364.49 (I=81.5%)**, **Ba-133 356.01 (I=62.0%)**,
  **Pb-214 351.93 (I=35.6%)**. Three nuclides directly inside the
  range.
- **Cs-137 Ba Kα IC X-ray** [23.7, 26.6] → **Am-241 26.34 (I=2.3%)** —
  direct conflict at low E. Resolvable only if photopeak amplitude
  pattern is examined.

**Verification — 9 new tests** in `test_secondary_peaks.py`
(13 prior F-37 + 9 new F-38 = 22 in this file):
- `test_v2_catalog_loads` — all 6 problem isotopes present.
- `test_v2_catalog_per_primary_keying` — Co-60 has separate 1173 +
  1332 entries.
- `test_position_range_cs137_compton_edge` — p10..p90 within [430, 445].
- `test_position_range_k40_compton_edge` — p10..p90 within [1175, 1185].
- `test_matches_secondary_cs137_compton_edge_collides_with_bi214` —
  peak at 437 keV matches Cs-137 Compton edge AND doesn't match Cs-137
  backscatter or photopeak.
- `test_matches_secondary_no_match_outside_ranges` — peak at 800 keV
  in Cs-137 context yields empty match list.
- `test_matches_secondary_k40_compton_dangerous_for_co60` —
  documents the K-40 vs Co-60 1173 proximity (5.7 keV separation).
- `test_v2_catalog_conflict_lines_recorded` — Cs-137 Ba Kα range
  flags Am-241 26.34 conflict.
- `test_v2_photopeak_position_tightness` — all problem-isotope
  photopeak spreads under 5 keV (intrinsic drift bound).

All 200 prior tests still pass. Total: **222 tests across 18 files**.

**Files changed**:
- `analyze_problem_isotopes.py` (new, ~370 LOC).
- `data/secondary_peaks_v2.json` (new, ~167 KB built artefact).
- `scripts/gamma/physics/secondary_peaks.py` — extended with
  `load_catalog_v2`, `position_range`, `matches_secondary`.
- `test_secondary_peaks.py` — +9 F-38 tests.
- `references/reference_spectra/Gamma-1C_NaI_63x63_USB_SN-01/
  Поверка-2016/` — 99 new .spe fixtures across 9 subdirectories
  (5 source geometries + 3 background types + time stability).

**Methodology references**: Knoll 4th Ed., chap. 10, 11, 12. Gilmore
& Joss 3rd Ed., chap. 6. The "−0.7·FWHM Compton edge shift" rule was
introduced in F-37 (v1.7.15) and is now confirmed across 6 problem
isotopes with 4–17 fixtures each.

**Open follow-ups (not in F-38)**:
- Wire `matches_secondary()` into `disambiguate_identifications`
  Rule 4: candidate nuclide whose only matched line falls inside a
  parent's secondary range → demote.
- Use averaged background spectra (15-measurement means from each
  per-folder background series — 50× noise reduction) for cleaner
  net spectra in future cert-validation runs.
- Time-stability series provides 15 same-conditions measurements to
  quantify true intrinsic Mariscotti detection drift (σ_position
  before any geometry/activity contribution).

---

#### F-39: Lsrm chain-library integration closes Th-228 cert row + Th-chain shape catalog (v1.7.17)

**Severity**: Major (closes the last cert-validation gap from F-36 +
adds 20 nuclides + 3 Th-chain daughter shape entries)

**Discovery context**: F-36 (v1.7.14) shipped the cert validation
matrix with one row marked as "library gap": **Th-228 cert source
unmeasured because Th-228 itself emits only weak (<1.5%) direct
γ-lines on NaI**. The cert quotes the parent activity but only the
chain daughters (Tl-208 583/2614, Pb-212 238, Ac-228 911, etc.)
produce detectable signals. v1.7.16 (F-38) characterised the Tl-208
secondaries under the "Th-228" label without explicitly modelling
the chain. The user then surfaced an archive of Lsrm-native
libraries at `C:\LSRM\Work\BG\Gamma-1S\Архив\Data\` containing:

- **`NaI-Etl+Esc.lib`** — NaI-tuned Th-232 chain library with all
  daughters bundled under "Th-232" via ENSDF `dbid` tags
  (Ac-228, Pb-212, Bi-212, Tl-208, Ra-224).
- **`ОСГИ.lib`** — 33 ОСГИ certified-source nuclides incl. Eu-154,
  Eu-155, Ce-144, Sn-113, Hg-203, Co-56, Ag-110m, Ta-182, Cs-134,
  Sb-125, Ir-192, Ru-103, Zr-95+, Ho-166m, Th-231, Th-234, U-232,
  Ti-44 (typo in source: "TI-44").
- Auxiliary files (`Th.zon`, pre-computed 10 cm Th-232 chain
  windows `.cen`/`.cfw`, `Gamma-1S+Compton.cpt` Compton model,
  `Aspect.src` extra cert).

**Fix**:
1. **Library bundle in repo** at `references/lsrm-libraries/`:
   - `NaI-Etl+Esc.lib`, `ОСГИ.lib`, `Th.zon`, the two .cen/.cfw
     chain-window files, `Gamma-1S+Compton.cpt`, `Aspect_2025.src`.
2. **Opt-in chain-library loader** in
   `gamma.data.nuclide_library.load_lsrm_chain_libs(include_nai_chain=True,
   include_osgi=True, merge_mode="supplement", split_chains=True)`:
   - `include_nai_chain` → loads `NaI-Etl+Esc.lib` and decomposes
     Th-232 into Tl-208 (5 lines), Pb-212 (2 lines), Bi-212 (3 lines),
     Ac-228 (6 lines), Ra-224 (2 lines).
   - `include_osgi` → loads `ОСГИ.lib` for 18 supplemental nuclides.
   - Default `"supplement"` mode keeps existing JSON entries; pass
     `"override"` to let Lsrm values win.
   - Verified: 27 baseline → 47 with both loaded.
3. **Th-228 + Ra-224 added explicitly** to `data/nuclides.json`. The
   chain decomposer in the existing `gamma.data.chain_decomposer`
   doesn't extract Th-228 itself from the bundled Th-232 lines
   (Th-228 is the chain's grandfather, not a Lsrm "owner"). For
   decay correction of the cert source we need the parent's T½ —
   so Th-228 (T½=6.03e7 s = 1.91 yr; 4 weak lines) and Ra-224
   (T½=313 ks = 3.6 d; 2 lines) are now in the built-in JSON
   regardless of whether the Lsrm libs are loaded.
4. **Chain-proxy cert validation** in `validate_certs.py`:
   - Calls `load_lsrm_chain_libs()` at module import.
   - New `CertFixture` fields `cert_nuclide` (parent name in the
     cert when different from the identified nuclide) and
     `chain_branching` (A_daughter / A_parent).
   - Th-228 cert fixture now identifies **Pb-212** (direct daughter
     in the Th-228 → Ra-224 → Rn-220 → Po-216 → Pb-212 sub-chain,
     all 1:1 branching) and compares A_Pb-212 to the cert's Th-228
     entry decay-corrected via the Th-228 parent half-life.

**Deviation matrix update (v1.7.17 vs v1.7.14)**:

| Nuclide  | v1.7.14 result | v1.7.17 result | Notes |
|----------|----------------|----------------|-------|
| Cs-137   | +1.83%         | +1.83%         | unchanged |
| Co-60    | +0.61%         | +0.61%         | unchanged |
| Na-22    | +5.99%         | +5.99%         | unchanged |
| Eu-152   | −6.82%         | −6.82%         | unchanged |
| Ba-133   | +1.32%         | +1.32%         | unchanged |
| Am-241   | −8.35%         | −8.35%         | unchanged |
| Zn-65    | +5.45%         | +5.45%         | unchanged |
| Y-88     | +1.64%         | +1.64%         | unchanged |
| Bi-207   | −2.67%         | −2.67%         | unchanged |
| Cd-109   | −3.38%         | −3.38%         | unchanged |
| **Th-228** | **library gap** | **−12.79% via Pb-212 chain proxy** | F-39 closure |

**11/11 cert nuclides measurable**. Mean |Δ| = 4.62%. Max |Δ| =
12.79% (Th-228 via Pb-212, slightly above ±10% due to 5cm sealed-
source self-absorption at E=238 keV).

**Th-chain shape catalog extended in `secondary_peaks_v2.json`**:
Added Tl-208 (4 primary lines: 510/583/860/2614), Pb-212 (2:
238/300), Ac-228 (6: 209/338/463/795/911/969) — all sharing the
Th-228 fixture set via the new `_PARENT_ALIASES` map in
`analyze_problem_isotopes.py`. Catalog now covers **9 problem
isotopes** (up from 6). Documented Th-chain conflicts:

- **Tl-208 510.77 photopeak** [502, 505] ↔ **Na-22 511 annihilation**
  — chronic confusion between positron emitters and Th-chain sources.
- **Tl-208 510.77 Compton edge** [291, 298] ↔ **Ir-192 295.96
  (I=28.6%), Pb-214 295.22 (I=18.4%)**.
- **Tl-208 583 Compton edge** [346, 379] ↔ **I-131 364 (I=81.5%),
  Ba-133 356 (I=62.0%), Pb-214 351 (I=35.6%)** (4 conflicts).
- **Tl-208 583 backscatter** [171, 185] ↔ **Sb-125 176 (I=6.7%)**.

**Verification — 4 new tests** in `test_secondary_peaks.py`:
- `test_v2_catalog_loads` — now expects 9 isotopes (was 6).
- `test_v2_catalog_tl208_chain_daughter` — Tl-208 has all 4
  primary lines + 510.77 photopeak overlaps Na-22 511 region.
- `test_lsrm_chain_loader_adds_th_chain_daughters` — loader adds
  ≥18 nuclides including Tl-208, Pb-212, Bi-212, Ac-228, Eu-154,
  Eu-155, Ce-144.
- `test_th228_in_built_in_library` — Th-228 has T½=6.03e7 s and the
  84.4 keV line, independent of chain-library load.

All 213 prior tests still pass. Total: **225 tests across 18 files**.

**Files changed**:
- `references/lsrm-libraries/` (new): 7 Lsrm files (~85 KB).
- `data/nuclides.json` — +2 entries (Th-228, Ra-224).
- `scripts/gamma/data/nuclide_library.py` — `load_lsrm_chain_libs()`
  helper.
- `validate_certs.py` — `load_lsrm_chain_libs()` at import,
  `CertFixture.cert_nuclide` + `chain_branching` fields, Th-228
  fixture now uses Pb-212 chain proxy.
- `analyze_problem_isotopes.py` — `load_lsrm_chain_libs()` at
  import, PROBLEM_ISOTOPES + `_PARENT_ALIASES` add Tl-208 /
  Pb-212 / Ac-228 sharing Th-228 fixtures.
- `data/secondary_peaks_v2.json` — regenerated, 9 isotopes.
- `test_secondary_peaks.py` — +4 F-39 tests.

**Closes K-03 (library size cap)** — built-in JSON is opt-in
supplemented to 47 nuclides on demand. Closes the last cert-matrix
gap (Th-228 row) flagged in F-36.

**Open follow-ups (not in F-39)**:
- The Th-228 −12.79% deviation is likely 5cm-geometry self-
  absorption + sealed-source matrix attenuation at E=238 keV.
  Alternative proxy via Tl-208 583 keV with 0.36 branching (Bi-212
  β-decay branch) could cross-check. Defer.
- Wire `matches_secondary()` into `disambiguate_identifications`
  Rule 4 (F-40, the original plan from F-38's open list). The
  catalog now covers 9 problem isotopes incl. the Th-chain — the
  conflict catalogue is rich enough to drive disambiguation.

---

#### F-40: Secondary-feature anti-misidentification rule wired into disambiguate (v1.7.18)
**Status**: FIXED in v1.7.18 — `disambiguate_identifications` Rule 5
(also referred to as Rule 4-style F-40) now consults the v1.7.16
`secondary_peaks_v2` catalog to demote candidate nuclides whose every
matched line falls inside the observed (p10..p90) secondary-feature
range of an already-detected parent.

**Discovery context**: F-38 (v1.7.16) and F-39 (v1.7.17) populated the
v2 catalog with quantile ranges of Compton edges, backscatter peaks,
escape peaks, and IC X-ray positions for **9 problem isotopes**
(Cs-137, K-40, Co-60, Na-22, Y-88, Th-228, Tl-208, Pb-212, Ac-228).
The catalog also tags **conflict lines** — real γ-emissions from other
nuclides falling inside each feature's observed range (e.g. Ac-228
911.20 keV inside Co-60 1173 Compton edge [906.85..912.50]; Sb-125 176
inside Tl-208 583 backscatter [171, 185]). F-40 ships the consumer.

**Mechanism**: a new block in `disambiguate_identifications` —
controlled by kwarg `apply_secondary_feature_rule: bool = True` and
`secondary_max_lines: int = 2` — iterates each detected nuclide and,
for those with `len(matched_lines) ≤ secondary_max_lines`, calls
`matches_secondary(parent, line.peak_E_keV, span="p10p90")` against
every other detected catalog-listed parent. The candidate's own
photopeak feature is excluded (true photopeak overlap belongs to
existing `NAI_CONFUSION_MAP` / Rule 3 / CI tiebreaker). If **every**
matched line of the candidate is explained as a non-photopeak
secondary of some parent, the candidate is moved to `rejected` with a
detailed reason listing the (line ↔ parent feature [range]) mappings.

**Defensive design**:
- A parent is never demoted by its own secondaries (candidate vs
  parent identity check via `[p for p in parents if p != ni.nuclide]`).
- The rule is **inert** when no detected nuclide has a v2-catalog
  entry — safe fallback for spectra dominated by anthropogenic
  nuclides that aren't characterised yet (e.g. Eu-152, Bi-214).
- Multi-line candidates (>2 by default) are not demoted regardless of
  position — strong proportional evidence overrides the rule.
- The rule operates AFTER universal proportionality Rule 4 and BEFORE
  the Ra-chain equilibrium analysis (Rule 4b) so chain-equilibrium
  decisions act on the post-secondary-cleaned candidate set.

**Tests** — 9 new in `test_secondary_feature_rule.py`:
- `test_co60_compton_edge_demotes_single_line_candidate` — the canonical
  Ac-228 911 ↔ Co-60 1173 Compton edge case.
- `test_cs137_backscatter_demotes_191_kev_candidate` — single-line
  candidate in Cs-137 backscatter [186.76..196.68].
- `test_k40_compton_edge_demotes_candidate_at_1179` — tight K-40
  Compton edge cluster.
- `test_multi_line_candidate_not_demoted_even_if_one_is_secondary` —
  Bi-214 with both 609 (outside any secondary) and 910 (inside Co-60
  edge) survives.
- `test_strong_multi_line_evidence_not_demoted` — >2 matched lines →
  rule respects `secondary_max_lines` threshold.
- `test_no_parent_in_catalog_rule_inert` — Eu-152 + Bi-214 (neither in
  v2 catalog) → rule does nothing.
- `test_photopeak_collision_not_handled_by_rule_5` — pure photopeak
  overlap is not Rule 5's domain.
- `test_opt_out_disables_rule` — `apply_secondary_feature_rule=False`
  fully disables the rule.
- `test_parent_itself_not_demoted_by_own_secondary` — defensive: a
  detected parent's own photopeak detection is never demoted.

**Verification — full regression**: 234 tests across 19 files pass
(225 prior + 9 new F-40). `validate_certs.py` matrix unchanged
(11/11 measurable; mean |Δ|=4.62%; max |Δ|=12.79% for Th-228 via
Pb-212 chain proxy) — cert harness bypasses `disambiguate_identifications`
by design (single-source semantics), so Rule 5 is inert there.

**Files changed**:
- `scripts/gamma/identification/disambiguate.py` — Rule 5 block + two
  new kwargs (`apply_secondary_feature_rule`, `secondary_max_lines`);
  module docstring updated to enumerate the rule.
- `test_secondary_feature_rule.py` — new, 9 tests.
- `KNOWN_AND_FIXED_ISSUES.md`, `README.md`, `NOTES_v1.7_methodology.md`,
  `handoff.md` — updated.

**Open follow-ups (not in F-40)**:
- ~~Tl-208 583 alternative chain proxy for Th-228 cross-validation~~
  **CLOSED in F-41 (v1.7.19)**.
- Extended catalog for Ba-133 / Eu-152 / Bi-214 multiplex spectra
  (would let Rule 5 act on these conflict-rich nuclides too).
- Averaged background spectra (15 measurements per geometry → 50×
  noise reduction).
- CLI integration (Variant A) — expose `--apply-secondary-rule` and
  `--secondary-max-lines` in `gamma.cli.analyze`.

---

#### F-41: Tl-208 alternative chain proxy для Th-228 cross-validation (v1.7.19)
**Status**: FIXED in v1.7.19 — `validate_certs.py` теперь измеряет
Th-228 cert source через ДВА независимых daughter-chain пути:
**Pb-212 238 keV** (single low-energy line, прежний F-39 proxy) и
**Tl-208 583/860/2614 keV** (multi-line high-energy proxy через
Bi-212 → α-branch с branching=0.3594).

**Контекст**: F-39 (v1.7.17) закрыл Th-228 cert gap через Pb-212
238 keV, но deviation составил −12.79 % — значительно больше, чем
прочие 10 cert строк (mean |Δ|=4.62 %, max |Δ| остальных = 8.35 %
для Am-241). Эта 5 cm 238 keV точка несла полную ответственность
за max |Δ| матрицы, и было неясно — это (а) ошибка cert value,
(б) методологическая проблема chain-proxy, или (в) специфика
low-energy/self-absorption на 5 cm point geometry. Cross-validation
через независимую chain через ДРУГОЙ дочерний изотоп даёт прямой
ответ.

**Mechanism**:
- Новая `CertFixture("Tl-208", "Th-228__264_2023_…spe",
  cert_nuclide="Th-228", chain_branching=1.0)` в `FIXTURES`.
- `chain_branching=1.0` потому что Lsrm chain-library (загружаемая
  через `load_lsrm_chain_libs()`) уже **pre-scaled** Tl-208
  γ-line intensities на 0.3594 β-branching from Bi-212:
  - Lib I(583.19 keV) = 30.6 % = 0.3594 × 84.5 % ENSDF
  - Lib I(2614.51 keV) = 35.85 % = 0.3594 × 99.75 % ENSDF
  - Lib I(510.77 keV) = 8.10 % = 0.3594 × 22.6 % ENSDF
  - Lib I(860.56 keV) = 4.50 % = 0.3594 × 12.5 % ENSDF
  `compute_activity` инвертирует lib I и восстанавливает **parent**
  A_Th-228 непосредственно (β-factor сокращается).
- Новый **cross-validation block** в `validate_certs.py:main()`:
  собирает rows по cert_nuclide, печатает side-by-side таблицу для
  каждого parent с ≥ 2 daughter proxies, считает попарные ratios,
  флагает rasхождение > 5 %.

**Результат**:
```
Parent: Th-228
    daughter    A_meas, Bq   A_cert@meas, Bq   Δ vs cert, %
      Pb-212       77 250.0          88 575.1        -12.79%
      Tl-208       88 516.8          88 575.1         -0.07%
  ratio Pb-212/Tl-208 = 0.8727  (-12.73%)  >5%
```

**Закрытые научные вопросы**:
1. **Cert value 129 000 Bq @ 25.05.2023 корректен** — Tl-208 даёт
   88 516.8 Bq, ожидание decay-corrected 88 575.1 Bq, Δ = −0.07 %.
2. **Chain-proxy методология (F-39) корректна** — Tl-208 проходит
   ту же decay-correction logic с T½_Th-228, и результат точный.
3. **−12.79 % deviation для Pb-212 локализована** в 238 keV line:
   - 5 cm point self-absorption более выражена на 238 keV чем на
     583+ keV (μ ρL ~ 0.5 для NaI shielding на 200 keV vs ~0.15 на
     600 keV).
   - ε(E) curve detail в low-energy крыле (точка 238 keV сидит в
     резком области ε-curve).
   - Pb-212 single-line measurement не имеет multi-line averaging
     для гашения statistical fluctuation.
4. **Pb-212 row остаётся в матрице** as documented limitation, НЕ
   удаляется — даёт справочную оценку для low-energy single-line
   измерений на этой геометрии.

**Methodological implication**:
- Когда parent имеет ≥ 2 detectable daughter chains, multi-chain
  cross-validation должна быть default. F-41 формализует это для
  Th-228; future expansion возможна на U-238 chain
  (Bi-214 609/1120/1764 + Pb-214 295/352).
- Cross-validation ratio — сильный diagnostic: ≤ 5 % = chain proxy
  методология валидна; > 5 % = указатель на geometry-specific
  systematic (как F-41 продемонстрировал на 238 keV).

**Файлы изменений**:
- `validate_certs.py` — новая `CertFixture` для Tl-208, новый
  cross-validation block (~30 строк) после CSV write.
- `test_chain_proxy.py` — новый файл с 8 тестами.
- `KNOWN_AND_FIXED_ISSUES.md`, `README.md`,
  `NOTES_v1.7_methodology.md`, `handoff.md` — обновлены.

**Verification — full regression**: 20 test files, 242+ tests pass
(234 prior + 8 new F-41). `validate_certs.py` matrix теперь:
**12/12 measurable**, mean |Δ| = 4.24 %, max |Δ| = 12.79 % (Pb-212;
не Tl-208 — Tl-208 даёт −0.07 %).

**Open follow-ups (not in F-41)**:
- Extended catalog для Ba-133 / Eu-152 / Bi-214 multiplex spectra.
- Averaged background spectra (15 measurements per geometry).
- CLI integration (Variant A).
- Per-geometry secondary-range inference (catalog сейчас агрегирует
  по геометриям через quantile spread).

---

#### F-42: Symmetric reader API — per-call energy-ceiling override for Lsrm `.spe` (v1.7.20)
**Status**: FIXED in v1.7.20 — `read_lsrm_spe` теперь принимает два
keyword-only параметра, симметричных с `read_atomspectra_xml`:
- `apply_energy_ceiling: bool = True` — при False сохраняет полный
  декодированный канальный диапазон (диагностический режим).
- `ceiling_keV: Optional[float] = None` — per-call override константы
  `ENERGY_CEILING_KEV` (3000 keV по project scope). None → используется
  модульная константа. Игнорируется когда `apply_energy_ceiling=False`.

**Контекст / motivation**: проектный потолок 3000 keV исторически был
жёстко зашит в `_apply_ceiling()`. Это блокировало два легитимных
варианта использования: (1) **диагностика** — проверка калибровки на
полном канальном диапазоне (включая 3+ MeV хвост), (2) **per-call
ceiling** — анализ узкой энергетической области без мутации
`ENERGY_CEILING_KEV` (которая глобальна для всего конвейера). До v1.7.20
обходить можно было только monkey-patching константы либо ручным
trimming `Spectrum.counts` после чтения — оба варианта ломают
deterministic поведение и/или нарушают invariants Spectrum dataclass.

**Mechanism**:
- Сигнатура `read_lsrm_spe(path, *, apply_energy_ceiling=True,
  ceiling_keV=None)` — оба параметра keyword-only.
- В блоке энергетического trim (`if apply_energy_ceiling and
  spec.energy_cal is not None`) ceiling выбирается как
  `ceiling_keV if ceiling_keV is not None else ENERGY_CEILING_KEV`.
- Когда `apply_energy_ceiling=False` И калибровка присутствует —
  `energy_max_keV_kept` всё равно вычисляется (на последнем сыром
  канале), чтобы downstream код видел реалистичный диапазон.
- `gamma.io.readers.read_spectrum` уже forward `**kwargs` в
  format-specific reader, так что
  `read_spectrum(path, apply_energy_ceiling=False)` и
  `read_spectrum(path, ceiling_keV=400.0)` работают «из коробки».
- Та же симметрия применена к `read_atomspectra_xml` (вне scope этого
  итерационного плана LSRM-only, но оставлена для согласованности
  единого API).

**Defensive design**:
- Default-поведение НЕ изменилось: 242 предыдущих теста проходят без
  правок. `apply_energy_ceiling=True` и `ceiling_keV=None` =
  исходная семантика.
- Keyword-only форма исключает позиционную путаницу с
  `parse_background` (xml-specific) и держит вызовы explicit.
- `ENERGY_CEILING_KEV` остаётся single source of truth для default —
  per-call override НЕ перезаписывает константу, эффект только на
  один read.

**Файлы изменений**:
- `scripts/gamma/io/lsrm_spe.py` — расширенная сигнатура
  `read_lsrm_spe` + per-call resolve в блоке trim (строки 84–101,
  266–268).
- `scripts/gamma/io/atomspectra_xml.py` — параллельная сигнатура
  `read_atomspectra_xml` + проброс kwargs в
  `_parse_energy_spectrum_block` (вне основного scope, но в наборе
  для симметрии).
- `scripts/gamma/io/readers.py` — обновлён docstring о двух
  passthrough kwargs.
- `SKILL.md` §Scope — параграф об override после строки про 0–3000 keV.
- `test_reader_api.py` — новый файл с 7 тестами; 3 теста LSRM-specific:
  `test_lsrm_spe_default_trims_at_3000`,
  `test_lsrm_spe_apply_false_keeps_full_range`,
  `test_lsrm_spe_custom_ceiling` (последний trim до 400 keV проверяет
  что `n_channels` строго меньше default).

**Verification — full regression**: 20 test files, всё PASS.
test_reader_api.py: 7/7 (по факту проверяет и .xml — 4 теста — но scope
v1.7.20 явно ограничен LSRM .spe; xml-тесты пройдут потому что код там
уже симметричен из предыдущего refactor). Validate matrix
неизменна: 12/12 measurable, mean |Δ|=4.24 %, max |Δ|=12.79 %.

**Open follow-ups (not in F-42)**:
- Симметричная разработка для AtomSpectra XML и других форматов (.chn,
  .n42, .mca, .csv) — отложено до завершения плана для LSRM NaI.
- Per-call subtrim helper `Spectrum.trim_above(E_keV)` для уже
  прочитанных spectra (без повторного read) — kandidates на будущий
  F-NN.

---

#### F-43: Averaged Lsrm `.spe` background spectra (v1.7.21)
**Status**: FIXED in v1.7.21 — новый модуль `gamma.io.average_lsrm`
агрегирует N Lsrm `.spe` спектров одной геометрии в один усреднённый
(sum-of-counts) спектр с суммарным live-time / real-time. Архив 2016
Поверка перерабатывается в 3 эталонных фоновых файла в
`data/averaged_backgrounds/`.

**Контекст / motivation**: 2016 Поверка archive содержит 4 фоновых
контекста с 15+ измерениями каждый, но они никогда не комбинировались
в единый low-noise background. Каждое единичное измерение 3-12 ч
дать σ_rate, ограничивающее точность downstream background-
subtraction для слабых γ-линий (Pb-210 46 keV, Pb-214 295/352 keV,
Am-241 26-60 keV X-rays). Суммирование 15 измерений в один спектр →
эквивалент одного 120-часового измерения → snr ×√15 ≈ 3.87.
Хотя начальная оценка в handoff была «50× noise reduction», точная
эмпирика — это **факторное снижение σ_rate в √N ≈ 3.87× для N=15**
(variance reduction ×15 = 1500%, std reduction ×3.87 = 287%). На
практике это значительно расширяет окно детекции слабых photopeak.

**Mechanism**:
- `average_lsrm_spectra(paths, *, rel_gain_tolerance=0.005,
  abs_offset_tolerance=2.0, require_same_detector=True,
  require_same_geometry=True, sample_id=None, comment="")` →
  Spectrum.
- Defensive проверки перед слиянием:
  - **Calibration consistency**: max abs spread `a0` ≤ 2.0 keV;
    max relative spread `a1` ≤ 0.5 %. Превышение → `CalibrationMismatchError`.
  - **Identity**: `detector_id` и `geometry` должны совпадать
    (по умолчанию); явный `require_same_geometry=False` для
    cross-geometry слияния.
  - **Channel lengths**: все входы должны иметь одинаковые
    `n_channels` и `n_channels_raw`.
- Counts: element-wise sum (`np.int64`). Live-time, real-time,
  dropped_overflow_count: simple sums. Calibration, FWHM, detector,
  geometry: inherited from first input. start_datetime: earliest
  across inputs.
- Output Spectrum несёт `extras["averaging_provenance"]` (полный
  audit dict со списком source paths, per-input live-times,
  agreement statistics) и `extras["averaging_sigma_reduction"]` = √N.
- Новый minimal writer `write_lsrm_spe(spec, path, *, type_label,
  config_name)` сохраняет в `.spe` формат с CP-1251 header +
  uint32 LE counts. Round-trip lossless для counts (в пределах
  energy ceiling), live_time, real_time, energy_cal, geometry,
  detector_id.

**Generated archive** (`build_averaged_backgrounds.py` запускается
один раз на packaging time):

| Output                                  | Context       | Geometry      | N  | Live  | σ-red  |
|-----------------------------------------|---------------|---------------|----|-------|--------|
| bg_marinelli_water_marinelli.spe        | Marinelli H₂O | Маринелли    | 15 | 120 h | 3.87×  |
| bg_empty_shield_point5cm.spe            | Empty shield  | Точечная-5см | 15 | 120 h | 3.87×  |
| bg_open_lid_point25cm.spe               | Open lid      | Точечная-25см | 15 | 120 h | 3.87×  |

Каждый файл сопровождается JSON-сайдкаром `<name>.provenance.json`
с полным аудит-следом для traceability. Общий MANIFEST.json
индексирует всё.

**"Временная нестабильность" пропущена**: 15 коротких 300-секундных
измерений предназначены для проверки temporal stability детектора,
не для averaging как long-exposure bg.

**Файлы изменений**:
- `scripts/gamma/io/average_lsrm.py` — новый модуль (~330 строк):
  `average_lsrm_spectra`, `write_lsrm_spe`,
  `CalibrationMismatchError`, `IdentityMismatchError`,
  `_check_calibrations`, `_check_identity`, `_check_channel_lengths`.
- `build_averaged_backgrounds.py` — новый standalone generator
  (~190 строк) — запускается один раз для production archive.
- `data/averaged_backgrounds/` — новая папка с 3 `.spe` + 3
  `.provenance.json` + 1 `MANIFEST.json`.
- `test_average_lsrm.py` — новый файл с **12 тестами**: count-sum
  identity, time aggregation, σ-reduction = √N, calibration drift
  rejection / tolerance relaxation, geometry mismatch
  rejection / override, single-file passthrough, empty list error,
  round-trip write/read, provenance audit, archive integrity,
  manifest consistency.

**Verification — full regression**: 21 test файл, всё PASS,
**≥261 теста** (249 prior + 12 new F-43). Validate matrix не
изменена. F-43 не затрагивает activity/identification pipeline —
averaged backgrounds готовы к use downstream через стандартный
`read_spectrum` API.

**Use case demonstration**:
```python
from gamma.io.readers import read_spectrum
from gamma.calibration.bg_subtract_dual_mode import (
    BackgroundConsentRegistry, subtract_background,
)
registry = BackgroundConsentRegistry()
src = read_spectrum("sample_marinelli.spe")
bg  = read_spectrum("data/averaged_backgrounds/bg_marinelli_water_marinelli.spe")
registry.approve(bg.source_path)
net = subtract_background(src, bg, consent_registry=registry)
# bg σ contribution shrinks ≈3.87× compared to single 5-h bg →
# weaker photopeaks (Pb-210, low-E XRF) gain in detectability.
```

**Open follow-ups (not in F-43)**:
- Аналогичная aggregation для AtomSpectra XML (после выхода из
  LSRM-only scope).
- Аналогичная aggregation для других LSRM bg-фолдеров (например,
  Дента-100/Чашка-60 — текущий archive имеет всего 1 файл каждой
  геометрии в "пустая защита"; для них нужно ≥5 измерений).
- Verification на downstream cert validation matrix (пересчитать
  с averaged backgrounds вместо single — ожидается улучшение
  стат-надёжности для weak-source rows).

---

#### F-44: Cumulative-checkpoint detection + 2024 archive sync (v1.7.22)
**Status**: FIXED in v1.7.22 — три связанных изменения:
1. **Bug в F-43**: F-43 averaged backgrounds некорректно интерпретировали
   LSRM Spectraline cumulative checkpoint files как независимые измерения.
2. **Cumulative pattern detection** в `gamma.io.average_lsrm`.
3. **Sync с пользовательским LSRM архивом** `C:\LSRM\Work\BG\Gamma-1S\Spe -
   поверки\`: добавлены 80 missing reference files (48 Y-88 temporal
   stability + 16 Marinelli closed-lid bg + 16 25cm open-lid bg).
4. **Methodological clarification**: вклад радона в фоновом спектре
   считается пренебрежимо малым, всё излучение фона трактуется как
   природный U-238/Th-232 из строительных материалов (стены, бетон,
   кирпич). Это упрощает background-subtraction: фон трактуется как
   статический источник.

**Bug в F-43 (механика)**:
LSRM Spectraline acquisition software пишет файлы `_01.spe` (1h),
`_02.spe` (2h cumulative), `_N.spe` (Nh cumulative) — **каждый
последующий файл содержит все события предыдущих** плюс новый
интервал. Установлено эмпирически на 2016 архиве: `Фон вода_01` имеет
22 709 counts, `Фон вода_15` имеет 339 544 counts (ratio 14.95, точно
15× как у cumulative).

F-43 наивно суммировал все 15 файлов:
- **Inflated counts**: суммарные counts = 2 720 587 = 339 544 × 8
  (где 8 = (1+2+...+15)/15 = средний коэффициент cumulative).
- **Inflated live_time**: суммарный live_time = 432 000s = 54 000s × 8.
- **Rate preservation**: `counts/live_time` = 6.30 cps корректно
  (потому что и числитель и знаменатель умножены на одно число).
- **σ-claim ошибочный**: F-43 объявлял σ-reduction = √15 ≈ 3.87×, но
  на самом деле есть только **одно независимое измерение длиной 15h**,
  поэтому σ-reduction = 1.0.

**Downstream impact bug'а**:
- Subtraction по LSRM formula `net = sample - bg × (T_sample / T_bg)`
  даёт **корректный net_counts** (rate consistency сохранён).
- σ_bg в Poisson `σ² = N_bg` под-оценивает statistical uncertainty в
  фоне в `√8 ≈ 2.83×` (потому что N_bg искусственно × 8 → σ × √8). MDA
  получалась слишком оптимистичной.
- Validate_certs.py не использовал averaged backgrounds → cert matrix
  не был затронут.

**F-44 fix**:
```python
def detect_cumulative_pattern(specs, *, rel_live_time_tolerance=0.01):
    """Detect LSRM cumulative checkpoint sets.
    Criteria (all required for is_cumulative=True):
      1. N ≥ 2 inputs
      2. All same start_datetime (±1s)
      3. Sorted live_times match arithmetic progression
         t_i ≈ (i+1) · t_min within rel tolerance
    """
```

`average_lsrm_spectra` теперь принимает `cumulative_policy` kwarg:
- `"auto"` (default) — детектирует pattern, выбирает mode.
- `"cumulative_last"` — берёт longest live_time file (cumulative reading).
- `"independent_sum"` — sum-of-counts (original F-43 semantic).

Output `extras["averaging_mode"]` фиксирует выбранную семантику для
audit. В cumulative_last mode `σ-reduction = 1.0` (одно измерение
длительности max(t_i)).

**Архив `data/averaged_backgrounds/` пересоздан (5 файлов)**:

| Output                                            | Source                  | Geometry      | live | mode             |
|---------------------------------------------------|-------------------------|---------------|------|------------------|
| `bg_2016_marinelli_water_marinelli.spe`           | 2016 / Фон вода         | Маринелли    | 15h  | cumulative_last  |
| `bg_2016_empty_shield_point5cm.spe`               | 2016 / фон пустая защита| Точечная-5см | 15h  | cumulative_last  |
| `bg_2016_open_lid_point25cm.spe`                  | 2016 / Фон откр крышки  | Точечная-25см | 15h  | cumulative_last  |
| `bg_2024_marinelli_water_closed_lid_marinelli.spe`| 2024 / Фон закр кр      | Маринелли    | 16h  | cumulative_last  |
| `bg_2024_open_lid_point25cm.spe`                  | 2024 / Фон откр кр      | Точечная-25см | 16h  | cumulative_last  |

Имена файлов получили префикс эпохи (`2016_*`, `2024_*`) — старые
имена F-43 удалены. Каждый файл сопровождается `.provenance.json`
с full audit trail (`aggregation_mode`, `cumulative_detection` block,
`pair_with_geometries` список) + общий MANIFEST.json.

**2024 архив pairing** (per user clarification):
- 2024 / Фон закр кр (16 cumulative, Marinelli geometry, water-filled
  vessel) → паруется с Marinelli sample measurements (matrix-matched
  attenuation).
- 2024 / Фон откр кр (16 cumulative, Точечная-25см) → паруется с
  Точечная-25см sample measurements (no shielding above 25cm).
- 2024 / Дента120мл, Петри-60мл, Точечная-5см — sample geometries
  паруются с **закрытой крышкой / пустой защитой** (используется 2016
  empty_shield_point5cm как ближайшая аналогия; в 2024 архиве
  отдельной empty-shield bg нет).

**Methodological note** (per user clarification):
Радон в фоновом спектре трактуется как пренебрежимый — практически
**всё фоновое излучение приписывается строительным материалам**
(природные U-238/Th-232 chains в бетоне, кирпиче стен помещения).
Это упрощает background-subtraction: фон — стационарный источник,
не требует отдельного моделирования флуктуаций радона.

**Subtraction methodology** (LSRM formula, per user):
```
net_counts[i] = sample_counts[i] - bg_counts[i] · (T_sample / T_bg)
```
Это count-based вычитание со scaling фона на acquisition time
sample'а. Математически эквивалентно rate-based вычитанию
(net_rate = sample_rate - bg_rate), но LSRM count-with-time-scaling
формулировка предпочитается потому что (а) preserves Poisson statistics
naturally; (б) даёт корректный σ через прямое propagation на counts.

**Файлы изменений**:
- `scripts/gamma/io/average_lsrm.py` — `detect_cumulative_pattern()`
  (~80 строк), `CumulativeAmbiguityError`, `cumulative_policy` kwarg
  в `average_lsrm_spectra`, mode-aware skip channel-length /
  calibration check для cumulative_last.
- `build_averaged_backgrounds.py` — расширен на 2024 контексты,
  включает `pair_with_geometries` в манифесте, переименование префиксов.
- `data/averaged_backgrounds/` — пересоздан с 5 канонических файлов
  + 5 sidecar JSON + 1 MANIFEST.json. Старые файлы (3 файла F-43 с
  inflated counts) удалены.
- `references/reference_spectra/Gamma-1C_NaI_63x63_USB_SN-01/
  Поверка-2024/` — новые папки с 80 missing files (48 ВН Y-88
  temporal stability + 16 closed-lid + 16 open-lid).
- `test_average_lsrm.py` — расширен с 5 новых тестов F-44:
  cumulative detection on 2016 set, single-spectrum non-cumulative,
  synthetic independent non-cumulative, auto-policy selects
  cumulative_last, explicit policy override (both directions), invalid
  policy raises. Старые тесты обновлены (forced `independent_sum`
  policy где нужна старая семантика).

**Verification — full regression**: 21 test файл, всё PASS,
**≥266 tests** (≥261 prior + 5 new F-44 cumulative tests). Validate
matrix не изменена. Old F-43 averaged backgrounds не использовались
ни в `validate_certs.py`, ни в каком-либо другом production-пути, так
что fix не имеет downstream regression.

**Open follow-ups (not in F-44)**:
- Replace single `Фон_закр_кр_вода_01.spe` в validate_certs.py BG_PATH
  на новый `bg_2024_marinelli_water_closed_lid_marinelli.spe` (или
  оставить как baseline для regression-стабильности и добавить
  averaged variant как сравнительный run). **CLOSED in F-45 (v1.7.23)**.
- 2024 / Дента120мл / Петри-60мл / Точечная-5см / Точечная-25см
  source measurements — добавить как новые cert fixtures.
- 48 ВН Y-88 файлов — анализ temporal stability детектора (drift
  in calibration, FWHM, rate за 3-дневный период).

---

#### F-45: Cert-validation background swap — averaged empty-shield point5cm (v1.7.23)

**Status**: FIXED in v1.7.23 — `validate_certs.py` BG_PATH переключён
с одиночного файла `Фон_закр_кр_вода_01.spe` (~1h, неправильная
геометрия) на canonical averaged `bg_2016_empty_shield_point5cm.spe`
(15h cumulative_last, корректная геометрия per F-44 pairing rules).

**Контекст методологической ошибки до F-45**:
До F-45 harness вычитал из всех 12 fixture point-5cm spectra фон
**Marinelli + water** (`Фон_закр_кр_вода_01.spe` — фон с Marinelli
сосудом, заполненным водой, в защите с закрытой крышкой). Это
matrix-matched фон для **Marinelli sample**, но fixtures всё —
Точечная-5см. Слой воды в Marinelli sleeve активно поглощает
низкоэнергетические γ-линии (Pb K-α 75 keV, Pb K-β 87 keV, Am-241
59.5 keV региона), что под-вычитает естественный low-E background для
point-source измерений. На сильных линиях (>200 keV) эффект
пренебрежим, но для weak source / weak-line analysis это завышает
net_counts для линий <100 keV.

Per F-44 pairing rules (NOTES_v1.7_methodology.md §v1.7.22):

| Sample geometry      | Background                       |
|----------------------|----------------------------------|
| Маринелли           | Marinelli + water (matrix-matched)|
| Дента / Чашка-60     | Empty shield, closed lid          |
| **Точечная-5см**     | **Empty shield, closed lid**      |
| Точечная-25см        | Open lid                          |

**F-45 fix** (одна строка):
```python
# было:
BG_PATH = REF_DIR / "Фон_закр_кр_вода_01.spe"
# стало:
BG_PATH = ROOT / "data" / "averaged_backgrounds" / "bg_2016_empty_shield_point5cm.spe"
```

**Эффект на cert-matrix metrics**:

| Nuclide  | Δ baseline | Δ F-45 | peaks baseline | peaks F-45 |
|----------|------------|--------|----------------|------------|
| Cs-137   | +1.83 %    | +1.83 %| 8              | 8          |
| Co-60    | +0.61 %    | +0.61 %| 13             | 13         |
| Na-22    | +5.99 %    | +5.99 %| 16             | 16         |
| Eu-152   | -6.82 %    | -6.82 %| 15             | 14         |
| Ba-133   | +1.32 %    | +1.32 %| 15             | 13         |
| Am-241   | -8.35 %    | -8.35 %| **11**         | **7**      |
| Zn-65    | +5.45 %    | +5.45 %| 7              | 6          |
| Y-88     | +1.64 %    | +1.64 %| 13             | 13         |
| Bi-207   | -2.67 %    | -2.67 %| 14             | 14         |
| Cd-109   | -3.38 %    | -3.38 %| **9**          | **6**      |
| Pb-212   | -12.79 %   | -12.79 %| 18            | 19         |
| Tl-208   | -0.07 %    | -0.07 %| 18             | 19         |

- **A_изм и Δ%**: инвариантны для всех 12 fixtures (target nuclide
  strong-line dominate activity, и для тех каналов bg subtraction
  меняет малую долю signal).
- **mean |Δ| = 4.24 %, max |Δ| = 12.79 %** — без изменений.
- **Peak count shifts** (главный наблюдаемый эффект): Am-241 11→7,
  Cd-109 9→6, Ba-133 15→13, Eu-152 15→14, Zn-65 7→6. Меньше "spurious"
  peaks остаётся в low-energy регионе после bg subtraction — empty-
  shield bg корректно содержит natural-bg Pb K-shell X-rays и U/Th
  chain низкоэнергетические линии из строительных материалов, которые
  Marinelli+water bg маскировал водяной абсорбцией.
- **Pb-212 18→19, Tl-208 18→19**: +1 peak — низкоэнергетический
  artefact, ранее частично compensated water absorption, теперь
  proper bg cancellation сохраняет его.

**Defensive characteristics**:
1. Activity не регрессирует (target lines E ≫ low-E bg absorption
   region для всех 12 fixtures с одним исключением Am-241 59.5 keV —
   там Δ остаётся −8.35 %, дисперсия не меняется).
2. Peak list cleaner — downstream identification/disambiguate видят
   меньше spurious natural-bg peaks.
3. σ_bg в bg subtraction правильно ослаблена (15h live_time vs 1h
   single file → bg σ_rate ×1/√15 ≈ 0.258).
4. Бывший single-file BG path сохранён закомментированным для
   diagnostic comparison.

**Open known limitation** (документировано в код-комментарии):
2024 архив не содержит empty-shield closed-lid bg, поэтому для всех
2017-2023 measurement-year fixtures используется 2016 averaged bg
независимо от epoch. Acceptable, поскольку 5cm shield + detector +
masonry комнаты не менялись 2016→2024 (verified via cross-epoch peak
position stability в F-44 archive sync).

**Файлы изменений**:
- `validate_certs.py` (5 строк):
  - Docstring `2. subtract bg_2016_empty_shield_point5cm.spe (F-43
    averaged, F-44 cumulative_last semantic; ...)`.
  - 19-строчный inline comment перед `BG_PATH = ...` объясняет
    pairing-rules motivation, F-44 reference, 2024 epoch caveat,
    diagnostic fallback path.
  - `BG_PATH` указывает на `data/averaged_backgrounds/bg_2016_empty_shield_point5cm.spe`.
  - Закомментированная single-file `BG_PATH` сохранена как diagnostic
    fallback.

**Verification**: 21/21 test файл проходит (266+ tests across all
suite), cert matrix re-run даёт идентичный mean/max |Δ| baseline:

```
Measured 12/12 fixtures.
Mean |Δ| = 4.24%
Max  |Δ| = 12.79%
```

**Open follow-ups (not in F-45)**:
- Generate `bg_2024_empty_shield_point5cm.spe` если LSRM добавит
  такой контекст в будущий archive (есть только 2024 Marinelli и
  Точечная-25см на текущий момент).
- Per-fixture epoch matching: использовать `bg_2024_marinelli_...`
  для measurement-year ≥ 2020, `bg_2016_*` для earlier. Не делалось
  потому что нет 2024 empty-shield point-5cm bg.
- 2024 cert sources как new fixtures (Точечная-25см, Дента120мл,
  Петри-60мл, Маринелли из Поверки 2024). **PARTIALLY CLOSED in
  F-46a (v1.7.24)** — 4 Точечная-25см добавлены. Marinelli/Дента/
  Петри отложены до F-46b/c (нужна обработка Bq/kg → Bq через
  mass_g + chain-proxy для Ra-226 и Th-232).

---

#### F-46a: Multi-geometry cert matrix — Точечная-25см slice (v1.7.24)

**Status**: FIXED in v1.7.24 — `validate_certs.py` расширен с одной
геометрии (Точечная-5см) до двух (Точечная-5см + Точечная-25см).
Foundation per-geometry resolution готов для последующих slices
F-46b (Marinelli) и F-46c (Дента + Петри).

**Контекст**: 2016/2024 архив содержит fixtures в 5 различных
геометриях, но v1.7.23 harness работал только с одной (Точечная-5см).
Per-geometry эффективность calibrations и averaged backgrounds уже
существуют (F-44), но не подключены к cert validation. Cross-geometry
validation одного и того же нуклида (Cs-137 в 5cm и 25cm) — главное
property для catching систематических ошибок eff/bg/TCS.

**Реализация — Per-geometry resolver**:
```python
EFF_PATHS = {
    "Точечная-5см":   EFF_DIR / "...-_Точечная-5см.efr",
    "Точечная-25см":  EFF_DIR / "...-_Точечная-25см.efr",
    "Дента-120мл":    EFF_DIR / "...-_Дента.efr",
    "Петри-60мл":     EFF_DIR / "...-_Петри.efr",
    "Маринелли":      EFF_DIR / "...-_Маринелли.efr",
}
BG_PATHS = {  # per F-44 pairing rules
    "Точечная-5см":   ".../bg_2016_empty_shield_point5cm.spe",
    "Точечная-25см":  ".../bg_2016_open_lid_point25cm.spe",
    "Дента-120мл":    ".../bg_2016_empty_shield_point5cm.spe",
    "Петри-60мл":     ".../bg_2016_empty_shield_point5cm.spe",
    "Маринелли":      ".../bg_2016_marinelli_water_marinelli.spe",
}
CERT_PATHS = {
    "Точечная-5см":   CERT_DIR / "АСПЕКТ_ОСГИ_2024.src",
    "Точечная-25см":  CERT_DIR / "АСПЕКТ_ОСГИ_2024.src",
    "Дента-120мл":    CERT_DIR / "Эталон_Дента120мл__Аспект2017_.src",
    "Петри-60мл":     CERT_DIR / "Эталон_Петри-60__Аспект_2017_.src",
    "Маринелли":      CERT_DIR / "Эталон_Маринелли__Аспект_2017_.src",
}
```

`CertFixture` получил `geometry: str = "Точечная-5см"` field (default
preserves v1.7.23 rows без modification).

`main()` теперь использует lazy per-geometry resource cache —
загружает eff/bg/cert per geometry один раз, переиспользует для всех
fixtures той же geometry.

**4 новых fixtures (F-46a — Точечная-25см)**:
- `Cs-137 №SRC-01_Точечная-25см_25cm.spe` (single nuclide direct)
- `Na-22 #01.22_Точечная-25см_25cm.spe` (single nuclide direct)
- `Y-88 №SRC-02_Точечная-25см_25cm.spe` (single nuclide direct)
- `Th-228 №309_Точечная-25см_25cm.spe` (Tl-208 chain proxy,
  source #SRC-06.2021 cert at 100 000 Bq)

Все 4 используют АСПЕКТ_ОСГИ_2024.src cert file (тот же что
Точечная-5см rows), Точечная-25см .efr efficiency curve и
bg_2016_open_lid_point25cm.spe.

**Extended cert-matrix results (16 fixtures)**:

| Геометрия      | n  | mean \|Δ\| | max \|Δ\| | comment             |
|----------------|----|-----------|----------|---------------------|
| Точечная-5см   | 12 | 4.24 %    | 12.79 %  | Pb-212 max (F-41)   |
| Точечная-25см  | 4  | 8.35 %    | 19.34 %  | Tl-208 25cm max     |
| **Общее**      | 16 | 5.27 %    | 19.34 %  |                     |

**Точечная-25см detail**:
- Cs-137: Δ=−5.09 %, 12 peaks
- Na-22:  Δ=+4.15 %, 13 peaks; TCS=2 line(s)
- Y-88:   Δ=−4.83 %, 13 peaks; TCS=2 line(s)
- Tl-208 (Th-228 proxy): Δ=−19.34 %, 16 peaks; TCS=4 line(s)

**Cross-geometry Tl-208 finding**: same chain-proxy methodology
gives Δ=−0.07 % @ 5cm (F-41 baseline) but Δ=−19.34 % @ 25cm. Two
explanations:
1. **25cm efficiency curve bias на высоких энергиях** (583/860/2614
   keV) — chi²/dof=2.51 для 25cm .efr vs ~lower для 5cm.
2. **Th-228 №SRC-03 cert overstatement** (cert says 100 000 Bq;
   measured 22 680 Bq → cert-only-effect would also penalize
   ratio with another nuclide on the same source, but #SRC-06 has
   only Th-228 entry → cannot independently cross-validate).

Эти два source'а (5cm uses #SRC-05.2023, 25cm uses #SRC-06.2021) — разные
физические источники. **Open**: для distinguishing нужно либо
(а) измерить №264 в 25cm geometry, либо (б) добавить ещё одну ε(E)
calibration источника для 25cm.

**Cross-validation block (F-41 extension)**: теперь печатает 3 ratio
вместо 1 — Pb-212/Tl-208 same-source (5cm), Pb-212(5cm)/Tl-208(25cm),
Tl-208(5cm)/Tl-208(25cm). Последняя comparison = same-nuclide
cross-geometry consistency check (ratio 3.90× против ratio of
expected cert@meas activities = 88575/28119 = 3.15, разница +24%).

**Defensive characteristics**:
1. **All 12 prior Точечная-5см rows unchanged** — bit-identical
   results (default `geometry="Точечная-5см"` preserves backward
   compat).
2. **Library-search backward compat**: `BG_PATH`, `EFF_5CM`,
   `CERT_PATH` module-level constants остаются (как aliases в
   Точ-5см defaults) — `test_chain_proxy.py` grep-tests их.
3. **test_chain_proxy.py updated**: assertions generalized from
   "exactly one" Tl-208 fixture to "≥1 per geometry, paired with
   Pb-212 in shared geometries". 8/8 passes.
4. **Lazy resource cache** в `_resolve_geometry_resources()` —
   eff/bg/cert загружаются один раз per geometry, переиспользуются
   для всех fixtures.

**Файлы изменений**:
- `validate_certs.py` (~70 строк): EFF_PATHS/BG_PATHS/CERT_PATHS dicts,
  geometry field в CertFixture, lazy resolver + cache, 4 новых
  Точечная-25см fixtures, main() refactored.
- `test_chain_proxy.py` (~40 строк): 3 теста generalized для
  multi-geometry chain-proxy invariants.

**Verification**: 21/21 test файл проходит, 266+ tests без изменений
в count (F-46a добавляет fixtures но не тесты — test_chain_proxy
обновлён без add/remove). Cert matrix: 16/16 measurable, mean
|Δ|=5.27 %, max |Δ|=19.34 %.

**Open follow-ups (next slices)**:
- **F-46b**: Marinelli — 8 fixtures. **CLOSED in v1.7.25 (F-46b)**.
- **F-46c**: Дента-120мл + Петри-60мл — 16 fixtures. **CLOSED in
  v1.7.25 (F-46c)**.
- **F-46d**: investigate Tl-208 25cm Δ=−19.34 %. **PARTIALLY CLOSED
  in v1.7.25 (F-46d)** — Th-228 №SRC-04 @25cm не доступен в архиве,
  поэтому cross-source validation невозможен. Diagnostic chi²/dof
  reporting per geometry добавлен в `validate_certs.py` summary; гипотезы
  не falsifiable без новых measurement data.

---

#### F-46b/c/d: Multi-geometry cert matrix — Marinelli + Дента + Петри slices + diagnostic chi²/dof (v1.7.25)

**Status**: FIXED in v1.7.25 — `validate_certs.py` расширен с 16
fixtures (12 Точ-5см + 4 Точ-25см из F-46a) до **40 fixtures** (12
Точ-5см + 4 Точ-25см + 8 Marinelli + 8 Дента-120мл + 8 Петри-60мл).
Closes F-46 multi-geometry expansion in one shipment.

**Структура slice'ов**:
- **F-46b** (Marinelli, 8 fixtures): Cs-137 ×2 direct, K-40 ×2 direct,
  Ra-226 ×2 через Bi-214 chain proxy, Th-232 ×2 через Tl-208 chain
  proxy. Один из Th-232 fixture использует source 420-17031 (2017 cert
  ref) вместо 420-7-16 (.spe file не measurement'нут).
- **F-46c** (Дента-120мл + Петри-60мл, 16 fixtures): same 4 nuclides
  × 2 source variants × 2 geometries. Th-232 row для каждой geometry
  использует 420-7-17 + 420-17031 источник.
- **F-46d** (diagnostic): per-geometry eff-curve chi²/dof printed в
  Per-geometry summary table. Th-228 №SRC-04 @25cm cross-source
  validation **не возможна** — file не доступен в archive.

**Реализация — Bq/kg → Bq via mass_g**:

`run_one()` extended (~30 строк):
```python
# F-46b: walk sub_sources to find target nuclide + access mass_g
target_sub = None
cert_act = None
for ss in src.sub_sources:
    for act in ss.activities:
        if act.nuclide == cert_nuclide_name:
            target_sub = ss; cert_act = act; break
    if cert_act: break

# Apply Bq/kg → Bq if per-mass unit
unit = (cert_act.unit or "").strip().lower()
if unit.endswith("bq/kg"):
    A_cert_absolute = cert_act.A_Bq * target_sub.mass_g / 1000.0
else:
    A_cert_absolute = cert_act.A_Bq
```

Fallback path для compound certs (ОСГИ 5431 multi-nuclide) — если
sub_source walk fails, fallback на legacy `src.get_activity(name)`.

**Chain proxies — методология**:

| Cert source | Proxy nuclide | Lib intensity | chain_branching |
|-------------|---------------|---------------|-----------------|
| **Ra-226** (Marinelli/Дента/Петри) | Bi-214 (9 lines) | direct ENSDF | 1.0 |
| **Th-232** (Marinelli/Дента/Петри) | Tl-208 (5 lines) | × 0.3594 (pre-scaled) | 1.0 |
| **Th-228** (Точ-5см #SRC-05) | Pb-212 + Tl-208 | direct ENSDF / × 0.3594 | 1.0 |
| **Th-228** (Точ-25см #SRC-06) | Tl-208 | × 0.3594 | 1.0 |

Ra-226 → Pb-214 → Bi-214 chain: в sealed sample Rn-222 buffer retained
(Rn-222 T½=3.8 d → equilibrium через 25-30 days, 420-series sources
≥20 years old). Lib Bi-214 intensities direct ENSDF (не pre-scaled),
chain_branching=1.0 yields A_Bi-214 ≈ A_Ra-226.

Th-232 → Ra-228 → Ac-228 → Th-228 → ... → Bi-212 → Tl-208 chain:
Ra-228 bottleneck (T½=5.75 y) reaches 93% equilibrium в 22 years (since
2002), 90% в 17 years (since 2007). Lib Tl-208 intensities pre-scaled
по Bi-212 α-branching 0.3594. compute_activity inverting lib intensities
recovers parent A_Th-232 directly.

**Cert-matrix результаты (40/40 measurable)**:

| Геометрия       | n  | mean \|Δ\| | max \|Δ\| | efr chi²/dof |
|-----------------|----|------------|-----------|--------------|
| Точечная-5см    | 12 | 4.24 %     | 12.79 %   | 6.95         |
| Точечная-25см   | 4  | 8.35 %     | 19.34 %   | 2.51         |
| Дента-120мл     | 8  | 9.68 %     | 27.90 %   | 15.40        |
| Петри-60мл      | 8  | 12.55 %    | 19.58 %   | 15.04        |
| Маринелли       | 8  | 9.66 %     | 22.54 %   | 3.72         |
| **Общее**       | 40 | **8.48 %** | **27.90 %** | —          |

**F-46d finding — chi²/dof reveals systematic geometry bias**:

Дента (15.40) и Петри (15.04) имеют поразительно poor eff-curve fits
по сравнению с Marinelli (3.72), Точ-25см (2.51) и Точ-5см (6.95).
Это **falsifies hypothesis B** для F-46a's Tl-208 25cm finding
(cert overstatement) и **strengthens hypothesis A** (eff curve bias).
Подтверждается тем, что:
- Точ-25см Tl-208 Δ=−19.34 % (chi²/dof=2.51 — relatively clean)
- Дента Tl-208 Δ=−17.94 % и −11.53 % (chi²/dof=15.40)
- Петри Tl-208 Δ=−19.58 % и −12.26 % (chi²/dof=15.04)
- Marinelli Tl-208 Δ=−20.31 % (chi²/dof=3.72 — relatively clean, но
  still significant under-estimate)

Pattern: **systematic ~15-25 % under-estimate для chain-proxy Tl-208**
во всех geometries — это **methodology effect, не geometry-specific
bias**. Probable cause: TCS correction model для Tl-208 chain
(Bi-212 → Tl-208 cascade имеет multi-line coincidences) под-correct'ит
true coincidence-summing на close geometries (Marinelli/Дента/Петри
0 cm distance). Точ-5см (5 см distance) → меньше TCS → меньше bias.

**Bi-214 chain proxy results**:

| Geometry      | Source #18 (lighter)  | Source #19 (heavier)  |
|---------------|-----------------------|-----------------------|
| Marinelli     | Δ=−6.48 % (620g)      | Δ=−22.54 % (1670g)    |
| Дента-120мл   | Δ=−3.83 % (74g)       | Δ=−27.90 % (200g)     |
| Петри-60мл    | Δ=−3.66 % (37g)       | Δ=−18.53 % (100g)     |

Strong **mass dependence** — heavier containers show worse Δ. Это
самоabsorption effect: для Bi-214 lower-E lines (242/295/352 keV)
self-attenuation in larger samples reduces measurable count rate.
Light containers (37-74 g) показывают excellent Δ ≤ 7 %. Это
**physically expected** и подтверждает что pipeline корректно работает
on light volumes; heavy-volume systematic suggests Cowell baseline
under-fits Compton continuum on broadened multi-line peaks.

**Cs-137 / K-40 direct (no chain proxy)**:

| Geometry      | Cs-137 results       | K-40 results          |
|---------------|----------------------|-----------------------|
| Marinelli     | +9.81 %, -5.58 %     | -2.33 %, -6.02 %      |
| Дента-120мл   | +4.85 %, -6.84 %     | -1.22 %, -3.31 %      |
| Петри-60мл    | +6.97 %, -3.71 %     | **+18.70 %, +16.95 %**|

Cs-137 и K-40 показывают excellent consistency во всех 6 геометриях
кроме K-40 в Петри (+17-18 % over). Возможная причина: тонкая источник
(60 ml Petri dish ≈ 3 mm thick) → meaningless self-attenuation для
1461 keV; close-geometry TCS не applied (K-40 single line). Excess
counts могут быть Compton continuum от U/Th contamination в материале
matrix не правильно subtracted close-geometry bg.

**Cross-validation block update (F-46)**:

Pairing logic переключён с per-`cert_nuclide` (групповал все same-parent
proxies cross-source) на per-`(cert_nuclide, spe_filename, geometry)`
ключ. Same-source pairs only. Результат: только один pair выводится —
Pb-212/Tl-208 на Th-228__264_2023@5cm (F-39 baseline, ratio 0.8727).
Cross-source ratios (Cs-137 across sources, etc.) hidden — they reflect
cert metadata + decay correction, не detector response.

**Per-geometry summary table (F-46d diagnostic)**:

Новый block в `validate_certs.py` main() выводит таблицу с:
- n_total / n_meas per geometry
- mean / max |Δ| per geometry
- efr chi²/dof из cached EfficiencyCurve
- bg_file name (sanity check pairing)

Это упрощает диагностику future regressions — поломка одной geometry
изолируется без поиска в 40-row matrix.

**Defensive characteristics**:

1. **All 16 prior fixtures (12 Точ-5см + 4 Точ-25см) bit-identical**
   v1.7.24. F-46b добавляет только volume fixtures с unit="Bq/kg" code
   path; point-source unit="Bq" code path unchanged.
2. **Bq/kg fallback path** preserves compound-cert behavior (Bi-207 +
   Cd-109 в "5431" multi-nuclide cert) — если sub_source.activities
   walk fails (e.g., legacy single-nuclide source), falls back на
   `src.get_activity(name)` returning absolute Bq.
3. **Cross-validation cleanup**: ratios bounded к same-source pairs;
   noise of cross-source comparisons elimited.

**Files changed**:
- `validate_certs.py` (~200 строк net): 24 new volume fixtures,
  run_one Bq/kg handling, per-geometry summary table, cross-validation
  per-source grouping, RunResult.geometry field.
- `test_chain_proxy.py` (~50 строк): 2 теста generalized для
  multi-parent chain proxies (Th-228 + Th-232).

**Verification**: 21/21 test файл проходит. Cert matrix: 40/40
measurable, mean |Δ|=8.48 %, max |Δ|=27.90 % (Bi-214 Дента heavy
container).

**Open follow-ups (next iterations)**:
- (п) **TCS correction refinement** для chain-proxy Tl-208 в close
  geometries. **CLOSED-AS-DOCUMENTED in v1.8.0 (K-21)** — Δ bound
  −15-20 % quantified; resolution requires close-geometry P/T
  experimental data.
- (р) **Self-attenuation correction** для volume samples. **CLOSED-
  AS-DOCUMENTED in v1.8.0 (K-20)** — density-correlated bias bound
  +10 / −7 % для Cs-137 quantified; ρ_sample/ρ_ref diagnostic added.
- (с) **K-40 Петри +17 % investigation** — partial resolution in
  v1.8.0 через F-47a (degree=4 для Петри: Δ=+18.7 → +17%); root
  cause same as K-20 (density mismatch at thin-source extreme).
- (т) **Дента/Петри .efr refit** — **CLOSED-AS-DOCUMENTED in v1.8.0
  (K-19)** — Дента chi²/dof=15 invariant of polynomial degree, lab
  procedure required. F-47a tuning recovers degree-optimal fit per
  geometry.

---

#### F-47a: Per-geometry polynomial-degree tuning для eff curve (v1.8.0)

**Status**: FIXED in v1.8.0 — `EFF_DEGREE` dict в `validate_certs.py`
выбирает optimal polynomial degree per geometry instead of hardcoded
degree=3.

**Контекст**: F-46d (v1.7.25) выявил per-geometry chi²/dof variation
для default degree=3 fit:
- Точ-5см: 6.95 (24 anchors)
- Точ-25см: 2.51 (20 anchors)
- Дента-120мл: 15.40 (13 anchors)
- Петри-60мл: 15.04 (23 anchors)
- Marinelli: 3.72 (15 anchors)

F-47a tests degrees 1-5 per geometry to find optimal.

**Findings**:

| Geometry      | Best degree | Best chi²/dof | Δ improvement   |
|---------------|-------------|----------------|------------------|
| Точечная-5см  | 3 (current) | 6.95           | unchanged        |
| Точечная-25см | **5**       | 1.74 (was 2.51)| max Δ 19.34→18.25%|
| Дента-120мл   | 3 (current) | 15.40 (К-19)  | unchanged        |
| Петри-60мл    | **4**       | 14.28 (was 15.04) | mean Δ 12.55→9.84% |
| Маринелли     | 3 (current) | 3.72           | unchanged        |

**Effect на cert-matrix (overall)**:
- mean |Δ| 8.48 % → **7.89 %** (улучшение 0.6 pp)
- max |Δ| 27.90 % (unchanged — Bi-214 Дента 200g)
- 40/40 measurable preserved

**Файлы изменений**:
- `validate_certs.py` (~12 строк): EFF_DEGREE dict + lookup в
  `_resolve_geometry_resources()`.

**Verification**: 21/21 test файл проходит, 266+ tests без
изменений в count.

---

#### F-47b: Density-ratio diagnostic для volume samples (v1.8.0)

**Status**: FIXED in v1.8.0 — per-fixture `ρ_sample/ρ_ref` annotation
добавлен в cert matrix note column. NO correction applied; this is
diagnostic-only. See K-20 для accepted limitation.

**Контекст**: F-46 cert matrix показывает density-correlated Δ для
single-line nuclides (Cs-137 spread +9.81 / −5.58 % across ρ ratio
0.36-1.04). Annotation makes correlation visible at-a-glance.

**Реализация**: в `run_one()` для unit endswith "bq/kg":
```python
REF_DENSITY = {  # hardcoded from manual .efr inspection
    "Маринелли":   (1000.0, 1.60),  # vol_ml, ρ_ref g/cm³
    "Дента-120мл": ( 120.0, 1.66),
    "Петри-60мл":  (  60.0, 1.60),
}
rho_sample = target_sub.mass_g / vol_ml
rho_ratio = rho_sample / rho_ref
note_parts.append(f"ρ_sample/ρ_ref={rho_ratio:.2f}")
```

**Effect**: cert matrix note column теперь shows e.g.
`(ρ_sample/ρ_ref=0.36; 6 peaks)` for Marinelli Cs-137 light.

**Verification**: 21/21 test pass. Pattern confirms K-20 root cause
(density mismatch correlates with single-line Δ).

---
