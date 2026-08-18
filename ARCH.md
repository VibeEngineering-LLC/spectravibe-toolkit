# ARCH.md — архитектура skill'а SpectraVibe

> **v1.17.9.2 (30.05.2026)** — сводный архитектурный документ.
> Обновлять при изменении публичного API модуля, добавлении новой
> подсистемы, изменении потока данных.
>
> **v1.17.9.2 изменения**: F-155 (все скрипты в `scripts/`),
> F-156 (тесты по step-папкам), TOOLS_INVENTORY обновлён.

---

## 1. Высокоуровневая структура

```
┌──────────────────────────────────────────────────────────────────┐
│  CLI (gamma)  +  PyAPI (gamma.*)  +  RAG (gamma.knowledge.*)     │
└────────────────┬─────────────────────────────────────────────────┘
                 │
        ┌────────▼──────────┐
        │  staged_pipeline  │   ◀── центральный orchestrator (Step 1..11)
        └────┬──────┬──┬────┘
             │      │  │
   ┌─────────┘      │  └──────────────┐
   │                │                 │
   ▼                ▼                 ▼
┌─────┐   ┌──────────────┐   ┌───────────────┐
│ I/O │   │ Calibration  │   │ Peaks / Multi │
└──┬──┘   └──────┬───────┘   └───────┬───────┘
   │            │                   │
   └────┬───────┴─────────┬─────────┘
        ▼                 ▼
  ┌──────────────┐  ┌────────────────┐
  │   Physics    │  │ Identification │
  │ TCS / dt /   │  │ + Residual cls │
  │ self-attn /  │  └──────┬─────────┘
  │ a-priori bg  │         │
  └─────┬────────┘         │
        │                  ▼
        │           ┌──────────────┐
        └──────────▶│   Activity   │
                    │   MDA / dose │
                    └──────┬───────┘
                           ▼
                  ┌──────────────────┐
                  │    Reporting     │
                  │ JSON / MD / HTML │
                  │ + PDF + PNG      │
                  └──────────────────┘
```

---

## 2. Подсистемы (модули `scripts/gamma/*`)

### 2.1. `gamma.io` — Format I/O (15 модулей)

**Назначение**: чтение/конверсия спектральных форматов, разрешение фона.

| Модуль | Контракт | F-rule |
|---|---|---|
| `readers.py` | Dispatch reader по расширению | — |
| `lsrm_spe.py` | LSRM `.spe` parser (KEY=VALUE, ZONES, Peaks) | F-130 |
| `lsrm_spe_text.py` | LSRM text writer | — |
| `atomspectra_xml.py` | AtomSpectra `.xml` parser | — |
| `becqmoni_xml.py` | BecqMoni `.bsp/.xml` parser/writer | F-148 (планируется) |
| `n42_2012.py` | ANSI N42.42-2012 parser | — |
| `convert.py` | Универсальный конвертор форматов | F-149 / F-154 |
| `background.py` | Resolve external bg (link → file → directory) | F-131 |
| `background_search.py` | Auto-discovery фонового кандидата | F-131 |
| `average_lsrm.py` | Усреднение N .spe в один | — |
| `format_registry.py` | Реестр поддерживаемых форматов | См. `FORMAT_REGISTRY.md` |
| `filename_hints.py` | Извлечение hints из имени файла | F-115 |
| `lsrm_efficiency.py` | LSRM `.efr/.efa` (efficiency) reader | — |
| `lsrm_library.py` | LSRM `.lib` nuclide library reader | — |
| `lsrm_src.py` | LSRM `.src` source spec reader | — |

**Контракты I/O**:
- Все readers возвращают `gamma.spectrum.Spectrum` dataclass (см. §2.3).
- Энергетическая калибровка хранится как полином `ENERGY=N,a0,a1,a2,...` с
  опциональными zones (F-FUTURE: F-189 per-zone orthopoly).
- Фоновый файл всегда resolves в `Spectrum.background_embedded`, тег
  `extras["background_source"]` принимает значения `{embedded,
  external_resolved, auto_resolved_from_directory, link_unresolved, none}`.

### 2.2. `gamma.calibration` — Энергия + FWHM + F-145 (15 модулей)

| Модуль | Контракт | F-rule |
|---|---|---|
| `energy_fit.py` | Полиномиальная подгонка E(N), deg≤4 | F-87 |
| `fwhm_fit.py` | FWHM(E) подгонка; **TODO F-168**: FWHM²=a+bE+cE² | F-125, F-168 (план) |
| `fwhm_provider.py` | Callable FWHM(E_keV) → keV | — |
| `bootstrap.py` | F-87 bootstrap refit при anchor disagreement | F-87 |
| `anchor_recalibration.py` | Stored cal sanity-check vs anchors | F-87 |
| `stored_check.py` | Adaptive stored-cal validation | — |
| `seven_line_check.py` | [LSRM-Algo-9] 7-line ЕРН проверка | LSRM-9 |
| `subcalibration.py` | Subzone refit | — |
| `multiplet_self_calibration.py` | F-145 Phase B+C: δ(N) coupling | F-145 |
| `efficiency.py` | ε(E) interpolation | — |
| `efficiency_autoload.py` | Auto-load ε(E) для детектора | — |
| `detector_type.py` | Detect detector type из spec | — |
| `bg_subtract_dual_mode.py` | Background subtraction, two modes (rate_normalized_channel / energy_aligned) + F-243 `ZERO_POINT_MATCH_THRESHOLD_KEV=30` safety gate + `BackgroundConsentRegistry/Required` + full uncertainty propagation. Used by **`scripts/validate_certs.py`** (standalone cert validation). | F-58, F-243 |
| `bg_subtract_energy.py` | Background subtraction, energy-domain `np.interp` (always energy_aligned, no channel-mode trap). Used by **`gamma.identification.staged_pipeline`** (production pipeline / `run_plan_a.py`). | F-58, F-160 |

**F-58 untangle (2026-06-21).** Two distinct `subtract_background` implementations live in this package by design, with separate consumers and separate guarantees:

- **`bg_subtract_dual_mode.py` — FULL.** Standalone cert validation (`scripts/validate_certs.py`). Two explicit modes auto-selected by gain match (`Δa₁/a₁ < 0.5%` → channel; else energy-aligned). F-243 zero-point guard (`|Δa₀| > 30 keV` forces energy-aligned even when gains match). `BackgroundConsentRegistry/Required` for per-file approval. Full σ propagation. Designed for cert/certification paths where the operator has explicitly approved the bg file.
- **`bg_subtract_energy.py` — LITE.** Production pipeline (`gamma.identification.staged_pipeline` → `run_plan_a.py`). Single-function `np.interp` in energy-domain, always energy-aligned. No safety gate (structurally inapplicable — there is no channel-mode path to mis-select). 4 diagnostic fields (Phase 2): `net_uncertainties`, `gain_mismatch_relative`, `zero_point_mismatch_keV`, `n_channels_clipped`. σ-formula corrected 2026-06-21 (k → k² Poisson propagation on scaled bg).

Naming history: the LITE file was `bg_subtraction.py` until Phase 1 (2026-06-21); it was renamed to `bg_subtract_energy.py` to eliminate the graphify label-collision "Background Subtraction" across two communities (c0 + c22). The function name `subtract_background` is shared by both modules; consumers import explicitly from the module that fits their guarantees.

**F-451 direction inversion — done 2026-06-22 (B1+V2: audit-only, no downstream override).** Both modules now scale **DOWN to the shorter live-time** ("к меньшему", operator-fixed direction), not UP to the sample. Three branches in `subtract_background`:

- `t_s > t_bg` → `scale_direction="sample_down"`, `applied_scale = t_bg / t_s`, sample counts scaled down; `effective_live_time = t_bg`.
- `t_s < t_bg` → `scale_direction="bg_down"`, `applied_scale = t_s / t_bg`, bg counts scaled down; `effective_live_time = t_s`.
- `t_s == t_bg` → `scale_direction="equal"`, `applied_scale = 1.0`; `effective_live_time = t_s`.

`BackgroundSubtractionResult` carries new fields `applied_scale`, `scale_direction`, `effective_live_time` (frozen dataclass). Legacy `scale_factor` = `t_s / t_bg` is preserved verbatim for backward compat / JSON schemas. σ-propagation uses `applied_scale` correctly in both directions: `σ²(net) = applied_scale²·N_sample + N_bg` (sample_down) or `σ²(net) = N_sample + applied_scale²·N_bg` (bg_down). The cps-invariant `net_cps = sample_rate − bg_rate` holds identically under all three branches (verified by snapshot tests).

Downstream contract (B1+V2, operator-locked 2026-06-22): `staged_pipeline` **does NOT** rebuild `spec.live_time` after subtraction. `spec.live_time` and `spec.real_time` remain at the original sample values so that all downstream rates / MDA / activity formulas continue to operate at the full sample count statistics (B2 «full breaking downstream» was prototyped, then explicitly reverted — propagating `effective_live_time` into `spec` produced inconsistent MDA across the dual_mode/lite split and broke the `activity = S_net / (t·ε·I)` contract). The F-451 «к меньшему» direction (`applied_scale`, `scale_direction`, `effective_live_time`) lives entirely inside `BackgroundSubtractionResult` for σ-propagation and is mirrored into `spec.extras["background_subtraction_*"]` for audit/diagnostics only. The `notes` string carries `"F-451 scale_direction=…, applied=…, t_s=…, t_bg=…, effective_live_time=…"` for downstream audit. Plan: `audit/_plans/F-451_bg_subtract_direction_invert.md`. Tests: `tests/snapshot/test_bg_subtract_energy.py` (15 inc. 6 F-451-dedicated) + `tests/snapshot/test_bg_subtract_dual_mode.py` (27).

**F-452 — FWHM model degree uplift (2026-06-22, LSRM poly-4 √E).** Legacy NNLS quadratic fit (`FwhmModelKind.legacy_quadratic`, `(a,b,c)` tuple) replaced for LSRM-calibrated detectors by `FwhmModel` callable with `kind="lsrm_poly_sqrt_E"` and 5-coefficient polynomial in √E. Ground truth: `references/lsrm_ground_truth/<detector_base>/fwhm_calibration_lsrm.json` (e.g. Gamma-1S anchors: 60.3 / 122.1 / 661.7 / 1332.5 / 2612.9 keV). `build_fwhm_model(spec)` returns `(FwhmModel, source_tag)` where source becomes `"lsrm_ground_truth_reference_poly4_sqrtE"` when an LSRM JSON is available; falls back to legacy quadratic otherwise. Backward-compat: `fwhm_keV_at_energy(model, e_keV)` and `_make_fwhm_at_channel` accept both `FwhmModel` and 3-tuple; `fwhm_model_legacy_abc(model)` re-fits poly-4 to `(a,b,c)` for the legacy JSON schema. Below-anchor extrapolation is clamped to a 0.1 keV floor (prevents negative FWHM near low-E edge). On Th-232 Marinelli the LSRM poly-4 gives FWHM(2614)=112.8 keV vs legacy quadratic FWHM≈116.9 keV; LSRM cert value is 112.80 keV. Only `build_fwhm_model` call site changed (line 2157); all FWHM consumers are model-agnostic via `fwhm_keV_at_energy`.

**F-452-FU2 — Currie L_C pre-MAD non-detection filter (2026-06-22).** F-452 (more accurate FWHM) **exposed** a latent bug in `compute_activity_for_nuclide` (BUG-38/39 MAD-outlier-rejection at `scripts/gamma/activity/compute.py:978-1038`): when most matched lines of a nuclide get `peak_area_source="deconvolved_coupled"` with numerical-noise values (A_i ~ 1e-21 Bq around zero), the MAD median itself falls into the nano-zone, and real high-count lines (e.g. Tl-208 2614 keV with A_i~2700 Bq) get rejected as 19σ outliers against a nano-consensus. Observed pre-fix: Tl-208 weighted activity collapsed from ~2925 Bq baseline to **9.8e-24 Bq** on Th-232 Marinelli. Fix: pre-MAD filter at `compute.py:906` drops any line with `A_i / σ_A_i < 1.0` (Currie L_C non-detection criterion, Currie 1968 / Gilmore §5.5) **before** MAD-rejection runs, with explicit `lines_skipped` provenance (`below_Currie_LC A_i=… σ=… A/σ=… < 1.0 (pre-MAD non-detection filter, F-452-FU2)`). Verified: Tl-208=2585, Ac-228=2710, Pb-212=2943 Bq; Th-232 chain equilibrium ratio=1.14 (`in_equilibrium=True`) in both NO-BG and WITH-BG paths. Side-effect on snapshot test `test_F389_th232_demo_v2_activity_parity_with_prod`: previously-exact 0.0% prod-vs-V2 parity broke to 15.3% on Ac-228 — this is **disclosure of real V2 path divergence** (V2-only lines 129.06 + 1630.6 keV at Ac-228 were never present in prod; previously masked by symmetric nano-line distortion across both pipelines). Tolerance widened from 5% to 18% with detailed inline justification. Diagnostic: `audit/_drafts/F-452-FU2_prod_vs_v2_diff.py`.

**F-452-FU3 — V2-vs-prod parity on close-multiplets (CLOSED not-a-bug, 2026-06-22).** Initial backlog assumption (расширить `_V2_ONLY_CHANNELS` coverage) was wrong. Investigation (`audit/_drafts/F-452-FU3_v2_only_channels_inspect.py`): the set autogenerates correctly from `peaks_v2 \ peaks_prod` per-channel (snapshot on Th-232 = {177, 251, 318, 378, 429, 778}, 6 channels — channel 539 is NOT in it because both Mariscotti and matched_filter find it). The actual root cause is **peak-search semantic divergence on close multiplets**: at channel 539 (~1588 keV) prod-Mariscotti resolves ONE peak, two Ac-228 library lines compete (1588.20 keV I=6.84% dominant + 1630.63 keV I=1.51% satellite), `shared_peak_dedupe` keeps the dominant, the satellite → `lines_skipped`. V2 dual-method resolves 1588+1630 as TWO separate peaks → both validly enter `lines_used`; 1630.60 (S=10407, A_i=6864 Bq) drags Ac-228 weighted mean down. The 129.06 keV survival in V2 is a derivative effect: with 1630.6 inflating MAD-spread, Ā falls 2740 → 2325 Bq and 129.06 (A_i=1430) is no longer a 12.8σ outlier. V2 super-resolution gives more correct physics; parity on close-multiplets is **impossible by design** (prod dedupes, V2 resolves). `PARITY_TOL_FRACTION = 0.18` retained as honest empirical ceiling; not-a-bug per operator decision (alternatives rejected: applying `shared_peak_dedupe` inside V2 path would kill V2 super-resolution value; extending `_V2_ONLY_CHANNELS` with a super-resolution-multiplet detector is premature optimization without operator signal).

**Контракты калибровки**:
- `Spectrum.energy_cal` = `tuple[float, ...]` (полиномиальные коэффициенты).
- `Spectrum.energy_cal_source` = провенанс (например `"F-145_multiplet_self_calibration"`).
- F-145 4-фазная самокалибровка:
  - Phase A (`peaks/coupled_multiplet.py::coupled_intensity_fit(free_centroids=True)`)
  - Phase B+C (`calibration/multiplet_self_calibration.py::recalibrate_from_multiplet_centroids`)
  - Phase D (`identification/staged_pipeline.py` повторный fit с `centroid_window_frac=0.15`)
- **TODO (F-171)**: после Phase D вызвать `seven_line_check.run_seven_line_check` как gate.

### 2.3. `gamma.spectrum` — Spectrum dataclass

```python
@dataclass(frozen=True, slots=True)
class Spectrum:
    counts: np.ndarray              # shape (n_channels,)
    energy_cal: tuple[float, ...]   # polynomial coefficients
    energy_cal_degree: int
    energy_cal_source: str
    livetime_s: float
    realtime_s: float
    detector_type: str              # "NaI" | "HPGe" | "LaBr3" | ...
    detector_profile: str           # e.g. "Gamma-1C"
    geometry_canonical: str         # "Marinelli" | "Petri" | "Dente" | ""
    sample_mass_kg: Optional[float]
    sample_density_g_cm3: Optional[float]
    sample_volume_cm3: Optional[float]
    fwhm_model: Optional[Callable[[float], float]]
    background_embedded: Optional["Spectrum"]
    background_link: Optional[str]
    extras: dict                    # provenance, hints, F-rule flags
    ...
```

Канонический ENERGY_CEILING_KEV = 3000 keV (sharded ceiling для NaI).

### 2.4. `gamma.peaks` — Поиск и подгонка пиков (6 модулей)

| Модуль | Контракт | F-rule |
|---|---|---|
| `search.py` | Mariscotti 2nd-derivative search | [GILMORE-9.3] |
| `convolution_search.py` | Matched-filter Gaussian search | F-124, F-129 |
| `area.py` | Net peak area: Σ_y − step+linear continuum | — |
| `peak_image.py` | Gauss + low-energy tail + Compton step модель | F-90, F-120, F-126, F-127, F-133 |
| `coupled_multiplet.py` | Intensity-coupled multiplet fit + F-145 Phase A | F-117, F-118, F-145 |
| `deconvolve.py` | Forced-multiplet decomposition chain | F-33, F-34, F-118 |

**Контракты peak shape**:
- Базис: `I(x) = A·G(x;x0,σ) + T·exp((x-x0)/β)·H(x0-x) + h_step·A·0.5·erfc((x-x0)/(σ√2))`
- `h_step` = 0.03 для NaI; **TODO F-183**: 0.003 для HPGe.
- `T(E)`: log-linear модель, slope=0.15, T_ref=0.7 @ 662 keV (NaI).
  **TODO F-190**: per-line T_k свободные.
- `G+tail` C¹-continuous в точке `x0-T·σ`.

**Контракты coupled multiplet**:
- Амплитуды связаны: `a_k = A_nuc · I_γ_k / 100`.
- Phase A `free_centroids=True`: `dE_k ∈ [±0.5·FWHM(E_k)]`.
- Phase D `centroid_window_frac=0.15`: soft-locked после калибровки.

### 2.5. `gamma.physics` — Физический слой (7 модулей)

| Модуль | Контракт | F-rule |
|---|---|---|
| `cascade_summing.py` | TCS `C(E)=1/(1−Σp·ε_T)` ([LSRM-Algo-10]) | F-128, K-17, K-18, F-186 (план) |
| `dead_time.py` | Effective dt `t_d=A·Σy+B·Σ(y·i)` ([LSRM-Algo-15]) | F-95, F-173 (план) |
| `bg_lines_apriori.py` | A-priori bg линии для subtraction | F-96, F-174 (план) |
| `self_attenuation.py` | Marinelli/Petri/Dente: `f=(1−exp(−μρd))/(μρd)` | F-122 |
| `pileup.py` | Pile-up correction (rate-dependent) | — |
| `secondary.py` | Secondary peak detection (Compton edge/BS) | F-141 |
| `secondary_peaks.py` | Selective CE/BS/SE/DE injection | F-141 |

### 2.6. `gamma.identification` — Pipeline + классификация (14 модулей)

| Модуль | Контракт | F-rule |
|---|---|---|
| `staged_pipeline.py` | **Главный orchestrator** Step-1..11 | — |
| `window.py` | Окно идентификации; **TODO F-167**: `±k·FWHM(E)` | F-167 (план) |
| `identify.py` | Базовая identification logic | — |
| `disambiguate.py` | 7-Д шаг (разрешение неоднозначностей) | F-115 |
| `cross_check.py` | Cross-check между линиями | — |
| `proportionality.py` | Intensity proportionality test | — |
| `completeness.py` | Chain completeness (Th-232 / U-238) | — |
| `chain_equilibrium.py` | Secular equilibrium check | F-119 |
| `anchor_ranks.py` | Anchor ranking | — |
| `confidence.py` | Confidence index CI | — |
| `ci_gating.py` | Confidence-based gating | — |
| `ern_set.py` | ЕРН-set (естественные радионуклиды) | LSRM-9 |
| `mda.py` | Currie L_C/L_D; **TODO F-169**: paired-blank | F-169 (план) |
| `residual_classifier.py` | chain_secondary / compton_residual / X-ray | F-143 |

**Pipeline Step-1..11**:
1. Read file + parse metadata + detect bg
2. Determine measurement environment (F-102/F-108)
3. Preliminary peak search (Mariscotti + convolution)
4. Detector type ID (R(662))
5. Bootstrap energy calibration (5α/5β/5γ + bg + F-145)
6. FWHM(E) calibration
7. Identification using isolated lines (7А-7Д)
8. Targeted multiplet deconvolution (F-117/F-118/F-126/F-145)
9. Peak areas + activities + ISO 11929 MDA
10. Secondary peaks + intrinsic detector activity
11. Report bundle (JSON+MD+HTML+PNG+PDF)

### 2.7. `gamma.activity` — Расчёт активностей (3 модуля)

| Модуль | Контракт | F-rule |
|---|---|---|
| `compute.py` | Базовый Bq/Bq·kg расчёт + F-91 σ propagation | F-91 |
| `template_method.py` | [LSRM-Algo-12] sensitivity matrix R | F-100, F-180 (план) |
| `quasitemplate.py` | [LSRM-Algo-13] quasi-template (full-spectrum fit) | F-98 |

### 2.8. `gamma.reporting` — Отчёты (12 модулей)

| Модуль | Контракт | F-rule |
|---|---|---|
| `wrapper.py` | `analyze_and_report` высокоуровневый wrapper | — |
| `build.py` | Сборка отчётного bundle | — |
| `json_report.py` | Канонический JSON отчёт (SKILL_VERSION) | — |
| `markdown_report.py` | Markdown отчёт (RU narrative, F-108) | F-108 |
| `html_report.py` | Базовый HTML | — |
| `interactive_html.py` | Интерактивный HTML с D3.js (F-113 iOS WebView) | F-113, F-134, F-141 |
| `pdf_export.py` | HTML → PDF через headless Edge | F-114 |
| `plots.py` | matplotlib PNG графики | — |
| `chat_summary.py` | 3-8 line chat summary | — |
| `environment.py` | Environment classification (F-102) | F-102 |
| `cost_estimator.py` | F-132 token cost estimate | F-132 |
| `anonymize.py` | F-115 анонимизация (operator/SN forbidden) | F-115 |

### 2.9. `gamma.knowledge` — RAG (3 модуля + 1 init)

| Модуль | Контракт | F-rule |
|---|---|---|
| `rag_extract.py` | PDF → JSON corpus (закреплённый opt-in tool) | F-151, F-154 |
| `rag_index.py` | BM25 индекс builder + serializer | F-151 |
| `rag_search.py` | `rag_query/explain/cite/verify` API + CLI | F-152, F-153 |

Двухслойная архитектура: **manually-curated** `knowledge_index.json`
(52 entries, всегда есть) + **optional corpus** `knowledge_corpus.json`
(собирается явно через `rag_extract`, ~10-20k chunks).

### 2.10. `gamma.detectors` — Detector profiles (1 модуль + 1 каталог)

| Артефакт | Назначение |
|---|---|
| `detectors/gamma1c.py` | Профиль NaI 63×63 Gamma-1C |
| `../detectors/Gamma-1C/` | Эталонные спектры, ε(E), averaged backgrounds, lsrm-libs |

### 2.11. `gamma.data` — Reference data accessors (6 модулей)

| Модуль | Контракт |
|---|---|
| `nuclide_library.py` | Доступ к `data/nuclides.json` |
| `aliases.py` | Алиасы нуклидов (`Cs137 → Cs-137`) |
| `anchors.py` | Anchor patterns для F-87 bootstrap |
| `chain_decomposer.py` | Decomposition Th-232 / U-238 chains |
| `iaea_fetcher.py` | IAEA Live Chart cache fetcher |
| `xrf_catalog.py` | XRF (X-ray fluorescence) линии |

---

## 3. Потоки данных

### 3.1. Базовый flow `gamma analyze --full-report`

```
.spe / .xml input
    ▼
[io/readers.dispatch] ─────────────▶ Spectrum dataclass
    ▼
[io/background.resolve]            (+ background_embedded)
    ▼
[identification/staged_pipeline.analyze_lsrm_spe]
    │
    ├─ Step 1-2: parse + environment
    ├─ Step 3:   peaks (Mariscotti + convolution)
    ├─ Step 4:   detector type ID
    ├─ Step 5:   energy_cal (5α/5β/5γ + bg + F-87 bootstrap)
    ├─ Step 6:   FWHM(E) calibration
    ├─ Step 7:   identification (7А-7Д isolated lines)
    ├─ Step 8:   multiplet deconvolution
    │           │
    │           ├─ F-117/F-118 coupled fit
    │           ├─ F-145 Phase A (free centroids)
    │           ├─ F-145 Phase B+C (recalibrate δ(N))
    │           └─ F-145 Phase D (soft-locked refit)
    │
    ├─ Step 9:   activities (compute + TCS + dt + self-attn) + MDA
    ├─ Step 10:  secondary peaks (CE/BS/SE/DE)
    └─ Step 11:  StagedAnalysisResult
                    │
                    ▼
        [reporting/wrapper.analyze_and_report]
                    │
                    ├─ json_report   ────▶ {stem}_report.json
                    ├─ markdown_report ──▶ {stem}_report.md
                    ├─ interactive_html ─▶ {stem}_report.html
                    ├─ pdf_export    ────▶ {stem}_report.pdf  (опц.)
                    ├─ plots         ────▶ {stem}_spectrum.png, cluster_*.png
                    └─ chat_summary  ────▶ stdout
```

### 3.2. RAG flow

```
references/books/*.pdf
    │
    ├─ (opt-in) rag_extract.py ──▶ knowledge_corpus.json
    │
    └─ (manually curated)
            │
            ▼
    references/knowledge_index.json (52 entries)
            │
            ▼
    rag_index.build_bm25_index ──▶ knowledge_bm25.json
                                       │
                          ┌────────────┼────────────┐
                          ▼            ▼            ▼
                       rag_query  rag_explain  rag_cite/verify
                          │            │            │
                          └────────────┼────────────┘
                                       ▼
                          gamma rag {query, explain, cite, verify, …}
```

---

## 4. Контракты / закреплённые правила

### Permanent contracts (ЗАКРЕПЛЕНО НАВСЕГДА)

| Контракт | Источник |
|---|---|
| F-91 σ = max(σ_weighted, σ_scatter) | [LSRM-Algo-7] (Birge ratio) |
| F-95 dt модель `A·Σy + B·Σ(y·i)` | [LSRM-Algo-15] |
| F-96 a-priori bg peak-list | [LSRM-Algo-9.4] |
| F-98/F-100 template + quasi-template | [LSRM-Algo-12]-13 |
| F-102/F-108 RU narrative + environment | внутренний |
| F-113 iOS WebView template fix | внутренний |
| F-115 анонимизация | внутренний |
| F-117/F-118 intensity-coupled multiplets | [LSRM-Algo-8.4] + extension |
| F-122 Marinelli self-attenuation | [GILMORE-10] |
| F-125 NaI 63×63 FWHM refit | Шендрик pt2 |
| F-126 nonlinear curve_fit | [LSRM-Algo-8.4] + scipy |
| F-127 T(E) per-line | [LSRM-Algo-8.4] + emp |
| F-128 Bi-212 TCS | [LSRM-Algo-10] |
| F-130 auto-density extraction | LSRM .spe format |
| F-131/F-135 auto-bg apply | F-131 + safety |
| F-133 per-line Compton step (ГОСТ) | [LSRM-Algo-8.4.4] + ГОСТ |
| F-138 knowledge library | references/books/ |
| F-141 selective CE/BS/SE/DE injection | [GILMORE-7] |
| F-142 Cs-Kα 32 keV check | nuclide data |
| F-143 chain_secondary / compton_residual | внутренний |
| F-145 4-фазная самокалибровка | [LSRM-Algo-8.4.4] extension |
| F-149 convert_spectrum.py preserved | F-154 contract |
| F-150 PDFs исключены из архива | maintenance |
| F-151..F-153 RAG protocol | knowledge library |
| F-154 tool preservation | maintenance |
| F-155 все Python в `scripts/` | maintenance |
| F-156 тесты по step-папкам | maintenance |

### Maintenance contracts

- **F-150**: Релизный архив исключает `references/books/*.pdf`, transient
  `_test_*/`, `_tmp/`, `__pycache__/`, `.git/`, `.vscode/`, `.idea/`.
- **F-153**: Новые методологические F-rules ОБЯЗАНЫ цитировать запись из
  `knowledge_index.json` через `rag_cite()`.
- **F-154**: Инструменты не удаляются; устаревшие → `LEGACY` в `TOOLS_INVENTORY.md`.
- **F-155**: Все Python-скрипты в `scripts/`; корень — только `*.md`.
- **F-156**: Все тесты в `tests/stepNN_*/` (Step 1..11) или
  `tests/{io,knowledge,smoke,snapshot}/`. Обновлять `tests/INDEX.md`
  при добавлении.
- **handoff.md ↔ handoff_ru.md**: синхронны (RU companion).
- **`KNOWN_AND_FIXED_ISSUES.md`**: обновляется в конце каждой сессии.
- **`INDEX.md`**: обновляется при добавлении/удалении файлов.
- **`ARCH.md`**: обновляется при изменении публичного API модуля.

---

## 5. Зависимости (Python packages)

| Package | Использование |
|---|---|
| `numpy` | Спектры, массивы, polynomial fit |
| `scipy` | `optimize.{least_squares,curve_fit,lsq_linear,nnls}` |
| `matplotlib` | PNG графики |
| `pypdf` | F-151 PDF extraction |
| `pytest` | Регрессионные тесты (9.0.3) |

**Опционально**:
- `reportlab` (планируется F-159 для technical PDF report)

**Не используется**:
- ML-зависимости (sentence-transformers, transformers) — RAG работает на
  чистом BM25 без эмбеддингов.

---

## 6. Точки расширения (будущая архитектура)

См. `ROADMAP_v1_17_8_plus.md`, разделы Tier 1..3:

- **F-167..F-175 (Tier 1)**: HIGH-priority фиксы аудита — точечные правки.
- **F-176..F-188 (Tier 2)**: MED-priority methodological refinements.
- **F-189..F-198 (Tier 3)**: архитектурные:
  - **F-189**: per-zone orthopoly calibration (SpectraLine ZONES support).
  - **F-190**: per-line shape parameters (free σ_k, T_k, h_step_k).
  - **F-193**: X-ray escape classifier для NaI iodine 28 keV.
  - **F-194**: dose-weighted MDA characteristic_line selection.

### Не входит в архитектуру (отвергнуто)

- Полнотекстовые embeddings (sentence-transformers): RAG достаточен через BM25.
- Отдельный database для history (SQLite/Postgres): анализ stateless.
- Cloud-side processing: skill локальный, off-line.

---

## 7. Версионирование

- `SKILL_VERSION` в `scripts/gamma/reporting/json_report.py:21`
- Bump при каждом релизе.
- Snapshot тесты (`test_v1_*.py`) фиксируют historic JSON структуру.
- Release archive: `SpectraVibe_v{N.N.N}.zip` в `../../1_Version/`.

История:
- v1.15.0 — F-86 step pipeline
- v1.16.0..v1.17.0 — refinement
- v1.17.2 — **reference contract** (не трогать references/demo_contract_v1_17_2)
- v1.17.6..v1.17.7 — алгоритмическая refinement, F-126..F-144
- v1.17.8 — F-145 + F-150
- v1.17.9 — F-151..F-154 RAG
- **v1.17.9.1** — RAG-аудит + tests/ + INDEX + ARCH (документация only)
- v1.17.10 — Tier 1 фиксы аудита (план)

---

Подробные F-rule контракты: `KNOWN_AND_FIXED_ISSUES.md`.
Точки расширения: `ROADMAP_v1_17_8_plus.md`.
Каталог всех файлов: `INDEX.md`.
Каталог тестов: `tests/INDEX.md`.
Инструменты: `scripts/TOOLS_INVENTORY.md`.
