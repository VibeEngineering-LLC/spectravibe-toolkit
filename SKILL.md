---
name: spectravibe
description: "Comprehensive analysis of gamma-ray spectrum files. Two-pass workflow: (1) metadata, environment detection, bootstrap calibration, background recognition; (2) strict identification from isolated lines + intensity ratios + decay chains, targeted multiplet deconvolution for confirmed nuclides, secondary peak classification, elemental XRF identification (Kα/Kβ patterns Z 13-95) with XRF-gamma disambiguation, intrinsic detector activity, dead-time correction, confidence index, completeness DC, MDA per ISO 11929. Use whenever user uploads .spe/.chn/.n42/.mca/.txt/.csv/.spm files, mentions gamma/HPGe/NaI/LaBr3/CeBr3/CZT spectra, or asks to identify radionuclides, calibrate, find peaks, deconvolve multiplets, compute MDA, identify elements from XRF, or diagnose shielding. Methods: Gilmore & Joss 3rd Ed., Lsrm Algorithmic Foundations."
---

# SpectraVibe

A full workflow for analyzing a gamma-ray spectrum from an arbitrary file. Methodology is grounded in Gilmore & Joss, *Practical Gamma-ray Spectrometry*, 3rd Ed. (Wiley, 2024), with several refinements adopted from Lsrm *Algorithmic Foundations* (2025): library-directed peak search, two-pass identification-first → multiplet deconvolution, dead-time correction model, Currie/ISO 11929 detection limits, intensity-ratio χ² check, confidence index (CI), and Dose-Contribution completeness metric (DC). Nuclear data from ENSDF/NuDat 3 (NNDC) and LNHB recommended values.

## Scope

**Energy range of routine interest: 0–3000 keV.** Rationale: ²⁰⁸Tl 2614.5 is the highest gamma line of routine interest from natural background; above 3 MeV mostly hosts overflow markers, pile-up, and cosmic-ray secondaries.

**Since v1.18.32 (BUG-9, 2026-06-03) the readers DO NOT trim at 3 MeV by default.** All decoded channels are returned. The previous default silently dropped the high-E tail — most painfully for low-channel-count NaI files where a 1024-channel Gamma-1C spectrum with a slightly negative `a0` (calibration drifted left) would shrink to 1003 channels at read time, surprising users who saw their channel count shrink during conversion. The 3 MeV ceiling is now opt-in.

*Per-call override*: both readers (`read_atomspectra_xml`, `read_lsrm_spe`) accept keyword-only `apply_energy_ceiling: bool = False` (set True to drop channels whose energy exceeds the ceiling) and `ceiling_keV: float | None = None` (use a different ceiling for a single call without changing the `ENERGY_CEILING_KEV` constant; ignored when `apply_energy_ceiling=False`). The same kwargs pass through `gamma.io.readers.read_spectrum`. For analysis-stage trimming there is also `gamma.spectrum.trim_to_working_energy(spec, max_keV=3000)`, which performs the same cut after the read so the action is visible at the call site instead of hidden in the reader.

**Not in scope (explicitly NOT computed):**
- Efficiency calibration ε(E) for a specific sample-detector geometry
- Sample mass / volume / density (treated as nominal placeholders even when present in file)
- Source-detector distance and geometry-dependent corrections
- Specific activity (Bq/kg or Bq/L)

Consequence: MDA is reported in **counts per second** (cps), not in becquerels. Per-nuclide reporting uses **peak count rate (cps)** instead of activity. Decay correction to a preparation date remains available because count rate decays with the same time constant.

**Supported file formats**: `.spe`, `.chn`, `.n42`, `.mca`, `.txt`/`.csv`, `.spm` (Lsrm multi-section binary), and **`.xml`** (AtomSpectra / Lsrm `ResultDataFile`). AtomSpectra XML is the primary tested format and supports embedded background spectrum (`<BackgroundEnergySpectrum>` block) with its own independent calibration.

**Active spectrometric complex (v1.12.0): `Gamma-1C`** — УДС-ГЦ-63×63 NaI(Tl) detector head by Aspect + Lsrm SpectraLine DAQ software (canonical alias resolves Колибри-1М, Гамма-1С, БДЭГ-63×63-USB and related vendor tokens). All detector-specific assets — `.efr` efficiency curves, `.spe` reference spectra, source passports, LSRM nuclide libraries, the empirical secondary-peak catalogues, averaged backgrounds, intrinsic-activity reference and dead-time A/B coefficients — live under **`detectors/Gamma-1C/`** and are reached exclusively through `gamma.detectors.gamma1c` (path resolver, single source of truth). Per the isolation policy adopted 2026-05-29, other complexes (AtomSpectra / AtomNano / RadiaCode) are recognized at the alias level only and will gain their own `detectors/<canonical>/` subtree + path resolver after the Gamma-1C pipeline stabilizes; algorithms in `gamma.peaks`, `gamma.identification`, `gamma.calibration`, `gamma.activity`, `gamma.physics` will be replicated rather than shared once that work begins.

**Round-5 quantification flags (v1.13.0)**: `analyze_lsrm_spe(...)` now accepts `apply_deconvolution: bool=False` (run multiplet deconvolution on identified clusters — canonical use case: the 600–680 keV Bi-214 / Cs-137 / Cs-134 cluster on Th-rich samples), `compute_activities: bool=False` (per-nuclide activity in Bq via `compute_activities_for_all` with cascade-summing dispatched per detected nuclide and K-21 close-geometry P/T scaling), `sample_mass_kg: float=None` (emits specific activity in Bq/kg), `compute_mda: bool=False` (per-line ISO 11929 / Lsrm §6.3 detection limits for every detected line plus the standard ЕРН + technogenic suite). All flags default OFF — backward-compatible with v1.12.0.

**Step 11 report + complete-workflow umbrella (v1.14.0)**: convenience flag `analyze_lsrm_spe(complete_workflow=True)` autonomously enables Steps 8/9/10 (deconvolution + activities + MDA, with `deconvolution_overlap_fwhm` bumped to 3·FWHM for the 600–680 keV cluster). New `gamma.reporting` package builds Step-11 artefacts from the resulting `StagedAnalysisResult`: `gamma.reporting.build_report(result, *, output_dir=..., write_json=True, write_markdown=False)` writes `{stem}_report.json` (machine-readable schema 0.1) and optionally `{stem}_report.md` (full 13-section long-form), and returns a 3–8 line chat summary. `gamma.reporting.classify_environment(result)` returns `"natural"` / `"low_background"` / `"unknown"` (Step 2). Intrinsic-activity signature for Gamma-1C NaI 63×63 is emitted in `diagnostics.intrinsic_activity_signature` — only I K-escape at E−28.6 / E−32.3 keV; no ¹³⁸La / ²²⁷Ac / Ge or CdTe K-escape (sanity check distinguishing NaI from other scintillators / semiconductors).

## Token economy rules

**Hard budget**: one full spectrum analysis must fit in **60k tokens** (30% of the 200k Pro context window). Includes SKILL.md + loaded reference modules + spectrum metadata + AI reasoning + report. The counts array itself is NEVER in the AI context — it lives only in Python memory.

> Production cost target (per-spectrum analysis): 20k Claude output-tokens
> alarm threshold for algorithmic-optimization review. Scope: production-run
> only (dev = MAXIMUM compute). See `~/.claude/skills/spectravibe-dev/SKILL.md` §«Production cost target».

### Where to spend tokens (high value)

The AI is the **arbiter** for:
- Ambiguous identifications near the rejection threshold (score borderline, χ²/ν between strict and lenient gates)
- Conflicting signals (filename says X, spectrum suggests Y; metadata says HPGe, FWHM says NaI)
- Mechanism assignment for XRF peaks when more than one is plausible (shield vs IC vs matrix)
- Composing the 3-8 line in-chat summary at the end
- Answering follow-up questions about an already-analyzed spectrum

### Where NOT to spend tokens (deterministic, do in Python)

The AI must NOT think through:
- Polynomial fits, residual checks, χ² calculations — `gamma.calibration.*`
- Peak search (Mariscotti), Currie L_C significance — `gamma.peaks.search`
- Library lookups by energy, intensity-ratio χ², DPR equilibrium status — `gamma.identification.*`
- CI and DC formula evaluation — `gamma.identification.{ci,dc}`
- Secondary-peak classification by energy and width — `gamma.physics.secondary`
- XRF Kα/Kβ doublet matching — `gamma.physics.xrf`
- L_C, L_D, MDA arithmetic — `gamma.physics.mda`
- JSON report assembly — `gamma.reporting.json_report`

The AI consumes the output of these modules (a small structured result, not the inputs and intermediates).

### Reference modules — lazy load only when needed

These references are **not** loaded by default. Read each one ONLY when its trigger condition is met:

| Module | Tokens | Load when |
|--------|--------|-----------|
| `01_metadata_calibration.md` | ~3.5k | Reading a new file format, or bootstrap calibration triggered |
| `02_peak_search_deconvolution.md` | ~2.5k | Step 8 (multiplet deconvolution) triggered |
| `03a_identification_algorithm.md` | ~2k | Step 7 (identification) starting — **frequent** |
| `03b_nuclide_library.md` | ~3k | Verifying a specific nuclide's full line set, OR DPR chain check |
| `03c_iso11929_mda.md` | ~0.5k | Step 9 (MDA) |
| `04_secondary_peaks.md` | ~2.5k | Step 10 secondary classification |
| `05_intrinsic_detector_activity.md` | ~2k | Detector is scintillator (LaBr₃, NaI, CeBr₃) OR has unexplained low-E peaks |
| `06_report_format.md` | ~1.5k | Step 11 report assembly |
| `07_dead_time_correction.md` | ~1.5k | Dead time > 5% OR pile-up signatures present |
| `08a_xrf_principles.md` | ~1.5k | Peaks in 5–110 keV range AND XRF possible |
| `08b_xrf_lines_catalog.md` | ~2.5k | Verifying a specific element's K/L lines |
| `09_xrf_gamma_interference.md` | ~3k | Ambiguous coincidence (e.g., 59.5 keV Am vs W) |

**Default load** at workflow start: only SKILL.md (this file) + `01` (for metadata) + `03a` (for identification algorithm). Other modules pulled in as needed.

### Output discipline

- **Default report format is JSON**, not Markdown. Write `report.json` to disk and return a 3–8 line summary to chat.
- **Full Markdown report** only on explicit request (`--full-report` or "give me a full / detailed report").
- **Never echo `counts` array** to chat. The numpy array stays in Python.
- **Never paste large library extracts** to chat. Reference by nuclide name + line energy.
- **One analysis, one pipeline run.** Follow-up questions about an already-analyzed spectrum are answered from `report.json`, not by re-running. Re-run only on explicit user request or new input file.

### Spectrum counts handling

The counts array is the dominant data structure (8k channels × 8 bytes ≈ 64 KB = ~16k tokens raw). It is never serialized into AI context. Only summary statistics (total, peak channel, range) are exposed. Anyone modifying `Spectrum.to_summary_dict()` must keep `counts` excluded.

## Working principles

1. **Autonomous staging.** Decide which steps to run, which to skip, and which to repeat — based on what the spectrum itself shows. Do not blindly execute all 11 steps. Concrete examples of valid skipping:
 - If the file's stored energy calibration agrees with bootstrap recovery within 0.3·FWHM at every anchor, do not refit; reuse stored coefficients.
 - If FWHM(E) from the file matches measured isolated-peak widths within ~5%, do not refit.
 - If dead time < 5% and no peaks at 2·E_strong are visible, skip pile-up analysis.
 - If detector is HPGe and there are no peaks below 120 keV, skip Ge K-escape search.
 - If no Pb/Sn/Cd/Cu/W K-X-ray signatures appear, skip shielding-diagnostic block.

 Before starting, briefly state which steps are planned and which are deferred/skipped based on the file's initial scan. The user can override; otherwise proceed.

2. **Resource economy.** Use the cheapest correct method first. Concrete rules:
 - **Display smoothing**: never on data points; OK to apply minimal cosmetic smoothing to *overlaid model curves* if the user explicitly asks for cleaner plots. Internal Savitzky–Golay for second-derivative search is allowed and does not affect displayed data.
 - **Polynomial degrees in calibrations are capped at 4** (energy, FWHM, efficiency). Higher orders fit noise. If residuals at degree 4 are still poor, segment the energy range and fit piecewise polynomials with continuity at boundaries rather than going to degree 5+.
 - **Multiplet deconvolution**: do not blindly deconvolve every multi-peak region. Only deconvolve multiplets whose components correspond to **already-identified** nuclides (see step 7 — the strict order is identification first).
 - **Library-directed peak search**: only for the working candidate list, not the entire NuDat database.

3. **Trust the spectrum, not the metadata.** Filename, embedded text, stored calibration are *priors*. Independent confirmation from the spectrum is required for every claim.

4. **Display conventions.**
 - Y-axis: counts per second (cps), log scale.
 - No smoothing on displayed data points.
 - Title contains filename + date/time + live/real time + sample ID.
 - Primary FEPs and secondary peaks in distinct colors with explicit labels.

5. **When information is insufficient, say so.** Do not guess.

6. **Local-first delegation (MANDATORY).** Per the user's global Local-First policy (~/.claude/CLAUDE.md), project AGENTS.md §5, and the LOCKED-2026-06-04 model policy (`audit/_drafts/_ollama_helpers/_context/ollama_models_2026-06-04.md`). Delegate to local Ollama (`http://127.0.0.1:11434`) BEFORE consuming Claude tokens on any of these. **Single generative default: `qwen3-coder:30b`** for all extraction / classification / templating / summarization. **No fallback chains** — fallback is itself a source of instability. Embeddings: `bge-m3:latest` only.
 - **Reading raw multi-MB dumps**: PDFs of references, multi-megabyte `.spe`/`.spm` text dumps, NuDat exports, long `.csv` peak lists — pre-process through `qwen3-coder:30b` (summarize ≤200 lines markdown, or `format='json'` for structured extraction). Claude reads the digest, not the source.
 - **Bulk classification**: hundreds of detected peaks to filter as primary/secondary/noise; large nuclide-line tables to filter by intensity — `qwen3-coder:30b` with `format='json'` (qwen3:4b disabled 2026-06-04 after 10/10 batch failure on O-2).
 - **Templated report sections**: per-nuclide table rows, per-multiplet summary blocks, standard prose ("Identification confidence: ..."), error/warning lists — `qwen3-coder:30b` generates skeleton; Claude finalizes the 3-8 line chat summary and arbitrates ambiguous cases (per §"Where to spend tokens").
 - **Long-document semantic search**: locate which page of Gilmore § discusses a topic, which ISO 11929 paragraph covers a formula — `bge-m3:latest` embeddings + cosine similarity. Don't grep through 600-page PDFs in Claude context.
 - Helper-scripts: `audit/_drafts/_ollama_helpers/` (project-level) or `scripts/ollama/` (per-skill-instance).
 - **NEVER delegate to Ollama**: Edit/Write/Bash on project files, git operations, final anti-hallucination check (every claim ↔ offset/line), user dialog, deciding "do or don't". Ollama is a text workhorse, not an agent.

## Strict order of operations

The workflow is a **two-pass identification-first scheme**:

```
Pass 1 (calibration & survey):
  Step 1  — Read file, parse metadata, detect background spectrum if present
  Step 2  — Determine measurement environment (low-background vs natural)
  Step 3  — Preliminary peak search in channel space
  Step 4  — Preliminary detector-type identification
  Step 5  — Bootstrap energy calibration (only if stored is unreliable)
  Step 6  — FWHM(E) calibration (only if needed)

Pass 2 (identification & analysis):
  Step 7  — Identification using ONLY isolated characteristic lines + intensity ratios + DPR
  Step 8  — Targeted multiplet deconvolution ONLY for confirmed nuclides
  Step 9  — Peak areas, count rates, MDA per ISO 11929
  Step 10 — Secondary-peak classification and intrinsic detector activity
  Step 11 — Report with confidence index (CI), completeness (DC), diagnostics
```

### Per-step instructions

The schema above is binding for ordering. Detailed instructions for each
step — sub-step breakdowns, decision rules, formulas, F-rule references,
nuclide tables — live in **`references/STEPS.md`**. Read that file ONCE
when you start an analysis (or per-step on demand) before executing
the step. The detail file mirrors the schema 1:1.

Quick map (with the reference files each step depends on):

| Step | Topic | Detail in `references/` |
|---|---|---|
| 1 | Read file, parse metadata, detect background | `STEPS.md` §Step 1 + `01_metadata_calibration.md` |
| 2 | Determine measurement environment | `STEPS.md` §Step 2 |
| 3 | Preliminary peak search in channel space | `STEPS.md` §Step 3 + `02_peak_search_deconvolution.md` |
| 4 | Preliminary detector-type identification | `STEPS.md` §Step 4 |
| 5 | Bootstrap energy calibration (conditional) | `STEPS.md` §Step 5 + `01_metadata_calibration.md` |
| 6 | FWHM(E) calibration (conditional) | `STEPS.md` §Step 6 |
| 7 | Identification from isolated lines + ratios + DPR + CI | `STEPS.md` §Step 7 + `03a_identification_algorithm.md` + `03b_nuclide_library.md` |
| 8 | Targeted multiplet deconvolution | `STEPS.md` §Step 8 + `02_peak_search_deconvolution.md` |
| 9 | Peak areas, count rates, MDA per ISO 11929 | `STEPS.md` §Step 9 + `03c_iso11929_mda.md` |
| 10 | Secondary peaks + intrinsic detector activity + XRF | `STEPS.md` §Step 10 + `04_secondary_peaks.md` + `08a/08b_xrf_*.md` + `09_xrf_gamma_interference.md` |
| 11 | Report with CI / DC / diagnostics + F-115 anonymisation | `STEPS.md` §Step 11 + `06_report_format.md` |

Critical contract reminders that apply across steps and which the agent
must NEVER skip without explicit user opt-out:

- **F-131 / F-135 — background subtraction.** When a background spectrum
  is available (paired file, auto-search hit, or `background_path`
  argument), it is ALWAYS subtracted before activity / MDA. Silent
  fall-through to gross-spectrum analysis is forbidden.
- **F-115 — anonymisation.** Every emitted JSON / MD / HTML / PDF passes
  through `anonymize_report_inplace()` before disk write. Operator /
  serial / absolute-path leakage is a contract violation.
- **F-256 / F-260 / F-319 — user-facing language.** Reports must use the
  two-layer citation scheme and bilingual narrator pipeline (
  `_f317_apply_user_facing_compliance` in `scripts/gamma/reporting/build.py`).
- **F-318 / F-324 / F-327 — references list.** «Список использованной
  литературы» is injected INSIDE the `.page` container; baseline
  refs {2, 7, 12, 19, 24} always included (ref 1 removed per F-337.4/v1.18.19.1).
- **F-326 — passport activity comparison.** Section is rendered for all
  non-background-only spectra; when `passport_activity_Bq` not provided,
  the deferred-state message + Python instruction are shown. For
  `measurement_environment == "background_only"` spectra the section is
  suppressed (F-UX-04 / 2026-06-04).

**File may contain a background spectrum** in addition to the sample spectrum. Specifically: Lsrm `.spm` format is a binary container with multiple spectra (each a section in standard Lsrm format), some marked as background measurements. N42 files may carry multiple `<RadMeasurement>` blocks with `measurementClassCode="Background"`. If a background spectrum is detected:
- Extract it separately.
- **Its energy calibration may differ from the sample spectrum's** even if the same detector — temperature drift, gain shift, time gap between measurements. Do not assume calibration consistency.
- Apply steps 5–6 (calibration recovery) independently to the background spectrum.
- Use the recalibrated background spectrum for subtraction or for background-peak listing as per step 10's continuum/peak handling.

**Display the sample spectrum**: cps, log scale, no smoothing on data, title with metadata. If background is present, display it as a second trace in the same plot with reduced opacity for visual comparison.

Read `references/01_metadata_calibration.md` for filename token parsing, format-specific metadata fields, date/time usage, and background-spectrum handling.

## When to ask the user vs decide autonomously

Decide autonomously:
- Whether to rebuild calibration (use residual test)
- Polynomial degree (start low, increase as needed)
- Which steps to skip (per the rules in §"Working principles")
- Which multiplets to deconvolve (after step 7 identifies what's in them)

Ask the user only when:
- No anchor peaks found (pure source, no natural background) AND no filename hint
- Filename strongly suggests nuclide X but spectrum shows nuclide Y (confirm interpretation: mislabeling vs unexpected source)
- The detector type cannot be determined unambiguously from the spectrum
- Dead time > 30% and A, B coefficients are unknown for this detector

## Source of methods

- Gilmore G., Joss D. *Practical Gamma-ray Spectrometry*, 3rd Ed., Wiley, 2024 — primary methodology (Ch. 6, 7, 8, 9).
- Lsrm Algorithmic Foundations (2025) — confidence index CI, completeness DC, dead-time correction, two-pass identification scheme, library-directed search, intensity-ratio χ², peak stripping, ISO 11929 implementation.
- ENSDF / NuDat 3 (NNDC, BNL), IAEA Live Chart of Nuclides — nuclear data.
- LNHB Recommended Data — precision I_γ values.
- NIST X-Ray Transition Energies Database — XRF lines.
- Currie L. A. (1968); ISO 11929 — L_C, L_D, MDA.
- Quarati et al., NIM A 729 (2013); Cámara et al., Appl. Radiat. Isot. 109 (2016) — LaBr₃ intrinsic.

---

## Knowledge library / RAG protocol  (F-151..F-153 / v1.17.9)

### Library at a glance

`references/books/` содержит 6 PDF (≈30 МБ, исключены из релизного архива
по F-150, описаны в `INDEX.md` + curated `references/knowledge_index.json`):

1. **LSRM Algorithmic Foundations 2022** — §7 σ, §8.4 peak shape, §8.4.4 Compton step, §9 cal, §9.4 a-priori bg, §9.7 continuum, §10 TCS, §11 ISO 11929, §12 template, §13 quasi-template, §14.2 dose, §15 dead time.
2. **LSRM .spe format spec** — KEY=VALUE header, ZONES, MATERIAL/Ro.
3. **Gilmore 2008 Practical γ-ray Spectrometry** — Ch.2 decay, Ch.4 detectors, Ch.5 efficiency + P/T, Ch.6 background, Ch.7 interactions, §8.5 TCS, §9.3 peak search, §9.7 continuum, Ch.10 geometry/self-attn, Ch.11 cal sources.
4. **Шендрик «Сцинтилляторы» pt.1** — физика NaI(Tl), peak shape.
5. **Шендрик «Сцинтилляторы» pt.2** — FWHM(E) модель, новые сцинтилляторы.
6. **«Анализ и представление результатов эксперимента»** — статистика, χ², weighted mean, ISO uncertainties, Currie MDA.

### RAG-протокол при принятии методологических решений (F-153)

ПЕРЕД тем как принять или предложить методологическое решение
(новая модель, изменение алгоритма, выбор параметра, изменение
границы применимости) — **обязательно** выполнить один из четырёх
RAG-паттернов:

```python
from gamma.knowledge import rag_query, rag_explain, rag_cite, rag_verify

# Pattern 1 — ASK: «что говорит библиотека по теме»
hits = rag_query("Compton step erfc NaI", k=5)

# Pattern 2 — EXPLAIN: связный ответ + цитата + формула
exp = rag_explain("каскадное суммирование TCS Co-60")

# Pattern 3 — CITE: каноническая цитата для F-rule docstring/отчёта
cite = rag_cite("Marinelli self-attenuation")  # → [Gilmore, Ch.10, p.335-355]

# Pattern 4 — VERIFY: guard ПЕРЕД сильным утверждением
verdict = rag_verify("h_step ≈ 0.03 для NaI peak shape")
if not verdict.supported:
    raise ValueError(f"методология не обоснована: {verdict.reason}")
```

**CLI shortcut** для интерактивной проверки в чате:
```bash
PYTHONPATH=scripts python -m gamma.cli rag explain "тема"
PYTHONPATH=scripts python -m gamma.cli rag cite "тема"
PYTHONPATH=scripts python -m gamma.cli rag verify "утверждение"
```

### Когда RAG обязателен

- **Новый методологический F-rule** (изменяет физическую модель / алгоритм / границы применимости) → `rag_cite()` в docstring И запись в `references/knowledge_index.json`.
- **Изменение параметра** (h_step, FWHM коэффициенты, окно идентификации) → `rag_verify()` показать что выбранное значение из источника.
- **Спор «как правильно»** (peak shape, continuum модель, MDA правило) → `rag_explain()` чтобы дать пользователю цитату.

### Когда RAG НЕ обязателен (но желателен)

- Maintenance / safety / process F-rules (как F-150, F-154, F-153).
- Чисто инфраструктурные правки (CLI, форматирование отчётов).
- Уже однократно процитированный источник, на который ссылается F-rule, есть в `knowledge_index.json` → можно ссылаться по `doc_id`.

### Расширение библиотеки

При добавлении нового источника:
1. Положить PDF в `references/books/`.
2. Обновить `references/books/INDEX.md` (биб. метаданные + SHA).
3. Добавить ≥1 запись в `references/knowledge_index.json`.
4. `python -m gamma.cli rag rebuild` для пересборки BM25.
5. Опционально — `python -m gamma.knowledge.rag_extract` для full-text корпуса.

Подробный контракт см. `KNOWN_AND_FIXED_ISSUES.md` §F-151..F-153.

---

## Citation scheme — двухслойная (F-256 / F-260 / F-265, v1.17.9.4)

### F-256 — формат ссылок на первоисточники

**Layer 1 (внутренний, FAST)** — формат `[RAG-ID]` для **моей работы**:
- Python docstrings, F-rule обоснования
- `KNOWN_AND_FIXED_ISSUES.md`, `AUDIT_*.md`, `ROADMAP*.md`
- JSON-отчёты, поле `code_citations`
- Audit-tables, RAG-search keys

Пример: `[LSRM-Algo-10]`, `[BUDYKA-7.5]`, `[LSRM-ACT-11]`, `[GILMORE-9.7]`.

**Layer 2 (внешний, ГОСТ Р 7.0.5–2008)** — `[№, локатор]` + затекстовая
библиография для **пользователя**:
- Чат-отчёты (сводные таблицы, итоговые report-сообщения в диалоге)
- HTML / PDF / MD regulated outputs `reporting/*.py`
- Печатная нормоконтрольная документация

Пример: `[5, Прил.5]`, `[7, §10]`, `[12, §7.5]` + список литературы.

**Автотрансляция**: `scripts/gamma/reporting/citation_translator.py` — функция
`translate_text(text)` заменяет Layer 1 → Layer 2 через mapping
prefix→GOST-№ в `references/REFERENCES.md` §0.

**Источник истины** библиографии: `references/REFERENCES.md` (21 запись).

### F-260 — двуязычный narrative enricher (ru/en)

`scripts/gamma/reporting/bilingual_narrator.py` — при ПЕРВОМ упоминании русского
термина в Markdown/HTML автоматически добавляет английский эквивалент в скобках
из словаря Будыка-2021 (`data/glossary_budyka_2021.json`, 255 терминов).

**Pipeline** для user-facing отчётов:
```python
text = enrich_text(raw_narrative)           # F-260 ru+en
text, _ = translate_text(text)              # F-256 Layer 1 → Layer 2
save(text)
```

### F-265 — глоссарий как pre-context (контракт для AI-агента)

ПЕРЕД написанием/изменением методологического кода (новый F-rule,
docstring, отчётный narrative) — **предварительно** свериться с
`data/glossary_budyka_2021.json` для канонической терминологии:

1. **Терминологический lookup**: «как Будыка называет X?» — не выдумывай
   синонимы, используй канон.
2. **Gap-discovery**: если термин из словаря НЕ покрыт кодом/F-rules —
   кандидат на расширение roadmap (например, «Болометр» стр.28,
   «Времяпролётный спектрометр» стр.31).
3. **Двуязычный API-naming**: имена функций и docstrings должны
   соответствовать паре `ru ↔ en` из словаря (избегаем диалектных
   вариантов «эффективность регистрации» vs «эффективность детектора»).

**Лимит**: глоссарий — справочник, не догма. Если ЛСРМ-источник [7] или
МВИ-Активность [5] даёт более точный термин — применяется F-157 (ЛСРМ-приоритет).

---

## Repository layout (F-155 / F-156, v1.17.9.2)

### F-155 — все Python-скрипты в `scripts/`

Корень skill'а содержит **только `*.md`** документы + каталоги. Любой
новый Python-файл (CLI, диагностика, build-time, validation) создаётся
в `scripts/` и регистрируется в `scripts/TOOLS_INVENTORY.md`. См.
`KNOWN_AND_FIXED_ISSUES.md` §F-155.

### F-156 — тесты по step-папкам

Все 70 тестов разложены по этапам pipeline:

```
tests/
├── conftest.py              (sys.path → scripts/, collection order)
├── INDEX.md
├── step01_io_and_metadata/         (Step 1: parse + bg detect)
├── step02_environment/             (Step 2: F-102/F-108)
├── step03_peak_search/             (Step 3: Mariscotti / matched filter)
├── step04_detector_type/           (Step 4: reserved)
├── step05_energy_calibration/      (Step 5: bootstrap + F-145)
├── step06_fwhm/                    (Step 6: FWHM(E))
├── step07_identification/          (Step 7: 7А-7Д identify)
├── step08_multiplets/              (Step 8: F-117/F-118/F-126/F-145)
├── step09_activity_mda/            (Step 9: compute + MDA + TCS)
├── step10_secondary_peaks/         (Step 10: F-141 CE/BS/SE/DE)
├── step11_reporting/               (Step 11: JSON/MD/HTML/PDF)
│
└── (aux) io/, knowledge/, smoke/, snapshot/
```

При добавлении нового теста:
1. Определить, какой Step pipeline он покрывает.
2. Положить в соответствующую `tests/stepNN_*/` папку.
3. Обновить `tests/INDEX.md` синхронно.

См. `KNOWN_AND_FIXED_ISSUES.md` §F-156.

### Контракт обновления документов в каждом релизе

При каждом релизе ОБЯЗАТЕЛЬНО обновляются:
- `KNOWN_AND_FIXED_ISSUES.md` — новые F-rules.
- `ROADMAP_v1_17_8_plus.md` — закрыть выполненные пункты, добавить новые.
- `handoff.md` + `handoff_ru.md` — release notes (синхронно).
- `INDEX.md` — каталог изменившихся артефактов.
- `ARCH.md` — изменения архитектуры (если есть).
- `TOOLS_INVENTORY.md` — новые/обновлённые скрипты.
- `tests/INDEX.md` — изменения раскладки тестов.
- `SKILL_VERSION` в `scripts/gamma/reporting/json_report.py`.
