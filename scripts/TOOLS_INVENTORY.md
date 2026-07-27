# Tools Inventory — gamma-spectrum-analysis

> **F-154 / v1.17.9 (закреплено навсегда)** — канонический реестр всех
> прикладных инструментов skill'а. Перед написанием нового конвертора /
> ридера / экстрактора **ОБЯЗАТЕЛЬНО** проверь этот файл; если подходящий
> инструмент существует — используй его (с доработкой при необходимости),
> а не переписывай с нуля.
>
> Инструменты НИКОГДА не удаляются. Устаревшие помечаются `LEGACY`, но
> файлы остаются как fallback. См. `KNOWN_AND_FIXED_ISSUES.md` §F-154.

---

## Точки входа верхнего уровня

### `scripts/analyze_spectrum.py`
- **Назначение**: legacy CLI обёртка вокруг `python -m gamma.cli analyze`.
- **Статус**: ACTIVE, используется в HTML/PDF отчётах.
- **Запуск**: `python scripts/analyze_spectrum.py <file.spe> [options]`.

### `scripts/build_release_archive.py`  (F-150)
- **Назначение**: канонический упаковщик релиза в `1_Version/`.
- **Контракт исключений**: фиксирован (PDF в `references/books/`,
  `__pycache__`, transient `_test_*`/`_smoke`/`_tmp`).
- **Запуск**: `python scripts/build_release_archive.py 1.17.9`.
- **Результат**: `../../1_Version/gamma-spectrum-analysis_v1.17.9.zip`.

### `scripts/convert_spectrum.py`  (F-149)
- **Назначение**: универсальный конвертор спектров в формат `.spe` LSRM.
  Поддерживает входы: `.xml` AtomSpectra, `.bsp/.txt` BecqMoni, `.csv`,
  `.spc` Canberra, generic 2-column ASCII.
- **Статус**: ACTIVE. ЗАКРЕПЛЁН — переписывать с нуля запрещено.
- **Запуск**: `python scripts/convert_spectrum.py <input> --output <out.spe>`.
- **Примечание**: для новых форматов добавлять дополнительные
  reader-функции в этот же файл, не создавать дубли.

### `scripts/analyze_problem_isotopes.py`  (F-38, перенесён v1.17.9.2)
- **Назначение**: характеризация диапазонов residual + intensity ratio +
  shape descriptors (FWHM, low-tail asymmetry, edge-width) для проблемных
  изотопов (Cs-137 Compton edge ↔ Bi-214 503; Be-7 478; ...).
- **Статус**: ACTIVE, диагностический.
- **Запуск**: `python scripts/analyze_problem_isotopes.py`.

### `scripts/analyze_secondaries.py`  (F-37, перенесён v1.17.9.2)
- **Назначение**: характеризация secondary peaks (CE/BS/SE/DE/xray_escape/
  ic_xray/k40_natural) на чистых Cs-137 и K-40 фикстурах.
- **Статус**: ACTIVE, диагностический.
- **Запуск**: `python scripts/analyze_secondaries.py`.
- **Выход**: `data/secondary_peaks.json`.

### `scripts/build_averaged_backgrounds.py`  (F-43/F-44, перенесён v1.17.9.2)
- **Назначение**: build-time генерация canonical long-exposure backgrounds
  из 2016+2024 Поверка архивов. F-44: cumulative-pattern detection через
  `average_lsrm_spectra::detect_cumulative_pattern`.
- **Статус**: ACTIVE, build-time (run-once при обновлении архива).
- **Запуск**: `PYTHONPATH=scripts python scripts/build_averaged_backgrounds.py [--dry-run]`.
- **Выход**: `data/averaged_backgrounds/`.

### `scripts/sync_knowledge_index_gost_refs.py`  (F-256 + F-154, v1.17.9.4)
- **Назначение**: каноничный (НАВСЕГДА сохраняемый по F-154) sync-инструмент
  для двухслойной схемы ссылок (см. `references/REFERENCES.md` §0). Добавляет
  поле `gost_ref_num` во все entries+books `references/knowledge_index.json`
  для авто-Layer 1 → Layer 2 трансляции в `citation_translator.py`.
- **Статус**: ACTIVE. ИДЕМПОТЕНТЕН — запускать при каждом расширении RAG.
- **Запуск**: `PYTHONIOENCODING=utf-8 python scripts/sync_knowledge_index_gost_refs.py [--dry-run]`.
- **Выход**: in-place patch `references/knowledge_index.json` (14 books + 89 entries в v1.17.9.4).

### `scripts/retrofit_citations_to_layer1.py`  (F-256 + F-154, v1.17.9.4)
- **Назначение**: каноничный регулярочный конвертор legacy ad-hoc ссылок
  (`LSRM §10`, `Gilmore §8.5`, `Будыка §7.5`, `Budyka §2.3`, `Shendrik pt.X гл.Y`)
  в Layer 1 формат `[RAG-ID]` (см. `references/REFERENCES.md` §0).
- **Статус**: ACTIVE. ИДЕМПОТЕНТЕН.
- **Запуск**: `PYTHONIOENCODING=utf-8 python scripts/retrofit_citations_to_layer1.py <files...> [--dry-run]`.
- **История**: в v1.17.9.4 применён к 10 mаster-документам — суммарно 173 замены.

### `scripts/validate_certs.py`  (F-31a/F-31b/F-35, перенесён v1.17.9.2)
- **Назначение**: extended multi-source cert validation на NaI 63×63
  point-source 5cm fixtures. Сравнивает измеренные активности vs
  паспорт (decay-corrected до measurement date).
- **Статус**: ACTIVE, validation.
- **Запуск**: `PYTHONPATH=scripts python scripts/validate_certs.py`.
- **Выход**: `cert_validation_matrix.csv` рядом со скриптом.

---

## RAG knowledge layer (F-151..F-153 / v1.17.9)

### `scripts/gamma/knowledge/rag_extract.py`  (F-151)
- **Назначение**: PDF → JSON корпус с постраничными чанками и
  распознаванием секций. Используется при первичной сборке полного
  корпуса знаний из библиотеки `references/books/`.
- **Статус**: ACTIVE, запускается по явной команде.
- **Запуск**: `python -m gamma.knowledge.rag_extract --books-dir references/books --out references/knowledge_corpus.json`.
- **Время**: ~30-60 с (6 книг, ~770 страниц суммарно).

### `scripts/gamma/knowledge/rag_index.py`  (F-151)
- **Назначение**: сборка BM25 индекса из manually-curated
  `knowledge_index.json` + (опционально) полнотекстового
  `knowledge_corpus.json`. Хранит обратный индекс term → doc_id.
- **Статус**: ACTIVE.
- **Запуск**: `python -m gamma.knowledge.rag_index --rebuild`.

### `scripts/gamma/knowledge/rag_search.py`  (F-152)
- **Назначение**: публичный API для query / explain / cite / verify.
- **Статус**: ACTIVE.
- **Python API**:
  ```python
  from gamma.knowledge import rag_query, rag_cite
  hits = rag_query("Compton step erfc NaI", k=5)
  for h in hits:
      print(h.book, h.section, h.score)
  ```
- **CLI**: `python -m gamma.cli rag query "<вопрос>"`.

---

## I/O layer (преобразование форматов спектров)

### `scripts/gamma/io/readers.py`
- **Назначение**: dispatch reader по расширению (`.spe` → LSRM,
  `.xml` → AtomSpectra, `.bsp` → BecqMoni, прочее).
- **Статус**: ACTIVE.

### `scripts/gamma/io/lsrm_spe.py`
- **Назначение**: парсер `.spe` LSRM (header KEY=VALUE + ZONES +
  channels). Поддерживает orthopoly зональные калибровки, FWHM модели,
  cascades.
- **Статус**: ACTIVE. ЗАКРЕПЛЁН.

### `scripts/gamma/io/atomspectra_xml.py`
- **Назначение**: парсер AtomSpectra Pro XML формата.
- **Статус**: ACTIVE.

### `scripts/gamma/io/background.py`
- **Назначение**: resolve external background по BackgroundSpectrumFile
  link или auto-discovery в типовых местах (F-131).
- **Статус**: ACTIVE.

---

## Reporting / визуализация

### `scripts/gamma/reporting/bilingual_narrator.py`  (F-260, v1.17.9.4)
- **Назначение**: двуязычный narrative enricher (ru → ru+en) на базе словаря
  Будыка-2021 (`data/glossary_budyka_2021.json`, 255 терминов, источник [13]).
  При ПЕРВОМ упоминании русского термина в Markdown/HTML автоматически
  добавляет английский эквивалент в скобках: «Активность» → «Активность (Activity)».
- **Pipeline**: `enrich_text(content) → translate_text(content) → save`.
  JSON-отчёты НЕ обогащаются.
- **Особенности**: поддержка аббревиатур `(ППП)` отдельной формой;
  protection от вставки в code-blocks/links/already-enriched fragments.
- **Doctests**: 3/3 pass.
- **CLI**: `python -m gamma.reporting.bilingual_narrator <file> [--in-place|--stats]`.

### `scripts/gamma/reporting/citation_translator.py`  (F-256, v1.17.9.4)
- **Назначение**: автоматическая трансляция Layer 1 ссылок `[RAG-ID]` →
  Layer 2 (ГОСТ Р 7.0.5–2008) `[№, локатор]` для user-facing отчётов
  (HTML/PDF/MD). Контракт: каждый `reporting/*.py`, генерирующий
  пользовательский вывод, обязан вызывать `translate_text(content)` перед
  сохранением. JSON-отчёты остаются на Layer 1.
- **Статус**: ACTIVE.
- **Public API**: `translate_rag_id(rag_id)`, `translate_text(text)`, `translate_file(path)`.
- **CLI**: `python -m gamma.reporting.citation_translator <file> [--in-place|--strict]`.
- **Источник истины маппинга**: `references/REFERENCES.md` §0 (20 prefix→№ записей).
- **Doctests**: 5 проходят.

### `scripts/gamma/reporting/json_report.py`
- **Назначение**: канонический JSON-отчёт результата анализа.
- **Статус**: ACTIVE.

### `scripts/gamma/reporting/markdown_report.py`
- **Назначение**: Markdown-отчёт (RU-narrative).
- **Статус**: ACTIVE.

### `scripts/gamma/reporting/interactive_html.py`
- **Назначение**: интерактивный HTML с D3.js-графиками.
- **Статус**: ACTIVE.

### `scripts/gamma/reporting/plots.py`
- **Назначение**: matplotlib PNG графики спектра / multiplet'ов.
- **Статус**: ACTIVE.

### `scripts/gamma/reporting/pdf_report.py`  (F-114)
- **Назначение**: HTML → PDF через headless Edge (Chromium).
- **Статус**: ACTIVE.

---

## Анализ / физика / калибровка

### `scripts/gamma/peaks/coupled_multiplet.py`  (F-117/F-118/F-145)
- **Назначение**: intensity-coupled multiplet fit + F-145 Phase A.
- **Статус**: ACTIVE.

### `scripts/gamma/peaks/peak_image.py`  (F-90/F-120/F-126/F-127/F-133)
- **Назначение**: peak-image модель (Gaussian + low-energy tail +
  per-line Compton step). ГОСТ.
- **Статус**: ACTIVE.

### `scripts/gamma/peaks/deconvolve.py`
- **Назначение**: forced-multiplet deconvolution chain.
- **Статус**: ACTIVE.

### `scripts/gamma/peaks/convolution_search.py`  (F-124/F-129)
- **Назначение**: matched-filter поиск пиков (Gilmore §9.3).
- **Статус**: ACTIVE.

### `scripts/gamma/calibration/multiplet_self_calibration.py`  (F-145)
- **Назначение**: Phase B+C — δ(N) коррекция энергетической шкалы.
- **Статус**: ACTIVE.

### `scripts/gamma/physics/cascade_summing.py`  (K-17/K-18/F-128)
- **Назначение**: TCS coincidence summing correction (LSRM §10).
- **Статус**: ACTIVE.

### `scripts/gamma/physics/dead_time.py`  (F-95)
- **Назначение**: эффективное мёртвое время `t_d = A·Σy + B·Σ(y·i)`.
- **Статус**: ACTIVE.

### `scripts/gamma/physics/bg_lines_apriori.py`  (F-96)
- **Назначение**: a-priori фоновые линии (K-40, Tl-208, Bi-214, …).
- **Статус**: ACTIVE.

### `scripts/gamma/activity/template_method.py`  (F-100)
- **Назначение**: template-метод (sensitivity matrix R).
- **Статус**: ACTIVE.

### `scripts/gamma/activity/quasitemplate.py`  (F-98)
- **Назначение**: квази-template метод (full-spectrum simultaneous fit).
- **Статус**: ACTIVE.

### `scripts/gamma/activity/mda.py`
- **Назначение**: ISO 11929 MDA, Currie L_C/L_D.
- **Статус**: ACTIVE.

---

## Identification pipeline

### `scripts/gamma/identification/staged_pipeline.py`
- **Назначение**: главный orchestrator Step-1..11.
- **Статус**: ACTIVE.

### `scripts/gamma/identification/residual_classifier.py`  (F-143)
- **Назначение**: классификация residual peaks (chain / Compton / binding).
- **Статус**: ACTIVE.

---

---

## Activity accuracy slice — v1.17.20 (F-294..F-298)

### `scripts/gamma/activity/self_absorption.py`  (F-294 / T-002+T-024, v1.17.20)
- **Назначение**: Cutshall-style analytic self-absorption correction для Marinelli (low-E ±15% при ρ_sample ≠ ρ_calibration). Встроенная NIST XCOM water μ/ρ table (30..3000 keV).
- **Public API**: `SelfAbsorptionInputs`, `self_absorption_factor()`, `mu_over_rho_water()`, `correct_activity_for_self_absorption()`, `batch_self_absorption_factors()`.
- **Статус**: ACTIVE. Stdlib-only. Integration: pending v1.18.1 (`compute.py` wire-up).

### `scripts/gamma/activity/pt_ratio_nai.py`  (F-295 / T-011, v1.17.20)
- **Назначение**: Peak-to-total ratio P/T(E) для NaI по Gilmore Table 8.4 (3"×3" + 4"×4" + diameter linear interp). Используется для `total_efficiency_from_fep` (нужно TCS).
- **Public API**: `GILMORE_TABLE_8_4_3IN3`, `GILMORE_TABLE_8_4_4IN4`, `pt_ratio_nai()`, `pt_ratio_for_detector()`, `total_efficiency_from_fep()`.
- **Статус**: ACTIVE. Stdlib-only. Integration: pending v1.18.1.

### `scripts/gamma/activity/tcs_close_geometry.py`  (F-296 / T-007+T-082, v1.17.20)
- **Назначение**: TCS coincidence-summing correction `C_TCS = 1/(1−Σ p·ε_T)` (Gilmore §8.6). Preset cascades для Co-60, Eu-152, Cs-137 (no cascade), Ba-133.
- **Public API**: `CascadeLine`, `CascadePair`, `CASCADE_PRESETS`, `compute_tcs_correction()`, `tcs_correction_for_detector()`.
- **Статус**: ACTIVE. Stdlib-only. Integration: pending v1.18.1.

### `scripts/gamma/activity/matrix_method_chi2.py`  (F-297 / T-027, v1.17.20)
- **Назначение**: Stdlib WLS multi-nuclide deconvolution через Gauss-Jordan inversion. Acceptance χ²_red ≤ 1.5 OK, >3.0 missing nuclide.
- **Public API**: `PeakObservation`, `NuclideContribution`, `MatrixMethodResult`, `solve_matrix_method()`.
- **Статус**: ACTIVE. Stdlib-only. Integration: pending v1.18.2 (`--method=matrix` opt-in).

### `scripts/gamma/activity/bg_lines_builder.py`  (F-298 / T-013, v1.17.20)
- **Назначение**: F-96 canonical bg library (27 lines: K-40, Th-232/U-238/U-235 chains, Bi-207, Co-60 contamination) → F-131 deconvolution input. `classify_chain_dominance` для NORM-классификации.
- **Public API**: `F96_BG_LIBRARY`, `filter_bg_lines_in_window()`, `build_f131_input()`, `get_anchor_candidates()`, `classify_chain_dominance()`.
- **Статус**: ACTIVE. Stdlib-only. Integration: pending v1.18.2.

## Peak-image .cpt layer — v1.17.21 (F-299..F-301)

### `scripts/gamma/peaks/peak_image_tabulated.py`  (F-299 / T-021a, v1.17.21)
- **Назначение**: Anchor-based tabulated peak shape (LSRM §8.4.1). `PeakShapeAnchor(E_keV, fwhm_keV, tail_fraction, tail_slope_inv_keV, step_height_frac, asymmetry, weight)`. Preset-defaults по detector_class (NaI: tail=0.03/step=0.05; HPGe: 0.001/0.005).
- **Public API**: `PeakShapeAnchor`, `TabulatedPeakImage`, `build_anchors_from_calibration()`.
- **Статус**: ACTIVE. Stdlib-only. Integration: pending v1.18.3 (`peak_image.py` shape source).

### `scripts/gamma/peaks/peak_image_logspline.py`  (F-300 / T-021b, v1.17.21)
- **Назначение**: Log-log interp FWHM/tail/step + linear asymmetry (LSRM §8.4.3). Bracket-search между anchors; `was_extrapolated` flag для E вне диапазона.
- **Public API**: `InterpolatedPeakShape`, `interpolate_peak_shape()`, `fwhm_at_E()`, `batch_interpolate()`.
- **Статус**: ACTIVE. Stdlib-only. Integration: pending v1.18.3.

### `scripts/gamma/io/cpt_io.py`  (F-301 / T-021c, v1.17.21)
- **Назначение**: LSRM .cpt XML I/O (Calibrated Peak Template). Round-trip preserves anchors+metadata. Forward-compat: unknown tags warn (non-strict) или raise (strict).
- **Public API**: `CPT_SCHEMA_VERSION`, `build_cpt_xml()`, `parse_cpt_xml()`, `write_cpt_file()`, `read_cpt_file()`.
- **Статус**: ACTIVE. Stdlib-only. Integration: pending v1.18.3 (CLI import/export tool).

## Quasi-template (ЛСРМ §13) — v1.18.0 (F-302..F-304)

### `scripts/gamma/activity/quasi_template_ppp.py`  (F-302 / T-022a, v1.18.0)
- **Назначение**: Per-nuclide PPP (Peak-Plus-Pedestal) sum-spectrum builder. Erf-analytic gaussian integration по channel edges (math.erf, без numpy). Continuum injection через F-304 factory. **Additive** — старый `quasitemplate.py` (numpy) сохранён.
- **Public API**: `NuclideLine`, `NuclideDef`, `PPPTemplate`, `build_nuclide_template()`, `build_templates_for_library()`, `validate_template_collection()`.
- **Статус**: ACTIVE. Stdlib-only. Integration: pending v1.18.4 (`--solver=quasi-template` CLI).

### `scripts/gamma/activity/quasi_template_fit.py`  (F-303 / T-022b, v1.18.0)
- **Назначение**: Full-spectrum WLS simultaneous fit (LSRM §13). Gauss-Jordan normal equations (как F-297 matrix_method). Poisson weights w=1/max(obs,1). χ²_red≤1.5 acceptance.
- **Public API**: `QuasiTemplateFitResult`, `solve_quasi_template_fit()`.
- **Статус**: ACTIVE. Stdlib-only. Integration: pending v1.18.4.

### `scripts/gamma/activity/compton_continuum.py`  (F-304 / T-022c, v1.18.0)
- **Назначение**: Analytic Klein-Nishina Compton continuum для NaI. `E_C = 2E²/(m_e c² + 2E)`. Erf-smoothed step + backscatter peak (10% area) + plateau. Area-preserving.
- **Public API**: `ELECTRON_REST_MASS_KEV`, `compton_edge_keV()`, `backscatter_peak_keV()`, `compton_continuum_for_line()`, `make_continuum_func()`.
- **Статус**: ACTIVE. Stdlib-only. Integration: pending v1.18.4 (drop-in F-302 pedestal).

## Production QA gates — v1.17.18 (F-289..F-292)

### `scripts/gamma/calibration/spectrometer_compliance.py`  (F-289 / T-083, v1.17.18)
- **Назначение**: Gate соответствия спектрометра паспорту (FWHM, channels, ИНЛ, drift). PASS/WARNING/FAIL per field. Default spec для Gamma-1S NaI 63×63.
- **Статус**: ACTIVE. Integration: ready (вызывается из QA-skripta при необходимости).

### `scripts/gamma/activity/sample_gates.py`  (F-290 / T-084, v1.17.18)
- **Назначение**: Gate геометрии и плотности образца vs validity-range калибровки. 4-level status; `requires_self_absorption_correction` flag для extrap-зоны.
- **Статус**: ACTIVE. Integration: pending v1.18.1 (вместе с self_absorption).

### `scripts/gamma/calibration/bg_drift.py`  (F-291 / T-054, v1.17.18)
- **Назначение**: F-test variance ratio + z-score на mean для bg-ROI. Pre-computed F-critical (no scipy hard dep).
- **Статус**: ACTIVE. Integration: ready (periodic-QA hook).

### `scripts/gamma/calibration/sensitivity_drift.py`  (F-292 / T-053, v1.17.18)
- **Назначение**: Quarterly ε(E)-drift check vs reference. Per-line threshold (2%/5% per LSRM §14.4) + calendar window + Pearson r для monotonic shift.
- **Статус**: ACTIVE. Integration: ready (periodic-QA hook).

---

## Books library tooling — v1.17.19 (F-293)

### `scripts/verify_books_inventory.py`  (F-293, v1.17.19)
- **Назначение**: Sanity-check `books_library/` inventory vs `references/books/INDEX.md`. Warning (не блок) при MISMATCH. Env override `GAMMA_BOOKS_LIBRARY_DIR`.
- **Запуск**: `python scripts/verify_books_inventory.py` или `--strict-books-inventory` через build_release_archive.

### `scripts/build_books_archive.py`  (F-293, v1.17.19)
- **Назначение**: Standalone packager `books_library/` → `1_Version/books_library/gamma-books_vYYYY-MM-DD.zip` с MANIFEST.sha256 и rotation (.prev.zip).
- **Запуск**: `python scripts/build_books_archive.py`.

---

## Политика добавления нового инструмента

При добавлении нового конвертора / ридера / экстрактора / визуализатора:

1. Прочитать этот файл — возможно подходящий уже существует.
2. Если нужен новый — добавить запись СЮДА в той же коммите (имя файла,
   назначение, статус, запуск, F-rule если применимо).
3. Зарегистрировать в `KNOWN_AND_FIXED_ISSUES.md` как новый F-rule.
4. Если новый инструмент дублирует функцию старого — пометить старый
   `LEGACY` (но не удалять).

См. F-154 в `KNOWN_AND_FIXED_ISSUES.md` для полного контракта.
