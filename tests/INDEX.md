# tests/ — каталог тестов (70 файлов в 14 step/aux папках)

> **v1.17.9.2 (F-156 контракт)**: тесты разложены по этапам pipeline.
> `conftest.py` автоматически настраивает `sys.path` и сохраняет
> flat-alphabetical collection order. **680 / 681 PASS**
> (1 deselect: TD-1 pytest 9 LocalPath).

---

## Структура

```
tests/
├── conftest.py                              PYTHONPATH=scripts + order
├── INDEX.md                                 (этот файл)
│
├── step01_io_and_metadata/        (2)       parse + bg detection
├── step02_environment/            (2)       F-102/F-108 environment
├── step03_peak_search/            (3)       Mariscotti + matched filter
├── step04_detector_type/          (—)       reserved
├── step05_energy_calibration/     (4)       bootstrap + F-145
├── step06_fwhm/                   (1)       F-125 FWHM(E)
├── step07_identification/         (14)      7А-7Д isolated lines
├── step08_multiplets/             (9)       F-117/F-118/F-126/F-145
├── step09_activity_mda/           (8)       compute + MDA + TCS + dt
├── step10_secondary_peaks/        (2)       F-141 CE/BS/SE/DE
├── step11_reporting/              (8)       JSON/MD/HTML/PDF/anonymize
│
├── io/                            (7)       formats, readers, conversion
├── knowledge/                     (1)       F-151..F-154 RAG
├── smoke/                         (4)       Cs-137/K-40/Ra-226 production
└── snapshot/                      (5)       v1.15.x..v1.17.0 fixations
```

---

## Запуск

```bash
# Полная регрессия
PYTHONIOENCODING=utf-8 python -m pytest tests/ -q

# Один этап
PYTHONIOENCODING=utf-8 python -m pytest tests/step08_multiplets/ -v

# Один тест
PYTHONIOENCODING=utf-8 python -m pytest tests/knowledge/test_knowledge_rag.py -v

# Параллельно (если pytest-xdist установлен)
PYTHONIOENCODING=utf-8 python -m pytest tests/ -n auto

# Исключить slow
PYTHONIOENCODING=utf-8 python -m pytest tests/ -m "not slow"
```

---

## Step 1: I/O + Metadata (`tests/step01_io_and_metadata/`)

| Файл | Что покрывает |
|---|---|
| `test_background_subtraction.py` | F-58 energy-rebinned bg subtraction |
| `test_bg_only_environment.py` | Bg-only environment processing |

## Step 2: Environment (`tests/step02_environment/`)

| Файл | Что покрывает |
|---|---|
| `test_bg_apply_and_F136.py` | F-135 default apply + F-136 binding suppression |
| `test_background_auto_search.py` | F-131 auto bg search |

## Step 3: Peak search (`tests/step03_peak_search/`)

| Файл | Что покрывает |
|---|---|
| `test_adaptive_mariscotti.py` | Адаптивный Mariscotti |
| `test_convolution_peak_search.py` | F-124 matched-filter Гауссиан |
| `test_convolution_cli.py` | F-129 CLI `--peak-search-method` |

## Step 4: Detector type (`tests/step04_detector_type/`)

(reserved, нет тестов на v1.17.9.2)

## Step 5: Energy calibration (`tests/step05_energy_calibration/`)

| Файл | Что покрывает |
|---|---|
| `test_stored_check_adaptive.py` | Adaptive stored-cal check |
| `test_subcalibration.py` | Subzone refit |
| `test_self_calibration_F145.py` | F-145 Phase A→B→C→D (5 cases, 1 slow) |
| `test_t_at_e_calibration.py` | F-127 T(E) per-line tail |

## Step 6: FWHM (`tests/step06_fwhm/`)

| Файл | Что покрывает |
|---|---|
| `test_nai_fwhm_refit.py` | F-125 NaI 63×63 refit |

## Step 7: Identification (`tests/step07_identification/`)

| Файл | Что покрывает |
|---|---|
| `test_identification.py` | Базовая identification logic |
| `test_staged_identification.py` | StagedAnalysisResult |
| `test_disambiguation.py` | 7-Д disambiguation |
| `test_filename_binding.py` | F-115 filename binding |
| `test_chain_completeness.py` | Th-232 / U-238 chain completeness |
| `test_chain_equilibrium_guard.py` | F-119 equilibrium guard |
| `test_chain_proxy.py` | Цепь как proxy для родителя |
| `test_k22_chain_equilibrium.py` | K-22 chain equilibrium |
| `test_pb212_in_th_chain.py` | F-123 Pb-212 238 keV |
| `test_th_chain_M1_M2_present.py` | Th-232 M1, M2 |
| `test_th_composite_present.py` | Th composite |
| `test_u238_chain_multiplets.py` | U-238 chain multiplets |
| `test_u238_equilibrium_guard.py` | U-238 equilibrium guard |
| `test_priority_express.py` | Priority express режим |

## Step 8: Multiplets (`tests/step08_multiplets/`)

| Файл | Что покрывает |
|---|---|
| `test_coupled_fit_M1.py` | F-117 Ra-226 M1 |
| `test_coupled_fit_M2.py` | F-117 M2 |
| `test_deconvolve.py` | Forced-multiplet деконволюция |
| `test_nonlinear_curve_fit.py` | F-126 nonlinear curve_fit |
| `test_peak_image_step_F133.py` | F-133 per-line Compton step |
| `test_peak_image_wiring.py` | F-120 peak_image интеграция |
| `test_overlay_arrays_F134.py` | F-134 fit overlay |
| `test_peak_area.py` | Площадь пика |
| `test_round5_pipeline.py` | Round 5 pipeline |

## Step 9: Activity + MDA + Physics (`tests/step09_activity_mda/`)

| Файл | Что покрывает |
|---|---|
| `test_activity.py` | Bq + Bq/kg + σ propagation |
| `test_ac228_activity_recovered.py` | Ac-228 recovery |
| `test_bi212_tcs.py` | F-128 Bi-212 TCS |
| `test_cascade_summing.py` | K-17/K-18 cascade summing |
| `test_tcs_method_scale.py` | TCS scaling |
| `test_self_attenuation.py` | F-122 self-attenuation |
| `test_self_attenuation_marinelli.py` | F-122 Marinelli |
| `test_auto_density.py` | F-130 auto-density |

## Step 10: Secondary peaks (`tests/step10_secondary_peaks/`)

| Файл | Что покрывает |
|---|---|
| `test_secondary_peaks.py` | Вторичные пики |
| `test_secondary_feature_rule.py` | Правило вторичных features |

## Step 11: Reporting (`tests/step11_reporting/`)

| Файл | Что покрывает |
|---|---|
| `test_step11_report.py` | Step 11 финальный отчёт |
| `test_interactive_report.py` | F-113 iOS WebView interactive HTML |
| `test_pdf_artefact.py` | F-114 PDF через Edge |
| `test_cost_footer.py` | Cost footer |
| `test_cost_estimate_per_stage.py` | F-132 per-stage cost |
| `test_rows_sorted_ascending.py` | Sorted table rows |
| `test_no_en_leak.py` | F-108 narrative leak |
| `test_anonymization.py` | F-115 анонимизация |

## I/O (`tests/io/`)

| Файл | Что покрывает |
|---|---|
| `test_average_lsrm.py` | Усреднение LSRM .spe |
| `test_format_conversion.py` | convert_spectrum.py (TD-1 pytest 9 fail) |
| `test_lsrm_spe_extended.py` | LSRM .spe extended fields |
| `test_lsrm_spe_reader.py` | LSRM .spe reader базовый |
| `test_lsrm_src.py` | LSRM source detection |
| `test_reader_api.py` | Reader dispatcher |
| `test_efficiency.py` | Detector efficiency |

## Knowledge / RAG (`tests/knowledge/`)

| Файл | Что покрывает |
|---|---|
| `test_knowledge_rag.py` | F-151..F-154 RAG (13 cases) |

## Smoke (`tests/smoke/`)

| Файл | Что покрывает |
|---|---|
| `test_cs137_smoke.py` | Cs-137 production demo |
| `test_k40_smoke.py` | K-40 production demo |
| `test_ra226_smoke.py` | Ra-226 production demo |
| `test_phase21b_infra.py` | Phase 21b infrastructure |

## Snapshot (`tests/snapshot/`)

| Файл | Что покрывает |
|---|---|
| `test_v1_15_delivery.py` | v1.15 snapshot |
| `test_v1_16_0.py` | v1.16.0 snapshot |
| `test_v1_16_1.py` | v1.16.1 snapshot |
| `test_v1_16_2.py` | v1.16.2 snapshot |
| `test_v1_17_0.py` | v1.17.0 snapshot |

---

## Контракты

- **`conftest.py`** добавляет `scripts/` в `sys.path`, устанавливает
  `PYTHONIOENCODING=utf-8`, регистрирует marker `slow`, сохраняет
  flat-alphabetical collection order через `pytest_collection_modifyitems`.
- **Test paths**: `Path(__file__).parent.parent.parent` → корень skill
  (после раскладки на 3 уровня: tests/stepNN/test_*.py).
- **Slow tests**: `@pytest.mark.slow`, исключаются дефолтом через `--deselect`.
- **TD-1**: `test_format_conversion.py::test_convert_api_explicit_formats`
  pre-existing pytest 9 LocalPath issue.
- **TD-2**: numerical-тесты в step07/step09 order-dependent — collection
  ordering в conftest спасает; полноценный фикс требует teardown.

## F-156 контракт

При добавлении нового теста:
1. Определить этап pipeline (Step 1..11) или aux-категорию.
2. Положить в соответствующую папку.
3. Обновить эту таблицу в `tests/INDEX.md`.
4. Прогнать полную регрессию.
