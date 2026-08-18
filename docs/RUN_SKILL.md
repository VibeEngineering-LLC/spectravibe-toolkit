# `run_skill.py` — оператор-гайд

> F-398 / v1.18.28. End-to-end оркестратор skill pipeline без AI babysitting.
> Один CLI/Python entry-point: сырой `.spe` → полный отчётный bundle.
> Спецификация: `AGENTS.md` §5.

---

## Что делает

Один запуск:

1. **Pre-flight** — валидация входов, демо-root, версия скилла.
2. **Load & calibrate** — чтение спектра, проверка stored calibration.
3. **PROD analyze** — `gamma.reporting.analyze_and_report(...)` (полный pipeline Step-1..11).
4. **PROD artefacts** — гейт presence (`*_report.json/md/html/.pdf`, plots, BecqMoni XML).
5. **V2 analyze** *(опц)* — `gamma.experimental.v2_integration.analyze_and_report_v2(...)`.
6. **V2 artefacts** *(опц)* — гейт V2-выхода.
7. **Compare** *(опц)* — `gen_v2_compare_th232.py` поверх обоих JSON.
8. **Index** — рендер `bundle/index.html` с навигацией.
9. **Finalize** — `run_skill_summary.json` со всеми результатами фаз.

Каждая фаза независима, идемпотентна, и помечается маркером `bundle/.phases/phase_N.done`. Прерывание не теряет прогресс.

---

## Установка / окружение

```bash
# Windows + Git Bash (рекомендуется)
export GAMMA_DEMO_REPORTS_DIR="<WORKDIR>/demo_reports"
export PYTHONIOENCODING=utf-8
export PYTHONPATH=scripts                # для `from gamma...` импортов

cd "0_Work/gamma-spectrum-analysis"

# Python 3.14 (вотfix Project-state); pip install -r scripts/requirements.txt
```

**Зачем `PYTHONIOENCODING=utf-8`:** Windows cp1251 консоль крашится на ` →` символах в `print()` диагностики.

**`GAMMA_DEMO_REPORTS_DIR`:** контракт F-384 — bundle лежит ВНЕ репозитория. Если не задано — `run_skill.py` берёт `~/demo_reports/` или путь из `gamma.data.demo_reports_root.ensure_demo_reports_root()`.

---

## Сценарии

### Минимум — один спектр

```bash
python scripts/run_skill.py path/to/sample.spe
```

Что произойдёт:
- bundle → `$GAMMA_DEMO_REPORTS_DIR/sample/`
- background — auto-detect (F-131/F-135 `background_auto=apply`)
- mass — выводится из filename token (`*kg`/`*г`) или geometry default
- production-only (без V2)
- генерация всех артефактов: JSON/MD/HTML/Technical PDF/PNG plots/BecqMoni XML

### С явными параметрами

```bash
python scripts/run_skill.py path/to/Th232.spe \
  --background path/to/Фон.spe \
  --mass 0.5 \
  --output-dir "$GAMMA_DEMO_REPORTS_DIR/v1_18_28_th232" \
  --include-v2
```

`--include-v2` добавляет:
- `bundle/sample_v2/` — V2 dual-method peak search
- `bundle/v2_compare/v2_compare_report.html` — 2-column сравнение PROD vs V2

### Батч-режим

```bash
python scripts/run_skill.py \
  --batch "evals/fixtures/*.spe" \
  --output-dir "$GAMMA_DEMO_REPORTS_DIR/auto_$(date +%Y%m%d)"
```

Каждый файл получает свой sub-bundle. Финальный `manifest.csv` в root содержит `stem,bundle,exit_code,phases_ok,phases_failed,elapsed_total_s`.

### Возобновление после краша

```bash
python scripts/run_skill.py --resume "$GAMMA_DEMO_REPORTS_DIR/v1_18_28_th232"
```

Фазы с маркером `.phases/phase_N.done` пропускаются (`status=resumed`). Безопасно дёргать сколько угодно раз — повторный успешный запуск идемпотентен.

---

## Флаги CLI

| Флаг | Назначение |
|---|---|
| `spectrum` *(positional)* | Путь к спектру. Не нужен с `--batch` / `--resume`. |
| `--background <path>` | Явный фон. Иначе auto-resolve (см. background priority chain ниже). |
| `--auto-detect-bg` | Дополнительно искать sibling-файл по filename-эвристике (`фон`/`bkg`/`закр`). |
| `--mass <kg>` | Масса образца. Иначе filename token / geometry default. |
| `--output-dir <path>` | Целевая директория bundle. По умолчанию `$GAMMA_DEMO_REPORTS_DIR/<stem>/`. |
| `--include-v2` | Запустить V2 ветку + compare. |
| `--config <path>` | YAML/JSON config с override defaults. |
| `--resume <bundle>` | Возобновить прерванный прогон. |
| `--batch <glob\|dir>` | Батч-режим. |
| `--no-pdf` / `--no-html` / `--no-plots` / `--no-markdown` / `--no-xml` | Отключить отдельный артефакт. |
| `--quiet` / `-q` | Без stderr (файл лога пишется всегда). |
| `--verbose` / `-v` | DEBUG в stderr. |
| `--version` | Версия скилла и выход. |

---

## Конфигурация (`config.yaml` или `.json`)

Все секции опциональны — отсутствующие ключи берутся из defaults. YAML требует `pip install pyyaml`; JSON работает из коробки.

```yaml
# config.yaml — пример override
analyze:
  sample_mass_kg: 0.5
  sample_density_g_cm3: 1.6
  export_becqmoni: both         # off | sample | bg | both
  full_report: true
  complete_workflow: true
  background_auto: apply        # off | suggest | apply
  background_auto_max_days: 90
  peak_search_method: mariscotti
  filter_narrow_peaks: null     # null = auto per detector (F-316)
  narrow_peak_fwhm_ratio: 0.3
  allow_stage2: true
  allow_stage3: true

multiplet:
  # F-387.1 NaI calibration — defaults owned by Agent A, mirror only.
  unresolved_separation_fwhm_factor: 1.1
  max_components_per_cluster: 3
  min_significance_snr: 3.0           # F-391
  min_significance_snr_singleton: 5.0 # F-391

v2:
  enabled: false                 # CLI --include-v2 = override true
  compare: true

artefacts:
  json: true
  markdown: true
  html: true
  technical_pdf: true
  plots: true
  xml_bq: true

output:
  base_dir: null                 # null = derive from env
  subdir_sample: sample
  subdir_sample_v2: sample_v2
  subdir_compare: v2_compare
```

---

## Background priority chain

`run_skill.py` решает источник фона детерминированно. Первый матч выигрывает:

1. **Explicit** `--background <path>` или `analyze.background_path` в config → используется как есть.
2. **`--auto-detect-bg`** — sibling-файл (`*фон*`/`*bkg*`/`*закр*`) в той же директории; берётся ближайший по mtime.
3. **F-397.1 embedded extraction** — если у спектра-файла есть встроенный фоновый блок:
   • AtomSpectra `<BackgroundEnergySpectrum>`
   • N42-2012 `<RadMeasurement measurementClassCode="Background">`
   В Phase 1 `run_skill.py` извлекает этот блок в `<bundle>/.embedded_bg/<stem>_embedded_bg.spe` (LSRM формат) и подаёт его в pipeline как обычный `--background-path`. F-397 (bg peak detection + HTML toggle «Фон») отрабатывает штатно — без изменений в `staged_pipeline`.
4. **F-131/F-135 pipeline auto-resolve** (`background_auto=apply`) — если предыдущие источники не сработали, сам `analyze_lsrm_spe` ищет совместимый фон в стандартных папках (`data/averaged_backgrounds`, `*Фон*/`); фильтры — тот же детектор, совместимая геометрия, `|Δt| ≤ 90 дн`.
5. **Без фона** — пайплайн пишет `background_status="absent_no_subtraction"`.

Источник фиксируется в `phase_1.detail.background_source`:
- `external` — flag 1 или 2 сработал
- `embedded_extracted` — flag 3 сработал (с путём к извлечённому файлу в `background_extracted_path`)
- `none` — flag 4 решит сам pipeline (F-131); status увидите в `sample/*_report.json`

## Bundle layout

```
<bundle>/
├── sample/                                    # Phase 3 PROD
│   ├── <stem>_report.json
│   ├── <stem>_report.md
│   ├── <stem>_report.html
│   ├── <stem>_technical_report.pdf
│   ├── <stem>_plots/spectrum.png
│   ├── <stem>_plots/multiplets/multiplet_*.png
│   ├── <stem>_calibrated.bq.xml
│   └── <bg_stem>_calibrated.bq.xml
├── sample_v2/                                 # Phase 5 V2 (opt)
│   └── ... (та же структура что sample/)
├── v2_compare/                                # Phase 6 (opt)
│   ├── compare_data.json
│   └── v2_compare_report.html
├── index.html                                 # Phase 7 навигация
├── run_skill_summary.json                     # Phase 8 манифест
├── run_skill.log                              # полный лог прогона
├── .embedded_bg/                              # F-397.1 (если применимо)
│   └── <stem>_embedded_bg.spe                 # извлечённый из spec.background_embedded
└── .phases/
    ├── phase_0.done
    ├── phase_1.done
    └── ...
```

---

## Exit codes

| Код | Условие |
|---|---|
| 0 | Все фазы успешны (или штатно skipped). |
| 1 | Partial — одна или несколько non-fatal фаз упали; bundle частично пригоден. |
| 2 | Fatal — Phase 0 не прошла (входы плохие, env сломано). |

Не fatal:
- V2 analyze упал — Phase 5/6 пропускаются, PROD bundle остаётся валидным.
- Compare упал — index.html всё равно сгенерируется.
- Technical PDF падает в редких случаях reportlab — JSON/HTML/MD остаются.

Fatal только Phase 0 — отсутствие спектра, недоступный demo_root, malformed `--resume` маркер.

---

## `run_skill_summary.json` — что внутри

```json
{
  "skill_version": "v1.18.28",
  "spectrum": "...",
  "background": "...",
  "bundle_base": "...",
  "include_v2": false,
  "started_at_unix": 1748880000.0,
  "elapsed_total_s": 87.4,
  "phases": {
    "0": {"phase": 0, "name": "preflight", "status": "ok",
          "elapsed_s": 0.05, "detail": {...}},
    "1": {...},
    "...": {...},
    "8": {...}
  },
  "metadata": {"stem": "...", "detector_hint": "...", "geometry_hint": "...",
               "distance_cm": 0, "sample_mass_kg": 0.5},
  "config_resolved": {...}
}
```

Этот файл — single-source-of-truth прогона. Удобно грепать в батч-режиме:

```bash
for d in batch_out/*/; do
  jq '.phases | to_entries | map(select(.value.status=="failed")) | length' \
    "$d/run_skill_summary.json"
done
```

---

## Защитные контракты

- **F-384** `demo_reports` ВНЕ скилла → bundle не пишется в репо если есть `GAMMA_DEMO_REPORTS_DIR`.
- **F-386** «пик вылета», не «ускользание» — index.html соблюдает.
- **F-372** IAEA hard-lock — nuclides library не модифицируется orchestrator-ом.
- **F-387.1** defaults `factor=1.1`, `max_K=3` — отражены в config, but **не меняются** без physics-обоснования (см. AGENTS.md §3.1).
- **No-data-loss** — каждая фаза в своём try/except; частичный bundle сохраняется при крахе.
- **No-unsanctioned** — orchestrator НЕ трогает physics/identification/calibration. Только параметризация + reporting hooks.

---

## Troubleshooting

| Симптом | Причина | Фикс |
|---|---|---|
| `UnicodeEncodeError: 'charmap' codec` | Windows cp1251 console + `→` в diagnostic strings | `export PYTHONIOENCODING=utf-8` |
| `ModuleNotFoundError: No module named 'gamma'` | scripts/ не в PYTHONPATH | `export PYTHONPATH=scripts` |
| `ERROR: spectrum not found` | rel path резолвится от CWD | передавать absolute path или `cd` в repo root |
| Phase 6 `gen_v2_compare exit=1` | sample/ или sample_v2/ нет нужного JSON | проверить Phase 2/4 завершились (`.phases/phase_2.done`, `phase_4.done`) |
| `--resume bundle not found` | bundle создан в другой кодировке cwd | передавать absolute path |
| `sample_mass_kg не задан…` warning | filename не содержит mass token | передать `--mass` явно или `analyze.sample_mass_kg` в config |
| Все фазы `resumed` после крашa | штатное поведение `--resume` | пересоздать bundle с нуля, если нужно — `rm -rf bundle/.phases/` |

---

## Тесты

```bash
# Fast unit tests (всегда быстро, <1s)
pytest tests/snapshot/test_run_skill_orchestration.py -m "not slow"

# Slow integration (full pipeline, ~30s)
pytest tests/snapshot/test_run_skill_orchestration.py -m "slow"

# Полный комплект
pytest tests/snapshot/test_run_skill_orchestration.py
```

---

## Связанные документы

- `AGENTS.md` — distribution Math vs Reports, §5 спецификация этого скрипта.
- `HANDOFF_v1_18_*.md` — текущая verаsion snapshot.
- `SKILL.md` — методология skill (Pass 1/2, autonomous staging).
- `scripts/gamma/cli.py` — низкоуровневый CLI, который оборачивает `run_skill.py`.
