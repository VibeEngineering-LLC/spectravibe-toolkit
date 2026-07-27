# README_RESEARCH_PLAN — bilingual latte-style landing pages

> Phase-1 deliverable (research + plan only). No code, README.md, or `*.py`
> were changed to produce this file. Anti-hallucination: every non-trivial
> claim cites a concrete file:line. Uncertain items are marked **TODO/Open**,
> not invented.

---

## Часть A — Находки (ответы 1–9)

### 1. Что это за проект (1–2 предложения, для оператора)

**SpectraVibe** — инструмент анализа файлов гамма-спектров: он читает спектр с
NaI/HPGe-детектора, находит пики, идентифицирует радионуклиды, считает
активность и формирует отчёт (JSON / Markdown / HTML / PDF).
Источники: `CLAUDE.md:6-8` («Gamma-spectroscopy methodology skill: peak search,
identification, calibration, deconvolution, activity calculation»),
`README.md:85` («Skill for the comprehensive analysis of gamma-ray spectrum files»).

### 2. Целевой пользователь (оператор vs разработчик)

Два контура, явно разделены в `CLAUDE.md:10-12`:
- **Оператор** — запускает анализ спектра, **не** грузит dev-skill.
- **Разработчик** — грузит `spectravibe-dev` skill (`~/.claude/skills/spectravibe-dev/SKILL.md`)
  с правилами разработки/релиза (`CLAUDE.md:10-12`, mirror `docs/spectravibe-dev_skill_mirror_2026-06-06.md:14`).

Лендинги ориентируем в первую очередь на **оператора** (точка входа), с короткой
секцией «для разработчиков» со ссылками.

### 3. Ключевые особенности (Features) — только реальные, с источником

1. **Мульти-форматный ввод/вывод + N-to-N конвертер.** 4 нативных формата:
   LSRM SpectraLine `.spe` (binary + ASCII), AtomSpectra/BecqMoni XML,
   ANSI/IEEE N42.42-2012. — `FORMAT_REGISTRY.md:20-27`; утилита
   `scripts/convert_spectrum.py`.
2. **Полный конвейер анализа Step-1..11** (поиск пиков → идентификация →
   деконволюция → активность → MDA) одним вызовом `analyze_and_report(...)`. —
   `scripts/gamma/cli.py:11-15, 79-199`; changelog `README.md:69-73` (v1.14.0/v1.15.0).
3. **Поиск пиков двумя методами**: Mariscotti (2-я производная) и matched-filter
   свёртка, плюс режим `compare`. — `scripts/gamma/cli.py:310-317`;
   `README.md:170` (Phase 1.3 DONE).
4. **Стадийная идентификация нуклидов** (Stage-1 ЕРН по умолчанию; Stage-2/3
   opt-in) + chain-dominance + priority-express порядок. —
   `scripts/gamma/cli.py:356-365`; changelog `README.md:65-67` (v1.15.1/v1.15.2).
5. **Активность и MDA по ISO 11929**, опц. TCS-коррекция, самопоглощение,
   matrix-method. — `scripts/gamma/cli.py:120-156, 424-463`;
   changelog `README.md:63, 73` (v1.13.0/v1.16.x).
6. **Комплект отчётов**: JSON + Markdown + HTML + PNG-графики + Technical PDF +
   переносимый XML-экспорт. — `scripts/gamma/cli.py:182-230`;
   `docs/RUN_SKILL.md:181-208` (bundle layout).
7. **End-to-end оркестратор без AI-babysitting** `run_skill.py`: сырой `.spe` →
   полный bundle, с `--batch` и `--resume`. — `docs/RUN_SKILL.md:9-91`.
8. **Библиотека детекторов** (11 папок: `AtomSpectra`, `GP_HPGe20`, `Gamma-1C`,
   `Handy_HPGe/LaBr/NaI`, `Simple_*`). — `ls detectors/`; changelog `README.md:27-29` (v1.26.0).
9. **RAG-поиск по библиотеке знаний** (`references/books/`, BM25):
   `query/explain/cite/verify/rebuild/stats`. — `scripts/gamma/cli.py:472-512`.
10. **Реестр канонических имён/синонимов** для детекторов/геометрий
    (`data/aliases.json`). — `CLAUDE.md:25-32`.

> ⚠ Anti-hallucination нюанс: базовый Scope (`README.md:87-91`) гласит «efficiency
> не вычисляется, результат в cps». Активность в Bq/kg появилась позже как
> **opt-in** (`--sample-mass-kg`, `scripts/gamma/cli.py:293-299`). В лендинге
> подаём это как «активность — опционально, при наличии калибровки эффективности».
> Анализ-конвейер **scope-locked** на LSRM NaI `Gamma-1C` (`README.md:75`).

### 4. Структура репозитория (фактическая, по `ls`, НЕ по `README.md:93-143`)

> Внимание: секция `## Layout` в `README.md:93-143` описывает *вложенное*
> skill-дерево `gamma-spectrum-analysis/`, а **не** фактический корень репозитория.
> Лендинг описывает реальный корень (проверено `ls`).

| Путь | Назначение | Источник |
|---|---|---|
| `scripts/gamma/` | Python-пакет: `activity, calibration, data, detectors, diagnostics, experimental, identification, io, knowledge, math, peaks, physics, reporting` + `cli.py`, `spectrum.py` | `ls scripts/gamma/` |
| `scripts/*.py` | Утилиты (`run_skill.py`, `convert_spectrum.py`, `version_check.py`, …) | `ls scripts/*.py` |
| `references/` | Методология, таблицы линий, спецификации алгоритмов, книги | `ls references/` |
| `detectors/` | Ассеты по детекторам (11 папок) | `ls detectors/` |
| `data/` | `aliases.json`, `nuclides.json`, `xrf_lines.json`, `anchor_patterns.json`, … | `ls data/` |
| `docs/` | `RUN_SKILL.md`, `METHODOLOGY_LESSONS.md`, dev-skill mirror, `archive/` | `ls docs/` |
| `tests/` | Регрессионные тесты | `ls tests/` |
| `evals/` | Фикстуры + сценарии приёмки | `README.md:134-142` |
| `1_Version/` | Релизные архивы (zip), включая `books_library/` | `CLAUDE.md:36-40` |
| `audit/` | Аудиты + RAG-индекс `audit/_rag/RAG_INDEX.json` | `CLAUDE.md:71` |

**НЕ в релизе / local-only** (важно для секции Structure):
- `books_library/` — рабочая копия, **полностью исключается** из релизных
  архивов (`CLAUDE.md:36-40`, F-150/F-293).
- `detectors/Gamma-1C/raw_lsrm/` — сырые `.spe`, никогда не в git (`.gitignore:32-38`).
- `archive/`, снапшоты `gamma-spectrum-analysis_*` — local-only (`.gitignore:40-47`).

### 5. Установка

- **Зависимости**: `scripts/requirements.txt` — `numpy>=1.24,<2.4`,
  `scipy>=1.10,<1.18`, `matplotlib>=3.7`, `reportlab>=4.0`, `pytest>=7.4`,
  `pytest-xdist>=3.5`, `defusedxml>=0.7.1` (`scripts/requirements.txt:14-37`).
- **Запуск пакета**: `PYTHONPATH=scripts` (`README.md:149`, `docs/RUN_SKILL.md:33`).
- **Python-версия**: локально `python3 --version` → **3.9.6**; `docs/RUN_SKILL.md:37`
  упоминает **Python 3.14** (Windows-окружение оператора). → **Open Q2** ниже.
- **`pip install -e .`**: `scripts/requirements.txt:8` упоминает dev-install,
  но `pyproject.toml` содержит **только** `[tool.pytest.ini_options]` — нет
  `[project]` / `[build-system]`, поэтому editable-install, скорее всего, **не
  работает**. Рекомендуемый путь — `pip install -r scripts/requirements.txt` +
  `PYTHONPATH=scripts`. → **Open Q3**.

### 6. Использование (минимальный пример + что на выходе)

- **Быстрый разбор файла** (Phase-0): `PYTHONPATH=scripts python -m gamma.cli
  analyze "<file>"` → JSON-сводка (метаданные, калибровка, dead-time, фон,
  потолок энергии). — `README.md:147-152`, `scripts/gamma/cli.py:42-72`.
- **Полный отчёт**: `... analyze "<file>" --full-report --output-dir ./out` →
  bundle JSON + Markdown + HTML + PNG (+ opt PDF/XML). — `scripts/gamma/cli.py:11-15, 182-230`.
- **Оркестратор «всё-в-одном»**: `python scripts/run_skill.py sample.spe` →
  bundle в `$GAMMA_DEMO_REPORTS_DIR/<stem>/`. — `docs/RUN_SKILL.md:50-59`.

### 7. Конфигурация

- **Детекторы**: 11 папок в `detectors/` (`ls detectors/`); канонические имена и
  синонимы — `data/aliases.json` (`CLAUDE.md:30`).
- **Геометрии**: Marinelli / Дента / Петри / Точечная — через alias-реестр
  (`CLAUDE.md:30`, упоминания в changelog `README.md:77`).
- **Ключевые флаги CLI**: `--background-auto {off,suggest,apply}` (default
  `apply`, `cli.py:330-342`), `--sample-mass-kg` (`cli.py:293-299`),
  `--sample-density-g-cm3` (`cli.py:300-309`), `--allow-stage2/3`
  (`cli.py:356-365`), `--peak-search-method` (`cli.py:310-317`),
  `--enable-tcs-correction` (`cli.py:425-430`).
- **`run_skill.py` config** (`config.yaml`/`.json`) + env `GAMMA_DEMO_REPORTS_DIR`
  (`docs/RUN_SKILL.md:115-159, 42`).

### 8. Требования / ограничения (known limitations)

- Энергодиапазон **0–3000 keV**; каналы выше потолка отбрасываются
  (`README.md:89`).
- В базовом scope efficiency ε(E) / масса / объём / расстояние **не вычисляются**;
  результат в **cps** (`README.md:90-91`); активность — opt-in (см. §3 нюанс).
- Анализ-конвейер **scope-locked** на LSRM NaI `Gamma-1C`; AtomSpectra/RadiaCode
  распознаются только на уровне alias (`README.md:75`).
- Bootstrap-калибровка имеет известные проблемы на NaI 50×50 фонах с Pb-XRF
  (`README.md:169`).
- Полный журнал — `KNOWN_AND_FIXED_ISSUES.md` (`README.md:179-188`).
- Жизненный цикл: **Phase 1 — Development & Validation** (не GA),
  `docs/spectravibe-dev_skill_mirror_2026-06-06.md:37-39`.

### 9. Дизайн механизма «автоподдержки» лендинга (требование B)

- **Канонический дом правил (в этом репо)**:
  - `.claude/skills/readme-sync/SKILL.md` — project-skill (сейчас в `.claude/`
    только `agents/`, папки `skills/` нет — создаём; `ls .claude/`).
  - `docs/readme-maintenance/README_CONTRACT.md` + `UPDATE_CHECKLIST.md`.
  - User-level `spectravibe-dev` — **отдельный репозиторий**, сюда не коммитим
    (`CLAUDE.md:4`, mirror header `:3-5`).
- **Встраивание в `agent-c-docs`** (владелец `*.md`, бампит `SKILL_VERSION`,
  `.claude/agents/agent-c-docs.md:30-38`): добавить `readme_en.md` + `readme_ru.md`
  в зону «Can edit» и триггер «на merge цикла / новой фиче обновить оба лендинга
  по `docs/readme-maintenance/UPDATE_CHECKLIST.md`»; явно отметить, что
  **`README.md` — pipeline-owned** (правится штатным релизным процессом,
  не landing-процессом).
- **Указатель в `CLAUDE.md`**: 1–3 строки на skill+contract. Mirror dev-skill
  (`docs/spectravibe-dev_skill_mirror_*.md`) помечен «do not edit»
  (mirror header `:3`) — **не трогаем его**.

---

## Часть B — Карта секций лендинга (EN и RU, по образцу latte)

Образец `latte` (`raw.githubusercontent.com/arinltte/latte/main/README.md`):
центрированная HTML-шапка (лого + слоган + бейджи + переключатель языка
`English | 中文文档`), далее emoji-секции `🏗️ Features`, `Requirements`,
`🚀 Installation`, `Getting Started`, `📂 Data`, `🌐 Supported Sites`,
`📜 License`.

Структура **обоих** лендингов (`readme_en.md` / `readme_ru.md`), одинаковая по
секциям, эквивалентная по содержанию:

1. **Header (centered)** — название `SpectraVibe` + лого-плейсхолдер (TODO: реального
   ассета нет) + слоган.
2. **Бейджи** — версия `v1.27.0` (`README.md:3` / `json_report.py:30`), Python,
   License (**TODO/Open Q7** — `LICENSE` в корне не найден).
3. **Переключатель языка** — `English | Русский`, ссылки на соседний файл
   (`readme_en.md` ↔ `readme_ru.md`).
4. **Описание** (§1).
5. **✨ Особенности / Features** (§3, 8–10 пунктов).
6. **Требования / Requirements** (§5, §8).
7. **🚀 Установка / Installation** (§5).
8. **▶️ Быстрый старт / Usage** (§6).
9. **⚙️ Конфигурация / Configuration** (§7).
10. **📂 Структура / Structure** (§4, с пометкой «не в релизе»).
11. **🧩 Форматы / Formats** (§3.1, из `FORMAT_REGISTRY.md`).
12. **📚 Документация / Docs & links** — `SKILL.md`, `ARCH.md`,
    `FORMAT_REGISTRY.md`, `INDEX.md`, `docs/RUN_SKILL.md`.
13. **🆕 What's new** — короткий блок → ссылка на `README.md` как полный
    changelog / история версий.
14. **📜 Лицензия / License** (**TODO/Open Q7**).

---

## Часть C — План файлов

**Создаём:**
- `readme_en.md` — EN-лендинг (стартовая страница).
- `readme_ru.md` — RU-лендинг (стартовая страница).
- `docs/readme-maintenance/README_CONTRACT.md` — контракт лендинга.
- `docs/readme-maintenance/UPDATE_CHECKLIST.md` — чеклист обновления.
- `.claude/skills/readme-sync/SKILL.md` — project-skill автоподдержки.

**Правим (только `.md`):**
- `.claude/agents/agent-c-docs.md` — зона ответственности + триггер.
- `CLAUDE.md` — короткий указатель (1–3 строки).

**НЕ трогаем:** `README.md`, любые `scripts/**` и `*.py` (вкл.
`version_check.py`), git-теги, dev-skill mirror.

**Коммиты (Phase 3, conventional, без push):**
- `docs(readme): bilingual latte-style landing pages (readme_en + readme_ru)`
- `docs(readme-sync): maintenance skill + rules + checklist for landing pages`

---

## Часть D — Открытые вопросы ко мне (оператору)

1. **[RESOLVED 2026-06-06]** version_check / SKILL_VERSION / tag mismatch при
   написании-обновлении `readme_*` **допустим и игнорируется**: репозиторий
   обновляется постоянно, поэтому расхождение версий (`git tag='v1.27.6' !=
   SKILL_VERSION='v1.27.0'`) — ожидаемое состояние, не блокер для landing-работы.
   Landing-файлы не меняют `README.md`/`json_report.py`/теги, поэтому
   `version_check --allow-no-tag` останется красным по этой pre-existing причине —
   это нормально. Правило вносится в `docs/readme-maintenance/` (контракт +
   чеклист) и в `.claude/skills/readme-sync/SKILL.md` на Phase 3.
2. **Python-версия в Install**: писать `3.14` (по `docs/RUN_SKILL.md:37`) или
   нейтрально `Python 3.9+` с пометкой? (локально 3.9.6).
3. **`pip install -e .`**: подтверждаете, что документируем путь
   `pip install -r scripts/requirements.txt` + `PYTHONPATH=scripts` (без editable,
   т.к. в `pyproject.toml` нет `[project]`)?
4. **License**: `LICENSE` в корне репо **не найден** (`ls` корня). Что писать в
   секции License — `TODO: выбрать лицензию`, либо указать конкретную (какую)?
5. **Лого/бейджи**: реального лого-ассета нет — ставлю текстовый/emoji
   плейсхолдер + `TODO: logo`. Бейджи — version + Python + (license после Q4).
   Подойдёт?
6. **Структура в лендинге**: подтверждаете, что описываем **фактический** корень
   (scripts/gamma, detectors/, …), а не устаревшую `## Layout` из `README.md`?

**СТОП. Жду «ОК» перед Phase 2 (генерация `readme_en.md` + `readme_ru.md`).**
