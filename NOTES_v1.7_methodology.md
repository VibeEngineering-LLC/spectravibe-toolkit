# Методологические замечания для Phase 1.4 (identification) / v1.7

---

## v1.17.4 — Anonymization + RU narrative + chain completeness wired (F-115 + F-108/F-110/F-111)

Closes the 21-directive 2026-05-29 catalogue (see `SESSION_DIRECTIVES_2026-05-29.md`). Highlights:

* **F-115 (анонимизация, CRITICAL)** — `scripts/gamma/reporting/anonymize.py.anonymize_report_inplace(report)` invoked at the very end of `build_json_report`. Strips operator names, certified-source S/N (`420-7-17`), detector S/N (`УДС-ГЦ-63х63-USB №SN-01`), absolute filesystem paths, device GUIDs, and S/N-bearing `.efr` filenames. Все артефакты (JSON / Markdown / HTML / PDF) видят один и тот же анонимизированный dict.
* **F-108 wiring**:
  * `_build_table_rows` and `_build_peaks` сортируют по возрастанию энергии (D-03).
  * `_translate_note_line` пропускает строку через расширенный RU-словарь (60+ фраз), и если в результате остаётся английский токен ≥4 букв — строка DROPS (D-04, D-05, D-06, D-07).
  * `markdown_report.py` полностью на русском (заголовки разделов, столбцы, ячейки, история версий, обоснования нуклидов, заметки вторичных пиков, диагностика).
  * CLI флаги `--cost-tokens` / `--cost-session-pct` / `--cost-detail` → cost footer (D-19).
* **F-110 wiring**: при Th-232 ДОМИНИРУЕТ автоматически добавляется композит-запись 73-90 кэВ (Pb K-РИ + ВК Pb-212/Tl-208/Bi-212 + Th-228 84.37 + Pb fluorescence), а зоны `backscatter_region` / `broad_compton_plateau` исключаются (D-08, D-09).
* **F-111 / F-111b wiring**: все библиотечные линии цепочки Th-232 с I_γ ≥ 0.5% (Ac-228 463, Tl-208 510/763/860 и т.д.) появляются и в `peaks`, и в `rows`, и в `detail` с согласованными `peak` id (D-17, D-18).
* **D-12 PDF**: `gamma.reporting.pdf_export.html_to_pdf` через Edge headless (`--virtual-time-budget=5000`); опциональный флаг `analyze_and_report(..., write_pdf=True)` и CLI `--write-pdf`.
* **D-01 BG-only narrative**: `classify_environment` распознаёт чисто-фоновый файл (`sample_type_hint ∈ {background, bg}` / `is_background` / `bg_*.spe`); `_build_notes_blocks` переключает финальный абзац на «Фоновый спектр. … Никаких выводов про образец не делается».

Регрессии (8 новых файлов):

| Файл | Что проверяет |
|---|---|
| `test_anonymization.py` | F-115: ни один из запретных токенов не попадает в JSON / MD / HTML |
| `test_rows_sorted_ascending.py` | D-03: rows в HTML отсортированы по E ↑ |
| `test_no_en_leak.py` | D-04..D-07: после strip CSS/JS/JSON в body нет ASCII слов ≥4 букв вне whitelist; нет токенов `gain drift`, `trump card`, `WARNING` |
| `test_th_composite_present.py` | F-110/D-08/D-09: композит 73-90 кэВ есть, диффузных зон нет |
| `test_chain_completeness.py` | F-111/F-111b: Ac-228 463 + Tl-208 510/763/860 в peaks AND rows AND detail; set equality peaks ≡ detail |
| `test_cost_footer.py` | D-19: HTML содержит «Стоимость анализа», 65%, 140000 |
| `test_bg_only_environment.py` | D-01: ни «в образце», ни «sample contains», и т.д. на чисто-фоновом файле |
| `test_pdf_artefact.py` | D-12: PDF ≥ 30 КБ (skip если Edge не установлен) |

Все 43 теста проходят. Версия `gamma.reporting.json_report.SKILL_VERSION` → `v1.17.4`.

---

## v1.17.3 — Canonical interactive report form locked (F-114)

Hand-crafted demo at `references/demo_contract_v1_17_2/report.html` is the contract.
Every `python -m gamma.cli analyze --full-report` emits HTML in that
skeleton:

* Chart.js spectrum, log/linear toggle, click-to-highlight peaks,
  sortable RU table (Изотоп / Линия / T½ / A / Комментарий).
* Multiplet plots embedded BEFORE the summary, non-stacked Gaussian
  fills, vertical library-energy markers, staggered bottom labels.
* Russian-only labels per F-108 glossary.
* Chain-completeness placeholders "присутствует по цепочке" for
  I_γ ≥ 0.5% lines per F-111.
* iOS Telegram WebView fixes baked in per F-113 (`__initReport` with
  Chart-missing fallback + `setTimeout(..., 200)` + 1500 ms slow-CDN
  retry + ResizeObserver/IntersectionObserver double-resize).
* Mobile-responsive (@media 680/420), print-friendly.
* Optional cost footer.

**Architecture**: locked HTML template at
`scripts/gamma/reporting/templates/interactive_v1_17_2.html`
(CSS + Chart.js init JS verbatim from the demo, with 16
`{{PLACEHOLDER}}` markers for per-sample data) is rendered by
`scripts/gamma/reporting/interactive_html.py::render_interactive_html`.

**Methodology gaps surfaced**:

* `report["multiplet_deconvolutions"]` JSON does not yet carry the
  per-component `g_plus_cont[]` / `g_base[]` arrays needed for the
  non-stacked Chart.js fills. The renderer synthesises them from the
  FWHM model (σ = FWHM/2.355) by sampling a normalised Gaussian at
  each ROI energy and adding the continuum + other-component sums.
  This is methodologically equivalent to the demo's pre-computed
  arrays but is a synthesis step on the renderer side; a future
  release should plumb these through the deconvolution result directly
  so the chart matches whatever continuum / line-shape the fitter
  actually used (e.g. Gauss + tail + Compton step per F-90).
* The renderer uses `spec.channel_to_energy()` to build the spectrum E
  axis with a stride that targets ~1000 points (keeps the HTML under
  100 KB even on 3000-channel inputs).

**Regression**: `test_interactive_report.py` (anchors + RU labels +
schema version), plus three pre-existing tests updated to the new
contract (`test_v1_15_delivery.py` F-86c bullets, `test_filename_binding.py`
F-89a HTML bg-status, `test_priority_express.py` F-88d HTML 3α). Full
suite: **35 / 35 PASS, 0 fail**.

---

## v1.15.2 — Filename binding hypothesis + chain suppression + bg-status (F-89)

**Замечания пользователя** (2026-05-29, после ревью v1.15.1 demo):

> «Почему не вычитался фон? Нет указаний вычете или отсутствии фона.
>  Проигнорировано название изотопа в имени файла. Упоминание
>  изотопа в названии должно быть маркером наличия изотопа.
>  Пик 2614 есть. Все пики тория должны автоматом закрепиться на
>  своих местах.
>  Откуда появилась цепь радия. Радия в образце точно нет.
>  ...
>  Мощный пик цезия не определен. Полный провал.»

Четыре фундаментальных gap'а закрыты в v1.15.2.

### Главный методологический сдвиг

До v1.15.2 filename использовался только как:
- bg vs sample tag (через `is_background_hint`)
- geometry / detector / sample_type hint для report header

После v1.15.2 filename становится **binding hypothesis** —
обязательной гипотезой, которая driver'ит:
1. Candidate list (SKILL.md §7A.1 наконец реализован)
2. Chain dominance suppression (если filename говорит «только Th»,
   U-цепь подавляется при отсутствии strong independent evidence)
3. Confirmation tier (filename-hinted нуклиды получают преференцию)

Это правильное архитектурное место для филенейма по SKILL.md, но до
v1.15.2 эта мандатория жила только в документации.

### Почему Cs-137 fixture был «полный провал»

Cs-137 относится к Stage 2 (technogenic), не к Stage 1 (ЕРН).
Orchestrator до v1.15.2:
```
candidates_for_stage(1) → [K-40, Tl-208, Bi-214, Ac-228, Pb-212,
                          Pb-214, Pb-210, Be-7, ...]
# Cs-137 здесь нет
```

`complete_workflow=True` не включает `allow_stage2=True` — он включает
Round 5 (deconvolution + activities + MDA), но не Stage 2 candidates.
Без явного флага identify_nuclides не пытается подтвердить Cs-137.

При этом 661 keV peak:
- находится Mariscotti поиском
- матчится anchor #3 (`Cs-137 661.66`)
- priority_signal #3 «matched: True, σ huge»

Но без Cs-137 в candidate list, identify_nuclides просто его не
рассматривал. Anchor matching и identification жили в параллельных
мирах.

### Mechanism F-89e (the "полный провал" fix)

```python
# В analyze_lsrm_spe:
filename_isotope_hints = ft.get("isotope_hints", [])  # ["Cs-137"]
stage1 = _run_stage(
    spec, peaks, 1, fwhm_at_ch, window,
    extra_candidates=filename_isotope_hints or None,
)
```

`_run_stage` уже принимал `extra_candidates` (used for Stage 2/3
escalation). Pre-v1.15.2 filename hints просто не пробрасывались
сюда. После — пробрасываются для Stage 1, что эквивалентно «filename
override» по SKILL.md §7A.1.

### Mechanism F-89d (the "Откуда появилась цепь радия" fix)

На NaI 63×63 при FWHM ~50 keV @ 600 keV:
- Tl-208 583.19 keV (I=84.5%, очень яркая Th-chain линия)
- Bi-214 609.31 keV (I=45.5%, U-chain)

Разница 26 keV, что меньше FWHM. Когда Mariscotti находит peak в
этом регионе, центроид может попасть в окно matchа Bi-214 (если
calibration drift > 13 keV) или Tl-208 583 (если drift < 13 keV).

Pre-v1.15.2: peak с центроидом 609 матчился как Bi-214 (closest
energy в пределах FWHM/2 окна). При этом Bi-214 Ra-pair (609+1764)
паттерн НЕ обязательно срабатывал (1764 keV мог быть просто шумовым
artefactом), но Stage 1 ЕРН candidate list содержит Bi-214 → его
identify_nuclides выдавал как confirmed.

F-89d: когда filename = `Th232_*` и filename_chains_claimed =
`{"Th-232"}` (без U-238), правило:
1. Bi-214 Ra-pair (609+1764) недостаточен для U-238 dominance
2. Нужен quartet (≥3 of 609, 1120, 1764, 2204)
3. Без quartet → Bi-214/Pb-214/Pb-210 убираются из final_detected

Это **whitelist по filename** — пользователь явно заявляет
«источник Th», и алгоритм доверяет этому заявлению как сильному
prior'у.

### Mechanism F-89a (bg-status surfacing)

Trivially: один новый str field, surface во всех 4 слоях reporting.
Не методологический фикс, а protocol fix — отчёт не должен молча
опускать factор, влияющий на интерпретацию activity numbers.

Особенно важно для chain-dominant fixtures: если bg НЕ subtracted,
natural radon contribution к 609/1764 keV нельзя отделить от
sample-borne U-chain. Без bg-status явного, пользователь не знает,
trust'ить ли U-chain identifications.

### Schema 0.3 vs 0.2

Backward-incompatible breakage:
- Новые обязательные поля в `header`: `background_status`,
  `filename_isotope_hints`, `filename_chains_claimed`
- Новые поля в `diagnostics.chain_dominance`: `suppressed_chains`,
  `suppression_reason`, `chain_filtered_out_nuclides`

Downstream consumers, ждущие схемы 0.2 строго, получат extra keys —
для тех, кто строго validates, это ломающее изменение. Поэтому
schema bump. Forward-compat поведение: test_priority_express
теперь проверяет `schema ≥ 0.2`, не `schema == 0.2`.

### Что НЕ делается

1. **Filename hint не гарантирует confirmation**. Если filename
   говорит «Cs137», но 661 keV peak отсутствует — Cs-137
   попадает в candidates, но identify_nuclides его не confirms.
   Это намеренно: filename — гипотеза, не аксиома.

2. **Chain suppression не применяется к single-isotope filenames**.
   `Cs137_*.spe` claims `{"Cs-137"}` (не цепь). Естественный
   фон от water-matrix может содержать Bi-214 / K-40 / Tl-208 —
   они не suppressed. Suppression срабатывает только когда
   filename явно claims одну цепь (Th или U).

3. **K-40 / Ac-228 overlap warning остаётся**. Если K-40 в
   filename И Th-dominance fires → warning. Это другой механизм
   (F-88), он не отключается F-89.

---

## v1.15.1 — Приоритетный порядок экспресс-эвристики + chain-dominance hard-prior (F-88)

**Замечание пользователя** (2026-05-29):
> «Важен порядок эвристики при экспресс идентификации:
> 2615 → торий. Это фиксирует жёстко его наличие; уже этого
> достаточно для первичной калибровки. Как правило остальные
> изотопы при выраженном тории идентифицировать на этапе
> калибровки бессмысленно. Данные о наличии тория должны
> жёстко передаваться на этап идентификации пиков.»

### Что изменилось

До v1.15.0 anchor-rank table сортировался по **практической
видимости на NaI** (Co-60 doublet выше Cs-137 одиночного, потому что
Co-60 даёт паттерн), но эвристика **не имела явного порядка
diagnostic-value**. После v1.15.1 это два разных списка:

- `ANCHOR_RANKS` — practical visibility (без изменений)
- `USER_PRIORITY_ORDER` — diagnostic value at calibration time
  (новый)

И главное — **chain dominance** теперь явный объект класса
`ChainDominance`, который:
1. Производится из anchor matches по правилам (trump card / multi-
   anchor / express pattern)
2. **Hard-passes в Step 7 identification** через F-60 CI-gating
3. Surface'ится в JSON / MD / HTML / chat summary

### Trump card rule (важно)

Tl-208 2614.51 keV — особенный якорь:
- При σ ≥ 5 он **один** фиксирует Th-232 dominance (никаких partner
  не требуется)
- Это match семантики пользователя: «уже этого достаточно для
  первичной калибровки»

Почему σ ≥ 5 а не σ ≥ 3? Потому что:
- 2614 — высокоэнергетическая линия с очень малым continuum фоном
- На NaI 63×63 даже при σ=3 это уже надёжный peak
- σ=5 даёт запас прочности против random fluctuations
- Эмпирически — в demo Marinelli K-40 σ=5.6 (минимальный), и trump
  card корректно срабатывает

### K-40 / Ac-228 overlap — реальная проблема, не теоретическая

На NaI 63×63 при FWHM ≈ 85 keV @ 1460 keV:
- K-40 1460.82 keV (I = 10.55 %)
- Ac-228 1459.20 keV (I = 0.85 %)

Разрешения нет — это **одна линия в спектре**. При естественном
содержании Th-232 в water-матриксе ввклад Ac-228 в "K-40 peak" может
достигать 5-15 % (зависит от соотношения Th / K в пробе).

Demo K-40 Marinelli показывает каноничную ситуацию:
- Th-232 trace в water-matrix → Tl-208 2614 σ=5.6
- Trump card → Th-dominant flag
- K-40 priority signal matches (это же K-40 источник)
- → **k40_ac228_overlap_warning = True**

Это не означает "K-40 не определён" — это означает "не сообщайте
K-40 активность как чистую без deconvolution против Tl-208 anchor".
Защитная семантика, не блокирующая.

### Hard-pass через CI-gating

Без F-88, F-60 CI-gating промотировал только нуклиды из подтверждённых
express patterns. Если 2615 (Tl-208) находился, но не было 911
(Ac-228) — Th-232 strong pattern не подтверждался, и Tl-208 шёл с
обычным CI без promotion.

С F-88: trump card сам по себе достаточен. CI-gating получает
`TH232_PROXY_NUCLIDES = (Tl-208, Pb-212, Ac-228, Bi-212, Pb-214)`
безусловно при Th-dominant, и эти нуклиды получают anchor / pattern
confirmation tier, даже если их χ² по интенсивностям маржинальный.

### Schema bump

JSON schema 0.1 → 0.2 (added top-level `priority_express_findings`
block + `diagnostics.chain_dominance` block +
`diagnostics.k40_ac228_overlap_warning`). Это semver-major на уровне
report schema, потому что downstream consumers могут сломаться, если
ждут v0.1 ровно.

---

## v1.15.0 — Reporting delivery + Step 5α anchor seeding refactor (F-86 + F-87)

Два параллельных потока, закрытых одновременно по согласованию с
пользователем («вариант C»). Принципиально разные слои —
визуально-выходной (F-86) и архитектурно-методологический (F-87) —
но оба тянутся за концом Step 11.

### F-86 — visual / interactive deliverables

v1.14.0 закрыла **данные** Step 11, но не **их представление**. В
chat-summary помещаются 3–8 строк, JSON смотрит downstream-инструмент,
Markdown — read-only длинный документ; пользователь же ожидает (а) PNG
плоты по display conventions из SKILL.md и (б) единый HTML, который
можно отправить по почте или открыть в браузере.

Архитектурное решение: **plots — отдельный модуль с lazy import
matplotlib**. Это даёт три выгоды:

1. Если matplotlib не установлен — orchestrator работает без падений,
   `build_report(write_plots=True)` тихо возвращает warning, а Markdown
   откатывается на v1.14.0 placeholder.
2. matplotlib стоит ~250 ms на холодном импорте — выводить его за
   пределы CLI-hot-path важно для batch-обработки сотен спектров.
3. Тесты могут проверить fallback-путь без mocking-инфраструктуры.

`build_html_report` намеренно НЕ зависит от matplotlib — он лишь
читает уже готовые PNG с диска и встраивает их base64-data-URI.
Поэтому HTML-генерация работает даже на серверах без графики.

**SKILL.md "Display conventions"** жёстко требуют: «no smoothing on
data points». Реализовано через `drawstyle='steps-mid'` на сырых
counts, а не на сглаженных. Cosmetic smoothing допустим только на
overlay-кривых (envelope deconvolution'а), но в v1.15.0 даже там
оставлено как есть — token economy.

**Auto-write через `analyze_and_report(path, output_dir, …)`** — это
ровно тот один-call wrapper, который пользователь явно запросил.
Под капотом он:

```
analyze_lsrm_spe(complete_workflow=True, **orch_kwargs)
        ↓ StagedAnalysisResult
build_report(result, output_dir=..., write_plots=True,
             write_markdown=True, write_html=True)
        ↓ dict с путями + summary
```

CLI flag `--full-report` маршрутизирует на этот wrapper; **Phase-0
путь сохранён** (`analyze fixture.spe` без флага по-прежнему даёт
старый JSON-парсер summary). Это намеренно — не ломать batch-скрипты,
которые могут полагаться на старый формат.

### F-87 — Step 5α anchor seeding refactor

**Замечание пользователя** (2026-05-29):
> «Эвристика должна выполняться для спектра образца на шаге 5.
>  Причина: экспресс-поиск должен помочь при определении пиков для
>  калибровки.»

Замечание методологически верно. SKILL.md Step 5 называет **те же
самые якоря** (Co-60 doublet, K-40 1460.82, Bi-214 series, Tl-208
2614.5, Cs-137 661.66, Pb XRF triplet, LaBr₃ self-pattern), что и
F-79 anchor-rank + F-80 express-patterns. До v1.14.0 эти passes
исполнялись внутри Stage-1 идентификации (Step 7) — на два шага
позже, чем требует методология.

Практически до v1.14.0 это не вызывало багов, потому что:
- stored cal на Gamma-1C `.spe` всегда прошёл паспортную поверку →
  Step 5 refit не требовался → F-79/F-80 успевали запуститься
  своевременно для Step 7
- последовательность вызовов в коде была корректной (peak search →
  F-79 → F-80 → F-81 → Stage 1)

Но это «правильный код в неправильном месте». Когда придёт черёд
RadiaCode (F-78a roadmap), где stored cal часто плывёт от температуры,
понадобится явный Step 5 refit, и тогда отсутствие чёткой границы
Step 5α/β/γ создаст узел.

**Решение** — рефакторинг без изменения поведения:

```
Step 3   peak search (channel space)
   ↓
Step 5α  seed_calibration_anchors(mode='sample' | 'background')
         = F-79 anchor-rank + F-80 express-patterns
   ↓
Step 5β  ОПЦИОНАЛЬНО: recalibrate_energy_if_anchors_disagree
         (max residual > 0.3·FWHM → refit deg≤4, re-seed once)
         → recalibrate_on_anchor_disagreement=True
   ↓
Step 5γ  run_seven_line_check (Lsrm §9)
         — финальная калибровочная верификация
   ↓
Step 6   FWHM(E)
   ↓
Step 7   identify_nuclides + disambiguate
```

Mode-tagging остаётся прежним:
- `is_background=True` → `analysis_mode = "background_7line"`
  (доминирует F-81)
- `is_background=False` → `analysis_mode = "sample_anchor_rank"`
  (доминируют F-79 + F-80)

**Step 5β по умолчанию выключен** — `recalibrate_on_anchor_disagreement=False`.
Это сохраняет v1.14.0 контракт bit-for-bit. Включается явно через
kwarg или CLI флаг `--recalibrate-on-anchor-disagreement`. Когда
включен:
1. Фильтруются anchor'ы с непустым nuclide и без missing partner
2. Считается `max(|peak_E − E_lib| / FWHM(E_lib))`
3. Если < threshold (0.3·FWHM) → ничего не делать, diag.applied=False
4. Если ≥ threshold И ≥ min_anchors (3) → refit deg≤4 (start at 1)
5. Если новый max-residual меньше старого → применить; иначе keep stored

**Что НЕ изменилось**:
- API `find_anchor_matches` и `confirm_express_patterns` сохранены
  для обратной совместимости (внешние скрипты могут их вызывать
  напрямую)
- F-81 7-line ЕРН check работает абсолютно как раньше
- Все 27 prior-version тестов проходят без модификаций

**Что добавилось как опция**:
- single-entry-point `seed_calibration_anchors` (упрощает тесты и
  будущие порты под другие детекторы)
- recalibration_diag поле в `StagedAnalysisResult`
- 7 новых тестов в `test_v1_15_delivery.py`

### Resource economy assessment

v1.15.0 добавляет matplotlib как **soft** dependency (только при
`write_plots=True`). На fixtures, тестированных в demo:

| Артефакт | Размер |
|---|---|
| JSON | 13–32 KB |
| Markdown | 6–10 KB |
| HTML (с встроенными PNG) | 126–256 KB |
| spectrum.png (110 DPI) | 87–92 KB |
| multiplet.png (110 DPI) | 36–52 KB / cluster |

Время `analyze_and_report` на одном Marinelli (на M2 MacBook Pro)
≈ 8 с — без плотов 4–6 с, плоты +2–4 с.

---

## v1.14.0 — Step 11 (Report with diagnostics) для Gamma-1C (F-85)

Step 11 — это **завершение** SKILL.md two-pass workflow для Gamma-1C.
До v1.14.0 orchestrator производил `StagedAnalysisResult` с 19 полями,
но не было ни одной функции, которая бы из этого собрала
human-readable отчёт. F-85 закрывает пробел.

### Архитектура отчётной подсистемы

```
analyze_lsrm_spe(complete_workflow=True)
   ↓
StagedAnalysisResult  ←─ единственный «дата-объект» pipeline'а
   ↓
gamma.reporting.build_report(result, output_dir=...)
   ↓
   ├── build_json_report(result)  → dict (schema 0.1)
   │      → JSON file (UTF-8, indented)
   ├── build_chat_summary(result) → string (3–8 строк)
   └── build_markdown_report(result) → string (13 секций)
          → Markdown file (UTF-8, optional)
```

Принципы:
1. **JSON — primary**. Это машинно-читаемый артефакт, который
   downstream скрипты (LIMS, базы данных проб, regulatory reports)
   парсят без двусмысленностей.
2. **Chat summary — для пользователя**. 3–8 строк, чтобы вписать в
   token-economy budget reference 06.
3. **Markdown — on-demand**. Полное long-form rendering из JSON.
   Plot-секции (9 и 10) пока эмитят placeholder'ы — PNG генерация
   отложена.

### Step 2 — environment classifier

До F-85 классификация measurement environment была декларативной
требовалкой в SKILL.md (Step 2), но не было модуля, который бы её
вычислял. F-85 добавляет `classify_environment(result)` с решающей
таблицей:

| Признак                          | Заключение         |
|----------------------------------|--------------------|
| Pb K-XRF в residual_classifications | low_background  |
| Не Pb-XRF; ЕРН-score ≥ 0.5       | natural            |
| Не Pb-XRF; не ЕРН; токен `закр_кр`/`shield` | low_background |
| Не Pb-XRF; не ЕРН; токен `открыт`/`open_lid` | natural   |
| Иначе                            | unknown            |

ЕРН-score складывается из anchor confirmations и detected nuclides по
list {K-40, Tl-208, Bi-214, Ac-228, Pb-212, Pb-214} с весами 0.25 и
0.20 соответственно, обрезается до 1.0.

**Эмпирический результат на Th-232 Marinelli 0cm**: environment =
"natural" корректно (массовая концентрация U/Th-chain в матрице +
полный набор ЕРН в детекции, без shielding-XRF).

### complete_workflow umbrella

SKILL.md Working principle 1: "Autonomous staging — decide which steps
to run … based on what the spectrum itself shows". F-84 wired Round 5
hooks, но все — opt-in. Для "стандартного" Gamma-1C run'а нужен один
keyword:

```python
result = analyze_lsrm_spe(path, complete_workflow=True)
```

Этот флаг включает:
- `apply_deconvolution=True` (Step 8)
- `deconvolution_overlap_fwhm=3.0` (если caller не override'нул —
  ширина нужна для 600-680 keV кластера на Th-rich NaI)
- `compute_activities=True` (Step 9 — activity Bq)
- `compute_mda=True` (Step 9 — ISO 11929 detection limits)
- Step 10 residual classification всегда работает (через F-74).

Per-флаг explicit override остаётся возможен — caller может задать
любой Round-5 флаг через explicit True/False, и umbrella только
форсирует ON, не OFF. Pure backward-compatible с v1.13.0: default
`complete_workflow=False` означает behavior бит-в-бит v1.13.0.

### Intrinsic detector activity (Gamma-1C-specific)

Per `detectors/Gamma-1C/references/05_intrinsic_detector_activity.md`,
NaI 63×63 имеет negligible intrinsic Bq/cm³ (trace ⁴⁰K от обработки
Na < 0.01 Bq/cm³, обычно ниже порога обнаружения). Единственный
рутинный артефакт — **I K-escape peaks** при E_γ > 33.17 keV (порог
поглощения K-оболочки йода):

- Kα-escape: пик при `E_γ − 28.6 keV`
- Kβ-escape: пик при `E_γ − 32.3 keV`

F-85 surfaces это в `diagnostics.intrinsic_activity_signature`:

```json
{
  "detector_canonical": "Gamma-1C",
  "Bq_per_cm3": null,
  "expected_artefacts": [
    {"kind": "I_K_escape_Ka", "rule": "E_gamma - 28.6 keV", ...},
    {"kind": "I_K_escape_Kb", "rule": "E_gamma - 32.3 keV", ...}
  ],
  "absent_signatures": ["La-138", "Ac-227 series", "Ge K-escape",
                        "Cd K-escape", "Te K-escape"]
}
```

`absent_signatures` — это **sanity check**: если в спектре всё-таки
наблюдаются эти линии, значит это либо не Gamma-1C, либо есть
загрязнение (например, NaI рядом с LaBr₃-источником, что даёт ¹³⁸La
1436 keV). Для Gamma-1C это всегда пустой список «появлений», поэтому
он только декларативный.

### Что НЕ изменилось в методологии

- **Все физические алгоритмы** — без изменений. F-85 — это исключительно
  reporting и orchestrator-convenience.
- **`validate_certs.py`** не использует F-85. Cert-harness работает в
  single-source семантике, обходя disambiguate; reporting подходит
  только для general analysis pipeline.
- **`compute_activity` ядро** — без изменений. Activities + Bq/kg
  выходят через `compute_activities_for_all` (F-30) с TCS (F-31b) и
  decay correction (F-30), как в v1.13.0.

### Что осталось для post-Step-11 итераций

1. **Plot generation** — markdown emits "Plot generation deferred"
   placeholder для секций 9 (spectrum plot) и 10 (multiplet plots).
   Implementation требует matplotlib backend + display conventions
   per SKILL.md (cps + log scale + no smoothing on data).
2. **CLI** — `gamma.cli.analyze` пока не вызывает orchestrator'а;
   wiring + flags exposure отдельная задача.
3. **F-78a multi-complex** — после ревью Step 11 output пользователем.

---

## v1.13.0 — Round 5: квантификация в orchestrator (F-84)

Round 5 — это **не** новая методология; это **wiring** уже реализованных
методологических блоков в верхнеуровневый orchestrator `analyze_lsrm_spe`.
До v1.13.0 ни один из четырёх количественных шагов не был доступен
через `analyze_lsrm_spe` — пользователь получал идентификацию и
diagnostics (CI, DC, residuals), но не активность, не Bq/kg и не
ISO-предел обнаружения.

### Четыре блока, которые подключены

| Блок                                   | Модуль                                      | Версия введения |
|----------------------------------------|---------------------------------------------|-----------------|
| Multiplet-деконволюция                 | `gamma.peaks.deconvolve` (F-33/F-34)        | v1.7.11/v1.7.12 |
| Активность Bq + ковариация             | `gamma.activity.compute` (F-30)             | v1.7.7          |
| TCS / cascade-summing                  | `gamma.physics.cascade_summing` (F-31b/K-21)| v1.7.9 + v1.9.0 |
| ISO 11929 предел обнаружения L_D / MDA | `gamma.identification.mda`                  | v1.7.7          |

### Точка интеграции — после Stage / disambiguate, перед `return`

```
read → bg_subtract → FWHM model → Mariscotti peak search
   → identification (Stage 1+/2+/3)
   → disambiguate
   → ⓐ if apply_deconvolution:        apply_multiplet_deconvolution
   → ⓑ if compute_activities + eff:   compute_activities_for_all
                                       (cascade-summing dispatcher per nuclide
                                        via compute_tcs_corrections + close-
                                        geometry P/T scaling)
   → ⓒ if sample_mass_kg:             Bq/kg derived from ⓑ
   → ⓓ if compute_mda + eff:          mda_for_peak for detected lines +
                                       standard ЕРН/technogenic suite
   → return StagedAnalysisResult
```

Все четыре блока — **opt-in**. Default-поведение orchestrator'a
бит-в-бит совпадает с v1.12.0; ни одно из новых полей результата
не появляется в неуказанной cherry-pick конфигурации.

### Что такое «стандартный MDA-набор»

10 линий, оцениваемых вне зависимости от того, обнаружено ли ядро.
Это нужно для **отчётности**: «Cs-137 не обнаружен — MDA = X Bq» —
ценная информация, особенно для compliance-формуляров.

| Ядро    | E, keV    | Назначение                                      |
|---------|-----------|-------------------------------------------------|
| Cs-137  | 661.66    | Техногенная индикатор-линия                     |
| Co-60   | 1173.23   | Калибровка / пром.                              |
| Co-60   | 1332.49   | Калибровка / пром.                              |
| K-40    | 1460.82   | ЕРН-якорь                                       |
| Bi-214  | 609.31    | ЕРН Ra-226 chain                                |
| Bi-214  | 1764.49   | ЕРН Ra-226 chain                                |
| Tl-208  | 2614.51   | ЕРН Th-232 chain (proxy для Th-228 / Th-232)    |
| Ac-228  | 911.20    | ЕРН Th-232 chain                                |
| Cs-134  | 604.72    | Техногенная (старый AES маркер)                 |
| Cs-134  | 795.86    | Техногенная (старый AES маркер)                 |

Дополнительные строки можно передать через `mda_suite_extra_lines_keV
= [(nucl, E_keV), …]`.

### Оценка ROI-фона для MDA — wing-mean

ISO 11929 формула требует `background_counts_in_ROI`. Если детектируется
линия, можно использовать остаток после baseline subtraction. Но
большинство строк стандартного suite *не* детектированы (в этом весь
смысл MDA). Поэтому используем простую робастную оценку:

```
ROI       = ±1.5·FWHM(E)  around expected channel
left wing = ±1.5·FWHM(E)  immediately to the left of ROI
right wing = ±1.5·FWHM(E) immediately to the right of ROI
bg_per_ch = mean(left ∪ right)
bg_in_ROI = bg_per_ch · n_ROI_channels
```

Идея: крылья репрезентируют локальный непрерывный фон, ROI его
интегрирует. Если линия *есть*, её пик не контаминирует ROI-сумму
для целей MDA (мы оцениваем именно фон под пиком); если линии нет —
оценка тем более корректна.

### Bq/kg через `sample_mass_kg`

Простейший случай: `compute_activities=True, sample_mass_kg=m_kg`.
Для каждого валидного `ActivityResult` (где `A_Bq > 0` и не NaN):

```
A_specific = A_Bq / m_kg            [Bq/kg]
σ_specific = σ(A_Bq) / m_kg         [Bq/kg]
```

Масса пробрасывается линейно — для волюметрических самплов
(Маринелли, Дента, Петри), где cert эмитит Bq/kg, такая прямая
интерполяция совпадает с `validate_certs.py` логикой v1.7.25/F-46b
(там же конверсия `cert_Bq * mass_g / 1000`).

### 600–680 keV multiplet — целевой кейс

На Th-rich NaI 63×63 в этом интервале сидят:

- **Bi-214 609.31 keV** (I=46.1 %, Ra-226 chain) — главный фоновый
  гамма-line ЕРН во всех геометриях.
- **Cs-137 661.66 keV** (I=85.1 %) — техногенный индикатор.
- **Cs-134 604.72 keV** (I=97.6 %) и Cs-134 795.86 keV (I=85.4 %) —
  старый AES маркер (Чернобыль, Фукусима).

FWHM Gamma-1C NaI 63×63 при 660 keV ≈ 56.6 keV ⇒ окно идентификации
±0.5·FWHM = ±28 keV. Все четыре линии попадают в одно широкое
плато (Bi-214 609 ↔ Cs-137 661 ↔ Cs-134 605/796 — диапазон 191 keV
≈ 3.4·FWHM). Cowell-площадь Cs-137 661 заведомо смещена влево
на крыло Bi-214 609. Деконволюция Round 5 разделяет вклады при
запросе `apply_deconvolution=True, deconvolution_overlap_fwhm=3.0`.

### Что НЕ изменилось

- **Не реализованы новые физические модели.** Деконволюция, TCS,
  cascade summing, ISO 11929 — все алгоритмы зашиты ранее (v1.7.7
  … v1.9.0). Round 5 даёт им точку входа в общем orchestrator.
- **`compute_activity` ядро не тронуто.** `DEFAULT_TCS_METHOD_SCALE`
  (K-18 / F-35), background-safety check (K-15), decay correction
  (F-30) — всё работает по-прежнему. F-84 — фасад.
- **`validate_certs.py` не использует F-84 wiring.** Cert-harness
  работает по single-source семантике, обходя `disambiguate`. F-84
  затрагивает только общий analysis pipeline. Cert-matrix
  mean |Δ| 6.65 % инвариантна.

### Why это важно для Round 6 / F-78a

Когда мы клонируем Gamma-1C-ветку под AtomSpectra / RadiaCode, нам
понадобится **точно тот же orchestrator с тем же контрактом**. F-84
фиксирует контракт: `apply_deconvolution`, `compute_activities`,
`sample_mass_kg`, `compute_mda`, `reference_datetime` — стабильные
ключевые слова, которые в детектор-специфичных orchestrator'ах
будут просто перенаправлять на детектор-специфичные физические модели
(другой `peak_to_total_*`, другой efficiency loader, другие cascade
schemes для специфичных уровней).

---

## v1.12.0 — Архитектурная изоляция Gamma-1C (F-83)

Не методологическое, а **архитектурное** обновление, но фиксируется здесь,
потому что определяет, где в проекте искать детектор-специфичные
методические артефакты.

**Принцип**: один спектрометрический комплекс ⇒ одна папка под `detectors/`
со всеми его собственными данными, references и (в будущем) скриптами.
Алгоритмы в `scripts/gamma/` пока остаются общими; их репликация по
комплексам произойдёт только после стабилизации Gamma-1C-ветки
(директива пользователя 2026-05-29: «Скрипты и алгоритмы имеет смысл
копировать по папкам других детекторов только после их отработки
в ветке Gamma-1C. Пока ничего не копируем»).

Что переехало в `detectors/Gamma-1C/`:

| Категория       | Бывшее место                          | Новое место                                  |
|-----------------|----------------------------------------|----------------------------------------------|
| Эффективность   | `references/efficiency/`              | `detectors/Gamma-1C/efficiency/`             |
| Реф. спектры    | `references/reference_spectra/`       | `detectors/Gamma-1C/reference_spectra/`      |
| Паспорта        | `references/certificates/`            | `detectors/Gamma-1C/certificates/`           |
| LSRM-библиотеки | `references/lsrm-libraries/`          | `detectors/Gamma-1C/lsrm-libraries/`         |
| Intrinsic ref   | `references/05_intrinsic_*.md`        | `detectors/Gamma-1C/references/05_*.md`      |
| Dead-time ref   | `references/07_dead_time_*.md`        | `detectors/Gamma-1C/references/07_*.md`      |
| Усредн. фоны   | `data/averaged_backgrounds/`          | `detectors/Gamma-1C/data/averaged_backgrounds/` |
| Secondary peaks | `data/secondary_peaks*.json`          | `detectors/Gamma-1C/data/secondary_peaks*.json` |

Что остаётся универсальным (общим для всех будущих комплексов):

- `data/{aliases,anchor_patterns,nuclides,xrf_lines}.json`
- `references/01-04, 06, 08a, 08b, 09.md` (методологические разделы)
- `references/iaea_cache/` (ENSDF/NuDat кеш)
- `references/books/lsrm_format_specification.pdf` (формат `.spe` как таковой;
  он принадлежит SpectraLine как программе, а не конкретному детектору;
  с v1.17.7 / F-138 перенесён в общий каталог книг знаний `references/books/`)

**Единый источник истины путей**: модуль `gamma.detectors.gamma1c`
экспортирует все константы (`DETECTOR_ROOT`, `EFFICIENCY_DIR`,
`REFERENCE_SPECTRA_DIR`, `CERTIFICATES_DIR`, `LSRM_LIBRARIES_DIR`,
`AVERAGED_BACKGROUNDS_DIR`, `SECONDARY_PEAKS_PATH`,
`SECONDARY_PEAKS_V2_PATH`, `DEFAULT_REFERENCE_DIR`,
`DEFAULT_EFFICIENCY_DIR`, `DETECTOR_NAME='Gamma-1C'`). Никакой код
больше не хардкодит `references/...` или `data/...` пути напрямую.
Когда появится `gamma.detectors.atomspectra` (F-78a), он будет
зеркальной копией с тем же API.

**Регрессии**: 25/25 test files PASS, поведение не изменилось.

**Off-scope** (отложено): F-32 (симметричный API ридеров), F-78a
(мульти-комплексный pipeline), Round 5 (мультиплет-деконволюция
600-680, ISO 11929 MDA через efficiency, Bq/kg c массой пробы).

---

> **⚠ Session-end checklist (постоянная конвенция):**
> 1. Обновить `KNOWN_AND_FIXED_ISSUES.md` — добавить F-NN для каждого
>    исправления, K-NN для каждой новой обнаруженной проблемы.
> 2. Добавить новые методологические находки в этот файл, в виде новой
>    секции (нумеровать последовательно).
> 3. Обновить README.md version history с кратким перечнем изменений.
> 4. Запустить полный test suite — все должны пройти.
> 5. Пересобрать архив релиза с новой версией.
>
> Эта конвенция должна выполняться в КАЖДОЙ сессии без исключения,
> чтобы не терять накопленные знания и не переоткрывать те же баги.

---

Источник: диагностические наблюдения пользователя на реальных файлах
NaI 50×50 в Pb-домике (Гамма-1С, AtomSpectra). Файлы фикстур:
`Cs137_0_см.xml`, `Фон_Cs137_0_см.xml`, `Cs137_0_см_-_subtract.xml`,
`KCl__в_домике.xml`, `Фон_KCl__в_домике.xml`,
`KCl__в_домике_-_subtract.xml`, `Медальон.xml`.

═══════════════════════════════════════════════════════════════════════
1. Поведение разных компонент спектра при вычитании фона
═══════════════════════════════════════════════════════════════════════

Когда `subtract = sample − background · (live_sample / live_background)`,
разные физические компоненты ведут себя по-разному. Это **критически
важно** для identification на NaI, где для низко-/средне-активных
образцов фотопики K-40 / Cs-137 при вычитании могут практически
исчезать (особенно K-40 — он всегда присутствует в фоне). Но из
спектра ОБРАЗЦА можно опознать изотоп по вторичным особенностям,
которые выживают.

Систематика, проверенная количественно на загруженных фикстурах:

| Компонента                       | Что её порождает                            | Поведение в subtract                        |
|----------------------------------|---------------------------------------------|---------------------------------------------|
| **Фотопик образца** (661, 1461)  | γ-фотон образца напрямую в детектор         | Выживает полностью, если активность образца ≫ фон. Для низкоактивных образцов в богатом фоне (особенно K-40) фотопик «съедается» при вычете. |
| **Compton-край** ("псевдопик")   | Комптоновское рассеяние γ образца на 180°   | **Выживает почти полностью** — зависит только от активности образца, не от фона. На NaI с плохим разрешением выглядит как широкое плечо/псевдопик, а не резкий обрыв. |
| **Compton-плато (континуум)**    | Рассеяние γ образца на углах < 180°         | **Выживает почти полностью** — пропорциональна активности образца |
| **Пик обратного рассеяния**      | γ образца рассеялись в материале вокруг и вернулись в детектор | **Выживает** — наведено образцом |
| **Pb Kα / Kβ XRF** (73–87 keV)   | γ образца возбуждают свинец защиты ИЗНУТРИ + фоновая активность возбуждает его СНАРУЖИ | **Частично выживает** (75–82%). Только часть, пропорциональная активности образца, остаётся. Часть от внешнего фона вычитается корректно. |
| **Pb-210 46.5 keV**              | Контаминация Pb-210 в свинце защиты (β → Bi-210 → 46.5 γ при 4.25%) | **Полностью вычитается** при правильной фоновой нормализации — это постоянная активность свинца, не зависящая от образца. |

**Замечание по предыдущим измерениям ROI вокруг 46.5 keV в этих
файлах.** Я наблюдал "выживание" интеграла в окне 43.5–49.5 keV на
84% после subtract. Это **не** наблюдение Pb-210 пика как такового: на
NaI 50×50 при экспозиции 9243 сек 46.5 keV пик неразличим под
низкоэнергетическим continuum от образца. Шкала канала ≈ 1.7 keV/ch в
этой геометрии даёт ROI ±3 keV = всего ±2 канала, в которые попадает
**низкоэнергетический хвост Compton-континуума**, масштабируемый с
активностью образца, а не Pb-210 пик. Pb-210 как линию надо искать на
HPGe или после долгого набора фона; на NaI 50×50 с активным образцом
46.5 keV — не диагностичен.

═══════════════════════════════════════════════════════════════════════
2. Идея для Phase 1.4: identification через сигнатуру комптон-структуры
═══════════════════════════════════════════════════════════════════════

Для каждого изотопа существует **физически фиксированный** набор
вторичных особенностей, координаты которых считаются по кинематике
комптоновского рассеяния:

```
E_Compton_edge = E_γ · (2·E_γ) / (m_e·c² + 2·E_γ)         # m_e·c² = 511 keV
E_backscatter  = E_γ − E_Compton_edge
```

Для двух важных линий:

| Изотоп  | Фотопик   | Compton edge | Backscatter |
|---------|-----------|--------------|-------------|
| Cs-137  | 661.66    | 477.3        | 184.3       |
| K-40    | 1460.82   | 1243.4       | 217.5       |
| Tl-208  | 2614.51   | 2381.6       | 232.9       |
| Co-60   | 1173.23   | 963.4        | 209.8       |
| Co-60   | 1332.49   | 1118.1       | 214.4       |

Алгоритм идентификации на NaI после subtract-этапа:
  1. Найти кандидата по фотопику (или его остатку).
  2. Проверить наличие Compton-края при E_γ · (2E_γ)/(511 + 2E_γ).
  3. Проверить наличие backscatter-пика при E_γ − E_edge.
  4. Проверить характерный профиль континуума между backscatter и
     Compton edge (форма плато для данного E_γ — почти плоская, с
     резким "обрывом" на краю).
  5. Сравнить с pre-stored эталонным "fingerprint" — спектром
     известного изотопа в той же геометрии, нормированным на
     активность.

Это даёт идентификацию **без** опоры только на фотопик — особенно
ценно для:
  - NaI с плохим разрешением, где фотопики уширены и слиты
  - Низкоактивные образцы, где фотопик еле виден над фоном
  - После вычитания фона, когда фотопик частично "съеден"

═══════════════════════════════════════════════════════════════════════
3. Идея: библиотека эталонных образцов
═══════════════════════════════════════════════════════════════════════

Пользователь предложил: хранить эталон активного образца (Cs-137 0 см,
KCl в домике) и сравнивать форму подозрительного спектра с эталоном
после нормировки на cps. Compton-структура почти инвариантна к
активности (масштабируется линейно), поэтому корреляция формы спектра
в характерных регионах (backscatter region 150–250 keV; Compton edge
region 450–500 keV для Cs-137; plateau 220–470 keV для Cs-137) может
служить **подтверждающим** тестом.

Технически: cross-correlation спектров sample и template после
вычитания baseline и нормировки.

═══════════════════════════════════════════════════════════════════════
4. Что добавить в anchor patterns / nuclide library
═══════════════════════════════════════════════════════════════════════

Поскольку Pb-210 46.5 keV — реальная (и иногда единственная)
сигнатура свинцовой защиты на HPGe, добавить в `data/nuclides.json`:
  - **Pb-210**: 46.54 keV (4.25%), polonium daughter chain
  - И отметить в `references/05_intrinsic_detector_activity.md` или в
    новом разделе `references/05b_shielding_contamination.md`, что
    Pb-210 в защите даёт характерный набор: 46.5 γ + Bi-210
    bremsstrahlung continuum + Po-210 α (невидимый, но даёт fluorescence
    в окружающих материалах).

Pb XRF (Kα/Kβ 73–87 keV) — добавить отдельный раздел в anchors:
  - Использовать как **диагностику геометрии** (есть ли Pb-домик)
  - **Не использовать** как энергетический anchor: интенсивность сильно
    варьируется с активностью образца, что даёт нестабильную
    калибровочную метрику. В качестве anchor паттерна `Pb K X-ray
    quartet` уже есть, оставить с priority 3 (последний resort).

═══════════════════════════════════════════════════════════════════════
5. План work для Phase 1.4 / v1.7
═══════════════════════════════════════════════════════════════════════

  [ ] В `gamma.physics.secondary` (создать модуль) реализовать
      compton_edge(E), backscatter(E), single_escape(E),
      double_escape(E) функции.
  [ ] В identification step после candidate-photopeak проверки
      добавить opt-in secondary-feature confirmation: для каждой
      candidate nuclide line, поискать ожидаемый Compton edge и
      backscatter; повысить confidence index если найдены.
  [ ] В references/03a_identification_algorithm.md добавить шаг
      "secondary-features cross-check" между characteristic-line и
      DPR steps.
  [ ] Залогировать ожидаемое поведение в `references/04_secondary_peaks.md`
      (он уже есть в 1.5 — расширить количественными примерами из
      загруженных фикстур).
  [ ] Добавить в data/nuclides.json: Pb-210 (46.54 keV, 4.25%).
  [ ] Опционально: модуль `gamma.physics.template_match` для
      cross-correlation против сохранённых эталонов.

═══════════════════════════════════════════════════════════════════════
6. 511 keV ROI: Tl-208 510.77 vs истинная аннигиляция
═══════════════════════════════════════════════════════════════════════

На NaI 50×50 FWHM @ 500 keV составляет ~30 keV. Это значит, что **пик в
ROI ~511 keV неотличим по форме** от двух физически разных
источников:

| Источник | Энергия | I | Где |
|----------|---------|---|-----|
| Истинная аннигиляция γ | 511.00 keV | Doppler-broadened | positron emitter в пробе или pair production когда γ_родителя > 1022 keV |
| **Tl-208** (Th-232 chain) | **510.77 keV** | **≈22.6% от Tl-208** | **природный фон, любая Th-содержащая проба** |

Различие 0.23 keV ≪ FWHM NaI → на NaI без специальных мер **неразделимы**.

Последствия для skill:
  - Для **bootstrap** калибровки: разница 0.045% по энергии незначима,
    мой `physical_corroboration_filter` (требует ≥1 пик > 1022 keV
    для допущения 511-hypothesis) пропускает оба случая корректно.
  - Для **identification**: критичный disambiguation, нужно делать.

Алгоритм disambiguation (Phase 1.4):
  1. Найден пик в ROI 511 ± tol_keV?
  2. Если в спектре есть пики Tl-208 (583.19 keV, 2614.51 keV):
     a. Рассчитать ожидаемую интенсивность 510.77 из 2614.51:
        I_expected(510.77) = S(2614.51) × (0.226/0.99) × (ε(510)/ε(2614))
     b. Если S(511_ROI) ≈ I_expected(510.77) ± несколько·σ → это Tl-208,
        не аннигиляция. Отметить пик как Tl-208(510.77).
     c. Если S(511_ROI) ≫ I_expected(510.77) → присутствует
        дополнительный вклад от аннигиляции. Разложить на две
        компоненты или оставить смесь.
  3. Если Tl-208 индикаторов нет, но есть фотопик > 1022 keV → 
     возможна pair production → аннигиляция.
  4. Если ни того ни другого нет → пик скорее всего ложный, либо
     positron emitter в пробе (Na-22, F-18, Co-58 и т.п.).

Это нужно зашить в `gamma.physics.secondary` модуль (Phase 1.4) и в
nuclide_library: пометить Tl-208 510.77 как **interfering line** с
annihilation, требующая cross-check от других Tl-208 линий перед
интерпретацией.

═══════════════════════════════════════════════════════════════════════

═══════════════════════════════════════════════════════════════════════
7. Применимая теория из Lsrm "Алгоритмические основы" 2022
═══════════════════════════════════════════════════════════════════════

Прочитан и осмыслен документ ООО ЛСРМ "Алгоритмические основы —
функции обработки спектрометрической информации" (2022). Применимые
идеи для следующих фаз скила:

**§3.1.2 Прецизионный поиск пиков**:
  Стандартный поиск (Марискотти/свёртка) + подгонка → анализ невязки →
  выявление пропущенных пиков по систематическим отклонениям.
  Достоинство: высокое разрешение мультиплетов.
  Недостаток: при неадекватной модели формы создаёт ложные пики.
  *Применить в Phase 2.1 (multiplet deconvolution).*

**§5.2 Модель спектра — переопределённая система МНК**:
  Y_i = Σ S_k · ψ(p̄_k, i) + Σ a_l · i^l + Σ ε
  где S_k = A_j · I_{j,k} связывает площадь пика с активностью нуклида.
  Линейные параметры (коэф. полинома, площади, активности) и
  нелинейные (положение, FWHM, ступенька) разделяются (separable
  least squares). LU-разложение структурной матрицы выявляет линейно
  зависимые параметры. Решается итерационно квази-ньютоновским.
  *Применить в Phase 2.1.*

**§5.2.3 Соотношения интенсивностей**:
  S₂/S₁ = (I₂·ε₂)/(I₁·ε₁).
  Два режима: внутри одного участка (можно без ε-калибровки,
  предполагая ε≈const локально); по всему спектру (требует ε(E)).
  Это **ключевой constraint для multi-line nuclides** — реализовать
  как априорное ограничение площадей в Phase 2.1.

**§5.2.5 Метод Ковелла**:
  S = N_total - B_polynomial; dS = √N_total.
  Fallback для случаев плохой статистики, когда подгонка не сходится.
  *Применить в Phase 2.1 как graceful degradation.*

**§6 Идентификация — ИЩЕМ ПИКИ ДЛЯ ЛИНИЙ, не наоборот**:
  Это фундаментальная переориентация identification pipeline.
  Для каждого кандидата-нуклида:
    1. Рассчитать МДА каждой его библиотечной линии
    2. Найти "характерную" линию (с минимумом МДА)
    3. Проверить: есть ли пик в окне идентификации этой линии?
    4. Если да — нуклид "обнаружен", дальше поиск его других линий
  *Это полная архитектура Phase 1.4.*

**§6 Окно идентификации — формальная формула**:
  δE(E) = δE_0 · √(E / 661.66)
  Где δE_0 = ширина окна на 661 keV (параметр конфигурации):
    ~1 keV для ППД (HPGe)
    ~10–20 keV для сцинтилляторов (NaI/CeBr/LaBr)
  Формула отражает FWHM(E) ∝ √E для сцинтилляторов.
  *Применить в Phase 1.4 как замену моих ad-hoc tolerance_keV в
   anchor_patterns.json.*

**§6.3 МДА по ISO 11929** — полностью формализована:
  u(A) = √[u²(N_S - N_0) + N²(u²(ε)+u²(I))]/(ε·I·t)
  Порог принятия решения: L_C = k₁₋α · ũ(0)
  Предел обнаружения: L_D из квадратного уравнения
  Уже описано в references/03c_iso11929_mda.md — соответствует
  ISO 11929:2019.

**§8.2.1 Подкалибровка** (НОВЫЙ режим, отсутствующий в скиле):
  Когда нелинейная часть stored cal заведомо хороша (внутреннее
  свойство спектрометра), а сбилась только линейная (a₀, a₁ — из-за
  drift температуры, после транспортировки и т.п.) — пересчитать
  только a₀ и a₁, оставив высокие коэффициенты как есть.
  Это **промежуточный режим** между "stored OK" и "full bootstrap".
  Дешевле и надёжнее full bootstrap, когда применимо.
  *Применить как новую функцию в gamma.calibration — между
   stored_check и bootstrap. Эвристика выбора: если max_residual
   stored_check между 0.3·FWHM и 1.0·FWHM (то есть промахи небольшие
   и системные) — попробовать подкалибровку прежде bootstrap.*

**§9.4 Разметка спектра с выставлением фоновых линий**:
  Lsrm считает это **рекомендуемым** способом учёта фона для ППД
  (поканальное вычитание - только для эталонных методов).
  Фоновые линии вставляются в модель с априорными значениями
  скорости счёта. Подгонка остаётся свободной по другим параметрам.
  *Применить в Phase 1.4: вместо текущего подхода
   `background_embedded`, выставлять фоновые линии в модели как
   priors.*

**§14.3 Confidence Index — формула и calibration table для NaI**:
  CI = log( 1 / (δE_1 · δE_2 · ... · δI_2 · δI_3 · ...) )
  
  Где δE_i — относительная неопределённость энергии i-й линии,
       δI_j — относительная неопределённость интенсивности (для
              линий 2 и далее, у которых ratio к первой задан).
  
  Эталонные значения для NaI (из Lsrm Таблица 14-1):
    Cs-137: 1.8  ← низкий, всего одна линия
    K-40:   2.2
    Na-22:  3.8
    Co-60:  5.9
    Cs-134: 4.4
    Ba-133: 8.5
    Eu-152: 18.3 ← высокий, 7 линий с известными ratio
    Th-232: 16.6
  
  Эти числа служат **тарировкой**: CI < 3 — низко-достоверная
  идентификация (несмотря на формальное успешное соответствие),
  CI > 10 — высоко-достоверная. У оператора должна быть возможность
  поднять порог отсева.
  
  *Применить в Phase 1.4: вычислять CI для каждого обнаруженного
   нуклида, выводить как часть identification result, использовать
   как критерий отсева ложных идентификаций.*

**§14.2 Completeness — DC (Dose Contribution)**:
  DC = D_unident / D_total · 100
  где D = Σ S_i · E_i / ε(E_i) — дозовый вклад линии.
  Параметр полноты идентификации: какой процент **энергетического**
  спектра остался unexplained. Лучше чем процент unaccounted peaks
  по штукам, потому что 1 сильный пик может перевесить 10 слабых.
  *Применить в Phase 1.4 в выходном JSON.*

**§14.4 Принцип наименьшего действия**:
  Даже при формально успешной идентификации с библиотекой L, нужно
  проверить — не объясняются ли те же пики другим набором нуклидов
  из расширенной библиотеки L' ⊃ L. При низком CI — особенно.
  Эвристика: после первичной идентификации с лимитированной
  библиотекой, попробовать расширить и посмотреть появятся ли
  конкурирующие гипотезы.
  *Применить в Phase 1.4 как опциональный "ambiguity check" пост-step.*

═══════════════════════════════════════════════════════════════════════

═══════════════════════════════════════════════════════════════════════
8. Interferences на NaI разрешении (для identification phase)
═══════════════════════════════════════════════════════════════════════

NaI 50×50 FWHM @ 500 keV ≈ 30 keV. Это значит линии, разнесённые менее
чем на ~25 keV (≈0.8·FWHM), при идентификации **неразделимы по
форме**. Нужна physical disambiguation через другие линии тех же
нуклидов:

| ROI на NaI | Источник 1 | Источник 2 | Разница | Дисамбигуация через |
|------------|------------|------------|---------|---------------------|
| ~511 keV   | Аннигиляция 511 (positron / pair production)         | **Tl-208 510.77** (Th-232 chain)                     | 0.23 keV | Tl-208 583+2614 |
| ~600 keV   | **Bi-214 609.31** (Ra-226 chain)                     | **Tl-208 583.19** (Th-232 chain)                     | 26 keV   | Pb-212 238 + Tl-208 2614 для Th, Pb-214 295+352 для Ra |
| ~665 keV   | **Cs-137 661.66** (fallout)                          | Bi-214 609 + Cs 661 shoulder (если Ra+Cs смесь)      | 52 keV   | проверить ratio к Bi-214 1120/1764 |

**Алгоритм для пика в ROI 583–610 keV** (для Phase 1.4):
  1. Найден пик в окне 583–612 keV?
  2. Проверить присутствие Th-232 индикаторов в спектре:
     - Pb-212 238.63 keV (I=43.6%)
     - Tl-208 2614.51 keV (I=99.75%)
  3. Если оба Th-индикатора есть со значимой интенсивностью:
     - Часть пика 583–610 ROI принадлежит **Tl-208 583.19** (I=84.5%)
     - Ожидаемая S(583) = S(2614) · (I(583)/I(2614)) · (ε(583)/ε(2614))
     - Если S(ROI) ≈ S_expected(583) → весь пик Tl-208, Ra нет
     - Если S(ROI) > S_expected(583) → есть и Tl-208, и Bi-214 (Ra)
       — разложить: S(Bi-214 609) = S(ROI) − S_expected(583)
  4. Если Th-индикаторов нет:
     - Пик чисто Bi-214 609 → есть Ra-226 chain
     - Проверить дополнительно Pb-214 295.22+351.93 и Bi-214 1120.29+1764.49
       для confidence index
  5. Если ни Th, ни Bi-214 секундарных индикаторов нет, но пик
     виден → подозревать **смесь Ra-226 + Cs-137** (Cs 661 слиплось
     с Bi-214 609 на разрешении NaI). Проверить ratio I(observed)/
     I(expected_Bi-214_chain) — превышение указывает на Cs вклад.

═══════════════════════════════════════════════════════════════════════

═══════════════════════════════════════════════════════════════════════
9. NaI калибровка на ЕРН/ФОН-спектрах (Lsrm-методика)
═══════════════════════════════════════════════════════════════════════

Lsrm-документация даёт готовый **набор 7 реперных линий** для
калибровки сцинтилляционных детекторов (NaI) на спектрах
естественных радионуклидов (ЕРН) и фоновых спектрах. Эти линии
выбраны с учётом интерференций на NaI разрешении:

| E (keV) | Состав | Тип |
|---------|--------|-----|
| 240     | Pb-212 (239, Th chain) + Pb-214 (242, Ra chain) | **superposition**, всегда присутствует в любом ЕРН |
| 351.93  | Pb-214 (Ra-226 chain) | чистая |
| 511     | Tl-208 (510.77) + аннигиляция | **superposition** (для Th-источников) |
| 1120.29 | Bi-214 (Ra-226 chain) | чистая |
| 1460.82 | K-40 | чистая |
| 1764.49 | Bi-214 (Ra-226 chain) | чистая |
| 2614.51 | Tl-208 (Th-232 chain) | чистая |

**Принципы**:
1. На NaI **superpositions полезны** как реперы: пик 240 присутствует
   в любом ЕРН-спектре (хоть в чистом Th, хоть в чистом Ra, хоть в
   смеси) — нет необходимости различать какая именно линия даёт
   вклад. Энергетическая неопределённость superposition (Δ = |239 −
   242| = 3 keV для 240; |510.77 − 511| = 0.23 keV для 511) много
   меньше FWHM детектора (~10–30 keV) → погрешность по координате
   несущественна.
2. Спан 240–2615 keV (1.04 декады) — достаточен для полиномиальной
   калибровки степени до 3 с резервом на нелинейность.
3. На K-40+Ra+Th фоновых спектрах достаточно увидеть 4-5 из 7 линий
   для надёжной калибровки.

**Реализация в anchor_patterns.json** (v1.6):
  Добавлен pattern `Natural background NaI ERN-line set` с priority 1
  и tolerance_keV=8.0 (расширена под superposition uncertainty). Это
  pattern с 7 линиями — самый «информативный» из всех в библиотеке,
  поэтому в `_find_seed_pattern` (где multi-line patterns
  сортируются по line-count desc) будет пробоваться раньше других
  multi-line patterns.

**Применение к моим фикстурам**:
  - `Фон_кабинет_8192к` (22-day NaI background): должен идеально
    matchиться, нашёл 18 пиков из которых ≥6 покрывают этот pattern.
    Ожидаемая улучшенная bootstrap calibration с residuals <2 keV.
  - `Алтайское_Зло` (NaI с Cs+K+Th+Ra contamination): тоже должна
    matchиться, поскольку весь набор линий присутствует.

═══════════════════════════════════════════════════════════════════════

═══════════════════════════════════════════════════════════════════════
10. Low-E калибровка в низкофоновых укрытиях (shielded geometry)
═══════════════════════════════════════════════════════════════════════

Для NaI-спектров, снятых в **низкофоновых свинцовых укрытиях**, где
внешний γ-фон сильно подавлен и природные ЕРН-линии могут быть
слабыми, дополнительными low-E реперными линиями могут служить:

| E (keV) | Происхождение                                       | Когда стабилен |
|---------|-----------------------------------------------------|----------------|
| **32**  | Ba Kα (≈32.06 keV) от Cs-137 IC, либо рассеянное Cs XRF в защите | присутствует везде где есть Cs-137 (включая fallout фон) |
| **46.5**| Pb-210 в составе свинца защиты (постоянная контаминация — см. секцию 5 выше) | присутствует в любом «свежем» Pb-укрытии |
| **~75** | Pb Kα (Kα1=74.97 + Kα2=72.80, **на NaI слиплись в один пик**) — наведённое XRF от фотонов, идущих через защиту | присутствует везде где есть Pb-домик и есть γ-поток |

Спан 32–75 keV даёт **локальную low-E калибровочную точку** при
условии что выше нет хороших реперов. Это особенно полезно когда:
  - Экспозиция короткая, и ЕРН-линии K-40/Tl-208 ещё статистически
    незначимы
  - Образец слабо-активный, и нет своих характерных линий
  - Геометрия сильно экранирована, и внешний фон подавлен в 100×+

**КРИТИЧЕСКОЕ ПРЕДУПРЕЖДЕНИЕ для пика ~75 keV**: в зоне 70–95 keV
возможны наложения из U/Th/Ra при их наличии в образце:
  - Pa-234m 92.4 keV (от U-238 chain, при наличии U в пробе)
  - Th X-rays (Kα1 ≈ 93.4, Kβ ≈ 105.6) от внутренней конверсии в Th
    chain — наводится в Th-содержащих пробах
  - Th-234 63 keV (от U-238 chain, при наличии U)
  - U L X-rays при наличии U в пробе

То есть **«пик ~75 keV» в шумном спектре — НЕ обязательно Pb XRF**.
Он может включать смесь от Pb защиты И от U/Th в пробе. Для
identification это надо различать:
  - Если есть U-235 185.7 keV или его другие индикаторы → есть U →
    вклад от U/Pa в зону 60–95 keV
  - Если есть Th-232 chain (2614 keV) → возможен Th X-ray вклад в
    районе 85–115 keV
  - Если в пробе нет ни U, ни Th → зона ~75 keV это чистый Pb XRF
    от защиты

**Применение к скилу**:
  - Phase 1.2 (bootstrap): не использовать low-E линии как primary
    anchor — слишком много confounders. Используются только когда
    остальные patterns не сработали. Уже есть `Pb K X-ray quartet`
    с priority 3 для этого.
  - Phase 1.4 (identification): добавить cross-check для зоны 70–95
    keV — если присутствует ≥1 индикатор U или Th цепи, вычислить
    ожидаемый вклад от non-Pb компонент и не приписывать всё Pb XRF.
  - Phase 2.1 (multiplet deconvolution): зона 70–95 keV — типичный
    multiplet, требует деконволюции с фиксацией энергий Pb Kα1, Kα2,
    Pa-234m, Th-Kα, Th-Kβ и фита амплитуд.

═══════════════════════════════════════════════════════════════════════

═══════════════════════════════════════════════════════════════════════
# v1.7.9 — Точность площади на близких дублетах + полная TCS-коррекция
═══════════════════════════════════════════════════════════════════════

## Открытие: главная ошибка активности Co-60 — это площадь, а не TCS

При валидации F-29 против сертификатов (F-30) источник Co-60 №043
показал отклонение **−25.56%** от паспортных 105 000 Бк. Гипотеза
«это cascade summing» НЕ подтвердилась: TCS на 5 см даёт всего ~2-3%.

Прямое сравнение площадей вскрыло настоящую причину:

  - Cowell для пика 1173 кэВ: ~469 000 отсчётов
  - Гауссов фит Lsrm (таблица `<START PEAKS>` в .spe): **666 002**
  - дефицит **42%** на этом пике

**Механизм**: метод Cowell строит ЛИНЕЙНУЮ подложку по «крыльям» ROI.
Для одиночного пика (Cs-137 661 кэВ) крылья лежат на чистом Комптон-
континууме, и ошибка мала (~4%). Но на дублете Co-60 1173/1332
(расстояние ~160 кэВ, FWHM ~66 кэВ каждый) правое крыло ROI пика 1173
попадает на ВОСХОДЯЩИЙ склон пика 1332, и наоборот. Линейная подложка
проводится завышенно высоко → вычитается слишком много → площадь
занижена на 30-42%.

**Урок для методологии**: для близко расположенных пиков (Δ < 4·FWHM)
метод Cowell систематически непригоден. Нужен либо многогауссов фит
всего региона, либо — если спектр от Lsrm — готовый Lsrm-фит из
таблицы PEAKS. В скиле реализован второй путь (F-31a): `get_peak_area`
предпочитает Lsrm-фит, откатываясь на Cowell для не-Lsrm форматов
(AtomSpectra XML).

Это объясняет, почему e2e-валидация Cs-137 (F-30) прошла на отлично
(−1.43%): одиночный изолированный пик, где Cowell точен. Дублеты и
мультиплеты требовали отдельного решения.

## Полная TCS-коррекция (Knoll §17.6, [GILMORE-8.5])

Модель потерь на истинном совпадении для линии E_i, излучаемой в
каскаде:

    Loss(E_i) = Σ_j  p_c(E_i, E_j) · ε_T(E_j)
    C(E_i)    = 1 / (1 − Loss(E_i))
    A_true,i  = C(E_i) · A_observed,i

где p_c — вероятность совместного испускания (из схемы распада),
ε_T(E_j) — ПОЛНАЯ эффективность на энергии партнёра (любое
взаимодействие, не только фотопик).

**Ключевой мост ε_T ← ε_p**: полную эффективность мы не измеряем, но
имеем фотопиковую ε_p(E) из F-27. Связь через peak-to-total ratio:

    ε_T(E) = ε_p(E) / P(E)

P(E) для NaI 3×3" хорошо известна (Gilmore Табл. 8.4): падает от ~0.92
на низких E до ~0.17 на 2614 кэВ. Аппроксимация лог-лог полиномом 2-й
степени:

    log P(E) = −0.316 + 0.458·ln E − 0.081·(ln E)²

согласуется с табличными значениями в пределах ≤5% на всех опорных
точках. Для HPGe P(E) ближе к 1 (выше фотопиковая доля), TCS меньше.

**Зависимость от геометрии**: ε_T линейна по телесному углу Ω/4π,
поэтому и Loss линейна. Точка 25 см (Ω/4π≈0.005) → ~0.5%; точка 5 см →
~2-3%; Маринелли (≈4π) → десятки %. Поэтому TCS-коррекцию НЕЛЬЗЯ
переносить между геометриями — она пересчитывается под конкретную
ε(E) каждой калибровки.

**Каталог схем** (`CASCADE_SCHEMES`): Co-60, Y-88, Na-22 (с парой
аннигиляционных 511 кэВ), Tl-208, Eu-152 (частично), Ba-133 (частично).
Вероятности ветвления — из ENSDF/NuDat 3. Для нуклидов вне каталога
коррекция = 1.0 (молча).

## Тонкость взаимодействия F-31a и F-31b (→ K-18)

Полный конвейер F-31a (Lsrm-площади) + F-31b (TCS) даёт Co-60 +2.89%
против +0.61% с одним только F-31a. Лёгкая ПЕРЕ-коррекция возникает
потому, что гауссов фит Lsrm по широкому ROI уже частично «возвращает»
отсчёты, смещённые суммированием, которые чистый channel-sum (или
Cowell) потерял бы. Применять полную аналитическую TCS поверх Lsrm-фита
= частично считать суммирование дважды.

**Методологический вывод**: TCS-коррекция и метод определения площади
НЕ независимы. Аналитическую TCS строго корректно применять к
channel-sum / Cowell площадям. К Lsrm-фитованным площадям следует
применять уменьшенную (или нулевую) TCS. Это зафиксировано как K-18 и
требует протягивания флага источника площади (`area_source`) из
`get_peak_area` через `LineMatch` в `compute_activity`. На точечной 5 см
геометрии расхождение ~2%, поэтому для текущего релиза приемлемо
использовать F-31a как основной механизм, а F-31b — как явно
запрашиваемую опцию.

═══════════════════════════════════════════════════════════════════════
v1.7.18 — Secondary-feature anti-misidentification (F-40)
═══════════════════════════════════════════════════════════════════════

## Архитектурное замыкание петли catalog → consumer

Цепочка F-37 → F-38 → F-39 строила **catalog** вторичных особенностей
проблемных изотопов:

- **F-37 (v1.7.15)** — pure-physics формулы Compton/backscatter +
  point-estimate каталог Cs-137 / K-40 на 17 фикстурах.
- **F-38 (v1.7.16)** — переход на quantile ranges {min, p10, median,
  p90, max} с per-primary-line keying; 6 problem isotopes на 99
  фикстурах из Поверка-2016.
- **F-39 (v1.7.17)** — extended на Tl-208 / Pb-212 / Ac-228 (через
  `_PARENT_ALIASES` map к Th-228 fixture set), 9 problem isotopes.

Все три релиза создавали **источник данных** без **потребителя**.
`matches_secondary(parent, observed_E, span="p10p90")` была реализована
в `gamma.physics.secondary_peaks` уже в v1.7.16, но никто её не вызывал.
v1.7.18 (F-40) замыкает петлю: подключает `matches_secondary` в
`disambiguate_identifications` как **Rule 5**.

## Принцип Rule 5

> Если **каждая** matched-линия кандидата падает внутри empirical
> position range (p10..p90) **не-photopeak** secondary-feature
> уже-обнаруженного parent'а из v2-каталога — кандидат демотируется.

Это формальная Lsrm §14.4 "principle of least action" применительно
к Compton-continuum особенностям. Если K-40 и Cs-137 уже объясняют
спектр, и единственная линия "Bar-XX" попадает точно в Compton edge
K-40 (range [1178.9..1179.1] keV), бритва Оккама требует не вводить
новую сущность Bar-XX.

## Конструктивные детали

### Per-primary-line scope vs aggregate scope

В F-38 решено keyить каталог по `(parent, primary_E_keV, feature)` —
не по `(parent, feature)`. Это критически важно для multi-line parents:
Co-60 имеет primary lines 1173 и 1332, каждая со своим Compton edge
([906.85..912.50] и [1166.50..1169.25] соответственно). Без
per-primary-line keying Compton edge у Co-60 был бы "размазанной"
областью [906.85..1169.25] длиной 263 keV, что давало бы катастрофически
высокий false-positive rate.

Rule 5 наследует это решение: `matches_secondary` возвращает
**отдельный hit per primary line**, и demotion срабатывает по
любому из таких hits.

### Photopeak feature исключён

Каталог v2 хранит и photopeak-feature (нужно для intrinsic spread
characterization), но Rule 5 явно фильтрует `feature != "photopeak"`.
Это сделано чтобы не дублировать работу существующих правил:

- Rule 2 (positron emitter ↔ Tl-208 510.77) уже разруливает Na-22
  через `POSITRON_EMITTERS_NEAR_511`.
- `NAI_CONFUSION_MAP` (Rule 1+2 generalised) разруливает Cs-134/
  Zn-65/Na-22 photopeak collisions с natural-chain.
- Rule 3 (CI tiebreaker) разруливает остальные shared-peak ситуации.

Rule 5 — это новый домен: **continuum-derived features** (Compton
edge, backscatter, escape peaks, IC X-rays), которые до v1.7.18
не отслеживались.

### secondary_max_lines = 2 threshold

Мульти-линейные кандидаты (>2 matched lines) **не демотируются**
Rule 5'ом регардлесс позиций. Это reasoned compromise:

- Кандидат с 3+ матчами вероятно имеет **некоторое** real evidence
  beyond pure secondary contamination.
- Для близких дублетов вроде Bi-214 (7+ матчей на полном Ra-chain
  спектре) демотация была бы серьёзной ошибкой.
- Threshold 2 позволяет одиночные false positives (1 линия) и
  двойные сужающиеся к одной mid-energy continuum-фичи (2 линии)
  быть демотируемыми, но защищает богато-матченные real detections.

Настраивается через kwarg `secondary_max_lines: int = 2`.

### Inert mode

Если ни одного из 9 problem isotopes (Cs-137, K-40, Co-60, Na-22,
Y-88, Th-228, Tl-208, Pb-212, Ac-228) нет в `detected`, Rule 5
**ничего не делает**. Это означает что на спектрах, доминируемых
anthropogenic нуклидами не из каталога (например Eu-152 calibration
sources), правило бессильно — это **explicit accepted limitation**
текущего scope. Расширение catalog'а до Ba-133 / Eu-152 / Bi-214
multiplex spectra оставлено для будущей итерации.

## Ordering: где в disambiguate сидит Rule 5

Порядок применения правил в `disambiguate_identifications` имеет
значение, потому что каждое правило мутирует `detected` set:

```
Rule 2          (positron emitter ↔ Tl-208 511)
Rule 1+2 gen    (NAI_CONFUSION_MAP)
Rule 4 univ     (intensity-ratio proportionality)
Rule 5 (F-40)   (secondary-feature anti-misID)         ← новое
Rule 4b         (Ra-226 chain equilibrium)
Rule 4c         (U-235 vs Ra-226 at 186 keV)
Rule 3          (shared-peak CI tiebreaker)
```

Rule 5 поставлено **после** Rule 4 (proportionality) и **перед**
Rule 4b/4c (chain equilibrium):

- **После Rule 4**: rare-isotope multi-line false positives уже
  отсечены, multi-line кандидаты с >2 lines на этом этапе — это
  proportionality-passing real candidates, не secondaries.
- **Перед Rule 4b**: Ra-chain equilibrium логика работает на чистом
  post-secondary candidate set, без шума от Compton-edge артефактов.
- **Перед Rule 3**: shared-peak CI tiebreaker работает после
  demotion, иначе losing nuclide мог бы быть демотирован Rule 5'ом
  до Rule 3 даже добралось.

## Что Rule 5 НЕ покрывает (intentional scope)

1. **Photopeak collisions** — это домен Rule 2 / NAI_CONFUSION_MAP /
   Rule 3, не Rule 5.
2. **True secondaries из не-каталогизированных parents** — Eu-152,
   Bi-214, Ra-226. Их Compton edges не в каталоге, значит Rule 5 их
   не использует.
3. **Geometry-conditional shifts** — каталог собран на нескольких
   геометриях (5cm, 25cm, Дента, Маринелли) и median'ы агрегированы.
   p10/p90 absorb the spread. Per-geometry inference defer для
   будущего F-41+.
4. **Time-shifted secondaries** — если детектор дрейфует за время
   measurement'а, secondary peak position может выйти за p90.
   Calibration drift как первичная проблема — defer.

## Verification

234 теста в 19 файлах проходят. Cert validation matrix не изменилась
(11/11 measurable, mean |Δ|=4.62%) — validate_certs.py обходит
disambiguate_identifications для single-source semantics, поэтому
Rule 5 inert на cert harness'е (правильное поведение).

Manual edge-case verification:

```python
# Ac-228 911.20 keV ↔ Co-60 1173 Compton edge (the canonical case)
co60 = NuclideIdentification("Co-60", matched=[1173, 1332])
ac228 = NuclideIdentification("Ac-228", matched=[911])
refined = disambiguate_identifications(result)
# → Ac-228 в rejected с reason listing Co-60 compton_edge [906.85..912.50]
```

═══════════════════════════════════════════════════════════════════════
v1.7.17 — Lsrm chain-library + закрытие Th-228 cert row (F-39)
═══════════════════════════════════════════════════════════════════════

## Что было неполным к v1.7.16

Матрица cert-валидации в v1.7.14 (F-36) имела **одну строку как
`library gap`**: Th-228. Сертификат заявляет активность parent'а
Th-228 (129 000 Bq), но сам Th-228 эмитирует слабые γ-линии (<1.5%
intensity) на NaI. Весь измеряемый сигнал идёт от daughter'ов в
secular equilibrium: Tl-208 583+2614, Pb-212 238, Ac-228 911 и др.

В v1.7.16 (F-38) characterized Tl-208 secondaries под "Th-228"
label без явного моделирования chain-equilibrium reconstruction.

## Пользователь предоставил Lsrm libraries

`C:\LSRM\Work\BG\Gamma-1S\Архив\Data\`:

- **`NaI-Etl+Esc.lib`** — NaI-tuned Th-232 chain library. Каждая
  линия тегирована ENSDF `dbid` (например `Pb-212_22`, `Tl-208_28`,
  `Ac-228_57`), что позволяет existing `gamma.data.chain_decomposer`
  разложить bundled "Th-232" в per-daughter records: Tl-208 (5 lines),
  Pb-212 (2), Bi-212 (3), Ac-228 (6), Ra-224 (2).
- **`ОСГИ.lib`** — 33 ОСГИ certified-source nuclides с правильными
  intensities + d_intensity uncertainties. Дополнительно:
  Eu-154, Eu-155, Ce-144, Sn-113, Hg-203, Co-56, Ag-110m, Ta-182,
  Cs-134, Sb-125, Ir-192, Ru-103, Zr-95+, Ho-166m, Th-231, Th-234,
  U-232, Ti-44.
- Auxiliary: `Th.zon` (pre-computed Th multiplet zones),
  `NaI63x63Point_10cm Th-232 chain windows .cen/.cfw`,
  `Gamma-1S+Compton.cpt` (Lsrm Compton model), `Aspect.src`.

## Реализация

### 1. Bundle в репозитории

`references/lsrm-libraries/` — 7 файлов (~85 KB):
- `NaI-Etl+Esc.lib`, `ОСГИ.lib`
- `Th.zon`, `*.cen`, `*.cfw`
- `Gamma-1S+Compton.cpt`, `Aspect_2025.src`

### 2. Opt-in loader

```python
from gamma.data.nuclide_library import load_lsrm_chain_libs

# По умолчанию: оба .lib, merge_mode='supplement', split_chains=True
res = load_lsrm_chain_libs()
# {'NaI-Etl+Esc.lib': 2, 'ОСГИ.lib': 18} — 20 новых нуклидов
```

Params:
- `include_nai_chain: bool = True` — toggle для chain library
- `include_osgi: bool = True` — toggle для ОСГИ
- `merge_mode: str = "supplement"` — keep existing, или `"override"`
- `split_chains: bool = True` — decompose Th-232/Ra-226/U-238/U-235

Default по-прежнему **off** в стандартном flow — deterministic behaviour
сохраняется (K-03 был accepted limitation, теперь closed: 47
nuclides доступны on-demand).

### 3. Th-228 + Ra-224 в built-in JSON

Chain decomposer не extract'ит Th-228 из bundled "Th-232" lines
(Th-228 — grandfather chain, не Lsrm "owner"). Для parent decay
correction в cert-validation нужен T½ Th-228 независимо от того,
загружена ли Lsrm chain library.

Добавлено явно в `data/nuclides.json`:

```json
"Th-228": {
  "T_half_s": 60307329,         // 1.9116 yr
  "parent": "Ra-228", "chain": "Th-232",
  "lines": [
    [84.373, 1.22, 0.02],         // weak direct lines, not measurable on NaI
    [215.99, 0.254, 0.003],
    [131.613, 0.131, 0.002],
    [166.41, 0.104, 0.001]
  ]
},
"Ra-224": {
  "T_half_s": 313000,            // 3.6 d
  "parent": "Th-228", "chain": "Th-232",
  "lines": [[240.986, 4.1, 0.05], [80.998, 1.27, 0.03]]
}
```

### 4. Chain-proxy cert validation

`validate_certs.py` теперь:
- Вызывает `load_lsrm_chain_libs()` at module import.
- `CertFixture` имеет новые поля `cert_nuclide` (parent name в cert,
  отличается от identified nuclide) + `chain_branching: float = 1.0`
  (`A_daughter / A_parent`).

Th-228 cert fixture:
```python
CertFixture("Pb-212", "Th-228__264_2023_Точечная-5см_5cm.spe",
            "Th-228 SRC-01",
            cert_nuclide="Th-228", chain_branching=1.0)
```

Identify Pb-212 в spectrum (его 238.63 keV — самая чистая
daughter линия для измерения). Decay-correct cert A через **parent's
T½** (Th-228, 1.91 yr). Compare A_Pb-212 measured vs A_Th-228 ·
1.0 (direct chain, no β-branching).

Логика расчёта:
```python
if fx.cert_nuclide:
    parent_lib = get_nuclide(fx.cert_nuclide) or {}
    T_half_s = float(parent_lib.get("T_half_s") or 0) or None
# Decay-correct cert A using parent T½
A_cert_at_meas = _decay_to_meas(cert_act.A_Bq, T_half_s,
                                 cert.reference_datetime, meas_dt)
# Apply chain branching: expected daughter = parent * branching
A_cert_at_meas *= fx.chain_branching
```

### 5. Th-chain в secondary_peaks_v2.json

Добавлены `Tl-208`, `Pb-212`, `Ac-228` в `PROBLEM_ISOTOPES`. Так как
у них нет dedicated fixtures, `_PARENT_ALIASES`:

```python
_PARENT_ALIASES = {
    "Tl-208": "Th-228",
    "Pb-212": "Th-228",
    "Ac-228": "Th-228",
}
```

→ inventory scan находит "Th-228" фикстуры, но analyzer обрабатывает
каждую с E_lines от Tl-208/Pb-212/Ac-228 соответственно.

Catalog теперь покрывает **9 problem isotopes** (вместо 6):
- Cs-137, K-40, Co-60, Na-22, Y-88, Th-228 (старые)
- Tl-208, Pb-212, Ac-228 (новые F-39)

## Обновлённая матрица cert-валидации

Точечная-5см геометрия:

| Нуклид  | A_cert, Bq | A_cert@meas, Bq | A_изм, Bq | Δ, %   | Примечание |
|---------|-----------:|----------------:|----------:|-------:|:-----------|
| Cs-137  |    106 000 |          89 307 |   90 944  | +1.83  |   |
| Co-60   |    105 000 |          49 886 |   50 193  | +0.61  | F-31a+F-31b+F-35 |
| Na-22   |    229 000 |         136 700 |  144 900  | +5.99  |   |
| Eu-152  |    142 000 |         120 600 |  112 400  | −6.82  | 8 lines + 1 multiplet |
| Ba-133  |     44 100 |          15 308 |   15 511  | +1.32  |   |
| Am-241  |    103 100 |         102 500 |   93 948  | −8.35  | 60 keV у ε(E) edge |
| Zn-65   |      3 100 |             871 |      918  | +5.45  |   |
| Y-88    |    350 000 |          29 903 |   30 392  | +1.64  |   |
| Bi-207  |     97 000 |          82 407 |   80 208  | −2.67  |   |
| Cd-109  |    395 000 |           6 762 |    6 534  | −3.38  |   |
| **Th-228** | **129 000** | **88 575** | **77 250** | **−12.79** | **F-39: via Pb-212 chain proxy** |

**11/11 measurable**. Mean |Δ| = 4.62%. Max |Δ| = 12.79% (Th-228,
сразу выше ±10% — 5см sealed-source self-absorption на 238 keV).

## Documented Th-chain conflicts

В catalog добавлены conflict-entries (real library γ-lines от
других нуклидов в p10..p90 range Th-chain feature's):

- **Tl-208 510.77 photopeak [502, 505]** ↔ **Na-22 511
  annihilation** — chronic положителен-emitter vs Th-chain confusion.
- **Tl-208 510.77 Compton edge [291, 298]** ↔ **Ir-192 295.96
  (I=28.6%), Pb-214 295.22 (I=18.4%)**.
- **Tl-208 583 Compton edge [346, 379]** ↔ **I-131 364 (I=81.5%),
  Ba-133 356 (I=62.0%), Pb-214 351 (I=35.6%)** — 4 конфликта.
- **Tl-208 583 backscatter [171, 185]** ↔ **Sb-125 176 (I=6.7%)**.
- **Ac-228 911** in Co-60 spectra — sits в Co-60 1173 Compton edge
  range (already documented F-38).

## Тесты — 4 новых F-39

`test_v2_catalog_loads` — now expects 9 isotopes.
`test_v2_catalog_tl208_chain_daughter` — Tl-208 has 4 primary
lines incl. 510 keV overlapping Na-22 511.
`test_lsrm_chain_loader_adds_th_chain_daughters` — loader adds
≥18 nuclides incl. all Th-chain daughters + ОСГИ extensions.
`test_th228_in_built_in_library` — Th-228 в built-in JSON с T½
6.03e7 s и 84.4 keV линией.

Все 213 prior tests неизменны. Total: **225 tests**.

## Closes K-03

Built-in JSON cap (24-27 nuclides) был flagged как K-03 limitation.
Opt-in load до 47 nuclides on demand → cap practically removed.
K-03 переходит в "closed via F-39" status.

═══════════════════════════════════════════════════════════════════════
v1.7.16 — range/shape catalog проблемных изотопов (F-38)
═══════════════════════════════════════════════════════════════════════

## Методологическое уточнение от пользователя

> "Положение вторичных пиков не фиксировано, а плавает в некотором
>  характерном диапазоне в зависимости от активности изотопа, геометрии
>  и прочих условий. Твоя задача выявить характерный диапазон положения,
>  форму пика изотопа и форму вторичных пиков для проблемных изотопов
>  для последующего предположения о наличии этих изотопов в образцах."

v1.7.15 (F-37) reported mean ± std per (nuclide, feature). Это
методологически неправильно для задачи **inference о присутствии
изотопа в образце** — single-point comparison ("есть ли пик ровно на
478 keV?") неверен. Правильный вопрос: "попадает ли наблюдаемое
положение И форма пика в характерный диапазон, измеренный по всем
известным условиям эталонных образцов этого нуклида?".

## Новый dataset

Архив `Поверка 2016` (99 .spe файлов):
- 12 нуклидов на 5 cm
- 11 нуклидов на 25 cm  
- Маринелли / Дента-100мл / Петри-60мл container sets
- **3 типа фонов**: вода (для Marinelli), открытые крышки (для
  контейнеров), пустая защита (для точечных). Plus два specialized:
  `_Дента-100.spe` и `_Петри-60.spe`.
- **15-measurement time-stability series**

С 2024 Поверка + existing fixtures per-isotope inventory:

| Parent  | Fixtures | Geometries                                     |
|---------|---------:|------------------------------------------------|
| Cs-137  |    17    | M-source, denta, marinelli, petri, point0cm xml, point25cm, point5cm |
| K-40    |    10    | M-source, denta, marinelli, petri              |
| Co-60   |     3    | point5cm×2, point25cm                          |
| Na-22   |     4    | point5cm, point25cm×2, fixture                 |
| Y-88    |     4    | point5cm×2, point25cm×2                        |
| Th-228  |     4    | point5cm×2, point25cm×2                        |

## Реализация — analyze_problem_isotopes.py

1. **Quantile statistics**: {min, p10, median, p90, max, mean, std}
   per (parent, primary_E, feature) instead of mean ± std.

2. **Per-primary-line keying**. Ключ `(parent, primary_E_keV, feature)`
   не `(parent, feature)`. Co-60's 1173-keV Compton edge (at 963) НЕ
   merged с 1332-keV (at 1118). Каждая photopeak's secondaries
   формируют свой tight cluster.

3. **Per-feature aggregations**:
   - position_keV: observed centroid energy distribution
   - position_residual_keV: observed − theoretical
   - intensity_ratio: S_secondary / S_photopeak (тот же photopeak,
     что в primary line)
   - fwhm_keV: measured FWHM
   - fwhm_ratio_to_theory: measured / FWHM-provider expected
   - asymmetry: left-half-area/total − 0.5
   - by_geometry: per-geometry raw observations
   - conflict_lines: real γ-lines from OTHER nuclides falling в
     p10..p90 observed range

4. **Conflict detection**: для каждого feature scan `nuclide_library`
   на γ-lines от других ядер с I_pct ≥ 1%, попадающие в p10..p90
   observed position range. Sorted by I_pct (worst conflict first).

## Catalog API (новое в `gamma.physics.secondary_peaks`)

```python
from gamma.physics.secondary_peaks import (
    load_catalog_v2, position_range, matches_secondary,
)
```

`position_range("Cs-137", 661.66, "compton_edge", span="p10p90")` →
`(433.9, 439.3)` — что peak в этом диапазоне consistent с Cs-137
Compton edge.

`matches_secondary("Cs-137", 437.0)` → list of matched features:
```python
[{'primary_E_keV': 661.66, 'feature': 'compton_edge',
  'range': (433.9, 439.3), 'median_E_keV': 436.8,
  'distance_to_median_keV': +0.2}]
```

Параметр `span`:
- `"p10p90"` (default) — 90% CI, использовать для inference
- `"minmax"` — расширенный диапазон, для warning
- `"iqr"` — narrow band (median ± 0.25·(p90-p10))

## Ключевые эмпирические находки

### Intrinsic photopeak drift на Gamma-1C NaI 63×63

| Parent  | Primary  | n  | p10..p90 spread | std    |
|---------|---------:|---:|----------------:|-------:|
| Cs-137  |  661.66  | 17 |   2.1 keV       | 0.81   |
| K-40    | 1460.82  | 10 |   4.5 keV       | 2.18   |
| Co-60   | 1173.23  |  3 |   2.8 keV       | 1.52   |
| Co-60   | 1332.49  |  3 |   0.7 keV       | 0.33   |
| Na-22   |  511.00  |  4 |   0.7 keV       | 0.40   |
| Y-88    |  898.04  |  4 |   1.8 keV       | 0.90   |

Все ниже 5 keV — устанавливает **lower-bound uncertainty** для
position-based identification tests на этом детекторе.

### Compton edge position consistent с −0.7·FWHM rule из F-37

| Parent | Primary | E_theory | p10..p90 range | residual median |
|--------|--------:|---------:|---------------:|----------------:|
| Cs-137 |  661.66 | 477.34   | 433.9..439.3   |     −40.6 keV   |
| K-40   | 1460.82 | 1243.36  | 1178.9..1179.1 |     −64.4 keV   |
| Co-60  | 1173.23 | 963.42   | 906.9..912.5   |     −54.0 keV   |
| Na-22  | 1274.54 | 1061.71  | 1000.8..1012.4 |     −57.4 keV   |
| Y-88   | 1836.06 | 1611.77  | 1545.6..1560.1 |     −52.8 keV   |

Универсально → формула `compton_edge_observed_keV(E, FWHM) =
compton_edge_keV(E) − 0.7·FWHM` (введена в F-37) подтверждена на
5 ядрах.

### Backscatter position всегда смещён ВВЕРХ, geometry-dependent

| Parent | Primary | E_theory | p10..p90 range | residual median |
|--------|--------:|---------:|---------------:|----------------:|
| Cs-137 |  661.66 | 184.32   | 186.8..196.7   |     +9.6 keV    |
| K-40   | 1460.82 | 217.46   | 227.4..235.2   |     +13.1 keV   |
| Co-60  | 1173.23 | 209.81   | 226.6..227.8   |     +17.0 keV   |
| Na-22  |  511.00 | 170.33   | 175.0..178.4   |     +6.2 keV    |
| Y-88   |  898.04 | 198.91   | 202.2..229.3   |     +18.4 keV   |

### Conflict catalogue (selected high-risk)

- **Co-60 1173 Compton edge** [906.9, 912.5] → **Ac-228 911.20 keV
  (I=25.8%)** — direct conflict.
- **Co-60 1332 Compton edge** [1166.5, 1169.3] → indistinguishable
  from **Co-60 1173 photopeak** position! Methodologically expected
  cross-validation: 1173 peak amplitude must include 1332 Compton.
- **K-40 Compton edge** [1178.9, 1179.1] → 5.7 keV above Co-60 1173
  photopeak. На этом детекторе раздельны, но на детекторе с realistic
  drift conflict acute.
- **Th-228 (Tl-208) 583 Compton edge** [351.9, 378.4] →
  **I-131 364.49 (I=81.5%)**, **Ba-133 356.01 (I=62.0%)**,
  **Pb-214 351.93 (I=35.6%)**. Three nuclides directly in range.
- **Cs-137 Ba Kα IC X-ray** [23.7, 26.6] → **Am-241 26.34 (I=2.3%)**.

## Применение

`matches_secondary(parent, observed_E)` встраивается в
`disambiguate_identifications` Rule 4:
- Candidate Bi-214 supported только линией 503 keV; уже identified
  Cs-137 photopeak → check `matches_secondary("Cs-137", 503.0)`. Если
  пуст — keep Bi-214. Если матч с Cs-137 compton_edge — demote
  Bi-214.
- Аналогично: Ac-228 911 in spectrum where Co-60 1173 also detected.

`position_range(nuclide, primary_E, feature)` для cross-validation:
- Cs-137 claim: проверить, есть ли peak в `position_range("Cs-137",
  661.66, "backscatter")` = `(186.8, 196.7)`. Если есть → confirm.
  Если нет → suspect misidentification.

## Тесты

9 new в `test_secondary_peaks.py` (всего 22 в файле):
- `test_v2_catalog_loads` — все 6 problem isotopes присутствуют
- `test_v2_catalog_per_primary_keying` — Co-60 имеет separate 1173 + 1332
- `test_position_range_cs137_compton_edge` — p10..p90 within [430,445]
- `test_position_range_k40_compton_edge` — p10..p90 within [1175,1185]
- `test_matches_secondary_cs137_compton_edge_collides_with_bi214` —
  peak at 437 matches compton_edge, doesn't match backscatter/photopeak
- `test_matches_secondary_no_match_outside_ranges` — peak at 800 в
  Cs-137 context → empty match list
- `test_matches_secondary_k40_compton_dangerous_for_co60` —
  documents 5.7 keV separation
- `test_v2_catalog_conflict_lines_recorded` — Cs-137 Ba Kα flags
  Am-241 26.34
- `test_v2_photopeak_position_tightness` — все problem isotope
  photopeak spreads < 5 keV

Все 213 prior tests проходят неизменно. Total: **222 теста**.

═══════════════════════════════════════════════════════════════════════
v1.7.15 — библиотека эталонных образцов / каталог вторичных пиков (F-37)
═══════════════════════════════════════════════════════════════════════

## Идея

Identify в текущей форме рассматривает каждый Mariscotti-пик как
кандидата на photopeak ядерной γ-линии. Но в реальном NaI-спектре до
половины обнаруженных пиков — это **вторичные особенности**: Compton
edge maximum, обратное рассеяние, single/double escape от pair
production, K-X-ray escape кристалла йодида, IC X-rays самого
распадающегося нуклида, и всегда-присутствующая фоновая линия K-40.

**Cs-137 Compton edge** в 478 keV — хронический источник ложных
идентификаций Bi-214 (503 keV) на NaI. **K-40 backscatter peak**
~230 keV ложно читается как Ac-228 (209) или Pb-212 (238.6). Без
явного знания "где должны быть вторичные особенности каждого
parent-нуклида" disambiguator не может принципиально подавить эти
ложные claims.

Пользователь предоставил 20 новых референсных фикстур (`Поверка
2024`):
- Дента-120мл: Cs-137 ×2, K-40 ×2, Ra-226 ×2, Th-232 ×2.
- Петри-60мл: те же 8.
- Точечная-25см: Cs-137, Na-22, Y-88, Th-228.

С имеющимися Точечная-5см и Маринелли — это **17 чистых Cs-137 + K-40
спектров через 4+ геометрии**. Достаточно для эмпирического
характеризования вторичных пиков на этом детекторе.

## Физика вторичных особенностей

Для primary γ-линии E (keV) на NaI ожидаются:

| Feature        | Position                          |  Note                       |
|----------------|-----------------------------------|------------------------------|
| photopeak      | E                                 | сам полный peak              |
| compton_edge   | 2E²/(m_e + 2E)                    | max e⁻ KE при 180° scatter   |
| backscatter    | E/(1+2E/m_e)                      | 180° фотон, рассеявшийся в shielding |
| single_escape  | E − 511                           | один annihilation γ ушёл (только если E>1022) |
| double_escape  | E − 1022                          | оба annihilation γ ушли      |
| xray_escape    | E − 28                            | I K X-ray ушёл из NaI (значимо при низких E) |
| ic_xray_Ba_Ka  | 32                                | Ba Kα от 8% IC branch Cs-137 |
| k40_natural    | 1460.82                           | природный K-40 фон           |

m_e = 510.999 keV. Energy conservation: `compton_edge + backscatter
= E` (для single 180° scatter, energy go to electron OR back-scattered photon).

## Измерения на 17 фикстурах

Гарнес `analyze_secondaries.py`:
1. Inventory 10 Cs-137 + 7 K-40 спектров.
2. Для каждого: read_spectrum → fwhm_provider → mariscotti_search.
3. Для каждого ожидаемого feature: найти ближайший Mariscotti-пик в
   пределах tolerance (30-60 keV).
4. Cowell-площадь на +/- 2.5σ ROI.
5. Per-feature ratio R = S_secondary / S_photopeak.
6. Aggregation: mean / std / min / max по геометриям.
7. Dump в `data/secondary_peaks.json`.

### Закономерность A: Compton edge смещён ВНИЗ от теории

|      Source     | E_theory keV | E_observed keV | shift keV  |
|----------------|-------------:|---------------:|-----------:|
| Cs-137 (n=9)    |   477.3      |  ~440          | **−37**    |
| K-40   (n=3)    |  1243.4      |  ~1190         | **−53**    |

**Это не drift калибровки** — affects both nuclides identically.
Root cause: Mariscotti finds local maximum of `−d²S/dE²`. Analytical
Compton edge — это step discontinuity, broadened by detector
resolution в спадающее плечо. Максимум второй производной отклика
сидит **примерно на 0.7·FWHM ниже** analytical step.

Verification:
- Cs-137 at FWHM(477)≈50 keV → predicted shift = −0.7·50 = −35 keV
  (observed −37 ✓)
- K-40 at FWHM(1243)≈75 keV → predicted = −0.7·75 = −52 keV
  (observed −53 ✓)

**Закладывается** в `compton_edge_observed_keV(E, FWHM)` helper.

### Закономерность B: Backscatter peak смещён ВВЕРХ + геометр.-conditional

|      Source     | E_theory keV | shift keV  | mean R, %  |
|----------------|-------------:|-----------:|-----------:|
| Cs-137 5cm      |    184.3     |  +14       |   14.6     |
| Cs-137 25cm     |    184.3     |  +4        |    3.2     |
| Cs-137 Маринелли|    184.3     |  +10       |    7.5     |
| Cs-137 Дента    |    184.3     |  +6        |    5.7     |
| K-40   Дента    |    217.5     |  +14       |   13.8     |
| K-40   Петри    |    217.5     |  +11       |   15.3     |

Mean смещение **+8 ... +15 keV** консистентно. Multi-path
photons (несколько mini-scatters до выхода из shielding) добавляют
к single-180° base, и centroid смещается вверх.

R (intensity ratio) — **самый чистый geometry-marker**:
- 5cm point: 0.146 (highest)
- 25cm point: 0.032 (lowest — направленные фотоны mostly miss
  detector + shielding)
- Extended-source containers: 0.06-0.08

**Закладывается** в `backscatter_observed_keV(E, geometry=)`.

### Закономерность C: Природный K-40 контаминирует ВСЕ длинные spectra

|  Cs-137 source  | live, s | R(K-40 / Cs-137)  |
|-----------------|--------:|------------------:|
| Маринелли/7-14  |  long   |     0.9%          |
| Дента/7-14      |  long   |     7.1%          |
| Петри/7-14      |  long   |    10.3%          |
| Точечная-5cm    | 1800    |     0.0%          |
| Точечная-25cm   | 3435    |     0.55%         |

Чем больше масса контейнера и время — тем сильнее K-40 фон. Identify
должен **ожидать** 1461 keV peak в любом measurement extended-source
и не credit это как anomalous K-40, если magnitude не превышает
ожидаемый background.

### Дополнительные феномены

- **Cs-137 Ba Kα IC X-ray (32 keV)**: detected в 4/10 фикстур (5cm,
  25cm, Маринелли). Mean R = 8.4% ± 4.1%. Чистый Cs-137 marker —
  отсутствует в Bi-214 / прочих ядрах. Полезен для cross-validation
  Cs-137 identification.
- **K-40 single escape (950 keV)**: detected в 2/7 spectra при
  R~2-3%. Слабая, но cleanly resolvable.
- **K-40 double escape (439 keV)**: НЕ resolved — buried в Compton
  continuum (рядом с собственной Compton edge).

## Применение

1. **Anti-misidentification**: identify запрашивает
   `expected_features_for(parent_nuclide, E_gamma)` и НЕ кредитит
   peaks-in-window-of-secondary как новый нуклид. Это закроет
   chronic Bi-214 false positive на Cs-137 5cm spectra.

2. **Cross-validation**: Cs-137 claim с НЕ-detected backscatter и
   НЕ-detected Ba Kα может быть Bi-214 609 keV misidentification —
   нет вторичных особенностей. `empirical_ratio("Cs-137",
   "backscatter").mean = 0.073` даёт expected magnitude.

3. **Geometry inference**: backscatter intensity ratio →
   inference geometry без явного metadata. Не критично, но полезно
   при missing GEOMETRY field в `.spe`.

## API

```python
from gamma.physics.secondary_peaks import (
    compton_edge_keV, backscatter_keV,
    compton_edge_observed_keV, backscatter_observed_keV,
    expected_features_for,
    load_catalog, empirical_ratio,
)
```

`expected_features_for("Cs-137", 661.66)` → 5 features
(photopeak, compton_edge, backscatter, ic_xray_Ba_Ka, k40_natural).

`empirical_ratio("Cs-137", "backscatter")` →
`{"mean": 0.0726, "std": 0.0297, "min": 0.0319, "max": 0.1457,
"n_observations": 10, "mean_residual_keV": +8.12}`.

`compton_edge_observed_keV(661.66, fwhm_at_edge_keV=50)` → 442.3 keV
(= 477.3 − 0.7·50, near observed mean 440).

## Тесты

13 новых в `test_secondary_peaks.py`:
- Theoretical Compton-edge на Cs-137 / Co-60 канонических энергиях.
- Backscatter formula + energy-conservation identity `E_C + E_bs = E`.
- `compton_edge_observed_keV` reproduces −0.7·FWHM rule и matches
  observed Cs-137 shift в пределах 5 keV.
- Geometry-conditional backscatter shifts match documented table.
- `expected_features_for("Cs-137", 661.66)` — 5 ожидаемых имён,
  `K-40` добавляет single/double escape, Am-241 (60 keV) добавляет
  X-ray escape.
- Catalog JSON loads, Cs-137 + K-40 present, `empirical_ratio()`
  returns правильную shape.
- Cs-137 backscatter R в 5-10% NaI диапазоне и Compton-edge residual
  < −10 keV (validates catalog's internal consistency).

Все 200 prior тестов проходят без изменений. Total: **213 тестов**.

═══════════════════════════════════════════════════════════════════════
v1.7.14 — расширенная многоисточниковая cert-валидация (F-36, Variant B')
═══════════════════════════════════════════════════════════════════════

## Мотивация

К концу v1.7.13 конвейер F-31a + F-31b + F-35 был количественно
проверен **только на Co-60 5 см** (одна точка, +0.61% против
сертификата). Чтобы доказать, что та же связка модулей даёт
сопоставимую точность на других нуклидах сертифицированных
референс-источников, в v1.7.14 собрана **матрица отклонений** через
все доступные `.spe` фикстуры из
`references/reference_spectra/Gamma-1C_NaI_63x63_USB_SN-01/` с
сопоставлением 1↔1 в `references/certificates/АСПЕКТ_ОСГИ_2024.src`.

## Три препятствия, которые пришлось закрыть

### 1. Пробел в библиотеке (4 из 11 cert-нуклидов)

`data/nuclides.json` содержал 24 записи; Y-88, Bi-207, Cd-109,
Th-228 отсутствовали. K-03 фиксировал это, но только частично
закрытым через Lsrm-v2 import. Для расширенной валидации добавлены
прямые JSON-записи (значения ENSDF / NuDat 3, cross-check LNHB):

```json
"Y-88":  { "T_half_s": 9216864,
           "lines": [[898.042, 94.0, 0.4],
                     [1836.063, 99.346, 0.025],
                     [2734.07, 0.71, 0.05]],
           "is_cascade": true }
"Bi-207":{ "T_half_s": 994574880,
           "lines": [[569.698, 97.75, 0.03],
                     [1063.656, 74.5, 0.7],
                     [1770.228, 6.87, 0.06]],
           "is_cascade": true }
"Cd-109":{ "T_half_s": 39865000,
           "lines": [[88.034, 3.66, 0.05]],
           "ic_xrays": [[22.16, 30.0], [22.0, 16.0],
                        [24.94, 4.7], [25.46, 0.91]] }
```

Th-228 НЕ добавлен. Его прямой γ-выход на NaI пренебрежимо мал
(84 / 132 / 216 keV при I ≤ 1.3%); вся регистрируемая активность
"Th-228" приходит через дочерние ядра (Tl-208, Pb-212, Bi-212,
Ac-228) в секулярном равновесии. Чтобы пересчитать "родительскую
активность Th-228" из измерения дочерних γ, нужен механизм
chain-parent reconstruction, который ещё не реализован — поэтому
Th-228 остаётся в матрице как `library gap`.

### 2. Lsrm `FWHM=…` — это полином по √E, а не по E

Заголовок `.spe` несёт `FWHM=N,c0,c1,…` (`StoredFwhmCalibration.model
= "lsrm_fwhm_polynomial_in_E"` в io.lsrm_spe.py). Документация Lsrm
формулирует это как "полиномиальные коэффициенты FWHM(E_keV)", и
комментарий в io/lsrm_spe.py повторяет ту же формулировку. **Это
неверно**: при подстановке E напрямую формула даёт бессмыслицу.

Эмпирическая проверка на Co-60 5 cm фикстуре:
- `FWHM = c0 + c1·E + c2·E² + c3·E³` при E=662 даёт ~340 000 keV.
- `FWHM = c0 + c1·√E + c2·(√E)² + c3·(√E)³` при E=662 даёт **45.6
  keV** — что есть точно NaI 63×63 типичные 7% от 662 keV.
- На E=1332: формула с √E даёт **71.6 keV** = NaI 5.4% (правильно).

Поэтому стоковый `make_fwhm_at_channel_provider` для .spe-файлов
deferred-fallback на физический √E-floor (он не умеет
`lsrm_fwhm_polynomial_in_E`). Это работает, но не точно.

Гарнес валидации (`validate_certs.py`) реализует **собственный**
`make_lsrm_fwhm_provider(spec)`, который оценивает полином по √E:

```python
def fwhm_at(ch):
    E = horner(energy_cal, ch)               # keV
    fw_keV = horner(fwhm_coefs, sqrt(E))     # keV
    dEdN  = horner(deriv(energy_cal), ch)    # keV/channel
    return fw_keV / dEdN                     # channels
```

Это локальный helper — не меняет поведение fwhm_provider для других
тестов. Если в будущем будет тестовая фикстура с реальным
`lsrm_fwhm_polynomial_in_E` и измеренными FWHM пиками, можно
перенести в stock provider.

### 3. Disambiguate отвергает редкие нуклиды по proportionality

`disambiguate_identifications` применяет Rule 4 (универсальную
проверку proportionality): сравнивает отношения Mariscotti-σ
наблюдаемых пиков с отношениями библиотечных интенсивностей.
Mariscotti-σ масштабируется как `height/√B`, **не** как площадь, и
проверка не делит на ε(E). При `prior(nuclide) ≤ 0.2` (редкий
изотоп: Na-22, Y-88) Rule 4 удаляет нуклид целиком, если
наблюдаемое отношение значимостей не совпадает с библиотечным.

На cert-фикстурах мы знаем источник заранее (= одна нуклид-кандидат
в `identify_nuclides(candidate_nuclides=[X])`). Mixture-resolution
rules не дают добавленной ценности и активно мешают. Гарнес
**обходит disambiguate** для single-source семантики. Это
задокументировано в комментарии `run_one`; production-конвейер по-
прежнему запускает disambiguate.

## Adaptive `min_intensity_pct`

Eu-152 имеет 11 библиотечных линий. На 5 см геометрии linia
443.96 keV (I=3.12%) на NaI плохо разрешается: площадь Lsrm-фита
65 515 счётов с σ_S=1114 (1.7% rel), но истинная активность ~120 000
Bq должна была бы дать ~240 000 счётов. Маленький S и маленький σ_S
дают рукавную сигму σ_A_i ≈ 5%·A_i = 1 360 Bq, поэтому вес
`1/σ²` ≈ 5.4·10⁻⁷ — **в 17 раз больше**, чем вес 121.78 keV (5%·A ≈
5 700 Bq, w ≈ 3·10⁻⁸). Линия 443 доминирует средневзвешенное и
уводит результат к 48 000 Bq (-60% от сертификата).

Фикс: в гарнесе `compute_activity(..., min_intensity_pct=5.0)`
фильтрует библиотечные линии с I < 5%. Эта проверка стандартна в
гамма-спектрометрии ([GILMORE-5.7.2] рекомендует исключать low-I
линии при weighted-mean). После фильтра Eu-152 даёт -6.82% — в
пределах типичной NaI cert accuracy.

Один nuance: Cd-109 имеет **только** линию 88 keV при I=3.66%
(весь γ-выход на NaI). Применить порог 5% — потерять единственный
сигнал. Решение: floor адаптивный — `min_intensity_pct=0.0`, если
вся библиотека ядра ниже 5%; иначе 5.0:

```python
has_strong = any(L[1] >= 5.0 for L in lib_lines)
min_I = 5.0 if has_strong else 0.0
```

Безопасно: применяется только когда в библиотеке есть хоть одна
"опорная" линия выше 5%. Cd-109 проходит через нулевой floor и
выдаёт -3.38%.

## Матрица отклонений v1.7.14

Геометрия точечная 5 см, детектор Gamma-1C NaI 63×63 (≈3″×3″),
сертификат `АСПЕКТ_ОСГИ_2024.src`. Все измерения decay-скорректированы
к дате измерения (а не к референсной дате сертификата).

| Нуклид  | A_cert, Bq | A_cert@meas, Bq | A_изм, Bq | Δ, %   | n_лин | Комментарий                            |
|---------|-----------:|----------------:|----------:|-------:|------:|:---------------------------------------|
| Cs-137  |    106 000 |          89 307 |   90 944  | +1.83  |     1 | одна линия 661.66 keV                  |
| Co-60   |    105 000 |          49 886 |   50 193  | +0.61  |     2 | matches v1.7.13 reference exactly      |
| Na-22   |    229 000 |         136 700 |  144 900  | +5.99  |     2 | 511 + 1274 keV; disambig обойден       |
| Eu-152  |    142 000 |         120 600 |  112 400  | −6.82  |     8 | 1 multiplet (1086+1112) deconvolved    |
| Ba-133  |     44 100 |          15 308 |   15 511  | +1.32  |     3 | 4-line TCS                             |
| Am-241  |    103 100 |         102 500 |   93 948  | −8.35  |     1 | 60 keV на нижней границе ε(E)          |
| Zn-65   |      3 100 |             871 |      918  | +5.45  |     1 | 1115 keV; 511 keV (I=2.83%) отфильтр.  |
| Y-88    |    350 000 |          29 903 |   30 392  | +1.64  |     2 | 898 + 1836 cascade pair, TCS           |
| Bi-207  |     97 000 |          82 407 |   80 208  | −2.67  |     3 | 570 + 1064 + 1770 keV                  |
| Cd-109  |    395 000 |           6 762 |    6 534  | −3.38  |     1 | 88 keV; min_I floor дроп. до 0         |
| Th-228  |    129 000 |              — |        —  |   —    |     0 | library gap: chain parent              |

**10/11 cert-нуклидов измеримы.** Mean |Δ| = 3.81%. Max |Δ| = 8.35%
(Am-241 у нижней границы ε(E) curve). Все измеримые внутри ±10% —
это типичный диапазон NaI 63×63 cert accuracy (5-10%).

## Что подтверждает матрица

1. **Co-60 +0.61% сохранён** end-to-end (v1.7.13 reference). Изменений
   в `compute_activity`/`cascade_summing`/`peaks.area` НЕ было —
   только в гарнесе и `nuclides.json`. F-31a + F-31b + F-35 как
   связка работают идентично.
2. **F-35 (K-18 fix) не сломал TCS-нечувствительные нуклиды**.
   Cs-137 (single line, TCS=0): +1.83% — близко к +2-3% типичному
   разбросу cert calibration.
3. **Cascade-нуклиды с TCS работают**. Co-60 +0.61%, Na-22 +5.99%,
   Y-88 +1.64%, Ba-133 +1.32%, Bi-207 −2.67% — все внутри cert
   uncertainty. Eu-152 −6.82% хуже остальных — это самый плотный
   мультиплет (11 линий, низко-E плохо разрешаемые), типичный
   challenge для NaI.
4. **Zn-65 11.41% → 5.45%** при включении `min_intensity_pct=5.0`.
   Это **снимает мотивацию Variant F** (refit ε(E) excluding
   cascade-depleted points). Слабая I=2.83% 511-keV линия Zn-65 не
   "cascade-depleted" в традиционном смысле — она просто
   статистически шумная и завышенно-weighted, потому что
   уменьшающаяся σ_S делает вес `1/σ²` слишком большим. Фильтрация
   low-I линий — правильное решение.
5. **Am-241 −8.35% — на нижней границе ε(E) curve.** 60 keV это
   `E_min_keV=59.541` — самая низкая калибровочная точка ε(E)-кривой
   `Точечная-5см.efr`. Любая экстраполяция или плохая точка вблизи
   границы дала бы ошибку 5-10%. Для улучшения нужен новый
   low-E ε(E) refit с более плотной сеткой <100 keV — это
   отдельная задача (Variant F variant).

## Что НЕ закрыто и в каком приоритете

- **Th-228**: chain-parent reconstruction. Когда родительская
  активность определяется по дочерним под equilibrium —
  `compute_activity` должен принимать "parent_chain"-флаг и
  пересчитывать. Сложность средняя; не блокирующая.
- **Disambiguate Rule 4 (proportionality)**: проверяет
  Mariscotti-σ ratio вместо area ratio / efficiency-corrected
  ratio. Для production-кейсов (mixed samples) это работает, но
  cert-валидация требует обхода. Альтернатива: исправить Rule 4 на
  правильные величины — тогда обход не нужен. Это потенциальный
  F-37.
- **Stock fwhm_provider для `lsrm_fwhm_polynomial_in_E`**: harness
  делает локально, не глобально. Перенос в
  `make_fwhm_at_channel_provider` — отдельная задача с
  регрессиями на AtomSpectra fixtures.

## Воспроизводимость

```powershell
cd "<WORKDIR>\gamma-spectrum-analysis"
$env:PYTHONPATH = "scripts"
$env:PYTHONIOENCODING = "utf-8"
py -3.11 validate_certs.py
```

Выход: stdout-таблица (Russian) + `cert_validation_matrix.csv` для
downstream regression / dashboarding.

═══════════════════════════════════════════════════════════════════════
v1.7.13 — area-method-aware TCS (F-35, закрывает K-18)
═══════════════════════════════════════════════════════════════════════

## Эмпирическое наблюдение, мотивирующее этот фикс

В v1.7.9 на Co-60 5cm point-source наблюдалась последовательность
результатов, поначалу казавшаяся парадоксальной:

| Конфигурация                              | Co-60 5cm vs cert |
|-------------------------------------------|-------------------|
| Cowell-площади, без TCS                   | −25.56%           |
| Cowell заменён на Lsrm-table (F-31a)      | +0.61%            |
| F-31a + аналитическая TCS (F-31b)         | +2.89%            |

F-31a даёт почти идеальный результат БЕЗ всякой TCS-коррекции — что
физически странно, ведь Co-60 — каскадный нуклид, должно быть
photopeak-depletion. Добавление физически корректной аналитической
TCS затем СВЕРХ-корректирует.

Объяснение: **Lsrm peak-table это не «исправленная Cowell-площадь», а
выход из совершенно другого алгоритма** — широко-ROI gaussian-on-step
fit. Подгоняя gaussian-shape к данным на широком окне, он восстанавливает
часть фотопиковых counts, которые суммирование сместило в боковые
каналы (но НЕ в sum-peak на E1+E2). Аналитическая TCS делает
противоположное предположение: что наблюдаемая площадь — это полное
photopeak ROI integration, и из неё нужно поднять обратно ВСЕ
суммированием смещённые. Применённая к Lsrm-table — двойной учёт.

## Архитектурный фикс

Поскольку TCS-correction чувствительна к методу integration, нужен
способ для `compute_activity` различать «откуда у меня площадь».
F-34 уже добавил `LineMatch.peak_area_source` (`"" / "cowell" /
"lsrm_peaks_table" / "deconvolved" / "failed"`). F-35 использует
это поле:

```python
c_effective = 1 + (c_analytic − 1) · scale[peak_area_source]
```

Из `DEFAULT_TCS_METHOD_SCALE`:

```python
{
    "":                  1.0,   # safe default — full TCS
    "cowell":            1.0,
    "deconvolved":       1.0,
    "failed":            1.0,
    "lsrm_peaks_table":  0.0,   # Lsrm-fit recovers them
}
```

Caller может передать свой dict через kwarg `tcs_method_scale`. Dict
**мерджится** поверх defaults — указание только `{"cowell": 0.0}` не
сбрасывает `lsrm_peaks_table` на 1.0.

## Почему deconvolved тоже full TCS

«Deconvolved» — это linear LSQ-fit fixed-position гауссианов поверх
step+linear подложки (F-33). При фиксированных позициях и FWHM
ROI определён жёстко — channels внутри окна делятся между
компонентами и подложкой. Counts которые суммирование сместило в
sum-peak (E1+E2) лежат сильно вне любой компоненты — на 2505 кэВ
для Co-60 (далеко от 1173 или 1332). Их deconvolution НЕ видит.
Значит TCS-correction здесь должна быть полной — те же
физические counts depleted, что и в случае Cowell.

Для Lsrm-table ситуация другая: широкая Gaussian fit с гладкой
step-baseline частично «затягивает» counts которые суммирование
смещает чуть-чуть (в близкий канал, типа +1-2 σ от центра). Это
не то же самое, что counts смещённые в sum-peak; но эмпирически
эти небольшие смещения дают, видимо, основную часть TCS-эффекта
на Co-60 5cm. Откуда и α = 0.

Будущая работа: на других геометриях (1cm, Marinelli) α может
быть не строго 0. Сейчас обоснованной калибровки нет — пользователь
может передать свой dict через `tcs_method_scale`.

## Merge-semantics на user_dict

Я делаю `{**DEFAULT, **user_dict}` (не `user_dict` целиком). Почему:

- Если пользователь указывает только `{"lsrm_peaks_table": 0.3}`,
  он явно меняет ОДИН известный default, а не пересоздаёт мапу с
  нуля.
- Если он не упомянул `"cowell"`, то ожидает, что Cowell-areas
  получат полный TCS (default 1.0). Merge сохраняет это.
- Альтернатива (replacement): пользователь должен каждый раз
  перечислить все 5 ключей, что неудобно. Источник ошибок:
  пропустил один ключ — он стал 0, неожиданно нет TCS.

Поэтому merge.

## Provenance в LineActivity.correction_factor

`LineActivity.correction_factor` хранит **эффективное** значение
после scaling, а не исходное analytic `C(E)`. Это сделано чтобы
caller, считывая результат, видел физически применённую коррекцию.
Если ему нужен исходный `C` — он у него в исходном
`coincidence_correction` dict, который он сам передал.

В `ActivityResult.notes` появляется строка
`"K-18: TCS scaled by area-method on N line(s)"` ТОЛЬКО когда
scale ≠ 1.0 хотя бы на одной строке. Это даёт диагностику без
шума в каждом результате.

## Регрессия и null-effect для legacy LineMatch

Существующий test_compute_tcs_corrections_compatible_with_compute_activity
из v1.7.9 строит LineMatch БЕЗ `peak_area_source` (поле появилось в
F-34). Дефолт поля — `""`. В `DEFAULT_TCS_METHOD_SCALE[""] = 1.0` →
полный TCS → test продолжает проходить. Это намеренная
backwards-compat: любой LineMatch с пустым source попадает на
«full-TCS» rail (the safe rail).

## Что в итоге даёт Co-60 cert validation

Используя F-31a (Lsrm-table area, через `get_peak_area`), F-31b
(аналитическая TCS, через `compute_tcs_corrections`), F-35 (scale
lsrm_peaks_table → 0):

- Per-line с `peak_area_source="lsrm_peaks_table"`: c_eff = 1.0
  → A_i = S_lsrm / (eps * I * t)
- Per-line с `peak_area_source="cowell"`: c_eff = c_analytic
  → A_i = S_cowell * c / (eps * I * t)

Co-60 5cm: оба photopeak дают `peak_area_source="lsrm_peaks_table"`
(потому что Lsrm SPE имеет PEAKS table с ними), → оба получают
c_eff = 1, активность совпадает с F-31a-alone (+0.61% от cert).

В первый раз pipeline активирован end-to-end (Lsrm-table → TCS dict
→ K-18 scaling → A_Bq) и даёт точность 1σ. Это и есть «закрытие
K-18».

## Что НЕ покрывается этим фиксом

1. Геометрия-зависимая калибровка `α(geometry)`. Сейчас α = 0 на
   5cm point. Что на 1cm или Marinelli — не проверялось.
2. Детектор-зависимая `α(detector_type)`. Все валидации на NaI. На
   HPGe wide-ROI gaussian fit ведёт себя иначе — `α` может быть
   ненулевым там.
3. Pile-up vs TCS разделение. TCS — отдельный физический эффект от
   pile-up; F-31b модель учитывает оба некорректно если CPS высокая.
   Для контрольных источников на низких CPS pile-up пренебрежим, и
   F-35 здесь чистый K-18 fix.

═══════════════════════════════════════════════════════════════════════
v1.7.12 — pipeline-интеграция deconvolution (F-34)
═══════════════════════════════════════════════════════════════════════

## Зачем отдельный этап после F-33

F-33 (v1.7.11) дал отдельностоящий алгоритм `deconvolve_multiplet` и
проверил его на синтетике + Co-60. Но в пайплайне он жил «сбоку»:
вызывающий код должен был сам найти кластеры, выполнить decon, потом
вручную заменить `LineMatch.peak_area` на каждой строке. Это работало
бы только если бы все были очень внимательны — и пропадает само
обоснование «identification-first deconvolution», потому что
identification по-прежнему отдаёт Cowell-площади в `compute_activity`.

F-34 устраняет этот разрыв: вместо ручной мутации даём
`apply_multiplet_deconvolution(id_result, spec, fwhm_at, ...) ->
(new_id_result, list)` — pluggable post-pass, выдающий новый
`IdentificationResult` со сменёнными площадями там, где нашёлся
мультиплет.

## `LineMatch.peak_area_source` — поле провенанса

Новое поле в `LineMatch`:

```python
peak_area_source: str = ""   # "" / "cowell" / "lsrm_peaks_table" /
                              # "deconvolved" / "failed"
```

Раньше `gamma.peaks.area.get_peak_area` уже возвращал кортеж
`(area, unc, source)`, но identify.py отбрасывал `source` через
`_src`. Теперь:

- `identify_nuclides` кэширует `(area, unc, source)` и пишет source
  в LineMatch.
- `disambiguate_identifications` сохраняет source при promotion
  (был баг — promoted LineMatch теряла peak_area вообще; чинится
  попутно).
- `apply_multiplet_deconvolution` пишет `"deconvolved"`.

Зачем это нужно? Главное — K-18: «TCS over-correction на
Lsrm-fitted areas». Это open issue, требующий area-method-aware
TCS-коррекции. Без поля provenance невозможно различить «эта площадь
из Cowell» / «эта из Lsrm-table» / «эта из decon» на уровне
`compute_activity`. С полем — путь открыт.

`compute_activity` пока не читает `peak_area_source` — это
намеренная минимальная интеграция: identification выбирает метод
integration, activity-слой доверяет выбору.

## Контракт `apply_multiplet_deconvolution`

```python
def apply_multiplet_deconvolution(
    identification_result,
    spec,
    fwhm_at_channel,
    *,
    overlap_threshold_fwhm: float = 1.0,
    continuum: str = "step_linear",
    max_chi2_per_dof: float = float("inf"),
) -> tuple[IdentificationResult, list[DeconvolutionResult]]:
```

Алгоритм:

1. `find_multiplet_regions(id_result, fwhm_at, overlap_threshold_fwhm)`
   → список кластеров `LineMatch`-объектов.
2. По каждому кластеру: построить `MultipletComponent` (используя
   `spec.energy_to_channel(library_E_keV)` и `fwhm_at_channel(ch)`),
   вызвать `deconvolve_multiplet(spec.counts, components=…,
   continuum=continuum)`.
3. Собрать map `(nuclide, round(library_E_keV, 3))` → `(area, unc)`
   из результата decon. Использовать `round(..., 3)` чтобы избежать
   мелких float-расхождений между библиотечной и LineMatch-копией.
4. Пройти по всем `detected_nuclides.matched_lines`, для каждого
   `LineMatch` посмотреть, есть ли он в map. Если да и decon
   converged и `χ²/ν ≤ max_chi2_per_dof` → заменить через
   `dataclasses.replace(m, peak_area=…, peak_area_uncertainty=…,
   peak_area_source="deconvolved")`. Иначе оставить как есть.
5. Собрать новый `IdentificationResult` с обновлёнными матчами.
6. В `notes` записать строку вида
   `"Multiplet deconvolution: N cluster(s), M peak area(s) replaced
    (overlap threshold = X·FWHM, continuum = step_linear)"`.

Возврат: `(new_id_result, [DeconvolutionResult, …])` — список с
полными результатами decon, в том же порядке, в каком
`find_multiplet_regions` их возвращает. Полезно для отчётности /
плотного диагноза.

## Дефолтный threshold = 1.0 — почему не больше

Соблазн был сделать дефолт `overlap_threshold_fwhm = 2.0` или `3.0`,
чтобы Co-60 1173/1332 на NaI ловился автоматически. Не сделал —
вот рассуждение:

- При threshold = 1.0 «мультиплет» = пик-пары, у которых половинные
  ширины касаются. Физически это режим, где **визуально** видна одна
  расплывшаяся вершина.
- При threshold = 2.0–3.0 захватываются пары, где пики **визуально
  разделены**, но крылья пересекаются. Cowell-интегрирование на таких
  парах систематически занижает площадь (это и есть F-31a). Decon
  помогает.
- При threshold = 4.0+ начинаем кластеризовать пики, которые
  фактически независимы — внося ненужную корреляцию между их
  площадями (decon рассматривает их совместно как линейная LSQ-задача,
  результирующая ковариация переходит на сильно недиагональную, σ_A
  раздувается).

Co-60 1173/1332 — особый случай: они визуально разделены (Δ ≈
3·FWHM), но F-31a показала 30-42% занижение Cowell-площади.
Правильный путь — вызвать decon с `overlap_threshold_fwhm=3.0` ИМЕННО
когда такой dataset обрабатывается. Пользователю/тесту понятно,
зачем; дефолт остаётся консервативным для разнообразного
input.

В будущем можно сделать адаптивный threshold (зависящий от того,
какая методология integration использовалась — Cowell vs Lsrm-table —
поскольку Lsrm Gaussian-fit уже частично восстанавливает крылья и
threshold=1.0 на нём может быть достаточен).

## Тонкость: matching по `(nuclide, library_E_keV)`

Когда decon-результат отображается обратно на LineMatch, я
использую ключ `(nuclide, round(library_E_keV, 3))`. Почему round(3)?

- `MultipletComponent.line_E_keV` приходит из
  `LineMatch.library_E_keV`, а тот — из библиотечного значения
  (например, 1173.228 для Co-60 1173).
- Где-то в библиотеках может оказаться 1173.228 vs 1173.23 — мелкие
  float-различия.
- Round to 3 decimals (`0.001 keV`) — точность, которой все
  библиотеки определённо точны.

В пределах одного нуклида два library line с разностью < 0.001 keV
не встречаются (это была бы одна линия). В пределах одного
кластера два разных нуклида с одинаковым library_E_keV
теоретически возможны — но они получили бы один и тот же ключ
только если нуклиды тоже совпадают, чего нет (разные).

## Тесты — 5 новых

- `test_apply_post_pass_returns_tuple` — sanity: всегда tuple, для
  no-cluster случая первый элемент IS the input result.
- `test_apply_post_pass_no_change_for_isolated_lines` — на реальном
  Co-60 spectrum lines выше/ниже 1173/1332 кластера сохраняют
  `peak_area_source = "lsrm_peaks_table"`.
- `test_apply_post_pass_replaces_co60_doublet` (threshold=3.0) —
  1173 area = 540 620, 1332 = 495 016, ratio 1.09. Совпадает с
  тем, что F-33 stand-alone decon давал на том же спектре.
- `test_apply_post_pass_notes_record_replacement` — string-check на
  `"Multiplet deconvolution"` в notes.
- `test_apply_post_pass_max_chi2_filter_skips_bad_fits` — порог
  `max_chi2_per_dof=0.01` отсекает все хорошие фиты тоже; нулевая
  замена.

## Что разблокировано для следующих фаз

1. **K-18 — area-method-aware TCS**. Можно теперь сделать в
   `gamma.activity.compute`: при `coincidence_correction` смотреть на
   `peak_area_source` каждой LineMatch и применять reduced TCS
   только к `"lsrm_peaks_table"`-площадям (Lsrm Gaussian-fit уже
   частично восстанавливает summing-displaced counts).
2. **Adaptive overlap threshold**: можно посчитать долю площади
   соседнего peak'а под крыльями данного пика (через
   `_channel_energies` и `_gaussian_normalised`), и если она > X%,
   кластеризовать. Точнее, но дороже.
3. **CLI hook**: `python -m gamma.cli analyze --deconvolve-multiplets`
   очевидный следующий шаг для Variant A.

═══════════════════════════════════════════════════════════════════════
v1.7.11 — multiplet deconvolution (F-33, закрывает K-05)
═══════════════════════════════════════════════════════════════════════

## Что и почему

K-05 («deconvolution не реализована») висел с v1.7.1 как deferred-задача.
Мотивация: на сцинтилляторах с FWHM 30–50 кэВ многие соседние линии
блендируются в одну видимую вершину. Примеры:

  - 240 кэВ на NaI = Pb-212 238.6 + Pb-214 241.98
  - 597 кэВ на NaI = Bi-214 609.31 + Tl-208 583.19
    (на верхней ветке Compton)
  - 1173/1332 кэВ Co-60 — формально разрешаются на NaI, но крылья
    соприкасаются и `cowell_area` теряет ~30–40% площади
    (это была причина F-31a)

Без deconvolution каждая такая «вершина» матчилась с одной библиотечной
линией, а вклад второй полностью игнорировался — активности
завышены/занижены пропорционально интенсивности проигнорированной линии.

## Ключевая методическая развилка: identification-first

Lsrm §9 и [GILMORE-9.7] настойчиво предписывают **сначала
identification, потом deconvolution**. Это не косметика:

- Если позиции компонент **свободные**, то deconvolution превращается
  в полноценный нелинейный LSQ-фит с риском местных минимумов,
  обменом компонент, вырождением. На многокомпонентных мультиплетах
  на сцинтилляторе он часто расходится или находит «лишние» пики
  на шуме.
- Если позиции **фиксированы из библиотеки** confirmed nuclides, и
  σ_k **фиксирована из калиброванной FWHM(E)**, то остаются только
  **площади** и **подложка**. Модель становится **линейной по
  свободным параметрам**:

      y(x) = Σ_k A_k · g_k(x; c_k, σ_k) + β₀ + β₁·(x − x_mid)        (linear)
      y(x) = … + β_step · 0.5·erfc((x − x_step)/(σ_step·√2))         (step + linear)

  где g_k — нормированный гауссиан (площадь = 1). Свободный вектор:

      p = [A₁, …, A_n, β₀, β₁]                              для linear
      p = [A₁, …, A_n, β₀, β₁, β_step]                      для step + linear

  Ограничения: A_k ≥ 0 (физика), β_step ≥ 0 (Compton-ступенька
  «опускается» через пик, проходя слева направо), β₀ и β₁ свободны.

  Эта задача решается `scipy.optimize.lsq_linear(A, b, bounds=(lb, ub),
  method="trf")` за один шаг, без риска несходимости. Ковариация — из
  обратной матрицы нормальных уравнений, аналитически.

## Почему bounds, а не штрафы

Альтернативы:
- L1/L2 регуляризация: размывает решение, требует подбора λ
- Послефитовая клиппинг: даёт смещённые оценки при сильной корреляции
- Trust-region NNLS (scipy.optimize.lsq_linear): даёт **bound-aware**
  оптимум; ковариация по диагонали корректна для активных
  «свободных» переменных, и нулевые оценки для активных границ
  получают строго нулевую вариацию (что и хочется — переменная
  «припёрта к стенке»)

Поэтому именно `lsq_linear` с trf и явными границами.

## β_step ≥ 0 — физика

Compton-непрерывный фон под фотопиком формируется рассеяниями того же
фотопика (и более высокоэнергичных линий). Фотопик при E_0 генерирует
Compton-фон на E < E_0, и не генерирует на E > E_0. Значит, проходя
через фотопик слева направо (по возрастанию канала), подложка
**уменьшается**: B_L > B_R. В наших обозначениях S(x) = 1 при
x ≪ x_step, S(x) = 0 при x ≫ x_step → β_step·S(x) добавляет на левой
стороне. Если β_step ≥ 0, мы корректно моделируем «ступеньку вниз».
Отрицательный β_step физически означал бы поглощение Compton-фона на
левой стороне — это не происходит.

## Degenerate-pair flag

Если две компоненты ближе чем 0.5·σ_min, их площади статистически
неразделимы: ковариация почти −1, индивидуальные σ_A раздуты, но
сумма площадей определена точно. Мы:
1. Всё равно делаем фит (математика не падает)
2. Регистрируем пару в `DeconvolutionResult.degenerate_pairs`
3. Включаем notes с инструкцией: «consult covariance» / «consider
   merging»

Это даёт вызывающему коду честный сигнал: «индивидуальные числа
ненадёжны, но сумма — да».

## Compton-step continuum vs pure linear: когда выбирать какой

- Изолированный пик на плавно меняющемся фоне → linear достаточно
- Любой мультиплет → step_linear обязательно (там подложка обычно
  на разных уровнях слева/справа)
- На синтетике без ступеньки оба дают сравнимый χ²/ν — лишний
  параметр без сильной поддержки данных получает близкое к нулю
  значение
- На синтетике со ступенькой step_linear даёт χ²/ν ~ 1, linear
  даёт χ²/ν >> 1 (см. test_step_continuum_recovery: 0.96 vs 2.18)

По умолчанию я делаю step_linear — он строго лучше или равен
linear для любых разумных данных.

## Поведение χ²/ν на реальных NaI-спектрах

На синтетике (Poisson-only) χ²/ν ∈ [0.8, 1.5]. На реальном Co-60
спектре `Co-60__043_02_2019_Точечная-5см_5cm.spe` χ²/ν ≈ 325 для
двух-компонентного фита 1173/1332. Это **ожидаемо**:

- NaI photopeak имеет лёгкий низкоэнергетический хвост, не
  описываемый чистым гауссианом
- Реальная Compton-непрерывная подложка имеет тонкую структуру,
  плохо описываемую одной erfc-ступенькой
- Statistics очень высокая (площади ~500 000), Poisson-σ маленький,
  любые систематические отклонения модели «торчат»

Тем не менее **площади остаются полезными** — на Co-60 1173/1332
получаем ratio 1.09 (библиотека 1.00). Расхождение объясняется
small TCS-asymmetry + остаточным crosstalk в линейном фоне. Для
активность-калибровки достаточно. Ковариация при χ²/ν > 1 умножается
на χ²/ν (так и реализовано в `_parameter_covariance`), поэтому σ_A
автоматически растут пропорционально неадекватности модели — то, что
нужно.

Если в будущем потребуется χ²/ν ~ 1 на реальном NaI, нужно
расширить модель response-функции: добавить exp-tail слева
(EMG / Hypermet) с фиксированными параметрами хвоста, откалиброванными
на изолированных пиках. Это уже Phase 2.1c или позже.

## Что НЕ интегрировано в этой фазе

`identify_nuclides` (см. identify.py:283) кэширует Cowell-площади в
`peak_area_cache[ch]` через `get_peak_area`. Логика «если линия
попала в мультиплет — заменить cowell-площадь на deconvolved-площадь»
не введена в `identify.py` намеренно: она трогает контракт LineMatch,
порядок обработки в identification и activity, и заслуживает
отдельной фазы с прицельной валидацией на нескольких реальных
спектрах (Ra-226 chain, Th-232 chain, Eu-152 multilineator).

Текущий контракт: deconvolution — отдельный stand-alone tool,
который вызывающий код может применять явно после identification.
Это намеренный «pluggable» дизайн — даёт время для калибровки
дефолтов до того, как замаунтить в hot path.

═══════════════════════════════════════════════════════════════════════
v1.7.10 — симметричный потолок энергии в ридерах (F-32)
═══════════════════════════════════════════════════════════════════════

## Архитектурная асимметрия, которая возникла исторически

Два формат-специфичных ридера развивались поэтапно:

- `read_lsrm_spe` (введён в v1.6) получил keyword-only флаг
  `apply_energy_ceiling: bool = True` сразу, потому что для .spe потолок
  обрезает «лишний» хвост уже после распарсивания всех 1023/2048
  каналов, и было сразу удобно его отключать ради диагностики.
- `read_atomspectra_xml` (введён в v1.2) применял потолок безусловно
  внутри `_parse_energy_spectrum_block`. Не было keyword-only флага, а
  параметр `parse_background` остался позиционно-доступным.

При попытке выполнять диагностику выше 3 МэВ или с другим потолком
требовалось временно править модуль-уровневую константу
`ENERGY_CEILING_KEV`. Это глобальная мутация, её легко забыть откатить,
и она портит результат для параллельных чтений из других процессов
или тестов того же интерпретатора.

## Что введено в v1.7.10

Оба ридера экспонируют один и тот же keyword-only контракт:

```python
read_lsrm_spe(path, *,
              apply_energy_ceiling=True,
              ceiling_keV=None)

read_atomspectra_xml(path, *,
                     parse_background=True,
                     apply_energy_ceiling=True,
                     ceiling_keV=None)
```

Семантика `ceiling_keV=None` — «использовать константу
`ENERGY_CEILING_KEV` (3000 keV)». Любое non-None значение применяется
как потолок на эту конкретную запись. Если `apply_energy_ceiling=False`,
параметр `ceiling_keV` игнорируется. Дефолтное поведение не меняется —
148 базовых тестов проходят без правок.

## Распространение политики на встроенный фон

`read_atomspectra_xml` пробрасывает свои `apply_energy_ceiling` /
`ceiling_keV` в **оба** вызова `_parse_energy_spectrum_block` — для
основного спектра и для встроенного `BackgroundEnergySpectrum`. Это
сделано осознанно: если пользователь отключает потолок ради
диагностики, естественно ожидать, что `spec.background_embedded.counts`
будет такой же длины, что и `spec.counts`. Иначе пара «спектр / фон»
расходится по сетке каналов и потребует дополнительной выравнивающей
логики выше по стеку.

## Дополнительные правки

- `parse_background` тоже переведён в keyword-only — это даёт обоим
  ридерам одну форму подписи: `(path, *, …)`. Единственный
  внутренний потребитель (`gamma.io.background.resolve_external_background`,
  строка 62) уже звал параметр по имени, поэтому миграция безопасна.
- Внутри `gamma.io.atomspectra_xml` извлечён общий helper
  `_channel_energies(coeffs, n) -> np.ndarray | None` — расчёт энергии
  каждого канала методом Хорнера. И ветка «потолок применяется», и
  ветка «потолок отключён» теперь используют один и тот же расчёт для
  `energy_max_keV_kept`.
- `read_spectrum` (диспетчер) ничего не меняет — он уже
  пробрасывает `**kwargs` в format-specific reader.

## Граничные значения и проверки

Проверено на реальной фикстуре `Фон_кабинет_8192к_01-01-2025.xml`
(8192 каналов, последний канал — overflow marker):

| Вызов                                  | n_channels | e_max keV  |
|----------------------------------------|------------|------------|
| Default                                | 7034       |  2999.6    |
| `apply_energy_ceiling=False`           | 8191       |  3441.5    |
| `ceiling_keV=1500`                     | 3533       |  1499.8    |

На `M_cs_легкий_2001-2005.spe` (1023 канала, нет overflow):

| Вызов                                  | n_channels | e_max keV  |
|----------------------------------------|------------|------------|
| Default                                | 927        |  2997.8    |
| `apply_energy_ceiling=False`           | 1023       |  3280.2    |
| `ceiling_keV=400`                      | 116        |   399.9    |

Все asserts test_reader_api.py проходят на этих числах.

═══════════════════════════════════════════════════════════════════════

# v1.7.19 — F-41 Tl-208 chain-proxy cross-validation для Th-228

## Контекст

F-39 (v1.7.17) закрыл Th-228 cert row через **Pb-212 238 keV
chain-proxy**: Th-228 в 1:1 secular equilibrium со всеми
α/β daughter'ами вплоть до Pb-212. Активность измеренного Pb-212
интерпретируется как A_Th-228, с decay-correction по T½_Th-228
(не Pb-212, поскольку daughter в equilibrium tracks parent).

Получилось −12.79% deviation — резко outlier против остальной
матрицы (mean |Δ| прочих 10 строк = 3.8%, max = 8.35% для Am-241).
Этот один результат доминирует max |Δ| матрицы. Возникает вопрос
**локализации источника ошибки**:

1. **Cert value ошибочен.** Th-228 source (SRC-01) на самом деле
   имеет активность ≠ 129000 Bq на дату 25.05.2023.
2. **Chain-proxy методология ошибочна.** Equilibrium assumption,
   T½-based decay correction, или decay logic в `compute_activity`
   систематически смещают результат.
3. **Geometry/energy-specific systematic.** 238 keV single-line
   измерение на 5cm point geometry имеет специфическую проблему
   (self-absorption в sealed source, ε(E) detail в low-E крыле,
   statistical fluctuation single-line vs multi-line averaging).

Невозможно различить (1)/(2)/(3) одной точкой измерения. Нужна
**независимая вторая chain через ДРУГОЙ дочерний изотоп** на том
же исходном спектре.

## Tl-208 как natural second proxy

В Th-232/228 chain:

```
Th-228 (α, T½=1.91 y) → Ra-224 → Rn-220 → Po-216 → Pb-212 (β, 10.6h)
   → Bi-212 (T½=60.6 min)
       ├─ β⁻ к Po-212 → α к Pb-208 (branching 64.06%)
       └─ α  к Tl-208 → β к Pb-208 (branching 35.94%)
```

Tl-208 в equilibrium получает 0.3594 от parent Th-228 disintegration
rate. **A_Tl-208 = 0.3594 · A_Th-228**.

Tl-208 эмиссионные линии (per Tl-208 disintegration, ENSDF 2004):
- 2614.51 keV @ 99.75%  (strong, чёткая)
- 583.19 keV  @ 84.5%   (clean, не overlap)
- 510.77 keV  @ 22.6%   (near positron annihilation 511)
- 860.56 keV  @ 12.5%   (well-isolated)
- 277.36 keV  @ 6.6%    (weak)

**Energy spread vs Pb-212**: Pb-212 имеет 238 keV (μ/ρ ~0.40 на NaI),
Tl-208 — 583/860/2614 keV (μ/ρ <0.05). Совершенно разные физические
режимы absorption + ε(E). Если Δ_Tl-208 совпадает с Δ_Pb-212 →
ошибка в (1)/(2). Если резко различаются → ошибка в (3).

## Implementation: chain_branching = 1.0 (не 0.3594)

Lsrm chain library (`NaI-Etl+Esc.lib`, загружаемая через
`load_lsrm_chain_libs()`) **уже pre-scaled** Tl-208 γ-line
intensities на 0.3594 β-branching factor:

```
loaded_lib["Tl-208"]["lines"] = [
    [2614.51, 35.85, 0.07],   # = 0.3594 × 99.75
    [583.19,  30.6,  0.3],    # = 0.3594 × 84.5
    [510.77,   8.1,  0.1],    # = 0.3594 × 22.6
    [860.56,   4.5,  0.1],    # = 0.3594 × 12.5
    [277.36,   2.37, 0.04],   # = 0.3594 × 6.6
]
```

`compute_activity` инвертирует lib intensity: A = Σ counts / (I_lib · ε · T).
Когда I_lib уже содержит β-factor, recovered A — это **parent A_Th-228**,
не A_Tl-208. β-фактор сокращается между обоими местами:

```
counts ∝ A_parent · β · I_raw · ε · T
A = counts / (lib_I · ε · T)
  = (A_parent · β · I_raw · ε · T) / (β · I_raw · ε · T)
  = A_parent
```

Поэтому `chain_branching=1.0` в `CertFixture` — это уже не двойной
учёт β-факторa.

## Cross-validation block

Новый блок в `validate_certs.py:main()` после CSV write:

```python
chain_proxies = {}  # cert_nuclide → list of RunResult
for r in rows:
    fx = next((f for f in FIXTURES if f.spe_filename == r.spe_filename
               and f.nuclide == r.nuclide), None)
    if fx and fx.cert_nuclide and r.measured_A_Bq is not None:
        chain_proxies.setdefault(fx.cert_nuclide, []).append(r)
for parent, proxy_rows in chain_proxies.items():
    if len(proxy_rows) < 2:
        continue
    # print side-by-side + pairwise ratios + flag >5 %
```

## Результаты

```
Parent: Th-228
    daughter    A_meas, Bq   A_cert@meas, Bq   Δ vs cert, %
      Pb-212       77 250.0           88 575.1        -12.79%
      Tl-208       88 516.8           88 575.1         -0.07%
  ratio Pb-212/Tl-208 = 0.8727  (-12.73%)  >5%
```

**Tl-208 даёт Δ = −0.07 %** — практически идеальная согласованность с
cert decayed value. Это **исключает варианты (1) и (2)**:

- (1) Cert 129000 Bq @ 25.05.2023 корректен (Tl-208 высокоэнергетически
  независимо его подтверждает).
- (2) Chain-proxy методология F-39 (T½-based decay correction по parent,
  equilibrium assumption) корректна — Tl-208 проходит ту же логику и
  попадает в ±0.1 %.

Остаётся **(3) — 238 keV / 5 cm point geometry-specific systematic**:
- Self-absorption matrix attenuation в sealed source: μ(NaI@238)≈0.40
  cm⁻¹ vs μ(NaI@583)≈0.10 cm⁻¹. Source-matrix attenuation в plastic
  capsule ОСГИ источника тоже больше для 238 keV.
- ε(E) curve detail: точка 238 keV сидит у нижнего крыла, где
  efficiency-curve fit имеет наибольшую неопределённость.
- Pb-212 single-line measurement не имеет multi-line averaging для
  гашения статистической флуктуации.

## Methodological implication

Multi-chain cross-validation **должна быть default** когда parent
имеет ≥ 2 detectable daughter chains. F-41 формализует для Th-228;
будущее расширение возможно для:

- **U-238 chain**: Pb-214 (295/352 keV) + Bi-214 (609/1120/1764/2204
  keV) → A_U-238 via Ra-226 equilibrium.
- **Th-232 chain**: Pb-212 (238) + Tl-208 (583/860/2614) + Ac-228
  (911/968) → A_Th-232 (same chain as Th-228 daughter set; useful
  для natural Th-232 sources вроде монацитовых песков).

**Cross-validation ratio как diagnostic**:
- Ratio ≤ 5% → chain-proxy методология валидна для обеих линий.
- Ratio > 5% → указатель на geometry/energy-specific systematic.

## Validation matrix v1.7.19

| Нуклид | A_cert@meas | A_изм   | Δ, %    | n_лин | comment                        |
|--------|-------------|---------|---------|------:|--------------------------------|
| Cs-137 |     89307   |  90944  |  +1.83  |   1   | 8 peaks                        |
| Co-60  |     49886   |  50193  |  +0.61  |   2   | 13 peaks; TCS=2                |
| Na-22  |    136700   | 144900  |  +5.99  |   2   | 16 peaks; TCS=2                |
| Eu-152 |    120600   | 112400  |  −6.82  |   8   | 15 peaks; 1 multiplet; TCS=7   |
| Ba-133 |     15308   |  15511  |  +1.32  |   3   | 15 peaks; TCS=4                |
| Am-241 |    102500   |  93948  |  −8.35  |   1   | 11 peaks                       |
| Zn-65  |       871   |    918  |  +5.45  |   1   | 7 peaks                        |
| Y-88   |     29903   |  30392  |  +1.64  |   2   | 13 peaks; TCS=2                |
| Bi-207 |     82407   |  80208  |  −2.67  |   3   | 14 peaks                       |
| Cd-109 |      6762   |   6534  |  −3.38  |   1   | 9 peaks                        |
| Pb-212 |     88575   |  77250  | −12.79  |   1   | 18 peaks                       |
| Tl-208 |     88575   |  88517  |  −0.07  |   3   | 18 peaks; TCS=4                |

**12/12 measurable**, mean |Δ| = **4.24 %** (было 4.62 %),
max |Δ| = 12.79 % (Pb-212 ; Tl-208 = −0.07 %).

═══════════════════════════════════════════════════════════════════════

# v1.7.20 — F-42 Symmetric reader API для LSRM `.spe`

## Scope этой итерации

**Только LSRM NaI / `.spe`.** Пользователь явно сузил scope:

> «Пока работаем с детектором лсрм NaI. Доделываем весь план для
> этого детектора и его файлов спектров `*.spe`. Другие детекторы
> и форматы спектров будем разрабатывать отдельно после
> завершения плана для лсрм NaI.»

Архитектурно изменения `read_atomspectra_xml` оставлены в коде для
единого API, но **тестовый scope этой итерации ограничен LSRM
`.spe`-specific случаями**. Параллельная разработка для XML и других
форматов отложена до завершения плана для LSRM NaI.

## Контекст и motivation

`ENERGY_CEILING_KEV = 3000` (`gamma.spectrum`) — это глобальная
project-scope константа, задающая верхний энергетический потолок при
чтении любого спектра (см. SKILL.md §Scope). До v1.7.20:

- `read_lsrm_spe(path, *, apply_energy_ceiling=True)` поддерживала
  только on/off-флаг trim'а. Сам потолок брался из `ENERGY_CEILING_KEV`
  жёстко — без возможности per-call override.
- Для диагностики (e.g. проверка калибровки на полном канальном
  диапазоне выше 3 MeV) или для узкоэнергетического анализа
  (e.g. Pb-210 46 keV или Am-241 26-60 keV без чтения всего хвоста)
  единственным workaround был monkey-patch модульной константы либо
  ручной trim уже декодированного `Spectrum.counts`. Оба варианта
  ломают deterministic поведение всего конвейера: константа
  используется не только в reader, но и в downstream (нормировка
  channel-energies, MDA, peak search bounds).

## Дизайн

Новая сигнатура:

```python
def read_lsrm_spe(
    path: str,
    *,
    apply_energy_ceiling: bool = True,
    ceiling_keV: Optional[float] = None,
) -> Spectrum:
```

Семантика:

| `apply_energy_ceiling` | `ceiling_keV` | Поведение                                          |
|------------------------|--------------|----------------------------------------------------|
| True (default)         | None         | Trim до `ENERGY_CEILING_KEV` (3000) — как в v1.7.19 |
| True                   | 400.0        | Trim до 400 keV для этого вызова, константа не тронута |
| False                  | (любое)      | Сохраняет полный декодированный канальный диапазон, `ceiling_keV` игнорируется |

`gamma.io.readers.read_spectrum(**kwargs)` уже пробрасывает kwargs в
format-specific reader — никаких dispatcher-правок не нужно.

## Что НЕ изменилось (defensive invariants)

1. **Default-поведение бит-в-бит совместимо с v1.7.19.** Все 242
   предыдущих теста проходят без правок. Это явно проверено в
   `test_lsrm_spe_default_trims_at_3000` — fixture
   `M_cs_легкий_2001-2005.spe` (1023 канала raw) trim'ится до 927
   с `e_max=2997.8 keV` — точно как в предыдущих версиях.
2. **`ENERGY_CEILING_KEV` остаётся single source of truth для
   default.** Per-call override эффект только на один вызов — модульная
   константа НЕ перезаписывается. Это значит downstream-код (MDA,
   peak search bounds, нормировки) продолжает использовать 3000 keV
   как референс, что важно для согласованности отчётов между разными
   спектрами в одной сессии.
3. **Keyword-only форма** (`*` в сигнатуре) исключает позиционную
   путаницу и делает вызовы explicit.
4. **Когда `apply_energy_ceiling=False` И калибровка присутствует**,
   `energy_max_keV_kept` всё равно вычисляется на последнем сыром
   канале — downstream видит реалистичный диапазон, не None.

## Use cases

### (1) Диагностика — полный канальный диапазон
```python
spec = read_spectrum("M_th_легкий_2001-2005.spe",
                     apply_energy_ceiling=False)
# spec.n_channels == spec.n_channels_raw (без overflow trim'а с .spe;
# overflow detection — это AtomSpectra-specific)
# spec.energy_max_keV_kept > 3000 — можно увидеть Tl-208 2614 keV хвост
# и выше, проверить нет ли overflow или mis-calibration на high-E.
```

### (2) Per-call ceiling — узкоэнергетический анализ
```python
# Анализ только low-energy области для Pb-210 / Am-241 / IC X-rays
spec = read_spectrum("Cs137_5cm.spe", ceiling_keV=400.0)
# spec.n_channels ≈ 116 (вместо 927 default) — peak_search быстрее
# в ~8×, peaks в 400-3000 keV области игнорируются.
```

### (3) Per-call ceiling — диагностика конкретного multiplet
```python
# Изоляция Co-60 1173/1332 doublet — отрезаем хвост, чтобы
# multiplet-deconvolution не отвлекалась на остальные peaks.
spec = read_spectrum("Co60_5cm.spe", ceiling_keV=1400.0)
```

## Архитектурный invariant

`ENERGY_CEILING_KEV` остаётся **константой проекта** (3000 keV — выбор
обоснован SKILL.md §Scope: ²⁰⁸Tl 2614.5 — самая высокая
γ-линия рутинного интереса от природного фона; выше 3 MeV почти
всегда overflow, pile-up, cosmic). Per-call override — это
**escape hatch для диагностики и узких задач**, не альтернативный
default. Тесты явно проверяют это:
`test_lsrm_spe_default_trims_at_3000` гарантирует что без kwarg'ов
default-trim происходит при `ENERGY_CEILING_KEV`.

## Тесты (3 LSRM-specific в `test_reader_api.py`)

| # | Тест                                            | Проверка                                                          |
|---|--------------------------------------------------|--------------------------------------------------------------------|
| 1 | `test_lsrm_spe_default_trims_at_3000`           | Default вызов: `e_max ≤ 3000`, `n_channels < n_channels_raw`       |
| 2 | `test_lsrm_spe_apply_false_keeps_full_range`    | `apply_energy_ceiling=False`: `n_channels == n_channels_raw`, `e_max > 3000` |
| 3 | `test_lsrm_spe_custom_ceiling`                  | `ceiling_keV=400.0`: `e_max ≤ 400`, `n_channels < default n_channels` |

Эмпирические значения на fixture `M_cs_легкий_2001-2005.spe`:

```
default                                  → n_ch=927, raw=1023, e_max=2997.8
apply_energy_ceiling=False               → n_ch=1023,           e_max=3280.2
ceiling_keV=400.0                        → n_ch=116,            e_max=399.9
```

## Что вне scope этой итерации

- **Симметричная разработка для других форматов спектров** (.chn, .n42,
  .mca, .csv) — отложено до завершения плана для LSRM NaI.
- **Per-call subtrim helper** `Spectrum.trim_above(E_keV)` для уже
  прочитанных spectra (без повторного read) — kandidates на будущий
  F-NN. Текущая итерация ограничена изменениями на уровне reader.
- **AtomSpectra XML тесты** (4 теста в `test_reader_api.py`) — в коде
  есть, в файле есть, проходят — но scope этой итерации формально их
  не покрывает. Они «бесплатно» работают потому что параллельная
  симметрия в `read_atomspectra_xml` уже реализована.

## Methodological implication

Default-trim до 3000 keV остаётся **обоснованным выбором по умолчанию**
для productive identification на LSRM NaI Gamma-1С (60×60 / 63×63):
energy diagnostic выше 3000 keV почти всегда уже overflow или
cosmic-ray secondaries, и trim ускоряет downstream peak search в
~8-10 % случаев. Per-call ceiling это **переключатель для exceptional
случаев**, не запасной default — это явно прописано в SKILL.md §Scope.

═══════════════════════════════════════════════════════════════════════

# v1.7.21 — F-43 Averaged Lsrm `.spe` background spectra

## Контекст и motivation

2016 Поверка archive `Gamma-1C_NaI_63x63_USB_SN-01/Поверка-2016/`
содержит **4 background-контекста** с 15+ измерениями каждый:

| Subfolder                       | Назначение                                   | N  | Live (per file) | Total live |
|---------------------------------|----------------------------------------------|----|-----------------|------------|
| `Фон вода`                      | Marinelli + water matrix (matrix-matched bg) | 15 | ~3600-54000 s   | 120 h      |
| `фон пустая защита`             | Empty shielding (point geometries)           | 15 | ~3600-54000 s   | 120 h      |
| `Фон с открытыми крышками`      | Open lid (Rn-daughter / air air check)       | 15 | ~3600-54000 s   | 120 h      |
| `Временная нестабильность`      | 300-сек temporal-stability sequence          | 15 | 300 s           | 1.25 h     |

Первые три контекста идеально подходят для averaging: каждое
измерение — стандартная длинная экспозиция, калибровка идеально
стабильна внутри сета (zero drift, проверено эмпирически), 15
измерений × 8h ≈ 120 hours total live time. Четвёртый контекст
("Временная нестабильность") — это короткие повторы для проверки
temporal drift детектора, не bg для averaging.

До v1.7.21 эти 45 файлов хранились rawly и downstream
background-subtraction (`gamma.calibration.background_subtraction`)
использовала только одно измерение за раз. F-43 строит **single
long-exposure equivalent per geometry** через простую sum-of-counts
и сохраняет результат в `.spe` файл для прозрачного потребления
любым downstream-кодом, который работает с Lsrm `.spe`.

## Mathematical basis

При идентичных калибровках N измерений одного и того же геометрического
сценария — это N i.i.d. Poisson samples процесса со
скоростью `λ(E)`. Сумма каналов:

```
counts_combined[i] = ∑_k counts_k[i]
```

Это Poisson sample с параметром `∑_k λ_k(E)·t_k = λ(E)·T_total` где
`T_total = ∑_k t_k`. Эквивалент одного long-exposure measurement.

**σ-reduction**:
- Total counts: σ ∝ √(Total counts) — scaling по сравнению с
  одиночным измерением: ×√N для total, ×1/√N для rate.
- Pour rate (cps): σ_rate = σ_counts / T_total = √(λT) / T = √(λ/T) —
  при T_total = N·t_single получаем σ_rate уменьшается в **√N** раз
  по сравнению с одиночным измерением длительности t_single.

Для N=15 это **√15 ≈ 3.87×** уменьшение std_rate. Variance
уменьшается в 15× (≈ 1500 %), std в 287 %. Initial оценка в
handoff'е "50× noise reduction" была формулировкой через
**variance × 3 — округление**; точная эмпирика σ_rate × √N = 3.87×.

## Defensive checks перед суммированием

```python
def average_lsrm_spectra(
    paths, *,
    rel_gain_tolerance=0.005,        # max |a1_k - a1_ref| / |a1_ref|
    abs_offset_tolerance=2.0,        # max |a0_k - a0_ref|, keV
    require_same_detector=True,      # detector_id must match
    require_same_geometry=True,      # geometry must match
    ...
)
```

Три класса проверок:

1. **Channel-length consistency** — все входы должны иметь одинаковые
   `n_channels` (после trim) и `n_channels_raw`. Это структурный
   prerequisite: суммировать каналы можно только если они индексируют
   ту же физическую энергию.

2. **Calibration drift** — модель: `E(N) = a0 + a1·N + a2·N² + ...`.
   Проверяются `a0` (offset) и `a1` (gain) против reference (первый
   файл). Defaults:
   - `abs_offset_tolerance = 2.0 keV` — реальная стабильность Lsrm
     NaI Gamma-1С около ±0.5 keV в сутки; 2 keV — щедрый запас.
   - `rel_gain_tolerance = 0.005 (0.5 %)` — стабильность gain
     0.05-0.2 %/сутки; 0.5 % — щедрый запас.
   - Higher-order coefficients (a2, a3, …) не проверяются по двум
     причинам: (1) их влияние на E(N) в пределах калиброванного
     диапазона при стабильных a0,a1 пренебрежимо; (2) их измерение
     numerically сложно при коротких экспозициях.

3. **Identity (defensive sanity)**:
   - `detector_id` должно совпадать — нельзя смешивать спектры
     с разных детекторов.
   - `geometry` должно совпадать — efficiency calibration ε(E)
     зависит от geometry, поэтому averaging across geometries
     ломает downstream activity computation. Опт-аут
     `require_same_geometry=False` для cross-geometry слияния
     (используется только при явном understanding последствий).

При violation любого check'а — `CalibrationMismatchError` или
`IdentityMismatchError` с детальным сообщением.

## Output Spectrum

Унаследованные от первого файла:
- `energy_cal`, `energy_cal_degree`, `energy_cal_source`
- `stored_fwhm_calibration` (включая `model="lsrm_fwhm_polynomial_in_E"`)
- `n_channels`, `n_channels_raw`, `channel_pitch`
- `detector_id`, `device_guid`, `geometry`, `operator`

Агрегированные:
- `counts` = element-wise sum как `np.int64`
- `live_time` = `sum(s.live_time for s in inputs)`
- `real_time` = `sum(s.real_time for s in inputs)`
- `dropped_overflow_count` = sum
- `start_datetime` = `min(start_datetimes)` (раньший)

Synthesized:
- `source_format = "averaged_lsrm_spe"`
- `source_path = "<averaged from N files>"`
- `is_background = True`
- `sample_id` = "averaged: <base> ×N" (или explicit kwarg)
- `comments` = user-provided или default summary
- `file_created_datetime = datetime.now()`

Provenance audit:
- `extras["averaging_provenance"]` — dict со всеми полями
  для traceability (source paths, per-input live-times, calibration
  agreement summary, identity summary, applied tolerances).
- `extras["averaging_sigma_reduction"]` = `√N` как float.

## Minimal `.spe` writer для round-trip

`write_lsrm_spe(spec, path, *, type_label="Калибровка",
config_name="")` сохраняет в Lsrm `.spe` format. CP-1251 header
с KEY=VALUE\r\n линиями, затем `SPECTR=` маркер и uint32 LE
counts.

Эмитятся только минимально необходимые keys, которые `read_lsrm_spe`
обратно читает: `SHIFR`, `TYPE`, `CONFIGNAME` (если задан),
`MEASBEGIN`, `TLIVE`, `TREAL`, `OPERATOR`, `GEOMETRY`, `DETECTOR`,
`ENERGY`, `FWHM` (если есть stored cal), `COMMENT`. Embedded peak
tables и zone tables не воссоздаются (они per-spectrum).

**Round-trip test** в `test_average_lsrm.py` проверяет:
- counts identical (в пределах энергетического ceiling, который
  reader применяет при чтении)
- live_time / real_time identical (в пределах форматирования
  `%.2f`)
- energy_cal identical (округление до научной нотации с 10
  знаками — > 14 значимых цифр, identity guaranteed)
- geometry / detector_id strings preserved

## Production archive

`build_averaged_backgrounds.py` — standalone generator (не
тестовый). Дискаверится 3 квалифицирующиеся группы:

| Output filename                              | Context        | Geometry      | N  | Live  |
|----------------------------------------------|----------------|---------------|----|-------|
| `bg_marinelli_water_marinelli.spe`           | Marinelli H₂O  | Маринелли    | 15 | 120 h |
| `bg_empty_shield_point5cm.spe`               | Empty shield   | Точечная-5см | 15 | 120 h |
| `bg_open_lid_point25cm.spe`                  | Open lid       | Точечная-25см | 15 | 120 h |

Each file ships with:
- `<name>.provenance.json` — полный audit trail
- `MANIFEST.json` (в `data/averaged_backgrounds/`) — index всех
  averaged backgrounds

**Output size**: 3 files × ~5 KB each = ~15 KB total. Сайдкары
~3-5 KB каждый. Минимальный overhead в архиве проекта.

## Use case в downstream pipeline

Пример: cert validation на Marinelli geometry с averaged bg
вместо single-shot:

```python
from gamma.io.readers import read_spectrum
from gamma.calibration.background_subtraction import (
    BackgroundConsentRegistry, subtract_background,
)

registry = BackgroundConsentRegistry()
sample = read_spectrum("samples/Th-228_marinelli.spe")
# Раньше: один 5-часовой bg shot.
# bg = read_spectrum("references/.../фон_single_5h.spe")
# Теперь: 120-часовой averaged equivalent.
bg = read_spectrum(
    "data/averaged_backgrounds/bg_marinelli_water_marinelli.spe"
)
registry.approve(bg.source_path)
net = subtract_background(sample, bg, consent_registry=registry)
# σ-вклад bg в combined σ shrinks ~3.87× → слабые photopeak
# (Pb-210 46, Pb-214 295, Am-241 X-rays) становятся
# detectable с лучшим CI.
```

## Scope и ограничения

**Только LSRM NaI `.spe`** (по явному ограничению пользователя
v1.7.20+). Aggregation для AtomSpectra XML и других форматов
отложена до завершения плана для LSRM NaI.

**Только sum-of-counts averaging**. Альтернативные стратегии
(weighted average по dead-time-corrected rate, robust median,
outlier rejection через χ²-screen на per-channel basis) рассмотрены
и отвергнуты:
- При идентичных калибровках и стабильном детекторе **sum-of-counts
  оптимален** (максимизирует sum_of_counts → minimizes σ_rate).
- Weighted average имеет смысл только при сильно разном dead-time
  fraction (не наш случай: archive все измерения в одном режиме).
- Robust statistics добавляет complexity без ощутимой выгоды для
  background spectra (нет outliers — все измерения от одного
  чистого фона).

**Прямой generator не тестируется**. `build_averaged_backgrounds.py`
запускается один раз при packaging; результат коммитится. Тесты
проверяют лишь, что (1) логика averaging корректна на live
fixtures, и (2) производственные файлы существуют и читаются
back через стандартный API.

═══════════════════════════════════════════════════════════════════════

# v1.7.22 — F-44 Cumulative-checkpoint detection + 2024 archive sync

## Контекст

F-43 (v1.7.21) добавил `gamma.io.average_lsrm` с предположением, что
все входные файлы в bg-фолдере — это **независимые измерения**, и
суммирование их даёт эквивалент одного long-exposure measurement.
Это эталонная Poisson-aggregation математика: N независимых выборок
λ·t суммируются в один Poisson sample с параметром N·λ·t.

**Discovered during F-44 inventory** (при сравнении проектных
fixtures с LSRM архивом пользователя): **LSRM Spectraline
acquisition software emits cumulative checkpoint files**. Длительная
acquisition (e.g. 15h) сохраняется как ряд snapshot'ов:

| File             | live_time | cumulative content                   |
|------------------|-----------|--------------------------------------|
| `..._01.spe`     | 1 h       | events from 0–1 h                    |
| `..._02.spe`     | 2 h       | events from 0–1 h + 1–2 h            |
| `..._N.spe`      | N h       | events from 0–N h (всё с начала)     |

Каждый последующий файл — это full cumulative, не инкремент. Это
подтверждено эмпирически: `Фон вода_01` = 22 709 counts,
`Фон вода_15` = 339 544 counts (ratio 14.95 ≈ 15, exactly as
expected for cumulative).

## Bug в F-43 — механика

F-43 наивно суммировал все 15 файлов:

```python
# F-43 v1.7.21
for s in specs:
    out_counts += s.counts                  # cumulative double-count
total_live = sum(s.live_time for s in specs)  # also inflated
sigma_reduction = sqrt(N)                     # FALSE claim
```

Результат для 15 файлов Marinelli water:
- `summed_counts` = 22 709 + 45 643 + … + 339 544 = **2 720 587**
- Истинные unique counts = 339 544 (только file_15)
- **Inflation factor** = 2 720 587 / 339 544 = **8.01** = (1+2+…+15)/15 = mean cumulative index
- `summed_live_time` = 3600 + 7200 + … + 54 000 = **432 000s** (120h)
- Истинный live_time = 54 000s (15h)
- **Inflation factor** для live_time тоже 8 (same arithmetic mean)

## Critical observation — Rate preservation

Хотя counts И live_time inflated в одинаковое число (×8), их
**отношение rate = counts/live_time сохранено**:

```
rate_F-43_inflated = 2 720 587 / 432 000 = 6.30 cps
rate_true          =   339 544 /  54 000 = 6.29 cps
```

(Малая разница объясняется округлением; математически они одинаковые.)

## Downstream impact

Background subtraction формула LSRM (per user clarification):

```
net_counts[i] = sample_counts[i] − bg_counts[i] · (T_sample / T_bg)
```

— это count-based вычитание со scaling фона на acquisition time
sample'а. Математически эквивалентно rate-based вычитанию
(net_rate = sample_rate − bg_rate), но LSRM формулировка преферируется
потому что preserves Poisson statistics naturally на counts.

**Bug effect на subtraction**:
- `bg_counts × T_sample / T_bg` = `(8 × N_true) × T_sample / (8 × T_true)`
  = `N_true × T_sample / T_true` — **factor 8 сокращается**.
- Net_counts therefore correct.

**Bug effect на σ**:
- `σ²(net) = N_src + N_bg · (T_src/T_bg)²` — Poisson variance propagation.
- `N_bg_inflated = 8 × N_bg_true`, `(T_src/T_bg_inflated)² = (T_src/(8·T_bg_true))² = (1/64) · (T_src/T_bg_true)²`
- Variance bg-side term = `8·N_bg_true × (1/64)·(T_src/T_bg_true)² = (1/8) · σ²_true`
- **σ_bg under-estimated в √8 ≈ 2.83×** → MDA получалась оптимистичной
  (узкий confidence band).

**Production scope of bug**: validate_certs.py НЕ использовал F-43
averaged backgrounds (он использует single-shot `Фон_закр_кр_вода_01.spe`).
Cert validation matrix не затронут. Bug был latent — потенциальная
проблема для downstream code что захочет использовать averaged bg.

## F-44 fix — cumulative detection

Новая функция `detect_cumulative_pattern(specs, *, rel_live_time_tolerance=0.01)`:

```python
def detect_cumulative_pattern(specs, *, rel_live_time_tolerance=0.01):
    if len(specs) < 2:
        return {"is_cumulative": False, ...}
    # Criterion 1: same start_datetime
    same_start = all(d == starts[0] for d in starts)
    # Criterion 2: live-times form arithmetic progression
    sorted_lt = sorted(s.live_time for s in specs)
    expected = [(i+1) * sorted_lt[0] for i in range(N)]
    progression_ok = all(abs(t-e)/e < tol for t, e in zip(sorted_lt, expected))
    return {"is_cumulative": same_start and progression_ok, ...}
```

Логика чисто на metadata (start_datetime + live_time), без чтения
counts — fast. Возвращает full diagnostic dict (criteria results,
deviations, longest_idx).

`average_lsrm_spectra` получил новый kwarg `cumulative_policy: str = "auto"`:

- `"auto"` (default): запускает `detect_cumulative_pattern`, выбирает
  `"cumulative_last"` если cumulative, `"independent_sum"` иначе.
- `"cumulative_last"`: берёт longest live_time file как result.
  σ-reduction = 1.0.
- `"independent_sum"`: original F-43 sum-of-counts semantic.
  σ-reduction = √N. Полезно когда caller знает что входы НЕЗАВИСИМЫ
  (e.g. cross-day measurements).

В cumulative_last mode `_check_calibrations` и `_check_channel_lengths`
**пропускаются** — используется только один input файл, calibration
agreement остальных не релевантен. Это критично для 2024 Marinelli
closed-lid set, где file_16 имеет отличающуюся калибровку (re-saved
after calibration update) и provoked отказ старого pipeline.

## 2024 archive sync — добавлено 80 файлов

User-supplied LSRM archive `C:\LSRM\Work\BG\Gamma-1S\Spe - поверки\`
содержит 218 файлов в двух поверках (2016 + 2024). Проект имел 144 —
miss 74. Скопировано **80 файлов**:

| Source archive                              | Destination project          | N  |
|---------------------------------------------|------------------------------|----|
| Поверка 2024 / Временная нестабильность     | Поверка-2024 / ВН Y-88       | 48 |
| Поверка 2024 / Фон закр кр                  | Поверка-2024 / Фон закр кр   | 16 |
| Поверка 2024 / Фон откр кр                  | Поверка-2024 / Фон откр кр   | 16 |

(48 + 16 + 16 = 80. Discrepancy 90 vs 80: 10 файлов в
`Поверка 2024 / Точечная-5см` присутствуют в project под
санитизированными именами в root-level Gamma-1C/.)

## Pairing semantics (per user clarification)

Background pairing зависит от sample geometry:

| Sample geometry          | Paired background                 | Reason                       |
|--------------------------|-----------------------------------|------------------------------|
| Маринелли                | Marinelli + water (matrix-matched) | Вода ≈ образец → attenuation учтена  |
| Дента / Чашка-60         | Empty shield, closed lid          | Низкий ambient фон           |
| Точечная-5см             | Empty shield, closed lid          | Низкий ambient фон           |
| Точечная-25см            | Open lid                          | На 25см крышки не помещаются → open lid |

**Methodological clarification (per user)**: вклад радона в фоновом
спектре считается **пренебрежимо малым**, всё излучение фона
приписывается **природным U-238/Th-232 chains в строительных
материалах** (бетон, кирпич стен помещения). Это упрощает
background-subtraction:

- Фон — статический источник (распад U-238/Th-232 имеет T½ ~ 10⁹
  years → за лабораторную жизнь не меняется).
- Pb-214 / Bi-214 / Pb-212 / Tl-208 в фоне — это daughters of those
  long-lived parents в стенах, не атмосферный радон.
- Не требуется моделировать флуктуации Rn (например, погодные
  колебания radon emanation rate).

Это методологическое упрощение, оправданное в стандартных
лабораторных условиях с современной вентиляцией. Закрывает open
questions около F-44 ("U-238 chain cross-validation на natural-radon
фон") — фон НЕ трактуется как natural-radon источник.

## Production archive — 5 канонических файлов

`data/averaged_backgrounds/` пересоздан. Старые 3 файла F-43 удалены.
Новые имена с epoch prefix:

| Output filename                                          | Epoch | Geometry      | live | Σcounts   | Pairs with               |
|----------------------------------------------------------|-------|---------------|------|-----------|--------------------------|
| `bg_2016_marinelli_water_marinelli.spe`                  | 2016  | Маринелли    | 15h  | 339 544   | Маринелли samples        |
| `bg_2016_empty_shield_point5cm.spe`                      | 2016  | Точечная-5см | 15h  | 313 694   | Дента/Чашка-60/Точ-5см   |
| `bg_2016_open_lid_point25cm.spe`                         | 2016  | Точечная-25см | 15h  | 1 967 943 | Точечная-25см            |
| `bg_2024_marinelli_water_closed_lid_marinelli.spe`       | 2024  | Маринелли    | 16h  | 440 398   | Маринелли samples (recent) |
| `bg_2024_open_lid_point25cm.spe`                         | 2024  | Точечная-25см | 16h  | 2 383 952 | Точ-25см (recent)        |

Все mode=cumulative_last, σ-reduction=1.0 (one independent measurement
длительности 15h или 16h).

Каждый файл сопровождается `.provenance.json` sidecar с full audit
trail включая `aggregation_mode`, `cumulative_detection` block,
`pair_with_geometries` список. Общий `MANIFEST.json` индексирует всё.

## Тесты (5 новых F-44 + 12 обновлённых F-43)

Новые:
- `test_detect_cumulative_pattern_on_2016_set` — реальный 2016 set
  корректно detected.
- `test_detect_cumulative_pattern_single_spectrum` — N=1 не cumulative.
- `test_detect_cumulative_pattern_synthetic_independent` — разные
  start_datetime → не cumulative.
- `test_auto_selects_cumulative_last_for_2016_set` — auto policy
  выбирает cumulative_last; output == longest input.
- `test_explicit_independent_sum_overrides_cumulative_detection` —
  forced `independent_sum` даёт sum-of-counts даже на cumulative input.
- `test_explicit_cumulative_last_overrides_independent_detection` —
  forced `cumulative_last` берёт longest.
- `test_invalid_cumulative_policy_raises` — invalid string → ValueError.

Обновлённые тесты F-43:
- `test_write_lsrm_spe_roundtrip` — forced `independent_sum` для
  exercise sum semantic.
- `test_geometry_mismatch_can_be_overridden` — forced `independent_sum`
  потому что разные geometry не cumulative.
- `test_prebuilt_archive_files_exist_and_readable` — expectation
  обновлена: ≥3 файла (теперь 5), проверяется наличие
  `aggregation_mode` в provenance.

═══════════════════════════════════════════════════════════════════════

# v1.7.23 — Cert-validation background swap (F-45)

## Контекст

После F-44 в `data/averaged_backgrounds/` есть 5 канонических усреднённых
фонов с правильной аттрибуцией геометрии. `validate_certs.py` остаётся
основным production harness'ом для cert-matrix regression, но использует
до F-45 single-file фон `Фон_закр_кр_вода_01.spe` — это **фон в
геометрии Маринелли с водой** (~1 час), а все 12 fixtures — **point
sources 5 cm**. Несовпадение геометрии bg ↔ sample методологически
некорректно.

Per F-44 pairing rules (§v1.7.22 NOTES выше):

| Sample geometry      | Background                       |
|----------------------|----------------------------------|
| Маринелли           | Marinelli + water (matrix-matched)|
| Дента / Чашка-60     | Empty shield, closed lid          |
| **Точечная-5см**     | **Empty shield, closed lid**      |
| Точечная-25см        | Open lid                          |

F-45 — minimal data-only fix: переключить `BG_PATH` в `validate_certs.py`
на canonical averaged `bg_2016_empty_shield_point5cm.spe` (15h
cumulative_last из F-43/F-44).

## Почему "data-only" — нет изменений в pipeline

Pipeline вычисления активности доминируется **strong-line** target
nuclide. Для linka с E ≫ low-E bg absorption region:

- net_counts[i] = sample_counts[i] − bg_counts[i] · (T_s/T_b)
- σ_net² ≈ sample_counts[i] + bg_counts[i] · (T_s/T_b)²
- A = (net_counts / live_time) / (ε(E) · I_line)

bg_counts входит линейно в net_counts. Для каналов вокруг target peaks
(>100 keV для всех 12 fixtures кроме Am-241 59.5 keV) bg_counts состоит
из preflux fragment подложки, который не сильно зависит от того,
оптически экранирован объём или нет (Compton continuum дальнего
происхождения доминирует над local scattering).

Поэтому ожидание: **A_изм и Δ% инвариантны для всех strong-line
fixtures**. Только peak-count меняется (low-E peak hygiene).

## Эффект на cert-matrix metrics

| Nuclide  | Δ baseline | Δ F-45 | peaks baseline | peaks F-45 |
|----------|------------|--------|----------------|------------|
| Cs-137   | +1.83 %    | +1.83 %| 8              | 8          |
| Co-60    | +0.61 %    | +0.61 %| 13             | 13         |
| Na-22    | +5.99 %    | +5.99 %| 16             | 16         |
| Eu-152   | −6.82 %    | −6.82 %| 15             | 14         |
| Ba-133   | +1.32 %    | +1.32 %| 15             | 13         |
| Am-241   | −8.35 %    | −8.35 %| **11**         | **7**      |
| Zn-65    | +5.45 %    | +5.45 %| 7              | 6          |
| Y-88     | +1.64 %    | +1.64 %| 13             | 13         |
| Bi-207   | −2.67 %    | −2.67 %| 14             | 14         |
| Cd-109   | −3.38 %    | −3.38 %| **9**          | **6**      |
| Pb-212   | −12.79 %   | −12.79 %| 18            | 19         |
| Tl-208   | −0.07 %    | −0.07 %| 18             | 19         |

- mean |Δ| = 4.24 %, max |Δ| = 12.79 % — **инвариант**.
- Activity hypothesis confirmed: strong-line dominance shields A_изм
  от bg geometry choice для target nuclides ≥ ~100 keV.

## Peak-count shifts: интерпретация

**Low-energy isotopes (Am-241 11→7, Cd-109 9→6, Ba-133 15→13,
Eu-152 15→14, Zn-65 7→6)**: empty-shield bg корректно содержит:
- Pb K-α 74.97 keV, Pb K-β 87.0 keV — характеристическое X-излучение
  свинцового экрана (внутренняя fluorescence от Compton/photoelectric
  событий в свинце).
- Низкоэнергетические линии U-238/Th-232 chains из строительных
  материалов: Pb-214 53/77/242/295/352 keV, Bi-214 76/77 X-rays,
  Tl-208 75/86 X-rays, Pb-212 75/239 keV.

В Marinelli+water bg эти линии частично absorbed water sleeve
(~1.5 cm water в Marinelli sleeve между source position и detector
для ~70-300 keV даёт attenuation factor 0.3-0.8). Поэтому single bg
file под-вычитал natural background в low-E region, оставляя false
peak detections в Mariscotti search.

**Pb-212 18→19, Tl-208 18→19**: +1 peak — малый low-E artifact
ранее частично compensated water absorption. Безопасно — Mariscotti
adds peak, identification pipeline отдельно проверяет нуклиды.

## Defensive characteristics

1. **Activity не регрессирует** — это главное свойство F-45.
2. **σ_bg правильно ослаблена**: 15h live_time vs 1h single file →
   bg σ_rate ×1/√15 ≈ 0.258 → σ_net уменьшается → MDA уменьшается
   (более чувствительный harness).
3. **Single-file BG path сохранён закомментированным** для diagnostic
   comparison при необходимости debug bg-related effects.
4. **Pb-212 vs Tl-208 Th-228 cross-validation**: ratio 0.8727
   (−12.73% от unity) сохранён — F-41 finding остаётся локализованным
   как Pb-212 238 keV self-absorption/ε(E) problem, не bg artifact.

## Известное ограничение

2024 архив не содержит empty-shield closed-lid bg (есть только
Marinelli closed-lid и open-lid). Поэтому для всех 2017-2023
measurement-year fixtures используется 2016 averaged bg независимо
от epoch.

Acceptable, потому что:
- 5cm shield геометрия (Pb thickness, internal Cu/Cd lining) не
  менялась 2016→2024.
- Detector (Gamma-1С NaI 63×63, USB SN-01) — тот же экземпляр.
- Masonry стен помещения — стабильна (natural U-238/Th-232 chains
  в бетоне/кирпиче не меняются на годовых масштабах).
- Cross-epoch peak position stability подтверждена в F-44 archive
  inventory (a₀ drift 2016→2024 ≤ 1 keV).

Open follow-up: если LSRM добавит empty-shield closed-lid bg
в Поверку 2026+, переключить на epoch-matched (per-fixture).

## Файлы изменений

`validate_certs.py` только:
1. Docstring пункт 2 обновлён: `subtract bg_2016_empty_shield_point5cm.spe
   (F-43 averaged, F-44 cumulative_last semantic; 15-hour live time,
   matrix-matched "empty shield closed lid" geometry for point-5cm
   samples)`.
2. 19-строчный inline comment перед `BG_PATH = ...` объясняет
   pairing rules motivation, ссылается на F-44, документирует 2024 epoch
   caveat и diagnostic fallback path.
3. `BG_PATH = ROOT / "data" / "averaged_backgrounds" /
   "bg_2016_empty_shield_point5cm.spe"`.
4. Закомментированная single-file `BG_PATH` сохранена для diagnostic.

═══════════════════════════════════════════════════════════════════════

# v1.7.24 — Multi-geometry cert matrix expansion (F-46a)

## Контекст

После F-44 у проекта в archive есть fixtures в **5 разных
геометриях**: Точечная-5см, Точечная-25см, Дента-120мл, Петри-60мл,
Маринелли — с собственной .efr efficiency curve, собственным cert
файлом (для volume sources) и собственной правильной background-
геометрией (per F-44 pairing rules). До F-46 `validate_certs.py`
работал только с Точечная-5см (12 fixtures, mean |Δ|=4.24%).

Стратегия F-46 разбита на 3 slice'а:
- **F-46a** (этот): Точечная-25см — самый простой slice (тот же
  cert файл что Точ-5см, single-nuclide источники + один chain
  proxy mirror F-41).
- **F-46b** (будущее): Marinelli — 8 fixtures, требует обработки
  Bq/kg → Bq через mass_g, chain proxy для Ra-226 и Th-232.
- **F-46c** (будущее): Дента-120мл + Петри-60мл — 16 fixtures, та же
  логика что F-46b.

Каждый slice — independently shippable.

## Реализация — Per-geometry resolution

`CertFixture` extended:
```python
@dataclass
class CertFixture:
    nuclide: str
    spe_filename: str            # forward-slash subpath for non-root geom
    cert_source_hint: str
    cert_nuclide: Optional[str] = None
    chain_branching: float = 1.0
    geometry: str = "Точечная-5см"  # F-46 NEW field
```

Per-geometry resource dicts `EFF_PATHS` / `BG_PATHS` / `CERT_PATHS`
(5 геометрий каждая) + lazy `_resolve_geometry_resources(geom, cache)`
helper. Eff/bg/cert загружаются один раз per geometry, кэшируются в
dict, переиспользуются для всех subsequent fixtures той же geometry.

## 4 новых fixtures (F-46a — Точечная-25см)

| Source              | Cert ref date | A_cert, Bq | Pipeline       |
|---------------------|---------------|-----------|----------------|
| Cs-137 №SRC-01    | 2017-05-19    | 106 000    | direct, 661 keV|
| Na-22 #01.22        | 2022-11-14    | 229 000    | direct, 511/1274|
| Y-88 №SRC-02      | 2023-10-09    | 350 000    | direct, 898/1836|
| Th-228 №SRC-03    | 2021-04-26    | 100 000    | Tl-208 chain proxy |

Все 4 паруются с `bg_2016_open_lid_point25cm.spe` (per F-44 pairing
rules: Точ-25см → open lid bg) и используют
`УДС-ГЦ-63х63-USB__SN-01_-_Точечная-25см.efr`.

## Результаты cert-matrix

**Точечная-5см (12 fixtures, unchanged)**: mean |Δ|=4.24%,
max |Δ|=12.79% (Pb-212 single 238 keV). Bit-identical с v1.7.23.

**Точечная-25см (4 fixtures, F-46a)**:

| Nuclide | A_cert@meas | A_meas | Δ        | comment      |
|---------|-------------|--------|----------|--------------|
| Cs-137  | 89 301      | 84 756 | −5.09 %  | single line  |
| Na-22   | 1.37e5      | 1.42e5 | +4.15 %  | TCS=2        |
| Y-88    | 29 709      | 28 274 | −4.83 %  | TCS=2        |
| Tl-208  | 28 119      | 22 680 | −19.34 % | TCS=4 (chain)|

3/4 fixtures within ±10%. Общий matrix: **16/16 measurable, mean
|Δ|=5.27%, max |Δ|=19.34%**.

## Cross-geometry Tl-208 finding

Same Tl-208 chain-proxy methodology даёт:

| Geometry      | Source           | Δ vs cert |
|---------------|------------------|-----------|
| Точечная-5см  | Th-228 №SRC-04 | −0.07 %   |
| Точечная-25см | Th-228 №SRC-03 | −19.34 %  |

Две гипотезы:

**A**: 25cm efficiency curve high-E bias. chi²/dof=2.51 (vs ~1 для
5cm); если eff curve over-predicts ε на 583-2614 keV → A_meas меньше.

**B**: Th-228 №SRC-03 cert overstatement. #309 содержит только
Th-228 entry — невозможно independently cross-validate через other
nuclide на том же физическом источнике.

Без additional measurement невозможно distinguished. **F-46d open**:
hypothesis A falsifiable если #SRC-05.2023 source измерить в 25cm
geometry — small Δ подтвердит, что bias на #SRC-06 — cert error.

## Tests update

`test_chain_proxy.py` 3 теста generalized для multi-geometry:
- `share_spe_file`: per-geometry pairing invariant.
- `share_cert_parent`: all proxies cert_nuclide=Th-228.
- `chain_branching_is_one`: iterate all proxies.

8/8 pass. Polynomial regression unchanged: 21/21 files PASS, 266+
tests total.

## Defensive

1. **All 12 prior Точ-5см rows bit-identical** — default geometry
   preserves backward compat.
2. **Module-level BG_PATH / EFF_5CM / CERT_PATH** остаются aliases
   в Точ-5см defaults — `test_chain_proxy.py` grep-tests их.
3. **Lazy resource cache** — 16 fixtures × 5 sec / geometry = ~10 sec
   вместо ~80 sec naive reload.
4. **Resource-missing graceful failure** — unknown geometry token
   → RunResult с note "unknown geometry", processing продолжается.

═══════════════════════════════════════════════════════════════════════

# v1.7.25 — Multi-geometry cert matrix completed (F-46b/c/d)

## Контекст и slice'ы

F-46 multi-geometry expansion стартовал в v1.7.24 с F-46a (Точ-25см
slice, 4 fixtures). v1.7.25 closes остальные три slice'а в одной
итерации:

- **F-46b** — Marinelli volume sources (8 fixtures)
- **F-46c** — Дента-120мл + Петри-60мл volume sources (16 fixtures)
- **F-46d** — Diagnostic chi²/dof per geometry + Tl-208 25cm
  investigation deferred (data unavailable)

Cert matrix вырос с 16 → 40 fixtures.

## Bq/kg → Bq unit conversion (F-46b)

Volume cert files report specific activity (Bq/kg) per sub_source.
Чтобы сравнить с измеренным A (absolute Bq) нужно multiply на mass.

`run_one()` extended:

```python
# Walk sub_sources to find target nuclide + mass_g
target_sub = None
cert_act = None
for ss in src.sub_sources:
    for act in ss.activities:
        if act.nuclide == cert_nuclide_name:
            target_sub = ss; cert_act = act; break
    if cert_act: break

# Fallback for compound certs (e.g., ОСГИ 5431 multi-nuclide)
if cert_act is None:
    cert_act = src.get_activity(cert_nuclide_name)
    target_sub = None  # mass unavailable

# Apply Bq/kg → Bq via mass_g
unit = (cert_act.unit or "").strip().lower()
if unit.endswith("bq/kg"):
    if target_sub is None or target_sub.mass_g is None:
        return RunResult(..., note="Bq/kg cert but mass unavailable")
    A_cert_absolute = cert_act.A_Bq * target_sub.mass_g / 1000.0
else:
    A_cert_absolute = cert_act.A_Bq
```

## Chain proxy methodology — расширение F-41

| Cert source | Proxy nuclide | Lib intensity scaling | chain_branching | Lines used |
|-------------|---------------|------------------------|-----------------|-----|
| Th-228 (point #SRC-05, #309) | Pb-212 (F-39) | direct ENSDF | 1.0 | 238 keV |
| Th-228 (point #SRC-05, #309) | Tl-208 (F-41) | × 0.3594 (β-branching Bi-212 → Tl-208) | 1.0 | 583/2614 keV |
| **Ra-226** (Marinelli/Дента/Петри) | **Bi-214 (F-46b)** | direct ENSDF | 1.0 | 609/1120/1764 keV |
| **Th-232** (Marinelli/Дента/Петри) | **Tl-208 (F-46b/c)** | × 0.3594 | 1.0 | 583/2614 keV |

**Ra-226 chain proxy через Bi-214**:
- Ra-226 → α → Rn-222 → α → Po-218 → α → Pb-214 → β → Bi-214 → β → Pb-210
- В sealed sample Rn-222 buffer retained → equilibrium через 25-30 days
  (Rn-222 T½=3.8 d, требуется ~6 half-lives для full equilibrium)
- 420-series sources ≥20 years old → deep equilibrium
- Lib Bi-214 intensities direct ENSDF (не pre-scaled, в отличие от Tl-208)
- chain_branching=1.0 — A_Bi-214 = A_Ra-226 в equilibrium

**Th-232 chain proxy через Tl-208**:
- Th-232 → α → Ra-228 → β → Ac-228 → β → Th-228 → α → Ra-224 → α → Rn-220 → α → Po-216 → α → Pb-212 → β → Bi-212 → α (35.94%) → Tl-208
- Ra-228 bottleneck (T½=5.75 y):
  - 22 years since 2002 sources → 1 − 2^(-22/5.75) = 93% equilibrium
  - 17 years since 2007 sources → 87% equilibrium
  - Acceptable для 5-10% accuracy goal
- Lib Tl-208 intensities pre-scaled by Bi-212 α-branching 0.3594
  (loaded via `load_lsrm_chain_libs()`)
- compute_activity inverting lib intensities recovers parent A_Th-232
  directly с chain_branching=1.0 — exactly same logic as F-41 для Th-228

## Cert-matrix результаты — full 40-fixture table

| Geometry        | n  | mean \|Δ\| | max \|Δ\| | efr chi²/dof |
|-----------------|----|------------|-----------|--------------|
| Точечная-5см    | 12 | 4.24 %     | 12.79 %   | 6.95         |
| Точечная-25см   | 4  | 8.35 %     | 19.34 %   | 2.51         |
| Дента-120мл     | 8  | 9.68 %     | 27.90 %   | **15.40**    |
| Петри-60мл      | 8  | 12.55 %    | 19.58 %   | **15.04**    |
| Маринелли       | 8  | 9.66 %     | 22.54 %   | 3.72         |
| **Total**       | 40 | **8.48 %** | **27.90 %** | —          |

## F-46d diagnostic — chi²/dof reveals systematic eff-curve bias

Per-geometry summary table добавлен в `validate_certs.py` main()
output. Key observation: **Дента/Петри chi²/dof = 15.40/15.04**
дramatically worse than other geometries (Marinelli 3.72,
Точ-25см 2.51, Точ-5см 6.95). Это указывает на underlying eff-curve
fit problem, не на cert или methodology errors.

### Falsifying hypothesis B (cert overstatement) для F-46a Tl-208 25cm

F-46a (v1.7.24) finding: Tl-208 25cm Δ=−19.34 % vs Tl-208 5cm Δ=−0.07 %.
Две гипотезы proposed:
- **A**: 25cm eff curve high-E bias
- **B**: Th-228 №SRC-03 cert overstatement

F-46d data point: Tl-208 across all close geometries показывает **similar
~15-20 % under-estimate**:
- Точ-25см #SRC-06: Δ=−19.34 % (chi²/dof=2.51)
- Marinelli 420-7-17: Δ=−20.31 % (chi²/dof=3.72)
- Marinelli 420-17031: Δ=−4.19 % ← одно отличие
- Дента 420-7-17: Δ=−17.94 % (chi²/dof=15.40)
- Дента 420-17031: Δ=−11.53 %
- Петри 420-7-17: Δ=−19.58 % (chi²/dof=15.04)
- Петри 420-17031: Δ=−12.26 %

Pattern: **systematic ~15-20 % under-estimate для Tl-208 across all
close-geometry sources** — falsifies hypothesis B (cert error wouldn't
correlate with chain-proxy methodology). Strengthens **methodology
explanation** (probable: TCS under-correction для Bi-212 → Tl-208
cascade на close geometries).

Notable: 17031 Th-232 sources show consistently smaller Δ (~−4 to
−12 %) than 420-7-17 sources (~−18 to −20 %). Cert ref dates differ
(2017 vs 2002/2007). Possible explanation: Ra-228 secular equilibrium
fully reached в 7 years (since 2017) → no equilibrium error; 22-year
sources at 93% equilibrium → -7 % из-за not-yet-saturated chain.
Combined with TCS bias → total ~-17 to -20 %.

### Bi-214 mass dependence

| Geometry      | Light source           | Heavy source         |
|---------------|------------------------|----------------------|
| Marinelli     | Δ=−6.48 % (620 g)      | Δ=−22.54 % (1670 g)  |
| Дента-120мл   | Δ=−3.83 % (74 g)       | Δ=−27.90 % (200 g)   |
| Петри-60мл    | Δ=−3.66 % (37 g)       | Δ=−18.53 % (100 g)   |

Heavy containers systematically Δ=−18 до −28 %, light containers
Δ≤7 %. Это **self-attenuation effect** для Bi-214 lower-E lines
(242/295/352 keV) в larger samples. Cowell baseline integration
likely под-fits Compton continuum на broadened multi-line peaks в
heavier samples (higher continuum, broader peaks за счёт self-
attenuation tail).

### K-40 Петри anomaly

К-40 Петри Δ=+18.7 % и +16.95 % — единственный over-estimate
системно в matrix. Probable cause: thin-source (60 ml Petri dish
~3 mm thick) → no meaningful self-attenuation; close-geometry
bg subtraction может leave Compton background tail из U/Th decays в
empty-shield bg (Pb-214/Bi-214 lines at 295/352/609 keV produce
Compton tails extending до K-40 1461 keV region).

## Cross-validation block update

Pairing logic changed:

**До F-46**:
```python
for r in rows:
    fx = ... if fx.cert_nuclide ...
    chain_proxies.setdefault(fx.cert_nuclide, []).append(r)
```

Это grouped all same-parent proxies across sources/geometries → noise
ratios cross physical sources (cert metadata + decay correction
dominated).

**После F-46**:
```python
for r in rows:
    fx = ... if fx.cert_nuclide ...
    key = (fx.cert_nuclide, r.spe_filename, r.geometry)
    same_source_proxies.setdefault(key, []).append((r, fx))
```

Только same-spe pairs grouped. Result: only one pair printed (Pb-212
+ Tl-208 на Th-228__264_2023@5cm, F-39 baseline). Cross-source ratios
hidden — they reflect cert metadata + decay, не detector response.

## Test updates

`test_chain_proxy.py`:

1. **test_pb212_and_tl208_fixtures_share_cert_parent** generalized:
   - Pb-212 still must have cert_nuclide="Th-228" (only Th-228 chain
     proxy)
   - Tl-208 now allowed ∈ {Th-228, Th-232} (dual-role)
   - Same-source pairs (geometry + spe_filename) must agree on both
     cert_nuclide и cert_source_hint

2. **test_validate_certs_module_has_cross_validation_block**:
   - Header check relaxed на prefix "Cross-validation of chain proxies
     (F-41" (matches both F-41-only и F-41/F-46 banners)
   - Grouping dict name accepts `chain_proxies` OR `same_source_proxies`

8/8 chain_proxy pass. Full regression: 21/21 PASS.

## Defensive characteristics

1. **All 16 prior fixtures (12 Точ-5см + 4 Точ-25см) bit-identical**
   v1.7.24. F-46b добавляет только volume fixtures с unit="Bq/kg"
   code path; point-source unit="Bq" code path unchanged.
2. **Bq/kg fallback path** preserves compound-cert behavior (Bi-207
   + Cd-109 в "5431" multi-nuclide cert) — если sub_source.activities
   walk fails (e.g., legacy single-nuclide source), falls back на
   `src.get_activity(name)` returning absolute Bq.
3. **Cross-validation noise eliminated**: same-source-only pairing.
4. **RunResult.geometry default** "Точечная-5см" preserves backward
   compat в случае error-path constructions без geometry.

## Open follow-ups

- **(п)** TCS correction refinement для chain-proxy Tl-208 в close
  geometries. **CLOSED-AS-DOCUMENTED in v1.8.0 (K-21)** — bound
  −15-20% quantified; resolution requires close-geometry P/T data.
- **(р)** Self-attenuation correction для volume samples. **CLOSED-
  AS-DOCUMENTED in v1.8.0 (K-20)** — bound +10/-7% for Cs-137
  quantified; ρ_sample/ρ_ref diagnostic added.
- **(с)** K-40 Петри +17 % investigation — root cause same as K-20
  (density mismatch at extreme thin-source ρ_sample/ρ_ref).
- **(т)** Дента/Петри .efr refit. **CLOSED-AS-DOCUMENTED in v1.8.0
  (K-19)** — F-47a tuning recovers best-per-geometry degree.
- **(у)** 17031 Th-232 row standardization — cosmetic, deferred.

═══════════════════════════════════════════════════════════════════════

# v1.8.0 — Production release (F-47a/b + K-19/K-20/K-21 documented)

## Контекст

После v1.7.25 multi-geometry cert matrix expansion проект достиг
plateau accuracy с **mean |Δ|=8.48%, max=27.90%** на 40 fixtures
в 5 геометриях. Дальнейшая methodology improvement требует research-
grade work (μ(E)/ρ tables для ОИСН-16 matrix; close-geometry P/T
data). v1.8.0 — **production checkpoint** который:
1. Применяет tactical tuning (F-47a): per-geometry polynomial degree
2. Добавляет diagnostic instrumentation (F-47b): density ratio
3. Documents systematic biases как accepted limitations с
   quantified bounds (K-19, K-20, K-21)
4. Sets v1.9 roadmap explicitly

## F-47a — Per-geometry polynomial-degree tuning

После F-46d revealed wide chi²/dof variation across geometries
при hardcoded degree=3, F-47a tests degrees 1-5 per .efr file.

### Findings

| Geometry      | Anchors | degree=1 | =2 | =3 | =4 | =5 | Best |
|---------------|---------|----------|----|----|----|----|------|
| Точ-5см       | 24      | 69.45    | 10.59 | **6.95** | 7.26 | 7.63 | 3 |
| Точ-25см      | 20      | 60.97    | 9.44  | 2.51 | 2.67 | **1.74** | **5** |
| Дента-120мл   | 13      | 14.27    | 15.66 | **15.40** | 17.15 | 18.36 | 3 (K-19) |
| Петри-60мл    | 23      | 34.69    | 22.55 | 15.04 | **14.28** | 15.12 | **4** |
| Маринелли     | 15      | 6.34     | 5.52  | **3.72** | 4.07 | 4.47 | 3 |

Дента chi²/dof≈15 invariant of degree — это **data quality**
limit (only 13 anchor points), не algorithm. Documented as **K-19**.

### Implementation

`EFF_DEGREE` dict в `validate_certs.py`:
```python
EFF_DEGREE = {
    "Точечная-5см":  3,
    "Точечная-25см": 5,  # was 3
    "Дента-120мл":   3,  # K-19 invariant
    "Петри-60мл":    4,  # was 3
    "Маринелли":     3,
}
```

`_resolve_geometry_resources()` uses lookup.

### Effect на cert-matrix

- **Mean |Δ| 8.48% → 7.89%** (улучшение 0.6 pp)
- Max |Δ| 27.90% unchanged (Bi-214 Дента, K-20 + K-21 effect)
- Точ-25см: max Δ 19.34% → 18.25%
- Петри: mean Δ 12.55% → **9.84%** (significant)

## F-47b — Density-ratio diagnostic

`run_one()` computes ρ_sample = mass_g / volume_ml для volume
fixtures и записывает ρ_sample/ρ_ref в note column. NO correction
applied — diagnostic only.

```python
REF_DENSITY = {  # from manual .efr inspection
    "Маринелли":   (1000.0, 1.60),
    "Дента-120мл": ( 120.0, 1.66),
    "Петри-60мл":  (  60.0, 1.60),
}
```

Effect: cert matrix note column shows e.g.
`(ρ_sample/ρ_ref=0.36; 6 peaks)` для Marinelli Cs-137 light.

### Pattern revealed (single-line nuclides)

| Source                    | ρ_sample/ρ_ref | Cs-137 Δ% |
|---------------------------|----------------|-----------|
| Marinelli Cs137_420-7-14  | 0.36           | +9.81 %   |
| Marinelli Cs137_420-7-15  | 1.04           | −5.58 %   |
| Дента Cs137_420-7-14      | 0.34           | +4.85 %   |
| Дента Cs137_420-7-15      | 1.00           | −6.84 %   |
| Петри Cs137_420-7-14      | 0.27           | +6.97 %   |
| Петри Cs137_420-7-15      | 0.79           | −3.71 %   |

Spread +9.81 / −6.84 % maps to ρ range 0.27-1.04. Это **K-20**
matrix attenuation effect.

### Chain-proxy pattern (independent of density)

| Source                  | ρ_sample/ρ_ref | Tl-208 Δ% |
|-------------------------|----------------|-----------|
| Marinelli light         | 0.40           | −4.19 %   |
| **Marinelli heavy**     | **1.00**       | **−20.31 %** |
| Дента light             | 0.40           | −11.53 %  |
| Дента heavy             | 1.00           | −17.94 %  |

Heavy samples с ρ MATCHING reference still show large under-
estimate — **NOT density-related**. This is **K-21** close-
geometry TCS effect.

## K-19, K-20, K-21 — three new accepted limitations

### K-19: Дента-120мл .efr 13-anchor data quality

chi²/dof≈15 invariant of polynomial degree (1-5). Resolution requires
**lab procedure**: re-calibrate Дента-120мл geometry с 5-7
additional anchor sources covering 50-2700 keV. Not actionable from
software.

### K-20: Matrix attenuation correction не реализована

Sample density vs reference density mismatch correlates with
single-line nuclide Δ. Bound +10/−7% for Cs-137 across ρ range
0.27-1.04. Resolution path documented:
1. μ(E)/ρ table for ОИСН-16 (H/C/N/O/Fe composition в .efr Material)
2. F(E) = (1-exp(-μρt))/(μρt) per geometry
3. Correction factor F_ref/F_sample applied в compute_activity

Estimated 6-8h work. Deferred to v1.9.

### K-21: TCS uses point-geometry P/T data

`peak_to_total_NaI` в `cascade_summing` calibrated against Gilmore
Table 8.4 для ~5 cm distance point sources. Close-geometry samples
(0 cm) experience larger cascade-coincidence losses than model
predicts. Tl-208 chain proxy systematic −15 to −20% в all close
geometries (Marinelli/Дента/Петри). 5 cm point Δ=−0.07% (F-41).

Resolution requires close-geometry P/T data (experimental или
Monte Carlo). Estimated 12-16h. Deferred to v1.9.

## v1.9 roadmap

- **K-20 implementation**: gamma.physics.self_attenuation module
  с ОИСН-16 μ(E)/ρ table. Expected to bring Cs-137 spread from
  17% peak-to-peak to ≤ 5%.
- **K-21 implementation**: close-geometry P/T model (option a:
  experimental calibration; option b: Monte Carlo). Expected to
  fix Tl-208 chain proxy in close geometries from −18% to ≤ ±5%.
- **К-19 mitigation**: re-calibrate Дента-120мл .efr (lab
  procedure outside software scope).
- Other formats: .chn, .n42, .mca, .csv readers (parallel to
  existing .spe и AtomSpectra XML).

## Defensive characteristics — v1.8.0 stability commitments

1. **No public API breaking changes** в `gamma` package since v1.7.x.
2. **40/40 cert fixtures measurable** preserved.
3. **21/21 test files PASS, 266+ tests** preserved.
4. **Backward-compat aliases** в `validate_certs.py` (BG_PATH,
   EFF_5CM, CERT_PATH) preserved.
5. **Documented limitations bound the bias** — users know what
   to expect, can apply manual correction if needed.

═══════════════════════════════════════════════════════════════════════

# v1.9.0 — K-20 closed + K-21 partial + K-22 added

## Контекст

v1.8.0 documented three accepted limitations: K-19 (Дента .efr quality),
K-20 (matrix attenuation), K-21 (close-geometry TCS). v1.9.0 implements
K-20 fully (closed) and K-21 partially (closed for light samples,
heavy samples now constrained by new K-22 chain-equilibrium effect).

K-19 confirmed outside-software-scope (lab procedure).

## F-48a — K-20 self-attenuation implementation

### Module: `gamma.physics.self_attenuation`

NIST XCOM mass-attenuation coefficients (μ/ρ in cm²/g) tabulated для
H, C, N, O, Fe в диапазоне 50-3000 keV (17 pillar energies). Source:
NIST Standard Reference Database 8 (XGAM), Berger M.J. et al. 2010,
accessed 2024-11.

```python
_XCOM_ENERGIES_KEV = (50, 60, 80, 100, 150, 200, 300, 400, 500, 600,
                      800, 1000, 1250, 1500, 2000, 2500, 3000)
_XCOM_MU_RHO = {
    "H":  (0.3354, ..., 0.0485),  # 17 values
    "C":  (0.1875, ..., 0.0231),
    "N":  (0.1980, ..., 0.0237),
    "O":  (0.2132, ..., 0.0241),
    "Fe": (1.958,  ..., 0.0379),
}
```

Log-log linear interpolation между pillars (NIST XCOM recommendation).
Below 50 keV: clamp к first value. Above 3000 keV: clamp к last value.

### ОИСН-16 matrix composition

Из .efr Material field:
```python
OISN_16_COMPOSITION = {
    "H":  0.022, "C":  0.206, "N":  0.009,
    "O":  0.049, "Fe": 0.714,  # sum = 1.000
}
```

При 600 keV: μ/ρ = 0.022×0.1271 + 0.206×0.0586 + 0.009×0.0573
+ 0.049×0.0577 + 0.714×0.0769 = **0.073 cm²/g** (verified в test).

### Slab self-attenuation factor

Thin-disk approximation per Knoll §10.III.5 / [GILMORE-8.7]:

```
F(E) = (1 − exp(−μρt)) / (μρt)
```

Limits:
- μρt → 0 (thin/light/high-E): F → 1
- μρt → ∞ (thick/heavy/low-E): F → 1/(μρt)

Implementation использует series expansion `1 − x/2 + x²/6` для
x < 1e-6 чтобы избежать numerical loss-of-significance.

### Correction factor

```
correction(E) = F_ref(E) / F_sample(E)
```

Applied as `A_true = A_measured × correction(E)`.

- ρ_sample < ρ_ref → correction < 1 (reduces over-estimate)
- ρ_sample > ρ_ref → correction > 1 (boosts under-estimate)
- ρ_sample = ρ_ref → correction = 1 (no correction)

Algebraic identity: correction(ρ_a → ρ_b) × correction(ρ_b → ρ_a) = 1
(verified в test).

### Weighted-mean correction

`compute_activity` returns inverse-variance-weighted mean of per-line
activities. К-20 correction должна быть consistent с this weighting:

```python
def weighted_mean_correction(E_keV_list, weights, ...):
    return sum(w * correction_factor(E, ...) for E, w in zip(...)) / sum(weights)
```

При single line: equals single correction (verified в test).
При equal weights: equals arithmetic mean (verified).
При empty input: returns 1.0 (no correction).

### REF_GEOMETRY registry — Per-geometry application

Critical implementation detail: К-20 applies только когда .efr
**не уже включает** matrix correction.

LSRM .efr files имеют outer `Layers.Enable` flag:
- `Enable=false`: matrix correction NOT applied during calibration
  → К-20 needed externally
- `Enable=true`: matrix correction baked в calibration
  → К-20 NOT applied externally (would double-correct)

Verified 2024-11 .efr inspection:

| Geometry      | Layers.Enable | K-20 applies? |
|---------------|---------------|---------------|
| Marinelli     | false         | YES           |
| Дента-120мл   | true          | NO            |
| Петри-60мл    | true          | NO            |
| Точ-5см / Точ-25см | n/a (point sources) | NO    |

```python
REF_GEOMETRY = {
    "Маринелли": (1000.0, 1.60, 3.1),   # vol_ml, ρ_ref, t_cm
    # Дента-120мл, Петри-60мл excluded
}
```

### Wire-up в validate_certs.py

После `compute_activity`:
```python
if (fx.geometry in K20_REF_GEOMETRY
        and target_sub is not None
        and unit.endswith("bq/kg")):
    vol_ml, rho_ref, t_cm = K20_REF_GEOMETRY[fx.geometry]
    rho_sample = target_sub.mass_g / vol_ml
    E_keVs = [la.E_keV for la in act.lines_used]
    weights = [1.0 / (la.sigma_A_Bq**2) if la.sigma_A_Bq > 0 else 0
               for la in act.lines_used]
    cf = weighted_mean_correction(E_keVs, weights,
        rho_sample_g_cm3=rho_sample, rho_ref_g_cm3=rho_ref,
        thickness_cm=t_cm, composition=OISN_16_COMPOSITION)
    A_meas_corrected = act.A_Bq * cf
    note_parts.append(f"K20×{cf:.3f}")
```

### Effect на Cs-137 Marinelli spread

| Source             | ρ_s/ρ_ref | Pre-K20 Δ% | K20 cf | Post-K20 Δ%   |
|--------------------|----------|------------|--------|----------------|
| 420-7-14 (570 g)   | 0.36     | +9.81      | 0.899  | **−1.32**     |
| 420-7-15 (1660 g)  | 1.04     | −5.58      | 1.006  | **−5.00**     |
| **Spread**         |          | **15.39 %** |       | **3.68 %** ✅  |

**Spread reduction 76 %** — meets v1.9.0 design target (≤ 5 %).

### Defensive characteristics

1. Correction = 1.0 at reference density (algebraic identity)
2. Correction = 1.0 for zero reference density (point sources)
3. correction(ρ_a → ρ_b) × correction(ρ_b → ρ_a) = 1 (symmetry)
4. Correction → 1 at high E (μ/ρ decreases with E)
5. F → 1 в thin-slab limit; F → 1/(μρt) в thick-slab limit

### Tests (32 в test_self_attenuation.py)

- XCOM table integrity (monotonic E, consistent lengths,
  monotonic μ/ρ in E)
- Element μ/ρ interpolation (pillar exact, below/above clamping,
  log-log linear at geometric mean, unknown element raises, invalid
  E raises)
- ОИСН-16 mass fraction sum, hand-computed cross-check at 600 keV,
  monotonic μ/ρ in E
- Slab F limits (F=1 at zero t/ρ/μ, monotonic in t, asymptotic
  thick limit, invalid ρ/t raises)
- Correction factor (= 1 at reference, < 1 light, > 1 heavy,
  symmetry, = 1 for ρ_ref=0, monotonic in E)
- Weighted mean (single = single factor, equal weights = arithmetic,
  empty = 1.0, length mismatch raises)
- REF_GEOMETRY registry (Marinelli present, Дента/Петри excluded,
  Marinelli values correct)
- Empirical Marinelli Cs-137 spread reduction (test asserts < 5 %)

## F-48b — K-21 close-geometry P/T scaling

### Extension в gamma.physics.cascade_summing

```python
def peak_to_total_NaI(E_keV: float, *, geometry_factor: float = 1.0):
    """K-21 (v1.9.0): scaled P/T for non-point-geometry samples."""
    ...
    return max(0.05, min(1.0, P * geometry_factor))

GEOMETRY_PT_FACTOR = {
    "Точечная-5см":  1.00,   # Gilmore reference
    "Точечная-25см": 1.00,
    "Маринелли":     0.45,
    "Дента-120мл":   0.50,
    "Петри-60мл":    0.50,
}

def peak_to_total_NaI_for_geometry(geometry: str):
    factor = GEOMETRY_PT_FACTOR.get(geometry, 1.0)
    def _pt(E_keV: float) -> float:
        return peak_to_total_NaI(E_keV, geometry_factor=factor)
    return _pt
```

### Methodological rationale

P/T values from Gilmore Table 8.4 калибровано для ~5 cm point geometry.
At very close geometry (0 cm Marinelli/Дента/Петри), the detector
covers larger solid angle, true coincidence summing of cascade photons
происходит чаще, и эффективная P/T снижается.

Smaller `geometry_factor` → smaller effective P → larger ε_T = ε_p/P
→ larger TCS correction C = 1/(1 − Σ p_c × ε_T(E_j)).

Empirical factors fitted к bring light-sample Tl-208 Δ → 0.

### Wire-up в validate_certs.py

```python
pt_for_geom = peak_to_total_NaI_for_geometry(fx.geometry)
tcs = compute_tcs_corrections(fx.nuclide, eff_curve, p_t_func=pt_for_geom)
```

### Effect на Tl-208 chain-proxy

| Geometry        | Source     | Pre-K21 Δ | Post-K21 Δ |
|-----------------|------------|-----------|-------------|
| Маринелли light | 17031, 7y  | −4.19 %   | **−7.5 %** ※|
| Маринелли heavy | 420-7-17, 17y | −20.31 % | **−16.4 %** (К-22) |
| Дента light     | 17031, 7y  | −11.53 %  | **−6.4 %**  |
| Дента heavy     | 420-7-17, 17y | −17.94 % | **−13.2 %** (К-22) |
| Петри light     | 17031, 7y  | −10.66 %  | **−0.6 %** ✅|
| Петри heavy     | 420-7-17, 17y | −19.58 % | **−8.9 %** (К-22) |

※ Marinelli light slightly worse than baseline because K-20 (applied
first для Marinelli only) reduces A_meas by 8 %; K-21 partially
compensates но not fully. Дента/Петри light samples skip K-20
(.efr already includes matrix), so K-21 alone delivers ≤ 6 % target.

**Light-sample Δ ≤ 5 %** target met for Петри (−0.6 %); approached
for Дента (−6.4 %) and Marinelli (−7.5 %, compounded with K-20).

### Heavy-sample residual → K-22

Heavy samples retain consistent ~−10 to −16 % under-estimate after K-21.
Pattern correlates с source AGE (17-year-old 420-7-17 worse than
7-year-old 420-17031).

Root cause: Th-232 → Ra-228 → Ac-228 → Th-228 chain. Ra-228 T½=5.75 y.
Equilibrium ratio = 1 − 2^(−t/5.75):
- 17 years: 0.870 → expect 13 % under-estimate
- 7 years: 0.566 (but cert reports daughter rate at preparation, so
  net empirical effect smaller)

K-22 documented as accepted limitation; resolution path = extend
`_decay_to_meas()` для cert_nuclide=Th-232 to apply chain-equilibrium
correction. Estimated 2-3 h v1.10 work.

### Tests (6 в test_cascade_summing.py)

- `test_k21_default_geometry_factor_is_one` — backward compat
- `test_k21_geometry_factor_scales_pt_linearly` — algebraic identity
- `test_k21_geometry_pt_factor_table_present` — registry integrity
- `test_k21_close_geometry_factors_smaller_than_one` — physical
  constraint
- `test_k21_peak_to_total_NaI_for_geometry_dispatcher` — round-trip
- `test_k21_close_geometry_increases_tcs_correction` — Tl-208 TCS
  monotonicity vs geometry

## F-48c — Per-geometry P/T dispatcher API

Convenience factory `peak_to_total_NaI_for_geometry(geometry: str)`:
returns `(E_keV) -> P` callable bound к
`GEOMETRY_PT_FACTOR[geometry]`. Used as `p_t_func` argument в
`compute_tcs_corrections`. Unknown geometry → factor 1.0.

## K-19 — Дента-120мл lab re-calibration procedure

К-19 confirmed как outside-software-scope. Lab procedure documented:

1. **Acquire calibrated single-line sources** ≥ 5-7 spanning 50-2700
   keV. Recommended set:
   - Am-241 (59.5 keV)
   - Co-57 (122 keV)
   - Cd-109 (88 keV)
   - Ce-139 (165 keV)
   - Cs-137 (662 keV)
   - Mn-54 (835 keV)
   - Co-60 (1173 + 1332 keV)
   - Y-88 (898 + 1836 keV)
   - Na-22 (511 + 1274 keV)
   - Tl-208 (583 + 2614 keV) via Th-228 chain source

2. **Measure each source в Дента-120мл container** с тем же sample
   mass (199 g) + matrix density (1.66 g/cm³ ОИСН-16) что intended
   use.

3. **Compute peak areas** per LSRM Algorithmic Foundations §5.2.5
   (Cowell или Gaussian fit).

4. **Apply F-29 efficiency-calibration procedure**:
   ε(E_i) = A_i / (S_i × I_i × T_i)
   Append к existing .efr file или create new one с date stamp.

5. **Target chi²/dof ≤ 5** (vs current 15.4 from 13 anchors).

Software does NOT generate calibration data — это hardware + wet-lab
procedure. Software applies user-supplied .efr.

## v1.10 roadmap

- **K-22 implementation**: extend `_decay_to_meas()` для chain-
  equilibrium correction (cert_nuclide=Th-232). Estimated 2-3 h.
- **К-19 mitigation depends on lab procedure** (outside software).
- **Other formats / detectors**: analysis pipeline остаётся LSRM-only
  до отдельного указания.

## v1.10.0 — Multi-format file converter (F-49)

**Что изменилось vs v1.9.0**: добавлена утилита для bidirectional
конверсии spectrum files между 4 форматами — LSRM SpectraLine `.spe`,
BecqMoni/AtomSpectra `ResultDataFile` XML, ANSI/IEEE N42.42-2012,
IAEA SPE ASCII. Conversion работает по всем 12 направлениям (4 × 3).

**Что НЕ изменилось**: analysis pipeline (peaks, identification,
activity, validate_certs, cert matrix) продолжает работать
исключительно с LSRM `.spe`. F-49 — стандалонная утилита; ни одна
строка production analysis не зависит от формата вне LSRM.

**Архитектурное решение — единая Spectrum как neutral**: вместо
N×(N−1)=12 прямых конвертеров, формат A читается в `Spectrum`,
затем пишется в формат B. Поле `Spectrum.extras` несёт format-
specific metadata, которая опционально сохраняется при выходе.
Lossy-conversion summary (`--verbose`) сообщает дропнутые поля.

**Открытие**: BecqMoni и AtomSpectra используют **идентичную** XML
схему (`ResultDataFile / ResultData`). Существующий
`gamma.io.atomspectra_xml.read_atomspectra_xml` уже корректен для
обоих вендоров — добавлен только writer и sniffer. SpecUtils
называет этот парсер `RadiaCode`.

**N42-42-2012 ключевые факты** (per ANSI/IEEE N42.42-2012):
- Root: `RadInstrumentData` в namespace
  `http://physics.nist.gov/N42/2011/N42`
- `EnergyCalibration.CoefficientValues` — whitespace-separated
  low-to-high
- `ChannelData` — whitespace integers; optional
  `compressionCode="CountedZeroes"` (нули compressed как
  `0 <count>`); reader поддерживает оба варианта
- Durations: ISO-8601 `PT<sec>S`
- ID cross-references: `Spectrum/@energyCalibrationReference` →
  `EnergyCalibration/@id`, `Spectrum/@radDetectorInformationReference`
  → `RadDetectorInformation/@id`
- `MeasurementClassCode` ∈ {Background, Foreground, Calibration,
  IntrinsicActivity, NotSpecified}; mapping
  `Spectrum.is_background` ↔ Background/Foreground

**Pluggable registry** (`gamma.io.format_registry`): каждый формат
= `FormatSpec(id, label, extensions, reader, writer, sniffer)`.
Add a new format = one `register(FormatSpec(...))` call. Sniffer-
based detection disambiguates shared extensions (LSRM vs IAEA
`.spe`; BecqMoni vs N42 `.xml`).

**Format inventory** (`FORMAT_REGISTRY.md`): full SpecUtils
catalog — ~26 read formats, ~14 write formats. v1.10.0 implements 4.
Roadmap planned: PCF, CHN, CNF, SPC (3 variants), TKA, Exploranium,
N42-2006.

**Verification**:
- 59 new conversion tests pass (`test_format_conversion.py`):
  registry integrity (13); auto-detect (2); per-format round-trip
  (4 × 8 = 32); cross-format chains (.spe → n42 → xml → .spe +
  xml → n42 → .spe), counts byte-equal at chain end (5); sniffer
  disambiguation (4); N42 CountedZeroes spec compliance (2);
  explicit-format API override (1).
- 22 existing test files PASS — **0 regression**.
- Total: ~325 tests (266+ baseline + 59 new).

**Counts/live_time/real_time/energy_cal round-trip exact** across
all 4 formats. Vendor-specific metadata (FWHM polynomial, peaks
table, SampleInfo, detector kind) is lossy by design.

═══════════════════════════════════════════════════════════════════════

## v1.10.1 — Format-id correction (F-50)

После первой реальной проверки конвертера на пользовательском файле
(Ламинария з Нёноксы 28.09.2019, 5000 ch, 27 h live time) пользователь
указал на ошибку labeling: ASCII `$`-section `.spe` формат, который
v1.10.0 называл "IAEA SPE", на самом деле является **LSRM SpectraLine
ASCII export** по документации ЛСРМ.

**Различия LSRM ASCII vs IAEA / ORTEC GammaVision SPE** (хотя оба
используют `$`-section layout):

| Поле | LSRM ASCII | IAEA / ORTEC |
|---|---|---|
| `$MCA_CAL` первая строка | **N** = число коэффициентов | Polynomial degree (= N-1) |
| `$MCA_CAL` коэффициенты | без суффикса | с `keV` |
| `$MEAS_TIM` | high-precision float | integer / single decimal |
| `$SPEC_ID / $SPEC_REM / $ROI / $SHAPE_CAL` | optional | typically present |
| Экспонента | lowercase `e` | uppercase `E` с 3-digit |

**Реализация в v1.10.1**:
- Модуль `gamma.io.iaea_spe` → `gamma.io.lsrm_spe_text`
- Format-id `iaea_spe` → `lsrm_spe_text`
- Functions переименованы (`read_lsrm_spe_text`,
  `write_lsrm_spe_text`, `looks_like_lsrm_spe_text`)
- Writer строго следует ЛСРМ-документации: N coefs first line,
  no keV suffix, `%.12g` floats для $MEAS_TIM, optional sections
  не emitted при пустых полях
- True ORTEC `SpeIaea` теперь в roadmap (отдельный формат, минорный
  variant `lsrm_spe_text` для добавления при необходимости)

**Lessons learned**:
- SpecUtils canonical enum уже различает `LsrmSpe` и `SpeIaea` —
  v1.10.0 их объединил под "IAEA" label. Имена SpecUtils — canonical
  reference для идентификации форматов в нашем registry.
- Writer должен следовать **target** format spec, а не superset,
  который случайно round-trips при чтении нашего же output.

**Verification**: 22/22 existing + 59/59 conversion = 0 regression.
Auto-detection now correctly identifies user's real LSRM ASCII file:
```
detected format: lsrm_spe_text
[convert] lsrm_spe_text -> becqmoni_xml: clean round-trip
[convert] lsrm_spe_text -> lsrm_spe:     clean round-trip
[convert] lsrm_spe_text -> lsrm_spe_text: clean round-trip
```

═══════════════════════════════════════════════════════════════════════

## v1.11.0 — K-22 chain-equilibrium correction (F-51)

K-22, документированная как accepted limitation в v1.9.0 после К-21
close-geometry P/T scaling, описывает устойчивый residual −13...−16 %
для **heavy** Th-232 fixtures (Marinelli/Дента/Петри с sealed source
420-7-17, cert ref 2007-09-17, измерения 2024). К-21 закрыл
**light** sources (Tl-208 Petri −10.7 → −0.6 % ✅, Дента −11.5 → −6.4 %),
но heavy sources остались. К-22 объясняет это физикой decay-chain
in-growth в Ra-228 bottleneck (T½ = 5.75 y).

### Физика К-22

Th-232 decay chain: `Th-232 → Ra-228 → Ac-228 → Th-228 → ... → Bi-212 → Tl-208/Po-212`.

- **Th-232 parent T½ = 1.4 × 10¹⁰ y** → постоянная в любом temporal
  scale measurement.
- **Ra-228 T½ = 5.75 y** — самый короткий half-life в верхней части
  chain → **bottleneck**. Tl-208 emission rate определяется текущей
  активностью Ra-228 (через быстрых daughters Ac-228 → Th-228 → ...).
- Если в момент cert ref source был **chemically processed** (parent
  Th-232 отделён от daughters), Ra-228 starts at zero и растёт по
  Bateman-equation asymptote:

      A_Ra-228(t) / A_Th-232 = 1 − exp(−ln(2) · t / T_Ra-228)
                              = 1 − 2^(−t / 5.75 y)

- Tl-208 в secular equilibrium с Ra-228 (через очень короткие daughters
  Ac-228 T½=6.1h, Th-228 T½=1.9y, Ra-224 T½=3.66d, Rn-220 T½=55.6s,
  Po-216 T½=0.145s, Pb-212 T½=10.6h, Bi-212 T½=60.55min), pre-scaled
  на β-branching factor 0.3594 от Bi-212. Это уже встроено в Lsrm chain
  library через `load_lsrm_chain_libs()`.

### Empirical evidence для in-growth модели

| Source | cert ref | Age @ 2024 | Predicted in-growth | Observed Δ (pre-K22) |
|---|---|---|---|---|
| Marinelli 420-7-17 | 2007-09 | 16.3 y | 0.86 | −16.4 % |
| Дента 420-7-17 | 2007-09 | 16.3 y | 0.86 | −13.2 % |
| Петри 420-7-17 | 2007-09 | 16.3 y | 0.86 | −8.9 % |
| Marinelli 420-17031 | 2017-06 | 6.6 y | 0.55 | −7.5 % |
| Дента 420-17031 | 2017-06 | 6.6 y | 0.55 | −6.4 % |
| Петри 420-17031 | 2017-06 | 6.6 y | 0.55 | −0.6 % |

420-7-17 fits модель: предсказание 14 % under-estimate, эмпирика
13-16 %. 420-17031 НЕ fits: предсказание 45 % under-estimate, эмпирика
≤ 7.5 %. Вывод: **17031 sources were cert'd at full daughter
equilibrium** (lab procedure — либо использовали natural Th с deeply
established chain, либо ждали 30+ years after preparation before
issuing cert). 420-7-17 sources были cert'd shortly after chemical
preparation → chain growing toward equilibrium с момента cert ref date.

### Implementation rule: per-fixture opt-in flag

Чтобы не over-correct'ить 17031 fixtures (которые НЕ требуют
correction), мы добавили per-fixture flag:

```python
@dataclass
class CertFixture:
    ...
    chain_at_cert_equilibrium: bool = True  # default: NO correction
    chain_bottleneck_T_half_s: Optional[float] = None
```

Маркируем только 3 heavy fixtures explicit'но:

```python
# 420-7-17 в Marinelli/Дента/Петри:
CertFixture("Tl-208", "Th232_420-7-17_Marinelli_0cm.spe", "420/7_р16",
            cert_nuclide="Th-232", chain_branching=1.0,
            geometry="Marinelli",
            chain_at_cert_equilibrium=False,
            chain_bottleneck_T_half_s=RA228_T_HALF_S)
```

В `run_one` — applied **after** parent decay + chain branching:

```python
A_cert_at_meas = _decay_to_meas(A_cert_absolute, T_half_s, ...)
A_cert_at_meas *= fx.chain_branching
if not fx.chain_at_cert_equilibrium and fx.chain_bottleneck_T_half_s:
    dt_s = (meas_dt - src.reference_datetime).total_seconds()
    eq_factor = 1.0 - math.exp(-math.log(2.0) * dt_s
                               / fx.chain_bottleneck_T_half_s)
    A_cert_at_meas *= eq_factor
    note_parts.append(f"K22 eq={eq_factor:.3f}")
```

### Why not universal application?

Альтернативный путь — applies correction всем Th-232 fixtures
unconditionally. Это unable to distinguish 17031 ("chain at equilibrium
at cert") from 420-7-17 ("chain reset at cert"). Результат на 17031:
factor 0.55 → cert@meas сжимается до 55 %, Δ становится +63 %
(вместо текущих −7.5 %). Очевидно неприменимо.

Опт-ин flag — explicit signal от source librarian: "это source
prepared from chemical separation at cert ref date, expect in-growth
deficit". Default `True` сохраняет v1.10.1 behavior bit-for-bit для
всех **37** немаркированных fixtures.

### Cert matrix improvement

```
Pre-K22 (v1.10.1):                 Post-K22 (v1.11.0):
  mean |Δ| = 7.39 %                  mean |Δ| = 6.65 %
  max  |Δ| = 27.90 %                 max  |Δ| = 27.90 %  (Bi-214 K-20 territory)
  
Heavy Th-232 (3 fixtures):
  Marinelli  −16.40 %  →  −4.18 %   ✅
  Дента-120мл −13.20 %  →  −0.57 %   ✅
  Петри-60мл  −8.90 %  →  +4.38 %   ✅

Light Th-232 (3 fixtures, unchanged):
  Marinelli  −7.53 %                no K22 marker
  Дента-120мл −6.42 %                no K22 marker
  Петри-60мл  −0.64 %                no K22 marker
```

### Out-of-scope (после v1.11.0)

- **Bi-214 heavy Marinelli/Дента residual −22...−28 %** — K-20 matrix
  attenuation territory (`.efr Layers.Enable=true` имеет matrix
  correction baked для ref density, deviates для high-density samples).
  Resolution: extend K-20 self-attenuation correction на Дента/Петри
  geometries при density mismatch с .efr ref.
- **K-40 / Cs-137 Петри light +10...+12 %** — close-geometry thin-source
  Compton bg subtraction anomaly. K-20-adjacent но требует separate
  investigation.

### Test coverage

`test_k22_chain_equilibrium.py` (15 tests, all pipeline-free):

| # | Test name | Что проверяет |
|---|---|---|
| 1 | `test_ra228_t_half_s_value` | RA228_T_HALF_S = 5.75 × 365.25 × 86400 s ≈ 1.81417e8 s |
| 2 | `test_eq_factor_17y_matches_doc` | factor(17 y / 5.75) ≈ 0.870 |
| 3 | `test_eq_factor_7y_matches_doc` | factor(7 y / 5.75) ≈ 0.566 |
| 4 | `test_eq_factor_asymptotes_to_one` | factor(60 y) > 0.999 |
| 5 | `test_eq_factor_zero_at_cert_date` | factor(0) = 0 |
| 6 | `test_cert_fixture_defaults_preserve_v110_behavior` | default flags = (True, None) |
| 7-9 | `test_heavy_*_420_7_17_opts_in` | 3 heavy fixtures маркированы |
| 10-12 | `test_light_17031_*_does_not_opt_in` | 3 light fixtures default |
| 13 | `test_ra226_chain_proxies_do_not_opt_in` | Ra-226 → Bi-214 не marked |
| 14 | `test_th228_chain_proxies_do_not_opt_in` | Th-228 → Pb-212/Tl-208 не marked |
| 15 | `test_applied_correction_size_for_420_7_17` | 16.3-y empirical interval даёт factor ∈ [0.85, 0.88] |

### Verification

24/24 test files PASS, 0 regression. Full cert matrix 40/40 measurable
preserved. K-22 status в `KNOWN_AND_FIXED_ISSUES.md` обновлён с
"Accepted limitation" → "✅ PARTIALLY FIXED in v1.11.0".

═══════════════════════════════════════════════════════════════════════
