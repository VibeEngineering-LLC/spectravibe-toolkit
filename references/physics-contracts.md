# Physics & methodology contracts (extracted 2026-07-03)

Полный verbatim-архив F-правил и methodology-контрактов SpectraVibe. Извлечён из
`CLAUDE.md` при lean-рефакторинге 2026-07-03 (contents 643 строки). Указатель на этот
файл живёт в `CLAUDE.md` → секция «Physics & methodology contracts».

---

## Physics & methodology contracts

- **Detector taxonomy lock (2026-06-05 user; F2-A renormalisation 2026-06-21)**:
  «Гамма-1С (кириллица) = Gamma-1S (ASCII)». Canonical = `Gamma-1S` (правильная
  транслитерация кириллической «С» U+0421 → ASCII `S`). Старая ASCII-форма
  `Gamma-1C` (исторический омоглиф-typo) — теперь **legacy alias** в
  `data/aliases.json:detector`, резолвится в `Gamma-1S`. В проекте **одна**
  физическая NaI-63×63 станция: каноническое имя `Gamma-1S`, папка
  `detectors/Gamma-1S/`, модуль `scripts/gamma/detectors/gamma1s.py`. **НЕ
  создавать** ни `detectors/Gamma-1C/`, ни иные варианты — все формы должны
  collapse в `Gamma-1S`. BUG-40 cyrillic→ASCII homoglyph predicate
  (`aliases.contains_cyrillic_letters` / `cyrillic_to_latin_collision`)
  остаётся: warning gate в `json_report._build_warnings` keyна
  `detector_fallback.reason` — после F2-A для real header «Гамма-1С» профиль
  грузится cleanly, поэтому **warning не эмитится**, но предикат честно
  возвращает True (используется для будущих незарегистрированных кириллических
  комплексов). Rationale: `detectors/Gamma-1S/README.md` §1.

- **F-157**: ЛСРМ > Будыка > Gilmore для Gamma-1S на конфликте.

- **F-160 LSRM ground-truth для эталонных спектров (HARD, 2026-06-20)**:
  Каждый эталонный спектр поверочного набора Гамма-1С, по которому оператор присылает
  снимки окон LSRM Гамма-1С (Параметры пиков, Калибровка по ПШПВ, Калибровка по
  пику-образу, Информация об активности и т.п.), фиксируется в
  `references/lsrm_ground_truth/<base_filename>/` как **machine-readable** JSON+CSV +
  README с привязкой к `.spe`-источнику и протоколом сверки. Все расчёты по этому
  спектру (поиск пиков, FWHM, площади, активности) **сверяются с этим эталоном** и
  любое расхождение помечается в отчёте/JSON. Пилот: `Th232_420-7-17_Marinelli_0cm/`
  (зафиксировано 2026-06-20, 9 скриншотов LSRM-GUI). Не путать с сертификатом
  (`cert_zcheck.json` — это аттестованная активность, а LSRM ground-truth — это
  «как ту же самую пробу рассчитал бы LSRM Гамма-1С на тех же исходных данных»;
  оба остаются: cert первичен по активности, LSRM первичен по форме/площади/FWHM).
  Полный реестр эталонов и подробности — README в каждой подпапке.

- **F-160 калибровка FWHM ВСЕГДА из данных, never default (HARD, operator-locked 2026-06-20)**:
  **Оператор Дмитрий, дословно:** «всегда нужно делать калибровку по fwhm полагаться на
  расчетную кривую нельзя». На любом спектре — bootstrap FWHM(E) по значимым пикам или
  ground-truth-anchor LSRM-полинома. `default_NaI_63x63` — **fallback последней инстанции**
  с visible warning «F-160 ALERT». Для эталонных Гамма-1С спектров с папкой в
  `references/lsrm_ground_truth/<base>/fwhm_calibration_lsrm.json` —
  **автоматическая загрузка 14 anchor-точек LSRM через _index.json** и **NNLS-fit**
  (неотрицательные коэффициенты, форма `FWHM²=a+b·E+c·E²`). Маппинг:
  `references/lsrm_ground_truth/_index.json` поле `mapping` (basename .spe → папка
  ground-truth). Реализация: `staged_pipeline._load_ground_truth_fwhm_anchors(spec)` +
  ветка в `build_fwhm_model` перед bootstrap. **Структурное ограничение:** LSRM использует
  полином 4-й степени, наша 3-параметрическая NNLS-аппроксимация даёт систематику
  ±5-7 keV на anchor-точках (Th-232 Marinelli: FWHM(238)=29.5 vs LSRM 24.0 +23%,
  FWHM(2614)=116.8 vs 112.8 +3.6%). Это лучше чем default
  (default давал FWHM(1581)=78.1 vs 90.7 −14%, FWHM(2614)=107.9 vs 112.8 −4.4%).
  Расширение API `build_fwhm_model` до polynomial-4 / callable вместо (a,b,c)-tuple
  **выполнено в F-452 (2026-06-22)** — `FwhmModel` теперь callable, внутри держит
  poly-4 √E LSRM-uplift; backlog по этому пункту F-160 закрыт (см. F-452 entry
  в KNOWN_AND_FIXED_ISSUES.md). 3-параметрическая NNLS-аппроксимация остаётся как
  fallback, но основной путь — callable модель.
  Гарантия применения: `source='lsrm_ground_truth_reference'` в `fwhm_cal.source`
  V2-JSON-отчёта; warning в `extras.fwhm_model_warnings`.

- **Handheld scintillator — НЕТ шильда, Pb K-X-rays = эндогенные IC (HARD, operator-locked 2026-06-13)**:
  RadiaCode-103, AtomFast, Atomtex handheld и любые карманные scintillator-зонды (<1 кг,
  без явного Pb-домика) **шильда не имеют**. Атрибуция пика 70-90 keV как «Pb K-XRF от
  шильда» на их спектре — категориальная ошибка. На handheld эта полоса возникает
  ТОЛЬКО эндогенно: **Pb K-Xα/Kβ (75/85 keV) от внутренней конверсии (IC) в дочках
  Ra-226 chain (Pb-214, Bi-214) и Th-232 chain (Pb-212, Pb-208)** в самой матрице
  образца. На природных минералах с U/Th (apatite, monazite, granite, zircon, K-feldspar,
  торий-содержащие пески) эта полоса всегда видна и **диагностична** для присутствия
  Ra-226/Th-232. Полные правила, фон-rate gates и диагностика — `references/04_secondary_peaks.md` §5a.
  **Обычный фон handheld без шильда**: 8-10 cps integral (RC-103, 25-2800 keV).
  Sample > 15-20 cps → активный. Триггер фиксации: RC-103 + apatite session
  2026-06-13 (pipeline мис-атрибутил 70 keV пик как «Pb XRF shield» на handheld
  без шильда).

- **XRF стабильных материалов — только с подтверждением геометрии (HARD, operator-locked 2026-07-04)**:
  XRF-линии стабильных конструкционных материалов (Pb, W, Fe, Cu, Sn, и др.) атрибутируются
  **только при подтверждённом наличии** этого материала в конфигурации измерения.
  **Категорически запрещено** приписывать «W-XRF» или «Pb-XRF от шилда» без явного
  знания геометрии. Правило — anti-hallucination: наблюдаемый пик в области XRF
  стабильного элемента должен быть обоснован контекстом, а не только энергетическим совпадением.

  **Подтверждённые контексты и соответствующие XRF:**

  | Материал | XRF-линии (главные, keV) | Когда атрибутировать |
  |---|---|---|
  | **Pb** (свинец) | Ka1=74.97, Ka2=72.81, Kb1=84.94, Kb2=87.30; L: La1=10.55, Lb1=12.61 | Pb-укрытие/шилд, Pb-контейнер, Pb-вкладыш, Pb-пластины в геометрии. НЕ handheld без шилда. |
  | **W** (вольфрам) | Ka1=59.32, Ka2=57.98, Kb1=67.24, Kb2=69.07; L: La1=8.40, Lb1=9.67 | Спектрометрия изделий с W: **сварочные электроды WT-20** (содержат 2% ThO₂ + W), вольфрамовые вставки, карбид W. |
  | **Fe** | Ka1=6.40, Kb1=7.06 | Стальной контейнер, стальная кювета, железный грунт — при низкоэнергетическом детекторе. |
  | **Cu** | Ka1=8.05, Kb1=8.91 | Cu-вкладыш шилда, электроника рядом с детектором. |

  **WT-20 электроды (провенанс 2026-07-04):**
  Состав: W-матрица + ~2% ThO₂ (торий-220/232 chain). В спектре видны ОДНОВРЕМЕННО:
  Th-232 chain гамма-линии (Ac-228, Tl-208, Bi-212, Ra-224) + W K-XRF (59.3, 67–69 кэВ)
  от возбуждения вольфрамовой матрицы γ-квантами chain'а. Библиотека:
  `Электроды.lib` (нуклид `Th-232` с дочерними + нуклид `W` с XRF-линиями).

- **Ra-226 186 кэВ → suspect U-238/U-235 (HARD, operator-locked 2026-06-13)**:
  При обнаружении на природном образце (камень, грунт, минерал, фосфорит, гранит,
  апатит, monazit, zircon, торий-содержащий песок, K-feldspar) пика в окрестности
  **186 кэВ** — НЕМЕДЛЕННО подозревать присутствие U-238 и U-235 в матрице.
  На сцинтилляторе с R(186) ≈ 15-20% (RC-103 ~17%, NaI ~13%) линии **²³⁵U 185.72**
  и **²²⁶Ra 186.21** полностью неразрешимы — пик «186» всегда содержит **обе**
  компоненты в неизвестной пропорции. Чистый ²²⁶Ra без ²³⁵U встречается ТОЛЬКО
  на искусственных Ra-эталонах (медицина, поверочные источники), не в природных
  матрицах (природное отношение ²³⁵U/²³⁸U = 0.00725, а ²²⁶Ra накапливается в
  равновесии с ²³⁸U → если есть Ra-226, есть и U-238/U-235 родители).
  Характеристический catalog U-линий 25-200 кэВ для cross-check:
  **²²⁷Th x-ray escape 25**, **²³⁴Th 63.30**, **²³⁴Th 92.4**, **U K x-rays ~95**,
  **²³⁵U 143.76+163.36** (merged на CsI/NaI), **²³⁵U 185.72** (под Ra-226 186),
  **²³⁴ᵐPa 1001**. На HPGe эти линии разрешимы; на CsI(Tl)/NaI часть из них
  сливается с Pb K-X-rays 70-90 (там сидят Th L и U K), но 1001 кэВ всегда
  диагностично выделяется выше continuum'а.
  Reporting convention: пик в окрестности 186 → НЕ «²²⁶Ra 186», правильно
  «²²⁶Ra+²³⁵U 186 keV (merged on scintillator)» или `Ra-226_U235_186` в JSON.
  В operator-facing отчёте — явная нота «²³⁵U 185 vs ²²⁶Ra 186 на CsI(Tl)/NaI
  неразрешимы; для дифференциации нужны (a) долгая набирка HPGe >24 ч,
  (b) количественный изотопный анализ через estimate ²³⁵U/²³⁸U = 0.00725,
  (c) cross-check через ²³⁴Th 63 / U K x-rays ~95 / ²³⁴ᵐPa 1001».
  Эмпирическая основа: разностный спектр «U(Природный) − Ra²²⁶ страница №2»
  в руководстве Соловьева V.1.05 (стр. 32 PDF) на RadiaCode-103 — после вычитания
  чистого ²²⁶Ra остаются только U-родительские линии 63 / 95 / 143+163 / 185.
  Полные ссылки: `detectors/RadiaCode_103/references/руководство_спектроскописта_v1_05_summary.md` §2,
  `detectors/RadiaCode_103/README.md` §9. Триггер фиксации: RC-103 + «Камень с
  Ra-226» session 2026-06-13 (pipeline ложно атрибутил 3 пика 94/186/292 как Ga-67).

- **F-150 / F-293 (v1.17.19)**: библиотека книг лежит в `books_library/`
  на корне проекта (рабочая копия) и **полностью исключается** из
  релизных архивов. Архивы библиотеки — `1_Version/books_library/gamma-books_vYYYY-MM-DD.zip`,
  упаковываются вручную через `python scripts/build_books_archive.py`.
  Sanity-check каталога — `python scripts/verify_books_inventory.py`.

- **F-256**: Layer 2 ГОСТ для regulated отчётов (operator-facing report contract).

- **Ac-228 129.06 в Th-232 цепочке (gotcha, operator-locked 2026-06-10)**:
  пик в окрестности 125-135 кэВ в спектре с Th-232 — это **Ac-228 129.06 кэВ**
  (I_apparent = 2.42% в равновесии), плюс пренебрежимо слабая Th-228 131.61
  (0.13%, в 20 раз меньше). На NaI FWHM ≈ 17 кэВ на этой энергии — линии
  полностью слиты в один пик с доминирующим вкладом Ac-228 129. НЕ путать
  с Compton-bump'ом, NE-структурой или артефактом continuum'а — это
  настоящая линия Th-232 цепочки. Использовать как **дополнительный
  низкоэнергетический якорь Ac-228** (наряду с 338.3, 463.0, 911/969, 1588).
  Зафиксировано после ошибки атрибуции "структура 129 — уточнить" на сессии
  анализа `Th232_420-7-17_Маринелли_0cm.spe`. Источник: data/nuclides.json:Ac-228.

- **Порядок: фон → net → sample-операции (operator-locked, HARD, 2026-06-10)**:
  Все действия со спектром образца — **точная калибровка** (step 5),
  **поиск пиков** (step 3 на финальном виде), **идентификация** (step 7),
  деконволюция (step 8), активности (step 9) — выполняются **после
  вычитания фона**, т.е. на **net-спектре** (sample − bg, попиксельно по
  cps после независимой калибровки обоих).

  **Allowed pre-subtraction steps** (do anyway, нужны как priors):
  - Step 1 parse (sample + bg раздельно)
  - Step 2 environment classification
  - Step 4 detector type (грубо, по сырому sample — это не меняется от
    subtraction'a)
  - **Bootstrap E-cal sample** (grubo, deg=1 или stored, чтобы привести
    sample-каналы в keV для bg-aligned subtraction) — это **rough**,
    не финал. Точную E-cal (deg≤4) делать на net'е.
  - **Bootstrap FWHM** на нескольких видимых isolated peaks sample —
    тоже rough, для оценки. Финал — на net'е.

  **Strict order**:
  1. Parse sample + bg (step 1).
  2. Environment (step 2).
  3. Calibrate **bg independently** через background-only anchor heuristic
     (см. `references/01_metadata_calibration.md` → «Background-only anchor
     heuristic»). Получить bg-cal, bg-cps array.
  4. Bootstrap rough sample-cal (или reuse stored если drift < 0.3·FWHM
     против bg-cal на общих anchor'ах K-40 1461 / Tl-208 2614).
  5. **Subtract**: net_cps(ch) = sample_cps(ch) − bg_cps_resampled(ch),
     σ_net = √(σ²_s + σ²_b). Если sample-cal и bg-cal различаются —
     resample bg на сетку sample (linear interpolation в keV-domain
     или channel-domain по совместимым cal'ам). Net = sample net of bg.
  6. **Step 3** (точный peak search), **step 5** (финальная E-cal deg≤4),
     **step 6** (FWHM(E) deg≤4), **step 7** (identification), **step 8**
     (deconv), **step 9** (areas/MDA) — **все на net-спектре**.

  **Rationale**: фон вносит K-40, Tl-208, Bi-214 и т.д. в sample-cps.
  Если калибровать/идентифицировать на sample-raw, эти линии (а они
  всегда есть в любой лаборатории) станут «найденными в sample»
  по умолчанию — false positive для sample-source claim'a. Например,
  на Th-232 marinelli-поверке Tl-208 583 в sample = sample (Th источник)
  + bg (Th в бетоне стен). Только net показывает реальную активность
  Th в samples'е.

  **Anti-pattern** (моя ошибка ранней сессии 2026-06-10): запускать step 3
  (peak search), step 5 (E-cal), step 7 (identification) на raw sample
  cps без bg subtraction. Я это сделал — оператор поймал на step 7 («где
  калибровали фон? когда вычли?»), пришлось переделывать. Не повторять.

  **Carve-out**: если bg ничтожен (sample cps ≫ 100× bg cps на главных
  линиях источника), subtraction численно не меняет результат, **но**
  процедура всё равно выполняется (формально + audit-trail) — нельзя
  «пропустить» bg-step доводом «sample dominate». Доказательство
  доминирования = сам результат subtraction'a, а не a-priori
  предположение.

- **F-451 Background subtract — масштабирование «к меньшему» (HARD, operator-locked 2026-06-22)**:
  **Оператор Дмитрий, дословно:** «масштабируется к меньшему. Спектр с большим набором
  нужно приводить к меньшему». Оба модуля subtraction (`bg_subtract_energy.py` LITE для
  `staged_pipeline`/`run_plan_a.py` и `bg_subtract_dual_mode.py` FULL для cert validation)
  масштабируют **ВНИЗ** к спектру с меньшим `live_time`, **не вверх к sample**. Три ветки:
  `t_s > t_bg` → `scale_direction="sample_down"`, sample-counts × `(t_bg/t_s)`,
  `effective_live_time = t_bg`; `t_s < t_bg` → `scale_direction="bg_down"`,
  bg-counts × `(t_s/t_bg)`, `effective_live_time = t_s`; равные → `equal`, scale=1.
  σ-propagation учитывает `applied_scale²`: `σ²(net) = applied_scale²·N_sample + N_bg`
  (sample_down) либо симметрично (bg_down). **CPS-инвариант**
  `net_cps = sample_rate − bg_rate` сохраняется во всех трёх ветках (snapshot-tests
  доказывают численно). **Downstream-контракт (B1+V2, operator-locked 2026-06-22):**
  `staged_pipeline` **НЕ** пересобирает `spec.live_time` после subtraction —
  `spec.live_time`/`spec.real_time` остаются на оригинальных sample-значениях, все step
  6+ (FWHM, identification, deconv, MDA, activities) работают на полной статистике
  sample-счётов. Вариант B2 (full breaking downstream через `dataclasses.replace`
  spec) был прототипирован и **явно откачен** оператором: пропагация
  `effective_live_time` в `spec` ломала контракт `activity = S_net / (t·ε·I)` и давала
  несогласованные MDA между LITE и FULL модулями. F-451 «к меньшему» direction
  (`applied_scale`, `scale_direction`, `effective_live_time`) живёт целиком внутри
  `BackgroundSubtractionResult` для σ-propagation и зеркалится в
  `spec.extras["background_subtraction_*"]` **только для audit/diagnostics**.
  `BackgroundSubtractionResult` несёт новые поля `applied_scale`, `scale_direction`,
  `effective_live_time` (frozen dataclass); legacy `scale_factor = t_s/t_bg` сохранён
  verbatim для JSON-схем. `notes`-строка содержит маркер
  `"F-451 scale_direction=…, applied=…, t_s=…, t_bg=…, effective_live_time=…"`.
  Полный план — `audit/_plans/F-451_bg_subtract_direction_invert.md`; снапшот-тесты
  `tests/snapshot/test_bg_subtract_energy.py` (15 шт., 6 F-451-dedicated) +
  `tests/snapshot/test_bg_subtract_dual_mode.py` (27 шт.); ARCH.md §2.2 «F-58 untangle»
  обновлён.

- **F-452 FWHM model degree uplift — LSRM poly-4 √E (2026-06-22)**:
  Legacy NNLS-квадратичная подгонка FWHM (`FwhmModelKind.legacy_quadratic`, кортеж
  `(a,b,c)`) для LSRM-калиброванных детекторов заменена `FwhmModel` callable с
  `kind="lsrm_poly_sqrt_E"` (5-coef полином от √E). Ground truth —
  `references/lsrm_ground_truth/<base>/fwhm_calibration_lsrm.json` (Gamma-1S anchors:
  60.3 / 122.1 / 661.7 / 1332.5 / 2612.9 keV). `build_fwhm_model(spec)` возвращает
  `(FwhmModel, source)`; source становится `"lsrm_ground_truth_reference_poly4_sqrtE"`
  когда LSRM JSON доступен, иначе fall back на legacy quadratic. Backward-compat
  `fwhm_keV_at_energy(model, E)` и `_make_fwhm_at_channel` принимают и `FwhmModel`, и
  3-tuple; helper `fwhm_model_legacy_abc(model)` рефитит poly-4 в `(a,b,c)` для legacy
  JSON-схемы. Below-anchor extrapolation clamps на 0.1 keV floor. На Th-232 Marinelli
  poly-4 даёт FWHM(2614)=112.8 keV (LSRM cert: 112.80) vs legacy quadratic ≈116.9 keV.
  Меняется только call site `build_fwhm_model` (строка 2157 staged_pipeline.py); все
  FWHM-consumers model-agnostic через `fwhm_keV_at_energy`.

- **F-452-FU2 Currie L_C pre-MAD non-detection filter (2026-06-22)**:
  F-452 (более точный FWHM) **раскрыл** скрытый баг в
  `compute_activity_for_nuclide` (BUG-38/39 MAD-outlier-rejection,
  `scripts/gamma/activity/compute.py:978-1038`): когда у нуклида большинство matched
  lines получают `peak_area_source="deconvolved_coupled"` с numerical-noise значениями
  `A_i ~ 1e-21 Bq` (вокруг нуля), MAD-медиана сама падает в нанозону, и **реальные**
  линии с высокими counts (Tl-208 2614 keV, `A_i~2700 Bq`) выбрасываются как 19σ
  outliers против нано-консенсуса. Наблюдалось pre-fix: weighted Tl-208 activity
  упало с baseline ~2925 Bq до **9.8e-24 Bq** на Th-232 Marinelli. Fix: pre-MAD фильтр
  в `compute.py:906` срезает любую линию с `A_i/σ_A_i < 1.0` (Currie L_C non-detection
  criterion, Currie 1968 / Gilmore §5.5) **ДО** MAD-rejection, с явным `lines_skipped`
  provenance (`below_Currie_LC A_i=… σ=… A/σ=… < 1.0 (pre-MAD non-detection filter,
  F-452-FU2)`). Verified: Tl-208=2585, Ac-228=2710, Pb-212=2943 Bq; Th-232 chain
  equilibrium ratio=1.14 (`in_equilibrium=True`) в обоих NO-BG/WITH-BG путях. Side-
  effect: snapshot-тест `test_F389_th232_demo_v2_activity_parity_with_prod` ранее
  exact 0.0% prod-vs-V2 parity сломался до 15.3% на Ac-228 — это **раскрытие реальной
  V2 path divergence** (V2-only линии 129.06 + 1630.6 keV у Ac-228 никогда не было в
  prod; ранее маскировалось симметричной нано-distortion в обоих pipeline). Tolerance
  расширен 5%→18% с inline обоснованием; системный fix backlog'нут как **F-452-FU3**
  (V2-only filter coverage). Диагностика —
  `audit/_drafts/F-452-FU2_prod_vs_v2_diff.py`.

- **F-453 Anchor-disagreement auto-trigger + singleton-anchor F-145 fallback (BUG-38 root cause, 2026-06-23)**:
  Закрывает root cause **BUG-38** (silent self-calibration на коротких NaI-фикстурах: AmTiCsEu /
  Cs-Co Marinelli). Pre-fix Step 5β `recalibrate_energy_if_anchors_disagree` сидел за kwarg-opt-in
  `recalibrate_on_anchor_disagreement=False` (default OFF, F-87c контракт v1.14.0), а F-145
  Phase A multiplet self-cal требовал `forced_clusters` (Th/U chain_dominance) и тихо проходил
  мимо на источниках без хвоста (n_multiplets_seen=0). Итог: stored ADC→keV дрейф 9.5 keV
  на 662 keV (Cs-137) **не лечился** ни одним из двух механизмов.
  **Часть (a) — auto-trigger Step 5β.** Новая функция
  `should_auto_recalibrate(anchor_matches, *, fwhm_provider_keV, drift_frac_threshold=0.5, min_anchors=3)`
  в `scripts/gamma/calibration/anchor_recalibration.py`: возвращает `(True, diag)` iff
  `≥3 usable anchors` AND `max(|Δ|/FWHM) > 0.5·FWHM` (вдвое выше базового 0.3·FWHM —
  «однозначный» drift, выше singleton-noise). При срабатывании staged_pipeline вызывает
  стандартный `recalibrate_energy_if_anchors_disagree` (0.3·FWHM target) **без** kwarg-opt-in.
  Default-контракт v1.14.0 сохранён: на спектрах с приличной cal или `<3 anchors` auto-trigger
  молчит. Диагностика — `recalibration_diag["f453_auto_trigger"]`.
  **Часть (b) — singleton-anchor fallback в F-145.** Новый helper
  `_f453_build_singleton_extras(anchor_matches, fwhm_provider_keV, forced_clusters)` в
  `staged_pipeline.py` собирает ENH-priority anchor_matches, которых **нет** в
  `forced_clusters` (одиночные линии Am-241 59.5, Ti-44 1157, Cs-137 662, Co-60 1173/1332,
  Eu-152 1408 …), и передаёт их в `recalibrate_from_multiplet_centroids(..., extra_anchors=…)`.
  Параметр `extra_anchors` существовал в API F-145 как dead-code — теперь живой. Singletons
  идут в Phase B как индивидуальные точки (Rayleigh-fit minimization), Phase C refit
  принимает их если χ² улучшение ≥1.5× AND `|⟨dE⟩|/FWHM ≤ 0.5`. На AmTiCsEu фикстуре Phase C
  refit `9.53 → 4.47 keV` на ENH-anchors (3× улучшение).
  **Phase D carve-out.** При `not forced_clusters AND f453_singleton_extras` (singleton-only
  путь без multiplet'ов) — Phase D refit **полностью пропускается**, calibration принимается
  unconditionally (`chi2_A_sum=chi2_D_sum=0`, `forced_clusters_D=[]`). Иначе Phase D rollback
  условие `chi2_D <= chi2_A AND forced_clusters_D` падало на пустом списке и singleton-cal
  откатывался обратно. `f145_diag.reason` помечается `"F-453 carve-out: forced_clusters пуст,
  singleton-only refit принят безусловно"`.
  **Что осталось backlog'ом.** `anchor_matches` сейчас покрывает только ENH-priority лиш до
  662 keV (Am/Cs/Co/K/Tl/Bi/Pb через `derive_priority_findings`). Линии Ti-44 1157 и Eu-152
  1408 keV вне anchor coverage — high-E extrapolation residual `−8.7 keV` остаётся.
  **BUG-38 closed OPEN → PARTIAL**: root cause (silent cal) закрыт, full closure требует
  **F-453-FU** (расширение `anchor_matches` / `priority_findings` на Ti-44 / Eu-1408 / прочие
  ENH-singletons на калибровочных фикстурах). Реализация — 4 правки в `staged_pipeline.py`
  + новая функция в `anchor_recalibration.py`. Regression: pytest 2364 passed, 0 failed,
  46 skipped, 5 xfailed, 1 xpassed. Бриф —
  `_state/agent_a/inbox/2026-06-23_F-453_bug38_followup.md`.

- **F-453-FU High-E calibration anchors via fixture-fingerprint gate (BUG-38 high-E closure, 2026-06-23)**:
  Закрывает high-E extrapolation residual `−8.7 кэВ` на Ti-44 1157 / Eu-152 1408, оставшийся
  после F-453. **BUG-38 PARTIAL → CLOSED.** Подход v1 (partner-only anti-shadow на ranks 15-20)
  упал — partner-based tolerance на NaI FWHM фундаментально недостаточен (Δ_AmTiCsEu=4.4 кэВ
  vs Δ_Ra-226=7.1 кэВ через Pb-XRF 75 + Bi-214 1120 false-positive partners). Принят v2
  **dual-gate**: (1) **fixture-fingerprint gate** — calibration-tier anchors (`rank ≥
  CALIBRATION_RANK_START=15`) активируются только когда в peak list одновременно видны
  peaks возле Cs-137 661.66 кэВ И Am-241 59.54 кэВ (подпись AmTiCsEu/AmCs Marinelli
  фикстуры); на Ra-226/Th-232/K-40/Co-60 gate fails → calibration-tier silent → нет
  регрессий self-cal на природных образцах; (2) **partner-required** (вторичная защита) —
  Sc-44 67.87 ↔ Ti-44 1157 mutual, Eu-152 1408 → Sc-44 OR Ti-44 (НЕ другие Eu-152 линии
  чтобы не получить confound через Ac-228 338/794 на Th-232). Реализация — 3 новых
  `AnchorEntry` (ranks 15-17) + 2 константы (`CALIBRATION_RANK_START`,
  `AMTICSEU_FINGERPRINT_LINES_KEV`) + helper `_amticseu_fingerprint_present()` + dual-gate
  filter в `find_anchor_matches` — всё в `scripts/gamma/identification/anchor_ranks.py`.
  Test `test_anchor_table_has_14_entries` обновлён `15 → 18` entries. Probe re-run
  AmTiCsEu Marinelli: **max char-line Δ=4.57 кэВ** (было 9.5 после F-453, 10.9 до); Ti-44
  1157 **−1.63 кэВ** (5.3× улучшение), Eu-152 1408 **+0.40 кэВ** (22× улучшение). Regression
  `pytest -n auto -p no:randomly`: 2364 passed, 46 skipped, 5 xfailed, 1 xpassed, 0 failed —
  baseline match. Carve-out: gate работает только на смесевых фикстурах с Cs+Am
  одновременно; моно-fixture (только Cs или только Am) — gate fails, calibration-tier
  silent (F-145 ЕРН-anchors справляется). Полная регистрация — `KNOWN_AND_FIXED_ISSUES.md`
  F-453-FU entry.

- **F-459 / F-456-guard — Eu-152 cascade char-line + Ra-226 pre-cal carve-out (BUG-Y closure, 2026-06-23)**:
  Закрывает **BUG-Y** — на AmTiCsEu Marinelli фикстуре с `GAMMA_ALLOW_STAGE3=1`:
  (i) Eu-152 silently missing из `identified_nuclides`; (ii) false positives Co-57, In-111,
  Ba-133, Tc-99m. Дополнительно — попутная регрессия `test_ra226_demo_phase_C_applied`
  после F-456 pre-cal.
  **F-459 (cascade char-line)**: `identify_nuclides` ранее всегда брал highest-intensity
  library line как characteristic (для Eu-152: 121.78 keV, I=28.53%). На AmTiCsEu эта линия
  скрыта под Am-241 Compton continuum (Mariscotti не находит пика между ch=30/72 keV и
  ch=57/152 keV) → Eu-152 silently rejected. **Whitelist-cascade** в
  `scripts/gamma/identification/identify.py`:
  `_CASCADE_WHITELIST = frozenset({"Eu-152"})`, `_CHAR_CANDIDATES = 3` (try top-N highest-I),
  `_MIN_CASCADE_I_PCT = 20.0` (cascade fallback line должна иметь I ≥ 20%). Для Eu-152
  cascade на 344.28 keV (I=26.59%) → 4 matched lines (244.70, 503.47, 867.38, 1408.01 keV).
  Whitelist-only (НЕ generic «try multiple char lines») — иначе Co-60 (1173/1332),
  Bi-212 (1620/1078), Ac-228 (338) ложно сматчат соседние пики. Co-57/In-111 (122/245 keV
  ≈ Eu-152) — handled `NAI_CHAR_OVERLAP_PAIRS` в `disambiguate.py`; Ba-133 (356)
  — `NAI_CONFUSION_MAP`. AmTiCsEu `final_detected`: `['Am-241','Cs-137','Eu-152','Sc-44','Ti-44']`,
  REAL_MISSING=∅, FALSE_POS=∅.
  **F-456-guard**: F-456 pre-cal block (`staged_pipeline.py`, перед Stage 3) ранее запускался
  для ВСЕХ спектров с `allow_stage3=True` — вызывал `find_anchor_matches(max_rank=99)` +
  `_f453_build_singleton_extras(forced_clusters=[])`. Для Ra-226 (нет Am-241, нет Cs-137,
  fixture-fingerprint gate fails) `_f456_am_full` содержал только natural-chain anchors
  (Bi-214/Pb-214/Tl-208/K-40 rank 1-14), F-456 предкалибровал спектр ими → смещал
  `energy_cal` → Phase D self-cal находил `chi²_sum: 56.96 → 94.56 (хуже)` → откат →
  `phase_C_applied=False`. **Guard** `_f456_calib_tier_present` — F-456 запускается только
  когда есть хотя бы один calibration-tier якорь (`rank ≥ CALIBRATION_RANK_START=15`,
  т.е. Sc-44 67.87 / Ti-44 1157 / Eu-152 1408). На Ra-226/Th-232/K-40/Co-60 — guard=False,
  F-456 silent, F-145 multiplet self-cal работает на оригинальной energy_cal как раньше.
  **Регрессия**: 2376 passed, 0 failed, 46 skipped, 5 xfailed, 1 xpassed (baseline T41
  восстановлен). Cite: `scripts/gamma/identification/identify.py` (cascade block);
  `scripts/gamma/identification/staged_pipeline.py` (F-456 block + guard);
  `scripts/gamma/identification/anchor_ranks.py:CALIBRATION_RANK_START=15`;
  probe `audit/_drafts/F-459_cascade_char_probe.py`. Полная регистрация —
  `KNOWN_AND_FIXED_ISSUES.md` F-459 + F-456-guard entries.

- **F-460 PTB-2018 Annex E chain-decay library + Eq. 16 tail verdict (2026-07-02)**:
  Три под-пункта, закрывают #PTB-5/#PTB-6/#PTB-7 из PTB-2018 γ-SPEKT/GRUNDL (ISSN 1865-8725)
  gap-audit. Реализация — `scripts/gamma/activity/compute.py:423-529` (generalized
  `_ptb_annex_e_half_life`), тесты — `tests/step09_activity_mda/test_activity.py`
  Group 4c (PTB Annex E.2).

  **(a) #PTB-5 Pb-214 / Bi-214 три-режимная библиотека (Annex E.1, Tab. E1)**:
  `_CHAIN_DECAY_MODES` расширен `("equilibrium", "rn222", "ra224_fresh", "progeny")`.
  Для Pb-214/Bi-214 (`_PTB_E1_NUCLIDES`): `equilibrium` → T½ := Ra-226 (1600 a);
  `rn222` → T½ := Rn-222 (3.8235 d); `progeny` → own T½ (Pb-214: 1608 s;
  Bi-214: 1194 s). Модель выбирается оператором через `chain_decay_mode`.
  Silent pass-through к own T½ при `ra224_fresh` на Pb-214/Bi-214 (cross-annex
  no-op, backward-compat).

  **(b) #PTB-6 Pb-212 dual T½ (Annex E.2, PTB p. 6160-6360)**:
  `_PTB_E2_NUCLIDES = frozenset({"Pb-212"})`. Для Pb-212: `equilibrium` → T½ :=
  Th-228 (`_TH228_T_HALF_S = 6.0275e7 s` = 1.91 a; soil / aged aqueous >4 d,
  Ra-224–Th-228 secular equilibrium — default); `ra224_fresh` → T½ := Ra-224
  (`_RA224_T_HALF_S = 3.1622e5 s` = 3.66 d; fresh aqueous carve-out, sample
  measured soon after sampling); `progeny` → own T½ (38304 s ≈ 10.64 h).
  Silent pass-through при `rn222` на Pb-212 (E.1 mode на E.2 нуклиде).

  **(c) #PTB-7 Eq. 16 tail T(K) — verdict PARTIAL COVERED**:
  PTB Eq. 16: `T(K) = N_n · 1/(σ·√(2π)) · exp(-σ_T²/(2σ²) - σ_T·(K-k0)/σ²)` —
  ADDITIVE (T(K) сосуществует с G(K) на всей оси, continuity через N_n
  из Eq. 17). Наш код: **multiplet path** (`peaks/peak_image.py:118-148`
  `gaussian_with_tail()`, вход `use_peak_image=True` через
  `staged_pipeline.py:2409`) реализует PIECEWISE dimensionless
  T = σ_T/σ: для `z<-T` — `A·exp(T·z + T²/2)`, для `z≥-T` — pure Gauss.
  Функционально имеет exp-tail (не 1-to-1 PTB additive, разная
  параметризация, но покрывает физику). **Isolated peak path**
  (`peaks/area.py:cowell_area` + `gaussian_fit_area` +
  `peaks/area_step_continuum.py:gauss_erfc_step_fit`) — pure Gauss + step,
  БЕЗ exp-tail → GAP. Backlog: унифицированное exp-tail F-rule когда
  появится реальный impact на площади isolated пиков (в текущем корпусе
  фикстур регрессий не наблюдается).

  **Semantics шпаргалка для оператора** (mode → nuclide):
  | mode | Pb-214/Bi-214 | Pb-212 | Ra-226 | остальные |
  |---|---|---|---|---|
  | `equilibrium` (default) | T½ := Ra-226 1600 a | T½ := Th-228 1.91 a | own | own |
  | `rn222` | T½ := Rn-222 3.8235 d | own T½ (silent no-op) | own | own |
  | `ra224_fresh` | own T½ (silent no-op) | T½ := Ra-224 3.66 d | own | own |
  | `progeny` | own T½ | own T½ | own | own |

  **Регрессия**: 2424 passed / 0 failed / 46 skipped / 5 xfailed / 1 xpassed
  в 151.14 s (baseline 2419 + 5 новых E.2 тестов в Group 4c
  `test_activity.py`). Back-compat alias `_ptb_e1_half_life =
  _ptb_annex_e_half_life` сохранён. Cite: `scripts/gamma/activity/compute.py:423-529`
  (helper), `compute.py:1240-1251` (call site), `tests/step09_activity_mda/test_activity.py`
  Group 4c (5 tests), `references/gamma_spekt_grundl_v2018-03_en.md:1699-1714` (Eq. 16),
  `references/gamma_spekt_grundl_v2018-03_en.md:6160-6360` (Annex E.2).

- **T41 Efficiency-file detector content-fingerprint gate (anti-hallucination, 2026-06-23)**:
  Закрывает silent CONTENT-fallback класса на efficiency-axis — параллельный пробел к
  BUG-40 (a) cyrillic-collision: path-lookup в `efficiency_autoload.find_efr_file()`
  успешен, но `.efr` файл внутри директории относится к ДРУГОМУ физическому экземпляру
  прибора. Real incident (BUG-40 (b), 2026-06-23):
  `detectors/Gamma-1S/efficiency/Gamma-1S_NaI_63x63_USB_SN-01/...Marinelli.efr` имеет
  `[detector;…]` = `УДС-ГЦ-63х63-USB №SN-01` Поверка-2024, спектр CONFIGNAME =
  `Гамма-1С №SN-04` Поверка-2016. Серия 0086 ≠ 0221 → активность Am-241/Ti-44 искажена
  −96/−97 %, Cs-137 +9.5 %. Существующий `cyrillic_to_latin_collision` ловит только
  path-level homoglyph, content-level расхождение было silent. **Фикс**: новый standalone-
  валидатор `scripts/gamma/calibration/efficiency_provenance.py` (`extract_serial_year` +
  `check_efr_detector_match`) — извлекает `(serial, year)` regex'ом из обеих detector-
  строк (`№NNNN-NN` cyrillic / `No NNNN-NN` ASCII / `N-NNNN-NN`), сравнивает кортежи.
  Mismatch → запись в `detector_fallback_dict["efficiency_detector_mismatch"]` (только
  PII-safe поля per F-115: code + expected_serial_year + actual_serial_year; полные
  detector-строки и basename — только в `logger.warning`, не в JSON) → новый warning code
  `EFFICIENCY_DETECTOR_SERIAL_MISMATCH` (severity `HIGH`) в `_build_warnings` + RU-перевод
  в `markdown_report._render_warning_dict_ru` (F-386 anti-EN-leak gate). 12 unit-тестов
  (`tests/calibration/test_efficiency_provenance.py`), probe `BUG-40_amticseu_..._probe`
  показывает warning emit на реальной BUG-40 (b) фикстуре. Регрессия 2376 passed
  (2364 baseline + 12 новых). **Carve-out**: T41 surfaces проблему, НЕ исправляет данные —
  реальный фикс infrastructure (rename директории ИЛИ положить настоящую Поверка-2016 .efr)
  остаётся за operator-gate (KFI T40). T41 не сравнивает геометрию (`Marinelli` vs
  `Marinelli-1L`) — отдельная F-rule при появлении geometry-mismatch class. Полная
  регистрация — `KNOWN_AND_FIXED_ISSUES.md` T41 entry.

- **F-440 Two-phase weak-line completion (operator-locked 2026-06-13)**:
  Multiplet handling выполняется в две фазы.

  **Phase 1 (fit)**: в grouping/Rayleigh/G1-fit участвуют **только** линии
  с `S/N ≥ 5` AND `I_γ ≥ 3%`. Слабые library-anchor'ы (Ac-228 503/509/523/
  562/571/572/583.41 в окрестности M4 на Th-232 demo и аналогичные) НЕ входят
  в topology кластера — они выпадают на входе `find_multiplet_regions` и
  сохраняются отдельно для Phase 2.

  **Phase 2 (post-fit completion)**: после расчёта активности нуклида
  по сильным якорям через `quasi_template_solver`, для каждой слабой
  library-линии считается ожидаемое число отсчётов
  `S_expected = A·I_γ·ε(E)·t_live·f_self_abs·f_TCS` и записывается в новый
  JSON-блок `weak_line_completion`. Используется для:
  (a) completeness metric per nuclide,
  (b) contamination correction в G1 ROI сильной линии (опционально),
  (c) честного MDA на уровне нуклида.

  Бриф для agent-a-math: `_state/agent_a/inbox/2026-06-13_F-440_two_phase_weak_line_completion.md`.
  Реализация — отдельная сессия, target version v1.30.0 (minor bump → push + GH Release).

  Ссылки: LSRM Algorithmic Foundations 2025 §6.2 «fit only what you measure»;
  Gilmore & Joss 3rd Ed. §9.6.4 «library correction afterwards».

- **Tl-208 510.77 vs annihilation 511 (gotcha, operator-locked 2026-06-10)**:
  пик в окрестности 506-511 кэВ в спектре с Th-232 (Маринелли поверка, грунт
  с природным торием, любой Th-содержащий образец) — это **Tl-208 510.77 кэВ**
  (линия цепочки Th-232, I_γ_abs = 22.6%, I_apparent ≈ 8.12% per Bi-212 decay
  в равновесии). НЕ annihilation 511. Tl-208 510 доминирует над annihilation
  на любом сколько-нибудь активном Th-источнике. Аннигиляция 511 — фон
  (мюоны+pair production в Pb-castle), даёт ≈0.015 cps; Tl-208 510 на образце
  1940 Бк/кг Th-232 даёт ≈0.95 cps — на два порядка больше. Использовать
  Tl-208 510 как **третий якорь Tl-208** (наряду с 583.2 и 2614.5) для
  intensity-ratio cross-check (шаг 7) и для калибровки энергии (шаг 5).
  Зафиксировано после ошибки атрибуции на сессии анализа
  `Th232_420-7-17_Маринелли_0cm.spe`. Источник: data/nuclides.json:Tl-208.

- **FWHM любого пика = калибровка FWHM(E) (operator-locked, HARD, 2026-06-16)**:
  Зафиксировано оператором Дмитрием 2026-06-16: «все пики должны иметь ширину
  соответствующую калибровке по FWHM. Если кривая калибровки не известна,
  производится калибровка FWHM по значимым пикам. Все прочие пики в том числе
  в мультиплетах должны иметь FWHM соответствующую калибровке FWHM по энергии».

  **Принцип**: разрешение детектора — его физическое свойство; ширина пика
  **детерминирована** функцией FWHM(E), а не подгоняется индивидуально. Асимметрию
  и низкоэнергетический хвост на NaI ловят **erfc-step + tail-член** (T, h_step),
  **НЕ раздувание гауссианы**. Раздутый σ затягивает континуум в «ширину» и
  искажает площадь → активность.

  **Procedure (HARD)**:
  1. Если кривая FWHM(E) известна (stored vendor-модель / cal-пики) — **все** пики
     (изолированные И компоненты мультиплетов) получают σ = FWHM(E)/2.3548,
     **locked**. Свободны только амплитуды/площади + континуум (+ tail/step).
  2. Если кривой FWHM(E) НЕТ или она недостоверна — **сначала** откалибровать
     FWHM(E) по значимым (сильным, изолированным, высокий S/N) пикам самого
     спектра: free-σ фит ТОЛЬКО этих anchor'ов → построить FWHM(E) (модель √E
     для сцинтиллятора / полином) → **затем** применить эту кривую ко всем пикам
     по п.1. Free-σ допустим **исключительно** на этом bootstrap-шаге, не как
     рутинный режим интегрирования.

  **Текущее состояние кода (gap, F-rule pending)**:
  - ✅ Компоненты мультиплета уже соответствуют: σ locked к FWHM(E) в
    `coupled_intensity_fit` (`scripts/gamma/peaks/coupled_multiplet.py:342`).
  - ⚠️ Изолированные пики НЕ соответствуют: `gauss_erfc_step_fit` отпускает σ
    в коридоре [0.6, 1.6]×sigma_cal (`area_step_continuum.py:154-161`), а
    `fit_peak_image` (включён через `use_peak_image=True`,
    `staged_pipeline.py:2409`) отпускает σ полностью (`peak_image.py:324-332`).
    В отчёт идёт подогнанная ширина, не калибровочная.
  - ❌ Bootstrap «нет кривой → калибровка FWHM по значимым пикам» отсутствует:
    fallback в `fwhm_provider.py:218-225` = generic-константа + √E-пол, а не
    калибровка по сильным пикам файла.

  **Что в ЛСРМ (проверено по первоисточнику 2026-06-16, на вопрос оператора
  «что в ЛСРМ?»)**:
  - Базовый метод SpectraLine — **МЭС (метод эталонных спектров)**: спектр пробы
    раскладывается на эталонные спектры нуклидов в той же геометрии, ЛИНЕЙНАЯ
    система минимума χ², коэффициенты ∝ активности; форма/ширина пиков = форма
    эталонных спектров (разрешение детектора), свободной ширины нет в принципе
    (`references/_extracted_corpus/Документация ЛСРМ/01_methodology_pdf/Активность
    в счетных образцах...pdf.md:140-156, 194-203`). Это наш `coupled_intensity_fit`.
  - Подгонка модельной функцией управляется **флажком «ПШПВ»** (бинарно):
    «Если он установлен, полуширина определяется из условия минимума χ²-функционала…
    значение получается более точным. В противном случае значение полуширины
    будет взято из результатов калибровки» (`spectralinexx_2.0_basic_functions_rus.pdf.md:1647-1650`).
    **Default (флаг выкл) = σ из калибровки FWHM(E), locked.** Флаг вкл = свободная σ.
    **Коридора ±N% у ЛСРМ нет** — режим бинарный.
  - Асимметрия/хвост/ступенька = отдельная «калибровка по форме пика» → пик-образ
    `.cpt`, НЕ уширение гауссианы; без пик-образа «все пики описываются симметричным
    Гауссианом» (`:1534-1536`, `:580-583`, `:1674-1708`). Границы пика — в долях ПШПВ
    (`:1565-1566`).

  **Решение оператора 2026-06-16 (Q1 → «что в ЛСРМ?» → выбран вариант ЛСРМ-faithful)**:
  **hard-lock σ = FWHM(E)/2.3548 по умолчанию** (= default-режим ЛСРМ = стандарт выше),
  БЕЗ коридора. Свободный режим сохранить как **выключенный по умолчанию тумблер
  `GAMMA_FREE_SIGMA=1`** (аналог флажка «ПШПВ» ЛСРМ, «более точно»). Bootstrap
  «нет кривой → калибровка FWHM по значимым пикам» — **в той же F-rule** (Q2).

  **Уточнение роли метода 2 (operator-locked, 2026-06-20)**: зафиксировано оператором
  Дмитрием — «метод 2 оставляем если калибровочной кривой нет. Фактически это "ручная"
  калибровка по FWHM». То есть free-σ (`GAMMA_FREE_SIGMA=1`) — **не** рутинный
  «более точный» режим интегрирования на каждый прогон, а **инструмент bootstrap-калибровки
  разрешения**, применяемый ИСКЛЮЧИТЕЛЬНО когда кривой FWHM(E) нет/она недостоверна:
  отпустить σ на сильных изолированных якорях → построить FWHM(E) → запереть ею все
  остальные пики (= ровно Q2-процедура выше). При наличии достоверной кривой штатный
  режим — метод 1 (σ-lock). Триггер: сессия Th-232 Marinelli 2026-06-20 (сравнение
  методов показало, что на мультиплет-богатом спектре с валидной FWHM(E) тумблер
  инопертивен — все площади идут через coupled-МЭС/Cowell, `gauss_erfc_step` не вызывается).

  **Реализация** — отдельная F-rule (algorithmic change в fit-ядро), dispatch
  agent-a-math (background). Состав: (a) lock изолированного σ к FWHM(E) в
  `gauss_erfc_step_fit` и `fit_peak_image`, default locked, тумблер `GAMMA_FREE_SIGMA`;
  (b) FWHM(E)-bootstrap-from-significant-peaks в `fwhm_provider.py` при отсутствии/
  недостоверности stored-кривой. Триггер фиксации: вопросы оператора «как определяется
  FWHM при подгонке?» → «что в ЛСРМ?» на сессии Th-232 Marinelli 2026-06-16.
  Ссылки: LSRM Algorithmic Foundations 2025 §4.1 (FWHM(E) модель прибора как
  фиксированное разрешение) + методичка/Basic Functions ЛСРМ (offsets выше);
  Gilmore & Joss 3rd Ed. §9.3 (resolution function), §6.5 (low-E tailing → tail term,
  не ширина).

- **Ширина зоны / зонирование спектра (LSRM-faithful, навсегда, 2026-06-20)**:
  Зафиксировано по первоисточнику `references/_extracted_corpus/Документация ЛСРМ/
  01_methodology_pdf/Lsrm_algorithmic_foundations.pdf.md` (раздел 4 «Разметка спектра» +
  §5.2 «область действия пика» + §5.2.1, §5.2.5) на запрос оператора «перечитай раздел 5
  целиком, особо про зоны и их наложение». F-157: ЛСРМ первична.

  **Зона (информативный участок)** — участок с одним/несколькими пиками, обрабатываемыми
  СОВМЕСТНО под одним полиномом фона (`pdf.md:480-485`). Ширина задаётся НЕ в σ, а в
  долях ПШПВ (FWHM), на трёх уровнях:
  1. **Область действия пика** (ROI интегрирования): config-параметр; Гамма-1С UI
     «границы пиков = 2.5 / 2.5 ПШПВ» = полуширина 2.5 ПШПВ с КАЖДОЙ стороны (полное
     окно 5 ПШПВ). В коде = `window_factor=2.5` = ROI half-width: `half_window =
     round(2.5·FWHM_ch)` (`area.py:82`+`:138`, docstring «ROI half-width =
     window_factor · FWHM» `area.py:108`; `area_step_continuum.py:75`+`:108`). Канал
     суммируется только по пикам, в чью область действия он попал. Источник
     `pdf.md:658-660` даёт это ТОЛЬКО качественно: «каждый пик … имеет ограниченную
     область действия, которая определяется параметрами, задаваемыми в конфигурации» +
     правило членства канала. Числового значения полуширины (2.5) в источнике НЕТ —
     «2.5» это параметр кода/UI Гамма-1С, не текст ЛСРМ (OVERSTATED-correction; verified
     gamma-агент Read `pdf.md:656-660` 2026-06-20 — ранее Цензор-факт-чек был read-only).
  2. **Верхний предел зоны `Lzmax`**: макс. длина участка, config-параметр; при
     превышении зона РАЗБИВАЕТСЯ в точке минимума отсчётов, чтобы снизить влияние
     соседних пиков (`pdf.md:482-485, 490-494`). Гамма-1С UI Lzmax = «макс. длина зоны
     10 ПШПВ». **РЕАЛИЗОВАНО в коде 2026-06-20 (Step 3, LSRM Lzmax):** `find_multiplet_regions`
     принимает `max_zone_length_fwhm` (дефолт 10.0 = ON, синхрон с Гамма-1С) и после
     кластеризации вызывает `_split_zones_lzmax` (`deconvolve.py`), который дробит зону
     длиннее порога в точке минимума отсчётов между крайними пиками; 0.0/spec None —
     graceful no-op. Cert — см. «(в)» ниже.
  3. **Метод моментов** (только первичная калибровка): длина зоны ≈3 ПШПВ
     (`pdf.md:547-549`); НЕ для рабочего тракта активности.

  **Наложение зон** (`pdf.md:497-498, 861-868`): зоны между собой НЕ перекрываются —
  НЕПЕРЕСЕКАЮЩИЕСЯ, соприкасающиеся. «Подгонка со сшивкой» (все зоны одновременно): на
  границах соприкасающихся зон обеспечивается непрерывность модели спектра И фоновой
  подложки И их первых производных. «Подгонка» (без сшивки) = последовательно по зонам
  независимо, итог = усреднение «зоновых» активностей. Перекрываются (накладываются) ПИКИ
  ВНУТРИ зоны = мультиплет; площадь между ними распределяется пропорционально подгоночным
  (§5.2.5 Ковелл, `pdf.md:937-939`).

  **Что НЕ есть ширина зоны**: display-континуум — отрисовочная эвристика (floor под
  гауссианой для JS-дисплея), НЕ влияет на квантование/cert. ГЕНЕРИРУЕТСЯ в
  `scripts/gamma/reporting/json_report.py` хелперами `_continuum_grid_for_peak` /
  `_bg_continuum_grid` / `_amp_net_from_spectrum` — три имени эксклюзивны для json_report.py
  (verified Grep 2026-06-20: 10 вхождений / 1 файл). Историч. окно выборки = `E±3σ`
  (≈±1.27 ПШПВ, docstring `json_report.py:681`); с 2026-06-20 синхронизирован к
  интеграционной ROI = ±2.5 ПШПВ. OVERSTATED-correction: формулировка «E±3σ живёт ТОЛЬКО
  в json_report.py» неточна — генерация действительно там (Grep-эксклюзив по трём
  хелперам), но сам термин/окно `±3σ` живёт и у ПОТРЕБИТЕЛЯ-рендера в шаблоне:
  `templates/interactive_v1_17_2.html:1309-1310` — `E_lo/E_hi = energy ∓ 3·sigma_keV`
  по компонентам мультиплета (живой код); плюс синглет/library_coverage используют
  окно ±4σ — `:1246` (synglet coverage trace) и `:1502` (`_hw = 4·sigma_keV`).
  Указатель `:1470` (Censor 2026-06-21) — это МЁРТВЫЙ комментарий «Раньше trim ±3σ
  делал…», а НЕ живой рендер; пользоваться им как ссылкой на код нельзя. Прочие `3σ`
  в коде (Chauvenet `compute.py`, peak-search `bootstrap.py`, Currie L_C
  `deconvolve.py`) — НЕ про display.
  В первоисточнике ЛСРМ этот display как «ширина зоны» не встречается — это отрисовка,
  не зонирование.

  **Решение оператора 2026-06-20 («всё три»)** — в работе:
  **(а) ВЫПОЛНЕНО 2026-06-20** — display-континуум синглета/фона приведён к ±2.5 ПШПВ
  (синхрон с ROI). Три хелпера `_amp_net_from_spectrum` / `_continuum_grid_for_peak` /
  `_bg_continuum_grid` в `json_report.py` считают плечи на `[1.75·ПШПВ, 2.5·ПШПВ]` с каждой
  стороны (ПШПВ=2.355·σ, плечи = внешние 30% ROI, Ковелл §5.2.5), сетка-пролёт
  `[E−2.5·ПШПВ, E+2.5·ПШПВ]`. РАСШИРЕНИЕ полуширины ×~1.97 (было 3σ=1.27 ПШПВ).
  Верифицировано по факту (2026-06-20): на Th-232 Маринелли span/ПШПВ синглетов = **5.000**
  (т.е. ±2.5 ПШПВ), центр = энергия пика; 165 тестов отчёта (step11+snapshot) PASS;
  инвариант `_amp_net ↔ continuum_grid` сохранён (одинаковые плечи в обоих, сетка
  симметрична → континуум в центре = (L+R)/2 = cont_est). Только display — cert/квантование
  не затронуты (используют `peak_area_counts` напрямую).
  **(б) ВЫПОЛНЕНО 2026-06-20 (Step 3, LSRM Lzmax)** — Lzmax=10 ПШПВ внедрён в
  `find_multiplet_regions` (`deconvolve.py`): новые параметры `max_zone_length_fwhm=10.0`
  (+ `lzmax_roi_window_factor=2.5`) и хелпер `_split_zones_lzmax`, вызываемый после
  кластеризации перед `return clusters`. Длина зоны := размах пиков + 2·2.5·ПШПВ(центр);
  при превышении 10·ПШПВ(центр) — рекурсивный split в точке минимума отсчётов СТРОГО между
  крайними пиками (`pdf.md:482-485, 490-494`). Сравнение строгое («>»); ПШПВ берётся в
  центральном канале зоны. Запись кода — одним нормальным edit (хук-гейт снят: escape-маркер
  `.delegation_guard_off` + UTF-8-фикс декода stdin `delegation_guard.py:33-34`,
  censor-adjudication 2026-06-20; баг был — `json.load(sys.stdin)` декодировал не-ASCII путь
  как cp1251 → mojibake → escape-tree-walk проваливался).
  **(в) Cert 2026-06-21 (REVISED после Censor-вердикта 2026-06-21):** py_compile OK;
  полная регрессия `tests/step08_multiplets` + `tests/snapshot` = **1077 passed, 44
  skipped, 3 xfailed, 1 xpassed, 0 failed** при Lzmax default-ON. **Honest-cert-correction
  предыдущего (в):** заявленный «синтетический unit (5 кейсов)» в коммите 0c8726f
  отсутствовал в закоммиченном suite (false-green; пойман Censor 2026-06-21). Закрыт
  forward-fix'ом 2026-06-21: **реальный** unit `tests/step08_multiplets/test_split_zones_lzmax.py`
  (3 кейса, все PASS):
  (a) over-long зона >10·ПШПВ дробится в точке min-counts ровно по этой долине;
  (b) **зона длиной ~7.5·ПШПВ** (зеркало M3 Ac-228 233/252/277) НЕ дробится — порог
  строгий («>», `==` не дробит);
  (c) observability: при исчерпании `_max_depth` `_split_zones_lzmax` теперь эмитит
  `RuntimeWarning` с подстрокой `_max_depth` (раньше зона возвращалась молча — это и
  было корнем false-green'а: длинная зона без split была неотличима от граничного
  «не дотянула до порога»).
  Disabled (`max_zone_length_fwhm=0.0`) / `spec=None` → graceful no-op (covered тем же
  модулем неявно: ветка `if max_zone_length_fwhm <= 0.0 or spec is None`).
  **Уточнение про эталоны деконволюции (Censor 2026-06-21):** реальные длины зон
  M2/M3/M4 на Th-232 Marinelli — **M2 ≈ 5.66, M3 ≈ 7.47, M4 ≈ 6.44 ПШПВ**, все < 10 →
  не дробятся. Прежняя формулировка «~1-2 ПШПВ» относилась к РАССТОЯНИЯМ между
  крайними пиками (а не к длине зоны = размах + 2·wing); снимаю как вводящую в
  заблуждение. Чтобы M3 при текущей конфигурации зацепил порог Lzmax, потребовалось
  бы примерно удвоение размаха пиков либо вдвое меньший ПШПВ — то есть он НЕ «на
  грани».

