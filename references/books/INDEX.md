# Knowledge Library Index — gamma-spectrum-analysis

> **F-138 / v1.17.7 (закреплено навсегда)** — централизованная библиотека
> справочных PDF-документов. Все алгоритмические/методологические/нормативные
> ссылки в коде должны указывать на файлы из этой библиотеки.
>
> **F-293 / v1.17.19 (закреплено навсегда)** — **файлы библиотеки физически
> перенесены** из `references/books/` в папку
> `books_library/` **в корне проекта** (рабочая копия):
> `D:\...\0_Work\gamma-spectrum-analysis\books_library\`. В
> `references/books/` остался только этот `INDEX.md` — каталог с
> метаданными.
>
> **Архивы библиотеки** хранятся отдельно от релизных архивов проекта в
> `D:\...\1_Version\books_library\gamma-books_vYYYY-MM-DD.zip`,
> упаковываются вручную через `python scripts/build_books_archive.py`.
>
> См. также `scripts/verify_books_inventory.py` (auto-sanity-check перед
> релизом). Путь настраивается через env `GAMMA_BOOKS_LIBRARY_DIR`
> (default: `<root>/books_library`).
>
> Структура: `books_library/<short-id>.pdf` (или `.ppt/.pptx/.docx`) +
> аннотации в этом файле. При добавлении новой книги: положить файл в
> `books_library/`, добавить запись ниже (id / автор / год / тема /
> страничный диапазон ключевых разделов).
>
> **F-150 / v1.17.8 (расширено в v1.17.19, F-333 в v1.18.18.6)** — папка
> `books_library/` **полностью исключена** из релизных архивов
> `gamma-spectrum-analysis_v*.zip`. Дополнительно из архива исключены
> `references/_converted_tmp/` (build-cache PPT→PPTX, ≈40 МБ).
>
> **F-333 / v1.18.18.6** — rendered PNG страниц для multimodal-чтения
> перенесены из `references/_extracted_corpus/_pages/` (legacy) в
> `books_library/_corpus_pages/{lsrm_activity_2014,vartanov}/` (≈63 МБ,
> 301 файл). PNG генерируются скриптами `scripts/_render_scanned.py`,
> `scripts/_ocr_activity.py`, `scripts/_ocr_vartanov.py` и более не
> засоряют дерево `references/`. Текстовые экстракты
> `references/_extracted_corpus/*.md` сохраняются для RAG. Список
> исключений закреплён в `scripts/build_release_archive.py`.
>
> **F-157 / v1.17.9.4 (закреплено навсегда)** — для Gamma-1C
> литература ЛСРМ имеет повышенный приоритет. На конфликте:
> ЛСРМ > Будыка > Gilmore (числа) / ЛСРМ-алгоритм без числа → Будыка/Шендрик.
>
> **v1.17.9.4 структура** — `Документация ЛСРМ/` разложена по 5 тематическим
> категориям (см. раздел B ниже).

---

## Список

> **Inventory note (DEEP-10, v1.26.2)** — три записи ниже (Lsrm_algorithmic_foundations.pdf, Мощность дозы..pdf, Прецизионные измерения.pdf) физически расположены в `books_library/Документация ЛСРМ/01_methodology_pdf/` и `02_topical_pdf/` (не в корне `books_library/`). Записи вынесены из заголовков `### N.` чтобы не включаться в top-level inventory (verify_books_inventory.py). Это **исторический artefact** перемещения F-293. Файлы проиндексированы в `books_library/Документация ЛСРМ/INDEX.md`. `verify_books_inventory --strict` должен проходить без MISMATCH.
>
> **Top-level inventory count** (`references/books/INDEX.md` vs `books_library/` root): **13 файлов** (entries #2–#16 ниже, без учёта вложенных LSRM-документов).

### 1-nested. Lsrm_algorithmic_foundations.pdf (in Документация ЛСРМ/01_methodology_pdf/)
- **Полное название**: «Алгоритмические основы программ обработки спектрометрической информации SpectraLine», 2022
- **Расположение (фактическое)**: `books_library/Документация ЛСРМ/01_methodology_pdf/Lsrm_algorithmic_foundations.pdf`
- **Издатель**: ООО «ЛСРМ» (Москва)
- **Ключевые разделы**:
  - §7 (стр. 7-1…7-8) — обработка результатов измерений, формула σ = max(scatter, weighted-mean)
  - §8.4 (стр. 8-3…8-7) — калибровка по форме линии, **peak-image модель** (Гаусс + tail + Compton-step)
  - §8.4.4 — Compton ступенька, формула `step(x) = (A_step/2)·erfc(...)`, доля `h_step` ≈ 0.03 для NaI
  - §9 — калибровка по энергии, ёмкость окна идентификации, седьмая 7-line проверка ЕРН
  - §9.4 — peak-list a-priori representation фоновых линий (F-96)
  - §9.7 — Compton continuum, step+linear модель (Gilmore §9.7 ссылается на ту же физику)
  - §10 — каскадное суммирование (TCS) — формула `C(E) = 1 / (1 − Σ p·ε_T)`
  - §11 — ISO 11929 правило upper limit при σ_A/A > 50%
  - §12 — шаблонный метод (template method, sensitivity matrix R)
  - §13 — квази-шаблонный метод (full-spectrum simultaneous fit)
  - §14.2 — Dose Contribution metric (F-61)
  - §15 — эффективное мёртвое время `t_d = A·Σy + B·Σ(y·i)`
- **Цитируется в коде**:
  - `scripts/gamma/peaks/coupled_multiplet.py` — F-117, F-127, F-133
  - `scripts/gamma/peaks/peak_image.py` — F-90
  - `scripts/gamma/physics/cascade_summing.py` — K-17, K-18
  - `scripts/gamma/physics/dead_time.py` — F-95
  - `scripts/gamma/physics/bg_lines_apriori.py` — F-96
  - `scripts/gamma/activity/template_method.py` — F-100
  - `scripts/gamma/activity/quasitemplate.py` — F-98

### 2. `lsrm_format_specification.pdf`
- **Полное название**: «Описание формата файлов SpectraLine `.spe`»
- **Издатель**: ООО «ЛСРМ»
- **Ключевые разделы**:
  - §7.5.2.1 — структура header (KEY=VALUE\r\n), таблица пиков, ZONES
  - Поля метаданных: SAMPLEMASS, SAMPLEVOLUME, MATERIAL (с `Ro` плотностью)
  - Энергетическая калибровка `ENERGY=N,a0,a1,a2,...` и orthopoly зоны
  - FWHM модель `FWHM=...`
- **Цитируется в коде**:
  - `scripts/gamma/io/lsrm_spe.py` — весь модуль
  - F-130 (auto-density extraction)

### 3. `pgs-gilmore-2008.pdf`
- **Полное название**: «Practical Gamma-ray Spectrometry», Gordon R. Gilmore, 2nd Edition, 2008
- **Издатель**: John Wiley & Sons
- **Ключевые разделы**:
  - Ch. 2 — Радиоактивный распад и фундаментальные основы (T½, branching ratios)
  - Ch. 4 — Сцинтилляционные детекторы (NaI(Tl), CsI, BGO, LaBr, CeBr — характеристики)
  - Ch. 5 — Эффективность детектирования, peak-to-total ratio (Table 8.4 — данные для NaI 3×3″)
  - Ch. 6 — Background, environmental gamma
  - Ch. 7 — Спектральные взаимодействия γ-фотонов
  - Ch. 8.5 — Каскадное суммирование (TCS), Table 8.4 P/T ratio
  - Ch. 9.3 — Методы поиска пиков (Mariscotti, matched-filter convolution)
  - Ch. 9.7 — Континуум, step+linear модель под фотопиком
  - Ch. 10 — Source preparation, geometry effects
  - Ch. 11 — Calibration sources & techniques
- **Цитируется в коде**:
  - `scripts/gamma/peaks/convolution_search.py` — F-124 (matched-filter §9.3)
  - `scripts/gamma/peaks/deconvolve.py` — F-33, F-34, F-117 (step+linear §9.7)
  - `scripts/gamma/physics/cascade_summing.py` — Table 8.4 (P/T NaI)
  - `scripts/gamma/activity/compute.py` — σ propagation §5.7.2

### 4. `Shendrik_Scintillators_pt1.pdf`
- **Полное название**: «Введение в физику сцинтилляторов. Часть 1», Шендрик
- **Тема**: основы сцинтилляции, механизмы возбуждения и излучения, NaI(Tl), CsI(Tl), световой выход, временные характеристики
- **Релевантно для skill'а**: понимание формы пика (Гаусс + low-energy tail) — связано с потерей заряда / световыхода и малоугловым рассеянием. Подтверждает физическую основу F-133 (per-line step) и F-127 (T(E) tail).

### 5. `Shendrik_Scintillators_pt2.pdf`
- **Полное название**: «Введение в физику сцинтилляторов. Часть 2», Шендрик
- **Тема**: продвинутая физика сцинтилляторов, разрешающая способность, дозиметрические применения, новые материалы (LaBr, CeBr, SrI, CLYC)
- **Релевантно для skill'а**: модель FWHM(E) для NaI — `FWHM²(E) = a + b·E + c·E²` (F-125 калибровка), сравнение с другими сцинтилляторами для будущих детектор-расширений.

### 7. `Budyka_Spektrometriya_ioniziruyushchikh_izlucheniy._Osnovnye_ponyatiya_2021.pdf` (v1.17.9.3)
- **Полное название**: «Спектрометрия ионизирующих излучений. Основные понятия и терминология» (учебно-методическое пособие)
- **Автор**: А.К. Будыка
- **Издатель**: НИЯУ МИФИ, 2021
- **ISBN**: 978-5-7262-2794-8
- **Объём**: 144 стр.
- **Тема**: глоссарий 500+ терминов γ/β/нейтронного спектрометрии с английскими эквивалентами (алфавитный порядок).
- **Использование**: lookup-словарь для RU/EN narrative (F-108) + точные определения для F-rule docstrings.

### 8. `Budyka_Spektrometriya_ioniziruyushchikh_izlucheniy_2021.pdf` (v1.17.9.3)
- **Полное название**: «Спектрометрия ионизирующих излучений. Гамма-спектрометрия» (учебное пособие)
- **Автор**: А.К. Будыка
- **Издатель**: НИЯУ МИФИ, 2021
- **Объём**: 225 стр.
- **Ключевые разделы**:
  - Гл.2 §2.3 (стр. 22-24) — Энергетическое разрешение R(E), формула R²=R²_stat+R²_intr+R²_PMT
  - Гл.2 §2.4 (стр. 25-30) — Функция отклика и аппаратурная форма линии
  - Гл.3 §3.2 (стр. 35-52) — Аппаратурная форма линии γ-излучения (полная декомпозиция)
  - Гл.3 §3.3 (стр. 53-57) — Эталонные источники (ОСГИ/ИИГИ маркировки)
  - Гл.4 §4.4 (стр. 70-82) — Сцинтилляторы + PMT параметры (QE, M, LY=38000 ф/МэВ)
  - Гл.5 §5.6 (стр. 136-141) — МКА: ADC, DNL ≤1%, INL ≤0.05%
  - Гл.5 §5.7 (стр. 141-142) — Стабилизаторы спектра (gain drift correction)
  - Гл.6 (стр. 155-165) — Anti-coincidence, Compton suppression, парные спектрометры
  - Гл.7 §7.3-7.4 (стр. 178-183) — Поиск пиков: статистические, Mariscotti, correlation
  - Гл.7 §7.5 (стр. 183-188) — Центроид: σ(centroid)=FWHM/(2.355·√N)
  - Гл.7 §7.7 (стр. 193-196) — Площадь: Currie net, Covell, Sterlinski TPA
  - Гл.7 §7.8 (стр. 196-203) — Разрешение мультиплетов: peak stripping vs nonlinear fit
  - Гл.8 §8.1 (стр. 205-212) — High count rate: pile-up, dead time
  - Гл.8 §8.2 (стр. 212-223) — Low-activity: MDA, low-level counting

### 9. `Минимальная детектируемая активность. Основные понятия и определения.pdf` (v1.17.9.3)
- **Объём**: 6 стр. (статья)
- **Тема**: 4 типа MDA по Currie: L_C (critical), L_D (detection), L_Q (determination ≈ 10·σ_0), quantification level. Paired-blank `L_D=2.71+4.65·√n_0` vs well-known-blank.
- **Использование**: обоснование для F-169 Currie paired-blank в плане Tier 1 v1.17.10.

### 10-nested. Мощность дозы. Методика расчета из спектра гамма-излучения..pdf (in Документация ЛСРМ/01_methodology_pdf/)
- **Полное название**: «Мощность дозы. Методика расчёта из спектра гамма-излучения»
- **Издатель**: ВНИИФТРИ + ООО «ЛСРМ»
- **Год**: 2000
- **Объём**: 10 стр.
- **Тема**: расчёт Ḣ*(10) из спектра — `Ḣ*(10) = Σ N_k·g(E_k)/(ε·t·K)`, где g(E_k) из ICRP-74. Нормы погрешности: ≤30% для прецизионных установок, ≤50% для полевых.
- **Использование**: обоснование для **F-177 (NEW dose_rate.py)** в плане Tier 1 v1.17.10.

### 11-nested. Прецизионные измерения.pdf (in Документация ЛСРМ/02_topical_pdf/)
- **Полное название**: «Прецизионные измерения. Образцовые и калибровочные источники»
- **Объём**: 46 стр.
- **Ключевые разделы**:
  - Калибровочные источники: рабочие δ=10-15%, эталон 2-го разряда δ=4-6%, 1-го разряда δ=3-4%, первичный <1%.
  - Поправка на самопоглощение: plane-slab ≤5% при d/D<0.5; Monte-Carlo для δ<3%.
  - TCS + random summing: `C_TCS = 1/(1-Σp·ε_T)`, `n_RSC ≈ 2·n1·n2·τ`.
  - K-фактор combined uncertainty: σ_K_attn ±3%, σ_K_TCS ±2%, σ_K_decay <0.1%.
- **Использование**: обоснование F-179 (K-budget) и F-180 (cert grade) в Tier 1.

### 6. `Experiment_results_analysis.pdf`
- **Полное название**: «Анализ и представление результатов эксперимента»
- **Тема**: статистическая обработка данных, неопределённости измерений, нормальное и пуассоновское распределение, ISO формулировки погрешностей, χ² критерии согласия
- **Релевантно для skill'а**:
  - σ propagation в `compute.py` (F-91 LSRM §7 max(scatter, weighted-mean))
  - χ²/ν критерий качества fit'а в `coupled_multiplet.py`
  - ISO 11929 правила upper limit (`mda.py`)
  - Currie L_C/L_D в `convolution_search.py`

### 12. `ortec_gammavision_v9_a66.pdf` (v1.17.9.5 ← новый)
- **Полное название**: «GammaVision® / Maestro-PRO® — Gamma-Ray Spectrum Analysis and MCA Emulators for Microsoft Windows. Software User's Manual» (A66-BW, A66SV-BW, A66MP-BW, Software Version 9)
- **Издатель**: ORTEC® (Advanced Measurement Technology, Inc. / AMETEK), Oak Ridge, TN, USA
- **Год / редакция**: 2020 / Manual Revision **M** (ORTEC Part No. 783620, 0220)
- **Объём**: 83 стр. (с обширными приложениями)
- **Категория**: vendor-software methodology — третий международный канон одного уровня с Gilmore (USA), параллельно ЛСРМ (RU)
- **Релевантность для skill'а**: содержит **6 движков анализа**, в том числе `NAI32` — единственный коммерческий движок специально для NaI/CsI/LaBr (низкое разрешение), что делает его прямым аналогом скилу для Gamma-1C.
- **Ключевые разделы (главы)**:
  - Ch. 6.2 (стр. 231-237) — **6 движков анализа**: WAN32, GAM32, NPP32, ENV32, **NAI32**, ROI32. Decision matrix (Table 6.2.2). Library Reduction (key-line + fraction-limit).
  - Ch. 6.3.1 (стр. 238-242) — **Background methods**: Automatic (5/3/1-point adaptive), X-Point Average, **X.X·FWHM** (по умолчанию ширина = ceil(X·FWHM_кан)).
  - Ch. 6.3.2 (стр. 243-250) — Peak Area Singlets: **Total Summation**, **Directed Fit**, **ISO NORM Singlet**.
  - Ch. 6.3.4 (стр. 251) — Peak Uncertainty (включая **ZDT-спектры**).
  - Ch. 6.3.7 (стр. 255-258) — **Mariscotti Peak Search**: smoothed second-difference, weights k_i, j=4 (regular) / j=9 (wide), wide-filter порог 0.15 keV/канал; PEAK SEARCH SENSITIVITY (S, F=1.35/1.0). Equations 53-57.
  - Ch. 6.5 (стр. 259-267) — **Multiplets**: deconvolution width = 3.08·FWHM между пиками; multiplet region = ±1.5·FWHM; **3 типа background**: stepped (slope test), parabolic (E<200 кэВ + 3 точки <line), straight-line (default). Library-Based Peak Stripping (auto + manual).
  - Ch. 6.7 (стр. 268-271) — **Nuclide Activity** (eq. 59-60): A = N_E·TDC·Mult·RSF / (LT·ε_E·Br·GeoFac·AttCorr·Div·s). Weighted-mean (eq. 61) + **Activity Range Test** (eq. 60a, F=MIN(r·σ₁, F_max)).
  - Ch. 6.9 (стр. 273-282) — **19 методов MDA**:
    1. Traditional ORTEC, 2. Critical Level ORTEC, 3. Suppress, 4. **KTA Rule**, 5. Japan 2σ, 6. Japan 3σ, 7. **Currie Limit**, 8. RISO MDA, 9. LLD ORTEC, 10. Peak Area, 11. GIMRAD/DIN 25 482 Air, 12. **Regulatory Guide 4.16 (USA)**, 13. Counting Lab USA, 14. **DIN 25 482.5 Erkennungsgrenze (CL)**, 15. **DIN 25 482.5 Nachweisgrenze (MDA ≈ 2·Method 14)**, 16. **EDF France** (CEA-R-5506), 17. **NUREG 0472**, 18. **ISO Decision Threshold (CL)**, 19. **ISO Detection Limit (MDA, 11929)**.
  - Ch. 6.10 (стр. 283-294) — Corrections: DDA (eq. 83), Decay Correction (eq. 84), Decay During Collection, **PBC** (Peaked Background — by Nuclide / by Energy, eq. 86-87), Geometry (eq. 88, .UFO-based table), Absorption (External eq. 89 / Internal — ratio + table methods, ASTM E181-82).
  - Ch. 6.11 (стр. 294-296) — **Random Summing** (RSF, eq. 91; experimental slope from 88Y+137Cs pair).
  - Ch. 6.12 (стр. 296-308) — **Reported Uncertainty** (quadrature): Counting, Additional Normal, Random Summing, Absorption, Nuclide, Efficiency (TCC-poly, Interpolative, Linear/Quadratic/Polynomial, Matrix), Geometry, Uniform, User-Defined, Sample Size, PBC (single + multi).
  - Ch. 6.16 (стр. 313) — **True Coincidence Correction (TCC)** — full peak/total efficiency model (ref. Knoll 3rd ed., Gilmore-Hemingway, ANSI N42.14-1991).
  - Ch. 6.17 (стр. 314-326) — **ISO NORM Implementation (ISO 11929:2010)**: модель y=N·κ; Critical Level (eq. 143-144); Peak MDA (eq. 145-153, Special Cases f=1); MDA-to-CL Ratio (eq. 154-163, R_max ограничение); Nuclide Activity / Best Estimated Activity / Lower-Upper Confidence Limits.
  - Ch. 6.18 (стр. 326-328) — EDF Gamma Total Analysis (Cs-equivalence).
  - Ch. 8 (стр. 359-374) — **QA (ANSI N13.30 + N42.14)**: Total Background, Total Activity, Average FWHM ratio, Average FWTM ratio, Average Library Peak Energy Shift; warning/acceptance limits, control charts.
- **Цитируется в коде** (после реализации Tier 1, v1.17.10+):
  - F-167 ID window — cross-check с ORTEC Library Match Width default 0.5·FWHM (ID matching) vs Mariscotti search vs deconvolution 3.08·FWHM
  - F-169-REVISED МИА — добавить ORTEC Method 18/19 (ISO 11929) и Method 7 (Currie) как secondary validators
  - F-170 PBC — выровнять с ORTEC PBC By Energy/By Nuclide (eq. 86-87)
  - F-176 Random Summing — ORTEC eq. 91 канон
  - F-178/F-179 TCC, ISO NORM — ORTEC §6.16-6.17
- **Конфликты с ЛСРМ** (новые, требуют разрешения в audit/v1_17_9_5_ortec/):
  - ORTEC default Library Match Width = **0.5·FWHM** (NaI32) vs ЛСРМ k=**1.5·FWHM** ([LSRM-Algo-9]) — разные контексты: ORTEC — ID match; ЛСРМ — peak deconvolution окно. Требует чёткой документации semantics.
  - ORTEC multiplet criterion = **3.08·FWHM** (Section 6.5.1) vs Gilmore default 3.0 vs ЛСРМ peak-overlap-range >= 1.0. Близкие, но не идентичные.
  - ORTEC `Activity Range Test` (eq. 60a) — отсутствует прямой аналог в ЛСРМ методике; кандидат на NEW F-rule.

### 13. `GammaDet.docx` (v1.17.22 ← новый top-level каталог)
- **Полное название**: «GammaDet — конструкторская и пользовательская документация на Гамма-1С/Колибри-семейство»
- **Формат**: DOCX (extracted MD доступен)
- **Тема**: техническое описание Gamma-1C сцинтилляционного блока (NaI(Tl) 63×63, PMT, USB-инструментальная оболочка)
- **Использование**: справка по конструкции детектора и базовым тех-характеристикам (FWHM model, активная зона, мёртвая область).

### 14. `spectralinexx_2.0_basic_functions_rus.pdf` (v1.17.22 ← новый top-level каталог)
- **Полное название**: «SpectraLine XX 2.0 — Основные функции (руководство пользователя)»
- **Издатель**: ООО «ЛСРМ»
- **Формат**: PDF (русский)
- **Тема**: пользовательский ман SpectraLine — описание UI, форматов файлов (.spe / .cpt / .efr / .lib), сценариев анализа
- **Использование**: справка по interoperability — поведение референсного коммерческого тула для сверки результатов skill'а.

### 15. `Идентификация и расчет активности.pptx` (v1.17.22 ← новый top-level каталог)
- **Полное название**: «Идентификация и расчёт активности» (презентация-методичка ЛСРМ)
- **Формат**: PPTX
- **Тема**: workflow — отбор пиков → совмещение с библиотекой → расчёт активности (LSRM canonical flow)
- **Использование**: cross-check описанной workflow vs реализованной в `staged_pipeline.py` (Step-1..11).

### 16. `2_Спектрометрия_2012_05 Илья.ppt` (v1.17.22 ← новый top-level каталог)
- **Полное название**: «Спектрометрия — 2012, лекция (Илья)»
- **Формат**: PPT (legacy Office)
- **Тема**: учебная презентация по спектрометрии (введение, energy/efficiency calibration на примерах)
- **Использование**: исторический материал — НЕ источник методологии, информационный пласт RAG.

---

### B1. `Документация ЛСРМ/01_methodology_pdf/` — официальные методики ЛСРМ
**Приоритет 1 по F-157** — обязательная справка для Gamma-1C.

| Файл | Год | Что |
|---|---|---|
| `Lsrm_algorithmic_foundations.pdf` | 2022 | Алгоритмические основы SpectraLine (дубль — основной в корне `books/`) |
| `Активность в счетных образцах. Методика измерений на гамма-спектрометрах с использоваонием ПО СпектраЛайн.pdf` | 2024 | **Основная методика ЛСРМ для активности в счётных образцах** (8.5 MB) |
| `Мощность дозы. Методика расчета из спектра гамма-излучения..pdf` | 2000 | Методика расчёта Ḣ\*(10) из спектра (ВНИИФТРИ + ЛСРМ) |

### B2. `Документация ЛСРМ/02_topical_pdf/` — тематические работы (PDF/DJVU)

| Файл | Что |
|---|---|
| `5_1_Прецизионные методы.pdf` | LSRM презентация: прецизионные методы (тезисы) |
| `5_2_Практическая спектрометрия-ядерные материалы.pdf` | LSRM презентация: практ. спектрометрия — ядерные материалы |
| `Прецизионные измерения.pdf` (4.5 MB, full) | Полная версия — прецизионные методы, K-budget, δ-разряды источников |
| `Практические_методы_сцинтилляционной_гамма_спектрометрии_Вартанов.djvu` | Вартанов — классический учебник по практ. сцинтилляционной γ-спектрометрии (NaI) |

### B3. `Документация ЛСРМ/03_lectures/` — лекции (PPT/PPTX/DOC, 22 файла)
Учебные курсы Кувыкина, ЛСРМ, презентации экспертной системы и пр.
**Не в RAG** (визуальный материал, ручная справка).

### B4. `Документация ЛСРМ/04_lab_works_docs/` — лабораторные работы (PDF/DOCX, 94 файла)
Подпапки:
- `EffMaker/` — Лабораторная работа по построению ε(E) (EffMaker)
- `GammaLab/` — Лабораторная работа по GammaLab
- `NuclideMaster/` — Обзор NuclideMaster + составление .lib (важно для F-220 candidate)
- `Калибровки/Калибровки/` — лабораторные по энергии/полуширине/эффективности + первичная калибровка
- `Измерение активности_учет фона/` — измерение активности и учёт фона
- `Сценарии_протоколы/` — SpectraLine scripting / автоматизация
- `Учет самопоглощения/` — прецизионная обработка спектров (самопоглощение)
- `Файлы паспортов/` — создание и использование .cpt паспортов LSRM
- `Расчет эффективности EffCalcMC/` — Monte-Carlo расчёт ε(E)
- `_текст_misc/` — текстовые версии презентаций ЛСРМ

### B5. `Документация ЛСРМ/05_lab_data/` — реальные данные ЛСРМ (236 файлов)
`3_Lab_sums-2012.03.01/` — Lab_Sums + Lab_Гамма-1П №SN-07:
- `Data/` — `.efr` (efficiency), `.lib` (nuclide), `.cpt` (passport), `.efa` (eff matrix)
- `Spe/` — реальные `.spe` спектры (Самопоглощение, Относительная эффективность, HighLoad)
- `Scenario/` — сценарии SpectraLine
- `*.ini`, `*.cnf` — конфигурация

**Использование**: эталонные данные для будущей регрессии (test fixtures).

---

---

## Унифицированные библиографические ссылки (v1.17.9.4)

Полные записи по **ГОСТ Р 7.0.5–2008** вынесены в отдельный файл:
[`../REFERENCES.md`](../REFERENCES.md) (21 запись + соглашения проекта).

**Формат ссылки в коде/F-rules/аудитах**: `[RAG-id, §раздел]`
(например `[LSRM-Algo, §10]`, `[LSRM-Activity-2014, Прил.5]`, `[Будыка-2021, §7.5]`).

**Формат ссылки в пользовательских отчётах (HTML/PDF)**: `[№, с. N]` с раскрытием
в разделе «Список использованных источников» (затекстовая ссылка по §7.5
ГОСТ Р 7.0.5–2008).

---

## Правила использования (закреплено F-138)

1. **Любая методологическая ссылка** в коде должна указывать на конкретный раздел из этой библиотеки (например `LSRM §8.4.4`, `Gilmore §9.7`, `Shendrik pt.1 гл.3`).
2. **При добавлении нового F-rule**, методологическое обоснование должно цитировать конкретный раздел из библиотеки, иначе F-rule отклоняется как «не методологически обоснованный».
3. **Книги в `references/books/`** — read-only архив; редактирование запрещено.
4. **Аннотации** в этом файле обновляются при каждом добавлении/удалении PDF.

## История перемещений

- **2026-05-29 (F-138, v1.17.7)** — создана `references/books/`, перенесены 2 LSRM PDF из `references/`, добавлены 4 новые книги (Gilmore 2008, Шендрик pt1+pt2, Анализ и представление результатов эксперимента).
