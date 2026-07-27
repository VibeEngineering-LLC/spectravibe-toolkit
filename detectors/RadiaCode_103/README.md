# detectors/RadiaCode_103/

Detector profile for **RadiaCode-103** — handheld CsI(Tl) scintillator-spectrometer.

- **Detector head:** RadiaCode-103, CsI(Tl) кристалл, ФЭУ-блок интегрирован в корпус
- **DAQ software:** RadiaCode Android app + RadiaCode Studio (export to XML/RC3)
- **Form factor:** карманный, 25×123×123 мм, ~150 г
- **Serial seen:** RC-103-012013 (оператор Дмитрий, 2026-06-13)
- **Shielding (HARD, operator-locked 2026-06-13):** ШИЛЬДА НЕТ. RC-103 — handheld
  без свинцового домика. Любая интерпретация пиков 73-90 keV как «Pb K-XRF от шильда»
  для этого детектора — **ошибка**. Pb K-X-rays на спектре RC-103 могут возникать
  ТОЛЬКО эндогенно: из внутренней конверсии (IC) в дочках Ra-226 / Th-232 в образце
  (Pb-214, Bi-214 → Pb K-Xα 74.97, Kα 72.80, Kβ 84.9) или как Th/U L-X-rays
  (Th Lα1 12.97, Pb Lα 10.55). НЕ ассоциировать с shielding XRF никогда.
- **Background rate (operator-locked 2026-06-13):** обычный фон RC-103 без шильда —
  **8-10 cps** в полном спектре (intergral 25-2800 keV). Если sample cps в 2-5 раз
  выше — образец активен. Sample 5.3 ч на 2026-06-13 показал ~12-15 cps integral
  → слабо активный (apatite Th + Ra + K).
- **First imported XML:** `Spectrum 13-06-2026 фон.xml` + `Spectrum 13-06-2026 (2).xml`
  (run folder `1_Version/v1.30.0/analysis_runs/2026-06-13_RC103_unknown_vs_bg20d/data/`)
- **Reference spectra archive:** `reference_spectra/<date>_<label>/` — копии sample+bg
  XML каждой пройденной сессии для regression и обучения. Текущие: апатит 2026-06-13,
  камень с Ra-226 2026-06-13.
- **Primary reference catalog (operator-locked 2026-06-13):**
  `references/Руководство_спектроскописта_V1.05_Соловьев_2024.pdf` + summary
  `references/руководство_спектроскописта_v1_05_summary.md` — справочник 14 нуклидов
  на RC-103 CsI(Tl) 1см³ + таблица L/K x-rays для всех элементов. Используется как
  empirical reference для шагов 4/7/10 pipeline.

---

## 1. Crystal-class & expected resolution

- **Crystal:** CsI(Tl)
- **Expected R(662 keV):** 8–10% FWHM (typical для CsI(Tl) handheld; хуже NaI 6–8%
  из-за non-prop)
- **FWHM(E) модель:** `FWHM(E) = k·√(E + α·E²)` с non-prop поправкой; деградировать
  до `FWHM = a + b·√E` если статистики мало
- **Non-prop:** CsI(Tl) имеет заметную нелинейность отклика на низких E (<200 keV);
  паспортная калибровка стораджа учитывает квадратичный член, но для тонкой E-cal
  на multi-line источниках рекомендуется degree ≤ 4 polynom

---

## 2. File format — RadiaCode XML

### Структура (root и ключевые узлы)

- Root: `ResultDataFile` (UTF-8 XML)
- Путь к спектру: `ResultDataFile / ResultDataList / ResultData / EnergySpectrum`
- Калибровка: `EnergySpectrum / EnergyCalibration` содержит:
  - `PolynomialOrder` (обычно 2)
  - `Coefficients / Coefficient` (массив, длина = PolynomialOrder + 1)
- Спектр: `EnergySpectrum / Spectrum / DataPoint` (ровно 1024 элементов)
- Time live: `EnergySpectrum / MeasurementTime` (seconds)
- Time real: вычисляется как `EndTime − StartTime` (на уровне `ResultData`)
- Background marker: `ResultData / BackgroundSpectrumFile` (пустой self-closing тег
  если raw, непустой если спектр уже net)
- Sample: `ResultData / SampleInfo / Name` и `Note`

### Поля

| XPath | Содержимое | Use |
|---|---|---|
| `//DeviceConfigReference/Name` | `RadiaCode-103` | детектор-ID |
| `//StartTime`, `//EndTime` | ISO 8601 без TZ | real time = EndTime − StartTime |
| `//NumberOfChannels` | `1024` (всегда для RC-103) | sanity-check |
| `//PolynomialOrder` | обычно `2` | количество ожидаемых coefficient-ов |
| `//Coefficients/Coefficient` | `a₀`, `a₁`, `a₂` | stored E-cal: `E(n) = a₀ + a₁·n + a₂·n²` |
| `//MeasurementTime` | seconds (live time) | live time (НЕ real) |
| `//Spectrum/DataPoint` | counts per channel | spectrum array (ровно 1024) |
| `//BackgroundSpectrumFile` | пустое self-closing — raw | если непустое — спектр уже net |
| `//SampleInfo/Name` | произвольное | sample-ID |
| `//SampleInfo/Note` | произвольное | metadata, часто пустое |

### Live vs real time

- **MeasurementTime** = live time (секунды) — что использовать для `cps`
- **EndTime − StartTime** = real time
- **Dead time fraction** = 1 − live/real

Пример (фон 2026-06-13): real = 20.4 сут = 1,762,690 s; live = 1,724,890 s →
dead = 2.14%. Для образцов с низкой загрузкой dead ≈ 0%.

---

## 3. Stored energy calibration — known instance

Серийный номер **RC-103-012013** (2026-06-13):

```
E(n) = 1.4309082 + 2.3888474·n + 0.00036652226·n²
```

Это **одинаковая калибровка** для всех записей этого прибора в одной серии
измерений (фон и образец, разнесённые на недели — без drift-а). Это упрощает
анализ: bg-subtract в channel-domain без resample.

**Crystal-class invariant**: формула одинакова у всех RC-103, но коэффициенты
индивидуальны (gain-stretching при производстве). НЕ переносить эти конкретные
коэффициенты на другой прибор.

### Channel ↔ Energy map (этот прибор)

| ch | E [keV] | Что там |
|---|---|---|
| 0 | 1.43 | threshold |
| 5 | 13.4 | electronic noise spike (маскировать) |
| 10 | 25.3 | ниже маскируется как noise |
| 30 | 73.4 | Pb K-Xα **endogenous** (IC от Pb-214/Bi-214 в Ra-226 chain) или Th L-X-rays |
| 35 | 85.4 | Pb K-Xβ endogenous (та же физика) |
| 97 | 237.3 | Pb-212 238.6 keV (Th-232 chain) |
| 213 | 510.9 | annihilation 511 + Tl-208 510.77 (Th-232) |
| 243 | 583.6 | Tl-208 583.19 keV |
| 254 | 609.4 | Bi-214 609.32 |
| 277 | 661.9 | Cs-137 661.66 |
| 380 | 911.2 | Ac-228 911.20 |
| 405 | 968.9 | Ac-228 968.97 |
| 562 | 1460.9 | K-40 1460.82 |
| 657 | 1764.4 | Bi-214 1764.49 |
| 953 | 2614.5 | Tl-208 2614.51 (Th-232 chain end) |
| 1023 | 2829 | end of range, overflow bin |

**Tl-208 2614.5 keV → ch ≈ 953** при stored cal (и ch 953 при пересборке
по якорям Pb K / K-40 / Tl-208 2614 на серии 2026-06-13). В диапазоне
RC-103, но **сильно размыт**: FWHM(2614) ≈ 165 keV (по модели
`FWHM = a·√E + b` с a=3.23, b=−7.89) → пик растянут на ~165 каналов
≈ 16% от ширины спектра. На большинстве источников **не выделяется**
из continuum'a peak-search'ем; даже когда виден, S/σ низкий. Поэтому
**primary anchor для Th-232 chain — Tl-208 583** (CI=0.95, S/σ=44),
**secondary** — Tl-208 510.77 (CLAUDE.md gotcha: Tl-208 510 доминирует
над annihilation 511 на любом Th-источнике), 2614.5 — только как
support при достаточной статистике.

---

## 4. Channel range gotchas

### Threshold spike (ch 0–6)

Низкоэнергетический шум усиления. На длинных live (20 сут фона) виден как
острый пик с count ~10^5–10^6. **НЕ интерпретировать как реальный пик.**

**Правило**: маскировать E < 25 keV (ch < 10) как electronic noise.

### Overflow bin (ch 1023)

Последний канал — накопитель «всё что > E_max». На фоне за 20 сут видел
`DP[1023] = 35989` — это не пик, это накопленные overflow events.
**НЕ интерпретировать как высоко-E линию.**

### Tl-208 2614.5 keV — в диапазоне, но размыт

`E(1023) = 1.43 + 2444.18 + 383.66 = 2829.27 keV` — формальный max RC-103.
Tl-208 2614.5 → **ch ≈ 953**, в диапазоне (не «вне», как было записано в
первой версии этого README). Но FWHM(2614) ≈ 165 keV растягивает пик на
~165 каналов с низкой амплитудой над continuum. На образце 5.3 ч
(2026-06-13) peak-search не выделил его как самостоятельный пик —
сливается с continuum'ом. Использовать **Tl-208 583** как primary anchor
для Th-232 цепочки; 2614 — support при долгих набирках (фон 20 сут
покажет его явно).

### Низкое разрешение на низких E

CsI(Tl) FWHM ≈ 8–10 keV в районе 80 keV → Pb K-Xα (73 keV) и K-Xβ (85 keV)
**сливаются** в один пик (как и на NaI). Не пытаться разрешить.

CsI(Tl) FWHM ≈ 50–80 keV в районе 600 keV → большинство «мультиплетов»
естественно слиты в один пик. Не плодить произвольные free-component
deconvolutions.

---

## 5. Workflow notes (HARD operator rule)

См. главный CLAUDE.md → «Порядок: фон → net → sample-операции». Для RC-103
конкретно:

1. Parse XML фона и образца раздельно
2. Сверить калибровки: для одного прибора в одной серии — одинаковы → можно
   subtract в channel-domain без resample. **Cross-check на 1+ якорной линии**
   (любой видимой K-X-ray или K-40)
3. `net_cps(ch) = sample_cps(ch) − bg_cps(ch)`
4. `σ_net = √(σ²_s + σ²_b)`
5. Все subsequent шаги (peak search, identification) — на net

**Если фон НЕ предоставлен**: спросить оператора, не запускать identification
на raw sample (false positives от K-40, Tl-208, Bi-214 окружения).

---

## 6. Parser implementation notes

Базовый Python-парсер RadiaCode XML использует `xml.etree.ElementTree`:
найти узел `.//ResultData`, оттуда `EnergySpectrum`, прочитать
`EnergyCalibration/Coefficients` (массив float), `NumberOfChannels` (int),
`MeasurementTime` (float, секунды live), и массив `Spectrum/DataPoint`
(int counts, ровно 1024 элементов). Метаданные: `StartTime`, `EndTime`,
`DeviceConfigReference/Name`, опционально `SampleInfo/Name` и `Note`.

`PYTHONIOENCODING=utf-8` обязателен (CLAUDE.md global). `ElementTree` парсит
UTF-8 XML напрямую — кодировка в файловом `open()` не нужна, она в XML-prologue.

---

## 7. Operator preview budget

См. CLAUDE.md → «Operator preview-запросы — бюджет ~20k токенов». На RC-103
preview-запросы («посмотри спектр», «полный анализ без активности») — один
скрипт solo, без `agent-a-math`, без 11 отдельных шагов. Полный пайплайн —
только по явной просьбе («ГОСТ Layer 2», «как для поверки»).

---

## 8. Known runs

| Дата | Run folder | Что измерялось | Вывод |
|---|---|---|---|
| 2026-06-13 | `1_Version/v1.30.0/analysis_runs/2026-06-13_RC103_unknown_vs_bg20d/` + `reference_spectra/2026-06-13_apatite/` | Образец **минерал апатит** (оператор, post-hoc) 5.3 ч vs фон 20 сут; без расчёта активности. Статистика слабая, R(662)≈11.4% | Th-232 chain в равновесии: Pb-212 238, Ac-228 338/(911+969 merged)/463 tent., Tl-208 510+583. Pb K-X-rays 70-75 keV — **эндогенные** (IC от Ra-226 chain в матрице апатита, шильда у RC-103 нет). Ожидается также Ra/U-238 chain (Bi-214 609/1764) и K-40 1460 — pipeline на 5.3 ч их не выделил из-за слабой статистики |
| 2026-06-13 | `demo_reports/2026-06-13_22-58_RC103_Stone_Ra_2026-05-30/` + `reference_spectra/2026-06-13_stone_Ra226/` | **Камень с Ra-226 (природный минерал)** 30 мин vs фон 199 сут; без расчёта активности. R(662)≈11% | **Ra-226 цепочка в равновесии**: ²²⁶Ra 186 (+ возможно ²³⁵U 185 — неразделимы на CsI), ²¹⁴Pb 295/352, ²¹⁴Bi 609, ²¹⁴Bi 1764 (слабо), ²⁰⁸Tl 583 (Th-232 background contribution от стен). **CRITICAL — operator HARD-LOCK 2026-06-13**: при обнаружении ²²⁶Ra 186 на природном образце ОБЯЗАТЕЛЬНО подозревать U-238/U-235 (см. §9 ниже). Pipeline ложно атрибутил 3 пика (94, 186, 292 кэВ) как **Ga-67** — это process-bug environment-context library (медицинские нуклиды не должны матчиться на природный камень). Ga-67 fix отложен. |

### Apatite (Ca₅(PO₄)₃(F,Cl,OH)) — что ожидать (operator-locked 2026-06-13)

Природный апатит — фосфатный минерал, regularly содержит:
- **Th-232 chain** (через изоморфное замещение Ca²⁺ → Th⁴⁺): Pb-212 238, Ac-228 338/463/911/969, Tl-208 510/583/2614
- **U-238 / Ra-226 chain** (тоже Ca²⁺ → U⁴⁺): Pb-214 295/352, Bi-214 609/1120/1764, Pb-210 46
- **K-40** (Ca²⁺ → K⁺): 1460.82 keV (всегда в природных Ca-минералах)
- **Pb K-X-rays 70-90 keV** — эндогенные от IC в дочках Ra-226 (Pb-214/Bi-214). НЕ shielding.
- Возможны Th L-X-rays 12-20 keV (ниже threshold RC-103)

На RC-103 с R(662)~11% многие линии **сливаются**:
- Pb-214 295 + 352 → один пик в районе 320-330 keV (но конкурирует с Ac-228 338)
- Bi-214 1120 + Ac-228 1110 + K-40 sum effects
- Tl-208 510 + annihilation 511 → один пик
- Ac-228 911 + 969 → один merged пик при ch ~370-400

Слабая статистика 5 ч на handheld → видны только сильнейшие лидеры (Pb-212 238, Ac-228 338/(911+969), Tl-208 583). K-40 1460, Bi-214 609/1764 на короткой набирке могут не выделиться выше continuum'а — это **не** означает их отсутствие в образце, только недостаточную статистику.

### Шаблоны результатов

- HTML preview-отчёт: `<run>/reports/report_full.html` (Arial, embedded PNG, без активности)
- State-файлы: `<run>/state/*.json` (parsed, net, ecal_final, fwhm_model, peaks_significant, final_attribution)
- Скрипты прогона: `<run>/scripts/01_*.py … 10_render_report.py`

---

## 9. HARD-LOCK: Ra-226 186 → suspect U-238/U-235 (operator-locked 2026-06-13, навсегда)

**Зафиксировано оператором Дмитрием 2026-06-13** при анализе спектра «Камень с Ra-226»
на RC-103: «когда ты нашел 186 Ra226 ты сразу должен подозревать наличие урана 238/235
с характерными сильными линиями в районе 63-186. Запомни навсегда».

### Правило

Если на природном образце (камень, грунт, минерал, фосфорит, гранит, апатит, monazit,
zircon, торий-содержащий песок, K-feldspar) обнаружен пик в окрестности **186 кэВ** —
**немедленно** подозревать присутствие U-238 и U-235 в матрице и искать
характеристические U-линии в диапазоне **25-200 кэВ**:

| E, кэВ | Эмиттер | Цепочка | Что искать на RC-103 |
|---|---|---|---|
| **25** | ²²⁷Th | x-ray escape от U K x-rays | малый пик у edge'а спектра |
| **63.30** | ²³⁴Th | ²³⁸U → ²³⁴Th β⁻ | плечо/пик 50-70 кэВ — часто слит с Pb K-X-rays 73 |
| **92.4** | ²³⁴Th | ²³⁸U → ²³⁴Th | под U K x-rays |
| **~95** | U K x-rays | IC в actinide chain | сильный пик 90-100 кэВ на природном U |
| **143.76** | ²³⁵U | ²³⁵U прямой | merged с 163 в широкий пик 140-170 |
| **163.36** | ²³⁵U | ²³⁵U прямой | та же merge |
| **185.72** | ²³⁵U | ²³⁵U прямой (I_γ_per_decay 57.2%) | **сливается с ²²⁶Ra 186.21** при R(186)~17% на RC-103 |
| **1001** | ²³⁴ᵐPa | ²³⁸U → ²³⁴ᵐPa | слабый пик ~1000 кэВ выше continuum'а |

### Эмпирическая основа

См. `references/руководство_спектроскописта_v1_05_summary.md` §1, стр. **32 PDF**
(«Разница U(Природный) − Ra²²⁶, стр. №2») — на RadiaCode-103 разностный спектр явно
показывает 63 / ~95 / 143+163 / 185 как U-only сигнатуру, остающуюся после вычитания
чистого ²²⁶Ra.

### Reporting convention на RC-103 при подозрении

- Пик в окрестности 186 кэВ → НЕ называть «²²⁶Ra 186» автоматически.
  Правильно: «**²²⁶Ra+²³⁵U 186 keV (merged on CsI(Tl) R~11-17%)**» или
  `Ra-226_U235_186` в JSON-полях.
- В operator-facing отчёте (`<run>/reports/report_full.html`) — пометить, что
  разрешить ²³⁵U 185 vs ²²⁶Ra 186 на RC-103 невозможно. Для дифференциации нужны:
  (a) долгая набирка HPGe (>24 ч), (b) количественный изотопный анализ через
  estimate ²³⁵U/²³⁸U = 0.00725 (природное отношение), (c) cross-check через
  ²³⁴Th 63 / U K x-rays ~95 / ²³⁴ᵐPa 1001 — если эти линии видны выше bg,
  U-238 родитель **подтверждён**.

### Anti-pattern (что было до фиксации правила)

❌ В сессии 2026-06-13 на «Камень с Ra-226» pipeline на 30-мин набирке атрибутил пик
**~186 кэВ** как **Ga-67 184.6** (медицинский изотоп, T½ 3.26 сут, артефактно
матчился из-за отсутствия environment-context фильтра). Категориальная ошибка:
медицинский короткоживущий изотоп на природном минерале с Ra-226+ДПР в равновесии.
Ga-67 атрибуция — process-bug pipeline, fix в `gamma/id/candidates.py` или
`scripts/run_plan_a.py` отложен.

### Связь с другими правилами

- **handheld no-shield rule** (CLAUDE.md): Pb K-X-rays 70-90 кэВ на RC-103
  всегда эндогенные → не отбрасывать как «фон шильда», они часть U/Ra-226 сигнатуры.
- **Operator preview budget ~20k**: при подозрении на U — НЕ запускать полный
  11-шаговый пайплайн через `agent-a-math`. Один скрипт solo + проверка 4-5 U-линий
  через `references/руководство_спектроскописта_v1_05_summary.md` §2.