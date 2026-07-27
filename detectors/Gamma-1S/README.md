# detectors/Gamma-1C/

Detector-specific assets for the **Gamma-1C** spectrometric complex:

- **Detector head:** УДС-ГЦ-63×63 (NaI(Tl) 63×63 mm crystal, БДЭГ-63×63-USB) by Aspect
- **DAQ software:** Lsrm SpectraLine
- **Aliases:** Колибри-1М, Гамма-1С, БДЭГ-63×63 (see `data/aliases.json` → `detector.Gamma-1C`)

---

## 1. Taxonomy lock — Crystal vs Station (2026-06-05, user)

Two distinct levels — never conflate:

| Level | What it is | Carries… |
|---|---|---|
| **Crystal-class** | физика сцинтиллятора (тип + размеры) | универсальные паттерны линий нуклида (относительные позиции, BR-нормированные интенсивности после FWHM-нормализации). Шаблоны переносимы между всеми экземплярами этого класса. |
| **Station-instance** | конкретный экземпляр детектора + электроника + serial | индивидуальные FWHM-полином, efficiency-кривая, drift во времени, dead-time. Калибровка не переносима. |
| **Resolution (FWHM-кривая)** | свойство конкретного экземпляра | у каждого NaI 63×63 разная (разные кристаллы — разный световыход, ФЭУ-токи). Влияет на ширину окон поиска пиков, но не на отношения интенсивностей. |

**Gamma-1C как сущность**: единственная физическая **station-instance** этого
проекта (LSRM, BДЭГ-63×63-USB №SN-01 + БДЭГ-63×63-USB №SN-04 — оба
serial-числа относятся к **одной** станции в разные эпохи / Поверка-циклы),
а её **crystal-class = NaI 63×63**. См. `audit/_rag/visual_templates/SCHEMA.md:10`
и `scripts/gamma/io/lsrm_spe.py:29` («DETECTOR=… — e.g. "Гамма-1С" (NaI 63×63)»).

### Canonical name fixation (HARD-LOCK 2026-06-05, user)

> «Гамма-1с (кириллица) = Gamma-1C (латиница)»

**Canonical = `Gamma-1C`** (омоглиф-mapping: кириллическая «С» → латинская «C»).
`detectors/Gamma-1C/` — **единственная** папка для этой станции.

**Do NOT** create `detectors/Gamma-1S/`, `detectors/Гамма-1С/`, или любые иные
варианты — это duplicates одной и той же sub-tree.

`data/aliases.json:detector` содержит *два* canonical токена (`Gamma-1C` и
`Gamma-1S`) — это legacy v1.11.1 registry decision (`_meta.version`), сохранённый
для backwards compat с .spe-headers, где LSRM-conv транслитерирует Cyrillic-«С»
по звуку → latin-«S». **Для нашей станции authoritative canonical = `Gamma-1C`**.
BUG-40 (`KNOWN_AND_FIXED_ISSUES.md:1401-1422`) — defensive warning, когда
canonicalizer ловит cyrillic-latin homoglyph ambiguity на .spe header.

**Следствие методологии**: visual-templates для NaI-63×63 crystal-class
переносимы на любую future NaI-63×63 station (не только нашу); в каждой
template-записи провенанс указывает station-instance (для downstream
efficiency/FWHM-attribution). Сейчас все 24 канонических templates observed-on
= `Gamma-1C`.

---

## 2. Shielding configuration

- **Outer shield:** **Pb 50 mm** — primary attenuation of cosmic background and external γ
- **Inner liner (graded shield):** **Cd + Cu cup** — suppresses Pb K-X-rays (73–87 keV)
  generated in the Pb shield by cosmic-ray excitation, cascading through Cd
  (K-XR 23 keV) → Cu (K-XR 8 keV)
- **Observable consequence:** without the Cd+Cu liner, the Pb K-XR peak in BG
  dominates the 73–87 keV region (~1400 counts in maximum vs ~700 continuum in
  a side-by-side comparison). With the liner, residual Pb K-XR in BG is
  ~0.1–0.3 cps — at the level of natural background continuum.
- **Cd Kα 23 keV and Cu Kα 8 keV are not seen in the BG spectrum** — both sit
  below the ADC threshold of this NaI 63×63 + USB front-end configuration
  (effective cutoff ~10–15 keV). Absence is not evidence against the liner.
- **For sample analysis:** the Pb K-XR cluster (73–87 keV) seen in any sample
  spectrum is dominated by **internal conversion X-rays from the source decays**
  (Pb-212 238 keV IC, Bi-212 / Tl-208 transitions IC) plus direct chain γ
  (Th-228 84.37 keV) — *not* by Pb fluorescence of the shielding.

The graded liner is what makes Gamma-1C a **low-background** instrument suitable
for trace-level identification at LSRM/ISO 11929 sensitivities.

---

## 3. Vessel canonicalization (ЛСРМ source-of-truth, стр. 11)

**Vessel canonicalization rule** (ЛСРМ source-of-truth «Прецизионные
измерения. Образцовые и калибровочные источники.», стр. 11 +
user-locks 2026-06-05):

### Constructive vs operational — три отдельных поля, не один tag

Ранее ошибочно conflated в один `geometry` string. Правильное разделение:

| Поле | Что | Источник | Per-detector? |
|---|---|---|---|
| `vessel_class` | Конструктивная спецификация сосуда (capacity, габариты, d_эфф) | ЛСРМ-source-of-truth | НЕТ — переносим между детекторами одного crystal-class |
| `effective_thickness_mm` | d_эфф для self-attenuation correction | ЛСРМ-source-of-truth | НЕТ |
| `useful_sample_volume_ml` | Реально загруженный объём пробы для данной измерительной геометрии | Operator / passport / .spe header | **ДА** — может различаться |
| `placement_distance_cm` | Sample-to-detector дистанция | Operator / passport | **ДА** |

### Canonical vessel-classes (по ЛСРМ-source стр. 11)

| `vessel_class` | Capacity | Габариты (мм) | `effective_thickness_mm` |
|---|---|---|---|
| `marinelli_0.5L` | 0.5 L | ⌀125, H=100 | `[15, 2]` |
| `marinelli_1L`   | 1.0 L | ⌀150, H=110 | `[26, 2]` |
| `marinelli_3L`   | 3.0 L | ⌀180, H=200 | `[60, 5]` |
| `denta_120ml`    | 0.12 L | ⌀75, H=35 | `[36, 2]` |
| `petri_75ml`     | 0.075 L | ⌀88, H=14 | `[15, 2]` |
| `point_source`   | — | non-table | — |
| `other_<spec>`   | non-ЛСРМ | passport | passport |

### Operator-data canonicalization (этот FS layout → vessel tags)

| Folder name | `vessel_class` | `useful_sample_volume_ml` | Detector scope |
|---|---|---|---|
| `Маринелли` (89 records) | `marinelli_1L` (lab convention: default = 1L) | по умолчанию = vessel capacity 1000 мл (fill до ring mark) | universal |
| `Маринелли 1л` (37 records) | `marinelli_1L` | 1000 мл | universal |
| `MARINELLI` (11) / `Marinelli` (1) | `marinelli_1L` (presumed, verify per .spe) | 1000 мл (presumed) | universal |
| `marinelli_0cm` | `marinelli_1L` (suffix = placement, не volume; +`placement_distance_cm: 0`) | 1000 мл | universal |
| `Дента-120мл` | `denta_120ml` | 120 мл | universal |
| `Дента-120` (typo) | `denta_120ml` | 120 мл | universal |
| `Дента-100` (3 records, Поверка-2016) | **PENDING source-reconciliation** (см. ниже) | — | — |
| `Петри-60мл` | `petri_75ml` (vessel = 75 мл ЛСРМ) | **60 мл (useful)** — **ТОЛЬКО Gamma-1C** | Gamma-1C **only** |

### `Дента-100` — pending reconciliation

User lock 2026-06-05 = «реальная отдельная геометрия 100 мл», но
ЛСРМ-source стр. 11 не упоминает 100 мл vessel. Три кандидата
интерпретации:
- (a) non-ЛСРМ vessel 100 мл → `vessel_class: "other_denta_100ml"`, флаг `lsrm_standard: false`;
- (b) ЛСРМ-Дента 120 мл с partial-fill → `vessel_class: "denta_120ml"` + `useful_sample_volume_ml: 100`;
- (c) operator typo. Решение отложено до explicit user-instruction.

### Petri 60 мл = useful, не vessel — **ТОЛЬКО для Gamma-1C** (user lock 2026-06-05)

> «Это общий объем, полезный — 60 мл его и используем для Гамма-1С.
>  Для других типов детекторов надо смотреть индивидуально.» (user)

⇒ Для других crystal-class (HPGe-coaxial-20pct, LaBr3, CdZnTe, Si(Li),
Si-surface-barrier) с тем же `vessel_class: "petri_75ml"` —
`useful_sample_volume_ml` распаковывается из passport / .spe header
конкретной записи, **не наследуется** от Gamma-1C default 60 мл.

### Архитектурное правило

Visual-template shape определяется кортежем (`crystal_class` ×
`vessel_class` × `useful_sample_volume_ml`). Constructive part
(vessel + d_эфф) переносим между станциями того же crystal-class.
Operational part (useful_volume + placement_distance) — per-station.

---

## 4. Runtime conflict-resolution rule (user lock 2026-06-05)

> «При конкретных замерах параметры геометрий могут изменяться
>  оператором. Поэтому при анализе спектра в случае разночтений,
>  ты должен уточнять.» (user)

**Geometry parameters at measurement time may diverge from any nominal
spec** — оператор может загрузить partial-fill, использовать non-LSRM
vessel substitute, переместить пробу на нестандартную дистанцию, и
т.д. Folder name / .spe header / passport / operator metadata —
**четыре независимых provenance layer'a**, и в общем случае они могут
расходиться.

**Precedence для logging** (не для silent decision):

| Rank | Layer | Authority |
|---|---|---|
| 1 | `operator_explicit_metadata` | Оператор явно указал в session metadata |
| 2 | `passport_pdf_block` | Паспорт источника / certificate |
| 3 | `spe_header_sample_geometry` / `spe_header_sample_volume` | .spe header field |
| 4 | `folder_name_hint` | Имя папки (легко переименовывается, не authoritative) |
| 5 | `lsrm_source_default` | ЛСРМ-source constructive default (только для vessel_class fallback, **никогда** для operational useful_volume) |

**Conflict handling**:

- **Offline RAG-build pipeline (W4)**: при detected разночтении в
  vessel/volume/distance между layers — log в `__geometry_conflicts`,
  set `__needs_operator_review: true`, template уходит в
  `_pending_review/` (НЕ в canonical pool, **не** eligible для
  similarity API).
- **Online analyzer (production runtime)**: при analyze конкретного
  спектра, если detected разночтение между provenance layers —
  analyzer **ОБЯЗАН prompt оператора**: «Folder name suggests X,
  but .spe header reports Y. Какие фактические параметры замера?»
  Efficiency calculation **gated** до явного operator-decision.
  **НЕ resolve по precedence silent** — precedence rank это лишь
  для structured logging.

Это применимо ко всем трём operational полям (`vessel_class`
substitution / `useful_sample_volume_ml` / `placement_distance_cm`)
и развивает anti-hallucination правило проекта CLAUDE.md «каждое
утверждение ссылается на конкретный offset/строку/таблицу в
исходнике» — расширяя его на multi-layer conflict case.

---

## 5. Layout

```
detectors/Gamma-1C/
├── README.md                          # this file
├── certificates/                      # passports of standard sources (.xls / .pdf / .src)
├── data/                              # secondary peaks catalogs + aliases overrides
│   ├── averaged_backgrounds/          # 5 averaged .spe background files per geometry
│   ├── secondary_peaks.json           # Cs-137 + K-40 secondary peak catalog
│   └── secondary_peaks_v2.json        # 9-isotope rich catalog (incl. chain proxies)
├── efficiency/                        # .efr efficiency curves per geometry
│   └── Gamma-1C_NaI_63x63_USB_SN-01/
├── lsrm-libraries/                    # LSRM SpectraLine nuclide libraries
├── reference_spectra/                 # .spe reference spectra (verification campaigns)
│   └── Gamma-1C_NaI_63x63_USB_SN-01/
├── references/
│   ├── 05_intrinsic_detector_activity.md   # NaI(Tl) 63×63-specific intrinsic signatures
│   └── 07_dead_time_correction.md          # A, B coefficients for the УДС-ГЦ
└── raw_lsrm/                          # LOCAL-ONLY operator LSRM-tree, gitignored (F-115)
    ├── Work/
    │   ├── BG/Gamma-1C/Spe/           ← рабочий ствол LSRM-станции
    │   │   ├── Маринелли/             ┐ обе подпапки → один канонический
    │   │   ├── Маринелли 1л/          ┘ tag `marinelli_1L` (default = 1 L,
    │   │   │                            user lock 2026-06-05)
    │   │   ├── Точечная-25см/
    │   │   ├── Background/
    │   │   └── Spe — поверки/Поверка YYYY/
    │   ├── Calibration/
    │   └── …
    └── passports/                     ← (опционально) .pdf/.txt паспорта источников
```

---

## 6. Path resolver

Python code accesses these assets via the resolver module:

```python
from gamma.detectors.gamma1c import (
    DETECTOR_ROOT,
    CERTIFICATES_DIR,
    EFFICIENCY_DIR,
    REFERENCE_SPECTRA_DIR,
    LSRM_LIBRARIES_DIR,
    AVERAGED_BACKGROUNDS_DIR,
    SECONDARY_PEAKS_PATH,
    SECONDARY_PEAKS_V2_PATH,
    DEFAULT_REFERENCE_DIR,
    DEFAULT_EFFICIENCY_DIR,
    DETECTOR_NAME,
)
```

Never hardcode `detectors/Gamma-1C/...` paths in calling code — always go through the
resolver so future detectors can swap to their own subtrees without invasive churn.

---

## 7. Local-only `raw_lsrm/` working copy (НЕ коммитится)

Полная LSRM-tree оператора кладётся в подпапку `raw_lsrm/`, и в git
**никогда не попадает** (паттерн F-150/F-293 «books_library»).

### Exclusion guarantees (defence-in-depth)

`raw_lsrm/` исключается из артефактов на **трёх** уровнях:

1. **`.gitignore`** — паттерн `detectors/Gamma-1C/raw_lsrm/` →
   git track никогда не подхватит .spe / паспорта / сертификаты оператора.
2. **`scripts/build_release_archive.py:EXCLUDE_DIRS`** — basename
   `raw_lsrm` исключён → даже если случайно файл попадёт в working tree
   между gitignore-passes, release-archive его пропустит.
3. **F-115 anonymizer** (`scripts/gamma/reporting/anonymize.py`) — любой
   output-артефакт (JSON / Markdown / HTML / PDF), который сошлётся
   на путь внутри `raw_lsrm/`, будет скраббить absolute path до basename
   до записи на диск.

---

## 8. Isolation policy (v1.12.0)

This folder is **isolated** from any other detector. Algorithms in `scripts/gamma/`
are shared, but data, certificates, .efr curves, intrinsic-activity references and
secondary-peak catalogues here are valid **only** for the Gamma-1C complex.

When the AtomSpectra / AtomNano / RadiaCode pipelines are added (deferred), each
gets its own `detectors/<canonical>/` folder. Scripts will be copied across only
after they are stabilized in the Gamma-1C branch (per user policy 2026-05-29).

---

## 9. Crystal-class map (зафиксировано на 2026-06-05)

Что точно известно из проектных source-комментариев — для future stations:

| Station-instance | Crystal-class (точно) | Source-pin |
|---|---|---|
| `Gamma-1C` (the only one in this project) | NaI 63×63 | `audit/_rag/visual_templates/SCHEMA.md:10`, `scripts/gamma/io/lsrm_spe.py:29` |
| `GP_HPGe20` | HPGe coaxial 20% | имя scope-glob `Work\GP\HPGe(20%)` (`build_spectra_index.py:60`) |
| `NM_HPGe20` | HPGe coaxial 20% | имя scope-glob `Work\NM\HPGe(20%)` (`build_spectra_index.py:66`) |
| `Handy_NaI` | NaI (размер TBD) | имя scope-glob `Work\Handy\Handy(NaI)` |
| `Handy_HPGe` | HPGe (размер TBD) | имя scope-glob `Work\Handy\Handy(HPGe)` |
| `Handy_LaBr` | LaBr3 (размер TBD) | имя scope-glob `Work\Handy\Handy(LaBr)` |
| `Simple_NaI` | NaI (Demo, размер TBD) | имя scope-glob `Work\Simple\NaI(Demo)` |
| `Simple_HPGe` | HPGe (Demo, размер TBD) | имя scope-glob `Work\Simple\HPGe(Demo)` |
| `Simple_TeCd` | CdZnTe / CZT (размер TBD) | имя scope-glob `Work\Simple\TeCd(Demo)` |
| `Simple_SiLi` | Si(Li) (размер TBD) | имя scope-glob `Work\Simple\SiLi(Demo)` |
| `Simple_Alpha` | Si surface barrier (alpha) | имя scope-glob `Work\Simple\Alpha(Demo)` |

«размер TBD» — выяснить у оператора (или из паспорта DETECTOR= поля
конкретных .spe) ДО W5+ харнесса; в build-скрипт добавить
`station → crystal_class` lookup-table.

**Note**: «Gamma-1S» — НЕ separate station. Если оператор привёз новую
LSRM-выгрузку с DETECTOR-header «Гамма-1С №NNNN-NN» — это **тот же**
Gamma-1C в другой Поверка-эпохе. Кладётся в `raw_lsrm/Work/BG/Gamma-1C/`
(см. §7); cross-epoch / cross-Поверка drift study обрабатывается через
`_drift_study/` mirror (см. SCHEMA.md «drift-study isolation»).

---

## 10. F-070 W4 Use-case (extended Gamma-1C ingest)

См. `audit/_plans/F-070_W4_gamma1c_visual_templates_TODO.md` (renamed from
`_gamma1s_` after 2026-06-05 user lock).

Краткий контракт:
1. Оператор копирует LSRM-tree в `detectors/Gamma-1C/raw_lsrm/Work/...`
2. Build-скрипт `scripts/rag/build_visual_templates_nai63x63.py` (TODO)
   читает `raw_lsrm/...` через **относительные пути**.
3. Каждый emitted VT-*.json:
   - провенанс basename only (F-115 anonymizer)
   - `detector_id` → `УДС-ГЦ-63×63-USB` (S/N `№NNNN-NN` срезан)
   - `sample_id` cert-S/N паттерны (`420-7-XX`) → `None`
   - `crystal_class: "NaI-63x63"` поле + `station_observed_on: "Gamma-1C"`
     поле (новые в schema 0.2 — добавляются в S0 retrofit).
4. Результат — JSON-templates в
   `audit/_rag/visual_templates/<class>/VT-<NUC>-<GEOM>-<EPOCH>.json`
   (с примечанием station-instance в провенансе).

---

## 11. Cross-refs

- F-78 / F-78a (aliases): `data/aliases.json`, `scripts/gamma/data/aliases.py`
- F-83 (detector isolation): this README §8
- F-115 (анонимизация): `scripts/gamma/reporting/anonymize.py`
- F-150 / F-293 (паттерн external-data working-copy): `books_library/` precedent
- F-155 (allow-list корневых папок включает `detectors/`)
- BUG-40 (cyrillic-latin homoglyph warning): `KNOWN_AND_FIXED_ISSUES.md:1401-1422`,
  `tests/step04_detector_type/test_bug40_cyrillic_latin_warning.py`
- Crystal-vs-station distinction (LOCK 2026-06-05): this file §1
- W1 / W2 / W3 visual templates harness for NaI-63×63 / Gamma-1C station:
  `audit/_rag/visual_templates/SCHEMA.md`, `SIMILARITY_POLICY.md`
- SPECTRA_INDEX Gamma-1S метадата-записи (~394, исторический artifact
  v1.11.1 aliases): `audit/_rag/SPECTRA_INDEX.json` (filter `by_detector.Gamma-1S`)
  — это та же физическая Gamma-1C, transliterated by LSRM-header parser
  via cyrillic-latin omoglyph mapping; см. §1 «Canonical name fixation»

---

## 12. Last release that touched this folder

- `v1.12.0` — initial isolation of Gamma-1C-specific assets (F-83).
- `v1.25.0` (in progress) — F-070 W4 schema 0.2 retrofit (crystal-class
  abstraction + vessel taxonomy + runtime conflict-resolution rule).
