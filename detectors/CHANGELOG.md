# Журнал изменений по спектрам детекторов

Подробная хроника того, что происходило с файлами в `detectors/`. Общая сводка
по количеству спектров — в [`SPECTRA_MAP.md`](../SPECTRA_MAP.md), она
генерируется по индексу git и здесь не дублируется.

---

## 2026-07-27 — вшитый фон 5-й степени портил энергетическую шкалу

**Затронуто:** 51 XML в `GP_HPGe20` (29) и `NM_HPGe20` (22).
**Статус:** исправлено, фон из этих файлов вынут.

### Что случилось

BecqMoni проверяет `<EnergySpectrum>` и `<BackgroundEnergySpectrum>` как две
**независимые** калибровки. `PolynomialEnergyCalibration.CheckCalibration()`
отвергает любой порядок выше четвёртого:

```csharp
if (this.polynomialOrder > 4 || this.polynomialOrder < 1) { return false; }
```

`DocumentManager.CheckDocument()` пробует спасти такую калибровку понижением
степени, но понижает её только за счёт нулевых старших коэффициентов. У
настоящего полинома 5-й степени нулей нет, срабатывает ветка `zerosCount == 0`,
и калибровка заменяется на **дефолтную**. Файл открывается, ошибки нет, шкала
потеряна.

Проверка фона идёт отдельной веткой (`DocumentManager.cs:204`), поэтому фон со
степенью 5 губил документ **даже когда сам спектр был линейным**. Именно так и
оказалось: у 47 из 51 файла собственная калибровка была 1-й или 2-й степени —
негодным их делал исключительно подшитый фон.

### Виновник

Один файл на класс:

| Класс | Фон | Степень | Скольким файлам подшит |
|---|---|---:|---:|
| `GP_HPGe20` | `Bckg_1.spe` | 5 | 29 |
| `NM_HPGe20` | `Bckg_1.spe` | 5 | 22 |

`Bckg_15.spe` в обоих классах тоже несёт 5-ю степень, но ни одному измерению
подшит не был.

### Как обнаружено

Полный скан 835 XML под `detectors/` с разбором степени полинома в обоих
блоках. До правки:

| | файлов |
|---|---:|
| степень >4 где-либо | 56 |
| — виноват только фон | 47 |
| — виноват только спектр | 5 |
| — оба | 4 |

Прежний счётчик в `SPECTRA_MAP.md` показывал 9: он читал флаг
`becqmoni_reads_as_linear` из `INDEX.json`, а тот ставился только по калибровке
самого спектра. Фон не разбирался вовсе.

### Что сделано

1. `<BackgroundEnergySpectrum>` вынут из 51 XML. Ссылка
   `<BackgroundSpectrumFile>` сохранена — видно, какой фон относится к
   измерению.
2. Оба конвертера больше не вшивают фон со степенью выше 4:
   `scripts/convert/lsrm_tree_publish.py`, `scripts/convert/reference_kits_to_becqmoni.py`.
3. В описи добавлены поля `background_cal_degree`, `background_embedded`,
   `background_dropped_high_order`.
4. `scripts/convert/build_spectra_map.py` считает степень по самим XML, а не по
   флагу в описи, и разбирает оба блока.

Коэффициенты калибровки нигде не пересчитывались — сохранённые данные не
трогали.

### Осталось как есть

9 спектров, у которых **сам** спектр записан с полиномом 5-й степени. Они
перечислены в `SPECTRA_MAP.md`; коэффициенты лежат в `INDEX.json` полем
`energy_cal`. Пересчёт к 4-й степени был реализован и отклонён: на HPGe он
ложится в пределах 0,21 кэВ, но на изогнутой NaI-калибровке промахивается на
57 кэВ.

### Полный перечень

**`GP_HPGe20` — 29 файлов**, пути от `reference_spectra/becqmoni/`:

```
Work/GP/HPGe(20_)/Spe/joint processing/Barrel_10cm_D1.xml
Work/GP/HPGe(20_)/Spe/joint processing/Barrel_10cm_D2.xml
Work/GP/HPGe(20_)/Spe/joint processing/Barrel_10cm_D3.xml
Work/GP/HPGe(20_)/Spe/joint processing/NZK_50cm_D1.xml
Work/GP/HPGe(20_)/Spe/joint processing/NZK_50cm_D2.xml
Work/GP/HPGe(20_)/Spe/joint processing/NZK_50cm_D3.xml
Work/GP/HPGe(20_)/Spe/joint processing/NZK_50cm_D4.xml
Work/GP/HPGe(20_)/Spe/Marinelli/m08085_mix09.xml
Work/GP/HPGe(20_)/Spe/Marinelli/m_cs06.xml
Work/GP/HPGe(20_)/Spe/Marinelli/m_cs16.xml
Work/GP/HPGe(20_)/Spe/Marinelli/m_k06.xml
Work/GP/HPGe(20_)/Spe/Marinelli/m_k16.xml
Work/GP/HPGe(20_)/Spe/Marinelli/m_ra06.xml
Work/GP/HPGe(20_)/Spe/Marinelli/m_ra16.xml
Work/GP/HPGe(20_)/Spe/Marinelli/m_th06.xml
Work/GP/HPGe(20_)/Spe/Marinelli/m_th16.xml
Work/GP/HPGe(20_)/Spe/Point25/Am241-SRC-01-25cm.xml
Work/GP/HPGe(20_)/Spe/Point25/Ba133-SRC-01-25cm.xml
Work/GP/HPGe(20_)/Spe/Point25/Cd109-SRC-01-25cm.xml
Work/GP/HPGe(20_)/Spe/Point25/Ce139-SRC-01-25cm.xml
Work/GP/HPGe(20_)/Spe/Point25/Co57-SRC-01-25cm.xml
Work/GP/HPGe(20_)/Spe/Point25/Co60-SRC-01-25cm.xml
Work/GP/HPGe(20_)/Spe/Point25/Cs137-SRC-01-25cm.xml
Work/GP/HPGe(20_)/Spe/Point25/Eu152-SRC-01-25cm.xml
Work/GP/HPGe(20_)/Spe/Point25/Mn54-SRC-01-25cm.xml
Work/GP/HPGe(20_)/Spe/Point25/Sn113-SRC-01-25cm.xml
Work/GP/HPGe(20_)/Spe/Point25/Th228-SRC-01-25cm.xml
Work/GP/HPGe(20_)/Spe/Point25/Y88-SRC-01-25cm.xml
Work/GP/HPGe(20_)/Spe/Th228--BadCalibre.xml
```

**`NM_HPGe20` — 22 файла:**

```
Work/NM/HPGe(20_)/Spe/Marinelli/m08085_mix09.xml
Work/NM/HPGe(20_)/Spe/Marinelli/m_cs06.xml
Work/NM/HPGe(20_)/Spe/Marinelli/m_cs16.xml
Work/NM/HPGe(20_)/Spe/Marinelli/m_k06.xml
Work/NM/HPGe(20_)/Spe/Marinelli/m_k16.xml
Work/NM/HPGe(20_)/Spe/Marinelli/m_ra06.xml
Work/NM/HPGe(20_)/Spe/Marinelli/m_ra16.xml
Work/NM/HPGe(20_)/Spe/Marinelli/m_th06.xml
Work/NM/HPGe(20_)/Spe/Marinelli/m_th16.xml
Work/NM/HPGe(20_)/Spe/Point25/Am241-SRC-01-25cm.xml
Work/NM/HPGe(20_)/Spe/Point25/Ba133-SRC-01-25cm.xml
Work/NM/HPGe(20_)/Spe/Point25/Cd109-SRC-01-25cm.xml
Work/NM/HPGe(20_)/Spe/Point25/Ce139-SRC-01-25cm.xml
Work/NM/HPGe(20_)/Spe/Point25/Co57-SRC-01-25cm.xml
Work/NM/HPGe(20_)/Spe/Point25/Co60-SRC-01-25cm.xml
Work/NM/HPGe(20_)/Spe/Point25/Cs137-SRC-01-25cm.xml
Work/NM/HPGe(20_)/Spe/Point25/Eu152-SRC-01-25cm.xml
Work/NM/HPGe(20_)/Spe/Point25/Mn54-SRC-01-25cm.xml
Work/NM/HPGe(20_)/Spe/Point25/Sn113-SRC-01-25cm.xml
Work/NM/HPGe(20_)/Spe/Point25/Th228-SRC-01-25cm.xml
Work/NM/HPGe(20_)/Spe/Point25/Y88-SRC-01-25cm.xml
Work/NM/HPGe(20_)/Spe/Th228--BadCalibre.xml
```

---

## 2026-07-27 — пары «измерение + фон» подбирались наугад

**Затронуто:** 87 измерений в четырёх классах.
**Статус:** данные оставлены как есть, пары помечены в описях как непроверенные.

### Что случилось

`pick_background()` брал фон из папки измерения, а если его там не было —
выбирал из всех фонов класса тот, у которого больше совпавших имён каталогов в
пути. Это не «ближайший по дереву», а пересечение множеств компонентов пути:
мера, ничего не говорящая об условиях измерения. При равенстве брался первый.

У классов, где фон вынесен в отдельный подкаталог, срабатывал именно запасной
вариант, и всем измерениям доставался один и тот же файл:

| Класс | Пар из папки измерения | Подобрано автоматически | Что доставалось |
|---|---:|---:|---|
| `GP_HPGe20` | 0 | 29 | `Bckg_1.spe` |
| `NM_HPGe20` | 0 | 26 | `Bckg_1.spe`, `Background_LEGe_…` |
| `Handy_LaBr` | 0 | 24 | `Background_1.spe` |
| `Simple_HPGe` | 0 | 8 | `Фон_2MeV.spe` |
| `Handy_HPGe` | 37 | 0 | — |
| `Handy_NaI` | 31 | 0 | — |

`Bckg_1.spe` — это **маринелли с дистиллированной водой**. Он оказался подшит к
точечным измерениям на 25 см, к бочке Barrel-SKG(Water) и к НЗК-150. Четыре
демо-спектра U-Coaxial получили фон планарного LEGe, то есть другого детектора.

### Признак верной пары — курирование, а не совпадение полей

Первым решением была сверка `GEOMETRY` спектра и фона. Правило негодное. В
`reference_kits/Denta_120mL/Cs-137/` фон записан как `GEOMETRY=Точечная-5см`
при образце `Дента-120мл`, а в его комментарии сказано прямо: снят под Дента,
Чашку-60 и Точечную-5см. Сверка полей забраковала бы верную пару.

Достоверна пара, которую сложил человек: фон, положенный рядом с образцом,
снят под эту задачу, и заголовок это подтверждать не обязан. Обратное тоже
верно — пара, которую никто не складывал, ничем не подтверждается.

### Что сделано

1. Запасное правило из `pick_background()` убрано. Фон берётся только из папки
   измерения; нет — значит пары нет. Ложная пара хуже отсутствующей: вычитание
   по ней считается молча.
2. В описи добавлено поле `background_curated`, в корень — список
   `backgrounds_available` с геометрией и примечанием каждого фона класса, и
   `background_pairing_note` там, где непроверенные пары есть.
3. **Существующие XML не тронуты.** Автоматически подобранный фон остался
   вшитым там, где был; курирование — отдельная работа.

⚠️ Следствие: повторный прогон `lsrm_tree_publish.py` по `GP_HPGe20`,
`NM_HPGe20`, `Handy_LaBr` или `Simple_HPGe` пересоберёт их XML **без фона** —
87 пар исчезнут. Это и есть намеренное поведение, но знать о нём надо заранее.

`Gamma-1S` устроен иначе и не затронут: его собирает
`reference_kits_to_becqmoni.py`, который разбирает каждую листовую папку как
готовую пару «образец + фон». Проверено — у всех 40 XML фон взят из своего же
листа.

---

## 2026-07-27 — фон вынесен в отдельную папку каждого детектора

Фоновые измерения больше не рассыпаны по дереву измерений. Для каждого
детектора они собраны в одну ветку:

```
detectors/<класс>/reference_spectra/background/lsrm/       исходные .spe
detectors/<класс>/reference_spectra/background/becqmoni/   XML
```

Исходная иерархия внутри ветки сохранена. Перенесено 60 файлов:

| Класс | `.spe` | XML |
|---|---:|---:|
| `Handy_LaBr` | 17 | 17 |
| `NM_HPGe20` | 6 | 4 |
| `GP_HPGe20` | 3 | 3 |
| `Handy_HPGe` | 2 | 2 |
| `Handy_NaI` | 2 | 2 |
| `Simple_HPGe` | 1 | 1 |

Пути в `INDEX.json` обновлены — записи фона теперь начинаются с `background/`.

**`Gamma-1S` не затронут ничем из перечисленного.** Класс перечислен в
`FROZEN_LAYOUT` (`scripts/convert/lsrm_tree_publish.py`): его фон остаётся там,
где лежит в дереве измерений, описи не переписываются, и повторный прогон
конвертера ничего не переставит. Курированные наборы `reference_kits/` тоже
нетронуты — там фон намеренно лежит в одной папке с образцом, парой на
геометрию.

---

## 2026-07-03 — единицы `SampleInfo` роняли BecqMoni

`DCSampleInfoView` трактует `Weight` как килограммы, `Volume` как литры, с
контролем диапазона 0,001…100. ЛСРМ пишет граммы и миллилитры, поэтому без
деления на 1000 приложение падало с `ArgumentOutOfRangeException` ещё до
отрисовки документа. Оба конвертера приводят единицы перед записью.

⚠️ Сама библиотека `scripts/gamma/io/becqmoni_xml.py` этой правки не содержит —
она пишет значение из заголовка как есть. Вызывать `write_becqmoni_xml()`
напрямую на данных ЛСРМ нельзя, только через конвертеры.
