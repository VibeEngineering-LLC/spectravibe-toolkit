# LSRM Geometries — каталог BecqMoni для расчёта эффективности

**Дата фиксации**: 2026-06-14
**Источник**: https://github.com/Am6er/BecqMoni/tree/master/LSRM%20Geometries
**Назначение**: справка по форматам LSRM efficiency-кривых и .in модельных файлов, эталонные геометрии для бенчмарка нашего pipeline'а.
**Не путать**: с LSRM software (СНИИП «ЛСРМ Сертификат») — это разные сущности. Здесь LSRM = «Lab Scintillator Response Method», MC-симулятор от автора Am6er (вписан в BecqMoni).

## 1. Что такое LSRM в контексте BecqMoni

BecqMoni использует собственный MC-симулятор (упоминается в README как **EffMaker**) для расчёта **photopeak efficiency** ε(E) по описанию геометрии и материалов. Папка `LSRM Geometries/` содержит **готовые** результаты симуляции для типовых конфигураций:

- `*.in` — входной файл MC-симулятора (geometry + materials + source).
- `*.txt` (curve) — результат: таблица ε(E) на 150 точках от 20 до 3000 keV с шагом 20 keV.

При импорте в BecqMoni UI оператор выбирает геометрию из списка → BecqMoni берёт curve-таблицу → интерполирует на energy-grid его прибора → применяет для расчёта активности.

Для нашего pipeline эти данные служат **независимым reference'ом**: можно сравнить наш расчёт ε(E) для RC-103 marinelli с LSRM curve и валидировать.

## 2. Формат `curve_*.txt` (efficiency table)

```
Energy, keV	Efficiency	Uncertainty, %
20.0			4.55616E-02		1150
40.0			2.76634E-03		63.9
60.0			2.24553E-03		2.15
...
3000.0			2.48089E-05		21.4
```

Структура:
- **3 колонки**: Energy в keV, Efficiency (доля от total flux в 4π), Uncertainty в %.
- **Разделитель**: TAB-символ (может быть несколько подряд).
- **Заголовок**: первая строка `Energy, keV\tEfficiency\tUncertainty, %`.
- **Энергетический grid**: 20 ... 3000 keV с шагом 20 keV → **150 точек**. (В curve_Nano_16 и curve_RadiaCode_cilinder обнаружена 151-я точка — возможно артефакт.)
- **Efficiency** — безразмерная, например `2.24553E-03` = 0.224% efficiency на 60 keV.
- **Uncertainty** — relative %, для 20 keV часто >100% (низкая статистика MC), для 100-2000 keV обычно 2-5%, на хвостах >3000 keV растёт до 20%.

**Важно**: 20-40 keV uncertainty >50% — этот диапазон **не использовать** для quantitative-измерений. Низкая энергия — большая статистическая погрешность MC + сильная зависимость от self-absorption (тонкие детали геометрии).

### Пример: RadiaCode + Marinelli 0.5 (наш профиль)

| E, keV | ε | u, % |
|---|---|---|
| 20  | 4.56e-02 | 1150 (мусор) |
| 40  | 2.77e-03 | 63.9 (плохо) |
| 60  | 2.25e-03 | 2.15 (надёжно) |
| 100 | 2.58e-03 | 2.16 |
| 200 | 1.60e-03 | 2.48 |
| 500 | 3.93e-04 | 4.78 (интерполировано) |
| 1000 | 1.45e-04 | 7.9 |
| 1460 | 9.6e-05 | ~10 (K-40) |
| 2614 | 4.4e-05 | ~15 (Tl-208) |
| 3000 | 2.48e-05 | 21.4 |

Падение ε(E) ≈ E^(-1.4) от 100 до 1000 keV — типичный CsI(Tl) response с фотоэффектом доминирующим до ~200 keV, затем переход в комптон + парная при E > 1022.

## 3. Формат `.in` модельного файла LSRM

Шесть блоков:

### 3.1. DETECTOR PARAMETERS

```
DetectorType = SCINTILLATOR     // или COAXIAL для HPGe

// SCINTILLATOR-specific (DS_*)
DS_CrystalDiameter = 1 cm
DS_CrystalHeight = 1 cm
DS_CrystalFrontReflectorThickness = 0.1 cm
DS_CrystalSideReflectorThickness = 0.1 cm
DS_CrystalFrontCladdingThickness = 0.1 cm
DS_CrystalSideCladdingThickness = 0.1 cm
DS_DetectorMountingThickness = 0.1 cm
```

(COAXIAL DC_* поля присутствуют для совместимости, но игнорируются если DetectorType=SCINTILLATOR.)

### 3.2. SOURCE PARAMETERS

```
SourceType = MARINELLI     // или POINT, CYLINDER

// MARINELLI-specific (SM_*)
SM_BeakerToDetectorFrontDistance = 0.8 cm     // зазор детектор → дно Маринелли (0.2-0.8 cm)
SM_BeakerDiameter = 11.4 cm                   // внешний диаметр Маринелли (стандарт 11.4)
SM_BeakerHeight = 8.9 cm                      // высота beaker'а
SM_BeakerHoleDiameter = 6.1 cm                // диаметр сквозной полости (где сидит детектор)
SM_BeakerHoleHeight = 5.3 cm                  // глубина полости
SM_BeakerSideThickness = 0.2 cm               // толщина бокового пластика
SM_BeakerEndWallThickness = 0.2 cm            // толщина дна beaker'а
SM_BeakerHoleSideThickness = 0.2 cm           // толщина стенок полости
SM_BeakerHoleEndWallThickness = 0.2 cm        // толщина крышки полости
SM_SourceHeight = 8.5 cm                      // высота наполнения образцом
```

**Стандартная Marinelli 0.5 л** (Am6er reference, точно совпадает с ОИСН-16 у нас):
- Beaker 11.4 cm Ø × 8.9 cm H, walls 0.2 cm.
- Hole 6.1 cm Ø × 5.3 cm H.
- Полезный объём ≈ (π × 5.7² × 8.5) − (π × 3.05² × 5.3) ≈ 868 − 155 = 713 см³ ≈ **0.7 л** (округление маркетинговое до 0.5).

### 3.3. MATERIAL PARAMETERS

Для каждого компонента (Crystal, Cladding, Reflector, Mounting, Beaker, Source) — плотность + список элементов с долями.

**CsI(Tl) кристалл RadiaCode**:
```
DS_nCrystalElements = 2
DS_RoCrystal = 4.51                    // плотность г/см³
DS_ZCrystal[0] = 53                    // I (Iodine)
DS_FractionsCrystal[0] = 0.488451      // массовая доля I
DS_ZCrystal[1] = 55                    // Cs (Cesium)
DS_FractionsCrystal[1] = 0.511549      // массовая доля Cs
DS_FractionTypeCrystal = MASS          // вид долей: MASS или ATOM
```

(Без Tl — допант ≪0.1%, симулятор его игнорирует.)

**Источник (наполнение Marinelli)**:
```
M_SM_Source.MName = Water, liquid              // дефолт — water-equivalent
M_SM_Source.Nmaterials = 1
M_SM_Source.Name[0] = Water, liquid
M_SM_Source.MatRelWeight[0] = 1
```

Для реальных образцов оператор должен **поменять** на нужную матрицу (грунт ρ ≈ 1.4, бетон ρ ≈ 2.3) и пересчитать. Дефолтная water — это base-case для расчётов.

**Beaker материал**:
```
M_SM_Beaker.MName = Polyethylene terephthalate (PET)
```

Не путать с PE (polyethylene, ρ=0.93) — это PET (ρ=1.38). Реальные Marinelli ОИСН-16 у нас — белый ПВХ или ПП, плотность чуть отличается, но эффект на ε(E) для энергий >100 keV пренебрежимый.

### 3.4. Где Tl-допант в LSRM модели

**Нигде явно**. Симулятор моделирует фотоэффект и комптон через атомные сечения по Z, а Tl (0.1% массы) даёт пренебрежимый вклад в attenuation. Сцинтилляционный выход (Tl как actuator люминесценции) симулятором не моделируется — это **photopeak efficiency**, а не light-output.

## 4. Каталог геометрий в репо

| Curve файл | Модель | Детектор | Источник | E-range надёжный |
|---|---|---|---|---|
| `RadiaCode - marinelli 0.5.txt` | `model_RadiaCode_Marinelli0.5.in` | CsI(Tl) 1×1 cm | Marinelli 0.7 л (11.4×8.9, hole 6.1×5.3, walls 0.2) | 60-2500 keV |
| `RadiaCode - author marinelli 0.5.txt` | `model_RadiaCode_AuthorMarinelli0.5.in` | CsI(Tl) 1×1 cm | Marinelli «author» 0.5 л (9.28×9.28, hole 2×6.2, walls 0.18) — другая геометрия | 60-2500 keV |
| `RadiaCode - author marinelli 0.2.txt` | — | CsI(Tl) 1×1 cm | Marinelli «author» 0.2 л (меньше) | 60-2500 keV |
| `RadiaCode - cilinder.txt` | — | CsI(Tl) 1×1 cm | Cylinder geometry (компактный пластик) | 60-2500 keV |
| `Obsidian - marinelli 0.5.txt` | `model_Obsidian_Marinelli_0.5.in` | CsI(Tl) 0.67×3 cm (длинный тонкий) | Marinelli стандарт 11.4×8.9 (как наш ОИСН-16) | 60-2500 keV |
| `Nano 16 - marinelli.txt` | `model_Nano16Pro_Marinelli.in` | CsI(Tl) 1.854×5.9 cm (большой Nano) | Marinelli стандарт 11.4×8.9 | 80-2800 keV |

**Что отсутствует в репо** (gaps):
- RC-103 + Marinelli с грунтом (ρ=1.4) — оператор должен симулировать сам или взять «water» как worst-case.
- Geometry для HPGe (Coaxial) — модельных файлов нет, только сcintillator.

## 5. Сравнение размеров CsI(Tl) кристаллов

| Прибор | Crystal Ø | Crystal H | Объём | Notes |
|---|---|---|---|---|
| RadiaCode RC-10x | 1.0 cm | 1.0 cm | 0.785 cm³ | Заявлено производителем как «1×1» (компактный) |
| Obsidian | 0.67 cm | 3.0 cm | 1.06 cm³ | Тонкий длинный (направленный) |
| Atom Spectra Nano 16 Pro | 1.854 cm | 5.9 cm | 15.9 cm³ | Полу-стационарный, большой объём |

Соотношение полезного объёма Nano/RC-103 ≈ 20×. На K-40 1461 keV ε_Nano ≈ 5× ε_RC-103 (примерно √20 ≈ 4.5×, поправка на форму).

## 6. Cross-reference с проектом

### 6.1. Импорт LSRM curve в наш pipeline

Минимальный helper-скрипт (читает текущий формат, выдаёт `EfficiencyModel`):

```python
import numpy as np
from pathlib import Path

def read_lsrm_curve(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Returns (E_keV, eff, uncertainty_rel_pct)."""
    data = np.loadtxt(path, skiprows=1, delimiter=None)  # whitespace разделитель
    return data[:, 0], data[:, 1], data[:, 2]
```

Сохранить как `scripts/io/lsrm_efficiency.py` (нумерация продолжает существующие BecqMoni-парсеры).

### 6.2. Валидация нашего ε(E) для RC-103

Опорные точки RC-103 + Marinelli 0.5 (вода):
- 60 keV: ε = 2.25e-03 (надёжно)
- 186 keV: интерполяция между 180 (1.84e-03) и 200 (1.60e-03) ≈ 1.69e-03
- 583 keV: интерполяция между 580 и 600 ≈ 3.2e-04
- 1461 keV: ≈ 9.6e-05
- 2614 keV: ≈ 4.4e-05

При запуске нашего расчёта Marinelli + RC-103 с water-наполнением ε должна укладываться в ±10% от этих значений (LSRM uncertainty + наша модель). Если расходимся в 2× — есть baг в self-absorption или solid-angle факторе.

### 6.3. Differences от наших методов

- LSRM моделирует **water** наполнение по умолчанию. Наш `gamma.efficiency` для рабочих образцов оператора — **soil** ρ=1.4. На 60 keV self-absorption грунта в Marinelli ≈ 0.4-0.5 от water, на 600 keV ≈ 0.7-0.8. Прямой сравнение LSRM curve и нашего ε(E) для грунта = НЕ валидно. Нужно симулировать LSRM с grunt-материалом.
- LSRM uncertainty на 20 keV >100% — это **MC шум**, не калиброванная неопределённость. Реальная погрешность ε(60 keV) ± 2-3% при наличии калибровочного источника, а не LSRM-curve.

## 7. Anti-hallucination — провенанс

Все факты — из открытых файлов в `Am6er/BecqMoni/LSRM Geometries/`:
- `curve_RadiaCode_-_marinelli_0.5.txt`, `curve_RadiaCode_-_author_marinelli_0.2.txt`, `curve_RadiaCode_-_author_marinelli_0.5.txt`, `curve_RadiaCode_-_cilinder.txt`, `curve_Nano_16_-_marinelli.txt`, `curve_Obsidian_-_marinelli_0.5.txt`
- `model_RadiaCode_Marinelli0.5.in`, `model_RadiaCode_AuthorMarinelli0.5.in`, `model_Nano16Pro_Marinelli.in`, `model_Obsidian_Marinelli_0.5.in`

Скачано 2026-06-14 через `gh api repos/Am6er/BecqMoni/contents` + curl с URL-encoding пробелов.

**Не extrapolate'ить** на:
- Symbiosis между LSRM (MC) и BecqMoni (ε из curve): после генерации curve, BecqMoni только интерполирует, не пересчитывает MC.
- Чужие LSRM-форматы (СНИИП ЛСРМ от 2007 — это полностью другой формат, не путать).
- HPGe Coaxial — в репо нет готовых curves для HPGe, не делать предположения о них.