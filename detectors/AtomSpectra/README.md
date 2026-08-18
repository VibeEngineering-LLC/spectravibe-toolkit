# AtomSpectra / AtomNano — detector subtree (F-307 isolation)

> **F-307 / v1.18.7 (закреплено 2026-05-31)** — изолированная папка-заготовка для
> AtomSpectra / AtomNano family детекторов per **F-157 isolation policy**
> (CLAUDE.md). Все ресурсы, связанные с этим детектором, живут здесь,
> отдельно от `detectors/Gamma-1C/`.

## Текущий статус

**Pre-Phase 1** — детектор-семейство **AtomSpectra Pro / AtomNano** распознаётся
проектом только на alias-уровне (через `gamma.io.atomspectra_xml` reader). Полная
методологическая ветка для этого детектора будет создана **после Phase 3 GA**
ветки Gamma-1C (см. CLAUDE.md §«Lifecycle phases»).

## Что лежит здесь сейчас

```
detectors/AtomSpectra/
├── README.md                       (этот файл)
└── data/
    └── fixtures/                   (тестовые xml-спектры, 11 шт.)
        ├── AtomSpectra Pro device:
        │   ├── Bq-2024-08-03_13-38-15-Жуковка_очищенная.xml
        │   ├── Cs137_0_см.xml + ..._-_subtract.xml
        │   ├── KCl__в_домике.xml + ..._-_subtract.xml
        │   ├── Алтайское_Зло_в_домике_маринелли_294_6г.xml
        │   ├── Дски_лаба_GS5050_8k.xml
        │   ├── Фон_Cs137_0_см.xml
        │   ├── Фон_KCl__в_домике.xml
        │   └── Фон_типовой_8192к_01-01-2025.xml
        └── AtomNano device:
            └── Радон.xml
```

## История перемещения

- **v1.18.7 / F-307 (2026-05-31)**: 11 xml-файлов перенесены из
  `evals/fixtures/` сюда. До переноса оставались в общей папке fixtures
  вместе с .spe файлами Gamma-1C, что нарушало F-157 (per-detector isolation).
  В `evals/fixtures/` остались только `M_cs/M_k/M_ra/M_th/*.spe`
  (Gamma-1C LSRM .spe формат).

## Roadmap для этого детектора (отложено)

| Релиз | Содержание | Готовность к старту |
|---|---|---|
| TBD | `detectors/AtomSpectra/efficiency/` — ε(E) калибровки | требует passport-данных от вендора |
| TBD | `gamma.detectors.atomspectra` — path resolver (аналогично `gamma.detectors.gamma1c`) | требует Phase 3 GA Gamma-1C |
| TBD | Алгоритмическая ветка: peak_image / identification / activity replicated from Gamma-1C tree | требует cert-fixtures для regression |

## Запреты

- НЕ удалять файлы из `data/fixtures/` (используются existing tests
  `tests/io/test_reader_api.py`, `tests/step05_energy_calibration/*.py`,
  `tests/step07_identification/test_identification.py`,
  `tests/step08_multiplets/test_peak_area.py`)
- НЕ перемещать обратно в `evals/fixtures/` (нарушает F-157)
- НЕ копировать алгоритмический код из `gamma.peaks/identification/calibration/...`
  в `gamma.detectors.atomspectra` до Phase 3 GA закрытия Gamma-1C ветки

См. CLAUDE.md §«Lifecycle phases» + SKILL.md §«Active spectrometric complex».
