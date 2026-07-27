# AtomSpectra / AtomNano test fixtures

> Перенесено сюда из `evals/fixtures/` в v1.18.7 (F-307) per F-157 detector
> isolation. Содержимое — реальные спектры от AtomSpectra Pro (10 файлов)
> и AtomNano (1 файл — `Радон.xml`), не от Gamma-1C.

## Использование в тестах

Эти файлы используются как fixtures в следующих тестах:

- `tests/io/test_format_conversion.py` — формат-конверсия
- `tests/io/test_reader_api.py` — embedded-bg parsing (Cs137_0_см.xml)
- `tests/step05_energy_calibration/test_stored_check_adaptive.py` — Алтайское_Зло (известный +48 keV сдвиг stored cal на 2614 keV)
- `tests/step05_energy_calibration/test_subcalibration.py` — Фон_типовой
- `tests/step07_identification/test_identification.py` — Фон_типовой
- `tests/step08_multiplets/test_peak_area.py` — Фон_типовой

При перемещении/переименовании файлов обновлять path в test-источниках
синхронно.
