# NM_HPGe20 — спектры LSRM и экспорт в BecqMoni

Спектров всего: **427**. Исходные `.spe` — в `lsrm/`, структура папок
сохранена. Экспорт в BecqMoni — в `becqmoni/`.

Из них 397 спектров опубликованы **только как `.spe`**: это объёмные
исследовательские коллекции (Pu/U LNHB, INR, IPPPE, kpti). У HPGe 16 384 канала,
поэтому один BecqMoni-XML весит около 590 КБ против 33 КБ исходника — держать их
в репозитории неоправданно. Собираются по требованию:

```bash
python scripts/convert/lsrm_tree_publish.py --class NM_HPGe20 --xml-all
```

Сконвертировано: 30 (26 образцов, 4 фоновых).

## Очистка персональных данных

Из заголовков удалено поле `OPERATOR` — оно было заполнено в 38 файлах.
Остальные поля не изменялись: `COMMENT` содержит только паспортные активности
источников, персональных данных в нём нет.

## Калибровка, которую BecqMoni прочитает неверно

`PolynomialEnergyCalibration.ChannelToEnergy()` разбирает только степени 2, 3 и 4;
всё, что выше, молча уходит в линейную ветку `c[1]·n + c[0]`. У перечисленных
файлов LSRM записал полином **5-й степени** — коэффициенты сохранены как в
источнике, но энергетическая шкала в BecqMoni будет неверной. Для расчётов
берите калибровку из `INDEX.json` (поле `energy_cal`).

- `Work/NM/HPGe(20_)/Spe/Background/Bckg_1.spe` — степень 5
- `Work/NM/HPGe(20_)/Spe/Background/Bckg_15.spe` — степень 5
- `Work/NM/HPGe(20_)/Spe/Marinelli/m08085_mix09.spe` — степень 5
- `Work/NM/HPGe(20_)/Spe/Point25/Y88-SRC-01-25cm.spe` — степень 5

## Перечень

| Файл | Источник | Паспорт | Геометрия | Детектор | Дата | Живое, с | Каналов | Фон |
|---|---|---|---|---|---|---:|---:|---|
| `Work/NM/HPGe(20_)/Spe/Marinelli/m08085_mix09.spe` | m08085_mix09 | Ti-44 A=3530 Bq /kg dA=6% 26-11-08 | Marinelli | GEM20P4-70 #SN-02 | 2008-12-19 | 3516 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Marinelli/m_cs06.spe` | m_cs06 | Cs-137 A=1760 Bk/kg dA=5% 24-05-02 | MARINELLI | GEM20P4-70 #SN-02 | 2008-12-19 | 1800 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Marinelli/m_cs16.spe` | m_cs16 | Cs-137 A=1650 Bk/kg dA=5% 17-09-07 | MARINELLI | GEM20P4-70 #SN-02 | 2008-12-17 | 1800 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Marinelli/m_k06.spe` | m_k06 | K-40 A=2530 Bk/kg dA=6% 24-05-02 | MARINELLI | GEM20P4-70 #SN-02 | 2008-12-19 | 3600 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Marinelli/m_k16.spe` | m_k16 | K-40 A=1950 Bk/kg dA=6% 17-09-07 | MARINELLI | GEM20P4-70 #SN-02 | 2008-12-17 | 3600 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Marinelli/m_ra06.spe` | m_ra06 | Ra-226 A=1780 Бк/кг dA=6% 24-05-02 | MARINELLI | GEM20P4-70 #SN-02 | 2008-12-19 | 3600 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Marinelli/m_ra16.spe` | m_ra16 | Ra-226 A=1790 Bk/kg dA=6% 17-09-07 | MARINELLI | GEM20P4-70 #SN-02 | 2008-12-17 | 1800 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Marinelli/m_th06.spe` | m_th06 | Th-232 A=2240 Bk/kg dA=6% 24-05-02 | MARINELLI | GEM20P4-70 #SN-02 | 2008-12-19 | 3585 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Marinelli/m_th16.spe` | m_th16 | Th-232 A=1940 Bk/kg dA=6% 17-09-07 | MARINELLI | GEM20P4-70 #SN-02 | 2008-12-17 | 1800 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Point25/Am241-SRC-01-25cm.spe` | Am241-SRC-01-25cm | Am-241 A=116600 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 300 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Point25/Ba133-SRC-01-25cm.spe` | Ba133-SRC-01-25cm | Ba-133 A=44100 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 600 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Point25/Cd109-SRC-01-25cm.spe` | Cd109-SRC-01-25cm | Cd-109 A=1.033E6 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 300 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Point25/Ce139-SRC-01-25cm.spe` | Ce139-SRC-01-25cm | Ce-139 A=116200 Bk dA=2% 01-10-2008 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 300 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Point25/Co57-SRC-01-25cm.spe` | Co57-SRC-01-25cm | Co-57 A=99500 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 600 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Point25/Co60-SRC-01-25cm.spe` | Co60-SRC-01-25cm | Co-60 A=107800 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 600 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Point25/Cs137-SRC-01-25cm.spe` | Cs137-SRC-01-25cm | Cs-137 A=94200 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 600 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Point25/Eu152-SRC-01-25cm.spe` | Eu152-SRC-01-25cm | Eu-152 A=46700 Bq dA=2% 01-10-2008 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 1800 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Point25/Mn54-SRC-01-25cm.spe` | Mn54-SRC-01-25cm | Mn-54 A=83600 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 600 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Point25/Sn113-SRC-01-25cm.spe` | Sn113-SRC-01-25cm | Sn-113 A=120000 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 600 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Point25/Th228-SRC-01-25cm.spe` | Th228-SRC-01-25cm | Th-228 A=37700 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 2700 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Point25/Y88-SRC-01-25cm.spe` | Y88-SRC-01-25cm | Y-88 A=134300 Bq dA=2% 01-10-2008 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 900 | 8192 | Bckg_1.spe |
| `Work/NM/HPGe(20_)/Spe/Th228--BadCalibre.spe` | Th228-SRC-01-25cm | Th-228 A=37700 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 2700 | 8192 | Bckg_1.spe |
| `Work/NM/U/U-Coaxial(Demo)/Spe/U_NBL_0043_93_CorrTlive.spe` | U_NBL_0043_93 | U-234 A=0.98 % (mass fraction) dA=0.29% 12-05-1987 | CRM146_9_8cm | GC6020 | 2008-09-29 | 106219 | 8192 | Background_LEGe_1_day_room_104.spe |
| `Work/NM/U/U-Coaxial(Demo)/Spe/ГСО-87-16000_Point_25cm.spe` | ГСО-87-16000 | — | Point | GEM30 | 2021-02-17 | 1564 | 16384 | Background_LEGe_1_day_room_104.spe |
| `Work/NM/U/U-Coaxial(Demo)/Spe/ГСО-87_25cm U-coaxial.spe` | ГСО-87 | — | Temp | HPGe(Coaxial) | 2021-02-17 | 1617 | 8192 | Background_LEGe_1_day_room_104.spe |
| `Work/NM/U/U-Coaxial(Demo)/Spe/ГСО-87_50-2614keV_25cm.spe` | ГСО-87 | — | 50-2614keV | HPGe(Coaxial) | 2021-02-17 | 1563 | 8192 | Background_LEGe_1_day_room_104.spe |

## Фоновые спектры

| Файл | Живое, с | Каналов | Дата |
|---|---:|---:|---|
| `Work/NM/HPGe(20_)/Spe/Background/Bckg_1.spe` | 3600 | 8192 | 2008-12-17 |
| `Work/NM/HPGe(20_)/Spe/Background/Bckg_15.spe` | 53994 | 8192 | 2008-12-17 |
| `Work/NM/HPGe(20_)/Spe/Background/Bckg_5.spe` | 17998 | 8192 | 2008-12-17 |
| `Work/NM/U/Planar/Data/Background_LEGe_1_day_room_104_D.spe` | 86400 | 16384 | 2013-01-24 |

`INDEX.json` рядом — те же данные плюс калибровка, реальное время и разобранный паспорт.
