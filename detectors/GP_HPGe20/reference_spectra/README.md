# GP_HPGe20 — спектры LSRM и экспорт в BecqMoni

Спектров всего: **32**. Исходные `.spe` — в `lsrm/`, структура папок
сохранена. Экспорт в BecqMoni — в `becqmoni/`.

Сконвертировано: 32 (29 образцов, 3 фоновых).

## Очистка персональных данных

Из заголовков удалено поле `OPERATOR` — оно было заполнено в 25 файлах.
Остальные поля не изменялись: `COMMENT` содержит только паспортные активности
источников, персональных данных в нём нет.

## Калибровка, которую BecqMoni прочитает неверно

`PolynomialEnergyCalibration.ChannelToEnergy()` разбирает только степени 2, 3 и 4;
всё, что выше, молча уходит в линейную ветку `c[1]·n + c[0]`. У перечисленных
файлов LSRM записал полином **5-й степени** — коэффициенты сохранены как в
источнике, но энергетическая шкала в BecqMoni будет неверной. Для расчётов
берите калибровку из `INDEX.json` (поле `energy_cal`).

- `Work/GP/HPGe(20_)/Spe/Background/Bckg_1.spe` — степень 5
- `Work/GP/HPGe(20_)/Spe/Background/Bckg_15.spe` — степень 5
- `Work/GP/HPGe(20_)/Spe/Marinelli/m08085_mix09.spe` — степень 5
- `Work/GP/HPGe(20_)/Spe/Point25/Y88-SRC-01-25cm.spe` — степень 5

## Перечень

| Файл | Источник | Паспорт | Геометрия | Детектор | Дата | Живое, с | Каналов | Фон |
|---|---|---|---|---|---|---:|---:|---|
| `Work/GP/HPGe(20_)/Spe/Marinelli/m08085_mix09.spe` | m08085_mix09 | Ti-44 A=3530 Bq /kg dA=6% 26-11-08 | Marinelli | GEM20P4-70 #SN-02 | 2008-12-19 | 3516 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Marinelli/m_cs06.spe` | m_cs06 | Cs-137 A=1760 Bk/kg dA=5% 24-05-02 | MARINELLI | GEM20P4-70 #SN-02 | 2008-12-19 | 1800 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Marinelli/m_cs16.spe` | m_cs16 | Cs-137 A=1650 Bk/kg dA=5% 17-09-07 | MARINELLI | GEM20P4-70 #SN-02 | 2008-12-17 | 1800 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Marinelli/m_k06.spe` | m_k06 | K-40 A=2530 Bk/kg dA=6% 24-05-02 | MARINELLI | GEM20P4-70 #SN-02 | 2008-12-19 | 3600 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Marinelli/m_k16.spe` | m_k16 | K-40 A=1950 Bk/kg dA=6% 17-09-07 | MARINELLI | GEM20P4-70 #SN-02 | 2008-12-17 | 3600 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Marinelli/m_ra06.spe` | m_ra06 | Ra-226 A=1780 Бк/кг dA=6% 24-05-02 | MARINELLI | GEM20P4-70 #SN-02 | 2008-12-19 | 3600 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Marinelli/m_ra16.spe` | m_ra16 | Ra-226 A=1790 Bk/kg dA=6% 17-09-07 | MARINELLI | GEM20P4-70 #SN-02 | 2008-12-17 | 1800 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Marinelli/m_th06.spe` | m_th06 | Th-232 A=2240 Bk/kg dA=6% 24-05-02 | MARINELLI | GEM20P4-70 #SN-02 | 2008-12-19 | 3585 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Marinelli/m_th16.spe` | m_th16 | Th-232 A=1940 Bk/kg dA=6% 17-09-07 | MARINELLI | GEM20P4-70 #SN-02 | 2008-12-17 | 1800 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Point25/Am241-SRC-01-25cm.spe` | Am241-SRC-01-25cm | Am-241 A=116600 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 300 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Point25/Ba133-SRC-01-25cm.spe` | Ba133-SRC-01-25cm | Ba-133 A=44100 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 600 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Point25/Cd109-SRC-01-25cm.spe` | Cd109-SRC-01-25cm | Cd-109 A=1.033E6 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 300 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Point25/Ce139-SRC-01-25cm.spe` | Ce139-SRC-01-25cm | Ce-139 A=116200 Bk dA=2% 01-10-2008 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 300 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Point25/Co57-SRC-01-25cm.spe` | Co57-SRC-01-25cm | Co-57 A=99500 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 600 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Point25/Co60-SRC-01-25cm.spe` | Co60-SRC-01-25cm | Co-60 A=107800 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 600 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Point25/Cs137-SRC-01-25cm.spe` | Cs137-SRC-01-25cm | Cs-137 A=94200 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 600 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Point25/Eu152-SRC-01-25cm.spe` | Eu152-SRC-01-25cm | Eu-152 A=46700 Bq dA=2% 01-10-2008 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 1800 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Point25/Mn54-SRC-01-25cm.spe` | Mn54-SRC-01-25cm | Mn-54 A=83600 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 600 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Point25/Sn113-SRC-01-25cm.spe` | Sn113-SRC-01-25cm | Sn-113 A=120000 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 600 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Point25/Th228-SRC-01-25cm.spe` | Th228-SRC-01-25cm | Th-228 A=37700 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 2700 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Point25/Y88-SRC-01-25cm.spe` | Y88-SRC-01-25cm | Y-88 A=134300 Bq dA=2% 01-10-2008 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 900 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/Th228--BadCalibre.spe` | Th228-SRC-01-25cm | Th-228 A=37700 Bk dA=2% 01-10-08 | Point25 | GEM20P4-70 #SN-02 | 2008-12-17 | 2700 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/joint processing/Barrel_10cm_D1.spe` | Barrel-SKG(Water) Co60 e07 D2 10cm | — | Barrel-SKG(Water) | GEM20P4-70 #SN-02 | 2021-07-23 | 12 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/joint processing/Barrel_10cm_D2.spe` | Barrel-SKG(Water) Co60 e07 D2 10cm | — | Barrel-SKG(Water) | GEM20P4-70 #SN-02 | 2021-07-23 | 12 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/joint processing/Barrel_10cm_D3.spe` | Barrel-SKG(Water) Co60 e07 D2 10cm | — | Barrel-SKG(Water) | GEM20P4-70 #SN-02 | 2021-07-23 | 12 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/joint processing/NZK_50cm_D1.spe` | GEM-30_D1_Front50_50 | — | NZK-150 | GEM20P4-70 #SN-02 | 2021-07-23 | 129 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/joint processing/NZK_50cm_D2.spe` | GEM-30_D1_Front50_50 | — | NZK-150 | GEM20P4-70 #SN-02 | 2021-07-23 | 129 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/joint processing/NZK_50cm_D3.spe` | GEM-30_D1_Front50_50 | — | NZK-150 | GEM20P4-70 #SN-02 | 2021-07-23 | 129 | 8192 | Bckg_1.spe |
| `Work/GP/HPGe(20_)/Spe/joint processing/NZK_50cm_D4.spe` | GEM-30_D1_Front50_50 | — | NZK-150 | GEM20P4-70 #SN-02 | 2021-07-23 | 129 | 8192 | Bckg_1.spe |

## Фоновые спектры

| Файл | Живое, с | Каналов | Дата |
|---|---:|---:|---|
| `Work/GP/HPGe(20_)/Spe/Background/Bckg_1.spe` | 3600 | 8192 | 2008-12-17 |
| `Work/GP/HPGe(20_)/Spe/Background/Bckg_15.spe` | 53994 | 8192 | 2008-12-17 |
| `Work/GP/HPGe(20_)/Spe/Background/Bckg_5.spe` | 17998 | 8192 | 2008-12-17 |

`INDEX.json` рядом — те же данные плюс калибровка, реальное время и разобранный паспорт.
