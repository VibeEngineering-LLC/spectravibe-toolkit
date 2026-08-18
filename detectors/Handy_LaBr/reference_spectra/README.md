# Handy_LaBr — спектры LSRM и экспорт в BecqMoni

Спектров всего: **41**. Исходные `.spe` — в `lsrm/`, структура папок
сохранена. Экспорт в BecqMoni — в `becqmoni/`.

Сконвертировано: 41 (24 образцов, 17 фоновых).

## Очистка персональных данных

Из заголовков удалено поле `OPERATOR` — оно было заполнено в 41 файлах.
Остальные поля не изменялись: `COMMENT` содержит только паспортные активности
источников, персональных данных в нём нет.

## Перечень

| Файл | Источник | Паспорт | Геометрия | Детектор | Дата | Живое, с | Каналов | Фон |
|---|---|---|---|---|---|---:|---:|---|
| `Work/Handy/Handy(LaBr)/Spe/Am241_#14-05_24sm.spe` | Am241_#14-05_24sm | A=100000 Bq dA=3% 06-06-2005 | Point24 | BrilLanCe380 | 2008-07-03 | 3600 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Ba133_#SRC-13_24sm.spe` | Ba133_#SRC-13_24sm | A=92700 Bq dA=3% 13-04-2002 | Point24 | BrilLanCe380 | 2008-07-04 | 3600 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Cd109_#SRC-14_24sm.spe` | Cd109_#SRC-14_24sm | Cd-109 A=348500 Bq dA=5% 27-12-2006 | Point24 | BrilLanCe380 | 2008-07-02 | 3600 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Ce139_#SRC-15_24sm_1.spe` | Ce139_#SRC-15_24sm_1 | A=149500 Bq dA=2% 31-03-2006 | Point24 | BrilLanCe380 | 2008-07-03 | 3600 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Ce139_#SRC-15_24sm_10.spe` | Ce139_#SRC-15_24sm_10 | A=149500 Bq dA=2% 31-03-2006 | Point24 | BrilLanCe380 | 2008-07-03 | 36000 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Ce139_#SRC-15_24sm_2.spe` | Ce139_#SRC-15_24sm_2 | подложка для источника на расстоянии 24 cм от крышки детектора  Ce-139 A=149500 Бк dA=2% 31-03-06 | Точечная | BrilLanCe380 | 2008-07-03 | 7200 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Ce139_#SRC-15_24sm_3.spe` | Ce139_#SRC-15_24sm_3 | подложка для источника на расстоянии 24 cм от крышки детектора  Ce-139 A=149500 Бк dA=2% 31-03-06 | Точечная | BrilLanCe380 | 2008-07-03 | 10800 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Ce139_#SRC-15_24sm_4.spe` | Ce139_#SRC-15_24sm_4 | подложка для источника на расстоянии 24 cм от крышки детектора  Ce-139 A=149500 Бк dA=2% 31-03-06 | Точечная | BrilLanCe380 | 2008-07-03 | 14400 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Ce139_#SRC-15_24sm_5.spe` | Ce139_#SRC-15_24sm_5 | Ce-139 A=149500 Bq dA=2% 31-03-2006 | Point24 | BrilLanCe380 | 2008-07-03 | 18000 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Ce139_#SRC-15_24sm_6.spe` | Ce139_#SRC-15_24sm_6 | подложка для источника на расстоянии 24 cм от крышки детектора  Ce-139 A=149500 Бк dA=2% 31-03-06 | Точечная | BrilLanCe380 | 2008-07-03 | 21600 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Ce139_#SRC-15_24sm_7.spe` | Ce139_#SRC-15_24sm_7 | подложка для источника на расстоянии 24 cм от крышки детектора  Ce-139 A=149500 Бк dA=2% 31-03-06 | Точечная | BrilLanCe380 | 2008-07-03 | 25200 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Ce139_#SRC-15_24sm_8.spe` | Ce139_#SRC-15_24sm_8 | подложка для источника на расстоянии 24 cм от крышки детектора  Ce-139 A=149500 Бк dA=2% 31-03-06 | Точечная | BrilLanCe380 | 2008-07-03 | 28800 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Ce139_#SRC-15_24sm_9.spe` | Ce139_#SRC-15_24sm_9 | подложка для источника на расстоянии 24 cм от крышки детектора  Ce-139 A=149500 Бк dA=2% 31-03-06 | Точечная | BrilLanCe380 | 2008-07-03 | 32400 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Ce139_#SRC-16_24sm.spe` | Ce139_#SRC-16_24sm | подложка для источника на расстоянии 1 мм от крышки детектора #SRC-16 Ce-139 A=350600 Бк dA=5% 01-11-2003 | Точечная | BrilLanCe380 | 2008-07-03 | 2364 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Co57_#SRC-15_24sm.spe` | Co57_#SRC-15_24sm | A=230500 Bq dA=2% 31-03-2006 | Point24 | BrilLanCe380 | 2008-07-02 | 1800 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Co60_#SRC-13_24sm.spe` | Co60_#SRC-13_24sm | A=105000 Bq dA=3% 13-04-2002 | Point24 | BrilLanCe380 | 2008-07-03 | 1800 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Cs137_#SRC-13_24sm.spe` | Cs137_#SRC-13_24sm | Cs-137 A=98950 Bq dA=3% 13-04-2002 | Point24 | BrilLanCe380 | 2008-07-03 | 1800 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Eu152_#SRC-13_24sm.spe` | Eu152_#SRC-13_24sm | A=84400 Bq dA=3% 13-04-2002 | Point24 | BrilLanCe380 | 2008-07-04 | 3600 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Mn54_#SRC-15_24sm.spe` | Mn54_#SRC-15_24sm | Mn-54 A=210000 Bq dA=2% 31-03-2006 | Point24 | BrilLanCe380 | 2008-07-02 | 1800 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Na22_#SRC-15_24sm.spe` | Na22_#SRC-15_24sm | Na-22 A=119700 Bq dA=2% 31-03-2006 | Point24 | BrilLanCe380 | 2008-07-02 | 1800 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Na22_#SRC-13_24sm.spe` | Na22_#SRC-13_24sm | Na-22 A=42500 Bq dA=3% 13-04-2002 | Point24 | BrilLanCe380 | 2008-07-04 | 3600 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Th228_#SRC-17_24sm.spe` | Th228_#SRC-17_24sm | Th-228 A=28900 Bq dA=3% 17-06-2005 | Point24 | BrilLanCe380 | 2008-07-03 | 3600 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Ti44_#SRC-18_24sm.spe` | Ti44_#SRC-18_24sm | Ti-44 A=10000 Bq dA=2% 02-04-2004 | Point24 | BrilLanCe380 | 2008-07-03 | 3600 | 1024 | Background_1.spe |
| `Work/Handy/Handy(LaBr)/Spe/Zn65_#SRC-15_24sm.spe` | Zn65_#SRC-15_24sm | Zn-65 A=156300 Bq dA=2% 31-03-2006 | Point24 | BrilLanCe380 | 2008-07-03 | 3600 | 1024 | Background_1.spe |

## Фоновые спектры

| Файл | Живое, с | Каналов | Дата |
|---|---:|---:|---|
| `Work/Handy/Handy(LaBr)/Spe/Background/Background_1.spe` | 3600 | 1024 | 2008-07-02 |
| `Work/Handy/Handy(LaBr)/Spe/Background/Background_10.spe` | 36000 | 1024 | 2008-07-02 |
| `Work/Handy/Handy(LaBr)/Spe/Background/Background_2.spe` | 7200 | 1024 | 2008-07-02 |
| `Work/Handy/Handy(LaBr)/Spe/Background/Background_3.spe` | 10800 | 1024 | 2008-07-02 |
| `Work/Handy/Handy(LaBr)/Spe/Background/Background_4.spe` | 14400 | 1024 | 2008-07-02 |
| `Work/Handy/Handy(LaBr)/Spe/Background/Background_5.spe` | 18000 | 1024 | 2008-07-02 |
| `Work/Handy/Handy(LaBr)/Spe/Background/Bckg_10.spe` | 36000 | 1024 | 2008-07-02 |
| `Work/Handy/Handy(LaBr)/Spe/Background/MarBckg_14.spe` | 51650 | 1024 | 2008-06-30 |
| `Work/Handy/Handy(LaBr)/Spe/Background/Фон без защиты_1.spe` | 3600 | 1024 | 2008-07-02 |
| `Work/Handy/Handy(LaBr)/Spe/Background/Фон без защиты_2.spe` | 7200 | 1024 | 2008-07-02 |
| `Work/Handy/Handy(LaBr)/Spe/Background/Фон без защиты_3.spe` | 10800 | 1024 | 2008-07-02 |
| `Work/Handy/Handy(LaBr)/Spe/Background/Фон без защиты_4.spe` | 14400 | 1024 | 2008-07-02 |
| `Work/Handy/Handy(LaBr)/Spe/Background/Фон без защиты_5.spe` | 18000 | 1024 | 2008-07-02 |
| `Work/Handy/Handy(LaBr)/Spe/Background/Фон без защиты_6.spe` | 21600 | 1024 | 2008-07-02 |
| `Work/Handy/Handy(LaBr)/Spe/Background/Фон без защиты_7.spe` | 25200 | 1024 | 2008-07-02 |
| `Work/Handy/Handy(LaBr)/Spe/Background/Фон без защиты_8.spe` | 28800 | 1024 | 2008-07-02 |
| `Work/Handy/Handy(LaBr)/Spe/Background/Фон без защиты_9.spe` | 32400 | 1024 | 2008-07-02 |

`INDEX.json` рядом — те же данные плюс калибровка, реальное время и разобранный паспорт.
