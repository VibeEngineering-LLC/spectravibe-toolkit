# Handy_NaI — спектры LSRM и экспорт в BecqMoni

Спектров всего: **33**. Исходные `.spe` — в `lsrm/`, структура папок
сохранена. Экспорт в BecqMoni — в `becqmoni/`.

Сконвертировано: 33 (31 образцов, 2 фоновых).

## Очистка персональных данных

Из заголовков удалено поле `OPERATOR` — оно было заполнено в 2 файлах.
Остальные поля не изменялись: `COMMENT` содержит только паспортные активности
источников, персональных данных в нём нет.

## Перечень

| Файл | Источник | Паспорт | Геометрия | Детектор | Дата | Живое, с | Каналов | Фон |
|---|---|---|---|---|---|---:|---:|---|
| `Work/Handy/Handy(NaI)/Spe/Am-241  A=321200 Bq.spe` | Am-241  A=321200 Bq | Am-241 без контейнера  A=321200 Бк на 23.01.2007 г. (три источника) на расстоянии 15 см | Point-15cm | Gamma-1S-NB1 | 2007-01-23 | 297 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Ba-133   A=813600 Bq.spe` | Ba-133   A=813600 Bq | Ba-133   A=813600 Bq | Point-15cm | Gamma-1S-NB1 | 2007-01-23 | 292 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Ba-133 A=813000 Bq  in KT1-5.spe` | Ba-133 A=813000 Bq  in KT1-5 | Ва-133 A=813000 Бк  в KT1-5   (все 6 источников) дата аттестации источника 22.01.2006г. | Point-15cm | Gamma-1S-NB1 | 2007-01-31 | 594 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Ba-133 A=813600 Bq  in KT1-10.spe` | Ba-133 A=813600 Bq  in KT1-10 | Ba-133 A=813600 Бк  в KT1-10   (все 6 источников)  дата аттестации источника 22.01.2007 г. | Point-15cm | Gamma-1S-NB1 | 2007-01-31 | 595 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Co-60   A=197600 Bq.spe` | Co-60   A=197600 Bq | Co-60 без контейнера  A=197600 Бк на 22.01.2007 г. (три источника) и расстоянии 15 см | Point-15cm | Gamma-1S-NB1 | 2007-01-23 | 296 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Co-60 A=197600 Bq  in KT1-10.spe` | Co-60 A=197600 Bq  in KT1-10 | Co-60 A=197600 Бк  в KT1-10   (все 3 источника) дата аттестации источника 1.11.2006 г. | Point-15cm | Gamma-1S-NB1 | 2007-01-31 | 593 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Co-60 A=197600 Bq  in KT1-20.spe` | Co-60 A=197600 Bq  in KT1-20 | Co-60 A=197600 Бк  в KT1-20    дата аттестации источникa 22.01.2007 г. | Point-15cm | Gamma-1S-NB1 | 2007-02-02 | 594 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Co-60 A=197600 Bq  in KT1-5.spe` | Co-60 A=197600 Bq  in KT1-5 | Co-60 A=197600 Бк  в KT1-5   (все 3 источника) дата аттестации источника 22.01.2006г. | Point-15cm | Gamma-1S-NB1 | 2007-01-31 | 593 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Cs-134   A=95100 Bq.spe` | Cs-134   A=95100 Bq | Cs-134 без контейнера  A=95100 Бк на 1.11.2006 г.  на расстоянии 15 см | Point-15cm | Gamma-1S-NB1 | 2007-01-23 | 297 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Cs-134 A=95100 Bq  in KT1-10.spe` | Cs-134 A=95100 Bq  in KT1-10 | Cs-134 A=95100  Бк  в KT1-10     дата аттестации источника 1.11.2006 г. | Point-15cm | Gamma-1S-NB1 | 2007-01-31 | 595 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Cs-134 A=95100 Bq  in KT1-20.spe` | Cs-134 A=95100 Bq  in KT1-20 | Cs-134 A=95100 Бк  в KT1-20  дата аттестации источникa 1.11.2006 г. | Point-15cm | Gamma-1S-NB1 | 2007-02-02 | 595 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Cs-134 A=95100 Bq  in KT1-5.spe` | Cs-134 A=95100 Bq  in KT1-5 | Cs-134  A=95100 Бк  в KT1-5    дата аттестации источника 1.11.2006 г. | Point-15cm | Gamma-1S-NB1 | 2007-01-31 | 594 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Cs-137  A=74800 Bq.spe` | Cs-137  A=74800 Bq | Cs-137 без контейнера  A=74800 Бк на 1.11.2006 г. (четыре источника) на расстоянии 15 см | Point-15cm | Gamma-1S-NB1 | 2007-01-23 | 297 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/EU-152 A=243600 Bq  in KT1-20.spe` | EU-152 A=243600 Bq  in KT1-20 | EU-152 A=243600 Bq  in KT1-20 | Point-15cm | Gamma-1S-NB1 | 2007-02-02 | 595 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Eu-152  A=243600 Bq  in KT1-10.spe` | Eu-152  A=243600 Bq  in KT1-10 | Eu-152  A=243600  Бк  в KT1-10   (все 6 источников)  дата аттестации источников 22.01.2007 г. | Point-15cm | Gamma-1S-NB1 | 2007-01-31 | 594 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Eu-152  A=243600 Bq  in KT1-5.spe` | Eu-152  A=243600 Bq  in KT1-5 | Eu-152  A=243600 Бк  в KT1-5   (все 6 источников) дата аттестации источника 22.01.2007 г. | Point-15cm | Gamma-1S-NB1 | 2007-01-31 | 594 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Eu-152  A=243600 Bq.spe` | Eu-152  A=243600 Bq | Eu-152 без контейнера  A=243600 Бк на 23.01.2007 г. (шесть источника) на расстоянии 15 см | Point-15cm | Gamma-1S-NB1 | 2007-01-23 | 295 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Eu-152 Calibr.spe` | Eu-152 Calibr | — | Point-15cm | Gamma-1S-NB1 | 2007-02-14 | 576 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Na-22   A=91900 Bq.spe` | Na-22   A=919000 Bq | Na=22  A=91900 Бк, без контейнера, на расстоянии 15 см и дату аттестации 1.11.2006 | Point-15cm | Gamma-1S-NB1 | 2007-01-23 | 296 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Na-22  A=91900 Bq  in KT1-10.spe` | Na-22  A=91900 Bq  in KT1-10 | Na-22  A=91900 Бк  в KT1-10   дата аттестации источника 1.11.2006 г. | Point-15cm | Gamma-1S-NB1 | 2007-01-31 | 595 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Na-22  A=91900 Bq  in KT1-20.spe` | Na-22  A=91900 Bq  in KT1-20 | Na-22  A=91900 Бк  в KT1-20    дата аттестации источникa 1.11.2006г. | Point-15cm | Gamma-1S-NB1 | 2007-02-02 | 297 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Na-22 A=91900 Bq  in KT1-5.spe` | Na-22 A=91900 Bq  in KT1-5 | Na-22 A=91900 Бк в KT1-5   дата аттестации источника 1.11.2006 г. | Point-15cm | Gamma-1S-NB1 | 2007-01-31 | 594 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Th-228  A=67700 Bq  in KT1-5.spe` | Th-228  A=67700 Bq  in KT1-5 | Th-228  A=67700 Бк  в KT1-5   дата аттестации источника 1.11.2006 г. | Point-15cm | Gamma-1S-NB1 | 2007-01-31 | 595 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Th-228  A=67700 Bq in KT1-10.spe` | Th-228  A=67700 Bq in KT1-10 | Th-228  A=67700 Бк  в КТ1-10  на расстоянии 15 см дата аттестации источника 1.11.2006 г. | Point-15cm | Gamma-1S-NB1 | 2007-01-31 | 595 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Th-232  A=16700 Bq  in KT1-5.spe` | Th-232  A=16700 Bq  in KT1-5 | Th-232  A=16700 Бк  в KT1-5   дата аттестации источника 1.11.2006 г. | Point-15cm | Gamma-1S-NB1 | 2007-01-31 | 595 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Th-232  A=16700 Bq.spe` | Th-232  A=16700 Bq | Th-232 без контейнера  A=16700 Бк на 23.01.2007 г. (четыре источника) на расстоянии 15 см | Point-15cm | Gamma-1S-NB1 | 2007-01-23 | 297 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Th-232Calibr.spe` | Th-232  A=67700 Bq in KT1-10 | Th-232  A=67700 Бк  в КТ1-10  на расстоянии 15 см дата аттестации источника 1.11.2006 г. | Точечный | Handy_NaI | — | 592 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Ti-44   A=72400 Bq.spe` | Ti-44   A=72400 Bq | Ti-44 без контейнера  A=72400 Бк на 1.11.2006 г. и расстоянии 15 см | Point-15cm | Gamma-1S-NB1 | 2007-01-23 | 296 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Ti-44 A=72400 Bq  in KT1-10.spe` | Ti-44 A=72400 Bq  in KT1-10 | Ti-44 A=72400 Бк  в KT1-10   дата аттестации источника 1.11.2006 г. | Point-15cm | Gamma-1S-NB1 | 2007-01-31 | 595 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Ti-44 A=72400 Bq  in KT1-20.spe` | Ti-44 A=72400 Bq  in KT1-20 | Ti-44 A=72400 Бк  в KT1-20    дата аттестации источникa 1.11.2006г. | Point-15cm | Gamma-1S-NB1 | 2007-02-02 | 595 | 1024 | Background140207.spe |
| `Work/Handy/Handy(NaI)/Spe/Ti-44 A=72400 Bq  in KT1-5.spe` | Ti-44 A=72400 Bq  in KT1-5 | Ti-44  A=72400 Бк  в KT1-5   дата аттестации источника 1.11.2006 г. | Point-15cm | Gamma-1S-NB1 | 2007-01-31 | 594 | 1024 | Background140207.spe |

## Фоновые спектры

| Файл | Живое, с | Каналов | Дата |
|---|---:|---:|---|
| `Work/Handy/Handy(NaI)/Data/Background140207.spe` | 297 | 1024 | 2007-02-14 |
| `Work/Handy/Handy(NaI)/Spe/Background140207.spe` | 297 | 1024 | 2007-02-14 |

`INDEX.json` рядом — те же данные плюс калибровка, реальное время и разобранный паспорт.
