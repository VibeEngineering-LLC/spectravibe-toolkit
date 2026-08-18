# Reference Kits — Gamma-1C NaI 63×63 USB №SN-01

> F-331 / v1.18.18.5 — Канонические комплекты «образец + фон»
> для использования в регрессионных и приёмочных тестах.

## Структура

```
reference_kits/
├── {geometry}/
│   └── {nuclide}/
│       ├── sample_{...}.spe        — спектр эталонного источника
│       └── background_{...}.spe    — релевантный усреднённый фон
```

Все остальные исторические спектры из легаси-папки
`Gamma-1C_NaI_63x63_USB_SN-01/` перенесены в
`detectors/Gamma-1C/reference_spectra/archive/` с сохранением
относительной структуры.

## Манифест комплектов

| Geometry | Nuclide | Sample (md5) | Background (md5) |
|---|---|---|---|
| Marinelli_1L | Cs-137 | `sample_M_cs_легкий_2001-2005.spe` (3a5accc6e7) | `background_bg_2016_marinelli_water_marinelli.spe` (f792057a8a) |
| Marinelli_1L | K-40 | `sample_M_k_легкий_2001-2005.spe` (efa302f064) | `background_bg_2016_marinelli_water_marinelli.spe` (f792057a8a) |
| Marinelli_1L | Ra-226 | `sample_M_ra_легкий_2001-2007.spe` (c8c15cee5d) | `background_bg_2016_marinelli_water_marinelli.spe` (f792057a8a) |
| Marinelli_1L | Th-232 | `sample_M_th_легкий_2001-2005.spe` (f447ff9916) | `background_bg_2016_marinelli_water_marinelli.spe` (f792057a8a) |
| Point_5cm | Am-241 | `sample_Am-241 42.13_Точечная-5см_5cm.spe` (156ba48f09) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Point_5cm | Ba-133 | `sample_Ba-133 #SRC-07_Точечная-5см_5cm.spe` (ad5bad5caf) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Point_5cm | Bi-207 | `sample_Bi-207__176_04_2017_Точечная-5см_5cm.spe` (04763ab36a) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Point_5cm | Cd-109 | `sample_Cd-109 #SRC-07_Точечная-5см_5cm.spe` (4f4107782b) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Point_5cm | Ce-139 | `sample_Ce-139_591_Точечная-5см_5cm.spe` (20ae97dde7) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Point_5cm | Co-57 | `sample_Co-57 #SRC-07_Точечная-5см_5cm.spe` (b5ed893095) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Point_5cm | Co-60 | `sample_Co-60 #SRC-07_Точечная-5см_5cm.spe` (0b1cc83e5b) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Point_5cm | Cs-137 | `sample_Cs-137 #SRC-07_Точечная-5см_5cm.spe` (9feedbe1af) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Point_5cm | Eu-152 | `sample_Eu-152 #SRC-07_Точечная-5см_5cm.spe` (dbe63ba336) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Point_5cm | Mn-54 | `sample_Mn-54_587_Точечная-5см_5cm.spe` (59d411fc13) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Point_5cm | Na-22 | `sample_Na-22_585_Точечная-5см_5cm.spe` (7e443018a4) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Point_5cm | Th-228 | `sample_Th-228 #SRC-07_Точечная-5см_5cm.spe` (50c79cc070) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Point_5cm | Y-88 | `sample_Y-88_589_Точечная-5см_5cm.spe` (047b75cfb8) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Point_5cm | Zn-65 | `sample_Zn-65__342_2019_Точечная-5см_5cm.spe` (6b98afa2ae) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Point_25cm | Am-241 | `sample_Am-241 42.13_Точечная-25см_25cm.spe` (6f033ecdbc) | `background_bg_2016_open_lid_point25cm.spe` (5d6a9751da) |
| Point_25cm | Ba-133 | `sample_Ba-133 #SRC-07_Точечная-25см_25cm.spe` (e1218aa5f9) | `background_bg_2016_open_lid_point25cm.spe` (5d6a9751da) |
| Point_25cm | Cd-109 | `sample_Cd-109 #SRC-07_Точечная-25см_25cm.spe` (bbc7ece1e7) | `background_bg_2016_open_lid_point25cm.spe` (5d6a9751da) |
| Point_25cm | Ce-139 | `sample_Ce-139_591_Точечная-25см_25cm.spe` (b6dc5d9dae) | `background_bg_2016_open_lid_point25cm.spe` (5d6a9751da) |
| Point_25cm | Co-60 | `sample_Co-60 #SRC-07_Точечная-25см_25cm.spe` (a4818d5bf4) | `background_bg_2016_open_lid_point25cm.spe` (5d6a9751da) |
| Point_25cm | Cs-137 | `sample_Cs-137 №SRC-01_Точечная-25см_25cm.spe` (5fb8d60b4e) | `background_bg_2016_open_lid_point25cm.spe` (5d6a9751da) |
| Point_25cm | Eu-152 | `sample_Eu-152 #SRC-07_Точечная-25см_25cm.spe` (0e46173dd8) | `background_bg_2016_open_lid_point25cm.spe` (5d6a9751da) |
| Point_25cm | Mn-54 | `sample_Mn-54_587_Точечная-25см_25cm.spe` (678e321c70) | `background_bg_2016_open_lid_point25cm.spe` (5d6a9751da) |
| Point_25cm | Na-22 | `sample_Na-22 #01.22_Точечная-25см_25cm.spe` (0fa53e2ecb) | `background_bg_2016_open_lid_point25cm.spe` (5d6a9751da) |
| Point_25cm | Th-228 | `sample_Th-228 №309_Точечная-25см_25cm.spe` (7a00a12d61) | `background_bg_2016_open_lid_point25cm.spe` (5d6a9751da) |
| Point_25cm | Y-88 | `sample_Y-88 №SRC-02_Точечная-25см_25cm.spe` (7701a3f5ad) | `background_bg_2016_open_lid_point25cm.spe` (5d6a9751da) |
| Petri_60mL | Cs-137 | `sample_Cs137_420-7-14_Петри-60мл_0cm.spe` (93a2fcac71) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Petri_60mL | K-40 | `sample_K40_420-7-20_Петри-60мл_0cm.spe` (18fb2a7cf1) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Petri_60mL | Ra-226 | `sample_Ra226_420-7-18_Петри-60мл_0cm.spe` (bbc841a573) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Petri_60mL | Th-232 | `sample_Th232_420-7-17_Петри-60мл_0cm.spe` (13063ba24d) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Denta_120mL | Cs-137 | `sample_Cs137_420-7-14_Дента-120мл_0cm.spe` (3b95078775) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Denta_120mL | K-40 | `sample_K40_420-7-20_Дента-120мл_0cm.spe` (f832d0dbae) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Denta_120mL | Ra-226 | `sample_Ra226_420-7-18_Дента-120мл_0cm.spe` (271d501ba8) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |
| Denta_120mL | Th-232 | `sample_Th232_420-7-17_Дента-120мл_0cm.spe` (fd06b486f1) | `background_bg_2016_empty_shield_point5cm.spe` (fd7bde8eb0) |

## Использование в тестах

```python
from pathlib import Path
KIT = Path('detectors/Gamma-1C/reference_spectra/reference_kits')
from gamma.reporting import analyze_and_report

# Marinelli 1L Cs-137 kit
kit = KIT / 'Marinelli_1L' / 'Cs-137'
sample = next(kit.glob('sample_*.spe'))
bg = next(kit.glob('background_*.spe'))

artefacts = analyze_and_report(
    str(sample),
    output_dir='./out',
    background_path=str(bg),
    sample_mass_kg=0.570,
)
```

Passport activities для Marinelli автоматически читаются из
COMMENT-поля .spe (F-330 v1.18.18.4) — `passport_activity_Bq`
передавать не нужно.