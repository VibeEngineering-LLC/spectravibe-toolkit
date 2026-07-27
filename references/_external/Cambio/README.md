# Cambio — Spectroscopic file conversion tool (Sandia National Labs)

**Источник:** https://github.com/sandialabs/cambio
**Provenance:** `gh repo view sandialabs/cambio` + `gh api repos/sandialabs/cambio/contents/README.md`
**Дата fetch:** 2026-06-24

## Метаданные репо

| Поле | Значение |
|---|---|
| Описание | Spectroscopic file conversion tool |
| Default branch | `master` |
| Последний push | 2025-11-05 (активен) |
| Создан | 2019-03-19 |
| Архивирован | нет |
| Лицензия | LGPL v2.1 |
| Звёзды | 20 |
| Языки | C++ (783 740), C (1 154 010), CMake (34 786), Obj-C++ (1 272), Java (2 541), QMake (3 319), Batch (6 276) |
| Topics | snl-applications, snl-visualization, scr-988 |
| Автор | William Johnson |
| Copyright | 2018 NTESS / Sandia, Contract DE-NA0003525, SCR #SRC-19 |

## Что делает (verbatim из README)

> Cambio converts spectrum files from nearly all common handheld and lab-based
> spectroscopic gamma radiation detectors, radiation portal monitors, or search
> systems to a format of your choice (N42, PCF, CSV, TXT, CHN, SPC, HTML, and more).
>
> Cambio can be built either as a GUI that allows previewing the spectrum files,
> editing the meta-data, and interactively choosing the output format options,
> or Cambio can be built as a command line only utility. The command line utility
> is useful for calling from batch scripts, or calling from other programs to
> take care of reading the hundreds of potential formats, making it so your
> program only needs to read in a single format.

## Статус (verbatim, важно)

> The GUI version of cambio is no longer actively maintained — please see
> [InterSpec](https://github.com/sandialabs/InterSpec) as a far more capable
> replacement.
>
> The command line version will continue to be maintained.

Т.е. **GUI deprecated**, CLI поддерживается.

## Родственные проекты Sandia

| Репо | Что | URL |
|---|---|---|
| **SpecUtils** | C++ библиотека парсинга спектров (база Cambio); биндинги для нескольких языков; Python — `pip install SandiaSpecUtils` | https://github.com/sandialabs/SpecUtils |
| **InterSpec** | Современный полноценный аналитический пакет; рекомендованная замена GUI Cambio | https://github.com/sandialabs/InterSpec |

## Релевантность проекту `gamma-spectrum-analysis`

1. **Конвертация форматов.** Cambio (CLI) и `SandiaSpecUtils` (Python) умеют читать
   "сотни" форматов гамма-спектров (handheld, lab-based, RPM, search systems) и
   писать в N42/PCF/CSV/TXT/CHN/SPC/HTML. Полезно как референс для расширения
   собственных парсеров спектров (`scripts/gamma/io/` или эквивалент) и как
   ground-truth конвертер при тестах форматов.
2. **N42 reference impl.** Cambio/SpecUtils — авторитетная C++ реализация ANSI/IEEE
   N42.42 (стандарт). Можно сверять собственный N42-парсер/райтер с поведением
   SpecUtils как с эталоном.
3. **PCF (Sandia format).** PCF mention здесь — если когда-то понадобится
   читать/писать PCF, у Sandia есть и спецификация
   (https://prod-ng.sandia.gov/techlib-noauth/access-control.cgi/2017/179107.pdf),
   и эталонная реализация в SpecUtils.
4. **InterSpec — рекомендованный GUI.** Если оператору нужен внешний
   просмотрщик/анализатор спектров с GUI для cross-check работы нашего пайплайна
   (peak search, identification, MDA) — InterSpec, а не Cambio-GUI.

## Build / install (verbatim)

Pre-compiled exe — в `releases/` репо. Для исходной сборки:

- C++14 compiler
- Boost 1.65 – 1.84
- Qt 5.15 (только для GUI)
- CMake
- SpecUtils (клонируется рядом)

CLI-only: `-DBUILD_CAMBIO_COMMAND_LINE=ON -DBUILD_CAMBIO_GUI=OFF`.

Установка — копия бинаря, инсталлятор не нужен.

## Privacy (verbatim, важно для security/compliance)

> Cambio does not collect any user information or statistics, and it does not
> send or receive any information over the network (i.e., nothing is downloaded
> to, or leaves from your computer).

Без сетевой телеметрии — можно использовать на машине оператора без опсек-рисков
(в отличие от инструментов, требующих cloud-вызовов).

## Когда вспомнить про Cambio в нашей работе

- если возникает спектр в формате, который наш собственный парсер не поддерживает,
  и нужна быстрая конвертация в наш входной формат (.spe / .n42 / CSV);
- если нужен независимый source-of-truth для парсинга N42/PCF/CHN/SPC при отладке
  собственного I/O;
- если оператор хочет визуальный просмотр спектра внешним инструментом —
  направлять на **InterSpec**, не на Cambio-GUI (deprecated).

## Не использовать для

- анализа (peak search / identification / efficiency / MDA) — это не аналитический
  пакет, а конвертер. Для анализа — InterSpec.
- замены нашего пайплайна.

---
SCR #SRC-19. LGPL-2.1.