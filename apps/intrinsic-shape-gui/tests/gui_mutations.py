#!/usr/bin/env python3
"""Мутационная поверка тестов окна (tests/test_gui_detectors.py): каждая точечная мутация логики GUI портит одну
проверяемую величину в копии src/, тесты прогоняются на копии (переменная ISGUI_SRC); мутация убита, если покраснел хотя бы
один тест. Исправная версия обязана быть зелёной. Итог -- таблица в файл отчёта (аргумент 1, по умолчанию
<корень репозитория>/audit/gui_4detectors_mutations.md). Запуск: python apps/intrinsic-shape-gui/tests/gui_mutations.py.
Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-11-16."""

import os
import re
import sys
import shutil
import pathlib
import subprocess
import tempfile
import time

sys.stdout.reconfigure(encoding="utf-8")

ROOT = pathlib.Path(__file__).resolve().parents[1]  # apps/intrinsic-shape-gui
TK = ROOT.parents[1]  # корень репозитория
TEST = ROOT / "tests" / "test_gui_detectors.py"

# Список мутаций: (имя, цель на русском, файл в src/isgui/, старый фрагмент, новый фрагмент)
MUTS = [
    ("M1_profiles", "список профилей: только La-138", "backend_nai.py",
     "profiles.list_profiles(all_methods=True)", "profiles.list_profiles()"),
    ("M2_tabs", "скрытие вкладок чужого метода", "main_window.py",
     "            self.tabs.setTabVisible(i, is_la)\n            self.tabs.setTabEnabled(i, is_la)",
     "            self.tabs.setTabVisible(i, True)\n            self.tabs.setTabEnabled(i, True)"),
    ("M3_banner", "пометка «экспериментально» на вкладке", "nai_panel.py",
     'EXPERIMENTAL_BANNER = "ЭКСПЕРИМЕНТАЛЬНО: ', 'EXPERIMENTAL_BANNER = "Метод NaI: '),
    ("M3b_summary", "пометка «экспериментально» в итоге", "backend_nai.py",
     "<b>ЭКСПЕРИМЕНТАЛЬНО: шкала не проверяется</b>", "<b>шкала не проверяется</b>"),
    ("M4_guard", "выбор метода: проверка метода команды", "main_window.py",
     "        if method != cur:", "        if False:"),
    ("M4b_naimethod", "выбор метода: запуск NaI как La-138", "main_window.py",
     "self.nai.show_result, method=backend_nai.NAI_METHOD)", "self.nai.show_result)"),
    ("M5_detect", "автоопределение: прибор не распознаётся", "backend_nai.py",
     "        pid, msg = profiles.detect(dev)", '        pid, msg = None, ""'),
    ("M5b_autoswitch", "автоопределение: профиль не переключается", "main_window.py",
     'if self.auto_det.isChecked() and det["pid"] and det["pid"] != cur:', "if False:"),
    ("M6_settings", "сохранение выбора профиля", "main_window.py",
     'self.settings.setValue("detector", pid)', 'self.settings.setValue("detector_off", pid)'),
    ("M6b_autoflag", "сохранение флага автовыбора", "main_window.py",
     '"true" if on else "false"', '"true"'),
    ("M7_cancel", "отмена", "backend_nai.py",
     '    if robust.CANCEL["flag"]:', "    if False:"),
    ("M8_mismatch", "предупреждение о несовпадении (распознан другой)", "backend_nai.py",
     '    if det["pid"] is not None and det["pid"] != selected_pid:', "    if False:"),
    ("M8b_mismatch_win", "предупреждение о несовпадении в окне", "main_window.py",
     'self.det_warning = backend_nai.mismatch_text(cur, det)', 'self.det_warning = ""'),
    ("M9_peaks_nuclides", "режим по линиям: явный список нуклидов", "backend_nai.py",
     "            nuclides = nuclides or sorted(n for n in nuclide_names(pid) if n in anchors)",
     "            pass"),
]


def make_copy(tmp):
    """Создаёт копию каталога src/ (без __pycache__) в tmp/src и возвращает путь."""
    dst = tmp / "src"
    shutil.copytree(ROOT / "src", dst, ignore=shutil.ignore_patterns("__pycache__"))
    return dst


def apply(src, rel, old, new):
    """Применяет точечную мутацию: заменяет old на new в файле src/isgui/rel."""
    p = src / "isgui" / rel
    with open(p, encoding="utf-8", newline="") as f:
        text = f.read()
    n = text.count(old)
    if n == 0:
        # Повторная попытка с CRLF-окончаниями (main_window.py использует CRLF)
        old2 = old.replace("\n", "\r\n")
        new2 = new.replace("\n", "\r\n")
        n = text.count(old2)
        old, new = old2, new2
    if n != 1:
        raise SystemExit("мутация не применима: %s в %s найдено %d раз" % (old, rel, n))
    with open(p, "w", encoding="utf-8", newline="") as f:
        f.write(text.replace(old, new))


def run_tests(src):
    """Прогоняет pytest на копии src/ и возвращает (rc, список упавших тестов, строка итога, секунды)."""
    env = dict(os.environ, ISGUI_SRC=str(src), QT_QPA_PLATFORM="offscreen", PYTHONIOENCODING="utf-8")
    cmd = [sys.executable, "-X", "utf8", "-m", "pytest", str(TEST), "-q", "-p", "no:cacheprovider", "-rf"]
    t0 = time.monotonic()
    proc = subprocess.run(
        cmd, cwd=str(TK), env=env, capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=900
    )
    out = proc.stdout + proc.stderr
    failed = sorted(
        set(re.findall(r"^FAILED \S+::(\w+)", out, re.M))
        | set(re.findall(r"^ERROR \S+::(\w+)", out, re.M))
    )
    summary = ""
    for line in reversed(out.splitlines()):
        if re.search(r"\d+ (passed|failed)", line):
            summary = line.strip()
            break
    return proc.returncode, failed, summary, time.monotonic() - t0


def main():
    """Основная логика: базовый прогон, прогон с каждой мутацией, запись отчёта."""
    report = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else TK / "audit" / "gui_4detectors_mutations.md"
    rows = []

    with tempfile.TemporaryDirectory(prefix="isgui_mut_") as td:
        td_path = pathlib.Path(td)

        # --- Базовый прогон (исправная версия) ---
        src = make_copy(td_path / "base")
        rc, failed, summ, sec = run_tests(src)
        print("[BASE] rc=%d %s (%.0f с)" % (rc, summ, sec))
        base_summ = summ
        if rc != 0:
            print("исправная версия красная -- мутации не прогоняются")
            for name in failed:
                print("  ", name)
            sys.exit(2)

        # --- Прогон с каждой мутацией ---
        for name, purpose, rel, old, new in MUTS:
            src = make_copy(td_path / name)
            apply(src, rel, old, new)
            rc, failed, summ, sec = run_tests(src)
            killed = rc != 0 and len(failed) > 0
            rows.append((name, purpose, rel, killed, failed, summ))
            print("[%s] %s: %s; покраснели: %s" % (
                name, "УБИТА" if killed else "ВЫЖИЛА", summ, ", ".join(failed) or "-"))

    # --- Запись отчёта в Markdown ---
    report.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    lines.append("# Мутационная поверка тестов окна: 4 детектора\n\n")
    lines.append("Тесты: `apps/intrinsic-shape-gui/tests/test_gui_detectors.py`; прогон: `python apps/intrinsic-shape-gui/tests/gui_mutations.py`.\n\n")
    lines.append("Исправная версия: %s\n\n" % base_summ)
    lines.append("| Мутация | Что портит | Файл | Итог | Покраснели тесты |\n")
    lines.append("|---|---|---|---|---|\n")
    for name, purpose, rel, killed, failed, summ in rows:
        verdict = "убита" if killed else "ВЫЖИЛА"
        failed_str = ", ".join("`%s`" % t for t in failed) if failed else "—"
        lines.append("| %s | %s | %s | %s | %s |\n" % (name, purpose, rel, verdict, failed_str))
    killed_count = sum(1 for r in rows if r[3])
    total = len(rows)
    final_line = "Итого: мутаций %d, убито %d, выжило %d.\n" % (total, killed_count, total - killed_count)
    lines.append("\n" + final_line)

    with open(report, "w", encoding="utf-8") as f:
        f.write("".join(lines))

    print(final_line)
    sys.exit(0 if killed_count == total else 1)


if __name__ == "__main__":
    main()
