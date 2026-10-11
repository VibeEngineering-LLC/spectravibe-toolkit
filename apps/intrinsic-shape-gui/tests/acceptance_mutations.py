# Мутационная поверка приёмки: для каждой точечной мутации прогоняет tests/acceptance.py на нужном подмножестве шагов и показывает, что
# краснеет ровно целевая проверка, а исправная версия зелёная. Сгенерировано gen_code.py (qwen3.8:27b) по спецификации,
# 2026-10-10

import os
import re
import sys
import json
import time
import pathlib
import subprocess

sys.stdout.reconfigure(encoding="utf-8")

ROOT = pathlib.Path(__file__).resolve().parents[1]

MUTS = [("", "FGH", None),
        ("F3_k", "F", "F3"), ("F6_onesided", "F", "F6"), ("F7_nested", "F", "F7"), ("F8_hide", "F", "F8"),
        ("G2_plain", "G", "G2"), ("G3_maskon", "G", "G3"), ("G5_lu", "G", "G5"), ("H3_stale", "H", "H3")]


def run_one(mut, only, log_dir):
    """
    Запускает tests/acceptance.py с заданной мутацией и подмножеством шагов.
    Возвращает словарь результатов проверок и строку критериев.
    """
    env = dict(os.environ,
               ISGUI_ACC_MUT=mut,
               ISGUI_ACC_ONLY=only,
               PYTHONIOENCODING="utf-8",
               PYTHONUNBUFFERED="1",
               QT_QPA_PLATFORM="offscreen")

    try:
        proc = subprocess.run(
            [sys.executable, str(ROOT / "tests" / "acceptance.py")],
            env=env,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=5400
        )
        stdout = proc.stdout
        stderr = proc.stderr
    except subprocess.TimeoutExpired:
        stdout = ""
        stderr = "TIMEOUT"

    # Запись полного вывода в лог
    log_file = log_dir / ("mut_%s.log" % (mut or "clean"))
    with open(log_file, "w", encoding="utf-8") as f:
        f.write("=== STDOUT ===\n")
        f.write(stdout)
        f.write("\n=== STDERR ===\n")
        f.write(stderr)

    # Разбор результатов
    results = {}
    criteria_line = ""

    all_lines = (stdout + "\n" + stderr).splitlines()
    for line in all_lines:
        if line.startswith("[OK] "):
            # Формат: [OK] Имя проверки: описание
            rest = line[5:]  # после "[OK] "
            # Имя проверки — текст до первого ": "
            if ": " in rest:
                name = rest.split(": ")[0].strip()
                results[name] = True
        elif line.startswith("[FAIL] "):
            rest = line[6:]  # после "[FAIL] "
            if ": " in rest:
                name = rest.split(": ")[0].strip()
                results[name] = False
        elif line.startswith("[CRITERIA]"):
            criteria_line = line

    return results, criteria_line


def main():
    # Определение каталога логов
    if len(sys.argv) > 1:
        log_dir = pathlib.Path(sys.argv[1])
    else:
        log_dir = ROOT / "tests" / "mutlogs"
    log_dir.mkdir(parents=True, exist_ok=True)

    # Определение списка мутаций
    if len(sys.argv) > 2 and sys.argv[2]:
        selected_muts = [m.strip() for m in sys.argv[2].split(",") if m.strip()]
    else:
        selected_muts = [m[0] for m in MUTS]

    # Словарь для быстрого поиска параметров мутации
    mut_params = {m[0]: (m[1], m[2]) for m in MUTS}

    base = None
    all_criteria = ""
    mut_count = 0
    target_red_count = 0
    not_red_muts = []
    clean_failures = []

    for mut in selected_muts:
        if mut not in mut_params:
            print("[WARN] Неизвестная мутация: %s" % mut)
            continue

        only, prefix = mut_params[mut]
        print("[RUN] %s шаги %s" % (mut or "clean", only))
        start_time = time.time()

        res, criteria = run_one(mut, only, log_dir)
        elapsed = time.time() - start_time

        if criteria:
            all_criteria = criteria

        if mut == "":
            # Чистый прогон
            base = res
            failures = [n for n, ok in res.items() if not ok]
            fail_str = "; ".join(failures) if failures else "нет"
            print("[CLEAN] проверок %d, провалов %d: %s" % (len(res), len(failures), fail_str))
            clean_failures = failures
        else:
            # Мутация
            mut_count += 1
            # Поиск целевых проверок
            target = [n for n in res if n.startswith(prefix + " ")]
            t_red = target and all(not res[n] for n in target)

            # Определение прочих покрасневших
            if base is not None:
                others_red = [n for n, ok in res.items()
                              if not ok and not n.startswith(prefix + " ") and base.get(n, True)]
            else:
                # Если чистого прогона не было, считаем все прочие красные как "прочие"
                others_red = [n for n, ok in res.items()
                              if not ok and not n.startswith(prefix + " ")]

            others_str = "; ".join(others_red) if others_red else "нет"
            status = "КРАСНЕЕТ" if t_red else "НЕ КРАСНЕЕТ"

            if t_red:
                target_red_count += 1
            else:
                not_red_muts.append(mut)

            print("[MUT] %-12s цель %-4s %s; прочие покрасневшие: %s; %.0f с" %
                  (mut, prefix, status, others_str, elapsed))

    # Итог
    not_red_str = "; ".join(not_red_muts) if not_red_muts else "нет"
    print("[MUT-ИТОГ] мутаций %d, целевая краснеет у %d, не краснеет: %s" %
          (mut_count, target_red_count, not_red_str))

    if all_criteria:
        print(all_criteria)

    # Код возврата
    success = True
    if clean_failures:
        success = False
    if not_red_muts:
        success = False

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
