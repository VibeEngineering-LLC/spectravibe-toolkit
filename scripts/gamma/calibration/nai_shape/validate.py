"""Модуль валидации и суррогатного тестирования формы НАИ (пакет gamma.calibration.nai_shape)."""

import numpy as np

from gamma.calibration.intrinsic_shape.core import load_spectrum, poly_eval
from . import pipeline, surrogate


def segment(path_b, path_a=None):
    """Извлечь независимый отрезок вложенного (кумулятивного) накопления.

    Если задан path_a — вычитается из path_b. Отрицательные каналы означают
    ошибку вложенности.
    """
    B = load_spectrum(path_b)
    if path_a is not None:
        A = load_spectrum(path_a)
        n = len(B["counts"])
        a_counts = A["counts"][:n]
        if len(a_counts) < n:
            a_counts = np.concatenate([a_counts, np.zeros(n - len(a_counts))])
        B["counts"] = B["counts"] - a_counts
        B["live"] -= A["live"]
        B["real"] -= A["real"]
        B["name"] = B["name"] + " минус " + A["name"]
        if np.any(B["counts"] < 0):
            raise ValueError("файлы не вложены")
    return B


def nested_status(B, A, tol_s=1.0):
    """Ревизия 6 (criteria_v6.md П6.5): A вложен в B по времени и поканально. Возвращает (bool, причина)."""
    import datetime as _dt
    try:
        sa, sb = (_dt.datetime.fromisoformat(str(X["start"])) for X in (A, B))
    except ValueError:
        return False, "невложенные: нет времени начала"
    if (sa - sb).total_seconds() < -tol_s:
        return False, "невложенные: A начат раньше B"
    if (sa - sb).total_seconds() + A["real"] > B["real"] + tol_s:
        return False, "невложенные: A закончен позже B"
    n = min(len(A["counts"]), len(B["counts"]))
    if np.any(np.asarray(A["counts"][:n]) > np.asarray(B["counts"][:n])):
        return False, "невложенные: поканально A > B"
    return True, ""


def surrogate_cases(bg, profile, cases, seeds=(1, 2, 3), frac=0.5, thin=1.0, free=("d", "gN")):
    """Прогнать набор суррогатных случаев на реальном фоне и собрать расхождения.

    Параметры:
        bg       — dict спектра фона (результат load_spectrum)
        profile  — профиль тестирования (содержит "test_kev")
        cases    — список dict с ключами "d", "g", необязательно "c2", "sigma_extra"
        seeds    — кортеж/список целых сидов
        frac     — доля каналов, отводимая под шаблон
        thin     — коэффициент прореживания
        free     — кортеж свободных параметров

    Возвращает список dict-строк с результатами.
    """
    n = len(bg["counts"])
    ref = bg["counts"].copy()
    ref[-1] = 0

    rows = []
    for ci, case in enumerate(cases):
        c_true = [float(case["d"]), float(case["g"]), float(case.get("c2", 0.0))]
        for seed in seeds:
            # ручная правка (VERIFY.md п. 3): независимый поток на каждую пару (случай, реализация);
            # было default_rng(seed) -- одни и те же шаблон и события во всех случаях
            rng = np.random.default_rng([ci, int(seed)])
            tmpl, base = surrogate.split_poisson(ref, frac, rng)
            smp = surrogate.drift_events(
                base, c_true, n, rng,
                thin=thin, sigma_extra=float(case.get("sigma_extra", 0.0)),
            )
            B = dict(bg, counts=tmpl, live=bg["live"] * frac, real=bg["real"] * frac, name="шаблон")
            S = dict(bg, counts=smp, live=bg["live"] * (1 - frac) * thin,
                     real=bg["real"] * (1 - frac) * thin, name="суррогат")
            r = pipeline._run_pair(B, S, profile, mode="background",
                                  free=list(free), verbose=False)

            xs = np.array(r["x_test"])
            h_true = poly_eval(c_true, xs)
            h_fit = poly_eval(r["c"], xs)
            sig = np.array(r["sigma_h_test"])
            pull = (h_fit - h_true) / sig

            rows.append({
                "case": case,
                "seed": seed,
                "test_kev": list(profile["test_kev"]),
                "h_true": [float(v) for v in h_true],
                "h_fit": [float(v) for v in h_fit],
                "sigma": [float(v) for v in sig],
                "pull": [float(v) for v in pull],
                "c_true": c_true,
                "c_fit": [float(v) for v in r["c"]],
                "chi2ndf": float(r["chi2ndf"]),
            })
    return rows


def surrogate_summary(rows):
    """Сводная статистика по pull-расхождениям всех суррогатных строк.

    Возвращает dict с n, max_abs_pull, rms_pull, rms/mean по контрольным точкам.
    """
    pulls = np.array([r["pull"] for r in rows])
    n = len(rows)
    n_points = pulls.shape[1] if pulls.ndim > 1 else 1

    rms_by_point = [float(np.sqrt(np.mean(pulls[:, i] ** 2))) for i in range(n_points)]
    mean_by_point = [float(np.mean(pulls[:, i])) for i in range(n_points)]

    return {
        "n": n,
        "max_abs_pull": float(np.max(np.abs(pulls))),
        "rms_pull": float(np.sqrt(np.mean(pulls ** 2))),
        "rms_pull_by_point": rms_by_point,
        "mean_pull_by_point": mean_by_point,
    }


def null_pair(bg, other, profile, free=("d", "gN")):
    """Нулевая пара: сравнение двух реальных спектров в режиме фона.

    Возвращает dict с именем пар, коэффициентами, chi2ndf, тестами и сдвигом.
    """
    r = pipeline._run_pair(bg, other, profile, mode="background",
                          free=list(free), verbose=False)
    shift_test = [float(h - x) for h, x in zip(r["h_test"], r["x_test"])]
    return {
        "bg": bg["name"],
        "other": other["name"],
        "c": [float(v) for v in r["c"]],
        "chi2ndf": float(r["chi2ndf"]),
        "tests": r["tests"],
        "verdict": r["verdict"],
        # ревизия 6: внутренний результат по умолчанию не проверен; флаг снимает только api.run
        "scale_verified": False,
        "flag": "шкала не проверена",
        "status": "experimental",
        "shift_test": shift_test,
        "sigma_h_test": [float(v) for v in r["sigma_h_test"]],
    }


def format_surrogate_table(rows):
    """Сформировать Markdown-таблицу суррогатных результатов по контрольным точкам."""
    header = "| случай (d, g) | seed | E, кэВ | h задано | h получено | σ | расхождение, σ |"
    sep = "|---|---|---|---|---|---|---|"
    lines = [header, sep]
    for row in rows:
        case = row["case"]
        case_label = "d={0:g}, g={1:g}".format(case["d"], case["g"])
        seed = row["seed"]
        kev = row["test_kev"]
        for i in range(len(kev)):
            line = "| {0} | {1} | {2:.1f} | {3:.2f} | {4:.2f} | {5:.3f} | {6:+.2f} |".format(
                case_label,
                seed,
                kev[i],
                row["h_true"][i],
                row["h_fit"][i],
                row["sigma"][i],
                row["pull"][i],
            )
            lines.append(line)
    return "\n".join(lines)


def format_null_table(items):
    """Сформировать Markdown-таблицу тестовых окон для нулевых пар."""
    header = "| пара | E, кэВ | метка | z_sum | z_dip | исключено |"
    sep = "|---|---|---|---|---|---|"
    lines = [header, sep]
    for item in items:
        pair = "{0} / {1}".format(item["bg"][:40], item["other"][:40])
        for t in item["tests"]:
            line = "| {0} | {1:.1f} | {2} | {3:+.2f} | {4:+.2f} | {5} |".format(
                pair,
                t["E"],
                t["tag"],
                t["z_sum"],
                t["z_dip"],
                "да" if t["excluded"] else "нет",
            )
            lines.append(line)
    return "\n".join(lines)
