"""
GUI-обработчики для пары спектров «фон + образец» детектора LaBr3.
Заменяют три обработчика isgui.backend, сохраняя ключи результата и добавляя новые.
Вся вычислительная работа делегируется модулю gamma.calibration.intrinsic_shape.pipeline.
Сгенерировано gen_code.py (qwen3.8:27b) по specs/backend_pair.md, 2026-10-10
"""

import re
import json
import pathlib
import datetime
import numpy as np
from .backend import UserError, capture, pair_info, parse_interp, explain, file_check


def parse_manual(text):
    """
    Разбирает строку с ручными опорами.
    Формат: "a37 sum" или "a37 sum edge".
    Возвращает словарь или None, если текст пуст.
    """
    if not text or not text.strip():
        return None

    # Разбиваем по пробелам, запятым и точкам с запятой
    parts = re.split(r'[\s,;]+', text.strip())
    parts = [p for p in parts if p]

    if len(parts) not in (2, 3):
        raise UserError(
            "Ручные опоры: нужно 2 или 3 числа -- каналы фона пика ~37 кэВ, "
            "сумм-пика La-138 и (по желанию) β-края, например «110.4 3183.7»"
        )

    try:
        vals = [float(p) for p in parts]
    except ValueError:
        raise UserError(
            "Ручные опоры: нужно 2 или 3 числа -- каналы фона пика ~37 кэВ, "
            "сумм-пика La-138 и (по желанию) β-края, например «110.4 3183.7»"
        )

    a37, s = vals[0], vals[1]
    edge = vals[2] if len(vals) == 3 else None

    # Проверка положительных значений и порядка
    if a37 <= 0 or s <= 0:
        raise UserError(
            "Ручные опоры: нужно 2 или 3 числа -- каналы фона пика ~37 кэВ, "
            "сумм-пика La-138 и (по желанию) β-края, например «110.4 3183.7»"
        )

    if a37 >= s:
        raise UserError(
            "Ручные опоры: нужно 2 или 3 числа -- каналы фона пика ~37 кэВ, "
            "сумм-пика La-138 и (по желанию) β-края, например «110.4 3183.7»"
        )

    if edge is not None:
        if edge <= 0 or a37 >= edge or edge >= s:
            raise UserError(
                "Ручные опоры: нужно 2 или 3 числа -- каналы фона пика ~37 кэВ, "
                "сумм-пика La-138 и (по желанию) β-края, например «110.4 3183.7»"
            )

    result = {"a37": a37, "sum": s}
    if edge is not None:
        result["edge"] = edge

    return result


def run_register(bg_path, sample_path, manual=None, robust="auto", out_json=None, log=None):
    """
    Выполняет совместную регистрацию фона и образца.
    Возвращает словарь с результатами регистрации и информацией для GUI.
    """
    # Валидация файлов
    pair_info(bg_path, sample_path)

    # Импортируем toolkit внутри функции
    from gamma.calibration.intrinsic_shape import pipeline as pl
    from gamma.calibration.intrinsic_shape.core import poly_eval

    # Запускаем регистрацию
    res, text = capture(
        pl.run_register,
        bg_path,
        sample_path,
        manual=manual,
        robust=robust,
        profiles=False,
        out_json=out_json,
        log=log
    )

    out = res["out"]
    items = parse_interp(text)
    c = out["c"]
    anchor = out["anchor"]
    anchor_z = (poly_eval(c, anchor[0]) - anchor[1]) / anchor[2]

    return {
        "h": res["h"],
        "perr": out["perr"],
        "chi2": out["chi2"],
        "ndf": out["ndf"],
        "chi2ndf": out["chi2ndf"],
        "k": out["k"],
        "n_minima": len(out["minima"]),
        "n_close": out["spread"]["n"],
        "grid_on_edge": out["grid_on_edge"],
        "anchor_z": float(anchor_z),
        "interp": items,
        "interp_explained": [(it, explain(it)) for it in items],
        "h_json": out,
        "auto": out["auto"],
        "robust": {k: v for k, v in out["robust"].items() if k != "model"},
        "structure": out["structure"],
        "structure_plain": out.get("structure_plain", []),
        "structure_summary": out.get("structure_summary"),
        "nested": res["ctx"]["nested"],
        "n_ch": res["ctx"]["n"],
        "log": text
    }


def run_ecal(bg_path, sample_path, h, manual=None, source=None, source_nuclide="Am-241",
             e37=None, use37=True, out_json=None, log=None, reg_json=None):
    """
    Строит энергетическую шкалу по паре спектров.
    Возвращает словарь с коэффициентами шкалы и информацией для GUI.
    """
    if not isinstance(h, list) or len(h) != 3:
        raise UserError("Сначала выполните «Совместить»: нет коэффициентов h")

    # Импортируем toolkit внутри функции
    from gamma.calibration.intrinsic_shape import pipeline as pl

    # Запускаем калибровку
    res, text = capture(
        pl.run_ecal,
        bg_path,
        sample_path,
        list(h),
        manual=manual,
        source=source or None,
        source_nuclide=source_nuclide,
        e37=e37,
        use37=use37,
        out_json=out_json,
        robust_json=reg_json,
        log=log
    )

    result = dict(res)
    result["checks"] = [line.strip() for line in text.splitlines() if line.startswith("[CHK]")]
    result["log"] = text

    # Формируем строки таблицы для GUI
    rows = []
    for p in res["points"]:
        rows.append((
            p["name"],
            "основная",
            float(p["ch"]),
            float(p["ch_err"]),
            float(p["E"]),
            float(p["residual"]),
            float(p["loo"])
        ))

    if "low" in res:
        for r in res["low"]["residuals"]:
            loo_val = res["low"]["loo"].get(r["name"], float('nan'))
            rows.append((
                r["name"],
                "ниже 100 кэВ",
                float(r["ch"]),
                float('nan'),
                float(r["E"]),
                float(r["after"]),
                float(loo_val)
            ))

    result["branch_rows"] = rows

    return result


def run_subtract(bg_path, sample_path, h, ecal, fw_points, out_zeroed, out_signed,
                 name=None, date=None, monotone_to=1.5, sum_window=None, z_limit=2.0,
                 reg_json=None, manual=None, xml_scale="composite", log=None, lt_out="sample"):
    """
    Выполняет вычитание фона из образца с проверками.
    Возвращает словарь с результатами вычитания и информацией для GUI.
    """
    # Валидация файлов
    pair = pair_info(bg_path, sample_path)

    # Проверки входных данных
    if not isinstance(h, list) or len(h) != 3:
        raise UserError("Сначала выполните «Совместить»: нет коэффициентов h")

    if not ecal:
        raise UserError("Сначала выполните «Шкала»: нет энергетической шкалы")

    if not out_zeroed or not out_signed:
        raise UserError("Укажите файлы для записи результата")

    out_zeroed_path = pathlib.Path(out_zeroed).resolve()
    out_signed_path = pathlib.Path(out_signed).resolve()
    bg_path_res = pathlib.Path(bg_path).resolve()
    sample_path_res = pathlib.Path(sample_path).resolve()

    if out_zeroed_path == out_signed_path:
        raise UserError("Файлы результата не должны совпадать")

    if out_zeroed_path == bg_path_res or out_zeroed_path == sample_path_res or \
       out_signed_path == bg_path_res or out_signed_path == sample_path_res:
        raise UserError("Нельзя записывать результат поверх исходных файлов")

    # Дата по умолчанию
    if date is None:
        date = datetime.date.today().isoformat()

    # Импортируем toolkit внутри функции
    from gamma.calibration.intrinsic_shape import pipeline as pl
    from gamma.calibration.intrinsic_shape import fwhm as fwm
    from gamma.calibration.intrinsic_shape.core import poly_eval, rebin_conserve_var
    from gamma.calibration.intrinsic_shape import register as reg

    # Запускаем вычитание
    res, text = capture(
        pl.run_subtract,
        bg_path,
        sample_path,
        list(h),
        ecal,
        list(fw_points),
        out_zeroed,
        out_signed,
        name=name,
        date=date,
        monotone_to=monotone_to,
        sum_window=sum_window,
        manual=manual,
        reg_json=reg_json,
        xml_scale=xml_scale,
        lt_out=lt_out,
        log=log
    )

    # Вычисления для GUI
    resid = np.asarray(res["resid"], float)
    var = np.asarray(res["var"], float)
    lo, hi = res["sum_window"]
    tot = resid[lo:hi+1].sum()
    sig = np.sqrt(var[lo:hi+1].sum())
    z = tot / sig if sig > 0 else 0.0

    # Проверка файлов
    files = [
        file_check(out_zeroed, res["zeroed"], False),
        file_check(out_signed, res["signed"], True)
    ]

    # Проверка ПШПВ
    fw_ok, fw_vals, fw_vertex, fw_margin = fwm.check(res["fwhm"], len(resid))

    # Вычисление энергетической шкалы и данных для графика
    S, B, k, _, _ = reg.load_inputs(bg_path, sample_path)
    edges = poly_eval(list(h), np.arange(len(S) + 1) - 0.5)
    s_on_bg, _ = rebin_conserve_var(np.asarray(S, float), edges, len(B))
    kb = s_on_bg - resid
    energy = poly_eval(list(res["coeffs_xml"]), np.arange(len(B), dtype=float))

    # Проверка структур
    # правка Claude (раунд 4): вердикт окон -- verdict.summarize (статус «не оценено» при перекрытии линиями образца,
    # двусторонний z_χ²); несошедшаяся устойчивая регистрация -- all_ok=False
    checks = res["checks"]
    summ = res["checks_summary"]
    computed = [c for c in checks if c["evaluated"]]
    struct_ok = summ["struct_ok"]
    converged = bool(res.get("converged", True))

    # Проверка файлов
    files_ok = all(
        f["match"] and f["all_int"] and f["has_fwhm"] and f["min_ok"]
        for f in files
    )

    all_ok = bool(fw_ok and files_ok and converged)   # правка Claude: вердикт окон справочный, в all_ok не входит

    return {
        "k": float(res["k"]),
        "sum": float(tot),
        "sigma": float(sig),
        "z": float(z),
        "z_ok": abs(z) <= z_limit,
        "z_limit": z_limit,
        "z_info_only": True,
        "files": files,
        "fwhm": [float(v) for v in res["fwhm"]],
        "fwhm_ok": bool(fw_ok),
        "fwhm_values": fw_vals,
        "fwhm_margin": float(fw_margin),
        "struct": checks,
        "struct_ok": bool(struct_ok),
        "n_struct_computed": len(computed),
        "struct_summary": summ,
        "converged": converged,
        "lt_out": res.get("lt_out", "sample"),
        "live_s": float(res.get("live_s", 0.0)),
        "live_b": float(res.get("live_b", 0.0)),
        "live_out": float(res.get("live_out", 0.0)),
        "excluded": res["excluded"],
        "all_ok": all_ok,
        "channels": np.arange(len(B)),
        "energy": energy,
        "sample": s_on_bg,
        "bg_k": kb,
        "resid": resid,
        "zeroed": np.asarray(res["zeroed"]),
        "signed": np.asarray(res["signed"]),
        "n_negative": int((resid < 0).sum()),
        "coeffs_xml": list(res["coeffs_xml"]),
        "sum_window": list(res["sum_window"]),
        "log": text
    }


def _k40_row(c):
    k = c.get("k40")
    if not k:
        return []
    return ['<tr><td colspan="8">K-40 в окне: площадь %.0f ± %.0f отсчётов, σ %.2f кан</td></tr>' % (
        k["area"], k["area_err"], k["sigma"])]


def struct_table_html(checks, title="Окна собственных структур"):
    """Правка Claude (раунд 4): таблица окон -- энергии, доля каналов, скрытых линиями образца, статус «не оценено»."""
    import html
    rows = ["<h3>%s</h3>" % html.escape(title),
            '<p style="color:#555"><i>справочно: самопроверка неполная, см. ограничения</i></p>', '<table border="1" cellpadding="3">',
            "<tr><th>Окно</th><th>кэВ</th><th>Каналы</th><th>Скрыто линиями образца</th><th>z суммы</th>"
            "<th>z диполя</th><th>z χ²</th><th>Итог</th></tr>"]
    for c in checks:
        if not c.get("evaluated"):
            st, col = "не оценено", "#b8860b"   # жёлтый: окно скрыто, а не «в норме»
        else:
            st, col = ("в норме", "#333") if c["ok"] else ("вне нормы", "#c62828")   # нейтральный, не зелёный
        why = (" (" + c["reason"] + ")") if c.get("reason") else ""
        rows.append("<tr><td>%s</td><td>%.0f–%.0f</td><td>%d–%d</td><td>%.0f %%</td><td>%+.2f</td><td>%+.2f</td>"
                    '<td>%+.2f</td><td style="color:%s">%s</td></tr>' % (
                        html.escape(c["name"]), c["e_lo"], c["e_hi"], c["lo"], c["hi"], 100.0 * c["hidden_frac"],
                        c["z_sum"], c["z_dip"], c["z_chi2"], col, html.escape(st + why)))
        rows.extend(_k40_row(c))
    rows.append("</table>")
    return "\n".join(rows)


def register_html(r):
    """
    Формирует HTML-фрагмент для отображения результатов регистрации.
    """
    import html

    parts = []
    parts.append("<h3>Опоры автопресета</h3>")
    parts.append('<table border="1" cellpadding="3">')
    parts.append("<tr><th>Опора</th><th>Канал фона</th><th>±</th><th>Статус</th></tr>")

    anchors = r["auto"]["anchors"]
    anchor_names = {
        "a37": "пик ~37 кэВ",
        "sum": "сумм-пик La-138",
        "edge": "β-край La-138"
    }

    for key in ["a37", "sum", "edge"]:
        if key in anchors:
            a = anchors[key]
            name = anchor_names[key]
            ch = a["x"]
            err = a["err"]
            ok = a["ok"]
            manual = a["manual"]

            if manual:
                status = "вручную"
            elif ok:
                status = "уверенно"
            else:
                status = "НЕУВЕРЕННО"

            parts.append(
                f"<tr><td>{html.escape(name)}</td>"
                f"<td>{ch:.2f}</td>"
                f"<td>{err:.2f}</td>"
                f"<td>{html.escape(status)}</td></tr>"
            )

    parts.append("</table>")

    # Отображение каналов
    n_ch = r["n_ch"]
    auto = r["auto"]
    if auto["identity"]:
        map_text = "тождественное -- пресеты без пересчёта"
    else:
        m = auto["map"]
        map_text = f"x = {m[0]:.5g} + {m[1]:.7g}·x + {m[2]:.3g}·x²"

    parts.append(
        f"<p>Отображение каналов эталона (2048) на прибор ({n_ch} каналов): "
        f"{html.escape(map_text)}, масштаб {auto['r']:.5f}</p>"
    )

    # Предупреждения
    for w in auto.get("warnings", []):
        parts.append(f'<p style="color:#c62828">{html.escape(w)}</p>')

    # Вложенность
    if r["nested"]:
        parts.append(
            f'<p style="color:#c62828"><b>Внимание:</b> {html.escape(r["nested"])}</p>'
        )

    # Линии образца
    parts.append("<h3>Линии образца</h3>")
    robust = r["robust"]
    if robust["mode"] == "plain":
        parts.append("<p>простая регистрация (образец не добавляет структур)</p>")
    else:
        parts.append("<p>устойчивая регистрация</p>")

    parts.append(
        f"<p>χ²/ndf простой регистрации {robust['chi2ndf_plain']:.3f} → "
        f"итог {robust['chi2ndf']:.3f} (разные каналы и модели -- справочно)</p>"
    )
    sb = robust.get("same_basis")   # правка Claude (раунд 4, п.6): сравнение на одних каналах и одной модели
    if sb:
        parts.append(f"<p>На одних каналах ({sb['n_ch']}): без модели образца χ²/ndf простая {sb['no_model'][0]:.3f} → "
                     f"устойчивая {sb['no_model'][1]:.3f}; с моделью {sb['with_model'][0]:.3f} → {sb['with_model'][1]:.3f}</p>")
    if robust["mode"] == "robust" and not robust.get("converged", True):
        parts.append('<p style="background:#fff3b0; color:#6b4e00"><b>Устойчивая регистрация не сошлась</b>: образец с '
                     'множеством линий по всей шкале%s — итог не принимается</p>' % (
                         (" (" + html.escape(robust.get("stop_reason")) + ")") if robust.get("stop_reason") else ""))

    excluded = robust.get("excluded", [])
    if excluded:
        excl_text = ", ".join(
            f"кан {e[0]:d}–{e[1]:d} ({e[2]:.0f}–{e[3]:.0f} кэВ)"
            for e in excluded
        )
        parts.append(f"<p>Исключённые окна: {html.escape(excl_text)}</p>")
    else:
        parts.append("<p>окна не исключались</p>")

    nuclides = robust.get("nuclides", [])
    if nuclides:
        nucl_text = ", ".join(nuclides)
        parts.append(f"<p>Нуклиды: {html.escape(nucl_text)}</p>")

    k40 = robust.get("k40")
    if k40:
        parts.append(
            f"<p>площадь K-40 образца под сумм-пиком: "
            f"{k40['area']:.0f} ± {k40['area_err']:.0f} отсчётов</p>"
        )

    parts.append(struct_table_html(r["structure"], "Окна собственных структур (итоговая регистрация)"))
    if r.get("structure_plain") and r["robust"]["mode"] == "robust":
        parts.append(struct_table_html(r["structure_plain"], "Те же окна, простая регистрация (та же основа)"))

    return "\n".join(parts)


def ecal_html(r):
    """
    Формирует HTML-фрагмент для отображения результатов калибровки.
    """
    import html

    parts = []

    # Основная шкала
    coeffs = r["coeffs"]   # правка Claude: любая степень; отвергнутая шкала -- красной строкой, в XML калибровка файла
    scale_text = "Основная шкала E(x) = " + " + ".join("%.6g·x^%d" % (c, i) if i else "%.6g" % c for i, c in enumerate(coeffs))
    if r.get("rejected"):
        parts.append('<p style="color:#c62828"><b>Шкала по реперам отвергнута:</b> %s. В XML остаётся калибровка '
                     'файла фона (показана ниже).</p>' % html.escape(r["rejected"]))
    parts.append(f"<p>{html.escape(scale_text)}</p>")
    parts.append(
        f"<p>χ²/ndf = {r['chi2ndf']:.2f}, rms невязок {r['rms']:.2f} кэВ, "
        f"LOO rms {r['loo_rms']:.2f} кэВ</p>"
    )

    # Реперы образца
    if r["lu"]:
        parts.append("<p>Реперы образца: из пресета Lu-176</p>")
    else:
        lines = r.get("lines", [])
        if lines:
            lines_text = ", ".join(f"{l['nuclide']} {l['E']:.1f} кэВ" for l in lines)   # правка Claude: ключ nuclide
            parts.append(
                f"<p>Реперы образца: линии, найденные в разности: "
                f"{html.escape(lines_text)}</p>"
            )
        else:
            parts.append(
                "<p>линий образца не найдено -- только собственные структуры фона</p>"
            )

    # Шкала ниже 100 кэВ
    if "low" in r:
        low = r["low"]
        parts.append("<h3>Шкала ниже 100 кэВ</h3>")
        parts.append(
            f"<p>s = {low['s']:.4f} ± {low['sigma_s']:.4f} (форма {html.escape(low['shape'])}), "
            f"ступень на K-крае La {low['dK']:.3f}, "
            f"гамма-эквивалент пика ~37 кэВ {low['e37']:.2f} кэВ</p>"
        )
        parts.append(f"<p>порог прибора ≈ {low['e_threshold']:.1f} кэВ</p>")

        if low.get("warn_extrapolation"):
            parts.append(
                '<p style="color:#c62828">Ниже 37 кэВ — экстраполяция: '
                'нужен репер Am-241 / Ba-133</p>'
            )

        xml = low.get("xml", {})
        if xml.get("fallback"):
            dev = xml.get("damage", float("nan"))   # правка Claude: порча шкалы выше 100 кэВ
            parts.append(
                f"<p>В XML записывается: основная шкала "
                f"(составная портит шкалу выше 100 кэВ на {dev:.2f} кэВ)</p>"
            )
        else:
            deg = xml.get("deg", 2)
            max_dev = xml.get("max_dev", {})
            if max_dev:
                dev_text = ", ".join(f"{k}: {v:.2f} кэВ" for k, v in max_dev.items())
                parts.append(
                    f"<p>В XML записывается: полином степени {deg} "
                    f"(макс. отклонение от модели: {html.escape(dev_text)})</p>"
                )
            else:
                parts.append(
                    f"<p>В XML записывается: полином степени {deg}</p>"
                )

    # Строки [CHK]
    for line in r.get("checks", []):
        parts.append(f"<p>{html.escape(line)}</p>")

    return "\n".join(parts)


def subtract_html(r):
    """
    Формирует HTML-фрагмент для отображения результатов вычитания.
    """
    import html

    def status(ok, text):
        color = "#1a7f37" if ok else "#c62828"
        symbol = "✔" if ok else "✘"
        return f'<span style="color:{color}; font-weight:bold">{symbol} {html.escape(text)}</span>'

    parts = []

    # Общий статус
    if r["all_ok"]:
        parts.append(status(True, "Файлы, ПШПВ и сходимость в порядке"))
    else:
        parts.append(status(False, "Есть проблемы"))

    parts.append(f"<p>k = {r['k']:.6f}</p>")
    parts.append(f"<p>Разность в cps по каналам: net_cps = S'/LT_s − B̃/LT_b, LT_s = {r['live_s']:.1f} с, LT_b = {r['live_b']:.1f} с; "
                 f"в XML записано net_cps·LT_out (целые), LiveTime файла = {r['live_out']:.1f} с "
                 f"({'образца' if r['lt_out'] == 'sample' else 'фона'}) — BecqMoni показывает ровно разность cps</p>")

    parts.append(struct_table_html(r["struct"]))
    sm_ = r.get("struct_summary") or {}
    parts.append('<p style="color:#555">Итог окон (справочно, в общий итог не входит): оценено %s из %s, в норме %s, не оценено %s</p>' % (sm_.get("n_eval"), sm_.get("n"), sm_.get("n_ok"), sm_.get("n_not_eval")))

    # Остаток разности
    parts.append(
        f"<p>Остаток разности по всей области (справочно, у образца со своими линиями "
        f"не равен нулю): {r['sum']:.1f} ± {r['sigma']:.1f} (z = {r['z']:+.2f})</p>"
    )

    # ПШПВ
    parts.append(
        f"<p>Модель ПШПВ годна для BecqMoni (запас за последним каналом "
        f"{r['fwhm_margin']:.0f} кан.)</p>"
    )

    # Файлы
    for f in r["files"]:
        all_int = "да" if f["all_int"] else "нет"
        match = "да" if f["match"] else "нет"
        has_fwhm = "есть" if f["has_fwhm"] else "нет"
        sum_text = f"{f['sum']:.0f}" if f["sum"] is not None else "N/A"

        parts.append(
            f"<p>{html.escape(f['name'])}: точек {f['n_points']}, "
            f"все целые: {all_int}, совпадает после перечитывания: {match}, "
            f"блок ПШПВ: {has_fwhm}, сумма {sum_text}</p>"
        )
        parts.append(f'<p style="font-size:small">{html.escape(f["path"])}</p>')

    # Отрицательные каналы
    parts.append(f"<p>Отрицательных каналов в разности до обнуления: {r['n_negative']}</p>")

    return "\n".join(parts)


def markers_from(ecal_result):
    """
    Извлекает маркеры из результата калибровки.
    Возвращает список словарей {"name", "ch", "E"}.
    """
    markers = []
    seen_names = set()

    # Основные точки
    for p in ecal_result.get("points", []):
        markers.append({
            "name": p["name"],
            "ch": float(p["ch"]),
            "E": float(p["E"])
        })
        seen_names.add(p["name"])

    # Точки низкой энергии
    if "low" in ecal_result:
        for r in ecal_result["low"].get("residuals", []):
            if r["name"] not in seen_names:
                markers.append({
                    "name": r["name"],
                    "ch": float(r["ch"]),
                    "E": float(r["E"])
                })
                seen_names.add(r["name"])

    return markers

