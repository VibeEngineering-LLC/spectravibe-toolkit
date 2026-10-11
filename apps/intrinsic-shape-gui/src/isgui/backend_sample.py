"""Обработчики режима «калибровка по пробе» (разность на шкале пробы): шкала и вычитание через
gamma.calibration.intrinsic_shape.samplecal_net; результат -- в формате вкладки «Вычесть» и графика.
Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-10-11 (ответ повторил шапку спецификации -- взят код; ключи
decision_html исправлены Claude); run_sample_cal и decision_html переписаны по VERIFY7
(gen_code.py, qwen3.8:27b, 2026-10-11)."""
import datetime
import pathlib
import numpy as np
from .backend import UserError, capture, pair_info, file_check


def run_sample_cal(bg_path, sample_path, repers, log=None, use_reper_scale=False):
    """Калибровка по пробе (VERIFY7): шкала XML по умолчанию -- калибровка файла пробы; use_reper_scale=True -- записать шкалу по
реперам (только если рекомендована, иначе ValueError из toolkit). branch_rows -- таблица «файл / реперы»."""
    pair_info(bg_path, sample_path)
    from gamma.calibration.intrinsic_shape import samplecal_net as sn
    res, text = capture(sn.run_sample_cal, bg_path, sample_path, list(repers), log=log, use_reper_scale=use_reper_scale)
    result = dict(res, log=text)
    d = res["decision"]
    result["branch_headers"] = ["Линия", "Роль", "Канал пробы", "E, кэВ", "Невязка: файл, кэВ", "Невязка: реперы, кэВ", "LOO, кэВ"]
    result["branch_rows"] = [(name, "репер", ch, E, rf, rn, loo) for name, E, ch, rf, rn, loo in d.get("rows_repers", [])] + \
                            [("%s %.1f" % (nuc, E), "независимая", ch, E, rf, rn, float("nan")) for nuc, E, ch, rf, rn in d.get("rows_ind", [])]
    return result


def run_sample_subtract(bg_path, sample_path, cal_res, fw_points, out_zeroed, out_signed, 
                        name=None, date=None, lt_out="sample", log=None):
    """
    Запуск вычитания фона по пробе.
    
    Args:
        bg_path: Путь к файлу фона.
        sample_path: Путь к файлу пробы.
        cal_res: Результат калибровки.
        fw_points: Точки для проверки FWHM.
        out_zeroed: Путь для записи результата с нулевым фоном.
        out_signed: Путь для записи знакового результата.
        name: Имя (опционально).
        date: Дата (опционально).
        lt_out: Тип вывода (по умолчанию "sample").
        log: Логгер (опционально).
        
    Returns:
        Словарь с результатами вычитания.
    """
    # Проверки путей
    if not out_zeroed or not out_signed:
        raise UserError("Укажите файлы для записи результата")
    
    # Проверка совпадения путей
    out_zeroed_path = pathlib.Path(out_zeroed).resolve()
    out_signed_path = pathlib.Path(out_signed).resolve()
    bg_path_resolved = pathlib.Path(bg_path).resolve()
    sample_path_resolved = pathlib.Path(sample_path).resolve()
    
    if out_zeroed_path == out_signed_path:
        raise UserError("Файлы результата не должны совпадать")
    
    if out_zeroed_path == bg_path_resolved or out_zeroed_path == sample_path_resolved:
        raise UserError("Нельзя записывать результат поверх исходных файлов")
    
    if out_signed_path == bg_path_resolved or out_signed_path == sample_path_resolved:
        raise UserError("Нельзя записывать результат поверх исходных файлов")
    
    date = date or datetime.date.today().isoformat()
    
    from gamma.calibration.intrinsic_shape import samplecal_net as sn
    from gamma.calibration.intrinsic_shape import fwhm as fwm
    from gamma.calibration.intrinsic_shape.core import poly_eval
    
    res, text = capture(sn.run_sample_subtract, bg_path, sample_path, cal_res, 
                        list(fw_points), out_zeroed, out_signed, 
                        name=name, date=date, lt_out=lt_out, log=log)
    
    resid = np.asarray(res["resid"], float)
    var = np.asarray(res["var"], float)
    n = len(resid)
    tot = float(resid.sum())
    sig = float(np.sqrt(var.sum()))
    z = tot / sig if sig > 0 else 0.0
    
    files = [
        file_check(out_zeroed, res["zeroed"], False),
        file_check(out_signed, res["signed"], True)
    ]
    
    fw_ok, fw_vals, fw_vertex, fw_margin = fwm.check(res["fwhm"], n)
    
    files_ok = all(f["match"] and f["all_int"] and f["has_fwhm"] and f["min_ok"] for f in files)
    
    S = resid + res["k"] * np.asarray(res["bg_on_sample"], float)
    
    return {
        "k": float(res["k"]),
        "sum": tot,
        "sigma": sig,
        "z": float(z),
        "z_ok": True,
        "z_limit": 2.0,
        "z_info_only": True,
        "files": files,
        "fwhm": [float(v) for v in res["fwhm"]],
        "fwhm_ok": bool(fw_ok),
        "fwhm_values": fw_vals,
        "fwhm_margin": float(fw_margin),
        "struct": [],
        "struct_ok": True,
        "n_struct_computed": 0,
        "struct_summary": {"struct_ok": True, "n": 0, "n_eval": 0, "n_ok": 0, "n_not_eval": 0},
        "converged": True,
        "lt_out": res["lt_out"],
        "live_s": float(res["live_s"]),
        "live_b": float(res["live_b"]),
        "live_out": float(res["live_out"]),
        "excluded": [],
        "all_ok": bool(fw_ok and files_ok),
        "channels": np.arange(n),
        "energy": poly_eval(list(res["coeffs_xml"]), np.arange(n, dtype=float)),
        "sample": S,
        "bg_k": res["k"] * np.asarray(res["bg_on_sample"], float),
        "resid": resid,
        "zeroed": np.asarray(res["zeroed"]),
        "signed": np.asarray(res["signed"]),
        "n_negative": int((resid < 0).sum()),
        "coeffs_xml": list(res["coeffs_xml"]),
        "sum_window": [0, n - 1],
        "log": text,
        "sample_grid": True
    }


def decision_html(r):
    """HTML итога: рекомендация, обе шкалы на реперах и независимых линиях, перевод фона, предупреждения."""
    import html
    d = r["decision"]
    lines = ["<h3>Шкала XML</h3>"]
    fmt = lambda v: "%.2f" % v if v is not None and np.isfinite(v) else "--"
    
    if d["replaced"]:
        lines.append("<p><b>Записана шкала по реперам</b>: степень %d, LOO %s кэВ по всем реперам.</p>" % (d["deg"], fmt(d["loo_rms"])))
    elif d.get("recommend"):
        lines.append("<p><b>Рекомендована шкала по реперам</b> (степень %d, LOO %s кэВ): СКО калибровки файла на реперах %s кэВ -- больше LOO с запасом (x1.5 и +1 кэВ). В XML по умолчанию -- калибровка файла пробы; применить -- кнопкой «Записать шкалу по реперам».</p>" % (d["deg"], fmt(d["loo_rms"]), fmt(d.get("rms_file_repers"))))
    else:
        lines.append("<p><b>Калибровка файла пробы</b> -- шкала по реперам не рекомендована: %s</p>" % html.escape(d.get("why", "")))
        if d.get("need"):   # правка Claude (VERIFY8): что нужно для замены и что есть сейчас
            lines.append("<p><i>%s</i></p>" % html.escape(d["need"]))
    
    lines.append('<table border="1" cellspacing="0" cellpadding="3"><tr><th></th><th>Калибровка файла</th><th>Шкала по реперам</th></tr>')
    lines.append("<tr><td>СКО на реперах, кэВ</td><td>%s</td><td>%s (LOO)</td></tr>" % (fmt(d.get("rms_file_repers")), fmt(d.get("loo_rms"))))
    lines.append("<tr><td>СКО на %d независимых линиях, кэВ</td><td>%s</td><td>%s</td></tr>" % (len(d.get("independent", [])), fmt(d.get("rms_file")), fmt(d.get("rms_new"))))
    lines.append("<tr><td>Степень</td><td>%d</td><td>%d</td></tr></table>" % (len(d.get("coeffs_file", [])) - 1, d.get("deg", 0)))
    
    if r["mode"] == "visible":
        lines.append("<p>Фон переведён на сетку пробы: по видимым собственным структурам LaBr3</p>")
    else:
        lines.append("<p>Фон переведён на сетку пробы: через энергии (%s)</p>" % html.escape(r["transfer"].get("bg_scale_how", "")))
        
    for w in r["warnings"]:
        lines.append('<p style="background:#fff3cd; color:#664d03; padding:2px">⚠ %s</p>' % html.escape(w))
        
    return "\n".join(lines)
