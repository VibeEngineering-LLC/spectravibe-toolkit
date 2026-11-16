"""Обработчики слоя детекторов и метода NaI без интерфейса: список профилей, автоопределение прибора, запуск nai_shape.api в фоновом потоке, справка. Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-11-16."""

import html
import pathlib
import subprocess

from . import backend

LA_METHOD = "la138_intrinsic"
NAI_METHOD = "nai"
METHOD_TITLES = {"la138_intrinsic": "по форме La-138", "nai": "по линиям ЕРН (nai_shape)"}
STATUS_TITLES = {"": "рабочий", "unverified": "не проверен (UNVERIFIED)", "experimental": "экспериментально"}
METHODOLOGY = "references/intrinsic_shape_calibration.md"
BASE_VERSION = "0.2.0"
KIND_CHOICES = [("sample", "Образец минус фон"), ("background", "Проверка пары фон-фон")]
FORCE_CHOICES = [(None, "Авто (по доле фона: форма или по линиям)"), ("shape", "Форма"), ("peaks", "По линиям (peaks; при одной линии — peaks1)")]


def profile_rows() -> list[dict]:
    from gamma.calibration.intrinsic_shape import profiles

    rows = []
    for p in profiles.list_profiles(all_methods=True):
        pid = p["id"]
        title = p.get("title", pid)
        method = p.get("method", "")
        method_title = METHOD_TITLES.get(method, method)
        status = p.get("status", "")
        status_title = STATUS_TITLES.get(status, status)
        limitations = list(p.get("limitations", []))
        aliases = list(p.get("aliases", []))

        label = title
        if "[экспериментально]" not in label and status == "experimental":
            label += " [экспериментально]"
        label += " — метод " + method_title + " — " + status_title + " (" + pid + ")"

        rows.append({
            "id": pid,
            "title": title,
            "method": method,
            "method_title": method_title,
            "status": status,
            "status_title": status_title,
            "experimental": status == "experimental",
            "status_note": p.get("status_note", ""),
            "limitations": limitations,
            "aliases": aliases,
            "has_match": bool(p.get("match")),
            "label": label,
        })
    return rows


def method_of(pid) -> str:
    from gamma.calibration.intrinsic_shape import profiles
    return profiles.get(pid).get("method", "")


def device_name(path) -> str:
    from gamma.io.readers import read_spectrum
    try:
        obj = read_spectrum(str(path))
        return (getattr(obj, "detector_id", "") or "").strip()
    except Exception:
        return ""


def detect_for(path) -> dict:
    from gamma.calibration.intrinsic_shape import profiles
    dev = device_name(path)
    if dev:
        pid, msg = profiles.detect(dev)
    else:
        pid, msg = None, ""
    return {"device": dev, "pid": pid, "message": msg}


def mismatch_text(selected_pid, det) -> str:
    dev = det["device"]
    if dev == "":
        return ""
    if det["pid"] is not None and det["pid"] != selected_pid:
        return "Прибор в файле «%s» распознан как профиль «%s», а выбран «%s» — проверьте выбор детектора" % (dev, det["pid"], selected_pid)
    if det["pid"] is None and method_of(selected_pid) == LA_METHOD:
        return "Прибор в файле «%s» не распознан ни одним профилем La-138; выбран «%s» — проверьте выбор детектора" % (dev, selected_pid)
    return ""


def hint_text(det) -> str:
    if det["pid"] is None:
        return "Прибор в файле не распознан автоматически (у профилей NaI шаблона имени прибора нет) — выберите профиль детектора вручную"
    return "Прибор в файле: «%s» → профиль «%s»" % (det["device"], det["pid"])


def passport_html(pid) -> str:
    rows = profile_rows()
    row = None
    for r in rows:
        if r["id"] == pid:
            row = r
            break
    if row is None:
        return ""
    e = html.escape
    title = e(row["title"])
    method_title = e(row["method_title"])
    status_title = e(row["status_title"])
    status = row["status"]
    status_note = e(row["status_note"])
    limitations = row["limitations"]

    parts = ["<b>Паспорт профиля:</b> " + title + "; метод: " + method_title + "; статус: "]
    if status != "":
        parts.append('<span style="color:#b00000"><b>' + status_title + "</b></span>")
    else:
        parts.append(status_title)
    if row["status_note"]:
        parts.append("<br>" + status_note)
    for lim in limitations:
        parts.append("<br>• " + e(lim))
    return "".join(parts)


def nuclide_names(pid) -> list[str]:
    from gamma.calibration.nai_shape import api as nai_api
    return sorted(nai_api.load_profile(pid).get("sample_lines", {}).keys())


def _check_cancel(where):
    from gamma.calibration.intrinsic_shape import robust
    if robust.CANCEL["flag"]:
        raise robust.Cancelled("остановлено: " + where)


def run_nai(bg_path, sm_path, pid, kind="sample", force=None, nuclides=None, out_dir="", stem="", log=None) -> dict:
    from gamma.calibration.intrinsic_shape import profiles
    from gamma.calibration.intrinsic_shape.core import load_spectrum, poly_eval
    from gamma.calibration.nai_shape import api as nai_api
    import numpy as np

    if method_of(pid) != NAI_METHOD:
        raise backend.UserError("Профиль «%s» не относится к методу NaI" % pid)

    bgp = backend.check_input(bg_path)
    smp = backend.check_input(sm_path)

    if not str(out_dir).strip():
        raise backend.UserError("Укажите папку для результата")
    outp = pathlib.Path(out_dir)
    if not outp.is_dir():
        raise backend.UserError("Папка не существует: %s" % outp)

    stem_val = (stem or "").strip() or smp.stem + "_минус_фон_NaI"
    if stem_val.lower().endswith(".xml"):
        stem_val = stem_val[:-4]

    def _load(p):
        name = p.name   # в сообщении -- имя файла, не полный путь
        try:
            data = load_spectrum(p)
        except Exception as e:
            raise backend.UserError("Не удалось прочитать спектр %s: %s" % (name, e))
        if data["counts"].sum() <= 0:
            raise backend.UserError("В файле %s нет отсчётов спектра" % name)
        if data["live"] <= 0:
            raise backend.UserError("В файле %s нет живого времени (live time)" % name)
        return data

    bg = _load(bgp)
    sm = _load(smp)

    if kind == "background":
        force = None
    profile = nai_api.load_profile(pid)
    if kind == "background" or not nuclides:
        nuclides = None
    if kind == "sample":   # режим по линиям (peaks) берёт якорные линии нуклидов: None и нуклид без якорей там падают
        f = nai_api.live_check(bg, sm, profile, nuclides)["f"]
        if nai_api.choose_mode(f, profile["shape_min_bg_fraction"], force) == "peaks":
            anchors = profile.get("anchor_lines", {})
            miss = [n for n in (nuclides or []) if n not in anchors]
            if miss:
                raise backend.UserError("Режим по линиям: у нуклидов %s нет якорных линий в профиле" % ", ".join(miss))
            nuclides = nuclides or sorted(n for n in nuclide_names(pid) if n in anchors)
    profiles.activate(pid)

    if log:
        log("[NAI] ЭКСПЕРИМЕНТАЛЬНО: " + nai_api.EXPERIMENTAL_NOTE + "\n")

    _check_cancel("до расчёта")

    try:
        res, _txt = backend.capture(nai_api.run, bg, sm, profile, nuclides, kind=kind, force=force, log=log)
    except ValueError as exc:
        raise backend.UserError("Метод NaI отказал: %s" % exc)

    _check_cancel("после расчёта, файлы не записаны")

    paths = nai_api.write_outputs(res, sm, str(outp), stem_val)

    a = res.get("arrays", {})
    Eb = a.get("Eb")
    net = a.get("net")
    fit = a.get("fit_lin")
    n_s = len(sm["counts"])
    E_sm = poly_eval(res["ecal_bg"], poly_eval(res["c"], np.arange(n_s, dtype=float)))
    plot = {
        "E_bg": Eb,
        "bg_k": np.asarray(bg["counts"], float) * float(res.get("k", 1.0)),
        "E_sm": E_sm,
        "sample": np.asarray(sm["counts"], float),
        "net": net,
        "resid": (np.asarray(net, float) - np.asarray(fit, float)) if (net is not None and fit is not None) else None,
    }

    return {
        "pid": pid,
        "profile_title": profiles.get(pid).get("title", pid),
        "kind": kind,
        "res": res,
        "paths": paths,
        "out_dir": str(outp),
        "plot": plot,
    }


def _fmt(v, fmt) -> str:
    try:
        return fmt % v
    except (TypeError, ValueError):
        return str(v)


def summary_html(r) -> str:
    e = html.escape
    res = r["res"]
    lines = []
    lines.append('<span style="color:#b00000"><b>ЭКСПЕРИМЕНТАЛЬНО: шкала не проверяется</b></span>')
    lines.append("Профиль: <b>%s</b>; режим: %s (%s)" % (e(r["profile_title"]), e(str(res.get("mode_auto", ""))), e(str(res.get("mode", "")))))
    verdict = res.get("verdict", {})
    lines.append("Вердикт окон: ok=%s, окон %s; %s" % (e(str(verdict.get("ok", ""))), e(str(verdict.get("n", ""))), e(str(verdict.get("note", "")))))
    c = res.get("c", [0, 1, 0])
    c0 = c[0] if len(c) > 0 else 0.0
    c1 = c[1] if len(c) > 1 else 1.0
    c2 = c[2] if len(c) > 2 else 0.0
    lines.append("Отображение шкалы образца h(x) = c0 + c1·x + c2·x²: смещение c0 = " + _fmt(c0, "%.4f") + " кан, масштаб c1 = " + _fmt(c1, "%.6f") + ", c2 = " + _fmt(c2, "%.3e"))
    lines.append("Масштаб фона k = " + _fmt(res.get("k", float("nan")), "%.6g") + "; доля фона f = " + _fmt(res.get("bg_fraction", float("nan")), "%.3f"))
    lines.append("χ²/ndf = " + _fmt(res.get("chi2ndf", float("nan")), "%.3f"))
    lines.append("Шкала проверена: <b>%s</b>" % ("да" if res.get("scale_verified") else "нет"))
    if res.get("scale_reasons"):
        lines.append("Причины: " + "; ".join(e(str(x)) for x in res["scale_reasons"]))
    lines.append("Флаг: %s" % e(str(res.get("flag") or "нет")))
    lines.append("Файлы результата:")
    for key in ("xml_zeroed", "xml_signed", "csv", "json"):
        if key in r["paths"]:
            lines.append("• " + e(pathlib.Path(r["paths"][key]).name))
    return "<br>".join(lines)


def build_sha() -> str:
    f = pathlib.Path(__file__).with_name("_build_sha.txt")
    if f.exists():
        return f.read_text(encoding="utf-8").strip()
    try:
        r = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=pathlib.Path(__file__).parent,
            capture_output=True, text=True, timeout=5,
        )
        if r.returncode == 0 and r.stdout.strip():
            return r.stdout.strip()
    except Exception:
        pass
    return "dev"


def version() -> str:
    return BASE_VERSION + "-beta." + build_sha()


def about_html() -> str:
    e = html.escape
    rows = profile_rows()
    parts = []
    parts.append("<h3>IntrinsicShapeGUI %s</h3>" % e(version()))
    parts.append("<p>Калибровка гамма-спектров по форме собственных структур кристалла (LaBr3, La-138) и экспериментальный метод NaI по линиям ЕРН. Тестовая версия: для проверки метода, не для отчётных измерений.</p>")
    parts.append("<p><b>Детекторы:</b></p><ul>")
    for row in rows:
        parts.append("<li>%s — метод %s — статус: %s</li>" % (e(row["title"]), e(row["method_title"]), e(row["status_title"])))
    parts.append("</ul>")
    parts.append("<p><b>Ограничения:</b></p><ul>")
    for row in rows:
        for lim in row["limitations"]:
            parts.append("<li>%s: %s</li>" % (e(row["title"]), e(lim)))
        if row["status_note"]:
            parts.append("<li>%s: %s</li>" % (e(row["title"]), e(row["status_note"])))
    parts.append("</ul>")
    parts.append("<p>Методика: %s (в репозитории spectravibe-toolkit).</p>" % e(METHODOLOGY))
    return "".join(parts)
