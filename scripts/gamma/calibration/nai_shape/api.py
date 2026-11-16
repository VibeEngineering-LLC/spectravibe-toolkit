"""
Единственная публичная точка входа метода NaI «форма + вычет фона».
СТАТУС МЕТОДА — ЭКСПЕРИМЕНТАЛЬНЫЙ: каждый результат и каждый выход несёт пометку
«экспериментально, шкала образца не проверяется». Флаг «шкала проверена» ставится
только при всех условиях (режим «форма», f ≥ порога, зелёный вердикт по ≥ min_windows
окнам, живое время подтверждено, для пар фон-фон — подтверждено, что пара фон-фон).
Правила — `audit/nai_runs/criteria_v6.md`.
Ревизия 7: печать HMAC убрана; `write_outputs` пересчитывает флаг из содержимого результата и пишет только согласованный.
"""

import json
import math
import os
import pathlib
import re
import warnings
from xml.sax.saxutils import escape

import numpy as np

from gamma.calibration.intrinsic_shape.core import poly_eval, poly_deriv
from . import pipeline, peakreg, checks, bgcal, strongcheck

# ---------------------------------------------------------------------------
# Константы
# ---------------------------------------------------------------------------

__all__ = [
    "run", "write_outputs", "live_check", "live_decision", "header_plausible",
    "choose_mode", "scale_status", "apply_status", "background_verdict",
    "bgpair_check", "natural_lines", "profile_path", "load_profile",
    "KINDS", "FLAG_UNVERIFIED", "STATUS", "EXPERIMENTAL_NOTE",
]

KINDS = ("sample", "background")
FLAG_UNVERIFIED = "шкала не проверена"
STATUS = "experimental"
EXPERIMENTAL_NOTE = "экспериментально, шкала образца не проверяется"

PROFILE_DIR = pathlib.Path(__file__).resolve().parent / "profiles"
PROFILE_ALIASES = {"gamma1s": "nai_gamma1s", "gs2020": "nai_gs2020"}

# 9. (ревизия 7) Печать HMAC убрана: не защищала от умысла и мешала записи сохранённого результата.



# ---------------------------------------------------------------------------
# 1. Пути профилей
# ---------------------------------------------------------------------------

def profile_path(name: str) -> pathlib.Path:
    """Возвращает путь к JSON-профилю, учитывая алиасы."""
    if name in PROFILE_ALIASES:
        new = PROFILE_ALIASES[name]
        warnings.warn(
            "имя профиля %r устарело, используйте %r" % (name, new),
            DeprecationWarning,
            stacklevel=2,
        )
        name = new
    return PROFILE_DIR / (name + ".json")


def load_profile(name: str) -> dict:
    """Загружает JSON-профиль по имени."""
    return pipeline.load_profile(profile_path(name))


# ---------------------------------------------------------------------------
# 2. Функция ФПШ
# ---------------------------------------------------------------------------

def _fwhm_kev(profile: dict, E: float) -> float:
    """Оценка ФПШ (кэВ) по энергии E на основе относительного значения при 662 кэВ."""
    return profile["fwhm_rel662"] * 661.7 * math.sqrt(max(E, 1.0) / 661.7)


# ---------------------------------------------------------------------------
# 3. Проверка правдоподобности заголовка
# ---------------------------------------------------------------------------

def header_plausible(X: dict, dead_max: float) -> bool:
    """Проверяет, что live ≤ real и мёртвая доля не превышает dead_max."""
    return bool(
        X["real"] > 0
        and X["live"] <= X["real"] * (1.0 + 1e-6)
        and (1.0 - X["live"] / X["real"]) <= dead_max
    )


# ---------------------------------------------------------------------------
# 4. Чистое решение по живому времени
# ---------------------------------------------------------------------------

def live_decision(
    ratio, sigma, nS, nB, header_ok: bool, profile: dict
) -> tuple[bool, str]:
    """Чистая функция решения: возвращает (ok, reason) на первом провале."""
    min_counts = profile.get("live_check_min_counts", 400)
    sigma_max = profile.get("live_check_sigma_max", 0.005)
    tol = profile.get("live_check_tol", 0.0)
    nsig = profile.get("live_check_nsigma", 3.0)

    # 1. Правдоподобность заголовка
    if not header_ok:
        return (False, "заголовок живого времени неправдоподобен (live > real или мёртвая доля велика)")
    # 2. Достаточность отсчётов
    if nS < min_counts or nB < min_counts:
        return (False, "мало отсчётов в области проверки живого времени")
    # 3. Достаточность точности
    if ratio is None or sigma is None or sigma > sigma_max:
        return (False, "точность проверки живого времени недостаточна")
    # 4. Согласование ratio с 1.0
    if abs(ratio - 1.0) > nsig * sigma + tol:
        return (False, "живое время заголовка не согласуется с данными")
    return (True, "")


# ---------------------------------------------------------------------------
# 5. Проверка живого времени
# ---------------------------------------------------------------------------

def live_check(bg: dict, sm: dict, profile: dict, nuclides, E_sm=None, E_bg=None) -> dict:
    """
    Сравнивает живое время образца и фона через отношение сумм в области выше линий.
    nuclides: список имён или None (все нуклиды профиля), пустой список — линий нет.
    """
    lo, hi = profile["fit_range_kev"]

    # Собираем линии
    if nuclides is None or len(nuclides) > 0:
        lines = pipeline.sample_lines(profile, nuclides)
    else:
        lines = []

    # Граница области: выше максимальной линии + 3*ФПШ
    if lines:
        e_start = max(lines) + 3.0 * _fwhm_kev(profile, max(lines))
    else:
        e_start = lo
    e_start = max(e_start, lo)

    # Все попарные суммы (включая себя-с-собой)
    sums = [a + b for i, a in enumerate(lines) for b in lines[i:]]

    # Суммы в области для образца и фона
    def _region_sum(X, E_in=None):
        """Число отсчётов в области [e_start, hi] без окон ±2 ПШПВ у сумм линий; поправка 6.1: интегрирование по
границам каналов (доли каналов на краях), E_in — энергии центров каналов на общей шкале."""
        C = np.asarray(X["counts"], dtype=float).copy()
        C[-1] = 0.0
        n = len(C)

        E = np.asarray(E_in, dtype=float) if E_in is not None else poly_eval(X["cal"], np.arange(n, dtype=float))

        edges = np.empty(n + 1)
        edges[1:-1] = 0.5 * (E[:-1] + E[1:])
        edges[0] = E[0] - 0.5 * (E[1] - E[0])
        edges[-1] = E[-1] + 0.5 * (E[-1] - E[-2])

        cum = np.concatenate([[0.0], np.cumsum(C)])
        N = lambda a, b: float(np.interp(b, edges, cum) - np.interp(a, edges, cum)) if b > a else 0.0

        iv = sorted(
            (max(Es - 2.0 * _fwhm_kev(profile, Es), e_start),
             min(Es + 2.0 * _fwhm_kev(profile, Es), hi))
            for Es in sums
        )
        merged = []
        for a, b in iv:
            if a >= b:
                continue
            if merged and a <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))

        total = N(e_start, hi) - sum(N(a, b) for a, b in merged)
        return float(total)
    nS = _region_sum(sm, E_sm)
    nB = _region_sum(bg, E_bg)

    # Показатели из заголовков
    k_header = sm["live"] / bg["live"]
    f_header = peakreg.bg_fraction(bg, sm, profile)

    # Проверка заголовков
    dead_max = profile.get("dead_fraction_max", 0.002)
    header_ok = header_plausible(sm, dead_max) and header_plausible(bg, dead_max)

    # Данные
    if nS > 0 and nB > 0:
        k_data = nS / nB
        ratio = k_data / k_header
        sigma = ratio * math.sqrt(1.0 / nS + 1.0 / nB)
        k_eff = min(k_header, k_data)
    else:
        k_data = None
        ratio = None
        sigma = None
        k_eff = k_header

    ok, reason = live_decision(ratio, sigma, nS, nB, header_ok, profile)
    f = f_header * k_eff / k_header

    return {
        "k_header": float(k_header),
        "k_data": float(k_data) if k_data is not None else None,
        "ratio": float(ratio) if ratio is not None else None,
        "sigma": float(sigma) if sigma is not None else None,
        "n_region": int(min(nS, nB)),
        "e_start": float(e_start),
        "header_ok": bool(header_ok),
        "ok": bool(ok),
        "reason": reason,
        "f_header": float(f_header),
        "f": float(f),
    }


# ---------------------------------------------------------------------------
# 6. Выбор режима
# ---------------------------------------------------------------------------

def choose_mode(f: float, thr: float, force=None) -> str:
    """Выбирает режим: 'shape' если f >= thr, иначе 'peaks'. force переопределяет."""
    if force is not None:
        if force not in ("shape", "peaks"):
            raise ValueError("force должен быть 'shape', 'peaks' или None, получено: %r" % force)
        return force
    return "shape" if f >= thr else "peaks"


# ---------------------------------------------------------------------------
# 7. Статус шкалы
# ---------------------------------------------------------------------------

def scale_status(
    mode_auto: str,
    f: float,
    thr: float,
    verdict_ok,
    live_ok,
    n_windows=None,
    min_windows: int = 3,
    bgpair_ok: bool = True,
) -> tuple[bool, list[str]]:
    """
    Возвращает (verified, reasons). reasons — по строке на каждое невыполненное условие.
    """
    reasons: list[str] = []

    if mode_auto != "shape":
        reasons.append("режим по линиям (peaks/peaks1)")
    if not (f >= thr):
        reasons.append("доля фона f=%.3f меньше порога %.3f" % (f, thr))
    if n_windows is not None and n_windows < min_windows:
        reasons.append("нечем проверить: окон %d < %d" % (n_windows, min_windows))
    if verdict_ok is not True:
        reasons.append("вердикт окон не зелёный")
    if live_ok is not True:
        reasons.append("живое время не подтверждено данными")
    if bgpair_ok is not True:
        reasons.append("пара не подтверждена как фон-фон")

    verified = len(reasons) == 0
    return (verified, reasons)


# ---------------------------------------------------------------------------
# 8. Применение статуса к результату
# ---------------------------------------------------------------------------

def apply_status(
    r: dict,
    f: float,
    thr: float,
    live_ok: bool,
    min_windows: int = 3,
    bgpair_ok: bool = True,
) -> dict:
    """Применяет статус шкалы, запечатывает результат HMAC и возвращает r."""
    ok_raw = r["verdict"].get("ok")
    n_w = r["verdict"].get("n")

    verified, reasons = scale_status(
        r["mode_auto"], f, thr, ok_raw, live_ok, n_w, min_windows, bgpair_ok
    )

    r["scale_verified"] = bool(verified)
    r["scale_reasons"] = reasons
    r["flag"] = "" if verified else FLAG_UNVERIFIED
    r["status"] = STATUS
    r["status_note"] = EXPERIMENTAL_NOTE

    if verified:
        r["verdict"] = dict(r["verdict"], ok_raw=ok_raw, counted=True)
    else:
        r["verdict"] = dict(
            r["verdict"],
            ok_raw=ok_raw,
            ok=None,
            counted=False,
            note=FLAG_UNVERIFIED + ", в зачёт не идёт: " + "; ".join(reasons),
        )

    r["thresholds"] = {"thr": float(thr), "min_windows": int(min_windows)}   # ревизия 7: для пересчёта флага при записи
    return r


# ---------------------------------------------------------------------------
# 9. (ревизия 7) Печать HMAC убрана: не защищала от умысла и мешала записи сохранённого результата.
# ---------------------------------------------------------------------------
# 10. Вердикт по окнам (фон)
# ---------------------------------------------------------------------------

def background_verdict(r: dict, profile: dict) -> dict:
    """
    Проверяет отсутствие пиков в окнах нулевого фона (правило z_dip с боковыми полосами).
    """
    net = r["arrays"]["net"]
    var = r["arrays"]["var"]
    cal = {"ecal": r["ecal_bg"], "fwc": r["fwc"], "R662": r["R662"], "n": len(net)}

    fw = bgcal.fwhm_ch_fn(cal)
    e2ch = lambda E: float(bgcal.energy_to_ch(cal, E))
    fwhm_kev = lambda E: float(fw(e2ch(E)) * poly_deriv(cal["ecal"], e2ch(E)))

    wl = [w["E"] for w in profile["null_windows"]]
    tags = {w["E"]: w["tag"] for w in profile["null_windows"]}

    T = checks.window_tests(net, var, e2ch, fwhm_kev, wl, sideband=True, tags=tags)

    z_thr = profile.get("z_thr", 3.0)
    rms_max = profile.get("null_rms_max", 1.5)
    min_w = profile.get("min_windows", 3)

    zd = [t["z_dip"] for t in T if not t["excluded"]]
    n = len(zd)
    rms = math.sqrt(sum(z * z for z in zd) / n) if n else 0.0
    bad = [t["E"] for t in T if (not t["excluded"]) and abs(t["z_dip"]) > z_thr]
    ok = bool(n >= min_w and len(bad) == 0 and rms <= rms_max)

    max_abs_z_dip = max(abs(z) for z in zd) if zd else 0.0
    z_sums = [t["z_sum"] for t in T if not t["excluded"]]
    max_abs_z_sum = max(abs(z) for z in z_sums) if z_sums else 0.0

    return {
        "rule": "z_dip по окнам с боковыми полосами",
        "n": n,
        "n_ok": n - len(bad),
        "max_abs_z_dip": max_abs_z_dip,
        "rms_z_dip": rms,
        "max_abs_z_sum": max_abs_z_sum,
        "bad": bad,
        "ok": ok,
        "tests": T,
    }


# ---------------------------------------------------------------------------
# 11. Естественные линии фона
# ---------------------------------------------------------------------------

def natural_lines(profile: dict) -> list[float]:
    """Список энергий (кэВ) естественных линий, которые не должны давать пик в разности."""
    imin = profile.get("bgpair_natural_imin", 5.0)   # поправка 6.1 criteria_v6.md (было 2.0)
    E = [
        e for e, _i in strongcheck.sample_lines(["K-40", "Ra-226", "Th-232"], imin=imin)
    ] + [511.0] + [w["E"] for w in profile["null_windows"]]
    return sorted(set(float(e) for e in E))


# ---------------------------------------------------------------------------
# 12. Проверка пары фон-фон
# ---------------------------------------------------------------------------

def bgpair_check(r: dict, profile: dict) -> dict:
    """
    Правило П6.2: скан |z_sum| по сетке вне естественных линий.
    """
    net = r["arrays"]["net"]
    var = r["arrays"]["var"]
    cal = {"ecal": r["ecal_bg"], "fwc": r["fwc"], "R662": r["R662"], "n": len(net)}

    fw = bgcal.fwhm_ch_fn(cal)
    e2ch = lambda E: float(bgcal.energy_to_ch(cal, E))
    fwhm_kev = lambda E: float(fw(e2ch(E)) * poly_deriv(cal["ecal"], e2ch(E)))

    lo, hi = profile["fit_range_kev"]
    z_max = profile.get("bgpair_scan_z", 5.0)
    n_min = profile.get("bgpair_min_windows", 10)
    excl = profile.get("bgpair_excl_fwhm", 1.0)
    nat = natural_lines(profile)

    # Построение сетки
    grid: list[float] = []
    E = lo + 2.0 * fwhm_kev(lo)
    hi_limit = hi - 2.0 * fwhm_kev(hi)
    while E <= hi_limit:
        if all(abs(E - X) >= excl * fwhm_kev(X) for X in nat):
            grid.append(E)
        E += 0.5 * fwhm_kev(E)

    # Тесты окон
    if grid:
        T = checks.window_tests(net, var, e2ch, fwhm_kev, grid, sideband=True)
    else:
        T = []

    zs = [abs(t["z_sum"]) for t in T]
    n = len(zs)
    zmax = max(zs) if zs else 0.0
    E_worst = T[int(np.argmax(zs))]["E"] if zs else None

    # Решение
    if n < n_min:
        ok = False
        reason = "нечем проверить пару фон-фон: окон %d < %d" % (n, n_min)
    elif zmax > z_max:
        ok = False
        reason = "в разности линия вне окон фона: E=%.1f кэВ, |z|=%.1f" % (E_worst, zmax)
    else:
        ok = True
        reason = ""

    return {
        "rule": "скан |z_sum| вне линий естественного фона",
        "n": n,
        "max_abs_z": zmax,
        "E_worst": E_worst,
        "z_thr": z_max,
        "ok": ok,
        "reason": reason,
    }


# ---------------------------------------------------------------------------
# 13. Главный вход
# ---------------------------------------------------------------------------

def run(
    bg: dict,
    sm: dict,
    profile: dict,
    nuclides=None,
    kind: str = "sample",
    force=None,
    verbose: bool = False,
) -> dict:
    """
    Единственная публичная точка входа метода NaI «форма + вычет фона».
    kind: 'sample' или 'background'.
    nuclides: список имён нуклидов, None — все из профиля, [] — не заданы.
    force: None, 'shape' или 'peaks' (только для sample).
    """
    # --- Валидация входных параметров ---
    if kind not in KINDS:
        raise ValueError("kind должен быть из %r, получено: %r" % (KINDS, kind))
    if kind == "background" and force is not None:
        raise ValueError("для kind='background' параметр force не допускается")
    if kind == "sample" and nuclides is not None and len(nuclides) == 0:
        raise ValueError(
            "kind='sample' требует хотя бы один нуклид образца (или None — все нуклиды профиля)"
        )

    # --- Живое время и доля фона ---
    lc = live_check(bg, sm, profile, nuclides if kind == "sample" else [])
    f = lc["f"]
    thr = profile["shape_min_bg_fraction"]
    min_w = profile.get("min_windows", 3)
    bgpair_ok = True

    # --- Выполнение пайплайна ---
    if kind == "background":
        r = pipeline._run_pair(bg, sm, profile, mode="background", verbose=verbose)
        r["mode_auto"] = "shape"
        r["verdict_fit"] = r["verdict"]
        r["verdict"] = background_verdict(r, profile)
        r["bgpair"] = bgpair_check(r, profile)
        bgpair_ok = r["bgpair"]["ok"]
    else:
        mode = choose_mode(f, thr, force)
        if mode == "shape":
            r = pipeline._run_pair(
                bg, sm, profile, mode="sample", nuclides=nuclides, verbose=verbose
            )
        else:
            r = peakreg._register_by_lines(bg, sm, profile, nuclides, verbose=verbose)
        r["mode_auto"] = mode

    # поправка 6.1: живое время проверяется на общей шкале (шкала фона + отображение h образца), не на шкалах файлов
    E_bg = poly_eval(r["ecal_bg"], np.arange(len(bg["counts"]), dtype=float))
    E_sm = poly_eval(r["ecal_bg"], poly_eval(r["c"], np.arange(len(sm["counts"]), dtype=float)))
    lc_map = live_check(bg, sm, profile, nuclides if kind == "sample" else [], E_sm, E_bg)
    lc = dict(lc_map, f=lc["f"], f_header=lc["f_header"], ok_prefit=lc["ok"])

    # --- Метаданные ---
    r["kind"] = kind
    r["bg_fraction"] = f
    r["bg_fraction_header"] = lc["f_header"]
    r["live_check"] = lc
    r["sample_file_cal"] = [float(c) for c in sm["cal"]]
    r["profile_device"] = profile.get("device", "")

    # --- Статус и печать ---
    return apply_status(r, f, thr, lc["ok"], min_w, bgpair_ok)


# ---------------------------------------------------------------------------
# 14. Проверка согласованности результата (ревизия 7: печать HMAC убрана)
# ---------------------------------------------------------------------------

def _check_consistent(res):
    """Пересчитывает флаг шкалы из содержимого
    результата и сверяет с записанным; бросает ValueError при расхождении или отсутствии полей api.run. Защищает от
    небрежной правки словаря, не от умысла."""
    need = ("scale_verified", "scale_reasons", "flag", "verdict", "mode_auto", "bg_fraction", "live_check", "thresholds", "status")
    missing = [k for k in need if k not in res]
    if missing:
        raise ValueError("результат без статуса api.run (нет полей: %s)" % ", ".join(missing))
    th = res["thresholds"]
    v = res["verdict"]
    verified, reasons = scale_status(res["mode_auto"], res["bg_fraction"], th["thr"], v.get("ok_raw"), (res["live_check"] or {}).get("ok"), v.get("n"), th["min_windows"], (res.get("bgpair") or {}).get("ok", True))
    bad = []
    if bool(res["scale_verified"]) is not bool(verified):
        bad.append("scale_verified не следует из содержимого")
    if list(res["scale_reasons"]) != list(reasons):
        bad.append("причины не совпадают с пересчётом")
    if res["flag"] != ("" if verified else FLAG_UNVERIFIED):
        bad.append("флаг не совпадает с пересчётом")
    if v.get("counted") is not bool(verified):
        bad.append("counted не совпадает с пересчётом")
    if verified and v.get("ok") is not True:
        bad.append("ok вердикта не True при проверенной шкале")
    if (not verified) and v.get("ok") is not None:
        bad.append("ok вердикта не None при непроверенной шкале")
    if res["status"] != STATUS:
        bad.append("статус не experimental")
    if bad:
        raise ValueError("результат несогласован: " + "; ".join(bad))
# ---------------------------------------------------------------------------
# 15. Текстовое представление статуса
# ---------------------------------------------------------------------------

def _status_text(res: dict) -> str:
    """Собирает строку статуса для вставки в выходной XML."""
    parts: list[str] = []
    parts.append(EXPERIMENTAL_NOTE + "; ")
    if res["scale_verified"]:
        parts.append("шкала проверена")
    else:
        parts.append(FLAG_UNVERIFIED + ": " + "; ".join(res["scale_reasons"]))
    parts.append(
        "; вердикт окон ok=%s (до флага %s)"
        % (res["verdict"].get("ok"), res["verdict"].get("ok_raw"))
    )
    cal_str = ", ".join("%.7g" % c for c in res.get("sample_file_cal", []))
    parts.append("; исходная шкала файла образца: " + cal_str)
    return "".join(parts)


# ---------------------------------------------------------------------------
# 15a. Пометка SampleInfo в XML (ревизия 7: CDATA, отсутствующие элементы, контроль результата)
# ---------------------------------------------------------------------------
_CDATA_OR_TEXT = r"(?:<!\[CDATA\[(.*?)\]\]>|([^<]*))"


def _mark_sampleinfo(text, name, note_add):
    """Записать имя и дописать примечание в первый блок
    SampleInfo: понимает CDATA и пустые элементы, создаёт отсутствующие Name/Note. Возвращает новый текст; нет SampleInfo —
    ValueError."""
    from xml.sax.saxutils import unescape

    m = re.search(r"<SampleInfo>(.*?)</SampleInfo>", text, flags=re.S)
    if m is None:
        raise ValueError("в XML нет блока SampleInfo")
    block = m.group(1)

    new_name = "<Name>" + escape(name) + "</Name>"
    block, k = re.subn(r"<Name>" + _CDATA_OR_TEXT + r"</Name>|<Name\s*/>", lambda mm: new_name, block, count=1, flags=re.S)
    if k == 0:
        block = new_name + block

    def _note(mm):
        old = mm.group(1) if mm.group(1) is not None else unescape(mm.group(2) or "")
        txt = (old + "; " if old else "") + note_add
        return "<Note>" + escape(txt) + "</Note>"

    block, k = re.subn(r"<Note>" + _CDATA_OR_TEXT + r"</Note>|<Note\s*/>", _note, block, count=1, flags=re.S)
    if k == 0:
        block = block + "<Note>" + escape(note_add) + "</Note>"

    return text[:m.start(1)] + block + text[m.end(1):]


def _verify_marks(path, need_name, need_note):
    """Прочитать записанный XML читателем пакета и убедиться, что
    пометки есть в имени и примечании; иначе RuntimeError."""
    from gamma.io.atomspectra_xml import read_atomspectra_xml

    s = read_atomspectra_xml(str(path))
    miss = [t for t in need_name if t not in (s.sample_id or "")] + \
           [t for t in need_note if t not in (s.comments or "")]
    if miss:
        raise RuntimeError("пометка не записана в %s: %s" % (path, miss))
# ---------------------------------------------------------------------------
# 16. Запись выходных файлов
# ---------------------------------------------------------------------------

def write_outputs(
    res: dict,
    sm: dict,
    out_dir,
    stem: str,
    date: str = "",
) -> dict:
    """
    Записывает CSV, XML (нулевой и подписанный), JSON.
    Перед записью проверяет печать и согласованность результата.
    """
    _check_consistent(res)
    paths = pipeline._write_outputs(res, sm, out_dir, stem, date)
    # --- Дополнение CSV ---
    csv_path = paths["csv"]
    with open(csv_path, encoding="utf-8") as fh:
        lines = fh.read().rstrip("\n").split("\n")

    csv_header_extra = ",status,scale_verified,flag,verdict_ok,verdict_ok_raw"
    csv_row_extra = (
        ",%s,%d,%s,%s,%s"
        % (
            STATUS,
            int(res["scale_verified"]),
            res["flag"],
            res["verdict"].get("ok"),
            res["verdict"].get("ok_raw"),
        )
    )

    new_lines: list[str] = []
    for idx, line in enumerate(lines):
        if not line.strip():
            new_lines.append(line)
        elif idx == 0:
            new_lines.append(line + csv_header_extra)
        else:
            new_lines.append(line + csv_row_extra)

    with open(csv_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(new_lines) + "\n")

    # --- Дополнение XML (ревизия 7: через _mark_sampleinfo, с проверкой записанного) ---
    tag = " [экспериментально]" + ((" [" + res["flag"] + "]") if res["flag"] else "")
    need_name = ["[экспериментально]"] + (["[" + res["flag"] + "]"] if res["flag"] else [])
    status_text = _status_text(res)
    for key, suffix in (("xml_zeroed", " (обнулённый)"), ("xml_signed", " (со знаком)")):
        p = paths[key]
        with open(p, encoding="utf-8", newline="") as f:
            content = f.read()
        content = _mark_sampleinfo(content, stem + suffix + tag, status_text)
        with open(p, "w", encoding="utf-8", newline="") as f:
            f.write(content)
        _verify_marks(p, need_name, [EXPERIMENTAL_NOTE])
    return paths

