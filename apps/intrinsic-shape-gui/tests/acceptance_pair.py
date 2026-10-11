"""Приёмка пар спектров для isgui:
пара 8192 (вложенная и независимая), грунт с линиями образца, кнопка «Новый спектр». Критерии -- из
tests/acceptance_criteria.json (заданы до прогона); точечные мутации -- переменная ISGUI_ACC_MUT.
Сгенерировано gen_code.py (qwen3.8:27b) по specs/acceptance_pair.md, 2026-10-10 (раунд 4)"""

import os, re, json, time, pathlib, traceback
import numpy as np
import pyqtgraph as pg
from isgui import backend_pair
from gamma.calibration.intrinsic_shape import pipeline as pl, verdict as vd, robust as rb, register as reg

ROOT = pathlib.Path(__file__).resolve().parents[1]
FIX = ROOT.parents[1] / "tests" / "fixtures" / "intrinsic_shape"
BG = FIX / "labr3_background.xml"
SM = FIX / "labr3_lu176_sample.xml"
# реальные спектры -- только через переменные окружения, иначе шаги F/G пропускаются
def _env_path(name):
    v = os.environ.get(name, "")
    return pathlib.Path(v) if v else None
BG8, SM8 = _env_path("ISGUI_PAIR8_BG"), _env_path("ISGUI_PAIR8_SAMPLE")      # вложенная пара 8192
BG8_IND = _env_path("ISGUI_PAIR8_BG_IND")                                     # независимый фон 8192
BG_SOIL, SM_SOIL = _env_path("ISGUI_SOIL_BG"), _env_path("ISGUI_SOIL_SAMPLE")  # грунт с линиями образца
CRIT_PATH = ROOT / "tests" / "acceptance_criteria.json"
MUT = os.environ.get("ISGUI_ACC_MUT", "")
_LOAD_INPUTS = reg.load_inputs          # исходные функции -- до любых мутаций
_REGISTER_ROBUST = rb.register_robust
_WINDOW_VERDICT = vd.window_verdict
_NESTED = pl.nested_pair
_MARKERS = backend_pair.markers_from

_MASK_SWITCH_BROKEN = False

def load_criteria():
    return json.loads(CRIT_PATH.read_text("utf-8"))

def apply_mutation(name):
    if not name:
        return lambda: None
    
    if name == "F3_k":
        def _k_scaled(t):
            return (t[0], t[1], t[2] * 1.05, t[3], t[4])
        reg.load_inputs = lambda *a, **kw: _k_scaled(_LOAD_INPUTS(*a, **kw))
        return lambda: setattr(reg, "load_inputs", _LOAD_INPUTS)
        
    elif name == "F6_onesided":
        def wrapper(*a, **kw):
            res = _WINDOW_VERDICT(*a, **kw)
            if res.get("dependent"):
                res["dependent"] = False
                if "reason" in res and res["reason"]:
                    parts = res["reason"].split("; ")
                    new_parts = [p for p in parts if not p.startswith("χ² ниже статистики")]
                    if not new_parts:
                        res["ok"] = True
                        res["status"] = "OK"
                    res["reason"] = "; ".join(new_parts)
            return res
        vd.window_verdict = wrapper
        return lambda: setattr(vd, "window_verdict", _WINDOW_VERDICT)
        
    elif name == "F7_nested":
        pl.nested_pair = lambda *a, **kw: ""
        return lambda: setattr(pl, "nested_pair", _NESTED)
        
    elif name == "F8_hide":
        vd.window_verdict = lambda *a, **kw: _WINDOW_VERDICT(*a, **dict(kw, hide_max=1.0))
        return lambda: setattr(vd, "window_verdict", _WINDOW_VERDICT)
        
    elif name == "G2_plain":
        rb.register_robust = lambda *a, **kw: _REGISTER_ROBUST(*a, **dict(kw, mode="off"))
        return lambda: setattr(rb, "register_robust", _REGISTER_ROBUST)
        
    elif name == "G3_maskon":
        global _MASK_SWITCH_BROKEN
        _MASK_SWITCH_BROKEN = True
        def _undo():   # откат флага модуля
            global _MASK_SWITCH_BROKEN
            _MASK_SWITCH_BROKEN = False
        return _undo
        
    elif name == "G5_lu":
        backend_pair.markers_from = lambda r: _MARKERS(r) + [{"name": "Lu176_202", "ch": 142.9, "E": 201.83}]
        return lambda: setattr(backend_pair, "markers_from", _MARKERS)
        
    elif name == "H3_stale":
        from isgui import panels
        orig = panels.AnchorBox.reset
        def reset_keep(self):
            t = self.manual_edit.text()
            orig(self)
            self.manual_edit.setText(t)
        panels.AnchorBox.reset = reset_keep
        return lambda: setattr(panels.AnchorBox, "reset", orig)
        
    else:
        raise ValueError("неизвестная мутация " + name)

def masking_off():
    if _MASK_SWITCH_BROKEN:
        return lambda: None
    cur = rb.register_robust
    rb.register_robust = lambda *a, **kw: cur(*a, **dict(kw, exclude=False, extra_fn=None))
    return lambda: setattr(rb, "register_robust", cur)

def datapoints(path):
    text = pathlib.Path(path).read_text("utf-8")
    return re.findall(r"<DataPoint>([^<]*)</DataPoint>", text)

def chain(bg, sm, out_dir, robust="auto"):
    r = backend_pair.run_register(str(bg), str(sm), robust=robust)
    e = backend_pair.run_ecal(str(bg), str(sm), r["h"])
    out_dir.mkdir(parents=True, exist_ok=True)
    s = backend_pair.run_subtract(
        str(bg), str(sm), r["h"], e, 
        [tuple(p) for p in e["fwhm_points"]], 
        str(out_dir / "net.xml"), 
        str(out_dir / "net_signed.xml"), 
        reg_json=r["h_json"]
    )
    return (r, e, s)

def fmt_windows(checks):
    parts = []
    for c in checks:
        parts.append("%s %s скрыто %.0f%% z_sum %+.2f z_dip %+.2f z_χ² %+.2f" % (
            c["name"], c["status"], c["hidden_frac"] * 100, c["z_sum"], c["z_dip"], c["z_chi2"]
        ))
    return "; ".join(parts)

def print_table(tag, checks):
    for c in checks:
        print("[%s] %-13s %-11s скрыто %5.1f %% z_sum %+7.2f z_dip %+7.2f z_χ² %+8.2f %s" % (
            tag, c["name"], c["status"], c["hidden_frac"] * 100, c["z_sum"], c["z_dip"], c["z_chi2"], c["reason"]
        ))

def start_time(path):
    text = pathlib.Path(path).read_text("utf-8")
    m = re.search(r"<StartTime>([^<]*)</StartTime>", text)
    return m.group(1) if m else ""

def step_f(tmp, record, interp, crit):
    if not BG8.exists() or not SM8.exists():
        record("F пара 8192", False, "нет файлов")
        return
        
    r, e, s = chain(BG8, SM8, tmp / "f8192")
    
    # F1
    ok_f1 = r["n_ch"] == 8192
    record("F1 8192 каналы", ok_f1, "n_ch=%d" % r["n_ch"])
    
    # F2
    ok_f2 = all(a["ok"] for a in r["auto"]["anchors"].values())
    record("F2 якоря", ok_f2, "все якоря ok" if ok_f2 else "есть не ok якоря")
    
    # F4
    ok_f4 = all(f["n_points"] == 8192 and f["all_int"] and f["match"] for f in s["files"])
    record("F4 файлы", ok_f4, "все файлы корректны" if ok_f4 else "есть некорректные файлы")
    
    # F5
    ok_f5 = s["fwhm_ok"]
    record("F5 FWHM", ok_f5, "fwhm_ok=%s" % s["fwhm_ok"])
    
    # F3
    S, B, k, _, _ = _LOAD_INPUTS(BG8, SM8)
    lo, hi = r["h_json"]["auto"]["register_kw"]["mask"]
    m = slice(int(lo), int(hi) + 1)
    X = S - B
    pred = float(np.mean((X[m] + (1 - k) ** 2 * B[m]) / np.maximum(S[m] + k ** 2 * B[m], 1.0)))
    ratio = r["chi2ndf"] / pred
    ok_f3 = crit["F3"]["ratio_lo"] <= ratio <= crit["F3"]["ratio_hi"]
    detail_f3 = "χ²/ndf %.4f, ожидание вложенной пары %.4f, отношение %.3f (допуск %.2f..%.2f)" % (
        r["chi2ndf"], pred, ratio, crit["F3"]["ratio_lo"], crit["F3"]["ratio_hi"]
    )
    record("F3 8192 вложенная: χ²/ndf = ожидание для вложенной пары", ok_f3, detail_f3)
    
    # F6
    n_dep = sum(1 for c in s["struct"] if c["dependent"])
    ok_f6 = (n_dep >= crit["F6"]["min_dependent_windows"] and 
             not s["struct_summary"]["struct_ok"] and 
             not s["all_ok"])
    detail_f6 = "окон с флагом %d, struct_ok %s, all_ok %s; %s" % (
        n_dep, s["struct_summary"]["struct_ok"], s["all_ok"], fmt_windows(s["struct"])
    )
    record("F6 8192 вложенная: вердикт окон -- флаг зависимой пары, итог красный", ok_f6, detail_f6)
    
    # F7
    ind = pl.nested_pair(str(BG8_IND), str(SM8)) if BG8_IND.exists() else ""
    same_start = start_time(BG8) != "" and start_time(BG8) == start_time(SM8)
    ok_f7 = bool(r["nested"]) and same_start and ind == ""
    detail_f7 = "вложенная: %s; StartTime обоих %s; независимая пара: %s" % (
        r["nested"][:80] or "нет", start_time(BG8), ind[:60] or "не помечена"
    )
    record("F7 8192: вложенная пара распознана, независимая -- нет", ok_f7, detail_f7)
    
    interp.append("пара 8192: χ²/ndf %.3f при ожидании вложенной пары %.3f" % (r["chi2ndf"], pred))
    
    # F8
    if BG8_IND.exists():
        r2, e2, s2 = chain(BG8_IND, SM8, tmp / "f8192ind")
        hm = crit["F8"]["hide_max"]
        viol = [c["name"] for c in list(s2["struct"]) + list(r2["structure"]) 
                if c["evaluated"] and max(c["hidden_frac_l"], c["hidden_frac_r"]) > hm]
        conv = s2["converged"]
        ok_f8 = (not viol and 
                 all("hidden_frac" in c for c in s2["struct"]) and 
                 (conv or not s2["all_ok"]))
        detail_f8 = "сходимость %s, all_ok %s, оценено %d из %d, нарушений правила перекрытия %s; %s" % (
            conv, s2["all_ok"], s2["struct_summary"]["n_eval"], s2["struct_summary"]["n"], 
            viol or "нет", fmt_windows(s2["struct"])
        )
        record("F8 8192 независимый фон: окна, скрытые линиями образца, не оцениваются; несошедшаяся регистрация не зелёная", 
               ok_f8, detail_f8)
        print_table("F8-ОКНА", s2["struct"])
        if not conv:
            interp.append("независимая пара 8192: устойчивая регистрация не сошлась")

def step_g(tmp, record, interp, crit, app=None, w=None, wait_job=None):
    if not BG_SOIL.exists() or not SM_SOIL.exists():
        record("G грунт", False, "нет файлов")
        return
        
    r, e, s = chain(BG_SOIL, SM_SOIL, tmp / "soil")
    
    print_table("G-РОБАСТ", r["structure"])
    print_table("G-ПРОСТАЯ", r["structure_plain"])
    
    sb = r["robust"].get("same_basis")
    if sb:
        print("[G] на одних каналах (%d): без модели простая %.3f -> устойчивая %.3f; с моделью %.3f -> %.3f" % (
            sb["n_ch"], sb["no_model"][0], sb["no_model"][1], sb["with_model"][0], sb["with_model"][1]
        ))
        
    # G1
    ok_g1 = r["robust"]["mode"] == "robust" and r["robust"].get("converged", False)
    excluded = r["robust"]["excluded"]
    detail_g1 = "χ²/ndf простой %.3f -> итог %.3f; исключено %s" % (
        r["robust"]["chi2ndf_plain"], r["robust"]["chi2ndf"],
        ", ".join("(%d, %d)" % (x[0], x[1]) for x in excluded) if excluded else "нет"
    )
    record("G1 грунт: устойчивая регистрация включилась и сошлась", ok_g1, detail_g1)
    
    # G2
    def n_ok(cs):
        return sum(1 for c in cs if c["evaluated"] and c["ok"])
        
    def win(cs, nm):
        return next((c for c in cs if c["name"] == nm), None)
        
    kr = win(r["structure"], "La_sum_K40")
    kp = win(r["structure_plain"], "La_sum_K40")
    
    ok_g2 = (n_ok(r["structure"]) - n_ok(r["structure_plain"]) >= crit["G2"]["min_gain"] and 
             kr is not None and kp is not None and 
             kr["evaluated"] and kr["ok"] and 
             kp["evaluated"] and not kp["ok"])
             
    detail_g2 = "в норме: устойчивая %d, простая %d; La_sum_K40 устойчивая %s, простая %s" % (
        n_ok(r["structure"]), n_ok(r["structure_plain"]), 
        kr["status"] if kr else "нет", kp["status"] if kp else "нет"
    )
    record("G2 грунт, одна основа: окон в норме у устойчивой больше, сумм-пик в норме только у неё", ok_g2, detail_g2)
    
    # G3
    undo = masking_off()
    try:
        rm, em, sm2 = chain(BG_SOIL, SM_SOIL, tmp / "soil_nomask", robust="on")
    finally:
        undo()
        
    ok_g3 = not rm["structure_summary"]["struct_ok"]
    detail_g3 = "без исключения линий: исключено окон %d, оценено %d, в норме %d, struct_ok %s" % (
        len(rm["robust"]["excluded"]), rm["structure_summary"]["n_eval"], 
        rm["structure_summary"]["n_ok"], rm["structure_summary"]["struct_ok"]
    )
    record("G3 мутация «линии образца не учитывать» красит итог окон", ok_g3, detail_g3)
    
    # G4
    kk = kr["k40"] if kr else None
    ok_g4 = bool(kk) and np.isfinite(kk["area"]) and kk["area_err"] > 0
    if ok_g4:
        detail_g4 = "площадь %.0f ± %.0f, σ %.2f кан" % (kk["area"], kk["area_err"], kk["sigma"])
    else:
        detail_g4 = "нет"
    record("G4 грунт: K-40 под сумм-пиком -- отдельная гауссиана в окне, площадь ± σ", ok_g4, detail_g4)
    
    # G5
    if app is None or w is None or wait_job is None:
        record("G5 грунт: подписи на графике окна -- реперы пары, не лютециевые", False, "нет окна")
        return
        
    w.bg_edit.setText(str(BG_SOIL))
    w.sample_edit.setText(str(SM_SOIL))
    w._update_pair()
    
    out = tmp / "ui_soil"
    out.mkdir(parents=True, exist_ok=True)
    w.out_dir.setText(str(out))
    w.out_name.setText("soil")
    
    w._do_register()
    wait_job()
    w._do_ecal()
    wait_job()
    w._do_subtract()
    wait_job()
    
    drawn = [it.label.toPlainText() for it in w.plot.plot.getPlotItem().items 
             if isinstance(it, pg.InfiniteLine) and getattr(it, "label", None) is not None]
             
    expected = set()
    if w.ecal_res:
        expected |= {p["name"] for p in w.ecal_res.get("points", [])}
        expected |= {q["name"] for q in w.ecal_res.get("low", {}).get("residuals", [])}
        
    bad = [d for d in drawn if d.startswith(tuple(crit["G5"]["forbidden_prefixes"]))]
    
    ok_g5 = len(drawn) > 0 and set(drawn) <= expected and not bad
    detail_g5 = "нарисовано: %s; ожидалось из реперов пары: %s; запрещённых: %s" % (
        ", ".join(drawn), ", ".join(sorted(expected)), bad or "нет"
    )
    record("G5 грунт: подписи на графике окна -- реперы пары, не лютециевые", ok_g5, detail_g5)

def step_h(tmp, record, interp, app, w, wait_job):
    """
    Тест H: сценарий «Новый спектр» в UI.
    Проверяет диалог сохранения, сброс состояния и работу новой пары.
    """
    from PySide6 import QtWidgets

    # Предусловие: w.h должен быть установлен
    if w.h is None:
        record("H «Новый спектр»", False, "предусловие не выполнено: w.h is None")
        return

    # (a) Несохранённый результат, пользователь отвечает «Нет»
    w.saved = False
    orig_question = QtWidgets.QMessageBox.question
    QtWidgets.QMessageBox.question = staticmethod(lambda *a, **k: QtWidgets.QMessageBox.No)
    try:
        res = w._new_spectrum()
    finally:
        QtWidgets.QMessageBox.question = orig_question

    record(
        "H1 «Новый спектр»: при несохранённом результате спрашивает и при «Нет» ничего не сбрасывает",
        res is False and w.h is not None and w.bg_edit.text() != "",
        "res=%s, h=%s, bg=%s" % (res, w.h is not None, w.bg_edit.text() != "")
    )

    # (b) Пользователь отвечает «Да»; в поле ручной поправки -- устаревшее значение прошлой пары (, раунд 4)
    w.anchor_box.manual_edit.setText(load_criteria()["H3"]["stale_manual"])
    QtWidgets.QMessageBox.question = staticmethod(lambda *a, **k: QtWidgets.QMessageBox.Yes)
    try:
        res = w._new_spectrum()
    finally:
        QtWidgets.QMessageBox.question = orig_question
        app.processEvents()

    # Собираем состояние после сброса
    state = {
        "bg": w.bg_edit.text(),
        "sm": w.sample_edit.text(),
        "h": w.h,
        "ecal": w.ecal,
        "ecal_res": w.ecal_res,
        "net": w.net,
        "tab": w.tabs.currentIndex(),
        "ecal_rows": w.ecal_table.rowCount(),
        "plot_items": len(w.plot.plot.listDataItems()),
        "reg_html": w.reg_view.toPlainText().strip(),
        "pair": w.pair_label.text(),
        "manual_text": w.anchor_box.manual_edit.text()
    }

    record(
        "H2 «Новый спектр»: файлы и результаты сброшены, вкладки в исходном",
        res is True and
        state["bg"] == "" and
        state["sm"] == "" and
        state["h"] is None and
        state["ecal"] is None and
        state["ecal_res"] is None and
        state["net"] is None and
        state["tab"] == 0 and
        state["ecal_rows"] == 0 and
        state["plot_items"] == 0 and
        state["reg_html"] == "" and
        state["manual_text"] == "",
        str({k: (v if not isinstance(v, (list, dict)) else type(v).__name__) for k, v in state.items()})
    )

    # (c) Новая пара работает без перезапуска
    w.bg_edit.setText(str(BG))
    w.sample_edit.setText(str(SM))
    w._update_pair()
    w._do_register()
    wait_job()
    h_ref = backend_pair.run_register(str(BG), str(SM))["h"]   # эталон -- чистый запуск

    record(
        "H3 после «Новый спектр» новая пара совмещается так же, как в чистом запуске",
        w.h is not None and list(w.h) == list(h_ref),
        "h окна %s, h чистого запуска %s" % (w.h, h_ref)
    )


def run(tmp, record, interp, app=None, w=None, wait_job=None, only="FGH"):
    crit = load_criteria()
    undo = apply_mutation(MUT)
    if MUT:
        print("[MUT] " + MUT)
        
    try:
        for letter in "FGH":
            if letter not in only:
                continue
                
            if letter == "F" and not (BG8 and SM8 and BG8_IND):
                print("[SKIP] F: не заданы ISGUI_PAIR8_BG, ISGUI_PAIR8_SAMPLE, ISGUI_PAIR8_BG_IND")
            elif letter == "G" and not (BG_SOIL and SM_SOIL):
                print("[SKIP] G: не заданы ISGUI_SOIL_BG, ISGUI_SOIL_SAMPLE")
            elif letter == "F":
                try:
                    step_f(tmp, record, interp, crit)
                except Exception:
                    record("F выполнен", False, traceback.format_exc())
                    
            elif letter == "G":
                try:
                    step_g(tmp, record, interp, crit, app, w, wait_job)
                except Exception:
                    record("G выполнен", False, traceback.format_exc())
                    
            elif letter == "H":
                if app and w and wait_job:
                    try:
                        step_h(tmp, record, interp, app, w, wait_job)
                    except Exception:
                        record("H выполнен", False, traceback.format_exc())
    finally:
        undo()
