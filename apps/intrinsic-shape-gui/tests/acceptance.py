import os, sys, re, json, time, shutil, pathlib, tempfile, subprocess, traceback
sys.stdout.reconfigure(encoding="utf-8")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path[0:0] = [str(ROOT / "src"), str(ROOT.parents[1] / "scripts")]
import numpy as np
from isgui import backend, backend_pair
from PySide6 import QtWidgets as _QW   # модальные окна ошибок в приёмке не блокируют, а пишутся в журнал
_QW.QMessageBox.critical = staticmethod(lambda *a, **k: print("[DIALOG] %s" % (a[1:3],)))
_QW.QMessageBox.warning = staticmethod(lambda *a, **k: print("[DIALOG] %s" % (a[1:3],)))
FW = backend.default_fwhm_points()
RESULTS = []   # list of (name, ok_bool, detail_str)
INTERP = []    # strings for the «требует толкования» block
def record(name, ok, detail): RESULTS.append((name, bool(ok), detail)); print("[%s] %s: %s" % ("OK" if ok else "FAIL", name, detail))
def datapoints(path): return [v for v in re.findall(r"<DataPoint>([^<]*)</DataPoint>", pathlib.Path(path).read_text("utf-8"))]

FIX = ROOT.parents[1] / "tests" / "fixtures" / "intrinsic_shape"
BG = FIX / "labr3_background.xml"
SM = FIX / "labr3_lu176_sample.xml"

def step_b(tmp):
    ref = tmp / "ref"
    ref.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, PYTHONPATH=str(ROOT.parents[1] / "scripts"), PYTHONIOENCODING="utf-8")
    
    def cli(*args):
        proc = subprocess.run([sys.executable, "-m", "gamma.calibration.intrinsic_shape", *map(str, args)], env=env, capture_output=True, text=True, encoding="utf-8")
        if proc.returncode != 0:
            raise RuntimeError("CLI %s: rc=%d\n%s" % (args[0], proc.returncode, proc.stderr))
        return proc.stdout

    cli("register", "--bg", BG, "--sample", SM, "--out", ref / "h.json")
    cli("ecal", "--bg", BG, "--sample", SM, "--h", ref / "h.json", "--out", ref / "ecal.json")
    cli("subtract", "--bg", BG, "--sample", SM, "--h", ref / "h.json", "--ecal", ref / "ecal.json", "--fwhm-points", *["%.2f:%.2f" % p for p in FW], "--out-zeroed", ref / "net.xml", "--out-signed", ref / "net_signed.xml")

    gui = tmp / "gui"
    gui.mkdir(parents=True, exist_ok=True)
    
    # B идёт через backend_pair -- тот же обработчик, что у окна
    r = backend_pair.run_register(str(BG), str(SM))
    e = backend_pair.run_ecal(str(BG), str(SM), r["h"])
    s = backend_pair.run_subtract(str(BG), str(SM), r["h"], e, [tuple(p) for p in e["fwhm_points"]],
                                  str(gui / "net.xml"), str(gui / "net_signed.xml"), reg_json=r["h_json"])
    ref_dir = os.environ.get("ISGUI_REF_DIR", "")   # эталон -- вне репозитория, через переменную
    if not ref_dir:
        print("[SKIP] B0: не задана ISGUI_REF_DIR (папка замороженного выхода net.xml/net_signed.xml)")
    for nm in (("net.xml", "net_signed.xml") if ref_dir else ()):   # регрессия против замороженного выхода e5498c7 (до раунда 4)
        a, b = datapoints(pathlib.Path(ref_dir) / nm), datapoints(gui / nm)
        nd = sum(1 for x, y in zip(a, b) if x != y) + abs(len(a) - len(b))
        record("B0 регрессия DataPoint %s = e5498c7" % nm, nd == 0 and len(a) > 0, "точек %d / %d, отличий %d" % (len(a), len(b), nd))

    href = json.loads((ref / "h.json").read_text("utf-8"))["c"]
    diff_h = max(abs(x - y) for x, y in zip(r["h"], href))
    record("B1 h совпадает с CLI", diff_h == 0.0, "max|Δh| = %g" % diff_h)

    eref = json.loads((ref / "ecal.json").read_text("utf-8"))["coeffs"]
    diff_e = max(abs(x - y) for x, y in zip(e["coeffs"], eref))
    record("B2 шкала совпадает с CLI", diff_e == 0.0, "max|Δe| = %g" % diff_e)

    for ref_name, gui_name in [("net.xml", "net.xml"), ("net_signed.xml", "net_signed.xml")]:
        ref_path = ref / ref_name
        gui_path = gui / gui_name
        a = datapoints(ref_path)
        b = datapoints(gui_path)
        n_diff = sum(1 for x, y in zip(a, b) if x.strip() != y.strip()) + abs(len(a) - len(b))
        record("B3 DataPoint %s = CLI" % gui_path.name, n_diff == 0 and len(a) > 0, "точек %d / %d, отличий %d" % (len(a), len(b), n_diff))
        non_int = sum(1 for x in b if not re.fullmatch(r"-?\d+", x.strip()))
        record("B4 все DataPoint целые: %s" % gui_path.name, non_int == 0, "нецелых %d из %d" % (non_int, len(b)))

    record("B5 блок ПШПВ годен", s["fwhm_ok"] and all(f["has_fwhm"] for f in s["files"]), "запас за N-1 = %.1f кан, блок в файлах: %s" % (s["fwhm_margin"], [f["has_fwhm"] for f in s["files"]]))
    record("B6 остаток z в норме", s["z_ok"], "z = %+.3f (порог |z| <= %g), сумма %.1f +- %.1f" % (s["z"], s["z_limit"], s["sum"], s["sigma"]))
    record("B7 совпадение после перечитывания", all(f["match"] for f in s["files"]), str([(f["name"], f["match"]) for f in s["files"]]))
    record("B8 итог обработчика all_ok", s["all_ok"], "all_ok=%s" % s["all_ok"])

    INTERP.extend("совмещение: " + it for it in r["interp"])
    if abs(s["z"]) > 1.5:
        INTERP.append("z = %+.2f близко к порогу" % s["z"])

    return (r, e)

def step_c(tmp, r, e):
    from gamma.calibration.intrinsic_shape import ops
    orig = ops.reg_apply
    def mutated(c, S, B, k, *a, **kw):
        return orig(c, S, B, k * 1.01, *a, **kw)
    ops.reg_apply = mutated
    try:
        mut_dir = tmp / "mut"
        mut_dir.mkdir(parents=True, exist_ok=True)
        m = backend_pair.run_subtract(str(BG), str(SM), r["h"], e, [tuple(p) for p in e["fwhm_points"]],
                                      str(mut_dir / "net.xml"), str(mut_dir / "net_signed.xml"), reg_json=r["h_json"])
    finally:
        ops.reg_apply = orig

    record("C1 мутация k*1.01: проверка остатка краснеет", not m["z_ok"], "z = %+.2f (без мутации было в норме)" % m["z"])
    record("C2 мутация k*1.01: итог all_ok краснеет", not m["all_ok"], "all_ok=%s" % m["all_ok"])

    ctl_dir = tmp / "ctl"
    ctl_dir.mkdir(parents=True, exist_ok=True)
    c = backend_pair.run_subtract(str(BG), str(SM), r["h"], e, [tuple(p) for p in e["fwhm_points"]],
                                  str(ctl_dir / "net.xml"), str(ctl_dir / "net_signed.xml"), reg_json=r["h_json"])
    record("C3 контроль: после снятия мутации снова зелёно", c["z_ok"], "z = %+.3f" % c["z"])

def step_d(tmp):
    bad = tmp / "bad"
    bad.mkdir(parents=True, exist_ok=True)
    
    (bad / "notxml.xml").write_text("это не спектр, просто текст\n" * 10, encoding="utf-8")
    bg_text = BG.read_text("utf-8")
    (bad / "broken.xml").write_text(bg_text[:3000], encoding="utf-8")
    (bad / "empty.xml").write_text("", encoding="utf-8")
    shutil.copy(str(BG), str(bad / "spectrum.txt"))
    missing_path = bad / "missing.xml"

    bad_paths = [bad / "notxml.xml", bad / "broken.xml", bad / "empty.xml", bad / "spectrum.txt", missing_path]
    for p in bad_paths:
        try:
            backend.pair_info(str(p), str(SM))
            record("D %s" % p.name, False, "исключения нет")
        except Exception as exc:
            msg = backend.friendly(exc)
            ok = (isinstance(exc, backend.UserError) or not msg.startswith("Внутренняя ошибка")) and "Traceback" not in msg
            record("D понятное сообщение: %s" % p.name, ok, msg)

    from PySide6 import QtCore, QtWidgets
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from isgui.widgets import Job
    got = {}
    job = Job(backend.run_register, str(bad / "broken.xml"), str(SM))
    job.failed.connect(lambda msg, tb: got.update(msg=msg))
    job.done.connect(lambda res: got.update(res=res))
    job.start()
    t0 = time.time()
    while ("msg" not in got and "res" not in got) and time.time() - t0 < 60:
        app.processEvents()
        time.sleep(0.02)
    job.wait(5000)
    record("D фоновый поток: ошибка доставлена в окно", "msg" in got and not got["msg"].startswith("Внутренняя ошибка"), got.get("msg", "нет сообщения"))

    from isgui.main_window import MainWindow
    w = MainWindow()
    w.bg_edit.setText(str(bad / "notxml.xml"))
    w.sample_edit.setText(str(SM))
    w._update_pair()
    text = w.pair_label.text()
    record("D окно: сообщение в строке файлов", "✘" in text and "Внутренняя" not in text, re.sub("<[^>]+>", "", text))
    
    return (app, w)

def step_e(tmp, app, w):
    def wait_job(timeout=600):
        t0 = time.time()
        while time.time() - t0 < timeout:
            app.processEvents()
            if w.job is None or (w.job.isFinished() and not w.progress.isVisible()):
                return True
            time.sleep(0.02)
        return False

    w.bg_edit.setText(str(BG))
    w.sample_edit.setText(str(SM))
    w._update_pair()
    out = tmp / "ui"
    out.mkdir(parents=True, exist_ok=True)
    w.out_dir.setText(str(out))
    w.out_name.setText("net")

    w._do_register()
    ok1 = wait_job()
    w._do_ecal()
    ok2 = wait_job()
    w._do_subtract()
    ok3 = wait_job()

    record("E1 цепочка кнопок отработала", ok1 and ok2 and ok3 and w.net is not None and w.net["all_ok"], "h=%s, шкала=%s, all_ok=%s" % (w.h is not None, w.ecal is not None, w.net["all_ok"] if w.net else None))
    
    a = datapoints(out / "net.xml")
    b = datapoints(tmp / "gui" / "net.xml")
    record("E2 файл из окна = файл обработчика", a == b and len(a) > 0, "точек %d, отличий %d" % (len(a), sum(1 for x, y in zip(a, b) if x != y)))
    
    n_items = len(w.plot.plot.listDataItems())
    record("E3 график построен", n_items >= 3, "кривых %d" % n_items)

def wait_ui(app, w, timeout=900):
    """дождаться конца фоновой задачи окна и доставки сигнала done в UI-поток."""
    t0 = time.time()
    while time.time() - t0 < timeout:
        app.processEvents()
        if w.job is None or (w.job.isFinished() and not w.progress.isVisible()):
            break
        time.sleep(0.02)
    for _ in range(50):
        app.processEvents()
        time.sleep(0.02)
    return time.time() - t0 < timeout


def main():
    tmp = pathlib.Path(tempfile.mkdtemp(prefix="isgui_acc_"))
    print("[TMP]", tmp)
    t0 = time.time()
    
    import hashlib   # критерии -- файлом с sha до прогона; ISGUI_ACC_ONLY -- подмножество шагов
    print("[CRITERIA] tests/acceptance_criteria.json sha256 %s" % hashlib.sha256((ROOT / "tests" / "acceptance_criteria.json").read_bytes()).hexdigest())
    only = os.environ.get("ISGUI_ACC_ONLY", "BCDEFGH")
    r, e = None, None
    if "B" in only or "C" in only:
        try:
            r, e = step_b(tmp)
        except Exception:
            record("B выполнен", False, traceback.format_exc())
    if "C" in only:
        if r is not None and e is not None:
            try:
                step_c(tmp, r, e)
            except Exception:
                record("C выполнен", False, traceback.format_exc())
        else:
            record("C выполнен", False, "пропущено из-за ошибки в B")

    app, w = None, None
    if any(x in only for x in "DEGH"):
        try:
            app, w = step_d(tmp)
        except Exception:
            record("D выполнен", False, traceback.format_exc())
    if any(x in only for x in "EH"):
        if app is not None and w is not None:
            try:
                step_e(tmp, app, w)
            except Exception:
                record("E выполнен", False, traceback.format_exc())
        else:
            record("E выполнен", False, "пропущено из-за ошибки в D")

    import acceptance_pair   # пара 8192, грунт (линии образца), кнопка «Новый спектр»
    acceptance_pair.run(tmp, record, INTERP, app, w, (lambda: wait_ui(app, w)) if w is not None else None,
                        only="".join(x for x in "FGH" if x in only))

    print("[ТРЕБУЕТ ТОЛКОВАНИЯ]")
    if INTERP:
        for line in INTERP:
            print("  " + line)
    else:
        print("  (нет пунктов)")
        
    n_fail = sum(1 for x in RESULTS if not x[1])
    print("ИТОГ ПРИЁМКИ: проверок %d, провалов %d, время %.1f с" % (len(RESULTS), n_fail, time.time() - t0))
    return 0 if n_fail == 0 else 1

if __name__ == "__main__":
    sys.exit(main())
