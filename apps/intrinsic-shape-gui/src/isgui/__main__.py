"""Точка входа: python -m isgui (из src) или собранный exe. ISGUI_SMOKE_FILE -- приёмка: время до показа окна, выход."""
import os, sys, time, pathlib
_T0 = time.time()
if not getattr(sys, "frozen", False):   # разработка: модуль toolkit -- scripts/ корня репозитория
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[4] / "scripts"))
else:   # ядерные данные toolkit (nuclides.json, X-ray CSV) лежат в сборке в data/ рядом с модулями
    import gamma.data
    gamma.data.DATA_DIR = pathlib.Path(getattr(sys, "_MEIPASS", pathlib.Path(sys.executable).parent)) / "data"
if os.environ.get("ISGUI_SMOKE_CALC"):   # приёмка сборки: "фон|образец|папка" -> расчёт без окна, итог в папка/calc.txt
    from isgui import backend_pair as bp; from gamma.calibration.intrinsic_shape import selftest as st
    bg, sm, out = os.environ["ISGUI_SMOKE_CALC"].split("|"); pathlib.Path(out).mkdir(parents=True, exist_ok=True)
    r = bp.run_register(bg, sm); e = bp.run_ecal(bg, sm, r["h"])
    s = bp.run_subtract(bg, sm, r["h"], e, e["fwhm_points"], out + "/net.xml", out + "/net_signed.xml", reg_json=r["h_json"])
    try:   # калибровка по пробе в сборке (реперы -- только видимые структуры LaBr3)
        from isgui import backend_sample as bs; c = bs.run_sample_cal(bg, sm, []); sc_txt = "samplecal_mode=%s visible=%s replaced=%s recommend=%s" % (c["mode"], ",".join(c["visible_sample"]), c["decision"]["replaced"], c["decision"]["recommend"])
    except Exception as exc:
        sc_txt = "samplecal_error=%s" % exc
    pathlib.Path(out, "calc.txt").write_text("all_ok=%s struct_ok=%s n_ch=%d low=%s selftest_preset=%s %s" % (
        s["all_ok"], s["struct_ok"], r["n_ch"], "low" in e, st.DEFAULT_CFG.is_file(), sc_txt), encoding="utf-8")
    sys.exit(0)
from PySide6 import QtCore, QtWidgets
from isgui.main_window import MainWindow
app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
win = MainWindow()
win.show()
smoke = os.environ.get("ISGUI_SMOKE_FILE")
if smoke:
    QtCore.QTimer.singleShot(0, lambda: pathlib.Path(smoke).write_text("shown %.3f visible=%s title=%s" % (
        time.time() - _T0, win.isVisible(), win.windowTitle()), encoding="utf-8"))
    if os.environ.get("ISGUI_SMOKE_PNG"):   # снимок окна для приёмки
        QtCore.QTimer.singleShot(150, lambda: win.grab().save(os.environ["ISGUI_SMOKE_PNG"]))
    QtCore.QTimer.singleShot(int(os.environ.get("ISGUI_SMOKE_HOLD_MS", "300")), app.quit)
sys.exit(app.exec())
