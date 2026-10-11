"""Счётчик «прошло mm:ss»: растёт во время операции 3 с, после остановки фиксируется (правка Claude, 2026-10-10).
Запуск: python tests/test_op_timer.py (код 0 -- зелёный) или через pytest."""
import os, sys, time, pathlib
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path[0:0] = [str(ROOT / "src"), str(ROOT.parents[1] / "scripts")]
from PySide6 import QtWidgets
from isgui import widgets


def pump(app, lbl, sec):
    """Мок длительной операции: UI-поток крутит события sec секунд, собирает показанные тексты метки."""
    t0, seen = time.monotonic(), []
    while time.monotonic() - t0 < sec:
        app.processEvents()
        time.sleep(0.05)
        seen.append(lbl.text())
    return seen


def test_timer_grows_then_freezes():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    lbl = QtWidgets.QLabel("")
    t = widgets.OpTimer(lbl)
    t.start("Совместить")
    seen = pump(app, lbl, 3.3)
    el = t.stop()
    final, after = lbl.text(), set(pump(app, lbl, 1.5))
    print("показано: %s; итог: %s; после: %s; прошло %.2f с" % (sorted(set(seen)), final, after, el))
    assert {"Совместить: прошло 00:01", "Совместить: прошло 00:02", "Совместить: прошло 00:03"} <= set(seen)
    assert final == "Совместить: прошло 00:03" and after == {final} and 3.0 <= el < 4.0


if __name__ == "__main__":
    test_timer_grows_then_freezes(); print("OK")
