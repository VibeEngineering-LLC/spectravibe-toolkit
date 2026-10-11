"""Вкладка «Шкала по пробе»: щелчок по пику → подгонка → диалог (линия из ядерных данных / ручная энергия) →
таблица; правка и удаление репера; «Шкала по пробе» в главном окне заполняет совмещение и шкалу. Сгенерировано
gen_code.py (qwen3.8:27b) по спецификации, 2026-10-10. Запуск: python tests/test_peakpicker.py (код 0 -- зелёный) или pytest."""

import os, sys, re, math, pathlib, tempfile
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path[0:0] = [str(ROOT / "src"), str(ROOT.parents[1] / "scripts")]
import numpy as np
from PySide6 import QtCore, QtWidgets
from isgui import peakpicker as pp
FIX = ROOT.parents[1] / "tests" / "fixtures" / "intrinsic_shape"
BG = FIX / "labr3_background.xml"
SAMPLE = FIX / "labr3_lu176_sample.xml"


def _app():
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


class _Ev:
    def __init__(self, pos):
        self._pos = pos

    def button(self):
        return QtCore.Qt.LeftButton

    def scenePos(self):
        return self._pos


def _click(panel, ch):
    vb = panel.plot.getPlotItem().getViewBox()
    panel.plot.resize(900, 400)
    panel.plot.show()
    _app().processEvents()
    vb.setRange(xRange=(0, panel.n), padding=0)
    _app().processEvents()
    y = vb.viewRange()[1]
    pos = vb.mapViewToScene(QtCore.QPointF(float(ch), (y[0] + y[1]) / 2.0))
    panel._on_click(_Ev(pos))


INITIAL_ROWS = []


def _dialog_auto(choice):
    def exec_(self):
        INITIAL_ROWS.append(self.table.currentRow())   # строка не выбирается молча
        if choice == "lib":
            self.radio_lib.setChecked(True)
            self.table.selectRow(0)
        else:
            self.e_edit.setText(choice)
        self._accept()
        return QtWidgets.QDialog.Accepted if self.result_value else QtWidgets.QDialog.Rejected
    return exec_


def _sample_with_cs(tmp):
    from gamma.calibration.intrinsic_shape.core import load_spectrum, poly_inverse
    d = load_spectrum(BG)
    B = d["counts"].astype(float)
    x = np.arange(len(B), dtype=float)
    c = float(poly_inverse(list(d["cal"]), 661.657))
    S = B + 60000.0 / (math.sqrt(2 * math.pi) * 4.0) * np.exp(-0.5 * ((x - c) / 4.0) ** 2)
    t = BG.read_text("utf-8")
    vals = iter(np.round(S).astype(int))
    p = pathlib.Path(tmp) / "sample_cs.xml"
    p.write_text(re.sub(r"<DataPoint>[^<]*</DataPoint>", lambda m: "<DataPoint>%d</DataPoint>" % next(vals), t), encoding="utf-8")
    return (p, c)


def test_click_fit_dialog_table_edit_delete(monkeypatch=None):
    _app()
    tmp = tempfile.mkdtemp()
    p, c = _sample_with_cs(tmp)
    panel = pp.PeakPickerPanel()
    panel.set_sample(str(p))
    msgs = []
    panel.message.connect(msgs.append)
    old = pp.ReperDialog.exec
    try:
        pp.ReperDialog.exec = _dialog_auto("lib")
        _click(panel, c + 1.5)
        assert len(panel.repers) == 1
        r = panel.repers[0]
        assert r["nuclide"] == "Cs-137"
        assert abs(r["E"] - 661.657) < 0.01
        assert r["source"] == "библиотека"
        assert abs(r["ch"] - c) < 0.5
        assert panel.table.rowCount() == 1
        assert len(panel.marks) == 1
        _click(panel, c + 1.0)   # повтор той же линии не добавляется
        assert len(panel.repers) == 1 and any("повтор" in m for m in msgs)
        assert INITIAL_ROWS and all(r == -1 for r in INITIAL_ROWS)

        pp.ReperDialog.exec = _dialog_auto("700.5")
        panel._edit_row(0)
        assert panel.repers[0]["source"] == "ручная"
        assert abs(panel.repers[0]["E"] - 700.5) < 1e-9
        assert abs(panel.repers[0]["ch"] - r["ch"]) < 1e-9
        assert "ручных 1" in panel.summary.text()

        panel.table.setCurrentCell(0, 1)   # как щелчок пользователя по ячейке (не selectRow)
        panel._delete()
        assert panel.repers == []
        assert panel.table.rowCount() == 0
        assert panel.marks == []
        import pyqtgraph as pg   # линии реально сняты с графика, а не только из списка
        assert not any(isinstance(i, pg.InfiniteLine) for i in panel.plot.getPlotItem().items)

        pp.ReperDialog.exec = _dialog_auto("lib")
        _click(panel, 5.0)
        assert panel.repers == []
        assert any("не подогнан" in m for m in msgs)
    finally:
        pp.ReperDialog.exec = old


def test_main_window_sample_scale_fills_h_and_ecal():
    _app()
    from isgui.main_window import MainWindow
    from isgui import backend_sample   # архитектура «разность на шкале пробы»
    w = MainWindow()
    names = [w.tabs.tabText(i) for i in range(w.tabs.count())]
    assert "2б. Шкала по пробе" in names
    w.bg_edit.setText(str(BG))
    w.sample_edit.setText(str(SAMPLE))
    w._update_pair()
    assert w.picker.counts is not None and w.picker.n == 2048
    r = backend_sample.run_sample_cal(str(BG), str(SAMPLE), [])
    w._show_sample_cal(r)
    assert w.sample_res is r and w.h == [0.0, 1.0, 0.0]
    assert w.ecal == r["coeffs"] and not r["decision"]["replaced"]   # без реперов -- калибровка файла пробы
    assert w.ecal_table.rowCount() == len(r["branch_rows"]) == 0
    txt = w.sc_view.toPlainText()
    assert "Шкала XML" in txt and "Калибровка файла пробы" in txt and "Фон переведён на сетку пробы" in txt
    assert w.fw_points_edit.text().strip() != ""
    assert all(b in w.buttons for b in w.picker.buttons)
    assert not w.apply_scale_btn.isEnabled() and w.sc_table.rowCount() == 0   # без рекомендации -- файл
    rec = dict(r, decision=dict(r["decision"], recommend=True, why=""), branch_rows=[("Tl208_2615", "репер", 4978.3, 2614.5, 15.0, 0.4, 1.1)])
    w._show_sample_cal(rec)
    assert w.apply_scale_btn.isEnabled() and w.ecal == r["decision"]["coeffs_file"]
    assert w.sc_table.rowCount() == 1 and w.sc_table.columnCount() == len(r["branch_headers"])
    assert "Рекомендована шкала по реперам" in w.sc_view.toPlainText()
    w._show_sample_cal(r)
    from isgui import peakpicker as pp2
    w.picker.repers = [{"ch": 100.0, "E": 150.0}, {"ch": 500.0, "E": 750.0}]
    assert abs(w.picker._energy(300.0) - 450.0) < 1e-6   # диалог: шкала по 2 реперам, не калибровка файла
    w.picker.repers = []
    h = backend_sample.decision_html(dict(r, warnings=["точки ПШПВ пробы не дают годной модели"]))
    assert "#fff3cd" in h and "точки ПШПВ пробы" in h
    w.sample_edit.setText(str(BG)); w._update_pair()   # другая пара
    assert w.sample_res is None   # смена пары сбрасывает результат режима «по пробе»

def test_reper_scale_button_state():
    """кнопка гаснет после применения, при смене пары, «Новый спектр», запуске задачи; без
    рекомендации запись не запускается; «что нужно для замены» -- в окне."""
    _app()
    from isgui.main_window import MainWindow
    from isgui import backend_sample
    w = MainWindow(); w.bg_edit.setText(str(BG)); w.sample_edit.setText(str(SAMPLE)); w._update_pair()
    r = backend_sample.run_sample_cal(str(BG), str(SAMPLE), [])
    assert "для замены нужно" in backend_sample.decision_html(r)
    rec = dict(r, decision=dict(r["decision"], recommend=True, why=""), branch_rows=[("x", "репер", 1.0, 2.0, 0.1, 0.2, 0.3)])
    w._show_sample_cal(dict(rec, decision=dict(rec["decision"], replaced=True)))
    assert not w.apply_scale_btn.isEnabled()                          # уже записана
    calls, errs = [], []
    w._run, w._error = (lambda *a, **k: calls.append(a)), errs.append
    w._show_sample_cal(r); w._do_sample_cal(use_reper_scale=True)
    assert not calls and errs                                          # без рекомендации -- не запускается
    w._show_sample_cal(rec); assert w.apply_scale_btn.isEnabled()
    w.sample_edit.setText(str(BG)); w._update_pair()
    assert not w.apply_scale_btn.isEnabled() and w.sc_table.rowCount() == 0 and w.sc_view.toPlainText() == ""
    w._show_sample_cal(rec); w._new_spectrum(confirm=False); assert not w.apply_scale_btn.isEnabled()
    del w._run
    w._show_sample_cal(rec); w._run("тест", lambda: None, (), {}, lambda res: None)
    assert not w.apply_scale_btn.isEnabled(); w.job.wait()


if __name__ == "__main__":
    test_click_fit_dialog_table_edit_delete()
    test_main_window_sample_scale_fills_h_and_ecal()
    test_reper_scale_button_state()
    print("OK")
