"""Разность на лог-шкале: обе ветви (>0 и |<0|) и предупреждение об обнулении (правка Claude, 2026-10-10).
Запуск: python tests/test_net_plot.py (код 0 -- зелёный) или через pytest."""
import os, sys, pathlib
import numpy as np
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path[0:0] = [str(ROOT / "src"), str(ROOT.parents[1] / "scripts")]
from PySide6 import QtWidgets
from isgui import widgets


def test_log_plot_draws_negative_branch_and_warns():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    ch = np.arange(2048)
    S, bgk = np.full(2048, 100.0), np.where(ch % 2 == 0, 150.0, 60.0)   # фон·k выше образца в половине каналов
    data = {"channels": ch, "energy": ch * 1.4, "sample": S, "bg_k": bgk, "signed": S - bgk,
            "zeroed": np.maximum(S - bgk, 0), "live_s": 100.0, "live_out": 100.0}
    p = widgets.SpectrumPlot()
    p.log_check.setChecked(True); p.mode_combo.setCurrentIndex(0); p.unit_combo.setCurrentIndex(0)
    p.set_data(data, [])
    items = {it.name(): it for it in p.plot.listDataItems()}
    n_pos, n_neg = (int(np.isfinite(items[k].yData).sum()) for k in ("Разность", widgets.NEG_NAME))
    print("кривых: %s; точек >0: %d, |<0|: %d; предупреждение: %s" % (list(items), n_pos, n_neg, p.zeroed_warning))
    assert n_pos > 0 and n_neg > 0 and "50 %" in p.zeroed_warning and p.note.text() == p.zeroed_warning


if __name__ == "__main__":
    test_log_plot_draws_negative_branch_and_warns(); print("OK")
