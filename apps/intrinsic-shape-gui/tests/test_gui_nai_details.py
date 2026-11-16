# -*- coding: utf-8 -*-
"""Окно, метод NaI: нулевое живое время, кривая фон·k на графике, режим фон-фон с выбранным режимом совмещения,
русское сообщение о битом XML (offscreen). Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-11-16."""

import os, sys, time, pathlib
import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT = pathlib.Path(__file__).resolve().parents[1]
TK = ROOT.parents[1]
SRC = pathlib.Path(os.environ.get("ISGUI_SRC", str(ROOT / "src")))
sys.path[0:0] = [str(SRC), str(TK / "scripts")]

FIX = TK / "tests" / "fixtures" / "nai_shape"
NAI_BG = FIX / "nai_background.xml"
NAI_SM = FIX / "nai_sample_cs137.xml"

from PySide6 import QtWidgets
from gamma.calibration.intrinsic_shape import profiles
from gamma.calibration.intrinsic_shape.core import load_spectrum


# ─── Фикстуры ────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def app():
    """Оффскрин QApplication."""
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture
def msgs(monkeypatch):
    """Перехват QMessageBox, сбор текстов в список."""
    box: list[str] = []

    def fake(*a, **k):
        box.append(str(a[2]) if len(a) > 2 else "")
        return QtWidgets.QMessageBox.Ok

    def question(*a, **k):
        box.append(str(a[2]) if len(a) > 2 else "")
        return QtWidgets.QMessageBox.Yes

    monkeypatch.setattr(QtWidgets.QMessageBox, "warning", fake)
    monkeypatch.setattr(QtWidgets.QMessageBox, "critical", fake)
    monkeypatch.setattr(QtWidgets.QMessageBox, "information", fake)
    monkeypatch.setattr(QtWidgets.QMessageBox, "about", fake)
    monkeypatch.setattr(QtWidgets.QMessageBox, "question", question)
    return box


@pytest.fixture
def win(app, msgs, tmp_path, monkeypatch):
    """MainWindow с профилем nai_gamma1s."""
    monkeypatch.setenv("ISGUI_SETTINGS_INI", str(tmp_path / "settings.ini"))
    from isgui.main_window import MainWindow
    w = MainWindow()
    w.det_combo.setCurrentIndex(w.det_combo.findData("nai_gamma1s"))
    yield w
    w.close()
    w.deleteLater()
    app.processEvents()
    profiles.activate(None)


# ─── Вспомогательные функции ─────────────────────────────────────────────────

def load_pair(win, bg, sm):
    """Загрузить пару файлов (фон + образец) в окно."""
    win.bg_edit.setText(str(bg))
    win.sample_edit.setText(str(sm))
    win._pair_key = None
    win._update_pair()


def wait_job(app, win, timeout=120):
    """Ожидание завершения job."""
    t0 = time.monotonic()
    while True:
        app.processEvents()
        if win.job is None:
            break
        if win.job.isFinished():
            for _ in range(20):
                app.processEvents()
                time.sleep(0.01)
            break
        if time.monotonic() - t0 > timeout:
            raise AssertionError("операция не завершилась")
        time.sleep(0.02)


def run_window(app, win, out):
    """Запуск NaI-расчёта через окно в каталог out."""
    out.mkdir(exist_ok=True)
    win.nai.out_dir.setText(str(out))
    win._do_nai()
    wait_job(app, win)


# ─── Тесты ───────────────────────────────────────────────────────────────────

def test_zero_live_time_is_user_error(tmp_path):
    """Нулевое живое время → UserError с русским текстом."""
    text = NAI_BG.read_text(encoding="utf-8").replace(">20000<", ">0<")
    bad = tmp_path / "live0.xml"
    bad.write_text(text, encoding="utf-8")
    from isgui import backend, backend_nai
    with pytest.raises(backend.UserError) as ei:
        backend_nai.run_nai(str(bad), str(NAI_SM), "nai_gamma1s", out_dir=str(tmp_path))
    assert "живого времени" in str(ei.value)
    assert "live0.xml" in str(ei.value)


def test_bg_k_curve_is_background_times_k(app, win, msgs, tmp_path):
    """Кривая «Фон·k» равна фону, умноженному на k."""
    load_pair(win, NAI_BG, NAI_SM)
    win.nai.log_check.setChecked(False)
    run_window(app, win, tmp_path / "out")
    assert msgs == [], str(msgs)
    r = win.nai.result
    assert r is not None
    k = float(r["res"]["k"])
    assert k != 1.0
    items = {it.name(): it for it in win.nai.plot.getPlotItem().listDataItems()}
    y = np.asarray(items["Фон·k (масштаб k)"].yData, float)
    bg = np.asarray(load_spectrum(NAI_BG)["counts"], float)
    assert y.shape == bg.shape
    assert np.allclose(y, bg * k, rtol=1e-9, atol=0)
    assert not np.allclose(y, bg)


def test_background_kind_ignores_chosen_mode(app, win, msgs, tmp_path):
    """Режим фон-фон: kind=background, выбранное совмещение не влияет."""
    load_pair(win, NAI_BG, NAI_SM)
    win.nai.kind_combo.setCurrentIndex(win.nai.kind_combo.findData("background"))
    win.nai.force_combo.setCurrentIndex(win.nai.force_combo.findData("shape"))
    assert win.nai.params()["force"] == "shape"
    run_window(app, win, tmp_path / "bgbg")
    assert msgs == [], str(msgs)
    r = win.nai.result
    assert r is not None
    assert r["kind"] == "background"
    assert r["res"]["kind"] == "background"


def test_truncated_xml_message_is_russian(app, win, msgs, tmp_path):
    """Обрезанный XML → русское сообщение без английских XML-терминов."""
    raw = NAI_BG.read_text(encoding="utf-8")
    cut = tmp_path / "cut.xml"
    cut.write_text(raw[: len(raw) // 3], encoding="utf-8")
    load_pair(win, cut, NAI_SM)
    assert "обрезан или повреждён" in win.pair_label.text()
    run_window(app, win, tmp_path / "o4")
    assert msgs, "нет сообщения об ошибке"
    m = msgs[-1]
    assert "обрезан или повреждён" in m
    for word in ("line", "column", "token", "unclosed", "not well-formed", "no element"):
        assert word not in m.lower(), m
    assert "Traceback" not in m
