"""Окно с четырьмя профилями детектора (offscreen): выбор из слоя профилей, скрытие вкладок чужого метода, сохранение
выбора, запуск метода NaI на фикстурах репозитория, предупреждение о несовпадении прибора, отмена. Запуск: python -m pytest
apps/intrinsic-shape-gui/tests/test_gui_detectors.py. Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-11-16."""

import os, sys, time, pathlib
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT = pathlib.Path(__file__).resolve().parents[1]  # apps/intrinsic-shape-gui
TK = ROOT.parents[1]  # корень репозитория
SRC = pathlib.Path(os.environ.get("ISGUI_SRC", str(ROOT / "src")))
sys.path[0:0] = [str(SRC), str(TK / "scripts")]

FIX = TK / "tests" / "fixtures"
NAI_BG = FIX / "nai_shape" / "nai_background.xml"
NAI_SM = FIX / "nai_shape" / "nai_sample_cs137.xml"
LA_BG = FIX / "intrinsic_shape" / "labr3_background.xml"
LA_SM = FIX / "intrinsic_shape" / "labr3_lu176_sample.xml"

from PySide6 import QtWidgets
from gamma.calibration.intrinsic_shape import profiles
from gamma.io.atomspectra_xml import read_atomspectra_xml


# ── фикстуры ──────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def app():
    """Общий QApplication (offscreen)."""
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture
def msgs(monkeypatch):
    """Перехват QMessageBox: собирает тексты в список."""
    box = []

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
def ini(tmp_path, monkeypatch):
    """Изолированный файл настроек."""
    p = tmp_path / "settings.ini"
    monkeypatch.setenv("ISGUI_SETTINGS_INI", str(p))
    return p


@pytest.fixture
def win(app, msgs, ini):
    """Окно MainWindow на каждый тест."""
    from isgui.main_window import MainWindow

    w = MainWindow()
    yield w
    w.close()
    w.deleteLater()
    app.processEvents()
    profiles.activate(None)


# ── вспомогательные функции ────────────────────────────────────────────────


def select(win, pid):
    """Выбрать профиль по id в комбобоксе."""
    i = win.det_combo.findData(pid)
    assert i >= 0, pid
    win.det_combo.setCurrentIndex(i)


def load_pair(win, bg, sm):
    """Записать пути фона/образца и обновить пару."""
    win.bg_edit.setText(str(bg))
    win.sample_edit.setText(str(sm))
    win._update_pair()


def wait_job(app, win, timeout=120):
    """Дождаться завершения фоновой задачи."""
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
            raise AssertionError("операция не завершилась за %d с" % timeout)
        time.sleep(0.02)


def combo_ids(win):
    """Список id из комбобокса."""
    return [win.det_combo.itemData(i) for i in range(win.det_combo.count())]


def tab_state(win, w):
    """Видимость и доступность вкладки."""
    i = win.tabs.indexOf(w)
    return (win.tabs.isTabVisible(i), win.tabs.isTabEnabled(i))


# ── тесты ──────────────────────────────────────────────────────────────────


def test_combo_filled_from_profile_layer(win):
    """Комбобокс заполнен из слоя профилей, подписи корректны."""
    ids = combo_ids(win)
    assert ids == profiles.ids(all_methods=True)
    assert len(ids) == 4
    assert set(ids) == {"labr_atomnano5pro", "labr_as2pro", "nai_gamma1s", "nai_gs2020"}

    labels = {win.det_combo.itemData(i): win.det_combo.itemText(i) for i in range(win.det_combo.count())}

    for pid in ("nai_gamma1s", "nai_gs2020"):
        assert "[экспериментально]" in labels[pid]
        assert "метод по линиям ЕРН" in labels[pid]

    for pid in ("labr_atomnano5pro", "labr_as2pro"):
        assert "[экспериментально]" not in labels[pid]
        assert "метод по форме La-138" in labels[pid]

    assert "UNVERIFIED" in labels["labr_as2pro"]
    assert "UNVERIFIED" not in labels["labr_atomnano5pro"]


def test_foreign_method_tabs_hidden(win):
    """Вкладки чужого метода скрыты и неактивны."""
    select(win, "nai_gamma1s")
    for w in win.la_tabs:
        assert tab_state(win, w) == (False, False)
    assert tab_state(win, win.nai) == (True, True)
    assert win.tabs.currentWidget() is win.nai

    select(win, "labr_as2pro")
    for w in win.la_tabs:
        assert tab_state(win, w) == (True, True)
    assert tab_state(win, win.nai) == (False, False)
    assert len(win.la_tabs) == 8

    txt = win.passport.text()
    assert "UNVERIFIED" in txt and "4 парах" in txt and "La_beta_edge" in txt


def test_no_path_to_foreign_command(win, msgs):
    """Запуск чужого метода без путей → предупреждение, задача не стартует."""
    select(win, "nai_gamma1s")
    load_pair(win, NAI_BG, NAI_SM)
    n0 = len(msgs)
    win._do_register()
    assert win.job is None
    assert any("недоступна" in m for m in msgs[n0:])

    select(win, "labr_atomnano5pro")
    win.auto_det.setChecked(False)
    n1 = len(msgs)
    win._do_nai()
    assert win.job is None
    assert any("недоступна" in m for m in msgs[n1:])


def test_selection_persisted(app, msgs, ini):
    """Выбор профиля и авто-режим сохраняются в settings.ini."""
    from isgui.main_window import MainWindow

    w1 = MainWindow()
    select(w1, "nai_gs2020")
    w1.auto_det.setChecked(False)
    w1.settings.sync()
    w1.close()

    w2 = MainWindow()
    assert w2.det_combo.currentData() == "nai_gs2020"
    assert w2.auto_det.isChecked() is False
    assert tab_state(w2, w2.nai) == (True, True)
    w2.close()
    profiles.activate(None)


def test_nai_run_writes_experimental_xml(app, win, msgs, tmp_path):
    """Запуск NaI: результат, XML-разметка, HTML, график, баннер, лог."""
    select(win, "nai_gamma1s")
    load_pair(win, NAI_BG, NAI_SM)
    out = tmp_path / "out"
    out.mkdir()
    win.nai.out_dir.setText(str(out))
    win._do_nai()
    assert win.job is not None
    wait_job(app, win)

    assert msgs == [], "ошибки: %s" % msgs

    r = win.nai.result
    assert r is not None
    for key in ("xml_zeroed", "xml_signed"):
        s = read_atomspectra_xml(r["paths"][key])
        assert "[экспериментально]" in s.sample_id
        assert "экспериментально, шкала образца не проверяется" in s.comments

    html_txt = win.nai.view.toPlainText()
    assert "ЭКСПЕРИМЕНТАЛЬНО" in html_txt and "Шкала проверена: нет" in html_txt
    assert win.nai.banner.text().startswith("ЭКСПЕРИМЕНТАЛЬНО")

    names = [it.name() for it in win.nai.plot.getPlotItem().listDataItems()]
    assert names == ["Образец", "Фон·k (масштаб k)", "Остаток (разность − модель линий)"]
    assert "NaI" in win.nai.plot.getPlotItem().titleLabel.text
    assert win.nai.open_btn.isEnabled()
    assert "[TIME] Метод NaI" in win.log.toPlainText()


def test_autodetect_switches_to_la_profile(win):
    """Авто-определение переключает на профиль La при паре La."""
    select(win, "nai_gamma1s")
    win.auto_det.setChecked(True)
    load_pair(win, LA_BG, LA_SM)
    assert win.det_combo.currentData() == "labr_atomnano5pro"
    assert tab_state(win, win.la_tabs[0]) == (True, True)
    assert win.det_warning == ""


def test_mismatch_is_warning_not_refusal(win, msgs):
    """Несовпадение профиля — предупреждение, а не отказ."""
    win.auto_det.setChecked(False)
    select(win, "nai_gamma1s")
    load_pair(win, LA_BG, LA_SM)
    assert win.det_combo.currentData() == "nai_gamma1s"
    assert "распознан как" in win.det_warning
    assert "labr_atomnano5pro" in win.det_warning
    assert "распознан как" in win.pair_label.text()
    assert msgs == []

    select(win, "labr_atomnano5pro")
    load_pair(win, NAI_BG, NAI_SM)
    assert "не распознан ни одним профилем La-138" in win.det_warning
    assert msgs == []


def test_nai_hint_without_match(win):
    """Авто-режим без совпадения подсказывает выбрать вручную."""
    win.auto_det.setChecked(True)
    select(win, "nai_gs2020")
    load_pair(win, NAI_BG, NAI_SM)
    assert win.det_combo.currentData() == "nai_gs2020"
    assert "выберите профиль детектора вручную" in win.pair_label.text()


def test_cancel_nai_run(app, win, msgs, tmp_path):
    """Отмена запуска NaI: задача прекращена, файлов нет."""
    select(win, "nai_gamma1s")
    load_pair(win, NAI_BG, NAI_SM)
    out = tmp_path / "cancel"
    out.mkdir()
    win.nai.out_dir.setText(str(out))
    win._do_nai()
    win._cancel()
    wait_job(app, win)

    assert any(m.startswith("Отменено") for m in msgs), str(msgs)
    assert list(out.iterdir()) == []
    assert win.nai.result is None
    assert win.nai.run_btn.isEnabled()
    assert not win.cancel_btn.isVisible()


def test_broken_xml_friendly_error(app, win, msgs, tmp_path):
    """Повреждённый XML → дружелюбная ошибка на русском, без трейсбека."""
    select(win, "nai_gamma1s")
    bad = tmp_path / "bad.xml"
    bad.write_text("<ResultDataFile><broken", encoding="utf-8")
    empty = tmp_path / "empty.xml"
    empty.write_text("", encoding="utf-8")

    for f in (bad, empty):
        msgs.clear()
        win.bg_edit.setText(str(f))
        win.sample_edit.setText(str(NAI_SM))
        win._pair_key = None
        win._update_pair()
        win.nai.out_dir.setText(str(tmp_path))
        win._do_nai()
        wait_job(app, win)
        assert msgs and "Traceback" not in msgs[-1]
        assert any(ch in msgs[-1] for ch in "абвгдеёжзийклмнопрстуфхцчшщыэюя")


def test_about_and_version(win):
    """Диапазон about и версия в заголовке."""
    from isgui import backend_nai

    a = backend_nai.about_html()
    assert backend_nai.version().startswith("0.2.0-beta.")
    assert backend_nai.version() in win.windowTitle()

    for pid in profiles.ids(all_methods=True):
        assert profiles.get(pid)["title"].split(" [")[0] in a

    assert "UNVERIFIED" in a and "экспериментально" in a and "intrinsic_shape_calibration.md" in a
