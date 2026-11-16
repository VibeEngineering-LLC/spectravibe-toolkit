"""Вкладка метода NaI по линиям ЕРН (ЭКСПЕРИМЕНТАЛЬНО). Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-11-16;
правка: Qt-флаги из QtCore, QDesktopServices из QtGui, «Открыть папку» вне списка блокируемых кнопок."""
import pathlib
import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets
from . import backend_nai

EXPERIMENTAL_BANNER = "ЭКСПЕРИМЕНТАЛЬНО: шкала не проверяется. Метод NaI по линиям ЕРН — для проверки метода, не для отчётных измерений; в каждом выходном файле пометка «[экспериментально]» и флаг шкалы."


class NaiPanel(QtWidgets.QWidget):
    """Вкладка метода NaI по линиям ЕРН (ЭКСПЕРИМЕНТАЛЬНО): параметры запуска, итог, график образец / фон·k / остаток.
Расчёт запускает главное окно (фоновый поток, отмена, счётчик времени). Сгенерировано gen_code.py (qwen3.8:27b) по
спецификации, 2026-11-16."""

    run_requested = QtCore.Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.pid = None
        self.result = None
        self.buttons = []

        layout = QtWidgets.QVBoxLayout(self)

        self.banner = QtWidgets.QLabel(EXPERIMENTAL_BANNER)
        self.banner.setWordWrap(True)
        self.banner.setStyleSheet(
            "background:#fff3cd; color:#8a0000; font-weight:bold; padding:6px; border:1px solid #d39e00;"
        )
        layout.addWidget(self.banner)

        info_lbl = QtWidgets.QLabel(
            "Файлы — строки «Фон» и «Образец» вверху окна (BecqMoni XML). Профиль — строка «Детектор»."
        )
        info_lbl.setWordWrap(True)
        layout.addWidget(info_lbl)

        grid = QtWidgets.QGridLayout()

        grid.addWidget(QtWidgets.QLabel("Вид расчёта:"), 0, 0)
        self.kind_combo = QtWidgets.QComboBox()
        for val, label in backend_nai.KIND_CHOICES:
            self.kind_combo.addItem(label, val)
        grid.addWidget(self.kind_combo, 0, 1)

        grid.addWidget(QtWidgets.QLabel("Режим совмещения:"), 1, 0)
        self.force_combo = QtWidgets.QComboBox()
        for val, label in backend_nai.FORCE_CHOICES:
            self.force_combo.addItem(label, val)
        grid.addWidget(self.force_combo, 1, 1)

        grid.addWidget(QtWidgets.QLabel("Нуклиды образца:"), 2, 0)
        self.nuc_list = QtWidgets.QListWidget()
        self.nuc_list.setMaximumHeight(80)
        self.nuc_list.setToolTip("Все отмечены — все нуклиды профиля")
        grid.addWidget(self.nuc_list, 2, 1)

        grid.addWidget(QtWidgets.QLabel("Папка для результата:"), 3, 0)
        self.out_dir = QtWidgets.QLineEdit()
        btn_browse = QtWidgets.QPushButton("Обзор…")
        btn_browse.clicked.connect(self._browse_dir)
        grid.addWidget(self.out_dir, 3, 1)
        grid.addWidget(btn_browse, 3, 2)

        grid.addWidget(QtWidgets.QLabel("Основа имени файлов:"), 4, 0)
        self.stem = QtWidgets.QLineEdit()
        self.stem.setPlaceholderText("по умолчанию — имя образца + _минус_фон_NaI")
        grid.addWidget(self.stem, 4, 1)

        self.kind_combo.currentIndexChanged.connect(self._on_kind)

        layout.addLayout(grid)

        btn_row = QtWidgets.QHBoxLayout()
        self.run_btn = QtWidgets.QPushButton("Запустить метод NaI")
        self.run_btn.clicked.connect(self.run_requested.emit)
        self.buttons.append(self.run_btn)
        btn_row.addWidget(self.run_btn)

        self.open_btn = QtWidgets.QPushButton("Открыть папку результата")
        self.open_btn.clicked.connect(self.open_out_dir)
        self.open_btn.setEnabled(False)
        btn_row.addWidget(self.open_btn)

        self.log_check = QtWidgets.QCheckBox("Лог. шкала Y")
        self.log_check.setChecked(True)
        self.log_check.toggled.connect(self.redraw)
        btn_row.addWidget(self.log_check)

        btn_row.addStretch()
        layout.addLayout(btn_row)

        self.view = QtWidgets.QTextBrowser()
        self.view.setOpenExternalLinks(False)
        self.view.setMaximumHeight(170)   # итог -- компактно, основное место графику
        layout.addWidget(self.view)

        self.plot = pg.PlotWidget(background="w")
        self.plot.showGrid(x=True, y=True, alpha=0.3)
        self.plot.addLegend()
        self.plot.setLabel("left", "отсчёты на канал")
        self.plot.setLabel("bottom", "Энергия, кэВ")
        self.plot.setMinimumHeight(260)
        layout.addWidget(self.plot, 1)

    # ------------------------------------------------------------------ API

    def set_profile(self, pid, title):
        self.pid = pid
        t = title if "экспериментально" in title else title + " [экспериментально]"
        self.plot.setTitle("Профиль: " + t)
        try:
            names = backend_nai.nuclide_names(pid)
        except Exception:
            names = []
        self.nuc_list.clear()
        for n in names:
            item = QtWidgets.QListWidgetItem(n)
            item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
            item.setCheckState(QtCore.Qt.Checked)
            self.nuc_list.addItem(item)

    def _on_kind(self, *_):
        s = self.kind_combo.currentData() == "sample"
        self.force_combo.setEnabled(s)
        self.nuc_list.setEnabled(s)

    def _browse_dir(self):
        d = QtWidgets.QFileDialog.getExistingDirectory(self, "Папка для результата", self.out_dir.text())
        if d:
            self.out_dir.setText(d)

    def params(self):
        names = [
            self.nuc_list.item(i).text()
            for i in range(self.nuc_list.count())
            if self.nuc_list.item(i).checkState() == QtCore.Qt.Checked
        ]
        all_count = self.nuc_list.count()
        nuclides = None if (len(names) == all_count) else names
        return {
            "kind": self.kind_combo.currentData(),
            "force": self.force_combo.currentData(),
            "nuclides": nuclides,
            "out_dir": self.out_dir.text().strip(),
            "stem": self.stem.text().strip(),
        }

    def default_out(self, sample_path):
        if not self.out_dir.text().strip() and sample_path:
            self.out_dir.setText(str(pathlib.Path(sample_path).parent))

    def show_result(self, r):
        self.result = r
        self.view.setHtml(backend_nai.summary_html(r))
        self.open_btn.setEnabled(True)
        self.redraw()

    def clear(self):
        self.result = None
        self.view.setHtml("")
        self.open_btn.setEnabled(False)
        self.plot.clear()

    def open_out_dir(self):
        d = (self.result or {}).get("out_dir") or self.out_dir.text().strip()
        if d and pathlib.Path(d).is_dir():
            QtGui.QDesktopServices.openUrl(QtCore.QUrl.fromLocalFile(d))

    def _pos(self, y):
        y = np.asarray(y, dtype=float)
        if self.log_check.isChecked():
            return np.where(y > 0, y, np.nan)
        return y

    def redraw(self, *_):
        self.plot.clear()
        if not self.result:
            return
        p = self.result["plot"]
        logy = self.log_check.isChecked()
        self.plot.setLogMode(x=False, y=logy)

        curves = [
            ("Образец", p["E_sm"], p["sample"], pg.mkPen("#1f77b4", width=1)),
            ("Фон·k (масштаб k)", p["E_bg"], p["bg_k"], pg.mkPen("#7f7f7f", width=1)),
            ("Остаток (разность − модель линий)", p["E_bg"], p["resid"], pg.mkPen("#d62728", width=1)),
        ]

        for name, x, y, pen in curves:
            if x is None or y is None:
                continue
            x = np.asarray(x, dtype=float)
            if name.startswith("Остаток") and logy:
                y_plot = np.abs(y)
            else:
                y_plot = y
            y_plot = np.asarray(y_plot, dtype=float)
            if len(x) != len(y_plot):
                continue
            x_use = x
            y_use = self._pos(y_plot)
            self.plot.plot(x_use, y_use, pen=pen, name=name, connect="finite")
