"""Фоновый поток, журнал, графики. Сгенерировано gen_code.py (qwen3.8:27b) по specs/widgets.md, 2026-10-10.
Job -- по образцу LoaderThread (поток загрузки с сигналами)."""
import time
import traceback
import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets, QtGui
from . import backend

NEG_NAME = "Разность < 0 (|net|)"
# правка Claude: доля обнулённых каналов разности в области 60-1700 (шкала 2048, пересчёт по числу каналов),
# выше которой вычитание фактически обнуляет образец (фон снят в другой обстановке); 30 % -- порог из запроса
ZEROED_WARN_FRAC = 0.30


def zeroed_warning(signed):
    """Текст жёлтого предупреждения или "" -- по доле каналов разности < 0 в области 60-1700."""
    sg = np.asarray(signed, float)
    n = len(sg)
    part = sg[int(60 * n / 2048):min(n, int(1700 * n / 2048) + 1)]
    frac = float(np.mean(part < 0)) if len(part) else 0.0
    if frac <= ZEROED_WARN_FRAC:
        return ""
    return ("Фон·k выше образца в %.0f %% каналов: внешний фон измерения фона выше, чем у образца (другая обстановка/"
            "защита) — вычитание обнуляет образец; смотрите разность со знаком" % (100 * frac))


class OpTimer(QtCore.QObject):
    """Правка Claude: счётчик «прошло mm:ss» длительной операции; обновляется таймером Qt в UI-потоке раз в секунду."""
    def __init__(self, label, parent=None):
        super().__init__(parent)
        self.label, self.title, self.t0, self.elapsed = label, "", None, 0.0
        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(1000)
        self.timer.timeout.connect(self.tick)

    def start(self, title):
        self.title, self.t0 = title, time.monotonic()
        self.tick()
        self.timer.start()

    def tick(self):
        self.elapsed = time.monotonic() - self.t0
        self.label.setText("%s: прошло %02d:%02d" % (self.title, int(self.elapsed) // 60, int(self.elapsed) % 60))

    def stop(self):
        self.timer.stop()
        self.tick()
        return self.elapsed


class Job(QtCore.QThread):
    """Фоновый исполнитель задач с передачей ошибок в UI-поток."""
    line = QtCore.Signal(str)
    done = QtCore.Signal(object)
    failed = QtCore.Signal(str, str)

    def __init__(self, fn, *args, parent=None, **kwargs):
        super().__init__(parent)
        self._fn = fn
        self._args = args
        self._kwargs = kwargs

    def _emit(self, s):
        self.line.emit(s)

    def run(self):
        try:
            res = self._fn(*self._args, log=self._emit, **self._kwargs)
            self.done.emit(res)
        except Exception as exc:
            self.failed.emit(backend.friendly(exc), traceback.format_exc())


class LogView(QtWidgets.QPlainTextEdit):
    """Только-для-чтения журнал вывода."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setReadOnly(True)
        self.setMaximumBlockCount(20000)
        self.setStyleSheet("font-family: Consolas, monospace; font-size: 9pt;")

    def append_text(self, s):
        self.moveCursor(QtGui.QTextCursor.End)
        self.insertPlainText(s)
        self.moveCursor(QtGui.QTextCursor.End)


class SpectrumPlot(QtWidgets.QWidget):
    """Отображение спектра с элементами управления."""

    def __init__(self, parent=None):
        super().__init__(parent)
        pg.setConfigOptions(antialias=True)
        self.plot = pg.PlotWidget(background="w")
        self.plot.showGrid(x=True, y=True, alpha=0.3)
        self.legend = self.plot.addLegend()

        controls = QtWidgets.QHBoxLayout()
        self.axis_combo = QtWidgets.QComboBox()
        self.axis_combo.addItems(["Каналы (шкала фона)", "Энергия, кэВ"])
        self.mode_combo = QtWidgets.QComboBox()
        self.mode_combo.addItems(["Разность: обнулённая", "Разность: со знаком"])
        self.unit_combo = QtWidgets.QComboBox()   # правка Claude (раунд 4, п.11): ось Y -- cps (по умолчанию) или отсчёты
        self.unit_combo.addItems(["cps на канал", "Отсчёты на канал"])
        self.log_check = QtWidgets.QCheckBox("Лог. шкала Y")
        self.log_check.setChecked(True)
        self.markers_check = QtWidgets.QCheckBox("Метки реперов")
        self.markers_check.setChecked(True)
        self.show_sample = QtWidgets.QCheckBox("Образец")
        self.show_sample.setChecked(True)
        self.show_bg = QtWidgets.QCheckBox("Фон·k")
        self.show_bg.setChecked(True)
        self.show_net = QtWidgets.QCheckBox("Разность")
        self.show_net.setChecked(True)

        controls.addWidget(self.axis_combo)
        controls.addWidget(self.mode_combo)
        controls.addWidget(self.unit_combo)
        controls.addWidget(self.log_check)
        controls.addWidget(self.markers_check)
        controls.addWidget(self.show_sample)
        controls.addWidget(self.show_bg)
        controls.addWidget(self.show_net)
        controls.addStretch()

        self.axis_combo.currentIndexChanged.connect(lambda *a: self.redraw())
        self.mode_combo.currentIndexChanged.connect(lambda *a: self.redraw())
        self.log_check.toggled.connect(lambda *a: self.redraw())
        self.unit_combo.currentIndexChanged.connect(lambda *a: self.redraw())
        self.markers_check.toggled.connect(lambda *a: self.redraw())
        self.show_sample.toggled.connect(lambda *a: self.redraw())
        self.show_bg.toggled.connect(lambda *a: self.redraw())
        self.show_net.toggled.connect(lambda *a: self.redraw())

        self.note = QtWidgets.QLabel("Нет данных: выполните «Вычесть»")
        self.note.setWordWrap(True)
        self.note.setStyleSheet("color: #666;")

        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(controls)
        layout.addWidget(self.plot, 1)
        layout.addWidget(self.note)

        self._data = None
        self._markers = []
        self._items = []

    def set_data(self, result, markers):
        self._data = result
        self._markers = markers
        self.redraw()

    def redraw(self):
        self.plot.clear()
        if self.legend:
            self.legend.clear()

        if self._data is None:
            return

        data = self._data
        energy = self.axis_combo.currentIndex() == 1
        x = data["energy"] if energy else data["channels"].astype(float)

        logy = self.log_check.isChecked()
        self.plot.setLogMode(x=False, y=logy)

        def prep(y):
            y = np.asarray(y, float)
            if logy:
                return np.where(y > 0, y, np.nan)
            return y

        net = data["signed"] if self.mode_combo.currentIndex() == 1 else data["zeroed"]
        cps = self.unit_combo.currentIndex() == 0 and data.get("live_s")
        if cps:   # образец и фон·k -- отсчёты за LT_s; разность в файле -- за LT_out
            data = dict(data, sample=np.asarray(data["sample"], float) / data["live_s"],
                        bg_k=np.asarray(data["bg_k"], float) / data["live_s"])
            net = np.asarray(net, float) / data["live_out"]

        if self.show_sample.isChecked():
            self.plot.plot(x, prep(data["sample"]), pen=pg.mkPen((30, 90, 200), width=1),
                           name="Образец (на шкале фона)", connect="finite")
        if self.show_bg.isChecked():
            self.plot.plot(x, prep(data["bg_k"]), pen=pg.mkPen((230, 120, 20), width=1),
                           name="Фон·k (совмещён)", connect="finite")
        if self.show_net.isChecked() and logy:   # правка Claude: лог-шкала -- разность со знаком, |net<0| отдельно
            sg = np.asarray(data["signed"], float) / (data["live_out"] if cps else 1.0)
            self.plot.plot(x, np.where(sg > 0, sg, np.nan), pen=pg.mkPen((20, 140, 40), width=1.5),
                           name="Разность", connect="finite")
            self.plot.plot(x, np.where(sg < 0, -sg, np.nan), connect="finite", name=NEG_NAME,
                           pen=pg.mkPen((200, 30, 30), width=1.2, style=QtCore.Qt.DashLine))
        elif self.show_net.isChecked():
            self.plot.plot(x, prep(net), pen=pg.mkPen((20, 140, 40), width=1.5),
                           name="Разность", connect="finite")

        if self.markers_check.isChecked() and self._markers:
            for m in self._markers:
                pos = m["E"] if energy else m["ch"]
                if pos is None or not np.isfinite(pos):
                    continue
                line = pg.InfiniteLine(pos=pos, angle=90, movable=False,
                                       pen=pg.mkPen((150, 0, 150), style=QtCore.Qt.DashLine),
                                       label=m["name"],
                                       labelOpts={"position": 0.92, "rotateAxis": (1, 0),
                                                   "color": (120, 0, 120),
                                                   "anchors": [(0, 0), (0, 0)]})
                self.plot.addItem(line)

        self.plot.setLabel("bottom", "Энергия, кэВ" if energy else "Канал (шкала фона)")
        self.plot.setLabel("left", "cps на канал (S/LT_s - B/LT_b)" if cps else "Отсчёты на канал")

        self.zeroed_warning = zeroed_warning(data["signed"]) if self.mode_combo.currentIndex() == 0 else ""
        self.note.setStyleSheet("background: #fff3b0; color: #6b4e00;" if self.zeroed_warning else "color: #666;")
        self.note.setText(self.zeroed_warning)

    def clear_data(self):
        self._data = None
        self._markers = []
        self.plot.clear()
        self.note.setText("Нет данных: выполните «Вычесть»")


class FwhmPlot(QtWidgets.QWidget):
    """Отображение модели ПШПВ и измеренных точек."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.plot = pg.PlotWidget(background="w")
        self.plot.showGrid(x=True, y=True, alpha=0.3)
        self.plot.addLegend()
        self.plot.setLabel("bottom", "Канал")
        self.plot.setLabel("left", "ПШПВ, каналов")
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self.plot)

    def set_data(self, points, x, y):
        self.plot.clear()
        self.plot.plot(x, y, pen=pg.mkPen((30, 90, 200), width=2), name="Модель")
        self.plot.plot([p[0] for p in points], [p[1] for p in points], pen=None,
                       symbol="o", symbolSize=9, symbolBrush=(220, 50, 50),
                       name="Измеренные точки")


def fill_table(table, headers, rows):
    """Заполняет QTableWidget заголовками и строками."""
    table.clear()
    table.setColumnCount(len(headers))
    table.setHorizontalHeaderLabels(headers)
    table.setRowCount(len(rows))
    for r, row in enumerate(rows):
        for c, v in enumerate(row):
            text = ("%.3f" % v) if isinstance(v, float) else str(v)
            item = QtWidgets.QTableWidgetItem(text)
            item.setFlags(item.flags() & ~QtCore.Qt.ItemIsEditable)
            table.setItem(r, c, item)
    table.resizeColumnsToContents()


def status_html(ok, text):
    """Возвращает HTML-строку статуса."""
    return '<span style="color:%s; font-weight:bold">%s %s</span>' % (
        "#1a7f37" if ok else "#c62828", "✔" if ok else "✘", text)
