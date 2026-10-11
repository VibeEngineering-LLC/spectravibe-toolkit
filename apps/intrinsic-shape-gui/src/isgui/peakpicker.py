"""Вкладка «Калибровка по пробе»: выбор пиков пробы щелчком по графику, диалог энергии (линия из ядерных данных или
ручной ввод), таблица реперов, набор реперов в JSON. Сгенерировано gen_code.py (qwen3.8:27b) по specs/peakpicker.md,
2026-10-10 (ответ модели повторил спецификацию перед кодом -- взят блок кода, )."""
import math
import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets
from .widgets import fill_table


class ReperDialog(QtWidgets.QDialog):
    def __init__(self, fit, E_est, candidates, parent=None, current=None):
        super().__init__(parent)
        self.setWindowTitle("Репер пробы")
        self.result_value = None
        self.candidates = candidates

        layout = QtWidgets.QVBoxLayout(self)

        # Метка с информацией о пике
        info_text = (
            "Пик: канал %.2f ± %.3f, ПШПВ %.2f кан, значимость %.1f; "
            "оценка энергии по текущей шкале %.1f кэВ"
        ) % (fit["ch"], fit["ch_err"], fit["fwhm"], fit["signif"], E_est)
        info_label = QtWidgets.QLabel(info_text)
        info_label.setWordWrap(True)
        info_label.setTextFormat(QtCore.Qt.RichText)
        layout.addWidget(info_label)

        # Предупреждение, если есть
        warn_list = fit.get("warn", [])
        if warn_list:
            warn_text = "Предупреждение: " + "; ".join(warn_list) + " — репер можно добавить"
            warn_label = QtWidgets.QLabel(warn_text)
            warn_label.setWordWrap(True)
            warn_label.setTextFormat(QtCore.Qt.RichText)
            warn_label.setStyleSheet("color: #c62828;")
            layout.addWidget(warn_label)

        # Радио-кнопки и таблица кандидатов
        self.radio_lib = QtWidgets.QRadioButton("Линия из ядерных данных")
        self.radio_lib.setChecked(True)
        layout.addWidget(self.radio_lib)

        self.table = QtWidgets.QTableWidget()
        self.table.setColumnCount(5)   # примечание -- дублет, чужой нуклид, ряд виден
        self.table.setHorizontalHeaderLabels(["Нуклид", "E, кэВ", "Интенсивность, %", "E − оценка, кэВ", "Примечание"])
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.table.verticalHeader().setVisible(False)

        if candidates:
            self.table.setRowCount(len(candidates))
            for i, c in enumerate(candidates):
                self.table.setItem(i, 0, QtWidgets.QTableWidgetItem(c["nuclide"]))
                self.table.setItem(i, 1, QtWidgets.QTableWidgetItem("%.3f" % c["E"]))
                self.table.setItem(i, 2, QtWidgets.QTableWidgetItem("%.2f" % c["I"]))
                self.table.setItem(i, 3, QtWidgets.QTableWidgetItem("%+.1f" % c["dE_est"]))
                self.table.setItem(i, 4, QtWidgets.QTableWidgetItem(c.get("note", "")))
                if c.get("foreign"):   # чужие нуклиды -- серым, строка 0 молча не выбирается
                    for j in range(5):
                        self.table.item(i, j).setForeground(QtGui.QColor("#888888"))
            self.table.clearSelection()
            self.table.setCurrentCell(-1, -1)
            self.table.itemDoubleClicked.connect(lambda item: self._accept())
        else:
            self.table.setRowCount(0)

        layout.addWidget(self.table)

        # Ручной ввод
        self.radio_man = QtWidgets.QRadioButton("Ввести энергию вручную (помечается «ручная»)")
        layout.addWidget(self.radio_man)

        manual_layout = QtWidgets.QHBoxLayout()
        self.e_edit = QtWidgets.QLineEdit()
        self.e_edit.setPlaceholderText("кэВ")
        self.e_edit.textChanged.connect(lambda: self.radio_man.setChecked(True))
        manual_layout.addWidget(self.e_edit)

        self.nuc_edit = QtWidgets.QLineEdit()
        self.nuc_edit.setPlaceholderText("нуклид (необязательно)")
        manual_layout.addWidget(self.nuc_edit)
        layout.addLayout(manual_layout)

        # Обработка текущего репера (правка)
        if current is not None:
            if current.get("source") == "ручная":
                self.radio_man.setChecked(True)
                self.e_edit.setText("%.3f" % current["E"])
                self.nuc_edit.setText(current.get("nuclide", ""))
            else:
                # Найти строку с тем же нуклидом и близкой энергией
                for i, c in enumerate(candidates):
                    if c["nuclide"] == current.get("nuclide") and abs(c["E"] - current["E"]) < 0.01:
                        self.table.selectRow(i)
                        break

        # Кнопки
        buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _accept(self):
        if self.radio_lib.isChecked():
            row = self.table.currentRow()
            if row < 0 or row >= len(self.candidates):
                QtWidgets.QMessageBox.warning(self, "Репер", "Выберите линию или введите энергию")
                return
            c = self.candidates[row]
            self.result_value = (float(c["E"]), c["nuclide"])
        else:
            try:
                E = float(self.e_edit.text().replace(",", "."))
                if E <= 0:
                    raise ValueError
            except (ValueError, TypeError):
                QtWidgets.QMessageBox.warning(self, "Репер", "Энергия — число больше нуля, кэВ")
                return
            nuc = self.nuc_edit.text().strip() or None
            self.result_value = (E, nuc)
        self.accept()


class PeakPickerPanel(QtWidgets.QWidget):
    scale_requested = QtCore.Signal()
    message = QtCore.Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.counts = None
        self.cal = None
        self.n = 0
        self.path = ""
        self.repers = []
        self.scale = None

        layout = QtWidgets.QVBoxLayout(self)

        # Метка с инструкцией
        instr_label = QtWidgets.QLabel(
            "Щёлкните по пику на графике пробы: центроид подгоняется в окне ±2 ПШПВ, "
            "затем выберите линию из ядерных данных или введите энергию вручную. "
            "К реперам пробы добавляются видимые в пробе собственные структуры LaBr3 (37, 789, 1472 кэВ). "
            "Двойной щелчок по строке таблицы — изменить энергию."
        )
        instr_label.setWordWrap(True)
        layout.addWidget(instr_label)

        # График
        self.plot = pg.PlotWidget(background="w")
        self.plot.showGrid(x=True, y=True)
        self.plot.setLabel("left", "Отсчёты")
        self.plot.setLabel("bottom", "Канал пробы")
        self.plot.getPlotItem().setLogMode(y=True)
        self.plot.setMinimumHeight(260)
        self.curve = self.plot.plot([], [], pen=pg.mkPen("#1565c0", width=1))
        self.marks = []
        self.plot.scene().sigMouseClicked.connect(self._on_click)
        layout.addWidget(self.plot)

        # Таблица реперов
        self.table = QtWidgets.QTableWidget()
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)   # иначе selectedRows() пуст
        self.table.itemDoubleClicked.connect(self._on_table_double_click)
        layout.addWidget(self.table)

        # Кнопки
        btn_layout = QtWidgets.QHBoxLayout()
        self.buttons = []

        btn_delete = QtWidgets.QPushButton("Удалить репер")
        btn_delete.clicked.connect(self._delete)
        btn_layout.addWidget(btn_delete)
        self.buttons.append(btn_delete)

        btn_clear = QtWidgets.QPushButton("Очистить")
        btn_clear.clicked.connect(self._clear)
        btn_layout.addWidget(btn_clear)
        self.buttons.append(btn_clear)

        btn_save = QtWidgets.QPushButton("Сохранить набор…")
        btn_save.clicked.connect(self._save)
        btn_layout.addWidget(btn_save)
        self.buttons.append(btn_save)

        btn_load = QtWidgets.QPushButton("Загрузить набор…")
        btn_load.clicked.connect(self._load)
        btn_layout.addWidget(btn_load)
        self.buttons.append(btn_load)

        btn_scale = QtWidgets.QPushButton("Шкала по пробе")
        btn_scale.clicked.connect(lambda *a: self.scale_requested.emit())   # clicked передаёт bool
        btn_layout.addWidget(btn_scale)
        self.buttons.append(btn_scale)

        layout.addLayout(btn_layout)

        # Сводка
        self.summary = QtWidgets.QLabel()
        self.summary.setWordWrap(True)
        self.summary.setTextFormat(QtCore.Qt.RichText)
        layout.addWidget(self.summary)

    def _on_table_double_click(self, item):
        row = item.row()
        self._edit_row(row)

    def set_sample(self, path):
        from gamma.calibration.intrinsic_shape.core import load_spectrum
        d = load_spectrum(path)
        self.counts = np.asarray(d["counts"], float)
        self.cal = list(d["cal"])
        self.n = len(self.counts)
        self.path = str(path)
        self.scale = None
        self.curve.setData(np.arange(self.n, dtype=float), np.maximum(self.counts, 0.5))
        self._refresh()

    def _coeffs(self):
        """шкала для диалога -- принятая шкала, иначе текущая шкала по реперам пробы (>= 2),
        иначе калибровка файла."""
        from gamma.calibration.intrinsic_shape import samplecal_q as q
        if self.scale is not None:
            return self.scale
        return q.current_scale(self.repers) or self.cal

    def _fwhm_keV(self, fit):
        """ПШПВ пика в кэВ по текущей шкале -- ширина слияния дублетов в списке линий."""
        from gamma.calibration.intrinsic_shape.core import poly_deriv
        coeffs = self._coeffs()
        return float(fit["fwhm"]) * abs(float(poly_deriv(coeffs, fit["ch"]))) if np.isfinite(fit.get("fwhm", float("nan"))) else None

    def _energy(self, x):
        from gamma.calibration.intrinsic_shape.core import poly_eval
        coeffs = self._coeffs()
        return float(poly_eval(coeffs, x))

    def _on_click(self, ev):
        if self.counts is None or ev.button() != QtCore.Qt.LeftButton:
            return

        vb = self.plot.getPlotItem().getViewBox()
        if not vb.sceneBoundingRect().contains(ev.scenePos()):
            return

        x = float(vb.mapSceneToView(ev.scenePos()).x())
        if not (0 <= x < self.n):
            return

        from gamma.calibration.intrinsic_shape import samplecal as sc
        fw = sc.fwhm_guess(x, self.n)
        fit = sc.fit_click(self.counts, x, fw, fwhm_expected=fw)

        if fit is None:
            self.message.emit("Пик у канала %.1f не подогнан (центроид не определён) — репер не добавлен\n" % x)
            return

        if any("мультиплет" in w for w in fit["warn"]):   # предложение двойного фита
            ans = QtWidgets.QMessageBox.question(self, "Вероятно мультиплет", "Пик шире ожидаемого -- вероятно мультиплет.\n"
                                                 "Выполнить двойной фит и взять компонент, ближайший к щелчку?")
            if ans == QtWidgets.QMessageBox.Yes:
                dd = sc.q.fit_doublet(self.counts, x, fw)
                if dd:
                    fit = dict(min(dd, key=lambda c: abs(c["ch"] - x)), warn=["двойной фит: компоненты %.1f и %.1f кан" % (dd[0]["ch"], dd[1]["ch"])])
                else:
                    self.message.emit("Двойной фит не сошёлся -- оставлен одиночный\n")
        E_est = self._energy(fit["ch"])
        fwk = self._fwhm_keV(fit)
        cands = sc.candidate_lines(E_est, fwhm_keV=fwk, visible_series=sc.q.series_seen(self.repers))
        dlg = ReperDialog(fit, E_est, cands, self)
        if dlg.exec() == QtWidgets.QDialog.Accepted and dlg.result_value:
            E, nuc = dlg.result_value
            r = sc.reper_from_choice(fit, E, nuc, fwhm_keV=fwk)
            dup = sc.duplicate_of(r, self.repers)   # одна линия -- один репер
            if dup is not None:
                self.message.emit("Репер %s не добавлен: повтор репера %s (та же линия -- в пределах 0.5 ПШПВ по энергии или каналу)\n"
                                  % (r["name"], dup["name"]))
                return
            self.repers.append(r)
            warn_str = ("; " + "; ".join(r["warn"])) if r.get("warn") else ""
            self.message.emit(
                "Репер %s: канал %.2f ± %.3f, E = %.3f кэВ (%s)%s\n"
                % (r["name"], r["ch"], r["ch_err"], r["E"], r["source"], warn_str)
            )
            self._refresh()

    def _edit_row(self, row):
        if row < 0 or row >= len(self.repers):
            return

        from gamma.calibration.intrinsic_shape import samplecal as sc
        r = self.repers[row]
        fit = {
            "ch": r["ch"],
            "ch_err": r["ch_err"],
            "fwhm": r.get("fwhm", float("nan")),
            "signif": r.get("signif", float("nan")),
            "warn": r.get("warn", []),
        }
        E_est = self._energy(r["ch"])
        fwk = self._fwhm_keV(fit)
        cands = sc.candidate_lines(E_est, fwhm_keV=fwk, visible_series=sc.q.series_seen(self.repers))
        dlg = ReperDialog(fit, E_est, cands, self, current=r)
        if dlg.exec() == QtWidgets.QDialog.Accepted and dlg.result_value:
            E, nuc = dlg.result_value
            new_r = sc.reper_from_choice(fit, E, nuc, fwhm_keV=fwk)
            dup = sc.duplicate_of(new_r, [q for i, q in enumerate(self.repers) if i != row])
            if dup is not None:
                self.message.emit("Репер не изменён: энергия %.1f кэВ -- повтор репера %s\n" % (new_r["E"], dup["name"]))
                return
            self.repers[row] = new_r
            self._refresh()

    def _delete(self):
        selected_rows = self.table.selectionModel().selectedRows()
        indices = sorted([idx.row() for idx in selected_rows], reverse=True)
        for i in indices:
            if 0 <= i < len(self.repers):
                del self.repers[i]
        self._refresh()

    def _clear(self):
        self.repers = []
        self._refresh()

    def _save(self):
        if not self.repers:
            self.message.emit("Нет реперов для сохранения\n")
            return

        fn, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Сохранить набор реперов", "", "Набор реперов (*.json)"
        )
        if not fn:
            return

        from gamma.calibration.intrinsic_shape import samplecal as sc
        anc = sc.anchors_for_save(self.counts, self.cal, self.repers, verbose=False)   # ложные опоры не сохраняются
        sc.save_json(fn, self.repers, self.path, self.n, anc, cal=self.cal)
        self.message.emit("Набор реперов сохранён: %s\n" % fn)

    def _load(self):
        if self.counts is None:
            self.message.emit("Сначала выберите файл образца\n")
            return

        fn, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Загрузить набор реперов", "", "Набор реперов (*.json)"
        )
        if not fn:
            return

        from gamma.calibration.intrinsic_shape import samplecal as sc
        try:
            data = sc.load_json(fn)
            if int(data["n_ch"]) != self.n:
                raise ValueError("набор для %d каналов, спектр — %d" % (data["n_ch"], self.n))
            reps, m, dropped = sc.apply_to(data, self.counts, self.cal, verbose=False)
            for msg in dropped:
                self.message.emit("ПРЕДУПРЕЖДЕНИЕ: " + msg + "\n")
            self.repers = [r for r in reps if r["source"] != sc.SRC_LABR]
            self.message.emit(
                "Набор реперов применён: отображение опор x' = %.3f + %.5f·x, реперов %d\n"
                % (m[0], m[1], len(self.repers))
            )
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, "Набор реперов", str(exc))
        self._refresh()

    def _refresh(self):
        # Удалить старые метки
        for mark in self.marks:
            self.plot.removeItem(mark)
        self.marks = []

        # Добавить новые метки
        for r in self.repers:
            color = "#c62828" if r["source"] == "ручная" else "#2e7d32"
            line = pg.InfiniteLine(
                pos=r["ch"],
                angle=90,
                pen=pg.mkPen(color, style=QtCore.Qt.DashLine),
                label="%.1f" % r["E"],
                labelOpts={"position": 0.9, "color": "#333"},
            )
            self.plot.addItem(line)
            self.marks.append(line)

        # Заполнить таблицу
        headers = ["Репер", "Канал пробы", "± канал", "ПШПВ, кан", "Значимость", "E, кэВ", "Источник", "Предупреждение"]
        rows = [
            (
                r["name"],
                float(r["ch"]),
                float(r["ch_err"]),
                float(r.get("fwhm", float("nan"))),
                float(r.get("signif", float("nan"))),
                float(r["E"]),
                r["source"],
                "; ".join(r.get("warn", [])),
            )
            for r in self.repers
        ]
        fill_table(self.table, headers, rows)

        # Сводка
        n_manual = sum(1 for r in self.repers if r["source"] == "ручная")
        self.summary.setText("Реперов пробы: %d (ручных %d)" % (len(self.repers), n_manual))

    def set_scale(self, coeffs):
        self.scale = list(coeffs) if coeffs else None

    def set_enabled(self, on):
        for btn in self.buttons:
            btn.setEnabled(on)
