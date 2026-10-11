"""Главное окно. Сгенерировано gen_code.py (qwen3.8:27b) по specs/main_window.md, 2026-10-10.
Расчёты -- только через isgui.backend в фоновом потоке (isgui.widgets.Job)."""
import json
import pathlib
import html

import numpy as np

from PySide6 import QtCore, QtGui, QtWidgets

from . import backend, backend_pair, backend_sample   # обработчики пары (автопресет, линии образца, ветвь шкалы)
from .panels import LowScalePlot, AnchorBox
from .peakpicker import PeakPickerPanel   # калибровка по пробе
from .widgets import Job, LogView, SpectrumPlot, FwhmPlot, fill_table, status_html, OpTimer


class MainWindow(QtWidgets.QMainWindow):
    """Главное окно программы калибровки и вычитания фона."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Калибровка по форме собственных структур (LaBr3) — вычитание фона для BecqMoni — тестовая версия")
        self.resize(1250, 900)

        self.settings = QtCore.QSettings("SpectraVibe", "IntrinsicShapeGUI")

        # Состояние
        self.h = None
        self.h_json = None
        self.ecal = None
        self.ecal_points = []
        self.net = None
        self.job = None
        self.ecal_res = None   # полный результат «Шкалы» (составная шкала, точки ПШПВ пары)
        self.manual = None
        self.n_ch = 2048
        self.saved = False
        self.sample_res = None   # результат режима «по пробе» (разность на сетке пробы)
        self.buttons = []

        # Центральное виджет
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        main_layout = QtWidgets.QVBoxLayout(central)

        # 1. Файлы
        files_group = self._build_files()
        main_layout.addWidget(files_group)

        # 2. Сплиттер: вкладки сверху, журнал снизу
        splitter = QtWidgets.QSplitter(QtCore.Qt.Vertical)

        self.tabs = QtWidgets.QTabWidget()
        self.tabs.addTab(self._build_register(), "1. Совместить")
        self.tabs.addTab(self._build_ecal(), "2. Шкала")
        self.tabs.addTab(self._build_samplecal(), "2б. Шкала по пробе")
        self.tabs.addTab(self._build_subtract(), "3. Вычесть")
        self.tabs.addTab(self._build_plot(), "График")
        self.tabs.addTab(self._build_fwhm(), "ПШПВ")
        self.tabs.addTab(self._build_backgrounds(), "Фоны")
        self.tabs.addTab(self._build_selftest(), "Самопроверка")
        splitter.addWidget(self.tabs)

        log_widget = QtWidgets.QWidget()
        log_layout = QtWidgets.QVBoxLayout(log_widget)
        log_layout.addWidget(QtWidgets.QLabel("Журнал работы"))
        self.log = LogView()
        log_layout.addWidget(self.log)
        splitter.addWidget(log_widget)

        splitter.setSizes([620, 220])
        main_layout.addWidget(splitter)

        # 3. Статус-бар
        self.progress = QtWidgets.QProgressBar()
        self.progress.setMaximumWidth(220)
        self.progress.setVisible(False)
        self.statusBar().addPermanentWidget(self.progress)
        self.elapsed_label = QtWidgets.QLabel("")   # счётчик времени и отмена длительной операции
        self.statusBar().addPermanentWidget(self.elapsed_label)
        self.op_timer = OpTimer(self.elapsed_label, self)
        self.cancel_btn = QtWidgets.QPushButton("Отмена")
        self.cancel_btn.setVisible(False)
        self.cancel_btn.clicked.connect(self._cancel)
        self.statusBar().addPermanentWidget(self.cancel_btn)

        self.statusBar().showMessage("Выберите файл фона и файл образца")

    # --- Вспомогательные методы ---

    def _btn(self, text, slot):
        """Создать кнопку, добавить в список управляемых кнопок."""
        b = QtWidgets.QPushButton(text)
        b.clicked.connect(lambda *a: slot())
        self.buttons.append(b)
        return b

    def _html_view(self):
        """Создать QTextBrowser для отображения HTML."""
        v = QtWidgets.QTextBrowser()
        v.setOpenExternalLinks(False)
        return v

    def _last_dir(self):
        """Получить последнюю рабочую папку из настроек."""
        return str(self.settings.value("last_dir", str(pathlib.Path.home())))

    def _remember_dir(self, path):
        """Запомнить папку файла в настройках."""
        self.settings.setValue("last_dir", str(pathlib.Path(path).parent))

    def _pick_file(self, title):
        """Выбрать один файл."""
        fn, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, title, self._last_dir(), "Спектры BecqMoni (*.xml);;Все файлы (*)"
        )
        if fn:
            self._remember_dir(fn)
            return fn
        return ""

    def _pick_files(self, title):
        """Выбрать несколько файлов."""
        fns, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self, title, self._last_dir(), "Спектры BecqMoni (*.xml);;Все файлы (*)"
        )
        if fns:
            self._remember_dir(fns[0])
            return fns
        return []

    def _pick_dir(self, title):
        """Выбрать папку."""
        d = QtWidgets.QFileDialog.getExistingDirectory(self, title, self._last_dir())
        if d:
            self.settings.setValue("last_dir", d)
            return d
        return ""

    def _error(self, msg):
        """Показать сообщение об ошибке."""
        QtWidgets.QMessageBox.warning(self, "Ошибка", msg)
        self.statusBar().showMessage("Ошибка: " + msg.splitlines()[0])

    def _paths(self):
        """Получить пути к файлам фона и образца."""
        return (self.bg_edit.text().strip(), self.sample_edit.text().strip())

    def _run(self, title, fn, args, kwargs, on_done):
        """Запустить фоновую задачу."""
        if self.job is not None and self.job.isRunning():
            self._error("Дождитесь окончания текущей операции")
            return

        self.log.append_text("\n=== %s ===\n" % title)
        for b in self.buttons:
            b.setEnabled(False)
        self.apply_scale_btn.setEnabled(False)   # во время любой задачи -- недоступна
        self.progress.setRange(0, 0)
        self.progress.setVisible(True)
        self.statusBar().showMessage(title + " — выполняется…")
        from gamma.calibration.intrinsic_shape import robust as _rb
        _rb.CANCEL["flag"] = False
        self.op_timer.start(title)
        self.cancel_btn.setVisible(True)

        job = Job(fn, *args, **kwargs)
        job.line.connect(self.log.append_text)
        job.done.connect(lambda res: self._finish(title, res, on_done))
        job.failed.connect(lambda msg, tb: self._fail(title, msg, tb))
        self.job = job
        job.start()

    def _end_job(self):
        """Завершить состояние выполнения задачи."""
        for b in self.buttons:
            b.setEnabled(True)
        self.progress.setVisible(False)
        self.cancel_btn.setVisible(False)
        if self.op_timer.t0 is not None and self.op_timer.timer.isActive():
            self.log.append_text("[TIME] %s: %.0f с\n" % (self.op_timer.title, self.op_timer.stop()))

    def _cancel(self):
        """запросить отмену; рабочий поток проверяет флаг между итерациями."""
        from gamma.calibration.intrinsic_shape import robust as _rb
        _rb.CANCEL["flag"] = True
        self.log.append_text("Отмена запрошена — операция остановится после текущего шага\n")

    def _finish(self, title, res, on_done):
        """Обработка успешного завершения задачи."""
        self._end_job()
        try:
            on_done(res)
        except Exception as exc:
            self._error(backend.friendly(exc))
            return
        self.statusBar().showMessage(title + " — готово")

    def _fail(self, title, msg, tb):
        """Обработка ошибки задачи."""
        self._end_job()
        if msg.startswith("Внутренняя ошибка"):   # трассировка -- только для внутренних ошибок,
            self.log.append_text(tb)               # сообщение показывается ВСЕГДА (было внутри if)
        self.log.append_text("ОШИБКА: " + msg + "\n")
        self._error(msg)

    # --- Вкладка: Файлы ---

    def _build_files(self):
        """Построить группу выбора файлов."""
        group = QtWidgets.QGroupBox("Файлы спектров (BecqMoni XML)")
        layout = QtWidgets.QGridLayout(group)

        # Строка 0: Фон
        layout.addWidget(QtWidgets.QLabel("Фон:"), 0, 0)
        self.bg_edit = QtWidgets.QLineEdit()
        self.bg_edit.setPlaceholderText("файл фона .xml")
        self.bg_edit.editingFinished.connect(self._update_pair)
        layout.addWidget(self.bg_edit, 0, 1)
        btn_bg = QtWidgets.QPushButton("Обзор…")
        btn_bg.clicked.connect(self._choose_bg)
        layout.addWidget(btn_bg, 0, 2)

        # Строка 1: Образец
        layout.addWidget(QtWidgets.QLabel("Образец:"), 1, 0)
        self.sample_edit = QtWidgets.QLineEdit()
        self.sample_edit.setPlaceholderText("файл образца .xml")
        self.sample_edit.editingFinished.connect(self._update_pair)
        layout.addWidget(self.sample_edit, 1, 1)
        btn_sm = QtWidgets.QPushButton("Обзор…")
        btn_sm.clicked.connect(self._choose_sample)
        layout.addWidget(btn_sm, 1, 2)
        btn_new = QtWidgets.QPushButton("Новый спектр")   # сброс файлов и результатов без перезапуска
        btn_new.setToolTip("Сбросить выбранные файлы и все результаты, начать новую пару")
        btn_new.clicked.connect(self._new_spectrum)
        layout.addWidget(btn_new, 0, 3, 2, 1)
        self.buttons.append(btn_new)

        # Строка 2: Информация о паре
        self.pair_label = QtWidgets.QLabel("Файлы не выбраны")
        self.pair_label.setWordWrap(True)
        self.pair_label.setTextFormat(QtCore.Qt.RichText)
        layout.addWidget(self.pair_label, 2, 0, 1, 3)

        return group

    def _choose_bg(self):
        """Выбрать файл фона."""
        fn = self._pick_file("Файл фона")
        if fn:
            self.bg_edit.setText(fn)
            self._update_pair()

    def _choose_sample(self):
        """Выбрать файл образца."""
        fn = self._pick_file("Файл образца")
        if fn:
            self.sample_edit.setText(fn)
            self._update_pair()

    def _update_pair(self):
        """Обновить информацию о выбранной паре файлов."""
        # editingFinished срабатывает и при потере фокуса -- сбрасывать результаты только при смене файлов
        if self._paths() == getattr(self, "_pair_key", None):
            return
        self._pair_key = self._paths()
        # Сброс состояния
        self.h = None
        self.h_json = None
        self.ecal = None
        self.net = None
        self.ecal_points = []
        self.ecal_res = None; self.saved = False; self.sample_res = None
        self._reset_sample_scale()
        self.plot.clear_data()

        bg, sm = self._paths()
        if not bg or not sm:
            self.pair_label.setText("Выберите оба файла")
            return

        try:
            p = backend.pair_info(bg, sm)
        except Exception as exc:
            self.pair_label.setText(status_html(False, html.escape(backend.friendly(exc))))
            return

        self.n_ch = p["bg"]["n_ch"]   # число каналов пары (ПШПВ, параметры прибора)
        # Формирование HTML-информации
        parts = []
        bg_info = p["bg"]
        sm_info = p["sample"]
        parts.append("Фон: <b>%s</b> — живое время %.1f с, каналов %d, отсчётов %.0f" % (
            html.escape(bg_info["name"]), bg_info["live"], bg_info["n_ch"], bg_info["total"]
        ))
        parts.append("Образец: <b>%s</b> — живое время %.1f с, каналов %d, отсчётов %.0f" % (
            html.escape(sm_info["name"]), sm_info["live"], sm_info["n_ch"], sm_info["total"]
        ))
        parts.append("Отношение живых времён k = %.6f" % p["k"])
        if p["warning"]:
            parts.append(status_html(False, html.escape(p["warning"])))

        self.pair_label.setText("<br>".join(parts))
        try:   # спектр пробы -- на вкладку выбора пиков
            self.picker.set_sample(sm)
        except Exception as exc:
            self.log.append_text("Вкладка «Шкала по пробе»: спектр не загружен — %s\n" % backend.friendly(exc))

        # Предзаполнение имени и папки результата
        if not self.out_name.text() or self.out_name.text().endswith("_минус_фон"):
            self.out_name.setText(pathlib.Path(sm).stem + "_минус_фон")
        if not self.out_dir.text():
            self.out_dir.setText(str(pathlib.Path(sm).parent))

    def _new_spectrum(self, confirm=True):
        """«Новый спектр» -- сброс файлов, результатов и вкладок без перезапуска."""
        if confirm and (self.h is not None or self.ecal is not None) and not self.saved:
            ans = QtWidgets.QMessageBox.question(self, "Новый спектр", "Результат совмещения/шкалы не сохранён в файл. Сбросить всё?")
            if ans != QtWidgets.QMessageBox.Yes:
                return False
        self.bg_edit.clear(); self.sample_edit.clear(); self.out_name.clear(); self.out_dir.clear()
        self._pair_key = None
        self.h = self.h_json = self.ecal = self.ecal_res = self.net = self.manual = None
        self.ecal_points = []; self.n_ch = 2048; self.saved = False; self.sample_res = None
        self._reset_views()
        self._reset_sample_scale()
        self.log.append_text("\n=== Новый спектр: файлы и результаты сброшены ===\n")
        self.statusBar().showMessage("Выберите файл фона и файл образца")
        return True

    def _reset_views(self):
        """все вкладки -- в исходное состояние."""
        self.pair_label.setText("Файлы не выбраны")
        self.anchor_box.reset(); self.low_plot.clear_data(); self.plot.clear_data(); self.lt_combo.setCurrentIndex(0)
        self.picker._clear(); self.picker.curve.setData([], []); self.picker.counts = None
        for v in (self.reg_view, self.sub_view, self.dec_view, self.st_view, self.sc_view):
            v.setHtml("")
        for tb in (self.ecal_table, self.fw_table, self.inv_table):
            tb.clear(); tb.setRowCount(0); tb.setColumnCount(0)
        for lb in (self.ecal_summary, self.fw_view, self.inv_view):
            lb.setText("")
        self.fw_points_edit.setText(" ".join("%.2f:%.2f" % p for p in backend.default_fwhm_points()))
        self.fw_plot.plot.clear(); self._device_text = ""
        self.tabs.setCurrentIndex(0)

    # --- Вкладка: Совместить ---

    def _build_register(self):
        """Построить вкладку совмещения."""
        w = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(w)

        label = QtWidgets.QLabel(
            "Совмещение шкалы образца со шкалой фона по форме собственных структур кристалла. "
            "Занимает несколько секунд."
        )
        label.setWordWrap(True)
        layout.addWidget(label)
        self.anchor_box = AnchorBox()   # опоры автопресета, ручная поправка, линии образца
        layout.addWidget(self.anchor_box)

        layout.addWidget(self._btn("Совместить", self._do_register))

        self.reg_view = self._html_view()
        layout.addWidget(self.reg_view)

        return w

    def _do_register(self):
        """Запустить совмещение."""
        bg, sm = self._paths()
        try:   # ручная поправка опор и режим линий образца
            self.manual = backend_pair.parse_manual(self.anchor_box.manual_text())
        except backend.UserError as exc:
            self._error(str(exc))
            return
        self._run("Совмещение", backend_pair.run_register, (bg, sm),
                  {"manual": self.manual, "robust": self.anchor_box.robust_mode()}, self._show_register)

    def _show_register(self, r):
        """Показать результат совмещения."""
        self.h = r["h"]
        self.h_json = r["h_json"]
        self.ecal = None
        self.sample_res = None   # «Совместить» возвращает режим пары
        self.net = None

        h = r["h"]
        perr = r["perr"]
        x_sum = r["auto"]["anchors"]["sum"]["x"]   # канал сумм-пика этого прибора
        h846 = h[0] + h[1] * x_sum + h[2] * x_sum ** 2

        parts = []
        parts.append("<h3>Результат совмещения</h3>")
        parts.append(
            "h(x) = c0 + c1·x + c2·x²: c0 = %.5f ± %.5f, c1 = %.7f ± %.7f, c2 = %.4e ± %.1e" % (
                h[0], perr[0], h[1], perr[1], h[2], perr[2]
            )
        )
        parts.append("Сдвиг на канале %.0f (сумм-пик La-138): %+.2f канала" % (x_sum, h846 - x_sum))
        parts.append("Качество совпадения форм χ²/ndf = %.3f (ожидается около 1)" % r["chi2ndf"])
        parts.append(
            "Различных минимумов найдено: %d; из них почти равноценных (Δχ² ≤ 4): %d" % (
                r["n_minima"], r["n_close"]
            )
        )
        z_ok = abs(r["anchor_z"]) <= 2
        parts.append(
            "Якорь ~37 кэВ: отклонение z = %+.2f %s" % (
                r["anchor_z"], status_html(z_ok, "в норме" if z_ok else "тянет решение")
            )
        )

        parts.append("<h3>Требует толкования</h3>")
        if not r["interp_explained"]:
            parts.append(status_html(True, "Замечаний нет"))
        else:
            for item, expl in r["interp_explained"]:
                parts.append("<p><b>%s</b><br>%s</p>" % (html.escape(item), html.escape(expl)))

        parts.append("<p>Дальше: вкладка «2. Шкала».</p>")
        self.reg_view.setHtml("<br>".join(parts) + backend_pair.register_html(r))
        self.anchor_box.show_found(r["auto"])
        self.ecal_res = None; self.saved = False

    # --- Вкладка: Шкала ---

    def _build_ecal(self):
        """Построить вкладку энергетической шкалы."""
        w = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(w)

        layout.addWidget(QtWidgets.QLabel(   # реперы пары и ветвь ниже 100 кэВ
            "Энергетическая шкала фона по реперам пары: собственные структуры LaBr3 и линии, найденные в разности "
            "(Lu-176 -- по пресету). Ниже 100 кэВ -- отдельная ветвь (непропорциональность, K-край La)."
        ))
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Внешний источник (необязательно):"))
        self.src_edit = QtWidgets.QLineEdit(); self.src_edit.setPlaceholderText("спектр Am-241 или Ba-133 этого прибора")
        row.addWidget(self.src_edit)
        btn_src = QtWidgets.QPushButton("Обзор…")
        btn_src.clicked.connect(lambda: (fn := self._pick_file("Спектр источника")) and self.src_edit.setText(fn))
        row.addWidget(btn_src)
        self.src_combo = QtWidgets.QComboBox(); self.src_combo.addItems(["Am-241", "Ba-133"])
        row.addWidget(self.src_combo)
        layout.addLayout(row)

        layout.addWidget(self._btn("Построить шкалу", self._do_ecal))

        self.ecal_table = QtWidgets.QTableWidget()
        layout.addWidget(self.ecal_table)

        self.ecal_summary = QtWidgets.QLabel()
        self.ecal_summary.setWordWrap(True)
        self.ecal_summary.setTextFormat(QtCore.Qt.RichText)
        layout.addWidget(self.ecal_summary)
        self.low_plot = LowScalePlot()   # ветвь шкалы ниже 100 кэВ и полоса неопределённости
        layout.addWidget(self.low_plot, 1)

        return w

    def _do_ecal(self):
        """Запустить построение шкалы."""
        if self.h is None:
            self._error("Сначала выполните «Совместить»")
            return
        bg, sm = self._paths()
        self._run("Шкала", backend_pair.run_ecal, (bg, sm, self.h),
                  {"manual": self.manual, "reg_json": self.h_json, "source": self.src_edit.text().strip() or None,
                   "source_nuclide": self.src_combo.currentText()}, self._show_ecal)

    def _show_ecal(self, r):
        """Показать результат шкалы. реперы пары, ветвь ниже 100 кэВ, точки ПШПВ пары."""
        self.ecal = r["coeffs"]
        self.ecal_res = r
        self.sample_res = None   # «Шкала» пары -- режим пары
        self.ecal_points = backend_pair.markers_from(r)
        self.net = None
        fill_table(self.ecal_table, ["Репер", "Ветвь", "Канал фона", "± канал", "E, кэВ", "Невязка, кэВ", "LOO, кэВ"],
                   r["branch_rows"])
        self.ecal_summary.setText(backend_pair.ecal_html(r))
        self.low_plot.set_data(r)
        if r.get("fwhm_points"):
            self.fw_points_edit.setText(" ".join("%.2f:%.3f" % (p[0], p[1]) for p in r["fwhm_points"]))

    def _build_samplecal(self):
        """вкладка калибровки по пробе -- выбор пиков, реперы, шкала, совмещение фона."""
        w = QtWidgets.QWidget(); layout = QtWidgets.QVBoxLayout(w)
        self.picker = PeakPickerPanel()
        self.picker.message.connect(lambda m: self.log.append_text(m))   # журнал создаётся после вкладок
        self.picker.scale_requested.connect(self._do_sample_cal)
        self.buttons.extend(self.picker.buttons)
        layout.addWidget(self.picker, 3)
        self.sc_view = self._html_view()
        layout.addWidget(self.sc_view, 1)
        # таблица «файл / реперы» и применение шкалы по реперам решением пользователя
        self.sc_table = QtWidgets.QTableWidget()
        layout.addWidget(self.sc_table, 1)
        self.apply_scale_btn = QtWidgets.QPushButton("Записать шкалу по реперам")
        self.apply_scale_btn.setEnabled(False)
        self.apply_scale_btn.setToolTip("Доступно, только если шкала по реперам рекомендована; по умолчанию в XML -- калибровка файла пробы")
        self.apply_scale_btn.clicked.connect(lambda *a: self._do_sample_cal(use_reper_scale=True))
        layout.addWidget(self.apply_scale_btn)
        return w

    def _reset_sample_scale(self):
        """смена пары / «Новый спектр» -- кнопка гаснет, таблица и рекомендация сбрасываются."""
        self.apply_scale_btn.setEnabled(False)
        self.sc_table.setRowCount(0)
        self.sc_view.clear()

    def _do_sample_cal(self, use_reper_scale=False):
        """шкала по реперам пробы и видимым структурам LaBr3; результат заменяет «Совместить»+«Шкала»."""
        bg, sm = self._paths()
        if not bg or not sm:
            self._error("Выберите файл фона и файл образца")
            return
        if use_reper_scale and not (self.sample_res and self.sample_res["decision"].get("recommend")):
            self._error("Шкала по реперам не рекомендована -- сначала «Шкала по пробе»")
            return
        self._run("Шкала по пробе", backend_sample.run_sample_cal, (bg, sm, list(self.picker.repers)),
                  {"use_reper_scale": bool(use_reper_scale)}, self._show_sample_cal)

    def _show_sample_cal(self, r):
        """шкала XML и перевод фона на пробу; «Вычесть» строит разность на сетке пробы."""
        self.sample_res, self.h, self.h_json = r, [0.0, 1.0, 0.0], None
        self.ecal, self.ecal_res, self.net = r["coeffs"], r, None
        self.ecal_points = backend_pair.markers_from(r)
        fill_table(self.ecal_table, r["branch_headers"], r["branch_rows"])   # «файл / реперы»
        fill_table(self.sc_table, r["branch_headers"], r["branch_rows"])
        d = r["decision"]
        self.apply_scale_btn.setEnabled(bool(d.get("recommend")) and not d.get("replaced"))
        self.ecal_summary.setText(backend_sample.decision_html(r))
        self.sc_view.setHtml(backend_sample.decision_html(r))
        if r.get("fwhm_points"):
            self.fw_points_edit.setText(" ".join("%.2f:%.3f" % (p[0], p[1]) for p in r["fwhm_points"]))
        self.picker.set_scale(r["coeffs"] if r["decision"]["replaced"] else None)

    # --- Вкладка: Вычесть ---

    def _build_subtract(self):
        """Построить вкладку вычитания."""
        w = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(w)

        grid = QtWidgets.QGridLayout()

        # Строка 0: Папка
        grid.addWidget(QtWidgets.QLabel("Папка для результата:"), 0, 0)
        self.out_dir = QtWidgets.QLineEdit()
        grid.addWidget(self.out_dir, 0, 1)
        btn_dir = QtWidgets.QPushButton("Обзор…")
        btn_dir.clicked.connect(lambda: (d := self._pick_dir("Папка для результата")) and self.out_dir.setText(d))
        grid.addWidget(btn_dir, 0, 2)

        # Строка 1: Имя
        grid.addWidget(QtWidgets.QLabel("Имя файла (без .xml):"), 1, 0)
        self.out_name = QtWidgets.QLineEdit()
        grid.addWidget(self.out_name, 1, 1)

        # Строка 2: Описание
        desc = QtWidgets.QLabel(
            "Будут созданы два файла: «имя.xml» — разность с обнулёнными отрицательными каналами "
            "(для просмотра в BecqMoni) и «имя_со_знаком.xml» — разность со знаком (для проверок). "
            "В оба вписан блок ПШПВ (модель разрешения) — точки берутся с вкладки «ПШПВ»."
        )
        desc.setWordWrap(True)
        grid.addWidget(desc, 2, 0, 1, 3)
        # живое время выходного файла -- выбор, нормировку не отключает
        grid.addWidget(QtWidgets.QLabel("Живое время файла:"), 4, 0)
        self.lt_combo = QtWidgets.QComboBox()
        self.lt_combo.addItems(["образца (по умолчанию)", "фона"])
        self.lt_combo.setToolTip("Разность всегда в cps: net_cps = S/LT_s − B/LT_b. В файл пишется net_cps·LT выбранного "
                                 "спектра; cps в BecqMoni от выбора не меняется.")
        grid.addWidget(self.lt_combo, 4, 1)

        # Строка 3: Кнопки
        btn_sub = self._btn("Вычесть и сохранить", self._do_subtract)
        grid.addWidget(btn_sub, 3, 0)
        btn_open = QtWidgets.QPushButton("Открыть папку")
        btn_open.clicked.connect(self._open_out_dir)
        grid.addWidget(btn_open, 3, 1)

        layout.addLayout(grid)

        self.sub_view = self._html_view()
        layout.addWidget(self.sub_view)

        return w

    def _out_paths(self):
        """Получить пути к файлам результата."""
        d = self.out_dir.text().strip()
        n = self.out_name.text().strip()
        if not d or not n:
            raise backend.UserError("Укажите папку и имя файла результата")
        if n.lower().endswith(".xml"):
            n = n[:-4]
        return (str(pathlib.Path(d) / (n + ".xml")), str(pathlib.Path(d) / (n + "_со_знаком.xml")))

    def _do_subtract(self):
        """Запустить вычитание."""
        if self.h is None:
            self._error("Сначала выполните «Совместить»")
            return
        if self.ecal is None:
            self._error("Сначала выполните «Шкала»")
            return

        try:
            points = backend.parse_points(self.fw_points_edit.text())
            z_path, s_path = self._out_paths()
        except backend.UserError as exc:
            self._error(str(exc))
            return

        d = self.out_dir.text().strip()
        if not pathlib.Path(d).is_dir():
            self._error("Папка не существует: " + d)
            return

        existing = [p for p in (z_path, s_path) if pathlib.Path(p).exists()]
        if existing:
            ans = QtWidgets.QMessageBox.question(
                self, "Файлы существуют", "Перезаписать?\n" + "\n".join(existing)
            )
            if ans != QtWidgets.QMessageBox.Yes:
                return

        bg, sm = self._paths()
        lt = "bg" if self.lt_combo.currentIndex() == 1 else "sample"
        if self.sample_res is not None:   # режим «по пробе» -- разность на сетке пробы
            self._run("Вычитание (на шкале пробы)", backend_sample.run_sample_subtract,
                      (bg, sm, self.sample_res, points, z_path, s_path), {"lt_out": lt}, self._show_subtract)
            return
        self._run(
            "Вычитание", backend_pair.run_subtract,   # составная шкала, модель линий образца
            (bg, sm, self.h, self.ecal_res or self.ecal, points, z_path, s_path),
            {"reg_json": self.h_json, "manual": self.manual,
             "lt_out": "bg" if self.lt_combo.currentIndex() == 1 else "sample"}, self._show_subtract
        )

    def _show_subtract(self, r):
        """Показать результат вычитания. проверки окон собственных структур вместо «сумма = 0»."""
        self.net = r
        self.saved = True
        self.plot.set_data(r, self.ecal_points)
        first_name = r["files"][0]["name"] if r["files"] else "результат"
        self.sub_view.setHtml(backend_pair.subtract_html(r) +
                              "<p>Откройте «%s» в BecqMoni. График — на вкладке «График».</p>" % html.escape(first_name))
        self.tabs.setCurrentWidget(self.plot_tab)

    def _open_out_dir(self):
        """Открыть папку результата в проводнике."""
        d = self.out_dir.text().strip()
        if d and pathlib.Path(d).is_dir():
            QtGui.QDesktopServices.openUrl(QtCore.QUrl.fromLocalFile(d))

    # --- Вкладка: График ---

    def _build_plot(self):
        """Построить вкладку графика."""
        self.plot = SpectrumPlot()
        self.plot_tab = self.plot
        return self.plot

    # --- Вкладка: ПШПВ ---

    def _build_fwhm(self):
        """Построить вкладку ПШПВ."""
        w = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(w)

        label = QtWidgets.QLabel(
            "ПШПВ (ширина пика на полувысоте) — модель разрешения прибора. "
            "Точки «канал:ширина в каналах» через пробел. "
            "По умолчанию — точки, измеренные на этом приборе (AtomNano 5 PRO LaBr3, 2048 каналов)."
        )
        label.setWordWrap(True)
        layout.addWidget(label)

        self.fw_points_edit = QtWidgets.QLineEdit()
        self.fw_points_edit.setText(" ".join("%.2f:%.2f" % p for p in backend.default_fwhm_points()))
        layout.addWidget(self.fw_points_edit)

        btn_row = QtWidgets.QHBoxLayout()
        btn_row.addWidget(self._btn("Построить модель", self._do_fwhm))
        btn_row.addWidget(self._btn("Значения для настроек прибора BecqMoni", self._do_device))
        btn_copy = QtWidgets.QPushButton("Копировать в буфер")
        btn_copy.clicked.connect(self._copy_device)
        btn_row.addWidget(btn_copy)
        btn_row.addStretch()
        layout.addLayout(btn_row)

        self.fw_plot = FwhmPlot()
        layout.addWidget(self.fw_plot, 1)

        self.fw_table = QtWidgets.QTableWidget()
        self.fw_table.setMaximumHeight(180)
        layout.addWidget(self.fw_table)

        self.fw_view = QtWidgets.QLabel()
        self.fw_view.setWordWrap(True)
        self.fw_view.setTextFormat(QtCore.Qt.RichText)
        self.fw_view.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        layout.addWidget(self.fw_view)

        self._device_text = ""

        return w

    def _points_or_error(self):
        """Парсить точки ПШПВ или показать ошибку."""
        try:
            return backend.parse_points(self.fw_points_edit.text())
        except backend.UserError as exc:
            self._error(str(exc))
            return None

    def _do_fwhm(self):
        """Запустить построение модели ПШПВ."""
        pts = self._points_or_error()
        if pts is None:
            return
        self._run("Модель ПШПВ", backend.run_fwhm_model, (pts,), {"n_ch": self.n_ch}, lambda r: self._show_fwhm(pts, r))

    def _show_fwhm(self, pts, r):
        """Показать результат модели ПШПВ."""
        self.fw_plot.set_data(pts, r["x"], r["y"])
        fill_table(
            self.fw_table,
            ["Канал", "Измерено, кан", "Модель, кан", "Невязка, %"],
            r["rows"]
        )
        self.fw_view.setText(
            status_html(r["ok"], "Модель годна для BecqMoni" if r["ok"] else "Модель НЕ годна") +
            "<br>Коэффициенты c0, c1, c2 = " + ", ".join("%.6g" % c for c in r["coeffs"]) +
            "<br>ПШПВ на каналах 0, 500, 1309, N−1: " + ", ".join("%.2f" % v for v in r["values"])
        )

    def _do_device(self):
        """Запустить расчёт параметров прибора."""
        pts = self._points_or_error()
        if pts is None:
            return
        self._run("Параметры прибора", backend.run_device, (pts,),   # канал отсчёта -- в каналах прибора
                  {"ch_ref": int(round(1309 * self.n_ch / 2048)), "n_ch": self.n_ch}, self._show_device)

    def _show_device(self, r):
        """Показать параметры прибора."""
        self._device_text = r["clipboard"]
        self.fw_view.setText(
            "<b>Для настроек прибора BecqMoni:</b><br>"
            "FWHM_AT_0 = <b>%s</b><br>"
            "Ch_Fwhm = <b>%s</b><br>"
            "Width_Fwhm = <b>%s</b><br><br>"
            "%s" % (
                r["FWHM_AT_0"], r["Ch_Fwhm"], r["Width_Fwhm"],
                html.escape(r["hint"]).replace("\n", "<br>")
            )
        )

    def _copy_device(self):
        """Скопировать параметры в буфер обмена."""
        if not self._device_text:
            self._error("Сначала нажмите «Значения для настроек прибора BecqMoni»")
            return
        QtWidgets.QApplication.clipboard().setText(self._device_text)
        self.statusBar().showMessage("Скопировано в буфер обмена")

    # --- Вкладка: Фоны ---

    def _build_backgrounds(self):
        """Построить вкладку фонов."""
        w = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(w)

        # Группа 1: Проверка папки
        group1 = QtWidgets.QGroupBox("Проверка папки фонов (дубли, группы)")
        g1_layout = QtWidgets.QVBoxLayout(group1)

        row1 = QtWidgets.QHBoxLayout()
        self.inv_dir = QtWidgets.QLineEdit()
        row1.addWidget(self.inv_dir)
        btn_inv_dir = QtWidgets.QPushButton("Обзор…")
        btn_inv_dir.clicked.connect(lambda: (d := self._pick_dir("Папка фонов")) and self.inv_dir.setText(d))
        row1.addWidget(btn_inv_dir)
        self.inv_rec = QtWidgets.QCheckBox("с подпапками")
        row1.addWidget(self.inv_rec)
        row1.addWidget(self._btn("Проверить папку", self._do_inventory))
        g1_layout.addLayout(row1)

        self.inv_table = QtWidgets.QTableWidget()
        g1_layout.addWidget(self.inv_table)

        self.inv_view = QtWidgets.QLabel()
        self.inv_view.setWordWrap(True)
        self.inv_view.setTextFormat(QtCore.Qt.RichText)
        g1_layout.addWidget(self.inv_view)

        layout.addWidget(group1)

        # Группа 2: Разложение
        group2 = QtWidgets.QGroupBox("Разложение фона: доля K-40 под сумм-пиком (по желанию)")
        g2_layout = QtWidgets.QVBoxLayout(group2)

        row2 = QtWidgets.QHBoxLayout()
        row2.addWidget(QtWidgets.QLabel("Экранированный фон:"))
        self.dec_ref = QtWidgets.QLineEdit()
        row2.addWidget(self.dec_ref)
        btn_dec_ref = QtWidgets.QPushButton("Обзор…")
        btn_dec_ref.clicked.connect(lambda: (fn := self._pick_file("Экранированный фон")) and self.dec_ref.setText(fn))
        row2.addWidget(btn_dec_ref)
        g2_layout.addLayout(row2)

        row3 = QtWidgets.QHBoxLayout()
        row3.addWidget(QtWidgets.QLabel("Неэкранированные фоны:"))
        self.dec_list = QtWidgets.QLineEdit()
        row3.addWidget(self.dec_list)
        btn_dec_list = QtWidgets.QPushButton("Обзор…")
        btn_dec_list.clicked.connect(lambda: (fns := self._pick_files("Неэкранированные фоны")) and self.dec_list.setText(" | ".join(fns)))
        row3.addWidget(btn_dec_list)
        g2_layout.addLayout(row3)

        g2_layout.addWidget(self._btn("Разложить", self._do_decompose))

        self.dec_view = self._html_view()
        g2_layout.addWidget(self.dec_view)

        layout.addWidget(group2)

        return w

    def _do_inventory(self):
        """Запустить проверку папки фонов."""
        d = self.inv_dir.text().strip()
        self._run(
            "Проверка фонов", backend.run_inventory, (d,),
            {"recursive": self.inv_rec.isChecked()}, self._show_inventory
        )

    def _show_inventory(self, r):
        """Показать результат проверки фонов."""
        rows = [
            (pathlib.Path(x["path"]).name, x["n_raw"], x["live"], x["start"], x["rate"], x["md5"])
            for x in r["rows"]
        ]
        fill_table(
            self.inv_table,
            ["Файл", "Каналов", "Живое время, с", "Начало", "Скорость счёта, 1/с", "Отпечаток"],
            rows
        )

        parts = ["Файлов: %d, различных измерений: %d, групп дублей: %d" % (
            r["n_files"], r["n_unique"], len(r["dup_groups"])
        )]
        for g in r["dup_groups"]:
            parts.append("Одинаковый счёт: " + " = ".join(html.escape(pathlib.Path(p).name) for p in g))
        self.inv_view.setText("<br>".join(parts))

    def _do_decompose(self):
        """Запустить разложение фона."""
        ref = self.dec_ref.text().strip()
        lst = [s.strip() for s in self.dec_list.text().split("|") if s.strip()]
        self._run(
            "Разложение фона", backend.run_decompose, (ref, lst),
            # шкала «Шкалы» относится к фону пары -- передаётся, только если эталон и есть этот фон
            {"ecal": self.ecal if ref and ref == self._paths()[0] else None}, self._show_decompose
        )

    def _show_decompose(self, r):
        """Показать результат разложения."""
        self.dec_view.setHtml(
            "<h3>Итог разложения</h3><pre>" +
            html.escape(json.dumps(r["summaries"], ensure_ascii=False, indent=1, default=str)) +
            "</pre><p>Подробности — в журнале ниже.</p>"
        )

    # --- Вкладка: Самопроверка ---

    def _build_selftest(self):
        """Построить вкладку самопроверки."""
        w = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(w)

        label = QtWidgets.QLabel(
            "Мутационная самопроверка метода на выбранной паре: метод намеренно портится разными способами, "
            "и каждая порча должна быть поймана. Занимает несколько минут; окно при этом не замирает."
        )
        label.setWordWrap(True)
        layout.addWidget(label)

        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Дополнительные фоны (необязательно):"))
        self.st_alt = QtWidgets.QLineEdit()
        row.addWidget(self.st_alt)
        btn_alt = QtWidgets.QPushButton("Обзор…")
        btn_alt.clicked.connect(lambda: (fns := self._pick_files("Дополнительные фоны")) and self.st_alt.setText(" | ".join(fns)))
        row.addWidget(btn_alt)
        layout.addLayout(row)

        layout.addWidget(self._btn("Запустить самопроверку", self._do_selftest))

        self.st_view = self._html_view()
        layout.addWidget(self.st_view)

        return w

    def _do_selftest(self):
        """Запустить самопроверку."""
        bg, sm = self._paths()
        alt = [s.strip() for s in self.st_alt.text().split("|") if s.strip()]
        self._run(
            "Самопроверка", backend.run_selftest, (bg, sm),
            {"alt_bg": alt, "manual": self.manual}, self._show_selftest
        )

    def _show_selftest(self, r):
        """Показать результат самопроверки."""
        parts = [
            "<h2>" + status_html(
                r["green"],
                "ЗЕЛЁНЫЙ: метод на этой паре работает, все порчи пойманы" if r["green"]
                else "КРАСНЫЙ: есть непойманные порчи или неожиданные провалы"
            ) + "</h2>"
        ]
        parts.extend(html.escape(x) for x in r["itog"])
        if r["interp"]:
            parts.append("<h3>Требует толкования</h3>")
            parts.extend(html.escape(x) for x in r["interp"])
        self.st_view.setHtml("<br>".join(parts))

    # --- Закрытие окна ---

    def closeEvent(self, event):
        """Обработка закрытия окна."""
        if self.job is not None and self.job.isRunning():
            ans = QtWidgets.QMessageBox.question(
                self, "Операция выполняется", "Идёт расчёт. Закрыть программу всё равно?"
            )
            if ans != QtWidgets.QMessageBox.Yes:
                event.ignore()
                return
        event.accept()
