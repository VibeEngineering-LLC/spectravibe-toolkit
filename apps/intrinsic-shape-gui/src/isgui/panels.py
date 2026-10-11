"""Панели для вкладки «Шкала»: график низкой шкалы и блок опор автопресета.
Сгенерировано gen_code.py (qwen3.8:27b) по specs/panels.md, 2026-10-10"""

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets


class LowScalePlot(QtWidgets.QWidget):
    """График энергетической шкалы в области низких энергий."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.plot = pg.PlotWidget(background="w")
        self.plot.showGrid(x=True, y=True, alpha=0.3)
        self.plot.addLegend()
        self.plot.setLabel("bottom", "Канал (шкала фона)")
        self.plot.setLabel("left", "Энергия, кэВ")

        self.note = QtWidgets.QLabel()
        self.note.setWordWrap(True)

        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self.plot, 1)
        layout.addWidget(self.note)
        self.setMinimumHeight(260)

    def set_data(self, ecal):
        """Заполняет график данными из словаря ecal."""
        self.plot.clear()
        if self.plot.plotItem.legend is not None:   # легенда не очищается plot.clear()
            self.plot.plotItem.legend.clear()

        if "low" not in ecal or not ecal["low"]["band"]["ch"]:
            self.note.setText("Ветвь шкалы ниже 100 кэВ не построена")
            return

        low = ecal["low"]
        ch = np.asarray(low["band"]["ch"], dtype=float)
        E = np.asarray(low["band"]["E"], dtype=float)
        sg = np.asarray(low["band"]["sigma"], dtype=float)

        # Оставляем каналы, где энергия <= 130 кэВ
        mask = E <= 130
        ch = ch[mask]
        E = E[mask]
        sg = sg[mask]

        # Вычисление основной шкалы (полином, коэффициенты по возрастанию степени)
        coeffs = ecal["coeffs"]
        main = np.zeros_like(ch)
        for i, c in enumerate(coeffs):
            main += c * (ch ** i)

        # Полоса неопределенности
        e_upper = E + sg
        e_lower = E - sg
        curve_upper = self.plot.plot(ch, e_upper, pen=pg.mkPen(None))
        curve_lower = self.plot.plot(ch, e_lower, pen=pg.mkPen(None))
        fill = pg.FillBetweenItem(curve_upper, curve_lower, brush=(100, 150, 230, 70))
        self.plot.addItem(fill)
        self.plot.plotItem.legend.addItem(curve_upper, "Полоса неопределённости")

        # Композитная шкала E
        curve_e = self.plot.plot(ch, E, pen=pg.mkPen(color=(0, 0, 255), width=2))
        self.plot.plotItem.legend.addItem(curve_e, "Шкала с ветвью ниже 100 кэВ")

        # Основная шкала (полином)
        curve_main = self.plot.plot(ch, main, pen=pg.mkPen(color=(128, 128, 128), style=QtCore.Qt.DashLine))
        self.plot.plotItem.legend.addItem(curve_main, "Основная шкала (полином)")

        # Кривая из XML (если есть)
        if "coeffs_xml" in ecal:
            coeffs_xml = ecal["coeffs_xml"]
            xml_curve = np.zeros_like(ch)
            for i, c in enumerate(coeffs_xml):
                xml_curve += c * (ch ** i)
            curve_xml = self.plot.plot(ch, xml_curve, pen=pg.mkPen(color=(255, 165, 0), width=1))
            self.plot.plotItem.legend.addItem(curve_xml, "Записывается в XML")

        # Реперы (точки)
        points = low.get("points", [])
        if points:
            px = np.array([p["ch"] for p in points], dtype=float)
            py = np.array([p["E"] for p in points], dtype=float)
            curve_points = self.plot.plot(px, py, pen=None, symbol="o", symbolSize=9, symbolBrush=(220, 50, 50))
            self.plot.plotItem.legend.addItem(curve_points, "Реперы")

            # Текстовые метки для каждого репера
            for p in points:
                text_item = pg.TextItem(p["name"], color=(120, 0, 120), anchor=(0, 1))
                self.plot.addItem(text_item)
                text_item.setPos(p["ch"], p["E"])

        # Вертикальная линия на 37 кэВ (если есть предупреждение об экстраполяции)
        if low.get("warn_extrapolation"):
            # Находим канал, где E пересекает 37 кэВ (предполагаем, что E возрастает)
            ch_37 = np.interp(37, E, ch)
            line = pg.InfiniteLine(pos=ch_37, angle=90, pen=pg.mkPen(color=(220, 50, 50), style=QtCore.Qt.DashLine),
                                   label="37 кэВ: ниже -- экстраполяция", labelOpts={"position": 0.9})
            self.plot.addItem(line)

        self.note.setText(
            "s = %.4f ± %.4f; полоса: стат. + форма + K-край; порог ≈ %.1f кэВ"
            % (low["s"], low["sigma_s"], low["e_threshold"])
        )

    def clear_data(self):
        """Очищает график и заметку."""
        self.plot.clear()
        if self.plot.plotItem.legend is not None:
            self.plot.plotItem.legend.clear()
        self.note.setText("")


class AnchorBox(QtWidgets.QGroupBox):
    """Блок управления опорами автопресета и линиями образца."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setTitle("Опоры автопресета (каналы фона) и линии образца")

        grid = QtWidgets.QGridLayout(self)

        # Строка 0: Ручная поправка
        grid.addWidget(QtWidgets.QLabel("Ручная поправка опор:"), 0, 0)
        self.manual_edit = QtWidgets.QLineEdit()
        self.manual_edit.setPlaceholderText("пусто -- автоматически; или «A37 SUM [EDGE]», например 110.4 3183.7")
        grid.addWidget(self.manual_edit, 0, 1)

        # Строка 1: Линии образца
        grid.addWidget(QtWidgets.QLabel("Линии образца:"), 1, 0)
        self.robust_combo = QtWidgets.QComboBox()
        self.robust_combo.addItem("автоматически", "auto")
        self.robust_combo.addItem("учитывать (устойчивая регистрация)", "on")
        self.robust_combo.addItem("не учитывать (простая регистрация)", "off")
        grid.addWidget(self.robust_combo, 1, 1)

        # Строка 2: Найдено
        self.found_label = QtWidgets.QLabel("Опоры будут найдены при совмещении")
        self.found_label.setWordWrap(True)
        self.found_label.setTextFormat(QtCore.Qt.RichText)
        grid.addWidget(self.found_label, 2, 0, 1, 2)

    def manual_text(self):
        """Возвращает текст ручной поправки."""
        return self.manual_edit.text().strip()

    def robust_mode(self):
        """Возвращает текущий режим устойчивой регистрации."""
        return self.robust_combo.currentData()

    def show_found(self, auto):
        """Отображает найденные опоры и параметры совмещения."""
        anchors = auto.get("anchors", {})
        identity = auto.get("identity", False)
        r = auto.get("r", 0.0)

        # Извлекаем значения для текста
        a37 = anchors.get("a37", {}).get("x", 0.0)
        sum_peak = anchors.get("sum", {}).get("x", 0.0)
        beta_edge = anchors.get("edge", {}).get("x", 0.0)

        text = "Найдено: пик ~37 кэВ — %.2f, сумм-пик La-138 — %.2f, β-край — %.2f; масштаб к эталону %.4f" % (
            a37, sum_peak, beta_edge, r
        )

        if identity:
            text += " (тождественно эталону)"

        # Проверяем неуверенные опоры
        uncertain_names = []
        for name, data in anchors.items():
            if not data.get("ok", True):
                uncertain_names.append(name)

        if uncertain_names:
            names_str = ", ".join(uncertain_names)
            text += ' <span style="color:#c62828"> — НЕУВЕРЕННО: %s</span>' % names_str

        self.found_label.setText(text)

    def reset(self):
        """Сбрасывает поля к значениям по умолчанию."""
        self.manual_edit.clear()
        self.robust_combo.setCurrentIndex(0)
        self.found_label.setText("Опоры будут найдены при совмещении")
