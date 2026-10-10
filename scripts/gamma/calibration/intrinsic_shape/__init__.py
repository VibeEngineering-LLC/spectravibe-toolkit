"""Калибровка по форме собственных структур спектра (intrinsic-shape calibration).

Совмещение шкалы образца со шкалой фона по форме собственной активности детектора
(LaBr3: La-138, Ac-227; LYSO/LuAG: Lu-176) с последующим вычитанием фона, моделью ПШПВ
и энергетической шкалой по реперам. Командная строка:

    python -m gamma.calibration.intrinsic_shape <подкоманда> ...

Методика: references/intrinsic_shape_calibration.md; скилл: skills/intrinsic-shape-cal/SKILL.md.
"""
