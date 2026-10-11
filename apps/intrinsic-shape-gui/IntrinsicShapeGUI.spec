# PyInstaller spec -- IntrinsicShapeGUI (one-folder). Сборка: build.bat. Модуль toolkit берётся из scripts/ корня репозитория.
import os
from PyInstaller.utils.hooks import collect_submodules
ROOT = os.path.abspath(SPECPATH)
TK = os.path.abspath(os.path.join(ROOT, "..", ".."))   # корень репозитория toolkit
SCRIPTS = os.path.join(TK, "scripts")
PRESETS = os.path.join(SCRIPTS, "gamma", "calibration", "intrinsic_shape", "presets")
IS = "gamma.calibration.intrinsic_shape."
a = Analysis([os.path.join(ROOT, "src", "isgui", "__main__.py")], pathex=[os.path.join(ROOT, "src"), SCRIPTS],
    datas=[(PRESETS, "gamma/calibration/intrinsic_shape/presets"),   # ecal/selftest ищут пресеты рядом с модулем
           (os.path.join(TK, "data", "nuclides.json"), "data"),   # ядерные данные toolkit (реперы пары)
           (os.path.join(TK, "data", "xray_lines_lsrm.csv"), "data")],
    hiddenimports=[IS + m for m in ("core", "register", "ops", "ecal", "peakfit", "fwhm", "netxml", "background",
                                    "selftest", "selftest_boot", "autopreset", "robust", "pairref", "lowscale", "pipeline", "verdict",
                                    "samplecal", "samplecal_q", "samplecal_net")]   # модули режима «по пробе» (ленивый импорт)
                  + ["gamma.data", "gamma.data.nuclide_library"]
                  + collect_submodules("pyqtgraph", filter=lambda n: not any(s in n for s in ("opengl", "examples", "jupyter"))),
    excludes=["tkinter", "matplotlib", "PIL", "pandas", "PyQt5", "PyQt6", "PySide2", "torch", "IPython",
              "pyqtgraph.opengl", "OpenGL", "PySide6.QtWebEngineCore", "PySide6.Qt3DCore"],
    noarchive=False)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="IntrinsicShapeGUI", console=False, upx=False)
coll = COLLECT(exe, a.binaries, a.datas, upx=False, name="IntrinsicShapeGUI")
