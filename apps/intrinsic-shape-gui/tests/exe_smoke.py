"""Приёмка exe: запуск, окно показано (файл-метка из __main__), время старта, снимок окна, размер папки сборки."""
import os, sys, time, pathlib, subprocess, tempfile
sys.stdout.reconfigure(encoding="utf-8")
DIST = pathlib.Path(__file__).resolve().parents[1] / "dist" / "IntrinsicShapeGUI"
EXE = DIST / "IntrinsicShapeGUI.exe"
tmp = pathlib.Path(tempfile.mkdtemp(prefix="isgui_exe_"))
for i in range(int(sys.argv[1]) if len(sys.argv) > 1 else 3):
    mark, png = tmp / ("shown%d.txt" % i), tmp / ("window%d.png" % i)
    env = dict(os.environ, ISGUI_SMOKE_FILE=str(mark), ISGUI_SMOKE_PNG=str(png), ISGUI_SMOKE_HOLD_MS="1500")
    t0 = time.time(); proc = subprocess.Popen([str(EXE)], env=env)
    while not mark.exists() and time.time() - t0 < 60 and proc.poll() is None:
        time.sleep(0.02)
    t_shown = time.time() - t0
    rc = proc.wait(60)
    print("[RUN %d] окно показано: %s, от запуска exe до показа %.2f с, метка: %s, rc=%s, снимок: %s" % (
        i, mark.exists(), t_shown, mark.read_text("utf-8") if mark.exists() else "-", rc, png if png.exists() else "нет"))
size = sum(f.stat().st_size for f in DIST.rglob("*") if f.is_file())
print("[SIZE] папка сборки %s: %.1f МБ, файлов %d" % (DIST, size / 2**20, sum(1 for f in DIST.rglob("*") if f.is_file())))
