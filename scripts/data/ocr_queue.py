"""
OCR queue: запускает UL-OCR последовательно на очереди документов.
Ждёт завершения предыдущего прогона перед запуском следующего.

Запуск: python ocr_queue.py
"""
import subprocess
import sys
import time
from pathlib import Path

PY = Path(r"<USERPROFILE>\unlimited-ocr-eval\.venv\Scripts\python.exe")
RUN_OCR = Path(r"<USERPROFILE>\.claude\skills\unlimited-ocr\scripts\run_ocr.py")
BASE = Path(r"<WORKDIR>\gamma-spectrum-analysis")
OUT_ROOT = BASE / "references" / "_ocr_output"

QUEUE = [
    # (pdf_path, out_dir, sentinel_page, pages)
    # sentinel_page: имя подпапки последней страницы (page_NN), если None — ещё не завершён
    (
        BASE / "books_library" / "Документация ЛСРМ" / "01_methodology_pdf"
        / "Активность в счетных образцах. Методика измерений на гамма-спектрометрах с использоваонием ПО СпектраЛайн.pdf",
        OUT_ROOT / "Aktivnost_v_schetnyh_obrazcah_ocr",
        "all",
    ),
    (
        BASE / "books_library" / "Документация ЛСРМ" / "02_topical_pdf"
        / "Прецизионные измерения.pdf",
        OUT_ROOT / "Pretsizionnye_izmerenia_ocr",
        "all",
    ),
    (
        BASE / "books_library" / "Документация ЛСРМ" / "02_topical_pdf"
        / "5_2_Практическая спектрометрия-ядерные материалы.pdf",
        OUT_ROOT / "5_2_Prakticheskaya_spektrometria_ocr",
        "all",
    ),
    (
        BASE / "detectors" / "RadiaCode_103" / "references"
        / "Руководство_спектроскописта_V1.05_Соловьев_2024.pdf",
        OUT_ROOT / "Rukovodstvo_spektroskopista_ocr",
        "all",
    ),
]

# Ждать завершения первого прогона (Lsrm_algorithmic_foundations)
FIRST_SENTINEL = BASE / "references" / "Lsrm_algorithmic_foundations_ocr" / "page_54" / "result.md"
WAIT_TIMEOUT = 90 * 60  # 90 минут макс
POLL_INTERVAL = 60       # проверять раз в минуту


def wait_for_file(path: Path, timeout: int, interval: int) -> bool:
    elapsed = 0
    while elapsed < timeout:
        if path.exists():
            return True
        print(f"[queue] Жду завершения первого OCR... {elapsed//60} мин прошло", flush=True)
        time.sleep(interval)
        elapsed += interval
    return False


def run_ocr(pdf: Path, out: Path, pages: str) -> int:
    out.mkdir(parents=True, exist_ok=True)
    print(f"\n[queue] → {pdf.name}", flush=True)
    result = subprocess.run(
        [str(PY), str(RUN_OCR), str(pdf), "--pages", pages, "--dpi", "200", "--out", str(out)],
        env={**__import__("os").environ, "PYTHONIOENCODING": "utf-8"},
    )
    return result.returncode


def main():
    print("[queue] Ждём page_54/result.md от первого OCR...", flush=True)
    if not wait_for_file(FIRST_SENTINEL, WAIT_TIMEOUT, POLL_INTERVAL):
        print("[queue] TIMEOUT: первый OCR не завершился за 90 мин. Продолжаю всё равно.", flush=True)

    print("[queue] Первый OCR завершён. Запускаю очередь.", flush=True)

    for pdf, out, pages in QUEUE:
        if not pdf.exists():
            print(f"[queue] SKIP (не найден): {pdf.name}", flush=True)
            continue
        rc = run_ocr(pdf, out, pages)
        if rc != 0:
            print(f"[queue] ERROR rc={rc}: {pdf.name}", flush=True)
        else:
            print(f"[queue] OK: {pdf.name} → {out}", flush=True)

    print("[queue] Очередь завершена.", flush=True)


if __name__ == "__main__":
    main()
