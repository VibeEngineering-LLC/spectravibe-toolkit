from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""OCR the 26 rendered pages of 'Активность в счетных образцах' via Tesseract.
Parallel across all CPU cores.
"""
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# F-333 / v1.18.18.6: PNG originals → books_library/_corpus_pages/
# (out of release archive via F-293 exclude list).
PAGES_DIR = ROOT / "books_library" / "_corpus_pages" / "lsrm_activity_2014"
OUT_MD = ROOT / "references" / "_extracted_corpus" / "Документация ЛСРМ" / "01_methodology_pdf" / "Активность в счетных образцах. Методика измерений на гамма-спектрометрах с использоваонием ПО СпектраЛайн.pdf.md"
TESSERACT = r"C:\Program Files\Tesseract-OCR\tesseract.exe"
TESSDATA = os.environ.get("USERPROFILE", "") + r"\tessdata"


def ocr_one(png_path: str) -> tuple[str, str]:
    name = Path(png_path).name
    proc = subprocess.run(
        [TESSERACT, png_path, "-", "-l", "rus+eng", "--tessdata-dir", TESSDATA, "--psm", "1"],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    return name, proc.stdout or f"[OCR ERROR: {proc.stderr}]"


def main() -> int:
    pages = sorted(PAGES_DIR.glob("page_*.png"))
    print(f"OCR {len(pages)} pages via Tesseract rus+eng…", flush=True)
    workers = max(1, (os.cpu_count() or 4) - 1)
    t0 = time.time()
    results: dict[str, str] = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(ocr_one, str(p)): p for p in pages}
        for fut in as_completed(futs):
            name, txt = fut.result()
            results[name] = txt
            print(f"  OK  {name}  {len(txt)} chars", flush=True)
    elapsed = round(time.time() - t0, 1)
    parts = ["<!-- src: Документация ЛСРМ/01_methodology_pdf/Активность в счетных образцах… -->",
             "<!-- ext: .pdf (OCR via Tesseract rus+eng) -->",
             "<!-- meta: scanned PDF, OCR-extracted -->",
             "",
             "# Активность в счетных образцах. Методика измерений на гамма-спектрометрах с использованием ПО СпектраЛайн",
             ""]
    for name in sorted(results):
        i = int(name.split("_")[1].split(".")[0])
        parts.append(f"\n## PAGE {i}\n\n{results[name].strip()}")
    OUT_MD.write_text("\n".join(parts), encoding="utf-8")
    total = sum(len(t) for t in results.values())
    print(f"\nDONE in {elapsed}s — {total:,} chars total")
    return 0


if __name__ == "__main__":
    sys.exit(main())
