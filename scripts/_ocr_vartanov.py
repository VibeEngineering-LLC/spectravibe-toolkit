from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""Render + OCR all 275 pages of Вартанов DJVU (converted to PDF earlier).
Pages are rendered at 2x zoom, OCR'd via Tesseract rus+eng in parallel,
then concatenated into a single corpus markdown file.
"""
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PDF = ROOT / "references" / "_converted_tmp" / "Vartanov.pdf"
OUT_MD = ROOT / "references" / "_extracted_corpus" / "Документация ЛСРМ" / "02_topical_pdf" / "Vartanov_Prakticheskie_metody_scintillyatsionnoy_gamma-spektrometrii.djvu.md"
# F-333 / v1.18.18.6: PNG originals → books_library/_corpus_pages/
# (out of release archive via F-293 exclude list).
TMP_PAGES = ROOT / "books_library" / "_corpus_pages" / "vartanov"
TESSERACT = r"C:\Program Files\Tesseract-OCR\tesseract.exe"
TESSDATA = os.environ.get("USERPROFILE", "") + r"\tessdata"


def render_page(args) -> tuple[int, str]:
    pno, pdf_path, out_dir = args
    import fitz
    doc = fitz.open(pdf_path)
    page = doc[pno]
    pix = page.get_pixmap(matrix=fitz.Matrix(1.5, 1.5), alpha=False)  # smaller zoom for speed
    png = Path(out_dir) / f"page_{pno+1:03d}.png"
    pix.save(str(png))
    doc.close()
    return pno + 1, str(png)


def ocr_page(args) -> tuple[int, str]:
    pno, png_path = args
    proc = subprocess.run(
        [TESSERACT, png_path, "-", "-l", "rus+eng", "--tessdata-dir", TESSDATA, "--psm", "1"],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=60,
    )
    return pno, proc.stdout or ""


def main() -> int:
    import fitz
    doc = fitz.open(str(PDF))
    n = doc.page_count
    doc.close()
    TMP_PAGES.mkdir(parents=True, exist_ok=True)
    workers = max(1, (os.cpu_count() or 4) - 1)
    print(f"Render+OCR {n} pages on {workers} workers", flush=True)

    t0 = time.time()
    # Render stage
    args_list = [(i, str(PDF), str(TMP_PAGES)) for i in range(n)]
    pngs: list[tuple[int, str]] = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for fut in as_completed([ex.submit(render_page, a) for a in args_list]):
            pngs.append(fut.result())
    print(f"  Render done in {round(time.time()-t0,1)}s")

    # OCR stage
    t1 = time.time()
    results: dict[int, str] = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        ocr_args = [(p, png) for p, png in pngs]
        for fut in as_completed([ex.submit(ocr_page, a) for a in ocr_args]):
            pno, txt = fut.result()
            results[pno] = txt
            if pno % 20 == 0:
                print(f"  OCR p{pno:3d}  {len(txt)} chars", flush=True)
    print(f"  OCR done in {round(time.time()-t1,1)}s")

    parts = ["<!-- src: Вартанов DJVU 275 стр., OCR via Tesseract -->", "", "# Практические методы сцинтилляционной гамма-спектрометрии (Вартанов)", ""]
    for i in sorted(results):
        parts.append(f"\n## PAGE {i}\n\n{results[i].strip()}")
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(parts), encoding="utf-8")
    total_chars = sum(len(t) for t in results.values())
    print(f"\nDONE total {round(time.time()-t0,1)}s, {total_chars:,} chars, {OUT_MD.stat().st_size//1024} KB md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
