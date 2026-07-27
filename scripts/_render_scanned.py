from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""Render every page of the LSRM 'Активность в счетных образцах' (scanned PDF)
into PNG for multimodal sub-agent reading.

Output: books_library/_corpus_pages/lsrm_activity_2014/page_NN.png

F-333 / v1.18.18.6 — PNG originals relocated out of references/_extracted_corpus/
into books_library/_corpus_pages/ (excluded from release archives via F-293).
F-293 also moved the source PDF: references/books/ → books_library/.
"""
import fitz
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# F-293: source PDF lives in books_library/ (was references/books/).
SRC = ROOT / "books_library" / "Документация ЛСРМ" / "01_methodology_pdf" / "Активность в счетных образцах. Методика измерений на гамма-спектрометрах с использоваонием ПО СпектраЛайн.pdf"
# F-333: PNG originals → books_library/_corpus_pages/ (out of release archive).
OUT = ROOT / "books_library" / "_corpus_pages" / "lsrm_activity_2014"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    doc = fitz.open(str(SRC))
    print(f"Rendering {doc.page_count} pages")
    zoom = 2.0   # ~144 dpi — good balance for OCR/vision
    mat = fitz.Matrix(zoom, zoom)
    for i, page in enumerate(doc, 1):
        pix = page.get_pixmap(matrix=mat, alpha=False)
        out_png = OUT / f"page_{i:02d}.png"
        pix.save(str(out_png))
        sz_kb = out_png.stat().st_size // 1024
        print(f"  page_{i:02d}.png  {pix.width}x{pix.height}  {sz_kb} KB")
    doc.close()
    print(f"\nDONE: {doc.page_count} PNGs in {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
