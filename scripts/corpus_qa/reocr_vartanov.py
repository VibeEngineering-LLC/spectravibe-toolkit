# -*- coding: utf-8 -*-
"""#LIB-1 Ф5 класс 4 — ре-OCR Vartanov.pdf через Tesseract rus+eng (замена мусорного djvu-OCR)."""
import fitz, os, sys, io
sys.stdout.reconfigure(encoding="utf-8")
import pytesseract
from PIL import Image
pytesseract.pytesseract.tesseract_cmd = r"C:\Program Files\Tesseract-OCR\tesseract.exe"

ROOT = r"<WORKDIR>\gamma-spectrum-analysis"
SRC = os.path.join(ROOT, r"books_library\Документация ЛСРМ\02_topical_pdf\Vartanov.pdf")
OUT = os.path.join(ROOT, r"references\_extracted_corpus\Документация ЛСРМ\02_topical_pdf\Vartanov_Prakticheskie_metody_scintillyatsionnoy_gamma-spektrometrii.reocr.md")
LOG = os.path.join(ROOT, r"audit\_drafts\_lib1\vartanov_reocr.log")

doc = fitz.open(SRC)
n = doc.page_count
parts = ["# Vartanov — Практические методы сцинтилляционной гамма-спектрометрии",
         "", "> Ре-OCR Tesseract rus+eng (#LIB-1 Ф5, 2026-07-12) — замена мусорного djvu-OCR.", ""]
with open(LOG, "w", encoding="utf-8") as lg:
    for i in range(n):
        pix = doc[i].get_pixmap(dpi=300)
        img = Image.open(io.BytesIO(pix.tobytes("png")))
        try:
            txt = pytesseract.image_to_string(img, lang="rus+eng")
        except Exception as e:
            txt = f"[OCR-ERROR: {e}]"
        parts.append(f"## PAGE {i+1}")
        parts.append(txt.strip())
        parts.append("")
        if (i+1) % 25 == 0:
            lg.write(f"{i+1}/{n} done\n"); lg.flush()
doc.close()
open(OUT, "w", encoding="utf-8").write("\n".join(parts))
with open(LOG, "a", encoding="utf-8") as lg:
    lg.write(f"DONE {n} pages -> {OUT}\n")
print(f"DONE {n} pages -> {OUT}")