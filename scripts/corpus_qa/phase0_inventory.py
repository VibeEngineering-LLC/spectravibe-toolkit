# -*- coding: utf-8 -*-
"""Ф0 #LIB-1 — инвентаризация корпуса: реестр источник<->md, класс A/B/C/office/scan."""
import fitz, os, sys, json, re
sys.stdout.reconfigure(encoding="utf-8")

ROOT = r"<WORKDIR>\gamma-spectrum-analysis"
BOOKS = os.path.join(ROOT, "books_library")
CORPUS = os.path.join(ROOT, "references", "_extracted_corpus")
REFS = os.path.join(ROOT, "references")
OUT = os.path.join(ROOT, "audit", "_drafts", "_lib1", "manifest.json")

SCAN_PRODUCERS = ("finereader", "abbyy", "acrobat distiller", "scan", "turboscan", "camscanner", "piksoft")
DIGITAL_PRODUCERS = ("word", "latex", "tex", "indesign", "pdflatex", "dvips", "libreoffice", "openoffice", "writer", "powerpoint", "excel")
# born-digital office/text по расширению источника
OFFICE_EXT = (".pptx", ".ppt", ".docx", ".doc", ".odt", ".odp", ".xlsx", ".rtf")
SCAN_EXT = (".djvu",)

def classify_pdf(path):
    try:
        doc = fitz.open(path)
    except Exception as e:
        return {"error": str(e), "class": "pdf-err"}
    meta = doc.metadata or {}
    producer = (meta.get("producer") or "").lower()
    creator = (meta.get("creator") or "").lower()
    npages = doc.page_count
    idxs = sorted(set(int(i * (npages - 1) / 7) for i in range(8))) if npages > 1 else [0]
    total, sampled = 0, 0
    for i in idxs:
        try:
            total += len((doc[i].get_text() or "").strip()); sampled += 1
        except Exception:
            pass
    doc.close()
    per_page = total / sampled if sampled else 0
    prod = producer + " " + creator
    is_digital = any(k in prod for k in DIGITAL_PRODUCERS)
    is_scan = any(k in prod for k in SCAN_PRODUCERS)
    if per_page < 30:
        cls = "C"
    elif is_scan:
        cls = "B"
    elif is_digital or per_page > 200:
        cls = "A"
    else:
        cls = "B"
    return {"producer": producer, "creator": creator, "pages": npages,
            "layer_chars_per_page": round(per_page, 1), "class": cls}

# индекс всех источников в books_library по имени
src_index = {}
for r, d, f in os.walk(BOOKS):
    for fn in f:
        src_index.setdefault(fn.lower(), os.path.join(r, fn))

# все md/txt корпуса
mds = []
for r, d, f in os.walk(CORPUS):
    for fn in f:
        if fn.lower().endswith((".md", ".txt")):
            mds.append(os.path.join(r, fn))
for fn in os.listdir(REFS):
    fp = os.path.join(REFS, fn)
    if os.path.isfile(fp) and fn.lower().endswith(".md"):
        low = fn.lower()
        if any(k in low for k in ("ocr", "grundl", "green1995", "gulamrazul", "isma2014", "rjmcmc", "spe_format", "multiplet")):
            mds.append(fp)

records = []
src_used = set()
for md in mds:
    base = os.path.basename(md)
    # "<name>.<ext>.md" -> "<name>.<ext>"
    m = re.match(r"(.+\.[a-z0-9]+)\.md$", base, re.I)
    src_name = m.group(1).lower() if m else None
    cand = src_index.get(src_name) if src_name else None
    rec = {"md": os.path.relpath(md, ROOT), "size_bytes": os.path.getsize(md),
           "source": os.path.relpath(cand, ROOT) if cand else None}
    if cand:
        src_used.add(cand)
        ext = os.path.splitext(cand)[1].lower()
        if ext == ".pdf":
            rec.update(classify_pdf(cand))
        elif ext in OFFICE_EXT:
            rec["class"] = "office"     # born-digital, сверка нативной либой
            rec["src_ext"] = ext
        elif ext in SCAN_EXT:
            rec["class"] = "scan-djvu"  # скан, только рендер
            rec["src_ext"] = ext
        else:
            rec["class"] = "other"; rec["src_ext"] = ext
    else:
        rec["class"] = "no-source"
    records.append(rec)

orphan_src = [os.path.relpath(p, ROOT) for p in set(src_index.values()) if p not in src_used]

summary = {}
for rec in records:
    summary[rec.get("class", "?")] = summary.get(rec.get("class", "?"), 0) + 1

out = {"generated": "2026-07-12", "n_md": len(records), "n_source_files": len(set(src_index.values())),
       "class_summary": summary, "n_orphan_sources": len(orphan_src),
       "orphan_sources": sorted(orphan_src), "records": records}
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)

print(f"md: {len(records)}  источников в books_library: {len(set(src_index.values()))}")
print(f"классы: {summary}")
print(f"источники без извлечения (сироты): {len(orphan_src)}")