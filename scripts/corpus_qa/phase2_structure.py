# -*- coding: utf-8 -*-
"""Ф2 #LIB-1 — структура: PAGE-маркеры vs page_count, плотность, баланс LaTeX.
Ф2b — для класса A: word-set diff md <-> text-layer (ловит OCR-искажения)."""
import fitz, os, sys, json, re
sys.stdout.reconfigure(encoding="utf-8")
ROOT = r"<WORKDIR>\gamma-spectrum-analysis"
MANIFEST = os.path.join(ROOT, "audit", "_drafts", "_lib1", "manifest.json")
OUT = os.path.join(ROOT, "audit", "_drafts", "_lib1", "findings_structure.json")
man = json.load(open(MANIFEST, encoding="utf-8"))

def latex_balance(text):
    issues = []
    if text.count(r"\(") != text.count(r"\)"): issues.append(f"inline \\(={text.count(chr(92)+'(')} \\)={text.count(chr(92)+')')}")
    if text.count(r"\[") != text.count(r"\]"): issues.append(f"display \\[={text.count(chr(92)+'[')} \\]={text.count(chr(92)+']')}")
    if text.count("$$") % 2: issues.append("нечётное $$")
    return issues

def page_markers(text):
    return [int(m.group(1)) for m in re.finditer(r"(?im)^#+\s*PAGE\s+(\d+)", text)]

WORD = re.compile(r"[A-Za-zА-Яа-яЁё]{4,}")
def words(text):
    return set(w.lower() for w in WORD.findall(text))

results = []
for rec in man["records"]:
    cls = rec.get("class")
    md_path = os.path.join(ROOT, rec["md"])
    try:
        md = open(md_path, encoding="utf-8", errors="replace").read()
    except Exception as e:
        continue
    entry = {"md": rec["md"], "class": cls}
    # структура
    pm = page_markers(md)
    if pm:
        non_mono = sum(1 for a, b in zip(pm, pm[1:]) if b <= a)
        entry["page_markers"] = {"n": len(pm), "min": min(pm), "max": max(pm), "non_monotonic": non_mono}
    lb = latex_balance(md)
    if lb: entry["latex_imbalance"] = lb
    # vs PDF page_count + Ф2b diff (класс A/B/C с PDF)
    src = rec.get("source")
    if src and src.lower().endswith(".pdf"):
        try:
            doc = fitz.open(os.path.join(ROOT, src))
            npages = doc.page_count
            if pm:
                cov = len(set(pm)) / npages if npages else 0
                entry["page_coverage"] = {"pdf_pages": npages, "md_pages": len(set(pm)), "coverage": round(cov, 3)}
            # Ф2b: word diff только для класса A (born-digital, layer надёжен)
            if cls == "A":
                layer = "\n".join(doc[i].get_text() for i in range(npages))
                w_md, w_pdf = words(md), words(layer)
                only_md = w_md - w_pdf
                # подозрительные: содержат не-кир/лат (порча) или очень длинные
                susp = [w for w in only_md if any(ord(c) > 0x04FF for c in w) or len(w) > 22]
                entry["diff_A"] = {"words_md": len(w_md), "words_pdf": len(w_pdf),
                                   "only_in_md": len(only_md), "suspicious": sorted(susp)[:30]}
            doc.close()
        except Exception as e:
            entry["pdf_error"] = str(e)
    # копим только если есть что показать
    flag = ("latex_imbalance" in entry or
            (entry.get("page_coverage", {}).get("coverage", 1) < 0.9) or
            (entry.get("page_markers", {}).get("non_monotonic", 0) > 0) or
            (entry.get("diff_A", {}).get("suspicious")))
    entry["_flagged"] = bool(flag)
    results.append(entry)

flagged = [r for r in results if r["_flagged"]]
json.dump({"generated": "2026-07-12", "n_records": len(results), "n_flagged": len(flagged),
           "results": results}, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=2)

print(f"обработано записей: {len(results)}  с флагами: {len(flagged)}")
print()
print("=== LaTeX-дисбаланс ===")
for r in results:
    if "latex_imbalance" in r:
        print(f"  {os.path.basename(r['md'])[:55]:55} {r['latex_imbalance']}")
print()
print("=== покрытие страниц < 0.9 (класс A/B/C с PDF) ===")
for r in results:
    pc = r.get("page_coverage")
    if pc and pc["coverage"] < 0.9:
        print(f"  {os.path.basename(r['md'])[:50]:50} {pc['md_pages']}/{pc['pdf_pages']} = {pc['coverage']}")
print()
print("=== немонотонные PAGE-маркеры ===")
for r in results:
    nm = r.get("page_markers", {}).get("non_monotonic", 0)
    if nm: print(f"  {os.path.basename(r['md'])[:55]:55} {nm} обратных переходов")
print()
print("=== Ф2b: класс A — подозрительные слова (в md, нет в text-layer) ===")
for r in results:
    d = r.get("diff_A")
    if d and d["suspicious"]:
        print(f"  {os.path.basename(r['md'])[:50]:50} only_md={d['only_in_md']}: {d['suspicious']}")