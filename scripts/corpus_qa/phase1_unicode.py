# -*- coding: utf-8 -*-
"""Ф1 #LIB-1 — Unicode-скан корпуса. S1 чужие скрипты, S2 гомоглифы, S3 mojibake."""
import os, sys, json, re, unicodedata
sys.stdout.reconfigure(encoding="utf-8")

ROOT = r"<WORKDIR>\gamma-spectrum-analysis"
CORPUS = os.path.join(ROOT, "references", "_extracted_corpus")
REFS = os.path.join(ROOT, "references")
RAG = os.path.join(ROOT, "audit", "_rag")
OUT = os.path.join(ROOT, "audit", "_drafts", "_lib1", "findings_unicode.json")

# --- разрешённые скрипты/блоки
def allowed(ch):
    o = ord(ch)
    if o < 0x0250: return True                 # ASCII + Latin-1 + Latin Extended-A
    if 0x0250 <= o <= 0x02FF: return True       # IPA + spacing modifiers
    if 0x0300 <= o <= 0x036F: return True       # комбинируемая диакритика
    if 0x0370 <= o <= 0x03FF: return True       # греческий
    if 0x0400 <= o <= 0x04FF: return True       # кириллица
    if 0x2000 <= o <= 0x206F: return True       # пунктуация (тире, кавычки)
    if 0x2070 <= o <= 0x209F: return True       # суб/суперскрипты
    if 0x20A0 <= o <= 0x20CF: return True       # валютные
    if 0x2100 <= o <= 0x214F: return True       # letterlike (№, ℮)
    if 0x2190 <= o <= 0x21FF: return True       # стрелки
    if 0x2200 <= o <= 0x22FF: return True       # мат.операторы
    if 0x2300 <= o <= 0x23FF: return True       # тех.символы
    if 0x25A0 <= o <= 0x25FF: return True       # геом.фигуры
    if 0x2600 <= o <= 0x26FF: return True       # разное
    if 0x2A00 <= o <= 0x2AFF: return True       # доп.мат.операторы
    if o in (0x00B0, 0x00B5, 0x00B1, 0x00D7, 0x00F7): return True  # ° µ ± × ÷
    return False

CYR = lambda c: 0x0400 <= ord(c) <= 0x04FF
LAT = lambda c: (0x41 <= ord(c) <= 0x5A) or (0x61 <= ord(c) <= 0x7A)

def ctx(text, i, w=40):
    return text[max(0, i - w):i + w + 1].replace("\n", " ")

def scan_text(text):
    f = {"S1": [], "S2": [], "S3": []}
    # S1: недопустимые символы
    for i, ch in enumerate(text):
        if not allowed(ch) and not ch.isspace():
            try:
                name = unicodedata.name(ch)
            except ValueError:
                name = "?"
            f["S1"].append({"pos": i, "char": ch, "cp": f"U+{ord(ch):04X}", "name": name, "ctx": ctx(text, i)})
    # S2: гомоглифы — смешение кир/лат в одном слове (буквы, без пробелов/цифр)
    for m in re.finditer(r"[A-Za-zЀ-ӿ]{2,}", text):
        w = m.group()
        has_cyr = any(CYR(c) for c in w)
        has_lat = any(LAT(c) for c in w)
        if has_cyr and has_lat:
            f["S2"].append({"pos": m.start(), "word": w, "ctx": ctx(text, m.start())})
    # S3: mojibake
    for i, ch in enumerate(text):
        if ch == "�":
            f["S3"].append({"pos": i, "cp": "U+FFFD", "ctx": ctx(text, i)})
    return f

targets = []
for base in (CORPUS,):
    for r, d, fs in os.walk(base):
        for fn in fs:
            if fn.lower().endswith((".md", ".txt")):
                targets.append(os.path.join(r, fn))
for fn in os.listdir(REFS):
    fp = os.path.join(REFS, fn)
    if os.path.isfile(fp) and fn.lower().endswith((".md", ".json")):
        targets.append(fp)
for r, d, fs in os.walk(RAG):
    for fn in fs:
        if fn.lower().endswith(".json"):
            targets.append(os.path.join(r, fn))

results = []
tot = {"S1": 0, "S2": 0, "S3": 0}
for t in targets:
    try:
        with open(t, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
    except Exception as e:
        results.append({"file": os.path.relpath(t, ROOT), "read_error": str(e)}); continue
    f = scan_text(text)
    n = {k: len(v) for k, v in f.items()}
    if any(n.values()):
        results.append({"file": os.path.relpath(t, ROOT), "counts": n, "findings": f})
        for k in tot: tot[k] += n[k]

out = {"generated": "2026-07-12", "n_files_scanned": len(targets),
       "n_files_with_findings": len(results), "totals": tot, "results": results}
with open(OUT, "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)

print(f"файлов просканировано: {len(targets)}")
print(f"файлов с находками: {len(results)}")
print(f"итого S1(чужой скрипт)={tot['S1']}  S2(гомоглиф)={tot['S2']}  S3(mojibake)={tot['S3']}")
print(f"-> {OUT}")