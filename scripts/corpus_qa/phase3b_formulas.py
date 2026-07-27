# -*- coding: utf-8 -*-
"""Ф3b #LIB-1 — формулы: синтаксис (баланс/CJK внутри math) + парсируемость (latexwalker)."""
import os, sys, json, re
sys.stdout.reconfigure(encoding="utf-8")
from pylatexenc.latexwalker import LatexWalker, LatexWalkerParseError
ROOT = r"<WORKDIR>\gamma-spectrum-analysis"
CORPUS = os.path.join(ROOT, "references", "_extracted_corpus")
REFS = os.path.join(ROOT, "references")
OUT = os.path.join(ROOT, "audit", "_drafts", "_lib1", "findings_formulas.json")

# извлечь формулы: \( .. \), \[ .. \], $$ .. $$, $ .. $
FORMULA_RE = re.compile(r"\\\((.+?)\\\)|\\\[(.+?)\\\]|\$\$(.+?)\$\$|(?<!\$)\$(?!\$)(.+?)(?<!\$)\$(?!\$)", re.S)
FOREIGN = re.compile(r"[㐀-鿿가-힯֐-׿؀-ۿ]")  # CJK/хангыль/ивр/араб

targets = []
for r, d, fs in os.walk(CORPUS):
    for fn in fs:
        if fn.lower().endswith((".md", ".txt")): targets.append(os.path.join(r, fn))
for fn in os.listdir(REFS):
    fp = os.path.join(REFS, fn)
    if os.path.isfile(fp) and fn.lower().endswith(".md"): targets.append(fp)

findings = []
n_formulas = 0
for t in targets:
    try: text = open(t, encoding="utf-8", errors="replace").read()
    except Exception: continue
    rel = os.path.relpath(t, ROOT)
    for m in FORMULA_RE.finditer(text):
        body = next((g for g in m.groups() if g is not None), "")
        if not body.strip(): continue
        n_formulas += 1
        issues = []
        if body.count("{") != body.count("}"): issues.append(f"скобки {{}}: {body.count('{')}/{body.count('}')}")
        fr = FOREIGN.search(body)
        if fr: issues.append(f"чужой символ '{fr.group()}' U+{ord(fr.group()):04X}")
        try:
            LatexWalker(body).get_latex_nodes()
        except LatexWalkerParseError as e:
            issues.append(f"parse: {str(e)[:50]}")
        except Exception:
            pass
        if issues:
            findings.append({"file": rel, "pos": m.start(), "formula": body.strip()[:80], "issues": issues})

json.dump({"generated":"2026-07-12","n_formulas":n_formulas,"n_flagged":len(findings),
           "findings":findings}, open(OUT,"w",encoding="utf-8"), ensure_ascii=False, indent=2)
print(f"формул извлечено: {n_formulas}   с проблемами: {len(findings)}")
print()
# группировка по типу
from collections import Counter
by_type = Counter()
for f in findings:
    for i in f["issues"]:
        by_type[i.split(":")[0]] += 1
print("=== типы проблем ===")
for k, n in by_type.most_common(): print(f"  {k:20} {n}")
print()
print("=== формулы с чужими символами (порча) ===")
for f in findings:
    if any("чужой" in i for i in f["issues"]):
        print(f"  [{os.path.basename(f['file'])[:30]}] {f['issues']}  «{f['formula'][:50]}»")
print()
print("=== топ-15 формул с дисбалансом скобок ===")
cnt = 0
for f in findings:
    if any("скобки" in i for i in f["issues"]) and cnt < 15:
        print(f"  [{os.path.basename(f['file'])[:28]}] {[i for i in f['issues'] if 'скобки' in i]}  «{f['formula'][:45]}»")
        cnt += 1