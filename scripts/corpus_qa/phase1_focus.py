# -*- coding: utf-8 -*-
"""Ф1 фокус: CJK по файлам + гомоглифы уникальные + гомоглифы в критичных файлах."""
import json, os, sys
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8")
ROOT = r"<WORKDIR>\gamma-spectrum-analysis"
data = json.load(open(os.path.join(ROOT, "audit", "_drafts", "_lib1", "findings_unicode.json"), encoding="utf-8"))

def is_cjk(cp):
    o = int(cp[2:], 16); return 0x3400 <= o <= 0x9FFF or 0xF900 <= o <= 0xFAFF

# критичные файлы (идут в расчёт/идентификацию)
CRIT = ("nuclide", "physics-contracts", "identification", "xrf", "glossary", "metadata_calibration",
        "peak_search", "STEPS", "reports-procedure", "lsrm", "dose", "spe_format", "aspect")

print("=== CJK по файлам ===")
for r in data["results"]:
    cjk = [it for it in r.get("findings", {}).get("S1", []) if is_cjk(it["cp"])]
    if cjk:
        chars = " ".join(f"{it['char']}@{it['pos']}" for it in cjk)
        print(f"  {os.path.basename(r['file'])[:55]:55} {len(cjk)}  {chars}")

print()
print("=== S2 гомоглифы: уникальные слова с частотой (весь корпус, топ-40) ===")
words = Counter()
for r in data["results"]:
    for it in r.get("findings", {}).get("S2", []):
        words[it["word"]] += 1
for w, n in words.most_common(40):
    print(f"  {n:5}  {w}")

print()
print("=== S2 в критичных файлах (не учебники, не JSON-индексы) ===")
for r in data["results"]:
    fn = os.path.basename(r["file"]).lower()
    if any(c in fn for c in CRIT) and not fn.endswith(".json"):
        s2 = r.get("findings", {}).get("S2", [])
        if s2:
            uniq = sorted(set(it["word"] for it in s2))
            print(f"  {os.path.basename(r['file'])[:50]:50} {len(s2)}: {', '.join(uniq[:20])}")