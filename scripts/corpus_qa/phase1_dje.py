# -*- coding: utf-8 -*-
import json, os, sys
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8")
ROOT = r"<WORKDIR>\gamma-spectrum-analysis"
data = json.load(open(os.path.join(ROOT,"audit","_drafts","_lib1","findings_unicode.json"), encoding="utf-8"))
# Џ (U+040F) — где, в каких словах
dje_files = Counter()
dje_words = Counter()
for r in data["results"]:
    for it in r.get("findings",{}).get("S2",[]):
        if "Џ" in it["word"]:
            dje_files[os.path.basename(r["file"])]+=1
            dje_words[it["word"]]+=1
print("=== Џ (U+040F) в словах ===")
for w,n in dje_words.most_common(20): print(f"  {n:4}  {w}")
print("--- файлы ---")
for f,n in dje_files.most_common(): print(f"  {n:4}  {f[:60]}")
# новый CJK в rjmcmc_synthesis — контекст
print()
print("=== CJK контексты в источниках ===")
for r in data["results"]:
    fn=os.path.basename(r["file"])
    if fn in ("rjmcmc_synthesis_green_gulamrazul_isma2014.md",):
        for it in r.get("findings",{}).get("S1",[]):
            o=int(it["cp"][2:],16)
            if 0x3400<=o<=0x9FFF:
                print(f"  {fn}: {it['char']} ({it['cp']}) ctx: ...{it['ctx'].strip()}...")