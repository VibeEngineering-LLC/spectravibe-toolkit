# -*- coding: utf-8 -*-
"""Ф1 агрегатор: разбить S1 по Unicode-блокам, S2 по типу гомоглифа."""
import json, os, sys, unicodedata
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8")

ROOT = r"<WORKDIR>\gamma-spectrum-analysis"
F = os.path.join(ROOT, "audit", "_drafts", "_lib1", "findings_unicode.json")
data = json.load(open(F, encoding="utf-8"))

def block(cp):
    o = int(cp[2:], 16)
    if 0x4E00 <= o <= 0x9FFF: return "CJK-иероглифы"
    if 0x3400 <= o <= 0x4DBF: return "CJK-ext-A"
    if 0xAC00 <= o <= 0xD7AF: return "хангыль"
    if 0x0600 <= o <= 0x06FF: return "арабский"
    if 0x0590 <= o <= 0x05FF: return "иврит"
    if 0x3040 <= o <= 0x30FF: return "кана"
    if 0xE000 <= o <= 0xF8FF: return "PUA(office-шрифты)"
    if 0xF900 <= o <= 0xFAFF: return "CJK-compat"
    if 0xFB00 <= o <= 0xFB4F: return "лат/ивр-лигатуры"
    if 0xFF00 <= o <= 0xFFEF: return "fullwidth-формы"
    return f"прочее({cp})"

s1_blocks = Counter()
s1_by_file = Counter()
s1_examples = {}
s2_type = Counter()
s2_by_file = Counter()
s2_examples = []
for r in data["results"]:
    fn = os.path.basename(r["file"])
    for it in r.get("findings", {}).get("S1", []):
        b = block(it["cp"])
        s1_blocks[b] += 1
        s1_by_file[(b, fn)] += 1
        if b not in s1_examples:
            s1_examples[b] = f"{it['char']} ({it['cp']} {it.get('name','')}) в «{it['ctx'].strip()}»"
    for it in r.get("findings", {}).get("S2", []):
        w = it["word"]
        # тип: cyr-dominant с латинской вставкой или наоборот
        ncyr = sum(1 for c in w if 0x0400 <= ord(c) <= 0x04FF)
        nlat = sum(1 for c in w if c.isascii() and c.isalpha())
        kind = "кир-слово+лат-буквы" if ncyr >= nlat else "лат-слово+кир-буквы"
        s2_type[kind] += 1
        s2_by_file[fn] += 1
        if len(s2_examples) < 25:
            s2_examples.append(f"{w}  [{fn}]")

print("=== S1: чужие скрипты по блокам Unicode ===")
for b, n in s1_blocks.most_common():
    print(f"  {b:28} {n:5}   пример: {s1_examples.get(b,'')[:80]}")
print()
print("=== S1: топ файлов с не-PUA находками ===")
for (b, fn), n in s1_by_file.most_common():
    if not b.startswith("PUA"):
        print(f"  {n:4}  {b:20} {fn[:60]}")
print()
print("=== S2: типы гомоглифов ===")
for k, n in s2_type.most_common():
    print(f"  {k:24} {n}")
print()
print("=== S2: топ файлов ===")
for fn, n in s2_by_file.most_common(12):
    print(f"  {n:4}  {fn[:65]}")
print()
print("=== S2 примеры (первые 25) ===")
for e in s2_examples:
    print("  " + e)