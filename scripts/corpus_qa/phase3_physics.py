# -*- coding: utf-8 -*-
"""Ф3 v2 #LIB-1 — физика: чистые границы изотопов + цепочечные линии + сужённое окно."""
import os, sys, json, re, csv, glob
from collections import defaultdict
sys.stdout.reconfigure(encoding="utf-8")
ROOT = r"<WORKDIR>\gamma-spectrum-analysis"
CORPUS = os.path.join(ROOT, "references", "_extracted_corpus")
REFS = os.path.join(ROOT, "references")
OUT = os.path.join(ROOT, "audit", "_drafts", "_lib1", "findings_physics.json")

ELEMENTS = set("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr".split())

valid = set()
energies = defaultdict(set)
chain_of = {}
nj = json.load(open(os.path.join(ROOT, "data", "nuclides.json"), encoding="utf-8"))
for k, v in nj.items():
    if k == "_schema" or not isinstance(v, dict): continue
    valid.add(k)
    if v.get("chain"): chain_of[k] = v["chain"]
    for ln in v.get("gamma_lines", []) or []:
        e = ln.get("energy_keV") or ln.get("energy")
        if isinstance(e, (int, float)): energies[k].add(round(float(e), 2))
for f in glob.glob(os.path.join(REFS, "iaea_cache", "*_g.csv")):
    m = re.match(r"(\d+)([a-z]+)(m?)$", os.path.basename(f)[:-6])
    if not m: continue
    sym = m.group(2).capitalize()
    if sym not in ELEMENTS: continue
    iso = f"{sym}-{m.group(1)}"; valid.add(iso)
    try:
        for row in csv.DictReader(open(f, encoding="utf-8")):
            try:
                e = float(row["energy"]); inten = float(row.get("intensity") or 0)
                if inten >= 0.3: energies[iso].add(round(e, 2))
            except (ValueError, KeyError, TypeError): pass
    except Exception: pass

# цепочечные линии: изотоп в равновесии "видит" линии всей своей цепочки
chain_lines = defaultdict(set)
for iso, ch in chain_of.items():
    chain_lines[ch] |= energies.get(iso, set())
def lib_lines(iso):
    s = set(energies.get(iso, set()))
    if iso in chain_of: s |= chain_lines[chain_of[iso]]
    return s

# чистые границы: перед символом — не буква/цифра/дефис; после массового — не буква
ISO_RE = re.compile(r"(?<![A-Za-z0-9\-])([A-Z][a-z]?)-(\d{1,3})m?(?![A-Za-z])|\^\{?(\d{1,3})\}?\s*([A-Z][a-z]?)(?![a-z])")
ENERGY_RE = re.compile(r"(\d{2,4}(?:[.,]\d{1,2})?)\s*(?:кэВ|keV|кэв)")
STOP = re.compile(r"[,;.]|\bи\b|\band\b|[A-Z][a-z]?-\d")  # обрыв окна изотоп->энергия

def find_isos(text):
    out = []
    for m in re.finditer(ISO_RE, text):
        g = m.groups()
        sym, A = (g[0], g[1]) if g[0] else (g[3], g[2])
        if sym in ELEMENTS and 1 <= int(A) <= 300:
            out.append((f"{sym}-{int(A)}", m.start(), m.end()))
    return out

targets = []
for r, d, fs in os.walk(CORPUS):
    for fn in fs:
        if fn.lower().endswith((".md", ".txt")): targets.append(os.path.join(r, fn))
for fn in os.listdir(REFS):
    fp = os.path.join(REFS, fn)
    if os.path.isfile(fp) and fn.lower().endswith(".md"): targets.append(fp)

nonexist = defaultdict(list); mismatch = []
for t in targets:
    try: text = open(t, encoding="utf-8", errors="replace").read()
    except Exception: continue
    rel = os.path.relpath(t, ROOT)
    for iso, s, e in find_isos(text):
        if iso not in valid:
            if len(nonexist[iso]) < 4:
                nonexist[iso].append({"file": rel, "ctx": text[max(0,s-25):e+15].replace("\n"," ")})
            continue
        lib = lib_lines(iso)
        if not lib: continue
        win = text[e:e+30]
        stop = STOP.search(win)
        win = win[:stop.start()] if stop else win   # обрыв на разделителе/новом изотопе
        em = ENERGY_RE.search(win)
        if em:
            E = float(em.group(1).replace(",", "."))
            if 20 <= E <= 3000:
                nearest = min(lib, key=lambda x: abs(x-E))
                if abs(nearest-E) > 2.0:
                    mismatch.append({"file": rel, "iso": iso, "text_E": E, "nearest": round(nearest,2),
                                     "delta": round(E-nearest,2), "ctx": text[max(0,s-10):e+30].replace("\n"," ")})

json.dump({"generated":"2026-07-12","n_valid":len(valid),
           "nonexistent":{k:v for k,v in sorted(nonexist.items(),key=lambda x:-len(x[1]))},
           "mismatches":mismatch}, open(OUT,"w",encoding="utf-8"), ensure_ascii=False, indent=2)

print(f"валидных изотопов: {len(valid)}   вне баз: {len(nonexist)}   энерг.расхожд: {len(mismatch)}")
print()
print("=== изотопы вне баз (после чистки границ) ===")
for iso, occ in sorted(nonexist.items(), key=lambda x:-len(x[1]))[:25]:
    print(f"  {iso:8} x{len(occ)}  «{occ[0]['ctx'].strip()[:50]}» [{os.path.basename(occ[0]['file'])[:26]}]")
print()
print("=== энергия != библиотека (с цепочечными линиями, обрыв окна) ===")
for m in sorted(mismatch, key=lambda x:-abs(x['delta'])):
    print(f"  {m['iso']:8} txt={m['text_E']:7} биб={m['nearest']:7} Δ={m['delta']:+8}  «{m['ctx'].strip()[:44]}» [{os.path.basename(m['file'])[:22]}]")