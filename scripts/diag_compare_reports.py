"""Diagnostic: deep-compare two report.json files (A vs B).
Usage: python scripts/diag_compare_reports.py <A.json> <B.json>
Prints top-level keys of A, total diff count, and first 50 differing paths.
Float tolerance 1e-9. Reusable for any two-run comparison (e.g. FWHM sigma-mode).
"""
import json, sys

def load(p):
    with open(p, encoding="utf-8") as f:
        return json.load(f)

def walk(a, b, path=""):
    d = []
    if type(a) is not type(b):
        return [(path, f"type {type(a).__name__} vs {type(b).__name__}")]
    if isinstance(a, dict):
        for k in sorted(set(a) | set(b), key=str):
            if k not in a:
                d.append((path + "/" + str(k), "only in B"))
            elif k not in b:
                d.append((path + "/" + str(k), "only in A"))
            else:
                d += walk(a[k], b[k], path + "/" + str(k))
    elif isinstance(a, list):
        if len(a) != len(b):
            d.append((path, f"len {len(a)} vs {len(b)}"))
        for i, (x, y) in enumerate(zip(a, b)):
            d += walk(x, y, f"{path}[{i}]")
    else:
        if isinstance(a, (int, float)) and isinstance(b, (int, float)):
            if abs(float(a) - float(b)) > 1e-9:
                d.append((path, f"{a} vs {b}"))
        elif a != b:
            d.append((path, f"{a!r} vs {b!r}"))
    return d

A = load(sys.argv[1])
B = load(sys.argv[2])
print("A top-level keys:", list(A.keys()))
diffs = walk(A, B)
print(f"TOTAL DIFFS: {len(diffs)}")
for p, v in diffs[:50]:
    print("  ", p, "=>", v)