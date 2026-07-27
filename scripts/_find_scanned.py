import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent / "references" / "_extracted_corpus"
report = json.load(open(ROOT / "_extraction_report.json", encoding="utf-8"))

scanned = [r for r in report["results"] if r.get("n_pages") and not r.get("has_text_layer") and "error" not in r]
print(f"{len(scanned)} fully scanned PDFs (no text layer):")
for r in scanned:
    print(f"  {r['n_pages']:3d} p  {r['src']}")

print()

sparse = [
    r for r in report["results"]
    if r.get("n_pages") and r.get("has_text_layer") and r.get("n_chars", 0) / max(1, r["n_pages"]) < 200
]
print(f"{len(sparse)} sparse-text PDFs (<200 chars/page — partial scan):")
for r in sparse:
    print(f"  {r['n_pages']:3d} p, {r['n_chars']:5d} ch  {r['src']}")
