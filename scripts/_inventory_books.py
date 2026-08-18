import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parent.parent / "references" / "books"
by_ext = defaultdict(list)
for p in ROOT.rglob("*"):
    if p.is_file():
        by_ext[p.suffix.lower() or "<no-ext>"].append(p.relative_to(ROOT).as_posix())

for ext in sorted(by_ext):
    print(f"{ext}: {len(by_ext[ext])}")

print()
print("Unusual extensions:")
for ext in sorted(by_ext):
    if ext not in {".pdf", ".docx", ".pptx", ".doc", ".ppt"}:
        for f in by_ext[ext][:5]:
            print(f"  {ext}: {f}")
