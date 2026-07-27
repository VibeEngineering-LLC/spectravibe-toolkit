"""Extract LibreOffice-converted ODT->DOCX into corpus."""
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from _corpus_extractor import extract_docx, OUT

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "references" / "_converted_tmp"
DST = OUT / "_converted_odt"

DST.mkdir(parents=True, exist_ok=True)
for f in sorted(SRC.glob("*.docx")):
    if "Vartanov" in f.name:
        continue
    if any(s in f.name for s in ["Прецизион", "Порядок"]):
        text, meta = extract_docx(f)
        out = DST / (f.stem + ".odt.md")
        out.write_text(
            f"<!-- src: odt:{f.name} -->\n<!-- meta: {meta} -->\n\n# {f.name}\n\n{text}",
            encoding="utf-8",
        )
        print(f"OK  {f.name}  {meta}")
