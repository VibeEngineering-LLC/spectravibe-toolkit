from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""Second pass: extract the LibreOffice-converted .pptx/.docx into the same
_extracted_corpus tree (under a `_converted/` subfolder).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _corpus_extractor import extract_docx, extract_pptx, OUT  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "references" / "_converted_tmp"
DST = OUT / "_converted_legacy"


def main() -> int:
    DST.mkdir(parents=True, exist_ok=True)
    files = sorted(p for p in SRC.iterdir() if p.suffix.lower() in {".pptx", ".docx"})
    print(f"Found {len(files)} converted files")
    total_chars = 0
    for f in files:
        if f.suffix.lower() == ".pptx":
            text, meta = extract_pptx(f)
        else:
            text, meta = extract_docx(f)
        out = DST / (f.stem + f.suffix + ".md")
        out.write_text(
            f"<!-- src: legacy:{f.name} -->\n<!-- meta: {meta} -->\n\n# {f.name}\n\n{text}",
            encoding="utf-8",
        )
        total_chars += meta.get("n_chars", 0)
        print(f"OK  {f.name}  {meta}")
    print(f"\nDONE: {total_chars:,} chars")
    return 0


if __name__ == "__main__":
    sys.exit(main())
