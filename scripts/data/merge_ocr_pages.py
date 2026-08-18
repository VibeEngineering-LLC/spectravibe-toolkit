"""
Объединяет page_NN/result.md из папки OCR-вывода в единый .md файл.

Использование:
    python merge_ocr_pages.py <ocr_dir> <out_file> [--title "..."]
"""
import sys
import argparse
from pathlib import Path


def merge(ocr_dir: Path, out_file: Path, title: str):
    page_dirs = sorted(
        [d for d in ocr_dir.iterdir() if d.is_dir() and d.name.startswith("page_")],
        key=lambda d: int(d.name.split("_")[1])
    )
    if not page_dirs:
        print(f"ERROR: нет page_NN/ в {ocr_dir}", file=sys.stderr)
        sys.exit(1)

    out_file.parent.mkdir(parents=True, exist_ok=True)
    parts = []
    parts.append(f"<!-- src: {ocr_dir.name} -->")
    parts.append(f"<!-- ext: .pdf (OCR via baidu/Unlimited-OCR VLM, 2026-07-04) -->")
    parts.append(f"<!-- meta: pages={len(page_dirs)}, dpi=200 -->")
    parts.append("")
    parts.append(f"# {title}")
    parts.append("")

    for pd in page_dirs:
        page_num = int(pd.name.split("_")[1])
        result = pd / "result.md"
        if not result.exists():
            parts.append(f"## PAGE {page_num}")
            parts.append("")
            parts.append("_(страница не распознана)_")
            parts.append("")
            continue
        content = result.read_text(encoding="utf-8").strip()
        parts.append(f"## PAGE {page_num}")
        parts.append("")
        parts.append(content)
        parts.append("")

    out_file.write_text("\n".join(parts), encoding="utf-8")
    print(f"OK: {len(page_dirs)} стр. → {out_file}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ocr_dir")
    ap.add_argument("out_file")
    ap.add_argument("--title", default="")
    args = ap.parse_args()
    merge(Path(args.ocr_dir), Path(args.out_file), args.title)


if __name__ == "__main__":
    main()
