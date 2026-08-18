"""Extract embedded <BackgroundEnergySpectrum> from a RadiaCode/AtomSpectra XML
into a standalone BecqMoni/AtomSpectra XML file.

Зачем: run_plan_a.py требует отдельный GAMMA_BG-файл, а RadiaCode пишет фон
встроенным блоком <BackgroundEnergySpectrum> внутри файла образца. Этот хелпер
вытаскивает встроенный фон в самостоятельный файл, пригодный как --background /
GAMMA_BG для канонического pipeline.

Usage:
    python scripts/extract_embedded_bg.py <sample.xml> <out_bg.xml>

Exit codes: 0 OK, 1 нет встроенного фона, 2 неверные аргументы.
"""
from __future__ import annotations
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from gamma.io.atomspectra_xml import read_atomspectra_xml
from gamma.io.becqmoni_xml import write_becqmoni_xml


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: extract_embedded_bg.py <sample.xml> <out_bg.xml>")
        return 2
    src, out = sys.argv[1], sys.argv[2]
    spec = read_atomspectra_xml(src)
    bg = spec.background_embedded
    if bg is None:
        print(f"NO embedded background in {src}")
        return 1
    write_becqmoni_xml(bg, out)
    import numpy as np
    print(f"OK: embedded bg -> {out}")
    print(f"   channels={len(bg.counts)} live_time={bg.live_time:.1f}s "
          f"real_time={bg.real_time:.1f}s sum_counts={int(np.asarray(bg.counts).sum())}")
    return 0


if __name__ == "__main__":
    sys.exit(main())