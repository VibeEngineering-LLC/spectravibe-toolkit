"""D-01 — environment classifier must NOT generate sample-narrative
text when the analysed file is a pure background spectrum.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, "scripts")

from gamma.reporting import analyze_and_report  # noqa: E402


FORBIDDEN_SUBSTRINGS = [
    "в образце",
    "в пробе",
    "источник содержит",
    "sample contains",
]


def test_bg_only_environment(tmp_path):
    # P1-3c: per-test tmp_path prevents xdist concurrent-write race on
    # fixed shared dir (was "demo_reports/_test_bg_only").
    out = str(tmp_path)
    bg = (
        "detectors/Gamma-1S/data/averaged_backgrounds/"
        "bg_2016_marinelli_water_marinelli.spe"
    )
    assert os.path.exists(bg), f"fixture missing: {bg}"

    res = analyze_and_report(
        bg,
        output_dir=out,
        write_html=True,
        write_markdown=True,
        write_plots=False,
    )

    html = open(res["html"], encoding="utf-8").read()
    md = open(res["markdown"], encoding="utf-8").read()

    for text, where in ((html, "HTML"), (md, "Markdown")):
        lower = text.lower()
        for needle in FORBIDDEN_SUBSTRINGS:
            assert needle.lower() not in lower, (
                f"D-01: forbidden phrase {needle!r} present in {where} "
                f"for pure-background file"
            )


if __name__ == "__main__":
    import tempfile, pathlib
    test_bg_only_environment(pathlib.Path(tempfile.mkdtemp(prefix="_test_bg_only_")))
    print("OK")
