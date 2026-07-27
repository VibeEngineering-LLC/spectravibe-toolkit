"""F-115 (D-10) — anonymization regression.

Runs ``analyze_and_report`` on the Th-232 fixture and verifies that
no operator name, certified-source S/N, absolute path, or detector
S/N leaks into the JSON / Markdown / HTML artefacts.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, "scripts")

from gamma.reporting import analyze_and_report  # noqa: E402


FORBIDDEN = [
    "Нарчук",                  # operator name
    "SN-01",                 # detector S/N
    "420-7-17",                # source certificate S/N
    "D:\\",                    # absolute Windows path
    "D:/",                     # absolute Windows path (fwd slash)
    "УДС-ГЦ-63х63-USB №",      # detector type + serial token
]


def test_anonymization(tmp_path):
    # P1-3c: per-test tmp_path prevents xdist concurrent-write race on
    # fixed shared dir (was "demo_reports/_test_anonymize").
    out = str(tmp_path)
    sp = (
        "detectors/Gamma-1S/reference_spectra/"
        "archive/"
        "Th232_420-7-17_Маринелли_0cm.spe"
    )
    bg = (
        "detectors/Gamma-1S/data/averaged_backgrounds/"
        "bg_2016_marinelli_water_marinelli.spe"
    )
    res = analyze_and_report(
        sp,
        output_dir=out,
        write_html=True,
        write_markdown=True,
        write_plots=False,
        sample_mass_kg=0.5,
        background_path=bg,
    )

    for kind in ("json", "markdown", "html"):
        path = res.get(kind)
        assert path and os.path.exists(path), f"missing {kind}: {path}"
        text = open(path, encoding="utf-8").read()
        for token in FORBIDDEN:
            assert token not in text, (
                f"F-115 leak in {kind}: {token!r} appears in {path}"
            )


if __name__ == "__main__":
    import tempfile, pathlib
    test_anonymization(pathlib.Path(tempfile.mkdtemp(prefix="_test_anonymize_")))
    print("OK")
