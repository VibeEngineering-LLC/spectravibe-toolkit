"""F-108 (D-03) — table rows must be sorted by ascending energy.

Generates the Th-232 interactive HTML report and parses the JS
``const rows=[…]`` array.  Energies extracted from each row's
``line`` field must be monotonically non-decreasing.
"""
from __future__ import annotations

import json
import os
import re
import sys

sys.path.insert(0, "scripts")

from gamma.reporting import analyze_and_report  # noqa: E402


def test_rows_sorted_ascending(tmp_path):
    # P1-3c: per-test tmp_path prevents xdist concurrent-write race on
    # fixed shared dir (was "demo_reports/_test_sorted").
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
        write_markdown=False,
        write_plots=False,
        sample_mass_kg=0.5,
        background_path=bg,
    )

    html = open(res["html"], encoding="utf-8").read()

    # F-141 / v1.17.7 — HTML JS теперь может содержать символы ';'
    # внутри строк (например в `note`). Парсим через JSON.loads с
    # balanced bracket extraction.
    m = re.search(r"const\s+rows\s*=\s*(\[)", html)
    assert m, "could not locate `const rows = [...]` in HTML"
    start = m.end() - 1  # позиция '['
    depth, i = 0, start
    while i < len(html):
        c = html[i]
        if c == "[": depth += 1
        elif c == "]":
            depth -= 1
            if depth == 0:
                break
        i += 1
    rows = json.loads(html[start:i+1])
    assert rows, "rows array is empty"

    energies = []
    for r in rows:
        line = r.get("line") or ""
        em = re.search(r"\d+(?:\.\d+)?", line)
        if em:
            energies.append(float(em.group(0)))
    assert len(energies) >= 5, f"too few energies parsed: {energies}"

    for i in range(1, len(energies)):
        assert energies[i] >= energies[i-1] - 1e-6, (
            f"rows not sorted ascending at i={i}: "
            f"{energies[i-1]} > {energies[i]}"
        )


if __name__ == "__main__":
    import tempfile, pathlib
    test_rows_sorted_ascending(pathlib.Path(tempfile.mkdtemp(prefix="_test_sorted_")))
    print("OK")
