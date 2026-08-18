"""F-108 (D-19) — cost footer regression."""
from __future__ import annotations

import os
import sys

sys.path.insert(0, "scripts")

from gamma.reporting import analyze_and_report  # noqa: E402


def test_cost_footer(tmp_path):
    # P1-3c: per-test tmp_path prevents xdist concurrent-write race on
    # fixed shared dir (was "demo_reports/_test_cost_footer").
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
        cost_estimate={
            "tokens": 140000,
            "session_pct": "65%",
            "detail": "regression",
        },
    )
    html = open(res["html"], encoding="utf-8").read()
    assert "Стоимость анализа" in html, "cost footer label missing"
    assert "65%" in html, "session_pct missing in footer"
    # Token count should appear in the footer (we render it as ~140000)
    assert "140000" in html or "140 000" in html, (
        "tokens value missing in footer"
    )


if __name__ == "__main__":
    import tempfile, pathlib
    test_cost_footer(pathlib.Path(tempfile.mkdtemp(prefix="_test_cost_footer_")))
    print("OK")
