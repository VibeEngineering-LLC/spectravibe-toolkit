"""F-110 (D-08, D-09) — 73-90 кэВ composite must be present when
Th-232 chain is dominant, AND diffuse zones (backscatter_region,
broad_compton_plateau) must be filtered out.
"""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, "scripts")

from gamma.reporting import analyze_and_report  # noqa: E402


def test_th_composite_present(tmp_path):
    # P1-3c: per-test tmp_path prevents xdist concurrent-write race on
    # fixed shared dir (was "demo_reports/_test_th_composite").
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
        write_html=False,
        write_markdown=False,
        write_plots=False,
        sample_mass_kg=0.5,
        background_path=bg,
    )

    report = json.load(open(res["json"], encoding="utf-8"))
    sp_list = report.get("secondary_peaks") or []

    has_cluster = any(
        (e.get("type") == "cluster"
         and e.get("energy_keV") is not None
         and 73.0 <= float(e["energy_keV"]) <= 90.0)
        for e in sp_list
    )
    assert has_cluster, (
        "F-110: 73-90 кэВ cluster missing from secondary_peaks "
        f"(found types: {[e.get('type') for e in sp_list]})"
    )

    for e in sp_list:
        t = e.get("type") or ""
        assert t not in ("backscatter_region", "broad_compton_plateau"), (
            f"F-110: diffuse zone {t!r} should be filtered out"
        )


if __name__ == "__main__":
    import tempfile, pathlib
    test_th_composite_present(pathlib.Path(tempfile.mkdtemp(prefix="_test_th_composite_")))
    print("OK")
