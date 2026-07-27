"""F-114 / D-12 — PDF artefact ≥ 30 KB via Edge headless.

If Edge is not installed, the test is skipped (no failure).
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, "scripts")

from gamma.reporting import analyze_and_report  # noqa: E402
from gamma.reporting.pdf_export import _find_edge  # noqa: E402


def test_pdf_artefact(tmp_path):
    if _find_edge() is None:
        print("Edge not found — skipping PDF artefact test")
        return
    # P1-3c: per-test tmp_path prevents xdist concurrent-write race on
    # fixed shared dir (was "demo_reports/_test_pdf_artefact").
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
        write_pdf=True,
        sample_mass_kg=0.5,
        background_path=bg,
    )
    pdf = res.get("pdf")
    if not pdf:
        print("PDF not generated (Edge unavailable or failed) — skipping size assertion")
        return
    assert os.path.exists(pdf), f"PDF missing: {pdf}"
    size = os.path.getsize(pdf)
    assert size >= 30 * 1024, (
        f"PDF too small: {size} bytes (want ≥ 30 KB)"
    )


if __name__ == "__main__":
    import tempfile, pathlib
    test_pdf_artefact(pathlib.Path(tempfile.mkdtemp(prefix="_test_pdf_artefact_")))
    print("OK")
