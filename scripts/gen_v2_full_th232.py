# -*- coding: utf-8 -*-
from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""F-367 / v1.18.24.2 — Th-232 V2-full demo generator.

Запускает production `analyze_and_report` дважды на одном Th-232 датасете:
1. Через `gamma.experimental.v2_integration.analyze_and_report_v2`
   (peak search = V2 dual-method, всё остальное — production)
2. Output → demo_reports/v1_18_24_th232_full/sample_v2/

Результат содержит ТОТ ЖЕ комплект артефактов что и production sample:
* `Th232_*_report.json`
* `Th232_*_report.md`
* `Th232_*_report.html` (Chart.js spectrum + multiplet PNGs + tables)
* `Th232_*_technical_report.pdf`
* `Th232_*_plots/spectrum.png`
* `Th232_*_plots/multiplets/multiplet_*.png`
* `Th232_*_calibrated.bq.xml` + bg.bq.xml (BecqMoni)

Отличается от production sample отчёта ТОЛЬКО peak search-методом —
V2 dual-method даёт обычно +N пиков (типично 11→16 для Th-232).
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from gamma.experimental.v2_integration import analyze_and_report_v2


def main():
    KIT = (
        REPO / "detectors" / "Gamma-1S" / "reference_spectra"
        / "reference_kits" / "Marinelli_1L" / "Th-232"
    )
    sample = KIT / "Th232_420-7-17_Маринелли_0cm.spe"
    bg = KIT / "Фон закр кр вода_13.spe"
    output_dir = REPO / "demo_reports" / "v1_18_24_th232_full" / "sample_v2"
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[V2-full] analyze_and_report_v2 на {sample.name}")
    print(f"[V2-full] bg={bg.name}")
    print(f"[V2-full] output={output_dir}")

    artefacts = analyze_and_report_v2(
        str(sample),
        background_path=str(bg),
        output_dir=str(output_dir),
        sample_mass_kg=0.5,
        write_json=True,
        write_markdown=True,
        write_plots=True,
        write_html=True,
        # F-RPT-03 / v1.18.29 — Technical PDF OFF by default. Flip обратно
        # вручную, если нужен PDF в demo bundle.
        write_technical_pdf=False,
        # F-RPT-04 / v1.18.29 — BecqMoni XML export OFF by default.
        export_becqmoni="off",
    )

    print("\n[V2-full] artefacts:")
    for k, v in artefacts.items():
        if k in ("result", "warnings", "summary", "html_text"):
            continue
        print(f"  {k:14s} {v}")
    if artefacts.get("warnings"):
        print("\n[V2-full] warnings:")
        for w in artefacts["warnings"]:
            print(f"  - {w}")
    print(f"\n✓ V2 full report: {output_dir}")
    return artefacts


if __name__ == "__main__":
    main()
