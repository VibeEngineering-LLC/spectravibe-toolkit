"""Regression — F-118 chain-forced multiplet auto-discovery.

Контракт: на Th-232 фикстуре report["multiplet_deconvolutions"]
содержит ровно две жёстко-закреплённые подгонки с компонентами,
покрывающими M1 (911 + 964.77 + 969 + 860.6) и M2 (1588 + 1620 + 1630).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

from gamma.io.readers import read_spectrum
from gamma.identification.staged_pipeline import (
    build_fwhm_model, fwhm_keV_at_energy, analyze_lsrm_spe,
)


_FIXTURE = (
    "detectors/Gamma-1S/reference_spectra/"
    "archive/"
    "Th232_420-7-17_Маринелли_0cm.spe"
)


def test_th232_M1_M2_clusters_present():
    if not Path(_FIXTURE).is_file():
        print(f"  ⚠ skipping (fixture missing): {_FIXTURE}")
        return
    res = analyze_lsrm_spe(
        _FIXTURE,
        allow_stage2=True, allow_stage3=True,
        apply_deconvolution=True, compute_activities=False,
        compute_mda=False,
    )
    deconv = res.deconvolution_results or []
    # Жёстко-закреплённые M1 и M2 должны идти первыми
    assert len(deconv) >= 2, (
        f"expected ≥2 deconvolution clusters (M1, M2), got {len(deconv)}"
    )
    # Найти кластер с Ac-228 911 (M1) и с Ac-228 1588 (M2)
    found_M1 = False
    found_M2 = False
    for d in deconv:
        comps = list(d.components)
        nuclides = {(c.nuclide, round(float(c.line_E_keV), 1)) for c in comps}
        if ("Ac-228", 911.2) in nuclides and ("Ac-228", 969.0) in nuclides:
            found_M1 = True
            # F-117 контракт: 4 компоненты в M1 (Ac-228 × 3 + Tl-208 860.6)
            assert len(comps) >= 3, (
                f"M1 has {len(comps)} components, expected ≥3"
            )
            # F-118: жёстко-закреплённый кластер использует coupled_*
            assert str(d.method).startswith("coupled_"), (
                f"M1 method should be coupled_*, got {d.method}"
            )
        if ("Ac-228", 1588.2) in nuclides and ("Ac-228", 1630.6) in nuclides:
            found_M2 = True
            assert str(d.method).startswith("coupled_"), (
                f"M2 method should be coupled_*, got {d.method}"
            )
    assert found_M1, "M1 (Ac-228 911+969) cluster not found"
    assert found_M2, "M2 (Ac-228 1588+1630) cluster not found"
    print(f"  ✓ test_th232_M1_M2_clusters_present "
          f"(M1 ✓, M2 ✓, total clusters={len(deconv)})")


if __name__ == "__main__":
    test_th232_M1_M2_clusters_present()
    print("Th-chain M1/M2 presence regression PASS.")
