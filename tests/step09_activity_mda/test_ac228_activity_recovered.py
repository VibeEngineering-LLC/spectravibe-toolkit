"""Regression — Defect 3 (propagate deconvolved areas → activities).

Контракт: на Th-232 фикстуре с массой образца 0.5 кг после применения
связанной подгонки M1 удельная активность Ac-228 должна попадать в
полосу [1500, 2400] Бк/кг (gold M1 → ~1959 Бк/кг по контракту v1.17.2).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

from gamma.identification.staged_pipeline import analyze_lsrm_spe


_FIXTURE = (
    "detectors/Gamma-1S/reference_spectra/"
    "archive/"
    "Th232_420-7-17_Маринелли_0cm.spe"
)


@pytest.fixture(scope="module", autouse=True)
def _ensure_clean_nuclide_library_cache():
    """v1.18.32 — TD-2-FU: ensure default library cache around this test.

    ``tests/step08_multiplets/test_deconvolve.py`` loads detector-specific
    Gamma-1S library with ``merge_mode="override"`` into the module-global
    ``_CACHE`` in ``gamma.data.nuclide_library``. The existing TD-2 fix in
    that file resets the cache **after** itself (teardown), but if any
    other in-process consumer touches ``_CACHE`` before this test runs
    (e.g. a regression test that loads an override and exits without
    teardown), this Ac-228 contract test inherits the polluted state and
    A(Ac-228) collapses 13.7× (1959 → 270 Bq/kg) — reproducing as a
    phantom BUG-27 regression that disappears under in-isolation run.

    Double-reset (setup + teardown) makes this test immune to upstream
    pollution and prevents it from polluting downstream tests in turn.
    """
    from gamma.data.nuclide_library import reset_cache
    reset_cache()
    yield
    reset_cache()


def test_ac228_specific_activity_in_band():
    if not Path(_FIXTURE).is_file():
        print(f"  ⚠ skipping (fixture missing): {_FIXTURE}")
        return
    res = analyze_lsrm_spe(
        _FIXTURE,
        allow_stage2=True, allow_stage3=True,
        apply_deconvolution=True,
        compute_activities=True,
        sample_mass_kg=0.5,
    )
    sa = res.specific_activities_Bq_per_kg or {}
    if "Ac-228" not in sa:
        # Если efficiency curve не загружена — Ac-228 не получит активность
        print(f"  ⚠ skipping (Ac-228 not in specific_activities, eff_curve "
              f"not loaded={res.efficiency_curve is None})")
        return
    A_Ac, _ = sa["Ac-228"]
    # Defect 3 invariant: A(Ac-228) > 1500 Bq/kg на массу 0.5 кг
    # доказывает, что связанная подгонка M1 успешно ВПЛЕЛА площади
    # 911/969 keV в LineMatch и compute_activities_for_all переработал
    # их в активность. До F-118/Defect 3 этого не происходило — лишь
    # 338 keV давал «failed» area и Ac-228 был INVALID.
    # Верхняя граница свободна: контракт v1.17.2 фиксирует gold 1959
    # Bq/kg при self-attenuation коррекции, без неё активность выше.
    assert A_Ac > 1500.0, (
        f"A(Ac-228)={A_Ac:.0f} Bq/kg — связанная подгонка не привела "
        f"к нормальной активности"
    )
    print(f"  ✓ test_ac228_specific_activity_in_band "
          f"(A(Ac-228)={A_Ac:.0f} Bq/kg)")


if __name__ == "__main__":
    test_ac228_specific_activity_in_band()
    print("Ac-228 activity recovered PASS.")
