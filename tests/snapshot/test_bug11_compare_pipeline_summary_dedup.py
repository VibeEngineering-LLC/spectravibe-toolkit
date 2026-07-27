"""BUG-11 / v1.18.30+ (Agent B) — pipeline-summary V2/Production rows
in compare report must NOT contain duplicate (nuclide, peak_channel) pairs.

Symptom: multiplet fitter (post-BUG-3) assigns >1 library line к одному
channel'у. До фикса `_peak_rows` дедуплицировал по (nuclide,
round(peak_E_keV, 0)) — для Th-232 это работало случайно (все компоненты
на ch=203 имеют peak_E_keV=583.01), но не было инвариантом. Для канала
с разнесённой `peak_E_keV` (что встречается при coupled_multiplet fit'е
с гибким positioning) старый ключ оставлял дубли.

Fix: ключ дедупа = (nuclide, peak_channel). Tiebreak: max(S/σ), затем
max(library_I_pct).

Этот тест — синтетический: 5 raw `primary_feps` на ch=200, все Ac-228,
с разными library_E_keV/library_I_pct и разной peak_E_keV (что ломает
старый E_keV-based ключ). Ожидание: после `_peak_rows` остаётся ровно
ОДНА строка с максимальным S/σ.

Также добавлено пересечение с реальным Th-232 demo-снимком (если есть):
ни одна (nuclide, channel) не должна встретиться дважды в выводе.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from gen_v2_compare_th232 import _peak_rows  # noqa: E402


def _make_peak(channel, peak_E, nuclide, lib_E, lib_I, area, sigma_area):
    """Build a primary_feps-shaped dict matching F-380 contract."""
    return {
        "peak_channel": channel,
        "peak_E_keV": peak_E,
        "peak_area_counts": area,
        "peak_area_sigma": sigma_area,
        "peak_area_source": "deconvolved",
        "nuclide": nuclide,
        "library_E_keV": lib_E,
        "library_I_pct": lib_I,
    }


def test_bug11_synthetic_three_lines_one_channel_collapse_to_one():
    """Три library-линии Ac-228 на одном канале — ОДНА строка в выводе,
    с максимальным S/σ. Это база bug-report acceptance criterion.

    Setup: ch=200, peak_E_keV jitter'ит (583.0 / 583.5 / 583.9) — что
    эмулирует случай, когда coupled_multiplet возвращает разнесённую
    E_keV для каждого компонента. Старый round(E_keV) ключ дал бы
    (Ac-228, 583) + (Ac-228, 584) — 2 строки. Новый
    (nuclide, peak_channel) ключ → 1 строка.
    """
    peaks = [
        # 3 line-evidence для Ac-228 на ch=200, разные σ
        _make_peak(200, 583.0, "Ac-228", 562.50, 0.87, 100.0, 50.0),   # S/σ=2.0
        _make_peak(200, 583.5, "Ac-228", 570.91, 0.18, 1000.0, 50.0),  # S/σ=20.0 (winner)
        _make_peak(200, 583.9, "Ac-228", 583.41, 0.11, 50.0, 50.0),    # S/σ=1.0
    ]
    rows = _peak_rows(peaks)

    # ровно 1 строка
    assert len(rows) == 1, (
        f"expected 1 dedup'ed row for (Ac-228, ch=200), got {len(rows)}: "
        f"{[(r['nuclide'], r['channel'], r['library_E_keV']) for r in rows]}"
    )
    r = rows[0]
    assert r["nuclide"] == "Ac-228"
    assert r["channel"] == 200
    # Победитель — строка с max S/σ (1000/50=20.0 → lib_E=570.91)
    assert r["library_E_keV"] == 570.91, (
        f"expected library_E_keV=570.91 (highest S/σ), got {r['library_E_keV']}"
    )
    assert r["sigma"] == 20.0


def test_bug11_tiebreak_library_I_pct():
    """При равных S/σ выбирается строка с max library_I_pct."""
    peaks = [
        # 3 строки с одинаковым area/sigma_area (S/σ identical),
        # разный library_I_pct → выбирается максимальный
        _make_peak(150, 400.0, "Pb-212", 238.63, 1.0, 100.0, 10.0),    # S/σ=10, I=1.0
        _make_peak(150, 400.0, "Pb-212", 300.09, 99.5, 100.0, 10.0),   # S/σ=10, I=99.5 (winner)
        _make_peak(150, 400.0, "Pb-212", 415.00, 50.0, 100.0, 10.0),   # S/σ=10, I=50.0
    ]
    rows = _peak_rows(peaks)
    assert len(rows) == 1
    assert rows[0]["library_I_pct"] == 99.5
    assert rows[0]["library_E_keV"] == 300.09


def test_bug11_different_nuclides_same_channel_kept():
    """Разные нуклиды на одном канале — НЕ дедуплицируются (это
    physically valid случай: Tl-208 233.4 ↔ Pb-212 238.6 на NaI
    с FWHM ~30 keV → один канал, две библиотечных линии).
    """
    peaks = [
        _make_peak(85, 234.6, "Tl-208", 233.36, 0.11, 1000.0, 100.0),
        _make_peak(85, 234.6, "Pb-212", 238.63, 43.6, 500.0, 100.0),
    ]
    rows = _peak_rows(peaks)
    assert len(rows) == 2, (
        f"expected 2 rows (different nuclides), got {len(rows)}: "
        f"{[(r['nuclide'], r['channel']) for r in rows]}"
    )
    nucs = sorted(r["nuclide"] for r in rows)
    assert nucs == ["Pb-212", "Tl-208"]


def test_bug11_phantom_anchors_filtered_before_dedup():
    """Phantom anchors (peak_area_source=library_anchor*) фильтруются ДО
    дедупа. На ch=300 — 1 real + 2 phantom → 1 row (real)."""
    peaks = [
        _make_peak(300, 800.0, "Cs-137", 661.66, 85.1, 5000.0, 100.0),
        # 2 phantom anchors на том же канале
        {
            "peak_channel": 300,
            "peak_E_keV": 800.0,
            "peak_area_counts": None,
            "peak_area_sigma": None,
            "peak_area_source": "library_anchor",
            "nuclide": "Cs-137",
            "library_E_keV": 32.0,
            "library_I_pct": 5.6,
        },
        {
            "peak_channel": 300,
            "peak_E_keV": 800.0,
            "peak_area_counts": None,
            "peak_area_sigma": None,
            "peak_area_source": "library_anchor_phantom",
            "nuclide": "Cs-137",
            "library_E_keV": 283.5,
            "library_I_pct": 0.0058,
        },
    ]
    rows = _peak_rows(peaks)
    assert len(rows) == 1
    assert rows[0]["library_E_keV"] == 661.66


def test_bug11_real_th232_demo_no_dups():
    """Regression: реальный Th-232 demo-snapshot не должен иметь
    duplicate (nuclide, peak_channel) в выводе `_peak_rows`. Если
    GAMMA_DEMO_REPORTS_DIR не задан и sibling-папка отсутствует — skip.
    """
    env = os.environ.get("GAMMA_DEMO_REPORTS_DIR")
    demo_root: Path | None = None
    if env and Path(env).is_dir():
        demo_root = Path(env)
    else:
        sibling = REPO_ROOT.parent / "demo_reports"
        if sibling.is_dir():
            demo_root = sibling
        elif (REPO_ROOT / "demo_reports").is_dir():
            demo_root = REPO_ROOT / "demo_reports"

    if demo_root is None:
        pytest.skip("demo_reports/ не найден")

    # ищем любой Th-232 run с sample/*_report.json
    candidates = [
        c for c in demo_root.iterdir()
        if c.is_dir()
        and (c / "sample").is_dir()
        and list((c / "sample").glob("*_report.json"))
    ]
    if not candidates:
        pytest.skip("нет demo-run с sample/*_report.json")

    report_path = next(iter((candidates[0] / "sample").glob("*_report.json")))
    with report_path.open("r", encoding="utf-8") as f:
        report = json.load(f)
    primary_feps = report.get("primary_feps") or []
    if not primary_feps:
        pytest.skip("primary_feps пуст")

    rows = _peak_rows(primary_feps)
    seen = set()
    for r in rows:
        key = (r["nuclide"], r["channel"])
        assert key not in seen, (
            f"duplicate (nuclide, peak_channel) in real Th-232 dedup output: {key}"
        )
        seen.add(key)
