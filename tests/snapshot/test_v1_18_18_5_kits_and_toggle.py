# -*- coding: utf-8 -*-
"""F-331 + F-332 / v1.18.18.5 — Reference kits + HTML 4-way spectrum toggle.

F-331 контракт: canonical {geometry}/{nuclide}/{sample,background}.spe
структура в `detectors/Gamma-1S/reference_spectra/reference_kits/`.
F-332 контракт: HTML report включает chart-toggle payload (gross + bg
+ net на единой энергетической оси).
"""
from __future__ import annotations
import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

KITS_ROOT = ROOT / "detectors" / "Gamma-1S" / "reference_spectra" / "reference_kits"


def _manifest_path(rel: str) -> Path:
    """v1.27.9 lane A — normalise a path string read from MANIFEST.json.

    The manifest writer (`scripts/build_reference_kits.py`) currently stores
    paths with `\\` separators (Windows-native at write-time). On a Linux
    CI runner `ROOT / "detectors\\Gamma-1S\\..."` produces a path with a
    literal `\\` inside the *filename* component → IsADirectoryError /
    FileNotFoundError. Splitting on either separator and joining via
    Path/`/` operator gives a cross-platform safe result without touching
    the data file itself.
    """
    parts = rel.replace("\\", "/").split("/")
    return ROOT.joinpath(*parts)


# ─── F-331 ──────────────────────────────────────────────────────────

def test_F331_kits_root_present():
    """reference_kits/ существует и содержит README + MANIFEST + 5 geom."""
    assert KITS_ROOT.is_dir(), "reference_kits/ отсутствует"
    assert (KITS_ROOT / "README.md").is_file()
    assert (KITS_ROOT / "MANIFEST.json").is_file()
    expected_geoms = {"Marinelli_1L", "Point_5cm", "Point_25cm",
                      "Petri_60mL", "Denta_120mL"}
    actual = {p.name for p in KITS_ROOT.iterdir() if p.is_dir()}
    assert expected_geoms.issubset(actual), (
        f"Missing geometries: {expected_geoms - actual}"
    )


def test_F331_marinelli_kits_have_sample_and_background():
    """Каждый Marinelli/{nuclide}/ leaf содержит ровно один sample + один
    background. v1.18.24.0: проверка через MANIFEST.json как single source
    of truth (вместо хардкода naming-convention `sample_*.spe` /
    `background_*.spe`) — позволяет user replace kit с произвольным
    naming. Для Th-232: replaced 2024-10-24 → Th232_420-7-17_*.spe +
    Фон закр кр вода_*.spe."""
    m = json.loads((KITS_ROOT / "MANIFEST.json").read_text(encoding="utf-8"))
    found_by_nuc = {}
    for row in m["kits"]:
        if row.get("geometry") != "Marinelli_1L":
            continue
        nuc = row.get("nuclide")
        found_by_nuc.setdefault(nuc, []).append(row)
    for nuc in ("Cs-137", "K-40", "Ra-226", "Th-232"):
        leaf = KITS_ROOT / "Marinelli_1L" / nuc
        assert leaf.is_dir(), f"missing kit leaf: {leaf}"
        rows = found_by_nuc.get(nuc, [])
        assert len(rows) == 1, (
            f"{nuc}: expected 1 manifest row, got {len(rows)}"
        )
        row = rows[0]
        sample = _manifest_path(row["sample_kit_path"])
        bg = _manifest_path(row["background_kit_path"])
        assert sample.is_file(), f"{nuc}: sample missing: {sample}"
        assert bg.is_file(), f"{nuc}: bg missing: {bg}"


def test_F331_manifest_matches_filesystem():
    """MANIFEST.json content lines up with actual files on disk."""
    m = json.loads((KITS_ROOT / "MANIFEST.json").read_text(encoding="utf-8"))
    assert "kits" in m
    assert m["version"] == "v1.18.18.5"
    assert len(m["kits"]) >= 30, "expected ≥30 kits"
    for row in m["kits"]:
        sample = _manifest_path(row["sample_kit_path"])
        bg = _manifest_path(row["background_kit_path"])
        assert sample.is_file(), f"sample missing: {sample}"
        assert bg.is_file(), f"bg missing: {bg}"


def test_F331_legacy_subdir_gone():
    """Legacy flat folder retired (replaced by archive/ + reference_kits/)."""
    legacy = (ROOT / "detectors" / "Gamma-1S" / "reference_spectra"
              / "Gamma-1S_NaI_63x63_USB_SN-01")
    assert not legacy.exists(), (
        f"Legacy folder still present: {legacy}"
    )


def test_F331_archive_contains_leftovers():
    """archive/ has the historical spectra (≥150 .spe files)."""
    archive = ROOT / "detectors" / "Gamma-1S" / "reference_spectra" / "archive"
    assert archive.is_dir()
    n_spe = sum(1 for _ in archive.rglob("*.spe"))
    assert n_spe >= 150, f"expected ≥150 archived .spe, got {n_spe}"


# ─── F-332 — HTML chart toggle ──────────────────────────────────────

DEMO_HTML = (
    ROOT / "demo_reports" / "v1_18_18" / "M_cs_легкий_report.html"
)


def _read_demo():
    if not DEMO_HTML.is_file():
        pytest.skip(f"demo HTML missing: {DEMO_HTML}")
    return DEMO_HTML.read_text(encoding="utf-8")


def test_F332_toggle_ui_present():
    """HTML report содержит 4-way toggle UI элементы."""
    html = _read_demo()
    assert "fp-view-toggle" in html
    for view in ("sample", "bg", "overlay", "net"):
        assert f'data-view="{view}"' in html, f"missing button data-view={view}"


def test_F332_toggle_button_labels():
    """Russian labels for 4 view modes are rendered."""
    html = _read_demo()
    # F-335.7 / v1.18.18.13 — «Чистый (выч.)» переименовано в «Вычет».
    for label in ("Образец", "Фон", "Образец + Фон", "Вычет"):
        assert label in html, f"missing button label: {label}"


def test_F332_chart_payload_serialised():
    """`const CHART={...}` injected with has_background=true (M_cs bg
    был применён)."""
    html = _read_demo()
    m = re.search(r"const CHART=(\{[^;]*?\});", html, re.DOTALL)
    assert m, "DATA_CHART not injected"
    payload = json.loads(m.group(1))
    assert payload.get("has_background") is True
    for k in ("E", "C_net", "C_gross", "C_bg", "t_sample", "t_bg",
              "bg_scale"):
        assert k in payload, f"missing key in payload: {k}"
    # Все 4 массива должны иметь одинаковую длину (= общая E-сетка).
    n = len(payload["E"])
    assert n > 100
    assert len(payload["C_net"]) == n
    assert len(payload["C_gross"]) == n
    assert len(payload["C_bg"]) == n


def test_F332_bg_scale_is_live_time_ratio():
    """`bg_scale` ≈ t_sample / t_bg (как в subtract_background)."""
    html = _read_demo()
    m = re.search(r"const CHART=(\{[^;]*?\});", html, re.DOTALL)
    payload = json.loads(m.group(1))
    expected = payload["t_sample"] / payload["t_bg"]
    assert payload["bg_scale"] == pytest.approx(expected, rel=1e-4)


def test_F332_setView_function_present():
    """JS handler `setView()` injected to drive the toggle.

    F-334 / v1.18.18.7 — переход с `function setView(view){...}` на
    `let setView = function(view){...}` чтобы можно было определить
    no-op fallback до if-branch с has_background. Тест принимает оба
    варианта.
    """
    html = _read_demo()
    assert ("function setView(" in html) or ("setView = function(" in html), \
        "setView handler missing in HTML"
    # And the dispatch on click for each button.
    assert "fp-view-btn" in html
    # And the legend toggle (overlay shows legend).
    assert "plugins.legend.display" in html


def test_F332_payload_calibrated_to_sample_energy_axis():
    """Все 3 серии (gross / bg / net) выровнены на одной E-сетке
    (sample's channel→keV mapping). E ascending; 47–3000 keV ranges
    для NaI Gamma-1S 3000-keV ceiling."""
    html = _read_demo()
    m = re.search(r"const CHART=(\{[^;]*?\});", html, re.DOTALL)
    payload = json.loads(m.group(1))
    E = payload["E"]
    # Monotonic ascending
    for i in range(1, len(E)):
        assert E[i] >= E[i - 1], f"E not monotone at {i}"
    # Range sanity: Gamma-1S 47-3000 keV
    assert E[0] >= 30 and E[0] <= 100
    assert E[-1] >= 2500 and E[-1] <= 3100


# ─── F-332 — non-bg case ────────────────────────────────────────────

def test_F332_toggle_default_hidden_via_inline_css():
    """Toggle UI стартует с inline `display:none`; JS показывает его
    только когда `CHART.has_background === true`. Это гарантирует, что
    spectrum-only reports не показывают broken toggle."""
    html = _read_demo()
    assert 'id="fp-view-toggle" style="display:none;"' in html, (
        "Toggle div must start hidden; JS unhides only if "
        "has_background=true"
    )
