"""BUG-24 / v1.18.31+ (Agent B) — каждый компонент мультиплета с S>0
должен иметь видимую Гауссову форму на графике.

Symptom (Th-232 demo M5 «Tl-208 510.8 + Ac-228 463.0 + Tl-208 583.2 keV»):
  Ac-228 463.00 keV (I=4.4%, S=4569) показывал вертикальную линию +
  label `463.00`, но НЕ имел видимой заполненной Гауссовой формы.

Root cause: в has_overlay-ветке `_build_multiplets_data` для компонентов
без entry в `overlay_components` (k >= len(overlays)) `g_plus_cont`
ставился равным `cont_arr`, что давало `g_plus_cont - g_base = 0`.
В JS-рендере (interactive_v1_17_2.html:1369) это даёт peakY ≡ continuum
→ Chart.js fill rectangle нулевой высоты — заливки не видно.
HTML-нота legend (`_build_multiplet_blocks`) при этом перечисляла
компонент → асимметрия «легенда есть, заливки нет».

Fix:
  • interactive_html.py:_build_multiplets_data — fallback: при S>0 и
    вырожденном overlay (max(g_plus_cont - g_base) ≈ 0) синтезируем
    Гауссиану из FWHM-модели → видимая форма.
  • interactive_html.py:_build_multiplet_blocks — фильтр S=0 из note
    легенды → симметрия с JS-render gate (c.area === 0 пропуск).
"""
from __future__ import annotations

import sys
import types

sys.path.insert(0, "scripts")

from gamma.reporting.interactive_html import (  # noqa: E402
    _build_multiplets_data,
    _build_multiplet_blocks,
)


# ──── Minimal stubs ───────────────────────────────────────────────────


class _FakeComponent:
    def __init__(self, nuclide, line_E_keV, library_I_pct,
                 center_channel=0, fwhm_channels=4.0):
        self.nuclide = nuclide
        self.line_E_keV = line_E_keV
        self.library_I_pct = library_I_pct
        self.center_channel = center_channel
        self.fwhm_channels = fwhm_channels


class _FakeDeconResult:
    def __init__(self, *, components, areas, roi_low_ch, roi_high_ch,
                 overlay_E_keV=None, overlay_data=None,
                 overlay_continuum=None, overlay_total=None,
                 overlay_components=None, chi2_per_dof=2.0,
                 continuum_model="step_linear", continuum_params=(0.0, 0.0)):
        self.components = tuple(components)
        self.areas = tuple(areas)
        self.roi_low_ch = roi_low_ch
        self.roi_high_ch = roi_high_ch
        self.overlay_E_keV = overlay_E_keV
        self.overlay_data = overlay_data
        self.overlay_continuum = overlay_continuum
        self.overlay_total = overlay_total
        self.overlay_components = overlay_components
        self.chi2_per_dof = chi2_per_dof
        self.continuum_model = continuum_model
        self.continuum_params = continuum_params
        self.cluster_id = ""
        self.phase_A_chi2_per_dof = None


class _FakeSpec:
    """Простой канал-энергия mapping: E = 0.5 * ch (kev/ch)."""
    def __init__(self, n_channels=2048):
        self.counts = [10.0] * n_channels

    def channel_to_energy(self, ch):
        return 0.5 * ch


class _FakeAnalysisResult:
    def __init__(self, decon_results):
        self.deconvolution_results = list(decon_results)
        self.spec = _FakeSpec()
        # FWHM model: const ≈ 5 keV² → FWHM ≈ 2.24 keV (achievable on Gamma-1S)
        self.fwhm_model = (5.0, 0.0, 0.0)


# ──── Tests ──────────────────────────────────────────────────────────


def test_bug24_small_S_component_renders_visible_shape():
    """Компонент с малой ненулевой S и ОТСУТСТВУЮЩИМ overlay → renderer
    должен синтезировать Гауссиану из FWHM модели чтобы заливка была
    видна (g_plus_cont > g_base где-то на интервале)."""
    # ROI 800-1400 ch ≈ 400-700 keV (Tl-208 510 + Ac-228 463 + Tl-208 583).
    roi_lo, roi_hi = 800, 1400
    n_pts = roi_hi - roi_lo
    E_arr = [0.5 * ch for ch in range(roi_lo, roi_hi)]
    cont = [50.0] * n_pts
    # Overlay arrays — фиктивные total/data чтобы has_overlay = True.
    # overlay_components даём только для 2 компонентов из 3 — третий
    # компонент (Ac-228 463) должен попасть в fallback-синтез.
    comp1_overlay = [c + 100.0 for c in cont]  # Tl-208 510
    comp2_overlay = [c + 200.0 for c in cont]  # Tl-208 583
    overlays = (tuple(comp1_overlay), tuple(comp2_overlay))
    total = [a + b - c for a, b, c in zip(comp1_overlay, comp2_overlay, cont)]

    d = _FakeDeconResult(
        components=[
            _FakeComponent("Tl-208", 510.77, 8.12, center_channel=1021),
            _FakeComponent("Tl-208", 583.19, 30.55, center_channel=1166),
            _FakeComponent("Ac-228", 463.00, 4.40, center_channel=926),
        ],
        areas=[54751.0, 205922.0, 4569.0],
        roi_low_ch=roi_lo, roi_high_ch=roi_hi,
        overlay_E_keV=tuple(E_arr),
        overlay_data=tuple(cont),  # placeholder data
        overlay_continuum=tuple(cont),
        overlay_total=tuple(total),
        overlay_components=overlays,  # ТОЛЬКО для k=0,1 — k=2 будет fallback
        chi2_per_dof=15.0,
    )

    out = _build_multiplets_data(
        {"multiplet_deconvolutions": [d.__dict__]},
        _FakeAnalysisResult([d]),
    )
    assert len(out) == 1
    comp_payload = out[0]["components"]
    assert len(comp_payload) == 3

    # Ac-228 463 (index 2 in our setup) must have visible peak above continuum.
    ac228 = next(c for c in comp_payload if c["nuclide"] == "Ac-228")
    assert ac228["area"] == 4569.0
    deltas = [g - b for g, b in zip(ac228["g_plus_cont"], ac228["g_base"])]
    max_delta = max(deltas)
    assert max_delta > 0.5, (
        f"BUG-24: Ac-228 463 keV компонент с S=4569 должен дать видимый "
        f"подъём над continuum (max_delta={max_delta:.3f}). Renderer "
        f"должен синтезировать Гауссиану из FWHM-модели, когда overlay "
        f"отсутствует/вырожден."
    )


def test_bug24_zero_S_component_skipped_from_legend_note():
    """HTML-note legend в `_build_multiplet_blocks` должна СКРЫВАТЬ
    компоненты с S=0 (симметрия с JS-render gate c.area === 0)."""
    mp = {
        "id": "M1",
        "title": "тест",
        "chi2_per_dof": 2.0,
        "closure_pct": 0.5,
        "n_channels": 100,
        "continuum_model": "linear",
        "components": [
            {"nuclide": "Tl-208", "E_keV": 510.77, "I_pct": 8.12,
             "area": 54751.0},
            {"nuclide": "Ac-228", "E_keV": 463.0, "I_pct": 4.40,
             "area": 0.0},                              # phantom S=0
            {"nuclide": "Tl-208", "E_keV": 583.19, "I_pct": 30.55,
             "area": 205922.0},
        ],
    }
    html = _build_multiplet_blocks([mp])
    # Phantom-S=0 компонент НЕ должен попасть в HTML note (легенду).
    assert "463" not in html.split("<p class=\"fp-mp-note\">")[1].split("</p>")[0], (
        "BUG-24: компонент с S=0 не должен фигурировать в HTML note-легенде "
        "мультиплета — JS не рендерит его как dataset (gate c.area === 0), "
        "так что легенда обязана быть симметрична."
    )
    # Контроль: ненулевые компоненты — в легенде есть.
    assert "510.8" in html or "510" in html
    assert "583.2" in html or "583" in html


def test_bug24_nonzero_S_component_listed_in_legend_note():
    """Малое НЕНУЛЕВОЕ S (S=4569) должно фигурировать в HTML note
    легенде — это пара к видимой Гауссиане (test_bug24_small_S above)."""
    mp = {
        "id": "M1",
        "title": "тест",
        "chi2_per_dof": 2.0,
        "closure_pct": 0.5,
        "n_channels": 100,
        "continuum_model": "linear",
        "components": [
            {"nuclide": "Tl-208", "E_keV": 510.77, "I_pct": 8.12,
             "area": 54751.0},
            {"nuclide": "Ac-228", "E_keV": 463.0, "I_pct": 4.40,
             "area": 4569.0},
            {"nuclide": "Tl-208", "E_keV": 583.19, "I_pct": 30.55,
             "area": 205922.0},
        ],
    }
    html = _build_multiplet_blocks([mp])
    note_block = html.split('<p class="fp-mp-note">')[1].split("</p>")[0]
    assert "463" in note_block and "4569" in note_block, (
        "BUG-24: компонент Ac-228 463.0 с S=4569 ДОЛЖЕН быть в HTML note "
        "(он же должен рендериться на canvas — см. парный тест)."
    )


if __name__ == "__main__":
    test_bug24_small_S_component_renders_visible_shape()
    test_bug24_zero_S_component_skipped_from_legend_note()
    test_bug24_nonzero_S_component_listed_in_legend_note()
    print("OK")
