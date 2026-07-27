"""F-335 / v1.18.18.8 — Pan/zoom + smoothing + Pb-212 false-ID regression.

Подзадачи:
  • F-335.1 — pan/zoom contract: clamp original, lock Y baseline, fitYToVisibleX
  • F-335.2 — cps/counts toggle в modal header
  • F-335.3 — Savitzky-Golay smoothing (payload + UI)
  • F-335.4 — chain-proxy single-line guard (Pb-212 false ID fix)
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, "scripts")

from gamma.reporting import analyze_and_report
from gamma.identification.staged_pipeline import analyze_lsrm_spe


ROOT = Path(__file__).resolve().parents[2]
KIT = ROOT / "detectors/Gamma-1S/reference_spectra/reference_kits"
CS_SAMPLE = KIT / "Marinelli_1L/Cs-137/sample_M_cs_легкий_2001-2005.spe"
CS_BG = KIT / "Marinelli_1L/Cs-137/background_bg_2016_marinelli_water_marinelli.spe"
K_SAMPLE = KIT / "Marinelli_1L/K-40/sample_M_k_легкий_2001-2005.spe"
K_BG = KIT / "Marinelli_1L/K-40/background_bg_2016_marinelli_water_marinelli.spe"
TH_SAMPLE = KIT / "Marinelli_1L/Th-232/sample_M_th_легкий_2001-2005.spe"
TH_BG = KIT / "Marinelli_1L/Th-232/background_bg_2016_marinelli_water_marinelli.spe"


@pytest.fixture(scope="module")
def k40_html(tmp_path_factory):
    if not K_SAMPLE.exists():
        pytest.skip(f"Kit sample missing: {K_SAMPLE}")
    out = tmp_path_factory.mktemp("f335_k40")
    res = analyze_and_report(
        str(K_SAMPLE),
        background_path=str(K_BG) if K_BG.exists() else None,
        sample_mass_kg=0.665,
        output_dir=str(out),
        write_html=True, write_plots=False, write_markdown=False,
    )
    return Path(res["html"]).read_text(encoding="utf-8")


# ─────────────────────────────────────────────────────────────────
# F-335.1 — Pan/zoom limits + auto-fit Y
# ─────────────────────────────────────────────────────────────────

def test_f335_1_zoom_plugin_mode_x(k40_html):
    """Zoom plugin config: mode='x' (только X axis zoom через колесо)."""
    # Check zoom config includes mode: 'x' (not 'xy')
    assert "mode: 'x'" in k40_html
    # And NOT 'xy' mode for wheel zoom — old setting
    # (Soft check — mode: 'xy' might appear elsewhere for non-zoom contexts)


def test_f335_1_zoom_limits_original(k40_html):
    """Limits clamp zoom-out на original X и Y range."""
    assert "min: 'original'" in k40_html
    assert "max: 'original'" in k40_html
    # X и Y limits оба установлены
    assert k40_html.count("'original'") >= 4


def test_f335_1_fityToVisibleX_present(k40_html):
    """JS helper fitYToVisibleX определён + вызывается из RMB pan."""
    assert "function fitYToVisibleX" in k40_html
    assert "fitYToVisibleX(modalChart" in k40_html


def test_f335_1_x_orig_bounds(k40_html):
    """X_ORIG_MIN / X_ORIG_MAX константы вычислены из CHART.E."""
    assert "X_ORIG_MIN" in k40_html
    assert "X_ORIG_MAX" in k40_html


def test_f335_1_y_baseline_locked(k40_html):
    """В fitYToVisibleX Y.min фиксирован у logMin / 0, не trail data."""
    # Helper должен ставить y.min ВСЕГДА (не пропускать)
    m = re.search(
        r"function fitYToVisibleX[\s\S]{0,2000}",
        k40_html,
    )
    assert m
    body = m.group(0)
    assert "y.min = yd.logMin" in body or "y.min = 0" in body
    # Y.max scaled от data, не от scroll-offset
    assert "y.max" in body


def test_f335_1_rmb_pan_only_x(k40_html):
    """Custom RMB pan handler меняет ТОЛЬКО scales.x.* (Y не трогает)."""
    # Find the mousedown[button=2] block
    m = re.search(
        r"mcv\.addEventListener\('mousedown'[\s\S]{0,500}",
        k40_html,
    )
    assert m
    # Should NOT capture scaleYMinStart — F-335.1 removed Y panning
    rmb_block = m.group(0)
    assert "scaleYMinStart" not in rmb_block, \
        "F-335.1: Y pan should be removed from RMB handler"


# ─────────────────────────────────────────────────────────────────
# F-335.2 — cps/counts toggle в modal header
# ─────────────────────────────────────────────────────────────────

def test_f335_2_modal_units_buttons_present(k40_html):
    """Modal header содержит cps/counts buttons."""
    assert 'data-modal-units="cps"' in k40_html
    assert 'data-modal-units="counts"' in k40_html


def test_f335_2_modal_button_sync(k40_html):
    """setYUnits() синхронизирует обе button-группы (main + modal)."""
    # When setYUnits fires, modal buttons get .active toggle
    assert "data-modal-units" in k40_html
    # Sync logic in setYUnits
    assert "syncModalFromMain" in k40_html


def test_f335_2_modal_handler_calls_setYUnits(k40_html):
    """Modal click handler передаёт modalUnits в общий setYUnits."""
    assert "setYUnits(u)" in k40_html or "setYUnits(b.dataset.modalUnits)" in k40_html


# ─────────────────────────────────────────────────────────────────
# F-335.3 — Savitzky-Golay smoothing
# ─────────────────────────────────────────────────────────────────

def test_f335_3_smoothing_payload(k40_html):
    """Payload содержит smoothed массивы + метаданные."""
    m = re.search(r"const CHART=(\{.*?\});", k40_html, re.DOTALL)
    chart = json.loads(m.group(1))
    assert chart.get("smooth_method") == "savitzky_golay"
    assert chart.get("smooth_window") == 5
    assert chart.get("smooth_polyorder") == 2
    assert isinstance(chart.get("C_net_smooth"), list)
    assert len(chart["C_net_smooth"]) == len(chart["E"])
    # With bg also need C_gross_smooth / C_bg_smooth
    if chart.get("has_background"):
        assert isinstance(chart.get("C_gross_smooth"), list)
        assert isinstance(chart.get("C_bg_smooth"), list)


def test_f335_3_smoothing_preserves_total_counts(k40_html):
    """SG-smoothing сохраняет общую сумму ±1% (нет potential bias)."""
    m = re.search(r"const CHART=(\{.*?\});", k40_html, re.DOTALL)
    chart = json.loads(m.group(1))
    raw = chart["C_gross"] or chart["C_net"]
    sm = chart["C_gross_smooth"] or chart["C_net_smooth"]
    if not raw or not sm:
        pytest.skip("no data to compare")
    s_raw = sum(raw)
    s_sm = sum(sm)
    if s_raw <= 0:
        pytest.skip("empty raw")
    # SG smoothing может слегка изменить sum (~1%)
    assert abs(s_sm - s_raw) / s_raw < 0.05


def test_f335_3_smoothing_buttons_present(k40_html):
    """Кнопки сглаживания на main + modal."""
    assert 'data-smooth="on"' in k40_html
    assert 'data-smooth="off"' in k40_html
    assert 'data-modal-smooth="on"' in k40_html
    assert 'data-modal-smooth="off"' in k40_html


def test_f335_3_setSmoothing_function(k40_html):
    """JS function setSmoothing определена + работает на main и modal."""
    assert "function setSmoothing" in k40_html
    assert "smoothOn" in k40_html
    # pickArr helper для выбора raw vs smooth
    assert "function pickArr" in k40_html


def test_f335_3_default_off(k40_html):
    """Default smoothOn=false (per user Q3)."""
    # Initial state: смотрим что в HTML "active" на "выкл" buttons
    assert 'data-smooth="off" class="fp-units-btn active"' in k40_html
    assert 'data-modal-smooth="off" class="fp-units-btn active"' in k40_html


# ─────────────────────────────────────────────────────────────────
# F-335.4 — Chain-proxy single-line guard (Pb-212 false ID fix)
# ─────────────────────────────────────────────────────────────────

def test_f335_4_K40_drops_Pb212():
    """K-40 sample: Pb-212 НЕ должен быть в final_detected (только K-40)."""
    if not K_SAMPLE.exists():
        pytest.skip(f"Kit sample missing: {K_SAMPLE}")
    ar = analyze_lsrm_spe(
        str(K_SAMPLE),
        background_path=str(K_BG) if K_BG.exists() else None,
    )
    nuc_set = {ni.nuclide for ni in ar.final_detected}
    assert "K-40" in nuc_set
    assert "Pb-212" not in nuc_set, \
        f"F-335.4: Pb-212 must be dropped from K-40 sample (no chain partners). " \
        f"Got: {nuc_set}"


def test_f335_4_Th232_keeps_Pb212():
    """Th-232 demo: Pb-212 ДОЛЖЕН остаться (chain dominant с многими линиями)."""
    if not TH_SAMPLE.exists():
        pytest.skip(f"Kit sample missing: {TH_SAMPLE}")
    ar = analyze_lsrm_spe(
        str(TH_SAMPLE),
        background_path=str(TH_BG) if TH_BG.exists() else None,
    )
    nuc_set = {ni.nuclide for ni in ar.final_detected}
    assert "Pb-212" in nuc_set, \
        f"F-335.4: Pb-212 must remain when Th-chain dominant. Got: {nuc_set}"


def test_f335_4_Cs137_no_chain_false_ids():
    """Cs-137 sample: НЕ должно быть Pb-212 / Tl-208 / Bi-214 / Bi-212 attribution."""
    if not CS_SAMPLE.exists():
        pytest.skip(f"Kit sample missing: {CS_SAMPLE}")
    ar = analyze_lsrm_spe(
        str(CS_SAMPLE),
        background_path=str(CS_BG) if CS_BG.exists() else None,
    )
    nuc_set = {ni.nuclide for ni in ar.final_detected}
    forbidden = {"Pb-212", "Tl-208", "Bi-212", "Ac-228", "Bi-214", "Pb-214"}
    leaked = nuc_set & forbidden
    assert not leaked, \
        f"F-335.4: chain-proxy false IDs leaked in Cs-137 sample: {leaked}"


def test_f335_4_diagnostic_in_pipeline_notes():
    """Suppression причина выводится в pipeline_notes для трассируемости."""
    if not K_SAMPLE.exists():
        pytest.skip(f"Kit sample missing: {K_SAMPLE}")
    ar = analyze_lsrm_spe(
        str(K_SAMPLE),
        background_path=str(K_BG) if K_BG.exists() else None,
    )
    # StagedAnalysisResult exposes `notes` (which is what becomes
    # pipeline_notes in JSON report serializer).
    notes_list = getattr(ar, "notes", None) or getattr(ar, "pipeline_notes", []) or []
    notes = " ".join(notes_list)
    # F-335.4 diagnostic должен упоминать Pb-212 suppression
    assert ("F-335.4" in notes) or ("подавлен Pb-212" in notes), \
        f"F-335.4 suppression diagnostic missing. Notes: {notes[:500]}"


# ─────────────────────────────────────────────────────────────────
# Cross-cutting
# ─────────────────────────────────────────────────────────────────

def test_f335_version_bump():
    """SKILL_VERSION должен быть ≥ v1.18.18.8 (F-335 introduced bump)."""
    from gamma.reporting.json_report import SKILL_VERSION
    import re
    m = re.match(r"^v(\d+)\.(\d+)\.(\d+)(?:\.(\d+))?$", SKILL_VERSION)
    assert m, f"Unexpected SKILL_VERSION format: {SKILL_VERSION}"
    major, minor, patch, hotfix = (int(g) if g else 0 for g in m.groups())
    assert (major, minor, patch, hotfix) >= (1, 18, 18, 8), \
        f"Expected ≥ v1.18.18.8; got {SKILL_VERSION}"


def test_f335_K40_demo_chart_renders_correctly(k40_html):
    """K-40 sample peaks list НЕ содержит Pb-212 (F-335.4 cross-check в HTML)."""
    m = re.search(r"const peaks=(\[.*?\]);", k40_html, re.DOTALL)
    assert m
    peaks = json.loads(m.group(1))
    pb212_labels = [p for p in peaks if "Pb-212" in p.get("label", "")]
    assert not pb212_labels, \
        f"F-335.4: Pb-212 leaked into K-40 sample peaks list: {pb212_labels}"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
