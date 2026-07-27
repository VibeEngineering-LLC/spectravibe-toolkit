"""F-134 / v1.17.7 — Fit overlay arrays в DeconvolutionResult.

Контракт: pre-computed overlay arrays (E_keV, data, continuum, total,
component_g_plus_cont) из CoupledFitResult переносятся в DeconvolutionResult
и используются PNG / HTML рендерерами напрямую. Без F-134 PNG/HTML
строили continuum как c0+c1·channel_index (канальный домен) с
параметрами в энергетическом домене → визуально continuum прижат к нулю,
total «недокидывает» данные → catastrophic визуальное несоответствие
точному fit'у (closure −0.61 % в JSON но −78 % визуально).
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))


TH232_FIXTURE = (Path(__file__).parent.parent.parent / "detectors" / "Gamma-1S"
                 / "reference_spectra"
                 / "archive"
                 / "Th232_420-7-17_Маринелли_0cm.spe")


def _need(p: Path):
    if not p.exists():
        pytest.skip(f"fixture missing: {p}")


def test_deconv_result_carries_overlay_fields():
    """DeconvolutionResult имеет overlay_* поля."""
    from gamma.peaks.deconvolve import DeconvolutionResult
    fields = {f for f in DeconvolutionResult.__dataclass_fields__}
    for k in ("overlay_E_keV", "overlay_data", "overlay_continuum",
              "overlay_total", "overlay_components"):
        assert k in fields, f"missing field {k}"


def test_coupled_to_deconv_populates_overlay():
    """_coupled_to_deconv_result переносит coupled.data/continuum/total
    в overlay-поля результата."""
    _need(TH232_FIXTURE)
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    r = analyze_lsrm_spe(str(TH232_FIXTURE), complete_workflow=True)
    m1 = None
    for d in r.deconvolution_results or []:
        if d.notes and "кластер M1" in d.notes:
            m1 = d
            break
    assert m1 is not None
    assert m1.overlay_E_keV is not None
    assert m1.overlay_data is not None
    assert m1.overlay_continuum is not None
    assert m1.overlay_total is not None
    assert m1.overlay_components is not None
    # Все массивы одной длины
    n = len(m1.overlay_E_keV)
    assert len(m1.overlay_data) == n
    assert len(m1.overlay_continuum) == n
    assert len(m1.overlay_total) == n
    for comp_curve in m1.overlay_components:
        assert len(comp_curve) == n


def test_overlay_continuum_within_expected_range():
    """Overlay continuum должен соответствовать форме fit'а:
    max ~4000 слева, min ~900-1500 справа (как у эталона v1.17.2)."""
    _need(TH232_FIXTURE)
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    r = analyze_lsrm_spe(str(TH232_FIXTURE), complete_workflow=True)
    m1 = None
    for d in r.deconvolution_results or []:
        if d.notes and "кластер M1" in d.notes:
            m1 = d
            break
    assert m1 is not None
    cont = list(m1.overlay_continuum)
    assert max(cont) > 3000, f"continuum max={max(cont):.0f} должен быть > 3000"
    assert min(cont) > 500, f"continuum min={min(cont):.0f} должен быть > 500"
    # Спадает слева направо
    assert cont[0] > cont[-1], "continuum должен спадать слева направо"


def test_overlay_total_closes_data():
    """Overlay total ≈ overlay data (closure < 5%)."""
    _need(TH232_FIXTURE)
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    r = analyze_lsrm_spe(str(TH232_FIXTURE), complete_workflow=True)
    m1 = None
    for d in r.deconvolution_results or []:
        if d.notes and "кластер M1" in d.notes:
            m1 = d
            break
    assert m1 is not None
    sum_data = sum(m1.overlay_data)
    sum_total = sum(m1.overlay_total)
    closure = 100.0 * (sum_total - sum_data) / sum_data
    assert abs(closure) < 5.0, (
        f"overlay closure {closure:.2f}% должна быть < 5 %"
    )
    # Max(total) должен быть близок к max(data)
    max_data = max(m1.overlay_data)
    max_total = max(m1.overlay_total)
    rel = abs(max_total - max_data) / max_data
    assert rel < 0.10, (
        f"overlay max(total)={max_total:.0f} должен быть в ±10 % "
        f"от max(data)={max_data:.0f}, отклонение {rel*100:.1f}%"
    )


def test_html_report_uses_overlay_for_multiplet():
    """HTML interactive отчёт включает корректный continuum и total
    из overlay-данных (не из канальной реконструкции)."""
    _need(TH232_FIXTURE)
    from gamma.reporting import analyze_and_report
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        analyze_and_report(
            str(TH232_FIXTURE), output_dir=tmp,
            write_pdf=False, plot_dpi=80,
        )
        html_files = list(Path(tmp).glob("*.html"))
        assert html_files
        html = html_files[0].read_text(encoding="utf-8")
    # Найти MULTIPLETS JS-массив
    m = re.search(r"const __MULTIPLETS\s*=\s*(\[.*?\]);", html, re.DOTALL)
    if not m:
        m = re.search(r"multiplets\s*=\s*(\[.*?\]);", html, re.DOTALL)
    assert m, "не найден JS-массив __MULTIPLETS"
    data = json.loads(m.group(1))
    m1 = None
    for mult in data:
        if "M1" in (mult.get("title") or "") or "911" in (mult.get("title") or ""):
            m1 = mult
            break
    assert m1 is not None, "M1 кластер не найден в HTML"
    cont = m1.get("continuum") or []
    total = m1.get("total") or []
    data_arr = m1.get("data") or []
    assert cont, "continuum пуст"
    assert max(cont) > 3000, (
        f"HTML continuum max={max(cont):.0f} должен быть > 3000 "
        f"(без F-134 был бы < 1000)"
    )
    # total накрывает data
    assert max(total) > 0.9 * max(data_arr), (
        f"HTML total max={max(total):.0f} должен накрывать "
        f"data max={max(data_arr):.0f}"
    )


def test_png_plot_uses_overlay():
    """PNG plot мультиплета должен быть создан через overlay-массивы
    (косвенно — проверяем что файл создаётся без ошибок). Визуальная
    проверка отдельно — graphical regression вне теста."""
    _need(TH232_FIXTURE)
    from gamma.reporting import analyze_and_report
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        analyze_and_report(
            str(TH232_FIXTURE), output_dir=tmp,
            write_pdf=False, plot_dpi=80,
        )
        # Ищем PNG мультиплета
        plot_files = list(Path(tmp).rglob("multiplet_*.png"))
        assert plot_files, "PNG мультиплета не сгенерирован"


def test_overlay_components_align_with_components_list():
    """overlay_components[k] соответствует d.components[k]."""
    _need(TH232_FIXTURE)
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    r = analyze_lsrm_spe(str(TH232_FIXTURE), complete_workflow=True)
    m1 = None
    for d in r.deconvolution_results or []:
        if d.notes and "кластер M1" in d.notes:
            m1 = d
            break
    assert m1 is not None
    assert len(m1.overlay_components) == len(m1.components)


def test_legacy_path_still_works_without_overlay():
    """DeconvolutionResult без overlay (через apply_multiplet_deconvolution
    с free NNLS-fit) сохраняет рабочий legacy fallback в PNG / HTML."""
    from gamma.peaks.deconvolve import (
        DeconvolutionResult, MultipletComponent,
    )
    # Создадим минимальный legacy DeconvolutionResult без overlay
    comp = MultipletComponent(
        nuclide="Test", line_E_keV=911.0, library_I_pct=25.0,
        center_channel=911.0, fwhm_channels=15.0,
    )
    d = DeconvolutionResult(
        components=(comp,), areas=(1000.0,), area_uncertainties=(30.0,),
        continuum_params=(100.0, -0.1), continuum_model="linear",
        chi2_per_dof=1.0, n_dof=10, roi_low_ch=850, roi_high_ch=970,
        gross_counts=2000.0, converged=True, method="lstsq",
    )
    # overlay поля должны быть None
    assert d.overlay_E_keV is None
    assert d.overlay_data is None
    assert d.overlay_continuum is None
    assert d.overlay_total is None
    assert d.overlay_components is None
