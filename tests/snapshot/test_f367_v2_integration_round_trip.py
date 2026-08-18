# -*- coding: utf-8 -*-
"""F-367 / v1.18.24.2 — V2 production-pipeline integration guard.

Контракт: `analyze_and_report_v2` через `v2_peak_search_patched()` context
manager даёт sample-отчёт **строго идентичный** production sample-отчёту по
составу артефактов и JSON-схеме; отличается ТОЛЬКО результатами peak search
(V2 dual-method = Mariscotti ∪ matched filter), всё downstream — без изменений.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))


def test_F367_v2_mariscotti_replacement_signature_compat():
    """V2 replacement принимает все kwargs Mariscotti без TypeError."""
    from gamma.experimental.v2_integration import v2_mariscotti_replacement
    import numpy as np
    counts = np.zeros(500, dtype=float)
    # синтетический пик на канале 250
    for i, val in enumerate([5, 15, 35, 60, 80, 60, 35, 15, 5]):
        counts[246 + i] = val
    # Вызвать со всеми kwargs production-Mariscotti
    result = v2_mariscotti_replacement(
        counts=counts,
        fwhm_channels=lambda ch: 4.0,
        sigma_threshold=3.0,
        min_separation_factor=1.0,
        edge_margin=10,
        band_ratio=1.2,
        filter_narrow_peaks=False,
        min_fwhm_ratio=0.3,
    )
    assert isinstance(result, list), "V2 replacement должен вернуть list"
    # Каждый элемент — FoundPeak с production-полями
    from gamma.peaks.search import FoundPeak
    for p in result:
        assert isinstance(p, FoundPeak), (
            f"expected FoundPeak, got {type(p).__name__}"
        )
        # Production-fields обязаны быть заполнены
        assert p.channel >= 0
        assert p.fwhm_channels > 0
        assert p.significance >= 0


def test_F367_context_manager_restores_on_exception():
    """`v2_peak_search_patched()` гарантированно восстанавливает оригинал
    даже при исключении внутри блока."""
    from gamma.experimental.v2_integration import v2_peak_search_patched
    from gamma.identification import staged_pipeline as _sp
    original = _sp.mariscotti_search

    with pytest.raises(RuntimeError, match="forced"):
        with v2_peak_search_patched():
            # внутри патча, attribute заменён
            assert _sp.mariscotti_search is not original
            raise RuntimeError("forced")
    # после exception — оригинал восстановлен
    assert _sp.mariscotti_search is original, (
        "context manager не восстановил оригинал после exception"
    )


SAMPLE_DIR = REPO / "demo_reports" / "v1_18_24_th232_full" / "sample"
V2_DIR = REPO / "demo_reports" / "v1_18_24_th232_full" / "sample_v2"


def test_F367_th232_demo_v2_sample_strict_composition():
    """V2 sample-отчёт должен содержать ТОТ ЖЕ набор файлов что production
    sample (JSON + MD + HTML + Technical PDF + spectrum.png + ≥1 multiplet
    PNG + 2 BecqMoni XML образца и фона)."""
    if not SAMPLE_DIR.is_dir() or not V2_DIR.is_dir():
        pytest.skip("demo dirs not present")

    def file_categories(d: Path) -> dict:
        files = list(d.iterdir())
        return {
            "json": sorted([f.name for f in files if f.suffix == ".json"]),
            "md": sorted([f.name for f in files if f.suffix == ".md"]),
            "html": sorted([f.name for f in files if f.suffix == ".html"]),
            "pdf": sorted([f.name for f in files if f.suffix == ".pdf"]),
            "xml": sorted([f.name for f in files if f.suffix == ".xml"]),
        }

    prod_cat = file_categories(SAMPLE_DIR)
    v2_cat = file_categories(V2_DIR)

    # Каждая категория должна иметь одинаковое количество файлов
    for cat in ("json", "md", "html", "pdf", "xml"):
        assert len(prod_cat[cat]) == len(v2_cat[cat]), (
            f"category {cat}: prod={prod_cat[cat]} vs v2={v2_cat[cat]}"
        )
        assert len(prod_cat[cat]) > 0, (
            f"production не содержит файлов категории {cat}; "
            "регенерируйте через `python -m gamma.cli analyze ... --full-report "
            "--export-becqmoni both`"
        )

    # Должно быть ровно 2 BecqMoni XML (sample + bg)
    assert len(prod_cat["xml"]) == 2, (
        f"expected 2 BecqMoni XML (sample+bg), got {prod_cat['xml']}"
    )
    assert len(v2_cat["xml"]) == 2, (
        f"expected 2 BecqMoni XML (sample+bg), got {v2_cat['xml']}"
    )

    # Plots: spectrum + multiplets
    prod_plots = SAMPLE_DIR / "Th232_Маринелли_0cm_plots"
    v2_plots = V2_DIR / "Th232_Маринелли_0cm_plots"
    assert (prod_plots / "spectrum.png").is_file()
    assert (v2_plots / "spectrum.png").is_file()
    prod_multi = list((prod_plots / "multiplets").glob("*.png"))
    v2_multi = list((v2_plots / "multiplets").glob("*.png"))
    assert len(prod_multi) >= 1, "production должен иметь ≥1 multiplet PNG"
    assert len(v2_multi) >= 1, "V2 должен иметь ≥1 multiplet PNG"


def test_F367_th232_demo_v2_json_schema_matches_production():
    """V2 JSON top-level keys строго равны production keys (различия
    только в data, не в schema)."""
    if not SAMPLE_DIR.is_dir() or not V2_DIR.is_dir():
        pytest.skip("demo dirs not present")
    prod_json = next(SAMPLE_DIR.glob("*_report.json"), None)
    v2_json = next(V2_DIR.glob("*_report.json"), None)
    assert prod_json is not None, "production JSON отсутствует"
    assert v2_json is not None, "V2 JSON отсутствует"
    prod_data = json.loads(prod_json.read_text(encoding="utf-8"))
    v2_data = json.loads(v2_json.read_text(encoding="utf-8"))
    prod_keys = set(prod_data.keys())
    v2_keys = set(v2_data.keys())
    assert prod_keys == v2_keys, (
        f"JSON-schema mismatch: only in prod={prod_keys-v2_keys}, "
        f"only in v2={v2_keys-prod_keys}"
    )
    # Версии должны совпадать
    assert prod_data["skill_version"] == v2_data["skill_version"], (
        f"version mismatch: prod={prod_data['skill_version']} "
        f"v2={v2_data['skill_version']}"
    )
    # V2 ожидаемо имеет >= peaks чем production (dual-method finds more)
    assert (
        len(v2_data["primary_feps"]) >= len(prod_data["primary_feps"])
    ), (
        f"V2 dual-search должен найти ≥ production-Mariscotti: "
        f"prod={len(prod_data['primary_feps'])} v2={len(v2_data['primary_feps'])}"
    )
