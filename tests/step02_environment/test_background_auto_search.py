"""F-131 / v1.17.7 — Auto-background search heuristic."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.io.background_search import (
    BackgroundCandidate,
    BG_FILENAME_MARKERS,
    BG_TYPE_MARKERS,
    find_background_candidates,
    render_suggestion_note,
    render_applied_note,
    _filename_looks_like_bg,
    _type_field_says_bg,
    _geometry_compatible,
    _detector_compatible,
    _recency_score,
)


FIXTURE_DIR = (Path(__file__).parent.parent.parent / "detectors" / "Gamma-1S"
               / "reference_spectra" / "archive")
TH232_FIXTURE = FIXTURE_DIR / "Th232_420-7-17_Маринелли_0cm.spe"


# ─── Unit tests ──────────────────────────────────────────────────────────

def test_filename_markers_detected():
    """Файлы с маркерами фона в имени распознаются."""
    assert _filename_looks_like_bg("bg_2024_marinelli.spe")
    assert _filename_looks_like_bg("BKG_open_lid.spe")
    assert _filename_looks_like_bg("Фон_закр_кр_вода_01.spe")
    assert _filename_looks_like_bg("background_v1.spe")


def test_filename_non_bg_rejected():
    """Обычные sample-файлы НЕ помечаются как BG."""
    assert not _filename_looks_like_bg("Th232_420-7-17.spe")
    assert not _filename_looks_like_bg("Cs137_фикстура.spe")
    assert not _filename_looks_like_bg("K40_Маринелли_0cm.spe")


def test_type_field_bg_markers():
    """TYPE=Фон/Background → True."""
    assert _type_field_says_bg("Фон")
    assert _type_field_says_bg("Background")
    assert _type_field_says_bg("BG")
    assert not _type_field_says_bg("Калибровка")
    assert not _type_field_says_bg("")


def test_geometry_exact_match():
    """Одинаковая геометрия → score=1.0."""
    ok, score = _geometry_compatible("Маринелли", "Маринелли")
    assert ok
    assert score == 1.0


def test_geometry_partial_match():
    """Точечная BG для Marinelli sample допустима с пониженной уверенностью."""
    ok, score = _geometry_compatible("Маринелли", "Точечная-5см")
    assert ok
    assert 0.4 < score < 0.7


def test_geometry_incompatible():
    """Несовместимые геометрии → отброс."""
    ok, _ = _geometry_compatible("Marinelli", "Дента-120мл")
    assert not ok


def test_detector_exact_match():
    """Точное совпадение детектора → score=1.0."""
    ok, score = _detector_compatible("Гамма-1С", "Гамма-1С")
    assert ok
    assert score == 1.0


def test_detector_partial_match():
    """Partial-совпадение (общий токен) → score 0.7-0.85.

    F-313 / v1.18.12: одинаковый Gamma-1S alias в обоих → score 0.85
    (раньше 0.7). Разные aliases из набора → 0.7. Расширено для
    canonical alias normalization (см. KNOWN_AND_FIXED_ISSUES F-313).
    """
    ok, score = _detector_compatible(
        "Гамма-1С NaI 63х63", "Гамма-1С USB SN-01",
    )
    assert ok
    assert 0.5 < score <= 0.9     # was 0.8; F-313 расширил scoring до 0.85


def test_F313_detector_alias_match_vendor_tokens():
    """F-313: vendor-токены УДС-ГЦ / БДЭГ / Колибри matched as Gamma-1S alias."""
    # Sample detector_id из M_cs (LSRM SpectraLine standard)
    # Background detector_id из bg_2016_marinelli (Aspect УДС-ГЦ vendor)
    ok, score = _detector_compatible(
        "Гамма-1С", "УДС-ГЦ-63х63-USB №SN-01",
    )
    assert ok, "Gamma-1С vs УДС-ГЦ — canonical alias match failed"
    assert score >= 0.5


def test_detector_unknown_compatible():
    """Если у sample/кандидата детектор пуст — допустимо но с малым score."""
    ok, score = _detector_compatible("", "Гамма-1С")
    assert ok
    assert score < 0.5


def test_recency_score_monotone():
    """Чем ближе по дате, тем выше score."""
    s_close = _recency_score(0.0, 90.0)
    s_mid = _recency_score(45.0, 90.0)
    s_far = _recency_score(89.0, 90.0)
    assert s_close > s_mid > s_far
    assert s_close == 1.0


def test_recency_score_out_of_range():
    """Δt > max_days → score 0 (фактически кандидат отсеется выше)."""
    assert _recency_score(150.0, 90.0) == 0.0


def test_recency_score_unknown_date():
    """Неизвестная дата → нейтральный bonus."""
    s = _recency_score(None, 90.0)
    assert 0.2 < s < 0.5


# ─── E2E tests на реальной фикстуре ──────────────────────────────────

def test_find_candidates_for_th232_fixture():
    """Th-232 Marinelli фикстура должна найти ≥1 BG-кандидата
    в той же папке (есть «Фон_закр_кр_вода_*.spe») или в
    data/averaged_backgrounds/."""
    if not TH232_FIXTURE.exists():
        pytest.skip(f"fixture missing: {TH232_FIXTURE}")
    from gamma.io.readers import read_spectrum
    spec = read_spectrum(str(TH232_FIXTURE))
    candidates = find_background_candidates(spec, str(TH232_FIXTURE))
    assert len(candidates) >= 1
    # Лучший кандидат должен иметь высокую уверенность (≥3.0)
    best = candidates[0]
    assert best.confidence_score >= 3.0
    # И быть в той же папке либо averaged_backgrounds
    assert ("Фон" in best.path.parent.name
            or "averaged_backgrounds" in best.path.parts
            or _filename_looks_like_bg(best.path.name))


def test_candidate_dict_serializable():
    """BackgroundCandidate.to_dict даёт JSON-serializable словарь."""
    if not TH232_FIXTURE.exists():
        pytest.skip(f"fixture missing: {TH232_FIXTURE}")
    from gamma.io.readers import read_spectrum
    spec = read_spectrum(str(TH232_FIXTURE))
    candidates = find_background_candidates(spec, str(TH232_FIXTURE))
    assert candidates
    d = candidates[0].to_dict()
    import json
    # Должен сериализоваться без ошибок
    s = json.dumps(d, ensure_ascii=False)
    assert "confidence_score" in s
    assert "filename" in s


def test_render_suggestion_note_ru():
    """render_suggestion_note возвращает RU-нарратив с маркером F-131."""
    if not TH232_FIXTURE.exists():
        pytest.skip(f"fixture missing: {TH232_FIXTURE}")
    from gamma.io.readers import read_spectrum
    spec = read_spectrum(str(TH232_FIXTURE))
    candidates = find_background_candidates(spec, str(TH232_FIXTURE))
    note = render_suggestion_note(candidates[0])
    assert "F-131" in note
    assert "background-path" in note or "background-auto" in note
    # RU фразы
    assert "предложен" in note or "фоновый" in note


# ─── Pipeline integration ─────────────────────────────────────────────

def test_pipeline_suggest_mode_does_not_apply():
    """suggest-режим не вычитает фон, только предлагает."""
    if not TH232_FIXTURE.exists():
        pytest.skip(f"fixture missing: {TH232_FIXTURE}")
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    r = analyze_lsrm_spe(
        str(TH232_FIXTURE), complete_workflow=True,
        background_auto="suggest",
    )
    assert r.auto_background_mode == "suggest"
    assert r.auto_background_applied_path is None
    assert r.background_status != "auto_resolved_from_directory"
    # Кандидаты найдены
    assert r.auto_background_candidates
    # И F-131 в notes
    notes_str = " ".join(r.notes or [])
    assert "F-131" in notes_str


def test_pipeline_apply_mode_subtracts_background():
    """apply-режим автоматически вычитает лучшего кандидата."""
    if not TH232_FIXTURE.exists():
        pytest.skip(f"fixture missing: {TH232_FIXTURE}")
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    r = analyze_lsrm_spe(
        str(TH232_FIXTURE), complete_workflow=True,
        background_auto="apply",
    )
    assert r.auto_background_mode == "apply"
    assert r.auto_background_applied_path is not None
    assert r.background_status == "auto_resolved_from_directory"
    assert r.background_subtraction is not None


def test_pipeline_off_mode_skips_search():
    """off-режим вообще не запускает поиск."""
    if not TH232_FIXTURE.exists():
        pytest.skip(f"fixture missing: {TH232_FIXTURE}")
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    r = analyze_lsrm_spe(
        str(TH232_FIXTURE), complete_workflow=True,
        background_auto="off",
    )
    assert r.auto_background_mode == "off"
    assert r.auto_background_candidates is None
    assert r.background_status == "absent_no_subtraction"


def test_explicit_background_path_overrides_auto():
    """Если --background-path задан явно, авто-поиск не должен
    подменить выбор пользователя."""
    if not TH232_FIXTURE.exists():
        pytest.skip(f"fixture missing: {TH232_FIXTURE}")
    bg_path = (FIXTURE_DIR / "Фон_закр_кр_вода_01.spe")
    if not bg_path.exists():
        pytest.skip(f"bg fixture missing: {bg_path}")
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    r = analyze_lsrm_spe(
        str(TH232_FIXTURE), complete_workflow=True,
        background_path=str(bg_path),
        background_auto="apply",
    )
    # background_status — обычный subtracted_from_external_file,
    # НЕ auto_resolved_from_directory.
    assert r.background_status == "subtracted_from_external_file"
    # F-325 / v1.18.18.1 — applied_path теперь записывается ВСЕГДА когда
    # фон был вычтен (auto или explicit). Это нужно чтобы reports
    # surface'или имя bg-файла. Контракт «explicit overrides auto»
    # проверяется через background_status (не "auto_resolved_..."),
    # а не через None-ность applied_path.
    assert r.auto_background_applied_path == str(bg_path)


def test_json_report_includes_auto_background_block():
    """JSON-отчёт должен содержать diagnostics.auto_background_search."""
    if not TH232_FIXTURE.exists():
        pytest.skip(f"fixture missing: {TH232_FIXTURE}")
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    from gamma.reporting.json_report import build_json_report
    r = analyze_lsrm_spe(
        str(TH232_FIXTURE), complete_workflow=True,
        background_auto="suggest",
    )
    rep = build_json_report(r)
    block = rep.get("diagnostics", {}).get("auto_background_search")
    assert block is not None
    assert block.get("mode") == "suggest"
    assert block.get("candidates")


def test_max_days_filter():
    """Кандидаты с |Δt| > max_days отбрасываются."""
    if not TH232_FIXTURE.exists():
        pytest.skip(f"fixture missing: {TH232_FIXTURE}")
    from gamma.io.readers import read_spectrum
    spec = read_spectrum(str(TH232_FIXTURE))
    # max_days=0 — только если ровно в тот же день
    candidates_strict = find_background_candidates(
        spec, str(TH232_FIXTURE), max_days_apart=0,
    )
    candidates_loose = find_background_candidates(
        spec, str(TH232_FIXTURE), max_days_apart=365,
    )
    assert len(candidates_loose) >= len(candidates_strict)


def test_does_not_suggest_sample_as_its_own_background():
    """Sample-файл НЕ должен попасть в список кандидатов фона для самого себя."""
    if not TH232_FIXTURE.exists():
        pytest.skip(f"fixture missing: {TH232_FIXTURE}")
    from gamma.io.readers import read_spectrum
    spec = read_spectrum(str(TH232_FIXTURE))
    candidates = find_background_candidates(spec, str(TH232_FIXTURE))
    sample_resolved = Path(str(TH232_FIXTURE)).resolve()
    for c in candidates:
        assert c.path.resolve() != sample_resolved
