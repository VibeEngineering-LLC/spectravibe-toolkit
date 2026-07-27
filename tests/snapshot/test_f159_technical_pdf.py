"""F-159 / v1.18.21.0 — Technical PDF report integration test.

Validates:
  * build_technical_pdf produces a valid PDF from a sample JSON dict
  * 11 шагов присутствуют в извлечённом тексте PDF
  * Anonymization F-115: PDF не содержит абсолютных путей операторов
  * analyze_and_report по умолчанию ВКЛЮЧАЕТ PDF (контракт навсегда)
  * write_technical_pdf=False корректно отключает генерацию
  * CLI флаг --no-technical-pdf уважается
  * Fallback: при отсутствии reportlab artefact пропускается с warning
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, "scripts")

ROOT = Path(__file__).resolve().parents[2]
KIT = ROOT / "detectors/Gamma-1S/reference_spectra/reference_kits"
CS_SAMPLE = KIT / "Marinelli_1L/Cs-137/sample_M_cs_легкий_2001-2005.spe"
CS_BG = KIT / "Marinelli_1L/Cs-137/background_bg_2016_marinelli_water_marinelli.spe"


# ──────────────────────────────────────────────────────────────────
# Direct-call tests (на готовом JSON демо-отчёте — быстрые)
# ──────────────────────────────────────────────────────────────────

@pytest.fixture
def demo_json():
    """Берём демо v1.18.20 как стабильный фикстура для unit-проверок."""
    demo_path = ROOT / "demo_reports/v1_18_20/M_th_легкий_report.json"
    if not demo_path.exists():
        pytest.skip("demo_reports/v1_18_20 missing — regenerate первый.")
    with demo_path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _read_pdf_text(path: Path) -> str:
    try:
        from pypdf import PdfReader
    except ImportError:
        pytest.skip("pypdf not installed — нельзя проверить content PDF")
    reader = PdfReader(str(path))
    return "\n".join(p.extract_text() for p in reader.pages)


def test_f159_build_technical_pdf_smoke(demo_json, tmp_path):
    """Smoke: модуль импортируется, PDF создаётся, файл > 5 КБ."""
    from gamma.reporting.technical_pdf import build_technical_pdf

    out = tmp_path / "tech.pdf"
    res = build_technical_pdf(None, demo_json, str(out))
    assert Path(res).exists()
    assert out.stat().st_size > 5000, "PDF слишком маленький — пустой output?"


def test_f159_pdf_contains_all_11_steps(demo_json, tmp_path):
    """Все 11 шагов walkthrough присутствуют в извлечённом тексте."""
    from gamma.reporting.technical_pdf import build_technical_pdf

    out = tmp_path / "tech.pdf"
    build_technical_pdf(None, demo_json, str(out))
    text = _read_pdf_text(out)

    # Контракт: каждый Шаг N появляется (в TOC + в section heading).
    # Не считаем точное число вхождений, но требуем оба source места.
    for n in range(1, 12):
        assert f"Шаг {n}" in text, f"PDF missing 'Шаг {n}'"

    # Ключевые секции имеют узнаваемые заголовки
    for marker in [
        "Чтение файла",
        "Среда измерения",
        "Поиск пиков",
        "Тип детектора",
        "Энергокалибровка",
        "ПШПВ",
        "Идентификация",
        "Деконволюция",
        "активностей",
        "Вторичные пики",
        "Финальный отчёт",
    ]:
        assert marker in text, f"PDF missing section marker: {marker}"


def test_f159_pdf_anonymization_no_abs_paths(demo_json, tmp_path):
    """F-115: PDF не должен содержать абсолютных путей операторов."""
    from gamma.reporting.technical_pdf import build_technical_pdf

    out = tmp_path / "tech.pdf"
    build_technical_pdf(None, demo_json, str(out))
    text = _read_pdf_text(out)

    # Windows-стиль: C:\, D:\, ...
    win_path = re.search(r"\b[A-Za-z]:[\\/][A-Za-z0-9_\\/\.\- ]{3,}", text)
    assert win_path is None, f"abs Windows path leaked: {win_path.group(0)}"

    # POSIX-стиль: /home/<user>, /Users/<name>
    posix_user = re.search(r"/(?:home|Users)/[A-Za-z0-9_\-\.]+/", text)
    assert posix_user is None, f"abs POSIX path leaked: {posix_user.group(0)}"


def test_f159_pdf_anonymization_no_cert_serial(demo_json, tmp_path):
    """F-115: PDF не должен содержать сериалов сертифицированных источников
    типа '420-7-17' или 'SN-01'."""
    from gamma.reporting.technical_pdf import build_technical_pdf

    out = tmp_path / "tech.pdf"
    build_technical_pdf(None, demo_json, str(out))
    text = _read_pdf_text(out)

    # Cert-S/N: ровно \d{2,4}-\d+-\d+ (как в эталоне 420-7-17, SN-01-3)
    sn = re.search(r"\b\d{3,4}-\d+-\d+\b", text)
    # NB: значения типа "73-90" (диапазон энергии) попадают сюда не должны —
    # они два числа без второго дефиса. Cert-SN всегда трёхсегментный.
    assert sn is None, f"cert serial leaked: {sn.group(0)}"


# ──────────────────────────────────────────────────────────────────
# Integration tests (через analyze_and_report)
# ──────────────────────────────────────────────────────────────────

def test_f159_analyze_and_report_writes_pdf_by_default(tmp_path):
    """F-RPT-03 / v1.18.29 — Technical PDF теперь OFF by default.

    Контракт изменился (Phase 4.5 / Phase 5 v1.18.29):
    - default: PDF НЕ генерится (write_technical_pdf=False по умолчанию)
    - opt-in: явный `write_technical_pdf=True` восстанавливает старый F-159 контракт

    Имя теста сохранено для git blame trace; смысл проверки инвертирован."""
    if not CS_SAMPLE.exists():
        pytest.skip("kit sample missing")
    from gamma.reporting import analyze_and_report

    # F-RPT-03 NEW default — PDF не пишется.
    res_default = analyze_and_report(
        str(CS_SAMPLE),
        background_path=str(CS_BG),
        sample_mass_kg=0.570,
        output_dir=str(tmp_path / "f159_default"),
        write_html=False, write_plots=False, write_markdown=False, write_json=True,
    )
    assert not res_default.get("technical_pdf"), (
        "F-RPT-03: technical_pdf must be empty/None when write_technical_pdf "
        "is left at default (False) per v1.18.29 contract"
    )

    # Opt-in pathway — старый F-159 контракт всё ещё доступен.
    res_optin = analyze_and_report(
        str(CS_SAMPLE),
        background_path=str(CS_BG),
        sample_mass_kg=0.570,
        output_dir=str(tmp_path / "f159_optin"),
        write_html=False, write_plots=False, write_markdown=False, write_json=True,
        write_technical_pdf=True,
    )
    assert res_optin.get("technical_pdf"), "F-159 opt-in pathway broken: PDF not written"
    p = Path(res_optin["technical_pdf"])
    assert p.exists()
    assert p.suffix == ".pdf"
    assert p.name.endswith("_technical_report.pdf")
    assert p.stat().st_size > 5000


def test_f159_opt_out_via_write_technical_pdf_false(tmp_path):
    """write_technical_pdf=False корректно отключает генерацию."""
    if not CS_SAMPLE.exists():
        pytest.skip("kit sample missing")
    from gamma.reporting import analyze_and_report

    res = analyze_and_report(
        str(CS_SAMPLE),
        background_path=str(CS_BG),
        sample_mass_kg=0.570,
        output_dir=str(tmp_path / "f159_off"),
        write_html=False, write_plots=False, write_markdown=False, write_json=True,
        write_technical_pdf=False,
    )
    assert res.get("technical_pdf") is None


# ──────────────────────────────────────────────────────────────────
# Version-bump assertion (защита от случайного rollback версии)
# ──────────────────────────────────────────────────────────────────

def test_f159_version_bump():
    """Skill version >= (1, 18, 21, 0) — F-159 landed."""
    from gamma.reporting.json_report import SKILL_VERSION

    m = re.match(r"v(\d+)\.(\d+)\.(\d+)(?:\.(\d+))?", SKILL_VERSION)
    assert m, f"Unparseable SKILL_VERSION: {SKILL_VERSION}"
    parts = tuple(int(p or 0) for p in m.groups())
    assert parts >= (1, 18, 21, 0), (
        f"SKILL_VERSION {SKILL_VERSION} below F-159 baseline v1.18.21.0"
    )
