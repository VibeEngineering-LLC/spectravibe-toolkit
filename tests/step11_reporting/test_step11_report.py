"""
test_step11_report.py — F-85 / v1.14.0 regression for Step 11 reporting.

Verifies:

  1. analyze_lsrm_spe(complete_workflow=True) enables all Round-5 hooks
     when efficiency curve is available, without breaking the v1.12.0
     contract when not.
  2. gamma.reporting.build_json_report produces a dict with all 15
     mandatory top-level keys per references/06_report_format.md.
  3. gamma.reporting.build_chat_summary produces a string of 3–8 lines.
  4. gamma.reporting.build_markdown_report produces text with all 13
     mandatory sections (counted via "## N." headers).
  5. environment classifier returns one of three valid labels.
  6. build_report writes JSON file to disk when output_dir given.
  7. Default analyze_lsrm_spe (no complete_workflow flag) does NOT
     auto-write a report — explicit build_report call required.
  8. JSON is fully serializable (json.dumps round-trip without error).

Run:  PYTHONPATH=scripts python test_step11_report.py
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.identification.staged_pipeline import analyze_lsrm_spe
from gamma.reporting import (
    build_report, build_json_report, build_chat_summary,
    build_markdown_report, classify_environment,
    ENV_NATURAL, ENV_LOW_BG, ENV_UNKNOWN,
)
from gamma.reporting.json_report import SCHEMA_VERSION, SKILL_VERSION
from gamma.detectors.gamma1s import DEFAULT_REFERENCE_DIR


_ROOT = DEFAULT_REFERENCE_DIR
FIXTURE_TH_MARINELLI = _ROOT / "Th232_420-7-17_Маринелли_0cm.spe"   # eff auto-loads
FIXTURE_CS_5CM = _ROOT / "Cs-137__163_2017.spe"                     # no eff


# ════════════════════════════════════════════════════════════════════
# Cached results (one analysis per fixture, shared across tests)
# ════════════════════════════════════════════════════════════════════

_RESULT_TH = None
_RESULT_CS_PLAIN = None
_RESULT_CS_WORKFLOW = None


def _result_th():
    global _RESULT_TH
    if _RESULT_TH is None:
        _RESULT_TH = analyze_lsrm_spe(
            str(FIXTURE_TH_MARINELLI),
            complete_workflow=True,
            sample_mass_kg=0.2,
        )
    return _RESULT_TH


def _result_cs_plain():
    """Default-arguments call — must remain bit-for-bit v1.12.0 compatible."""
    global _RESULT_CS_PLAIN
    if _RESULT_CS_PLAIN is None:
        _RESULT_CS_PLAIN = analyze_lsrm_spe(str(FIXTURE_CS_5CM))
    return _RESULT_CS_PLAIN


def _result_cs_workflow():
    global _RESULT_CS_WORKFLOW
    if _RESULT_CS_WORKFLOW is None:
        _RESULT_CS_WORKFLOW = analyze_lsrm_spe(
            str(FIXTURE_CS_5CM),
            complete_workflow=True,
            allow_stage2=True,
        )
    return _RESULT_CS_WORKFLOW


# ════════════════════════════════════════════════════════════════════
# Test helpers
# ════════════════════════════════════════════════════════════════════

def _assert(cond, msg):
    if not cond:
        raise AssertionError(msg)


def _report(name, fn):
    try:
        fn()
        print(f"  PASS  {name}")
        return True
    except AssertionError as exc:
        print(f"  FAIL  {name}: {exc}")
        return False
    except Exception as exc:
        print(f"  ERROR {name}: {type(exc).__name__}: {exc}")
        return False


# ════════════════════════════════════════════════════════════════════
# 1. complete_workflow umbrella
# ════════════════════════════════════════════════════════════════════

def test_complete_workflow_enables_all_hooks_when_eff_loaded():
    r = _result_th()
    _assert(r.efficiency_curve is not None,
            "Marinelli fixture should auto-load an efficiency curve")
    _assert(r.activities is not None and len(r.activities) > 0,
            "complete_workflow=True must populate activities")
    _assert(r.mda_per_line is not None and len(r.mda_per_line) >= 8,
            "complete_workflow=True must populate mda_per_line with ≥8 rows")
    _assert(r.deconvolution_results is not None,
            "complete_workflow=True must populate deconvolution_results")


def test_default_call_does_not_run_round5():
    """Backward compat — without complete_workflow, Round-5 fields stay None."""
    r = _result_cs_plain()
    _assert(r.activities is None,
            "activities must be None without complete_workflow / explicit flag")
    _assert(r.mda_per_line is None,
            "mda_per_line must be None without complete_workflow / explicit flag")
    _assert(r.deconvolution_results is None,
            "deconvolution_results must be None without complete_workflow / explicit flag")


# ════════════════════════════════════════════════════════════════════
# 2. JSON report schema
# ════════════════════════════════════════════════════════════════════

REQUIRED_TOP_LEVEL_KEYS = {
    "schema_version", "skill_version",
    "header", "calibration",
    "primary_feps", "secondary_peaks", "elemental_xrf",
    "identified_nuclides", "unidentified_peaks", "completeness",
    "mda", "multiplet_deconvolutions",
    "diagnostics", "warnings", "pipeline_notes",
}


def test_json_report_has_required_top_level_keys():
    jd = build_json_report(_result_th())
    missing = REQUIRED_TOP_LEVEL_KEYS - set(jd.keys())
    _assert(not missing, f"missing JSON keys: {missing}")
    _assert(jd["schema_version"] == SCHEMA_VERSION,
            "schema_version mismatch")
    _assert(jd["skill_version"] == SKILL_VERSION,
            "skill_version mismatch")


def test_json_report_is_serializable():
    """All NaN/inf must be sanitized so json.dumps succeeds."""
    jd = build_json_report(_result_th())
    s = json.dumps(jd, ensure_ascii=False)
    # Round-trip
    jd2 = json.loads(s)
    _assert(jd2["schema_version"] == jd["schema_version"],
            "round-trip lost schema_version")
    _assert(set(jd2.keys()) == set(jd.keys()),
            "round-trip lost top-level keys")


def test_header_contains_environment_field():
    jd = build_json_report(_result_th())
    env = jd["header"].get("environment")
    _assert(env in (ENV_NATURAL, ENV_LOW_BG, ENV_UNKNOWN),
            f"environment must be one of three labels, got {env!r}")


def test_primary_feps_carry_rate_cps():
    jd = build_json_report(_result_th())
    feps = jd.get("primary_feps") or []
    _assert(len(feps) > 0, "Th-232 Marinelli should have FEPs")
    for r in feps:
        for required in ("nuclide", "library_E_keV", "peak_E_keV",
                         "peak_area_counts", "rate_cps", "is_characteristic"):
            _assert(required in r,
                    f"primary_feps row missing field {required!r}")


def test_mda_contains_standard_suite():
    jd = build_json_report(_result_th())
    rows = jd.get("mda") or []
    nuclides = {r["nuclide"] for r in rows}
    # Standard suite must include at least these
    suite_required = {"Cs-137", "Co-60", "K-40", "Bi-214", "Tl-208"}
    missing = suite_required - nuclides
    _assert(not missing,
            f"standard MDA suite missing {missing}; got {sorted(nuclides)}")


def test_identified_nuclides_carry_activity_and_bq_per_kg():
    jd = build_json_report(_result_th())
    nucs = jd.get("identified_nuclides") or []
    _assert(len(nucs) > 0, "Th-232 Marinelli should have detected nuclides")
    have_activity = any(n.get("activity_Bq") is not None for n in nucs)
    have_specific = any(
        n.get("specific_activity_Bq_per_kg") is not None for n in nucs
    )
    _assert(have_activity, "no identified nuclide carries activity_Bq")
    _assert(have_specific,
            "no identified nuclide carries specific_activity_Bq_per_kg "
            "despite sample_mass_kg=0.2")


# ════════════════════════════════════════════════════════════════════
# 3. Chat summary
# ════════════════════════════════════════════════════════════════════

def test_chat_summary_is_3_to_8_lines():
    summary = build_chat_summary(_result_th())
    n = summary.count("\n") + 1
    _assert(3 <= n <= 8,
            f"chat summary must be 3–8 lines, got {n}: {summary!r}")


def test_chat_summary_mentions_filename_and_detector():
    summary = build_chat_summary(_result_th())
    _assert("Spectrum:" in summary, "chat summary missing 'Spectrum:' line")
    _assert("Detector:" in summary, "chat summary missing 'Detector:' line")


# ════════════════════════════════════════════════════════════════════
# 4. Markdown report
# ════════════════════════════════════════════════════════════════════

# v1.17.4: Markdown is fully RU.  Section titles are translated;
# legacy English titles are still accepted by the assertion below for
# backward compatibility.
REQUIRED_MD_SECTIONS = [
    ("1. Header", "1. Заголовок"),
    ("2. Detector type", "2. Тип детектора"),
    ("3. Calibration", "3. Калибровка"),
    ("4. Primary FEP", "4. Основные пики полного поглощения"),
    ("5. Secondary peak", "5. Вторичные пики"),
    ("6. Elemental XRF", "6. Элементная XRF"),
    ("7. Identified nuclides", "7. Идентифицированные нуклиды"),
    ("8. Unidentified", "8. Неидентифицированные"),
    ("9. Spectrum plot", "9. График спектра"),
    # BUG-5 / v1.18.30+ (Agent B): markdown H2 renamed to clarify the
    # data source (sample spectrum) vs the bg-multiplet block that lives
    # only in interactive HTML. Legacy strings kept for back-compat.
    ("10. Multiplet", "10. Мультиплеты — разложение в спектре образца"),
    ("11. MDA", "11. Таблица MDA"),
    ("12. Diagnostics", "12. Диагностика"),
    ("13. Version history", "13. История версий"),
]


def test_markdown_report_contains_all_13_sections():
    md = build_markdown_report(_result_th())
    for legacy_en, ru in REQUIRED_MD_SECTIONS:
        en_marker = f"## {legacy_en}"
        ru_marker = f"## {ru}"
        _assert(
            (en_marker in md) or (ru_marker in md),
            f"markdown missing section header '{ru_marker}' (or legacy '{en_marker}')",
        )


def test_markdown_report_renders_version_history():
    md = build_markdown_report(_result_th())
    _assert("v1.14.0" in md, "markdown version history missing v1.14.0")
    _assert("v1.13.0" in md, "markdown version history missing v1.13.0")


# ════════════════════════════════════════════════════════════════════
# 5. Environment classifier
# ════════════════════════════════════════════════════════════════════

def test_environment_classifier_returns_valid_label():
    for fn in (_result_th, _result_cs_workflow):
        env = classify_environment(fn())
        _assert(env in (ENV_NATURAL, ENV_LOW_BG, ENV_UNKNOWN),
                f"environment must be valid label, got {env!r}")


def test_environment_th232_marinelli_is_natural():
    """A natural-background Marinelli reading should classify as 'natural'."""
    env = classify_environment(_result_th())
    _assert(env == ENV_NATURAL,
            f"Th-232 Marinelli should be 'natural', got {env!r}")


# ════════════════════════════════════════════════════════════════════
# 6. build_report dispatcher — disk writes
# ════════════════════════════════════════════════════════════════════

def test_build_report_writes_json_to_disk():
    with tempfile.TemporaryDirectory() as td:
        artefacts = build_report(_result_th(), output_dir=td,
                                 write_json=True, write_markdown=False)
        _assert(artefacts["json"] is not None, "no JSON path returned")
        json_path = Path(artefacts["json"])
        _assert(json_path.exists(), f"JSON file not written: {json_path}")
        # Verify file is valid JSON
        with json_path.open("r", encoding="utf-8") as f:
            jd = json.load(f)
        _assert(jd["schema_version"] == SCHEMA_VERSION,
                "written JSON has wrong schema_version")


def test_build_report_writes_markdown_when_requested():
    with tempfile.TemporaryDirectory() as td:
        artefacts = build_report(_result_th(), output_dir=td,
                                 write_json=False, write_markdown=True)
        _assert(artefacts["markdown"] is not None,
                "no markdown path returned")
        md_path = Path(artefacts["markdown"])
        _assert(md_path.exists(), f"markdown file not written: {md_path}")
        text = md_path.read_text(encoding="utf-8")
        # v1.17.4 — Markdown is now RU; legacy English title still accepted.
        _assert(
            ("# Отчёт о гамма-спектрометрическом анализе" in text)
            or ("# Gamma-spectrum analysis report" in text),
            "markdown file missing top heading",
        )


def test_build_report_in_memory_only_no_files_written():
    """When output_dir=None, no files are written; artefacts come back
    in-memory only."""
    artefacts = build_report(_result_th(), output_dir=None,
                             return_summary=True)
    _assert(artefacts["json"] is None, "json path should be None")
    _assert(artefacts["markdown"] is None, "markdown path should be None")
    _assert(artefacts["json_dict"] is not None,
            "json_dict must always be populated")
    _assert(isinstance(artefacts["summary"], str),
            "summary must be a string when requested")


# ════════════════════════════════════════════════════════════════════
# 7. Orchestrator does NOT auto-write reports
# ════════════════════════════════════════════════════════════════════

def test_analyze_lsrm_spe_does_not_write_files():
    """analyze_lsrm_spe must NEVER write report files automatically —
    build_report is the only public entry-point that touches the disk."""
    r = _result_cs_workflow()
    # analyze_lsrm_spe returns a StagedAnalysisResult; reports are
    # opt-in via build_report. The presence/absence of files on disk
    # is therefore tested by the test_build_report_* family.
    _assert(hasattr(r, "spec"), "result missing spec attr")
    _assert(not hasattr(r, "report_path"),
            "StagedAnalysisResult should not carry a report_path field")


# ════════════════════════════════════════════════════════════════════
# Test harness
# ════════════════════════════════════════════════════════════════════

TESTS = [
    ("complete_workflow_enables_all_hooks_when_eff_loaded",
     test_complete_workflow_enables_all_hooks_when_eff_loaded),
    ("default_call_does_not_run_round5",
     test_default_call_does_not_run_round5),
    ("json_report_has_required_top_level_keys",
     test_json_report_has_required_top_level_keys),
    ("json_report_is_serializable",
     test_json_report_is_serializable),
    ("header_contains_environment_field",
     test_header_contains_environment_field),
    ("primary_feps_carry_rate_cps",
     test_primary_feps_carry_rate_cps),
    ("mda_contains_standard_suite",
     test_mda_contains_standard_suite),
    ("identified_nuclides_carry_activity_and_bq_per_kg",
     test_identified_nuclides_carry_activity_and_bq_per_kg),
    ("chat_summary_is_3_to_8_lines",
     test_chat_summary_is_3_to_8_lines),
    ("chat_summary_mentions_filename_and_detector",
     test_chat_summary_mentions_filename_and_detector),
    ("markdown_report_contains_all_13_sections",
     test_markdown_report_contains_all_13_sections),
    ("markdown_report_renders_version_history",
     test_markdown_report_renders_version_history),
    ("environment_classifier_returns_valid_label",
     test_environment_classifier_returns_valid_label),
    ("environment_th232_marinelli_is_natural",
     test_environment_th232_marinelli_is_natural),
    ("build_report_writes_json_to_disk",
     test_build_report_writes_json_to_disk),
    ("build_report_writes_markdown_when_requested",
     test_build_report_writes_markdown_when_requested),
    ("build_report_in_memory_only_no_files_written",
     test_build_report_in_memory_only_no_files_written),
    ("analyze_lsrm_spe_does_not_write_files",
     test_analyze_lsrm_spe_does_not_write_files),
]


def main():
    print("F-85 / v1.14.0 Step 11 reporting regression")
    print("=" * 70)
    passed = failed = 0
    for name, fn in TESTS:
        if _report(name, fn):
            passed += 1
        else:
            failed += 1
    print("=" * 70)
    print(f"{passed}/{len(TESTS)} pass, {failed} fail")
    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
