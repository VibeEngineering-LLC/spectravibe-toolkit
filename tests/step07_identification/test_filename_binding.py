"""
test_filename_binding.py — F-89 / v1.15.2 regression.

Covers the four user-reported gaps from the v1.15.1 demo review
(2026-05-29):

  1. F-89a — background-subtraction status MUST be explicitly
     reported in the JSON header, Markdown table, HTML table, and
     chat summary. No silent omission.
  2. F-89b — filename tokens like "Cs137", "Th232", "K40" must be
     canonicalised to standard labels ("Cs-137", "Th-232", "K-40")
     and surfaced as a binding hypothesis.
  3. F-89c/d — Th-232 source: U-chain identifications (Bi-214,
     Pb-214, Pb-210) must be suppressed because on NaI 63×63 the
     609 keV peak is often Tl-208 583 shifted by Compton overlap,
     and the user has bound the spectrum to Th-only.
  4. F-89e — Cs-137 source: Cs-137 MUST appear in the identified
     list (the v1.15.1 "полный провал" bug). The filename hint
     drives the Stage-1 candidate list per SKILL.md §7A.1.

Run:  PYTHONPATH=scripts python test_filename_binding.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.identification.staged_pipeline import analyze_lsrm_spe
from gamma.identification.anchor_ranks import (
    derive_chain_dominance, AnchorMatch, AnchorEntry,
    ChainDominance, U238_PROXY_NUCLIDES, TH232_PROXY_NUCLIDES,
)
from gamma.io.filename_hints import (
    extract_isotope_hints, chains_claimed_by_isotope_hints,
    parse_filename,
)
from gamma.reporting import (
    build_json_report, build_chat_summary,
    build_markdown_report, build_html_report,
)
from gamma.detectors.gamma1s import DEFAULT_REFERENCE_DIR


_ROOT = DEFAULT_REFERENCE_DIR
FIX_TH = _ROOT / "Th232_420-7-17_Маринелли_0cm.spe"
FIX_CS = _ROOT / "Cs137_420-7-14_Маринелли_0cm.spe"
FIX_K  = _ROOT / "K40_420-7-20_Маринелли_0cm.spe"


_R_TH = _R_CS = _R_K = None


def _result_th():
    global _R_TH
    if _R_TH is None:
        _R_TH = analyze_lsrm_spe(str(FIX_TH), complete_workflow=True,
                                sample_mass_kg=0.2)
    return _R_TH


def _result_cs():
    global _R_CS
    if _R_CS is None:
        _R_CS = analyze_lsrm_spe(str(FIX_CS), complete_workflow=True,
                                sample_mass_kg=0.2)
    return _R_CS


def _result_k():
    global _R_K
    if _R_K is None:
        _R_K = analyze_lsrm_spe(str(FIX_K), complete_workflow=True,
                               sample_mass_kg=0.2)
    return _R_K


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
        import traceback
        traceback.print_exc()
        print(f"  ERROR {name}: {type(exc).__name__}: {exc}")
        return False


# ════════════════════════════════════════════════════════════════════
# F-89b — filename isotope extraction
# ════════════════════════════════════════════════════════════════════

def test_extract_isotope_hints_canonical_form():
    """Canonical labels: Cs137 → Cs-137, Th232 → Th-232, etc."""
    cases = [
        ("Cs137_420-7-14_Маринелли_0cm.spe", ["Cs-137"]),
        ("Th232_420-7-17_Маринелли_0cm.spe", ["Th-232"]),
        ("K40_420-7-20_Маринелли_0cm.spe", ["K-40"]),
        ("Co60_check_source.spe", ["Co-60"]),
        ("Am241_test.spe", ["Am-241"]),
        ("Ra226_marinelli.spe", ["Ra-226"]),
        ("Eu152_calib.spe", ["Eu-152"]),
        ("Фон_лаба_2024.spe", []),
        ("noise_spectrum.spe", []),
    ]
    for fn, expected in cases:
        got = extract_isotope_hints(fn)
        _assert(got == expected,
                f"{fn!r}: expected {expected}, got {got}")


def test_chains_claimed_by_isotope_hints():
    """Chain mapping: Th-232 → {"Th-232"}, Bi-214 → {"U-238"}, etc."""
    _assert(chains_claimed_by_isotope_hints(["Th-232"]) == {"Th-232"},
            "Th-232 must claim Th-232 chain")
    _assert(chains_claimed_by_isotope_hints(["Tl-208"]) == {"Th-232"},
            "Tl-208 is Th-232 chain member")
    _assert(chains_claimed_by_isotope_hints(["Bi-214"]) == {"U-238"},
            "Bi-214 is U-238 chain member")
    _assert(chains_claimed_by_isotope_hints(["Cs-137"]) == {"Cs-137"},
            "Cs-137 has no chain — returns its own label")
    _assert(chains_claimed_by_isotope_hints([]) == set(),
            "empty hint → empty set")


def test_parse_filename_carries_hints_and_chains():
    pt = parse_filename("Th232_420-7-17_Маринелли_0cm.spe")
    _assert(pt["isotope_hints"] == ["Th-232"],
            f"isotope_hints wrong: {pt['isotope_hints']}")
    _assert(pt["chains_claimed"] == ["Th-232"],
            f"chains_claimed wrong: {pt['chains_claimed']}")


# ════════════════════════════════════════════════════════════════════
# F-89e — filename hints DRIVE candidate list (the "полный провал" bug)
# ════════════════════════════════════════════════════════════════════

def test_cs137_fixture_identifies_cs137():
    """CRITICAL: Cs-137 must be in the identified list on a Cs-137 fixture.

    Pre-v1.15.2: Cs-137 was in Stage 2 (technogenic); without explicit
    allow_stage2=True, it never entered the candidate list — the
    bright 661 keV peak ended up unidentified. v1.15.2 / F-89e
    auto-adds filename-hinted isotopes to Stage 1.
    """
    r = _result_cs()
    nuclide_names = {n.nuclide for n in r.final_detected}
    _assert("Cs-137" in nuclide_names,
            f"Cs-137 missing from identified list! Got: {nuclide_names}")


def test_k40_fixture_identifies_k40():
    """K-40 fixture → K-40 in identified list."""
    r = _result_k()
    nuclide_names = {n.nuclide for n in r.final_detected}
    _assert("K-40" in nuclide_names,
            f"K-40 missing from identified list! Got: {nuclide_names}")


def test_orchestrator_carries_filename_isotope_hints():
    r_cs = _result_cs()
    _assert(r_cs.filename_isotope_hints == ["Cs-137"],
            f"Cs result hints: {r_cs.filename_isotope_hints}")
    r_th = _result_th()
    _assert(r_th.filename_isotope_hints == ["Th-232"],
            f"Th result hints: {r_th.filename_isotope_hints}")
    r_k = _result_k()
    _assert(r_k.filename_isotope_hints == ["K-40"],
            f"K result hints: {r_k.filename_isotope_hints}")


# ════════════════════════════════════════════════════════════════════
# F-89d — U-chain suppression on Th-only filename
# ════════════════════════════════════════════════════════════════════

def test_th232_fixture_suppresses_u238_chain():
    """Th-232 fixture must NOT identify Bi-214, Pb-214, Pb-210.

    User feedback: "Радия в образце точно нет. Его следы могут быть
    только от фона, но cps фона пренебрежимо мал."
    """
    r = _result_th()
    nuclide_names = {n.nuclide for n in r.final_detected}
    forbidden = {"Bi-214", "Pb-214", "Pb-210", "Ra-226"}
    intersection = nuclide_names & forbidden
    _assert(not intersection,
            f"Th-232 source still has U-chain identifications: "
            f"{intersection}. Identified: {nuclide_names}")


def test_th232_fixture_keeps_th_chain_nuclides():
    """Th-232 fixture must keep Tl-208 / Bi-212 / Ac-228 etc."""
    r = _result_th()
    nuclide_names = {n.nuclide for n in r.final_detected}
    expected_present = {"Tl-208", "Ac-228"}
    missing = expected_present - nuclide_names
    _assert(not missing,
            f"Th-chain nuclides missing: {missing}. Got: {nuclide_names}")


def test_th232_fixture_records_suppression():
    """suppressed_chains populated; U-chain не попадает в final.

    F-89d удаляет U-chain нуклиды post-Stage; F-123 + F-125 (v1.17.6)
    могут предотвратить ложную идентификацию ещё на Stage. Контракт
    тот же: в результате не должно быть U-chain.
    """
    r = _result_th()
    _assert("U-238" in r.chain_dominance.suppressed_chains,
            f"suppressed_chains missing U-238: "
            f"{r.chain_dominance.suppressed_chains}")
    # Либо chain_filtered_out содержит U-chain (классический F-89d путь),
    # либо ни одного U-chain нуклида нет в final_detected (v1.17.6 путь
    # через узкое окно идентификации после F-125 refit).
    u_chain_nuclides = {"Bi-214", "Pb-214", "Pb-210", "Ra-226", "U-238"}
    u_in_final = {n.nuclide for n in r.final_detected} & u_chain_nuclides
    u_in_filtered = set(r.chain_filtered_out) & u_chain_nuclides
    _assert(not u_in_final,
            f"U-chain нуклиды не должны быть в final_detected, "
            f"но обнаружены: {u_in_final}")
    # Допустимы оба пути; ниже регистрируем какой именно сработал.
    _assert(not u_in_final or u_in_filtered,
            f"либо chain_filtered_out содержит U-chain, либо они "
            f"вообще не попали в final_detected. Filtered={u_in_filtered}, "
            f"final U-chain={u_in_final}")


def test_derive_chain_dominance_filename_arg_works():
    """Synthetic test: passing filename_chains_claimed={'Th-232'}
    suppresses U-238 verdict even with Ra-pair pattern."""
    from gamma.identification.anchor_ranks import (
        PatternConfirmation, ExpressPattern,
    )
    pat = ExpressPattern(
        name="Bi-214 Ra-chain pair",
        nuclide="Bi-214",
        required_lines_keV=(609.31, 1764.49),
        minimum_required=2,
        description="",
    )
    pc = PatternConfirmation(
        pattern=pat, matched_lines_keV=[609.31, 1764.49],
        missing_lines_keV=[], confirmed=True,
    )
    # Without filename binding: U-238 fires
    cd_default = derive_chain_dominance([], [pc])
    _assert(cd_default.u238 is True,
            "without filename binding, Ra-pair fires U-238 dominance")
    # With Th-only filename: U-238 suppressed
    cd_bound = derive_chain_dominance(
        [], [pc], filename_chains_claimed={"Th-232"},
    )
    _assert(cd_bound.u238 is False,
            f"filename binding must suppress U-238; got: {cd_bound.reason}")
    _assert("U-238" in cd_bound.suppressed_chains,
            "suppressed_chains must include U-238")
    _assert("filename" in (cd_bound.suppression_reason or "").lower(),
            "suppression_reason must mention filename binding")


def test_k40_fixture_no_suppression():
    """K-40 fixture: filename doesn't bind a chain, so no suppression."""
    r = _result_k()
    _assert(len(r.chain_dominance.suppressed_chains) == 0,
            f"K-40 fixture should not suppress chains; "
            f"got: {r.chain_dominance.suppressed_chains}")


# ════════════════════════════════════════════════════════════════════
# F-89a — background status surfacing
# ════════════════════════════════════════════════════════════════════

_VALID_BG_STATUS = {
    "absent_no_subtraction",
    "embedded_present_not_subtracted",
    "subtracted_from_external_file",
    # F-135 / v1.17.7 — добавилось 4-е значение при default=apply
    "auto_resolved_from_directory",
}


def test_background_status_field_populated_on_all_fixtures():
    """Все 3 demo фикстуры получают одно из валидных значений.
    F-135 (default apply): обычно `auto_resolved_from_directory`,
    если в той же папке найден подходящий фон."""
    for label, fn in [("Cs", _result_cs),
                      ("Th", _result_th),
                      ("K",  _result_k)]:
        r = fn()
        _assert(r.background_status in _VALID_BG_STATUS,
                f"{label} fixture background_status: {r.background_status!r}")


def test_background_status_in_json_header():
    r = _result_th()
    j = build_json_report(r)
    h = j["header"]
    _assert("background_status" in h,
            "JSON header missing background_status field")
    _assert(h["background_status"] in _VALID_BG_STATUS,
            f"unexpected status: {h['background_status']}")


def test_background_status_in_markdown_table():
    r = _result_th()
    md = build_markdown_report(r)
    # v1.17.4: Markdown is fully RU.  The status label is one of
    # "Фон ВЫЧТЕН (внешний файл)" / "Фон НЕ вычтен (...)".
    _assert("Фон" in md and ("НЕ вычтен" in md or "ВЫЧТЕН" in md),
            f"Markdown must show RU background status label")


def test_background_status_in_html_header():
    r = _result_th()
    html = build_html_report(r, plots=None)
    # F-114 / v1.17.3 — canonical interactive form: background status
    # is surfaced in the RU subtitle as "без вычитания фона" /
    # "фон вычтен" / "фон встроен".
    _assert(
        ("без вычитания фона" in html)
        or ("фон вычтен" in html)
        or ("фон встроен" in html)
        or ("NOT subtracted" in html),
        "HTML must show background status label"
    )


def test_chat_summary_inputs_line_present():
    """Chat summary must include 'Inputs: filename → X · bg ...' line."""
    r = _result_th()
    j = build_json_report(r)
    summ = build_chat_summary(None, json_dict=j)
    _assert("Inputs:" in summ and "filename" in summ.lower(),
            f"Chat summary must include filename input row; got:\n{summ}")
    _assert("bg" in summ.lower(),
            "Chat summary must mention background status")


# ════════════════════════════════════════════════════════════════════
# F-89b — JSON / Markdown / HTML carry filename hints
# ════════════════════════════════════════════════════════════════════

def test_json_header_has_filename_isotope_hints():
    r = _result_cs()
    j = build_json_report(r)
    h = j["header"]
    _assert("filename_isotope_hints" in h,
            "JSON header missing filename_isotope_hints field")
    _assert(h["filename_isotope_hints"] == ["Cs-137"],
            f"Cs-137 hint wrong: {h['filename_isotope_hints']}")


def test_json_diagnostics_has_chain_suppression():
    """JSON.diagnostics.chain_dominance carries suppressed_chains;
    chain_filtered_out_nuclides поле присутствует (может быть пустым,
    если v1.17.6 окно идентификации не пустило U-chain).
    """
    r = _result_th()
    j = build_json_report(r)
    cd = (j.get("diagnostics", {}) or {}).get("chain_dominance") or {}
    _assert("suppressed_chains" in cd,
            "diagnostics.chain_dominance missing suppressed_chains")
    _assert("U-238" in cd["suppressed_chains"],
            f"U-238 not in suppressed_chains: {cd['suppressed_chains']}")
    _assert("chain_filtered_out_nuclides" in cd,
            "diagnostics.chain_dominance missing chain_filtered_out_nuclides")
    # v1.17.6: chain_filtered_out_nuclides может быть пустым (предотвращено
    # на Stage). Главное — U-chain отсутствует в результатах.
    u_chain = {"Bi-214", "Pb-214", "Pb-210", "Ra-226"}
    detected_names = {n["nuclide"] for n in j.get("nuclides", []) or []}
    _assert(not (detected_names & u_chain),
            f"U-chain не должен попадать в результат; "
            f"detected={detected_names & u_chain}")


def test_markdown_shows_chain_suppression_box():
    # F-388 / v1.18.26 — conditional render: блок показывается ⟺ хотя бы
    # один primary_fep принадлежит подавляемой цепочке. На Th-232 fixture
    # U-chain полностью отсутствует в primary_feps (Stage-1 отбраковка
    # по filename binding), поэтому блок СКРЫТ. F-89d info-канал смещён
    # в JSON diagnostics (test_json_diagnostics_has_chain_suppression).
    r = _result_th()
    md = build_markdown_report(r)
    feps = build_json_report(r).get("primary_feps") or []
    u_chain = {"Pb-214", "Bi-214", "Pb-210", "Bi-210", "Po-214", "Po-218",
               "Ra-226", "Rn-222", "U-238", "U-234", "Th-234", "Pa-234"}
    has_u_pickup = any((p.get("nuclide") or "") in u_chain for p in feps)
    if has_u_pickup:
        _assert("Подавление цепочки" in md,
                "Markdown must show 'Подавление цепочки' block when U-pickup present")
        _assert("U-238" in md,
                "Markdown must mention U-238 in suppression block")
    else:
        _assert("Подавление цепочки" not in md,
                "F-388: Markdown must NOT show suppression block when no U-pickup in primary_feps")


def test_html_shows_chain_suppression_box():
    # F-388 / v1.18.26 — conditional render симметрично markdown.
    r = _result_th()
    html = build_html_report(r, plots=None)
    feps = build_json_report(r).get("primary_feps") or []
    u_chain = {"Pb-214", "Bi-214", "Pb-210", "Bi-210", "Po-214", "Po-218",
               "Ra-226", "Rn-222", "U-238", "U-234", "Th-234", "Pa-234"}
    has_u_pickup = any((p.get("nuclide") or "") in u_chain for p in feps)
    if has_u_pickup:
        _assert(
            ("не определяется" in html) or ("Chain suppression" in html),
            "HTML must surface chain-suppression narrative when U-pickup present"
        )
        _assert("U-238" in html,
            "HTML must mention U-238 in the suppression narrative")
    else:
        # F-388: блок отсутствует в HTML notes-block (но slovar остаётся в JSON diagnostics)
        _assert("Почему U-238" not in html,
                "F-388: HTML must NOT show '«Почему U-238 ...»' block when no U-pickup")


def test_schema_version_bumped_to_0_3():
    # Forward-compat: F-89 bumped schema to 0.3; later releases may bump higher.
    # We require schema >= 0.3 and skill_version >= v1.15.2 (lexical
    # comparison works on the "vMAJOR.MINOR.PATCH" prefix).
    r = _result_th()
    j = build_json_report(r)
    sv = j["schema_version"]
    parts = tuple(int(p) for p in sv.split("."))
    _assert(parts >= (0, 3),
            f"expected schema ≥ 0.3, got {sv}")
    _assert(j["skill_version"] >= "v1.15.2",
            f"expected skill ≥ v1.15.2, got {j['skill_version']}")


# ════════════════════════════════════════════════════════════════════
# Runner
# ════════════════════════════════════════════════════════════════════

ALL_TESTS = [
    # F-89b — filename → canonical
    ("F-89b   extract_isotope_hints canonical form (Cs137 → Cs-137 etc.)",
     test_extract_isotope_hints_canonical_form),
    ("F-89b   chains_claimed_by_isotope_hints chain mapping",
     test_chains_claimed_by_isotope_hints),
    ("F-89b   parse_filename carries isotope_hints + chains_claimed",
     test_parse_filename_carries_hints_and_chains),
    # F-89e — filename drives candidates (the main "полный провал" fix)
    ("F-89e   Cs-137 fixture → Cs-137 in identified list (NEW)",
     test_cs137_fixture_identifies_cs137),
    ("F-89e   K-40 fixture → K-40 in identified list",
     test_k40_fixture_identifies_k40),
    ("F-89e   StagedAnalysisResult carries filename_isotope_hints",
     test_orchestrator_carries_filename_isotope_hints),
    # F-89d — chain suppression
    ("F-89d   Th-232 fixture suppresses U-chain (no Bi-214/Pb-214/Pb-210)",
     test_th232_fixture_suppresses_u238_chain),
    ("F-89d   Th-232 fixture keeps Th-chain nuclides (Tl-208/Ac-228)",
     test_th232_fixture_keeps_th_chain_nuclides),
    ("F-89d   Th-232 fixture records suppression in result",
     test_th232_fixture_records_suppression),
    ("F-89d   derive_chain_dominance filename_chains_claimed kwarg works",
     test_derive_chain_dominance_filename_arg_works),
    ("F-89d   K-40 fixture (no chain claim) → no suppression",
     test_k40_fixture_no_suppression),
    # F-89a — bg status
    ("F-89a   background_status field populated on all fixtures",
     test_background_status_field_populated_on_all_fixtures),
    ("F-89a   JSON header has background_status",
     test_background_status_in_json_header),
    ("F-89a   Markdown header table shows bg status label",
     test_background_status_in_markdown_table),
    ("F-89a   HTML header shows bg status label",
     test_background_status_in_html_header),
    ("F-89a   Chat summary has 'Inputs:' row with filename + bg",
     test_chat_summary_inputs_line_present),
    # F-89b/d surfacing
    ("F-89b   JSON header has filename_isotope_hints field",
     test_json_header_has_filename_isotope_hints),
    ("F-89d   JSON diagnostics has chain suppression block",
     test_json_diagnostics_has_chain_suppression),
    ("F-89d   Markdown shows chain-suppression notice",
     test_markdown_shows_chain_suppression_box),
    ("F-89d   HTML shows chain-suppression block",
     test_html_shows_chain_suppression_box),
    ("F-89    Schema bumped 0.2 → 0.3; skill v1.15.2",
     test_schema_version_bumped_to_0_3),
]


def main() -> int:
    print("=" * 72)
    print("test_filename_binding.py — F-89 / v1.15.2")
    print("=" * 72)
    for f in (FIX_TH, FIX_CS, FIX_K):
        if not f.exists():
            print(f"FATAL: fixture missing: {f}", file=sys.stderr)
            return 1
    passed = failed = 0
    for name, fn in ALL_TESTS:
        if _report(name, fn):
            passed += 1
        else:
            failed += 1
    print()
    print("=" * 72)
    print(f"F-89 filename binding: {passed}/{passed + failed} pass, "
          f"{failed} fail")
    print("=" * 72)
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
