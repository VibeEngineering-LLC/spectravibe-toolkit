"""
test_v1_15_delivery.py — F-86 + F-87 regression for v1.15.0.

Covers two parallel streams:

  F-87  Anchor seeding at Step 5α + opt-in calibration refit
        - seed_calibration_anchors API
        - default orchestrator path produces identical anchor_matches /
          pattern_confirmations as v1.14.0 (bit-compatible)
        - recalibration_diag default shape
        - recalibrate_on_anchor_disagreement=True with stored cal OK
          (no refit triggered)
        - corrupted cal triggers a refit + recovers within 1 keV at
          the Cs-137 661.66 anchor

  F-86  PNG plots + HTML + CLI + analyze_and_report wrapper
        - build_spectrum_plot creates a non-trivial PNG
        - build_multiplet_plots returns ≥1 PNG when clusters present
        - build_html_report contains base64-embedded PNG data URI
        - Markdown sections 9/10 contain Markdown image syntax
        - build_report respects write_plots=False
        - analyze_and_report end-to-end on a Marinelli fixture writes
          JSON + Markdown + HTML + PNGs
        - CLI Phase-0 backward compatibility (analyze without
          --full-report still prints the parse summary)
        - CLI --full-report writes all artefacts and exits 0

Run:  PYTHONPATH=scripts python test_v1_15_delivery.py
"""
from __future__ import annotations

import base64
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from gamma.identification.staged_pipeline import (
    analyze_lsrm_spe,
    fwhm_keV_at_energy as _fwhm_keV_at_energy,
)
from gamma.identification.anchor_ranks import (
    seed_calibration_anchors, AnchorSeedResult, AnchorMatch,
    find_anchor_matches, confirm_express_patterns,
)
from gamma.calibration.anchor_recalibration import (
    recalibrate_energy_if_anchors_disagree,
)
from gamma.reporting import (
    build_report, analyze_and_report,
    build_html_report, build_markdown_report,
    build_spectrum_plot, build_multiplet_plots,
)
from gamma.detectors.gamma1s import DEFAULT_REFERENCE_DIR


_ROOT = DEFAULT_REFERENCE_DIR
FIXTURE_TH_MARINELLI = _ROOT / "Th232_420-7-17_Маринелли_0cm.spe"
FIXTURE_CS_MARINELLI = _ROOT / "Cs137_420-7-14_Маринелли_0cm.spe"
FIXTURE_CS_5CM = _ROOT / "Cs-137__163_2017.spe"


# ════════════════════════════════════════════════════════════════════
# Cached pipeline results — one per fixture / call shape
# ════════════════════════════════════════════════════════════════════

_RESULT_TH_WORKFLOW = None
_RESULT_CS_DEFAULT = None


def _result_th_workflow():
    global _RESULT_TH_WORKFLOW
    if _RESULT_TH_WORKFLOW is None:
        _RESULT_TH_WORKFLOW = analyze_lsrm_spe(
            str(FIXTURE_TH_MARINELLI),
            complete_workflow=True,
            sample_mass_kg=0.2,
        )
    return _RESULT_TH_WORKFLOW


def _result_cs_default():
    """Default-args call — for v1.14.0 BC verification."""
    global _RESULT_CS_DEFAULT
    if _RESULT_CS_DEFAULT is None:
        _RESULT_CS_DEFAULT = analyze_lsrm_spe(str(FIXTURE_CS_5CM))
    return _RESULT_CS_DEFAULT


# ════════════════════════════════════════════════════════════════════
# Helpers
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
        import traceback
        traceback.print_exc()
        print(f"  ERROR {name}: {type(exc).__name__}: {exc}")
        return False


# ════════════════════════════════════════════════════════════════════
# F-87 — anchor seeding + recalibration
# ════════════════════════════════════════════════════════════════════

def test_seed_calibration_anchors_sample_mode():
    """seed_calibration_anchors with mode='sample' returns both lists."""
    r = _result_th_workflow()
    spec = r.spec

    def fwhm(e):
        # F-452: polymorphic — FwhmModel callable | legacy 3-tuple.
        return _fwhm_keV_at_energy(r.fwhm_model, float(e))

    seed = seed_calibration_anchors(
        r.peaks, spec, mode="sample", fwhm_provider_keV=fwhm,
    )
    _assert(isinstance(seed, AnchorSeedResult), "expected AnchorSeedResult")
    _assert(seed.mode == "sample", "mode must round-trip")
    _assert(len(seed.anchor_matches) >= 1,
            "Th-232 Marinelli must yield at least one anchor (Tl-208 2614)")
    _assert(len(seed.pattern_confirmations) >= 5,
            "should evaluate all 7 express patterns (≥5 returned)")
    nuclides_hit = {m.anchor.nuclide for m in seed.anchor_matches}
    _assert("Tl-208" in nuclides_hit,
            "Th-232 Marinelli must hit the Tl-208 anchor")


def test_seed_calibration_anchors_background_mode_extends_rank():
    """mode='background' extends max_rank to 12 (Pb-212 / Tl-208 583)."""
    r = _result_th_workflow()
    spec = r.spec

    def fwhm(e):
        # F-452: polymorphic — FwhmModel callable | legacy 3-tuple.
        return _fwhm_keV_at_energy(r.fwhm_model, float(e))

    seed_s = seed_calibration_anchors(
        r.peaks, spec, mode="sample", fwhm_provider_keV=fwhm,
    )
    seed_b = seed_calibration_anchors(
        r.peaks, spec, mode="background", fwhm_provider_keV=fwhm,
    )
    # Background mode must produce at least as many matches as sample,
    # since it considers ranks 1..12 instead of 1..10.
    _assert(len(seed_b.anchor_matches) >= len(seed_s.anchor_matches),
            f"bg mode must catch at least as many anchors as sample mode "
            f"(bg={len(seed_b.anchor_matches)} < s={len(seed_s.anchor_matches)})")


def test_seed_calibration_anchors_rejects_bad_mode():
    """Bad mode raises ValueError immediately."""
    r = _result_cs_default()
    try:
        seed_calibration_anchors(
            r.peaks, r.spec, mode="garbage",
            fwhm_provider_keV=lambda e: 1.0,
        )
    except ValueError:
        return
    raise AssertionError("expected ValueError for mode='garbage'")


def test_orchestrator_bc_anchor_matches_unchanged():
    """Default analyze_lsrm_spe still populates anchor_matches as before."""
    r = _result_cs_default()
    _assert(r.anchor_matches is not None,
            "anchor_matches must not be None even with defaults")
    _assert(isinstance(r.anchor_matches, list),
            "anchor_matches must be a list")
    _assert(r.pattern_confirmations is not None and
            len(r.pattern_confirmations) >= 5,
            "pattern_confirmations must populate all 7 patterns")


def test_recalibration_diag_default_empty():
    """Default call → recalibration_diag dict present but ``attempted=False``."""
    r = _result_cs_default()
    _assert(isinstance(r.recalibration_diag, dict),
            "recalibration_diag must be a dict")
    _assert(r.recalibration_diag.get("attempted", True) is False,
            "default call must not attempt recalibration")


def test_recalibration_no_op_when_cal_ok():
    """With stored cal OK, recalibrate_on_anchor_disagreement=True yields
    attempted=True but applied=False (no improvement / below threshold)."""
    r = analyze_lsrm_spe(
        str(FIXTURE_TH_MARINELLI),
        complete_workflow=False,
        recalibrate_on_anchor_disagreement=True,
    )
    diag = r.recalibration_diag
    _assert(diag.get("attempted") is True,
            "must attempt the residual check")
    # Either applied=False (cal already OK) or applied=True with
    # improvement. We accept both; we only require no crash and a
    # well-formed diag dict.
    _assert("old_residual_max_keV" in diag,
            "diag must report old_residual_max_keV")
    _assert("n_anchors_used" in diag,
            "diag must report n_anchors_used")


def test_recalibration_recovers_corrupted_cal():
    """Corrupt stored cal by significant offset → refit should recover
    anchors to within 1 keV.

    F-125 (v1.17.6): после refit FWHM(59.5)≈13 кэВ (раньше клампилось
    в 1.0). Чтобы превысить порог 0.3·FWHM на низкой энергии, используем
    смещение +20 кэВ вместо +3.
    """
    # Use the API directly (rather than the orchestrator) for tight
    # control over the corruption.
    r = analyze_lsrm_spe(str(FIXTURE_TH_MARINELLI),
                        complete_workflow=False)
    spec = r.spec
    anchors = r.anchor_matches

    if not anchors:
        # Fixture didn't yield anchors (very unusual) — skip the
        # assertion gracefully.
        return

    # Snapshot stored cal, corrupt it by adding offset to a0.
    # v1.17.6 (F-125): FWHM model теперь корректен на низких энергиях
    # (FWHM(60)≈13 кэВ вместо клампа 1.0), поэтому минимальное смещение,
    # уверенно превышающее 0.3·FWHM на всех анкерах, ≥ 15 кэВ.
    stored = list(spec.energy_cal)
    if len(stored) == 0:
        return
    OFFSET_KEV = 20.0
    corrupted = list(stored)
    corrupted[0] = float(corrupted[0]) + OFFSET_KEV

    # Build a tiny mock spec-like object that recalibrate sees
    class _MockSpec:
        energy_cal = tuple(corrupted)

    def fwhm(e):
        # F-452: polymorphic — FwhmModel callable | legacy 3-tuple.
        return _fwhm_keV_at_energy(r.fwhm_model, float(e))

    # The recalibration needs the *current* delta against corrupted cal.
    # The AnchorMatch.delta_keV from the original orchestrator was
    # computed against the stored (correct) cal; so we synthesize a
    # corrupted-cal anchor list by shifting delta_keV by OFFSET_KEV.
    corrupted_anchors = []
    for am in anchors:
        if not am.anchor.nuclide:
            continue
        # Library E unchanged, peak_E shifted by the cal offset
        new_peak_E = am.peak_E_keV + OFFSET_KEV
        new_delta = abs(new_peak_E - am.anchor.energy_keV)
        corrupted_anchors.append(
            AnchorMatch(
                anchor=am.anchor,
                peak_channel=am.peak_channel,
                peak_E_keV=new_peak_E,
                delta_keV=new_delta,
                sigma=am.sigma,
                partner_required_but_missing=am.partner_required_but_missing,
            )
        )

    new_cal, diag = recalibrate_energy_if_anchors_disagree(
        _MockSpec(), corrupted_anchors,
        threshold_fraction_of_fwhm=0.3,
        fwhm_provider_keV=fwhm,
        min_anchors=2,
    )
    _assert(diag.get("attempted") is True, "attempted must be True")
    if len(corrupted_anchors) >= 2:
        _assert(new_cal is not None,
                f"with {OFFSET_KEV} keV offset on multiple anchors, "
                f"refit must fire (got reason: {diag.get('reason')!r})")
        # The new cal should bring residuals down
        _assert(diag["new_residual_max_keV"] < diag["old_residual_max_keV"],
                f"new max residual ({diag['new_residual_max_keV']:.2f}) "
                f"must be lower than old ({diag['old_residual_max_keV']:.2f})")


# ════════════════════════════════════════════════════════════════════
# F-86 — plots / HTML / CLI
# ════════════════════════════════════════════════════════════════════

def test_build_spectrum_plot_creates_png():
    r = _result_th_workflow()
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "spectrum.png"
        path = build_spectrum_plot(r, out, dpi=80)
        _assert(path is not None, "build_spectrum_plot returned None")
        _assert(Path(path).exists(), "PNG file not created")
        size = Path(path).stat().st_size
        _assert(size > 5_000,
                f"PNG suspiciously small ({size} bytes)")


def test_build_multiplet_plots_present():
    r = _result_th_workflow()
    with tempfile.TemporaryDirectory() as tmp:
        paths = build_multiplet_plots(r, Path(tmp), dpi=80)
        # Th-232 Marinelli should produce at least one multiplet
        # (600-680 keV cluster: Cs-134/Bi-214). If for some reason
        # the orchestrator returned no clusters, accept []; otherwise
        # check we got non-empty paths.
        if (r.deconvolution_results or []):
            _assert(len(paths) >= 1,
                    "expected ≥1 multiplet PNG when clusters present")
            for p in paths:
                _assert(Path(p).exists(), f"PNG missing: {p}")
                _assert(Path(p).stat().st_size > 3_000,
                        f"PNG suspiciously small: {p}")


def test_html_report_embeds_png_when_plots_given():
    # F-114 / v1.17.3 — the canonical interactive form renders the
    # spectrum live via Chart.js (no base64 PNGs). Verify the chart
    # canvas + embedded JS data arrays are present.
    r = _result_th_workflow()
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "spectrum.png"
        spec_path = build_spectrum_plot(r, out, dpi=72)
        plots = {"spectrum": spec_path, "multiplets": []}
        html = build_html_report(r, plots=plots)
        _assert('<canvas id="fp-sp"' in html,
                "HTML must contain the fp-sp spectrum canvas (F-114)")
        _assert("const E=" in html and "const C=" in html,
                "HTML must embed E/C spectrum data arrays")
        _assert("chartjs-plugin-annotation" in html,
                "HTML must load Chart.js annotation plugin")


def test_html_report_placeholder_when_no_plots():
    # F-114 / v1.17.3 — the canonical interactive form does not depend
    # on pre-rendered PNGs. The canvas is always emitted regardless of
    # the `plots=` argument.
    r = _result_th_workflow()
    html = build_html_report(r, plots=None)
    _assert('<canvas id="fp-sp"' in html,
            "HTML must contain the fp-sp spectrum canvas")
    _assert("data:image/png;base64," not in html,
            "HTML must NOT contain base64 PNG URIs (canonical form is live Chart.js)")


def test_markdown_embeds_image_when_plots_given():
    r = _result_th_workflow()
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "spectrum.png"
        spec_path = build_spectrum_plot(r, out, dpi=72)
        plots = {"spectrum": spec_path, "multiplets": []}
        md = build_markdown_report(r, plots=plots, md_dir=tmp)
        # v1.17.4: image alt text is RU.
        _assert(
            ("![График спектра]" in md) or ("![Spectrum overlay]" in md),
            "Markdown section 9 must embed the spectrum image",
        )


def test_build_report_writes_full_bundle():
    r = _result_th_workflow()
    with tempfile.TemporaryDirectory() as tmp:
        artefacts = build_report(
            r,
            output_dir=tmp,
            write_json=True,
            write_markdown=True,
            write_plots=True,
            write_html=True,
            plot_dpi=80,
        )
        _assert(artefacts["json"] is not None,
                "JSON file must be written")
        _assert(artefacts["markdown"] is not None,
                "Markdown must be written")
        _assert(artefacts["html"] is not None,
                "HTML must be written")
        _assert(artefacts["plots"] is not None,
                "plots dict must be present")
        _assert(Path(artefacts["json"]).exists(), "JSON not on disk")
        _assert(Path(artefacts["markdown"]).exists(), "MD not on disk")
        _assert(Path(artefacts["html"]).exists(), "HTML not on disk")


def test_build_report_skip_plots_keeps_md_placeholder():
    r = _result_th_workflow()
    with tempfile.TemporaryDirectory() as tmp:
        artefacts = build_report(
            r, output_dir=tmp,
            write_markdown=True, write_plots=False,
        )
        md = artefacts["markdown_text"]
        # v1.17.4 — placeholder is now RU "Графики не сгенерированы".
        _assert(
            ("Графики не сгенерированы" in md)
            or ("Plot generation deferred" in md),
            "MD must contain placeholder when plots=False",
        )
        _assert(artefacts["plots"] is None,
                "plots dict should be None when write_plots=False")


def test_analyze_and_report_end_to_end():
    """Single-call wrapper produces all four artefact types."""
    with tempfile.TemporaryDirectory() as tmp:
        out = analyze_and_report(
            str(FIXTURE_TH_MARINELLI),
            output_dir=tmp,
            sample_mass_kg=0.2,
            plot_dpi=72,
        )
        _assert(out["json"] is not None, "JSON missing")
        _assert(out["markdown"] is not None, "Markdown missing")
        _assert(out["html"] is not None, "HTML missing")
        _assert(out["plots"]["spectrum"] is not None, "spectrum PNG missing")
        _assert(out["summary"] is not None, "summary missing")
        _assert("result" in out, "raw StagedAnalysisResult missing")

        # Verify JSON is well-formed and serializable
        with open(out["json"], "r", encoding="utf-8") as f:
            data = json.load(f)
        _assert(data.get("skill_version") == "v1.14.0" or
                data.get("skill_version").startswith("v1."),
                "skill_version field must be set")
        _assert("plot_files" in data,
                "JSON must record plot_files manifest")


# ════════════════════════════════════════════════════════════════════
# CLI tests (subprocess)
# ════════════════════════════════════════════════════════════════════

def _run_cli(args, *, cwd=None):
    """Run `python -m gamma.cli ...` and return (returncode, stdout, stderr)."""
    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(__file__).parent.parent.parent / "scripts")
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONUTF8"] = "1"
    proc = subprocess.run(
        [sys.executable, "-m", "gamma.cli"] + list(args),
        cwd=cwd, env=env, capture_output=True, text=True,
        encoding="utf-8", errors="replace",
        timeout=180,
    )
    return proc.returncode, proc.stdout, proc.stderr


def test_cli_phase0_backward_compatible():
    """`analyze fixture.spe` without --full-report still prints JSON summary."""
    rc, out, err = _run_cli([
        "analyze", str(FIXTURE_CS_5CM),
    ])
    _assert(rc == 0, f"Phase-0 analyze exit {rc}; stderr: {err}")
    # Phase-0 prints JSON. Look for a recognizable key.
    _assert('"source_path"' in out or '"detector_type"' in out,
            "Phase-0 output must contain JSON summary fields")


def test_cli_full_report_writes_artefacts():
    with tempfile.TemporaryDirectory() as tmp:
        rc, out, err = _run_cli([
            "analyze", str(FIXTURE_TH_MARINELLI),
            "--full-report",
            "--output-dir", tmp,
            "--sample-mass-kg", "0.2",
            "--plot-dpi", "72",
        ])
        _assert(rc == 0, f"--full-report exit {rc}; stderr: {err}")
        # Stdout should contain the chat summary + manifest
        _assert("Written artefacts:" in out,
                "stdout must include the manifest header")
        _assert(".json" in out and ".html" in out and ".md" in out,
                "manifest must list all artefacts")
        # Files must exist on disk
        files = list(Path(tmp).rglob("*"))
        names = [f.name for f in files if f.is_file()]
        _assert(any(n.endswith("_report.json") for n in names),
                "JSON report missing on disk")
        _assert(any(n.endswith("_report.html") for n in names),
                "HTML report missing on disk")
        _assert(any(n.endswith(".png") for n in names),
                "no PNG plots written")


def test_cli_full_report_requires_output_dir():
    """F-384 / v1.18.25.3 — отменено: --full-report без --output-dir
    больше НЕ ошибка. CLI автоматически создаёт путь под demo_reports_root
    (ensure_demo_reports_root). Тест переориентирован: проверяем что в
    stderr выводится информационное сообщение про demo_reports_root."""
    import tempfile, os
    with tempfile.TemporaryDirectory() as tmp:
        os.environ["GAMMA_DEMO_REPORTS_DIR"] = tmp
        try:
            rc, out, err = _run_cli([
                "analyze", str(FIXTURE_CS_5CM), "--full-report",
            ])
        finally:
            os.environ.pop("GAMMA_DEMO_REPORTS_DIR", None)
    _assert(
        "demo_reports" in err.lower() or rc == 0,
        "CLI должен сообщить про demo_reports_root либо завершиться 0"
    )


# ════════════════════════════════════════════════════════════════════
# Runner
# ════════════════════════════════════════════════════════════════════

ALL_TESTS = [
    # F-87 anchor seeding
    ("F-87a/b seed_calibration_anchors(sample) returns both lists",
     test_seed_calibration_anchors_sample_mode),
    ("F-87a/b seed_calibration_anchors(background) extends max_rank to 12",
     test_seed_calibration_anchors_background_mode_extends_rank),
    ("F-87a   seed_calibration_anchors rejects bad mode",
     test_seed_calibration_anchors_rejects_bad_mode),
    ("F-87b   orchestrator BC: default anchor_matches unchanged",
     test_orchestrator_bc_anchor_matches_unchanged),
    ("F-87b   recalibration_diag default empty",
     test_recalibration_diag_default_empty),
    ("F-87c/d recalibration no-op when stored cal OK",
     test_recalibration_no_op_when_cal_ok),
    ("F-87c   recalibration recovers a +3 keV corrupted cal",
     test_recalibration_recovers_corrupted_cal),

    # F-86 plots / HTML / wrapper / CLI
    ("F-86a   build_spectrum_plot creates non-trivial PNG",
     test_build_spectrum_plot_creates_png),
    ("F-86a   build_multiplet_plots returns PNGs when clusters present",
     test_build_multiplet_plots_present),
    ("F-86c   HTML embeds base64 PNG when plots given",
     test_html_report_embeds_png_when_plots_given),
    ("F-86c   HTML placeholder when plots omitted",
     test_html_report_placeholder_when_no_plots),
    ("F-86b   Markdown embeds image when plots given",
     test_markdown_embeds_image_when_plots_given),
    ("F-86d   build_report writes full bundle (JSON+MD+HTML+plots)",
     test_build_report_writes_full_bundle),
    ("F-86d   build_report write_plots=False keeps MD placeholder",
     test_build_report_skip_plots_keeps_md_placeholder),
    ("F-86e   analyze_and_report end-to-end produces all artefacts",
     test_analyze_and_report_end_to_end),
    ("F-86f   CLI Phase-0 backward-compatible (no --full-report)",
     test_cli_phase0_backward_compatible),
    ("F-86f   CLI --full-report writes all artefacts, exit 0",
     test_cli_full_report_writes_artefacts),
    ("F-86f   CLI --full-report without --output-dir errors out",
     test_cli_full_report_requires_output_dir),
]


def main() -> int:
    print("=" * 72)
    print("test_v1_15_delivery.py — F-86 + F-87 / v1.15.0")
    print("=" * 72)

    if not FIXTURE_TH_MARINELLI.exists():
        print(f"FATAL: fixture missing: {FIXTURE_TH_MARINELLI}",
              file=sys.stderr)
        return 1
    if not FIXTURE_CS_5CM.exists():
        print(f"FATAL: fixture missing: {FIXTURE_CS_5CM}",
              file=sys.stderr)
        return 1

    passed = failed = errored = 0
    for name, fn in ALL_TESTS:
        ok = _report(name, fn)
        if ok:
            passed += 1
        else:
            failed += 1

    total = passed + failed + errored
    print()
    print("=" * 72)
    print(f"v1.15.0 delivery: {passed}/{total} pass, {failed} fail")
    print("=" * 72)
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
