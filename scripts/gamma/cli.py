"""
CLI entry point.

Commands:

    python -m gamma.cli analyze <spectrum_file>
        Phase-0 default: parse the file, also try to auto-resolve the
        external background if a BackgroundSpectrumFile link is set,
        and print a JSON summary to stdout.

    python -m gamma.cli analyze <spectrum_file> --full-report
        F-86 (v1.15.0): run the full Step-1..11 pipeline (Round 5 on
        by default — multiplet deconvolution + activities + ISO 11929
        MDA) and write the Step-11 report bundle to ``--output-dir``.
        Prints the 3–8 line chat summary to stdout.

    python -m gamma.cli summarize <spectrum_file>
        Alias for the Phase-0 analyze command.

The Phase-0 path is preserved so existing callers continue to work.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from gamma.io.readers import read_spectrum
from gamma.io.background import resolve_external_background
from gamma.data.demo_reports_root import (
    ensure_demo_reports_root,
    get_demo_reports_root_default,
)


# ──────────────────────────────────────────────────────────────────
# Phase-0 path (parse only)
# ──────────────────────────────────────────────────────────────────

def _cmd_analyze_phase0(args) -> int:
    path = Path(args.spectrum)
    if not path.exists():
        print(f"ERROR: file not found: {path}", file=sys.stderr)
        return 2

    try:
        spec = read_spectrum(str(path))
    except (ValueError, NotImplementedError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 3

    # Try to resolve the external background file if its hint is set
    # and no embedded background is already present
    if spec.background_link and spec.background_embedded is None:
        external_bg = resolve_external_background(spec, args.search_dir)
        if external_bg is not None:
            spec.background_embedded = external_bg
            spec.extras["background_source"] = "external_resolved"
        else:
            spec.extras["background_source"] = "link_unresolved"
    elif spec.background_embedded is not None:
        spec.extras["background_source"] = "embedded"
    else:
        spec.extras["background_source"] = "none"

    summary = spec.to_summary_dict()
    summary["extras"] = spec.extras

    print(json.dumps(summary, indent=2, ensure_ascii=False, default=str))
    return 0


# ──────────────────────────────────────────────────────────────────
# F-86 full report path (v1.15.0)
# ──────────────────────────────────────────────────────────────────

def _cmd_analyze_full_report(args) -> int:
    path = Path(args.spectrum)
    if not path.exists():
        print(f"ERROR: file not found: {path}", file=sys.stderr)
        return 2

    # F-384 / v1.18.25.3 — demo_reports first-run check.
    # СРАБАТЫВАЕТ ТОЛЬКО когда --output-dir НЕ задан явно:
    # создаётся autopath <demo_root>/<stem>/ через
    # ensure_demo_reports_root() (interactive prompt при первом запуске
    # для выбора demo_root, поддержка env var GAMMA_DEMO_REPORTS_DIR).
    # Если пользователь задал --output-dir явно (любой путь) — используем
    # его как есть, никаких side-effects.
    if not args.output_dir:
        demo_root = ensure_demo_reports_root()
        out_dir = demo_root / path.stem
        args.output_dir = str(out_dir)
        print(
            f"[gamma] --output-dir не задан, использую "
            f"{args.output_dir} (под demo_reports_root)",
            file=sys.stderr,
        )

    # Lazy import — keeps Phase-0 startup fast
    try:
        from gamma.reporting import analyze_and_report
    except ImportError as e:
        print(f"ERROR: gamma.reporting unavailable: {e}", file=sys.stderr)
        return 5

    kwargs = {}
    if args.sample_mass_kg is not None:
        kwargs["sample_mass_kg"] = float(args.sample_mass_kg)
    if args.background_path:
        kwargs["background_path"] = args.background_path
    if args.allow_stage2:
        kwargs["allow_stage2"] = True
    if args.allow_stage3:
        kwargs["allow_stage3"] = True
    if args.recalibrate_on_anchor_disagreement:
        kwargs["recalibrate_on_anchor_disagreement"] = True
    # F-122 / v1.17.6 — self-attenuation для Marinelli / Дента / Петри.
    if getattr(args, "sample_density_g_cm3", None) is not None:
        kwargs["sample_density_g_cm3"] = float(args.sample_density_g_cm3)
    # F-129 / v1.17.7 — выбор метода поиска пиков.
    if getattr(args, "peak_search_method", None):
        kwargs["peak_search_method"] = str(args.peak_search_method)
    # F-139 / v1.17.7 — отбраковка узких пиков.
    if getattr(args, "filter_narrow_peaks", False):
        kwargs["filter_narrow_peaks"] = True
    if getattr(args, "narrow_peak_fwhm_ratio", None) is not None:
        kwargs["narrow_peak_fwhm_ratio"] = float(args.narrow_peak_fwhm_ratio)
    # F-131 / v1.17.7 — авто-поиск фона.
    if getattr(args, "background_auto", None):
        kwargs["background_auto"] = str(args.background_auto)
    if getattr(args, "background_auto_max_days", None) is not None:
        kwargs["background_auto_max_days"] = int(args.background_auto_max_days)

    # F-309 / v1.18.8 — opt-in флаги activity integrations v1.18.1..v1.18.4.
    # Прокидываются через analyze_lsrm_spe (F-308) в compute_activities_for_all.
    if getattr(args, "enable_tcs_correction", False):
        kwargs["enable_tcs_correction"] = True
    if getattr(args, "tcs_detector_id", None):
        kwargs["tcs_detector_id"] = str(args.tcs_detector_id)
    if getattr(args, "enable_cutshall_self_abs", False):
        kwargs["enable_cutshall_self_abs"] = True
    if getattr(args, "cutshall_path_cm", None) is not None:
        kwargs["cutshall_path_cm"] = float(args.cutshall_path_cm)
    if getattr(args, "cutshall_calib_density_g_cm3", None) is not None:
        kwargs["cutshall_calib_density_g_cm3"] = float(
            args.cutshall_calib_density_g_cm3
        )
    if getattr(args, "enable_matrix_method", False):
        kwargs["enable_matrix_method"] = True
    if getattr(args, "matrix_method_energy_tolerance_keV", None) is not None:
        kwargs["matrix_method_energy_tolerance_keV"] = float(
            args.matrix_method_energy_tolerance_keV
        )

    # F-108 / D-19 cost-estimate footer (tokens, session %, detail).
    # F-132 / v1.17.7 — cost footer теперь ОБЯЗАТЕЛЬНЫЙ; авто-оценка
    # включается всегда. CLI флаги работают как override:
    #   --cost-tokens N             → переопределить итог
    #   --cost-session-pct "65%"    → не используется (% теперь авто)
    #   --cost-detail "v1.17.7"     → строка-комментарий
    #   --cost-session-token-budget → бюджет 5-час. сессии (default 200k)
    cost_estimate = {}
    if args.cost_tokens is not None:
        cost_estimate["tokens"] = int(args.cost_tokens)
    if args.cost_detail:
        cost_estimate["detail"] = str(args.cost_detail)
    if getattr(args, "cost_session_token_budget", None):
        cost_estimate["session_token_budget"] = int(
            args.cost_session_token_budget
        )
    # Передаём даже пустой dict — build_report сам подставит авто-оценку.
    kwargs["cost_estimate"] = cost_estimate or {}
    if getattr(args, "write_pdf", False):
        kwargs["write_pdf"] = True
    # F-160 / v1.18.19.0 — BecqMoni XML export опция
    if getattr(args, "export_becqmoni", "off") != "off":
        kwargs["export_becqmoni"] = str(args.export_becqmoni)

    try:
        artefacts = analyze_and_report(
            str(path),
            output_dir=args.output_dir,
            write_json=not args.no_json,
            write_markdown=not args.no_markdown,
            write_plots=not args.no_plots,
            write_html=not args.no_html,
            # F-RPT-03 / v1.18.29 — Technical PDF default OFF; вкл через
            # --technical-pdf (новый opt-in). --no-technical-pdf оставлен
            # как back-compat shim (always wins as off если задан).
            write_technical_pdf=(
                bool(getattr(args, "technical_pdf", False))
                and not bool(getattr(args, "no_technical_pdf", False))
            ),
            plot_dpi=int(args.plot_dpi),
            **kwargs,
        )
    except (ValueError, NotImplementedError, FileNotFoundError) as e:
        print(f"ERROR: {type(e).__name__}: {e}", file=sys.stderr)
        return 6

    # Stdout: chat summary, then a single-line manifest of written paths
    summary = artefacts.get("summary") or "(no summary generated)"
    print(summary)
    print()
    print("Written artefacts:")
    for k in ("json", "markdown", "html"):
        v = artefacts.get(k)
        if v:
            print(f"  {k:10s} {v}")
    plots = artefacts.get("plots") or {}
    if plots.get("spectrum"):
        print(f"  spectrum   {plots['spectrum']}")
    for i, mp in enumerate(plots.get("multiplets") or [], start=1):
        print(f"  cluster {i:<2d} {mp}")
    # F-160 / v1.18.19.0 — surface XML round-trip outputs (sample + bg)
    # F-159 / v1.18.21.0 — surface technical PDF artefact
    # F-375: user-facing labels — нейтральные "xml_sample"/"xml_bg" вместо
    # бренда. Ключи в artefacts dict остаются исходными (контракт wrapper).
    _label_map = {
        "becqmoni_sample": "xml_sample",
        "becqmoni_bg":     "xml_bg",
    }
    for k in ("becqmoni_sample", "becqmoni_bg", "pdf", "technical_pdf"):
        v = artefacts.get(k)
        if v:
            label = _label_map.get(k, k)
            print(f"  {label:16s} {v}")

    for w in artefacts.get("warnings", []) or []:
        # BUG-40: warnings may be heterogeneous (str | dict).
        if isinstance(w, dict):
            w = w.get("message") or w.get("code") or str(w)
        print(f"WARNING: {w}", file=sys.stderr)

    return 0


# ──────────────────────────────────────────────────────────────────
# Dispatcher
# ──────────────────────────────────────────────────────────────────

def _cmd_analyze(args) -> int:
    if getattr(args, "full_report", False):
        return _cmd_analyze_full_report(args)
    return _cmd_analyze_phase0(args)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="gamma",
        description=(
            "Gamma-ray spectrum analysis. Default `analyze` is the "
            "Phase-0 file-parse summary; pass --full-report for the "
            "Step-1..11 pipeline with the Step-11 report bundle "
            "(JSON + Markdown + HTML + PNGs)."
        ),
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    # ── analyze (with --full-report opt-in) ──────────────────────
    p_analyze = sub.add_parser(
        "analyze",
        help="Parse a spectrum file (Phase-0) or run the full "
             "Step-1..11 pipeline with --full-report.",
    )
    p_analyze.add_argument(
        "spectrum",
        help="Path to the spectrum file (.spe LSRM or .xml AtomSpectra).",
    )
    p_analyze.add_argument(
        "--search-dir",
        action="append",
        default=[],
        help="Extra directory to search for the linked background file "
             "(can be given multiple times).",
    )
    # Full-report opt-in
    p_analyze.add_argument(
        "--full-report",
        action="store_true",
        help="Run the full Step-1..11 pipeline (Round 5 on by default) "
             "and write the Step-11 report bundle to --output-dir.",
    )
    p_analyze.add_argument(
        "--output-dir",
        default=None,
        help="Output directory for the report bundle "
             "(required with --full-report).",
    )
    p_analyze.add_argument(
        "--sample-mass-kg",
        type=float,
        default=None,
        help="Sample mass in kg — enables Bq/kg specific-activity "
             "derivation when an efficiency curve is available.",
    )
    p_analyze.add_argument(
        "--sample-density-g-cm3",
        type=float,
        default=None,
        help="(F-122 / v1.17.6) Плотность образца в г/см³ для "
             "коррекции самопоглощения в Marinelli / Дента / Петри "
             "геометриях. По умолчанию ОИСН-16 (ρ=1.60 г/см³). "
             "v1.17.7: если не задано, авто-определяется из "
             "MATERIAL.Ro / SAMPLEMASS+SAMPLEVOLUME заголовка .spe (F-130).",
    )
    p_analyze.add_argument(
        "--peak-search-method",
        choices=["mariscotti", "convolution", "compare"],
        default="mariscotti",
        help="(F-129 / v1.17.7) Метод поиска пиков: 'mariscotti' "
             "(default, second-derivative), 'convolution' (matched-filter "
             "Гауссиан, F-124), 'compare' (оба + сравнение в diagnostics).",
    )
    p_analyze.add_argument(
        "--filter-narrow-peaks", action="store_true",
        help="(F-139 / v1.17.7) Отбраковка узких шумовых пиков: измерять "
             "фактическую FWHM пика по полувысоте и отсекать кандидатов "
             "с FWHM < ratio·FWHM(E). Полезно для NaI 63×63 с реальными "
             "спектрами для удаления 1-2 канальных спайков.",
    )
    p_analyze.add_argument(
        "--narrow-peak-fwhm-ratio", type=float, default=0.3,
        help="(F-139) Пороговое соотношение measured/expected FWHM для "
             "отбраковки. Default 0.3.",
    )
    p_analyze.add_argument(
        "--background-auto",
        choices=["off", "suggest", "apply"],
        default="apply",
        help="(F-131 + F-135 / v1.17.7) Авто-поиск и применение "
             "подходящего фонового .spe в той же папке + типовых местах "
             "(data/averaged_backgrounds, *Фон*/) когда --background-path "
             "не задан. F-135: default=apply (фон ВСЕГДА вычитается при "
             "наличии подходящего кандидата). 'off' — отключить; "
             "'suggest' — только записать предложение в notes; "
             "'apply' (default) — вычесть лучшего кандидата. "
             "Эвристика: тот же детектор, совместимая геометрия, |Δt| ≤ 90 дн.",
    )
    p_analyze.add_argument(
        "--background-auto-max-days",
        type=int,
        default=90,
        help="(F-131) Максимальный временной разрыв между sample и фоновым "
             "кандидатом, дней. По умолчанию 90.",
    )
    p_analyze.add_argument(
        "--background-path",
        default=None,
        help="Path to a paired background spectrum (.spe) for "
             "energy-rebinned subtraction (F-58).",
    )
    p_analyze.add_argument(
        "--allow-stage2",
        action="store_true",
        help="Allow Stage-2 identification (Cs-137/Cs-134/Co-60/I-131).",
    )
    p_analyze.add_argument(
        "--allow-stage3",
        action="store_true",
        help="Allow Stage-3 identification (Na-22, Be-7, Am-241, ...).",
    )
    p_analyze.add_argument(
        "--recalibrate-on-anchor-disagreement",
        action="store_true",
        help="Opt into F-87 Step-5β bootstrap refit when the seeded "
             "anchors disagree with the stored energy calibration.",
    )
    p_analyze.add_argument(
        "--plot-dpi", type=int, default=120,
        help="DPI for the rendered PNG plots (default 120).",
    )
    p_analyze.add_argument("--no-json", action="store_true",
                           help="Skip writing the JSON report.")
    p_analyze.add_argument("--no-markdown", action="store_true",
                           help="Skip writing the Markdown report.")
    p_analyze.add_argument("--no-plots", action="store_true",
                           help="Skip writing the PNG plots.")
    p_analyze.add_argument("--no-html", action="store_true",
                           help="Skip writing the HTML report.")
    p_analyze.add_argument("--write-pdf", action="store_true",
                           help="Render the HTML report to PDF via Edge headless (F-114).")
    # F-159 / v1.18.21.0 — Technical PDF report (контракт сохранён как opt-in)
    # F-RPT-03 / v1.18.29 — default OFF; opt-in через --technical-pdf.
    p_analyze.add_argument(
        "--technical-pdf", action="store_true",
        help="(F-159, F-RPT-03) Сгенерировать Technical PDF "
             "(пошаговый walkthrough 11 шагов через reportlab). "
             "По умолчанию ВЫКЛЮЧЕНО (F-RPT-03 v1.18.29) — "
             "включается явным флагом.",
    )
    p_analyze.add_argument(
        "--no-technical-pdf", action="store_true",
        help="(back-compat) Explicit suppression of Technical PDF. "
             "Default уже OFF (F-RPT-03), флаг сохранён для скриптов.",
    )
    # F-160 / v1.18.19.0 — XML round-trip export (флаг сохранён как контракт)
    p_analyze.add_argument(
        "--export-becqmoni",
        choices=["off", "sample", "bg", "both"],
        default="off",
        help="(F-160) Экспортировать spectra в переносимый XML "
             "(ResultDataFile schema) после анализа. 'sample' — только "
             "образец (с F-145 калибровкой), 'bg' — только фон (re-read от "
             "--background), 'both' — оба, 'off' (default) — пропустить. "
             "Files: <output-dir>/<stem>_calibrated.bq.xml.",
    )
    # F-108 / D-19 — cost footer
    p_analyze.add_argument("--cost-tokens", type=int, default=None,
                           help="Estimated tokens consumed by the analysis "
                                "(displayed in the HTML footer).")
    p_analyze.add_argument("--cost-session-pct", type=str, default=None,
                           help='Estimated share of the free session (e.g. "65%%").')
    p_analyze.add_argument("--cost-detail", type=str, default=None,
                           help='Free-form cost-line detail (e.g. "v1.17.4 prod run").')
    p_analyze.add_argument(
        "--cost-session-token-budget", type=int, default=None,
        help="(F-132 / v1.17.7) Бюджет 5-часовой сессии в токенах "
             "(default 200000) для расчёта %% использования.",
    )
    # ── F-309 / v1.18.8 — opt-in флаги activity integrations v1.18.1..v1.18.4
    p_analyze.add_argument(
        "--enable-tcs-correction", action="store_true",
        help="(F-296 / v1.18.1) Auto-TCS correction для cascade нуклидов "
             "из CASCADE_PRESETS (Co-60, Eu-152, Ba-133, Bi-214, Bi-212). "
             "Использует F-295 P/T ratio (Gilmore Table 8.4). Default OFF.",
    )
    p_analyze.add_argument(
        "--tcs-detector-id", default="Gamma-1S",
        choices=["Gamma-1S", "3in3", "4in4", "NaI_63x63"],
        help="(F-295 / v1.18.1) Детектор для P/T ratio lookup. "
             "Default 'Gamma-1S' (63×63 NaI).",
    )
    p_analyze.add_argument(
        "--enable-cutshall-self-abs", action="store_true",
        help="(F-294 / v1.18.1) Cutshall analytic self-absorption fallback "
             "когда REF_GEOMETRY (F-122) не сработал. Использует NIST XCOM "
             "water μ/ρ. Требует --sample-density-g-cm3. Default OFF.",
    )
    p_analyze.add_argument(
        "--cutshall-path-cm", type=float, default=None,
        help="(F-294) Average γ-path в образце для Cutshall, см. "
             "Default — Marinelli 0.5L (1.75 cm).",
    )
    p_analyze.add_argument(
        "--cutshall-calib-density-g-cm3", type=float, default=1.0,
        help="(F-294) Плотность калибровочного source-материала. Default 1.0 "
             "(water equivalent).",
    )
    p_analyze.add_argument(
        "--enable-matrix-method", action="store_true",
        help="(F-297 / v1.18.2) Matrix-method simultaneous χ² solver для "
             "multi-nuclide activity (alternative к per-line weighted-mean). "
             "Требует ≥2 нуклидов и ≥2 пиков. Default OFF.",
    )
    p_analyze.add_argument(
        "--matrix-method-energy-tolerance-keV", type=float, default=1.0,
        help="(F-297) Допуск парирования peak ↔ library line для matrix "
             "method, кэВ. Default 1.0.",
    )
    p_analyze.set_defaults(func=_cmd_analyze)

    # ── summarize (alias) ────────────────────────────────────────
    p_summarize = sub.add_parser("summarize", help="Alias for `analyze`.")
    p_summarize.add_argument("spectrum")
    p_summarize.add_argument("--search-dir", action="append", default=[])
    p_summarize.set_defaults(func=_cmd_analyze_phase0)

    # ── rag (F-151..F-153 / v1.17.9) ─────────────────────────────
    p_rag = sub.add_parser(
        "rag",
        help="(F-151..F-153 / v1.17.9) RAG поиск по библиотеке знаний "
             "references/books/ — query/explain/cite/verify/rebuild.",
    )
    rag_sub = p_rag.add_subparsers(dest="rag_cmd", required=True)

    p_rag_q = rag_sub.add_parser("query", help="top-k чанков по запросу")
    p_rag_q.add_argument("query", help="свободный текст запроса (RU/EN)")
    p_rag_q.add_argument("-k", type=int, default=5,
                         help="число хитов (default 5)")
    p_rag_q.add_argument("--json", action="store_true",
                         help="вывод JSON")
    p_rag_q.set_defaults(func=_cmd_rag_query)

    p_rag_e = rag_sub.add_parser("explain", help="связное объяснение темы")
    p_rag_e.add_argument("topic")
    p_rag_e.add_argument("-k", type=int, default=3)
    p_rag_e.set_defaults(func=_cmd_rag_explain)

    p_rag_c = rag_sub.add_parser("cite",
                                 help="каноническая цитата top-1")
    p_rag_c.add_argument("topic")
    p_rag_c.set_defaults(func=_cmd_rag_cite)

    p_rag_v = rag_sub.add_parser("verify",
                                 help="проверка утверждения по библиотеке")
    p_rag_v.add_argument("claim")
    p_rag_v.add_argument("--min-score", type=float, default=1.5)
    p_rag_v.set_defaults(func=_cmd_rag_verify)

    p_rag_r = rag_sub.add_parser("rebuild",
                                 help="пересборка BM25 индекса")
    p_rag_r.add_argument("--no-corpus", action="store_true",
                         help="не подключать knowledge_corpus.json даже если есть")
    p_rag_r.set_defaults(func=_cmd_rag_rebuild)

    p_rag_s = rag_sub.add_parser("stats",
                                 help="статистика загруженного индекса")
    p_rag_s.set_defaults(func=_cmd_rag_stats)

    args = parser.parse_args(argv)
    return args.func(args)


# ──────────────────────────────────────────────────────────────────
# F-152 — RAG subcommand handlers
# ──────────────────────────────────────────────────────────────────

def _cmd_rag_query(args) -> int:
    from gamma.knowledge.rag_search import cli_query as _q
    # Reuse the rag_search CLI handler
    class _Args: pass
    a = _Args()
    a.query = args.query
    a.k = args.k
    a.json = args.json
    return _q(a)


def _cmd_rag_explain(args) -> int:
    from gamma.knowledge.rag_search import cli_explain as _e
    class _Args: pass
    a = _Args()
    a.topic = args.topic
    a.k = args.k
    return _e(a)


def _cmd_rag_cite(args) -> int:
    from gamma.knowledge.rag_search import cli_cite as _c
    class _Args: pass
    a = _Args()
    a.topic = args.topic
    return _c(a)


def _cmd_rag_verify(args) -> int:
    from gamma.knowledge.rag_search import cli_verify as _v
    class _Args: pass
    a = _Args()
    a.claim = args.claim
    a.min_score = args.min_score
    return _v(a)


def _cmd_rag_rebuild(args) -> int:
    from gamma.knowledge.rag_index import rebuild_index
    try:
        index, out_path = rebuild_index(
            include_corpus=not args.no_corpus,
        )
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    n_cur = sum(1 for d in index.docs if d.source_layer == "curated")
    n_cor = sum(1 for d in index.docs if d.source_layer == "corpus")
    print(f"OK: BM25 index → {out_path}")
    print(f"  total: {index.n_docs}  curated: {n_cur}  corpus: {n_cor}")
    print(f"  avgdl: {index.avgdl:.1f}  vocab: {len(index.df)}")
    return 0


def _cmd_rag_stats(args) -> int:
    from gamma.knowledge.rag_search import load_index
    try:
        index = load_index()
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    n_cur = sum(1 for d in index.docs if d.source_layer == "curated")
    n_cor = sum(1 for d in index.docs if d.source_layer == "corpus")
    print(f"docs:    {index.n_docs}")
    print(f"curated: {n_cur}")
    print(f"corpus:  {n_cor}")
    print(f"avgdl:   {index.avgdl:.1f}")
    print(f"vocab:   {len(index.df)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
