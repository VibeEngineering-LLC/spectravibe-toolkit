"""
One-call wrapper — analyze_and_report (F-86e / v1.15.0).

The full Step-1..11 pipeline + Step-11 report assembly behind a
single entry point. Intended for CLI use and notebook one-liners.

    from gamma.reporting import analyze_and_report

    artefacts = analyze_and_report(
        "data/Cs137_Marinelli.spe",
        output_dir="./out",
        sample_mass_kg=0.500,
    )
    print(artefacts["summary"])

The wrapper picks sensible defaults for typical Gamma-1S operation:
* ``complete_workflow=True`` (Round 5: deconvolution, activities, MDA)
* ``write_plots=True`` (PNG spectrum + multiplet overlays)
* ``write_markdown=True`` (so the Markdown picks up the plots)
* ``write_html=True`` (one self-contained file to email)

Override any of these through ``**kwargs`` — they pass through to
``analyze_lsrm_spe`` and ``build_report`` respectively.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, Optional, Tuple


# Keywords routed to analyze_lsrm_spe (the orchestrator). Anything
# else in kwargs goes to build_report.
_ORCHESTRATOR_KEYS = {
    "detector_type",
    "sigma_threshold",
    "fwhm_window_multiple",
    "allow_stage2",
    "allow_stage3",
    "auto_escalate",
    "user_confirmed_stage2_nuclides",
    "user_confirmed_stage3_nuclides",
    "background_path",
    "apply_deconvolution",
    "deconvolution_overlap_fwhm",
    "compute_activities",
    "sample_mass_kg",
    # F-122 / v1.17.6 — self-attenuation kwargs
    "sample_density_g_cm3",
    "matrix_composition",
    # F-129 / v1.17.7 — peak search method dispatch
    "peak_search_method",
    # F-139 / v1.17.7 — отбраковка узких пиков
    "filter_narrow_peaks",
    "narrow_peak_fwhm_ratio",
    # F-131 / v1.17.7 — auto-background search
    "background_auto",
    "background_auto_max_days",
    "compute_mda",
    "mda_suite_extra_lines_keV",
    "reference_datetime",
    "complete_workflow",
    "recalibrate_on_anchor_disagreement",
    "recalibration_threshold_fwhm",
    # F-309 / v1.18.8 — opt-in activity integration flags v1.18.1..v1.18.4.
    # Pass-through через analyze_lsrm_spe (F-308) → compute_activities_for_all.
    "enable_tcs_correction",
    "tcs_detector_id",
    "enable_cutshall_self_abs",
    "cutshall_path_cm",
    "cutshall_calib_density_g_cm3",
    "enable_matrix_method",
    "matrix_method_energy_tolerance_keV",
    # F-322 / v1.18.16 — opt-in F-96 bg-anchors в multiplet deconvolution.
    "enable_f96_bg_anchors",
}


def analyze_and_report(
    path: str,
    *,
    output_dir: Optional[str] = None,
    sample_mass_kg: Optional[float] = None,
    write_json: bool = True,
    write_markdown: bool = True,
    write_plots: bool = True,
    write_html: bool = True,
    write_pdf: bool = False,
    write_technical_pdf: bool = False,
    export_becqmoni: str = "off",
    return_summary: bool = True,
    complete_workflow: bool = True,
    plot_dpi: int = 120,
    report_stem: Optional[str] = None,
    cost_estimate: Optional[Dict[str, Any]] = None,
    passport_activity_Bq: Optional[Dict[str, float]] = None,
    **kwargs,
) -> Dict[str, Any]:
    """Run the full pipeline + assemble the Step-11 report in one call.

    Parameters
    ----------
    path : str
        Path to the spectrum file (.spe LSRM or .xml AtomSpectra).
    output_dir : str, optional
        Where to write the report artefacts. When ``None``, nothing is
        written and only the in-memory result + summary are returned.
    sample_mass_kg : float, optional
        Enables Bq/kg specific-activity derivation in the orchestrator.
    write_json, write_markdown, write_plots, write_html : bool
        Toggle individual artefacts. By default all four are on
        (``write_html=True`` is the noticeable change vs ``build_report``).
    return_summary : bool, default True
        Include the 3–8 line chat summary in the return dict.
    complete_workflow : bool, default True
        Turn on Round 5 in the orchestrator (deconvolution + activities
        + MDA per SKILL.md autonomous defaults).
    plot_dpi : int, default 120
    report_stem : str, optional
        Override the file stem for output filenames.
    **kwargs
        Routed to either ``analyze_lsrm_spe`` (if the key is recognized
        as an orchestrator kwarg) or to ``build_report``. Unknown keys
        cause a TypeError on the receiving function — keep them honest.

    Returns
    -------
    dict
        ``build_report`` output extended with one extra key:
        * ``result`` — the raw ``StagedAnalysisResult``.
    """
    # Lazy imports to keep wrapper module light
    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    from gamma.reporting.build import build_report

    # Split kwargs
    orch_kw: Dict[str, Any] = {}
    report_kw: Dict[str, Any] = {}
    for k, v in kwargs.items():
        if k in _ORCHESTRATOR_KEYS:
            orch_kw[k] = v
        else:
            report_kw[k] = v

    # Wrapper-level conveniences override kwargs if user passed them
    orch_kw.setdefault("complete_workflow", complete_workflow)
    if sample_mass_kg is not None:
        orch_kw.setdefault("sample_mass_kg", sample_mass_kg)
    # F-135 / v1.17.7 — production отчёты ВСЕГДА в режиме background
    # auto-apply. analyze_and_report — public API для генерации отчётов
    # (CLI / Jupyter / scripts). Тестовая инфраструктура использует
    # analyze_lsrm_spe напрямую с pipeline default "suggest" (back-compat
    # с synthetic-тестами v1.17.6).
    orch_kw.setdefault("background_auto", "apply")

    result = analyze_lsrm_spe(path, **orch_kw)

    # F-114 / v1.17.3 — stash sample_mass_kg on spec.extras so the
    # interactive renderer can surface it in the grid card without
    # round-tripping through the JSON.
    if sample_mass_kg is not None:
        try:
            result.spec.extras["sample_mass_kg"] = float(sample_mass_kg)
        except Exception:
            pass

    if cost_estimate is not None:
        report_kw.setdefault("cost_estimate", cost_estimate)
    # F-326 / v1.18.18.1 — passport activity comparison (всегда передаётся
    # в build_report, даже если None — тогда блок рендерится с explicit
    # deferred message + инструкцией как включить).
    # F-330 / v1.18.18.4 — если пользователь НЕ передал passport явно,
    # пробуем auto-извлечь из LSRM .spe COMMENT (entries собраны в
    # spec.extras["lsrm_passport"] reader-ом). Конвертируем Bq/kg → Bq
    # по массе образца + decay correction до даты измерения.
    auto_passport, auto_meta = _auto_passport_from_spec(
        result, explicit=passport_activity_Bq,
    )
    if passport_activity_Bq is None and auto_passport:
        report_kw.setdefault("passport_activity_Bq", auto_passport)
        report_kw.setdefault("passport_meta", auto_meta)
    else:
        report_kw.setdefault("passport_activity_Bq", passport_activity_Bq)
        if passport_activity_Bq is not None:
            report_kw.setdefault(
                "passport_meta",
                {"source": "explicit", "nuclides": list(passport_activity_Bq.keys())},
            )

    artefacts = build_report(
        result,
        output_dir=output_dir,
        write_json=write_json,
        write_markdown=write_markdown,
        write_plots=write_plots,
        write_html=write_html,
        write_technical_pdf=write_technical_pdf,
        return_summary=return_summary,
        report_stem=report_stem,
        plot_dpi=plot_dpi,
        **report_kw,
    )
    artefacts["result"] = result

    # F-114 / D-12 — optional PDF artefact via Edge headless.
    if write_pdf and artefacts.get("html"):
        try:
            from gamma.reporting.pdf_export import html_to_pdf
            pdf_path = html_to_pdf(artefacts["html"])
            if pdf_path:
                artefacts["pdf"] = pdf_path
        except Exception as e:
            artefacts.setdefault("warnings", []).append(
                f"pdf export failed: {type(e).__name__}: {e}"
            )

    # F-160 / v1.18.19.0 — Export откалиброванных спектров sample / bg в
    # BecqMoni/AtomSpectra XML формат для round-trip с AtomSpectra PRO и
    # сторонними BecqMoni-совместимыми инструментами. Writer уже реализован
    # в gamma.io.becqmoni_xml.write_becqmoni_xml; задача wrapper — связать
    # его с pipeline output.
    #
    # CLI: --export-becqmoni {off,sample,bg,both}  (default off)
    # • sample: записывает spec (включая F-145 калибровку если применилась).
    # • bg:     re-reads bg от background_path → пишет в BecqMoni XML.
    # • both:   оба файла (default для пользовательских пакетов отчёта).
    #
    # Output:
    #   <output-dir>/<sample-stem>_calibrated.bq.xml
    #   <output-dir>/<bg-stem>_calibrated.bq.xml
    #
    # Anonymisation contract (F-115): writer не вытаскивает абсолютные пути.
    if output_dir and export_becqmoni and export_becqmoni != "off":
        from pathlib import Path as _Path
        mode = str(export_becqmoni).lower().strip()
        valid_modes = {"sample", "bg", "both"}
        if mode not in valid_modes:
            artefacts.setdefault("warnings", []).append(
                f"export_becqmoni: invalid mode '{export_becqmoni}' — "
                f"ожидалось off/sample/bg/both, пропущено."
            )
        else:
            try:
                from gamma.io.becqmoni_xml import write_becqmoni_xml
                stem = report_stem or _Path(path).stem
                if mode in ("sample", "both"):
                    out = _Path(output_dir) / f"{stem}_calibrated.bq.xml"
                    write_becqmoni_xml(result.spec, str(out))
                    artefacts["becqmoni_sample"] = str(out)
                bg_path = orch_kw.get("background_path")
                if mode in ("bg", "both") and bg_path:
                    from gamma.io.readers import read_spectrum
                    bg_spec = read_spectrum(str(bg_path))
                    bg_stem = _Path(str(bg_path)).stem
                    out_bg = _Path(output_dir) / f"{bg_stem}_calibrated.bq.xml"
                    write_becqmoni_xml(bg_spec, str(out_bg))
                    artefacts["becqmoni_bg"] = str(out_bg)
                elif mode in ("bg", "both"):
                    # bg запрошен, но background_path не передан — warning
                    artefacts.setdefault("warnings", []).append(
                        "export_becqmoni=bg/both: background_path не передан, "
                        "bg XML не записан."
                    )
            except Exception as e:
                artefacts.setdefault("warnings", []).append(
                    f"becqmoni export failed: {type(e).__name__}: {e}"
                )

    return artefacts


def _auto_passport_from_spec(
    result, *, explicit: Optional[Dict[str, float]],
):
    """F-330 / v1.18.18.4 — auto-build passport dict из spec.extras.

    Returns (passport_dict_or_None, meta_dict).

    meta = {
        "source": "spe_comment" | "explicit" | "none",
        "raw_entries": [...],          # original parsed entries
        "mass_kg": 0.570 | None,
        "ref_dates": {nuc: "1997-05-30"},
        "decay_corrected": {nuc: True/False},
        "meas_date": "1999-08-04",
        "notes": "...",                # human-readable provenance
    }
    """
    meta = {
        "source": "none",
        "raw_entries": [],
        "mass_kg": None,
        "ref_dates": {},
        "decay_corrected": {},
        "meas_date": None,
        "notes": "",
    }
    if explicit is not None:
        meta["source"] = "explicit"
        return None, meta

    try:
        extras = getattr(result.spec, "extras", {}) or {}
        raw_entries = extras.get("lsrm_passport") or []
    except AttributeError:
        extras = {}
        raw_entries = []
    if not raw_entries:
        # F-369 / v1.18.24.4 — нет inline COMMENT-passport; попробуем
        # certificate .src по serial из filename.
        cert_data, cert_meta = _passport_from_certificate(result)
        if cert_data:
            return cert_data, cert_meta
        meta["notes"] = cert_meta.get("notes", "") or meta["notes"]
        return None, meta

    # Mass for Бк/кг → Бк conversion
    try:
        mass_kg = float(extras.get("lsrm_sample_mass_kg") or 0.0)
    except (TypeError, ValueError):
        mass_kg = 0.0
    meta["mass_kg"] = mass_kg if mass_kg > 0 else None

    # Measurement date for decay correction
    try:
        sd = result.spec.start_datetime
        meas_date = sd.date() if sd else None
    except AttributeError:
        meas_date = None
    meta["meas_date"] = meas_date.isoformat() if meas_date else None

    from gamma.io.lsrm_passport import decay_correct, half_life_seconds
    from datetime import date as _date

    out: Dict[str, float] = {}
    notes_lines = []
    for entry in raw_entries:
        nuc = entry.get("nuclide")
        if not nuc:
            continue
        is_spec = bool(entry.get("is_specific_activity"))
        val = float(entry.get("value") or 0.0)
        unit = str(entry.get("unit") or "")
        # Apply kBq → Bq prefix multiplier
        kilo_mult = 1000.0 if unit.lower().startswith(("kbq", "кбк")) else 1.0
        val_abs = val * kilo_mult

        # Convert specific activity → total Bq
        if is_spec:
            if mass_kg <= 0:
                notes_lines.append(
                    f"{nuc}: масса образца не известна — конверсия Бк/кг → Бк пропущена."
                )
                continue
            A0_Bq = val_abs * mass_kg
        else:
            A0_Bq = val_abs

        # Decay correction если есть ref_date и known half-life
        ref_iso = entry.get("reference_date")
        ref_dt = None
        if ref_iso:
            try:
                ref_dt = _date.fromisoformat(ref_iso)
            except ValueError:
                ref_dt = None
        decay_done = False
        A_final = A0_Bq
        if ref_dt and meas_date and half_life_seconds(nuc) is not None:
            corr = decay_correct(A0_Bq, nuc, ref_dt, meas_date)
            if corr is not None:
                A_final = corr
                decay_done = True

        meta["ref_dates"][nuc] = ref_iso
        meta["decay_corrected"][nuc] = decay_done
        out[nuc] = float(A_final)

    # Notes are always emitted (e.g. mass-missing reason) regardless of
    # whether passport dict came out non-empty — they explain WHY auto-
    # routing skipped или partial.
    if notes_lines:
        meta["notes"] = " ".join(notes_lines)

    if not out:
        # F-369 / v1.18.24.4 — fallback к certificate (.src) если все
        # COMMENT-entries оказались пустыми (например mass отсутствует
        # и Бк/кг→Бк конверсия пропущена).
        cert_data, cert_meta = _passport_from_certificate(result)
        if cert_data:
            return cert_data, cert_meta
        return None, meta

    meta["source"] = "spe_comment"
    meta["raw_entries"] = raw_entries
    return out, meta


# F-369 / v1.18.24.4 — .src certificate auto-discovery
# ────────────────────────────────────────────────────
# Pattern для serial number в LSRM source names: "420-7-17", "420-7-18", etc.
# Также допускаем "SRC-02" вариант.
import re as _re
_CERT_SERIAL_RE = _re.compile(r"(\d{3}-\d{1,2}-\d{1,3})")

# Mapping LSRM geometry-names (как в .src `[<source>]Geometry=`) →
# подстрока в filename (как в `Th232_420-7-17_Маринелли_0cm.spe`).
_CERT_GEOM_HINTS = {
    "маринелли": ("маринелли", "marinelli"),
    "петри": ("петри", "petri"),
    "дента": ("дента", "denta"),
    "точечный": ("точечный", "point", "5cm", "25cm"),
}


def _extract_serial_from_filename(path_str: str) -> Optional[str]:
    """Извлечь serial номер источника из имени .spe файла.
    Например `Th232_420-7-17_Маринелли_0cm.spe` → `420-7-17`.
    """
    name = os.path.basename(path_str or "")
    m = _CERT_SERIAL_RE.search(name)
    return m.group(1) if m else None


def _filename_geometry_hint(path_str: str) -> str:
    """Определить геометрию по подстроке в filename."""
    name = (os.path.basename(path_str or "")).lower()
    for canon, aliases in _CERT_GEOM_HINTS.items():
        for a in aliases:
            if a in name:
                return canon
    return ""


def _geometry_matches(cert_geom: str, hint: str) -> bool:
    """Проверка соответствия cert source.geometry → filename hint.
    Обе строки lowercase-normalized; hint содержит canonical alias."""
    if not hint:
        return True  # No hint → match any
    cg = (cert_geom or "").lower()
    aliases = _CERT_GEOM_HINTS.get(hint, (hint,))
    return any(a in cg for a in aliases)


def _passport_from_certificate(result) -> Tuple[Optional[Dict[str, float]], Dict]:
    """F-369 / v1.18.24.4 — fallback passport-извлечение из `.src`
    сертификатов LSRM Аспект, если в .spe COMMENT нет inline-данных.

    Алгоритм:
      1. Из filename .spe извлекается serial (например `420-7-17`).
      2. Из filename определяется geometry hint (Маринелли / Петри / ...).
      3. Сканируется `detectors/Gamma-1S/certificates/*.src`.
      4. В каждом сертификате — ищется source с подходящей геометрией
         и sub_source чьё имя содержит serial.
      5. Активности sub_source конвертируются Бк/кг → Бк по mass,
         корректируются на распад от reference_datetime до даты измерения.

    Returns (passport_dict_or_None, meta).
    """
    meta = {
        "source": "none",
        "raw_entries": [],
        "mass_kg": None,
        "ref_dates": {},
        "decay_corrected": {},
        "meas_date": None,
        "notes": "",
        "cert_file": None,
        "cert_source": None,
        "cert_subsource": None,
    }

    # Spec path. Spectrum class использует `source_path` атрибут;
    # legacy reader-ы могут использовать `path`. Пробуем оба.
    try:
        spec_path = (
            getattr(result.spec, "source_path", None)
            or getattr(result.spec, "path", None)
            or ""
        )
    except AttributeError:
        return None, meta
    spec_path = str(spec_path) if spec_path else ""
    if not spec_path:
        return None, meta

    serial = _extract_serial_from_filename(spec_path)
    if not serial:
        meta["notes"] = "serial-номер не найден в имени файла"
        return None, meta
    geom_hint = _filename_geometry_hint(spec_path)

    # Find repo root via spec.path → walk up to find detectors/Gamma-1S/certificates
    spec_p = Path(spec_path)
    cert_dir = None
    for parent in [spec_p, *spec_p.parents]:
        cand = parent / "detectors" / "Gamma-1S" / "certificates"
        if cand.is_dir():
            cert_dir = cand
            break
    if cert_dir is None:
        meta["notes"] = "папка certificates/ не найдена"
        return None, meta

    cert_files = sorted(cert_dir.glob("*.src"))
    if not cert_files:
        meta["notes"] = "В certificates/ нет .src файлов"
        return None, meta

    # Measurement date for decay correction
    try:
        sd = result.spec.start_datetime
        meas_dt = sd if sd else None
    except AttributeError:
        meas_dt = None
    if meas_dt is not None:
        meta["meas_date"] = meas_dt.date().isoformat()

    # Search
    from gamma.io.lsrm_src import read_certificate_file
    from gamma.io.lsrm_passport import decay_correct, half_life_seconds

    serial_norm = serial.replace("-", "").replace("_", "").replace(" ", "").lower()

    found = None
    for cert_path in cert_files:
        try:
            cert = read_certificate_file(cert_path)
        except Exception:
            continue
        for src_name, src in cert.sources.items():
            if not _geometry_matches(src.geometry, geom_hint):
                continue
            for sub in src.sub_sources:
                sub_norm = (
                    sub.name.replace("-", "").replace("_", "").lower()
                )
                if serial_norm in sub_norm:
                    found = (cert, src, sub)
                    break
            if found:
                break
        if found:
            break

    if not found:
        meta["notes"] = (
            f"в .src не найден sub-source с serial {serial!r} "
            f"и geometry-hint {geom_hint!r}"
        )
        return None, meta

    cert, src, sub = found
    meta["cert_file"] = Path(cert.file_path).name
    meta["cert_source"] = src.name
    meta["cert_subsource"] = sub.name
    # Mass for Бк/кг → Бк
    mass_g = sub.mass_g if sub.mass_g else src.mass_g
    mass_kg = (mass_g / 1000.0) if mass_g else 0.0
    meta["mass_kg"] = mass_kg if mass_kg > 0 else None

    unit_l = (src.activity_unit or "").lower()
    is_per_kg = "kg" in unit_l or "кг" in unit_l

    raw_entries = []
    out: Dict[str, float] = {}
    notes_lines = []

    for act in sub.activities:
        nuc = act.nuclide
        A_passport = float(act.A_Bq)
        if is_per_kg:
            if mass_kg <= 0:
                notes_lines.append(
                    f"{nuc}: масса не известна в .src — конверсия "
                    f"Бк/кг → Бк пропущена."
                )
                continue
            A0_Bq = A_passport * mass_kg
        else:
            A0_Bq = A_passport

        # Decay correction
        decay_done = False
        A_final = A0_Bq
        ref_dt = src.reference_datetime
        if ref_dt and meas_dt and half_life_seconds(nuc) is not None:
            corr = decay_correct(A0_Bq, nuc, ref_dt.date(), meas_dt.date())
            if corr is not None:
                A_final = corr
                decay_done = True

        meta["ref_dates"][nuc] = ref_dt.date().isoformat() if ref_dt else None
        meta["decay_corrected"][nuc] = decay_done
        out[nuc] = float(A_final)
        raw_entries.append({
            "nuclide": nuc,
            "value": A_passport,
            "unit": src.activity_unit or "Bq",
            "is_specific_activity": is_per_kg,
            "reference_date": (
                ref_dt.date().isoformat() if ref_dt else None
            ),
            "sigma_pct": act.sigma_pct,
            "confidence_sigma": act.confidence_sigma,
        })

    if notes_lines:
        meta["notes"] = " ".join(notes_lines)

    if not out:
        return None, meta

    meta["source"] = "cert_src"
    meta["raw_entries"] = raw_entries
    return out, meta


__all__ = ["analyze_and_report"]
