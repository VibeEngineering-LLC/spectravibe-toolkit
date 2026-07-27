# -*- coding: utf-8 -*-
from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""F-398 / v1.18.28 — run_skill.py: end-to-end orchestrator (no AI babysitting).

Один CLI/Python entry-point, прогоняющий весь skill pipeline для оператора:
от сырого .spe → полного report bundle (production + опционально V2 +
compare + bundle index).

См. AGENTS.md §5 для спецификации (фазы 0-8, defensive design, config,
тесты, документация).

Usage
-----
    # Минимум — один спектр
    python scripts/run_skill.py <spectrum.spe>

    # С фоном, массой и V2
    python scripts/run_skill.py <sample.spe> \\
        --background <bg.spe> \\
        --mass 0.5 \\
        --include-v2 \\
        --output-dir "$GAMMA_DEMO_REPORTS_DIR/v1_18_28_th232"

    # Батч
    python scripts/run_skill.py --batch "data/*.spe" \\
        --output-dir "$GAMMA_DEMO_REPORTS_DIR/auto"

    # Возобновить прерванный прогон
    python scripts/run_skill.py --resume "$GAMMA_DEMO_REPORTS_DIR/v1_18_28_th232"

Exit codes
----------
    0  success (all requested phases ok)
    1  partial — non-fatal phase failure(s); bundle usable но неполный
    2  fatal — input invalid / setup error; bundle не создан

Bundle layout
-------------
    <output-dir>/
    ├── sample/                            # Phase 3 PROD artefacts
    │   ├── <stem>_report.json
    │   ├── <stem>_report.md
    │   ├── <stem>_report.html
    │   ├── <stem>_technical_report.pdf
    │   ├── <stem>_plots/spectrum.png
    │   ├── <stem>_plots/multiplets/multiplet_*.png
    │   ├── <stem>_calibrated.bq.xml
    │   └── <bg_stem>_calibrated.bq.xml
    ├── sample_v2/                         # Phase 5 V2 artefacts (opt)
    │   └── ... (same shape as sample/)
    ├── v2_compare/                        # Phase 6 compare (opt)
    │   ├── compare_data.json
    │   └── v2_compare_report.html
    ├── index.html                         # Phase 7 navigation
    ├── run_skill_summary.json             # Phase 8 manifest
    ├── run_skill.log                      # structured log (all phases)
    └── .phases/                           # phase markers (resume)
        ├── phase_0.done
        ├── phase_1.done
        └── ...
"""

import argparse
import concurrent.futures as _futures
import datetime as _dt
import glob as _glob
import json
import logging
import os
import re
import sys
import time
import traceback
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

# ──────────────────────────────────────────────────────────────────
# Bootstrap: ensure gamma package is importable when run from anywhere
# ──────────────────────────────────────────────────────────────────

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

# ──────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────

SKILL_VERSION_FALLBACK = "unknown"
PHASE_NAMES = {
    0: "preflight",
    1: "load_calibrate",
    2: "prod_analyze",
    3: "prod_artefacts",
    4: "v2_analyze",
    5: "v2_artefacts",
    6: "compare",
    7: "index",
    8: "finalize",
}

# Regex set for filename token extraction.
# F-386 терминология не затрагивает имена файлов; здесь чисто parse.
_RX_MASS_KG = re.compile(r"(?P<n>\d+(?:[.,]\d+)?)\s*(?:kg|кг)", re.IGNORECASE)
_RX_MASS_G = re.compile(r"(?P<n>\d+(?:[.,]\d+)?)\s*(?:g|г)(?![а-яА-Яa-z])", re.IGNORECASE)
_RX_VOLUME_L = re.compile(r"(?P<n>\d+(?:[.,]\d+)?)\s*(?:l|л)(?![а-яА-Я])", re.IGNORECASE)
_RX_DISTANCE = re.compile(r"(?P<n>\d+)\s*cm", re.IGNORECASE)
_RX_BG_TOKENS = re.compile(
    r"(фон|bkg|background|bg_|_bg|закр|пустой)", re.IGNORECASE,
)
_RX_GEOMETRY = re.compile(
    r"(маринелли|marinelli|чашка|чаша|петри|petri|шар|дента|denta)",
    re.IGNORECASE,
)
_RX_DETECTOR = re.compile(
    # F2-A (2026-06-21): canonical NaI station is Gamma-1S (transliteration
    # of cyrillic «Гамма-1С»). Legacy ASCII «Gamma-1C» still matched for
    # backward-compat with archived filenames; both lower to "gamma-1s"/"gamma-1c"
    # respectively (post-merge canonicalisation is via aliases.canonicalize()).
    r"(gamma[-_]?1[cs]|hpge|labr3?|nai|cebr3?|czt|cdznte)",
    re.IGNORECASE,
)

# Geometry → default sample mass (kg) when filename token is absent.
# Empirically chosen to match HANDOFF / regen_demo_reports defaults.
GEOMETRY_DEFAULT_MASS = {
    "marinelli_1l": 0.5,
    "marinelli_0_5l": 0.25,
    "petri": 0.05,
    "chashka": 0.1,
    "sphere": 1.0,
    "unknown": 1.0,
}


# ──────────────────────────────────────────────────────────────────
# Skill version (best-effort import; tolerates missing module)
# ──────────────────────────────────────────────────────────────────

def _read_skill_version() -> str:
    """Прочитать SKILL_VERSION из json_report.py (line 21).

    Fail-safe: вернуть 'unknown' если файл не найден / не парсится.
    """
    candidate = REPO_ROOT / "scripts" / "gamma" / "reporting" / "json_report.py"
    if not candidate.exists():
        return SKILL_VERSION_FALLBACK
    try:
        text = candidate.read_text(encoding="utf-8", errors="replace")
        m = re.search(r'SKILL_VERSION\s*=\s*["\']([^"\']+)["\']', text)
        return m.group(1) if m else SKILL_VERSION_FALLBACK
    except OSError:
        return SKILL_VERSION_FALLBACK


# ──────────────────────────────────────────────────────────────────
# Config (YAML if available, else dict-defaults)
# ──────────────────────────────────────────────────────────────────

DEFAULT_CONFIG: Dict[str, Any] = {
    "analyze": {
        "sample_mass_kg": None,              # if None — auto-derive
        "sample_density_g_cm3": None,
        "export_becqmoni": "off",
        "full_report": True,
        "complete_workflow": True,
        "background_auto": "apply",
        "background_auto_max_days": 90,
        "peak_search_method": "mariscotti",
        "filter_narrow_peaks": None,
        "narrow_peak_fwhm_ratio": 0.3,
        "allow_stage2": True,
        # F-442 / v1.30.2 — Stage 3 (EXOTIC: Ga-67, Tc-99m, In-111, Mn-54,
        # Na-22, Be-7, Am-241, Eu-152, Ba-133, Co-57, Ti-44, Sc-44) **DEFAULT OFF**.
        # Doctrine (ern_set.py:30-38): Stage 3 nuclides «should NEVER be
        # auto-proposed; require explicit user confirmation». Override:
        # env var GAMMA_ALLOW_STAGE3=1 (или --allow-stage3 CLI flag).
        # Trigger: Ga-67 false-positive на природном «Камне с Ra-226» 2026-06-13.
        "allow_stage3": False,
    },
    "multiplet": {
        # Defaults are owned by Agent A (deconvolve.py). These keys mirror
        # the lock so config can override per-detector. NOT changing the
        # defaults below — they reproduce the Agent A contract.
        "unresolved_separation_fwhm_factor": 1.1,   # F-387.1 NaI
        "max_components_per_cluster": 3,            # F-387.1
        "min_significance_snr": 3.0,                # F-391
        "min_significance_snr_singleton": 5.0,      # F-391
    },
    "v2": {
        "enabled": False,
        "compare": True,
    },
    "artefacts": {
        "json": True,
        "markdown": True,
        "html": True,
        # F-RPT-03 / v1.18.29 — Technical PDF (F-159) OFF by default.
        # Re-enable manually via --config or by editing user config.
        "technical_pdf": False,
        "plots": True,
        # F-RPT-04 / v1.18.29 — BecqMoni XML export OFF by default
        # (некоторые входные .spe файлы конвертируются с потерей
        # n_channels — audit отложен на v1.18.30). Re-enable вручную.
        "xml_bq": False,
    },
    "output": {
        "base_dir": None,        # None = $GAMMA_DEMO_REPORTS_DIR/<stem>
        "subdir_sample": "sample",
        "subdir_sample_v2": "sample_v2",
        "subdir_compare": "v2_compare",
    },
}


def _deep_merge(dst: Dict[str, Any], src: Dict[str, Any]) -> Dict[str, Any]:
    """Mutating deep-merge of src into dst (dict values recurse, others overwrite)."""
    for k, v in src.items():
        if k in dst and isinstance(dst[k], dict) and isinstance(v, dict):
            _deep_merge(dst[k], v)
        else:
            dst[k] = v
    return dst


def _load_config(path: Optional[Path]) -> Dict[str, Any]:
    """Load user config and merge over DEFAULT_CONFIG.

    PyYAML необязателен — если файл с расширением .json, парсим как JSON.
    Если PyYAML нет и расширение .yaml/.yml — кидаем явную ошибку.
    """
    cfg = json.loads(json.dumps(DEFAULT_CONFIG))  # deep clone
    if path is None:
        return cfg
    if not path.exists():
        raise FileNotFoundError(f"Config not found: {path}")
    text = path.read_text(encoding="utf-8")
    user_cfg: Dict[str, Any]
    if path.suffix.lower() == ".json":
        user_cfg = json.loads(text)
    else:
        try:
            import yaml  # type: ignore
        except ImportError as e:
            raise RuntimeError(
                f"PyYAML required for {path.suffix} config; "
                f"install `pip install pyyaml` or convert to JSON."
            ) from e
        user_cfg = yaml.safe_load(text) or {}
    if not isinstance(user_cfg, dict):
        raise ValueError(f"Config root must be a mapping, got {type(user_cfg).__name__}")
    return _deep_merge(cfg, user_cfg)


# ──────────────────────────────────────────────────────────────────
# Bundle layout + phase markers
# ──────────────────────────────────────────────────────────────────

@dataclass
class BundleLayout:
    """Все пути bundle, выведенные из base_dir."""
    base: Path
    sample: Path
    sample_v2: Path
    compare: Path
    index_html: Path
    summary_json: Path
    log_file: Path
    phases_dir: Path

    @classmethod
    def from_base(cls, base: Path, cfg: Dict[str, Any]) -> "BundleLayout":
        out = cfg["output"]
        return cls(
            base=base,
            sample=base / out["subdir_sample"],
            sample_v2=base / out["subdir_sample_v2"],
            compare=base / out["subdir_compare"],
            index_html=base / "index.html",
            summary_json=base / "run_skill_summary.json",
            log_file=base / "run_skill.log",
            phases_dir=base / ".phases",
        )

    def ensure_dirs(self) -> None:
        for p in (self.base, self.phases_dir):
            p.mkdir(parents=True, exist_ok=True)

    def marker(self, phase: int) -> Path:
        return self.phases_dir / f"phase_{phase}.done"


# ──────────────────────────────────────────────────────────────────
# Logging setup (per-bundle file + console)
# ──────────────────────────────────────────────────────────────────

def _setup_logger(layout: BundleLayout, verbose: bool, quiet: bool) -> logging.Logger:
    """Per-bundle logger. Файл всегда DEBUG, консоль — INFO/WARNING по флагам.

    Каждый прогон получает свежий logger (per-bundle); хендлеры предыдущих
    bundle закрываются. Это позволяет батч-режиму не путать логи.
    """
    logger = logging.getLogger(f"run_skill.{layout.base.name}")
    # Очистка от прошлого прогона (важно для --resume и для batch).
    for h in list(logger.handlers):
        logger.removeHandler(h)
        try:
            h.close()
        except Exception:
            pass
    logger.setLevel(logging.DEBUG)
    logger.propagate = False

    layout.ensure_dirs()
    fh = logging.FileHandler(layout.log_file, mode="a", encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    ))
    logger.addHandler(fh)

    if not quiet:
        ch = logging.StreamHandler(sys.stderr)
        ch.setLevel(logging.DEBUG if verbose else logging.INFO)
        ch.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
        logger.addHandler(ch)
    return logger


# ──────────────────────────────────────────────────────────────────
# Metadata extraction from filename + spectrum header
# ──────────────────────────────────────────────────────────────────

@dataclass
class SpectrumMetadata:
    stem: str
    detector_hint: Optional[str] = None
    geometry_hint: Optional[str] = None
    distance_cm: Optional[int] = None
    sample_mass_kg: Optional[float] = None
    is_background_candidate: bool = False
    source: str = "filename"

    @classmethod
    def from_path(cls, path: Path) -> "SpectrumMetadata":
        stem = path.stem
        name = stem.lower()

        meta = cls(stem=stem)

        if (m := _RX_DETECTOR.search(name)):
            meta.detector_hint = m.group(1).lower()
        if (m := _RX_GEOMETRY.search(name)):
            meta.geometry_hint = m.group(1).lower()
        if (m := _RX_DISTANCE.search(name)):
            try:
                meta.distance_cm = int(m.group("n"))
            except ValueError:
                pass

        # Mass: kg explicit → kg; g → /1000; else None (fall through to
        # geometry default later).
        if (m := _RX_MASS_KG.search(name)):
            meta.sample_mass_kg = _parse_float(m.group("n"))
        elif (m := _RX_MASS_G.search(name)):
            g = _parse_float(m.group("n"))
            if g is not None:
                meta.sample_mass_kg = g / 1000.0

        if _RX_BG_TOKENS.search(name):
            meta.is_background_candidate = True
        return meta


def _parse_float(s: str) -> Optional[float]:
    try:
        return float(s.replace(",", "."))
    except (ValueError, AttributeError):
        return None


def _default_mass_for_geometry(geo_hint: Optional[str]) -> float:
    """Best-effort геометрия → масса. Используется только при отсутствии
    явного token в имени и при отсутствии --mass в CLI."""
    if geo_hint is None:
        return GEOMETRY_DEFAULT_MASS["unknown"]
    g = geo_hint.lower()
    if "маринелли" in g or "marinelli" in g:
        return GEOMETRY_DEFAULT_MASS["marinelli_1l"]
    if "петри" in g or "petri" in g:
        return GEOMETRY_DEFAULT_MASS["petri"]
    if "чашка" in g or "чаша" in g:
        return GEOMETRY_DEFAULT_MASS["chashka"]
    if "шар" in g:
        return GEOMETRY_DEFAULT_MASS["sphere"]
    return GEOMETRY_DEFAULT_MASS["unknown"]


def _read_spec_mass_kg(spectrum_path: Path) -> Optional[tuple]:
    """BUG-1 / 2026-06-02 — light-read sample mass from the spectrum file.

    Returns (mass_kg, uncertainty_kg_or_None) if the underlying format
    exposes SAMPLEMASS (currently only LSRM .spe via the typed
    `Spectrum.sample_mass_kg` field); returns None on any error or when
    the field is absent. Pure read; no caching, no side effects beyond
    the reader itself.
    """
    try:
        from gamma.io.readers import read_spectrum
        spec = read_spectrum(str(spectrum_path))
    except Exception:
        return None
    m = getattr(spec, "sample_mass_kg", None)
    if m is None:
        return None
    try:
        m_f = float(m)
    except (TypeError, ValueError):
        return None
    if not (m_f > 0.0):
        return None
    u = getattr(spec, "sample_mass_uncertainty_kg", None)
    try:
        u_f = float(u) if u is not None else None
    except (TypeError, ValueError):
        u_f = None
    return (m_f, u_f)


def _extract_embedded_background(
    spec: Any,
    layout: "BundleLayout",
    stem: str,
    logger: logging.Logger,
) -> Optional[Path]:
    """F-397.1: вытащить embedded background spectrum из исходного файла
    в отдельный `.spe` под bundle и вернуть путь.

    Используется, когда:
    • spec.background_embedded is not None (N42 RadMeasurement
      `measurementClassCode="Background"` или AtomSpectra
      `<BackgroundEnergySpectrum>`)
    • ctx.background не задан (ни explicit --background, ни auto-detect)

    После извлечения путь подаётся в analyze_and_report как обычный
    --background-path; F-397 в staged_pipeline отрабатывает штатно
    (peak detection + bg view в HTML), без изменений в Agent A зоне.

    Контракт F-397.1 файла:
    • Пишется в `<bundle>/.embedded_bg/<stem>_embedded_bg.spe` (скрытая
      подпапка, чтобы не путать с production артефактами).
    • LSRM format (через write_lsrm_spe) — переносимый, читается
      reader-ом без дополнительных hooks.
    • Существующий файл перезаписывается — идемпотентно при resume.
    """
    if spec is None or getattr(spec, "background_embedded", None) is None:
        return None

    bg_spec = spec.background_embedded
    out_dir = layout.base / ".embedded_bg"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{stem}_embedded_bg.spe"
    try:
        from gamma.io.lsrm_spe import write_lsrm_spe
        # F-397.1 marker in extra_header to make extracted files traceable.
        write_lsrm_spe(
            bg_spec,
            str(out_path),
            extra_header={"COMMENT_F397_1": "extracted from embedded section by run_skill"},
        )
    except Exception as e:
        logger.warning(
            f"F-397.1: embedded bg extraction failed "
            f"({type(e).__name__}: {e}); fallback to no-bg"
        )
        return None
    logger.info(f"F-397.1: extracted embedded background → {out_path.name}")
    return out_path


def _auto_detect_background(
    sample_path: Path,
    logger: logging.Logger,
    max_days: int = 90,
) -> Optional[Path]:
    """Auto-detect background sibling. Используется только если --auto-detect-bg
    и не передан --background. Эвристика мягкая — алгоритм-приоритет всё равно
    отдаётся background_auto=apply в pipeline (F-131/F-135).
    """
    parent = sample_path.parent
    if not parent.is_dir():
        return None
    candidates = []
    for p in parent.iterdir():
        if not p.is_file():
            continue
        if p == sample_path:
            continue
        if p.suffix.lower() not in (".spe", ".chn", ".n42", ".mca"):
            continue
        if _RX_BG_TOKENS.search(p.name.lower()):
            candidates.append(p)
    if not candidates:
        logger.debug(f"auto-detect bg: no sibling candidates in {parent}")
        return None

    # Pick closest mtime to sample's mtime.
    try:
        sample_mtime = sample_path.stat().st_mtime
    except OSError:
        return candidates[0]
    candidates.sort(key=lambda p: abs(p.stat().st_mtime - sample_mtime))
    pick = candidates[0]
    logger.info(f"auto-detected background: {pick.name}")
    return pick


# ──────────────────────────────────────────────────────────────────
# Run context (passed through phases)
# ──────────────────────────────────────────────────────────────────

@dataclass
class RunContext:
    spectrum: Path
    background: Optional[Path]
    metadata: SpectrumMetadata
    cfg: Dict[str, Any]
    layout: BundleLayout
    logger: logging.Logger
    skill_version: str
    include_v2: bool
    skip_artefacts: Dict[str, bool] = field(default_factory=dict)
    phase_results: Dict[int, Dict[str, Any]] = field(default_factory=dict)
    started_at: float = field(default_factory=time.monotonic)


@dataclass
class PhaseResult:
    phase: int
    name: str
    status: str           # "ok" | "skipped" | "failed" | "resumed"
    elapsed_s: float
    detail: Dict[str, Any] = field(default_factory=dict)
    error: Optional[str] = None


@contextmanager
def _phase_timer():
    t0 = time.monotonic()
    yield lambda: time.monotonic() - t0


def _record_phase(
    ctx: RunContext,
    phase: int,
    result: PhaseResult,
) -> None:
    ctx.phase_results[phase] = asdict(result)
    if result.status in ("ok", "resumed"):
        marker = ctx.layout.marker(phase)
        marker.write_text(
            json.dumps(asdict(result), ensure_ascii=False, indent=2),
            encoding="utf-8",
        )


def _phase_already_done(ctx: RunContext, phase: int) -> bool:
    return ctx.layout.marker(phase).exists()


def _load_prior_phase(ctx: RunContext, phase: int) -> Optional[Dict[str, Any]]:
    p = ctx.layout.marker(phase)
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


# ──────────────────────────────────────────────────────────────────
# Phase implementations
# ──────────────────────────────────────────────────────────────────

def phase_0_preflight(ctx: RunContext, resume: bool) -> PhaseResult:
    """Validate inputs, env, demo_root, skill version.

    F-384 contract: если output-dir не передан, используется
    GAMMA_DEMO_REPORTS_DIR/<stem>/.
    """
    name = PHASE_NAMES[0]
    with _phase_timer() as t:
        if resume and _phase_already_done(ctx, 0):
            prior = _load_prior_phase(ctx, 0) or {}
            ctx.logger.info(f"phase 0 ({name}): resumed from marker")
            return PhaseResult(
                phase=0, name=name, status="resumed",
                elapsed_s=t(), detail=prior.get("detail", {}),
            )

        if not ctx.spectrum.exists():
            return PhaseResult(
                phase=0, name=name, status="failed",
                elapsed_s=t(),
                error=f"spectrum not found: {ctx.spectrum}",
            )

        if ctx.background is not None and not ctx.background.exists():
            return PhaseResult(
                phase=0, name=name, status="failed",
                elapsed_s=t(),
                error=f"background not found: {ctx.background}",
            )

        # Demo root visibility — info only; analyze_and_report uses
        # output_dir explicitly anyway.
        env_root = os.environ.get("GAMMA_DEMO_REPORTS_DIR")
        detail = {
            "spectrum": str(ctx.spectrum),
            "background": str(ctx.background) if ctx.background else None,
            "metadata": asdict(ctx.metadata),
            "skill_version": ctx.skill_version,
            "include_v2": ctx.include_v2,
            "GAMMA_DEMO_REPORTS_DIR": env_root,
            "bundle_base": str(ctx.layout.base),
            "python": sys.version.split()[0],
        }
        ctx.logger.info(
            f"phase 0 ({name}): ok — "
            f"spectrum={ctx.spectrum.name}, "
            f"bg={ctx.background.name if ctx.background else 'none'}, "
            f"v={ctx.skill_version}"
        )
        return PhaseResult(
            phase=0, name=name, status="ok",
            elapsed_s=t(), detail=detail,
        )


def phase_1_load_calibrate(ctx: RunContext, resume: bool) -> PhaseResult:
    """Light spectrum read for header/calibration sanity check.

    Не вызываем тяжёлый pipeline — это делает Phase 2. Здесь только проверка,
    что файл читается и spec.live_time/.real_time/.channels валидны.
    """
    name = PHASE_NAMES[1]
    with _phase_timer() as t:
        if resume and _phase_already_done(ctx, 1):
            prior = _load_prior_phase(ctx, 1) or {}
            ctx.logger.info(f"phase 1 ({name}): resumed from marker")
            return PhaseResult(
                phase=1, name=name, status="resumed",
                elapsed_s=t(), detail=prior.get("detail", {}),
            )

        try:
            from gamma.io.readers import read_spectrum
            spec = read_spectrum(str(ctx.spectrum))
        except Exception as e:
            ctx.logger.error(f"phase 1: read failed: {type(e).__name__}: {e}")
            return PhaseResult(
                phase=1, name=name, status="failed",
                elapsed_s=t(),
                error=f"read_spectrum: {type(e).__name__}: {e}",
            )

        detail = {
            "channels": int(getattr(spec, "n_channels", 0) or 0),
            "live_time_s": float(getattr(spec, "live_time", 0) or 0.0),
            "real_time_s": float(getattr(spec, "real_time", 0) or 0.0),
            "stored_calibration_present": bool(
                getattr(spec, "energy_calibration", None) is not None
            ),
            "background_link": getattr(spec, "background_link", None),
            "background_embedded_present": (
                getattr(spec, "background_embedded", None) is not None
            ),
            "background_source": "external" if ctx.background else "none",
        }
        # F-397.1 — extract embedded bg into a sibling file if present
        # AND no external bg supplied; downstream pipeline then uses it
        # exactly like a normal --background-path (F-397 trigger fires).
        if (detail["background_embedded_present"]
                and ctx.background is None):
            extracted = _extract_embedded_background(
                spec, ctx.layout, ctx.metadata.stem, ctx.logger,
            )
            if extracted is not None:
                ctx.background = extracted
                detail["background_source"] = "embedded_extracted"
                detail["background_extracted_path"] = str(extracted)
        # Dead time hint, surfaced for the diagnostic block.
        if detail["real_time_s"] > 0:
            dt_pct = 100.0 * (1.0 - detail["live_time_s"] / detail["real_time_s"])
            detail["dead_time_pct"] = round(dt_pct, 2)
        ctx.logger.info(
            f"phase 1 ({name}): ok — channels={detail['channels']}, "
            f"live={detail['live_time_s']:.0f}s, real={detail['real_time_s']:.0f}s, "
            f"bg_source={detail['background_source']}"
        )
        return PhaseResult(
            phase=1, name=name, status="ok",
            elapsed_s=t(), detail=detail,
        )


def _build_orch_kwargs(ctx: RunContext) -> Dict[str, Any]:
    """Свернуть config + metadata в kwargs для analyze_and_report."""
    cfg = ctx.cfg
    a = cfg["analyze"]
    art = cfg["artefacts"]

    # Mass: priority — CLI/config explicit > filename token > .spe
    # SAMPLEMASS field > geometry default. BUG-1 / 2026-06-02 fix:
    # previously SAMPLEMASS from the .spe header was ignored at this
    # layer, so geometry-default (0.5 kg for Маринелли) was fed to
    # analyze_and_report — inflating Бк/кг by SAMPLEMASS / 0.5×.
    # Now we light-read the file via _read_spec_mass_kg() so the typed
    # `Spectrum.sample_mass_kg` (populated from LSRM SAMPLEMASS) wins
    # over the geometry default. Only when both the CLI flag AND the
    # .spe field are absent does the F-378 default-warning fire.
    mass = a.get("sample_mass_kg")
    if mass is None:
        mass = ctx.metadata.sample_mass_kg
    if mass is None:
        spec_pair = _read_spec_mass_kg(ctx.spectrum)
        if spec_pair is not None:
            mass, u_kg = spec_pair
            u_txt = (
                f", uncertainty {u_kg:.3f} кг" if u_kg is not None else ""
            )
            ctx.logger.info(
                f"масса из .spe SAMPLEMASS = {mass:.3f} кг{u_txt} "
                f"(приоритет: CLI > имя файла > .spe SAMPLEMASS > default)"
            )
    if mass is None:
        mass = _default_mass_for_geometry(ctx.metadata.geometry_hint)
        ctx.logger.warning(
            f"F-378: sample_mass_kg не задан (--sample-mass-kg отсутствует, "
            f"имя файла не содержит массу, .spe SAMPLEMASS не найден); "
            f"принят default {mass} kg для геометрии "
            f"'{ctx.metadata.geometry_hint or 'unknown'}'. "
            f"Удельная активность может быть искажена. "
            f"Передайте --sample-mass-kg или используйте .spe-файл с "
            f"корректным SAMPLEMASS."
        )

    export = "off"
    # F-RPT-04 / v1.18.29 — BecqMoni export OFF by default. Both
    # xml_bq artefact toggle AND analyze.export_becqmoni must be
    # opt-in for the wrapper to receive a non-off mode.
    if art.get("xml_bq", False):
        export = a.get("export_becqmoni", "off")

    kwargs: Dict[str, Any] = dict(
        sample_mass_kg=float(mass),
        write_json=bool(art.get("json", True)),
        write_markdown=bool(art.get("markdown", True)),
        write_html=bool(art.get("html", True)),
        write_plots=bool(art.get("plots", True)),
        # F-RPT-03 / v1.18.29 — default-off если ключ отсутствует.
        write_technical_pdf=bool(art.get("technical_pdf", False)),
        export_becqmoni=export,
        complete_workflow=bool(a.get("complete_workflow", True)),
        background_auto=str(a.get("background_auto", "apply")),
        background_auto_max_days=int(a.get("background_auto_max_days", 90)),
        peak_search_method=str(a.get("peak_search_method", "mariscotti")),
        narrow_peak_fwhm_ratio=float(a.get("narrow_peak_fwhm_ratio", 0.3)),
        allow_stage2=bool(a.get("allow_stage2", True)),
        allow_stage3=bool(a.get("allow_stage3", False)),
        # post-#118 — same gate как phase_7_index. Не orchestrator-кей →
        # wrapper маршрутизирует в build_report → build_html_report →
        # render_interactive_html, чтобы back-nav emit'ился только когда
        # bundle landing генерируется.
        bundle_index=bool(
            cfg.get("reports", {}).get("bundle_index", False)
        ),
    )
    if a.get("sample_density_g_cm3") is not None:
        kwargs["sample_density_g_cm3"] = float(a["sample_density_g_cm3"])
    if a.get("filter_narrow_peaks") is not None:
        kwargs["filter_narrow_peaks"] = bool(a["filter_narrow_peaks"])
    if ctx.background is not None:
        kwargs["background_path"] = str(ctx.background)
    return kwargs


def phase_2_prod_analyze(ctx: RunContext, resume: bool) -> PhaseResult:
    """PROD: analyze_and_report → bundle/sample/."""
    name = PHASE_NAMES[2]
    out_dir = ctx.layout.sample
    with _phase_timer() as t:
        if resume and _phase_already_done(ctx, 2):
            ctx.logger.info(f"phase 2 ({name}): resumed (marker present)")
            prior = _load_prior_phase(ctx, 2) or {}
            return PhaseResult(
                phase=2, name=name, status="resumed",
                elapsed_s=t(), detail=prior.get("detail", {}),
            )

        out_dir.mkdir(parents=True, exist_ok=True)
        try:
            from gamma.reporting import analyze_and_report
        except ImportError as e:
            return PhaseResult(
                phase=2, name=name, status="failed",
                elapsed_s=t(),
                error=f"gamma.reporting import failed: {e}",
            )

        kwargs = _build_orch_kwargs(ctx)
        ctx.logger.info(
            f"phase 2 ({name}): analyze_and_report → {out_dir}"
        )
        try:
            artefacts = analyze_and_report(
                str(ctx.spectrum),
                output_dir=str(out_dir),
                **kwargs,
            )
        except Exception as e:
            tb = traceback.format_exc()
            ctx.logger.error(f"phase 2: {type(e).__name__}: {e}\n{tb}")
            return PhaseResult(
                phase=2, name=name, status="failed",
                elapsed_s=t(),
                error=f"{type(e).__name__}: {e}",
                detail={"traceback_tail": tb.splitlines()[-3:]},
            )

        # Detach the heavy 'result' object before serializing.
        detail = {k: v for k, v in artefacts.items()
                  if k not in ("result", "html_text")}
        # Convert Path-like values to str for JSON serialization.
        detail = _stringify_paths(detail)
        ctx.logger.info(
            f"phase 2 ({name}): ok — "
            f"artefacts: {sorted(k for k in detail if not k.startswith('_'))}"
        )
        return PhaseResult(
            phase=2, name=name, status="ok",
            elapsed_s=t(), detail=detail,
        )


def phase_3_prod_artefacts(ctx: RunContext, resume: bool) -> PhaseResult:
    """Verify PROD artefacts on disk. analyze_and_report уже пишет их в Phase 2;
    эта фаза — гейт, проверяющий presence обязательных файлов."""
    name = PHASE_NAMES[3]
    with _phase_timer() as t:
        if resume and _phase_already_done(ctx, 3):
            ctx.logger.info(f"phase 3 ({name}): resumed (marker present)")
            return PhaseResult(
                phase=3, name=name, status="resumed",
                elapsed_s=t(), detail=(_load_prior_phase(ctx, 3) or {}).get("detail", {}),
            )

        missing = _check_required_artefacts(ctx.layout.sample, ctx.cfg)
        if missing:
            ctx.logger.warning(f"phase 3: missing artefacts: {missing}")
            return PhaseResult(
                phase=3, name=name, status="failed",
                elapsed_s=t(),
                error=f"missing artefacts: {missing}",
                detail={"missing": missing},
            )
        files = _list_artefacts(ctx.layout.sample)
        ctx.logger.info(f"phase 3 ({name}): ok — {len(files)} files present")
        return PhaseResult(
            phase=3, name=name, status="ok",
            elapsed_s=t(), detail={"files": files},
        )


def phase_4_v2_analyze(ctx: RunContext, resume: bool) -> PhaseResult:
    """V2: analyze_and_report_v2 → bundle/sample_v2/. Skipped if --include-v2 off."""
    name = PHASE_NAMES[4]
    with _phase_timer() as t:
        if not ctx.include_v2:
            ctx.logger.info(f"phase 4 ({name}): skipped (V2 disabled)")
            return PhaseResult(
                phase=4, name=name, status="skipped",
                elapsed_s=t(), detail={"reason": "include_v2=False"},
            )
        if resume and _phase_already_done(ctx, 4):
            ctx.logger.info(f"phase 4 ({name}): resumed (marker present)")
            return PhaseResult(
                phase=4, name=name, status="resumed",
                elapsed_s=t(),
                detail=(_load_prior_phase(ctx, 4) or {}).get("detail", {}),
            )

        out_dir = ctx.layout.sample_v2
        out_dir.mkdir(parents=True, exist_ok=True)
        try:
            from gamma.experimental.v2_integration import analyze_and_report_v2
        except ImportError as e:
            return PhaseResult(
                phase=4, name=name, status="failed",
                elapsed_s=t(),
                error=f"v2_integration import failed: {e}",
            )

        kwargs = _build_orch_kwargs(ctx)
        ctx.logger.info(
            f"phase 4 ({name}): analyze_and_report_v2 → {out_dir}"
        )
        try:
            artefacts = analyze_and_report_v2(
                str(ctx.spectrum),
                output_dir=str(out_dir),
                **kwargs,
            )
        except Exception as e:
            tb = traceback.format_exc()
            ctx.logger.error(f"phase 4: {type(e).__name__}: {e}\n{tb}")
            return PhaseResult(
                phase=4, name=name, status="failed",
                elapsed_s=t(),
                error=f"{type(e).__name__}: {e}",
                detail={"traceback_tail": tb.splitlines()[-3:]},
            )

        detail = _stringify_paths(
            {k: v for k, v in artefacts.items()
             if k not in ("result", "html_text")}
        )
        ctx.logger.info(
            f"phase 4 ({name}): ok — "
            f"artefacts: {sorted(k for k in detail if not k.startswith('_'))}"
        )
        return PhaseResult(
            phase=4, name=name, status="ok",
            elapsed_s=t(), detail=detail,
        )


def phase_5_v2_artefacts(ctx: RunContext, resume: bool) -> PhaseResult:
    """Verify V2 artefacts on disk."""
    name = PHASE_NAMES[5]
    with _phase_timer() as t:
        if not ctx.include_v2:
            return PhaseResult(
                phase=5, name=name, status="skipped",
                elapsed_s=t(), detail={"reason": "include_v2=False"},
            )
        if resume and _phase_already_done(ctx, 5):
            ctx.logger.info(f"phase 5 ({name}): resumed (marker present)")
            return PhaseResult(
                phase=5, name=name, status="resumed",
                elapsed_s=t(),
                detail=(_load_prior_phase(ctx, 5) or {}).get("detail", {}),
            )
        missing = _check_required_artefacts(ctx.layout.sample_v2, ctx.cfg)
        if missing:
            return PhaseResult(
                phase=5, name=name, status="failed",
                elapsed_s=t(),
                error=f"missing V2 artefacts: {missing}",
                detail={"missing": missing},
            )
        files = _list_artefacts(ctx.layout.sample_v2)
        ctx.logger.info(f"phase 5 ({name}): ok — {len(files)} files")
        return PhaseResult(
            phase=5, name=name, status="ok",
            elapsed_s=t(), detail={"files": files},
        )


def phase_6_compare(ctx: RunContext, resume: bool) -> PhaseResult:
    """V2 vs PROD compare. Calls gen_v2_compare_th232.py on bundle base."""
    name = PHASE_NAMES[6]
    with _phase_timer() as t:
        if not (ctx.include_v2 and ctx.cfg["v2"].get("compare", True)):
            return PhaseResult(
                phase=6, name=name, status="skipped",
                elapsed_s=t(),
                detail={"reason": "compare disabled or v2 not run"},
            )
        if resume and _phase_already_done(ctx, 6):
            return PhaseResult(
                phase=6, name=name, status="resumed",
                elapsed_s=t(),
                detail=(_load_prior_phase(ctx, 6) or {}).get("detail", {}),
            )

        compare_script = REPO_ROOT / "scripts" / "gen_v2_compare_th232.py"
        if not compare_script.exists():
            return PhaseResult(
                phase=6, name=name, status="failed",
                elapsed_s=t(),
                error=f"compare script not found: {compare_script}",
            )

        # gen_v2_compare_th232 reads <run_dir>/sample/ and <run_dir>/sample_v2/
        # and writes <run_dir>/v2_compare/.
        ctx.layout.compare.mkdir(parents=True, exist_ok=True)
        ctx.logger.info(f"phase 6 ({name}): compare in {ctx.layout.base}")
        try:
            import subprocess
            env = dict(os.environ)
            env.setdefault("PYTHONPATH", str(SCRIPTS_DIR))
            env.setdefault("PYTHONIOENCODING", "utf-8")
            res = subprocess.run(
                [sys.executable, str(compare_script), str(ctx.layout.base)],
                capture_output=True, text=True, env=env, check=False,
                encoding="utf-8", errors="replace",
            )
            if res.returncode != 0:
                ctx.logger.error(
                    f"phase 6: compare exit={res.returncode}\n"
                    f"STDOUT:\n{res.stdout[-800:]}\n"
                    f"STDERR:\n{res.stderr[-800:]}"
                )
                return PhaseResult(
                    phase=6, name=name, status="failed",
                    elapsed_s=t(),
                    error=f"gen_v2_compare exit={res.returncode}",
                    detail={"stderr_tail": res.stderr.splitlines()[-5:]},
                )
        except Exception as e:
            return PhaseResult(
                phase=6, name=name, status="failed",
                elapsed_s=t(), error=f"{type(e).__name__}: {e}",
            )

        out_files = sorted(p.name for p in ctx.layout.compare.iterdir()
                           if p.is_file())
        ctx.logger.info(f"phase 6 ({name}): ok — {out_files}")
        return PhaseResult(
            phase=6, name=name, status="ok",
            elapsed_s=t(), detail={"files": out_files},
        )


def phase_7_index(ctx: RunContext, resume: bool) -> PhaseResult:
    """Render bundle/index.html — minimal navigation.

    Gated by cfg.reports.bundle_index (default False since V2 became
    canonical — task #114). The renderer `_render_index_html` is kept
    intact as a template for future V2-vs-Production compare runs;
    set `cfg.reports.bundle_index = true` to re-enable.
    """
    name = PHASE_NAMES[7]
    with _phase_timer() as t:
        # Gate: bundle landing page is reserved for compare-mode runs.
        # Single-variant V2-canonical runs skip it (user request 2026-06-03).
        if not ctx.cfg.get("reports", {}).get("bundle_index", False):
            return PhaseResult(
                phase=7, name=name, status="skipped",
                elapsed_s=t(),
                detail={
                    "reason": (
                        "bundle_index disabled; enable via "
                        "cfg.reports.bundle_index=true for compare-mode runs"
                    ),
                },
            )

        if resume and _phase_already_done(ctx, 7):
            return PhaseResult(
                phase=7, name=name, status="resumed",
                elapsed_s=t(),
                detail=(_load_prior_phase(ctx, 7) or {}).get("detail", {}),
            )

        try:
            html = _render_index_html(ctx)
            ctx.layout.index_html.write_text(html, encoding="utf-8")
        except Exception as e:
            tb = traceback.format_exc()
            ctx.logger.error(f"phase 7: {type(e).__name__}: {e}\n{tb}")
            return PhaseResult(
                phase=7, name=name, status="failed",
                elapsed_s=t(), error=f"{type(e).__name__}: {e}",
            )
        ctx.logger.info(f"phase 7 ({name}): ok — {ctx.layout.index_html}")
        return PhaseResult(
            phase=7, name=name, status="ok",
            elapsed_s=t(), detail={"path": str(ctx.layout.index_html)},
        )


def phase_8_finalize(ctx: RunContext, resume: bool) -> PhaseResult:
    """Write run_skill_summary.json with phase results + manifest."""
    name = PHASE_NAMES[8]
    with _phase_timer() as t:
        summary = {
            "skill_version": ctx.skill_version,
            "spectrum": str(ctx.spectrum),
            "background": str(ctx.background) if ctx.background else None,
            "bundle_base": str(ctx.layout.base),
            "include_v2": ctx.include_v2,
            "started_at_unix": time.time() - (time.monotonic() - ctx.started_at),
            "elapsed_total_s": round(time.monotonic() - ctx.started_at, 3),
            "phases": ctx.phase_results,
            "metadata": asdict(ctx.metadata),
            "config_resolved": ctx.cfg,
        }
        try:
            ctx.layout.summary_json.write_text(
                json.dumps(summary, indent=2, ensure_ascii=False, default=str),
                encoding="utf-8",
            )
        except OSError as e:
            return PhaseResult(
                phase=8, name=name, status="failed",
                elapsed_s=t(), error=f"write summary: {e}",
            )
        ctx.logger.info(
            f"phase 8 ({name}): ok — {ctx.layout.summary_json.name} written"
        )
        return PhaseResult(
            phase=8, name=name, status="ok",
            elapsed_s=t(), detail={"summary_path": str(ctx.layout.summary_json)},
        )


# ──────────────────────────────────────────────────────────────────
# Helpers: artefact verification, path stringify, index rendering
# ──────────────────────────────────────────────────────────────────

def _stringify_paths(d: Dict[str, Any]) -> Dict[str, Any]:
    """Convert Path-like values nested inside dict/list to str for JSON."""
    def conv(v):
        if isinstance(v, Path):
            return str(v)
        if isinstance(v, dict):
            return {k: conv(x) for k, x in v.items()}
        if isinstance(v, list):
            return [conv(x) for x in v]
        return v
    return conv(d)  # type: ignore[return-value]


def _check_required_artefacts(out: Path, cfg: Dict[str, Any]) -> List[str]:
    """Return list of missing artefact-categories under `out`. Empty list = ok."""
    art = cfg["artefacts"]
    missing: List[str] = []
    if not out.exists():
        return ["<dir not found>"]
    have_files = [p for p in out.iterdir() if p.is_file()]
    have_names = [p.name for p in have_files]

    def _any(suffix: str, *needles: str) -> bool:
        for n in have_names:
            lo = n.lower()
            if lo.endswith(suffix) and (not needles or any(s in lo for s in needles)):
                return True
        return False

    if art.get("json", True) and not _any("_report.json"):
        missing.append("json")
    if art.get("markdown", True) and not _any("_report.md"):
        missing.append("markdown")
    if art.get("html", True) and not _any("_report.html"):
        missing.append("html")
    # F-RPT-03 / v1.18.29 — default-off если ключ отсутствует.
    if art.get("technical_pdf", False) and not _any(".pdf", "technical"):
        missing.append("technical_pdf")
    if art.get("plots", True):
        # PNG лежат в подпапке _plots/
        plots_subdir = next(
            (p for p in out.iterdir() if p.is_dir() and p.name.endswith("_plots")),
            None,
        )
        if plots_subdir is None or not list(plots_subdir.glob("*.png")):
            missing.append("plots")
    # F-RPT-04 / v1.18.29 — default-off если ключ отсутствует.
    if art.get("xml_bq", False) and not _any(".bq.xml"):
        missing.append("xml_bq")
    return missing


def _list_artefacts(out: Path) -> List[str]:
    """Flat relative paths of files under out (depth 2)."""
    if not out.exists():
        return []
    items: List[str] = []
    for p in sorted(out.rglob("*")):
        if p.is_file():
            items.append(str(p.relative_to(out)).replace(os.sep, "/"))
    return items


def _render_index_html(ctx: RunContext) -> str:
    """F-RPT-01 / v1.18.29 — cards-style bundle landing-page.

    Reference: `0_Work/demo_reports/v1_18_25_0_th232/index.html` (gold).
    Structure:
      * header: title + sub
      * kit-block: Sample / Background / Соотношение / Масса / Library
      * cards (3): production sample / V2 sample / V2 vs prod compare
      * footer: F-rules list + дата запуска
    Inline CSS, light/dark via `prefers-color-scheme`. F-386 терминология
    соблюдена — нет слов «ускользание».
    """
    import html as _html

    def _esc(s: Any) -> str:
        return _html.escape(str(s)) if s is not None else ""

    sample_stem = ctx.metadata.stem
    p1 = ctx.phase_results.get(1, {}).get("detail", {})
    p2 = ctx.phase_results.get(2, {}).get("detail", {}) or {}

    # --- Kit-block fields -------------------------------------------------
    sample_real = float(p1.get("real_time_s", 0.0) or 0.0)
    sample_live = float(p1.get("live_time_s", 0.0) or 0.0)
    sample_dt = float(p1.get("dead_time_pct", 0.0) or 0.0)
    sample_total = int(p1.get("channels", 0) or 0)
    # cps estimate если данных хватает; иначе пропускаем
    try:
        sample_cps = (
            (float(p1.get("counts_total", 0)) / sample_live)
            if (sample_live > 0 and p1.get("counts_total"))
            else None
        )
    except Exception:
        sample_cps = None

    bg_name = ctx.background.name if ctx.background else "—"
    sample_name = ctx.spectrum.name

    # Sample mass — приоритет: cfg → metadata → .spe SAMPLEMASS → geometry
    # default. Mirror precedence in `_build_orch_kwargs` (BUG-1 fix). The
    # display-side fall-through keeps the index page consistent with what
    # the pipeline actually consumed.
    mass_kg = (
        ctx.cfg.get("analyze", {}).get("sample_mass_kg")
        or ctx.metadata.sample_mass_kg
    )
    if mass_kg is None:
        spec_pair = _read_spec_mass_kg(ctx.spectrum)
        if spec_pair is not None:
            mass_kg = spec_pair[0]
    if mass_kg is None:
        mass_kg = _default_mass_for_geometry(ctx.metadata.geometry_hint)

    # --- File existence helpers -------------------------------------------
    def _exists(rel: str) -> bool:
        return (ctx.layout.base / rel).exists()

    sample_rel = ctx.layout.sample.relative_to(ctx.layout.base).as_posix()
    v2_rel = ctx.layout.sample_v2.relative_to(ctx.layout.base).as_posix()
    cmp_rel = ctx.layout.compare.relative_to(ctx.layout.base).as_posix()

    # BUG-10 / 2026-06-02 — index links must point at the actual report
    # filenames on disk. The pipeline's report stem is produced by
    # `gamma.reporting.build._safe_filename_stem`, which scrubs embedded
    # S/N tokens (F-115). E.g. input `Th232_420-7-17_Маринелли_0cm.spe`
    # is stored on disk as `Th232_Маринелли_0cm_report.html`, while
    # `ctx.metadata.stem` retains the raw filename stem (S/N intact)
    # because the bundle directory is keyed off the original spectrum
    # name. Previously this function naïvely joined `{ctx.metadata.stem}
    # _report.html` and the resulting <a href="…"> 404'd whenever the
    # source filename carried an S/N token.
    #
    # Resolution order per sub-dir:
    #   1. Phase artefact path (phase 2/4 detail['html']) if present —
    #      this is the authoritative path from the writer.
    #   2. Glob `*_report.html` in the sub-dir — picks up resumed runs
    #      or freshly-rendered pages without phase detail.
    #   3. Fall back to the metadata-stem path (preserves dim-card
    #      behaviour when nothing has been written yet, e.g. unit tests).
    def _phase_html(detail: Dict[str, Any]) -> Optional[str]:
        v = detail.get("html") if isinstance(detail, dict) else None
        return v if isinstance(v, str) and v else None

    def _resolve_report_rel(
        sub_dir: Path,
        sub_rel: str,
        phase_detail: Dict[str, Any],
    ) -> str:
        """Return repo-relative href for the `_report.html` in `sub_dir`.

        Always returns *some* path (even if the file is absent) so
        `available=False` cards still carry a stable href attribute.
        """
        # 1) explicit artefact path from phase result
        p = _phase_html(phase_detail)
        if p:
            try:
                rel = Path(p).resolve().relative_to(ctx.layout.base.resolve())
                return rel.as_posix()
            except (ValueError, OSError):
                pass
        # 2) glob inside the sub-dir
        if sub_dir.exists():
            for candidate in sorted(sub_dir.glob("*_report.html")):
                return f"{sub_rel}/{candidate.name}"
        # 3) metadata-stem fallback (legacy behaviour)
        return f"{sub_rel}/{sample_stem}_report.html"

    p4 = ctx.phase_results.get(4, {}).get("detail", {}) or {}

    sample_report = _resolve_report_rel(ctx.layout.sample, sample_rel, p2)
    v2_report = _resolve_report_rel(ctx.layout.sample_v2, v2_rel, p4)
    compare_report = f"{cmp_rel}/v2_compare_report.html"

    # Effective stem used to discover plots/PDF siblings of the report
    # (these always share the report's filename prefix). We strip the
    # trailing `_report.html` to obtain the writer-side stem.
    def _stem_of(rel: str) -> str:
        leaf = rel.rsplit("/", 1)[-1]
        return leaf[: -len("_report.html")] if leaf.endswith("_report.html") else sample_stem

    sample_disk_stem = _stem_of(sample_report)
    v2_disk_stem = _stem_of(v2_report)

    has_sample = _exists(sample_report)
    has_v2 = ctx.include_v2 and _exists(v2_report)
    has_compare = (
        ctx.include_v2
        and ctx.cfg.get("v2", {}).get("compare", True)
        and _exists(compare_report)
    )

    # --- CSS (verbatim из gold v1.18.25.0 reference, slight extensions) ----
    css = """
:root {
  --bg:#fbfaf3; --bg-secondary:#f5f4ee;
  --text:#1a1a1a; --text-secondary:#5f5e5a;
  --border-secondary:rgba(0,0,0,.20); --radius-md:6px;
  --accent:#37a;
}
@media (prefers-color-scheme:dark) {
  :root { --bg:#1c1c1b; --bg-secondary:#262625;
           --text:#ece9d8; --text-secondary:#b4b2a9;
           --border-secondary:rgba(255,255,255,.22); --accent:#7af; }
}
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
     background:var(--bg);color:var(--text);margin:0;padding:0;line-height:1.5;}
.page{max-width:920px;margin:0 auto;padding:32px 24px;}
h1{font-size:24px;margin:0 0 6px;}
.sub{font-size:14px;color:var(--text-secondary);margin:0 0 28px;}
.card{display:block;border:0.5px solid var(--border-secondary);border-radius:var(--radius-md);
      padding:18px 22px;margin-bottom:14px;background:var(--bg-secondary);
      text-decoration:none;color:inherit;transition:transform .12s ease;}
.card:hover{transform:translateY(-1px);border-color:var(--accent);}
.card-title{font-size:17px;font-weight:600;margin:0 0 6px;color:var(--accent);}
.card-desc{font-size:13px;color:var(--text-secondary);margin:0 0 6px;}
.card-meta{font-size:11px;color:var(--text-secondary);margin-top:8px;display:flex;
           gap:14px;flex-wrap:wrap;}
.card-meta span{display:inline-flex;align-items:center;gap:4px;}
.card.dim{opacity:.55;cursor:not-allowed;}
.tag{display:inline-block;font-size:10px;padding:2px 8px;border-radius:10px;
     background:rgba(55,170,55,.2);color:#3a7;letter-spacing:.5px;
     text-transform:uppercase;font-weight:600;margin-right:6px;}
.tag.exp{background:rgba(140,90,180,.2);color:#85f;}
.tag.bg{background:rgba(170,140,55,.2);color:#a73;}
.kit{border:0.5px solid var(--border-secondary);border-radius:var(--radius-md);
     padding:14px 18px;margin:18px 0 28px;background:var(--bg-secondary);
     font-size:13px;}
.kit b{color:var(--text);}
.kit-row{display:flex;gap:18px;align-items:baseline;flex-wrap:wrap;
         padding:4px 0;border-bottom:0.5px solid var(--border-secondary);}
.kit-row:last-child{border-bottom:none;}
.kit-row .label{color:var(--text-secondary);min-width:120px;font-size:12px;}
code{background:rgba(0,0,0,.06);padding:1px 6px;border-radius:3px;font-size:11px;}
@media (prefers-color-scheme:dark){
  code{background:rgba(255,255,255,.1);}
}
.footer{margin-top:38px;padding-top:16px;border-top:0.5px solid var(--border-secondary);
        font-size:11px;color:var(--text-secondary);}
""".strip()

    title = f"{sample_stem} — комплект отчётов"
    geometry_label = ctx.metadata.geometry_hint or "—"
    detector_label = ctx.metadata.detector_hint or "Gamma-1S NaI 63×63"
    timestamp = _dt.datetime.now().strftime("%Y-%m-%d %H:%M")

    parts: List[str] = []
    parts.append(
        "<!doctype html>\n"
        f"<html lang=\"ru\"><head>\n"
        f"<meta charset=\"utf-8\">\n"
        f"<meta name=\"viewport\" content=\"width=device-width,initial-scale=1.0\">\n"
        f"<title>{_esc(title)} (SpectraVibe {_esc(ctx.skill_version)})</title>\n"
        f"<style>\n{css}\n</style></head>\n"
        "<body><div class=\"page\">\n"
    )
    parts.append(f"<h1>{_esc(title)}</h1>\n")
    parts.append(
        f"<p class=\"sub\">{_esc(geometry_label)} · {_esc(detector_label)} · "
        f"SpectraVibe {_esc(ctx.skill_version)} · bundle "
        f"<code>{_esc(ctx.layout.base.name)}</code> · запуск {timestamp}</p>\n"
    )

    # --- Kit-block (Sample / Background / Соотношение / Масса / Library) ---
    parts.append("<div class=\"kit\">\n")
    sample_meta_bits = []
    if sample_real > 0:
        sample_meta_bits.append(f"{sample_real:.0f} с real")
    if sample_live > 0:
        sample_meta_bits.append(f"{sample_live:.0f} с live")
    if sample_dt > 0:
        sample_meta_bits.append(f"{sample_dt:.2f}% dead")
    if sample_cps is not None:
        sample_meta_bits.append(f"{sample_cps:.2g} cps")
    if sample_total > 0:
        sample_meta_bits.append(f"{sample_total} каналов")
    parts.append(
        "<div class=\"kit-row\">\n"
        "  <span class=\"label\">Sample</span>\n"
        f"  <code>{_esc(sample_name)}</code>\n"
        f"  <span style=\"color:var(--text-secondary)\">{_esc(' · '.join(sample_meta_bits))}</span>\n"
        "</div>\n"
    )
    parts.append(
        "<div class=\"kit-row\">\n"
        "  <span class=\"label\">Background</span>\n"
        f"  <code>{_esc(bg_name)}</code>\n"
        "</div>\n"
    )
    # Соотношение / Масса — известны только если pipeline вычислил;
    # делаем soft-output, иначе пропускаем строку.
    if mass_kg is not None:
        parts.append(
            "<div class=\"kit-row\">\n"
            "  <span class=\"label\">Масса образца</span>\n"
            f"  <b>{mass_kg:.3f} кг</b>\n"
            "  <span style=\"color:var(--text-secondary)\">"
            f"(гипотеза по геометрии {_esc(geometry_label)})</span>\n"
            "</div>\n"
        )
    parts.append(
        "<div class=\"kit-row\">\n"
        "  <span class=\"label\">Библиотека нуклидов</span>\n"
        "  <b>IAEA Live Chart</b>\n"
        "  <span style=\"color:var(--text-secondary)\">F-372 hard-lock: "
        "<code>data/nuclides.json</code></span>\n"
        "</div>\n"
    )
    parts.append("</div>\n")

    # --- Cards -----------------------------------------------------------
    def _card(
        href: str,
        tag_text: str,
        tag_class: str,
        title_text: str,
        desc: str,
        meta_items: List[str],
        available: bool,
    ) -> str:
        if not available:
            return (
                f"<div class=\"card dim\" title=\"артефакт отсутствует — фаза не отработала\">"
                f"<span class=\"tag {tag_class}\">{_esc(tag_text)}</span>"
                f"<div class=\"card-title\">{_esc(title_text)}</div>"
                f"<p class=\"card-desc\">{_esc(desc)} <em>(отсутствует)</em></p>"
                f"</div>\n"
            )
        meta_html = "".join(f"<span>{_esc(m)}</span>" for m in meta_items)
        return (
            f"<a class=\"card\" href=\"{_esc(href)}\">"
            f"<span class=\"tag {tag_class}\">{_esc(tag_text)}</span>"
            f"<div class=\"card-title\">{_esc(title_text)}</div>"
            f"<p class=\"card-desc\">{_esc(desc)}</p>"
            f"<div class=\"card-meta\">{meta_html}</div>"
            f"</a>\n"
        )

    # Card 1: production sample
    sample_meta = [
        "JSON + MD + HTML",
    ]
    if _exists(f"{sample_rel}/{sample_disk_stem}_plots/spectrum.png"):
        sample_meta.append("Spectrum + multiplet PNG")
    if _exists(f"{sample_rel}/{sample_disk_stem}_technical_report.pdf"):
        sample_meta.append("Technical PDF")
    parts.append(_card(
        href=sample_report,
        tag_text="production",
        tag_class="",
        title_text="1. Полный отчёт образца (production)",
        desc=(
            "Step-1..11 pipeline через analyze_lsrm_spe на образце с вычетом фона. "
            "Идентификация по IAEA библиотеке, расчёт активностей в Бк/кг, "
            "мультиплеты, secondary peaks."
        ),
        meta_items=sample_meta,
        available=has_sample,
    ))

    # Card 2: V2 sample (only if include_v2)
    if ctx.include_v2:
        v2_meta = ["JSON + MD + HTML"]
        if _exists(f"{v2_rel}/{v2_disk_stem}_plots/spectrum.png"):
            v2_meta.append("Spectrum + multiplet PNG")
        if _exists(f"{v2_rel}/{v2_disk_stem}_technical_report.pdf"):
            v2_meta.append("Technical PDF")
        parts.append(_card(
            href=v2_report,
            tag_text="experimental v2",
            tag_class="exp",
            title_text="2. Полный отчёт образца (V2-метод)",
            desc=(
                "Тот же analyze_and_report на образце с вычетом фона — "
                "идентичен Пункту 1 по составу артефактов. Отличается "
                "ТОЛЬКО peak search-методом: V2 dual-search "
                "(Mariscotti ∪ matched filter)."
            ),
            meta_items=v2_meta,
            available=has_v2,
        ))

        # Card 3: V2 vs PROD compare
        cmp_meta = ["compare_data.json + HTML", "R&D, не production"]
        parts.append(_card(
            href=compare_report,
            tag_text="experimental F-354",
            tag_class="exp",
            title_text="3. Production vs Experimental v2 — сравнение",
            desc=(
                "KPI-карточки + peaks diff side-by-side между Пунктами 1 и 2. "
                "V2 находит дополнительные пики (matched filter ловит широкие "
                "на continuum)."
            ),
            meta_items=cmp_meta,
            available=has_compare,
        ))

    # --- Footer ------------------------------------------------------------
    parts.append(
        "<div class=\"footer\">\n"
        f"SpectraVibe {_esc(ctx.skill_version)} · "
        "F-RPT-01..05 (v1.18.29) cards landing-page + back-nav + "
        "Technical PDF/BecqMoni opt-in · "
        "F-372 IAEA hard-lock · F-384 demo_reports вне скилла · "
        "F-386 «вылет» терминология · F-389 V2 contract · "
        "F-391 S/N thresholds · F-397 bg-secondary analysis · F-398 run_skill.py · "
        f"запуск {timestamp}\n"
        "</div>\n"
    )

    parts.append("</div></body></html>\n")
    return "".join(parts)


# ──────────────────────────────────────────────────────────────────
# Run-one orchestration
# ──────────────────────────────────────────────────────────────────

PHASE_FNS = [
    phase_0_preflight,
    phase_1_load_calibrate,
    phase_2_prod_analyze,
    phase_3_prod_artefacts,
    phase_4_v2_analyze,
    phase_5_v2_artefacts,
    phase_6_compare,
    phase_7_index,
    phase_8_finalize,
]


def run_one(
    spectrum: Path,
    background: Optional[Path],
    out_base: Path,
    cfg: Dict[str, Any],
    *,
    include_v2: bool,
    resume: bool = False,
    verbose: bool = False,
    quiet: bool = False,
) -> Tuple[int, RunContext]:
    """Run all phases on a single spectrum.

    Returns (exit_code, ctx). Exit codes:
        0 — all ok or all skipped intentionally
        1 — at least one non-fatal phase failed
        2 — fatal: phase 0 failed (input invalid)
    """
    layout = BundleLayout.from_base(out_base, cfg)
    layout.ensure_dirs()
    logger = _setup_logger(layout, verbose=verbose, quiet=quiet)

    skill_version = _read_skill_version()
    metadata = SpectrumMetadata.from_path(spectrum)

    ctx = RunContext(
        spectrum=spectrum,
        background=background,
        metadata=metadata,
        cfg=cfg,
        layout=layout,
        logger=logger,
        skill_version=skill_version,
        include_v2=include_v2,
    )
    logger.info(
        f"=== run_skill start — bundle={layout.base.name}, "
        f"resume={resume}, v={skill_version} ==="
    )

    any_failed = False
    fatal = False
    for phase_id, fn in enumerate(PHASE_FNS):
        try:
            res = fn(ctx, resume)
        except Exception as e:
            tb = traceback.format_exc()
            logger.error(
                f"phase {phase_id}: unhandled exception "
                f"{type(e).__name__}: {e}\n{tb}"
            )
            res = PhaseResult(
                phase=phase_id, name=PHASE_NAMES[phase_id],
                status="failed", elapsed_s=0.0,
                error=f"unhandled: {type(e).__name__}: {e}",
            )
        _record_phase(ctx, phase_id, res)
        if res.status == "failed":
            # Always surface the error reason — otherwise stderr only shows
            # "phase N failed" without why.
            logger.error(
                f"phase {phase_id} ({res.name}) failed: {res.error or '(no detail)'}"
            )
            if phase_id == 0:
                fatal = True
                logger.error("phase 0 fatal — aborting")
                break
            any_failed = True
            # Non-fatal: continue, downstream phases may still produce useful
            # partial output (e.g. compare can fail but index still useful).

    # Always try to finalize summary even on partial failure.
    if fatal:
        # Best-effort: still emit a summary so users can see why we stopped.
        try:
            phase_8_finalize(ctx, resume=False)
        except Exception:
            pass
        return 2, ctx
    if any_failed:
        return 1, ctx
    return 0, ctx


# ──────────────────────────────────────────────────────────────────
# Batch
# ──────────────────────────────────────────────────────────────────

def _batch_row_from_ctx(spectrum: Path, bundle: Path, code: int,
                         include_v2: bool, ctx) -> Dict[str, Any]:
    """Build manifest row from RunContext — kept identical to legacy serial path."""
    return {
        "stem": spectrum.stem,
        "bundle": str(bundle),
        "exit_code": code,
        "include_v2": include_v2,
        "phases_ok": sum(
            1 for r in ctx.phase_results.values()
            if r.get("status") in ("ok", "resumed", "skipped")
        ),
        "phases_failed": sum(
            1 for r in ctx.phase_results.values()
            if r.get("status") == "failed"
        ),
        "elapsed_total_s": round(
            time.monotonic() - ctx.started_at, 2
        ),
    }


def _batch_worker(spectrum_str: str, bundle_str: str, cfg: Dict[str, Any],
                   include_v2: bool, resume: bool, verbose: bool, quiet: bool,
                   ) -> Dict[str, Any]:
    """ProcessPoolExecutor worker — runs run_one and returns serialisable row.

    RunContext is not pickle-safe (carries logging.Logger / BundleLayout),
    so we extract the manifest row inside the child and only return primitives.
    """
    spectrum = Path(spectrum_str)
    bundle = Path(bundle_str)
    code, ctx = run_one(
        spectrum, None, bundle, cfg,
        include_v2=include_v2, resume=resume,
        verbose=verbose, quiet=quiet,
    )
    return _batch_row_from_ctx(spectrum, bundle, code, include_v2, ctx)


def run_batch(
    inputs: Iterable[Path],
    out_root: Path,
    cfg: Dict[str, Any],
    *,
    include_v2: bool,
    resume: bool,
    verbose: bool,
    quiet: bool,
    jobs: int = 1,
) -> int:
    """Run pipeline on multiple spectra, write a top-level manifest.csv.

    AUDIT-F1 (2026-06-25): ``jobs>1`` enables ``ProcessPoolExecutor`` fan-out;
    ``jobs<=1`` keeps the legacy serial path bit-identical for regression gate.
    Output (manifest.csv rows) is always written in input order.
    """
    rows: List[Dict[str, Any]] = []
    worst_code = 0
    inputs = list(inputs)
    n = len(inputs)
    out_root.mkdir(parents=True, exist_ok=True)

    if jobs <= 1:
        print(f"[run_skill] batch: {n} spectra → {out_root}", file=sys.stderr)
        for i, spectrum in enumerate(inputs, start=1):
            # bg=None ⇒ pipeline auto-resolves через background_auto=apply.
            bundle = out_root / spectrum.stem
            print(f"[run_skill] [{i}/{n}] {spectrum.name}", file=sys.stderr)
            code, ctx = run_one(
                spectrum, None, bundle, cfg,
                include_v2=include_v2, resume=resume,
                verbose=verbose, quiet=quiet,
            )
            worst_code = max(worst_code, code)
            rows.append(_batch_row_from_ctx(
                spectrum, bundle, code, include_v2, ctx
            ))
    else:
        max_workers = min(int(jobs), n) if n > 0 else 1
        print(
            f"[run_skill] batch: {n} spectra → {out_root} "
            f"(jobs={max_workers})",
            file=sys.stderr,
        )
        # Submit in input order, collect results into an index→row map so the
        # manifest stays deterministic regardless of completion order.
        results: Dict[int, Dict[str, Any]] = {}
        with _futures.ProcessPoolExecutor(max_workers=max_workers) as pool:
            future_index: Dict[Any, int] = {}
            for i, spectrum in enumerate(inputs):
                bundle = out_root / spectrum.stem
                fut = pool.submit(
                    _batch_worker,
                    str(spectrum), str(bundle), cfg,
                    include_v2, resume, verbose, quiet,
                )
                future_index[fut] = i
            done_count = 0
            for fut in _futures.as_completed(future_index):
                i = future_index[fut]
                spectrum = inputs[i]
                try:
                    row = fut.result()
                except Exception as e:
                    # Surface but do not abort the batch: record a failure
                    # row so the manifest still covers every input.
                    tb = traceback.format_exc()
                    print(
                        f"[run_skill] worker exception on {spectrum.name}: "
                        f"{type(e).__name__}: {e}\n{tb}",
                        file=sys.stderr,
                    )
                    row = {
                        "stem": spectrum.stem,
                        "bundle": str(out_root / spectrum.stem),
                        "exit_code": 2,
                        "include_v2": include_v2,
                        "phases_ok": 0,
                        "phases_failed": 0,
                        "elapsed_total_s": 0.0,
                    }
                results[i] = row
                worst_code = max(worst_code, int(row.get("exit_code", 0)))
                done_count += 1
                print(
                    f"[run_skill] [{done_count}/{n}] done: {spectrum.name} "
                    f"(exit={row.get('exit_code')})",
                    file=sys.stderr,
                )
        # Re-order strictly by input index so manifest.csv is deterministic.
        rows = [results[i] for i in range(n)]

    manifest_path = out_root / "manifest.csv"
    _write_csv(manifest_path, rows)
    print(f"[run_skill] manifest → {manifest_path}", file=sys.stderr)
    return worst_code


def _write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    import csv
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


# ──────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="run_skill",
        description=(
            "End-to-end SpectraVibe orchestrator. "
            "Один прогон: спектр → полный bundle (PROD + опц V2 + compare + index)."
        ),
    )
    p.add_argument(
        "spectrum", nargs="?",
        help="Путь к .spe/.chn/.n42 спектру. Не требуется при --batch или --resume.",
    )
    p.add_argument(
        "--background", default=None,
        help="Путь к фоновому спектру (.spe). Если не задан, "
             "auto-detect via background_auto=apply.",
    )
    p.add_argument(
        "--auto-detect-bg", action="store_true",
        help="Дополнительно искать sibling-файл по эвристике filename "
             "(фон/bkg/закр) и передать как --background, если найден.",
    )
    p.add_argument(
        "--mass", "--sample-mass-kg", dest="mass", type=float, default=None,
        help="Масса образца в кг. Если не задано — выводится из имени файла "
             "(токены kg/г) или из геометрии (default 1.0 кг для unknown).",
    )
    p.add_argument(
        "--output-dir", "-o", default=None,
        help="Куда положить bundle. Default: $GAMMA_DEMO_REPORTS_DIR/<spectrum-stem>/",
    )
    p.add_argument(
        "--include-v2", action="store_true",
        help="Дополнительно прогнать V2 ветку + compare.",
    )
    p.add_argument(
        "--config", default=None, type=Path,
        help="YAML/JSON config файл. Override через deep-merge поверх дефолтов.",
    )
    p.add_argument(
        "--resume", default=None,
        help="Возобновить прерванный прогон по пути к bundle-директории. "
             "Завершённые фазы (.phases/phase_N.done) пропускаются.",
    )
    p.add_argument(
        "--batch", default=None,
        help="Glob или директория с .spe файлами для батч-режима. "
             "Каждый файл получает свой bundle под <output-dir>/<stem>/.",
    )
    p.add_argument(
        "--jobs", default=1, type=int,
        help="AUDIT-F1: процессов в батче (default 1 — serial, bit-identical "
             "к легаси-пути). >1 — concurrent.futures.ProcessPoolExecutor; "
             "manifest.csv остаётся в порядке входов.",
    )
    p.add_argument(
        "--no-pdf", action="store_true", help="Skip technical PDF.",
    )
    p.add_argument(
        "--no-html", action="store_true", help="Skip HTML report.",
    )
    p.add_argument(
        "--no-plots", action="store_true", help="Skip PNG plots.",
    )
    p.add_argument(
        "--no-markdown", action="store_true", help="Skip Markdown.",
    )
    p.add_argument(
        "--no-xml", action="store_true", help="Skip BecqMoni XML round-trip.",
    )
    p.add_argument(
        "--allow-stage3", action="store_true",
        help=(
            "F-442: разрешить Stage 3 identification (EXOTIC nuclide pool: "
            "Ga-67, Tc-99m, In-111, Mn-54, Na-22, Be-7, Am-241, Eu-152, "
            "Ba-133, Co-57, Ti-44, Sc-44). По умолчанию OFF — на природных "
            "samples Stage 3 даёт medical false-positives (Ga-67 vs Ra-226+U-235 "
            "186 keV). Включать ТОЛЬКО для low-background castle calibration "
            "с sealed AmTiCsEu Marinelli-source. Эквивалент: GAMMA_ALLOW_STAGE3=1."
        ),
    )
    p.add_argument("--quiet", "-q", action="store_true",
                   help="Suppress console output (file log unaffected).")
    p.add_argument("--verbose", "-v", action="store_true",
                   help="Verbose console output (DEBUG level).")
    p.add_argument("--version", action="store_true",
                   help="Show skill version and exit.")
    return p


def _resolve_out_base(args, stem: str) -> Path:
    """Implement F-384 + HARD-LOCK v1.2.6 — timestamped auto-default output dir."""
    if args.output_dir:
        return Path(args.output_dir).resolve()
    # HARD-LOCK v1.2.6: auto-default folder must carry timestamp prefix
    ts = _dt.datetime.now().strftime("%Y-%m-%d_%H-%M")
    ts_stem = f"{ts}_{stem}"
    env = os.environ.get("GAMMA_DEMO_REPORTS_DIR")
    if env:
        return Path(env).resolve() / ts_stem
    # Defer to demo_reports_root helper (creates default at first run).
    try:
        from gamma.data.demo_reports_root import ensure_demo_reports_root
        root = ensure_demo_reports_root(interactive=False)
    except Exception:
        root = REPO_ROOT / "demo_reports"
        root.mkdir(parents=True, exist_ok=True)
    return root / ts_stem


def _apply_cli_overrides(cfg: Dict[str, Any], args) -> None:
    """Reflect CLI flags into resolved config."""
    if args.no_pdf:
        cfg["artefacts"]["technical_pdf"] = False
    if args.no_html:
        cfg["artefacts"]["html"] = False
    if args.no_plots:
        cfg["artefacts"]["plots"] = False
    if args.no_markdown:
        cfg["artefacts"]["markdown"] = False
    if args.no_xml:
        cfg["artefacts"]["xml_bq"] = False
    if args.mass is not None:
        cfg["analyze"]["sample_mass_kg"] = float(args.mass)
    if args.include_v2:
        cfg["v2"]["enabled"] = True
    # F-442 / v1.30.2 — Stage 3 opt-in (CLI flag OR env var).
    if getattr(args, "allow_stage3", False):
        cfg["analyze"]["allow_stage3"] = True
    env_s3 = os.environ.get("GAMMA_ALLOW_STAGE3", "").strip().lower()
    if env_s3 in ("1", "true", "yes", "on"):
        cfg["analyze"]["allow_stage3"] = True
    elif env_s3 in ("0", "false", "no", "off"):
        cfg["analyze"]["allow_stage3"] = False


def _resolve_batch_inputs(pattern: str) -> List[Path]:
    p = Path(pattern)
    if p.is_dir():
        return sorted(p.glob("*.spe"))
    return sorted(Path(x) for x in _glob.glob(pattern, recursive=True))


def main(argv: Optional[List[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.version:
        print(_read_skill_version())
        return 0

    # Load config (or default-only)
    try:
        cfg = _load_config(args.config)
    except Exception as e:
        print(f"ERROR: config load: {e}", file=sys.stderr)
        return 2
    _apply_cli_overrides(cfg, args)
    include_v2 = bool(cfg["v2"].get("enabled") or args.include_v2)

    # Resume path
    if args.resume:
        bundle = Path(args.resume).resolve()
        if not bundle.exists():
            print(f"ERROR: --resume bundle not found: {bundle}", file=sys.stderr)
            return 2
        # Recover input paths from prior phase 0 marker if available.
        marker0 = bundle / ".phases" / "phase_0.done"
        if not marker0.exists():
            print(
                f"ERROR: cannot resume — {marker0} missing. "
                f"Run без --resume для нового прогона.",
                file=sys.stderr,
            )
            return 2
        try:
            prior = json.loads(marker0.read_text(encoding="utf-8"))
            detail = prior.get("detail", {})
            spectrum = Path(detail["spectrum"])
            bg_str = detail.get("background")
            background = Path(bg_str) if bg_str else None
        except (KeyError, json.JSONDecodeError) as e:
            print(f"ERROR: resume marker malformed: {e}", file=sys.stderr)
            return 2
        code, _ = run_one(
            spectrum, background, bundle, cfg,
            include_v2=include_v2, resume=True,
            verbose=args.verbose, quiet=args.quiet,
        )
        return code

    # Batch path
    if args.batch:
        inputs = _resolve_batch_inputs(args.batch)
        if not inputs:
            print(f"ERROR: --batch pattern matched 0 files: {args.batch}",
                  file=sys.stderr)
            return 2
        out_root = (Path(args.output_dir).resolve() if args.output_dir
                    else _resolve_out_base(args, "batch"))
        out_root.mkdir(parents=True, exist_ok=True)
        return run_batch(
            inputs, out_root, cfg,
            include_v2=include_v2, resume=False,
            verbose=args.verbose, quiet=args.quiet,
            jobs=int(getattr(args, "jobs", 1) or 1),
        )

    # Single path
    if not args.spectrum:
        parser.print_usage(sys.stderr)
        print(
            "ERROR: spectrum argument required (or use --batch / --resume).",
            file=sys.stderr,
        )
        return 2
    spectrum = Path(args.spectrum).resolve()
    background = Path(args.background).resolve() if args.background else None
    if background is None and args.auto_detect_bg:
        # Use minimal logger for pre-bundle detection.
        tmp_log = logging.getLogger("run_skill.bootstrap")
        tmp_log.addHandler(logging.NullHandler())
        background = _auto_detect_background(spectrum, tmp_log)

    out_base = _resolve_out_base(args, spectrum.stem)
    code, _ = run_one(
        spectrum, background, out_base, cfg,
        include_v2=include_v2, resume=False,
        verbose=args.verbose, quiet=args.quiet,
    )
    return code


if __name__ == "__main__":
    sys.exit(main())
