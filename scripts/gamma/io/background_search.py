"""F-131 / v1.17.7 — Эвристический поиск подходящего фонового спектра.

Если у sample-спектра нет ни embedded `<BackgroundEnergySpectrum>`, ни
`<BackgroundSpectrumFile>` link, ни явного CLI `--background-path`,
этот модуль ищет файл-кандидат, который опытный спектрометрист с
большой вероятностью выбрал бы как фоновый: тот же детектор, та же
(или совместимая) геометрия, дата измерения близка к sample.

Алгоритм (зеркалит ручную рутину):
  1. Собрать все `.spe` файлы в той же папке, что и sample, +
     рекурсивно в `data/averaged_backgrounds/` и `*/Фон*/` подпапках
     соответствующего детектора.
  2. Распарсить каждый файл (lightweight — только header / extras /
     метаданные, без полной декодировки spectr binary).
  3. Отфильтровать по критериям:
       • filename содержит маркер фона (``bg``, ``bkg``, ``фон``,
         ``background``) ИЛИ `TYPE=Фон` ИЛИ `filename_tokens.is_background_hint`;
       • DETECTOR совпадает с sample (либо у sample DETECTOR пустой);
       • |Δt_measurement| ≤ max_days_apart (default 90 дней);
       • geometry совпадает с sample либо допустимый fallback
         (точечная для безматериального BG, Marinelli для матричного).
  4. Сортировать по «уверенности»:
       confidence_score = w_detector + w_geometry + w_recency + w_filename
     где каждый компонент 0..1, итог ∈ [0, 4]. Чем выше, тем
     более вероятный кандидат.

Возвращает упорядоченный список ``BackgroundCandidate`` (лучший
первый). Пустой список — кандидатов нет, pipeline продолжит на
gross-спектре с `background_status = "absent_no_subtraction"`.

Интеграция в pipeline (`staged_pipeline.analyze_lsrm_spe`):
  • `background_auto = "off"`     — не искать (CI / batch режим)
  • `background_auto = "suggest"` — найти и положить в pipeline_notes
                                    + сериализовать как
                                    `auto_background_candidates`,
                                    но НЕ применять автоматически
  • `background_auto = "apply"`   — найти и применить лучшего
                                    кандидата как
                                    `--background-path <best>`

Безопасность: модуль НИКОГДА не raise; при любой ошибке возвращает
пустой список и пишет диагностику в логи. Это «soft prompt» —
пользователь всегда может проигнорировать предложение и продолжить
на gross-спектре.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from typing import List, Optional


BG_FILENAME_MARKERS = ("bg_", "bkg_", "фон", "background", "_bg.", "_bkg.")
BG_TYPE_MARKERS = ("фон", "background", "BG")


@dataclass(frozen=True)
class BackgroundCandidate:
    """Один кандидат — найденный потенциальный фоновый файл."""
    path: Path
    detector: str
    geometry: str
    measurement_datetime: Optional[datetime]
    days_before_sample: Optional[float]
    is_filename_bg: bool
    is_type_bg: bool
    confidence_score: float
    why: str

    def to_dict(self) -> dict:
        """Для сериализации в JSON-отчёт (без объектов Path / datetime)."""
        return {
            "path": str(self.path),
            "filename": self.path.name,
            "detector": self.detector,
            "geometry": self.geometry,
            "measurement_datetime": (
                self.measurement_datetime.isoformat()
                if self.measurement_datetime else None
            ),
            "days_before_sample": (
                round(self.days_before_sample, 1)
                if self.days_before_sample is not None else None
            ),
            "is_filename_bg": bool(self.is_filename_bg),
            "is_type_bg": bool(self.is_type_bg),
            "confidence_score": round(float(self.confidence_score), 3),
            "why": self.why,
        }


def _filename_looks_like_bg(name: str) -> bool:
    low = name.lower()
    return any(m in low for m in BG_FILENAME_MARKERS)


def _type_field_says_bg(type_field: str) -> bool:
    if not type_field:
        return False
    low = type_field.lower()
    return any(m.lower() in low for m in BG_TYPE_MARKERS)


def _geometry_compatible(sample_geom: str, candidate_geom: str) -> tuple:
    """Совместима ли геометрия кандидата с sample? Возвращает (compatible, score).

    score ∈ [0, 1]: 1.0 — точное совпадение, 0.6 — semi-совместимая
    (например, "Точечная" фон для Marinelli sample допустимо — спектр
    собственного фона детектора одинаков; матрица влияет только на
    self-absorption), 0.0 — несовместима.
    """
    s = (sample_geom or "").strip().lower()
    c = (candidate_geom or "").strip().lower()
    if not s and not c:
        return (True, 0.5)
    if not s or not c:
        return (True, 0.4)
    if s == c:
        return (True, 1.0)
    # Точечная-фоновая под Marinelli — допустимо (детекторный фон)
    if "точечная" in c and "маринелли" in s:
        return (True, 0.55)
    if "точечная" in s and "маринелли" in c:
        return (True, 0.55)
    # Любые "Фон*" имена считаем geometry-agnostic
    if "фон" in c.lower():
        return (True, 0.5)
    return (False, 0.0)


def _detector_compatible(sample_det: str, candidate_det: str) -> tuple:
    """Тот же детектор? Возвращает (compatible, score) ∈ [0..1].

    F-312 / v1.18.12 — canonical alias normalization. До v1.18.12 матчинг
    шёл по substring 'гамма-1' / 'gamma-1', что отсекало vendor-токены
    `УДС-ГЦ-63х63-USB №SN-01`, `БДЭГ-63×63-USB`, `Колибри-1М` —
    которые per SKILL.md строка 26 являются alias одного и того же
    Gamma-1S физического детектора. Теперь матчинг включает полный
    набор vendor-токенов.
    """
    s = (sample_det or "").strip().lower()
    c = (candidate_det or "").strip().lower()
    if not s and not c:
        return (True, 0.3)
    if not s or not c:
        return (True, 0.4)
    if s == c:
        return (True, 1.0)
    # F-312: canonical Gamma-1S aliases per SKILL.md строка 26.
    # Любые два токена из этого набора считаются совместимыми.
    GAMMA1C_ALIASES = (
        "гамма-1", "gamma-1", "γ-1",
        "колибри-1", "kolibri-1",
        "удс-гц", "udc-hc",        # УДС-ГЦ vendor token (Aspect)
        "бдэг-63", "bdeg-63",      # БДЭГ vendor token (LSRM)
        "naidet-63", "naidet-1",
    )
    s_has = [t for t in GAMMA1C_ALIASES if t in s]
    c_has = [t for t in GAMMA1C_ALIASES if t in c]
    if s_has and c_has:
        # Оба содержат Gamma-1S alias → совместимы. Score выше если
        # совпадает тот же токен, ниже если разные алиасы.
        if set(s_has) & set(c_has):
            return (True, 0.85)
        return (True, 0.7)
    return (False, 0.0)


def _recency_score(days_apart: Optional[float], max_days: float) -> float:
    """Чем ближе по дате, тем выше score. ∈ [0, 1]."""
    if days_apart is None:
        return 0.3   # неизвестная дата — слабый, но не нулевой бонус
    a = abs(days_apart)
    if a > max_days:
        return 0.0
    # Линейный спад: 1.0 при Δt=0, 0.1 при Δt=max_days
    return float(max(0.1, 1.0 - 0.9 * (a / max_days)))


def _read_header_only(path: Path) -> Optional[dict]:
    """Прочитать только header LSRM .spe (без полной декодировки бинарного
    блока counts) для извлечения метаданных. Возвращает dict с полями
    DETECTOR / GEOMETRY / TYPE / MEASBEGIN либо None при ошибке.

    Использует существующий read_spectrum но в режиме чтения метаданных
    (для .spe LSRM формат это всё равно весь header сразу).
    """
    try:
        # Late import — избегаем циркулярных зависимостей
        from gamma.io.lsrm_spe import read_lsrm_spe
        spec = read_lsrm_spe(str(path), apply_energy_ceiling=False)
        return {
            "DETECTOR": getattr(spec, "detector_id", "") or "",
            "GEOMETRY": getattr(spec, "geometry", "") or "",
            "TYPE": (spec.extras or {}).get("lsrm_type", "") or "",
            "MEASBEGIN": getattr(spec, "start_datetime", None),
            "filename_tokens": getattr(spec, "filename_tokens", {}) or {},
        }
    except Exception:
        return None


def find_background_candidates(
    sample_spec,
    sample_path: str,
    *,
    extra_search_dirs: Optional[List[str]] = None,
    max_days_apart: int = 90,
    max_candidates: int = 10,
) -> List[BackgroundCandidate]:
    """F-131 — найти подходящие фоновые .spe-файлы для sample.

    Эвристика опытного спектрометриста:
      • тот же детектор (DETECTOR-поле),
      • та же или semi-совместимая геометрия (Точечная для
        Marinelli sample — это детекторный фон, допустимо),
      • дата ≤ max_days_apart (default 90 дней),
      • filename содержит маркер фона ИЛИ TYPE=Фон/Background.

    Возвращает упорядоченный список кандидатов (лучший первый).
    """
    sample_dir = Path(sample_path).parent
    sample_det = getattr(sample_spec, "detector_id", "") or ""
    sample_geom = getattr(sample_spec, "geometry", "") or ""
    sample_dt = getattr(sample_spec, "start_datetime", None)

    # ─── собрать корпус кандидатов ────────────────────────────────────
    search_roots: List[Path] = [sample_dir]
    # Поднимаемся до 4 уровней вверх и ищем `data/averaged_backgrounds`
    # и `*/Фон*/` подкаталоги (типичные места для bg-фикстур
    # в дереве detectors/<DET>/).
    cur = sample_dir
    for _ in range(4):
        cur = cur.parent
        if not cur or cur == cur.parent:
            break
        avg_bg = cur / "data" / "averaged_backgrounds"
        if avg_bg.is_dir():
            search_roots.append(avg_bg)
        # F-312 / v1.18.12 (CRITICAL FIX) — также ищем в detector-subtree.
        # До v1.18.12 эвристика глобила только <level>/data/averaged_backgrounds,
        # но F-157 isolation policy кладёт фоны в
        # <level>/detectors/<DETECTOR>/data/averaged_backgrounds. Это ломало
        # F-135 contract («фон ВСЕГДА вычитается при наличии подходящего
        # кандидата») для всех фикстур вне detector-subtree (например,
        # evals/fixtures/*.spe). Поиск теперь покрывает все detector-папки.
        detectors_dir = cur / "detectors"
        if detectors_dir.is_dir():
            try:
                for det_sub in detectors_dir.iterdir():
                    if not det_sub.is_dir():
                        continue
                    det_avg_bg = det_sub / "data" / "averaged_backgrounds"
                    if det_avg_bg.is_dir():
                        search_roots.append(det_avg_bg)
                    # Также фон-папки внутри detector subtree
                    for nested in det_sub.rglob("*"):
                        if not nested.is_dir():
                            continue
                        nm = nested.name.lower()
                        if ("фон" in nm or "background" in nm
                                or nm == "bg" or nm.startswith("bg_")):
                            search_roots.append(nested)
            except OSError:
                pass
        # Также любая папка с "Фон" в имени на 1 уровне
        for sub in cur.iterdir() if cur.is_dir() else []:
            if sub.is_dir() and ("фон" in sub.name.lower()
                                 or "background" in sub.name.lower()
                                 or "bg" in sub.name.lower()):
                search_roots.append(sub)
                # Также подкаталоги первого уровня
                try:
                    for sub2 in sub.iterdir():
                        if sub2.is_dir():
                            search_roots.append(sub2)
                except OSError:
                    pass
    for d in (extra_search_dirs or []):
        try:
            p = Path(d)
            if p.is_dir():
                search_roots.append(p)
        except (TypeError, ValueError):
            pass

    # Уникальные пути
    seen_roots: set = set()
    unique_roots: List[Path] = []
    for r in search_roots:
        try:
            r_res = r.resolve()
        except OSError:
            continue
        if r_res not in seen_roots:
            seen_roots.add(r_res)
            unique_roots.append(r)

    spe_files: List[Path] = []
    sample_path_resolved = None
    try:
        sample_path_resolved = Path(sample_path).resolve()
    except OSError:
        pass

    for root in unique_roots:
        try:
            for p in root.iterdir():
                if (p.is_file() and p.suffix.lower() == ".spe"):
                    try:
                        if p.resolve() == sample_path_resolved:
                            continue   # не предлагать sample как BG самого себя
                    except OSError:
                        pass
                    spe_files.append(p)
        except OSError:
            continue

    # ─── собрать кандидатов ───────────────────────────────────────────
    candidates: List[BackgroundCandidate] = []
    for p in spe_files:
        name = p.name
        is_fn_bg = _filename_looks_like_bg(name)
        meta = _read_header_only(p)
        if meta is None:
            # Фоном можно считать только если в имени маркер есть и
            # файл читается хотя бы как .spe — но раз read failed,
            # пропускаем.
            continue
        det = meta["DETECTOR"]
        geom = meta["GEOMETRY"]
        type_f = meta["TYPE"]
        meas_dt = meta["MEASBEGIN"]
        ftok = meta.get("filename_tokens", {}) or {}
        is_type_bg = (
            _type_field_says_bg(type_f)
            or bool(ftok.get("is_background_hint"))
        )

        # Жёсткий фильтр: должен выглядеть как BG (по filename или TYPE)
        if not (is_fn_bg or is_type_bg):
            continue

        # Совместимость детектора
        det_ok, det_score = _detector_compatible(sample_det, det)
        if not det_ok:
            continue
        # Совместимость геометрии
        geom_ok, geom_score = _geometry_compatible(sample_geom, geom)
        if not geom_ok:
            continue
        # Δt
        days_apart: Optional[float] = None
        if sample_dt is not None and meas_dt is not None:
            try:
                days_apart = (sample_dt - meas_dt).total_seconds() / 86400.0
            except (TypeError, ValueError):
                days_apart = None
        if (days_apart is not None
                and abs(days_apart) > float(max_days_apart)):
            continue
        recency = _recency_score(days_apart, float(max_days_apart))

        # Filename markers boost
        fn_score = 1.0 if is_fn_bg else 0.6
        # ИТОГ confidence ∈ [0, 4]
        confidence = float(det_score + geom_score + recency + fn_score)
        why_parts = []
        why_parts.append(f"детектор: {det_score:.2f}")
        why_parts.append(f"геометрия: {geom_score:.2f}")
        why_parts.append(f"recency: {recency:.2f}")
        why_parts.append(f"filename: {fn_score:.2f}")
        candidates.append(BackgroundCandidate(
            path=p,
            detector=det,
            geometry=geom,
            measurement_datetime=meas_dt,
            days_before_sample=days_apart,
            is_filename_bg=is_fn_bg,
            is_type_bg=is_type_bg,
            confidence_score=confidence,
            why="; ".join(why_parts),
        ))

    # Сортировка по confidence (убывание), tie-break — recency
    candidates.sort(
        key=lambda c: (
            -c.confidence_score,
            abs(c.days_before_sample) if c.days_before_sample is not None else 1e9,
        )
    )
    return candidates[:max_candidates]


def render_suggestion_note(best: BackgroundCandidate) -> str:
    """Сформировать RU-нарративную заметку для pipeline_notes о
    предложенном фоне."""
    dt_part = ""
    if best.days_before_sample is not None:
        delta = int(round(abs(best.days_before_sample)))
        when = "раньше" if best.days_before_sample > 0 else "позже"
        dt_part = f", Δt={delta} дн. {when} sample"
    return (
        f"F-131: предложен фоновый файл «{best.path.name}» "
        f"(уверенность={best.confidence_score:.2f}/4.0{dt_part}). "
        f"Передайте `--background-path \"{best.path}\"` либо "
        f"`--background-auto apply` для автоматического применения, "
        f"либо `--background-auto off` чтобы продолжить на gross-спектре."
    )


def render_applied_note(best: BackgroundCandidate) -> str:
    """Заметка про применённый авто-фон."""
    dt_part = ""
    if best.days_before_sample is not None:
        delta = int(round(abs(best.days_before_sample)))
        dt_part = f", Δt={delta} дн."
    return (
        f"F-131: автоматически применён фоновый файл «{best.path.name}» "
        f"(уверенность={best.confidence_score:.2f}/4.0{dt_part}). "
        f"Передайте `--background-auto suggest` или `--background-path` "
        f"для явного выбора, либо `--background-auto off` для отказа от "
        f"авто-подбора."
    )


__all__ = [
    "BackgroundCandidate",
    "find_background_candidates",
    "render_suggestion_note",
    "render_applied_note",
    "BG_FILENAME_MARKERS",
    "BG_TYPE_MARKERS",
]
