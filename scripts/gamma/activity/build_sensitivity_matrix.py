# -*- coding: utf-8 -*-
"""F-448 — построение матрицы чувствительности МЭС (LSRM §12) из эталонов.

Модуль связывает уже существующий fitter/builder ``template_method.py``
(``build_R_matrix`` / ``template_assay``) с реальными эталонными спектрами
Gamma-1S (Маринелли, Поверка 2024). Здесь:

1. Таблица энергетических окон ``MES_WINDOWS`` — характеристические линии
   Th-232 chain / Cs-137 / K-40 / Ra-226 chain с провенансом по
   ``data/nuclides.json`` (каждое окно подписано линией и её энергией).
2. ``load_etalon_calibration_specs`` — читает .spe эталонов через
   существующий ``gamma.io.lsrm_spe.read_lsrm_spe`` (НЕ дублирует reader),
   достаёт сертифицированную активность из COMMENT через
   ``parse_lsrm_passport_comment`` (НЕ дублирует COMMENT-парсер), строит
   ``CalibrationSpec`` на нуклид.
3. ``build_sensitivity_matrix_from_etalons`` — собирает ``BackgroundSpec``
   из фона и вызывает ``build_R_matrix`` → ``SensitivityMatrix``.

ВАЖНО (anti-hallucination): каждый эталон имеет СВОЮ энергетическую
калибровку (gain MCA слегка различается между файлами — проверено на
Поверке 2024: sample cal a0=-8.29, Th-эталон 420-17031 a0=-12.37, фон
a0=-5.04). Поэтому ``sum_in_windows`` для каждого спектра вызывается с
ЕГО собственной осью энергий (``spec.channel_to_energy``), а не с общей
сеткой. Окна заданы в keV, что делает суммирование инвариантным к
небольшим различиям gain между файлами.

Циркулярность: эталон Th-232 — это файл ``420-17031`` (A=860 Бк/кг,
независимый источник), НЕ образец ``420-7-17`` (A=1940 Бк/кг). Раннер
``run_template_assay.py`` ассертит, что путь образца != путь Th-эталона.

Reference: LSRM Algorithmic Foundations 2022 §12 «Шаблонный метод».
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from gamma.io.lsrm_spe import read_lsrm_spe
from gamma.activity.template_method import (
    CalibrationSpec,
    BackgroundSpec,
    SensitivityMatrix,
    build_R_matrix,
)


# ---------------------------------------------------------------------------
# Окна МЭС с провенансом (характеристические линии, data/nuclides.json)
# ---------------------------------------------------------------------------
#
# Формат: (E_low_keV, E_high_keV, "линия-провенанс").
# Окна шире FWHM(E) NaI (~7% на 662 keV → ~46 keV полная), чтобы вместить
# полный photopeak с запасом на низкоэнергетический хвост, но достаточно
# узкие, чтобы изолировать характеристическую линию нуклида. Перекрытия
# Tl-208 583 / Bi-214 609 разнесены границей 600 keV; cross-talk между
# окнами штатно учитывается МЭС через off-diagonal члены матрицы R
# (каждый эталон вносит вклад во ВСЕ окна, LSQ решает совместно).
#
# Th-232 chain (один эталон 420-17031, цепочка в равновесии → одна A_k):
#   238.63 keV  Pb-212  (I=43.6%)   data/nuclides.json:Pb-212
#   338.32 keV  Ac-228  (I=11.27%)  data/nuclides.json:Ac-228
#   583.19 keV  Tl-208  (I=30.55%)  data/nuclides.json:Tl-208
#   911.20 keV  Ac-228  (I=25.8%)   data/nuclides.json:Ac-228 (+968.97 I=15.8%)
#   2614.51 keV Tl-208  (I=35.85%)  data/nuclides.json:Tl-208
# Cs-137 (эталон 420-7-14):
#   661.66 keV  Cs-137  (I=85.1%)   data/nuclides.json:Cs-137
# K-40 (эталон 420-7-20):
#   1460.82 keV K-40    (I=10.66%)  data/nuclides.json:K-40
# Ra-226 chain (эталон 420-7-18, цепочка в равновесии → одна A_k):
#   351.93 keV  Pb-214  (I=35.72%)  data/nuclides.json:Pb-214
#   609.32 keV  Bi-214  (I=45.44%)  data/nuclides.json:Bi-214
#   1764.49 keV Bi-214  (I=15.29%)  data/nuclides.json:Bi-214

MES_WINDOWS: Tuple[Tuple[float, float, str], ...] = (
    (228.0, 252.0, "Pb-212 238.63 (Th-232 chain)"),
    (322.0, 358.0, "Ac-228 338.32 (Th-232 chain) | Pb-214 351.93 (Ra-226 chain)"),
    (560.0, 600.0, "Tl-208 583.19 (Th-232 chain)"),
    (600.0, 632.0, "Bi-214 609.32 (Ra-226 chain)"),
    (645.0, 685.0, "Cs-137 661.66"),
    (885.0, 990.0, "Ac-228 911.20 + 968.97 (Th-232 chain)"),
    (1380.0, 1540.0, "K-40 1460.82"),
    (1700.0, 1830.0, "Bi-214 1764.49 (Ra-226 chain)"),
    (2540.0, 2690.0, "Tl-208 2614.51 (Th-232 chain)"),
)


def windows_keV_only(
    windows: Sequence[Tuple[float, float, str]] = MES_WINDOWS,
) -> Tuple[Tuple[float, float], ...]:
    """Strip the provenance label → (lo, hi) tuples for build_R_matrix."""
    return tuple((float(lo), float(hi)) for lo, hi, _ in windows)


# ---------------------------------------------------------------------------
# Реестр эталонов по нуклиду (Поверка 2024, Маринелли 0cm)
# ---------------------------------------------------------------------------
# Относительные пути от корня detectors/Gamma-1S/raw_lsrm/Work/...
# Ключ — канонический нуклид (для chain — родитель/представитель цепочки).

_MARINELLI_DIR = (
    "detectors/Gamma-1S/raw_lsrm/Work/BG/Gamma-1S/"
    "Spe - поверки/Поверка 2024/Маринелли"
)

ETALON_FILES: Dict[str, str] = {
    "Th-232": f"{_MARINELLI_DIR}/Th-232_420-17031_Маринелли_0cm.spe",
    "Cs-137": f"{_MARINELLI_DIR}/Cs137_420-7-14_Маринелли_0cm.spe",
    "K-40":   f"{_MARINELLI_DIR}/K40_420-7-20_Маринелли_0cm.spe",
    "Ra-226": f"{_MARINELLI_DIR}/Ra226_420-7-18_Маринелли_0cm.spe",
}

# Файл-образец (НЕ эталон) — для assert циркулярности в раннере.
SAMPLE_TH232_FILE = (
    f"{_MARINELLI_DIR}/Th232_420-7-17_Маринелли_0cm.spe"
)

DEFAULT_BG_FILE = (
    "detectors/Gamma-1S/raw_lsrm/Work/BG/Gamma-1S/"
    "Spe - поверки/Поверка 2024/Фон закр кр/Фон закр кр вода_13.spe"
)


@dataclass
class EtalonInfo:
    """Что прочитано из одного эталонного .spe (для аудита/провенанса)."""
    nuclide: str
    path: str
    t_live_s: float
    A_certified_Bq_per_kg: float
    sigma_A_rel: float
    cert_raw: str
    energy_cal: tuple
    matrix_name: str = ""          # F-448: имя матрицы наполнения (ОИСН-06 и т.п.)
    density_g_cm3: Optional[float] = None  # F-448: плотность Маринелли


def matrix_descriptor(spec) -> Tuple[str, Optional[float]]:
    """F-448: извлечь (имя_матрицы, плотность г/см3) из spec.extras.

    МЭС требует ОДНУ И ТУ ЖЕ матрицу/геометрию у эталона и пробы (LSRM
    §12.1): счёт детектора пропорционален ОБЪЁМНОЙ активности
    (Бк/л = Бк/кг · ρ) с поправкой на самопоглощение, которая зависит от
    плотности и состава. Нормировка R на УДЕЛЬНУЮ активность (Бк/кг)
    корректна только когда ρ и состав совпадают. Дескриптор используется
    guard'ом для предупреждения о несовпадении.

    Источники (anti-hallucination):
        spec.extras['lsrm_sample_density_g_cm3'] — плотность (reader),
        spec.extras['lsrm_material'] — JSON {"Name","Ro","Compound"}.
    Возвращает ("", None) если данных нет (не угадываем).
    """
    extras = getattr(spec, "extras", {}) or {}
    rho = extras.get("lsrm_sample_density_g_cm3")
    try:
        rho = float(rho) if rho is not None else None
    except (TypeError, ValueError):
        rho = None
    name = ""
    mat = extras.get("lsrm_material")
    if isinstance(mat, str) and mat.strip():
        try:
            md = json.loads(mat)
            name = str(md.get("Name", "") or "")
            if rho is None and md.get("Ro") is not None:
                rho = float(md["Ro"])
        except (json.JSONDecodeError, TypeError, ValueError):
            name = ""
    return name, rho


def _passport_for_nuclide(spec, nuclide: str) -> Optional[dict]:
    """Найти запись паспорта (из COMMENT) для нужного нуклида.

    Использует уже распарсенные ``spec.extras['lsrm_passport']`` (их
    кладёт reader через ``parse_lsrm_passport_comment`` — НЕ дублируем
    парсер). Возвращает первую запись с совпадающим нуклидом.
    """
    pp = spec.extras.get("lsrm_passport") or []
    for e in pp:
        if e.get("nuclide", "").upper() == nuclide.upper():
            return e
    # fallback: единственная запись
    if len(pp) == 1:
        return pp[0]
    return None


def load_etalon_calibration_specs(
    *,
    repo_root: Path,
    nuclides: Optional[Sequence[str]] = None,
    etalon_files: Optional[Dict[str, str]] = None,
    windows: Sequence[Tuple[float, float, str]] = MES_WINDOWS,
) -> Tuple[Dict[str, CalibrationSpec], List[EtalonInfo]]:
    """Загрузить эталоны → {nuclide: CalibrationSpec} + список EtalonInfo.

    Активность каждого эталона берётся из COMMENT (паспорт, Бк/кг). Ось
    энергий каждого эталона — собственная (``channel_to_energy``), окна
    в keV.

    Возвращает (calibration_specs, etalon_infos).
    Бросает FileNotFoundError, если файл эталона отсутствует, и
    ValueError, если в COMMENT нет сертифицированной активности для
    нужного нуклида (не угадываем — anti-hallucination).
    """
    files = dict(etalon_files or ETALON_FILES)
    keys = list(nuclides) if nuclides is not None else list(files.keys())

    specs: Dict[str, CalibrationSpec] = {}
    infos: List[EtalonInfo] = []

    for nuc in keys:
        rel = files.get(nuc)
        if rel is None:
            raise ValueError(f"Нет эталонного файла в реестре для {nuc!r}")
        path = (repo_root / rel).resolve()
        if not path.exists():
            raise FileNotFoundError(f"Эталон не найден: {path}")
        spec = read_lsrm_spe(str(path))

        entry = _passport_for_nuclide(spec, nuc)
        if entry is None:
            raise ValueError(
                f"В COMMENT эталона {path.name} нет паспортной активности "
                f"для {nuc!r}; COMMENT={spec.comments!r}"
            )
        if not entry.get("is_specific_activity"):
            raise ValueError(
                f"Эталон {path.name}: активность {nuc!r} не удельная "
                f"(Бк/кг ожидается для Маринелли); raw={entry.get('raw_match')!r}"
            )
        A_cert = float(entry["value"])  # Бк/кг (value_Bq_per_kg для is_specific)
        sigma_rel = float(entry.get("uncertainty_pct", 5.0)) / 100.0

        # Ось энергий ЭТОГО эталона (собственная калибровка).
        n = int(spec.n_channels)
        energies = np.array(
            [spec.channel_to_energy(ch) for ch in range(n)], dtype=float,
        )

        specs[nuc] = CalibrationSpec(
            nuclide=nuc,
            counts=tuple(float(c) for c in np.asarray(spec.counts, dtype=float)),
            t_live_s=float(spec.live_time),
            A_certified_Bq_per_kg=A_cert,
            sigma_A_rel=sigma_rel,
        )
        # CalibrationSpec не хранит ось энергий — она нужна для R.
        # Прокинем через атрибут на info; build использует общий путь
        # ниже (per-spec energies).
        mat_name, mat_rho = matrix_descriptor(spec)
        infos.append(EtalonInfo(
            nuclide=nuc, path=str(path), t_live_s=float(spec.live_time),
            A_certified_Bq_per_kg=A_cert, sigma_A_rel=sigma_rel,
            cert_raw=str(entry.get("raw_match", "")),
            energy_cal=tuple(spec.energy_cal) if spec.energy_cal else (),
            matrix_name=mat_name, density_g_cm3=mat_rho,
        ))

    return specs, infos


def _spec_energies(spec) -> np.ndarray:
    n = int(spec.n_channels)
    return np.array([spec.channel_to_energy(ch) for ch in range(n)], dtype=float)


def build_R_matrix_per_spec(
    calibration_specs: Dict[str, CalibrationSpec],
    spec_energies: Dict[str, np.ndarray],
    bg_counts: Sequence[float],
    bg_t_live_s: float,
    bg_energies: Sequence[float],
    *,
    windows_keV: Sequence[Tuple[float, float]],
) -> SensitivityMatrix:
    """build_R_matrix, но с СОБСТВЕННОЙ осью энергий на каждый спектр.

    ``template_method.build_R_matrix`` принимает единую ``energies_per_ch``
    для всех спектров. В реальных данных Поверки 2024 калибровки эталонов
    различаются (проверено), поэтому здесь R считается тем же формулам
    LSRM §12.2:

        R_ki = (S_kr_i / t_kr − B_i / t_b) / A_k

    но S_kr_i суммируется по окнам в ОСИ ЭНЕРГИЙ k-го эталона, а B_i — в
    оси энергий фона. Окна одни и те же (keV). Это устраняет циркулярность
    общей сетки и численно эквивалентно при совпадающих калибровках.
    """
    from gamma.activity.template_method import sum_in_windows

    nuclides = list(calibration_specs.keys())
    K = len(nuclides)
    M = len(windows_keV)
    R = np.zeros((M, K), dtype=float)
    sigma_R = np.zeros((M, K), dtype=float)

    B_counts = sum_in_windows(bg_counts, bg_energies, windows_keV)
    B_rate = B_counts / bg_t_live_s
    sigma_B_rate = np.sqrt(np.maximum(B_counts, 1.0)) / bg_t_live_s

    for k, nuc in enumerate(nuclides):
        cspec = calibration_specs[nuc]
        en = spec_energies[nuc]
        S_counts = sum_in_windows(cspec.counts, en, windows_keV)
        S_rate = S_counts / cspec.t_live_s
        sigma_S_rate = np.sqrt(np.maximum(S_counts, 1.0)) / cspec.t_live_s

        if cspec.A_certified_Bq_per_kg <= 0:
            R[:, k] = np.nan
            sigma_R[:, k] = np.nan
            continue
        net_rate = S_rate - B_rate
        R[:, k] = net_rate / cspec.A_certified_Bq_per_kg
        sigma_net = np.sqrt(sigma_S_rate ** 2 + sigma_B_rate ** 2)
        rel_cert = cspec.sigma_A_rel
        sigma_R[:, k] = np.abs(R[:, k]) * np.sqrt(
            (sigma_net / np.where(net_rate != 0, net_rate, 1.0)) ** 2
            + rel_cert ** 2
        )

    return SensitivityMatrix(
        R=R, sigma_R=sigma_R,
        nuclides=tuple(nuclides),
        windows_keV=tuple(tuple(w) for w in windows_keV),
        notes="per-spec energy axes (F-448)",
    )


def build_sensitivity_matrix_from_etalons(
    *,
    repo_root: Path,
    bg_file: Optional[str] = None,
    nuclides: Optional[Sequence[str]] = None,
    etalon_files: Optional[Dict[str, str]] = None,
    windows: Sequence[Tuple[float, float, str]] = MES_WINDOWS,
) -> Tuple[SensitivityMatrix, List[EtalonInfo], dict]:
    """Полный билд R из эталонов + фона. Возвращает (R, infos, bg_meta)."""
    specs, infos = load_etalon_calibration_specs(
        repo_root=repo_root, nuclides=nuclides,
        etalon_files=etalon_files, windows=windows,
    )
    # пере-читаем эталоны ещё раз для осей энергий — дешево (1024 ch).
    files = dict(etalon_files or ETALON_FILES)
    spec_energies: Dict[str, np.ndarray] = {}
    for nuc in specs.keys():
        sp = read_lsrm_spe(str((repo_root / files[nuc]).resolve()))
        spec_energies[nuc] = _spec_energies(sp)

    bg_rel = bg_file or DEFAULT_BG_FILE
    bg_path = (repo_root / bg_rel).resolve()
    if not bg_path.exists():
        raise FileNotFoundError(f"Фон не найден: {bg_path}")
    bgspec = read_lsrm_spe(str(bg_path))
    bg_counts = tuple(float(c) for c in np.asarray(bgspec.counts, dtype=float))
    bg_energies = _spec_energies(bgspec)

    wins = windows_keV_only(windows)
    R = build_R_matrix_per_spec(
        specs, spec_energies, bg_counts, float(bgspec.live_time), bg_energies,
        windows_keV=wins,
    )
    bg_meta = {
        "path": str(bg_path),
        "t_live_s": float(bgspec.live_time),
        "comment": bgspec.comments,
    }
    return R, infos, bg_meta


__all__ = [
    "MES_WINDOWS",
    "windows_keV_only",
    "ETALON_FILES",
    "SAMPLE_TH232_FILE",
    "DEFAULT_BG_FILE",
    "EtalonInfo",
    "matrix_descriptor",
    "load_etalon_calibration_specs",
    "build_R_matrix_per_spec",
    "build_sensitivity_matrix_from_etalons",
]