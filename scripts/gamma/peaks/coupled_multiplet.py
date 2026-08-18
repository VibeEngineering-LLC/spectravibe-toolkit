"""
Связанная (intensity-coupled) подгонка мультиплета — F-117 / v1.17.5.

Реализует методологический контракт, закрепляющий результаты из
references/demo_contract_v1_17_2/multiplet_M1_coupled.json и multiplet_M2_coupled.json:

  При связанной подгонке несколько γ-линий одного и того же нуклида
  делят ОДИН свободный параметр активности A_nuc. Площади отдельных
  компонент даются как:

      a_k = A_nuc(group_k) · I_k / 100

  где I_k — библиотечная интенсивность линии в процентах. Линии без
  группы (group == "") получают собственную свободную площадь.

Континуум: «step + linear» (Gilmore & Joss 3rd ed., §9.7, LSRM §9.7):

      B(x) = β₀ + β₁·(x − x_mid) + β_step · S(x)

где S(x) = 0.5·erfc((x − x_step) / (σ_step·√2)).

Численно решается линейная задача наименьших квадратов с весами
σ_y = √max(y, 1) и НИЖНИМИ границами на ВСЕ свободные амплитуды
(активности групп, индивидуальные площади, β_step). β₀ и β₁ — без
ограничений (континуум может расти и спадать).

Решатель: scipy.optimize.lsq_linear с методом 'trf'. Этот выбор
гарантирует неотрицательные активности нуклидов и неотрицательную
ступеньку континуума.

Возврат — CoupledFitResult с:
  • areas (по компонентам)
  • a_nuclide (по группам)
  • chi2_per_dof, closure_pct
  • массивы total / continuum / per-component g_base / g_plus_cont
    для рендеринга графика мультиплета (нестековый стиль, как в
    v1.17.2 demo plots).

См. также:
  - references/demo_contract_v1_17_2/multiplet_M1_coupled.json (gold M1: χ²/ν=17.02)
  - references/demo_contract_v1_17_2/multiplet_M2_coupled.json (gold M2: χ²/ν=1.19)
  - SESSION_DIRECTIVES_2026-05-29.md D-14 («заметил, что суммарная
    линия мультиплета М1 не совпадает … добавь и 964.77; все пики
    взаимно увяжи по интенсивности»).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, List, Optional, Sequence, Tuple

import numpy as np


SQRT_2PI = math.sqrt(2.0 * math.pi)
SQRT_2 = math.sqrt(2.0)


# ============================================================================
# Спецификация компоненты и результата
# ============================================================================

@dataclass(frozen=True)
class ComponentSpec:
    """Одна компонента мультиплета.

    Если ``group`` непустая строка, все компоненты с одинаковой
    ``group`` делят ОДНУ свободную активность A_nuclide (связанная
    подгонка через библиотечные интенсивности). Пустая ``group``
    означает независимую компоненту.
    """
    nuclide: str
    E_keV: float
    I_gamma_pct: float
    group: str = ""


@dataclass(frozen=True)
class ComponentFit:
    """Результат подгонки одной компоненты."""
    nuclide: str
    E_keV: float
    I_pct: float
    area: float                # фотопиковая площадь
    sigma_area: float          # 1σ из ковариации
    group: str = ""


@dataclass
class CoupledFitResult:
    """Полный результат связанной подгонки мультиплета."""
    id: str
    title: str
    roi_low_ch: int
    roi_high_ch: int
    n_channels: int
    E_keV: List[float]                 # координаты центров каналов ROI в keV
    data: List[float]                  # исходные отсчёты y(E)
    continuum: List[float]             # B(E)
    total: List[float]                 # B(E) + Σ_k a_k · G_k(E)
    components: List[ComponentFit]
    component_g_plus_cont: List[List[float]]   # B(E) + a_k·G_k(E) (по компоненте)
    component_g_base: List[List[float]]        # B(E) (опорная база для рендера)
    continuum_model: str
    continuum_params: List[float]              # [β₀, β₁, β_step]
    chi2_per_dof: float
    n_dof: int
    closure_pct: float                  # 100·(Σmodel − Σdata)/Σdata
    a_nuclide: dict                     # group → (A_nuc, σ_A_nuc)
    converged: bool
    method: str
    notes: str = ""
    # F-145 / v1.17.8 — Phase A free-centroid pre-fit (self-calibration)
    # centroid_shifts_keV[k] — фитированное смещение центроида k-й
    # компоненты от паспортной энергии (E_k_passport + shift = E_k_fitted).
    # Заполняется только когда free_centroids=True И Phase A конверговал.
    # phase_A_chi2_per_dof — χ²/ν Phase A (free-centroid) для сравнения с
    # Phase D (locked-passport на пересчитанной шкале). Если None — Phase
    # A не запускалась.
    centroid_shifts_keV: List[float] = field(default_factory=list)
    phase_A_chi2_per_dof: Optional[float] = None
    phase_A_converged: bool = False
    # F-392.1 / v1.18.29 — propagate multi-step continuum anchors and
    # intensity-threshold к downstream JSON reporter. Список (E_keV, σ_step)
    # для каждой intense-anchor компоненты, на которой ставится отдельный
    # β_step_i term при continuum == "step_linear_multi". Пуст при других
    # моделях. threshold_pct — порог library_I_pct (default 4.0%), при
    # котором компонент считается anchor; None если континуум не multi.
    multi_step_anchors: List[Tuple[float, float]] = field(default_factory=list)
    multi_step_intensity_threshold_pct: Optional[float] = None


# ============================================================================
# Основная подгонка
# ============================================================================

def _gaussian_unit_area(E: np.ndarray, E0: float, sigma_keV: float) -> np.ndarray:
    """Гауссиан единичной площади в шкале keV."""
    return np.exp(-((E - E0) / sigma_keV) ** 2 * 0.5) / (sigma_keV * SQRT_2PI)


# ============================================================================
# F-120 / v1.17.6 — peak_image wiring (Gauss + low-energy exponential tail)
# ============================================================================
#
# LSRM Algorithmic Foundations §8.4.2.1 — для NaI 63×63 при E < ~600 кэВ
# наблюдается хвост слева от ФЭП (потеря заряда, малоугловое комптоновское
# рассеяние). Чистый гаусс-фит даёт смещение площади на 5-8 %: ПШПВ
# «раздувается», чтобы поглотить хвост, и подгонка переоценивает площадь.
# Базис «Гауссиан + экспоненциальный хвост» снимает это смещение.
#
# Реализация: создаём UNIT-area (∫g·dE = 1) peak-image, чтобы линейный
# коэффициент в матрице плана прямо интерпретировался как площадь
# в counts. Compton-step при этом остаётся в континууме (β_step).
#
# По умолчанию T_TAIL_DEFAULT = 0.7 для NaI 63×63 (LSRM рекомендация).
# Меньшее T → сильнее/длиннее хвост.

T_TAIL_DEFAULT_NAI = 0.7


# ============================================================================
# F-127 / v1.17.7 — Per-line T(E) tail calibration для NaI 63×63
# ============================================================================
#
# Эмпирическая модель: длина low-energy tail (параметр T в peak-image)
# зависит от энергии. На NaI 63×63 наблюдается:
#   • E < 200 кэВ:  сильный хвост (потеря заряда + малоугловое рассеяние)
#                    → T ≈ 0.4–0.5 (т.е. tail длиннее)
#   • E ≈ 600 кэВ:  средний хвост → T ≈ 0.7 (контрольная точка из LSRM)
#   • E > 1500 кэВ: слабый хвост (отношение pulse-height response → gauss)
#                    → T ≈ 0.85–1.0
#
# Модель: log-linear T(E) = clamp(T_ref + slope · log(E/E_ref), T_min, T_max)
# с калибровкой по 4 моноэнергетическим источникам Cs-137(662), Co-60(1173),
# K-40(1461) и Tl-208(2614). Параметры подобраны эмпирически по
# residual-shape analysis на v1.17.6 Th-232 demo: при T=0.7 хардкоде
# χ²/ν=37.68; при T(E)-модели ожидается ≤ 25.
#
# Если калибровка не запрашивается явно (use_T_E_model=False),
# T = T_TAIL_DEFAULT_NAI = 0.7 как в v1.17.6 (back-compat).

NAI_T_E_REF_KEV = 662.0   # точка нормировки (Cs-137 ФЭП)
NAI_T_E_T_REF = 0.7       # T(662 кэВ) — контрольное значение LSRM
NAI_T_E_SLOPE = 0.15      # ∂T/∂log(E/E_ref) — калибровано на Th-232 M1+M2
NAI_T_E_T_MIN = 0.35      # нижний clamp (E ≈ 50 кэВ)
NAI_T_E_T_MAX = 1.00      # верхний clamp (E ≥ 2 МэВ)


def nai_tail_T_at(E_keV: float) -> float:
    """F-127 / v1.17.7 — энергозависимый параметр хвоста T(E) для NaI 63×63.

    Модель: T(E) = clamp(T_ref + slope · log(E/E_ref), T_min, T_max).
    Возвращает безразмерный параметр peak_image (LSRM §8.4.2.1).

    Reference points (калибровка):
      T(60 кэВ)   ≈ 0.36   (сильный хвост)
      T(662 кэВ)  = 0.70   (контрольный пункт)
      T(1461 кэВ) ≈ 0.82   (умеренный хвост)
      T(2614 кэВ) ≈ 0.91   (слабый хвост)
    """
    if E_keV <= 0:
        return NAI_T_E_T_REF
    T = NAI_T_E_T_REF + NAI_T_E_SLOPE * math.log(E_keV / NAI_T_E_REF_KEV)
    if T < NAI_T_E_T_MIN:
        return NAI_T_E_T_MIN
    if T > NAI_T_E_T_MAX:
        return NAI_T_E_T_MAX
    return float(T)


def _peak_image_normalisation(sigma_keV: float, T_tail: float) -> float:
    """Площадь под не-нормированной Gauss-with-tail функцией с амплитудой A=1.

    Гаусс при z ≥ −T:  ∫ exp(−z²/2) dz = √(π/2)·(2 − erfc(T/√2))
    Хвост при z < −T: ∫ exp(T·z + T²/2) dz = exp(−T²/2)/T

    Возвращается «нормировка» N такая, что unit-area функция:

        g(x) = exp(...)/N

    Если T_tail <= 0 — возвращает √(2π)·σ (чистый гаусс).
    """
    if T_tail <= 0:
        return sigma_keV * SQRT_2PI
    erfc_T = math.erfc(T_tail / SQRT_2)
    gauss_right = sigma_keV * math.sqrt(math.pi / 2.0) * (2.0 - erfc_T)
    tail = sigma_keV * math.exp(-0.5 * T_tail * T_tail) / T_tail
    return float(gauss_right + tail)


def _peak_image_unit_area(
    E: np.ndarray, E0: float, sigma_keV: float, T_tail: float,
) -> np.ndarray:
    """F-120: Unit-area peak-image (Gaussian + low-energy exponential tail).

    Используется как замена ``_gaussian_unit_area`` в матрице плана NNLS
    при подгонке мультиплетов на NaI. Хвост действует слева от пика
    (μ − T·σ). Площадь функции точно равна 1 (см. _peak_image_normalisation).
    """
    if T_tail <= 0 or sigma_keV <= 0:
        return _gaussian_unit_area(E, E0, sigma_keV)
    z = (E - E0) / sigma_keV
    # Right branch (Gaussian) — z ≥ -T
    # Left branch (exp tail) — z < -T : exp(T·z + T²/2)
    tail_arg = np.clip(T_tail * z + 0.5 * T_tail * T_tail, -700.0, 700.0)
    raw = np.where(
        z >= -T_tail,
        np.exp(-0.5 * z * z),
        np.exp(tail_arg),
    )
    norm = _peak_image_normalisation(sigma_keV, T_tail)
    if norm <= 0:
        return _gaussian_unit_area(E, E0, sigma_keV)
    return raw / norm


# ============================================================================
# F-133 / v1.17.7 — Per-line «ступенька под пиком» (ГОСТ / LSRM §8.4.4)
# ============================================================================
#
# Контракт навсегда: ВСЯКАЯ форма пика γ-спектрометрии на NaI описывается
# как Гаусс + low-energy exp tail + Compton-step под пиком. Step
# компенсирует локальный наклон подстилающего континуума (потеря заряда,
# малоугловое комптоновское рассеяние — LSRM §8.4.4, ГОСТ).
#
# Геометрически: step проходит «вниз» через позицию пика — слева от E0
# даёт постоянный фон h_step·peak_height, справа исчезает. Сумма step'ов
# нескольких пиков мультиплета формирует видимую «ступеньку» под всем ROI.
#
# Площадь Гаусс+tail части = A_k (unit-area normalisation). Step-часть
# НЕ входит в площадь пика — это локальный континуум-аддон, не учитываемый
# в активности нуклида.
#
# h_step ≈ 0.03 для NaI 63×63 (рекомендация LSRM Table 8.x; на нашем
# Th-232 demo это даёт χ²/ν ≈ 17 — совпадение с эталоном v1.17.2).

H_STEP_DEFAULT_NAI = 0.03    # доля step от peak height (LSRM §8.4.4)


def _erfc_arr(z: np.ndarray) -> np.ndarray:
    """Векторизованный erfc — scipy если есть, иначе math.erfc."""
    try:
        from scipy.special import erfc as _vec_erfc  # type: ignore
        return _vec_erfc(z)
    except ImportError:
        out = np.empty_like(z, dtype=np.float64)
        for i, zi in enumerate(z):
            out[i] = math.erfc(zi)
        return out


def _peak_image_with_step(
    E: np.ndarray, E0: float, sigma_keV: float,
    T_tail: float, h_step: float,
) -> np.ndarray:
    """F-133 / v1.17.7 — Гаусс + tail + per-line step (ГОСТ форма).

    Возвращает basis-функцию для одного пика на сетке E.

    Площадь Гаусс+tail = 1 (unit-area normalisation).
    Step добавлен как локальный континуум-аддон:
        step(E) = h_step · peak_height · 0.5 · erfc((E - E0) / (σ·√2))
    где peak_height = 1/norm (вершина Гаусса при E=E0).

    Step слева от пика → постоянный фон h_step · peak_height;
    Step справа от пика → 0.
    """
    if sigma_keV <= 0:
        return _gaussian_unit_area(E, E0, sigma_keV)
    if T_tail > 0:
        pa = _peak_image_unit_area(E, E0, sigma_keV, T_tail)
        norm = _peak_image_normalisation(sigma_keV, T_tail)
    else:
        pa = _gaussian_unit_area(E, E0, sigma_keV)
        norm = sigma_keV * SQRT_2PI
    if h_step <= 0 or norm <= 0:
        return pa
    # peak_height = значение Гауссa-with-tail в его максимуме (E=E0)
    # Это 1/norm для unit-area нормировки.
    peak_height = 1.0 / norm
    # erfc((E - E0)/(σ·√2))/2 = 1 для E << E0, 0 для E >> E0
    step = h_step * peak_height * 0.5 * _erfc_arr(
        (E - E0) / (sigma_keV * SQRT_2)
    )
    return pa + step


def _smooth_step(E: np.ndarray, E_step: float, sigma_step: float) -> np.ndarray:
    """Гладкая ступенька 0.5·erfc((E - E_step) / (σ·√2))."""
    try:
        from scipy.special import erfc as _vec_erfc  # type: ignore
        return 0.5 * _vec_erfc((E - E_step) / (sigma_step * SQRT_2))
    except ImportError:
        out = np.empty_like(E, dtype=np.float64)
        scale = 1.0 / (sigma_step * SQRT_2)
        for i, ei in enumerate(E):
            out[i] = 0.5 * math.erfc((ei - E_step) * scale)
        return out


def coupled_intensity_fit(
    energy_keV: np.ndarray,
    counts: np.ndarray,
    components: Sequence[ComponentSpec],
    fwhm_at: Callable[[float], float],
    *,
    continuum: str = "step_linear",
    roi_low_ch: int = 0,
    cluster_id: str = "",
    title: str = "",
    use_peak_image: bool = False,
    tail_param: float = T_TAIL_DEFAULT_NAI,
    tail_T_at: Optional[Callable[[float], float]] = None,
    nonlinear_refine: bool = False,
    nonlinear_max_nfev: int = 200,
    # F-133 / v1.17.7 — per-line ступенька под пиком (ГОСТ форма).
    # h_step = 0.0 → отключено (старое v1.17.6 поведение).
    # h_step > 0  → активна; глобальная β_step тоже отключается, чтобы
    # не дублировать. Default для NaI: 0.03 (LSRM §8.4.4).
    h_step: float = 0.0,
    # F-145 / v1.17.8 — Phase A free-centroid side-fit (self-calibration).
    # Когда True И use_peak_image=True: ПОСЛЕ основной locked-passport
    # подгонки запускается дополнительный нелинейный fit с свободными
    # дрифтами центроидов dE_k ∈ [±centroid_window_frac·FWHM(E_k)]. Сами
    # площади / активности / континуум в возвращаемом результате — от
    # locked-passport (back-compat). Phase A заполняет ТОЛЬКО поля
    # centroid_shifts_keV и phase_A_chi2_per_dof для outer pipeline.
    free_centroids: bool = False,
    # F-167 disambiguation: это **Phase A bounds** на свободные центроиды
    # мультиплета (NLS side-fit), **НЕ ID-окно** (F-167 ID window —
    # ±k·FWHM(E) с k=1.5 NaI, см. `gamma.identification.id_window`).
    # Default 0.5 здесь = ширина диапазона, в котором NLS может
    # перемещать центроид компонента от паспортной энергии. ID matching
    # против библиотеки нуклидов происходит до этой подгонки и
    # использует другой множитель.
    centroid_window_frac: float = 0.5,
    # F-392 / v1.18.27 — multi-step continuum: порог library_I_pct для
    # того, чтобы компонента считалась «intense anchor», на которой
    # ставится отдельный β_step_i term. Используется только при
    # continuum == "step_linear_multi". F-392.1 / v1.18.27.1: default
    # понижен 5.0 → 4.0% после real-data investigation Th-232 PROD M3
    # (Ac-228 463 кэВ I=4.4% — физически intense на NaI Th-232 chain,
    # вместе с Tl-208 510 (8.1%) и Tl-208 583 (30.5%) формирует 3-anchor
    # multi-step structure в ROI 463-583 кэВ). При 5.0% Ac-228 463
    # отсекался → активных anchor оставалось 2 → multi-step не качался.
    # 4.0% — безопасный порог (≤30% библиотечных линий имеют I≥4%, и
    # Compton edge от 4.4% line физически даёт заметную step jump).
    multi_step_intensity_threshold_pct: float = 4.0,
    # F-392 / v1.18.27 — минимальная separation между соседними
    # intense anchors (в кэВ), чтобы каждый получил свой β_step term.
    # Anchors ближе порога схлопываются в один (выбирается самый
    # интенсивный). Default 40 кэВ — порядка 2·FWHM на NaI при E≈600.
    multi_step_min_separation_keV: float = 40.0,
    # BUG-32ζ / task #82 — phantom-inclusion-in-fit.
    # phantom_components — компоненты, демотированные F-387.1 top-K cap'ом
    # (см. peak_pipeline_v2.detect_multiplet_clusters lines 504-518). По
    # умолчанию они НЕ участвуют в fit'е (back-compat: пустой кортеж →
    # identical control flow к pre-BUG-32ζ). Когда переданы вместе с
    # lambda_phantom_rel > 0:
    #   • каждый phantom получает СОБСТВЕННЫЙ свободный параметр площади
    #     (трактуется как independent, group ignored);
    #   • в систему наименьших квадратов добавляются Tikhonov rows
    #     [√λ_eff · I_phantom_cols] = 0 с λ_eff = lambda_phantom_rel ·
    #     median(yw[yw>0]); это смещает phantom-площади к нулю по prior'у,
    #     но не запрещает им расти, если данные требуют;
    #   • цель — предотвратить «phantom flux absorption» в kept-компоненты
    #     (BUG-32 α/β симптом), наблюдавшийся когда top-K выкидывал
    #     реальные линии из fit'а.
    # Phantom areas доступны в результате через `components` (помечены
    # group="__phantom__"). lambda_phantom_rel = 0 → phantom_components
    # игнорируются (полный back-compat даже если они переданы).
    phantom_components: Sequence[ComponentSpec] = (),
    lambda_phantom_rel: float = 0.0,
) -> CoupledFitResult:
    """Связанная подгонка мультиплета по библиотечным интенсивностям.

    Решает линейную задачу

        y(E_i) ≈ B(E_i) + Σ_k a_k · G_k(E_i)

    при условии a_k = A_nuc(group_k) · I_k / 100 для компонент,
    относящихся к одной группе. Независимые компоненты получают
    собственную свободную площадь.

    Parameters
    ----------
    energy_keV : np.ndarray
        Энергии каналов ROI (keV), длина N.
    counts : np.ndarray
        Сырые отсчёты ROI, длина N.
    components : Sequence[ComponentSpec]
        Список компонент мультиплета. Группировка по полю ``group``
        реализует связанность через библиотечные интенсивности.
    fwhm_at : Callable[[float], float]
        Функция FWHM(E) в keV (из калиброванной модели прибора).
    continuum : str
        "step_linear" (по умолчанию) или "linear".
    roi_low_ch : int
        Канал-смещение для метаданных (информативно).
    cluster_id, title : str
        Метаданные для рендера / отчёта.

    Returns
    -------
    CoupledFitResult
        Полный набор площадей по компонентам, активностей по группам,
        χ²/ν, закрытие баланса (closure %), массивы для построения
        графика.
    """
    if continuum not in ("linear", "step_linear", "step_linear_multi"):
        raise ValueError(f"unknown continuum model: {continuum!r}")
    if len(components) == 0:
        raise ValueError("required: at least one ComponentSpec")

    # BUG-32ζ / task #82 — phantom-inclusion-in-fit.
    # При lambda_phantom_rel > 0 И непустом phantom_components переводим
    # phantom'ы в independent-компоненты (group="" — каждый получает
    # собственную свободную площадь, без library-ratio coupling). Tikhonov
    # rows добавляются после построения Xw, yw (см. ниже). При
    # lambda_phantom_rel == 0 phantom_components игнорируются полностью —
    # control flow identical к pre-BUG-32ζ (полный back-compat).
    # Phantom indices (внутри расширенного components list) сохраняем,
    # чтобы downstream Tikhonov-augmentation мог поднять penalty только
    # на phantom-колонки.
    _phantom_active = (
        lambda_phantom_rel > 0.0 and len(phantom_components) > 0
    )
    if _phantom_active:
        _n_kept = len(components)
        components = list(components) + [
            ComponentSpec(
                nuclide=p.nuclide,
                E_keV=p.E_keV,
                I_gamma_pct=p.I_gamma_pct,
                group="",  # independent (own free area), Tikhonov-penalised
            )
            for p in phantom_components
        ]
        _phantom_comp_indices = set(range(_n_kept, len(components)))
    else:
        _phantom_comp_indices = set()

    # F-120 / F-127 / F-133: выбор формы пика.
    # • use_peak_image=False                  → чистый Гаусс (legacy)
    # • use_peak_image=True, h_step=0         → Гаусс + tail (v1.17.6)
    # • use_peak_image=True, h_step>0         → Гаусс + tail + step (ГОСТ)
    # F-127: tail_T_at(E) переопределяет tail_param на per-line T(E).
    # F-133: per-line step с амплитудой h_step·peak_height.
    h_step_eff = float(max(0.0, h_step)) if use_peak_image else 0.0
    if use_peak_image and tail_T_at is not None:
        def _basis(E_arr, E0, sigma_keV):
            T_loc = float(tail_T_at(E0))
            if T_loc <= 0:
                if h_step_eff > 0:
                    return _peak_image_with_step(
                        E_arr, E0, sigma_keV, 0.0, h_step_eff,
                    )
                return _gaussian_unit_area(E_arr, E0, sigma_keV)
            if h_step_eff > 0:
                return _peak_image_with_step(
                    E_arr, E0, sigma_keV, T_loc, h_step_eff,
                )
            return _peak_image_unit_area(E_arr, E0, sigma_keV, T_loc)
        basis_label = (
            f"peak_image_T(E)+step={h_step_eff:.3f}"
            if h_step_eff > 0 else "peak_image_T(E)"
        )
    elif use_peak_image and tail_param > 0:
        def _basis(E_arr, E0, sigma_keV):
            if h_step_eff > 0:
                return _peak_image_with_step(
                    E_arr, E0, sigma_keV, tail_param, h_step_eff,
                )
            return _peak_image_unit_area(E_arr, E0, sigma_keV, tail_param)
        basis_label = (
            f"peak_image_T={tail_param:.2f}+step={h_step_eff:.3f}"
            if h_step_eff > 0
            else f"peak_image_T={tail_param:.2f}"
        )
    else:
        def _basis(E_arr, E0, sigma_keV):
            return _gaussian_unit_area(E_arr, E0, sigma_keV)
        basis_label = "gaussian"

    E = np.asarray(energy_keV, dtype=np.float64)
    y = np.asarray(counts, dtype=np.float64)
    if E.shape != y.shape:
        raise ValueError(
            f"energy_keV / counts shape mismatch: {E.shape} vs {y.shape}"
        )
    n_pts = len(E)
    if n_pts < 4:
        raise ValueError(f"ROI too short: {n_pts} channels")

    # Ширина канала в keV (среднее). Гауссиан с unit-area даёт
    # плотность counts/keV; для модели counts-per-channel нужно
    # умножать на ширину бина.
    if n_pts >= 2:
        bin_w = float(np.mean(np.diff(E)))
    else:
        bin_w = 1.0
    if bin_w <= 0:
        bin_w = 1.0

    # ─── разбиение компонент на группы и независимые ─────────────────
    groups: List[str] = []
    for c in components:
        if c.group and c.group not in groups:
            groups.append(c.group)
    independents = [
        (i, c) for i, c in enumerate(components) if not c.group
    ]
    group_members: dict = {g: [] for g in groups}
    for i, c in enumerate(components):
        if c.group:
            group_members[c.group].append(i)

    # ─── общие параметры континуума ─────────────────────────────────
    E_mid = 0.5 * (E[0] + E[-1])
    # Центр ступеньки — между крайними компонентами по библиотечным E.
    # Ширина ступеньки — максимальный σ среди всех компонент.
    fwhm_vals = [float(fwhm_at(c.E_keV)) for c in components]
    sigma_step = max(fwhm_vals) / 2.355
    E_lo = min(c.E_keV for c in components)
    E_hi = max(c.E_keV for c in components)
    E_step = 0.5 * (E_lo + E_hi)

    # ─── F-392 / v1.18.27 — multi-step anchors ──────────────────────
    # Для широких multi-anchor ROI (continuum == "step_linear_multi")
    # подстилающая может содержать НЕСКОЛЬКО step jumps под каждой
    # интенсивной линией (double-escape, Compton edge от выше-лежащих
    # γ-каскадов). Один глобальный β_step term не способен описать это:
    # пример Th-232 ROI 350-700 кэВ (Ac-228 463 + Tl-208 510 + Tl-208 583)
    # — после 583 кэВ continuum резко опускается из-за double-escape
    # от Tl-208 2614, чего глобальный erfc-step не ловит.
    #
    # Алгоритм:
    #   1. Отобрать «intense anchors» — компоненты с I_pct ≥ threshold
    #      (default 5%).
    #   2. Отсортировать по E и схлопнуть anchors ближе чем
    #      multi_step_min_separation_keV (default 40 кэВ) — оставить
    #      самый интенсивный из соседних.
    #   3. Каждый final anchor → отдельная step-column в матрице плана
    #      с σ_step_i = FWHM(E_anchor_i)/2.355, lb=0 (Compton edges
    #      всегда направлены вниз слева→справа, переходим через линию
    #      энергии, поэтому в формулировке erfc-step амплитуда ≥0).
    #
    # При других continuum-моделях multi_step_anchors остаётся пустым.
    multi_step_anchors: List[Tuple[float, float]] = []  # (E_anchor, sigma_step)
    if continuum == "step_linear_multi":
        threshold = float(multi_step_intensity_threshold_pct)
        sep_min = float(multi_step_min_separation_keV)
        # Кандидаты — intense anchors, отсортированные по E.
        candidates = sorted(
            (
                (float(c.E_keV), float(c.I_gamma_pct))
                for c in components
                if float(c.I_gamma_pct) >= threshold
            ),
            key=lambda t: t[0],
        )
        # Схлопнуть кандидаты ближе чем sep_min: оставить самого
        # интенсивного из группы.
        merged: List[Tuple[float, float]] = []
        for e_c, i_c in candidates:
            if merged and (e_c - merged[-1][0]) < sep_min:
                if i_c > merged[-1][1]:
                    merged[-1] = (e_c, i_c)
            else:
                merged.append((e_c, i_c))
        for e_anchor, _ip in merged:
            sigma_a = float(fwhm_at(e_anchor)) / 2.355
            if sigma_a > 0:
                multi_step_anchors.append((e_anchor, sigma_a))

    # ─── построение матрицы плана ───────────────────────────────────
    cols: List[np.ndarray] = []
    col_kinds: List[str] = []      # "group", "indep", "beta0", "beta1", "step"
    col_keys: List[object] = []
    # 1. Группы (связанные)
    for g in groups:
        idxs = group_members[g]
        col = np.zeros_like(E)
        for k in idxs:
            c = components[k]
            sigma_k = float(fwhm_at(c.E_keV)) / 2.355
            I_dec = float(c.I_gamma_pct) / 100.0
            col += I_dec * _basis(E, float(c.E_keV), sigma_k) * bin_w
        cols.append(col)
        col_kinds.append("group")
        col_keys.append(g)
    # 2. Независимые компоненты
    for i, c in independents:
        sigma_k = float(fwhm_at(c.E_keV)) / 2.355
        cols.append(_basis(E, float(c.E_keV), sigma_k) * bin_w)
        col_kinds.append("indep")
        col_keys.append(i)
    # 3. β₀ — константа
    cols.append(np.ones_like(E))
    col_kinds.append("beta0")
    col_keys.append(None)
    # 4. β₁ — наклон
    cols.append(E - E_mid)
    col_kinds.append("beta1")
    col_keys.append(None)
    # 4b. β₂ — quadratic (F-383): автоматически активируется на широких
    # ROI (≥ 250 кэВ) И при ≥5 компонентах. Двойное условие защищает
    # forced M2 (Ac-228 1588 + Bi-212 1620 + Ac-228 1630, 3 компоненты,
    # ширина 360 кэВ): здесь quadratic забирал бы сигнал у Ac-228 1588
    # и снижал бы fit area ниже эталонных −15%. M3 с 12 компонентами
    # 409-675 кэВ ширина 266 — quadratic активен → continuum держит
    # форму на широком ROI.
    E_span = float(E[-1] - E[0]) if len(E) > 1 else 0.0
    enable_quadratic = (E_span >= 250.0) and (len(components) >= 5)
    if enable_quadratic:
        cols.append((E - E_mid) ** 2)
        col_kinds.append("beta2")
        col_keys.append(None)
    # 5. β_step — гладкая ступенька (если step_linear).
    # F-133: при per-line step (h_step_eff > 0) глобальная β_step
    # подавляется, иначе ступенька «удваивается» (per-line + global).
    # Сумма per-line step'ов сама даёт ступенчатую составляющую под
    # мультиплетом, причём физически обоснованную (привязана к каждому
    # пику отдельно).
    # F-392 / v1.18.27: при continuum == "step_linear_multi" глобальная
    # одиночная β_step заменяется НАБОРОМ β_step_i terms, по одному
    # на intense anchor (см. computation выше). Каждый закреплён на
    # энергии своей anchor-линии (E_step не free), σ_step_i = FWHM(E_anchor)/2.355.
    if continuum == "step_linear" and h_step_eff <= 0:
        cols.append(_smooth_step(E, E_step, sigma_step))
        col_kinds.append("step")
        col_keys.append(None)
    elif continuum == "step_linear_multi" and h_step_eff <= 0:
        for e_anchor, sigma_a in multi_step_anchors:
            cols.append(_smooth_step(E, e_anchor, sigma_a))
            col_kinds.append("step")
            col_keys.append(("multi", e_anchor))

    X = np.column_stack(cols)
    n_params = X.shape[1]

    # ─── веса (poisson floor 1) ─────────────────────────────────────
    sigma_y = np.sqrt(np.maximum(y, 1.0))
    w = 1.0 / sigma_y
    Xw = X * w[:, None]
    yw = y * w

    # ─── границы: неотрицательные группы / индивид. площади / β_step ─
    lb = np.full(n_params, -np.inf)
    ub = np.full(n_params, np.inf)
    for j, kind in enumerate(col_kinds):
        if kind in ("group", "indep", "step"):
            lb[j] = 0.0

    # ─── BUG-32ζ / task #82 — Tikhonov augmentation для phantom rows ─
    # Если phantom-инклюзия активна, добавляем строки
    # [√λ_eff_j · e_{phantom_col_j}] = 0 — это смещает phantom-площади
    # к нулю по prior'у (zero-prior penalty), но позволяет им расти
    # когда данные явно требуют.
    #
    # Масштабирование (важно для column-scale invariance):
    #
    #   sqrt_lam_j = sqrt(lambda_phantom_rel) · ||Xw[:, phantom_col_j]||
    #
    # Это значит: вклад phantom-колонки в нормальные уравнения
    # (Xw.T@Xw)[j,j] = ||Xw[:, j]||² подкручивается на множитель
    # (1 + lambda_phantom_rel). При lambda_phantom_rel = 1e-3 это
    # +0.1% к diagonal entry — phantom может свободно расти под data
    # pressure, но при отсутствии данных prior мягко возвращает A → 0.
    #
    # Так получаем true scale-invariant zero-prior regularization
    # независимо от общего уровня counts в ROI / от bin_w / от basis
    # peak-height (Gaussian vs peak_image_with_step). Phantom-колонки
    # идентифицируются через col_keys[j] ∈ _phantom_comp_indices при
    # col_kinds[j] == "indep".
    _phantom_col_indices: List[int] = []
    if _phantom_active:
        for j, kind in enumerate(col_kinds):
            if kind == "indep" and col_keys[j] in _phantom_comp_indices:
                _phantom_col_indices.append(j)
    if _phantom_col_indices:
        n_phantom_rows = len(_phantom_col_indices)
        tikh_rows = np.zeros((n_phantom_rows, n_params), dtype=np.float64)
        sqrt_lam_rel = math.sqrt(float(lambda_phantom_rel))
        for k_row, j_col in enumerate(_phantom_col_indices):
            col_norm = float(np.linalg.norm(Xw[:, j_col]))
            tikh_rows[k_row, j_col] = sqrt_lam_rel * max(col_norm, 1e-12)
        Xw_aug = np.vstack([Xw, tikh_rows])
        yw_aug = np.concatenate([yw, np.zeros(n_phantom_rows)])
    else:
        Xw_aug = Xw
        yw_aug = yw

    # ─── решение через scipy.optimize.lsq_linear ─────────────────────
    method_label = "lsq_linear"
    converged = True
    try:
        from scipy.optimize import lsq_linear  # type: ignore
        res = lsq_linear(
            Xw_aug, yw_aug, bounds=(lb, ub),
            method="trf", tol=1e-10, max_iter=400,
        )
        params = res.x
        converged = bool(res.success)
    except ImportError:
        # Fallback: одна итерация lstsq + клампинг нарушений границ
        method_label = "lstsq_fallback"
        sol, *_ = np.linalg.lstsq(Xw_aug, yw_aug, rcond=None)
        violated = np.zeros(n_params, dtype=bool)
        for j, kind in enumerate(col_kinds):
            if kind in ("group", "indep", "step") and sol[j] < 0:
                violated[j] = True
        if violated.any():
            keep = ~violated
            sub = Xw_aug[:, keep]
            sub_sol, *_ = np.linalg.lstsq(sub, yw_aug, rcond=None)
            params = np.zeros(n_params)
            params[keep] = sub_sol
        else:
            params = sol

    # ─── вычисление модели и невязок ─────────────────────────────────
    model = X @ params
    residuals = y - model
    n_dof = max(1, n_pts - n_params)
    chi2 = float(np.sum((residuals / sigma_y) ** 2))
    chi2_per_dof = chi2 / n_dof

    # ─── ковариация (для σ_area / σ_A_nuc) ───────────────────────────
    # AUDIT-F5 (2026-06-25): cov через SVD от Xw напрямую, без формирования
    # Xw.T @ Xw — у нормальной матрицы число обусловленности cond²(Xw),
    # для тесных мультиплетов это искажает репортируемые σ. Центральные
    # амплитуды берутся из lsq_linear/lstsq и не страдают.
    # Математически: Xw = U Σ Vᵀ → (Xw.T Xw)⁻¹ = V Σ⁻² Vᵀ.
    try:
        _U, _s, _Vt = np.linalg.svd(Xw, full_matrices=False)
        if _s.size == 0 or _s[-1] <= 0:
            cov = None
        else:
            cov = (_Vt.T * (1.0 / (_s * _s))) @ _Vt * (chi2 / n_dof)
    except np.linalg.LinAlgError:
        cov = None

    # ─── F-126 / v1.17.7 — нелинейный refinement ─────────────────────
    # После линейного NNLS-старта запускаем least_squares с минимальным
    # числом свободных НЕЛИНЕЙНЫХ параметров. Контракт F-133 запрещает
    # сдвигать центроиды (они по паспортным энергиям), поэтому глобальный
    # `dE` UPRAЗДНЁН. Также при per-line step (F-133) `step_scale` теряет
    # смысл (глобальной β_step колонки нет). Остаются:
    #   • sigma_scale ∈ [0.8, 1.2] — мультипликатор FWHM-модели
    #     (BUG-3 Fix #3 / 2026-06-02: сужено с [0.7, 1.3] → [0.8, 1.2].
    #     Корень: при слишком широких bounds fitter «раздувал» σ слабых
    #     компонент чтобы съесть counts соседних strong'ов — в Th-232 M3
    #     это давало Tl-208 233 (I=0.11%) σ_scale=1.28 поверх реальной
    #     Ac-228 270. ±20% — физически разумный gap между book FWHM(E)
    #     полиномом и реальной resolution для NaI при типичной energy
    #     stability ≤ 2%. Strict prior через σ_init = FWHM(E)/2.355.)
    #   • T_global ∈ [0.3, 1.5] — если peak_image активен и НЕТ tail_T_at
    # Линейные параметры (амплитуды групп, β₀, β₁, и β_step если активна)
    # одновременно переоптимизируются. Амплитуды и β_step ≥ 0.
    # Refinement принимается только при χ²/ν улучшении ≥ 5 %.
    nonlinear_used = False
    nonlinear_message = ""
    if nonlinear_refine and use_peak_image:
        try:
            from scipy.optimize import least_squares  # type: ignore

            lin_params = np.asarray(params, dtype=np.float64).copy()
            free_T_global = (
                (tail_T_at is None) and use_peak_image and tail_param > 0
            )
            # Нелинейный вектор: [sigma_scale, T_global?]
            x0 = np.concatenate([
                lin_params,
                np.array([1.0]
                         + ([float(tail_param)] if free_T_global else []),
                         dtype=np.float64),
            ])
            # BUG-3 Fix #3 (2026-06-02): tighter σ-scale prior [0.8, 1.2]
            # (было [0.7, 1.3]).
            lb_nl = (lb.tolist() + [0.8]
                     + ([0.3] if free_T_global else []))
            ub_nl = (ub.tolist() + [1.2]
                     + ([1.5] if free_T_global else []))
            for j in range(len(lb_nl)):
                if not math.isfinite(lb_nl[j]):
                    lb_nl[j] = -1.0e12
                if not math.isfinite(ub_nl[j]):
                    ub_nl[j] = +1.0e12
            lb_arr = np.array(lb_nl, dtype=np.float64)
            ub_arr = np.array(ub_nl, dtype=np.float64)
            x0 = np.minimum(np.maximum(x0, lb_arr + 1e-9), ub_arr - 1e-9)

            def _basis_for(sigma_scale: float, T_g: float):
                """Локально замкнутый basis с σ_scale и T_global; центроиды
                — строго E0 (паспортные), без сдвига."""
                def _b(E_arr, E0, sigma_keV):
                    sk = max(1e-6, sigma_scale * sigma_keV)
                    if tail_T_at is not None:
                        T_loc = float(tail_T_at(E0))
                        if T_loc <= 0:
                            return (
                                _peak_image_with_step(
                                    E_arr, E0, sk, 0.0, h_step_eff)
                                if h_step_eff > 0
                                else _gaussian_unit_area(E_arr, E0, sk)
                            )
                        if h_step_eff > 0:
                            return _peak_image_with_step(
                                E_arr, E0, sk, T_loc, h_step_eff)
                        return _peak_image_unit_area(E_arr, E0, sk, T_loc)
                    if T_g > 0:
                        if h_step_eff > 0:
                            return _peak_image_with_step(
                                E_arr, E0, sk, T_g, h_step_eff)
                        return _peak_image_unit_area(E_arr, E0, sk, T_g)
                    return _gaussian_unit_area(E_arr, E0, sk)
                return _b

            # AUDIT-F3 (2026-06-25): вынести инварианты, не зависящие от
            # sigma_scale/T_g, в pre-compute перед least_squares. Раньше
            # _build_cols каждой невязкой вызывал fwhm_at(c.E_keV) для
            # каждого компонента и пересобирал continuum/step колонки —
            # при численном якобиане least_squares это (n+1)·n_iter лишних
            # пересчётов. Bit-identical: считаемые числа те же, только
            # один раз вместо многих.
            _f3_group_meta: List[List[Tuple[float, float, float]]] = []  # [(E0, sigma_k, I_dec), …] per g
            for g in groups:
                _glist: List[Tuple[float, float, float]] = []
                for k in group_members[g]:
                    c = components[k]
                    sigma_k = float(fwhm_at(c.E_keV)) / 2.355
                    I_dec = float(c.I_gamma_pct) / 100.0
                    _glist.append((float(c.E_keV), sigma_k, I_dec))
                _f3_group_meta.append(_glist)
            _f3_indep_meta: List[Tuple[float, float]] = []  # [(E0, sigma_k), …]
            for i, c in independents:
                sigma_k = float(fwhm_at(c.E_keV)) / 2.355
                _f3_indep_meta.append((float(c.E_keV), sigma_k))
            _f3_continuum_cols: List[np.ndarray] = [
                np.ones_like(E),  # β₀
                E - E_mid,        # β₁
            ]
            if enable_quadratic:
                _f3_continuum_cols.append((E - E_mid) ** 2)  # β₂ (F-383)
            if continuum == "step_linear" and h_step_eff <= 0:
                _f3_continuum_cols.append(_smooth_step(E, E_step, sigma_step))
            elif continuum == "step_linear_multi" and h_step_eff <= 0:
                # F-392: каждый intense anchor — отдельный step term.
                for e_anchor, sigma_a in multi_step_anchors:
                    _f3_continuum_cols.append(_smooth_step(E, e_anchor, sigma_a))

            def _build_cols(basis_fn) -> List[np.ndarray]:
                """Полная матрица плана с заданным basis и текущими
                энергиями (центроиды по паспортным E_keV).
                AUDIT-F3: sigma_k / I_dec / continuum-колонки берутся
                из pre-compute (см. _f3_*_meta / _f3_continuum_cols выше);
                пересчитываем ТОЛЬКО peak-колонки через basis_fn."""
                cs: List[np.ndarray] = []
                for _glist in _f3_group_meta:
                    col = np.zeros_like(E)
                    for E0, sigma_k, I_dec in _glist:
                        col += I_dec * basis_fn(E, E0, sigma_k) * bin_w
                    cs.append(col)
                for E0, sigma_k in _f3_indep_meta:
                    cs.append(basis_fn(E, E0, sigma_k) * bin_w)
                cs.extend(_f3_continuum_cols)
                return cs

            def _model_at(p_vec: np.ndarray) -> np.ndarray:
                lin = p_vec[:n_params]
                sigma_scale = float(p_vec[n_params + 0])
                T_g = (float(p_vec[n_params + 1])
                       if free_T_global else tail_param)
                basis_fn = _basis_for(sigma_scale, T_g)
                cs = _build_cols(basis_fn)
                Xn = np.column_stack(cs)
                return Xn @ lin

            def _residuals_nl(p_vec: np.ndarray) -> np.ndarray:
                m = _model_at(p_vec)
                return (y - m) / sigma_y

            res_nl = least_squares(
                _residuals_nl, x0, bounds=(lb_arr, ub_arr),
                method="trf", max_nfev=int(nonlinear_max_nfev),
                xtol=1e-10, ftol=1e-10, gtol=1e-10,
            )
            chi2_nl = float(np.sum(res_nl.fun ** 2))
            n_dof_nl = max(1, n_pts - len(x0))
            chi2_per_dof_nl = chi2_nl / n_dof_nl

            # BUG-3 Fix #3 (2026-06-02): tighter σ-scale acceptance window
            # [0.8, 1.2] (было [0.7, 1.3]).
            if (chi2_per_dof_nl < 0.95 * chi2_per_dof
                    and 0.8 <= float(res_nl.x[n_params + 0]) <= 1.2):
                params = res_nl.x[:n_params]
                sigma_scale_used = float(res_nl.x[n_params + 0])
                T_g_used = (float(res_nl.x[n_params + 1])
                            if free_T_global else tail_param)
                # Подменяем _basis для downstream рендера на refined версию.
                _basis = _basis_for(sigma_scale_used, T_g_used)
                cols = _build_cols(_basis)
                X = np.column_stack(cols)
                Xw = X * w[:, None]
                model = X @ params
                residuals = y - model
                chi2 = chi2_nl
                chi2_per_dof = chi2_per_dof_nl
                n_dof = n_dof_nl
                # AUDIT-F5 (2026-06-25): см. комментарий выше у первой cov-ветки.
                try:
                    _U, _s, _Vt = np.linalg.svd(Xw, full_matrices=False)
                    if _s.size == 0 or _s[-1] <= 0:
                        cov = None
                    else:
                        cov = (_Vt.T * (1.0 / (_s * _s))) @ _Vt * (chi2 / n_dof)
                except np.linalg.LinAlgError:
                    cov = None
                method_label = "lsq_linear+nl_refine"
                nonlinear_used = True
                nonlinear_message = (
                    f"F-126 nonlinear refine: σ_scale={sigma_scale_used:.3f}"
                    + (f", T={T_g_used:.2f}" if free_T_global else "")
                    + f", χ²/ν={chi2_per_dof_nl:.2f}"
                )
            else:
                nonlinear_message = (
                    "F-126 nonlinear refine отброшен "
                    f"(χ²/ν_lin={chi2_per_dof:.2f}, "
                    f"χ²/ν_nl={chi2_per_dof_nl:.2f})"
                )
        except ImportError:
            nonlinear_message = "F-126: scipy.optimize.least_squares недоступен"
        except Exception as exc:
            nonlinear_message = f"F-126 refinement не удался: {exc!r}"

    # ─── F-145 / v1.17.8 — Phase A free-centroid side-fit ───────────
    # Запускается ТОЛЬКО после основной locked-passport подгонки.
    # Сами params / model / continuum НЕ изменяются — это диагностический
    # side-fit для self-calibration. Возвращает в результате:
    #   • centroid_shifts_keV[k] — фитированный сдвиг центра k-й компоненты
    #   • phase_A_chi2_per_dof — χ²/ν Phase A (для сравнения с locked)
    #   • phase_A_converged — bool
    # Outer pipeline использует эти данные для refit'a E(N) калибровки и
    # последующего Phase D (locked-passport на новой шкале).
    pA_shifts: List[float] = [0.0] * len(components)
    pA_chi2_per_dof: Optional[float] = None
    pA_converged: bool = False
    if free_centroids and use_peak_image:
        try:
            from scipy.optimize import least_squares as _ls_pA  # type: ignore
            from scipy.optimize import lsq_linear as _lsq_pA   # type: ignore
            n_comp = len(components)
            dE_lb_arr = np.array(
                [-float(centroid_window_frac) * float(fwhm_at(c.E_keV))
                 for c in components], dtype=np.float64)
            dE_ub_arr = np.array(
                [+float(centroid_window_frac) * float(fwhm_at(c.E_keV))
                 for c in components], dtype=np.float64)
            free_T_pA = (tail_T_at is None) and use_peak_image and tail_param > 0
            # BUG-3 Fix #3 (2026-06-02): tighter σ-scale prior [0.8, 1.2].
            extra_lb = [0.8] + ([0.3] if free_T_pA else [])
            extra_ub = [1.2] + ([1.5] if free_T_pA else [])
            extra_x0 = [1.0] + ([float(tail_param)] if free_T_pA else [])
            n_extra_pA = len(extra_x0)
            x0_pA = np.concatenate([
                np.array(extra_x0, dtype=np.float64),
                np.zeros(n_comp, dtype=np.float64),
            ])
            lb_pA = np.concatenate([
                np.array(extra_lb, dtype=np.float64), dE_lb_arr])
            ub_pA = np.concatenate([
                np.array(extra_ub, dtype=np.float64), dE_ub_arr])

            def _basis_with_shift_pA(E_arr, E0_eff, sk, T_loc):
                if h_step_eff > 0 and T_loc > 0:
                    return _peak_image_with_step(
                        E_arr, E0_eff, sk, T_loc, h_step_eff)
                if h_step_eff > 0 and T_loc <= 0:
                    return _peak_image_with_step(
                        E_arr, E0_eff, sk, 0.0, h_step_eff)
                if T_loc > 0:
                    return _peak_image_unit_area(E_arr, E0_eff, sk, T_loc)
                return _gaussian_unit_area(E_arr, E0_eff, sk)

            def _build_cols_pA(sigma_scale_pA, T_g_pA, dE_arr):
                cs_local: List[np.ndarray] = []
                for g in groups:
                    col = np.zeros_like(E)
                    for k_idx in group_members[g]:
                        c = components[k_idx]
                        E0_eff = float(c.E_keV) + float(dE_arr[k_idx])
                        sigma_k = float(fwhm_at(c.E_keV)) / 2.355
                        sk_eff = max(1e-6, sigma_scale_pA * sigma_k)
                        T_loc = (float(tail_T_at(c.E_keV))
                                 if tail_T_at is not None else T_g_pA)
                        I_dec = float(c.I_gamma_pct) / 100.0
                        col += I_dec * _basis_with_shift_pA(
                            E, E0_eff, sk_eff, T_loc) * bin_w
                    cs_local.append(col)
                for i_idx, c in independents:
                    E0_eff = float(c.E_keV) + float(dE_arr[i_idx])
                    sigma_k = float(fwhm_at(c.E_keV)) / 2.355
                    sk_eff = max(1e-6, sigma_scale_pA * sigma_k)
                    T_loc = (float(tail_T_at(c.E_keV))
                             if tail_T_at is not None else T_g_pA)
                    cs_local.append(_basis_with_shift_pA(
                        E, E0_eff, sk_eff, T_loc) * bin_w)
                cs_local.append(np.ones_like(E))
                cs_local.append(E - E_mid)
                if enable_quadratic:
                    cs_local.append((E - E_mid) ** 2)        # β₂ (F-383)
                if continuum == "step_linear" and h_step_eff <= 0:
                    mean_shift = float(np.mean(dE_arr))
                    cs_local.append(_smooth_step(
                        E, E_step + mean_shift, sigma_step))
                elif continuum == "step_linear_multi" and h_step_eff <= 0:
                    # F-392: каждый anchor сдвигается на собственный dE_arr
                    # для своей компоненты; phantom anchors используют
                    # mean shift (анкер не входит в components list).
                    mean_shift = float(np.mean(dE_arr))
                    for e_anchor, sigma_a in multi_step_anchors:
                        # Найти компоненту с E_keV ≈ e_anchor чтобы взять
                        # её индивидуальный shift; fallback на mean_shift.
                        local_shift = mean_shift
                        for k_idx, c in enumerate(components):
                            if abs(float(c.E_keV) - e_anchor) < 1e-6:
                                local_shift = float(dE_arr[k_idx])
                                break
                        cs_local.append(_smooth_step(
                            E, e_anchor + local_shift, sigma_a))
                return cs_local

            def _residuals_pA(p_vec):
                sigma_scale_pA = float(p_vec[0])
                T_g_pA = (float(p_vec[1]) if free_T_pA else float(tail_param))
                dE_arr = p_vec[n_extra_pA:]
                cs_l = _build_cols_pA(sigma_scale_pA, T_g_pA, dE_arr)
                X_l = np.column_stack(cs_l)
                Xw_l = X_l * w[:, None]
                # Solve linear amplitudes at this configuration with same bounds
                try:
                    res_lin = _lsq_pA(
                        Xw_l, yw, bounds=(lb, ub),
                        method="trf", tol=1e-8, max_iter=200)
                    lin = res_lin.x
                except Exception:
                    lin, *_ = np.linalg.lstsq(Xw_l, yw, rcond=None)
                model_pA = X_l @ lin
                return (y - model_pA) / sigma_y

            res_pA = _ls_pA(
                _residuals_pA, x0_pA, bounds=(lb_pA, ub_pA),
                method="trf", max_nfev=300,
                xtol=1e-8, ftol=1e-8, gtol=1e-8,
            )
            chi2_pA = float(np.sum(res_pA.fun ** 2))
            n_dof_pA = max(1, n_pts - len(x0_pA) - n_params)
            pA_chi2_per_dof = chi2_pA / n_dof_pA
            pA_converged = bool(res_pA.success)
            pA_shifts = [float(v) for v in res_pA.x[n_extra_pA:]]
        except ImportError:
            pA_chi2_per_dof = None
        except Exception:
            pA_chi2_per_dof = None

    # ─── распаковка площадей по компонентам и активностей групп ──────
    a_nuc: dict = {}
    areas_by_comp: List[float] = [0.0] * len(components)
    sigma_by_comp: List[float] = [0.0] * len(components)

    # 1. группы → площади a_k = A_nuc · I_k/100
    for j, kind in enumerate(col_kinds):
        if kind != "group":
            continue
        g = col_keys[j]
        A_g = float(max(0.0, params[j]))
        sA = float(math.sqrt(max(0.0, cov[j, j]))) if cov is not None else float("nan")
        a_nuc[g] = (A_g, sA)
        for k in group_members[g]:
            I_dec = float(components[k].I_gamma_pct) / 100.0
            areas_by_comp[k] = A_g * I_dec
            sigma_by_comp[k] = sA * I_dec

    # 2. независимые → собственная площадь
    indep_iter = iter(independents)
    for j, kind in enumerate(col_kinds):
        if kind != "indep":
            continue
        try:
            i, c = next(indep_iter)
        except StopIteration:
            break
        a_i = float(max(0.0, params[j]))
        s_i = float(math.sqrt(max(0.0, cov[j, j]))) if cov is not None else float("nan")
        areas_by_comp[i] = a_i
        sigma_by_comp[i] = s_i

    # ─── континуум и модельные кривые для рендера ─────────────────────
    # Континуум — только параметры континуума, без вкладов компонент.
    cont = np.zeros_like(E)
    cont_params_list: List[float] = []
    for j, kind in enumerate(col_kinds):
        if kind in ("beta0", "beta1", "beta2", "step"):
            cont += cols[j] * params[j]
            cont_params_list.append(float(params[j]))

    # F-373 — clamp continuum to ≥0 (counts can't be negative). On long
    # ROIs with a negative β₁ slope, a pure linear continuum can dive
    # below zero and the renderer then showed a fit that went into the
    # negatives. NNLS bounds didn't prevent this because β₀ and β₁ are
    # unconstrained — only β_step has lb=0.
    #
    # BUG-3 Fix #4 (2026-06-02) — hard continuum-non-negative constraint
    # diagnostics. Если хотя бы один channel дал negative continuum (до
    # clamp), это сигнал что fit натянул β₀/β₁ так чтобы compensate за
    # strong-line counts которые fit не смог распределить по компонентам
    # (например: strong line «съедена» в baseline). Эмитим RuntimeWarning
    # и помечаем converged=False — downstream может (1) показать badge,
    # (2) исключить такой fit из activity calculation, (3) trigger
    # diagnostic re-фит. Сам clamp оставляем (renderer-safety).
    cont_pre_clamp_min = float(np.min(cont)) if cont.size else 0.0
    if cont_pre_clamp_min < 0.0:
        import warnings as _warnings_bug3
        _warnings_bug3.warn(
            f"BUG-3 Fix #4: continuum went negative "
            f"(min={cont_pre_clamp_min:.3g}) — fit pathology suspected. "
            f"Marking converged=False.",
            RuntimeWarning,
            stacklevel=2,
        )
        converged = False
    cont = np.maximum(cont, 0.0)

    # per-component «g_plus_cont» (континуум + только данная компонента)
    # Самплируем в counts/channel: density × bin_width.
    comp_g_plus_cont: List[List[float]] = []
    comp_g_base: List[List[float]] = []
    for k, c in enumerate(components):
        sigma_k = float(fwhm_at(c.E_keV)) / 2.355
        gauss_k = _basis(E, float(c.E_keV), sigma_k) * bin_w
        contrib_k = gauss_k * areas_by_comp[k]
        comp_g_plus_cont.append([float(v) for v in (cont + contrib_k)])
        comp_g_base.append([float(v) for v in cont])

    total = cont.copy()
    for k in range(len(components)):
        sigma_k = float(fwhm_at(components[k].E_keV)) / 2.355
        gauss_k = _basis(E, float(components[k].E_keV), sigma_k) * bin_w
        total += gauss_k * areas_by_comp[k]

    # closure %: 100 · (Σ model − Σ data) / Σ data
    sum_data = float(np.sum(y))
    sum_model = float(np.sum(total))
    closure_pct = 100.0 * (sum_model - sum_data) / max(sum_data, 1.0)

    # Сборка ComponentFit
    comp_fits: List[ComponentFit] = []
    for k, c in enumerate(components):
        comp_fits.append(ComponentFit(
            nuclide=c.nuclide,
            E_keV=float(c.E_keV),
            I_pct=float(c.I_gamma_pct),
            area=float(areas_by_comp[k]),
            sigma_area=float(sigma_by_comp[k]),
            group=c.group,
        ))

    return CoupledFitResult(
        id=cluster_id,
        title=title,
        roi_low_ch=int(roi_low_ch),
        roi_high_ch=int(roi_low_ch + n_pts),
        n_channels=int(n_pts),
        E_keV=[float(v) for v in E],
        data=[float(v) for v in y],
        continuum=[float(v) for v in cont],
        total=[float(v) for v in total],
        components=comp_fits,
        component_g_plus_cont=comp_g_plus_cont,
        component_g_base=comp_g_base,
        continuum_model=continuum,
        continuum_params=cont_params_list,
        chi2_per_dof=float(chi2_per_dof),
        n_dof=int(n_dof),
        closure_pct=float(closure_pct),
        a_nuclide=a_nuc,
        converged=bool(converged),
        method=method_label,
        notes=(
            f"basis={basis_label}"
            + (
                f"; F-392 multi-step anchors={len(multi_step_anchors)}"
                if continuum == "step_linear_multi" else ""
            )
            + (f"; {nonlinear_message}" if nonlinear_message else "")
            + (
                "; F-145 Phase A: "
                f"χ²/ν={pA_chi2_per_dof:.2f}, max|dE|="
                f"{max((abs(v) for v in pA_shifts), default=0.0):.3f} кэВ"
                if pA_chi2_per_dof is not None else ""
            )
        ),
        centroid_shifts_keV=pA_shifts,
        phase_A_chi2_per_dof=pA_chi2_per_dof,
        phase_A_converged=bool(pA_converged),
        # F-392.1 / v1.18.29 — surface для downstream JSON. Anchors как
        # list of (E_keV, σ_step); threshold_pct — None при non-multi
        # continuum (чтобы JSON показывал null вместо misleading 4.0).
        multi_step_anchors=list(multi_step_anchors),
        multi_step_intensity_threshold_pct=(
            float(multi_step_intensity_threshold_pct)
            if continuum == "step_linear_multi" else None
        ),
    )


__all__ = [
    "ComponentSpec", "ComponentFit", "CoupledFitResult",
    "coupled_intensity_fit",
    # F-120
    "T_TAIL_DEFAULT_NAI",
    # F-127 / v1.17.7
    "nai_tail_T_at",
    "NAI_T_E_REF_KEV", "NAI_T_E_T_REF", "NAI_T_E_SLOPE",
    "NAI_T_E_T_MIN", "NAI_T_E_T_MAX",
    # F-133 / v1.17.7 — per-line step (ГОСТ форма пика)
    "H_STEP_DEFAULT_NAI",
    "_peak_image_with_step",
]
