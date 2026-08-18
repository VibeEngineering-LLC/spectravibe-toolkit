"""
F-346 / v1.18.23.0 — ГОСТ 26874-86 § 3.3.2-3.3.4 канонические методы.

Реализует три **формальных** ГОСТ-метода определения центроиды пика
и два метода вычитания фонового пьедестала под пиком, как
ОПЦИОНАЛЬНУЮ альтернативу нашему основному peak_image (Gauss+tail+step)
fit'у. Полезно когда:

  * нужен **независимый** sanity-check центроиды из peak_image
    (особенно после F-145 self-calibration);
  * выпускается формальный отчёт о соответствии ГОСТ-методу;
  * в diagnostics для образовательных walkthrough-отчётов
    (Technical PDF F-159 Шаг 5/6) нужно явно показать ГОСТ-формулы.

Источник: ГОСТ 26874-86 «Спектрометры энергий ионизирующих излучений.
Методы измерения основных параметров» — § 3.3.2, § 3.3.3, § 3.3.4.
Layer 2: [24, §3.3.3] и [24, §3.3.4]. Layer 1: [GOST-26874-3],
[GOST-26874-3.3.2].

Контракт F-346:
  1. Pedestal subtraction:
       - symmetric: горизонтальная линия через средние Nl, Nh слева/справа
       - asymmetric: 9-канальное усреднение на расстоянии ≥4·FWHM,
         линейная интерполяция под пиком (формула 4)
       - threshold: фон < 2% от N_max → не вычитать (§ 3.3.2 первый абзац)

  2. Centroid determination (после вычитания фона):
       - graphical (§3.3.3.1): пересечения с огибающей на 1/2 и 3/4 высоты
       - weighted_mean (§3.3.3.2, формула 5): n_c = Σ N_i·n_i / Σ N_i
         только по точкам ВЫШЕ полувысоты
       - graphoanalytic (§3.3.3.3, формулы 6-10): линейная регрессия
         ln(N_i / N_{i+1}) = A·n_i − B → n_c = B/A + 1/2

  3. FWHM (§3.3.4):
       - linear interpolation (§3.3.4.2, формула 12): соседние точки
         на полувысоте слева и справа
       - graphoanalytic (§3.3.4.3, формула 13): Δn = 2·√(ln2/A)
         ⚠️ ВНИМАНИЕ: формула в ГОСТ даёт FWHM/√2 от истинного значения
         (см. F-347 sanity-check); используем **исправленную форму**
         Δn = 2·√(2·ln2/A) = 2σ·√(2·ln2), но возвращаем ОБА значения
         для прозрачности.

  4. ROI selection: peak_channel ± k·FWHM, k=0.5 для интегрирования,
     pedestal-окна на расстоянии ≥4·FWHM от максимума.

Все функции read-only по отношению к входному спектру — возвращают
diagnostic dict + новый массив counts с вычтенным фоном (или None,
если фон не нужен).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Sequence, Tuple

import numpy as np


# ──────────────────────────────────────────────────────────────────
# Константы
# ──────────────────────────────────────────────────────────────────

# Стандартный фактор связи σ ↔ FWHM для Гаусса: FWHM = σ · 2√(2·ln2).
_SIGMA_TO_FWHM = 2.0 * math.sqrt(2.0 * math.log(2.0))  # ≈ 2.354820045

# Порог вычитания пьедестала по § 3.3.2: при N_фон < 2% от N_max
# пьедестал НЕ вычитают.
PEDESTAL_THRESHOLD_FRAC = 0.02

# Множитель FWHM для расстояния от максимума до 9-канального окна
# (§ 3.3.2.2: «не менее четырёхкратного значения ширины пика»).
DEFAULT_PEDESTAL_GAP_FWHM = 4.0

# Ширина окна усреднения пьедестала (§ 3.3.2.2: «9 каналов»).
DEFAULT_PEDESTAL_WINDOW_CHANNELS = 9


# ──────────────────────────────────────────────────────────────────
# Result dataclasses
# ──────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class PedestalResult:
    """Результат вычитания фонового пьедестала под пиком."""
    method: str                          # 'symmetric' | 'asymmetric' | 'none'
    pedestal: np.ndarray                 # массив фона на ROI (counts)
    counts_net: np.ndarray               # counts − pedestal на ROI
    roi_lo: int
    roi_hi: int                          # exclusive
    N_l_avg: Optional[float] = None      # среднее по 9-кан. окну слева
    N_h_avg: Optional[float] = None      # ... справа
    pedestal_gap_channels: Optional[int] = None  # отступ от max до окна
    skipped_reason: str = ""             # если method='none', почему


@dataclass(frozen=True)
class CentroidResult:
    """Результат определения центроиды одним из ГОСТ-методов."""
    method: str                          # 'graphical' | 'weighted_mean' | 'graphoanalytic'
    n_c: float                           # центроида (канальный номер)
    sigma_n_c: Optional[float] = None    # σ_центроиды если оценимо
    A: Optional[float] = None            # параметр Гаусс ln-fit'а (1/σ²); только для graphoanalytic
    B: Optional[float] = None            # сдвиг; только для graphoanalytic
    n_points_used: int = 0
    notes: str = ""


@dataclass(frozen=True)
class FwhmResult:
    """Результат определения FWHM одним из ГОСТ-методов."""
    method: str                          # 'linear_interp' | 'graphoanalytic' | 'graphoanalytic_corrected'
    fwhm_channels: float
    # Для graphoanalytic: A из ln-fit (см. § 3.3.3.3 формула 7).
    A: Optional[float] = None
    # Сырое ГОСТ-значение (формула 13). Может отличаться от corrected на √2.
    raw_gost_fwhm: Optional[float] = None
    # Математически корректная FWHM Гаусса: 2·√(2·ln2/A).
    corrected_gauss_fwhm: Optional[float] = None
    n_max: Optional[float] = None
    half_max: Optional[float] = None
    n_left: Optional[float] = None
    n_right: Optional[float] = None
    notes: str = ""


# ──────────────────────────────────────────────────────────────────
# § 3.3.2 — Pedestal (фоновый пьедестал)
# ──────────────────────────────────────────────────────────────────

def should_subtract_pedestal(
    counts: Sequence[float],
    peak_channel: int,
    fwhm_channels: Optional[float] = None,
    *,
    bg_gap_fwhm: float = DEFAULT_PEDESTAL_GAP_FWHM,
    bg_window: int = DEFAULT_PEDESTAL_WINDOW_CHANNELS,
) -> Tuple[bool, float, float]:
    """§ 3.3.2 первый абзац: вычитание не проводят при N_фон < 2% от N_max.

    Оценка фона: 9-канальное усреднение на расстоянии ≥4·FWHM от
    максимума пика (как в § 3.3.2.2). Если fwhm_channels не задан —
    используется дефолтный gap=20 каналов (грубая оценка).

    Parameters
    ----------
    counts : sequence
        Полный спектр (counts на канал).
    peak_channel : int
        Положение максимума пика.
    fwhm_channels : float, optional
        Оценочная FWHM пика в каналах. Если None — gap=20.
    bg_gap_fwhm : float, default 4.0
        Множитель FWHM для отступа от пика до фон-окна.
    bg_window : int, default 9
        Ширина окна усреднения (§ 3.3.2.2).

    Returns
    -------
    (do_subtract, n_max, n_bg_est)
        do_subtract — True если N_bg/N_max ≥ 0.02 (§ 3.3.2 порог)
        n_max — значение в пике
        n_bg_est — оценочный фон (среднее двух 9-кан. окон)
    """
    c = np.asarray(counts, dtype=float)
    if not (0 <= peak_channel < len(c)):
        raise ValueError(f"peak_channel {peak_channel} out of range [0, {len(c)})")
    n_max = float(c[peak_channel])
    if n_max <= 0:
        return False, n_max, 0.0

    if fwhm_channels and fwhm_channels > 0:
        gap = max(1, int(round(bg_gap_fwhm * fwhm_channels)))
    else:
        gap = 20

    half_w = bg_window // 2
    l_center = peak_channel - gap
    h_center = peak_channel + gap

    if l_center - half_w < 0:
        l_center = half_w
    if h_center + half_w >= len(c):
        h_center = len(c) - half_w - 1

    if l_center + half_w >= peak_channel or h_center - half_w <= peak_channel:
        # Окно вылазит за пик — fallback на края
        edges = np.concatenate([c[:max(1, peak_channel - 3)], c[peak_channel + 4:]])
        n_bg = float(np.mean(edges)) if len(edges) else 0.0
    else:
        S_l = float(np.sum(c[l_center - half_w:l_center + half_w + 1]))
        S_h = float(np.sum(c[h_center - half_w:h_center + half_w + 1]))
        n_bg = 0.5 * (S_l + S_h) / bg_window

    return (n_bg / n_max >= PEDESTAL_THRESHOLD_FRAC), n_max, n_bg


def gost_pedestal_symmetric(
    counts: Sequence[float],
    peak_channel: int,
    fwhm_channels: float,
    *,
    roi_half_fwhm: float = 2.5,
    pedestal_gap_fwhm: float = DEFAULT_PEDESTAL_GAP_FWHM,
    pedestal_window: int = DEFAULT_PEDESTAL_WINDOW_CHANNELS,
) -> PedestalResult:
    """§ 3.3.2.1 — симметричный пьедестал (Черт. 3, а).

    Через средние N̅_l и N̅_h слева/справа от пика на расстоянии
    ≥ pedestal_gap_fwhm·FWHM проводят ГОРИЗОНТАЛЬНУЮ линию (среднее
    из обеих сторон). Используется когда |N̅_l − N̅_h| мало по сравнению
    со статистической погрешностью σ_stat.

    Возвращает PedestalResult с уже вычисленным counts_net на ROI.
    """
    return _gost_pedestal_impl(
        counts, peak_channel, fwhm_channels,
        roi_half_fwhm=roi_half_fwhm,
        pedestal_gap_fwhm=pedestal_gap_fwhm,
        pedestal_window=pedestal_window,
        force_method="symmetric",
    )


def gost_pedestal_asymmetric(
    counts: Sequence[float],
    peak_channel: int,
    fwhm_channels: float,
    *,
    roi_half_fwhm: float = 2.5,
    pedestal_gap_fwhm: float = DEFAULT_PEDESTAL_GAP_FWHM,
    pedestal_window: int = DEFAULT_PEDESTAL_WINDOW_CHANNELS,
) -> PedestalResult:
    """§ 3.3.2.2 — асимметричный пьедестал (Черт. 3, б) + формула (4).

    На расстоянии ≥ pedestal_gap_fwhm·FWHM от максимума выбирают
    участки по pedestal_window каналов; вычисляют:
        S_l = Σ_{i=l-4..l+4} N_i,   N̅_l = S_l / 9
        S_h = Σ_{i=h-4..h+4} N_i,   N̅_h = S_h / 9
    Затем линейная интерполяция в каждом канале i под пиком:
        N_li = N̅_l + (N̅_h − N̅_l)·(i − l)/(h − l)   ... (формула 4)
    """
    return _gost_pedestal_impl(
        counts, peak_channel, fwhm_channels,
        roi_half_fwhm=roi_half_fwhm,
        pedestal_gap_fwhm=pedestal_gap_fwhm,
        pedestal_window=pedestal_window,
        force_method="asymmetric",
    )


def gost_select_pedestal_method(
    counts: Sequence[float],
    peak_channel: int,
    fwhm_channels: float,
    *,
    symmetry_tol_sigma: float = 1.0,
    **kwargs,
) -> PedestalResult:
    """§ 3.3.2 — авто-выбор симметричного / асимметричного метода.

    Сначала считаем N̅_l, N̅_h по 9-канальному усреднению. Если
    |N̅_l − N̅_h| ≤ symmetry_tol_sigma · max(σ_l, σ_h) (статистически
    одинаковы) → используем симметричный (горизонтальную линию из их
    среднего). Иначе — асимметричный (формула 4).

    Также проверяем порог 2% (см. should_subtract_pedestal). Если фон
    меньше — возвращаем PedestalResult с method='none'.
    """
    do_sub, n_max, n_bg = should_subtract_pedestal(
        counts, peak_channel, fwhm_channels,
        bg_gap_fwhm=kwargs.get("pedestal_gap_fwhm", DEFAULT_PEDESTAL_GAP_FWHM),
        bg_window=kwargs.get("pedestal_window", DEFAULT_PEDESTAL_WINDOW_CHANNELS),
    )
    if not do_sub:
        # Возвращаем «пустую» подложку — counts_net = counts на ROI.
        roi_lo, roi_hi = _roi_around_peak(
            len(counts), peak_channel, fwhm_channels,
            kwargs.get("roi_half_fwhm", 2.5),
        )
        c = np.asarray(counts, dtype=float)
        zero_ped = np.zeros(roi_hi - roi_lo)
        return PedestalResult(
            method="none",
            pedestal=zero_ped,
            counts_net=c[roi_lo:roi_hi].copy(),
            roi_lo=roi_lo, roi_hi=roi_hi,
            N_l_avg=n_bg, N_h_avg=n_bg,
            skipped_reason=(
                f"N_bg/N_max = {n_bg/n_max if n_max > 0 else 0:.3f} "
                f"< {PEDESTAL_THRESHOLD_FRAC} (§ 3.3.2)"
            ),
        )

    # Сначала пробуем асимметричный (он всегда корректен), потом
    # проверяем симметрию N̅_l ↔ N̅_h.
    asym = _gost_pedestal_impl(
        counts, peak_channel, fwhm_channels,
        force_method="asymmetric", **kwargs,
    )
    if asym.N_l_avg is None or asym.N_h_avg is None:
        return asym

    # Статистическая погрешность среднего ≈ √(N̅ / window).
    window = kwargs.get("pedestal_window", DEFAULT_PEDESTAL_WINDOW_CHANNELS)
    sigma_l = math.sqrt(max(asym.N_l_avg, 1.0) / window)
    sigma_h = math.sqrt(max(asym.N_h_avg, 1.0) / window)
    diff = abs(asym.N_l_avg - asym.N_h_avg)
    if diff <= symmetry_tol_sigma * max(sigma_l, sigma_h):
        return _gost_pedestal_impl(
            counts, peak_channel, fwhm_channels,
            force_method="symmetric", **kwargs,
        )
    return asym


def _gost_pedestal_impl(
    counts: Sequence[float],
    peak_channel: int,
    fwhm_channels: float,
    *,
    roi_half_fwhm: float = 2.5,
    pedestal_gap_fwhm: float = DEFAULT_PEDESTAL_GAP_FWHM,
    pedestal_window: int = DEFAULT_PEDESTAL_WINDOW_CHANNELS,
    force_method: str = "asymmetric",
) -> PedestalResult:
    """Общая реализация § 3.3.2.1 + § 3.3.2.2."""
    c = np.asarray(counts, dtype=float)
    n_ch = len(c)
    if not (0 <= peak_channel < n_ch):
        raise ValueError(f"peak_channel {peak_channel} out of range")
    if fwhm_channels <= 0:
        raise ValueError(f"fwhm_channels must be > 0, got {fwhm_channels}")
    if pedestal_window < 3:
        raise ValueError(f"pedestal_window must be ≥ 3, got {pedestal_window}")

    # ROI пика для возврата counts_net.
    roi_lo, roi_hi = _roi_around_peak(n_ch, peak_channel, fwhm_channels, roi_half_fwhm)

    # Центры 9-канальных окон слева/справа.
    gap_channels = max(1, int(round(pedestal_gap_fwhm * fwhm_channels)))
    half_window = pedestal_window // 2

    l_center = peak_channel - gap_channels
    h_center = peak_channel + gap_channels

    # Сжатие если окна вылазят за края спектра.
    if l_center - half_window < 0:
        l_center = half_window
    if h_center + half_window >= n_ch:
        h_center = n_ch - half_window - 1

    # Если после сжатия окна перекрывают пик — pedestal не вычитаем.
    if l_center + half_window >= peak_channel or h_center - half_window <= peak_channel:
        zero_ped = np.zeros(roi_hi - roi_lo)
        return PedestalResult(
            method="none",
            pedestal=zero_ped,
            counts_net=c[roi_lo:roi_hi].copy(),
            roi_lo=roi_lo, roi_hi=roi_hi,
            skipped_reason="pedestal windows overlap peak (insufficient spectrum margin)",
        )

    S_l = float(np.sum(c[l_center - half_window:l_center + half_window + 1]))
    S_h = float(np.sum(c[h_center - half_window:h_center + half_window + 1]))
    N_l_avg = S_l / pedestal_window
    N_h_avg = S_h / pedestal_window

    pedestal = np.empty(roi_hi - roi_lo)
    if force_method == "symmetric":
        # § 3.3.2.1 — горизонтальная линия среднего двух пьедесталов.
        N_const = 0.5 * (N_l_avg + N_h_avg)
        pedestal[:] = N_const
    elif force_method == "asymmetric":
        # § 3.3.2.2 формула (4) — линейная интерполяция между центрами.
        slope = (N_h_avg - N_l_avg) / float(h_center - l_center)
        for k, i in enumerate(range(roi_lo, roi_hi)):
            pedestal[k] = N_l_avg + slope * (i - l_center)
    else:
        raise ValueError(f"unknown method: {force_method}")

    counts_net = np.maximum(0.0, c[roi_lo:roi_hi] - pedestal)
    return PedestalResult(
        method=force_method,
        pedestal=pedestal,
        counts_net=counts_net,
        roi_lo=roi_lo, roi_hi=roi_hi,
        N_l_avg=N_l_avg, N_h_avg=N_h_avg,
        pedestal_gap_channels=gap_channels,
    )


# ──────────────────────────────────────────────────────────────────
# § 3.3.3 — Centroid (центроида пика)
# ──────────────────────────────────────────────────────────────────

def gost_centroid_graphical(
    counts_net: Sequence[float],
    *,
    channel_offset: int = 0,
) -> CentroidResult:
    """§ 3.3.3.1 — графический метод (Черт. 4).

    На полувысоте проводим горизонтальную линию; через середину
    пересечения с огибающей пика — нормаль. Численная аппроксимация:
    находим n_left и n_right на полувысоте (линейная интерполяция),
    центроида = (n_left + n_right) / 2.

    channel_offset — добавляется к результату чтобы вернуть абсолютный
    номер канала (если counts_net — это ROI, не полный спектр).
    """
    c = np.asarray(counts_net, dtype=float)
    if len(c) < 3:
        raise ValueError("centroid: need ≥3 points in ROI")
    n_max = float(c.max())
    if n_max <= 0:
        raise ValueError("centroid: empty peak (max ≤ 0)")
    peak_idx = int(np.argmax(c))
    half = 0.5 * n_max

    # Левый склон: первая точка слева с counts < half.
    n_left = _find_crossing_left(c, peak_idx, half)
    n_right = _find_crossing_right(c, peak_idx, half)
    if n_left is None or n_right is None:
        raise ValueError(
            "centroid: failed to find half-max crossings "
            f"(left={n_left}, right={n_right})"
        )

    n_c = 0.5 * (n_left + n_right) + channel_offset
    return CentroidResult(
        method="graphical",
        n_c=float(n_c),
        n_points_used=int(np.sum(c >= half)),
        notes=f"left={n_left:.3f}, right={n_right:.3f}, half={half:.3f}",
    )


def gost_centroid_weighted_mean(
    counts_net: Sequence[float],
    *,
    channel_offset: int = 0,
) -> CentroidResult:
    """§ 3.3.3.2 формула (5) — средневзвешенное.

        n_c = Σ N_i · n_i / Σ N_i

    ВАЖНО: суммирование только по точкам ВЫШЕ полувысоты (симметричная
    часть пика по ГОСТ), чтобы исключить асимметричный хвост.
    """
    c = np.asarray(counts_net, dtype=float)
    if len(c) < 3:
        raise ValueError("centroid: need ≥3 points in ROI")
    n_max = float(c.max())
    if n_max <= 0:
        raise ValueError("centroid: empty peak (max ≤ 0)")

    half = 0.5 * n_max
    mask = c >= half
    if int(np.sum(mask)) < 3:
        raise ValueError(
            f"centroid weighted_mean: only {int(np.sum(mask))} points above half-max"
        )

    idx = np.arange(len(c))
    n_above = c[mask]
    i_above = idx[mask]
    sum_w = float(np.sum(n_above))
    sum_wi = float(np.sum(n_above * i_above))
    n_c = sum_wi / sum_w

    # σ_n_c ≈ sqrt(Σ N_i·(n_i − n_c)²) / Σ N_i ≈ σ_peak / √N_total
    var = float(np.sum(n_above * (i_above - n_c) ** 2) / sum_w)
    sigma_n_c = math.sqrt(var) / math.sqrt(max(sum_w, 1.0))

    return CentroidResult(
        method="weighted_mean",
        n_c=float(n_c + channel_offset),
        sigma_n_c=float(sigma_n_c),
        n_points_used=int(np.sum(mask)),
        notes=f"half={half:.3f}, n_above={int(np.sum(mask))}, sum_w={sum_w:.1f}",
    )


def gost_centroid_graphoanalytic(
    counts_net: Sequence[float],
    *,
    channel_offset: int = 0,
    min_count_floor: float = 1.0,
) -> CentroidResult:
    """§ 3.3.3.3 — графоаналитический (Черт. 5, формулы 6-10).

    Предполагаем Гауссиан N_i = N_max·exp[-(n_0 − n_i)²/(2σ²)].
    Тогда:
        ln(N_i / N_{i+1}) = A·n_i − B
    где
        A = 1/σ²
        B = (2·n_0 − 1) / (2·σ²)
    → n_c = B/A + 1/2

    Точки берём ТОЛЬКО выше полувысоты (как требует § 3.3.3.3).
    Малые значения N_i (< min_count_floor) исключаем, чтобы ln не
    взорвался от Пуассон-шума.

    Возвращает CentroidResult с A, B, σ_n_c из МНК.
    """
    c = np.asarray(counts_net, dtype=float)
    if len(c) < 4:
        raise ValueError("centroid graphoanalytic: need ≥4 points in ROI")
    n_max = float(c.max())
    if n_max <= 0:
        raise ValueError("centroid graphoanalytic: empty peak")

    half = 0.5 * n_max
    mask = (c >= half) & (c >= min_count_floor)
    # Также нужно чтобы N_{i+1} тоже было ≥ floor для ln-отношения.
    idx = np.arange(len(c))
    pairs = []
    for i in idx[mask]:
        if i + 1 >= len(c):
            continue
        if c[i + 1] < min_count_floor:
            continue
        pairs.append((float(i), math.log(c[i] / c[i + 1])))

    if len(pairs) < 3:
        raise ValueError(
            f"centroid graphoanalytic: only {len(pairs)} pairs above half-max"
        )

    x = np.array([p[0] for p in pairs])
    y = np.array([p[1] for p in pairs])

    # ГОСТ Прил. 2: МНК с равными весами.
    # y = A·x − B → линейная регрессия.
    m = len(x)
    sum_x = float(np.sum(x))
    sum_y = float(np.sum(y))
    sum_xx = float(np.sum(x * x))
    sum_xy = float(np.sum(x * y))
    denom = m * sum_xx - sum_x ** 2
    if abs(denom) < 1e-12:
        raise ValueError("centroid graphoanalytic: singular MNK system")

    A = (m * sum_xy - sum_x * sum_y) / denom
    minus_B = (sum_y * sum_xx - sum_x * sum_xy) / denom  # это (−B) в форме y = A·x + (−B)
    B = -minus_B

    if A <= 0:
        raise ValueError(
            f"centroid graphoanalytic: A={A:.6e} ≤ 0 — fit failed "
            "(peak shape not Gaussian above half-max)"
        )

    n_c = B / A + 0.5
    sigma_squared = 1.0 / A
    sigma_n_c = math.sqrt(sigma_squared / m)  # грубая оценка

    return CentroidResult(
        method="graphoanalytic",
        n_c=float(n_c + channel_offset),
        sigma_n_c=float(sigma_n_c),
        A=float(A),
        B=float(B),
        n_points_used=m,
        notes=f"σ²=1/A={sigma_squared:.3f} ch², σ={math.sqrt(sigma_squared):.3f} ch",
    )


# ──────────────────────────────────────────────────────────────────
# § 3.3.4 — FWHM
# ──────────────────────────────────────────────────────────────────

def gost_fwhm_linear_interp(
    counts_net: Sequence[float],
) -> FwhmResult:
    """§ 3.3.4.2 формула (12) — линейная интерполяция точек пересечения
    уровня полувысоты слева и справа от максимума.
    """
    c = np.asarray(counts_net, dtype=float)
    if len(c) < 3:
        raise ValueError("fwhm linear_interp: need ≥3 points")
    n_max = float(c.max())
    if n_max <= 0:
        raise ValueError("fwhm linear_interp: empty peak")
    peak_idx = int(np.argmax(c))
    half = 0.5 * n_max

    n_left = _find_crossing_left(c, peak_idx, half)
    n_right = _find_crossing_right(c, peak_idx, half)
    if n_left is None or n_right is None:
        raise ValueError(
            "fwhm linear_interp: half-max crossings not found "
            f"(left={n_left}, right={n_right})"
        )

    fwhm = float(n_right - n_left)
    return FwhmResult(
        method="linear_interp",
        fwhm_channels=fwhm,
        n_max=n_max,
        half_max=half,
        n_left=float(n_left),
        n_right=float(n_right),
        notes=f"peak_idx={peak_idx}",
    )


def gost_fwhm_graphoanalytic(
    A: float,
    *,
    use_corrected_formula: bool = True,
) -> FwhmResult:
    """§ 3.3.4.3 формула (13) — FWHM из параметра A графоаналитического
    fit'а (§ 3.3.3.3).

    ГОСТ формула (13):   Δn = 2·√(ln2 / A)        — даёт FWHM/√2 (ошибка!)
    Стандарт Gauss:      Δn = 2·√(2·ln2 / A)      — настоящая FWHM

    Расхождение происходит из того, что для Гауссиана
    f(x) = N_max·exp[-(x−x₀)²/(2σ²)] FWHM связана с σ через
    FWHM = 2σ·√(2·ln2), а не 2σ·√(ln2). Параметр A = 1/σ², поэтому
    корректная формула — 2·√(2·ln2/A).

    Возвращаем ОБА значения:
      - raw_gost_fwhm = ГОСТ-формула как есть (для формального отчёта)
      - corrected_gauss_fwhm = математически правильная FWHM Гаусса

    Если use_corrected_formula=True (по умолчанию) — `fwhm_channels`
    возвращает corrected (рекомендуется); False — сырое ГОСТ-значение.
    """
    if A <= 0:
        raise ValueError(f"fwhm graphoanalytic: A={A:.6e} ≤ 0 — invalid")

    raw_gost = 2.0 * math.sqrt(math.log(2.0) / A)
    corrected = 2.0 * math.sqrt(2.0 * math.log(2.0) / A)

    fwhm = corrected if use_corrected_formula else raw_gost
    method = "graphoanalytic_corrected" if use_corrected_formula else "graphoanalytic_gost_raw"
    return FwhmResult(
        method=method,
        fwhm_channels=float(fwhm),
        A=float(A),
        raw_gost_fwhm=float(raw_gost),
        corrected_gauss_fwhm=float(corrected),
        notes=(
            f"ratio corrected/raw_gost = {corrected/raw_gost:.4f} (= √2 ≈ 1.4142). "
            f"ГОСТ формула (13) даёт FWHM/√2 — вероятно опечатка в стандарте; "
            f"используем corrected по умолчанию."
        ),
    )


# ──────────────────────────────────────────────────────────────────
# ROI helper
# ──────────────────────────────────────────────────────────────────

def gost_roi_from_fwhm(
    peak_channel: int,
    fwhm_channels: float,
    *,
    half_fwhm: float = 2.5,
    n_total: Optional[int] = None,
) -> Tuple[int, int]:
    """ROI пика ± half_fwhm·FWHM (полуоткрытый интервал [lo, hi))."""
    return _roi_around_peak(
        n_total if n_total is not None else (peak_channel + 1000),
        peak_channel, fwhm_channels, half_fwhm,
    )


def _roi_around_peak(
    n_total: int,
    peak_channel: int,
    fwhm_channels: float,
    half_fwhm: float,
) -> Tuple[int, int]:
    half = max(2, int(round(half_fwhm * fwhm_channels)))
    lo = max(0, peak_channel - half)
    hi = min(n_total, peak_channel + half + 1)
    return lo, hi


def _find_crossing_left(c: np.ndarray, peak_idx: int, level: float) -> Optional[float]:
    """Линейная интерполяция пересечения c[i] = level слева от peak_idx."""
    for i in range(peak_idx - 1, -1, -1):
        if c[i] < level <= c[i + 1]:
            denom = c[i + 1] - c[i]
            if abs(denom) < 1e-12:
                return float(i + 0.5)
            return float(i + (level - c[i]) / denom)
    return None


def _find_crossing_right(c: np.ndarray, peak_idx: int, level: float) -> Optional[float]:
    """Линейная интерполяция пересечения c[i] = level справа от peak_idx."""
    for i in range(peak_idx, len(c) - 1):
        if c[i] >= level > c[i + 1]:
            denom = c[i] - c[i + 1]
            if abs(denom) < 1e-12:
                return float(i + 0.5)
            return float(i + (c[i] - level) / denom)
    return None


__all__ = [
    "PedestalResult", "CentroidResult", "FwhmResult",
    "PEDESTAL_THRESHOLD_FRAC",
    "should_subtract_pedestal",
    "gost_pedestal_symmetric", "gost_pedestal_asymmetric",
    "gost_select_pedestal_method",
    "gost_centroid_graphical", "gost_centroid_weighted_mean",
    "gost_centroid_graphoanalytic",
    "gost_fwhm_linear_interp", "gost_fwhm_graphoanalytic",
    "gost_roi_from_fwhm",
]
