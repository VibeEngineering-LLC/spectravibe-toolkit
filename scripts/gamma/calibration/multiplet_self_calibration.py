"""
F-145 / v1.17.8 — Two-phase multiplet self-calibration.

Контракт «ЗАКРЕПЛЕНО НАВСЕГДА»:

  Когда мультиплет M фитируется со связанными интенсивностями (F-117/F-118)
  на NaI-спектре с дрейфом ADC-калибровки, контракт «центроиды locked на
  паспортные» переносит дрейф напрямую в χ². В качестве диагностики
  `coupled_intensity_fit(..., free_centroids=True)` запускает Phase A
  side-fit: per-component dE_k ∈ [±0.5·FWHM(E_k)] свободны, амплитуды
  пересчитываются на каждом шаге.

  Этот модуль реализует Phase B (convergence test) + Phase C (calibration
  refit) — собирает все Phase A результаты с успешной конвергенцией,
  отбирает достоверные сдвиги (max|dE|/FWHM < 0.5), формирует список
  anchor-точек {(channel_fitted, E_passport)}, докалибровывает спектр
  через :func:`gamma.calibration.energy_fit.polynomial_energy_fit` степени
  ≤ 2 (для NaI 63×63 деградация выше степени-2 несостоятельна), и
  откатывает калибровку, если новый residual хуже старого.

  Phase D (final locked-passport fit на пересчитанной шкале) выполняется
  outer pipeline'ом (staged_pipeline) повторным вызовом
  ``run_chain_forced_multiplets`` с ``free_centroids=False``.

См. также:
  - ROADMAP_v1_17_8_plus.md F-145 — обоснование и метрики
  - gamma.calibration.anchor_recalibration — F-87 Step 5β сестринская
    система (anchor-based recalibration по AnchorMatch'ам)
  - LSRM Algorithmic Foundations §8.4.4, Gilmore §6.4
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Iterable, List, Optional, Tuple, Dict, Any

import numpy as np

from gamma.calibration.energy_fit import (
    polynomial_energy_fit, MAX_POLYNOMIAL_DEGREE,
)


# F-145 пороги — закреплены навсегда
PHASE_B_MAX_DRIFT_FRACTION_OF_FWHM = 0.5    # |dE|/FWHM > 0.5 → отвергнуть anchor
# 2026-06-11 attempt to relax to 0.0 (operator methodology "fit-first → recal")
# деградировал cert_zcheck 3/3 → 1/3 на Th-232 Marinelli из-за:
#  (a) component swapping в Phase A (Ac-228 964.77/968.97 поменялись местами),
#  (b) deg=4 cal refit на 2 anchor'ах → over-fit.
# Возврат к 1.5 на момент починки degeneracy + degree-cap в Phase C.
# Рендер "fit-sum приклеена к спектру" реализован отдельно: json_report
# перерисовывает overlay на fitted centroids когда Phase A конвергнулась.
PHASE_B_MIN_CHI2_IMPROVEMENT = 1.5
PHASE_C_MIN_ANCHORS = 2                     # backward-compat: 2 anchor'а допустимо
                                            # на real NaI данных, где мультиплетов
                                            # часто только 2 (Pb-214 + Bi-214 ИЛИ
                                            # K-40 + Tl-208 в зависимости от типа).
PHASE_C_RECOMMENDED_MIN_ANCHORS = 3         # F-270 (v1.17.11, T-019) — рекомендация
                                            # для нового кода: ≥3 anchor'ов даёт
                                            # 1 избыточную точку для валидации.
                                            # Если real data позволяет — использовать
                                            # min_anchors=3 в production-вызовах.
# F-270 (v1.17.11, T-019) — range-gate: anchor'ы не должны кучковаться
# в узкой полосе. Refit deg≥1 с anchor'ами в 200 кэВ полосе будет
# экстраполировать на остальные ~2000 кэВ — это опаснее чем ничего.
# Порог 250 кэВ: близкие к границам мультиплеты (Pb-214 295+352 =
# 57 кэВ) комбинируются с дальними (Bi-214 609, 1120, 1764, 2614 и
# K-40 1461) — типичная картина даёт диапазон ≥ 250 кэВ.
PHASE_C_MIN_ANCHOR_RANGE_KEV = 250.0
PHASE_C_MAX_DEGREE_NAI = 2                  # для NaI 63x63 степень > 2 несостоятельна

# F-446 (v1.30.3 continuation) — adaptive Phase C polynomial degree для
# δ(N) refit. На малом числе anchor'ов polyfit deg ≥ 1 ЭКСТРАПОЛИРУЕТ
# линейно за пределы [ch_min, ch_max] anchor channels — что нарушает
# operator contract «расстояние между центроидами соответствует library
# spacing». На Th-232 ОИСН-16 cluster-Δ собрал ровно 2 anchor'а (M1 ~911,
# M2 ~1588), F-445 deg=1 polyfit ломал spacing на Pb-212 238 (ch≈80) и
# Tl-208 2614 (ch≈890) → Phase D откатывал из-за роста χ² мультиплетов.
#
# F-446 заменяет линейный trend на constant shift при n_anchors < 3.
# Constant shift = единый Δ_const ко всему диапазону каналов, spacing
# тривиально сохраняется везде. При n_anchors == 3 — линейный (1 anchor
# избыточен для валидации), при n_anchors >= 4 — parabola до cap NAI.
#
# Пороги hardcoded (а не отдельные константы): adaptive policy зашита
# в код, не настраивается через kwarg — это часть контракта F-446.
PHASE_C_ADAPTIVE_DEGREE = True              # F-446 marker (используется в diag/tests)
PHASE_C_MIN_ANCHORS_FOR_LINEAR = 3          # n<3 → deg=0
PHASE_C_MIN_ANCHORS_FOR_PARABOLA = 4        # n<4 → deg ≤ 1

# F-145 Phase D — смягчение locked-passport.
# После refit'a E(N) на 2-3 anchor'ах остаточные нелинейности шкалы
# приводят к ~0.1-0.3·FWHM нелокальным сдвигам в областях между anchor'ами.
# Жёсткий lock на паспортных приведёт к избыточному χ². Поэтому в Phase D
# центроиды свободны в малом окне ±PHASE_D_CENTROID_TOLERANCE_FRAC·FWHM(E_k).
# Это допускает «мягкое» прижатие к паспорту с эффектом регуляризации
# через ограниченную свободу. Параметр настраиваемый через kwarg
# `phase_D_tolerance_frac` в staged_pipeline.
#
# Значения по дефолту (для NaI 63×63, NaI Гамма-1С):
#   0.0  — жёсткий lock на паспортные (старый F-117/F-118 контракт)
#   0.10 — узкое окно, ~0.1·FWHM (для хорошо откалиброванных спектров)
#   0.15 — компромисс default, ~0.15·FWHM
#   0.25 — широкое (для спектров с заметным остаточным drift'ом)
#   0.50 — это уже Phase A bounds (полностью свободные центроиды)
#
# F-167 (2026-05-30) — переименован в `PHASE_D_REGULARIZATION_K_FWHM`
# и перенесён в каноничный справочник `gamma.identification.id_window`.
# Имя `PHASE_D_CENTROID_TOLERANCE_FRAC` сохраняется здесь как
# **deprecated alias** для backward-compat (импорты в `staged_pipeline.py`
# и legacy-вызовах). Использовать новое имя в новом коде.
from gamma.identification.id_window import (
    PHASE_D_REGULARIZATION_K_FWHM as _PHASE_D_REGULARIZATION_K_FWHM,
)
PHASE_D_CENTROID_TOLERANCE_FRAC = _PHASE_D_REGULARIZATION_K_FWHM


@dataclass
class CentroidAnchor:
    """Anchor-точка для refit E(N): фитированный канал ↔ паспортная энергия."""
    nuclide: str
    E_passport_keV: float
    E_fitted_keV: float      # = E_passport + centroid_shift
    channel_fitted: float    # ⟵ channel где данные показывают этот пик
    fwhm_keV: float
    drift_fraction_of_fwhm: float
    source: str              # "multiplet_M1", "Cs-Kα", "binding_singleton", ...


@dataclass
class SelfCalibrationDiag:
    """Полный диагностический блок для F-145."""
    attempted: bool = False
    phase_A_run: bool = False
    phase_B_passed: bool = False
    phase_C_applied: bool = False
    n_multiplets_seen: int = 0
    n_multiplets_phase_A_converged: int = 0
    n_anchors_collected: int = 0
    n_anchors_after_filter: int = 0
    anchors_used: List[Dict[str, Any]] = field(default_factory=list)
    old_energy_cal: Optional[list] = None
    new_energy_cal: Optional[list] = None
    old_residual_max_keV: Optional[float] = None
    new_residual_max_keV: Optional[float] = None
    degree_used: Optional[int] = None
    reason: str = ""
    phase_A_chi2_per_mult: Dict[str, float] = field(default_factory=dict)
    # F-446 diagnostics — adaptive Phase C degree policy
    delta_degree_used: Optional[int] = None        # actual deg used for δ(N): 0/1/2
    delta_const_keV: Optional[float] = None        # mean Δ при deg=0 (None иначе)
    degree_choice_reason: str = ""                 # «n_anchors=2 < 3 → constant shift»
    accepted_cluster_ids: List[Any] = field(default_factory=list)  # cid list (F-446
    # render-override guard support — set by caller when Phase D accepts cluster cal)


def _collect_anchors_from_multiplets(
    fitted_multiplets: Iterable,
    fwhm_provider_keV: Callable[[float], float],
    energy_to_channel: Callable[[float], float],
    diag: SelfCalibrationDiag,
) -> List[CentroidAnchor]:
    """Phase B implementation: фильтр Phase A результатов мультиплетов.

    Принимает iterable объектов с полями:
      • id или cluster_id
      • chi2_per_dof, phase_A_chi2_per_dof, phase_A_converged
      • components (List[ComponentFit] — у каждого .E_keV, .nuclide, .I_pct)
      • centroid_shifts_keV (List[float], выровнен по components)

    Для каждого мультиплета формируется ОДИН anchor через
    I_pct-взвешенное усреднение сдвигов центроидов:

        ⟨dE⟩ = Σ_k I_k · dE_k / Σ_k I_k
        ⟨E⟩  = Σ_k I_k · E_k  / Σ_k I_k

    Это снимает шум от слабых линий (где Phase A фиттер слабо
    ограничен) и сохраняет общий drift калибровки, под который
    подстраивается весь мультиплет.

    Фильтр Phase B:
      1. phase_A_converged == True
      2. phase_A_chi2_per_dof * PHASE_B_MIN_CHI2_IMPROVEMENT ≤ chi2_per_dof
         (Phase A улучшил χ² ≥ 1.5× относительно locked-passport)
      3. |⟨dE⟩| / FWHM(⟨E⟩) ≤ PHASE_B_MAX_DRIFT_FRACTION_OF_FWHM
         (взвешенный сдвиг в пределах ±0.5·FWHM)
    """
    anchors: List[CentroidAnchor] = []
    for m in fitted_multiplets:
        diag.n_multiplets_seen += 1
        m_id = str(getattr(m, "id", getattr(m, "cluster_id", "M?")))
        comps = list(getattr(m, "components", []) or [])
        shifts = list(getattr(m, "centroid_shifts_keV", []) or [])
        pA_chi2 = getattr(m, "phase_A_chi2_per_dof", None)
        pA_conv = bool(getattr(m, "phase_A_converged", False))
        chi2_locked = float(getattr(m, "chi2_per_dof", 0.0))
        if pA_chi2 is not None:
            diag.phase_A_chi2_per_mult[m_id] = float(pA_chi2)
        if not pA_conv or pA_chi2 is None:
            continue
        diag.n_multiplets_phase_A_converged += 1
        # Phase B condition 2: χ² improvement
        if not (pA_chi2 * PHASE_B_MIN_CHI2_IMPROVEMENT <= chi2_locked):
            continue
        if len(shifts) != len(comps) or not comps:
            continue
        # I_pct-weighted average shift and energy
        # Берём ТОЛЬКО компоненты, у которых индивидуальный |dE|/FWHM ≤ 1.0
        # (отсечка диких выбросов от ложных свободных параметров для очень
        # слабых линий — типа Bi-214 665 при 1.51% I_γ). Из оставшихся
        # формируем I_pct-взвешенное среднее.
        weights = []
        weighted_E = 0.0
        weighted_dE = 0.0
        sum_w = 0.0
        for k, comp in enumerate(comps):
            # MultipletComponent → line_E_keV; ComponentFit → E_keV.
            E_lib_raw = (
                getattr(comp, "line_E_keV", None)
                if getattr(comp, "line_E_keV", None) is not None
                else getattr(comp, "E_keV", 0.0)
            )
            E_lib = float(E_lib_raw or 0.0)
            if E_lib <= 0:
                continue
            dE = float(shifts[k])
            fwhm_k = max(1e-6, float(fwhm_provider_keV(E_lib)))
            # Outlier cut: per-component |dE|/FWHM > 1.0 → исключить
            if abs(dE) > 1.0 * fwhm_k:
                continue
            # Веса = I_pct (минимум 0.5 чтобы слабые не выпали полностью);
            # для одиночных групп independent — используем I_pct = 100 как fallback.
            # MultipletComponent → library_I_pct; ComponentFit → I_pct.
            I_pct = float(
                getattr(comp, "library_I_pct", None)
                if getattr(comp, "library_I_pct", None) is not None
                else getattr(comp, "I_pct", 100.0)
            ) or 100.0
            w_k = max(0.5, I_pct)
            weighted_E += w_k * E_lib
            weighted_dE += w_k * dE
            sum_w += w_k
            weights.append((E_lib, dE, w_k, fwhm_k))
        if sum_w <= 0 or not weights:
            continue
        E_avg = weighted_E / sum_w
        dE_avg = weighted_dE / sum_w
        fwhm_avg = max(1e-6, float(fwhm_provider_keV(E_avg)))
        frac_avg = abs(dE_avg) / fwhm_avg
        # Phase B condition 3: weighted drift cap
        if frac_avg > PHASE_B_MAX_DRIFT_FRACTION_OF_FWHM:
            continue
        E_fit_avg = E_avg + dE_avg
        try:
            ch_fit = float(energy_to_channel(E_fit_avg))
        except Exception:
            continue
        # Имя нуклида = первый non-empty из компонент
        nuc = ""
        for c in comps:
            n = str(getattr(c, "nuclide", ""))
            if n:
                nuc = n
                break
        anchors.append(CentroidAnchor(
            nuclide=nuc,
            E_passport_keV=E_avg,
            E_fitted_keV=E_fit_avg,
            channel_fitted=ch_fit,
            fwhm_keV=fwhm_avg,
            drift_fraction_of_fwhm=frac_avg,
            source=f"multiplet_{m_id}_I_pct_weighted",
        ))
        diag.n_anchors_collected += 1
    return anchors


def _f445_try_cluster_anchors(
    spec, fitted_multiplets, fwhm_provider_keV, diag,
    counts_arr, continuum_arrays, E_arr, use_cluster_global,
):
    """F-445: try cluster-Δ anchors. Returns [] on any failure."""
    if not (use_cluster_global and counts_arr is not None):
        return []
    try:
        from gamma.calibration.cluster_shift_anchors import (
            collect_cluster_global_anchors,
        )
        out = collect_cluster_global_anchors(
            spec, fitted_multiplets, counts_arr,
            continuum_arrays, E_arr, fwhm_provider_keV,
        )
        if out:
            diag.n_anchors_collected += len(out)
        return list(out)
    except Exception as exc:
        diag.reason = "F-445 cluster collector exception: " + repr(exc)
        return []


def _f446_compute_constant_delta(anchors, delta_targets) -> float:
    """F-446: constant Δ = uniform mean of Δ across anchors.

    Cluster-Δ + per-component смешиваются как равноценные. Per-component
    уже несут I_pct усреднение на уровне коллектора (см.
    _collect_anchors_from_multiplets), повторное взвешивание не требуется.
    Использует CentroidAnchor.weight если есть (future-proof), иначе
    uniform mean (weight=1.0).
    """
    weights = []
    for a in anchors:
        w = float(getattr(a, "weight", 1.0)) or 1.0
        weights.append(max(0.5, w))
    sum_w = sum(weights)
    if sum_w > 0:
        return sum(d * w for d, w in zip(delta_targets, weights)) / sum_w
    return sum(delta_targets) / max(1, len(delta_targets))


def _f446_choose_phase_c_degree(n_anchors: int, max_degree: int) -> Tuple[int, str]:
    """F-446: adaptive Phase C degree policy for δ(N) refit.

    n<3  → deg=0 (constant shift, preserves spacing everywhere).
    n==3 → deg=1 (linear, 1 redundant point for validation).
    n>=4 → deg=min(max_degree, n-1) (parabola up to NaI cap).
    """
    if n_anchors < PHASE_C_MIN_ANCHORS_FOR_LINEAR:
        return 0, (
            f"n_anchors={n_anchors} < {PHASE_C_MIN_ANCHORS_FOR_LINEAR} "
            f"→ constant shift (deg=0)"
        )
    if n_anchors < PHASE_C_MIN_ANCHORS_FOR_PARABOLA:
        return 1, f"n_anchors={n_anchors} == 3 → linear (deg=1)"
    deg = min(int(max_degree), n_anchors - 1)
    return deg, (
        f"n_anchors={n_anchors} ≥ {PHASE_C_MIN_ANCHORS_FOR_PARABOLA} "
        f"→ polyfit deg={deg} (cap max_degree={max_degree})"
    )


def recalibrate_from_multiplet_centroids(
    spec,
    fitted_multiplets: Iterable,
    *,
    fwhm_provider_keV: Callable[[float], float],
    extra_anchors: Optional[Iterable[Tuple[float, float, str]]] = None,
    max_degree: int = PHASE_C_MAX_DEGREE_NAI,
    min_anchors: int = PHASE_C_MIN_ANCHORS,
    counts_arr=None,
    continuum_arrays: Optional[dict] = None,
    E_arr=None,
    use_cluster_global: bool = True,
) -> Tuple[Optional[Tuple[float, ...]], SelfCalibrationDiag]:
    """F-145 Phase B + C — собирает anchor'ы и refit'ит E(N).

    Parameters
    ----------
    spec : Spectrum
        Исходный спектр со старой ``energy_cal``. Используются методы
        ``energy_to_channel(E)`` и поле ``energy_cal``.
    fitted_multiplets : iterable
        Результаты Phase A — объекты типа CoupledFitResult или
        DeconvolutionResult с полями ``components``, ``centroid_shifts_keV``,
        ``chi2_per_dof``, ``phase_A_chi2_per_dof``, ``phase_A_converged``.
    fwhm_provider_keV : Callable[[float], float]
        F(E) в кэВ.
    extra_anchors : iterable of (E_passport_keV, channel_fitted, source_label), optional
        Дополнительные anchor'ы из других источников (Cs-Kα 32 кэВ из
        F-142, binding-good singletons и т.п.).
    max_degree : int
        Cap для polyfit степени. Для NaI 63×63 — 2 (LSRM рекомендация).
    min_anchors : int
        Минимум anchor'ов для безопасного refit. < 3 → skip.

    Returns
    -------
    (new_cal, diag) : tuple
        ``new_cal`` — refitted коэффициенты (a0, a1, ...) или None при
        отказе. ``diag`` — полная диагностика F-145.
    """
    diag = SelfCalibrationDiag()
    diag.attempted = True
    diag.phase_A_run = True
    stored_cal = tuple(getattr(spec, "energy_cal", ()) or ())
    diag.old_energy_cal = list(stored_cal) if stored_cal else None

    def _e2c(E):
        if hasattr(spec, "energy_to_channel"):
            return spec.energy_to_channel(E)
        # Fallback: linear inverse
        if stored_cal and len(stored_cal) >= 2 and stored_cal[1] != 0:
            return (E - stored_cal[0]) / stored_cal[1]
        raise ValueError("spec has no energy_to_channel and no linear cal")

    # F-445 / v1.30.3: collect cluster-level Δ_cluster anchors AND
    # legacy per-component anchors, then merge unique by (nuclide,
    # E_passport). On conflict cluster-Δ wins (matches operator
    # contract: cluster-spacing preserved). Cluster-Δ ALONE was tried
    # in early F-445 iteration but regressed Ra-226 M1 χ²/ν 61.69 → 66.56
    # because Phase D ran with fewer anchors. Merge restores Phase D
    # statistical power on Ra-226 while still unlocking Th-232 cases
    # where per-comp filter rejected all anchors.
    cluster_anchors = _f445_try_cluster_anchors(
        spec, fitted_multiplets, fwhm_provider_keV, diag,
        counts_arr, continuum_arrays, E_arr, use_cluster_global,
    )
    legacy_anchors = _collect_anchors_from_multiplets(
        fitted_multiplets, fwhm_provider_keV, _e2c, diag,
    )
    seen_keys = set()
    anchors = []
    for a in list(cluster_anchors) + list(legacy_anchors):
        key = (str(a.nuclide), round(float(a.E_passport_keV), 2))
        if key in seen_keys:
            continue
        seen_keys.add(key)
        anchors.append(a)
    # Add extra anchors (Cs-Kα, binding singletons)
    if extra_anchors:
        for E_lib, ch_fit, src in extra_anchors:
            fwhm = max(1e-6, float(fwhm_provider_keV(E_lib)))
            anchors.append(CentroidAnchor(
                nuclide="",
                E_passport_keV=float(E_lib),
                E_fitted_keV=float(E_lib),  # для extra мы знаем passport
                channel_fitted=float(ch_fit),
                fwhm_keV=fwhm,
                drift_fraction_of_fwhm=0.0,
                source=str(src),
            ))
            diag.n_anchors_collected += 1

    diag.n_anchors_after_filter = len(anchors)
    diag.anchors_used = [
        {
            "nuclide": a.nuclide,
            "E_passport_keV": round(a.E_passport_keV, 3),
            "E_fitted_keV": round(a.E_fitted_keV, 3),
            "channel_fitted": round(a.channel_fitted, 3),
            "drift_keV": round(a.E_fitted_keV - a.E_passport_keV, 3),
            "drift_frac_of_FWHM": round(a.drift_fraction_of_fwhm, 3),
            "source": a.source,
        }
        for a in anchors
    ]

    if not anchors:
        diag.reason = "ни один мультиплет с конвергнувшей фазой А не прошёл фильтр"
        return None, diag

    if len(anchors) < min_anchors:
        diag.reason = (
            f"недостаточно anchor'ов для safe refit "
            f"({len(anchors)} < {min_anchors})"
        )
        diag.phase_B_passed = False
        return None, diag

    # F-270 (v1.17.11, T-019) — диапазон покрытия anchor'ов.
    # Даже 3 anchor'а в узкой полосе (например, все в 600–700 кэВ) не
    # обеспечивают safe refit deg-2 — это приведёт к экстраполяции.
    anchor_energies = [a.E_passport_keV for a in anchors]
    e_min = min(anchor_energies)
    e_max = max(anchor_energies)
    if (e_max - e_min) < PHASE_C_MIN_ANCHOR_RANGE_KEV:
        diag.reason = (
            f"anchor'ы покрывают только {e_max - e_min:.0f} кэВ "
            f"(E_min={e_min:.0f}, E_max={e_max:.0f}; "
            f"требуется ≥ {PHASE_C_MIN_ANCHOR_RANGE_KEV:.0f} кэВ) — "
            f"refit отклонён, риск экстраполяции"
        )
        diag.phase_B_passed = False
        return None, diag

    diag.phase_B_passed = True

    # Phase C: подгонка СДВИГА к stored E(N), а не замена полинома.
    # Сохраняем deg(stored_cal) старшие степени, добавляем low-degree
    # коррекцию δ(N) (a0 + a1·N для 2 anchor'ов; +a2·N² для ≥3).
    #
    #   E_corrected(N) = E_stored(N) + δ(N)
    #   δ(N) fitted via polyfit на (channel, E_passport - E_stored(channel))
    channels = [a.channel_fitted for a in anchors]
    energies_passport = [a.E_passport_keV for a in anchors]

    if stored_cal:
        def _e_old(N):
            E = 0.0
            for c in reversed(stored_cal):
                E = E * N + c
            return E
        old_resids = [abs(_e_old(ch) - E) for ch, E in zip(channels, energies_passport)]
        diag.old_residual_max_keV = float(max(old_resids))
        # Residuals to fit: E_passport - E_stored(channel)
        delta_targets = [
            float(E) - float(_e_old(ch))
            for ch, E in zip(channels, energies_passport)
        ]
        # F-446: adaptive Phase C degree policy.
        # n<3 → deg=0 constant; n==3 → deg=1; n>=4 → deg=min(max_degree, n-1).
        n_anchors_for_deg = len(anchors)
        delta_deg, choice_reason = _f446_choose_phase_c_degree(
            n_anchors_for_deg, int(max_degree)
        )
        diag.degree_choice_reason = choice_reason
        diag.delta_degree_used = int(delta_deg)
        if delta_deg == 0:
            # F-446: constant shift = mean Δ across anchor channels.
            # Cluster-Δ anchor'ы trivially равноценны (cluster-global centroid),
            # per-component anchor уже несёт I_pct-weighted dE внутри. Здесь
            # берём uniform mean чтобы не накладывать веса дважды.
            delta_const = _f446_compute_constant_delta(anchors, delta_targets)
            diag.delta_const_keV = float(delta_const)
            delta_lo_to_hi = [float(delta_const)]
        else:
            # Подгонка δ(N) — numpy.polyfit, high-to-low → переворачиваем
            delta_hi_to_lo = np.polyfit(
                np.asarray(channels, dtype=np.float64),
                np.asarray(delta_targets, dtype=np.float64),
                delta_deg,
            )
            delta_lo_to_hi = list(delta_hi_to_lo[::-1])
        # E_new(N) = E_stored(N) + δ(N) поэлементно
        n_max = max(len(stored_cal), len(delta_lo_to_hi))
        stored_padded = list(stored_cal) + [0.0] * (n_max - len(stored_cal))
        delta_padded = delta_lo_to_hi + [0.0] * (n_max - len(delta_lo_to_hi))
        new_cal_list = [
            float(stored_padded[i] + delta_padded[i]) for i in range(n_max)
        ]
        new_cal = tuple(new_cal_list)
        diag.new_energy_cal = list(new_cal)
        diag.degree_used = int(max(len(stored_cal), delta_deg + 1) - 1)
        # New residual
        def _e_new(N):
            E = 0.0
            for c in reversed(new_cal):
                E = E * N + c
            return E
        new_max = float(max(abs(_e_new(ch) - float(E))
                            for ch, E in zip(channels, energies_passport)))
        diag.new_residual_max_keV = new_max
    else:
        # Нет stored_cal — fallback на полный polyfit
        diag.old_residual_max_keV = float("inf")
        target = 0.3 * min(
            max(1e-6, fwhm_provider_keV(E)) for E in energies_passport
        )
        fit = polynomial_energy_fit(
            channels=channels,
            energies=energies_passport,
            max_degree=max_degree,
            target_residual_keV=target,
            min_degree=1,
        )
        new_cal = tuple(fit.coefficients)
        diag.new_energy_cal = list(new_cal)
        diag.degree_used = int(fit.degree)
        new_resids_pred = fit.predict(channels)
        new_max = float(max(abs(float(p) - float(E))
                            for p, E in zip(new_resids_pred, energies_passport)))
        diag.new_residual_max_keV = new_max

    # Решение: применять ли новую калибровку
    if diag.old_residual_max_keV is None:
        diag.phase_C_applied = True
        diag.reason = (
            f"калибровка применена (старая отсутствует, "
            f"new_residual_max={new_max:.2f} кэВ)"
        )
        return new_cal, diag

    if new_max < diag.old_residual_max_keV:
        diag.phase_C_applied = True
        diag.reason = (
            f"подгонка улучшила невязки: "
            f"{diag.old_residual_max_keV:.2f} → {new_max:.2f} кэВ "
            f"(степень {diag.degree_used}, {len(anchors)} опорных точек)"
        )
        return new_cal, diag

    diag.phase_C_applied = False
    diag.reason = (
        f"подгонка НЕ улучшила невязки "
        f"({diag.old_residual_max_keV:.2f} → {new_max:.2f} кэВ), "
        f"откат к сохранённой калибровке"
    )
    return None, diag


__all__ = [
    "CentroidAnchor",
    "SelfCalibrationDiag",
    "recalibrate_from_multiplet_centroids",
    "PHASE_B_MAX_DRIFT_FRACTION_OF_FWHM",
    "PHASE_B_MIN_CHI2_IMPROVEMENT",
    "PHASE_C_MIN_ANCHORS",
    "PHASE_C_RECOMMENDED_MIN_ANCHORS",
    "PHASE_C_MIN_ANCHOR_RANGE_KEV",
    "PHASE_C_MAX_DEGREE_NAI",
    "PHASE_C_ADAPTIVE_DEGREE",
    "PHASE_C_MIN_ANCHORS_FOR_LINEAR",
    "PHASE_C_MIN_ANCHORS_FOR_PARABOLA",
    "_f446_choose_phase_c_degree",
    "_f446_compute_constant_delta",
]
