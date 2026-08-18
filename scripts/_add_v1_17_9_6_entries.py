from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
"""Append v1.17.9.6 deep-read entries to knowledge_index.json.

Idempotent: skips entries whose id already exists.
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
IDX = ROOT / "references" / "knowledge_index.json"

NEW_ENTRIES = [
    {
        "id": "LSRM-Algo-3.1.1-PeakSearchSigma",
        "book": "lsrm_algorithmic_foundations",
        "section": "§3.1.1",
        "page_from": 6, "page_to": 7,
        "topic_ru": "Convolution peak search detection threshold (default 3σ)",
        "topic_en": "Convolution peak search detection threshold",
        "keywords": ["convolution", "peak search", "threshold", "3 sigma", "Mariscotti"],
        "summary_ru": "Свёрточный быстрый поиск пиков. Рекомендуемый порог детекции = 3σ. Увеличение порога — пропуски пиков; уменьшение — статистические false-positive. Калибровочная ошибка FWHM >2× делает результаты ненадёжными.",
        "code_citations": ["peak_search.py default sigma_threshold"],
        "f157_priority": "PRIORITY-1-LSRM-OFFICIAL",
        "gost_ref_num": 7
    },
    {
        "id": "LSRM-Algo-5.2-ExclusionThresholds",
        "book": "lsrm_algorithmic_foundations",
        "section": "§5.2",
        "page_from": 14, "page_to": 16,
        "topic_ru": "Sequential nonlinear-parameter exclusion thresholds (dS/S 0.05/0.1/1.0)",
        "topic_en": "Sequential nonlinear parameter exclusion in fit",
        "keywords": ["nonlinear fit", "parameter exclusion", "Compton step", "FWHM", "convergence"],
        "summary_ru": "При separable nonlinear MNK поэтапно исключаются параметры на пиках низкой статистики: dS/S>0.05 → исключить Compton-step; >0.1 → исключить FWHM; >1.0 → исключить ВСЕ нелинейные. Предотвращает расхождение фита на слабых пиках.",
        "code_citations": ["coupled_multiplet.py может игнорировать; T-052"],
        "f157_priority": "PRIORITY-1-LSRM-OFFICIAL",
        "gost_ref_num": 7
    },
    {
        "id": "LSRM-Algo-6-IDwindow-sqrtE",
        "book": "lsrm_algorithmic_foundations",
        "section": "§6",
        "page_from": 18, "page_to": 19,
        "topic_ru": "ID-window scaling через √E-якорь Cs-137: Δ_E=(Δ_E^Cs/√E_Cs)·√E",
        "topic_en": "Identification window energy scaling rooted at Cs-137",
        "keywords": ["ID window", "identification", "√E scaling", "NaI", "HPGe", "Cs-137"],
        "summary_ru": "Окно идентификации зависит от энергии через √E-якорь Cs-137. Δ_E^Cs ≈ 10-20 keV для NaI/CsI, ≈1 keV для HPGe. Формула: Δ_E = (Δ_E^Cs/√E_Cs)·√E. Это ДРУГАЯ формула чем k·FWHM(E) — они расходятся на низких/высоких E. Кандидат на пересмотр F-167.",
        "code_citations": ["F-167 id_window.py; T-001"],
        "f157_priority": "PRIORITY-1-LSRM-OFFICIAL",
        "gost_ref_num": 7
    },
    {
        "id": "LSRM-Algo-6.3-ISO11929-quad",
        "book": "lsrm_algorithmic_foundations",
        "section": "§6.3",
        "page_from": 20, "page_to": 23,
        "topic_ru": "ISO 11929 decision threshold a* и detection limit a# (quadratic form)",
        "topic_en": "ISO 11929 detection limit full quadratic formula",
        "keywords": ["ISO 11929", "MDA", "decision threshold", "detection limit", "Type-B"],
        "summary_ru": "ISO 11929 формулы 6.3-1...6.3-9. Decision threshold a* = k_α·u(0). Detection limit a# из квадратного уравнения из-за a-dependent variance. Отдельные формы для paired-blank (sample n_b/t_b + чистый bg n_g/t_g) и area-domain S*, S#. Включает Type-B u_rel(g) термин для efficiency/intensity uncertainty.",
        "code_citations": ["mda.py; T-021, T-057"],
        "f157_priority": "PRIORITY-1-LSRM-OFFICIAL",
        "gost_ref_num": 7
    },
    {
        "id": "LSRM-Algo-8.4.1-PeakImage",
        "book": "lsrm_algorithmic_foundations",
        "section": "§8.4.1",
        "page_from": 33, "page_to": 35,
        "topic_ru": "Tabulated peak-image (smoothed-spline) calibration",
        "topic_en": "Tabular peak-image calibration",
        "keywords": ["peak image", "spline", "log-spline", "tabulated", "shape calibration"],
        "summary_ru": "Реальный отклик спектрометра хранится как сглаженно-сплайновая табулированная функция. Bg-полиномы в двух боковых окнах ~3 FWHM. Сплайн (или log-spline для лучшей dynamics) интерполирует. Процедура масштабирования к другой position/FWHM. Skill: peak-image layer ОТСУТСТВУЕТ; добавление улучшит multiplet deconvolution.",
        "code_citations": ["T-023"],
        "f157_priority": "PRIORITY-1-LSRM-OFFICIAL",
        "gost_ref_num": 7
    },
    {
        "id": "LSRM-Algo-8.4.2.1-NaIpureGauss",
        "book": "lsrm_algorithmic_foundations",
        "section": "§8.4.2.1",
        "page_from": 36, "page_to": 37,
        "topic_ru": "NaI peak shape = pure Gaussian (NO tail; tail только HPGe/alpha)",
        "topic_en": "NaI peak shape is pure Gaussian without tail",
        "keywords": ["NaI", "peak shape", "Gaussian", "no tail", "HPGe", "alpha"],
        "summary_ru": "Для гамма-пиков NaI рекомендуется чистая гауссиана (без tail). Для гамма-пиков HPGe и alpha-пиков — модифицированная гауссиана с low-energy экспоненциальным хвостом; tail-параметр T выбирается golden-section χ²-минимизацией. Skill: убедиться, что не добавляется tail безусловно (T-003).",
        "code_citations": ["coupled_multiplet.py peak shape selection; T-003"],
        "f157_priority": "PRIORITY-1-LSRM-OFFICIAL",
        "gost_ref_num": 7
    },
    {
        "id": "LSRM-Algo-8.5.2-NewSelfAbs",
        "book": "lsrm_algorithmic_foundations",
        "section": "§8.5.2",
        "page_from": 44, "page_to": 47,
        "topic_ru": "NEW self-absorption thin-plate model + d_eff table 8-2",
        "topic_en": "New thin-plate self-absorption with d_eff lookup table",
        "keywords": ["self-absorption", "thin plate", "d_eff", "Marinelli", "Petri", "Дента", "SpectraLine 1.4"],
        "summary_ru": "Используется с SpectraLine 1.4: g(E)=(1−exp(−μρd))/(μρd) вместо exp(−μρd_eff). Новая table d_eff (vs старая): Marinelli 0.5L=15±2mm (vs 8mm), 1.0L=26±2mm (vs 17mm), 3.0L=60±5mm (vs 30mm), Petri 0.075L=15±2mm, Дента 0.12L=36±2mm. Корректирует до десятков % bias для Am-241 59 keV в плотных пробах.",
        "code_citations": ["self_attenuation.py; T-002, T-024"],
        "f157_priority": "PRIORITY-1-LSRM-OFFICIAL",
        "gost_ref_num": 7
    },
    {
        "id": "LSRM-Algo-11-UpperLimitGate",
        "book": "lsrm_algorithmic_foundations",
        "section": "§11",
        "page_from": 49, "page_to": 49,
        "topic_ru": "Activity uncertainty >50% → upper-limit report (A+ΔA only)",
        "topic_en": "Upper limit reporting when relative uncertainty >50%",
        "keywords": ["upper limit", "50% threshold", "insignificant", "activity"],
        "summary_ru": "Если относительная неопределённость активности >50% — сообщать только верхнюю границу (A_upp = A + ΔA), нуклид маркируется 'insignificant'. Порог конфигурируется (ErrorUpperLevel=50% default в SpectraLine).",
        "code_citations": ["json_report.py; T-016"],
        "f157_priority": "PRIORITY-1-LSRM-OFFICIAL",
        "gost_ref_num": 7
    },
    {
        "id": "LSRM-Algo-14.3-ConfidenceIndex",
        "book": "lsrm_algorithmic_foundations",
        "section": "§14.3",
        "page_from": 51, "page_to": 52,
        "topic_ru": "Confidence Index CI=log10(1/Π(δE·δI)) + per-detector Table 14-1 thresholds",
        "topic_en": "Confidence Index for identification quality + thresholds table",
        "keywords": ["confidence index", "CI", "identification", "NaI thresholds", "Table 14-1"],
        "summary_ru": "Эвристический показатель уникальности сигнатуры нуклида. CI = log10(1/(δE_1·δE_2·…·δI_1·δI_2·…)). Higher CI = более надёжная идентификация. Таблица 14-1 для NaI: Cs-137=1.8, K-40=2.2, Na-22=3.8, Cs-134=4.4, Co-60=5.9, Ba-133=8.5, Th-232=16.6, Eu-152=18.3. Skill: модуль ОТСУТСТВУЕТ — T-022.",
        "code_citations": ["NEW: confidence_index.py; T-022"],
        "f157_priority": "PRIORITY-1-LSRM-OFFICIAL",
        "gost_ref_num": 7
    },
    {
        "id": "LSRM-Algo-15.1-DeadTimeRecipe",
        "book": "lsrm_algorithmic_foundations",
        "section": "§15.1",
        "page_from": 53, "page_to": 54,
        "topic_ru": "Dead-time empirical calibration: Co-60 ref + 2 high-load sources",
        "topic_en": "Dead-time A,B calibration recipe",
        "keywords": ["dead time", "τ_d", "A B coefficients", "Co-60", "high load", "calibration"],
        "summary_ru": "Effective dead-time τ_d = A·Σy_i + B·Σ(y_i·i). Коэффициенты A,B определяются эмпирически: источник 1 = Co-60 ≤500 cps для reference rate n₁; источники 2 и 3 дают r₂, r₃ при ≤5·10⁴ cps. Методика валидна до 5·10⁴ cps с uncertainty ≤5%. Skill: коэффициенты для Gamma-1S = None → F-95 не применяется.",
        "code_citations": ["dead_time.py; T-012"],
        "f157_priority": "PRIORITY-1-LSRM-OFFICIAL",
        "gost_ref_num": 7
    },
    {
        "id": "LSRM-Algo-16-EffectiveCenter",
        "book": "lsrm_algorithmic_foundations",
        "section": "§16",
        "page_from": 55, "page_to": 55,
        "topic_ru": "Effective detector center d_eff(E) two-distance method",
        "topic_en": "Energy-dependent effective center via two-distance source",
        "keywords": ["effective center", "d_eff", "two-distance", "1/R²", "geometry"],
        "summary_ru": "Замена реального детектора точкой на d_eff за крышкой; энергозависима. Определяется измерением одного источника на двух дистанциях R₁, R₂ → счётности n₁, n₂. d_eff(E) = (R2·√n2 − R1·√n1) / (√n1 − √n2). 1/R² breaks для d<25 cm в лабораторной геометрии. .efd файл = 2-col text (E_MeV, distance_cm).",
        "code_citations": ["NEW: efd_loader.py; T-081"],
        "f157_priority": "PRIORITY-1-LSRM-OFFICIAL",
        "gost_ref_num": 7
    },
    {
        "id": "LSRM-Act-2014-4.1-AULranges",
        "book": "lsrm_activity_counting_samples_2024",
        "section": "§4.1",
        "page_from": 4, "page_to": 5,
        "topic_ru": "Range AUL Cs-137/K-40/Th-232/Ra-226 для Marinelli 1L NaI 63×63 Gamma-1S",
        "topic_en": "AUL ranges per nuclide for Gamma-1S Marinelli 1L",
        "keywords": ["Gamma-1S", "Marinelli 1L", "AUL", "Cs-137", "K-40", "Th-232", "Ra-226", "Bq/kg"],
        "summary_ru": "Диапазоны АУ (1h count, 50% uncertainty, no interferents): Cs-137 = 2-10⁵ Bq/kg, K-40 = 40-10⁵, Th-232 = 4-10⁵, Ra-226 = 4-10⁵. Это методические целевые пороги для Gamma-1S 63×63 NaI с Marinelli 1L. Skill может валидировать что достигаемая чувствительность лежит в этом диапазоне.",
        "code_citations": ["regression tests benchmark; T-008"],
        "f157_priority": "PRIORITY-1-LSRM-OFFICIAL",
        "gost_ref_num": 17
    },
    {
        "id": "LSRM-Act-2014-5.1-SpectrometerCompliance",
        "book": "lsrm_activity_counting_samples_2024",
        "section": "§5.1",
        "page_from": 6, "page_to": 6,
        "topic_ru": "Spectrometer compliance: NaI ≤8% @ 661 keV, INL ≤1%, stab ≤1%/24h",
        "topic_en": "Spectrometer compliance requirements per methodology",
        "keywords": ["NaI compliance", "energy resolution", "INL", "temporal stability", "Gamma-1S"],
        "summary_ru": "Requirements: NaI ≥40×40 мм OR LaBr ≥25×25 мм; разрешение @ 661 keV ≤8% NaI / ≤3.5 keV HPGe; интегральная нелинейность ±1% NaI / ±0.05% HPGe; временная стабильность 24h ≤1% NaI / ≤0.1% HPGe. Skill должен refuse применять методику если spec не соблюден.",
        "code_citations": ["compliance gate; T-083"],
        "f157_priority": "PRIORITY-1-LSRM-OFFICIAL",
        "gost_ref_num": 17
    },
    {
        "id": "LSRM-Act-2014-App5-MIA",
        "book": "lsrm_activity_counting_samples_2024",
        "section": "App 5",
        "page_from": 21, "page_to": 22,
        "topic_ru": "МИА (Minimum Measurable Activity) = 4.42·√n_f/(η·√t), 4h cap",
        "topic_en": "Minimum measurable activity methodical form",
        "keywords": ["MIA", "МИА", "minimum activity", "background rate", "sensitivity", "4 hours"],
        "summary_ru": "Single-isolated-nuclide МИА: 4.42·√n_f/(η·√t), где n_f = bg count rate (cps), η = sensitivity (cps/Bq). Scaling: MIA(t) = MIA(t₀)·√(t₀/t). Time selection: T > t_0·MIA²(t_0)/A_min². Рекомендуемое время 10-60 мин; >4h impractical (Type-A < Type-B уже).",
        "code_citations": ["mda.py; T-014"],
        "f157_priority": "PRIORITY-1-LSRM-OFFICIAL",
        "gost_ref_num": 17
    },
    {
        "id": "LSRM-Act-2014-App4-BgDriftGate",
        "book": "lsrm_activity_counting_samples_2024",
        "section": "App 4",
        "page_from": 18, "page_to": 19,
        "topic_ru": "Background F-statistic drift gate: |F-F₀|/F₀ ≤ 0.1, F=Σy_n/t_f",
        "topic_en": "Background drift control gate via integral F-statistic",
        "keywords": ["background drift", "F statistic", "10% gate", "operational check"],
        "summary_ru": "Ежедневный integral F = Σy_n/t_f по полному energy-range. Comparison |F−F₀|/F₀ ≤ 0.1 (10% порог). 10 мин достаточно если Σy>1000. >0.1 → investigate / re-measure background. Skill: pre-flight check перед использованием prior bg.",
        "code_citations": ["bg_lines_apriori.py; T-054"],
        "f157_priority": "PRIORITY-1-LSRM-OFFICIAL",
        "gost_ref_num": 17
    },
    {
        "id": "LSRM-Dose-2000-SensMatrix",
        "book": "dose_rate_lsrm_2000",
        "section": "§12.1-12.3",
        "page_from": 8, "page_to": 10,
        "topic_ru": "Sensitivity matrix T_ik (12 intervals × 16 windows) + μ_en/ρ + f_k(10) tables",
        "topic_en": "Dose-rate sensitivity matrix and conversion tables",
        "keywords": ["dose rate", "sensitivity matrix", "T_ik", "μ_en/ρ", "f_k(10)", "12 intervals", "16 windows"],
        "summary_ru": "Восстановление истинного спектра G_k из аппаратурных окон S_i: S_i=Σ T_ik·G_k. Ẋ_k = D·G_k·(μ_en/ρ)·Ē_k. H*(10) = Σ Ẋ_k·f_k(10). Энергетические интервалы (m=12): 40,70,100,140,180,250,420,600,750,950,1400,2000,3000 keV. Окна (n=16): 40,70,100,140,180,250,300,420,600,700,800,920,1050,1400,1800,2200,3000 keV. Полная таблица μ_en/ρ для air ρ=1.205 kg/m³.",
        "code_citations": ["NEW: dose_rate.py; T-008, T-009"],
        "f157_priority": "PRIORITY-1-LSRM-OFFICIAL",
        "gost_ref_num": 11
    },
    {
        "id": "Vartanov-6-MultipletResolvability",
        "book": "vartanov_practical_scint_djvu",
        "section": "§6 Fig.20",
        "page_from": 35, "page_to": 37,
        "topic_ru": "Multiplet resolvability NaI: n2 ≥ 1.81·(FWHM/2) equal-intensity doublet",
        "topic_en": "Multiplet resolvability criterion for NaI",
        "keywords": ["multiplet", "resolvability", "NaI", "Vartanov", "Fig.20", "deconvolution criterion"],
        "summary_ru": "Для близких гауссианов: для равной интенсивности (ratio=1) два пика разрешаются когда n2 = 1.81·(FWHM/2). При intensity ratio = 1.13 второй пик исчезает (только inflection остаётся). Это hard upper bound на multiplet detectability в NaI. Phase-D multiplet finder НЕ должен пытаться разрешать ниже этих порогов.",
        "code_citations": ["coupled_multiplet.py multiplet trigger"],
        "f157_priority": "PRIORITY-2-RUSSIAN-CANON",
        "gost_ref_num": 20
    },
    {
        "id": "Vartanov-3-IodineKescape",
        "book": "vartanov_practical_scint_djvu",
        "section": "§3, §12",
        "page_from": 13, "page_to": 16,
        "topic_ru": "Iodine K-escape peak NaI ~28.6 keV (3% area >200 keV, 40% @ 30-40 keV)",
        "topic_en": "Iodine K-escape peak NaI quantitative",
        "keywords": ["iodine K-escape", "NaI", "28 keV", "Am-241", "low-energy satellite"],
        "summary_ru": "Для NaI photopeak <200 keV, I-Kα/Kβ при 28 keV может escape и создать satellite 28 keV ниже main peak. Формула: P_X/P_y = (0.875)·(0.84)·(τ_K)·[…]. Доля: ~3% при 200 keV, ~40% при 30-40 keV. При E_γ ≥ 200 keV escape можно lump в photopeak. Для Am-241 (59 keV) — must model как отдельный Gaussian.",
        "code_citations": ["NEW: residual_classifier добавить LBL_XRAY_ESCAPE; T-025"],
        "f157_priority": "PRIORITY-2-RUSSIAN-CANON",
        "gost_ref_num": 20
    },
    {
        "id": "Vartanov-12-13-TCS",
        "book": "vartanov_practical_scint_djvu",
        "section": "§12-13",
        "page_from": 77, "page_to": 84,
        "topic_ru": "TCS explicit formulas NaI: 4π well + cylindrical geometry",
        "topic_en": "True coincidence summing explicit formulas",
        "keywords": ["TCS", "cascade summing", "4pi well", "Co-60", "Eu-152", "close geometry"],
        "summary_ru": "Source-near-detector (S/D <5 cm), cascade γ-emitters: true coincidence summing уменьшает S1, S2 и создаёт sum-peak в E1+E2. S1_obs = N0·ε(E1)·P(E1)·Ω·[1 − ε(E2)·Ω·W(0°)]. Sum-peak: S12 = N0·ε(E1)·P(E1)·ε(E2)·P(E2)·Ω²·W(0°). 4π well-crystal: ε(E_i) extracted из sum-peak ratio without external calibration. Effect scales Ω² (~1/d⁴), negligible @ d>20 cm для типичного 70×70 NaI.",
        "code_citations": ["cascade_summing.py; T-007"],
        "f157_priority": "PRIORITY-2-RUSSIAN-CANON",
        "gost_ref_num": 20
    },
    {
        "id": "Vartanov-14-Zimmermann",
        "book": "vartanov_practical_scint_djvu",
        "section": "§14",
        "page_from": 85, "page_to": 88,
        "topic_ru": "Zimmermann linearization centroid: ln g(n)=α−(n−n0)²/(2σ²)",
        "topic_en": "Zimmermann linearization for centroid",
        "keywords": ["centroid", "Zimmermann", "linearization", "Gauss fit", "Vartanov"],
        "summary_ru": "Plot ln g(n) vs n в peak region — для Gaussian парабола. Tangent intersections с abscissa дают n0 и σ. Vartanov Eq.(59): ln g(n) = α − n²/(2σ²). Более точно чем naive max-channel. Statistical error: σ_z = ±√(g/g(n−1)²). Area: S = σ·√(2π)·g(n0) = 1.064·FWHM·H (±2.65% error).",
        "code_citations": ["centroid_uncertainty.py optional verification mode"],
        "f157_priority": "PRIORITY-2-RUSSIAN-CANON",
        "gost_ref_num": 20
    },
    {
        "id": "Gilmore-16.3.5-186keV-NORM",
        "book": "pgs_gilmore_2008",
        "section": "Ch.16.3.5",
        "page_from": 319, "page_to": 322,
        "topic_ru": "186 keV NORM apportionment 235U/226Ra на NaI: 0.5709/0.02662",
        "topic_en": "186 keV NORM math apportionment",
        "keywords": ["186 keV", "235U", "226Ra", "NORM", "secular equilibrium", "NaI", "apportionment"],
        "summary_ru": "При natural isotopic ratio + secular equilibrium: Corrected 226Ra = 0.5709·Apparent_186keV (±0.86%). Estimated 235U = 0.02662·Apparent_186keV (±2.16%). Включает correction для small 230Th (186.05 keV, 0.0088%). Без apportionment NORM 226Ra на NaI завышается на ~43%. Mathematical apportionment ТОЛЬКО при natural ratio + equilibrium.",
        "code_citations": ["NEW: nuclide_disambiguation.py; T-005"],
        "f157_priority": "PRIORITY-3-INTERNATIONAL",
        "gost_ref_num": 19
    },
    {
        "id": "Gilmore-6.5-FWHM-sqrtquad",
        "book": "pgs_gilmore_2008",
        "section": "Ch.6.5",
        "page_from": 138, "page_to": 141,
        "topic_ru": "Sqrt-quadratic FWHM: √(e²+p²·E+c²·E²) (empirical F=0.108)",
        "topic_en": "Square-root quadratic FWHM(E) calibration",
        "keywords": ["FWHM model", "sqrt-quadratic", "Fano factor", "HPGe", "calibration"],
        "summary_ru": "Gilmore pooled 22 HPGe calibrations. RMS-difference vs energy: sqrt-quadratic 0.0057 (best), simple quadratic 0.0065, linear 0.020, Debertin-Helmer 0.026, Genie-2000 (a+b√E) 0.055 (worst). Эмпирические параметры: e=0.956 keV, p=0.0422, c=5.29e-4. Implied Fano F=0.108. Для NaI: W(%) = a + b/√E (Knoll alternative).",
        "code_citations": ["fwhm_fit.py; F-168 confirmed; T-021"],
        "f157_priority": "PRIORITY-3-INTERNATIONAL",
        "gost_ref_num": 19
    },
    {
        "id": "Gilmore-10.6.7-LaBrIntrinsic",
        "book": "pgs_gilmore_2008",
        "section": "Ch.10.6.7",
        "page_from": 213, "page_to": 215,
        "topic_ru": "LaBr3(Ce) 138La intrinsic activity makes it UNUSABLE for low-bg",
        "topic_en": "LaBr3 138La intrinsic activity rules out low-background use",
        "keywords": ["LaBr3", "138La", "intrinsic activity", "low background", "scintillator selection"],
        "summary_ru": "LaBr3(Ce) имеет ~1 Bq/cm³ intrinsic 138La activity. Bremsstrahlung до 255 keV; 788.74 keV summed with beta; 1435.80 keV summed с Ba X-rays. UNUSABLE для low-background work. Несмотря на 2.5× лучшее разрешение vs NaI — НЕ рекомендуется как замена NaI для low-activity NORM.",
        "code_citations": ["doc-level rule; T-034"],
        "f157_priority": "PRIORITY-3-INTERNATIONAL",
        "gost_ref_num": 19
    },
    {
        "id": "Budyka-7.6-h2-12",
        "book": "budyka_textbook",
        "section": "§7.6 eq.7.22",
        "page_from": 190, "page_to": 191,
        "topic_ru": "Channel-discretization FWHM correction: σ²_actual = σ²_gauss − h²/12",
        "topic_en": "Channel-bin discretization variance correction",
        "keywords": ["FWHM", "discretization", "channel bin", "correction", "Gauss"],
        "summary_ru": "Поправка на discretization: σ²_actual = σ²_gauss − h²/12; FWHM_G² = FWHM_A² + 0.462·h². Table 7.1: при FWHM=3 канала погрешность 3%; ≥15 каналов → <0.1%. Для Gamma-1S @ 661 keV FWHM≈40 keV/0.6 keV/ch ≈ 67 каналов — погрешность <0.05%. Для низких E может стать значимой.",
        "code_citations": ["fwhm calibration QC"],
        "f157_priority": "PRIORITY-2-RUSSIAN-CANON",
        "gost_ref_num": 12
    },
    {
        "id": "Budyka-7.7-eq7.29-SmoothStepBg",
        "book": "budyka_textbook",
        "section": "§7.7 eq.7.29",
        "page_from": 193, "page_to": 196,
        "topic_ru": "Smoothed-step background under FEP (sigmoid)",
        "topic_en": "Sigmoid step background under photopeak",
        "keywords": ["background", "sigmoid step", "FEP", "Compton step", "forward Compton"],
        "summary_ru": "B_n = C_L − (C_L−C_U)·Σ_{i=L+1..n}C_i / Σ_{i=L+1..U}C_i. Эмпирически учитывает forward-Compton edge + incomplete charge collection. Более физично чем pure linear для broad NaI peaks на Compton step. Net peak areas ~1-3% larger чем при linear background.",
        "code_citations": ["deconvolve.py continuum model; T-040, T-074"],
        "f157_priority": "PRIORITY-2-RUSSIAN-CANON",
        "gost_ref_num": 12
    },
    {
        "id": "Budyka-8.2-LQ",
        "book": "budyka_textbook",
        "section": "§8.2 eq.8.16-17",
        "page_from": 217, "page_to": 220,
        "topic_ru": "L_Q (limit of quantitation): 50·(1+√(1+B·(n+n²/m)/25))",
        "topic_en": "Limit of quantitation for 10% target uncertainty",
        "keywords": ["L_Q", "quantitation limit", "Currie", "Будыка", "third threshold"],
        "summary_ru": "Third statistical threshold above L_D. L_Q = k_Q²/2·(1+√(1+4σ_0²/k_Q²)). С peak area: L_Q ≈ 50·(1+√(1+B·(n+n²/m)/25)). Для типичного B-limited: L_Q ≈ 3·L_D. Reportable vs detectable threshold. Skill пока stop at L_D; добавление L_Q даёт regulator-friendly третий тиер.",
        "code_citations": ["mda.py extension; T-036"],
        "f157_priority": "PRIORITY-2-RUSSIAN-CANON",
        "gost_ref_num": 12
    },
    {
        "id": "ORTEC-A.2.2-n30winds-NaI",
        "book": "ortec_gammavision_v9_a66",
        "section": "App A.2.2 (n30winds.ini)",
        "page_from": 380, "page_to": 385,
        "topic_ru": "NaI32 engine defaults DIFFERENT from HPGe (LibReduction OFF, PeakOverlap=2.0)",
        "topic_en": "ORTEC NaI32 engine defaults differ from HPGe",
        "keywords": ["NaI32", "n30winds.ini", "library reduction", "peak overlap", "ORTEC", "NaI defaults"],
        "summary_ru": "n30winds.ini для NaI32: Library Nuclide Reduction Flag = FALSE (vs HPGe TRUE); Library Peak Critical Level Test Flag = FALSE (vs HPGe TRUE); Peak overlap range = 2.0·FWHM (vs HPGe 3.5). Skill, если использует HPGe defaults для NaI/Gamma-1S, это методологическая ошибка.",
        "code_citations": ["staged_pipeline.py default config; T-004"],
        "f157_priority": "PRIORITY-3-INTERNATIONAL",
        "gost_ref_num": 22
    },
    {
        "id": "ORTEC-6.17.2.2-RMax",
        "book": "ortec_gammavision_v9_a66",
        "section": "§6.17.2.2 eq.157-163",
        "page_from": 322, "page_to": 325,
        "topic_ru": "ISO NORM MDA Ratio RMax cap (default 3.0); β-risk dynamic adjustment",
        "topic_en": "ORTEC RMax cap on ISO NORM MDA/CL ratio",
        "keywords": ["RMax", "ISO NORM", "MDA cap", "ORTEC extension", "β risk"],
        "summary_ru": "Когда uncertainty очень большая, MDA correction factor f → ∞. ORTEC clamps R = MDA/CL at RMax (default 3.0, range 2-1000 via n30winds.ini). При превышении RMax β-risk dynamically RAISED пока ratio не вернётся к RMax. Asymptotic k_{1-β}=0 ↔ β=50%. Flag 'I' рядом с isotope name при capping. ISO 11929 не описывает; это ORTEC extension.",
        "code_citations": ["mda.py if ISO NORM mode; T-066"],
        "f157_priority": "PRIORITY-3-INTERNATIONAL",
        "gost_ref_num": 22
    },
    {
        "id": "MDA-paper-OneBestLine",
        "book": "mda_basics_ru",
        "section": "Заключение",
        "page_from": 5, "page_to": 5,
        "topic_ru": "MDA reported для ONE best line only — НЕ averaged across nuclide lines",
        "topic_en": "MDA reporting rule: single best line",
        "keywords": ["MDA", "one best line", "averaging", "Currie", "Lochamy", "Zimmer"],
        "summary_ru": "Если активность рассчитывается как weighted mean по N γ-линиям, MDA сообщается для ОДНОЙ best линии — highest quantum yield × ε / background, cleanest isolation. НЕ average. Skill: проверить что не averaged. MDA result ALWAYS prefixed '<' (less than) sign.",
        "code_citations": ["mda.py reporting; T-080"],
        "f157_priority": "PRIORITY-2-RUSSIAN-CANON",
        "gost_ref_num": 5
    },
    {
        "id": "Kuvkin-Ident-Table14-1",
        "book": "lsrm_algorithmic_foundations",
        "section": "§14.3 Table 14-1 (репродуцировано в 4 лекциях Кувыкина)",
        "page_from": 52, "page_to": 52,
        "topic_ru": "Per-detector CI thresholds Table 14-1: HPGe / NaI / LaBr3 per nuclide",
        "topic_en": "CI thresholds per detector and nuclide",
        "keywords": ["CI", "confidence index", "Table 14-1", "HPGe", "NaI", "LaBr3", "thresholds", "Кувыкин"],
        "summary_ru": "Минимальные CI для confident identification per detector class. NaI (для Gamma-1S): Cs-137=1.8, K-40=2.2, Na-22=3.8, Cs-134=4.4, Co-60=5.9, Ba-133=8.5, Th-232=16.6, Eu-152=18.3. HPGe значения 3-100× выше из-за лучшего разрешения. LaBr3 промежуточно. Skill: использовать как minimum confidence thresholds для NaI ID.",
        "code_citations": ["NEW: confidence_index.py thresholds table; T-022, NEW-LECT-001"],
        "f157_priority": "PRIORITY-1-LSRM-OFFICIAL",
        "gost_ref_num": 7
    }
]


def main() -> int:
    data = json.loads(IDX.read_text(encoding="utf-8"))
    existing_ids = {e["id"] for e in data["entries"]}
    added = 0
    skipped = 0
    for new in NEW_ENTRIES:
        if new["id"] in existing_ids:
            skipped += 1
            continue
        data["entries"].append(new)
        added += 1
    data["version"] = "1.4-v1.17.9.6"
    data["generated_at"] = "2026-05-30"
    data["source"] = data.get("source", "") + " v1.17.9.6 += 30 new entries from full-corpus deep-read."
    IDX.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Added: {added}, Skipped (already present): {skipped}")
    print(f"Total entries now: {len(data['entries'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
