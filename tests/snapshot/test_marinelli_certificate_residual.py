# -*- coding: utf-8 -*-
"""Wave 2 / 2026-06-04 — Real-data certificate residual validation
для Marinelli 1L reference kits (Phase 1 exit criterion #1).

Проверяет, что измеренная удельная активность (Bq/kg) для эталонных
Marinelli-источников из `detectors/Gamma-1S/reference_spectra/
reference_kits/Marinelli_1L/` сходится с паспортной активностью
(сертификат) в допустимом диапазоне.

# Паспортные активности (источник: .spe COMMENT-поле + .src certificate)

| Kit            | Cert SA  | Ref. date  | Mass    | Source                          |
|----------------|----------|------------|---------|---------------------------------|
| Cs-137         | 1890 Bq/kg ±5%  | 30.05.1997 | 0.570 kg | .spe line 29 COMMENT          |
| K-40           | 2540 Bq/kg ±10% | n/a        | 0.665 kg | .spe line 30 COMMENT          |
| Ra-226         | 1850 Bq/kg ±10% | n/a        | 0.622 kg | .spe line 29 COMMENT          |
| Th-232         | 1940 Bq/kg ±6%  | 17.09.2007 | 1.600 kg | .spe line 42 COMMENT + .src line 32 |

# Ground-truth nuclide per kit

NaI с разрешением ~7 % @662 кэВ плохо разделяет Ac-228 multiplets
(338+340+911+964+969 keV), поэтому daughter-line residuals варьируются.
Для Phase 1 exit criterion #1 (residual <5%) пригодны строки, дающие
clean singlet с высокой intensity:

- Ra-226 chain: Pb-214 295/352 keV (in secular equilibrium with Ra-226 →
  same Bq/kg) — currently -1.2% on real data (PASS for <5% target).
- Th-232 chain: Pb-212 238 keV (clean singlet, BR=43.6%, in secular
  equilibrium с Th-232) — -11.5% (NaI realistic, не достигает 5%).
- K-40: прямая 1460.8 keV строка — -12.8%.

Daughter-line «known weaker proxies» (Ac-228, Bi-214):
- Ac-228 911 keV — -40% (multiplet smearing; carry-over для Wave 3
  multiplet improvement).
- Bi-214 609 keV — -31% (overlap с другими Ra-226 chain lines).

# Tolerance policy (документировано как deliberate)

- TIER A (Phase 1 exit-criterion #1 «<5%»): Pb-214 / Ra-226 — ENFORCED.
- TIER B («realistic NaI ±15%»): K-40 direct, Pb-212 / Th-232 —
  ENFORCED (regression guard, не должно деградировать).
- TIER C («known underestimate»): Ac-228, Bi-214 — captured as data,
  asserted within ±50% (catches catastrophic regressions only).

# Anti-hallucination provenance (cited offsets)

- Cs-137 cert 1890 Bq/kg ref 30.05.97 → `detectors/Gamma-1S/
  reference_spectra/reference_kits/Marinelli_1L/Cs-137/
  sample_M_cs_легкий_2001-2005.spe` line 29 (COMMENT=).
- K-40 cert 2540 Bq/kg → same path /K-40/sample_M_k_легкий_2001-2005.spe
  line 30 (COMMENT=).
- Ra-226 cert 1850 Bq/kg → same /Ra-226/sample_M_ra_легкий_2001-2007.spe
  line 29 (COMMENT=).
- Th-232 cert 1940 Bq/kg ref 17.09.2007 → /Th-232/Th232_420-7-17_
  Маринелли_0cm.spe line 42 (COMMENT=Th-232 A=1940 Бк/кг dA=6%
  17-09-2007) + `detectors/Gamma-1S/certificates/
  Эталон_Маринелли__Аспект_2017_.src` lines 22-32
  (section `[..., Th232_420-7-17, Act]` → `Th-232=1940,6`).
- Sample masses из SAMPLEMASS-поля .spe (line 19/20).

# Carry-over (Phase 1 exit-criterion gap)

Тест документирует, что <5% residual goal достигнут ТОЛЬКО для одной
строки (Pb-214 / Ra-226). Остальные ground-truth strings показывают
-11..-13% (физически объяснимо — NaI 7% FWHM, Compton continuum,
self-absorption в matrix `Грунт-16`). Phase 2 RC может потребовать
HPGe data или re-fit efficiency curve для достижения <5% на всех
4 nuclides.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

KIT = REPO / "detectors/Gamma-1S/reference_spectra/reference_kits/Marinelli_1L"


# ─────────────────────────────────────────────────────────────────────
# Certificate ground truth (cited above)
# ─────────────────────────────────────────────────────────────────────

CERT = {
    "Cs-137": {
        "sample_spe": "Cs-137/sample_M_cs_легкий_2001-2005.spe",
        "bg_spe": "Cs-137/background_bg_2016_marinelli_water_marinelli.spe",
        "mass_kg": 0.570,
        "cert_SA_Bq_per_kg": 1890.0,
        "cert_uncertainty_pct": 5.0,
        # Decay 1997-05-30 → 1999-08-04 (MEASBEGIN) = 2.18 y, T½=30.17 y
        # factor = 0.5^(2.18/30.17) = 0.9506
        "decay_factor_at_meas": 0.9506,
        "ground_truth_line": "Cs-137",  # 661.66 keV — direct
        "tier": "C",  # currently not detected on this fixture (NaI limit)
    },
    "K-40": {
        "sample_spe": "K-40/sample_M_k_легкий_2001-2005.spe",
        "bg_spe": "K-40/background_bg_2016_marinelli_water_marinelli.spe",
        "mass_kg": 0.665,
        "cert_SA_Bq_per_kg": 2540.0,
        "cert_uncertainty_pct": 10.0,
        # K-40 T½=1.25e9 y → no decay correction over ~25 y
        "decay_factor_at_meas": 1.0,
        "ground_truth_line": "K-40",  # 1460.83 keV — direct
        "tier": "B",  # ~13% empirical
    },
    "Ra-226": {
        "sample_spe": "Ra-226/sample_M_ra_легкий_2001-2007.spe",
        "bg_spe": "Ra-226/background_bg_2016_marinelli_water_marinelli.spe",
        "mass_kg": 0.622,
        "cert_SA_Bq_per_kg": 1850.0,
        "cert_uncertainty_pct": 10.0,
        # Ra-226 T½=1600 y → no significant decay
        "decay_factor_at_meas": 1.0,
        "ground_truth_line": "Pb-214",  # 295.22/351.93 keV secular eq. proxy
        "tier": "A",  # <5% empirical
    },
    "Th-232": {
        "sample_spe": "Th-232/Th232_420-7-17_Маринелли_0cm.spe",
        "bg_spe": "Th-232/Фон закр кр вода_13.spe",
        "mass_kg": 1.600,
        "cert_SA_Bq_per_kg": 1940.0,
        "cert_uncertainty_pct": 6.0,
        # Th-232 T½=1.4e10 y → no decay
        "decay_factor_at_meas": 1.0,
        "ground_truth_line": "Pb-212",  # 238.63 keV (BR 43.6%) secular eq.
        "tier": "B",  # ~12% empirical
    },
}

# Tolerance tiers (deliberately documented above)
TIER_A_RESIDUAL_MAX = 0.05   # Phase 1 exit criterion #1 — <5%
TIER_B_RESIDUAL_MAX = 0.15   # «realistic NaI ±15%»
TIER_C_RESIDUAL_MAX = 0.50   # regression catch — catastrophic only


# ─────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────


def _expected_SA(entry: dict) -> float:
    """Cert SA скорректировано на decay до MEASBEGIN."""
    return entry["cert_SA_Bq_per_kg"] * entry["decay_factor_at_meas"]


def _run_analysis(entry: dict, out_dir: Path) -> dict:
    """Запускает analyze_and_report, возвращает report.json как dict."""
    from gamma.reporting import analyze_and_report

    sample = KIT / entry["sample_spe"]
    bg = KIT / entry["bg_spe"]
    out_dir.mkdir(parents=True, exist_ok=True)
    analyze_and_report(
        str(sample),
        background_path=str(bg),
        output_dir=str(out_dir),
        sample_mass_kg=entry["mass_kg"],
        write_json=True,
        write_markdown=False,
        write_html=False,
        write_plots=False,
        write_technical_pdf=False,
    )
    rep_path = next(out_dir.glob("*_report.json"), None)
    assert rep_path is not None, f"report.json не создан в {out_dir}"
    return json.loads(rep_path.read_text(encoding="utf-8"))


def _specific_activity_for(rep: dict, nuclide: str) -> float | None:
    """Извлекает specific_activity_Bq_per_kg для конкретного nuclide."""
    for n in rep.get("identified_nuclides", []):
        if n.get("nuclide") == nuclide:
            sa = n.get("specific_activity_Bq_per_kg")
            if sa is not None and sa > 0:
                return float(sa)
    return None


def _residual(measured: float, expected: float) -> float:
    """Relative residual = (measured − expected) / expected."""
    return (measured - expected) / expected


def _all_specific_activities(rep: dict) -> dict[str, float]:
    out = {}
    for n in rep.get("identified_nuclides", []):
        nm = n.get("nuclide")
        sa = n.get("specific_activity_Bq_per_kg")
        if nm and sa is not None and sa > 0:
            out[nm] = float(sa)
    return out


# ─────────────────────────────────────────────────────────────────────
# Fixtures availability gate
# ─────────────────────────────────────────────────────────────────────


@pytest.fixture(params=list(CERT.keys()), scope="module")
def kit_entry(request):
    """Per-kit fixture, кидает skip если файлов нет."""
    name = request.param
    entry = CERT[name]
    sample = KIT / entry["sample_spe"]
    bg = KIT / entry["bg_spe"]
    if not sample.is_file() or not bg.is_file():
        pytest.skip(f"{name} fixture missing: {sample} or {bg}")
    return name, entry


# ─────────────────────────────────────────────────────────────────────
# Tests
# ─────────────────────────────────────────────────────────────────────


def test_marinelli_certificates_loaded():
    """Sanity: все 4 паспортных значения присутствуют в CERT с положительной активностью."""
    assert set(CERT.keys()) == {"Cs-137", "K-40", "Ra-226", "Th-232"}
    for nm, e in CERT.items():
        assert e["cert_SA_Bq_per_kg"] > 0, f"{nm} cert SA нулевая"
        assert 0.5 <= e["mass_kg"] <= 2.0, f"{nm} mass вне Marinelli диапазона"
        assert 0.5 <= e["decay_factor_at_meas"] <= 1.0, (
            f"{nm} decay factor вне (0.5..1.0]"
        )


def test_marinelli_th232_pb212_residual_within_tier_B(tmp_path):
    """Th-232 Marinelli: Pb-212 (238 keV, secular eq.) residual ≤ ±15%.

    Phase 1 exit criterion #1 цель <5% НЕ достигается на NaI (текущий
    factual residual ~−11.5%). Tier-B contract: catch регресс > ±15%.
    """
    entry = CERT["Th-232"]
    sample = KIT / entry["sample_spe"]
    bg = KIT / entry["bg_spe"]
    if not sample.is_file() or not bg.is_file():
        pytest.skip(f"Th-232 fixture missing: {sample}")

    rep = _run_analysis(entry, tmp_path / "out")
    measured = _specific_activity_for(rep, "Pb-212")
    assert measured is not None, (
        "Pb-212 не обнаружен на Th-232 Marinelli — регрессия "
        f"identification. Обнаружено: {list(_all_specific_activities(rep))}"
    )
    expected = _expected_SA(entry)
    residual = _residual(measured, expected)
    assert abs(residual) <= TIER_B_RESIDUAL_MAX, (
        f"Th-232 (proxy Pb-212): measured={measured:.1f} Bq/kg vs "
        f"cert={expected:.1f} Bq/kg → residual={residual:+.1%} "
        f"(tier-B limit ±{TIER_B_RESIDUAL_MAX:.0%})"
    )


def test_marinelli_ra226_pb214_residual_within_tier_A(tmp_path):
    """Ra-226 Marinelli: Pb-214 (295/352 keV, secular eq.) ≤ ±5%.

    ЭТО ЕДИНСТВЕННАЯ строка, удовлетворяющая Phase 1 exit-criterion #1
    (<5% residual) на текущей prod-сборке. Регрессия здесь = серьёзная
    деградация identification/efficiency.
    """
    entry = CERT["Ra-226"]
    sample = KIT / entry["sample_spe"]
    bg = KIT / entry["bg_spe"]
    if not sample.is_file() or not bg.is_file():
        pytest.skip(f"Ra-226 fixture missing: {sample}")

    rep = _run_analysis(entry, tmp_path / "out")
    measured = _specific_activity_for(rep, "Pb-214")
    assert measured is not None, (
        "Pb-214 (Ra-226 secular-equilibrium proxy) не найден — регрессия. "
        f"Обнаружено: {list(_all_specific_activities(rep))}"
    )
    expected = _expected_SA(entry)
    residual = _residual(measured, expected)
    assert abs(residual) <= TIER_A_RESIDUAL_MAX, (
        f"Ra-226 (proxy Pb-214): measured={measured:.1f} Bq/kg vs "
        f"cert={expected:.1f} Bq/kg → residual={residual:+.1%} "
        f"(Phase 1 exit-criterion #1 ±{TIER_A_RESIDUAL_MAX:.0%})"
    )


def test_marinelli_k40_direct_residual_within_tier_B(tmp_path):
    """K-40 Marinelli: K-40 1460.83 keV ≤ ±15%.

    Direct line — нет цепочки распадов. Текущий empirical residual ~−13%.
    """
    entry = CERT["K-40"]
    sample = KIT / entry["sample_spe"]
    bg = KIT / entry["bg_spe"]
    if not sample.is_file() or not bg.is_file():
        pytest.skip(f"K-40 fixture missing: {sample}")

    rep = _run_analysis(entry, tmp_path / "out")
    measured = _specific_activity_for(rep, "K-40")
    assert measured is not None, (
        "K-40 не обнаружен на K-40 Marinelli — регрессия identification."
    )
    expected = _expected_SA(entry)
    residual = _residual(measured, expected)
    assert abs(residual) <= TIER_B_RESIDUAL_MAX, (
        f"K-40 direct: measured={measured:.1f} Bq/kg vs "
        f"cert={expected:.1f} Bq/kg → residual={residual:+.1%} "
        f"(tier-B ±{TIER_B_RESIDUAL_MAX:.0%})"
    )


def test_marinelli_th232_ac228_underestimate_tier_C(tmp_path):
    """Th-232 Marinelli: Ac-228 (911 keV proxy) — known underestimate.

    Ac-228 на NaI находится в multiplet с 904, 964, 969 keV — частично
    смешивается; в текущей сборке residual ≈ −40%. Tier-C: assert ≤ 50%
    (catches catastrophic regression только). Carry-over для Wave 3
    multiplet/Ac-228 improvement.
    """
    entry = CERT["Th-232"]
    sample = KIT / entry["sample_spe"]
    bg = KIT / entry["bg_spe"]
    if not sample.is_file() or not bg.is_file():
        pytest.skip(f"Th-232 fixture missing: {sample}")

    rep = _run_analysis(entry, tmp_path / "out")
    measured = _specific_activity_for(rep, "Ac-228")
    if measured is None:
        pytest.skip(
            "Ac-228 не идентифицирован в этой сборке — tier-C метрика "
            "недоступна (но не fail)."
        )
    expected = _expected_SA(entry)
    residual = _residual(measured, expected)
    assert abs(residual) <= TIER_C_RESIDUAL_MAX, (
        f"Ac-228 catastrophic regression: measured={measured:.1f} Bq/kg "
        f"vs cert={expected:.1f} Bq/kg → residual={residual:+.1%} "
        f"(tier-C catastrophic limit ±{TIER_C_RESIDUAL_MAX:.0%})"
    )


def test_marinelli_ra226_bi214_underestimate_tier_C(tmp_path):
    """Ra-226 Marinelli: Bi-214 (609 keV) — known overlap underestimate.

    Bi-214 609 keV в NaI overlaps с другими U/Th chain lines; tier-C.
    """
    entry = CERT["Ra-226"]
    sample = KIT / entry["sample_spe"]
    bg = KIT / entry["bg_spe"]
    if not sample.is_file() or not bg.is_file():
        pytest.skip(f"Ra-226 fixture missing: {sample}")

    rep = _run_analysis(entry, tmp_path / "out")
    measured = _specific_activity_for(rep, "Bi-214")
    if measured is None:
        pytest.skip("Bi-214 не идентифицирован — tier-C метрика недоступна.")
    expected = _expected_SA(entry)
    residual = _residual(measured, expected)
    assert abs(residual) <= TIER_C_RESIDUAL_MAX, (
        f"Bi-214 catastrophic regression: measured={measured:.1f} Bq/kg "
        f"vs cert={expected:.1f} Bq/kg → residual={residual:+.1%} "
        f"(tier-C ±{TIER_C_RESIDUAL_MAX:.0%})"
    )


def test_marinelli_phase1_exit_criterion_at_least_one_under_5pct(tmp_path):
    """Phase 1 exit-criterion #1: «certificate residual < 5%».

    Контракт: ХОТЯ БЫ ОДИН nuclide ground-truth line в Marinelli kits
    достигает |residual| ≤ 5%. Сегодня этому удовлетворяет Pb-214/Ra-226.
    Если в будущем (efficiency re-fit / better matrix correction)
    больше nuclides сходятся <5%, тест автоматически продолжит
    проходить.

    NOTE: НЕ закрывает Phase 1 exit-criterion полностью — для GA
    необходимо <5% на КАЖДОМ ground-truth nuclide. Эта точка покрытия —
    регресс-guard.
    """
    hits = []
    misses = []
    for name, entry in CERT.items():
        sample = KIT / entry["sample_spe"]
        bg = KIT / entry["bg_spe"]
        if not sample.is_file() or not bg.is_file():
            continue
        rep = _run_analysis(entry, tmp_path / f"{name}_out")
        measured = _specific_activity_for(rep, entry["ground_truth_line"])
        expected = _expected_SA(entry)
        if measured is None:
            misses.append(f"{name}/{entry['ground_truth_line']}=not_detected")
            continue
        residual = _residual(measured, expected)
        record = (
            f"{name} ({entry['ground_truth_line']}): "
            f"meas={measured:.1f} cert={expected:.1f} "
            f"res={residual:+.1%}"
        )
        if abs(residual) <= TIER_A_RESIDUAL_MAX:
            hits.append(record)
        else:
            misses.append(record)

    assert hits, (
        "Phase 1 exit-criterion #1: НЕТ ни одного nuclide с residual "
        f"≤ {TIER_A_RESIDUAL_MAX:.0%} на Marinelli kits.\n"
        f"Hits: {hits}\nMisses: {misses}"
    )
