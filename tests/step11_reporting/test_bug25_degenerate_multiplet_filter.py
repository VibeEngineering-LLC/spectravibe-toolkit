"""BUG-25 / v1.18.31+ (Agent B) — отсечение вырожденных мультиплетов
(все компоненты S=0 ИЛИ катастрофический χ²/ν) на reporting layer.

Symptom (Th-232 demo M6 «Ac-228 674.8 + Ac-228 726.9 + Bi-212 727.3 keV»):
  χ²/ν = 958.48, все три компонента S=0 (Ac-228 674.8, Ac-228 726.9,
  Bi-212 727.3 → S=0). График показывал bleeding от Tl-208 583 + малый
  bump near 720 + sloped continuum, но реального fit-сигнала не было.

Root cause: deconvolution возвращает результат NNLS даже когда
оптимум — все коэффициенты 0 (только континуум). Reporting слой
рендерил эти phantom-мультиплеты без квалитативного фильтра.

Fix: новый `_is_meaningful_multiplet` (interactive_html.py) применяется
в `render_interactive_html` к sample- и bg-сериям ДО `_build_multiplet_blocks`
и сериализации в DATA_MULTIPLETS. Гейт «meaningful»:
  • хотя бы один компонент с S > 0 (если σ_S задан — требуется S/σ ≥ 3);
  • chi2_per_dof < MAX_ACCEPTABLE_CHI2_DOF (= 1000.0).

JSON multiplet_deconvolutions НЕ трогается — full audit-trail сохраняется,
фильтр только на представлении.
"""
from __future__ import annotations

import sys

sys.path.insert(0, "scripts")

from gamma.reporting.interactive_html import (  # noqa: E402
    _is_meaningful_multiplet,
    _filter_meaningful_multiplets,
    MAX_ACCEPTABLE_CHI2_DOF,
)


def _mp(*, id_, chi2, areas):
    """Helper: построить multiplet-dict в формате `_build_multiplets_data`."""
    return {
        "id": id_,
        "chi2_per_dof": chi2,
        "components": [
            {"nuclide": f"X-{i}", "E_keV": 100.0 + i, "I_pct": 1.0,
             "area": float(a)}
            for i, a in enumerate(areas)
        ],
    }


def test_bug25_all_S_zero_multiplet_filtered():
    """Все компоненты S=0 → multiplet отбрасывается."""
    mp = _mp(id_="M6", chi2=958.48, areas=[0.0, 0.0, 0.0])
    assert not _is_meaningful_multiplet(mp), (
        "BUG-25: multiplet с χ²/ν=958.48 и всеми S=0 должен быть "
        "квалифицирован как вырожденный (has_signal=False)."
    )


def test_bug25_one_significant_component_passes():
    """Один компонент с большой S (даже при χ²/ν~150) → multiplet остаётся."""
    mp = _mp(id_="M4", chi2=180.05, areas=[54751.0, 4569.0, 205922.0])
    assert _is_meaningful_multiplet(mp), (
        "BUG-25: multiplet с реальным сигналом (Tl-208 510 S=54751 и др.) "
        "не должен отбрасываться даже при χ²/ν~180 — это качество подгонки, "
        "а не отсутствие сигнала."
    )


def test_bug25_catastrophic_chi2_filtered_even_with_signal():
    """χ²/ν выше потолка → multiplet отбрасывается (catastrophic fit)."""
    mp = _mp(id_="M99", chi2=MAX_ACCEPTABLE_CHI2_DOF + 1.0, areas=[1000.0])
    assert not _is_meaningful_multiplet(mp), (
        "BUG-25: catastrophic chi2/dof > MAX_ACCEPTABLE_CHI2_DOF должен "
        "отбрасываться даже при наличии non-zero S — fit не доверяется."
    )


def test_bug25_existing_demo_multiplets_pass():
    """Существующие Th-232 demo мультиплеты M1-M4 должны проходить фильтр:
    M1 χ²=29.23 (4/4 nonzero), M2 χ²=1.05 (3/3), M3 χ²=931.19 (6/9 nonzero),
    M4 χ²=180.05 (3/10 nonzero) — все имеют реальный сигнал."""
    cases = [
        ("M1", 29.23, [116355.0, 22504.0, 71256.0, 0.0]),
        ("M2", 1.05, [20661.0, 0.0, 10266.0]),
        ("M3", 931.19, [570049.0, 0.0, 0.0, 0.0, 0.0, 43153.0]),
        ("M4", 180.05, [54751.0, 4569.0, 205922.0, 0.0, 0.0, 0.0]),
    ]
    for id_, chi, areas in cases:
        mp = _mp(id_=id_, chi2=chi, areas=areas)
        assert _is_meaningful_multiplet(mp), (
            f"BUG-25 regression: existing Th-232 demo {id_} (χ²/ν={chi}) "
            f"должен оставаться видимым — фильтр слишком агрессивный."
        )


def test_bug25_filter_drops_only_degenerate():
    """Списочный helper: должен оставить только meaningful multiplets."""
    mps = [
        _mp(id_="M1", chi2=2.0, areas=[100.0, 200.0]),   # keep
        _mp(id_="M2", chi2=958.0, areas=[0.0, 0.0]),     # drop (all S=0)
        _mp(id_="M3", chi2=5.0, areas=[50.0]),           # keep
        _mp(id_="M4", chi2=10000.0, areas=[100.0]),      # drop (catastrophic)
    ]
    kept = _filter_meaningful_multiplets(mps)
    kept_ids = [m["id"] for m in kept]
    assert kept_ids == ["M1", "M3"], (
        f"BUG-25: filter должен оставить только M1 и M3 (meaningful), "
        f"actual={kept_ids}."
    )


def test_bug25_snr_threshold_filters_noise_components():
    """Компонент с S>0 но S/σ_S < 3 (low SNR) НЕ считается значимым.
    Если ВСЕ компоненты ниже SNR-порога → multiplet отбрасывается."""
    mp = {
        "id": "Mnoise",
        "chi2_per_dof": 2.0,
        "components": [
            {"nuclide": "X", "E_keV": 100.0, "I_pct": 1.0,
             "area": 10.0, "area_sigma": 50.0},   # SNR=0.2 → noise
        ],
    }
    assert not _is_meaningful_multiplet(mp), (
        "BUG-25: компонент с S/σ_S=0.2 — статистический шум, multiplet "
        "не должен квалифицироваться как meaningful."
    )
    # Поднимем S до значимого уровня → должен пройти.
    mp["components"][0]["area"] = 200.0  # SNR = 200/50 = 4 ≥ 3
    assert _is_meaningful_multiplet(mp), (
        "BUG-25 regression: S/σ_S=4 ≥ MIN_COMPONENT_SNR=3 — компонент "
        "должен квалифицироваться как значимый."
    )


if __name__ == "__main__":
    test_bug25_all_S_zero_multiplet_filtered()
    test_bug25_one_significant_component_passes()
    test_bug25_catastrophic_chi2_filtered_even_with_signal()
    test_bug25_existing_demo_multiplets_pass()
    test_bug25_filter_drops_only_degenerate()
    test_bug25_snr_threshold_filters_noise_components()
    print("OK")
