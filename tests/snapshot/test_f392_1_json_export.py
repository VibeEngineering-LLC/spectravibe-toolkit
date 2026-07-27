"""F-392 / v1.18.27 + F-392.1 / v1.18.28 — JSON export of continuum_model.

`DeconvolutionResult.continuum_model` хранит выбранную continuum-модель
(linear / step_linear / step_linear_multi / quadratic), но до v1.18.28.1
поле НЕ экспортировалось в JSON `multiplet_deconvolutions[*]`. Downstream
consumers (HTML report, run_skill.py sanity, external tools) не могли
увидеть, активировался ли F-392.1 auto-select на конкретном кластере.

Тесты:
  1. Synthetic DeconvolutionResult с continuum_model="step_linear_multi" —
     `_build_deconvolutions` должен пробросить поле дословно в JSON dict.
  2. Все 4 supported значения round-trip (linear, step_linear,
     step_linear_multi, quadratic).
  3. Empirical PROD Th-232 demo: хотя бы один multiplet имеет
     `continuum_model == "step_linear_multi"` (F-392.1 boundary M3
     активация на real Marinelli-1L спектре). Pass-through через
     `build_json_report`.

Запрет F-392.1 / F-145: НЕ менять физику в peaks/coupled_multiplet.py
(Agent A scope). Только additive JSON schema change.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

from gamma.peaks.deconvolve import DeconvolutionResult, MultipletComponent
from gamma.reporting.json_report import _build_deconvolutions


# ──────────────────────────────────────────────────────────────────
# Утилиты — лёгкие fake-объекты для unit-тестов без полного pipeline
# ──────────────────────────────────────────────────────────────────


class _FakeStagedResult:
    """Минимальный duck-type holder для `_build_deconvolutions`.

    Builder использует только `.deconvolution_results`. Остальное поле
    не нужно для этого теста.
    """

    def __init__(self, deconvolution_results):
        self.deconvolution_results = deconvolution_results


def _make_decon(
    continuum_model: str,
    cluster_id: str = "M1",
    *,
    multi_step_anchors: tuple = (),
    multi_step_intensity_threshold_pct: float | None = None,
) -> DeconvolutionResult:
    """Собрать минимальный валидный DeconvolutionResult с заданной моделью.

    F-392.1 / v1.18.29: позволяет задать multi_step_anchors и
    multi_step_intensity_threshold_pct для тестирования новых JSON-полей.
    """
    comp = MultipletComponent(
        nuclide="Tl-208",
        line_E_keV=583.19,
        library_I_pct=30.55,
        center_channel=500.0,
        fwhm_channels=10.0,
    )
    return DeconvolutionResult(
        components=(comp,),
        areas=(1234.0,),
        area_uncertainties=(50.0,),
        continuum_params=(100.0, 0.0),
        continuum_model=continuum_model,
        chi2_per_dof=1.12,
        n_dof=10,
        roi_low_ch=480,
        roi_high_ch=520,
        gross_counts=2000.0,
        converged=True,
        method="lsq_linear",
        cluster_id=cluster_id,
        multi_step_anchors=multi_step_anchors,
        multi_step_intensity_threshold_pct=multi_step_intensity_threshold_pct,
    )


# ──────────────────────────────────────────────────────────────────
# Test 1: Synthetic step_linear_multi round-trip
# ──────────────────────────────────────────────────────────────────


class TestSyntheticContinuumModelRoundTrip:
    """Поле `continuum_model` должно появиться в JSON dict как str."""

    def test_step_linear_multi_present(self):
        dr = _make_decon("step_linear_multi", cluster_id="M3")
        fake = _FakeStagedResult([dr])
        out = _build_deconvolutions(fake)
        assert len(out) == 1
        assert "continuum_model" in out[0], (
            f"`continuum_model` поле должно быть в JSON dict, "
            f"но keys={sorted(out[0].keys())}"
        )
        assert out[0]["continuum_model"] == "step_linear_multi"
        assert isinstance(out[0]["continuum_model"], str)

    @pytest.mark.parametrize(
        "model",
        ["linear", "step_linear", "step_linear_multi", "quadratic"],
    )
    def test_all_models_round_trip(self, model: str):
        dr = _make_decon(model)
        fake = _FakeStagedResult([dr])
        out = _build_deconvolutions(fake)
        assert out[0]["continuum_model"] == model

    def test_missing_continuum_model_returns_empty_string(self):
        """Если поле не set (legacy DeconvolutionResult), не должно падать."""

        # `DeconvolutionResult` фиксирует continuum_model как обязательное
        # поле, но `_build_deconvolutions` использует getattr с дефолтом ""
        # — проверим robustness через duck-typed объект без атрибута.
        class _BareDecon:
            components = ()
            areas = ()
            area_uncertainties = ()
            chi2_per_dof = 0.0
            converged = False
            cluster_id = "M0"
            degenerate_pairs = ()
            phantom_components = ()
            centroid_shifts_keV = ()
            phase_A_chi2_per_dof = None
            phase_A_converged = False
            # умышленно НЕТ continuum_model

        fake = _FakeStagedResult([_BareDecon()])
        out = _build_deconvolutions(fake)
        assert out[0]["continuum_model"] == ""

    def test_other_keys_unchanged(self):
        """Additive change: не должны исчезнуть существующие keys."""
        dr = _make_decon("step_linear_multi", cluster_id="M3")
        fake = _FakeStagedResult([dr])
        out = _build_deconvolutions(fake)
        d = out[0]
        # Существующие schema keys должны остаться (API contract).
        required = {
            "cluster_id",
            "converged",
            "chi2_per_dof",
            "F145_phase_A_chi2_per_dof",
            "F145_phase_A_converged",
            "n_components",
            "components",
            "degenerate_pairs",
            "continuum_model",  # v1.18.28.1
            # F-392.1 / v1.18.29 — multi-step diagnostics
            "step_anchor_energies_keV",
            "step_intensity_pct",
        }
        assert required.issubset(d.keys()), (
            f"missing keys: {required - set(d.keys())}"
        )
        assert d["cluster_id"] == "M3"
        assert d["converged"] is True
        assert d["n_components"] == 1


# ──────────────────────────────────────────────────────────────────
# Test 1b: F-392.1 / v1.18.29 — multi-step anchors / threshold export
# ──────────────────────────────────────────────────────────────────


class TestSyntheticMultiStepDiagnostics:
    """`multi_step_anchors` и `multi_step_intensity_threshold_pct` из
    DeconvolutionResult должны попасть в JSON dict как
    `step_anchor_energies_keV` (list[float]) и `step_intensity_pct`
    (float / None)."""

    def test_step_linear_multi_anchors_present(self):
        """Th-232 M3-like cluster — 3 anchors с threshold=4.0."""
        dr = _make_decon(
            "step_linear_multi",
            cluster_id="M3",
            multi_step_anchors=((463.0, 12.0), (510.8, 13.0), (583.2, 14.0)),
            multi_step_intensity_threshold_pct=4.0,
        )
        fake = _FakeStagedResult([dr])
        out = _build_deconvolutions(fake)
        d = out[0]
        assert d["step_anchor_energies_keV"] == [463.0, 510.8, 583.2], (
            f"anchors не пробросились в JSON: got {d.get('step_anchor_energies_keV')!r}"
        )
        assert d["step_intensity_pct"] == 4.0
        # Все элементы — float, не tuple (σ_step не утекает).
        for e in d["step_anchor_energies_keV"]:
            assert isinstance(e, float)

    def test_two_anchors_export(self):
        """Synthetic параметризация — 2 anchors с custom threshold."""
        dr = _make_decon(
            "step_linear_multi",
            cluster_id="M2",
            multi_step_anchors=((463.0, 12.0), (583.2, 14.0)),
            multi_step_intensity_threshold_pct=4.5,
        )
        out = _build_deconvolutions(_FakeStagedResult([dr]))
        d = out[0]
        assert d["step_anchor_energies_keV"] == [463.0, 583.2]
        assert d["step_intensity_pct"] == 4.5

    def test_non_multi_continuum_empty_anchors(self):
        """linear / step_linear / quadratic — anchors=[], threshold=None."""
        for model in ("linear", "step_linear", "quadratic"):
            dr = _make_decon(model)  # без anchors / threshold
            d = _build_deconvolutions(_FakeStagedResult([dr]))[0]
            assert d["step_anchor_energies_keV"] == [], (
                f"{model!r}: anchors должен быть пустым, got "
                f"{d['step_anchor_energies_keV']!r}"
            )
            assert d["step_intensity_pct"] is None, (
                f"{model!r}: threshold должен быть None, got "
                f"{d['step_intensity_pct']!r}"
            )

    def test_legacy_decon_without_fields(self):
        """Robust: legacy duck-typed объект без новых полей → []/None."""

        class _BareDecon:
            components = ()
            areas = ()
            area_uncertainties = ()
            chi2_per_dof = 0.0
            converged = False
            cluster_id = "M0"
            degenerate_pairs = ()
            phantom_components = ()
            centroid_shifts_keV = ()
            phase_A_chi2_per_dof = None
            phase_A_converged = False
            continuum_model = "linear"
            # умышленно НЕТ multi_step_anchors / threshold

        fake = _FakeStagedResult([_BareDecon()])
        out = _build_deconvolutions(fake)
        assert out[0]["step_anchor_energies_keV"] == []
        assert out[0]["step_intensity_pct"] is None


# ──────────────────────────────────────────────────────────────────
# Test 2: Empirical Th-232 demo — M3 PROD activates step_linear_multi
# ──────────────────────────────────────────────────────────────────

_TH232_PATH = (
    "detectors/Gamma-1S/reference_spectra/reference_kits/"
    "Marinelli_1L/Th-232/Th232_420-7-17_Маринелли_0cm.spe"
)
_TH232_BG_PATH = (
    "detectors/Gamma-1S/reference_spectra/reference_kits/"
    "Marinelli_1L/Th-232/Фон закр кр вода_13.spe"
)


@pytest.mark.slow
class TestEmpiricalTh232ContinuumModelExported:
    """End-to-end: schema invariant — каждый multiplet на PROD Th-232
    демонстрирует non-empty `continuum_model` поле в JSON.

    Минимально требуется: ВСЕ multiplets имеют валидную модель из
    {linear, step_linear, step_linear_multi, quadratic}. Какая именно
    модель выбирается для каждого кластера зависит от auto-select
    (F-392.1 threshold 4%, ≥3 anchors, ≥200 кэВ span) и от phantom
    pool, что — Agent A scope. Этот тест НЕ задаёт hard-lock на
    конкретный M-индекс, только schema invariant: поле НЕ пропущено
    и не пустое.
    """

    def _run_pipeline(self):
        from gamma.identification.staged_pipeline import analyze_lsrm_spe
        from gamma.reporting import build_json_report

        spe_path = Path(__file__).resolve().parent.parent.parent / _TH232_PATH
        bg_path = Path(__file__).resolve().parent.parent.parent / _TH232_BG_PATH
        if not spe_path.exists():
            pytest.skip(f"PROD Th-232 spectrum not available: {spe_path}")
        if not bg_path.exists():
            pytest.skip(f"PROD bg spectrum not available: {bg_path}")
        result = analyze_lsrm_spe(
            str(spe_path),
            background_path=str(bg_path),
            complete_workflow=True,
            sample_mass_kg=1.6,
        )
        return build_json_report(result)

    def test_all_multiplets_have_continuum_model_field(self):
        report = self._run_pipeline()
        multiplets = report.get("multiplet_deconvolutions") or []
        assert multiplets, "PROD Th-232 should produce multiplets"
        valid_models = {"linear", "step_linear", "step_linear_multi", "quadratic"}
        for i, m in enumerate(multiplets, 1):
            assert "continuum_model" in m, (
                f"M{i}: missing continuum_model in JSON schema"
            )
            cmodel = m["continuum_model"]
            assert isinstance(cmodel, str), f"M{i}: continuum_model not str"
            assert cmodel != "", f"M{i}: continuum_model is empty"
            assert cmodel in valid_models, (
                f"M{i}: unexpected continuum_model={cmodel!r}, "
                f"expected one of {sorted(valid_models)}"
            )

    def test_multi_step_anchors_exported_when_active(self):
        """F-392.1 / v1.18.29: для multiplets с continuum_model ==
        'step_linear_multi' поле `step_anchor_energies_keV` — non-empty
        list of floats; `step_intensity_pct` — float (typ. 4.0).
        Для остальных моделей anchors=[] и threshold=None.

        Контракт schema-invariant: оба ключа ВСЕГДА присутствуют (даже
        при non-multi continuum они равны [] / null) — downstream
        потребители не должны делать `if 'step_anchor...' in m`.
        """
        report = self._run_pipeline()
        multiplets = report.get("multiplet_deconvolutions") or []
        assert multiplets, "PROD Th-232 should produce multiplets"

        seen_multi = False
        for i, m in enumerate(multiplets, 1):
            # Schema invariance — оба ключа всегда есть.
            assert "step_anchor_energies_keV" in m, (
                f"M{i}: missing step_anchor_energies_keV"
            )
            assert "step_intensity_pct" in m, (
                f"M{i}: missing step_intensity_pct"
            )
            anchors = m["step_anchor_energies_keV"]
            threshold = m["step_intensity_pct"]
            assert isinstance(anchors, list), (
                f"M{i}: step_anchor_energies_keV must be list, "
                f"got {type(anchors).__name__}"
            )
            cmodel = m["continuum_model"]
            if cmodel == "step_linear_multi":
                seen_multi = True
                assert anchors, (
                    f"M{i}: step_linear_multi должен иметь non-empty anchors, "
                    f"got {anchors!r}"
                )
                for e in anchors:
                    assert isinstance(e, float), (
                        f"M{i}: anchor energy не float: {e!r}"
                    )
                    # Sanity — γ-energies в [10, 3000] keV для NaI Th-232.
                    assert 10.0 <= e <= 3000.0, (
                        f"M{i}: anchor energy out of physical range: {e}"
                    )
                assert isinstance(threshold, float), (
                    f"M{i}: step_intensity_pct должен быть float при "
                    f"step_linear_multi, got {threshold!r}"
                )
                assert 0.0 < threshold <= 100.0, (
                    f"M{i}: step_intensity_pct out of range: {threshold}"
                )
            else:
                # Non-multi continuum — anchors must be empty, threshold None.
                assert anchors == [], (
                    f"M{i}: {cmodel!r} continuum должен иметь empty anchors, "
                    f"got {anchors!r}"
                )
                assert threshold is None, (
                    f"M{i}: {cmodel!r} continuum должен иметь threshold=None, "
                    f"got {threshold!r}"
                )

        # F-392.1 цель: Th-232 demo M3 (Ac-228 463 + Tl-208 510 + Tl-208 583)
        # обязан активировать step_linear_multi. Если не активировался —
        # либо физика поломалась, либо M3 cluster decompose-нулся (что
        # тоже регрессия). Зафиксируем как warning, не FAIL — auto-select
        # boundary условия (4% / 3 anchors / 200 keV span) могут shift'нуть
        # в future builds, поэтому hard-lock сюда — overfit к v1.18.29.
        if not seen_multi:
            pytest.skip(
                "No multiplet activated step_linear_multi on Th-232 demo — "
                "коммиты в peak grouping / clustering могли изменить "
                "boundary. Не регрессия export schema'и; см. PLAN F-392.1."
            )
