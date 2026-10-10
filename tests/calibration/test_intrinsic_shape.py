"""Тесты калибровки по форме собственных структур на авторской паре LaBr3 (фон + образец с Lu-176).
Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-10-10. Два теста selftest медленные (минуты), помечены slow."""

import re
import json
import pathlib
import numpy as np
import pytest

from gamma.calibration.intrinsic_shape.core import (
    load_spectrum,
    rebin_conserve_var,
    poly_eval,
    poly_inverse,
)
from gamma.calibration.intrinsic_shape import fwhm as F
from gamma.calibration.intrinsic_shape import register as reg
from gamma.calibration.intrinsic_shape import ops
from gamma.calibration.intrinsic_shape import selftest
from gamma.calibration.intrinsic_shape.cli import main as cli_main

FIX = pathlib.Path(__file__).resolve().parents[1] / "fixtures" / "intrinsic_shape"
BG = FIX / "labr3_background.xml"
SAMPLE = FIX / "labr3_lu176_sample.xml"
H_REF = [0.06290719, 1.00669051, -4.00234e-06]
ECAL_REF = [-0.48496698, 1.35007262, 4.58971597e-04]
FW_POINTS = [(29.20, 2.47), (142.87, 6.14), (211.97, 7.29), (338.14, 8.62), (847.10, 11.63)]
PRESETS = pathlib.Path(reg.__file__).resolve().parent / "presets"


def test_fixtures_guid_zeroed_and_live_time():
    """Проверка нулевого GUID и времени жизни в фикстурах."""
    for path in (BG, SAMPLE):
        text = path.read_text("utf-8")
        assert "<Guid>00000000-0000-0000-0000-000000000000</Guid>" in text
        uuids = re.findall(
            r"[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}",
            text,
        )
        assert all(u == "00000000-0000-0000-0000-000000000000" for u in uuids)

    bg = load_spectrum(BG)
    sm = load_spectrum(SAMPLE)
    assert len(bg["counts"]) == 2048
    assert len(sm["counts"]) == 2048
    assert bg["live"] == pytest.approx(291077.867571)
    assert sm["live"] == pytest.approx(244689.14413721429)


def test_rebin_conserves_counts():
    """Проверка сохранения суммарного числа событий при ребининге."""
    rng = np.random.default_rng(1)
    counts = rng.poisson(50, 300).astype(float)
    edges = np.arange(301) - 0.5 + 0.37
    values, var = rebin_conserve_var(counts, edges, 320)
    assert values.sum() == pytest.approx(counts.sum(), rel=1e-12)
    assert np.all(var <= values + 1e-9)


def test_poly_inverse_roundtrip():
    """Проверка обратимости полиномиального преобразования энергии."""
    for E in (40.0, 500.0, 1460.0):
        x = poly_inverse(ECAL_REF, E)
        assert poly_eval(ECAL_REF, x) == pytest.approx(E, abs=1e-6)


def test_fwhm_vertex_margin_fixed():
    """Проверка ограничения на вершину параболы FWHM."""
    x = [p[0] for p in FW_POINTS]
    f = [p[1] for p in FW_POINTS]
    c_new = F.fit_fwhm(2048, x, f)
    ok, vals, vertex, margin = F.check(c_new, 2048)
    assert ok
    assert margin >= 1000

    c_old = F.fit_fwhm(2048, x, f, monotone_to=1.0)
    assert F.check(c_old, 2048)[3] < 1.0

    c_c2 = F.fit_fwhm(2048, x, f, monotone_to="c2>=0")
    assert c_c2[2] >= -1e-15


def test_fwhm_block_channel_integer_energy_consistent():
    """Проверка целочисленности каналов и согласованности энергий в блоке FWHM."""
    block = F.fwhm_block("\n", [1.0, 0.2, 0.0], [(29.2, 2.47), (847.6, 11.6)], [0.0, 1.5, 0.0])
    channels = [int(v) for v in re.findall(r"<Channel>(-?\d+)</Channel>", block)]
    assert channels == [29, 848]
    energies = [float(v) for v in re.findall(r"<Energy>([^<]+)</Energy>", block)]
    assert energies == pytest.approx([43.5, 1272.0])
    assert "<SqrtFwhmCalibration>" in block


def test_make_net_xml_products(tmp_path):
    """Проверка XML-продуктов make_net: целые значения, суммы, FWHM."""
    r = ops.make_net(
        BG,
        SAMPLE,
        H_REF,
        ECAL_REF,
        FW_POINTS,
        tmp_path / "zeroed.xml",
        tmp_path / "signed.xml",
        date="2026-10-10",
    )

    for fname in ("zeroed.xml", "signed.xml"):
        text = (tmp_path / fname).read_text("utf-8")
        pts = re.findall(r"<DataPoint>([^<]*)</DataPoint>", text)
        assert len(pts) == 2048
        assert all(re.fullmatch(r"-?\d+", p) for p in pts)

        valid = int(re.search(r"<ValidPulseCount>(-?\d+)</ValidPulseCount>", text).group(1))
        total = int(re.search(r"<TotalPulseCount>(-?\d+)</TotalPulseCount>", text).group(1))

        values = [int(p) for p in pts]

        if fname == "zeroed.xml":
            assert all(v >= 0 for v in values)
            assert valid == total
            assert valid == sum(values)
        else:
            assert any(v < 0 for v in values)
            assert valid == sum(values)
            assert total - valid == 2351887 - 2349796

        assert "<SqrtFwhmCalibration>" in text
        m = re.search(r"<SqrtFwhmCalibration>.*?<Coefficients>(.*?)</Coefficients>", text, re.S)
        assert m is not None
        coeffs = [float(c) for c in re.findall(r"<Coefficient>([^<]+)</Coefficient>", m.group(1))]
        assert len(coeffs) == 3
        assert F.check(coeffs, 2048)[0] is True

    signed_values = [int(p) for p in re.findall(r"<DataPoint>([^<]*)</DataPoint>", (tmp_path / "signed.xml").read_text("utf-8"))]
    assert abs(sum(signed_values[352:2047])) < 3 * 500


def test_register_reproduces_reference():
    """Проверка воспроизведения эталонного преобразования h."""
    S, B, k, ls, lb = reg.load_inputs(BG, SAMPLE)
    anchor = reg.anchor37(B, S, verbose=False)
    mask = np.zeros(len(B), bool)
    mask[430:1701] = True
    r = reg.register(S, B, k, mask, anchor, n_ms=10, verbose=False)
    assert poly_eval(r["c"], 845.0) == pytest.approx(poly_eval(H_REF, 845.0), abs=0.05)
    assert 0.90 < r["chi2ndf"] < 1.00


def test_ecal_preset_reproduces_scale_and_flags_sum_peak_loo():
    """Проверка воспроизведения энергетической шкалы и LOO для суммарного пика."""
    cfg = json.loads((PRESETS / "labr3_ecal.json").read_text("utf-8"))
    S, B, k, ls, lb = reg.load_inputs(BG, SAMPLE)
    res = ops.calibrate(cfg, B, S, H_REF)
    assert res["coeffs"] == pytest.approx(ECAL_REF, rel=1e-5)
    sum_peak = next(p for p in res["points"] if p["name"] == "La138_sum1473")
    assert sum_peak["loo"] < -15


def test_cli_help_and_fwhm_device():
    """Проверка CLI: --help и подкоманда fwhm."""
    with pytest.raises(SystemExit) as e:
        cli_main(["--help"])
    assert e.value.code == 0

    argv = ["fwhm", "--device", "--points"] + ["%g:%g" % p for p in FW_POINTS]
    assert cli_main(argv) == 0


@pytest.mark.slow
def test_selftest_green_by_new_failures(capsys):
    """Проверка зелёного selftest по новым провалам."""
    rc = selftest.main(BG, SAMPLE)
    out = capsys.readouterr().out
    assert rc == 0
    assert "неожиданных 0" in out
    assert "мутаций с >=1 новым провалом 8 из 8" in out
    assert "fitB->A|alpha3:z_dip" in out
    assert "FAIL full|reg:dchi2" not in out
    assert "[INFO вне зачёта: согласованность с аргминимумом" in out


@pytest.mark.slow
def test_selftest_turns_red_without_known_list(tmp_path):
    """Проверка красного selftest при пустом списке известных провалов."""
    cfg = json.loads(selftest.DEFAULT_CFG.read_text("utf-8"))
    cfg["known_base_failures"] = []
    cfg["mutations"] = cfg["mutations"][:1]
    (tmp_path / "cfg.json").write_text(json.dumps(cfg, ensure_ascii=False), "utf-8")
    assert selftest.main(BG, SAMPLE, tmp_path / "cfg.json") == 1


def test_fwhm_bad_point_rejected_or_clear_error(tmp_path, capsys):
    """Плохая точка ПШПВ (нет пика в разности у канала 847): отбраковка при измерении, понятная ошибка при подгонке."""
    bad = [(143.0, 6.75), (212.0, 8.09), (338.0, 10.55), (847.0, 2.75)]
    with pytest.raises(ValueError, match="CheckCalibration"):
        F.fit_fwhm(2048, [p[0] for p in bad], [p[1] for p in bad])
    r = ops.make_net(BG, SAMPLE, H_REF, ECAL_REF, FW_POINTS, tmp_path / "z.xml", tmp_path / "s.xml")
    capsys.readouterr()
    res = ops.measure_fwhm(r["zeroed"], [143, 212, 338, 847], 8.0)
    out = capsys.readouterr().out
    assert [round(c) for c, f, e in res] == [143, 212, 339]
    assert "ch 847: ОТБРОШЕНА" in out
    argv = ["fwhm", "--points"] + ["%g:%g" % p for p in bad]
    assert cli_main(argv) == 2
