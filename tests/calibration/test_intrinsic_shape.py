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


# --- автопресет (#LABR-8), устойчивая регистрация, ветвь шкалы ниже 100 кэВ (#LABR-9), мелочи (#LABR-3) ---
# Сгенерировано gen_code.py (qwen3.8:27b) по спецификации, 2026-10-10.
from gamma.calibration.intrinsic_shape import autopreset as ap
from gamma.calibration.intrinsic_shape import robust as rb
from gamma.calibration.intrinsic_shape import lowscale as lw
from gamma.calibration.intrinsic_shape import pipeline as pl
from gamma.calibration.intrinsic_shape import background as bgm
from gamma.calibration.intrinsic_shape.core import rebin_conserve


def _stretch(counts, a, b, n_dst):
    edges = a + b * (np.arange(len(counts) + 1) - 0.5)
    if not np.all(np.diff(edges) > 0):
        raise ValueError("edges not strictly increasing")
    return rebin_conserve(np.asarray(counts, float), edges, n_dst)


def test_autopreset_identity_on_reference_pair():
    """Проверка, что autopreset возвращает идентичность для эталонной пары спектров."""
    S, B, k, _, _ = reg.load_inputs(BG, SAMPLE)
    r = ap.autopreset(B, S, load_spectrum(BG)["cal"], verbose=False)
    assert r["identity"] is True
    assert r["ok"] is True
    assert r["register_kw"]["mask"] == [430, 1700]


def test_pipeline_register_equals_plain_register():
    """Регрессионный тест: pipeline.run_register совпадает с reg.register."""
    S, B, k, _, _ = reg.load_inputs(BG, SAMPLE)
    anchor = reg.anchor37(B, S, verbose=False)
    mask = np.zeros(len(B), dtype=bool)
    mask[430:1701] = True
    r0 = reg.register(S, B, k, mask, anchor, n_ms=10, verbose=False)
    r1 = pl.run_register(BG, SAMPLE, n_ms=10, verbose=False)
    assert max(abs(x - y) for x, y in zip(r0["c"], r1["h"])) == 0.0
    assert r1["out"]["robust"]["mode"] == "plain"


def test_pipeline_subtract_datapoints_equal_make_net(tmp_path):
    """Проверка, что pipeline.run_subtract генерирует те же DataPoint, что и ops.make_net."""
    a = ops.make_net(BG, SAMPLE, H_REF, ECAL_REF, FW_POINTS, tmp_path / "a.xml", tmp_path / "a_s.xml", date="2026-10-10")
    b = pl.run_subtract(BG, SAMPLE, H_REF, list(ECAL_REF), FW_POINTS, str(tmp_path / "b.xml"), str(tmp_path / "b_s.xml"), date="2026-10-10", xml_scale="main", verbose=False)
    for (x, y) in (("a.xml", "b.xml"), ("a_s.xml", "b_s.xml")):
        text_a = (tmp_path / x).read_text("utf-8")
        text_b = (tmp_path / y).read_text("utf-8")
        pts_a = re.findall(r"<DataPoint>([^<]*)</DataPoint>", text_a)
        pts_b = re.findall(r"<DataPoint>([^<]*)</DataPoint>", text_b)
        assert len(pts_a) == 2048
        assert pts_a == pts_b


@pytest.mark.parametrize("a,b", [(0.0, 4.0), (3.0, 4.1)])
def test_autopreset_recovers_synthetic_map(a, b):
    """Проверка восстановления синтетического отображения каналов autopreset'ом."""
    S, B, k, _, _ = reg.load_inputs(BG, SAMPLE)
    Bd = _stretch(B, a, b, 8192)
    Sd = _stretch(S, a, b, 8192)
    r = ap.autopreset(Bd, Sd, None, verbose=False)
    x = np.linspace(0, 1700, 50)
    err = np.max(np.abs(poly_eval(r["map"], x) - (a + b * x)))
    print("ошибка отображения, кан прибора:", err)
    assert err < 0.5 * b
    assert r["identity"] is False
    mapped_mask = r["register_kw"]["mask"]
    assert abs(mapped_mask[0] - (a + b * 430)) <= 1 + b
    assert abs(mapped_mask[1] - (a + b * 1700)) <= 1 + b


@pytest.mark.parametrize("name", ["a37", "sum", "edge"])
def test_autopreset_spoiled_anchor_turns_red(name):
    """Проверка, что испорченный якорь приводит к неудаче autopreset."""
    S, B, k, _, _ = reg.load_inputs(BG, SAMPLE)
    good = ap.autopreset(B, S, load_spectrum(BG)["cal"], verbose=False)
    assert good["ok"]
    x = good["anchors"][name]["x"]
    bad = ap.autopreset(B, S, load_spectrum(BG)["cal"], manual={name: x * 1.03}, verbose=False)
    assert bad["ok"] is False
    failed = [c["name"] for c in bad["checks"] if not c["ok"]]
    assert len(failed) >= 1


def test_robust_registration_with_sample_lines():
    """Проверка робастной регистрации на синтетическом образце с линиями."""
    S, B, k, _, _ = reg.load_inputs(BG, SAMPLE)
    n = len(B)
    x = np.arange(n, dtype=float)
    lines = (
        2.0e5 * np.exp(-0.5 * ((x - 425) / 4.0) ** 2) +
        3.0e4 * np.exp(-0.5 * ((x - 640) / 5.5) ** 2) +
        1.5e4 * np.exp(-0.5 * ((x - 1300) / 7.0) ** 2)
    )
    continuum = 30.0 * np.exp(-x / 700.0) * (x > 20)
    rng = np.random.default_rng(7)
    S2 = S + rng.poisson(lines + continuum).astype(float)   # правка Claude: шум -- только у добавки (S уже пуассоновский)
    S2[-1] = 0

    anchor = reg.anchor37(B, S2, verbose=False)
    mask = np.zeros(n, dtype=bool)
    mask[430:1701] = True
    fw = lambda c: float(np.sqrt(max(2.1 + 0.153 * c, 1.0)))

    plain = rb.register_robust(S2, B, k, mask, anchor, fw, k40_ch=None, mode="off", verbose=False, n_ms=10)
    rob = rb.register_robust(S2, B, k, mask, anchor, fw, k40_ch=None, mode="on", verbose=False, n_ms=10)

    d_plain = abs(poly_eval(plain["c"], 845.0) - poly_eval(H_REF, 845.0))
    d_rob = abs(poly_eval(rob["c"], 845.0) - poly_eval(H_REF, 845.0))
    print("d_plain:", d_plain, "d_rob:", d_rob)

    windows = json.loads((PRESETS / "labr3_selftest.json").read_text("utf-8"))["windows"]
    net_p, var_p = reg.apply(plain["c"], S2, B, k)
    net_r, var_r = reg.apply(rob["c"], S2, B, k)
    chk_p = rb.structure_checks(net_p, var_p, None, windows)
    chk_r = rb.structure_checks(net_r, var_r, rob["model"], windows, rob["excluded"])
    n_bad_p = sum(not c["ok"] for c in chk_p)
    n_bad_r = sum(not c["ok"] and c.get("n_used", 0) >= 5 for c in chk_r)
    print("n_bad_p:", n_bad_p, "n_bad_r:", n_bad_r)

    assert rob["mode"] == "robust"
    assert len(rob["excluded"]) >= 2
    assert d_rob < 0.3
    assert d_rob < d_plain
    assert n_bad_p >= 3
    assert n_bad_r <= 1


def test_nested_pair_detected():
    """Проверка детектирования вложенной пары (фон как образец)."""
    assert pl.nested_pair(BG, BG) != ""
    assert pl.nested_pair(BG, SAMPLE) == ""


def test_lowscale_recovers_known_s():
    """Проверка восстановления известного s в lowscale.fit_s."""
    C = list(ECAL_REF)
    s_true = -0.05
    energies = [35.75, 55.36, 63.5, 75.2, 88.34]
    pts = []
    for E in energies:
        L = lw.g_of_e(E, s_true)
        ch = poly_inverse(C, L)
        if E == 35.75:
            pts.append({"name": "intrinsic37", "E": 35.75, "dE": 0.5, "ch": ch, "ch_err": 0.01, "kind": "intrinsic37"})
        else:
            pts.append({"name": f"photon_{E}", "E": E, "dE": 0.05, "ch": ch, "ch_err": 0.01, "kind": "photon"})
    fit = lw.fit_s(C, pts)
    assert abs(fit["s"] - s_true) < 2e-3
    fit0 = lw.fit_s(C, pts, dK=0.0)
    # правка Claude: ступень 1.7 % на 37 кэВ = 0.61 кэВ при sigma репера 0.5 кэВ -- мутация видна на уровне ~1.3 sigma
    print("chi2 с изломом %.3g, без излома %.3g" % (fit["chi2"], fit0["chi2"]))
    assert fit0["chi2"] > fit["chi2"] + 1.0


def test_lowscale_continuity_and_becqmoni_check():
    """Проверка непрерывности composite_energy и проверки becqmoni."""
    C = list(ECAL_REF)
    fit = {"s": -0.05, "shape": "quad", "dK": 0.017}
    ch_j = poly_inverse(C, 100.0)
    e1 = lw.composite_energy(C, ch_j - 1e-3, fit)
    e2 = lw.composite_energy(C, ch_j + 1e-3, fit)
    assert abs(e2 - e1) < 0.01
    assert lw.becqmoni_ok(C, 2048)[0] is True
    assert lw.becqmoni_ok([0.0, 1.0, 0.0, 0.0, 0.0, 1e-20], 2048)[0] is False
    assert lw.becqmoni_ok([0.0, 1.0, -1e-3], 2048)[0] is False


def test_decompose_degenerate_ecal_clear_error(capsys):
    """Проверка явной ошибки при вырожденном ecal в decompose."""
    with pytest.raises(ValueError, match="сумм-пик"):
        bgm.check_ecal([0.0, 1.0, 0.0], bgm.load_any(BG)["counts"], "--ecal")
    rc = cli_main(["decompose", "--ref", str(BG), "--unshielded", str(SAMPLE), "--ecal", "0", "1", "0"])
    err = capsys.readouterr().err
    assert rc == 2
    assert "Ошибка" in err
    assert "Traceback" not in err
    bgm.check_ecal(list(load_spectrum(BG)["cal"]), bgm.load_any(BG)["counts"], "файл")


def test_fwhm_single_point_warns(capsys):
    """Проверка предупреждения при одной точке FWHM."""
    rc = cli_main(["fwhm", "--points", "100:5"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "ПРЕДУПРЕЖДЕНИЕ" in out
    assert "< 3" in out
    cli_main(["fwhm", "--points"] + ["%g:%g" % p for p in FW_POINTS])
    out2 = capsys.readouterr().out
    assert "ПРЕДУПРЕЖДЕНИЕ: точек" not in out2
