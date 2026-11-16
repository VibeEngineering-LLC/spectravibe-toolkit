"""Тесты параметров профилей детекторов, которые не требуют данных AS2PRO (идут на CI), и регрессия AtomNano 5 PRO
против выхода main 2ca4d15 на обеих парах фикстур."""

import hashlib
import json
import pathlib

import pytest

from gamma.calibration.intrinsic_shape import profiles
from gamma.calibration.intrinsic_shape import verdict as vd
from gamma.calibration.intrinsic_shape.cli import main as cli_main

FIX = pathlib.Path(__file__).resolve().parents[1] / "fixtures" / "intrinsic_shape"
REGRESS = json.loads((FIX / "regress_main_2ca4d15.json").read_text("utf-8"))
CS137_KEV = 661.66


@pytest.fixture(autouse=True)
def _reset_profile():
    profiles.activate(None)
    yield
    profiles.activate(None)


def _cs137_limit_ch() -> float:
    """Эталон AS2PRO: канал + 2*ПШПВ у точки, ближайшей к 1736 (Cs-137 661.7 кэВ)."""
    ref = profiles.load_preset("reference", "labr_as2pro")
    fwhm_points = ref["fwhm_points"]
    best = min(fwhm_points, key=lambda pair: abs(pair[0] - 1736))
    ch, fwhm = best
    return ch + 2 * fwhm


def test_as2pro_register_mask_above_cs137():
    """Нижняя граница маски register не ниже предела Cs-137 и ниже окна La_beta_edge selftest."""
    ref = profiles.load_preset("reference", "labr_as2pro")
    lo = ref["register"]["mask"][0]
    assert lo >= _cs137_limit_ch()
    st = profiles.load_preset("selftest", "labr_as2pro")
    assert lo < st["windows"]["La_beta_edge"][0]


def test_as2pro_selftest_mask_and_integrals_above_cs137():
    """full_mask и все интегралы selftest начинаются не ниже предела Cs-137."""
    st = profiles.load_preset("selftest", "labr_as2pro")
    lim = _cs137_limit_ch()
    assert st["full_mask"][0] >= lim
    for pair in st["integrals"]:
        assert pair[0] >= lim, pair


def test_profile_ids_and_aliases_unique():
    """id и алиасы всех профилей не пересекаются."""
    names: list[str] = []
    for p in profiles.list_profiles():
        names.append(p["id"])
        names.extend(p.get("aliases", []))
    assert len(names) == len(set(names)), names


def test_as2pro_alpha_windows_follow_sum_peak():
    """Окна α-линий не зависят от шкалы, а только от позиции сумм-пика."""
    with profiles.using("labr_as2pro"):
        wins = vd.load_windows_kev()
        a = vd.windows_to_channels(wins, [0.0, 0.377], 8192, sum_ch=3816.84)
        b = vd.windows_to_channels(wins, [0.0, 0.30], 8192, sum_ch=3816.84)
    for name in ("alpha1", "alpha2", "alpha3", "alpha4"):
        assert name in a and name in b
        assert (a[name]["lo"], a[name]["hi"]) == (b[name]["lo"], b[name]["hi"]), name


def test_as2pro_alpha_windows_match_selftest():
    """Вычисленные окна α1–α3 совпадают с окнами selftest в пределах 3 каналов."""
    with profiles.using("labr_as2pro"):
        ch = vd.windows_to_channels(vd.load_windows_kev(), [0.0, 0.377], 8192, sum_ch=3816.84)
    st = profiles.load_preset("selftest", "labr_as2pro")
    for name in ("alpha1", "alpha2", "alpha3"):
        lo, hi = ch[name]["lo"], ch[name]["hi"]
        slo, shi = st["windows"][name]
        assert abs(lo - slo) <= 3 and abs(hi - shi) <= 3, (name, lo, hi, slo, shi)


def test_scan_centroid_tolerance_by_profile(monkeypatch):
    """Допуск по каналам зависит от активного профиля."""
    from gamma.calibration.intrinsic_shape import register as reg

    # профиль по умолчанию (AtomNano): допуск 3 кан, отброс центроидов на 4 кан
    it1 = iter([100.0, 100.0, 104.0, 104.0, 104.0])
    def fake1(*a, **k):
        return {"centroid": next(it1), "centroid_err": 0.5}
    monkeypatch.setattr(reg, "fit_peak", fake1)
    profiles.activate(None)
    r = reg.scan_centroid(None, 100.0, 5.0, 10)
    assert r["n_ok"] == 3 and r["mean"] == pytest.approx(104.0)

    # профиль labr_as2pro: допуск 12 кан
    it2 = iter([100.0, 100.0, 104.0, 104.0, 104.0])
    def fake2(*a, **k):
        return {"centroid": next(it2), "centroid_err": 0.5}
    monkeypatch.setattr(reg, "fit_peak", fake2)
    with profiles.using("labr_as2pro"):
        r = reg.scan_centroid(None, 100.0, 5.0, 10)
    assert r["n_ok"] == 5 and r["mean"] == pytest.approx(102.4)


def test_ecal_scan_tolerance_by_profile(monkeypatch):
    """ecal.scan применяет допуск активного профиля."""
    from gamma.calibration.intrinsic_shape import ecal

    a = {"x0": 100.0, "w0": 5.0, "hw": 10}

    # профиль по умолчанию: допуск 3 кан
    it1 = iter([100.0, 100.0, 104.0, 104.0, 104.0])
    def fake1(*a, **k):
        return {"centroid": next(it1), "centroid_err": 0.5}
    monkeypatch.setattr(ecal, "fit_peak", fake1)
    mean, sig, cs = ecal.scan(None, a, "peak")
    assert len(cs) == 3 and mean == pytest.approx(104.0)

    # профиль labr_as2pro: допуск 12 кан
    it2 = iter([100.0, 100.0, 104.0, 104.0, 104.0])
    def fake2(*a, **k):
        return {"centroid": next(it2), "centroid_err": 0.5}
    monkeypatch.setattr(ecal, "fit_peak", fake2)
    with profiles.using("labr_as2pro"):
        mean, sig, cs = ecal.scan(None, a, "peak")
    assert len(cs) == 5 and mean == pytest.approx(102.4)


def test_scan_tolerance_default_without_profile_field(monkeypatch):
    """Профиль без поля scan_tol_ch (будущие профили) -- допуск по умолчанию 3 кан в register и ecal."""
    from gamma.calibration.intrinsic_shape import register as reg, ecal
    monkeypatch.setattr(profiles, "get", lambda pid=None: {"id": "без_допуска"})
    for mod, call in ((reg, lambda: reg.scan_centroid(None, 100.0, 5.0, 10)["n_ok"]),
                      (ecal, lambda: len(ecal.scan(None, {"x0": 100.0, "w0": 5.0, "hw": 10}, "peak")[2]))):
        it = iter([100.0, 100.0, 104.0, 104.0, 104.0])
        monkeypatch.setattr(mod, "fit_peak", lambda *a, **k: {"centroid": next(it), "centroid_err": 0.5})
        assert call() == 3, mod.__name__


def test_find_sum_window_ignores_decoy():
    """Полуокно поиска 6 % не захватывает более высокий пик на 7.5 % ниже прогноза."""
    import numpy as np
    from gamma.calibration.intrinsic_shape import autopreset as ap

    ref = ap.load_ref()
    xs = ref["anchors"]["sum"]["x"]
    a37 = ref["anchors"]["a37"]["x"]
    ch = np.arange(2048, dtype=float)
    counts = (
        10.0
        + 1000.0 * np.exp(-0.5 * ((ch - xs) / 4.6) ** 2)
        + 3000.0 * np.exp(-0.5 * ((ch - 0.925 * xs) / 4.6) ** 2)
    )
    # калибровки нет: прогноз по отношению к пику 37 кэВ = xs
    r = ap.find_sum(counts, a37, ref, None)
    assert abs(r["x"] - xs) < 3.0, r["x"]
@pytest.mark.parametrize("n_ch", ["2048", "8192"])
def test_atomnano_cli_outputs_match_main(n_ch, tmp_path):
    """CLI AtomNano 5 PRO даёт те же SHA-256, что и main 2ca4d15."""
    cfg = REGRESS["pairs"][n_ch]
    bg = str(FIX / cfg["bg"])
    sm = str(FIX / cfg["sample"])
    o = lambda name: str(tmp_path / name)

    assert cli_main(["register", "--bg", bg, "--sample", sm, "--out", o("reg.json")]) == 0
    assert cli_main(["ecal", "--bg", bg, "--sample", sm, "--h", o("reg.json"), "--out", o("ecal.json")]) == 0
    assert cli_main(["subtract", "--bg", bg, "--sample", sm, "--h", o("reg.json"),
                     "--ecal", o("ecal.json"), "--out-zeroed", o("net0.xml"),
                     "--out-signed", o("nets.xml"), "--date", "2026-11-15"]) == 0
    assert cli_main(["sample-cal", "--bg", bg, "--sample", sm,
                     "--out-zeroed", o("sc0.xml"), "--out-signed", o("scs.xml"),
                     "--date", "2026-11-15"]) == 0

    bad = []
    for name, want in cfg["sha256_lf"].items():
        got = hashlib.sha256((tmp_path / name).read_bytes().replace(b"\r\n", b"\n")).hexdigest()
        if got != want:
            bad.append(name)
    assert bad == [], bad
